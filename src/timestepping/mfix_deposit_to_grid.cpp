#include <mfix.H>
#include <mfix_utils.H>

#include <mfix_deposition_op.H>
#include <mfix_diffusion_op.H>

using namespace amrex;

Real mfix::
deposit_volume_to_grid ( Real const a_time )
{
  BL_PROFILE("mfix::deposit_volume_to_grid()");

  // Print what if we use the regular factory?
  MFIXDepOpVolAvg volDepOp(bc_list, pc, m_eb->particle_factory()[0].get(), leveldata().epf());

  volDepOp.Deposit();
  volDepOp.Redistribute();
  volDepOp.CopyToDest();

  // At this point, we have the particle volume on the fluid grid (ep_s).
  // We will diffuse it first, then convert it to epf.
  if ( m_coupling.FilterDeposition() ) {

    DepositionFilter const* depop_filter = m_coupling.getDepositionFilter();
    if ( !depop_filter->ConstantSize() ) {
      ComputeVariableFilter( leveldata().epf_const() );
    }

    diffOpVoidFrac()->solve( leveldata().epf() );
  }

  for (int lev(0); lev<nlev(); ++lev) {
    // Now define this epf = (1 - particle_vol)
    leveldata().epf(lev)->mult(-1.0, leveldata().epf(lev)->nGrow());
    leveldata().epf(lev)->plus( 1.0, leveldata().epf(lev)->nGrow());

    // We set epf to 1 rather than 0 in covered cells so that when we divide
    // by epf, we don't have to protect against division by 0.
    EB_set_covered(*(leveldata().epf(lev)),1.0);

    bcs().set_epf_bcs(lev, a_time, leveldata().epf(lev), BCType::ext_dir);
  }


  // Sum up all the values of epf[lev], weighted by each cell's EB volfrac
  // Note epf = 1 - particle_volume / this_cell_volume where
  //    this_cell_volume = (volfrac * dx * dy * dz)
  // When we define the sum we add up (epf * volfrac) so that the total sum
  //    does not depend on whether a particle is in a full or cut cell.

  // Compute the monitor on coarsest level
  const int lev = 0;
  EulerianMonitor::VolumeIntegral monitor(leveldata(), m_eb->factory(), geom[lev].ProbDomain());

  // we compute the monitor on the only component of the epf Multifab
  const int comp = 0;
  auto domain_vol = monitor.volume_integral(lev, {leveldata().epf(lev)}, {comp});

  // sum_vol_orig is a vector with only one component
  const int var = 0;
  return domain_vol[var];
}


void mfix::
deposit_forces_to_grid ( Real const a_dt )
{
  BL_PROFILE("mfix::deposit_forces_to_grid()");

  MFIXDepOpTxfr txfrDepOp(bc_list, pc, m_eb->particle_factory()[0].get(), leveldata().txfr(), a_dt);

  if( m_coupling.include_virtual_mass() ) { txfrDepOp.setVirtualMass(); }

  // Deposit the interphase transfer forces to the grid and reduce to level-0.
  txfrDepOp.Deposit();
  txfrDepOp.Redistribute();
  txfrDepOp.CopyToDest();

  // Apply mean field diffusion
  if( m_coupling.FilterDeposition() ) {
    diffOpTxfr()->solve( leveldata().txfr() );
  }

}


void mfix::
ComputeVariableFilter ( Vector< MultiFab const*> const& a_eps )
{

  DepositionFilter const* depop_filter = m_coupling.getDepositionFilter();

  Real const sample_size( depop_filter->size() );
  Real const mineps( depop_filter->mineps() );

  // Compute an averaged particle radius
  Vector<std::unique_ptr<MultiFab>> diffCoeff(nlev());
  Vector<std::unique_ptr<MultiFab>> avgRadius(nlev());

  for (int lev(0); lev<nlev(); ++lev) {
    diffCoeff[lev].reset( new MultiFab(grids[lev], dmap[lev], 1, 1, MFInfo(), m_eb->Factory(lev)));
    avgRadius[lev].reset( new MultiFab(grids[lev], dmap[lev], 1, 1, MFInfo(), m_eb->Factory(lev)));

    diffCoeff[lev]->setVal(0.);
    avgRadius[lev]->setVal(0.);
  }

  MFIXDepOpVolAvg radDepOp(bc_list, pc, m_eb->particle_factory()[0].get(), GetVecOfPtrs(avgRadius));

  // Set flag to deposit radius
  radDepOp.setDepComp( SoArealData::radius );

  // Deposit, compute the average, and copy into avgRadius. We skip redistributione
  // because there is no need to conserver "radius"

  radDepOp.Deposit();
  radDepOp.Average();
  radDepOp.CopyToDest();

  // Compute the variable filter size
  for (int lev(0); lev<nlev(); ++lev) {

    Real min_length = std::numeric_limits<Real>::max();

    for (int idir(0); idir < AMREX_SPACEDIM; ++idir)
    { min_length = amrex::min(min_length, geom[lev].ProbLength(idir)); }

    Real const third = 1.0/3.0;
    Real const inv_16ln2 = 1.0 / (16.0*std::log(2.0));

    for (MFIter mfi(*diffCoeff[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& tbox = mfi.tilebox();

      Array4<Real      > const& coeff  = diffCoeff[lev]->array(mfi);
      Array4<Real const> const& avg_rp = avgRadius[lev]->const_array(mfi);
      Array4<Real const> const& eps    = a_eps[lev]->const_array(mfi);

      ParallelFor(tbox, [coeff, avg_rp, eps, sample_size, inv_16ln2, mineps, min_length, third]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real const dp =  (2.0*avg_rp(i,j,k)) + std::numeric_limits<Real>::min();
        Real const eps_ = amrex::max(eps(i,j,k), mineps);

        Real const delta = std::pow(sample_size*(dp*dp*dp)/eps_, third);

        coeff(i,j,k) = amrex::min( delta*delta*inv_16ln2, min_length);

      });

    } // MFIter

    diffCoeff[lev]->FillBoundary(geom[lev].periodicity());

  } // lev

  diffOpVoidFrac()->setDiffCoeff( GetVecOfConstPtrs(diffCoeff) );
  diffOpTxfr()->setDiffCoeff( GetVecOfConstPtrs(diffCoeff) );

}

void mfix::
average_pc_data_to_grid ( Vector< MultiFab* > const& a_avg_data )
{
  for (int lev(0); lev < nlev(); ++lev) {
    a_avg_data[lev]->setVal(0.);
  }

  MFIXDepOpVolAvg depOp(bc_list, pc, (eb()->pc_factory_const())[0], a_avg_data);

  // Set flags to deposit
  depOp.setDepComp( avg_pc_parms().comps() );

  // Deposit, compute the average, and copy into a_avg_data.
  // We skip redistributione because there is no need to conservere
  depOp.Deposit();
  depOp.Average();
  //depOp.Redistribute();
  depOp.CopyToDest();

}
