#include <mfix.H>
#include <mfix_filcc.H>

#include <AMReX_FillPatchUtil.H>
#include <mfix_pc.H>
#include <mfix_mf_helpers.H>

#include <mfix_dem.H>
#include <mfix_pic.H>

#include <mfix_deposition_op.H>

void mfix::pic_iteration (const bool apply_forces,
                          const bool update_parcels,
                          const bool use_taylor_approx,
                          const Real advance_vel_p,
                          Real dt_in,
                          RealVect& gravity_in,
                          Vector< Array<MultiFab*,3> >& vel_s_in,
                          Vector< MultiFab* >& ep_s_out,
                          EBFArrayBoxFactory* ebfactory_in,
                          const int ls_refinement,
                          const MultiFab* ls_phi)
{
  BL_PROFILE("mfix::pic_iteration()");

  if (nlev() > 2)
    amrex::Abort("For right now mfix::pic_iteration can only handle up to 2 levels");

  Vector< Array<MultiFab*,3> > vel_s_out(nlev());

  for (int lev(0); lev < nlev(); lev ++ ){

    const BoxArray&            pba = pc->ParticleBoxArray(lev);
    const DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

    const int pic_ncomp  = vel_s_in[lev][0]->nComp();
    const int pic_nghost = vel_s_in[lev][0]->nGrow();

    // Averaged particle velocity
    vel_s_out[lev][0] = new MultiFab(BoxArray(pba).surroundingNodes(0), pdm, pic_ncomp,
                          pic_nghost, MFInfo(), *(m_eb->particle_factory()[lev]));
    vel_s_out[lev][1] = new MultiFab(BoxArray(pba).surroundingNodes(1), pdm, pic_ncomp,
                          pic_nghost, MFInfo(), *(m_eb->particle_factory()[lev]));
    vel_s_out[lev][2] = new MultiFab(BoxArray(pba).surroundingNodes(2), pdm, pic_ncomp,
                          pic_nghost, MFInfo(), *(m_eb->particle_factory()[lev]));

    vel_s_out[lev][0]->setVal(0.0, 0, pic_ncomp, pic_nghost);
    vel_s_out[lev][1]->setVal(0.0, 0, pic_ncomp, pic_nghost);
    vel_s_out[lev][2]->setVal(0.0, 0, pic_ncomp, pic_nghost);

    // We must have ghost cells for each FAB so that a particle in one grid can spread
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (ep_s_out[lev]->nGrow() < 1)
      amrex::Error("Must have at least one ghost cell when in CalcVolumeFraction");

    ep_s_out[lev]->setVal(0.0,0,1,ep_s_out[lev]->nGrow());
  }

  const Geometry& gm  = Geom(0);
  const FabArray<EBCellFlagFab>* flags = nullptr;
  const MultiFab* volfrac = nullptr;

  for (int lev = 0; lev < nlev(); lev++) {

    // Use level 0 to define the EB factory
    if (lev == 0) {

      flags   = &(m_eb->particle_factory()[lev]->getMultiEBCellFlagFab());
      volfrac = &(m_eb->particle_factory()[lev]->getVolFrac());

    } else {

      Vector<int> ngrow = {1,1,1};
      EBFArrayBoxFactory* crse_factory;

      crse_factory = (makeEBFabFactory(gm, ep_s_out[lev]->boxArray(),
                                      ep_s_out[lev]->DistributionMap(),
                                      ngrow, EBSupport::volume)).release();

      flags   = &(crse_factory->getMultiEBCellFlagFab());
      volfrac = &(crse_factory->getVolFrac());

      delete crse_factory;
    }

    // Deposit particle volume to the grid
    pc->PICHydroStep(lev, apply_forces, update_parcels, use_taylor_approx,
                     advance_vel_p, dt_in, gravity_in, vel_s_in,
                     *ep_s_out[lev], vel_s_out[lev], volfrac, flags,
                     ebfactory_in, ls_refinement, ls_phi);
  }

  {
    // The deposition occurred on level 0, thus the next few operations
    // only need to be carried out on level 0.
    int lev0(0);

    MFIXDepositionOp DepOp(bc_list, pc);

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    DepOp.ApplyBC(lev0, ep_s_out[lev0]);

    // print_state(*ep_s_out[0],{1,23,1},0);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    ep_s_out[lev0]->SumBoundary(gm.periodicity());

    // const Real max_eps(1.1*m_pic.ep_cp());
    const Real max_eps(0.95);

    for (MFIter mfi(*ep_s_out[lev0],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // We don't want to do this for ghost cells
      const Box& bx = mfi.tilebox();

      if ( (*flags)[mfi].getType(amrex::grow(bx,0)) != FabType::covered &&
           (*flags)[mfi].getType(amrex::grow(bx,0)) != FabType::regular ) {

        Array4<Real> const& eps_arr = ep_s_out[lev0]->array(mfi);
        Array4<const EBCellFlag> const& flags_arr = (*flags)[mfi].array();

        // Clean up overfiled cells.
        amrex::ParallelFor(bx, [flags_arr, eps_arr, max_eps]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if(flags_arr(i,j,k).isSingleValued() && eps_arr(i,j,k) > max_eps) {
            eps_arr(i,j,k) = max_eps;
          }
        });

         }// Non-regular, Non-covered grid
       }// MFIter

    ep_s_out[lev0]->FillBoundary(gm.periodicity());

  }

  bool do_deposition = false;
  MFIX_CalcAvgSolidsVel(vel_s_out, do_deposition);

  for (int lev = 0; lev < nlev(); lev++) {
    // Swap old and new variables
    std::swap(vel_s_in[lev][0], vel_s_out[lev][0]);
    std::swap(vel_s_in[lev][1], vel_s_out[lev][1]);
    std::swap(vel_s_in[lev][2], vel_s_out[lev][2]);
    delete vel_s_out[lev][0];
    delete vel_s_out[lev][1];
    delete vel_s_out[lev][2];
  }

}
