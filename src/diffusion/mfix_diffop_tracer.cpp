#include <mfix_diffusion_op.H>

using namespace amrex;

MFIXDiffOpTracer::
MFIXDiffOpTracer ( int const a_nlev,
                   Vector<Geometry>                   const& a_geom,
                   Vector<BoxArray>                   const& a_grids,
                   Vector<DistributionMapping>        const& a_dmap,
                   Vector< const EBFArrayBoxFactory* >const& a_ebfactory,
                   int const a_ncomp, int const a_matrix_ncomp,
                   MFIXFluidPhase& a_fluid,
                   MFIXBoundaryConditions& a_bcs,
                   Vector<BCRec> const& a_bcrec)
  : MFIXDiffusionOp(a_nlev, a_geom, a_ncomp, a_fluid,  a_bcs, a_bcrec)
{
  readParameters();

  define( a_grids, a_dmap, a_ebfactory, a_matrix_ncomp);

  LPInfo info;

  info.setMaxCoarseningLevel(m_mg_max_coarsening_level);
  info.setAgglomerationGridSize(m_mg_agg_grid_size);

  m_matrix.reset( new MLEBABecLap(m_geom, a_grids, a_dmap, info,
      a_ebfactory, a_matrix_ncomp));

  // It is essential that we set MaxOrder to 2 if we want to use the standard
  // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
  // The solver's default order is 3 and this uses three points for the gradient.
  m_matrix->setMaxOrder(2);

  // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
  m_matrix->setDomainBC(a_bcs.diff_scal_lobc(), a_bcs.diff_scal_hibc());
}


void MFIXDiffOpTracer::
computeLap ( Vector< MultiFab      *> const& a_laps,
             Vector< MultiFab      *> const& a_tracer,
             Vector< MultiFab const*> const& /*a_epf*/,
             Vector< MultiFab const*> const& /*a_rho*/)
{
  BL_PROFILE("MFIXDiffOpTracer::computeLap");
  AMREX_ASSERT( m_nlev  == a_tracer.size()      );
  AMREX_ASSERT( m_ncomp == a_tracer[0]->nComp() );
  AMREX_ASSERT( m_b_defined == 1 );

  Vector< MultiFab> tracer(m_nlev);

  for(int lev(0); lev<m_nlev; lev++) {

    a_laps[lev]->setVal(0.0);

    const BoxArray            ba = a_tracer[lev]->boxArray();
    const DistributionMapping dm = a_tracer[lev]->DistributionMap();
    const FabFactory<FArrayBox>& ebfac = a_tracer[lev]->Factory();

    tracer[lev].define(ba, dm, m_ncomp, 1, MFInfo(), ebfac);
    MultiFab::Copy(tracer[lev], *a_tracer[lev], 0, 0, m_ncomp, 1);

    EB_set_covered(tracer[lev], 0, 1, m_ncomp, 0.);
  }

  // We want to return div (mu grad)) phi
  m_matrix->setScalars(0.0, -1.0);

  // Compute the coefficients
  for (int lev(0); lev<m_nlev; lev++) {

    if ( m_const_coeff ) {

      AMREX_ASSERT( m_b_const.size() == m_ncomp );

      for(int idim(0); idim < 3; idim++) {
        for(int n(0); n < m_ncomp; n++) {
          m_b[lev][idim]->setVal(m_b_const[n],n,1);
        }
      } // idim

    } else {

      int const b_ncomp(1);  // Other changes are needed if we want m_ncomp
      AMREX_ASSERT( m_b_cc.size() == b_ncomp );

      Vector<BCRec> dummy_bc(AMREX_SPACEDIM*b_ncomp);

      EB_interp_CellCentroid_to_FaceCentroid( *m_b_cc[lev], GetArrOfPtrs(m_b[lev]),
          0, 0, b_ncomp, m_geom[lev], /*m_bcrec*/ dummy_bc);
    }

#if 0
    if (eb_is_dirichlet) {
      m_matrix->setEBDirichlet(lev, *phi_eb[lev], mu_s);
    }
#endif

    m_matrix->setBCoeffs(lev, GetArrOfConstPtrs(m_b[lev]),MLMG::Location::FaceCentroid);
    m_matrix->setLevelBC(lev, &tracer[lev]);
  }

  MLMG solver(*m_matrix);

  solver.apply(a_laps, GetVecOfPtrs(tracer));

  for(int lev(0); lev<m_nlev; lev++) {
    EB_set_covered(*a_laps[lev], 0, a_laps[lev]->nComp(), a_laps[lev]->nGrow(), 0.);
  }

}


//
// Implicit solve for tracer diffusion
//
void MFIXDiffOpTracer::
solve ( Vector< MultiFab      * > const& a_tracer,
        Vector< MultiFab const* > const& a_epf,
        Vector< MultiFab const* > const& a_rho,
        Real const a_dt)
{
  BL_PROFILE("MFIXDiffOpTracer::solve");
  AMREX_ASSERT( m_nlev  == a_tracer.size()      );
  AMREX_ASSERT( m_ncomp == a_tracer[0]->nComp() );

  // Update the coefficients of the matrix going into the solve based on the current state of the
  // simulation. Recall that the relevant matrix is
  //
  //      alpha a - beta div ( b grad )   <--->   rho - dt div ( mu_s grad )
  //
  // So the constants and variable coefficients are:
  //
  //      alpha: 1
  //      beta: dt
  //      a: ro
  //      b: eta

  // Set alpha and beta
  m_matrix->setScalars(1.0, a_dt);

  for(int lev(0); lev<m_nlev; lev++) {

    MultiFab epf_rho(a_epf[lev]->boxArray(), a_epf[lev]->DistributionMap(),
        1, 1, MFInfo(), a_epf[lev]->Factory());

    MultiFab::Copy(epf_rho, *a_epf[lev], 0, 0, 1, epf_rho.nGrow());
    MultiFab::Multiply(epf_rho, *a_rho[lev], 0, 0, 1, epf_rho.nGrow());

    // This sets the coefficients
    m_matrix->setACoeffs(lev, epf_rho);

    if ( m_const_coeff ) {

      AMREX_ASSERT( m_b_const.size() == m_ncomp );

      for(int idim(0); idim < 3; idim++) {
        for(int n(0); n < m_ncomp; n++) {
          m_b[lev][idim]->setVal(m_b_const[n],n,1);
        }
      } // idim

    } else {

      int const b_ncomp(1);  // Other changes are needed if we want m_ncomp
      AMREX_ASSERT( m_b_cc.size() == b_ncomp );

      Vector<BCRec> dummy_bc(AMREX_SPACEDIM*b_ncomp);

      EB_interp_CellCentroid_to_FaceCentroid( *m_b_cc[lev], GetArrOfPtrs(m_b[lev]),
          0, 0, b_ncomp, m_geom[lev], /*m_bcrec*/ dummy_bc);
    }

    m_matrix->setBCoeffs(lev, GetArrOfConstPtrs(m_b[lev]), MLMG::Location::FaceCentroid);

    // Zero these out just to have a clean start because they have 3 components
    //      (due to reuse with velocity solve)
    m_phi[lev]->setVal(0.0);
    m_rhs[lev]->setVal(0.0);

    // Set the right hand side to equal rhs
    MultiFab::Copy((*m_rhs[lev]),(*a_tracer[lev]), 0, 0, m_ncomp, 0);

    // Multiply rhs by epf*rho  -- we are solving
    //
    //  epf rho s_star = epf rho s_old + dt ( - div (epf rho U s) + epf rho div (mu (grad s)) )
    //
    for(int n(0); n < m_ncomp; n++) {
      MultiFab::Multiply((*m_rhs[lev]), epf_rho, 0, n, 1, 0);
    }

    MultiFab::Copy(*m_phi[lev],*a_tracer[lev], 0, 0, m_ncomp, 1);
    m_matrix->setLevelBC(lev, GetVecOfConstPtrs(m_phi)[lev]);
  }

  if(m_verbose > 0) { Print() << "Diffusing tracers one at a time ...\n"; }

  MLMG solver(*m_matrix);
  setSolverSettings(solver);

  // This ensures that ghost cells of sol are correctly filled when returned from the solver
  solver.setFinalFillBC(true);

  solver.solve(GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs),
      m_mg_rtol, m_mg_atol);

  for(int lev(0); lev<m_nlev; lev++) {
    MultiFab::Copy(*a_tracer[lev], *m_phi[lev], 0, 0, m_ncomp, 1);
  }
}
