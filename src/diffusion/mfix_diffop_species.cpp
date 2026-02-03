#include <mfix_diffusion_op.H>
#include <mfix_fluid.H>
#include <mfix_run_on.H>

using namespace amrex;

MFIXDiffOpSpecies::
MFIXDiffOpSpecies ( int const a_nlev,
                    Vector<Geometry>                   const& a_geom,
                    Vector<BoxArray>                   const& a_grids,
                    Vector<DistributionMapping>        const& a_dmap,
                    Vector< const EBFArrayBoxFactory* >const& a_ebfactory,
                    int const a_ncomp,
                    MFIXFluidPhase& a_fluid,
                    MFIXBoundaryConditions& a_bcs,
                    Vector<BCRec> const& a_bcrec)
  : MFIXDiffusionOp(a_nlev, a_geom, a_ncomp, a_fluid, a_bcs, a_bcrec)
{
  readParameters();

  define( a_grids, a_dmap, a_ebfactory, a_ncomp );

  LPInfo info;

  info.setMaxCoarseningLevel(m_mg_max_coarsening_level);
  info.setAgglomerationGridSize(m_mg_agg_grid_size);

  m_matrix.reset( new  MLEBABecLap(m_geom, a_grids, a_dmap, info,
      a_ebfactory, a_ncomp));

  // It is essential that we set MaxOrder to 2 if we want to use the standard
  // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
  // The solver's default order is 3 and this uses three points for the gradient.
  m_matrix->setMaxOrder(2);

  // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
  m_matrix->setDomainBC( a_bcs.diff_species_lobc(), a_bcs.diff_species_hibc());

  // Mass diffusion coefficient:
  AMREX_ASSERT(m_b_defined == 0);
  m_b_cc.resize(a_nlev);
  for (int lev(0); lev < a_nlev; ++lev) {
    m_b_cc[lev].reset( new MultiFab( a_grids[lev], a_dmap[lev], a_ncomp, 1,
        MFInfo(), *a_ebfactory[lev]) );
  }
  m_b_defined = 1;

}


/***************************************************************************
 * Compute species diffusive fluxes: returns -Jk                           *
 *                                                                         *
 *  -Jk = -(ερDk)∇)Xk   <-- minus sign comes from AMReX formulation        *
 *                                                                         *
 *  Jk := diffusive flux of k-th species                                   *
 *   ε := volume fraction                                                  *
 *   ρ := density                                                          *
 *  Dk := diffusion coefficient of k-th species                            *
 *  Xk := k-th species mass fraction                                       *
 *                                                                         *
 ***************************************************************************/
void MFIXDiffOpSpecies::
computeFlux ( Vector< Array< MultiFab*, AMREX_SPACEDIM> > const& a_Jk,
              Vector< MultiFab*      > const& a_X_f,
              Vector< MultiFab const*> const& a_rho,
              Vector< MultiFab const*> const& a_epf)
{
  BL_PROFILE("MFIXDiffOpSpecies::ComputeFlux");
  AMREX_ASSERT( m_nlev  == a_X_f.size()       );

  AMREX_ASSERT( m_ncomp == m_fluid.nspecies()  );
  AMREX_ASSERT( m_ncomp == a_X_f[0]->nComp()  );

  // Set lapX_out to 0
  for(int lev(0); lev < m_nlev; lev++) {
    for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
      a_Jk[lev][dir]->setVal(0.0);
    }
  }

  // We want to return div (epf rho D_k grad)) phi
  m_matrix->setScalars(0.0, -1.0);

  // Compute the coefficients
  for (int lev(0); lev<m_nlev; lev++) {

    // Compute cell-center diffusion coefficient: ερDk
    setDiffCoeff(lev, a_rho, a_epf);

    // Set diffusion coefficient at MI to zero to prevent additional
    // diffusion into the domain from the mass inflow.
    avgDiffCoeffToFaces(lev);

    // Set BCoeffs
    m_matrix->setBCoeffs(lev, GetArrOfConstPtrs(m_b[lev]), MLMG::Location::FaceCentroid);

    // Set LevelBC
    m_matrix->setLevelBC(lev, GetVecOfConstPtrs(a_X_f)[lev]);
  }

  MLMG solver(*m_matrix);
  setSolverSettings(solver);

  // This ensures that ghost cells of sol are correctly filled when returned
  // from the solver
  solver.setFinalFillBC(true);

  solver.getFluxes(a_Jk, a_X_f, MLMG::Location::FaceCentroid);
}


/***************************************************************************
 * Compute divergence of species diffusive fluxes: returns -∇·Jk           *
 *                                                                         *
 *  -∇·Jk = ∇·(-Jk) = -∇·(ερDk)∇)Xk  <-- minus sign comes from AMReX       *
 *                                       returning -Jk from getFluxes()    *
 *                                                                         *
 *  Jk := diffusive flux of k-th species                                   *
 *   ε := volume fraction                                                  *
 *   ρ := density                                                          *
 *  Dk := diffusion coefficient of k-th species                            *
 *  Xk := k-th species mass fraction                                       *
 *                                                                         *
 ***************************************************************************/
void MFIXDiffOpSpecies::
computeDivJ ( Vector< MultiFab*      > const& divJ_out,
              Vector< Array< MultiFab*, AMREX_SPACEDIM> > const& a_Jk)
{
  BL_PROFILE("MFIXDiffOpSpecies::computeDivJ");

  // Set divJ_out to 0
  for(int lev(0); lev<m_nlev; lev++) {
    divJ_out[lev]->setVal(0.0);
  }

  // Species flux correction has been removed from this point.
  // Last develop commit before species flux correction was
  // b74b06253008954a5eb6b67640f477ce572ac113

  for (int lev(0); lev<m_nlev; lev++) {

    const auto& ebfact =
        dynamic_cast<EBFArrayBoxFactory const&>(divJ_out[lev]->Factory());

    if (ebfact.isAllRegular()) {

      amrex::computeDivergence(*divJ_out[lev], GetArrOfConstPtrs(a_Jk[lev]), m_geom[lev]);

    } else {

      amrex::EB_computeDivergence(*divJ_out[lev], GetArrOfConstPtrs(a_Jk[lev]),
          m_geom[lev], /*already_on_centroids=*/true);
    }
  }

  for(int lev(0); lev <m_nlev; lev++) {
    EB_set_covered(*divJ_out[lev], 0, divJ_out[lev]->nComp(), divJ_out[lev]->nGrow(), 0.);
  }
}


/***************************************************************************
 * Implicit species diffusion solve:                                       *
 *                                                                         *
 *                       (ερ - dt∇·(ερD)∇)X = ερX*                         *
 *                                                                         *
 *   ε := volume fraction @ t^(n+1)                                        *
 *   ρ := density @ t^(n+1)                                                *
 *   D := diffusion coefficient                                            *
 *   X := species mass fraction @ t^(n+1)  <--  (unknown)                  *
 *                                                                         *
 * ....................................................................... *
 *                                                                         *
 * The AMReX MLABecLaplacian linear solver class uses the canonical form:  *
 *                                                                         *
 *  (αA -β∇·B∇)φ = RHS     where     φ := unknown                          *
 *                                                                         *
 *  α :=  1;   A := ερ;                                                    *
 *  β := dt;   B := ερD;                                                   *
 *                                                                         *
 *  RHS = ερX*  <--  X* is the current value of X @ t^(n+1,*)              *
 *                                                                         *
 ***************************************************************************/
void MFIXDiffOpSpecies::
solve ( Vector< MultiFab      * > const& a_X_f,
        Vector< Array< MultiFab*, AMREX_SPACEDIM>> const& a_Jk,
        Vector< MultiFab const* > const& a_epf,
        Vector< MultiFab const* > const& a_rho,
        Real const a_dt)
{
  BL_PROFILE("MFIXDiffOpSpecies::solve");
  AMREX_ASSERT( m_nlev  == a_X_f.size() );

  if( m_verbose > 0) { Print() << "Diffusing species mass fractions ...\n"; }

  // Number of fluid species
  int const nspecies( m_ncomp );
  AMREX_ASSERT( nspecies == m_fluid.nspecies() );
  AMREX_ASSERT( nspecies == a_X_f[0]->nComp() );

  // Set alpha and beta
  m_matrix->setScalars(1.0, a_dt);

  for(int lev(0); lev<m_nlev; lev++) {

    // Create MultiFab to hold 'bulk density.'
    MultiFab epf_rho(a_epf[lev]->boxArray(), a_epf[lev]->DistributionMap(),
        1, 1, MFInfo(), a_epf[lev]->Factory());

    MultiFab::Copy(epf_rho,     *a_epf[lev], 0, 0, 1, epf_rho.nGrow());
    MultiFab::Multiply(epf_rho, *a_rho[lev], 0, 0, 1, epf_rho.nGrow());

    // This sets the coefficients
    m_matrix->setACoeffs (lev, epf_rho);

    // Compute cell-center diffusion coefficient: ερDk
    setDiffCoeff(lev, a_rho, a_epf);

    // Set diffusion coefficient at MI to zero to prevent additional
    // diffusion into the domain from the mass inflow.
    avgDiffCoeffToFaces(lev);

    m_matrix->setBCoeffs (lev, GetArrOfConstPtrs(m_b[lev]),
        MLMG::Location::FaceCentroid);

    // Zero these out just to have a clean start
    m_phi[lev]->setVal(0.0);
    m_rhs[lev]->setVal(0.0);

    for (MFIter mfi(*a_X_f[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      Box const& bx = mfi.growntilebox(IntVect(0));

      if (bx.ok()) {

        Array4<Real const> const& rhop  = epf_rho.const_array(mfi);
        Array4<Real const> const& X_f   = a_X_f[lev]->const_array(mfi);
        Array4<Real      > const& s_rhs = m_rhs[lev]->array(mfi);

        amrex::ParallelFor(bx, [rhop, X_f, s_rhs, nspecies]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for (int n(0); n < nspecies; ++n) {
            s_rhs(i,j,k,n) = rhop(i,j,k)*X_f(i,j,k,n);
          }
        });
      }
    }

    MultiFab::Copy(*m_phi[lev], *a_X_f[lev], 0, 0, nspecies, 1);
    m_matrix->setLevelBC(lev, GetVecOfConstPtrs(m_phi)[lev]);
  }

  MLMG solver(*m_matrix);
  setSolverSettings(solver);

  // This ensures that ghost cells of sol are correctly filled when returned from the solver
  solver.setFinalFillBC(true);

  solver.solve(GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs),
      m_mg_rtol, m_mg_atol);

  // Species mass fractions normalization and clamping
  for(int lev(0); lev<m_nlev; lev++) {

    const auto& factory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(a_epf[lev]->Factory());
    const auto& flags_fab = factory.getMultiEBCellFlagFab();

    for (MFIter mfi(*m_phi[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      Box const& bx = mfi.growntilebox(IntVect(0));

      if (bx.ok()) {

        Array4<Real> const& phi = m_phi[lev]->array(mfi);
        auto const& flags = flags_fab.const_array(mfi);

        amrex::ParallelFor(bx, [phi, nspecies, flags]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (!flags(i,j,k).isCovered()) {
            amrex::Real sum(0);
            for (int n(0); n < nspecies; ++n) {
              phi(i,j,k,n) = amrex::Clamp(phi(i,j,k,n), 0., 1.);
              sum += phi(i,j,k,n);
            }

            for (int n(0); n < nspecies; ++n) {
              phi(i,j,k,n) /= sum;
            }
          }
        });
      }
    }
  }

  solver.getFluxes(a_Jk, GetVecOfPtrs(m_phi), MLMG::Location::FaceCentroid);

  for(int lev(0); lev <m_nlev; lev++) {
    MultiFab::Copy(*a_X_f[lev], *m_phi[lev], 0, 0, nspecies, 1);
  }
}


/***************************************************************************
 * Compute species mass diffusion coefficient.                             *
 *                                                                         *
 *  Full diffusion term := (∇·(ερDk)∇)Xk                                   *
 *                                                                         *
 *  Diffusion coefficient := ερDk                                          *
 *                                                                         *
 *   ε := volume fraction                                                  *
 *   ρ := density                                                          *
 *  Dk := diffusion coefficient of k-th species                            *
 *                                                                         *
 ***************************************************************************/
void MFIXDiffOpSpecies::
setDiffCoeff ( int const a_lev,
               Vector< MultiFab const*> const& a_epf,
               Vector< MultiFab const*> const& a_rho)
{
  const auto& fluid_parms = m_fluid.parameters<run_on>();

  int const nspecies(m_ncomp);

  Box domain(m_geom[a_lev].Domain());

  const auto dlo = lbound(domain);
  const auto dhi = ubound(domain);

  const auto& bc_list = m_bcs.get_bc_list();

  m_b_cc[a_lev]->setVal(0.);

  int const nghost( m_b_cc[a_lev]->nGrow() );

  for (MFIter mfi(*m_phi[a_lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    Box const& bx = mfi.growntilebox(nghost);

    Array4<Real      > const& b_cc = m_b_cc[a_lev]->array(mfi);
    Array4<Real const> const& epf  = a_epf[a_lev]->const_array(mfi);
    Array4<Real const> const& rho  = a_rho[a_lev]->const_array(mfi);

    ParallelFor(bx, nspecies, [b_cc, epf, rho, fluid_parms]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      b_cc(i,j,k,n) = epf(i,j,k)*rho(i,j,k)*fluid_parms.get_D_g();
    });

    // Set boundary diffusion coefficients

    { Array4<int> const& bct_lo = bc_list.bc_ilo[a_lev]->array();
      Array4<int> const& bct_hi = bc_list.bc_ihi[a_lev]->array();

      ParallelFor(bx, m_ncomp, [b_cc, epf, rho, fluid_parms, dlo, dhi, bct_lo, bct_hi]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        if (i == dlo.x) {
          if (bct_lo(dlo.x-1,j,k,0) == BCList::pout) {
            b_cc(i-1,j,k,n) = epf(i,j,k)*rho(i,j,k)*fluid_parms.get_D_g();
          } else {
            b_cc(i-1,j,k,n) = 0.;
          }
        }
        if (i == dhi.x) {
          if (bct_hi(dhi.x+1,j,k,0) == BCList::pout) {
            b_cc(i+1,j,k,n) = epf(i,j,k)*rho(i,j,k)*fluid_parms.get_D_g();
          } else {
            b_cc(i+1,j,k,n) = 0.;
          }
        }
      });
    } // x-dir

    { Array4<int> const& bct_lo = bc_list.bc_jlo[a_lev]->array();
      Array4<int> const& bct_hi = bc_list.bc_jhi[a_lev]->array();

      ParallelFor(bx, m_ncomp, [b_cc, epf, rho, fluid_parms, dlo, dhi, bct_lo, bct_hi]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        if (j == dlo.y) {
          if (bct_lo(i,dlo.y-1,k,0) == BCList::pout) {
            b_cc(i,j-1,k,n) = epf(i,j,k)*rho(i,j,k)*fluid_parms.get_D_g();
          } else {
            b_cc(i,j-1,k,n) = 0.;
          }
        }
        if (j == dhi.y) {
          if (bct_hi(i,dhi.y+1,k,0) == BCList::pout) {
            b_cc(i,j+1,k,n) = epf(i,j,k)*rho(i,j,k)*fluid_parms.get_D_g();
          } else {
            b_cc(i,j+1,k,n) = 0.;
          }
        }
      });
    } // y-dir

    { Array4<int> const& bct_lo = bc_list.bc_klo[a_lev]->array();
      Array4<int> const& bct_hi = bc_list.bc_khi[a_lev]->array();

      ParallelFor(bx, m_ncomp, [b_cc, epf, rho, fluid_parms, dlo, dhi, bct_lo, bct_hi]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        if (k == dlo.z) {
          if (bct_lo(i,j,dlo.z-1,0) == BCList::pout) {
            b_cc(i,j,k-1,n) = epf(i,j,k)*rho(i,j,k)*fluid_parms.get_D_g();
          } else {
            b_cc(i,j,k-1,n) = 0.;
          }
        }
        if (k == dhi.z) {
          if (bct_hi(i,j,dhi.z+1,0) == BCList::pout) {
            b_cc(i,j,k+1,n) = epf(i,j,k)*rho(i,j,k)*fluid_parms.get_D_g();
          } else {
            b_cc(i,j,k+1,n) = 0.;
          }
        }
      });
    } // z-dir
  }
  m_b_cc[a_lev]->FillBoundary(m_geom[a_lev].periodicity());
}
