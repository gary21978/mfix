#include <mfix_diffusion_op.H>

#include <mfix_fluid.H>
#include <mfix_solvers.H>

#include <AMReX_EB_utils.H>

using namespace amrex;
using namespace Solvers;

MFIXDiffOpEnergy::
MFIXDiffOpEnergy ( int const a_nlev,
                   Vector<Geometry>                   const& a_geom,
                   Vector<BoxArray>                   const& a_grids,
                   Vector<DistributionMapping>        const& a_dmap,
                   Vector< const EBFArrayBoxFactory* >const& a_ebfactory,
                   MFIXFluidPhase& a_fluid,
                   MFIXBoundaryConditions& a_bcs,
                   Vector<BCRec> const& a_bcrec)
  : MFIXDiffusionOp(a_nlev, a_geom, /*ncomp=*/1, a_fluid, a_bcs, a_bcrec)
{
  readParameters();

  int const matrix_ncomp(1);
  define( a_grids, a_dmap, a_ebfactory, matrix_ncomp );

  LPInfo info;

  // Turning off multigrid coarsening for temperature matrix
  // since its only called with "apply" not "solve"
  info.setMaxCoarseningLevel(0);
  info.setAgglomerationGridSize(m_mg_agg_grid_size);

  m_matrix.reset( new MLEBABecLap(m_geom, a_grids, a_dmap, info,
      a_ebfactory, matrix_ncomp) );

  // It is essential that we set MaxOrder to 2 if we want to use the standard
  // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
  // The solver's default order is 3 and this uses three points for the gradient.
  m_matrix->setMaxOrder(2);

  // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
  m_matrix->setDomainBC(a_bcs.diff_temperature_lobc(), a_bcs.diff_temperature_hibc());

  // Thermal diffusion coefficient: εκ
  AMREX_ASSERT(m_b_defined == 0);
  m_b_cc.resize(a_nlev);
  for (int lev(0); lev < a_nlev; ++lev) {
    m_b_cc[lev].reset( new MultiFab( a_grids[lev], a_dmap[lev], 1, 1,
        MFInfo(), *a_ebfactory[lev]) );
  }
  m_b_defined = 1;

}

/***************************************************************************
 * Compute explicit thermal conduction:  (∇·εκ∇)T                          *
 *                                                                         *
 *   ε := volume fraction                                                  *
 *   κ := thermal conductivity (diffusion coefficient)                     *
 *   T := temperature                                                      *
 *                                                                         *
 ***************************************************************************/
void MFIXDiffOpEnergy::
computeLapT ( Vector< MultiFab*      > const& a_lapT,
              Vector< MultiFab*      > const& a_Tf,
              Vector< MultiFab const*> const& a_epf,
              Vector< MultiFab const*> const& a_Teb)
{
  BL_PROFILE("MFIXDiffnOEnergyp::computeLapT");

  Vector< MultiFab> temperature(m_nlev);

  for(int lev = 0; lev<m_nlev; lev++) {

    a_lapT[lev]->setVal(0.0);

    const BoxArray            ba = a_epf[lev]->boxArray();
    const DistributionMapping dm = a_epf[lev]->DistributionMap();
    const FabFactory<FArrayBox>& ebfac = a_epf[lev]->Factory();

    temperature[lev].define(ba, dm, 1, 1, MFInfo(), ebfac);
    MultiFab::Copy(temperature[lev], *a_Tf[lev], 0, 0, 1, 1);

    EB_set_covered(temperature[lev], 0, 1, 1, 0.);
  }

  // We want to return div (a_epf k_g grad)) a_Tf
  m_matrix->setScalars(0.0, -1.0);

  // Compute the coefficients
  for (int lev = 0; lev<m_nlev; lev++) {

    // Compute thermal diffusion coefficient: epf*kappa
    setDiffCoeff(lev, a_epf, GetVecOfConstPtrs(a_Tf));

    // Set diffusion coefficient at MI to zero to prevent additional
    // diffusion of heat into the domain from the mass inflow.
    avgDiffCoeffToFaces(lev);

    if (m_set_eb_dirichlet) {

      AMREX_ASSERT( a_Teb[lev] != nullptr );
      AMREX_ASSERT( a_Teb[lev]->nComp() == m_ncomp);

      // The following is a WIP in AMReX
      //m_matrix->setPhiOnCentroid();

      setEBDiffCoeff(lev, a_Teb);

      m_matrix->setEBDirichlet(lev, *a_Teb[lev], *m_b_eb[lev]);
    }

    m_matrix->setBCoeffs(lev, GetArrOfConstPtrs(m_b[lev]),
        MLMG::Location::FaceCentroid);

    m_matrix->setLevelBC(lev, &temperature[lev]);
  }

  MLMG solver(*m_matrix);

  solver.apply(a_lapT, GetVecOfPtrs(temperature));

  for(int lev = 0; lev<m_nlev; lev++) {
    EB_set_covered(*a_lapT[lev], 0, a_lapT[lev]->nComp(), a_lapT[lev]->nGrow(), 0.);
  }
}


/***************************************************************************
 * Compute the energy change from interdiffusion of species:  -∇·hkJk      *
 *                                                                         *
 *  -∇·hkJk = ∇·(-hkJk) = -(∇·hk(ερDk)∇)Xk <-- minus sign comes from AMReX *
 *                                             AMReX getFluxes()           *
 *  Jk := diffusive flux of k-th species                                   *
 *  hk := specific enthalpy of k-th species                                *
 *   ε := volume fraction                                                  *
 *   ρ := density                                                          *
 *  Dk := diffusion coefficient of k-th species                            *
 *  Xk := k-th species mass fraction                                       *
 *                                                                         *
 ***************************************************************************/
void MFIXDiffOpEnergy::
computeDivhJ ( Vector< MultiFab* > const& divhJ_out,
               Vector< Array< MultiFab*, AMREX_SPACEDIM> > const& a_hf_fc,
               Vector< Array< MultiFab*, AMREX_SPACEDIM> > const& a_Jk,
               Vector< MultiFab const* > const& a_Tf,
               int const update_enthalpies)
{
  BL_PROFILE("MFIXDiffOpEnergy::ComputeDivhJ");

  // Set divhJ_out to 0
  for(int lev = 0; lev<m_nlev; lev++) {
    divhJ_out[lev]->setVal(0.0);
  }

  // Number of fluid species
  const int nspecies = m_fluid.nspecies();
  const auto fluid_props = m_fluid.props.data<run_on>();

  for (int lev = 0; lev<m_nlev; ++lev) {

    Array<MultiFab, AMREX_SPACEDIM> hJ_fc;

    for (int dim(0); dim < AMREX_SPACEDIM; ++dim) {

      const IntVect dimVec = IntVect::TheDimensionVector(dim);
      const BoxArray& ba = amrex::convert(m_rhs[lev]->boxArray(), dimVec);

      hJ_fc[dim].define(ba, m_rhs[lev]->DistributionMap(),
          1, 0, MFInfo(), m_rhs[lev]->Factory());

      hJ_fc[dim].setVal(0.);

      for (MFIter mfi(hJ_fc[dim],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real      > const& hJ = hJ_fc[dim].array(mfi);
        Array4<Real      > const& h  = a_hf_fc[lev][dim]->array(mfi);

        Array4<Real const> const& Tf = a_Tf[lev]->const_array(mfi);
        Array4<Real const> const& J  = a_Jk[lev][dim]->const_array(mfi);

        int const update_h( update_enthalpies );

        ParallelFor(bx, [nspecies, dimVec, update_h, Tf, h, J, hJ, fluid_props]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          hJ(i,j,k) = 0.0;

          for (int n(0); n < nspecies; ++n) {

            if (update_h) {

              IntVect ijk = IntVect{i,j,k};

              // "upwind" enthalpy to cell face
              if (J(i,j,k,n) > 0.) { ijk -= dimVec; }

              h(i,j,k,n) = fluid_props.enthalpy(n,Tf(ijk));
            }

            hJ(i,j,k) += (h(i,j,k,n) * J(i,j,k,n));
          }
        });
      }
    }

    divhJ_out[lev]->setVal(0.);

    const auto& ebfact =
      static_cast<amrex::EBFArrayBoxFactory const&>(a_Tf[lev]->Factory());

    if (ebfact.isAllRegular()) {

      computeDivergence( *divhJ_out[lev], GetArrOfConstPtrs(hJ_fc), m_geom[lev]);

    } else {

      const bool already_on_centroids = true;
      EB_computeDivergence(*divhJ_out[lev], GetArrOfConstPtrs(hJ_fc),
          m_geom[lev], already_on_centroids);
    }
  }

  for(int lev(0); lev<m_nlev; ++lev) {
    EB_set_covered(*divhJ_out[lev], 0, divhJ_out[lev]->nComp(), divhJ_out[lev]->nGrow(), 0.);
  }
}

/***************************************************************************
 * Implicit energy diffusion solve:                                        *
 *                                                                         *
 *                    ερh(T) - dt∇·(εκ)∇)T = ερh*                          *
 *                                                                         *
 *     ε := volume fraction @ t^(n+1)                                      *
 *     ρ := density @ t^(n+1)                                              *
 *    h* := specific enthalpy @ t^(n+1,*)                                  *
 *     κ := thermal conductivity (diffusion coefficient)                   *
 *     T := temperature @ t^(n+1) <-- (unknown)                            *
 *  h(T) := specific enthalpy @ t^(n+1)    <-- (function of unknown T)     *
 *                                                                         *
 *  The equation is solved iteratively for temperature with convergence    *
 *  determined when max(abs(T^(k+1) - T^k)) < max(|T*|)*tolerance          *
 *                                                                         *
 *    T^(k=0) = T*  and the 'new' enthalpy approximated by                 *
 *                                                                         *
 *    h(T^(k+1)) ≈ h(T^k) + Cp(T^k)*[T^(k+1) - T^k]                        *
 *                                                                         *
 *  Substituting these definitions into the equation above,                *
 *                                                                         *
 *    ερh(T^k) + ερCp(T^k)*[T^(k+1) - T^k] - dt∇·(εκ)∇)T^(k+1) = ερh(T*)   *
 *                                                                         *
 *  Rearranging:                                                           *
 *                                                                         *
 *    [ερCp(T^k)- dt∇·(εκ)∇)]T^(k+1) = ερh* - ερh(T^k) + ερCp(T^k)*T^k     *
 *                                                                         *
 * ....................................................................... *
 *                                                                         *
 * The AMReX MLABecLaplacian linear solver class uses the canonical form:  *
 *                                                                         *
 *  [αA - β∇·B∇]φ = RHS     where    φ := unknown                          *
 *                                                                         *
 *  α :=  1;   A := ερCp(T^k);                                             *
 *  β := dt;   B := εκ;                                                    *
 *                                                                         *
 *  RHS = ερh* - ερh(T^k) + ερCp(T^k)*T^k                                  *
 *                                                                         *
 ***************************************************************************/
void MFIXDiffOpEnergy::
solve ( Vector< MultiFab      * > const& a_Tf,
        Vector< MultiFab      * > const& a_hf,
        Vector< MultiFab const* > const& a_epf,
        Vector< MultiFab const* > const& a_rho,
        Vector< MultiFab const* > const& X_gk,
        Vector< MultiFab const* > const& a_Teb,
        Real const dt,
        Real const /*abstol*/,
        Real const reltol,
        int const maxiter)
{
  BL_PROFILE("MFIXDiffOpEnergy::solve");

  if(m_verbose > 0) { Print() << "Diffusing temperature ...\n"; }

  int const fluid_is_a_mixture = m_fluid.isMixture();

  const auto fluid_props = m_fluid.props.data<run_on>();

  amrex::Vector<amrex::MultiFab*> A(m_nlev, nullptr);

  for (int lev(0); lev<m_nlev; ++lev) {

    A[lev] = new MultiFab(m_rhs[lev]->boxArray(), m_rhs[lev]->DistributionMap(),
        1, 1, MFInfo(), m_rhs[lev]->Factory());

    A[lev]->setVal(0.);
  }

  // **************************************************************************
  // Define the norm function
  // **************************************************************************

  auto norm0 = [&] (const Vector<MultiFab*>& vec_of_MFs) -> Real
  {
    Vector<Real> vec_of_norms(vec_of_MFs.size(), 0.);
    std::transform(vec_of_MFs.begin(), vec_of_MFs.end(), vec_of_norms.begin(),
      [](MultiFab* MF) { return MF->norm0(0, 0, false, true); });

    Real norm(0);
    for (auto value: vec_of_norms)
      norm = amrex::max(norm, value);

    ParallelDescriptor::ReduceRealMax(norm);

    return norm;
  };

  // **************************************************************************
  // Do the Newton iterations
  // **************************************************************************

  amrex::Vector<amrex::MultiFab*> update(m_nlev, nullptr);
//  amrex::Vector<amrex::MultiFab*> residue(m_nlev, nullptr);

  for (int lev(0); lev<m_nlev; ++lev) {
    update[lev] = new MultiFab(m_rhs[lev]->boxArray(), m_rhs[lev]->DistributionMap(),
        1, 0, MFInfo(), a_epf[lev]->Factory());

    update[lev]->setVal(0.);

//    residue[lev] = new MultiFab(m_rhs[lev]->boxArray(), m_rhs[lev]->DistributionMap(),
//        1, 0, MFInfo(), m_rhs[lev]->Factory());
//    residue[lev]->setVal(0.);
  }

  int iter(0);
  const amrex::Real update_rel_tol = reltol*norm0(a_Tf);
  //const amrex::Real residue_rel_tol = reltol*norm0(residue);

  if (m_set_eb_dirichlet) {

    for (int lev(0); lev<m_nlev; ++lev) {

      AMREX_ASSERT( a_Teb[lev] != nullptr );
      AMREX_ASSERT( a_Teb[lev]->nComp() == m_ncomp);

      // The following is a WIP in AMReX
      //m_matrix->setPhiOnCentroid();
      setEBDiffCoeff(lev, a_Teb);
    }
  }

  do {

    // Set alpha and beta
    m_matrix->setScalars(1.0, dt);

    for(int lev = 0; lev<m_nlev; lev++) {

      A[lev]->setVal(0.);

      for (MFIter mfi(*a_epf[lev]); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*a_epf[lev])[mfi]);
        const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

        if (bx.ok()) {

          Array4<Real const> dummy_arr;

          Array4<Real      > const& A_array       = A[lev]->array(mfi);
          Array4<Real const> const& epf    = a_epf[lev]->const_array(mfi);
          Array4<Real const> const& rho    = a_rho[lev]->const_array(mfi);
          Array4<Real const> const& Tf     = a_Tf[lev]->const_array(mfi);
          Array4<Real const> const& X_gk_array    = fluid_is_a_mixture ?
            X_gk[lev]->const_array(mfi) : dummy_arr;

          auto const& flags_arr = flags.const_array();

          amrex::ParallelFor(bx, [epf,Tf,rho,fluid_props,
              A_array,X_gk_array,dt,flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            if (!cell_is_covered) {

              Real cp_loc = fluid_props.specificHeat(IntVect(i,j,k), Tf, X_gk_array);

              A_array(i,j,k) = epf(i,j,k)*rho(i,j,k)*cp_loc;
            }
          });
        }
      }

      // Compute thermal diffusion coefficient: epf*kappa
      setDiffCoeff(lev, a_epf, GetVecOfConstPtrs(a_Tf));

      // Set diffusion coefficient at MI to zero to prevent additional
      // diffusion of heat into the domain from the mass inflow.
      avgDiffCoeffToFaces(lev);

      if (m_set_eb_dirichlet) {
        m_matrix->setEBDirichlet(lev, *a_Teb[lev], *m_b_eb[lev]);
      }

      m_matrix->setACoeffs(lev, *A[lev]);
      m_matrix->setBCoeffs(lev, GetArrOfConstPtrs(m_b[lev]),
          MLMG::Location::FaceCentroid);

      m_phi[lev]->setVal(0.);
      m_rhs[lev]->setVal(0.);

      for (MFIter mfi(*a_epf[lev]); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.tilebox();

        const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*a_epf[lev])[mfi]);
        const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

        if (bx.ok())
        {
          Array4<Real const> dummy_arr;

          Array4<Real      > const& rhs    = m_rhs[lev]->array(mfi);
          Array4<Real const> const& epf    = a_epf[lev]->const_array(mfi);
          Array4<Real const> const& rho    = a_rho[lev]->const_array(mfi);
          Array4<Real const> const& Tf     = a_Tf[lev]->const_array(mfi);
          Array4<Real const> const& h_f     = a_hf[lev]->const_array(mfi);
          Array4<Real const> const& X_gk_array    = fluid_is_a_mixture ?
            X_gk[lev]->const_array(mfi) : dummy_arr;

          auto const& flags_arr = flags.const_array();

          amrex::ParallelFor(bx, [epf,Tf,rho,fluid_props,
              X_gk_array,h_f,dt,rhs,flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            if (!cell_is_covered) {
              const Real Tf_loc = Tf(i,j,k);
              const Real epf_rho = epf(i,j,k)*rho(i,j,k);

              Real cp_loc = fluid_props.specificHeat(IntVect(i,j,k), Tf, X_gk_array);
              Real hf_loc = fluid_props.enthalpy(IntVect(i,j,k), Tf, X_gk_array);

              rhs(i,j,k) = epf_rho*(cp_loc*Tf_loc - hf_loc + h_f(i,j,k));
            }
          });
        }
      }

      MultiFab::Copy(*m_phi[lev], *a_Tf[lev], 0, 0, 1, 1);
      m_matrix->setLevelBC(lev, GetVecOfConstPtrs(m_phi)[lev]);
    } // end of loop on lev

    MLMG solver(*m_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned
    // from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs),
        m_mg_rtol, m_mg_atol);

    for(int lev = 0; lev<m_nlev; lev++) {
      m_phi[lev]->FillBoundary(m_geom[lev].periodicity());
    }

    for(int lev = 0; lev<m_nlev; lev++) {

      for (MFIter mfi(*a_epf[lev]); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*a_epf[lev])[mfi]);
        const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

        if (bx.ok())
        {
          Array4<Real      > const& update_array = update[lev]->array(mfi);
          Array4<Real      > const& Tf    = a_Tf[lev]->array(mfi);
          Array4<Real const> const& phi   = m_phi[lev]->const_array(mfi);

          auto const& flags_arr = flags.const_array();

          amrex::ParallelFor(bx, [update_array,Tf,phi,flags_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

            if (!cell_is_covered) {
              update_array(i,j,k) = phi(i,j,k) - Tf(i,j,k);
              Tf(i,j,k) = phi(i,j,k);
            }
          });
        }
      }
    } // end of loop on lev

    iter++;
    if (iter > maxiter) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "MFIXDiffOpEnergy::solve !!! Newton solver failed to converge!!!\n"
          << "Newton solver iterations = " << iter << '\n'
          << "Newton solver update = " << norm0(update);
    }

  } while (//(norm0(residue) > residue_rel_tol) ||
           (norm0(update) > update_rel_tol));

  for (int lev(0); lev<m_nlev; ++lev) {
    delete A[lev];
    delete update[lev];
//    delete residue[lev];
  }


  // **************************************************************************
  // **************************************************************************
  // **************************************************************************

  for(int lev = 0; lev<m_nlev; lev++) {

    for (MFIter mfi(*a_epf[lev]); mfi.isValid(); ++mfi) {

      const EBFArrayBox& epg_fab = static_cast<EBFArrayBox const&>((*a_epf[lev])[mfi]);
      const EBCellFlagFab& flags = epg_fab.getEBCellFlagFab();

      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      Array4<Real const> dummy_arr;

      Array4<Real      > const& h_f  = a_hf[lev]->array(mfi);
      Array4<Real const> const& Tf  = a_Tf[lev]->const_array(mfi);
      Array4<Real const> const& X_gk_array = fluid_is_a_mixture ?
        X_gk[lev]->const_array(mfi) : dummy_arr;

      auto const& flags_arr = flags.const_array();

      amrex::ParallelFor(bx, [h_f,Tf,X_gk_array,fluid_props,flags_arr]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

        Real hg = fluid_props.enthalpy(IntVect(i,j,k), Tf, X_gk_array, cell_is_covered);

        h_f(i,j,k) = hg;
      });
    }
  }
}


void MFIXDiffOpEnergy::
setDiffCoeff ( int const a_lev,
               Vector< MultiFab const*> const& a_epf,
               Vector< MultiFab const*> const& a_Tf )
{
  const auto& fluid_parms = m_fluid.parameters<run_on>();

  Box domain(m_geom[a_lev].Domain());

  const auto dlo = lbound(domain);
  const auto dhi = ubound(domain);

  const auto& bc_list = m_bcs.get_bc_list();

  m_b_cc[a_lev]->setVal(0.);

  int const nghost( m_b_cc[a_lev]->nGrow() );

  for (MFIter mfi(*a_epf[a_lev]); mfi.isValid(); ++mfi) {

    Box const& bx = mfi.growntilebox(nghost);

    Array4<Real      > const& b_cc = m_b_cc[a_lev]->array(mfi);
    Array4<Real const> const& epf  = a_epf[a_lev]->const_array(mfi);
    Array4<Real const> const& Tf   = a_Tf[a_lev]->const_array(mfi);

    ParallelFor(bx, [b_cc,epf,Tf,fluid_parms]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      b_cc(i,j,k) = epf(i,j,k)*fluid_parms.calc_k_g(Tf(i,j,k));
    });

    // Set boundary diffusion coefficients
    { Array4<int> const& bct_lo = bc_list.bc_ilo[a_lev]->array();
      Array4<int> const& bct_hi = bc_list.bc_ihi[a_lev]->array();

      ParallelFor(bx, [b_cc, epf, Tf, fluid_parms, dlo, dhi, bct_lo, bct_hi]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (i == dlo.x) {
          if (bct_lo(dlo.x-1,j,k,0) == BCList::pout) {
            b_cc(i-1,j,k) = epf(i,j,k)*fluid_parms.calc_k_g(Tf(i,j,k));
          } else {
            b_cc(i-1,j,k) = 0.;
          }
        }
        if (i == dhi.x) {
          if (bct_hi(dhi.x+1,j,k,0) == BCList::pout) {
            b_cc(i+1,j,k) = epf(i,j,k)*fluid_parms.calc_k_g(Tf(i,j,k));
          } else {
            b_cc(i+1,j,k) = 0.;
          }
        }
      });
    } // x-dir

    { Array4<int> const& bct_lo = bc_list.bc_jlo[a_lev]->array();
      Array4<int> const& bct_hi = bc_list.bc_jhi[a_lev]->array();

      ParallelFor(bx, [b_cc, epf, Tf, fluid_parms, dlo, dhi, bct_lo, bct_hi]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (j == dlo.y) {
          if (bct_lo(i,dlo.y-1,k,0) == BCList::pout) {
            b_cc(i,j-1,k) = epf(i,j,k)*fluid_parms.calc_k_g(Tf(i,j,k));
          } else {
            b_cc(i,j-1,k) = 0.;
          }
        }
        if (j == dhi.y) {
          if (bct_hi(i,dhi.y+1,k,0) == BCList::pout) {
            b_cc(i,j+1,k) = epf(i,j,k)*fluid_parms.calc_k_g(Tf(i,j,k));
          } else {
            b_cc(i,j+1,k) = 0.;
          }
        }
      });
    } // y-dir

    { Array4<int> const& bct_lo = bc_list.bc_klo[a_lev]->array();
      Array4<int> const& bct_hi = bc_list.bc_khi[a_lev]->array();

      ParallelFor(bx, [b_cc, epf, Tf, fluid_parms, dlo, dhi, bct_lo, bct_hi]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (k == dlo.z) {
          if (bct_lo(i,j,dlo.z-1,0) == BCList::pout) {
            b_cc(i,j,k-1) = epf(i,j,k)*fluid_parms.calc_k_g(Tf(i,j,k));
          } else {
            b_cc(i,j,k-1) = 0.;
          }
        }
        if (k == dhi.z) {
          if (bct_hi(i,j,dhi.z+1,0) == BCList::pout) {
            b_cc(i,j,k+1) = epf(i,j,k)*fluid_parms.calc_k_g(Tf(i,j,k));
          } else {
            b_cc(i,j,k+1) = 0.;
          }
        }
      });
    } // z-dir

  } // MFIter

  m_b_cc[a_lev]->FillBoundary(m_geom[a_lev].periodicity());
}



void MFIXDiffOpEnergy::
setEBDiffCoeff ( int const a_lev,
                 Vector< MultiFab const*> const& a_Teb )
{
  const auto& fluid_parms = m_fluid.parameters<run_on>();

  for (MFIter mfi(*a_Teb[a_lev]); mfi.isValid(); ++mfi) {

    Box const& bx = mfi.growntilebox(IntVect(1,1,1));

    Array4<Real      > const& eb_k = m_b_eb[a_lev]->array(mfi);
    Array4<Real const> const& Teb  = a_Teb[a_lev]->const_array(mfi);

    amrex::ParallelFor(bx, [eb_k,Teb,fluid_parms]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      if (Teb(i,j,k) > 0.) {
        eb_k(i,j,k) = fluid_parms.calc_k_g(Teb(i,j,k));
      }
    });
  }
}
