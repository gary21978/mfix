#include <AMReX_MultiFabUtil.H>
#include <AMReX_EB_Redistribution.H>

#include <mfix_diffusion_op.H>
#include <mfix_eb_parms.H>
#include <mfix_run_on.H>
#include <mfix_solids.H>
#include <mfix_diffop_velocity_K.H>

using namespace amrex;

MFIXDiffOpVelocity::
MFIXDiffOpVelocity ( int const a_nlev,
                     Vector<Geometry> const& a_geom,
                     Vector<BoxArray> const& a_grids,
                     Vector<DistributionMapping> const& a_dmap,
                     Vector< const EBFArrayBoxFactory* >const& a_ebfactory,
                     MFIXFluidPhase& a_fluid,
                     MFIXAvgParticleParms& a_avg_pc_parms,
                     MFIXBoundaryConditions& a_bcs,
                     Vector<BCRec> const& a_bcrec)
  : MFIXDiffusionOp(a_nlev, a_geom, AMREX_SPACEDIM, a_fluid, a_bcs, a_bcrec)
  , m_avg_pc_parms(a_avg_pc_parms)
{
  readParameters();

  int const matrix_ncomp(1);
  define( a_grids, a_dmap, a_ebfactory, matrix_ncomp );

  LPInfo info;

  info.setMaxCoarseningLevel(m_mg_max_coarsening_level);
  info.setAgglomerationGridSize(m_mg_agg_grid_size);

  m_matrix.reset( new MLEBTensorOp(m_geom, a_grids, a_dmap, info, a_ebfactory));

  // It is essential that we set MaxOrder to 2 if we want to use the standard
  // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
  // The solver's default order is 3 and this uses three points for the gradient.
  m_matrix->setMaxOrder(2);

  // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
  m_matrix->setDomainBC(a_bcs.diff_vel_lobc(), a_bcs.diff_vel_hibc());

  // Effective viscosity stored in m_b_cc
  AMREX_ASSERT(m_b_defined == 0);
  m_b_cc.resize(a_nlev);
  for (int lev(0); lev < a_nlev; ++lev) {
    m_b_cc[lev].reset( new MultiFab(a_grids[lev], a_dmap[lev], 1, 1,
        MFInfo(), *a_ebfactory[lev]) );
  }
  m_b_defined = 1;

}


//
// Explicit evaluation of velocity diffusion (div(tau))
//
void
MFIXDiffOpVelocity::
computeDivTau ( Vector< MultiFab      * > const& divtau_out,
                Vector< MultiFab      * > const& a_vel,
                Vector< MultiFab const* > const& a_epf,
                Vector< MultiFab const* > const& a_rho,
                Vector< MultiFab const* > const& a_Tf,
                Vector< MultiFab const* > const& a_avg_particle_data,
                Vector< MultiFab const* > const& a_Xfk,
                Vector< MultiFab const* > const& a_eb_inflow,
                int const include_eddy_viscosity/*= true*/)
{
  BL_PROFILE("MFIXDiffOpVelocity::ComputeDivTau");
  AMREX_ASSERT( m_nlev  == a_vel.size() );

  Vector< MultiFab* > divtau_aux(m_nlev);
  Vector< MultiFab> velocity(m_nlev);

  for(int lev(0); lev<m_nlev; ++lev) {

    const BoxArray            ba = a_epf[lev]->boxArray();
    const DistributionMapping dm = a_epf[lev]->DistributionMap();
    const FabFactory<FArrayBox>& ebfac = a_epf[lev]->Factory();

    divtau_aux[lev] = new MultiFab(ba, dm, divtau_out[lev]->nComp(),
        a_epf[lev]->nGrow(), MFInfo(), ebfac);

    divtau_aux[lev]->setVal(0.0);

    velocity[lev].define(ba, dm, AMREX_SPACEDIM, 1, MFInfo(), ebfac);
    MultiFab::Copy(velocity[lev], *a_vel[lev], 0, 0, AMREX_SPACEDIM, 1);

    EB_set_covered(velocity[lev], 0, AMREX_SPACEDIM, 1, 0.);
  }

  setDiffCoeff( GetVecOfConstPtrs(a_vel), a_epf, a_rho, a_Tf,
      a_avg_particle_data, a_Xfk, include_eddy_viscosity);

  // We want to return div (eff_visc grad)) phi
  m_matrix->setScalars(0.0, -1.0);

  // Compute the coefficients
  for (int lev(0); lev<m_nlev; ++lev) {

    Vector<BCRec> dummy_bc(AMREX_SPACEDIM);

    // average_cellcenter_to_face( GetArrOfPtrs(b[lev]), *eta_in[lev], geom[lev] );
    EB_interp_CellCentroid_to_FaceCentroid( *m_b_cc[lev],
        GetArrOfPtrs(m_b[lev]), 0, 0, 1, m_geom[lev], /*m_bcrec*/ dummy_bc);

    m_matrix->setShearViscosity(lev, GetArrOfConstPtrs(m_b[lev]),
        MLMG::Location::FaceCentroid);

    if (m_set_eb_dirichlet) {

      m_matrix->setEBShearViscosityWithInflow(lev,
          *m_b_cc[lev], *a_eb_inflow[lev]);

    } else {

      m_matrix->setEBShearViscosity(lev, *m_b_cc[lev]);

    }

    m_matrix->setLevelBC(lev, &velocity[lev]);
  }

  MLMG solver(*m_matrix);

  solver.apply(divtau_aux, GetVecOfPtrs(velocity));

  for(int lev(0); lev<m_nlev; ++lev) {

    single_level_weighted_redistribute(*divtau_aux[lev], *divtau_out[lev],
        *a_epf[lev], 0, AMREX_SPACEDIM, m_geom[lev], true);

    EB_set_covered(*divtau_out[lev], 0, divtau_out[lev]->nComp(), divtau_out[lev]->nGrow(), 0.);
  }

  for(int lev(0); lev<m_nlev; ++lev) {
    if (divtau_aux[lev] != nullptr) { delete divtau_aux[lev]; }
  }
}


//
// Implicit tensor solve for velocity diffusion
//
/***************************************************************************
 * Implicit tensor solve for momentum diffusion:                           *
 *                                                                         *
 *               (ερ + dt*Sp)u - dt[ε∇·τ] = ερu* + dt*Sc                   *
 *                                                                         *
 * where the viscous stress is given by                                    *
 *                                                                         *
 *   τ = μ [∇u + (∇u)^T] + λ(∇·u)I                                         *
 *                                                                         *
 *   ε := volume fraction @ t^(n+1/2)                                      *
 *   ρ := density @ t^(n+1/2)                                              *
 *   μ := diffusion coefficient (effective viscoity)                       *
 *   u := velocity @ t^(n+1)  <--  (unknown)                               *
 *                                                                         *
 * and dt(Sc - uSp) is the general linearized source term.                 *
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
void
MFIXDiffOpVelocity::
solve ( Real const a_dt, int const a_drag_includes_divtau,
        Vector< MultiFab      *> const& a_vel,
        Vector< MultiFab const*> const& a_epf,
        Vector< MultiFab const*> const& a_rho,
        Vector< MultiFab const*> const& a_Tf,
        Vector< MultiFab const*> const& a_avg_particle_data,
        Vector< MultiFab const*> const& a_Xfk,
        Real const a_drag_dt,
        Vector< MultiFab const* >const& S_p_in,
        Vector< MultiFab const* >const& S_c_in,
        Vector< MultiFab const*> const& a_eb_inflow)
{
  BL_PROFILE("MFIXDiffOpVelocity::solve");
  AMREX_ASSERT( m_nlev  == a_vel.size() );

  Vector< MultiFab*> Acoeff(m_nlev);

  for(int lev(0); lev<m_nlev; ++lev) {

    const BoxArray                  ba = a_epf[lev]->boxArray();
    const DistributionMapping       dm = a_epf[lev]->DistributionMap();
    const FabFactory<FArrayBox>& ebfac = a_epf[lev]->Factory();

    Acoeff[lev] = new MultiFab(ba, dm, 1, 1, MFInfo(), ebfac);

    // Acoeff = epg*rho
    MultiFab::Copy(*Acoeff[lev], *a_rho[lev], 0, 0, 1, 1);
    MultiFab::Multiply(*Acoeff[lev], *a_epf[lev], 0, 0, 1, 1);

    // rhs = epg*rho*vel
    MultiFab::Copy((*m_rhs[lev]),(*a_vel[lev]), 0, 0, AMREX_SPACEDIM, 0);
    for (int idim(0); idim < AMREX_SPACEDIM; ++idim)
    { MultiFab::Multiply((*m_rhs[lev]), (*Acoeff[lev]), 0, idim, 1, 0); }

    // Include drag terms if needed.
    if (a_drag_dt > 0.) {

      // Acoeff = epg*rho + dt*Sp
      MultiFab::Saxpy(*Acoeff[lev], a_drag_dt, *S_p_in[lev], 0, 0, 1, 1);

      // rhs = epg*rho*vel + dt*Sc
      MultiFab::Saxpy(*m_rhs[lev], a_drag_dt, *S_c_in[lev], 0, 0, 3, 0);
    }

    // Remove epg from Acoeff and rhs if drag includes shear stress
    if(a_drag_includes_divtau) {

      // Acoeff = (epg*rho + dt*Sp) / epg
      MultiFab::Divide(*Acoeff[lev], *a_epf[lev], 0, 0, 1, 1);

      // rhs = (epg*rho*vel + dt*Sc) / epg
      for (int idim(0); idim < AMREX_SPACEDIM; ++idim)
      { MultiFab::Divide(*m_rhs[lev], *a_epf[lev], 0, idim, 1, 0); }
    }
  }

  setDiffCoeff(GetVecOfConstPtrs(a_vel), a_epf, a_rho, a_Tf,
      a_avg_particle_data, a_Xfk, /*include_eddy_viscosity=*/1);

  // Update the coefficients of the matrix going into the solve based on the
  // current state of the simulation. Recall that the relevant matrix is
  //
  //      alpha A - beta div ( b grad )   <--->   rho - dt div ( mu grad )
  //
  //      alpha: 1
  //      beta: dt
  //      A: epg*rho or rho if drag includes div(tau)
  //      b: eff_visc

  // Set alpha and beta
  m_matrix->setScalars(1.0, a_dt);

  for(int lev(0); lev<m_nlev; ++lev) {

    Vector<BCRec> dummy_bc(AMREX_SPACEDIM);

    // Compute the spatially varying b coefficients (on faces) to equal the apparent viscosity
    EB_interp_CellCentroid_to_FaceCentroid( *m_b_cc[lev],
      GetArrOfPtrs(m_b[lev]), 0, 0, 1, m_geom[lev], /*m_bcrec*/ dummy_bc);

    // This sets the coefficients
    m_matrix->setACoeffs(lev, *Acoeff[lev]);
    m_matrix->setShearViscosity(lev, GetArrOfConstPtrs(m_b[lev]),
        MLMG::Location::FaceCentroid);

    m_matrix->setEBShearViscosity(lev, *m_b_cc[lev]);

    if (m_set_eb_dirichlet) {

      m_matrix->setEBShearViscosityWithInflow(lev,
          *m_b_cc[lev], *a_eb_inflow[lev]);

    }
  }

  if(m_verbose > 0) { Print() << "Diffusing velocity components all together\n"; }

  for(int lev(0); lev<m_nlev; ++lev) {

    // By this point we must have filled the Dirichlet values of phi stored in ghost cells
    MultiFab::Copy(*m_phi[lev],*a_vel[lev], 0, 0, AMREX_SPACEDIM, 1);
    m_matrix->setLevelBC(lev, GetVecOfConstPtrs(m_phi)[lev]);

  }

  MLMG solver(*m_matrix);
  setSolverSettings(solver);

  // This ensures that ghost cells of phi are filled by the solver
  solver.setFinalFillBC(true);

  solver.solve(GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs),
      m_mg_rtol, m_mg_atol);

  for(int lev(0); lev<m_nlev; ++lev) {
    MultiFab::Copy(*a_vel[lev], *m_phi[lev], 0, 0, AMREX_SPACEDIM, 1);
  }

  for (int lev(0); lev<m_nlev; ++lev) {
    if (Acoeff[lev] != nullptr ) { delete Acoeff[lev]; }
  }

  if(m_verbose > 0) { Print() << " Done diffusing all velocity components\n"; }
}

void
MFIXDiffOpVelocity::
setDiffCoeff ( Vector<MultiFab const* >const& a_vel,
               Vector<MultiFab const* >const& a_epf,
               Vector<MultiFab const* >const& a_rho,
               Vector<MultiFab const* >const& a_Tf,
               Vector<MultiFab const* >const& a_avg_particle_data,
               Vector<MultiFab const* >const& a_Xfk,
               int const a_include_eddy_viscosity)
{
  const auto fluid_props = m_fluid.props.data<run_on>();

  AMREX_ASSERT( m_nlev  == a_vel.size() );
  for (int lev(0); lev < m_nlev; ++lev) {
    m_b_cc[lev]->setVal(0.);

    const FabFactory<FArrayBox>& lev_ebfac = a_epf[lev]->Factory();
    const auto& ebfac = static_cast<EBFArrayBoxFactory const&>(lev_ebfac);

    const RealVect dx(AMREX_D_DECL(m_geom[lev].CellSize()[0],
                                   m_geom[lev].CellSize()[1],
                                   m_geom[lev].CellSize()[2]));

    for (MFIter mfi(*a_vel[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      Box const& bx = mfi.growntilebox(IntVect(1));

      Array4<Real      > const& mf = m_b_cc[lev]->array(mfi);

      Array4<Real const> const& Tf = m_fluid.solve_enthalpy()
                                   ? a_Tf[lev]->const_array(mfi)
                                   : Array4<const Real>();

      Array4<Real const> const& Xfk = m_fluid.solve_species()
                                    ? a_Xfk[lev]->const_array(mfi)
                                    : Array4<const Real>();

      ParallelFor(bx, [fluid_props,mf,Tf,Xfk]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        mf(i,j,k) = fluid_props.molViscosity(i, j, k, Tf, Xfk);
      });

      if ( !m_fluid.eddy_visc_model( EddyViscModel::None) ||
           !m_fluid.susp_visc_model( SuspViscModel::None) )
      {

        EBCellFlagFab const& flags_fab = ebfac.getMultiEBCellFlagFab()[mfi];
        const auto& flags = flags_fab.const_array();
        const auto& vfrac = ebfac.getVolFrac().const_array(mfi);

        Array4<Real const> const& vel = a_vel[lev]->const_array(mfi);
        Array4<Real const> const& rho = a_rho[lev]->const_array(mfi);
        Array4<Real const> const& epf = a_epf[lev]->const_array(mfi);

        int const avg_radius_comp = m_avg_pc_parms.comp(SoArealData::radius);
        int const avg_velx_comp = m_avg_pc_parms.comp(SoArealData::velx);
        int const avg_vely_comp = m_avg_pc_parms.comp(SoArealData::vely);
        int const avg_velz_comp = m_avg_pc_parms.comp(SoArealData::velz);

        if ( m_fluid.susp_visc_model( SuspViscModel::Sato ) ||
             m_fluid.susp_visc_model( SuspViscModel::Subramaniam ) ) {
          AMREX_ASSERT( avg_radius_comp >= 0 );
          AMREX_ASSERT( avg_velx_comp >= 0 );
          AMREX_ASSERT( avg_vely_comp >= 0 );
          AMREX_ASSERT( avg_velz_comp >= 0 );
        }

        Array4<Real const> const& avg_radius =
          (avg_radius_comp < 0) ? Array4<const Real>()
          : a_avg_particle_data[lev]->const_array(mfi,avg_radius_comp);

        Array4<Real const> const& avg_velp_x =
          (avg_velx_comp < 0) ? Array4<const Real>()
          : a_avg_particle_data[lev]->const_array(mfi,avg_velx_comp);

        Array4<Real const> const& avg_velp_y =
          (avg_vely_comp < 0) ? Array4<const Real>()
          : a_avg_particle_data[lev]->const_array(mfi, avg_vely_comp);

        Array4<Real const> const& avg_velp_z =
          (avg_velz_comp < 0) ? Array4<const Real>()
          : a_avg_particle_data[lev]->const_array(mfi, avg_velz_comp);

        ParallelFor(bx, [fluid_props,flags,vfrac,mf,vel,rho,epf,Tf,dx,
          avg_radius,avg_velp_x,avg_velp_y,avg_velp_z,
          a_include_eddy_viscosity,
          eddy_visc_model=m_fluid.eddy_visc_model(),susp_visc_model=m_fluid.susp_visc_model(),
          smagorinsky_constant=m_fluid.smagorinsky_constant(),
          wale_constant=m_fluid.wale_constant(),
          brinkman_constant=m_fluid.brinkman_constant(),
          roscoe_c1=m_fluid.roscoe_c1(),roscoe_c2=m_fluid.roscoe_c2(),
          chenglaw_constant=m_fluid.chenglaw_constant(),
          sato_constant = m_fluid.sato_constant(),
          subramaniam_constant = m_fluid.subramaniam_constant(),
          max_effective_viscosity_factor=m_fluid.max_effective_viscosity_factor()]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
           Real mol_visc = mf(i,j,k);
           Real edd_visc = 0.;
           Real susp_visc = 0.;

           if (eddy_visc_model == EddyViscModel::Smagorinsky) {

             edd_visc = calc_edd_visc_smagorinsky(i, j, k, vel, rho,
                 flags, dx, vfrac, smagorinsky_constant);

           } else if (eddy_visc_model == EddyViscModel::WALE) {

             edd_visc = calc_edd_visc_wale(i, j, k, vel, rho,
                 flags, dx, vfrac, wale_constant);
           }

           if (susp_visc_model == SuspViscModel::Einstein ) {

             susp_visc = calc_susp_visc_einstein(i, j, k, epf,
                 mol_visc);

           } else if (susp_visc_model == SuspViscModel::Brinkman ) {

             susp_visc = calc_susp_visc_brinkman(i, j, k, epf,

                 mol_visc, brinkman_constant);
           } else if (susp_visc_model == SuspViscModel::Roscoe ) {

             susp_visc = calc_susp_visc_roscoe(i, j, k, epf,
                 mol_visc, roscoe_c1, roscoe_c2);

           } else if (susp_visc_model == SuspViscModel::ChengLaw ) {

             susp_visc = calc_susp_visc_chenglaw(i, j, k, epf,
                 mol_visc, chenglaw_constant);

           } else if (susp_visc_model == SuspViscModel::Sato ) {

             susp_visc = calc_susp_visc_sato(i, j, k, epf,
                vel, rho, avg_radius, avg_velp_x, avg_velp_y, avg_velp_z,
                sato_constant);

           } else if (susp_visc_model == SuspViscModel::Subramaniam ) {

             susp_visc = calc_susp_visc_subramaniam(i, j, k, epf,
                mol_visc, vel, rho, avg_radius, avg_velp_x, avg_velp_y, avg_velp_z,
                subramaniam_constant);
           }

           Real eff_visc = mol_visc + susp_visc;
           if (a_include_eddy_viscosity) { eff_visc += edd_visc; }

           mf(i,j,k) = min(eff_visc, max_effective_viscosity_factor*mol_visc);
        });
      }
    }

    m_b_cc[lev]->FillBoundary(m_geom[lev].periodicity());
    EB_set_covered(*m_b_cc[lev], 1.e40);
  }
}
