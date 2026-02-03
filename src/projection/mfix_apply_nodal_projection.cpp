#include <mfix.H>
#include <hydro_NodalProjector.H>

using namespace amrex;


/***************************************************************************
 *  Apply the nodal projection                                             *
 *                                                                         *
 *                ∇·(ε/ρ)∇φ =  ∇·(εu^*) - ∇·(εu)                           *
 *                                                                         *
 *   ε := volume fraction      ( )^* denotes explicit update velocity      *
 *   ρ := density                                                          *
 *   u := velocity                                                         *
 *   φ := unknown                                                          *
 *                                                                         *
 * The AMReX MLNodeLaplacian linear solver class uses the canonical form:  *
 *                                                                         *
 *  ∇·σ∇φ = RHS    where    σ := ε/ρ   and   RHS = -(∇·(εu) - ∇·(εu)^*)    *
 *                                                                         *
 *  We add the pressure gradient back to velocity before the projection    *
 *                                                                         *
 *    u^** = u^* +  (Δt/ρ) ∇p                                              *
 *                                                                         *
 *    /---------------------------------------------------------------\    *
 *    | If we take the implicit update for drag:                      |    *
 *    |                                                               |    *
 *    |   u^** = u^* + Δt * ( ε^(n+1/2) ∇p ) / ( (ερ)^(n+1) + Δt*fp ) |    *
 *    |                                                               |    *
 *    | where fp is the drag coefficient.                             |    *
 *    \---------------------------------------------------------------/    *
 *                                                                         *
 *  Update velocity and perturbational pressure:                           *
 *                                                                         *
 *     u^(n+1) = u^** - (1/ρ)∇φ                                            *
 *                                                                         *
 *     pert_p = (1/Δt) ∇φ                                                  *
 *                                                                         *
 ***************************************************************************/
void mfix::
apply_nodal_projection ( Real a_time,
                         Real a_dt,
                         Real a_prev_dt,
                         bool proj_2,
                         Vector<MultiFab      *> const& a_rhs,
                         Vector<MultiFab      *> const& a_vel_old,
                         Vector<MultiFab      *> const& a_vel_new,
                         Vector<MultiFab      *> const& a_pert_p,
                         Vector<MultiFab      *> const& a_grad_p,
                         Vector<MultiFab const*> const& a_epf,
                         Vector<MultiFab const*> const& a_rho,
                         Vector<MultiFab const*> const& a_eb_vel)
{
  BL_PROFILE("mfix::apply_nodal_projection");

  // If we have dropped the dt substantially for whatever reason, use a different form of the approximate
  // projection that projects (U^*-U^n + dt grad_p) rather than (U^* + dt grad_p)
  bool proj_for_small_dt = (a_time > 0. && a_dt < 0.1 * a_prev_dt);

  // Create sigma
  Vector<MultiFab> sigma_mf(nlev());

  for (int lev(0); lev < nlev(); ++lev) {

    sigma_mf[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), *(m_eb->factory()[lev]));

    for (MFIter mfi(*a_vel_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // Tilebox
      Box bx = mfi.tilebox();

      Array4<Real      > const&  vel   = a_vel_new[lev]->array(mfi);
      Array4<Real      > const&  sigma = sigma_mf[lev].array(mfi);

      Array4<Real const> const&  rho = a_rho[lev]->const_array(mfi);
      Array4<Real const> const&  epf = a_epf[lev]->const_array(mfi);
      Array4<Real const> const&  grad_p  = a_grad_p[lev]->const_array(mfi);

      amrex::ParallelFor(bx,[a_dt, proj_2, vel, sigma, rho, epf, grad_p]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        sigma(i,j,k) = 1.0/rho(i,j,k);

        if (proj_2) {
          vel(i,j,k,0) += sigma(i,j,k)*a_dt*grad_p(i,j,k,0);
          vel(i,j,k,1) += sigma(i,j,k)*a_dt*grad_p(i,j,k,1);
          vel(i,j,k,2) += sigma(i,j,k)*a_dt*grad_p(i,j,k,2);
        }

        sigma(i,j,k) *= epf(i,j,k);

      });
    }

    // Print level infos
    if (proj_for_small_dt) {
      Print() << "Before projection (with small dt modification):\n";
    } else {
      Print() << "Before projection:\n";
    }
    m_rw->mfix_print_max_vel(lev, a_vel_new, a_pert_p);
    m_rw->mfix_print_max_gp(lev, a_grad_p);
    Print() << "Min and Max of epf " << a_epf[lev]->min(0) << " " << a_epf[lev]->max(0) << '\n';

    // Set velocities BC before projection
    bcs().set_velocity_bcs(lev, a_time, a_vel_new[lev]);

  } // lev

  // Define "vel" to be U^* - U^n rather than U^*
  if (proj_for_small_dt || (!proj_2)) {

    for(int lev = 0; lev < nlev(); lev++) {

      bcs().set_velocity_bcs(lev, a_time, a_vel_old[lev]);

      MultiFab::Saxpy(*a_vel_new[lev], -1.0, *a_vel_old[lev],
          0, 0, 3, a_vel_new[lev]->nGrow());
    }
  }

  //
  // Compute epu = epf*vel
  //
  Vector< std::unique_ptr<MultiFab>> epu(nlev());

  for (int lev(0); lev < nlev(); ++lev) {
    // We only need one ghost cell here -- so no need to make it bigger
    epu[lev].reset( new MultiFab(grids[lev], dmap[lev], 3, 1,
        MFInfo(), *(m_eb->factory()[lev])) );

    epu[lev]->setVal(1.e200);

    MultiFab::Copy(*epu[lev], *a_vel_new[lev], 0, 0, 3, epu[lev]->nGrow());

    for (int n(0); n < 3; n++) {
      MultiFab::Multiply(*epu[lev], *a_epf[lev], 0, n, 1, epu[lev]->nGrow());
    }
    epu[lev]->setBndry(0.0);
  }

  // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
  //  no-slip walls are treated exactly like slip walls --
  // Note that this routine is essential to impose the correct inflow bc's on
  //  the product ep  * vel
  //
  //  TODO
  // Why are we using this instead of simply multiplying vel and ep with their BCs in place
  // already?
  //
  // For incremental dt and initial iterations, keep boundary values as zero
  if (!proj_for_small_dt && proj_2) {
    for (int lev(0); lev < nlev(); ++lev) {
      bcs().set_vec_bcs(lev, a_time, epu[lev].get());
    }
  }

  for (int lev(0); lev < nlev(); ++lev) {
    // We set these to zero because if the values in the covered cells are undefined,
    //   even though they are multiplied by zero in the divu computation, we can still get NaNs
    EB_set_covered(*epu[lev].get(), 0, epu[lev]->nComp(), 1, 0.0);
  }

  for (int lev(0); lev < nlev(); lev++) {
    EB_set_covered(*a_rhs[lev], 0, a_rhs[lev]->nComp(), 1, 0.0);
  }

  if ( m_redistribute_before_nodal_proj ) {
    PreProjectionRedistribution(a_time);
  }

  //
  // Setup the nodal projector
  //

  LPInfo info;
  info.setMaxCoarseningLevel(nodalproj_options->max_coarsening_level);
  info.setAgglomerationGridSize(agg_grid_size);

  nodal_projector = std::make_unique<Hydro::NodalProjector>(a_vel_new,
      GetVecOfConstPtrs(sigma_mf), Geom(0,nlev()-1), info);

  nodalproj_options->apply(*nodal_projector);
  nodal_projector->setDomainBC(bcs().ppe_lobc(), bcs().ppe_hibc());

  // By setting alpha = epf, the nodal projection will correct the velocity by
  // (sigma / alpha) grad(phi) rather than sigma grad phi
  nodal_projector->setAlpha(GetVecOfConstPtrs(a_epf));

  if (m_embedded_boundaries.has_flow()) {
    for (int lev(0); lev < nlev(); ++lev) {
      nodal_projector->getLinOp().setEBInflowVelocity(lev, *a_eb_vel[lev]);
    }
  }

  Vector< std::unique_ptr<MultiFab>> diveu(nlev());
  for (int lev(0); lev < nlev(); ++lev) {
    diveu[lev].reset( new MultiFab( amrex::convert(grids[lev], IntVect{1,1,1}), dmap[lev],
      /*ncomp =*/1, /*nghost =*/1, MFInfo(), *(m_eb->factory()[lev])));
  }

  nodal_projector->computeRHS(GetVecOfPtrs(diveu), GetVecOfPtrs(epu), a_rhs);
  nodal_projector->setCustomRHS(GetVecOfConstPtrs(diveu));

  nodal_projector->project(nodalproj_options->mg_rtol, nodalproj_options->mg_atol);

  // Define "vel" to be U^{n+1} rather than (U^{n+1}-U^n)
  if (proj_for_small_dt || (!proj_2)) {
    for(int lev = 0; lev < nlev(); lev++) {
      MultiFab::Saxpy(*a_vel_new[lev], 1.0, *a_vel_old[lev], 0, 0, 3, a_vel_new[lev]->nGrow());
    }
  }

  // Get phi and fluxes
  Vector< const MultiFab* > phi(nlev());
  Vector< const MultiFab* > gradphi(nlev());

  phi     = nodal_projector->getPhiConst();
  gradphi = nodal_projector->getGradPhiConst();

  // Since I did not pass dt, I have to normalize here
  Real qdt(1.0/a_dt);
  for (int lev(0); lev < nlev(); ++lev) {
    if (proj_2) {
      // p := phi/dt
      MultiFab::Copy(*a_pert_p[lev], *phi[lev], 0, 0, 1, phi[lev]->nGrow());
      MultiFab::Copy(*a_grad_p[lev], *gradphi[lev], 0, 0, 3, gradphi[lev]->nGrow());
      a_pert_p[lev]->mult(qdt);
      a_grad_p[lev]->mult(qdt);

    } else {
      // p := p + phi/dt
      MultiFab::Saxpy(*a_pert_p[lev], qdt, *phi[lev], 0, 0, 1, phi[lev]->nGrow());
      MultiFab::Saxpy(*a_grad_p[lev], qdt, *gradphi[lev], 0, 0, 3, gradphi[lev]->nGrow());
    }
  }

  for (int lev = finest_level-1; lev >= 0; --lev) {
    amrex::EB_average_down( *a_grad_p[lev+1], *a_grad_p[lev], 0,
        AMREX_SPACEDIM, refRatio(lev));
  }

  // Perform the redistribution operation on the updated (projected) velocity
  // field and update grad_p to maintain consistency
  if ( m_redistribute_nodal_proj ) {
    PostProjectionRedistribution(a_time, a_dt, GetVecOfPtrs(sigma_mf));
  }

  PostProjectionDiagnostics(a_time, a_vel_new, a_pert_p, a_grad_p, proj_for_small_dt);
}

void mfix::
PostProjectionDiagnostics ( Real a_time,
                            Vector<MultiFab      *> const& a_vel_new,
                            Vector<MultiFab      *> const& a_pert_p,
                            Vector<MultiFab      *> const& a_grad_p,
                            bool proj_for_small_dt)
{
  // Extrapolate Dirichlet values to ghost cells -- but do it differently in that
  //  no-slip walls are treated exactly like slip walls --
  // Note that this routine is essential to impose the correct inflow bc's on
  //  the product ep  * vel

  for (int lev = nlev()-1; lev > 0; lev--) {
    avgDown(lev-1, *a_vel_new[lev], *a_vel_new[lev-1]);
    avgDown(lev-1, *a_grad_p[lev], *a_grad_p[lev-1]);
  }

  // Print level info after projection
  for (int lev(0); lev < nlev(); lev++) {

    // Swap ghost cells and apply BCs to velocity
    bcs().set_velocity_bcs(lev, a_time, a_vel_new[lev]);

    if (proj_for_small_dt) {
      Print() << "After  projection (with small dt modification):\n";
    } else {
      Print() << "After  projection:\n";
    }
    m_rw->mfix_print_max_vel(lev, a_vel_new, a_pert_p);
  }
}
