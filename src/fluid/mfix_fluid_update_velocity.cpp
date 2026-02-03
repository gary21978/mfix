#include <mfix_fluid_update.H>

using namespace amrex;

void FluidUpdate::
Velocity ( int const a_advect_momentum,
           int const a_include_src,
           int const a_drag_includes_divtau,
           MFIXDiffOpVelocity* a_diffOp)
{
  BL_PROFILE("FluidUpdate::Velocity");

  int const explicit_diffusion( m_explicit_diffusion );

  int const advect_momentum( a_advect_momentum );
  int const drag_includes_divtau( a_drag_includes_divtau );

  for (int lev(0); lev < nlev(); ++lev) {

    if ( m_step_type == StepType::Predictor ) {

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real const> const& vel_old = leveldata().vel_old_const(lev,mfi);
        Array4<Real      > const& vel_new = leveldata().vel(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);
        Array4<Real const> const& rho_nph = stepdata().rho_nph_const(lev,mfi);

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);
        Array4<Real const> const& epf_nph = stepdata().epf_nph_const(lev,mfi);

        Array4<Real const> const& dudt_old = predictor().conv_velocity_const(lev,mfi);

        Array4<Real const> const& divtau_old = leveldata().div_tau_const(lev,mfi);
        Array4<Real const> const& forces = stepdata().vel_forces_const(lev,mfi);

        if (!a_advect_momentum) {// convective differencing

          ParallelFor(bx, [dt=m_dt, explicit_diffusion, drag_includes_divtau,
          epf_nph, rho_nph, vel_old, vel_new, dudt_old, divtau_old, forces ]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Real vel_nx = epf_nph(i,j,k)*vel_old(i,j,k,0);
            Real vel_ny = epf_nph(i,j,k)*vel_old(i,j,k,1);
            Real vel_nz = epf_nph(i,j,k)*vel_old(i,j,k,2);

            vel_nx += dt*dudt_old(i,j,k,0);
            vel_ny += dt*dudt_old(i,j,k,1);
            vel_nz += dt*dudt_old(i,j,k,2);

            vel_nx /= epf_nph(i,j,k);
            vel_ny /= epf_nph(i,j,k);
            vel_nz /= epf_nph(i,j,k);

            if (explicit_diffusion) {

              amrex::Real const rop_nph = rho_nph(i,j,k) *
                (drag_includes_divtau ? 1.0 : epf_nph(i,j,k));

                vel_nx += dt * divtau_old(i,j,k,0)/rop_nph;
                vel_ny += dt * divtau_old(i,j,k,1)/rop_nph;
                vel_nz += dt * divtau_old(i,j,k,2)/rop_nph;
            }

            // gravity and pressure grad
            vel_new(i,j,k,0) = vel_nx + dt * forces(i,j,k,0);
            vel_new(i,j,k,1) = vel_ny + dt * forces(i,j,k,1);
            vel_new(i,j,k,2) = vel_nz + dt * forces(i,j,k,2);

          });

        } else {// advect momentum

          ParallelFor(bx, [dt=m_dt, explicit_diffusion, drag_includes_divtau,
            epf_old, epf_new, rho_old, rho_new, vel_old, vel_new,
            dudt_old, divtau_old, forces]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            amrex::Real vel_nx = epf_old(i,j,k)*rho_old(i,j,k)*vel_old(i,j,k,0);
            amrex::Real vel_ny = epf_old(i,j,k)*rho_old(i,j,k)*vel_old(i,j,k,1);
            amrex::Real vel_nz = epf_old(i,j,k)*rho_old(i,j,k)*vel_old(i,j,k,2);

            vel_nx += dt*dudt_old(i,j,k,0);
            vel_ny += dt*dudt_old(i,j,k,1);
            vel_nz += dt*dudt_old(i,j,k,2);

            if (explicit_diffusion) {

              amrex::Real const epf( drag_includes_divtau ? epf_new(i,j,k) : 1.0 );

              vel_nx += dt * epf * divtau_old(i,j,k,0);
              vel_ny += dt * epf * divtau_old(i,j,k,1);
              vel_nz += dt * epf * divtau_old(i,j,k,2);
            }

            vel_nx/= (epf_new(i,j,k)*rho_new(i,j,k));
            vel_ny/= (epf_new(i,j,k)*rho_new(i,j,k));
            vel_nz/= (epf_new(i,j,k)*rho_new(i,j,k));

            vel_new(i,j,k,0) = vel_nx + dt * forces(i,j,k,0);
            vel_new(i,j,k,1) = vel_ny + dt * forces(i,j,k,1);
            vel_new(i,j,k,2) = vel_nz + dt * forces(i,j,k,2);

          });
        }
      } // mfi

    } else { AMREX_ASSERT( m_step_type == StepType::Corrector );

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real const> const& vel_old = leveldata().vel_old_const(lev,mfi);
        Array4<Real      > const& vel_new = leveldata().vel(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);
        Array4<Real const> const& rho_nph = stepdata().rho_nph_const(lev,mfi);

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);
        Array4<Real const> const& epf_nph = stepdata().epf_nph_const(lev,mfi);

        Array4<Real const> const& dudt_old = predictor().conv_velocity_const(lev,mfi);
        Array4<Real const> const& dudt_new = corrector().conv_velocity_const(lev,mfi);

        Array4<Real const> const& divtau_old = leveldata().div_tau_const(lev,mfi);
        Array4<Real const> const& forces = stepdata().vel_forces_const(lev,mfi);

        if (!advect_momentum) { // convective differencing

          ParallelFor(bx, [dt=m_dt, drag_includes_divtau, epf_nph, rho_nph,
          vel_old, vel_new, dudt_old, dudt_new, divtau_old, forces]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            amrex::Real vel_nx = epf_nph(i,j,k)*vel_old(i,j,k,0);
            amrex::Real vel_ny = epf_nph(i,j,k)*vel_old(i,j,k,1);
            amrex::Real vel_nz = epf_nph(i,j,k)*vel_old(i,j,k,2);

            vel_nx += 0.5 * dt * (dudt_old(i,j,k,0) + dudt_new(i,j,k,0));
            vel_ny += 0.5 * dt * (dudt_old(i,j,k,1) + dudt_new(i,j,k,1));
            vel_nz += 0.5 * dt * (dudt_old(i,j,k,2) + dudt_new(i,j,k,2));

            vel_nx /= epf_nph(i,j,k);
            vel_ny /= epf_nph(i,j,k);
            vel_nz /= epf_nph(i,j,k);

            amrex::Real const rop_nph = rho_nph(i,j,k) *
              (drag_includes_divtau ? 1.0 : epf_nph(i,j,k));

            // Crank-Nicolson so we only add half of the explicit term.

            vel_nx += dt * 0.5 * divtau_old(i,j,k,0) / rop_nph;
            vel_ny += dt * 0.5 * divtau_old(i,j,k,1) / rop_nph;
            vel_nz += dt * 0.5 * divtau_old(i,j,k,2) / rop_nph;

            vel_new(i,j,k,0) = vel_nx + dt*forces(i,j,k,0);
            vel_new(i,j,k,1) = vel_ny + dt*forces(i,j,k,1);
            vel_new(i,j,k,2) = vel_nz + dt*forces(i,j,k,2);

          });

        } else { // advect_momentum

          amrex::ParallelFor(bx, [dt=m_dt, drag_includes_divtau, rho_old, rho_new,
          epf_old, epf_new, vel_old, vel_new, dudt_old, dudt_new, divtau_old, forces]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            amrex::Real vel_nx = epf_old(i,j,k)*rho_old(i,j,k)*vel_old(i,j,k,0);
            amrex::Real vel_ny = epf_old(i,j,k)*rho_old(i,j,k)*vel_old(i,j,k,1);
            amrex::Real vel_nz = epf_old(i,j,k)*rho_old(i,j,k)*vel_old(i,j,k,2);

            vel_nx += dt * 0.5 * (dudt_old(i,j,k,0) + dudt_new(i,j,k,0));
            vel_ny += dt * 0.5 * (dudt_old(i,j,k,1) + dudt_new(i,j,k,1));
            vel_nz += dt * 0.5 * (dudt_old(i,j,k,2) + dudt_new(i,j,k,2));

            amrex::Real const epf( drag_includes_divtau ? epf_new(i,j,k) : 1.0 );

            // Crank-Nicolson so we only add half of the explicit term.

            vel_nx += 0.5 * dt * epf * divtau_old(i,j,k,0);
            vel_ny += 0.5 * dt * epf * divtau_old(i,j,k,1);
            vel_nz += 0.5 * dt * epf * divtau_old(i,j,k,2);

            vel_nx/= (epf_new(i,j,k)*rho_new(i,j,k));
            vel_ny/= (epf_new(i,j,k)*rho_new(i,j,k));
            vel_nz/= (epf_new(i,j,k)*rho_new(i,j,k));

            vel_new(i,j,k,0) = vel_nx + dt * forces(i,j,k,0);
            vel_new(i,j,k,1) = vel_ny + dt * forces(i,j,k,1);
            vel_new(i,j,k,2) = vel_nz + dt * forces(i,j,k,2);

          });
        }

      } // mfi
    } // step_type

    // ****************************************************************************
    // Add additional sources implicitly:  u^(n+1,*) = Sc - Sp*u^(n+1,*)
    // ****************************************************************************

    if ( explicit_diffusion && a_include_src ) {

      VelocityImplicitUpdate( leveldata().vel(lev),
        (!a_advect_momentum ? stepdata().epf_nph(lev) : leveldata().epf(lev)),
        (!a_advect_momentum ? stepdata().rho_nph(lev) : leveldata().rho(lev)),
        stepdata().vel_S_p_const(lev), stepdata().vel_S_c_const(lev));
    }

    leveldata().vel(lev)->FillBoundary(m_geom[lev].periodicity());

  } // lev

  // TODO: This should only be needed for the diffusive solve.
  for (int lev(0); lev < nlev(); lev++) {
    m_bcs.set_velocity_bcs(lev, m_time + m_dt, leveldata().vel(lev));
  }

  if ( !explicit_diffusion ) {

    // drag_dt > 0 indicates that we include drag in the implicit solve.
    Real const drag_dt( a_include_src ? m_dt : 0.0);
    Real const diff_dt( m_dt*(m_step_type == StepType::Predictor ? 1.0 : 0.5) );

    // Diffuse velocity
    a_diffOp->solve( diff_dt, a_drag_includes_divtau, leveldata().vel(),
      (!a_advect_momentum ? stepdata().epf_nph_const() : leveldata().epf_const()),
      (!a_advect_momentum ? stepdata().rho_nph_const() : leveldata().rho_const()),
      leveldata().T_const(), stepdata().avg_particle_data_const(),
      leveldata().X_const(), drag_dt, stepdata().vel_S_p_const(),
      stepdata().vel_S_c_const(), stepdata().eb_flow_velocity_const());

    // TODO: This shouldn't be needed.
    for (int lev(0); lev < nlev(); lev++) {
      m_bcs.set_velocity_bcs(lev, m_time + m_dt, leveldata().vel(lev));
    }
  }

}

/*****************************************************************
* Incorporate velocity source terms                              *
*                                                                *
*    u^(n+1,*) = Sc - Sp*u^(n+1,*)                               *
*                                                                *
******************************************************************/
void FluidUpdate::
VelocityImplicitUpdate ( MultiFab*       a_vel,
                         MultiFab const* a_epf,
                         MultiFab const* a_rho,
                         MultiFab const* a_S_p,
                         MultiFab const* a_S_c)
{
  BL_PROFILE("FluidUpdate::VelocityImplicitUpdate");

  for (MFIter mfi(*a_vel,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    // Tilebox
    Box bx = mfi.tilebox();

    Array4<Real      > const& vel = a_vel->array(mfi);

    Array4<Real const> const& epf = a_epf->const_array(mfi);
    Array4<Real const> const& rho = a_rho->const_array(mfi);

    Array4<Real const> const& S_p = a_S_p->const_array(mfi);
    Array4<Real const> const& S_c = a_S_c->const_array(mfi);

    amrex::ParallelFor(bx,[dt=m_dt, vel, epf, rho, S_c, S_p]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      Real rop  = epf(i,j,k)*rho(i,j,k);

      for (int idim(0); idim<AMREX_SPACEDIM; ++idim) {
        vel(i,j,k,idim) = (rop*vel(i,j,k,idim) + dt*S_c(i,j,k,idim))
                        / (rop + dt*S_p(i,j,k));
      }
    });
  }
}


/*****************************************************************
* Incorporate velocity source terms                              *
*                                                                *
*    u^(n+1,*) = Sc - Sp*u^(n+1,*)                               *
*                                                                *
******************************************************************/
void FluidUpdate::
ComputeVelTxfrSrc ( AdvectionType const a_adv_type,
                    int const a_advect_momentum,
                    int const a_include_virtual_mass)
{
  BL_PROFILE("FluidUpdate::ComputeVelTxfrSrc");

  InterphaseTxfrIndexes txfr_comps;
  int const drag_comp( txfr_comps.drag_coeff);
  int const vm_comp( txfr_comps.vm_coeff);
  int const Sc_comp( txfr_comps.vel_src_c);

  int const nspecies( m_fluid.nspecies() );

  AMREX_ASSERT( (nspecies <= 1 && !use_species_advection() ) ||
                (nspecies >  1 &&  use_species_advection()  ));

  int const advect_momentum( a_advect_momentum );
  int const include_virtual_mass( a_include_virtual_mass );

  for (int lev(0); lev < nlev(); ++lev) {

    if ( m_step_type == StepType::Predictor ) {

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real      > const& S_p = stepdata().vel_S_p(lev,mfi);
        Array4<Real      > const& S_c = stepdata().vel_S_c(lev,mfi);

        Array4<Real const> const& drag_coeff = leveldata().txfr_const(lev,mfi,drag_comp);
        Array4<Real const> const& vm_coeff   = leveldata().txfr_const(lev,mfi,vm_comp);
        Array4<Real const> const& vel_src    = leveldata().txfr_const(lev,mfi,Sc_comp);

        Array4<Real const> const& epf = ( a_adv_type == AdvectionType::MOL ) ?
            leveldata().epf_old_const(lev,mfi) : stepdata().epf_nph_const(lev,mfi);

        Array4<Real const> const& rho = ( a_adv_type == AdvectionType::MOL ) ?
            leveldata().rho_old_const(lev,mfi) : stepdata().rho_nph_const(lev,mfi);

        Array4<Real const> const& vel_old = leveldata().vel_old_const(lev,mfi);

        Array4<Real const> const& dudt_old = predictor().conv_velocity_const(lev,mfi);

        Array4<Real const> const& drdt_old = !use_species_advection()
            ? predictor().conv_scalars_const(lev,mfi,AdvectionComp::density)
            : predictor().conv_species_const(lev,mfi);

        ParallelFor(bx, [dt=m_dt, nspecies, advect_momentum, include_virtual_mass,
        S_p, S_c, drag_coeff, vm_coeff, vel_src, epf, rho,
        vel_old, dudt_old, drdt_old ]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          S_p(i,j,k) = drag_coeff(i,j,k);

          S_c(i,j,k,0) = vel_src(i,j,k,0);
          S_c(i,j,k,1) = vel_src(i,j,k,1);
          S_c(i,j,k,2) = vel_src(i,j,k,2);

          if (include_virtual_mass) {

            RealVect epfu_divu(dudt_old(i,j,k,0), dudt_old(i,j,k,1), dudt_old(i,j,k,2));

            if (advect_momentum) {

              Real A_rho = drdt_old(i,j,k,0);
              for (int n(1); n < nspecies; ++n)
              {  A_rho += drdt_old(i,j,k,n); }

              epfu_divu[0] = (epfu_divu[0] - vel_old(i,j,k,0)*A_rho) / rho(i,j,k);
              epfu_divu[1] = (epfu_divu[1] - vel_old(i,j,k,1)*A_rho) / rho(i,j,k);
              epfu_divu[2] = (epfu_divu[2] - vel_old(i,j,k,2)*A_rho) / rho(i,j,k);
            }

            S_p(i,j,k) += (vm_coeff(i,j,k) / dt);

            // Cvm*dupdt is included in vel_src deposition so we only need to account
            // for the explicit component from the fluid material derivative.
            S_c(i,j,k,0) += vm_coeff(i,j,k)*(vel_old(i,j,k,0)/dt + epfu_divu[0]/epf(i,j,k) );
            S_c(i,j,k,1) += vm_coeff(i,j,k)*(vel_old(i,j,k,1)/dt + epfu_divu[1]/epf(i,j,k) );
            S_c(i,j,k,2) += vm_coeff(i,j,k)*(vel_old(i,j,k,2)/dt + epfu_divu[2]/epf(i,j,k) );
          }
        });
      }

    } else { AMREX_ASSERT( m_step_type == StepType::Corrector );

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real      > const& S_p = stepdata().vel_S_p(lev,mfi);
        Array4<Real      > const& S_c = stepdata().vel_S_c(lev,mfi);

        Array4<Real const> const& drag_coeff = leveldata().txfr_const(lev,mfi,drag_comp);
        Array4<Real const> const& vm_coeff   = leveldata().txfr_const(lev,mfi,vm_comp);
        Array4<Real const> const& vel_src    = leveldata().txfr_const(lev,mfi,Sc_comp);

        Array4<Real const> const& vel_old = leveldata().vel_old_const(lev,mfi);
        Array4<Real const> const& vel_new = leveldata().vel_const(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

        Array4<Real const> const& dudt_old = predictor().conv_velocity_const(lev,mfi);
        Array4<Real const> const& dudt_new = corrector().conv_velocity_const(lev,mfi);

        Array4<Real const> const& drdt_old = !use_species_advection()
            ? predictor().conv_scalars_const(lev,mfi,AdvectionComp::density)
            : predictor().conv_species_const(lev,mfi);

        Array4<Real const> const& drdt_new = !use_species_advection()
            ? corrector().conv_scalars_const(lev,mfi,AdvectionComp::density)
            : corrector().conv_species_const(lev,mfi);

        ParallelFor(bx, [dt=m_dt, nspecies, advect_momentum, include_virtual_mass,
          S_p, S_c, drag_coeff, vm_coeff, vel_src, epf_old, epf_new,
          rho_old, rho_new,  vel_old, vel_new, dudt_old, dudt_new,
          drdt_old, drdt_new ]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          S_p(i,j,k) = drag_coeff(i,j,k);

          S_c(i,j,k,0) = vel_src(i,j,k,0);
          S_c(i,j,k,1) = vel_src(i,j,k,1);
          S_c(i,j,k,2) = vel_src(i,j,k,2);

          if (include_virtual_mass) {

            RealVect epfu_divu_old(dudt_old(i,j,k,0), dudt_old(i,j,k,1), dudt_old(i,j,k,2));
            RealVect epfu_divu_new(dudt_new(i,j,k,0), dudt_new(i,j,k,1), dudt_new(i,j,k,2));

            if(advect_momentum) {

              Real A_rho_old = drdt_old(i,j,k,0);
              for (int n(1); n < nspecies; ++n)
              {  A_rho_old += drdt_old(i,j,k,n); }

              epfu_divu_old[0] = (epfu_divu_old[0] - vel_old(i,j,k,0)*A_rho_old) / rho_old(i,j,k);
              epfu_divu_old[1] = (epfu_divu_old[1] - vel_old(i,j,k,1)*A_rho_old) / rho_old(i,j,k);
              epfu_divu_old[2] = (epfu_divu_old[2] - vel_old(i,j,k,2)*A_rho_old) / rho_old(i,j,k);

              Real A_rho_new = drdt_new(i,j,k,0);
              for (int n(1); n < nspecies; ++n)
              {  A_rho_new += drdt_new(i,j,k,n); }

              epfu_divu_new[0] = (epfu_divu_new[0] - vel_new(i,j,k,0)*A_rho_new) / rho_new(i,j,k);
              epfu_divu_new[1] = (epfu_divu_new[1] - vel_new(i,j,k,1)*A_rho_new) / rho_new(i,j,k);
              epfu_divu_new[2] = (epfu_divu_new[2] - vel_new(i,j,k,2)*A_rho_new) / rho_new(i,j,k);

            }

            S_p(i,j,k) += vm_coeff(i,j,k);

            // dt*dupdt is included in vel_src deposition
            S_c(i,j,k,0) += vm_coeff(i,j,k)*vel_old(i,j,k,0);
            S_c(i,j,k,1) += vm_coeff(i,j,k)*vel_old(i,j,k,1);
            S_c(i,j,k,2) += vm_coeff(i,j,k)*vel_old(i,j,k,2);

            S_c(i,j,k,0) += 0.5*dt*vm_coeff(i,j,k)*( epfu_divu_old[0]/epf_old(i,j,k) );
            S_c(i,j,k,1) += 0.5*dt*vm_coeff(i,j,k)*( epfu_divu_old[1]/epf_old(i,j,k) );
            S_c(i,j,k,2) += 0.5*dt*vm_coeff(i,j,k)*( epfu_divu_old[2]/epf_old(i,j,k) );

            S_c(i,j,k,0) += 0.5*dt*vm_coeff(i,j,k)*( epfu_divu_new[0]/epf_new(i,j,k) );
            S_c(i,j,k,1) += 0.5*dt*vm_coeff(i,j,k)*( epfu_divu_new[1]/epf_new(i,j,k) );
            S_c(i,j,k,2) += 0.5*dt*vm_coeff(i,j,k)*( epfu_divu_new[2]/epf_new(i,j,k) );

          }
        });
      } // mfi
    } // step_type
  } //lev
}



/*****************************************************************
* Incorporate velocity source terms                              *
*                                                                *
*    u^(n+1,*) = Sc - Sp*u^(n+1,*)                               *
*                                                                *
******************************************************************/
void FluidUpdate::
ComputeAcceleration ( AdvectionType const a_adv_type,
                      int const a_advect_momentum)
{
  BL_PROFILE("FluidUpdate::ComputeAcceleration");

  int const nspecies( m_fluid.nspecies() );

  AMREX_ASSERT( (nspecies == 0 && !use_species_advection() ) ||
                (nspecies >  0 &&  use_species_advection()  ));

  int const advect_momentum( a_advect_momentum );

  for (int lev(0); lev < nlev(); ++lev) {

    if ( m_step_type == StepType::Predictor ) {

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real      > const& DufDt = stepdata().acceleration(lev,mfi);

        Array4<Real const> const& epf = ( a_adv_type == AdvectionType::MOL ) ?
            leveldata().epf_old_const(lev,mfi) : stepdata().epf_nph_const(lev,mfi);

        Array4<Real const> const& rho = ( a_adv_type == AdvectionType::MOL ) ?
            leveldata().rho_old_const(lev,mfi) : stepdata().rho_nph_const(lev,mfi);

        Array4<Real const> const& vel_old = leveldata().vel_old_const(lev,mfi);
        Array4<Real const> const& vel_new = leveldata().vel_const(lev,mfi);

        Array4<Real const> const& dudt_old = predictor().conv_velocity_const(lev,mfi);

        Array4<Real const> const& drdt_old = !use_species_advection()
            ? predictor().conv_scalars_const(lev,mfi,AdvectionComp::density)
            : predictor().conv_species_const(lev,mfi);

        ParallelFor(bx, [dt=m_dt, advect_momentum, nspecies, DufDt, epf, rho,
        vel_old, vel_new, dudt_old, drdt_old]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

          DufDt(i,j,k,0) = (vel_new(i,j,k,0) - vel_old(i,j,k,0)) / dt;
          DufDt(i,j,k,1) = (vel_new(i,j,k,1) - vel_old(i,j,k,1)) / dt;
          DufDt(i,j,k,2) = (vel_new(i,j,k,2) - vel_old(i,j,k,2)) / dt;

          if (!advect_momentum) {// convective differencing

            DufDt(i,j,k,0) -= (dudt_old(i,j,k,0) / epf(i,j,k));
            DufDt(i,j,k,1) -= (dudt_old(i,j,k,1) / epf(i,j,k));
            DufDt(i,j,k,2) -= (dudt_old(i,j,k,2) / epf(i,j,k));

          } else {

            Real A_rho = drdt_old(i,j,k,0);
            for (int n(1); n < nspecies; ++n)
            {  A_rho += drdt_old(i,j,k,n); }

            Real const rop(epf(i,j,k)*rho(i,j,k));

            DufDt(i,j,k,0) -= (dudt_old(i,j,k,0) - A_rho*vel_old(i,j,k,0)) / rop;
            DufDt(i,j,k,1) -= (dudt_old(i,j,k,1) - A_rho*vel_old(i,j,k,1)) / rop;
            DufDt(i,j,k,2) -= (dudt_old(i,j,k,2) - A_rho*vel_old(i,j,k,2)) / rop;

          }
        });

      }// MFIter

    } else { AMREX_ASSERT( m_step_type == StepType::Corrector );

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real      > const& DufDt = stepdata().acceleration(lev,mfi);

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);

        Array4<Real const> const& vel_old = leveldata().vel_old_const(lev,mfi);
        Array4<Real const> const& vel_new = leveldata().vel_const(lev,mfi);

        Array4<Real const> const& dudt_old = predictor().conv_velocity_const(lev,mfi);
        Array4<Real const> const& dudt_new = corrector().conv_velocity_const(lev,mfi);

        Array4<Real const> const& drdt_old = !use_species_advection()
            ? predictor().conv_scalars_const(lev,mfi,AdvectionComp::density)
            : predictor().conv_species_const(lev,mfi);

        Array4<Real const> const& drdt_new = !use_species_advection()
            ? corrector().conv_scalars_const(lev,mfi,AdvectionComp::density)
            : corrector().conv_species_const(lev,mfi);

        ParallelFor(bx, [dt=m_dt, advect_momentum, nspecies, epf_old, epf_new,
        rho_old, rho_new, vel_old, vel_new, dudt_old, dudt_new, drdt_old, drdt_new,
        DufDt]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

          DufDt(i,j,k,0) = (vel_new(i,j,k,0) - vel_old(i,j,k,0)) / dt;
          DufDt(i,j,k,1) = (vel_new(i,j,k,1) - vel_old(i,j,k,1)) / dt;
          DufDt(i,j,k,2) = (vel_new(i,j,k,2) - vel_old(i,j,k,2)) / dt;

          if (!advect_momentum) {// convective differencing

            DufDt(i,j,k,0) -= 0.5*dt*(dudt_old(i,j,k,0) / epf_old(i,j,k));
            DufDt(i,j,k,1) -= 0.5*dt*(dudt_old(i,j,k,1) / epf_old(i,j,k));
            DufDt(i,j,k,2) -= 0.5*dt*(dudt_old(i,j,k,2) / epf_old(i,j,k));

            DufDt(i,j,k,0) -= 0.5*dt*(dudt_new(i,j,k,0) / epf_new(i,j,k));
            DufDt(i,j,k,1) -= 0.5*dt*(dudt_new(i,j,k,1) / epf_new(i,j,k));
            DufDt(i,j,k,2) -= 0.5*dt*(dudt_new(i,j,k,2) / epf_new(i,j,k));

          } else { // advect momentum

            { Real A_rho_old = drdt_old(i,j,k,0);
              for (int n(1); n < nspecies; ++n)
              { A_rho_old += drdt_old(i,j,k,n); }

              amrex::Real const rop_old(epf_old(i,j,k)*rho_old(i,j,k));

              DufDt(i,j,k,0) -= 0.5*(dudt_old(i,j,k,0) - A_rho_old*vel_old(i,j,k,0)) / rop_old;
              DufDt(i,j,k,1) -= 0.5*(dudt_old(i,j,k,1) - A_rho_old*vel_old(i,j,k,1)) / rop_old;
              DufDt(i,j,k,2) -= 0.5*(dudt_old(i,j,k,2) - A_rho_old*vel_old(i,j,k,2)) / rop_old;
            }

            { Real A_rho_new = drdt_new(i,j,k,0);
              for (int n(1); n < nspecies; ++n)
              { A_rho_new += drdt_new(i,j,k,n); }

              amrex::Real const rop_new(epf_new(i,j,k)*rho_new(i,j,k));

              DufDt(i,j,k,0) -= 0.5*(dudt_new(i,j,k,0) - A_rho_new*vel_new(i,j,k,0)) / rop_new;
              DufDt(i,j,k,1) -= 0.5*(dudt_new(i,j,k,1) - A_rho_new*vel_new(i,j,k,1)) / rop_new;
              DufDt(i,j,k,2) -= 0.5*(dudt_new(i,j,k,2) - A_rho_new*vel_new(i,j,k,2)) / rop_new;

            }
          }
        });


      } // mfi
    } // step_type
  } //lev
}
