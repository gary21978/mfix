#include <mfix_fluid_update.H>

using namespace amrex;

void FluidUpdate::
Tracers ( int const a_nscalars,
          MFIXDiffOpTracer* a_diffOp)
{
  int const nscalars( a_nscalars );
  int const adv_comp ( AdvectionComp::tracer );

  int const explicit_diffusion( m_explicit_diffusion );

  for (int lev(0); lev < nlev(); ++lev) {

    if ( m_step_type == StepType::Predictor ) {

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);

        Array4<Real const> const& S_old = leveldata().tracer_old_const(lev,mfi);
        Array4<Real      > const& S_new = leveldata().tracer(lev,mfi);

        Array4<Real const> const& dSdt_old = predictor().conv_scalars_const(lev,mfi,adv_comp);

        Array4<Real const> const& lapS_old = predictor().lap_tracer_const(lev,mfi);

        ParallelFor(bx, nscalars, [dt=m_dt, rho_old, rho_new, epf_old, epf_new,
        S_old, S_new, dSdt_old, lapS_old, explicit_diffusion ]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          Real S = epf_old(i,j,k)*rho_old(i,j,k)*S_old(i,j,k,n);

          S += dt * dSdt_old(i,j,k,n);

          if (explicit_diffusion) { S += dt * lapS_old(i,j,k,n); }

          S_new(i,j,k,n) = S / (epf_new(i,j,k)*rho_new(i,j,k));
        });

      } // mfi

    } else { AMREX_ASSERT( m_step_type == StepType::Corrector );

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);

        Array4<Real const> const& S_old = leveldata().tracer_old_const(lev,mfi);
        Array4<Real      > const& S_new = leveldata().tracer(lev,mfi);

        Array4<Real const> const& dSdt_old = predictor().conv_scalars_const(lev,mfi,adv_comp);
        Array4<Real const> const& dSdt_new = corrector().conv_scalars_const(lev,mfi,adv_comp);

        Array4<Real const> const& lapS_old = predictor().lap_tracer_const(lev,mfi);
        Array4<Real const> const& lapS_new = corrector().lap_tracer_const(lev,mfi);

        ParallelFor(bx, nscalars, [dt=m_dt, rho_old, rho_new, epf_old, epf_new,
        S_old, S_new, dSdt_old, dSdt_new, lapS_old, lapS_new, explicit_diffusion ]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          Real S = epf_old(i,j,k)*rho_old(i,j,k)*S_old(i,j,k,n);
          S += dt * 0.5*(dSdt_old(i,j,k,n) + dSdt_new(i,j,k,n));

          // Crank-Nicolson: add half of the explicit old lap
          S += dt* (0.5*lapS_old(i,j,k,n));

          // Full explicit add half of new lap too
          if (explicit_diffusion) { S += dt*(0.5*lapS_new(i,j,k,n)); }

          S_new(i,j,k,n) = S / (epf_new(i,j,k)*rho_new(i,j,k));
        });

      } // mfi
    } // step_type

    leveldata().tracer(lev)->FillBoundary(m_geom[lev].periodicity());

  } // lev

  for (int lev(0); lev < nlev(); lev++) {
    m_bcs.set_tracer_bcs(lev, m_time + m_dt, m_fluid, leveldata().tracer(lev));
  }

  if ( !explicit_diffusion ) {

    Real const dt( m_dt*(m_step_type == StepType::Predictor ? 1.0 : 0.5) );

    a_diffOp->solve( leveldata().tracer(), leveldata().epf_const(),
        leveldata().rho_const(), dt);

    for (int lev(0); lev < nlev(); lev++) {
      m_bcs.set_tracer_bcs(lev, m_time + m_dt, m_fluid, leveldata().tracer(lev));
    }

  } else {

   // Need to average down tracers when diffusion is explicit because
   // the diffusion solver didn't do it for us.
    for (int lev( finestLevel()-1 ); lev >= 0; --lev) {
      average_down(lev, leveldata().tracer());
    }
  }
}
