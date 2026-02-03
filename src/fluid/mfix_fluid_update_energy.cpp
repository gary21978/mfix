#include <mfix_fluid_update.H>
#include <mfix_solvers.H>

using namespace amrex;

void FluidUpdate::
Energy (
         Real const& a_pT_old,
         Real      & a_pT_new,
         MFIXDiffOpEnergy* a_diffOp,
         Real const a_abstol, Real const a_reltol, int const a_maxiter )
{
  int const is_mixture( m_fluid.isMixture() );
  int const nspecies( m_fluid.nspecies() );
  int const include_chem( m_include_chem );

  int const explicit_diffusion( m_explicit_diffusion );
  int const species_diffusion( m_fluid.solve_species() );

  int const adv_comp(AdvectionComp::enthalpy);
  int const chem_comp( m_include_chem ? m_chem_idxs.chem_h : 0 );

  const auto fluid_props = m_fluid.props.data<run_on>();

  Real const atol(a_abstol);
  Real const rtol(a_reltol);

  int  const maxiter(a_maxiter);

  // Update thermodynamic pressure
  Real dpdt(0.);
  if (m_fluid.constraint.isIdealGasClosedSystem()) {
    dpdt = ((m_step_type == StepType::Predictor) ? predictor().RHS_therm_p
        : 0.5*(predictor().RHS_therm_p + corrector().RHS_therm_p));
  }
  a_pT_new = a_pT_old + m_dt*dpdt;

  // Compute enthalpy flux from species diffusion
  if (!explicit_diffusion && m_fluid.solve_species()) {

    if ( m_step_type == StepType::Predictor ) {

      a_diffOp->computeDivhJ(predictor().div_hJ(), stepdata().h_k(),
          stepdata().J_k(), leveldata().T_old_const() );

    } else { AMREX_ASSERT( m_step_type == StepType::Corrector );

      a_diffOp->computeDivhJ(corrector().div_hJ(), stepdata().h_k(),
          stepdata().J_k(), leveldata().T_const() );
    }
  }

  for (int lev(0); lev < nlev(); ++lev) {

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(leveldata().rho(lev)->Factory());

    if ( m_step_type == StepType::Predictor ) {

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);

        Array4<Real const> const& h_old = leveldata().h_old_const(lev,mfi);
        Array4<Real      > const& h_new = leveldata().h(lev,mfi);

        Array4<Real const> const& T_old = leveldata().T_old_const(lev,mfi);
        Array4<Real      > const& T_new = leveldata().T(lev,mfi);

        Array4<Real const> const& dhdt_old = predictor().conv_scalars_const(lev,mfi,adv_comp);

        Array4<Real const> const& lapT_old = predictor().lap_T_const(lev,mfi);

        Array4<Real const> const& divhJ_old = predictor().div_hJ_const(lev,mfi);

        Array4<Real const> const& chem_src = leveldata().chem_txfr_const(lev,mfi,chem_comp);

        Array4<Real const> const& X_new = leveldata().X_const(lev,mfi);

        EBCellFlagFab const& flagfab = factory.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flags = flagfab.const_array();

        ParallelFor(bx, [dt=m_dt, explicit_diffusion, species_diffusion, nspecies,
        epf_old, epf_new, rho_old, rho_new, X_new, h_old, h_new, T_old, T_new,
        dhdt_old, lapT_old, divhJ_old, dpdt, chem_src, flags, fluid_props,
        include_chem, is_mixture, atol, rtol, maxiter]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          int const is_covered = static_cast<int>(flags(i,j,k).isCovered());

          if (!is_covered) {

            Real h = epf_old(i,j,k)*rho_old(i,j,k)* h_old(i,j,k);

            h += dt * dhdt_old(i,j,k);

            if (explicit_diffusion) { h += dt * lapT_old(i,j,k); }
            if (species_diffusion)  { h -= dt * divhJ_old(i,j,k); }
            if (include_chem)       { h += dt * chem_src(i,j,k); }

            h += dt * epf_new(i,j,k) * dpdt;

            h_new(i,j,k) = h / (epf_new(i,j,k)*rho_new(i,j,k));

            // Newton-Raphson solver for solving implicit equation for
            // temperature
            Solvers::Newton::FluidEnthalpy::Residue residue(i, j, k,
                is_covered, fluid_props, X_new, h_new(i,j,k));

            Solvers::Newton::FluidEnthalpy::Gradient gradient(i, j, k, fluid_props, X_new);

            amrex::Real T(T_old(i,j,k));

            auto output = Solvers::Newton::solve(T, residue, gradient,
                atol, rtol, maxiter);

            if (output.iterations == -1) {
              amrex::Abort("FluidUpdate::Energy StepType::Predictor\n!!!Newton solver did not converge!!!");
            }

            T_new(i,j,k) = T;
          }

        });

      } // mfi

    } else { AMREX_ASSERT( m_step_type == StepType::Corrector );

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);

        Array4<Real const> const& h_old = leveldata().h_old_const(lev,mfi);
        Array4<Real      > const& h_new = leveldata().h(lev,mfi);

        Array4<Real const> const& T_old = leveldata().T_old_const(lev,mfi);
        Array4<Real      > const& T_new = leveldata().T(lev,mfi);

        Array4<Real const> const& X_new = leveldata().X_const(lev,mfi);

        Array4<Real const> const& dhdt_old = predictor().conv_scalars_const(lev,mfi,adv_comp);
        Array4<Real const> const& dhdt_new = corrector().conv_scalars_const(lev,mfi,adv_comp);

        Array4<Real const> const& lapT_old = predictor().lap_T_const(lev,mfi);
        Array4<Real const> const& lapT_new = corrector().lap_T_const(lev,mfi);

        Array4<Real const> const& chem_src = leveldata().chem_txfr_const(lev,mfi,chem_comp);

        Array4<Real const> const& divhJ_old = predictor().div_hJ_const(lev,mfi);
        Array4<Real const> const& divhJ_new = corrector().div_hJ_const(lev,mfi);

        EBCellFlagFab const& flagfab = factory.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flags = flagfab.const_array();

        ParallelFor(bx, [dt=m_dt, explicit_diffusion, species_diffusion, nspecies,
        epf_old, epf_new, rho_old, rho_new, X_new, h_old, h_new, T_old, T_new,
        dhdt_old, dhdt_new, lapT_old, lapT_new, divhJ_old, divhJ_new, dpdt, chem_src,
        flags, fluid_props, include_chem, is_mixture, atol, rtol, maxiter]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          int const is_covered = static_cast<int>(flags(i,j,k).isCovered());

          if (!is_covered) {

            Real h = epf_old(i,j,k)*rho_old(i,j,k)*h_old(i,j,k);

            h += dt*0.5*(dhdt_old(i,j,k) + dhdt_new(i,j,k));

            if (include_chem) { h += dt*chem_src(i,j,k); }

            if (species_diffusion) { h -= dt*0.5*(divhJ_old(i,j,k) + divhJ_new(i,j,k)); }

            // Crank-Nicolson so we only add half of the diffusive term
            h += dt*0.5*lapT_old(i,j,k);
            if (explicit_diffusion) { h += dt*0.5*lapT_new(i,j,k); }

            h += dt * epf_new(i,j,k)*dpdt;

            h_new(i,j,k) = h / (epf_new(i,j,k)*rho_new(i,j,k));

            // Newton-Raphson solver for solving implicit equation for
            // temperature
            Solvers::Newton::FluidEnthalpy::Residue residue(i, j, k,
                is_covered, fluid_props, X_new, h_new(i,j,k));

            Solvers::Newton::FluidEnthalpy::Gradient gradient(i, j, k, fluid_props, X_new);

            amrex::Real T(T_old(i,j,k));

            auto output = Solvers::Newton::solve(T, residue, gradient, atol, rtol, maxiter);

            if (output.iterations == -1) {
              amrex::Abort("FluidUpdate::Energy StepType::Corrector\n!!!Newton solver did not converge!!!");
            }
            T_new(i,j,k) = T;

          }
        });

      } // mfi
    } // step_type
  } // lev

  // We only diffuse temperature so we need to average_down enthalpy ourselves.
  for (int lev( finestLevel()-1 ); lev >= 0; --lev) {
    average_down(lev, leveldata().T());
  }

  if ( !explicit_diffusion ) {

    Real const new_time = m_time+m_dt;

    for (int lev(0); lev < nlev(); lev++) {
      m_bcs.set_temperature_bcs(lev, new_time, m_fluid, leveldata().T(lev));
      m_bcs.set_enthalpy_bcs(lev, new_time, m_fluid, leveldata().h(lev));
    }

    Real const dt( m_dt*(m_step_type == StepType::Predictor ? 1.0 : 0.5) );

    // Diffuse temperature
    a_diffOp->solve(leveldata().T(), leveldata().h(), leveldata().epf_const(),
        leveldata().rho_const(), leveldata().X_const(),
        stepdata().eb_temperature_const(), dt,
        a_abstol, a_reltol, a_maxiter);

    // TODO: Do we need to rest the BC values here?
    for (int lev(0); lev < nlev(); lev++) {
      m_bcs.set_temperature_bcs(lev, new_time, m_fluid, leveldata().T(lev));
      m_bcs.set_enthalpy_bcs(lev, new_time, m_fluid, leveldata().h(lev));
    }

  } else {

    // Average down the temperature when diffusion is explicit because
    // the diffusion solver didn't do it for us.
    for (int lev( finestLevel()-1 ); lev >= 0; --lev) {
      average_down(lev, leveldata().T());
    }
  }

}
