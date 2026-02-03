#include <mfix_fluid_update.H>

using namespace amrex;

void FluidUpdate::
Species ( MFIXDiffOpSpecies* a_diffOp )
{
  int const nspecies( m_fluid.nspecies() );
  AMREX_ASSERT( nspecies > 0 );

  int const explicit_diffusion( m_explicit_diffusion );

  int const include_chem( m_include_chem );
  int const chem_comp( m_include_chem ? m_chem_idxs.chem_ro_gk : -1 );

  for (int lev(0); lev < nlev(); ++lev) {

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(
        leveldata().X(lev)->Factory() );

    if ( m_step_type == StepType::Predictor ) {

      for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
        Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

        Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
        Array4<Real const> const& rho_new = leveldata().rho_const(lev,mfi);

        Array4<Real const> const& X_old = leveldata().X_old_const(lev,mfi);
        Array4<Real      > const& X_new = leveldata().X(lev,mfi);

        Array4<Real const> const& dXdt_old = predictor().conv_species_const(lev,mfi);

        Array4<Real const> const& divJ_old = predictor().div_J_const(lev,mfi);

        Array4<Real const> const& chem_src = leveldata().chem_txfr_const(lev, mfi, chem_comp);

        EBCellFlagFab const& flagfab = factory.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flags = flagfab.const_array();

        ParallelFor(bx, [dt=m_dt, explicit_diffusion, nspecies, include_chem,
        rho_old, rho_new, epf_old, epf_new, X_old, X_new, dXdt_old,
        divJ_old, chem_src, flags]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real sum_X(0.);

          for (int n(0); n < nspecies; ++n) {

            Real Xn = epf_old(i,j,k)*rho_old(i,j,k)*X_old(i,j,k,n);

            Xn += dt * dXdt_old(i,j,k,n);

            if (include_chem) { Xn += dt * chem_src(i,j,k,n); }

            if (explicit_diffusion) {
              Xn -= dt * divJ_old(i,j,k,n);
            }

            Xn /= (epf_new(i,j,k)*rho_new(i,j,k));

            Xn = amrex::Clamp(Xn, 0., 1.);
            sum_X += Xn;

            X_new(i,j,k,n) = Xn;
          }

          if (!flags(i,j,k).isCovered()) {
            for (int n = 0; n < nspecies; ++n) {
              X_new(i,j,k,n) /= sum_X;
            }
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

        Array4<Real const> const& X_old = leveldata().X_old_const(lev,mfi);
        Array4<Real      > const& X_new = leveldata().X(lev,mfi);

        Array4<Real const> const& dXdt_old = predictor().conv_species_const(lev,mfi);
        Array4<Real const> const& dXdt_new = corrector().conv_species_const(lev,mfi);

        Array4<Real const> const& divJ_old = predictor().div_J_const(lev,mfi);
        Array4<Real const> const& divJ_new = corrector().div_J_const(lev,mfi);

        Array4<Real const> const& chem_src = leveldata().chem_txfr_const(lev,mfi,chem_comp);

        EBCellFlagFab const& flagfab = factory.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flags = flagfab.const_array();

        ParallelFor(bx, [dt=m_dt, nspecies, explicit_diffusion,
        rho_old, rho_new, epf_old, epf_new, X_old, X_new, dXdt_old, dXdt_new,
        divJ_old, divJ_new, include_chem, chem_src, flags]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real sum_X(0.);

          for (int n(0); n < nspecies; ++n) {

            Real Xn = epf_old(i,j,k)*rho_old(i,j,k)*X_old(i,j,k,n);

            Xn += 0.5*dt*(dXdt_old(i,j,k,n) + dXdt_new(i,j,k,n));

            if (include_chem) { Xn += dt * chem_src(i,j,k,n); }

            if (explicit_diffusion) {
              Xn -= 0.5*dt*(divJ_old(i,j,k,n)+divJ_new(i,j,k,n));
            } else {
              // Crank-Nicolson so we only add half of the diffusive term
              Xn -= 0.5*dt * divJ_old(i,j,k,n);
            }

            Xn /= (epf_new(i,j,k)*rho_new(i,j,k));

            Xn = amrex::Clamp(Xn, 0.,1.);
            sum_X += Xn;

            X_new(i,j,k,n) = Xn;
          }

          if (!flags(i,j,k).isCovered()) {
            for (int n(0); n < nspecies; ++n) {
              X_new(i,j,k,n) /= sum_X;
            }
          }
        });

      } // mfi
    } // step_type

    //leveldata().X(lev)->FillBoundary(m_geom[lev].periodicity());

  } // lev


  if ( !explicit_diffusion ) {

    // When using implicit diffusion for species, we "Add" (subtract) the
    // correction term computed at time t^{star,star} to the RHS before
    // doing the implicit diffusion
    for (int lev(0); lev < nlev(); lev++) {
      m_bcs.set_species_bcs(lev, m_time + m_dt, nspecies, leveldata().X(lev));
    }

    Real const dt( m_dt*(m_step_type == StepType::Predictor ? 1.0 : 0.5) );

    // Diffuse species mass fractions
    a_diffOp->solve(leveldata().X(), stepdata().J_k(),
        leveldata().epf_const(), leveldata().rho_const(), dt);

    for (int lev(0); lev < nlev(); lev++) {
      m_bcs.set_species_bcs(lev, m_time + m_dt, nspecies, leveldata().X(lev));
    }

  } else {

   // Need to average down species when diffusion is explicit because
   // the diffusion solver didn't do it for us.
    for (int lev( finestLevel()-1 ); lev >= 0; --lev) {
      average_down(lev, leveldata().X());
    }

  }

  // Update ghost cells
  for (int lev(0); lev < nlev(); lev++) {
    leveldata().X(lev)->FillBoundary(m_geom[lev].periodicity());
  }
}
