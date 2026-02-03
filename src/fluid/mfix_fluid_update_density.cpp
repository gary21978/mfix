#include <mfix_fluid_update.H>

using namespace amrex;

void FluidUpdate::Density ( )
{
  int const nspecies( m_fluid.nspecies() );
  int const include_chem( m_include_chem );

  AMREX_ASSERT( (nspecies <= 1 && !use_species_advection() ) ||
                (nspecies >  1 &&  use_species_advection()  ));

  if (m_fluid.solve_density() ) {

    for (int lev(0); lev < nlev(); ++lev) {

      if ( m_step_type == StepType::Predictor ) {

        for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          Box const& bx = mfi.tilebox();

          Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
          Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

          Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
          Array4<Real      > const& rho_new = leveldata().rho(lev,mfi);

          Array4<Real const> const& drdt_old = !use_species_advection()
              ? predictor().conv_scalars_const(lev,mfi,AdvectionComp::density)
              : predictor().conv_species_const(lev,mfi);

          Array4<Real const> const& chem_src = m_include_chem
              ? leveldata().chem_txfr_const(lev, mfi, m_chem_idxs.chem_ro_gk)
              : Array4<Real const>{};

          ParallelFor(bx, [dt=m_dt, nspecies, rho_old, rho_new,
          epf_old, epf_new, drdt_old, include_chem, chem_src ]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Real rho = epf_old(i,j,k)*rho_old(i,j,k);

            Real A_rho = drdt_old(i,j,k,0);
            for (int n(1); n < nspecies; ++n)
            {  A_rho += drdt_old(i,j,k,n); }

            rho += dt * A_rho;

            if (include_chem) {
              for (int n(0); n < nspecies; ++n)
              {  rho += dt * chem_src(i,j,k,n); }
            }

            rho /= epf_new(i,j,k);

            rho_new(i,j,k) = rho;

          });

        } // mfi

      } else { AMREX_ASSERT( m_step_type == StepType::Corrector );

        for (MFIter mfi(*leveldata().rho(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          Box const& bx = mfi.tilebox();

          Array4<Real const> const& epf_old = leveldata().epf_old_const(lev,mfi);
          Array4<Real const> const& epf_new = leveldata().epf_const(lev,mfi);

          Array4<Real const> const& rho_old = leveldata().rho_old_const(lev,mfi);
          Array4<Real      > const& rho_new = leveldata().rho(lev,mfi);

          Array4<Real const> const& drdt_old = !use_species_advection()
              ? predictor().conv_scalars_const(lev,mfi,AdvectionComp::density)
              : predictor().conv_species_const(lev,mfi);

          Array4<Real const> const& drdt_new = !use_species_advection()
              ? corrector().conv_scalars_const(lev,mfi,AdvectionComp::density)
              : corrector().conv_species_const(lev,mfi);

          Array4<Real const> const& chem_src = m_include_chem
              ? leveldata().chem_txfr_const(lev,mfi, m_chem_idxs.chem_ro_gk)
              : Array4<Real const>{};

          ParallelFor(bx, [dt=m_dt, nspecies, rho_old, rho_new,
            epf_old, epf_new, drdt_old, drdt_new, include_chem, chem_src ]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Real rho = epf_old(i,j,k)*rho_old(i,j,k);

            Real A_rho = 0.5*(drdt_old(i,j,k,0) + drdt_new(i,j,k,0));
            for (int n(1); n < nspecies; ++n)
            { A_rho += 0.5*(drdt_old(i,j,k,n) + drdt_new(i,j,k,n)); }

            rho += dt*A_rho;

            if (include_chem) {
              for (int n(0); n < nspecies; ++n)
              {  rho += dt * chem_src(i,j,k,n); }
            }

            rho /= epf_new(i,j,k);

            rho_new(i,j,k) = rho;
          });

        } // mfi
      } // step_type

    } // lev

    // Fill coarse with solution from fine.
    for (int lev( finestLevel()-1 ); lev >= 0; --lev) {
      average_down(lev, leveldata().rho());
    }

    Real const new_time = m_time + m_dt;
    int const ng = (m_step_type == StepType::Corrector) ? 0 : 1;

    for (int lev(0); lev < nlev(); lev++) {

      // Fill ghost cells of new-time density if needed
      // (ghost cells of old density should already  be filled)
      m_bcs.fillpatch(lev, new_time, BCFillVar::rho, leveldata().rho(), ng);

      MultiFab::LinComb( *stepdata().rho_nph(lev), Real(0.5), *leveldata().rho(lev), 0,
          Real(0.5), *leveldata().rho_old(lev), 0, 0, 1, stepdata().rho_nph(lev)->nGrow());
    }

  } else {

    for (int lev(0); lev < nlev(); ++lev) {

      MultiFab::Copy(*(stepdata().rho_nph(lev)), *(leveldata().rho_old(lev)),
          0, 0, 1, stepdata().rho_nph(lev)->nGrow());
    }

  }

}
