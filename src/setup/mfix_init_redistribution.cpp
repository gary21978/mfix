#include <mfix.H>
#include <mfix_fluid.H>

#include <AMReX_EB_Redistribution.H>

void mfix::
InitialRedistribution (Real a_time)
{
  // Next we must redistribute the initial solution if we are going to use
  // StateRedist redistribution scheme
  if ( m_redistribution_type == "StateRedist") {

    // We use the "old" data as the input here
    // We must fill internal ghost values before calling redistribution
    // We also need any physical boundary conditions imposed if we are
    //    calling state redistribution (because that calls the slope routine)

    for (int lev(0); lev < nlev(); lev++) {

      bcs().fillpatch(lev, a_time, BCFillVar::vel,
          leveldata().vel(), leveldata().vel_old(lev), 3);

      bcs().fillpatch(lev, a_time, BCFillVar::rho,
          leveldata().rho(), leveldata().rho_old(lev), 3);

      if (leveldata(lev)->has_enthalpy()) {
        bcs().fillpatch(lev, a_time, BCFillVar::h,
            leveldata().h(), leveldata().h_old(lev), 3);
      }
      if (leveldata(lev)->has_temperature()) {
        bcs().fillpatch(lev, a_time, BCFillVar::T,
            leveldata().T(), leveldata().T_old(lev), 3);
      }
      if (leveldata(lev)->has_species()) {
        bcs().fillpatch(lev, a_time, BCFillVar::X,
            leveldata().X(), leveldata().X_old(lev), 3);
      }
      if (leveldata(lev)->has_tracer()) {
        bcs().fillpatch(lev, a_time, BCFillVar::tracer,
            leveldata().tracer(), leveldata().tracer_old(lev), 3);
      }
    }

    for (int lev(0); lev < nlev(); lev++) {

      for (MFIter mfi(*(leveldata().rho(lev)),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();
        auto const& fact = m_eb->Factory(lev);

        EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flag = flagfab.const_array();

        if ( (flagfab.getType(bx)                != FabType::covered) &&
             (flagfab.getType(amrex::grow(bx,4)) != FabType::regular) ) {

          Array4<Real const> fcx, fcy, fcz, ccc, vfrac, apx, apy, apz;
          fcx = fact.getFaceCent()[0]->const_array(mfi);
          fcy = fact.getFaceCent()[1]->const_array(mfi);
          fcz = fact.getFaceCent()[2]->const_array(mfi);
          ccc   = fact.getCentroid().const_array(mfi);
          apx = fact.getAreaFrac()[0]->const_array(mfi);
          apy = fact.getAreaFrac()[1]->const_array(mfi);
          apz = fact.getAreaFrac()[2]->const_array(mfi);
          vfrac = fact.getVolFrac().const_array(mfi);

          int ncomp = AMREX_SPACEDIM;

          auto const& bc_vel = bcs().get_hydro_velocity_bcrec_device_ptr();
          ApplyInitialRedistribution(bx,ncomp,
              leveldata().vel(lev,mfi), leveldata().vel_old(lev,mfi),
              flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
              bc_vel, geom[lev], m_redistribution_type);

          if (fluid.solve_density()) {

            ncomp = 1;

            auto const& bc_den = bcs().get_density_bcrec_device_ptr();
            ApplyInitialRedistribution( bx, ncomp,
                leveldata().rho(lev,mfi), leveldata().rho_old(lev,mfi),
                flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                bc_den, geom[lev], m_redistribution_type);
          }
          if (fluid.solve_enthalpy()) {

            ncomp = 1;

            auto const& bc_h = bcs().get_enthalpy_bcrec_device_ptr();
            ApplyInitialRedistribution(bx,ncomp,
                leveldata().h(lev,mfi), leveldata().h_old(lev,mfi),
                flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                bc_h, geom[lev],
                m_redistribution_type);
          }
          if (fluid.solve_tracer()) {

            ncomp = fluid.ntracer();

            auto const& bc_t = bcs().get_tracer_bcrec_device_ptr();
            ApplyInitialRedistribution(bx,ncomp,
                leveldata().tracer(lev,mfi), leveldata().tracer_old(lev,mfi),
                flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                bc_t, geom[lev], m_redistribution_type);
          }
          if (fluid.solve_species()) {

            ncomp = fluid.nspecies();

            auto const& bc_X = bcs().get_species_bcrec_device_ptr();
            ApplyInitialRedistribution(bx,ncomp,
                leveldata().X(lev,mfi), leveldata().X_old(lev,mfi),
                flag, apx, apy, apz, vfrac, fcx, fcy, fcz, ccc,
                bc_X, geom[lev], m_redistribution_type);
          }
        }
      }

      // We fill internal ghost values after calling redistribution
      leveldata().vel(lev)->FillBoundary();
      leveldata().rho(lev)->FillBoundary();

      // Since we redistributed enthalpy should we recompute temperature?
      if (fluid.solve_enthalpy()) { leveldata().h(lev)->FillBoundary(); }
      if (fluid.solve_tracer()) { leveldata().tracer(lev)->FillBoundary(); }
      if (fluid.solve_species()) { leveldata().X(lev)->FillBoundary(); }
    }
  }
}
