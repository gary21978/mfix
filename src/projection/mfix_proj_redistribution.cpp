#include <mfix.H>
#include <mfix_fluid.H>

#include <AMReX_EB_Redistribution.H>

void
mfix::PostProjectionRedistribution (Real l_time, Real l_dt,
                                    const Vector<MultiFab*>& sigma)
{
    if ( m_redistribution_type == "StateRedist") {

      // We must fill internal ghost values before calling redistribution
      // We also need any physical boundary conditions imposed if we are
      //    calling state redistribution (because that calls the slope routine)

      for (int lev(0); lev < nlev(); lev++) {
        bcs().fillpatch(lev, l_time,  BCFillVar::vel, leveldata().vel(), 3);
      }

      for (int lev = 0; lev < nlev(); lev++)
      {
        int ncomp = AMREX_SPACEDIM;

        // Make a temporary to the redistributed velocity into
        MultiFab new_vel(grids[lev], dmap[lev], ncomp, 0);
        new_vel.setVal(0.);

        auto const& bc_vel = bcs().get_hydro_velocity_bcrec_device_ptr();

        for (MFIter mfi(*(leveldata().vel(lev)),TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& fact = m_eb->Factory(lev);

            EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            Array4<Real> const& vel_redist = new_vel.array(mfi);
            Array4<Real> const& vel_orig   = leveldata().vel(lev,mfi);
            Array4<Real> const& gradp = leveldata().grad_p(lev,mfi);
            Array4<Real> const& epf   = leveldata().epf(lev,mfi);
            Array4<Real> const& sig   = sigma[lev]->array(mfi);

            if ( (flagfab.getType(amrex::grow(bx,1)) != FabType::covered) &&
                 (flagfab.getType(amrex::grow(bx,1)) != FabType::regular) )
            {
                Array4<Real const> fcx, fcy, fcz, ccc, vfrac, apx, apy, apz;
                fcx = fact.getFaceCent()[0]->const_array(mfi);
                fcy = fact.getFaceCent()[1]->const_array(mfi);
                fcz = fact.getFaceCent()[2]->const_array(mfi);
                ccc   = fact.getCentroid().const_array(mfi);
                apx = fact.getAreaFrac()[0]->const_array(mfi);
                apy = fact.getAreaFrac()[1]->const_array(mfi);
                apz = fact.getAreaFrac()[2]->const_array(mfi);
                vfrac = fact.getVolFrac().const_array(mfi);

                ApplyInitialRedistribution(bx, ncomp, vel_redist, vel_orig,
                    flag, apx, apy, apz, vfrac,
                    fcx, fcy, fcz, ccc, bc_vel,
                    geom[lev], m_redistribution_type);

                // We update gradp so that (vel_redist + dt gradp_redistnew/rho) == (vel_orig + dt gradp_orig/rho)
                // Note that we do not change rho in the redistribution
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    Real delta_vel = vel_redist(i,j,k,n) - vel_orig(i,j,k,n);
                    gradp(i,j,k,n) -= delta_vel * epf(i,j,k) / (l_dt * sig(i,j,k));
                });

            } else {
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    vel_redist(i,j,k,n) = vel_orig(i,j,k,n);
                });
            }
        }

        // Copy back into the original MultiFab
        MultiFab::Copy(*(leveldata().vel(lev)), new_vel, 0,0, AMREX_SPACEDIM, 0);
      }
    }
}

void
mfix::PreProjectionRedistribution (Real l_time)
{
    if ( m_redistribution_type == "StateRedist") {

      // We must fill internal ghost values before calling redistribution
      // We also need any physical boundary conditions imposed if we are
      //    calling state redistribution (because that calls the slope routine)

      for (int lev(0); lev < nlev(); lev++) {
        bcs().fillpatch(lev, l_time,  BCFillVar::vel, leveldata().vel(), 3);
      }

      for (int lev = 0; lev < nlev(); lev++) {

        int ncomp = AMREX_SPACEDIM;

        // Make a temporary to the redistributed velocity into
        MultiFab new_vel(grids[lev], dmap[lev], ncomp, 0);
        new_vel.setVal(0.);

        auto const& bc_vel = bcs().get_hydro_velocity_bcrec_device_ptr();

        for (MFIter mfi(*leveldata().vel(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& fact = m_eb->Factory(lev);

            EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
            Array4<EBCellFlag const> const& flag = flagfab.const_array();

            Array4<Real> const& vel_redist = new_vel.array(mfi);
            Array4<Real> const& vel_orig   = leveldata().vel(lev,mfi);

            if ( (flagfab.getType(amrex::grow(bx,1)) != FabType::covered) &&
                 (flagfab.getType(amrex::grow(bx,1)) != FabType::regular) )
            {
                Array4<Real const> fcx, fcy, fcz, ccc, vfrac, apx, apy, apz;
                fcx = fact.getFaceCent()[0]->const_array(mfi);
                fcy = fact.getFaceCent()[1]->const_array(mfi);
                fcz = fact.getFaceCent()[2]->const_array(mfi);
                ccc   = fact.getCentroid().const_array(mfi);
                apx = fact.getAreaFrac()[0]->const_array(mfi);
                apy = fact.getAreaFrac()[1]->const_array(mfi);
                apz = fact.getAreaFrac()[2]->const_array(mfi);
                vfrac = fact.getVolFrac().const_array(mfi);

                ApplyInitialRedistribution(bx, ncomp, vel_redist, vel_orig,
                    flag, apx, apy, apz, vfrac,
                    fcx, fcy, fcz, ccc, bc_vel,
                    geom[lev],
                    m_redistribution_type);

            } else {
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    vel_redist(i,j,k,n) = vel_orig(i,j,k,n);
                });
            }
        }

        // Copy back into the original MultiFab
        MultiFab::Copy(*(leveldata().vel(lev)), new_vel, 0,0, AMREX_SPACEDIM, 0);
      }
    }
}
