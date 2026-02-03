#include <mfix.H>
#include <mfix_bc_fillpatch_K.H>
#include <mfix_compute_cell_grads.H>

#include <AMReX_FillPatchUtil.H>

void
MFIXReadWrite::FillPatchVelForVort (int lev, MultiFab& vel)
{
   vel.setVal(mfix::covered_val);

   const int minf = BCList::minf;
   const int cover = BCList::cover;
   const int ng = 2;
   const int max_level = nlev - 1;
   const Real dummy_time = 0.;

   if (lev == 0) {
      PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > physbc
        (geom[lev], m_bcrec_velocity, MFIXVelFill{minf, cover, max_level,
            m_boundary_conditions.bc_u_g().data(),
            m_boundary_conditions.bc_v_g().data(),
            m_boundary_conditions.bc_w_g().data(),
            bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
            bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
            bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
            });

      FillPatchSingleLevel(vel, IntVect(ng), dummy_time,
                           {leveldata().vel(lev)},
                           {dummy_time}, 0, 0, vel.nComp(), geom[lev],
                           physbc, 0);
   }
   else
   {
     PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > cphysbc
       (geom[lev-1],m_bcrec_velocity, MFIXVelFill{minf, cover, max_level,
           m_boundary_conditions.bc_u_g().data(),
           m_boundary_conditions.bc_v_g().data(),
           m_boundary_conditions.bc_w_g().data(),
           bc_list.bc_ilo[lev-1]->array(), bc_list.bc_ihi[lev-1]->array(),
           bc_list.bc_jlo[lev-1]->array(), bc_list.bc_jhi[lev-1]->array(),
           bc_list.bc_klo[lev-1]->array(), bc_list.bc_khi[lev-1]->array()
           });

     PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > fphysbc
       (geom[lev], m_bcrec_velocity, MFIXVelFill{minf, cover, max_level,
           m_boundary_conditions.bc_u_g().data(),
           m_boundary_conditions.bc_v_g().data(),
           m_boundary_conditions.bc_w_g().data(),
           bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
           bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
           bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
           });

     Interpolater* mapper = (Interpolater*)(&eb_cell_cons_interp);

     FillPatchTwoLevels(vel, IntVect(ng), dummy_time,
                        {leveldata().vel(lev-1)},
                        {dummy_time},
                        {leveldata().vel(lev)},
                        {dummy_time},
                        0, 0, vel.nComp(), geom[lev-1], geom[lev],
                        cphysbc, 0, fphysbc, 0,
                        ref_ratio[lev-1], mapper, m_bcrec_velocity, 0);

   }
}

void
MFIXReadWrite::ComputeVort ()
{
    BL_PROFILE("mfix::mfix_compute_vort");

    Vector<std::unique_ptr<MultiFab> > vel(nlev);

    if (nlev > 1) {
       for (int lev = 0; lev < nlev; ++lev) {
          // Make a copy to do fillpatch
          vel[lev] = std::make_unique<MultiFab>(leveldata().vel(lev)->boxArray(),
                leveldata().vel(lev)->DistributionMap(),
                leveldata().vel(lev)->nComp(), leveldata().vel(lev)->nGrow(),
                MFInfo(), leveldata().vel(lev)->Factory());
       }

       for (int lev = 0; lev < nlev; ++lev) {
          FillPatchVelForVort(lev, *vel[lev]);
       }
    }

    for (int lev = 0; lev < nlev; lev++)
    {
      leveldata().vorticity(lev)->setVal(0.0);

       for (MFIter mfi(*(leveldata().vel(lev)),TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          // Tilebox
          Box bx = mfi.tilebox ();
          const RealVect dx(geom[lev].CellSize()[0],
                geom[lev].CellSize()[1], geom[lev].CellSize()[2]);

          Array4<Real> const& vorticity = leveldata().vorticity(lev,mfi);
          Array4<Real> const& velocity_g = (nlev > 1) ? vel[lev]->array(mfi) : leveldata().vel(lev,mfi);

          const Real odx(1./dx[0]), ody(1./dx[1]), odz(1./dx[2]);

          // This is to check efficiently if this tile contains any eb stuff
          const EBFArrayBox& vel_fab =
            static_cast<EBFArrayBox const&>((*(leveldata().vel(lev)))[mfi]);

          const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

          if (flags.getType(amrex::grow(bx,0)) == FabType::regular)
          {
            amrex::ParallelFor(bx, [odx,ody,odz,velocity_g,vorticity]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              Real uy = .5*ody*(velocity_g(i,j+1,k,0) - velocity_g(i,j-1,k,0));
              Real uz = .5*odz*(velocity_g(i,j,k+1,0) - velocity_g(i,j,k-1,0));
              Real vx = .5*odx*(velocity_g(i+1,j,k,1) - velocity_g(i-1,j,k,1));
              Real vz = .5*odz*(velocity_g(i,j,k+1,1) - velocity_g(i,j,k-1,1));
              Real wx = .5*odx*(velocity_g(i+1,j,k,2) - velocity_g(i-1,j,k,2));
              Real wy = .5*ody*(velocity_g(i,j+1,k,2) - velocity_g(i,j-1,k,2));

              vorticity(i,j,k) = std::sqrt((wy-vz)*(wy-vz) +
                                           (uz-wx)*(uz-wx) +
                                           (vx-uy)*(vx-uy));
            });
          }
          // COPIED FROM INCFLO
          else if (flags.getType(amrex::grow(bx,0)) == FabType::singlevalued)
          {
            const auto& flag_fab = flags.const_array();
            amrex::ParallelFor(bx, [dx,velocity_g,vorticity,flag_fab]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (flag_fab(i,j,k).isCovered())
                {
                    vorticity(i,j,k) = Real(0.0);
                }
                else
                {
                    Real ux, vx, wx, uy, vy, wy, uz, vz, wz;
                    mfix_comp_cell_grads(i, j, k, ux, uy, uz, vx, vy, vz,
                          wx, wy, wz, velocity_g, flag_fab, dx);

                    vorticity(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
                }
            });
          }
          else if (flags.getType(amrex::grow(bx,0)) == FabType::covered)
          {
            amrex::ParallelFor(bx, [vorticity]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               vorticity(i,j,k) = Real(0.0);
            });
          }
       }
    }
}
