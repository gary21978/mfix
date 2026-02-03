#include <AMReX_FillPatchUtil.H>

#include <mfix.H>
#include <mfix_bc_fillpatch_K.H>

//
// Set the BCs for face-centroid-based velocity components only
//
void
mfix::set_MAC_velocity_bcs ( int a_lev,
                             MultiFab const* a_mac_rhs,
                             Array<MultiFab*,3> a_mac_vel)
{
  BL_PROFILE("MacProjection::set_MAC_velocity_bcs()");

  Box domain(geom[a_lev].Domain());

  for (MFIter mfi(*a_mac_rhs, false); mfi.isValid(); ++mfi)
  {
    const Box& ubx = (*a_mac_vel[0])[mfi].box();
    IntVect ubx_lo(ubx.loVect());
    IntVect ubx_hi(ubx.hiVect());

    const Box& vbx = (*a_mac_vel[1])[mfi].box();
    IntVect vbx_lo(vbx.loVect());
    IntVect vbx_hi(vbx.hiVect());

    const Box& wbx = (*a_mac_vel[2])[mfi].box();
    IntVect wbx_lo(wbx.loVect());
    IntVect wbx_hi(wbx.hiVect());

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    Array4<Real> const& ep_u = a_mac_vel[0]->array(mfi);
    Array4<Real> const& ep_v = a_mac_vel[1]->array(mfi);
    Array4<Real> const& ep_w = a_mac_vel[2]->array(mfi);

    Array4<int> const& bct_ilo = bc_list.bc_ilo[a_lev]->array();
    Array4<int> const& bct_ihi = bc_list.bc_ihi[a_lev]->array();
    Array4<int> const& bct_jlo = bc_list.bc_jlo[a_lev]->array();
    Array4<int> const& bct_jhi = bc_list.bc_jhi[a_lev]->array();
    Array4<int> const& bct_klo = bc_list.bc_klo[a_lev]->array();
    Array4<int> const& bct_khi = bc_list.bc_khi[a_lev]->array();

    const int nlft = amrex::max(0, dom_lo[0]-ubx_lo[0]);
    const int nbot = amrex::max(0, dom_lo[1]-vbx_lo[1]);
    const int ndwn = amrex::max(0, dom_lo[2]-wbx_lo[2]);

    const int nrgt = amrex::max(0, ubx_hi[0]-dom_hi[0]);
    const int ntop = amrex::max(0, vbx_hi[1]-dom_hi[1]);
    const int nup  = amrex::max(0, wbx_hi[2]-dom_hi[2]);

    amrex::Real* p_bc_u_g = m_boundary_conditions.bc_u_g().data();
    amrex::Real* p_bc_v_g = m_boundary_conditions.bc_v_g().data();
    amrex::Real* p_bc_w_g = m_boundary_conditions.bc_w_g().data();
    amrex::Real* p_bc_epf = m_boundary_conditions.bc_epf().data();

    // NOTE - we only call this for MAC velocities which are only defined on normal faces

    if (nlft > 0)
    {
      // Create InVects for following Box
      IntVect ulo_bx_yz_lo(ubx_lo);
      IntVect ulo_bx_yz_hi(ubx_hi);

      // Fix lo and hi limits
      ulo_bx_yz_lo[0] = dom_lo[0];
      ulo_bx_yz_hi[0] = dom_lo[0];

      const Box ulo_bx_yz(ulo_bx_yz_lo, ulo_bx_yz_hi);

      amrex::ParallelFor(ulo_bx_yz,
        [bct_ilo,dom_lo,p_bc_u_g,p_bc_epf,ep_u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_ilo(dom_lo[0]-1,j,k,1);
        const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
        if(bct == BCList::minf) ep_u(i,j,k) = p_bc_u_g[bcv] * p_bc_epf[bcv];
      });
    }

    if (nrgt > 0)
    {
      // Create InVects for following Box
      IntVect uhi_bx_yz_lo(ubx_lo);
      IntVect uhi_bx_yz_hi(ubx_hi);

      // Fix lo and hi limits
      uhi_bx_yz_lo[0] = dom_hi[0]+1;
      uhi_bx_yz_hi[0] = dom_hi[0]+1;

      const Box uhi_bx_yz(uhi_bx_yz_lo, uhi_bx_yz_hi);

      amrex::ParallelFor(uhi_bx_yz,
        [bct_ihi,dom_hi,p_bc_u_g,p_bc_epf,ep_u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_ihi(dom_hi[0]+1,j,k,1);
        const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
        if (bct == BCList::minf) ep_u(i,j,k) = p_bc_u_g[bcv] * p_bc_epf[bcv];
      });
    }

    if (nbot > 0)
    {
      // Create InVects for following Box
      IntVect vlo_bx_xz_lo(vbx_lo);
      IntVect vlo_bx_xz_hi(vbx_hi);

      // Fix lo and hi limits
      vlo_bx_xz_lo[1] = dom_lo[1];
      vlo_bx_xz_hi[1] = dom_lo[1];

      const Box vlo_bx_xz(vlo_bx_xz_lo, vlo_bx_xz_hi);

      amrex::ParallelFor(vlo_bx_xz,
        [bct_jlo,dom_lo,p_bc_v_g,p_bc_epf,ep_v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_jlo(i,dom_lo[1]-1,k,1);
        const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
        if (bct == BCList::minf) ep_v(i,j,k) = p_bc_v_g[bcv] * p_bc_epf[bcv];
      });
    }

    if (ntop > 0)
    {
      // Create InVects for following Box
      IntVect vhi_bx_xz_lo(vbx_lo);
      IntVect vhi_bx_xz_hi(vbx_hi);

      // Fix lo and hi limits
      vhi_bx_xz_lo[1] = dom_hi[1]+1;
      vhi_bx_xz_hi[1] = dom_hi[1]+1;

      const Box vhi_bx_xz(vhi_bx_xz_lo, vhi_bx_xz_hi);

      amrex::ParallelFor(vhi_bx_xz,
        [bct_jhi,dom_hi,p_bc_v_g,p_bc_epf,ep_v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_jhi(i,dom_hi[1]+1,k,1);
        const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
        if (bct == BCList::minf) ep_v(i,j,k) = p_bc_v_g[bcv] * p_bc_epf[bcv];
      });
    }

    if (ndwn > 0)
    {
      // Create InVects for following Boxes
      IntVect wlo_bx_xy_lo(wbx_lo);
      IntVect wlo_bx_xy_hi(wbx_hi);

      // Fix lo and hi limits
      wlo_bx_xy_lo[2] = dom_lo[2];
      wlo_bx_xy_hi[2] = dom_lo[2];

      const Box wlo_bx_xy(wlo_bx_xy_lo, wlo_bx_xy_hi);

      amrex::ParallelFor(wlo_bx_xy,
        [bct_klo,dom_lo,p_bc_w_g,p_bc_epf,ep_w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_klo(i,j,dom_lo[2]-1,1);
        const int bct = bct_klo(i,j,dom_lo[2]-1,0);
        if (bct == BCList::minf) ep_w(i,j,k) = p_bc_w_g[bcv] * p_bc_epf[bcv];
      });
    }

    if (nup > 0)
    {
      // Create InVects for following Boxes
      IntVect whi_bx_xy_lo(wbx_lo);
      IntVect whi_bx_xy_hi(wbx_hi);

      // Fix lo and hi limits
      whi_bx_xy_lo[2] = dom_hi[2]+1;
      whi_bx_xy_hi[2] = dom_hi[2]+1;

      const Box whi_bx_xy(whi_bx_xy_lo, whi_bx_xy_hi);

      amrex::ParallelFor(whi_bx_xy,
        [bct_khi,dom_hi,p_bc_w_g,p_bc_epf,ep_w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int bcv = bct_khi(i,j,dom_hi[2]+1,1);
        const int bct = bct_khi(i,j,dom_hi[2]+1,0);
        if (bct == BCList::minf) ep_w(i,j,k) = p_bc_w_g[bcv] * p_bc_epf[bcv];
      });
    }
  }
}


void
mfix::fillpatch_mac (Vector< MultiFab* > const& ep_u_mac,
                     Vector< MultiFab* > const& ep_v_mac,
                     Vector< MultiFab* > const& ep_w_mac)
{
   Real fake_time = 0;

   for (int lev = 0; lev < nlev(); ++lev)
   {
      const int minf = BCList::minf;
      const int cover = BCList::cover;

      Array<MultiFab*, AMREX_SPACEDIM> u_fine;
      AMREX_D_TERM(u_fine[0] = ep_u_mac[lev];,
                   u_fine[1] = ep_v_mac[lev];,
                   u_fine[2] = ep_w_mac[lev];);


      PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > fphysbc_x
        (geom[lev], bcs().get_velocity_bcrec(), MFIXScalarFill{minf, cover, maxLevel(),
            /*stride=*/0, m_boundary_conditions.bc_u_g().data(), /*covered_val=*/0.,
            bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
            bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
            bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
            });

      PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > fphysbc_y
        (geom[lev], bcs().get_velocity_bcrec(), MFIXScalarFill{minf, cover, maxLevel(),
            /*stride=*/0, m_boundary_conditions.bc_v_g().data(), /*covered_val=*/0.,
            bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
            bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
            bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
            });

      PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > fphysbc_z
        (geom[lev], bcs().get_velocity_bcrec(), MFIXScalarFill{minf, cover, maxLevel(),
            /*stride=*/0, m_boundary_conditions.bc_w_g().data(), /*covered_val=*/0.,
            bc_list.bc_ilo[lev]->array(), bc_list.bc_ihi[lev]->array(),
            bc_list.bc_jlo[lev]->array(), bc_list.bc_jhi[lev]->array(),
            bc_list.bc_klo[lev]->array(), bc_list.bc_khi[lev]->array()
            });

      Array<PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill>>,AMREX_SPACEDIM>
         fbndyFuncArr = {AMREX_D_DECL(fphysbc_x,fphysbc_y,fphysbc_z)};

      if (lev == 0) {
          for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
          {
             amrex::FillPatchSingleLevel(*u_fine[idim], fake_time,
                                         {u_fine[idim]}, {fake_time}, 0, 0, 1,
                                         geom[lev], fbndyFuncArr[idim], idim);
          }
      }
      else
      {
         Array<MultiFab*, AMREX_SPACEDIM> u_crse;
         AMREX_D_TERM(u_crse[0] = ep_u_mac[lev-1];,
                      u_crse[1] = ep_v_mac[lev-1];,
                      u_crse[2] = ep_w_mac[lev-1];);

         // Divergence preserving interp
         Interpolater* mapper = &face_divfree_interp;

         const Array<Vector<BCRec>,AMREX_SPACEDIM> bcrecArr = {AMREX_D_DECL(bcs().get_velocity_bcrec(),
                                                                            bcs().get_velocity_bcrec(),
                                                                            bcs().get_velocity_bcrec())};

         PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > cphysbc_x
           (geom[lev-1], bcs().get_velocity_bcrec(), MFIXScalarFill{minf, cover, maxLevel(),
               /*stride=*/0, m_boundary_conditions.bc_u_g().data(), /*covered_val=*/0.,
               bc_list.bc_ilo[lev-1]->array(), bc_list.bc_ihi[lev-1]->array(),
               bc_list.bc_jlo[lev-1]->array(), bc_list.bc_jhi[lev-1]->array(),
               bc_list.bc_klo[lev-1]->array(), bc_list.bc_khi[lev-1]->array()
               });


         PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > cphysbc_y
           (geom[lev-1], bcs().get_velocity_bcrec(), MFIXScalarFill{minf, cover, maxLevel(),
               /*stride=*/0, m_boundary_conditions.bc_v_g().data(), /*covered_val=*/0.,
               bc_list.bc_ilo[lev-1]->array(), bc_list.bc_ihi[lev-1]->array(),
               bc_list.bc_jlo[lev-1]->array(), bc_list.bc_jhi[lev-1]->array(),
               bc_list.bc_klo[lev-1]->array(), bc_list.bc_khi[lev-1]->array()
               });

          PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > cphysbc_z
           (geom[lev-1], bcs().get_velocity_bcrec(), MFIXScalarFill{minf, cover, maxLevel(),
               /*stride=*/0, m_boundary_conditions.bc_w_g().data(), /*covered_val=*/0.,
               bc_list.bc_ilo[lev-1]->array(), bc_list.bc_ihi[lev-1]->array(),
               bc_list.bc_jlo[lev-1]->array(), bc_list.bc_jhi[lev-1]->array(),
               bc_list.bc_klo[lev-1]->array(), bc_list.bc_khi[lev-1]->array()
               });

         Array<PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill>>,AMREX_SPACEDIM>
            cbndyFuncArr = {AMREX_D_DECL(cphysbc_x,cphysbc_y,cphysbc_z)};

         Array<int, AMREX_SPACEDIM> idx = {AMREX_D_DECL(0,1,2)};
         int ngmac = 1; //nghost_mac();
         amrex::FillPatchTwoLevels(u_fine, IntVect(ngmac), fake_time,
                                   {u_crse}, {fake_time},
                                   {u_fine}, {fake_time},
                                   0, 0, 1, geom[lev-1], geom[lev],
                                   cbndyFuncArr, idx, fbndyFuncArr, idx,
                                   refRatio(lev-1), mapper, bcrecArr, idx);
      }
   }
}
