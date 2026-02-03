#include <AMReX.H>
#include <mfix.H>
#include <mfix_mf_helpers.H>

using namespace amrex;

void mfix::MFIX_CalcAvgSolidsVel (Vector< Array<MultiFab*,3> >& vel_s,
                                  const bool do_deposition)
{
  BL_PROFILE("MFIX_CalcAvgSolidsProp()");

  const amrex::FabArray<EBCellFlagFab>* flags = nullptr;
  //Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac; UNUSED

  for (int lev = 0; lev < nlev(); lev ++ )
  {

    const amrex::BoxArray&            pba = pc->ParticleBoxArray(lev);
    const amrex::DistributionMapping& pdm = pc->ParticleDistributionMap(lev);

    // Use level 0 to define the EB factory. If we are not on level 0
    // then create a copy of the coarse factory to use.
   if (lev == 0) {
      flags  = &(m_eb->particle_factory()[lev]->getMultiEBCellFlagFab());
      // areafrac = particle_ebfactory[lev]->getAreaFrac(); UNUSED

    } else {

      Vector<int> ngrow = {1,1,1};
      EBFArrayBoxFactory* crse_factory;

      crse_factory = (amrex::makeEBFabFactory(geom[0], pba, pdm, ngrow, EBSupport::volume)).release();

      flags  = &(crse_factory->getMultiEBCellFlagFab());
      delete crse_factory;
    }

    // We deposit vel*pmass and pmass. Later on we need to divide out
    // the accumulated mass to get the mass-averaged grid velocities.
   if ( do_deposition ) {
     pc->MFIX_PC_SolidsVelocityDeposition(lev, vel_s[lev], flags);
   }
  }


  {
    // The deposition occurred on level 0, thus the next few operations
    // only need to be carried out on level 0.
    int lev(0);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    vel_s[lev][0]->SumBoundary(geom[lev].periodicity());
    vel_s[lev][1]->SumBoundary(geom[lev].periodicity());
    vel_s[lev][2]->SumBoundary(geom[lev].periodicity());

    constexpr Real tolerance = std::numeric_limits<double>::min();

    Box domain(geom[lev].Domain());

    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    Array4<const int> const& bct_ilo = bc_list.bc_ilo[lev]->array();
    Array4<const int> const& bct_ihi = bc_list.bc_ihi[lev]->array();
    Array4<const int> const& bct_jlo = bc_list.bc_jlo[lev]->array();
    Array4<const int> const& bct_jhi = bc_list.bc_jhi[lev]->array();
    Array4<const int> const& bct_klo = bc_list.bc_klo[lev]->array();
    Array4<const int> const& bct_khi = bc_list.bc_khi[lev]->array();

    for (MFIter mfi(*vel_s[lev][0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& xbx = mfi.growntilebox(IntVect(0,1,1));

      Array4<Real> const& u_s = vel_s[lev][0]->array(mfi);

      // A mix of regular and cut-cells.
      if ( (*flags)[mfi].getType(xbx) != FabType::covered) {
        amrex::ParallelFor(xbx, [u_s] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if ( u_s(i,j,k,1) > tolerance ) {
            u_s(i,j,k,0) /= u_s(i,j,k,1);
          } else {
            u_s(i,j,k,0) = 0.0;
          }
        });

      } else {
        amrex::ParallelFor(xbx, [u_s] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u_s(i,j,k,0) = 0.0; });
      }

      IntVect u_lo((*vel_s[lev][0])[mfi].loVect());
      IntVect u_hi((*vel_s[lev][0])[mfi].hiVect());

      const int nlft = amrex::max(0, dom_lo[0]-u_lo[0]);
      const int nrgt = amrex::max(0, u_hi[0]-dom_hi[0]);

      if (nlft > 0) {

        IntVect bx_yz_lo_lo_2D(u_lo), bx_yz_lo_hi_2D(u_hi);
        bx_yz_lo_hi_2D[0] = dom_lo[0];
        const Box bx_yz_lo_2D(bx_yz_lo_lo_2D, bx_yz_lo_hi_2D);

        amrex::ParallelFor(bx_yz_lo_2D, [bct_ilo,dom_lo,u_s]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          const int bct = bct_ilo(dom_lo[0]-1,j,k,0);
          if((bct == BCList::pinf) || (bct == BCList::minf)) u_s(i,j,k,0) = 0.0;
        });
      } // nlft


      if (nrgt > 0 ) {

        IntVect bx_yz_hi_lo_2D(u_lo), bx_yz_hi_hi_2D(u_hi);
        bx_yz_hi_lo_2D[0] = dom_hi[0]+1;
        const Box bx_yz_hi_2D(bx_yz_hi_lo_2D, bx_yz_hi_hi_2D);

        amrex::ParallelFor(bx_yz_hi_2D, [bct_ihi,dom_hi,u_s]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          const int bct = bct_ihi(dom_hi[0]+1,j,k,0);
          if((bct == BCList::pinf) || (bct == BCList::minf)) u_s(i,j,k,0) = 0.0;
        });
      } // nrgt
    }// end vel[0] MFIter loop


    for (MFIter mfi(*vel_s[lev][1], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& ybx = mfi.growntilebox(IntVect(1,0,1));

      Array4<Real> const& v_s = vel_s[lev][1]->array(mfi);

      // A mix of regular and cut-cells.
      if ( (*flags)[mfi].getType(ybx) != FabType::covered) {
        amrex::ParallelFor(ybx, [v_s] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if ( v_s(i,j,k,1) > tolerance ) {
            v_s(i,j,k,0) /= v_s(i,j,k,1);
          } else {
            v_s(i,j,k,0) = 0.0;
          }
        });

      } else {
        amrex::ParallelFor(ybx, [v_s] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { v_s(i,j,k,0) = 0.0; });
      }

      IntVect v_lo((*vel_s[lev][1])[mfi].loVect());
      IntVect v_hi((*vel_s[lev][1])[mfi].hiVect());

      const int nbot = amrex::max(0, dom_lo[1]-v_lo[1]);
      const int ntop = amrex::max(0, v_hi[1]-dom_hi[1]);

      if (nbot > 0) {

        IntVect bx_xz_lo_lo_2D(v_lo), bx_xz_lo_hi_2D(v_hi);
        bx_xz_lo_hi_2D[1] = dom_lo[1];
        const Box bx_xz_lo_2D(bx_xz_lo_lo_2D, bx_xz_lo_hi_2D);

        amrex::ParallelFor(bx_xz_lo_2D, [bct_jlo,dom_lo,v_s]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          const int bct = bct_jlo(i,dom_lo[1]-1,k,0);
          if((bct == BCList::pinf) || (bct == BCList::minf)) v_s(i,j,k,0) = 0.0;
        });
      } // nbot

      if (ntop > 0 ) {

        IntVect bx_xz_hi_lo_2D(v_lo), bx_xz_hi_hi_2D(v_hi);
        bx_xz_hi_lo_2D[1] = dom_hi[1]+1;
        const Box bx_xz_hi_2D(bx_xz_hi_lo_2D, bx_xz_hi_hi_2D);

        amrex::ParallelFor(bx_xz_hi_2D, [bct_jhi,dom_hi,v_s]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          const int bct = bct_jhi(i,dom_hi[1]+1,k,0);
          if((bct == BCList::pinf) || (bct == BCList::minf)) v_s(i,j,k,0) = 0.0;
        });
      } // ntop
    }// end vel[1] MFIter loop

    for (MFIter mfi(*vel_s[lev][2], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& zbx = mfi.growntilebox(IntVect(1,1,0));

      Array4<Real> const& w_s = vel_s[lev][2]->array(mfi);

      // A mix of regular and cut-cells.
      if ( (*flags)[mfi].getType(zbx) != FabType::covered) {

        amrex::ParallelFor(zbx, [w_s] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if ( w_s(i,j,k,1) > tolerance ) {
            w_s(i,j,k,0) /= w_s(i,j,k,1);
          } else {
            w_s(i,j,k,0) = 0.0;
          }
        });

      } else {
        amrex::ParallelFor(zbx, [w_s] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { w_s(i,j,k,0) = 0.0; });
      }
      IntVect w_lo((*vel_s[lev][2])[mfi].loVect());
      IntVect w_hi((*vel_s[lev][2])[mfi].hiVect());

      const int ndwn = amrex::max(0, dom_lo[2]-w_lo[2]);
      const int nup  = amrex::max(0, w_hi[2]-dom_hi[2]);

      if (ndwn > 0) {

        IntVect bx_xy_lo_lo_2D(w_lo), bx_xy_lo_hi_2D(w_hi);
        bx_xy_lo_hi_2D[2] = dom_lo[2];
        const Box bx_xy_lo_2D(bx_xy_lo_lo_2D, bx_xy_lo_hi_2D);

        amrex::ParallelFor(bx_xy_lo_2D, [bct_klo,dom_lo,w_s]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          const int bct = bct_klo(i,j,dom_lo[2]-1,0);
          if((bct == BCList::pinf) || (bct == BCList::minf)) w_s(i,j,k,0) = 0.0;

        });
      } // ndwn

      if (nup  > 0 ) {

        IntVect bx_xy_hi_lo_2D(w_lo), bx_xy_hi_hi_2D(w_hi);
        bx_xy_hi_lo_2D[2] = dom_hi[2]+1;
        const Box bx_xy_hi_2D(bx_xy_hi_lo_2D, bx_xy_hi_hi_2D);

        amrex::ParallelFor(bx_xy_hi_2D, [bct_khi,dom_hi,w_s]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          const int bct = bct_khi(i,j,dom_hi[2]+1,0);
          if((bct == BCList::pinf) || (bct == BCList::minf)) w_s(i,j,k,0) = 0.0;

        });
      } // nup
    }// end vel[2] MFIter loop

    vel_s[lev][0]->FillBoundary(geom[lev].periodicity());
    vel_s[lev][1]->FillBoundary(geom[lev].periodicity());
    vel_s[lev][2]->FillBoundary(geom[lev].periodicity());
  }

  // Interpolate from the coarse grid to define the fine grid arrays
  if (nlev() > 1) {

    amrex::Abort(" HACK: For right now the PIC implementation only works"
                 " with a single level. If we want to make sure it works"
                 " for more, then we need to fix this.");
  }

}
