#include <mfix_bc.H>

using namespace amrex;

void
MFIXBoundaryConditions::set_gradp_bcs (const Box& /*bx*/,
                                       const int lev,
                                       FArrayBox& gp_fab,
                                       Box& domain)
{
  // Extract the lower and upper boundaries of Domain
  const IntVect dom_lo(domain.loVect()), dom_hi(domain.hiVect());

  Array4<Real> const& gp_arr = gp_fab.array();
  const IntVect gp_lo(gp_fab.loVect()), gp_hi(gp_fab.hiVect());

  Array4<int> const& bct_ilo = m_bc_list.bc_ilo[lev]->array();
  Array4<int> const& bct_ihi = m_bc_list.bc_ihi[lev]->array();
  Array4<int> const& bct_jlo = m_bc_list.bc_jlo[lev]->array();
  Array4<int> const& bct_jhi = m_bc_list.bc_jhi[lev]->array();
  Array4<int> const& bct_klo = m_bc_list.bc_klo[lev]->array();
  Array4<int> const& bct_khi = m_bc_list.bc_khi[lev]->array();

  if(gp_lo[0] <= dom_lo[0])
  {
    IntVect bx_lo_yz_lo(gp_lo), bx_lo_yz_hi(gp_hi);
    bx_lo_yz_lo[0] = dom_lo[0]-1;
    bx_lo_yz_hi[0] = dom_lo[0]-1;
    const Box bx_yz_lo(bx_lo_yz_lo, bx_lo_yz_hi);

    amrex::ParallelFor(bx_yz_lo, [bct_ilo,dom_lo,gp_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ilo(dom_lo[0]-1,j,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        gp_arr(i,j,k,0) = gp_arr(i+1,j,k,0);
        gp_arr(i,j,k,1) = gp_arr(i+1,j,k,1);
        gp_arr(i,j,k,2) = gp_arr(i+1,j,k,2);
      }
      else if(bct == BCList::minf)
      {
        gp_arr(i,j,k,0) = gp_arr(i+1,j,k,0);
        gp_arr(i,j,k,1) = 0;
        gp_arr(i,j,k,2) = 0;
      }
    });
  }

  if(gp_hi[0] >= dom_hi[0]+1)
  {
    IntVect bx_hi_yz_lo(gp_lo), bx_hi_yz_hi(gp_hi);
    bx_hi_yz_lo[0] = dom_hi[0]+1;
    bx_hi_yz_hi[0] = dom_hi[0]+1;
    const Box bx_yz_hi(bx_hi_yz_lo, bx_hi_yz_hi);

    amrex::ParallelFor(bx_yz_hi, [bct_ihi,dom_hi,gp_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_ihi(dom_hi[0]+1,j,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        gp_arr(i,j,k,0) = gp_arr(i-1,j,k,0);
        gp_arr(i,j,k,1) = gp_arr(i-1,j,k,1);
        gp_arr(i,j,k,2) = gp_arr(i-1,j,k,2);
      }
      else if(bct == BCList::minf)
      {
        gp_arr(i,j,k,0) = gp_arr(i-1,j,k,0);
        gp_arr(i,j,k,1) = 0;
        gp_arr(i,j,k,2) = 0;
      }
    });
  }

  if(gp_lo[1] <= dom_lo[1])
  {
    IntVect bx_lo_xz_lo(gp_lo), bx_lo_xz_hi(gp_hi);
    bx_lo_xz_lo[1] = dom_lo[1]-1;
    bx_lo_xz_hi[1] = dom_lo[1]-1;
    const Box bx_xz_lo(bx_lo_xz_lo, bx_lo_xz_hi);

    amrex::ParallelFor(bx_xz_lo, [bct_jlo,dom_lo,gp_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jlo(i,dom_lo[1]-1,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        gp_arr(i,j,k,0) = gp_arr(i,j+1,k,0);
        gp_arr(i,j,k,1) = gp_arr(i,j+1,k,1);
        gp_arr(i,j,k,2) = gp_arr(i,j+1,k,2);
      }
      else if(bct == BCList::minf)
      {
        gp_arr(i,j,k,0) = 0;
        gp_arr(i,j,k,1) = gp_arr(i,j+1,k,1);
        gp_arr(i,j,k,2) = 0;
      }
    });
  }

  if(gp_hi[1] >= dom_hi[1]+1)
  {
    IntVect bx_hi_xz_lo(gp_lo), bx_hi_xz_hi(gp_hi);
    bx_hi_xz_lo[1] = dom_hi[1]+1;
    bx_hi_xz_hi[1] = dom_hi[1]+1;
    const Box bx_xz_hi(bx_hi_xz_lo, bx_hi_xz_hi);

    amrex::ParallelFor(bx_xz_hi, [bct_jhi,dom_hi,gp_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_jhi(i,dom_hi[1]+1,k,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        gp_arr(i,j,k,0) = gp_arr(i,j-1,k,0);
        gp_arr(i,j,k,1) = gp_arr(i,j-1,k,1);
        gp_arr(i,j,k,2) = gp_arr(i,j-1,k,2);
      }
      else if(bct == BCList::minf)
      {
        gp_arr(i,j,k,0) = 0;
        gp_arr(i,j,k,1) = gp_arr(i,j-1,k,1);
        gp_arr(i,j,k,2) = 0;
      }
    });
  }

  if(gp_lo[2] <= dom_lo[2])
  {
    IntVect bx_lo_xy_lo(gp_lo), bx_lo_xy_hi(gp_hi);
    bx_lo_xy_lo[2] = dom_lo[2]-1;
    bx_lo_xy_hi[2] = dom_lo[2]-1;
    const Box bx_xy_lo(bx_lo_xy_lo, bx_lo_xy_hi);

    amrex::ParallelFor(bx_xy_lo, [bct_klo,dom_lo,gp_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_klo(i,j,dom_lo[2]-1,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        gp_arr(i,j,k,0) = gp_arr(i,j,k+1,0);
        gp_arr(i,j,k,1) = gp_arr(i,j,k+1,1);
        gp_arr(i,j,k,2) = gp_arr(i,j,k+1,2);
      }
      else if(bct == BCList::minf)
      {
        gp_arr(i,j,k,0) = 0;
        gp_arr(i,j,k,1) = 0;
        gp_arr(i,j,k,2) = gp_arr(i,j,k+1,2);
      }
    });
  }

  if(gp_hi[2] >= dom_hi[2]+1)
  {
    IntVect bx_hi_xy_lo(gp_lo), bx_hi_xy_hi(gp_hi);
    bx_hi_xy_lo[2] = dom_hi[2]+1;
    bx_hi_xy_hi[2] = dom_hi[2]+1;
    const Box bx_xy_hi(bx_hi_xy_lo, bx_hi_xy_hi);

    amrex::ParallelFor(bx_xy_hi, [bct_khi,dom_hi,gp_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int bct = bct_khi(i,j,dom_hi[2]+1,0);

      if((bct == BCList::pinf) || (bct == BCList::pout))
      {
        gp_arr(i,j,k,0) = gp_arr(i,j,k-1,0);
        gp_arr(i,j,k,1) = gp_arr(i,j,k-1,1);
        gp_arr(i,j,k,2) = gp_arr(i,j,k-1,2);
      }
      else if(bct == BCList::minf)
      {
        gp_arr(i,j,k,0) = 0;
        gp_arr(i,j,k,1) = 0;
        gp_arr(i,j,k,2) = gp_arr(i,j,k-1,2);
      }
    });
  }
}
