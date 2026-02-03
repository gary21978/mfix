#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>

#include <mfix_deposition_op.H>

using namespace amrex;

void
MFIXDepositionOp::
ApplyBC (int const a_lev,
         MultiFab* a_MF)
{
  BL_PROFILE("MFIXDepositionOp::ApplyBC()");

  Box domain(m_pc->Geom(a_lev).Domain());

  for (MFIter mfi(*a_MF, false); mfi.isValid(); ++mfi) {

    const IntVect sbx_lo( ((*a_MF)[mfi]).loVect());
    const IntVect sbx_hi( ((*a_MF)[mfi]).hiVect());

    // Extract the lower and upper boundaries of Box sbx and Domain
    //const IntVect sbx_lo(sbx.loVect()), sbx_hi(sbx.hiVect());
    const amrex::Dim3& dom_lo = amrex::lbound(domain);
    const amrex::Dim3& dom_hi = amrex::ubound(domain);

    // Create a 2D Box collapsing sbx on x-direction
    IntVect sbx_yz_hi(sbx_hi);
    sbx_yz_hi[0] = sbx_lo[0];
    const Box sbx_yz(sbx_lo, sbx_yz_hi);

    // Create a 2D Box collapsing sbx on y-direction
    IntVect sbx_xz_hi(sbx_hi);
    sbx_xz_hi[1] = sbx_lo[1];
    const Box sbx_xz(sbx_lo, sbx_xz_hi);

    // Create a 2D Box collapsing sbx on z-direction
    IntVect sbx_xy_hi(sbx_hi);
    sbx_xy_hi[2] = sbx_lo[2];
    const Box sbx_xy(sbx_lo, sbx_xy_hi);

    Array4<int> const& bc_ilo_type = m_bc_list.bc_ilo[a_lev]->array();
    Array4<int> const& bc_ihi_type = m_bc_list.bc_ihi[a_lev]->array();
    Array4<int> const& bc_jlo_type = m_bc_list.bc_jlo[a_lev]->array();
    Array4<int> const& bc_jhi_type = m_bc_list.bc_jhi[a_lev]->array();
    Array4<int> const& bc_klo_type = m_bc_list.bc_klo[a_lev]->array();
    Array4<int> const& bc_khi_type = m_bc_list.bc_khi[a_lev]->array();

    const int ncomp = a_MF->nComp();

    Array4<Real> const& vol = a_MF->array(mfi);

    if(sbx_lo[0] < dom_lo.x || sbx_hi[0] > dom_hi.x ||
       sbx_lo[1] < dom_lo.y || sbx_hi[1] > dom_hi.y ||
       sbx_lo[2] < dom_lo.z || sbx_hi[2] > dom_hi.z) {

      const int sbx_yz_npoints = sbx_yz.numPts();
      const auto sbx_yz_lo = amrex::lbound(sbx_yz);
      const auto sbx_yz_len = amrex::length(sbx_yz);

      const int sbx_xz_npoints = sbx_xz.numPts();
      const auto sbx_xz_lo = amrex::lbound(sbx_xz);
      const auto sbx_xz_len = amrex::length(sbx_xz);

      const int sbx_xy_npoints = sbx_xy.numPts();
      const auto sbx_xy_lo = amrex::lbound(sbx_xy);
      const auto sbx_xy_len = amrex::length(sbx_xy);

      const int npoints = amrex::max(sbx_yz_npoints,sbx_xz_npoints,sbx_xy_npoints);

      ParallelFor(npoints, [sbx_xy_npoints,sbx_yz_npoints,sbx_xz_npoints,
          sbx_xy_len,sbx_yz_len,sbx_xz_len,sbx_xy_lo,sbx_yz_lo,sbx_xz_lo,
          dom_lo,dom_hi,bc_ilo_type,bc_ihi_type,bc_jlo_type,bc_jhi_type,
          bc_klo_type,bc_khi_type,vol,ncomp,sbx_hi,sbx_lo]
        AMREX_GPU_DEVICE (int idx) noexcept
      {
        if(idx < sbx_yz_npoints)
        {
          int k = idx / (sbx_yz_len.y);
          int j = idx - k*(sbx_yz_len.y);

          j += sbx_yz_lo.y;
          k += sbx_yz_lo.z;

          if(sbx_lo[0] < dom_lo.x) {
            int i = dom_lo.x;

            const int bct_ilo = bc_ilo_type(i-1,j,k,0);
            if(bct_ilo == BCList::pinf || bct_ilo == BCList::minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,  j,k,n) += vol(i-1,j,k,n);
                vol(i-1,j,k,n) = 0;
              }
            }
          }
          if(sbx_hi[0] > dom_hi.x) {
            int i = dom_hi.x;

            const int bct_ihi = bc_ihi_type(i+1,j,k,0);
            if(bct_ihi == BCList::pinf || bct_ihi == BCList::minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,  j,k,n) += vol(i+1,j,k,n);
                vol(i+1,j,k,n) = 0;
              }
            }
          }
        }

        if(idx < sbx_xz_npoints)
        {
          int k = idx / (sbx_xz_len.x);
          int i = idx - k*(sbx_xz_len.x);

          i += sbx_xz_lo.x;
          k += sbx_xz_lo.z;

          if(sbx_lo[1] < dom_lo.y) {
            int j = dom_lo.y;

            const int bct_jlo = bc_jlo_type(i,j-1,k,0);
            if(bct_jlo == BCList::pinf || bct_jlo == BCList::minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,j  ,k,n) += vol(i,j-1,k,n);
                vol(i,j-1,k,n) = 0;
              }
            }
          }
          if(sbx_hi[1] > dom_hi.y) {
            int j = dom_hi.y;

            const int bct_jhi = bc_jhi_type(i,j+1,k,0);
            if(bct_jhi == BCList::pinf || bct_jhi == BCList::minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,j  ,k,n) += vol(i,j+1,k,n);
                vol(i,j+1,k,n) = 0;
              }
            }
          }
        }

        if(idx < sbx_xy_npoints)
        {
          int j = idx / (sbx_xy_len.x);
          int i = idx - j*(sbx_xy_len.x);

          i += sbx_xy_lo.x;
          j += sbx_xy_lo.y;

          if(sbx_lo[2] < dom_lo.z) {
            int k = dom_lo.z;

            const int bct_klo = bc_klo_type(i,j,k-1,0);
            if(bct_klo == BCList::pinf || bct_klo == BCList::minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,j,k  ,n) += vol(i,j,k-1,n);
                vol(i,j,k-1,n) = 0;
              }
            }
          }
          if(sbx_hi[2] > dom_hi.z) {
            int k = dom_hi.z;

            const int bct_khi = bc_khi_type(i,j,k+1,0);
            if(bct_khi == BCList::pinf || bct_khi == BCList::minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,j,k  ,n) += vol(i,j,k+1,n);
                vol(i,j,k+1,n) = 0;
              }
            }
          }
        }
      });
    }
  }

}
