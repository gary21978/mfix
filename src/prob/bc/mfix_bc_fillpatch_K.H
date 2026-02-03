#ifndef MFIX_FILLPATCH_BC_H_
#define MFIX_FILLPATCH_BC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>


struct MFIXForFill
{
  int probtype;

  AMREX_GPU_HOST
  constexpr MFIXForFill (int a_probtype) : probtype(a_probtype) {}

  AMREX_GPU_DEVICE
  void operator() (const amrex::IntVect& /*iv*/, amrex::Array4<amrex::Real> const& /*vel*/,
                   const int /*dcomp*/, const int /*numcomp*/,
                   amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                   const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                   const int /*orig_comp*/) const
  {}
};


struct MFIXVelFill
{
  const int minf;
  const int cover;
  const int max_level;

  const amrex::Real* bc_u_g;
  const amrex::Real* bc_v_g;
  const amrex::Real* bc_w_g;

  amrex::Array4<const int> bct_ilo;
  amrex::Array4<const int> bct_ihi;
  amrex::Array4<const int> bct_jlo;
  amrex::Array4<const int> bct_jhi;
  amrex::Array4<const int> bct_klo;
  amrex::Array4<const int> bct_khi;

  AMREX_GPU_HOST
  constexpr MFIXVelFill (const int a_minf,
                         const int a_cover,
                         const int a_max_level,
                         const amrex::Real* a_bc_u_g,
                         const amrex::Real* a_bc_v_g,
                         const amrex::Real* a_bc_w_g,
                         amrex::Array4<const int> const& a_bct_ilo,
                         amrex::Array4<const int> const& a_bct_ihi,
                         amrex::Array4<const int> const& a_bct_jlo,
                         amrex::Array4<const int> const& a_bct_jhi,
                         amrex::Array4<const int> const& a_bct_klo,
                         amrex::Array4<const int> const& a_bct_khi)
    : minf(a_minf), cover(a_cover), max_level(a_max_level),
      bc_u_g(a_bc_u_g), bc_v_g(a_bc_v_g), bc_w_g(a_bc_w_g),
      bct_ilo(a_bct_ilo), bct_ihi(a_bct_ihi),
      bct_jlo(a_bct_jlo), bct_jhi(a_bct_jhi),
      bct_klo(a_bct_klo), bct_khi(a_bct_khi)
  {}

  AMREX_GPU_DEVICE
  void operator() (const amrex::IntVect& iv,
                   amrex::Array4<amrex::Real> const& vel,
                   const int /*dcomp*/,
                   const int /*numcomp*/,
                   amrex::GeometryData const& geom,
                   const amrex::Real /*time*/,
                   const amrex::BCRec* /*bcr*/,
                   const int /*bcomp*/,
                   const int /*orig_comp*/) const
  {
    using namespace amrex;

    // We are hard-wiring this fillpatch routine to define the Dirichlet values
    //    at the faces (not the ghost cell center)
    const int i = iv[0];
    const int j = iv[1];
    const int k = iv[2];

    // const int bct_test = bct_ilo(-1,0,0,0);
    // printf("(%d,%d,%d) bct %d,  bcv %d\n", i,j,k, bct, bcv);

    const Box& domain_box = geom.Domain();

    if (i < domain_box.smallEnd(0)) {

      const int bct = bct_ilo(domain_box.smallEnd(0)-1,j,k,0);
      const int bcv = bct_ilo(domain_box.smallEnd(0)-1,j,k,1);

      if (bct == minf) {
        vel(i,j,k,0) = bc_u_g[bcv];
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = 0;
      } else if (max_level > 0 && bct == cover) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = 0;
      }

    } else if (i > domain_box.bigEnd(0)) {

      const int bct = bct_ihi(domain_box.bigEnd(0)+1,j,k,0);
      const int bcv = bct_ihi(domain_box.bigEnd(0)+1,j,k,1);

      if (bct == minf) {
        vel(i,j,k,0) = bc_u_g[bcv];
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = 0;
      } else if (max_level > 0 && bct == cover) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = 0;
      }

    }

    if (j < domain_box.smallEnd(1)) {

      const int bct = bct_jlo(i,domain_box.smallEnd(1)-1,k,0);
      const int bcv = bct_jlo(i,domain_box.smallEnd(1)-1,k,1);

      if (bct == minf) {
          vel(i,j,k,0) = 0;
          vel(i,j,k,1) = bc_v_g[bcv];
          vel(i,j,k,2) = 0;
      } else if (max_level > 0 && bct == cover) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = 0;
      }

    } else if (j > domain_box.bigEnd(1)) {

      const int bcv = bct_jhi(i,domain_box.bigEnd(1)+1,k,1);
      const int bct = bct_jhi(i,domain_box.bigEnd(1)+1,k,0);

      if (bct == minf) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = bc_v_g[bcv];
        vel(i,j,k,2) = 0;
      } else if (max_level > 0 && bct == cover) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = 0;
      }

    }

    if (k < domain_box.smallEnd(2)) {

      const int bcv = bct_klo(i,j,domain_box.smallEnd(2)-1,1);
      const int bct = bct_klo(i,j,domain_box.smallEnd(2)-1,0);

      if (bct == minf) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = bc_w_g[bcv];
      } else if (max_level > 0 && bct == cover) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = 0;
      }

    } else if (k > domain_box.bigEnd(2)) {

      const int bcv = bct_khi(i,j,domain_box.bigEnd(2)+1,1);
      const int bct = bct_khi(i,j,domain_box.bigEnd(2)+1,0);

      if (bct == minf) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = bc_w_g[bcv];
      } else if (max_level > 0 && bct == cover) {
        vel(i,j,k,0) = 0;
        vel(i,j,k,1) = 0;
        vel(i,j,k,2) = 0;
      }
    }
  }
};


struct MFIXScalarFill
{
  const int minf;
  const int cover;
  const int max_level;

  const int stride;

  const amrex::Real* bc_scalar;

  amrex::Real const m_covered_val;

  amrex::Array4<const int> bct_ilo;
  amrex::Array4<const int> bct_ihi;
  amrex::Array4<const int> bct_jlo;
  amrex::Array4<const int> bct_jhi;
  amrex::Array4<const int> bct_klo;
  amrex::Array4<const int> bct_khi;

  AMREX_GPU_HOST
  constexpr MFIXScalarFill ( int const a_minf,
                             int const a_cover,
                             int const a_max_level,
                             int const a_stride,
                             const amrex::Real* a_bc_scalar,
                             amrex::Real const a_covered_val,
                             amrex::Array4<const int> const& a_bct_ilo,
                             amrex::Array4<const int> const& a_bct_ihi,
                             amrex::Array4<const int> const& a_bct_jlo,
                             amrex::Array4<const int> const& a_bct_jhi,
                             amrex::Array4<const int> const& a_bct_klo,
                             amrex::Array4<const int> const& a_bct_khi)
    : minf(a_minf), cover(a_cover), max_level(a_max_level)
    , stride(a_stride), bc_scalar(a_bc_scalar)
    , m_covered_val(a_covered_val)
    , bct_ilo(a_bct_ilo), bct_ihi(a_bct_ihi)
    , bct_jlo(a_bct_jlo), bct_jhi(a_bct_jhi)
    , bct_klo(a_bct_klo), bct_khi(a_bct_khi)
  {}

  AMREX_GPU_DEVICE
  void operator() (const amrex::IntVect& iv,
                   amrex::Array4<amrex::Real> const& scalar,
                   const int/*dcomp*/,
                   const int numcomp,
                   amrex::GeometryData const& geom,
                   const amrex::Real /*time*/,
                   const amrex::BCRec* /*bcr*/,
                   const int /*bcomp*/,
                   const int /*orig_comp*/) const
  {
    using namespace amrex;

    // We are hard-wiring this fillpatch routine to define the Dirichlet values
    //    at the faces (not the ghost cell center)
    const int i = iv[0];
    const int j = iv[1];
    const int k = iv[2];

    const Box& domain_box = geom.Domain();

    for(int n(0); n<numcomp; n++) {

      if (i < domain_box.smallEnd(0)) {

        const int bct = bct_ilo(domain_box.smallEnd(0)-1,j,k,0);
        const int bcv = bct_ilo(domain_box.smallEnd(0)-1,j,k,1);

        if (bct == minf) {
          scalar(i,j,k,n) = bc_scalar[bcv + n*stride];
        } else if (max_level > 0 && bct == cover) {
          scalar(i,j,k,n) = m_covered_val;
        }

      } else if (i > domain_box.bigEnd(0)) {

        const int bct = bct_ihi(domain_box.bigEnd(0)+1,j,k,0);
        const int bcv = bct_ihi(domain_box.bigEnd(0)+1,j,k,1);

        if (bct == minf) {
          scalar(i,j,k,n) = bc_scalar[bcv + n*stride];
        } else if (max_level > 0 && bct == cover) {
          scalar(i,j,k,n) = m_covered_val;
        }
      }

      if (j < domain_box.smallEnd(1)) {

        const int bct = bct_jlo(i,domain_box.smallEnd(1)-1,k,0);
        const int bcv = bct_jlo(i,domain_box.smallEnd(1)-1,k,1);

        if (bct == minf) {
          scalar(i,j,k,n) = bc_scalar[bcv + n*stride];
        } else if (max_level > 0 && bct == cover) {
          scalar(i,j,k,n) = m_covered_val;
        }

      } else if (j > domain_box.bigEnd(1)) {

        const int bcv = bct_jhi(i,domain_box.bigEnd(1)+1,k,1);
        const int bct = bct_jhi(i,domain_box.bigEnd(1)+1,k,0);

        if (bct == minf) {
          scalar(i,j,k,n) = bc_scalar[bcv + n*stride];
        } else if (max_level > 0 && bct == cover) {
          scalar(i,j,k,n) = m_covered_val;
        }
      }

      if (k < domain_box.smallEnd(2)) {

        const int bcv = bct_klo(i,j,domain_box.smallEnd(2)-1,1);
        const int bct = bct_klo(i,j,domain_box.smallEnd(2)-1,0);

        if (bct == minf) {
          scalar(i,j,k,n) = bc_scalar[bcv + n*stride];
        } else if (max_level > 0 && bct == cover) {
          scalar(i,j,k,n) = m_covered_val;
        }

      } else if (k > domain_box.bigEnd(2)) {

        const int bcv = bct_khi(i,j,domain_box.bigEnd(2)+1,1);
        const int bct = bct_khi(i,j,domain_box.bigEnd(2)+1,0);

        if (bct == minf) {
          scalar(i,j,k,n) = bc_scalar[bcv + n*stride];
        } else if (max_level > 0 && bct == cover) {
          scalar(i,j,k,n) = m_covered_val;
        }
      }

    } // ncomp
  } // operator
};
#endif
