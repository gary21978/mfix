#include <mfix_run_on.H>
#include <mfix_reporter.H>
#include <mfix_bc.H>
#include <mfix_fluid.H>

using namespace amrex;

void
MFIXBoundaryConditions::
set_bc_list (int const a_lev, int const a_nghost_bc)
{
  Real dx = m_geom[a_lev].CellSize(0);
  Real dy = m_geom[a_lev].CellSize(1);
  Real dz = m_geom[a_lev].CellSize(2);

  const GpuArray<Real, 3> plo = m_geom[a_lev].ProbLoArray();

  { // begin x direction

    const int init_x = m_geom[a_lev].isPeriodic(0) ? BCList::undefined : BCList::cover;

    Array4<int> const& bc_ilo_type = m_bc_list.bc_ilo[a_lev]->array();
    Array4<int> const& bc_ihi_type = m_bc_list.bc_ihi[a_lev]->array();

    Box domainx(m_geom[a_lev].Domain());

    domainx.grow(1,a_nghost_bc);  // Add ghost cells to y
    domainx.grow(2,a_nghost_bc);  // Add ghost cells to z

    { // x-lo side of the domain

      Box box_ilo = amrex::adjCellLo(domainx,0,1);

      IntVect ibx_lo(box_ilo.loVect());
      IntVect ibx_hi(box_ilo.hiVect());

      // Initialize x-lo domain extent.
      ParallelFor(box_ilo, [bc_ilo_type, init_x]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {bc_ilo_type(i,j,k,0) = init_x;});

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_xlo().size(); ++lc) {

        const int bcv  = bc_xlo(lc);
        const int type = m_bc[bcv].type;

        if (lc > 0){
          ibx_lo[1] = m_bc[bcv].Lo(1,plo[1],dy);
          ibx_lo[2] = m_bc[bcv].Lo(2,plo[2],dz);

          ibx_hi[1] = m_bc[bcv].Hi(1,plo[1],dy);
          ibx_hi[2] = m_bc[bcv].Hi(2,plo[2],dz);
        }

        const Box lo_box(ibx_lo, ibx_hi);
        ParallelFor(lo_box, [bc_ilo_type, type, bcv]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          bc_ilo_type(i,j,k,0) = type;
          bc_ilo_type(i,j,k,1) = bcv;
        });
      }
    } // end x-lo side of the domain

    { // x-hi side of the domain

      Box box_ihi = amrex::adjCellHi(domainx,0,1);

      IntVect ibx_lo(box_ihi.loVect());
      IntVect ibx_hi(box_ihi.hiVect());

      // Initialize x-lo domain extent.
      ParallelFor(box_ihi, [bc_ihi_type, init_x]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {bc_ihi_type(i,j,k,0) = init_x;});

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_xhi().size(); ++lc) {

        const int bcv  = bc_xhi(lc);
        const int type = m_bc[bcv].type;

        if (lc > 0){
          ibx_lo[1] = m_bc[bcv].Lo(1,plo[1],dy);
          ibx_lo[2] = m_bc[bcv].Lo(2,plo[2],dz);

          ibx_hi[1] = m_bc[bcv].Hi(1,plo[1],dy);
          ibx_hi[2] = m_bc[bcv].Hi(2,plo[2],dz);
        }

        const Box hi_box(ibx_lo, ibx_hi);
        ParallelFor(hi_box, [bc_ihi_type, type, bcv]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          bc_ihi_type(i,j,k,0) = type;
          bc_ihi_type(i,j,k,1) = bcv;
        });
      }

    } // end x-hi side of the domain
  } // end x-direction


  { // begin y-direction

    const int init_y = m_geom[a_lev].isPeriodic(1) ? BCList::undefined : BCList::cover;

    Array4<int> const& bc_jlo_type = m_bc_list.bc_jlo[a_lev]->array();
    Array4<int> const& bc_jhi_type = m_bc_list.bc_jhi[a_lev]->array();

    Box domainy(m_geom[a_lev].Domain());
    domainy.grow(0,a_nghost_bc);  // Add ghost cells to x
    domainy.grow(2,a_nghost_bc);  // Add ghost cells to z

    { // y-lo side of the domain

      Box box_jlo = amrex::adjCellLo(domainy,1,1);
      IntVect jbx_lo(box_jlo.loVect());
      IntVect jbx_hi(box_jlo.hiVect());

      // Initialize y-lo domain extent.
      ParallelFor(box_jlo, [bc_jlo_type, init_y]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {bc_jlo_type(i,j,k,0) = init_y;});

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_ylo().size(); ++lc) {

        const int bcv  = bc_ylo(lc);
        const int type = m_bc[bcv].type;

        if (lc > 0){
          jbx_lo[0] = m_bc[bcv].Lo(0,plo[0],dx);
          jbx_lo[2] = m_bc[bcv].Lo(2,plo[2],dz);

          jbx_hi[0] = m_bc[bcv].Hi(0,plo[0],dx);
          jbx_hi[2] = m_bc[bcv].Hi(2,plo[2],dz);
        }

        const Box lo_box(jbx_lo, jbx_hi);
        ParallelFor(lo_box, [bc_jlo_type, type, bcv]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          bc_jlo_type(i,j,k,0) = type;
          bc_jlo_type(i,j,k,1) = bcv;
        });
      }
    }// end y-lo side of the domain


    { // y-hi side of the domain

      Box box_jhi = amrex::adjCellHi(domainy,1,1);

      IntVect jbx_lo(box_jhi.loVect());
      IntVect jbx_hi(box_jhi.hiVect());

      // Initialize y-hi domain extent.
      ParallelFor(box_jhi, [bc_jhi_type, init_y]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {bc_jhi_type(i,j,k,0) = init_y;});

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_yhi().size(); ++lc) {

        const int bcv  = bc_yhi(lc);
        const int type = m_bc[bcv].type;

        if (lc > 0){
          jbx_lo[0] = m_bc[bcv].Lo(0,plo[0],dx);
          jbx_lo[2] = m_bc[bcv].Lo(2,plo[2],dz);

          jbx_hi[0] = m_bc[bcv].Hi(0,plo[0],dx);
          jbx_hi[2] = m_bc[bcv].Hi(2,plo[2],dz);
        }

        const Box hi_box(jbx_lo, jbx_hi);
        ParallelFor(hi_box, [bc_jhi_type, type, bcv]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          bc_jhi_type(i,j,k,0) = type;
          bc_jhi_type(i,j,k,1) = bcv;
        });
      }
    } // end y-hi side of the domain
  } // end y-direction


  { // begin z-direction

    const int init_z = m_geom[a_lev].isPeriodic(2) ? BCList::undefined : BCList::cover;

    Array4<int> const& bc_klo_type = m_bc_list.bc_klo[a_lev]->array();
    Array4<int> const& bc_khi_type = m_bc_list.bc_khi[a_lev]->array();

    Box domainz(m_geom[a_lev].Domain());
    domainz.grow(0,a_nghost_bc);  // Add ghost cells to x
    domainz.grow(1,a_nghost_bc);  // Add ghost cells to y

    { // z-lo side of the domain

      Box box_klo = amrex::adjCellLo(domainz,2,1);

      IntVect kbx_lo(box_klo.loVect());
      IntVect kbx_hi(box_klo.hiVect());

      // Initialize z-lo domain extent.
      ParallelFor(box_klo, [bc_klo_type, init_z]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {bc_klo_type(i,j,k,0) = init_z;});

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_zlo().size(); ++lc) {

        const int bcv  = bc_zlo(lc);
        const int type = m_bc[bcv].type;

        if (lc > 0){
          kbx_lo[0] = m_bc[bcv].Lo(0,plo[0],dx);
          kbx_lo[1] = m_bc[bcv].Lo(1,plo[1],dy);

          kbx_hi[0] = m_bc[bcv].Hi(0,plo[0],dx);
          kbx_hi[1] = m_bc[bcv].Hi(1,plo[1],dy);
        }

        const Box lo_box(kbx_lo, kbx_hi);
        ParallelFor(lo_box, [bc_klo_type, type, bcv]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          bc_klo_type(i,j,k,0) = type;
          bc_klo_type(i,j,k,1) = bcv;
        });
      }
    } // end z-lo side of the domain


    { // z-hi side of the domain

      Box box_khi = amrex::adjCellHi(domainz,2,1);

      IntVect kbx_lo(box_khi.loVect());
      IntVect kbx_hi(box_khi.hiVect());

      // Initialize z-hi domain extent.
      ParallelFor(box_khi, [bc_khi_type, init_z]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {bc_khi_type(i,j,k,0) = init_z;});

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_zhi().size(); ++lc) {

        const int bcv  = bc_zhi(lc);
        const int type = m_bc[bcv].type;

        if (lc > 0){
          kbx_lo[0] = m_bc[bcv].Lo(0,plo[0],dx);
          kbx_lo[1] = m_bc[bcv].Lo(1,plo[1],dy);

          kbx_hi[0] = m_bc[bcv].Hi(0,plo[0],dx);
          kbx_hi[1] = m_bc[bcv].Hi(1,plo[1],dy);
        }

        const Box hi_box(kbx_lo, kbx_hi);
        ParallelFor(hi_box, [bc_khi_type, type, bcv]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          bc_khi_type(i,j,k,0) = type;
          bc_khi_type(i,j,k,1) = bcv;
        });
      }
    } // end z-hi side of the domain
  }// end z-direction

}
