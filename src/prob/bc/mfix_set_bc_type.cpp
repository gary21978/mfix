#include <mfix_reporter.H>
#include <mfix_bc.H>
#include <mfix_fluid.H>

using namespace amrex;

void
MFIXBoundaryConditions::
set_bc_type ( int a_lev, int a_nghost,
              MFIXFluidPhase& a_fluid )
{
    Real dx = m_geom[a_lev].CellSize(0);
    Real dy = m_geom[a_lev].CellSize(1);
    Real dz = m_geom[a_lev].CellSize(2);

    const GpuArray<Real, 3> plo = m_geom[a_lev].ProbLoArray();

    const int l_species = a_fluid.nspecies();
    const int l_ntrac = a_fluid.ntracer();

    // Set the defaults for BCRecs
    m_bcrec_velocity.resize(AMREX_SPACEDIM);
    m_bcrec_hydro_velocity.resize(AMREX_SPACEDIM);
    m_bcrec_volfrac.resize(1);
    m_bcrec_density.resize(1);
    m_bcrec_enthalpy.resize(1);
    m_bcrec_tracer.resize(l_ntrac);
    m_bcrec_species.resize(l_species);

    { // begin x direction

      const int dir = 0;

      Array4<int> const& bc_ilo_type = m_bc_list.bc_ilo[a_lev]->array();
      Array4<int> const& bc_ihi_type = m_bc_list.bc_ihi[a_lev]->array();

      const int init_x = m_geom[a_lev].isPeriodic(0) ? BCList::undefined : BCList::cover;

      Box domainx(m_geom[a_lev].Domain());
      domainx.grow(1,a_nghost);  // Add ghost cells to y
      domainx.grow(2,a_nghost);  // Add ghost cells to z

      { // x-lo side of the domain

        Box box_ilo = amrex::adjCellLo(domainx,0,1);

        IntVect ibx_lo(box_ilo.loVect());
        IntVect ibx_hi(box_ilo.hiVect());

        int xlo_type = init_x;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_ilo, [bc_ilo_type, init_x]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_ilo_type(i,j,k,0) = init_x;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_xlo().size(); ++lc) {

          const int bcv  = bc_xlo(lc);
          const int type = bc(bcv).type;

          const BC_t& bct = bc(bcv);

          xlo_type = type;

          if (lc > 0){
            ibx_lo[1] = static_cast<int>(amrex::Math::floor((bct.region->lo(1)-plo[1])/dy + 0.5));
            ibx_lo[2] = static_cast<int>(amrex::Math::floor((bct.region->lo(2)-plo[2])/dz + 0.5));

            ibx_hi[1] = static_cast<int>(amrex::Math::floor((bct.region->hi(1)-plo[1])/dy + 0.5)-1);
            ibx_hi[2] = static_cast<int>(amrex::Math::floor((bct.region->hi(2)-plo[2])/dz + 0.5)-1);
          }

          const Box lo_box(ibx_lo, ibx_hi);

          amrex::ParallelFor(lo_box, [bc_ilo_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_ilo_type(i,j,k,0) = type;
               bc_ilo_type(i,j,k,1) = bcv;
             });
        }

        set_bcrec_lo(a_lev, dir, xlo_type);

      } // end x-lo side of the domain


      { // x-hi side of the domain

        Box box_ihi = amrex::adjCellHi(domainx,0,1);

        IntVect ibx_lo(box_ihi.loVect());
        IntVect ibx_hi(box_ihi.hiVect());

        int xhi_type = init_x;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_ihi, [bc_ihi_type, init_x]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_ihi_type(i,j,k,0) = init_x;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_xhi().size(); ++lc) {

          const int bcv  = bc_xhi(lc);
          const int type = bc(bcv).type;

          const BC_t& bct = bc(bcv);

          xhi_type = type;

          if (lc > 0){
            ibx_lo[1] = static_cast<int>(amrex::Math::floor((bct.region->lo(1)-plo[1])/dy + 0.5));
            ibx_lo[2] = static_cast<int>(amrex::Math::floor((bct.region->lo(2)-plo[2])/dz + 0.5));

            ibx_hi[1] = static_cast<int>(amrex::Math::floor((bct.region->hi(1)-plo[1])/dy + 0.5)-1);
            ibx_hi[2] = static_cast<int>(amrex::Math::floor((bct.region->hi(2)-plo[2])/dz + 0.5)-1);
          }

          const Box hi_box(ibx_lo, ibx_hi);

          amrex::ParallelFor(hi_box, [bc_ihi_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_ihi_type(i,j,k,0) = type;
               bc_ihi_type(i,j,k,1) = bcv;
             });

        }

        set_bcrec_hi(a_lev, dir, xhi_type);
      } // end x-hi side of the domain
    } // end x-direction


    { // begin y-direction

      const int dir = 1;

      const int init_y = m_geom[a_lev].isPeriodic(1) ? BCList::undefined : BCList::cover;

      Array4<int> const& bc_jlo_type = m_bc_list.bc_jlo[a_lev]->array();
      Array4<int> const& bc_jhi_type = m_bc_list.bc_jhi[a_lev]->array();

      Box domainy(m_geom[a_lev].Domain());
      domainy.grow(0,a_nghost);  // Add ghost cells to x
      domainy.grow(2,a_nghost);  // Add ghost cells to z

      { // y-lo side of the domain

        Box box_jlo = amrex::adjCellLo(domainy,1,1);

        IntVect jbx_lo(box_jlo.loVect());
        IntVect jbx_hi(box_jlo.hiVect());

        int ylo_type = init_y;

        // Initialize y-lo domain extent.
        amrex::ParallelFor(box_jlo, [bc_jlo_type, init_y]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_jlo_type(i,j,k,0) = init_y;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_ylo().size(); ++lc) {

          const int bcv  = bc_ylo(lc);
          const int type = bc(bcv).type;

          const BC_t& bct = bc(bcv);

          ylo_type = type;

          if (lc > 0){
            jbx_lo[0] = static_cast<int>(amrex::Math::floor((bct.region->lo(0)-plo[0])/dx + 0.5));
            jbx_lo[2] = static_cast<int>(amrex::Math::floor((bct.region->lo(2)-plo[2])/dz + 0.5));

            jbx_hi[0] = static_cast<int>(amrex::Math::floor((bct.region->hi(0)-plo[0])/dx + 0.5)-1);
            jbx_hi[2] = static_cast<int>(amrex::Math::floor((bct.region->hi(2)-plo[2])/dz + 0.5)-1);
          }

          const Box lo_box(jbx_lo, jbx_hi);

          amrex::ParallelFor(lo_box, [bc_jlo_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_jlo_type(i,j,k,0) = type;
               bc_jlo_type(i,j,k,1) = bcv;
             });

        }

        set_bcrec_lo(a_lev, dir, ylo_type);

      }// end y-lo side of the domain


      { // y-hi side of the domain

        Box box_jhi = amrex::adjCellHi(domainy,1,1);

        IntVect jbx_lo(box_jhi.loVect());
        IntVect jbx_hi(box_jhi.hiVect());

        int yhi_type = init_y;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_jhi, [bc_jhi_type, init_y]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_jhi_type(i,j,k,0) = init_y;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_yhi().size(); ++lc) {

          const int bcv  = bc_yhi(lc);
          const int type = bc(bcv).type;

          const BC_t& bct = bc(bcv);

          yhi_type = type;

          if (lc > 0){
            jbx_lo[0] = static_cast<int>(amrex::Math::floor((bct.region->lo(0)-plo[0])/dx + 0.5));
            jbx_lo[2] = static_cast<int>(amrex::Math::floor((bct.region->lo(2)-plo[2])/dz + 0.5));

            jbx_hi[0] = static_cast<int>(amrex::Math::floor((bct.region->hi(0)-plo[0])/dx + 0.5)-1);
            jbx_hi[2] = static_cast<int>(amrex::Math::floor((bct.region->hi(2)-plo[2])/dz + 0.5)-1);
          }

          const Box hi_box(jbx_lo, jbx_hi);

          amrex::ParallelFor(hi_box, [bc_jhi_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_jhi_type(i,j,k,0) = type;
               bc_jhi_type(i,j,k,1) = bcv;
             });
        }

        set_bcrec_hi(a_lev, dir, yhi_type);

      } // end y-hi side of the domain
    } // end y-direction

    { // begin z-direction

      const int dir = 2;

      const int init_z = m_geom[a_lev].isPeriodic(2) ? BCList::undefined : BCList::cover;

      Array4<int> const& bc_klo_type = m_bc_list.bc_klo[a_lev]->array();
      Array4<int> const& bc_khi_type = m_bc_list.bc_khi[a_lev]->array();

      Box domainz(m_geom[a_lev].Domain());
      domainz.grow(0,a_nghost);  // Add ghost cells to x
      domainz.grow(1,a_nghost);  // Add ghost cells to y

      { // z-lo side of the domain

        Box box_klo = amrex::adjCellLo(domainz,2,1);

        IntVect kbx_lo(box_klo.loVect());
        IntVect kbx_hi(box_klo.hiVect());

        int zlo_type = init_z;

        // Initialize y-lo domain extent.
        amrex::ParallelFor(box_klo, [bc_klo_type, init_z]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_klo_type(i,j,k,0) = init_z;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_zlo().size(); ++lc) {

          const int bcv  = bc_zlo(lc);
          const int type = bc(bcv).type;

          const BC_t& bct = bc(bcv);

          zlo_type = type;

          if (lc > 0){
            kbx_lo[0] = static_cast<int>(amrex::Math::floor((bct.region->lo(0)-plo[0])/dx + 0.5));
            kbx_lo[1] = static_cast<int>(amrex::Math::floor((bct.region->lo(1)-plo[1])/dy + 0.5));

            kbx_hi[0] = static_cast<int>(amrex::Math::floor((bct.region->hi(0)-plo[0])/dx + 0.5)-1);
            kbx_hi[1] = static_cast<int>(amrex::Math::floor((bct.region->hi(1)-plo[1])/dy + 0.5)-1);
          }

          const Box lo_box(kbx_lo, kbx_hi);

          amrex::ParallelFor(lo_box, [bc_klo_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_klo_type(i,j,k,0) = type;
               bc_klo_type(i,j,k,1) = bcv;
             });

        }

        set_bcrec_lo(a_lev, dir, zlo_type);

      } // end z-lo side of the domain


      { // z-hi side of the domain

        Box box_khi = amrex::adjCellHi(domainz,2,1);

        IntVect kbx_lo(box_khi.loVect());
        IntVect kbx_hi(box_khi.hiVect());

        int zhi_type = init_z;

        // Initialize x-lo domain extent.
        amrex::ParallelFor(box_khi, [bc_khi_type, init_z]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {bc_khi_type(i,j,k,0) = init_z;});

        // Define specific BC conditions from inputs
        for (int lc(0); lc < bc_zhi().size(); ++lc) {

          const int bcv  = bc_zhi(lc);
          const int type = bc(bcv).type;

          const BC_t& bct = bc(bcv);

          zhi_type = type;

          if (lc > 0){
            kbx_lo[0] = static_cast<int>(amrex::Math::floor((bct.region->lo(0)-plo[0])/dx + 0.5));
            kbx_lo[1] = static_cast<int>(amrex::Math::floor((bct.region->lo(1)-plo[1])/dy + 0.5));

            kbx_hi[0] = static_cast<int>(amrex::Math::floor((bct.region->hi(0)-plo[0])/dx + 0.5)-1);
            kbx_hi[1] = static_cast<int>(amrex::Math::floor((bct.region->hi(1)-plo[1])/dy + 0.5)-1);
          }

          const Box hi_box(kbx_lo, kbx_hi);

          amrex::ParallelFor(hi_box, [bc_khi_type, type, bcv]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
               bc_khi_type(i,j,k,0) = type;
               bc_khi_type(i,j,k,1) = bcv;
             });

        }

        set_bcrec_hi(a_lev, dir, zhi_type);

      } // end z-hi side of the domain
    }// end z-direction


    {
      m_bcrec_velocity_d.resize(AMREX_SPACEDIM);
      m_bcrec_hydro_velocity_d.resize(AMREX_SPACEDIM);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_velocity_d.data(), m_bcrec_velocity.data(), sizeof(BCRec)*AMREX_SPACEDIM);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_hydro_velocity_d.data(), m_bcrec_hydro_velocity.data(), sizeof(BCRec)*AMREX_SPACEDIM);
    }

    { // Always with volume fraction (epf)
      m_bcrec_volfrac_d.resize(1);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_volfrac_d.data(), m_bcrec_volfrac.data(), sizeof(BCRec));
    }

    { // Always include rho for fillpatch calls
      m_bcrec_density_d.resize(1);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_density_d.data(), m_bcrec_density.data(), sizeof(BCRec));
    }

    if (a_fluid.solve_enthalpy() || a_fluid.constraint.isIdealGas()) {
      m_bcrec_enthalpy_d.resize(1);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_enthalpy_d.data(), m_bcrec_enthalpy.data(), sizeof(BCRec));
    }

    if (l_ntrac > 0) {
      m_bcrec_tracer_d.resize(l_ntrac);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_tracer_d.data(), m_bcrec_tracer.data(), sizeof(BCRec)*l_ntrac);
    }

    if (l_species > 0) {
      m_bcrec_species_d.resize(l_species);
#ifdef AMREX_USE_GPU
      Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (m_bcrec_species_d.data(), m_bcrec_species.data(), sizeof(BCRec)*l_species);
    }


    if (a_fluid.solve()) {
      Real ltime(0.);

      set_velocity_bc_values(ltime);

      set_density_bc_values(ltime);

      if (a_fluid.solve_tracer()) {
        set_tracer_bc_values(ltime, a_fluid);
      }

      if (a_fluid.solve_species()) {
        set_species_bc_values(ltime, a_fluid.nspecies());
      }

      if (a_fluid.solve_enthalpy() || a_fluid.constraint.isIdealGas() ) {
        set_energy_bc_values(ltime, a_fluid);
      }

    }

    set_epf_values( a_fluid.solve() );
    set_pressure_values( a_fluid.solve() );

    Gpu::synchronize();
}


void MFIXBoundaryConditions::
set_bcrec_lo ( const int a_lev, const int dir, const int l_type)
{

  // Velocity BC Recs
  if (l_type == BCList::pinf) {

    m_bcrec_velocity[0].setLo(dir, BCType::foextrap);
    m_bcrec_velocity[1].setLo(dir, BCType::foextrap);
    m_bcrec_velocity[2].setLo(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[0].setLo(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[1].setLo(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[2].setLo(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[dir].setLo(dir, BCType::ext_dir);

  } else if (l_type == BCList::pout) {

    m_bcrec_velocity[0].setLo(dir, BCType::foextrap);
    m_bcrec_velocity[1].setLo(dir, BCType::foextrap);
    m_bcrec_velocity[2].setLo(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[0].setLo(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[1].setLo(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[2].setLo(dir, BCType::foextrap);

  } else if (l_type == BCList::minf ) {

    m_bcrec_velocity[0].setLo(dir, BCType::ext_dir);
    m_bcrec_velocity[1].setLo(dir, BCType::ext_dir);
    m_bcrec_velocity[2].setLo(dir, BCType::ext_dir);

    m_bcrec_hydro_velocity[0].setLo(dir, BCType::ext_dir);
    m_bcrec_hydro_velocity[1].setLo(dir, BCType::ext_dir);
    m_bcrec_hydro_velocity[2].setLo(dir, BCType::ext_dir);

  } else if (m_geom[a_lev].isPeriodic(dir)) {

    m_bcrec_velocity[0].setLo(dir, BCType::int_dir);
    m_bcrec_velocity[1].setLo(dir, BCType::int_dir);
    m_bcrec_velocity[2].setLo(dir, BCType::int_dir);

    m_bcrec_hydro_velocity[0].setLo(dir, BCType::int_dir);
    m_bcrec_hydro_velocity[1].setLo(dir, BCType::int_dir);
    m_bcrec_hydro_velocity[2].setLo(dir, BCType::int_dir);
  }

  // Scalar BC Recs
  if (l_type == BCList::pinf || l_type == BCList::pout ) {

    m_bcrec_volfrac[0].setLo(dir, BCType::ext_dir);
    m_bcrec_density[0].setLo(dir, BCType::foextrap);
    m_bcrec_enthalpy[0].setLo(dir, BCType::foextrap);
    for (auto& b : m_bcrec_tracer) b.setLo(dir, BCType::foextrap);
    for (auto& b : m_bcrec_species) b.setLo(dir, BCType::foextrap);

  } else if (l_type == BCList::minf) {

    m_bcrec_volfrac[0].setLo(dir, BCType::ext_dir);
    m_bcrec_density[0].setLo(dir, BCType::ext_dir);
    m_bcrec_enthalpy[0].setLo(dir, BCType::ext_dir);
    for (auto& b : m_bcrec_tracer) b.setLo(dir, BCType::ext_dir);
    for (auto& b : m_bcrec_species) b.setLo(dir, BCType::ext_dir);

  } else if (m_geom[a_lev].isPeriodic(dir)) {

    m_bcrec_volfrac[0].setLo(dir, BCType::int_dir);
    m_bcrec_density[0].setLo(dir, BCType::int_dir);
    m_bcrec_enthalpy[0].setLo(dir, BCType::int_dir);
    for (auto& b : m_bcrec_tracer) b.setLo(dir, BCType::int_dir);
    for (auto& b : m_bcrec_species) b.setLo(dir, BCType::int_dir);
  }

}



void MFIXBoundaryConditions::
set_bcrec_hi (const int a_lev, const int dir, const int l_type)
{

  // Velocity BC Recs
  if (l_type == BCList::pinf) {

    m_bcrec_velocity[0].setHi(dir, BCType::foextrap);
    m_bcrec_velocity[1].setHi(dir, BCType::foextrap);
    m_bcrec_velocity[2].setHi(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[0].setHi(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[1].setHi(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[2].setHi(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[dir].setHi(dir, BCType::ext_dir);

  } else if (l_type == BCList::pout) {

    m_bcrec_velocity[0].setHi(dir, BCType::foextrap);
    m_bcrec_velocity[1].setHi(dir, BCType::foextrap);
    m_bcrec_velocity[2].setHi(dir, BCType::foextrap);

    m_bcrec_hydro_velocity[0].setHi(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[1].setHi(dir, BCType::foextrap);
    m_bcrec_hydro_velocity[2].setHi(dir, BCType::foextrap);

  } else if (l_type == BCList::minf ) {

    m_bcrec_velocity[0].setHi(dir, BCType::ext_dir);
    m_bcrec_velocity[1].setHi(dir, BCType::ext_dir);
    m_bcrec_velocity[2].setHi(dir, BCType::ext_dir);

    m_bcrec_hydro_velocity[0].setHi(dir, BCType::ext_dir);
    m_bcrec_hydro_velocity[1].setHi(dir, BCType::ext_dir);
    m_bcrec_hydro_velocity[2].setHi(dir, BCType::ext_dir);

  } else if (m_geom[a_lev].isPeriodic(dir)) {

    m_bcrec_velocity[0].setHi(dir, BCType::int_dir);
    m_bcrec_velocity[1].setHi(dir, BCType::int_dir);
    m_bcrec_velocity[2].setHi(dir, BCType::int_dir);

    m_bcrec_hydro_velocity[0].setHi(dir, BCType::int_dir);
    m_bcrec_hydro_velocity[1].setHi(dir, BCType::int_dir);
    m_bcrec_hydro_velocity[2].setHi(dir, BCType::int_dir);
  }

  // Scalar BC Recs
  if (l_type == BCList::pinf || l_type == BCList::pout ) {

    m_bcrec_volfrac[0].setHi(dir, BCType::ext_dir);
    m_bcrec_density[0].setHi(dir, BCType::foextrap);
    m_bcrec_enthalpy[0].setHi(dir, BCType::foextrap);
    for (auto& b : m_bcrec_tracer) b.setHi(dir, BCType::foextrap);
    for (auto& b : m_bcrec_species) b.setHi(dir, BCType::foextrap);

  } else if (l_type == BCList::minf) {

    m_bcrec_volfrac[0].setHi(dir, BCType::ext_dir);
    m_bcrec_density[0].setHi(dir, BCType::ext_dir);
    m_bcrec_enthalpy[0].setHi(dir, BCType::ext_dir);
    for (auto& b : m_bcrec_tracer) b.setHi(dir, BCType::ext_dir);
    for (auto& b : m_bcrec_species) b.setHi(dir, BCType::ext_dir);

  } else if (m_geom[a_lev].isPeriodic(dir)) {

    m_bcrec_volfrac[0].setHi(dir, BCType::int_dir);
    m_bcrec_density[0].setHi(dir, BCType::int_dir);
    m_bcrec_enthalpy[0].setHi(dir, BCType::int_dir);
    for (auto& b : m_bcrec_tracer) b.setHi(dir, BCType::int_dir);
    for (auto& b : m_bcrec_species) b.setHi(dir, BCType::int_dir);
  }

}
