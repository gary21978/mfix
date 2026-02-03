#ifndef _MFIX_PC_GEN_HEX_CLOSE_PACK_K_H_
#define _MFIX_PC_GEN_HEX_CLOSE_PACK_K_H_

#include <limits>

#include <mfix_ic.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_PODVector.H>

// This is a very small number we add to the length to catch
// any particles that are "on the edge"
#define small_number 10.*std::numeric_limits<amrex::Real>::epsilon()

#define loc_tolerance std::numeric_limits<amrex::Real>::epsilon()

#define sqrt6x2o3 std::sqrt(6.0)*(2.0/3.0)
#define sqrt3o3 std::sqrt(3.0)/3.0
#define sqrt3 std::sqrt(3.0)

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: hex_close_pack                                          !
//                                                                      !
//  Purpose: Generate initial solids packing based on hexagonal close   !
//           packing of mono-sized spheres.                             !
//                                                                      !
//  TODO: * generalize fill direction to follow gravity.                !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
class Hex_ClosePack
{

  public:

    // Constructor
    AMREX_GPU_HOST_DEVICE
    Hex_ClosePack ()
      : m_is_initialized(0)
      , m_plo(0.)
      , m_dx(0.)
    { initialize_member_variables(); }

    AMREX_GPU_HOST_DEVICE
    Hex_ClosePack ( const amrex::RealVect& plo,
                    const amrex::RealVect& dx )
      : m_is_initialized(1)
      , m_plo(plo)
      , m_dx(dx)
    { initialize_member_variables(); }

  private:

    int m_is_initialized;

    amrex::RealVect m_plo;

    amrex::RealVect m_dx;

    int m_np;

    // Number of particles in plane type 1; row type 1
    int m_p1r1_np;

    // Number of particles in plane type 2; row type 1
    int m_p2r1_np;

    // Number of particles in plane type 1
    int m_p1_np;

    // Number of particles in combined type 1 + type 2 plane
    int m_p1p2_np;

    // Number of particles in combined type 1 + type 2 plane
    int m_r1r2_np;

    int m_klo;

    int m_p1_jlo;

    int m_p1r1_ilo;
    int m_p1r2_ilo;

    int m_p2_jlo;
    int m_p2r1_ilo;
    int m_p2r2_ilo;

    amrex::RealVect m_ic_dlo;

    amrex::Real m_eff_rad;

    amrex::RealVect m_rbx_lo;
    amrex::RealVect m_rbx_hi;

  public:

    AMREX_GPU_HOST_DEVICE
    void initialize ( const amrex::RealVect& plo,
                      const amrex::RealVect& dx )
    { m_is_initialized = 1;

      for (int dir(0); dir < 3; ++dir) {
        m_plo[dir] = plo[dir];
        m_dx[dir] = dx[dir];
      }
    }

    AMREX_GPU_HOST_DEVICE
    void setup (amrex::Box const& a_bx,
                const amrex::RealBox& region,
                const amrex::Real mean,
                const amrex::Real eps,
                const amrex::Real eps_tolerance = 1.e-15)
    {

     const amrex::IntVect lo(a_bx.loVect());
     const amrex::IntVect hi(a_bx.hiVect());

      // physical volume of IC region
      const amrex::Real vol = region.volume();

      m_np = 0;

      if (!m_is_initialized) { return; }
      if (!(eps > eps_tolerance && eps < 1.)) { return; }

      // Total particle count is based on mean particle size
      const long seed =
        static_cast<long>(vol * eps / ((M_PI/6.0)*mean*mean*mean));

      if (seed == 0) { return; }

      m_ic_dlo = amrex::RealVect(region.lo(0), region.lo(1), region.lo(2));

      amrex::RealVect ic_dhi(region.hi(0), region.hi(1), region.hi(2));

      const amrex::Real xlen = ic_dhi[0] - m_ic_dlo[0];
      const amrex::Real ylen = ic_dhi[1] - m_ic_dlo[1];
      const amrex::Real zlen = ic_dhi[2] - m_ic_dlo[2];

      amrex::RealVect rbx_lo(m_plo[0] + m_dx[0]*static_cast<amrex::Real>(lo[0]),
                             m_plo[1] + m_dx[1]*static_cast<amrex::Real>(lo[1]),
                             m_plo[2] + m_dx[2]*static_cast<amrex::Real>(lo[2]));

      amrex::RealVect rbx_hi(m_plo[0] + m_dx[0]*static_cast<amrex::Real>(hi[0]+1),
                             m_plo[1] + m_dx[1]*static_cast<amrex::Real>(hi[1]+1),
                             m_plo[2] + m_dx[2]*static_cast<amrex::Real>(hi[2]+1));

      const amrex::Real mean_rad = 0.5*mean;

      // Verify that the IC region overlaps with the IC region.
      if (ic_dhi[0] <= rbx_lo[0] + mean_rad || rbx_hi[0] <= m_ic_dlo[0] + mean_rad ||
          ic_dhi[1] <= rbx_lo[1] + mean_rad || rbx_hi[1] <= m_ic_dlo[1] + mean_rad ||
          ic_dhi[2] <= rbx_lo[2] + mean_rad || rbx_hi[2] <= m_ic_dlo[2] + mean_rad) {

        return;
      }

      // This is the radius we need for spacing particles to get the
      // desired volume fraction and particle count.
      m_eff_rad = hcp_rad(mean_rad, xlen, ylen, zlen, seed);
      //printf ("\nPacking radius: %21.18f\n", m_eff_rad);

      // Pull back just less than one particle radius to keep particles
      // fully inside the IC region.
      ic_dhi -= (m_eff_rad - small_number);

      rbx_hi[0] = amrex::min(rbx_hi[0], ic_dhi[0]);
      rbx_hi[1] = amrex::min(rbx_hi[1], ic_dhi[1]);
      rbx_hi[2] = amrex::min(rbx_hi[2], ic_dhi[2]);

      m_rbx_hi = rbx_hi;
      m_rbx_lo = rbx_lo;

      // Transition into relative coordinates with origin at (0,0,0)
      rbx_lo -= m_ic_dlo;
      rbx_hi -= m_ic_dlo;

      // Compute the first plane.

      m_klo = hcp_k(m_eff_rad, rbx_lo[2]);        // k index of first plane
      const int p_cnt = hcp_k(m_eff_rad, rbx_hi[2], m_klo) - m_klo; // Total number of planes


      const int p2_cnt = (p_cnt - (p_cnt%2))/2;

      const int p1_cnt = p_cnt - p2_cnt;

      m_p1_jlo = hcp_j(m_eff_rad, rbx_lo[1], (m_klo%2));
      const int p1_rows  = hcp_j(m_eff_rad, rbx_hi[1], (m_klo%2), m_p1_jlo) - m_p1_jlo;


      const int p1r2_cnt = (p1_rows - (p1_rows%2))/2;
      const int p1r1_cnt = p1_rows - p1r2_cnt;

      m_p1r1_ilo = hcp_i(m_eff_rad, rbx_lo[0], ((m_klo+m_p1_jlo  )%2));
      m_p1r2_ilo = hcp_i(m_eff_rad, rbx_lo[0], ((m_klo+m_p1_jlo+1)%2));

      m_p1r1_np  = hcp_i(m_eff_rad, rbx_hi[0], ((m_klo+m_p1_jlo  )%2), m_p1r1_ilo) - m_p1r1_ilo;
      const int p1r2_np  = hcp_i(m_eff_rad, rbx_hi[0], ((m_klo+m_p1_jlo+1)%2), m_p1r2_ilo) - m_p1r2_ilo;

      m_p1_np = m_p1r1_np*p1r1_cnt + p1r2_np*p1r2_cnt;

      // Build out plane 2
      m_p2_jlo = hcp_j(m_eff_rad, rbx_lo[1], ((m_klo+1)%2));
      const int p2_rows  = hcp_j(m_eff_rad, rbx_hi[1], ((m_klo+1)%2), m_p2_jlo) - m_p2_jlo;


      const int p2r2_cnt = (p2_rows - (p2_rows%2))/2;
      const int p2r1_cnt = p2_rows - p2r2_cnt;

      m_p2r1_ilo = hcp_i(m_eff_rad, rbx_lo[0], ((m_klo+1+m_p2_jlo  )%2));
      m_p2r2_ilo = hcp_i(m_eff_rad, rbx_lo[0], ((m_klo+1+m_p2_jlo+1)%2));

      m_p2r1_np  = hcp_i(m_eff_rad, rbx_hi[0], ((m_klo+1+m_p2_jlo  )%2), m_p2r1_ilo) - m_p2r1_ilo;
      const int p2r2_np  = hcp_i(m_eff_rad, rbx_hi[0], ((m_klo+1+m_p2_jlo+1)%2), m_p2r2_ilo) - m_p2r2_ilo;

      const int p2_np = m_p2r1_np*p2r1_cnt + p2r2_np*p2r2_cnt;

      m_np = m_p1_np*p1_cnt + p2_np*p2_cnt;

      m_p1p2_np = m_p1_np + p2_np;

      m_r1r2_np = m_p1r1_np + p1r2_np;

    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int get_particles_number () const
    { return m_np; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_particles_radius () const
    { return m_eff_rad; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect get_position (const int p_id,
                                  amrex::RandomEngine const&) const
    {

      // Compute the grouped (type-1 + type-2) plane containing the particle.
      const int lk0 = p_id/m_p1p2_np; // (int/int) no need for floor()
      // Compute the particle index in the grouped plane.
      const int lk1 = p_id - lk0*m_p1p2_np;
      // 0: particle is in type-1 plane
      // 1: particle is in type-2 plane
      const int lk2 = (lk1 < m_p1_np ? 0 : 1);

      // Particle index in the hcp on this grid.
      const int local_k = lk2 + 2*lk0;
      // Particle index in the full hcp.
      const int k = m_klo + local_k;

      // Particle index in the k-th plane.
      const int local_ko2 = local_k/2; // (int/int) no need for floor()
      const int lj0 = p_id - ((local_k%2)*m_p1_np + local_ko2*m_p1p2_np);
      // Index of grouped (type-1 + type-2) row in the k-th plane.
      const int lj1 = lj0/m_r1r2_np; // (int/int) not need for floor()
      // Particle index within the grouped row
      const int lj2 = lj0 - lj1*m_r1r2_np;
      // Number of particles in first row or this plane
      const int lj3 = (1-(local_k%2))*m_p1r1_np + (local_k%2)*m_p2r1_np;

      // Particle index in the hcp on this grid.
      // 0: particle is in the first row of the group
      // 1: particle is in the second row of the group
      const int local_j = (lj2 < lj3 ? 0 : 1) + 2*lj1;
      // Particle index in the full hcp.
      const int j = local_j + (1-local_k%2)*m_p1_jlo
                            + (  local_k%2)*m_p2_jlo;


      // Particle index in the hcp on this grid.
      const int li0 = (lj2 < lj3 ? 0 : 1)*lj3;
      const int local_i = lj2 - li0;
      // Particle index in the full hcp.
      const int i = local_i + (1-(local_k%2))*((1-(local_j%2))*m_p1r1_ilo  +
                                               (  (local_j%2))*m_p1r2_ilo) +
                              (  (local_k%2))*((1-(local_j%2))*m_p2r1_ilo  +
                                               (  (local_j%2))*m_p2r2_ilo);

      const amrex::Real x = m_ic_dlo[0] + m_eff_rad * (1.0 + static_cast<amrex::Real>(2*i + (j+k)%2));
      const amrex::Real y = m_ic_dlo[1] + m_eff_rad * (1.0 + sqrt3o3*(static_cast<amrex::Real>(3*j + k%2)));
      const amrex::Real z = m_ic_dlo[2] + m_eff_rad * (1.0 + sqrt6x2o3*(static_cast<amrex::Real>(k)));

      AMREX_ALWAYS_ASSERT(m_rbx_lo[0] <= x && x < m_rbx_hi[0]);
      AMREX_ALWAYS_ASSERT(m_rbx_lo[1] <= y && y < m_rbx_hi[1]);
      AMREX_ALWAYS_ASSERT(m_rbx_lo[2] <= z && z < m_rbx_hi[2]);

      return amrex::RealVect(x, y, z);
    }


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void initialize_member_variables ()
    {
      m_np = 0;

      m_p1r1_np  = 0;
      m_p2r1_np  = 0;
      m_p1_np    = 0;
      m_p1p2_np  = 0;
      m_r1r2_np  = 0;
      m_klo      = 0;
      m_p1_jlo   = 0;
      m_p1r1_ilo = 0;
      m_p1r2_ilo = 0;
      m_p2_jlo   = 0;
      m_p2r1_ilo = 0;
      m_p2r2_ilo = 0;

      m_eff_rad = 0.;

      m_ic_dlo = amrex::RealVect::TheZeroVector();
      m_rbx_lo = amrex::RealVect::TheZeroVector();
      m_rbx_hi = amrex::RealVect::TheZeroVector();
    }


  private:

    // Auxiliary functions
    AMREX_GPU_HOST_DEVICE
    int hcp_i (amrex::Real rad,
               amrex::Real xend,
               int offset_int,
               int start_i = 0) const
    {
      const amrex::Real offset = static_cast<amrex::Real>(offset_int);
      int i(start_i);
      while ( loc_tolerance < xend - (1.0 + 2.0*static_cast<amrex::Real>(i) + offset)*rad ) {
        i++;
      }
      return i;
    }

    AMREX_GPU_HOST_DEVICE
    int hcp_j (amrex::Real rad,
               amrex::Real yend,
               int offset_int,
               int start_j = 0) const
    {
      const amrex::Real offset = static_cast<amrex::Real>(offset_int);
      int j(start_j);
      while ( loc_tolerance < yend - (1.0 + sqrt3o3 * (3.0 * static_cast<amrex::Real>(j) + offset))*rad ) {
        j++;
      }
      return j;
    }


    AMREX_GPU_HOST_DEVICE
    int hcp_k (amrex::Real rad,
               amrex::Real zend,
               int start_k = 0) const
    {
      int k(start_k);
      while ( loc_tolerance < zend - (1.0 + sqrt6x2o3*static_cast<amrex::Real>(k))*rad ) {
        k++;
      }
      return k;
    }

    AMREX_GPU_HOST_DEVICE
    long hcp_count (amrex::Real rad,
                    amrex::Real xlen,
                    amrex::Real ylen,
                    amrex::Real zlen) const
    {


      // Compute the first plane.

      const long p_cnt = static_cast<long>(hcp_k(rad, zlen)); // Total number of planes

      const long p2_cnt = (p_cnt - (p_cnt%2))/2;
      const long p1_cnt = p_cnt - p2_cnt;

      const long p1_rows = static_cast<long>(hcp_j(rad, ylen, 0));

      const long p1r2_cnt = (p1_rows - (p1_rows%2))/2;
      const long p1r1_cnt = p1_rows - p1r2_cnt;

      const long p1r1_np = static_cast<long>(hcp_i(rad, xlen, 0));
      const long p1r2_np = static_cast<long>(hcp_i(rad, xlen, 1));

      const long p1_np = p1r1_np*p1r1_cnt + p1r2_np*p1r2_cnt;

      // Build out plane 2
      const long p2_rows = static_cast<long>(hcp_j(rad, ylen, 1));

      const long p2r2_cnt = (p2_rows - (p2_rows%2))/2;
      const long p2r1_cnt = p2_rows - p2r2_cnt;

      const long p2r1_np  = static_cast<long>(hcp_i(rad, xlen, 1));
      const long p2r2_np  = static_cast<long>(hcp_i(rad, xlen, 0));

      const int p2_np = p2r1_np*p2r1_cnt + p2r2_np*p2r2_cnt;

      // Total number of particles
      const long np = p1_np*p1_cnt + p2_np*p2_cnt;

      return np;

    }

    AMREX_GPU_HOST_DEVICE
    amrex::Real hcp_rad (amrex::Real rad_in,
                         amrex::Real xlen,
                         amrex::Real ylen,
                         amrex::Real zlen,
                         long target) const
    {
      amrex::Real high = 1.0;
      amrex::Real low  = 1.0;

      amrex::Real radius(rad_in);

      const int max_iter = 15;
      int iter(0);

      long pc0 = hcp_count(rad_in, xlen, ylen, zlen);
      while ((pc0 > target) && (iter < max_iter)) {
        low = high;
        high *= 2.;
        radius = high*rad_in;
        pc0 = hcp_count(radius, xlen, ylen, zlen);

        iter++;
      }

      amrex::Real scale=0.5*(low + high);

      iter = 0;
      int keep_searching(1);
      while (keep_searching) {

        radius = scale*rad_in;
        long pc = hcp_count(radius, xlen, ylen, zlen);

        if(pc == target || (pc == pc0 && iter > 0)) {
          scale = 1.0;

        } else if (pc < target) {
          high  = scale;
          scale = 0.5*(low + scale);

        } else {
          high *= ( (scale >= high) ? 2. : 1.);
          low   = amrex::max(scale, 1.0);
          scale = 0.5*(high + scale);

        }

        pc0 = pc;
        iter++;
        keep_searching = ((scale != 1.0) && (iter < max_iter)) ? 1 : 0;
      }

      return radius;
    }
};

#endif
