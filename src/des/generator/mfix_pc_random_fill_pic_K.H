#ifndef _MFIX_PC_GEN_RANDOM_FILL_PIC_K_H_
#define _MFIX_PC_GEN_RANDOM_FILL_PIC_K_H_

#include <limits>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>

#include <mfix_random.H>
#include <mfix_distributions.H>

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: random fill                                             !
//                                                                      !
//  Purpose: Generate initial solids packing based on randomly placing  !
//           particles in the ic region.                                !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
class RandomFill_PIC
{
  public:
    // Constructor
    RandomFill_PIC (const amrex::RealVect& plo,
                    const amrex::RealVect& dx,
                    const INPUT_DIST_t& a_diameter,
                    bool const a_generate_rand_seed = false)
      : m_plo(plo)
      , m_dx(dx)
      , m_diameter(a_diameter)
      , m_generate_rand_seed(a_generate_rand_seed)
      , m_np(0)
      , m_ordered(0)
      , m_max_i(0.)
      , m_max_ij(0.)
      , m_inv_wppc(std::numeric_limits<amrex::Real>::infinity())
      , m_inv_max_i(std::numeric_limits<amrex::Real>::infinity())
      , m_inv_max_ij(std::numeric_limits<amrex::Real>::infinity())
    {}

  private:

    const amrex::RealVect m_plo;

    const amrex::RealVect m_dx;

    INPUT_DIST_t m_diameter;

    int const m_generate_rand_seed;

    int m_np;

    int m_ordered;

    amrex::Real m_max_i;
    amrex::Real m_max_ij;

    // whole parcels per cell
    amrex::Real m_inv_wppc;

    amrex::Real m_inv_max_i;
    amrex::Real m_inv_max_ij;

    // sub-region of initial condition region
    amrex::RealBox m_rbx;

  public:

    void setup ( const amrex::Box& a_bx,
                 const SOLIDS_t& a_solid,
                 amrex::Real const a_multi_particle_volume )
    {

      int const rand_seed = m_generate_rand_seed ? aux::calc_random_int() :
        amrex::ParallelDescriptor::MyProc()+1;

      amrex::ResetRandomSeed(rand_seed);

      m_rbx = amrex::RealBox(a_bx, m_dx.begin(), m_plo.begin());

      // Number of particles to generate in this box within this subregion of the IC.
      m_np = a_solid.calc_particles_from_vol(m_rbx.volume(), a_multi_particle_volume);

      const int cells = static_cast<int>(a_bx.numPts());

      // number of particles that will be placed in a somewhat ordered
      // fashion to get good coverage of the initial condition region.
      m_ordered = m_np - (m_np % cells);

      // whole parcels per cell
      int const wppc = m_ordered / cells;

      if ( wppc > 0 ) {

        m_inv_wppc = 1.0/static_cast<amrex::Real>(wppc);

        m_max_i  = static_cast<amrex::Real>(wppc*a_bx.length(0));
        m_max_ij = m_max_i*static_cast<amrex::Real>(a_bx.length(1));

        m_inv_max_i  = 1.0/m_max_i;
        m_inv_max_ij = 1.0/m_max_ij;
      }
    }

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect get_position (const int a_pid,
                                amrex::RandomEngine const& engine) const
    {

      amrex::RealVect position(0.);

      if ( a_pid < m_ordered) {

        amrex::Real const pid_r(a_pid);

        const int k = static_cast<int>(pid_r * m_inv_max_ij);
        amrex::Real const k_r(static_cast<amrex::Real>(k));

        const int j = static_cast<int>( (pid_r - k_r*m_max_ij) * m_inv_max_i );
        amrex::Real const j_r(static_cast<amrex::Real>(j));

        const int i = static_cast<int>( (pid_r - k_r*m_max_ij - j_r*m_max_i) * m_inv_wppc );
        amrex::Real const i_r(static_cast<amrex::Real>(i));

        position = { m_rbx.lo(0) + (i_r + amrex::Random(engine)) * m_dx[0],
                     m_rbx.lo(1) + (j_r + amrex::Random(engine)) * m_dx[1],
                     m_rbx.lo(2) + (k_r + amrex::Random(engine)) * m_dx[2]};

      } else {

        position = { m_rbx.lo(0) + m_rbx.length(0)*amrex::Random(engine),
                     m_rbx.lo(1) + m_rbx.length(1)*amrex::Random(engine),
                     m_rbx.lo(2) + m_rbx.length(2)*amrex::Random(engine)};

      }

      return position;
    }

};

#endif
