#ifndef _MFIX_PC_GEN_RANDOM_FILL_DEM_K_H_
#define _MFIX_PC_GEN_RANDOM_FILL_DEM_K_H_

#define DEBUG_RANDOM_FILL_DEM 0

#include <limits>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_PODVector.H>
#include <AMReX_GpuAsyncArray.H>

#include <mfix_random.H>
#include <mfix_distributions.H>
#include <mfix_run_on.H>

class RandomFill_DEM
{

  public:

    // Constructor
    RandomFill_DEM (const amrex::RealVect& a_plo,
                    const amrex::RealVect& a_dx,
                    const INPUT_DIST_t& a_diameter,
                    bool const a_generate_rand_seed = false)
      : m_max_fails(1000)
      , m_max_seq_fails(100)
      , m_plo(a_plo)
      , m_dx(a_dx)
      , m_diameter(a_diameter)
      , m_is_const_rad(a_diameter.is_constant())
      , m_mean_rad(0.5*a_diameter.get_mean())
      , m_generate_rand_seed(a_generate_rand_seed)
      , m_np(0)
      , m_offset(0.)
    {}

  private:

    int const m_max_fails;

    int const m_max_seq_fails;

    amrex::RealVect const m_plo;

    amrex::RealVect const m_dx;

    INPUT_DIST_t m_diameter;

    int const m_is_const_rad;

    amrex::Real m_mean_rad;

    int const m_generate_rand_seed;

    int m_np;

    amrex::RealVect m_offset;

    amrex::Real* m_h_pos_ptr;
    amrex::Real* m_d_pos_ptr;

  public:

    void setup ( const amrex::Box& a_bx,
                 const SOLIDS_t& a_solid,
                 amrex::Real const a_multi_particle_volume,
                 amrex::Gpu::HostVector<amrex::Real>& h_pos,
                 amrex::Gpu::HostVector<amrex::Real>& h_rad,
                 amrex::Gpu::DeviceVector<amrex::Real>& d_pos,
                 amrex::Gpu::DeviceVector<amrex::Real>& d_rad)
    {
      h_pos.clear(); h_rad.clear();
      d_pos.clear(); d_rad.clear();

      amrex::RealBox const r_bx(a_bx, m_dx.begin(), m_plo.begin());

      m_offset = {r_bx.lo(0), r_bx.lo(1), r_bx.lo(2)};

      amrex::Real const lbin = m_diameter.get_max()*(m_is_const_rad ? 2. : 1.);

      int const max_i( static_cast<int>(r_bx.length(0)/lbin) );
      int const max_j( static_cast<int>(r_bx.length(1)/lbin) );
      int const max_k( static_cast<int>(r_bx.length(2)/lbin) );

      int const max_ij(max_i*max_j);
      int const max_ijk(max_ij*max_k);

#if DEBUG_RANDOM_FILL_DEM
      amrex::Print(amrex::Print::AllProcs)
                     << "Cells " << amrex::IntVect(max_i, max_j, max_k) << '\n'
                     << "   max_ij " <<  max_ij << '\n'
                     << "  max_ijk " << max_ijk << '\n';
#endif

      amrex::RealVect const idx =
        { static_cast<amrex::Real>(max_i)/r_bx.length(0),
          static_cast<amrex::Real>(max_j)/r_bx.length(1),
          static_cast<amrex::Real>(max_k)/r_bx.length(2)};

      // Number of particles to generate in this box within this
      // subregion of the IC.
      int const max_np(a_solid.calc_particles_from_vol(
                  r_bx.volume(), a_multi_particle_volume));

#if DEBUG_RANDOM_FILL_DEM
      amrex::Print(amrex::Print::AllProcs) << "Number to generate: " << max_np << '\n';
#endif

      // Running tally of number of particles in this cell.
      amrex::Vector<int> pinc(max_ijk, 0);

      // local ids of particles in this cell
      int nb(amrex::max(8, 2 + max_np/max_ijk));
      amrex::Vector<int> pbin(max_ijk*nb, -1);

      h_pos.reserve(max_np*3);
      if ( !m_is_const_rad ) { h_rad.reserve(max_np); }

      amrex::Real const const_rad_min_dist = 4.0*m_mean_rad*m_mean_rad;

      int const rand_seed = m_generate_rand_seed ? aux::calc_random_int() :
        amrex::ParallelDescriptor::MyProc()+1;

#ifdef AMREX_USE_GPU
      amrex::RandomEngine engine(nullptr);
#else
      amrex::RandomEngine engine;
#endif
      amrex::ResetRandomSeed(rand_seed);

      m_np = 0;

      // Number of sequential fails.
      int seq_fails(0);

      while(m_np < max_np && seq_fails < m_max_seq_fails) {

        int overlaps, fails(0);

        do {

          overlaps = 0;

          amrex::Real const rad_p = m_is_const_rad ? m_mean_rad :
            0.5*m_diameter.sample<amrex::RunOn::Host>(engine);

          // We don't want the particle to extend past the edge
          // of the current box because this could result in overlaps
          // with particles generated on other grids.
          amrex::RealVect const pos =
            { rad_p + (r_bx.length(0)-2.0*rad_p)*amrex::Random(engine),
              rad_p + (r_bx.length(1)-2.0*rad_p)*amrex::Random(engine),
              rad_p + (r_bx.length(2)-2.0*rad_p)*amrex::Random(engine) };

          // cell containing the new particle
          amrex::IntVect const cell =
            { static_cast<int>(amrex::Math::floor(pos[0]*idx[0])),
              static_cast<int>(amrex::Math::floor(pos[1]*idx[1])),
              static_cast<int>(amrex::Math::floor(pos[2]*idx[2])) };

          // Local grid search for collisions.
          for(int k = amrex::max(cell[2]-1, 0); k<=amrex::min(cell[2]+1, max_k-1); ++k) {
            for(int j = amrex::max(cell[1]-1, 0); j<=amrex::min(cell[1]+1, max_j-1); ++j) {
              for(int i = amrex::max(cell[0]-1, 0); i<=amrex::min(cell[0]+1, max_i-1); ++i) {

                int const ijk = i + j*(max_i) + k*(max_ij);

                for (int p(0); p<pinc[ijk]; ++p) {

                  // neighbor particle index
                  int const nbr = pbin[nb*ijk + p];

                  // particle data locally stored as:
                  // p0(x,y,z), p1(x,y,z), p2(x,y,z), ...
                  const amrex::Real dx = h_pos[nbr*3 + 0] - pos[0];
                  const amrex::Real dy = h_pos[nbr*3 + 1] - pos[1];
                  const amrex::Real dz = h_pos[nbr*3 + 2] - pos[2];

                  const amrex::Real dist = dx*dx + dy*dy + dz*dz;

                  if (m_is_const_rad) {
                    if (dist <= const_rad_min_dist)  { overlaps++; }
                  } else {
                    const amrex::Real sum_rp = rad_p + h_rad[nbr];
                    if(dist <= (sum_rp*sum_rp)) { overlaps++; }
                  }
                } // p
              } // i
            } // j
          } // k

          if (overlaps == 0) {

            h_pos.push_back(pos[0]);
            h_pos.push_back(pos[1]);
            h_pos.push_back(pos[2]);

            if (!m_is_const_rad) { h_rad.push_back(rad_p); }

            int const ijk = cell[0] + cell[1]*(max_i) + cell[2]*(max_ij);

            if(pinc[ijk] == nb - 1) {

              int const old_nb(nb);
              nb += 4;

              amrex::Vector<int> old_pbin(max_ijk*nb, -1);
              std::swap(pbin, old_pbin);

              for(int k = 0; k < max_k; ++k) {
                for(int j = 0; j < max_j; ++j) {
                  for(int i = 0; i < max_i; ++i) {

                    int const ijk_1d = i + j*(max_i) + k*(max_ij);

                    for (int p(0); p<pinc[ijk_1d]; ++p) {

                      pbin[nb*ijk_1d + p] = old_pbin[old_nb*ijk_1d + p];

                    }
                  }
                }
              }
            }

            pbin[nb*ijk + pinc[ijk]] = m_np;

            pinc[ijk]++;

            m_np++;

            seq_fails = 0;

          } else {

            fails++;

          }

        } while(overlaps > 0 && fails < m_max_fails);

        if (overlaps > 0) { seq_fails++; }
#if DEBUG_RANDOM_FILL_DEM
        if ( seq_fails != 0 && (seq_fails%50 == 0)) {
          amrex::Print(amrex::Print::AllProcs)
            << "seq fails " << seq_fails << '\n';
        }
#endif
      }


      d_pos.resize(3*m_np);

      amrex::Gpu::copy(amrex::Gpu::hostToDevice,
        h_pos.begin(), h_pos.end(), d_pos.begin());

      m_h_pos_ptr = h_pos.data();
      m_d_pos_ptr = d_pos.data();

      if (!m_is_const_rad) {
        m_diameter.set_sample_pointers( m_np, h_rad, d_rad);
      }

    }

    int get_particles_number () const
    { return m_np; }


    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect get_position (const int p_id,
                                  amrex::RandomEngine const&) const
    {
      amrex::Real x;
      amrex::Real y;
      amrex::Real z;

      if (run_on == amrex::RunOn::Gpu) {
        x = m_d_pos_ptr[p_id*3+0];
        y = m_d_pos_ptr[p_id*3+1];
        z = m_d_pos_ptr[p_id*3+2];
      } else if (run_on == amrex::RunOn::Cpu) {
        x = m_h_pos_ptr[p_id*3+0];
        y = m_h_pos_ptr[p_id*3+1];
        z = m_h_pos_ptr[p_id*3+2];
      }

      return amrex::RealVect(x, y, z) + m_offset;
    }

};

#endif
