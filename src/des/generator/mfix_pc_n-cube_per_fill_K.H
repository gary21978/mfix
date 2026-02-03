#ifndef _MFIX_PC_GEN_N_CUBE_PER_K_H_
#define _MFIX_PC_GEN_N_CUBE_PER_K_H_

#include <limits>

#include <mfix_ic.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_PODVector.H>

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: n_cube_per_fill                                         !
//                                                                      !
//  Purpose: Generate initial solids packing based on putting one       !
//           per cell                                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
class nCubePer_Fill
{
  public:
    // Constructor
    nCubePer_Fill ( const amrex::RealVect& a_plo,
                    const amrex::RealVect& a_dx,
                    const INPUT_DIST_t& a_diameter,
                    const int a_base)
      : m_plo(a_plo)
      , m_dx(a_dx)
      , m_diameter(a_diameter)
      , m_base(a_base)
      , m_np(0)
    {
      // This assertion should not be needed.
      AMREX_ASSERT(amrex::almostEqual(m_dx[0],m_dx[1]) &&
                   amrex::almostEqual(m_dx[0],m_dx[2]));
    }

  private:

    const amrex::RealVect m_plo;

    const amrex::RealVect m_dx;

    INPUT_DIST_t m_diameter;

    const int m_base;

    int m_np;

    amrex::IntVect m_delta_bx;

    amrex::IntVect m_seed_lo;

    amrex::Real m_dxn;

  public:

    void setup ( amrex::Box const& a_bx)
    {
      amrex::Real const max_dp(m_diameter.get_max());

      // Scaled m_dx[0] for particle spacing
      m_dxn = m_dx[0]/static_cast<amrex::Real>(m_base);

      if (max_dp > m_dxn) {

        std::vector<std::string> regions;
        amrex::ParmParse pp("ic");
        pp.queryarr("regions", regions);
        amrex::Print() << "\n\n";
        amrex::Print() << "**************************************************************\n";
        amrex::Print() << "  Invalid particle initialization. The spacing for " << m_base << "-cube\n";
        amrex::Print() << "  packing is less than the maximum particle diameter. \n";
        amrex::Print() << "  Fix the inputs file.\n";
        amrex::Print() << "**************************************************************\n";
        amrex::Print() << "\n\n";
        amrex::Abort("Fix the inputs file.");
      }

      // indices
      amrex::IntVect seed_hi;

      m_seed_lo[0] = a_bx.smallEnd(0)*m_base;
      m_seed_lo[1] = a_bx.smallEnd(1)*m_base;
      m_seed_lo[2] = a_bx.smallEnd(2)*m_base;

      seed_hi[0] = (a_bx.bigEnd(0) + 1)*m_base - 1;
      seed_hi[1] = (a_bx.bigEnd(1) + 1)*m_base - 1;
      seed_hi[2] = (a_bx.bigEnd(2) + 1)*m_base - 1;

      const amrex::Box bx(m_seed_lo, seed_hi);

      m_delta_bx[0] = amrex::max(0, seed_hi[0] - m_seed_lo[0] + 1);
      m_delta_bx[1] = amrex::max(0, seed_hi[1] - m_seed_lo[1] + 1);
      m_delta_bx[2] = amrex::max(0, seed_hi[2] - m_seed_lo[2] + 1);

      m_np = m_delta_bx[0] * m_delta_bx[1] * m_delta_bx[2];
    }

    int get_particles_number () const
    { return m_np; }

    template <amrex::RunOn >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect get_position (const int p_id,
                                amrex::RandomEngine const&) const
    {
      const int local_j = int(amrex::Math::floor(amrex::Real(p_id) / amrex::Real(m_delta_bx[0]*m_delta_bx[2])));
      const int local_k = int(amrex::Math::floor(amrex::Real(p_id - local_j*m_delta_bx[0]*m_delta_bx[2]) / amrex::Real(m_delta_bx[0])));
      const int local_i = p_id - local_j*m_delta_bx[0]*m_delta_bx[2] - local_k*m_delta_bx[0];

      const int i = local_i + m_seed_lo[0];
      const int j = local_j + m_seed_lo[1];
      const int k = local_k + m_seed_lo[2];

      const amrex::Real x = m_plo[0] + (i + 0.5) * m_dxn;
      const amrex::Real y = m_plo[1] + (j + 0.5) * m_dxn;
      const amrex::Real z = m_plo[2] + (k + 0.5) * m_dxn;

      return amrex::RealVect(x, y, z);
    }


};

#endif
