#ifndef _MFIX_AVG_PARTICLE_PARMS_H_
#define _MFIX_AVG_PARTICLE_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_GpuContainers.H>

#include <mfix_eb_parms.H>
#include <mfix_fluid.H>


/**
 * @brief Struct to identify particle data averaged to the grid
 *
 * This structure contains indexes for particle quantities that are
 * model dependent. The data is stored as a struct of arrays (SoA).
 */
struct MFIXAvgParticleParms {

  public:


    void Initialize ( int const a_pc_var_count,
                      MFIXFluidPhase& a_fluid,
                      MFIXEmbeddedBoundaries& a_eb_parms);


    int nComp () const { return m_ncomp; }

    amrex::Gpu::HostVector<int> const& comps () const {
      return m_h_depComp;
    }

    int comp ( int const a_comp ) const
    { return (m_pc_var_count > 0) ? m_comp_map.at(a_comp) : -1;}

  private:

    int m_ncomp = 0;
    int m_pc_var_count = -1;

    amrex::Gpu::HostVector<int> m_h_depComp;
    amrex::Vector<int> m_comp_map;

    void setDepComp ( int const a_comp ) {
      AMREX_ASSERT( m_pc_var_count > 0 );
      AMREX_ASSERT( a_comp < m_pc_var_count );

      if (!m_h_depComp[a_comp]) {
        m_h_depComp[a_comp] = 1;
        m_comp_map[a_comp] = m_ncomp;
        m_ncomp++;
      }
    }

};

#endif
