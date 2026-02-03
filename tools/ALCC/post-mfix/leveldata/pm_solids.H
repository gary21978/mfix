#ifndef POST_MFIX_SOLIDS_
#define POST_MFIX_SOLIDS_

#include <AMReX_AmrMesh.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <pm_particles.H>

class pm_solids {

  public:

    pm_solids ( int const a_max_level,
                amrex::Vector<amrex::Geometry>            const & a_geom,
                amrex::Vector<amrex::DistributionMapping> const & a_dmap,
                amrex::Vector<amrex::BoxArray>            const & a_ba,
                pm_particles* a_particles);

    ~pm_solids ();

    void reset ();

    int contains ( std::string a_var ) noexcept
    { return (m_variable_map.count(a_var)) ? 1 : 0; }

    int index ( std::string a_var ) noexcept {
      if (m_variable_map.count(a_var))
      { return m_variable_map[a_var]; }
      return -1;
    }

    int nComp () const noexcept { return m_comps; }

    std::string varName (int const a_comp) const noexcept
    { return m_variables[a_comp]; }

    amrex::Vector<std::string> varName () const noexcept
    { return m_variables; }

    amrex::MultiFab* get_data (const int a_lev) noexcept
    { return m_data[a_lev].get(); }


    amrex::MultiFab const* get_const_data (const int a_lev) noexcept
    { return m_data[a_lev].get(); }


    amrex::Vector< amrex::MultiFab* > get_data () const noexcept {
      amrex::Vector<amrex::MultiFab*> r;
      r.reserve(m_data.size());
      for (auto& lev_data : m_data) { r.push_back( lev_data.get() ); }
      return r;
    }

  private:

    int const m_ngrow = 1;

    int m_max_level;

    amrex::Vector<amrex::Geometry>            const m_geom;
    amrex::Vector<amrex::DistributionMapping> const m_dmap;
    amrex::Vector<amrex::BoxArray>            const m_grids;

    std::map<std::string, std::string> m_LP_to_EL_map;

    int m_comps;
    amrex::Vector<std::string> m_variables;

    std::map<std::string, int> m_variable_map;

    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_data;

    int get_nlev ( ) noexcept { return m_max_level + 1; }
};

#endif
