#ifndef FILTER_ML_DERIVED_
#define FILTER_ML_DERIVED_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

class fml_derived {

  public:

    fml_derived ( int const a_max_level,
                 amrex::Vector<amrex::Geometry>            const & a_geom,
                 amrex::Vector<amrex::DistributionMapping> const & a_dmap,
                 amrex::Vector<amrex::BoxArray>            const & a_grids);

    ~fml_derived ();

    int contains ( std::string a_var ) noexcept
    { return (m_variable_map.count(a_var)) ? 1 : 0; }

    int index ( std::string a_var ) noexcept {
      if (m_variable_map.count(a_var))
      { return m_variable_map[a_var]; }
      return -1;
    }

    int nComp () const noexcept { return m_comps; }

    // Access to variable names
    std::string varName (int const a_comp) const noexcept
    { return m_variables[a_comp]; }

    amrex::Vector<std::string> varNames () const noexcept
    { return m_variables; }

    // Access to MultiFab data
    amrex::MultiFab* get_data (const int a_lev) noexcept
    { return m_data[a_lev].get(); }

    amrex::Vector< amrex::MultiFab* > get_data () const noexcept {
      amrex::Vector<amrex::MultiFab*> r;
      r.reserve(m_data.size());
      for (auto& lev_data : m_data) { r.push_back( lev_data.get() ); }
      return r;
    }

    amrex::MultiFab const* get_const_data (const int a_lev) noexcept
    { return m_data[a_lev].get(); }

    amrex::Vector< amrex::MultiFab const* > get_const_data () const noexcept {
      amrex::Vector<amrex::MultiFab const*> r;
      r.reserve(m_data.size());
      for (auto& lev_data : m_data) { r.push_back( lev_data.get() ); }
      return r;
    }

  protected:

    int const m_ngrow = 1;

    int m_max_level;

    amrex::Vector<amrex::Geometry>            const m_geom;
    amrex::Vector<amrex::DistributionMapping> const m_dmap;
    amrex::Vector<amrex::BoxArray>            const m_grids;

    int m_comps;
    amrex::Vector<std::string> m_variables;

    std::map<std::string, int> m_variable_map;

    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_data;

    int get_nlev ( ) noexcept { return m_max_level + 1; }

    int var_index ( std::string a_var ) noexcept {
      if (m_variable_map.count(a_var))
      { return m_variable_map[a_var]; }
      return -1;
    }

    amrex::MultiFab* define (const int a_lev, const int a_ncomp) noexcept
    { AMREX_ALWAYS_ASSERT( m_data[a_lev] == nullptr );
      m_data[a_lev].reset(new amrex::MultiFab(m_grids[a_lev], m_dmap[a_lev], a_ncomp, 1));
      m_data[a_lev]->setVal(0.0);
      return m_data[a_lev].get();
    }

    void add_variable ( std::string );

  private:

    void add_single_variable ( std::string a_var );
};

#endif
