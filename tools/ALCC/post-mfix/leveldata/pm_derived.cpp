#include <pm_derived.H>

using namespace amrex;

pm_derived::
~pm_derived ()
{
  for (auto& lev_data : m_data) { lev_data.reset( nullptr ); }
}

pm_derived::
pm_derived ( int const a_max_level,
             amrex::Vector<amrex::Geometry>            const & a_geom,
             amrex::Vector<amrex::DistributionMapping> const & a_dmap,
             amrex::Vector<amrex::BoxArray>            const & a_grids,
             pm_fluid* /*a_fuild*/, pm_solids* /*a_solids*/)
  : m_max_level(a_max_level)
  , m_geom(a_geom)
  , m_dmap(a_dmap)
  , m_grids(a_grids)
{
  m_data.resize(get_nlev());
}


void
pm_derived::
add_variable ( std::string a_variable, int const a_comp ) {

  m_variables.push_back(a_variable);
  m_variable_map[a_variable] = a_comp;
}
