#include <fml_derived.H>

using namespace amrex;

fml_derived::
~fml_derived ()
{
  for (auto& lev_data : m_data) { lev_data.reset( nullptr ); }
}

fml_derived::
fml_derived ( int const a_max_level,
             amrex::Vector<amrex::Geometry>            const & a_geom,
             amrex::Vector<amrex::DistributionMapping> const & a_dmap,
             amrex::Vector<amrex::BoxArray>            const & a_grids)
  : m_max_level(a_max_level)
  , m_geom(a_geom)
  , m_dmap(a_dmap)
  , m_grids(a_grids)
{
  m_data.resize(get_nlev());
}

void
fml_derived::
add_variable ( std::string a_var ) {

  // Only a single component.
  if( a_var.find(",") == std::string::npos) {

    add_single_variable(a_var);

  // Comma separated list of variables.
  } else {

    std::replace(a_var.begin(), a_var.end(), ',', ' ');

    std::istringstream var_iss(a_var);
    std::istream_iterator<std::string> var_comp_it(var_iss);

    amrex::Vector<std::string> vars;

    std::copy(var_comp_it, std::istream_iterator<std::string>(),
               std::back_inserter(vars));

    for (int comp(0); comp<vars.size(); ++comp) {
      add_single_variable(vars[comp]);
    }
  }
}


void
fml_derived::
add_single_variable ( std::string a_var ) {
  if (std::find(m_variables.begin(), m_variables.end(), a_var) == m_variables.end() ) {
    m_variables.push_back(a_var);
    m_variable_map[a_var] = m_variables.size()-1;
  }
}
