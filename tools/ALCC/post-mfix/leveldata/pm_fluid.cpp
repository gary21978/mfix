#include<AMReX.H>

#include <pm_fluid.H>

using namespace amrex;

pm_fluid::
~pm_fluid ()
{
  for (auto& lev_data : m_data) { lev_data.reset( nullptr ); }
}


pm_fluid::
pm_fluid ( int const a_max_level,
           amrex::Vector<amrex::Geometry>            const & a_geom,
           amrex::Vector<amrex::DistributionMapping> const & a_dmap,
           amrex::Vector<amrex::BoxArray>            const & a_grids,
           pm_plotfile* a_plotfile)
  : m_max_level(a_max_level)
  , m_geom(a_geom)
  , m_dmap(a_dmap)
  , m_grids(a_grids)
{

  m_variables = {"Uf", "Vf", "Wf", "Pf",
                 "grad_x(Pf)", "grad_y(Pf)", "grad_z(Pf)" };

  m_comps = static_cast<int>(m_variables.size());

  // create a map between variable name and component
  for (int lc(0); lc < m_comps; ++lc) {
    m_variable_map[m_variables[lc]] = lc;
  }

  // get variable names from plotfile
  Vector<std::string> plt_variables(a_plotfile->fluid_variables());
  int const plt_comps = static_cast<int>(plt_variables.size());

  // create a map between plot file variable name and component
  std::map<std::string, int> plt_variable_map;
  for (int lc(0); lc < plt_comps; ++lc) {
    plt_variable_map[plt_variables[lc]] = lc;
  }

  // create a map linking post-mfix variable names
  // and names in the plot file
  std::map<std::string, std::string> fluid_to_plt_map;

  fluid_to_plt_map["Uf"] = "u_g";
  fluid_to_plt_map["Vf"] = "v_g";
  fluid_to_plt_map["Wf"] = "w_g";
  fluid_to_plt_map["Pf"] = "p_g";
  fluid_to_plt_map["grad_x(Pf)"] = "gpx";
  fluid_to_plt_map["grad_y(Pf)"] = "gpy";
  fluid_to_plt_map["grad_z(Pf)"] = "gpz";

  m_data.resize(get_nlev());

  for (int lev(0); lev < get_nlev(); ++lev) {

    m_data[lev].reset(new MultiFab(m_grids[lev], m_dmap[lev], m_comps, m_ngrow));

    const MultiFab pf_data = a_plotfile->get_mf_data(lev);

    for (int comp(0); comp<m_comps; ++comp) {

      std::string f_var = m_variables[comp]; // variable to copy

      // Make sure we have that variable
      if ( fluid_to_plt_map.count(f_var) ) {

        std::string p_var = fluid_to_plt_map[f_var];

        int const srccomp = plt_variable_map[p_var];
        int const dstcomp = m_variable_map[f_var];

        // Copy mf into m_data copy 1 component from pf_data. No ghosts.
        (m_data[lev].get())->ParallelCopy(pf_data, srccomp, dstcomp, 1, 0, 0);
      }
      else {
        Print() << "Failed to find " << f_var << " in fluid_to_plt_map\n";
      }
    }
    m_data[lev]->FillBoundary(m_geom[lev].periodicity());

  }

}
