#include <AMReX_ParmParse.H>

#include <fml_solids.H>
#include <fml_diffusion.H>

using namespace amrex;

fml_solids::
fml_solids ( int const a_max_level,
            Vector<Geometry>            const & a_geom,
            Vector<DistributionMapping> const & a_dmap,
            Vector<BoxArray>            const & a_grids,
            fml_particles* a_particles,
            Vector<std::string> a_vars)
  : fml_derived(a_max_level, a_geom, a_dmap, a_grids)
{

  m_LP_to_EL_map["alpha_p"] = "volume";
  m_LP_to_EL_map["Up"] = "velx";
  m_LP_to_EL_map["Vp"] = "vely";
  m_LP_to_EL_map["Wp"] = "velz";

  if ( std::find(a_vars.begin(), a_vars.end(), "beta") != a_vars.end() )
  { add_variable("alpha_p,Up,Vp,Wp"); }

  if ( std::find(a_vars.begin(), a_vars.end(), "slip_vel") != a_vars.end() ) {
    add_variable("Up,Vp,Wp");
  } else {
    if ( std::find(a_vars.begin(), a_vars.end(), "Up") != a_vars.end()) { add_variable("Up"); }
    if ( std::find(a_vars.begin(), a_vars.end(), "Vp") != a_vars.end()) { add_variable("Vp"); }
    if ( std::find(a_vars.begin(), a_vars.end(), "Wp") != a_vars.end()) { add_variable("Wp"); }
  }

  if ( std::find(a_vars.begin(), a_vars.end(), "alpha_p") != a_vars.end() ||
       std::find(a_vars.begin(), a_vars.end(), "alpha_f") != a_vars.end() )
  { add_variable("alpha_p"); }

  // Number of 'base' variables
  m_comps = static_cast<int>(m_variables.size());

  m_data.resize(get_nlev());
  for (int lev(0); lev < get_nlev(); ++lev) {

    m_data[lev].reset(new MultiFab(m_grids[lev], m_dmap[lev], m_comps, m_ngrow));
    m_data[lev]->setVal(0.0, 0, m_comps, m_ngrow);

    for (int comp(0); comp<m_comps; ++comp) {

      if ( m_variables[comp].find('*') != std::string::npos) {

          auto const pos = m_variables[comp].find('*');

          std::string const EL_var1 = m_variables[comp].substr(0, pos);
          std::string const LP_var1 = m_LP_to_EL_map[EL_var1];
          AMREX_ALWAYS_ASSERT( a_particles->contains(LP_var1) );

          std::string const EL_var2 = m_variables[comp].substr(pos+1);
          std::string const LP_var2 = m_LP_to_EL_map[EL_var2];
          AMREX_ALWAYS_ASSERT( a_particles->contains(LP_var2) );

          Print() << "  depositing " << LP_var1+"*"+LP_var2
            << " to " << EL_var1+"*"+EL_var2 << "\n";

          a_particles->deposit_mult(lev, m_data[lev].get(), comp, LP_var1, LP_var2);

      } else if ( m_variables[comp].find('^') != std::string::npos) {

          auto const pos = m_variables[comp].find('^');

          std::string const EL_var = m_variables[comp].substr(0, pos);
          std::string const LP_var = m_LP_to_EL_map[EL_var];
          AMREX_ALWAYS_ASSERT( a_particles->contains(LP_var) );

          int power = std::stoi( m_variables[comp].substr(pos+1) );

          Print() << "  depositing std::pow(" << LP_var+","+std::to_string(power)+")"
            << " to " << EL_var+"^"+std::to_string(power) << "\n";

          a_particles->deposit_pow(lev, m_data[lev].get(), comp, LP_var, 1.0, power);


      } else {

          std::string const EL_var = m_variables[comp];
          std::string const LP_var = m_LP_to_EL_map[EL_var];

          Print() << "  depositing " << LP_var << " to " << EL_var << "\n";

          if (toLower(LP_var).compare("volume") )
          { AMREX_ALWAYS_ASSERT( a_particles->contains(LP_var) ); }

          a_particles->deposit(lev, m_data[lev].get(), comp, LP_var);
      }

    }
    m_data[lev]->FillBoundary(m_geom[lev].periodicity());
  }

  int const alpha_p_comp = m_variable_map["alpha"];

  // Initialize the diffuison operator.
  fml_diffusion diffusion( m_max_level, m_geom, m_dmap, m_grids);

  diffusion.reset( GetVecOfConstPtrs(m_data), alpha_p_comp,
    a_particles->get_mean_diameter());

  // Smooth (diffuse) deposited particle data
  diffusion.smooth( GetVecOfPtrs(m_data), 0, m_comps, m_variables);
}


void fml_solids::
reset ()
{
  int const alpha_p_comp = m_variable_map["alpha"];

  for (int lev(0); lev<get_nlev(); ++lev) {

    for (MFIter mfi(*m_data[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& tbox = mfi.tilebox();

      Array4<Real const> const& alpha = m_data[lev]->const_array(mfi,alpha_p_comp);
      Array4<Real      > const& sdata = m_data[lev]->array(mfi);

      ParallelFor(tbox, m_comps, [alpha, alpha_p_comp, sdata]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { if ((n != alpha_p_comp) && alpha(i,j,k) > 0.)
        { sdata(i,j,k,n) /= alpha(i,j,k); }
      });
    }
  }
}
