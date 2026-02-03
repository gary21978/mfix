#ifndef FILTER_ML_
#define FILTER_ML_

#include <regex>

#include <fml_plotfile.H>
#include <fml_particles.H>
#include <fml_fluid.H>
#include <fml_solids.H>

#include <fml_derived.H>

class FilterML
{
  public:

    FilterML (std::string a_plotifle);
    ~FilterML ();

    fml_particles* get_particles () noexcept { return m_particles; }

    void create_Eulerian_solids ( amrex::Vector<std::string> a_vars);

    void compute_Eulerian_filter (int const a_filter_size);

    amrex::MultiFab const* get_const_alpha_p (const int a_lev) noexcept
    { return m_alpha_p[a_lev].get(); }

    amrex::Vector< amrex::MultiFab const* > get_const_alpha_p () const noexcept {
      amrex::Vector<amrex::MultiFab const*> r;
      r.reserve(m_alpha_p.size());
      for (auto& lev_alpha_p : m_alpha_p) { r.push_back( lev_alpha_p.get() ); }
      return r;
    }

    amrex::MultiFab const* get_const_alpha_f (const int a_lev) noexcept
    { return m_alpha_f[a_lev].get(); }

    amrex::Vector< amrex::MultiFab const* > get_const_alpha_f () const noexcept {
      amrex::Vector<amrex::MultiFab const*> r;
      r.reserve(m_alpha_f.size());
      for (auto& lev_alpha_f : m_alpha_f) { r.push_back( lev_alpha_f.get() ); }
      return r;
    }

  private:

    int const m_verbose;

    std::string m_plotfile;

    int m_max_level;

    amrex::Vector<amrex::Geometry> m_geom;
    amrex::Vector<amrex::BoxArray> m_grids;
    amrex::Vector<amrex::DistributionMapping> m_dmap;

    fml_fluid* m_fluid;

    fml_particles* m_particles;

    fml_solids* m_solids;

    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_alpha_p;

    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_alpha_f;

    amrex::Vector<int> m_ref_ratio;

    int get_nlev () noexcept { return m_max_level + 1; }

    void fill_var_MF ( int const a_lev, std::string const a_var,
                       amrex::MultiFab* a_MF, int const a_dstcomp=0);

    amrex::Vector<std::string> get_var_comps ( std::string const a_var ) noexcept
    {
      std::string var(a_var);

      amrex::Vector<std::string> var_comps;

      if (var.size() > 0) {
        // Only a single component.
        if( var.find("*") == std::string::npos) {

          var_comps.push_back(var);

        } else { // we need to split this into pieces.

          std::replace(var.begin(), var.end(), '*', ' ');

          std::istringstream var_iss(var);
          std::istream_iterator<std::string> var_comp_it(var_iss);

          std::copy(var_comp_it, std::istream_iterator<std::string>(),
                     std::back_inserter(var_comps));
        }
      }
      return var_comps;
    }

    std::string csv_filename (int const a_lev, int const a_filter_size) noexcept {

      // Find the start of the plotfile name (e.g., plt)
      size_t last = m_plotfile.find_last_of("/\\");
      last = (last == std::string::npos) ? 0 : last + 1;

      // Name of plot file with step index but excluding path
      std::string base = m_plotfile.substr(last,m_plotfile.size()-last);

      std::string filename = "fml_" + base +
        amrex::Concatenate("_lev",a_lev,1) +
        amrex::Concatenate("_fs",a_filter_size,4) + ".csv";

      return filename;
    }

};
#endif
