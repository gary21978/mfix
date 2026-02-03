#ifndef POST_MFIX_
#define POST_MFIX_

#include <regex>
#include <AMReX_AmrCore.H>

#include <pm_plotfile.H>
#include <pm_particles.H>
#include <pm_fluid.H>
#include <pm_solids.H>

#include <pm_derived.H>

class post_mfix : public amrex::AmrCore
{
  public:

    //post_mfix (std::string a_plot_file_name);
    post_mfix (amrex::Geometry const& a_level_0_geom,
               amrex::AmrInfo const& a_amr_info,
               pm_plotfile* m_p_plotfile);

    ~post_mfix ();

    pm_particles* get_particles () noexcept { return m_particles; }

    void compute_derived_quantites ();

    amrex::Vector< amrex::Vector<amrex::Real> >
      compute_averages (amrex::Vector<std::string> a_vars);

    amrex::Vector<amrex::Vector<amrex::Real>>
      compute_fluctuations (amrex::Vector<std::string> a_vars,
              amrex::Vector<amrex::Vector<amrex::Real>> a_avgs);

    void output_hdf5 ( std::string a_plotfile );

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

    void compute_grad ( int const a_lev, int const a_dir,
                        std::string a_var, amrex::MultiFab* a_MF);

    void compute_stress ( int const a_lev, std::string a_stress,
                          amrex::MultiFab* a_MF);

    void compute_granular_viscosity ( int const a_lev,
                          std::string a_term,
                          amrex::MultiFab* a_MF);

    amrex::Real calc_fluct ( amrex::MultiFab const* a_MF,
                             int const a_comp, amrex::Real const a_avg );

    void fill_var_MF ( int const a_lev, std::string const a_var,
                       amrex::MultiFab* a_MF, int const a_dstcomp=0);
  private:

    int const m_verbose;

    pm_fluid* m_fluid;

    pm_particles* m_particles;

    pm_solids* m_solids;

    //pm_derived* m_gradients;

    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_alpha_p;

    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_alpha_f;

    int get_nlev ( ) noexcept { return maxLevel() + 1; }

    int grad_dir ( std::string a_var) noexcept
    {
      if ( a_var.find("_x") != std::string::npos ) { return 0; }
      if ( a_var.find("_y") != std::string::npos ) { return 1; }
      if ( a_var.find("_z") != std::string::npos ) { return 2; }
      return -1;
    }

    int is_gradient ( std::string a_var ) noexcept
    { // Skip if its a known gradient (e.g., fluid pressure)
      if (m_fluid->contains(a_var)) { return 0; }
      else if (m_solids->contains(a_var)) { return 0; }

      // For grad_x( ), grad_y( ), grad_z( ) format
      std::regex grad_comp("(grad_[xyz]\\(.*?\\))");
      std::smatch match;

      std::regex_search(a_var, match, grad_comp);
      return match.size();
    }

    amrex::Vector<std::string> get_var_comps ( std::string const a_var ) noexcept
    {
      std::string var(a_var);

      amrex::Vector<std::string> var_comps;

      // For grad_x( ), grad_y( ), grad_z( ) format
      std::regex grad_comp("(grad_[xyz]\\(.*?\\))");
      std::smatch match;

      // This will pull out the derivates
      while (std::regex_search(var, match, grad_comp)) {
        var_comps.push_back(match[0]);
        var = match.suffix().str();
      }

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

  protected:

    //! Tagging cells for refinement
    virtual void ErrorEst (int /*lev*/, amrex::TagBoxArray & /*tags*/,
                           amrex::Real /*time*/, int /*ngrow*/) override
    { amrex::Abort("ErrorEst: To be implemented"); }

    //! Make a new level using provided BoxArray and DistributionMapping and
    //! fill with interpolated coarse level data. Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int /*lev*/, amrex::Real /*time*/,
                                         const amrex::BoxArray & /*ba*/,
                                         const amrex::DistributionMapping & /*dm*/) override
    { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse data. Called
    //! by AmrCore::regrid.
    virtual void RemakeLevel (int /*lev*/, amrex::Real /*time*/,
                              const amrex::BoxArray & /*ba*/,
                              const amrex::DistributionMapping & /*dm*/) override
    { amrex::Abort("RemakeLevel: To be implemented"); }

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int /*lev*/) override
    { amrex::Abort("ClearLevel: To be implemented"); }

    void MakeNewLevelFromScratch (int lev, amrex::Real time,
                                  const amrex::BoxArray & new_grids,
                                  const amrex::DistributionMapping & new_dmap) override;

};
#endif
