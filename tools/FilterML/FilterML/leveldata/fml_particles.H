#ifndef FILTER_ML_PARTICLE_CONTAINER_
#define FILTER_ML_PARTICLE_CONTAINER_

#include <AMReX_Particles.H>
#include <AMReX_ParIter.H>

#include <fml_plotfile.H>

class fml_particles : public amrex::ParticleContainer<0,0>
{

  public:

    fml_particles ( int const a_max_level,
                   amrex::Vector<amrex::Geometry>            const & a_geom,
                   amrex::Vector<amrex::DistributionMapping> const & a_dmap,
                   amrex::Vector<amrex::BoxArray>            const & a_ba,
                   amrex::Vector<int>                        const & a_rr,
                   fml_plotfile* a_plotfile);

    int nComps_int () const noexcept { return m_int_comps; }

    int nComps_real () const noexcept { return m_real_comps; }

    amrex::Vector<std::string> VarNames_int () const noexcept
    { return m_int_comp_names; }

    amrex::Vector<std::string> VarNames_real () const noexcept
    { return m_real_comp_names; }

    amrex::Real max ( std::string const a_var,
                      bool const a_local = false);

    amrex::Real sum ( std::string const a_var,
                      bool const a_local = false);

    amrex::Real fluct ( std::string const a_var,
                        bool a_local = false )
    {
      amrex::Real avg = sum(a_var, 0);
      avg /= static_cast<amrex::Real>(TotalNumberOfParticles());

      amrex::Real loc_fluct = fluct(a_var, avg, a_local);
      return loc_fluct;
    }

    amrex::Real fluct ( std::string const a_var,
                        amrex::Real const a_avg,
                        bool a_local = false );


    int contains ( std::string a_var ) noexcept
    { return (m_variable_map_r.count(a_var)) ? 1 : 0; }

    int var_index ( std::string a_var ) noexcept {
      if (m_variable_map_r.count(a_var))
      { return m_variable_map_r[a_var]; }
      return -1;
    }

    amrex::Real get_mean_diameter () noexcept { return m_mean_diameter; }

    void deposit ( int const a_lev, amrex::MultiFab* a_MF, int const a_dstcomp,
                   std::string const a_var);

    void deposit_mult ( int const a_lev, amrex::MultiFab* a_MF, int const a_dstcomp,
                        std::string const a_var1, std::string const a_var2);

    void deposit_pow ( int const a_lev, amrex::MultiFab* a_MF, int const a_dstcomp,
                       std::string const a_var, amrex::Real const a_a, int const a_b);

  private:

    int m_max_level;

    int m_verbose;

    int const m_real_comps;
    int const m_int_comps;

    amrex::Vector<std::string> m_int_comp_names;
    amrex::Vector<std::string> m_real_comp_names;

    std::map<std::string, int> m_variable_map_r;
    std::map<std::string, int> m_variable_map_i;

    amrex::Real m_mean_diameter;

    int get_nlev () { return m_max_level + 1; }

};
#endif
