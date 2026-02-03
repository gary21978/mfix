#ifndef POST_MFIX_PARTICLE_CONTAINER_
#define POST_MFIX_PARTICLE_CONTAINER_

#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_ParIter.H>

#include <pm_plotfile.H>

class PMParIter : public amrex::ParIter<0,0,0,0>
{
  public:

    using amrex::ParIter<0,0,0,0>::ParIter;
    using RealVector = amrex::ParIter<0,0,0,0>::ContainerType::RealVector;

};

struct PMCheckPair
{
  amrex::Real neighborhood;
  PMCheckPair (amrex::Real neighborhood_) : neighborhood(neighborhood_) {}

  template <int NSR, int NSI, int NAR, int NAI>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  bool operator() (const amrex::ConstParticleTileData<
    amrex::Particle<NSR, NSI>, NAR, NAI>& ptile, int i, int j) const
  {
    const auto p1 = ptile.getSuperParticle(i);
    const auto p2 = ptile.getSuperParticle(j);

    if (i < j)
      if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
        return true;

    return false;
  }
};


class pm_particles : public amrex::NeighborParticleContainer<0,0,0,0>
{

  public:

    pm_particles ( int const a_max_level,
                   amrex::Vector<amrex::Geometry>            const & a_geom,
                   amrex::Vector<amrex::DistributionMapping> const & a_dmap,
                   amrex::Vector<amrex::BoxArray>            const & a_ba,
                   pm_plotfile* a_plotfile);

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

      //else if (m_variable_map_i.count(a_var))
      //{ return m_variable_map_i[a_var]; }

      else {
        std::string message;
        message = "Variable " + a_var + " not found in particle plot file!\n";
        amrex::Abort(message);
      }
      return -1;
    }

    amrex::Real get_mean_diameter () noexcept { return m_mean_diameter; }

    void deposit ( int const a_lev, amrex::MultiFab* a_MF, int const a_dstcomp,
                   std::string const a_var);

    void deposit_mult ( int const a_lev, amrex::MultiFab* a_MF, int const a_dstcomp,
                        std::string const a_var1, std::string const a_var2);

    void deposit_pow ( int const a_lev, amrex::MultiFab* a_MF, int const a_dstcomp,
                       std::string const a_var, amrex::Real const a_a, int const a_b);

  protected:

    int m_max_level;

    int m_verbose;

  private:

    int get_nlev () { return m_max_level + 1; }


    int const m_real_comps;
    int const m_int_comps;

    amrex::Vector<std::string> m_int_comp_names;
    amrex::Vector<std::string> m_real_comp_names;

    std::map<std::string, int> m_variable_map_r;
    std::map<std::string, int> m_variable_map_i;

    amrex::Real m_mean_diameter;

};
#endif
