#ifndef FILTER_ML_SOLIDS_
#define FILTER_ML_SOLIDS_

#include <fml_derived.H>
#include <fml_particles.H>

class fml_solids : public fml_derived {

  public:

    fml_solids ( int const a_max_level,
                amrex::Vector<amrex::Geometry>            const & a_geom,
                amrex::Vector<amrex::DistributionMapping> const & a_dmap,
                amrex::Vector<amrex::BoxArray>            const & a_ba,
                fml_particles*              a_particles,
                amrex::Vector<std::string>  a_vars);

    void reset ();

  private:

    std::map<std::string, std::string> m_LP_to_EL_map;

};

#endif
