#ifndef FILTER_ML_FLUID_
#define FILTER_ML_FLUID_

#include <fml_derived.H>
#include <fml_plotfile.H>

class fml_fluid : public fml_derived {

  public:

    fml_fluid ( int const a_max_level,
               amrex::Vector<amrex::Geometry>            const & a_geom,
               amrex::Vector<amrex::DistributionMapping> const & a_dmap,
               amrex::Vector<amrex::BoxArray>            const & a_grids,
               fml_plotfile* a_plotfile);

};
#endif
