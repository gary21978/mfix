#include <AMReX_ParmParse.H>

#include <mfix_tag.H>
#include <mfix_reporter.H>

using namespace amrex;


void MFIXTag::
Initialize ( int const a_max_level,
             int const a_solve_fluid,
             int const a_has_particles,
             const MFIXRegions& a_regions )
{

  ParmParse pp_tag("mfix.tag");

  { Vector<Real> values;
    if ( pp_tag.queryarr("rho", values ) ) {
      Real last = values.back();
      values.resize(a_max_level+1, last);

      bool less_than = false;
      if ( pp_tag.queryAdd("rho.less_than", less_than) )

      rho.set(true, less_than, values);
      m_tags++;

    } else {
      values.resize(a_max_level+1, std::numeric_limits<Real>::max());
      rho.set(false, false, values);
    }
  }

  { Vector<Real> values;
    if ( pp_tag.queryarr("grad_rho", values ) ) {
      Real last = values.back();
      values.resize(a_max_level+1, last);

      bool less_than = false;
      if ( pp_tag.queryAdd("grad_rho.less_than", less_than) )

      grad_rho.set(true, less_than, values);
      m_tags++;
    } else {
      values.resize(a_max_level+1, std::numeric_limits<Real>::max());
      grad_rho.set(false, false, values);
    }
  }

  { Vector<Real> values;
    if ( pp_tag.queryarr("vel", values ) ) {
      Real last = values.back();
      values.resize(a_max_level+1, last);

      bool less_than = false;
      if ( pp_tag.queryAdd("vel.less_than", less_than) )

      vel.set(true, less_than, values);
      m_tags++;
    } else {
      values.resize(a_max_level+1, std::numeric_limits<Real>::max());
      vel.set(false, false, values);
    }
  }

  { Vector<Real> values;
    if ( pp_tag.queryarr("grad_vel", values ) ) {
      Real last = values.back();
      values.resize(a_max_level+1, last);

      bool less_than = false;
      if ( pp_tag.queryAdd("grad_vel.less_than", less_than) )

      grad_vel.set(true, less_than, values);
      m_tags++;
    } else {
      values.resize(a_max_level+1, std::numeric_limits<Real>::max());
      grad_vel.set(false, false, values);
    }
  }

#if 0
  { Vector<Real> values;
    if ( pp_tag.queryarr("vorticity", values ) ) {
      Real last = values.back();
      values.resize(a_max_level+1, last);

      bool less_than = false;
      if ( pp_tag.queryAdd("vorticity.less_than", less_than) )

      vorticity.set(true, less_than, values);
      m_tags++;
    } else {
      values.resize(a_max_level+1, std::numeric_limits<Real>::max());
      vorticity.set(false, false, values);
    }
  }
#endif

  { Vector<std::string> tag_regions;
    if ( pp_tag.queryarr("regions", tag_regions ) ) {

      m_num_regions = tag_regions.size();
      h_regions.resize(num_regions());

      for (int rcv(0); rcv < num_regions(); rcv++) {

        RealBox const* rbox = a_regions.getRegion(tag_regions[rcv]);

        if ( rbox == nullptr ) {
          reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Unable to find tag region in defined regions!\n"
              << "Tag region name: " << tag_regions[rcv];
        } else {
          h_regions[rcv] = *rbox;
        }
        m_tags++;
      }
    }
  }


  if (m_tags > 0) {

    if (a_max_level == 0) {
     reporter::Log(reporter::Error,__FILE__, __LINE__)
         << "AMR tags are defined but amr.max_level = 0!";
    }
    if (!a_solve_fluid) {
     reporter::Log(reporter::Error,__FILE__, __LINE__)
         << "AMR tags are defined but the fluid is not solved!";
    }
  }

  if (a_max_level > 0 && a_has_particles) {
   reporter::Log(reporter::Error,__FILE__, __LINE__)
       << "AMR + particle support is not yet supported!";
  }
}
