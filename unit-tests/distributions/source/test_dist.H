#ifndef MFIX_TEST_DIST_H_
#define MFIX_TEST_DIST_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>

#include <mfix_distributions.H>

static int constexpr PASS(0);
static int constexpr FAIL(1);

class test_dist {

  public:

    test_dist () = default;
    ~test_dist () = default;

    std::string const Nw_dtype = "number-weighted";
    std::string const Vw_dtype = "volume-weighted";

    bool compare (amrex::Real const x, amrex::Real const y, amrex::Real const scale=10.)
    { return (std::abs(x-y) <= scale*std::numeric_limits<amrex::Real>::epsilon()); }

    int compute_mean_and_std ( INPUT_DIST_t & a_diameter,
                               unsigned long const a_samples,
                               amrex::Real& a_mean,
                               amrex::Real& a_stddev,
                               int const a_print_pdf = 0,
                               int const a_bins = -1,
                               int const a_use_midpoint = 1);


    amrex::Real normal_mean_Nw_to_Vw ( amrex::Real a_Nw_mu, amrex::Real a_Nw_sigma )
    { amrex::Real const m(a_Nw_mu), s(a_Nw_sigma);
      return (3.*s*s*s*s + 6.*s*s*m*m + m*m*m*m) / (3.*s*s*m + m*m*m);
    }

    amrex::Real normal_stdev_Nw_to_Vw ( amrex::Real a_Nw_mu, amrex::Real a_Nw_sigma )
    { amrex::Real const Vw_mu(normal_mean_Nw_to_Vw(a_Nw_mu, a_Nw_sigma));
      amrex::Real const m(a_Nw_mu), s(a_Nw_sigma);
      amrex::Real const dm = (m - Vw_mu);
      return std::sqrt( ((3*s*s*s*s)*(3*m + 2*dm) + (m*m*m*dm*dm) +
      (s*s)*(m*m*m + 6.*m*m*dm + 3*m*dm*dm) ) / (3.*s*s*m + m*m*m));
    }
};

#endif

