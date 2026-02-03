#ifndef MFIX_TEST_CUSTOM_DIST_H_
#define MFIX_TEST_CUSTOM_DIST_H_

#include <string>

#include <AMReX_REAL.H>

#include <test_dist.H>

class test_custom_dist : public test_dist
{

  public:

    test_custom_dist () = default;
    ~test_custom_dist () = default;

    int test ();

  private:

    void create_custom ( std::string const a_name, int const is_cdf, int const a_entries,
                         amrex::Real const* a_bins, amrex::Real const* a_cdf );

    void create_custom_normal ( std::string const a_name, int const a_is_cdf,
                                int const a_bins, int const a_reset_base,
                                amrex::Real a_min, amrex::Real a_max,
                                amrex::Real a_mean, amrex::Real a_stdev);

};

#endif
