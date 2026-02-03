#include <test_constant_dist.H>

using namespace amrex;

int test_constant_dist::test () {

  ParmParse pp("constant");

  std::string constant = "constant";

  Real const mu = 1000.0e-6;

  // Check for an invalid constant value.
  { std::string field = "constant.solid0";

    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      pp.add("solid0.diameter", constant);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      pp.add("solid0.diameter.constant", mu);

      // Base inputs are all in the database and valid -- make sure that holds.
      if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }

      if (!diameter.is_constant()   ) { return FAIL; }
      if ( diameter.is_normal()     ) { return FAIL; }
      if ( diameter.is_log_normal() ) { return FAIL; }
      if ( diameter.is_uniform()    ) { return FAIL; }
      if ( diameter.is_custom()     ) { return FAIL; }
    }

    // Flag error if the mean is invalid
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      pp.remove("solid0.diameter.constant");
      pp.add("solid0.diameter.constant", 0.);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }

      pp.remove("solid0.diameter.constant");
      pp.add("solid0.diameter.constant", -mu);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }

      pp.remove("solid0.diameter.constant");
      pp.add("solid0.diameter.constant", mu);

      if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }

      if (!compare(diameter.get_mean(),    mu, 0.)) { return FAIL; }
      if (!compare(diameter.get_min(),     mu, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(),     mu, 0.)) { return FAIL; }
      if (!compare(diameter.get_stddev(), 0.0, 0.)) { return FAIL; }

      { unsigned long const samples( 100000 );
        Real mean(0.), stddev(0.);

        if (compute_mean_and_std(diameter, samples, mean, stddev) == FAIL) { return FAIL; }

        if (amrex::Math::abs(mu - mean) > 1.0e-18) { return FAIL; }
        if (amrex::Math::abs(stddev)    > 1.0e-18) { return FAIL; }
      }
    }
  }

  return PASS;
}
