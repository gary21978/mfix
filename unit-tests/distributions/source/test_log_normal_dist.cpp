#include <test_log_normal_dist.H>

using namespace amrex;

int
test_log_normal_dist::test () {
  ParmParse pp("log-normal");

  std::string log_normal_dist = "log-normal";

  Print() << "\nLog-normal distribution tests:\n";

  // Checks on the validity of input values
  { std::string field = "log-normal.solid0";

    pp.add("solid0.diameter",    log_normal_dist);

    Real const valid_mu    = std::log(1000.0e-6);
    Real const valid_sigma = 1.04;

    Real const valid_dmin =   50.0e-6;
    Real const valid_dmax = 5000.0e-6;

    // Base inputs are valid -- make sure that holds.
    Print() << "Log-normal distribution test 01";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      pp.add("solid0.diameter.mean",   valid_mu);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      pp.add("solid0.diameter.std", valid_sigma);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      pp.add("solid0.diameter.min",  valid_dmin);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      pp.add("solid0.diameter.max",  valid_dmax);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      pp.add("solid0.diameter.type", Vw_dtype);

      // Base inputs are all in the database and valid -- make sure that holds.
      if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }

      if ( diameter.is_constant()   ) { return FAIL; }
      if ( diameter.is_normal()     ) { return FAIL; }
      if (!diameter.is_log_normal() ) { return FAIL; }
      if ( diameter.is_uniform()    ) { return FAIL; }
      if ( diameter.is_custom()     ) { return FAIL; }
    }
    Print() << " PASSED!\n";

    // Flag error if the number of bins is invalid
    Print() << "Log-normal distribution test 02";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      pp.add("solid0.diameter.bins", 0);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }

      pp.remove("solid0.diameter.bins");
      pp.add("solid0.diameter.bins", 64);
      if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }
    }
    Print() << " PASSED!\n";

    // Error for invalid min/max
    Print() << "Log-normal distribution test 03";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      // (max < min) invalid
      pp.remove("solid0.diameter.min");
      pp.add("solid0.diameter.min", 2.*valid_dmax);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }

      // (max == min) invalid
      pp.remove("solid0.diameter.min");
      pp.add("solid0.diameter.min", valid_dmax);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }

      // (min < 0.) invalid
      pp.remove("solid0.diameter.min");
      pp.add("solid0.diameter.min", -valid_dmin);
      if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }

      // Rest back to good state.
      pp.remove("solid0.diameter.min");
      pp.add("solid0.diameter.min", valid_dmin);
      if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }
    }
    Print() << " PASSED!\n";
  }

  Real const Nw_mu    = std::log(1000.0e-6);
  Real const Nw_sigma = 1.04;

  Real const Nw_min = std::exp(Nw_mu - 5.*Nw_sigma);
  Real const Nw_max = std::exp(Nw_mu + 5.*Nw_sigma);

  Real const Vw_mu    = Nw_mu + 3.*Nw_sigma*Nw_sigma;
  Real const Vw_sigma = Nw_sigma;

  Real const Vw_min = std::exp(Vw_mu - 5.*Vw_sigma);
  Real const Vw_max = std::exp(Vw_mu + 5.*Vw_sigma);

  { std::string field = "log-normal.solid1";

    Print() << "Log-normal distribution test 10";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(1);

      pp.add("solid1.diameter.bins", 2048);

      pp.add("solid1.diameter", log_normal_dist);
      pp.add("solid1.diameter.mean",      Nw_mu);
      pp.add("solid1.diameter.std",    Nw_sigma);
      pp.add("solid1.diameter.min",      Nw_min);
      pp.add("solid1.diameter.max",      Nw_max);
      pp.add("solid1.diameter.type",   Nw_dtype);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      Real const m = std::exp(Nw_mu +0.5*Nw_sigma*Nw_sigma);
      Real const s = std::sqrt(m*m*(std::exp(Nw_sigma*Nw_sigma)-1.));

      if (!compare(diameter.get_mean(),     m, 0.)) { return FAIL; }
      if (!compare(diameter.get_stddev(),   s, 0.)) { return FAIL; }
      if (!compare(diameter.get_min(), Nw_min, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), Nw_max, 0.)) { return FAIL; }

      if (!diameter.is_number_weighted()) { return FAIL; }
      if ( diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << m << "  " << mean
                << "\nstdd: " << s << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(m - mean) > 5.0e-4*m) { return FAIL; }
        if (amrex::Math::abs(s - stdev) > 1.0e-2*s) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    // Reprocess the previous inputs with volume weighted sampling enabled.
    // The min and max values need to be increased given the volume-weighted
    // distribution is much wider.
    Print() << "Log-normal distribution test 11";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);

      pp.remove("solid1.diameter.min");
      pp.add("solid1.diameter.min", Vw_min);

      pp.remove("solid1.diameter.max");
      pp.add("solid1.diameter.max", Vw_max);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      Real const m = std::exp(Vw_mu +0.5*Vw_sigma*Vw_sigma);
      Real const s = std::sqrt(m*m*(std::exp(Vw_sigma*Vw_sigma)-1.));

      if (!compare(diameter.get_mean(),     m, 1.)) { return FAIL; }
      if (!compare(diameter.get_stddev(),   s, 1.)) { return FAIL; }
      if (!compare(diameter.get_min(), Vw_min, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), Vw_max, 0.)) { return FAIL; }

      if (!diameter.is_number_weighted()) { return FAIL; }
      if ( diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << m << "  " << mean
                << "\nstdd: " << s << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(m - mean) > 5.0e-3*m) { return FAIL; }
        if (amrex::Math::abs(s - stdev) > 1.0e-2*s) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";
  } // end number-weighted checks

  // volume-weighted checks
  { std::string field = "log-normal.solid2";

    Print() << "Log-normal distribution test 20";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);

      pp.add("solid2.diameter.bins", 2048);

      pp.add("solid2.diameter", log_normal_dist);
      pp.add("solid2.diameter.mean",      Vw_mu);
      pp.add("solid2.diameter.std",    Vw_sigma);
      pp.add("solid2.diameter.min",      Vw_min);
      pp.add("solid2.diameter.max",      Vw_max);
      pp.add("solid2.diameter.type",   Vw_dtype);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      Real const m = std::exp(Vw_mu +0.5*Vw_sigma*Vw_sigma);
      Real const s = std::sqrt(m*m*(std::exp(Vw_sigma*Vw_sigma)-1.));

      if (!compare(diameter.get_mean(),     m, 0.)) { return FAIL; }
      if (!compare(diameter.get_stddev(),   s, 0.)) { return FAIL; }
      if (!compare(diameter.get_min(), Vw_min, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), Vw_max, 0.)) { return FAIL; }

      if ( diameter.is_number_weighted()) { return FAIL; }
      if (!diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << m << "  " << mean
                << "\nstdd: " << s << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(m - mean) > 5.0e-4*m) { return FAIL; }
        if (amrex::Math::abs(s - stdev) > 1.0e-2*s) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    // Reprocess the previous inputs with number-weighted sampling.
    Print() << "Log-normal distribution test 21";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(1);

      pp.remove("solid2.diameter.min");
      pp.add("solid2.diameter.min", Nw_min);

      pp.remove("solid2.diameter.max");
      pp.add("solid2.diameter.max", Nw_max);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      Real const m = std::exp(Nw_mu +0.5*Nw_sigma*Nw_sigma);
      Real const s = std::sqrt(m*m*(std::exp(Nw_sigma*Nw_sigma)-1.));

      if (!compare(diameter.get_mean(),     m, 1.)) { return FAIL; }
      if (!compare(diameter.get_stddev(),   s, 1.)) { return FAIL; }
      if (!compare(diameter.get_min(), Nw_min, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), Nw_max, 0.)) { return FAIL; }

      if ( diameter.is_number_weighted()) { return FAIL; }
      if (!diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << m << "  " << mean
                << "\nstdd: " << s << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(m - mean) > 5.0e-4*m) { return FAIL; }
        if (amrex::Math::abs(s - stdev) > 1.0e-2*s) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";
  } // end volume-weighted checks

  return PASS;
}
