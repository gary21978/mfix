#include <test_normal_dist.H>

using namespace amrex;

int test_normal_dist::test () {

  ParmParse pp("normal");

  std::string normal_dist = "normal";

  Print() << "\nNormal distribution tests:\n";

  // Checks on the validity of input values
  { std::string field = "normal.solid0";

    pp.add("solid0.diameter",    normal_dist);

    Real const valid_mu    = 1000.0e-6;
    Real const valid_sigma =  250.0e-6;

    Real const valid_dmin =  500.0e-6;
    Real const valid_dmax = 1500.0e-6;

    Print() << "Normal distribution test 01";
    // Base inputs are valid -- make sure that holds.
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
      pp.add("solid0.diameter.type",   Vw_dtype);

      // Base inputs are all in the database and valid -- make sure that holds.
      if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }

      if ( diameter.is_constant()   ) { return FAIL; }
      if (!diameter.is_normal()     ) { return FAIL; }
      if ( diameter.is_log_normal() ) { return FAIL; }
      if ( diameter.is_uniform()    ) { return FAIL; }
      if ( diameter.is_custom()     ) { return FAIL; }
    }
    Print() << " PASSED!\n";

    // Flag error if the number of bins is invalid
    Print() << "Normal distribution test 02";
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
    Print() << "Normal distribution test 03";
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

      // Rest back to good state.
      pp.remove("solid0.diameter.min");
      pp.add("solid0.diameter.min", valid_dmin);
      if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }
    }
    Print() << " PASSED!\n";

    // Error for valid distribution type
    Print() << "Normal distribution test 04";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      { std::string dtype = "number weighted";
        pp.remove("solid0.diameter.type");
        pp.add("solid0.diameter.type", dtype);
        if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      }
      { std::string dtype = "numberweighted";
        pp.remove("solid0.diameter.type");
        pp.add("solid0.diameter.type", dtype);
        if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      }

      { std::string dtype = "NuMbEr-WeIgHtEd";
        pp.remove("solid0.diameter.type");
        pp.add("solid0.diameter.type", dtype);
        if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }
      }
    }
    Print() << " PASSED!\n";

    // Error for valid distribution type
    Print() << "Normal distribution test 05";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      { std::string dtype = "volume weighted";
        pp.remove("solid0.diameter.type");
        pp.add("solid0.diameter.type", dtype);
        if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      }
      { std::string dtype = "volumeweighted";
        pp.remove("solid0.diameter.type");
        pp.add("solid0.diameter.type", dtype);
        if (diameter.define(field, "diameter", &dist_data, 0) == PASS) { return FAIL; }
      }

      { std::string dtype = "vOlUmE-wEiGhTeD";
        pp.remove("solid0.diameter.type");
        pp.add("solid0.diameter.type", dtype);
        if (diameter.define(field, "diameter", &dist_data, 0) == FAIL) { return FAIL; }
      }
    }
    Print() << " PASSED!\n";
  }


  Real const Nw_mu    = 1000.0e-6;
  Real const Nw_sigma =  250.0e-6;

  Real const dmin = Nw_mu - 50.*Nw_sigma;
  Real const dmax = Nw_mu + 50.*Nw_sigma;

  Real const Vw_mu = normal_mean_Nw_to_Vw(Nw_mu, Nw_sigma);
  Real const Vw_sigma = normal_stdev_Nw_to_Vw(Nw_mu, Nw_sigma);

  // Check number-weighted distributions
  { std::string field = "normal.solid1";

    Print() << "Normal distribution test 10";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(1);

      pp.add("solid1.diameter.bins", 2048);

      pp.add("solid1.diameter",    normal_dist);
      pp.add("solid1.diameter.mean",     Nw_mu);
      pp.add("solid1.diameter.std",   Nw_sigma);
      pp.add("solid1.diameter.min",       dmin);
      pp.add("solid1.diameter.max",       dmax);
      pp.add("solid1.diameter.type",  Nw_dtype);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      if (!compare(diameter.get_mean(),      Nw_mu, 0.)) { return FAIL; }
      if (!compare(diameter.get_stddev(), Nw_sigma, 0.)) { return FAIL; }
      if (!compare(diameter.get_min(),        dmin, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(),        dmax, 0.)) { return FAIL; }

      if (!diameter.is_number_weighted()) { return FAIL; }
      if ( diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << Nw_mu    << "  " << mean
                << "\nstdd: " << Nw_sigma << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(Nw_mu - mean) > 1.0e-4*Nw_mu) { return FAIL; }
        if (amrex::Math::abs(Nw_sigma - stdev) > 1.0e-4*Nw_sigma) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    // Reprocess the same data, but enable volume weighted sampling
    Print() << "Normal distribution test 11";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      if (!compare(diameter.get_mean(),      Vw_mu, 0.1)) { return FAIL; }
      if (!compare(diameter.get_stddev(), Vw_sigma, 0.1)) { return FAIL; }
      if (!compare(diameter.get_min(),        dmin, 0.0)) { return FAIL; }
      if (!compare(diameter.get_max(),        dmax, 0.0)) { return FAIL; }

      if (!diameter.is_number_weighted()) { return FAIL; }
      if ( diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << Vw_mu    << "  " << mean
                << "\nstdd: " << Vw_sigma << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(Vw_mu - mean) > 1.0e-4*Vw_mu) { return FAIL; }
        if (amrex::Math::abs(Vw_sigma - stdev) > 1.0e-4*Vw_sigma) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";
  } // end number-weighted checks

  // volume-weighted distribution checks
  { std::string field = "normal.solid2";

    Print() << "Normal distribution test 20";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);

      pp.add("solid2.diameter.bins", 2048);

      pp.add("solid2.diameter", normal_dist);
      pp.add("solid2.diameter.mean",    Vw_mu);
      pp.add("solid2.diameter.std",  Vw_sigma);
      pp.add("solid2.diameter.min",      dmin);
      pp.add("solid2.diameter.max",      dmax);
      pp.add("solid2.diameter.type", Vw_dtype);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      if (!compare(diameter.get_mean(),      Vw_mu, 0.)) { return FAIL; }
      if (!compare(diameter.get_stddev(), Vw_sigma, 0.)) { return FAIL; }
      if (!compare(diameter.get_min(),        dmin, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(),        dmax, 0.)) { return FAIL; }

      if ( diameter.is_number_weighted()) { return FAIL; }
      if (!diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << Vw_mu    << "  " << mean
                << "\nstdd: " << Vw_sigma << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(Vw_mu - mean) > 1.0e-4*Vw_mu) { return FAIL; }
        if (amrex::Math::abs(Vw_sigma - stdev) > 1.0e-4*Vw_sigma) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    Print() << "Normal distribution test 21";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(1);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      if (!compare(diameter.get_mean(),      Nw_mu, 5.)) { return FAIL; }
      if (!compare(diameter.get_stddev(), Nw_sigma, 5.)) { return FAIL; }
      if (!compare(diameter.get_min(),        dmin, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(),        dmax, 0.)) { return FAIL; }

      if ( diameter.is_number_weighted()) { return FAIL; }
      if (!diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << Nw_mu    << "  " << mean
                << "\nstdd: " << Nw_sigma << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(Nw_mu - mean) > 1.0e-4*Nw_mu) { return FAIL; }
        if (amrex::Math::abs(Nw_sigma - stdev) > 1.0e-4*Nw_sigma) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";
  } // end volume-weighted checks

  return PASS;
}
