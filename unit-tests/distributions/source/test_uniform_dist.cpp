#include <test_uniform_dist.H>

using namespace amrex;

int test_uniform_dist::test () {

  ParmParse pp("uniform");

  std::string uniform_dist = "uniform";

  Print() << "\nUniform distribution tests:\n";

  // Checks on the validity of input values
  { std::string field = "uniform.solid0";

    pp.add("solid0.diameter", uniform_dist);

    Real const valid_dmin =  250.0e-6;
    Real const valid_dmax = 2500.0e-6;

    // Base inputs are valid -- make sure that holds.
    Print() << "Uniform distribution test 01";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

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
      if ( diameter.is_log_normal() ) { return FAIL; }
      if (!diameter.is_uniform()    ) { return FAIL; }
      if ( diameter.is_custom()     ) { return FAIL; }
    }
    Print() << " PASSED!\n";

    // Flag error if the number of bins is invalid
    Print() << "Uniform distribution test 02";
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
    Print() << "Uniform distribution test 03";
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
  }

  Real const dmin =  250.e-6;
  Real const dmax = 2500.e-6;

  // number-weighted checks
  { std::string field = "uniform.solid1";

    Print() << "Uniform distribution test 10";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(1);

      pp.add("solid1.diameter.bins", 2048);

      pp.add("solid1.diameter",  uniform_dist);
      pp.add("solid1.diameter.min",      dmin);
      pp.add("solid1.diameter.max",      dmax);
      pp.add("solid1.diameter.type", Nw_dtype);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      Real const m( 0.5*(dmax + dmin));
      Real const s( (dmax - dmin)/ sqrt(12.) );

      if (!compare(diameter.get_mean(),   m, 0.)) { return FAIL; }
      if (!compare(diameter.get_stddev(), s, 0.)) { return FAIL; }
      if (!compare(diameter.get_min(), dmin, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), dmax, 0.)) { return FAIL; }

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
        if (amrex::Math::abs(m - mean) > 1.0e-3*m) { return FAIL; }
        if (amrex::Math::abs(s - stdev) > 1.0e-4*s) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    // Reprocess the previous inputs with volume weighted sampling enabled.
    Print() << "Uniform distribution test 11";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      Real m, s;

      { Real const a(dmin),     b(dmax);

        Real const a4(a*a*a*a), b4(b*b*b*b);
        Real const a5(a*a4),    b5(b*b4);
        Real const a6(a*a5),    b6(b*b5);

        m = (4.0*(b5 - a5)) / (5.0*(b4 - a4));

        Real const v(4.0*((b6-a6)/6. - 2.*m*(b5-a5)/5. + m*m*(b4-a4)/4. ) / (b4-a4));

        s = std::sqrt(v);
      }

      if (!compare(diameter.get_mean(),   m, .1)) { return FAIL; }
      if (!compare(diameter.get_stddev(), s, .1)) { return FAIL; }
      if (!compare(diameter.get_min(), dmin, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), dmax, 0.)) { return FAIL; }

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
        if (amrex::Math::abs(m - mean) > 1.0e-3*m) { return FAIL; }
        if (amrex::Math::abs(s - stdev) > 5.0e-3*s) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";
  } // end number-weighted checks


  // volume-weighted checks
  { std::string field = "uniform.solid2";

    Print() << "Uniform distribution test 20";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);

      pp.add("solid2.diameter.bins", 2048);

      pp.add("solid2.diameter",  uniform_dist);
      pp.add("solid2.diameter.min",      dmin);
      pp.add("solid2.diameter.max",      dmax);
      pp.add("solid2.diameter.type", Vw_dtype);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      Real const m( 0.5*(dmax + dmin));
      Real const s( (dmax - dmin)/ sqrt(12.) );

      if (!compare(diameter.get_mean(),   m, 0.)) { return FAIL; }
      if (!compare(diameter.get_stddev(), s, 0.)) { return FAIL; }
      if (!compare(diameter.get_min(), dmin, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), dmax, 0.)) { return FAIL; }

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
        if (amrex::Math::abs(m - mean) > 1.0e-3*m) { return FAIL; }
        if (amrex::Math::abs(s - stdev) > 1.0e-4*s) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    // Reprocess the previous inputs with number weighted sampling enabled.
    Print() << "Uniform distribution test 21";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(1);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }

      Real m, s;
      { Real const a(dmin), b(dmax);
        Real const a2(a*a), b2(b*b);

        m = (2.*(a*b)) / (b + a);

        Real const v(2.*(a2*b2*std::log(b/a) - 2.*m*a*b*(b-a)) / (b2-a2) + m*m);

        s = std::sqrt(v);
      }

      if (!compare(diameter.get_mean(),   m, .1)) { return FAIL; }
      if (!compare(diameter.get_stddev(), s, .1)) { return FAIL; }
      if (!compare(diameter.get_min(), dmin, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), dmax, 0.)) { return FAIL; }

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
        if (amrex::Math::abs(m - mean) > 1.0e-3*m) { return FAIL; }
        if (amrex::Math::abs(s - stdev) > 5.0e-3*s) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";
  } // volume-weighted

  return PASS;
}
