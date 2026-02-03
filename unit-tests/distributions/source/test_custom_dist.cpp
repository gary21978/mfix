#include <fstream>
#include <iostream>

#include <AMReX_ParallelDescriptor.H>

#include <test_custom_dist.H>

using namespace amrex;

int test_custom_dist::test()
{
  ParmParse pp("custom");

  std::string custom = "custom";

  Print() << "\nCustom distribution tests:\n";

  Real const Nw_mu    = 1000.0e-6;
  Real const Nw_sigma =  250.0e-6;

  Real const Nw_min = Nw_mu - 3.*Nw_sigma;
  Real const Nw_max = Nw_mu + 3.*Nw_sigma;

  Real const Vw_mu(normal_mean_Nw_to_Vw(Nw_mu, Nw_sigma));
  Real const Vw_sigma(normal_stdev_Nw_to_Vw(Nw_mu, Nw_sigma));

  Real const Vw_min = Vw_mu - 3.*Vw_sigma;
  Real const Vw_max = Vw_mu + 3.*Vw_sigma;

  // Check for an invalid inputs
  { std::string field = "custom.solid0";

    int const dist_samples = 32;

    // ---> number-weighted sampling of number-weighted distribution
    // ---> do NOT rebase CDF
    // ---> do NOT interpolate between bins
    Print() << "Custom distribution test 01";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      // Create a CDF by taking dist_samples from a normal distribution
      // with mean Nw_mu and standard deivation Nw_sigma.
      std::string filename = "custom_cdf_01.txt";

      int const is_cdf(1);
      int const rebase(0);
      create_custom_normal(filename, is_cdf, dist_samples, rebase,
        Nw_min, Nw_max, Nw_mu, Nw_sigma);

      int const use_number_weighted(1);          // number-weighted sampling

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == PASS) { return FAIL; }
      pp.add("solid0.diameter", custom);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == PASS) { return FAIL; }
      pp.add("solid0.diameter.type",   Nw_dtype);

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == PASS) { return FAIL; }
      pp.add("solid0.diameter.custom", filename);

      pp.add("solid0.diameter.interpolate", 1);  // turn on  interpolation <-- should fail

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == PASS) { return FAIL; }

      pp.remove("solid0.diameter.interpolate");
      pp.add("solid0.diameter.interpolate", 0);  // turn off interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      if (!compare(diameter.get_min(), Nw_min, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), Nw_max, 0.)) { return FAIL; }

      if (!diameter.is_number_weighted()) { return FAIL; }
      if ( diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(0);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean - dist_mean) > 5.0e-4*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 1.0e-3*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    // Reprocess the same data, but
    // ---> number-weighted sampling of number-weighted distribution
    // ===> rebase CDF
    // ---> do NOT interpolate between bins
    Print() << "Custom distribution test 02";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      // Create a CDF by taking dist_samples from a normal distribution
      // with mean Nw_mu and standard deivation Nw_sigma.
      std::string filename = "custom_cdf_02.txt";

      int const is_cdf(1);
      int const rebase(1);
      create_custom_normal(filename, is_cdf, dist_samples, rebase,
        Nw_min, Nw_max, Nw_mu, Nw_sigma);

      pp.remove("solid0.diameter.custom");
      pp.add("solid0.diameter.custom", filename);

      int const use_number_weighted(1);          // number-weighted sampling

      pp.remove("solid0.diameter.interpolate");
      pp.add("solid0.diameter.interpolate", 0);  // turn off interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      if (!compare(diameter.get_min(), Nw_min, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), Nw_max, 0.)) { return FAIL; }

      if (!diameter.is_number_weighted()) { return FAIL; }
      if ( diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(0);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean - dist_mean) > 5.0e-4*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 1.0e-3*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    // Reprocess the same data, but
    // ---> number-weighted sampling of number-weighted distribution
    // ===> rebase CDF
    // ===> interpolate between bins
    Print() << "Custom distribution test 03";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(1);          // number-weighted sampling

      pp.remove("solid0.diameter.interpolate");
      pp.add("solid0.diameter.interpolate", 1);  // turn on  interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      if (!compare(diameter.get_min(), Nw_min, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(), Nw_max, 0.)) { return FAIL; }

      if (!diameter.is_number_weighted()) { return FAIL; }
      if ( diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        // Increase the number of bins to use bin midpoint and compare with
        // the mean and sigma of base distribution.
        Real dist_mean = Nw_mu;
        Real dist_stdev = Nw_sigma;

        int const print_pdf(0);
        int const use_midpoint(1);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
        { return FAIL; }
#else
        if (amrex::Math::abs(mean - dist_mean) > 5.0e-3*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 5.0e-2*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

    // Reprocess the same data, but
    // ===> volume-weighted sampling of number-weighted distribution
    // ===> rebase CDF
    // ---> do NOT interpolate between bins
    Print() << "Custom distribution test 04";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);          // use volume weighted

      pp.remove("solid0.diameter.interpolate");
      pp.add("solid0.diameter.interpolate", 0);  // turn off interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      // Error comparing to analytical values
      if (std::abs(diameter.get_mean() - Vw_mu)      > 2.5e-2*Vw_mu   ) { return FAIL; }
      if (std::abs(diameter.get_stddev() - Vw_sigma) > 2.5e-2*Vw_sigma) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(0);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << mean  << "  " << dist_mean
                << "\nstdd: " << stdev << "  " << dist_stdev
                << "\n";
#else
        if (amrex::Math::abs(mean - dist_mean) > 5.0e-4*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 5.0e-4*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";


    // Reprocess the same data, but
    // ===> volume-weighted sampling of number-weighted distribution
    // ===> rebase CDF
    // ===> do NOT interpolate between bins
    Print() << "Custom distribution test 05";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);          // use volume weighted

      pp.remove("solid0.diameter.interpolate");
      pp.add("solid0.diameter.interpolate", 1);  // turn on  interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(1);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean - dist_mean) > 2.5e-2*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 2.5e-2*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

  } // custom.solid0


  { std::string field = "custom.solid1";

    int const dist_samples = 27;


    // ---> volume-weighted sampling of volume-weighted distribution
    // ---> do NOT rebase PDF
    // ---> do NOT interpolate between bins
    Print() << "Custom distribution test 11";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      // Create a PDF by taking dist_samples from a normal distribution
      // with mean Vw_mu and standard deivation Vw_sigma.
      std::string filename = "custom_pdf_11.txt";

      int const is_cdf(0);
      int const rebase(0);
      create_custom_normal(filename, is_cdf, dist_samples, rebase,
        Vw_min, Vw_max, Vw_mu, Vw_sigma);

      int const use_number_weighted(0);          // volume-weighted sampling

      pp.add("solid1.diameter", custom);
      pp.add("solid1.diameter.type",   Vw_dtype);
      pp.add("solid1.diameter.custom", filename);

      pp.add("solid1.diameter.interpolate", 1);  // turn on  interpolation <-- should fail

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == PASS) { return FAIL; }

      pp.remove("solid1.diameter.interpolate");
      pp.add("solid1.diameter.interpolate", 0);  // turn off interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      // Error comparing to analytical values
      if (std::abs(diameter.get_mean() - Vw_mu)      > 1.0e-6*Vw_mu   ) { return FAIL; }
      if (std::abs(diameter.get_stddev() - Vw_sigma) > 1.0e-2*Vw_sigma) { return FAIL; }

      if ( diameter.is_number_weighted()) { return FAIL; }
      if (!diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(0);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << mean  << "  " << dist_mean  << "  diff: " << std::abs(mean-dist_mean)
                << "\nstdd: " << stdev << "  " << dist_stdev << "  diff: " << std::abs(stdev-dist_stdev)
                << "\n";
#else
        if (amrex::Math::abs(mean-dist_mean) > 2.5e-2*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 5.0e-2*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";


    // ---> volume-weighted sampling of volume-weighted distribution
    // ===> rebase PDF
    // ---> do NOT interpolate between bins
    Print() << "Custom distribution test 12";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      // Create a PDF by taking dist_samples from a normal distribution
      // with mean Vw_mu and standard deivation Vw_sigma.
      std::string filename = "custom_pdf_12.txt";

      int const is_cdf(0);
      int const rebase(1);
      create_custom_normal(filename, is_cdf, dist_samples, rebase,
        Vw_min, Vw_max, Vw_mu, Vw_sigma);

      int const use_number_weighted(0);          // volume-weighted sampling

      pp.add("solid1.diameter.custom", filename);
      pp.add("solid1.diameter.interpolate", 0);  // turn off interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      // Error comparing to analytical values
      if (std::abs(diameter.get_mean() - Vw_mu)      > 1.0e-6*Vw_mu   ) { return FAIL; }
      if (std::abs(diameter.get_stddev() - Vw_sigma) > 5.0e-2*Vw_sigma) { return FAIL; }

      if ( diameter.is_number_weighted()) { return FAIL; }
      if (!diameter.is_volume_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(0);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean-dist_mean) > 5.0e-3*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 5.0e-3*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";


    // ---> volume-weighted sampling of volume-weighted distribution
    // ===> rebase PDF
    // ===> interpolate between bins
    Print() << "Custom distribution test 13";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      // Create a PDF by taking dist_samples from a normal distribution
      // with mean Vw_mu and standard deivation Vw_sigma.

      int const use_number_weighted(0);          // volume-weighted sampling

      pp.remove("solid1.diameter.interpolate");
      pp.add("solid1.diameter.interpolate", 1);  // turn on  interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(1);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean-dist_mean) > 5.0e-2*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 5.0e-2*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";


    // ===> number-weighted sampling of volume-weighted distribution
    // ===> rebase PDF
    // ---> do NOT interpolate between bins
    Print() << "Custom distribution test 14";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      // Create a PDF by taking dist_samples from a normal distribution
      // with mean Vw_mu and standard deivation Vw_sigma.

      int const use_number_weighted(1);          // number-weighted sampling

      pp.remove("solid1.diameter.interpolate");
      pp.add("solid1.diameter.interpolate", 0);  // turn off interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(0);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean-dist_mean) > 5.0e-4*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 5.0e-5*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";


    // ===> number-weighted sampling of volume-weighted distribution
    // ===> rebase PDF
    // ===> interpolate between bins
    Print() << "Custom distribution test 15";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(1);          // number-weighted sampling

      pp.remove("solid1.diameter.interpolate");
      pp.add("solid1.diameter.interpolate", 1);  // turn on  interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(1);
        if (compute_mean_and_std(diameter, samples, mean, stdev,
              print_pdf, diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean-dist_mean) > 5.0e-2*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 5.0e-2*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

  } // custom.solid1


  { std::string field = "custom.solid2";

    // Create a PDF of bi-disperse mixture -- v0
    // ---> volume-weighted sampling of volume-weighted distribution
    // ---> do NOT interpolate between bins
    std::string filename = "custom_pdf_21.txt";

    int const is_cdf(0);

    int const bins(3);
    amrex::Vector<amrex::Real> bin_dp( {0.5, 1.0, 2.0} );
    amrex::Vector<amrex::Real> dist(   {0.0, 0.5, 0.5} );

    create_custom(filename, is_cdf, bins, bin_dp.data(), dist.data());

    // do NOT interpolate between bins
    Print() << "Custom distribution test 21";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);          // volume-weighted sampling

      pp.add("solid2.diameter", custom);
      pp.add("solid2.diameter.type",   Vw_dtype);
      pp.add("solid2.diameter.custom", filename);

      pp.remove("solid2.diameter.interpolate");
      pp.add("solid2.diameter.interpolate", 0);  // turn off interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      if (!compare(diameter.get_mean(),   1.5, 0.)) { return FAIL; }
      if (!compare(diameter.get_stddev(), 0.5, 0.)) { return FAIL; }
      if (!compare(diameter.get_min(),    0.5, 0.)) { return FAIL; }
      if (!compare(diameter.get_max(),    2.0, 0.)) { return FAIL; }

      if (diameter.is_number_weighted()) { return FAIL; }

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);
        Real dist_mean = diameter.get_mean();
        Real dist_stdev = diameter.get_stddev();

        int const print_pdf(0);
        int const use_midpoint(0);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean-dist_mean) > 5.0e-4*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 5.0e-6*dist_stdev) { return FAIL; }
#endif
      }

    }
    Print() << " PASSED!\n";

    // Reprocess the same data, but
    // ---> volume-weighted sampling of volume-weighted distribution
    // ===> interpolate between bins
    Print() << "Custom distribution test 22";
    { INPUT_DIST_t diameter;
      INPUT_DIST_DATA_t dist_data;

      int const use_number_weighted(0);          // volume-weighted sampling

      pp.remove("solid2.diameter.interpolate");
      pp.add("solid2.diameter.interpolate", 1);  // turn on  interpolation

      if (diameter.define(field, "diameter", &dist_data, use_number_weighted) == FAIL) { return FAIL; }
      diameter.copyAsync( &dist_data );

      { unsigned long const samples( 500000 );
        Real mean(0.), stdev(0.);

        // Analytically computed values
        Real dist_mean(1.125);
        Real dist_stdev(0.375);

        int const print_pdf(0);
        int const use_midpoint(1);

        if (compute_mean_and_std(diameter, samples, mean, stdev, print_pdf,
              diameter.get_bins(), use_midpoint) == FAIL) { return FAIL; }
#if 0
        Print() << "\nmean: " << dist_mean  << "  " << mean
                << "\nstdd: " << dist_stdev << "  " << stdev
                << "\n";
#else
        if (amrex::Math::abs(mean-dist_mean) > 5.0e-4*dist_mean) { return FAIL; }
        if (amrex::Math::abs(stdev - dist_stdev) > 1.0e-6*dist_stdev) { return FAIL; }
#endif
      }
    }
    Print() << " PASSED!\n";

  } // custom.solid2
  return PASS;


}


void
test_custom_dist::
create_custom_normal ( std::string const a_name, int const a_is_cdf,
                       int const a_bins, int const a_rebase,
                       amrex::Real a_min, amrex::Real a_max,
                       amrex::Real a_mean, amrex::Real a_stdev)
{
  Vector<Real> bin_dp(a_bins, 0.), dist(a_bins, 0.);

  Real const db((a_max - a_min)/static_cast<Real>(a_bins-1));

  Real const inv_sqrt2( 1.0/std::sqrt(2.));
  Real const inv_stdev( 1.0/a_stdev );
  Real const inv_stdev_sqrt2pi( 1.0/(a_stdev*std::sqrt(2.*M_PI)) );

  Real shift(0.);
  Real scale(1.);

  if (a_rebase) {
    Real x( (a_min - a_mean)*inv_stdev );
    shift = (a_is_cdf) ? 0.5*(1.0 + std::erf( x*inv_sqrt2 ))
                       : std::exp(-0.5*x*x)*inv_stdev_sqrt2pi;
    if (a_is_cdf) {
      x = (a_max - a_mean)*inv_stdev;
      scale = 1.0/(0.5*(1.0 + std::erf( x*inv_sqrt2 )) - shift);
    }
  }

  for (int n(0); n<a_bins; ++n) {

    bin_dp[n] = a_min + db*static_cast<amrex::Real>(n);

    amrex::Real const x( (bin_dp[n] - a_mean)*inv_stdev );

    dist[n] = (a_is_cdf) ? 0.5*(1.0 + std::erf( x*inv_sqrt2 ))
                         : std::exp(-0.5*x*x)*inv_stdev_sqrt2pi;

    dist[n] -= shift;
    dist[n] *= scale;
  }

  create_custom(a_name, a_is_cdf, a_bins, bin_dp.data(), dist.data());
}


void
test_custom_dist::
create_custom ( std::string const a_name, int const is_cdf, int const a_entries,
                Real const* a_bins, Real const* a_cdf ) {

  if (!ParallelDescriptor::IOProcessor()) { return; }

  std::remove(a_name.c_str());

  std::ofstream cfile;
  cfile.open(a_name.c_str());

  cfile << a_entries << "," << 1 << '\n';
  cfile << (is_cdf ? "CDF" : "PDF") << '\n';
  cfile << '\n';

  for (int entry(0); entry< a_entries; ++entry) {
    cfile << a_bins[entry] << "   " << a_cdf[entry] << '\n';
  }
}
