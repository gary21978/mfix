#include <test_dist.H>

using namespace amrex;

int test_dist::
compute_mean_and_std ( INPUT_DIST_t & a_diameter,
                       unsigned long const a_samples,
                       Real& a_mean,
                       Real& a_stddev,
                       int const a_print_pdf,
                       int const a_bins,
                       int const a_use_midpoint)
{
  amrex::RandomEngine engine;
  amrex::ResetRandomSeed(0);

  int const bins = (a_bins > 0) ? a_bins : a_diameter.get_bins();
  if (bins < 1) { return FAIL; }

  amrex::Vector<int> counts(bins, 0);

  Real const dmin = a_diameter.get_min();
  Real const dmax = a_diameter.get_max();

  amrex::Real db(dmax - dmin);
  if (bins > 1){ db = (dmax-dmin)/static_cast<amrex::Real>(bins); }
  else { db = 1.0; }//amrex::max(db, std::numeric_limits<amrex::Real>::epsilon()); }

  amrex::Real idb = 1.0/db;

  amrex::Real dlo = dmin - ((bins > 1) ? 0.0 : std::numeric_limits<amrex::Real>::epsilon() );
  amrex::Real dhi = dmax + ((bins > 1) ? 0.0 : std::numeric_limits<amrex::Real>::epsilon() );

  Vector<Real> custom_bins;
  if ( a_diameter.is_custom() ) {
    custom_bins.resize(a_diameter.get_bins());
    for (int bin(0); bin<bins; ++bin) {
      auto const cbin = a_diameter.get_custom_bin(bin);
      custom_bins[bin] = std::get<0>(cbin);
    }
  }

  int const a_print_samples(0);

  if (a_print_samples) { amrex::Print() <<  "\n\nGenerate samples:\n"; }
  for (unsigned long i(0); i<a_samples; ++i) {
    Real const dp  = a_diameter.sample<amrex::RunOn::Host>(engine);
    if ( dp > dhi ) { return FAIL; }
    if ( dp < dlo ) { return FAIL; }
    if (a_print_samples) { amrex::Print() << "  " << i << "  " << dp << "\n"; }

    int bin;
    if ( !a_diameter.is_custom() ) {
      bin = static_cast<int>((dp - dlo)*idb);
    } else {
      auto ub = std::upper_bound(custom_bins.begin(), custom_bins.end(), dp);
      if (ub == custom_bins.end() ) { bin = bins-1; }
      else { bin = std::distance(custom_bins.begin(), ub)-1; }
    }
    if (bin < 0 || bin > bins-1) { return FAIL; }
    counts[bin]++;
  }
  if (a_print_samples) { Print() << "\n\n"; }

  a_mean = 0.;
  a_stddev = 0.;

  if (a_print_pdf) { Print() << "\n   bin       PDF\n"; }

  Real const midpoint ((a_use_midpoint && (bins > 1)) ? 0.5 : 0.0);

  int const use_uniform_bins = (!a_diameter.is_custom() ||
                                 a_diameter.get_bins() != bins);

  // Compute the mean and standard deviation
  for (int bin(0); bin<bins; ++bin) {

    Real dp;
    if (use_uniform_bins) {
      dp = dmin + db*(static_cast<Real>(bin) + midpoint);
    } else {
      if (a_use_midpoint && (bin+1 < bins)) {
        dp = 0.5*(custom_bins[bin] + custom_bins[bin+1]);
      } else {
        dp = custom_bins[bin];
      }
    }

    Real const p = static_cast<Real>(counts[bin])/static_cast<Real>(a_samples);

    a_mean += dp*p;
    a_stddev += dp*dp*p;

    if (a_print_pdf) {
      Print() << std::setw(6) << bin << "  "
              << std::setw(8) << counts[bin] << "  "
        << std::scientific << std::setw(8) << std::setprecision(2) << dp << "   "
        << std::scientific << std::setw(12) << std::setprecision(6) << p << '\n';
    }
  }
  if (a_print_pdf) { Print() << "\n\n"; }

  a_stddev = std::sqrt(a_stddev - a_mean*a_mean);

  return 0;
}

