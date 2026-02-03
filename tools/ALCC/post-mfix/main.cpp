//
// Average a variable with optional density weighting (Favre average)

#include <regex>
#include <algorithm>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <post_mfix.H>

using namespace amrex;

void PrintHelp ()
{
  Print() << "\n"
          << " Compute fluid avarges for variable in a plotfile\n"
          << " usage: \n"
          << "    post-mfix inputs\n"
          << "\n"
          << "   inputs         : file containing options for post-processing\n"
          << "\n"
          << "\n"
          << " Example: \n"
          << " ./post-mfix3d.gnu.ex inputs\n"
          << std::endl;
}


int get_plotfile ( std::string& a_pfile );

void write_averages ( Vector<std::string> a_vars,
                      Vector<int> a_steps, Vector<Real> a_times,
                      Vector<Vector< Vector<Real>>> a_avgs,
                      Vector<Vector< Vector<Real>>> a_flcts);

int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, true);

  {
    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    if( amrex::command_argument_count() == 0) {
      PrintHelp();
      return  0;
    }

    amrex::ParmParse pp;

    bool write_hdf5 = false;
    pp.query("hdf5", write_hdf5);
#ifndef AMREX_USE_HDF5
    if( write_hdf5) {
      amrex::Print() << "HDF5 output enabled but post-mfix was not built with HDF5 support!\n";
      return  -1;
    }
#endif

    Vector<std::string> vars;
    pp.queryarr("averages",vars);

    bool calc_avgs = false;
    if (vars.size() > 0) { calc_avgs = amrex::toLower(vars[0]).compare("none"); }

    std::string pfile;
    get_plotfile(pfile);

    Vector<int> steps;
    Vector<amrex::Real> times;
    Vector<Vector< Vector<Real>>> avgs, flcts;
    do {

      pm_plotfile plotfile(pfile);

      steps.push_back(plotfile.get_step());
      times.push_back(plotfile.get_time());

      post_mfix post(plotfile.get_level_0_geom(),
                     plotfile.get_amr_info(),
                     &plotfile);

      if (calc_avgs) {
        avgs.push_back( post.compute_averages(vars) );
        flcts.push_back( post.compute_fluctuations(vars, avgs.back()) );
      }

      if (write_hdf5) { post.output_hdf5(pfile); }

    }
    while ( get_plotfile(pfile) );

    if (calc_avgs) { write_averages( vars, steps, times, avgs, flcts ); }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

  }
  amrex::Finalize();
}

int get_plotfile ( std::string& a_pfile ) {

  amrex::ParmParse pp;
  if (a_pfile.size() == 0) {
    pp.get("plotfile", a_pfile);
    //Print() << "This is our first plot file. Done.\n";
    return 1;
  } else {

    int max_index = -1;
    if ( !pp.query("plotfile.max_index", max_index) ) {
      // We are only processing one plot file and we're done.
      return 0;
    } else {
      // Find the start of the plotfile time (e.g., plt)
      size_t last = a_pfile.find_last_of("/\\");
      if (last == std::string::npos) { last = 0; }

      // Get the index of the current plot file.
      int const first = a_pfile.find_first_of("0123456789",last);
      int step = std::stoi( a_pfile.substr(first) );

      // Name of plot file (including path) without step index
      std::string base = a_pfile.substr(0,first);

      // Find the next plot file that exists
      do { a_pfile = amrex::Concatenate(base,++step);
      } while (!FileExists(a_pfile) && step <= max_index);

      return ( step <= max_index ? 1 : 0);
    }
  }
  Print() << "Call for an adult. Something when very wrong.\n";
  return 0;
}


void write_averages ( Vector<std::string> a_vars,
                      Vector<int> a_steps, Vector<Real> a_times,
                      Vector<Vector< Vector<Real>>> a_avgs,
                      Vector<Vector< Vector<Real>>> a_flcts){

  int const sample_size = static_cast<int>(a_steps.size());

  int const nlev = static_cast<int>(a_avgs[0].size());

  int const ncomp = static_cast<int>(a_vars.size());

  if (sample_size == 1 ) {

    for (int lev(0); lev<nlev; ++lev) {

      Print() << "\nResults on level " << lev
              << " at step " << a_steps[0]
              << " and time " << a_times[0] << ":\n";
      Print() << "========================================================================\n"
              << "                       variable         average             fluctuations\n";

      for (int comp(0); comp < ncomp; ++comp) {

        Print() << std::setw(32) << std::setfill(' ') << std::right << "<"+a_vars[comp]+">"
          << std::setw(18) << std::setprecision(6) << std::scientific
                             << std::right << a_avgs[0][lev][comp]
          << "    "
          << std::setw(18) << std::setprecision(6) << std::scientific
                             << std::right << a_flcts[0][lev][comp]
          << "\n";
      }
      Print() << "========================================================================\n\n\n";
    }
  } else {

    int const fill_size = 2 + 32 + 18*sample_size;

    for (int lev(0); lev<nlev; ++lev) {

      Print() << "\nAverages on level " << lev << ":\n";
      Print() << std::setw(fill_size) << std::setfill('=') << "\n";

      Print() << std::setw(32) << std::setfill(' ') << std::right << "step";
      for (int sample(0); sample<sample_size; ++sample) {
          Print() << std::setw(18) << std::right << a_steps[sample];
      }
      Print() << "\n";

      Print() << std::setw(32) << std::setfill(' ') << std::right << "time";
      for (int sample(0); sample<sample_size; ++sample) {
          Print() << std::setw(18) << std::right << a_times[sample];
      }
      Print() << "\n\n";

      for (int comp(0); comp < ncomp; ++comp) {

        Print() << std::setw(32) << std::setfill(' ') << std::right << "<"+a_vars[comp]+">";

        for (int sample(0); sample<sample_size; ++sample) {
          Print() << std::setw(18) << std::setprecision(6) << std::scientific
                             << std::right << a_avgs[sample][lev][comp];
        }
        Print() << "\n";
      }
      Print() << std::setw(fill_size) << std::setfill('=') << "\n";
      Print() << "\n\n";

      Print() << "\nFluctuations on level " << lev << ":\n";
      Print() << std::setw(fill_size) << std::setfill('=') << "\n";

      Print() << std::setw(32) << std::setfill(' ') << std::right << "step";
      for (int sample(0); sample<sample_size; ++sample) {
          Print() << std::setw(18) << std::right << a_steps[sample];
      }
      Print() << "\n";

      Print() << std::setw(32) << std::setfill(' ') << std::right << "time";
      for (int sample(0); sample<sample_size; ++sample) {
          Print() << std::setw(18) << std::right << a_times[sample];
      }
      Print() << "\n\n";

      for (int comp(0); comp < ncomp; ++comp) {

        Print() << std::setw(32) << std::setfill(' ') << std::right << "<"+a_vars[comp]+">";

        for (int sample(0); sample<sample_size; ++sample) {
          Print() << std::setw(18) << std::setprecision(6) << std::scientific
                             << std::right << a_flcts[sample][lev][comp];
        }
        Print() << "\n";
      }
      Print() << std::setw(fill_size) << std::setfill('=') << "\n";
      Print() << "\n\n";

    }
  }



};
