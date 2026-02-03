//
// Average a variable with optional density weighting (Favre average)

#include <regex>
#include <algorithm>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <FilterML.H>

using namespace amrex;

void PrintHelp ()
{
  Print() << "\n"
          << " Compute and report average volume fraction and fluid-solids\n"
          << " relative velocity magnitude and (fine/coarse) heterogeneous\n"
          << " index for drag.\n"
          << "\n"
          << "\n"
          << " Example: \n"
          << " mpirun -n 4 ./FilterML3d.gnu.MPI.ex inputs\n"
          << "\n"
          << "   inputs: text file containing options for post-processing\n"
          << std::endl;
}


int get_plotfile ( std::string& a_pfile );

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

    bool calc_fml_Eulerian = false;
    pp.query("filter.Eulerian", calc_fml_Eulerian);

    Vector<std::string> fml_Eulerian_vars {"alpha_f", "slip_vel", "beta"};

    int filter_size(1024);
    int min_filter_size(1024);

    if (calc_fml_Eulerian) {
      amrex::ParmParse pp_filter("filter");
      pp_filter.get("Eulerian.size", filter_size);
      pp_filter.query("Eulerian.min", min_filter_size);
    }

    std::string pfile;
    get_plotfile(pfile);

    do {

      FilterML filterML(pfile);

      int fml_size = filter_size;

      if (calc_fml_Eulerian) {

        filterML.create_Eulerian_solids(fml_Eulerian_vars);

        do { Print() << "Computing filter size: " << fml_size << '\n';

          filterML.compute_Eulerian_filter(fml_size);

          fml_size /= 2;

        } while (fml_size >= min_filter_size);
      }
    }
    while ( get_plotfile(pfile) );

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

      // Length of the numeric part
      int numLength = a_pfile.size() - first;

      // Name of plot file (including path) without step index
      std::string base = a_pfile.substr(0,first);

      // Find the next plot file that exists
      do {
        std::ostringstream ss;
        ss << base << std::setw(numLength) << std::setfill('0') << ++step;
        a_pfile = ss.str();
      } while (!FileExists(a_pfile) && step <= max_index);

      return ( step <= max_index ? 1 : 0);
    }
  }
  Print() << "Call for an adult. Something when very wrong.\n";
  return 0;
}
