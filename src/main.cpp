#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_buildInfo.H>

#include <mfix.H>

#ifdef MFIX_MPMD
#include <AMReX_MPMD.H>
#endif

#include "build_info.H"
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_rw.H>

using namespace amrex;

// Set defaults that are different that what ARMeX uses.  We only
// add them if they are not already specified in the inputs file.
void add_par () {

  ParmParse pp;

  // Set the extend domain flag by default
  bool extend_face(true);
  pp.queryAdd("eb2.extend_domain_face", extend_face);

  // Disable managed memory for GPUs
  int not_managed(0);
  pp.queryAdd("amrex.the_arena_is_managed", not_managed);

  // Disable GPU aware MPI by default (Joule3 workaround)
  bool not_gpu_aware(false);
  pp.queryAdd("amrex.use_gpu_aware_mpi", not_gpu_aware);

  // Set default output to NATIVE
  std::string format = "NATIVE";
  pp.queryAdd("fab.format", format);
}


const char* HypreVersion ();
void writeBuildInfo ();

int main (int argc, char* argv[])
{
    // check to see if it contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                writeBuildInfo();
                return 0;
            }
        }
    }

    // Issue an error if AMR input file is not given
    if ( argc < 2 ) {
       std::cerr << "AMReX input file missing" << std::endl << std::endl;
       std::cerr << "Usage:  " << argv[0] << " inputs [--describe]" << std::endl;
       return -1;
    }

    // AMReX will now read the inputs file and the command line arguments, but the
    //        command line arguments are in mfix-format so it will just ignore them.
#ifdef MFIX_MPMD
    MPI_Comm app_comm = amrex::MPMD::Initialize(argc,argv);
    amrex::Initialize(argc, argv, true, app_comm, add_par);
#else
    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, add_par);
#endif

    { // This start bracket and the end bracket before Finalize are essential so
      // that the mfix object is deleted before Finalize
    BL_PROFILE_VAR("main()", pmain)
    BL_PROFILE_REGION_START("mfix::main()");

    // Write out the MFIX git hash (the AMReX git hash is already written)
    const char* githash_mfix = buildInfoGetGitHash(1);
    Print() << "   MFIX git describe: " << githash_mfix<< "\n";
    Print() << "AMReX-Hydro git hash: " << HydroGitHash() << "\n";
#ifdef CSG_EB
    Print() << "     CSG-EB git hash: " << CsgEbGitHash() << "\n";
#endif
#ifdef AMREX_USE_HYPRE
    Print() << "       HYPRE Version: " << HypreVersion() << "\n";
#endif

    // Default constructor. Note inheritance: mfix : AmrCore : AmrMesh
    //                                                             |
    //  => Geometry is constructed here: (constructs Geometry) ----+
    mfix mfix;

    auto& rw = *(mfix.m_rw);
    auto& timer = mfix.timer();
    auto& eb = *(mfix.m_eb);

    // Initialize internals from ParamParse database
    mfix.InitParams();

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    eb.make_geometry(mfix.Geom(), rw.restart_file);

    // Initialize derived internals
    mfix.Init(timer.time());

    Real init_time = timer.elapsed_runtime();
    ParallelDescriptor::ReduceRealMax(init_time, ParallelDescriptor::IOProcessorNumber());

    Print() << "Time spent in init      " << init_time << "\n\n";

    if (!rw.only_print_grid_report) {

       if ( ParmParse::QueryUnusedInputs() && mfix.run_type(RunType::Standard) ) {
         if (rw.stop_for_unused_inputs ) {
           reporter::Log(reporter::Error,__FILE__, __LINE__)
               << "Aborting due to unused inputs.";
         }
       }

       rw.writeNow(timer, /*first=*/true, /*last=*/false);
       rw.writeStaticPlotFiles();

       mfix.m_mass_balance->ComputeMassAccum(rw.report_mass_balance, 0);

       { // Start profiling solve here

         BL_PROFILE("mfix_solve");

         if (timer.SteadyState()) {

           mfix.EvolveSteadyState();
           timer.nstep() = 1;

         } else {

           while (timer.ok()) {

             timer.step_start();

             mfix.Regrid( timer.nstep() );

             mfix.Evolve();

             Print() << "   Time per step        " << timer.step_end() << '\n';

             timer.advance_time();
             timer.advance_nstep();

             rw.writeNow(timer);
           }
         }
       }

       if (timer.run_status_type() != MFIXRunStatusType::RuntimeIsOver ||
           timer.run_status_type() == MFIXRunStatusType::UserStop)
       { rw.writeNow(timer, /*first=*/false, /*last=*/true); }


       Real end_time = timer.elapsed_runtime();
       ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

       Print() << "Time spent in main (after init) " << end_time-init_time << '\n'
               << "Time spent in main      " << end_time << "\n\n";

       BL_PROFILE_REGION_STOP("mfix::main()");
       BL_PROFILE_VAR_STOP(pmain);

    }

    } // This end bracket and the start bracket after Initialize are essential so
      // that the mfix object is deleted before Finalize

    amrex::Finalize();
#ifdef MFIX_MPMD
    amrex::MPMD::Finalize();
#endif
    return 0;
}
