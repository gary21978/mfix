#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_buildInfo.H>

#include <mfix.H>
#include <restarter.H>

#include <mfix_dem.H>
#include <mfix_fluid.H>
#include <mfix_rw.H>


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

  // Set default output to NATIVE
  std::string format = "NATIVE";
  pp.queryAdd("fab.format", format);
}


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
  amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, add_par);

  {
    // Write out the MFIX git hash (the AMReX git hash is already written)
    const char* githash_mfix = buildInfoGetGitHash(1);
    amrex::Print() << "   MFIX git describe: " << githash_mfix << "\n";

    // Setting format to NATIVE rather than default of NATIVE_32
    FArrayBox::setFormat(FABio::FAB_NATIVE);


    // ************************************************************************
    // Create the coarse-mesh mfix object for reading the PIC chckpt file
    // ************************************************************************

    mfix* mfix_coarse = new mfix;

    auto& rw_coarse = *(mfix_coarse->m_rw);
    auto& timer_coarse = mfix_coarse->timer();
    auto& eb_coarse = *(mfix_coarse->m_eb);

    // Initialize internals from ParamParse database
    mfix_coarse->InitParams();

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    eb_coarse.make_geometry(mfix_coarse->Geom(),
                            rw_coarse.restart_file);

    // Initialize derived internals
    mfix_coarse->set_run_type(RunType::Restart);
    mfix_coarse->Init(timer_coarse.time());


    // Fill level-sets on each level
    if (mfix_coarse->has_particles()) {
      eb_coarse.fill_levelsets(mfix_coarse->Geom(),
                               mfix_coarse->pc,
                               mfix_coarse->m_boundary_conditions,
                               mfix_coarse->m_porous_media);
    }

    // ************************************************************************
    // Create the PIC2DEM restarter object
    // ************************************************************************

    MFIXRestarter::change_inputs_table();

    MFIXRestarter mfix_restarter(mfix_coarse->nlev());

    mfix_restarter.allocate_avgdPIC_coarse(mfix_coarse);

    mfix_restarter.deposit_PIC(mfix_coarse);

    // ************************************************************************
    // Free up unnecessary memory
    // ************************************************************************

    // After having deposited PIC quantities on the coarse mesh, we can free up
    // memory space by deleting the PIC particles container of the coarse mfix
    delete mfix_coarse->pc;
    mfix_coarse->pc = nullptr;


    // ************************************************************************
    // Create the fine-mesh mfix object
    // ************************************************************************

    mfix mfix_fine;

    auto& rw_fine = *(mfix_fine.m_rw);
    auto& timer_fine = mfix_fine.timer();
    auto& eb_fine = *(mfix_fine.m_eb);

    mfix_fine.set_run_type(RunType::PIC2DEM);
    mfix_fine.InitParams();

    // Set fine mesh objects by refinement from the coarse mesh
    mfix_restarter.set_fine_grids_from_coarse(&mfix_fine, mfix_coarse);

    // Initialize EB geometry. This needs to be done before grid creation (in
    // mfix::Init), as the grids are created using each EB-level's volfrac.
    eb_fine.make_geometry(mfix_fine.Geom(),
                          rw_fine.restart_file);

    // Reset the timer of the fine mfix object from the one in the coarse mfix
    timer_fine.reset(timer_coarse);

    for ( int lev(0); lev < mfix_coarse->nlev(); ++lev) {

      Real time = timer_fine.time();
      BoxArray const& ba = mfix_fine.boxArray(lev);
      DistributionMapping const& dm = mfix_fine.DistributionMap(lev);

      mfix_fine.MakeNewLevelFromScratch(lev, time, ba, dm);
    }

    mfix_fine.InitParticlesFromScratch();

    // Finalize initialization
    mfix_fine.Init(timer_fine.time());

    // Sanity check
    AMREX_ALWAYS_ASSERT(mfix_fine.m_pic.solve() == 0);

    eb_fine.write_surface(mfix_fine.Geom(),
                          mfix_fine.DistributionMap(),
                          mfix_fine.boxArray());

    if (mfix_fine.m_dem.solve()) {
      eb_fine.fill_levelsets(mfix_fine.Geom(),
                             mfix_fine.pc,
                             mfix_fine.m_boundary_conditions,
                             mfix_fine.m_porous_media);
    }


    // ************************************************************************
    // Convert coarse-mesh data to fine-mesh
    // ************************************************************************

    mfix_restarter.allocate_avgdPIC_fine(&mfix_fine);

    // Here starts the part with new stuff
    // Set mfix_fine dem_solve to false so we initialize only the fluid data
    // add here the copy of fluid's coarse to fine variables in here
    mfix_restarter.convert_coarse_data(mfix_coarse, &mfix_fine);

    // Get a copy of the coarse solids volume fraction and coarse geometry
    // before we delete the coarse mfix object to free up memory
    mfix_restarter.get_eps_coarse(mfix_coarse);
    const Geometry geom_coarse(mfix_coarse->Geom(0));

    // Free up memory that is not needed anymore
    delete mfix_coarse;

    // Free up memory that is not needed anymore
    for (int lev(0); lev < mfix_coarse->nlev(); ++lev) {
      // Coarse txfr data is not needed anymore
      delete mfix_restarter.avgdPIC_coarse[lev];
      mfix_restarter.avgdPIC_coarse[lev] = nullptr;
    }

    // ************************************************************************
    // Generate DEM particles
    // ************************************************************************

    mfix_restarter.get_particles_radius(&mfix_fine);
    mfix_restarter.generate_particles(geom_coarse, &mfix_fine);
    mfix_restarter.init_particles_data(&mfix_fine);

    // Free memory allocated in get_eps_coarse
    mfix_restarter.free_eps_coarse();

    std::string check_filename = rw_fine.check_filename();

    rw_fine.WriteCheckPointFile(check_filename,
                                timer_fine.nstep(),
                                timer_fine.dt(),
                                timer_fine.time());

    rw_fine.reportGridStats();

    if (rw_fine.stop_for_unused_inputs && ParmParse::QueryUnusedInputs())
      amrex::Warning("there were unused inputs");

    Real end_time = timer_coarse.elapsed_runtime();
    ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
      std::cout << "Time spent in PIC2DEM conversion app: " << end_time << std::endl;

    amrex::Print() << " " << std::endl;

  } // This end bracket and the start bracket after Initialize are essential so
    // that the mfix object is deleted before Finalize

  amrex::Finalize();

  return 0;
}
