#include <AMReX_ParmParse.H>
#include <AMReX_EBAmrUtil.H>

#include <mfix.H>
#include <mfix_reporter.H>
#include <mfix_init_fluid.H>
#include <mfix_regions.H>
#include <mfix_utils.H>

using namespace amrex;

void
mfix::InitParams ()
{
  timer().Initialize();
  regions.Initialize();

  // Read and process species, fluid and DEM particle model options.
  species.Initialize();
  reactions.Initialize(species);
  fluid.Initialize(species, reactions);
  solids.Initialize(species, reactions);

  m_dem.Initialize();
  m_pic.Initialize();

  // Read in regions, initial and boundary conditions. Note that
  // regions need to be processed first as they define the
  // physical extents of ICs and BCs.
  bcs().Initialize(regions, fluid, solids, m_dem, m_pic);

  ics().Initialize(regions, fluid, solids, m_dem, m_pic);

  pms().Initialize(regions, geom);

  tag().Initialize(max_level, fluid.solve(), has_particles(), regions);

  m_rw->Initialize();

  m_mass_balance->Initialize();

  m_loadbalance.reset( new LoadBalance(max_level, fluid.solve(), has_particles()));

  // Verbosity and MLMG parameters are now ParmParse with "nodal_proj" in the
  nodalproj_options = std::make_unique<MfixUtil::MLMGOptions>("nodal_proj");

  { ParmParse pp("fluid.newton_solver");

    pp.query("absolute_tol", newton_abstol);
    pp.query("relative_tol", newton_reltol);
    pp.query("max_iterations", newton_maxiter);
  }

  { ParmParse pp("mfix");

    // Flag to set verbosity
    m_verbose = 0;
    pp.query("verbose", m_verbose);

    pp.query("ooo_debug", ooo_debug);

    // Temporary keyword for testing new depdt algorithm
    pp.query("use_new_depdt_algo", m_use_new_depdt_algo);

    Array<Real,3> gravity_in{0.0, 0.0, 0.0};
    pp.get("gravity", gravity_in);
    for (int dir = 0; dir < 3; dir++)
    { gravity[dir] = gravity_in[dir]; }

    // frequency and bin size for sorting particles
    pp.query("particle_sorting_bin", particle_sorting_bin);
    sort_particle_int = -1;
    pp.query("sort_particle_int", sort_particle_int);

    // Options to control initial projections (mostly we use these for
    // debugging)
    pp.query("initial_iterations", m_initial_iterations);
    pp.query("do_initial_proj", do_initial_proj);

    pp.query("test_tracer_conservation", test_tracer_conservation);
    if (test_tracer_conservation && !fluid.solve_tracer()) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Flag to test tracer conservation is true but tracers are not solved.";
    }

    // Are we advecting velocity or momentum (default is velocity)
    pp.query("advect_momentum"                  , m_advect_momentum);

    pp.query("grid_pruning", m_grid_pruning);

    pp.query("regrid_int", m_regrid_int);
    if ( m_regrid_int <= 0 ) {

      m_regrid_int = 0;

    } else if ( (ParallelDescriptor::NProcs() == 1) && (max_level == 0) ) {
      m_regrid_int = 0;
      reporter::Log(reporter::Warning)
          << "Regridding disabled for single level, single rank run.";
    }

    { // This option should not have been added to main code.
      FixInputs fix("Mar. 2025");
      fix.remove("mfix.use_drag_coeff_in_proj_gp");
    }

    { // This option should not have been added to main code.
      FixInputs fix("Jul. 2025");
      fix.remove("mfix.use_drag_in_godunov");
    }

    {
      FixInputs fix("Nov. 2025");
      // I can't stand seeing "axe"
      fix.swap<std::string>("cylinder.rotation_axe", "cylinder.rotation_axis");
      // Consistency
      fix.swap<std::string>("mfix.use_mac_phi_in_godunov", "mfix.godunov_use_mac_phi");
      fix.swap<std::string>("mfix.godunov_ppm", "mfix.godunov_use_ppm");
    }

    // Redistribute before the nodal projection
    pp.query("redistribute_before_nodal_proj"   , m_redistribute_before_nodal_proj);

    // Redistribute after the nodal projection
    pp.query("redistribute_nodal_proj"          , m_redistribute_nodal_proj);

    // Redistribute after the initial nodal projection
    pp.query("redistribute_after_initial_nodal_proj", m_redistribute_after_initial_nodal_proj);

    // Threshold volfrac for correcting small cell velocity in the predictor and corrector
    pp.query("correction_small_volfrac"         , m_correction_small_volfrac);

    // Are we using MOL or Godunov?
    std::string l_advection_type = "Godunov";
    pp.query("advection_type"                   , l_advection_type);
    pp.query("godunov_use_ppm"                  , m_godunov_ppm);
    pp.query("godunov_use_forces_in_trans"      , m_godunov_use_forces_in_trans);
    pp.query("godunov_include_diff_in_forcing"  , m_godunov_include_diff_in_forcing);
    pp.query("godunov_use_mac_phi"              , m_godunov_use_mac_phi);

    // agglomeration for GMG coarse levels
    pp.query("agg_grid_size", agg_grid_size);

    pp.query("redistribution_type", m_redistribution_type);
    if (m_redistribution_type != "NoRedist" &&
        m_redistribution_type != "FluxRedist" &&
        m_redistribution_type != "StateRedist") {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Invalid redistribution type. Valid options include:\n"
          << "   FluxRedist, NoRedist or StateRedist";
    }

    // Default to Godunov
    if(amrex::toLower(l_advection_type).compare("mol") == 0) {
      m_advection_type = AdvectionType::MOL;
    } else if(amrex::toLower(l_advection_type).compare("godunov") == 0) {
      m_advection_type = AdvectionType::Godunov;
    } else {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Invalid advection scheme. Valid options include:\n   MOL, Godunov";
    }

    if (advection_type() == AdvectionType::MOL)
    { m_godunov_include_diff_in_forcing = false; }

    // MOL: Explicit predictor / Crank_Nicolson corrector
    // Godunov: Implicit predictor / No corrector
    if (advection_type() == AdvectionType::MOL) {
      m_predictor_diff_type = DiffusionType::Explicit;
      m_corrector_diff_type = DiffusionType::Crank_Nicolson;
    } else {
      m_predictor_diff_type = DiffusionType::Implicit;
      m_corrector_diff_type = DiffusionType::Undefined;
    }


    // Options to control time stepping
    // MOL: default CFl = 0.5
    // Godunov: default CFL = 0.9
    m_cfl = (advection_type() == AdvectionType::MOL) ? 0.5 : 0.9;
    pp.query("cfl", m_cfl);

    if (advection_type() == AdvectionType::MOL && m_cfl > 0.5) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "The provide cfl value " << m_cfl << " is invalid.\n"
          << "Method-of-Lines advection scheme requires cfl <= 0.5";
    }
    if (advection_type() == AdvectionType::Godunov && m_cfl > 1.0) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "The provide cfl value " << m_cfl << " is invalid.\n"
          << "Godunov advection scheme requires cfl <= 1.0";
    }

    // This will multiply the time-step in the very first step only
    pp.query("scale_init_dt", m_scale_init_dt);
    if (m_scale_init_dt > 1.0) {
      amrex::Abort("We require scale_init_dt <= 1.0");
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Setting scale_init_dt > 1. is invalid.";
    }

    // Verbosity and MLMG parameters are now ParmParse with "mac_proj" in the
    // inputs file
    // Examples: mac_proj.verbose = 1
    //           mac_proj.bottom_verbose = 1
    //           mac_proj.maxiter
    //           mac_proj.bottom_maxiter
    //           mac_proj.bottom_rtol
    //           mac_proj.bottom_atol
    //           mac_proj.bottom_solver
    // More info at "AMReX-Hydro/Projections/hydro_MacProjector.cpp"
    macproj_options = std::make_unique<MfixUtil::MLMGOptions>("mac_proj");

    // Checks for hypre namespace
    if (nodalproj_options->bottom_solver_type == "hypre" &&
          macproj_options->bottom_solver_type == "hypre") {

      std::string nodal_ns = nodalproj_options->hypre_namespace;
      std::string mac_ns = macproj_options->hypre_namespace;

      if (nodal_ns == "hypre" && mac_ns != "hypre") {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "hypre namespace required for nodal projection";
      }

      if (nodal_ns != "hypre" && mac_ns == "hypre") {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "hypre namespace required for MAC projection";
      }

      if ((nodal_ns == mac_ns) && nodal_ns != "hypre") {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "same hypre namespace other than \"hypre\" not allowed\n"
            << "for nodal and MAC projections";
      }
    }
  }

  if (has_particles()) {

    ParmParse pp_particles("particles");

    // Keep particles that are initially touching the wall. Used by DEM tests.
    pp_particles.query("removeOutOfRange", removeOutOfRange);
    pp_particles.query("reduceGhostParticles", reduceGhostParticles);

    if (!fluid.solve()) {

      if (m_timer.timestep_type() != MFIXTimer::TimestepType::Fixed) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Particle-only runs require a positive fixed_dt to be specified";
      }
    }
  }


  const int solve_fluid_and_solids =
      static_cast<int>( fluid.solve() && has_particles() );

  m_coupling.Initialize(fluid.solve_enthalpy(), fluid.nspecies(),
      fluid.isMixture(), solve_fluid_and_solids);

  {
    ParmParse reports_pp("mfix.reports");

    reports_pp.query("mass_balance_int", m_rw->mass_balance_report_int);
    reports_pp.query("mass_balance_per_approx", m_rw->mass_balance_report_per_approx);

    if ((m_rw->mass_balance_report_int > 0 && m_rw->mass_balance_report_per_approx > 0) ) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Invalid mass balance write frequency. Specify only one:\n"
          << "mass_balance_int or mass_balance_report_per_approx";
    }
    // Add check to turn off report if not solving species
    if (fluid.solve_species() && fluid.nspecies() >= 1) {
      m_rw->report_mass_balance = (m_rw->mass_balance_report_int > 0 ||
                                     m_rw->mass_balance_report_per_approx > 0);
    } else {
      m_rw->report_mass_balance = 0;
    }
  }

}



void
mfix::Init (Real a_time)
{
  if (ooo_debug) { Print() << "Init\n"; }

  m_rw->InitIOPltData();

  if ( run_type(RunType::Standard) ) {

    finest_level = 0;

    /****************************************************************************
     *                                                                          *
     * Generate levels using ErrorEst tagging.                                  *
     *                                                                          *
     ***************************************************************************/

    // This tells the AmrMesh class not to iterate when creating the initial
    // grid hierarchy
    SetIterateToFalse();

    // This tells the Cluster routine to use the new chopping routine which
    // rejects cuts if they don't improve the efficiency
    SetUseNewChop();

    // InitFromScratch (AmrCore)
    // ├── MakeNewGrids (AmrMesh)
    // │   ├── ba = MakeBaseGrids
    // │   │   ├── ChopGrids(0, ba, ParallelDescriptor::NProcs());
    // │   │   └── PostProcessBaseGrids(ba);
    // │   └── dm = MakeDistributionMap(0,ba)
    // └── MakeNewLevelFromScratch (MFIX-Exa)

    InitFromScratch(a_time);

    therm_p = therm_po = fluid.thermodynamic_pressure();

    /****************************************************************************
    *                                                                           *
    * Particle data creation                                                    *
    *                                                                           *
    ****************************************************************************/
    if (has_particles() && !m_rw->only_print_grid_report ) {
      InitParticlesFromScratch();
    }

  // Setup for restart
  } else if ( run_type(RunType::Restart) ) {

    // NOTE: mfix::levelset_restart == true loading level-set from a
    // checkpoint file. However, if this is a replicating restart,
    // mfix::levelset_restart is set to false again (so that the level-sets
    // are recomputed for the replicated system).
    levelset_restart = true;

    // NOTE: during replication 1) this also re-builds ebfactories and
    // level-set 2) this can change the grids
    IntVect Nrep(m_rw->repl_x, m_rw->repl_y, m_rw->repl_z);
    Restart( m_rw->restart_file, m_timer.nstep(), m_timer.dt(), m_timer.time(), Nrep);

  }

  // Skip the remaining initialization if we are only writing the grid report.
  m_rw->reportGridStats();

  // Write out EB surface
  eb()->write_surface(geom, dmap, grids);

  if (m_rw->only_print_grid_report) { return; }

  m_rw->set_nlev(finest_level+1);
  m_rw->setReportTime(timer().time());

  if ( has_particles() ) {

    AMREX_ALWAYS_ASSERT( pc );

    int const pc_var_count( SoArealData::count + pc->m_runtimeRealData.count);
    avg_pc_parms().Initialize( pc_var_count, fluid, eb_parms() );

    // Updating m_rw pc pointer is needed since mfix pc has changed
    m_mass_balance->set_pc(pc);
    m_rw->set_pc(pc);
    m_rw->set_monitors_pc(pc);

    eb()->fill_levelsets(geom, pc, bcs(), pms());

    // Auto generated particles may be out of the domain. This call will
    // remove them. Note that this has to occur after the EB geometry is
    // created. if (particle_init_type == "Auto" && not restarting &&
    // particle_ebfactory[nlev()-1])

    if ( run_type(RunType::Standard) ) {

      if (removeOutOfRange) {

        int const is_pic2dem = run_type(RunType::PIC2DEM) ? 1 : 0;

        for ( int lev(0); lev <= pc->finestLevel(); ++lev) {

          if ( m_eb->particle_factory()[nlev()-1]) {
              //___________________________________________________________________
              // Only 1 refined level-set

              Print() << "Clean up auto-generated particles.\n";

              const MultiFab * ls_data = m_eb->level_sets()[1].get();
              iMultiFab ls_valid(ls_data->boxArray(), ls_data->DistributionMap(),
                                 ls_data->nComp(), ls_data->nGrow());

              pc->RemoveOutOfRange(nlev()-1, m_eb->particle_factory()[nlev()-1].get(),
                  ls_data, m_eb->levelset_refinement(), is_pic2dem, m_porous_media);

          } else if (!run_type(RunType::Restart) && m_eb->particle_factory()[nlev()-1]) {

            //___________________________________________________________________
            // Multi-level everything

            Print() << "Clean up auto-generated particles.\n";

            for (int ilev = 0; ilev < nlev(); ilev ++) {

              const MultiFab * ls_data = m_eb->level_sets()[ilev].get();
              iMultiFab ls_valid(ls_data->boxArray(), ls_data->DistributionMap(),
                                 ls_data->nComp(), ls_data->nGrow());

              pc->RemoveOutOfRange(ilev, m_eb->particle_factory()[ilev].get(),
                                   ls_data, 1, is_pic2dem, m_porous_media);
            }
          }
        }
      } // remove out of range particles
    }  // standard run

    pc->setSortingBinSizes(IntVect(particle_sorting_bin));
    pc->setSortInt(sort_particle_int);

    if (m_dem.solve() && !run_type(RunType::PIC2DEM) ) {
        pc->MFIX_PC_InitCollisionParams();
        pc->setReduceGhostParticles(reduceGhostParticles);
    }

    if (run_type(RunType::Standard)) {
      pc->InitParticlesRuntimeVariables(fluid.solve_enthalpy());
    }

  } // has_particles


  // Initialize weights for load balancing. This loop would have some
  // issues with AMR+particles because the number of levels may not
  // be the same.
  for (int lev(0); lev < nlev(); lev++) {

    if ( m_loadbalance->DualGrid() || !fluid.solve() ) {

      m_loadbalance->ResetWeights(lev, pc->ParticleBoxArray(lev),
          pc->ParticleDistributionMap(lev));

    } else { AMREX_ASSERT( m_loadbalance->SingleGrid() );

      m_loadbalance->ResetWeights(lev, grids[lev], dmap[lev]);
    }

    if (m_regrid_int > 0 && has_particles() &&
       (m_loadbalance->DualGrid()              ||
        m_loadbalance->WeightByParticleCount() ||
        m_loadbalance->WeightByParticleRunTime()) ) {

      using MFIXParIter = MFIXParticleContainer::MFIXParIter;
      for (MFIXParIter pti(*pc, lev); pti.isValid(); ++pti) {
        Real const weight(pti.numParticles() );
        m_loadbalance->UpdateWeights(lev, pti.index(), weight);
      }

      // Regrid here if it won't happen in first pass through main. This
      // should only happen for restarts where the restart step does
      // not align with the regrid interval.
      if ( m_timer.nstep()%m_regrid_int != 0 ) { Regrid(); }

    // Not sure if this logical branch is even needed.
    } else if (has_particles() && m_loadbalance->DualGrid() ) {

      bcs().calc_bc_areas(lev, pc->ParticleBoxArray(lev),
          pc->ParticleDistributionMap(lev),
          eb()->particle_factory()[lev].get());
    }
  }

  if (fluid.solve() ) {

    set_gp0();

    init_advection();

    Real const time( m_timer.time() );

    // Compute the initial volume fraction from particles
    // and any porous media
    if ( run_type(RunType::Standard) ) {

      if (has_particles()) {

        sum_vol_orig = deposit_volume_to_grid(time);
        Print() << "Setting original sum_vol to " << sum_vol_orig << '\n';

        for( int lev(0); lev < nlev(); ++lev) {
          MultiFab::Copy(*(leveldata().epf_old(lev)), *(leveldata().epf(lev)),
            0, 0, 1, leveldata().epf(lev)->nGrow());
        }
      }

      // Add in porous media contributions
      if ( pms().nregions() > 0) {
        pms().impose_vfrac(leveldata().epf());
        pms().impose_vfrac(leveldata().epf_old());
      }
    }


    for( int lev(0); lev < nlev(); ++lev) {

      int const ng = nghost_state();

      // We must fill internal ghost values before calling redistribution
      leveldata().epf(lev)->FillBoundary(geom[lev].periodicity());
      bcs().fillpatch(lev, time, BCFillVar::epf, leveldata().epf(), ng);

      leveldata().rho(lev)->FillBoundary(geom[lev].periodicity());
      bcs().fillpatch(lev, time, BCFillVar::rho, leveldata().rho(), ng);

      leveldata().vel(lev)->FillBoundary(geom[lev].periodicity());
      bcs().fillpatch(lev, time, BCFillVar::vel, leveldata().vel(), ng);

      if (leveldata(lev)->has_tracer()) {
        leveldata().tracer(lev)->FillBoundary(geom[lev].periodicity());
        bcs().fillpatch(lev, time, BCFillVar::tracer, leveldata().tracer(), ng);
      }
      if (leveldata(lev)->has_enthalpy()) {
        leveldata().h(lev)->FillBoundary(geom[lev].periodicity());
        bcs().fillpatch(lev, time, BCFillVar::h, leveldata().h(), ng);
      }
      if (leveldata(lev)->has_temperature()) {
        leveldata().T(lev)->FillBoundary(geom[lev].periodicity());
        bcs().fillpatch(lev, time, BCFillVar::T, leveldata().T(), ng);
      }
      if (leveldata(lev)->has_species()) {
        leveldata().X(lev)->FillBoundary(geom[lev].periodicity());
        bcs().fillpatch(lev, time, BCFillVar::X, leveldata().X(), ng);
      }

      // Make sure to fill the "old state" before we start.
      leveldata(lev)->resetOldWithNew();

      // resetOldWithNew does not include epf so we do it here
      if ( has_particles() ) {
        MultiFab::Copy(*(leveldata().epf_old(lev)), *(leveldata().epf(lev)),
          0, 0, 1, leveldata().epf(lev)->nGrow());
      }
    }


    if ( run_type(RunType::Standard) ) {

      InitialRedistribution(time);

      // Project the initial velocity field
      if (do_initial_proj) { project_velocity(); }

      // Iterate to compute the initial pressure
      if (m_initial_iterations > 0) { initial_iterations(); }

      if (m_run_type != RunType::PIC2DEM) {
        m_mass_balance->InitMassBalance(m_rw->report_mass_balance, fluid.nspecies());
      }
    }

    if ( has_particles() ) {

      if ( !run_type(RunType::PIC2DEM) ) {

        EulerianMonitor::VolumeIntegral monitor( leveldata(),
            m_eb->factory(), geom[0].ProbDomain());

        // we compute the monitor on the only component of the epf Multifab
        const int comp = 0;
        auto domain_vol = monitor.volume_integral(0, {leveldata().epf(0)}, {comp});

        // domain_vol is a vector, but only has one variable. We get the value
        // at index 0
        const int var = 0;

        Print() << "Enclosed domain volume is   " << domain_vol[var] << '\n';

      }
    }
  } // solving fluid
}



void mfix::
PruneBaseGrids (BoxArray &ba) const
{
  if (m_verbose > 0) { Print() << "Pruning base grids\n"; }

  // Use 1 ghost layer
  EBDataCollection ebdc(*(m_eb->levels()[0]), geom[0],
      ba, DistributionMapping{ba}, {1}, EBSupport::basic);

  const auto &cflag = ebdc.getMultiEBCellFlagFab();
  Vector<Box> uncovered;

  for (MFIter mfi(cflag); mfi.isValid(); ++mfi) {

    FabType t = cflag[mfi].getType();
    const Box& vbx = mfi.validbox();
    if (t != FabType::covered) {
      uncovered.push_back(vbx);
    }
  }

  amrex::AllGatherBoxes(uncovered);
  ba = BoxArray(BoxList(std::move(uncovered)));
}


void
mfix::ChopGrids (const Box& domain,
                 BoxArray& ba,
                 int target_size) const
{
  if (m_verbose > 0) { Print() << "Chopping grids\n"; }

  reporter::Log(reporter::Warning)
      << "Using max_grid_size did not make enough grids for the\n"
      << "number of processors";

  // Here we hard-wire the maximum number of times we divide the boxes.
  int n = 10;

  // Here we hard-wire the minimum size in any one direction the boxes can be
  int min_grid_size = 4;

  IntVect chunk(domain.length(0),domain.length(1),domain.length(2));

  int j = -1;

  for (int cnt = 1; cnt <= n; ++cnt) {

    if (chunk[0] >= chunk[1] && chunk[0] >= chunk[2])      { j = 0; }
    else if (chunk[1] >= chunk[0] && chunk[1] >= chunk[2]) { j = 1; }
    else if (chunk[2] >= chunk[0] && chunk[2] >= chunk[1]) { j = 2; }

    chunk[j] /= 2;

    if (chunk[j] >= min_grid_size) {

        ba.maxSize(chunk);

    } else {

      reporter::Log(reporter::Warning)
          << "ChopGrids was unable to make enough grids for the\n"
          << "number of processors";

      return;
    }

    // Test if we now have enough grids
    if (ba.size() >= target_size) { return; }
  }
}


void mfix::
PostProcessBaseGrids( BoxArray& a_grids) const
{
  if ( m_verbose > 0 ) { Print() << "PostProcessBaseGrids\n"; }
  if (m_grid_pruning) { PruneBaseGrids(a_grids); }
}


void mfix::
MakeNewLevelFromScratch ( int a_lev, Real /*a_time*/,
                          const BoxArray& a_new_grids,
                          const DistributionMapping& a_new_dmap)
{
  if ( m_verbose > 0 ) {
    Print() << "Making new level " << a_lev << " from scratch\n";
    if (m_verbose > 1) { Print() << "with box array " << a_new_grids << '\n'; }
  }

  SetBoxArray(a_lev, a_new_grids);
  SetDistributionMap(a_lev, a_new_dmap);

  m_eb->make_factory(a_lev, geom[a_lev], grids[a_lev], dmap[a_lev]);

  // Skip allocating and initializing grid data if we are only
  // generating the grid report.
  if ( m_rw->only_print_grid_report ) { return; }

  // Initialize the level data container.
  m_level_data.m_level_data[a_lev] =
      std::make_unique<LevelData>(fluid, reactions);

  leveldata().define(a_lev, nghost_state(), grids[a_lev], dmap[a_lev],
      *(m_eb->factory()[a_lev]) );

  leveldata(a_lev)->initVals( init_value );

  // Create the BC data structures
  bc_list.MakeBCArrays(a_lev, nghost_state(), geom[a_lev], ooo_debug);

  bcs().calc_bc_areas(a_lev, grids[a_lev], dmap[a_lev], m_eb->factory()[a_lev].get());

  bcs().set_bc_type(a_lev, nghost_state(), fluid);
  bcs().set_bc_list(a_lev, nghost_state());

  // Impose the initial condition on the new level. These routines need reworked.
  if ( run_type(RunType::Standard) && fluid.solve() ) {

    Box domain(geom[a_lev].Domain());

    MultiFab& epf = *(leveldata().epf(a_lev));

    Real dx = geom[a_lev].CellSize(0);
    Real dy = geom[a_lev].CellSize(1);
    Real dz = geom[a_lev].CellSize(2);

    const GpuArray<Real,3> plo = geom[a_lev].ProbLoArray();

    // We deliberately don't tile this loop since we will be looping
    //    over bc's on faces and it makes more sense to do this one grid at a time
    for (MFIter mfi(epf, false); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.validbox();
      const Box& sbx = epf[mfi].box();

      set_ic_fluid(a_lev, sbx, bx, domain, mfi, leveldata(), dx, dy, dz, plo,
          test_tracer_conservation, m_initial_conditions, fluid);
    }
  } // not restarting
}



void mfix::
InitParticlesFromScratch ()
{
  Real strt_init_part = ParallelDescriptor::second();

  BoxList             pbl;
  BoxArray            pba;
  DistributionMapping pdm;

  if ( m_loadbalance->SingleGrid() ) {

    pbl = boxArray(0).boxList();
    pba = boxArray(0);
    pdm = DistributionMap(0);

  } else { // dual grid is enabled

    ParmParse pp("particles");

    IntVect particle_max_grid_size;

    const int contains_x = pp.query("max_grid_size_x", particle_max_grid_size[0]);
    const int contains_y = pp.query("max_grid_size_y", particle_max_grid_size[1]);
    const int contains_z = pp.query("max_grid_size_z", particle_max_grid_size[2]);

    if (!contains_x && !contains_y && !contains_z) {

      reporter::Log(reporter::Warning)
          << "Setting particle max grid size to fluid max grid size";

      particle_max_grid_size = max_grid_size[0];

    } else if (!(contains_x && contains_y && contains_z)) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "The max grid size for particles must be specified for x, y and z!\n"
          << "Please correct the input deck.";

    }

    // Define coarse level BoxArray and DistributionMap
    //const BoxArray& ba = MakeBaseGrids(geom[0].Domain(), particle_max_grid_size);
    {
      if (ooo_debug) { Print() << "MakeBaseGrids - Particles\n"; }

      pba = BoxArray(geom[0].Domain());
      pba.maxSize(max_grid_size[0]);

      if (m_grid_pruning) { PruneBaseGrids(pba); }

      // We only call ChopGrids if dividing up the grid using max_grid_size didn't
      //    create enough grids to have at least one grid per processor.
      // This option is controlled by "refine_grid_layout" which defaults to true.
      if (refine_grid_layout && pba.size() < ParallelDescriptor::NProcs())
      { ChopGrids(geom[0].Domain(), pba, ParallelDescriptor::NProcs()); }

      // avoid duplicates
      if (pba == grids[0]) { pba = grids[0]; }

      Print() << "In MakeBaseGrids: BA HAS " << pba.size() << " GRIDS\n";
    }

    pdm.define(pba, ParallelDescriptor::NProcs());
    pbl = pba.boxList();

    // chop grids to use all the gpus for particle generation
    if (pba.size() < ParallelDescriptor::NProcs()) {

      IntVect reduced_size = particle_max_grid_size;

      while (pbl.size() < ParallelDescriptor::NProcs()) {

        pbl.clear();

        int maxdir = reduced_size.maxDir(false);
        reduced_size[maxdir] /= 2;

        for (auto i=0; i<pba.size(); i++) {

          BoxArray tmpba{pba[i]};
          tmpba.maxSize(reduced_size);
          pbl.join(tmpba.boxList());
        }
      }

      pba.define(pbl);
      pdm.define(pba, ParallelDescriptor::NProcs());
    }
  }

  pc = new MFIXParticleContainer(geom, pdm, pba, ics(), bcs(),
      solids, m_dem, m_pic, fluid, reactions,
      m_coupling.include_virtual_mass());

  pc->AllocData();

  if ( m_loadbalance->DualGrid() ) {
    for ( int lev(0); lev <= pc->finestLevel(); ++lev) {
      m_eb->make_particle_factory(lev, geom[lev], pc->ParticleBoxArray(lev),
        pc->ParticleDistributionMap(lev));
    }
  }

  if (run_type(RunType::Standard)) {
    if (ics().AutoParticleInit()) {
      Print() << "Auto generating particles ...\n";
      pc->InitParticlesAuto(m_eb->particle_factory()[0].get());
    } else {
      Print() << "Reading particles from particle_input.dat ...\n";
      pc->InitParticlesAscii("particle_input.dat");
    }

    pc->Redistribute(0, 0, 0, 0, false);
  }

  Real end_init_part = ParallelDescriptor::second() - strt_init_part;

  ParallelDescriptor::ReduceRealMax(end_init_part,
      ParallelDescriptor::IOProcessorNumber());

  Print() << "Time spent in initializing particles " << end_init_part << '\n';

}
