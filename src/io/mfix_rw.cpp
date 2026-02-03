#include <mfix.H>
#include <mfix_rw.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

#include <AMReX_ParmParse.H>
#include <AMReX_EBFArrayBox.H>

#include <string>
#include <sstream>
#include <unordered_set>
#include <mfix_fix_inputs.H>

#ifdef AMREX_USE_CATALYST
#include "catalyst.hpp"
#include "AMReX_Conduit_Blueprint.H"
#endif


using namespace amrex;


MFIXReadWrite::MFIXReadWrite (int nlev_in,
                              amrex::Vector<amrex::BoxArray>& grids_in,
                              amrex::Vector<amrex::Geometry>& geom_in,
                              MFIXParticleContainer* pc_in,
                              MFIXFluidPhase& fluid_in,
                              MFIXLevelData& a_leveldata,
                              amrex::Real const& a_therm_p,
                              amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory_in,
                              amrex::Vector<amrex::DistributionMapping>& dmap_in,
                              bool ooo_debug_in,
                              amrex::Vector<std::unique_ptr<amrex::MultiFab>> const& level_sets_in,
                              int levelset_refinement_in,
                              int levelset_pad_in,
                              int levelset_eb_refinement_in,
                              int levelset_eb_pad_in,
                              MFIXSolidsPhase& solids_in,
                              MFIXDEM& dem,
                              MFIXPIC& pic,
                              MFIXReactions& reactions_in,
                              const amrex::Vector<amrex::IntVect>& ref_ratio_in,
                              BCList& bc_list_in,
                              Vector<std::shared_ptr<EBFArrayBoxFactory>> const& particle_ebfactory_in,
                              std::shared_ptr<MFIXMassBalance> a_mass_balance,
                              MFIXRegions& regions_in,
                              MFIXBoundaryConditions& boundary_conditions_in,
                              const amrex::Vector<amrex::BCRec>& bcrec_velocity_in)
  : nlev(nlev_in)
  , grids(grids_in)
  , geom(geom_in)
  , pc(pc_in)
  , fluid(fluid_in)
  , m_leveldata(a_leveldata)
  , m_therm_p(a_therm_p)
  , ebfactory(ebfactory_in)
  , dmap(dmap_in)
  , ooo_debug(ooo_debug_in)
  , level_sets(level_sets_in)
  , levelset_refinement(levelset_refinement_in)
  , levelset_pad(levelset_pad_in)
  , levelset_eb_refinement(levelset_eb_refinement_in)
  , levelset_eb_pad(levelset_eb_pad_in)
  , solids(solids_in)
  , m_dem(dem)
  , m_pic(pic)
  , reactions(reactions_in)
  , ref_ratio(ref_ratio_in)
  , bc_list(bc_list_in)
  , particle_ebfactory(particle_ebfactory_in)
  , m_mass_balance(a_mass_balance)
  , regions(regions_in)
  , m_boundary_conditions(boundary_conditions_in)
  , m_bcrec_velocity(bcrec_velocity_in)
  , m_ascent_actions_yaml("")
{

  // These are all fixes to get around regtest that sets them.
  { FixInputs fix(0);

    fix.swap<int>("amr.plt_regtest", "mfix.plt_regtest", 1);
    fix.swap<int>("amr.plot_int","mfix.plot_int",1);
    fix.swap<std::string>("amr.plot_file","mfix.plot_file",1);

    fix.swap<bool>("amr.checkpoint_files_output","mfix.checkpoint_files_output",1);
    fix.swap<int>("amr.check_int", "mfix.check_int",1);
    fix.swap<std::string>("amr.check_file", "mfix.check_file",1);

    fix.swap<std::string>("amr.restart", "mfix.restart",1);
  }

  readParameters();

#ifdef AMREX_USE_CATALYST
  if (catalyst_enabled) {

    conduit_cpp::Node params;

    params["catalyst/scripts/script0"].set_string(catalyst_script);
    params["catalyst_load/implementation"].set_string(catalyst_impl);
    params["catalyst_load/search_paths/paraview"].set_string(catalyst_library_path);

    catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&params));

    if (err != catalyst_status_ok) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Catalyst failed to initialize";
    } else {
      reporter::Log(reporter::Status)
        << "Catalyst was successfully initialized.";
    }
  }
#endif
}


MFIXReadWrite::~MFIXReadWrite ()
{
#ifdef AMREX_USE_CATALYST
  if (catalyst_enabled) {
      conduit_node* f_params = conduit_node_create();
      catalyst_finalize(f_params);
  }
#endif
}



void MFIXReadWrite::readParameters ()
{
  {
     ParmParse pp("mfix");

     // Minimum number of digits in all filenames
     pp.query("file_digits", m_min_digits);

     // Checkpoint output control
     pp.query("checkpoint_files_output", checkpoint_files_output);
     pp.query("check_file", check_file);
     pp.query("check_int", check_int);

     //pp.query("check_per_exact", check_per_exact);
     pp.query("check_per_approx", check_per_approx);

     std::string walltime_interval_in;
     int has_walltime_interval = pp.query("check_walltime_interval", walltime_interval_in);

     if (has_walltime_interval) {
       int HH(0), MM(0), SS(0);
       if (sscanf(walltime_interval_in.c_str(), "%d:%d:%d", &HH, &MM, &SS) >= 2) {
         check_walltime_interval = static_cast<Real>(HH*3600 + MM*60 + SS);
       } else {
         std::string message =
           " Error: Unable to correctly parse checkpoint walltime interval "
           + walltime_interval_in + "\n" + " The correct format is HH:MM:SS\n";
         amrex::Print() << message;
         amrex::Abort(message);
       }
     }

     if (/*(check_int        > 0 && check_per_exact  > 0) ||*/
         (check_int        > 0 && check_per_approx        > 0) ||
         (check_int        > 0 && check_walltime_interval > 0) ||
         (check_per_approx > 0 && check_walltime_interval > 0) /*||
         (check_per_exact  > 0 && check_per_approx > 0) */ )
       amrex::Abort("Must choose only one of check_int, check_per_approx, or "
           "check_walltime_interval");

     // Plot output control
     pp.query("plot_file", plot_file);
     pp.query("plotfile_on_restart", plotfile_on_restart);
     pp.query("plot_int", plot_int);

     //pp.query("plot_per_exact", plot_per_exact);
     pp.query("plot_per_approx", plot_per_approx);

     if (/*(plot_int       > 0 && plot_per_exact  > 0) ||*/
         (plot_int       > 0 && plot_per_approx > 0) /*||
         (plot_per_exact > 0 && plot_per_approx > 0) */ )
       amrex::Abort("Must choose only one of plot_int or plot_per_approx");

     // Plot solids only in specific regions
     {
       ParmParse ppSolids("mfix.solids");

       std::vector<std::string> solids_regions;
       ppSolids.queryarr("regions", solids_regions);

       if (solids_regions.size() > 0 &&
           amrex::toLower(solids_regions[0]).compare("none") != 0) {

         m_solids_plot_regions.resize(solids_regions.size());

         // Loop over input plot regions
         for (size_t n(0); n < solids_regions.size(); n++) {

           auto& plot_region = m_solids_plot_regions[n];

           const std::string& region_name = solids_regions[n];
           plot_region.m_region_name = region_name;

           ppSolids.queryarr(region_name, plot_region.m_plot_names);

           const std::string solids_region_prefix = "mfix.solids." + region_name;
           ParmParse ppSolidsRegion(solids_region_prefix);

           plot_region.m_pp_string = solids_region_prefix;

           int ppint = ppSolidsRegion.query("plot_int", plot_region.m_plot_int);
           int ppapprox = ppSolidsRegion.query("plot_per_approx", plot_region.m_plot_per_approx);

           AMREX_ALWAYS_ASSERT(ppint || ppapprox);

           ppSolidsRegion.queryarr("plt_fluid_vars", plot_region.m_plot_fluid_vars);
         }
       }
     }

     // Monitors
     m_monitors.initialize(regions, m_leveldata, ebfactory, fluid,
         pc, particle_ebfactory, solids);

     // Ascent output control
     pp.query("ascent_on_restart", ascent_on_restart);

     pp.query("ascent_int", ascent_int);
     pp.query("ascent_per_approx", ascent_per_approx);

     if ((ascent_int > 0 && ascent_per_approx > 0) )
       amrex::Abort("Must choose only one of ascent_int or ascent_per_approx");

     pp.query("par_ascii_file", par_ascii_file);
     pp.query("par_ascii_int", par_ascii_int);
     pp.query("par_ascii_per_approx", par_ascii_per_approx);

     pp.query("restart", restart_file);

     ParmParse("pic2dem").query("convert", restart_file);

     pp.query("repl_x", repl_x);
     pp.query("repl_y", repl_y);
     pp.query("repl_z", repl_z);

  }

  {
     ParmParse pp("mfix");

     pp.query("write_ls", write_ls);
     pp.query("stop_for_unused_inputs", stop_for_unused_inputs);
     pp.query("only_print_grid_report", only_print_grid_report);
     pp.query("plt_geom", plt_geom);
  }

  {
     ParmParse pp("ascent");

     pp.query("actions", m_ascent_actions_yaml);
  }

#ifdef AMREX_USE_CATALYST
  {
    ParmParse pp("catalyst");
    pp.query("catalyst_on_restart", catalyst_on_restart);
    pp.query("script", catalyst_script);
    pp.query("implementation", catalyst_impl);
    pp.query("library_path", catalyst_library_path);
    pp.query("enabled", catalyst_enabled);
  }
#endif

#ifdef MFIX_MPMD
  {
    ParmParse pp("mfix");

    pp.query("mpmd_int", m_mpmd_int);
    pp.query("mpmd_per_approx", m_mpmd_per_approx);

    std::unordered_set<std::string> unique_names;
    pp.queryarr("mpmd_static_mfs", m_mpmd_static_mf_names);
    for (const auto& name : m_mpmd_static_mf_names) {
      if ((name != "volfrac") && (name != "centroid")) {
        std::string message = "Unknown input '" + name + "' specified in mpmd_static_mfs";
        amrex::Abort(message);
      } else if (unique_names.find(name) != unique_names.end()) {
        std::string message = "Input '" + name + "' specified more than once in mpmd_static_mfs";
        amrex::Abort(message);
      } else {
        unique_names.insert(name);
      }
    }

    unique_names.clear();
    pp.queryarr("mpmd_mfs", m_mpmd_mf_names);
    for (const auto& name : m_mpmd_mf_names) {
      if ((name != "vel_g") && (name != "ep_g") && (name != "T_g")
          && (name != "X_gk")) {
        std::string message = "Unknown input '" + name + "' specified in mpmd_mfs";
        amrex::Abort(message);
      } else if (unique_names.find(name) != unique_names.end()) {
        std::string message = "Input '" + name + "' specified more than once in mpmd_mfs";
        amrex::Abort(message);
      } else {
        unique_names.insert(name);
      }
    }
  }
#endif

}


void
MFIXReadWrite::Initialize ()
{
  // Finalize initialization of solids plot regions
  for (int n(0); n < m_solids_plot_regions.size(); ++n) {

    auto& plot_region = m_solids_plot_regions[n];

    const int N_names = plot_region.m_plot_names.size();
    if (N_names > 0) {

      plot_region.m_h_plot_types.clear();
      plot_region.m_h_plot_types.resize(N_names);

      for (int i(0); i < N_names; ++i) {
        const std::string& name = plot_region.m_plot_names[i];

        int found(0);

        for (int j(0); j < solids.names().size(); ++j) {
          const std::string& solids_name = solids.names(j);

          if (name.compare(solids_name) == 0) {
            plot_region.m_h_plot_types[i] = j+1;
            found = 1;
            break;
          }
        }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(found, "Solid type not found in solids");
      }

      plot_region.m_d_plot_types.resize(N_names);
      Gpu::copy(Gpu::hostToDevice, plot_region.m_h_plot_types.begin(),
                plot_region.m_h_plot_types.end(), plot_region.m_d_plot_types.begin());
    }

    const RealBox* region_extents = regions.getRegion(plot_region.m_region_name);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(region_extents != nullptr, "Invalid solids plot region!");

    plot_region.m_region_extents = *region_extents;
  }

  m_p_reports_diameter = new reports::diameter(m_pic.solve(), regions);
}


void
MFIXReadWrite::writeNow (MFIXTimer& timer,
                         bool first,
                         bool last)
{
/*--------------------------------------------------------------------------------------------------
 *
 *                                     AMReX Plot File Output Control
 *
 *------------------------------------------------------------------------------------------------*/
    int plot_test = 0;

    if ( first ) {
        if ( (restart_file.empty() || plotfile_on_restart) &&
          (plot_int > 0 /*|| plot_per_exact > 0*/ || plot_per_approx > 0) )
          plot_test = 1;
    }

    else if ( last && plot_int > 0 ) {
        plot_test = 1;
    }

    else if (plot_per_approx > 0.0)
    {
        plot_test = test_per_approx(timer.time(), timer.dt(), plot_per_approx);

    }/*
    else if ( plot_per_exact  > 0 && (amrex::Math::abs(remainder(timer.time(), plot_per_exact)) < 1.e-12) )
    {
        plot_test = 1;
    }*/


    if ( (plot_test == 1) || ( ( plot_int > 0) && ( timer.nstep() %  plot_int == 0 ) ) )
    {
        if (fluid.solve() && mfix::m_run_type != RunType::PIC2DEM)
          ComputeVort();

        WritePlotFile(plot_file, timer.nstep(), timer.time());
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                     AMReX Solids Plot File Output Control
 *
 *------------------------------------------------------------------------------------------------*/
    if ((m_dem.solve() || m_pic.solve()) && (solids_plot_regions() == true)) {

      BL_PROFILE("mfix::WriteSolidsPlotFile()");

      for (int n(0); n <m_solids_plot_regions.size(); ++n) {

        auto& plot_region = m_solids_plot_regions[n];

        plot_test = 0;

        if (first) {
          if (//(restart_file.empty() || plotfile_on_restart) &&
              (plot_region.m_plot_int > 0 /*|| plot_region.m_plot_per_exact > 0*/ ||
               plot_region.m_plot_per_approx > 0))

            plot_test = 1;
        }

        else if (last && plot_region.m_plot_int > 0) {
          plot_test = 1;
        }

        else if (plot_region.m_plot_per_approx > 0.0)
        {
          plot_test = test_per_approx(timer.time(), timer.dt(), plot_region.m_plot_per_approx);

        }/*
        else if (plot_region.m_plot_per_exact  > 0 &&
                 (amrex::Math::abs(remainder(timer.time(), plot_region.m_plot_per_exact)) < 1.e-12) )
        {
            plot_test = 1;
        }*/


        if ((plot_test == 1) || ((plot_region.m_plot_int > 0) &&
            (timer.nstep() % plot_region.m_plot_int == 0)))
        {
          WriteSolidsPlotFile(plot_region, timer.nstep(), timer.time());
        }
      }
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                     MFIX Monitors Plot File Output Control
 *
 *------------------------------------------------------------------------------------------------*/
    for (int i(0); i < m_monitors.size(); ++i) {

      auto& monitor = m_monitors.get(i);

      int monitor_test = 0;

      if ( first || (timer.nstep() == 0) ) {
        if ( (monitor.plot_int() > 0 || monitor.plot_per_approx() > 0) )
          monitor_test = 1;
      }

      else if ( last ) {
        monitor_test = 0;
      }

      else if (monitor.plot_per_approx() > 0.0) {
        monitor_test = test_per_approx(timer.time(), timer.dt(), monitor.plot_per_approx());
      }

      else if ((monitor.plot_int() > 0) && (timer.nstep() % monitor.plot_int() == 0)) {
        monitor_test = 1;
      }

      if ( monitor_test == 1 ) {

        monitor.write_csv(timer.time(), timer.dt());
      }
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                       Ascent Output Control
 *
 *------------------------------------------------------------------------------------------------*/

#ifdef AMREX_USE_ASCENT
    int ascent_test = 0;

    if ( first ) {
        if ((restart_file.empty() || ascent_on_restart) &&
            (ascent_int > 0 || ascent_per_approx > 0) )
            ascent_test = 1;

    } else if (ascent_per_approx > 0.0) {
      ascent_test = test_per_approx(timer.time(), timer.dt(), ascent_per_approx);
    }

    if ( (ascent_test == 1) || ( ( ascent_int > 0) && ( timer.nstep() %  ascent_int == 0 ) ) )
    {
        WriteAscentFile(timer.nstep(), timer.time());
    }
#endif


/*--------------------------------------------------------------------------------------------------
 *
 *                               AMReX checkpoint file output control
 *
 *------------------------------------------------------------------------------------------------*/

    if (checkpoint_files_output) {

      int check_test = 0;

      if ( check_int > 0 /* || check_per_exact > 0*/ ||
           check_per_approx > 0. || check_walltime_interval > 0.) {

        // We automatically write checkpoint files with the initial data
        if ( first ) {
          check_test = (restart_file.empty()) ? 1 : 0;
        }
        // We automatically write checkpoint files with the final data
        else if (last) {
          check_test = (timer.nstep() != last_chk) ? 1 : 0;
        }
        else if (check_per_approx > 0) {
          check_test = test_per_approx(timer.time(), timer.dt(), check_per_approx);
        }/*
        else if (check_per_exact > 0 &&
                 (Math::abs(remainder(timer.time(), check_per_exact)) < 1.e-12)) {
          check_test = 1;
        }*/
        else if (check_int > 0) {
          check_test = (timer.nstep() % check_int == 0) ? 1 : 0;
        }
        else if (check_walltime_interval > 0) {
          if (ParallelDescriptor::IOProcessor()) {
            check_test = test_walltime_interval(timer);
          }
          ParallelDescriptor::Bcast(&check_test, 1, ParallelDescriptor::IOProcessorNumber());
        }

        if (check_test == 1) {

          Real time_start = timer.system_time();
          WriteCheckPointFile(check_file, timer.nstep(), timer.dt(), timer.time());
          Real chkpt_write_time = timer.elapsed_runtime(time_start);

          if (timer.walltime_limit() > 0) {
            timer.update_max_write_chkpt_time(chkpt_write_time);
          }

          last_chk = timer.nstep();
        }

        if ( timer.walltime_limit() > 0. && check_test == 0 ) {
          int sufficient_time(1);
          if (ParallelDescriptor::IOProcessor()) {
            sufficient_time = timer.runtime_left_is_sufficient();
          }
          ParallelDescriptor::Bcast(&sufficient_time, 1, ParallelDescriptor::IOProcessorNumber());

          if ( !sufficient_time ) {
            WriteCheckPointFile(check_file, timer.nstep(), timer.dt(), timer.time());
          }
        }
      }

    }

/*--------------------------------------------------------------------------------------------------
 *
 *                               AMReX particle ASCII output control
 *
 *------------------------------------------------------------------------------------------------*/
    int par_ascii_test = 0;

    if ( par_ascii_int > 0) {
      if ( first || last ) {
        par_ascii_test = 1;
      } else if ( timer.nstep() %  par_ascii_int == 0 ) {
        par_ascii_test = 1;
      }

    } else if (par_ascii_per_approx > 0.0) {
      par_ascii_test = test_per_approx(timer.time(), timer.dt(),
          par_ascii_per_approx);

    }

    if( par_ascii_test == 1) {
      WriteParticleAscii(par_ascii_file, timer.nstep());
      last_par_ascii = timer.nstep();
    }


/*--------------------------------------------------------------------------------------------------
 *
 *                                  MFIX mass balance output control
 *
 *------------------------------------------------------------------------------------------------*/
    int mass_balance_report_test = 0;

    if ( mass_balance_report_int > 0 ) {

      if ( first ) { // Never write the first.
        mass_balance_report_test = 0;

      } else if (last) { // Always write the last.
        mass_balance_report_test = (timer.nstep() != last_mb_report) ? 1 :0;

      } else {
        mass_balance_report_test = (timer.nstep() % mass_balance_report_int == 0) ? 1 : 0;
      }
    }

    else if (mass_balance_report_per_approx > 0.0) {
      mass_balance_report_test = test_per_approx(timer.time(),
          timer.dt(), mass_balance_report_per_approx);
    }

    if ( mass_balance_report_test == 1) {
      WriteMassBalanceReport(timer.time());
      last_mb_report = timer.nstep();
    }

    m_p_reports_diameter->write_now(timer, timer.dt(), first, last, pc);

/*--------------------------------------------------------------------------------------------------
 *
 *                                  MFIX call to catalyst
 *
 *------------------------------------------------------------------------------------------------*/

#ifdef AMREX_USE_CATALYST
    if (catalyst_enabled) {
      RunCatalystAdaptor(timer.nstep(), timer.time());
    }
#endif

/*--------------------------------------------------------------------------------------------------
 *
 *                                  MFIX send to MPMD
 *
 *------------------------------------------------------------------------------------------------*/

#ifdef MFIX_MPMD
    bool enable_mpmd = ((m_mpmd_int > 0) || (m_mpmd_per_approx > 0)) ? true : false;

    bool test_mpmd = false;
    if ((m_mpmd_int > 0) && (timer.nstep() % m_mpmd_int == 0)) {
      test_mpmd = true;
    } else if ((m_mpmd_per_approx > 0) &&
        test_per_approx(timer.time(), timer.dt(), m_mpmd_per_approx)) {
      test_mpmd = true;
    }

    if (first && enable_mpmd) {
      MPMDInit();
    } else if (test_mpmd && (!last)) {
      MPMDSend(timer.time());
    } else if (last && enable_mpmd) {
      MPMDFinal();
    }
#endif
}

void MFIXReadWrite::writeStaticPlotFiles() const
{
   if ((m_dem.solve() || m_pic.solve()) && write_ls) {
      WriteStaticPlotFileParticleLevelSet(static_plt_file_ls);
   }

   if (plt_geom) {
      WriteStaticPlotFileEBGeometry(static_plt_file_geom);
   }
}

void MFIXReadWrite::reportGridStats() const
{
   if (fluid.solve())
     ReportGridStats();
}

//
// Print the maximum values of the velocity components
//
void
MFIXReadWrite::mfix_print_max_vel (int lev,
                                   const Vector<MultiFab*>& vel_g_in,
                                   const Vector<MultiFab*>& p_g_in)
{
    amrex::Print() << "   max(abs(u/v/w/p))  = "
                   << vel_g_in[lev]->norm0(0,0,false,true) << "  "
                   << vel_g_in[lev]->norm0(1,0,false,true) << "  "
                   << vel_g_in[lev]->norm0(2,0,false,true) << "  "
                   << p_g_in[lev]->norm0(0,0,false,true) << std::endl;
}

//
// Print the maximum values of the pressure gradient components
//
void
MFIXReadWrite::mfix_print_max_gp (int lev,
                                  const Vector<MultiFab*>& gp_g_in)
{
    amrex::Print() << "   max(abs(gpx/gpy/gpz))  = "
                   << gp_g_in[lev]->norm0(0,0,false,true) << "  "
                   << gp_g_in[lev]->norm0(1,0,false,true) << "  "
                   << gp_g_in[lev]->norm0(2,0,false,true) <<  std::endl;
}


//
// Determine if it is time to write based on approximate interval
//
int
MFIXReadWrite::test_per_approx (const Real time,
                                const Real dt,
                                const Real per_approx)
{
  // Check to see if we've crossed a _per_approx interval by comparing
  // the number of intervals that have elapsed for both the current
  // time and the time at the beginning of this timestep.

  int num_per_old = static_cast<int>( (time-dt) / per_approx );
  int num_per_new = static_cast<int>( (time   ) / per_approx );

  // Before using these, however, we must test for the case where we're
  // within machine epsilon of the next interval. In that case, increment
  // the counter, because we have indeed reached the next par_ascii_per_approx interval
  // at this point.

  const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
  const Real next_time = (num_per_old + 1) * per_approx;

  if ((num_per_new == num_per_old) && amrex::Math::abs(time - next_time) <= eps) {
      num_per_new += 1;
  }

  // Similarly, we have to account for the case where the old time is within
  // machine epsilon of the beginning of this interval, so that we don't double
  // count that time threshold -- we already plotted at that time on the last timestep.

  if ((num_per_new != num_per_old) && amrex::Math::abs((time - dt) - next_time) <= eps)
      num_per_old += 1;

  return (num_per_old != num_per_new) ? 1 : 0;
}


//
// Determine if it is time to write before job is killed
//
int
MFIXReadWrite::test_walltime_interval (const MFIXTimer& timer)
{
  const Real walltime = m_interval_nb * check_walltime_interval;

  if (timer.elapsed_runtime() > walltime) {
    m_interval_nb++;
    return 1;
  }

  return 0;
}


//
//
//
void
MFIXReadWrite::ReportGridStats () const
{
  std::vector<long> counts(6,0);

  int lev = 0;

  const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());

  // Count the number of regular cells
  counts[0] = static_cast<int>(amrex::ReduceSum(*volfrac, 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,&dm] (int i, int j, int k) noexcept
      {if(vfrc(i,j,k)==1.0) dm += 1;});

      return dm;
    }));

  // Count the number of covered cells
  counts[1] = static_cast<int>(amrex::ReduceSum( *volfrac, 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,&dm] (int i, int j, int k) noexcept
      {if(vfrc(i,j,k)==0.0) dm += 1;});

      return dm;
    }));

  // Count the number of cut cells
  counts[2] = static_cast<int>(amrex::ReduceSum( *volfrac, 0,
    [=] AMREX_GPU_HOST_DEVICE (Box const & bx,
                               Array4<const Real> const & vfrc) -> int
    {
      int dm = 0;

      amrex::Loop(bx, [vfrc,&dm] (int i, int j, int k) noexcept
      {if(0.0 < vfrc(i,j,k) && vfrc(i,j,k) < 1.0) dm += 1;});

      return dm;
    }));

  int regular(0), covered(0), cut(0);

  const FabArray<EBCellFlagFab>* flags = &(ebfactory[lev]->getMultiEBCellFlagFab());

  for (MFIter mfi(*volfrac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.growntilebox(IntVect(1,1,1));
    FabType t = (*flags)[mfi].getType(bx);

    // Count number of regular grids
    if ( t == FabType::regular ) {
      regular += 1;
    } else if ( t == FabType::covered ) {
      covered += 1;
    } else {
      cut += 1;
    }
  }

  counts[3] = regular;
  counts[4] = covered;
  counts[5] = cut;

  ParallelDescriptor::ReduceLongSum(counts.data(), 6);

  if(ParallelDescriptor::IOProcessor()){
    printf("\n\n****************************************\n");
    printf("  Coverage report:  Grids        Cells\n");
    printf("          regular:  %5ld   %10ld\n", counts[3], counts[0]);
    printf("          covered:  %5ld   %10ld\n", counts[4], counts[1]);
    printf("              cut:  %5ld   %10ld\n", counts[5], counts[2]);
    printf("****************************************\n\n");
  }
}

//
// Print the minimum volume fraction and cell location.
//
IntVect
MFIXReadWrite::mfix_print_min_epg ()
{

  ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
  ReduceData<int, int, int, int> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (int lev = 0; lev < nlev; lev++) {

    const Real tolerance = std::numeric_limits<Real>::epsilon();
    const Real min_epf = leveldata().epf(lev)->min(0);

    for (MFIter mfi(*leveldata().epf(lev),false); mfi.isValid(); ++mfi) {
      Box const& bx = mfi.tilebox();
      Array4<Real const> const& epf = leveldata().epf_const(lev,mfi);

      reduce_op.eval(bx, reduce_data, [epf, min_epf, tolerance]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        int found(0);
        int iloc(0);
        int jloc(0);
        int kloc(0);
        if( amrex::Math::abs(epf(i,j,k) - min_epf) < tolerance ){
          iloc = i;
          jloc = j;
          kloc = k;
          found = 1;
        }
        return {found, iloc, jloc, kloc};
      });

      ReduceTuple htuple = reduce_data.value();

      const int found(amrex::get<0>(htuple));
      if(found > 0){

        IntVect epf_cell = {amrex::get<1>(htuple),
                            amrex::get<2>(htuple),
                            amrex::get<3>(htuple)};

        amrex::Print(Print::AllProcs)
          << std::endl << std::endl << "min epf "  << min_epf
          << "  at " << epf_cell[0] << "  " << epf_cell[1] << "  " << epf_cell[2]
          << "   total found " << found << std::endl << std::endl;

        return epf_cell;

      }

      //AMREX_ALWAYS_ASSERT(min_epg > 0.275);

    } // mfi
  } // lev

  IntVect fake = {0,0,0};
  return fake;
}


void
MFIXReadWrite::WriteMassBalanceReport (const Real new_time)
{

  if (!report_mass_balance) {
    return;
  }

  // Compute current mass in system
  m_mass_balance->ComputeMassAccum(report_mass_balance, 1);

  long const total_np = (!m_dem.solve() && !m_pic.solve()) ? 0 :
      pc->TotalNumberOfParticles(/*only valid=*/true,/*local=*/false);

  const int nspecies_g = fluid.nspecies();
  const int nspecies_s = solids.nspecies();

  if(ParallelDescriptor::IOProcessor()) {
    printf("\n**********");
    for (int col=0; col < 6; ++col) {
      printf("**************");
    }
    printf("****************\n");

    printf("  Species mass balance for interval:  %12.6f  to %12.6f\n",
           mass_balance_report_time, new_time);

    mass_balance_report_time = new_time;

    printf("\n  %-8s%14s%14s%14s%14s%14s%14s%14s\n", "Species",
           "mass(t+dt)", "mass(t) ", "production",
           "mass(in)", "mass(out)","net accu","flux   ");

    amrex::Array<amrex::Real,7> totals{0.};

    if (fluid.solve()) {

      auto& mass_accum = m_mass_balance->fluid_mass_accum();
      auto& mass_prod = m_mass_balance->fluid_mass_prod();
      auto& mass_inflow = m_mass_balance->fluid_mass_inflow();
      auto& mass_outflow = m_mass_balance->fluid_mass_outflow();

      for (int n=0; n < nspecies_g; ++n) {
        Real delta_accum(mass_accum[n+nspecies_g] - mass_accum[n] - mass_prod[n]);
        Real bc_flux(mass_inflow[n] - mass_outflow[n]);

        //net_acc  += delta_accum;
        //net_flux += bc_flux;

        printf("  %-8s%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e\n",
            fluid.species_names(n).c_str(),
            mass_accum[n+nspecies_g], mass_accum[n], mass_prod[n],
            mass_inflow[n], mass_outflow[n], delta_accum, bc_flux);

        totals[0] += mass_accum[n+nspecies_g];
        totals[1] += mass_accum[n];
        totals[2] += mass_prod[n];
        totals[3] += mass_inflow[n];
        totals[4] += mass_outflow[n];
        totals[5] += delta_accum;
        totals[6] += bc_flux;

        mass_accum[n] = mass_accum[n+nspecies_g];
        mass_inflow[n] = 0.;
        mass_outflow[n] = 0.;
        mass_prod[n] = 0.;
      }
    }

    if (m_dem.solve() || m_pic.solve()) {

      auto& mass_accum = m_mass_balance->solids_mass_accum();
      auto& mass_prod = m_mass_balance->solids_mass_prod();
      auto& mass_inflow = m_mass_balance->solids_mass_inflow();
      auto& mass_outflow = m_mass_balance->solids_mass_outflow();

      for (int n=0; n < nspecies_s; ++n) {
        Real delta_accum(mass_accum[n+nspecies_s] - mass_accum[n] - mass_prod[n]);
        Real bc_flux(mass_inflow[n] - mass_outflow[n]);

        //net_acc += delta_accum;
        //net_flux += bc_flux;

        printf("  %-8s%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e\n",
            solids.species_names(n).c_str(),
            mass_accum[n+nspecies_s], mass_accum[n], mass_prod[n],
            mass_inflow[n], mass_outflow[n],
            delta_accum, bc_flux);

        totals[0] += mass_accum[n+nspecies_s];
        totals[1] += mass_accum[n];
        totals[2] += mass_prod[n];
        totals[3] += mass_inflow[n];
        totals[4] += mass_outflow[n];
        totals[5] += delta_accum;
        totals[6] += bc_flux;

        m_mass_balance->ResetMassBalance(n);
      }
    }

    printf("----------");
    for (int col=0; col < 6; ++col) {
      printf("--------------");
    }
    printf("----------------\n");

    printf("  %-8s%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e\n", "Total",
           totals[0], totals[1], totals[2], totals[3], totals[4], totals[5], totals[6]);

    //printf("\n  total accu := %14.4e %77s\n", net_acc,"sum(mass(t+dt) - mass(t) - production)");
    //printf("    net flux := %14.4e %77s\n", net_flux,"sum(mass(in) - mass(out))");
    //printf("  difference := %14.4e %77s\n", totals[5]-totals[6],"(total acc) - (net flux))");

    if (m_dem.solve() || m_pic.solve()) {
      if (m_dem.solve() ) { printf("\n  total number of particles = "); }
      if (m_pic.solve() ) { printf("\n  total number of parcels =   "); }
      printf("%16s\n", reporter::FormatWithCommas(total_np).c_str() ) ;
    }

    Real diff = std::abs(totals[5]-totals[6]);
    printf("  |(total acc) - (net flux)| =  %14.4e\n", diff);

    printf("**********");
    for (int col=0; col < 6; ++col) {
      printf("**************");
    }
    printf("****************\n\n");
  }
}
