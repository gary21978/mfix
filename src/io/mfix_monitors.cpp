#include <mfix_monitors.H>
#include <mfix_interp_K.H>
#include <mfix_eb_interp_shepard_K.H>
#include <mfix.H>

#include <AMReX_FabArrayUtility.H>
#include <AMReX_Loop.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>


using namespace amrex;
using MFIXParIter = MFIXParticleContainer::MFIXParIter;


Monitors::~Monitors ()
{
  for (int i(0); i < m_monitors.size(); ++i) {
    delete m_monitors[i];
    m_monitors[i] = nullptr;
  }
}


void
Monitors::initialize (MFIXRegions const& regions,
                      MFIXLevelData& a_leveldata,
                      amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& ebfactory,
                      MFIXFluidPhase& fluid,
                      MFIXParticleContainer* pc,
                      amrex::Vector<std::shared_ptr<amrex::EBFArrayBoxFactory>> const& particle_ebfactory,
                      MFIXSolidsPhase& solids)
{
  ParmParse ppRoot("mfix");
  ParmParse ppMonitors("mfix.monitors");

  Vector<std::string> monitor_names;
  ppRoot.queryarr("monitors", monitor_names);

  for (const auto& monitor_name: monitor_names) {

    std::string monitor_type;
    ppMonitors.get(monitor_name, monitor_type);

    std::replace(monitor_type.begin(), monitor_type.end(), ':', ' ');

    std::istringstream iss(monitor_type);
    std::vector<std::string> monitor_specs((std::istream_iterator<std::string>(iss)),
                                           std::istream_iterator<std::string>());

    if (monitor_specs.size() != 3) {
      Print() << "Monitor " << monitor_name << " specs are ill-defined.\n";
      amrex::Abort("Inputs error");
    }

    std::transform(monitor_specs.begin(), monitor_specs.end(),
                   monitor_specs.begin(), amrex::toLower);
    std::string monitor_prefix = "mfix.monitors." + monitor_name;

    std::array<std::string, 2> specs = {monitor_specs[1], monitor_specs[2]};

    std::string mspec0(toLower(monitor_specs[0]));
    std::string mspec1(toLower(monitor_specs[1]));

    if (mspec0 == "eulerian") {

      if (mspec1 ==  "uniformscalarfield") {

        m_monitors.push_back(new EulerianMonitor::UniformScalarField(specs,
              monitor_prefix, a_leveldata, ebfactory, fluid, regions));

      } else if (mspec1 == "pointregion") {

        m_monitors.push_back(new EulerianMonitor::PointRegion(specs,
              monitor_prefix, a_leveldata, ebfactory, fluid, regions));

      } else if (mspec1 ==  "arearegion") {

        m_monitors.push_back(new EulerianMonitor::AreaRegion(specs,
              monitor_prefix, a_leveldata, ebfactory, fluid, regions));

      } else if (mspec1 == "volumeregion") {

        m_monitors.push_back(new EulerianMonitor::VolumeRegion(specs,
              monitor_prefix, a_leveldata, ebfactory, fluid, regions));

      } else if (mspec1 == "surfaceintegral") {

        m_monitors.push_back(new EulerianMonitor::SurfaceIntegral(specs,
              monitor_prefix, a_leveldata, ebfactory, fluid, regions));

      } else if (mspec1 ==  "volumeintegral") {

        m_monitors.push_back(new EulerianMonitor::VolumeIntegral(specs,
              monitor_prefix, a_leveldata, ebfactory, fluid, regions));

      } else {

        amrex::Abort("Unknown Eulerian monitor type " + monitor_specs[1]);
      }

    } else if (mspec0 == "lagrangian") {

      if (mspec1 == "generalproperty") {

        m_monitors.push_back(new LagrangianMonitor::GeneralProperty(specs,
              monitor_prefix, pc, particle_ebfactory, solids, regions));

      } else if (mspec1 == "averagedproperty") {

        m_monitors.push_back(new LagrangianMonitor::AveragedProperty(specs,
              monitor_prefix, pc, particle_ebfactory, solids, regions));

      } else if (mspec1 == "flowrate") {

        m_monitors.push_back(new LagrangianMonitor::FlowRate(specs,
              monitor_prefix, pc, particle_ebfactory, solids, regions));

      } else {

        amrex::Abort("Unknown Lagrangian monitor type " + monitor_specs[1]);
      }

    } else {

      amrex::Abort("Unknown Monitor type " + monitor_specs[0]);
    }
  }
}


void
Monitors::set_pc (MFIXParticleContainer* pc)
{
  for (Monitor* monitor: m_monitors)
    monitor->set_pc(pc);
}


Monitor::Monitor (Vector<std::shared_ptr<EBFArrayBoxFactory>> const& ebfactory,
                  RealBox const& region,
                  RealBox const* plane)
  : m_nlev(0)
  , m_monitoring_results(0)
  , m_ebfactory(ebfactory)
  , m_region(&region)
  , m_boxes(0)
  , m_plane(plane)
{
  m_nlev = m_ebfactory.size();

  AMREX_ASSERT(m_nlev > 0);

  for (int lev(0); lev < m_nlev; ++lev) {

    Box box = calc_box(m_ebfactory[lev]->Geom(), *m_region);
    AMREX_ASSERT(box.ok());

    m_boxes.push_back(box);
  }
}


Monitor::Monitor (const std::array<std::string,2> specs,
                  const std::string monitor_prefix,
                  Vector<std::shared_ptr<EBFArrayBoxFactory>> const& ebfactory,
                  MFIXRegions const& regions)
  : is_initialized(0)
  , m_nlev(0)
  , m_specs(specs)
  , m_monitoring_results(0)
  , m_count_flag(0)
  , m_counts(0)
  , m_ebfactory(ebfactory)
  , m_regions(&regions)
  , m_filename(std::string())
  , m_plot_int(-1)
  , m_plot_per_approx(0.)
  , m_openmode("app")
  , m_setw(0)
  , m_setfill("")
  , m_setprecision(0)
  , m_formatflag("")
  , m_region_name(std::string())
  , m_region(nullptr)
  , m_variables(0)
  , m_input_variables(0)
  , m_boxes(0)
  , m_direction(-1)
  , m_plane_name(std::string())
  , m_plane(nullptr)
{
  ParmParse ppMonitor(monitor_prefix);

  std::string filename;
  ppMonitor.get("plot_file", filename);
  std::filesystem::path full_path_filename(filename);

  m_inputs_full_path = full_path_filename.parent_path();
  m_filename = full_path_filename.filename().string();
  m_filename.append(".csv");

  ppMonitor.getarr("variables", m_input_variables);

  const int plot_int = ppMonitor.query("plot_int", m_plot_int);
  const int per_approx = ppMonitor.query("plot_per_approx", m_plot_per_approx);

  // XOR operation
  AMREX_ALWAYS_ASSERT(plot_int ^ per_approx);

  // Get count flag
  ppMonitor.query("count", m_count_flag);

  ppMonitor.query("output.openmode", m_openmode);

  ppMonitor.query("output.setw", m_setw);
  ppMonitor.query("output.setfill", m_setfill);
  ppMonitor.query("output.setprecision", m_setprecision);
  ppMonitor.query("output.format", m_formatflag);

  if (!(m_setfill.empty() || m_setfill.length() == 1)) {
    amrex::Abort("setfill must be either empty or a single character");
  }

  if (!(m_formatflag.empty() ||
        toLower(m_formatflag).compare("defaultfloat") == 0 ||
        toLower(m_formatflag).compare("fixed") == 0 ||
        toLower(m_formatflag).compare("scientific") == 0)) {
    amrex::Abort("format must be either empty or fixed/scientific");
  }

  ppMonitor.get("region", m_region_name);

  if (m_specs[0].compare("flowrate") == 0) {
    ppMonitor.get("plane", m_plane_name);
  }
}


void
Monitor::initialize ()
{
  is_initialized = 1;

  m_nlev = m_ebfactory.size();

  AMREX_ASSERT(m_nlev > 0);

  AMREX_ASSERT(m_regions->isInitialized());
  m_region = m_regions->getRegion(m_region_name);

  AMREX_ASSERT_WITH_MESSAGE(m_region != nullptr, "Region does not exist");

  AMREX_ASSERT(m_region->ok());

  for (int lev(0); lev < m_nlev; ++lev) {

    Box box = calc_box(m_ebfactory[lev]->Geom(), *m_region);
    AMREX_ASSERT(box.ok());

    m_boxes.push_back(box);
  }

  // Setup variables names
  this->setup_variables();
}


std::ostream&
Monitor::setformat(std::ostream& output) const
{
  std::ostream& result(output);

  if (m_setw > 0)
    result << std::setw(m_setw);

  if (!m_setfill.empty())
    result << std::setfill(m_setfill[0]);

  if (m_setprecision > 0)
    result << std::setprecision(m_setprecision);

  if (toLower(m_formatflag).compare("fixed") == 0)
    result << std::fixed;

  if (toLower(m_formatflag).compare("scientific") == 0)
    result << std::scientific;

  return result;
}


void
Monitor::write_csv (const Real& time,
                    const Real& dt)
// NOTE: right now we only print out lev 0
{
  if (!is_initialized) {
    this->initialize();


    //
//    m_full_paths.resize(m_nlev);
    m_full_paths.resize(1);

//    if (m_nlev > 1) {
//      for (int lev(0); lev < m_nlev; lev++) {
//        std::stringstream lev_ss;
//        lev_ss << std::setw(2) << std::setfill('0') << lev;
//        m_full_paths[lev] = "monitors_results" / m_inputs_full_path / ("level_"+lev_ss.str());
//      }
//    } else {
//      m_full_paths[0] = "monitors_output" / m_inputs_full_path;
      m_full_paths[0] = m_inputs_full_path;
//    }


    // Write csv output file headers
    if (ParallelDescriptor::IOProcessor()) {

//      for (int lev(0); lev < m_nlev; ++lev) {
      for (int lev(0); lev < 1; ++lev) {
        if (!m_full_paths[lev].empty()) {
          if (!std::filesystem::exists(m_full_paths[lev])) {
            std::filesystem::create_directories(m_full_paths[lev]);
          }
        }
      }

//      for (int lev(0); lev < m_nlev; ++lev) {
      for (int lev(0); lev < 1; ++lev) {

        std::string filename = m_full_paths[lev].string() + m_filename;

        std::ofstream output_file;

        if (toLower(m_openmode).compare("app") == 0) {
          output_file.open(filename, std::ios::out | std::ios::app);
        } else if (toLower(m_openmode).compare("trunc") == 0) {
          output_file.open(filename, std::ios::out | std::ios::trunc);
        } else {
          amrex::Abort("Unknown openmode for monitor output file");
        }

        AMREX_ALWAYS_ASSERT(output_file.good());

        output_file << "time";

        if (m_count_flag) {
          output_file << "," << "N";
        }

        for (const std::string& variable: m_variables)
          output_file << "," << variable;
        output_file << std::endl;

        output_file.close();
      }
    }
  }


  // Do the Monitoring
  this->monitor(dt);

  // Count
  if (m_count_flag) {
    m_counts.clear();
    m_counts.resize(m_nlev, 0);
    this->count();
  }


  // MPI synchronization
  ParallelDescriptor::Barrier();


  // Write csv output file monitored values
  if (ParallelDescriptor::IOProcessor()) {

//    for (int lev(0); lev < m_nlev; ++lev) {
    for (int lev(0); lev < 1; ++lev) {

      std::string filename = m_full_paths[lev].string() + m_filename;

      std::ofstream output_file;
      output_file.open(filename, std::ios::out | std::ios::app);

      setformat(output_file) << time;

      if (m_count_flag) {
        auto output_value = m_counts[lev];
        output_file << ",";
        output_file.unsetf(std::ios::floatfield);
        output_file << static_cast<int>(output_value);
      }

      for (int var(0); var < m_monitoring_results[lev].size(); ++var) {

        auto output_value = m_monitoring_results[lev][var];
        output_value = std::abs(output_value) < FLT_MIN ? 0 : output_value;

        output_file << ",";
        setformat(output_file) << output_value;
      }
      output_file << std::endl;

      output_file.close();
    }
  }
}


Box
Monitor::calc_box (const Geometry& geometry,
                   const RealBox& realbox) const
{
  const GpuArray<Real,3> dxi = geometry.InvCellSizeArray();

  const GpuArray<Real,3> prob_lo = geometry.ProbLoArray();
  const GpuArray<Real,3> prob_hi = geometry.ProbHiArray();

  const Real* realbox_lo = realbox.lo();
  const Real* realbox_hi = realbox.hi();

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
    AMREX_ALWAYS_ASSERT(!(realbox_lo[dir] < prob_lo[dir]));
    AMREX_ALWAYS_ASSERT(!(realbox_hi[dir] > prob_hi[dir]));
  }

  const IntVect box_type(AMREX_D_DECL(
        static_cast<int>(std::abs(realbox_hi[0] - realbox_lo[0]) < 1.e-15),
        static_cast<int>(std::abs(realbox_hi[1] - realbox_lo[1]) < 1.e-15),
        static_cast<int>(std::abs(realbox_hi[2] - realbox_lo[2]) < 1.e-15)));

  const int sum = AMREX_D_TERM(box_type[0], + box_type[1], + box_type[2]);
  const int prod = AMREX_D_TERM(box_type[0], * box_type[1], * box_type[2]);

  AMREX_ALWAYS_ASSERT(((sum <= 1) && (prod == 0)) ||
      (amrex::toLower(m_specs[0]).compare("pointregion") == 0 && sum == 3));

  IntVect box_lo, box_hi;

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {

    if (box_type[dir] == 0) {

      box_lo[dir] = static_cast<int>(Math::ceil((realbox_lo[dir]-prob_lo[dir])*dxi[dir]));
      box_hi[dir] = static_cast<int>(Math::floor((realbox_hi[dir]-prob_lo[dir])*dxi[dir]));

    } else if (box_type[dir] == 1) {

      box_lo[dir] = static_cast<int>(std::round((realbox_lo[dir]-prob_lo[dir])*dxi[dir]));
      box_hi[dir] = box_lo[dir];

    } else {
      amrex::Abort("How did we arrive here?");
    }
  }

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
    AMREX_ASSERT(box_lo[dir] <= box_hi[dir]);
  }

  return Box(box_lo, box_hi, box_type);
}


RealBox
Monitor::calc_realbox (const Geometry& geometry,
                       const Box& box) const
{
  const GpuArray<Real,3> dx = geometry.CellSizeArray();

  const GpuArray<Real,3> prob_lo = geometry.ProbLoArray();

  const IntVect box_lo = box.smallEnd();
  const IntVect box_hi = box.bigEnd();

  RealBox realbox;

  realbox.setLo(Vector<Real>{AMREX_D_DECL(prob_lo[0]+box_lo[0]*dx[0],
                                          prob_lo[1]+box_lo[1]*dx[1],
                                          prob_lo[2]+box_lo[2]*dx[2])});

  realbox.setHi(Vector<Real>{AMREX_D_DECL(prob_lo[0]+(box_hi[0]+1)*dx[0],
                                          prob_lo[1]+(box_hi[1]+1)*dx[1],
                                          prob_lo[2]+(box_hi[2]+1)*dx[2])});

#ifdef AMREX_DEBUG
  const GpuArray<Real,3> prob_hi = geometry.ProbHiArray();

  const Real* realbox_lo = realbox.lo();
  const Real* realbox_hi = realbox.hi();

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
    AMREX_ASSERT(!(realbox_lo[dir] < (prob_lo[dir]-std::numeric_limits<Real>::epsilon())));
    AMREX_ASSERT(!(realbox_hi[dir] > (prob_hi[dir]+std::numeric_limits<Real>::epsilon())));
  }
#endif

  return realbox;
}


RealBox
Monitor::realboxes_intersection (const RealBox& realbox_a,
                                 const RealBox& realbox_b) const
{
  RealBox realbox;

  const Real* a_lo = realbox_a.lo();
  const Real* b_lo = realbox_b.lo();

  realbox.setLo(Vector<Real>{AMREX_D_DECL(amrex::max(a_lo[0], b_lo[0]),
                                          amrex::max(a_lo[1], b_lo[1]),
                                          amrex::max(a_lo[2], b_lo[2]))});

  const Real* a_hi = realbox_a.hi();
  const Real* b_hi = realbox_b.hi();

  realbox.setHi(Vector<Real>{AMREX_D_DECL(amrex::min(a_hi[0], b_hi[0]),
                                          amrex::min(a_hi[1], b_hi[1]),
                                          amrex::min(a_hi[2], b_hi[2]))});

  return realbox;
}


namespace EulerianMonitor {

void
BaseMonitor::setup_variables ()
{
  m_mf.clear();
  m_mf.resize(m_nlev, Vector<const MultiFab*>());

  m_scalars.clear();
  m_scalars.resize(m_nlev, Vector<const Real*>());

  m_components.clear();

  Vector<std::string> variable_names;
  variable_names.clear();

  InterphaseTxfrIndexes txfr_idxs;
  InterphaseChemTxfrIndexes chem_txfr_idxs(leveldata(0)->fluid().nspecies(),
                                           leveldata(0)->reactions().solve());

  for (int i(0); i < m_input_variables.size(); ++i) {

    const std::string var = m_input_variables[i];

    if (amrex::toLower(var).compare("ones") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      // resize the vector of MultiFab if it's empty
      if (m_ones.size() != m_nlev)
        m_ones.resize(m_nlev);

      for (int lev(0); lev < m_nlev; lev++) {
        // reference to the epf MultiFab to be used to define Ones MultiFabs
        const MultiFab& epf = *leveldata().epf(lev);

        // define the Ones MultiFab if it was not already defined
        if (!m_ones[lev].isDefined()) {
          // define the Ones MultiFab to be used in monitors
          m_ones[lev].define(epf.boxArray(), epf.DistributionMap(), 1, epf.nGrow(), MFInfo(), epf.Factory());
          // set the value of current Ones MultiFab to 1
          m_ones[lev].setVal(1.);
          EB_set_covered(m_ones[lev], 0, 1, epf.nGrow(), mfix::covered_val);
        }
      }

      // assign Ones MultiFab to the list of multifabs to be monitores
      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(&(m_ones[lev]));

    } else if (var.compare("ep_g") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().epf(lev));

    } else if (var.compare("p_g") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().pert_p(lev));

    } else if (var.compare("thermo_p_g") == 0) {

      reporter::Log(reporter::Warning)
        << "The thermodynamic pressure (thermo_p_g)"
        << " is currently disabled as a monitor variable.";
    } else if (var.compare("ro_g") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().rho(lev));

    } else if (var.compare("trac") == 0) {

      const int ntrac = leveldata().tracer(0)->nComp();

      if (ntrac == 1) {

        variable_names.push_back(var);
        m_components.push_back(0);

        for (int lev(0); lev < m_nlev; lev++)
          m_mf[lev].push_back( leveldata().tracer(lev) );

      } else {

        for (int n(0); n < ntrac; ++n) {

          variable_names.push_back("trac_"+std::to_string(n+1));
          m_components.push_back(n);

          for (int lev(0); lev < m_nlev; lev++)
            m_mf[lev].push_back( leveldata().tracer(lev));

        }
      }

    } else if (var.compare("vel_g") == 0) {

      variable_names.push_back("vel_g_x");
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().vel(lev));

      variable_names.push_back("vel_g_y");
      m_components.push_back(1);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().vel(lev));

      variable_names.push_back("vel_g_z");
      m_components.push_back(2);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().vel(lev));

    } else if (var.compare("vel_g_x") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().vel(lev));

    } else if (var.compare("vel_g_y") == 0) {

      variable_names.push_back(var);
      m_components.push_back(1);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().vel(lev));

    } else if (var.compare("vel_g_z") == 0) {

      variable_names.push_back(var);
      m_components.push_back(2);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().vel(lev));

    } else if (var.compare("gp") == 0) {

      variable_names.push_back("gp_x");
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().grad_p(lev));

      variable_names.push_back("gp_y");
      m_components.push_back(1);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().grad_p(lev));

      variable_names.push_back("gp_z");
      m_components.push_back(2);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().grad_p(lev));

    } else if (var.compare("gp_x") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().grad_p(lev));

    } else if (var.compare("gp_y") == 0) {

      variable_names.push_back(var);
      m_components.push_back(1);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().grad_p(lev));

    } else if (var.compare("gp_z") == 0) {

      variable_names.push_back(var);
      m_components.push_back(2);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().grad_p(lev));

    } else if (var.compare("T_g") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().T(lev));

    } else if (var.compare("h_g") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().h(lev));

    } else if (var.compare("X_gk") == 0) {

      for (int n_g(0); n_g < m_fluid->nspecies(); ++n_g) {

        variable_names.push_back("X_gk_"+m_fluid->species_names(n_g));
        m_components.push_back(n_g);

        for (int lev(0); lev < m_nlev; lev++)
          m_mf[lev].push_back(leveldata().X(lev));

      }

    } else if (var.substr(0,5).compare("X_gk_") == 0) {

      const int var_name_size = var.size();
      const std::string var_species = var.substr(5,var_name_size-5);

      for (int n_g(0); n_g < m_fluid->nspecies(); ++n_g) {
        if (var_species.compare(m_fluid->species_names(n_g)) == 0) {

          variable_names.push_back("X_gk_"+m_fluid->species_names(n_g));
          m_components.push_back(n_g);

          for (int lev(0); lev < m_nlev; lev++)
            m_mf[lev].push_back(leveldata().X(lev));

          break;
        }
      }

    } else if (var.compare("vort") == 0) {

      variable_names.push_back("vort");
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().vorticity(lev));

    } else if (var.compare("txfr_velocity") == 0) {

      variable_names.push_back("txfr_vel_x");
      m_components.push_back(txfr_idxs.vel_src_c+0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

      variable_names.push_back("txfr_vel_y");
      m_components.push_back(txfr_idxs.vel_src_c+1);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

      variable_names.push_back("txfr_vel_z");
      m_components.push_back(txfr_idxs.vel_src_c+2);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

    } else if (var.compare("txfr_vel_x") == 0) {

      variable_names.push_back(var);
      m_components.push_back(txfr_idxs.vel_src_c+0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

    } else if (var.compare("txfr_vel_y") == 0) {

      variable_names.push_back(var);
      m_components.push_back(txfr_idxs.vel_src_c+1);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

    } else if (var.compare("txfr_vel_z") == 0) {

      variable_names.push_back(var);
      m_components.push_back(txfr_idxs.vel_src_c+2);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

    } else if (var.compare("txfr_beta") == 0) {

      variable_names.push_back(var);
      m_components.push_back(txfr_idxs.drag_coeff);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

    } else if (var.compare("txfr_energy_source") == 0) {

      variable_names.push_back(var);
      m_components.push_back(txfr_idxs.energy_source);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

    } else if (var.compare("txfr_gamma") == 0) {

      variable_names.push_back(var);
      m_components.push_back(txfr_idxs.convection_coeff);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().txfr(lev));

    } else if (var.compare("chem_txfr_X_gk") == 0) {

      for (int n_g(0); n_g < m_fluid->nspecies(); ++n_g) {

        variable_names.push_back("chem_txfr_X_gk_"+m_fluid->species_names(n_g));
        m_components.push_back(chem_txfr_idxs.chem_ro_gk+n_g);

        for (int lev(0); lev < m_nlev; lev++)
          m_mf[lev].push_back(leveldata().chem_txfr(lev));

      }

    } else if (var.substr(0,15).compare("chem_txfr_X_gk_") == 0) {

      const int var_name_size = var.size();
      const std::string var_species = var.substr(15,var_name_size-15);

      for (int n_g(0); n_g < m_fluid->nspecies(); ++n_g) {
        if (var_species.compare(m_fluid->species_names(n_g)) == 0) {

          variable_names.push_back("chem_txfr_X_gk_"+m_fluid->species_names(n_g));
          m_components.push_back(chem_txfr_idxs.chem_ro_gk+n_g);

          for (int lev(0); lev < m_nlev; lev++)
            m_mf[lev].push_back(leveldata().chem_txfr(lev));

          break;
        }
      }

    } else if (var.compare("chem_txfr_h") == 0) {

      variable_names.push_back(var);
      m_components.push_back(chem_txfr_idxs.chem_h);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().chem_txfr(lev));

    } else if (var.compare("mac_phi") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

      for (int lev(0); lev < m_nlev; lev++)
        m_mf[lev].push_back(leveldata().mac_phi(lev));

    } else if (toLower(var).compare("none") == 0) {

      if (m_specs[1].compare("area") == 0) {

        variable_names.push_back("area");
        m_components.push_back(0);

        for (int lev(0); lev < m_nlev; lev++)
          m_mf[lev].push_back(nullptr);

      } else if (m_specs[1].compare("massflowrate") == 0) {

        variable_names.push_back("mass_flow_rate");
        m_components.push_back(0);

        for (int lev(0); lev < m_nlev; lev++)
          m_mf[lev].push_back(nullptr);

      } else if (m_specs[1].compare("volumeflowrate") == 0) {
        variable_names.push_back("volume_flow_rate");
        m_components.push_back(0);

        for (int lev(0); lev < m_nlev; lev++)
          m_mf[lev].push_back(nullptr);

      } else if (m_specs[1].compare("volume") == 0) {

        variable_names.push_back("volume");
        m_components.push_back(0);

        for (int lev(0); lev < m_nlev; lev++)
          m_mf[lev].push_back(nullptr);

      } else {

        amrex::Abort("inputs error");
      }

    } else {

      Print() << "var = " << var << "\n";
      amrex::Abort("Unrecognized variable to monitor. Fix inputs");
    }
  }

  m_variables.swap(variable_names);
}


void
BaseMonitor::count ()
{
  for (int lev(0); lev < m_nlev; lev++) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<int> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& /*ijk*/) -> ReduceTuple
    { return 1; };

    ReduceTuple host_tuple = apply(lev, nullptr, 0, reduce_data, reduce_op, R, default_value);

    int l_count = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceIntSum(&l_count, 1);

    m_counts[lev] = l_count;
  }
}


void
UniformScalarField::monitor (const Real& /*dt*/)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  if (m_specs[1].compare("value") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {

      m_monitoring_results[lev] = value(m_scalars[lev]);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
UniformScalarField::value (const Vector<const Real*>& scalars)
{
  const int var_nb = scalars.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    if (scalars[var] != nullptr)
      result[var] = *(scalars[var]);
  }

  return result;
}


void
PointRegion::set_point_coords ()
{
  AMREX_ASSERT(AMREX_D_TERM(std::abs(m_region->length(0)) < 1.e-15, &&
                            std::abs(m_region->length(1)) < 1.e-15, &&
                            std::abs(m_region->length(2)) < 1.e-15));

  m_point = RealVect(m_region->lo());
}


void
PointRegion::check_boxes_are_ok () const
{
  for (int lev(0); lev < m_nlev; ++lev) {

    AMREX_ASSERT(AMREX_D_TERM(m_boxes[lev].length(0) == 1, &&
                              m_boxes[lev].length(1) == 1, &&
                              m_boxes[lev].length(2) == 1));
  }
}


bool
PointRegion::check_mf_is_ok (const MultiFab& mf) const
{
  const BoxArray& box_array = mf.boxArray();
  const IndexType& idx_type = box_array.ixType();

  if (!(idx_type.cellCentered() || idx_type.nodeCentered())) {
    return false;
  }

  return true;
}


void
PointRegion::monitor (const Real& /*dt*/)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  if (m_specs[1].compare("value") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {

      m_monitoring_results[lev] = value(lev, m_mf[lev], m_components);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
PointRegion::value (const int lev,
                    const Vector<const MultiFab*>& mf,
                    const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<const MultiFab*> MFs(mf);
  Vector<int> ConversionFlags(var_nb, 0);

  // TODO TODO TODO
  // Be sure if we need to reset mf to its original value after we have deleted
  // the converted MF
//  convert_mf_if_needed(MFs, ConversionFlags, lev);

  for (int i(0); i < MFs.size(); ++i) {
    if (MFs[i] != nullptr) {
      AMREX_ASSERT(check_mf_is_ok(*MFs[i]));
    }
  }

  Vector<Real> result(var_nb, 0.);

  const Vector<Real> default_value(var_nb, std::numeric_limits<Real>::min());

  const auto& geom = m_ebfactory[lev]->Geom();

  const auto dxi = geom.InvCellSizeArray();
  const auto dx  = geom.CellSizeArray();
  const auto plo = geom.ProbLoArray();

  const auto& cellcent = m_ebfactory[lev]->getCentroid();
  const auto& bndrycent = m_ebfactory[lev]->getBndryCent();
  const auto& areafrac = m_ebfactory[lev]->getAreaFrac();

  amrex::BoxArray box_array = this->get_box_array(lev);
  const amrex::DistributionMapping& d_map = this->get_d_map(lev);

  const amrex::FabArray<amrex::EBCellFlagFab>& flags = m_ebfactory[lev]->getMultiEBCellFlagFab();

  MultiFab* interp_mf;
  interp_mf = new MultiFab(box_array, d_map, var_nb, 1, MFInfo(), *m_ebfactory[lev]);

  for (int var(0); var < var_nb; ++var) {

    AMREX_ASSERT(MFs[var] != nullptr);

    MultiFab::Copy(*interp_mf, *MFs[var], components[var], var, 1, 1);
  }

  interp_mf->FillBoundary(geom.periodicity());

  for (amrex::MFIter mfi(box_array, d_map, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    const int iloc = static_cast<int>(amrex::Math::floor((m_point[0] - plo[0])*dxi[0]));
    const int jloc = static_cast<int>(amrex::Math::floor((m_point[1] - plo[1])*dxi[1]));
    const int kloc = static_cast<int>(amrex::Math::floor((m_point[2] - plo[2])*dxi[2]));

    if (bx.contains(iloc, jloc, kloc)) {

      const amrex::EBCellFlagFab& flags_fab = flags[mfi];

      if (flags_fab.getType(amrex::grow(bx,0)) != FabType::covered) {

        const int grown_bx_is_regular = (flags_fab.getType(amrex::grow(bx,1)) == FabType::regular);

        const auto& flags_array = flags_fab.array();

        const Array4<const Real> empty_array;

        // mf array
        const auto& interp_mf_array = interp_mf->const_array(mfi);

        // Cell centroids
        const auto& ccent_fab = grown_bx_is_regular ? empty_array : cellcent.const_array(mfi);
        // Centroid of EB
        const auto& bcent_fab = grown_bx_is_regular ? empty_array : bndrycent.const_array(mfi);
        // Area fractions
        const auto& apx_fab = grown_bx_is_regular ? empty_array : areafrac[0]->const_array(mfi);
        const auto& apy_fab = grown_bx_is_regular ? empty_array : areafrac[1]->const_array(mfi);
        const auto& apz_fab = grown_bx_is_regular ? empty_array : areafrac[2]->const_array(mfi);

        if (grown_bx_is_regular) {

          trilinear_interp(m_point, result.dataPtr(), interp_mf_array, plo, dxi, var_nb);

        } else {

          if (flags_array(iloc,jloc,kloc).isCovered()) {

            result = default_value;

            break;

          } else {

            const int i = static_cast<int>(amrex::Math::floor((m_point[0] - plo[0])*dxi[0] + 0.5));
            const int j = static_cast<int>(amrex::Math::floor((m_point[1] - plo[1])*dxi[1] + 0.5));
            const int k = static_cast<int>(amrex::Math::floor((m_point[2] - plo[2])*dxi[2] + 0.5));

            // All cells in the stencil are regular. Use traditional trilinear
            // interpolation
            if (flags_array(i-1,j-1,k-1).isRegular() &&
                flags_array(i  ,j-1,k-1).isRegular() &&
                flags_array(i-1,j  ,k-1).isRegular() &&
                flags_array(i  ,j  ,k-1).isRegular() &&
                flags_array(i-1,j-1,k  ).isRegular() &&
                flags_array(i  ,j-1,k  ).isRegular() &&
                flags_array(i-1,j  ,k  ).isRegular() &&
                flags_array(i  ,j  ,k  ).isRegular()) {

              trilinear_interp(m_point, result.dataPtr(), interp_mf_array, plo, dxi, var_nb);

            } else {

              shepard_interp(m_point, iloc, jloc, kloc, dx, dxi, plo, flags_array,
                             ccent_fab, bcent_fab, apx_fab, apy_fab, apz_fab,
                             interp_mf_array, result.dataPtr(), 0, 0, var_nb);
            }
          }

        } // grown box is not regular
      } // FAB is not covered
    } // this Box contains m_point
  } // MFIter loop

  ParallelDescriptor::ReduceRealMax(result.dataPtr(), var_nb);

  // Free up memory
  delete interp_mf;

  for (int i(0); i < MFs.size(); ++i) {
    if (ConversionFlags[i] == 1) {
      delete MFs[i];
    }
  }

  return result;
}


int
PointRegion::convert_mf_if_needed (const MultiFab*& mf,
                                   const int lev)
{
  if (mf == nullptr)
    return 0;

  int conversion_flag(0);

  BoxArray box_array = this->get_box_array(lev);
  const DistributionMapping& d_map = this->get_d_map(lev);

  if (mf->boxArray().ixType().nodeCentered()) {

    MultiFab* mf_cc = new MultiFab(box_array, d_map, mf->nComp(), 1, MFInfo(), *m_ebfactory[lev]);
    amrex::average_node_to_cellcenter(*mf_cc, 0, *mf, 0, mf->nComp(), 1);
    mf = mf_cc;
    conversion_flag = 1;

  } else if (!mf->boxArray().ixType().cellCentered()) {

    amrex::Abort("Face centered to cell centered conversion not yet implemented here");
  }

  AMREX_ASSERT(check_mf_is_ok(*mf));

  return conversion_flag;
}


void
AreaMonitor::set_direction ()
{
  const IndexType& idx_type = m_boxes[0].ixType();

  for (int dir(0); dir < AMREX_SPACEDIM; ++dir)
    if (idx_type.nodeCentered(dir)) {
      m_direction = dir;
      break;
    }

  AMREX_ASSERT(AMREX_D_TERM(m_boxes[0].length(m_direction) == 1, &&
                            m_boxes[0].length((m_direction+1)%AMREX_SPACEDIM) > 1, &&
                            m_boxes[0].length((m_direction+2)%AMREX_SPACEDIM) > 1));
}


void
AreaMonitor::check_boxes_are_ok () const
{
  for (int lev(0); lev < m_nlev; ++lev) {

    const int min_box_length = amrex::min(AMREX_D_DECL(m_boxes[lev].length(0),
                                                       m_boxes[lev].length(1),
                                                       m_boxes[lev].length(2)));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(min_box_length == 1,
        "AreaMonitor Box is over-dimensionated");
  }
}


bool
AreaMonitor::check_mf_is_ok (const MultiFab& mf) const
{
  const BoxArray& box_array = mf.boxArray();
  const IndexType& idx_type = box_array.ixType();

  if (!idx_type.nodeCentered(m_direction)) {
    return false;
  }

  for (int dir(1); dir < AMREX_SPACEDIM; ++dir) {
    if (!idx_type.cellCentered((m_direction+dir)%AMREX_SPACEDIM)) {
      return false;
    }
  }

  return true;
}


int
AreaMonitor::convert_mf_if_needed (const MultiFab*& mf,
                                   const int lev)
{
  if (mf == nullptr)
    return 0;

  int conversion_flag(0);

  BoxArray box_array = m_ebfactory[lev]->boxArray();
  const DistributionMapping& d_map = this->get_d_map(lev);

  if (mf->boxArray().ixType().cellCentered()) {

    Array<MultiFab*, AMREX_SPACEDIM> mf_fc;

    for(int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      BoxArray edge_ba = box_array;
      edge_ba.surroundingNodes(dir);
      mf_fc[dir] = new MultiFab(edge_ba, d_map, mf->nComp(), 1, MFInfo(), *m_ebfactory[lev]);
      mf_fc[dir]->setVal(0.);
    }

    const auto& geom = m_ebfactory[lev]->Geom();
    Vector<BCRec> bc_rec(mf->nComp(), BCRec());

    EB_interp_CellCentroid_to_FaceCentroid(*mf, mf_fc, 0, 0, mf->nComp(), geom, bc_rec);

    mf = mf_fc[m_direction];
    conversion_flag = 1;

    for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
      if (dir != m_direction)
        delete mf_fc[dir];
    }

  } else if (mf->boxArray().ixType().nodeCentered()) {

    // first convert to cell centered
    MultiFab* mf_cc = new MultiFab(box_array, d_map, mf->nComp(), 1, MFInfo(), *m_ebfactory[lev]);
    amrex::average_node_to_cellcenter(*mf_cc, 0, *mf, 0, mf->nComp(), 1);

    // then interpolate to facecentroid
    Array<MultiFab*, AMREX_SPACEDIM> mf_fc;

    for(int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      BoxArray edge_ba = box_array;
      edge_ba.surroundingNodes(dir);
      mf_fc[dir] = new MultiFab(edge_ba, d_map, mf->nComp(), 1, MFInfo(), *m_ebfactory[lev]);
      mf_fc[dir]->setVal(0.);
    }

    const auto& geom = m_ebfactory[lev]->Geom();
    Vector<BCRec> bc_rec(1, BCRec());

    EB_interp_CellCentroid_to_FaceCentroid(*mf_cc, mf_fc, 0, 0, mf->nComp(), geom, bc_rec);

    mf = mf_fc[m_direction];
    conversion_flag = 1;

    delete mf_cc;

    for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
      if (dir != m_direction)
        delete mf_fc[dir];
    }
  }

  AMREX_ASSERT(check_mf_is_ok(*mf));

  return conversion_flag;
}


void
AreaRegion::monitor (const Real& /*dt*/)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  m_epsilon.clear();
  m_epsilon.resize(m_nlev, nullptr);

  m_density.clear();
  m_density.resize(m_nlev, nullptr);

  m_velocity.clear();
  m_velocity.resize(m_nlev, nullptr);

  // Select the monitor
  if (m_specs[1].compare("sum") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = sum(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("min") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = min(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("max") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = max(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("average") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = average(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("standarddeviation") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = stddev(lev, m_mf[lev], m_components);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
AreaRegion::sum (const int lev,
                 const Vector<const MultiFab*>& mf,
                 const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*areafrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_sum = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_sum, 1);

    result[var] = l_sum;
  }

  return result;
}


Vector<Real>
AreaRegion::min (const int lev,
                 const Vector<const MultiFab*>& mf,
                 const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::max()};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*areafrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_min = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMin(&l_min, 1);

    result[var] = l_min;
  }

  return result;
}


Vector<Real>
AreaRegion::max (const int lev,
                 const Vector<const MultiFab*>& mf,
                 const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::min()};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*areafrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_max = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMax(&l_max, 1);

    result[var] = l_max;
  }

  return result;
}


Vector<Real>
AreaRegion::average (const int lev,
                     const Vector<const MultiFab*>& mf,
                     const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,long> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0., 0};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*areafrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk), 1}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    long denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceLongSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / Real(denominator);
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


Vector<Real>
AreaRegion::stddev (const int lev,
                    const Vector<const MultiFab*>& mf,
                    const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  const Vector<Real> avg = this->average(lev, mf, components);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,long> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0., 0};

    const Real d_avg = avg[var];

    auto R = [d_avg] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                       Array4<Real const> const& /*areafrac_arr*/,
                                       Array4<Real const> const& /*epsilon_arr*/,
                                       Array4<Real const> const& /*density_arr*/,
                                       Array4<Real const> const& /*velocity_arr*/,
                                       const IntVect& ijk) -> ReduceTuple
    {
      Real diff = mf_arr(ijk) - d_avg;
      return {diff*diff, 1};
    };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    long denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceLongSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = std::sqrt(numerator / Real(denominator));
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


void
VolumeMonitor::check_boxes_are_ok () const
{
  for (int lev(0); lev < m_nlev; ++lev) {

    const IndexType& idx_type = m_boxes[lev].ixType();

    if (!idx_type.cellCentered())
      amrex::Abort("Box is not cell-centered");

    for (int dir(0); dir < AMREX_SPACEDIM; ++dir) {
      if (m_boxes[lev].length(dir) < 1)
        amrex::Abort("VolumeMonitor monitor Box is not made of one single cell");
    }
  }
}


bool
VolumeMonitor::check_mf_is_ok (const MultiFab& mf) const
{
  const BoxArray& box_array = mf.boxArray();
  const IndexType& idx_type = box_array.ixType();

  if (!(idx_type.cellCentered() || idx_type.nodeCentered()))
    return false;

  return true;
}


int
VolumeMonitor::convert_mf_if_needed (const MultiFab*& mf,
                                     const int lev)
{
  if (mf == nullptr)
    return 0;

  int conversion_flag = 0;

  BoxArray box_array = this->get_box_array(lev);
  const DistributionMapping& d_map = this->get_d_map(lev);

  if (mf->boxArray().ixType().nodeCentered()) {

    MultiFab* mf_cc = new MultiFab(box_array, d_map, mf->nComp(), 1, MFInfo(), *m_ebfactory[lev]);
    amrex::average_node_to_cellcenter(*mf_cc, 0, *mf, 0, mf->nComp(), 1);
    mf = mf_cc;
    conversion_flag = 1;

  } else if (!mf->boxArray().ixType().cellCentered()) {

    amrex::Abort("Face centered to cell centered conversion not yet implemented here");
  }

  AMREX_ASSERT(check_mf_is_ok(*mf));

  return conversion_flag;
}


void
VolumeRegion::monitor (const Real& /*dt*/)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  m_epsilon.clear();
  m_epsilon.resize(m_nlev, nullptr);

  m_density.clear();
  m_density.resize(m_nlev, nullptr);

  m_velocity.clear();
  m_velocity.resize(m_nlev, nullptr);

  // Select the monitor
  if (m_specs[1].compare("sum") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = sum(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("min") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = min(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("max") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = max(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("average") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = average(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("standarddeviation") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = stddev(lev, m_mf[lev], m_components);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
VolumeRegion::sum (const int lev,
                   const Vector<const MultiFab*>& mf,
                   const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_sum = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_sum, 1);

    result[var] = l_sum;
  }

  return result;
}


Vector<Real>
VolumeRegion::min (const int lev,
                   const Vector<const MultiFab*>& mf,
                   const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::max()};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_min = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMin(&l_min, 1);

    result[var] = l_min;
  }

  return result;
}


Vector<Real>
VolumeRegion::max (const int lev,
                   const Vector<const MultiFab*>& mf,
                   const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::min()};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return mf_arr(ijk); };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_max = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMax(&l_max, 1);

    result[var] = l_max;
  }

  return result;
}


Vector<Real>
VolumeRegion::average (const int lev,
                       const Vector<const MultiFab*>& mf,
                       const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,long> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0., 0};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& /*volfrac_arr*/,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk), 1}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    long denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceLongSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / Real(denominator);
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


Vector<Real>
VolumeRegion::stddev (const int lev,
                      const Vector<const MultiFab*>& mf,
                      const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  const Vector<Real> avg = this->average(lev, mf, components);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,long> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0., 0};

    const Real d_avg = avg[var];

    auto R = [d_avg] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                       Array4<Real const> const& /*volfrac_arr*/,
                                       Array4<Real const> const& /*epsilon_arr*/,
                                       Array4<Real const> const& /*density_arr*/,
                                       Array4<Real const> const& /*velocity_arr*/,
                                       const IntVect& ijk) -> ReduceTuple
    {
      Real diff = mf_arr(ijk) - d_avg;
      return {diff*diff, 1};
    };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    long denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceLongSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = std::sqrt(numerator / Real(denominator));
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


void
SurfaceIntegral::set_dA ()
{
  m_dA.resize(m_nlev);

  for (int lev(0); lev < m_nlev; ++lev) {

    const Real* dx = m_ebfactory[lev]->Geom().CellSize();

    m_dA[lev] = 1.;

    for (int dir(1); dir < AMREX_SPACEDIM; ++dir) {
      const int i = (m_direction + dir) % AMREX_SPACEDIM;
      m_dA[lev] *= dx[i];
    }
  }
}


void
SurfaceIntegral::monitor (const Real& /*dt*/)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  m_epsilon.clear();
  m_epsilon.resize(m_nlev, nullptr);

  m_density.clear();
  m_density.resize(m_nlev, nullptr);

  m_velocity.clear();
  m_velocity.resize(m_nlev, nullptr);

  // Select the monitor
  if (m_specs[1].compare("area") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = area(lev);
    }

  } else if (m_specs[1].compare("areaweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = area_weighted_average(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("flowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_epsilon[lev] = leveldata().epf(lev);
      m_density[lev] = leveldata().rho(lev);
      m_velocity[lev] = leveldata().vel(lev);

      m_monitoring_results[lev] = flow_rate(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("massflowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_epsilon[lev] = leveldata().epf(lev);
      m_density[lev] = leveldata().rho(lev);
      m_velocity[lev] = leveldata().vel(lev);

      m_monitoring_results[lev] = mass_flow_rate(lev);
    }

  } else if (m_specs[1].compare("massweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_epsilon[lev] = leveldata().epf(lev);
      m_density[lev] = leveldata().rho(lev);
      m_velocity[lev] = leveldata().vel(lev);

      m_monitoring_results[lev] = mass_weighted_average(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("volumeflowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_epsilon[lev] = leveldata().epf(lev);
      m_velocity[lev] = leveldata().vel(lev);

      m_monitoring_results[lev] = volume_flow_rate(lev);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
SurfaceIntegral::area (const int lev)
{
  Vector<Real> result(1, 0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);

  using ReduceTuple = typename decltype(reduce_data)::Type;

  ReduceTuple default_value = {0.};

  const Real dA = m_dA[lev];

  auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& areafrac_arr,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
  { return areafrac_arr(ijk)*dA; };

  ReduceTuple host_tuple = apply(lev, nullptr, 0, reduce_data, reduce_op, R, default_value);

  Real l_area = amrex::get<0>(host_tuple);

  ParallelDescriptor::ReduceRealSum(&l_area, 1);

  result[0] = l_area;

  return result;
}


Vector<Real>
SurfaceIntegral::area_weighted_average (const int lev,
                                        const Vector<const MultiFab*>& mf,
                                        const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0., 0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& areafrac_arr,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk)*areafrac_arr(ijk), areafrac_arr(ijk)}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / denominator;
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


Vector<Real>
SurfaceIntegral::flow_rate (const int lev,
                            const Vector<const MultiFab*>& mf,
                            const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    const Real dA = m_dA[lev];

    auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& areafrac_arr,
                                    Array4<Real const> const& epsilon_arr,
                                    Array4<Real const> const& density_arr,
                                    Array4<Real const> const& velocity_arr,
                                    const IntVect& ijk) -> ReduceTuple
    { return {epsilon_arr(ijk)*density_arr(ijk)*mf_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_flow_rate = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_flow_rate, 1);

    result[var] = l_flow_rate;
  }

  return result;
}


Vector<Real>
SurfaceIntegral::mass_flow_rate (const int lev)
{
  Vector<Real> result(1, 0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);

  using ReduceTuple = typename decltype(reduce_data)::Type;

  ReduceTuple default_value = {0.};

  const Real dA = m_dA[lev];

  auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& areafrac_arr,
                                  Array4<Real const> const& epsilon_arr,
                                  Array4<Real const> const& density_arr,
                                  Array4<Real const> const& velocity_arr,
                                  const IntVect& ijk) -> ReduceTuple
  { return {epsilon_arr(ijk)*density_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA}; };

  ReduceTuple host_tuple = apply(lev, nullptr, 0, reduce_data, reduce_op, R, default_value);

  Real l_mass_flow_rate = amrex::get<0>(host_tuple);

  ParallelDescriptor::ReduceRealSum(&l_mass_flow_rate, 1);

  result[0] = l_mass_flow_rate;

  return result;
}


Vector<Real>
SurfaceIntegral::mass_weighted_average (const int lev,
                                        const Vector<const MultiFab*>& mf,
                                        const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0., 0.};

    const Real dA = m_dA[lev];

    auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& areafrac_arr,
                                    Array4<Real const> const& epsilon_arr,
                                    Array4<Real const> const& density_arr,
                                    Array4<Real const> const& velocity_arr,
                                    const IntVect& ijk) -> ReduceTuple
    { return {epsilon_arr(ijk)*density_arr(ijk)*mf_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA,
              epsilon_arr(ijk)*density_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / denominator;
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


Vector<Real>
SurfaceIntegral::volume_flow_rate (const int lev)
{
  Vector<Real> result(1, 0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);

  using ReduceTuple = typename decltype(reduce_data)::Type;

  ReduceTuple default_value = {0.};

  const Real dA = m_dA[lev];

  auto R = [dA] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& areafrac_arr,
                                  Array4<Real const> const& epsilon_arr,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& velocity_arr,
                                  const IntVect& ijk) -> ReduceTuple
  { return {epsilon_arr(ijk)*velocity_arr(ijk)*areafrac_arr(ijk)*dA}; };

  ReduceTuple host_tuple = apply(lev, nullptr, 0, reduce_data, reduce_op, R, default_value);

  Real l_volume_flow_rate = amrex::get<0>(host_tuple);

  ParallelDescriptor::ReduceRealSum(&l_volume_flow_rate, 1);

  result[0] = l_volume_flow_rate;

  return result;
}


void
VolumeIntegral::set_dV ()
{
  m_dV.resize(m_nlev);

  for (int lev(0); lev < m_nlev; ++lev) {

    m_dV[lev] = 1.;

    for (int dir(0); dir < AMREX_SPACEDIM; ++dir)
      m_dV[lev] *= m_ebfactory[lev]->Geom().CellSize(dir);
  }
}


void
VolumeIntegral::monitor (const Real& /*dt*/)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  m_epsilon.clear();
  m_epsilon.resize(m_nlev, nullptr);

  m_density.clear();
  m_density.resize(m_nlev, nullptr);

  m_velocity.clear();
  m_velocity.resize(m_nlev, nullptr);

  // Select the monitor
  if (m_specs[1].compare("volume") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {

      m_monitoring_results[lev] = volume(lev);
    }

  } else if (m_specs[1].compare("volumeintegral") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = volume_integral(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("volumeweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = volume_weighted_average(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("massweightedintegral") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_epsilon[lev] = leveldata().epf(lev);
      m_density[lev] = leveldata().rho(lev);

      m_monitoring_results[lev] = mass_weighted_integral(lev, m_mf[lev], m_components);
    }

  } else if (m_specs[1].compare("massweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_epsilon[lev] = leveldata().epf(lev);
      m_density[lev] = leveldata().rho(lev);

      m_monitoring_results[lev] = mass_weighted_average(lev, m_mf[lev], m_components);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
VolumeIntegral::volume (const int lev)
{
  Vector<Real> result(1, 0.);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);

  using ReduceTuple = typename decltype(reduce_data)::Type;

  ReduceTuple default_value = {0.};

  const Real dV = m_dV[lev];

  auto R = [dV] AMREX_GPU_DEVICE (Array4<Real const> const& /*mf_arr*/,
                                  Array4<Real const> const& volfrac_arr,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
  { return volfrac_arr(ijk)*dV; };

  ReduceTuple host_tuple = apply(lev, nullptr, 0, reduce_data, reduce_op, R, default_value);

  Real l_volume = amrex::get<0>(host_tuple);

  ParallelDescriptor::ReduceRealSum(&l_volume, 1);

  result[0] = l_volume;

  return result;
}


Vector<Real>
VolumeIntegral::volume_integral (const int lev,
                                 const Vector<const MultiFab*>& mf,
                                 const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    const Real dV = m_dV[lev];

    auto R = [dV] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& volfrac_arr,
                                    Array4<Real const> const& /*epsilon_arr*/,
                                    Array4<Real const> const& /*density_arr*/,
                                    Array4<Real const> const& /*velocity_arr*/,
                                    const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk)*volfrac_arr(ijk)*dV}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_volume_integral = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_volume_integral, 1);

    result[var] = l_volume_integral;
  }

  return result;
}


Vector<Real>
VolumeIntegral::volume_weighted_integral (const int lev,
                                          const Vector<const MultiFab*>& mf,
                                          const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    const Real dV = m_dV[lev];

    auto R = [dV] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& volfrac_arr,
                                    Array4<Real const> const& epsilon_arr,
                                    Array4<Real const> const& /*density_arr*/,
                                    Array4<Real const> const& /*velocity_arr*/,
                                    const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk)*epsilon_arr(ijk)*volfrac_arr(ijk)*dV}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_volume_integral = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_volume_integral, 1);

    result[var] = l_volume_integral;
  }

  return result;
}


Vector<Real>
VolumeIntegral::volume_weighted_average (const int lev,
                                         const Vector<const MultiFab*>& mf,
                                         const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0., 0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& volfrac_arr,
                                  Array4<Real const> const& /*epsilon_arr*/,
                                  Array4<Real const> const& /*density_arr*/,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {mf_arr(ijk)*volfrac_arr(ijk),
              volfrac_arr(ijk)}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / denominator;
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


Vector<Real>
VolumeIntegral::mass_weighted_integral (const int lev,
                                        const Vector<const MultiFab*>& mf,
                                        const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0.};

    const Real dV = m_dV[lev];

    auto R = [dV] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                    Array4<Real const> const& volfrac_arr,
                                    Array4<Real const> const& epsilon_arr,
                                    Array4<Real const> const& density_arr,
                                    Array4<Real const> const& /*velocity_arr*/,
                                    const IntVect& ijk) -> ReduceTuple
    { return epsilon_arr(ijk)*density_arr(ijk)*mf_arr(ijk)*volfrac_arr(ijk)*dV; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real l_integral = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_integral, 1);

    result[var] = l_integral;
  }

  return result;
}


Vector<Real>
VolumeIntegral::mass_weighted_average (const int lev,
                                       const Vector<const MultiFab*>& mf,
                                       const Vector<int>& components)
{
  const int var_nb = mf.size();

  AMREX_ASSERT(components.size() == var_nb);

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0., 0.};

    auto R = [] AMREX_GPU_DEVICE (Array4<Real const> const& mf_arr,
                                  Array4<Real const> const& volfrac_arr,
                                  Array4<Real const> const& epsilon_arr,
                                  Array4<Real const> const& density_arr,
                                  Array4<Real const> const& /*velocity_arr*/,
                                  const IntVect& ijk) -> ReduceTuple
    { return {epsilon_arr(ijk)*density_arr(ijk)*mf_arr(ijk)*volfrac_arr(ijk),
              epsilon_arr(ijk)*density_arr(ijk)*volfrac_arr(ijk)}; };

    ReduceTuple host_tuple = apply(lev, mf[var], components[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / denominator;
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}

} // end namespace EulerianMonitor


namespace LagrangianMonitor {

void
BaseMonitor::setup_variables ()
{
  m_components.clear();

  Vector<std::string> variable_names;
  variable_names.clear();

  for (int i(0); i < m_input_variables.size(); ++i) {

    const std::string& var = m_input_variables[i];

    if (amrex::toLower(var).compare("ones") == 0) {

      variable_names.push_back(var);
      m_components.push_back(GetParticleValue::ones);

    } else if (amrex::toLower(var).compare("volume") == 0) {

      variable_names.push_back(var);
      m_components.push_back(GetParticleValue::volume);

    } else if (amrex::toLower(var).compare("mass") == 0) {

      variable_names.push_back(var);
      m_components.push_back(GetParticleValue::mass);

    } else if (amrex::toLower(var).compare("k_energy") == 0) {

      variable_names.push_back(var);
      m_components.push_back(GetParticleValue::k_energy);

    } else if (amrex::toLower(var).compare("position") == 0) {

      variable_names.push_back("pos_x");
      m_components.push_back(0);

      variable_names.push_back("pos_y");
      m_components.push_back(1);

      variable_names.push_back("pos_z");
      m_components.push_back(2);

    } else if (amrex::toLower(var).compare("pos_x") == 0) {

      variable_names.push_back(var);
      m_components.push_back(0);

    } else if (amrex::toLower(var).compare("pos_y") == 0) {

      variable_names.push_back(var);
      m_components.push_back(1);

    } else if (amrex::toLower(var).compare("pos_z") == 0) {

      variable_names.push_back(var);
      m_components.push_back(2);

    } else if (amrex::toLower(var).compare("id") == 0) {

      variable_names.push_back(var);
      m_components.push_back(3);

    } else if (amrex::toLower(var).compare("cpu") == 0) {

      variable_names.push_back(var);
      m_components.push_back(4);

    } else if (amrex::toLower(var).compare("radius") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::radius);

    } else if (amrex::toLower(var).compare("density") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::density);

    } else if (amrex::toLower(var).compare("velocity") == 0) {

      variable_names.push_back("vel_x");
      m_components.push_back(5+SoArealData::velx);

      variable_names.push_back("vel_y");
      m_components.push_back(5+SoArealData::vely);

      variable_names.push_back("vel_z");
      m_components.push_back(5+SoArealData::velz);

    } else if (amrex::toLower(var).compare("vel_x") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::velx);

    } else if (amrex::toLower(var).compare("vel_y") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::vely);

    } else if (amrex::toLower(var).compare("vel_z") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::velz);

    } else if (amrex::toLower(var).compare("omega") == 0) {

      variable_names.push_back("omega_x");
      m_components.push_back(5+SoArealData::omegax);

      variable_names.push_back("omega_y");
      m_components.push_back(5+SoArealData::omegay);

      variable_names.push_back("omega_z");
      m_components.push_back(5+SoArealData::omegaz);

    } else if (amrex::toLower(var).compare("omega_x") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::omegax);

    } else if (amrex::toLower(var).compare("omega_y") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::omegay);

    } else if (amrex::toLower(var).compare("omega_z") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::omegaz);

    } else if (amrex::toLower(var).compare("statwt") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::statwt);

    } else if (amrex::toLower(var).compare("dragcoeff") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::drag_coeff);

    } else if (amrex::toLower(var).compare("drag") == 0) {

      variable_names.push_back("drag_x");
      m_components.push_back(5+SoArealData::vel_source_x);

      variable_names.push_back("drag_y");
      m_components.push_back(5+SoArealData::vel_source_y);

      variable_names.push_back("drag_z");
      m_components.push_back(5+SoArealData::vel_source_z);

    } else if (amrex::toLower(var).compare("drag_x") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::vel_source_x);

    } else if (amrex::toLower(var).compare("drag_y") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::vel_source_y);

    } else if (amrex::toLower(var).compare("drag_z") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::vel_source_z);

    } else if (amrex::toLower(var).compare("t_s") == 0) {

      variable_names.push_back(var);
      m_components.push_back(5+SoArealData::temperature);

    } else if (amrex::toLower(var).compare("phase") == 0) {

      variable_names.push_back(var);
      const int idx = 5+SoArealData::count;
      m_components.push_back(idx+SoAintData::phase);

    } else if (amrex::toLower(var).compare("state") == 0) {

      variable_names.push_back(var);
      const int idx = 5+SoArealData::count;
      m_components.push_back(idx+SoAintData::state);

    } else if (amrex::toLower(var).compare("x_sn") == 0) {

      for (int n_s(0); n_s < m_solids->nspecies(); ++n_s) {
        variable_names.push_back("X_sn_"+m_solids->species_names(n_s));
        const int idx = 5+SoArealData::count+SoAintData::count;
        m_components.push_back(idx+n_s);
      }

    } else if (amrex::toLower(var).substr(0,5).compare("x_sn_") == 0) {

      const int var_name_size = var.size();
      const std::string var_species = var.substr(5,var_name_size-5);

      for (int n_s(0); n_s < m_solids->nspecies(); ++n_s) {
        if (var_species.compare(m_solids->species_names(n_s)) == 0) {
          variable_names.push_back("X_sn_"+m_solids->species_names(n_s));
          const int idx = 5+SoArealData::count+SoAintData::count;
          m_components.push_back(idx+n_s);
          break;
        }
      }

    } else {

      amrex::Abort("Unrecognized variable to monitor. Fix inputs");
    }
  }

  m_variables.swap(variable_names);
}


void
BaseMonitor::count ()
{
  for (int lev(0); lev < m_nlev; lev++) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<int> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {0};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& /*ptile_data*/,
                                           const int /*index*/,
                                           const int /*i*/) -> ReduceTuple
    {
      return 1;
    };

    ReduceTuple host_tuple = apply(lev, INT_MAX, reduce_data, reduce_op, R, default_value);

    int l_count = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceIntSum(&l_count, 1);

    m_counts[lev] = l_count;
  }
}


void
BaseMonitor::check_realbox_is_ok () const
{
  AMREX_ASSERT_WITH_MESSAGE(m_region->volume() > 0., "RealBox has non-positive volume");
}


void
GeneralProperty::monitor (const Real& /*dt*/)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  if (m_specs[1].compare("sum") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = sum(lev, m_components);
    }

  } else if (m_specs[1].compare("min") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = min(lev, m_components);
    }

  } else if (m_specs[1].compare("max") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = max(lev, m_components);
    }

  } else if (m_specs[1].compare("massweightedsum") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = massweightedsum(lev, m_components);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
GeneralProperty::sum (const int lev,
                      const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);

      return statwt*p_value;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_sum = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_sum, 1);

    result[var] = l_sum;
  }

  return result;
}


Vector<Real>
GeneralProperty::min (const int lev,
                      const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::max()};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const Real p_value = get_value(ptile_data, index, i);
      return p_value;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_value);

    Real l_min = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMin(&l_min, 1);

    result[var] = l_min;
  }

  return result;
}


Vector<Real>
GeneralProperty::max (const int lev,
                      const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    ReduceTuple default_value = {std::numeric_limits<Real>::min()};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const Real p_value = get_value(ptile_data, index, i);
      return p_value;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_value);

    Real l_max = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealMax(&l_max, 1);

    result[var] = l_max;
  }

  return result;
}


Vector<Real>
GeneralProperty::massweightedsum (const int lev,
                                  const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real volume = SoArealData::volume(particle.rdata(SoArealData::radius));
      const Real mass = particle.rdata(SoArealData::density)*volume;
      const Real statwt = particle.rdata(SoArealData::statwt);

      return statwt*mass*p_value;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_sum = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_sum, 1);

    result[var] = l_sum;
  }

  return result;
}


void
AveragedProperty::monitor (const Real& /*dt*/)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  if (m_specs[1].compare("average") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = average(lev, m_components);
    }

  } else if (m_specs[1].compare("standarddeviation") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = stddev(lev, m_components);
    }

  } else if (m_specs[1].compare("massweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = mass_weighted_average(lev, m_components);
    }

  } else if (m_specs[1].compare("volumeweightedaverage") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = volume_weighted_average(lev, m_components);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
AveragedProperty::average (const int lev,
                           const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = ReduceData<Real,Real>::Type;

    ReduceTuple default_value = {0., 0.};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);

      return {statwt*p_value, statwt};
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / denominator;
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


Vector<Real>
AveragedProperty::stddev (const int lev,
                          const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  const Vector<Real> avg = average(lev, indexes);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = ReduceData<Real,Real>::Type;

    ReduceTuple default_value = {0., 0.};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    const Real d_avg = avg[var];

    auto R = [get_value,d_avg] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                                 const int index,
                                                 const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);

      const Real diff = p_value - d_avg;

      return {statwt*diff*diff, statwt};
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = std::sqrt(numerator / denominator);
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


Vector<Real>
AveragedProperty::mass_weighted_average (const int lev,
                                         const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = ReduceData<Real,Real>::Type;

    ReduceTuple default_value = {0., 0.};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);
      const Real volume = SoArealData::volume(particle.rdata(SoArealData::radius));
      const Real mass = particle.rdata(SoArealData::density)*volume;

      return {statwt*mass*p_value, statwt*mass};
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / denominator;
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


Vector<Real>
AveragedProperty::volume_weighted_average (const int lev,
                                           const Vector<int>& indexes)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<Real,Real> reduce_data(reduce_op);

    using ReduceTuple = ReduceData<Real,Real>::Type;

    ReduceTuple default_value = {0., 0.};

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value] AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                                           const int index,
                                           const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);
      const Real volume = SoArealData::volume(particle.rdata(SoArealData::radius));

      return {statwt*volume*p_value, statwt*volume};
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, default_value);

    Real numerator = amrex::get<0>(host_tuple);
    Real denominator = amrex::get<1>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&numerator, 1);
    ParallelDescriptor::ReduceRealSum(&denominator, 1);

    if (std::abs(denominator) > 1.e-15)
      result[var] = numerator / denominator;
    else {
      AMREX_ASSERT(std::abs(numerator) < 1.e-15);
      result[var] = 0.;
    }
  }

  return result;
}


void
FlowRate::set_direction ()
{
  for (int dir(0); dir < AMREX_SPACEDIM; ++dir)
    if(std::abs(m_plane->lo(dir) - m_plane->hi(dir)) < 1.e-15) {
      m_direction = dir;
      break;
    }

  AMREX_ASSERT(AMREX_D_TERM(std::abs(m_plane->length(m_direction)) < 1.e-15, &&
                                     m_plane->length((m_direction+1)%AMREX_SPACEDIM) > 1.e-15, &&
                                     m_plane->length((m_direction+2)%AMREX_SPACEDIM) > 1.e-15));
}


void
FlowRate::check_plane_is_ok () const
{
  AMREX_ASSERT(std::abs(m_plane->volume()) < 1.e-15);
}


void
FlowRate::set_coordinate ()
{
  AMREX_ASSERT(std::abs(m_plane->lo(m_direction) - m_plane->hi(m_direction)) < 1.e-15);

  m_coordinate = m_plane->lo(m_direction);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
FlowRate::CrossesFlowPlane::operator() (const Real& position,
                                        const Real& velocity,
                                        const Real& plane_coordinate,
                                        const Real& dt) const
{
  // If dt is not strictly positive, return false
  if (dt < 1.e-15) {
    return 0;
  }

  const Real min_velocity = (plane_coordinate - position) / dt;

  const Real abs_min_velocity = std::abs(min_velocity);
  const Real abs_velocity = std::abs(velocity);

  const int it_crosses = (min_velocity*velocity > 0.) &&
    (std::abs(abs_min_velocity - abs_velocity) < 1.e-15 || abs_velocity > abs_min_velocity);

  return it_crosses;
}


void
FlowRate::monitor (const Real& dt)
{
  m_monitoring_results.clear();
  m_monitoring_results.resize(m_nlev, Vector<Real>());

  if (m_specs[1].compare("flowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = flow_rate(lev, m_components, dt);
    }

  } else if (m_specs[1].compare("massweightedflowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = mass_weighted_flow_rate(lev, m_components, dt);
    }

  } else if (m_specs[1].compare("volumeweightedflowrate") == 0) {

    for (int lev(0); lev < m_nlev; ++lev) {
      m_monitoring_results[lev] = volume_weighted_flow_rate(lev, m_components, dt);
    }

  } else {

    Print() << "Unrecognized monitor type: " << m_specs[1] << "\n";
    amrex::Abort("Inputs error");
  }
}


Vector<Real>
FlowRate::flow_rate (const int lev,
                     const Vector<int>& indexes,
                     const Real dt)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    const int direction = m_direction;
    const Real plane_coordinate = m_coordinate;

    CrossesFlowPlane crossing_check;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value,crossing_check,direction,plane_coordinate,dt]
      AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                        const int index,
                        const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);

      const Real pos = particle.pos(direction);
      const Real vel = particle.rdata(SoArealData::velx+direction);

      const Real sign = std::abs(vel) < 1.e-15 ? 0. : (vel > 0. ? 1. : -1.);

      if (crossing_check(pos, vel, plane_coordinate, dt))
        return statwt*p_value*sign;
      else
        return 0.;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_flow_rate = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_flow_rate, 1);

    result[var] = l_flow_rate;
  }

  return result;
}


Vector<Real>
FlowRate::mass_weighted_flow_rate (const int lev,
                                   const Vector<int>& indexes,
                                   const Real dt)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    const int direction = m_direction;
    const Real plane_coordinate = m_coordinate;

    CrossesFlowPlane crossing_check;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value,crossing_check,direction,plane_coordinate,dt]
      AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                        const int index,
                        const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);
      const Real volume = SoArealData::volume(particle.rdata(SoArealData::radius));
      const Real mass = particle.rdata(SoArealData::density)*volume;

      const Real pos = particle.pos(direction);
      const Real vel = particle.rdata(SoArealData::velx+direction);

      const Real sign = std::abs(vel) < 1.e-15 ? 0. : (vel > 0. ? 1. : -1.);

      if (crossing_check(pos, vel, plane_coordinate, dt))
        return statwt*mass*p_value*sign;
      else
        return 0.;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_flow_rate = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_flow_rate, 1);

    result[var] = l_flow_rate;
  }

  return result;
}


Vector<Real>
FlowRate::volume_weighted_flow_rate (const int lev,
                                     const Vector<int>& indexes,
                                     const Real dt)
{
  const int var_nb = indexes.size();

  Vector<Real> result(var_nb, 0.);

  for (int var(0); var < var_nb; ++var) {

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    const int direction = m_direction;
    const Real plane_coordinate = m_coordinate;

    CrossesFlowPlane crossing_check;

    GetParticleValue get_value(m_pc->m_runtimeRealData);

    auto R = [get_value,crossing_check,direction,plane_coordinate,dt]
      AMREX_GPU_DEVICE (const ParticleTileData& ptile_data,
                        const int index,
                        const int i) -> ReduceTuple
    {
      const auto& particle = ptile_data.getSuperParticle(i);

      const Real p_value = get_value(ptile_data, index, i);
      const Real statwt = particle.rdata(SoArealData::statwt);
      const Real volume = SoArealData::volume(particle.rdata(SoArealData::radius));

      const Real pos = particle.pos(direction);
      const Real vel = particle.rdata(SoArealData::velx+direction);

      const Real sign = std::abs(vel) < 1.e-15 ? 0. : (vel > 0. ? 1. : -1.);

      if (crossing_check(pos, vel, plane_coordinate, dt))
        return statwt*volume*p_value*sign;
      else
        return 0.;
    };

    ReduceTuple host_tuple = apply(lev, indexes[var], reduce_data, reduce_op, R, 0.);

    Real l_flow_rate = amrex::get<0>(host_tuple);

    ParallelDescriptor::ReduceRealSum(&l_flow_rate, 1);

    result[var] =  l_flow_rate;
  }

  return result;
}

} // end namespace LagrangianMonitor
