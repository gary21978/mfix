#include <mfix_timer.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>
#include <fstream>

#include <mfix_reporter.H>

using namespace amrex;

void
MFIXTimer::Initialize ()
{
  ParmParse pp("mfix");

  pp.query("fixed_dt", m_dt);

  if(m_dt > 0.) { m_timestep_type = TimestepType::Fixed; }
  else { m_timestep_type = TimestepType::Dynamic; }

  pp.query("dt_min", m_dt_min);
  pp.query("dt_max", m_dt_max);

  std::string walltime_limit_in;
  if (pp.query("walltime_limit", walltime_limit_in)) {
    int HH(0), MM(0), SS(0);
    if (sscanf(walltime_limit_in.c_str(), "%d:%d:%d", &HH, &MM, &SS) >= 2) {
      m_walltime_limit = static_cast<Real>(HH*3600 + MM*60 + SS);
    } else {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Unable to parse walltime limit " << walltime_limit_in << "\n"
        << "The expected dformat is HH:MM:SS. Please correct the input deck.";
    }
  }

  // Is this a steady-state calculation
  pp.query("steady_state", m_steady_state);

  int stop_time_set = pp.query("stop_time", m_stop_time);
  int max_step_set  = pp.query("max_step", m_max_step);

  pp.query("steady_state_tol", m_steady_state_tol);
  pp.query("steady_state_maxiter", m_steady_state_maxiter);

  if (!m_steady_state) {

    reporter::Log(reporter::Status) << "Simulation is transient";

    if ( stop_time_set && max_step_set) {

      reporter::Log(reporter::Warning)
        << "mfix.stop_time and mfix.max_step are both defined.\n"
        << "The simulation will end when either condition is met.";

    } else if (!stop_time_set && !max_step_set) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid inputs - Either mfix.stop_time or mfix.max_step must be defined!\n"
        << "Please correct the input deck.";

    } else if (m_steady_state_tol > 0.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid inputs - mfix.steady_state_tol is defined for transient simulation.\n"
        << "Please correct the input deck.";
    }
  } else { // steady_state

    if (m_steady_state_tol < 0.) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid inputs - mfix.steady_state_tol is required for steady-state simulation!\n"
        << "Please correct the input deck.";
    }

    reporter::Log(reporter::Status)
      << "Simulation is steady-state."
      << "\n  Max iterations: " << m_steady_state_maxiter
      << "\n  Tolerance " << m_steady_state_tol;

    if ( stop_time_set || max_step_set) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid inputs - Neither mfix.stop_time nor mfix.max_step should be defined!\n"
        << "Please correct the input deck.";
    }
  }

  pp.query("overstep_end_time", m_overstep_end_time);
  m_usr_max_step = m_max_step; // save the original max step selected by the user

  pp.query("clean_exit", m_clean_exit);
}


void
MFIXTimer::reset (const MFIXTimer& other)
{
  m_runtime_start = other.runtime_start();
  m_walltime_limit = other.walltime_limit();
  m_avg_step_runtime = other.avg_step_runtime();
  m_start_time = other.start_time();
  m_time = other.time();
  m_stop_time = other.stop_time();
  m_steady_state = other.SteadyState();
  m_timestep_type = other.timestep_type();
  m_dt = other.dt();
  m_dt_min = other.dt_min();
  m_dt_max = other.dt_max();
  m_first_step = other.first_step();
  m_nstep = other.nstep();
  m_max_step = other.max_step();
  m_usr_max_step = other.usr_max_step();
  m_clean_exit = other.clean_exit();
  m_run_status_type = other.run_status_type();
}


int
MFIXTimer::ok ()
{
  if ((m_stop_time >= 0.) && ((m_time+0.1*m_dt) >= m_stop_time))
  { m_run_status_type = MFIXRunStatusType::TimeIsOver; }

  if ((m_max_step >= 0) && (m_nstep >= m_max_step))
  { m_run_status_type = MFIXRunStatusType::IsFinalStep; }

  if ( m_steady_state ) { m_run_status_type = MFIXRunStatusType::SteadyState; }

  if (ParallelDescriptor::IOProcessor()) {

    if ((m_walltime_limit > 0.) && (!runtime_left_is_sufficient()))
    { m_run_status_type = MFIXRunStatusType::RuntimeIsOver; }

    if((m_clean_exit != "") && std::ifstream(m_clean_exit).good())
    { m_run_status_type = MFIXRunStatusType::UserStop; }
  }

  if (((m_walltime_limit > 0.)) || (m_clean_exit != "")) {
    ParallelDescriptor::Bcast(&m_run_status_type, 1, ParallelDescriptor::IOProcessorNumber());
  }

  return (m_run_status_type == MFIXRunStatusType::OK);
}


int
MFIXTimer::runtime_left_is_sufficient ()
{
  const Real needed_time = 1.2*(m_max_write_chkpt_time + m_avg_step_runtime);

  const Real missing_time = m_walltime_limit - elapsed_runtime();

  if (needed_time < missing_time)
    return 1;

  return 0;
}
