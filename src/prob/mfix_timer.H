#ifndef _MFIX_TIMER_H_
#define _MFIX_TIMER_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_ParallelDescriptor.H>

#include <limits>
#include <string>

struct MFIXRunStatusType { enum {OK, TimeIsOver, IsFinalStep, RuntimeIsOver, UserStop, SteadyState}; };

class MFIXTimer
{
  public:
    struct TimestepType { enum {Undefined=-1, Fixed, Dynamic}; };

    // Constructor
    MFIXTimer ()
      : m_runtime_start(0.)
      , m_walltime_limit(-1.)
      , m_avg_step_runtime(0.)
      , m_max_write_chkpt_time(0.)
      , m_start_time(0.)
      , m_time(0.)
      , m_stop_time(std::numeric_limits<amrex::Real>::infinity())
      , m_steady_state(0)
      , m_steady_state_tol(-1.)
      , m_steady_state_maxiter(std::numeric_limits<int>::max())
      , m_timestep_type(TimestepType::Undefined)
      , m_dt(-1.)
      , m_dt_min(0.)
      , m_dt_max(1.e14)
      , m_prev_dt(-1.)
      , m_first_step(0)
      , m_nstep(0)
      , m_max_step(-1)
      , m_usr_max_step(-1)
      , m_clean_exit("")
      , m_run_status_type(MFIXRunStatusType::OK)
      , m_step_timing(std::numeric_limits<amrex::Real>::max())
    {
      m_start_time = system_time();
    }

    void Initialize ();

    void reset (const MFIXTimer& other);

    amrex::Real runtime_start () const { return m_runtime_start; }
    amrex::Real walltime_limit () const { return m_walltime_limit; }
    amrex::Real& walltime_limit () { return m_walltime_limit; }
    amrex::Real avg_step_runtime () const { return m_avg_step_runtime; }
    amrex::Real max_write_chkpt_time () const { return m_max_write_chkpt_time; }

    AMREX_FORCE_INLINE
    void update_avg_step_runtime (const amrex::Real step_time)
    {
      m_avg_step_runtime = (m_avg_step_runtime*m_nstep + step_time) / (m_nstep+1);
      amrex::ParallelDescriptor::ReduceRealMax(&m_avg_step_runtime, 1);

      if ((m_max_step > 0) && (m_walltime_limit > 0.))
        this->update_max_step();
    }

    AMREX_FORCE_INLINE
    void update_max_write_chkpt_time (const amrex::Real write_time)
    {
      m_max_write_chkpt_time = amrex::max(m_max_write_chkpt_time, write_time);
      amrex::ParallelDescriptor::ReduceRealMax(&m_max_write_chkpt_time, 1);
    }

    AMREX_FORCE_INLINE
    amrex::Real system_time () const { return amrex::ParallelDescriptor::second(); }

    AMREX_FORCE_INLINE
    amrex::Real elapsed_runtime () const
    { return elapsed_runtime(m_start_time); }

    AMREX_FORCE_INLINE
    amrex::Real elapsed_runtime (const amrex::Real start) const
    { return (system_time() - start); }

    amrex::Real start_time () const { return m_start_time; }
    amrex::Real stop_time () const { return m_stop_time; }

    int overstep_end_time () const {
      return (m_overstep_end_time && (m_stop_time > 0.));
    }

    int SteadyState () const { return m_steady_state; }

    amrex::Real SteadyStateTol () const { return m_steady_state_tol; }

    int SteadyStateMaxIter () const { return m_steady_state_maxiter; }

    amrex::Real  time () const { return m_time; }
    amrex::Real& time () { return m_time; }

    amrex::Real  new_time () const { return m_time + m_dt; }

    AMREX_FORCE_INLINE
    void advance_time () { m_time += m_dt; }

    int timestep_type () const { return m_timestep_type; }
    amrex::Real  dt () const { return m_dt; }
    amrex::Real& dt () { return m_dt; }

    amrex::Real  prev_dt () const { return m_prev_dt; }
    amrex::Real& prev_dt () { return m_prev_dt; }

    amrex::Real dt_min () const { return m_dt_min; }
    amrex::Real dt_max () const { return m_dt_max; }

    int first_step () const { return m_first_step; }
    int nstep () const { return m_nstep; }
    int& nstep () { return m_nstep; }
    int max_step () const { return m_max_step; }
    int usr_max_step () const { return m_usr_max_step; }

    AMREX_FORCE_INLINE
    void advance_nstep () { m_nstep++; }

    const std::string& clean_exit () const { return m_clean_exit; }

    int run_status_type () const { return m_run_status_type; }

    int ok ();
    int runtime_left_is_sufficient ();

    void step_start () { m_step_timing = system_time(); }

    amrex::Real step_end () {

      amrex::Real step_time = system_time() - m_step_timing;

      amrex::ParallelDescriptor::ReduceRealMax(step_time,
          amrex::ParallelDescriptor::IOProcessorNumber());

      if (walltime_limit() > 0.) {
        update_avg_step_runtime(step_time);
      }

      return step_time;
    }

    AMREX_FORCE_INLINE
    int test_per_approx ( amrex::Real a_dt, amrex::Real a_per_approx)
    {
      // Check to see if we've crossed a _per_approx interval by comparing
      // the number of intervals that have elapsed for both the current
      // time and the time at the beginning of this timestep.

      int num_per_old = static_cast<int>( (m_time-a_dt) / a_per_approx );
      int num_per_new = static_cast<int>( (m_time     ) / a_per_approx );

      // Before using these, however, we must test for the case where we're
      // within machine epsilon of the next interval. In that case, increment
      // the counter, because we have indeed reached the next par_ascii_per_approx interval
      // at this point.

      const amrex::Real sn = 10. * std::numeric_limits<amrex::Real>::epsilon();
      const amrex::Real eps = sn * amrex::Math::abs(m_time);

      const amrex::Real next_time = (num_per_old + 1) * a_per_approx;

      if ((num_per_new == num_per_old) && amrex::Math::abs(m_time - next_time) <= eps) {
          num_per_new += 1;
      }

      // Similarly, we have to account for the case where the old time is within
      // machine epsilon of the beginning of this interval, so that we don't double
      // count that time threshold -- we already plotted at that time on the last timestep.

      if ((num_per_new != num_per_old) && amrex::Math::abs((m_time - a_dt) - next_time) <= eps)
          num_per_old += 1;

      return (num_per_old != num_per_new) ? 1 : 0;
    }

  private:
    amrex::Real m_runtime_start;
    amrex::Real m_walltime_limit;
    amrex::Real m_avg_step_runtime;
    amrex::Real m_max_write_chkpt_time;

    amrex::Real m_start_time;
    amrex::Real m_time;
    amrex::Real m_stop_time;
    int m_overstep_end_time = 0;

    int m_steady_state;

    /// Tolerance to check for steady state.
    amrex::Real m_steady_state_tol;

    /// Maximum number of iterations to steady state --
    int m_steady_state_maxiter;

    int m_timestep_type;
    amrex::Real m_dt;
    amrex::Real m_dt_min;
    amrex::Real m_dt_max;

    amrex::Real m_prev_dt;

    int m_first_step;
    int m_nstep;
    int m_max_step;
    int m_usr_max_step;

    std::string m_clean_exit;

    int m_run_status_type;

    amrex::Real m_step_timing;

    AMREX_FORCE_INLINE
    void update_max_step ()
    {
      if (m_walltime_limit > 0. && m_avg_step_runtime > 0.) {
        m_max_step = std::min(m_usr_max_step, int(amrex::Math::floor(m_walltime_limit / m_avg_step_runtime)));
      }
    }
};

#endif
