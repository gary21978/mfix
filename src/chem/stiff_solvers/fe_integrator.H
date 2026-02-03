#ifndef FE_INTEGRATOR_H
#define FE_INTEGRATOR_H

#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>

#include <integrator.H>
#include <mfix_actual_rhs.H>
#include <mfix_reporter.H>


/**
 * \brief This class integrates the chemistry ODE problem using the StiffSolver
 * Forward Euler integrator from AMReX-Astro/Microphysics
 */
class ForwardEuler: public Integrator
{
  public:
    /**
     * \brief This class defines the amount of auxiliary memory needed by the
     * StiffSolver ForwardEuler integrator during the ODE integration. The
     * auxiliary memory is used for temporary vectors and variables and for GPU
     * runs is either allocated on the device global or shared memory. The
     * aux_memory struct is needed to know how much auxiliary memory to be
     * allocated on the device memory.
     */
    struct aux_memory
    {
      /**
       * \brief Construct the aux_memory class
       *
       * This constructor gets the size of the ODE system, aka the number of
       * variables for the problem to integrate, and allocates the indexes for
       * auxiliary memory variables for this integrator by calling the
       * corresponding constructors.
       *
       * \param N is the number of variables (unknowns) for the ODE system
       */
      AMREX_GPU_HOST_DEVICE
      aux_memory (const int N,
                  [[maybe_unused]] const int jacobian_type = JacobianType::Numerical)
        : m_short1D(N)
        , m_real1D(N)
        , m_real2D(N)
      {}

      /**
       * \brief This class defines the indexes for the auxiliary variables of
       * type vector<short>
       */
      struct short1D
      {
        /**
         * \brief Constructs the short1D class
         *
         * No arguments are needed as the ForwardEuler integrator does NOT need
         * any vector<short> auxiliary variable
         */
        AMREX_GPU_HOST_DEVICE
        short1D (const int)
          : count(0)
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<short> auxiliary
         * memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(short); }

        // The number of auxiliary vector<short> auxiliary variables
        int count;
      } m_short1D;

      /**
       * \brief This class defines the indexes for the auxiliary variables of
       * type vector<Real>
       */
      struct real1D
      {
        /**
         * \brief Constructs the real1D class
         *
         * No arguments are needed as the ForwardEuler integrator does NOT need
         * any vector<Real> auxiliary variable
         */
        AMREX_GPU_HOST_DEVICE
        real1D (const int)
          : count(0)
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<Real> auxiliary
         * memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(amrex::Real); }

        int count;
      } m_real1D;

      /**
       * \brief This class defines the indexes for the auxiliary variables of
       * type vector<vector<Real>>
       */
      struct real2D
      {
        /**
         * \brief Constructs the real2D class
         *
         * No arguments are needed as the ForwardEuler integrator does NOT need
         * any vector<vector<Real>> auxiliary variable
         */
        AMREX_GPU_HOST_DEVICE
        real2D (const int)
          : count(0)
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<vector<Real>>
         * auxiliary memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(amrex::Real); }

        // The number of auxiliary vector<vector<Real>> auxiliary variables
        int count;
      } m_real2D;

      /**
       * \brief Get the amount of bytes needed for the auxiliary memory
       */
      AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
      size_t mem_bytes () const { return m_short1D.mem_bytes() +
                                         m_real1D.mem_bytes() +
                                         m_real2D.mem_bytes(); }

      /**
       * \brief Get the number of auxiliary vector<short> variables
       */
      AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
      int short_size () const { return m_short1D.count; }

      /**
       * \brief Get the number of auxiliary vector<Real> variables
       */
      AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
      int real1D_size () const { return m_real1D.count; }

      /**
       * \brief Get the number of auxiliary vector<vector<Real>> variables
       */
      AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
      int real2D_size () const { return m_real2D.count; }

    } m_aux_memory;


    /**
     * \brief Constructs the ForwardEuler class
     *
     * \param params is the reference to the IntegratorParms class which contain
     * all the user inputs and integrator settings for the constructor
     */
    AMREX_GPU_HOST_DEVICE
    ForwardEuler (const IntegratorParms& params)
      : Integrator(params)
      , m_aux_memory(0)
      , m_max_dt_change_factor(params.m_max_dt_change_factor)
    {}


    /**
     * \brief Defines the auxiliary memory number of variables and assigns the
     * auxiliary memory pointers to the corresponding allocated memory space
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup_aux_memory ()
    {
      // Sanity check that the current problem size, aka the ODE number of
      // variables, is strictly positive
      AMREX_ASSERT(m_N > 0);

      // Assign the correct size for the auxiliary variables
      m_aux_memory = aux_memory(m_N);

      const int stride = m_mem_handler->stride();

      ChemData<amrex::Real> aux_real1D(m_mem_handler->aux_real1D_mem(), m_aux_memory.real1D_size(), stride);
    }


    /**
     * \brief Integrate the ODE system by doing one step of the ForwardEuler
     * class
     *
     * \param dt is the timestep size
     * \param is_retry is the flag for retrying the ODE integration
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_aux_variables ()
    {}


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void actual_integrator (const amrex::Real& dt,
                            bool is_retry = false)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("ForwardEuler::actual_integrator()");
#endif

      reset_aux_variables();

      // Call the Integrator::integrator_setup function
      integrator_setup(dt, is_retry);

      // Call the Integrator::integrator_backup function
      integrator_backup();

      // Integrate the ODE system
      fe_integrator(dt);

      // Call the Integrator::integrator_cleanup function
      integrator_cleanup(dt);
    }


  private:
    /**
     * \brief Compute a timestep that allows the solution to do not change too
     * much
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calculate_dt ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("ForwardEuler::calculate_dt()");
#endif

      // Calculate the timestep.
      amrex::Real dt = 1.0e200;

      // Our timestepping strategy is to prevent any quantity from changing by
      // more than a certain factor in any timestep. We ignore this below atol.

      for (uint64_t n(1); n <= m_N; ++n) {

        if (m_y(n) >= m_atol) {

          amrex::Real target_dX;
          if (m_y_rhs(n) > 0.0) {
            target_dX = (m_max_dt_change_factor - 1.0) * m_y(n);
          } else {
            target_dX = (1.0 - 1.0 / m_max_dt_change_factor) * m_y(n);
          }

          amrex::Real dXdt = amrex::max(std::abs(m_y_rhs(n)), 1.0e-30);

          dt = amrex::min(dt, target_dX / dXdt);
        }
      }

      dt = amrex::min(dt, m_ode_max_dt);

      return dt;
    }


    /**
     * \brief The ForwardEuler integration algorithm
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void fe_integrator (amrex::Real dt)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("BackwardEuler::fe_integrator()");
#endif

      m_n_rhs = 0;
      m_n_jac = 0;
      m_n_step = 0;

      m_integration_status = IntegrationStatus::IERR_SUCCESS;

      while ((m_integration_status == IntegrationStatus::IERR_SUCCESS) &&
             (m_time < (1.0 - m_timestep_safety_factor) * dt) &&
             (m_n_step < m_ode_max_steps))
      {
        rhs(m_time, m_y, m_y_rhs);

        m_n_rhs += 1;

        amrex::Real dt_sub = calculate_dt();

        // Prevent the timestep from overshooting the final time.
        if (m_time + dt_sub > dt) {
          dt_sub = dt - m_time;
        }

        // Apply the update.
        for (uint64_t n(1); n <= m_N; ++n) {
          m_y(n) += m_y_rhs(n) * dt_sub;
        }

        m_time += dt_sub;
        ++m_n_step;
      }

      if (m_n_step >= m_ode_max_steps) {
        m_integration_status = IntegrationStatus::IERR_TOO_MANY_STEPS;
      }

      m_n_jac = 0;

      return;
    }


  private:
    amrex::Real m_max_dt_change_factor;

    static constexpr amrex::Real m_timestep_safety_factor = 1.0e-12;
};


#endif
