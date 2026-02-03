#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>

#include <mfix_actual_rhs.H>
#include <mfix_reporter.H>
#include <integrator_parameters.H>
#include <mfix_chem_memory_handler.H>
#include <mfix_fix_inputs.H>


/**
 * \brief This is the Base ODE integrator class. All the ODE integrators derive
 * from this base class
 *
 * The Integrator class contains all the data structures and information that is
 * common to all the ODE integrators in MFIX-Exa. It contains debug information
 * such as the integration success/error code, cell indexes, start and stop
 * times, tolerances, verbosity, etc... It also contains the pointers to the ODE
 * solution and ODE right-hand side.
 */
class Integrator
{
  public:
    /**
     * \brief Constructs the Integrator class with the inputs values from the
     * integrator parameters object
     *
     * This class constructor assigns default values to the time variables,
     * debug information, and integration flags. It also assigns the user's
     * inputs values to the solver settings and tolerances
     *
     * \param params is the object that carries the user's inputs parameters
     * information
     */
    AMREX_GPU_HOST_DEVICE
    Integrator (const IntegratorParms& params)
      : m_integration_status(IntegrationStatus::IERR_SUCCESS)
      , m_burner_verbose(params.m_burner_verbose)
      , m_ode_max_steps(params.m_ode_max_steps)
      , m_jacobian_type(params.m_jacobian_type)
      , m_jacobian_inputs(params.m_jacobian_type)
      , m_atol(params.m_atol)
      , m_rtol(params.m_rtol)
      , m_retry_atol(params.m_retry_atol)
      , m_retry_rtol(params.m_retry_rtol)
      , m_retry_swap_jacobian(params.m_retry_swap_jacobian)
      , m_ode_max_dt(params.m_ode_max_dt)
    {}


    /**
     * \brief Defines the integrator variables
     *
     * This member function defines the ODE problem size, aka the number of
     * variables, assigns the ODE solution and right-hand side pointers of this
     * class' members, assigns the debug information, and sets the current time
     * for the ODE integration
     *
     * \param memory_handler is the object that contains the information about
     * the memory allocated for the ODE integration
     * \param ijk are the indexes of the cell where the chemistry is being
     * solved
     * \param time is the current simulation time
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void define (DeviceChemMemoryHandler& memory_handler,
                 const amrex::IntVect& ijk,
                 const amrex::Real time)
    {
      // Sanity check that the memory_handler argument is correctly sized
      AMREX_ASSERT(memory_handler.nvariables() > 0);

      // Assign the memory handler
      m_mem_handler = &memory_handler;

      // Assign the problem size
      m_N = m_mem_handler->nvariables();

      // Assign the memory pointers
      m_y.define(m_mem_handler->variables_address(), m_N, m_mem_handler->stride());
      m_y_rhs.define(m_mem_handler->rhs_address(), m_N, m_mem_handler->stride());
      m_y_old.define(m_mem_handler->y_old_address(), m_N, m_mem_handler->stride());

      // Assign the debug information for the current ODE integration cell
      AMREX_D_DECL(m_i = ijk[0], m_j = ijk[1], m_k = ijk[2]);

      // Assign current time
      m_time = time;
    }


    /*
     * \brief Set the status of the integrator
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set_status (const IntegrationStatus status) { m_integration_status = status; }


    /*
     * \brief Get a reference to the integrator status variable
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntegrationStatus& integration_status () { return m_integration_status; }


    /*
     * \brief Get how to compute the jacobian of the ODE's rhs
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int jacobian_type () const { return m_jacobian_inputs; }


    /**
     * \brief Change the ODE problem size
     *
     * This member function can be useful when switching between eulerian and
     * lagrangian chemical reactions, to reduce the size of the ODE problem that
     * will be integrated. Indeed, in case of eulerian reactions, we can ignore
     * the solids variables in the ODE solution and right-hand side, thus
     * reducing the size of the problem to be solved, and make the integrator
     * more efficient
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void change_size (const int N)
    {
      AMREX_ASSERT(N > 0);
      m_N = N;
    }


    /**
     * \brief Calls the right-hand side function evaluation
     *
     * \param time is the current integration time
     * \param a_y is the vector of variables that will be used to compute the
     * right-hand side
     * \param a_y_rhs is the vector where the right-hand side computed quantities
     * will be stored
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void rhs (const amrex::Real time,
              ChemVector<amrex::Real>& a_y,
              ChemVector<amrex::Real>& a_y_rhs)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("Integrator::rhs()");
#endif

      m_time = time;

      // We are integrating a system of
      // dy/dt(1:N) = y_rhs(1:N)

      // we come in with (integrator type) state having the current
      // solution (or predicted) state and we need to copy this over to the
      // state to interface with the actual reaction network
      //
      // Fix the state as necessary -- this ensures that the mass fractions that
      // enter are valid (and optionally normalized)

      // Before the right-hand side computation, clean the state by constraining
      // the mass fractions within the [0,1] interval and by normalizing them
      clean_state(a_y);

      // at this point, the state is synchronized with the integrator and is
      // thermodynamically consistent.

      // Call the specific network routine to get the RHS.
      actual_rhs(a_y, a_y_rhs);
    }


    /**
     * \brief Define some of the integrator parameters
     *
     * \param dt is the integration timestep
     * \param is_retry is the flag to specify whether we're retrying the
     * integration or not
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void integrator_setup (const amrex::Real dt,
                           bool is_retry)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("Integrator::integrator_setup()");
#endif

      // Set the integration tolerances
      if (is_retry) {
        m_atol = m_retry_atol;
        m_rtol = m_retry_rtol;
      }

      // set the Jacobian type
      if (is_retry && m_retry_swap_jacobian) {
        AMREX_ASSERT(m_jacobian_type == JacobianType::Analytical ||
                     m_jacobian_type == JacobianType::Numerical ||
                     m_jacobian_type == JacobianType::Broyden);

        if (m_jacobian_type == JacobianType::Analytical) {
          if (m_jacobian_inputs == JacobianType::Numerical) {
            m_jacobian_type = JacobianType::Numerical;
          } else {
            m_jacobian_type = JacobianType::Broyden;
          }
        } else {
          m_jacobian_type = JacobianType::Analytical;
        }
      }

      // Start off by assuming a successful burn.
      m_integration_status = IntegrationStatus::IERR_SUCCESS;

      // Initialize the integration time.
      m_time = 0.0;
      m_tout = dt;
    }


    /**
     * \brief Store the initial state into a backup state for debugging and
     * verbosity purposes
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void integrator_backup ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("Integrator::integrator_backup()");
#endif

      for (uint64_t n(1); n <= m_N; ++n) {
        m_y_old(n) = m_y(n);
      }
    }


    /**
     * \brief Check whether the computed solution is physically feasible or not
     *
     * \param dt is the integration timestep
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void integrator_cleanup (amrex::Real dt)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("Integrator::integrator_cleanup()");
#endif

      // VODE does not always fail even though it can lead to unphysical states.
      // Add some checks that indicate a burn fail even if VODE thinks the
      // integration was successful.

      // Check whether the mass fractions have physically feasible values
      {
        auto const& indexes = m_mem_handler->indexes();

        if (indexes.fluid_indexes().count() > 0)
        {
          constexpr int species_begin = ChemMemoryIndexes::fluid_species_begin();
          const int species_end = indexes.fluid_species_end();

          amrex::Real X_tot(0.);
          for (int n(species_begin); n < species_end; ++n) {
            X_tot += m_y[n];
          }

          for (int n(species_begin); n < species_end; ++n) {
            amrex::Real X_sn = m_y[n] / X_tot;

            if (X_sn < -m_species_failure_tolerance ||
                X_sn > (1. + m_species_failure_tolerance)) {
              m_integration_status = IntegrationStatus::IERR_INVALID_MASS_FRACTIONS;
            }
          }
        }

        if (indexes.solids_indexes().count() > 0)
        {
          const int species_begin = indexes.solids_species_begin();
          const int species_end = indexes.solids_species_end();

          amrex::Real X_tot(0.);
          for (int n(species_begin); n < species_end; ++n) {
            X_tot += m_y[n];
          }

          if (m_y[indexes.solids_density()] > 0 &&
              m_y[indexes.solids_radius()] > 0 &&
              !amrex::almostEqual(m_y[indexes.solids_density()], 0.) &&
              !amrex::almostEqual(m_y[indexes.solids_radius()], 0.)) {

            for (int n(species_begin); n < species_end; ++n) {
              amrex::Real X_sn = m_y[n] / X_tot;

              if (X_sn < -m_species_failure_tolerance ||
                  X_sn > (1. + m_species_failure_tolerance)) {
                m_integration_status = IntegrationStatus::IERR_INVALID_MASS_FRACTIONS;
              }
            }
          } else {
            m_integration_status = IntegrationStatus::IERR_INVALID_SOLIDS_MASS;
          }
        }
      }

#ifndef AMREX_USE_GPU
      if (m_burner_verbose) {
        // Print out some integration statistics, if desired.
        printf("integration summary:\n"
               "number of steps taken: %d\n"
               "number of rhs evaluations: %d\n", m_n_step, m_n_rhs);
      }
#endif

      // If we failed, print out the current state of the integration.
      if (m_integration_status != IntegrationStatus::IERR_SUCCESS &&
          m_integration_status != IntegrationStatus::IERR_INVALID_SOLIDS_MASS) {
#ifndef AMREX_USE_GPU
        std::cout << amrex::Font::Bold << amrex::FGColor::Red
                  << "[ERROR] integration failed in net" << amrex::ResetDisplay << std::endl;
        std::cout << "error description: " << IntegrationErrors[static_cast<int>(m_integration_status)] << std::endl;
        if (m_integration_status == IntegrationStatus::IERR_INVALID_MASS_FRACTIONS) {
          std::cout << "  Integrator exited succefully, but a check on the data values failed" << std::endl;
        }
        std::cout << "zone = (" << m_i << ", " << m_j << ", " << m_k << ")" << std::endl;
        std::cout << "time = " << m_time << std::endl;
        std::cout << "dt = " << std::setprecision(16) << dt << std::endl;
        std::cout << "y start = ";
        for (uint64_t n(1); n <= m_N; ++n) {
          std::cout << std::setprecision(16) << m_y_old(n) << " ";
        }
        std::cout << std::endl;
        std::cout << "y current = ";
        for (uint64_t n(1); n <= m_N; ++n) {
          std::cout << std::setprecision(16) << m_y(n) << " ";
        }
        std::cout << std::endl;
#elif defined(ALLOW_GPU_PRINTF)
        AMREX_DEVICE_PRINTF("[ERROR] integration failed in net\n error description = %s, time = %g, dt = %g\n",
            IntegrationErrors[static_cast<int>(m_integration_status)], m_time, dt);
#endif
      }
    }


    /**
     * \brief Clean the state by constraining the mass fractions into the [0,1]
     * interval and normalizing the fluid and solids phases mass fractions
     *
     * \param a_y is the vector containing the thermodynamic state to clean
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void clean_state (ChemVector<amrex::Real>& a_y)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("Integrator::clean_state()");
#endif

      auto const& indexes = m_mem_handler->indexes();

      // Normalize fluid species mass fractions
      if (indexes.fluid_indexes().count() > 0)
      {
        constexpr int species_begin = ChemMemoryIndexes::fluid_species_begin();
        const int species_end = indexes.fluid_species_end();

        // Auxiliary variable computed by summing up species mass fractions
        // constrained within the [0,1] interval
        amrex::Real X_tot(0);

        for (int n(species_begin); n < species_end; ++n) {
          a_y[n] = amrex::Clamp(a_y[n], 0., 1.);
          X_tot += a_y[n];
        }

        //
        if (!amrex::almostEqual(X_tot, 1.))
        {
          for (int n(species_begin); n < species_end; ++n) {
            a_y[n] = a_y[n] / X_tot;
          }
        }
      }

      // Normalize solids species mass fractions
      if (indexes.solids_indexes().count() > 0 &&
          m_y[indexes.solids_density()] > 0 &&
          m_y[indexes.solids_radius()] > 0 &&
          !amrex::almostEqual(m_y[indexes.solids_density()], 0.) &&
          !amrex::almostEqual(m_y[indexes.solids_radius()], 0.))
      {
        const int species_begin = indexes.solids_species_begin();
        const int species_end = indexes.solids_species_end();

        // Auxiliary variable computed by summing up species mass fractions
        // constrained within the [0,1] interval
        amrex::Real X_tot(0);

        for (int n(species_begin); n < species_end; ++n) {
          a_y[n] = amrex::Clamp(a_y[n], 0., 1.);
          X_tot += a_y[n];
        }

        //
        if (!amrex::almostEqual(X_tot, 1.))
        {
          for (int n(species_begin); n < species_end; ++n) {
            a_y[n] = a_y[n] / X_tot;
          }
        }
      }
    }


    // The class that computes the righ-hand side
    MFIXActualRHS actual_rhs;


    /**
     * \brief Compute the analytical jacobian
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void actual_jac (ChemVector<amrex::Real>& /*xn*/,
                     ChemMatrix<amrex::Real>& /*jac*/)
    {
      amrex::Abort("Actual Jacobian function not implemented");
    }


  protected:
    static constexpr amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();

    // The data containing the solution
    ChemVector<amrex::Real> m_y;
    ChemVector<amrex::Real> m_y_rhs;
    ChemVector<amrex::Real> m_y_old;

    //
    DeviceChemMemoryHandler* m_mem_handler = nullptr;

    // The size of the problem to be solved
    uint64_t m_N = 0;

    // On exit, this will be set to the integration time we achieved
    amrex::Real m_time = 0.;
    amrex::Real m_tout = 0.;

    // for diagnostics / error reporting
    int m_i = -1;
    int m_j = -1;
    int m_k = -1;

    int m_n_rhs = 0;
    int m_n_jac = 0;
    int m_n_step = 0;

    // Was the integration successful?
    IntegrationStatus m_integration_status;

    int m_burner_verbose;

    // Maximum number of integration steps
    int m_ode_max_steps;

    // Switching between the analytical jacobian or the numerically approximated
    // one
    short m_jacobian_type;
    short m_jacobian_inputs;

    // Relative and absolute tolerances
    amrex::Real m_atol;
    amrex::Real m_rtol;

    amrex::Real m_retry_atol;
    amrex::Real m_retry_rtol;
    int m_retry_swap_jacobian;

    // ODE maximum timestep
    amrex::Real m_ode_max_dt;

    static constexpr amrex::Real m_species_failure_tolerance = 1.e-2;
};


#endif
