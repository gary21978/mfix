#ifndef BE_INTEGRATOR_H
#define BE_INTEGRATOR_H

#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>

#include <integrator.H>
#include <linear_algebra.H>
#include <mfix_actual_rhs.H>
#include <mfix_reporter.H>
#include <mfix_chem_memory_handler.H>


/**
 * \brief This class integrates the chemistry ODE problem using the StiffSolver
 * Backward Euler integrator from AMReX-Astro/Microphysics
 */
class BackwardEuler: public Integrator
{
  public:
    /**
     * \brief This class defines the amount of auxiliary memory needed by the
     * StiffSolver BackwardEuler integrator during the ODE integration. The
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
                  const int jacobian_type = JacobianType::Numerical)
        : m_short1D(N)
        , m_real1D(N, jacobian_type)
        , m_real2D(N, jacobian_type)
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
         * \param N is the number of variables (unknowns) for the ODE system
         */
        AMREX_GPU_HOST_DEVICE
        short1D (const int N)
          : pivot(0)
          , count(pivot+N)
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<short> auxiliary
         * memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(short); }

        int pivot;
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
        real1D (const int N,
                const int jacobian_type)
          : y_init(0)
          , y_copy(y_init+N)
          , b(y_copy+N)
          , y_fine(b+N)
          , y_rhs_plus(y_fine+N)
          , y_delta(y_rhs_plus+N)
          , y_rhs_delta(y_delta+N*(jacobian_type == JacobianType::Broyden))
          , count(y_rhs_delta+N*(jacobian_type == JacobianType::Broyden))
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<Real> auxiliary
         * memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(amrex::Real); }

        int y_init;
        int y_copy;
        int b;
        int y_fine;
        int y_rhs_plus;
        int y_delta;
        int y_rhs_delta;
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
         * \param N is the number of variables (unknowns) for the ODE system
         */
        AMREX_GPU_HOST_DEVICE
        real2D (const int N,
                const int jacobian_type)
          : jac(0)
          , jac_old(jac+N*N)
          , count(jac_old+N*N*(jacobian_type == JacobianType::Broyden))
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<vector<Real>>
         * auxiliary memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(amrex::Real); }

        int jac;
        int jac_old;
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
     * \brief Constructs the BackwardEuler class
     *
     * \param params is the reference to the IntegratorParms class which contain
     * all the user inputs and integrator settings for the constructor
     */
    AMREX_GPU_HOST_DEVICE
    BackwardEuler (const IntegratorParms& params)
      : Integrator(params)
      , m_aux_memory(0)
      , m_max_iter(params.m_max_iter)
      , m_linalg_do_pivoting(1)
      , m_rp_rtol(params.m_rp_rtol)
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
      m_aux_memory = aux_memory(m_N, m_jacobian_type);

      const int stride = m_mem_handler->stride();

      ChemData<short> aux_short1D(m_mem_handler->aux_short1D_mem(), m_aux_memory.short_size(), stride);
      ChemData<amrex::Real> aux_real1D(m_mem_handler->aux_real1D_mem(), m_aux_memory.real1D_size(), stride);
      ChemData<amrex::Real> aux_real2D(m_mem_handler->aux_real2D_mem(), m_aux_memory.real2D_size(), stride);

      m_pivot.define(aux_short1D.ptr(m_aux_memory.m_short1D.pivot), m_N, stride);

      m_y_init.define(aux_real1D.ptr(m_aux_memory.m_real1D.y_init), m_N, stride);
      m_y_copy.define(aux_real1D.ptr(m_aux_memory.m_real1D.y_copy), m_N, stride);
      m_b.define(aux_real1D.ptr(m_aux_memory.m_real1D.b), m_N, stride);
      m_y_fine.define(aux_real1D.ptr(m_aux_memory.m_real1D.y_fine), m_N, stride);
      m_y_rhs_plus.define(aux_real1D.ptr(m_aux_memory.m_real1D.y_rhs_plus), m_N, stride);

      if (m_jacobian_type == JacobianType::Broyden) {
        m_y_delta.define(aux_real1D.ptr(m_aux_memory.m_real1D.y_delta), m_N, stride);
        m_y_rhs_delta.define(aux_real1D.ptr(m_aux_memory.m_real1D.y_rhs_delta), m_N, stride);
      }

      m_jac.define(aux_real2D.ptr(m_aux_memory.m_real2D.jac), m_N, m_N, stride);

      if (m_jacobian_type == JacobianType::Broyden) {
        m_jac_old.define(aux_real2D.ptr(m_aux_memory.m_real2D.jac_old), m_N, m_N, stride);
      }
    }


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_aux_variables ()
    {
      m_pivot.fill(0.);

      m_y_init.fill(0.);
      m_y_copy.fill(0.);
      m_b.fill(0.);
      m_y_fine.fill(0.);
      m_y_rhs_plus.fill(0.);

      if (m_jacobian_type == JacobianType::Broyden) {
        m_y_delta.fill(0.);
        m_y_rhs_delta.fill(0.);
      }

      m_jac.fill(0.);

      if (m_jacobian_type == JacobianType::Broyden) {
        m_jac_old.fill(0.);
      }
    }


//    struct jac_info_t {
//      amrex::Real h;
//    };


    /**
     * \brief Integrate the ODE system by doing one step of the BackwardEuler
     * class
     *
     * \param dt is the timestep size
     * \param is_retry is the flag for retrying the ODE integration
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void actual_integrator (const amrex::Real& dt,
                            bool is_retry = false)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("BackwardEuler::actual_integrator()");
#endif

      reset_aux_variables();

      // Call the Integrator::integrator_setup function
      integrator_setup(dt, is_retry);
      // Call the Integrator::integrator_backup function
      integrator_backup();

      // Integrate the ODE system
      be_integrator();

      // Call the Integrator::integrator_cleanup function
      integrator_cleanup(dt);
    }


  private:
    /**
     * \brief Itial timestep estimator for VODE integrator
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real initial_react_dt ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("BackwardEuler::initial_react_dt()");
#endif

      // this is a version of the initial timestep estimator from VODE
      // we come in with y() storing the initial solution

      // initial lower and upper bounds on the timestep
      const amrex::Real hL = 100.0 * std::numeric_limits<amrex::Real>::epsilon() * m_tout;
      const amrex::Real hU = 0.1 * m_tout;

      // initial guess for the iteration
      amrex::Real h = std::sqrt(hL * hU);
      amrex::Real h_old = 10.0 * h;

      // save the old state and update the stored state
      for (uint64_t n(1); n <= m_N; ++n) {
        m_y_init(n) = m_y(n);
      }

      for (uint64_t iter = 1; iter <= 4; iter++) {

        h_old = h;

        // Construct the trial point.
        for (uint64_t n(1); n <= m_N; n++) {
          m_y(n) = m_y_init(n) + h * m_y_rhs(n);
        }

        // Call the RHS, then estimate the finite difference.

        rhs(m_time, m_y, m_y_rhs_plus);

        amrex::Real yddnorm = 0.0;
        for (uint64_t n(1); n <= m_N; n++) {
          // Get the error weighting -- this is similar to VODE's dewset routine
          const amrex::Real ewt = m_rtol * std::abs(m_y_init(n)) + m_atol;

          // Iterate on ddydtt = (RHS(t + h, y + h * dydt) - dydt) / h
          const amrex::Real ddydtt = (m_y_rhs_plus(n) - m_y_rhs(n)) / h;

          yddnorm += amrex::Math::powi<2>(ddydtt * ewt);
        }
        yddnorm = std::sqrt(yddnorm / m_N);

        if (yddnorm*hU*hU > 2.0) {
          h = std::sqrt(2.0 / yddnorm);
        } else {
          h = std::sqrt(h * hU);
        }

        if (h_old < 2.0 * h && h_old > 0.5 * h) {
          break;
        }
      }

      // Save the final timestep, with a bias factor.
      amrex::Real dt = h / 2.0;
      dt = amrex::Clamp(h, hL, hU);

      dt = amrex::min(dt, m_ode_max_dt);

      // restore the old time solution
      for (uint64_t n(1); n <= m_N; ++n) {
        m_y(n) = m_y_init(n);
      }

      return dt;
    }


    /**
     * \brief Creates current guess for the solution
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void single_step (const amrex::Real& dt)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("BackwardEuler::single_step()");
#endif

      m_integration_status = IntegrationStatus::IERR_SUCCESS;
      bool converged = false;

      // create our current guess for the solution -- just as a first order
      // explicit prediction
      rhs(m_time, m_y, m_y_rhs);

      m_n_rhs += 1;

      for (uint64_t n(1); n <= m_N; n++) {
        m_y_copy(n) = m_y(n);
        m_y(n) = m_y(n) + dt * m_y_rhs(n);
      }

      // Newton loop
      for (int iter = 1; iter <= m_max_iter; iter++) {

        // work with the current guess, and get the y_rhss for our current guess of y
        rhs(m_time, m_y, m_y_rhs);
        m_n_rhs += 1;

        // construct the Jacobian
        if (m_jacobian_type == JacobianType::Analytical) {
          actual_jac(m_y, m_jac);
        }
        else if (m_jacobian_type == JacobianType::Numerical) {
//          jac_info_t jac_info;
//          jac_info.h = dt;
//          numerical_jac(jac_info);
          numerical_jac();
          m_n_rhs += m_N;
        }
        else if (m_jacobian_type == JacobianType::Broyden) {
          if (iter == 1) {
//            jac_info_t jac_info;
//            jac_info.h = dt;
//            numerical_jac(jac_info);
            numerical_jac();
            m_n_rhs += m_N;
          } else {
            broyden_jac();
          }
        }

        m_n_jac++;

        // construct the matrix for the linear system
        // (I - dt J) dy^{n+1} = rhs
        for (uint64_t n(1); n <= m_N; n++) {
          for (uint64_t m(1); m <= m_N; m++) {
            m_jac(m,n) *= -dt;
          }
          m_jac(n,n) += 1.;
        }

        // construct the RHS of our linear system
        for (uint64_t n(1); n <= m_N; n++) {
          m_b(n) = m_y_copy(n) - m_y(n) + dt*m_y_rhs(n);
        }

        // solve the linear system
        int ierr_linpack(0);
        m_pivot.fill(0);

        if (m_linalg_do_pivoting == 1) {
          constexpr bool allow_pivot{true};
          dgefa<allow_pivot>(m_jac, m_N, m_pivot, ierr_linpack);
        } else {
          constexpr bool allow_pivot{false};
          dgefa<allow_pivot>(m_jac, m_N, m_pivot, ierr_linpack);
        }

        if (ierr_linpack != 0) {
          m_integration_status = IntegrationStatus::IERR_LU_DECOMPOSITION_ERROR;
          break;
        }

        if (m_linalg_do_pivoting == 1) {
          constexpr bool allow_pivot{true};
          dgesl<allow_pivot>(m_jac, m_N, m_pivot, m_b);
        } else {
          constexpr bool allow_pivot{false};
          dgesl<allow_pivot>(m_jac, m_N, m_pivot, m_b);
        }

        if (m_jacobian_type == JacobianType::Broyden) {
          for (uint64_t n(1); n <= m_N; n++) {
            m_y_delta(n) = m_y(n) - m_y_delta(n);
            m_y_rhs_delta(n) = m_y_rhs(n) - m_y_rhs_delta(n);
            for (uint64_t m(1); m <= m_N; m++) {
              m_jac_old(n,m) = m_jac(n,m);
            }
          }
        }

        // update our current guess for the solution
        for (uint64_t n(1); n <= m_N; n++) {
          m_y(n) += m_b(n);
        }

        // check to see if we converged
        // we compute the norms
        amrex::Real y_norm = 0.0;
        amrex::Real b_norm = 0.0;
        for (uint64_t n(1); n <= m_N; n++) {
          y_norm += m_y(n) * m_y(n);
          b_norm += m_b(n) * m_b(n);
        }
        y_norm = std::sqrt(y_norm / m_N);
        b_norm = std::sqrt(b_norm / m_N);

        if (b_norm < (m_rp_rtol * y_norm)) {
          converged = true;
          break;
        }
      }

      // we are done iterating -- did we converge?
      if (!converged) {
        if (m_integration_status == IntegrationStatus::IERR_SUCCESS) {

          // if we didn't set another error, then we probably ran
          // out of iterations, so set nonconvergence
          m_integration_status = IntegrationStatus::IERR_CORRECTOR_CONVERGENCE;

          // reset the solution to the original
          for (uint64_t n(1); n <= m_N; n++) {
            m_y(n) = m_y_copy(n);
          }
        }
      }

      return;
    }


    /**
     * \brief The BackwardEuler integration algorithm
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void be_integrator ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("BackwardEuler::be_integrator()");
#endif

      m_n_rhs = 0;
      m_n_jac = 0;
      m_n_step = 0;

      m_integration_status = IntegrationStatus::IERR_SUCCESS;

      // estimate the timestep
      rhs(m_time, m_y, m_y_rhs);

      m_n_rhs += 1;

      amrex::Real dt_sub = initial_react_dt();

      // main timestepping loop
      while ((m_integration_status == IntegrationStatus::IERR_SUCCESS) &&
             (m_time < (1.0 - m_timestep_safety_factor) * m_tout) &&
             (m_n_step < m_ode_max_steps))
      {
        // store the current solution -- we'll revert to this if a step fails
        for (uint64_t n(1); n <= m_N; ++n) {
          m_y_init(n) = m_y(n);
        }

        // don't go too far
        if (m_time + dt_sub > m_tout) {
          dt_sub = m_tout - m_time;
        }

        // our strategy is to take 2 steps at dt/2 and one at dt and
        // to compute the error from those

        // try to take a step dt
        // first do 2 (fine) dt/2 steps
        m_y_fine.fill(0.);

        single_step(dt_sub/2);

        if (m_integration_status == IntegrationStatus::IERR_SUCCESS) {

          single_step(dt_sub/2);

          for (uint64_t n(1); n <= m_N; ++n) {
            // store the fine dt solution
            m_y_fine(n) = m_y(n);

            // now that single (coarse) dt step
            // first reset the solution
            m_y(n) = m_y_init(n);
          }

          single_step(dt_sub);
        }

        amrex::Real rel_error = 0.0;

        for (uint64_t n(1); n <= m_N; n++) {
          // define a weight for each variable to use in checking the error
          const amrex::Real w = 1.0 / (m_rtol * std::abs(m_y_fine(n)) + m_atol);

          // now look for w |y_fine - y_coarse| < 1
          rel_error = amrex::max(rel_error, w * std::abs(m_y_fine(n) - m_y(n)));
        }

        bool step_success = false;
        if (rel_error < 1.0) {
          step_success = true;
        }

        if (m_integration_status == IntegrationStatus::IERR_SUCCESS && step_success) {
          // there were no problems with the mechanics of the step (linear algebra,
          // etc.) and we met our error goals.

          // y_fine has the current best solution
          m_time += dt_sub;

          for (uint64_t n(1); n <= m_N; ++n) {
            m_y(n) = m_y_fine(n);
          }

          // can we potentially increase the timestep?
          // backward-Euler has a local truncation error of dt**2
          amrex::Real dt_new = dt_sub * std::pow(1.0 / rel_error, 0.5);
          dt_sub = amrex::Clamp(dt_new, dt_sub/2.0, 2.0*dt_sub);

        } else {

          // roll back the solution
          for (uint64_t n(1); n <= m_N; ++n) {
            m_y(n) = m_y_init(n);
          }

          // adjust the timestep and try again
          dt_sub /= 2;

        }

        ++m_n_step;
        ++m_n_step;
      }

      if (m_n_step >= m_ode_max_steps) {
        m_integration_status = IntegrationStatus::IERR_TOO_MANY_STEPS;
      }

      return;
    }


    /**
     * \brief Computes a numerical first-order approximation of the Jacobian of
     * the right-hand side
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
//    void numerical_jac (const jac_info_t& jac_info)
    void numerical_jac ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("BackwardEuler::numerical_jac()");
#endif

//      // start by computing |f|, which is the norm of the RHS / weight
//      amrex::Real fac = 0.0;
//      amrex::Real weight(0.0);
//      for (uint64_t n(1); n <= m_N; n++) {
//        // the weight is 1 / EWT defined in VODE
//        weight = 1.0 / (m_rtol * std::abs(m_y(n)) + m_atol);
//        fac += (m_y_rhsm(n) * weight) * (m_y_rhsm(n) * weight);
//      }
//
//      fac = std::sqrt(fac / m_N);
//
//      // now compute r0, defined in the LSODE paper, Eq. 3.34
//      amrex::Real r0 = 1000.0 * std::abs(jac_info.h) * UROUND * m_N * fac;
//      if (r0 == 0.0) {
//        r0 = 1.0;
//      }

      // now we compute the Jacobian elements

      for (uint64_t n(1); n <= m_N; n++) {

//        weight = m_rtol * std::abs(yj) + m_atol;
//
//        // the incremement we use in the derivative is defined in the LSODE paper, Eq. 3.35
//        amrex::Real dy = amrex::max(std::sqrt(UROUND) * std::abs(yj), r0 * weight);
        const amrex::Real dy = std::sqrt(std::numeric_limits<amrex::Real>::epsilon()) * (std::abs(m_y(n)) + 1.e-8);

        m_y(n) += dy;

        actual_rhs(m_y, m_y_rhs_plus);

        // now fill in all of the rows for this column X_n
        for (uint64_t m(1); m <= m_N; m++) {
          m_jac(m,n) = (m_y_rhs_plus(m) - m_y_rhs(m)) / dy;
        }

        m_y(n) -= dy;
      }
    }


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void broyden_jac ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("BackwardEuler::broyden_jac()");
#endif

      amrex::Real Dx_sqr_norm2(0);
      for (uint64_t i(1); i <= m_N; ++i) {
        Dx_sqr_norm2 += m_y_delta(i)*m_y_delta(i);
      }

      for (uint64_t i(1); i <= m_N; ++i) {
        amrex::Real JDx(0.);
        for (uint64_t j(1); j <= m_N; ++j) {
          JDx += m_jac_old(i,j)*m_y_delta(j);
        }

        for (uint64_t j(1); j <= m_N; ++j) {
          m_jac(i,j) = m_jac_old(i,j) + ((m_y_rhs_delta(i) - JDx)/Dx_sqr_norm2)*m_y_delta(j);
        }
      }
    }


  private:
    int m_max_iter;

    int m_linalg_do_pivoting;

    amrex::Real m_rp_rtol;

    static constexpr amrex::Real m_timestep_safety_factor = 1.0e-12;

    ChemVector<short> m_pivot;

    ChemVector<amrex::Real> m_y_init;
    ChemVector<amrex::Real> m_y_copy;
    ChemVector<amrex::Real> m_b;
    ChemVector<amrex::Real> m_y_fine;
    ChemVector<amrex::Real> m_y_rhs_plus;
    ChemVector<amrex::Real> m_y_delta;
    ChemVector<amrex::Real> m_y_rhs_delta;

    ChemMatrix<amrex::Real> m_jac;
    ChemMatrix<amrex::Real> m_jac_old;
};


#endif
