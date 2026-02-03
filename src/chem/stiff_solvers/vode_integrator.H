#ifndef VODE_INTEGRATOR_H
#define VODE_INTEGRATOR_H

#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>

#include <integrator.H>
#include <linear_algebra.H>
#include <mfix_actual_rhs.H>
#include <mfix_reporter.H>

#include <mfix_chem_memory_structs.H>


// For the backward differentiation formula (BDF) integration
// the maximum order should be no greater than 5.
#define VODE_MAXORD 5
#define VODE_LMAX   (VODE_MAXORD+1)


/**
 * \brief This class integrates the chemistry ODE problem using the StiffSolver
 * VODE integrator from AMReX-Astro/Microphysics
 */
class VODE: public Integrator
{
  public:
    /**
     * \brief This class defines the amount of auxiliary memory needed by the
     * StiffSolver VODE integrator during the ODE integration. The auxiliary
     * memory is used for temporary vectors and variables and for GPU runs is
     * either allocated on the device global or shared memory. The aux_memory
     * struct is needed to know how much auxiliary memory to be allocated on the
     * device memory.
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
      aux_memory (const uint64_t N,
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
         * \param N is the number of variables (unknowns) for the ODE system
         */
        AMREX_GPU_HOST_DEVICE
        short1D (const uint64_t N)
          : pivot(0)
          , count(pivot+N)
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<short> auxiliary
         * memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(short); }

        uint64_t pivot;
        uint64_t count;
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
        real1D (const uint64_t N)
          : ewt(0)
          , savf(ewt+N)
          , acor(savf+N)
          , count(acor+N)
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<Real> auxiliary
         * memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(amrex::Real); }

        uint64_t ewt;
        uint64_t savf;
        uint64_t acor;
        uint64_t count;
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
        real2D (const uint64_t N)
          : jac(0)
          , yh(jac+N*N)
#ifdef VODE_JACOBIAN_CACHING
          , jac_save(yh+N*VODE_LMAX)
          , count(jac_save+N*N)
#else
          , count(yh+N*VODE_LMAX)
#endif
        {}

        /**
         * \brief Get the amount of bytes needed for the vector<vector<Real>>
         * auxiliary memory
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        size_t mem_bytes () const { return count*sizeof(amrex::Real); }

        uint64_t jac;
        uint64_t yh;
#ifdef VODE_JACOBIAN_CACHING
        uint64_t jac_save;
#endif
        uint64_t count;
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
      uint64_t short_size () const { return m_short1D.count; }

      /**
       * \brief Get the number of auxiliary vector<Real> variables
       */
      AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
      uint64_t real1D_size () const { return m_real1D.count; }

      /**
       * \brief Get the number of auxiliary vector<vector<Real>> variables
       */
      AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
      uint64_t real2D_size () const { return m_real2D.count; }

    } m_aux_memory;


    /**
     * \brief Constructs the VODE class
     *
     * \param params is the reference to the IntegratorParms class which contain
     * all the user inputs and integrator settings for the constructor
     */
    AMREX_GPU_HOST_DEVICE
    VODE (const IntegratorParms& params)
      : Integrator(params)
      , m_aux_memory(0)
#ifdef VODE_JACOBIAN_CACHING
      , m_use_jacobian_caching(params.m_use_jacobian_caching)
#endif
      , m_linalg_do_pivoting(params.m_linalg_do_pivoting)
      , m_X_reject_buffer(params.m_X_reject_buffer)
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

      const uint64_t stride = m_mem_handler->stride();

      ChemData<short> aux_short1D(m_mem_handler->aux_short1D_mem(), m_aux_memory.short_size(), stride);
      ChemData<amrex::Real> aux_real1D(m_mem_handler->aux_real1D_mem(), m_aux_memory.real1D_size(), stride);
      ChemData<amrex::Real> aux_real2D(m_mem_handler->aux_real2D_mem(), m_aux_memory.real2D_size(), stride);

      m_pivot.define(aux_short1D.ptr(m_aux_memory.m_short1D.pivot), m_N, stride);

      m_ewt.define(aux_real1D.ptr(m_aux_memory.m_real1D.ewt), m_N, stride);
      m_savf.define(aux_real1D.ptr(m_aux_memory.m_real1D.savf), m_N, stride);
      m_acor.define(aux_real1D.ptr(m_aux_memory.m_real1D.acor), m_N, stride);

      m_jac.define(aux_real2D.ptr(m_aux_memory.m_real2D.jac), m_N, m_N, stride);
      m_yh.define(aux_real2D.ptr(m_aux_memory.m_real2D.yh), m_N, VODE_LMAX, stride);
#ifdef VODE_JACOBIAN_CACHING
      m_jac_save.define(aux_real2D.ptr(m_aux_memory.m_real2D.jac_save), m_N, m_N, stride);
#endif
    }


    /**
     * \brief Integrate the ODE system by doing one step of the BackwardEuler
     * class
     *
     * \param dt is the timestep size
     * \param is_retry is the flag for retrying the ODE integration
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_aux_variables ()
    {
      m_pivot.fill(0.);

      m_ewt.fill(0.);
      m_savf.fill(0.);
      m_acor.fill(0.);

      m_jac.fill(0.);
      m_yh.fill(0.);
#ifdef VODE_JACOBIAN_CACHING
      m_jac_save.fill(0.);
#endif
    }


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void actual_integrator (const amrex::Real& dt,
                            bool is_retry = false)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::actual_integrator()");
#endif

      reset_aux_variables();

      // Call the Integrator::integrator_setup function
      integrator_setup(dt, is_retry);

      // Call the Integrator::integrator_backup function
      integrator_backup();

      // Integrate the ODE system
      dvode();

      // Call the Integrator::integrator_cleanup function
      integrator_cleanup(dt);
    }


  private:
    /**
     * \brief This routine computes the step size, H0, to be attempted on the
     * first step, when the user has not supplied a value for this.
     *
     * \param H0 is the computed timestep size
     * \param NITER is the number of iterations
     * \param IER is the error code
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void dvhin (amrex::Real& H0,
                int& NITER,
                int& IER)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::dvhin()");
#endif

      // This routine computes the step size, H0, to be attempted on the
      // first step, when the user has not supplied a value for this.

      // First we check that tout - T0 differs significantly from zero. Then
      // an iteration is done to approximate the initial second derivative
      // and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
      // A bias factor of 1/2 is applied to the resulting h.
      // The sign of H0 is inferred from the initial values of tout and T0.

      constexpr amrex::Real PT1 = 0.1e0;

      NITER = 0;
      amrex::Real TDIST = std::abs(m_tout - m_time);
      amrex::Real TROUND = UROUND * amrex::max(std::abs(m_time), std::abs(m_tout));

      if (TDIST < 2.0 * TROUND) {
        // Error return for tout - time too small.
        IER = -1;
        return;
      }

      // Set a lower bound on h based on the roundoff level in time and tout.
      const amrex::Real HLB = 100.0 * TROUND;

      // Set an upper bound on h based on tout-time and the initial Y and YDOT.
      amrex::Real HUB = PT1 * TDIST;

      for (uint64_t i(1); i <= m_N; ++i) {
        amrex::Real DELYI = PT1 * std::abs(m_yh(i,1)) + m_atol;
        amrex::Real AFI = std::abs(m_yh(i,2));
        if (AFI * HUB > DELYI) {
          HUB = DELYI / AFI;
        }
      }

      // Set initial guess for h as geometric mean of upper and lower bounds.
      int iter = 0;
      amrex::Real HG = std::sqrt(HLB * HUB);

      // If the bounds have crossed, exit with the mean value.
      bool do_iterations = true;
      if (HUB < HLB) {
        H0 = HG;
        do_iterations = false;
      }

      amrex::Real hnew{};

      if (do_iterations) {

        while (true) {

          // Estimate the second derivative as a difference quotient in f.
          const amrex::Real H = std::copysign(HG, m_tout - m_time);
          for (uint64_t i(1); i <= m_N; ++i) {
            m_y(i) = m_yh(i,1) + H * m_yh(i,2);
          }

          const amrex::Real t1 = m_time + H;

          rhs(t1, m_y, m_acor);

          amrex::Real YDDNRM = 0.0;

          for (uint64_t i(1); i <= m_N; ++i) {
            m_acor(i) = (m_acor(i) - m_yh(i,2)) / H;

            YDDNRM += (m_acor(i) * m_ewt(i)) * (m_acor(i) * m_ewt(i));
          }
          YDDNRM = std::sqrt(YDDNRM / m_N);

          // Get the corresponding new value of h.
          if (YDDNRM * HUB * HUB > 2.0) {
            hnew = std::sqrt(2.0 / YDDNRM);
          }
          else {
            hnew = std::sqrt(HG * HUB);
          }
          iter += 1;

          // Test the stopping conditions.
          // Stop if the new and previous h values differ by a factor of .lt. 2.
          // Stop if four iterations have been done.  Also, stop with previous h
          // if hnew/HG .gt. 2 after first iteration, as this probably means that
          // the second derivative value is bad because of cancellation error.

          if (iter >= 4) {
            break;
          }

          amrex::Real HRAT = hnew / HG;
          if ((HRAT > 0.5) && (HRAT < 2.0)) {
            break;
          }

          if ((iter >= 2) && (hnew > 2.0 * HG)) {
            hnew = HG;
            break;
          }

          HG = hnew;
        }

        // Iteration done. Apply bounds, bias factor, and sign. Then exit.
        H0 = hnew * 0.5;
        if (H0 < HLB) {
          H0 = HLB;
        }
        if (H0 > HUB) {
          H0 = HUB;
        }
      }

      H0 = std::copysign(H0, m_tout - m_time);
      NITER = iter;
      IER = 0;
    }


    /**
     * \brief Compute and process the matrix P = I - h*rl1*J , where J is an
     * approximation to the Jacobian that we obtain either through direct
     * evaluation or caching from a previous evaluation. P is then subjected to
     * LU decomposition in preparation for later solution of linear systems with
     * P as coefficient matrix. This is done by DGEFA.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void dvjac (int& IERPJ)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::dvjac()");
#endif

      // dvjac is called by dvnlsd to compute and process the matrix
      // P = I - h*rl1*J , where J is an approximation to the Jacobian
      // that we obtain either through direct evaluation or caching from
      // a previous evaluation. P is then subjected to LU decomposition
      // in preparation for later solution of linear systems with P as
      // coefficient matrix. This is done by DGEFA.

      IERPJ = 0;

#ifdef VODE_JACOBIAN_CACHING
      // See whether the Jacobian should be evaluated. Start by basing
      // the decision on whether we're caching the Jacobian.

      int evaluate_jacobian = 1;

      if (m_use_jacobian_caching) {
        evaluate_jacobian = 0;
      }

      if (m_use_jacobian_caching) {
        // Now evaluate the cases where we're caching the Jacobian but aren't
        // going to be using the cached Jacobian.

        // On the first step we don't have a cached Jacobian. Also, after enough
        // steps, we consider the cached Jacobian too old and will want to re-evaluate
        // it, so we look at whether the step of the last Jacobian evaluation (m_NSLJ)
        // is more than m_max_steps_between_jacobian_evals steps in the past.
        if (m_n_step == 0 || m_n_step > m_NSLJ + m_max_steps_between_jacobian_evals) {
          evaluate_jacobian = 1;
        }

        // See the non-linear solver for details on these conditions.
        if (m_ICF == 1 && m_DRC < m_CCMXJ) {
          evaluate_jacobian = 1;
        }

        if (m_ICF == 2) {
          evaluate_jacobian = 1;
        }
      }

      if (evaluate_jacobian == 1) {
#endif

        // We want to evaluate the Jacobian -- now the path depends on
        // whether we're using the numerical or analytic Jacobian.

        if (m_jacobian_type == JacobianType::Analytical) {

          // For the analytic Jacobian, call the user-supplied function.

          // Increment the Jacobian evaluation counter.
          m_n_jac += 1;

          // Refresh the timestep marker for the last Jacobian evaluation.
          m_NSLJ = m_n_step;

          // Indicate that the Jacobian is current for this solve.
          m_JCUR = 1;

          actual_jac(m_y, m_jac);

#ifdef VODE_JACOBIAN_CACHING
          // Store the Jacobian if we're caching.
          if (m_use_jacobian_caching == 1) {
            for (uint64_t i(1); i <= m_N; ++i)
              for (uint64_t j(1); j <= m_N; ++j)
                m_jac_save(i,j) = m_jac(i,j);
          }
#endif

        }
        else {

          // For the numerical Jacobian, make N calls to the RHS to approximate
          // it.

          // Increment the Jacobian evaluation counter.
          m_n_jac += 1;

          // Refresh the timestep marker for the last Jacobian evaluation.
          m_NSLJ = m_n_step;

          // Indicate that the Jacobian is current for this solve.
          m_JCUR = 1;

//          amrex::Real fac = 0.0;
//          for (uint64_t i(1); i <= m_N; ++i) {
//            fac += (m_savf(i) * m_ewt(i)) * (m_savf(i) * m_ewt(i));
//          }
//          fac = std::sqrt(fac / m_N);
//
//          amrex::Real R0 = 1000.0 * std::abs(m_H) * UROUND * m_N * fac;
//          if (R0 == 0.0) {
//            R0 = 1.0;
//          }

          for (uint64_t j(1); j <= m_N; ++j) {
            amrex::Real yj = m_y(j);

//            const amrex::Real R = amrex::max(std::sqrt(UROUND) * std::abs(yj), R0 / m_ewt(j));
            amrex::Real R = std::sqrt(std::numeric_limits<amrex::Real>::epsilon()) * (std::abs(yj) + 1.e-8);
            m_y(j) += R;
            R = 1.0 / R;

            rhs(m_tn, m_y, m_acor);
            for (uint64_t i(1); i <= m_N; ++i) {
              m_jac(i,j) = (m_acor(i) - m_savf(i)) * R;
            }

            m_y(j) = yj;
          }

          // Increment the RHS evaluation counter by N.
          m_n_rhs += m_N;

#ifdef VODE_JACOBIAN_CACHING
          // Store the Jacobian if we're caching.
          if (m_use_jacobian_caching == 1) {
            for (uint64_t i(1); i <= m_N; ++i)
              for (uint64_t j(1); j <= m_N; ++j)
                m_jac_save(i,j) = m_jac(i,j);
          }
#endif
        }

#ifdef VODE_JACOBIAN_CACHING
      }
      else {
        // Load the cached Jacobian.

        // Indicate the Jacobian is not current for this step.
        m_JCUR = 0;

        for (uint64_t i(1); i <= m_N; ++i)
          for (uint64_t j(1); j <= m_N; ++j)
            m_jac(i,j) = m_jac_save(i,j);
      }
#endif

      // Multiply Jacobian by a scalar, add the identity matrix
      // (along the diagonal), and do LU decomposition.

      const amrex::Real hrl1 = m_H * m_RL1;
      const amrex::Real con = -hrl1;

      for (uint64_t i(1); i <= m_N; ++i) {
        for (uint64_t j(1); j <= m_N; ++j) {
          m_jac(i,j) *= con;
        }
        m_jac(i,i) += 1.0;
      }

      int IER{};

      if (m_linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgefa<allow_pivot>(m_jac, m_N, m_pivot, IER);
      } else {
        constexpr bool allow_pivot{false};
        dgefa<allow_pivot>(m_jac, m_N, m_pivot, IER);
      }

      if (IER != 0) {
        IERPJ = 1;
      }
    }


    /**
     * \brief Adjusts the YH array on reduction of order, and also when the
     * order is increased.
     *
     * \param IORD is the integer flag used to indicate an order increase (IORD
     * = +1) or an order decrease (IORD = -1).
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void dvjust (int IORD)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::dvjust()");
#endif

      // This subroutine adjusts the YH array on reduction of order,
      // and also when the order is increased.

      //  IORD  = An integer flag used to indicate an order
      //          increase (IORD = +1) or an order decrease (IORD = -1).

      if ((m_NQ == 2) && (IORD != 1)) {
        return;
      }

      int NQM1 = m_NQ - 1;
      int NQM2 = m_NQ - 2;

      amrex::Real HSUM{}, XI{}, ALPH0{}, ALPH1{}, PROD{}, XIOLD{}, T1{};

      // Check to see if the order is being increased or decreased.

      if (IORD != 1) {

        // Order decrease.
        for (int j(1); j <= VODE_LMAX; ++j) {
          m_el(j) = 0.0;
        }
        m_el(3) = 1.0;
        HSUM = 0.0;
        for (int j(1); j <= NQM2; ++j) {
          // Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
          HSUM += m_tau(j);
          XI = HSUM / m_HSCAL;
          for (int iback(1); iback <= j+1; ++iback) {
            int i = (j + 4) - iback;
            m_el(i) = m_el(i) * XI + m_el(i-1);
          }
        }

        // Subtract correction terms from YH array.
        for (int j(3); j <= m_NQ; ++j) {
          for (uint64_t i(1); i <= m_N; ++i) {
            m_yh(i,j) -= m_yh(i,m_L) * m_el(j);
          }
        }

      }
      else {

        // Order increase.
        for (int j(1); j <= VODE_LMAX; ++j) {
          m_el(j) = 0.0;
        }

        m_el(3) = 1.0;
        ALPH0 = -1.0;
        ALPH1 = 1.0;
        PROD = 1.0;
        XIOLD = 1.0;
        HSUM = m_HSCAL;

        if (m_NQ != 1) {
          for (int j(1); j <= NQM1; ++j) {
            // Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
            HSUM += m_tau(j+1);
            XI = HSUM / m_HSCAL;
            PROD *= XI;
            ALPH0 -= 1.0 / (j + 1);
            ALPH1 += 1.0 / XI;
            for (int iback(1); iback <= j+1; ++iback) {
              int i = (j + 4) - iback;
              m_el(i) = m_el(i) * XIOLD + m_el(i-1);
            }
            XIOLD = XI;
          }
        }

        T1 = (-ALPH0 - ALPH1) / PROD;
        // Load column L + 1 in YH array.
        for (uint64_t i(1); i <= m_N; ++i) {
          m_yh(i,m_L+1) = T1 * m_yh(i,VODE_LMAX);
        }

        // Add correction terms to YH array.
        for (int j(3); j <= m_NQ + 1; ++j) {
          for (uint64_t i(1); i <= m_N; ++i) {
            m_yh(i,j) += m_el(j) * m_yh(i,m_L+1);
          }
        }
      }
    }


    /**
     * \brief Nonlinear system solver that uses a chord (modified Newton) method
     * with direct linear algebraic system solvers. It then handles the
     * corrector phase of this integration package.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real dvnlsd (int& NFLAG)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::dvnlsd()");
#endif

      amrex::Real ACNRM = 1.e10;

      // dvnlsd is a nonlinear system solver that uses a chord (modified
      // Newton) method with direct linear algebraic system solvers.
      // It then handles the corrector phase of this integration package.

      // Parameter declarations
      constexpr amrex::Real CCMAX = 0.3e0;
      constexpr amrex::Real CRDOWN = 0.3e0;
      constexpr amrex::Real RDIV  = 2.0e0;
      constexpr int MAXCOR = 3;
      constexpr int MSBP = 20;

      amrex::Real DEL{};
      int M{};

      // On the first step, on a change of method order, or after a
      // nonlinear convergence failure with NFLAG = -2, set IPUP = 1
      // to force a Jacobian update.

      if (NFLAG == 0) {
        m_ICF = 0;
      }
      else if (NFLAG == -2) {
        m_IPUP = 1;
      }

      // RC is the ratio of new to old values of the coefficient H / EL(2) = h / l1.
      // When RC differs from 1 by more than CCMAX, IPUP is set to 1
      // to force DVJAC to be called, if a Jacobian is involved.
      // In any case, DVJAC is called at least every MSBP steps.

      m_DRC = std::abs(m_RC - 1.0);
      if (m_DRC > CCMAX || m_n_step >= m_NSLP + MSBP) {
        m_IPUP = 1;
      }

      // Up to MAXCOR corrector iterations are taken.  A convergence test is
      // made on the r.m.s. norm of each correction, weighted by the error
      // weight array ewt.  The sum of the corrections is accumulated in the
      // array acor.  The YH array is not altered in the corrector loop.

      bool converged = false;

      while (true) {

        M = 0;
        amrex::Real DELP = 0.0;

        for (uint64_t i(1); i <= m_N; ++i) {
          m_y(i) = m_yh(i,1);
        }

        rhs(m_tn, m_y, m_savf);
        m_n_rhs += 1;

        if (m_IPUP == 1) {

          // If indicated, the matrix P = I - h*rl1*J is reevaluated and
          // preprocessed before starting the corrector iteration.  IPUP is set
          // to 0 as an indicator that this has been done.

          int IERPJ{};
          dvjac(IERPJ);

          m_IPUP = 0;
          m_RC = 1.0;
          m_DRC = 0.0;
          m_CRATE = 1.0;
          m_NSLP = m_n_step;

          // If matrix is singular, take error return to force cut in step size.
          if (IERPJ != 0) {
            NFLAG = -1;
            m_ICF = 2;
            m_IPUP = 1;
            return ACNRM;
          }

        }

        for (uint64_t i(1); i <= m_N; ++i) {
          m_acor(i) = 0.0;
        }

        // Corrector iteration loop.
        while (true) {

          // Compute the corrector error, and solve the linear system with
          // that as right-hand side and  P as coefficient matrix. The
          // correction is scaled by the factor 2/(1+RC) to account for
          // changes in h*rl1 since the last dvjac call.

          for (uint64_t i(1); i <= m_N; ++i) {
            m_y(i) = (m_RL1 * m_H) * m_savf(i) -
                     (m_RL1 * m_yh(i,2) + m_acor(i));
          }

          if (m_linalg_do_pivoting == 1) {
            constexpr bool allow_pivot{true};
            dgesl<allow_pivot>(m_jac, m_N, m_pivot, m_y);
          } else {
            constexpr bool allow_pivot{false};
            dgesl<allow_pivot>(m_jac, m_N, m_pivot, m_y);
          }

          if (m_RC != 1.0) {
            amrex::Real CSCALE = 2.0 / (1.0 + m_RC);
            for (uint64_t i(1); i <= m_N; ++i) {
              m_y(i) *= CSCALE;
            }
          }

          DEL = 0.0;
          for (uint64_t i(1); i <= m_N; ++i) {
            DEL += (m_y(i) * m_ewt(i)) * (m_y(i) * m_ewt(i));

            m_acor(i) += m_y(i);

            m_y(i) = m_yh(i,1) + m_acor(i);
          }
          DEL = std::sqrt(DEL / m_N);

          // Test for convergence. If M > 0, an estimate of the convergence
          // rate constant is stored in CRATE, and this is used in the test.

          if (M != 0) {
            m_CRATE = amrex::max(CRDOWN * m_CRATE, DEL / DELP);
          }

          amrex::Real DCON = DEL * amrex::min(1.0, m_CRATE) / m_tq(4);
          if (DCON <= 1.0) {
            // we converged, exit the outer loop
            converged = true;
            break;
          }

          M += 1;
          if (M == MAXCOR) {
            //! exit the inner correction iteration
            break;
          }

          if (M >= 2 && DEL > RDIV * DELP) {
            // exit the inner correction iteration
            break;
          }

          DELP = DEL;
          rhs(m_tn, m_y, m_savf);
          m_n_rhs += 1;

        }

        if (converged) {
          break;
        }

        if (m_JCUR == 1) {
          NFLAG = -1;
          m_ICF = 2;
          m_IPUP = 1;
          return ACNRM;
        }

        m_ICF = 1;
        m_IPUP = 1;
      }

      // Return for successful step.
      NFLAG = 0;
      m_JCUR = 0;
      m_ICF = 0;
      if (M == 0) {
        ACNRM = DEL;
      }

      if (M > 0) {
        ACNRM = 0.0;
        for (uint64_t i(1); i <= m_N; ++i) {
          ACNRM += (m_acor(i) * m_ewt(i)) * (m_acor(i) * m_ewt(i));
        }
        ACNRM = std::sqrt(ACNRM / m_N);
      }

      return ACNRM;
    }


    /**
     * \brief The VODE integration algorithm
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void dvode ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::dvode()");
#endif

      // Local variables
      amrex::Real H0{}, S{};
      int IER{}, NITER{};

      // Flag determining if we were successful.
      m_integration_status = IntegrationStatus::IERR_SUCCESS;

      // Set the (inverse of the) timestep limiter.
      m_HMXI = 1.0 / m_ode_max_dt;

      // Return if the final time matches the starting time.
      if (m_tout == m_time) {
        return;
      }

      // All remaining initializations, the initial call to F,
      // and the calculation of the initial step size.
      // The error weights in m_ewt are inverted after being loaded.
      m_tn = m_time;

      m_n_step = 0;
      m_n_jac = 0;
      m_NSLJ = 0;

      // Initial call to the RHS.

      rhs(m_time, m_y, m_y_rhs);

      m_n_rhs = 1;

      // H is temporarily set to 1.0
      m_NQ = 1;
      m_H = 1.0;

      for (uint64_t i(1); i <= m_N; ++i) {
        m_yh(i,2) = m_y_rhs(i);

        // Load the initial value array in yh.
        m_yh(i,1) = m_y(i);

        // Load and invert the ewt array.
        m_ewt(i) = m_rtol * std::abs(m_yh(i,1)) + m_atol;
        m_ewt(i) = 1.0 / m_ewt(i);
      }

      // Call DVHIN to set initial step size H0 to be attempted.
      H0 = 0.0;
      dvhin(H0, NITER, IER);
      m_n_rhs += NITER;

      if (IER != 0) {
#ifndef AMREX_USE_GPU
        std::cout << "DVODE: TOUT too close to T to start integration" << std::endl;
#endif
        m_integration_status = IntegrationStatus::IERR_SPRAD_CONVERGENCE;
        return;
      }

      // Load H with H0 and scale yh(:,2) by H0.
      m_H = H0;
      for (uint64_t i(1); i <= m_N; ++i) {
        m_yh(i,2) *= H0;
      }

      // Start with the order set to 1, and initialize other variables.  ETAMAX is
      // the maximum ratio by which H can be increased in a single step. It is
      // normally 10, but is larger during the first step to compensate for the
      // small initial H. If a failure occurs (in corrector convergence or error
      // test), ETAMAX is set to 1 for the next increase.

      m_NQ = 1;
      m_NEWQ = 1;
      m_L = 2;
      m_tau(1) = m_H;
      m_PRL1 = 1.0;
      m_RC = 0.0;
      m_ETAMAX = 1.0e4;
      m_NQWAIT = 2;
      m_HSCAL = m_H;
      m_NEWH = 0;
      m_NSLP = 0;
      m_IPUP = 1;

      bool skip_loop_start = true;

      // Now do the actual integration as a loop over dvstep.

      while (m_integration_status == IntegrationStatus::IERR_SUCCESS) {

        if (!skip_loop_start) {

          // First check for too many steps being taken, update ewt (if not at
          // start of problem), check for too much accuracy being requested, and
          // check for H below the roundoff level in T.

          if (m_n_step >= m_ode_max_steps) {
            // The maximum number of steps was taken before reaching TOUT.
#ifndef AMREX_USE_GPU
            std::cout << amrex::Font::Bold << amrex::FGColor::Red
                      << "DVODE: maximum number of steps taken before reaching TOUT"
                      << amrex::ResetDisplay << std::endl;
#endif
            for (uint64_t i(1); i <= m_N; ++i) {
              m_y(i) = m_yh(i,1);
            }

            m_time = m_tn;

            m_integration_status = IntegrationStatus::IERR_TOO_MANY_STEPS;
            return;
          }

          for (uint64_t i(1); i <= m_N; ++i) {
            m_ewt(i) = m_rtol * std::abs(m_yh(i,1)) + m_atol;
            m_ewt(i) = 1.0 / m_ewt(i);
          }

        }
        else {
          skip_loop_start = false;
        }

        amrex::Real TOLSF = 0.0;
        for (uint64_t i(1); i <= m_N; ++i) {
          TOLSF += (m_yh(i,1) * m_ewt(i)) * (m_yh(i,1) * m_ewt(i));
        }
        TOLSF = UROUND * std::sqrt(TOLSF / m_N);

        if (TOLSF > 1.0) {

          if (m_n_step == 0) {
#ifndef AMREX_USE_GPU
            std::cout << amrex::Font::Bold << amrex::FGColor::Red
                      << "DVODE: too much accuracy requested at start of integration"
                      << amrex::ResetDisplay << std::endl;
#endif
            m_integration_status = IntegrationStatus::IERR_TOO_MUCH_ACCURACY_REQUESTED;
            return;
          }

          // Too much accuracy requested for machine precision.
#ifndef AMREX_USE_GPU
          std::cout << "DVODE: too much accuracy requested" << std::endl;
#endif
          for (uint64_t i(1); i <= m_N; ++i) {
            m_y(i) = m_yh(i,1);
          }

          m_time = m_tn;

          m_integration_status = IntegrationStatus::IERR_TOO_MUCH_ACCURACY_REQUESTED;
          return;
        }

        int kflag = dvstep();

        // Branch on KFLAG. KFLAG can be 0, -1, or -2.
        if (kflag == -1) {
          // Error test failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
          std::cout << amrex::Font::Bold << amrex::FGColor::Red
                    << "DVODE: error test failed repeatedly or with abs(H) = HMIN"
                    << amrex::ResetDisplay << std::endl;
#endif
          // Set Y array, T, and optional output.
          for (uint64_t i(1); i <= m_N; ++i) {
            m_y(i) = m_yh(i,1);
          }

          m_time = m_tn;

          m_integration_status = IntegrationStatus::IERR_DT_UNDERFLOW;
          return;

        }
        else if (kflag == -2) {
          // Convergence failed repeatedly or with ABS(H) = HMIN.
#ifndef AMREX_USE_GPU
          std::cout << amrex::Font::Bold << amrex::FGColor::Red
                    << "DVODE: corrector convergence failed repeatedly or with abs(H) = HMIN"
                    << amrex::ResetDisplay << std::endl;
#endif
          // Set Y array, T, and optional output.
          for (uint64_t i(1); i <= m_N; ++i) {
            m_y(i) = m_yh(i,1);
          }

          m_time = m_tn;

          m_integration_status = IntegrationStatus::IERR_CORRECTOR_CONVERGENCE;
          return;

        }

        // Otherwise, we've had a successful return from the integrator (kflag = 0).
        // Test for our stopping condition.

        if ((m_tn - m_tout) * m_H < 0.0)
          continue;

        // If TOUT has been reached, interpolate.
        for (uint64_t i(1); i <= m_N; ++i) {
          m_y(i) = m_yh(i,m_L);
        }

        S = (m_tout - m_tn) / m_H;

        for (int jb(1); jb <= m_NQ; ++jb) {
          int j = m_NQ - jb;
          for (uint64_t i(1); i <= m_N; ++i) {
            m_y(i) = m_yh(i,j+1) + S * m_y(i);
          }
        }

        m_time = m_tout;

        m_integration_status = IntegrationStatus::IERR_SUCCESS;
        return;
      }
    }


    /**
     * \brief Sets coefficients for use in VODE::dvstep()
     *
     * For each order m_NQ, the coefficients in EL are calculated by use of the
     * generating polynomial lambda(x), with coefficients EL(i).
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void dvset ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::dvset()");
#endif

      // dvset is called by dvstep and sets coefficients for use there.
      //
      // For each order m_NQ, the coefficients in EL are calculated by use of
      // the generating polynomial lambda(x), with coefficients EL(i).
      //      lambda(x) = EL(1) + EL(2)*x + ... + EL(m_NQ+1)*(x**NQ).
      // Since we are using the BDF method, this is
      //                                     NQ-1
      //     lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
      //                                     i = 1
      // xi(i) is defined by
      // H*xi(i) = t sub n  -  t sub (n-i)
      //         = H + TAU(1) + TAU(2) + ... TAU(i-1).

      constexpr amrex::Real CORTES = 0.1e0;

      amrex::Real FLOTL = m_L;
      int NQM1 = m_NQ - 1;
      int NQM2 = m_NQ - 2;

      for (int i(3); i <= m_L; ++i) {
        m_el(i) = 0.0;
      }

      m_el(1) = 1.0;
      m_el(2) = 1.0;
      amrex::Real ALPH0 = -1.0;
      amrex::Real AHATN0 = -1.0;
      amrex::Real HSUM = m_H;
      amrex::Real RXI = 1.0;
      amrex::Real RXIS = 1.0;

      if (m_NQ != 1) {
        for (int j(1); j <= NQM2; ++j) {
          // In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)).
          HSUM += m_tau(j);
          RXI = m_H / HSUM;
          ALPH0 -= 1.0 / (j+1);
          for (int iback(1); iback <= j+1; ++iback) {
            const int i = (j + 3) - iback;
            m_el(i) += m_el(i-1) * RXI;
          }
        }

        ALPH0 -= 1.0 / m_NQ;
        RXIS = -m_el(2) - ALPH0;
        HSUM += m_tau(NQM1);
        RXI = m_H / HSUM;
        AHATN0 = -m_el(2) - RXI;
        for (int iback(1); iback <= m_NQ; ++iback) {
          const int i = (m_NQ + 2) - iback;
          m_el(i) += m_el(i-1) * RXIS;
        }
      }

      const amrex::Real T1 = 1.0 - AHATN0 + ALPH0;
      const amrex::Real T2 = 1.0 + m_NQ * T1;
      m_tq(2) = std::abs(ALPH0 * T2 / T1);
      m_tq(5) = std::abs(T2 / (m_el(m_L) * RXI / RXIS));

      if (m_NQWAIT == 1) {
        const amrex::Real CNQM1 = RXIS / m_el(m_L);
        const amrex::Real T3 = ALPH0 + 1.0 / m_NQ;
        const amrex::Real T4 = AHATN0 + RXI;
        amrex::Real ELP = T3 / (1.0 - T4 + T3);
        m_tq(1) = std::abs(ELP / CNQM1);
        HSUM += m_tau(m_NQ);
        RXI = m_H / HSUM;
        const amrex::Real T5 = ALPH0 - 1.0 / (m_NQ+1);
        const amrex::Real T6 = AHATN0 - RXI;
        ELP = T2 / (1.0 - T6 + T5);
        m_tq(3) = std::abs(ELP * RXI * (FLOTL + 1.0) * T5);
      }

      m_tq(4) = CORTES * m_tq(2);
    }


    /**
     * \brief Effectively multiplies the Nordsieck history array by the Pascal
     * triangle matrix.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void advance_nordsieck ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::advance_nordsieck()");
#endif

      // Effectively multiplies the Nordsieck history
      // array by the Pascal triangle matrix.

      for (int k(m_NQ); k >= 1; --k) {
        for (int j(k); j <= m_NQ; ++j) {
          for (uint64_t i(1); i <= m_N; ++i) {
            m_yh(i,j) += m_yh(i,j+1);
          }
        }
      }
    }


    /**
     * \brief Undoes the Pascal triangle matrix multiplication implemented in
     * subroutine advance_nordsieck.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void retract_nordsieck ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::retract_nordsieck()");
#endif

      // Undoes the Pascal triangle matrix multiplication
      // implemented in subroutine advance_nordsieck.

      for (int k(m_NQ); k >= 1; --k) {
        for (int j(k); j <= m_NQ; ++j) {
          for (uint64_t i(1); i <= m_N; ++i) {
            m_yh(i,j) -= m_yh(i,j+1);
          }
        }
      }
    }


    /**
     * \brief performs one step of the integration of an initial value problem
     * for a system of ordinary differential equations. calls subroutine dvnlsd
     * for the solution of the nonlinear system arising in the time step.  Thus
     * it is independent of the problem Jacobian structure and the type of
     * nonlinear system solution method. returns a completion flag kflag.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int dvstep ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("VODE::dvstep()");
#endif

      // dvstep performs one step of the integration of an initial value
      // problem for a system of ordinary differential equations.
      // dvstep calls subroutine dvnlsd for the solution of the nonlinear system
      // arising in the time step.  Thus it is independent of the problem
      // Jacobian structure and the type of nonlinear system solution method.
      // dvstep returns a completion flag kflag.

      // A return with kflag = -1 or -2 means either abs(H) = HMIN or 10
      // consecutive failures occurred. On a return with kflag negative,
      // the values of TN and the m_yh array are as of the beginning of the last
      // step, and m_H is the last step size attempted.

      // On a successful return, ETAMAX is reset and ACOR is scaled.

      // Parameter declarations
      constexpr int KFC = -3;
      constexpr int KFH = -7;
      constexpr int MXNCF = 10;
      constexpr amrex::Real ADDON = 1.0e-6;
      constexpr amrex::Real BIAS1 = 6.0e0;
      constexpr amrex::Real BIAS2 = 6.0e0;
      constexpr amrex::Real BIAS3 = 10.0e0;
      constexpr amrex::Real ETACF = 0.25e0;
      constexpr amrex::Real ETAMIN = 0.1e0;
      constexpr amrex::Real ETAMXF = 0.2e0;
      constexpr amrex::Real ETAMX2 = 10.0e0;
      constexpr amrex::Real ETAMX3 = 10.0e0;
      constexpr amrex::Real ONEPSM = 1.00001e0;
      constexpr amrex::Real THRESH = 1.5e0;

      amrex::Real CNQUOT{}, DDN{}, DSM{}, DUP{}, TOLD{};
      amrex::Real FLOTL{}, R{};
      int NCF{}, NFLAG{};

      int kflag = 0;
      TOLD = m_tn;
      NCF = 0;
      m_JCUR = 0;
      NFLAG = 0;

      // Take preliminary actions if the driver changed H. ETA must be
      // reset and the history array yh is rescaled.

      // If a change of order was dictated on the previous step, then
      // it is done here and appropriate adjustments in the history are made.

      if (m_NEWH != 0) {

        if (m_NEWQ < m_NQ) {
          dvjust(-1);
          m_NQ = m_NEWQ;
          m_L = static_cast<short>(m_NQ + 1);
          m_NQWAIT = m_L;
        }
        else if (m_NEWQ > m_NQ) {
          dvjust(1);
          m_NQ = m_NEWQ;
          m_L = static_cast<short>(m_NQ + 1);
          m_NQWAIT = m_L;
        }

        // Rescale the history array for a change in H by a factor of ETA.

        R = 1.0;

        for (int j(2); j <= m_L; ++j) {
          R *= m_ETA;
          for (uint64_t i(1); i <= m_N; ++i) {
            m_yh(i,j) *= R;
          }
        }

        m_H = m_HSCAL * m_ETA;
        m_HSCAL = m_H;
        m_RC = m_RC * m_ETA;

      }

      // Compute the predicted values by effectively
      // multiplying the yh array by the Pascal triangle matrix.
      // dvset is called to calculate all integration coefficients.
      // RC is the ratio of new to old values of the coefficient H / EL(2) = h / l1.

      // Save the value of y before calling the solver.
      // We will use this saved value to determine if any
      // constraints were violated in the update.

      for (uint64_t n(1); n <= m_N; ++n)
        m_y_old(n) = m_y(n);

      while (true) {

        m_tn += m_H;

        advance_nordsieck();

        dvset();

        m_RL1 = 1.0 / m_el(2);
        m_RC *= (m_RL1 / m_PRL1);
        m_PRL1 = m_RL1;


        // Call the nonlinear system solver.

        const amrex::Real ACNRM = dvnlsd(NFLAG);

        if (NFLAG != 0) {

          // The VNLS routine failed to achieve convergence (NFLAG != 0).
          // The yh array is retracted to its values before prediction.
          // The step size H is reduced and the step is retried, if possible.
          // Otherwise, an error exit is taken.

          NCF += 1;
          m_ETAMAX = 1.0;
          m_tn = TOLD;

          retract_nordsieck();

          if (std::abs(m_H) <= m_HMIN * ONEPSM) {
            kflag = -2;
            return kflag;
          }

          if (NCF == MXNCF) {
            kflag = -2;
            return kflag;
          }

          m_ETA = ETACF;
          m_ETA = amrex::max(m_ETA, m_HMIN / std::abs(m_H));

          // Rescale the history array for a change in H by a factor of ETA.
          R = 1.0;

          for (int j(2); j <= m_L; ++j) {
            R *= m_ETA;
            for (uint64_t i(1); i <= m_N; ++i) {
              m_yh(i,j) *= R;
            }
          }

          m_H = m_HSCAL * m_ETA;
          m_HSCAL = m_H;
          m_RC *= m_ETA;

          continue;
        }

        // We add to VODE some constraints. If these constraints are violated,
        // we treat it the same way we will treat an error test failure (below).

        m_integration_status = IntegrationStatus::IERR_SUCCESS;

        {
          auto const& indexes = m_mem_handler->indexes();

          if (indexes.fluid_indexes().count() > 0)
          {
            constexpr uint64_t species_begin = ChemMemoryIndexes::fluid_species_begin();
            uint64_t const species_end = indexes.fluid_species_end();

            amrex::Real X_tot(0.);
            amrex::Real X_tot_save(0.);
            for (uint64_t n(species_begin); n < species_end; ++n) {
              X_tot += m_y[n];
              X_tot_save += m_y_old[n];
            }

            for (uint64_t n(species_begin); n < species_end; ++n) {
              amrex::Real X_gk = m_y[n] / X_tot;
              amrex::Real X_gk_save = m_y_old[n] / X_tot_save;

              // Constrain abundances such that they don't change by more than
              // a certain factor in  a given step (as long as their initial
              // abundance was large enough for us to worry about)

              if (std::abs(X_gk_save) > m_X_reject_buffer * m_atol &&
                  std::abs(X_gk) > m_X_reject_buffer * m_atol &&
                  (std::abs(X_gk) > m_vode_increase_change_factor * std::abs(X_gk_save) ||
                   std::abs(X_gk) < m_vode_decrease_change_factor * std::abs(X_gk_save))) {
#ifdef AMREX_DEBUG
#ifndef AMREX_USE_GPU
                std::cout << "rejecting step based on fluid species " << n+1 << " from "
                          << X_gk_save << " to " << X_gk << std::endl;
#endif
#endif
                m_integration_status = IntegrationStatus::IERR_INVALID_MASS_FRACTIONS;
                break;
              }

              // Constrain abundances such that they are not negative (within
              // a tolerance) or greater than one (within a tolerance)

              if (X_gk < -m_species_failure_tolerance) {
                m_integration_status = IntegrationStatus::IERR_INVALID_MASS_FRACTIONS;
                break;
              }

              if (X_gk > 1.0 + m_species_failure_tolerance) {
                m_integration_status = IntegrationStatus::IERR_INVALID_MASS_FRACTIONS;
                break;
              }
            }
          }

          if ((indexes.solids_indexes().count() > 0) && (m_N >= indexes.count()))
          {
            const uint64_t species_begin = indexes.solids_species_begin();
            uint64_t const species_end = indexes.solids_species_end();

            amrex::Real X_tot(0.);
            amrex::Real X_tot_save(0.);
            for (uint64_t n(species_begin); n < species_end; ++n) {
              X_tot += m_y[n];
              X_tot_save += m_y_old[n];
            }

            if (m_y[indexes.solids_density()] > 0 &&
                m_y[indexes.solids_radius()] > 0 &&
                !amrex::almostEqual(m_y[indexes.solids_density()], 0.) &&
                !amrex::almostEqual(m_y[indexes.solids_radius()], 0.)) {

              for (uint64_t n(species_begin); n < species_end; ++n) {
                amrex::Real X_sn = m_y[n] / X_tot;
                amrex::Real X_sn_save = m_y_old[n] / X_tot_save;

                // Constrain abundances such that they don't change by more than
                // a certain factor in  a given step (as long as their initial
                // abundance was large enough for us to worry about)

                if (std::abs(X_sn_save) > m_X_reject_buffer * m_atol &&
                    std::abs(X_sn) > m_X_reject_buffer * m_atol &&
                    (std::abs(X_sn) > m_vode_increase_change_factor * std::abs(X_sn_save) ||
                     std::abs(X_sn) < m_vode_decrease_change_factor * std::abs(X_sn_save))) {
#ifdef AMREX_DEBUG
#ifndef AMREX_USE_GPU
                  std::cout << "rejecting step based on solids species " << n+1 << " from "
                            << X_sn_save << " to " << X_sn << std::endl;
#endif
#endif
                  m_integration_status = IntegrationStatus::IERR_INVALID_MASS_FRACTIONS;
                  break;
                }

                // Constrain abundances such that they are not negative (within
                // a tolerance) or greater than one (within a tolerance)

                if (X_sn < -m_species_failure_tolerance) {
                  m_integration_status = IntegrationStatus::IERR_INVALID_MASS_FRACTIONS;
                  break;
                }

                if (X_sn > 1.0 + m_species_failure_tolerance) {
                  m_integration_status = IntegrationStatus::IERR_INVALID_MASS_FRACTIONS;
                  break;
                }
              }
            } else {
              m_integration_status = IntegrationStatus::IERR_INVALID_SOLIDS_MASS;
              break;
            }
          }
        }

        // The corrector has converged (NFLAG = 0). The local error test is
        // made and control passes to statement 500 if it fails.

        DSM = ACNRM / m_tq(2);
        if ((DSM <= 1.0) && (m_integration_status == IntegrationStatus::IERR_SUCCESS))
        {

          // After a successful step, update the yh and TAU arrays and decrement
          // NQWAIT. If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
          // for use in a possible order increase on the next step.
          // If ETAMAX = 1 (a failure occurred this step), keep NQWAIT >= 2.

          kflag = 0;
          m_n_step += 1;
          for (int iback(1); iback <= m_NQ; ++iback) {
            int i = m_L - iback;
            m_tau(i+1) = m_tau(i);
          }

          m_tau(1) = m_H;
          for (int j(1); j <= m_L; ++j) {
            for (uint64_t i(1); i <= m_N; ++i) {
              m_yh(i,j) += m_el(j) * m_acor(i);
            }
          }

          m_NQWAIT -= 1;
          if ((m_L != VODE_LMAX) && (m_NQWAIT == 1)) {
            for (uint64_t i(1); i <= m_N; ++i) {
              m_yh(i,VODE_LMAX) = m_acor(i);
            }
            m_CONP = m_tq(5);
          }

          if (m_ETAMAX != 1.0) {
            break;
          }

          if (m_NQWAIT < 2) {
            m_NQWAIT = 2;
          }

          m_NEWQ = m_NQ;
          m_NEWH = 0;
          m_ETA = 1.0;
          m_ETAMAX = ETAMX3;
          if (m_n_step <= 10) {
            m_ETAMAX = ETAMX2;
          }

          R = 1.0 / m_tq(2);
          for (uint64_t i(1); i <= m_N; ++i) {
            m_acor(i) *= R;
          }

          return kflag;
        }

        // The error test failed (or our constraints on the species failed).
        // kflag keeps track of multiple failures.
        // Restore m_tn and the yh array to their previous values, and prepare
        // to try the step again. Compute the optimum step size for the
        // same order. After repeated failures, H is forced to decrease
        // more rapidly.

        kflag -= 1;
        NFLAG = -2;
        m_tn = TOLD;

        retract_nordsieck();

        if (std::abs(m_H) <= m_HMIN * ONEPSM) {
          kflag = -1;
          return kflag;
        }

        m_ETAMAX = 1.0;

        if (kflag > KFC) {

          // Compute ratio of new H to current H at the current order.
          FLOTL = m_L;
          m_ETA = 1.0 / (std::pow(BIAS2 * DSM, 1.0 / FLOTL) + ADDON);
          m_ETA = amrex::max(m_ETA, m_HMIN / std::abs(m_H), ETAMIN);
          if ((kflag <= -2) && (m_ETA > ETAMXF)) {
            m_ETA = ETAMXF;
          }

          // Rescale the history array for a change in H by a factor of ETA.
          R = 1.0;

          for (int j(2); j <= m_L; ++j) {
            R *= m_ETA;
            for (uint64_t i(1); i <= m_N; ++i) {
              m_yh(i,j) *= R;
            }
          }

          m_H = m_HSCAL * m_ETA;
          m_HSCAL = m_H;
          m_RC = m_RC * m_ETA;

          continue;

        }

        // Control reaches this section if 3 or more consecutive failures
        // have occurred. It is assumed that the elements of the yh array
        // have accumulated errors of the wrong order. The order is reduced
        // by one, if possible. Then H is reduced by a factor of 0.1 and
        // the step is retried.  After a total of 7 consecutive failures,
        // an exit is taken with kflag = -1.

        if (kflag == KFH) {
          kflag = -1;
          return kflag;
        }

        if (m_NQ != 1) {
          m_ETA = amrex::max(ETAMIN, m_HMIN / std::abs(m_H));
          dvjust(-1);
          m_L = m_NQ;
          m_NQ -= 1;
          m_NQWAIT = m_L;

          // Rescale the history array for a change in H by a factor of ETA.
          R = 1.0;

          for (int j(2); j <= m_L; ++j) {
            R *= m_ETA;
            for (uint64_t i(1); i <= m_N; ++i) {
              m_yh(i,j) *= R;
            }
          }

          m_H = m_HSCAL * m_ETA;
          m_HSCAL = m_H;
          m_RC = m_RC * m_ETA;

          continue;
        }

        m_ETA = amrex::max(ETAMIN, m_HMIN / std::abs(m_H));
        m_H *= m_ETA;
        m_HSCAL = m_H;
        m_tau(1) = m_H;
        rhs(m_tn, m_y, m_savf);
        m_n_rhs += 1;
        for (uint64_t i(1); i <= m_N; ++i) {
          m_yh(i,2) = m_H * m_savf(i);
        }

        m_NQWAIT = 10;
      }

      // If NQWAIT = 0, an increase or decrease in order by one is considered.
      // Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
      // be multiplied at order q, q-1, or q+1, respectively.
      // The largest of these is determined, and the new order and
      // step size set accordingly.
      // A change of H or NQ is made only if H increases by at least a
      // factor of THRESH.  If an order change is considered and rejected,
      // then NQWAIT is set to 2 (reconsider it after 2 steps).

      bool already_set_eta = false;

      // Compute ratio of new H to current H at the current order.
      FLOTL = m_L;
      amrex::Real ETAQ = 1.0 / (std::pow(BIAS2 * DSM, 1.0 / FLOTL) + ADDON);
      if (m_NQWAIT == 0) {
        m_NQWAIT = 2;

        amrex::Real ETAQM1{0.0};

        if (m_NQ != 1) {
          // Compute ratio of new H to current H at the current order less one.
          DDN = 0.0;
          for (uint64_t i(1); i <= m_N; ++i) {
            DDN += (m_yh(i,m_L) * m_ewt(i)) * (m_yh(i,m_L) * m_ewt(i));
          }
          DDN = std::sqrt(DDN / m_N) / m_tq(1);
          ETAQM1 = 1.0 / (std::pow(BIAS1 * DDN, 1.0 / (FLOTL - 1.0)) + ADDON);
        }

        amrex::Real ETAQP1{0.0};

        if (m_L != VODE_LMAX) {
          // Compute ratio of new H to current H at current order plus one.
          CNQUOT = (m_tq(5) / m_CONP) * std::pow(m_H / m_tau(2), m_L);
          DUP = 0.0;
          for (uint64_t i(1); i <= m_N; ++i) {
            m_savf(i) = m_acor(i) - CNQUOT * m_yh(i,VODE_LMAX);

            DUP += (m_savf(i) * m_ewt(i)) * (m_savf(i) * m_ewt(i));
          }
          DUP = std::sqrt(DUP / m_N) / m_tq(3);
          ETAQP1 = 1.0 / (std::pow(BIAS3 * DUP, 1.0 / (FLOTL + 1.0)) + ADDON);
        }

        if (ETAQ < ETAQP1) {
          if (ETAQP1 > ETAQM1) {
            m_ETA = ETAQP1;
            m_NEWQ = static_cast<short>(m_NQ + 1);
            for (uint64_t i(1); i <= m_N; ++i) {
              m_yh(i,VODE_LMAX) = m_acor(i);
            }
          }
          else {
            m_ETA = ETAQM1;
            m_NEWQ = static_cast<short>(m_NQ - 1);
          }
          already_set_eta = true;
        }

        if (ETAQ < ETAQM1 && !already_set_eta) {
          m_ETA = ETAQM1;
          m_NEWQ = static_cast<short>(m_NQ - 1);
          already_set_eta = true;
        }

      }

      if (!already_set_eta) {
        m_ETA = ETAQ;
        m_NEWQ = m_NQ;
      }

      // Test tentative new H against THRESH and ETAMAX, and HMXI, then exit.
      if (m_ETA >= THRESH && m_ETAMAX != 1.0) {
        m_ETA = amrex::min(m_ETA, m_ETAMAX);
        m_ETA = m_ETA / amrex::max(1.0, std::abs(m_H) * m_HMXI * m_ETA);
        m_NEWH = 1;
        m_ETAMAX = ETAMX3;
        if (m_n_step <= 10) {
          m_ETAMAX = ETAMX2;
        }
        R = 1.0 / m_tq(2);
        for (uint64_t i(1); i <= m_N; ++i) {
          m_acor(i) *= R;
        }
        return kflag;
      }

      m_NEWQ = m_NQ;
      m_NEWH = 0;
      m_ETA = 1.0;
      m_ETAMAX = ETAMX3;
      if (m_n_step <= 10) {
        m_ETAMAX = ETAMX2;
      }
      R = 1.0 / m_tq(2);
      for (uint64_t i(1); i <= m_N; ++i) {
        m_acor(i) *= R;
      }

      return kflag;
    }


  private:
#ifdef VODE_JACOBIAN_CACHING
    int m_use_jacobian_caching;

    // CCMXJ  = Threshold on m_DRC for updating the Jacobian
    static constexpr amrex::Real m_CCMXJ = 0.2e0;

    // How many timesteps should pass before refreshing the Jacobian
    static constexpr int m_max_steps_between_jacobian_evals = 50;
#endif

    int m_linalg_do_pivoting;

    amrex::Real m_X_reject_buffer;

    static constexpr amrex::Real m_vode_increase_change_factor = 4.0;

    static constexpr amrex::Real m_vode_decrease_change_factor = 0.25;

    static constexpr amrex::Real m_HMIN = 0.0;

    // CONP   = The saved value of TQ(5)
    amrex::Real m_CONP = 0.;

    // CRATE  = Estimated corrector convergence rate constant
    amrex::Real m_CRATE = 0.;

    // Relative change in m_H*m_RL1 since last DVJAC call
    amrex::Real m_DRC = 0.;

    // ETA    = Saved tentative ratio of new to old m_H
    amrex::Real m_ETA = 0.;

    // ETAMAX = Saved maximum value of m_ETA to be allowed
    amrex::Real m_ETAMAX = 0.;

    // H      = The step size
    amrex::Real m_H = 0.;

    // HSCAL  = Stepsize in scaling of YH array
    amrex::Real m_HSCAL = 0.;

    // PRL1   = The saved value of m_RL1
    amrex::Real m_PRL1 = 0.;

    // HMXI   = Inverse of the maximum absolute value of H to be used.
    //          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
    amrex::Real m_HMXI = 0.;

    // RC     = Ratio of current H*RL1 to value on last DVJAC call
    amrex::Real m_RC = 0.;

    // RL1    = The reciprocal of the coefficient EL(1)
    amrex::Real m_RL1 = 0.;

    // TN     = The independent variable, updated on each step taken
    amrex::Real m_tn = 0.;

    // ICF    = Integer flag for convergence failure in DVNLSD:
    //            0 means no failures
    //            1 means convergence failure with out of date Jacobian
    //                   (recoverable error)
    //            2 means convergence failure with current Jacobian or
    //                   singular matrix (unrecoverable error)
    short m_ICF = 0;

    // IPUP   = Saved flag to signal updating of Newton matrix
    short m_IPUP = 0;

    // JCUR   = Output flag from DVJAC showing Jacobian status:
    //            JCUR = 0 means J is not current
    //            JCUR = 1 means J is current
    short m_JCUR = 0;

    // L      = Integer variable, NQ + 1, current order plus one
    short m_L = 0;

    // NEWH   = Saved integer to flag change of H
    short m_NEWH = 0;

    // NEWQ   = The method order to be used on the next step
    short m_NEWQ = 0;

    // NQ     = Integer variable, the current integration method order
    short m_NQ = 0;

    // NQWAIT = A counter controlling the frequency of order changes.
    //          An order change is about to be considered if NQWAIT = 1.
    short m_NQWAIT = 0;

    // NSLJ   = The number of steps taken as of the last Jacobian update
    int m_NSLJ = 0;

    // NSLP   = Saved value of NST as of last Newton matrix update
    int m_NSLP = 0;

    // EL     = amrex::Real array of integration coefficients.  See DVSET
    amrex::Array1D<amrex::Real, 1, VODE_LMAX> m_el;

    // TAU    = amrex::Real vector of past NQ step sizes, length 13
    amrex::Array1D<amrex::Real, 1, VODE_LMAX> m_tau;

    // TQ     = A real vector of length 5 in which DVSET stores constants
    //          used for the convergence test, the error test, and the
    //          selection of H at a new order.
    amrex::Array1D<amrex::Real, 1, 5> m_tq;

    // Jacobian
    ChemMatrix<amrex::Real> m_jac;

#ifdef VODE_JACOBIAN_CACHING
    // Saved Jacobian
    ChemMatrix<amrex::Real> m_jac_save;
#endif

    // the Nordsieck history array
    ChemMatrix<amrex::Real> m_yh;

    ChemVector<amrex::Real> m_ewt;

    ChemVector<amrex::Real> m_savf;

    // Array of size NEQ used for the accumulated corrections on each
    // step, scaled in the output to represent the estimated local
    // error in Y on the last step.  This is the vector e in the
    // description of the error control.  It is defined only on a
    // successful return from DVODE.
    ChemVector<amrex::Real> m_acor;

    ChemVector<short> m_pivot;
};


#endif
