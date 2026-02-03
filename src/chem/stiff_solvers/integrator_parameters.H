#ifndef INTEGRATOR_PARAMETERS_H
#define INTEGRATOR_PARAMETERS_H

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#include <mfix_reporter.H>
#include <mfix_solids.H>
#include <mfix_fix_inputs.H>

using amrex::toLower;

/**
 * \brief Define the integrator type categories to switch between the different
 * available ODE integrators
 */
struct IntegratorType
{
  // The MFIX-Exa native integrator
  enum { ForwardEuler = 0 };

  // The AMReX-Astro/Microphysics integrators labeled here as "StiffSolvers"
  struct StiffSolver {
    // FE = ForwardEuler
    // BE = BackwardEuler
    // VODE = VODE (Variable-coefficient ODE solver)
    enum { FE = 1, BE, VODE };
  };
};


/**
 *
 */
struct JacobianType { enum { Analytical = 0, Numerical, Broyden }; };


/**
 *
 */
struct SolidsUpdateType { enum { ConstantVolume = 0, ConstantDensity, ConstantSpeciesDensities }; };


/**
 * \brief The different types of integration errors that can occur during the
 * ODE integration algorithms
 */
enum struct IntegrationStatus {
  IERR_SUCCESS = 0,
  IERR_BAD_INPUTS = 1,
  IERR_DT_UNDERFLOW = 2,
  IERR_SPRAD_CONVERGENCE = 3,
  IERR_TOO_MANY_STEPS = 4,
  IERR_TOO_MUCH_ACCURACY_REQUESTED = 5,
  IERR_CORRECTOR_CONVERGENCE = 6,
  IERR_LU_DECOMPOSITION_ERROR = 7,
  IERR_NORMALIZATION_FAILURE = 8,
  IERR_INVALID_SOLIDS_MASS = 9,
  IERR_INVALID_MASS_FRACTIONS = 10,
};

static amrex::Vector<std::string> IntegrationErrors {
  "integrator success",
  "bad inputs",
  "dt underflow",
  "sprad convergence",
  "too many steps",
  "too much accuracy requested",
  "corrector convergence",
  "LU decomposition error",
  "normalization failure",
  "invalid solids mass",
  "invalid mass fractions",
};


/**
 * \brief This class contains all the possible parameters for the ODE
 * integrators
 *
 * This class contains only trivially-constructible members, so it can be copied
 * to the GPU memory and then it can be used to define the ODE integrators
 * within a GPU kernel
 */
class IntegratorParms
{
  public:
    /**
     * \brief Constructs the default-valued IntegratorParms class
     */
    AMREX_GPU_HOST_DEVICE
    IntegratorParms ()
      : m_burner_verbose(0)
      , m_ode_max_steps(150000)
      , m_jacobian_type(JacobianType::Numerical)
      , m_atol(1.e-6)
      , m_rtol(1.e-6)
      , m_ode_max_dt(1.e30)
      , m_max_dt_change_factor(1.001)
      , m_max_iter(100)
      , m_rp_rtol(1.e-6)
#ifdef VODE_JACOBIAN_CACHING
      , m_use_jacobian_caching(1)
#endif
      , m_linalg_do_pivoting(1)
      , m_X_reject_buffer(1.0)
      , m_retry_atol(-1)
      , m_retry_rtol(-1)
      , m_retry_swap_jacobian(0)
    {}

    int m_burner_verbose;
    int m_ode_max_steps;
    int m_jacobian_type;
    amrex::Real m_atol;
    amrex::Real m_rtol;
    amrex::Real m_ode_max_dt;

    amrex::Real m_max_dt_change_factor;

    int m_max_iter;
    amrex::Real m_rp_rtol;

#ifdef VODE_JACOBIAN_CACHING
    int m_use_jacobian_caching;
#endif

    int m_linalg_do_pivoting;
    amrex::Real m_X_reject_buffer;

    amrex::Real m_retry_atol;
    amrex::Real m_retry_rtol;
    int m_retry_swap_jacobian;
};


/**
 * \brief This class parses the inputs file for importing ODE integrator user's
 * settings and parameters
 */
class IntegratorInputs
{
  public:
    /**
     * \brief Constructs the IntegratorInputs class
     *
     * The class constructor parses the inputs file and assigns the user's
     * inputs and parameters to the IntegratorParms object that is member of
     * this class
     */
    IntegratorInputs ()
      : m_solids_update_type(SolidsUpdateType::ConstantVolume)
      , m_solids_thresholds()
      , m_params()
    {
      FixInputs fix("Dec. 2024");
      fix.swap<std::string>("mfix.integrator", "chemistry.integrator");

      // Integrator
      amrex::ParmParse pp("chemistry");

      std::string solids_update_type("constant_volume");
      pp.query("solids.update_type", solids_update_type);

      std::string ulower(toLower(solids_update_type));
      if (ulower == "constant_volume")
        m_solids_update_type = SolidsUpdateType::ConstantVolume;
      else if (ulower == "constant_density")
        m_solids_update_type = SolidsUpdateType::ConstantDensity;
      else if (ulower == "constant_species_densities") {
        m_solids_update_type = SolidsUpdateType::ConstantSpeciesDensities;

        std::string density_model("None");
        pp.query("solids.density", density_model);

        if (amrex::toLower(density_model) != "mixture") {
          reporter::Log(reporter::Error, __FILE__, __LINE__)
            << "Solids update type mismatch:\n"
            << "Solids density model should be set to mixture.\n"
            << "Fix the input deck.";
        }
      } else {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Solids update type " << solids_update_type <<" unknown.\nPlease correct the input deck.";
      }

      pp.query("solids.mass_threshold", m_solids_thresholds.mass);
      pp.query("solids.radius_threshold", m_solids_thresholds.radius);
      pp.query("solids.density_threshold", m_solids_thresholds.density);

      // By default we set the MFIX-Exa native ForwardEuler integrator as the
      // ODE integrator
      std::string integrator("ForwardEuler");

      // We query the user's inputs for the ODE integrator
      pp.query("integrator", integrator);
      std::string ilower(toLower(integrator));

      // If the ODE integrator is one of the StiffSolvers
      if (ilower == "stiffsolver::forwardeuler" ||
          ilower == "stiffsolver::backwardeuler" ||
          ilower == "stiffsolver::vode")
      {
        amrex::ParmParse ppIntegrator("chemistry.integrator");

        ppIntegrator.query("burner_verbose", m_params.m_burner_verbose);
        ppIntegrator.query("ode_max_steps", m_params.m_ode_max_steps);
        ppIntegrator.query("atol", m_params.m_atol);
        ppIntegrator.query("rtol", m_params.m_rtol);
        ppIntegrator.query("retry_atol", m_params.m_retry_atol);
        ppIntegrator.query("retry_rtol", m_params.m_retry_rtol);
        ppIntegrator.query("retry_swap_jacobian", m_params.m_retry_swap_jacobian);
        ppIntegrator.query("ode_max_dt", m_params.m_ode_max_dt);

        // If the StiffSolver is ForwardEuler
        if (ilower == "stiffsolver::forwardeuler")
        {
          m_type = IntegratorType::StiffSolver::FE;

          ppIntegrator.query("max_dt_change_factor",
              m_params.m_max_dt_change_factor);
        }

        // If the StiffSolver is BackwardEuler
        if (ilower == "stiffsolver::backwardeuler")
        {
          m_type = IntegratorType::StiffSolver::BE;

          ppIntegrator.query("max_iter", m_params.m_max_iter);
          ppIntegrator.query("linalg_do_pivoting", m_params.m_linalg_do_pivoting);
          ppIntegrator.query("rp_rtol", m_params.m_rp_rtol);

          std::string jacobian("numerical");
          ppIntegrator.query("jacobian_type", jacobian);

          if (toLower(jacobian) == "analytical")
            m_params.m_jacobian_type = JacobianType::Analytical;
          else if (toLower(jacobian) == "numerical")
            m_params.m_jacobian_type = JacobianType::Numerical;
          else if (toLower(jacobian) == "broyden")
            m_params.m_jacobian_type = JacobianType::Broyden;
          else {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Wrong Jacobian type" << jacobian << "\n"
              << "Options are: Analytical, Numerical, Broyden\n"
              << "Please correct the input deck.";
          }
        }

        // If the StiffSolver is VODE
        if (ilower == "stiffsolver::vode")
        {
          m_type = IntegratorType::StiffSolver::VODE;

#ifdef VODE_JACOBIAN_CACHING
          ppIntegrator.query("use_jacobian_caching", m_params.m_use_jacobian_caching);
#endif

          ppIntegrator.query("linalg_do_pivoting", m_params.m_linalg_do_pivoting);
          ppIntegrator.query("X_reject_buffer", m_params.m_X_reject_buffer);

          std::string jacobian("numerical");
          ppIntegrator.query("jacobian_type", jacobian);

          if (toLower(jacobian) == "analytical")
            m_params.m_jacobian_type = JacobianType::Analytical;
          else if (toLower(jacobian) == "numerical")
            m_params.m_jacobian_type = JacobianType::Numerical;
          else {
            reporter::Log(reporter::Error,__FILE__, __LINE__)
              << "Jacobian type " << jacobian << " unknown.\n"
              << "Options are: Analytical, Numerical\n"
              << "Please correct the input deck.";
          }
        }

        // If the ODE integrator is not any of the StiffSolvers, but it is the
        // MFIX-Exa native ForwardEuler
      } else if (ilower == "forwardeuler") {

        m_type = IntegratorType::ForwardEuler;

      } else {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Integrator type " << integrator << " unknown. Options are:\n"
          << "1) ForwardEuler\n"
          << "2) StiffSolver::ForwardEuler\n"
          << "3) StiffSolver::BackwardEuler\n"
          << "4) StiffSolver::VODE\n\n"
          << "Please correct the input deck.";
      }

    }

    int m_type;

    int m_solids_update_type;

    struct solids_thresholds {
      solids_thresholds () : mass(0), radius(0), density(0) {}
      amrex::Real mass;
      amrex::Real radius;
      amrex::Real density;
    } m_solids_thresholds;

    IntegratorParms m_params;
};


#endif
