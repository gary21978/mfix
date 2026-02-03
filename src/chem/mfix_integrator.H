#ifndef MFIX_INTEGRATOR_H
#define MFIX_INTEGRATOR_H

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_Vector.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Gpu.H>

#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_reactions.H>
#include <mfix_solvers.H>
#include <mfix_chem_memory_handler.H>
#include <integrator.H>


/**
 * \brief This class integrates the chemistry ODE problem using the Forward
 * Euler approach
 */
class MFIXIntegrator: public Integrator
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
                                         m_real2D.mem_bytes();}

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
     * \brief Constructs the MFIXIntegrator class
     *
     * \param params is the reference to the IntegratorParms class which contain
     * all the user inputs and integrator settings for the constructor
     */
    AMREX_GPU_HOST_DEVICE
    MFIXIntegrator (const IntegratorParms& params)
      : Integrator(params)
      , m_aux_memory(0)
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
    }

    /**
     * \brief Integrate the ODE system by doing one step of the Forward Euler
     * method
     *
     * \param dt is the timestep size
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void actual_integrator (const amrex::Real& dt)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXIntegrator::actual_integrator()");
#endif

      // compute the rhs
      actual_rhs(m_y, m_y_rhs);

      // update the solution
      for (uint64_t n(1); n <= m_N; ++n)
        m_y(n) += dt*m_y_rhs(n);

      // clean the new state by restricting the species mass fractions into the
      // [0,1] interval and normalizing them
      clean_state(m_y);
    }
};


#endif
