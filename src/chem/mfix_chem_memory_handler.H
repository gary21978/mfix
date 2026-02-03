#ifndef MFIX_CHEM_MEMORY_HANDLER_H
#define MFIX_CHEM_MEMORY_HANDLER_H

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Gpu.H>

#include <mfix_chem_memory_structs.H>


/**
 * \brief Class that allocates and manages the host and device memory to be used
 * during the chemistry update
 *
 * This class allocates, stores, and deallocates all the host and device memory
 * that is needed for the whole chamistry computation, including the chemistry
 * ODE integration. The device memory is allocated either on the GPU global
 * memory or on its shared memory if the total memory request is within the GPU
 * resources limits.
 */
class ChemMemoryHandler
{
  public:
    /**
     * \brief Constructs a ChemMemoryHandler object
     *
     * The constructor assigns the class members for a correct allocation of the
     * host and device memory for the ODE integration. In addition, it allocates
     * the host and device memory and defines a flag to discriminate whether the
     * device will use its global or shared memory. In case the total amount of
     * memory necessary for the ODE integration is less than the total available
     * share dmemory space on the device, the GPU global memory is not allocated
     * as the GPU will use only its shared memory space.
     *
     * \param size is the number of Box cells where chemistry is solved
     * \param nreactions is the maximum number between lagrangian and eulerian
     * chemical reactions
     * \param nspecies_g is the number of fluid species
     * \param nspecies_s is the number of solids species
     * \param solve_fluid is the flag for solving for the fluid phase  in
     * chemical reactions
     * \param solve_solids is the flag for solving for the solids phase in
     * chemical reactions
     * \param lagrangian_reactions is the flag for solving for lagrangian
     * chemical reactions
     * \param integrator is the integrator type chosen by the user to solve the
     * chemistry ODE problem
     */
    template <class T>
    ChemMemoryHandler ([[maybe_unused]] const int size,
                       const int nreactions,
                       const int nspecies_g,
                       const int nspecies_s,
                       const int solve_fluid,
                       const int solve_solids,
                       const int lagrangian_reactions,
                       const T& integrator)
      : m_nblocks(0)
      , m_needed_sm_bytes(0)
      , m_use_shared_mem(1)
      , m_nreactions(nreactions)
      , m_nvariables(0)
      , m_rates_mem(0)
      , m_variables_mem(0)
      , m_rhs_mem(0)
      , m_aux_short1D_mem(0)
      , m_aux_real1D_mem(0)
      , m_aux_real2D_mem(0)
    {
      // Setup the memory indexes to properly allocate and access the memory
      // spaces
      m_indexes.setup_fluid(nspecies_g, solve_fluid);
      m_indexes.setup_solids(nspecies_s, solve_solids, lagrangian_reactions);

      // Get the total number of variables for the ODE to integrate
      m_nvariables = m_indexes.count();

      // Get the auxiliary memory class for the integrator passed as argument.
      // This class is needed to define the amount of auxiliary variables to be
      // allocated
      auto aux_memory = typename T::aux_memory(m_nvariables, integrator.jacobian_type());

      // Get the number of auxiliary variables of type short
      m_naux_short1D = aux_memory.m_short1D.count;

      // Get the number of auxiliary variables of type real
      m_naux_real1D = aux_memory.m_real1D.count;

      // Get the number of auxiliary bidimensional variables of type real
      m_naux_real2D = aux_memory.m_real2D.count;

#ifdef AMREX_USE_GPU
      // Compute the number of blocks of threads needed to run on the device.
      // The number of threads in each block is a member of this class. The cail
      // operation is needed because we need an integer number of blocks and we
      // need to be sure the number of blocks is sufficient to cover the whole
      // vector with size m_size.
      m_nblocks = static_cast<int>(amrex::Math::ceil(amrex::Real(size) / m_nthreads_per_block));

      // Compute the amount of shared memory (in bytes and for a single block of
      // threads) that is necessary for the ODE integration. For each thread in
      // each block of threads we need enough shared memory to store
      // m_nreactions variables for the reactions rates, m_nvariables variables
      // for the ODE solution, and m_nvariables variables for the ODE right-hand
      // side.
      m_needed_sm_bytes = m_nreactions*m_nthreads_per_block*sizeof(amrex::Real); // rates
      m_needed_sm_bytes += m_nvariables*m_nthreads_per_block*sizeof(amrex::Real); // ODE y
      m_needed_sm_bytes += m_nvariables*m_nthreads_per_block*sizeof(amrex::Real); // ODE y_rhs
      m_needed_sm_bytes += m_nvariables*m_nthreads_per_block*sizeof(amrex::Real); // ODE y_old

      // Also, add the amount of shared memory (in bytes and for a single block
      // of threads) that is necessary to allocate the auxiliary variables
      // needed during the ODE itnegration by the chosen ODE integrator
      m_needed_sm_bytes += m_naux_short1D*m_nthreads_per_block*sizeof(short);
      m_needed_sm_bytes += m_naux_real1D*m_nthreads_per_block*sizeof(amrex::Real);
      m_needed_sm_bytes += m_naux_real2D*m_nthreads_per_block*sizeof(amrex::Real);

      // Get the amount of available GPU shared memory on the current hardware
      std::size_t available_sm_per_block = amrex::Gpu::Device::sharedMemPerBlock();

      // Check whether the needed amount of shared memory exceeds the hardware
      // resources limits. If true, use the GPU global memory; if false, use the
      // GPU shared memory
      if (m_needed_sm_bytes > available_sm_per_block) {
        m_use_shared_mem = 0;
        m_needed_sm_bytes = 0;

        // Allocate the GPU global memory for the reactions rates, the ODE
        // variables, and the ODE right-hand side
        m_rates_mem.resize(size*m_nreactions, 0);
        m_variables_mem.resize(size*m_nvariables, 0);
        m_rhs_mem.resize(size*m_nvariables, 0);
        m_y_old_mem.resize(size*m_nvariables, 0);

        // Allocate the GPU global memory for the integrator auxiliary variables
        m_aux_short1D_mem.resize(size*m_naux_short1D, 0);
        m_aux_real1D_mem.resize(size*m_naux_real1D, 0);
        m_aux_real2D_mem.resize(size*m_naux_real2D, 0);
      }
#else
      // When running on the host, set the shared memory flag to false just for
      // sanity.
      m_use_shared_mem = 0;

      // On the host we need to allocate one vector for the reactions rates, one
      // for the ODE solution, and one for the ODE right-hand side.
      m_rates_mem.resize(m_nreactions, 0);
      m_variables_mem.resize(m_nvariables, 0);
      m_rhs_mem.resize(m_nvariables, 0);
      m_y_old_mem.resize(m_nvariables, 0);

      // On the host we also need to allocate one vector for each of the
      // different types of auxiliary variables for the chosen ODE integrator
      m_aux_short1D_mem.resize(m_naux_short1D, 0);
      m_aux_real1D_mem.resize(m_naux_real1D, 0);
      m_aux_real2D_mem.resize(m_naux_real2D, 0);
#endif
    }

    /**
     * \brief Get the number of GPU blocks of threads.
     */
    AMREX_FORCE_INLINE
    int nblocks () const { return m_nblocks; }

    /**
     * \brief Get the number of GPU threads for each block of threads.
     */
    AMREX_FORCE_INLINE
    static constexpr int nthreads_per_block () { return m_nthreads_per_block; }

    /**
     * \brief Get the value of the flag which switches between using or not
     * using the GPU shared memory
     */
    AMREX_FORCE_INLINE
    int use_shared_mem () const { return m_use_shared_mem; }

    /**
     * \brief Get the number of variables for the ODE system
     */
    AMREX_FORCE_INLINE
    int nvariables () const { return m_nvariables; }

    /**
     * \brief Get the amount in bytes of device shared memory needed for each
     * block of threads
     */
    AMREX_FORCE_INLINE
    int needed_sm_bytes () const { return m_needed_sm_bytes; }

    /**
     * \brief Get the size of the total memory allocated.
     */
    AMREX_FORCE_INLINE
    int glob_mem_size () const
    {
      return m_rates_mem.size() + m_variables_mem.size() + m_rhs_mem.size() +
        m_aux_short1D_mem.size() + m_aux_real1D_mem.size() + m_aux_real2D_mem.size();
    }

    /**
     * \brief Get the pointer to the memory space allocated for the reactions
     * rates
     */
    AMREX_FORCE_INLINE
    amrex::Real* rates_mem () { return m_rates_mem.dataPtr(); }

    /**
     * \brief Get the pointer to the memory space allocated for the ODE solution
     */
    AMREX_FORCE_INLINE
    amrex::Real* variables_mem () { return m_variables_mem.dataPtr(); }

    /**
     * \brief Get the pointer to the memory space allocated for the ODE
     * right-hand side
     */
    AMREX_FORCE_INLINE
    amrex::Real* rhs_mem () { return m_rhs_mem.dataPtr(); }

    /**
     * \brief
     */
    AMREX_FORCE_INLINE
    amrex::Real* y_old_mem () { return m_y_old_mem.dataPtr(); }

    /**
     * \brief Get the number of auxiliary variables of type short for the chosen
     * ODE integrator
     */
    AMREX_FORCE_INLINE
    int naux_short1D () const { return m_naux_short1D; }

    /**
     * \brief Get the number of auxiliary variables of type Real for the chosen
     * ODE integrator
     */
    AMREX_FORCE_INLINE
    int naux_real1D () const { return m_naux_real1D; }

    /**
     * \brief Get the number of auxiliary bidimensional variables of type Real
     * for the chosen ODE integrator
     */
    AMREX_FORCE_INLINE
    int naux_real2D () const { return m_naux_real2D; }

    /**
     * \brief Get the pointer to the auxiliary variables of type short for the
     * chosen ODE integrator
     */
    AMREX_FORCE_INLINE
    short* aux_short1D_mem ()
    {
      return (m_naux_short1D > 0)? m_aux_short1D_mem.dataPtr() : nullptr;
    }

    /**
     * \brief Get the pointer to the auxiliary variables of type Real for the
     * chosen ODE integrator
     */
    AMREX_FORCE_INLINE
    amrex::Real* aux_real1D_mem ()
    {
      return (m_naux_real1D > 0)? m_aux_real1D_mem.dataPtr() : nullptr;
    }

    /**
     * \brief Get the pointer to the auxiliary bidimensional variables of type
     * Real for the chosen ODE integrator
     */
    AMREX_FORCE_INLINE
    amrex::Real* aux_real2D_mem ()
    {
      return (m_naux_real2D > 0)? m_aux_real2D_mem.dataPtr() : nullptr;
    }

    /**
     * \brief
     */
    AMREX_FORCE_INLINE
    ChemMemoryIndexes indexes () const { return m_indexes; }

  private:
    static constexpr int m_nthreads_per_block = AMREX_GPU_MAX_THREADS;

    int m_nblocks;
    std::size_t m_needed_sm_bytes;
    int m_use_shared_mem;

    int m_nreactions;
    int m_nvariables;

    ChemMemoryIndexes m_indexes;

    amrex::Gpu::DeviceVector<amrex::Real> m_rates_mem;
    amrex::Gpu::DeviceVector<amrex::Real> m_variables_mem;
    amrex::Gpu::DeviceVector<amrex::Real> m_rhs_mem;
    amrex::Gpu::DeviceVector<amrex::Real> m_y_old_mem;

    int m_naux_short1D;
    int m_naux_real1D;
    int m_naux_real2D;

    amrex::Gpu::DeviceVector<short> m_aux_short1D_mem;
    amrex::Gpu::DeviceVector<amrex::Real> m_aux_real1D_mem;
    amrex::Gpu::DeviceVector<amrex::Real> m_aux_real2D_mem;
};


/**
 * \brief Class that manages the memory to be used during the chemistry update
 *
 * This class stores the pointers to the memory that will be used in the
 * chemistry ODE integration. When running on the device, this class correctly
 * sets the pointers to either the GPU global memory or to the GPU shared
 * memory. This class must be set up within a GPU kernel as it is unique for
 * each GPU thread. Indeed, each GPU thread accesses a different element in the
 * vector of the problem variables, and in this class this is achieved by using
 * GPU directives such as blockIdx, blockDim , and threadIdx.
 */
class DeviceChemMemoryHandler
{
  public:
    /**
     * \brief Constructs an empty DeviceChemMemoryHandler
     *
     * Data can be defined at a later time using the different define member
     * functions.
     */
    AMREX_GPU_HOST_DEVICE
    DeviceChemMemoryHandler ()
      : m_use_shared_mem(0)
      , m_rates_mem_ptr(nullptr)
      , m_variables_mem_ptr(nullptr)
      , m_rhs_mem_ptr(nullptr)
      , m_y_old_mem_ptr(nullptr)
      , m_aux_short1D_mem_ptr(nullptr)
      , m_aux_real1D_mem_ptr(nullptr)
      , m_aux_real2D_mem_ptr(nullptr)
    {}

    /**
     * \brief Defines the sizes of the various memory areas which the
     * DeviceChemMemoryHandler class will point to.
     *
     * This member function assigns the memory sizes class members, such as the
     * number of chemical reactions, number of variables, and the size of the
     * problem (aka the number of fluid cells). It also initializes the indexes
     * class with the solids and fluid phase numbers of species. These indexes
     * are necproperly access variables such as the ODE solution and right-hand
     * side.
     *
     * \param size is the number of Box cells where chemistry is solved
     * \param nreactions is the maximum number between lagrangian and eulerian
     * chemical reactions
     * \param nspecies_g is the number of fluid species
     * \param nspecies_s is the number of solids species
     * \param solve_fluid is the flag for solving for the fluid phase  in
     * chemical reactions
     * \param solve_solids is the flag for solving for the solids phase in
     * chemical reactions
     * \param lagrangian_reactions is the flag for solving for lagrangian
     * chemical reactions
     * \param use_shared_mem is the flag for using GPU shared memory
     * \param n_aux_short1D is the number of auxiliary variables of type short
     * \param n_aux_real1D is the number of auxiliary variables of type Real
     * \param n_aux_real2D is the number of bidimensional auxiliary variables of
     * type Real
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup ([[maybe_unused]] const int size,
                const int nreactions,
                const int nspecies_g,
                const int nspecies_s,
                const int solve_fluid,
                const int solve_solids,
                const int lagrangian_reactions,
                const int use_shared_mem,
                const int n_aux_short1D,
                const int n_aux_real1D,
                const int n_aux_real2D)
    {
#ifdef AMREX_USE_GPU
      // Define the size of the problem and the number of chemical reactions
      m_size = size;
#endif
      m_nreactions = nreactions;

      // Setup the fluid and solids variables indexes
      m_indexes.setup_fluid(nspecies_g, solve_fluid);
      m_indexes.setup_solids(nspecies_s, solve_solids, lagrangian_reactions);

      // Get the number of variables of the chemistry ODE problem
      m_nvariables = m_indexes.count();

      // Set the flag for the GPU shared memory
      m_use_shared_mem = use_shared_mem;

      // Define the numbers of auxiliary variables for each different type
      m_naux_short1D = n_aux_short1D;
      m_naux_real1D = n_aux_real1D;
      m_naux_real2D = n_aux_real2D;
    }

    /**
     * \brief Set up the memory pointers
     *
     * This member function is GPU thread specific, and assigns the proper
     * pointers to the various memory areas that are specific to this thread. In
     * case GPU shared memory can be used, the pointers point to the shared
     * memory, otherwise they point to the GPU global memory.
     *
     * \param aux_class is the auxiliary object needed for setting up the shared
     * memory when using the GPU shared memory
     * \param rates_mem_ptr is the pointer to the memory area to be used for
     * storing the reaction rates
     * \param variables_mem_ptr is the pointer to the memory area to be used for
     * the solution of the ODE problem
     * \param rhs_mem_ptr is the pointer to the memory area to be used for the
     * right-hand side of the ODE problem
     * \param aux_short1D_mem_ptr is the pointer to the memory area to be used
     * for the ODE integrator auxiliary variables of type short
     * \param aux_real1D_mem_ptr is the pointer to the memory area to be used
     * for the ODE integrator auxiliary variables of type Real
     * \param aux_real2D_mem_ptr is the pointer to the memory area to be used
     * for the ODE integrator auxiliary bidimensional variables of type Real
     */
    template <class U>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void setup (U& aux_class,
                amrex::Real* rates_mem_ptr,
                amrex::Real* variables_mem_ptr,
                amrex::Real* rhs_mem_ptr,
                amrex::Real* y_old_mem_ptr,
                short* aux_short1D_mem_ptr,
                amrex::Real* aux_real1D_mem_ptr,
                amrex::Real* aux_real2D_mem_ptr)
    {
      if (m_use_shared_mem) {
        // The shared memory areas is shared among all the threads within the
        // same block of threads. It is a vectorized memory space, and we want
        // each warp of threads to access contiguous addresses of the shared
        // memory area to avoid bank conflicts. Thus, each thread will store
        // the different variables of the problem in addresses that are not
        // contiguous, but strided, with a stride between each variable that
        // is equal to the dimension of the block of threads. Thus we set the
        // stride equal to the block dimension.
      } else {
        // The global memory areas are vectorized memory spaces, and we want
        // each thread of a given block of threads to access contiguous
        // addresses in the GPU global memory to ensure that the reading and
        // writing accesses to the GPU blobal memory are coalesced. To do so,
        // we store the different variables of the problem in addresses that
        // are not contiguous, but strided, with a stride between each
        // variable that is equal to the size of the problem (aka the number
        // of fluid cells).  Thus we set the stride equal to the size of the
        // problem.
      }

#ifdef AMREX_USE_GPU
      m_glob_idx = blockDim.x*blockIdx.x+threadIdx.x;
      const int thread_idx = threadIdx.x;
      m_stride = m_use_shared_mem ? blockDim.x : m_size;
      amrex::Real* real_sm = m_use_shared_mem ? aux_class.dataPtr() : nullptr;

      if (m_use_shared_mem) {

        // The first variables that we store in the shared memory area are the
        // reactions rates. Each thread's pointer points to an address that is
        // indexed by the thread's index within the block of threads
        int gap(0);
        m_rates_mem_ptr = &(real_sm[gap*m_stride + thread_idx]);

        // Then, we store the ODE problem solution in the shared memory, after
        // the reactions rates. We start from the first element after all the
        // elements needed for the reactions rates (m_nreactions*m_stride),
        // and then each thread's pointer is different from the others by the
        // thread's index within the block of threads.
        gap += m_nreactions;
        m_variables_mem_ptr = &(real_sm[gap*m_stride + thread_idx]);

        // Then, we store the ODE problem right-hand side in the shared
        // memory. We start from the first element after all the elements
        // needed for the reactions rates and the ODE solution ((m_nreactions
        // + m_nvariables) * m_stride), and then each thread's pointer is
        // different from the others by the thread's index within the block of
        // threads.
        gap += m_nvariables;
        m_rhs_mem_ptr = &(real_sm[gap*m_stride + thread_idx]);

        //
        gap += m_nvariables;
        m_y_old_mem_ptr = &(real_sm[gap*m_stride + thread_idx]);

        // Then, we store the ODE integrator auxiliary variables of type Real
        // in the shared memory. We start from the first element after all the
        // elements needed for the reactions rates and the ODE solution and
        // right-hand side ((m_nreactions + 2*m_nvariables) * m_stride), and
        // then each thread's pointer is different from the others by the
        // thread's index within the block of threads.
        gap += m_nvariables;
        m_aux_real1D_mem_ptr = (m_naux_real1D > 0)?
          &(real_sm[gap*m_stride + thread_idx]) : nullptr;

        // Then, we store the ODE integrator auxiliary bidimensional variables
        // of type Real in the shared memory. We start from the first element
        // after all the elements needed for the reactions rates and the ODE
        // solution and right-hand side ((m_nreactions + 2*m_nvariables
        // + m_naux_real1D) * m_stride), and then each thread's pointer is
        // different from the others by the thread's index within the block of
        // threads.
        gap += m_naux_real1D;
        m_aux_real2D_mem_ptr = (m_naux_real2D > 0)?
          &(real_sm[gap*m_stride + thread_idx]) : nullptr;

        // Finally, we store the ODE integrator auxiliary variables of type
        // short in the shared memory. We start from the first element after all
        // the elements needed for the reactions rates and the ODE solution and
        // right-hand side ((m_nreactions + 2*m_nvariables + m_naux_real1D
        // + m_naux_real2D) * m_stride), and then each thread's pointer is
        // different from the others by the thread's index within the block of
        // threads.
        gap += m_naux_real2D;
        short* short_sm = (short*)(&real_sm[gap*m_stride]);
        m_aux_short1D_mem_ptr = (m_naux_short1D > 0)? &(short_sm[thread_idx]) : nullptr;

      } else { // GPU shared memory can't be used, so we use the global memory

        // Only set up pointers if the thread index is within the limits of the
        // allocated memory.
        if (m_glob_idx < m_size) {

          // This GPU thread's pointer to the reactions rates points to the
          // element with index m_glob_idx of the corresponding GPU global
          // memory space
          m_rates_mem_ptr = &(rates_mem_ptr[m_glob_idx]);

          // This GPU thread's pointer to the ODE solution points to the element
          // with index m_glob_idx of the corresponding GPU global memory space
          m_variables_mem_ptr = &(variables_mem_ptr[m_glob_idx]);

          // This GPU thread's pointer to the ODE right-hand side points to the
          // element with index m_glob_idx of the corresponding GPU global
          // memory space
          m_rhs_mem_ptr = &(rhs_mem_ptr[m_glob_idx]);

          //
          m_y_old_mem_ptr = &(y_old_mem_ptr[m_glob_idx]);

          // This GPU thread's pointer to the ODE integrator auxiliary variables
          // of type short points to the element with index m_glob_idx of the
          // corresponding GPU global memory space
          m_aux_short1D_mem_ptr = (m_naux_short1D > 0)?
            &(aux_short1D_mem_ptr[m_glob_idx]) : nullptr;

          // This GPU thread's pointer to the ODE integrator auxiliary variables
          // of type Real points to the element with index m_glob_idx of the
          // corresponding GPU global memory space
          m_aux_real1D_mem_ptr = (m_naux_real1D > 0)?
            &(aux_real1D_mem_ptr[m_glob_idx]) : nullptr;

          // This GPU thread's pointer to the ODE integrator auxiliary
          // bidimensional variables of type Real points to the element with
          // index m_glob_idx of the corresponding GPU global memory space
          m_aux_real2D_mem_ptr = (m_naux_real2D > 0)?
            &(aux_real2D_mem_ptr[m_glob_idx]) : nullptr;
        }
      }

#else

      // The global index is simply the index of the ParallelFor loop
      m_glob_idx = aux_class;
      // For non-GPU runs, the stride between different variables is just equal
      // to 1
      m_stride = 1;

      // For non-GPU runs, the memory pointers are simply equal to the pointers
      // passed as arguments.
      m_rates_mem_ptr = rates_mem_ptr;
      m_variables_mem_ptr = variables_mem_ptr;
      m_rhs_mem_ptr = rhs_mem_ptr;
      m_y_old_mem_ptr = y_old_mem_ptr;

      // For non-GPU runs, the memory pointers to the ODE integrator auxiliary
      // variables are simply equal to the pointers passed as arguments.
      m_aux_short1D_mem_ptr = (m_naux_short1D > 0)? aux_short1D_mem_ptr : nullptr;
      m_aux_real1D_mem_ptr = (m_naux_real1D > 0)? aux_real1D_mem_ptr : nullptr;
      m_aux_real2D_mem_ptr = (m_naux_real2D > 0)? aux_real2D_mem_ptr : nullptr;
#endif
    }


    /**
     * \brief Get the memory address to the reactions rates variable with given
     * index q
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* rates_address (const int q = 0) const
    {
      AMREX_ASSERT(m_nreactions != 0);
      AMREX_ASSERT(q >= 0 && q < m_nreactions);
      return &(m_rates_mem_ptr[q*m_stride]);
    }

    /**
     * \brief Get the memory address to the ODE solution variable with given
     * index n
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* variables_address (const uint64_t n = 0) const
    {
      AMREX_ASSERT(m_indexes.count() != 0);

      if (n == m_indexes.count())
        return nullptr;

      AMREX_ASSERT(n < m_indexes.count());
      return &(m_variables_mem_ptr[n*m_stride]);
    }

    /**
     * \brief Get the memory address to the ODE right-hand side variable with
     * given index n
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* rhs_address (const uint64_t n = 0) const
    {
      AMREX_ASSERT(m_indexes.count() != 0);

      if (n == m_indexes.count())
        return nullptr;

      AMREX_ASSERT(n < m_indexes.count());
      return &(m_rhs_mem_ptr[n*m_stride]);
    }

    /**
     * \brief
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* y_old_address (const uint64_t n = 0) const
    {
      AMREX_ASSERT(m_indexes.count() != 0);

      if (n == m_indexes.count())
        return nullptr;

      AMREX_ASSERT(n < m_indexes.count());
      return &(m_y_old_mem_ptr[n*m_stride]);
    }

    /**
     * \brief Get the memory address to the ODE solution fluid-phase variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* fluid_vars_address () const
    {
      if (m_indexes.fluid_indexes().count() == 0)
        return nullptr;

      return variables_address(m_indexes.fluid_first_idx());
    }

    /**
     * \brief Get the memory address to the ODE solution solid-phase variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* solids_vars_address () const
    {
      if (m_indexes.solids_indexes().count() == 0)
        return nullptr;

      return variables_address(m_indexes.solids_first_idx());
    }

    /**
     * \brief Get the memory address to the ODE right-hand side fluid-phase
     * variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* fluid_rhs_address () const
    {
      if (m_indexes.fluid_indexes().count() == 0)
        return nullptr;

      return rhs_address(m_indexes.fluid_first_idx());
    }

    /**
     * \brief Get the memory address to the ODE right-hand side solid-phase
     * variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* solids_rhs_address () const
    {
      if (m_indexes.solids_indexes().count() == 0)
        return nullptr;

      return rhs_address(m_indexes.solids_first_idx());
    }

    /**
     * \brief Get the memory stride
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int const& stride () const { return m_stride; }

    /**
     * \brief Get the number of variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int const& nvariables () const { return m_nvariables; }

    /**
     * \brief Get the global index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int const& glob_idx () const { return m_glob_idx; }

    /**
     * \brief Get a reference to the ChemMemoryIndexes class
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const ChemMemoryIndexes& indexes () const { return m_indexes; }

    /**
     * \brief Get the memory address to the auxiliary variables of type short
     * for the chosen integrator
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    short* aux_short1D_mem () const
    { return m_aux_short1D_mem_ptr; }

    /**
     * \brief Get the memory address to the auxiliary variables of type Real
     * for the chosen integrator
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* aux_real1D_mem () const
    { return m_aux_real1D_mem_ptr; }

    /**
     * \brief Get the memory address to the auxiliary bidimensional variables of
     * type Real for the chosen integrator
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real* aux_real2D_mem () const
    { return m_aux_real2D_mem_ptr; }

  private:
    int m_glob_idx;
    int m_stride;
    int m_use_shared_mem;

#ifdef AMREX_USE_GPU
    int m_size;
#endif
    int m_nreactions;
    int m_nvariables;

    ChemMemoryIndexes m_indexes;

    amrex::Real* m_rates_mem_ptr;
    amrex::Real* m_variables_mem_ptr;
    amrex::Real* m_rhs_mem_ptr;
    amrex::Real* m_y_old_mem_ptr;

    int m_naux_short1D;
    int m_naux_real1D;
    int m_naux_real2D;

    short* m_aux_short1D_mem_ptr;
    amrex::Real* m_aux_real1D_mem_ptr;
    amrex::Real* m_aux_real2D_mem_ptr;
};


#endif
