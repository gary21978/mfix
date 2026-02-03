#ifndef MFIX_CHEM_MEMORY_STRUCTS_H
#define MFIX_CHEM_MEMORY_STRUCTS_H

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Gpu.H>


/**
 * \brief This class stores the indexes for accessing the memory variables used
 * in the ODE integration
 *
 * The class is a wrapper for other two classes, one for the fluid phase and
 * another for the solids phase, which contain the indexes to access the
 * variables for each single phase. This class combines the two classes indexes
 * to provide the global indexing for the whole ODE system to be solved.
 */
class ChemMemoryIndexes
{
  public:
    /**
     * \brief Constructs an empty ChemMemoryIndexes
     */
    AMREX_GPU_HOST_DEVICE
    ChemMemoryIndexes ()
      : m_fluid_is_initialized(0)
      , m_solids_is_initialized(0)
      , m_fluid_indexes()
      , m_solids_indexes()
    {}

    /**
     * \brief This class stores the indexes for accessing the memory variables
     * for the fluid phase.
     *
     * The indexes of the variables for the ODE solution and right-hand side for
     * the fluid phase are in the following order:
     * - from 0 to N_g-1: the fluid species densities \f$ \rho_g X_{gk} \f$;
     * - N_g: the fluid enthalpy density \f$ \rho_g h_g \f$.
     * where \f$ \rho_g \f$, \f$ h_g \f$, and \f$ X_{gk} \f$ are respectively
     * the fluid density, specific enthalpy, and mass fractions. \f$ N_g \f$ is
     * the number of fluid species.
     */
    class ChemFluidIndexes
    {
      public:
        static constexpr uint64_t m_density = 0;
        static constexpr uint64_t m_temperature = 1;
        static constexpr uint64_t m_specific_enthalpy = 2;
        static constexpr uint64_t m_mass_fractions = 3;

        /**
         * \brief Constructs an empty ChemFluidIndexes
         *
         * Data can be defined at a later time using the define member function.
         */
        AMREX_GPU_HOST_DEVICE
        ChemFluidIndexes ()
          : m_nspecies(0)
          , m_count(0)
        {}

        /**
         * \brief Defines the ChemFluidIndexes class
         *
         * This member function initializes the fluid-phase indexes in the
         * following order:
         * - from 0 to N_g-1: the fluid species densities \f$ \rho_g X_{gk} \f$;
         * - N_g: the fluid enthalpy density \f$ \rho_g h_g \f$.
         * where \f$ \rho_g \f$, \f$ h_g \f$, and \f$ X_{gk} \f$ are
         * respectively the fluid density, specific enthalpy, and mass
         * fractions. \f$ N_g \f$ is the number of fluid species.
         *
         * \param nspecies_g is the number of fluid species
         * \param solve_fluid is the flag for solving for the fluid phase in
         * chemical reactions
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void define (const uint64_t nspecies_g,
                     const uint64_t solve_fluid)
        {
          if (solve_fluid) {
            m_nspecies = nspecies_g;
            m_count = nspecies_g + m_mass_fractions;
          }
        }

        /**
         * \brief Get the value of the first fluid index
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t first_idx () { return m_density; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t density () { return m_density; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t temperature () { return m_temperature; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t specific_enthalpy () { return m_specific_enthalpy; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t mass_fractions () { return m_mass_fractions; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        uint64_t const& nspecies () const { return m_nspecies; }

        /**
         * \brief Get the number of variables
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        uint64_t const& count () const { return m_count; }

      private:
        uint64_t m_nspecies;
        uint64_t m_count;
    };

    /**
     * \brief This class stores the indexes for accessing the memory variables
     * for the solids phase.
     *
     * The indexes of the variables for the ODE solution and right-hand side for
     * the solids phase are in the following order:
     * - from 0 to N_s-1: the solids species densities \f$ \m_p X_{pk} \f$;
     * - N_s: the solids enthalpy \f$ m_p h_p \f$.
     * where \f$ \m_p \f$, \f$ h_p \f$, and \f$ X_{pk} \f$ are respectively
     * the particle mass, specific enthalpy, and mass fractions. \f$ N_s \f$ is
     * the number of solids species.
     */
    class ChemSolidsIndexes
    {
      public:
        static constexpr uint64_t m_density = 0;
        static constexpr uint64_t m_radius = 1;
        static constexpr uint64_t m_temperature = 2;
        static constexpr uint64_t m_specific_enthalpy = 3;
        static constexpr uint64_t m_mass_fractions = 4;

        /**
         * \brief Constructs an empty ChemSolidsIndexes
         *
         * Data can be defined at a later time using the define member function.
         */
        AMREX_GPU_HOST_DEVICE
        ChemSolidsIndexes ()
          : m_nspecies(0)
          , m_count(0)
        {}

        /**
         * \brief Defines the ChemSolidsIndexes class
         *
         * This member function initializes the solids-phase indexes  in the
         * following order:
         * - from 0 to N_s-1: the solids species densities \f$ \m_p X_{pk} \f$;
         * - N_s: the solids enthalpy \f$ m_p h_p \f$.
         * where \f$ \m_p \f$, \f$ h_p \f$, and \f$ X_{pk} \f$ are respectively
         * the particle mass, specific enthalpy, and mass fractions. \f$ N_s \f$
         * is the number of solids species.
         *
         * \param nspecies_s is the number of solids species
         * \param solve_solids is the flag for solving for the solids phase in
         * chemical reactions
         * \param lagrangian_reactions is the flag for solving for the
         * lagrangian reactions
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void define (const uint64_t nspecies_s,
                     const uint64_t solve_solids,
                     const uint64_t lagrangian_reactions)
        {
          if (solve_solids && lagrangian_reactions) {
            m_nspecies = nspecies_s;
            m_count = nspecies_s + m_mass_fractions;
          }
        }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t first_idx () { return m_density; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t density () { return m_density; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t radius () { return m_radius; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t temperature () { return m_temperature; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t specific_enthalpy () { return m_specific_enthalpy; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr uint64_t mass_fractions () { return m_mass_fractions; }

        /**
         * \brief
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        uint64_t const& nspecies () const { return m_nspecies; }

        /**
         * \brief Get the number of variables
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        uint64_t const& count () const { return m_count; }

      private:
        uint64_t m_nspecies;
        uint64_t m_count;
    };

    /**
     * \brief Call the define member function of the fluid indexes class and
     * sets the fluid initialization flag to true
     *
     * \param nspecies_g is the number of fluid species
     * \param solve_fluid is the flag for solving for the fluid phase in
     * chemical reactions
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup_fluid (const uint64_t nspecies_g,
                      const int solve_fluid)
    {
      m_fluid_is_initialized = 1;
      m_fluid_indexes.define(nspecies_g, solve_fluid);
    }

    /**
     * \brief Call the define member function of the solids indexes class and
     * sets the initialization flag to true
     *
     * \param nspecies_s is the number of solids species
     * \param solve_solids is the flag for solving for the solids phase in
     * chemical reactions
     * \param lagrangian_reactions is the flag for solving for lagrangian
     * reactions
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup_solids (const uint64_t nspecies_s,
                       const int solve_solids,
                       const int lagrangian_reactions)
    {
      m_solids_is_initialized = 1;
      m_solids_indexes.define(nspecies_s, solve_solids, lagrangian_reactions);
    }

    /**
     * \brief Get the first index for the fluid variables in the whole set of
     * variables for the ODE solution
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr uint64_t fluid_first_idx ()
    {
      return ChemFluidIndexes::first_idx();
    }

    /**
     * \brief Get the first index for the solids variables in the whole set of
     * variables for the ODE solution
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t solids_first_idx () const
    {
      AMREX_ASSERT(m_solids_is_initialized);
      return m_fluid_indexes.count() + ChemSolidsIndexes::first_idx();
    }

    /**
     * \brief Get the first index for the fluid species densities variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr uint64_t fluid_species_begin ()
    {
      return ChemFluidIndexes::mass_fractions();
    }

    /**
     * \brief Get the last index for the fluid species densities variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t fluid_species_end () const
    {
      AMREX_ASSERT(m_fluid_is_initialized);
      return ChemFluidIndexes::mass_fractions() + m_fluid_indexes.nspecies();
    }

    /**
     * \brief
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t solids_density () const
    {
      AMREX_ASSERT(m_solids_is_initialized);
      return m_fluid_indexes.count() + ChemSolidsIndexes::density();
    }

    /**
     * \brief
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t solids_radius () const
    {
      AMREX_ASSERT(m_solids_is_initialized);
      return m_fluid_indexes.count() + ChemSolidsIndexes::radius();
    }

    /**
     * \brief Get the first index for the solids species masses variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t solids_species_begin () const
    {
      AMREX_ASSERT(m_solids_is_initialized);
      return m_fluid_indexes.count() + ChemSolidsIndexes::mass_fractions();
    }

    /**
     * \brief Get the last index for the solids species masses variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t solids_species_end () const
    {
      AMREX_ASSERT(m_solids_is_initialized);
      return m_fluid_indexes.count() + ChemSolidsIndexes::mass_fractions() + m_solids_indexes.nspecies();
    }

    /**
     * \brief Get the total number of ODE variables
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    uint64_t count () const
    {
      AMREX_ASSERT(m_fluid_is_initialized || m_solids_is_initialized);
      return m_fluid_indexes.count() + m_solids_indexes.count();
    }

    /**
     * \brief Get a reference to the fluid indexes class
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const ChemFluidIndexes& fluid_indexes () const
    {
      AMREX_ASSERT(m_fluid_is_initialized);
      return m_fluid_indexes;
    }

    /**
     * \brief Get a reference to the solids indexes class
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const ChemSolidsIndexes& solids_indexes () const
    {
      AMREX_ASSERT(m_solids_is_initialized);
      return m_solids_indexes;
    }

  private:
    int m_fluid_is_initialized;
    int m_solids_is_initialized;
    ChemFluidIndexes m_fluid_indexes;
    ChemSolidsIndexes m_solids_indexes;
};


/**
 * \brief This class allows to properly address the variables of the ODE problem
 * and ODE integrator on the host or device memory during the ODE integration
 * and for the chemistry update
 *
 * The class is a wrapper for a pointer to a given memory area, and it allows to
 * access the elements of the pointer by a given stride. This is particularly
 * important when accessing GPU memory, as it allows to avoid GPU uncoalesced
 * accesses to its global memory and/or bank conflicts when accessing the GPU
 * shared memory.
 */
template <class T>
class ChemData
{
  public:
    /**
     * \brief Class constructor
     *
     * This class constructor gets the pointer to the memory area to be accessed
     * and the stride that will be used for accessing its elements.
     *
     * \param data_ptr is the pointer to the memory area to be accessed
     * \param size is the number of Box cells where chemistry is solved
     * \param stride is the stride to use when accessing the pointed memory
     * elements
     */
    AMREX_GPU_HOST_DEVICE
    ChemData (T* data_ptr,
              const uint64_t size,
              const uint64_t stride)
      : m_data_ptr(data_ptr)
      , m_size(size)
      , m_stride(stride)
    {}

    /**
     * \brief Get a copy of this class' pointer to the n-th variable in the
     * pointed memory space, which is at address n*stride
     *
     * \param n is the index of the variable
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* ptr (const uint64_t& n)
    {
      AMREX_ASSERT(n < m_size);
      return &(m_data_ptr[n*m_stride]);
    }

  private:
    T* m_data_ptr;
    uint64_t m_size;
    uint64_t m_stride;
};


/**
 * \brief This class allows to properly address the elements of the vectors on
 * the host or device during the ODE integration and for the chemistry update
 *
 * The class is a wrapper for a pointer to a given memory area, and it allows to
 * access the elements of the pointer by a given stride. This is particularly
 * important when accessing GPU memory, as it allows to avoid GPU uncoalesced
 * accesses to its global memory and/or bank conflicts when accessing the GPU
 * shared memory.
 */
template <class T>
class ChemVector
{
  public:
    /**
     * \brief Constructs an empty ChemVector
     *
     * Data can be defined at a later time using the define member function.
     */
    AMREX_GPU_HOST_DEVICE
    ChemVector () = default;

    /**
     * \brief Class constructor
     *
     * This class constructor gets the pointer to the memory area to be accessed
     * and the stride that will be used for accessing its elements.
     *
     * \param data_ptr is the pointer to the memory area to be accessed
     * \param size is the number of Box cells where chemistry is solved
     * \param stride is the stride to use when accessing the pointed memory
     * elements
     */
    AMREX_GPU_HOST_DEVICE
    ChemVector (T* data_ptr,
                const uint64_t size,
                const uint64_t stride)
      : m_data_ptr(data_ptr)
      , m_size(size)
      , m_stride(stride)
    {}

    /**
     * \brief Defines the ChemVector class
     *
     * This class constructor gets the pointer to the memory area to be accessed
     * and the stride that will be used for accessing its elements.
     *
     * \param data_ptr is the pointer to the memory area to be accessed
     * \param size is the number of Box cells where chemistry is solved
     * \param stride is the stride to use when accessing the pointed memory
     * elements
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void define (T* data_ptr,
                 const uint64_t size,
                 const uint64_t stride)
    {
      m_data_ptr = data_ptr;
      m_size = size;
      m_stride = stride;
    }

    /**
     * \brief Fill the vector with a given value
     *
     * \param value is the value to be used to fill up the vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void fill (const T& value)
    {
      for (uint64_t n(1); n <= m_size; ++n)
        m_data_ptr[(n-1)*m_stride] = value;
    }

    /**
     * \brief Get a constant reference to the (one-based) n-th component, which
     * is at address (n-1)*stride
     *
     * \param n is the one-based index of the component for this vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T const& operator() (const uint64_t& n) const
    {
      AMREX_ASSERT(n <= m_size);
      return m_data_ptr[(n-1)*m_stride];
    }

    /**
     * \brief Get a reference to the (one-based) n-th component, which is at
     * address (n-1)*stride
     *
     * \param n is the one-based index of the component for this vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator() (const uint64_t& n)
    {
      AMREX_ASSERT(n <= m_size);
      return m_data_ptr[(n-1)*m_stride];
    }

    /**
     * \brief Get a constant reference to the (zero-based) n-th component, which
     * is at address n*stride
     *
     * \param n is the zero-based index of the component for this vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T const& operator[] (const uint64_t& n) const
    {
      AMREX_ASSERT(n < m_size);
      return m_data_ptr[n*m_stride];
    }

    /**
     * \brief Get a reference to the (zero-based) n-th component, which is at
     * address n*stride
     *
     * \param n is the zero-based index of the component for this vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator[] (const uint64_t& n)
    {
      AMREX_ASSERT(n < m_size);
      return m_data_ptr[n*m_stride];
    }

    /**
     * \brief Get a reference to this class' pointer to the memory space
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T*& dataPtr ()
    {
      return m_data_ptr;
    }

    /**
     * \brief Get a copy of this class' pointer to the memory space
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* const& dataPtr () const
    {
      return m_data_ptr;
    }

    /**
     * \brief Get a copy of this class' pointer to the n-th component in the
     * pointed memory space, which is at address n*stride
     *
     * \param n is the zero-based index of the component for this vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* ptr (const uint64_t& n)
    {
      AMREX_ASSERT(n < m_size);
      return &(m_data_ptr[n*m_stride]);
    }

  private:
    T* m_data_ptr = nullptr;
    uint64_t m_size = 0;
    uint64_t m_stride = 0;
};


/**
 * \brief This class allows to properly address the elements of the matrices on
 * the host or device during the ODE integration and for the chemistry update
 *
 * The class is a wrapper for a pointer to a given memory area, and it allows to
 * access the elements of the pointer by a given stride. This is particularly
 * important when accessing GPU memory, as it allows to avoid GPU uncoalesced
 * accesses to its global memory and/or bank conflicts when accessing the GPU
 * shared memory.
 */
template <class T>
class ChemMatrix
{
  public:
    /**
     * \brief Constructs an empty ChemVector
     *
     * Data can be defined at a later time using the define member function.
     */
    AMREX_GPU_HOST_DEVICE
    ChemMatrix () = default;

    /**
     * \brief Class constructor
     *
     * This class constructor gets the pointer to the memory area to be accessed
     * and the stride that will be used for accessing its elements.
     *
     * \param data_ptr is the pointer to the memory area to be accessed
     * \param nrows is the number of rows of the bidimensional variable
     * \param ncols is the number of columns of the bidimensional variable
     * \param stride is the stride to use when accessing the pointed memory
     * elements
     */
    AMREX_GPU_HOST_DEVICE
    ChemMatrix (T* data_ptr,
                const uint64_t nrows,
                const uint64_t ncols,
                const uint64_t stride)
      : m_data_ptr(data_ptr)
      , m_nrows(nrows)
      , m_ncols(ncols)
      , m_stride(stride)
    {}

    /**
     * \brief Defines the ChemVector class
     *
     * This class constructor gets the pointer to the memory area to be accessed
     * and the stride that will be used for accessing its elements.
     *
     * \param data_ptr is the pointer to the memory area to be accessed
     * \param nrows is the number of rows of the bidimensional variable
     * \param ncols is the number of columns of the bidimensional variable
     * \param stride is the stride to use when accessing the pointed memory
     * elements
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void define (T* data_ptr,
                 const uint64_t nrows,
                 const uint64_t ncols,
                 const uint64_t stride)
    {
      m_data_ptr = data_ptr;
      m_nrows = nrows;
      m_ncols = ncols;
      m_stride = stride;

      this->fill(static_cast<T>(0));
    }

    /**
     * \brief Fill the vector with a given value
     *
     * \param value is the value to be used to fill up the matrix
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void fill (const T& value)
    {
      for (uint64_t m(1); m <= m_nrows; ++m)
        for (uint64_t n(1); n <= m_ncols; ++n)
          m_data_ptr[((m-1)*m_ncols + (n-1))*m_stride] = value;
    }


    /**
     * \brief Get a constant reference to the (one-based) m-th and n-th
     * component, which is at address ((m-1)*m_ncols + (n-1))*stride
     *
     * \param m is the one-based index of the row component for this matrix
     * \param n is the one-based index of the column component for this matrix
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T const& operator() (const uint64_t& m,
                         const uint64_t& n) const
    {
      AMREX_ASSERT(m <= m_nrows);
      AMREX_ASSERT(n <= m_ncols);

      return m_data_ptr[((m-1)*m_ncols + (n-1))*m_stride];
    }

    /**
     * \brief Get a reference to the (one-based) m-th and n-th component, which
     * is at address ((m-1)*m_ncols + (n-1))*stride
     *
     * \param m is the one-based index of the row component for this matrix
     * \param n is the one-based index of the column component for this matrix
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator() (const uint64_t& m,
                   const uint64_t& n)
    {
      AMREX_ASSERT(m <= m_nrows);
      AMREX_ASSERT(n <= m_ncols);

      return m_data_ptr[((m-1)*m_ncols + (n-1))*m_stride];
    }

  private:
    T* m_data_ptr = nullptr;
    uint64_t m_nrows = 0;
    uint64_t m_ncols = 0;
    uint64_t m_stride = 0;
};


#endif
