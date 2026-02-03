#ifndef MFIX_MATERIAL_INTERFACE_H
#define MFIX_MATERIAL_INTERFACE_H

#include <AMReX_REAL.H>

#include <mfix_actual_rhs.H>
#include <mfix_reactions.H>
#include <mfix_chem_data_structs.H>
#include <vode_integrator.H>
#include <mfix_chem_memory_handler.H>
#include <mfix_usr_reactions_rates_K.H>


/**
 * \brief This class manages the chemistry ODE integration for each fluid cell
 * inside a ParallelFor loop.
 *
 * This class takes care of putting together all the pieces that are necessary
 * for the ODE system integration, and then performans the actual ODE
 * integration in the fluid cell in which it has been defined. For each fluid
 * cell we need to define one MaterialInterface class. This class takes the
 * fluid-phase data, the solids data, the reactions data, the integrator, and
 * the memory pointers for the ODE solution, ODE right-hand side, and reactions
 * rates. Using all the gathered information, it performs the actual time
 * integration of the chemical problem.
 */
class MaterialInterface {

  public:
    /**
     * \brief Constructs an empty MaterialInterface class
     *
     * Data can be defined at a later time using the define member function.
     */
    AMREX_GPU_HOST_DEVICE
    MaterialInterface ()
      : m_fluid_data(nullptr)
      , m_solids_data(nullptr)
      , m_eulerian_reactions_data(nullptr)
      , m_lagrangian_reactions_data(nullptr)
      , m_memory_handler()
    {}

    /**
     * \brief Define the memory handler of this class
     *
     * This member function defines the memory handler object by setting up its
     * memory pointers and the corresponding sizes for properly accessing the
     * allocated memory
     *
     * \param size is the number of Box cells where chemistry is solved
     * \param nreactions is the number of chemical reactions
     * \param nspecies_g is the number of species for the fluid phase
     * \param nspecies_s is the number of species for the solids phase
     * \param use_shared_mem is the flag to switch on/off the device shared
     * memory
     * \param mem_aux is the auxiliary template object: admissible types are
     * device shared memory for GPU compilations or the ParallelFor loop index
     * for non-GPU compilations.
     * \param rates_mem_ptr is the pointer to the memory area dedicated to the
     * reactions rates
     * \param variables_mem_ptr is the pointer to the memory area dedicated to
     * the ODE solution
     * \param rhs_mem_ptr is the pointer to the memory are dedicated to the ODE
     * right-hand side
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup (const int size,
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
      m_memory_handler.setup(size, nreactions, nspecies_g, nspecies_s, solve_fluid, solve_solids,
          lagrangian_reactions, use_shared_mem, n_aux_short1D, n_aux_real1D, n_aux_real2D);
    }

    /**
     * \brief Get a reference to the memory handler
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    DeviceChemMemoryHandler& memory_handler () { return m_memory_handler; }

    /**
     * \brief Calls the setup function for the memory handler object in this
     * class
     */
    template <class U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup_memory (U& mem_aux,
                       amrex::Real* rates_mem_ptr,
                       amrex::Real* variables_mem_ptr,
                       amrex::Real* rhs_mem_ptr,
                       amrex::Real* y_old_mem_ptr,
                       short* aux_short1D_mem_ptr,
                       amrex::Real* aux_real1D_mem_ptr,
                       amrex::Real* aux_real2D_mem_ptr)
    {
      m_memory_handler.setup(mem_aux, rates_mem_ptr, variables_mem_ptr, rhs_mem_ptr,
          y_old_mem_ptr, aux_short1D_mem_ptr, aux_real1D_mem_ptr, aux_real2D_mem_ptr);
    }

    /**
     * \brief Define the eulerian reactions
     *
     * \param reactions_data is the pointer to the eulerian reactions class
     * \param solve_eulerian_reactions is the flag for solving for eulerian
     * chemical reactions
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup_eulerian_reactions (MFIXReactionsData* reactions_data,
                                   const int solve_eulerian_reactions)
    {
      m_solve_eulerian_reactions = solve_eulerian_reactions;
      m_eulerian_reactions_data = solve_eulerian_reactions? reactions_data : nullptr;
    }

    /**
     * \brief Define the lagrangian reactions
     *
     * \param reactions_data is the pointer to the lagrangian reactions class
     * \param solve_lagrangian_reactions is the flag for solving for lagrangian
     * chemical reactions
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup_lagrangian_reactions (MFIXReactionsData* reactions_data,
                                     const int solve_lagrangian_reactions)
    {
      m_solve_lagrangian_reactions = solve_lagrangian_reactions;
      m_lagrangian_reactions_data = solve_lagrangian_reactions? reactions_data : nullptr;
    }

    /**
     * \brief Define the fluid phase
     *
     * \param fluid_data is the pointer to the fluid data class
     * \param solve_fluid is the flag for solving for the fluid phase in
     * chemical reactions
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup_fluid (MFIXFluidData* fluid_data,
                      const int solve_fluid)
    {
      m_solve_fluid = solve_fluid;
      m_fluid_data = solve_fluid? fluid_data : nullptr;
    }

    /**
     * \brief Define the solids phase
     *
     * \param solids_data is the pointer to the solids data class
     * \param solve_solids is the flag for solving for the solids phase in
     * chemical reactions
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setup_solids (MFIXSolidsData* solids_data,
                       const int solve_solids)
    {
      m_solve_solids = solve_solids;
      m_solids_data = solve_solids? solids_data : nullptr;
    }

    /**
     * \brief Compute the chemical reactions by solving the ODE system and
     * updating the ODE solution
     */
    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void compute_reactions (T* integrator,
                            const amrex::Real dt) const
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("mfix::compute_reactions()");
#endif

      // Get a reference to the memory indexes to properly access the ODE memory
      // variables
      const ChemMemoryIndexes& indexes = m_memory_handler.indexes();

      // Compute the eulerian reactions if the user defined any
      if (m_solve_eulerian_reactions && m_solve_fluid) {

        // Initialize the fluid ODE solution by copying the fluid quantities
        // into the ODE variables
        m_fluid_data->init_vars();

        if constexpr (!std::is_same<T, VODE>::value) {
          integrator->change_size(indexes.fluid_indexes().count());
        }

        // Define the ODE right-hand side by setting the reactions type to
        // eulerian and passing the information that is necessary to compute the
        // right-hand side
        integrator->actual_rhs.setup(ReactionType::Eulerian, m_solve_fluid,
            0, m_fluid_data, nullptr, m_eulerian_reactions_data, indexes);

        // Reset to zero the eulerian and lagrangian reactions rates
        m_eulerian_reactions_data->reset_rates();

        m_fluid_data->reset_rhs();

        integrator->set_status(IntegrationStatus::IERR_SUCCESS);

        // Perform the actual ODE integration
        integrator->actual_integrator(dt);
      }

      // Compute the lagrangian reactions if the user defined any
      if (m_solve_lagrangian_reactions) {

        if (!m_solve_eulerian_reactions && m_solve_fluid) {
          m_fluid_data->init_vars();
        }

        if constexpr (!std::is_same<T, VODE>::value) {
          integrator->change_size(indexes.count());
        }

        // Define the ODE right-hand side by setting the reactions type to
        // lagrangian and passing the information that is necessary to compute
        // the right-hand side
        integrator->actual_rhs.setup(ReactionType::Lagrangian, m_solve_fluid,
            m_solve_solids, m_fluid_data, m_solids_data, m_lagrangian_reactions_data,
            indexes);

        // Get this fluid cell indexes
        const amrex::IntVect ijk = m_fluid_data->cell_ijk();
        // Get the number of particles in the current fluid cell
        const amrex::Long cell_parts_nb = m_solids_data->particles_nb_in_cell(ijk);

        // TODO
        // consider computing the nb of particles in a cell by subtracting the
        // global indexes between the current ijk cell and the "previous" ijk-1
        // cell, paying attention to consider as 0 the global index in the "-1"
        // cell.
        // This would allow to save the memory space for the MultiFab containing
        // the nb of particles per cell

        // Loop over the particles in this fluid cell
        for (int n(cell_parts_nb); (n > 0) && m_solve_solids; n--) {

          // Update the ID index in the solids data class with current particle
          // global ID. This is important to correctly access the
          // NeighborParticleContainer data for the current particle
          m_solids_data->update_p_id(n);

          // Initialize the solids ODE solution by copying this particle's
          // quantities into the ODE variables
          m_solids_data->init_vars();

          // Reset the lagrangian reactions rates for the current particle
          m_lagrangian_reactions_data->reset_rates();

          m_fluid_data->reset_rhs();
          m_solids_data->reset_rhs();

          integrator->set_status(IntegrationStatus::IERR_SUCCESS);

          // Perform the actual ODE integration
          integrator->actual_integrator(dt);

          // Update this particle's quantities with the solids variables in the
          // ODE solution.
          m_solids_data->update_variables();
        }
      }

      // Update the fluid phase chamical transfer quantities
      if (m_solve_fluid) {
        m_fluid_data->update_chem_txfr(dt);
      }
    }

  private:
    MFIXFluidData* m_fluid_data;
    MFIXSolidsData* m_solids_data;
    MFIXReactionsData* m_eulerian_reactions_data;
    MFIXReactionsData* m_lagrangian_reactions_data;

    int m_solve_fluid;
    int m_solve_solids;
    int m_solve_eulerian_reactions;
    int m_solve_lagrangian_reactions;

    DeviceChemMemoryHandler m_memory_handler;
};


#endif
