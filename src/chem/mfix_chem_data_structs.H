#ifndef MFIX_CHEM_DATA_STRUCTS_H
#define MFIX_CHEM_DATA_STRUCTS_H

#include <limits>

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
#include <integrator_parameters.H>
#include <mfix_chem_memory_handler.H>
#include <mfix_chem_memory_structs.H>


class MFIXActualRHS;
class MaterialInterface;


/**
 * \brief Class that stores all the fluid phase data needed in chemistry
 *
 * This class contains all the fluid phase data that can be useful to define the
 * reactions rates functions, to compute the chemistry ODE right-hand side, and
 * to integrate the chemistry ODE. In addition, this class stores the pointers
 * to the memory areas that will be used to save the ODE unknowns and right-hand
 * side.
 */
class MFIXFluidData
{
  // Define classes MFIXActualRHS and MaterialInterface friends so they can
  // access the private members of this class
  friend class MFIXActualRHS;
  friend class MaterialInterface;

  public:
    using ChemFluidIndexes = ChemMemoryIndexes::ChemFluidIndexes;

    /**
     * \brief Constructs the MFIXFluidData class
     *
     * Constructor of the class
     *
     * \param memory_handler is the class that contains all the information
     * about the memory areas that this MFIXFluidData class can use to store ODE
     * unknowns and right-hand side data
     * \param fluid_parms is the MFIXFluidParms class reference
     * \param ijk are the indexes for current fluid cell
     * \param ep_g is the fluid volume fraction
     * \param ro_g is the fluid density
     * \param vel_g is the fluid velocity
     * \param T_g is the fluid temperature
     * \param h_g is the fluid specific enthalpy
     * \param X_gk is the fluid species mass fractions
     * \param p_g is the fluid thermodynamic pressure
     * \param chem_txfr_g is the fluid transfer quantities due to chemical
     * reactions.  This is where the computed transfer quantities of the fluid
     * phase will be stored
     * \param chem_txfr_idxs is the class containing the indexes to properly
     * access the "chem_txfr_g" Array4 components
     * \param fluid_is_mixture is the flag to discriminate when the fluid is
     * a mixture of species or not
     */
    AMREX_GPU_HOST_DEVICE
    MFIXFluidData (DeviceChemMemoryHandler& memory_handler,
                   const MFIXFluidParms& fluid_parms,
                   const ThermoPropertyData& a_fluid_props,
                   const amrex::IntVect& ijk,
                   amrex::Array4<const amrex::Real> const& ep_g,
                   amrex::Array4<const amrex::Real> const& ro_g,
                   amrex::Array4<const amrex::Real> const& vel_g,
                   amrex::Array4<const amrex::Real> const& T_g,
                   amrex::Array4<const amrex::Real> const& h_g,
                   amrex::Array4<const amrex::Real> const& X_gk,
                   const amrex::Real& /*p_g*/,
                   const amrex::Real& cell_volume,
                   const amrex::Real& vfrac,
                   amrex::Array4<amrex::Real> const& chem_txfr_g,
                   const InterphaseChemTxfrIndexes& chem_txfr_idxs,
                   const int report_mass_balance,
                   amrex::Real* fluid_prod_ptr)
      : params(fluid_parms)
      , m_props(a_fluid_props)
      , m_y(memory_handler.fluid_vars_address(), memory_handler.nvariables(), memory_handler.stride())
      , m_persistent_y(memory_handler.fluid_vars_address(), memory_handler.nvariables(), memory_handler.stride())
      , m_rhs(memory_handler.fluid_rhs_address(), memory_handler.nvariables(), memory_handler.stride())
      , m_persistent_rhs(memory_handler.fluid_rhs_address(), memory_handler.nvariables(), memory_handler.stride())
      , m_ijk(ijk)
      , m_indexes(memory_handler.indexes().fluid_indexes())
      , m_ep_g(ep_g)
      , m_ro_g(ro_g)
      , m_vel_g(vel_g)
      , m_T_g(T_g)
      , m_h_g(h_g)
      , m_X_gk(X_gk)
      //, m_p_g(p_g)
      , m_cell_volume(cell_volume)
      , m_vfrac(vfrac)
      , m_chem_txfr_g(chem_txfr_g)
      , m_chem_txfr_idxs(chem_txfr_idxs)
      , m_report_mass_balance(report_mass_balance)
      , m_fluid_prod_ptr(fluid_prod_ptr)
    {}

    /**
     * \brief Get indexes of current fluid cell
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::IntVect const& cell_ijk () const { return m_ijk; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_rhs ()
    {
      m_rhs[ChemFluidIndexes::density()] = 0.;
      m_rhs[ChemFluidIndexes::temperature()] = 0.;
      m_rhs[ChemFluidIndexes::specific_enthalpy()] = 0.;

      for (int n_g(0); n_g < params.nspecies(); n_g++)
        m_rhs[ChemFluidIndexes::mass_fractions()+n_g] = 0.;
    }

  private:
    /**
     * Change this class pointer to the vector of unknowns. This is useful, for
     * example, when integrating the chemistry ODE system with iterative
     * integrators that use temporary variables which are passed to the
     * MFIXActualRHS class.
     *
     * \param y is the pointer to the memory area to be copied
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set_y_ptr (amrex::Real* const& y) { m_y.dataPtr() = y; }

    /**
     * Change this class pointer to the vector of the ODE right-hand side. This
     * is useful, for example, when integrating the chemistry ODE system with
     * iterative integrators that use temporary variables which are passed to
     * the MFIXActualRHS class.
     *
     * \param y_rhs is the pointer to the memory area to be copied
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set_y_rhs_ptr (amrex::Real* const& y_rhs) { m_rhs.dataPtr() = y_rhs; }

    /**
     * Reset the pointer to the vector of unknowns to the correct memory
     * space in case it was modified.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_y_ptr () { m_y.dataPtr() = m_persistent_y.dataPtr(); }

    /**
     * Reset the pointer to the vector of unknowns to the correct memory
     * space in case it was modified.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_y_rhs_ptr () { m_rhs.dataPtr() = m_persistent_rhs.dataPtr(); }

  public:
    /**
     * \brief Get the fluid shear viscosity at the current fluid cell.
     *
     * \param T_g is the fluid temperature
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real shear_viscosity (const amrex::Real& T_g) const
    {
      return m_props.molViscosity(T_g, m_ijk, m_X_gk);
    }

    /**
     * \brief Get the fluid shear viscosity at the current fluid cell.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real shear_viscosity () const
    {
      return shear_viscosity(temperature());
    }

    /**
     * \brief Get the fluid volume fraction at the current fluid cell.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real volume_fraction () const
    {
      return m_ep_g(m_ijk);
    }

    /**
     * \brief Get the fluid density at the current fluid cell.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real density () const
    {
      AMREX_ASSERT(m_y[ChemFluidIndexes::density()] > 0.);
      return m_y[ChemFluidIndexes::density()];
    }

    /**
     * \brief Get the fluid velocity at the current fluid cell.
     */
    amrex::RealVect velocity () const
    {
      return {AMREX_D_DECL(velocity(0), velocity(1), velocity(2))};
    }

    /**
     * \brief Get a component of the fluid velocity at the current fluid cell.
     *
     * \param n is the component of the fluid velocity
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real velocity (const int& n) const
    {
      AMREX_ASSERT(n >= 0 && n < AMREX_SPACEDIM);
      return m_vel_g(m_ijk, n);
    }

    /**
     * \brief Get the fluid thermodynamic pressure at the current fluid cell.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real thermodynamic_pressure () const
    {
      return MFIXFluidPhase::R*temperature()*density()*average_molar_mass();
    }

    /**
     * \brief Get the fluid temperature at the current fluid cell
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real temperature () const
    {
      return m_y[ChemFluidIndexes::temperature()];
    }

    /**
     * \brief Get the fluid specific enthalpy at the current fluid cell.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real specific_enthalpy () const
    {
      return m_y[ChemFluidIndexes::specific_enthalpy()];
    }

  public:
    /**
     * \brief Get the unique index associated to the given fluid species name
     *
     * Returns an index (between 0 and the number of fluid phase species minus
     * one) that is uniquely associated to the fluid species corresponding
     * names. If the requested species name is not found among those provided in
     * the inputs file, then it returns -1
     *
     * \param species_name is the string defining the fluid species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int species_idx (const char* const& species_name) const
    {
      const int n_g = params.species_idx(species_name);
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return n_g;
    }

    /**
     * \brief Get the fluid mass fraction at the current fluid cell for the
     * given species index.
     *
     * \param n_g is the species index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mass_fraction (const int& n_g) const
    {
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return m_y[ChemFluidIndexes::mass_fractions()+n_g];
    }

    /**
     * \brief Get the fluid mass fraction at the current fluid cell for the
     * given species name.
     *
     * \param species_name is the species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mass_fraction (const char* const& species_name) const
    {
      const int n_g = species_idx(species_name);
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return m_y[ChemFluidIndexes::mass_fractions()+n_g];
    }

    /**
     * \brief Get the fluid molar concentration at the current fluid cell for
     * the given species index.
     *
     * \param n_g is the species index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molar_concentration (const int& n_g) const
    {
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return density() * mass_fraction(n_g) / molar_mass(n_g);
    }

    /**
     * \brief Get the fluid molar concentration at the current fluid cell for
     * the given species name.
     *
     * \param species_name is the species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molar_concentration (const char* const& species_name) const
    {
      const int n_g = species_idx(species_name);
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return molar_concentration(n_g);
    }

    /**
     * \brief Get the fluid molar mass at the current fluid cell for the given
     * species index.
     *
     * \param n_g is the species index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molar_mass (const int& n_g) const
    {
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return params.get_MW_gk(n_g);
    }

    /**
     * \brief Get the fluid molar mass at the current fluid cell for the given
     * species name.
     *
     * \param species_name is the species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molar_mass (const char* const& species_name) const
    {
      const int n_g = species_idx(species_name);
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return molar_mass(n_g);
    }

    /**
     * \brief Get the fluid average molar mass at the current fluid cell
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real average_molar_mass () const
    {
      amrex::Real avg_MW_g(0.);
      for(int n_g(0); n_g < params.nspecies(); ++n_g) {
        avg_MW_g += mass_fraction(n_g) / molar_mass(n_g);
      }
      return 1./avg_MW_g;
    }

    /**
     * \brief Get the fluid mole fraction at the current fluid cell for the
     * given species index.
     *
     * \param n_g is the species index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mole_fraction (const int& n_g) const
    {
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return (mass_fraction(n_g) / molar_mass(n_g)) * average_molar_mass();
    }

    /**
     * \brief Get the fluid mass fraction at the current fluid cell for the
     * given species name.
     *
     * \param species_name is the species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mole_fraction (const char* const& species_name) const
    {
      const int n_g = species_idx(species_name);
      AMREX_ASSERT(n_g >= 0 && n_g < params.nspecies());
      return mole_fraction(n_g);
    }

  private:
    /**
     * \brief Initialize the vector of unknowns of the ODE to be integrated
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void init_vars ()
    {
      // copy the species densities
      m_y[ChemFluidIndexes::density()] = m_ro_g(m_ijk);

      for (int n_g(0); n_g < params.nspecies(); n_g++) {
        m_y[ChemFluidIndexes::mass_fractions()+n_g] = m_X_gk(m_ijk, n_g);
      }

      // copy the temperature
      m_y[ChemFluidIndexes::temperature()] = m_T_g(m_ijk);

      // copy the enthalpy
      m_y[ChemFluidIndexes::specific_enthalpy()] = m_h_g(m_ijk);
    }

    /**
     * \brief Use the new fluid variables to compute zero-order fluid chemical
     * transfer quantities.
     *
     * Update the Array4 containing the fluid transfer quantities due to
     * chemical reactions by doing a zero-order finite difference in time
     * between the new fluid variables computed through the ODE integration and
     * the old fluid quantities before the ODE integration.
     *
     * \param dt is the timestep used for the ODE integration
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void update_chem_txfr (const amrex::Real& dt)
    {
      // Get the fluid volume fraction
      const amrex::Real& ep_g = m_ep_g(m_ijk);
      // Get the fluid density before the ODE integration
      const amrex::Real& ro_g_old = m_ro_g(m_ijk);

      // Loop over the fluid species
      for (int n_g(0); n_g < params.nspecies(); n_g++) {
        // Get the fluid n-th species density before the ODE integration
        const amrex::Real ro_gk_old = m_X_gk(m_ijk, n_g) * ro_g_old;
        // Get the fluid n-th species density after the ODE integration
        const amrex::Real& ro_gk_new = mass_fraction(n_g) * density();

        // Compute the difference between the new and the old n-th species
        // density, and set it to zero if it is smaller than a given tolerance
        const amrex::Real diff = amrex::almostEqual(ro_gk_new, ro_gk_old) ?
          0. : (ro_gk_new - ro_gk_old);

        // Copy the zero-order transfer quantity approximation for the transfer
        // rate of the fluid n-th species density into the proper Array4 space
        m_chem_txfr_g(m_ijk, m_chem_txfr_idxs.chem_ro_gk+n_g) = (ep_g*diff) / dt;

        if (m_report_mass_balance) {
          const amrex::Real effective_cell_volume = m_vfrac*m_cell_volume;
          const amrex::Real fluid_mass_prod = effective_cell_volume*ep_g*diff;
          amrex::HostDevice::Atomic::Add(&m_fluid_prod_ptr[n_g], fluid_mass_prod);
        }
      }

      // Get the fluid enthalpy density before the ODE integration
      const amrex::Real ro_h_g_old = ro_g_old * m_h_g(m_ijk);
      // Get the fluid enthalpy density after the ODE integration
      const amrex::Real& ro_h_g_new = density() * specific_enthalpy();

      // Compute the difference between the new and the old fluid enthalpy
      // density, and set it to zero if it is smaller than a given tolerance
      const amrex::Real diff = amrex::almostEqual(ro_h_g_new, ro_h_g_old)?
        0. : (ro_h_g_new - ro_h_g_old);

      // Copy the zero-order transfer quantity approximation for the transfer
      // rate of the fluid enthalpy density into the proper Array4 space
      m_chem_txfr_g(m_ijk, m_chem_txfr_idxs.chem_h) = (ep_g*diff) / dt;
    }

    // A reference to the MFIXFluidParms class that will be used to compute
    // necessary fluid phase quantities
    const MFIXFluidParms& params;
    const ThermoPropertyData& m_props;

    ChemVector<amrex::Real> m_y;
    ChemVector<amrex::Real> const m_persistent_y;
    ChemVector<amrex::Real> m_rhs;
    ChemVector<amrex::Real> const m_persistent_rhs;
    const amrex::IntVect& m_ijk;
    const ChemMemoryIndexes::ChemFluidIndexes& m_indexes;

    amrex::Array4<const amrex::Real> const& m_ep_g;
    amrex::Array4<const amrex::Real> const& m_ro_g;
    amrex::Array4<const amrex::Real> const& m_vel_g;
    amrex::Array4<const amrex::Real> const& m_T_g;
    amrex::Array4<const amrex::Real> const& m_h_g;
    amrex::Array4<const amrex::Real> const& m_X_gk;
    //const amrex::Real& m_p_g;
    const amrex::Real& m_cell_volume;
    const amrex::Real& m_vfrac;
    amrex::Array4<amrex::Real> const& m_chem_txfr_g;
    const InterphaseChemTxfrIndexes& m_chem_txfr_idxs;

    const int m_report_mass_balance;
    amrex::Real* m_fluid_prod_ptr;
};


/**
 * \brief Class that stores all the solids phase data needed in chemistry
 *
 * This class contains all the solids phase data that can be useful to
 * define the reactions rates functions, to compute the chemistry ODE right-hand
 * side, and to integrate the chemistry ODE. In addition, this class stores the
 * pointers to the memory areas that will be used to save the ODE unknowns and
 * right-hand side.
 */
class MFIXSolidsData
{
  // Define classes MFIXActualRHS and MaterialInterface friends so they can
  // access the private members of this class
  friend class MFIXActualRHS;
  friend class MaterialInterface;

  public:
    // Aliases
    using ChemSolidsIndexes = ChemMemoryIndexes::ChemSolidsIndexes;
    using PType = amrex::Particle<0,0>;
    using PTileType = amrex::ParticleTileData<PType,SoArealData::count,SoAintData::count>;

    /**
     * \brief Constructs the MFIXSolidsData class
     *
     * Constructor of the class
     *
     * \param memory_handler is the class that contains all the information
     * about the memory areas that this MFIXSolidsData class can use to store
     * ODE unknowns and right-hand side data
     * \param solids_parms is the MFIXSolidsParms class reference
     * \param ijk are the indexes for current fluid cell where the particle is
     * contained
     * \param parts_nb_array is the number of particles in each fluid cell
     * \param parts_idxs_array is the first global index of particles in each
     * fluid cell
     * \param parts_idxs_ptr are the global indexes of particles vectorized in
     * order for each fluid cell
     * \param ptile_data is the MFIXParticleContainer particles data
     * \param runtime_indexes are the indexes to properly access particles
     * runtime variables
     */
    AMREX_GPU_HOST_DEVICE
    MFIXSolidsData (DeviceChemMemoryHandler& memory_handler,
                    const MFIXSolidsParms& solids_parms,
                    const ThermoPropertyData& solids_props,
                    const amrex::IntVect& ijk,
                    const amrex::Array4<const int>& parts_nb_array,
                    const amrex::Array4<const int>& parts_idxs_array,
                    const int* parts_idxs_ptr,
                    const PTileType& ptile_data,
                    const runtimeRealData& runtime_indexes,
                    const int solids_update_type,
                    const IntegratorInputs::solids_thresholds& solids_thresholds,
                    IntegrationStatus& integration_status,
                    const int report_mass_balance,
                    amrex::Real* solids_prod_ptr)
      : params(solids_parms)
      , m_props(solids_props)
      , m_y(memory_handler.solids_vars_address(), memory_handler.nvariables(), memory_handler.stride())
      , m_persistent_y(memory_handler.solids_vars_address(), memory_handler.nvariables(), memory_handler.stride())
      , m_rhs(memory_handler.solids_rhs_address(), memory_handler.nvariables(), memory_handler.stride())
      , m_persistent_rhs(memory_handler.solids_rhs_address(), memory_handler.nvariables(), memory_handler.stride())
      , m_p_id(-1)
      , m_ijk(ijk)
      , m_indexes(memory_handler.indexes().solids_indexes())
      , m_parts_nb_array(parts_nb_array)
      , m_parts_idxs_array(parts_idxs_array)
      , m_parts_idxs_ptr(parts_idxs_ptr)
      , m_ptile_data(ptile_data)
      , m_runtime_indexes(runtime_indexes)
      , m_solids_update_type(solids_update_type)
      , m_solids_thresholds(solids_thresholds)
      , m_integration_status(integration_status)
      , m_report_mass_balance(report_mass_balance)
      , m_solids_prod_ptr(solids_prod_ptr)
    {}

  private:
    /**
     * \brief Update current particle ID to access ptile_data information at the
     * correct position
     *
     * \param n is the current number of this particle within the fluid cell
     * where it is contained
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void update_p_id (const int& n)
    {
      // Get the particles global index for the current fluid cell
      const amrex::Long cell_idx = m_parts_idxs_array(m_ijk);
      // Get this particle global index
      const amrex::Long vect_idx = cell_idx - n;
      // Use this particle global index to obtain this particle global ID
      m_p_id = m_parts_idxs_ptr[vect_idx];
    }

  public:
    /**
     * \brief Get the total number of particles in the current fluid cell
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int particles_nb_in_cell () const
    {
      return m_parts_nb_array(m_ijk);
    }

    // TODO
    // consider computing the nb of particles in a cell by subtracting the
    // global indexes between the current ijk cell and the "previous" ijk-1
    // cell, paying attention to consider as 0 the global index in the "-1" cell
    // This would allow to save the memory space for the MultiFab containing the
    // nb of particles per cell

    /**
     * \brief Get the total number of particles in the fluid cell corresponding
     * to the gven indexes
     *
     * \param ijk are the cell indexes
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int particles_nb_in_cell (const amrex::IntVect& ijk) const
    {
      return m_parts_nb_array(ijk);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_rhs ()
    {
      m_rhs[ChemSolidsIndexes::density()] = 0.;
      m_rhs[ChemSolidsIndexes::radius()] = 0.;
      m_rhs[ChemSolidsIndexes::temperature()] = 0.;
      m_rhs[ChemSolidsIndexes::specific_enthalpy()] = 0.;

      for (int n_s(0); n_s < params.nspecies(); n_s++)
        m_rhs[ChemSolidsIndexes::mass_fractions()+n_s] = 0.;
    }

  private:
    /**
     * Change this class pointer to the vector of unknowns. This is useful, for
     * example, when integrating the chemistry ODE system with iterative
     * integrators that use temporary variables which are passed to the
     * MFIXActualRHS class.
     *
     * \param y is the pointer to be copied
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set_y_ptr (amrex::Real* const& y) { m_y.dataPtr() = y; }

    /**
     * Change this class pointer to the vector of the ODE right-hand side. This
     * is useful, for example, when integrating the chemistry ODE system with
     * iterative integrators that use temporary variables which are passed to
     * the MFIXActualRHS class.
     *
     * \param y_rhs is the pointer to be copied
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set_y_rhs_ptr (amrex::Real* const& y_rhs) { m_rhs.dataPtr() = y_rhs; }

    /**
     * Reset the pointer to the vector of unknowns to the correct memory space
     * in case it was modified.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_y_ptr () { m_y.dataPtr() = m_persistent_y.dataPtr(); }

    /**
     * Reset the pointer to the vector of unknowns to the correct memory space
     * in case it was modified.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_y_rhs_ptr () { m_rhs.dataPtr() = m_persistent_rhs.dataPtr(); }

  public:
    /**
     * \brief Get this particle's ID.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Long id () const
    {
      return amrex::Long(m_ptile_data.id(m_p_id));
    }

    /**
     * \brief Get this particle's cpu.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int cpu () const
    {
      return int(m_ptile_data.cpu(m_p_id));
    }

    /**
     * \brief Get this particle's position component.
     *
     * \brief n the component index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real position (const int& n) const
    {
      AMREX_ASSERT(n >= 0 && n < AMREX_SPACEDIM);
      return m_ptile_data.pos(n, m_p_id);
    }

    /**
     * \brief Get this particle's position.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect position () const
    {
      return {AMREX_D_DECL(position(0), position(1), position(2))};
    }

    /**
     * \brief Get this particle's radius.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real radius () const
    {
      if ((m_y[ChemSolidsIndexes::radius()] < 0.) ||
          (amrex::almostEqual(m_y[ChemSolidsIndexes::radius()], 0.))) {

        // set the integration status to invalid
        m_integration_status = IntegrationStatus::IERR_INVALID_SOLIDS_MASS;
        // return quiet NaN
        return std::numeric_limits<amrex::Real>::quiet_NaN();
      }

      return m_y[ChemSolidsIndexes::radius()];
    }

    /**
     * \brief Get this particle's volume.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real volume () const
    {
      return SoArealData::volume(radius());
    }

    /**
     * \brief Get this particle's mass.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mass () const
    {
      return density()*volume();
    }

    /**
     * \brief Get this particle's density.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real density () const
    {
      if ((m_y[ChemSolidsIndexes::density()] < 0.) ||
          (amrex::almostEqual(m_y[ChemSolidsIndexes::density()], 0.))) {

        // set the integration status to invalid
        m_integration_status = IntegrationStatus::IERR_INVALID_SOLIDS_MASS;
        // return quiet NaN
        return std::numeric_limits<amrex::Real>::quiet_NaN();
      }

      return m_y[ChemSolidsIndexes::density()];
    }

    /**
     * \brief Get this particle's inverse of the momentum of inertia.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real oneOverI () const
    {
      return SoArealData::oneOverI(radius(), mass());
    }

    /**
     * \brief Get this particle's statistical weight.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real statistical_weight () const
    {
      return m_ptile_data.rdata(SoArealData::statwt)[m_p_id];
    }

    /**
     * \brief Get this particle's velocity n-th component.
     *
     * \param n is the velocity component
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real velocity (const int& n) const
    {
      AMREX_ASSERT(n >= 0 && n < AMREX_SPACEDIM);
      return m_ptile_data.rdata(SoArealData::velx+n)[m_p_id];
    }

    /**
     * \brief Get this particle's velocity.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect velocity () const
    {
      return {AMREX_D_DECL(velocity(0), velocity(1), velocity(2))};
    }

    /**
     * \brief Get this particle's temperature.
     *
     * This particle's temperature is computed from the particle's specific
     * enthalpy using the Newton-Raphson algorithm because the particle's
     * specific enthalpy changes during the chemistry ODE integration and during
     * the computation of chemical reactions.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real temperature () const
    {
      return m_y[ChemSolidsIndexes::temperature()];
    }

    /**
     * \brief Get this particle's specific enthalpy.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real specific_enthalpy () const
    {
      return m_y[ChemSolidsIndexes::specific_enthalpy()];
    }

    /**
     * \brief Get the unique index associated to the given solids species name.
     *
     * Returns an index (between 0 and the number of solids phase species minus
     * one) that is uniquely associated to the solids species corresponding
     * names. If the requested species name is not found among those provided in
     * the inputs file, then it returns -1
     *
     * \param species_name is the string defining the solids species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int species_idx (const char* const& species_name) const
    {
      const int n_s = params.species_idx(species_name);
      AMREX_ASSERT(n_s >= 0 && n_s < params.nspecies());
      return n_s;
    }

    /**
     * \brief Get this particle's mass fraction for the given species index.
     *
     * \param n_s is the solids species index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mass_fraction (const int& n_s) const
    {
      AMREX_ASSERT(n_s >=0 && n_s < params.nspecies());
      return m_y[ChemSolidsIndexes::mass_fractions()+n_s];
    }

    /**
     * \brief Get this particle's mass fraction for the given species name.
     *
     * \param species_name is the solids species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mass_fraction (const char* const& species_name) const
    {
      const int n_s = species_idx(species_name);
      AMREX_ASSERT(n_s >=0 && n_s < params.nspecies());
      return mass_fraction(n_s);
    }

    /**
     * \brief Get this particle's molar concentration for the given species
     * index.
     *
     * \param n_s is the solids species index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molar_concentration (const int& n_s) const
    {
      AMREX_ASSERT(n_s >=0 && n_s < params.nspecies());
      return density() * mass_fraction(n_s) / molar_mass(n_s);
    }

    /**
     * \brief Get this particle's molar concentration for the given species
     * name.
     *
     * \param species_name is the solids species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molar_concentration (const char* const& species_name) const
    {
      const int n_s = species_idx(species_name);
      AMREX_ASSERT(n_s >=0 && n_s < params.nspecies());
      return molar_concentration(n_s);
    }

    /**
     * \brief Get this particle's molar mass for the given species index.
     *
     * \param n_s is the solids species index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molar_mass (const int& n_s) const
    {
      AMREX_ASSERT(n_s >= 0 && n_s < params.nspecies());
      return params.get_MW_sn(n_s);
    }

    /**
     * \brief Get the solids molar concentration at the current solids cell for
     * the given species name.
     *
     * \param species_name is the solids species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real molar_mass (const char* const& species_name) const
    {
      const int n_s = species_idx(species_name);
      AMREX_ASSERT(n_s >= 0 && n_s < params.nspecies());
      return molar_mass(n_s);
    }

    /**
     * \brief Get the particle's average molar mass
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real average_molar_mass () const
    {
      amrex::Real avg_MW_s(0.);
      for (int n_s(0); n_s < params.nspecies(); ++n_s) {
        avg_MW_s += mass_fraction(n_s) / molar_mass(n_s);
      }
      return 1./avg_MW_s;
    }

    /**
     * \brief Get the particle's mole fraction for the given species index.
     *
     * \param n_s is the solids species index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mole_fraction (const int& n_s) const
    {
      AMREX_ASSERT(n_s >= 0 && n_s < params.nspecies());
      return (mass_fraction(n_s) / molar_mass(n_s)) * average_molar_mass();
    }

    /**
     * \brief Get the particle's mole fraction for the given species name.
     *
     * \param species_name is the solids species name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real mole_fraction (const char* const& species_name) const
    {
      const int n_s = species_idx(species_name);
      AMREX_ASSERT(n_s >= 0 && n_s < params.nspecies());
      return mole_fraction(n_s);
    }

  private:

    /**
     * Initialize the vector of unknowns of the ODE to be integrated
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool constant_density ()
    {
      return (m_solids_update_type == SolidsUpdateType::ConstantDensity);
    }

    /**
     * Initialize the vector of unknowns of the ODE to be integrated
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool constant_volume ()
    {
      return (m_solids_update_type == SolidsUpdateType::ConstantVolume);
    }

    /**
     * Initialize the vector of unknowns of the ODE to be integrated
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void init_vars ()
    {
      m_y[ChemSolidsIndexes::density()] = m_ptile_data.rdata(SoArealData::density)[m_p_id];
      m_y[ChemSolidsIndexes::radius()] = m_ptile_data.rdata(SoArealData::radius)[m_p_id];

      // Initialize the solids species mass fractions
      for (int n_s(0); n_s < params.nspecies(); n_s++) {
        m_y[ChemSolidsIndexes::mass_fractions()+n_s] =
          m_ptile_data.m_runtime_rdata[m_runtime_indexes.X_sn+n_s][m_p_id];
      }

      // Initialize the solids temperature
      m_y[ChemSolidsIndexes::temperature()] = m_ptile_data.rdata(SoArealData::temperature)[m_p_id];

      // Initialize the solids specific enthalpy
      {
        amrex::Real h_s(0.);
        const amrex::Real T_s = m_ptile_data.rdata(SoArealData::temperature)[m_p_id];

        for (int n_s(0); n_s < params.nspecies(); ++n_s)
          h_s += m_ptile_data.m_runtime_rdata[m_runtime_indexes.X_sn+n_s][m_p_id] *
            m_props.enthalpy(n_s, T_s);

        m_y[ChemSolidsIndexes::specific_enthalpy()] = h_s;
      }
    }

    class SpecificEnthalpyResidue
    {
      public:
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SpecificEnthalpyResidue (ChemVector<amrex::Real> const& y,
                                 const ThermoPropertyData& props,
                                 const int nspecies)
          : m_y(y)
          , m_props(props)
          , m_nspecies(nspecies)
        {}

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::Real operator() (amrex::Real const& T_s) const
        {
          amrex::Real H_s(0.);
          for (int n_s(0); n_s < m_nspecies; ++n_s) {
            const amrex::Real h_sn = m_props.enthalpy(n_s, T_s);
            H_s += m_y[ChemSolidsIndexes::mass_fractions()+n_s]*h_sn;
          }

          return H_s - m_y[ChemSolidsIndexes::specific_enthalpy()];
        }

      private:
        ChemVector<amrex::Real> const& m_y;
        const ThermoPropertyData& m_props;
        const int m_nspecies;
    };

    class SpecificEnthalpyResidueGradient
    {
      public:
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        SpecificEnthalpyResidueGradient (ChemVector<amrex::Real> const& y,
                                         const ThermoPropertyData& props,
                                         const int nspecies)
          : m_y(y)
          , m_props(props)
          , m_nspecies(nspecies)
        {}

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::Real operator() (amrex::Real const& T_s) const
        {
          amrex::Real partial_H_s(0.);
          for (int n_s(0); n_s < m_nspecies; ++n_s) {
            const amrex::Real partial_h_sn = m_props.specificHeat(n_s, T_s);
            partial_H_s += m_y[ChemSolidsIndexes::mass_fractions()+n_s]*partial_h_sn;
          }

          return partial_H_s;
        }

      private:
        ChemVector<amrex::Real> const& m_y;
        const ThermoPropertyData& m_props;
        const int m_nspecies;
    };

    /**
     * \brief Update particle's variables from the vector of ODEs unknowns
     * after the ODE integration
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void update_variables ()
    {
      const amrex::Real p_radius_old = m_ptile_data.rdata(SoArealData::radius)[m_p_id];
      const amrex::Real p_volume_old = SoArealData::volume(p_radius_old);
      const amrex::Real p_density_old = m_ptile_data.rdata(SoArealData::density)[m_p_id];
      const amrex::Real p_mass_old = p_density_old * p_volume_old;

      const amrex::Real p_mass_new = mass();

      if (p_mass_new < m_solids_thresholds.mass ||
          m_integration_status == IntegrationStatus::IERR_INVALID_SOLIDS_MASS) {
        m_ptile_data.id(m_p_id).make_invalid();
        return;
      }

      if (m_solids_update_type == SolidsUpdateType::ConstantDensity) {

        amrex::Real p_radius_new = radius();

        if (p_radius_new < m_solids_thresholds.radius) {
          m_ptile_data.id(m_p_id).make_invalid();
          return;
        } else {
          m_ptile_data.rdata(SoArealData::radius)[m_p_id] = p_radius_new;
        }

      } else if (m_solids_update_type == SolidsUpdateType::ConstantVolume) {

        amrex::Real p_density_new = density();

        if (p_density_new < m_solids_thresholds.density) {
          m_ptile_data.id(m_p_id).make_invalid();
          return;
        } else {
          m_ptile_data.rdata(SoArealData::density)[m_p_id] = p_density_new;
        }

      } else if (m_solids_update_type == SolidsUpdateType::ConstantSpeciesDensities) {

        amrex::Real p_density_inv(0.);
        for (int n = 0; n < params.nspecies(); n++) {
          // note: if allowing infinite densities (denoted by negative),
          // uncomment the following line:
          //if (params.get_density(n) < -std::numeric_limits<amrex::Real>::min()) continue;
          p_density_inv += std::max(mass_fraction(n), 0.) / params.get_density(n);
        }

        amrex::Real p_density_new = 1./p_density_inv;
        amrex::Real p_volume_new = p_mass_new * p_density_inv;
        amrex::Real p_radius_new = std::cbrt((3.*p_volume_new)/(4.*M_PI));

        if (p_radius_new  < m_solids_thresholds.radius ||
            p_density_new < m_solids_thresholds.density) {
          m_ptile_data.id(m_p_id).make_invalid();
          return;
        } else {
          m_ptile_data.rdata(SoArealData::radius)[m_p_id] = p_radius_new;
          m_ptile_data.rdata(SoArealData::density)[m_p_id] = p_density_new;
        }

      } else {

        AMREX_ALWAYS_ASSERT(false);

      }

      // Compute temperature from specific enthalpy using a Newton-Raphson method
      {
        // Local variable to store final result for the particle's temperature. It
        // is initialized with "old" particle's temperature. This value is the
        // initial value used by the Newton algorithm.
        amrex::Real result(m_ptile_data.rdata(SoArealData::temperature)[m_p_id]);

        // Define the residue of the equation to be solved with the Newton
        // algorithm
        SpecificEnthalpyResidue residue(m_y, m_props, params.nspecies());

        // Define the gradient of the equation to be solved with the Newton
        // algorithm
        SpecificEnthalpyResidueGradient gradient(m_y, m_props, params.nspecies());

        // Call the Newton solver
        auto output = Solvers::Newton::solve(result, residue, gradient, 1.e-8, 1.e-8, 100);

        // Check if the Newton solver has converged
        if (output.iterations == -1)
          amrex::Abort("Newton solver did not converge");

        m_ptile_data.rdata(SoArealData::temperature)[m_p_id] = result;
      }

      //
      for (int n_s(0); n_s < params.nspecies(); n_s++) {
        const amrex::Real mass_sn_old =
          m_ptile_data.m_runtime_rdata[m_runtime_indexes.X_sn+n_s][m_p_id] * p_mass_old;

        m_ptile_data.m_runtime_rdata[m_runtime_indexes.X_sn+n_s][m_p_id] = mass_fraction(n_s);

        const amrex::Real mass_sn_new = mass_fraction(n_s) * p_mass_new;
        const amrex::Real diff = amrex::almostEqual(mass_sn_new, mass_sn_old) ?
          0 : (mass_sn_new - mass_sn_old);

        if (m_report_mass_balance) {
          amrex::HostDevice::Atomic::Add(&m_solids_prod_ptr[n_s], statistical_weight()*diff);
        }
      }
    }

    const MFIXSolidsParms& params;
    const ThermoPropertyData& m_props;

    ChemVector<amrex::Real> m_y;
    ChemVector<amrex::Real> const m_persistent_y;
    ChemVector<amrex::Real> m_rhs;
    ChemVector<amrex::Real> const m_persistent_rhs;
    int m_p_id;
    const amrex::IntVect& m_ijk;
    const ChemMemoryIndexes::ChemSolidsIndexes& m_indexes;

    const amrex::Array4<const int>& m_parts_nb_array;
    const amrex::Array4<const int>& m_parts_idxs_array;
    const int* m_parts_idxs_ptr;
    const PTileType& m_ptile_data;
    const runtimeRealData& m_runtime_indexes;

    const int m_solids_update_type;
    const IntegratorInputs::solids_thresholds& m_solids_thresholds;
    IntegrationStatus& m_integration_status;

    const int m_report_mass_balance;
    amrex::Real* m_solids_prod_ptr;
};


/**
 * \brief Class that stores all the reactions data needed in chemistry
 *
 * This class contains all the reactions data that can be useful to define the
 * reactions rates functions, to compute the chemistry ODE right-hand side, and
 * to integrate the chemistry ODE. In addition, this class stores the pointers
 * to the memory areas that will be used to save the eulerian and lagrangian
 * reactions rates.
 */
class MFIXReactionsData
{
    // Define class MaterialInterface friends so it can access the private
    // members of this class
    friend class MaterialInterface;
    friend class MFIXActualRHS;

  public:
    /**
     * \brief Constructs the MFIXReactionsData class
     *
     * Constructor of the class
     *
     * \param memory_handler is the class that contains all the information
     * about the memory areas that this MFIXReactionsData class can use to store
     * reactions rates
     * \param reactions_parms is the MFIXReactionsParms class reference
     * \param time is the current time
     * \param dt is the timestep used for ODE integration
     * \param cell_volume is the current fluid cell volume
     * \param vfrac is the current fluid cell volume fraction
     */
    AMREX_GPU_HOST_DEVICE
    MFIXReactionsData (DeviceChemMemoryHandler& memory_handler,
                       const MFIXReactionsParms& reactions_parms,
                       const amrex::Real& time,
                       const amrex::Real& dt,
                       const amrex::Real& cell_volume,
                       const amrex::Real& vfrac)
      : params(reactions_parms)
      , m_rates(memory_handler.rates_address(), reactions_parms.m_nreactions, memory_handler.stride())
      , m_time(time)
      , m_dt(dt)
      , m_cell_volume(cell_volume)
      , m_vfrac(vfrac)
    {
      this->reset_rates();
    }

  private:
    /**
     * \brief Member function to reset the reactions rates to zero.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset_rates ()
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXReactionsData::reset_rates()");
#endif

      for (int q(0); q < params.m_nreactions; ++q)
        m_rates[q] = 0.;
    }

  public:
    /**
     * \brief Get the number of reactions.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int const& reactions_number () const { return params.m_nreactions; }

    /**
     * \brief Get the unique index associated to the given chemical reaction
     * name.
     *
     * Returns an index (between 0 and the number of chemical reactions minus
     * one) that is uniquely associated to the given type of chemical reactions
     * corresponding names. If the requested reaction not found among those
     * provided in the inputs file, then it returns -1.
     *
     * \param reaction_name is the string defining the reaction name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int reaction_idx (const char* const& reaction_name) const
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXReactionsData::reaction_idx()");
#endif

      const int q = params.reaction_idx(reaction_name);
      AMREX_ASSERT(q >= 0 && q < params.m_nreactions);
      return q;
    }

    /**
     * \brief Returns a reference to the reaction rate with the given index.
     *
     * \param q is the reaction index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real& rate (const int& q)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXReactionsData::set_rate(int)");
#endif

      return m_rates[q];
    }

    /**
     * \brief Returns a reference to the reaction rate with the given name.
     *
     * \param reaction_name is the reaction name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real& rate (const char* const& reaction_name)
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXReactionsData::set_rate(char*)");
#endif

      const int q = reaction_idx(reaction_name);
      AMREX_ASSERT(q >= 0 && q < params.m_nreactions);
      return m_rates[q];
    }

    /**
     * \brief Returns the reaction rate with the given index.
     *
     * \param q is the reaction index
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real const& rate (const int& q) const
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXReactionsData::get_rate(int)");
#endif

      return m_rates[q];
    }

    /**
     * \brief Returns the reaction rate with the given name.
     *
     * \param reaction_name is the reaction name
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real const& rate (const char* const& reaction_name) const
    {
#ifndef AMREX_USE_GPU
      BL_PROFILE("MFIXReactionsData::get_rate(char*)");
#endif

      const int q = reaction_idx(reaction_name);
      AMREX_ASSERT(q >= 0 && q < params.m_nreactions);
      return m_rates[q];
    }

    /**
     * \brief Get the current time.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real const& time () const { return m_time; }

    /**
     * \brief Get the ODE integration timestep.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real const& dt () const { return m_dt; }

    /**
     * \brief Get the fluid cell geometric volume.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real const& cell_volume () const { return m_cell_volume; }

    /**
     * \brief Get the fluid cell volume fraction.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real const& vfrac () const { return m_vfrac; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int has_fluid_and_solids (const int q) const
    {
      return params.has_fluid_and_solids(q);
    }

  private:
    const MFIXReactionsParms& params;

    ChemVector<amrex::Real> m_rates;
    const amrex::Real& m_time;
    const amrex::Real& m_dt;
    const amrex::Real& m_cell_volume;
    const amrex::Real& m_vfrac;
};


#endif
