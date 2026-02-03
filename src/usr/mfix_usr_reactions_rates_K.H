#ifndef MFIX_USR_REACTIONS_RATES_K_H
#define MFIX_USR_REACTIONS_RATES_K_H

#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>


/**
 * \brief The class for user-defined lagrangian reaction rates
 */
class LagrangianReactionRates
{
  public:
    /**
     * \brief Default constructor
     */
    AMREX_GPU_HOST_DEVICE
    LagrangianReactionRates () = default;

    /**
     * \brief The class call operator that computes the lagrangian reaction
     * rates
     *
     * \param reactions the class that contains the reactions data for the
     * chemistry computation and storing the rates
     * \param fluid the class that contains the fluid data for the chemistry
     * computation
     * \param solids the class that contains the solids data for the
     * chemistrycomputation
     *
     *  Reactions class gives access to the following reactions quantities:
     *   - rate(q) with q in [0,reactions_nb]
     *   - rate("reaction_name")
     *   - time()
     *   - dt()
     *   - cell_volume()
     *
     *  Fluid class gives access to the following fluid quantities:
     *   - shear_viscosity()
     *   - volume_fraction()
     *   - density()
     *   - velocity() and velocity(dir) with dir in [0,AMREX_SPACEDIM]
     *   - thermodynamic_pressure()
     *   - temperature()
     *   - specific_enthalpy()
     *   - mass_fraction(n_g) with n_g in [0,fluid_nb_species]
     *   - mass_fraction("species_name")
     *   - molar_concentration(n_g) with n_g in [0,fluid_nb_species]
     *   - molar_concentration("species_name")
     *   - molar_mass(n_g) with n_g in [0,fluid_nb_species]
     *   - molar_mass("species_name")
     *   - average_molar_mass()
     *   - mole_fraction(n_g) with n_g in [0,fluid_nb_species]
     *   - mole_fraction("species_name")
     *
     *  Solids class gives access to the following solids quantities:
     *   - position() and position(dir) with dir in [0,AMREX_SPACEDIM]
     *   - radius()
     *   - volume()
     *   - density()
     *   - mass()
     *   - oneOverI()
     *   - statistical_weight()
     *   - velocity() and velocity(dir) with dir in [0,AMREX_SPACEDIM]
     *   - temperature()
     *   - specific_heat_capacity()
     *   - mass_fraction(n_s) with n_s in [0,solids_nb_species]
     *   - mass_fraction("species_name")
     *   - molar_concentration(n_g) with n_g in [0,solids_nb_species]
     *   - molar_concentration("species_name")
     *   - molar_mass(n_g) with n_g in [0,solids_nb_species]
     *   - molar_mass("species_name")
     *   - average_molar_mass()
     *   - mole_fraction(n_g) with n_g in [0,solids_nb_species]
     *   - mole_fraction("species_name")
     */
    template <class ReactionsData, class FluidData, class SolidsData>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void operator() (ReactionsData& reactions,
                     const FluidData& fluid,
                     const SolidsData& solids) const
    {
      amrex::ignore_unused(reactions);
      amrex::ignore_unused(fluid);
      amrex::ignore_unused(solids);
    }
};


/**
 * \brief The class for user-defined eulerian reaction rates
 */
class EulerianReactionRates
{
  public:
    /**
     * \brief Default constructor
     */
    AMREX_GPU_HOST_DEVICE
    EulerianReactionRates () = default;

    /**
     * \brief The class call operator that computes the eulerian reaction rates
     *
     * \param reactions the class that contains the reactions data for the
     * chemistry computation and storing the rates
     * \param fluid the class that contains the fluid data for the chemistry
     * computation
     *
     *  Fluid class gives access to the following fluid quantities:
     *   - shear_viscosity()
     *   - volume_fraction()
     *   - density()
     *   - velocity() and velocity(dir) with dir in [0,AMREX_SPACEDIM]
     *   - thermodynamic_pressure()
     *   - temperature()
     *   - specific_enthalpy()
     *   - mass_fraction(n_g) with n_g in [0,fluid_nb_species]
     *   - mass_fraction("species_name")
     *   - molar_concentration(n_g) with n_g in [0,fluid_nb_species]
     *   - molar_concentration("species_name")
     *   - molar_mass(n_g) with n_g in [0,fluid_nb_species]
     *   - molar_mass("species_name")
     *   - average_molar_mass()
     *   - mole_fraction(n_g) with n_g in [0,fluid_nb_species]
     *   - mole_fraction("species_name")
     *
     *  Solids class gives access to the following solids quantities:
     *   - position() and position(dir) with dir in [0,AMREX_SPACEDIM]
     *   - radius()
     *   - volume()
     *   - density()
     *   - mass()
     *   - oneOverI()
     *   - statistical_weight()
     *   - velocity() and velocity(dir) with dir in [0,AMREX_SPACEDIM]
     *   - temperature()
     *   - specific_heat_capacity()
     *   - mass_fraction(n_s) with n_s in [0,solids_nb_species]
     *   - mass_fraction("species_name")
     *   - molar_concentration(n_g) with n_g in [0,solids_nb_species]
     *   - molar_concentration("species_name")
     *   - molar_mass(n_g) with n_g in [0,solids_nb_species]
     *   - molar_mass("species_name")
     *   - average_molar_mass()
     *   - mole_fraction(n_g) with n_g in [0,solids_nb_species]
     *   - mole_fraction("species_name")
     */
    template <class ReactionsData, class FluidData>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void operator() (ReactionsData& reactions,
                     const FluidData& fluid) const
    {
      amrex::ignore_unused(reactions);
      amrex::ignore_unused(fluid);
    }
};

#endif
