#ifndef MFIX_USR_REACTIONS_RATES_K_H
#define MFIX_USR_REACTIONS_RATES_K_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>

#include <cmath>


class LagrangianReactionRates
{
  public:
    AMREX_GPU_HOST_DEVICE
    LagrangianReactionRates () = default;

    /**
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
    template <class Reactions, class Fluid, class Solids>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (Reactions& reactions,
                     const Fluid& fluid,
                     const Solids& solids) const
    {

#if defined(CHEM_TEST_LAGRANGIAN01) || \
    defined(CHEM_TEST_LAGRANGIAN02) || \
    defined(CHEM_TEST_LAGRANGIAN03)

      constexpr int A = 0;

      const amrex::Real time = reactions.time();
      const amrex::Real dt = reactions.dt();
      const amrex::Real t = time + dt/2;

      const amrex::Real y_dot = 0.15*M_PI*(     std::sin(6.*M_PI*t)*std::sin(3.*M_PI*t)
                                          - 2.0*std::cos(6.*M_PI*t)*std::cos(3.*M_PI*t));

      const amrex::Real vfrac_f = fluid.volume_fraction();
      const amrex::Real MW_A = fluid.molar_mass(A);
      const amrex::Real vol_c = fluid.volume_fraction()*reactions.cell_volume();
      const amrex::Real Np = solids.particles_nb_in_cell();

      const amrex::Real rate = y_dot * (vol_c * vfrac_f / MW_A) / Np;

      reactions.rate(0) = rate;

#elif defined(CHEM_TEST_LAGRANGIAN99)

      constexpr int A = 0;
      constexpr int AtoC = 0;

      const amrex::Real time = reactions.time();
      const amrex::Real dt = reactions.dt();
      const amrex::Real t = time + dt/2;

      const amrex::Real y_dot = 1.5*M_PI*std::cos(6.*M_PI*t) /
                             (0.5 + 0.25*std::sin(6.*M_PI*t));

      const amrex::Real vfrac_f = fluid.volume_fraction();
      const amrex::Real rho_f = fluid.density();
      const amrex::Real MW_A = fluid.molar_mass(A);
      const amrex::Real vol_c = fluid.volume_fraction()*reactions.cell_volume();
      const amrex::Real Np = solids.particles_nb_in_cell();

      reactions.rate(AtoC) = y_dot * (vol_c * vfrac_f * rho_f / MW_A) / Np;

#else

      amrex::ignore_unused(reactions);
      amrex::ignore_unused(fluid);
      amrex::ignore_unused(solids);

#endif
    }
};


class EulerianReactionRates
{
  public:
    AMREX_GPU_HOST_DEVICE
    EulerianReactionRates () = default;

    /**
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
    template <class Reactions, class Fluid>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (Reactions& reactions,
                     const Fluid& fluid) const
    {

#if defined(CHEM_TEST_EULERIAN01)

      constexpr int A = 0;

      constexpr int AtoB = 0;

      const amrex::Real time = reactions.time();
      const amrex::Real dt = reactions.dt();
      const amrex::Real t = time + dt/2;

      const amrex::Real y_dot = 1.5*M_PI*std::cos(6.*M_PI*t);

      const amrex::Real vfrac_f = fluid.volume_fraction();
      const amrex::Real rho_f = fluid.density();
      const amrex::Real MW_A = fluid.molar_mass(A);

      reactions.rate(AtoB) = y_dot * (vfrac_f * rho_f / MW_A);

#elif defined(CHEM_TEST_EULERIAN02)

      // Adiabatic flame:
      // CH4_Comb:  CH4 + 2O2 --> CO2 + 2H2O  (mol/cm^3.s)
      //---------------------------------------------------------------------//
      // Note: The CH4 combustion rate is artificial and used for the
      // adiabatic flame test case.

      amrex::Real const c_O2  = fluid.molar_concentration("O2");   // [mol/m^3]
      amrex::Real const c_CH4 = fluid.molar_concentration("CH4");  // [mol/m^3]

      amrex::Real k0(1.0e3); // reaction rate coefficient

      reactions.rate(0) = k0 * c_O2 * c_CH4;

#else

      amrex::ignore_unused(reactions);
      amrex::ignore_unused(fluid);

#endif
    }
};

#endif
