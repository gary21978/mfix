#ifndef MFIX_THERMO_PROPERTIES_H_
#define MFIX_THERMO_PROPERTIES_H_

#include <AMReX_REAL.H>

#include <mfix_prop_molecular_weight.H>
#include <mfix_prop_molecular_viscosity.H>
#include <mfix_prop_specific_heat.H>

struct ThermoPropertyData
{

  ThermoPropertyData ( const amrex::Real* a_MW, int const a_mV_model,
      int const a_mV_ncomp, const amrex::Real* a_mV_coeffs,
      int const a_nspecies, const amrex::Real* a_sH_coeffs,
      const int* a_sH_T_offset, const int* a_sH_coeff_offset, const amrex::Real* a_sH_T)
  : molViscosity( a_mV_model, a_mV_ncomp, a_MW, a_mV_coeffs)
  , enthalpy( a_nspecies, a_sH_coeffs, a_sH_T_offset, a_sH_coeff_offset, a_sH_T)
  , specificHeat( a_nspecies, a_sH_coeffs, a_sH_T_offset, a_sH_coeff_offset, a_sH_T)
  { }

  public:

    MolViscosityData molViscosity;
    EnthalpyData enthalpy;
    SpecificHeatData specificHeat;

};


class ThermoProperties
{

  public:

    template <amrex::RunOn run_on>
    [[nodiscard]] ThermoPropertyData data() const noexcept {
      if constexpr (run_on == amrex::RunOn::Cpu) {
        return { molWeight.m_h_MW.data(), molViscosity.m_model, molViscosity.m_ncomp, molViscosity.m_h_coeffs.data(),
                 molWeight.m_ncomp, specificHeat.m_h_coeffs.data(),
                 specificHeat.m_h_T_offset.data(), specificHeat.m_h_coeff_offset.data(), specificHeat.m_h_T.data() };
      } else {
        return { molWeight.m_d_MW.data(), molViscosity.m_model, molViscosity.m_ncomp, molViscosity.m_d_coeffs.data(),
                 molWeight.m_ncomp, specificHeat.m_d_coeffs.data(),
                 specificHeat.m_d_T_offset.data(), specificHeat.m_d_coeff_offset.data(), specificHeat.m_d_T.data() };
      }
    }

    /// Data containers

    MolecularWeight molWeight;
    MolViscosity molViscosity;
    SpecificHeat specificHeat;

};

#endif

