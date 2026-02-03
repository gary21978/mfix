#ifndef MFIX_DES_TXFR_COEFF_H_
#define MFIX_DES_TXFR_COEFF_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Geometry.H>

#include <mfix_pc.H>

#include <cmath>

struct ComputeTxfrCoeff
{

  /**
     *  \brief Returns: the calculated drag coefficient.
     *
     *  Inputs:
     *    eps_f     - gas volume fraction
     *    Mug     - gas laminar viscosity
     *    ROpg    - gas density * EP_g
     *    vrel    - magnitude of gas-solids relative velocity
     *    DPM     - particle diameter of solids phase M
     *    DPA     - average particle diameter
     *    PHIS    - solids volume fraction of solids phases
     *    fvelx   - x component of the fluid velocity at the particle position
     *    fvely   - y component of the fluid velocity at the particle position
     *    fvelz   - z component of the fluid velocity at the particle position
     *    i, j, k - particle cell indices
     *    pid     - particle id number
     */
  ComputeTxfrCoeff ( runtimeRealData a_runtimeReal_idxs,
                     DualGridAuxIndexes& a_interp_idxs,
                     const MFIXFluidParms& a_fluid_parms,
                     const ThermoPropertyData& a_fluid_props,
                     const amrex::GeometryData& a_geomData )
    : m_runtimeReal_idxs( a_runtimeReal_idxs )
    , m_fluid_idxs( a_interp_idxs )
    , m_fluid_params( a_fluid_parms)
    , m_fluid_props( a_fluid_props)
    , m_geomData( a_geomData)
    , contains_Tf (m_fluid_idxs.include_energy)
  { }

  protected:

    amrex::Real const m_epsilon = std::numeric_limits<amrex::Real>::epsilon();

    amrex::Real const m_small_number = 1.0e-15;
    amrex::Real const m_large_number = 1.0e+32;

    const runtimeRealData m_runtimeReal_idxs;
    const DualGridAuxIndexes m_fluid_idxs;

    const MFIXFluidParms m_fluid_params;
    const ThermoPropertyData m_fluid_props;
    const amrex::GeometryData m_geomData;

    int const contains_Tf;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Schiller_Naumann ( amrex::Real const a_Re ) const
    { return 24.0*(1.0 + 0.15*std::pow(a_Re, 0.687))/(a_Re+m_small_number); }

};

#endif
