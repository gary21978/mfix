#ifndef _MFIX_DES_ADD_MASS_COEFF_K_H_
#define _MFIX_DES_ADD_MASS_COEFF_K_H_

#include <mfix_des_txfr_coeff.H>

/*******************************************************************************************\
*  Structs that return the added mass force coefficient.                                    *
*                                                                                           *
*  Inputs:                                                                                  *
*    RO_par  - particle density                                                             *
*    RO_g    - gas density                                                                  *
*    EPg     - gas volume fraction                                                          *
*                                                                                           *
\********************************************************************************************/

struct AddedMassCoeff_Null : public ComputeTxfrCoeff
{
  AddedMassCoeff_Null ( runtimeRealData a_runtimeReal_idxs,
                        DualGridAuxIndexes& a_interp_idxs,
                        const MFIXFluidParms& a_fluid_parms,
                        const ThermoPropertyData& a_fluid_props,
                        const amrex::GeometryData& a_geomData )
  : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
      a_fluid_props, a_geomData )
  { }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator()
    ( int const /*a_id*/, MFIXParticleContainer::ParticleType& /*a_particle*/,
      const amrex::GpuArray<int*,                 SoAintData::count>& /*a_intarray*/,
      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& /*a_realarray*/,
      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*a_tile_data*/,
      amrex::Array4<amrex::Real const> const& /*a_fluid*/, amrex::Real const* const /*a_fluid_interp*/ ) const
  { return 0.0; }
};


struct AddedMassCoeff_Constant : public ComputeTxfrCoeff
{
  AddedMassCoeff_Constant ( runtimeRealData a_runtimeReal_idxs,
                            DualGridAuxIndexes& a_interp_idxs,
                            const MFIXFluidParms& a_fluid_parms,
                            const ThermoPropertyData& a_fluid_props,
                            const amrex::GeometryData& a_geomData,
                            amrex::Real const a_const)
  : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
      a_fluid_props, a_geomData )
  , m_const(a_const)
  { }

  amrex::Real const m_const;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator()
    ( int const /*a_id*/, MFIXParticleContainer::ParticleType& /*a_particle*/,
      const amrex::GpuArray<int*,                 SoAintData::count>& /*a_intarray*/,
      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& /*a_realarray*/,
      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*a_tile_data*/,
      amrex::Array4<amrex::Real const> const& /*a_fluid*/, amrex::Real const* const /*a_fluid_interp*/ ) const
  { return m_const; }
};


struct AddedMassCoeff_Zuber : public ComputeTxfrCoeff
{
  AddedMassCoeff_Zuber ( runtimeRealData a_runtimeReal_idxs,
                         DualGridAuxIndexes& a_interp_idxs,
                         const MFIXFluidParms& a_fluid_parms,
                         const ThermoPropertyData& a_fluid_props,
                         const amrex::GeometryData& a_geomData )
  : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
      a_fluid_props, a_geomData )
  { }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator()
    ( int const /*a_id*/, MFIXParticleContainer::ParticleType& /*a_particle*/,
      const amrex::GpuArray<int*,                 SoAintData::count>& /*a_intarray*/,
      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& /*a_realarray*/,
      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*a_tile_data*/,
      amrex::Array4<amrex::Real const> const& /*a_fluid*/, amrex::Real const* const a_fluid_interp ) const
  {
    amrex::Real const eps_f = a_fluid_interp[m_fluid_idxs.ep_g];

    return (1.5 - eps_f)/eps_f; // (1/2)*(1+2a)/(1-a)  whera a = 1-eps_f
  }
};


/*******************************************************************************************\
* Literature/Document References:                                                           *
*                                                                                           *
* T.M.J. Nijssen, H.A.M. Kuipers and J. van der Stel et al., Complete liquid-solid momentum *
* coupling for unresolved CFD-DEM simulations, International Journal of Multiphase Flow     *
*  132 (2020) 103425; Eqs. 12 -13; Table 2.                                                 *
*                                                                                           *
\********************************************************************************************/
struct AddedMassCoeff_Nijssen : public ComputeTxfrCoeff
{
  AddedMassCoeff_Nijssen ( runtimeRealData a_runtimeReal_idxs,
                           DualGridAuxIndexes& a_interp_idxs,
                           const MFIXFluidParms& a_fluid_parms,
                           const ThermoPropertyData& a_fluid_props,
                           const amrex::GeometryData& a_geomData )
  : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
      a_fluid_props, a_geomData )
  { }

  amrex::Real const Cvm0 =  0.5000;
  amrex::Real const a10  =  0.1300;
  amrex::Real const a11  =  0.0470;
  amrex::Real const a20  = -0.5800;
  amrex::Real const a21  = -0.0660;
  amrex::Real const a30  =  1.4200;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator()
    ( int const a_id, MFIXParticleContainer::ParticleType& /*a_particle*/,
      const amrex::GpuArray<int*,                 SoAintData::count>& /*a_intarray*/,
      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& a_realarray,
      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*a_tile_data*/,
      amrex::Array4<amrex::Real const> const& /*a_fluid*/, amrex::Real const* const a_fluid_interp ) const
  {
    amrex::Real const rho_f = a_fluid_interp[m_fluid_idxs.ro_g];

    amrex::Real const rho_p = 2.0*a_realarray[SoArealData::density][a_id];

    amrex::Real Cv1 = (a10 + a11*std::log(rho_p/rho_f)) * (1.0 - rho_f);
    amrex::Real Cv2 = (a20 + a21*std::log(rho_p/rho_f)) * (1.0 - rho_f) * (1.0 - rho_f);
    amrex::Real Cv3 = (a30                            ) * (1.0 - rho_f) * (1.0 - rho_f) * (1.0 - rho_f);

    return (Cvm0 + Cv1 + Cv2 + Cv3);
  }
};

#endif
