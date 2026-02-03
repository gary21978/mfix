#ifndef _MFIX_DES_CONV_COEFF_K_H_
#define _MFIX_DES_CONV_COEFF_K_H_

#include <mfix_des_txfr_coeff.H>

struct NullConvectionCoeff : public ComputeTxfrCoeff
{
  // Constructor
  NullConvectionCoeff ( runtimeRealData a_runtimeReal_idxs,
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


struct ComputeConvRanzMarshall : public ComputeTxfrCoeff
{

  ComputeConvRanzMarshall ( runtimeRealData a_runtimeReal_idxs,
                            DualGridAuxIndexes& a_interp_idxs,
                            const MFIXFluidParms& a_fluid_parms,
                            const ThermoPropertyData& a_fluid_props,
                            const amrex::GeometryData& a_geomData )
    : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
        a_fluid_props, a_geomData )
  { }

  /*
   *  \brief Returns: the calculated heat transfer coefficient.
   *
   * Literature/Document References:
   * Ranz, W.E. and Marshall, W.R., "Friction and transfer coefficients for
   *     single particles and packed beds," Chemical Engineering Science,
   *     Vol. 48, No. 5, pp 247-253, 1925.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator()
    ( int const a_id, MFIXParticleContainer::ParticleType& /*a_particle*/,
      const amrex::GpuArray<int*,                 SoAintData::count>& /*a_intarray*/,
      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& a_realarray,
      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*a_tile_data*/,
      amrex::Array4<amrex::Real const> const& /*a_fluid*/, amrex::Real const* const a_fluid_interp ) const
  {

    amrex::Real const eps_f = a_fluid_interp[m_fluid_idxs.ep_g];
    amrex::Real const rho_f = a_fluid_interp[m_fluid_idxs.ro_g];
    amrex::Real const rop_f = rho_f * eps_f;

    amrex::Real const T_f  = a_fluid_interp[m_fluid_idxs.T_g];

    amrex::Real const mu_f = m_fluid_props.molViscosity(a_fluid_interp,
        m_fluid_idxs.T_g, m_fluid_idxs.X_gk);

    amrex::Real const k_f  = m_fluid_params.calc_k_g(T_f);

    amrex::Real cp_f = m_fluid_props.specificHeat(a_fluid_interp,
        m_fluid_idxs.T_g, m_fluid_idxs.X_gk);

    amrex::RealVect vel_f;
    vel_f[0] = a_fluid_interp[m_fluid_idxs.vel_g+0];
    vel_f[1] = a_fluid_interp[m_fluid_idxs.vel_g+1];
    vel_f[2] = a_fluid_interp[m_fluid_idxs.vel_g+2];

    amrex::RealVect vel_p;
    vel_p[0] = a_realarray[SoArealData::velx][a_id];
    vel_p[1] = a_realarray[SoArealData::vely][a_id];
    vel_p[2] = a_realarray[SoArealData::velz][a_id];

    amrex::Real const diam_p = 2.0*a_realarray[SoArealData::radius][a_id];

    amrex::RealVect vslp;
    vslp[0] = vel_f[0] - vel_p[0];
    vslp[1] = vel_f[1] - vel_p[1];
    vslp[2] = vel_f[2] - vel_p[2];

    amrex::Real const vrel = vslp.vectorLength();

    // Calculate the Reynolds Number
    amrex::Real N_Re = (mu_f > 0.0) ? diam_p*vrel*rop_f/mu_f : m_large_number;

    // Calculate the Prandtl Number
    amrex::Real N_Pr = (k_f  > 0.0) ? cp_f*mu_f/k_f : m_large_number;

    // Calculate the Nusselt Number
    amrex::Real N_Nu = 2.0 + 0.6 * std::sqrt(N_Re) * std::pow(N_Pr, 0.333);

    return N_Nu * k_f / diam_p;
  }
};


struct ComputeConvGunn : public ComputeTxfrCoeff
{

  ComputeConvGunn ( runtimeRealData a_runtimeReal_idxs,
                    DualGridAuxIndexes& a_interp_idxs,
                    const MFIXFluidParms& a_fluid_parms,
                    const ThermoPropertyData& a_fluid_props,
                    const amrex::GeometryData& a_geomData )
    : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
        a_fluid_props, a_geomData )
  { }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator()
    ( int const a_id, MFIXParticleContainer::ParticleType& /*a_particle*/,
      const amrex::GpuArray<int*,                 SoAintData::count>& /*a_intarray*/,
      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& a_realarray,
      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*a_tile_data*/,
      amrex::Array4<amrex::Real const> const& /*a_fluid*/, amrex::Real const* const a_fluid_interp ) const
  {

    amrex::Real const eps_f = a_fluid_interp[m_fluid_idxs.ep_g];
    amrex::Real const rho_f = a_fluid_interp[m_fluid_idxs.ro_g];
    amrex::Real const rop_f = rho_f * eps_f;

    amrex::Real const T_f  = a_fluid_interp[m_fluid_idxs.T_g];

    amrex::Real const mu_f = m_fluid_props.molViscosity(a_fluid_interp,
        m_fluid_idxs.T_g, m_fluid_idxs.X_gk);

    amrex::Real const k_f  = m_fluid_params.calc_k_g(T_f);

    amrex::Real cp_f = m_fluid_props.specificHeat(a_fluid_interp,
        m_fluid_idxs.T_g, m_fluid_idxs.X_gk);

    amrex::RealVect vel_f;
    vel_f[0] = a_fluid_interp[m_fluid_idxs.vel_g+0];
    vel_f[1] = a_fluid_interp[m_fluid_idxs.vel_g+1];
    vel_f[2] = a_fluid_interp[m_fluid_idxs.vel_g+2];

    amrex::RealVect vel_p;
    vel_p[0] = a_realarray[SoArealData::velx][a_id];
    vel_p[1] = a_realarray[SoArealData::vely][a_id];
    vel_p[2] = a_realarray[SoArealData::velz][a_id];

    amrex::Real const diam_p = 2.0*a_realarray[SoArealData::radius][a_id];

    amrex::RealVect vslp;
    vslp[0] = vel_f[0] - vel_p[0];
    vslp[1] = vel_f[1] - vel_p[1];
    vslp[2] = vel_f[2] - vel_p[2];

    amrex::Real const vrel = vslp.vectorLength();

    // Calculate the Reynolds Number
    amrex::Real N_Re = (mu_f > 0.0) ? diam_p*vrel*rop_f/mu_f : m_large_number;

    // Calculate the Prandtl Number
    amrex::Real N_Pr = (k_f  > 0.0) ? cp_f*mu_f/k_f : m_large_number;

    // Calculate the Nusselt Number
    amrex::Real N_Nu =
      (7 - 10*eps_f + 5*eps_f*eps_f)*(1 + .7*std::pow(N_Re, 0.2)*std::cbrt(N_Pr))
      + (1.33 - 2.4*eps_f + 1.2*eps_f*eps_f)*std::pow(N_Re, 0.7)*std::cbrt(N_Pr);

    return N_Nu * k_f / diam_p;
  }
};

#endif
