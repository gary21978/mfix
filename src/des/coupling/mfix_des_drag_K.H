#ifndef MFIX_DES_DRAG_K_H_
#define MFIX_DES_DRAG_K_H_

#include <mfix_des_txfr_coeff.H>

struct ComputeDragWenYu : public ComputeTxfrCoeff
{

  // Literature/Document References:                                     C
  // Wen CY, Yu YH (1966). Chemical Engineering Progress Symposium    C
  //    Series 62: 100-111.

  ComputeDragWenYu ( runtimeRealData a_runtimeReal_idxs,
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

    amrex::Real const mu_f = m_fluid_props.molViscosity(a_fluid_interp,
        m_fluid_idxs.T_g, m_fluid_idxs.X_gk);

    amrex::RealVect vel_f;
    vel_f[0] = a_fluid_interp[m_fluid_idxs.vel_g+0];
    vel_f[1] = a_fluid_interp[m_fluid_idxs.vel_g+1];
    vel_f[2] = a_fluid_interp[m_fluid_idxs.vel_g+2];

    amrex::RealVect vel_p;
    vel_p[0] = a_realarray[SoArealData::velx][a_id];
    vel_p[1] = a_realarray[SoArealData::vely][a_id];
    vel_p[2] = a_realarray[SoArealData::velz][a_id];

    amrex::RealVect vslp;
    vslp[0] = vel_f[0] - vel_p[0];
    vslp[1] = vel_f[1] - vel_p[1];
    vslp[2] = vel_f[2] - vel_p[2];

    amrex::Real const vrel = vslp.vectorLength();

    amrex::Real const diam_p = 2.0*a_realarray[SoArealData::radius][a_id];

    amrex::Real const Re = (mu_f > 0.0) ? diam_p*vrel*rop_f/mu_f : m_large_number;

    if (Re < m_epsilon) { return 0.0; }

    amrex::Real const C_d = (Re > 1000.0) ?  0.44 : Schiller_Naumann(Re);

    return 0.75 * C_d * vrel * rho_f * std::pow(eps_f, -2.65) / diam_p;
  }
};


struct ComputeDragGidaspow : public ComputeTxfrCoeff
{
  ComputeDragGidaspow ( runtimeRealData a_runtimeReal_idxs,
                        DualGridAuxIndexes& a_interp_idxs,
                        const MFIXFluidParms& a_fluid_parms,
                        const ThermoPropertyData& a_fluid_props,
                        const amrex::GeometryData& a_geomData )
    : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
        a_fluid_props, a_geomData )
  { }

  // Literature/Document References:
  //  Ding J, Gidaspow D (1990). AIChE Journal 36: 523-538.
  //  original source unknown:
  //  Lathouwers D, Bellan J (2000). Proceedings of the 2000 U.S. DOE
  //     Hydrogen Program Review NREL/CP-570-28890. Available from
  //  http://www.eere.energy.gov/hydrogenandfuelcells/pdfs/28890k.pdf.

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

    amrex::Real const mu_f = m_fluid_props.molViscosity(a_fluid_interp,
        m_fluid_idxs.T_g, m_fluid_idxs.X_gk);

    amrex::RealVect vel_f;
    vel_f[0] = a_fluid_interp[m_fluid_idxs.vel_g+0];
    vel_f[1] = a_fluid_interp[m_fluid_idxs.vel_g+1];
    vel_f[2] = a_fluid_interp[m_fluid_idxs.vel_g+2];

    amrex::RealVect vel_p;
    vel_p[0] = a_realarray[SoArealData::velx][a_id];
    vel_p[1] = a_realarray[SoArealData::vely][a_id];
    vel_p[2] = a_realarray[SoArealData::velz][a_id];

    amrex::RealVect vslp;
    vslp[0] = vel_f[0] - vel_p[0];
    vslp[1] = vel_f[1] - vel_p[1];
    vslp[2] = vel_f[2] - vel_p[2];

    amrex::Real const vrel = vslp.vectorLength();

    amrex::Real const diam_p = 2.0*a_realarray[SoArealData::radius][a_id];

    amrex::Real const Re = (mu_f > 0.0) ? diam_p*vrel*rop_f/mu_f : m_large_number;
    if (Re < m_epsilon) { return 0.0; }

    // Dense phase: eps_f <= 0.8
    amrex::Real const Ergun = 150.0*(1.0 - eps_f)*mu_f / (eps_f*diam_p*diam_p) + 1.75*rho_f*vrel/diam_p;

    // Dilute phase: eps_f > 0.8
    amrex::Real const C_d = (Re > 1000.0) ?  0.44 : Schiller_Naumann(Re);
    amrex::Real const WenYu = 0.75*C_d*vrel*rop_f*std::pow(eps_f, -2.65) / diam_p;

    // switch function
    amrex::Real const phi_fp = std::atan(150.0*1.75*(eps_f - 0.8))/M_PI + 0.5;

    // blend the models
    return (1.0 - phi_fp)*Ergun + phi_fp*WenYu;
  }
};


/**
   Literature/Document References:
   van der Hoef, Beetstra, Kuipers, J. Fluid. Mech., 528 (2005)

   Beetstra, van der Hoef, Kuipers, Chem. Eng. Science, 62
     (Jan 2007)

 */
struct ComputeDragBVK2 : public ComputeTxfrCoeff
{
    ComputeDragBVK2 ( runtimeRealData a_runtimeReal_idxs,
                      DualGridAuxIndexes& a_interp_idxs,
                      const MFIXFluidParms& a_fluid_parms,
                      const ThermoPropertyData& a_fluid_props,
                      const amrex::GeometryData& a_geomData )
    : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
        a_fluid_props, a_geomData )
  { }

  // Literature/Document References:
  //  Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62
  //  (Jan 2007)
  //
  //  Tang, Peters, Kuipers, Kriebitzsch, van der Hoef, AIChEJ,
  //  61(2) (Feb 2015)
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

    amrex::Real const mu_f = m_fluid_props.molViscosity(a_fluid_interp,
        m_fluid_idxs.T_g, m_fluid_idxs.X_gk);

    amrex::RealVect vel_f;
    vel_f[0] = a_fluid_interp[m_fluid_idxs.vel_g+0];
    vel_f[1] = a_fluid_interp[m_fluid_idxs.vel_g+1];
    vel_f[2] = a_fluid_interp[m_fluid_idxs.vel_g+2];

    amrex::RealVect vel_p;
    vel_p[0] = a_realarray[SoArealData::velx][a_id];
    vel_p[1] = a_realarray[SoArealData::vely][a_id];
    vel_p[2] = a_realarray[SoArealData::velz][a_id];

    amrex::RealVect vslp;
    vslp[0] = vel_f[0] - vel_p[0];
    vslp[1] = vel_f[1] - vel_p[1];
    vslp[2] = vel_f[2] - vel_p[2];

    amrex::Real const vrel = vslp.vectorLength();

    amrex::Real const diam_p = 2.0*a_realarray[SoArealData::radius][a_id];
    amrex::Real const eps_s = 1.0 - eps_f;

    amrex::Real const Re = (mu_f > 0.0) ? diam_p*vrel*rop_f/mu_f : m_large_number;
    if (Re < m_epsilon) { return 0.0; }

    amrex::Real const inv_pow_eps_f_4 = 1.0/(eps_f*eps_f*eps_f*eps_f);

    // (this F_Stokes is /= of Koch_Hill by a factor of eps_f)
    // see eq(5) J. Fluid. Mech. 528, 2005
    amrex::Real F_Stokes = 18.0*mu_f*eps_f/diam_p/diam_p;

    amrex::Real F(0.);

#if 0

    // Beetstra, van der Hoef, Kuipers (2007). Chem. Eng. Science, 62.

    F = 10.0*eps_s/eps_f/eps_f + eps_f*eps_f*(1.0 + 1.5*sqrt(eps_s));

    F += 0.413*Re/(24.0*eps_f*eps_f) *
        (1.0/eps_f + 3.0*eps_f*eps_s + 8.4/std::pow(Re, 0.343)) /
        (1.0 + std::pow(10.0, 3.0*eps_s)/std::pow(Re, 0.5 + 2.0*eps_s));

#elif 0

    // Tenneti, Garg, Subramaniam (2011). Int. J. Multiphase Flow, 37(9).

    F = (amrex::max(1. + 0.15*std::pow(Re, 0.687), 0.44*Re/24.)
        + 5.81*eps_s + 0.48*std::pow(eps_s, 1./3.)/eps_f)/eps_f/eps_f
        + eps_f*eps_s*eps_s*eps_s*Re*(0.95 + 0.61*eps_s*eps_s*eps_s/eps_f/eps_f);

#else

    // Tang, Peters, Kuipers, Kriebitzsch, van der Hoef (2015). AIChEJ, 61(2)

    F = 10.0*eps_s/eps_f/eps_f + eps_f*eps_f*(1.0 + 1.5*sqrt(eps_s));

    F += Re*(0.11*eps_s*(1.0+eps_s) - 4.56e-3*inv_pow_eps_f_4 +
         std::pow(Re, -0.343)*(0.169*eps_f + 6.44e-2*inv_pow_eps_f_4));

#endif
    return F*F_Stokes;
  }
};


/**********************************************************************************
*
*********************************************************************************/
struct ComputeDragSyamOBrien1988 : public ComputeTxfrCoeff
{
    amrex::Real const drag_c1;
    amrex::Real const drag_d1;

    ComputeDragSyamOBrien1988 ( runtimeRealData a_runtimeReal_idxs,
                                DualGridAuxIndexes& a_interp_idxs,
                                const MFIXFluidParms& a_fluid_parms,
                                const ThermoPropertyData& a_fluid_props,
                                const amrex::GeometryData& a_geomData,
                                amrex::Real a_c1, amrex::Real a_d1)
    : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
        a_fluid_props, a_geomData )
    , drag_c1(a_c1)
    , drag_d1(a_d1)
  { }

  // Literature/Document References:
  // Syamlal M, O'Brien TJ (1988). International Journal of Multiphase Flow 14: 473-481.
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

    amrex::Real const mu_f = m_fluid_props.molViscosity(a_fluid_interp,
        m_fluid_idxs.T_g, m_fluid_idxs.X_gk);

    amrex::RealVect vel_f;
    vel_f[0] = a_fluid_interp[m_fluid_idxs.vel_g+0];
    vel_f[1] = a_fluid_interp[m_fluid_idxs.vel_g+1];
    vel_f[2] = a_fluid_interp[m_fluid_idxs.vel_g+2];

    amrex::RealVect vel_p;
    vel_p[0] = a_realarray[SoArealData::velx][a_id];
    vel_p[1] = a_realarray[SoArealData::vely][a_id];
    vel_p[2] = a_realarray[SoArealData::velz][a_id];

    amrex::RealVect vslp;
    vslp[0] = vel_f[0] - vel_p[0];
    vslp[1] = vel_f[1] - vel_p[1];
    vslp[2] = vel_f[2] - vel_p[2];

    amrex::Real const vrel = vslp.vectorLength();

    amrex::Real const diam_p = 2.0*a_realarray[SoArealData::radius][a_id];

    amrex::Real const Re = (mu_f > 0.0) ? diam_p*vrel*rop_f/mu_f : m_large_number;
    if (Re < m_epsilon) { return 0.0; }

    // Calculate V_rm
    amrex::Real const A = std::pow(eps_f, 4.14);

    amrex::Real const B = (eps_f <= 0.85) ?
        drag_c1*std::pow(eps_f, 1.28) :
                std::pow(eps_f, drag_d1);

    amrex::Real const V_RM = 0.5*(A - 0.06*Re +
       std::sqrt(0.0036*Re*Re + 0.12*Re*(2.0*B - A) + A*A));

    // Calculate single sphere drag correlation of Dalla Valle (1948)
    // multiplied by the Reynolds number.
    amrex::Real C_dsXRe = (0.63*std::sqrt(Re/V_RM) + 4.8);
    C_dsXRe *= C_dsXRe;

    return 0.75*mu_f*eps_f*C_dsXRe/(V_RM*diam_p*diam_p);
  }
};


struct ComputeDragUser : public ComputeTxfrCoeff
{

  ComputeDragUser ( runtimeRealData a_runtimeReal_idxs,
                    DualGridAuxIndexes& a_interp_idxs,
                    const MFIXFluidParms& a_fluid_parms,
                    const ThermoPropertyData& a_fluid_props,
                    const amrex::GeometryData& a_geomData )
    : ComputeTxfrCoeff( a_runtimeReal_idxs, a_interp_idxs, a_fluid_parms,
        a_fluid_props, a_geomData )
  { }

  AMREX_GPU_HOST_DEVICE
#if defined(AMREX_USE_DPCPP)
  SYCL_EXTERNAL
#endif
  amrex::Real operator()
    ( int const a_id, MFIXParticleContainer::ParticleType& /*a_particle*/,
      const amrex::GpuArray<int*,                 SoAintData::count>& /*a_intarray*/,
      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& a_realarray,
      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*a_tile_data*/,
      amrex::Array4<amrex::Real const> const& /*a_fluid*/, amrex::Real const* const a_fluid_interp ) const;
};

#endif
