#include <AMReX_Utility.H>

#include <mfix_des_drag_K.H>

AMREX_GPU_HOST_DEVICE
amrex::Real
ComputeDragUser::operator()
    ( int const a_id, MFIXParticleContainer::ParticleType& /*a_particle*/,
      const amrex::GpuArray<int*,                 SoAintData::count>& /*a_intarray*/,
      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& a_realarray,
      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& /*a_tile_data*/,
      amrex::Array4<amrex::Real const> const& /*a_fluid*/, amrex::Real const* const a_fluid_interp ) const
{
    amrex::Real const rho_f = a_fluid_interp[m_fluid_idxs.ro_g];

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

    amrex::Real const Re = (mu_f > 0.0) ? diam_p*vrel*rho_f/mu_f : m_large_number;

    amrex::Real Cd = 0.0;
    if (Re > m_epsilon) { Cd = (24.0/Re)*(1.0 + 0.15*std::pow(Re, 0.687)); }

    return 0.75*(rho_f*vrel/diam_p)*Cd;
}
