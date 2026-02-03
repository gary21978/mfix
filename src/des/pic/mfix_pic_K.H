#include "mfix_pc.H"

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real solids_stress (const amrex::Real Ps,
                           const amrex::Real beta,
                           const amrex::Real ep_cp,
                           const amrex::Real small_number,
                           const amrex::Real eps)
{

  return Ps*std::pow(eps,beta) /
    amrex::max(ep_cp-eps, small_number*(1.-eps));
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::RealVect
updated_pic_velocity (const amrex::RealVect& pos,         // parcel position
                      const amrex::RealVect& vel_p_old,   // parcel velocity (current)
                      const amrex::Real dt,               // time step
                      const amrex::Real coeff,            // coeff of vel old in update
                      const amrex::Real scale,            // update denominator
                      const amrex::RealVect S_c,          // forces on particle
                      const amrex::RealVect grad_tau_p,   // particle normal stress grad
                      amrex::Array4<amrex::Real const> const& u_s,
                      amrex::Array4<amrex::Real const> const& v_s,
                      amrex::Array4<amrex::Real const> const& w_s,
                      const amrex::Real en,               // collision damping factor
                      const amrex::Real vol_p,            // particle volume
                      const amrex::Real eps_p,            // solids volume fraction
                      const amrex::Real ep_cp,            // close pack volume fraction
                      const amrex::Real vel_ref_frame,    // scale factor near packing
                      const amrex::Real mfp_vel_in,       // mfp limit on velocity
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& p_lo)
{

  // updated parcel velocity devoid of the particle normal stress
  RealVect vel = {scale*(coeff*vel_p_old[0] + dt*S_c[0]),
                  scale*(coeff*vel_p_old[1] + dt*S_c[1]),
                  scale*(coeff*vel_p_old[2] + dt*S_c[2])};

#if SCALE_MFP_DISABLED
  const Real mfp_vel = mfp_vel_in;
#elif SCALE_MFP_CHAPMAN
  const Real mfp_vel = mfp_vel_in / (1.0 + 2.5*eps_p);
#elif SCALE_MFP_CUTCHIS
  const Real mfp_vel = mfp_vel_in *
    amrex::max(0.0, ep_cp - eps_p) / (ep_cp + (2.5*ep_cp - 1.0)*eps_p);
#elif SCALE_MFP_MUSSER
  constexpr amrex::Real ep_cp_scale = 1.0/24.0;
  const Real eps_lim = ep_cp*(1.0 + ep_cp_scale);
  const Real mfp_vel = mfp_vel_in *
    amrex::max(0.0, eps_lim - eps_p) / (eps_lim + (2.5*eps_lim - 1.0)*eps_p);
#elif !MFP_DISABLED
  amrex::Abort("Invalid MFP set in mfix_pic_K.H");
#endif

  const Real lx = (pos[0] - p_lo[0])*dxi[0] + 0.5;
  const Real ly = (pos[1] - p_lo[1])*dxi[1] + 0.5;
  const Real lz = (pos[2] - p_lo[2])*dxi[2] + 0.5;

  const int i = static_cast<int>(Math::floor(lx));
  const int j = static_cast<int>(Math::floor(ly));
  const int k = static_cast<int>(Math::floor(lz));

  const Real wx_hi(lx - static_cast<Real>(i));
  const Real wy_hi(ly - static_cast<Real>(j));
  const Real wz_hi(lz - static_cast<Real>(k));

  const Real wx_lo(1 - wx_hi);
  const Real wy_lo(1 - wy_hi);
  const Real wz_lo(1 - wz_hi);

  const RealVect face_vel = {
    // u_s face velocity
    wy_lo*wz_lo*u_s(i  ,j-1,k-1,0) +
    wy_hi*wz_lo*u_s(i  ,j  ,k-1,0) +
    wy_lo*wz_hi*u_s(i  ,j-1,k  ,0) +
    wy_hi*wz_hi*u_s(i  ,j  ,k  ,0),
    // v_s face velocity
    wx_lo*wz_lo*v_s(i-1,j  ,k-1,0) +
    wx_hi*wz_lo*v_s(i  ,j  ,k-1,0) +
    wx_lo*wz_hi*v_s(i-1,j  ,k  ,0) +
    wx_hi*wz_hi*v_s(i  ,j  ,k  ,0),
    // w_s face velocity
    wx_lo*wy_lo*w_s(i-1,j-1,k  ,0) +
    wx_hi*wy_lo*w_s(i  ,j-1,k  ,0) +
    wx_lo*wy_hi*w_s(i-1,j  ,k  ,0) +
    wx_hi*wy_hi*w_s(i  ,j  ,k  ,0)};

  RealVect avg_vel(face_vel);

  const int ii = i + (vel[0]*grad_tau_p[0] <= 0.0 ? 0 :
     static_cast<int>(amrex::Math::copysign(1.0, vel[0])));

  const int jj = j + (vel[1]*grad_tau_p[1] <= 0.0 ? 0 :
     static_cast<int>(amrex::Math::copysign(1.0, vel[1])));

  const int kk = k + (vel[2]*grad_tau_p[2] <= 0.0 ? 0 :
     static_cast<int>(amrex::Math::copysign(1.0, vel[2])));

  if( ii != i ) {

    avg_vel[0] =
      wy_lo*wz_lo*u_s(ii ,j-1,k-1,0) +
      wy_hi*wz_lo*u_s(ii ,j  ,k-1,0) +
      wy_lo*wz_hi*u_s(ii ,j-1,k  ,0) +
      wy_hi*wz_hi*u_s(ii ,j  ,k  ,0);

    if( (vel[0] > 0.0 && wx_lo > 0.5) ||
        (vel[0] < 0.0 && wx_hi > 0.5)) {

      if( vel[0]*avg_vel[0] > 0.0 )
        avg_vel[0] = vel_ref_frame*face_vel[0];
    }
  }

  if( jj != j ) {

    avg_vel[1] =
      wx_lo*wz_lo*v_s(i-1,jj ,k-1,0) +
      wx_hi*wz_lo*v_s(i  ,jj ,k-1,0) +
      wx_lo*wz_hi*v_s(i-1,jj ,k  ,0) +
      wx_hi*wz_hi*v_s(i  ,jj ,k  ,0);

    if ((vel[1] > 0.0 && wy_lo > 0.5) ||
        (vel[1] < 0.0 && wy_hi > 0.5)) {

      if (vel[1]*avg_vel[1] > 0.0 )
        avg_vel[1] = vel_ref_frame*face_vel[1];
    }
  }

  if( kk != k ) {
    avg_vel[2] =
      wx_lo*wy_lo*w_s(i-1,j-1,kk ,0) +
      wx_hi*wy_lo*w_s(i  ,j-1,kk ,0) +
      wx_lo*wy_hi*w_s(i-1,j  ,kk ,0) +
      wx_hi*wy_hi*w_s(i  ,j  ,kk ,0);

    if ((vel[2] > 0.0 && wz_lo > 0.5) ||
        (vel[2] < 0.0 && wz_hi > 0.5)) {

      if( vel[2]*avg_vel[2] > 0.0 )
        avg_vel[2] = vel_ref_frame*face_vel[2];
    }
  }


  for (int dir(0); dir < 3; dir++) {

    // solids stress velocity contribution:
    // -dt*Vp*grad_tau_p / (ep_s * (mass + rho*Cvm*Vp + dt*beta/mass))
    const Real del_up = -((dt*vol_p)/eps_p)*grad_tau_p[dir] * scale;

    // Slip velocity between the parcel and the bulk. This needs to be the
    // tentative new velocity, not the current velocity.
    const Real slip_vel = avg_vel[dir] - vel[dir];

    // Add in contribution from stress gradient.
    // 1) Negative gradient. Push parcel in positive direction.
    // 2) Positive gradient. Push parcel in negative direction.
    // Both cases are limited limited by the slip velocity with the bulk.
    Real u_tau = (grad_tau_p[dir] <= 0.0 ) ?
      amrex::max(0.0, amrex::min(del_up, en*slip_vel)):
      amrex::min(0.0, amrex::max(del_up, en*slip_vel));

    vel[dir] += u_tau;

#if !MFP_DISABLED
    // Limit parcel motion when moving in the direction of normal stress.
    if(vel[dir]*del_up > 0.0) {
      if ( u_tau > std::numeric_limits<double>::min()) {
        // Parcel moving with positive stress
        vel[dir] = amrex::min( vel[dir], face_vel[dir] + mfp_vel);
        // Parcel moving with negative stress
      } else if ( u_tau < -std::numeric_limits<double>::min()) {
        vel[dir] = amrex::max( vel[dir], face_vel[dir] - mfp_vel);
      }
    }
#endif

  }

  return vel;
}
