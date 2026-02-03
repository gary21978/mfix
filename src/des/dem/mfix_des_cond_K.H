#ifndef MFIX_DES_COND_K_H_
#define MFIX_DES_COND_K_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_RealVect.H>
#include <AMReX_Array.H>

// Particle-Particle conduction calculation
//   REF: Batchelor & O'Brien 1977
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real des_pp_conduction(amrex::Real dist_mag,
                              const amrex::Real p1radius,
                              const amrex::Real p2radius,
                              const amrex::Real kp1,
                              const amrex::Real kp2,
                              const amrex::Real Tp1,
                              const amrex::Real Tp2)
{
  amrex::Real kp_eff = 2.0*(kp1*kp2)/(kp1 + kp2);
  amrex::Real Olap   = p1radius + p2radius - dist_mag;
  amrex::Real Q_dot  = 0.0;

  /*
  // To-do: Area correction
  amrex::Real Hertz_act = 1.0;
  amrex::Real Hertz_sim = 1.0;
  amrex::Real Kn_sim    = 1.0;
  if (false) // do correction?
  {
    if (false) // is this Hertzian?
    {
      Olap = pow(Hertz_sim/Hertz_act,2.0/3.0)*Olap;
    }
    else // it's LSD
    {
      Olap = pow(Kn_sim*Olap/Hertz_act,2.0/3.0);
    }
  }
  */

  amrex::Real Rp_min = std::min(p1radius,p2radius);
  amrex::Real y      = Rp_min - Olap;
  amrex::Real lRad   = std::sqrt(Rp_min*Rp_min - y*y); // Contact radius -> base of triangle
  Q_dot = 2.0 * kp_eff * lRad * (Tp2 - Tp1); // Heat rate from 1 -> 2

  return Q_dot;
}

// Particle-Fluid-Particle conduction calculation
//   REF: Rong & Horio 1999
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real des_pfp_conduction(const amrex::Real dist_mag_eff,
                               const amrex::Real Rp_eff,
                               const amrex::Real Rlens_eff,
                               const amrex::Real Rough,
                               const amrex::Real kf,
                               const amrex::Real Tp1,
                               const amrex::Real Tp2)
{
  amrex::Real delta  = dist_mag_eff - Rp_eff;
  amrex::Real r1_tmp = std::sqrt(Rlens_eff*Rlens_eff - dist_mag_eff*dist_mag_eff);
  amrex::Real r0     = 0.0;
  amrex::Real Q_dot  = 0.0;

  // Limit lower bound to surface roughness
  if ( delta < Rough)
    r0 = sqrt(Rp_eff*Rp_eff - (Rough - dist_mag_eff)*(Rough - dist_mag_eff));
  // Limit upper bound to Rp_eff
  amrex::Real r1 = std::min(r1_tmp,Rp_eff);

  // Components of the solution
  amrex::Real C1 = sqrt(Rp_eff*Rp_eff - r1*r1);
  amrex::Real C0 = sqrt(Rp_eff*Rp_eff - r0*r0);
  amrex::Real A1 = amrex::Math::abs(dist_mag_eff-C1);
  amrex::Real A0 = amrex::Math::abs(dist_mag_eff-C0);
  Q_dot = M_PI * kf * ( dist_mag_eff*log(A1/A0)
      + C1 - C0 )*(Tp2 - Tp1); // Heat rate from 1 -> 2

  return Q_dot;
}

#endif
