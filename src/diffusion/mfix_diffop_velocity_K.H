#ifndef MFIX_DIFFOP_VELOCITY_K_H_
#define MFIX_DIFFOP_VELOCITY_K_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Array.H>
#include <AMReX_EBCellFlag.H>

#include <mfix_compute_cell_grads.H>

//
// Eddy Viscosity
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_edd_visc_smagorinsky (int i, int j, int k,
                                       amrex::Array4<amrex::Real const> const& vel,
                                       amrex::Array4<amrex::Real const> const& rho,
                                       amrex::Array4<amrex::EBCellFlag const> const& flags,
                                       const amrex::RealVect& dx,
                                       amrex::Array4<amrex::Real const> const& vfrac,
                                       amrex::Real c)
{
  amrex::Real ux, uy, uz, vx, vy, vz, wx, wy, wz;
  mfix_comp_cell_grads(i, j, k, ux, uy, uz, vx, vy, vz,
        wx, wy, wz, vel, flags, dx);

  amrex::Real Sxx = ux;
  amrex::Real Sxy = 0.5 * (uy + vx);
  amrex::Real Sxz = 0.5 * (uz + wx);
  amrex::Real Syx = 0.5 * (vx + uy);
  amrex::Real Syy = vy;
  amrex::Real Syz = 0.5 * (vz + wy);
  amrex::Real Szx = 0.5 * (uz + wx);
  amrex::Real Szy = 0.5 * (wy + vz);
  amrex::Real Szz = wz;

  amrex::Real trS2 = Sxx*Sxx + Syy*Syy + Szz*Szz + 2.0*(Sxy*Syx + Sxz*Szx + Syz*Szy);
  amrex::Real Ig2D = 0.5 * trS2;

  amrex::Real delta = std::cbrt(dx[0]*dx[1]*dx[2]*vfrac(i,j,k));

  return (2.0 * rho(i,j,k) * c*c *delta*delta * std::sqrt(Ig2D));
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_edd_visc_wale (int i, int j, int k,
                                amrex::Array4<amrex::Real const> const& vel,
                                amrex::Array4<amrex::Real const> const& rho,
                                amrex::Array4<amrex::EBCellFlag const> const& flags,
                                const amrex::RealVect& dx,
                                amrex::Array4<amrex::Real const> const& vfrac,
                                amrex::Real c)
{
  amrex::Real ux, uy, uz, vx, vy, vz, wx, wy, wz;
  mfix_comp_cell_grads(i, j, k, ux, uy, uz, vx, vy, vz,
        wx, wy, wz, vel, flags, dx);

  amrex::Real Sxx = ux;
  amrex::Real Sxy = 0.5 * (uy + vx);
  amrex::Real Sxz = 0.5 * (uz + wx);
  amrex::Real Syx = 0.5 * (vx + uy);
  amrex::Real Syy = vy;
  amrex::Real Syz = 0.5 * (vz + wy);
  amrex::Real Szx = 0.5 * (uz + wx);
  amrex::Real Szy = 0.5 * (wy + vz);
  amrex::Real Szz = wz;

  amrex::Real Oxx = 0.;
  amrex::Real Oxy = 0.5 * (uy - vx);
  amrex::Real Oxz = 0.5 * (uz - wx);
  amrex::Real Oyx = 0.5 * (vx - uy);
  amrex::Real Oyy = 0.;
  amrex::Real Oyz = 0.5 * (vz - wy);
  amrex::Real Ozx = 0.5 * (wx - uz);
  amrex::Real Ozy = 0.5 * (wy - vz);
  amrex::Real Ozz = 0.;

  amrex::Real trS2 = Sxx*Sxx + Syy*Syy + Szz*Szz + 2.0*(Sxy*Syx + Sxz*Szx + Syz*Szy);
  amrex::Real trOOt = Oxy*Oxy + Oxz*Oxz + Oyx*Oyx + Oyz*Oyz + Ozx*Ozx + Ozy*Ozy;

  amrex::Real S2xx = Sxx*Sxx + Sxy*Syx + Sxz*Szx;
  amrex::Real S2xy = Sxx*Sxy + Sxy*Syy + Sxz*Szy;
  amrex::Real S2xz = Sxx*Sxz + Sxy*Syz + Sxz*Szz;
  amrex::Real S2yx = Syx*Sxx + Syy*Syx + Syz*Szx;
  amrex::Real S2yy = Syx*Sxy + Syy*Syy + Syz*Szy;
  amrex::Real S2yz = Syx*Sxz + Syy*Syz + Syz*Szz;
  amrex::Real S2zx = Szx*Sxx + Szy*Syx + Szz*Szx;
  amrex::Real S2zy = Szx*Sxy + Szy*Syy + Szz*Szy;
  amrex::Real S2zz = Szx*Sxz + Szy*Syz + Szz*Szz;

  amrex::Real O2xx = Oxx*Oxx + Oxy*Oyx + Oxz*Ozx;
  amrex::Real O2xy = Oxx*Oxy + Oxy*Oyy + Oxz*Ozy;
  amrex::Real O2xz = Oxx*Oxz + Oxy*Oyz + Oxz*Ozz;
  amrex::Real O2yx = Oyx*Oxx + Oyy*Oyx + Oyz*Ozx;
  amrex::Real O2yy = Oyx*Oxy + Oyy*Oyy + Oyz*Ozy;
  amrex::Real O2yz = Oyx*Oxz + Oyy*Oyz + Oyz*Ozz;
  amrex::Real O2zx = Ozx*Oxx + Ozy*Oyx + Ozz*Ozx;
  amrex::Real O2zy = Ozx*Oxy + Ozy*Oyy + Ozz*Ozy;
  amrex::Real O2zz = Ozx*Oxz + Ozy*Oyz + Ozz*Ozz;

  amrex::Real S2O2xx = S2xx*O2xx + S2xy*O2yx + S2xz*O2zx;
  //amrex::Real S2O2xy = S2xx*O2xy + S2xy*O2yy + S2xz*O2zy;
  //amrex::Real S2O2xz = S2xx*O2xz + S2xy*O2yz + S2xz*O2zz;
  //amrex::Real S2O2yx = S2yx*O2xx + S2yy*O2yx + S2yz*O2zx;
  amrex::Real S2O2yy = S2yx*O2xy + S2yy*O2yy + S2yz*O2zy;
  //amrex::Real S2O2yz = S2yx*O2xz + S2yy*O2yz + S2yz*O2zz;
  //amrex::Real S2O2zx = S2zx*O2xx + S2zy*O2yx + S2zz*O2zx;
  //amrex::Real S2O2zy = S2zx*O2xy + S2zy*O2yy + S2zz*O2zy;
  amrex::Real S2O2zz = S2zx*O2xz + S2zy*O2yz + S2zz*O2zz;

  amrex::Real trS2O2 = S2O2xx + S2O2yy + S2O2zz;

  amrex::Real Sd2 = (1./6.)*(trS2*trS2 + trOOt*trOOt) + (2./3.)*(trS2*trOOt) + 2.*trS2O2;

  amrex::Real delta = std::cbrt(dx[0]*dx[1]*dx[2]*vfrac(i,j,k));

  if (Sd2 <= 0. || trS2  <= 0.) {
     return 0.;
  }

  return (rho(i,j,k) * std::pow(delta*c, 2) * std::pow(Sd2,(3./2.))) /
    (std::pow(trS2,(5./2.)) + std::pow(Sd2,(5./4.)));
}

//
// Suspension Viscosity
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_susp_visc_einstein (int i, int j, int k,
                                     amrex::Array4<amrex::Real const> const& epf,
                                     amrex::Real cell_mol_visc)
{
  amrex::Real eps = 1. - epf(i,j,k);
  return (2.5 * eps * cell_mol_visc);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_susp_visc_brinkman (int i, int j, int k,
                                     amrex::Array4<amrex::Real const> const& epf,
                                     amrex::Real cell_mol_visc,
                                     amrex::Real c)
{
  amrex::Real mu_star = std::pow(epf(i,j,k), -c);

  return (mu_star - 1.)* cell_mol_visc;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_susp_visc_roscoe (int i, int j, int k,
                                   amrex::Array4<amrex::Real const> const& epf,
                                   amrex::Real cell_mol_visc,
                                   amrex::Real c1,
                                   amrex::Real c2)
{
  amrex::Real mu_star = std::pow(std::max(c1 - (1. - epf(i,j,k)), 1.0e-9) / c1, -c2);

  return (mu_star - 1.)* cell_mol_visc;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_susp_visc_chenglaw (int i, int j, int k,
                                     amrex::Array4<amrex::Real const> const& epf,
                                     amrex::Real cell_mol_visc,
                                     amrex::Real c)
{
  amrex::Real mu_star = std::exp(2.5*(std::pow(epf(i,j,k), -c) - 1.)/c);

  return (mu_star - 1.)* cell_mol_visc;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_susp_visc_sato (int i, int j, int k,
                                 amrex::Array4<amrex::Real const> const& epf,
                                 amrex::Array4<amrex::Real const> const& vel,
                                 amrex::Array4<amrex::Real const> const& rho,
                                 amrex::Array4<amrex::Real const> const& avg_radius,
                                 amrex::Array4<amrex::Real const> const& avg_velp_x,
                                 amrex::Array4<amrex::Real const> const& avg_velp_y,
                                 amrex::Array4<amrex::Real const> const& avg_velp_z,
                                 amrex::Real c)
{
  amrex::Real urel = std::sqrt(  std::pow(vel(i,j,k,0) - avg_velp_x(i,j,k), 2)
                               + std::pow(vel(i,j,k,1) - avg_velp_y(i,j,k), 2)
                               + std::pow(vel(i,j,k,2) - avg_velp_z(i,j,k), 2)
                              );

  return c * ( 1. - epf(i,j,k)) * 2 * avg_radius(i,j,k) * urel * rho(i,j,k);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real calc_susp_visc_subramaniam (int i, int j, int k,
                                        amrex::Array4<amrex::Real const> const& epf,
                                        amrex::Real cell_mol_visc,
                                        amrex::Array4<amrex::Real const> const& vel,
                                        amrex::Array4<amrex::Real const> const& rho,
                                        amrex::Array4<amrex::Real const> const& avg_radius,
                                        amrex::Array4<amrex::Real const> const& avg_velp_x,
                                        amrex::Array4<amrex::Real const> const& avg_velp_y,
                                        amrex::Array4<amrex::Real const> const& avg_velp_z,
                                        amrex::Real c)
{
  // if phi < 1.0e-8, neglect PIT
  if (epf(i,j,k) + 1.0e-8 > 1.0) { return 0.0; }

  // small number
  amrex::Real const m_epsilon = std::numeric_limits<amrex::Real>::epsilon();

  // relative velocity mag
  amrex::Real urel = std::sqrt(  std::pow(vel(i,j,k,0) - avg_velp_x(i,j,k), 2)
                               + std::pow(vel(i,j,k,1) - avg_velp_y(i,j,k), 2)
                               + std::pow(vel(i,j,k,2) - avg_velp_z(i,j,k), 2)
                              );

  // mean flow Reynolds number
  amrex::Real const Rem = rho(i,j,k)*epf(i,j,k)*urel*2.0*avg_radius(i,j,k)/std::max(cell_mol_visc, m_epsilon);

  // if Rem is small, neglect PIT
  if (Rem < m_epsilon) { return 0.0; }

  // pseudo-turbulent kinetic energy in the fluid, see (7.1) of Mehrabadi etal (2015)
  amrex::Real expval = (epf(i,j,k) - 1.0)*std::sqrt(Rem);
  amrex::Real ptke_f = 2.0*(1.0 - epf(i,j,k))
                     + 2.5*(1.0 - epf(i,j,k))*std::pow(epf(i,j,k),3.0)*std::exp(expval);
  ptke_f *= 0.5*urel*urel;

  // drag force / Stokes drag, see (10.3) of Mehrabadi etal (2015)
  // replaced TGS drag with Tang "BVK2", see wdf notes
  amrex::Real const inv_epf_pow_4 = 1.0/(epf(i,j,k)*epf(i,j,k)*epf(i,j,k)*epf(i,j,k));
  amrex::Real F = 10.0*(1.0 - epf(i,j,k))/epf(i,j,k)/epf(i,j,k)
                + epf(i,j,k)*epf(i,j,k)*(1.0 + 1.5*sqrt(1.0 - epf(i,j,k)))
                + Rem*(0.11*(1.0 - epf(i,j,k))*(2.0 - epf(i,j,k)) - 4.56e-3*inv_epf_pow_4
                       + (0.169*epf(i,j,k) + 6.44e-2*inv_epf_pow_4)*std::pow(Rem, -0.343));

  // mean dissipation rate, see Eq. (10.1) of Mehrabadi etal (2015)
  amrex::Real ptEps_m = 4.5*(1.0 - epf(i,j,k))*cell_mol_visc*F*urel*urel
                        /rho(i,j,k)/avg_radius(i,j,k)/avg_radius(i,j,k);

  // rho_f * epg * epsilon_f = Epsilon_m, see Eq. (10.6) of Mehrabadi etal (2015)

  // nu_PIT = C_mu * k^2 / epsilon, see pg 221 of Tenneti & Subramaniam (2014)
  return rho(i,j,k)*c*ptke_f*ptke_f/ptEps_m;
}

#endif
