#ifndef MFIX_DES_ROLL_FRICTION_K_H_
#define MFIX_DES_ROLL_FRICTION_K_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>


struct NoRollingFriction
{
  public:
    NoRollingFriction ()
    {}

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::RealVect operator() (const amrex::Real& /*ri*/,
                                const amrex::Real& /*rj*/,
                                const amrex::Real& /*normal_force_mag*/,
                                const amrex::RealVect& /*omegai*/,
                                const amrex::RealVect& /*omegaj*/) const
    {
      return amrex::RealVect(0.);
    }
};


struct ModelA
{
  public:
    ModelA (const amrex::Real mu_r)
      : m_mu_r(mu_r)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::RealVect operator() (const amrex::Real& ri,
                                const amrex::Real& rj,
                                const amrex::Real& normal_force_mag,
                                const amrex::RealVect& omegai,
                                const amrex::RealVect& omegaj) const
    {
      amrex::Real Rij = (rj > 0.0) ? ri*rj/(ri+rj) : ri;
      amrex::RealVect omegaij = (rj > 0.0) ? omegai - omegaj : omegai;
      amrex::Real mag_omegaij = sqrt(dot_product(omegaij, omegaij));
      if (mag_omegaij > 1.0e-15)
      {
          return -1.0 * m_mu_r * Rij * normal_force_mag * omegaij / mag_omegaij;
      } else {
          return 0.0 * omegaij;
      }
    }

  private:
    const amrex::Real m_mu_r;
};



struct ModelB
{
  public:
    ModelB (const amrex::Real mu_r)
      : m_mu_r(mu_r)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::RealVect operator() (const amrex::Real& ri,
                                const amrex::Real& rj,
                                const amrex::Real& normal_force_mag,
                                const amrex::RealVect& omegai,
                                const amrex::RealVect& omegaj) const
    {
      if (rj > 0.0)
      {
          amrex::Real Rij = ri*rj/(ri+rj);
          amrex::RealVect omegaij = ri*omegai - rj*omegaj;
          return -1.0 * m_mu_r * Rij * normal_force_mag * omegaij;
      } else {
          return -1.0 * m_mu_r * ri * normal_force_mag * ri * omegai;
      }
    }

  private:
    const amrex::Real m_mu_r;
};


#endif
