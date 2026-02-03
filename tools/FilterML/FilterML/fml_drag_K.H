#ifndef FML_DRAG_K_
#define FML_DRAG_K_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>

/**********************************************************************************
   Purpose: Calculate the gas-solids drag coefficient

   Literature/Document References:
   - van der Hoef, Beetstra, Kuipers, J. Fluid. Mech., 528 (2005)
   - Beetstra, van der Hoef, Kuipers, Chem. Eng. Science, 62 (Jan 2007)
   - Tenneti, Garg, Subramaniam, Int. J. Multiphase Flow, 37(9) (2011)
   - Tang, Peters, Kuipers, Kriebitzsch, van der Hoef, AIChEJ, 61(2) (Feb 2015)
**********************************************************************************/
struct compute_drag_Tang
{
  amrex::Real const m_small_number;
  amrex::Real const m_large_number;

  compute_drag_Tang ()
    : m_small_number(1.0e-15)
    , m_large_number(1.0e32)
  {}

  /**
   *  \brief Returns: the calculated drag coefficient.
   *
   *  Inputs:
   *    af     - fluid volume fraction (alpha_f)
   *    mu_f   - fluid laminar viscosity
   *    rho_f  - fluid density
   *    vrel   - magnitude of fluid-solids relative velocity
   *    dp     - average particle diameter
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (amrex::Real af, amrex::Real mu_f, amrex::Real rho_f,
                          amrex::Real relV, amrex::Real dp) const
  {
    amrex::Real const Re = (mu_f > 0.0) ? dp*relV*(af*rho_f)/mu_f : m_large_number;
    if (Re < m_small_number) return 0.0;

    // solids volume fraction (alpha_p)
    amrex::Real const ap(1.0 - af);

    amrex::Real inv_af4 = 1.0/(af*af*af*af);

    // F_Stokes is not equal to that in Koch_Hill by a factor of alpha_f
    // see eq(5) J. Fluid. Mech. 528, 2005
    amrex::Real F_Stokes = 18.0*mu_f*af/dp/dp;

    // BVK2:
    amrex::Real F = 10.0*ap/af/af + af*af*(1.0 + 1.5*sqrt(ap));

    F += Re*(0.11*ap*(1.0+ap) - 4.56e-3*inv_af4 +
         std::pow(Re, -0.343)*(0.169*af + 6.44e-2*inv_af4));

    return F*F_Stokes;
  }
};

#endif
