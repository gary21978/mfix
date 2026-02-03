#ifndef _MFIX_UPWIND_H_
#define _MFIX_UPWIND_H_

namespace aux {

  //
  // Compute upwind non-normal velocity
  //
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real
  upwind (const amrex::Real velocity_minus, const amrex::Real velocity_plus, const amrex::Real u_edge)
  {
    // Small value to protect against tiny velocities used in upwinding
    const amrex::Real small_velocity(1.e-10);

    if(amrex::Math::abs(u_edge) < small_velocity)
      return .5*(velocity_minus+velocity_plus);

    return u_edge > 0 ? velocity_minus : velocity_plus;
  }
} // end namespace
#endif
