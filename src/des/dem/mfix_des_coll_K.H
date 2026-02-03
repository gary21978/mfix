#ifndef MFIX_DES_COLL_K_H_
#define MFIX_DES_COLL_K_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_RealVect.H>
#include <AMReX_Array.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cfrelvel (const amrex::RealVect& vel1,
               const amrex::RealVect& vel2,
               const amrex::Real& radius1,
               const amrex::Real& radius2,
               const amrex::RealVect& omega1,
               const amrex::RealVect& omega2,
               amrex::Real& vrn,
               amrex::RealVect& vslip,
               const amrex::RealVect& normal,
               const amrex::Real dist_mag)
{
    // translational relative velocity
    amrex::RealVect vreltrans(0.);
    vreltrans[0] =  vel1[0] - vel2[0];
    vreltrans[1] =  vel1[1] - vel2[1];
    vreltrans[2] =  vel1[2] - vel2[2];

    // distance from the contact point to the particle centers
    amrex::Real dist_c1 = (dist_mag*dist_mag
                          + radius1*radius1 - radius2*radius2) / (2.0*dist_mag);

    amrex::Real dist_c2 = dist_mag - dist_c1;

    amrex::RealVect omega_sum(0.);
    omega_sum[0] = omega1[0] * dist_c1 + omega2[0] * dist_c2;
    omega_sum[1] = omega1[1] * dist_c1 + omega2[1] * dist_c2;
    omega_sum[2] = omega1[2] * dist_c1 + omega2[2] * dist_c2;

    // the rotational relative velocity
    amrex::RealVect v_rot(0.);
    cross_product(omega_sum, normal, v_rot);

    // total relative velocity
    vreltrans[0] += v_rot[0];
    vreltrans[1] += v_rot[1];
    vreltrans[2] += v_rot[2];

    // normal component of relative velocity
    vrn = dot_product(vreltrans, normal);

    vslip[0] = vreltrans[0] - vrn*normal[0];
    vslip[1] = vreltrans[1] - vrn*normal[1];
    vslip[2] = vreltrans[2] - vrn*normal[2];
}

#endif
