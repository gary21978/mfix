#ifndef _MFIX_EB_INTERP_K_H_
#define _MFIX_EB_INTERP_K_H_

#include "mfix_pc.H"

#include "mfix_eb_newton_solver_K.H"
#include "mfix_eb_node_connectivity_K.H"

/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void fe_interp(const amrex::RealVect& pos,
               const int ip,
               const int jp,
               const int kp,
               const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
               const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
               const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
               amrex::Array4<const EBCellFlag> const& flags,
               amrex::Array4<const amrex::Real> const& ccent,
               amrex::Array4<const amrex::Real> const& bcent,
               amrex::Array4<const amrex::Real> const& apx,
               amrex::Array4<const amrex::Real> const& apy,
               amrex::Array4<const amrex::Real> const& apz,
               amrex::Array4<amrex::Real const> const& interp_array,
               amrex::Real* interp_vals,
               const int ncomp,
               const int scomp)
{

    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();

    const amrex::Real cdist_x = pos[0] - (ip + 0.5 + ccent(ip,jp,kp,0))*dx[0] - plo[0];
    const amrex::Real cdist_y = pos[1] - (jp + 0.5 + ccent(ip,jp,kp,1))*dx[1] - plo[1];
    const amrex::Real cdist_z = pos[2] - (kp + 0.5 + ccent(ip,jp,kp,2))*dx[2] - plo[2];

    // Distance between particle and cell centoid.
    const amrex::Real cdist =
      std::sqrt(cdist_x*cdist_x + cdist_y*cdist_y + cdist_z*cdist_z);

    // Before doing anything fance, just check that the particle isn't overlapping
    // the cell center. If it is, use the cell value.
    if(cdist < tolerance){

        for(int n(0); n < scomp; n++)
            interp_vals[n] = interp_array(ip,jp,kp,n);

        for(int n=scomp; n<ncomp; n++)
            interp_vals[n] = interp_array(ip,jp,kp,n);

        return;
    }

    amrex::Real par_dot_EB(2.);
    amrex::Real cent_dot_EB(1.);

    // If the particle is inside a cut-cell, verify that it is on the correct
    // side of the EB before trying to interpolate.
    if(flags(ip,jp,kp).isSingleValued()){

        const Real adx = apx(ip+1,jp  ,kp  ) - apx(ip  ,jp  ,kp  );
        const Real ady = apy(ip  ,jp+1,kp  ) - apy(ip  ,jp  ,kp  );
        const Real adz = apz(ip  ,jp  ,kp+1) - apz(ip  ,jp  ,kp  );

        const Real apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);

        // Temporary sanity check.
        AMREX_ASSERT_WITH_MESSAGE( apnorm > std::numeric_limits<amrex::Real>::epsilon(),
            "Invalid EB normal. Not sure how to interpolate!");

        // EB normal
        const amrex::Real apnorminv = 1. / apnorm;
        const amrex::RealVect normal = {adx*apnorminv, ady*apnorminv, adz*apnorminv};

        // Projection of vector pointing from EB centroid to particle onto EB normal
        par_dot_EB = (pos[0] - (ip + 0.5 + bcent(ip,jp,kp,0))*dx[0] - plo[0])*normal[0] +
                     (pos[1] - (jp + 0.5 + bcent(ip,jp,kp,1))*dx[1] - plo[1])*normal[1] +
                     (pos[2] - (kp + 0.5 + bcent(ip,jp,kp,2))*dx[2] - plo[2])*normal[2];

        // Projection of vector pointing from EB centroid to cell centroid onto EB normal
        cent_dot_EB = (ccent(ip,jp,kp,0) - bcent(ip,jp,kp,0))*dx[0]*normal[0] +
                      (ccent(ip,jp,kp,1) - bcent(ip,jp,kp,1))*dx[1]*normal[1] +
                      (ccent(ip,jp,kp,2) - bcent(ip,jp,kp,2))*dx[2]*normal[2];

        // Temporary sanity check
        AMREX_ASSERT_WITH_MESSAGE( cent_dot_EB > tolerance,
            "cent_dot_EB < tolerance ... this makes no sense!");


    }
    // Particle position relative to cell center [-0.5, 0.5]
    amrex::Real gx = (pos[0]-plo[0])*dxi[0] - (ip + 0.5);
    amrex::Real gy = (pos[1]-plo[1])*dxi[1] - (jp + 0.5);
    amrex::Real gz = (pos[2]-plo[2])*dxi[2] - (kp + 0.5);

    // Use the centoid location of the cell containing the particle
    // to guess the interpolation stencil.
    const int i = (gx < ccent(ip,jp,kp,0)) ? ip : ip + 1;
    const int j = (gy < ccent(ip,jp,kp,1)) ? jp : jp + 1;
    const int k = (gz < ccent(ip,jp,kp,2)) ? kp : kp + 1;

    const int di = i - ip; // 0 or 1
    const int dj = j - jp; // 0 or 1
    const int dk = k - kp; // 0 or 1

    // Count the number of non-conntected cells in the stencil
    int covered = 0;
    for(int kk(-1); kk<1; kk++){
      for(int jj(-1); jj<1; jj++){
        for(int ii(-1); ii<1; ii++){
          if(!flags(ip,jp,kp).isConnected(di+ii,dj+jj,dk+kk))
            covered += 1;
        }
      }
    }

    // A negative value implies that the particle is 'behind' the EB (and
    // therefore inside the wall). Although this shouldn't occur often, it could
    // for fast moving particles.  Set the fluid velocity to zero and use the
    // cell value for remaining items.
    if( flags(ip,jp,kp).isSingleValued() && (par_dot_EB <= tolerance)){

        // Set the velocity to zero
        for (int n(0); n < scomp; n++)
            interp_vals[n] = 0.;

        for(int n=scomp; n<ncomp; n++){
            interp_vals[n] = interp_array(ip,jp,kp,n);
        }

    // The particle is near the EB. It is either
    // 1) between the cell centroid and the EB, or
    // 2) one or more cells in the 2^3 stencil is covered.
    // Either way, interpolating to the particle is not straight forward so
    // (for now) we do a 1D interpolation to the EB.

    } else if((covered > 0) || (par_dot_EB < cent_dot_EB)){

        amrex::Real scale(1.);

        for(int kk(-1); kk<1; kk++){
            for(int jj(-1); jj<1; jj++){
                for(int ii(-1); ii<1; ii++){
                    if(flags(i+ii,j+jj,k+kk).isSingleValued()) {

                        // EB normal
                        const amrex::Real adx = apx(i+ii+1,j+jj  ,k+kk  ) - apx(i+ii,j+jj,k+kk);
                        const amrex::Real ady = apy(i+ii  ,j+jj+1,k+kk  ) - apy(i+ii,j+jj,k+kk);
                        const amrex::Real adz = apz(i+ii  ,j+jj  ,k+kk+1) - apz(i+ii,j+jj,k+kk);
                        const amrex::Real apnorminv = 1. / std::sqrt(adx*adx + ady*ady + adz*adz);

                        const amrex::RealVect normal = {adx*apnorminv, ady*apnorminv, adz*apnorminv};

                        // Ratio of distance of particle-to-EB and cell-to-EB
                        amrex::Real test_scale = ((pos[0] - (i+ii + 0.5 + bcent(i+ii,j+jj,k+kk,0))*dx[0] - plo[0])*normal[0] +
                                                  (pos[1] - (j+jj + 0.5 + bcent(i+ii,j+jj,k+kk,1))*dx[1] - plo[1])*normal[1] +
                                                  (pos[2] - (k+kk + 0.5 + bcent(i+ii,j+jj,k+kk,2))*dx[2] - plo[2])*normal[2]) /
                                                 ((ccent(ip,jp,kp,0) - (di+ii + bcent(i+ii,j+jj,k+kk,0)))*dx[0]*normal[0] +
                                                  (ccent(ip,jp,kp,1) - (dj+jj + bcent(i+ii,j+jj,k+kk,1)))*dx[1]*normal[1] +
                                                  (ccent(ip,jp,kp,2) - (dk+kk + bcent(i+ii,j+jj,k+kk,2)))*dx[2]*normal[2]);

                        // Save the shortest positive distance.
                        scale = (test_scale > 0.) ? amrex::min(scale,test_scale) : scale;

                    }
                }
            }
        }

        // 1D interpolation from cell center to wall.

        for (int n(0); n < scomp; n++)
            interp_vals[n] = scale*interp_array(ip,jp,kp,n);

        for(int n=scomp; n<ncomp; n++)
            interp_vals[n] = interp_array(ip,jp,kp,n);

    } else {

        amrex::GpuArray<amrex::GpuArray<amrex::Real, 3>, 8> nodes;

        // Here we pack the local 2x2x2 stencil into a local array. Note that
        // the node ordering is not consistent with the previous version.
        // Nodes 2 and 3 have been swapped as well as nodes 6 and 7. This was
        // to allow for a more compact for-loop filling.

        int lc(0);

        for(int kk(-1); kk<1; kk++){
            for(int jj(-1); jj<1; jj++){
                for(int ii(-1); ii<1; ii++){
                    nodes[lc][0] = (i + ii + 0.5 + ccent(i+ii,j+jj,k+kk,0))*dx[0];
                    nodes[lc][1] = (j + jj + 0.5 + ccent(i+ii,j+jj,k+kk,1))*dx[1];
                    nodes[lc][2] = (k + kk + 0.5 + ccent(i+ii,j+jj,k+kk,2))*dx[2];
                    lc+=1;
                }
            }
        }

        amrex::Real xi   = (pos[0] - nodes[0][0] - plo[0])*dxi[0];
        amrex::Real eta  = (pos[1] - nodes[0][1] - plo[1])*dxi[1];
        amrex::Real zeta = (pos[2] - nodes[0][2] - plo[2])*dxi[2];

        const amrex::RealVect plo_vect(plo[0], plo[1], plo[2]);
        get_interp_mapping(pos-plo_vect, nodes, xi, eta, zeta);

        for(int n=0; n < ncomp; n++) {
            interp_vals[n] = (1 - xi - eta - zeta + xi*eta + xi*zeta + eta*zeta - xi*eta*zeta) * interp_array(i-1,j-1,k-1,n)
                           + (    xi              - xi*eta - xi*zeta            + xi*eta*zeta) * interp_array(i  ,j-1,k-1,n)
                           + (                      xi*eta                      - xi*eta*zeta) * interp_array(i  ,j  ,k-1,n)
                           + (         eta        - xi*eta           - eta*zeta + xi*eta*zeta) * interp_array(i-1,j  ,k-1,n)
                           + (               zeta          - xi*zeta - eta*zeta + xi*eta*zeta) * interp_array(i-1,j-1,k  ,n)
                           + (                               xi*zeta            - xi*eta*zeta) * interp_array(i  ,j-1,k  ,n)
                           + (                                                    xi*eta*zeta) * interp_array(i  ,j  ,k  ,n)
                           + (                                         eta*zeta - xi*eta*zeta) * interp_array(i-1,j  ,k  ,n);
        }
    }
}

#endif
