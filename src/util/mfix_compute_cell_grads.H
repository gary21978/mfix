#ifndef _MFIX_COMPUTE_CELL_GRADS_H_
#define _MFIX_COMPUTE_CELL_GRADS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Array.H>
#include <AMReX_EBCellFlag.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mfix_comp_cell_grads (int i, int j, int k,
      amrex::Real& ux, amrex::Real& uy, amrex::Real& uz,
      amrex::Real& vx, amrex::Real& vy, amrex::Real& vz,
      amrex::Real& wx, amrex::Real& wy, amrex::Real& wz,
      amrex::Array4<amrex::Real const> const& vel,
      amrex::Array4<amrex::EBCellFlag const> const& flag,
      const amrex::RealVect& dx)
{
   // If the cell is covered, set all gradients to zero
   if (flag(i,j,k).isCovered()) {
      ux = uy = uz = vx = vy = vz = wx = wy = wz = 0.;
      return;
   }

   constexpr amrex::Real c0 = -1.5;
   constexpr amrex::Real c1 =  2.0;
   constexpr amrex::Real c2 = -0.5;
   const amrex::Real odx = 1./dx[0];
   const amrex::Real ody = 1./dx[1];
   const amrex::Real odz = 1./dx[2];

   if (!flag(i,j,k).isConnected( 1,0,0))
   {
      // Covered cell to the right, go fish left

      if (!flag(i,j,k).isConnected(-1,0,0)) {
         // Covered cell to the left and right, set to zero
         ux = vx = wx = 0.;

      } else if (!flag(i-1,j,k).isConnected(-1,0,0)) {
         // Covered cell two spaces to left, need to settle for
         // first order approx
         ux = (vel(i,j,k,0) - vel(i-1,j,k,0)) * odx;
         vx = (vel(i,j,k,1) - vel(i-1,j,k,1)) * odx;
         wx = (vel(i,j,k,2) - vel(i-1,j,k,2)) * odx;

      } else {
         ux = - (c0 * vel(i  ,j,k,0)
               + c1 * vel(i-1,j,k,0)
               + c2 * vel(i-2,j,k,0)) * odx;
         vx = - (c0 * vel(i  ,j,k,1)
               + c1 * vel(i-1,j,k,1)
               + c2 * vel(i-2,j,k,1)) * odx;
         wx = - (c0 * vel(i  ,j,k,2)
               + c1 * vel(i-1,j,k,2)
               + c2 * vel(i-2,j,k,2)) * odx;
      }
   }
   else if (!flag(i,j,k).isConnected(-1,0,0))
   {
      // Covered cell to the left, go fish right

      if (!flag(i,j,k).isConnected(1,0,0)) {
         // Covered cell to the right and left, set to zero
         ux = vx = wx = 0.;

      } else if (!flag(i+1,j,k).isConnected(1,0,0)) {
         // Covered cell two spaces to the right, need to settle for
         // first order approx
         ux = (vel(i+1,j,k,0) - vel(i,j,k,0)) * odx;
         vx = (vel(i+1,j,k,1) - vel(i,j,k,1)) * odx;
         wx = (vel(i+1,j,k,2) - vel(i,j,k,2)) * odx;

      } else {
         ux = (c0 * vel(i  ,j,k,0)
             + c1 * vel(i+1,j,k,0)
             + c2 * vel(i+2,j,k,0)) * odx;
         vx = (c0 * vel(i  ,j,k,1)
             + c1 * vel(i+1,j,k,1)
             + c2 * vel(i+2,j,k,1)) * odx;
         wx = (c0 * vel(i  ,j,k,2)
             + c1 * vel(i+1,j,k,2)
             + c2 * vel(i+2,j,k,2)) * odx;
      }
   }
   else
   {
      // No covered cells right or left, use standard stencil
      ux = 0.5 * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * odx;
      vx = 0.5 * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * odx;
      wx = 0.5 * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * odx;
   }

   // Do the same in y-direction
   if (!flag(i,j,k).isConnected(0, 1,0))
   {
      if (!flag(i,j,k).isConnected(0,-1,0)) {
         uy = vy = wy = 0.;

      } else if (!flag(i,j-1,k).isConnected(0,-1,0)) {
         uy = (vel(i,j,k,0) - vel(i,j-1,k,0)) * ody;
         vy = (vel(i,j,k,1) - vel(i,j-1,k,1)) * ody;
         wy = (vel(i,j,k,2) - vel(i,j-1,k,2)) * ody;

      } else {
         uy = - (c0 * vel(i,j  ,k,0)
               + c1 * vel(i,j-1,k,0)
               + c2 * vel(i,j-2,k,0)) * ody;
         vy = - (c0 * vel(i,j  ,k,1)
               + c1 * vel(i,j-1,k,1)
               + c2 * vel(i,j-2,k,1)) * ody;
         wy = - (c0 * vel(i,j  ,k,2)
               + c1 * vel(i,j-1,k,2)
               + c2 * vel(i,j-2,k,2)) * ody;
      }
   }
   else if (!flag(i,j,k).isConnected(0,-1,0))
   {
      if (!flag(i,j,k).isConnected(0,1,0)) {
         uy = vy = wy = 0.;

      } else if (!flag(i,j+1,k).isConnected(0,1,0)) {
         uy = (vel(i,j+1,k,0) - vel(i,j,k,0)) * ody;
         vy = (vel(i,j+1,k,1) - vel(i,j,k,1)) * ody;
         wy = (vel(i,j+1,k,2) - vel(i,j,k,2)) * ody;

      } else {
         uy = (c0 * vel(i,j  ,k,0)
             + c1 * vel(i,j+1,k,0)
             + c2 * vel(i,j+2,k,0)) * ody;
         vy = (c0 * vel(i,j  ,k,1)
             + c1 * vel(i,j+1,k,1)
             + c2 * vel(i,j+2,k,1)) * ody;
         wy = (c0 * vel(i,j  ,k,2)
             + c1 * vel(i,j+1,k,2)
             + c2 * vel(i,j+2,k,2)) * ody;
      }
   }
   else
   {
      uy = 0.5 * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * ody;
      vy = 0.5 * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * ody;
      wy = 0.5 * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * ody;
   }
   // Do the same in z-direction
   if (!flag(i,j,k).isConnected(0,0, 1))
   {
      if (!flag(i,j,k).isConnected(0,0,-1)) {
         uz = vz = wz = 0.;

      } else if (!flag(i,j,k-1).isConnected(0,0,-1)) {
         uz = (vel(i,j,k,0) - vel(i,j,k-1,0)) * odz;
         vz = (vel(i,j,k,1) - vel(i,j,k-1,1)) * odz;
         wz = (vel(i,j,k,2) - vel(i,j,k-1,2)) * odz;

      } else {
         uz = - (c0 * vel(i,j,k  ,0)
               + c1 * vel(i,j,k-1,0)
               + c2 * vel(i,j,k-2,0)) * odz;
         vz = - (c0 * vel(i,j,k  ,1)
               + c1 * vel(i,j,k-1,1)
               + c2 * vel(i,j,k-2,1)) * odz;
         wz = - (c0 * vel(i,j,k  ,2)
               + c1 * vel(i,j,k-1,2)
               + c2 * vel(i,j,k-2,2)) * odz;
      }
   }
   else if (!flag(i,j,k).isConnected(0,0,-1))
   {
      if (!flag(i,j,k).isConnected(0,0,1)) {
         uz = vz = wz = 0.;

      } else if (!flag(i,j,k+1).isConnected(0,0,1)) {
         uz = (vel(i,j,k+1,0) - vel(i,j,k,0)) * odz;
         vz = (vel(i,j,k+1,1) - vel(i,j,k,1)) * odz;
         wz = (vel(i,j,k+1,2) - vel(i,j,k,2)) * odz;

      } else {
         uz = (c0 * vel(i,j,k  ,0)
             + c1 * vel(i,j,k+1,0)
             + c2 * vel(i,j,k+2,0)) * odz;
         vz = (c0 * vel(i,j,k  ,1)
             + c1 * vel(i,j,k+1,1)
             + c2 * vel(i,j,k+2,1)) * odz;
         wz = (c0 * vel(i,j,k  ,2)
             + c1 * vel(i,j,k+1,2)
             + c2 * vel(i,j,k+2,2)) * odz;
      }
   }
   else
   {
      uz = 0.5 * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * odz;
      vz = 0.5 * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * odz;
      wz = 0.5 * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * odz;
   }

}

#endif
