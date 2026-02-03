#ifndef MFIX_LEASTSQUARES_K_H_
#define MFIX_LEASTSQUARES_K_H_

#include <AMReX.H>
#include <AMReX_Config.H>
#include <cmath>
#include <AMReX_LUSolver.H>

#define SMALLERROR 1.e-10


//---------------------------------------------------------------------
template <int maxPts>
AMREX_GPU_DEVICE AMREX_INLINE
amrex::Array1D<amrex::Real,0,3> LeastSquares (
                  amrex::Array2D<amrex::Real,0,maxPts-1, 0,2> const& xi
                , amrex::Array1D<amrex::Real,0,maxPts-1> const& b
                , amrex::Array1D<amrex::Real,0,2> const& xt
                , int const& numPts)
//---------------------------------------------------------------------
{ // Four coeff polynomial
  int const n = 4; // a0 + a1*x + a2*y + a3*z

  amrex::Array1D<amrex::Real,0,n-1> value;

  // A.x = b
  amrex::Array2D<amrex::Real,0,n-1,0,n-1,amrex::Order::C> Alocal;
  amrex::Real blocal[n];

  for (int jj = 0; jj < n; ++jj) {
      for (int ii = 0; ii < n; ++ii) {
          Alocal(jj,ii) = 0.;
      }
      blocal[jj] = 0.;
  }

  for (int lc(0); lc < numPts; ++lc) {

    amrex::Real dx = xi(lc,0) - xt(0);
    amrex::Real dy = xi(lc,1) - xt(1);
    amrex::Real dz = xi(lc,2) - xt(2);

    Alocal(0,0) += 1.;
    Alocal(0,1) += dx;
    Alocal(0,2) += dy;
    Alocal(0,3) += dz;
    Alocal(1,1) += dx * dx;
    Alocal(1,2) += dx * dy;
    Alocal(1,3) += dx * dz;
    Alocal(2,2) += dy * dy;
    Alocal(2,3) += dy * dz;
    Alocal(3,3) += dz * dz;

    blocal[0] += b(lc);
    blocal[1] += b(lc) * dx;
    blocal[2] += b(lc) * dy;
    blocal[3] += b(lc) * dz;
  }

  for (int i(0); i < n; i++) {
    for (int j = i; j < n; j++) {
      if (i != j) { Alocal(j,i) = Alocal(i,j); }
    }
  }

  // Solve the system of equations
  amrex::LUSolver<n,amrex::Real> lusolver(Alocal);
  double xlocal[n];
  lusolver(xlocal, blocal);

  for (int lc(0); lc < n; ++lc) {
      value(lc) = xlocal[lc];
  }

  return (value);
}

#endif
