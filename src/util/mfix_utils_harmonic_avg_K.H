#ifndef MFIX_UTILS_HARMONIC_AVG_K_H_
#define MFIX_UTILS_HARMONIC_AVG_K_H_


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void harmonic_avg_cc2face_x ( amrex::Box const& xbx, int ncomp,
                              amrex::Array4<amrex::Real      > const& phi_x,
                              amrex::Array4<amrex::Real const> const& phi,
                              int const periodic_x,
                              const amrex::Box& domain ) noexcept
{
  const amrex::Dim3 dlo = amrex::lbound(domain);
  const amrex::Dim3 dhi = amrex::ubound(domain);

  amrex::Loop(xbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if ( !periodic_x && (i == dlo.x) ) {

      phi_x(i,j,k,n) = phi(dlo.x-1,j,k,n);

    } else if ( !periodic_x && (i == dhi.x+1) ) {

      phi_x(i,j,k,n) = phi(dhi.x+1,j,k,n);

    } else {

      amrex::Real const kP(phi(i-1,j,k,n));
      amrex::Real const kE(phi(i  ,j,k,n));

      if ( amrex::almostEqual(kP, amrex::Real(0.0)) ||
           amrex::almostEqual(kE, amrex::Real(0.0)) ) {

        phi_x(i,j,k,n) = amrex::Real(0.0);

      } else {

        phi_x(i,j,k,n) = amrex::Real(2.0)*(kP*kE)/(kP + kE);

      }
    }
  });
}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_harmonic_avg_cc2face_x ( amrex::Box const& xbx, int ncomp,
                                 amrex::Array4<amrex::Real      > const& phi_x,
                                 amrex::Array4<amrex::Real const> const& phi,
                                 amrex::Array4<amrex::Real const> const& apx,
                                 amrex::Array4<amrex::Real const> const& cvol,
                                 amrex::Array4<amrex::Real const> const& ccent,
                                 int const periodic_x,
                                 const amrex::Box& domain ) noexcept
{
  const amrex::Dim3 dlo = amrex::lbound(domain);
  const amrex::Dim3 dhi = amrex::ubound(domain);

  amrex::Loop(xbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (amrex::almostEqual(apx(i,j,k), amrex::Real(0.0))) {

      phi_x(i,j,k) = amrex::Real(0.0);

    } else if ( !periodic_x && (i == dlo.x) ) {

      phi_x(i,j,k,n) = phi(dlo.x-1,j,k,n);

    } else if ( !periodic_x && (i == dhi.x+1) ) {

      phi_x(i,j,k,n) = phi(dhi.x+1,j,k,n);

    } else {

      amrex::Real const kP(phi(i-1,j,k,n));
      amrex::Real const kE(phi(i  ,j,k,n));

      if ( amrex::almostEqual(kP, amrex::Real(0.0)) ||
           amrex::almostEqual(kE, amrex::Real(0.0)) ) {

        phi_x(i,j,k,n) = amrex::Real(0.0);

      } else {

        if ( amrex::almostEqual( apx(i  ,j,k), amrex::Real(1.0)) &&
             amrex::almostEqual(cvol(i  ,j,k), amrex::Real(1.0)) &&
             amrex::almostEqual(cvol(i-1,j,k), amrex::Real(1.0)) ) {

          phi_x(i,j,k,n) = amrex::Real(2.0)*(kP*kE)/(kP + kE);

        } else {

          amrex::Real const d0  = amrex::Real(0.5) - ccent(i-1,j,k,0);
          amrex::Real const d1  = amrex::Real(0.5) + ccent(i  ,j,k,0);
          amrex::Real const d01 = (d0 + d1);

          if (amrex::almostEqual(d01, amrex::Real(0.0))) {

            phi_x(i,j,k,n) = amrex::Real(0.5)*(kP + kE);

          } else {

            amrex::Real const fe = d1 / (d0 + d1);

            phi_x(i,j,k,n) = amrex::Real(1.0) / ((amrex::Real(1.0) - fe)/kP + fe/kE);

          }
        }
      }
    }

  });
} // end eb_harmonic_avg_x



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void harmonic_avg_cc2face_y ( amrex::Box const& ybx, int ncomp,
                              amrex::Array4<amrex::Real      > const& phi_y,
                              amrex::Array4<amrex::Real const> const& phi,
                              int const periodic_y,
                              const amrex::Box& domain ) noexcept
{
  const amrex::Dim3 dlo = amrex::lbound(domain);
  const amrex::Dim3 dhi = amrex::ubound(domain);

  amrex::Loop(ybx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {

    if ( !periodic_y && (j == dlo.y) ) {

      phi_y(i,j,k,n) = phi(i,dlo.y-1,k,n);

    } else if ( !periodic_y && (j == dhi.y+1) ) {

      phi_y(i,j,k,n) = phi(i,dhi.y+1,k,n);

    } else {

      amrex::Real const kP(phi(i,j-1,k,n));
      amrex::Real const kN(phi(i,j  ,k,n));

      if ( amrex::almostEqual(kP, amrex::Real(0.0)) ||
           amrex::almostEqual(kN, amrex::Real(0.0)) ) {

        phi_y(i,j,k,n) = amrex::Real(0.0);

      } else {

        phi_y(i,j,k,n) = amrex::Real(2.0)*(kP*kN)/(kP + kN);

      }
    }
  });
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_harmonic_avg_cc2face_y ( amrex::Box const& ybx, int ncomp,
                                 amrex::Array4<amrex::Real      > const& phi_y,
                                 amrex::Array4<amrex::Real const> const& phi,
                                 amrex::Array4<amrex::Real const> const& apy,
                                 amrex::Array4<amrex::Real const> const& cvol,
                                 amrex::Array4<amrex::Real const> const& ccent,
                                 int const periodic_y,
                                 const amrex::Box& domain ) noexcept
{
  const amrex::Dim3 dlo = amrex::lbound(domain);
  const amrex::Dim3 dhi = amrex::ubound(domain);

  amrex::Loop(ybx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (amrex::almostEqual(apy(i,j,k), amrex::Real(0.0))) {

      phi_y(i,j,k) = amrex::Real(0.0);

    } else if ( !periodic_y && (j == dlo.y) ) {

      phi_y(i,j,k,n) = phi(i,dlo.y-1,k,n);

    } else if ( !periodic_y && (j == dhi.y+1) ) {

      phi_y(i,j,k,n) = phi(i,dhi.y+1,k,n);

    } else {

      amrex::Real const kP(phi(i,j-1,k,n));
      amrex::Real const kN(phi(i,j  ,k,n));

      if ( amrex::almostEqual(kP, amrex::Real(0.0)) ||
           amrex::almostEqual(kN, amrex::Real(0.0)) ) {

        phi_y(i,j,k,n) = amrex::Real(0.0);

      } else {

        if ( amrex::almostEqual( apy(i,j  ,k), amrex::Real(1.0)) &&
             amrex::almostEqual(cvol(i,j  ,k), amrex::Real(1.0)) &&
             amrex::almostEqual(cvol(i,j-1,k), amrex::Real(1.0)) ) {

          phi_y(i,j,k,n) = amrex::Real(2.0)*(kP*kN)/(kP + kN);

        } else {

          amrex::Real const d0  = amrex::Real(0.5) - ccent(i,j-1,k,1);
          amrex::Real const d1  = amrex::Real(0.5) + ccent(i,j  ,k,1);
          amrex::Real const d01 = (d0 + d1);

          if (amrex::almostEqual(d01, amrex::Real(0.0))) {

            phi_y(i,j,k,n) = amrex::Real(0.5)*(kP + kN);

          } else {

            amrex::Real const fn = d1 / (d0 + d1);

            phi_y(i,j,k,n) = amrex::Real(1.0) / ((amrex::Real(1.0) - fn)/kP + fn/kN);

          }
        }
      }
    }

  });
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void harmonic_avg_cc2face_z ( amrex::Box const& zbx, int ncomp,
                              amrex::Array4<amrex::Real      > const& phi_z,
                              amrex::Array4<amrex::Real const> const& phi,
                              int const periodic_z,
                              const amrex::Box& domain ) noexcept
{
  const amrex::Dim3 dlo = amrex::lbound(domain);
  const amrex::Dim3 dhi = amrex::ubound(domain);

  amrex::Loop(zbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {

    if ( !periodic_z && (k == dlo.z) ) {

      phi_z(i,j,k,n) = phi(i,j,dlo.z-1,n);

    } else if ( !periodic_z && (j == dhi.z+1) ) {

      phi_z(i,j,k,n) = phi(i,j,dhi.z+1,n);

    } else {

      amrex::Real const kP(phi(i,j,k-1,n));
      amrex::Real const kT(phi(i,j,k  ,n));

      if ( amrex::almostEqual(kP, amrex::Real(0.0)) ||
           amrex::almostEqual(kT, amrex::Real(0.0)) ) {

        phi_z(i,j,k,n) = amrex::Real(0.0);

      } else {

        phi_z(i,j,k,n) = amrex::Real(2.0)*(kP*kT)/(kP + kT);

      }
    }
  });
}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_harmonic_avg_cc2face_z ( amrex::Box const& zbx, int ncomp,
                                 amrex::Array4<amrex::Real      > const& phi_z,
                                 amrex::Array4<amrex::Real const> const& phi,
                                 amrex::Array4<amrex::Real const> const& apz,
                                 amrex::Array4<amrex::Real const> const& cvol,
                                 amrex::Array4<amrex::Real const> const& ccent,
                                 int const periodic_z,
                                 const amrex::Box& domain ) noexcept
{
  const amrex::Dim3 dlo = amrex::lbound(domain);
  const amrex::Dim3 dhi = amrex::ubound(domain);

  amrex::Loop(zbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (amrex::almostEqual(apz(i,j,k), amrex::Real(0.0))) {

      phi_z(i,j,k) = amrex::Real(0.0);

    } else if ( !periodic_z && (k == dlo.z) ) {

      phi_z(i,j,k,n) = phi(i,j,dlo.z-1,n);

    } else if ( !periodic_z && (k == dhi.z+1) ) {

      phi_z(i,j,k,n) = phi(i,j,dhi.z+1,n);

    } else {

      amrex::Real const kP(phi(i,j,k-1,n));
      amrex::Real const kT(phi(i,j,k  ,n));

      if ( amrex::almostEqual(kP, amrex::Real(0.0)) ||
           amrex::almostEqual(kT, amrex::Real(0.0)) ) {

        phi_z(i,j,k,n) = amrex::Real(0.0);

      } else {

        if ( amrex::almostEqual( apz(i,j,k  ), amrex::Real(1.0)) &&
             amrex::almostEqual(cvol(i,j,k  ), amrex::Real(1.0)) &&
             amrex::almostEqual(cvol(i,j,k-1), amrex::Real(1.0)) ) {

          phi_z(i,j,k,n) = amrex::Real(2.0)*(kP*kT)/(kP + kT);

        } else {

          amrex::Real const d0  = amrex::Real(0.5) - ccent(i,j,k-1,2);
          amrex::Real const d1  = amrex::Real(0.5) + ccent(i,j,k  ,2);
          amrex::Real const d01 = (d0 + d1);

          if (amrex::almostEqual(d01, amrex::Real(0.0))) {

            phi_z(i,j,k,n) = amrex::Real(0.5)*(kP + kT);

          } else {

            amrex::Real const ft = d1 / (d0 + d1);

            phi_z(i,j,k,n) = amrex::Real(1.0) / ((amrex::Real(1.0) - ft)/kP + ft/kT);

          }
        }
      }
    }
  });
} // end eb_harmonic_avg_z

#endif
