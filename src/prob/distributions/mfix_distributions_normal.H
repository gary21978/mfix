#ifndef MFIX_DISTRIBUTIONS_NORMAL_H_
#define MFIX_DISTRIBUTIONS_NORMAL_H_

namespace distributions {

struct normal {

  public:

    normal ( amrex::Real const a_mean, amrex::Real const a_stddev,
             int const a_is_number_weighted)
    {
      amrex::Real const variance(a_stddev*a_stddev);
      if (a_is_number_weighted) { setup_Nw_distribution(a_mean, variance); }
      else { setup_Vw_distribution(a_mean, variance); }
    }

    amrex::Real Nw_mean () { return m_Nw_mean; }
    amrex::Real Nw_variance () { return m_Nw_variance; }
    amrex::Real Nw_stddev () { return std::sqrt(m_Nw_variance); }

    amrex::Real Vw_mean () { return m_Vw_mean; }
    amrex::Real Vw_variance () { return m_Vw_variance; }
    amrex::Real Vw_stddev () { return std::sqrt(m_Vw_variance); }

  private:

    amrex::Real m_Nw_mean;
    amrex::Real m_Nw_variance;

    amrex::Real m_Vw_mean;
    amrex::Real m_Vw_variance;


    void setup_Nw_distribution (amrex::Real const a_mean,
                                amrex::Real const a_variance)
    {
      m_Nw_mean = a_mean;
      m_Nw_variance = a_variance;

      m_Vw_mean = g1({m_Nw_mean, m_Nw_variance}) /
                  g0({m_Nw_mean, m_Nw_variance});

      m_m = m_Vw_mean;

      m_Vw_variance  = g2({m_Nw_mean, m_Nw_variance}) /
                       g0({m_Nw_mean, m_Nw_variance});
    }


    void setup_Vw_distribution ( amrex::Real const a_mean,
                                 amrex::Real const a_variance)
    {
      m_Vw_mean = a_mean;
      m_Vw_variance  = a_variance;

      m_m = a_mean;
      m_v = a_variance;

      // Compute initial guess for the mean and variance for the
      // number-weighted distribution.

      amrex::Real const m(m_Vw_mean), v(m_Vw_variance);
      amrex::Real const stddev(std::sqrt(v));

      // Pick the min and max bounds to integrate over.
      amrex::Real const smin((m - std::sqrt(m*m - 12.*v))*0.5);
      amrex::Real const smax(m + std::sqrt(m*m - 12.*v) - smin);

      int const bins = 1024;
      amrex::Real const db( (smax - smin) / (static_cast<amrex::Real>(bins)));

      amrex::Real const db_2sqrt3 = (0.5*db)/std::sqrt(3);

      amrex::Array<amrex::Real,2> x {0., 0.};

      amrex::Real sum_PDF(0.);

      for (int bin(0); bin<bins; ++bin) {

        amrex::Real const midpoint = smin + db*(static_cast<amrex::Real>(bin)+0.5);

        // Use two-point Gaussian quadrature to integrate over the interval.
        amrex::Real const xl = midpoint - db_2sqrt3;
        amrex::Real const xh = midpoint + db_2sqrt3;

        amrex::Real const al = (xl - m_Vw_mean)/stddev;
        amrex::Real const ah = (xh - m_Vw_mean)/stddev;

        amrex::Real const p = 0.5*db*(std::exp(-0.5*al*al)/(xl*xl*xl)
                             + std::exp(-0.5*ah*ah)/(xh*xh*xh));

        x[0] += p*midpoint;
        x[1] += p*midpoint*midpoint;

        sum_PDF += p;
      }
      AMREX_ALWAYS_ASSERT(sum_PDF > 0.);

      x[0] /= sum_PDF;
      x[1] = x[1]/sum_PDF - x[0]*x[0];

      // Using the mean and variance computed from sampling the distribution,
      // find the roots of the polynomial using a homotopy method for solving
      // non-linear systems.

      int const nsteps = 4;
      m_dt = 1./static_cast<amrex::Real>(nsteps);
      m_Fx0 = {f1(x), f2(x)};
      for( int n(0); n<nsteps; ++n) { x = step(x); }

      m_Nw_mean = x[0];
      m_Nw_variance = x[1];

    }

    /***********************************************************************
    *                                                                      *
    * The following functions comprise the integrated normal distribution  *
    * function multiplied by x^3:                                          *
    *                                                                      *
    *  n(x) = ((x^3)/(sigma * sqrt(2*pi))) * exp( -0.5*((x - mu)/sigma)^2) *
    *                                                                      *
    *  mu    : mean of number-weighted distribution                        *
    *  sigma : standard deviation of number-weighted distribution          *
    *                                                                      *
    *  g0(x) = int n(x) from  -/+ infinity                                 *
    *  g1(x) = int x*n(x) from -/+ infinity                                *
    *  g2(x) = int (x^2)*n(x) form -/+ infinity                            *
    *                                                                      *
    *  g0    : number density of volume-weighted distribution              *
    *  g1/g0 : mean of volume-weighted distribution                        *
    *  g2/g0 : variance of volume-weighted distribution                    *
    *                                                                      *
    ***********************************************************************/

    amrex::Real m_v, m_m;

    // g0 = h(x) integrated over -/+ infinity
    amrex::Real g0 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0]), x1(a_x[1]);
      return (3.*x0*x1 + x0*x0*x0);
    }

    // partial derivate of g0 with respect to x0 (mean)
    amrex::Real dg0dx0 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0]), x1(a_x[1]);
      return 3.*(x0*x0 + x1);
    }

    // partial derivate of g0 with respect to x1 (variance)
    amrex::Real dg0dx1 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0])/*, x1(a_x[1])*/;
      return 3.*x0;
    }

    // g1 = x*h(x) integrated over -/+ infinity
    amrex::Real g1 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0]), x1(a_x[1]);
      return (3.*x1*x1 + 6.*x0*x0*x1 + x0*x0*x0*x0);
    }

    // partial derivate of g1 with respect to x0 (mean)
    amrex::Real dg1dx0 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0]), x1(a_x[1]);
      return 4.*x0*(x0*x0 + 3.*x1);
    }

    // partial derivate of g1 with respect to x1 (variance)
    amrex::Real dg1dx1 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0]), x1(a_x[1]);
      return 6.*(x0*x0 + x1);
    }


    // g2 = x^2*h(x) integrated over -/+ infinity
    amrex::Real g2 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0]), x1(a_x[1]);
      amrex::Real const b(x0 - m_m);
      return 3.*x1*x1*(5.*x0 - 2.*m_m) + x0*x1*(x0*x0 + 6.*x0*b + 3.*b*b) + x0*x0*x0*b*b;
    }

    // partial derivate of g2 with respect to x0 (mean)
    amrex::Real dg2dx0 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0]), x1(a_x[1]);
      amrex::Real const b(x0 - m_m);
      return 15.*x1*x1 + x1*(9.*x0*x0 + 18.*x0*b + 3.*b*b) + x0*x0*b*(2.*x0 + 3.*b);
    }

    // partial derivate of g2 with respect to x1 (variance)
    amrex::Real dg2dx1 ( amrex::Array<amrex::Real,2> a_x )
    { amrex::Real const x0(a_x[0]), x1(a_x[1]);
      amrex::Real const b(x0 - m_m);
      return 6.*x1*(5.*x0 - 2.*m_m) + x0*(x0*x0 + 6.*x0*b + 3.*b*b);
    }


    // f1(x) = g1 - (volume-weighted mean)*g1 = 0
    amrex::Real f1 ( amrex::Array<amrex::Real,2> a_x )
    { return g1(a_x) - m_m*g0(a_x); }

    // partial derivate of f1 with respect to x0 (mean)
    amrex::Real df1dx0 ( amrex::Array<amrex::Real,2> a_x )
    { return dg1dx0(a_x) - m_m*dg0dx0(a_x); }

    // partial derivate of f1 with respect to x1 (variance)
    amrex::Real df1dx1 ( amrex::Array<amrex::Real,2> a_x )
    { return dg1dx1(a_x) - m_m*dg0dx1(a_x); }


    // f2(x) = g2 - (volume-weighted variance)*g2 = 0
    amrex::Real f2 ( amrex::Array<amrex::Real,2> a_x )
    { return g2(a_x) - m_v*g0(a_x); }

    // partial derivate of f2 with respect to x0 (mean)
    amrex::Real df2dx0 ( amrex::Array<amrex::Real,2> a_x )
    { return dg2dx0(a_x) - m_v*dg0dx0(a_x); }

    // partial derivate of f2 with respect to x1 (variance)
    amrex::Real df2dx1 ( amrex::Array<amrex::Real,2> a_x )
    { return dg2dx1(a_x) - m_v*dg0dx1(a_x); }


    /**********************************************************************
    *                                                                     *
    * The following methods implement a homotopy method to solve the non- *
    * linear system of equations. Find x0 and x1 (mean and variance) of   *
    * the number-weighted distribution that corresponds :                 *
    *                                                                     *
    *  f1(x0,x1) = g1(x0,x1) - m * g0(x0,x1) = 0                          *
    *  f2(x0,x1) = g2(x0,x1) - v * g0(x0,x1) = 0                          *
    *                                                                     *
    *  where m and v are the mean and variance of the volume-weighted     *
    *  distribution, and x0 and x1 are the mean and variance of the       *
    *  number-weighted distribution.                                      *
    *                                                                     *
    **********************************************************************/

    // Step size and value of F evaluate at initial x
    amrex::Real m_dt;
    amrex::Array<amrex::Real,2> m_Fx0;


    amrex::Array<amrex::Real,2> solve_leq ( amrex::Array<amrex::Real,2> const a_x)
    {
      amrex::Array<amrex::Array<amrex::Real,2>,2> J, invJ;

      J[0][0] = df1dx0(a_x); J[0][1] = df1dx1(a_x);
      J[1][0] = df2dx0(a_x); J[1][1] = df2dx1(a_x);

      amrex::Real const detJ = (J[0][0]*J[1][1] - J[0][1]*J[1][0]);

      if (amrex::almostEqual(detJ,0.))
      { amrex::Abort("Invalid matrix -- det invJ is zero!\n"); }

      invJ[0][0] =  (J[1][1])/detJ;
      invJ[0][1] = -(J[0][1])/detJ;

      invJ[1][0] = -(J[1][0])/detJ;
      invJ[1][1] =  (J[0][0])/detJ;

      return { -m_dt*(invJ[0][0]*m_Fx0[0] + invJ[0][1]*m_Fx0[1]),
               -m_dt*(invJ[1][0]*m_Fx0[0] + invJ[1][1]*m_Fx0[1]) };
    }


    amrex::Array<amrex::Real,2> step ( amrex::Array<amrex::Real,2> const a_x )
    {
      amrex::Array<amrex::Real,2> w, k1, k2, k3, k4;

      k1 = solve_leq(a_x);
      for (int i(0); i<2; ++i) { w[i] = a_x[i] + 0.5*k1[i]; }

      k2 = solve_leq(w);
      for (int i(0); i<2; ++i) { w[i] = a_x[i] + 0.5*k2[i]; }

      k3 = solve_leq(w);
      for (int i(0); i<2; ++i) { w[i] = a_x[i] + k3[i]; }

      k4 = solve_leq(w);

      amrex::Array<amrex::Real,2> sol;
      for( int i(0); i<2; ++i) {
        sol[i] = a_x[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
      }

      return sol;
    }


};

} // end namespace distribution

#endif
