#ifndef MFIX_DISTRIBUTIONS_UNIFORM_H_
#define MFIX_DISTRIBUTIONS_UNIFORM_H_

namespace distributions {

struct uniform {

  public:

    uniform ( amrex::Real const a_min, amrex::Real const a_max,
             int const a_is_number_weighted)
    {
      if (a_is_number_weighted) { setup_Nw_distribution(a_min, a_max); }
      else { setup_Vw_distribution(a_min, a_max); }
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

    void setup_Nw_distribution (amrex::Real const a_min,
                                amrex::Real const a_max)
    {
      amrex::Real const b(a_max);
      amrex::Real const a(a_min);

      m_Nw_mean = 0.5*(b + a);
      m_Nw_variance = (b - a)*(b - a)/12.0;

      amrex::Real const a4(a*a*a*a), b4(b*b*b*b);
      amrex::Real const a5(a*a4),    b5(b*b4);
      amrex::Real const a6(a*a5),    b6(b*b5);

      m_Vw_mean = (4.0*(b5 - a5)) / (5.0*(b4 - a4));

      amrex::Real const m(m_Vw_mean);

      m_Vw_variance = 4.0*( (b6-a6)/6. - 2.*m*(b5-a5)/5. + m*m*(b4-a4)/4. ) / (b4-a4);

    }


    void setup_Vw_distribution ( amrex::Real const a_min,
                                 amrex::Real const a_max)
    {
      amrex::Real const b(a_max);
      amrex::Real const a(a_min);

      m_Vw_mean = 0.5*(b + a);
      m_Vw_variance = (b - a)*(b - a)/12.0;

      amrex::Real const a2(a*a), b2(b*b);

      m_Nw_mean = (2.*(a*b)) / (b + a);

      amrex::Real const m(m_Nw_mean);

      m_Nw_variance = (2.*a2*b2*std::log(b/a) - 4.*m*a*b*(b-a)) / (b2-a2) + m*m;

    }

};

} // end namespace distribution

#endif
