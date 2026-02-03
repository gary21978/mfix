#ifndef MFIX_DISTRIBUTIONS_LOG_NORMAL_H_
#define MFIX_DISTRIBUTIONS_LOG_NORMAL_H_

namespace distributions {

struct log_normal {

  public:

    log_normal ( amrex::Real const a_mean, amrex::Real const a_stddev,
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

      m_Vw_mean = a_mean + 3.0*a_variance;
      m_Vw_variance = a_variance;
    }

    void setup_Vw_distribution ( amrex::Real const a_mean,
                                 amrex::Real const a_variance)
    {
      m_Vw_mean = a_mean;
      m_Vw_variance  = a_variance;

      m_Nw_mean = a_mean - 3.0*a_variance;
      m_Nw_variance = a_variance;
    }
};

} // end namespace distribution

#endif
