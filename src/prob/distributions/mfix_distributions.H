#ifndef MFIX_DISTRIBUTIONS_PARMS_H_
#define MFIX_DISTRIBUTIONS_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_RandomEngine.H>

#include <mfix_run_on.H>

struct INPUT_DIST_DATA_t {

  public:

    INPUT_DIST_DATA_t ()
    {
      m_h_diameters.clear();
      m_d_diameters.clear();

      m_h_NwCDF.clear();
      m_d_NwCDF.clear();

      m_h_VwCDF.clear();
      m_d_VwCDF.clear();
    }

    amrex::Gpu::HostVector<amrex::Real> m_h_diameters;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_diameters;

    // Number-weighted cumulative distribution function
    amrex::Gpu::HostVector<amrex::Real> m_h_NwCDF;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_NwCDF;

    // Volume-weighted cumulative distribution function
    amrex::Gpu::HostVector<amrex::Real> m_h_VwCDF;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_VwCDF;

};

struct INPUT_DIST_t {

  public:

    AMREX_GPU_HOST_DEVICE
    INPUT_DIST_t()
      : m_is_constant(0)
      , m_is_uniform(0)
      , m_is_normal(0)
      , m_is_log_normal(0)
      , m_is_custom(0)
      , m_is_number_weighted(0)
      , m_is_volume_weighted(0)
      , m_use_volume_weighted(0)
      , m_mean(-1.0)
      , m_stddev(-1.0)
      , m_min(-1.0)
      , m_max(-1.0)
      , m_log_min(-1.0)
      , m_log_max(-1.0)
      , m_bins(64)
      , m_interp_custom(false)
      , m_h_ptr_diameters(nullptr)
      , m_d_ptr_diameters(nullptr)
      , m_h_ptr_NwCDF(nullptr)
      , m_d_ptr_NwCDF(nullptr)
      , m_h_ptr_VwCDF(nullptr)
      , m_d_ptr_VwCDF(nullptr)
      , m_h_ptr_CDF(nullptr)
      , m_d_ptr_CDF(nullptr)
      , m_sample_size(0)
      , m_h_ptr_samples(nullptr)
      , m_d_ptr_samples(nullptr)
    {}

    int define (std::string const dist_prefix,
                std::string const prop);

    int define (std::string const dist_prefix,
                std::string const prop,
                INPUT_DIST_DATA_t* a_data,
                int const a_use_number_weighted);

    void copyAsync ( INPUT_DIST_DATA_t* a_data );

    void set_sample_pointers ( int const a_samples,
             amrex::Gpu::HostVector<amrex::Real>& a_h_samples,
             amrex::Gpu::DeviceVector<amrex::Real>& a_d_samples);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int is_constant () const { return m_is_constant; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int is_uniform () const { return m_is_uniform; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int is_normal () const { return m_is_normal; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int is_log_normal () const { return m_is_log_normal; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int is_custom () const { return m_is_custom; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_mean () const {
      return ((!m_is_log_normal) ? m_mean :
        std::exp(m_mean + 0.5*m_stddev*m_stddev));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int is_number_weighted () const { return m_is_number_weighted; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int is_volume_weighted () const { return m_is_volume_weighted; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_stddev () const {
      if (!m_is_log_normal) { return m_stddev; }
      else { amrex::Real const variance= (std::exp(m_stddev*m_stddev) - 1.)
          * std::exp(2.*m_mean + m_stddev*m_stddev);
          return std::sqrt(variance);
      }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_max ()    const { return m_max; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real get_min ()    const { return m_min; }

    int get_bins () const { return m_bins; }

    std::tuple<amrex::Real, amrex::Real> get_custom_bin ( int const a_bin) const
    { AMREX_ASSERT(a_bin >= 0 && a_bin < m_bins);
      if(run_on == amrex::RunOn::Gpu)
      {return {m_h_ptr_diameters[a_bin], m_h_ptr_CDF[a_bin]}; }
      else {return {m_d_ptr_diameters[a_bin], m_d_ptr_CDF[a_bin]}; }
    }


    void report_distribution ( amrex::Real const a_scale = 1.0 ) const noexcept;


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::tuple<amrex::Real, amrex::Real> get_Vw_probability ( int const a_bin) const
    { AMREX_ASSERT(a_bin >= 0 && a_bin < m_bins);

      if (m_is_constant) {

        return {m_mean, 1.0};

      } else {

        amrex::Real p(0.);

        amrex::Real dbin = (m_max - m_min) / static_cast<amrex::Real>(m_bins);

        amrex::Real blo = m_min + static_cast<amrex::Real>(a_bin)*dbin;
        amrex::Real bhi = blo + dbin;

        if (m_is_normal) {

          amrex::Real const inv_sq2_stddev(1.0 / (std::sqrt(2)*m_Vw_stddev));

          p = std::erf( (bhi - m_Vw_mean) * inv_sq2_stddev)
              - std::erf( (blo - m_Vw_mean) * inv_sq2_stddev);

          p /= amrex::max( 10.*std::numeric_limits<amrex::Real>::min(),
                 std::erf( (m_max - m_Vw_mean) * inv_sq2_stddev)
                 - std::erf( (m_min - m_Vw_mean) * inv_sq2_stddev));

        } else if (m_is_log_normal) {

          amrex::Real const inv_sq2_stddev(1.0 / (std::sqrt(2)*m_Vw_stddev));

          p = std::erf( (std::log(bhi) - m_Vw_mean) * inv_sq2_stddev)
              - std::erf( (std::log(blo) - m_Vw_mean) * inv_sq2_stddev);

          p /= amrex::max( 10.*std::numeric_limits<amrex::Real>::min(),
                 std::erf( (std::log(m_max) - m_Vw_mean) * inv_sq2_stddev)
                 - std::erf( (std::log(m_min) - m_Vw_mean) * inv_sq2_stddev));

        } else if (m_is_uniform) {

          if ( m_is_number_weighted ) {

            p = (bhi*bhi*bhi*bhi - blo*blo*blo*blo) /
                  (m_max*m_max*m_max*m_max - m_min*m_min*m_min*m_min);

          } else {
            p = dbin / (m_max - m_min);
          }

        } else if (m_is_custom) {

          AMREX_ASSERT(m_h_ptr_diameters != nullptr);
          AMREX_ASSERT(m_h_ptr_VwCDF != nullptr);

          blo = m_h_ptr_diameters[a_bin];
          p = m_h_ptr_VwCDF[a_bin];

          if (!m_interp_custom) { bhi = blo; }
          else {
            if ( a_bin+1 < m_bins ) {
              bhi = m_h_ptr_diameters[a_bin+1];
            } else {
              p = 0.;
            }
          }
       }
       return {0.5*(blo+bhi), p};
     }

    }


    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real sample ( int const a_pid,
                         amrex::RandomEngine const& a_engine ) const {
      if ( m_sample_size > 0 ) {

        AMREX_ASSERT( 0 <= a_pid && a_pid < m_sample_size);

        if (run_on == amrex::RunOn::Gpu) { return m_d_ptr_samples[a_pid]; }
        else { return m_h_ptr_samples[a_pid]; }

      } else {

        return sample<run_on>(a_engine);
      }
      return -1.;
    }

    template <amrex::RunOn run_on>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real sample ( amrex::RandomEngine const& a_engine ) const {

      if (m_is_constant) {

        return m_mean;

      } else if (m_is_normal) {

        amrex::Real observation;

        do { observation = amrex::RandomNormal(m_mean, m_stddev, a_engine); }
        while (!(m_min <= observation && observation <= m_max));

        return observation;

      } else if (m_is_log_normal) {

        amrex::Real observation;

        do { observation = amrex::RandomNormal(m_mean, m_stddev, a_engine); }
        while (!(m_log_min <= observation && observation <= m_log_max));

        return std::exp(observation);

      } else if (m_is_uniform) {

        amrex::Real rand = amrex::Random(a_engine);

        if ( m_is_number_weighted && m_use_volume_weighted) {

          amrex::Real const b(m_max);
          amrex::Real const a(m_min);

          amrex::Real const a4(a*a*a*a), b4(b*b*b*b);

          return std::pow(a4 + rand*(b4-a4),0.25);

        } else if ( m_is_volume_weighted && !m_use_volume_weighted) {

          amrex::Real const b(m_max);
          amrex::Real const a(m_min);

          amrex::Real const a2(a*a), b2(b*b);
          return (a*b)/std::sqrt(b2 - rand*(b2-a2));

        } else {

          return m_min + (m_max-m_min)*rand;

        }

      } else if ( m_is_custom ) {

        amrex::Real* ptr_diameters;
        amrex::Real* ptr_CDF;

        if (run_on == amrex::RunOn::Cpu) {
          AMREX_ASSERT(m_h_ptr_diameters != nullptr);
          AMREX_ASSERT(m_h_ptr_CDF != nullptr);
          ptr_diameters = m_h_ptr_diameters;
          ptr_CDF = m_h_ptr_CDF;
        } else {
          AMREX_ASSERT(m_d_ptr_diameters != nullptr);
          AMREX_ASSERT(m_d_ptr_CDF != nullptr);
          ptr_diameters = m_d_ptr_diameters;
          ptr_CDF = m_d_ptr_CDF;
        }

        amrex::Real observation;

        do {

          amrex::Real const randno = amrex::Random(a_engine);

          amrex::Real* ub = amrex::upper_bound(ptr_CDF, ptr_CDF + m_bins, randno);
          int const b( std::distance(ptr_CDF, ub)-m_interp_custom);
          AMREX_ASSERT( 0 <= b && b < m_bins );
          amrex::Real pl, ph, xl, xh;
          pl = ptr_CDF[b];
          xl = ptr_diameters[b];
          observation = xl;
          if (m_interp_custom) {
            AMREX_ASSERT( b+1 < m_bins );
            ph = ptr_CDF[b+1];
            xh = ptr_diameters[b+1];
            observation += (randno-pl)*(xh-xl)/(ph-pl);
          }
        }
        while (!(m_min <= observation && observation <= m_max));

        return observation;

      } else {

        amrex::Abort("Unrecognized distribution type.");
        return std::numeric_limits<amrex::Real>::infinity();

      }
    }

    //std::string error_msg () { return m_error; }

  private:

    //std::string m_error;

    int m_is_constant;   // single valued
    int m_is_uniform;    // uniform
    int m_is_normal;     // Gaussian distribution
    int m_is_log_normal; // Log-normal distribution
    int m_is_custom;     // user defined via discrete CDF or PDF

    int m_is_number_weighted;
    int m_is_volume_weighted;

    int m_use_volume_weighted;

    amrex::Real m_mean;      // Mean value of distriubtion (or constant value)
    amrex::Real m_Vw_mean;   // Volume-weighted

    amrex::Real m_stddev;     // Standard deviation
    amrex::Real m_Vw_stddev;  // Volume-weighted

    //amrex::Real m_Nw_mean;    // Number-weighted
    //amrex::Real m_Nw_stddev;  // Number-weighted

    amrex::Real m_min;         // Minimum value of distribution
    amrex::Real m_max;         // Max value of distribution

    amrex::Real m_log_min;     // Minimum value of distribution
    amrex::Real m_log_max;     // Max value of distribution

    int m_bins;

    bool m_interp_custom;  // Interpolate between bins

    amrex::Real* m_h_ptr_diameters;
    amrex::Real* m_d_ptr_diameters;

    amrex::Real* m_h_ptr_NwCDF;
    amrex::Real* m_d_ptr_NwCDF;

    amrex::Real* m_h_ptr_VwCDF;
    amrex::Real* m_d_ptr_VwCDF;

    amrex::Real* m_h_ptr_CDF;
    amrex::Real* m_d_ptr_CDF;

    // Distriubiton samples created during the particle generation. These
    // are pointers to host and device vectors containing the radii of
    // each randomly seeded particle.
    int m_sample_size;

    amrex::Real* m_h_ptr_samples;
    amrex::Real* m_d_ptr_samples;

};

#endif
