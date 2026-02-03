#include <AMReX_ParmParse.H>

#include <mfix_distributions.H>
#include <mfix_distributions_normal.H>
#include <mfix_distributions_log_normal.H>
#include <mfix_distributions_uniform.H>
#include <mfix_distributions_custom.H>

using namespace amrex;

int
INPUT_DIST_t::
define (std::string const dist_prefix,
        std::string const prop)
{
  amrex::ParmParse pp(dist_prefix);

  std::string dist_prop_prefix = dist_prefix + "." + prop;
  amrex::ParmParse ppDist(dist_prop_prefix);

  // Get the distribution type.
  std::string distribution;
  if (!pp.query(prop, distribution)) {
    //m_error = "Distribution type not specified.";
    return 1;
  }

  int fnd_mean = ppDist.query("mean", m_mean);
  int fnd_std  = ppDist.query("std",  m_stddev);
  int fnd_min  = ppDist.query("min",  m_min);
  int fnd_max  = ppDist.query("max",  m_max);

  ppDist.query("bins", m_bins);

  if( amrex::toLower(distribution) == "constant") {

    m_is_constant = 1;
    m_bins = 1;

    int fnd_const = ppDist.query("constant", m_mean);

    if (!fnd_const || (fnd_const && fnd_mean) || m_mean <= 0.) {

      //m_error = "Invalid constant distribution parameters.";
      return 1;

    }

    m_stddev = 0.;
    m_max = m_mean;
    m_min = m_mean;

  } else if (amrex::toLower(distribution) ==  "normal") {

    m_is_normal = 1;

    if (!fnd_mean || !fnd_std || !fnd_min || !fnd_max ||
       (m_max <= m_min) || (m_bins < 1)) {

      //m_error = "Invalid normal distribution parameters.";
      return 1;

    }

  } else if (amrex::toLower(distribution) == "log-normal") {

    m_is_log_normal = 1;

    if (!fnd_mean || !fnd_std || !fnd_min || !fnd_max ||
       (m_bins < 1) || (m_max <= m_min) || (m_min < 0.)) {

      //m_error = "Invalid log-normal distribution parameters.";
      return 1;

    } else {

      // min and max of normal distribution
      m_log_min = std::log(m_min);
      m_log_max = std::log(m_max);

    }

  } else if (amrex::toLower(distribution) ==  "uniform") {

    m_is_uniform = 1;

    if (!fnd_min || !fnd_max || (m_max <= m_min) || (m_bins < 1) ) {

      //m_error = "Invalid uniform distribution parameters.";
      return 1;

    }

  } else if (amrex::toLower(distribution) == "custom") {

    m_is_custom = 1;

    std::string custom_input;
    if ( !ppDist.query("custom", custom_input) || (m_bins < 1) ) {

      //m_error = "Invalid custom distribution parameters.";
      return 1;

    }

    ppDist.query("interpolate", m_interp_custom);

  } else {

    //m_error = "Invalid distribution type: " + distribution + ".";
    return 1;

  }

  return 0;
}


int
INPUT_DIST_t::
define (std::string const dist_prefix,
        std::string const prop,
        INPUT_DIST_DATA_t* a_data,
        int const a_use_number_weighted)
{
  if ( define (dist_prefix, prop)) { return 1; }

  m_use_volume_weighted = 1-a_use_number_weighted;

  std::string dist_prop_prefix = dist_prefix + "." + prop;
  amrex::ParmParse ppDist(dist_prop_prefix);

  std::string weighting;
  if (ppDist.query("type", weighting)) {
    m_is_number_weighted = (amrex::toLower(weighting).compare("number-weighted") == 0);
    m_is_volume_weighted = (amrex::toLower(weighting).compare("volume-weighted") == 0);
  }
  if (!m_is_number_weighted && !m_is_volume_weighted) {
    if (!m_is_constant) { return 1; }
  }

  if ( m_is_constant ) {

    m_Vw_mean = m_mean;
    m_Vw_stddev = 0.0;

  } else if ( m_is_normal ) {

    distributions::normal ndist(m_mean, m_stddev, m_is_number_weighted);

    m_Vw_mean = ndist.Vw_mean();
    m_Vw_stddev = ndist.Vw_stddev();

    if (m_use_volume_weighted) {
      m_mean = m_Vw_mean;
      m_stddev = m_Vw_stddev;
    } else {
      m_mean = ndist.Nw_mean();
      m_stddev = ndist.Nw_stddev();
    }

  } else if ( m_is_log_normal ) {

    distributions::log_normal lndist(m_mean, m_stddev, m_is_number_weighted);

    m_Vw_mean = lndist.Vw_mean();
    m_Vw_stddev = lndist.Vw_stddev();

    if (m_use_volume_weighted) {
      m_mean = m_Vw_mean;
      m_stddev = m_Vw_stddev;
    } else {
      m_mean = lndist.Nw_mean();
      m_stddev = lndist.Nw_stddev();
    }

  } else if ( m_is_uniform ) {

    distributions::uniform udist(m_min, m_max, m_is_number_weighted);

    m_Vw_mean = udist.Vw_mean();
    m_Vw_stddev = udist.Vw_stddev();

    if (m_use_volume_weighted) {
      m_mean = m_Vw_mean;
      m_stddev = m_Vw_stddev;
    } else {
      m_mean = udist.Nw_mean();
      m_stddev = udist.Nw_stddev();
    }

  } else if ( m_is_custom ) {

    std::string custom_input;
    ppDist.get("custom", custom_input);

    distributions::custom cdist(custom_input);
    if (!cdist.ok()) { return 1; }

    m_bins = cdist.bins();

    a_data->m_h_diameters.resize(m_bins, 0.);

    a_data->m_h_NwCDF.resize(m_bins, 0.);
    a_data->m_h_VwCDF.resize(m_bins, 0.);

    m_h_ptr_diameters = a_data->m_h_diameters.data();

    m_h_ptr_NwCDF = a_data->m_h_NwCDF.data();
    m_h_ptr_VwCDF = a_data->m_h_VwCDF.data();

    cdist.setup_distribution(custom_input,
          m_bins, m_is_number_weighted,
          a_data->m_h_diameters.data(),
          a_data->m_h_NwCDF.data(),
          a_data->m_h_VwCDF.data());

    if (!cdist.ok()) { return 1; }

    if (m_interp_custom && !cdist.can_interp() ) { return 1; }

    if (m_use_volume_weighted) {
      m_mean = cdist.Vw_mean();
      m_stddev = cdist.Vw_stddev();
    } else {
      m_mean = cdist.Nw_mean();
      m_stddev = cdist.Nw_stddev();
    }
    m_min = (m_min < 0.) ? cdist.min() : amrex::max(cdist.min(), m_min);
    m_max = (m_max < 0.) ? cdist.max() : amrex::min(cdist.max(), m_max);
  }

  return 0;
}


void
INPUT_DIST_t::
copyAsync ( INPUT_DIST_DATA_t* a_data )
{
  // The following only applied when using custom distributions.
  if (!is_custom()) { return; }

  // diameters: host -> device copy and set pointers
  // ............................................................
  a_data->m_d_diameters.resize(m_bins);

  Gpu::copyAsync(Gpu::hostToDevice, a_data->m_h_diameters.begin(),
    a_data->m_h_diameters.end(), a_data->m_d_diameters.begin());

  m_h_ptr_diameters = a_data->m_h_diameters.data();
  m_d_ptr_diameters = a_data->m_d_diameters.data();

  // number-weighted CDF: host->device copy and set pointers
  // ............................................................
  a_data->m_d_NwCDF.resize(m_bins);

  Gpu::copyAsync(Gpu::hostToDevice, a_data->m_h_NwCDF.begin(),
    a_data->m_h_NwCDF.end(), a_data->m_d_NwCDF.begin());

  m_h_ptr_NwCDF = a_data->m_h_NwCDF.data();
  m_d_ptr_NwCDF = a_data->m_d_NwCDF.data();


  // volume-weighted CDF: host->device copy and set pointers
  // ............................................................
  a_data->m_d_VwCDF.resize(m_bins);

  Gpu::copyAsync(Gpu::hostToDevice, a_data->m_h_VwCDF.begin(),
    a_data->m_h_VwCDF.end(), a_data->m_d_VwCDF.begin());

  m_h_ptr_VwCDF = a_data->m_h_VwCDF.data();
  m_d_ptr_VwCDF = a_data->m_d_VwCDF.data();

  if (m_use_volume_weighted) {
    m_h_ptr_CDF = m_h_ptr_VwCDF;
    m_d_ptr_CDF = m_d_ptr_VwCDF;
  } else {
    m_h_ptr_CDF = m_h_ptr_NwCDF;
    m_d_ptr_CDF = m_d_ptr_NwCDF;
  }
}


void
INPUT_DIST_t::
set_sample_pointers ( int const a_sample_size,
                      Gpu::HostVector<amrex::Real>& a_h_samples,
                      Gpu::DeviceVector<amrex::Real>& a_d_samples)
{
  AMREX_ALWAYS_ASSERT( a_sample_size > 0);
  m_sample_size = a_sample_size;

  a_d_samples.resize(m_sample_size);

  amrex::Gpu::copy(amrex::Gpu::hostToDevice,
    a_h_samples.begin(), a_h_samples.end(), a_d_samples.begin());

  m_h_ptr_samples = a_h_samples.data();
  m_d_ptr_samples = a_d_samples.data();
}


