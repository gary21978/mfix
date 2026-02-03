#ifndef MFIX_THERMO_SPECIFIC_HEAT_H_
#define MFIX_THERMO_SPECIFIC_HEAT_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Algorithm.H>

struct SpecificHeatModel { enum {None, Constant, NASA7Polynomials}; };

/**
 * \brief Functions for computing specific heat (J/kg.T).
 *
 *   Cp = a0 * a1*T + a2*T^2 + a3*T^3 + a4*T^4
 *      = a0 + T*(a1 + T*(a2 + T*(a3 + T*a4)))
 */
struct SpecificHeatData
{

  SpecificHeatData ( int const a_nspecies, amrex::Real const* const a_coeffs,
      int const* const a_T_offset, int const* const a_coeff_offset,
      amrex::Real const* const a_T)
    : m_nspecies(a_nspecies) , p_coeffs(a_coeffs), p_T_offset(a_T_offset)
    , p_coeff_offset(a_coeff_offset), p_T(a_T)
  { }

  /// Compute the specific heat of the k-th species at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( int const a_k, amrex::Real const a_T) const noexcept
  {
    // Locate temperature range
    amrex::Real const* ub = amrex::upper_bound(p_T + p_T_offset[a_k],  p_T + p_T_offset[a_k+1], a_T);
    int const t( std::distance(p_T + p_T_offset[a_k], ub));

    // Offset for species coefficients
    int const k(p_coeff_offset[a_k] + 7*t);

    // a0 * a1*T + a2*T^2 + a3*T^3 + a4*T^4
    return (p_coeffs[k+0] + a_T*(p_coeffs[k+1] + a_T*(p_coeffs[k+2]
                  + a_T*(p_coeffs[k+3] + a_T*p_coeffs[k+4]))));
  }

  /// Compute the specific heat of mixture at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::IntVect a_cell,
                           amrex::Array4<amrex::Real const> const& a_T,
                           amrex::Array4<amrex::Real const> const& a_X ) const noexcept
  {
    amrex::Real Cp = operator()(0, a_T(a_cell));
    if (m_nspecies > 1) { Cp *= a_X(a_cell,0); }
    for (int k(1); k < m_nspecies; ++k) {
      Cp += a_X(a_cell,k) * operator()(k,a_T(a_cell));
    }
    return Cp;
  }

  /// Compute the specific heat of mixture at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::Real const a_T, amrex::Real const* const a_X) const noexcept
  {
    amrex::Real Cp = operator()(0,a_T);
    if (m_nspecies > 1) { Cp *= a_X[0]; }
    for (int k(1); k < m_nspecies; ++k) {
      Cp += a_X[k] * operator()(k,a_T);
    }
    return Cp;
  }

  /// Compute the specific heat of mixture at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::Real const a_T,
                           amrex::IntVect a_cell,
                           amrex::Array4<amrex::Real const> const& a_X ) const noexcept
  {
    amrex::Real Cp = operator()(0, a_T);
    if (m_nspecies > 1) { Cp *= a_X(a_cell,0); }
    for (int k(1); k < m_nspecies; ++k) {
      Cp += a_X(a_cell,k) * operator()(k,a_T);
    }
    return Cp;
  }

  /// Compute the specific heat of mixture at temperature T for particle
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::Real const a_T,
                           int const a_i,
                           amrex::ParticleReal* AMREX_RESTRICT * AMREX_RESTRICT a_ptile_runtime_rdata,
                           int const a_X_idx ) const noexcept
  {
    amrex::Real Cp = operator()(0, a_T);
    if (m_nspecies > 1) { Cp *= a_ptile_runtime_rdata[a_X_idx][a_i]; }
    for (int k(1); k < m_nspecies; ++k) {
      Cp += a_ptile_runtime_rdata[a_X_idx+k][a_i] * operator()(k,a_T);
    }
    return Cp;
  }

  /// Compute the specific heat of mixture from interpolated fluid data
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( const amrex::Real* a_fluid_interp, int const a_T_idx, int const a_X_idx) const noexcept
  {
    amrex::Real Cp = operator()(0,a_fluid_interp[a_T_idx]);
    if (m_nspecies > 1) { Cp *= a_fluid_interp[a_X_idx]; }
    for (int k(1); k < m_nspecies; ++k) {
      Cp += a_fluid_interp[a_X_idx+k] * operator()(k,a_fluid_interp[a_T_idx]);
    }
    return Cp;
  }

  public:

    int const m_nspecies;

    const amrex::Real* p_coeffs;

    const int* p_T_offset;
    const int* p_coeff_offset;
    const amrex::Real* p_T;
};

/**
 * \brief Functions for computing enthalpy.
 *
 *   H = a0*T * a1*T^2 + a2*T^3 + a3*T^4 + a4*T^5 + a5
 *     = T*(a0 + T*(a1/2 + T*(a2/3 + T*(a3/4 + T*a4/5)))) + a5
 */
struct EnthalpyData
{

  EnthalpyData ( int const a_nspecies, amrex::Real const* const a_coeffs,
      int const* const a_T_offset, int const* const a_coeff_offset,
      amrex::Real const* const a_T)
    : m_nspecies(a_nspecies) , p_coeffs(a_coeffs), p_T_offset(a_T_offset)
    , p_coeff_offset(a_coeff_offset), p_T(a_T)
  { }

  /// Compute the specific enthalpy of the k-th species at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( int const a_k, amrex::Real const a_T, int const a_cell_is_covered = 0) const noexcept
  {
    if (a_cell_is_covered) { return 0.; }

    // Locate temperature range
    amrex::Real const* ub = amrex::upper_bound(p_T + p_T_offset[a_k],  p_T + p_T_offset[a_k+1], a_T);
    int const t( std::distance(p_T + p_T_offset[a_k], ub));

    // Offset for species coefficients
    int const k( p_coeff_offset[a_k] + 7*t);

    return (a_T*(p_coeffs[k+0] + a_T*(p_coeffs[k+1]/2. + a_T*(p_coeffs[k+2]/3.
               + a_T*(p_coeffs[k+3]/4. + a_T*p_coeffs[k+4]/5.)))) + p_coeffs[k+5]);
  }

  /// Compute the specific enthalpy of mixture at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::IntVect a_cell,
                           amrex::Array4<amrex::Real const> const& a_T,
                           amrex::Array4<amrex::Real const> const& a_X,
                           int const a_cell_is_covered = 0) const noexcept
  {
    if (a_cell_is_covered) { return 0.; }

    amrex::Real H = operator()(0, a_T(a_cell));
    if (m_nspecies > 1) {  H *= a_X(a_cell,0); }
    for (int k(1); k < m_nspecies; ++k) {
      H += a_X(a_cell,k) * operator()(k,a_T(a_cell));
    }
    return H;
  }

  /// Compute the specific enthalpy of mixture at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::Real const a_T, amrex::Real const* const a_X, int const a_cell_is_covered = 0) const noexcept
  {
    if (a_cell_is_covered) { return 0.; }

    amrex::Real H = operator()(0,a_T);
    if ((m_nspecies > 1) && (a_X != nullptr)) {  H *= a_X[0]; }
    for (int k(1); k < m_nspecies; ++k) {
      H += a_X[k] * operator()(k,a_T);
    }
    return H;
  }

  /// Compute the specific enthalpy of mixture at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::Real const a_T,
                           amrex::IntVect a_cell,
                           amrex::Array4<amrex::Real const> const& a_X,
                           int const a_cell_is_covered = 0) const noexcept
  {
    if (a_cell_is_covered) { return 0.; }

    amrex::Real H = operator()(0, a_T);
    if (m_nspecies > 1) {  H *= a_X(a_cell,0); }
    for (int k(1); k < m_nspecies; ++k) {
      H += a_X(a_cell,k) * operator()(k,a_T);
    }
    return H;
  }

  /// Compute the specific enthalpy of mixture at temperature T for particle
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::Real const a_T,
                           int const a_i,
                           amrex::ParticleReal* AMREX_RESTRICT * AMREX_RESTRICT a_ptile_runtime_rdata,
                           int const a_X_idx) const noexcept
  {
    amrex::Real H = operator()(0, a_T);
    if (m_nspecies > 1) {  H *= a_ptile_runtime_rdata[a_X_idx][a_i]; }
    for (int k(1); k < m_nspecies; ++k) {
      H += a_ptile_runtime_rdata[a_X_idx+k][a_i] * operator()(k,a_T);
    }
    return H;
  }

  public:

    int const m_nspecies;

    const amrex::Real* p_coeffs;

    const int* p_T_offset;
    const int* p_coeff_offset;
    const amrex::Real* p_T;
};


class SpecificHeat
{

  public:

    SpecificHeat() = default;

    int define ( std::string base );

    int ReadT ( std::string base,
                amrex::Vector<std::string> a_species );

    int ReadTEntry ( std::string base,
                     int const comp );

    int ReadConstant ( amrex::Real const a_Tref, std::string base,
                       amrex::Vector<std::string> a_species );

    int ReadConstantEntry ( amrex::Real const a_Tref, std::string base,
                            int const comp );

    int ReadNASA7 ( std::string base,
                    amrex::Vector<std::string> a_species );

    int ReadNASA7Entry ( std::string base, int const comp);

    int CheckEntry ( int const a_ignore, int const comp,
                     std::string base, EnthalpyData& enthalpy,
                     SpecificHeatData& specificHeat);

  public:


    int m_model = SpecificHeatModel::None;

    int m_nspecies = 0;
    int m_ncomp = 0;

    // NASA-7 species specific heat coefficient
    amrex::Gpu::HostVector<amrex::Real>   m_h_coeffs;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_coeffs;

    // NASA-7 species temperatures
    amrex::Gpu::HostVector<int>   m_h_T_offset;
    amrex::Gpu::DeviceVector<int> m_d_T_offset;

    amrex::Gpu::HostVector<int>   m_h_coeff_offset;
    amrex::Gpu::DeviceVector<int> m_d_coeff_offset;

    amrex::Gpu::HostVector<amrex::Real>   m_h_T;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_T;
};

#endif
