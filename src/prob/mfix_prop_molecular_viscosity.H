#ifndef MFIX_THERMO_VISCOSITY_HEAT_H_
#define MFIX_THERMO_VISCOSITY_HEAT_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>

struct MolViscosityModel { enum {Undefined, Constant, Sutherland, Reid, Usr}; };

struct MolViscosityData
{

  MolViscosityData ( int const a_model, int const a_ncomp,
                     amrex::Real const* const a_MW,
                     amrex::Real const* const a_coeffs )
    : m_model(a_model), m_ncomp(a_ncomp), p_MW(a_MW), p_coeffs(a_coeffs)
    , m_is_mixture(a_ncomp > 1)
    , m_has_temperature((a_model == MolViscosityModel::Sutherland) ||
                        (a_model == MolViscosityModel::Reid))
  { }

  ///
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real compute_species ( int const a_k, amrex::Real const a_T ) const noexcept
  {

    if (m_model == MolViscosityModel::Constant) {

      return p_coeffs[a_k];

    } else if (m_model == MolViscosityModel::Sutherland) {

      // Species offset for coeff index: C1 and S
      int const k( 2*a_k );

      // viscosity := (C1 * T^(3/2)) / (T + S)
      return (p_coeffs[k+0] * std::pow(a_T, (3./2.))) / (a_T + p_coeffs[k+1]);

    } else if (m_model == MolViscosityModel::Reid) {

      // Species offset for coeff index: A, B, C, D
      int const k( 4*a_k );

      // viscosity := A*exp( B/T + C*T + D*T^2)
      return (p_coeffs[k+0] * std::exp( (p_coeffs[k+1]/a_T) +
                 (p_coeffs[k+2]*a_T) + (p_coeffs[k+3]*a_T*a_T)));
    }
    amrex::Abort("Invalid molecular viscosity model.");
    return std::numeric_limits<amrex::Real>::infinity();
  }


  ///
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real compute ( const amrex::Real a_T,
                        const amrex::CellData<amrex::Real const>& a_X ) const noexcept
  {
    amrex::Real mu(0.);

    if (!m_is_mixture) {

      mu = compute_species(0, a_T);

    } else {

      amrex::Real MW_mix(0.);
      for (int n(0); n < m_ncomp; ++n) {
        MW_mix += a_X[n] / p_MW[n];
      }
      MW_mix = amrex::Real(1.) / MW_mix;


      for (int m(0); m < m_ncomp; ++m) {

        amrex::Real const mw_m = p_MW[m];
        amrex::Real const mu_m = compute_species(m, a_T);

        amrex::Real y_phi = 0.;

        for (int n(0); n < m_ncomp; ++n) {

          amrex::Real const mw_n = p_MW[n];
          amrex::Real const mu_n = compute_species(n, a_T);

          amrex::Real phi(1.);

          if (m != n) {

            amrex::Real const mu_m_mu_n_1o2 = std::sqrt(mu_m/mu_n);
            amrex::Real const mw_n_mu_m_1o4 = std::sqrt(std::sqrt(mw_n/mw_m));
            amrex::Real const mw_m_mw_n = mw_m / mw_n;

            phi = std::pow((1. + mu_m_mu_n_1o2 * mw_n_mu_m_1o4), 2) / std::sqrt(8.*(1. + mw_m_mw_n));
          }

          y_phi += (a_X[n] * MW_mix/mw_n) * phi;
        }
        mu += (a_X[m] * MW_mix/mw_m) * (mu_m / y_phi);
      }

    } // mixture formulation

    return mu;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( int i, int j, int k,
                           amrex::Array4<amrex::Real const> const& a_T,
                           amrex::Array4<amrex::Real const> const& a_X ) const noexcept
  {
    amrex::Real l_T = (m_has_temperature) ? a_T(i,j,k) : 0.;
    auto l_X = (m_is_mixture) ?  a_X.cellData(i,j,k) : amrex::CellData<amrex::Real const>{nullptr, 0, 0};

    return compute(l_T, l_X);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( const amrex::Real a_T,
                           amrex::IntVect const& a_cell,
                           amrex::Array4<amrex::Real const> const& a_X ) const noexcept
  {
    amrex::Real l_T = (m_has_temperature) ? a_T : 0.;
    auto l_X = (m_is_mixture) ?  a_X.cellData(a_cell[0],a_cell[1],a_cell[2]) :
      amrex::CellData<amrex::Real const>{nullptr, 0, 0};

    return compute(l_T, l_X);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( const amrex::Real* a_fluid_interp,
                           const int a_T_idx,
                           const int a_X_idx ) const noexcept
  {
    amrex::Real l_T = (m_has_temperature) ? a_fluid_interp[a_T_idx] : 0.;
    auto l_X = (m_is_mixture) ?
      amrex::CellData<amrex::Real const>{&a_fluid_interp[a_X_idx], 1, m_ncomp} :
      amrex::CellData<amrex::Real const>{nullptr, 0, 0};

    return compute(l_T, l_X);
  }


  public:

    int const m_model;

    int const m_ncomp;

    const amrex::Real* p_MW;
    const amrex::Real* p_coeffs;

    bool m_is_mixture;
    bool m_has_temperature;
};


struct MolViscosity
{

  public:

    MolViscosity() = default;

    int define ( const std::string& base);

  public:

    int m_model = MolViscosityModel::Undefined;

    int m_ncomp = 0;
    int m_defined = 0;

    // Model specific viscosity coefficients
    amrex::Gpu::HostVector<amrex::Real>   m_h_coeffs;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_coeffs;

};

#endif
