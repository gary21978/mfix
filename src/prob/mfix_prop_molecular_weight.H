#ifndef MFIX_THERMO_MOLECULAR_WEIGHT_H_
#define MFIX_THERMO_MOLECULAR_WEIGHT_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>

/**
 * \brief Functions for computing molecular weight (mol/kg).
 *
 */
struct MolecularWeightData
{

  MolecularWeightData ( int const a_ncomp, amrex::Real const* const a_MW)
    : m_ncomp(a_ncomp) , p_MW(a_MW)
  { }

  /// Compute the specific heat of the k-th species at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( int const a_k ) const noexcept
  { return (p_MW[a_k]); }

  /// Compute the specific heat of mixture at temperature T
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() ( amrex::IntVect a_cell,
                           amrex::Array4<amrex::Real const> const& a_X ) const noexcept
  {
    amrex::Real MW(0.);

    if (m_ncomp == 1) { MW = p_MW[0]; }
    else {

      for (int n(0); n < m_ncomp; ++n) {
        MW += a_X(a_cell,n) / p_MW[n];
      }
      MW = amrex::Real(1.) / MW;
    }
    return MW;
  }


  public:

    int const m_ncomp;

    const amrex::Real* p_MW;
};


class MolecularWeight
{

  public:

    MolecularWeight() = default;

    int define ( std::string base );

    int ReadEntry ( std::string base, int const comp);

    int defined () const { return m_defined; }

  public:

    int m_ncomp = 0;

    int m_defined = 0;

    amrex::Gpu::HostVector<amrex::Real>   m_h_MW;
    amrex::Gpu::DeviceVector<amrex::Real> m_d_MW;
};

#endif
