#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>

#include <mfix_actual_rhs.H>
#include <mfix_reporter.H>
#include <integrator_parameters.H>
#include <mfix_chem_memory_handler.H>
#include <mfix_fix_inputs.H>


/**
 * \brief Part of the solution of a linear algebra problem
 */
template <bool allow_pivot>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static
void dgesl (ChemMatrix<amrex::Real>& a,
            const uint64_t N,
            ChemVector<short>& pivot,
            ChemVector<amrex::Real>& b)
{
  // solve a * x = b
  // first solve l * m_y = b
  if (N > 1) {
    for (uint64_t k = 1; k < N; ++k) {

      amrex::Real t{};
      if constexpr (allow_pivot) {
        uint64_t l = pivot(k);
        t = b(l);
        if (l != k) {
          b(l) = b(k);
          b(k) = t;
        }
      } else {
        t = b(k);
      }

      for (uint64_t j = k+1; j <= N; ++j) {
        b(j) += t * a(j,k);
      }
    }
  }

  // now solve u * x = m_y
  for (uint64_t kb = 1; kb <= N; ++kb) {

    uint64_t k = N + 1 - kb;
    b(k) = b(k) / a(k,k);
    amrex::Real t = -b(k);
    for (uint64_t j = 1; j <= k-1; ++j) {
      b(j) += t * a(j,k);
    }
  }
}


/**
 * \brief Part of the solution of a linear algebra problem
 */
template <bool allow_pivot>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static
void dgefa (ChemMatrix<amrex::Real>& a,
            const uint64_t N,
            ChemVector<short>& pivot,
            int& info)
{
  // dgefa factors a matrix by gaussian elimination.
  // a is returned in the form a = l * u where
  // l is a product of permutation and unit lower
  // triangular matrices and u is upper triangular.

  // gaussian elimination with partial pivoting

  info = 0;

  amrex::Real t;

  if (N > 1) {

    for (uint64_t k = 1; k < N; ++k) {

      // find l = pivot index
      uint64_t l = k;

      if constexpr (allow_pivot) {
        amrex::Real dmax = std::abs(a(k,k));
        for (uint64_t i = k+1; i <= N; ++i) {
          if (std::abs(a(i,k)) > dmax) {
            l = i;
            dmax = std::abs(a(i,k));
          }
        }

        pivot(k) = static_cast<short>(l);
      }

      // zero pivot implies this column already triangularized
      if (a(l,k) != 0.0e0) {

        if constexpr (allow_pivot) {
          // interchange if necessary
          if (l != k) {
            t = a(l,k);
            a(l,k) = a(k,k);
            a(k,k) = t;
          }
        }

        // compute multipliers
        t = -1.0e0 / a(k,k);
        for (uint64_t j = k+1; j <= N; ++j) {
          a(j,k) *= t;
        }

        // row elimination with column indexing
        for (uint64_t j = k+1; j <= N; ++j) {
          t = a(l,j);

          if constexpr (allow_pivot) {
            if (l != k) {
              a(l,j) = a(k,j);
              a(k,j) = t;
            }
          }

          for (uint64_t i = k+1; i <= N; ++i) {
            a(i,j) += t * a(i,k);
          }
        }
      } else {
        info = k;
      }
    }
  }

  if constexpr (allow_pivot) {
    pivot(N) = static_cast<short>(N);
  }

  if (a(N,N) == 0.0e0) {
    info = N;
  }
}


#endif
