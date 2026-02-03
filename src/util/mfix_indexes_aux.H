#ifndef MFIX_INDEXES_AUX_H
#define MFIX_INDEXES_AUX_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>


/**
 * \brief Get the 1D index associated to the 3D index ijk of a Box
 *
 * Given a Box with dimensions \f$ \Delta x \f$, \f$ \Delta y \f$, and \f$
 * \Delta z \f$, and a set of tridimensional indexes i, j, and k, compute the
 * corresponding monodimensional index n by the formula:
 * \f[
 *   n = i + j \cdot \Delta x + k \cdot \Delta x \cdot \Delta y
 * \f]
 *
 * This function is the inverse of the function get_local_indexes.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int get_global_index (const amrex::IntVect& ijk,
                      const amrex::Box& bx)
{
  return (ijk[0]-bx.smallEnd(0)) +
         (ijk[1]-bx.smallEnd(1))*bx.length(0) +
         (ijk[2]-bx.smallEnd(2))*bx.length(0)*bx.length(1);
}


/**
 * \brief Get the 3D index ijk associated to the 1D index n of a Box
 *
 * Given a Box with dimensions \f$ \Delta x \f$, \f$ \Delta y \f$, and \f$
 * \Delta z \f$, and a monodimensional index n, compute the corresponding set of
 * tridimensional indexes i, j, and k by the formula:
 * \f[
 * \left\{ \begin{array}{l}
 *   k = \lfloor \dfrac{n}{\Delta x \cdot \Delta y} \rfloor \\ [1.5em]
 *   j = \lfloor \dfrac{n - k \cdot \Delta x \cdot \Delta y}
 *       {\Delta x} \rfloor \\ [1.5em]
 *   i = n - k \cdot \Delta x \cdot \Delta y - j \cdot \Delta x
 * \end{array} \right.
 * \f]
 *
 * This function is the inverse of the function get_global_index.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::IntVect get_local_indexes (const int& n,
                                  const amrex::Box& bx)
{
  amrex::IntVect ijk(-1);

  int k = int(amrex::Math::floor(n / (bx.length(0)*bx.length(1))));
  int j = int(amrex::Math::floor((n - k*bx.length(0)*bx.length(1)) / bx.length(0)));
  int i = n - k*bx.length(0)*bx.length(1) - j*bx.length(0);

  ijk[2] = k + bx.smallEnd(2);
  ijk[1] = j + bx.smallEnd(1);
  ijk[0] = i + bx.smallEnd(0);

  return ijk;
}


#endif
