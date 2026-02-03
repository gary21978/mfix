#ifndef _MFIX_FILCC_H_
#define _MFIX_FILCC_H_

#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>

namespace mfix_aux {

void filcc (amrex::Real* data,
            const int* lo,
            const int* hi,
            const int* dom_lo,
            const int* dom_hi,
            const amrex::Real* dx,
            const amrex::Real* grd_lo,
            const amrex::Real* time,
            const int* bc);

} // end namespace mfix_aux

#endif
