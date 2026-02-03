#ifndef _MFIX_CALCCELL_H_
#define _MFIX_CALCCELL_H_

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>

void calc_cell_ic (const amrex::Real dx,
                   const amrex::Real dy,
                   const amrex::Real dz,
                   const amrex::Real* lo,
                   const amrex::Real* hi,
                   const amrex::Real* plo,
                   int& i_w,
                   int& i_e,
                   int& j_s,
                   int& j_n,
                   int& k_b,
                   int& k_t);



const amrex::Box calc_ic_box (const amrex::Geometry& geom,
                              const amrex::RealBox* region);

const amrex::Box calc_bc_box (amrex::Geometry const& a_geom,
                              amrex::RealBox  const* a_region,
                              int             const  a_dir=-1);

#endif
