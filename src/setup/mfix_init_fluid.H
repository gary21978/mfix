#ifndef __MFIX_INIT_FLUID_HPP_
#define __MFIX_INIT_FLUID_HPP_

#include <mfix.H>
#include <mfix_ic.H>
#include <mfix_fluid.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_GpuContainers.H>

void set_ic_fluid ( int const a_lev,
                    const amrex::Box& sbx,
                    const amrex::Box& bx,
                    const amrex::Box& domain,
                    const amrex::MFIter& mfi,
                    MFIXLevelData& a_ld,
                    const amrex::Real dx,
                    const amrex::Real dy,
                    const amrex::Real dz,
                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                    bool test_tracer_conservation,
                    MFIXInitialConditions& initial_conditions,
                    MFIXFluidPhase& fluid);

void init_periodic_tracer (const Box& bx,
                           const Box& domain,
                           amrex::Array4<Real> const& a_vel,
                           amrex::Array4<Real> const& a_trac,
                           const Real dx,
                           const Real dy,
                           const Real dz);

void init_fluid_parameters ( int const a_lev,
                             const Box& bx,
                             const MFIter& mfi,
                             MFIXLevelData& a_ld,
                             MFIXFluidPhase& fluid);

void set_ic (const Box& sbx,
             const Box& domain,
             const Real dx,
             const Real dy,
             const Real dz,
             FArrayBox& vel_g);

#endif
