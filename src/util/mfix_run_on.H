#ifndef _MFIX_RUN_ON_H_
#define _MFIX_RUN_ON_H_

#include <AMReX_GpuControl.H>


#ifdef AMREX_USE_GPU
static constexpr amrex::RunOn run_on = amrex::RunOn::Device;
#else
static constexpr amrex::RunOn run_on = amrex::RunOn::Host;
#endif


#endif
