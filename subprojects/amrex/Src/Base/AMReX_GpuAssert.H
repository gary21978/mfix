#ifndef AMREX_GPU_ASSERT_H_
#define AMREX_GPU_ASSERT_H_
#include <AMReX_Config.H>

#if !defined(__APPLE__)

#include <cassert>

#if defined(AMREX_USE_CUDA)
#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#endif

#endif  // !defined(__APPLE__)
#endif  // AMREX_GPU_ASSERT_H_
