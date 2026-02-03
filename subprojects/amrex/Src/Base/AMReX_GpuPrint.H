#ifndef AMREX_GPU_PRINT_H_
#define AMREX_GPU_PRINT_H_
#include <AMReX_Config.H>

#include <cstdio>

#if defined(AMREX_USE_CUDA)
#  define AMREX_DEVICE_PRINTF(...) std::printf(__VA_ARGS__);
#else
#  define AMREX_DEVICE_PRINTF(format,...) { \
      std::printf(format, __VA_ARGS__); }
#endif

#endif  // AMREX_GPU_PRINT_H_
