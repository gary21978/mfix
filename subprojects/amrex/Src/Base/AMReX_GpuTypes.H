#ifndef AMREX_GPU_TYPES_H_
#define AMREX_GPU_TYPES_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>

#ifdef AMREX_USE_GPU

namespace amrex {

}

#endif

namespace amrex::Gpu {

#if defined(AMREX_USE_GPU)
struct Handler
{
    AMREX_GPU_HOST_DEVICE constexpr explicit Handler (int n_active_threads = -1)
        : numActiveThreads(n_active_threads) {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isFullBlock () const {
        return (numActiveThreads >= (int)blockDim.x ||
                numActiveThreads <= 0);
    }

    int numActiveThreads;
};
#else
struct Handler {};
#endif

}

#endif
