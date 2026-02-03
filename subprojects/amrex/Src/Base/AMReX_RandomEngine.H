#ifndef AMREX_RANDOM_ENGINE_H_
#define AMREX_RANDOM_ENGINE_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>

#if defined(AMREX_USE_CUDA)
#include <curand.h>
#include <curand_kernel.h>
#endif

namespace amrex
{
#ifdef AMREX_USE_GPU

    using randState_t = curandState_t;
    using randGenerator_t = curandGenerator_t;

    extern randState_t* gpu_rand_state;

    inline
    randState_t* getRandState () {
        return gpu_rand_state;
    }

    struct RandomEngine {
        randState_t* rand_state;
        AMREX_GPU_HOST_DEVICE RandomEngine (randState_t* rs) : rand_state(rs) {}
    };

#else

    struct RandomEngine {}; // CPU

#endif

    AMREX_FORCE_INLINE
    RandomEngine getInvalidRandomEngine () {
#ifdef AMREX_USE_GPU
        return RandomEngine{nullptr};
#else
        return RandomEngine{};
#endif
    }

}

#endif
