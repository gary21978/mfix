// Do not include this header anywhere other than AMReX_GpuLaunchMacrosG.H.
// The purpose of this file is to avoid clang-tidy.

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    { auto const& amrex_i_tn = TN; \
    if (!amrex::isEmpty(amrex_i_tn)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_tn); \
        AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS, amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const TI : amrex::Gpu::Range(amrex_i_tn)) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI : amrex::Gpu::Range(amrex_i_tn)) { \
            block \
        } \
    }}}

// two fused launches
#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    { auto const& amrex_i_tn1 = TN1; auto const& amrex_i_tn2 = TN2; \
    if (!amrex::isEmpty(amrex_i_tn1) || !amrex::isEmpty(amrex_i_tn2)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(amrex_i_tn1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(amrex_i_tn2); \
        dim3 amrex_i_nblocks = amrex::max(amrex_i_ec1.numBlocks.x, \
                                          amrex_i_ec2.numBlocks.x); \
        amrex_i_nblocks.y = 2; \
        AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS, amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
                        block2 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
            block2 \
        } \
    }}}

// three fused launches
#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    { auto const& amrex_i_tn1 = TN1; auto const& amrex_i_tn2 = TN2; auto const& amrex_i_tn3 = TN3; \
    if (!amrex::isEmpty(amrex_i_tn1) || !amrex::isEmpty(amrex_i_tn2) || !amrex::isEmpty(amrex_i_tn3)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(amrex_i_tn1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(amrex_i_tn2); \
        const auto amrex_i_ec3 = amrex::Gpu::ExecutionConfig(amrex_i_tn3); \
        dim3 amrex_i_nblocks = amrex::max(amrex::max(amrex_i_ec1.numBlocks.x, \
                                                     amrex_i_ec2.numBlocks.x), \
                                                     amrex_i_ec3.numBlocks.x); \
        amrex_i_nblocks.y = 3; \
        AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS, amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
                        block2 \
                    } \
                    break; \
            case 2: for (auto const TI3 : amrex::Gpu::Range(amrex_i_tn3)) { \
                        block3 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
            block2 \
        } \
        for (auto const TI3 : amrex::Gpu::Range(amrex_i_tn3)) { \
            block3 \
        } \
    }}}

#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    { auto const& amrex_i_tn = TN; \
    if (!amrex::isEmpty(amrex_i_tn)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_tn); \
        AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS, amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const TI : amrex::Gpu::Range(amrex_i_tn)) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE: cannot call device function from host"); \
    }}}

// two fused launches
#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    { auto const& amrex_i_tn1 = TN1; auto const& amrex_i_tn2 = TN2; \
    if (!amrex::isEmpty(amrex_i_tn1) || !amrex::isEmpty(amrex_i_tn2)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(amrex_i_tn1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(amrex_i_tn2); \
        dim3 amrex_i_nblocks = amrex::max(amrex_i_ec1.numBlocks.x, \
                                          amrex_i_ec2.numBlocks.x); \
        amrex_i_nblocks.y = 2; \
        AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS, amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
                        block2 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }}}

// three fused launches
#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    { auto const& amrex_i_tn1 = TN1; auto const& amrex_i_tn2 = TN2; auto const& amrex_i_tn3 = TN3; \
    if (!amrex::isEmpty(amrex_i_tn1) || !amrex::isEmpty(amrex_i_tn2) || !amrex::isEmpty(amrex_i_tn3)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(amrex_i_tn1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(amrex_i_tn2); \
        const auto amrex_i_ec3 = amrex::Gpu::ExecutionConfig(amrex_i_tn3); \
        dim3 amrex_i_nblocks = amrex::max(amrex::max(amrex_i_ec1.numBlocks.x, \
                                                     amrex_i_ec2.numBlocks.x), \
                                                     amrex_i_ec3.numBlocks.x); \
        amrex_i_nblocks.y = 3; \
        AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS, amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
                        block2 \
                    } \
            case 2: for (auto const TI3 : amrex::Gpu::Range(amrex_i_tn3)) { \
                        block3 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }}}

// FOR_1D

#define AMREX_HOST_DEVICE_FOR_1D(n,i,block) \
{ \
    auto const& amrex_i_n = n; \
    using amrex_i_inttype = std::remove_const_t<decltype(n)>; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_n,[=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept block); \
    } else { \
        auto amrex_i_lambda = [=] (amrex_i_inttype i) noexcept block; \
        for (amrex_i_inttype i = 0; i < amrex_i_n; ++i) amrex_i_lambda(i); \
    } \
}

#define AMREX_HOST_DEVICE_PARALLEL_FOR_1D(n,i,block) \
{ \
    auto const& amrex_i_n = n; \
    using amrex_i_inttype = std::remove_const_t<decltype(n)>; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_n,[=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept block); \
    } else { \
        auto amrex_i_lambda = [=] (amrex_i_inttype i) noexcept block; \
        for (amrex_i_inttype i = 0; i < amrex_i_n; ++i) amrex_i_lambda(i); \
    } \
}

#define AMREX_FOR_1D(n,i,block) \
{ \
    using amrex_i_inttype = std::remove_const_t<decltype(n)>; \
    amrex::ParallelFor(n,[=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept block); \
}

#define AMREX_PARALLEL_FOR_1D(n,i,block) \
{ \
    using amrex_i_inttype = std::remove_const_t<decltype(n)>; \
    amrex::ParallelFor(n,[=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept block); \
}

// FOR_3D

#define AMREX_HOST_DEVICE_FOR_3D(box,i,j,k,block) \
{ \
    auto const& amrex_i_box = box; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept block); \
    } else { \
        amrex::LoopConcurrentOnCpu(amrex_i_box,[=] (int i, int j, int k) noexcept block); \
    } \
}

#define AMREX_HOST_DEVICE_PARALLEL_FOR_3D(box,i,j,k,block) \
{ \
    auto const& amrex_i_box = box; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept block); \
    } else { \
        amrex::LoopConcurrentOnCpu(amrex_i_box,[=] (int i, int j, int k) noexcept block); \
    } \
}

#define AMREX_FOR_3D(box,i,j,k,block) \
{ \
    amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept block); \
}

#define AMREX_PARALLEL_FOR_3D(box,i,j,k,block) \
{ \
    amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept block); \
}

// FOR_4D

#define AMREX_HOST_DEVICE_FOR_4D(box,ncomp,i,j,k,n,block) \
{ \
    auto const& amrex_i_box = box; \
    auto const& amrex_i_ncomp = ncomp; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_box,amrex_i_ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept block); \
    } else { \
        amrex::LoopConcurrentOnCpu(amrex_i_box,amrex_i_ncomp,[=] (int i, int j, int k, int n) noexcept block); \
    } \
}

#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D(box,ncomp,i,j,k,n,block) \
{ \
    auto const& amrex_i_box = box; \
    auto const& amrex_i_ncomp = ncomp; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_box,amrex_i_ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept block); \
    } else { \
        amrex::LoopConcurrentOnCpu(amrex_i_box,amrex_i_ncomp,[=] (int i, int j, int k, int n) noexcept block); \
    } \
}

#define AMREX_FOR_4D(box,ncomp,i,j,k,n,block) \
{ \
    amrex::ParallelFor(box,ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept block); \
}

#define AMREX_PARALLEL_FOR_4D(box,ncomp,i,j,k,n,block) \
{ \
    amrex::ParallelFor(box,ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept block); \
}
