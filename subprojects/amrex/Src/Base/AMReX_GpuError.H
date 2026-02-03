#ifndef AMREX_GPU_ERROR_H_
#define AMREX_GPU_ERROR_H_
#include <AMReX_Config.H>

#include <AMReX_Utility.H>
#include <string>

#if defined(AMREX_USE_GPU) && !defined(AMREX_GPU_NO_ERROR_CHECK)

namespace amrex {

using gpuError_t = cudaError_t;
constexpr gpuError_t gpuSuccess = cudaSuccess;

inline gpuError_t
gpuGetLastError () {
    return cudaGetLastError();
}

inline const char*
gpuGetErrorString (gpuError_t error) {
    return cudaGetErrorString(error);
}

namespace Gpu {

    inline void ErrorCheck (const char* file, int line) noexcept
    {
        gpuError_t err = gpuGetLastError();
        if (gpuSuccess != err) {
            std::string errStr(std::string("GPU last error detected in file ") + file
                               + " line " + std::to_string(line)
                               + ": " + std::string(gpuGetErrorString(err)));
            amrex::Abort(errStr);
        }
    }
}
}

#define AMREX_GPU_SAFE_CALL(call) { \
    amrex::gpuError_t amrex_i_err = call;          \
    if (amrex::gpuSuccess != amrex_i_err) { \
        std::string errStr(std::string("GPU error in file ") + __FILE__ \
                           + " line " + std::to_string(__LINE__) \
                           + " " + amrex::gpuGetErrorString(amrex_i_err)); \
        amrex::Abort(errStr); \
    }}

#ifdef AMREX_USE_CUDA
#define AMREX_CUDA_SAFE_CALL(call) { \
    cudaError_t amrex_i_err = call;          \
    if (cudaSuccess != amrex_i_err) { \
        std::string errStr(std::string("CUDA error ") + std::to_string(amrex_i_err) \
                           + std::string(" in file ") + __FILE__ \
                           + " line " + std::to_string(__LINE__) \
                           + ": " + cudaGetErrorString(amrex_i_err)); \
        amrex::Abort(errStr); \
    }}

#define AMREX_CURAND_SAFE_CALL(call) { \
    curandStatus_t amrex_i_err = call;          \
    if (CURAND_STATUS_SUCCESS != amrex_i_err) { \
        std::string errStr(std::string("CURAND error ") + std::to_string(amrex_i_err) \
                           + std::string(" in file ") + __FILE__ \
                           + " line " + std::to_string(__LINE__)); \
        amrex::Abort(errStr); \
    }}

#define AMREX_CUFFT_SAFE_CALL(call) { \
    cufftResult_t amrex_i_err = call; \
    if (CUFFT_SUCCESS != amrex_i_err) { \
        std::string errStr(std::string("CUFFT error ")+std::to_string(amrex_i_err) \
                           + std::string(" in file ") + __FILE__        \
                           + " line " + std::to_string(__LINE__));      \
        amrex::Abort(errStr); \
    }}

#endif

#define AMREX_GPU_ERROR_CHECK() amrex::Gpu::ErrorCheck(__FILE__, __LINE__)

#else

#define AMREX_GPU_SAFE_CALL(call)  (call)
#define AMREX_GPU_ERROR_CHECK()    ((void)0)

#endif

#endif
