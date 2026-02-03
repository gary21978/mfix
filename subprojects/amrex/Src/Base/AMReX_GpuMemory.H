#ifndef AMREX_GPU_MEMORY_H_
#define AMREX_GPU_MEMORY_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Arena.H>
#include <cstdlib>

namespace amrex::Gpu {

struct Managed {

#ifdef AMREX_USE_GPU

    [[nodiscard]] void *operator new (std::size_t len)
    {
        return The_Managed_Arena()->alloc(len);
    }

    void operator delete (void *ptr)
    {
        if (ptr) { The_Managed_Arena()->free(ptr); }
    }

#endif
};

struct Pinned {

#ifdef AMREX_USE_GPU

    [[nodiscard]] void *operator new (std::size_t len)
    {
        return The_Pinned_Arena()->alloc(len);
    }

    void operator delete (void *ptr)
    {
        if (ptr) { The_Pinned_Arena()->free(ptr); }
    }

#endif
};

struct Deleter {
    Arena* m_arena = nullptr;
    Deleter (Arena* ar) noexcept : m_arena(ar) {}
    void operator() (void* pt) const noexcept { m_arena->free(pt); }
};

template <class T, std::enable_if_t<std::is_trivially_copyable_v<T>,int> = 0>
struct DeviceScalar
{
    DeviceScalar (DeviceScalar const&) = delete;
    DeviceScalar (DeviceScalar &&) = delete;
    void operator= (DeviceScalar const&) = delete;
    void operator= (DeviceScalar &&) = delete;

#ifdef AMREX_USE_GPU

    DeviceScalar () {
        if (Gpu::inLaunchRegion()) {
            dp = (T*)(The_Arena()->alloc(sizeof(T)));
        } else {
            dp = (T*)(std::malloc(sizeof(T)));
        }
    }

    explicit DeviceScalar (T init_val) {
        if (Gpu::inLaunchRegion()) {
            dp = (T*)(The_Arena()->alloc(sizeof(T)));
            Gpu::htod_memcpy(dp, &init_val, sizeof(T));
        } else {
            dp = (T*)(std::malloc(sizeof(T)));
            *dp = init_val;
        }
    }

    ~DeviceScalar () {
        if (Gpu::inLaunchRegion()) {
            The_Arena()->free(dp);
        } else {
            std::free(dp);
        }
    }

    [[nodiscard]] T* dataPtr () { return dp; }
    [[nodiscard]] T const* dataPtr () const { return dp; }
    [[nodiscard]] T dataValue () const {
        if (Gpu::inLaunchRegion()) {
            T r;
            Gpu::dtoh_memcpy(&r, dp, sizeof(T));
            return r;
        } else {
            return *dp;
        }
    }

private:
    T* dp;

#else

    DeviceScalar (T const& init_val) : d(init_val) {}
    DeviceScalar () = default;
    ~DeviceScalar () = default;

    [[nodiscard]] T* dataPtr () { return &d; }
    [[nodiscard]] T const* dataPtr () const { return &d; }
    [[nodiscard]] T dataValue () const { return d; }

private:
    T d;

#endif
};

#ifdef AMREX_USE_GPU

template <class T>
struct SharedMemory
{
    [[nodiscard]] AMREX_GPU_DEVICE T* dataPtr () noexcept {
        static_assert(sizeof(T) < 0, "We must specialize struct SharedMemory");
        return nullptr;
    }
};

template <>
struct SharedMemory<double>
{
    [[nodiscard]] AMREX_GPU_DEVICE double* dataPtr () noexcept {
        extern __shared__ double amrex_sm_double[];
        return amrex_sm_double;
    }
};

template <>
struct SharedMemory<float>
{
    [[nodiscard]] AMREX_GPU_DEVICE float* dataPtr () noexcept {
        extern __shared__ float amrex_sm_float[];
        return amrex_sm_float;
    }
};

template <>
struct SharedMemory<long>
{
    [[nodiscard]] AMREX_GPU_DEVICE long* dataPtr () noexcept {
        extern __shared__ long amrex_sm_long[];
        return amrex_sm_long;
    }
};

template <>
struct SharedMemory<long long>
{
    [[nodiscard]] AMREX_GPU_DEVICE long long* dataPtr () noexcept {
        extern __shared__ long long amrex_sm_long_long[];
        return amrex_sm_long_long;
    }
};

template <>
struct SharedMemory<int>
{
    [[nodiscard]] AMREX_GPU_DEVICE int* dataPtr () noexcept {
        extern __shared__ int amrex_sm_int[];
        return amrex_sm_int;
    }
};

template <>
struct SharedMemory<short>
{
    [[nodiscard]] AMREX_GPU_DEVICE short* dataPtr () noexcept {
        extern __shared__ short amrex_sm_short[];
        return amrex_sm_short;
    }
};

template <>
struct SharedMemory<char>
{
    [[nodiscard]] AMREX_GPU_DEVICE char* dataPtr () noexcept {
        extern __shared__ char amrex_sm_char[];
        return amrex_sm_char;
    }
};

template <>
struct SharedMemory<unsigned long>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned long* dataPtr () noexcept {
        extern __shared__ unsigned long amrex_sm_ulong[];
        return amrex_sm_ulong;
    }
};

template <>
struct SharedMemory<unsigned long long>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned long long* dataPtr () noexcept {
        extern __shared__ unsigned long long amrex_sm_ulonglong[];
        return amrex_sm_ulonglong;
    }
};

template <>
struct SharedMemory<unsigned int>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned int* dataPtr () noexcept {
        extern __shared__ unsigned int amrex_sm_uint[];
        return amrex_sm_uint;
    }
};

template <>
struct SharedMemory<unsigned short>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned short* dataPtr () noexcept {
        extern __shared__ unsigned short amrex_sm_ushort[];
        return amrex_sm_ushort;
    }
};

template <>
struct SharedMemory<unsigned char>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned char* dataPtr () noexcept {
        extern __shared__ unsigned char amrex_sm_uchar[];
        return amrex_sm_uchar;
    }
};

template <>
struct SharedMemory<bool>
{
    [[nodiscard]] AMREX_GPU_DEVICE bool* dataPtr () noexcept {
        extern __shared__ bool amrex_sm_bool[];
        return amrex_sm_bool;
    }
};

#endif

}  //namespace


#endif
