#include <AMReX_BaseFab.H>
#include <AMReX_BLFort.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <cstring>
#include <cstdlib>

namespace amrex {

std::atomic<Long> atomic_total_bytes_allocated_in_fabs     {0L};
std::atomic<Long> atomic_total_bytes_allocated_in_fabs_hwm {0L};
std::atomic<Long> atomic_total_cells_allocated_in_fabs     {0L};
std::atomic<Long> atomic_total_cells_allocated_in_fabs_hwm {0L};
Long private_total_bytes_allocated_in_fabs     = 0L;
Long private_total_bytes_allocated_in_fabs_hwm = 0L;
Long private_total_cells_allocated_in_fabs     = 0L;
Long private_total_cells_allocated_in_fabs_hwm = 0L;

namespace
{
    bool basefab_initialized = false;
}

void
BaseFab_Initialize ()
{
    if (!basefab_initialized)
    {
        basefab_initialized = true;


#ifdef AMREX_MEM_PROFILING
        MemProfiler::add("Fab", std::function<MemProfiler::MemInfo()>
                         ([] () -> MemProfiler::MemInfo {
                             return {amrex::TotalBytesAllocatedInFabs(),
                                     amrex::TotalBytesAllocatedInFabsHWM()};
                         }));
#endif
    }

    amrex::ExecOnFinalize(amrex::BaseFab_Finalize);
}

void
BaseFab_Finalize()
{
    basefab_initialized = false;
}


Long
TotalBytesAllocatedInFabs () noexcept
{
    return private_total_bytes_allocated_in_fabs
        + atomic_total_bytes_allocated_in_fabs.load(std::memory_order_relaxed);
}

Long
TotalBytesAllocatedInFabsHWM () noexcept
{
    return private_total_bytes_allocated_in_fabs_hwm
        + atomic_total_bytes_allocated_in_fabs_hwm.load(std::memory_order_relaxed);
}

Long
TotalCellsAllocatedInFabs () noexcept
{
    return private_total_cells_allocated_in_fabs
        + atomic_total_cells_allocated_in_fabs.load(std::memory_order_relaxed);
}

Long
TotalCellsAllocatedInFabsHWM () noexcept
{
    return private_total_cells_allocated_in_fabs_hwm
        + atomic_total_cells_allocated_in_fabs_hwm.load(std::memory_order_relaxed);
}

void
ResetTotalBytesAllocatedInFabsHWM () noexcept
{
    {
        private_total_bytes_allocated_in_fabs_hwm = 0;
    }
    atomic_total_bytes_allocated_in_fabs_hwm.store(0,std::memory_order_relaxed);
}

void
update_fab_stats (Long n, Long s, size_t szt) noexcept
{
    {
        Long tst = s * static_cast<Long>(szt);
        Long old_bytes = amrex::atomic_total_bytes_allocated_in_fabs.fetch_add
            (tst,std::memory_order_relaxed);
        Long new_bytes = old_bytes + tst;
        Long prev_bytes_hwm = amrex::atomic_total_bytes_allocated_in_fabs_hwm.load
            (std::memory_order_relaxed);
        while (prev_bytes_hwm < new_bytes) {
            if (amrex::atomic_total_bytes_allocated_in_fabs_hwm.compare_exchange_weak
                (prev_bytes_hwm, new_bytes, std::memory_order_release, std::memory_order_relaxed)) {
                break;
            }
        }

        if(szt == sizeof(Real)) {
            Long old_cells = amrex::atomic_total_cells_allocated_in_fabs.fetch_add
                (n,std::memory_order_relaxed);
            Long new_cells = old_cells + n;
            Long prev_cells_hwm = amrex::atomic_total_cells_allocated_in_fabs_hwm.load
                (std::memory_order_relaxed);
            while (prev_cells_hwm < new_cells) {
                if (amrex::atomic_total_cells_allocated_in_fabs_hwm.compare_exchange_weak
                    (prev_cells_hwm, new_cells, std::memory_order_release, std::memory_order_relaxed)) {
                    break;
                }
            }
        }
    }
}

}
