
#include <AMReX_GpuElixir.H>
#include <AMReX_GpuDevice.H>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>

namespace amrex::Gpu {

namespace {

#if defined(AMREX_USE_GPU)

extern "C" {
    void CUDART_CB amrex_elixir_delete (void* p)
    {
        auto p_pa = reinterpret_cast<Vector<std::pair<void*,Arena*> >*>(p);
        for (auto const& pa : *p_pa) {
            pa.second->free(pa.first);
        }
        delete p_pa;
    }
}

#endif

}

void
Elixir::clear () noexcept
{
#if defined(AMREX_USE_GPU)
    if (Gpu::inLaunchRegion())
    {
        if (!m_pa.empty()) {
#ifdef AMREX_USE_CUDA
            auto p = new Vector<std::pair<void*,Arena*> >(std::move(m_pa));
            AMREX_CUDA_SAFE_CALL(cudaLaunchHostFunc(Gpu::gpuStream(),
                                                    amrex_elixir_delete, (void*)p));
#endif
        }
    }
    else
#endif
    {
        for (auto const& pa : m_pa) {
            pa.second->free(pa.first);
        }
    }
    m_pa.clear();
}

}
