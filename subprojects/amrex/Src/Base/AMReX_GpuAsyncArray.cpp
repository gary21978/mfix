#include <AMReX_GpuAsyncArray.H>

#ifdef AMREX_USE_GPU

extern "C" {
    void CUDART_CB amrex_asyncarray_delete (void* p)
    {
        void** pp = (void**)p;
        void* dp = pp[0];
        void* hp = pp[1];
        std::free(p);
        amrex::The_Arena()->free(dp);
        amrex::The_Pinned_Arena()->free(hp);
    }
}

#endif
