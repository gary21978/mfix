#ifndef AMREX_OPENMP_H_
#define AMREX_OPENMP_H_
#include <AMReX_Config.H>


namespace amrex::OpenMP {

    constexpr int get_num_threads () { return 1; }
    constexpr int get_max_threads () { return 1; }
    constexpr int get_thread_num  () { return 0; }
    constexpr int in_parallel     () { return false; }
    constexpr void set_num_threads (int) { /* nothing */ }
}


namespace amrex {
    /** ... */
    int
    numUniquePhysicalCores();
}

#endif
