#ifndef MFIX_MFHELPERS_H
#define MFIX_MFHELPERS_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>

namespace MFHelpers
{

    // Creates an exact copy ("shape" + value) of "mold"
    std::unique_ptr<amrex::MultiFab> createFrom (amrex::MultiFab& mold);

    // Creates a copy of "mold" and initializes it to "val".
    std::unique_ptr<amrex::MultiFab> createFrom (amrex::MultiFab& mold,
                                                 amrex::Real val);

    // Creates a copy of "mold" with nGrow ghosts and initializes it to "val"
    std::unique_ptr<amrex::MultiFab> createFrom (amrex::MultiFab& mold,
                                                 amrex::Real val,
                                                 int nGrow);

    // Creates a copy of "mold" with nComp components, nGrow ghosts and
    // initializes it to "val"
    std::unique_ptr<amrex::MultiFab> createFrom (amrex::MultiFab& mold,
                                                 amrex::Real val,
                                                 int nGrow,
                                                 int nComp);
}

#endif
