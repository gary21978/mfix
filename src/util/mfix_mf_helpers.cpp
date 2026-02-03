#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <mfix_mf_helpers.H>

using namespace amrex;


//
// Creates an exact copy ("shape" + value) of "mold"
//
std::unique_ptr< MultiFab >
MFHelpers::createFrom (MultiFab& mold)
{
    AMREX_ASSERT(mold.ok());

    std::unique_ptr< MultiFab > mf;

    mf = std::make_unique<MultiFab>(mold.boxArray(), mold.DistributionMap(),
                          mold.nComp(), mold.nGrow(), MFInfo(),
                          mold.Factory() );

    MultiFab::Copy(*mf, mold, 0, 0, mold.nComp(), mold.nGrow());

    return mf;
}


//
// Creates a copy of "mold" and initializes it to "val".
//
std::unique_ptr< MultiFab >
MFHelpers::createFrom ( MultiFab& mold, Real val )
{
    AMREX_ASSERT(mold.ok());

    std::unique_ptr< MultiFab > mf;

    mf = std::make_unique<MultiFab>(mold.boxArray(), mold.DistributionMap(),
                          mold.nComp(), mold.nGrow(), MFInfo(),
                          mold.Factory() );
    mf -> setVal(val);

    return mf;
}

//
// Creates a copy of "mold" with nGrow ghosts and initializes it to "val".
//
std::unique_ptr< MultiFab >
MFHelpers::createFrom (MultiFab& mold, Real val, int nGrow )
{
    AMREX_ASSERT(mold.ok());

    std::unique_ptr< MultiFab > mf;

    mf = std::make_unique<MultiFab>(mold.boxArray(), mold.DistributionMap(),
                          mold.nComp(), nGrow, MFInfo(),
                          mold.Factory() );
    mf -> setVal(val);

    return mf;

}


//
// Creates a copy of "mold" with nComp components, nGrow ghosts and initializes
// it to "val".
//
std::unique_ptr< MultiFab >
MFHelpers::createFrom (MultiFab& mold, Real val, int nGrow, int nComp)
{
    AMREX_ASSERT(mold.ok());

    std::unique_ptr< MultiFab > mf;

    mf = std::make_unique<MultiFab>(mold.boxArray(), mold.DistributionMap(),
                          nComp, nGrow, MFInfo(), mold.Factory());

    mf->setVal(val);

    return mf;
}
