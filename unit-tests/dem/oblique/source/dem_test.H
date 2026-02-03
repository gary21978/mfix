#ifndef DEM_TEST_H_
#define DEM_TEST_H_

#include <iostream>
#include <memory>


#include <AMReX_AmrCore.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BLProfiler.H>

class dem_test : public amrex::AmrCore
{
  public:

    dem_test ();

    void InitFromScratch (amrex::Real a_time);

    void setup ( amrex::Vector< amrex::EBFArrayBoxFactory const*> a_factory);

    int m_levelset_refinement = 1;
    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_level_sets;

  protected:

    //! Tagging cells for refinement
    virtual void ErrorEst (int /*lev*/, amrex::TagBoxArray & /*tags*/,
                           amrex::Real /*time*/, int /*ngrow*/) override
    { amrex::Abort("ErrorEst: To be implemented"); }

    //! Make a new level using provided BoxArray and DistributionMapping and
    //! fill with interpolated coarse level data. Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int /*lev*/, amrex::Real /*time*/,
                                         const amrex::BoxArray & /*ba*/,
                                         const amrex::DistributionMapping & /*dm*/) override
    { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse data. Called
    //! by AmrCore::regrid.
    virtual void RemakeLevel (int /*lev*/, amrex::Real /*time*/,
                              const amrex::BoxArray & /*ba*/,
                              const amrex::DistributionMapping & /*dm*/) override
    { amrex::Abort("RemakeLevel: To be implemented"); }


    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int /*lev*/) override
    { amrex::Abort("ClearLevel: To be implemented"); }

    // Make a level 0 grids covering the whole domain.
    // It does NOT install the new grids.
    amrex::BoxArray MakeBaseGrids() const;

    void MakeNewLevelFromScratch(int lev, amrex::Real time,
                                 const amrex::BoxArray& new_grids,
                                 const amrex::DistributionMapping& new_dmap) override;

  private:

    int const m_verbose;
    int const m_nghost;
};

#endif
