#ifndef CHEM_TEST_H_
#define CHEM_TEST_H_

#include <iostream>
#include <memory>


#include <AMReX_AmrCore.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BLProfiler.H>

#include <mfix_eb.H>
#include <mfix_leveldata.H>
#include <mfix_fluid.H>
#include <mfix_reactions.H>
#include <mfix_ic.H>


#include <chem_test_prob.H>

class chem_test : public amrex::AmrCore
{
  public:

    ~chem_test ();

    chem_test ();

    void Init ( MFIXFluidPhase        const& a_fluid,
                MFIXReactions         const& a_reactions);

    amrex::Real thermo_p;
    amrex::Real thermo_po;

    void setup ( amrex::Vector< amrex::EBFArrayBoxFactory const*> a_factory,
                 MFIXFluidPhase        const& a_fluid,
                 MFIXReactions         const& a_reactions,
                 MFIXInitialConditions const& a_initial_conditions);

    void update_interp_data ( MFIXFluidPhase const& a_fluid );

    int nlev () const { return (finest_level+1); }

    MFIXEB& eb () noexcept { return m_eb; }

    MFIXLevelData& leveldata () noexcept { return m_level_data; }

    LevelData* leveldata ( int const a_lev ) noexcept {
      AMREX_ASSERT( a_lev < m_level_data.m_level_data.size() );
      return m_level_data.m_level_data[a_lev].get();
    }

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

    void MakeNewLevelFromScratch(int /*lev*/, amrex::Real /*ltime*/,
                                 const amrex::BoxArray& /*lnew_grids*/,
                                 const amrex::DistributionMapping& /*lnew_dmap*/) override
    { amrex::Abort("MakeNewLevelFromScratch: To be implemented"); }

  private:

    int const m_verbose;

    int const m_nghost;

    int m_newton_maxiter;
    amrex::Real m_newton_abstol;
    amrex::Real m_newton_reltol;

    DualGridAuxIndexes* m_p_indexes;

    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_interp_data;

  public:

    MFIXEB m_eb;

    MFIXLevelData m_level_data;

    int newton_maxiter () const { return m_newton_maxiter; }
    amrex::Real newton_abstol () const { return m_newton_abstol; }
    amrex::Real newton_reltol () const { return m_newton_reltol; }

    amrex::Vector< amrex::MultiFab* > interp_data () const noexcept
    { return GetVecOfPtrs(m_interp_data); }

    amrex::MultiFab* interp_data (int const a_lev) const noexcept
    { return m_interp_data[a_lev].get(); }

};


#endif
