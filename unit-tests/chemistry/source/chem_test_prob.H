#ifndef CHEM_TEST_PROB_H_
#define CHEM_TEST_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <mfix_pc.H>

#include <mfix_timer.H>
#include <mfix_species.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>
#include <mfix_reactions.H>

#include <mfix_regions.H>
#include <mfix_ic.H>
#include <mfix_bc.H>

//#include <mfix_eb.H>

class chem_test_prob {

  public:

    chem_test_prob (int const a_maxLevel,
                    amrex::Vector<amrex::Geometry> const& a_geom);

    ~chem_test_prob ();

    MFIXTimer& timer () noexcept { return m_timer; }
    MFIXTimer* p_timer () noexcept { return &m_timer; }

    MFIXSpecies& species () noexcept { return m_species; }
    MFIXSpecies* p_species () noexcept { return &m_species; }

    MFIXReactions& reactions () noexcept { return m_reactions; }
    MFIXReactions* p_reactions () noexcept { return &m_reactions; }

    MFIXFluidPhase& fluid ()  noexcept { return m_fluid; }
    MFIXFluidPhase* p_fluid ()  noexcept { return &m_fluid; }

    MFIXSolidsPhase& solids () noexcept { return m_solids; }
    MFIXSolidsPhase* p_solids () noexcept { return &m_solids; }

    MFIXDEM& dem () noexcept { return m_dem; }
    MFIXDEM* p_dem () noexcept { return &m_dem; }

    MFIXPIC& pic () noexcept { return m_pic; }
    MFIXPIC* p_pic () noexcept { return &m_pic; }

    MFIXInitialConditions& initial_conditions () noexcept { return m_initial_conditions; }
    MFIXInitialConditions* p_initial_conditions () noexcept { return &m_initial_conditions; }

    MFIXBoundaryConditions& boundary_conditions () noexcept { return m_boundary_conditions; }
    MFIXBoundaryConditions* p_boundary_conditions () noexcept { return &m_boundary_conditions; }

    BCList m_bc_list;

  private:

    MFIXTimer m_timer;

    MFIXSpecies m_species;
    MFIXReactions m_reactions;
    MFIXFluidPhase m_fluid;
    MFIXSolidsPhase m_solids;
    MFIXDEM m_dem;
    MFIXPIC m_pic;

    MFIXRegions m_regions;

    MFIXEmbeddedBoundaries m_eb_boundary_conditions;

    MFIXInitialConditions  m_initial_conditions;
    MFIXBoundaryConditions m_boundary_conditions;

};
#endif
