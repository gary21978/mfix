#include <AMReX_ParmParse.H>

#include <dem_test_prob.H>

using namespace amrex;

dem_test_prob::
~dem_test_prob ()
{ }

dem_test_prob::
dem_test_prob (int const a_maxLevel, Vector<Geometry> a_geom)
  : m_bc_list(a_maxLevel + 1)
  , m_boundary_conditions(a_maxLevel, a_geom, {IntVect(1,1,1)}
  , m_bc_list, m_eb_boundary_conditions)
{
  m_timer.Initialize();

  m_species.Initialize();
  m_reactions.Initialize(m_species);
  m_fluid.Initialize(m_species, m_reactions);
  m_solids.Initialize(m_species, m_reactions);
  m_dem.Initialize();
  m_pic.Initialize();

  m_regions.Initialize();
  m_initial_conditions.Initialize(m_regions, m_fluid, m_solids, m_dem, m_pic);
  m_boundary_conditions.Initialize(m_regions, m_fluid, m_solids, m_dem, m_pic);

  ParmParse pp("mfix");
  pp.queryarr("gravity", m_gravity);
}
