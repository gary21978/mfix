#include <mfix_fluid_update.H>

using namespace amrex;

FluidUpdate::
FluidUpdate ( int const a_nlev,
              StepType const a_step_type,
              Vector<Geometry> const& a_geom,
              Vector<BoxArray> const& a_ba,
              const MFIXFluidPhase& a_fluid,
              MFIXLevelData& a_leveldata,
              MFIXStepData& a_stepdata,
              MFIXBoundaryConditions& a_bcs,
              int const a_include_chem,
              Real const a_time, Real const a_dt)
  : m_nlev(a_nlev)
  , m_step_type(a_step_type)
  , m_geom(a_geom)
  , m_ba(a_ba)
  , m_fluid(a_fluid)
  , m_leveldata(a_leveldata)
  , m_stepdata(a_stepdata)
  , m_bcs(a_bcs)
  , m_include_chem(a_include_chem)
  , m_time(a_time)
  , m_dt(a_dt)
  , m_chem_idxs(a_fluid.nspecies(), a_include_chem)
{
}
