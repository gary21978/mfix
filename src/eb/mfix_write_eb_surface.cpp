#include <mfix_eb.H>

#include <AMReX_WriteEBSurface.H>

using namespace amrex;

void MFIXEB::write_surface (Vector<Geometry> const& a_geom,
                            Vector<DistributionMapping> const& a_dmap,
                            Vector<BoxArray> const& a_grids ) const
{
  if ((!m_write_surface) || (a_geom[0].isAllPeriodic()))
  { return; }

  // Only write at the finest level!
  int lev = m_max_level;

  WriteEBSurface(a_grids[lev], a_dmap[lev],
      a_geom[lev], m_ebfactory[lev].get());
}
