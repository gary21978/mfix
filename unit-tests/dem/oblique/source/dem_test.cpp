#include <AMReX_EB2.H>
#include <AMReX_EB_utils.H>

#include <dem_test.H>

using namespace amrex;

dem_test::
dem_test ( )
  : m_verbose(0)
  , m_nghost(0)
{
  // Do not to iterate when creating the initial grid hierarchy
  SetIterateToFalse();

  // Use the new chopping routine which rejects cuts if
  // they don't improve the efficiency
  SetUseNewChop();
}

BoxArray
dem_test::
MakeBaseGrids () const
{
  if (m_verbose > 1) {
    amrex::Print() << "MakeBaseGrids\n";
  }
  BoxArray ba(geom[0].Domain());

  ba.maxSize(max_grid_size[0]);

  // to avoid duplicates
  if (ba == grids[0]) { ba = grids[0]; }

  if (m_verbose) {
    amrex::Print() << "Initial BoxArray has " << ba.size() << " grids\n";
  }

  return ba;
}

void
dem_test::
MakeNewLevelFromScratch (int a_lev, Real /*a_time*/,
                         const BoxArray& a_new_grids,
                         const DistributionMapping& a_new_dmap)

{
  if (m_verbose) {
    amrex::Print() << "MakeNewLevelFromScratch" << std::endl;

    if (a_lev == 0) {
      std::cout << "Making level 0 with box array\n" << a_new_grids << "\n";
    } else {
      Print() << "Setting refined region at level " << a_lev
              << " to " << a_new_grids << "\n";
    }
  }

  SetBoxArray(a_lev, a_new_grids);
  SetDistributionMap(a_lev, a_new_dmap);

  if (m_verbose > 1) {
    Print() << "SETTING NEW GRIDS IN MAKE NEW LEVEL " << a_new_grids << std::endl;
    Print() << "SETTING NEW DMAP IN MAKE NEW LEVEL " << a_new_dmap << std::endl;
  }

}


void
dem_test::
InitFromScratch (Real a_time)
{
  //! EB level constructed from building GeometryShop

  const BoxArray& ba = MakeBaseGrids();

  DistributionMapping dm(ba, ParallelDescriptor::NProcs());

  MakeNewLevelFromScratch(0, a_time, ba, dm);

  for (int lev(1); lev < maxLevel()+1; lev++) {
     MakeNewLevelFromScratch(lev, a_time, grids[lev], dmap[lev]);
  }
}

void
dem_test::
setup ( Vector< EBFArrayBoxFactory const*> a_factory)
{
  //! Generate level sets
  m_level_sets.resize(2);
  int levelset_pad = 2;

  BoxArray ls_ba = convert(boxArray(0), IntVect::TheNodeVector());
  m_level_sets[0] = std::make_unique<MultiFab>(ls_ba, DistributionMap(0), 1,
      levelset_pad/m_levelset_refinement);
  m_level_sets[1] = std::make_unique<MultiFab>(ls_ba, DistributionMap(0), 1, levelset_pad);

  if (maxLevel() == 0) {
    amrex::FillSignedDistance(*m_level_sets[1], *a_factory[0]->getEBLevel(), *a_factory[0],
                              m_levelset_refinement);
    MultiFab::Copy(*m_level_sets[0], *m_level_sets[1], 0, 0, 1, m_level_sets[0]->nGrow());
  }

}

