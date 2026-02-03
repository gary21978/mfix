#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EB2_IF.H>

#include <dem_test_eb.H>

using namespace amrex;

dem_test_eb::
dem_test_eb ( Vector<Geometry> const& a_geom, int const a_max_lev)
  : m_support_level(EBSupport::full)
{
  std::string type;

  ParmParse pp;
  if (pp.query("mfix.geometry", type)) {
    if (type == "box") {
      Print() << "\n Building box geometry.\n";
      make_box(a_geom, a_max_lev);
    } else {
      Abort("Unknown geometry type specified");
    }
  } else {
    EB2::AllRegularIF regular;
    auto gshop = EB2::makeShop(regular);
    build_level(gshop, a_geom, a_max_lev);
  }
}

void dem_test_eb::make_box (Vector<Geometry> const& geom,
                            int const max_level) {
  // Get box information from inputs file
  ParmParse pp("box");

  if(geom[0].isAllPeriodic())
  {
    EB2::AllRegularIF regular;
    auto gshop = EB2::makeShop(regular);

    build_level(gshop, geom, max_level);
  }
  else
  {
    /************************************************************************
     *                                                                      *
     * Define Box geometry:                                                 *
     *        -> box.{Lo,Hi} vector storing box lo/hi                       *
     *        -> box.offset  vector storing box offset                      *
     * NOTE: walls are placed _outside_ domain for periodic directions.     *
     *                                                                      *
     ************************************************************************/

    Vector<Real> boxLo(3), boxHi(3);
    Real offset = 1.0e-15;
    bool inside = true;

    for(int i = 0; i < 3; i++)
    {
        boxLo[i] = geom[0].ProbLo(i);
        boxHi[i] = geom[0].ProbHi(i);
    }

    pp.queryarr("Lo", boxLo, 0, 3);
    pp.queryarr("Hi", boxHi, 0, 3);

    pp.query("offset", offset);
    pp.query("internal_flow", inside);

    Real xlo = boxLo[0] + offset;
    Real xhi = boxHi[0] - offset;

    // This ensures that the walls won't even touch the ghost cells. By
    // putting them one domain width away
    if(geom[0].isPeriodic(0))
    {
        xlo = 2.0 * geom[0].ProbLo(0) - geom[0].ProbHi(0);
        xhi = 2.0 * geom[0].ProbHi(0) - geom[0].ProbLo(0);
    }

    Real ylo = boxLo[1] + offset;
    Real yhi = boxHi[1] - offset;

    // This ensures that the walls won't even touch the ghost cells. By
    // putting them one domain width away
    if(geom[0].isPeriodic(1))
    {
        ylo = 2.0 * geom[0].ProbLo(1) - geom[0].ProbHi(1);
        yhi = 2.0 * geom[0].ProbHi(1) - geom[0].ProbLo(1);
    }

    Real zlo = boxLo[2] + offset;
    Real zhi = boxHi[2] - offset;

    // This ensures that the walls won't even touch the ghost cells. By
    // putting them one domain width away
    if(geom[0].isPeriodic(2))
    {
        zlo = 2.0 * geom[0].ProbLo(2) - geom[0].ProbHi(2);
        zhi = 2.0 * geom[0].ProbHi(2) - geom[0].ProbLo(2);
    }

    RealArray lo {xlo, ylo, zlo};
    RealArray hi {xhi, yhi, zhi};

    EB2::BoxIF my_box(lo, hi, inside);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(my_box);

    build_level(gshop, geom, max_level);
  }
}

void
dem_test_eb::
make_factory ( Vector<Geometry>            const& a_geom,
               Vector<DistributionMapping> const& a_dmap,
               Vector<BoxArray>            const& a_grids)
{
    int const nlev(a_geom.size());

  // Build the eb factory which has the cut-cell information
  m_factory.resize(nlev);

  for (int lev(0); lev<nlev; ++lev) {
    m_factory[lev] = std::make_unique<EBFArrayBoxFactory>(*m_eb_levels[lev],
      a_geom[lev], a_grids[lev], a_dmap[lev],
      Vector<int>{nghost_basic(),nghost_volume(),nghost_full()}, m_support_level);
  }
}

