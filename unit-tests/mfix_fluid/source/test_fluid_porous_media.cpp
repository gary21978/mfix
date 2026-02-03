#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <mfix_regions.H>
#include <mfix_porous_media.H>

#include <test_fluid_porous_media.H>

using namespace amrex;

int test_fluid_porous_media::
pipe_vfrac ()
{
  ParmParse pp;

  int nlev = 1;

  pp.addarr("mfix.regions", std::vector<std::string>{"bed"});
  pp.addarr("regions.bed.lo", std::vector<double>{0.0030,0.0000,0.0000});
  pp.addarr("regions.bed.hi", std::vector<double>{0.0064,0.0020,0.0020});
  pp.addarr("pm.regions", std::vector<std::string>{"bed"});
  pp.add("pm.bed.volfrac", 0.2);
  pp.add("pm.bed.c1", 1.e-5);
  pp.add("pm.bed.c2", 0.);
  pp.add("pm.bed.allow_particles", 0);

  pp.add("fluid.solve", std::string("fluid0"));

  Vector<int> n_cell{40,10,10};
  IntVect grid_sizes(n_cell);
  IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
  hi -= IntVect::TheUnitVector();
  Box index_domain(lo,hi);

  Geometry a_geom;
  a_geom.define(index_domain, RealBox({0.0,0.0,0.0},{0.008,0.002,0.002}),
      0, {0,1,1});

  Vector<Geometry> geom;
  geom.push_back(a_geom);

  MFIXSpecies species;
  MFIXReactions rxns;

  MFIXFluidPhase fluid;

  species.Initialize();
  rxns.Initialize(species);

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("constant"));
  pp.add("fluid0.viscosity.molecular.constant", 2.0e-5);

  // Fluid fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  MFIXRegions regions;
  regions.Initialize();

  MFIXPorousMedia pm;
  pm.Initialize(regions, geom);

  BoxArray ba(geom[0].Domain());
  ba.maxSize(grid_sizes);

  DistributionMapping dm;
  dm.define(ba, ParallelDescriptor::NProcs());

  Vector<MultiFab*> ep_g, level_set;
  ep_g.resize(nlev);
  level_set.resize(nlev);

  for (int lev = 0; lev < nlev; ++lev) {
    ep_g[lev] = new MultiFab(ba, dm, 1, 0, MFInfo());
    ep_g[lev]->setVal(1.);
    const auto nd_ba = amrex::convert(ba, IntVect::TheNodeVector());
    level_set[lev] = new MultiFab(nd_ba, dm, 1, 0, MFInfo());
    level_set[lev]->setVal(1.);
  }

  {
    // Check the volume fraction
    pm.impose_vfrac(ep_g);

    ReduceOps<ReduceOpLogicalOr> reduce_op;
    ReduceData<int> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(*ep_g[0]); mfi.isValid(); ++mfi) {
      const auto& bx = mfi.validbox();
      Array4<Real const> const& ep_g_arr = ep_g[0]->array(mfi);

      reduce_op.eval(bx, reduce_data, [lPASS=PASS, lFAIL=FAIL, ep_g_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        if ((i >= 15) && (i <= 31)) {
          return (ep_g_arr(i,j,k) != 0.8) ? lFAIL : lPASS;
        }

        return (ep_g_arr(i,j,k) == 1.0) ? lPASS : lFAIL;
      });
    }
    ReduceTuple host_tuple = reduce_data.value(reduce_op);
    if (get<0>(host_tuple) == FAIL) { return FAIL; }
  }

  {
    // Check the level set
    // Note that fluid has positive sign
    pm.block_particles(*level_set[0], geom[0]);

    ReduceOps<ReduceOpLogicalOr> reduce_op;
    ReduceData<int> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(*level_set[0]); mfi.isValid(); ++mfi) {
      const auto& bx = mfi.validbox();
      Array4<Real const> const& levset_arr = level_set[0]->array(mfi);

      reduce_op.eval(bx, reduce_data, [lPASS=PASS, lFAIL=FAIL, levset_arr]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        if ((i >= 15) && (i <= 32)) {
          return (levset_arr(i,j,k) > 0.) ? lFAIL : lPASS;
        }

        return (levset_arr(i,j,k) >= 0.) ? lPASS : lFAIL;
      });
    }
    ReduceTuple host_tuple = reduce_data.value(reduce_op);
    if (get<0>(host_tuple) == FAIL) { return FAIL; }
  }

  // Clear multifab
  for (int lev = 0; lev < nlev; ++lev)
  {
    delete ep_g[lev];
    delete level_set[lev];
  }

  pp.remove("mfix.regions");
  pp.remove("regions.bed.lo");
  pp.remove("regions.bed.hi");
  pp.remove("pm.regions");
  pp.remove("pm.bed.volfrac");
  pp.remove("pm.bed.c1");
  pp.remove("pm.bed.c2");
  pp.remove("pm.bed.allow_particles");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.constant");

  return PASS;
}
