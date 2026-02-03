#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>

#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_fluid.H>
#include <mfix_run_on.H>
#include <mfix_diffusion_op.H>
#include <mfix_solids.H>
#include <mfix_avg_pc_parms.H>

#include <test_fluid_suspension_viscosity.H>

#include <memory>

using namespace amrex;

test_fluid_suspension_viscosity::test_fluid_suspension_viscosity ()
  : vel_fab(nullptr)
  , epf_fab(nullptr)
  , rho_fab(nullptr)
  , Tf_fab(nullptr)
  , avg_part_fab(nullptr)
{
  Vector<int> n_cell{4,4,4};
  IntVect grid_sizes(n_cell);
  IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
  hi -= IntVect::TheUnitVector();
  Box index_domain(lo,hi);
  geom.define(index_domain, RealBox({0.0,0.0,0.0},{4.0,4.0,4.0}), 0, {0,0,0});

  ba.define(geom.Domain());
  ba.maxSize(grid_sizes);
  dm.define(ba, ParallelDescriptor::NProcs());

  EB2::AllRegularIF regular;
  Build(EB2::makeShop(regular), geom, 0, 100);
  ebfact = std::make_unique<EBFArrayBoxFactory>(EB2::IndexSpace::top().getLevel(geom),
      geom, ba, dm, Vector<int>{2,2,2}, EBSupport::full);

  vel_fab = new MultiFab(ba, dm, 3, 1, MFInfo(), *ebfact);
  epf_fab = new MultiFab(ba, dm, 1, 1, MFInfo(), *ebfact);
  rho_fab = new MultiFab(ba, dm, 1, 1, MFInfo(), *ebfact);
  Tf_fab = new MultiFab(ba, dm, 1, 1, MFInfo(), *ebfact);
  avg_part_fab = new MultiFab(ba, dm, 4, 1, MFInfo(), *ebfact);
}

test_fluid_suspension_viscosity::~test_fluid_suspension_viscosity()
{
  delete vel_fab;
  delete epf_fab;
  delete rho_fab;
  delete Tf_fab;
  delete avg_part_fab;
}

void test_fluid_suspension_viscosity::init_fabs (const Box& bx,
                                                 Array4<Real>& epg,
                                                 Array4<Real>&rho)
{
  Real const epg_lo(0.35), epg_hi(1.0);
  Real const d_epg( (epg_hi - epg_lo) / static_cast<Real>(bx.numPts()-1));

  Real const rho_lo(1.2), rho_hi(12.0);
  Real const d_rho( (rho_hi - rho_lo) / static_cast<Real>(bx.numPts()-1));

  Real const l0 = bx.length(0);
  Real const l1 = bx.length(1);

  ParallelFor(bx, [l0,l1,epg_lo,d_epg,rho_lo,d_rho,epg,rho]
  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    epg(i,j,k) = epg_lo + d_epg*static_cast<Real>(i)
                        + l0*d_epg*static_cast<Real>(j)
                        + (l0*l1)*d_epg*static_cast<Real>(k);

    rho(i,j,k) = rho_lo + d_rho*static_cast<Real>(i)
                        + l0*d_rho*static_cast<Real>(j)
                        + (l0*l1)*d_rho*static_cast<Real>(k);

  });
}

///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                            Einstein                                 ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_suspension_viscosity::
einstein ()
{
  ParmParse pp;

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   0);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXSpecies species;
  MFIXReactions rxns;

  MFIXFluidPhase fluid;

  species.Initialize();
  rxns.Initialize(species);

  MFIXAvgParticleParms avg_particle_parms;
  MFIXEmbeddedBoundaries eb_parms;

  avg_particle_parms.Initialize( SoArealData::count, fluid, eb_parms );

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("conSTANT"));

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }

  pp.add("fluid0.viscosity.molecular.constant", MOL_VISC_CONST);

  // Constant molecular viscosity is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Define einstein suspension viscosity model
  pp.add("fluid0.viscosity.suspension", std::string("eINStein"));

  pp.add("fluid0.viscosity.max_effective_factor", MAX_EFF_VISC_FAC);
  // model is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Sanity checks
  if (!fluid.isInitialized() ) { return FAIL; }
  if ( fluid.susp_visc_model() != SuspViscModel::Einstein ) { return FAIL; }
  if (avg_particle_parms.nComp() != 0) { return FAIL; }

  const auto props = fluid.props.data<run_on>();

  {
    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }
      Array4<Real> epf = epf_fab->array(mfi);
      Array4<Real> rho = rho_fab->array(mfi);

      init_fabs(ba[0], epf, rho);
    }

    BCList bc_list(1);
    MFIXEmbeddedBoundaries mfix_eb;
    Vector<Geometry> gm = {geom};
    MFIXBoundaryConditions mfix_bcs(1, gm, {IntVect(1,1,1)}, bc_list, mfix_eb);
    Vector<BCRec> bcrec(1);

    auto diffOpVel = MFIXDiffOpVelocity(1, {geom}, {ba}, {dm}, {ebfact.get()},
        fluid, avg_particle_parms, mfix_bcs, bcrec);

    diffOpVel.setDiffCoeff({vel_fab}, {epf_fab}, {rho_fab}, {Tf_fab}, {avg_part_fab}, {nullptr}, 1);

    Vector<MultiFab const*> eff_visc_fabs = diffOpVel.getDiffCoeff();

    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }

      const auto& eff_visc = eff_visc_fabs[0]->array(mfi);
      const auto& epf = epf_fab->array(mfi);
      const auto& rho = rho_fab->array(mfi);

      ReduceOps<ReduceOpLogicalOr> reduce_op;
      ReduceData<int> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      reduce_op.eval(ba[0], reduce_data, [lPASS=PASS, lFAIL=FAIL, eff_visc,
                                          epf, rho, mol_visc=MOL_VISC_CONST,
                                          dummy=DUMMY, props]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        /* Compare against the actual Einstein expression given by
            mu_susp = mu_mol * 2.5 * (1 - epf)
        */
        Real visc0 = mol_visc * 2.5*(1. - epf(i,j,k));

        Real mvisc = props.molViscosity(dummy, IntVect(), Array4<Real const>());
        Real eff_visc0 = visc0 + mvisc;

        Real err = std::abs(eff_visc0 - eff_visc(i,j,k));

        return (err <= 1.e-14) ? lPASS : lFAIL;
      });

      ReduceTuple host_tuple = reduce_data.value(reduce_op);
      if (get<0>(host_tuple) == FAIL) { return FAIL; }
    }
  }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.constant");
  pp.remove("fluid0.viscosity.suspension");
  pp.remove("fluid0.viscosity.max_effective_factor");

  return PASS;
}


///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                            Brinkman                                 ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_suspension_viscosity::
brinkman ()
{
  ParmParse pp;

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   0);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXSpecies species;
  MFIXReactions rxns;

  MFIXFluidPhase fluid;

  species.Initialize();
  rxns.Initialize(species);

  MFIXAvgParticleParms avg_particle_parms;
  MFIXEmbeddedBoundaries eb_parms;

  avg_particle_parms.Initialize( SoArealData::count, fluid, eb_parms );

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("conSTANT"));

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }

  pp.add("fluid0.viscosity.molecular.constant", MOL_VISC_CONST);

  // Constant molecular viscosity is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Define Brinkman suspension viscosity model
  pp.add("fluid0.viscosity.suspension", std::string("BrINKman"));

  // Brinkman viscosity model required inputs:
  //   "viscosity.suspension.Brinkman.constant"

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.suspension.Brinkman.constant",BRINKMAN_CONST);

  pp.add("fluid0.viscosity.max_effective_factor", MAX_EFF_VISC_FAC);
  // model is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Sanity checks
  if (!fluid.isInitialized() ) { return FAIL; }
  if ( fluid.susp_visc_model() != SuspViscModel::Brinkman ) { return FAIL; }
  if (avg_particle_parms.nComp() != 0) { return FAIL; }

  const auto props = fluid.props.data<run_on>();

  {
    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }
      Array4<Real> epf = epf_fab->array(mfi);
      Array4<Real> rho = rho_fab->array(mfi);

      init_fabs(ba[0], epf, rho);
    }

    BCList bc_list(1);
    MFIXEmbeddedBoundaries mfix_eb;
    Vector<Geometry> gm = {geom};
    MFIXBoundaryConditions mfix_bcs(1, gm, {IntVect(1,1,1)}, bc_list, mfix_eb);
    Vector<BCRec> bcrec(1);

    auto diffOpVel = MFIXDiffOpVelocity(1, {geom}, {ba}, {dm}, {ebfact.get()},
        fluid, avg_particle_parms, mfix_bcs, bcrec);

    diffOpVel.setDiffCoeff({vel_fab}, {epf_fab}, {rho_fab}, {Tf_fab}, {avg_part_fab}, {nullptr}, 1);

    Vector<MultiFab const*> eff_visc_fabs = diffOpVel.getDiffCoeff();

    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }

      const auto& eff_visc = eff_visc_fabs[0]->array(mfi);
      const auto& epf = epf_fab->array(mfi);
      const auto& rho = rho_fab->array(mfi);

      ReduceOps<ReduceOpLogicalOr> reduce_op;
      ReduceData<int> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      reduce_op.eval(ba[0], reduce_data, [lPASS=PASS, lFAIL=FAIL, eff_visc, rho,
                                          epf, props, mol_visc=MOL_VISC_CONST,
                                          dummy=DUMMY, brinkmanConst=BRINKMAN_CONST]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        /* Compare against the actual Brinkman expression given by
            mu_susp = mu_mol * (epf^-c - 1)
        */
        Real visc0 = mol_visc * (std::pow(epf(i,j,k), -brinkmanConst) - 1.);

        Real mvisc = props.molViscosity(dummy, IntVect(), Array4<Real const>());
        Real eff_visc0 = visc0 + mvisc;

        Real eff_visc1 = eff_visc(i,j,k);

        Real err = std::abs(eff_visc0 - eff_visc1);

        return (err <= 1.e-14) ? lPASS : lFAIL;
      });

      ReduceTuple host_tuple = reduce_data.value(reduce_op);
      if (get<0>(host_tuple) == FAIL) { return FAIL; }
    }
  }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.constant");
  pp.remove("fluid0.viscosity.suspension");
  pp.remove("fluid0.viscosity.suspension.Brinkman.constant");
  pp.remove("fluid0.viscosity.max_effective_factor");

  return PASS;
}

///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                            Roscoe                                   ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_suspension_viscosity::
roscoe ()
{
  ParmParse pp;

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   0);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXSpecies species;
  MFIXReactions rxns;

  MFIXFluidPhase fluid;

  species.Initialize();
  rxns.Initialize(species);

  MFIXAvgParticleParms avg_particle_parms;
  MFIXEmbeddedBoundaries eb_parms;

  avg_particle_parms.Initialize( SoArealData::count, fluid, eb_parms );

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("conSTANT"));

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }

  pp.add("fluid0.viscosity.molecular.constant", MOL_VISC_CONST);

  // Constant molecular viscosity is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Define Roscoe suspension viscosity model
  pp.add("fluid0.viscosity.suspension", std::string("rosCOE"));

  // Brinkman viscosity model required inputs:
  //   "viscosity.suspension.Roscoe.c1"
  //   "viscosity.suspension.Roscoe.c2"

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.suspension.Roscoe.c1", ROSCOE_C1);

  // Fails because one more constant is not specified
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.suspension.Roscoe.c2", ROSCOE_C2);

  pp.add("fluid0.viscosity.max_effective_factor", MAX_EFF_VISC_FAC);
  // model is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Sanity checks
  if (!fluid.isInitialized() ) { return FAIL; }
  if ( fluid.susp_visc_model() != SuspViscModel::Roscoe ) { return FAIL; }
  if (avg_particle_parms.nComp() != 0) { return FAIL; }

  const auto props = fluid.props.data<run_on>();

  {
    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }
      Array4<Real> epf = epf_fab->array(mfi);
      Array4<Real> rho = rho_fab->array(mfi);

      init_fabs(ba[0], epf, rho);
    }

    BCList bc_list(1);
    MFIXEmbeddedBoundaries mfix_eb;
    Vector<Geometry> gm = {geom};
    MFIXBoundaryConditions mfix_bcs(1, gm, {IntVect(1,1,1)}, bc_list, mfix_eb);
    Vector<BCRec> bcrec(1);

    auto diffOpVel = MFIXDiffOpVelocity(1, {geom}, {ba}, {dm}, {ebfact.get()},
        fluid, avg_particle_parms, mfix_bcs, bcrec);

    diffOpVel.setDiffCoeff({vel_fab}, {epf_fab}, {rho_fab}, {Tf_fab}, {avg_part_fab}, {nullptr}, 1);

    Vector<MultiFab const*> eff_visc_fabs = diffOpVel.getDiffCoeff();

    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }

      const auto& eff_visc = eff_visc_fabs[0]->array(mfi);
      const auto& epf = epf_fab->array(mfi);
      const auto& rho = rho_fab->array(mfi);

      ReduceOps<ReduceOpLogicalOr> reduce_op;
      ReduceData<int> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;
      reduce_op.eval(ba[0], reduce_data, [lPASS=PASS, lFAIL=FAIL, eff_visc, rho,
                                          epf, props, mol_visc=MOL_VISC_CONST,
                                          dummy=DUMMY, c1=ROSCOE_C1,
                                          c2=ROSCOE_C2, fac=MAX_EFF_VISC_FAC]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        /* Compare against the actual Roscoe expression given by
            mu_susp = mu_mol * [(1-eps/c1)^(-c2)-1]
        */
        Real eps = 1. - epf(i,j,k);
        Real visc0 = 0.;
        if (std::abs(eps-c1) > 1.e-9) {
          visc0 = mol_visc * (std::pow((1. - (eps/c1)), -c2) - 1.);
        } else {
          // Special handling when eps is very close to c1
          visc0 = mol_visc * (std::pow(1.e-9/c1, -c2) - 1.);
        }

        Real mvisc = props.molViscosity(dummy,IntVect(), Array4<Real const>());
        Real eff_visc0 = min(visc0 + mvisc, fac*mol_visc);

        Real eff_visc1 = eff_visc(i,j,k);

        Real err = std::abs(eff_visc0 - eff_visc1);

        return (err <= 1.e-14) ? lPASS : lFAIL;
      });

      ReduceTuple host_tuple = reduce_data.value(reduce_op);
      if (get<0>(host_tuple) == FAIL) { return FAIL; }
    }
  }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.constant");
  pp.remove("fluid0.viscosity.suspension");
  pp.remove("fluid0.viscosity.suspension.Roscoe.c1");
  pp.remove("fluid0.viscosity.suspension.Roscoe.c2");
  pp.remove("fluid0.viscosity.max_effective_factor");

  return PASS;
}

///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                            ChengLaw                                 ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_suspension_viscosity::
cheng_law ()
{
  ParmParse pp;

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   0);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXSpecies species;
  MFIXReactions rxns;

  MFIXFluidPhase fluid;

  species.Initialize();
  rxns.Initialize(species);

  MFIXAvgParticleParms avg_particle_parms;
  MFIXEmbeddedBoundaries eb_parms;

  avg_particle_parms.Initialize( SoArealData::count, fluid, eb_parms );

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("conSTANT"));

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }

  pp.add("fluid0.viscosity.molecular.constant", MOL_VISC_CONST);

  // Constant molecular viscosity is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Define ChengLaw suspension viscosity model
  pp.add("fluid0.viscosity.suspension", std::string("chenGlaW"));

  // Brinkman viscosity model required inputs:
  //   "viscosity.suspension.ChengLaw.constant"

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.suspension.ChengLaw.constant", CHENG_LAW_CONST);

  pp.add("fluid0.viscosity.max_effective_factor", MAX_EFF_VISC_FAC);
  // model is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Sanity checks
  if (!fluid.isInitialized() ) { return FAIL; }
  if ( fluid.susp_visc_model() != SuspViscModel::ChengLaw ) { return FAIL; }
  if (avg_particle_parms.nComp() != 0) { return FAIL; }

  const auto props = fluid.props.data<run_on>();

  {
    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }
      Array4<Real> epf = epf_fab->array(mfi);
      Array4<Real> rho = rho_fab->array(mfi);

      init_fabs(ba[0], epf, rho);
    }

    BCList bc_list(1);
    MFIXEmbeddedBoundaries mfix_eb;
    Vector<Geometry> gm = {geom};
    MFIXBoundaryConditions mfix_bcs(1, gm, {IntVect(1,1,1)}, bc_list, mfix_eb);
    Vector<BCRec> bcrec(1);

    auto diffOpVel = MFIXDiffOpVelocity(1, {geom}, {ba}, {dm}, {ebfact.get()},
        fluid, avg_particle_parms, mfix_bcs, bcrec);

    diffOpVel.setDiffCoeff({vel_fab}, {epf_fab}, {rho_fab}, {Tf_fab}, {avg_part_fab}, {nullptr}, 1);

    Vector<MultiFab const*> eff_visc_fabs = diffOpVel.getDiffCoeff();

    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }

      const auto& eff_visc = eff_visc_fabs[0]->array(mfi);
      const auto& epf = epf_fab->array(mfi);
      const auto& rho = rho_fab->array(mfi);

      ReduceOps<ReduceOpLogicalOr> reduce_op;
      ReduceData<int> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;
      reduce_op.eval(ba[0], reduce_data, [lPASS=PASS, lFAIL=FAIL, eff_visc, rho,
                                       epf, props, mol_visc=MOL_VISC_CONST,
                                       dummy=DUMMY, c=CHENG_LAW_CONST, fac=MAX_EFF_VISC_FAC]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        /* Compare against the actual ChengLaw expression given by
            mu_susp = mu_mol * {exp[2.5*(1/((1-eps)^c)-1)/c] - 1}
        */
        Real visc0 = mol_visc * (std::exp(2.5*(std::pow(epf(i,j,k), -c) - 1.)/c) - 1.);

        Real mvisc = props.molViscosity(dummy, IntVect(), Array4<Real const>());
        Real eff_visc0 = min(visc0 + mvisc, fac*mol_visc);

        Real eff_visc1 = eff_visc(i,j,k);

        Real err = std::abs(eff_visc0 - eff_visc1);

        return (err <= 1.e-14) ? lPASS : lFAIL;
      });

      ReduceTuple host_tuple = reduce_data.value(reduce_op);
      if (get<0>(host_tuple) == FAIL) { return FAIL; }
    }
  }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.constant");
  pp.remove("fluid0.viscosity.suspension");
  pp.remove("fluid0.viscosity.suspension.ChengLaw.constant");
  pp.remove("fluid0.viscosity.max_effective_factor");

  return PASS;
}

///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                            Sato                                     ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_suspension_viscosity::
sato ()
{
  ParmParse pp;

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   0);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXSpecies species;
  MFIXReactions rxns;

  MFIXFluidPhase fluid;

  species.Initialize();
  rxns.Initialize(species);

  MFIXAvgParticleParms avg_particle_parms;
  MFIXEmbeddedBoundaries eb_parms;

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("conSTANT"));

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }

  pp.add("fluid0.viscosity.molecular.constant", MOL_VISC_CONST);

  // Constant molecular viscosity is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Define Sato suspension viscosity model
  pp.add("fluid0.viscosity.suspension", std::string("saTO"));
  pp.add("fluid0.viscosity.suspension.Sato.constant", SATO_CONST);
  pp.add("fluid0.viscosity.max_effective_factor", MAX_EFF_VISC_FAC);

  // model is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Sanity checks
  if (!fluid.isInitialized() ) { return FAIL; }
  if ( fluid.susp_visc_model() != SuspViscModel::Sato ) { return FAIL; }

  avg_particle_parms.Initialize( SoArealData::count, fluid, eb_parms );
  if (avg_particle_parms.nComp() == 0) { return FAIL; }

  // Set deposition components
  int const rad_comp  = avg_particle_parms.comp(SoArealData::radius);
  if ( rad_comp  < 0 ) { return FAIL; }

  int const velx_comp = avg_particle_parms.comp(SoArealData::velx);
  if ( velx_comp < 0 ) { return FAIL; }

  int const vely_comp = avg_particle_parms.comp(SoArealData::vely);
  if ( vely_comp < 0 ) { return FAIL; }

  int const velz_comp = avg_particle_parms.comp(SoArealData::velz);
  if ( velz_comp < 0 ) { return FAIL; }

  const auto props = fluid.props.data<run_on>();

  {
    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }
      Array4<Real> epf = epf_fab->array(mfi);
      Array4<Real> rho = rho_fab->array(mfi);
      Array4<Real> vel = vel_fab->array(mfi);
      Array4<Real> rad_velp = avg_part_fab->array(mfi);

      init_fabs(ba[0], epf, rho);

      ParallelFor(ba[0], [vel, rad_velp, rad_comp,
      velx_comp, vely_comp, velz_comp]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        vel(i,j,k,0) = 2.;
        vel(i,j,k,1) = 0.;
        vel(i,j,k,2) = 0.;

        rad_velp(i,j,k,rad_comp ) = 0.001;
        rad_velp(i,j,k,velx_comp) = 1.;
        rad_velp(i,j,k,vely_comp) = 0.;
        rad_velp(i,j,k,velz_comp) = 0.;
      });
    }


    BCList bc_list(1);
    MFIXEmbeddedBoundaries mfix_eb;
    Vector<Geometry> gm = {geom};
    MFIXBoundaryConditions mfix_bcs(1, gm, {IntVect(1,1,1)}, bc_list, mfix_eb);
    Vector<BCRec> bcrec(1);

    auto diffOpVel = MFIXDiffOpVelocity(1, {geom}, {ba}, {dm}, {ebfact.get()},
        fluid, avg_particle_parms, mfix_bcs, bcrec);

    diffOpVel.setDiffCoeff({vel_fab}, {epf_fab}, {rho_fab}, {Tf_fab}, {avg_part_fab}, {nullptr}, 1);

    Vector<MultiFab const*> eff_visc_fabs = diffOpVel.getDiffCoeff();

    {
      MFIter mfi(*epf_fab);
      if (!mfi.isValid()) { return FAIL; }

      const auto& eff_visc = eff_visc_fabs[0]->array(mfi);
      const auto& epf = epf_fab->array(mfi);
      const auto& rho = rho_fab->array(mfi);
      const auto& vel = vel_fab->array(mfi);
      const auto& rad_velp = avg_part_fab->array(mfi);

      ReduceOps<ReduceOpLogicalOr> reduce_op;
      ReduceData<int> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;
      reduce_op.eval(ba[0], reduce_data, [lPASS=PASS, lFAIL=FAIL, eff_visc, rho,
                                          epf, props, mol_visc=MOL_VISC_CONST,
                                          dummy=DUMMY, fac=MAX_EFF_VISC_FAC,
                                          vel, rad_velp, satoConst=SATO_CONST]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        /* Compare against the actual Sato expression given by
            mu_susp = CONST * d_p * (1 - epf) * | u_g - u_p | * rho
        */
        Real visc0 = satoConst * 2.* 0.001 * (1. - epf(i,j,k)) * 1. * rho(i,j,k);

        Real mvisc = props.molViscosity(dummy, IntVect(), Array4<Real const>());
        Real eff_visc0 = min(visc0 + mvisc, fac*mol_visc);

        Real eff_visc1 = eff_visc(i,j,k);

        Real err = std::abs(eff_visc0 - eff_visc1);

        return (err <= 1.e-14) ? lPASS : lFAIL;
      });

      ReduceTuple host_tuple = reduce_data.value(reduce_op);
      if (get<0>(host_tuple) == FAIL) { return FAIL; }
    }
  }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.constant");
  pp.remove("fluid0.viscosity.suspension");
  pp.remove("fluid0.viscosity.suspension.Sato.constant");
  pp.remove("fluid0.viscosity.max_effective_factor");
  return PASS;
}
