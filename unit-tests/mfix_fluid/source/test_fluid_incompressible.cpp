#include <AMReX_ParmParse.H>

#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_fluid.H>

#include <test_fluid_incompressible.H>

using namespace amrex;

///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                    Test basic incompressible fluid                  ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_incompressible::
basic ()
{
  ParmParse pp;

  MFIXSpecies species;
  MFIXReactions rxns;

  // Test disabling fluid solve by setting name to none.
  { pp.add("fluid.solve", std::string("None"));

    MFIXFluidPhase fluid;
    if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }
    if ( fluid.solve() ) { return FAIL; }
  }

  // Test disabling fluid solve by not defining it
  { pp.remove("fluid.solve");

    MFIXFluidPhase fluid;
    if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }
    if ( fluid.solve() ) { return FAIL; }
  }

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   0);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXFluidPhase fluid;

  // Fails because species and rxns are not initialized
  if (fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  species.Initialize();

  // Fails because rxns is not initialized
  if (fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  rxns.Initialize(species);

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("constant"));

  // Fails because molecular viscosity constant is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular.constant", 1.8e-3);

  // Basic incompressible fluid is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }


  // cleanup added pp entries.
  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");

  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.constant");

  return PASS;
}
