#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_fluid.H>
#include <mfix_run_on.H>

#include <test_fluid_eddy_viscosity.H>

using namespace amrex;

///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                            WALE                                     ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_eddy_viscosity::
wale ()
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

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("conSTANT"));

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }

  pp.add("fluid0.viscosity.molecular.constant", MOL_VISC_CONST);

  // Constant molecular viscosity is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Define wale eddy viscosity model
  pp.add("fluid0.viscosity.eddy", std::string("waLE"));

  // model is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Sanity checks
  if (!fluid.isInitialized() ) { return FAIL; }
  if ( fluid.eddy_visc_model() != EddyViscModel::WALE ) { return FAIL; }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.constant");
  pp.remove("fluid0.viscosity.eddy");

  return PASS;
}
