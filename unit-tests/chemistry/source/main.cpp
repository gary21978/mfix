#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <mfix_timer.H>
#include <mfix_pc.H>

#include <mfix_stepdata.H>
#include <mfix_fluid_update.H>
#include <mfix_cmi.H>

#include <chem_test.H>
#include <chem_test_prob.H>
#include <chem_test_check.H>

static int constexpr PASS(0);
static int constexpr FAIL(1);

using namespace amrex;


int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, true);
  {
    // initialize AmrCore
    chem_test test;

    // initialize problem setup (process inputs)
    chem_test_prob prob(test.maxLevel(), test.Geom());

    // Create the box array and distribution map
    test.Init( prob.fluid(), prob.reactions());

    int const include_virtual_mass(0);

    // initialize particle container
    MFIXParticleContainer particles( test.Geom(), test.DistributionMap(0),
        test.boxArray(0),
        prob.initial_conditions(),
        prob.boundary_conditions(),
        prob.solids(), prob.dem(), prob.pic(),
        prob.fluid(),  prob.reactions(), include_virtual_mass);

    // Set initial values for field and particle data
    test.setup( test.eb().factory_const(), prob.fluid(), prob.reactions(),
                 prob.initial_conditions());

    prob.m_bc_list.MakeBCArrays(0, 2, test.Geom(0), false);
    prob.boundary_conditions().set_bc_list(0, 2);

    // Create particles and populate container
    if ( prob.dem().solve() || prob.pic().solve() ) {
      if ( prob.initial_conditions().AutoParticleInit() ) {
        Print() << "Auto generating particles ...\n";
        particles.InitParticlesAuto(test.eb().factory()[0].get());
      } else {
        Print() << "Reading particles from particle_input.dat ...\n";
        particles.InitParticlesAscii("particle_input.dat");
      }
      particles.Redistribute();
      particles.InitParticlesRuntimeVariables(prob.fluid().solve_enthalpy());
    }

    MFIXTimer& timer = prob.timer();

    int const verbose(0);

    Real const* const dx = test.Geom(0).CellSize();
    Real const vol( dx[0]*dx[1]*dx[2] );

    MFIXFluidPhase const& fluid = prob.fluid();

    chem_test_check check( fluid.solve_species(),
        fluid.solve_enthalpy(), fluid.nspecies(),
        fluid.isMixture(), prob.dem().solve() || prob.pic().solve(), vol,
        fluid.props.data<run_on>(), prob.solids().props.data<run_on>());

    // Copy level data into interp data
    test.update_interp_data( fluid );
    check.write(0., test.interp_data(), &particles);

    MFIXStepData stepData(test.nlev(), test.boxArray(), test.DistributionMap(),
        test.eb().factory_const(), 1, 1, prob.boundary_conditions() );

    stepData.define( 0., fluid.solve_density(), fluid.solve_enthalpy(),
        fluid.nspecies(), 0, 0, 0);

    FluidUpdate updateFluid( /*nlev=*/1, StepType::Predictor, test.Geom(), test.boxArray(),
        prob.fluid(), test.leveldata(), stepData, prob.boundary_conditions(),
        prob.reactions().solve(), 0.0, timer.dt() );

    updateFluid.setExplicitDiffusion();

    do {

      timer.advance_nstep();
      timer.advance_time();

      Real const time = timer.time();

      ChemistryManagementInterface cmi;

      const int report_mass_balance = 0;

      cmi.update_pointers(prob.p_reactions(), prob.p_solids(), prob.p_fluid(),
        prob.p_dem(), prob.p_pic(), &particles, prob.p_timer(), report_mass_balance);

      cmi.calc_chemistry(test.Geom(), test.leveldata().chem_txfr(), test.interp_data(),
        test.thermo_p, time, timer.dt(), verbose);

      updateFluid.Density( );

      updateFluid.Species( /*diffOpSpecies=*/nullptr );

      updateFluid.Energy( test.thermo_po, test.thermo_p, /*diffOpEnergy*/nullptr,
          test.newton_abstol(), test.newton_reltol(), test.newton_maxiter());

      // Copy level data into interp data
      test.update_interp_data( fluid );

      check.write(time, test.interp_data(), &particles);

#if defined(CHEM_TEST_EULERIAN01)
      check.eulerian01(time, test.interp_data());
#endif

    } while (timer.ok());

    ParmParse pp("test");

    Real abs_tol(0.);
    pp.query("abs_tol", abs_tol);

    if (check.L1norm() < abs_tol ) { Print() << "PASSED!\n"; }
    //else { Print() << "FAILED!\n"; return FAIL; }

    int print_L1norm(0);
    pp.query("print_L1norm", print_L1norm);
    if (print_L1norm) { printf("L1 Norm: %21.18e\n", check.L1norm()); }

    // Make sure there isn't unused inputs
    if ( ParmParse::QueryUnusedInputs() ) {
      Print() << "Unused inputs detected. Test failed.\n";
      return FAIL;
    }
  }

  amrex::Finalize();
}
