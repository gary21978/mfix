#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Extension.H>

#include <mfix_timer.H>
#include <mfix_pc.H>

#include <dem_test.H>
#include <dem_test_eb.H>
#include <dem_test_prob.H>
#include <dem_test_check.H>

static int constexpr PASS(0);
static int constexpr FAIL(1);

using namespace amrex;

int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, true);
  {
    // initialize AmrCore
    dem_test test;
    test.InitFromScratch(0.);

    // initialize EB geometry
    dem_test_eb eb( test.Geom(), test.maxLevel() );
    eb.make_factory( test.Geom(), test.DistributionMap(), test.boxArray());

    // initialize problem setup (process inputs)
    dem_test_prob prob(test.maxLevel(), test.Geom());

    // Check reading of inputs
    if (!prob.dem().tan_history()) {
      Print() << "Tangential history not read correctly from input file" << std::endl;
      return FAIL;
    }
    if (prob.dem().tan_history_max_contacts() != 1) {
      Print() << "Tangential history max contacts not read correctly from input file" << std::endl;
      return FAIL;
    }

    // initialize particle container
    int include_virtual_mass = 0;
    MFIXParticleContainer pc( test.Geom(), test.DistributionMap(0),
        test.boxArray(0),
        prob.initial_conditions(),
        prob.boundary_conditions(),
        prob.solids(), prob.dem(), prob.pic(),
        prob.fluid(),  prob.reactions(), include_virtual_mass);

    LoadBalance loadbalance( /*max_level*/0, /*solve fluid*/0, /*solve solids*/1);
    loadbalance.ResetWeights(0, pc.ParticleBoxArray(0),
            pc.ParticleDistributionMap(0));

    // Set levelset values
    test.setup(eb.const_factory());

    // Create particles and populate container
    if ( prob.initial_conditions().AutoParticleInit() ) {
      Print() << "Auto generating particles ...\n";
      pc.InitParticlesAuto(eb.factory(0));
    } else {
#ifdef PARTICLE_INPUT_FILE
      std::string particle_input_file = AMREX_TO_STRING (PARTICLE_INPUT_FILE);
#else
      std::string particle_input_file = "particle_input.dat";
#endif
      Print() << "Reading particles from " + particle_input_file + " ...\n";
      pc.InitParticlesAscii(particle_input_file);
    }
    pc.Redistribute();
    pc.InitParticlesRuntimeVariables(prob.fluid().solve_enthalpy());

    // Initialize the collision parameters
    pc.MFIX_PC_InitCollisionParams();

    MFIXTimer& timer = prob.timer();

    std::string knapsack_weight_type = "RunTimeCosts";

    int steps = 0;
    do {

      timer.advance_nstep();
      timer.advance_time();

      for (int lev = 0; lev <= test.maxLevel(); ++lev) {
        int nsubsteps;
        const MultiFab* ls_data = test.m_level_sets[lev].get();

        pc.EvolveParticles(lev, timer.dt(), prob.gravity(), eb.factory(0),
            eb.factory(0), ls_data, test.m_levelset_refinement,
            &loadbalance, /*a_T_eb=*/nullptr, nsubsteps);
      }

      ++steps;

    } while (timer.time() < timer.stop_time() - 1e-8);

    dem_test_check check;
    if ( !check.compare(pc, test.maxLevel()+1) ) {
      Print() << "Comparison failed!\n";
      return FAIL;
    }
    if (!check.rolling_friction_inputs(prob.dem())) {
      Print() << "Rolling friction inputs failed!\n";
      return FAIL;
    }
  }

  // Make sure there isn't unused inputs
  if ( ParmParse::QueryUnusedInputs() ) {
    Print() << "Unused inputs detected. Test failed.\n";
    return FAIL;
  }

  amrex::Finalize();
}

