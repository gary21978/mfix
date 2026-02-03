#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include <test_fluid_incompressible.H>
#include <test_fluid_molecular_viscosity.H>
#include <test_fluid_suspension_viscosity.H>
#include <test_fluid_eddy_viscosity.H>
#include <test_fluid_porous_media.H>
#include <test_fluid_eddy_viscosity.H>

using namespace amrex;

int main(int argc, char* argv[])
{

  int constexpr PASS(0);
  int constexpr FAIL(1);

  bool constexpr build_parm_parse(true);

  amrex::Initialize(argc, argv, build_parm_parse);
  {

    { test_fluid_incompressible test;

      Print() << "test_fluid_incompressible::basic .................................. ";
      if (test.basic() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }
    }

    // Tests for molecular viscosity models
    { test_fluid_molecular_viscosity test;

      Print() << "test_fluid_molecular_viscosity::sutherland ........................ ";
      if (test.sutherland() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }

      Print() << "test_fluid_molecular_viscosity::reid .............................. ";
      if (test.reid() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }

      Print() << "test_fluid_molecular_viscosity::mixture ........................... ";
      if (test.mixture_sutherland() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }
    }

    // Tests for suspension viscosity models
    { test_fluid_suspension_viscosity test;

      Print() << "test_fluid_suspension_viscosity::einstein ......................... ";
      if (test.einstein() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }

      Print() << "test_fluid_suspension_viscosity::brinkman ......................... ";
      if (test.brinkman() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }

      Print() << "test_fluid_suspension_viscosity::roscoe ........................... ";
      if (test.roscoe() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }

      Print() << "test_fluid_suspension_viscosity::cheng_law ........................ ";
      if (test.cheng_law() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }

      Print() << "test_fluid_suspension_viscosity::sato       ....................... ";
      if (test.sato() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }
    }

    // Tests for eddy viscosity models
    { test_fluid_eddy_viscosity test;

      Print() << "test_fluid_eddy_viscosity::wale ................................... ";
      if (test.wale() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }
    }

    // Tests for porous media model
    { test_fluid_porous_media test;

      Print() << "test_fluid_porous_media::pipe ..................................... ";
      if (test.pipe_vfrac() == PASS) { Print() << "PASSED!\n"; }
      else { Print() << "FAILED!\n"; return FAIL;
      }
    }

  }
  amrex::Finalize();

  return PASS;
}



