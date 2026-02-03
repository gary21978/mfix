#ifndef TEST_FLUID_EDDY_VISCOSITY_H_
#define TEST_FLUID_EDDY_VISCOSITY_H_

#include <AMReX_REAL.H>

class test_fluid_eddy_viscosity
{
  public:

    test_fluid_eddy_viscosity () = default;
    ~test_fluid_eddy_viscosity () = default;

    int wale ();

  private:

    int const PASS = 0;
    int const FAIL = 1;
    amrex::Real const MOL_VISC_CONST = 1.2e-5;

};
#endif
