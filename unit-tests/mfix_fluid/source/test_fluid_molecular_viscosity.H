#ifndef TEST_FLUID_MOLECULAR_VISCOSITY_H_
#define TEST_FLUID_MOLECULAR_VISCOSITY_H_

class test_fluid_molecular_viscosity
{
  public:

    test_fluid_molecular_viscosity () = default;
    ~test_fluid_molecular_viscosity () = default;

    int sutherland ();
    int reid ();
    int mixture_sutherland ();

  private:

    int const PASS = 0;
    int const FAIL = 1;

};
#endif
