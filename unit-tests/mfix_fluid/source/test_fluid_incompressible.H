#ifndef TEST_FLUID_INCOMPRESSIBLE_H_
#define TEST_FLUID_INCOMPRESSIBLE_H_

class test_fluid_incompressible
{
  public:

    test_fluid_incompressible () = default;
    ~test_fluid_incompressible () = default;

    int basic ();

  private:

    int const PASS = 0;
    int const FAIL = 1;

};
#endif
