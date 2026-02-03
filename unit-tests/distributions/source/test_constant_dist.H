#ifndef MFIX_TEST_CONSTANT_DIST_H_
#define MFIX_TEST_CONSTANT_DIST_H_

#include <string>

#include <AMReX_REAL.H>

#include <test_dist.H>

class test_constant_dist : public test_dist
{

  public:

    test_constant_dist () = default;
    ~test_constant_dist () = default;

    int test ();
};

#endif
