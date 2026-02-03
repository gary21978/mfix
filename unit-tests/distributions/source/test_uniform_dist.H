#ifndef MFIX_TEST_UNIFORM_DIST_H_
#define MFIX_TEST_UNIFORM_DIST_H_

#include <string>

#include <AMReX_REAL.H>

#include <test_dist.H>

class test_uniform_dist : public test_dist
{

  public:

    test_uniform_dist () = default;
    ~test_uniform_dist () = default;

    int test ();
};

#endif
