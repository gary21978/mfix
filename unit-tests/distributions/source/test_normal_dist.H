#ifndef MFIX_TEST_NORMAL_DIST_H_
#define MFIX_TEST_NORMAL_DIST_H_

#include <string>

#include <AMReX_REAL.H>

#include <test_dist.H>

class test_normal_dist : public test_dist
{

  public:

    test_normal_dist () = default;
    ~test_normal_dist () = default;

    int test ();
};

#endif
