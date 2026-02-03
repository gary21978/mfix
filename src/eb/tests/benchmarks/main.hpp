#include <csg.hpp>

const int NN = 100;

// Compare a CSG IF with some other EB IF at "all" points within the region of
// interest
template <typename FF>
void test_for_equivalence(FF eb2_if, csg::CsgIF csg_if, double x_lo,
                          double x_hi, double y_lo, double y_hi, double z_lo,
                          double z_hi) {
  for (int ii = 0; ii < NN; ii++) {
    auto xx = ((NN - ii) * x_lo + ii * x_hi) / NN;
    for (int jj = 0; jj < NN; jj++) {
      auto yy = ((NN - jj) * y_lo + jj * y_hi) / NN;
      for (int kk = 0; kk < NN; kk++) {
        auto zz = ((NN - kk) * z_lo + kk * z_hi) / NN;

        REQUIRE(std::signbit(csg_if(xx, yy, zz)) ==
                std::signbit(eb2_if(xx, yy, zz)));
      }
    }
  }
}
