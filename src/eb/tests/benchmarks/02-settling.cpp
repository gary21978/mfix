#include "catch2/catch_all.hpp"
using namespace Catch;

#include <AMReX_EB2_IF_Box.H>

#include "main.hpp"

namespace {

using namespace amrex;

void check(const csg::CsgIF &my_if, const EB2::BoxIF &their_box,
           std::tuple<double, double, double> point) {
  auto xx = std::get<0>(point);
  auto yy = std::get<1>(point);
  auto zz = std::get<2>(point);
  CHECK(Approx(my_if(xx, yy, zz)) == their_box(xx, yy, zz));
}

TEST_CASE("02-settling", "[benchmarks]") {
  EB2::BoxIF their_box({0, 0, 0}, {0.0015, 0.0015, 0.0015}, false);

  auto csg_str = std::string(R"(
   cube(size = [0.0015, 0.0015, 0.0015], center = false);
  )");
  auto maybe_tree = csg::parse_csg(csg_str);
  REQUIRE(maybe_tree);
  csg::CsgIF my_if(maybe_tree);

  SECTION("Outside") {
    check(my_if, their_box, {0.0016, 0, 0});
    check(my_if, their_box, {-0.001, 0, 0});
    check(my_if, their_box, {-0.001, 0.0017, 0});
    check(my_if, their_box, {-0.001, 0.0017, -0.0020});
  }
  SECTION("Inside") {
    check(my_if, their_box, {0.00075, 0, 0});
    check(my_if, their_box, {-0.0075, 0, 0});
    check(my_if, their_box, {0, 0.00075, 0});
    check(my_if, their_box, {0, 0.0008, 0.0008});
  }
}

} // namespace
