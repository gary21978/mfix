#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

// Tests that don't fit elsewhere

//  group() same as union()

TEST_CASE("group", "[csg]") {
  CHECK(csg::parse_csg(R"(
group() {
   cylinder($fn = 0, $fa = 5, $fs = 0.1, h = 20, r1 = 5, r2 = 5, center = true);
}
)"));
}

//  render() should be accepted, but ignored
TEST_CASE("render", "[csg]") {
  CHECK(csg::parse_csg(R"(
render(convexity=2) {
   cylinder($fn = 0, $fa = 5, $fs = 0.1, h = 20, r1 = 5, r2 = 5, center = true);
}
)"));
}

//  First Basic Example when running OpenSCAD

TEST_CASE("OpenSCAD Basic Example", "[csg]") {
  CHECK(csg::parse_csg(R"(
multmatrix([[1, 0, 0, -24], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) {
  union() {
    cube(size = [15, 15, 15], center = true);
    sphere($fn = 0, $fa = 12, $fs = 2, r = 10);
  }
}
intersection() {
  cube(size = [15, 15, 15], center = true);
  sphere($fn = 0, $fa = 12, $fs = 2, r = 10);
}
multmatrix([[1, 0, 0, 24], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) {
  difference() {
    cube(size = [15, 15, 15], center = true);
    sphere(r = 10);
  }
}
group();
)"));
}
