#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>
#include <iostream>

namespace {
TEST_CASE("hull of translated, centered cube and a sphere", "[csg]") {
  auto st = *csg::parse_csg(R"(
hull() {
	multmatrix([[1, 0, 0, 5], [0, 1, 0, 6], [0, 0, 1, 7], [0, 0, 0, 1]]) {
		cube(size = [1, 2, 3], center = true);
	}
	sphere($fn = 0, $fa = 12, $fs = 2, r = 4);
}
)");
  auto hull = std::get<csg::Hull>(st.top.objs.back());
  auto [XX, YY, ZZ] = hull.cube_size;
  auto [Cx, Cy, Cz] = hull.cube_center;

  CHECK(XX == 1);
  CHECK(YY == 2);
  CHECK(ZZ == 3);

  CHECK(Cx == 5);
  CHECK(Cy == 6);
  CHECK(Cz == 7);

  CHECK(hull.sphere_radius == 4);
}

TEST_CASE("hull with non-centered cube not supported", "[csg]") {
  CHECK_THROWS_WITH(csg::parse_csg(R"(
hull() {
	multmatrix([[1, 0, 0, 5], [0, 1, 0, 6], [0, 0, 1, 7], [0, 0, 0, 1]]) {
		cube(size = [1, 2, 3], center = false);
	}
	sphere($fn = 0, $fa = 12, $fs = 2, r = 4);
}
)"),
                    Matchers::ContainsSubstring("action<hull>"));
}

TEST_CASE("hull with three elements not supported", "[csg]") {
  CHECK_THROWS_WITH(csg::parse_csg(R"(
hull() {
	multmatrix([[1, 0, 0, 5], [0, 1, 0, 6], [0, 0, 1, 7], [0, 0, 0, 1]]) {
		cube(size = [1, 2, 3], center = true);
	}
	sphere($fn = 0, $fa = 12, $fs = 2, r = 4);
	sphere($fn = 0, $fa = 12, $fs = 2, r = 4);
}
)"),
                    Matchers::ContainsSubstring("action<hull>"));
}

TEST_CASE("hull with other shape combos not supported", "[csg]") {
  CHECK_THROWS_WITH(csg::parse_csg(R"(
hull() {
   cube(size = [1, 2, 3], center = true);
	sphere($fn = 0, $fa = 12, $fs = 2, r = 4);
}
)"),
                    Matchers::ContainsSubstring("action<hull>"));
}

} // namespace
