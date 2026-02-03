#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

TEST_CASE("matmul in matmul", "[csg]") {
  auto maybe_st = csg::parse_csg(R"(
multmatrix(
[
[1, 0, 0, 1],
[0, 1, 0, 2],
[0, 0, 1, 3],
[0, 0, 0, 1]
]
) {
    cylinder(h = 2, r = 10, center=true);
multmatrix(
[
[1, 0, 0, 4],
[0, 1, 0, 5],
[0, 0, 1, 6],
[0, 0, 0, 1]
]
) {
    cylinder(h = 2, r = 10, center=true);
}}
)");
  CHECK(maybe_st != nullptr);
  auto st = maybe_st.get();
  auto mat = std::get<csg::Mulmatrix3D>(st->top.objs.back());
  auto mat2 = std::get<csg::Mulmatrix3D>(mat.group.objs.back());
  auto cyl = std::get<csg::Cylinder>(mat2.group.objs.at(0));
  CHECK(cyl.height == 2);
  CHECK(cyl.radius == 10);
}
