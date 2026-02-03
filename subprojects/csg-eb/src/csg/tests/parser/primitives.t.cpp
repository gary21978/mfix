#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

// Tests for the primitives on their own

namespace {
TEST_CASE("cylinder", "[csg]") {
  auto st = *csg::parse_csg(R"(
cylinder($name="my_cyl", h = 2, r = 10, center=true);
)");
  auto cyl = std::get<csg::Cylinder>(st.top.objs.at(0));
  CHECK(cyl.name == "my_cyl");
  CHECK(cyl.radius == 10);
  CHECK(cyl.height == 2);
}

TEST_CASE("nameless cylinder", "[csg]") {
  auto st = *csg::parse_csg(R"(
cylinder(h = 2, r = 10, center=true);
)");
  auto cyl = std::get<csg::Cylinder>(st.top.objs.at(0));
  CHECK_FALSE(cyl.name.has_value());
  CHECK(cyl.radius == 10);
  CHECK(cyl.height == 2);
}

TEST_CASE("cube", "[csg]") {
  auto st = *csg::parse_csg(R"(
cube(size = [1,2,3], center=true, $name="my_cube");
)");
  auto cub = std::get<csg::Cube>(st.top.objs.at(0));
  auto [XX, YY, ZZ] = cub.size;
  CHECK(cub.name == "my_cube");
  CHECK(XX == 1);
  CHECK(YY == 2);
  CHECK(ZZ == 3);
}

TEST_CASE("sphere", "[csg]") {
  auto st = *csg::parse_csg(R"(
sphere(r = 10, $name="my_sphere");
)");
  auto sph = std::get<csg::Sphere>(st.top.objs.at(0));
  CHECK(sph.name == "my_sphere");
  CHECK(sph.radius == 10);
}

TEST_CASE("semicolon optional", "[csg]") {
  auto st = *csg::parse_csg(R"(
sphere(r = 10)
)");
  auto sph = std::get<csg::Sphere>(st.top.objs.at(0));
  CHECK(sph.radius == 10);
}

#if USE_CGAL
TEST_CASE("square pyramid polyhedron", "[csg]") {
  auto st = *csg::parse_csg(R"(
polyhedron(
points = [[10, 10, 0], [10, -10, 0], [-10, -10, 0], [-10, 10, 0], [0, 0, 10]], 
faces = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4], [1, 0, 3], [2, 1, 3]], 
$name="my_polyhedron");
)");
  auto polyh = std::get<csg::Polyhedron>(st.top.objs.at(0));
  CHECK(polyh.name == "my_polyhedron");
  CHECK(polyh.cgal_polyhedron()->size_of_vertices() == 5);
  CHECK(polyh.cgal_polyhedron()->size_of_facets() == 6);
}
#endif

#if USE_CGAL
TEST_CASE("bad vs good polyhedron openscad example", "[csg]") {
  // https://en.wikibooks.org/wiki/OpenSCAD_User_Manual/Primitive_Solids#polyhedron
  CHECK_THROWS(csg::parse_csg(R"(
polyhedron(
points = [[0, -10, 60], [0, 10, 60], [0, 10, 0], 
[0, -10, 0], [60, -10, 60], [60, 10, 60], [10, -10, 50], 
[10, 10, 50], [10, 10, 30], [10, -10, 30], [30, -10, 50], 
[30, 10, 50]], 

faces = [[0, 2, 3], [0, 1, 2], [0, 4, 5], [0, 5, 1], 
[5, 4, 2], [2, 4, 3], [6, 8, 9], [6, 7, 8], [6, 10, 11], 
[6, 11, 7], [10, 8, 11], [10, 9, 8], [0, 3, 9], [9, 0, 6], 
[10, 6, 0], [0, 4, 10], [3, 9, 10], [3, 10, 4], [1, 7, 11], 
[1, 11, 5], [1, 7, 8], [1, 8, 2], [2, 8, 11], [2, 11, 5]], 

convexity = 1);
)"));

  auto st = *csg::parse_csg(R"(
polyhedron(
points = [[0, -10, 60], [0, 10, 60], [0, 10, 0], 
[0, -10, 0], [60, -10, 60], [60, 10, 60], [10, -10, 50], 
[10, 10, 50], [10, 10, 30], [10, -10, 30], [30, -10, 50], 
[30, 10, 50]], 

faces = [[0, 3, 2], [0, 2, 1], [4, 0, 5], [5, 0, 1], 
[5, 2, 4], [4, 2, 3], [6, 8, 9], [6, 7, 8], [6, 10, 11], 
[6, 11, 7], [10, 8, 11], [10, 9, 8], [3, 0, 9], [9, 0, 6], 
[10, 6, 0], [0, 4, 10], [3, 9, 10], [3, 10, 4], [1, 7, 11], 
[1, 11, 5], [1, 8, 7], [2, 8, 1], [8, 2, 11], [5, 11, 2]], 

convexity = 1);
)");
  auto polyh = std::get<csg::Polyhedron>(st.top.objs.at(0));
  CHECK(polyh.cgal_polyhedron()->size_of_vertices() == 12);
}
#endif

#if !(USE_CGAL)
TEST_CASE("without cgal test", "[csg]") {
  CHECK_THROWS(*csg::parse_csg(R"(
polyhedron(
points = [[10, 10, 0], [10, -10, 0], [-10, -10, 0], [-10, 10, 0], [0, 0, 10]], 
faces = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4], [1, 0, 3], [2, 1, 3]], 
$name="my_polyhedron");
)"));
}
#endif

} // namespace
