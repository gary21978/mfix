#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

// Tests for CSG union(), intersection(), and difference()

TEST_CASE("two shapes", "[csg]") {
  auto st = *csg::parse_csg(R"(
sphere(r = 10);
cube(size = [1,2,3], center=true);
)");
  auto sph = std::get<csg::Sphere>(st.top.objs.at(0));
  CHECK(sph.radius == 10);
  auto cub = std::get<csg::Cube>(st.top.objs.at(1));

  auto [XX, YY, ZZ] = cub.size;
  CHECK(XX == 1);
  CHECK(YY == 2);
  CHECK(ZZ == 3);
}

TEST_CASE("union", "[csg]") {
  auto st = *csg::parse_csg(R"(
union() {
  cube(size = [12, 14, 15], center = true);
  sphere(r = 8);
}
)");

  auto un = std::get<csg::Union3D>(st.top.objs.back());
  CHECK(st.top.objs.size() == 1);

  auto cub = std::get<csg::Cube>(un.objs.at(0));
  auto sph = std::get<csg::Sphere>(un.objs.at(1));

  auto [XX, YY, ZZ] = cub.size;
  CHECK(XX == 12);
  CHECK(YY == 14);
  CHECK(ZZ == 15);
  CHECK(sph.radius == 8);
}

TEST_CASE("intersection", "[csg]") {
  auto st = *csg::parse_csg(R"(
intersection() {
  cube(size = [15, 15, 15], center = true);
  sphere(r = 10);
}
)");

  CHECK(st.top.objs.size() == 1);
  auto ints = std::get<csg::Intersection3D>(st.top.objs.back());
  CHECK(ints.objs.size() == 2);

  auto cub = std::get<csg::Cube>(ints.objs.at(0));

  auto sph = std::get<csg::Sphere>(ints.objs.at(1));

  auto [XX, YY, ZZ] = cub.size;
  CHECK(XX == 15);
  CHECK(YY == 15);
  CHECK(ZZ == 15);
  CHECK(sph.radius == 10);
}

TEST_CASE("difference", "[csg]") {
  auto st = *csg::parse_csg(R"(
difference() {
  cube(size = [12, 12, 12], center = true);
  sphere(r = 8);
}
)");

  auto diff = std::get<csg::Difference3D>(st.top.objs.back());
  CHECK(st.top.objs.size() == 1);

  auto cub = std::get<csg::Cube>(*diff.first_obj);
  auto sph = std::get<csg::Sphere>(diff.next_objs.objs.at(0));

  auto [XX, YY, ZZ] = cub.size;
  CHECK(XX == 12);
  CHECK(YY == 12);
  CHECK(ZZ == 12);
  CHECK(sph.radius == 8);
}

TEST_CASE("incompatible 2D shape", "[csg]") {
  auto st = csg::parse_csg(R"(
sphere(r = 10);
circle(size = [1,2], center=true);
)");
  CHECK(st == nullptr);
}

TEST_CASE("sphere and mulmat cylinder", "[csg]") {
  auto st = *csg::parse_csg(R"(
sphere(r = 0.1);
multmatrix([[1, 0, 0, 0.5], [0, 1, 0, 0.1], [0, 0, 1, 0], [0, 0, 0, 1]]) {
   cylinder(h = 0.2, r1 = 0.3, r2 = 0.4, center = false);
}
)");

  CHECK(st.top.objs.size() == 2);

  auto sph = std::get<csg::Sphere>(st.top.objs.at(0));
  auto mm = std::get<csg::Mulmatrix3D>(st.top.objs.at(1));

  CHECK(sph.radius == 0.1);
  CHECK(mm.group.objs.size() == 1);
  CHECK(mm.translation == std::array<double, 3>({0.5, 0.1, 0}));
  CHECK(mm.rotation()[0] == std::array<double, 3>({1, 0, 0}));
  CHECK(mm.rotation()[1] == std::array<double, 3>({0, 1, 0}));
  CHECK(mm.rotation()[2] == std::array<double, 3>({0, 0, 1}));

  auto cone = std::get<csg::Cone>(mm.group.objs.at(0));
  CHECK(cone.height == 0.2);
  CHECK(cone.radius1 == 0.3);
  CHECK(cone.radius2 == 0.4);
  CHECK(cone.center == false);
}
