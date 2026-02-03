#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

//  color() should be accepted, but ignored
TEST_CASE("color", "[csg]") {
  CHECK(csg::parse_csg(R"(
color([0, 1, 0, 1]) {
 cylinder($fn = 0, $fa = 5, $fs = 0.1, h = 20, r1 = 5, r2 = 5, center = true);
}
)"));
}

TEST_CASE("one shape in matmul", "[csg]") {
  auto st = *csg::parse_csg(R"(
multmatrix(
[
[1, 0, 0, 0.0020],
[0, 1, 0, 0.0005],
[0, 0, 1, 0.0005],
[0, 0, 0, 1]
]
) {
    cylinder(h = 2, r = 10, center=true);
}
)");
  auto mat = std::get<csg::Mulmatrix3D>(st.top.objs.back());
  CHECK(mat.group.objs.size() == 1);
  CHECK(mat.translation == std::array<double, 3>({0.0020, 0.0005, 0.0005}));
  CHECK(mat.rotation()[0] == std::array<double, 3>({1, 0, 0}));
  CHECK(mat.rotation()[1] == std::array<double, 3>({0, 1, 0}));
  CHECK(mat.rotation()[2] == std::array<double, 3>({0, 0, 1}));
  auto cyl = std::get<csg::Cylinder>(mat.group.objs.at(0));
  CHECK(cyl.height == 2);
  CHECK(cyl.radius == 10);
}

TEST_CASE("two shapes in matmul", "[csg]") {
  auto st = *csg::parse_csg(R"(
multmatrix(
[
[1, 0, 0, 0.0020],
[0, 1, 0, 0.0005],
[0, 0, 1, 0.0005],
[0, 0, 0, 1]
]
) {
    cylinder(h = 2, r = 10, center=true);
    sphere(r = 10);
}
)");
  auto mat = std::get<csg::Mulmatrix3D>(st.top.objs.back());
  CHECK(mat.translation == std::array<double, 3>({0.0020, 0.0005, 0.0005}));
  CHECK(mat.rotation()[0] == std::array<double, 3>({1, 0, 0}));
  CHECK(mat.rotation()[1] == std::array<double, 3>({0, 1, 0}));
  CHECK(mat.rotation()[2] == std::array<double, 3>({0, 0, 1}));
  auto cyl = std::get<csg::Cylinder>(mat.group.objs.at(0));
  auto sph = std::get<csg::Sphere>(mat.group.objs.at(1));
  CHECK(cyl.height == 2);
  CHECK(cyl.radius == 10);
  CHECK(sph.radius == 10);
}

TEST_CASE("two matmuls", "[csg]") {
  auto st = *csg::parse_csg(R"(
multmatrix(
[
[1, 0, 0, 0.0020],
[0, 1, 0, 0.0005],
[0, 0, 1, 0.0005],
[0, 0, 0, 1]
]
) {
    cylinder(h = 2, r = 10, center=true);
    sphere(r = 11);
}
multmatrix(
[
[1, 0, 0, 0.0020],
[0, 1, 0, 0.0005],
[0, 0, 1, 0.0005],
[0, 0, 0, 1]
]
) {
    cube(size = [1,2,3], center=true);
    cylinder(h=4, r1=1, r2=2, center=true);
}
)");
  auto mat = std::get<csg::Mulmatrix3D>(st.top.objs.at(0));
  CHECK(mat.translation == std::array<double, 3>({0.0020, 0.0005, 0.0005}));
  CHECK(mat.rotation()[0] == std::array<double, 3>({1, 0, 0}));
  CHECK(mat.rotation()[1] == std::array<double, 3>({0, 1, 0}));
  CHECK(mat.rotation()[2] == std::array<double, 3>({0, 0, 1}));
  auto cyl = std::get<csg::Cylinder>(mat.group.objs.at(0));
  auto sph = std::get<csg::Sphere>(mat.group.objs.at(1));
  CHECK(cyl.height == 2);
  CHECK(cyl.radius == 10);
  CHECK(sph.radius == 11);

  auto mat2 = std::get<csg::Mulmatrix3D>(st.top.objs.at(1));
  CHECK(mat2.translation == std::array<double, 3>({0.0020, 0.0005, 0.0005}));
  CHECK(mat2.rotation()[0] == std::array<double, 3>({1, 0, 0}));
  CHECK(mat2.rotation()[1] == std::array<double, 3>({0, 1, 0}));
  CHECK(mat2.rotation()[2] == std::array<double, 3>({0, 0, 1}));
  auto cube = std::get<csg::Cube>(mat2.group.objs.at(0));
  auto cone = std::get<csg::Cone>(mat2.group.objs.at(1));
  auto [XX, YY, ZZ] = cube.size;
  CHECK(XX == 1);
  CHECK(YY == 2);
  CHECK(ZZ == 3);
  CHECK(cone.height == 4);
  CHECK(cone.radius1 == 1);
  CHECK(cone.radius2 == 2);
}

TEST_CASE("two shapes one matmul", "[csg]") {
  auto st = *csg::parse_csg(R"(
cylinder(h = 2, r = 10, center=true);
sphere(r = 11);
multmatrix(
[
[1, 0, 0, 0.0020],
[0, 1, 0, 0.0005],
[0, 0, 1, 0.0005],
[0, 0, 0, 1]
]
) {
    cube(size = [1,2,3], center=true);
    cylinder(h=4, r1=1, r2=2, center=true);
}
)");
  auto cyl = std::get<csg::Cylinder>(st.top.objs.at(0));
  auto sph = std::get<csg::Sphere>(st.top.objs.at(1));
  auto mat = std::get<csg::Mulmatrix3D>(st.top.objs.at(2));
  CHECK(mat.translation == std::array<double, 3>({0.0020, 0.0005, 0.0005}));
  CHECK(mat.rotation()[0] == std::array<double, 3>({1, 0, 0}));
  CHECK(mat.rotation()[1] == std::array<double, 3>({0, 1, 0}));
  CHECK(mat.rotation()[2] == std::array<double, 3>({0, 0, 1}));
  CHECK(cyl.height == 2);
  CHECK(cyl.radius == 10);
  CHECK(sph.radius == 11);

  auto cube = std::get<csg::Cube>(mat.group.objs.at(0));
  auto cone = std::get<csg::Cone>(mat.group.objs.at(1));
  auto [XX, YY, ZZ] = cube.size;
  CHECK(XX == 1);
  CHECK(YY == 2);
  CHECK(ZZ == 3);
  CHECK(cone.height == 4);
  CHECK(cone.radius1 == 1);
  CHECK(cone.radius2 == 2);
}
