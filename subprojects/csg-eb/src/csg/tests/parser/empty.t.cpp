#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

// Tests for CSG union(), intersection(), and difference()

namespace {
TEST_CASE("mulmat and empty group", "[csg]") {
  auto st = *csg::parse_csg(R"(
multmatrix([[0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1]]) {
   cylinder($fn = 0, $fa = 12, $fs = 2, h = 200, r1 = 7, r2 = 7, center = false);
}
group();
)");
  auto un = std::get<csg::Union3D>(st.top.objs.at(1));

  CHECK(un.objs.size() == 0);
}

TEST_CASE("difference and empty group", "[csg]") {
  auto st = *csg::parse_csg(R"(
difference() {
	sphere(r = 2);
	sphere(r = 1);
}
group();
)");
  auto diff = std::get<csg::Difference3D>(st.top.objs.at(0));
  auto sph1 = std::get<csg::Sphere>(*diff.first_obj);
  auto sph2 = std::get<csg::Sphere>(diff.next_objs.objs.at(0));

  auto un = std::get<csg::Union3D>(st.top.objs.at(1));

  CHECK(sph1.radius == 2);
  CHECK(sph2.radius == 1);
  CHECK(un.objs.size() == 0);
}

TEST_CASE("intersection and empty group", "[csg]") {
  auto st = *csg::parse_csg(R"(
intersection() {
	sphere(r = 2);
	sphere(r = 3);
}
group();
)");
  auto un = std::get<csg::Union3D>(st.top.objs.at(1));

  CHECK(un.objs.size() == 0);
}

TEST_CASE("extrude and empty group", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(height = 10, center = false, scale = [5, 5]) {
   circle(r = 5);
}
group();
)");
  auto un = std::get<csg::Union3D>(st.top.objs.at(1));

  CHECK(un.objs.size() == 0);
}

TEST_CASE("empty boolean", "[csg]") {
  auto st = *csg::parse_csg(R"(
union();
sphere(r=0.1);
)");
  auto un = std::get<csg::Union3D>(st.top.objs.at(0));
  CHECK(un.objs.size() == 0);

  auto sph = std::get<csg::Sphere>(st.top.objs.at(1));
  CHECK(sph.radius == 0.1);
}

TEST_CASE("empty intersection", "[csg]") {
  auto st = *csg::parse_csg(R"(
intersection();
)");
  auto in = std::get<csg::Intersection3D>(st.top.objs.at(0));
  CHECK(in.objs.size() == 0);
}

TEST_CASE("empty difference", "[csg]") {
  auto st = *csg::parse_csg(R"(
cylinder(h = 200, r1 = 1.7, r2 = 2.7, center = false);
difference();
)");
  auto cone = std::get<csg::Cone>(st.top.objs.at(0));
  CHECK(cone.radius1 == 1.7);
  CHECK(cone.radius2 == 2.7);
  CHECK(cone.height == 200);

  auto diff = std::get<csg::Difference3D>(st.top.objs.at(1));
  CHECK(diff.first_obj == nullptr);
  CHECK(diff.next_objs.objs.size() == 0);
}

TEST_CASE("empty matmul", "[csg]") {
  auto st = *csg::parse_csg(R"(
multmatrix(
[
[1, 0, 0, 0.0020],
[0, 1, 0, 0.0005],
[0, 0, 1, 0.0005],
[0, 0, 0, 1]
]
);
)");
  auto mat = std::get<csg::Mulmatrix3D>(st.top.objs.back());
  CHECK(mat.group.objs.size() == 0);
  CHECK(mat.translation == std::array<double, 3>({0.0020, 0.0005, 0.0005}));
  CHECK(mat.rotation()[0] == std::array<double, 3>({1, 0, 0}));
  CHECK(mat.rotation()[1] == std::array<double, 3>({0, 1, 0}));
  CHECK(mat.rotation()[2] == std::array<double, 3>({0, 0, 1}));
}

TEST_CASE("empty intersection inside linear extrude", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(
height = 10,
center = true,
scale = [10,1]
) {
    circle(r=0.2);
    intersection();
}
)");
  auto lin_ext = std::get<csg::LinearExtrude>(st.top.objs.at(0));
  CHECK(lin_ext.group.objs.size() == 2);

  auto cir = std::get<csg::Circle>(lin_ext.group.objs.at(0));
  CHECK(cir.radius == 0.2);

  auto in = std::get<csg::Intersection2D>(lin_ext.group.objs.at(1));
  CHECK(in.objs.size() == 0);
}

TEST_CASE("empty linear_extrude", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(height = 10, center = true, scale = [10,1]);
)");
  auto lin_ext = std::get<csg::LinearExtrude>(st.top.objs.at(0));
  CHECK(lin_ext.height == 10);
  CHECK(lin_ext.center == true);

  auto [sx, sy] = lin_ext.scale;
  CHECK(sx == 10);
  CHECK(sy == 1);
}

TEST_CASE("group and empty multmatrix", "[csg]") {
  auto st = *csg::parse_csg(R"(
group() {
   multmatrix([[1, 0, 0, 0.2], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) {
      cylinder($fn = 50, $fa = 12, $fs = 2, h = 0.636, r1 = 0.05, r2 = 0.05, center = false);
   }
}
multmatrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0.636], [0, 0, 0, 1]]);
)");
  auto gp = std::get<csg::Union3D>(st.top.objs.at(0));
  auto mm1 = std::get<csg::Mulmatrix3D>(gp.objs.at(0));
  auto cone = std::get<csg::Cone>(mm1.group.objs.at(0));
  CHECK(cone.height == 0.636);

  auto mm2 = std::get<csg::Mulmatrix3D>(st.top.objs.at(1));
  CHECK(mm2.group.objs.size() == 0);
}

} // namespace
