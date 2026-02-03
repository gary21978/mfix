#include <memory>

#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_matrix_functions.hpp>
#include <csg_types.hpp>

namespace {

TEST_CASE("Translation", "[Levelset Transform]") {

  csg::Cylinder my_cyl;
  my_cyl.radius = 2;
  my_cyl.height = 20;
  my_cyl.center = true;

  auto my_mat = csg::Mulmatrix3D({{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}});
  my_mat.translation = {100, 200, 300};
  my_mat.group.objs.push_back(my_cyl);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_mat);
  csg::CsgIF my_levelset(my_tree);

  SECTION("Outside") {
    CHECK_FALSE(0 < my_levelset(100, 200, 280));
    CHECK_FALSE(0 < my_levelset(100, 200, 320));
    CHECK_FALSE(0 < my_levelset(103, 200, 300));
    CHECK_FALSE(0 < my_levelset(97, 200, 300));
    CHECK_FALSE(0 < my_levelset(100, 203, 300));
    CHECK_FALSE(0 < my_levelset(100, 197, 300));
  }

  SECTION("Inside") {
    CHECK(0 < my_levelset(100, 200, 295));
    CHECK(0 < my_levelset(100, 200, 305));
    CHECK(0 < my_levelset(101, 200, 300));
    CHECK(0 < my_levelset(99, 200, 300));
    CHECK(0 < my_levelset(100, 201, 300));
    CHECK(0 < my_levelset(100, 199, 300));
  }
}

TEST_CASE("90° Rotation around y-axis, rotating z-axis into x-axis",
          "[Levelset Transform]") {

  csg::Cylinder my_cyl;
  my_cyl.radius = 2;
  my_cyl.height = 20;
  my_cyl.center = false;

  auto my_mat = csg::Mulmatrix3D({{0, 0, 1}}, {{0, 1, 0}}, {{-1, 0, 0}});
  my_mat.translation = {0, 0, 0};
  my_mat.group.objs.push_back(my_cyl);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_mat);
  csg::CsgIF my_levelset(my_tree);

  SECTION("near origin end") {
    CHECK_FALSE(0 < my_levelset(-1, 0, 0));
    CHECK(0 < my_levelset(1, 0, 0));

    CHECK_FALSE(0 < my_levelset(3, 0, 2.1));
    CHECK(0 < my_levelset(3, 0, 1.9));

    CHECK_FALSE(0 < my_levelset(3, 0, -2.1));
    CHECK(0 < my_levelset(3, 0, -1.9));
  }
  SECTION("near middle") {
    CHECK_FALSE(0 < my_levelset(10, 2.1, 0));
    CHECK(0 < my_levelset(10, 1.9, 0));

    CHECK_FALSE(0 < my_levelset(10, -2.1, 0));
    CHECK(0 < my_levelset(10, -1.9, 0));
  }
  SECTION("near other end") {
    CHECK_FALSE(0 < my_levelset(21, 0, 0));
    CHECK(0 < my_levelset(19, 0, 0));

    CHECK_FALSE(0 < my_levelset(18, 0, 2.1));
    CHECK(0 < my_levelset(18, 0, 1.9));

    CHECK_FALSE(0 < my_levelset(18, 0, -2.1));
    CHECK(0 < my_levelset(18, 0, -1.9));
  }
}

TEST_CASE("90° rotation + translation of cylinder", "[Levelset Transform]") {

  double radius = 4.5, height = 10;
  double Cx = 10, Cy = 5, Cz = 5;

  csg::Cylinder my_cyl;
  my_cyl.radius = radius;
  my_cyl.height = height;
  my_cyl.center = true;

  auto my_mat = csg::Mulmatrix3D({{0, 0, 1}}, {{0, 1, 0}}, {{-1, 0, 0}});
  my_mat.translation = {Cx, Cy, Cz};
  my_mat.group.objs.push_back(my_cyl);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_mat);
  csg::CsgIF my_levelset(my_tree);

  SECTION("near origin end") {
    CHECK_FALSE(0 < my_levelset(Cx - (1.1 * height / 2), Cy, Cz));
    CHECK(0 < my_levelset(Cx - (0.9 * height / 2), Cy, Cz));

    CHECK_FALSE(0 <
                my_levelset(Cx - (0.9 * height / 2), Cy, Cz + 1.1 * radius));
    CHECK(0 < my_levelset(Cx - (0.9 * height / 2), Cy, 0.9 * radius));

    CHECK_FALSE(0 <
                my_levelset(Cx - (0.9 * height / 2), Cy, Cz - 1.1 * radius));
    CHECK(0 < my_levelset(Cx - (0.9 * height / 2), Cy, Cz - 0.9 * radius));
  }
  SECTION("near middle") {
    CHECK_FALSE(0 < my_levelset(Cx, Cy + 1.1 * radius, Cz));
    CHECK(0 < my_levelset(Cx, Cy + 0.9 * radius, Cz));

    CHECK_FALSE(0 < my_levelset(Cx, Cy - 1.1 * radius, Cz));
    CHECK(0 < my_levelset(Cx, Cy - 0.9 * radius, Cz));
  }
  SECTION("near other end") {
    CHECK_FALSE(0 < my_levelset(Cx + (1.1 * height / 2), Cy, Cz));
    CHECK(0 < my_levelset(Cx + (0.9 * height / 2), Cy, Cz));

    CHECK_FALSE(0 <
                my_levelset(Cx + (0.9 * height / 2), Cy, Cz + 1.1 * radius));
    CHECK(0 < my_levelset(Cx + (0.9 * height / 2), Cy, Cz + 0.9 * radius));

    CHECK_FALSE(0 <
                my_levelset(Cx + (0.9 * height / 2), Cy, Cz - 1.1 * radius));
    CHECK(0 < my_levelset(Cx + (0.9 * height / 2), Cy, Cz - 0.9 * radius));
  }
}

TEST_CASE("45° Rotation around y-axis", "[Levelset Transform]") {

  double side = 20;
  csg::Cube my_cub;
  my_cub.size = {side, side, side};
  my_cub.center = true;

  auto my_mat = csg::Mulmatrix3D({{0.707107, 0, 0.707107}}, {{0, 1, 0}},
                                 {{-0.707107, 0, 0.707107}});
  my_mat.translation = {0, 0, 0};
  my_mat.group.objs.push_back(my_cub);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_mat);
  csg::CsgIF my_levelset(my_tree);

  SECTION("along x axis") {
    CHECK_FALSE(0 < my_levelset(-14.15, 0, 0));
    CHECK(0 < my_levelset(-14.13, 0, 0));

    CHECK_FALSE(0 < my_levelset(14.15, 0, 0));
    CHECK(0 < my_levelset(14.13, 0, 0));

    CHECK_FALSE(0 < my_levelset(14.13, 0, 1));
    CHECK(0 < my_levelset(12, 0, 1));

    CHECK_FALSE(0 < my_levelset(-14.13, 0, 1));
    CHECK(0 < my_levelset(-12, 0, 1));

    CHECK_FALSE(0 < my_levelset(-14.13, 0, -1));
    CHECK(0 < my_levelset(-12, 0, -1));

    CHECK_FALSE(0 < my_levelset(14.15, 1, 0));
    CHECK(0 < my_levelset(14.13, 1, 0));
  }
  SECTION("along y axis") {
    CHECK_FALSE(0 < my_levelset(0, -10.01, 0));
    CHECK(0 < my_levelset(0, -9.99, 0));

    CHECK_FALSE(0 < my_levelset(0, 10.01, 0));
    CHECK(0 < my_levelset(0, 9.99, 0));

    CHECK_FALSE(0 < my_levelset(0, 10.01, 1));
    CHECK(0 < my_levelset(0, 9.99, 1));

    CHECK_FALSE(0 < my_levelset(0, -10.01, -1));
    CHECK(0 < my_levelset(0, -9.99, -1));

    CHECK_FALSE(0 < my_levelset(1, 10.01, 0));
    CHECK(0 < my_levelset(1, 9.99, 0));
  }
  SECTION("along z axis") {
    CHECK_FALSE(0 < my_levelset(0, 0, -14.15));
    CHECK(0 < my_levelset(0, 0, -14.13));

    CHECK_FALSE(0 < my_levelset(0, 0, 14.15));
    CHECK(0 < my_levelset(0, 0, 14.13));

    CHECK_FALSE(0 < my_levelset(0, 1, 14.15));
    CHECK(0 < my_levelset(0, 1, 14.13));

    CHECK_FALSE(0 < my_levelset(0, -1, -14.15));
    CHECK(0 < my_levelset(0, -1, -14.13));

    CHECK_FALSE(0 < my_levelset(1, 0, 14.13));
    CHECK(0 < my_levelset(1, 0, 12));
  }
}

TEST_CASE("Skewed with shear y along z", "[Levelset Transform]") {
  double radius = 10, height = 10;
  double skew_yz = 0.7;

  csg::Cylinder my_cyl;
  my_cyl.radius = radius;
  my_cyl.height = height,
  my_cyl.center = false;

  auto my_mat = csg::Mulmatrix3D({{1, 0, 0}}, {{0, 1, skew_yz}}, {{0, 0, 1}});
  my_mat.translation = {0, 0, 0};
  my_mat.group.objs.push_back(my_cyl);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_mat);
  csg::CsgIF my_levelset(my_tree);

  auto check = [radius, &my_levelset](double Cx, double Cy, double Cz) {
    CHECK_FALSE(0 < my_levelset(Cx + radius * 1.01, Cy, Cz));
    CHECK(0 < my_levelset(Cx + radius * 0.99, Cy, Cz));

    CHECK_FALSE(0 < my_levelset(Cx, Cy + 1.01 * radius, Cz));
    CHECK(0 < my_levelset(Cx, Cy + radius * 0.99, Cz));

    CHECK_FALSE(0 < my_levelset(Cx, Cy + 1.01 * radius, Cz));
    CHECK(0 < my_levelset(Cx, Cy + radius * 0.99, Cz));

    CHECK_FALSE(0 < my_levelset(Cx, Cy - 1.01 * radius, Cz));
    CHECK(0 < my_levelset(Cx, Cy - radius * 0.99, Cz));
  };

  SECTION("near lower end") { check(0, 0, 0.01); }
  SECTION("near middle") { check(0, skew_yz * 4.99, 4.99); }
  SECTION("near upper end") { check(0, skew_yz * 9.99, 9.99); }
}

} // namespace
