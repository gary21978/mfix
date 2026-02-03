#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

#include <memory>

namespace {

TEST_CASE("union internal flow", "[Levelset Internal Flow]") {
  // Create a union of a cube of size 12 and a sphere of radius 8
  csg::Cube my_cub;
  my_cub.size = {12, 12, 12};
  my_cub.center = true;

  csg::Sphere my_sph;
  my_sph.radius = 8;

  auto my_union = csg::Union3D();
  my_union.objs.push_back(my_cub);
  my_union.objs.push_back(my_sph);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_union);
  csg::CsgIF my_ls(my_tree);
  csg::CsgIF my_ls_comp(my_tree, true);

  auto check = [&my_ls, &my_ls_comp](double xx, double yy, double zz) {
    CHECK(my_ls(xx, yy, zz) == -my_ls_comp(xx, yy, zz));
  };

  SECTION("Inside") {
    check(0, 7, 0);
    check(0, 0, -7.5);
  }
  SECTION("Outside") {
    check(0, -9, 0);
    check(9, 0, 0);
  }
}

TEST_CASE("CsgIF copy constructor", "[Levelset Internal Flow]") {

  csg::Sphere my_sph;
  my_sph.radius = 10;

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_sph);
  csg::CsgIF my_ls1(my_tree);
  csg::CsgIF my_ls2(my_ls1);
  csg::CsgIF my_ls_comp1(my_tree, true);
  csg::CsgIF my_ls_comp2(my_ls_comp1);

  auto check = [&my_ls1, &my_ls2, &my_ls_comp1,
                &my_ls_comp2](double xx, double yy, double zz) {
    CHECK(my_ls1(xx, yy, zz) == my_ls2(xx, yy, zz));
    CHECK(my_ls_comp1(xx, yy, zz) == my_ls_comp2(xx, yy, zz));
    CHECK(my_ls1(xx, yy, zz) == -my_ls_comp1(xx, yy, zz));
    CHECK(my_ls2(xx, yy, zz) == -my_ls_comp2(xx, yy, zz));
  };

  SECTION("Outside") {
    check(0, 8, 8);
    check(-8, 0, 8);
    check(0, 8, -8);
  }
  SECTION("Inside") {
    check(0, 0, 0);
    check(0, 7, 7);
    check(-7, 0, 7);
    check(-7, 0, 7);
    check(0, 7, -7);
  }
}

} // namespace
