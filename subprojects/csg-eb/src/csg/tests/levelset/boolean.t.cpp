#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

#include <memory>

namespace {

TEST_CASE("Union3D", "[Levelset Boolean]") {
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
  csg::CsgIF my_levelset(my_tree);

  CHECK(0 < my_levelset(0, 7, 0));
  CHECK(Approx(8 * 8 - 7 * 7) == my_levelset(0, 7, 0));

  CHECK(0 < my_levelset(0, 0, -7.5));
  CHECK(Approx(8 * 8 - 7.5 * 7.5) == my_levelset(0, 0, -7.5));

  CHECK_FALSE(0 < my_levelset(0, -9, 0));
  CHECK(Approx(-3) == my_levelset(0, -9, 0));

  CHECK_FALSE(0 < my_levelset(9, 0, 0));
  CHECK(Approx(-3) == my_levelset(9, 0, 0));
}

TEST_CASE("Intersection3D", "[Levelset Boolean]") {
  // Create an of a cube of size 12 and a sphere of radius 8
  csg::Cube my_cub;
  my_cub.size = {12, 12, 12};
  my_cub.center = true;

  csg::Sphere my_sph;
  my_sph.radius = 8;

  auto my_in = csg::Intersection3D();
  my_in.objs.push_back(my_cub);
  my_in.objs.push_back(my_sph);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_in);
  csg::CsgIF my_levelset(my_tree);

  CHECK_FALSE(0 < my_levelset(0, 7, 0));
  CHECK(Approx(-1) == my_levelset(0, 7, 0));

  CHECK_FALSE(0 < my_levelset(0, 0, -7.5));
  CHECK(Approx(-1.5) == my_levelset(0, 0, -7.5));

  CHECK(0 < my_levelset(0, 5, 0));
  CHECK(Approx(1) == my_levelset(0, 5, 0));

  CHECK(0 < my_levelset(-5, 0, 0));
  CHECK(Approx(1) == my_levelset(-5, 0, 0));
}

TEST_CASE("Difference3D", "[Levelset Boolean]") {
  // Create an of a cube of size 12 and remove a sphere of radius 8
  csg::Cube my_cub;
  my_cub.size = {12, 12, 12};
  my_cub.center = true;

  csg::Sphere my_sph;
  my_sph.radius = 8;

  csg::Union3D my_union;
  my_union.objs.push_back(my_sph);

  auto my_diff = csg::Difference3D({
      std::make_shared<csg::Type3D>(my_cub),
      csg::Union3D(my_union),
  });

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_diff);
  csg::CsgIF my_levelset(my_tree);

  CHECK_FALSE(0 < my_levelset(0, 0, 0));
  CHECK(Approx(-(8 * 8)) == my_levelset(0, 0, 0));

  CHECK_FALSE(0 < my_levelset(-4, 0, 0));
  CHECK(Approx(-(8 * 8 - 4 * 4)) == my_levelset(-4, 0, 0));

  CHECK(0 < my_levelset(5.99, 5.99, 0));
  CHECK(Approx(0.01) == my_levelset(5.99, 5.99, 0));

  CHECK(0 < my_levelset(-5.99, 0, 5.99));
  CHECK(Approx(0.01) == my_levelset(-5.99, 0, 5.99));
}

} // namespace
