#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

#include <memory>

namespace {

TEST_CASE("Empty Difference3D", "[Levelset Boolean]") {
  auto my_diff = csg::Difference3D();

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_diff);
  csg::CsgIF my_levelset(my_tree);

  // Any point should lie outside!
  CHECK_FALSE(0 < my_levelset(0, 0, 0));
  CHECK_FALSE(0 < my_levelset(0, 7, 0));
  CHECK_FALSE(0 < my_levelset(0, 0, -7.5));
  CHECK_FALSE(0 < my_levelset(0, -9, 0));
  CHECK_FALSE(0 < my_levelset(9, 0, 0));
}

TEST_CASE("Empty Difference2D", "[Levelset Boolean]") {
  double height = 100, radius = 10;
  csg::LinearExtrude my_lin_ext;
  my_lin_ext.height = 100;
  my_lin_ext.center = false;
  my_lin_ext.scale = {1, 1};
  my_lin_ext.group = csg::Union2D();

  auto my_diff = csg::Difference2D();

  my_lin_ext.center = false;
  my_lin_ext.group.objs.push_back(my_diff);
  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_lin_ext);
  csg::CsgIF my_levelset(my_tree);

  // Any point should lie outside!
  CHECK_FALSE(0 < my_levelset(0, 0, 0));
  CHECK_FALSE(0 < my_levelset(0, 7, 0));
  CHECK_FALSE(0 < my_levelset(0, 0, -7.5));
  CHECK_FALSE(0 < my_levelset(0, -9, 0));
  CHECK_FALSE(0 < my_levelset(9, 0, 0));
}

} // namespace
