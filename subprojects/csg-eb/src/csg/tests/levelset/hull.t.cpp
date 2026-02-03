#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

namespace {

TEST_CASE("hull of a sphere and translated cube", "[Levelset Hull]") {
  csg::Hull my_hull;

  double XX = 1.0, YY = 2.0, ZZ = 3.0, Cx = 5.0, Cy = 5.0, Cz = 6.0;

  SECTION("Unit radius"){
     double R = 1.0;
     my_hull.cube_size = {XX, YY, ZZ};
     my_hull.cube_center = {Cx, Cy, Cz};
     my_hull.sphere_radius = R;
     my_hull.generate_polyhedron();

     auto my_tree = std::make_shared<csg::Tree>();
     my_tree->top.objs.push_back(my_hull);
     csg::CsgIF my_levelset(my_tree);

     SECTION("Outside") {
       CHECK_FALSE(0 < my_levelset(0, 0, Cz));
       CHECK_FALSE(0 < my_levelset(1.3 * R, 0, 0));
       CHECK_FALSE(0 < my_levelset(Cx + 1.1 * XX / 2, Cy, Cz));
       CHECK_FALSE(0 < my_levelset(1.01*(Cx + 0.5*XX), Cy, Cz));
     }
     SECTION("Inside") {
       CHECK(0 < my_levelset(0, 0, 0));
       CHECK(0 < my_levelset(Cx, Cy, Cz));
       CHECK(0 < my_levelset(1.1 * R, 0, 0));
       CHECK(0 < my_levelset(0.99*(Cx + 0.5*XX), Cy, Cz));
     }
  }
  SECTION("Radius 0.1"){
     double R = 0.1;
     my_hull.cube_size = {XX, YY, ZZ};
     my_hull.cube_center = {Cx, Cy, Cz};
     my_hull.sphere_radius = R;
     my_hull.generate_polyhedron();

     auto my_tree = std::make_shared<csg::Tree>();
     my_tree->top.objs.push_back(my_hull);
     csg::CsgIF my_levelset(my_tree);

     SECTION("Outside") {
       CHECK_FALSE(0 < my_levelset(0, 0, Cz));
       CHECK_FALSE(0 < my_levelset(1.4 * R, 0, 0));
       CHECK_FALSE(0 < my_levelset(Cx + 1.1 * XX / 2, Cy, Cz));
       CHECK_FALSE(0 < my_levelset(1.01*(Cx + 0.5*XX), Cy, Cz));
     }
     SECTION("Inside") {
       CHECK(0 < my_levelset(0, 0, 0));
       CHECK(0 < my_levelset(Cx, Cy, Cz));
       CHECK(0 < my_levelset(1.1 * R, 0, 0));
       CHECK(0 < my_levelset(0.99*(Cx + 0.5*XX), Cy, Cz));
     }
  }
}

} // namespace
