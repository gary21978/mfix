#include "catch2/catch_all.hpp"
using namespace Catch;

#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Union.H>

#include "main.hpp"

using namespace amrex;

TEST_CASE("cyclone", "[benchmarks][.]") {

  // EB2::PlaneIF hbody(point, normal);
  // auto hopper0 = EB2::lathe(EB2::makeUnion(hbody, funnel));
  // auto hopper1 = EB2::rotate(hopper0, 2. * std::atan(1.0), 1);
  // auto eb_hopper = EB2::translate(hopper1, center);

  // test_for_equivalence(my_hopper, my_csg, -1, 1, -1, 1, -1, 1);

  auto maybe_st = csg::parse_csg(" \
group() { \
  union() { \
    multmatrix( \
        [ [ 1, 0, 0, 17 ], [ 0, 1, 0, 20 ], [ 0, 0, 1, 5 ], [ 0, 0, 0, 1 ] ]) { \
      cube(size = [ 18, 40, 7 ], center = false); \
    } \
    multmatrix( \
        [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 20 ], [ 0, 0, 1, 20 ], [ 0, 0, 0, 1 ] ]) { \
      union() { \
        difference() { \
          multmatrix([ \
            [ 1, 0, 0, 15 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] \
          ]) { \
            multmatrix([ \
              [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ -1, 0, 0, 0 ], [ 0, 0, 0, 1 ] \
            ]) { \
              cylinder($fn = 0, $fa = 12, $fs = 2, h = 40, r1 = 15, r2 = 15, \
                       center = true); \
            } \
          } \
          multmatrix([ \
            [ 1, 0, 0, 35 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] \
          ]) { \
            multmatrix([ \
              [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ -1, 0, 0, 0 ], [ 0, 0, 0, 1 ] \
            ]) { \
              cylinder($fn = 0, $fa = 12, $fs = 2, h = 40, r1 = 8, r2 = 8, \
                       center = true); \
            } \
          } \
        } \
        multmatrix([ \
          [ 1, 0, 0, 35 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] \
        ]) { \
          multmatrix([ \
            [ 0, 0, 1, 0 ], [ 0, 1, 0, 0 ], [ -1, 0, 0, 0 ], [ 0, 0, 0, 1 ] \
          ]) { \
            cylinder($fn = 0, $fa = 12, $fs = 2, h = 40, r1 = 5.5, r2 = 5.5, \
                     center = true); \
          } \
        } \
      } \
    } \
  } \
} \
");

  CHECK(maybe_st != nullptr);

  csg::CsgIF my_cyclone(maybe_st);
}
