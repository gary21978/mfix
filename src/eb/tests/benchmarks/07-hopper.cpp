#include "catch2/catch_all.hpp"
using namespace Catch;

#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Union.H>

#include "main.hpp"

namespace {

using namespace amrex;

const Real orifice_r = 0.00025;
const Real funnel_r = 0.00050;
const Real funnel_h = 0.00105;
const Real Cx = 0.0;
const Real Cy = 0.0;
const Real Cz = 0.0;

template <class H> void check_hopper(const H &hopper) {
  {
    INFO("Distant");
    CHECK(0 < hopper(Cx - 100 * funnel_h, Cy, Cz));
    CHECK(0 < hopper(Cx, Cy + 100 * funnel_h, Cz));
    CHECK(0 < hopper(Cx, Cy - 100 * funnel_h, Cz));
    CHECK(0 < hopper(Cx, Cy, Cz + 100 * funnel_h));
    CHECK(0 < hopper(Cx, Cy, Cz - 100 * funnel_h));
    CHECK_FALSE(0 < hopper(Cx + 100 * funnel_h, Cy, Cz));
  }

  {
    INFO("Length of the Cylinder");
    CHECK(0 < hopper(Cx + funnel_h, Cy + 1.01 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h, Cy, Cz + 1.01 * funnel_r));
    CHECK_FALSE(0 < hopper(Cx + 1.01 * funnel_h, Cy + 0.99 * funnel_r, Cz));
    CHECK_FALSE(0 < hopper(Cx + 1.01 * funnel_h, Cy + 0.99 * funnel_r, Cz));
    // CHECK_FALSE(0 < hopper(Cx + funnel_h, Cy, Cz + 0.99 * funnel_r)); TODO:
    // should this pass? CHECK_FALSE(0 < hopper(Cx + funnel_h, Cy, Cz + 0.99 *
    // funnel_r)); TODO: should this pass?
    CHECK(0 < hopper(Cx + 100 * funnel_h, Cy + 1.01 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + 100 * funnel_h, Cy, Cz + 1.01 * funnel_r));
    CHECK_FALSE(0 < hopper(Cx + 100 * funnel_h, Cy + 0.99 * funnel_r, Cz));
    CHECK_FALSE(0 < hopper(Cx + 100 * funnel_h, Cy, Cz + 0.99 * funnel_r));
  }

  {
    INFO("Along the funnel");
    CHECK_FALSE(0 < hopper(Cx + funnel_h / 2, Cy + 0.5 * funnel_r, Cz));
    CHECK_FALSE(0 < hopper(Cx + funnel_h / 2, Cy + 0.6 * funnel_r, Cz));
    CHECK_FALSE(0 < hopper(Cx + funnel_h / 2, Cy + 0.7 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 0.8 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 0.9 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 1.0 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 1.1 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 1.2 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 1.3 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 1.4 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 1.5 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 1.1 * funnel_r, Cz));
    CHECK(0 < hopper(Cx + funnel_h / 2, Cy + 1.1 * funnel_r, Cz));
  }
}

TEST_CASE("hopper", "[benchmarks]") {
  SECTION("EB Hopper") {
    Array<Real, 3> center = {Cx, Cy, Cz};

    Array<Real, 3> point{0.0, 0.0, 0.0};
    point[0] = orifice_r;

    Array<Real, 3> normal{0.0, 0.0, 0.0};
    normal[0] = funnel_h;
    normal[1] = orifice_r - funnel_r;

    EB2::PlaneIF funnel(point, normal);

    point = {funnel_r, 0., 0.};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF hbody(point, normal);
    auto hopper0 = EB2::lathe(EB2::makeUnion(hbody, funnel));
    auto hopper1 = EB2::rotate(hopper0, 2. * std::atan(1.0), 1);
    auto eb_hopper = EB2::translate(hopper1, center);

    check_hopper(eb_hopper);
  }

  SECTION("CSG Hopper") {
    auto maybe_st = csg::parse_csg(R"(
         group() {
            difference() {
               cube(size = [2002, 2002, 2002], center = true);
               multmatrix([[0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1]]) {
                  union() {
                     cylinder(h = 0.00105, r1 = 0.00025, r2 = 0.0005, center = false);
                     multmatrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0.00105], [0, 0, 0, 1]]) {
                        cylinder(h = 1000, r = 0.0005, center = false);
                     }
                  }
               }
            }
         }
      )");
    CHECK(maybe_st != nullptr);
    csg::CsgIF my_hopper(maybe_st);

    check_hopper(my_hopper);
  }
}

} // namespace
