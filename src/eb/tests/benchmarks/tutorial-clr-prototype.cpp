#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch_all.hpp"
using namespace Catch;

#include <fstream>
#include <sstream>

#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Union.H>

#include "main.hpp"

namespace {

using namespace amrex;

// This is used a padding on for the Y and Z directions.
const Real offset = 0.16;

// Air-reactor parameters
const Real reactor_radius = 0.10;
const Real reactor_height = 0.75;

// Riser parameters
const Real riser_radius = 0.05;
const Real riser_height = 3.75;

// Height of the cone connecting the air-reactor and riser
const Real c2c_height = 0.10;

const Real crossover_radius = 0.05;
const Real crossover_height = 3.60;
const Real crossover_length = 0.636;

// Cyclone parameters
const Real cyclone_top = 3.80;
const Real cyclone_radius = 0.10;
const Real cyclone_height = 0.35;
const Real cyclone_dipleg_radius = 0.04;
const Real cyclone_dipleg_bottom = 0.15;
const Real cyclone_to_dipleg_height = 0.20;

// Loopseal parameters
const Real loopseal_radius = 0.12;
const Real loopseal_bottom = 2.24;
const Real loopseal_height = 0.50;
const Real loopseal_offset = 0.004;
const Real loopseal_top = loopseal_bottom + loopseal_height;

// Fuel reactor parameters
const Real fuelreactor_radius = 0.12;
const Real fuelreactor_bottom = 0.52;
const Real fuelreactor_height = 1.54;
const Real fuelreactor_top = fuelreactor_bottom + fuelreactor_height;

// Lvalve parameters
const Real lvalve_height = 0.2;

template <class H> void check_clr(const H &clr) {

  {
    auto Cx = 0.0, Cy = offset, Cz = offset;
    INFO("Air reactor");
    {
      INFO("Along the air reactor body (cylinder)");

      CHECK_FALSE(
          0 < clr(Cx + 0.99 * reactor_height, Cy + 0.99 * reactor_radius, Cz));
      CHECK(0 <
            clr(Cx + 0.99 * reactor_height, Cy + 1.01 * reactor_radius, Cz));
      CHECK_FALSE(
          0 < clr(Cx + 0.99 * reactor_height, Cy, Cz + 0.99 * reactor_radius));
      CHECK(0 <
            clr(Cx + 0.99 * reactor_height, Cy, Cz + 1.01 * reactor_radius));
      CHECK_FALSE(
          0 < clr(Cx + 0.01 * reactor_height, Cy - 0.99 * reactor_radius, Cz));
      CHECK(0 <
            clr(Cx + 0.01 * reactor_height, Cy - 1.01 * reactor_radius, Cz));
      CHECK_FALSE(
          0 < clr(Cx + 0.01 * reactor_height, Cy, Cz - 0.99 * reactor_radius));
      CHECK(0 <
            clr(Cx + 0.01 * reactor_height, Cy, Cz - 1.01 * reactor_radius));
    }

    {
      INFO("Along the air reactor cap (truncated cone)");
      CHECK_FALSE(0 < clr(Cx + reactor_height + 0.5 * c2c_height,
                          Cy + 0.99 * ((reactor_radius + riser_radius) / 2),
                          Cz));
      CHECK(0 < clr(Cx + reactor_height + 0.5 * c2c_height,
                    Cy + 1.01 * ((reactor_radius + riser_radius) / 2), Cz));
      CHECK_FALSE(0 < clr(Cx + reactor_height + 0.5 * c2c_height, Cy,
                          Cz - 0.99 * ((reactor_radius + riser_radius) / 2)));
      CHECK(0 < clr(Cx + reactor_height + 0.5 * c2c_height, Cy,
                    Cz - 1.01 * ((reactor_radius + riser_radius) / 2)));
    }

    {
      INFO("Along the riser above cap (cylinder)");
      auto height_above_cap = riser_height - reactor_height - c2c_height;
      CHECK_FALSE(
          0 < clr(Cx + reactor_height + c2c_height + height_above_cap * 0.01,
                  Cy + 0.99 * riser_radius, Cz));
      CHECK(0 < clr(Cx + reactor_height + c2c_height + height_above_cap * 0.01,
                    Cy + 1.01 * riser_radius, Cz));
      CHECK_FALSE(
          0 < clr(Cx + reactor_height + c2c_height + height_above_cap * 0.99,
                  Cy, Cz - 0.99 * riser_radius));
      CHECK(0 <
            clr(Cx + height_above_cap * 0.99, Cy, Cz - 1.01 * riser_radius));

      // top blocked off
      CHECK(0 < clr(Cx + reactor_height + c2c_height + height_above_cap * 1.01,
                    Cy, Cz));
    }
  }

  {
    INFO("Crossover");
    auto Cx = crossover_height, Cy = offset, Cz = offset;
    CHECK_FALSE(0 < clr(Cx, Cy, Cz + riser_radius + 0.01 * crossover_length));
    CHECK(0 < clr(Cx + 1.01 * crossover_radius, Cy,
                  Cz + riser_radius + 0.01 * crossover_length));

    CHECK_FALSE(0 < clr(Cx, Cy - 0.99 * crossover_radius,
                        Cz + riser_radius + 0.01 * crossover_length));
    CHECK(0 < clr(Cx, Cy - 1.01 * crossover_radius,
                  Cz + riser_radius + 0.01 * crossover_length));

    CHECK_FALSE(0 < clr(Cx + 0.99 * crossover_radius, Cy,
                        Cz + crossover_length - cyclone_radius -
                            0.01 * crossover_length));
    CHECK(0 < clr(Cx + 1.01 * crossover_radius, Cy,
                  Cz + crossover_length - cyclone_radius -
                      0.01 * crossover_length));
  }

  {
    INFO("Cyclone with dipleg");
    auto Cx = cyclone_top, Cy = offset, Cz = offset + crossover_length;
    {
      INFO("Along the cyclone body");
      CHECK_FALSE(
          0 < clr(Cx - 0.2 * cyclone_height, Cy + 0.99 * cyclone_radius, Cz));
      CHECK(0 < clr(Cx - 0.2 * cyclone_height, Cy + 1.01 * cyclone_radius, Cz));

      CHECK_FALSE(
          0 < clr(Cx - 0.2 * cyclone_height, Cy, Cz - 0.99 * cyclone_radius));
      CHECK(0 < clr(Cx - 0.2 * cyclone_height, Cy, Cz - 1.01 * cyclone_radius));
    }
    {
      INFO("Along the cyclone cap");
      CHECK_FALSE(
          0 < clr(Cx - cyclone_height - cyclone_to_dipleg_height / 2,
                  Cy + 0.99 * ((cyclone_radius + cyclone_dipleg_radius) / 2),
                  Cz));
      CHECK(0 < clr(Cx - cyclone_height - cyclone_to_dipleg_height / 2,
                    Cy + 1.01 * ((cyclone_radius + cyclone_dipleg_radius) / 2),
                    Cz));
    }
    {
      INFO("Along the dipleg");
      CHECK_FALSE(0 < clr(Cx - cyclone_height - 2 * cyclone_to_dipleg_height,
                          Cy + 0.99 * cyclone_dipleg_radius, Cz));
      CHECK(0 < clr(Cx - cyclone_height - 2 * cyclone_to_dipleg_height,
                    Cy + 1.01 * cyclone_dipleg_radius, Cz));

      CHECK_FALSE(0 < clr(Cx - cyclone_height - 2 * cyclone_to_dipleg_height,
                          Cy, Cz - 0.99 * cyclone_dipleg_radius));
      CHECK(0 < clr(Cx - cyclone_height - 2 * cyclone_to_dipleg_height, Cy,
                    Cz - 1.01 * cyclone_dipleg_radius));

      CHECK_FALSE(0 < clr(1.01 * cyclone_dipleg_bottom, Cy, Cz));
      CHECK(0 < clr(0.99 * cyclone_dipleg_bottom, Cy, Cz));

      CHECK_FALSE(0 < clr(Cx + 0.01 * cyclone_height,
                          Cy + 0.99 * cyclone_dipleg_radius, Cz));
      CHECK(0 < clr(Cx + 0.01 * cyclone_height,
                    Cy + 1.01 * cyclone_dipleg_radius, Cz));
    }
  }

  {
    INFO("Loopseal");
    auto Cx = loopseal_bottom, Cy = offset, Cz = offset + crossover_length + loopseal_offset;
    CHECK_FALSE(0 <
                clr(Cx + loopseal_height / 2, Cy + 0.99 * loopseal_radius, Cz));
    CHECK(0 < clr(Cx + loopseal_height / 2, Cy + 1.01 * loopseal_radius, Cz));

    CHECK_FALSE(0 <
                clr(Cx + loopseal_height / 2, Cy, Cz - 0.99 * loopseal_radius));
    CHECK(0 < clr(Cx + loopseal_height / 2, Cy, Cz - 1.01 * loopseal_radius));
  }

  {
    INFO("Fuel reactor");
    auto Cx = fuelreactor_bottom, Cy = offset, Cz = offset + crossover_length;
    CHECK_FALSE(
        0 < clr(Cx + fuelreactor_height / 2, Cy + 0.99 * loopseal_radius, Cz));
    CHECK(0 <
          clr(Cx + fuelreactor_height / 2, Cy + 1.01 * loopseal_radius, Cz));

    CHECK_FALSE(
        0 < clr(Cx + fuelreactor_height / 2, Cy, Cz - 0.99 * loopseal_radius));
    CHECK(0 <
          clr(Cx + fuelreactor_height / 2, Cy, Cz - 1.01 * loopseal_radius));
  }

  {
    INFO("Lvalve");
    auto Cx = lvalve_height, Cy = offset, Cz = offset;
    CHECK_FALSE(0 < clr(Cx + 0.99 * crossover_radius, Cy,
                        Cz + reactor_radius + 0.01 * crossover_length));
    CHECK(0 < clr(Cx + 1.01 * crossover_radius, Cy,
                  Cz + reactor_radius + 0.01 * crossover_length));

    CHECK_FALSE(0 < clr(Cx, Cy - 0.99 * crossover_radius,
                        Cz + reactor_radius + 0.01 * crossover_length));
    CHECK(0 < clr(Cx, Cy - 1.01 * crossover_radius,
                  Cz + reactor_radius + 0.01 * crossover_length));

    CHECK_FALSE(0 < clr(Cx + 0.99 * crossover_radius, Cy,
                        Cz + crossover_length - cyclone_dipleg_radius -
                            0.01 * crossover_length));
    CHECK(0 < clr(Cx + 1.01 * crossover_radius, Cy,
                  Cz + crossover_length - cyclone_dipleg_radius -
                      0.01 * crossover_length));
  }
}

template <class H> void time_clr(const H &clr, int N) {
  std::vector<double> XX(N), YY(N), ZZ(N);
  std::generate(XX.begin(), XX.end(),
                [i = 0, N]() mutable { return (++i) * (riser_height / N); });
  std::generate(YY.begin(), YY.end(), [i = 0, N]() mutable {
    return (++i) * (fuelreactor_radius / N);
  });
  std::generate(ZZ.begin(), ZZ.end(), [i = 0, N]() mutable {
    return (++i) * (crossover_length / N);
  });

  for (auto xx : XX) {
    for (auto yy : YY) {
      for (auto zz : ZZ) {
        (void)clr(xx, yy, zz);
      }
    }
  }
}

TEST_CASE("clr", "[benchmarks]") {
  SECTION("EB CLR") {
    // set up ebfactory

    Array<Real, 3> point, normal, center;

    /****************************************************************************
     * Build the air reactor's reactor section *
     ***************************************************************************/

    // Define point+normal to define a cylinder
    point = {reactor_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF reactor_body(point, normal);

    // Define point+normal to define the cone cap
    point = {reactor_radius, reactor_height, 0.0};
    normal = {c2c_height, reactor_radius - riser_radius, 0.0};

    EB2::PlaneIF reactor_cap(point, normal);

    // Create reactor with cap
    auto reactor_wcap = EB2::makeUnion(reactor_body, reactor_cap);

    /****************************************************************************
     * Build the air reactor's riser section *
     ***************************************************************************/

    // Define point+normal to define a cylinder
    point = {riser_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF riser_body(point, normal);

    // Define point+normal to define a reactor top
    point = {0.0, riser_height, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF riser_top(point, normal);

    // Create reactor by clipping off top
    auto riser = EB2::makeUnion(riser_body, riser_top);

    /****************************************************************************
     * Combine the reactor and riser sections and make 3D *
     ***************************************************************************/

    auto airreactor0 = EB2::lathe(EB2::makeIntersection(reactor_wcap, riser));

    // Align to correct axis
    auto airreactor1 = EB2::rotate(airreactor0, 2.0 * std::atan(1.0), 1);

    // Translate to correct location
    center = {0.0, offset, offset};
    auto airreactorIF = EB2::translate(airreactor1, center);

    /****************************************************************************
     * Build a horizontal pipe. This is used to connect the riser to the *
     * cyclone as well as a portion of the l-valve. *
     ***************************************************************************/

    // Riser parameters
    Real pipe_radius = crossover_radius;

    // Define point+normal to define a cylinder
    point = {pipe_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF pipeIF_body(point, normal);

    // Define point+normal to define a reactor top
    point = {0.0, offset, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF pipeIF_lo(point, normal);

    // Define point+normal to define a reactor top
    point = {0.0, offset + crossover_length, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF pipeIF_hi(point, normal);

    auto pipeIF = EB2::lathe(EB2::makeUnion(pipeIF_body, pipeIF_lo, pipeIF_hi));

    /****************************************************************************
     * Build crossover to connect air-reactor riser to the cyclone. *
     ***************************************************************************/

    // Translate to correct location
    center = {crossover_height, offset, 0.0};

    auto crossoverIF = EB2::translate(pipeIF, center);

    /****************************************************************************
     * Build the cyclone *
     ***************************************************************************/

    Real cyclone_bottom = cyclone_top - cyclone_height;

    // Define point+normal to define the top of the cyclone
    point = {0.0, cyclone_top, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF cyclone_cap(point, normal);

    // Define point+normal to define a cylinder
    point = {cyclone_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF cyclone_body(point, normal);

    // Define point+normal to define the cone cap
    point = {cyclone_radius, cyclone_bottom, 0.0};
    normal = {cyclone_to_dipleg_height, cyclone_dipleg_radius - cyclone_radius,
              0.0};

    EB2::PlaneIF cyclone_c2c(point, normal);

    // Create reactor with cap
    auto cyclone_wc2c = EB2::makeUnion(cyclone_body, cyclone_cap, cyclone_c2c);

    // Define point+normal to define the top of the cyclone dipleg
    point = {0.0, cyclone_bottom, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF cyclone_dipleg_top(point, normal);

    // Define point+normal to define the bottom of the cyclone dipleg
    point = {0.0, cyclone_dipleg_bottom, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF cyclone_dipleg_bottomIF(point, normal);

    // Define point+normal to define the radius of the cyclone dipleg
    point = {cyclone_dipleg_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF cyclone_dipleg_body(point, normal);

    // Create reactor with cap
    auto cyclone_dipleg =
        EB2::makeUnion(cyclone_dipleg_bottomIF, cyclone_dipleg_body);

    auto cyclone_wdipleg = EB2::makeIntersection(cyclone_wc2c, cyclone_dipleg);

    /****************************************************************************
     * Combine the reactor and riser sections and make 3D *
     ***************************************************************************/

    auto cyclone0 = EB2::lathe(cyclone_wdipleg);

    // Align to correct axis
    auto cyclone1 = EB2::rotate(cyclone0, 2.0 * std::atan(1.0), 1);

    // Translate to correct location
    center = {0.0, offset, offset + crossover_length};
    auto cyclone2 = EB2::translate(cyclone1, center);

    /****************************************************************************
     * Create the loop-seal *
     ***************************************************************************/

    // Define point+normal to define the top of the cyclone dipleg
    point = {0.0, loopseal_top, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF loopsealIF_top(point, normal);

    // Define point+normal to define the bottom of the cyclone dipleg
    point = {0.0, loopseal_bottom, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF loopsealIF_bottom(point, normal);

    // Define point+normal to define the radius of the cyclone dipleg
    point = {loopseal_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF loopsealIF_body(point, normal);

    // Create reactor with cap
    auto loopsealIF =
        EB2::makeUnion(loopsealIF_bottom, loopsealIF_top, loopsealIF_body);

    auto loopseal0 = EB2::lathe(loopsealIF);

    // Align to correct axis
    auto loopseal1 = EB2::rotate(loopseal0, 2.0 * std::atan(1.0), 1);

    // Translate to correct location
    center = {0.0, offset, offset + crossover_length + loopseal_offset};
    auto loopseal2 = EB2::translate(loopseal1, center);

    /****************************************************************************
     * Create the loop-seal *
     ***************************************************************************/

    // Define point+normal to define the top of the cyclone dipleg
    point = {0.0, fuelreactor_top, 0.0};
    normal = {0.0, 1.0, 0.0};

    EB2::PlaneIF fuelreactorIF_top(point, normal);

    // Define point+normal to define the bottom of the cyclone dipleg
    point = {0.0, fuelreactor_bottom, 0.0};
    normal = {0.0, -1.0, 0.0};

    EB2::PlaneIF fuelreactorIF_bottom(point, normal);

    // Define point+normal to define the radius of the cyclone dipleg
    point = {fuelreactor_radius, 0.0, 0.0};
    normal = {1.0, 0.0, 0.0};

    EB2::PlaneIF fuelreactorIF_body(point, normal);

    // Create reactor with cap
    auto fuelreactorIF = EB2::makeUnion(fuelreactorIF_bottom, fuelreactorIF_top,
                                        fuelreactorIF_body);

    auto fuelreactor0 = EB2::lathe(fuelreactorIF);

    // Align to correct axis
    auto fuelreactor1 = EB2::rotate(fuelreactor0, 2.0 * std::atan(1.0), 1);

    // Translate to correct location
    center = {0.0, offset, offset + 0.636};
    auto fuelreactor2 = EB2::translate(fuelreactor1, center);

    /****************************************************************************
     * Create the loop-seal *
     ***************************************************************************/

    center = {lvalve_height, offset, 0.0};
    auto lvalve = EB2::translate(pipeIF, center);

    auto full_clr = EB2::makeIntersection(airreactorIF, crossoverIF, cyclone2,
                                          loopseal2, fuelreactor2, lvalve);

    check_clr(full_clr);

    BENCHMARK("Build EB geometry for 3x3x3 points") { time_clr(full_clr, 3); };
    BENCHMARK("Build EB geometry for 10x10x10 points") {
      time_clr(full_clr, 10);
    };
  }

  SECTION("CSG CLR") {
    std::ifstream csg_file(TUTORIAL_CLR_CSG);
    std::stringstream buffer;
    buffer << csg_file.rdbuf();

    auto maybe_st = csg::parse_csg(buffer.str());
    CHECK(maybe_st != nullptr);

    csg::CsgIF my_clr(maybe_st, true);

    check_clr(my_clr);

    BENCHMARK("Build CSG geometry for 3x3x3 points") { time_clr(my_clr, 3); };
    BENCHMARK("Build CSG geometry for 10x10x10 points") {
      time_clr(my_clr, 10);
    };
  }
}

} // namespace
