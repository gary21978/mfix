#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>

namespace {

TEST_CASE("Cylinder", "[Levelset Primitives]") {

  csg::Cylinder my_cyl;
  my_cyl.radius = 1;
  my_cyl.height = 2;
  my_cyl.center = false;

  SECTION("Not Centered") {
    my_cyl.center = false;

    auto my_tree = std::make_shared<csg::Tree>();
    my_tree->top.objs.push_back(my_cyl);
    csg::CsgIF my_levelset(my_tree);

    SECTION("Axis") {
      CHECK_FALSE(0 < my_levelset(0, 0, -1));
      CHECK(Approx(-1) == my_levelset(0, 0, -1));

      CHECK(0 < my_levelset(0, 0, 1));
      CHECK(Approx(1) == my_levelset(0, 0, 1));

      CHECK_FALSE(0 < my_levelset(0, 0, 3));
      CHECK(Approx(-1) == my_levelset(0, 0, 3));
    }

    SECTION("Cross section") {
      CHECK_FALSE(0 < my_levelset(-0.8, 0.8, 1));
      CHECK(Approx(-0.28) == my_levelset(-0.8, 0.8, 1));

      CHECK_FALSE(0 < my_levelset(0.8, -0.8, 1));
      CHECK(Approx(-0.28) == my_levelset(0.8, -0.8, 1));

      CHECK(0 < my_levelset(0, 0.8, 1));
      CHECK(Approx(0.36) == my_levelset(0, 0.8, 1));

      CHECK(0 < my_levelset(0.8, 0, 1));
      CHECK(Approx(0.36) == my_levelset(0.8, 0, 1));
    }
  }
  SECTION("Centered") {
    my_cyl.center = true;

    auto my_tree = std::make_shared<csg::Tree>();
    my_tree->top.objs.push_back(my_cyl);
    csg::CsgIF my_levelset(my_tree);

    SECTION("Axis") {
      CHECK_FALSE(0 < my_levelset(0, 0, -2));
      CHECK(Approx(-1) == my_levelset(0, 0, -2));

      CHECK(0 < my_levelset(0, 0, 0));
      CHECK(Approx(1) == my_levelset(0, 0, 0));

      CHECK_FALSE(0 < my_levelset(0, 0, 2));
      CHECK(Approx(-1) == my_levelset(0, 0, 2));
    }

    SECTION("Cross section") {
      CHECK_FALSE(0 < my_levelset(-0.8, 0.8, 0));
      CHECK(Approx(-0.28) == my_levelset(-0.8, 0.8, 0));

      CHECK_FALSE(0 < my_levelset(0.8, -0.8, 0));
      CHECK(Approx(-0.28) == my_levelset(0.8, -0.8, 0));

      CHECK(0 < my_levelset(0, 0.8, 0));
      CHECK(Approx(0.36) == my_levelset(0, 0.8, 0));

      CHECK(0 < my_levelset(0.8, 0, 0));
      CHECK(Approx(0.36) == my_levelset(0.8, 0, 0));
    }
  }
}

TEST_CASE("Cube", "[Levelset Primitives]") {

  csg::Cube my_cub;
  my_cub.size = {2, 3, 5};
  my_cub.center = false;

  SECTION("Not Centered") {
    my_cub.center = false;

    auto my_tree = std::make_shared<csg::Tree>();
    my_tree->top.objs.push_back(my_cub);
    csg::CsgIF my_levelset(my_tree);

    SECTION("Outside") {
      CHECK_FALSE(0 < my_levelset(1, 1, -1));
      CHECK(Approx(-1) == my_levelset(1, 1, -1));

      CHECK_FALSE(0 < my_levelset(1, -1, 1));
      CHECK(Approx(-1) == my_levelset(1, -1, 1));

      CHECK_FALSE(0 < my_levelset(-1, 1, 1));
      CHECK(Approx(-1) == my_levelset(-1, 1, 1));

      CHECK_FALSE(0 < my_levelset(3, 1, 1));
      CHECK(Approx(-1) == my_levelset(3, 1, 1));

      CHECK_FALSE(0 < my_levelset(1, 4, 1));
      CHECK(Approx(-1) == my_levelset(1, 4, 1));

      CHECK_FALSE(0 < my_levelset(1, 1, 6));
      CHECK(Approx(-1) == my_levelset(1, 1, 6));
    }
    SECTION("Inside") {
      CHECK(0 < my_levelset(1, 1, 4));
      CHECK(Approx(1) == my_levelset(1, 1, 4));

      CHECK(0 < my_levelset(1, 2, 1));
      CHECK(Approx(1) == my_levelset(1, 2, 1));

      CHECK(0 < my_levelset(1, 1, 1));
      CHECK(Approx(1) == my_levelset(1, 1, 1));
    }
  }

  SECTION("Centered") {
    my_cub.center = true;

    auto my_tree = std::make_shared<csg::Tree>();
    my_tree->top.objs.push_back(my_cub);
    csg::CsgIF my_levelset(my_tree);

    SECTION("Outside") {
      CHECK_FALSE(0 < my_levelset(-0.5, 1, -3));
      CHECK(Approx(-0.5) == my_levelset(-0.5, 1, -3));

      CHECK_FALSE(0 < my_levelset(-0.5, 2, -2));
      CHECK(Approx(-0.5) == my_levelset(-0.5, 2, -2));

      CHECK_FALSE(0 < my_levelset(-1.5, 1, -2));
      CHECK(Approx(-0.5) == my_levelset(-1.5, 1, -2));

      CHECK_FALSE(0 < my_levelset(0.5, -1, 3));
      CHECK(Approx(-0.5) == my_levelset(0.5, -1, 3));

      CHECK_FALSE(0 < my_levelset(0.5, -2, 2));
      CHECK(Approx(-0.5) == my_levelset(0.5, -2, 2));

      CHECK_FALSE(0 < my_levelset(1.5, -1, 2));
      CHECK(Approx(-0.5) == my_levelset(1.5, -1, 2));
    }
    SECTION("Inside") {
      CHECK(0 < my_levelset(0.5, -1, 2));
      CHECK(Approx(0.5) == my_levelset(0.5, -1, 2));

      CHECK(0 < my_levelset(-0.5, 1, -2));
      CHECK(Approx(0.5) == my_levelset(-0.5, 1, -2));
    }
  }
}

TEST_CASE("Sphere", "[Levelset Primitives]") {

  csg::Sphere my_sph;
  my_sph.radius = 10;

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_sph);
  csg::CsgIF my_levelset(my_tree);

  SECTION("Outside") {
    CHECK_FALSE(0 < my_levelset(0, 8, 8));
    CHECK(Approx(-28) == my_levelset(0, 8, 8));

    CHECK_FALSE(0 < my_levelset(-8, 0, 8));
    CHECK(Approx(-28) == my_levelset(-8, 0, 8));

    CHECK_FALSE(0 < my_levelset(0, 8, -8));
    CHECK(Approx(-28) == my_levelset(0, 8, -8));
  }
  SECTION("Inside") {
    CHECK(0 < my_levelset(0, 0, 0));
    CHECK(Approx(100) == my_levelset(0, 0, 0));

    CHECK(0 < my_levelset(0, 7, 7));
    CHECK(Approx(2) == my_levelset(0, 7, 7));

    CHECK(0 < my_levelset(-7, 0, 7));
    CHECK(Approx(2) == my_levelset(-7, 0, 7));

    CHECK(0 < my_levelset(0, 7, -7));
    CHECK(Approx(2) == my_levelset(0, 7, -7));
  }
}
TEST_CASE("Cone", "[Levelset Primitives]") {

  SECTION("radius1 < radius 2; regular)") {

    csg::Cone my_cone;
    my_cone.radius1 = 0;
    my_cone.radius2 = 1;
    my_cone.height = 2;
    my_cone.center = false;

    SECTION("Untruncated") {
      my_cone.radius1 = 0;

      SECTION("Not Centered") {
        my_cone.center = false;

        auto my_tree = std::make_shared<csg::Tree>();
        my_tree->top.objs.push_back(my_cone);
        csg::CsgIF my_levelset(my_tree);

        SECTION("Outside") {
          CHECK_FALSE(0 < my_levelset(0.15, -0.15, 0.1));
          CHECK(
              Approx((0.5 * 0.1) * (0.5 * 0.1) - (0.15 * 0.15 + 0.15 * .15)) ==
              my_levelset(0.15, -0.15, 0.1));

          CHECK_FALSE(0 < my_levelset(-0.15, 0.15, 0.1));
          CHECK(
              Approx((0.5 * 0.1) * (0.5 * 0.1) - (0.15 * 0.15 + 0.15 * .15)) ==
              my_levelset(-0.15, 0.15, 0.1));

          CHECK_FALSE(0 < my_levelset(0, 0.99, 0.9));
          CHECK(Approx((0.5 * 0.9) * (0.5 * 0.9) - 0.99 * 0.99) ==
                my_levelset(0, 0.99, 0.9));

          CHECK_FALSE(0 < my_levelset(-0.99, 0, 0.9));
          CHECK(Approx((0.5 * 0.9) * (0.5 * 0.9) - 0.99 * 0.99) ==
                my_levelset(-0.99, 0, 0.9));

          CHECK_FALSE(0 < my_levelset(0, 0, -0.1));
          CHECK(Approx(-0.1) == my_levelset(0, 0, -0.1));

          CHECK_FALSE(0 < my_levelset(0, 0, 2.1));
          CHECK(Approx(-0.1) == my_levelset(0, 0, 2.1));
        }

        SECTION("Inside") {
          CHECK(0 < my_levelset(0, 0, 0.1));
          CHECK(Approx((0.5 * 0.1) * (0.5 * 0.1)) == my_levelset(0, 0, 0.1));

          CHECK(0 < my_levelset(0, 0, 1.9));
          CHECK(Approx(0.1) == my_levelset(0, 0, 1.9));

          CHECK(0 < my_levelset(0.01, 0, 0.1));
          CHECK(Approx((0.5 * 0.1) * (0.5 * 0.1) - 0.01 * 0.01) ==
                my_levelset(0.01, 0, 0.1));

          CHECK(0 < my_levelset(0, -0.01, 0.1));
          CHECK(Approx((0.5 * 0.1) * (0.5 * 0.1) - 0.01 * 0.01) ==
                my_levelset(0, -0.01, 0.1));

          CHECK(0 < my_levelset(0, 0.4, 0.99));
          CHECK(Approx((0.5 * 0.99) * (0.5 * 0.99) - 0.4 * 0.4) ==
                my_levelset(0, 0.4, 0.99));

          CHECK(0 < my_levelset(-0.4, 0, 0.9));
          CHECK(Approx((0.5 * 0.9) * (0.5 * 0.9) - 0.4 * 0.4) ==
                my_levelset(-0.4, 0, 0.9));
        }
      }

      SECTION("Centered") {
        my_cone.center = true;

        auto my_tree = std::make_shared<csg::Tree>();
        my_tree->top.objs.push_back(my_cone);
        csg::CsgIF my_levelset(my_tree);

        SECTION("Outside") {
          CHECK_FALSE(0 < my_levelset(0.15, -0.15, -0.9));
          CHECK(
              Approx((0.5 * 0.1) * (0.5 * 0.1) - (0.15 * 0.15 + 0.15 * 0.15)) ==
              my_levelset(0.15, -0.15, -0.9));

          CHECK_FALSE(0 < my_levelset(-0.15, 0.15, -0.9));
          CHECK(
              Approx((0.5 * 0.1) * (0.5 * 0.1) - (0.15 * 0.15 + 0.15 * 0.15)) ==
              my_levelset(-0.15, 0.15, -0.9));

          CHECK_FALSE(0 < my_levelset(0, 0.99, -0.1));
          CHECK(Approx((0.5 * 0.9) * (0.5 * 0.9) - 0.99 * 0.99) ==
                my_levelset(0, 0.99, -0.1));

          CHECK_FALSE(0 < my_levelset(-0.99, 0, -0.1));
          CHECK(Approx((0.5 * 0.9) * (0.5 * 0.9) - 0.99 * 0.99) ==
                my_levelset(-0.99, 0, -0.1));

          CHECK_FALSE(0 < my_levelset(0, 0, -1.1));
          CHECK(Approx(-0.1) == my_levelset(0, 0, -1.1));

          CHECK_FALSE(0 < my_levelset(0, 0, 1.1));
          CHECK(Approx(-0.1) == my_levelset(0, 0, 1.1));
        }

        SECTION("Inside") {
          CHECK(0 < my_levelset(0, 0, -0.9));
          CHECK(Approx((0.5 * 0.1) * (0.5 * 0.1)) == my_levelset(0, 0, -0.9));

          CHECK(0 < my_levelset(0, 0, 0.9));
          CHECK(Approx(0.1) == my_levelset(0, 0, 0.9));

          CHECK(0 < my_levelset(0.01, 0, -0.9));
          CHECK(Approx((0.5 * 0.1) * (0.5 * 0.1) - 0.01 * 0.01) ==
                my_levelset(0.01, 0, -0.9));

          CHECK(0 < my_levelset(0, -0.01, -0.9));
          CHECK(Approx((0.5 * 0.1) * (0.5 * 0.1) - 0.01 * 0.01) ==
                my_levelset(0, -0.01, -0.9));

          CHECK(0 < my_levelset(0, 0.4, -0.01));
          CHECK(Approx((0.5 * 0.99) * (0.5 * 0.99) - 0.4 * 0.4) ==
                my_levelset(0, 0.4, -0.01));

          CHECK(0 < my_levelset(-0.4, 0, -0.1));
          CHECK(Approx((0.5 * 0.9) * (0.5 * 0.9) - 0.4 * 0.4) ==
                my_levelset(-0.4, 0, -0.1));
        }
      }
    }

    SECTION("Truncated") {
      my_cone.radius1 = 0.6;

      SECTION("Not Centered") {
        my_cone.center = false;

        auto my_tree = std::make_shared<csg::Tree>();
        my_tree->top.objs.push_back(my_cone);
        csg::CsgIF my_levelset(my_tree);

        SECTION("smaller end") {
          CHECK_FALSE(0 < my_levelset(0, 0, -0.01));
          CHECK(Approx(-0.01) == my_levelset(0, 0, -0.01));

          CHECK(0 < my_levelset(0, 0, 0.01));
          CHECK(Approx(0.01) == my_levelset(0, 0, 0.01));
        }
        SECTION("curved surface near smaller end") {
          CHECK_FALSE(0 < my_levelset(0.61, 0, 0.01));
          CHECK(Approx((0.6 + 0.2 * 0.01) * (0.6 + 0.2 * 0.01) - 0.61 * 0.61) ==
                my_levelset(0.61, 0, 0.01));

          CHECK(0 < my_levelset(0.59, 0, 0.01));
          CHECK(Approx(0.01) == my_levelset(0.59, 0, 0.01));

          CHECK_FALSE(0 < my_levelset(-0.5, 0.5, 0.4));
          CHECK(Approx((0.6 + 0.2 * 0.4) * (0.6 + 0.2 * 0.4) -
                       (0.5 * 0.5 + 0.5 * 0.5)) == my_levelset(-0.5, 0.5, 0.4));

          CHECK(0 < my_levelset(-0.3, 0.5, 0.4));
          CHECK(Approx((0.6 + 0.2 * 0.4) * (0.6 + 0.2 * 0.4) -
                       (0.3 * 0.3 + 0.5 * 0.5)) == my_levelset(-0.3, 0.5, 0.4));
        }
        SECTION("curved surface near the middle") {
          CHECK_FALSE(0 < my_levelset(0.73, -0.3, 0.8));
          CHECK(Approx((0.6 + 0.2 * 0.8) * (0.6 + 0.2 * 0.8) -
                       (0.73 * 0.73 + 0.3 * 0.3)) ==
                my_levelset(0.73, -0.3, 0.8));

          CHECK(0 < my_levelset(0.63, -0.3, 0.8));
          CHECK(Approx((0.6 + 0.2 * 0.8) * (0.6 + 0.2 * 0.8) -
                       (0.63 * 0.63 + 0.3 * 0.3)) ==
                my_levelset(0.63, -0.3, 0.8));

          CHECK_FALSE(0 < my_levelset(-0.9, 0.3, 1.4));
          CHECK(Approx((0.6 + 0.2 * 1.4) * (0.6 + 0.2 * 1.4) -
                       (0.9 * 0.9 + 0.3 * 0.3)) == my_levelset(-0.9, 0.3, 1.4));

          CHECK(0 < my_levelset(-0.8, 0.3, 1.4));
          CHECK(Approx((0.6 + 0.2 * 1.4) * (0.6 + 0.2 * 1.4) -
                       (0.8 * 0.8 + 0.3 * 0.3)) == my_levelset(-0.8, 0.3, 1.4));
        }
        SECTION("larger end") {
          CHECK_FALSE(0 < my_levelset(0, 0, 2.01));
          CHECK(Approx(-0.01) == my_levelset(0, 0, 2.01));

          CHECK(0 < my_levelset(0, 0, 1.99));
          CHECK(Approx(0.01) == my_levelset(0, 0, 1.99));
        }
        SECTION("curved surface near larger end") {
          CHECK_FALSE(0 < my_levelset(1.1, 0, 1.9));
          CHECK(Approx((0.6 + 0.2 * 1.9) * (0.6 + 0.2 * 1.9) - 1.1 * 1.1) ==
                my_levelset(1.1, 0, 1.9));

          CHECK(0 < my_levelset(0.9, 0, 1.9));
          CHECK(Approx(0.1) == my_levelset(0.9, 0, 1.9));

          CHECK_FALSE(0 < my_levelset(1, -0.3, 1.9));
          CHECK(Approx((0.6 + 0.2 * 1.9) * (0.6 + 0.2 * 1.9) -
                       (0.3 * 0.3 + 1 * 1)) == my_levelset(1, -0.3, 1.9));

          CHECK(0 < my_levelset(0.9, -0.3, 1.9));
          CHECK(Approx((0.6 + 0.2 * 1.9) * (0.6 + 0.2 * 1.9) -
                       (0.3 * 0.3 + 0.9 * 0.9)) == my_levelset(0.9, -0.3, 1.9));
        }
      }

      SECTION("Centered") {
        my_cone.center = true;

        auto my_tree = std::make_shared<csg::Tree>();
        my_tree->top.objs.push_back(my_cone);
        csg::CsgIF my_levelset(my_tree);

        SECTION("smaller end") {
          CHECK_FALSE(0 < my_levelset(0, 0, -1.01));
          CHECK(Approx(-0.01) == my_levelset(0, 0, -1.01));

          CHECK(0 < my_levelset(0, 0, -0.99));
          CHECK(Approx(0.01) == my_levelset(0, 0, -0.99));

          CHECK(0 < my_levelset(0, 0, -0.99));
          CHECK(Approx(0.01) == my_levelset(0, 0, -0.99));
        }
        SECTION("curved surface near smaller end") {
          CHECK_FALSE(0 < my_levelset(0.61, 0, -0.99));
          CHECK(Approx((0.6 + 0.2 * 0.01) * (0.6 + 0.2 * 0.01) - 0.61 * 0.61) ==
                my_levelset(0.61, 0, -0.99));

          CHECK(0 < my_levelset(0.5, 0, -0.99));
          CHECK(Approx(0.01) == my_levelset(0.5, 0, -0.99));

          CHECK_FALSE(0 < my_levelset(0.65, -0.1, -0.9));
          CHECK(Approx((0.6 + 0.2 * 0.1) * (0.6 + 0.2 * 0.1) -
                       (0.65 * 0.65 + 0.1 * 0.1)) ==
                my_levelset(0.65, -0.1, -0.9));

          CHECK(0 < my_levelset(0.57, -0.1, -0.9));
          CHECK(Approx((0.6 + 0.2 * 0.1) * (0.6 + 0.2 * 0.1) -
                       (0.57 * 0.57 + 0.1 * 0.1)) ==
                my_levelset(0.57, -0.1, -0.9));
        }
        SECTION("curved surface near middle") {
          CHECK_FALSE(0 < my_levelset(-0.9, 0.3, 0.2));
          CHECK(Approx((0.6 + 0.2 * 1.2) * (0.6 + 0.2 * 1.2) -
                       (0.9 * 0.9 + 0.3 * 0.3)) == my_levelset(-0.9, 0.3, 0.2));

          CHECK(0 < my_levelset(-0.7, 0.3, 0.2));
          CHECK(Approx((0.6 + 0.2 * 1.2) * (0.6 + 0.2 * 1.2) -
                       (0.7 * 0.7 + 0.3 * 0.3)) == my_levelset(-0.7, 0.3, 0.2));

          CHECK_FALSE(0 < my_levelset(0.75, -0.3, -0.2));
          CHECK(Approx((0.6 + 0.2 * 0.8) * (0.6 + 0.2 * 0.8) -
                       (0.75 * 0.75 + 0.3 * 0.3)) ==
                my_levelset(0.75, -0.3, -0.2));

          CHECK(0 < my_levelset(0.6, -0.3, -0.2));
          CHECK(Approx((0.6 + 0.2 * 0.8) * (0.6 + 0.2 * 0.8) -
                       (0.6 * 0.6 + 0.3 * 0.3)) ==
                my_levelset(0.6, -0.3, -0.2));
        }
        SECTION("larger end") {
          CHECK_FALSE(0 < my_levelset(0, 0, 1.01));
          CHECK(Approx(-0.01) == my_levelset(0, 0, 1.01));

          CHECK(0 < my_levelset(0, 0, 0.99));
          CHECK(Approx(0.01) == my_levelset(0, 0, 0.99));
        }
        SECTION("curved surface near the larger end") {
          CHECK_FALSE(0 < my_levelset(1.1, 0, 0.95));
          CHECK(Approx((0.6 + 0.2 * 1.95) * (0.6 + 0.2 * 1.95) - 1.1 * 1.1) ==
                my_levelset(1.1, 0, 0.95));

          CHECK(0 < my_levelset(0.8, 0, 0.95));
          CHECK(Approx(0.05) == my_levelset(0.8, 0, 0.95));

          CHECK_FALSE(0 < my_levelset(1.0, -0.3, 0.8));
          CHECK(Approx((0.6 + 0.2 * 1.8) * (0.6 + 0.2 * 1.8) -
                       (0.3 * 0.3 + 1 * 1)) == my_levelset(1.0, -0.3, 0.8));

          CHECK(0 < my_levelset(0.73, -0.3, 0.8));
          CHECK(Approx(0.2) == my_levelset(0.73, -0.3, 0.8));
        }
      }
    }
  }

  SECTION("radius1 > radius2; inverted") {
    csg::Cone my_cone;
    my_cone.radius1 = 1;
    my_cone.radius2 = 0;
    my_cone.height = 2;
    my_cone.center = false;

    SECTION("Untruncated") {
      my_cone.radius2 = 0;

      SECTION("Not Centered") {
        my_cone.center = false;

        auto my_tree = std::make_shared<csg::Tree>();
        my_tree->top.objs.push_back(my_cone);
        csg::CsgIF my_levelset(my_tree);

        SECTION("Outside") {
          CHECK_FALSE(0 < my_levelset(0.15, -0.15, 1.9));
          CHECK(Approx((1 - 0.5 * 1.9) * (1 - 0.5 * 1.9) -
                       (0.15 * 0.15 + 0.15 * .15)) ==
                my_levelset(0.15, -0.15, 1.9));

          CHECK_FALSE(0 < my_levelset(-0.15, 0.15, 1.9));
          CHECK(Approx((1 - 0.5 * 1.9) * (1 - 0.5 * 1.9) -
                       (0.15 * 0.15 + 0.15 * .15)) ==
                my_levelset(-0.15, 0.15, 1.9));

          CHECK_FALSE(0 < my_levelset(0, 0.99, 1.1));
          CHECK(Approx((1 - 0.5 * 1.1) * (1 - 0.5 * 1.1) - 0.99 * 0.99) ==
                my_levelset(0, 0.99, 1.1));

          CHECK_FALSE(0 < my_levelset(-0.99, 0, 1.1));
          CHECK(Approx((1 - 0.5 * 1.1) * (1 - 0.5 * 1.1) - 0.99 * 0.99) ==
                my_levelset(-0.99, 0, 1.1));

          CHECK_FALSE(0 < my_levelset(0, 0, -0.1));
          CHECK(Approx(-0.1) == my_levelset(0, 0, -0.1));

          CHECK_FALSE(0 < my_levelset(0, 0, 2.1));
          CHECK(Approx(-0.1) == my_levelset(0, 0, 2.1));
        }

        SECTION("Inside") {
          CHECK(0 < my_levelset(0, 0, 1.9));
          CHECK(Approx((1 - 0.5 * 1.9) * (1 - 0.5 * 1.9)) ==
                my_levelset(0, 0, 1.9));

          CHECK(0 < my_levelset(0, 0, 0.1));
          CHECK(Approx(0.1) == my_levelset(0, 0, 0.1));

          CHECK(0 < my_levelset(0.01, 0, 1.9));
          CHECK(Approx((1 - 0.5 * 1.9) * (1 - 0.5 * 1.9) - 0.01 * 0.01) ==
                my_levelset(0.01, 0, 1.9));

          CHECK(0 < my_levelset(0, -0.01, 1.9));
          CHECK(Approx((1 - 0.5 * 1.9) * (1 - 0.5 * 1.9) - 0.01 * 0.01) ==
                my_levelset(0, -0.01, 1.9));

          CHECK(0 < my_levelset(0, 0.4, 1.01));
          CHECK(Approx((1 - 0.5 * 1.01) * (1 - 0.5 * 1.01) - 0.4 * 0.4) ==
                my_levelset(0, 0.4, 1.01));

          CHECK(0 < my_levelset(-0.4, 0, 1.1));
          CHECK(Approx((1 - 0.5 * 1.1) * (1 - 0.5 * 1.1) - 0.4 * 0.4) ==
                my_levelset(-0.4, 0, 1.1));
        }
      }
    }
  }
}

#if USE_CGAL
TEST_CASE("Polyhedron", "[Levelset Primitives]") {
  double Lx = 10.0, Ly = 7.0, Lz = 5.0;
  // A 10 x 7 x 5 cuboid in positive quandrant formed using polyhedron
  csg::Polyhedron my_cub({{0, 0, 0},
                          {Lx, 0, 0},
                          {Lx, Ly, 0},
                          {0, Ly, 0},
                          {0, 0, Lz},
                          {Lx, 0, Lz},
                          {Lx, Ly, Lz},
                          {0, Ly, Lz}},
                         {{0, 1, 2, 3},
                          {4, 5, 1, 0},
                          {7, 6, 5, 4},
                          {5, 6, 2, 1},
                          {6, 7, 3, 2},
                          {7, 4, 0, 3}});

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_cub);
  csg::CsgIF my_levelset(my_tree);

  SECTION("Outside") {
    CHECK_FALSE(0 < my_levelset(1.01 * Lx, 0.1 * Ly, 0.1 * Lz));
    CHECK_FALSE(0 < my_levelset(-0.01 * Lx, 0.5 * Ly, 0.5 * Lz));
    CHECK_FALSE(0 < my_levelset(0.3 * Lx, 1.01 * Ly, 0.2 * Lz));
    CHECK_FALSE(0 < my_levelset(0.7 * Lx, -0.01 * Ly, 0.7 * Lz));
    CHECK_FALSE(0 < my_levelset(0.3 * Lx, 0.9 * Ly, 1.01 * Lz));
    CHECK_FALSE(0 < my_levelset(0.7 * Lx, 0.1 * Ly, -0.01 * Lz));
  }
  SECTION("Inside") {
    CHECK(0 < my_levelset(0.99 * Lx, 0.1 * Ly, 0.1 * Lz));
    CHECK(0 < my_levelset(0.01 * Lx, 0.5 * Ly, 0.5 * Lz));
    CHECK(0 < my_levelset(0.3 * Lx, 0.99 * Ly, 0.2 * Lz));
    CHECK(0 < my_levelset(0.7 * Lx, 0.01 * Ly, 0.7 * Lz));
    CHECK(0 < my_levelset(0.3 * Lx, 0.9 * Ly, 0.99 * Lz));
    CHECK(0 < my_levelset(0.7 * Lx, 0.1 * Ly, 0.01 * Lz));
  }
}
#endif

} // namespace
