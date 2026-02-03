#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_matrix_functions.hpp>
#include <csg_types.hpp>

namespace {

TEST_CASE("linear", "[Levelset Extrude]") {
  double height = 100, radius = 10;
  csg::LinearExtrude my_lin_ext;
  my_lin_ext.height = 100;
  my_lin_ext.center = false;
  my_lin_ext.scale = {1, 1};
  my_lin_ext.group = csg::Union2D();

  csg::Circle my_cir;
  my_cir.radius = radius;

  SECTION("Not centered") {
    my_lin_ext.center = false;
    my_lin_ext.group.objs.push_back(my_cir);
    auto my_tree = std::make_shared<csg::Tree>();
    my_tree->top.objs.push_back(my_lin_ext);
    csg::CsgIF my_levelset(my_tree);

    SECTION("Outside") {
      CHECK_FALSE(0 < my_levelset(0, radius * 0.99, -0.01 * height));
      CHECK(Approx(-0.01 * height) ==
            my_levelset(0, radius * 0.99, -0.01 * height));

      CHECK_FALSE(0 < my_levelset(0, radius * 1.01, 0.01 * height));
      CHECK(Approx((radius) * (radius) - (1.01 * radius) * (1.01 * radius)) ==
            my_levelset(0, radius * 1.01, 0.01 * height));

      CHECK_FALSE(0 < my_levelset(0, radius * 0.99, 1.01 * height));
      CHECK(Approx(-0.01 * height) ==
            my_levelset(0, radius * 0.99, 1.01 * height));

      CHECK_FALSE(0 < my_levelset(radius * 0.99, 0, 1.01 * height));
      CHECK(Approx(-0.01 * height) ==
            my_levelset(radius * 0.99, 0, 1.01 * height));
    }
    SECTION("Inside") {
      CHECK(0 < my_levelset(0, radius * 0.99, 0.01 * height));
      CHECK(Approx(0.01 * height) ==
            my_levelset(0, radius * 0.99, 0.01 * height));

      CHECK(0 < my_levelset(radius * 0.99, 0, 0.01 * height));
      CHECK(Approx(0.01 * height) ==
            my_levelset(radius * 0.99, 0, 0.01 * height));

      CHECK(0 < my_levelset(0, radius * 0.99, 0.9 * height));
      CHECK(Approx((radius) * (radius) - (0.99 * radius) * (0.99 * radius)) ==
            my_levelset(0, radius * 0.99, 0.9 * height));

      CHECK(0 < my_levelset(radius * 0.99, 0, 0.9 * height));
      CHECK(Approx((radius) * (radius) - (0.99 * radius) * (0.99 * radius)) ==
            my_levelset(radius * 0.99, 0, 0.9 * height));
    }
  }

  SECTION("Centered") {
    my_lin_ext.center = true;
    my_lin_ext.group.objs.push_back(my_cir);
    auto my_tree = std::make_shared<csg::Tree>();
    my_tree->top.objs.push_back(my_lin_ext);
    csg::CsgIF my_levelset(my_tree);

    SECTION("Outside") {
      CHECK_FALSE(0 < my_levelset(0, radius * 0.99, -1.01 * height / 2));
      CHECK(Approx(-0.01 * height / 2) ==
            my_levelset(0, radius * 0.99, -1.01 * height / 2));

      CHECK_FALSE(0 < my_levelset(0, radius * 1.01, 0.01 * height));
      CHECK(Approx((radius) * (radius) - (1.01 * radius) * (1.01 * radius)) ==
            my_levelset(0, radius * 1.01, 0.01 * height));

      CHECK_FALSE(0 < my_levelset(0, radius * 0.99, 1.01 * height / 2));
      CHECK(Approx(-0.01 * height / 2) ==
            my_levelset(0, radius * 0.99, 1.01 * height / 2));

      CHECK_FALSE(0 < my_levelset(radius * 0.99, 0, 1.01 * height / 2));
      CHECK(Approx(-0.01 * height / 2) ==
            my_levelset(radius * 0.99, 0, 1.01 * height / 2));
    }
    SECTION("Inside") {
      CHECK(0 < my_levelset(0, radius * 0.99, -0.99 * height / 2));
      CHECK(Approx(0.01 * height / 2) ==
            my_levelset(0, radius * 0.99, -0.99 * height / 2));

      CHECK(0 < my_levelset(radius * 0.99, 0, 0.01 * height / 2));
      CHECK(Approx((radius) * (radius) - (0.99 * radius) * (0.99 * radius)) ==
            my_levelset(radius * 0.99, 0, 0.01 * height / 2));

      CHECK(0 < my_levelset(0, radius * 0.99, 0.99 * height / 2));
      CHECK(Approx(0.01 * height / 2) ==
            my_levelset(0, radius * 0.99, 0.99 * height / 2));

      CHECK(0 < my_levelset(radius * 0.99, 0, 0.9 * height / 2));
      CHECK(Approx((radius) * (radius) - (0.99 * radius) * (0.99 * radius)) ==
            my_levelset(radius * 0.99, 0, 0.9 * height / 2));
    }
  }
}

TEST_CASE("two shape linear extrude", "[Levelset Extrude]") {
  csg::LinearExtrude my_lin_ext;
  my_lin_ext.height = 100;
  my_lin_ext.center = true;
  my_lin_ext.scale = {1, 1};
  my_lin_ext.group = csg::Union2D();

  csg::Circle my_cir;
  my_cir.radius = 1;

  my_lin_ext.group.objs.push_back(my_cir);

  auto my_mat = csg::Mulmatrix2D({{1, 0}}, {{0, 1}});
  my_mat.translation = {4, 0};
  my_mat.group.objs.push_back(my_cir);
  my_lin_ext.group.objs.push_back(my_mat);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_lin_ext);
  csg::CsgIF my_levelset(my_tree);

  SECTION("Outside") {
    CHECK_FALSE(0 < my_levelset(0, 2, 0));
    CHECK_FALSE(0 < my_levelset(-4, 0, 0));
    CHECK_FALSE(0 < my_levelset(3.5, 0, 51));
    CHECK_FALSE(0 < my_levelset(3.5, 0, -51));
  }
  SECTION("Inside") {
    CHECK(0 < my_levelset(0, 0, 0));
    CHECK(0 < my_levelset(0.5, 0, 0));
    CHECK(0 < my_levelset(-0.5, 0, 0));
    CHECK(0 < my_levelset(4.5, 0, 0));
    CHECK(0 < my_levelset(3.5, 0, 0));
    CHECK(0 < my_levelset(3.5, 0, 49));
    CHECK(0 < my_levelset(3.5, 0, -49));
  }
}

TEST_CASE("simple torus", "[Levelset Extrude]") {
  csg::RotateExtrude my_rot_ext;
  my_rot_ext.angle = 360;
  my_rot_ext.group = csg::Union2D();

  csg::Circle my_cir;
  my_cir.radius = 1;

  auto my_mat = csg::Mulmatrix2D({{1, 0}}, {{0, 1}});
  my_mat.translation = {2, 0};
  my_mat.group.objs.push_back(my_cir);
  my_rot_ext.group.objs.push_back(my_mat);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_rot_ext);
  csg::CsgIF my_levelset(my_tree);

  SECTION("Outside") {
    CHECK_FALSE(0 < my_levelset(0, 0, 0));
    CHECK_FALSE(0 < my_levelset(0.9, 0, 0));
    CHECK_FALSE(0 < my_levelset(-0.9, 0, 0));
    CHECK_FALSE(0 < my_levelset(2, 0, 1.1));
    CHECK_FALSE(0 < my_levelset(-2, 0, -1.1));
  }
  SECTION("Inside") {
    CHECK(0 < my_levelset(2, 0, 0));
    CHECK(0 < my_levelset(1.1, 0, 0));
    CHECK(0 < my_levelset(-1.1, 0, 0));
    CHECK(0 < my_levelset(2, 0, 0.9));
    CHECK(0 < my_levelset(-2, 0, -0.9));
  }
}

#if USE_CGAL
TEST_CASE("Linear extrude of a triangle with hole", "[Levelset Primitives]") {
  double ht = 100.0, outer = 100.0, inner = 80.0;
  double thick = outer - inner;

  csg::LinearExtrude my_lin_ext;
  my_lin_ext.height = ht;
  my_lin_ext.center = true;
  my_lin_ext.scale = {1, 1};
  my_lin_ext.group = csg::Union2D();

  // A 100 x 100 right triangle with a 80 x 80 hole formed by using a polygon
  // and extrude
  csg::Polygon my_outer_tri({{0, 0},
                             {outer, 0},
                             {0, outer},
                             {thick, thick},
                             {inner, thick},
                             {thick, inner}},
                            {0, 1, 2});

  csg::Polygon my_inner_tri({{0, 0},
                             {outer, 0},
                             {0, outer},
                             {thick, thick},
                             {inner, thick},
                             {thick, inner}},
                            {3, 4, 5});

  csg::Union2D my_union;
  my_union.objs.push_back(my_inner_tri);

  csg::Difference2D my_diff;
  my_diff.first_obj = std::make_shared<csg::Type2D>(my_outer_tri);
  my_diff.next_objs = csg::Union2D(my_union);

  my_lin_ext.group.objs.push_back(my_diff);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_lin_ext);
  csg::CsgIF my_levelset(my_tree);

  SECTION("Outside") {
    CHECK_FALSE(0 < my_levelset(0.5 * thick, 0.5 * thick, 1.01 * (ht / 2)));
    CHECK(Approx(-0.01 * ht / 2) ==
          my_levelset(0.5 * thick, 0.5 * thick, 1.01 * (ht / 2)));

    CHECK_FALSE(0 < my_levelset(0.5 * thick, 0.5 * thick, -1.01 * (ht / 2)));
    CHECK(Approx(-0.01 * ht / 2) ==
          my_levelset(0.5 * thick, 0.5 * thick, -1.01 * (ht / 2)));

    CHECK_FALSE(0 < my_levelset(1.5 * thick, 1.5 * thick, 0.0));
    CHECK(Approx(-(0.5 * thick) * (0.5 * thick)) ==
          my_levelset(1.5 * thick, 1.5 * thick, 0.0));

    CHECK_FALSE(0 < my_levelset(1.01 * outer, 0.5 * thick, 0.0));
    CHECK(Approx(-(1.01 * outer) * (1.01 * outer)) ==
          my_levelset(1.01 * outer, 0.5 * thick, 0.0));

    CHECK_FALSE(0 < my_levelset(0.5 * thick, 1.01 * outer, 0.0));
    CHECK(Approx(-(1.01 * outer) * (1.01 * outer)) ==
          my_levelset(0.5 * thick, 1.01 * outer, 0.0));

    CHECK_FALSE(0 < my_levelset(inner, inner, 0.0));
    CHECK(Approx(-inner * inner) == my_levelset(inner, inner, 0.0));
  }
  SECTION("Inside") {
    CHECK(0 < my_levelset(0.5 * thick, 0.5 * thick, 0.99 * (ht / 2)));
    CHECK(Approx(0.01 * ht / 2) ==
          my_levelset(0.5 * thick, 0.5 * thick, 0.99 * (ht / 2)));

    CHECK(0 < my_levelset(0.5 * thick, 0.5 * thick, -0.99 * (ht / 2)));
    CHECK(Approx(0.01 * ht / 2) ==
          my_levelset(0.5 * thick, 0.5 * thick, -0.99 * (ht / 2)));

    CHECK(0 < my_levelset(0.1 * thick, 0.5 * thick, 0.0));
    CHECK(Approx((0.1 * thick) * (0.1 * thick)) ==
          my_levelset(0.1 * thick, 0.5 * thick, 0.0));

    CHECK(0 < my_levelset(0.7 * outer, 0.1 * thick, 0.0));
    CHECK(Approx((0.1 * thick) * (0.1 * thick)) ==
          my_levelset(0.7 * outer, 0.1 * thick, 0.0));

    CHECK(0 < my_levelset(0.5 * thick, 0.7 * outer, 0.0));
    CHECK(Approx(ht / 2) == my_levelset(0.5 * thick, 0.7 * outer, 0.0));
  }
}
#endif

TEST_CASE("sliced torus", "[Levelset Extrude]") {
  csg::RotateExtrude my_rot_ext;
  my_rot_ext.angle = 45;
  my_rot_ext.group = csg::Union2D();

  csg::Circle my_cir;
  my_cir.radius = 1;

  auto my_mat = csg::Mulmatrix2D({{1, 0}}, {{0, 1}});
  my_mat.translation = {2, 0};
  my_mat.group.objs.push_back(my_cir);
  my_rot_ext.group.objs.push_back(my_mat);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_rot_ext);
  csg::CsgIF my_levelset(my_tree);

  SECTION("Outside") {
    CHECK_FALSE(0 < my_levelset(0, 0, 0));
    CHECK_FALSE(0 < my_levelset(0.9, 0, 0));
    CHECK_FALSE(0 < my_levelset(-0.9, 0, 0));
    CHECK_FALSE(0 < my_levelset(2, 0, 1.1));
    CHECK_FALSE(0 < my_levelset(-2, 0, -1.1));
    CHECK_FALSE(0 < my_levelset(2, -0.01, 0));
    CHECK_FALSE(0 < my_levelset(1.5, 1.6, 0.9));
    CHECK_FALSE(0 < my_levelset(-1.1, 0.0, 0));
    CHECK_FALSE(0 < my_levelset(1.5, 1.6, -0.5));
    CHECK_FALSE(0 < my_levelset(-2, 0.0, -0.9));
  }
  SECTION("Inside") {
    CHECK(0 < my_levelset(2, 0.01, 0));
    CHECK(0 < my_levelset(1.1, 0.01, 0));
    CHECK(0 < my_levelset(2, 0.01, 0.9));
    CHECK(0 < my_levelset(1.5, 1.4, 0.9));
    CHECK(0 < my_levelset(1.5, 1.4, -0.5));
  }
}

TEST_CASE("linear with twist", "[Levelset Extrude]") {
  double height = 10, radius = 1;
  csg::LinearExtrude my_lin_ext;
  my_lin_ext.height = height;
  my_lin_ext.center = false;
  my_lin_ext.scale = {1, 1};
  my_lin_ext.twist = 320;
  my_lin_ext.group = csg::Union2D();

  csg::Circle my_cir;
  my_cir.radius = radius;

  auto my_mat = csg::Mulmatrix2D({{1, 0}}, {{0, 1}});
  my_mat.translation = {2, 0};
  my_mat.group.objs.push_back(my_cir);
  my_lin_ext.group.objs.push_back(my_mat);

  SECTION("Not Centered") {
    my_lin_ext.center = false;
    auto my_tree = std::make_shared<csg::Tree>();
    my_tree->top.objs.push_back(my_lin_ext);
    csg::CsgIF my_levelset(my_tree);

    SECTION("Outside") {
      CHECK_FALSE(0 < my_levelset(-1, 0, height / 2));
      CHECK_FALSE(0 < my_levelset(-2.9, 0, height / 2));
    }

    SECTION("Inside") {
      CHECK(0 < my_levelset(-1.3, 0, height / 2));
      CHECK_FALSE(0 < my_levelset(-2.7, 0, height / 2));
    }
  }

  SECTION("Centered") {
    my_lin_ext.center = true;
    auto my_tree = std::make_shared<csg::Tree>();
    my_tree->top.objs.push_back(my_lin_ext);
    csg::CsgIF my_levelset(my_tree);

    SECTION("Outside") {
      CHECK_FALSE(0 < my_levelset(-1, 0, 0));
      CHECK_FALSE(0 < my_levelset(-2.9, 0, 0));
    }

    SECTION("Inside") {
      CHECK(0 < my_levelset(-1.3, 0, 0));
      CHECK_FALSE(0 < my_levelset(-2.7, 0, 0));
    }
  }
}

TEST_CASE("linear with mulmatrix rotation", "[Levelset Extrude]") {
  double height = 1, size = 2;
  csg::LinearExtrude my_lin_ext;
  my_lin_ext.height = height;
  my_lin_ext.center = false;
  my_lin_ext.scale = {1, 1};
  my_lin_ext.group = csg::Union2D();

  csg::Square my_sq;
  my_sq.size = {2, 2};
  my_sq.center = true;

  // 30 degree rotation
  auto my_mat = csg::Mulmatrix2D({{0.866025, -0.5}}, {{0.5, 0.866025}});

  my_mat.translation = {0, 0};
  my_mat.group.objs.push_back(my_sq);
  my_lin_ext.group.objs.push_back(my_mat);

  auto my_tree = std::make_shared<csg::Tree>();
  my_tree->top.objs.push_back(my_lin_ext);
  csg::CsgIF my_levelset(my_tree);

  SECTION("Outside") {
    CHECK(0 < my_levelset(-1.2 * 0.5 * height, 0, 0.5 * height));
    CHECK(0 < my_levelset(0.5 * height, 0.25 * height, 0.5 * height));
  }

  SECTION("Inside") {
    CHECK(0 < my_levelset(0, 0, 0.5 * height));
    CHECK(0 < my_levelset(-1.1 * 0.5 * height, 0, 0.5 * height));
    CHECK(0 < my_levelset(0.5 * height, 0.5 * height, 0.5 * height));
  }
}

} // namespace
