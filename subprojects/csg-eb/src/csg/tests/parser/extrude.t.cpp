#include "catch2/catch_all.hpp"

using namespace Catch;

#include <csg.hpp>
#include <csg_types.hpp>
#include <iostream>

namespace {
TEST_CASE("linear extrude circle", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(
height = 10,
center = true,
scale = [10,1]
) {
    circle(r = 1);
}
)");
  auto lin_ext = std::get<csg::LinearExtrude>(st.top.objs.back());
  CHECK(lin_ext.group.objs.size() == 1);
  CHECK(lin_ext.height == 10);
  CHECK(lin_ext.center == true);
  CHECK(!lin_ext.twist.has_value());

  auto [Sx, Sy] = lin_ext.scale;
  CHECK(Sx == 10);
  CHECK(Sy == 1);

  auto cir = std::get<csg::Circle>(lin_ext.group.objs.at(0));
  CHECK(cir.radius == 1);
}

TEST_CASE("linear extrude circle with twist", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(
height = 10,
center = true,
scale = [10,1],
twist = -100
) {
    circle(r = 1);
}
)");
  auto lin_ext = std::get<csg::LinearExtrude>(st.top.objs.back());
  CHECK(lin_ext.group.objs.size() == 1);
  CHECK(lin_ext.height == 10);
  CHECK(lin_ext.center == true);
  CHECK(lin_ext.twist == -100);

  auto [Sx, Sy] = lin_ext.scale;
  CHECK(Sx == 10);
  CHECK(Sy == 1);

  auto cir = std::get<csg::Circle>(lin_ext.group.objs.at(0));
  CHECK(cir.radius == 1);
}

TEST_CASE("rotate extrude", "[csg]") {
  auto st = *csg::parse_csg(R"(
rotate_extrude(
angle = 90
) {
    square(size = [2, 3], center = false);
}
)");
  auto rot_ext = std::get<csg::RotateExtrude>(st.top.objs.back());
  CHECK(rot_ext.group.objs.size() == 1);
  CHECK(rot_ext.angle == 90);

  auto sq = std::get<csg::Square>(rot_ext.group.objs.at(0));
  auto [xx, yy] = sq.size;
  CHECK(xx == 2);
  CHECK(yy == 3);
}

TEST_CASE("rotate extrude with bad angle", "[csg]") {
  CHECK_THROWS_WITH(csg::parse_csg(R"(
rotate_extrude(
angle = 370
) {
    square(size = [2, 3], center = false);
}
)"),
                    Matchers::ContainsSubstring("action<extrude_rot>"));
}

TEST_CASE("extrude two shapes", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(
height = 10,
center = true,
scale = [5, 3]
) {
    circle(r = 1);
    square(size = [2, 2], center = false);
}
)");

  auto lin_ext = std::get<csg::LinearExtrude>(st.top.objs.back());
  CHECK(lin_ext.group.objs.size() == 2);
  CHECK(lin_ext.height == 10);
  CHECK(lin_ext.center == true);
  CHECK(!lin_ext.twist.has_value());

  auto [Sx, Sy] = lin_ext.scale;
  CHECK(Sx == 5);
  CHECK(Sy == 3);

  auto cir = std::get<csg::Circle>(lin_ext.group.objs.at(0));
  CHECK(cir.radius == 1);

  auto sq = std::get<csg::Square>(lin_ext.group.objs.at(1));
  auto [xx, yy] = sq.size;
  CHECK(xx == 2);
  CHECK(yy == 2);
}

TEST_CASE("extrude with mulmatrix", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(height = 10, center = false, scale = [5, 3]) {
	circle(r = 2);
	multmatrix([[1, 0, 0, 5], [0, 1, 0, 5], [0, 0, 1, 0], [0, 0, 0, 1]]) {
		square(size = [4, 5], center = false);
	}
}
)");

  auto lin_ext = std::get<csg::LinearExtrude>(st.top.objs.back());
  CHECK(lin_ext.group.objs.size() == 2);
  CHECK(lin_ext.height == 10);
  CHECK(lin_ext.center == false);
  CHECK(!lin_ext.twist.has_value());

  auto [Sx, Sy] = lin_ext.scale;
  CHECK(Sx == 5);
  CHECK(Sy == 3);

  auto cir = std::get<csg::Circle>(lin_ext.group.objs.at(0));
  CHECK(cir.radius == 2);

  auto mm = std::get<csg::Mulmatrix2D>(lin_ext.group.objs.at(1));
  auto [Cx, Cy] = mm.translation;
  CHECK(Cx == 5);
  CHECK(Cy == 5);
  CHECK(mm.group.objs.size() == 1);

  auto sq = std::get<csg::Square>(mm.group.objs.at(0));
  auto [xx, yy] = sq.size;
  CHECK(xx == 4);
  CHECK(yy == 5);
  CHECK(sq.center == false);
}

TEST_CASE("extrude with 3D object", "[csg]") {
  auto st = csg::parse_csg(R"(
linear_extrude(height = 10, center = false, scale = [5, 5]) {
	multmatrix([[1, 0, 0, 5], [0, 1, 0, 5], [0, 0, 1, 0], [0, 0, 0, 1]]) {
		sphere(r = 8);
	}
}
)");

  CHECK(st == nullptr);
}

#if USE_CGAL
TEST_CASE("linear extrude polygon without paths", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(
height = 100, 
center = true, 
scale = [1, 1]) {
polygon(
points = [[0, 0], [100, 0], [130, 50], [30, 50]], 
paths = undef);
}
)");
  auto lin_ext = std::get<csg::LinearExtrude>(st.top.objs.back());
  CHECK(lin_ext.group.objs.size() == 1);
  CHECK(lin_ext.height == 100);
  CHECK(lin_ext.center == true);
  CHECK(!lin_ext.twist.has_value());

  auto [Sx, Sy] = lin_ext.scale;
  CHECK(Sx == 1);
  CHECK(Sy == 1);

  auto diff = std::get<csg::Difference2D>(lin_ext.group.objs.at(0));
  auto outer = std::get<csg::Polygon>(*diff.first_obj);

  CHECK(diff.next_objs.objs.size() == 0);
  CHECK(outer.cgal_polygon().size() == 4);
}
#endif

#if USE_CGAL
TEST_CASE("linear extrude polygon with hole", "[csg]") {
  auto st = *csg::parse_csg(R"(
linear_extrude(
height = 100, 
center = true, 
scale = [1, 1]) {
polygon(
points = [[0, 0], [100, 0], [0, 100], 
   [10, 10], [80, 10], [10, 80]], 
paths = [[0, 1, 2], [3, 4, 5]]);
}
)");
  auto lin_ext = std::get<csg::LinearExtrude>(st.top.objs.back());
  CHECK(lin_ext.group.objs.size() == 1);
  CHECK(lin_ext.height == 100);
  CHECK(lin_ext.center == true);
  CHECK(!lin_ext.twist.has_value());

  auto [Sx, Sy] = lin_ext.scale;
  CHECK(Sx == 1);
  CHECK(Sy == 1);

  auto diff = std::get<csg::Difference2D>(lin_ext.group.objs.at(0));
  auto outer = std::get<csg::Polygon>(*diff.first_obj);
  auto inner = std::get<csg::Polygon>(diff.next_objs.objs.at(0));
  CHECK(outer.cgal_polygon().size() == 3);
  CHECK(inner.cgal_polygon().size() == 3);
}
#endif

} // namespace
