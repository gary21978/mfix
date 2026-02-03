#if USE_CGAL
#include "csg_cgal_helper.hpp"
#endif
#include "csg_types.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template <class... Ts>
// clang-format off
overloaded(Ts...) -> overloaded<Ts...>; // not needed as of C++20
// clang-format on

} // namespace

namespace csg {

double signed_distance_2d(const Square &, double, double);
double signed_distance_2d(const Circle &, double, double);
double signed_distance_2d(const Type2D &, double, double);
double signed_distance_2d(const Difference2D &, double, double);
double signed_distance_2d(const Intersection2D &, double, double);
double signed_distance_2d(const Mulmatrix2D &, double, double);
double signed_distance_2d(const Union2D &, double, double);
#if USE_CGAL
double signed_distance_2d(const Polygon &, double, double);
#endif

double signed_distance_2d(const Union2D &group, double xx, double yy) {
  auto sdist = -std::numeric_limits<double>::max();
  for (const auto &member : group.objs) {
    auto sd = signed_distance_2d(member, xx, yy);
    sdist = std::max(sdist, sd);
  };

  return sdist;
}

double signed_distance_2d(const Intersection2D &group, double xx, double yy) {
  auto sdist = std::numeric_limits<double>::max();
  for (const auto &member : group.objs) {
    auto sd = signed_distance_2d(member, xx, yy);
    sdist = std::min(sdist, sd);
  };
  return sdist;
}

double signed_distance_2d(const Difference2D &group, double xx, double yy) {
  if (group.first_obj == nullptr)
    return -std::numeric_limits<double>::max();

  auto sdist = signed_distance_2d(*group.first_obj, xx, yy);

  for (const auto &member : group.next_objs.objs) {
    auto sd = signed_distance_2d(member, xx, yy);
    sdist = std::min(sdist, -sd);
  };
  return sdist;
}

double signed_distance_2d(const Mulmatrix2D &mm, double xx, double yy) {
  auto XX = xx - mm.translation[0];
  auto YY = yy - mm.translation[1];
  auto ri = mm.rotation_inv();
  return signed_distance_2d(mm.group, ri[0][0] * XX + ri[0][1] * YY,
                            ri[1][0] * XX + ri[1][1] * YY);
}

double signed_distance_2d(const Square &sq, double xx, double yy) {

  auto [Lx, Ly] = sq.size;
  auto XX = sq.center ? xx + Lx / 2 : xx;
  double sign_x = (XX >= 0 && XX <= Lx) ? -1.0 : 1.0;
  auto dist_x = std::min(std::fabs(XX), std::fabs(XX - Lx));

  auto YY = sq.center ? yy + Ly / 2 : yy;
  double sign_y = (YY >= 0 && YY <= Ly) ? -1.0 : 1.0;
  auto dist_y = std::min(std::fabs(YY), std::fabs(YY - Ly));

  return EXTERNAL_FLOW * std::max(sign_x * dist_x, sign_y * dist_y);
}

double signed_distance_2d(const Circle &cir, double xx, double yy) {

  auto delta = xx * xx + yy * yy - cir.radius * cir.radius;
  double sign = delta <= 0 ? -1.0 : 1.0;
  auto dist = std::fabs(delta);

  return EXTERNAL_FLOW * sign * dist;
}

#if USE_CGAL
double signed_distance_2d(const Polygon &polygon, double xx, double yy) {

  bool inside = cgal_helper::inside(polygon.cgal_polygon(), xx, yy);
  double sign = inside ? -1.0 : 1.0;
  double dist =
      cgal_helper::abs_max_distance(polygon.cgal_polygon(), xx, yy, inside);

  return EXTERNAL_FLOW * sign * dist;
}
#endif

double signed_distance_2d(const Type2D &obj, double xx, double yy) {

  return std::visit(
      overloaded{
          [xx, yy](auto &&arg) { return signed_distance_2d(arg, xx, yy); },
      },
      obj);
}

} // namespace csg
