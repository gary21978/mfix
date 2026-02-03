#if USE_CGAL
#include "csg_cgal_helper.hpp"
#endif
#include "csg_types.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

const double DEGREES_PER_RADIAN = 45 / std::atan(1);

namespace {

template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template <class... Ts>
// clang-format off
overloaded(Ts...) -> overloaded<Ts...>; // not needed as of C++20
// clang-format on

} // namespace

namespace csg {

double signed_distance_3d(const Cone &, double, double, double);
double signed_distance_3d(const Cube &, double, double, double);
double signed_distance_3d(const Cylinder &, double, double, double);
double signed_distance_3d(const Difference3D &, double, double, double);
double signed_distance_3d(const Intersection3D &, double, double, double);
double signed_distance_3d(const Mulmatrix3D &, double, double, double);
double signed_distance_3d(const Sphere &, double, double, double);
double signed_distance_3d(const Type3D &, double, double, double);
double signed_distance_3d(const LinearExtrude &, double, double, double);
double signed_distance_3d(const RotateExtrude &, double, double, double);
#if USE_CGAL
double signed_distance_3d(const Polyhedron &, double, double, double);
double signed_distance_3d(const Hull &, double, double, double);
#endif

double signed_distance_2d(const Union2D &, double, double);

double signed_distance_3d(const Union3D &group, double xx, double yy,
                          double zz) {
  auto sdist = -std::numeric_limits<double>::max();
  for (const auto &member : group.objs) {
    auto sd = signed_distance_3d(member, xx, yy, zz);
    sdist = std::max(sdist, sd);
  };

  return sdist;
}

double signed_distance_3d(const Intersection3D &group, double xx, double yy,
                          double zz) {
  auto sdist = std::numeric_limits<double>::max();
  for (const auto &member : group.objs) {
    auto sd = signed_distance_3d(member, xx, yy, zz);
    sdist = std::min(sdist, sd);
  };
  return sdist;
}

double signed_distance_3d(const Difference3D &group, double xx, double yy,
                          double zz) {
  if (group.first_obj == nullptr)
    return -std::numeric_limits<double>::max();

  auto sdist = signed_distance_3d(*group.first_obj, xx, yy, zz);

  for (const auto &member : group.next_objs.objs) {
    auto sd = signed_distance_3d(member, xx, yy, zz);
    sdist = std::min(sdist, -sd);
  };
  return sdist;
}

double signed_distance_3d(const Mulmatrix3D &mm, double xx, double yy,
                          double zz) {
  auto XX = xx - mm.translation[0];
  auto YY = yy - mm.translation[1];
  auto ZZ = zz - mm.translation[2];
  auto ri = mm.rotation_inv();
  return signed_distance_3d(mm.group,
                            ri[0][0] * XX + ri[0][1] * YY + ri[0][2] * ZZ,
                            ri[1][0] * XX + ri[1][1] * YY + ri[1][2] * ZZ,
                            ri[2][0] * XX + ri[2][1] * YY + ri[2][2] * ZZ);
}

double signed_distance_3d(const Cone &cone, double xx, double yy, double zz) {
  double ZZ = cone.center ? zz + cone.height / 2 : zz;
  double sign_z = (ZZ >= 0 && ZZ <= cone.height) ? -1.0 : 1.0;
  auto dist_z = std::min(std::fabs(ZZ), std::fabs(ZZ - cone.height));

  auto rr = cone.radius1 + (ZZ * (cone.radius2 - cone.radius1) / cone.height);
  auto delta_r = xx * xx + yy * yy - rr * rr;
  double sign_r = delta_r <= 0 ? -1.0 : 1.0;
  auto dist_r = std::fabs(delta_r);

  return EXTERNAL_FLOW * std::max(sign_z * dist_z, sign_r * dist_r);
}

double signed_distance_3d(const Cube &cube, double xx, double yy, double zz) {

  auto [Lx, Ly, Lz] = cube.size;
  auto XX = cube.center ? xx + Lx / 2 : xx;
  double sign_x = (XX >= 0 && XX <= Lx) ? -1.0 : 1.0;
  auto dist_x = std::min(std::fabs(XX), std::fabs(XX - Lx));

  auto YY = cube.center ? yy + Ly / 2 : yy;
  double sign_y = (YY >= 0 && YY <= Ly) ? -1.0 : 1.0;
  auto dist_y = std::min(std::fabs(YY), std::fabs(YY - Ly));

  auto ZZ = cube.center ? zz + Lz / 2 : zz;
  double sign_z = (ZZ >= 0 && ZZ <= Lz) ? -1.0 : 1.0;
  auto dist_z = std::min(std::fabs(ZZ), std::fabs(ZZ - Lz));

  return EXTERNAL_FLOW *
         std::max({sign_x * dist_x, sign_y * dist_y, sign_z * dist_z});
}

double signed_distance_3d(const Sphere &sph, double xx, double yy, double zz) {
  auto delta_r = xx * xx + yy * yy + zz * zz - sph.radius * sph.radius;
  double sign = delta_r <= 0 ? -1.0 : 1.0;
  auto dist = std::fabs(delta_r);
  return EXTERNAL_FLOW * sign * dist;
}

double signed_distance_3d(const Cylinder &cyl, double xx, double yy,
                          double zz) {
  auto ZZ = cyl.center ? zz + cyl.height / 2 : zz;
  double sign_z = (ZZ >= 0 && ZZ <= cyl.height) ? -1.0 : 1.0;
  auto dist_z = std::min(std::fabs(ZZ), std::fabs(ZZ - cyl.height));

  auto rr = cyl.radius;
  auto delta_r = xx * xx + yy * yy - rr * rr;
  double sign_r = delta_r <= 0 ? -1.0 : 1.0;
  auto dist_r = std::fabs(delta_r);

  return EXTERNAL_FLOW * std::max(sign_z * dist_z, sign_r * dist_r);
}

double signed_distance_3d(const LinearExtrude &lin_ext, double xx, double yy,
                          double zz) {
  // TODO: Support scale
  double XX = xx;
  double YY = yy;
  double ZZ = lin_ext.center ? zz + lin_ext.height / 2 : zz;
  double sign_z = (ZZ >= 0 && ZZ <= lin_ext.height) ? -1.0 : 1.0;
  auto dist_z = std::min(std::fabs(ZZ), std::fabs(ZZ - lin_ext.height));

  if (lin_ext.twist.has_value()) {
    double tt =
        (lin_ext.twist.value() * ZZ / lin_ext.height) / DEGREES_PER_RADIAN;
    XX = std::hypot(xx, yy) * std::cos(std::atan2(yy, xx) - tt);
    YY = std::hypot(xx, yy) * std::sin(std::atan2(yy, xx) - tt);
  }

  auto sd_2d = signed_distance_2d(lin_ext.group, XX, YY);
  double sign_2d = sd_2d >= 0 ? -1.0 : 1.0;
  auto dist_2d = std::fabs(sd_2d);

  return EXTERNAL_FLOW * std::max(sign_z * dist_z, sign_2d * dist_2d);
}

double signed_distance_3d(const RotateExtrude &rot_ext, double xx, double yy,
                          double zz) {
  double tt = std::atan2(yy, xx) * DEGREES_PER_RADIAN;
  tt = tt >= 0 ? tt : 360 + tt;

  assert(tt >= 0 && tt <= 360); // tt must be between 0 and 360

  double sign_t = (tt <= rot_ext.angle) ? -1.0 : 1.0;
  double dist_t = std::min(std::fabs(rot_ext.angle - tt), (360 - tt));

  auto XX = std::hypot(xx, yy);
  auto YY = zz;
  auto sd_2d = signed_distance_2d(rot_ext.group, XX, YY);
  double sign_2d = sd_2d >= 0 ? -1.0 : 1.0;
  auto dist_2d = std::fabs(sd_2d);

  return EXTERNAL_FLOW * std::max(sign_t * dist_t, sign_2d * dist_2d);
}

#if USE_CGAL
double signed_distance_3d(const Polyhedron &polyhedron, double xx, double yy,
                          double zz) {
  // TODO: support signed distance instead of -1.0/1.0
  return cgal_helper::inside(polyhedron.cgal_aabb_tree(), xx, yy, zz) ? 1.0
                                                                      : -1.0;
}
#endif

#if USE_CGAL
double signed_distance_3d(const Hull &hull, double xx, double yy, double zz) {
  // TODO: support signed distance instead of -1.0/1.0
  return cgal_helper::inside(hull.cgal_aabb_tree(), xx, yy, zz) ? 1.0 : -1.0;
}
#endif

double signed_distance_3d(const Type3D &obj, double xx, double yy, double zz) {

  return std::visit(
      overloaded{
          [xx, yy, zz](auto &&arg) {
            return signed_distance_3d(arg, xx, yy, zz);
          },
      },
      obj);
}

} // namespace csg
