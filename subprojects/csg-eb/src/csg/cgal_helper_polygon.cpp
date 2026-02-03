#include "csg_cgal_helper.hpp"

#include <limits>

namespace {
typedef CGAL::Point_2<cgal_helper::CK> Point;
typedef cgal_helper::CK::Segment_2 Segment;
} // namespace

namespace cgal_helper {

Polygon create_polygon(const std::vector<std::tuple<double, double>> &points,
                       const std::vector<unsigned int> &path) {

  Polygon polyg;

  if (path.empty()) { // Include all points
    for (const auto &p : points) {
      auto [px, py] = p;
      polyg.push_back(Point(px, py));
    }
  } else {
    for (auto p_index : path) {
      auto [px, py] = points[p_index];
      polyg.push_back(Point(px, py));
    }
  }

  return polyg;
}

bool inside(const Polygon &polygon, double xx, double yy) {

  return polygon.has_on_bounded_side(Point(xx, yy));
}

double abs_max_distance(const Polygon &polygon, double xx, double yy,
                        bool inside) {
  double dist = -std::numeric_limits<double>::max();
  auto m = Point(xx, yy);

  double sign = inside ? -1.0 : 1.0;

  for (auto e = polygon.edges_begin(); e != polygon.edges_end(); ++e) {
    auto d = sign * CGAL::squared_distance(*e, m);
    dist = std::max(dist, d);
  }

  return std::fabs(dist);
}

} // namespace cgal_helper
