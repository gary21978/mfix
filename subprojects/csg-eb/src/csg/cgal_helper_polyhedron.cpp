#include "csg_cgal_helper.hpp"
#include "csg_exception.hpp"

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/point_generators_3.h>

#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(6,0,0)
#include <CGAL/AABB_traits_3.h>
#else
#include <CGAL/AABB_traits.h>
#endif

namespace {

typedef cgal_helper::Polyhedron::HalfedgeDS HalfedgeDS;
typedef typename HalfedgeDS::Vertex Vertex;
typedef typename Vertex::Point Point;
typedef CGAL::Side_of_triangle_mesh<cgal_helper::Polyhedron, cgal_helper::CK>
    Point_inside;
typedef CGAL::Creator_uniform_3<double, cgal_helper::CK::Point_3> Creator;

// A modifier creating a polyhedron with the incremental builder.
template <class HDS> class PolyhedronBuilder : public CGAL::Modifier_base<HDS> {
private:
  const std::vector<std::tuple<double, double, double>> &m_points;
  const std::vector<std::vector<unsigned int>> &m_faces;

public:
  PolyhedronBuilder(
      const std::vector<std::tuple<double, double, double>> &points,
      const std::vector<std::vector<unsigned int>> &faces)
      : m_points(points), m_faces(faces) {}

  void operator()(HDS &hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    B.begin_surface(m_points.size(), m_faces.size());

    // Add all the vertices first
    for (const auto &p : m_points) {
      auto [px, py, pz] = p;
      B.add_vertex(Point(px, py, pz));
    }

    // Add facets next
    for (const auto &face : m_faces) {
      B.begin_facet();
      if (!B.test_facet(face.begin(), face.end()))
        throw csg::Exception(
            "Polyhedron builder",
            "Unable to create the polyhedron. test_facet failed. Check input!");

      for (const auto &p_index : face) {
        B.add_vertex_to_facet(p_index);
      }
      B.end_facet();
    }

    B.end_surface();
  }
};

void points_on_a_sphere(std::list<cgal_helper::CK::Point_3> &points, int N,
                        double radius) {
  points.clear();

  double phi = CGAL_PI * (3.0 - std::sqrt(5.0));
  for (int i = 0; i < N; ++i) {
    double y = 1.0 - (i / double(N - 1)) * 2;
    double r = std::sqrt(1.0 - y * y);
    double theta = phi * i;

    double x = std::cos(theta) * r;
    double z = std::sin(theta) * r;
    cgal_helper::CK::Point_3 pt(radius * x, radius * y, radius * z);
    points.push_back(pt);
  }
}

} // namespace

namespace cgal_helper {

std::shared_ptr<Polyhedron>
create_polyhedron(const std::vector<std::tuple<double, double, double>> &points,
                  const std::vector<std::vector<unsigned int>> &faces) {
  auto p = std::make_shared<Polyhedron>();

  // Build incrementally
  PolyhedronBuilder<HalfedgeDS> bp(points, faces);
  p->delegate(bp);
  CGAL_assertion(p->is_valid());

  // Triangulate faces - needed for levelset
  CGAL::Polygon_mesh_processing::triangulate_faces(*p);
  CGAL_assertion(p->is_valid());

  return p;
}

std::shared_ptr<AABBTree>
create_aabb_tree(const std::shared_ptr<Polyhedron> &polyhedron) {
  // Construct AABB tree with a KdTree
  auto tree = std::make_shared<AABBTree>(
      faces(*polyhedron).first, faces(*polyhedron).second, *polyhedron);
  tree->accelerate_distance_queries();
  return tree;
}

bool inside(const std::shared_ptr<AABBTree> &tree, double xx, double yy,
            double zz) {
  cgal_helper::CK::Point_3 pt(xx, yy, zz);
  Point_inside inside_tester(*tree);

  // Determine the side and return true if inside!
  return inside_tester(pt) == CGAL::ON_BOUNDED_SIDE;
}

std::shared_ptr<Polyhedron>
create_polyhedron_from_hull(const std::tuple<double, double, double> &cube_size,
                            const std::array<double, 3> &cube_center,
                            double sphere_radius) {

  // This is currently chosen somewhat arbitrarily
  // based on manual testing and CGAL examples
  // TODO: Is this sufficient?
  const int NUM_SAMPLING_PTS = 343;

  auto p = std::make_shared<Polyhedron>();

  std::list<cgal_helper::CK::Point_3> points;
  points_on_a_sphere(points, NUM_SAMPLING_PTS, sphere_radius);

  std::list<cgal_helper::CK::Point_3> points_c;
  CGAL::points_on_cube_grid_3(1.0, (std::size_t)(NUM_SAMPLING_PTS),
                              std::back_inserter(points_c), Creator());

  auto [sa, sb, sc] = cube_size;

  const cgal_helper::CK::Aff_transformation_3 transf(
      sa / 2, 0.0, 0.0, cube_center[0], 0.0, sb / 2, 0.0, cube_center[1], 0.0,
      0.0, sc / 2, cube_center[2]);
  std::transform(points_c.begin(), points_c.end(), points_c.begin(), transf);

  points.insert(points.end(), points_c.begin(), points_c.end());

  // compute convex hull of non-collinear points
  CGAL::convex_hull_3(points.begin(), points.end(), *p);

  CGAL_assertion(p->is_valid());

  return p;
}

} // namespace cgal_helper
