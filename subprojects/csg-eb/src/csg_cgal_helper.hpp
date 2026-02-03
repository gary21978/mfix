#ifndef CGAL_HELPER_H_
#define CGAL_HELPER_H_

#include <CGAL/version.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <memory>

namespace cgal_helper {

typedef CGAL::Simple_cartesian<double> CK;
typedef CGAL::Polyhedron_3<CK> Polyhedron;
typedef CGAL::Polygon_2<CK> Polygon;
typedef CGAL::AABB_face_graph_triangle_primitive<cgal_helper::Polyhedron>
    Primitive;

#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(6,0,0)
  typedef CGAL::AABB_traits_3<cgal_helper::CK, Primitive> Traits;
#else
  typedef CGAL::AABB_traits<cgal_helper::CK, Primitive> Traits;
#endif

typedef CGAL::AABB_tree<Traits> AABBTree;


std::shared_ptr<Polyhedron>
create_polyhedron(const std::vector<std::tuple<double, double, double>> &points,
                  const std::vector<std::vector<unsigned int>> &faces);

std::shared_ptr<AABBTree>
create_aabb_tree(const std::shared_ptr<Polyhedron> &polyhedron);

Polygon create_polygon(const std::vector<std::tuple<double, double>> &points,
                       const std::vector<unsigned int> &path);

bool inside(const std::shared_ptr<AABBTree> &tree, double xx, double yy,
            double zz);

bool inside(const Polygon &polygon, double xx, double yy);

double abs_max_distance(const Polygon &ploygon, double xx, double yy,
                        bool inside);

std::shared_ptr<Polyhedron>
create_polyhedron_from_hull(const std::tuple<double, double, double> &cube_size,
                  const std::array<double,3> &cube_center,
                  double sphere_radius);

} // namespace cgal_helper

#endif
