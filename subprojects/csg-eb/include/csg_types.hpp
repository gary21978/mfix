#ifndef CSG_TYPES_H_
#define CSG_TYPES_H_

#if USE_CGAL
#include "csg_cgal_helper.hpp"
#endif
#include "csg_matrix_functions.hpp"

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace csg {

struct Circle {
  std::optional<std::string> name;
  double radius;
};

struct Square {
  std::optional<std::string> name;
  std::tuple<double, double> size;
  bool center;
};

struct Sphere {
  std::optional<std::string> name;
  double radius;
};

struct Cube {
  std::optional<std::string> name;
  std::tuple<double, double, double> size;
  bool center;
};

struct Cylinder {
  std::optional<std::string> name;
  double radius;
  double height;
  bool center;
};

struct Cone {
  std::optional<std::string> name;
  double radius1;
  double radius2;
  double height;
  bool center;
};

#if USE_CGAL
struct Polyhedron {
private:
  std::vector<std::tuple<double, double, double>> m_points;
  std::vector<std::vector<unsigned int>> m_faces;
  std::shared_ptr<cgal_helper::Polyhedron> m_cgal_polyhedron;
  std::shared_ptr<cgal_helper::AABBTree> m_cgal_aabb_tree;

public:
  std::optional<std::string> name;
  Polyhedron(const std::vector<std::tuple<double, double, double>> &points,
             const std::vector<std::vector<double>> &faces)
      : m_points(points) {
    for (const auto &f : faces) {
      m_faces.push_back(std::vector<unsigned int>(f.begin(), f.end()));
    }
    m_cgal_polyhedron = cgal_helper::create_polyhedron(m_points, m_faces);
    m_cgal_aabb_tree = cgal_helper::create_aabb_tree(m_cgal_polyhedron);
  }

  const std::shared_ptr<cgal_helper::Polyhedron> &cgal_polyhedron() const {
    return m_cgal_polyhedron;
  }

  const std::shared_ptr<cgal_helper::AABBTree> &cgal_aabb_tree() const {
    return m_cgal_aabb_tree;
  }
};
#endif

#if USE_CGAL
struct Polygon {
private:
  std::vector<std::tuple<double, double>> m_points;
  std::vector<unsigned int> m_path;
  cgal_helper::Polygon m_cgal_polygon;

public:
  std::optional<std::string> name;
  Polygon(const std::vector<std::tuple<double, double>> &points,
          const std::vector<double> &path)
      : m_points(points) {
    for (const auto &p : path) {
      m_path.push_back((unsigned int)(p));
    }
    m_cgal_polygon = cgal_helper::create_polygon(m_points, m_path);
  }

  const cgal_helper::Polygon &cgal_polygon() const { return m_cgal_polygon; }
};
#endif

enum Dimension { D2, D3 };

template <Dimension dim> struct Mulmatrix;
template <Dimension dim> struct Union;
template <Dimension dim> struct Intersection;
template <Dimension dim> struct Difference;
struct LinearExtrude;
struct RotateExtrude;
struct Hull;

template <Dimension dim> struct TypeHelper;

#if USE_CGAL
template <> struct TypeHelper<Dimension::D2> {
  using Type =
      std::variant<Circle, Square, Union<Dimension::D2>,
                   Intersection<Dimension::D2>, Difference<Dimension::D2>,
                   Mulmatrix<Dimension::D2>, Polygon>;
};
#else
template <> struct TypeHelper<Dimension::D2> {
  using Type =
      std::variant<Circle, Square, Union<Dimension::D2>,
                   Intersection<Dimension::D2>, Difference<Dimension::D2>,
                   Mulmatrix<Dimension::D2>>;
};
#endif

#if USE_CGAL
template <> struct TypeHelper<Dimension::D3> {
  using Type = std::variant<Sphere, Cube, Cylinder, Cone, Union<Dimension::D3>,
                            Intersection<Dimension::D3>,
                            Difference<Dimension::D3>, Mulmatrix<Dimension::D3>,
                            LinearExtrude, RotateExtrude, Polyhedron, Hull>;
};
#else
template <> struct TypeHelper<Dimension::D3> {
  using Type =
      std::variant<Sphere, Cube, Cylinder, Cone, Union<Dimension::D3>,
                   Intersection<Dimension::D3>, Difference<Dimension::D3>,
                   Mulmatrix<Dimension::D3>, LinearExtrude, RotateExtrude>;
};
#endif

template <Dimension dim> struct Union {
  std::vector<typename TypeHelper<dim>::Type> objs;
};

template <Dimension dim> struct Intersection {
  std::vector<typename TypeHelper<dim>::Type> objs;
};

template <Dimension dim> struct Difference {
  std::shared_ptr<typename TypeHelper<dim>::Type> first_obj;
  Union<dim> next_objs;
};

template <> struct Mulmatrix<Dimension::D3> {
private:
  matrix::Mat3d m_rotation;
  matrix::Mat3d m_rotation_inv;

public:
  std::array<double, 3> translation;
  Union<Dimension::D3> group;

  Mulmatrix(const std::array<double, 3> &rot_row0,
            const std::array<double, 3> &rot_row1,
            const std::array<double, 3> &rot_row2) {
    m_rotation[0] = rot_row0;
    m_rotation[1] = rot_row1;
    m_rotation[2] = rot_row2;
    m_rotation_inv = matrix::inverse(m_rotation);
  }

  const matrix::Mat3d &rotation() const { return m_rotation; }

  const matrix::Mat3d &rotation_inv() const { return m_rotation_inv; }
};

template <> struct Mulmatrix<Dimension::D2> {
private:
  matrix::Mat2d m_rotation;
  matrix::Mat2d m_rotation_inv;

public:
  std::array<double, 2> translation;
  Union<Dimension::D2> group;

  Mulmatrix(const std::array<double, 2> &rot_row0,
            const std::array<double, 2> &rot_row1) {
    m_rotation[0] = rot_row0;
    m_rotation[1] = rot_row1;
    m_rotation_inv = matrix::inverse(m_rotation);
  }

  const matrix::Mat2d &rotation() const { return m_rotation; }

  const matrix::Mat2d &rotation_inv() const { return m_rotation_inv; }
};

// TODO: Support twist & slices
struct LinearExtrude {
  double height;
  bool center;
  std::tuple<double, double> scale;
  std::optional<double> twist;
  Union<Dimension::D2> group;
};

struct RotateExtrude {
  double angle;
  Union<Dimension::D2> group;
};

#if USE_CGAL
struct Hull {
private:
  std::shared_ptr<cgal_helper::Polyhedron> m_cgal_polyhedron;
  std::shared_ptr<cgal_helper::AABBTree> m_cgal_aabb_tree;

public:
  std::optional<std::string> name;

  // TODO: Extend support beyond just sphere & cube
  std::tuple<double, double, double> cube_size;
  std::array<double, 3> cube_center;
  double sphere_radius;

  void generate_polyhedron() {
    m_cgal_polyhedron = cgal_helper::create_polyhedron_from_hull(
        cube_size, cube_center, sphere_radius);
    m_cgal_aabb_tree = cgal_helper::create_aabb_tree(m_cgal_polyhedron);
  }

  const std::shared_ptr<cgal_helper::Polyhedron> &cgal_polyhedron() const {
    return m_cgal_polyhedron;
  }

  const std::shared_ptr<cgal_helper::AABBTree> &cgal_aabb_tree() const {
    return m_cgal_aabb_tree;
  }
};
#endif

struct Tree {
  Union<Dimension::D3> top;
};

const double EXTERNAL_FLOW = -1.0;

double signed_distance_3d(const Union<Dimension::D3> &, double, double, double);

// defining some useful aliases
using Mulmatrix3D = Mulmatrix<Dimension::D3>;
using Mulmatrix2D = Mulmatrix<Dimension::D2>;
using Union3D = Union<Dimension::D3>;
using Union2D = Union<Dimension::D2>;
using Intersection3D = Intersection<Dimension::D3>;
using Intersection2D = Intersection<Dimension::D2>;
using Difference3D = Difference<Dimension::D3>;
using Difference2D = Difference<Dimension::D2>;
using Type2D = TypeHelper<Dimension::D2>::Type;
using Type3D = TypeHelper<Dimension::D3>::Type;

} // namespace csg

#endif
