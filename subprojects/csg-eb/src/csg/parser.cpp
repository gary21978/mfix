// standard includes
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>

// subproject includes
#include <tao/pegtl.hpp>

// includes
#include "csg_exception.hpp"
#include "csg_types.hpp"

using namespace tao::pegtl;

namespace {

const std::string UNDEFINED_STR = "undef";

using Attr = std::variant<double, bool, std::vector<double>, std::string,
                          std::vector<std::vector<double>>>;
using AttrMap = std::map<std::string, Attr>;

struct parser_state {
  std::vector<std::vector<csg::Type3D>> current_3d_objs;
  std::vector<std::vector<csg::Type2D>> current_2d_objs;
  std::vector<csg::Type3D> current_3d_group;
  std::vector<csg::Type2D> current_2d_group;
  std::string current_name;
  std::vector<double> current_vec;
  std::vector<std::vector<double>> current_matrix;
  std::vector<std::vector<std::vector<double>>> current_matrices;
  AttrMap curr_attr;
  std::vector<AttrMap> curr_attrs;
};

enum Dimension { D2, D3 };

} // namespace

namespace {

struct escaped_x : seq<one<'x'>, rep<2, must<xdigit>>> {};
struct escaped_u : seq<one<'u'>, rep<4, must<xdigit>>> {};
struct escaped_U : seq<one<'U'>, rep<8, must<xdigit>>> {};
struct escaped_c
    : one<'\'', '"', '?', '\\', 'a', 'b', 'f', 'n', 'r', 't', 'v'> {};

struct escaped : sor<escaped_x, escaped_u, escaped_U, escaped_c> {};

struct character
    : if_must_else<one<'\\'>, escaped, utf8::range<0x20, 0x10FFFF>> {};
struct string_literal : if_must<one<'"'>, until<one<'"'>, character>> {};

struct plus_minus : opt<one<'+', '-'>> {};
struct dot : one<'.'> {};

struct decimal
    : if_then_else<dot, plus<digit>, seq<plus<digit>, opt<dot, star<digit>>>> {
};
struct e : one<'e', 'E'> {};
struct p : one<'p', 'P'> {};
struct exponent : seq<plus_minus, plus<digit>> {};
struct double_ : seq<plus_minus, decimal, opt<e, exponent>> {};

struct padded_double : pad<double_, space> {};

struct L_FUN : pad<string<'('>, space> {};
struct R_FUN : pad<string<')'>, space> {};
struct L_ARR : pad<string<'['>, space> {};
struct R_ARR : pad<string<']'>, space> {};
template <Dimension dim> struct L_BLK : pad<string<'{'>, space> {};
template <Dimension dim> struct R_BLK : pad<string<'}'>, space> {};
struct S_CLN : pad<string<';'>, space> {};

struct true_literal : string<'t', 'r', 'u', 'e'> {};
struct false_literal : string<'f', 'a', 'l', 's', 'e'> {};
struct undef_literal : string<'u', 'n', 'd', 'e', 'f'> {};

struct boolean_literal : sor<true_literal, false_literal> {};

struct special_identifier : seq<string<'$'>, identifier> {};

struct name : sor<identifier, special_identifier> {};

struct vector_cell : padded_double {};

struct vector : seq<L_ARR, list<vector_cell, one<','>>, R_ARR> {};

struct vector_attr : vector {};

struct row : vector {};

struct matrix : seq<L_ARR, list<row, one<','>>, R_ARR> {};

struct matrix_attr : matrix {};

struct value : sor<padded_double, boolean_literal, string_literal, vector_attr,
                   matrix_attr, undef_literal> {};

struct keyval : seq<pad<name, space>, string<'='>, pad<value, space>> {};

struct attr_list : list<keyval, one<','>> {};

struct sphere
    : seq<string<'s', 'p', 'h', 'e', 'r', 'e'>, L_FUN, attr_list, R_FUN> {};

struct cube : seq<string<'c', 'u', 'b', 'e'>, L_FUN, attr_list, R_FUN> {};

struct cylinder : seq<string<'c', 'y', 'l', 'i', 'n', 'd', 'e', 'r'>, L_FUN,
                      attr_list, R_FUN> {};

struct polyhedron
    : seq<string<'p', 'o', 'l', 'y', 'h', 'e', 'd', 'r', 'o', 'n'>, L_FUN,
          attr_list, R_FUN> {};

struct circle
    : seq<string<'c', 'i', 'r', 'c', 'l', 'e'>, L_FUN, attr_list, R_FUN> {};

struct square
    : seq<string<'s', 'q', 'u', 'a', 'r', 'e'>, L_FUN, attr_list, R_FUN> {};

struct polygon
    : seq<string<'p', 'o', 'l', 'y', 'g', 'o', 'n'>, L_FUN, attr_list, R_FUN> {
};

template <Dimension dim> struct shape;

template <>
struct shape<Dimension::D3>
    : seq<sor<cube, cylinder, sphere, polyhedron>, opt<S_CLN>> {};

template <>
struct shape<Dimension::D2> : seq<sor<circle, square, polygon>, opt<S_CLN>> {};

template <Dimension dim> struct obj_list;

template <Dimension dim>
struct bool_union
    : seq<string<'u', 'n', 'i', 'o', 'n'>, L_FUN, R_FUN,
          opt<seq<L_BLK<dim>, obj_list<dim>, R_BLK<dim>>>, opt<S_CLN>> {};

template <Dimension dim>
struct bool_intersection
    : seq<string<'i', 'n', 't', 'e', 'r', 's', 'e', 'c', 't', 'i', 'o', 'n'>,
          L_FUN, R_FUN, opt<seq<L_BLK<dim>, obj_list<dim>, R_BLK<dim>>>,
          opt<S_CLN>> {};

template <Dimension dim>
struct bool_diff
    : seq<string<'d', 'i', 'f', 'f', 'e', 'r', 'e', 'n', 'c', 'e'>, L_FUN,
          R_FUN, opt<seq<L_BLK<dim>, obj_list<dim>, R_BLK<dim>>>, opt<S_CLN>> {
};

template <Dimension dim>
struct bool_exp : sor<bool_union<dim>, bool_intersection<dim>, bool_diff<dim>> {
};

template <Dimension dim>
struct mulmat
    : seq<string<'m', 'u', 'l', 't', 'm', 'a', 't', 'r', 'i', 'x'>, L_FUN,
          matrix, R_FUN, opt<seq<L_BLK<dim>, obj_list<dim>, R_BLK<dim>>>,
          opt<S_CLN>> {};

template <Dimension dim>
struct group
    : seq<string<'g', 'r', 'o', 'u', 'p'>, L_FUN, R_FUN,
          opt<seq<L_BLK<dim>, obj_list<dim>, R_BLK<dim>>>, opt<S_CLN>> {};

struct extrude_rot : seq<string<'r', 'o', 't', 'a', 't', 'e', '_', 'e', 'x',
                                't', 'r', 'u', 'd', 'e'>,
                         L_FUN, attr_list, R_FUN,
                         opt<seq<L_BLK<Dimension::D2>, obj_list<Dimension::D2>,
                                 R_BLK<Dimension::D2>>>,
                         opt<S_CLN>> {};

struct extrude_lin : seq<string<'l', 'i', 'n', 'e', 'a', 'r', '_', 'e', 'x',
                                't', 'r', 'u', 'd', 'e'>,
                         L_FUN, attr_list, R_FUN,
                         opt<seq<L_BLK<Dimension::D2>, obj_list<Dimension::D2>,
                                 R_BLK<Dimension::D2>>>,
                         opt<S_CLN>> {};

struct hull
    : seq<string<'h', 'u', 'l', 'l'>, L_FUN, R_FUN, L_BLK<Dimension::D3>,
          obj_list<Dimension::D3>, R_BLK<Dimension::D3>> {};

struct render
    : seq<string<'r', 'e', 'n', 'd', 'e', 'r'>, L_FUN, attr_list, R_FUN,
          opt<seq<L_BLK<Dimension::D3>, obj_list<Dimension::D3>,
                  R_BLK<Dimension::D3>>>,
          opt<S_CLN>> {};

struct colorgroup : seq<string<'c', 'o', 'l', 'o', 'r'>, L_FUN, vector, R_FUN,
                        opt<seq<L_BLK<Dimension::D3>, obj_list<Dimension::D3>,
                                R_BLK<Dimension::D3>>>,
                        opt<S_CLN>> {};

template <Dimension dim>
struct csg_obj : sor<shape<dim>, mulmat<dim>, bool_exp<dim>, group<dim>,
                     extrude_lin, extrude_rot, hull, render, colorgroup> {};

template <Dimension dim> struct obj_list : plus<pad<csg_obj<dim>, space>> {};

struct grammar : seq<obj_list<Dimension::D3>, eof> {};

template <typename Rule> struct action {};

template <> struct action<L_BLK<Dimension::D3>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    // std::cout << " { " << std::endl;
    std::vector<csg::Type3D> new_3d_group;
    st.current_3d_objs.push_back(new_3d_group);
  }
};

template <> struct action<L_BLK<Dimension::D2>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    std::vector<csg::Type2D> new_2d_group;
    st.current_2d_objs.push_back(new_2d_group);
  }
};

template <> struct action<R_BLK<Dimension::D3>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    // std::cout << " } " << std::endl;
    st.current_3d_group.clear();
    for (auto obj : st.current_3d_objs.back()) {
      st.current_3d_group.push_back(obj);
    }
    st.current_3d_objs.pop_back();
  }
};

template <> struct action<R_BLK<Dimension::D2>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    st.current_2d_group.clear();
    for (auto obj : st.current_2d_objs.back()) {
      st.current_2d_group.push_back(obj);
    }
    st.current_2d_objs.pop_back();
  }
};

template <> struct action<vector_cell> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    double cell;
    ss >> cell;
    st.current_vec.push_back(cell);
  }
};

template <> struct action<matrix> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    st.current_matrices.push_back(st.current_matrix);
    st.current_matrix.clear();
  }
};

template <> struct action<attr_list> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    st.curr_attrs.push_back(st.curr_attr);
    st.curr_attr.clear();
    st.current_name.clear();
  }
};

template <> struct action<row> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    st.current_matrix.push_back(st.current_vec);
    st.current_vec.clear();
  }
};

void add_group_2d(parser_state &st) {
  csg::Union2D group;
  for (const auto &curr_obj : st.current_2d_group) {
    group.objs.push_back(curr_obj);
  }
  st.current_2d_group.clear();
  st.current_2d_objs.back().push_back(group);
}

void add_group_3d(parser_state &st) {
  csg::Union3D group;
  for (const auto &curr_obj : st.current_3d_group) {
    group.objs.push_back(curr_obj);
  }
  st.current_3d_group.clear();
  st.current_3d_objs.back().push_back(group);
}

template <> struct action<group<Dimension::D3>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    add_group_3d(st);
  }
};

template <> struct action<group<Dimension::D2>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    add_group_2d(st);
  }
};

template <> struct action<render> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    add_group_3d(st);
  }
};

template <> struct action<colorgroup> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    add_group_3d(st);
  }
};

template <> struct action<bool_union<Dimension::D3>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    add_group_3d(st);
  }
};

template <> struct action<bool_union<Dimension::D2>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    add_group_2d(st);
  }
};

template <> struct action<bool_intersection<Dimension::D3>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto csg_in = csg::Intersection3D();
    for (const auto &curr_obj : st.current_3d_group) {
      csg_in.objs.push_back(curr_obj);
    }
    st.current_3d_group.clear();
    st.current_3d_objs.back().push_back(csg_in);
  }
};

template <> struct action<bool_intersection<Dimension::D2>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto csg_in = csg::Intersection2D();
    for (const auto &curr_obj : st.current_2d_group) {
      csg_in.objs.push_back(curr_obj);
    }
    st.current_2d_group.clear();
    st.current_2d_objs.back().push_back(csg_in);
  }
};

template <> struct action<bool_diff<Dimension::D3>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto csg_diff = csg::Difference3D();
    for (auto it = st.current_3d_group.begin(); it != st.current_3d_group.end();
         ++it) {
      if (it == st.current_3d_group.begin()) {
        csg_diff.first_obj = std::make_shared<csg::Type3D>(*it);
      } else {
        csg_diff.next_objs.objs.push_back(*it);
      }
    }
    st.current_3d_group.clear();
    st.current_3d_objs.back().push_back(csg_diff);
  }
};

template <> struct action<bool_diff<Dimension::D2>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto csg_diff = csg::Difference2D();
    for (auto it = st.current_2d_group.begin(); it != st.current_2d_group.end();
         ++it) {
      if (it == st.current_2d_group.begin()) {
        csg_diff.first_obj = std::make_shared<csg::Type2D>(*it);
      } else {
        csg_diff.next_objs.objs.push_back(*it);
      }
    }
    st.current_2d_group.clear();
    st.current_2d_objs.back().push_back(csg_diff);
  }
};

template <> struct action<mulmat<Dimension::D2>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto mat = st.current_matrices.back();
    auto mulmat =
        csg::Mulmatrix2D({{mat[0][0], mat[0][1]}}, {{mat[1][0], mat[1][1]}});

    mulmat.translation = {
        mat[0][3],
        mat[1][3],
    };

    for (const auto &curr_obj : st.current_2d_group) {
      mulmat.group.objs.push_back(curr_obj);
    }
    st.current_2d_group.clear();
    st.current_2d_objs.back().push_back(mulmat);
    st.current_matrices.pop_back();
  }
};

template <> struct action<mulmat<Dimension::D3>> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto mat = st.current_matrices.back();
    auto mulmat = csg::Mulmatrix3D({{mat[0][0], mat[0][1], mat[0][2]}},
                                   {{mat[1][0], mat[1][1], mat[1][2]}},
                                   {{mat[2][0], mat[2][1], mat[2][2]}});

    mulmat.translation = {
        mat[0][3],
        mat[1][3],
        mat[2][3],
    };

    for (const auto &curr_obj : st.current_3d_group) {
      mulmat.group.objs.push_back(curr_obj);
    }
    st.current_3d_group.clear();
    st.current_3d_objs.back().push_back(mulmat);
    st.current_matrices.pop_back();
  }
};

template <> struct action<extrude_lin> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto lin_ext = csg::LinearExtrude();
    auto &curr_attr = st.curr_attrs.back();
    lin_ext.height = std::get<double>(curr_attr["height"]);
    lin_ext.center = std::get<bool>(curr_attr["center"]);
    auto scale = std::get<std::vector<double>>(curr_attr["scale"]);
    lin_ext.scale = {scale[0], scale[1]};

    lin_ext.twist = std::nullopt;
    if (curr_attr.count("twist") != 0) {
      lin_ext.twist = std::get<double>(curr_attr["twist"]);
    }

    for (const auto &curr_obj : st.current_2d_group) {
      lin_ext.group.objs.push_back(curr_obj);
    }

    st.current_2d_group.clear();
    st.current_3d_objs.back().push_back(lin_ext);
    st.curr_attrs.pop_back();
  }
};

template <> struct action<extrude_rot> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto rot_ext = csg::RotateExtrude();
    auto &curr_attr = st.curr_attrs.back();
    rot_ext.angle = std::get<double>(curr_attr["angle"]);

    if (rot_ext.angle < 0 || rot_ext.angle > 360) {
      throw csg::Exception("action<extrude_rot>",
                           "angle=" + std::to_string(rot_ext.angle) +
                               "; must be between 0 and 360!");
    }

    for (const auto &curr_obj : st.current_2d_group) {
      rot_ext.group.objs.push_back(curr_obj);
    }

    st.current_3d_group.clear();
    st.current_3d_objs.back().push_back(rot_ext);
    st.curr_attrs.pop_back();
  }
};

template <> struct action<hull> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
#if USE_CGAL
    csg::Hull hull;

    // Only a very limited hull is currently supported
    if (st.current_3d_group.size() == 2) {
      auto &obj0 = st.current_3d_group[0];
      auto &obj1 = st.current_3d_group[1];
      if (std::holds_alternative<csg::Mulmatrix3D>(obj0) and
          std::holds_alternative<csg::Sphere>(obj1)) {
        auto &mm = std::get<csg::Mulmatrix3D>(obj0);
        if (mm.group.objs.size() == 1 and
            std::holds_alternative<csg::Cube>(mm.group.objs[0])) {
          auto &cub = std::get<csg::Cube>(mm.group.objs[0]);
          if (cub.center) {
            hull.cube_size = cub.size;
            hull.cube_center = mm.translation;
            hull.sphere_radius = std::get<csg::Sphere>(obj1).radius;
            hull.generate_polyhedron();
            st.current_3d_objs.back().push_back(hull);
            return;
          }
        }
      }
    }
    st.current_3d_group.clear();

    // If the conditions were not met and throw exception
    std::string except_src = "action<hull>";
    std::string except_msg =
        "hull() support is limited to a multmatrix of a centered cube followed "
        "by a sphere; added for the CLR support";
    throw csg::Exception(except_src, except_msg);
#else
    throw csg::Exception("action<hull>",
                         "Needs to be built with CGAL enabled!");
#endif
  }
};

std::optional<std::string> get_name(const AttrMap &curr_attr) {
  if (curr_attr.count("$name")) {
    auto quoted_name = std::get<std::string>(curr_attr.find("$name")->second);
    auto name = quoted_name.substr(1, quoted_name.length() - 2);
    return name;
  }
  return std::nullopt;
}

template <> struct action<cube> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    csg::Cube cub;
    auto &curr_attr = st.curr_attrs.back();
    auto size = std::get<std::vector<double>>(curr_attr["size"]);
    cub.name = get_name(curr_attr);
    cub.size = {size[0], size[1], size[2]};
    cub.center = std::get<bool>(curr_attr["center"]);

    st.current_3d_objs.back().push_back(cub);
    st.curr_attrs.pop_back();
  }
};

template <> struct action<cylinder> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    auto &curr_attr = st.curr_attrs.back();
    if (curr_attr.count("r")) {
      // proper cylinder
      if (curr_attr.count("r1") || curr_attr.count("r2")) {
        throw csg::Exception(
            "action<cylinder>",
            " ERROR: cannot specify both r and (r1 or r2); ambiguous ");
      }
      csg::Cylinder cyl;
      cyl.name = get_name(curr_attr);
      cyl.center = std::get<bool>(curr_attr["center"]);
      cyl.height = std::get<double>(curr_attr["h"]);
      cyl.radius = std::get<double>(curr_attr["r"]);
      st.current_3d_objs.back().push_back(cyl);
      st.curr_attrs.pop_back();

    } else {
      // conic "cylinder"
      csg::Cone cone;
      cone.name = get_name(curr_attr);
      cone.center = std::get<bool>(curr_attr["center"]);
      cone.height = std::get<double>(curr_attr["h"]);
      cone.radius1 = std::get<double>(curr_attr["r1"]);
      cone.radius2 = std::get<double>(curr_attr["r2"]);
      st.current_3d_objs.back().push_back(cone);
      st.curr_attrs.pop_back();
    }
  }
};

template <> struct action<sphere> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    csg::Sphere sph;
    auto &curr_attr = st.curr_attrs.back();
    sph.name = get_name(curr_attr);
    sph.radius = std::get<double>(curr_attr["r"]);

    st.current_3d_objs.back().push_back(sph);
    st.curr_attrs.pop_back();
  }
};

template <> struct action<circle> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    csg::Circle cir;
    auto &curr_attr = st.curr_attrs.back();
    cir.name = get_name(curr_attr);
    cir.radius = std::get<double>(curr_attr["r"]);

    st.current_2d_objs.back().push_back(cir);
    st.curr_attrs.pop_back();
  }
};

template <> struct action<square> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    csg::Square sq;
    auto &curr_attr = st.curr_attrs.back();
    auto size = std::get<std::vector<double>>(curr_attr["size"]);
    sq.name = get_name(curr_attr);
    sq.size = {size[0], size[1]};
    sq.center = std::get<bool>(curr_attr["center"]);

    st.current_2d_objs.back().push_back(sq);
    st.curr_attrs.pop_back();
  }
};

template <> struct action<polygon> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
#if USE_CGAL
    auto &curr_attr = st.curr_attrs.back();

    auto points_raw =
        std::get<std::vector<std::vector<double>>>(curr_attr["points"]);

    std::vector<std::tuple<double, double>> points;
    for (const auto &pt : points_raw) {
      if (pt.size() != 2)
        throw csg::Exception("action<polygon>", "pt.size() != 2");
      points.push_back({pt[0], pt[1]});
    }

    std::vector<std::vector<double>> paths = {{}};
    if (std::holds_alternative<std::string>(curr_attr["paths"])) {
      auto paths_str = std::get<std::string>(curr_attr["paths"]);
      if (paths_str != UNDEFINED_STR)
        throw csg::Exception("action<polygon>",
                             "paths = " + paths_str + " not recognized");
    } else {
      paths = std::get<std::vector<std::vector<double>>>(curr_attr["paths"]);
    }

    auto diff = csg::Difference2D();
    for (auto it = paths.begin(); it != paths.end(); ++it) {
      if (it == paths.begin()) {
        diff.first_obj =
            std::make_shared<csg::Type2D>(csg::Polygon(points, *it));
      } else {
        diff.next_objs.objs.push_back(csg::Polygon(points, *it));
      }
    }

    st.current_2d_objs.back().push_back(diff);
    st.curr_attrs.pop_back();
#else
    throw csg::Exception("action<polygon>",
                         "Needs to be built with CGAL enabled");
#endif
  }
};

template <> struct action<polyhedron> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
#if USE_CGAL
    auto &curr_attr = st.curr_attrs.back();

    auto points_raw =
        std::get<std::vector<std::vector<double>>>(curr_attr["points"]);
    auto faces = std::get<std::vector<std::vector<double>>>(curr_attr["faces"]);

    std::vector<std::tuple<double, double, double>> points;
    for (const auto &pt : points_raw) {
      if (pt.size() != 3)
        throw csg::Exception("action<polyhedron>", "pt.size() != 3");
      points.push_back({pt[0], pt[1], pt[2]});
    }

    csg::Polyhedron polyh(points, faces);
    polyh.name = get_name(curr_attr);

    st.current_3d_objs.back().push_back(polyh);
    st.curr_attrs.pop_back();
#else
    throw csg::Exception("action<polyhedron>",
                         "Needs to be build with CGAL enabled!");
#endif
  }
};

template <> struct action<name> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    ss >> st.current_name;
  }
};

template <> struct action<string_literal> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    std::string v;
    ss >> v;
    if (!st.current_name.empty())
      st.curr_attr[st.current_name] = v;
  }
};

template <> struct action<double_> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    double v;
    ss >> v;
    if (!st.current_name.empty())
      st.curr_attr[st.current_name] = v;
  }
};

template <> struct action<true_literal> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    if (!st.current_name.empty())
      st.curr_attr[st.current_name] = true;
  }
};

template <> struct action<false_literal> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    if (!st.current_name.empty())
      st.curr_attr[st.current_name] = false;
  }
};

template <> struct action<undef_literal> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    if (!st.current_name.empty())
      st.curr_attr[st.current_name] = UNDEFINED_STR;
  }
};

template <> struct action<vector_attr> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    if (!st.current_name.empty())
      st.curr_attr[st.current_name] = st.current_vec;
    st.current_vec.clear();
  }
};

template <> struct action<matrix_attr> {
  template <typename Input>
  static void apply(const Input &in, parser_state &st) {
    std::stringstream ss(in.string());
    if (!st.current_name.empty())
      st.curr_attr[st.current_name] = st.current_matrix;
    st.current_matrix.clear();
  }
};

std::optional<parser_state> do_parse(std::string str) {
  parser_state st;
  std::vector<csg::Type3D> new_3d_group;
  st.current_3d_objs.push_back(new_3d_group);
  memory_input in(str, "std::cin");
  if (!parse<grammar, action>(in, st)) {
    return std::nullopt;
  }
  return st;
}

} // namespace

namespace csg {

std::shared_ptr<Tree> parse_csg(std::string str) {
  auto maybe_state = do_parse(str);
  if (!maybe_state.has_value()) {
    return nullptr;
  }
  auto st = maybe_state.value();

  assert(st.current_3d_objs.size() == 1);

  Tree tree;
  for (auto obj : st.current_3d_objs.back()) {
    tree.top.objs.push_back(obj);
  }

  if (tree.top.objs.empty())
    throw csg::Exception("parse_csg",
                         "Empty csg file"); // Disallow empty .csg file

  return std::make_shared<Tree>(tree);
}
} // namespace csg
