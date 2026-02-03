#include <fstream>
#include <iostream>
#include <sstream>

#include "csg.hpp"
#include "csg_types.hpp"
#include "csg_exception.hpp"

namespace {

const std::string csg_ext = ".csg";

bool has_ext(std::string filename, std::string extension) {
  return filename.rfind(extension) == (filename.size() - extension.size());
}

std::string read_csg(std::string geom_file) {
  if (!has_ext(geom_file, csg_ext)) {
     std::string err_msg = "ERROR:  Filename  " + geom_file + "  must have extension .csg";
     std::cout << err_msg << std::endl;
     throw csg::Exception("read_csg", err_msg);
  }

  std::ifstream in_stream;
  in_stream.open(geom_file.c_str());

  if (!in_stream.good()) {
     std::string err_msg = "ERROR:  Cannot open file   " + geom_file;
     std::cout << err_msg << std::endl;
     throw csg::Exception("read_csg", err_msg);
  }

  std::stringstream str_stream;
  str_stream << in_stream.rdbuf();
  return str_stream.str();
}

} // namespace

namespace csg {

class CsgIF::Impl {
public:
  double call_signed_distance(const std::shared_ptr<Tree> &a_tree, double xx,
                              double yy, double zz) const {
    return signed_distance_3d(a_tree->top, xx, yy, zz);
  }
};

CsgIF::CsgIF(const CsgIF &rhs)
    : m_state(rhs.m_state), m_is_internal_flow(rhs.m_is_internal_flow),
      m_pimpl(std::make_unique<Impl>(*rhs.m_pimpl)) {}

CsgIF::CsgIF(std::shared_ptr<Tree> a_state, bool is_internal_flow)
    : m_state(a_state), m_is_internal_flow(is_internal_flow),
      m_pimpl(std::make_unique<Impl>()) {}

CsgIF::CsgIF(CsgIF &&rhs) noexcept = default;
CsgIF::~CsgIF() = default;

double CsgIF::operator()(double xx, double yy, double zz) const noexcept {
  auto sd = m_pimpl->call_signed_distance(m_state, xx, yy, zz);
  if (m_is_internal_flow) {
    return -sd;
  }
  return sd;
}

std::shared_ptr<csg::Tree> get_csgtree_from_filename(std::string geom_file) {
  auto csg_str = read_csg(geom_file);
  return csg::parse_csg(csg_str);
}

std::unique_ptr<CsgIF> get_csgif_from_filename(std::string geom_file) {
  auto csg_obj = get_csgtree_from_filename(geom_file);

  csg::CsgIF csg_if(csg_obj);
  return std::make_unique<CsgIF>(csg_if);
}

std::unique_ptr<CsgIF> get_csgif(const std::string &geom_file,
                                 bool is_internal_flow) {
  auto csg_obj = get_csgtree_from_filename(geom_file);

  csg::CsgIF csg_if(csg_obj, is_internal_flow);
  return std::make_unique<CsgIF>(csg_if);
}

} // namespace csg
