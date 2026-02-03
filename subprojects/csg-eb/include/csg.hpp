#ifndef CSG_H_
#define CSG_H_

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace csg {

struct Tree;
std::shared_ptr<Tree> parse_csg(std::string);

class CsgIF {
public:
  CsgIF(std::shared_ptr<Tree> a_state, bool is_internal_flow = false);
  CsgIF(const CsgIF &rhs);
  CsgIF(CsgIF &&rhs) noexcept;
  CsgIF &operator=(const CsgIF &rhs) = delete;
  CsgIF &operator=(CsgIF &&rhs) = delete;
  ~CsgIF();

  double operator()(double xx, double yy, double zz) const noexcept;

  inline double operator()(const std::array<double, 3> &p) const noexcept {
    return this->operator()(p[0], p[1], p[2]);
  }

private:
  std::shared_ptr<Tree> m_state;
  bool m_is_internal_flow;

  class Impl;
  std::unique_ptr<Impl> m_pimpl;
};

std::shared_ptr<csg::Tree> get_csgtree_from_filename(std::string);
std::unique_ptr<CsgIF>
get_csgif_from_filename(std::string geom_file); // TODO: Deprecate once mfix has
                                                // stopped using this call
std::unique_ptr<CsgIF> get_csgif(const std::string &geom_file,
                                 bool is_internal_flow = false);

} // namespace csg

#endif
