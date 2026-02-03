#ifndef CSG_EXCEPTION_H
#define CSG_EXCEPTION_H

#include <stdexcept>
#include <string>

namespace csg {

class Exception : public std::exception {
public:
  Exception();
  explicit Exception(const std::string &message);
  Exception(const std::string &source, const std::string &message);
  const char *what() const noexcept;

private:
  std::string make_message(const std::string &source,
                           const std::string &message);
  std::string m_msg;
};

} // namespace csg

#endif
