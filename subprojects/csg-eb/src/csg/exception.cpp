#include "csg_exception.hpp"

#include <sstream>
#include <string>

namespace csg {

Exception::Exception() : m_msg(make_message("No Source", "No message")) {}

Exception::Exception(const std::string &message)
    : m_msg(make_message("No Source", message)) {}

Exception::Exception(const std::string &source, const std::string &message)
    : m_msg(make_message(source, message)) {}

const char *Exception::what() const noexcept { return m_msg.c_str(); }

std::string Exception::make_message(const std::string &source,
                                    const std::string &message) {
  std::stringstream s;
  s << "Exception Data:" << std::endl;
  s << "Source  : " << source << std::endl;
  s << "Message : " << message << std::endl;
  return s.str();
}

} // namespace csg
