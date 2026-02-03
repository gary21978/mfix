#ifndef MFIX_VERSIONS_H_
#define MFIX_VERSIONS_H_

#include <iostream>
#include <string>
#include <sstream>


namespace version {


struct checkpoint {

  checkpoint () = default;

  checkpoint ( std::string& a_version ) {
    std::vector<int> version = split(a_version);
    std::get<0>(m_version) = version[0];
    std::get<1>(m_version) = version[1];
  }

  // Current version of checkpoint files: major.minor.patch
  static constexpr int major = 1;
  static constexpr int minor = 6;

  static constexpr std::tuple<int,int> current {major,minor};

  // Return the current checkpoint version as a string
  std::string current_str () {

    std::stringstream ss;
      ss << static_cast<int>(major) << "."
         << static_cast<int>(minor);

    return ss.str();
  }

  std::vector<int> split ( const std::string& a_version_str) {

    std::vector<int> version_ints;

    std::stringstream ss(a_version_str);
    std::string item;

    while (std::getline(ss, item, '.')) {
      version_ints.push_back(std::stoi(item));
    }

    // Ensure the vector has three elements
    while (version_ints.size() < 2) {
        version_ints.push_back(0); // Default to 0 if missing
    }
    AMREX_ALWAYS_ASSERT(version_ints.size() == 2);
    return version_ints;
  }

  std::tuple<int,int> const& operator() () const noexcept
  { return m_version; }

  bool is_current () const noexcept
  { return (m_version == current); }

  std::tuple<int,int> m_version;
};

} // end namespace version


#endif
