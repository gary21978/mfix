#ifndef POST_MFIX_REPORTER_H_
#define POST_MFIX_REPORTER_H_

#include <AMReX_Print.H>

namespace post_mfix {

AMREX_FORCE_INLINE
std::string format(int const a_width,
                   int const a_percision,
                   amrex::Real const a_number)
{
  std::stringstream ss;

  ss << std::setw(a_width) << std::setprecision(a_percision);

  amrex::Real minval = std::pow(1.0, 1-a_percision);
  amrex::Real maxval = std::pow(1.0, a_width - (4+a_percision));
  if (a_number < std::numeric_limits<amrex::Real>::min() ||
     (a_number >= minval && a_number < maxval) )
  { ss << std::fixed << a_number; }
  else { ss << std::scientific << a_number; }

  return ss.str();
}


// Error
// Warning
// Status
// Info

class Error
{

  public:

    explicit Error (std::string const a_file, int const a_line)
      : m_boarder("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
      , m_nl("\n")
      , m_file(a_file)
      , m_line(std::to_string(a_line))
    { }

    ~Error () {

      std::string message =
          m_nl + m_boarder +
          m_nl + " File: " + m_file +
          m_nl + " Line: " + m_line + m_nl +
          m_nl + ss.str() +
          m_nl + m_boarder + m_nl +
          m_nl;

      amrex::Print() << message;
      amrex::Abort();

    }

    template <typename T>
    Error& operator<< (const T& x) {
        ss << x;
        return *this;
    }




  private:

    std::ostringstream ss;

    std::string m_boarder;
    std::string m_nl;

    std::string m_file;
    std::string m_line;
};

} // namespace reporter

#endif
