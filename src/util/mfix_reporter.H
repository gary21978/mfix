#ifndef MFIX_REPORT_H_
#define MFIX_REPORT_H_

#include <map>
#include <AMReX_Print.H>

namespace reporter {

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


struct threes : std::numpunct<char>
{
    std::string do_grouping() const { return "\3"; }
};

template<class T>
static std::string FormatWithCommas(T value)
{
    std::stringstream ss;
    ss.imbue(std::locale(std::cout.getloc(), new threes));
    ss << std::fixed << value;
    return ss.str();
}


static constexpr int Error   = 1;
static constexpr int Warning = 2;
static constexpr int Info    = 3;
static constexpr int Status  = 4;

class Log
{

  public:

    explicit Log (int const a_level, std::string const a_file, int const a_line)
      : m_nl("\n")
      , m_level(a_level)
      , m_file_info(1)
      , m_file(a_file)
      , m_line(std::to_string(a_line))
    { }


    Log (int const a_level)
      : m_nl("\n")
      , m_level(a_level)
      , m_file_info(0)
    { }

    ~Log () {
#ifndef MFIX_SILENT_REPORTER
      std::string message =
          Border() +
          Preamble() + LogLevel() + ss.str() +
          Border() + m_nl;

      amrex::Print() << message;
#endif
#ifndef MFIX_NO_ABORT
      if (m_level == Error) { amrex::Abort(); }
#endif
    }

    template <typename T>
    Log& operator<< (const T& x) {
        ss << x;
        return *this;
    }

  private:

    std::ostringstream ss;

    std::string m_nl;

    int const m_level;
    int const m_file_info;

    std::string m_file;
    std::string m_line;


    static constexpr int m_border_width = 80;

    std::string LogLevel () const {
      if (m_level == Error) { return "Error: "; }
      else if (m_level == Warning) { return "Warning: "; }
      else if (m_level == Info   ) { return "Info: ";    }
      return "";
    }

    std::string Preamble () const {
      if (m_file_info && (m_level == Error || m_level == Warning)) {
        return "File: " + m_file + m_nl +
               "Line: " + m_line + m_nl + m_nl;
      }
      return "";
    }

    std::string Border () const {
      if (m_level == Status) { return ""; }
      char b = (m_level == Error) ? '!' : ((m_level == Warning) ? '=' : '-');
      return m_nl + std::string(m_border_width, b) + m_nl;
    }

};

} // namespace reporter

#endif
