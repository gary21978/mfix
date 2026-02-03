#ifndef _MFIX_FIX_INPUTS_H_
#define _MFIX_FIX_INPUTS_H_

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include <mfix_reporter.H>

class FixInputs {

  public:

    FixInputs ( )
      : m_verbose(1)
      , m_has_date(0)
    {}

    FixInputs ( int const a_verbose )
      : m_verbose(a_verbose)
      , m_has_date(0)
    {}

    FixInputs ( std::string const a_date )
      : m_verbose(1)
      , m_has_date(1)
      , m_date(a_date)
    {}


    template <class T>
    void swap ( const std::string& old_input,
                const std::string& new_input,
                int const a_use_old = 0 )
    {

      if ( check_old(old_input, new_input) ) {

        T input_value;
        m_pp.get(old_input, input_value);
        m_pp.remove(old_input);

        if (m_pp.contains(new_input) ) {

          if (a_use_old) {

            warn_overwrite( old_input, new_input);

            m_pp.remove(new_input);
            m_pp.add(new_input, input_value);

          } else if (m_verbose) {

             report_clash(old_input, new_input);

          }
        } else {
          m_pp.add(new_input, input_value);
        }
      }
    }

    template <class T>
    void swap_arr ( const std::string& old_input,
                    const std::string& new_input,
                    int const a_use_old = 0 )
    {

      if ( check_old( old_input, new_input) ) {

        T input_value;
        m_pp.getarr(old_input, input_value);
        m_pp.remove(old_input);

        if (m_pp.contains(new_input) ) {

          if (a_use_old) {

            warn_overwrite(old_input, new_input);

            m_pp.remove(new_input);
            m_pp.addarr(new_input, input_value);

          } else if (m_verbose) {

            report_clash(old_input, new_input);

          }
        } else {
          m_pp.addarr(new_input, input_value);
        }
      }

    }


    void remove ( const std::string& input_string )
    {
      if (m_pp.contains(input_string)) {
        if (m_verbose) {
          reporter::Log(reporter::Warning)
            << "Input " << std::string(input_string) << " was removed.\n"
            << add_date()
            << "This warning message will be removed in the future.\n"
            << "Please update the inputs file.";
        }
        m_pp.remove(input_string);
      }
    }


    int check_old ( const std::string& old_input,
                    const std::string& new_input)
    {

      int const contains_old( m_pp.contains(old_input) );

      if (contains_old && m_verbose) {

        reporter::Log(reporter::Info)
          << "Keyword " << std::string(old_input) << " was replaced by "
          << std::string(new_input) << ".\n"
          << add_date()
          << "Update the inputs file, because automatic conversion for this\n"
          << "keyword may be removed in the future.";
      }

      return contains_old;
    }


    void warn_overwrite ( const std::string& old_input,
                          const std::string& new_input)
    {
      reporter::Log(reporter::Warning)
       << "Overwriting " << old_input << " "
       << "with value from " << new_input;
    }

    void report_clash ( const std::string& old_input,
                        const std::string& new_input)
    {
      reporter::Log(reporter::Info)
        << "This is embarrassing; old and new keywords found in the inputs file.\n"
        << "Old: " << old_input << "  "
        << "New: " << new_input << "\n"
        << "Update the inputs file to remove the old keyword.";
    }

  private:

    int const m_verbose;
    int const m_has_date;
    std::string m_date;

    std::string add_date () const {
      if (!m_has_date) { return ""; }
      return "Date: " + m_date + "\n";
    }

    amrex::ParmParse m_pp;
};

#endif
