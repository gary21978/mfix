#ifndef _MFIX_FVAR_H_
#define _MFIX_FVAR_H_

#include <AMReX_ParmParse.H>

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>

struct FVAR_ {

  public:

    FVAR_ ( amrex::Vector<int> a_test_sizes )
      : m_entries(0)
      , m_defined(0)
      , m_constant(1)
      , m_size(a_test_sizes.size())
      , m_test_sizes(a_test_sizes)
    { }

    // simple methods to access information
    int is_defined () const { return m_defined; }
    int is_constant () const { return m_constant; }

    int size () const {
      if (m_values.size() == 0) { return 0; }
      int const s = m_values[0].size();
      return (m_constant ? s : s-1);
    }

    // consume variable data from ParmParse
    int get( amrex::ParmParse a_pp, std::string a_name)
    {
      if (m_defined) { return 1;}

      m_entries = a_pp.countname(a_name);

      m_defined = ((m_entries == 0) ? 0 : 1);
      if ( !m_defined ) { return 0; }

      // Multi valued inputs must be the same length
      if ( m_entries > 1 ) {
        int min_size(std::numeric_limits<int>::max());
        int max_size(std::numeric_limits<int>::lowest());
        for (int k(0); k<m_entries; ++k) {
          amrex::Vector<amrex::Real> input;
          if (a_pp.queryktharr(a_name, k, input)) {
            int const sz = input.size();
            min_size = amrex::min(min_size, sz);
            max_size = amrex::max(max_size, sz);
          }
        }
        if (min_size != max_size) { return 1; }
      }


      { amrex::Vector<amrex::Real> input;
        a_pp.queryarr(a_name, input);

        if (!is_invalid(input.size())) {
          m_entries = 1;
        } else {
          m_constant = 0;
        }
      }

      m_values.resize(m_entries);

      if (m_constant) {

          a_pp.queryarr(a_name, m_values[0]);

          if ( is_invalid(m_values[0].size()) ) { return 1; }

      } else {

        for (int k(0); k<m_entries; ++k) {

          if (a_pp.queryktharr(a_name, k, m_values[k])) {
            if ( is_invalid(m_values[k].size(), 1) ) { return 1; }
          }
          else { return 1; }

        } // k-th entry

        // Duplicate time entries are invalid
        for (int i(0); i<m_entries; ++i) {
          for (int j(i+1); j<m_entries; ++j) {
            if (amrex::almostEqual(m_values[i][0], m_values[j][0]))
            { return 1; }
          }
        }

        // Sort based on time.
        for (int i(0); i<m_entries; ++i) {
          for (int j(0); j < m_entries-i-1; ++j) {
            if ( m_values[j][0] > m_values[j+1][0] )
            { m_values[j].swap(m_values[j+1]); }
          } // j-loop
        } // i-loop
      }

      return 0;
    }

    // Getting the variable out of the struct
    amrex::Vector<amrex::Real> operator() (amrex::Real a_time = 0.) const
    {
      if (m_constant) { return m_values[0]; }

      else {

        int const vsize = m_values[0].size()-1;

        amrex::Vector<amrex::Real> result(vsize,0.);

        // There is only one entry (use it) OR the current time
        // is less than the first entry so return the first entry.
        // amrex::Print() << "Defaulting to the first entry " << tg_table[0][1] << "\n";
        if( m_entries == 1 || a_time <= m_values[0][0] ) {

          std::copy(m_values[0].begin()+1, m_values[0].end(), result.begin());
          return result;

        } else {

          //Interpolating from user provided values
          for (int entry(1); entry<m_entries; ++entry) {

            // Time of current entry
            const amrex::Real entry_time = m_values[entry][0];

            // Too close for any calculations
            if ( amrex::almostEqual(a_time, entry_time) ) {

              std::copy(m_values[entry].begin()+1, m_values[entry].end(), result.begin());
              return result;

            // The entries are sorted based on time, so if the current time is
            // less than this entries time, then we need to interpolate between
            // the previous time and this time.
            } else if (a_time < entry_time) {

              const amrex::Real prv_time = m_values[entry-1][0];
              if (amrex::almostEqual(a_time, prv_time)) {

                std::copy(m_values[entry-1].begin()+1, m_values[entry-1].end(), result.begin());
                return result;

              } else {
                const amrex::Real delta_time = entry_time - prv_time;
                const amrex::Real whi = (a_time - prv_time) / delta_time;
                const amrex::Real wlo = 1.0 - whi;

                for (int i(0); i<vsize; ++i) {
                  result[i] = whi*m_values[entry  ][i+1] +
                              wlo*m_values[entry-1][i+1];
                }
                return result;
              }
            }
          }
          // If we made it here, time is larger than all the entries in the
          // table, therefore, use the last value.
          std::copy(m_values[m_entries-1].begin()+1, m_values[m_entries-1].end(), result.begin());
          return result;
        }
      }
    }

  private:

    int m_entries;
    int m_defined;
    int m_constant;

    amrex::Vector<amrex::Vector<amrex::Real>> m_values;

    int m_size;
    amrex::Vector<int> m_test_sizes;

    // Test if the provided input has a valid length
    int is_invalid ( int const a_size, int const a_offset=0 ) const {
      for ( int s(0); s<m_test_sizes.size(); ++s)
      { if ( a_size == m_test_sizes[s] + a_offset ) { return 0; } }
      return 1;
    }
};

#endif
