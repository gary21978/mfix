#ifndef MFIX_DISTRIBUTIONS_CUSTOM_H_
#define MFIX_DISTRIBUTIONS_CUSTOM_H_

#include <AMReX_VisMF.H>

namespace distributions {

struct custom {

  public:

    custom ( std::string a_filename )
    : m_has_error(0)
    , m_is_pdf(0)
    , m_bins(-1)
    , m_can_interp(0)
    , m_min(0.)
    , m_max(0.)
    , m_Nw_mean(0.)
    , m_Nw_variance(0.)
    , m_Vw_mean(0.)
    , m_Vw_variance(0.)
    { read_and_bcast(a_filename); }

    int ok () { return !m_has_error; }

    //std::string error_msg () { return m_error; }
    int bins () { return m_bins; }

    amrex::Real min () { return m_min; }
    amrex::Real max () { return m_max; }

    amrex::Real Nw_mean () { return m_Nw_mean; }
    amrex::Real Nw_variance () { return m_Nw_variance; }
    amrex::Real Nw_stddev () { return std::sqrt(m_Nw_variance); }

    amrex::Real Vw_mean () { return m_Vw_mean; }
    amrex::Real Vw_variance () { return m_Vw_variance; }
    amrex::Real Vw_stddev () { return std::sqrt(m_Vw_variance); }

    int can_interp () { return m_can_interp; }

  private:

    int m_has_error;

    int m_is_pdf;
    int m_bins;

    int m_can_interp;

    amrex::Real m_min;
    amrex::Real m_max;

    amrex::Real m_Nw_mean;
    amrex::Real m_Nw_variance;

    amrex::Real m_Vw_mean;
    amrex::Real m_Vw_variance;

    void read_and_bcast ( std::string  const a_file )
    {

      if ( amrex::ParallelDescriptor::IOProcessor() ) {
        std::ifstream iss;
        iss.open(a_file, std::ios::in);
        m_has_error = !iss.good();
        iss.close();
      }
      amrex::ParallelDescriptor::Bcast(&m_has_error, 1,
        amrex::ParallelDescriptor::IOProcessorNumber());

      // Unable to open distribution file
      if (m_has_error) { return; }

      amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());

      amrex::Vector<char> fileCharPtr;
      amrex::ParallelDescriptor::ReadAndBcastFile(a_file, fileCharPtr);
      std::string fileCharPtrString(fileCharPtr.dataPtr());
      std::istringstream is(fileCharPtrString, std::istringstream::in);

      { std::string s_val;

        std::string line;
        std::getline(is, line);
        std::replace(line.begin(), line.end(), ',', ' ');

        std::istringstream lis(line);

        lis >> s_val;
        if (lis.bad()) {
          //m_error = "Unable to determine number of bins";
          m_has_error = 1;
          return;
        }

        m_bins = std::stoi(s_val);
      }

      int is_cdf(0);
      { std::string dist_type;

        std::string line;
        std::getline(is, dist_type); // PDF, CDF

        m_is_pdf = (amrex::toLower(dist_type).compare("pdf") == 0) ? 1 : 0;
        is_cdf = (amrex::toLower(dist_type).compare("cdf") == 0) ? 1 : 0;

        if (!(m_is_pdf || is_cdf)) {
          //m_error = "Unable to determine if distribution is CDF or PDF";
          m_has_error = 1;
          return;
        }
      }

      if (m_bins <= 1) {
        //m_error = "Invalid number of bins specifying distribution";
        m_has_error = 1;
        return;
      }
    }

  public:

    void setup_distribution ( std::string a_filename, int const a_bins,
             int const a_is_number_weighted, amrex::Real* a_diam,
             amrex::Real* a_Nw_CDF, amrex::Real* a_Vw_CDF )
    {

      amrex::Vector<char> fileCharPtr;
      amrex::ParallelDescriptor::ReadAndBcastFile(a_filename, fileCharPtr);
      std::string fileCharPtrString(fileCharPtr.dataPtr());
      std::istringstream is(fileCharPtrString, std::istringstream::in);

      { std::string line;
        std::getline(is, line); // skip 1: Number of bins
        std::getline(is, line); // skip 2: CDF or PDF
        std::getline(is, line); // skip 3: comment
      }

      amrex::Real* p_CDF = (a_is_number_weighted ? a_Nw_CDF : a_Vw_CDF);

      for (int n(0); n<m_bins; ++n) {

        std::string line;

        std::getline(is, line);
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream lis(line);

        std::string s_val, s_prob;
        amrex::Real bval(std::numeric_limits<amrex::Real>::max());
        amrex::Real bprob(std::numeric_limits<amrex::Real>::max());

        lis >> s_val;
        if (lis.bad()) { m_has_error = 1; }
        { try { bval = std::stod(s_val); }
          catch (...) {m_has_error = 1; }
        }

        lis >> s_prob;
        if (lis.bad()) { m_has_error = 1; }
        { try { bprob = std::stod(s_prob); }
          catch (...) {m_has_error = 1; }
        }

        if (m_has_error) { return; }

        a_diam[n] = bval;
        p_CDF[n] = bprob;

        if(m_is_pdf && n>0) { p_CDF[n] += p_CDF[n-1]; }

        // verify monotonic inputs
        if (n-1 > 0) {
          if ( a_diam[n] < a_diam[n-1] ) {
            m_has_error = 1;
            return;
          }
          if ( p_CDF[n] < p_CDF[n-1] ) {
            m_has_error = 1;
            return;
          }
        }
      }
      if ( amrex::almostEqual(p_CDF[0], 0.) ) { m_can_interp = 1; }

      // Normal the CDF to run from [0,1];
      {
        amrex::Real const scale = 1.0/p_CDF[(m_bins-1)];
        for (int b(0); b<(m_bins); ++b) { p_CDF[b] *= scale; }
      }

      m_min = a_diam[0];
      m_max = a_diam[(m_bins-1)];

      a_diam[0] = m_min;
      a_diam[a_bins-1] = m_max;

      if (m_has_error) { return; }

      if (a_is_number_weighted) {
        setup_Nw_distribution(a_bins, a_diam, a_Nw_CDF, a_Vw_CDF);
      } else {
        setup_Vw_distribution(a_bins, a_diam, a_Nw_CDF, a_Vw_CDF);
      }
    }

  private:

    void setup_Nw_distribution ( int const a_bins,
                                 amrex::Real* const a_diam,
                                 amrex::Real* const a_Nw_CDF,
                                 amrex::Real* const a_Vw_CDF)
    {
      //amrex::Print() << "\n\nmfix_distributions_custom -- setup_NW_distribution --> Nw\n";

      a_Vw_CDF[0] = 0.0;

      for (int bin(m_can_interp); bin<a_bins; ++bin) {

        amrex::Real d(a_diam[bin]);
        amrex::Real p(a_Nw_CDF[bin]);
        if (bin > 0) { p -= a_Nw_CDF[bin-1]; }

        m_Nw_mean += d*p;
        m_Nw_variance += d*d*p;

        a_Vw_CDF[bin] = p*d*d*d;
        if (bin > 0) { a_Vw_CDF[bin] += a_Vw_CDF[bin-1]; }

        //amrex::Print() << bin << "  " << d << "  " << p << "  " << Vw_PDF << "\n";
      }
      m_Nw_variance -= m_Nw_mean*m_Nw_mean;

      //amrex::Print() << " Mean: " << m_Nw_mean << "\n"
      //               << " Var:  " << m_Nw_variance << "\n";

      //amrex::Print() << "\nfix_distributions_custom -- setup_NW_distribution --> Vw\n";

      amrex::Real const scale(1./a_Vw_CDF[a_bins-1]);
      for (int bin(m_can_interp); bin<a_bins; ++bin) {
        a_Vw_CDF[bin] *= scale;

        amrex::Real const d(a_diam[bin]);
        amrex::Real p(a_Vw_CDF[bin]);
        if (bin > 0) { p -= a_Vw_CDF[bin-1]; }

        //amrex::Print() << bin << "  " << d << "  " << p << "\n";

        m_Vw_mean += d*p;
        m_Vw_variance += d*d*p;
      }
      m_Vw_variance -= m_Vw_mean*m_Vw_mean;

      //amrex::Print() << "Vw Mean: " << m_Vw_mean << "\n"
      //               << "Vw Var:  " << m_Vw_variance << "\n";
    }

    void setup_Vw_distribution ( int const a_bins,
                                 amrex::Real* const a_diam,
                                 amrex::Real* const a_Nw_CDF,
                                 amrex::Real* const a_Vw_CDF)
    {
      //amrex::Print() << "\n\nmfix_distributions_custom -- setup_Vw_distribution --> Vw\n";

      a_Nw_CDF[0] = 0.0;

      for (int bin(m_can_interp); bin<a_bins; ++bin) {

        amrex::Real const d(a_diam[bin]);
        amrex::Real p(a_Vw_CDF[bin]);
        if (bin > 0) { p -= a_Vw_CDF[bin-1]; }

        m_Vw_mean += d*p;
        m_Vw_variance += d*d*p;

        a_Nw_CDF[bin] = p/(d*d*d);
        if (bin > 0) { a_Nw_CDF[bin] += a_Nw_CDF[bin-1]; }

        //amrex::Print() << bin << "  " << d << "  " << a_Vw_CDF[bin] << "  " << p << "\n";
      }
      m_Vw_variance -= m_Vw_mean*m_Vw_mean;

      //amrex::Print() << "Vw Mean: " << m_Vw_mean << "\n"
      //               << "Vw Var:  " << m_Vw_variance << "\n";

      //amrex::Print() << "\n\nmfix_distributions_custom -- setup_Vw_distribution --> Nw\n";

      amrex::Real const scale(1./a_Nw_CDF[a_bins-1]);
      for (int bin(m_can_interp); bin<a_bins; ++bin) {
        a_Nw_CDF[bin] *= scale;

        amrex::Real const d(a_diam[bin]);
        amrex::Real p(a_Nw_CDF[bin]);
        if (bin > 0) { p -= a_Nw_CDF[bin-1]; }

        //amrex::Print() << bin << "  " << d << "  " << a_Nw_CDF[bin] << "  " << p << "\n";

        m_Nw_mean += d*p;
        m_Nw_variance += d*d*p;
      }
      m_Nw_variance -= m_Nw_mean*m_Nw_mean;

      //amrex::Print() << "Nw Mean: " << m_Nw_mean << "\n"
      //               << "Nw Var:  " << m_Nw_variance << "\n";
    }

};

} // end namespace distribution

#endif
