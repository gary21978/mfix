#ifndef MFIX_IO_REPORTS_DIAMETER_H_
#define MFIX_IO_REPORTS_DIAMETER_H_

#include <AMReX_ParmParse.H>

#include <mfix_regions.H>
#include <mfix_timer.H>
#include <mfix_pc.H>

namespace reports {

class diameter
{

  public:

    diameter ( int const a_is_pic, MFIXRegions& a_regions );

    void fill_bin_data ( int const a_rt, amrex::Vector<amrex::Real>& a_data,
                         MFIXParticleContainer* a_p_pc );

    void write_now ( MFIXTimer& a_timer, amrex::Real a_dt,
                     bool a_first, bool a_last,
                     MFIXParticleContainer* a_p_pc);

  private:

    int const m_is_pic;

    int m_count = 0;

    amrex::Vector<std::string> m_reports;

    amrex::Vector< const amrex::RealBox*> m_region;

    amrex::Vector<amrex::Real> m_min;
    amrex::Vector<amrex::Real> m_max;

    amrex::Vector<int> m_int;
    amrex::Vector<int> m_last_int;

    amrex::Vector<amrex::Real> m_per_approx;

    amrex::Vector<int> m_bins;

    void write_bin_data (int const a_rt, int const a_step,
                         amrex::Vector<amrex::Real>& a_data);

    std::tuple<amrex::Real, amrex::Real>
    get_bin_extents ( int const a_rt, int const a_bin, amrex::Real a_scale = 1.0) const
    { AMREX_ALWAYS_ASSERT(a_bin >= 0 && a_bin < m_bins[a_rt]);

      amrex::Real blo = std::numeric_limits<amrex::Real>::min();
      amrex::Real bhi = std::numeric_limits<amrex::Real>::max();

      amrex::Real const bin_r  = static_cast<amrex::Real>(a_bin);
      amrex::Real const bins_r = static_cast<amrex::Real>(m_bins[a_rt]);

      amrex::Real const bw = (m_max[a_rt] - m_min[a_rt]) / (bins_r - 2.0);

      if (a_bin > 0) { blo = m_min[a_rt] + bw*(bin_r - 1.0); }
      if (a_bin < m_bins[a_rt]-1) { bhi = m_min[a_rt] + bw*bin_r; }

      return {a_scale*blo, a_scale*bhi};
    }

    amrex::Real get_midpoint ( int const a_rt, int const a_bin) const
    { AMREX_ALWAYS_ASSERT(a_bin > 0 && a_bin < m_bins[a_rt]-1);

      amrex::Real const bin_r  = static_cast<amrex::Real>(a_bin);
      amrex::Real const bins_r = static_cast<amrex::Real>(m_bins[a_rt]);

      amrex::Real const bw = (m_max[a_rt] - m_min[a_rt]) / (bins_r - 2.0);

      return m_min[a_rt] + bw*(bin_r - 0.5);
    }

};

} // namespace reports

#endif
