#include <mfix_reporter.H>

#include <mfix_io_reports_diameter.H>

using namespace amrex;

reports::diameter::
diameter ( int const a_is_pic, MFIXRegions& a_regions)
  : m_is_pic(a_is_pic)
{
  ParmParse ppDiamReports("mfix.reports");

  if (!ppDiamReports.queryarr("diameter", m_reports)) { return; }

  m_count = m_reports.size();
  m_region.resize( m_count, nullptr);

  m_min.resize( m_count, std::numeric_limits<amrex::Real>::max());
  m_max.resize( m_count, std::numeric_limits<amrex::Real>::min());

  m_int.resize( m_count, std::numeric_limits<int>::max());
  m_last_int.resize( m_count, 0);

  m_per_approx.resize( m_count, std::numeric_limits<amrex::Real>::max());

  m_bins.resize( m_count, 0);

  for (int rt(0); rt<m_count; ++rt) {

    std::string report_diameter_report_prefix = "mfix.reports.diameter." + m_reports[rt];
    ParmParse ppReport(report_diameter_report_prefix);

    std::string region;
    if ( !ppReport.query("region", region) ) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input " << ppReport.ParserPrefix << ".region not found.";
    }

    m_region[rt] = a_regions.getRegion(region);
    if ( m_region[rt] == nullptr ) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Invalid diameter report region";
    }

    if ( !ppReport.query("bins", m_bins[rt]) ) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input " << ppReport.ParserPrefix << ".bins not found";
    }

    if ( !(m_bins[rt] > 0) ) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input " << ppReport.ParserPrefix << ".bins must be greater than zero";
    }

    // add bins to catch under and over flow.
    m_bins[rt] += 2;

    if ( !ppReport.query("min", m_min[rt]) ) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input " << ppReport.ParserPrefix << ".min not found";
    }

    if ( !ppReport.query("max", m_max[rt]) ) {
      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input " << ppReport.ParserPrefix << ".max not found";
    }

    int const has_A = ppReport.query("int", m_int[rt]);
    int const has_B = ppReport.query("per_approx", m_per_approx[rt]);

    if ( (has_A && has_B) || !(has_A || has_B) ) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Required input for " << ppReport.ParserPrefix << " not found!\n"
        << "Either " << ppReport.ParserPrefix << ".int or " << ppReport.ParserPrefix << ".per_approx\n"
        << "must be specified.";
    }

  }

}

void
reports::diameter::
write_now ( MFIXTimer& a_timer, Real a_dt,
            bool a_first, bool a_last,
            MFIXParticleContainer* a_p_pc)
{
  for (int rt(0); rt<m_count; ++rt) {

    int test = 0;

    if ( a_first ) {

      test = 1;

    } else if ( a_last && (a_timer.nstep() != m_last_int[rt]) ) {

      test = 1;

    } else if ( m_int[rt] > 0 ) {

      test = (a_timer.nstep() % m_int[rt] == 0) ? 1 : 0;

    } else if ( !amrex::almostEqual(m_per_approx[rt], 0.) ) {

      test = a_timer.test_per_approx(a_dt, m_per_approx[rt]);

    }

    if ( test == 1) {

      Vector<Real> bin_data(4*m_bins[rt], 0.);

      for (int bin(0); bin<4*m_bins[rt]; ++bin)
      { bin_data[bin] = 0.; }

      fill_bin_data(rt, bin_data, a_p_pc);

      write_bin_data(rt, a_timer.nstep(), bin_data);

      m_last_int[rt] = a_timer.nstep();

      bin_data.clear();
    }

  } // loop over report types

}



void
reports::diameter::
fill_bin_data (int const a_rt, Vector<Real>& a_data,
               MFIXParticleContainer* a_p_pc)
{
  int const bins = m_bins[a_rt];

  using MFIXParIter = MFIXParticleContainer::MFIXParIter;
  using PairIndex = MFIXParticleContainer::PairIndex;

  const RealBox region(m_region[a_rt]->lo(), m_region[a_rt]->hi());

  for (int bin(0); bin<bins; ++bin) {

    //auto const [bin_l, bin_h] = get_bin_extents(a_rt, bin, 0.5);
    auto const bin_exts = get_bin_extents(a_rt, bin, 0.5);

    Real const bin_lo = std::get<0>(bin_exts);
    Real const bin_hi = std::get<1>(bin_exts);

    // Reduce sum over new particles to take into account
    // that the diameters may be different.
    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Real, Real, Real, Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (int lev(0); lev<a_p_pc->numLevels(); ++lev) {

      for (MFIXParIter pti(*a_p_pc, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto& ptile = a_p_pc->GetParticles(lev)[index];

        const int np = ptile.numParticles();

        auto const& aos = ptile.GetArrayOfStructs();
        auto const* pstruct = aos().dataPtr();

        auto const& p_real = (ptile.GetStructOfArrays()).realarray();

        auto const& p_radius = p_real[SoArealData::radius];
        auto const& p_statwt = p_real[SoArealData::statwt];
        auto const& p_density = p_real[SoArealData::density];

        reduce_op.eval(np, reduce_data, [pstruct, region, bin_lo, bin_hi,
          is_pic=m_is_pic, p_radius, p_statwt, p_density]
        AMREX_GPU_DEVICE (int pid) -> ReduceTuple
        {
          if ( region.contains(pstruct[pid].pos()) &&
              (bin_lo <= p_radius[pid] && p_radius[pid] < bin_hi)) {

            Real const statwt = is_pic ? p_statwt[pid] : 1.0;
            Real const p_volume = SoArealData::volume(p_radius[pid]);
            Real const p_mass = p_density[pid]*p_volume;

            return { 1., statwt, statwt*p_volume, statwt*p_mass};

          }
          return {0., 0., 0., 0.};
        });

      } // MFIXParIter
    } // lev

    ReduceTuple tuple = reduce_data.value();

    a_data[bin + 0*bins] += amrex::get<0>(tuple); // parcel count
    a_data[bin + 1*bins] += amrex::get<1>(tuple); // particle count
    a_data[bin + 2*bins] += amrex::get<2>(tuple); // total volume
    a_data[bin + 3*bins] += amrex::get<3>(tuple); // total mass

  } // loop over bins
}


void
reports::diameter::
write_bin_data (int const a_rt, int const a_step,
                Vector<Real>& a_data)
{

  ParallelDescriptor::ReduceRealSum(a_data.dataPtr(), a_data.size(),
                                    ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor()) {

    Vector<Real> sum_data(4, 0.);
    Vector<Real> inv_sum_data(4, 0.);

    int const bins = m_bins[a_rt];

    // parcels, particles, mass, volume
    for (int n(0); n<4; ++n) {
      for (int bin(0); bin<bins; ++bin) {
        sum_data[n] += a_data[bin + n*bins];
      }
      if (sum_data[n] > 0) { inv_sum_data[n] = 1.0/sum_data[n]; }
    }

    std::string filename = m_reports[a_rt] + amrex::Concatenate(".diam",a_step) + ".csv";

    std::ofstream File;
    File.open(filename, std::ios::out|std::ios::trunc);

    File << "# bin minimum diameter: " << m_min[a_rt] << '\n';
    File << "# bin maximum diameter: " << m_max[a_rt] << '\n';
    File << "# Number of bins: " << m_bins[a_rt]-2 << '\n';
    File << "# \n";
    File << "# Row descriptions:\n";
    { int row(1);
      File << "#  " << std::setw(2) << row << ". bin mean diameter\n"; ++row;
      if (m_is_pic) { File << "#  " << std::setw(2) << row << ". Number of parcles in interval\n"; ++row; }
      File << "#  " << std::setw(2) << row << ". Number of particles in interval\n"; ++row;
      File << "#  " << std::setw(2) << row << ". Total particle mass in interval\n"; ++row;
      File << "#  " << std::setw(2) << row << ". Total particle volume in interval\n"; ++row;
    }

    File << "# \n";
    File << "# \n";
    File << "# Total below smallest interval: \n";
    if (m_is_pic) { File << "#   Parcels:   " << static_cast<int>(a_data[0*bins]) << '\n'; }
    File << "#   Particles: " << static_cast<int>(a_data[1*bins]) << '\n';
    File << "#   Mass:      " << a_data[2*bins] << '\n';
    File << "#   Volume:    " << a_data[3*bins] << '\n';
    File << "# \n";
    File << "# \n";
    File << "# Total above largest interval: \n";
    if (m_is_pic) { File << "#   Parcels:   " << static_cast<int>(a_data[1*bins-1]) << '\n'; }
    File << "#   Particles: " << static_cast<int>(a_data[2*bins-1]) << '\n';
    File << "#   Mass:      " << a_data[3*bins-1] << '\n';
    File << "#   Volume:    " << a_data[4*bins-1] << '\n';
    File << "# \n";
    File << "# \n";

    for (int bin(1); bin<bins-1; ++bin) {
      File << "  " << std::scientific << std::setw(12)
                << std::setprecision(5) << get_midpoint(a_rt, bin);
      File << (bin < bins-2 ? ',' : '\n');
    }

    for (int n(0); n<4; ++n) {
      if (!m_is_pic && (n==0) ) { continue; }
      for (int bin(1); bin<bins-1; ++bin) {
        if (n < 2) {
          File << "  " << std::setw(12) << static_cast<long>(a_data[bin + n*bins]);
        } else {
          File << "  " << std::scientific << std::setw(12)
                    << std::setprecision(5) << a_data[bin+n*bins];
        }
        File << (bin < bins-2 ? ',' : '\n');
      }
    }

    File.flush();
    File.close();
  }

  ParallelDescriptor::Barrier();

}
