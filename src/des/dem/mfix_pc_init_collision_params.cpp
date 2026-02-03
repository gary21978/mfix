#include <mfix_bc.H>
#include <mfix_pc.H>
#include <mfix_dem.H>
#include <mfix_reporter.H>

using namespace amrex;


void MFIXParticleContainer::MFIX_PC_InitCollisionParams ()
{
  Vector<Real> max_dp, max_ro;
  Vector<Real> avg_dp, avg_ro;

  max_dp.resize(m_solids.ntypes());
  max_ro.resize(m_solids.ntypes());
  avg_dp.resize(m_solids.ntypes());
  avg_ro.resize(m_solids.ntypes());

  // Cycle through the different phases, starting from 1
  for (const int& phase: m_solids.phases())
  {
    Real h_pnum  = 0;   //number of particle
    Real h_pdiam = 0.0; //particle diameters
    Real h_pdens = 0.0; //particle density
    Real h_maxdiam = -1.e32;
    Real h_maxdens = -1.e32;

    for (int lev(0); lev<nlev(); lev++)
    {
      // Reduce sum operation for np, pdiam, pdens, maxdiam, maxdens
      ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMax, ReduceOpMax> reduce_op;
      ReduceData<Real, Real, Real, Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        const int np = particles.size();

        auto& soa = pti.GetStructOfArrays();
        auto p_realarray = soa.realarray();
        auto p_intarray = soa.intarray();

        reduce_op.eval(np, reduce_data, [p_realarray,p_intarray,phase]
            AMREX_GPU_DEVICE (int p_id) -> ReduceTuple
        {
          Real l_pnum  = 0._rt;
          Real l_pdiam = 0._rt;
          Real l_pdens = 0._rt;
          Real l_maxdiam = -1.e32;
          Real l_maxdens = -1.e32;

          if (phase == p_intarray[SoAintData::phase][p_id])
          {
            const Real density  = p_realarray[SoArealData::density][p_id];
            const Real diameter = 2.0*p_realarray[SoArealData::radius][p_id];

            l_pnum  = 1._rt;
            l_pdiam = diameter;
            l_pdens = density;
            l_maxdiam = diameter;
            l_maxdens = density;
          }

          return {l_pnum, l_pdiam, l_pdens, l_maxdiam, l_maxdens};
        });
      }

      ReduceTuple host_tuple = reduce_data.value();
      h_pnum  += amrex::get<0>(host_tuple);
      h_pdiam += amrex::get<1>(host_tuple);
      h_pdens += amrex::get<2>(host_tuple);
      h_maxdiam = amrex::max(h_maxdiam, amrex::get<3>(host_tuple));
      h_maxdens = amrex::max(h_maxdens, amrex::get<4>(host_tuple));
    }

    // A single MPI call passes all three variables
    ParallelDescriptor::ReduceRealSum({h_pnum,h_pdiam,h_pdens});
    ParallelDescriptor::ReduceRealMax({h_maxdiam,h_maxdens});

    //calculate averages or set = zero if no particles of that phase
    const int phase_idx = MFIXSolidsPhase::phase_to_index(phase);

    if (h_pnum==0) {
       avg_dp[phase_idx] = 0.0;
       avg_ro[phase_idx] = 0.0;

       max_dp[phase_idx] = 0.0;
       max_ro[phase_idx] = 0.0;

    } else {
       avg_dp[phase_idx] = h_pdiam/h_pnum;
       avg_ro[phase_idx] = h_pdens/h_pnum;

       max_dp[phase_idx] = h_maxdiam;
       max_ro[phase_idx] = h_maxdens;
    }
  }

  // Loop over BCs
  for (int bcv(0); bcv < m_boundary_conditions.bc().size(); ++bcv) {

    // EB flow with at least one solid
    if (m_boundary_conditions.bc(bcv).type == BCList::eb &&
        m_boundary_conditions.bc(bcv).solids.size() > 0) {

      const Real tolerance = std::numeric_limits<Real>::epsilon();

      for(int lcs(0); lcs < m_boundary_conditions.bc(bcv).solids.size(); lcs++) {

        const SOLIDS_t& bc_solid = m_boundary_conditions.bc(bcv).solids[lcs];

        if(bc_solid.volfrac > tolerance) {

          const Real mean_dp_bc = bc_solid.diameter.get_mean();
          const Real  max_dp_bc = bc_solid.diameter.get_max();

          const Real mean_rhop_bc = bc_solid.density.get_mean();
          const Real  max_rhop_bc = bc_solid.density.get_max();

          const int phase_idx = MFIXSolidsPhase::phase_to_index(bc_solid.phase);

          if ( avg_dp[phase_idx] == 0.0 ) { avg_dp[phase_idx] = mean_dp_bc; }
          if ( avg_ro[phase_idx] == 0.0 ) { avg_ro[phase_idx] = mean_rhop_bc; }

          max_dp[phase_idx] = amrex::max(max_dp[phase_idx], max_dp_bc);
          max_ro[phase_idx] = amrex::max(max_ro[phase_idx], max_rhop_bc);

        }
      }
    }
  }


  Real max_max_dp(0.);
  for (const int& phase: m_solids.phases())
  {
    const int phase_idx = MFIXSolidsPhase::phase_to_index(phase);

    if ((avg_dp[phase_idx] <= 0.0) || (max_dp[phase_idx] <= 0.0) ||
        (avg_ro[phase_idx] <= 0.0) || (max_ro[phase_idx] <= 0.0) ) {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
        << "Unable to determine DEM particle collision properties for "
        << m_solids.names(phase_idx)
        << "\ndiameter --"
        << "\n  average: " << avg_dp[phase_idx]
        << "\n  max:     " << max_dp[phase_idx]
        << "\ndensity --"
        << "\n  average: " << avg_ro[phase_idx]
        << "\n  max:     " << max_ro[phase_idx]
        << "\nPlease correct simulation settings.";
    }
     max_max_dp = amrex::max(max_max_dp, max_dp[phase_idx]);
  }

  // (3*max_dp/2)^2
  m_dem.set_neighborhood(2.25*max_max_dp*max_max_dp);

  // Polydisperse
  if (m_dem.pneig_flag()) {
      // Set up a the bin vector with type sizes
      int* ref_p = m_dem.prefratdata();
      auto& pbin = m_dem.get_pbin();
      for (int n(0); n<m_dem.nptypes(); ++n) pbin.push_back(max_max_dp / (amrex::Real) ref_p[n] );
      // Recursive loop to get smallest->largest ordering (e.g., 0-0 0-1 0-2, 1-1 1-2, 2-2 )
      amrex::Real* bin_p  = m_dem.pbindata();
      auto& pneighborhood = m_dem.get_pneighborhood();
      for (int i(0); i<m_dem.nptypes(); ++i) {
          for (int j(i); j<m_dem.nptypes(); ++j) {
              Real dist = 0.75 * ( bin_p[i] + bin_p[j] ); // 1.5 * (Rp1 + Rp2)
              pneighborhood.push_back(dist*dist);
          }
      }

      // Overwrite ptype based upon bin sizes
      Gpu::DeviceVector<Real> pbin_d(m_dem.nptypes());
      Gpu::copy(Gpu::hostToDevice, m_dem.pbindata(), m_dem.pbindata() + m_dem.nptypes(), pbin_d.begin());
      Real* bin_d = pbin_d.data();
      for (int lev(0); lev<nlev(); lev++) {
          for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti) {
              PairIndex index(pti.index(), pti.LocalTileIndex());
              auto& soa        = pti.GetStructOfArrays();
              auto p_realarray = soa.realarray();
              auto p_intarray  = soa.intarray();
              const int nrp    = GetParticles(lev)[index].numRealParticles();

              amrex::ParallelFor(nrp,
              [p_realarray,p_intarray,bin_d,lnptypes=m_dem.nptypes()]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                  // Bins are smallest -> largest
                  for (int n(0); n<lnptypes; ++n) {
                      Real diameter = 2.0 * p_realarray[SoArealData::radius][i];
                      if (diameter <= bin_d[n]) {
#if MFIX_POLYDISPERSE
                          p_intarray[SoAintData::ptype][i] = n;
#else
                          amrex::Abort("MFIX not built with POLYDISPERSE support");
#endif
                          break;
                      }
                  }
              }); // ParallelFor
          } // ParIter
      } // lev
  } // PneighFlag

  Real tcoll(1.0);

  MFIXDEM::A2D::array_type host_etan, host_etat, host_en;
  MFIXDEM::A1D::array_type host_etan_w, host_etat_w, host_en_w;
#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy_async(&host_en, m_dem.en().arrayPtr(), sizeof(MFIXDEM::A2D::array_type));
  amrex::Gpu::dtoh_memcpy_async(&host_en_w, m_dem.en_w().arrayPtr(), sizeof(MFIXDEM::A1D::array_type));
  amrex::Gpu::synchronize();
#else
  host_en = *(m_dem.en().arrayPtr());
  host_en_w = *(m_dem.en_w().arrayPtr());
#endif

  for (const int& phase: m_solids.phases()) {

    const int phase_idx = MFIXSolidsPhase::phase_to_index(phase);

    const Real dp_m = avg_dp[phase_idx];
    const Real mass_m = (M_PI/6.0)*(dp_m*dp_m*dp_m) * avg_ro[phase_idx];

    // Collision parameters between phase_idx M and all other phase_idxs
    for (int other_phase_idx(phase_idx); other_phase_idx < m_solids.ntypes(); ++other_phase_idx) {

      const Real dp_l = avg_dp[other_phase_idx];
      const Real mass_l = (M_PI/6.0)*(dp_l*dp_l*dp_l) * avg_ro[other_phase_idx];

      const Real mass_eff = (mass_m * mass_l) / (mass_m + mass_l);

      //Calculate the M-L normal and tangential damping coefficients
      host_etan(phase_idx,other_phase_idx) = 2.0*std::sqrt(m_dem.kn()*mass_eff);
      if(amrex::Math::abs(host_en(phase_idx,other_phase_idx)) > 0.0){
        const Real log_en = std::log(host_en(phase_idx,other_phase_idx));
        host_etan(phase_idx,other_phase_idx) *= amrex::Math::abs(log_en)/(std::sqrt(M_PI*M_PI + log_en*log_en));
      }
      host_etat(phase_idx,other_phase_idx) = m_dem.eta_fac() * host_etan(phase_idx,other_phase_idx);

      // Store the symmetric components
      host_etan(other_phase_idx,phase_idx) = host_etan(phase_idx,other_phase_idx);
      host_etat(other_phase_idx,phase_idx) = host_etat(phase_idx,other_phase_idx);

      // Collision time scale for M-L interactions
      Real tcoll_ml = M_PI/sqrt(m_dem.kn()/mass_eff -
         0.25*(host_etan(phase_idx,other_phase_idx)/mass_eff)*(host_etan(phase_idx,other_phase_idx)/mass_eff));

      tcoll = amrex::min(tcoll, tcoll_ml);
    }

    // Collision parameters between phase_idx M and wall
    {
      const Real mass_eff = mass_m;

      //Calculate the M-L normal and tangential damping coefficients
      host_etan_w(phase_idx) = 2.0*std::sqrt(m_dem.kn_w()*mass_eff);
      if(amrex::Math::abs(host_en_w(phase_idx)) > 0.0){
         const Real log_en = std::log(host_en_w(phase_idx));
         host_etan_w(phase_idx) *= amrex::Math::abs(log_en)/(std::sqrt(M_PI*M_PI + log_en*log_en));
      }
      host_etat_w(phase_idx) = m_dem.eta_w_fac() * host_etan_w(phase_idx);

      // Calculate the collision time scale.
      Real tcoll_w = M_PI/std::sqrt(m_dem.kn_w()/mass_eff -
         0.25*(host_etan_w(phase_idx)/mass_eff)*(host_etan_w(phase_idx)/mass_eff));

      tcoll = amrex::min(tcoll, tcoll_w);
    }
  }

#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy_async(m_dem.etan().arrayPtr(), &host_etan, sizeof(MFIXDEM::A2D::array_type));
  amrex::Gpu::htod_memcpy_async(m_dem.etat().arrayPtr(), &host_etat, sizeof(MFIXDEM::A2D::array_type));
  amrex::Gpu::htod_memcpy_async(m_dem.etan_w().arrayPtr(), &host_etan_w, sizeof(MFIXDEM::A1D::array_type));
  amrex::Gpu::htod_memcpy_async(m_dem.etat_w().arrayPtr(), &host_etat_w, sizeof(MFIXDEM::A1D::array_type));
  amrex::Gpu::synchronize();
#else
  *(m_dem.etan().arrayPtr()) = host_etan;
  *(m_dem.etat().arrayPtr()) = host_etat;
  *(m_dem.etan_w().arrayPtr()) = host_etan_w;
  *(m_dem.etat_w().arrayPtr()) = host_etat_w;
#endif

  ParmParse pp("mfix");
  Real tcoll_ratio = 50.;
  pp.query("tcoll_ratio", tcoll_ratio);

  m_dem.set_dtsolid(tcoll / tcoll_ratio);
}
