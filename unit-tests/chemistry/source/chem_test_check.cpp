#include <AMReX_ParmParse.H>

#include <chem_test_check.H>

using namespace amrex;

chem_test_check::
chem_test_check ( int  const a_solve_species,
                  int  const a_solve_enthalpy,
                  int  const a_nspecies,
                  int  const a_is_mixture,
                  int  const a_write_particles,
                  Real const a_vol,
                  const ThermoPropertyData& a_fluid_props,
                  const ThermoPropertyData& a_solids_props)
  : m_solve_species(a_solve_species)
  , m_solve_enthalpy(a_solve_enthalpy)
  , m_nspecies(a_nspecies)
  , m_is_mixture(a_is_mixture)
  , m_write_fluid(1)
  , m_write_particles(a_write_particles)
  , m_vol(a_vol)
  , m_fluid_props(a_fluid_props)
  , m_solids_props(a_solids_props)
  , m_verbose(0)
  , m_steps(0)
  , m_sum_abs_error(0.)
  , m_max_abs_error(0.)
{

  ParmParse pp("test");
  pp.query("verbose", m_verbose);

#if defined(CHEM_TEST_EULERIAN01)
  Print() << "CHEM_TEST_EULERIAN01 .............................................. ";
#endif

  if (m_verbose) { Print() << '\n'; }

}

void
chem_test_check::
write ( Real const a_time,
        Vector< MultiFab* > const& a_field_data,
        MFIXParticleContainer* a_particles)
{
  amrex::ignore_unused(a_time);
  amrex::ignore_unused(a_field_data);
  amrex::ignore_unused(a_particles);

  if (!m_verbose) { return; }

  Print() << std::setprecision(8) << std::scientific << a_time << ',';

  DualGridAuxIndexes fluid_idxs(m_solve_enthalpy, m_nspecies, m_is_mixture);

  int const nlev( a_field_data.size());
  for (int lev(0); lev<nlev; ++lev) {

    if (m_write_fluid) {

      for (MFIter mfi(*a_field_data[lev],false); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();

        Array4<Real const> const& epg = a_field_data[lev]->const_array(mfi, fluid_idxs.ep_g);

        Array4<Real const> const& rho = a_field_data[lev]->const_array(mfi, fluid_idxs.ro_g);

        Array4<Real const> const& Xgk = (!m_solve_species ? Array4<Real const>()
                                         : a_field_data[lev]->array(mfi, fluid_idxs.X_gk));

        Array4<Real const> const& hg = (!m_solve_enthalpy) ? Array4<Real const>()
                                     : a_field_data[lev]->const_array(mfi, fluid_idxs.h_g);

        Array4<Real const> const& Tg = (!m_solve_enthalpy) ? Array4<Real const>()
                                     : a_field_data[lev]->const_array(mfi, fluid_idxs.T_g);

        int const nspecies(m_solve_species ? m_nspecies : 0);

        amrex::ParallelFor(bx, [nspecies, solve_enthalpy=m_solve_enthalpy,
        epg, rho, Xgk, hg, Tg, vol=m_vol, props=m_fluid_props]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          printf("   %f,",rho(i,j,k));

          for (int n(0); n<nspecies; ++n) {
            printf("   %f,",Xgk(i,j,k,n));
          }

          if (solve_enthalpy) {

            Real Hg(0.);
            for (int n(0); n < nspecies; ++n) {
              Real const mass_n = epg(i,j,k)*rho(i,j,k)*Xgk(i,j,k,n)*vol; // kg
              Hg += props.enthalpy(n,Tg(i,j,k))*mass_n; // (J/kg)*(kg)
            }

            printf("   %f,   %f,",hg(i,j,k), Tg(i,j,k));
          }
        });
      }
    }

    if (m_write_particles) {

      for (MFIXParticleContainer::MFIXParIter pti(*a_particles, lev); pti.isValid(); ++pti) {

        auto& plev = a_particles->GetParticles(lev);

        // SoA real variables
        auto& soa = pti.GetStructOfArrays();

        auto radius = soa.GetRealData(SoArealData::radius     ).data();
        auto rho    = soa.GetRealData(SoArealData::density    ).data();
        auto Tp     = soa.GetRealData(SoArealData::temperature).data();

        std::pair<int, int> index(pti.index(), pti.LocalTileIndex());
        auto ptile_data = plev[index].getParticleTileData();

        int const nspecies( a_particles->m_runtimeRealData.count );
        int const X_p( a_particles->m_runtimeRealData.X_sn );

        amrex::ParallelFor(pti.numParticles(), [nspecies, solve_enthalpy=m_solve_enthalpy,
        radius, rho, Tp, ptile_data, X_p, props=m_solids_props]
        AMREX_GPU_DEVICE (int p) noexcept
        {
          Real mass = rho[p] * SoArealData::volume(radius[p]);
          printf("   %e,   %f,",mass, rho[p]);

          for (int n(0); n<nspecies; ++n) {
            printf("   %f,",ptile_data.m_runtime_rdata[X_p+n][p]);
          }

          if (solve_enthalpy) {

            Real hp = props.enthalpy(Tp[p], p, ptile_data.m_runtime_rdata, X_p);
            printf("   %f,   %f",hp, Tp[p]);
          }

        });

      }
    }
  }
  printf("\n");

}

Real chem_test_check::
get_cell_value ( MultiFab const* const& mf,
                 const IntVect& cell,
                 int const comp)
{
  Real bob(-1.);
  for (MFIter mfi(*mf); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.validbox();

    if (bx.contains(cell)) {

      auto const& fab = mf->const_array(mfi);
      Gpu::PinnedVector<amrex::Real> pv(1);
      auto* dp = pv.data();
      auto f = [=] AMREX_GPU_HOST_DEVICE ()
      { *dp = fab(cell, comp); };

#ifdef AMREX_USE_GPU
      if (mf->arena()->isManaged() || mf->arena()->isDevice()) {
        amrex::single_task(f);
        Gpu::streamSynchronize();
      } else
#endif
      {
        f();
      }
      bob = *dp;
    }
  }
  return bob;
}
