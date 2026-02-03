#include <restarter.H>
#include <deposition/mfix_deposition_K.H>
#include <mfix_indexes_aux.H>

#include <AMReX_ParmParse.H>

void MFIXDepOpPIC2DEM::
Deposit ( )
{
  int const lev0(0);

  for (int lev(0); lev < m_pc->numLevels(); ++lev) {

    if (m_deposition_scheme == DepositionScheme::trilinear) {

      Deposit(lev, TrilinearDeposition());

    } else if (m_deposition_scheme == DepositionScheme::square_dpvm) {

      Deposit(lev, TrilinearDPVMSquareDeposition());

    } else if (m_deposition_scheme == DepositionScheme::true_dpvm) {

      Deposit(lev, TrueDPVMDeposition());

    } else if (m_deposition_scheme == DepositionScheme::centroid) {

      Deposit(lev, CentroidDeposition());

    } else {

      amrex::Abort("Don't know this deposition_scheme!");
    }

    // Copy data deposited outside the domain back inside the domain
    // when BC is either a pressure inlet or mass inflow.
    ApplyBC(lev, m_ptr[lev]);
    ApplyBC(lev, m_eps[lev]);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    m_ptr[lev]->SumBoundary(m_pc->Geom(lev).periodicity());
    m_eps[lev]->SumBoundary(m_pc->Geom(lev).periodicity());

    // Copy any data on fine grids to level-0
    if (lev > lev0) {

      if ( m_reduce_to_lev0) {

        m_ptr[lev0]->ParallelCopy(*m_ptr[lev], 0, 0, m_ptr[lev]->nComp(), 0, 0,
          m_pc->Geom(lev).periodicity(), FabArrayBase::ADD);

        m_eps[lev0]->ParallelCopy(*m_eps[lev], 0, 0, m_eps[lev]->nComp(), 0, 0,
          m_pc->Geom(lev).periodicity(), FabArrayBase::ADD);

      } else {

        m_ptr[lev]->setBndry(0.0);
        m_eps[lev]->FillBoundary(m_pc->Geom(lev0).periodicity());

      }
    } // lev > lev0
  } // lev

  m_ptr[lev0]->setBndry(0.0);
  m_eps[lev0]->FillBoundary(m_pc->Geom(lev0).periodicity());

}


template <typename F>
void MFIXDepOpPIC2DEM::
Deposit ( int a_lev, F WeightFunc)
{
  BL_PROFILE("MFIXDepOpPIC2DEM::deposit()");

  int const lev0(0);

  // We always use the coarse dx
  const auto plo = m_pc->Geom(lev0).ProbLoArray();
  const auto dx  = m_pc->Geom(lev0).CellSizeArray();
  const auto dxi = m_pc->Geom(lev0).InvCellSizeArray();
  const Real dV = dx[0]*dx[1]*dx[2];

  const auto& solids = m_pc->get_solids();

  const int nspecies_s = solids.nspecies();

  Transfer txfr_idxs(solids);
  const int idx_eps     = txfr_idxs.idx_eps;
  const int idx_density = txfr_idxs.idx_density;
  const int idx_vel     = txfr_idxs.idx_vel;
  const int idx_temp    = txfr_idxs.idx_temp;
  const int idx_species = txfr_idxs.idx_species;

  const int idx_X_sn = m_pc->m_runtimeRealData.X_sn;

  EBFArrayBoxFactory* crse_factory = nullptr;

  const MultiFab* volfrac = nullptr;
  const FabArray<EBCellFlagFab>* flags = nullptr;

  if (a_lev == 0) {

    volfrac = &(m_lev0_factory->getVolFrac());
    flags   = &(m_lev0_factory->getMultiEBCellFlagFab());

  } else {

    crse_factory = (makeEBFabFactory(m_pc->Geom(0), m_pc->ParticleBoxArray(a_lev),
      m_pc->ParticleDistributionMap(a_lev), {1,1,1}, EBSupport::volume)).release();

    volfrac = &(crse_factory->getVolFrac());
    flags   = &(crse_factory->getMultiEBCellFlagFab());
  }

  {
    FArrayBox local_txfr;

    for (MFIXParIter pti(*m_pc, a_lev); pti.isValid(); ++pti) {

      PairIndex index(pti.index(), pti.LocalTileIndex());
      auto& ptile = m_pc->GetParticles(a_lev)[index];

      //Access to added variables
      auto ptile_data = ptile.getParticleTileData();

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();

      FArrayBox dummy_fab;

      FArrayBox& txfr_fab = (*m_ptr[a_lev])[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered) {

        auto txfr_arr = txfr_fab.array();

        const auto& flags_arr = (*flags)[pti].const_array();
        const auto& vfrac_arr = (*volfrac)[pti].const_array();

        const Real deposition_scale_factor = m_deposition_scale_factor;


        const int solve_enthalpy = solids.solve_enthalpy();
        const int solve_species  = solids.solve_species();

        amrex::ParallelFor(nrp, [pstruct,p_realarray,plo,dx,dxi,ptile_data,dV,idx_eps,
            deposition_scale_factor,WeightFunc,flags_arr,txfr_arr,solve_enthalpy,
            idx_vel,idx_temp,idx_species,nspecies_s,idx_X_sn,solve_species,vfrac_arr,
            idx_density]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          int i; int j; int k;

          const Real pradius = p_realarray[SoArealData::radius][ip];
          const Real pdensity = p_realarray[SoArealData::density][ip];

          const Real statwt  = p_realarray[SoArealData::statwt][ip];

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          WeightFunc(plo, dx, dxi, flags_arr, p.pos(), pradius, i, j, k, weights,
            deposition_scale_factor);

          const Real pvol = SoArealData::volume(p_realarray[SoArealData::radius][ip]);

          const Real pvel_x = p_realarray[SoArealData::velx][ip];
          const Real pvel_y = p_realarray[SoArealData::vely][ip];
          const Real pvel_z = p_realarray[SoArealData::velz][ip];

          Real ptemp(0.);
          if (solve_enthalpy)
            ptemp = p_realarray[SoArealData::temperature][ip];

          // Deposition
          for (int ii = -1; ii <= 0; ++ii) {
            for (int jj = -1; jj <= 0; ++jj) {
              for (int kk = -1; kk <= 0; ++kk) {
                if (flags_arr(i+ii,j+jj,k+kk).isCovered())
                  continue;

                Real weight = weights[ii+1][jj+1][kk+1] * (pvol/(vfrac_arr(i+ii,j+jj,k+kk)*dV));

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_eps), statwt*weight);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_density), statwt*weight*pdensity);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel+0), statwt*weight*pvel_x);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel+1), statwt*weight*pvel_y);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel+2), statwt*weight*pvel_z);

                if (solve_enthalpy) {
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_temp), statwt*weight*ptemp);
                }

                if (solve_species) {
                  for (int n_s(0); n_s < nspecies_s; ++n_s) {
                    HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_species+n_s),
                        statwt*weight*ptile_data.m_runtime_rdata[idx_X_sn+n_s][ip]);
                  }
                }
              }
            }
          }
        });

      }
    }
  }

  MultiFab::Copy(*m_eps[a_lev], *m_ptr[a_lev], idx_eps, 0, 1, m_eps[a_lev]->nGrow());
}
