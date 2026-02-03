#include <AMReX.H>

#include <mfix_deposition_op.H>
#include <mfix_deposition_K.H>
#include <mfix_reporter.H>

#include <mfix_pc.H>

using namespace amrex;

MFIXDepOpVolAvg::
MFIXDepOpVolAvg ( BCList const& a_bc_list,
                  MFIXParticleContainer* a_pc,
                  const EBFArrayBoxFactory* a_lev0_factory,
                  amrex::Vector <amrex::MultiFab*> a_dst )
  : MFIXDepositionOp(a_bc_list, a_pc, a_lev0_factory, a_dst)
{
  m_h_depComp.resize(SoArealData::count + a_pc->m_runtimeRealData.count, 0 );
  m_d_depComp.resize(SoArealData::count + a_pc->m_runtimeRealData.count, 0 );
}


void
MFIXDepOpVolAvg::
setDepComp ( int const a_comp )
{ AMREX_ASSERT( a_comp < SoArealData::count + m_pc->m_runtimeRealData.count );

  if (!m_h_depComp[a_comp]) {
    m_h_depComp[a_comp] = 1;
    m_ncomp++;
  }
}


void
MFIXDepOpVolAvg::
setDepComp ( amrex::Gpu::HostVector<int> const& a_h_depComp )
{ AMREX_ASSERT( a_h_depComp.size() == m_h_depComp.size() );
  AMREX_ASSERT( m_ncomp == 0 );

  for ( size_t comp(0); comp < a_h_depComp.size(); ++comp) {
    AMREX_ASSERT( m_h_depComp[comp] == 0 );
    if (a_h_depComp[comp] == 1) {
      m_h_depComp[comp] = 1;
      m_ncomp++;
    }
  }
}


void
MFIXDepOpVolAvg::
Deposit ()
{

  int const lev0(0);

  for (int lev(0); lev<m_pc->numLevels(); ++lev) {

    m_eps[lev]->setVal(0.0, 0, m_eps[lev]->nComp(), m_eps[lev]->nGrow());
    m_ptr[lev]->setVal(0.0, 0, m_ptr[lev]->nComp(), m_eps[lev]->nGrow());

    if (m_deposition_scheme == DepositionScheme::trilinear) {

      Deposit( lev, TrilinearDeposition());

    } else if (m_deposition_scheme == DepositionScheme::square_dpvm) {

      Deposit( lev, TrilinearDPVMSquareDeposition());

    } else if (m_deposition_scheme == DepositionScheme::true_dpvm) {

      Deposit( lev, TrueDPVMDeposition());

    } else if (m_deposition_scheme == DepositionScheme::centroid) {

      Deposit( lev, CentroidDeposition());

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

    // Copy any data on fine grids to lev0
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
  }//lev

  m_ptr[lev0]->setBndry(0.0);
  m_eps[lev0]->FillBoundary(m_pc->Geom(lev0).periodicity());

}



template <typename F>
void MFIXDepOpVolAvg::
Deposit (int const a_lev, F WeightFunc )
{
  BL_PROFILE("MFIXDepOpVolAvg::Deposit()");

  if (m_ncomp == 0) {
    m_deposit_only_volume = 1;
    m_ncomp = 1;
  }

  Gpu::copyAsync(Gpu::hostToDevice, m_h_depComp.begin(),
      m_h_depComp.end(), m_d_depComp.begin());

  int const ncomp( m_ncomp );
  AMREX_ASSERT( ncomp <= m_ptr[a_lev]->nComp() );

  int const SoAcount( SoArealData::count );
  int const RTRcount( m_pc->m_runtimeRealData.count );

  int const lev0(0);

  // We always use the coarse dx
  const auto plo = m_pc->Geom(lev0).ProbLoArray();
  const auto dx  = m_pc->Geom(lev0).CellSizeArray();
  const auto dxi = m_pc->Geom(lev0).InvCellSizeArray();

  const auto reg_cell_vol = dx[0]*dx[1]*dx[2];

  auto& plev = m_pc->GetParticles(a_lev);

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

  // Make sure m_d_depComp is copied
  Gpu::synchronize();

  {
    FArrayBox local_eps;
    FArrayBox local_ptr;

    for (MFIXParIter pti(*m_pc, a_lev); pti.isValid(); ++pti) {

      MFIXParticleContainer::PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& ptile = plev[index];
      auto ptile_data = ptile.getParticleTileData();

      const auto& particles = pti.GetArrayOfStructs();
      const ParticleType* pstruct = particles().dataPtr();

      auto& soa = pti.GetStructOfArrays();
      auto p_realarray = soa.realarray();

      const long nrp = pti.numParticles();

      FArrayBox dummy_fab;

      FArrayBox& eps_fab = (*m_eps[a_lev])[pti];
      FArrayBox& ptr_fab = (*m_ptr[a_lev])[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered) {

        auto ptr_arr = ptr_fab.array();
        auto eps_arr = eps_fab.array();

        const auto& flagsarr = (*flags)[pti].const_array();
        const auto&    vfrac = (*volfrac)[pti].const_array();

        Real const deposition_scale_factor = m_deposition_scale_factor;

        int* depComp = m_d_depComp.data();
        int deposit_only_volume = m_deposit_only_volume;

        ParallelFor(nrp, [WeightFunc, pstruct, p_realarray, ptile_data, ptr_arr, eps_arr,
             plo, dx, dxi, reg_cell_vol, flagsarr, vfrac, deposition_scale_factor,
             ncomp, depComp, SoAcount, RTRcount, deposit_only_volume ]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          int i;
          int j;
          int k;

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          Real const radius = p_realarray[SoArealData::radius][ip];

          WeightFunc(plo, dx, dxi, flagsarr, p.pos(), radius,
              i, j, k, weights, deposition_scale_factor);

          Real const statwt = p_realarray[SoArealData::statwt][ip];
          Real const volume = SoArealData::volume(p_realarray[SoArealData::radius][ip]);

          Real const pvol = statwt * volume / reg_cell_vol;

          // Deposition
          for (int ii = -1; ii <= 0; ++ii) {
            for (int jj = -1; jj <= 0; ++jj) {
              for (int kk = -1; kk <= 0; ++kk) {

                if (flagsarr(i+ii,j+jj,k+kk).isCovered()) { continue; }

                Real const weight = (weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk));
                Real const weighted_vol = weight*pvol;

                HostDevice::Atomic::Add(&eps_arr(i+ii,j+jj,k+kk), weighted_vol);

                int dcomp(0);

                if (deposit_only_volume) {

                  HostDevice::Atomic::Add(&ptr_arr(i+ii,j+jj,k+kk, 0), weighted_vol);

                } else {

                  for ( int n(0); n<SoAcount+RTRcount; ++n) {

                    if ( depComp[n] ) {

                      Real attribute(weighted_vol);

                      if (n < SoAcount) {
                        attribute *= p_realarray[n][ip];
                      } else {
                        attribute *= ptile_data.m_runtime_rdata[n-SoAcount][ip];
                      }

                      HostDevice::Atomic::Add(&ptr_arr(i+ii,j+jj,k+kk, dcomp), attribute);

                      dcomp++;

                    } // if depComp
                  } // SoAcount + RTRcount
                } // deposit_only_volume
              }
            }
          }
        });

      }
    }//pti
  }

  if (crse_factory != nullptr) { delete crse_factory; }
}


void
MFIXDepOpVolAvg::
Average ( )
{
  for (int lev(0); lev < m_ptr.size(); ++lev) {
    Average(lev, 0, m_ptr[lev]->nComp() );
  }
}


void
MFIXDepOpVolAvg::
Average ( int const a_lev, int const a_scomp, int const a_ncomp)
{
  AMREX_ASSERT( a_scomp >= 0 );
  AMREX_ASSERT( m_ptr[a_lev]->nComp() >= a_scomp+a_ncomp );

  for (MFIter mfi(*m_ptr[a_lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    Box const& bx = mfi.tilebox();

    Array4<Real const> const& eps = m_eps[a_lev]->const_array(mfi);
    Array4<Real      > const& ptr = m_ptr[a_lev]->array(mfi);

    ParallelFor(bx, [eps, ptr, scomp=a_scomp, ncomp=a_ncomp]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      Real const inv_eps = (eps(i,j,k) > std::numeric_limits<Real>::min())
                         ? (1.0 / eps(i,j,k) ) : 0.0;

      for ( int n(scomp); n<ncomp; ++n) { ptr(i,j,k,n) *= inv_eps; }
    });
  }
}
