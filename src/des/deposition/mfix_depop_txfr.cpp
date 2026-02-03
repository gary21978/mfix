#include <AMReX.H>
//#include <AMReX_Particles.H>
#include <AMReX_ParmParse.H>

#include <mfix_deposition_op.H>
#include <mfix_deposition_K.H>
#include <mfix_reporter.H>

#include <mfix_pc.H>

using namespace amrex;

MFIXDepOpTxfr::
MFIXDepOpTxfr ( BCList const& a_bc_list,
                MFIXParticleContainer* a_pc,
                const EBFArrayBoxFactory* a_lev0_factory,
                amrex::Vector <amrex::MultiFab*> a_dst,
                amrex::Real const a_dt )
  : MFIXDepositionOp(a_bc_list, a_pc, a_lev0_factory, a_dst)
  , m_dt( a_dt )
{ }


void
MFIXDepOpTxfr::
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
void MFIXDepOpTxfr::
Deposit (int const a_lev, F WeightFunc )
{
  BL_PROFILE("MFIXDepOpTxfr::Deposit()");

  int const lev0(0);

  // We always use the coarse dx
  const auto plo = m_pc->Geom(lev0).ProbLoArray();
  const auto dx  = m_pc->Geom(lev0).CellSizeArray();
  const auto dxi = m_pc->Geom(lev0).InvCellSizeArray();

  const auto reg_cell_vol = dx[0]*dx[1]*dx[2];

  const int solve_enthalpy = m_pc->get_fluid().solve_enthalpy();

  InterphaseTxfrIndexes txfr_idxs;

  int const idx_drag_coeff = txfr_idxs.drag_coeff;
  int const idx_vm_coeff   = txfr_idxs.vm_coeff;

  int const idx_vel_src_c  = txfr_idxs.vel_src_c;

  int const idx_convection_coeff_txfr = txfr_idxs.convection_coeff;
  int const idx_energy_source_txfr = txfr_idxs.energy_source;

  int const include_vm = includeVirtualMass();

  // particle acceleration and virtual mass coefficient
  int const idx_pc_dUpdt( (include_vm ? m_pc->m_runtimeRealData.acceleration : -1) );
  int const idx_pc_vm_coeff( (include_vm ? m_pc->m_runtimeRealData.vm_coeff : -1) );

  // convection coefficient
  int const idx_pc_conv_coeff( (solve_enthalpy ? m_pc->m_runtimeRealData.conv_coeff : -1) );

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

  {
    FArrayBox local_eps;
    FArrayBox local_txfr;

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

      FArrayBox& eps_fab  = (*m_eps[a_lev])[pti];
      FArrayBox& txfr_fab = (*m_ptr[a_lev])[pti];

      const Box& box = pti.tilebox(); // I need a box without ghosts

      if ((*flags)[pti].getType(box) != FabType::covered) {

        auto        txfr_arr = txfr_fab.array();
        auto         eps_arr = eps_fab.array();
        const auto& flagsarr = (*flags)[pti].const_array();
        const auto&    vfrac = (*volfrac)[pti].const_array();

        const Real deposition_scale_factor = m_deposition_scale_factor;



        const auto local_cg_dem = m_pc->get_dem().cg_dem();

        ParallelFor(nrp, [WeightFunc, pstruct, p_realarray, ptile_data, txfr_arr,
             plo, dx, dxi, reg_cell_vol, flagsarr, vfrac, eps_arr, solve_enthalpy,
             deposition_scale_factor, include_vm, idx_drag_coeff, idx_vm_coeff,
             idx_vel_src_c, idx_energy_source_txfr, idx_convection_coeff_txfr,
             idx_pc_dUpdt, idx_pc_vm_coeff, idx_pc_conv_coeff, local_cg_dem, dt=m_dt]
          AMREX_GPU_DEVICE (int ip) noexcept
        {
          const ParticleType& p = pstruct[ip];

          int i;
          int j;
          int k;

          Real pgamma(0.);
          Real vm_coeff(0.);
          RealVect dUpdt = {0.,0.,0.};

          const Real statwt = p_realarray[SoArealData::statwt][ip];

          GpuArray<GpuArray<GpuArray<Real,2>,2>,2> weights;

          WeightFunc(plo, dx, dxi, flagsarr, p.pos(), p_realarray[SoArealData::radius][ip],
                     i, j, k, weights, deposition_scale_factor);

          Real pvol  = statwt * SoArealData::volume(p_realarray[SoArealData::radius][ip]) / reg_cell_vol;
          Real pbeta = statwt * p_realarray[SoArealData::drag_coeff][ip] / reg_cell_vol;

          if (solve_enthalpy) {
            pgamma = statwt * ptile_data.m_runtime_rdata[idx_pc_conv_coeff ][ip] / reg_cell_vol;
          }

          if ( include_vm ) {
            // (rho_f * C_vm * vol_p)
            vm_coeff = ptile_data.m_runtime_rdata[idx_pc_vm_coeff][ip];

            // (rho_f * C_vm * vol_p) / vol_f ==>  (rho_f * C_vm * ep_s)
            vm_coeff *= (statwt / reg_cell_vol);

            // (rho_f * C_vm * ep_s) * (dUp/dt)
            dUpdt[0] = vm_coeff*ptile_data.m_runtime_rdata[idx_pc_dUpdt  ][ip];
            dUpdt[1] = vm_coeff*ptile_data.m_runtime_rdata[idx_pc_dUpdt+1][ip];
            dUpdt[2] = vm_coeff*ptile_data.m_runtime_rdata[idx_pc_dUpdt+2][ip];
          }

          if (local_cg_dem) {
            pvol = pvol / statwt;
            pbeta = pbeta / statwt;
          }

          Real pvx = p_realarray[SoArealData::velx][ip] * pbeta;
          Real pvy = p_realarray[SoArealData::vely][ip] * pbeta;
          Real pvz = p_realarray[SoArealData::velz][ip] * pbeta;

          Real pTp(0.);

          if (solve_enthalpy)
            pTp = p_realarray[SoArealData::temperature][ip] * pgamma;

          // Deposition
          for (int ii = -1; ii <= 0; ++ii) {
            for (int jj = -1; jj <= 0; ++jj) {
              for (int kk = -1; kk <= 0; ++kk) {

                if (flagsarr(i+ii,j+jj,k+kk).isCovered()) { continue; }

                Real weight_vol = weights[ii+1][jj+1][kk+1] / vfrac(i+ii,j+jj,k+kk);

                HostDevice::Atomic::Add(&eps_arr(i+ii,j+jj,k+kk), weight_vol*pvol);

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_drag_coeff), weight_vol*pbeta);

                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel_src_c+0), weight_vol*pvx);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel_src_c+1), weight_vol*pvy);
                HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel_src_c+2), weight_vol*pvz);

                if (include_vm) {
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vm_coeff), weight_vol*vm_coeff);

                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel_src_c+0), weight_vol*dUpdt[0]);
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel_src_c+1), weight_vol*dUpdt[1]);
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_vel_src_c+2), weight_vol*dUpdt[2]);
                }

                if (solve_enthalpy) {
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_energy_source_txfr   ), weight_vol*pTp);
                  HostDevice::Atomic::Add(&txfr_arr(i+ii,j+jj,k+kk,idx_convection_coeff_txfr), weight_vol*pgamma);
                }
              }
            }
          }
        });

      }
    }//pti
  }

  if (crse_factory != nullptr) { delete crse_factory; }
}
