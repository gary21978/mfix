#include <mfix.H>
#include <mfix_pc_updates_K.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>

using namespace amrex;
using namespace Solvers;

void MFIXParticleContainer::
MFIX_PC_AdvanceParcels ( Real dt,
                         LoadBalance* const a_loadbalance)
{

  BL_PROFILE("MFIXParticleContainer::MFIX_PC_AdvanceParcels()");

  const Real abstol = newton_abstol;
  const Real reltol = newton_reltol;
  const int maxiter = newton_maxiter;

  for (int lev(0); lev<nlev(); lev ++ )
  {

    for (MFIXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

      // Timer used for load-balancing
      Real const timer_start( ParallelDescriptor::second() );

      PairIndex index(pti.index(), pti.LocalTileIndex());

      const int nrp = GetParticles(lev)[index].numRealParticles();

      auto& plev = GetParticles(lev);
      auto& ptile = plev[index];
      auto& aos   = ptile.GetArrayOfStructs();
      ParticleType* pstruct = aos().dataPtr();

      auto ptile_data = ptile.getParticleTileData();

      auto& soa = ptile.GetStructOfArrays();
      auto p_realarray = soa.realarray();
      auto p_intarray  = soa.intarray();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR("pic_time_march()", pic_time_march);
#endif
      /********************************************************************
       * Move particles based on collision forces and torques             *
       *******************************************************************/

      const int nspecies_s = m_solids.nspecies();

      // Particles SoA indices
      const int idx_X_sn = m_runtimeRealData.X_sn;

      const int idx_conv_coeff = m_runtimeRealData.conv_coeff;
      const int idx_energy_src = m_runtimeRealData.energy_source;

      const int solve_enthalpy = m_solids.solve_enthalpy();

      const Real enthalpy_source = m_solids.enthalpy_source();

      const int solid_is_a_mixture = m_solids.isMixture();

      const auto solids_props = m_solids.props.data<run_on>();

      amrex::ParallelFor(nrp,
          [pstruct,p_realarray,p_intarray,ptile_data,dt,nspecies_s,idx_X_sn,
           solid_is_a_mixture,solve_enthalpy,enthalpy_source,idx_conv_coeff,
           idx_energy_src,solids_props, abstol,reltol,maxiter]
        AMREX_GPU_DEVICE (int i) noexcept
      {
        //*********************************************************************
        // First step: update parcels' temperature
        //*********************************************************************
        if (solve_enthalpy) {

          part_enthalpy_update(ptile_data, p_realarray, i, idx_X_sn,
              idx_conv_coeff, idx_energy_src, solids_props, dt, nullptr,
              enthalpy_source, abstol, reltol, maxiter, 0);
        }
      });

      /********************************************************************
       * Update runtime cost (used in load-balancing)                     *
       *******************************************************************/
      if (a_loadbalance->WeightByParticleRunTime()) {

        Real const weight( ParallelDescriptor::second() - timer_start );
        a_loadbalance->UpdateWeights(lev, pti.index(), weight);

      } else if (a_loadbalance->WeightByParticleCount()) {

        Real const weight(pti.numParticles() );
        a_loadbalance->UpdateWeights(lev, pti.index(), weight);

      }

      Gpu::synchronize();

#ifndef AMREX_USE_GPU
      BL_PROFILE_VAR_STOP(pic_time_march);
#endif

    } // particle-tile iterator

  } // loop over levels

}
