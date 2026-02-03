#ifndef _MFIX_PC_UPDATES_K_H_
#define _MFIX_PC_UPDATES_K_H_

#include <mfix_des_K.H>
#include <mfix_pc.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_calc_cell.H>


//***************************************************************
// Second step: update particles' positions and velocities
//***************************************************************
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
part_momentum_update (MFIXParticleContainer::ParticleType& p,
                      const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& p_realarray,
                      const amrex::GpuArray<int*, SoAintData::count>& p_intarray,
                      const int& i,
                      const amrex::Real& dt,
                      const int& ntot,
                      amrex::Real* fc_ptr, amrex::Real* tow_ptr,
                      const amrex::RealVect& gravity,
                      const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& p_lo,
                      const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& p_hi,
                      const amrex::GpuArray<int,6>& lo_hi_bc,
                      const amrex::GpuArray<int,3>& periodic,
                      const int implicit_drag, const amrex::Real& fluid_dt,
                      int const include_vm,  int const idx_pc_vm_coeff,
                      int const include_acc, int const idx_pc_acc)
{
  const amrex::Real p_radius = p_realarray[SoArealData::radius][i];
  const amrex::Real p_mass = p_realarray[SoArealData::density][i]*SoArealData::volume(p_radius);

  const amrex::Real p_velx_old = p_realarray[SoArealData::velx][i];
  const amrex::Real p_vely_old = p_realarray[SoArealData::vely][i];
  const amrex::Real p_velz_old = p_realarray[SoArealData::velz][i];

  amrex::Real p_velx_new, p_vely_new, p_velz_new;

  amrex::Real Sc_x = p_mass*gravity[0] + fc_ptr[i+0*ntot] + p_realarray[SoArealData::vel_source_x][i];
  amrex::Real Sc_y = p_mass*gravity[1] + fc_ptr[i+1*ntot] + p_realarray[SoArealData::vel_source_y][i];
  amrex::Real Sc_z = p_mass*gravity[2] + fc_ptr[i+2*ntot] + p_realarray[SoArealData::vel_source_z][i];

  amrex::Real Sp(0.);
  amrex::Real coeff(p_mass);

  //  drag forces
  if (implicit_drag) {
    if (include_vm) { coeff += ptile_data.m_runtime_rdata[idx_pc_vm_coeff][i]; }
    Sp = p_realarray[SoArealData::drag_coeff][i];
  }

  p_velx_new = (coeff*p_velx_old + dt*Sc_x)/(coeff + dt*Sp);
  p_vely_new = (coeff*p_vely_old + dt*Sc_y)/(coeff + dt*Sp);
  p_velz_new = (coeff*p_velz_old + dt*Sc_z)/(coeff + dt*Sp);

  if (include_acc) {
      ptile_data.m_runtime_rdata[idx_pc_acc  ][i] = dt*(p_velx_new-p_velx_old)/fluid_dt;
      ptile_data.m_runtime_rdata[idx_pc_acc+1][i] = dt*(p_vely_new-p_vely_old)/fluid_dt;
      ptile_data.m_runtime_rdata[idx_pc_acc+2][i] = dt*(p_velz_new-p_velz_old)/fluid_dt;
  }

  const amrex::Real p_omegax_old = p_realarray[SoArealData::omegax][i];
  const amrex::Real p_omegay_old = p_realarray[SoArealData::omegay][i];
  const amrex::Real p_omegaz_old = p_realarray[SoArealData::omegaz][i];

  const amrex::Real p_oneOverI = SoArealData::oneOverI(p_radius, p_mass);
  amrex::Real p_omegax_new = p_omegax_old + dt * p_oneOverI * tow_ptr[i];
  amrex::Real p_omegay_new = p_omegay_old + dt * p_oneOverI * tow_ptr[i+ntot];
  amrex::Real p_omegaz_new = p_omegaz_old + dt * p_oneOverI * tow_ptr[i+2*ntot];

  const amrex::Real p_posx_old = p.pos(0);
  const amrex::Real p_posy_old = p.pos(1);
  const amrex::Real p_posz_old = p.pos(2);

  amrex::Real p_posx_new = p_posx_old + dt * p_velx_new;
  amrex::Real p_posy_new = p_posy_old + dt * p_vely_new;
  amrex::Real p_posz_new = p_posz_old + dt * p_velz_new;

  if (p_intarray[SoAintData::state][i] != state::entering) {

    if (p_posx_new < p_lo[0]) {
      if (lo_hi_bc[0]) {
          p_posx_new = p_lo[0] + std::numeric_limits<amrex::Real>::epsilon();
          p_velx_new = -p_velx_new;
      } else if (!periodic[0] ) {
        p.id().make_invalid();
      }
    }
    if ( p_posx_new > p_hi[0]) {
      if (lo_hi_bc[1]) {
          p_posx_new = p_hi[0] - std::numeric_limits<amrex::Real>::epsilon();
          p_velx_new = -p_velx_new;
      } else if (!periodic[0] ) {
        p.id().make_invalid();
      }
    }

    if (p_posy_new < p_lo[1]) {
      if (lo_hi_bc[2]) {
        p_posy_new = p_lo[1] + std::numeric_limits<amrex::Real>::epsilon();
        p_vely_new = -p_vely_new;
      } else if (!periodic[1] ) {
        p.id().make_invalid();
      }
    }

    if (p_posy_new > p_hi[1]) {
      if (lo_hi_bc[3]) {
        p_posy_new = p_hi[1] - std::numeric_limits<amrex::Real>::epsilon();
        p_vely_new = -p_vely_new;
      } else if (!periodic[1] ) {
        p.id().make_invalid();
      }
    }

    if (p_posz_new < p_lo[2]) {
      if (lo_hi_bc[4]) {
        p_posz_new = p_lo[2] + std::numeric_limits<amrex::Real>::epsilon();
        p_velz_new = -p_velz_new;
      } else if (!periodic[2] ) {
        p.id().make_invalid();
      }
    }
    if (p_posz_new > p_hi[2]) {
      if (lo_hi_bc[5]) {
        p_posz_new = p_hi[2] - std::numeric_limits<amrex::Real>::epsilon();
        p_velz_new = -p_velz_new;
      } else if (!periodic[2] ) {
        p.id().make_invalid();
      }
    }

  } else {

    p.pos(0) += dt * p_velx_old + fc_ptr[i         ];
    p.pos(1) += dt * p_vely_old + fc_ptr[i +   ntot];
    p.pos(2) += dt * p_velz_old + fc_ptr[i + 2*ntot];

  }

  if (p.id().is_valid()) {

    // Update positions
    p.pos(0) = p_posx_new;
    p.pos(1) = p_posy_new;
    p.pos(2) = p_posz_new;

    // Update velocities
    p_realarray[SoArealData::velx][i] = p_velx_new;
    p_realarray[SoArealData::vely][i] = p_vely_new;
    p_realarray[SoArealData::velz][i] = p_velz_new;

    // Update angular velocities
    p_realarray[SoArealData::omegax][i] = p_omegax_new;
    p_realarray[SoArealData::omegay][i] = p_omegay_new;
    p_realarray[SoArealData::omegaz][i] = p_omegaz_new;

  }
}


//***************************************************************
// Third step: update particles' temperature
//***************************************************************
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
part_enthalpy_update (const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                      const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& p_realarray,
                      const int& i,
                      const int idx_X_sn,
                      const int /*idx_conv_coeff*/,
                      const int idx_energy_src,
                      const ThermoPropertyData& solids_props,
                      const amrex::Real& dt,
                      amrex::Real* cond_ptr,
                      const amrex::Real& enthalpy_source,
                      const amrex::Real& abstol,
                      const amrex::Real& reltol,
                      const int& maxiter,
                      const int& solve_dem)
{
  const amrex::Real p_radius = p_realarray[SoArealData::radius][i];
  const amrex::Real p_mass = p_realarray[SoArealData::density][i]*SoArealData::volume(p_radius);

  const amrex::Real Tp = p_realarray[SoArealData::temperature][i];
  amrex::Real p_enthalpy_old = solids_props.enthalpy(Tp, i, ptile_data.m_runtime_rdata, idx_X_sn);

  amrex::Real p_enthalpy_new(p_enthalpy_old);

  // convective heat transfer coefficient
  //amrex::Real const gamma = ptile_data.m_runtime_rdata[idx_conv_coeff][i];

  // constant source terms: gamma*(Tf - Tp) + ...
  amrex::Real const Sc = ptile_data.m_runtime_rdata[idx_energy_src][i];

  if (solve_dem) {
    p_enthalpy_new += dt*((Sc + cond_ptr[i] + enthalpy_source) / p_mass);
  } else {
    p_enthalpy_new += dt*((Sc + enthalpy_source) / p_mass);
  }

  // ************************************************************
  // Solvers::Newton-Raphson solver for solving implicit equation for
  // temperature
  // ************************************************************
  Solvers::Newton::SolidsEnthalpy::Residue residue(i, solids_props, ptile_data, idx_X_sn, p_enthalpy_new);

  Solvers::Newton::SolidsEnthalpy::Gradient gradient(i, solids_props, ptile_data, idx_X_sn);

  amrex::Real Tp_new(p_realarray[SoArealData::temperature][i]);

  auto output = Solvers::Newton::solve(Tp_new, residue, gradient, abstol, reltol, maxiter);

  if (output.iterations == -1)
    amrex::Abort("Newton solver did not converge");

  p_realarray[SoArealData::temperature][i] = Tp_new;
}

#endif
