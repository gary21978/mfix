#ifndef MFIX_SOLVERS_H_
#define MFIX_SOLVERS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_ParticleTile.H>

#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_run_on.H>

#include <cmath>
#include <functional>
#include <exception>
#include <string>

namespace Solvers { namespace Newton {

namespace FluidEnthalpy {

struct Residue
{
  AMREX_GPU_HOST_DEVICE
  Residue (const int& i,
           const int& j,
           const int& k,
           const int& cell_is_covered,
           const ThermoPropertyData& fluid_props,
           const amrex::Array4<const amrex::Real>& X_gk,
           const amrex::Real& h_g)
  : m_i(i)
  , m_j(j)
  , m_k(k)
  , m_cell_is_covered(cell_is_covered)
  , m_fluid_props(fluid_props)
  , m_X_gk(X_gk)
  , m_h_g(h_g)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (const amrex::Real& Tg_arg)
  {
    amrex::Real hg_loc = m_fluid_props.enthalpy(Tg_arg, amrex::IntVect(m_i,m_j,m_k), m_X_gk, m_cell_is_covered);

    return hg_loc - m_h_g;
  }

  const int& m_i; const int& m_j; const int& m_k;
  const int& m_cell_is_covered;
  const ThermoPropertyData& m_fluid_props;
  const amrex::Array4<const amrex::Real>& m_X_gk;
  const amrex::Real& m_h_g;
};

struct Gradient
{
  AMREX_GPU_HOST_DEVICE
  Gradient (const int& i,
            const int& j,
            const int& k,
            const ThermoPropertyData& fluid_props,
            const amrex::Array4<const amrex::Real>& X_gk)
  : m_i(i)
  , m_j(j)
  , m_k(k)
  , m_fluid_props(fluid_props)
  , m_X_gk(X_gk)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (const amrex::Real& Tg_arg)
  {
    amrex::Real gradient = m_fluid_props.specificHeat(Tg_arg, amrex::IntVect(m_i,m_j,m_k), m_X_gk);

    return gradient;
  }

  const int& m_i; const int& m_j; const int& m_k;
  const ThermoPropertyData& m_fluid_props;
  const amrex::Array4<const amrex::Real>& m_X_gk;
};

} // end namespace FluidEnthalpy

namespace SolidsEnthalpy {

struct Residue
{
  AMREX_GPU_HOST_DEVICE
  Residue (const int& i,
           const ThermoPropertyData& solids_props,
           const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
           const int idx_X_sn,
           const amrex::Real& hp)
    : m_i(i)
    , m_solids_props(solids_props)
    , m_ptile_data(ptile_data)
    , m_idx_X_sn(idx_X_sn)
    , m_hp(hp)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (const amrex::Real& Tp_arg)
  {
    amrex::Real hp_loc = m_solids_props.enthalpy(Tp_arg, m_i, m_ptile_data.m_runtime_rdata, m_idx_X_sn);

    return hp_loc - m_hp;
  }

  const int& m_i;
  const ThermoPropertyData& m_solids_props;
  const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& m_ptile_data;
  const int m_idx_X_sn;
  const amrex::Real& m_hp;
};

struct Gradient
{
  AMREX_GPU_HOST_DEVICE
  Gradient (const int& i,
            const ThermoPropertyData& solids_props,
            const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
            const int idx_X_sn)
    : m_i(i)
    , m_solids_props(solids_props)
    , m_ptile_data(ptile_data)
    , m_idx_X_sn(idx_X_sn)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real operator() (const amrex::Real& Tp_arg)
  {
    amrex::Real gradient = m_solids_props.specificHeat(Tp_arg, m_i, m_ptile_data.m_runtime_rdata, m_idx_X_sn);

    return gradient;
  }

  const int& m_i;
  const ThermoPropertyData& m_solids_props;
  const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& m_ptile_data;
  const int m_idx_X_sn;
};

} // end namespace SolidsEnthalpy

// Newton solver output information
struct Output
{
  AMREX_GPU_HOST_DEVICE
  Output(const int iter,
         const amrex::Real update,
         const amrex::Real residue)
    : iterations(iter)
    , update_norm(update)
    , residue_norm(residue)
  {}

  int iterations;
  amrex::Real update_norm;
  amrex::Real residue_norm;
};

// Solve for scalars
template<typename ResidueFunc, typename GradientFunc>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Output solve (amrex::Real& solution,
              ResidueFunc& R,
              GradientFunc& partial_R,
              const amrex::Real /*abs_tol*/,
              const amrex::Real rel_tol,
              const int max_iterations)
{
  amrex::Real gradient = partial_R(solution);
  amrex::Real residue  = R(solution);
  amrex::Real update(0.);

  int iter(0);
  const amrex::Real update_rel_tol = rel_tol*std::abs(solution);
  //const amrex::Real residue_rel_tol = rel_tol*std::abs(residue);

  // Do not do anything if we're already on the solution
  if (std::abs(residue) < 1.e-15) {
    return Output(0, 0, std::abs(residue));
  }

  if (std::abs(gradient) < 1.e-15) {
    amrex::Abort("Gradient in Newton solver is almost zero");
  }

  do {
    update = residue / gradient;
    solution -= update;
    residue = R(solution);
    gradient = partial_R(solution);

    ++iter;

    if (iter > max_iterations)
      return Output(-1, std::abs(update), std::abs(residue));

  } while(//(std::abs(residue) > residue_rel_tol) ||
          (std::abs(update) > update_rel_tol));

  return Output(iter, std::abs(update), std::abs(residue));
}

} } // end namespace Solvers::Newton


#endif
