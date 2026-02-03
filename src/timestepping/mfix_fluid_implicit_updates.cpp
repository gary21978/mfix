#include <AMReX_VisMF.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <mfix.H>
#include <mfix_run_on.H>
#include <mfix_mf_helpers.H>
#include <mfix_dem.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_eb.H>
#include <mfix_pic.H>
#include <mfix_solvers.H>

namespace EnthalpyImplicitUpdate {

struct Residue
{
  AMREX_GPU_HOST_DEVICE
  Residue (const int& i,
           const int& j,
           const int& k,
           const ThermoPropertyData& fluid_props,
           const Array4<const Real>& Xgk_array,
           const Real& hg,
           const Real& ep_ro_g,
           const Real& convection_coeff,
           const Real& energy_source,
           const Real& dt)
    : m_i(i)
    , m_j(j)
    , m_k(k)
    , m_fluid_props(fluid_props)
    , m_Xgk_array(Xgk_array)
    , m_hg(hg)
    , m_ep_ro_g(ep_ro_g)
    , m_convection_coeff(convection_coeff)
    , m_energy_source(energy_source)
    , m_dt(dt)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real operator() (const Real& Tg_arg)
  {
    Real hg_loc = m_fluid_props.enthalpy(Tg_arg, amrex::IntVect(m_i,m_j,m_k), m_Xgk_array);

    return m_ep_ro_g*(hg_loc - m_hg) + m_dt*m_convection_coeff*Tg_arg - m_dt*m_energy_source;
  }

  const int& m_i; const int& m_j; const int& m_k;
  const ThermoPropertyData& m_fluid_props;
  const Array4<const Real>& m_Xgk_array;
  const Real& m_hg;
  const Real& m_ep_ro_g;
  const Real& m_convection_coeff;
  const Real& m_energy_source;
  const Real& m_dt;
};

struct Gradient
{
  AMREX_GPU_HOST_DEVICE
  Gradient (const int& i,
            const int& j,
            const int& k,
            const ThermoPropertyData& fluid_props,
            const Array4<const Real>& Xgk_array,
            const Real& ep_ro_g,
            const Real& convection_coeff,
            const Real& dt)
    : m_i(i)
    , m_j(j)
    , m_k(k)
    , m_fluid_props(fluid_props)
    , m_Xgk_array(Xgk_array)
    , m_ep_ro_g(ep_ro_g)
    , m_convection_coeff(convection_coeff)
    , m_dt(dt)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real operator() (const Real& Tg_arg)
  {
    Real gradient = m_fluid_props.specificHeat(Tg_arg, amrex::IntVect(m_i,m_j,m_k), m_Xgk_array);

    return m_ep_ro_g*gradient + m_dt*m_convection_coeff;
  }

  const int& m_i; const int& m_j; const int& m_k;
  const ThermoPropertyData& m_fluid_props;
  const Array4<const Real>& m_Xgk_array;
  const Real& m_ep_ro_g;
  const Real& m_convection_coeff;
  const Real& m_dt;
};

} // end namespace EnthalpyImplicitUpdate

using namespace EnthalpyImplicitUpdate;
using namespace Solvers;

//
// Implicit solve for the intermediate velocity.
// Currently this means accounting for the implicit part of the fluid/particle
// momentum exchange
//
void mfix::
add_vel_src_implicit (Real const a_dt,
                       Vector<MultiFab*      > const& a_vel,
                       Vector<MultiFab const*> const& a_epg,
                       Vector<MultiFab const*> const& a_rho,
                       Vector<MultiFab const*> const& a_S_p,
                       Vector<MultiFab const*> const& a_S_c)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::add_vel_src_implicit");

  for (int lev(0); lev<nlev(); ++lev) {

    for (MFIter mfi(*a_vel[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // Tilebox
      Box bx = mfi.tilebox();

      Array4<Real      > const& vel = a_vel[lev]->array(mfi);

      Array4<Real const> const& epg = a_epg[lev]->const_array(mfi);
      Array4<Real const> const& rho = a_rho[lev]->const_array(mfi);

      Array4<Real const> const& S_p = a_S_p[lev]->const_array(mfi);
      Array4<Real const> const& S_c = a_S_c[lev]->const_array(mfi);

      amrex::ParallelFor(bx,[dt=a_dt,vel,epg,rho,S_c,S_p]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real rop  = epg(i,j,k)*rho(i,j,k);

        for (int idim(0); idim<AMREX_SPACEDIM; ++idim) {
          vel(i,j,k,idim) = (rop*vel(i,j,k,idim) + dt*S_c(i,j,k,idim))
                          / (rop + dt*S_p(i,j,k));
        }
      });
    }
    a_vel[lev]->FillBoundary(geom[lev].periodicity());
  }
}

void mfix::
add_enthalpy_txfr_implicit ( Real dt,
                             Vector<MultiFab*      > const& h_g_in,
                             Vector<MultiFab*      > const& T_g_in,
                             Vector<MultiFab const*> const& X_gk_in,
                             Vector<MultiFab const*> const& txfr_in,
                             Vector<MultiFab const*> const& rho_in,
                             Vector<MultiFab const*> const& ep_g_in)
{
  /*
     This adds both components of the drag term
     So the drag term we add is beta * (particle_velocity - fluid_velocity)
                              = drag(0:2) - drag(3) * fluid_velocity
  */

  BL_PROFILE("mfix::mfix_add_energy_txfr_implicit");

  const auto fluid_props = fluid.props.data<run_on>();

  InterphaseTxfrIndexes txfr_idxs;

  const int idx_energy_source_txfr = txfr_idxs.energy_source;
  const int idx_convection_coeff_txfr = txfr_idxs.convection_coeff;

  for (int lev = 0; lev < nlev(); lev++) {
    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(ep_g_in[lev]->Factory());
    const auto& flags = factory.getMultiEBCellFlagFab();

    for (MFIter mfi(*h_g_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // Tilebox
      Box bx = mfi.tilebox();

      Array4<Real const> const& txfr_array = txfr_in[lev]->const_array(mfi);
      Array4<Real const> const& ro_array   = rho_in[lev]->const_array(mfi);
      Array4<Real const> const& ep_array   = ep_g_in[lev]->const_array(mfi);

      const int fluid_is_a_mixture = fluid.isMixture();

      Array4<Real const> dummy_arr;

      Array4<Real      > const& hg_array  = h_g_in[lev]->array(mfi);
      Array4<Real      > const& Tg_array  = T_g_in[lev]->array(mfi);
      Array4<Real const> const& Xgk_array = fluid_is_a_mixture ? X_gk_in[lev]->const_array(mfi) : dummy_arr;

      auto const& flags_arr = flags.const_array(mfi);

      amrex::ParallelFor(bx,[dt,hg_array,Tg_array,txfr_array,ro_array,ep_array,
          fluid_props,Xgk_array,flags_arr,
          idx_energy_source_txfr, idx_convection_coeff_txfr,abstol=newton_abstol,
          reltol=newton_reltol,maxiter=newton_maxiter]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const int cell_is_covered = static_cast<int>(flags_arr(i,j,k).isCovered());

        if (!cell_is_covered) {

          const Real hg = hg_array(i,j,k);

          const Real energy_source = txfr_array(i,j,k,idx_energy_source_txfr);
          const Real convection_coeff = txfr_array(i,j,k,idx_convection_coeff_txfr);

          const Real epg_loc = ep_array(i,j,k);

          const Real ep_ro_g = epg_loc*ro_array(i,j,k);

          // ************************************************************
          // Newton-Raphson solver for solving implicit equation for
          // temperature
          // ************************************************************
          Residue residue(i, j, k, fluid_props, Xgk_array, hg, ep_ro_g, convection_coeff, energy_source, dt);

          Gradient gradient(i, j, k, fluid_props, Xgk_array, ep_ro_g, convection_coeff, dt);

          Real Tg_old = Tg_array(i,j,k);

          Real Tg_new(Tg_old);

          auto output = Newton::solve(Tg_new, residue, gradient, abstol, reltol, maxiter);

          if (output.iterations == -1) {

            amrex::Abort("mfix::mfix_add_energy_txfr_implicit\n!!!Newton solver did not converge!!!");
          }
          Tg_array(i,j,k) = Tg_new;

          Real hg_new = fluid_props.enthalpy(IntVect(i,j,k), Tg_array, Xgk_array);

          hg_array(i,j,k) = hg_new;
        }
      });
    }
  }
}
