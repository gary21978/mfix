#include <mfix_bc.H>
#include <mfix_fluid.H>

using namespace amrex;


void
MFIXBoundaryConditions::
set_epf_bcs ( int const a_lev, Real a_time,
              MultiFab* a_epf,
              int const a_minf_type )
{
  BL_PROFILE("MFIXBoundaryConditions::set_epf_bcs()");

  AMREX_ASSERT( (a_minf_type == amrex::BCType::ext_dir)  ||
                (a_minf_type == amrex::BCType::foextrap) ||
                (a_minf_type == amrex::BCType::hoextrap) ||
                (a_minf_type == amrex::BCType::bogus) );

  // Set all values outside the domain to covered_val just to avoid use of
  // undefined
  a_epf->setDomainBndry(m_covered_val, m_geom[a_lev]);

  { Real* p_bc_epf = m_bc_epf.data();

    auto bcs_function = [p_bc_epf]
      AMREX_GPU_DEVICE (const int bct,
                        const int bcv,
                        const IntVect& ijk,
                        const IntVect& dom_ijk,
                        const int n,
                        const Array4<Real>& mf_arr,
                        const Array4<const EBCellFlag>& /*flags_arr*/,
                        const int /*dir*/)
    {
      if(bct == BCList::pinf || bct == BCList::pout) {
        mf_arr(ijk,n) = mf_arr(dom_ijk,n);
      }
      else if (bct == BCList::minf) {
        mf_arr(ijk,n) = p_bc_epf[bcv];
      }
    };

    set_bcs(a_lev, a_time, bcs_function, a_epf);
  }

  if ( a_minf_type != BCType::ext_dir ) {
    auto bcs_function = [a_minf_type]
      AMREX_GPU_DEVICE (const int bct,
                        const IntVect& ijk,
                        const IntVect& dom_ijk,
                        const IntVect& near_ijk,
                        const int n,
                        const Array4<Real>& mf_arr)
    {
      if (bct == BCList::minf) {
        if (a_minf_type == BCType::hoextrap) {
          mf_arr(ijk,n) = 2*mf_arr(ijk,n) - mf_arr(near_ijk,n);
        } else if (a_minf_type == BCType::foextrap) {
          mf_arr(ijk,n) = mf_arr(dom_ijk,n);
        } else if (a_minf_type == BCType::bogus) {
          mf_arr(ijk,n) = std::numeric_limits<Real>::max();
        }
      }
    };

    set_bcs_2D(a_lev, a_time, bcs_function, a_epf);
  }

  EB_set_covered(*a_epf, 0, a_epf->nComp(), a_epf->nGrow(), m_covered_val);

  // Do this after to pick up terms that got updated in the call above
  a_epf->FillBoundary(m_geom[a_lev].periodicity());
}


void
MFIXBoundaryConditions::
set_temperature_bcs ( int const a_lev, Real a_time,
                      const MFIXFluidPhase& fluid,
                      MultiFab* a_T)
{
  BL_PROFILE("MFIXBoundaryConditions::set_temperature_bcs()");

  set_energy_bc_values(a_time, fluid);
  Real* p_bc_T = m_bc_T.data();

  auto bcs_function = [p_bc_T]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int /*dir*/)
  {
    if(bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf || bct == BCList::pinf) {
      mf_arr(ijk,n) = p_bc_T[bcv];
    }
  };

  set_bcs(a_lev, a_time, bcs_function, a_T);
}


void
MFIXBoundaryConditions::
set_enthalpy_bcs ( int const a_lev, Real a_time,
                   const MFIXFluidPhase& a_fluid,
                   MultiFab* a_h)
{
  BL_PROFILE("MFIXBoundaryConditions::set_enthalpy_bcs()");

  const int nspecies = a_fluid.nspecies();
  const int is_mixture = a_fluid.isMixture();

  set_energy_bc_values(a_time, a_fluid);
  Real* p_bc_T = m_bc_T.data();

  Real* p_bc_X = nullptr;
  if (is_mixture) {
    set_species_bc_values(a_time, a_fluid.nspecies());
    p_bc_X = m_bc_Xk.data();
  }

  int const stride = m_bc.size();

  const auto fluid_props = a_fluid.props.data<run_on>();

  auto bcs_function = [p_bc_T, p_bc_X, fluid_props, nspecies, stride]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& flags_arr,
                      const int /*dir*/)
  {

    if(bct == BCList::pout) {

      mf_arr(ijk,n) = mf_arr(dom_ijk,n);

    } else if (bct == BCList::minf || bct == BCList::pinf) {

      const int is_covered = static_cast<int>(flags_arr(ijk).isCovered());

      if ( p_bc_X ) {

        Real h_sum(0.);

        for (int nn(0); nn < nspecies; nn++) {
          const Real h_k = fluid_props.enthalpy(nn, p_bc_T[bcv], is_covered);
          h_sum += p_bc_X[bcv + nn*stride]*h_k;
        }

        mf_arr(ijk,n) = h_sum;

      } else {

        mf_arr(ijk,n) = fluid_props.enthalpy(p_bc_T[bcv], nullptr, is_covered);

      }
    }
  };

  set_bcs(a_lev, a_time, bcs_function, a_h);
}


void
MFIXBoundaryConditions::
set_velocity_bcs ( int const a_lev, Real a_time,
                   MultiFab* a_vel,
                   int const a_minf_type )
{
  BL_PROFILE("MFIXBoundaryConditions::set_velocity_bcs()");

  AMREX_ASSERT( (a_minf_type == amrex::BCType::ext_dir) ||
                (a_minf_type == amrex::BCType::hoextrap) );

  // Set all values outside the domain to covered_val just to avoid use of
  // undefined
  a_vel->setDomainBndry(m_covered_val, m_geom[a_lev]);

  set_velocity_bc_values(a_time);

  const Real* p_bc_u_g = m_bc_u_g.data();
  const Real* p_bc_v_g = m_bc_v_g.data();
  const Real* p_bc_w_g = m_bc_w_g.data();

  {
    auto bcs_function = [p_bc_u_g,p_bc_v_g,p_bc_w_g]
      AMREX_GPU_DEVICE (const int bct,
                        const int bcv,
                        const IntVect& ijk,
                        const IntVect& dom_ijk,
                        const int n,
                        const Array4<Real>& mf_arr,
                        const Array4<const EBCellFlag>& /*flags_arr*/,
                        const int dir)
    {
      GpuArray<const Real*,3> p_bc_vel_g = {p_bc_u_g, p_bc_v_g, p_bc_w_g};

      if(bct == BCList::pinf || bct == BCList::pout) {
        mf_arr(ijk,n) = mf_arr(dom_ijk,n);
      } else if (bct == BCList::minf) {
        mf_arr(ijk,n) = (n == dir) ? p_bc_vel_g[dir][bcv] : 0.;
      }
    };

    set_bcs(a_lev, a_time, bcs_function, a_vel);
  }

  if ( a_minf_type == BCType::hoextrap ) {

    auto bcs_function = []
      AMREX_GPU_DEVICE (const int bct,
                        const IntVect& ijk,
                        const IntVect& /*dom_ijk*/,
                        const IntVect& near_ijk,
                        const int n,
                        const Array4<Real>& mf_arr)
    {
      if(bct == BCList::minf) {
        mf_arr(ijk,n) = 2*mf_arr(ijk,n) - mf_arr(near_ijk,n);
      }
    };

    set_bcs_2D(a_lev, a_time, bcs_function, a_vel);

  }
}


void
MFIXBoundaryConditions::
set_vec_bcs ( int const a_lev, Real a_time,
              MultiFab* mf_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_vec_bcs()");

  const Real* p_bc_u_g = m_bc_u_g.data();
  const Real* p_bc_v_g = m_bc_v_g.data();
  const Real* p_bc_w_g = m_bc_w_g.data();

  const Real* p_bc_epf = m_bc_epf.data();

  auto bcs_function = [p_bc_u_g,p_bc_v_g,p_bc_w_g,p_bc_epf]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int dir)
  {
    GpuArray<const Real*,3> p_bc_vel_g = {p_bc_u_g, p_bc_v_g, p_bc_w_g};

    if(bct == BCList::pinf || bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf) {
      mf_arr(ijk,n) = p_bc_epf[bcv] * p_bc_vel_g[dir][bcv];
    }
  };

  set_bcs(a_lev, a_time, bcs_function, mf_in);
}


void
MFIXBoundaryConditions::
set_tracer_bcs ( int const a_lev, Real a_time,
                 const MFIXFluidPhase& fluid,
                 MultiFab* trac_in)
{
  BL_PROFILE("MFIXBoundaryConditions::set_tracer_bcs()");

  if (!fluid.solve_tracer()) { return; }

//  set_tracer_bc_values(a_time, fluid);
  Real* p_bc_trac = m_bc_tracer.data();

  auto bcs_function = [p_bc_trac]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int /*dir*/)
  {
    if(bct == BCList::pinf || bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf) {
      mf_arr(ijk,n) = p_bc_trac[bcv];
    }
  };

  set_bcs(a_lev, a_time, bcs_function, trac_in);
}


void
MFIXBoundaryConditions::
set_species_bcs ( int const a_lev, Real a_time,
                  int const a_nspecies, MultiFab* a_Xk)
{
  BL_PROFILE("MFIXBoundaryConditions::set_species_bcs()");

  set_species_bc_values(a_time, a_nspecies);
  Real* p_bc_Xk = m_bc_Xk.data();

  int const stride = m_bc.size();

  auto bcs_function = [p_bc_Xk, stride]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int /*dir*/)
  {
    if(bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf || bct == BCList::pinf) {
      mf_arr(ijk,n) = p_bc_Xk[bcv + n*stride];
    }
  };

  set_bcs(a_lev, a_time, bcs_function, a_Xk);
}


void
MFIXBoundaryConditions::
set_density_bcs (int const a_lev, Real a_time,
                 MultiFab* a_rho)
{
  BL_PROFILE("MFIXBoundaryConditions::set_density_bcs()");

  set_density_bc_values(a_time);
  Real* p_bc_rho = m_bc_rho.data();

  auto bcs_function = [p_bc_rho]
    AMREX_GPU_DEVICE (const int bct,
                      const int bcv,
                      const IntVect& ijk,
                      const IntVect& dom_ijk,
                      const int n,
                      const Array4<Real>& mf_arr,
                      const Array4<const EBCellFlag>& /*flags_arr*/,
                      const int /*dir*/)
  {
    if(bct == BCList::pinf || bct == BCList::pout) {
      mf_arr(ijk,n) = mf_arr(dom_ijk,n);
    }
    else if (bct == BCList::minf) {
      mf_arr(ijk,n) = p_bc_rho[bcv];
    }
  };

  set_bcs(a_lev, a_time, bcs_function, a_rho);
}
