#include <limits>

#include <mfix.H>
#include <mfix_run_on.H>
#include <mfix_diffusion_op.H>
#include <mfix_reporter.H>


void
mfix::compute_dt (Vector< MultiFab const* > const& a_avg_particle_data)
{
    // Store the dt we've just used in the previous time step as prev_dt
    m_timer.prev_dt() = m_timer.dt();

    // dt is always computed even when fixed_dt is set,
    // so we can issue a warning if the value of fixed dt does not satisfy the CFL condition.

    Real dt_new;
    /*
       Compute new dt by using the formula derived in
       "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
       by Kang et al. (JCP).

       dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )

      where

      C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection

      V = 2 * max(mu/ro) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion

      Fx, Fy, Fz = net acceleration due to external forces

      WARNING: We use a slightly modified version of C in the implementation below

    */
    bool explicit_diffusion = (m_predictor_diff_type == DiffusionType::Explicit);

    InterphaseTxfrIndexes txfr_idxs;

    const int idx_drag_coeff = txfr_idxs.drag_coeff;
    const int idx_vel_src_c = txfr_idxs.vel_src_c;

    const auto& fluid_parms = fluid.parameters<run_on>();

    diffOpVel()->setDiffCoeff(leveldata().vel_const(), leveldata().epf_const(),
        leveldata().rho_const(), leveldata().T_const(),
        a_avg_particle_data, leveldata().X_const(),
        /*include_eddy_viscosity=*/1);

    Vector< MultiFab const* > eff_visc = diffOpVel()->getDiffCoeff();

    // Max convective, diffusive and forcing CFL for all levels
    Real c_cfl = Real(0.0);
    Real v_cfl = Real(0.0);
    Real f_cfl = Real(0.0);
    for (int lev(0); lev < nlev(); ++lev) {

      const Real* dx = geom[lev].CellSize();

      Real odx(1.0 / dx[0]);
      Real ody(1.0 / dx[1]);
      Real odz(1.0 / dx[2]);

      const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(leveldata().rho(lev)->Factory());

      // Reduce max operation for c_cfl
      ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
      ReduceData<Real, Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for (MFIter mfi(*(leveldata().vel(lev)),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const auto& vel    = leveldata().vel(lev,mfi);
        const auto& ep     = leveldata().epf(lev,mfi);
        const auto& ro     = leveldata().rho_const(lev,mfi);
        const auto& gradp  = leveldata().grad_p_const(lev,mfi);
        const auto& txfr   = leveldata().txfr_const(lev,mfi);

        Box bx(mfi.tilebox());

        EBCellFlagFab const& flagfab = factory.getMultiEBCellFlagFab()[mfi];
        Array4<EBCellFlag const> const& flags = flagfab.const_array();

        // ew need this until we remove static attribute from mfix::gp0
        const RealVect gp0_dev(gp0);
        const RealVect gravity_dev(gravity);

        const auto& mu_f_arr = eff_visc[lev]->array(mfi);

        if (flagfab.getType(bx) != FabType::covered) {

          reduce_op.eval(bx, reduce_data,
            [ro,ep,gp0_dev,gradp,txfr,gravity_dev,vel,odx,ody,odz, flags,
             fluid_parms,idx_drag_coeff,idx_vel_src_c,mu_f_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
          {
            Real l_c_cfl = Real(-1.0);
            Real l_v_cfl = Real(-1.0);
            Real l_f_cfl = Real(-1.0);

            if (!flags(i,j,k).isCovered())
            {
                RealVect acc(0.);
                Real qro  = 1.0/ro(i,j,k);
                Real qep  = 1.0/ep(i,j,k);

                Real mu_g = mu_f_arr(i,j,k);

                // Compute the three components of the net acceleration
                // Explicit particle forcing is given by
                for (int idim(0); idim < 3; ++idim) {

                  Real delp = gp0_dev[idim] + gradp(i,j,k,idim);

                  Real fp   = txfr(i,j,k,idx_vel_src_c+idim)
                    - txfr(i,j,k,idx_drag_coeff) * vel(i,j,k,idim);

                  acc[idim] = gravity_dev[idim] + qro * ( - delp + fp*qep );
                }

                l_c_cfl = amrex::max(amrex::Math::abs(vel(i,j,k,0))*odx,
                                     amrex::Math::abs(vel(i,j,k,1))*ody,
                                     amrex::Math::abs(vel(i,j,k,2))*odz,
                                     l_c_cfl);

                l_v_cfl = amrex::max(2.0 * mu_g * qro * (odx*odx + ody*ody + odz*odz),
                                     l_v_cfl);

                l_f_cfl = amrex::max(amrex::Math::abs(acc[0])*odx,
                                     amrex::Math::abs(acc[1])*ody,
                                     amrex::Math::abs(acc[2])*odz,
                                     l_f_cfl);
            }
            return {l_c_cfl, l_v_cfl, l_f_cfl};
          });
        }
      }

      ReduceTuple host_tuple = reduce_data.value();
      c_cfl = amrex::max(c_cfl, amrex::get<0>(host_tuple));
      v_cfl = amrex::max(v_cfl, amrex::get<1>(host_tuple));
      f_cfl = amrex::max(f_cfl, amrex::get<2>(host_tuple));

    } //lev

    // Do global max operation
    ParallelDescriptor::ReduceRealMax({c_cfl, v_cfl, f_cfl});

    Real cv_cfl(c_cfl);
    if (explicit_diffusion) { cv_cfl += v_cfl; }

    // Max CFL factor for all levels
    Real cfl_max = cv_cfl + std::sqrt(cv_cfl*cv_cfl + Real(4.0) * f_cfl);

    // New dt
    Real eps = std::numeric_limits<Real>::epsilon();
    if (cfl_max > eps) {
       dt_new = m_cfl * 2.0 / cfl_max;
    } else {
       // This is totally random but just a way to set a timestep
       // when the initial velocity is zero and the forcing term
       // is not a body force
       dt_new = m_timer.dt_max();
    }

    // Optionally reduce CFL for initial step
    if ( m_timer.nstep() < 1) {
        dt_new *= m_scale_init_dt;
    }

    // Protect against cfl_max very small
    // This may happen, for example, when the initial velocity field
    // is zero for an inviscid flow with no external forcing
    //
    // Note, the following two lines of code are problematic for zero-flow
    // simulations used for testing
//    Real eps = std::numeric_limits<Real>::epsilon();
//    if ( nstep > 1 && cfl_max <= eps ) dt_new = 0.5 * prev_dt;

    // Don't let the timestep grow by more than 10% per step.
    if ( m_timer.nstep() > 0 ) { dt_new = amrex::min( dt_new, 1.1*m_timer.prev_dt() ); }

    // dt_new is the step calculated with a cfl constraint; dt is the value set by fixed_dt
    // When the test was on dt > dt_new, there were cases where they were effectively equal
    //   but (dt > dt_new) was being set to true due to precision issues.
    Real ope(1.0 + 1.e-8);

    if (m_timer.timestep_type() == MFIXTimer::TimestepType::Fixed) {

      if (m_timer.dt() > dt_new*ope && m_cfl > 0) {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Fixed dt is too large to satisfy CFL condition:\n"
          << "   fixed dt =       "  << m_timer.dt() << "\n"
          << "   dt based on cfl: " << dt_new;
      }

    } else {

      m_timer.dt() = amrex::min( dt_new, m_timer.dt_max() );

      if ( m_timer.dt() < m_timer.dt_min() ) {
        reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Current dt is smaller than dt_min";
      }
    }

    // Don't overshoot the final time if not running to steady state
    if ( !m_timer.SteadyState() && !m_timer.overstep_end_time()) {
      if ( (m_timer.new_time() - m_timer.stop_time()) > 1.e-15) {
        m_timer.dt() = m_timer.stop_time() - m_timer.time();
      }
    }

}
