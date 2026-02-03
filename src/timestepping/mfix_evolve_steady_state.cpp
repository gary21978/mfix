#include <AMReX_BCUtil.H>

#include <mfix.H>
#include <mfix_utils.H>

// This subroutine is the driver for the whole time stepping (fluid + particles )
void mfix::
EvolveSteadyState ()
{
  BL_PROFILE("mfix::EvolveSteadyState()");


  ///////////////////////////////////////////////////////////////////////////////
  //                           PREPARE FLUID FOR STEP                          //
  ///////////////////////////////////////////////////////////////////////////////
  MFIXStepData stepData(nlev(), grids, dmap, eb()->factory_const(),
      nghost_state(), nghost_force(), m_boundary_conditions);

  prepareToStep( stepData );

  stepData.include_pressure();
  stepData.set_pressure( leveldata().pert_p_const() );


  ///////////////////////////////////////////////////////////////////////////////
  //                            FLUID PREDICTOR STEP                           //
  ///////////////////////////////////////////////////////////////////////////////

  int iter (0);
  do {

    // Update iteration count
    ++iter;

    compute_dt(stepData.avg_particle_data_const());

    Print() << "\n   Iteration " << iter << " with dt = " << m_timer.dt() << "\n";

    // NOTE: It is important to do fillpatch on the NEW variables,
    // before swapping pointers.
    for (int lev(0); lev < nlev(); lev++) {

      Real const time = timer().time();
      int const ng = nghost_state();

      bcs().fillpatch(lev, time, BCFillVar::vel, leveldata().vel(), ng);
      bcs().fillpatch(lev, time, BCFillVar::rho, leveldata().rho(), ng);

      if (leveldata(lev)->has_enthalpy()) {
        bcs().fillpatch(lev, time, BCFillVar::h, leveldata().h(), ng);
      }
      if (leveldata(lev)->has_temperature()) {
        bcs().fillpatch(lev, time, BCFillVar::T, leveldata().T(), ng);
      }
      if (leveldata(lev)->has_species()) {
        bcs().fillpatch(lev, time, BCFillVar::X, leveldata().X(), ng);
      }
      if (leveldata(lev)->has_tracer()) {
        bcs().fillpatch(lev, time, BCFillVar::tracer, leveldata().tracer(), ng);
      }

      leveldata(lev)->swapOldAndNew();
    }
    std::swap(therm_p, therm_po);

    // Save a copy of the current perturbational pressure so we can compute
    // the residual for steady-state runs.
    stepData.set_pressure( leveldata().pert_p_const() );

    // Predictor step
    apply_predictor(stepData, m_timer.time(), m_timer.dt(),
        m_timer.prev_dt(), /*proj_2_pred=*/true);

  } while ( !SteadyStateReached(m_timer.dt(), iter,
             stepData.pert_pressure_const()));

}


//
// Check if steady state has been reached by verifying that
//
//      max(abs( u^(n+1) - u^(n) )) < tol * dt
//      max(abs( v^(n+1) - v^(n) )) < tol * dt
//      max(abs( w^(n+1) - w^(n) )) < tol * dt
//

int mfix::
SteadyStateReached ( Real dt, int iter,
                     Vector< MultiFab const*> const& a_p_go)
{

  Vector<int> condition1(nlev());
  Vector<int> condition2(nlev());

  Real time = 0.;

  //
  // Make sure velocity is up to date
  //
  for (int lev(0); lev<nlev(); ++lev) {

    bcs().set_velocity_bcs(lev, time, leveldata().vel(lev));

   //
   // Use temporaries to store the difference
   // between current and previous solution
   //
   MultiFab temp_vel(leveldata().vel(lev)->boxArray(),
                     leveldata().vel(lev)->DistributionMap(), 3, 0);

   MultiFab::LinComb(temp_vel, 1.0, *(leveldata().vel(lev)), 0, -1.0,
                     *(leveldata().vel_old(lev)), 0, 0, 3, 0);

   MultiFab tmp;

   const BoxArray & nd_grid = amrex::convert(grids[lev], IntVect{1,1,1});
   tmp.define(nd_grid, dmap[lev], 1, 0);

   MultiFab::LinComb(tmp, 1.0, *(leveldata().pert_p(lev)), 0, -1.0, *a_p_go[lev], 0, 0, 1, 0);

   Real delta_u = temp_vel.norm0(0,0,false,true);
   Real delta_v = temp_vel.norm0(1,0,false,true);
   Real delta_w = temp_vel.norm0(2,0,false,true);
   Real delta_p = temp_vel.norm0(0,0,false,true);

   Real tol = m_timer.SteadyStateTol();

   condition1[lev] = (delta_u < tol*dt) && (delta_v < tol*dt ) && (delta_w < tol*dt);

   //
   // Second stop condition
   //
   Periodicity period = geom[lev].periodicity();

   Real du_n1 = temp_vel.norm1(0,period,true);
   Real dv_n1 = temp_vel.norm1(1,period,true);
   Real dw_n1 = temp_vel.norm1(2,period,true);
   Real dp_n1 =      tmp.norm1(0,period,true);
   Real uo_n1 = leveldata().vel_old(lev)->norm1(0,period,true);
   Real vo_n1 = leveldata().vel_old(lev)->norm1(1,period,true);
   Real wo_n1 = leveldata().vel_old(lev)->norm1(2,period,true);
   Real po_n1 = a_p_go[lev]->norm1(0,period,true);

   Real tmp1, tmp2, tmp3, tmp4;

   Real local_tol = 1.0e-8;

   tmp1 = ( uo_n1 < local_tol ) ? 0.0 : du_n1 / uo_n1;
   tmp2 = ( vo_n1 < local_tol ) ? 0.0 : dv_n1 / vo_n1;
   tmp3 = ( wo_n1 < local_tol ) ? 0.0 : dw_n1 / wo_n1;
   tmp4 = ( po_n1 < local_tol ) ? 0.0 : dp_n1 / po_n1;

   condition2[lev] = (tmp1 < tol) && (tmp2 < tol) && (tmp3 < tol); // && (tmp4 < tol);

   //
   // Print out info on steady state checks
   //
   Print() << "\nSteady state check at level " << lev << ":"
     << "\n||u-uo||/||uo|| = " << std::setw(12) << std::scientific << std::setprecision(4) << tmp1
     << "      du/dt = " << std::setw(12) << std::scientific << std::setprecision(4) << delta_u/dt
     << "\n||v-vo||/||vo|| = " << std::setw(12) << std::scientific << std::setprecision(4) << tmp2
     << "      dv/dt = " << std::setw(12) << std::scientific << std::setprecision(4) << delta_v/dt
     << "\n||w-wo||/||wo|| = " << std::setw(12) << std::scientific << std::setprecision(4) << tmp3
     << "      dw/dt = " << std::setw(12) << std::scientific << std::setprecision(4) << delta_w/dt
     << "\n||p-po||/||po|| = " << std::setw(12) << std::scientific << std::setprecision(4) << tmp4
     << "      dp/dt = " << std::setw(12) << std::scientific << std::setprecision(4) << delta_p/dt
       << "\n\n";
  }

  int reached(1);

  for (int lev(0); lev<nlev(); ++lev) {
   reached = reached && (condition1[lev] || condition2[lev]);
  }

  reached = reached || (iter >= m_timer.SteadyStateMaxIter());

  // Always return fail for first iteration to prevent a "zero-velocity"
  // field initial condition to signal a false-positive.
  return (( iter == 1 ) ? 0 : reached);
}
