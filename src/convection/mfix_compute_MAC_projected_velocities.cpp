#include <hydro_MacProjector.H>
#include <hydro_utils.H>

#include <mfix.H>
#include <mfix_bc.H>
#include <mfix_eb_parms.H>

/***************************************************************************
 *  Compute the MAC projected velocities:  (εu)_mac                        *
 *                                                                         *
 *  1) Extrapolate cell-centered velocites to cell faces. u_face           *
 *  2) Average volume fraction and densty to cell faces.  ε_face, ρ_face   *
 *  3) Compute 'beta' coefficient. b_coeff = (ε_face / ρ_face)             *
 *  4) Compute (εu)_face = (ε_face * u_face)                               *
 *  5) Apply MAC projection to get (εu)_mac                                *
 *                                                                         *
 *             -∇·(ε/ρ)∇φ = ∇·(εu)_mac - ∇·(εu)_face                       *
 *                                                                         *
 *   ε := volume fraction      ( )_face denotes face variables             *
 *   ρ := density              ( )_mac  denotes MAC velocities             *
 *   u := velocity                                                         *
 *   φ := 'mac phi'                                                        *
 *                                                                         *
 *  The AMReX MLABecLaplacian linear solver class uses the canonical form: *
 *                                                                         *
 *  Aα -B∇·β∇φ = RHS                                                       *
 *                                                                         *
 *  Aα = 0;   B = 1;   β := ε/ρ; and RHS = ∇·(εu)_mac - ∇·(εu)_face        *
 *                                                                         *
 *  The fluid constraint defines ∇·(εu)_mac := a_mac_rhs                   *
 *                                                                         *
 *  AMReX Hydro MACProjector class calculates -∇·(εu)_face and combines    *
 *  with a_mac_rhs to complete the RHS.                                    *
 *                                                                         *
 ***************************************************************************/
void mfix::
compute_MAC_projected_velocities ( Real const /*a_time*/, Real const a_dt,
                                   Vector< MultiFab const*> const& a_vel,
                                   Vector< MultiFab const*> const& a_mac_rhs,
                                   Vector< Array<MultiFab      *,3>> const& a_mac_vel,
                                   Vector< MultiFab      *> const& a_mac_phi,
                                   Vector< MultiFab const*> const& a_epf,
                                   Vector< MultiFab const*> const& a_rho,
                                   Vector< MultiFab const*> const& a_txfr,
                                   Vector< MultiFab const*> const& a_eb_vel,
                                   Vector< MultiFab      *> const& a_vel_forces,
                                   Vector< MultiFab const*> const& a_divtau)
{
  BL_PROFILE("mfix::compute_MAC_projected_velocities()");

  if (m_verbose) { Print() << "MAC Projection:\n"; }

  if (!macproj) {
    macproj = std::make_unique<Hydro::MacProjector>( Geom(0,nlev()-1),
        MLMG::Location::FaceCentroid,  // Location of mac_vec
        MLMG::Location::FaceCentroid,  // Location of beta
        MLMG::Location::CellCenter,    // Location of solution variable phi
        MLMG::Location::CellCentroid); // Location of MAC RHS
  }

  // We first compute the velocity forcing terms to be used in predicting
  //    to faces before the MAC projection
  if (advection_type() != AdvectionType::MOL) {

    bool include_pressure_gradient = !(m_godunov_use_mac_phi);

    compute_vel_forces(a_vel_forces, a_vel, a_rho, a_txfr,
        include_pressure_gradient);

    if (m_godunov_include_diff_in_forcing) {

      for (int lev(0); lev < nlev(); ++lev) {
        for (MFIter mfi(*(leveldata().div_tau(lev)),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          Box const& bx = mfi.tilebox();

          Array4<Real      > const& forces = a_vel_forces[lev]->array(mfi);
          Array4<Real const> const& rho    = a_rho[lev]->const_array(mfi);
          Array4<Real const> const& divtau = a_divtau[lev]->const_array(mfi);

          ParallelFor(bx, [forces, divtau, rho]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            forces(i,j,k,0) += divtau(i,j,k,0)/rho(i,j,k);
            forces(i,j,k,1) += divtau(i,j,k,1)/rho(i,j,k);
            forces(i,j,k,2) += divtau(i,j,k,2)/rho(i,j,k);
          });
        }
      }
    }

    if (nghost_force() > 0) {
      for (int lev(0); lev < nlev(); ++lev) {
        bcs().fillpatch(lev, 0., BCFillVar::none, a_vel_forces, nghost_force());
      }
    }
  } // !MOL

  // Predict normal velocity to faces -- note that the {u_mac, v_mac, w_mac}
  // arrays returned from this call are on face CENTROIDS.
  for (int lev(0); lev < nlev(); ++lev) {

    const EBFArrayBoxFactory* ebfact = &m_eb->Factory(lev);

    // We need this to avoid FPE
    a_mac_vel[lev][0]->setVal(covered_val);
    a_mac_vel[lev][1]->setVal(covered_val);
    a_mac_vel[lev][2]->setVal(covered_val);

    std::string advection_string =
        (advection_type() == AdvectionType::Godunov) ? "Godunov" : "MOL";

    if (m_embedded_boundaries.has_flow()) {

      HydroUtils::ExtrapVelToFaces(*a_vel[lev], *a_vel_forces[lev],
          *a_mac_vel[lev][0], *a_mac_vel[lev][1], *a_mac_vel[lev][2],
          bcs().get_hydro_velocity_bcrec(), bcs().get_hydro_velocity_bcrec_device_ptr(),
          geom[lev], a_dt, *ebfact, a_eb_vel[lev],
          m_godunov_ppm, m_godunov_use_forces_in_trans,
          advection_string);

    } else {

       HydroUtils::ExtrapVelToFaces(*a_vel[lev], *a_vel_forces[lev],
          *a_mac_vel[lev][0], *a_mac_vel[lev][1], *a_mac_vel[lev][2],
          bcs().get_hydro_velocity_bcrec(), bcs().get_hydro_velocity_bcrec_device_ptr(),
          geom[lev], a_dt, *ebfact,
          m_godunov_ppm, m_godunov_use_forces_in_trans,
          advection_string);
    }
  }

  // Face-based coefficients b in MAC projection and implicit diffusion solve
  Vector< Array< std::unique_ptr<MultiFab>,3> > bcoeff( nlev() );

  // ro_face and ep_face are temporary, no need to keep it outside this routine
  Vector< Array< std::unique_ptr<MultiFab>,3> > ro_face( nlev() );
  Vector< Array< std::unique_ptr<MultiFab>,3> > ep_face( nlev() );


  for (int lev(0); lev<nlev(); ++lev) {

    EBFArrayBoxFactory const& factory = m_eb->Factory(lev);

    const BoxArray uba = a_mac_vel[lev][0]->boxArray();
    const BoxArray vba = a_mac_vel[lev][1]->boxArray();
    const BoxArray wba = a_mac_vel[lev][2]->boxArray();

    bcoeff[lev] = {AMREX_D_DECL(
        std::make_unique<MultiFab>(uba,dmap[lev],1,nghost_state(),MFInfo(),factory),
        std::make_unique<MultiFab>(vba,dmap[lev],1,nghost_state(),MFInfo(),factory),
        std::make_unique<MultiFab>(wba,dmap[lev],1,nghost_state(),MFInfo(),factory)
    )};

    ep_face[lev] = {AMREX_D_DECL(
        std::make_unique<MultiFab>(uba,dmap[lev],1,0,MFInfo(),factory),
        std::make_unique<MultiFab>(vba,dmap[lev],1,0,MFInfo(),factory),
        std::make_unique<MultiFab>(wba,dmap[lev],1,0,MFInfo(),factory)
    )};

    ro_face[lev] = {AMREX_D_DECL(
        std::make_unique<MultiFab>(uba,dmap[lev],1,0,MFInfo(),factory),
        std::make_unique<MultiFab>(vba,dmap[lev],1,0,MFInfo(),factory),
        std::make_unique<MultiFab>(wba,dmap[lev],1,0,MFInfo(),factory)
    )};

    for (int idim(0); idim < AMREX_SPACEDIM; ++idim) {
      bcoeff[lev][idim]->setVal(0.);
      ep_face[lev][idim]->setVal(covered_val);
    }

    // Define ep and rho on face centroids (using interpolation from cell centroids)
    // The only use of bcs in this call is to test on whether a domain boundary is ext_dir
    EB_interp_CellCentroid_to_FaceCentroid (*a_rho[lev], GetArrOfPtrs(ro_face[lev]),
        /*scomp=*/0, /*dcomp=*/0, /*ncomp=*/1, geom[lev], bcs().get_density_bcrec());

    EB_interp_CellCentroid_to_FaceCentroid (*a_epf[lev], GetArrOfPtrs(ep_face[lev]),
        /*scomp=*/0, /*dcomp=*/0, /*ncomp=*/1, geom[lev], bcs().get_volfrac_bcrec());

    // Compute beta coefficients for div(beta*grad(phi)) = RHS:  beta = ep / ro
    for (int idim(0); idim < AMREX_SPACEDIM; ++idim) {
      MultiFab::Copy(*bcoeff[lev][idim], *(ep_face[lev][idim]), 0, 0, 1, 0);
      MultiFab::Divide(*bcoeff[lev][idim], *(ro_face[lev][idim]), 0, 0, 1, 0);
    }

  }


  for (int lev(0); lev < nlev(); ++lev) {

    for (MFIter mfi(*a_vel[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& ubx = mfi.nodaltilebox(0);
      const Box& vbx = mfi.nodaltilebox(1);
      const Box& wbx = mfi.nodaltilebox(2);

      // Face-centered velocity components
      Array4<Real      > const& umac = a_mac_vel[lev][0]->array(mfi);
      Array4<Real      > const& vmac = a_mac_vel[lev][1]->array(mfi);
      Array4<Real      > const& wmac = a_mac_vel[lev][2]->array(mfi);

      // Face-centroid volume fractions
      Array4<Real const> const& epx = ep_face[lev][0]->const_array(mfi);
      Array4<Real const> const& epy = ep_face[lev][1]->const_array(mfi);
      Array4<Real const> const& epz = ep_face[lev][2]->const_array(mfi);

      // Now we multiply the face velocities by the phasic volume fraction
      // so we have {ep * u_mac, ep * v_mac, ep * w_mac}.
      amrex::ParallelFor(ubx, vbx, wbx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { umac(i,j,k) *= epx(i,j,k); },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { vmac(i,j,k) *= epy(i,j,k); },
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { wmac(i,j,k) *= epz(i,j,k); });
    }
  }


  //
  // Initialize (or redefine the beta in) the MacProjector
  //

  if (macproj->needInitialization()) {

    LPInfo lp_info;
    lp_info.setMaxCoarseningLevel(macproj_options->max_coarsening_level);
    lp_info.setAgglomerationGridSize(agg_grid_size);

    macproj->initProjector(lp_info, GetVecOfArrOfConstPtrs(bcoeff) );

    macproj_options->apply(*macproj);

    macproj->setDomainBC(m_boundary_conditions.mac_lobc(),
                         m_boundary_conditions.mac_hibc());

  } else {

    macproj->updateBeta( GetVecOfArrOfConstPtrs(bcoeff) );

  }

  if (m_verbose) { Print() << " >> Before projection\n" ; }

  Vector< std::unique_ptr<MultiFab> > div_epu( nlev() );

  for ( int lev(0); lev < nlev(); ++lev ) {

    // Set bcs on (ep * u_mac)
    set_MAC_velocity_bcs(lev, a_mac_rhs[lev], a_mac_vel[lev]);

    if (m_verbose) {

      div_epu[lev] = std::make_unique<MultiFab>(grids[lev],dmap[lev],
          /*ncomp=*/1, /*nghost=*/1, MFInfo(), m_eb->Factory(lev));

      EB_computeDivergence(*div_epu[lev], GetArrOfConstPtrs(a_mac_vel[lev]),
          geom[lev], /*already_on_centroid=*/true);

      Print() << "  * On level "<< lev << " max(abs(diveu)) = "
              << div_epu[lev]->norm0(0,0,false,true) << "\n";
    }
  }

  macproj->setUMAC(a_mac_vel);
  macproj->setDivU(a_mac_rhs);

  if (m_embedded_boundaries.has_flow()) {
    for (int lev(0); lev < nlev(); ++lev) {
      macproj->setEBInflowVelocity(lev, *a_eb_vel[lev]);
    }
  }

  Real const atol = macproj_options->mg_atol;
  Real const rtol = macproj_options->mg_rtol;

  if (m_timer.SteadyState()) {

    // Solve using mac_phi as an initial guess -- note that mac_phi is
    //       stored from iteration to iteration
    macproj->project(a_mac_phi, rtol, atol);

  } else {

    if (m_godunov_use_mac_phi) {

      for (int lev(0); lev < nlev(); ++lev) {
        a_mac_phi[lev]->mult(a_dt/2.,0,1,1);
      }

      macproj->project(a_mac_phi, rtol, atol);

      for (int lev=0; lev < nlev(); ++lev) {
        a_mac_phi[lev]->mult(2./a_dt,0,1,1);
      }

    } else {

      for (int lev(0); lev < nlev(); ++lev) {
        a_mac_phi[lev]->setVal(0.);
      }

      // Solve with initial guess of zero
      macproj->project(rtol, atol);
    }
  }

  if (m_verbose) { Print() << " >> After projection\n" ; }

  for ( int lev(0); lev < nlev() ; ++lev ) {

    a_mac_vel[lev][0]->FillBoundary(geom[lev].periodicity());
    a_mac_vel[lev][1]->FillBoundary(geom[lev].periodicity());
    a_mac_vel[lev][2]->FillBoundary(geom[lev].periodicity());

    if (m_verbose) {

      EB_computeDivergence(*div_epu[lev], GetArrOfConstPtrs(a_mac_vel[lev]),
          geom[lev], /*already_on_centroid=*/true);

      Print() << "  * On level "<< lev << " max(abs(diveu)) = "
          << div_epu[lev]->norm0(0,0,false,true) << "\n";
    }
  }

}
