#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>

// For multigrid
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLNodeLaplacian.H>

#include <mfix.H>
#include <mfix_bc.H>
#include <mfix_mf_helpers.H>
#include <mfix_utils.H>
#include <mfix_monitors.H>
#include <mfix_run_on.H>

/***************************************************************************
 * Compute the right-hand side (RHS) for the MAC projection:               *
 *                                                                         *
 *             -∇·(ε/ρ)∇φ = ∇·(εu)_mac - ∇·(εu)_edge                       *
 *                         |                        |                      *
 *                         |<--- projection RHS --->|                      *
 *                                                                         *
 *   ε := volume fraction      ( )_edge denotes edge variables             *
 *   ρ := density              ( )_mac  denotes MAC variables              *
 *   u := velocity                                                         *
 *   φ := unknown                                                          *
 *                                                                         *
 * The AMReX MLABecLaplacian linear solver class uses the canonical form:  *
 *                                                                         *
 *  Aα -B∇·β∇φ = RHS                                                       *
 *                                                                         *
 *  Aα = 0;   B = 1;   β := ε/ρ; and RHS = ∇·(εu)_mac - ∇·(εu)_edge        *
 *                                                                         *
 *  This routine computes ∇·(εu)_mac while AMReX Hydro MACProjector class  *
 *  calculates -∇·(εu)_edge to complete the RHS.                           *
 *                                                                         *
 *  The fluid constraint defines ∇·(εu)_mac                                *
 *                                                                         *
 *                                                                         *
 * Incompressible fluid .................................................. *
 *                                                                         *
 *   ∇·(εu)_mac = - ∂ε/∂t + G/ρ                                            *
 *                                                                         *
 *   G := Net mass production (consumption) of fluid from reactions        *
 *                                                                         *
 *                                                                         *
 * Ideal gas open-system ................................................. *
 *                                                                         *
 *   ∇·(εu)_mac = - ∂ε/∂t + S                                              *
 *                                                                         *
 *   where  S = Sh / (ρ*cp*T) + (1/ρ)Σ_k[W/Wk - (hk/(c*T))*SXk]            *
 *                                                                         *
 *   Sh  := RHS of enthalpy equation                                       *
 *   SXk := RHS of k-th species equation                                   *
 *                                                                         *
 *   cp  := mixture specific heat                                          *
 *   W   := mixture molecular weight                                       *
 *   T   := Temperature                                                    *
 *   Wk  := k-th species molecular weight                                  *
 *   hk  := k-th species specific enthalpy                                 *
 *                                                                         *
 *                                                                         *
 * Ideal gas closed-system ............................................... *
 *                                                                         *
 *   ∇·(εu)_mac = - ∂ε/∂t + δS - δθ(<S>/<θ>)                               *
 *                                                                         *
 *   where S is the ideal gas open-system constraint RHS and               *
 *                                                                         *
 *   θ = ε/(Γ.pT) = (ε/pT)*(1 - 1/(cp*ρ*T))                                *
 *                                                                         *
 *   S = <S> + δS      <x> := average of x over domain (scalar)            *
 *   θ = <θ> + δθ      δx  := fluctuation of x (field)                     *
 *                                                                         *
 ***************************************************************************/
void mfix::
compute_MAC_proj_RHS( Vector< MultiFab      * > const& a_rhs,
                      Vector< MultiFab const* > const& a_depdt,
                      Vector< MultiFab const* > const& a_epf,
                      Vector< MultiFab const* > const& a_rho,
                      Vector< MultiFab const* > const& a_Xf,
                      Vector< MultiFab const* > const& a_Tf,
                      Vector< MultiFab const* > const& a_lap_Tf,
                      Vector< MultiFab const* > const& a_div_hJ,
                      Vector< MultiFab const* > const& a_div_J,
                      Vector< MultiFab      * > const& a_txfr,
                      Vector< MultiFab const* > const& a_chem,
                      Real const& a_therm_p,
                      Real      & a_RHS_therm_p)
{
  InterphaseChemTxfrIndexes chem_idxs(fluid.nspecies(), reactions.solve());

  for (int lev(0); lev < nlev(); ++lev) {
    a_rhs[lev]->setVal(0.);
  }

  int const ngrow(0);

  Vector <std::unique_ptr<MultiFab>> energy_rhs(nlev());
  Vector <std::unique_ptr<MultiFab>> species_rhs(nlev());

  if ( fluid.constraint.isIdealGas() ) {

    for (int lev(0); lev < nlev(); ++lev) {

      const BoxArray&            ba( a_rhs[lev]->boxArray() );
      const DistributionMapping& dm( a_rhs[lev]->DistributionMap() );

      // Compute the RHS of the energy equation:
      // ------------------------------------------------------------------------------//
      if (fluid.solve_enthalpy()) {

        energy_rhs[lev].reset( new MultiFab( ba, dm, 1, ngrow,
            MFInfo(), a_rhs[lev]->Factory() ));

        energy_rhs[lev]->setVal(0.);

        // Thermal diffusion (condution)
        MultiFab::Add(*energy_rhs[lev], *a_lap_Tf[lev], 0, 0, 1, ngrow);

        // Species enthalpy flux -- fluxes are computed with a negative
        // sign, therefore substrat to add.
        if (fluid.solve_species()) {
          MultiFab::Subtract(*energy_rhs[lev], *a_div_hJ[lev], 0, 0, 1, ngrow);
        }

        // Convective heat transfer: coeff*(Tp - Tf) = coeff*Tp - coeff*Tp
        { InterphaseTxfrIndexes txfr_idxs;

          const int idx_convection_coeff_txfr = txfr_idxs.convection_coeff;
          const int idx_energy_source_txfr    = txfr_idxs.energy_source;

          // Add in the explicit component, coeff*Tp
          MultiFab::Add(*energy_rhs[lev], *a_txfr[lev], idx_energy_source_txfr, 0, 1, ngrow);

          // Flip the sign of the heat transfer coefficient
          a_txfr[lev]->mult(-1.0, idx_convection_coeff_txfr, 1, ngrow);

          // Subtract the coefficient multiplied by fluid temperature += -1*coeff*Tf
          MultiFab::AddProduct(*energy_rhs[lev], *a_txfr[lev], idx_convection_coeff_txfr,
                *a_Tf[lev], 0, 0, 1, ngrow);

          // Flip the sign of the heat transfer coefficient back
          a_txfr[lev]->mult(-1.0, idx_convection_coeff_txfr, 1, ngrow);
        }

        // Enthalpy transfer from mass transfer
        if (reactions.solve()) {

          int const idx_chem_source(chem_idxs.chem_h);

          MultiFab::Add(*energy_rhs[lev], *a_chem[lev], idx_chem_source, 0, 1, ngrow);
        }
      } // energy_rhs

      // Compute the RHS of the species equations:
      // ------------------------------------------------------------------------------//
      if (fluid.solve_species()) {

        int const nspecies( fluid.nspecies() );

        species_rhs[lev].reset( new MultiFab( ba, dm, nspecies, ngrow,
            MFInfo(), a_rhs[lev]->Factory() ));

        species_rhs[lev]->setVal(0.);

        // Species diffusion
        MultiFab::Subtract(*species_rhs[lev], *a_div_J[lev], 0, 0, nspecies, ngrow);

        if (reactions.solve()) {

          int const idx_chem_source(chem_idxs.chem_ro_gk);

          MultiFab::Add(*species_rhs[lev], *a_chem[lev], 0, idx_chem_source, nspecies, ngrow);

        }
      } // species_rhs
    } // level

    // Compute: S = Sh / (ρ*cp*T) + (1/ρ)Σ_k[W/Wk - (hk/(c*T))*SXk]
    idealgas_opensystem_rhs(a_rhs, GetVecOfConstPtrs(energy_rhs),
        GetVecOfConstPtrs(species_rhs), a_rho, a_Tf, a_Xf);

    // Compute: δS - δθ(<S>/<θ>) where S is computed from open-system
    if ( fluid.solve_enthalpy() && fluid.constraint.isIdealGasClosedSystem() ) {

      idealgas_closedsystem_rhs(a_rhs, a_epf, a_rho, a_Tf, a_Xf,
          a_therm_p, a_RHS_therm_p);
    }
  }

  for (int lev(0); lev < nlev(); lev++) {

    //- ∂ε/∂t = -(epf_new - epf_old) / dt
    MultiFab::Subtract( *a_rhs[lev], *a_depdt[lev], 0, 0, 1, ngrow);

    EB_set_covered(*a_rhs[lev], 0, 1, 0, 0.0);

    a_rhs[lev]->FillBoundary(geom[lev].periodicity());

#if 0
    Print() << "MAC proj on level " << lev
        << "  min(rhs) = " << a_rhs[lev]->min(0)
        << "  max(rhs) = " << a_rhs[lev]->max(0) << '\n';

    if (a_rhs[0]->nGrow() > 0) {
      bcs().fillpatch(lev, 0., BCFillVar::none, a_rhs);
    }
#endif
  }
}


/***************************************************************************
 * Compute the right-hand side (RHS) for the Nodal projection:             *
 *                                                                         *
 *                ∇·(ε/ρ)∇φ =  ∇·(εu^**) - ∇·(εu)                          *
 *                           |                    |                        *
 *                           |<- projection RHS ->|                        *
 *                                                                         *
 *   ε := volume fraction      ( )^** denotes explicit update              *
 *   ρ := density                                                          *
 *   u := velocity                                                         *
 *   φ := unknown                                                          *
 *                                                                         *
 * The AMReX MLNodeLaplacian linear solver class uses the canonical form:  *
 *                                                                         *
 *  ∇·σ∇φ = RHS                                                            *
 *                                                                         *
 *  σ := ε/ρ   and   RHS = -(∇·(εu) - ∇·(εu)^**)                           *
 *                         ^                                               *
 *                         |                                               *
 *                                                                         *
 *     --->!!! RHS HAS OPPOSITE SIGN FROM FROM MAC PROJECTION !!!<---      *
 *                                                                         *
 *                                                                         *
 *  This routine computes -∇·(εu) while AMReX calculates ∇·(εu)_edge       *
 *                                                                         *
 ***************************************************************************/
void mfix::
compute_nodal_proj_RHS( Vector< MultiFab      * > const& a_rhs,
                        Vector< MultiFab const* > const& a_depdt,
                        Vector< MultiFab const* > const& a_epf_o,
                        Vector< MultiFab const* > const& a_epf_n,
                        Vector< MultiFab const* > const& a_rho_o,
                        Vector< MultiFab const* > const& a_rho_n,
                        Vector< MultiFab const* > const& a_hf_o,
                        Vector< MultiFab const* > const& a_hf_n,
                        Vector< MultiFab const* > const& a_Tf_n,
                        Vector< MultiFab const* > const& a_dhdt,
                        Vector< MultiFab const* > const& a_Xf_o,
                        Vector< MultiFab const* > const& a_Xf_n,
                        Vector< MultiFab const* > const& a_dXdt,
                        Real const& a_therm_p,
                        Real      & a_RHS_therm_p,
                        Real const& a_dt)
{
  InterphaseChemTxfrIndexes chem_idxs(fluid.nspecies(), reactions.solve());

  for (int lev(0); lev < nlev(); ++lev) {
    a_rhs[lev]->setVal(0.);
  }

  int const ngrow(0);

  Vector <std::unique_ptr<MultiFab>> energy_rhs(nlev());
  Vector <std::unique_ptr<MultiFab>> species_rhs(nlev());

  if ( fluid.constraint.isIdealGas() ) {

    Real const inv_dt( 1_rt / a_dt );

    for (int lev(0); lev < nlev(); ++lev) {

      const BoxArray&            ba( a_rhs[lev]->boxArray() );
      const DistributionMapping& dm( a_rhs[lev]->DistributionMap() );

      // Compute the RHS of the energy equation:
      // ------------------------------------------------------------------------------//
      if (fluid.solve_enthalpy()) {

        energy_rhs[lev].reset( new MultiFab( ba, dm, 1, ngrow,
            MFInfo(), a_rhs[lev]->Factory() ));

        energy_rhs[lev]->setVal(0.);

        for (MFIter mfi(*a_rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          Box const& bx = mfi.tilebox();

          Array4<Real      > const& h_rhs  = energy_rhs[lev]->array(mfi);

          Array4<Real const> const& hf_n   = a_hf_n[lev]->const_array(mfi);
          Array4<Real const> const& hf_o   = a_hf_o[lev]->const_array(mfi);
          Array4<Real const> const& rho_n  = a_rho_n[lev]->const_array(mfi);
          Array4<Real const> const& rho_o  = a_rho_o[lev]->const_array(mfi);
          Array4<Real const> const& epf_n  = a_epf_n[lev]->const_array(mfi);
          Array4<Real const> const& epf_o  = a_epf_o[lev]->const_array(mfi);

          int const h_comp(AdvectionComp::enthalpy);
          Array4<Real const> const& dhdt = a_dhdt[lev]->const_array(mfi,h_comp);

          Real const dpdt( a_RHS_therm_p );

          ParallelFor(bx, [h_rhs, epf_n, rho_n, hf_n, epf_o, rho_o, hf_o, inv_dt, dhdt, dpdt]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            h_rhs(i,j,k) = (epf_n(i,j,k)*rho_n(i,j,k)*hf_n(i,j,k)
                          - epf_o(i,j,k)*rho_o(i,j,k)*hf_o(i,j,k)) * inv_dt;

            h_rhs(i,j,k) -= dhdt(i,j,k);

            Real const avg_epf = 0.5*(epf_o(i,j,k) + epf_n(i,j,k));
            h_rhs(i,j,k) -= avg_epf*dpdt;
          });
        }

      } // energy_rhs

      // Compute the RHS of the species equations:
      // ------------------------------------------------------------------------------//
      if (fluid.solve_species()) {

        int const nspecies( fluid.nspecies() );

        species_rhs[lev].reset( new MultiFab( ba, dm, nspecies, ngrow,
            MFInfo(), a_rhs[lev]->Factory() ));

        species_rhs[lev]->setVal(0.);

        for (MFIter mfi(*a_rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          Box const& bx = mfi.tilebox();

          Array4<Real      > const& X_rhs  = species_rhs[lev]->array(mfi);

          Array4<Real const> const& Xf_n  = a_Xf_n[lev]->const_array(mfi);
          Array4<Real const> const& Xf_o  = a_Xf_o[lev]->const_array(mfi);
          Array4<Real const> const& rho_n = a_rho_n[lev]->const_array(mfi);
          Array4<Real const> const& rho_o = a_rho_o[lev]->const_array(mfi);
          Array4<Real const> const& epf_n = a_epf_n[lev]->const_array(mfi);
          Array4<Real const> const& epf_o = a_epf_o[lev]->const_array(mfi);

          Array4<Real const> const& dXdt  = a_dXdt[lev]->const_array(mfi);

          ParallelFor(bx, [X_rhs, epf_n, rho_n, Xf_n, epf_o, rho_o, Xf_o, inv_dt, dXdt, nspecies]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {

            for (int n(0); n < nspecies; ++n) {

              X_rhs(i,j,k,n) = ( epf_n(i,j,k) * rho_n(i,j,k) * Xf_n(i,j,k,n)
                              -  epf_o(i,j,k) * rho_o(i,j,k) * Xf_o(i,j,k,n) ) * inv_dt;

              X_rhs(i,j,k,n) -= dXdt(i,j,k,n);
            }
          });
        }
      } // species_rhs
    } // level


    // Compute: S = Sh / (ρ*cp*T) + (1/ρ)Σ_k[W/Wk - (hk/(c*T))*SXk]
    idealgas_opensystem_rhs(a_rhs, GetVecOfConstPtrs(energy_rhs),
        GetVecOfConstPtrs(species_rhs), a_rho_n, a_Tf_n, a_Xf_n);

    // Compute: δS - δθ(<S>/<θ>) where S is computed from open-system
    if ( fluid.solve_enthalpy() && fluid.constraint.isIdealGasClosedSystem() ) {

      idealgas_closedsystem_rhs(a_rhs, a_epf_n, a_rho_n, a_Tf_n, a_Xf_n,
          a_therm_p, a_RHS_therm_p);
    }
  }

  for (int lev(0); lev < nlev(); lev++) {

    //- ∂ε/∂t = -(epf_new - epf_old) / dt
    MultiFab::Subtract( *a_rhs[lev], *a_depdt[lev], 0, 0, 1, ngrow);

    // REVERSE THE SIGN
    a_rhs[lev]->mult(-1.0,0,1,0);

    EB_set_covered(*a_rhs[lev], 0, 1, 0, 0.0);

    // COMMENTED FOR TESTING BEFORE CHANGE
    // if (a_rhs[0]->nGrow() > 0) {
    //   bcs().fillpatch(lev, 0., BCFillVar::none, GetVecOfPtrs(a_rhs));
    // }
  }

}



/* Ideal gas open-system ***************************************************
 *                                                                         *
 * Compute S = Sh / (ρ*cp*T) + (1/ρ)Σ_k[W/Wk - (hk/(c*T))*SXk]             *
 *                                                                         *
 *   Sh  := RHS of enthalpy equation                                       *
 *   SXk := RHS of k-th species equation                                   *
 *                                                                         *
 *   cp  := mixture specific heat                                          *
 *   W   := mixture molecular weight                                       *
 *   T   := Temperature                                                    *
 *   Wk  := k-th species molecular weight                                  *
 *   hk  := k-th species specific enthalpy                                 *
 *                                                                         *
 ***************************************************************************/
void
mfix::
idealgas_opensystem_rhs ( Vector< MultiFab*       > const& a_rhs,
                          Vector< MultiFab const* > const& a_enthalpy_rhs,
                          Vector< MultiFab const* > const& a_species_rhs,
                          Vector< MultiFab const* > const& a_rho,
                          Vector< MultiFab const* > const& a_Tf,
                          Vector< MultiFab const* > const& a_Xf)
{
  int const solve_enthalpy( fluid.solve_enthalpy() );
  int const solve_species(  fluid.solve_species()  );

  int const nspecies = fluid.isMixture() ? fluid.nspecies() : 0;

  const auto& params = fluid.parameters<run_on>();
  const auto props = fluid.props.data<run_on>();

  for (int lev(0); lev < nlev(); lev++) {

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(a_rho[lev]->Factory());
    const auto& flagsFab = factory.getMultiEBCellFlagFab();

    for (MFIter mfi(*a_rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      Array4<const Real> empty_arr;

      auto const& rhs   = a_rhs[lev]->array(mfi);

      auto const& h_RHS = solve_enthalpy ? a_enthalpy_rhs[lev]->const_array(mfi) : empty_arr;
      auto const& X_RHS = solve_species  ? a_species_rhs[lev]->const_array(mfi)  : empty_arr;

      auto const& rho = a_rho[lev]->const_array(mfi);
      auto const& Tf  = a_Tf[lev]->const_array(mfi);
      auto const& Xf  = (nspecies > 0) ? a_Xf[lev]->const_array(mfi) : empty_arr;

      auto const& flags = flagsFab.const_array(mfi);

      ParallelFor(bx, [rhs, h_RHS, X_RHS, rho, Tf, Xf, flags, nspecies,
          solve_enthalpy, solve_species, params, props]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (!flags(i,j,k).isCovered()) {

          const Real rho_loc = rho(i,j,k);
          const Real Tf_loc  = Tf(i,j,k);

          // set initial fluid molecular weight
          Real MW_loc(0.);
          Real cp_loc = (!solve_enthalpy) ? 0.0 :
              props.specificHeat(Tf_loc, IntVect(i,j,k), Xf);

          if (nspecies > 0) {

            for (int n(0); n < nspecies; ++n) {

              const Real MW_k = params.get_MW_gk(n);

              MW_loc += Xf(i,j,k,n) / MW_k;
            }

            MW_loc = Real(1.) / MW_loc;

          }

          if (solve_enthalpy) {

            rhs(i,j,k) += h_RHS(i,j,k) / (rho_loc*cp_loc*Tf_loc);

          }

          if (nspecies > 0) {

            for (int n(0); n < nspecies; ++n) {

              const Real MW_k = params.get_MW_gk(n);
              Real coeff = (MW_loc / MW_k);

              if (solve_enthalpy) {
                const Real h_k = props.enthalpy(n, Tf_loc);
                coeff -= (h_k / (cp_loc*Tf_loc));
              }

              rhs(i,j,k) += (coeff / rho_loc) * X_RHS(i,j,k,n);
            }

          } else if (solve_species) {

            const Real h_loc = props.enthalpy(Tf_loc, nullptr);

            Real coeff = Real(1.0) - (h_loc / (cp_loc*Tf_loc));

            rhs(i,j,k) += (coeff / rho_loc) * X_RHS(i,j,k,0);
          }
        }
      });
    } // MFIter
  } // nlev

}


/* Ideal gas closed-system *************************************************
 *                                                                         *
 * Compute: δS - δθ(<S>/<θ>)                                               *
 *                                                                         *
 * S is the ideal gas open-system constraint RHS and should already be     *
 * computed and stored in a_rhs before this routine is called.             *
 *                                                                         *
 *   θ = ε/(Γ.pT) = (ε/pT)*(1 - 1/(cp*ρ*T))                                *
 *                                                                         *
 *   S = <S> + δS      <x> := average of x over domain (scalar)            *
 *   θ = <θ> + δθ       δx := fluctuation of x (field)                     *
 *                                                                         *
 *   cp := mixture specific heat                                           *
 *   T  := Temperature                                                     *
 *   ρ  := density                                                         *
 *                                                                         *
 ***************************************************************************/
void
mfix::
idealgas_closedsystem_rhs ( Vector< MultiFab*       > const& a_rhs,
                            Vector< MultiFab const* > const& a_epf,
                            Vector< MultiFab const* > const& a_rho,
                            Vector< MultiFab const* > const& a_Tf,
                            Vector< MultiFab const* > const& a_Xf,
                            Real const& a_therm_p,
                            Real      & a_RHS_therm_p)
{
  const int nspecies   = fluid.nspecies();
  const int is_mixture = fluid.isMixture();

  const auto& params = fluid.parameters<run_on>();
  const auto props = fluid.props.data<run_on>();

  Vector <std::unique_ptr<MultiFab>> Theta(nlev());

  ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
  ReduceData<Real,Real,Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  // Compute Theta
  for (int lev(0); lev < nlev(); ++lev) {

    int const ngrow(0);

    Theta[lev].reset( new MultiFab(grids[lev], dmap[lev], 1, ngrow,
         MFInfo(), a_rhs[lev]->Factory()));

    Theta[lev]->setVal(0.);

    // Mask coarser data on lower levels by skipping covered cells. This is
    // performed by setting the fine mask to 1 and the coarse mask to 0.
    iMultiFab finemask;
    if (lev < finest_level ) {

      finemask = makeFineMask(grids[lev-1], dmap[lev-1], grids[lev],
          ref_ratio[lev-1], 1, 0);

    } else {

      finemask.define(grids[lev], dmap[lev], 1, 0);
      finemask.setVal(1);

    }

    const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(a_rho[lev]->Factory());
    const auto& flags_fab = factory.getMultiEBCellFlagFab();

    const MultiFab& volfrac =  factory.getVolFrac();

    Real const& pT(a_therm_p);

    // Compute Theta and sum volume, and volume weighted Theta and Sigma
    // over the domain.

    for (MFIter mfi(*Theta[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      auto const& flags = flags_fab.const_array(mfi);
      Array4< int const > const& mask = finemask.const_array(mfi);

      Array4< Real const > const& vfrac = volfrac.const_array(mfi);
      Array4< Real const > const& epf   = a_epf[lev]->const_array(mfi);
      Array4< Real const > const& rho   = a_rho[lev]->const_array(mfi);
      Array4< Real const > const& Xf    = is_mixture ? a_Xf[lev]->const_array(mfi) : Array4< Real const>{};
      Array4< Real const > const& Tf    = a_Tf[lev]->const_array(mfi);

      Array4< Real const > const& sigma = a_rhs[lev]->const_array(mfi);
      Array4< Real       > const& theta = Theta[lev]->array(mfi);

      reduce_op.eval(bx, reduce_data, [sigma, theta, vfrac, epf,
          rho, Xf, Tf, pT, mask, flags, is_mixture, nspecies, params, props]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        Real sumVol(0.);
        Real sumVolTheta(0.);
        Real sumVolSigma(0.);

        if (mask(i,j,k) && !flags(i,j,k).isCovered()) {

          Real cp = props.specificHeat(IntVect(i,j,k), Tf, Xf);

          // set initial fluid molecular weight
          theta(i,j,k) = (epf(i,j,k)/pT)*(1_rt - (1_rt/(cp*rho(i,j,k)*Tf(i,j,k))));

          sumVol = vfrac(i,j,k)*epf(i,j,k);
          sumVolTheta = sumVol*theta(i,j,k);
          sumVolSigma = sumVol*sigma(i,j,k);
        }
        return {sumVol, sumVolTheta, sumVolSigma};
      });
    }
  }// lev

  ReduceTuple host_tuple = reduce_data.value(reduce_op);

  Vector<Real> volSums = {amrex::get<0>(host_tuple),
                          amrex::get<1>(host_tuple),
                          amrex::get<2>(host_tuple)};

  ParallelDescriptor::ReduceRealSum(volSums.dataPtr(), volSums.size());

  Real const avgTheta( volSums[1]/volSums[0] );
  Real const avgSigma( volSums[2]/volSums[0] );

  Real const dpdt( avgSigma/avgTheta );

  a_RHS_therm_p = dpdt;

  // Finalize RHS
  for (int lev(0); lev < nlev(); ++lev) {

    for (MFIter mfi(*a_rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      Array4< Real       > const& rhs   = a_rhs[lev]->array(mfi);
      Array4< Real const > const& theta = Theta[lev]->const_array(mfi);

      amrex::ParallelFor(bx, [rhs, theta, dpdt, avgSigma, avgTheta]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real const S(rhs(i,j,k)  );
        Real const T(theta(i,j,k));

        Real const dS(S - avgSigma);
        Real const dT(T - avgTheta);

        rhs(i,j,k) = dS - dpdt*dT;
      });
    }
  }
}
