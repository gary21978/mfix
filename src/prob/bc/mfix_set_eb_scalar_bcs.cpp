#include <mfix_calc_cell.H>
#include <mfix_bc.H>
#include <mfix_fluid.H>

using namespace amrex;

void
MFIXBoundaryConditions::
set_eb_scalar_bcs ( Vector< MultiFab* > const& eb_scalars,
                    Vector< MultiFab* > const& eb_species,
                    int const a_has_enthalpy,
                    int const a_nspecies)
{
  BL_PROFILE("MFIXBoundaryConditions::set_eb_scalar_bcs()");

  int const nlev = eb_scalars.size();
  int const stride( m_bc.size() );
  int const has_enthalpy(a_has_enthalpy);

  for (int lev = 0; lev < nlev; lev++)
  {
    eb_scalars[lev]->setVal(0.);

    const auto& factory =
      dynamic_cast<EBFArrayBoxFactory const&>(eb_scalars[lev]->Factory());

    const auto& flags = factory.getMultiEBCellFlagFab();

    for (MFIter mfi(*eb_scalars[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();
      FabType t = flags[mfi].getType(bx);

      // We update only cut-cells values
      if (t == FabType::singlevalued) {

        const auto &eb_scalars_arr = (*eb_scalars[lev])[mfi].array();
        const auto &eb_species_arr = (a_nspecies > 0) ? (*eb_species[lev])[mfi].array() : Array4<Real>{};

        const auto &flags_arr    = factory.getMultiEBCellFlagFab()[mfi].const_array();
        const auto &eb_norm_arr  = factory.getBndryNormal()[mfi].const_array();

        for (int bcv(0); bcv < m_bc.size(); ++bcv) {

          if (m_bc[bcv].type == BCList::eb && m_bc[bcv].fluid.flow_thru_eb) {

            const Box ic_bx = calc_ic_box(m_geom[lev], m_bc[bcv].region);

            if (ic_bx.intersects(bx)) {

              // Intersection of ic box and mfi box
              const Box bx_int = bx & ic_bx;

              const int  has_normal = m_bc[bcv].eb.has_normal;
              amrex::GpuArray<amrex::Real,3> normal{0.};
              if (has_normal) {
                normal[0] = m_bc[bcv].eb.normal[0];
                normal[1] = m_bc[bcv].eb.normal[1];
                normal[2] = m_bc[bcv].eb.normal[2];
              }

              const Real pad = std::numeric_limits<float>::epsilon();
              const Real normal_tol = m_bc[bcv].eb.normal_tol;

              const Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
              const Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

              // (Rho), (Rho*Enthalpy), (Rho*Tracer)
              Real* p_bc_rho   = m_bc_rho.data();
              Real* p_bc_h     = m_bc_h.data();
              Real* p_bc_trac  = m_bc_tracer.data();
              Real* p_bc_X     = m_bc_Xk.data();

              const int num_trac = eb_scalars[lev]->nComp()-2;
              const int nspecies = a_nspecies;

              ParallelFor(bx_int, [bcv, num_trac, flags_arr, eb_scalars_arr,
              nspecies, eb_species_arr, eb_norm_arr, has_enthalpy,
              has_normal, normal, norm_tol_lo, norm_tol_hi,
              p_bc_rho, p_bc_h, p_bc_trac, p_bc_X, stride]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
               if (flags_arr(i,j,k).isSingleValued()) {

                  Real mask = Real(1.0);

                  if (has_normal) {

                    const Real dotprod = eb_norm_arr(i,j,k,0)*normal[0]
                                       + eb_norm_arr(i,j,k,1)*normal[1]
                                       + eb_norm_arr(i,j,k,2)*normal[2];

                    mask = ((norm_tol_lo <= dotprod) &&
                            (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
                  }

                  eb_scalars_arr(i,j,k,0) = mask*p_bc_rho[bcv];

                  eb_scalars_arr(i,j,k,1) = has_enthalpy
                      ? mask*p_bc_rho[bcv]*p_bc_h[bcv] : 0.;

                  for(int n(0); n<num_trac; n++) {
                    eb_scalars_arr(i,j,k,2+n) = mask*p_bc_rho[bcv]*p_bc_trac[bcv];
                  }

                  for(int n(0); n<nspecies; n++) {
                    eb_species_arr(i,j,k,n) = mask*p_bc_rho[bcv]*p_bc_X[bcv + n*stride];
                  }
                }
              });
            } // the ic box intersects the mfi box
          } // bc type is EB
        } // loop over bcv

      } // single valued fabs only
    }

    // Do this after as well as before to pick up terms that got updated in the
    // call above
    eb_scalars[lev]->FillBoundary(m_geom[lev].periodicity());
  }
}
