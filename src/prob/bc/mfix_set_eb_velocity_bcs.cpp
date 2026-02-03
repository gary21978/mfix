#include <mfix_calc_cell.H>
#include <mfix_bc.H>
#include <mfix_fluid.H>

#include <AMReX_ParallelReduce.H>

using namespace amrex;

void
MFIXBoundaryConditions::
set_eb_velocity_bcs (Real time_in,
                     Vector< MultiFab* > const& eb_vel_g)
{
  BL_PROFILE("MFIXBoundaryConditions::set_eb_velocity_bcs()");

  const int nlev = eb_vel_g.size();

  for (int lev = 0; lev < nlev; lev++)
  {
     eb_vel_g[lev]->setVal(0.);

     const auto& factory =
       dynamic_cast<EBFArrayBoxFactory const&>(eb_vel_g[lev]->Factory());

     const auto& flags = factory.getMultiEBCellFlagFab();

     for (MFIter mfi(*eb_vel_g[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

       const Box& bx = mfi.tilebox();
       FabType t = flags[mfi].getType(bx);

       // We update only cut-cells values
       if (t == FabType::singlevalued) {

           const auto &eb_vel_g_arr = (*eb_vel_g[lev])[mfi].array();
           const auto &flags_arr    = factory.getMultiEBCellFlagFab()[mfi].const_array();

           const auto &eb_norm_arr  = factory.getBndryNormal()[mfi].const_array();

           for (int bcv(0); bcv < m_bc.size(); ++bcv) {

             if (m_bc[bcv].type == BCList::eb) {

               const Box ic_bx = calc_ic_box(m_geom[lev], m_bc[bcv].region);

               if (ic_bx.intersects(bx)) {

                 // Intersection of ic box and mfi box
                 const Box bx_int = bx&(ic_bx);

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

                 int has_comps = 0;

                 Real vel_mag(0.0);
                 GpuArray<amrex::Real,3> vel_comps{0.};

                 // Flow is specified as a velocity magnitude
                 if ( m_bc[bcv].fluid.velocity.size() == 3 ) {

                    has_comps = 1;
                    const auto& bc_vels = m_bc[bcv].fluid.get_velocity(time_in);
                    vel_comps[0] = bc_vels[0];
                    vel_comps[1] = bc_vels[1];
                    vel_comps[2] = bc_vels[2];

                 } else if ( m_bc[bcv].fluid.velocity.size() == 1 ) {

                    vel_mag = m_bc[bcv].fluid.get_velocity_mag(time_in);

                 } else {

                     Real volflow(0.);
                     if ( m_bc[bcv].fluid.massflow.is_defined() ) {

                       Real massflow = m_bc[bcv].fluid.get_massflow(time_in);
                       Real density = m_bc[bcv].fluid.get_density(time_in);

                       // check density
                       if (amrex::almostEqual(density, 0.)) {
                       }

                       volflow = massflow / density;

                     } else if ( m_bc[bcv].fluid.volflow.is_defined() ) {

                       volflow = m_bc[bcv].fluid.get_volflow(time_in);

                     } else {
                       // Error check
                     }

                    Real eb_area = get_bc_area(bcv);
                    eb_area *= m_bc[bcv].fluid.volfrac;

                    if (amrex::almostEqual(eb_area, 0.)) {
                      // Error check
                    }

                   vel_mag = volflow / eb_area;

                 }

                 ParallelFor(bx_int, [flags_arr,eb_vel_g_arr, has_normal, normal,
                 norm_tol_lo, norm_tol_hi, has_comps, vel_comps, vel_mag, eb_norm_arr]
                 AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 {
                   if (flags_arr(i,j,k).isSingleValued()) {

                     Real mask = Real(1.0);

                     if(has_normal) {
                       const Real dotprod = eb_norm_arr(i,j,k,0)*normal[0]
                                          + eb_norm_arr(i,j,k,1)*normal[1]
                                          + eb_norm_arr(i,j,k,2)*normal[2];

                       mask = ((norm_tol_lo <= dotprod) &&
                               (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
                     }


                     if(has_comps) {
                       eb_vel_g_arr(i,j,k,0) = mask*vel_comps[0];
                       eb_vel_g_arr(i,j,k,1) = mask*vel_comps[1];
                       eb_vel_g_arr(i,j,k,2) = mask*vel_comps[2];

                     } else {

                       // The EB normal points out of the domain so we need to flip the
                       // when using it to convert magnitude to velocity components so
                       // the resulting vector points into the domain.
                       eb_vel_g_arr(i,j,k,0) = -mask*eb_norm_arr(i,j,k,0)*vel_mag;
                       eb_vel_g_arr(i,j,k,1) = -mask*eb_norm_arr(i,j,k,1)*vel_mag;
                       eb_vel_g_arr(i,j,k,2) = -mask*eb_norm_arr(i,j,k,2)*vel_mag;
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
     eb_vel_g[lev]->FillBoundary(m_geom[lev].periodicity());
  }
}
