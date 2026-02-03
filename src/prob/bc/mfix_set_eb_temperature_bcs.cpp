#include <AMReX_FillPatchUtil.H>

#include <mfix_calc_cell.H>
#include <mfix_bc.H>
#include <mfix_fluid.H>
#include <mfix_least_squares_K.H>
#include <mfix_bc_fillpatch_K.H>

using namespace amrex;

void
MFIXBoundaryConditions::
set_eb_temperature_bcs ( int const a_lev, Real const a_time,
                         Vector<MultiFab      *> a_Teb,
                         MultiFab const* a_epf,
                         MultiFab const* a_Tf,
                         int const a_avg_Tp_comp,
                         MultiFab const* a_avg_Tp)
{
  BL_PROFILE("MFIXBoundaryConditions::set_eb_temperature_bcs()");

  if (a_lev == 0) {

    a_Teb[a_lev]->setVal(0.);

    const auto& factory =
      dynamic_cast<EBFArrayBoxFactory const&>(a_Teb[a_lev]->Factory());

    const auto& flags_fab = factory.getMultiEBCellFlagFab();

    const auto& cellcent  = factory.getCentroid();
    const auto& bndrycent = factory.getBndryCent();
    const auto& bndrynorm = factory.getBndryNormal();

    if (a_avg_Tp_comp >= 0) { AMREX_ASSERT( a_avg_Tp ); }

    for (MFIter mfi(*a_Teb[a_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();
      FabType t = flags_fab[mfi].getType(bx);

      // Only update cut-cell values
      if (t != FabType::singlevalued) { continue; }

      // EB arrays
      Array4<Real      > const& Teb = a_Teb[a_lev]->array(mfi);
      Array4<Real const> const& epf = a_epf->const_array(mfi);
      Array4<Real const> const& Tf  = a_Tf->const_array(mfi);

      Array4<Real const> const& avg_Tp = (a_avg_Tp_comp >= 0) ?
          a_avg_Tp->const_array(mfi, a_avg_Tp_comp) : Array4<Real const>{};

      Array4<const EBCellFlag> const& flags = flags_fab.const_array(mfi);

      // Cell centroid, EB centroid, EB Normal
      Array4<const Real> const& ccent = cellcent.const_array(mfi);
      Array4<const Real> const& bcent = bndrycent.const_array(mfi);
      Array4<const Real> const& bnorm = bndrynorm.const_array(mfi);

      // Loop over BCs
      for (int bcv(0); bcv < m_bc.size(); ++bcv) {

        // Only work on boundaries that are EB
        if (bc(bcv).type != BCList::eb) { continue; }

        const Box eb_box = calc_ic_box(m_geom[a_lev], bc(bcv).region);

        // Only work on where this box intersects the EB's box
        if (!eb_box.intersects(bx)) { continue; }

        // Intersection of eb box and mfi box
        const Box bx_int = bx & eb_box;

        const int  has_normal = bc(bcv).eb.has_normal;
        amrex::GpuArray<amrex::Real,3> normal{0.};
        if (has_normal) {
          normal[0] = bc(bcv).eb.normal[0];
          normal[1] = bc(bcv).eb.normal[1];
          normal[2] = bc(bcv).eb.normal[2];
        }

        const Real pad = std::numeric_limits<float>::epsilon();
        const Real normal_tol = bc(bcv).eb.normal_tol;

        const Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
        const Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

        if ( bc(bcv).eb.constant_temperature() ) {

          // User-defined constant wall temperature in this region.
          const Real eb_temperature = bc(bcv).eb.temperature(a_time);

          ParallelFor(bx_int, [flags, eb_temperature, Teb, has_normal,
            normal, bnorm, norm_tol_lo, norm_tol_hi]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            if (flags(i,j,k).isSingleValued()) {

              Real mask = Real(1.0);

              if (has_normal) {

                const Real dotprod = bnorm(i,j,k,0)*normal[0]
                                   + bnorm(i,j,k,1)*normal[1]
                                   + bnorm(i,j,k,2)*normal[2];

                mask = ((norm_tol_lo <= dotprod) &&
                        (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
              }

              Teb(i,j,k) = mask*eb_temperature;
            }
          });

        } else if ( bc(bcv).eb.conjugate_heat_transfer_1D() ) {

          Real thermal_resistance(0.);

          // convective heat transfer: fluid-wall (internal)
          thermal_resistance += 1./bc(bcv).eb.h_int;

          // conduction through the wall
          if (bc(bcv).eb.L_wall > 0.) {
            thermal_resistance += bc(bcv).eb.L_wall/bc(bcv).eb.k_wall;
          }

          // convective heat transfer: wall-environment (external)
          thermal_resistance += 1./bc(bcv).eb.h_ext;

          Real const R_inv = (almostEqual(thermal_resistance,0.))
              ? 0. : (1./thermal_resistance);

          Real const hf_inv = (1./bc(bcv).eb.h_int);
          Real const T_env = bc(bcv).eb.T_env;

          bool const phase_average_temperature = bc(bcv).eb.phase_averaged_temperature;
          if (phase_average_temperature ) { AMREX_ASSERT(a_avg_Tp_comp >= 0); }

          ParallelFor(bx_int, [flags, hf_inv, T_env, R_inv, Teb, epf, Tf, avg_Tp,
            phase_average_temperature, has_normal, normal, norm_tol_lo, norm_tol_hi,
            bnorm, bcent, ccent]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            if (flags(i,j,k).isSingleValued()) {

              Real mask = Real(1.);
              Real T0 = Real(0.);

              if (has_normal) {

                const Real dotprod = bnorm(i,j,k,0)*normal[0]
                                   + bnorm(i,j,k,1)*normal[1]
                                   + bnorm(i,j,k,2)*normal[2];

                mask = ((norm_tol_lo <= dotprod) &&
                        (dotprod <= norm_tol_hi)) ? Real(1.0) : Real(0.0);
              }

              if (mask > 0.) {

                T0 = !phase_average_temperature ? Tf(i,j,k) :
                  epf(i,j,k)*Tf(i,j,k) + (1.-epf(i,j,k))*avg_Tp(i,j,k);

                // Compute the location where we want to 'sample' the
                // bulk fluid temperature. This is 1 dx away frm the eb,
                // normal to the EB face.
                Array1D<Real,0,2> pt = { bcent(i,j,k,0) - bnorm(i,j,k,0),
                                         bcent(i,j,k,1) - bnorm(i,j,k,1),
                                         bcent(i,j,k,2) - bnorm(i,j,k,2) };

                // Compute the index to neighbor cell that contains the point.
                IntVect offset = { ((pt(0) >= 0.5) ? 1 : (pt(0) <= -0.5 ? -1 : 0))
                                 , ((pt(1) >= 0.5) ? 1 : (pt(1) <= -0.5 ? -1 : 0))
                                 , ((pt(2) >= 0.5) ? 1 : (pt(2) <= -0.5 ? -1 : 0)) };

                if (flags(i,j,k).isConnected(offset)) {

                  // Update the position relative to cell.
                  pt(0) -= Real(offset[0]);
                  pt(1) -= Real(offset[1]);
                  pt(2) -= Real(offset[2]);

                  // cell is the (i,j,k) of cell containing pt. We look to collect
                  // information from all the neighboring cells.
                  IntVect cell = IntVect(i,j,k) + offset;

                  int constexpr maxPts( 27 ); // stencil is (-1,0,1)
                  int constexpr minPts( 10 ); // minimum needed for least square solve

                  int numPts = 0;

                  Array2D<Real, 0,maxPts-1, 0,2> points;
                  Array1D<Real, 0,maxPts-1> vals;

                  for ( int kk=-1; kk<=1; ++kk ) {
                    for ( int jj=-1; jj<=1; ++jj ) {
                      for ( int ii=-1; ii<=1; ++ii ) {

                        IntVect ngh_cell = cell + IntVect(ii,jj,kk);

                        if ( flags(cell).isConnected(ii,jj,kk) ) {

                          points(numPts,0) = Real(ii) + ccent(ngh_cell,0);
                          points(numPts,1) = Real(jj) + ccent(ngh_cell,1);
                          points(numPts,2) = Real(kk) + ccent(ngh_cell,2);

                          vals(numPts) = !phase_average_temperature ? Tf(i,j,k) :
                            epf(i,j,k)*Tf(i,j,k) + (1.-epf(i,j,k))*avg_Tp(i,j,k);

                          numPts++;
                        }
                      }
                    }
                  }

                  AMREX_ASSERT( numPts <= maxPts );

                  if (numPts > minPts) {
                    auto VecVal = LeastSquares<maxPts>(points, vals, pt, numPts);
                    T0 = VecVal(0);
                  }

                } // connected
              } // !masked


              Real flux = (T0 - T_env) * R_inv;

              Teb(i,j,k) = mask * (T0 - flux * hf_inv);
            }
          });
        }

      }//bcv loop
    } // end MFIter loop

  } else {

    Vector<BCRec> bcrec(1);

    for (int idim(0); idim<AMREX_SPACEDIM; ++idim) {
      if ( (m_geom[0].isPeriodic(idim)) ) {
        for (auto& b : bcrec) b.setLo(idim, BCType::int_dir);
        for (auto& b : bcrec) b.setHi(idim, BCType::int_dir);
      } else {
        for (auto& b : bcrec) b.setLo(idim, BCType::foextrap);
        for (auto& b : bcrec) b.setHi(idim, BCType::foextrap);
      }
    }

    PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > cphysbc
        (m_geom[a_lev-1], bcrec, MFIXForFill{/*probtype=*/-1});

    PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > fphysbc
        (m_geom[a_lev  ], bcrec, MFIXForFill{/*probtype=*/-1});

    Interpolater* mapper = &pc_interp;

    InterpFromCoarseLevel( *a_Teb[a_lev], a_time, *a_Teb[a_lev-1],
        /*scomp=*/0, /*dcomp=*/0, /*ncomp=*/1,
        m_geom[a_lev-1], m_geom[a_lev], cphysbc, 0, fphysbc, 0,
        m_ref_ratio[a_lev-1], mapper, bcrec, /*physbc_comp=*/0);
  }

  // Do this after as well as before to pick up terms that got updated in the
  // call above
  a_Teb[a_lev]->FillBoundary(m_geom[a_lev].periodicity());
}
