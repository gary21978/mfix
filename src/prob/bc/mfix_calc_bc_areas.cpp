#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>

#include <mfix_calc_cell.H>
#include <mfix_bc.H>

using namespace amrex;

void MFIXBoundaryConditions::
calc_bc_areas (int const lev,
               BoxArray            const& a_grids,
               DistributionMapping const& a_dmap,
               EBFArrayBoxFactory*        a_factory)
{
  BL_PROFILE("MFIXBoundaryConditions::calc_bc_areas()");

  if (m_fab_area.size() > 0) { this->m_fab_area.clear(); }

  int const num_bcs(m_bc.size());

  // Create a temp MultiFab to use in MFIter loops.
  // We want to 'save' the individual fab areas if PIC or DEM are
  // enabled because we can use that information to compute the
  // actual particle volumetric flow rate.
  MultiFab tmpMF(a_grids, a_dmap, 1, 0, MFInfo().SetAlloc(false), *a_factory);

  this->m_area.resize(num_bcs);
  this->m_fab_area.resize(num_bcs);

  std::vector<Real> bc_areas(num_bcs, 0.);

  Box domain(m_geom[lev].Domain());

  if( !m_geom[lev].isPeriodic(0) ) { // begin x direction

    int const dir(0);

    { // x-lo side of the domain

      Box box_ilo = amrex::adjCellLo(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_xlo().size(); ++lc) {

        const int bcv  = this->bc_xlo(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, bcv, box_ilo, bc, tmpMF, a_factory);
      } // loop over bc_xlo

    } // end x-lo side of the domain

    { // x-hi side of the domain
      Box box_ihi = amrex::adjCellHi(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_xhi().size(); ++lc) {

        const int bcv  = this->bc_xhi(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, bcv, box_ihi, bc, tmpMF, a_factory);
      }
    } // end x-hi side of the domain
  } // end x-direction


  if( !m_geom[lev].isPeriodic(1) ) { // begin y-direction

    const int dir = 1;

    { // y-lo side of the domain

      Box box_jlo = amrex::adjCellLo(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_ylo().size(); ++lc) {

        const int bcv  = this->bc_ylo(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, bcv, box_jlo, bc, tmpMF, a_factory);
      }

    }// end y-lo side of the domain

    { // y-hi side of the domain

      Box box_jhi = amrex::adjCellHi(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_yhi().size(); ++lc) {

        const int bcv  = this->bc_yhi(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, bcv, box_jhi, bc, tmpMF, a_factory);
      }

    } // end y-hi side of the domain
  } // end y-direction

  { // begin z-direction

    const int dir(2);

    { // z-lo side of the domain

      Box box_klo = amrex::adjCellLo(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_zlo().size(); ++lc) {

        const int bcv  = this->bc_zlo(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, bcv, box_klo, bc, tmpMF, a_factory);
      }


    } // end z-lo side of the domain


    { // z-hi side of the domain

      Box box_khi = amrex::adjCellHi(domain,dir,1);

      // Define specific BC conditions from inputs
      for (int lc(0); lc < bc_zhi().size(); ++lc) {

        const int bcv  = this->bc_zhi(lc);
        const BC_t& bc = this->bc(bcv);

        bc_areas[bcv] = calc_bc_area(lev, dir, bcv, box_khi, bc, tmpMF, a_factory);
      }

    } // end z-hi side of the domain
  }// end z-direction


  for (int bcv(0); bcv<num_bcs; bcv++) {

    BC_t const& bc = this->bc(bcv);

    if (bc.type == BCList::eb) {
      bc_areas[bcv] = calc_eb_bc_area(lev, bcv, bc, tmpMF, a_factory);
    }
  }


  // Do a global reduce sum and store in data container
  if (num_bcs > 0)
    ParallelDescriptor::ReduceRealSum(bc_areas.data(), bc_areas.size());

  for(int bcv(0); bcv < num_bcs; bcv++) {
    amrex::Print() << "BC: " << bcv << "  area: " << bc_areas[bcv] << "\n";
    this->m_area[bcv] = bc_areas[bcv];
  }

}


Real
MFIXBoundaryConditions::
calc_bc_area( int      const       a_lev,
              int      const       a_dir,
              int      const       a_bcv,
              Box      const&      a_bx,
              BC_t     const&      a_bc,
              MultiFab const&      a_MF,
              EBFArrayBoxFactory*  a_factory)
{
  Real dx = m_geom[a_lev].CellSize(0);
  Real dy = m_geom[a_lev].CellSize(1);
  Real dz = m_geom[a_lev].CellSize(2);

  // This should be caught elsewhere but just in case...
  AMREX_ALWAYS_ASSERT(amrex::almostEqual(dx,dy) && amrex::almostEqual(dy,dz));
  Real const da( dx*dx );

  const GpuArray<Real, 3> plo = m_geom[a_lev].ProbLoArray();

  IntVect bx_lo(a_bx.loVect());
  IntVect bx_hi(a_bx.hiVect());

  if (a_dir != 0){
    bx_lo[0] = static_cast<int>(amrex::Math::floor((a_bc.region->lo(0)-plo[0])/dx + 0.5));
    bx_hi[0] = static_cast<int>(amrex::Math::floor((a_bc.region->hi(0)-plo[0])/dx + 0.5)-1);
  }
  if (a_dir != 1){
    bx_lo[1] = static_cast<int>(amrex::Math::floor((a_bc.region->lo(1)-plo[1])/dy + 0.5));
    bx_hi[1] = static_cast<int>(amrex::Math::floor((a_bc.region->hi(1)-plo[1])/dy + 0.5)-1);
  }
  if (a_dir != 2){
    bx_lo[2] = static_cast<int>(amrex::Math::floor((a_bc.region->lo(2)-plo[2])/dz + 0.5));
    bx_hi[2] = static_cast<int>(amrex::Math::floor((a_bc.region->hi(2)-plo[2])/dz + 0.5)-1);
  }

  const Box bc_bx(bx_lo, bx_hi);

  Real totalArea(0.);

  // No tiling
  for (MFIter mfi(a_MF, false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.growntilebox(IntVect::TheDimensionVector(a_dir));

    std::pair<int,int> index(mfi.index(), mfi.LocalTileIndex());
    this->m_fab_area[a_bcv][index] = 0.;

    // Ensure lo_box intersects bx
    if (bc_bx.intersects(bx)) {

      const Box bx_int = bx&(bc_bx);

      EBCellFlagFab const& flagfab = a_factory->getMultiEBCellFlagFab()[mfi];
      FabType t = flagfab.getType(bx);

      if (t == FabType::regular ) {

        this->m_fab_area[a_bcv][index] = da*static_cast<Real>(bx_int.numPts());

      } else if (t == FabType::singlevalued ) {

        Array4<Real const> const& afrac = (a_factory->getAreaFrac()[a_dir])->const_array(mfi);

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(bx_int, reduce_data, [afrac]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        { return afrac(i,j,k); });

        ReduceTuple host_tuple = reduce_data.value(reduce_op);
        this->m_fab_area[a_bcv][index] = da*amrex::get<0>(host_tuple);

      }
    } // a_bx and the MFIter box intersect

    totalArea += this->m_fab_area[a_bcv][index];
  } // MFIter loop

  return totalArea;
}


Real
MFIXBoundaryConditions::
calc_eb_bc_area( int      const       a_lev,
                 int      const       a_bcv,
                 BC_t     const&      a_bc,
                 MultiFab const&      a_MF,
                 EBFArrayBoxFactory*  a_factory)
{
  BL_PROFILE("MFIXBoundaryConditions::calc_eb_bc_areas()");

  Real dx = m_geom[a_lev].CellSize(0);
  Real dy = m_geom[a_lev].CellSize(1);
  Real dz = m_geom[a_lev].CellSize(2);

  // This should be caught elsewhere but just in case...
  AMREX_ALWAYS_ASSERT(amrex::almostEqual(dx,dy) && amrex::almostEqual(dy,dz));
  Real const da( dx*dx );

  const int  has_normal = a_bc.eb.has_normal;

  amrex::GpuArray<amrex::Real,3> normal{0.};
  if (has_normal) {
     normal[0] = a_bc.eb.normal[0];
     normal[1] = a_bc.eb.normal[1];
     normal[2] = a_bc.eb.normal[2];
  }

  const Real pad = std::numeric_limits<float>::epsilon();
  const Real normal_tol = a_bc.eb.normal_tol;

  const Real norm_tol_lo = Real(-1.) - (normal_tol + pad);
  const Real norm_tol_hi = Real(-1.) + (normal_tol + pad);

  const Box ic_bx = calc_ic_box(m_geom[a_lev], a_bc.region);

  Real totalArea(0.);

  for (MFIter mfi(a_MF, false); mfi.isValid(); ++mfi) {

    std::pair<int,int> index(mfi.index(), mfi.LocalTileIndex());
    this->m_fab_area[a_bcv][index] = 0.;

    const Box& bx = mfi.tilebox();

    // Ensure lo_box intersects bx
    if (ic_bx.intersects(bx)) {

      EBCellFlagFab const& flagfab = a_factory->getMultiEBCellFlagFab()[mfi];
      FabType t = flagfab.getType(bx);

      if (t == FabType::singlevalued) {

        // EB boundary area and normal
        Array4<Real const> const& eb_area = (a_factory->getBndryArea()).const_array(mfi);
        Array4<Real const> const& eb_norm = (a_factory->getBndryNormal()).const_array(mfi);
        Array4<EBCellFlag const> const& flags = flagfab.const_array();

        const Box bx_int = bx&(ic_bx);

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(bx_int, reduce_data, [flags, has_normal, normal,
          norm_tol_lo, norm_tol_hi, eb_area, eb_norm]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
          Real area(0.0);

          if(flags(i,j,k).isSingleValued()) {

            area = eb_area(i,j,k);

            if(has_normal) {

                const Real dotprod = eb_norm(i,j,k,0)*normal[0]
                                   + eb_norm(i,j,k,1)*normal[1]
                                   + eb_norm(i,j,k,2)*normal[2];

                area *= (norm_tol_lo <= dotprod) ? Real(1.0) : Real(0.0);
                area *= (dotprod <= norm_tol_hi) ? Real(1.0) : Real(0.0);
            }
          }
          return {area};
        });

        ReduceTuple host_tuple = reduce_data.value(reduce_op);
        this->m_fab_area[a_bcv][index] = da*amrex::get<0>(host_tuple);

      } // single valued
    } // box intersects
    totalArea += this->m_fab_area[a_bcv][index];
  } // loop over MFIter

  return totalArea;

}
