#include <mfix.H>

#include <mfix_fluid.H>
#include <mfix_bc.H>
#include <mfix_ic.H>

#include <AMReX_GpuDevice.H>

#include <string>
#include <climits>

void
mfix::set_gp0 ()
{
  gp0[0] = 0.0;
  gp0[1] = 0.0;
  gp0[2] = 0.0;

  int const lev(0);

  Box domain(geom[lev].Domain());

  const Real tolerance = std::numeric_limits<Real>::epsilon();

  const IntVect dom_lo(domain.loVect());
  const IntVect dom_hi(domain.hiVect());

  auto make_ifab_host = [] (IArrayBox const* dfab) -> IArrayBox
  {
      IArrayBox ifab_host(dfab->box(),dfab->nComp(),The_Pinned_Arena());
#ifdef AMREX_USE_GPU
      Gpu::dtoh_memcpy_async(ifab_host.dataPtr(), dfab->dataPtr(), ifab_host.nBytes());
#else
      std::memcpy(ifab_host.dataPtr(), dfab->dataPtr(), ifab_host.nBytes());
#endif
      return ifab_host;
  };

  IArrayBox const& bc_ilo_host = make_ifab_host(bc_list.bc_ilo[lev]);
  IArrayBox const& bc_ihi_host = make_ifab_host(bc_list.bc_ihi[lev]);
  IArrayBox const& bc_jlo_host = make_ifab_host(bc_list.bc_jlo[lev]);
  IArrayBox const& bc_jhi_host = make_ifab_host(bc_list.bc_jhi[lev]);
  IArrayBox const& bc_klo_host = make_ifab_host(bc_list.bc_klo[lev]);
  IArrayBox const& bc_khi_host = make_ifab_host(bc_list.bc_khi[lev]);
  Gpu::synchronize();

  Array4<const int> const& bct_ilo = bc_ilo_host.const_array();
  Array4<const int> const& bct_ihi = bc_ihi_host.const_array();
  Array4<const int> const& bct_jlo = bc_jlo_host.const_array();
  Array4<const int> const& bct_jhi = bc_jhi_host.const_array();
  Array4<const int> const& bct_klo = bc_klo_host.const_array();
  Array4<const int> const& bct_khi = bc_khi_host.const_array();

  int delp_dir_loc(m_boundary_conditions.delp_dir());

  Real delp_x = m_boundary_conditions.delp(0);
  Real delp_y = m_boundary_conditions.delp(1);
  Real delp_z = m_boundary_conditions.delp(2);

  // ---------------------------------------------------------------->>>
  //     If the bc's are pressure inflow/outflow then be sure to
  //     capture that in p0 and gp0.
  // ---------------------------------------------------------------->>>

  if ( ( (bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],0) == BCList::pinf) &&
         (bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],0) == BCList::pout) )
       ||
       ( (bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],0) == BCList::pinf) &&
         (bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],0) == BCList::pout) ) )
    {

      const int bcv_lo = bct_ilo(dom_lo[0]-1,dom_lo[1],dom_lo[2],1);
      const int bcv_hi = bct_ihi(dom_hi[0]+1,dom_lo[1],dom_lo[2],1);

      const Real  p_lo = m_boundary_conditions.h_bc_p(bcv_lo);
      const Real  p_hi = m_boundary_conditions.h_bc_p(bcv_hi);

      delp_dir_loc = 0;
      delp_x = p_lo - p_hi;
    }


  else if( ( (bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],0) == BCList::pinf) &&
             (bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],0) == BCList::pout) )
           ||
           ( (bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],0) == BCList::pinf) &&
             (bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],0) == BCList::pout) ) )
    {

      const int bcv_lo = bct_jlo(dom_lo[0],dom_lo[1]-1,dom_lo[2],1);
      const int bcv_hi = bct_jhi(dom_lo[0],dom_hi[1]+1,dom_lo[2],1);

      const Real  p_lo = m_boundary_conditions.h_bc_p(bcv_lo);
      const Real  p_hi = m_boundary_conditions.h_bc_p(bcv_hi);

      delp_dir_loc = 1;
      delp_y = p_lo - p_hi;
    }

  else if( ( (bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,0) == BCList::pinf) &&
             (bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,0) == BCList::pout) )
           ||
           ( (bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,0) == BCList::pinf) &&
             (bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,0) == BCList::pout) ) )
    {

      const int bcv_lo = bct_klo(dom_lo[0],dom_lo[1],dom_lo[2]-1,1);
      const Real p_lo  = m_boundary_conditions.h_bc_p(bcv_lo);

      const int bcv_hi = bct_khi(dom_lo[0],dom_lo[1],dom_hi[2]+1,1);
      const Real p_hi  = m_boundary_conditions.h_bc_p(bcv_hi);

      delp_dir_loc = 2;
      delp_z = p_lo - p_hi;
  }

  // ---------------------------------------------------------------->>>

  //  Make sure that ic_p_g is set if using delp pressure conditions
  for(int icv(0); icv < m_initial_conditions.ic().size(); ++icv)
  {

    if((delp_dir_loc >= 0) && (delp_dir_loc == m_boundary_conditions.delp_dir())) {

      if (!m_initial_conditions.ic(icv).fluid.pressure_defined) {
        std::cout << "MUST DEFINE ic_p_g if using the DELP pressure condition" << std::endl;
        exit(0);
      }

    }
    else if((delp_dir_loc >= 0) && (delp_dir_loc != m_boundary_conditions.delp_dir())) {
      if( m_initial_conditions.ic(icv).fluid.pressure_defined ) {
        std::cout << "MUST not define ic_p_g if setting p_inflowandp_outflow" << std::endl;
        exit(0);
      }

    } else {

      const Real gravity_square_module =
        gravity[0]*gravity[0] + gravity[1]*gravity[1] + gravity[2]*gravity[2];

      if( !m_initial_conditions.ic(icv).fluid.pressure_defined || gravity_square_module > tolerance)
        {

          // HACK: This should probably take into consideration
          // variable fluid density.
          const Real ro_g0 = m_initial_conditions.ic(icv).fluid.get_density();

          for (int dim=0; dim<3; dim++){
            if (gravity[dim]*gravity[dim] > tolerance ) {
              gp0[dim] = ro_g0 * gravity[dim];
            }
          }

        return;
      }
    }
  }

  // ---------------------------------------------------------------->>>

  // Here the pressure in each cell is determined from a specified pressure
  // drop across the domain le>. This section requires that the pressure
  // is already defined in all initial condition regions (otherwise this
  // section would be skipped)

  Real xlen = geom[0].ProbHi(0) - geom[0].ProbLo(0);
  Real ylen = geom[0].ProbHi(1) - geom[0].ProbLo(1);
  Real zlen = geom[0].ProbHi(2) - geom[0].ProbLo(2);

  if (amrex::Math::abs(delp_x) > tolerance){
    gp0[0] = -delp_x / xlen;
  }

  if (amrex::Math::abs(delp_y) > tolerance){
    gp0[1] = -delp_y / ylen;
  }

  if (amrex::Math::abs(delp_z) > tolerance){
    gp0[2] = -delp_z / zlen;
  }

  return;
}


