#include <mfix.H>

using namespace amrex;

void mfix::
compute_tra_forces (Vector<MultiFab      *> const& tra_forces,
                    Vector<MultiFab const*> const& a_rho)
{
  // NOTE: this routine must return the force term for the update of (rho s), NOT just s.
  if (fluid.solve_tracer()) {
    for (int lev(0); lev < nlev(); ++lev) {
      for (MFIter mfi(*tra_forces[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();
        Array4<Real>       const& tra_f = tra_forces[lev]->array(mfi);
        Array4<Real const> const& rho   =    a_rho[lev]->const_array(mfi);

        int const ncomp(fluid.ntracer());

        ParallelFor(bx, ncomp, [tra_f, rho]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          // For now we don't have any external forces on the scalars
          tra_f(i,j,k,n) = 0.0;

          // Return the force term for the update of (rho s), NOT just s.
          tra_f(i,j,k,n) *= rho(i,j,k);
        });
      }
    }
  }
}

void mfix::
compute_vel_forces ( Vector<MultiFab*      > const& a_forces,
                     Vector<MultiFab const*> const& a_vel,
                     Vector<MultiFab const*> const& a_rho,
                     Vector<MultiFab const*> const& a_txfr,
                     bool const a_include_grad_p)
{
  if ( m_verbose > 1 ) {
    if ( a_include_grad_p ) {
      Print() << "\nIncluding pressure gradient in vel forces\n";
    } else {
      Print() << "\nNOT including pressure gradient in vel forces\n";
    }
  }
  for (int lev(0); lev < nlev(); ++lev) {
    compute_vel_forces_on_level (lev, *a_forces[lev], *a_vel[lev],
        *a_rho[lev], *a_txfr[lev], a_include_grad_p);
  }
}

void mfix::
compute_vel_forces_on_level ( int lev,
                              MultiFab& a_forces,
                              const MultiFab& a_vel,
                              const MultiFab& a_rho,
                              const MultiFab& a_txfr,
                              bool const a_include_grad_p)
{
  GpuArray<Real,3> l_gravity{gravity[0],gravity[1],gravity[2]};
  GpuArray<Real,3> l_gp0{gp0[0], gp0[1], gp0[2]};

  for (MFIter mfi(a_forces,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    Box const& bx = mfi.tilebox();
    Array4<Real>       const& forces =  a_forces.array(mfi);
    Array4<Real const> const&   rho  =     a_rho.const_array(mfi);
    Array4<Real const> const& gradp  = leveldata().grad_p_const(lev,mfi);

    if( !a_include_grad_p ){
      amrex::ParallelFor(bx, [forces, rho, l_gravity, gradp, l_gp0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real rhoinv = 1.0/rho(i,j,k);

        forces(i,j,k,0) = l_gravity[0]-(               l_gp0[0])*rhoinv;
        forces(i,j,k,1) = l_gravity[1]-(               l_gp0[1])*rhoinv;
        forces(i,j,k,2) = l_gravity[2]-(               l_gp0[2])*rhoinv;

      });

    } else if ( a_include_grad_p ){
#if 1
      amrex::ignore_unused(a_vel, a_txfr);

      amrex::ParallelFor(bx, [forces, rho, l_gravity, gradp, l_gp0]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real rhoinv = 1.0/rho(i,j,k);

        forces(i,j,k,0) = l_gravity[0]-(gradp(i,j,k,0)+l_gp0[0])*rhoinv;
        forces(i,j,k,1) = l_gravity[1]-(gradp(i,j,k,1)+l_gp0[1])*rhoinv;
        forces(i,j,k,2) = l_gravity[2]-(gradp(i,j,k,2)+l_gp0[2])*rhoinv;
      });
#else
      Array4<Real const> const& vel  = a_vel.const_array(mfi);
      Array4<Real const> const& txfr = a_txfr.const_array(mfi);
      Array4<Real const> const& epf  = leveldata().epf_const(lev,mfi);

      InterphaseTxfrIndexes txfr_idxs;

      const int drag_comp(txfr_idxs.drag_coeff);
      const int velx_comp(txfr_idxs.vel_src_c+0);
      const int vely_comp(txfr_idxs.vel_src_c+1);
      const int velz_comp(txfr_idxs.vel_src_c+2);

      ParallelFor(bx, [forces, rho, l_gravity, gradp, l_gp0, vel, epf, txfr,
          drag_comp, velx_comp, vely_comp, velz_comp]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        const Real rhoinv = 1.0/rho(i,j,k);
        const Real epginv = 1.0/epf(i,j,k);

        const Real beta = txfr(i,j,k,drag_comp);

        const Real drag_x = (txfr(i,j,k,velx_comp) - beta*vel(i,j,k,0))*epginv;
        const Real drag_y = (txfr(i,j,k,vely_comp) - beta*vel(i,j,k,1))*epginv;
        const Real drag_z = (txfr(i,j,k,velz_comp) - beta*vel(i,j,k,2))*epginv;

        forces(i,j,k,0) = l_gravity[0]-(gradp(i,j,k,0)+l_gp0[0]+drag_x)*rhoinv;
        forces(i,j,k,1) = l_gravity[1]-(gradp(i,j,k,1)+l_gp0[1]+drag_y)*rhoinv;
        forces(i,j,k,2) = l_gravity[2]-(gradp(i,j,k,2)+l_gp0[2]+drag_z)*rhoinv;
      });
#endif

    } else {
      Abort("Bad combo of options in compute vel forces");
    }
  }
}
