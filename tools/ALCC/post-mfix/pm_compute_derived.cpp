#include <AMReX.H>
#include <AMReX_Print.H>

#include <post_mfix.H>

using namespace amrex;


void
post_mfix::
compute_grad ( int const a_lev, int const a_dir,
               std::string a_var, MultiFab* a_MF)
{
  // copy the component MF to our local MF
  MultiFab filled_MF(grids[a_lev], dmap[a_lev], 1, 2);

  fill_var_MF(a_lev, a_var, &filled_MF, 0);

  filled_MF.FillBoundary(geom[a_lev].periodicity());

  // compute gradient of var at cell center
  a_MF->setVal(0.0);

  const auto dxi = geom[a_lev].InvCellSizeArray();

  for (MFIter mfi(*a_MF, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& tbox = mfi.tilebox();

    Array4<Real const> const& var  = filled_MF.const_array(mfi);
    Array4<Real      > const& grad = a_MF->array(mfi);

     // d/dx -- 2nd order central difference
    if (a_dir == 0) {

      ParallelFor(tbox, [dxi,var,grad]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { grad(i,j,k) = dxi[0]*(var(i+1,j  ,k  ) - var(i-1,j  ,k  ))*(Real(0.5)); });

    } else if (a_dir == 1) {

      ParallelFor(tbox, [dxi,var,grad]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { grad(i,j,k) = dxi[1]*(var(i  ,j+1,k  ) - var(i  ,j-1,k  ))*(Real(0.5)); });

    } else if (a_dir == 2) {

      ParallelFor(tbox, [dxi,var,grad]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { grad(i,j,k) = dxi[2]*(var(i  ,j  ,k+1) - var(i  ,j  ,k-1))*(Real(0.5)); });

    } else {

      amrex::Abort("You should not be here!");

    }
  }// MFIter

  //VisMF::Write(*lev_grad_MF, "gradients");

}


void
post_mfix::
compute_stress ( int const a_lev, std::string a_stress, MultiFab* a_MF)
{

  // copy the component MF to our local MF
  MultiFab tmp_MF(grids[a_lev], dmap[a_lev], 1, 2);
  tmp_MF.setVal(0.0);


  if ( a_stress.compare("sigma_11") != 0 ||
       a_stress.compare("sigma_22") != 0 ||
       a_stress.compare("sigma_33") != 0 ) {

    a_MF->setVal(0.0);

    fill_var_MF(a_lev, "grad_x(Uf)", &tmp_MF, 0);

    if ( a_stress.compare("sigma_11") != 0)
    { MultiFab::Saxpy(*a_MF, 2.0, tmp_MF, 0, 0 , 1, 0);}

    MultiFab::Saxpy(*a_MF, -2.0/3.0, tmp_MF, 0, 0 , 1, 0);

    fill_var_MF(a_lev, "grad_y(Vf)", &tmp_MF, 0);

    if ( a_stress.compare("sigma_22") != 0)
    { MultiFab::Saxpy(*a_MF, 2.0, tmp_MF, 0, 0 , 1, 0);}

    MultiFab::Saxpy(*a_MF, -2.0/3.0, tmp_MF, 0, 0 , 1, 0);

    fill_var_MF(a_lev, "grad_z(Wf)", &tmp_MF, 0);

    if ( a_stress.compare("sigma_33") != 0)
    { MultiFab::Saxpy(*a_MF, 2.0, tmp_MF, 0, 0 , 1, 0);}

    MultiFab::Saxpy(*a_MF, -2.0/3.0, tmp_MF, 0, 0 , 1, 0);

  } else {

    if ( a_stress.compare("sigma_12") != 0 ||
         a_stress.compare("sigma_21") != 0 ) {

      fill_var_MF(a_lev, "grad_x(Vf)", a_MF,    0);
      fill_var_MF(a_lev, "grad_y(Uf)", &tmp_MF, 0);

    } else if ( a_stress.compare("sigma_13") != 0 ||
                a_stress.compare("sigma_31") != 0 ) {

      fill_var_MF(a_lev, "grad_x(Wf)", a_MF,    0);
      fill_var_MF(a_lev, "grad_z(Uf)", &tmp_MF, 0);

    } else if ( a_stress.compare("sigma_23") != 0 ||
                a_stress.compare("sigma_32") != 0 ) {

      fill_var_MF(a_lev, "grad_y(Wf)", a_MF,    0);
      fill_var_MF(a_lev, "grad_z(Vf)", &tmp_MF, 0);

    }

    MultiFab::Add(*a_MF, tmp_MF, 0, 0, 1, 0);

  }
}


void
post_mfix::
compute_granular_viscosity ( int const a_lev, std::string a_term, MultiFab* a_MF)
{

  a_MF->setVal(0.0);

  if ( a_term.compare("sigPp_11") != 0 ||
       a_term.compare("sigPp_22") != 0 ||
       a_term.compare("sigPp_33") != 0 ) {

    fill_var_MF(a_lev, "Theta", a_MF, 0);

  }

  MultiFab tmp_MF(grids[a_lev], dmap[a_lev], 1, 2);
  tmp_MF.setVal(0.0);

  AMREX_ALWAYS_ASSERT( a_term.find("sig") == 0 );

  std::string term = a_term.substr(3);

  fill_var_MF(a_lev, term, &tmp_MF, 0);

  MultiFab::Subtract(*a_MF, tmp_MF, 0, 0, 1, 0);

}
