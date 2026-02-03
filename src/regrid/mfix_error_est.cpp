#include <AMReX_ParmParse.H>
#include <AMReX_EBAmrUtil.H>

#include <mfix.H>
#include <mfix_reporter.H>

#include <mfix_regions.H>
#include <mfix_utils.H>

//! Tag using each EB level's volfrac. This requires that the `eb_levels` have
//! already been build.
void mfix::
ErrorEst (int a_lev, TagBoxArray& a_tags, Real a_time, int /*a_ngrow*/)
{
  BL_PROFILE("mfix::ErrorEst()");
  if (m_verbose > 0) { Print() << "ErrorErr on level " << a_lev << '\n'; }

  //Refine cut cells by default
  //___________________________________________________________________________
  // Tag all cells with volfrac \in (0, 1)
  MultiFab volfrac(grids[a_lev], dmap[a_lev], 1, 1);
  eb()->levels()[a_lev]->fillVolFrac(volfrac, geom[a_lev]);

  amrex::TagVolfrac(a_tags, volfrac);

  const auto tagval = TagBox::SET;

  //bool tag_vorticity = tag().vorticity();

  bool tag_rho = tag().rho();
  bool tag_rho_lt = tag().rho.less_than();

  bool tag_grad_rho = tag().grad_rho();
  bool tag_grad_rho_lt = tag().grad_rho.less_than();

  Real rho_err = tag().rho.value(a_lev);
  Real grad_rho_err = tag().grad_rho.value(a_lev);

  if (m_verbose > 0) {
    if ( tag().rho() ) {
      std::string op_str = (tag_rho_lt) ? " < " : " > ";
      Print() << "Tag for density" << op_str << rho_err << '\n';
    }
    if ( tag().grad_rho() ) {
      std::string op_str = (tag_rho_lt) ? " < " : " > ";
      Print() << "Tag for density gradient" << op_str << grad_rho_err << '\n';
    }
  }

  if (tag_grad_rho) {
    bcs().fillpatch(a_lev, a_time, BCFillVar::rho, leveldata().rho(), 1);
  }


  bool tag_vel = tag().vel();
  bool tag_vel_lt = tag().vel.less_than();

  bool tag_grad_vel = tag().grad_vel();
  bool tag_grad_vel_lt = tag().grad_vel.less_than();

  Real vel_err = tag().vel.value(a_lev);
  Real grad_vel_err = tag().grad_vel.value(a_lev);

  if (m_verbose > 0) {
    if ( tag().vel() ) {
      std::string op_str = (tag_vel_lt) ? " < " : " > ";
      Print() << "Tag for velocity" << op_str << vel_err << '\n';
    }
    if ( tag().grad_vel() ) {
      std::string op_str = (tag_vel_lt) ? " < " : " > ";
      Print() << "Tag for velocity gradient" << op_str << grad_vel_err << '\n';
    }
  }

  if (tag_grad_vel) {
    bcs().fillpatch(a_lev, a_time, BCFillVar::vel, leveldata().vel(), 1);
  }

  int const num_regions( tag().num_regions() );
  Gpu::AsyncArray<RealBox> d_regions(tag().h_regions.data(), tag().h_regions.size());

  AMREX_D_TERM(const Real dx = geom[a_lev].CellSize(0);,
               const Real dy = geom[a_lev].CellSize(1);,
               const Real dz = geom[a_lev].CellSize(2););

  RealBox* p_regions = (num_regions > 0) ? d_regions.data() : nullptr;


  for (MFIter mfi(*leveldata().rho(a_lev),TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    Box const& bx = mfi.tilebox();
    auto const& tag = a_tags.array(mfi);

    // Refine by density or the gradient of density
    //___________________________________________________________________________

    if (tag_rho || tag_grad_rho) {

      Array4<Real const> const& rho = leveldata().rho_const(a_lev,mfi);

      ParallelFor(bx, [tag, tagval, rho, tag_rho, tag_rho_lt, rho_err,
          tag_grad_rho, tag_grad_rho_lt, grad_rho_err]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if ( tag_rho ) {
          if ( tag_rho_lt ) { if ( rho(i,j,k) < rho_err ) { tag(i,j,k) = tagval; } }
          else /*tag_rho_gt*/ { if ( rho(i,j,k) > rho_err ) { tag(i,j,k) = tagval; } }
        }

        if (tag_grad_rho) {

          Real ax = amrex::Math::abs(rho(i+1,j,k) - rho(i,j,k));
          Real ay = amrex::Math::abs(rho(i,j+1,k) - rho(i,j,k));
          Real az = amrex::Math::abs(rho(i,j,k+1) - rho(i,j,k));

          ax = amrex::max(ax,amrex::Math::abs(rho(i,j,k) - rho(i-1,j,k)));
          ay = amrex::max(ay,amrex::Math::abs(rho(i,j,k) - rho(i,j-1,k)));
          az = amrex::max(az,amrex::Math::abs(rho(i,j,k) - rho(i,j,k-1)));

          Real dr = amrex::max(ax,ay,az);

          if ( tag_grad_rho_lt ) { if ( dr < grad_rho_err ) { tag(i,j,k) = tagval; } }
          else /*tag_grad_rho_gt*/ { if ( dr > grad_rho_err ) { tag(i,j,k) = tagval; } }
        }

      });
    } // tag rho or grad_rho


    // Refine by velocity or the gradient of velocity
    //___________________________________________________________________________

    if (tag_vel || tag_grad_vel) {

      Array4<Real const> const& vel = leveldata().vel_const(a_lev,mfi);

      ParallelFor(bx, [tag, tagval, vel, tag_vel, tag_vel_lt, vel_err,
          tag_grad_vel, tag_grad_vel_lt, grad_vel_err ]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if ( tag_vel ) {

          RealVect vec_mag(vel(i,j,k,0), vel(i,j,k,1), vel(i,j,k,2));
          Real vel_mag = vec_mag.vectorLength();

          if ( tag_vel_lt ) { if ( vel_mag < vel_err) { tag(i,j,k) = tagval; } }
          else /*tag_vel_gt*/ { if ( vel_mag > vel_err) { tag(i,j,k) = tagval; } }
        }

        if (tag_grad_vel) {

          for ( int idim(0); idim<AMREX_SPACEDIM; ++idim) {

            Real ax = amrex::Math::abs(vel(i+1,j,k,idim) - vel(i,j,k,idim));
            Real ay = amrex::Math::abs(vel(i,j+1,k,idim) - vel(i,j,k,idim));
            Real az = amrex::Math::abs(vel(i,j,k+1,idim) - vel(i,j,k,idim));

            ax = amrex::max(ax,amrex::Math::abs(vel(i,j,k,idim) - vel(i-1,j,k,idim)));
            ay = amrex::max(ay,amrex::Math::abs(vel(i,j,k,idim) - vel(i,j-1,k,idim)));
            az = amrex::max(az,amrex::Math::abs(vel(i,j,k,idim) - vel(i,j,k-1,idim)));

            Real dv = amrex::max(ax,ay,az);

            if ( tag_grad_vel_lt ) { if ( dv < grad_vel_err ) { tag(i,j,k) = tagval; } }
            else /*tag_grad_vel_gt*/ { if ( dv > grad_vel_err ) { tag(i,j,k) = tagval; } }

          } // idim
        }
      });
    } // tag vel or grad_vel


    // Refine by vorticity
    //___________________________________________________________________________

    if ( num_regions > 0 ) {

      auto const& problo = geom[a_lev].ProbLoArray();

      ParallelFor(bx, [num_regions, p_regions, problo, dx, dy, dz, tag, tagval]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        RealVect point = { problo[0] + Real(i+0.5)*dx,
                           problo[1] + Real(j+0.5)*dy,
                           problo[2] + Real(k+0.5)*dz};

        // Tag if we are inside the specified box
        for ( int rcv(0); rcv < num_regions; ++rcv) {
          if ( p_regions[rcv].contains(point) ) {
            tag(i,j,k) = tagval;
          }
        }
      });
    }



  } // MFIter
}
