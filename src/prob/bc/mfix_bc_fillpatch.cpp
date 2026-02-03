#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <mfix_bc.H>
#include <mfix_bc_fillpatch_K.H>

using namespace amrex;


void MFIXBoundaryConditions::
fillpatch ( int const a_lev, Real const a_time, BCFillVar a_var,
            Vector<MultiFab*> const& a_MF, int const a_nghost )
{
  fillpatch(a_lev, a_time, a_var, a_MF, a_MF[a_lev], a_nghost );
}

void MFIXBoundaryConditions::
fillpatch ( int const a_lev,  Real const a_time, BCFillVar a_var,
            Vector<MultiFab*> const& a_fill_data, MultiFab* a_lev_MF,
            int const a_nghost )
{
  BL_PROFILE("MFIXBoundaryConditions::fillpatch");

  if ( m_verbose > 0 ) {
    Print() << "On level " << a_lev << " fillpatch " << BCFillVarName[a_var] << '\n';
  }

  const int minf  = BCList::minf;
  const int cover = BCList::cover;

  IntVect nghost(a_nghost);

  BCList const& bc_list = get_bc_list();

  MultiFab* fmf = a_fill_data[a_lev];
  AMREX_ASSERT( fmf != nullptr );

  int const  ncomp( a_lev_MF->nComp() );

  const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(fmf->Factory());

  switch (a_var) {

    case (BCFillVar::vel ): {

      Vector<BCRec> bcrec = get_velocity_bcrec();

      if (a_lev == 0) {

        PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > physbc(m_geom[a_lev], bcrec,
            MFIXVelFill{minf, cover, m_max_level,
              bc_u_g().data(), bc_v_g().data(), bc_w_g().data(),
              bc_list.bc_ilo[a_lev]->array(), bc_list.bc_ihi[a_lev]->array(),
              bc_list.bc_jlo[a_lev]->array(), bc_list.bc_jhi[a_lev]->array(),
              bc_list.bc_klo[a_lev]->array(), bc_list.bc_khi[a_lev]->array()
            });

        FillPatchSingleLevel(*a_lev_MF, nghost, a_time, {fmf}, {a_time}, /*scomp=*/0,
            /*dcomp=*/0, ncomp, m_geom[a_lev], physbc, /*physbc_comp=*/0);

      } else {

        MultiFab* cmf = a_fill_data[a_lev-1];
        AMREX_ASSERT( fmf != nullptr );

        PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > cphysbc(m_geom[a_lev-1], bcrec,
            MFIXVelFill{minf, cover, m_max_level,
              bc_u_g().data(), bc_v_g().data(), bc_w_g().data(),
              bc_list.bc_ilo[a_lev-1]->array(), bc_list.bc_ihi[a_lev-1]->array(),
              bc_list.bc_jlo[a_lev-1]->array(), bc_list.bc_jhi[a_lev-1]->array(),
              bc_list.bc_klo[a_lev-1]->array(), bc_list.bc_khi[a_lev-1]->array()
            });

        PhysBCFunct<GpuBndryFuncFab<MFIXVelFill> > fphysbc(m_geom[a_lev  ], bcrec,
            MFIXVelFill{minf, cover, m_max_level,
              bc_u_g().data(), bc_v_g().data(), bc_w_g().data(),
              bc_list.bc_ilo[a_lev  ]->array(), bc_list.bc_ihi[a_lev  ]->array(),
              bc_list.bc_jlo[a_lev  ]->array(), bc_list.bc_jhi[a_lev  ]->array(),
              bc_list.bc_klo[a_lev  ]->array(), bc_list.bc_khi[a_lev  ]->array()
            });

        Interpolater* mapper = (factory.isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);

        FillPatchTwoLevels(*a_lev_MF, nghost, a_time, {cmf}, {a_time}, {fmf}, {a_time},
            /*scomp=*/0, /*dcomp=*/0, ncomp, m_geom[a_lev-1], m_geom[a_lev], cphysbc, 0,
            fphysbc, 0, m_ref_ratio[a_lev-1], mapper, bcrec, /*physbc_comp=*/0);
      }

      break;
    }

    case (BCFillVar::epf ):
    case (BCFillVar::rho ):
    case (BCFillVar::X ):
    case (BCFillVar::T ):
    case (BCFillVar::h ):
    case (BCFillVar::tracer ): {

      // We use this for species bcs
      int const stride = (ncomp == 1) ? 0 : this->count();

      Vector<BCRec> bcrec;
      Real* bc_values = nullptr;

      get_fill_bcs( a_var, ncomp, bcrec, bc_values );

      // We need to treat volume fraction differently
      Real const covered_val = (a_var == BCFillVar::epf) ? 1. : 0.;

      if ( a_lev == 0) {

        PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > physbc(m_geom[a_lev], bcrec,
            MFIXScalarFill{minf, cover, m_max_level, stride, bc_values, covered_val,
              bc_list.bc_ilo[a_lev]->array(), bc_list.bc_ihi[a_lev]->array(),
              bc_list.bc_jlo[a_lev]->array(), bc_list.bc_jhi[a_lev]->array(),
              bc_list.bc_klo[a_lev]->array(), bc_list.bc_khi[a_lev]->array()
            });

        FillPatchSingleLevel(*a_lev_MF, nghost, a_time, {fmf}, {a_time}, /*scomp=*/0,
            /*dcomp=*/0, ncomp, m_geom[a_lev], physbc, /*physbc_comp=*/0 );

      } else {

        MultiFab* cmf = a_fill_data[a_lev-1];

        PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > cphysbc(m_geom[a_lev-1], bcrec,
            MFIXScalarFill{minf, cover, m_max_level, stride, bc_values, covered_val,
              bc_list.bc_ilo[a_lev-1]->array(), bc_list.bc_ihi[a_lev-1]->array(),
              bc_list.bc_jlo[a_lev-1]->array(), bc_list.bc_jhi[a_lev-1]->array(),
              bc_list.bc_klo[a_lev-1]->array(), bc_list.bc_khi[a_lev-1]->array()
            });

        PhysBCFunct<GpuBndryFuncFab<MFIXScalarFill> > fphysbc(m_geom[a_lev], bcrec,
            MFIXScalarFill{minf, cover, m_max_level, stride, bc_values, covered_val,
              bc_list.bc_ilo[a_lev  ]->array(), bc_list.bc_ihi[a_lev  ]->array(),
              bc_list.bc_jlo[a_lev  ]->array(), bc_list.bc_jhi[a_lev  ]->array(),
              bc_list.bc_klo[a_lev  ]->array(), bc_list.bc_khi[a_lev  ]->array()
            });

        Interpolater* mapper = (factory.isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);

        FillPatchTwoLevels(*a_lev_MF, nghost, a_time, {cmf}, {a_time}, {fmf}, {a_time},
            /*scomp=*/0, /*dcomp=*/0, ncomp, m_geom[a_lev-1], m_geom[a_lev], cphysbc, 0,
            fphysbc, 0, m_ref_ratio[a_lev-1], mapper, bcrec, /*physbc_comp=*/0);

      }

      break;
    }


    default: {

      Vector<BCRec> bcrec;
      Real* bc_values = nullptr;

      get_fill_bcs( a_var, ncomp, bcrec, bc_values );

      if ( a_lev == 0 ){

        PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > physbc
            (m_geom[a_lev], bcrec, MFIXForFill{/*probtype=*/-1});

        FillPatchSingleLevel(*a_lev_MF, nghost, a_time,
            {fmf}, {a_time}, /*scomp=*/0, /*scomp=*/0, ncomp,
            m_geom[a_lev], physbc, /*physbc_comp=*/0);

      } else { // a_lev > 0

        MultiFab* cmf = a_fill_data[a_lev-1];

        PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > cphysbc
            (m_geom[a_lev-1], bcrec, MFIXForFill{/*probtype=*/-1});

        PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > fphysbc
            (m_geom[a_lev  ], bcrec, MFIXForFill{/*probtype=*/-1});

        Interpolater* mapper = &pc_interp;

        FillPatchTwoLevels(*a_lev_MF, nghost, a_time,
            {cmf}, {a_time}, {fmf}, {a_time}, /*scomp=*/0, /*dcomp=*/0,
            ncomp, m_geom[a_lev-1], m_geom[a_lev], cphysbc, 0, fphysbc, 0,
            m_ref_ratio[a_lev-1], mapper, bcrec, /*physbc_comp=*/0);
      }
    } // generic fill

  } // switch a_var
}


void MFIXBoundaryConditions::
get_fill_bcs ( BCFillVar a_var, int const a_ncomp,
               Vector<BCRec>& a_bcrec, Real*& a_bc_values )
{

  switch (a_var) {

    case (BCFillVar::epf ): {
      a_bcrec = get_volfrac_bcrec();
      a_bc_values = bc_epf().data();
      break;
    }

    case (BCFillVar::rho ): {
      a_bcrec = get_density_bcrec();
      a_bc_values = bc_rho().data();
      break;
    }

    case (BCFillVar::vel ): {
      // TODO: Velocity is managed separately.
      // TODO: Change BCs to match other scalars.
      break;
    }

    case (BCFillVar::T ): {
      a_bcrec = get_enthalpy_bcrec();
      a_bc_values = bc_T().data();
      break;
    }

    case (BCFillVar::h ): {
      a_bcrec = get_enthalpy_bcrec();
      a_bc_values = bc_h().data();
      break;
    }

    case (BCFillVar::X ): {
      a_bcrec = get_species_bcrec();
      a_bc_values = bc_Xk().data();
      break;
    }

    case (BCFillVar::tracer ): {
      a_bcrec = get_tracer_bcrec();
      a_bc_values = bc_tracer().data();
      break;
    }

    default: {

      a_bcrec.resize(a_ncomp);

      for (int idim(0); idim<a_ncomp; ++idim) {

        if ( (m_geom[0].isPeriodic(idim)) ) {

          for (auto& b : a_bcrec) b.setLo(idim, BCType::int_dir);
          for (auto& b : a_bcrec) b.setHi(idim, BCType::int_dir);

        } else {

          for (auto& b : a_bcrec) b.setLo(idim, BCType::foextrap);
          for (auto& b : a_bcrec) b.setHi(idim, BCType::foextrap);

        }
      }

      break;
    }

  } // switch
}
