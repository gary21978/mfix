#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <mfix_bc.H>
#include <mfix_bc_fillpatch_K.H>

using namespace amrex;

void MFIXBoundaryConditions::
fillcoarsepatch ( int const a_lev,  Real const a_time, BCFillVar a_var,
                  Vector<MultiFab*> const& a_fill_data, MultiFab* a_lev_MF,
                  int const a_nghost )
{
  BL_PROFILE("MFIXBoundaryConditions::fillcoarsepatch");

  if ( m_verbose > 0 ) {
    Print() << "On level " << a_lev << ": fillcoarsepatch " << BCFillVarName[a_var] << '\n';
  }

  const int minf  = BCList::minf;
  const int cover = BCList::cover;

  IntVect nghost(a_nghost);

  BCList const& bc_list = get_bc_list();

  int const  ncomp( a_lev_MF->nComp() );

  const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(a_lev_MF->Factory());

  switch (a_var) {

    case (BCFillVar::vel ): {

      Vector<BCRec> bcrec = get_velocity_bcrec();

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

      InterpFromCoarseLevel( *a_lev_MF, nghost, a_time, *a_fill_data[a_lev-1],
          /*scomp=*/0, /*dcomp=*/0, ncomp,
          m_geom[a_lev-1], m_geom[a_lev], cphysbc, 0, fphysbc, 0,
          m_ref_ratio[a_lev-1], mapper, bcrec, /*physbc_comp=*/0);

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

      InterpFromCoarseLevel( *a_lev_MF, nghost, a_time, *a_fill_data[a_lev-1],
          /*scomp=*/0, /*dcomp=*/0, ncomp,
          m_geom[a_lev-1], m_geom[a_lev], cphysbc, 0, fphysbc, 0,
          m_ref_ratio[a_lev-1], mapper, bcrec, /*physbc_comp=*/0);

      break;
    }

    default: {

      Vector<BCRec> bcrec;
      Real* bc_values = nullptr;

      get_fill_bcs( a_var, ncomp, bcrec, bc_values );

      PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > cphysbc
          (m_geom[a_lev-1], bcrec, MFIXForFill{/*probtype=*/-1});

      PhysBCFunct<GpuBndryFuncFab<MFIXForFill> > fphysbc
          (m_geom[a_lev  ], bcrec, MFIXForFill{/*probtype=*/-1});

      Interpolater* mapper = &pc_interp;

      InterpFromCoarseLevel( *a_lev_MF, nghost, a_time, *a_fill_data[a_lev-1],
          /*scomp=*/0, /*dcomp=*/0, ncomp,
          m_geom[a_lev-1], m_geom[a_lev], cphysbc, 0, fphysbc, 0,
          m_ref_ratio[a_lev-1], mapper, bcrec, /*physbc_comp=*/0);

    } // default

  } // switch a_var
}
