#include <AMReX_ParmParse.H>

#include <pm_diffusion.H>
#include <pm_reporter.H>

using namespace amrex;

pm_diffusion::
~pm_diffusion ()
{
  for (auto& lev_bcoef : m_bcoef) {
    lev_bcoef[0].reset( nullptr );
    lev_bcoef[1].reset( nullptr );
    lev_bcoef[2].reset( nullptr );
  }
  for (auto& lev_phi : m_phi) { lev_phi.reset( nullptr ); }
  for (auto& lev_rhs : m_rhs) { lev_rhs.reset( nullptr ); }
}

pm_diffusion::
pm_diffusion ( int const a_max_level,
               Vector<Geometry>            const & a_geom,
               Vector<DistributionMapping> const & a_dmap,
               Vector<BoxArray>            const & a_grids)
  : m_max_level(a_max_level)
  , m_geom(a_geom)
  , m_dmap(a_dmap)
  , m_grids(a_grids)
  , m_filter_type(FilterType::none)
  , m_verbose(0)
  , m_bottom_verbose(0)
  , m_max_iter(100)
  , m_bottom_max_iter(100)
  , m_max_fmg_iter(0)
  , m_linop_maxorder(2)
  , m_agglomeration(true)
  , m_agg_grid_size(-1)
  , m_consolidation(true)
  , m_semicoarsening(false)
  , m_max_coarsening_level(32)
  , m_pre_smooth_iter(2)
  , m_post_smooth_iter(2)
  , m_final_smooth_iter(8)
  , m_bottom_smooth_iter(0)
  , m_max_semicoarsening_level(32)
  , m_rtol(1.0e-11)
  , m_atol(1.0e-14)
  , m_bottom_solver("bicgstab")
  , m_hypre_namespace("hypre")
{

  { ParmParse pp("diffusion");

    pp.query("verbose", m_verbose);
    pp.query("bottom_verbose", m_bottom_verbose);

    pp.query("max_iter", m_max_iter);
    pp.query("bottom_max_iter", m_bottom_max_iter);
    pp.query("max_fmg_iter", m_max_fmg_iter);
    pp.query("linop_maxorder", m_linop_maxorder);
    pp.query("agglomeration", m_agglomeration);
    pp.query("agg_grid_size", m_agg_grid_size);
    pp.query("consolidation", m_consolidation);
    pp.query("semicoarsening", m_semicoarsening);
    pp.query("max_coarsening_level", m_max_coarsening_level);
    pp.query("max_semicoarsening_level", m_max_semicoarsening_level);

    pp.query("pre_smooth_iter",    m_pre_smooth_iter);
    pp.query("post_smooth_iter",   m_post_smooth_iter);
    pp.query("final_smooth_iter",  m_final_smooth_iter);
    pp.query("bottom_smooth_iter", m_bottom_smooth_iter);

    pp.query("rtol", m_rtol);
    pp.query("atol", m_atol);
    pp.query("bottom_solver", m_bottom_solver);
    pp.query("hypre_namespace", m_hypre_namespace);
  }
}

void
pm_diffusion::
reset ( Vector<MultiFab const*> const a_alpha,
        int const a_alpha_p_comp, Real a_mean_diameter)
{
  Real diff_coeff(-1.0);

  Real filter_width(7.0);
  Real filter_sample_size(10.0);

  { ParmParse pp_filter("filter");

    std::string filter_type_str = "none";
    pp_filter.query("type", filter_type_str);

    if (amrex::toLower(filter_type_str).compare("none") == 0 ) {
      Print() << "\nParticle filtering disabled.\n";

      m_filter_type = FilterType::none;

    } else if (amrex::toLower(filter_type_str).compare("constant") == 0 ) {
      Print() << "\nConstant particle filtering enabled.\n";

      m_filter_type = FilterType::constant;

      if ( pp_filter.query("diff_coeff", diff_coeff)) {

        Print() << "\nSpecified constant diffusion coefficient: " << diff_coeff << "\n";

      } else {
        assert_has_mean_diameter(__FILE__, __LINE__, a_mean_diameter);
        pp_filter.query("width", filter_width);

        Real const delta = filter_width * a_mean_diameter;
        Real const sigma = delta / (2.0*std::sqrt(2.0*std::log(2.0)));

        diff_coeff = 0.5*sigma*sigma;

        Print() << "Filter width: " << filter_width << "\n";
        Print() << "Computed constant diffusion coefficient: " << diff_coeff << "\n";
      }

    } else if (amrex::toLower(filter_type_str).compare("variable") == 0 ) {
      Print() << "\nVariable particle filtering enabled.\n";
      assert_has_mean_diameter(__FILE__, __LINE__, a_mean_diameter);

      m_filter_type = FilterType::variable;

      pp_filter.query("width", filter_width);
      Print() << "Filter width: " << filter_width << "\n";

      pp_filter.query("sample_size", filter_sample_size);
      Print() << "Filter sample size: " << filter_sample_size << "\n";

    } else {

      post_mfix::Error(__FILE__, __LINE__) << " Error:"
        << " Invalid filter type " << filter_type_str << "\n"
        << " Please correct the input deck.";

    }
  }

  // If we aren't filtering deposited data, there is no need to
  // setup the MLMG solver.
  if ( no_filter() ) { return; }

  m_bcoef.resize(get_nlev());
  m_phi.resize(get_nlev());
  m_rhs.resize(get_nlev());

  amrex::Vector<amrex::MultiFab*> cc_bcoef(get_nlev());

  for( int lev(0); lev < get_nlev(); ++lev) {

    // cell-centered diffusion coefficient
    cc_bcoef[lev] = new MultiFab(m_grids[lev], m_dmap[lev], 1, 1);
    cc_bcoef[lev]->setVal(0.);

    // face-centered diffusion coefficient
    for(int dir = 0; dir < AMREX_SPACEDIM; dir++) {

      BoxArray edge_ba = m_grids[lev];
      edge_ba.surroundingNodes(dir);

      // one component, one ghost layer
      m_bcoef[lev][dir] = std::make_unique<MultiFab>(edge_ba, m_dmap[lev], 1, 1);
      m_bcoef[lev][dir]->setVal(0.);
    }

    // one component, one ghost layer
    m_phi[lev] = std::make_unique<MultiFab>(m_grids[lev], m_dmap[lev], 1, 1);
    m_phi[lev]->setVal(0.);

    // one component, no ghost layer
    m_rhs[lev] = std::make_unique<MultiFab>(m_grids[lev], m_dmap[lev], 1, 0);
    m_rhs[lev]->setVal(0.);

  }


  for (int lev(0); lev<get_nlev(); ++lev) {

    Real min_length = std::numeric_limits<Real>::max();

    for (int idir(0); idir < AMREX_SPACEDIM; ++idir)
    { min_length = amrex::min(min_length, m_geom[lev].ProbLength(idir)); }

    if ( m_filter_type == FilterType::constant ) {

      cc_bcoef[lev]->setVal(diff_coeff);

    } else if ( m_filter_type == FilterType::variable ) {

      Real const third = 1.0/3.0;
      Real const inv_2sqrt_2ln2 = 1.0 / (2.0*std::sqrt(2.0*std::log(2.0)));

      Real const dp = a_mean_diameter;
      Real const filter_sample_size_dp3 = filter_sample_size*(dp*dp*dp);

      // Reduce max operation for c_cfl
      ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpSum> reduce_op;
      ReduceData<Real, Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for (MFIter mfi(*cc_bcoef[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& tbox = mfi.tilebox();

        Array4<Real      > const& bcoef = cc_bcoef[lev]->array(mfi);
        Array4<Real const> const& alpha = a_alpha[lev]->const_array(mfi, a_alpha_p_comp);

        reduce_op.eval(tbox, reduce_data,
        [bcoef, alpha, third, filter_sample_size_dp3, inv_2sqrt_2ln2, min_length]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {

          Real const eps = alpha(i,j,k) + std::numeric_limits<Real>::epsilon();
          Real const sigma = std::pow( filter_sample_size_dp3/eps, third) * inv_2sqrt_2ln2;

          bcoef(i,j,k) = amrex::min( 0.5*sigma*sigma, min_length);

          return {bcoef(i,j,k), bcoef(i,j,k), bcoef(i,j,k)};

        });

      } // MFIter

      ReduceTuple host_tuple = reduce_data.value();
      Real min_b = amrex::get<0>(host_tuple);
      ParallelDescriptor::ReduceRealMin({min_b});

      Real max_b = amrex::get<1>(host_tuple);
      ParallelDescriptor::ReduceRealMax({max_b});

      Real avg_b = amrex::get<2>(host_tuple);
      ParallelDescriptor::ReduceRealSum({avg_b});

      avg_b /= static_cast<Real>((m_geom[lev].Domain()).numPts());

      Print() << "\nVariable diffusion coefficient:\n"
        << "  min: " << min_b << "\n"
        << "  max: " << max_b << "\n"
        << "  avg: " << avg_b << "\n";

    } // FilterType

    cc_bcoef[lev]->FillBoundary(m_geom[lev].periodicity());

    amrex::average_cellcenter_to_face(GetArrOfPtrs(m_bcoef[lev]), (*cc_bcoef[lev]), m_geom[lev]);

    delete cc_bcoef[lev];

  } // lev


  LPInfo info;
  info.setAgglomeration(m_agglomeration);
  info.setConsolidation(m_consolidation);
  info.setMaxCoarseningLevel(m_max_coarsening_level);

  info.setSemicoarsening(m_semicoarsening);
  info.setMaxSemicoarseningLevel(m_max_semicoarsening_level);

  m_mlabec = std::make_unique<MLABecLaplacian>(m_geom, m_grids, m_dmap, info);

  m_mlabec->setMaxOrder(m_linop_maxorder);

  // This is a 3D problem with periodic BC
  m_mlabec->setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                      LinOpBCType::Periodic,
                                      LinOpBCType::Periodic)},
                        {AMREX_D_DECL(LinOpBCType::Periodic,
                                      LinOpBCType::Periodic,
                                      LinOpBCType::Periodic)});

  m_mlabec->setScalars(1.0, 1.0);

  for (int lev(0); lev<get_nlev(); ++lev) {
    m_mlabec->setLevelBC(lev, nullptr);
    m_mlabec->setACoeffs(lev, 1.0);
    m_mlabec->setBCoeffs(lev, GetArrOfConstPtrs(m_bcoef[lev]));
  }
}

void
pm_diffusion::
assert_has_mean_diameter (std::string const a_file,
                          int const a_line,
                          Real const a_mean_diameter) {

  if ( a_mean_diameter < 0.) {

    post_mfix::Error(a_file, a_line) << " Error:"
      << " The selected filter requires a mean particle diameter. One was not provide"
      << " in the inputs file (particle.mean_diameter) nor computed from the plot file.\n"
      << " Please correct the input deck.";
  }

}


void
pm_diffusion::
smooth ( const Vector< MultiFab* >& a_MF,
         int const a_MF_srccomp,
         int const a_MF_numcomp,
         Vector< std::string > a_variables)
{

  if ( m_filter_type != FilterType::none) {

    if (a_variables.size()) { Print() << "\n"; }

    for (int comp(a_MF_srccomp); comp < a_MF_numcomp; ++comp) {

      if (a_variables.size()) { Print() << "  smoothing " << a_variables[comp] << "\n"; }

      for (int lev(0); lev<get_nlev(); ++lev) {

        // copy source component of a_MF into phi
        m_rhs[lev]->setVal(0.0);
        { int const srccomp = comp;
          int const dstcomp = 0;
          int const numcomp = 1;
          int const nghost  = 0;
          MultiFab::Copy(*m_rhs[lev], *a_MF[lev], srccomp, dstcomp, numcomp, nghost);
        }

        // copy source component of a_MF into phi
        m_phi[lev]->setVal(0.0);
        { int const srccomp = comp;
          int const dstcomp = 0;
          int const numcomp = 1;
          int const nghost  = 1;
          MultiFab::Copy(*m_phi[lev], *a_MF[lev], srccomp, dstcomp, numcomp, nghost);
        }
      }

      MLMG mlmg(*m_mlabec);

      // The default bottom solver is BiCG
      if(m_bottom_solver == "bicgstab")
      { mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab); }
      else if(m_bottom_solver == "smoother")
      { mlmg.setBottomSolver(MLMG::BottomSolver::smoother); }
      else if(m_bottom_solver == "cg")
      { mlmg.setBottomSolver(MLMG::BottomSolver::cg); }
      else if(m_bottom_solver == "bicgcg")
      { mlmg.setBottomSolver(MLMG::BottomSolver::bicgcg); }
      else if(m_bottom_solver == "cgbicg")
      { mlmg.setBottomSolver(MLMG::BottomSolver::cgbicg); }
      else if(m_bottom_solver == "hypre") {
#ifdef AMREX_USE_HYPRE
        mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
        mlmg.setHypreOptionsNamespace(m_hypre_namespace);
        mlmg.setHypreInterface(amrex::Hypre::Interface::ij);
#else
        post_mfix::Error(__FILE__, __LINE__) << " Error:"
          << " diffusion.bottom_solve = hypre but post_mfix was not built with hypre support\n"
          << " Please correct the input deck.";
#endif
      }
      else {
        post_mfix::Error(__FILE__, __LINE__) << " Error:"
          << " Unknown diffusion bottom solve " << m_bottom_solver << "\n"
          << " Please correct the input deck.";
      }

      mlmg.setMaxIter(m_max_iter);
      mlmg.setMaxFmgIter(m_max_fmg_iter);
      mlmg.setVerbose(m_verbose);
      mlmg.setBottomVerbose(m_bottom_verbose);

      mlmg.setFinalFillBC(false);

      mlmg.setPreSmooth(m_pre_smooth_iter);
      mlmg.setPostSmooth(m_post_smooth_iter);
      mlmg.setFinalSmooth(m_final_smooth_iter);
      mlmg.setBottomSmooth(m_bottom_smooth_iter);

      mlmg.solve(GetVecOfPtrs(m_phi), GetVecOfConstPtrs(m_rhs), m_rtol, m_atol);

      for(int lev(0); lev<get_nlev(); ++lev) {

        m_phi[lev]->FillBoundary(m_geom[lev].periodicity());

        // copy phi into the destination component of a_MF
        { int const srccomp = 0;
          int const dstcomp = comp;
          int const numcomp = 1;
          int const nghost  = 1;
          MultiFab::Copy(*a_MF[lev], *m_phi[lev], srccomp, dstcomp, numcomp, nghost);
        }
      }
    } // diffusing a_MF source component
  }
}
