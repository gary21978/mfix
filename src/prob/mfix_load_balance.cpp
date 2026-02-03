#include <AMReX_ParmParse.H>

#include <mfix_fix_inputs.H>
#include <mfix_load_balance.H>
#include <mfix_reporter.H>

using namespace amrex;

LoadBalance::
LoadBalance ( int const a_max_level,
              int const a_solve_fluid,
              int const a_solve_solids)
{

  { FixInputs fix("Mar. 2025");

    ParmParse pp;

    // mfix.dual_grid = 0 / 1
    // --> mfix.load_balance = SingleGrid / DualGrid
    if ( fix.check_old( "mfix.dual_grid", "mfix.load_balance") ) {

      int dual_grid(0);
      pp.get("mfix.dual_grid", dual_grid);
      pp.remove("mfix.dual_grid");

      if ( dual_grid == 1 ) {
        pp.add("mfix.load_balance", std::string("DualGrid"));
      } else {
        pp.add("mfix.load_balance", std::string("SingleGrid"));
      }
    }

    fix.swap<std::string>("mfix.load_balance_type",
                          "mfix.load_balance.strategy", 1);

    fix.swap<std::string>("mfix.knapsack_weight_type",
                          "mfix.load_balance.weighting", 1);

    fix.swap<std::string>("mfix.load_balance_fluid",
                          "mfix.load_balance.DualGrid.fluid", 1);
  } // END cleanup old inputs


  ParmParse ppMFIX("mfix");

  { std::string type_str;
   if ( ppMFIX.query("load_balance", type_str) ) {

     if ( toLower(type_str).compare("singlegrid") == 0) {

       m_type = LoadBalance::Type::SingleGrid;

     } else if ( toLower(type_str).compare("dualgrid") == 0) {

       m_type = LoadBalance::Type::DualGrid;

     } else {
       reporter::Log(reporter::Error,__FILE__, __LINE__)
           << "Unknown or invalid load balance type: " << type_str << "\n"
           << "Valid options: SingleGrid and DualGrid\n"
           << "Please correct the input deck.";
     }
   }
 }


  if ( m_type == LoadBalance::Type::DualGrid ) {

    if (!a_solve_fluid) {

      reporter::Log(reporter::Warning)
          << "Dual grid disabled because fluid is not solved";

    } else if ( !a_solve_solids) {

      reporter::Log(reporter::Warning)
          << "Dual grid disabled because particles are not solved";

    } else if (ParallelDescriptor::NProcs() == 1) {
      reporter::Log(reporter::Warning)
          << "Dual grid disabled because NProcs() == 1";
    }
    m_type = Type::DualGrid;

    ppMFIX.query("load_balance.DualGrid.fluid", m_dual_grid_load_balance_fluid);

  }


  // Load balancing strategy: Defaults to SFC
  { std::string strategy_str;
    if (ppMFIX.query("load_balance.strategy", strategy_str)) {
      if ( toLower(strategy_str).compare("knapsack") == 0) {

        m_strategy = LoadBalance::Strategy::KnapSack;

      } else if ( toLower(strategy_str).compare("sfc") == 0) {

        m_strategy = LoadBalance::Strategy::SFC;

      } else {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Unknown or invalid load balance strategy: " << strategy_str << "\n"
            << "Valid options: SFC and KnapSack\n"
            << "Please correct the input deck.";
      }
    }
  }


  // Load balancing weighting
  { std::string weighting_str;
    if (ppMFIX.query("load_balance.weighting", weighting_str)) {

      if ( toLower(weighting_str).compare("particleruntime") == 0 ||
           toLower(weighting_str).compare("runtimecosts") == 0) {

        m_weighting = LoadBalance::Weighting::ParticleRunTime;

      } else if ( toLower(weighting_str).compare("particlecount") == 0 ||
                  toLower(weighting_str).compare("numparticles") == 0) {

        m_weighting = LoadBalance::Weighting::ParticleCount;

      } else if ( toLower(weighting_str).compare("cellcount") == 0) {

        m_weighting = LoadBalance::Weighting::CellCount;

      } else {

        reporter::Log(reporter::Error,__FILE__, __LINE__)
            << "Unknown or invalid load balance weighting: " << weighting_str << "\n"
            << "Valid options: CellCount, ParticleCount and ParticleRunTime\n"
            << "Please correct the input deck.";
      }

    } else {

      m_weighting = (a_solve_solids)
          ? LoadBalance::Weighting::ParticleRunTime
          : LoadBalance::Weighting::CellCount;
    }
  }

  if (m_strategy == LoadBalance::Strategy::KnapSack) {
    ppMFIX.query("knapsack_nmax", m_knapsack_nmax);
  }

  {
    Print() << "\nLoad balance setup:\n"
        << "   Type:      " << (SingleGrid() ? "SingleGrid" : "DualGrid") << "\n"
        << "   Strategy:  " << (SFC() ? "SFC" : "KnapSack") << "\n"
        << "   Weighting: " << (WeightByParticleCount()   ? "ParticleCount" :
        (WeightByParticleRunTime() ? "ParticleRunTime" : "CellCount")) << "\n\n";
  }

  m_weights.resize(a_max_level+1);
  m_grids.resize(a_max_level+1);

}


void LoadBalance::
ResetWeights ( int const a_lev, const BoxArray & a_grids,
               const DistributionMapping & a_dmap)
{
  AMREX_ALWAYS_ASSERT ( a_lev >= 0 );
  AMREX_ALWAYS_ASSERT ( a_lev < m_weights.size() );

  m_grids[a_lev] = a_grids;
  m_weights[a_lev].reset( new LayoutData<Real>(a_grids, a_dmap) );

  if ( WeightByCellCount() ) {
    for (MFIter mfi(a_grids, a_dmap); mfi.isValid(); ++mfi) {
      Box const& bx = mfi.validbox();
      (*m_weights[a_lev])[mfi] = static_cast<Real>(bx.numPts());
    }
  }
}


DistributionMapping
LoadBalance::
MakeDistMapping ( int const a_lev )
{
  AMREX_ALWAYS_ASSERT ( a_lev >= 0 );
  AMREX_ALWAYS_ASSERT ( a_lev < m_weights.size() );

  Vector<Real> weights;
  weights.resize( m_weights[a_lev]->size() );

  ParallelDescriptor::GatherLayoutDataToVector( *m_weights[a_lev], weights,
      ParallelContext::IOProcessorNumberSub());

  ParallelDescriptor::Bcast(weights.data(), weights.size(),
      ParallelContext::IOProcessorNumberSub());

  if ( KnapSack() ) {

    return DistributionMapping::makeKnapSack(weights, m_knapsack_nmax);

  } else { AMREX_ALWAYS_ASSERT ( SFC() );

    return DistributionMapping::makeSFC(weights, m_grids[a_lev], false);

  }
}


DistributionMapping
LoadBalance::
MakeDistMappingFromGrids ( const BoxArray& a_grids,
                           const DistributionMapping& a_dmap)
{
  amrex::LayoutData<amrex::Real> ld_weights(a_grids, a_dmap);

  for (MFIter mfi(a_grids, a_dmap); mfi.isValid(); ++mfi) {
    Box const& bx = mfi.validbox();
    ld_weights[mfi] = static_cast<Real>(bx.numPts());
  }

  Vector<Real> weights( a_grids.size() );

  ParallelDescriptor::GatherLayoutDataToVector( ld_weights, weights,
      ParallelContext::IOProcessorNumberSub());

  ParallelDescriptor::Bcast(weights.data(), weights.size(),
      ParallelContext::IOProcessorNumberSub());

  if ( KnapSack() ) {

    return DistributionMapping::makeKnapSack(weights, m_knapsack_nmax);

  } else { AMREX_ALWAYS_ASSERT ( SFC() );

    return DistributionMapping::makeSFC(weights, a_grids, false);

  }
}
