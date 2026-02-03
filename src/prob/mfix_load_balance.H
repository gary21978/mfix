#ifndef _MFIX_LOAD_BALANCE_H_
#define _MFIX_LOAD_BALANCE_H_

#include <AMReX_DistributionMapping.H>
#include <AMReX_LayoutData.H>

class LoadBalance {

  public:

    LoadBalance ( int const a_max_level,
                  int const a_solve_fluid,
                  int const a_solve_solids);

    ~LoadBalance () = default;

    enum struct Type      { SingleGrid, DualGrid };
    enum struct Strategy  { SFC, KnapSack };
    enum struct Weighting { CellCount, ParticleCount, ParticleRunTime };

    int SingleGrid () const
    { return (m_type == Type::SingleGrid); }

    bool DualGrid () const
    { return (m_type == Type::DualGrid); }

    bool SFC () const
    { return (m_strategy == Strategy::SFC); }

    bool KnapSack () const
    { return (m_strategy == Strategy::KnapSack); }

    bool WeightByCellCount () const
    { return (m_weighting == Weighting::CellCount); }

    bool WeightByParticleCount () const
    { return (m_weighting == Weighting::ParticleCount); }

    bool WeightByParticleRunTime () const
    { return (m_weighting == Weighting::ParticleRunTime); }

    bool UpdateFluid () const
    { return (SingleGrid() || (m_dual_grid_load_balance_fluid > 0)); }

    void UpdateWeights ( int const a_lev, int const a_gid,
                         amrex::Real const a_weight )
    { //AMREX_ASSERT ( a_lev < m_weights.size() );
      //AMREX_ASSERT ( a_index < m_weights[a_lev].size() );
      (*m_weights[a_lev])[a_gid] += a_weight;
    }

    void ResetWeights ( int const a_lev, const amrex::BoxArray & a_grids,
                        const amrex::DistributionMapping & a_dmap);

    amrex::DistributionMapping
    MakeDistMapping ( int const a_lev );

    amrex::DistributionMapping
    MakeDistMappingFromGrids ( const amrex::BoxArray & a_grids,
                               const amrex::DistributionMapping & a_dmap);

  private:

    Type m_type = Type::SingleGrid;
    Strategy m_strategy = Strategy::KnapSack;
    Weighting m_weighting = Weighting::ParticleRunTime;

    int m_knapsack_nmax = 128;
    int m_dual_grid_load_balance_fluid = -1;

    amrex::Vector< std::unique_ptr< amrex::LayoutData<amrex::Real>>> m_weights;
    amrex::Vector< amrex::BoxArray > m_grids;

};

#endif
