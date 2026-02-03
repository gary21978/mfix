#ifndef _MFIX_PC_H_
#define _MFIX_PC_H_

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_FabArray.H>
#include <AMReX_EBSupport.H>
#include <AMReX_EBMultiFabUtil.H>

#include <mfix_ic.H>
#include <mfix_bc.H>
#include <mfix_pic.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_porous_media.H>
#include <mfix_load_balance.H>

///
/// We add a particle to the neighbor list if its center is within 3 particle radii
///
struct MFIXCheckPair
{
  amrex::Real neighborhood;
  MFIXCheckPair (amrex::Real neighborhood_) : neighborhood(neighborhood_) {}

  template <int NSR, int NSI, int NAR, int NAI>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  bool operator() (const amrex::ConstParticleTileData<amrex::Particle<NSR, NSI>, NAR, NAI>& ptile, int i, int j) const
  {
    const auto p1 = ptile.getSuperParticle(i);
    const auto p2 = ptile.getSuperParticle(j);

    if (i < j)
      if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
        return true;

    return false;
  }
};


struct MFIXCheckFullPair
{
  amrex::Real neighborhood;
  MFIXCheckFullPair (amrex::Real neighborhood_) : neighborhood(neighborhood_) {}

  template <int NSR, int NSI, int NAR, int NAI>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  bool operator() (const amrex::ConstParticleTileData<amrex::Particle<NSR, NSI>, NAR, NAI>& ptile, int i, int j) const
  {
    const auto p1 = ptile.getSuperParticle(i);
    const auto p2 = ptile.getSuperParticle(j);

    if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
      return true;

    return false;
  }
};


///
/// Neighbor search with multi-grid for polydisperse flows
///
struct MFIXCheckPolyPair
{
    int typeInd;
    int num_bin_types;
    amrex::Real* neighborhood;
    MFIXCheckPolyPair (int typeInd_, int num_bin_types_, amrex::Real* neighborhood_)
        : typeInd(typeInd_), num_bin_types(num_bin_types_), neighborhood(neighborhood_) {}

  template <int NSR, int NSI, int NAR, int NAI>
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  bool operator() (const amrex::ConstParticleTileData<amrex::Particle<NSR, NSI>, NAR, NAI>& ptile,
                   int i, int j, int type, bool ghost_i, bool ghost_j) const
  {
    const auto p1 = ptile.getSuperParticle(i);
    const auto p2 = ptile.getSuperParticle(j);

    int type_i = p1.idata(typeInd);
    int type_j = p2.idata(typeInd);

    if (!(ghost_i && ghost_j)) {                    // One particle must be a real particle
        if ((type_i == type) || (type_j == type)) { // One of the particles matches the type being searched
            if ((i < j) || (type_i < type_j)) {     // Smaller index or type particle finds the neighbor
                // Get index for NEIGHBORHOOD from analytical partial sums
                int imin = amrex::min(type_i,type_j);
                int imax = amrex::max(type_i,type_j);
                int k    = imin - 1;
                // \Sum_{i=N-k}^{N} (i)
                int indx = num_bin_types*(k+1) - (k*(k+1))/2;
                // \Sum_{i=1}^{l} (i)
                indx += imax - imin;
                if ( (p1.pos() - p2.pos()).radSquared() <= neighborhood[indx] ) return true;
            }
        }
    }

    return false;
  }
};


class MFIXParticleContainer final
    : public amrex::NeighborParticleContainer<0,0,SoArealData::count,SoAintData::count>
{
  public:

    using MFIXParIter = amrex::ParIter<0,0,SoArealData::count,SoAintData::count>;
    using MFIXParConstIter = amrex::ParConstIter<0,0,SoArealData::count,SoAintData::count>;

    MFIXParticleContainer (amrex::AmrCore* a_amrcore,
                           MFIXInitialConditions& a_initial_conditions,
                           MFIXBoundaryConditions& a_boundary_conditions,
                           MFIXSolidsPhase& a_solids,
                           MFIXDEM& a_dem, MFIXPIC& a_pic,
                           MFIXFluidPhase& a_fluid,
                           MFIXReactions& a_reactions,
                           int const a_include_vm);

    MFIXParticleContainer (const amrex::Vector<amrex::Geometry>& a_geom,
                           const amrex::DistributionMapping& a_dmap,
                           const amrex::BoxArray& a_ba,
                           MFIXInitialConditions& a_initial_conditions,
                           MFIXBoundaryConditions& a_boundary_conditions,
                           MFIXSolidsPhase& a_solids,
                           MFIXDEM& a_dem,
                           MFIXPIC& a_pic,
                           MFIXFluidPhase& a_fluid,
                           MFIXReactions& a_reactions,
                           int const a_include_vm);

    void define();

    virtual ~MFIXParticleContainer () {}

    int nlev () { return this->finestLevel()+1; }

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void InitParticlesAuto (amrex::EBFArrayBoxFactory* particle_ebfactory);

    void InitParticlesRuntimeVariables (const int adv_enthalpy);

    void RemoveOutOfRange (int lev,
                           const amrex::EBFArrayBoxFactory * ebfactory,
                           const amrex::MultiFab * ls_phi,
                           int const ls_refinement,
                           int const a_is_pic2dem,
                           const MFIXPorousMedia& a_porous_media);

    void ReportParticleGenerationStats (int lev);

    void Replicate (amrex::IntVect& Nrep,
                    amrex::Geometry& geom,
                    amrex::DistributionMapping& dmap,
                    amrex::BoxArray& ba);

    void AllocData ();

    amrex::Vector<std::map<std::pair<int, int>, amrex::NeighborList<ParticleType> > > get_neighbor_list()
    { return m_neighbor_list; }

    void EvolveParticles (int const lev,
                          amrex::Real const dt,
                          amrex::RealVect& gravity,
                          amrex::EBFArrayBoxFactory * ebfactory,
                          amrex::EBFArrayBoxFactory * particle_ebfactory,
                          const amrex::MultiFab * ls_phi,
                          const int ls_refinement,
                          LoadBalance* const a_loadbalance,
                          const amrex::MultiFab* a_Teb,
                          int & nsubsteps);

    void PostEvolveParticles (int const lev);

    template <class F1>
    void EvolveParticles (int const lev,
                          amrex::Real const dt,
                          amrex::RealVect& gravity,
                          amrex::EBFArrayBoxFactory * ebfactory,
                          amrex::EBFArrayBoxFactory * particle_ebfactory,
                          const amrex::MultiFab * ls_phi,
                          const int ls_refinement,
                          LoadBalance* const a_loadbalance,
                          const amrex::MultiFab* a_Teb,
                          int & nsubsteps,
                          F1 RollingFriction);

    void mfix_pc_inflow (int const lev,
                         int const is_dem,
                         int const is_pic,
                         amrex::Real dt,
                         const int advect_enthalpy,
                         amrex::EBFArrayBoxFactory* factory);

    void MFIX_PC_InitCollisionParams ();

    void output (int estatus,
                 int finish,
                 int nstep,
                 amrex::Real dt,
                 amrex::Real time);

    const MFIXFluidPhase& get_fluid () const { return m_fluid; }

    const MFIXSolidsPhase& get_solids () const { return m_solids; }

    const MFIXDEM& get_dem () const { return m_dem; }

    const MFIXPIC& get_pic () const {return m_pic; }

    const MFIXReactions& get_reactions () const { return m_reactions; }

    ///
    /// Deposit parcel velocity to the grid
    ///
    void MFIX_PC_SolidsVelocityDeposition (int lev,
                                           amrex::Array<amrex::MultiFab*,3>& vel_s_mf,
                                           const amrex::FabArray<amrex::EBCellFlagFab>* flags);

    void PICHydroStep (int lev,
                       const bool apply_forces,
                       const bool update_parcels,
                       const bool use_taylor_approx,
                       const amrex::Real advance_vel_p,
                       amrex::Real dt,
                       amrex::RealVect& gravity,
                       amrex::Vector< amrex::Array<amrex::MultiFab*,3> >& vel_s_in,
                       amrex::MultiFab & ep_s_out,
                       amrex::Array<amrex::MultiFab*,3>& vel_s_out,
                       const amrex::MultiFab * volfrac,
                       const amrex::FabArray<amrex::EBCellFlagFab>* flags,
                       amrex::EBFArrayBoxFactory* ebfactory,
                       const int ls_refinement,
                       const amrex::MultiFab* ls_phi);

    void MFIX_PC_AdvanceParcels (amrex::Real dt,
                                 LoadBalance* const a_loadbalance);

    const runtimeRealData m_runtimeRealData;
    const runtimeIntData m_runtimeIntData;

    void setSortingBinSizes (amrex::IntVect sorting_bin) { m_sorting_bin = sorting_bin; }

    void sortNow (int const a_nstep) {
      if ( ((m_sort_int > 0) && (a_nstep % m_sort_int == 0)) ) {
        amrex::Print() << "   Sorting particles at step " << a_nstep << '\n';
        SortParticlesByBin( m_sorting_bin );
      }
    }

    void setSortInt(int a_sort_particle_int)
    { m_sort_int = a_sort_particle_int; }

    void setReduceGhostParticles(bool reduce)
    { reduceGhostParticles = reduce; }

    // return the average/max particle count on IO process
    amrex::Real particleImbalance();

    void ReadParameters();

    int UseConstraint ()  noexcept { return m_use_constraint.sum(); }

    void ImposeMean (int const a_lev);

    MFIXBoundaryConditions& bcs () { return m_boundary_conditions; }

    MFIXFluidPhase& fluid () { return m_fluid; }

    MFIXSolidsPhase& solids () { return m_solids; }

    MFIXDEM& dem () { return m_dem; }

    MFIXPIC& pic () { return m_pic; }

    MFIXReactions& reactions () { return m_reactions; }

  private:

    int m_sort_int = -1; // # steps to sort particles

#ifdef AMREX_USE_GPU
    bool reduceGhostParticles = true;
#else
    bool reduceGhostParticles = false;
#endif

    amrex::Real newton_abstol = 1.e-6;
    amrex::Real newton_reltol = 1.e-6;
    int newton_maxiter = 100;

    amrex::IntVect m_use_constraint = amrex::IntVect::TheZeroVector();
    amrex::RealVect m_constraint = amrex::RealVect::TheZeroVector();

    amrex::IntVect max_grid_size{-1};

    amrex::IntVect m_sorting_bin = amrex::IntVect::TheZeroVector();

    static void ReadStaticParameters ();

    void usr0_des ();
    void usr1_des ();

    void usr2_des (int np, ParticleTileType& particles);

    void usr3_des (int np,void*&);

    MFIXInitialConditions& m_initial_conditions;

    MFIXBoundaryConditions& m_boundary_conditions;

    MFIXFluidPhase& m_fluid;

    MFIXSolidsPhase& m_solids;

    MFIXDEM& m_dem;

    MFIXPIC& m_pic;

    MFIXReactions& m_reactions;
};

#endif
