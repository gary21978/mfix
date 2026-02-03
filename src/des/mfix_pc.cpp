#include <mfix_pc.H>
#include <mfix_des_K.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_reporter.H>

using namespace amrex;

MFIXParticleContainer::
MFIXParticleContainer ( AmrCore* a_amrcore,
                        MFIXInitialConditions& a_initial_conditions,
                        MFIXBoundaryConditions& a_boundary_conditions,
                        MFIXSolidsPhase& a_solids,
                        MFIXDEM& a_dem, MFIXPIC& a_pic,
                        MFIXFluidPhase& a_fluid,
                        MFIXReactions&  a_reactions,
                        int const a_include_vm)
    : NeighborParticleContainer<0,0,SoArealData::count,SoAintData::count>(a_amrcore->GetParGDB(), 1)
    , m_runtimeRealData(a_solids.nspecies(),
                        a_pic.solve(),
                        a_include_vm,
                        a_include_vm || a_solids.track_acceleration(),
                        a_solids.solve_enthalpy(),
                        a_dem.tan_history_max_contacts())
    , m_runtimeIntData(a_dem.tan_history_max_contacts())
    , m_initial_conditions(a_initial_conditions)
    , m_boundary_conditions(a_boundary_conditions)
    , m_fluid(a_fluid)
    , m_solids(a_solids)
    , m_dem(a_dem)
    , m_pic(a_pic)
    , m_reactions(a_reactions)
{
  define();
}


MFIXParticleContainer::
MFIXParticleContainer ( const Vector<Geometry>& a_geom,
                        const DistributionMapping& a_dmap,
                        const BoxArray& a_ba,
                        MFIXInitialConditions&  a_initial_conditions,
                        MFIXBoundaryConditions& a_boundary_conditions,
                        MFIXSolidsPhase& a_solids,
                        MFIXDEM& a_dem, MFIXPIC& a_pic,
                        MFIXFluidPhase& a_fluid,
                        MFIXReactions&  a_reactions,
                        int const a_include_vm)
    : NeighborParticleContainer<0,0,SoArealData::count,SoAintData::count>(a_geom[0], a_dmap, a_ba, 1)
    , m_runtimeRealData(a_solids.nspecies(),
                        a_pic.solve(),
                        a_include_vm,
                        a_include_vm || a_solids.track_acceleration(),
                        a_solids.solve_enthalpy(),
                        a_dem.tan_history_max_contacts())
    , m_runtimeIntData(a_dem.tan_history_max_contacts())
    , m_initial_conditions(a_initial_conditions)
    , m_boundary_conditions(a_boundary_conditions)
    , m_fluid(a_fluid)
    , m_solids(a_solids)
    , m_dem(a_dem)
    , m_pic(a_pic)
    , m_reactions(a_reactions)
{
  define();
}


void MFIXParticleContainer::define ()
{
    ReadStaticParameters();
    ReadParameters();

    this->SetVerbose(0);

    // turn off certain components for ghost particle communication
    setRealCommComp(0, true);   // posx
    setRealCommComp(1, true);   // posy
    setRealCommComp(2, true);   // posz

    setRealCommComp(AMREX_SPACEDIM + SoArealData::radius,       true);   // radius
    setRealCommComp(AMREX_SPACEDIM + SoArealData::velx,         true);   // velx
    setRealCommComp(AMREX_SPACEDIM + SoArealData::vely,         true);   // vely
    setRealCommComp(AMREX_SPACEDIM + SoArealData::velz,         true);   // velz
    setRealCommComp(AMREX_SPACEDIM + SoArealData::omegax,       true);   // omegax
    setRealCommComp(AMREX_SPACEDIM + SoArealData::omegay,       true);   // omegay
    setRealCommComp(AMREX_SPACEDIM + SoArealData::omegaz,       true);   // omegaz
    setRealCommComp(AMREX_SPACEDIM + SoArealData::temperature,  true);   // temperature

    setRealCommComp(AMREX_SPACEDIM + SoArealData::density,      false);  // density
    setRealCommComp(AMREX_SPACEDIM + SoArealData::statwt,       false);  // statwt
    setRealCommComp(AMREX_SPACEDIM + SoArealData::drag_coeff,   false);  // drag_coeff
    setRealCommComp(AMREX_SPACEDIM + SoArealData::vel_source_x, false);  // vel_source_x
    setRealCommComp(AMREX_SPACEDIM + SoArealData::vel_source_y, false);  // vel_source_y
    setRealCommComp(AMREX_SPACEDIM + SoArealData::vel_source_z, false);  // vel_source_z

    // Add solids nspecies components
    for (int n(0); n < m_runtimeRealData.count_no_tan_history; ++n) {
      AddRealComp(true); // Turn on comm for redistribute on ghosting
      setRealCommComp(AMREX_SPACEDIM + SoArealData::count + n, false); // turn off for ghosting
    }

    for (int n(m_runtimeRealData.count_no_tan_history);
        n < m_runtimeRealData.count; ++n) {
      AddRealComp(true); // Turn on comm for redistribute on ghosting
      setRealCommComp(AMREX_SPACEDIM + SoArealData::count + n, true);
    }

#if defined(AMREX_DEBUG) || defined(AMREX_USE_ASSERTION)
    setIntCommComp(0, true);  // id
    setIntCommComp(1, true);  // cpu
#else
    setIntCommComp(0, false); // id
    setIntCommComp(1, false); // cpu
#endif
    setIntCommComp(2, true);  // phase
    setIntCommComp(3, true);  // state
#if MFIX_POLYDISPERSE
    setIntCommComp(4, true);  // ptype
#endif

    for (int n(0); n < m_runtimeIntData.count; ++n) {
      AddIntComp(true);
      setIntCommComp((SoAintData::count+1)+n, true);
    }
}

void MFIXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}

void MFIXParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    if (!initialized)
        initialized = true;
}

void MFIXParticleContainer::ReadParameters ()
{
  {
    ParmParse pp("solids.newton_solver");

    pp.query("absolute_tol", newton_abstol);
    pp.query("relative_tol", newton_reltol);
    pp.query("max_iterations", newton_maxiter);
  }

  { ParmParse pp("particles");

    std::string constraint_in = "none";
    pp.queryAdd("constraint", constraint_in);
    constraint_in = toLower(constraint_in); // case insensitive

    if ( constraint_in == "none")  {

      ; // do nothing

    } else if ( constraint_in == "mean_velocity")  {

      m_use_constraint[0] = pp.query("constraint.mean_velocity_x", m_constraint[0]);
      m_use_constraint[1] = pp.query("constraint.mean_velocity_y", m_constraint[1]);
      m_use_constraint[2] = pp.query("constraint.mean_velocity_z", m_constraint[2]);

    } else {

      reporter::Log(reporter::Error,__FILE__, __LINE__)
          << "Unknown or invalid input:\n"
          << "    particles.constraint = " << constraint_in;
    }
  }
}


namespace {
  typedef std::pair<int, int> BidNp;

  struct PairCompare {
    bool inverse = false;

    PairCompare(const bool a_inverse=false) : inverse(a_inverse) {}

    bool operator() (const BidNp& lhs, const BidNp& rhs)
    {
      return inverse ? lhs.second > rhs.second : lhs.second < rhs.second;
    }
  };

  typedef std::priority_queue<BidNp, Vector<BidNp>, PairCompare> BidNpHeap;
}



Real MFIXParticleContainer::
particleImbalance()
{
  // # particles on this process
  long local_count = 0;
  for (MFIXParIter pti(*this, 0); pti.isValid(); ++pti) {
    local_count += static_cast<long>(pti.numParticles());
  }

  long total = TotalNumberOfParticles(/*only valid=*/true,/*local=*/false);

  // max # particles per process
  ParallelDescriptor::ReduceLongMax(local_count,
      ParallelDescriptor::IOProcessorNumber());

  return ( static_cast<Real>(total) / ParallelDescriptor::NProcs()
         / (static_cast<Real>(local_count) + 1e-10));
}
