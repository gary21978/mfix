#ifndef MFIX_IC_PARMS_H_
#define MFIX_IC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_Vector.H>

#include <mfix_regions.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>


// Forward declaration
struct IC_t;


class MFIXInitialConditions
{
  public:
    // Constructor
    MFIXInitialConditions ()
      : m_allow_overlap(1)
      , m_particle_init_type(ParticleInitType::Undefined)
    {}

    void Initialize (const MFIXRegions& regions,
                     MFIXFluidPhase& fluid,
                     MFIXSolidsPhase& solids,
                     MFIXDEM& dem,
                     MFIXPIC& pic);

    int allow_overlap () const { return m_allow_overlap; }

    const amrex::Vector<IC_t>& ic () const { return m_ic; }
    const IC_t& ic (const int icv) const { return m_ic[icv]; }

    enum struct ParticleInitType {Undefined, Auto, AsciiFile};

    int AutoParticleInit () const {
      return (m_particle_init_type == ParticleInitType::Auto) ? 1 : 0;
    }

    // Store the number of particles generated in the IC region.
    void set_particle_count ( const int icv, const long np_in ) {
      AMREX_ASSERT( icv < m_np.size() );
      m_np[icv] = np_in;
    }

    // Return the number of particles generated in the IC region.
    long get_particle_count ( const int icv ) const {
      AMREX_ASSERT( icv < m_np.size() );
      return m_np[icv];
    }

    struct ICRankingType { enum {Undefined=-1, Inputs, Volume, Priority}; };

    class ICRanking {
      public:
        ICRanking ()
          : m_ranking_type(ICRankingType::Undefined)
        {}

        ICRanking (const int ranking_type)
          : m_ranking_type(ranking_type)
        {
          AMREX_ALWAYS_ASSERT(ranking_type == ICRankingType::Inputs ||
                              ranking_type == ICRankingType::Volume ||
                              ranking_type == ICRankingType::Priority);
        }

        void set_type (const int ranking_type)
        {
          if (ranking_type == ICRankingType::Inputs)
            m_ranking_type = ICRankingType::Inputs;
          else if (ranking_type == ICRankingType::Volume)
            m_ranking_type = ICRankingType::Volume;
          else if (ranking_type == ICRankingType::Priority)
            m_ranking_type = ICRankingType::Priority;
          else
            amrex::Abort("Error");
        }

        int get_type () const
        { return m_ranking_type; }

        bool operator() (const IC_t& left,
                         const IC_t& right) const;

      private:
        int m_ranking_type = ICRankingType::Undefined;
    };

    int has_granular_temperature ( ) const {
      for (int icv(0); icv < m_granular_temperature.size(); icv++ ) {
        if (has_granular_temperature(icv)) { return 1; }
      }
      return 0;
    }

    int has_granular_temperature ( const int icv ) const {
      AMREX_ASSERT( icv < m_granular_temperature.size() );
      return (m_granular_temperature[icv] > 0. ? 1 : 0);
    }

    amrex::Real get_granular_temperature ( const int icv ) const {
      AMREX_ASSERT( icv < m_granular_temperature.size() );
      return m_granular_temperature[icv];
    }

  private:
    int m_allow_overlap;
    ICRanking m_ic_ranking;

    amrex::Vector<IC_t> m_ic;

    ParticleInitType m_particle_init_type = ParticleInitType::Undefined;

    amrex::Vector<long> m_np;

    amrex::Vector<amrex::Real> m_granular_temperature;
};


struct IC_t {

  // Constructor
  IC_t (MFIXFluidPhase* fluid_in)
    : inputs_order(std::numeric_limits<int>::max())
    , volume(std::numeric_limits<amrex::Real>::max())
    , priority(std::numeric_limits<int>::max())
    , fluid(fluid_in)
  {}

  int inputs_order;
  amrex::Real volume;
  int priority;

  const amrex::RealBox* region;

  FLUID_t fluid;

  amrex::Vector<SOLIDS_t> solids;

  std::string packing;

  const SOLIDS_t* get_solid (const int phase) const
  {
    for (int n(0); n < solids.size(); ++n)
      if (phase == solids[n].phase)
        return &(solids[n]);

    amrex::Abort("Error: solid type not found");
    return nullptr;
  }
};


class MFIXICRegions
{
  public:
    // Constructor
    MFIXICRegions () = default;

    bool intersect (const amrex::RealBox& realbox) const
    {
      for (const auto& m_region: m_regions)
        if (m_region.intersects(realbox))
          return true;
      return false;
    }

    int add (const amrex::RealBox& realbox,
             const int allow_overlap = 1)
    {
      if (allow_overlap || (!this->intersect(realbox))) {
        m_regions.push_back(realbox);
        return 1; // success
      }
      return 0; // fail
    }

    int size () const
    { return m_regions.size(); }

    amrex::Gpu::HostVector<amrex::RealBox>::iterator begin ()
    { return m_regions.begin(); }

    amrex::Gpu::HostVector<amrex::RealBox>::iterator end ()
    { return m_regions.end(); }

    amrex::Gpu::HostVector<amrex::RealBox>::const_iterator begin () const
    { return m_regions.begin(); }

    amrex::Gpu::HostVector<amrex::RealBox>::const_iterator end () const
    { return m_regions.end(); }

  private:
    amrex::Gpu::HostVector<amrex::RealBox> m_regions;
};

#endif
