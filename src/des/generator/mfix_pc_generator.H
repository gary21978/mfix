#ifndef MFIX_PARTICLES_GENERATOR_H_
#define MFIX_PARTICLES_GENERATOR_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_GpuContainers.H>

#include <mfix_pc.H>
#include <mfix_ic.H>
#include <mfix_dem.H>


class ParticlesGenerator {

  public:

    using ParticleType = MFIXParticleContainer::ParticleType;
    using ParticleTileType = MFIXParticleContainer::ParticleTileType;

    // Constructor
    ParticlesGenerator ( const amrex::RealVect& plo,
                         const amrex::RealVect& dx,
                         amrex::RealBox const* const a_p_regions,
                         int const a_regions_nb,
                         int const a_allow_overlap,
                         const amrex::RealBox& a_ic_region,
                         const SOLIDS_t& a_ic_solid,
                         std::string const a_ic_pack_type,
                         int const a_has_granular_temperature,
                         int const a_is_dem,
                         int const a_is_pic,
                         int const a_is_cg_dem,
                         amrex::Real const a_multi_particle_volume);

    // Particle generator
    int generate (amrex::Box const& a_bx,
                  int const a_id, int const a_cpu,
                  ParticleTileType& particles);

  private:

    const amrex::RealVect& m_plo;
    const amrex::RealVect& m_dx;

    const amrex::RealBox& m_ic_region;

    const SOLIDS_t& m_ic_solid;

    int const m_is_dem;
    int const m_is_pic;
    int const m_is_cg_dem;

    amrex::Real const m_multi_particle_volume;

    int const m_is_hex_close_pack;
    int const m_is_random_pack;
    int const m_is_pseudo_random_pack;
    int m_cube_base;

    const int m_has_granular_temperature;

    amrex::RealBox const* const m_p_regions;
    const int m_regions_nb;
    const int m_allow_overlap;

  public:

    // Particle generator
    template <typename F1>
    int generate ( int const a_id, int const a_cpu,
                   ParticleTileType& particles,
                   F1 positions_generator);
};

#endif
