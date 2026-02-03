#include <mfix_pc.H>

using ParticleType = MFIXParticleContainer::ParticleType;
using ParticleTileType = MFIXParticleContainer::ParticleTileType;

void MFIXParticleContainer::usr2_des (int,
                                      ParticleTileType& particles)
{
// Purpose: This routine is called within the discrete phase time loop
// after the source terms are applied and the time step updated.

    auto& aos = particles.GetArrayOfStructs();
    ParticleType* pstruct = aos().dataPtr();

    auto& soa = particles.GetStructOfArrays();
    auto p_realarray = soa.realarray();

    // Move particles (Fortran ordering) 63-93 below particles 32-62 to fake a wall.
    for (unsigned i = 62; i < 94; ++i)
    {
        pstruct[i].pos(0)   = 0.0475;
        pstruct[i].pos(1)   = pstruct[i-31].pos(1);
        pstruct[i].pos(2)   = pstruct[i-31].pos(2);

        p_realarray[SoArealData::velx][i]  = 0.0;
        p_realarray[SoArealData::vely][i]  = 0.0;
        p_realarray[SoArealData::velz][i]  = 0.0;

        p_realarray[SoArealData::omegax][i]  = 0.0;
        p_realarray[SoArealData::omegay][i]  = 0.0;
        p_realarray[SoArealData::omegaz][i]  = 0.0;
    }
}

void MFIXParticleContainer::usr3_des (int, void*&)
{
}
