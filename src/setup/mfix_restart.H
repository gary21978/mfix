#ifndef MFIX_RESTART_H_
#define MFIX_RESTART_H_

struct pcCheckpointPatch
{

  pcCheckpointPatch ( std::tuple<int,int> a_checkpointVersion
                    , int const a_runtimeReal_count
                    , int const a_runtimeInt_count
                    , int const a_ncomp_ep_s
                    , int const a_ncomp_energy )
  : runtimeReal_count(a_runtimeReal_count)
  , runtimeInt_count(a_runtimeInt_count)
  {

    // Known checkpoint file versions:
    std::tuple<int const, int const> v_1_3 {1,3};
    std::tuple<int const, int const> v_1_4 {1,4};
    std::tuple<int const, int const> v_1_5 {1,5};
    std::tuple<int const, int const> v_1_6 {1,6};


    // runtimeReal_count and runtimeInt_count are initialized
    // to the size of the current particle container. We need to
    // increment (down or up) to account for added / removed entries.


    // v1.3 changes
    if ( a_checkpointVersion < v_1_3 )
    { soa_has_volume = 1; }  // removed from SoArealData in v1.3

    // v1.4 changes
    if ( a_checkpointVersion < v_1_4 )
    { soa_has_mass   = 1; }  // removed from SoArealData in v1.4

    // v1.5 changes
    if ( a_checkpointVersion < v_1_5 ) {
      soa_has_ep_s   = 1;  // removed from SoArealData in v1.5
      soa_has_cp_s   = 1;  // removed from SoArealData in v1.5

      rrd_has_ep_s   = 0;  // added to runtimeRealData in v1.5
      runtimeReal_count -= a_ncomp_ep_s;
    }

    // v1.6 changes
    if ( a_checkpointVersion < v_1_6 ) {
      soa_has_convection = 1;  // removed from SoArealData in v1.6

      rrd_has_energy_src = 0;  // added to runtimeRealData in v1.6
      runtimeReal_count -= a_ncomp_energy;
    }
  }

  // Changes in v1.3
  int soa_has_volume = 0;     ///< Flag: SoArealData contains volume

  // Changes in v1.4
  int soa_has_mass   = 0;     ///< Flag: SoArealData contains mass

  // Changes in v1.5
  int soa_has_ep_s   = 0;     ///< Flag: SoArealData contains mass
  int soa_has_cp_s   = 0;     ///< Flag: SoArealData contains mass
  int rrd_has_ep_s   = 1;     ///< Flag: runtimeReaData contains eps

  // Changes v1.6
  int soa_has_convection    = 0; ///< Flag: SoArealData contains convection
  int rrd_has_energy_src    = 1; ///< Flag: runtimeReaData contains energy_src

  int runtimeReal_count;    ///< Number of runtimeRealData in checkpoint
  int runtimeInt_count;     ///< Number of runtimeIntData in checkpoint


  // The SoAintData is included in the base struct because it is the same
  // for all versions (no changes).
  enum SoAintData { phase = 0, state,
#if MFIX_POLYDISPERSE
    ptype,
#endif
    count
  };

};


struct pcPatch_v_1_2 : public pcCheckpointPatch
{
  std::tuple<int const,int const> patchVersion {1,2}; ///< Version {Major,Minor}

  pcPatch_v_1_2 ( std::tuple<int,int> a_checkpointVersion
                , int const a_runtimeReal_count
                , int const a_runtimeInt_count
                , int const a_ncomp_ep_s
                , int const a_ncomp_energy )
    : pcCheckpointPatch ( a_checkpointVersion, a_runtimeReal_count, a_runtimeInt_count,
                          a_ncomp_ep_s, a_ncomp_energy)
  { }

  enum SoArealData { radius = 0, volume, mass, density, ep_s, velx, vely, velz,
         omegax, omegay, omegaz, statwt, drag_coeff, vel_source_x,
         vel_source_y, vel_source_z, cp_s, temperature, convection, count };
};


struct pcPatch_v_1_3 : public pcCheckpointPatch
{
  std::tuple<int const,int const> patchVersion {1,3}; ///< Version {Major,Minor}

  pcPatch_v_1_3 ( std::tuple<int,int> a_checkpointVersion
                , int const a_runtimeReal_count
                , int const a_runtimeInt_count
                , int const a_ncomp_ep_s
                , int const a_ncomp_energy )
    : pcCheckpointPatch ( a_checkpointVersion, a_runtimeReal_count, a_runtimeInt_count,
                          a_ncomp_ep_s, a_ncomp_energy)
  { }


  enum SoArealData {radius = 0, mass, density, ep_s, velx, vely, velz, omegax, omegay, omegaz,
    statwt, drag_coeff, vel_source_x, vel_source_y, vel_source_z, cp_s, temperature,
    convection, count};
};



struct pcPatch_v_1_4 : public pcCheckpointPatch
{
  std::tuple<int const,int const> patchVersion {1,4}; ///< Version {Major,Minor}

  pcPatch_v_1_4 ( std::tuple<int,int> a_checkpointVersion
                , int const a_runtimeReal_count
                , int const a_runtimeInt_count
                , int const a_ncomp_ep_s
                , int const a_ncomp_energy )
    : pcCheckpointPatch ( a_checkpointVersion, a_runtimeReal_count, a_runtimeInt_count,
                          a_ncomp_ep_s, a_ncomp_energy)
  { }

  enum SoArealData {radius = 0, density, ep_s, velx, vely, velz, omegax, omegay, omegaz,
      statwt, drag_coeff, vel_source_x, vel_source_y, vel_source_z, cp_s, temperature,
      convection, count};
};


struct pcPatch_v_1_5 : public pcCheckpointPatch
{
  std::tuple<int const,int const> patchVersion {1,5}; ///< Version {Major,Minor}

  pcPatch_v_1_5 ( std::tuple<int,int> a_checkpointVersion
                , int const a_runtimeReal_count
                , int const a_runtimeInt_count
                , int const a_ncomp_ep_s
                , int const a_ncomp_energy )
    : pcCheckpointPatch ( a_checkpointVersion, a_runtimeReal_count, a_runtimeInt_count,
                          a_ncomp_ep_s, a_ncomp_energy)
  { }

  enum struct SoArealData { radius = 0, density
      , velx, vely, velz, omegax, omegay, omegaz, statwt
      , drag_coeff, vel_source_x, vel_source_y, vel_source_z
      , temperature, convection, count
      , ep_s = 0
 };


};
#endif
