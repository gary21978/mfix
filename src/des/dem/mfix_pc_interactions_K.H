#ifndef _MFIX_PC_INTERACTIONS_K_H_
#define _MFIX_PC_INTERACTIONS_K_H_

#include <mfix_des_K.H>
#include <mfix_des_cond_K.H>
#include <mfix_des_coll_K.H>
#include <mfix_pc.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_reactions.H>
#include <mfix_bc.H>
#include <mfix_solvers.H>
#include <mfix_calc_cell.H>


template <class F1>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
compute_pw_collisions (const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                       amrex::RealVect& local_ft, amrex::RealVect& overlap_t, amrex::RealVect& tow_roll,
                       bool* touch, const int& i,
                       const amrex::Real& rp,
                       const amrex::RealVect& p1omega,
                       const amrex::RealVect& local_fn, const amrex::RealVect& vrel_t, const amrex::RealVect& normal,
                       const amrex::Real& dt, const amrex::Real& mew_w,
                       const amrex::Real& kt_des_w, const amrex::Real& etat_des_w, F1 RollingFriction,
                       const int& max_contacts_tan_history, const int& pft_neighbor_idx, const int& cpu_id_idx,
                       const int& pid_idx, const int& wall_pid_idx)
{
   // variables for tangential history
   amrex::RealVect pft_neighbor(0.);
   int contact_idx = -1;

   // calculate the tangential overlap
   for(int n(0); n < max_contacts_tan_history; ++n) {
      if (ptile_data.m_runtime_idata[pid_idx + n][i] == wall_pid_idx) {
         contact_idx = n;
         pft_neighbor[0] = ptile_data.m_runtime_rdata[pft_neighbor_idx + n][i];
         pft_neighbor[1] = ptile_data.m_runtime_rdata[pft_neighbor_idx + n
            + max_contacts_tan_history][i];
         pft_neighbor[2] = ptile_data.m_runtime_rdata[pft_neighbor_idx + n
            + 2*max_contacts_tan_history][i];
         break;
      }
   }

   // calculate the tangential displacement, adjusting old spring to new frame of reference
   overlap_t[0] = dt*vrel_t[0] + pft_neighbor[0] - normal[0]*dot_product(normal, pft_neighbor);
   overlap_t[1] = dt*vrel_t[1] + pft_neighbor[1] - normal[1]*dot_product(normal, pft_neighbor);
   overlap_t[2] = dt*vrel_t[2] + pft_neighbor[2] - normal[2]*dot_product(normal, pft_neighbor);

   // static friction (non-Coulomb) tangential force
   amrex::Real fnmd = sqrt(dot_product(local_fn, local_fn));  // normal force mag
   local_ft[0] = -kt_des_w*overlap_t[0] - etat_des_w*vrel_t[0];
   local_ft[1] = -kt_des_w*overlap_t[1] - etat_des_w*vrel_t[1];
   local_ft[2] = -kt_des_w*overlap_t[2] - etat_des_w*vrel_t[2];
   amrex::Real ftsmd = sqrt(dot_product(local_ft, local_ft));  //f_t static mag
   amrex::Real ftdmd = mew_w*fnmd;                             //f_t dynamic mag

   if (ftsmd > ftdmd) {
      amrex::RealVect tangent(0.);
      tangent[0] = -local_ft[0]/ftsmd;
      tangent[1] = -local_ft[1]/ftsmd;
      tangent[2] = -local_ft[2]/ftsmd;

      // dynamic  tangential friction
      local_ft[0] = -ftdmd * tangent[0];
      local_ft[1] = -ftdmd * tangent[1];
      local_ft[2] = -ftdmd * tangent[2];

      // limit tangential spring so that static friction = dynamic friction
      overlap_t[0] = (ftdmd*tangent[0] - etat_des_w*vrel_t[0])/kt_des_w;
      overlap_t[1] = (ftdmd*tangent[1] - etat_des_w*vrel_t[1])/kt_des_w;
      overlap_t[2] = (ftdmd*tangent[2] - etat_des_w*vrel_t[2])/kt_des_w;
   }

   tow_roll = RollingFriction(rp, -1.0, fnmd, p1omega, p1omega);

   // update tangential history contact data
   if (contact_idx == -1) {
      for (int n(0); n < max_contacts_tan_history; ++n) {
         if ((ptile_data.m_runtime_idata[cpu_id_idx + n][i] == -1) &&
             (ptile_data.m_runtime_idata[pid_idx + n][i] == -1)) {
            contact_idx = n;
            break;
         }
      }
   }
   AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (contact_idx != -1), "Particle neighbor contacts exceeded");
   ptile_data.m_runtime_idata[cpu_id_idx + contact_idx][i] = -1;
   ptile_data.m_runtime_idata[pid_idx + contact_idx][i] = wall_pid_idx;
   touch[contact_idx] = true;
   ptile_data.m_runtime_rdata[pft_neighbor_idx + contact_idx][i] = overlap_t[0];
   ptile_data.m_runtime_rdata[pft_neighbor_idx + contact_idx
      + max_contacts_tan_history][i] = overlap_t[1];
   ptile_data.m_runtime_rdata[pft_neighbor_idx + contact_idx
      + 2*max_contacts_tan_history][i] = overlap_t[2];
}

template <class F1>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
compute_pp_collisions (const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                       amrex::RealVect& local_ft, amrex::RealVect& overlap_t, amrex::RealVect& tow_roll,
                       bool* touch,
                       const int& i, const int& p2_cpu, const int& p2_id,
                       const amrex::Real& p1radius, const amrex::Real& p2radius,
                       const amrex::RealVect& p1omega, const amrex::RealVect& p2omega,
                       const amrex::RealVect& local_fn, const amrex::RealVect& vrel_t, const amrex::RealVect& normal,
                       const amrex::Real& dt, const amrex::Real& mew, const amrex::Real& kt_des,
                       const amrex::Real& etat_des, F1 RollingFriction,
                       const int& max_contacts_tan_history, const int& pft_neighbor_idx, const int& cpu_id_idx,
                       const int& pid_idx)
{
   // variables for tangential history
   amrex::RealVect pft_neighbor(0.);
   int contact_idx = -1;

   // calculate the tangential overlap
   for(int n(0); n < max_contacts_tan_history; ++n) {
      if ((ptile_data.m_runtime_idata[cpu_id_idx + n][i] == p2_cpu) &&
          (ptile_data.m_runtime_idata[pid_idx + n][i] == p2_id)) {
         contact_idx = n;
         pft_neighbor[0] = ptile_data.m_runtime_rdata[pft_neighbor_idx + n][i];
         pft_neighbor[1] = ptile_data.m_runtime_rdata[pft_neighbor_idx + n
            + max_contacts_tan_history][i];
         pft_neighbor[2] = ptile_data.m_runtime_rdata[pft_neighbor_idx + n
            + 2*max_contacts_tan_history][i];
         break;
      }
   }

   // calculate the tangential displacement, adjusting old spring to new frame of reference
   overlap_t[0] = dt*vrel_t[0] + pft_neighbor[0] - normal[0]*dot_product(normal, pft_neighbor);
   overlap_t[1] = dt*vrel_t[1] + pft_neighbor[1] - normal[1]*dot_product(normal, pft_neighbor);
   overlap_t[2] = dt*vrel_t[2] + pft_neighbor[2] - normal[2]*dot_product(normal, pft_neighbor);

    // static friction (non-Coulomb tangential force)
   amrex::Real fnmd = sqrt(dot_product(local_fn, local_fn));  // normal force mag
   local_ft[0] = -kt_des*overlap_t[0] - etat_des*vrel_t[0];
   local_ft[1] = -kt_des*overlap_t[1] - etat_des*vrel_t[1];
   local_ft[2] = -kt_des*overlap_t[2] - etat_des*vrel_t[2];
   amrex::Real ftsmd = sqrt(dot_product(local_ft, local_ft));  //f_t static mag
   amrex::Real ftdmd = mew*fnmd;                               //f_t dynamic mag

   if (ftsmd > ftdmd) {
       amrex::RealVect tangent(0.);
       tangent[0] = -local_ft[0]/ftsmd;
       tangent[1] = -local_ft[1]/ftsmd;
       tangent[2] = -local_ft[2]/ftsmd;

       //dynamic (Coulomb)  tangential friction
       local_ft[0] = -ftdmd * tangent[0];
       local_ft[1] = -ftdmd * tangent[1];
       local_ft[2] = -ftdmd * tangent[2];

       // limit tangential spring so that static friction = dynamic friction
       overlap_t[0] = (ftdmd*tangent[0] - etat_des*vrel_t[0])/kt_des;
       overlap_t[1] = (ftdmd*tangent[1] - etat_des*vrel_t[1])/kt_des;
       overlap_t[2] = (ftdmd*tangent[2] - etat_des*vrel_t[2])/kt_des;
   }

   tow_roll = RollingFriction(p1radius, p2radius, fnmd, p1omega, p2omega);

   // update tangential history contact data
   if (contact_idx == -1) {
      for (int n(0); n < max_contacts_tan_history; ++n) {
         if ((ptile_data.m_runtime_idata[cpu_id_idx + n][i] == -1) &&
             (ptile_data.m_runtime_idata[pid_idx + n][i] == -1)) {
            contact_idx = n;
            break;
         }
      }
   }
   AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (contact_idx != -1), "Particle neighbor contacts exceeded");
   ptile_data.m_runtime_idata[cpu_id_idx + contact_idx][i] = p2_cpu;
   ptile_data.m_runtime_idata[pid_idx + contact_idx][i] = p2_id;
   touch[contact_idx] = true;
   ptile_data.m_runtime_rdata[pft_neighbor_idx + contact_idx][i] = overlap_t[0];
   ptile_data.m_runtime_rdata[pft_neighbor_idx + contact_idx
      + max_contacts_tan_history][i] = overlap_t[1];
   ptile_data.m_runtime_rdata[pft_neighbor_idx + contact_idx
      + 2*max_contacts_tan_history][i] = overlap_t[2];
}

//**********************************************************
// Particle-wall collisions
//**********************************************************
template <class F1>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
particle_walls (const amrex::GpuArray<amrex::ParticleReal*,SoArealData::count>& p_realarray,
                const amrex::GpuArray<int*,SoAintData::count>& p_intarray,
                const int& i,
                const int& ls_refinement,
                const amrex::Array4<const amrex::Real>& phiarr,
                const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& plo,
                const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& dxi,
                const amrex::Real& dt,
                const amrex::RealVect& pos1,
                const MFIXSolidsParms& solids_parms,
                const int a_include_wall_conduction,
                amrex::Array4<amrex::Real const> const& a_T_eb,
                const amrex::Real& k_g,
                const amrex::Real& kn_w,
                const amrex::Real& kt_w,
                const MFIXDEM::A1D& etan_w,
                const MFIXDEM::A1D& etat_w,
                const amrex::Real& mew_w,
                amrex::RealVect& total_force,
                amrex::RealVect& total_tow_force,
                amrex::Real* cond_ptr,
                const int& istate,
                F1 RollingFriction,
                bool tan_history,
                const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                const int& max_contacts_tan_history,
                const int& pft_neighbor_idx,
                const int& cpu_id_idx,
                const int& pid_idx,
                const int& wall_pid_idx,
                bool* touch_ptr)
{
  constexpr amrex::Real small_number = 1.e-15;

  amrex::Real rp = p_realarray[SoArealData::radius][i];

  amrex::RealVect pos(pos1);

  amrex::Real ls_value = interp_level_set(pos, ls_refinement, phiarr, plo, dxi);

  amrex::Real overlap_n = rp - ls_value;

  // PFW conduction
  amrex::Real Tp1(0);
  amrex::Real Tp2(0);

  if ( a_include_wall_conduction ) {

    amrex::Real const FLPC  = solids_parms.get_flpc();
    amrex::Real const Rlens = (1.0 + FLPC)*rp;

    if (ls_value < Rlens) {

      const amrex::Real Rough = solids_parms.get_min_cond();

      Tp1 = p_realarray[SoArealData::temperature][i];
      Tp2 = Tp1;

      // Construct a point inside the wall (to machine precision)
      amrex::RealVect normal(0.);

      level_set_normal(pos, ls_refinement, normal, phiarr, plo, dxi);
      normal[0] *= -1;
      normal[1] *= -1;
      normal[2] *= -1;

      amrex::RealVect posw = normal*(ls_value + small_number) + pos1;

      amrex::Real Wx = (posw[0] - plo[0]) * dxi[0];
      amrex::Real Wy = (posw[1] - plo[1]) * dxi[1];
      amrex::Real Wz = (posw[2] - plo[2]) * dxi[2];

      int ip = static_cast<int>(amrex::Math::floor(Wx));
      int jp = static_cast<int>(amrex::Math::floor(Wy));
      int kp = static_cast<int>(amrex::Math::floor(Wz));

      // Protect against 'adiabatic walls' that have no temperature.
      if ( a_T_eb(ip,jp,kp) > std::numeric_limits<amrex::Real>::min()) {

        Tp2 = a_T_eb(ip,jp,kp);

        amrex::Real Q_dot = 2.*des_pfp_conduction( ls_value,
            rp, Rlens, Rough, k_g, Tp1, Tp2);

        amrex::HostDevice::Atomic::Add(&cond_ptr[i], Q_dot);
      }
    }
  }

  if (tan_history) {
    // Keep track of which wall contacts were touched during the update
    bool* touch = &touch_ptr[i*max_contacts_tan_history];
    for(int n(0); n < max_contacts_tan_history; ++n) {
       if (ptile_data.m_runtime_idata[pid_idx + n][i] == wall_pid_idx) {
          touch[n] = false;
       } else {
          touch[n] = true;
       }
    }
  }

  if ( ls_value < rp ) {

    // PP conduction (Tw already found from PFW conduction hit)
    if ( a_include_wall_conduction ) {
      const amrex::Real kp1 = solids_parms.calc_kp_sn(Tp1, 0);
      const amrex::Real kp2 = solids_parms.calc_kp_sn(Tp2, 0);
      amrex::Real Q_dot = des_pp_conduction(ls_value+rp,rp,rp,kp1,kp2,Tp1,Tp2);
      amrex::HostDevice::Atomic::Add(&cond_ptr[i],Q_dot);
    }

    amrex::RealVect normal(0.);
    level_set_normal(pos, ls_refinement, normal, phiarr, plo, dxi);

    normal[0] *= -1;
    normal[1] *= -1;
    normal[2] *= -1;

    amrex::RealVect v_rot(0.);
    v_rot[0] = ls_value * p_realarray[SoArealData::omegax][i];
    v_rot[1] = ls_value * p_realarray[SoArealData::omegay][i];
    v_rot[2] = ls_value * p_realarray[SoArealData::omegaz][i];

    amrex::RealVect vreltrans(0.);
    amrex::RealVect cprod(0.);

    cross_product(v_rot, normal, cprod);
    vreltrans[0] = p_realarray[SoArealData::velx][i] + cprod[0];
    vreltrans[1] = p_realarray[SoArealData::vely][i] + cprod[1];
    vreltrans[2] = p_realarray[SoArealData::velz][i] + cprod[2];

    amrex::Real vreltrans_norm = dot_product(vreltrans, normal);

    amrex::RealVect vrel_t(0.);
    vrel_t[0] = vreltrans[0] - vreltrans_norm*normal[0];
    vrel_t[1] = vreltrans[1] - vreltrans_norm*normal[1];
    vrel_t[2] = vreltrans[2] - vreltrans_norm*normal[2];

    const int phase = p_intarray[SoAintData::phase][i];
    const int phase_idx = MFIXSolidsPhase::phase_to_index(phase);

    amrex::Real kn_des_w   = kn_w;
    amrex::Real etan_des_w = etan_w(phase_idx);

    amrex::Real kt_des_w   = kt_w;
    amrex::Real etat_des_w = etat_w(phase_idx);

    amrex::RealVect local_fn(0.);
    amrex::RealVect local_ft(0.);
    amrex::RealVect overlap_t(0.);

    // calculate the normal contact force
    local_fn[0] = -(kn_des_w*overlap_n*normal[0]
                  + etan_des_w*vreltrans_norm*normal[0]);
    local_fn[1] = -(kn_des_w*overlap_n*normal[1]
                  + etan_des_w*vreltrans_norm*normal[1]);
    local_fn[2] = -(kn_des_w*overlap_n*normal[2]
                  + etan_des_w*vreltrans_norm*normal[2]);

    // torque due to rolling friction
    amrex::RealVect tow_roll(0.);

    if (tan_history) {

      bool* touch = &touch_ptr[i*max_contacts_tan_history];
      amrex::RealVect p1omega(p_realarray[SoArealData::omegax][i],
                              p_realarray[SoArealData::omegay][i],
                              p_realarray[SoArealData::omegaz][i]);

      compute_pw_collisions(ptile_data, local_ft, overlap_t, tow_roll, touch, i,
            rp, p1omega, local_fn, vrel_t, normal, dt, mew_w, kt_des_w, etat_des_w,
            RollingFriction, max_contacts_tan_history, pft_neighbor_idx, cpu_id_idx, pid_idx, wall_pid_idx);

    } else {

      // calculate the tangential displacement
      amrex::Real mag_overlap_t(0.);
      overlap_t[0] = dt*vrel_t[0];
      overlap_t[1] = dt*vrel_t[1];
      overlap_t[2] = dt*vrel_t[2];

      mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

      if (mag_overlap_t > 0.0) {
          amrex::Real fnmd = mew_w * sqrt(dot_product(local_fn, local_fn));
          amrex::RealVect tangent(0.);
          tangent[0] = overlap_t[0]/mag_overlap_t;
          tangent[1] = overlap_t[1]/mag_overlap_t;
          tangent[2] = overlap_t[2]/mag_overlap_t;
          local_ft[0] = -fnmd * tangent[0];
          local_ft[1] = -fnmd * tangent[1];
          local_ft[2] = -fnmd * tangent[2];
      } else {
          local_ft[0] = 0.0;
          local_ft[1] = 0.0;
          local_ft[2] = 0.0;
      }
    }

    if ( istate > 0 ) { // normal particles

      total_force[0] += local_fn[0] + local_ft[0];
      total_force[1] += local_fn[1] + local_ft[1];
      total_force[2] += local_fn[2] + local_ft[2];

      amrex::RealVect tow_force(0.);

      cross_product(normal, local_ft, tow_force);

      total_tow_force[0] += ls_value*tow_force[0];
      total_tow_force[1] += ls_value*tow_force[1];
      total_tow_force[2] += ls_value*tow_force[2];

      total_tow_force += tow_roll;

    } else { // entering particles

      amrex::Real velx = p_realarray[SoArealData::velx][i];
      amrex::Real vely = p_realarray[SoArealData::vely][i];
      amrex::Real velz = p_realarray[SoArealData::velz][i];

      amrex::Real velmag = std::sqrt(velx*velx + vely*vely + velz*velz);

      amrex::Real dotprod = (normal[0] * velx +
                             normal[1] * vely +
                             normal[2] * velz)/velmag;

      // This is to catch particles that are not moving normal to
      // the levelset so that we can adjust their velocity and make sure
      // they fully enter the domain.
      if(amrex::Math::abs(1.0 + dotprod) > std::numeric_limits<amrex::Real>::epsilon()) {

        p_realarray[SoArealData::velx][i] = -velmag*normal[0];
        p_realarray[SoArealData::vely][i] = -velmag*normal[1];
        p_realarray[SoArealData::velz][i] = -velmag*normal[2];

      }

    }

  // An entering particle is no longer overlapping the wall.
  } else if(istate == 0) {
    //amrex::AllPrint() << "setting particle to normal\n";

    // Set the state to normal so it no longer ignores forces.
    p_intarray[SoAintData::state][i] = 1;
  }

  if (tan_history) {
    // Clean up tangential history contacts
    bool* touch = &touch_ptr[i*max_contacts_tan_history];
    for(int n(0); n < max_contacts_tan_history; ++n) {
       if ((touch[n] == false) &&
           (ptile_data.m_runtime_idata[pid_idx + n][i] == wall_pid_idx)) {
          ptile_data.m_runtime_idata[cpu_id_idx + n][i] = -1;
          ptile_data.m_runtime_idata[pid_idx + n][i] = -1;
       }
    }
  }

}


//**********************************************************
// Particle-particle collisions
//**********************************************************
template <class F1>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
particle_particles (const MFIXParticleContainer::ParticleType& particle,
                    const amrex::Neighbors<MFIXParticleContainer::ParticleType>& neighbs,
                    const amrex::GpuArray<amrex::ParticleReal*, SoArealData::count>& p_realarray,
                    const amrex::GpuArray<int*, SoAintData::count>& p_intarray,
                    const int& i,
                    const int& solve_enthalpy,
                    const amrex::Real& dt,
                    const int& nrp,
                    const MFIXSolidsParms& solids_parms,
                    const amrex::Real& k_g,
                    const amrex::Real& kn,
                    const amrex::Real& kt,
                    const MFIXDEM::A2D& etan,
                    const MFIXDEM::A2D& etat,
                    const amrex::Real& mew,
                    amrex::RealVect& total_force,
                    amrex::RealVect& total_tow_force,
                    amrex::Real* fc_ptr,
                    amrex::Real* tow_ptr,
                    amrex::Real* cond_ptr,
                    const int& ntot,
                    const int& istate,
                    F1 RollingFriction,
                    bool tan_history,
                    const amrex::ParticleTileData<amrex::Particle<0,0>,SoArealData::count,SoAintData::count>& ptile_data,
                    const int& max_contacts_tan_history,
                    const int& pft_neighbor_idx,
                    const int& cpu_id_idx,
                    const int& pid_idx,
                    const int& wall_pid_idx,
                    bool* touch_ptr)
{
  constexpr amrex::Real small_number = 1.e-15;

  const amrex::RealVect& pos1(particle.pos());

  if (tan_history) {
    // Keep track of which contacts were touched during the update
    bool* touch = &touch_ptr[i*max_contacts_tan_history];
    for(int n(0); n < max_contacts_tan_history; ++n) {
       if (ptile_data.m_runtime_idata[pid_idx + n][i] != wall_pid_idx) {
          touch[n] = false;
       } else {
          touch[n] = true;
       }
    }
  }

  for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
  {
      const auto p2 = *mit;
      const int j = mit.index();

      amrex::Real dist_x = p2.pos(0) - pos1[0];
      amrex::Real dist_y = p2.pos(1) - pos1[1];
      amrex::Real dist_z = p2.pos(2) - pos1[2];

      amrex::Real r2 = dist_x*dist_x +
                       dist_y*dist_y +
                       dist_z*dist_z;

      const amrex::Real p1radius = p_realarray[SoArealData::radius][i];
      const amrex::Real p2radius = p_realarray[SoArealData::radius][j];

      amrex::Real r_lm = p1radius + p2radius;

      AMREX_ASSERT_WITH_MESSAGE(
          !(particle.id() == p2.id() && particle.cpu() == p2.cpu()),
          "A particle should not be its own neighbor!");

      // PFP conduction
      if(solve_enthalpy && solids_parms.get_do_pfp_cond()) {
        const amrex::Real FLPC = solids_parms.get_flpc();
        amrex::Real Rp_eff     = 2.0*(p1radius*p2radius)/(p1radius + p2radius);
        amrex::Real Rlens_eff  = (1.0 + FLPC)*Rp_eff;
        amrex::Real lens_lm    = 2.0*Rlens_eff;
        if ( r2 <= (lens_lm - small_number)*(lens_lm - small_number) ) {
          const amrex::Real Rough  = solids_parms.get_min_cond();
          const amrex::Real Tp1    = p_realarray[SoArealData::temperature][i];
          const amrex::Real Tp2    = p_realarray[SoArealData::temperature][j];
          amrex::Real dist_mag_eff = sqrt(r2)/2.0; // Two particles with a midpoint wall
          amrex::Real Q_dot = des_pfp_conduction(dist_mag_eff,Rp_eff,Rlens_eff,
                                          Rough,k_g,Tp1,Tp2);
          if(i < nrp) amrex::HostDevice::Atomic::Add(&cond_ptr[i], Q_dot);
          if(j < nrp) amrex::HostDevice::Atomic::Add(&cond_ptr[j],-Q_dot);
        }
      }

      if ( r2 <= (r_lm - small_number)*(r_lm - small_number) ) {

          const int jstate = p_intarray[SoAintData::state][j];

          amrex::Real dist_mag = sqrt(r2);

          AMREX_ASSERT(dist_mag >= std::numeric_limits<amrex::Real>::epsilon());

          // PP conduction
          if(solve_enthalpy && solids_parms.get_do_pfp_cond()) {
            const amrex::Real Tp1 = p_realarray[SoArealData::temperature][i];
            const amrex::Real Tp2 = p_realarray[SoArealData::temperature][j];
            const amrex::Real kp1 = solids_parms.calc_kp_sn(Tp1,0);
            const amrex::Real kp2 = solids_parms.calc_kp_sn(Tp2,0);
            amrex::Real Q_dot = des_pp_conduction(dist_mag,p1radius,p2radius,kp1,kp2,Tp1,Tp2);
            if(i < nrp) amrex::HostDevice::Atomic::Add(&cond_ptr[i], Q_dot);
            if(j < nrp) amrex::HostDevice::Atomic::Add(&cond_ptr[j],-Q_dot);
          }

          amrex::Real dist_mag_inv = 1.e0/dist_mag;

          amrex::RealVect normal(0.);
          normal[0] = dist_x * dist_mag_inv;
          normal[1] = dist_y * dist_mag_inv;
          normal[2] = dist_z * dist_mag_inv;

          amrex::Real overlap_n(0.);

          overlap_n = r_lm - dist_mag;

          amrex::Real vrel_trans_norm;
          amrex::RealVect vrel_t(0.);

          amrex::RealVect p1vel(p_realarray[SoArealData::velx][i],
                                p_realarray[SoArealData::vely][i],
                                p_realarray[SoArealData::velz][i]);

          amrex::RealVect p2vel(p_realarray[SoArealData::velx][j],
                                p_realarray[SoArealData::vely][j],
                                p_realarray[SoArealData::velz][j]);

          amrex::RealVect p1omega(p_realarray[SoArealData::omegax][i],
                                  p_realarray[SoArealData::omegay][i],
                                  p_realarray[SoArealData::omegaz][i]);

          amrex::RealVect p2omega(p_realarray[SoArealData::omegax][j],
                                  p_realarray[SoArealData::omegay][j],
                                  p_realarray[SoArealData::omegaz][j]);

          cfrelvel(p1vel, p2vel, p1radius, p2radius, p1omega,
                   p2omega, vrel_trans_norm, vrel_t, normal, dist_mag);

          const int phase1 = p_intarray[SoAintData::phase][i];
          const int phase2 = p_intarray[SoAintData::phase][j];

          const int phase1_idx = MFIXSolidsPhase::phase_to_index(phase1);
          const int phase2_idx = MFIXSolidsPhase::phase_to_index(phase2);

          amrex::Real kn_des   = kn;
          amrex::Real kt_des   = kt;
          amrex::Real etan_des = etan(phase1_idx, phase2_idx);
          amrex::Real etat_des = etat(phase1_idx, phase2_idx);

          amrex::RealVect local_fn(0.);
          amrex::RealVect local_ft(0.);
          amrex::RealVect overlap_t(0.);

          // calculate the normal contact force
          local_fn[0] = -(kn_des*overlap_n*normal[0]
                        + etan_des*vrel_trans_norm*normal[0]);
          local_fn[1] = -(kn_des*overlap_n*normal[1]
                        + etan_des*vrel_trans_norm*normal[1]);
          local_fn[2] = -(kn_des*overlap_n*normal[2]
                        + etan_des*vrel_trans_norm*normal[2]);

          // torque due to rolling friction
          amrex::RealVect tow_roll(0.);

          if (tan_history) {

            bool* touch = &touch_ptr[i*max_contacts_tan_history];
            compute_pp_collisions(ptile_data, local_ft, overlap_t, tow_roll, touch, i, p2.cpu(), p2.id(),
                  p1radius, p2radius, p1omega, p2omega,
                  local_fn, vrel_t, normal, dt, mew, kt_des, etat_des,
                  RollingFriction, max_contacts_tan_history, pft_neighbor_idx, cpu_id_idx, pid_idx);

          } else {

            // calculate the tangential overlap
            amrex::Real mag_overlap_t(0.);
            overlap_t[0] = dt*vrel_t[0];
            overlap_t[1] = dt*vrel_t[1];
            overlap_t[2] = dt*vrel_t[2];
            mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

            if (mag_overlap_t > 0.0) {
                amrex::Real fnmd = mew * sqrt(dot_product(local_fn, local_fn));
                amrex::RealVect tangent(0.);
                tangent[0] = overlap_t[0]/mag_overlap_t;
                tangent[1] = overlap_t[1]/mag_overlap_t;
                tangent[2] = overlap_t[2]/mag_overlap_t;
                local_ft[0] = -fnmd * tangent[0];
                local_ft[1] = -fnmd * tangent[1];
                local_ft[2] = -fnmd * tangent[2];

            } else {
                local_ft[0] = 0.0;
                local_ft[1] = 0.0;
                local_ft[2] = 0.0;
            }
          }

          amrex::Real dist_cl1 = 0.5 * (dist_mag + (p1radius*p1radius - p2radius*p2radius) * dist_mag_inv);
          dist_cl1 = dist_mag - dist_cl1;

          amrex::Real dist_cl2 = 0.5 * (dist_mag + (p2radius*p2radius - p1radius*p1radius) * dist_mag_inv);
          dist_cl2 = dist_mag - dist_cl2;

          amrex::RealVect local_tow_force(0.);
          cross_product(normal, local_ft, local_tow_force);

          if ( istate == state::normal ) {
            total_force[0] += local_fn[0] + local_ft[0];
            total_force[1] += local_fn[1] + local_ft[1];
            total_force[2] += local_fn[2] + local_ft[2];

            total_tow_force[0] += dist_cl1*local_tow_force[0] + tow_roll[0];
            total_tow_force[1] += dist_cl1*local_tow_force[1] + tow_roll[1];
            total_tow_force[2] += dist_cl1*local_tow_force[2] + tow_roll[2];
          }

          if (j < nrp && jstate != state::entering ) {
            amrex::HostDevice::Atomic::Add(&fc_ptr[j         ], -(local_fn[0] + local_ft[0]));
            amrex::HostDevice::Atomic::Add(&fc_ptr[j + ntot  ], -(local_fn[1] + local_ft[1]));
            amrex::HostDevice::Atomic::Add(&fc_ptr[j + 2*ntot], -(local_fn[2] + local_ft[2]));

            amrex::HostDevice::Atomic::Add(&tow_ptr[j         ], dist_cl2*local_tow_force[0] - tow_roll[0]);
            amrex::HostDevice::Atomic::Add(&tow_ptr[j + ntot  ], dist_cl2*local_tow_force[1] - tow_roll[1]);
            amrex::HostDevice::Atomic::Add(&tow_ptr[j + 2*ntot], dist_cl2*local_tow_force[2] - tow_roll[2]);
          }
          // Special case of two entering particles having an overlap
          if (istate == state::entering && jstate == state::entering) {

            const amrex::Real shift = 1.0001*overlap_n;
            const amrex::RealVect sumvel(p1vel + p2vel);
            const int imove = (( sumvel[0]*normal[0]
                               + sumvel[1]*normal[1]
                               + sumvel[2]*normal[2]) > 0.) ? 1 : 0;

            if (imove) {
              total_force[0] -= shift * normal[0];
              total_force[1] -= shift * normal[1];
              total_force[2] -= shift * normal[2];

            } else if (j < nrp) {
              {
                amrex::HostDevice::Atomic::Add(&fc_ptr[j         ], shift * normal[0]);
                amrex::HostDevice::Atomic::Add(&fc_ptr[j +   ntot], shift * normal[1]);
                amrex::HostDevice::Atomic::Add(&fc_ptr[j + 2*ntot], shift * normal[2]);
              }
            }
          } // end overlap between entering particles

      } // end overlap
  } // end neighbor loop

  if(i < nrp) {
      amrex::HostDevice::Atomic::Add(&fc_ptr[i         ], total_force[0]);
      amrex::HostDevice::Atomic::Add(&fc_ptr[i + ntot  ], total_force[1]);
      amrex::HostDevice::Atomic::Add(&fc_ptr[i + 2*ntot], total_force[2]);

      amrex::HostDevice::Atomic::Add(&tow_ptr[i         ], total_tow_force[0]);
      amrex::HostDevice::Atomic::Add(&tow_ptr[i + ntot  ], total_tow_force[1]);
      amrex::HostDevice::Atomic::Add(&tow_ptr[i + 2*ntot], total_tow_force[2]);

  }

  if (tan_history) {
    // Clean up tangential history contacts
    bool* touch = &touch_ptr[i*max_contacts_tan_history];
    for(int n(0); n < max_contacts_tan_history; ++n) {
       if ((touch[n] == false) &&
           (ptile_data.m_runtime_idata[pid_idx + n][i] != wall_pid_idx)) {
          ptile_data.m_runtime_idata[cpu_id_idx + n][i] = -1;
          ptile_data.m_runtime_idata[pid_idx + n][i] = -1;
       }
    }
  }
}

#endif
