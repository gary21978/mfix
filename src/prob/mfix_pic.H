#ifndef MFIX_PIC_PARMS_H_
#define MFIX_PIC_PARMS_H_

#include <AMReX_REAL.H>

#define MFP_DISABLED 0
#define SCALE_MFP_DISABLED 0
#define SCALE_MFP_CHAPMAN  0
#define SCALE_MFP_CUTCHIS  0
#define SCALE_MFP_MUSSER   1


class MFIXPIC {

  public:

    struct InitialStepType { enum {Undefined=-1, zero_eps, nth_eps, taylor_approx}; };

    MFIXPIC () = default;
    ~MFIXPIC () = default;

    void Initialize ();

    int solve () const { return m_solve; }
    void set_solve (const int val) { m_solve = val; }

    int NPHASE () const { return m_NPHASE; }

    int verbose () const { return m_verbose; }

    // Coefficients for computing the solid stress
    amrex::Real Ps () const { return m_Ps; }
    amrex::Real beta () const { return m_beta; }
    amrex::Real ep_cp () const { return m_ep_cp; }

    amrex::Real parcels_per_cell_at_pack () const { return m_parcels_per_cell_at_pack; }

    amrex::Real vel_ref_frame () const { return m_vel_ref_frame; }

    amrex::Real small_number () const { return m_small_number; }
    amrex::Real damping_factor () const { return m_damping_factor; }

    amrex::Real damping_factor_wall_normal () const { return m_damping_factor_wall_normal; }
    amrex::Real damping_factor_wall_tangent () const { return m_damping_factor_wall_tangent; }

    amrex::Real advance_vel_p () const { return m_advance_vel_p; }

    int max_iter () const { return m_max_iter; }

    int initial_step () const { return m_initial_step; }

  private:

    int m_NPHASE;

    int m_solve = 0;

    int m_verbose = 0;

    // Coefficients for computing the solid stress
    amrex::Real m_Ps = 100.0;
    amrex::Real m_beta = 3.0;
    amrex::Real m_ep_cp = 0.55;

    amrex::Real m_parcels_per_cell_at_pack;

    amrex::Real m_vel_ref_frame = 0.5;

    amrex::Real m_small_number = 1.e-7;
    amrex::Real m_damping_factor = 0.4;

    amrex::Real m_damping_factor_wall_normal = 0.3;
    amrex::Real m_damping_factor_wall_tangent = 0.99;

    amrex::Real m_advance_vel_p = 0.5;

    int m_max_iter = 3;

    int m_initial_step = InitialStepType::nth_eps;
};

#endif
