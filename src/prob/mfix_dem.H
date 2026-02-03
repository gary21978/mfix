#ifndef MFIX_DEM_PARMS_H_
#define MFIX_DEM_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <mfix_solids.H>

class MFIXDEM
{
  public:
    static constexpr int NMAX = 10;

    // Constructor
    MFIXDEM ();

    // Destructor
    ~MFIXDEM ();

    struct CollisionModel { enum {Undefined=-1, LSD}; };

    template <typename ARRAY_TYPE>
    struct DeviceArray
    {
        typedef ARRAY_TYPE array_type;

        template <typename... Is>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::Real& operator() (Is... index) const noexcept {
            return (*a)(index...);
        }

        void alloc () {
#ifdef AMREX_USE_GPU
            auto p = (ARRAY_TYPE*)amrex::The_Device_Arena()->alloc(sizeof(ARRAY_TYPE));
            amrex::single_task(amrex::Gpu::gpuStream(), [p] AMREX_GPU_DEVICE
            {
                new (p) ARRAY_TYPE;
            });
            amrex::Gpu::synchronize();
            a = p;
#else
            a = new ARRAY_TYPE;
#endif
        }

        void free () {
#ifdef AMREX_USE_GPU
            amrex::single_task(amrex::Gpu::gpuStream(), [p=a] AMREX_GPU_DEVICE
            {
                p->~ARRAY_TYPE();
            });
            amrex::Gpu::synchronize();
            amrex::The_Device_Arena()->free(a);
#else
            delete a;
#endif
        }

        ARRAY_TYPE* arrayPtr () { return a; }

    private:
        ARRAY_TYPE * a = nullptr;
    };

    using A2D = DeviceArray<amrex::Array2D<amrex::Real,0,NMAX-1,0,NMAX-1,amrex::Order::C> >;
    using A1D = DeviceArray<amrex::Array1D<amrex::Real,0,NMAX-1> >;

    void Initialize ();

    int NPHASE () const { return m_NPHASE; }

    int collision_model () const { return m_collision_model; }

    int solve () const { return m_solve; }
    void set_solve (const int val) { m_solve = val; }

    // implicit drag?
    int implicit_drag () const { return m_implicit_drag; }

    amrex::Real dtsolid () const { return m_dtsolid; }
    void set_dtsolid (const amrex::Real val) { m_dtsolid = val; }

    // tangential component of spring constant
    amrex::Real kt () const { return m_kt; }
    amrex::Real kt_w () const { return m_kt_w; }

    // normal component of spring constant
    amrex::Real kn () const { return m_kn; }
    amrex::Real kn_w () const { return m_kn_w; }

    // friction coefficients for particle-particle and particle-wall collisions
    amrex::Real mew () const { return m_mew; }
    amrex::Real mew_w () const { return m_mew_w; }

    // assumed to be a constant factor relating these for linear spring-dashpot
    amrex::Real kt_fac () const { return m_kt_fac; }
    amrex::Real kt_w_fac () const { return m_kt_w_fac; }

    // fluid thermal conductivity (see issues/288)
    amrex::Real k_g_dem () const { return m_k_g_dem; }

    // normal and tangential components of the damping coefficients
    A2D etan () const { return m_etan; }
    A1D etan_w () const { return m_etan_w; }

    A2D etat () const { return m_etat; }
    A1D etat_w () const { return m_etat_w; }

    // coefficients of restitution, normal and tangential
    A2D en () const { return m_en; }
    A1D en_w () const { return m_en_w; }

    amrex::Real eta_fac () const { return m_eta_fac; }
    amrex::Real eta_w_fac () const { return m_eta_w_fac; }

    // Tangential damping factor ratio
    amrex::Real eta_fac_pp () const { return m_eta_fac_pp; }
    amrex::Real eta_fac_pw () const { return m_eta_fac_pw; }

    amrex::Real neighborhood () const { return m_neighborhood; }
    void set_neighborhood (const amrex::Real val) { m_neighborhood = val; }

    // Polydisperse neighbor search
    amrex::Gpu::DeviceVector<amrex::Real>& get_pneighborhood () { return m_pneighborhood; }
    amrex::Real* pneighdata () { return m_pneighborhood.data(); }

    amrex::Vector<amrex::Real>& get_pbin () { return m_pbin; }
    amrex::Real* pbindata  () { return m_pbin.data();  }

    int* prefratdata () { return m_prefrats.data(); }

    bool pneig_flag () { return m_pneigh_flag; }
    int nptypes ()     { return m_nptypes; }
    int nppairs ()     { return ( (m_nptypes*(m_nptypes+1)) / 2 ); }

    // Coarse-grain DEM
    int cg_dem () const { return m_cg_dem; }

    amrex::Real cg_particles_per_cell_at_pack () const
    { return m_cg_particles_per_cell_at_pack; }

    int rolling_friction () const { return m_rolling_friction; }

    amrex::Real rolling_friction_coeff () const { return m_rolling_friction_coefficient; }

    bool tan_history () const { return m_tan_history; }
    int tan_history_max_contacts () const { return m_tan_history_max_contacts; }

  private:
    int m_NPHASE;

    int m_collision_model;

    int m_solve;

    int m_implicit_drag;

    amrex::Real m_dtsolid;

    // tangential component of spring constant
    amrex::Real m_kt;
    amrex::Real m_kt_w;

    // normal component of spring constant
    amrex::Real m_kn;
    amrex::Real m_kn_w;

    // friction coefficients for particle-particle and particle-wall collisions
    amrex::Real m_mew;
    amrex::Real m_mew_w;

    // assumed to be a constant factor relating these for linear spring-dashpot
    amrex::Real m_kt_fac;
    amrex::Real m_kt_w_fac;

    // fluid thermal conductivity (see issues/288)
    amrex::Real m_k_g_dem;

    // normal and tangential components of the damping coefficients
    A2D m_etan;
    A1D m_etan_w;

    A2D m_etat;
    A1D m_etat_w;

    // coefficients of restitution, normal and tangential
    A2D m_en;
    A1D m_en_w;

    amrex::Real m_eta_fac;
    amrex::Real m_eta_w_fac;

    // Tangential damping factor ratio
    amrex::Real m_eta_fac_pp;
    amrex::Real m_eta_fac_pw;

    amrex::Real m_neighborhood;

    // Polydisperse neighbor search
    amrex::Gpu::DeviceVector<amrex::Real> m_pneighborhood;  // Neighborhood criteria for i-j type pair
    amrex::Vector<amrex::Real> m_pbin;                      // Bin size for particle `types'
    amrex::Vector<int> m_prefrats;                          // Refinement ratios that control bin sizes
    bool m_pneigh_flag{false};                              // Flag to use Poly-search algorithm
    int  m_nptypes{1};                                      // Number of types

    // Coarse-grain DEM
    int m_cg_dem;
    amrex::Real m_cg_particles_per_cell_at_pack;

    int m_rolling_friction;
    amrex::Real m_rolling_friction_coefficient;

    bool m_tan_history = false;
    int m_tan_history_max_contacts = 0;
};

struct RollingFrictionModel { enum {Unspecified=-1, None, ModelA, ModelB}; };

#endif
