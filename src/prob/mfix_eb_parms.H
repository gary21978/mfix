#ifndef MFIX_EB_PARMS_H_
#define MFIX_EB_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>

#include <AMReX_Vector.H>
#include <AMReX_LO_BCTYPES.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>

#include <mfix_fluid_var.H>

class MFIXEmbeddedBoundaries
{
  public:
    // Constructor
    MFIXEmbeddedBoundaries () = default;
    ~MFIXEmbeddedBoundaries () = default;

    void set_has_flow () { m_has_flow = true; }

    void set_has_temperature () { m_has_temperature = true; }

    void set_needs_avg_particle_temperature ()
    { m_needs_avg_particle_temperature = true; }

    bool has_temperature () const { return m_has_temperature; }

    bool needs_avg_particle_temperature () const
    { return m_needs_avg_particle_temperature; }

    bool has_flow () const { return m_has_flow; }

  private:

    // Flag to specify if we have Dirichlet BCs on the EB
    bool m_has_temperature = false;

    // Flag stating we need to average solids temperature
    bool m_needs_avg_particle_temperature = false;

    // Flag to specify that we have Dirichlet flow though EB face
    bool m_has_flow = false;
};


struct EB_t {

  EB_t ()
    : m_temperature({1})
    , has_normal(0)
    , normal(0.)
    , normal_tol(0.)
  {}

  bool constant_temperature ( ) const noexcept
  { return temperature_bc_type == temperature_bc::Constant; }

  bool conjugate_heat_transfer_1D ( ) const noexcept
  { return temperature_bc_type == temperature_bc::CHT_1D; }

  void set_constant_temperature ( ) noexcept {
    temperature_bc_type = temperature_bc::Constant;
  }

  void set_conjugate_heat_transfer_1D ( ) noexcept {
    temperature_bc_type = temperature_bc::CHT_1D;
  }

  amrex::Real temperature ( const amrex::Real a_time = 0.) noexcept {
    return (m_temperature(a_time).at(0));
  }

  // Inhomogeneous Dirichlet BCs for fluid temperature solver
  // Adiabatic : No heat flux through walls
  // Constant  : Walls have a user-defined constant temperature
  // CHT_1D    : 1D conjugate heat transfer defines wall temperature
  enum struct temperature_bc { Adiabatic, Constant, CHT_1D };

  temperature_bc temperature_bc_type = temperature_bc::Adiabatic;

  // Temperature to set on EB walls
  FVAR_ m_temperature;

  // 1D Conjugate heat transfer: Tw = T - q/h_int where
  //
  // q = (T - T_env) / (1/h_int + L_wall/k_wall + 1/h_ext)
  //
  // T : temperature (K) inside boundary cell. If the flag 'phase_averaged_temperature'
  //     is true, then T = ε*Tf + (1-ε)*Tp where Tp is a volume averaged particle
  //     temperature, otherwise T = Tf.
  //
  // T_env : environment temperature (K)
  amrex::Real T_env = 298.15;
  // h_int : convective heat transfer coefficient, fluid-wall (W/m2.K)
  amrex::Real h_int = std::numeric_limits<amrex::Real>::max();
  // h_ext : convective heat transfer coefficient, wall-environment (W/m2.K)
  amrex::Real h_ext = std::numeric_limits<amrex::Real>::max();
  // k_wall : wall thermal conductivity (W/m.K)
  amrex::Real k_wall = std::numeric_limits<amrex::Real>::max();
  // L_wall : wall thickness (m)
  amrex::Real L_wall = 0.;
  // flag to use phase averaged temperature or fluid temperature
  bool phase_averaged_temperature = false;

  // Normal to exclude EB faces. By default all faces within a specified
  // region will
  int has_normal;
  amrex::Vector<amrex::Real> normal;
  amrex::Real normal_tol;
};

#endif
