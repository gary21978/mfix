#include <mfix_init_fluid.H>

#include <mfix_calc_cell.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_ic.H>


using namespace amrex;


namespace init_fluid_aux {

// Forward declarations
void set_ic_vel (const Box& sbx, const Box& domain, const Real dx,
                 const Real dy, const Real dz, const GpuArray<Real,3>& plo,
                 MFIXInitialConditions& initial_conditions,
                 Array4<Real      > const& a_vel);

void set_ic_energy (const Box& sbx, const Box& domain, const Real dx,
                  const Real dy, const Real dz, const GpuArray<Real,3>& plo,
                  MFIXInitialConditions& initial_conditions,
                  Array4<EBCellFlag const> const& a_flags,
                  Array4<Real      > const& a_T,
                  Array4<Real      > const& a_h,
                  Array4<Real const> const& a_X,
                  MFIXFluidPhase& fluid);

void set_ic_species_g (const Box& sbx, const Box& domain, const Real dx,
                       const Real dy, const Real dz, const GpuArray<Real,3>& plo,
                       MFIXInitialConditions& initial_conditions,
                       Array4<Real      > const& a_X);

void set_ic_rho (const Box& sbx, const Box& domain, const Real dx,
                 const Real dy, const Real dz, const GpuArray<Real,3>& plo,
                  MFIXInitialConditions& initial_conditions, Array4<Real> const& a_rho);


} // end namespace init_fluid_aux


using namespace init_fluid_aux;


//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid                                              !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

void set_ic_fluid ( int const a_lev,
                    const Box& sbx,
                    const Box& bx,
                    const Box& domain,
                    const MFIter& mfi,
                    MFIXLevelData& a_ld,
                    const Real dx,
                    const Real dy,
                    const Real dz,
                    const GpuArray<Real, 3>& plo,
                    bool test_tracer_conservation,
                    MFIXInitialConditions& initial_conditions,
                    MFIXFluidPhase& fluid)
{
  // Set user specified initial conditions (IC)

  const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(a_ld.rho(a_lev)->Factory());

  EBCellFlagFab const& flagfab = factory.getMultiEBCellFlagFab()[mfi];
  Array4<EBCellFlag const> const& flags = flagfab.const_array();

  // **************************************************************************
  // Set initial fluid velocity
  // **************************************************************************
  set_ic_vel(sbx, domain, dx, dy, dz, plo, initial_conditions,
      a_ld.vel(a_lev, mfi));

  // **************************************************************************
  // Set initial fluid tracer
  // **************************************************************************
  if (fluid.solve_tracer()) {

    if (test_tracer_conservation) {

      init_periodic_tracer(bx, domain, a_ld.vel(a_lev,mfi), a_ld.tracer(a_lev, mfi), dx, dy, dz);

    } else {

      const Real trac_0 = fluid.tracer();

      Array4<Real> const& trac = a_ld.tracer(a_lev,mfi);

      int const ncomp( fluid.ntracer() );

      ParallelFor(sbx, ncomp, [trac,trac_0]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { trac(i,j,k,n) = trac_0; });
    }
  }

  // ************************************************************************
  // Set initial fluid density
  // ************************************************************************
  set_ic_rho(sbx, domain, dx, dy, dz, plo, initial_conditions,
      a_ld.rho(a_lev,mfi));

  // **************************************************************************
  // Set initial fluid species mass fractions
  // **************************************************************************
  if (fluid.solve_species()) {
    set_ic_species_g(sbx, domain, dx, dy, dz, plo, initial_conditions,
        a_ld.X(a_lev,mfi));
  }

  // **************************************************************************
  // Set initial fluid temperature
  // **************************************************************************
  if (fluid.solve_enthalpy() || fluid.constraint.isIdealGas()) {

    Array4<Real      > const& T = a_ld.T(a_lev,mfi);
    Array4<Real      > const& h = a_ld.h(a_lev,mfi);
    Array4<Real const> const& X = a_ld.X_const(a_lev,mfi);

    set_ic_energy( sbx, domain, dx, dy, dz, plo, initial_conditions,
        flags, T, h, X, fluid);

    //FArrayBox& eps_fab = (*m_eps[a_lev])[pti];

    if (!fluid.solve_enthalpy()) {

      Array4<Real      > const& T_old = a_ld.T_old(a_lev,mfi);

      ParallelFor(mfi.tilebox(), [T, T_old]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { T_old(i,j,k) = T(i,j,k); });
    }
  }
}


void init_periodic_tracer (const Box& bx,
                           const Box& domain,
                           Array4<Real> const& a_vel,
                           Array4<Real> const& a_tracer,
                           const Real dx,
                           const Real dy,
                           const Real dz)
{
    const Real twopi(2. * M_PI);

    int dir(0);
    const Real A(1.0);

    Real L(0); // domain size
    Real C(0); // sin coefficient

    dir += (domain.bigEnd(1) == domain.bigEnd(2)) ? 1 : 0;
    dir += (domain.bigEnd(0) == domain.bigEnd(2)) ? 2 : 0;
    dir += (domain.bigEnd(0) == domain.bigEnd(1)) ? 4 : 0;

    switch (dir)
      {
      case 1:  // x-direction

        L = Real(domain.bigEnd(0)+1) * dx;
        C = twopi / L;
        amrex::ParallelFor(bx,[A,C,dx,dy,dz,a_tracer,a_vel]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            Real x = (Real(i) + .5) * dx - .00037;
            Real y = (Real(j) + .5) * dy - .00073;
            Real z = (Real(k) + .5) * dz - .00123;

            a_tracer(i,j,k) = A*( std::sin(C*(y+z) - 0.00042) + 1.0) * exp(x);

            a_vel(i,j,k,1) += 0.1*( std::sin(C*(x+z) - 0.00042) + 1.0) * exp(y);
            a_vel(i,j,k,2) += 0.1*( std::sin(C*(x+y) - 0.00042) + 1.0) * exp(z);
        });
        break;

    case 2: // y-direction

        L = Real(domain.bigEnd(1)+1) * dy;
        C = twopi / L;
        amrex::ParallelFor(bx,[A,C,dx,dy,dz,a_tracer,a_vel]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            Real y = (Real(j) + .5) * dy - .00037;

            Real x = (Real(i) + .5) * dx - .00073;
            Real z = (Real(k) + .5) * dz - .00123;

            a_tracer(i,j,k) = A*( std::sin(C*(x+z) - 0.00042) + 1.0) * exp(y);

            a_vel(i,j,k,0) += 0.1*( std::sin(C*(y+z) - 0.00042) + 1.0) * exp(x);
            a_vel(i,j,k,2) += 0.1*( std::sin(C*(y+x) - 0.00042) + 1.0) * exp(z);

        });
        break;

    case 4: // z-direction

        L = Real(domain.bigEnd(2)+1) * dz;
        C = twopi / L;
        amrex::ParallelFor(bx,[A,C,dx,dy,dz,a_tracer,a_vel]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (Real(k) + .5) * dz - .00037;

            Real x = (Real(i) + .5) * dx - .00073;
            Real y = (Real(j) + .5) * dy - .00123;

            a_tracer(i,j,k) = A*( std::sin(C*(x+y) - 0.00042) + 1.0) * exp(z);

            a_vel(i,j,k,0) += 0.1*( std::sin(C*(z+y) - 0.00042) + 1.0) * exp(x);
            a_vel(i,j,k,1) += 0.1*( std::sin(C*(z+x) - 0.00042) + 1.0) * exp(y);
        });
        break;

    default:
        amrex::Abort("Error: wrong direction number");
        break;
    }
}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: init_fluid_parameters                                   !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void init_fluid_parameters ( int const a_lev,
                             const Box& bx,
                             const MFIter& mfi,
                             MFIXLevelData& a_ld,
                             MFIXFluidPhase& fluid)
{

  const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(a_ld.rho(a_lev)->Factory());

  EBCellFlagFab const& flagfab = factory.getMultiEBCellFlagFab()[mfi];
  Array4<EBCellFlag const> const& flags = flagfab.const_array();


  if ( fluid.solve_enthalpy() ) {

    Array4<Real const> const& X = a_ld.X_const(a_lev,mfi);
    Array4<Real const> const& T = a_ld.T_const(a_lev, mfi);

    Array4<Real      > const& h  = a_ld.h(a_lev,mfi);

    const auto fluid_props = fluid.props.data<run_on>();

    ParallelFor(bx, [T,h,X,fluid_props,flags]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int cell_is_covered = static_cast<int>(flags(i,j,k).isCovered());

      h(i,j,k) = fluid_props.enthalpy(IntVect(i,j,k), T, X, cell_is_covered);

    });

  }
}


namespace init_fluid_aux {

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Purpose: Set velocity initial conditions.                           !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_vel (const Box& sbx,
                 const Box& domain,
                 const Real dx,
                 const Real dy,
                 const Real dz,
                 const GpuArray<Real, 3>& plo,
                 MFIXInitialConditions& initial_conditions,
                 Array4<Real      > const& a_vel)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  // Set the initial conditions.
  for(int icv(0); icv < initial_conditions.ic().size(); ++icv)
  {

    int i_w(0), j_s(0), k_b(0);
    int i_e(0), j_n(0), k_t(0);

    calc_cell_ic(dx, dy, dz,
                 initial_conditions.ic(icv).region->lo(),
                 initial_conditions.ic(icv).region->hi(),
                 plo.data(),
                 i_w, i_e, j_s, j_n, k_b, k_t);

    // Use the volume fraction already calculated from particle data
    const RealVect ic_vel = initial_conditions.ic(icv).fluid.get_velocity();

    const int istart = amrex::max(slo[0], i_w);
    const int jstart = amrex::max(slo[1], j_s);
    const int kstart = amrex::max(slo[2], k_b);
    const int iend   = amrex::min(shi[0], i_e);
    const int jend   = amrex::min(shi[1], j_n);
    const int kend   = amrex::min(shi[2], k_t);

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      amrex::ParallelFor(box1, [a_vel, ugx=ic_vel[0]]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { a_vel(i,j,k,0) = ugx; });

      if(slo[0] < domlo[0] && domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [a_vel, ugx=ic_vel[0]]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { a_vel(i,j,k,0) = ugx; });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend)
      {
        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);
        amrex::ParallelFor(box3, [a_vel, ugx=ic_vel[0]]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { a_vel(i,j,k,0) = ugx; });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      amrex::ParallelFor(box1, [a_vel, vgx=ic_vel[1]]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { a_vel(i,j,k,1) = vgx; });

      if (slo[1] < domlo[1] && domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [a_vel, vgx=ic_vel[1]]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { a_vel(i,j,k,1) = vgx; });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend)
      {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);
        amrex::ParallelFor(box3, [a_vel, vgx=ic_vel[1]]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { a_vel(i,j,k,1) = vgx; });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);
      amrex::ParallelFor(box1, [a_vel, wgx=ic_vel[2]]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          { a_vel(i,j,k,2) = wgx; });

      if (slo[2] < domlo[2] && domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);
        amrex::ParallelFor(box2, [a_vel, wgx=ic_vel[2]]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { a_vel(i,j,k,2) = wgx; });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend)
      {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);
        amrex::ParallelFor(box3, [a_vel, wgx=ic_vel[2]]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { a_vel(i,j,k,2) = wgx; });
      }
    }
  }
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Purpose: Set fluid temperature initial conditions.                  !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_energy (const Box& sbx,
                  const Box& domain,
                  const Real dx,
                  const Real dy,
                  const Real dz,
                  const GpuArray<Real, 3>& plo,
                  MFIXInitialConditions& initial_conditions,
                  Array4<EBCellFlag const> const& a_flags,
                  Array4<Real      > const& a_T,
                  Array4<Real      > const& a_h,
                  Array4<Real const> const& a_X,
                  MFIXFluidPhase& fluid)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  const int solve_enthalpy = fluid.solve_enthalpy();

  const auto fluid_props = fluid.props.data<run_on>();

  // Set the initial conditions.
  for(int icv(0); icv < initial_conditions.ic().size(); ++icv) {

    int i_w(0), j_s(0), k_b(0);
    int i_e(0), j_n(0), k_t(0);

    calc_cell_ic(dx, dy, dz,
                 initial_conditions.ic(icv).region->lo(),
                 initial_conditions.ic(icv).region->hi(),
                 plo.data(),
                 i_w, i_e, j_s, j_n, k_b, k_t);

    const Real temperature = initial_conditions.ic(icv).fluid.get_temperature();

    const int istart = amrex::max(slo[0], i_w);
    const int jstart = amrex::max(slo[1], j_s);
    const int kstart = amrex::max(slo[2], k_b);
    const int iend   = amrex::min(shi[0], i_e);
    const int jend   = amrex::min(shi[1], j_n);
    const int kend   = amrex::min(shi[2], k_t);

    // Define the function to be used on the different Box-es
    auto set_quantities = [a_T,a_h, a_X,temperature,fluid_props,a_flags,solve_enthalpy]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      const int cell_is_covered = static_cast<int>(a_flags(i,j,k).isCovered());

      a_T(i,j,k) = temperature;

      if (solve_enthalpy) {
        a_h(i,j,k) = fluid_props.enthalpy(temperature, IntVect(i,j,k), a_X, cell_is_covered);
      }
    };

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_quantities(i,j,k); });

      if(slo[0] < domlo[0] && domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend)
      {
        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_quantities(i,j,k); });

      if (slo[1] < domlo[1] && domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend)
      {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);
      ParallelFor(box1, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_quantities(i,j,k); });

      if (slo[2] < domlo[2] && domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend)
      {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_quantities] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_quantities(i,j,k); });
      }
    }
  }
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Purpose: Set fluid species mass fractions initial conditions.       !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_species_g (const Box& sbx,
                       const Box& domain,
                       const Real dx,
                       const Real dy,
                       const Real dz,
                       const GpuArray<Real, 3>& plo,
                       MFIXInitialConditions& initial_conditions,
                       Array4<Real      > const& a_X)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  const int nspecies = a_X.nComp();

  // Set the initial conditions.
  for(int icv(0); icv < initial_conditions.ic().size(); ++icv)
  {
    int i_w(0), j_s(0), k_b(0);
    int i_e(0), j_n(0), k_t(0);

    calc_cell_ic(dx, dy, dz,
                 initial_conditions.ic(icv).region->lo(), initial_conditions.ic(icv).region->hi(),
                 plo.data(),
                 i_w, i_e, j_s, j_n, k_b, k_t);

    // Get the initial condition values
    Gpu::DeviceVector< Real> mass_fractions_d(nspecies);
    Gpu::HostVector  < Real> mass_fractions_h(nspecies);
    for (int n(0); n < nspecies; n++) {
      mass_fractions_h[n] = initial_conditions.ic(icv).fluid.get_species(n);
    }
    Gpu::copy(Gpu::hostToDevice, mass_fractions_h.begin(), mass_fractions_h.end(),
              mass_fractions_d.begin());

    Real* p_mass_fractions = mass_fractions_d.data();

    const int istart = std::max(slo[0], i_w);
    const int jstart = std::max(slo[1], j_s);
    const int kstart = std::max(slo[2], k_b);
    const int iend   = std::min(shi[0], i_e);
    const int jend   = std::min(shi[1], j_n);
    const int kend   = std::min(shi[2], k_t);

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, nspecies, [a_X,p_mass_fractions]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { a_X(i,j,k,n) = p_mass_fractions[n]; });

      if(slo[0] < domlo[0] && domlo[0] == istart)
      {
        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, nspecies, [a_X,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { a_X(i,j,k,n) = p_mass_fractions[n]; });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend)
      {
        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, nspecies, [a_X,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { a_X(i,j,k,n) = p_mass_fractions[n]; });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, nspecies, [a_X,p_mass_fractions]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { a_X(i,j,k,n) = p_mass_fractions[n]; });

      if (slo[1] < domlo[1] && domlo[1] == jstart)
      {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);
        ParallelFor(box2, nspecies, [a_X,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { a_X(i,j,k,n) = p_mass_fractions[n]; });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend)
      {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);
        ParallelFor(box3, nspecies, [a_X,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { a_X(i,j,k,n) = p_mass_fractions[n]; });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);
      ParallelFor(box1, nspecies, [a_X,p_mass_fractions]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { a_X(i,j,k,n) = p_mass_fractions[n]; });

      if (slo[2] < domlo[2] && domlo[2] == kstart)
      {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, nspecies, [a_X,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { a_X(i,j,k,n) = p_mass_fractions[n]; });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend)
      {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);

        ParallelFor(box3, nspecies, [a_X,p_mass_fractions]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        { a_X(i,j,k,n) = p_mass_fractions[n]; });
      }
    }

    Gpu::synchronize();

  }
}

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Purpose: Set fluid density                                          !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_rho ( const Box& sbx,
                  const Box& domain,
                  const Real dx,
                  const Real dy,
                  const Real dz,
                  const GpuArray<Real, 3>& plo,
                  MFIXInitialConditions& initial_conditions,
                  Array4<Real> const& a_rho)
{
  const IntVect slo(sbx.loVect());
  const IntVect shi(sbx.hiVect());

  const IntVect domlo(domain.loVect());
  const IntVect domhi(domain.hiVect());

  // Set the initial conditions.
  for(int icv(0); icv < initial_conditions.ic().size(); ++icv)
  {
    int i_w(0), j_s(0), k_b(0);
    int i_e(0), j_n(0), k_t(0);

    calc_cell_ic(dx, dy, dz,
                 initial_conditions.ic(icv).region->lo(), initial_conditions.ic(icv).region->hi(),
                 plo.data(),
                 i_w, i_e, j_s, j_n, k_b, k_t);

    const int istart = std::max(slo[0], i_w);
    const int jstart = std::max(slo[1], j_s);
    const int kstart = std::max(slo[2], k_b);
    const int iend   = std::min(shi[0], i_e);
    const int jend   = std::min(shi[1], j_n);
    const int kend   = std::min(shi[2], k_t);

    const Real density = initial_conditions.ic(icv).fluid.get_density();

    // Define the function
    auto set_density = [a_rho,density]
      AMREX_GPU_DEVICE (int i, int j, int k) -> void
    { a_rho(i,j,k) = density; };

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_density(i,j,k); });

      if(slo[0] < domlo[0] && domlo[0] == istart) {

        const IntVect low2(slo[0], jstart, kstart), hi2(istart-1, jend, kend);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }

      if(shi[0] > domhi[0] && domhi[0] == iend) {

        const IntVect low3(iend+1, jstart, kstart), hi3(shi[0], jend, kend);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_density(i,j,k); });

      if (slo[1] < domlo[1] && domlo[1] == jstart) {
        const IntVect low2(istart, slo[1], kstart), hi2(iend, jstart-1, kend);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }

      if (shi[1] > domhi[1] && domhi[1] == jend) {
        const IntVect low3(istart, jend+1, kstart), hi3(iend, shi[1], kend);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }
    }

    {
      const IntVect low1(istart, jstart, kstart), hi1(iend, jend, kend);
      const Box box1(low1, hi1);

      ParallelFor(box1, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { set_density(i,j,k); });

      if (slo[2] < domlo[2] && domlo[2] == kstart) {
        const IntVect low2(istart, jstart, slo[2]), hi2(iend, jend, kstart-1);
        const Box box2(low2, hi2);

        ParallelFor(box2, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }

      if (shi[2] > domhi[2] && domhi[2] == kend) {
        const IntVect low3(istart, jstart, kend+1), hi3(iend, jend, shi[2]);
        const Box box3(low3, hi3);

        ParallelFor(box3, [set_density] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { set_density(i,j,k); });
      }
    }
  }
}

} // end namespace init_fluid_aux
