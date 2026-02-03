#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_fluid.H>
#include <mfix_run_on.H>

#include <test_fluid_molecular_viscosity.H>

using namespace amrex;

///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                            Sutherland's Law                         ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_molecular_viscosity::
sutherland ()
{
  ParmParse pp;

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   0);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXSpecies species;
  MFIXReactions rxns;

  MFIXFluidPhase fluid;

  species.Initialize();
  rxns.Initialize(species);

  // Fails because molecular viscosity model is not defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("SuThErLaNd"));

  // Sutherland molecular viscosity model required inputs:
  //   "viscosity.molecular.Sutherland.T_ref"
  //   "viscosity.molecular.Sutherland.mu_ref"
  //   "viscosity.molecular.Sutherland.S"

  // Fails because no required inputs are defined.
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular.Sutherland.T_ref", 273.);

  // Fails because mu_ref and S are missing
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular.Sutherland.mu_ref", 1.716e-5);

  // Fails because S is missing
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular.Sutherland.S", 110.);

  // Sutherland is fully defined.
  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Sanity checks
  if (!fluid.isInitialized() ) { return FAIL; }
  if ( fluid.props.molViscosity.m_model != MolViscosityModel::Sutherland ) { return FAIL; }

  const auto props = fluid.props.data<run_on>();

  { Box bx({0,0,0}, {3,3,3});
    FArrayBox T_fab(bx, 1, The_Async_Arena());
    Array4<Real> T = T_fab.array();

    ParallelFor(bx, [bx,T]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      Real const Tlow(175.), Tlen(1325.);
      Real const dT( Tlen / static_cast<Real>(bx.numPts()-1));

      int const ijk(i+j*bx.length(0)+k*bx.length(0)*bx.length(1));
      T(i,j,k) = Tlow + dT*static_cast<Real>(ijk);
    });

    ReduceOps<ReduceOpLogicalOr> reduce_op;
    ReduceData<int> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    reduce_op.eval(bx, reduce_data, [lPASS=PASS, lFAIL=FAIL, T, props]
    AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      GpuArray<Real,1> interp{T(i,j,k)};
      Real visc0 = props.molViscosity(interp.data(), 0, 0);
      Real visc1 = props.molViscosity(i,j,k,T, Array4<Real const>());

      return (visc0 == visc1) ? lPASS : lFAIL;
    });
    ReduceTuple host_tuple = reduce_data.value(reduce_op);
    if (get<0>(host_tuple) == FAIL) { return FAIL; }
  }


  /* Compare against values archived values which compare reasonably well
     with reported values:
     Shpilrain, E. E. (2011). Air (properties of). THERMOPEDIA.
     https://www.thermopedia.com/content/553/

     T :  100    200    300    400.   500.   600.   700.   800.   900.  1000.  1100.  1200.  1300.
     visc: 71.1  132.5  184.6  230.1  270.1  305.8  338.8  369.8  398.1  424.4  449.0  473.0  496.0
  */
  { int constexpr count(13);

    Array<Real const, count> visc0 = {
      69.3829550, 132.9399119, 184.6588439, 228.5556164, 267.0528754,
     301.6064423, 333.1456355, 362.2977818, 389.5062918, 415.0965350,
     439.3147089, 462.3519146, 484.3597013
    };

    Real T(100.);
    for (int n(0); n<count; ++n) {
      Real visc = props.molViscosity(T, IntVect(), Array4<Real const>())*1.e7;
      Real err = std::abs(visc0[n]-visc);
      if (err > 1.e-7) { return FAIL; }
      T+= 100.;
    }
  }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");
  pp.remove("fluid0.viscosity.molecular.Sutherland.T_ref");
  pp.remove("fluid0.viscosity.molecular.Sutherland.mu_ref");
  pp.remove("fluid0.viscosity.molecular.Sutherland.S");

  return PASS;
}



///////////////////////////////////////////////////////////////////////////
////                                                                     //
///                            Reid four parameter                      ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_molecular_viscosity::
reid ()
{
  ParmParse pp;

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   0);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXSpecies species;
  MFIXReactions rxns;

  species.Initialize();
  rxns.Initialize(species);

  { MFIXFluidPhase fluid;

    // Fails because molecular viscosity model is not defined
    if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
    pp.add("fluid0.viscosity.molecular", std::string("rEiD"));

    // Reid molecular viscosity model required inputs:
    //   "viscosity.molecular.Reid.A"
    //   "viscosity.molecular.Reid.B"
    //   "viscosity.molecular.Reid.C"
    //   "viscosity.molecular.Reid.D"

    // Fails because no required inputs are defined.
    if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
    pp.add("fluid0.viscosity.molecular.Reid.A",5.76e-6);

    // Fails because B, C and D are missing
    if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
    pp.add("fluid0.viscosity.molecular.Reid.B", 1492.836337);

    // Fails because C and D are missing
    if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
    pp.add("fluid0.viscosity.molecular.Reid.C", 0.);

    // Fails because D is missing
    if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
    pp.add("fluid0.viscosity.molecular.Reid.D", 0.);

    // Reid is fully defined.
    if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

    // Sanity checks
    if (!fluid.isInitialized() ) { return FAIL; }
    if ( fluid.props.molViscosity.m_model != MolViscosityModel::Reid ) { return FAIL; }

    const auto props = fluid.props.data<run_on>();

    { Box bx({0,0,0}, {3,3,3});
      FArrayBox T_fab(bx, 1, The_Async_Arena());
      Array4<Real> T = T_fab.array();

      ParallelFor(bx, [bx,T]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Real const Tlow(275.), Tlen(350.);
        Real const dT( Tlen / static_cast<Real>(bx.numPts()-1));

        int const ijk(i+j*bx.length(0)+k*bx.length(0)*bx.length(1));
        T(i,j,k) = Tlow + dT*static_cast<Real>(ijk);
      });

      ReduceOps<ReduceOpSum> reduce_op;
      ReduceData<int> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      reduce_op.eval(bx, reduce_data, [lPASS=PASS, lFAIL=FAIL, T, props]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
        GpuArray<Real,1> interp{T(i,j,k)};
        Real visc0 = props.molViscosity(interp.data(), 0, 0);
        Real visc1 = props.molViscosity(i,j,k, T, Array4<Real const>());

        return (visc0 == visc1) ? lPASS : lFAIL;
      });
      ReduceTuple host_tuple = reduce_data.value(reduce_op);
      if (get<0>(host_tuple) == FAIL) { return FAIL; }

    }

    /* Compare against values archived values which compare reasonably well
       with reported values for water: Engineering Toolbox

       T (K)
       273.16   283.15   293.15   298.15   303.15   313.15   323.15   333.15   343.15
       353.15   363.15   373.15   383.15   393.15   413.15   433.15   453.15   473.15
       493.15   513.15   533.15   553.15   573.15   593.15   613.15   633.15

       visc [Pa s] (x10^6)
          1791.4   1306.0   1001.6   890.0   797.2   652.7   546.5   466.0   403.5
           354.0    314.2    281.6   254.6   232.0   196.6   170.4   150.4   134.6
           121.8    111.1    101.8    93.6    85.9    78.3    70.3    60.3
    */
    { int constexpr count(26);
      Array<Real const, count> visc0 = {
        1361.03331226, 1122.35345451, 937.61086644, 860.86274934, 792.62660186, 677.28899591,
         584.39460994,  508.72674272, 446.45003907, 394.70527024, 351.33338727, 314.68458864,
         283.48424419,  256.73744210, 213.62902121, 180.80200293, 155.28949455, 135.10322371,
         118.87613030,  105.64656130,  94.72406793,  85.60371203,  77.91002656,  71.35956797,
          65.73539564,   60.86936433
      };

      //std::printf("\n\n");
      Real T(273.16);
      for (int n(0); n<count; ++n) {
        Real visc = props.molViscosity(T, IntVect(), Array4<Real const>())*1.e6;
        Real err = std::abs(visc0[n]-visc);
        //std::printf("T: %6.2f   %12.8f   %e\n", T, visc, err);
        if (err > 1.e-8) { return FAIL; }
        T+= ((n==0) ? 9.99 : ((n==2 || n==3) ? 5. : (n>12 ? 20. : 10.)));
      }
      //std::printf("\n\n");
    }
  }

  pp.remove("fluid0.viscosity.molecular.Reid.A");
  pp.remove("fluid0.viscosity.molecular.Reid.B");
  pp.remove("fluid0.viscosity.molecular.Reid.C");
  pp.remove("fluid0.viscosity.molecular.Reid.D");

  // Repeat with four coefficients fit for water.
  { MFIXFluidPhase fluid;

    // Reid molecular viscosity model required inputs:
    //   "viscosity.molecular.Reid.A"
    //   "viscosity.molecular.Reid.B"
    //   "viscosity.molecular.Reid.C"
    //   "viscosity.molecular.Reid.D"

    // Fails because no required inputs are defined.
    pp.add("fluid0.viscosity.molecular.Reid.A", 1.856e-14);
    pp.add("fluid0.viscosity.molecular.Reid.B", 4209.);
    pp.add("fluid0.viscosity.molecular.Reid.C", 0.04527);
    pp.add("fluid0.viscosity.molecular.Reid.D", -3.376e-5);

    // Reid is fully defined -- perform basic sanity checks
    if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }
    if (!fluid.isInitialized() ) { return FAIL; }
    if ( fluid.props.molViscosity.m_model != MolViscosityModel::Reid ) { return FAIL; }

    const auto props = fluid.props.data<run_on>();

    // Again compare against values archived values.

    { int constexpr count(26);
      Array<Real const, count> visc0 = {
        1725.37324002, 1305.29501235, 1017.64882449, 907.78414120, 814.86115894, 667.99599184,
         559.04796744,  476.47060824,  412.65482653, 362.46062883, 322.33671206, 289.77828662,
         262.98494412,  240.63988273,  205.62076571, 179.39975930, 158.81085998, 141.89113167,
         127.38879750,  114.49619167,  102.69757782,  91.67755795,  81.26208670,  71.37725775,
          62.01799786,   53.22268838
      };

      //std::printf("\n\n");
      Real T(273.16);
      for (int n(0); n<count; ++n) {
        Real visc = props.molViscosity(T, IntVect(), Array4<Real const>())*1.e6;
        Real err = std::abs(visc0[n]-visc);
        //std::printf("T: %6.2f   %12.8f   %e\n", T, visc, err);
        if (err > 1.e-8) { return FAIL; }
        T+= ((n==0) ? 9.99 : ((n==2 || n==3) ? 5. : (n>12 ? 20. : 10.)));
      }
    }
    //std::printf("\n\n");
  }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.viscosity.molecular");

  pp.remove("fluid0.viscosity.molecular.Reid.A");
  pp.remove("fluid0.viscosity.molecular.Reid.B");
  pp.remove("fluid0.viscosity.molecular.Reid.C");
  pp.remove("fluid0.viscosity.molecular.Reid.D");

  return PASS;
}

///////////////////////////////////////////////////////////////////////////
////                                                                     //
///           Mixture of same species using Sutherland                  ///
//                                                                     ////
///////////////////////////////////////////////////////////////////////////
int test_fluid_molecular_viscosity::
mixture_sutherland ()
{
  ParmParse pp;

  pp.add("mfix.advect_density",  0);
  pp.add("mfix.advect_enthalpy", 0);
  pp.add("mfix.solve_species",   1);
  pp.add("mfix.advect_tracer",   0);

  pp.add("fluid.solve", std::string("fluid0"));

  MFIXSpecies species;
  MFIXReactions rxns;

  // Mixture inputs
  pp.addarr("fluid0.species", std::vector<std::string>{"A", "B", "C"});
  pp.addarr("species.solve", std::vector<std::string>{"A", "B", "C"});
  pp.add("species.viscosity.molecular", std::string("SuThERLanD"));
  pp.add("species.diffusivity", std::string("constant"));
  pp.add("species.diffusivity.constant", 1.9e-5);

  // Sutherland mixture molecular viscosity model required inputs:
  pp.add("species.A.molecular_weight", 28.96e-3);
  pp.add("species.A.viscosity.molecular.Sutherland.T_ref", 273.);
  pp.add("species.A.viscosity.molecular.Sutherland.mu_ref", 1.716e-5);
  pp.add("species.A.viscosity.molecular.Sutherland.S", 110.);

  pp.add("species.B.molecular_weight", 28.96e-3);
  pp.add("species.B.viscosity.molecular.Sutherland.T_ref", 273.);
  pp.add("species.B.viscosity.molecular.Sutherland.mu_ref", 1.716e-5);
  pp.add("species.B.viscosity.molecular.Sutherland.S", 110.);

  pp.add("species.C.molecular_weight", 28.96e-3);
  pp.add("species.C.viscosity.molecular.Sutherland.T_ref", 273.);
  pp.add("species.C.viscosity.molecular.Sutherland.mu_ref", 1.716e-5);
  pp.add("species.C.viscosity.molecular.Sutherland.S", 110.);

  species.Initialize();
  rxns.Initialize(species);

  MFIXFluidPhase fluid;


  fluid.Initialize(species, rxns);
  // Fails because no molecular viscosity model is defined
  if ( fluid.Initialize(species, rxns) == PASS ) { return FAIL; }
  pp.add("fluid0.viscosity.molecular", std::string("mIXturE"));

  if ( fluid.Initialize(species, rxns) == FAIL ) { return FAIL; }

  // Sanity checks
  if (!fluid.isInitialized() ) { return FAIL; }
  if ( fluid.props.molViscosity.m_model != MolViscosityModel::Sutherland ) { return FAIL; }

  const auto props = fluid.props.data<run_on>();

  { Box bx({0,0,0}, {3,3,3});
    FArrayBox T_fab(bx, 1, The_Async_Arena());
    FArrayBox Xk_fab(bx, 3, The_Async_Arena());
    Array4<Real> T = T_fab.array();
    Array4<Real> Xk = Xk_fab.array();

    ParallelFor(bx, [bx,T,Xk]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      Real const Tlow(175.), Tlen(1325.);
      Real const dT( Tlen / static_cast<Real>(bx.numPts()-1));

      int const ijk(i+j*bx.length(0)+k*bx.length(0)*bx.length(1));
      T(i,j,k) = Tlow + dT*static_cast<Real>(ijk);

      Real const Xk_lo(0.);

      amrex::Real Xk_a = Xk_lo + dT*static_cast<Real>(ijk);
      Xk(i,j,k,0) = Xk_a;
      Xk(i,j,k,1) = 0.5*(1. - Xk_a);
      Xk(i,j,k,2) = 0.5*(1. - Xk_a);
    });

    ReduceOps<ReduceOpLogicalOr> reduce_op;
    ReduceData<int> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    reduce_op.eval(bx, reduce_data, [lPASS=PASS, lFAIL=FAIL, T, Xk, props]
    AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      GpuArray<Real,4> interp{T(i,j,k), Xk(i,j,k,0), Xk(i,j,k,1), Xk(i,j,k,2)};
      Real visc0 = props.molViscosity(interp.data(), 0, 1);
      Real visc1 = props.molViscosity(i,j,k,T, Xk);

      return (visc0 == visc1) ? lPASS : lFAIL;
    });
    ReduceTuple host_tuple = reduce_data.value(reduce_op);
    if (get<0>(host_tuple) == FAIL) { return FAIL; }
  }

  /* Compare against values archived values which compare reasonably well
     with reported values:
     Shpilrain, E. E. (2011). Air (properties of). THERMOPEDIA.
     https://www.thermopedia.com/content/553/

     T :  100    200    300    400.   500.   600.   700.   800.   900.  1000.  1100.  1200.  1300.
     visc: 71.1  132.5  184.6  230.1  270.1  305.8  338.8  369.8  398.1  424.4  449.0  473.0  496.0
  */
  { int constexpr count(13);

    Array<Real const, count> visc0 = {
      69.3829550, 132.9399119, 184.6588439, 228.5556164, 267.0528754,
     301.6064423, 333.1456355, 362.2977818, 389.5062918, 415.0965350,
     439.3147089, 462.3519146, 484.3597013
    };

    Real T(100.);
    for (int n(0); n<count; ++n) {
      Box bx({0,0,0}, {0,0,0});
      FArrayBox Xk_fab(bx, 3, The_Async_Arena());
      Array4<Real> Xk = Xk_fab.array();

      ParallelFor(bx, [bx,Xk,n]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        Xk(i,j,k,0) = static_cast<Real>(n/(count-1));
        Xk(i,j,k,1) = 0.5*(1. - static_cast<Real>(n/(count-1)));
        Xk(i,j,k,2) = 0.5*(1. - static_cast<Real>(n/(count-1)));
      });

      Real visc = props.molViscosity(T, IntVect(0,0,0), Xk)*1.e7;
      Real err = std::abs(visc0[n]-visc);
      if (err > 1.e-7) { return FAIL; }
      T+= 100.;
    }
  }

  pp.remove("mfix.advect_density");
  pp.remove("mfix.advect_enthalpy");
  pp.remove("mfix.solve_species");
  pp.remove("mfix.advect_tracer");
  pp.remove("fluid.solve");
  pp.remove("fluid0.species");
  pp.remove("species.solve");
  pp.remove("species.viscosity.molecular");
  pp.remove("species.diffusivity");
  pp.remove("species.diffusivity.constant");
  pp.remove("species.A.molecular_weight");
  pp.remove("species.B.molecular_weight");
  pp.remove("species.A.viscosity.molecular.Sutherland.T_ref");
  pp.remove("species.A.viscosity.molecular.Sutherland.mu_ref");
  pp.remove("species.A.viscosity.molecular.Sutherland.S");
  pp.remove("species.B.viscosity.molecular.Sutherland.T_ref");
  pp.remove("species.B.viscosity.molecular.Sutherland.mu_ref");
  pp.remove("species.B.viscosity.molecular.Sutherland.S");
  pp.remove("fluid0.viscosity.molecular");

  return PASS;
}
