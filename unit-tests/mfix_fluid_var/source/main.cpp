#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include <mfix_fluid_var.H>

using namespace amrex;

int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, true);
  {

    int const PASS(0);
    int const FAIL(1);

    ///////////////////////////////////////////////////////////////////////////
    ////                                                                     //
    ///                       Tests for scalar inputs                       ///
    //                                                                     ////
    ///////////////////////////////////////////////////////////////////////////

    { ParmParse pp("test_000");

      // Constant scalar format
      pp.addarr("scalar_1D", Vector<Real>{1.0});
      FVAR_ scalar_1D({1});
      if (scalar_1D.get(pp, "scalar_1D") == FAIL) { return FAIL; }

      // Transient scalar format
      pp.addarr("scalar_2D", Vector<Real>{1.0, 2.0});
      FVAR_ scalar_2D({1});
      if (scalar_2D.get(pp, "scalar_2D") == FAIL) { return FAIL; }

      // Invalid format
      pp.addarr("scalar_0D", Vector<Real>{});
      FVAR_ scalar_0D({1});
      if (scalar_0D.get(pp, "scalar_0D") == PASS) { return FAIL; }

      // Invalid format
      pp.addarr("scalar_3D", Vector<Real>{1.0, 2.0, 3.0});
      FVAR_ scalar_3D({1});
      if (scalar_3D.get(pp, "scalar_3D") == PASS) { return FAIL; }

      // Invalid format
      pp.addarr("scalar_5D", Vector<Real>{1., 2., 3., 4., 5.});
      FVAR_ scalar_5D({1});
      if (scalar_5D.get(pp, "scalar_5D") == PASS) { return FAIL; }

    } Print() << "test_000 PASSED!\n";


    // Test constant scalar inputs
    { ParmParse pp("test_001");

      FVAR_ scalar({1});

      // Fail if is_defined returns true before we call get()
      if (scalar.is_defined()) { return FAIL; }
      if (scalar.get(pp, "scalar") == FAIL) { return FAIL; }

      Real const val(100.);

      pp.add("scalar", val);

      if (scalar.get(pp, "scalar") == FAIL) { return FAIL; }
      if (scalar.get(pp, "scalar") == PASS) { return FAIL; }

      if (!scalar.is_defined())  { return FAIL; }

      if (!scalar.is_constant())   { return FAIL; }
      if ((scalar())[0]    != val) { return FAIL; }
      if ((scalar(0.0))[0] != val) { return FAIL; }
      if ((scalar(1.0))[0] != val) { return FAIL; }

    } Print() << "test_001 PASSED!\n";


    // Make sure we can manage duplicate scalar entries and if
    // there are duplicates, we return the "last" value.
    { ParmParse pp("test_002");

      FVAR_ scalar({1});

      pp.add("scalar", 100.);
      pp.add("scalar", 110.);
      pp.add("scalar", 120.);
      pp.add("scalar", 130.);
      pp.add("scalar", 140.);
      pp.add("scalar", 150.);

      if (scalar.get(pp, "scalar") == FAIL) { return 1; }

      if (!scalar.is_defined())  { return FAIL; }
      if (!scalar.is_constant()) { return FAIL; }
      if ((scalar())[0] != 150.) { return FAIL; }

    } Print() << "test_002 PASSED!\n";


    // Transient scalar entries
    { ParmParse pp("test_003");

      FVAR_ scalar({1});

      // Order should not matter
      pp.addarr("scalar", Vector<Real>{2.0, 300.});
      pp.addarr("scalar", Vector<Real>{0.0, 100.});
      pp.addarr("scalar", Vector<Real>{1.0, 200.});

      if (scalar.get(pp, "scalar") == FAIL) { return 1; }
      if (!scalar.is_defined())  { return FAIL; }
      if ( scalar.is_constant()) { return FAIL; }

      if ((scalar())[0]    != 100.) { return FAIL; }
      if ((scalar(0.0))[0] != 100.) { return FAIL; }
      if ((scalar(0.5))[0] != 150.) { return FAIL; }
      if ((scalar(1.0))[0] != 200.) { return FAIL; }
      if ((scalar(1.5))[0] != 250.) { return FAIL; }
      if ((scalar(2.0))[0] != 300.) { return FAIL; }
      if ((scalar(5.0))[0] != 300.) { return FAIL; }

    } Print() << "test_003 PASSED!\n";


    // If the first time entry for a transient is not zero,
    // make sure the scalar is constant until the first entry.
    { ParmParse pp("test_004");

      FVAR_ scalar({1});

      pp.addarr("scalar", Vector<Real>{1.0, 500.});
      pp.addarr("scalar", Vector<Real>{6.0, 400.});

      if (scalar.get(pp, "scalar") == FAIL) { return 1; }
      if (!scalar.is_defined())  { return FAIL; }
      if ( scalar.is_constant()) { return FAIL; }

      if ((scalar())[0]    != 500.) { return FAIL; }
      if ((scalar(0.0))[0] != 500.) { return FAIL; }
      if ((scalar(1.0))[0] != 500.) { return FAIL; }
      if ((scalar(3.5))[0] != 450.) { return FAIL; }
      if ((scalar(6.0))[0] != 400.) { return FAIL; }
      if ((scalar(9.0))[0] != 400.) { return FAIL; }

    } Print() << "test_004 PASSED!\n";


    // It is valid to only have one time entry
    { ParmParse pp("test_005");

      FVAR_ scalar({1});

      pp.addarr("scalar", Vector<Real>{100., -10.});

      if (scalar.get(pp, "scalar") == FAIL) { return 1; }
      if (!scalar.is_defined())  { return FAIL; }
      if ( scalar.is_constant()) { return FAIL; }

      if ((scalar())[0]      != -10.) { return FAIL; }
      if ((scalar(  0.0))[0] != -10.) { return FAIL; }
      if ((scalar(100.0))[0] != -10.) { return FAIL; }
      if ((scalar(200.0))[0] != -10.) { return FAIL; }

    } Print() << "test_005 PASSED!\n";


    // Mixed (constant / transient) entries are invalid
    { ParmParse pp("test_006");

      // Forward
      { FVAR_ scalar({1});

        pp.addarr("scalar", Vector<Real>{1.0, 500.});
        pp.addarr("scalar", Vector<Real>{400.});

        if (scalar.get(pp, "scalar") == PASS) { return 1; }
      }

      // Backward
      { FVAR_ scalar({1});

        pp.addarr("scalar", Vector<Real>{400.});
        pp.addarr("scalar", Vector<Real>{1.0, 500.});

        if (scalar.get(pp, "scalar") == PASS) { return 1; }
      }

      // CTC Sandwich
      { FVAR_ scalar({1});

        pp.addarr("scalar", Vector<Real>{400.});
        pp.addarr("scalar", Vector<Real>{1.0, 500.});
        pp.addarr("scalar", Vector<Real>{400.});

        if (scalar.get(pp, "scalar") == PASS) { return 1; }
      }

      // TCT Sandwich
      { FVAR_ scalar({1});

        pp.addarr("scalar", Vector<Real>{1.0, 500.});
        pp.addarr("scalar", Vector<Real>{400.});
        pp.addarr("scalar", Vector<Real>{1.0, 500.});

        if (scalar.get(pp, "scalar") == PASS) { return 1; }
      }

    } Print() << "test_006 PASSED!\n";


    // Duplicate time entries are invalid
    { ParmParse pp("test_007");

      FVAR_ scalar({1});

      pp.addarr("scalar", Vector<Real>{1.0, 500.});
      pp.addarr("scalar", Vector<Real>{1.0, 400.});

      if (scalar.get(pp, "scalar") == PASS) { return 1; }

    } Print() << "test_007 PASSED!\n";

    ///////////////////////////////////////////////////////////////////////////
    ////                                                                     //
    ///                       Tests for vector inputs                       ///
    //                                                                     ////
    ///////////////////////////////////////////////////////////////////////////

    // Test input formats for for 3D vector.
    { ParmParse pp("test_100");

      // Constant format
      pp.addarr("vector_3D", Vector<Real>{1.0, 2.0, 3.0});
      FVAR_ vector_3D({3});
      if (vector_3D.get(pp, "vector_3D") == FAIL) { return FAIL; }

      // Transient format
      pp.addarr("vector_4D", Vector<Real>{1.0, 2.0, 3.0, 4.0});
      FVAR_ vector_4D({3});
      if (vector_4D.get(pp, "vector_3D") == FAIL) { return FAIL; }

      // Invalid format
      pp.addarr("vector_0D", Vector<Real>{});
      FVAR_ vector_0D({3});
      if (vector_0D.get(pp, "vector_0D") == PASS) { return FAIL; }

      // Invalid format
      pp.addarr("vector_1D", Vector<Real>{1.0});
      FVAR_ vector_1D({3});
      if (vector_1D.get(pp, "vector_1D") == PASS) { return FAIL; }

      // Invalid format
      pp.addarr("vector_2D", Vector<Real>{1.0, 2.0});
      FVAR_ vector_2D({3});
      if (vector_2D.get(pp, "vector_2D") == PASS) { return FAIL; }

      // Invalid format
      pp.addarr("vector_5D", Vector<Real>{1., 2., 3., 4., 5.});
      FVAR_ vector_5D({3});
      if (vector_5D.get(pp, "vector_5D") == PASS) { return FAIL; }

    } Print() << "test_100 PASSED!\n";


    // Test constant vector inputs
    { ParmParse pp("test_101");

      FVAR_ vector({3});

      if (vector.is_defined()) { return FAIL; }
      if (vector.get(pp, "vector") == FAIL) { return FAIL; }

      Vector<Real> const val{100., 200., -300.};

      pp.addarr("vector", val);

      if (vector.get(pp, "vector") == FAIL) { return FAIL; }
      if (vector.get(pp, "vector") == PASS) { return FAIL; }

      if (!vector.is_defined())  { return FAIL; }

      if (!vector.is_constant())   { return FAIL; }
      if ((vector())[0]    != val[0]) { return FAIL; }
      if ((vector())[1]    != val[1]) { return FAIL; }
      if ((vector())[2]    != val[2]) { return FAIL; }
      if ((vector(0.0))[0] != val[0]) { return FAIL; }
      if ((vector(0.0))[1] != val[1]) { return FAIL; }
      if ((vector(0.0))[2] != val[2]) { return FAIL; }
      if ((vector(1.0))[0] != val[0]) { return FAIL; }
      if ((vector(1.0))[1] != val[1]) { return FAIL; }
      if ((vector(1.0))[2] != val[2]) { return FAIL; }

    } Print() << "test_101 PASSED!\n";


    // Make sure we can manage duplicate entries and if
    // there are duplicates, we return the "last" value.
    { ParmParse pp("test_102");

      FVAR_ vector({3});

      if (vector.is_defined()) { return FAIL; }
      if (vector.get(pp, "vector") == FAIL) { return FAIL; }

      pp.addarr("vector", Vector<Real>{10., 200., 3000.});
      pp.addarr("vector", Vector<Real>{11., 210., 3100.});
      pp.addarr("vector", Vector<Real>{12., 220., 3200.});
      pp.addarr("vector", Vector<Real>{13., 230., 3300.});
      pp.addarr("vector", Vector<Real>{14., 240., 3400.});
      pp.addarr("vector", Vector<Real>{15., 250., 3500.});

      if (vector.get(pp, "vector") == FAIL) { return FAIL; }
      if (vector.get(pp, "vector") == PASS) { return FAIL; }

      if (!vector.is_defined())  { return FAIL; }

      if (!vector.is_constant())   { return FAIL; }
      if ((vector())[0]  !=   15.) { return FAIL; }
      if ((vector())[1]  !=  250.) { return FAIL; }
      if ((vector())[2]  != 3500.) { return FAIL; }

      if ((vector(5.0))[0]  !=   15.) { return FAIL; }
      if ((vector(5.0))[1]  !=  250.) { return FAIL; }
      if ((vector(5.0))[2]  != 3500.) { return FAIL; }

    } Print() << "test_102 PASSED!\n";


    // Transient vector entries
    { ParmParse pp("test_103");

      FVAR_ vector({3});

      if (vector.is_defined()) { return FAIL; }
      if (vector.get(pp, "vector") == FAIL) { return FAIL; }

      // Order should not matter
      pp.addarr("vector", Vector<Real>{ 2.0, 12., 220., 3200.});
      pp.addarr("vector", Vector<Real>{ 3.0, 13., 230., 3300.});
      pp.addarr("vector", Vector<Real>{ 1.0, 11., 210., 3100.});
      pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});
      pp.addarr("vector", Vector<Real>{ 4.0, 14., 240., 3400.});
      pp.addarr("vector", Vector<Real>{ 0.0, 10., 200., 3000.});

      if (vector.get(pp, "vector") == FAIL) { return FAIL; }
      if (vector.get(pp, "vector") == PASS) { return FAIL; }

      if (!vector.is_defined())  { return FAIL; }
      if ( vector.is_constant()) { return FAIL; }

      if ((vector())[0] !=   10.) {  return FAIL; }
      if ((vector())[1] !=  200.) {  return FAIL; }
      if ((vector())[2] != 3000.) {  return FAIL; }

      if ((vector(0.0))[0] !=   10.) { return FAIL; }
      if ((vector(0.0))[1] !=  200.) { return FAIL; }
      if ((vector(0.0))[2] != 3000.) { return FAIL; }

      if ((vector(15.0))[0] !=   15.) { return FAIL; }
      if ((vector(15.0))[1] !=  250.) { return FAIL; }
      if ((vector(15.0))[2] != 3500.) { return FAIL; }

      if ((vector(4.5))[0]  !=   14.5) { return FAIL; }
      if ((vector(4.5))[1]  !=  245.0) { return FAIL; }
      if ((vector(4.5))[2]  != 3450.0) { return FAIL; }

    } Print() << "test_103 PASSED!\n";


    // If the first time entry for a transient is not zero,
    // make sure the scalar is constant until the first entry.
    { ParmParse pp("test_104");

      FVAR_ vector({3});

      if (vector.is_defined()) { return FAIL; }
      if (vector.get(pp, "vector") == FAIL) { return FAIL; }

      // Order should not matter
      pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});
      pp.addarr("vector", Vector<Real>{ 3.0, 13., 230., 3300.});

      if (vector.get(pp, "vector") == FAIL) { return FAIL; }
      if (vector.get(pp, "vector") == PASS) { return FAIL; }

      if (!vector.is_defined())  { return FAIL; }
      if ( vector.is_constant()) { return FAIL; }

      if ((vector())[0] !=   13.) {  return FAIL; }
      if ((vector())[1] !=  230.) {  return FAIL; }
      if ((vector())[2] != 3300.) {  return FAIL; }

      if ((vector(3.0))[0] !=   13.) {  return FAIL; }
      if ((vector(3.0))[1] !=  230.) {  return FAIL; }
      if ((vector(3.0))[2] != 3300.) {  return FAIL; }

      if ((vector(4.0))[0] !=   14.) {  return FAIL; }
      if ((vector(4.0))[1] !=  240.) {  return FAIL; }
      if ((vector(4.0))[2] != 3400.) {  return FAIL; }

      if ((vector(5.0))[0] !=   15.) {  return FAIL; }
      if ((vector(5.0))[1] !=  250.) {  return FAIL; }
      if ((vector(5.0))[2] != 3500.) {  return FAIL; }

      if ((vector(50.))[0] !=   15.) {  return FAIL; }
      if ((vector(50.))[1] !=  250.) {  return FAIL; }
      if ((vector(50.))[2] != 3500.) {  return FAIL; }

    } Print() << "test_104 PASSED!\n";


    // It is valid to only have one time entry
    { ParmParse pp("test_105");

      FVAR_ vector({3});

      if (vector.is_defined()) { return FAIL; }
      if (vector.get(pp, "vector") == FAIL) { return FAIL; }

      pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});

      if (vector.get(pp, "vector") == FAIL) { return FAIL; }
      if (vector.get(pp, "vector") == PASS) { return FAIL; }

      if (!vector.is_defined())  { return FAIL; }
      if ( vector.is_constant()) { return FAIL; }

      if ((vector())[0] !=   15.) {  return FAIL; }
      if ((vector())[1] !=  250.) {  return FAIL; }
      if ((vector())[2] != 3500.) {  return FAIL; }

      if ((vector(0.5))[0] !=   15.) {  return FAIL; }
      if ((vector(0.5))[1] !=  250.) {  return FAIL; }
      if ((vector(0.5))[2] != 3500.) {  return FAIL; }

      if ((vector(5.0))[0] !=   15.) {  return FAIL; }
      if ((vector(5.0))[1] !=  250.) {  return FAIL; }
      if ((vector(5.0))[2] != 3500.) {  return FAIL; }

      if ((vector(50.))[0] !=   15.) {  return FAIL; }
      if ((vector(50.))[1] !=  250.) {  return FAIL; }
      if ((vector(50.))[2] != 3500.) {  return FAIL; }

    } Print() << "test_105 PASSED!\n";


    // Mixed (constant / transient) entries are invalid
    { ParmParse pp("test_106");

      // Forward
      { FVAR_ vector({3});

        pp.addarr("vector", Vector<Real>{ 15., 250., 3500.});
        pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});

        if (vector.get(pp, "vector") == PASS) { return FAIL; }
      }

      // Backward
      { FVAR_ vector({3});

        pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});
        pp.addarr("vector", Vector<Real>{ 15., 250., 3500.});

        if (vector.get(pp, "vector") == PASS) { return FAIL; }
      }

      // CTC Sandwich
      { FVAR_ vector({3});

        pp.addarr("vector", Vector<Real>{ 15., 250., 3500.});
        pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});
        pp.addarr("vector", Vector<Real>{ 15., 250., 3500.});

        if (vector.get(pp, "vector") == PASS) { return FAIL; }
      }

      // TCT Sandwich
      { FVAR_ vector({3});

        pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});
        pp.addarr("vector", Vector<Real>{ 15., 250., 3500.});
        pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});

        if (vector.get(pp, "vector") == PASS) { return FAIL; }
      }

    } Print() << "test_106 PASSED!\n";


    // Duplicate time entries are invalid
    { ParmParse pp("test_107");

      FVAR_ vector({3});

      pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});
      pp.addarr("vector", Vector<Real>{ 5.0, 15., 250., 3500.});

      if (vector.get(pp, "vector") == PASS) { return FAIL; }

    } Print() << "test_107 PASSED!\n";


    ///////////////////////////////////////////////////////////////////////////
    ////                                                                     //
    ///                        Tests for mixed inputs                       ///
    //                                                                     ////
    ///////////////////////////////////////////////////////////////////////////

    // Test input formats for for mixed inputs.
    { ParmParse pp("test_200");

      // Constant scalar
      FVAR_ scalar_const({1,3});
      pp.add("scalar_const", 1.0);
      if (scalar_const.get(pp, "scalar_const") == FAIL) { return FAIL; }

      // Transient scalar
      FVAR_ scalar_dt({1,3});
      pp.addarr("scalar_dt", Vector<Real>{1.0, 2.0});
      if (scalar_dt.get(pp, "scalar_dt") == FAIL) { return FAIL; }

      // Constant vector
      FVAR_ vector_const({1,3});
      pp.addarr("vector_const", Vector<Real>{1.,2.,3.});
      if (vector_const.get(pp, "vector_const") == FAIL) { return FAIL; }

      // Transient vector
      FVAR_ vector_dt({1,3});
      pp.addarr("vector_dt", Vector<Real>{0., 1.,2.,3.});
      if (vector_dt.get(pp, "vector_dt") == FAIL) { return FAIL; }

      // Invalid format
      pp.addarr("vector_0D", Vector<Real>{});
      FVAR_ vector_0D({1,3});
      if (vector_0D.get(pp, "vector_0D") == PASS) { return FAIL; }

      // Invalid format
      pp.addarr("vector_5D", Vector<Real>{1., 2., 3., 4., 5.});
      FVAR_ vector_5D({1,3});
      if (vector_5D.get(pp, "vector_5D") == PASS) { return FAIL; }

      // No mixing
      pp.addarr("vector_mixed1", Vector<Real>{1.});
      pp.addarr("vector_mixed1", Vector<Real>{2., 3., 4.});
      FVAR_ vector_mixed1({1,3});
      if (vector_mixed1.get(pp, "vector_mixed1") == PASS) { return FAIL; }

      // No mixing
      pp.addarr("vector_mixed2", Vector<Real>{1., 2.});
      pp.addarr("vector_mixed2", Vector<Real>{2., 3., 4.});
      FVAR_ vector_mixed2({1,3});
      if (vector_mixed2.get(pp, "vector_mixed2") == PASS) { return FAIL; }

      // No mixing
      pp.addarr("vector_mixed3", Vector<Real>{1.});
      pp.addarr("vector_mixed3", Vector<Real>{0., 2., 3., 4.});
      FVAR_ vector_mixed3({1,3});
      if (vector_mixed3.get(pp, "vector_mixed3") == PASS) { return FAIL; }

      // No mixing
      pp.addarr("vector_mixed4", Vector<Real>{3., 1.});
      pp.addarr("vector_mixed4", Vector<Real>{0., 2., 3., 4.});
      FVAR_ vector_mixed4({1,3});
      if (vector_mixed4.get(pp, "vector_mixed4") == PASS) { return FAIL; }

    } Print() << "test_200 PASSED!\n";


  }
  amrex::Finalize();
}



