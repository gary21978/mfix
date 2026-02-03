#include <mfix_distributions.H>
#include <mfix_distributions_normal.H>
#include <mfix_distributions_log_normal.H>
#include <mfix_distributions_uniform.H>
#include <mfix_distributions_custom.H>

#include <test_constant_dist.H>
#include <test_normal_dist.H>
#include <test_log_normal_dist.H>
#include <test_uniform_dist.H>

#include <test_custom_dist.H>

int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, true);
  {
    amrex::ParmParse pp("amrex");

    pp.add("fpe_trap_invalid",  1);
    pp.add("fpe_trap_zero",     1);
    pp.add("fpe_trap_overflow", 1);

    static_assert(std::is_trivially_copyable_v<INPUT_DIST_t> == true);

    static_assert(std::is_trivially_copyable_v<distributions::normal> == true);
    static_assert(std::is_trivially_copyable_v<distributions::log_normal> == true);
    static_assert(std::is_trivially_copyable_v<distributions::uniform> == true);
    static_assert(std::is_trivially_copyable_v<distributions::custom> == true);

    test_constant_dist constant_dist;
    if (constant_dist.test() == PASS) {
      amrex::Print() << "constant distribution tests PASS!\n";
    } else {
      amrex::Print() << "\nconstant distribution tests FAILED!\n";
      return FAIL;
    }

    test_normal_dist normal_dist;
    if (normal_dist.test() == PASS) {
      amrex::Print() << "normal distribution tests PASS!\n";
    } else {
      amrex::Print() << "\nnormal distribution tests FAILED!\n";
      return FAIL;
    }

    test_log_normal_dist log_normal_dist;
    if (log_normal_dist.test() == PASS) {
      amrex::Print() << "log-normal distribution tests PASS!\n";
    } else {
      amrex::Print() << "\nlog-normal distribution tests FAILED!\n";
      return FAIL;
    }

    test_uniform_dist uniform_dist;
    if (uniform_dist.test() == PASS) {
      amrex::Print() << "uniform distribution tests PASS!\n";
    } else {
      amrex::Print() << "\nuniform distribution tests FAILED!\n";
      return FAIL;
    }

    test_custom_dist custom_dist;
    if (custom_dist.test() == PASS) {
      amrex::Print() << "custom distribution tests PASS!\n";
    } else {
      amrex::Print() << "\ncustom distribution tests FAILED!\n";
      return FAIL;
    }

  }
  amrex::Finalize();
}

