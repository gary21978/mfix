from pathlib import Path
from unittest.mock import patch

from rtl2 import regtest

EXAMPLE_SH = """
mpiexec --mca opal_warn_on_missing_libcuda 0 -n 4 BASE_DIR/rt-MFIX-Exa/build/benchmarks/rtl2/fluidbed_cyl/fluidbed_cyl.ex BASE_DIR/subprojects/mfix/benchmarks/rtl2/fluidbed_cyl/inputs mfix.plot_file=fluidbed_cyl_plt mfix.check_file=fluidbed_cyl_chk mfix.checkpoint_files_output=0  mfix.max_step=3
mpiexec --mca opal_warn_on_missing_libcuda 0 -n 4 BASE_DIR/rt-MFIX-Exa/build/benchmarks/rtl2/fluidbed_periodic/fluidbed_periodic.ex BASE_DIR/subprojects/mfix/benchmarks/rtl2/fluidbed_periodic/inputs mfix.plot_file=fluidbed_periodic_plt mfix.check_file=fluidbed_periodic_chk mfix.checkpoint_files_output=0  mfix.max_step=3
mpiexec --mca opal_warn_on_missing_libcuda 0 -n 4 BASE_DIR/rt-MFIX-Exa/build/benchmarks/rtl2/fluidbed_sq/fluidbed_sq.ex BASE_DIR/subprojects/mfix/benchmarks/rtl2/fluidbed_sq/inputs mfix.plot_file=fluidbed_sq_plt mfix.check_file=fluidbed_sq_chk mfix.checkpoint_files_output=0  mfix.max_step=3
mpiexec --mca opal_warn_on_missing_libcuda 0 -n 4 BASE_DIR/rt-MFIX-Exa/build/benchmarks/rtl2/hopper/hopper.ex BASE_DIR/subprojects/mfix/benchmarks/rtl2/hopper/inputs mfix.plot_file=hopper_plt mfix.check_file=hopper_chk mfix.checkpoint_files_output=0  mfix.max_step=3
mpiexec --mca opal_warn_on_missing_libcuda 0 -n 4 BASE_DIR/rt-MFIX-Exa/build/benchmarks/rtl2/riser_cyl/riser_cyl.ex BASE_DIR/subprojects/mfix/benchmarks/rtl2/riser_cyl/inputs mfix.plot_file=riser_cyl_plt mfix.check_file=riser_cyl_chk mfix.checkpoint_files_output=0  mfix.max_step=3
mpiexec --mca opal_warn_on_missing_libcuda 0 -n 4 BASE_DIR/rt-MFIX-Exa/build/benchmarks/rtl2/riser_sq/riser_sq.ex BASE_DIR/subprojects/mfix/benchmarks/rtl2/riser_sq/inputs mfix.plot_file=riser_sq_plt mfix.check_file=riser_sq_chk mfix.checkpoint_files_output=0  mfix.max_step=3
"""


@patch("rtl2.params.Path.cwd")
def test_sbatch_gen(path, tmp_path):
    path.return_value = Path("BASE_DIR")
    test_ini = tmp_path / "test.ini"
    with open(test_ini, "w") as testf:
        with open("MFIX-tests.ini") as default_ini:
            testf.write(
                default_ini.read().replace(
                    "[main]",
                    """[main]
job_manager=local
MPIcommand = mpiexec -n @nprocs@ @command@ > @output@ 2> @error@""",
                )
            )

    (args, test_list, suite) = regtest.setup([test_ini.as_posix()])
    suite.workDir = tmp_path
    runtimes = {}
    runners = [regtest.TestRunner(test, args, test_list, suite, runtimes) for test in test_list]
    tasks, _srun_sh = regtest.srun_script(suite, runners, args)

    # assert EXAMPLE_SH.strip() == srun_sh.strip()
    assert tasks == 30
