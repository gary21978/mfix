from __future__ import print_function
from pathlib import Path

import argparse
import datetime
import os
import shlex
import subprocess
import sys

from typing import Tuple

USAGE = """
The test suite and its tests are defined through an input file in an INI
configuration file format.

The "main" block specifies the global test suite parameters:

  [main]

  testTopDir     = < path to test output directory >
  webTopDir      = < path to test web output directory >
  workDir        = < path to work directory (for building and running) >
  sourceDir      = < path to source directory >
  wallclockFile  = < name of json file for storing past runtimes, to which .json will be appended;
                     set to wallclock_history by default >

  suiteName = < descriptive name (i.e. MFIX-Exa) >

  reportActiveTestsOnly = <0: web shows every test ever run;
                           1: just current tests >

  goUpLink = <1: add "Go UP" link at top of the web page >

  purge_output = <0: leave all plotfiles in place;
                  1: delete plotfiles after compare >

  job_manager = <local: run on current host;
                 slurm: submit to SLURM with srun Multiple Program configuration file >

  partition = < queue/partition to submit SLURM job; ignored if job_manager is not "slurm" >

  slurm_command = <srun: submit to SLURM using srun (preferred for CPU tests);
                   salloc: submit to SLURM using salloc (required for GPU tests) >

  ntasks_per_node = < when the slurm_command is "salloc" this value must be set >

  ntasks_per_socket = < when the slurm_command is "salloc" this value must be set >

  MPIcommand = < MPI run command, with holders for host, # of proc, command >

     This should look something like:

          mpiexec -host @host@ -n @nprocs@ @command@ >

  MPIhost = < host for MPI job -- depends on MPI implementation >

  sendEmailWhenFail = < 1: send email when any tests fail >

  emailTo = < list of email addresses separated by commas, such as,
              foo@example.com, bar@example.com >

  emailBody = < email body >


Each test is given its own block, with the general form:

  [fluidbed_cyl]

  buildDir = < relative path (from sourceDir) for this problem >

  target = <name of the make target associated to the test (cmake only)>

  inputFile = < input file name >

  dim = < dimensionality: 1, 2, or 3 >

  aux?File = < name of additional file needed by the test >
  link?File = < name of additional file needed by the test >

      Here "?" is 1, 2, or 3, allowing for several files per test

  restartTest = < is this a restart test? 0 for no, 1 for yes >
  restartFileNum = < # of file to restart from (if restart test) >

  useMPI = <is this a parallel (MPI) job? 0 for no, 1 for yes) >
  numprocs = < # of processors to run on (if parallel job) >

  useOMP = <is this an OpenMP job? 0 for no, 1 for yes) >
  numthreads = < # of threads to us with OpenMP (if OpenMP job) >

  debug = < 0 for normal run, 1 if we want debugging options on >

  compileTest = < 0 for normal run, 1 if we just test compilation >

  stSuccessString = < string to find in output to determine success >

  analysisRoutine = < name of the script to run on the output >

      The script is run as:

        analysisRoutine [options] plotfile

  analysisMainArgs = < commandline arguments to pass to the analysisRoutine --
                       these should refer to options from the [main] block >

  analysisOutputImage = < name on analysis result image to show on web page >

  tolerance = < floating point number representing the largest relative error
                permitted >

  diffDir = < directory/file to do a plain text diff on (recursive, if dir) >

  diffOpts = < options to use with the diff command for the diffDir comparison >

  check_performance = < 1: compare run time of test to average of past runs >
  performance_threshold = < ratio of run time / running average above which a
                            a performance warning will be issued, default is
                            1.2 >
  runs_to_average = < number of past runs to include when computing the average,
                      default is 5 >

Getting started:

To set up a test suite, it is probably easiest to write the
testfile.ini as described above and then run the test routine.
"""


class Log:
    """a simple logging class to show information to the terminal"""

    def __init__(self, output_file=None):

        # http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
        # which in-turn cites the blender build scripts
        self.warn_color = "\033[33m"
        self.success_color = "\033[32m"
        self.fail_color = "\033[31m"
        self.bold_color = "\033[1m"
        self.end_color = "\033[0m"

        self.current_indent = 0
        self.indent_str = ""

        if output_file is not None:
            try:
                self.of = open(output_file, "w")
            except IOError as err:
                print("ERROR: unable to open output file")
                raise IOError from err
            else:
                self.have_log = True
        else:
            self.of = None
            self.have_log = False

        self.suite = None

    def indent(self):
        """indent the log output by one stop"""
        self.current_indent += 1
        self.indent_str = self.current_indent * "   "

    def flush(self):
        """flush the output file (if it exists)"""
        if self.have_log:
            self.of.flush()

    def outdent(self):
        """undo one level of indent"""
        self.current_indent -= 1
        self.current_indent = max(0, self.current_indent)
        self.indent_str = self.current_indent * "   "

    def fail(self, string):
        """output a failure message to the log"""
        self._write(self.fail_color + string + self.end_color)
        self.close_log()
        sys.exit()

    def testfail(self, string):
        """output a test failure to the log"""
        self._write(self.fail_color + string + self.end_color)

    def warn(self, warn_msg):
        """
        output a warning.  It is always prefix with 'WARNING:'
        For multi-line warnings, send in a list of strings
        """
        prefix = self.indent_str + "WARNING: "
        filler = self.indent_str + "         "
        if isinstance(warn_msg, list):
            msg = [prefix + warn_msg[0]] + [filler + x for x in warn_msg[1:]]
            omsg = "\n".join(msg).strip()
        else:
            omsg = prefix + warn_msg
        self._write(self.warn_color + omsg + self.end_color)

    def success(self, string):
        """output a success message to the log"""
        self._write(self.success_color + string + self.end_color)

    def log(self, string):
        """output some text to the log"""
        self._write(string)

    def skip(self):
        """introduce a newline in the log"""
        self._write("")

    def bold(self, string):
        """emphasize a log message"""
        self._write(self.bold_color + string + self.end_color)

    def _write(self, string):
        """emphasize a log message"""
        now = datetime.datetime.now().isoformat()
        str_ = f"{now}  {self.indent_str}{string}"
        print(str_)
        if self.have_log:
            self.of.write(f"{str_}\n")

    def close_log(self):
        """close the log"""
        if self.have_log:
            self.of.close()


def get_args(arg_string=None):
    """parse the commandline arguments.  If arg_string is present, we
    parse from there, otherwise we use the default (sys.argv)"""

    parser = argparse.ArgumentParser(
        description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    tests_group = parser.add_argument_group(
        "test selection", "options that determine which tests to run"
    )
    tests_group.add_argument(
        "--single_test",
        type=str,
        default="",
        metavar="test-name",
        help="name of a single test to run",
    )
    tests_group.add_argument(
        "--tests",
        type=str,
        default="",
        metavar="'test1 test2 test3'",
        help="a space-separated list of tests to run",
    )
    tests_group.add_argument(
        "-d", type=int, default=-1, help="restrict tests to a particular dimensionality"
    )
    tests_group.add_argument(
        "--redo_failed", action="store_true", help="only run the tests that failed last time"
    )
    tests_group.add_argument(
        "--keyword",
        type=str,
        default=None,
        help="run tests only with this keyword specified in their definitions",
    )

    git_group = parser.add_argument_group(
        "git options", "options that control how we interact the git repos"
    )
    git_group.add_argument(
        "--no_update",
        type=str,
        default="None",
        metavar="name",
        help="which codes to exclude from the git update?"
        + " (None, All, or a comma-separated list of codes)",
    )
    git_group.add_argument(
        "--source_branch",
        type=str,
        default=None,
        metavar="branch-name",
        help="what git branch to use for the source repo",
    )
    git_group.add_argument(
        "--source_pr",
        type=int,
        default=None,
        metavar="PR-number",
        help="what github pull request number to use for the source repo",
    )
    git_group.add_argument(
        "--amrex_pr",
        type=int,
        default=None,
        metavar="PR-number",
        help="what github pull request number to use for the amrex repo",
    )

    run_group = parser.add_argument_group(
        "test running options", "options that control how the tests are run"
    )
    run_group.add_argument(
        "--compile_only",
        action="store_true",
        help="test only that the code compiles, without running anything",
    )
    run_group.add_argument(
        "--post_only",
        type=str,
        default=None,
        metavar="RUN_NAME",
        help="post-process (regenerate web reports) for the specified run,"
        " without building or running anything. For use after editing post.py scripts."
        " Build time/cmd/output will be empty in regenerated web report.",
    )

    run_group.add_argument("--with_valgrind", action="store_true", help="run with valgrind")
    run_group.add_argument(
        "--valgrind_options",
        type=str,
        default="--leak-check=yes --log-file=vallog.%p",
        help="valgrind options",
        metavar="'valgrind options'",
    )

    suite_options = parser.add_argument_group(
        "suite options", "options that control the test suite operation"
    )
    suite_options.add_argument(
        "--do_temp_run",
        action="store_true",
        help="is this a temporary run? (output not stored or logged)",
    )
    suite_options.add_argument(
        "--send_no_email", action="store_true", help="do not send emails when tests fail"
    )
    suite_options.add_argument(
        "--with_coverage", action="store_true", help="report parameter coverage for this test run"
    )
    suite_options.add_argument(
        "--check_performance",
        nargs=2,
        metavar=("performance_threshold", "runs_to_average"),
        help="measure the performance of each test run against the last runs_to_average runs, "
        + "supplying a warning on a ratio greater than performance_threshold",
    )
    suite_options.add_argument(
        "--note", type=str, default="", help="a note on the resulting test webpages"
    )
    suite_options.add_argument(
        "--complete_report_from_crash",
        type=str,
        default="",
        metavar="testdir",
        help="complete report generation from a crashed test suite run named testdir",
    )
    suite_options.add_argument(
        "--log_file",
        type=str,
        default=None,
        metavar="logfile",
        help="log file to write output to (in addition to stdout",
    )

    comp_options = parser.add_argument_group(
        "comparison options", "options that control how the comparisons are done"
    )
    comp_options.add_argument(
        "--tolerance",
        type=float,
        default=None,
        metavar="value",
        help="largest relative error permitted",
    )

    parser.add_argument(
        "input_file",
        metavar="input-file",
        type=str,
        nargs=1,
        help="the input file (INI format) containing the suite and test parameters",
    )

    if arg_string is not None:
        args = parser.parse_args(arg_string)
    else:
        args = parser.parse_args()

    return args


def run(
    string,
    log,
    stdin=False,
    outfile=None,
    store_command=False,
    outfile_mode="a",
    errfile=None,
    cwd=None,
):
    # shlex.split will preserve inner quotes
    prog = shlex.split(string)
    sin = None
    if stdin:
        sin = subprocess.PIPE

    p0 = subprocess.Popen(
        prog,
        stdin=sin,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=None,
        cwd=cwd,
    )

    stdout0, stderr0 = p0.communicate()

    stdout0 = stdout0.decode("utf-8")
    stderr0 = stderr0.decode("utf-8")

    if outfile is not None:
        outfile.parent.mkdir(parents=True, exist_ok=True)
        if outfile_mode == "a":
            outfile.touch(exist_ok=True)
        try:
            with open(outfile, outfile_mode) as cf:
                if store_command:
                    cf.write(p0.args)
                    cf.write("\n")
                for line in stdout0:
                    cf.write(line)

                if errfile is None:
                    for line in stderr0:
                        cf.write(line)
        except IOError:
            log.fail(f"  ERROR: unable to open outfile for writing {outfile}")

    if errfile is not None and stderr0.strip():
        try:
            with open(errfile, outfile_mode) as cf:
                for line in stderr0:
                    cf.write(line)
        except IOError:
            log.fail(f"  ERROR: unable to open errfile for writing {errfile}")

    return stdout0, stderr0, p0.returncode


def get_recent_filename(fdir, base, extension):
    """find the most recent file matching the base and extension"""

    files = sorted(
        [f for f in Path(fdir).iterdir() if (f.name.startswith(base) and f.suffix == extension)],
        key=os.path.getmtime,
    )
    if files:
        return files.pop()
    return None


def git_commit(path: Path) -> Tuple[str, str]:
    """Return tuple of (branch, SHA) for the Git repo at path
    if unable to read Git repo, returns ("UNKNOWN", "UNKNOWN")
    """
    branch = "Unknown"
    sha = "Unknown"
    try:
        branch = subprocess.check_output(
            ["git", "-C", path.as_posix(), "rev-parse", "--abbrev-ref", "HEAD"]
        ).decode("utf-8")
    except subprocess.CalledProcessError:
        print("Could not find current branch for Git repo: ", path.as_posix())
    try:
        sha = subprocess.check_output(
            ["git", "-C", path.as_posix(), "rev-parse", "--verify", "HEAD"]
        ).decode("utf-8")
    except subprocess.CalledProcessError:
        print("Could not find current commit hash for Git repo: ", path.as_posix())
    return (branch.strip(), sha.strip())
