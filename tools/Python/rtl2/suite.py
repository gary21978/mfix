from __future__ import print_function

import datetime
import json
import os
import glob
import shutil
import sys
import time
from pathlib import Path
from typing import List, Optional

from rtl2 import test_util

DO_TIMINGS_PLOTS = True

SUITE: Optional["Suite"] = None


def get_suite() -> "Suite":
    if SUITE is None:
        raise RuntimeError("global variable SUITE not initialized")
    return SUITE


try:
    import bokeh
    from bokeh.plotting import figure, save, ColumnDataSource
    from bokeh.resources import CDN
    from bokeh.models import HoverTool
    from datetime import datetime as dt

except:
    try:
        import matplotlib
    except:
        DO_TIMINGS_PLOTS = False
    else:
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

    try:
        from matplotlib import dates
    except:
        DO_TIMINGS_PLOTS = False


class Test:
    def __init__(self, name: str):

        self.name = name

        self.log = None

        self.buildDir = ""
        self.extra_build_dir = ""

        self.target = ""

        self.testSrcTree = ""

        self.inputFile = Path()
        self.probinFile = ""
        self.auxFiles: List[Path] = []
        self.linkFiles: List[str] = []

        self.dim = -1

        self.script_args = ""
        self.return_code = None

        self.restartTest = 0
        self.restartFileNum = -1

        self._compileTest = 0

        self.selfTest = 0
        self.stSuccessString = ""

        self.debug = 0

        self.acc = 0

        self.useMPI = 0
        self.numprocs = -1

        self.useOMP = 0
        self.numthreads = -1

        self._tolerance = None

        self.analysisRoutine = ""
        self.analysisMainArgs = ""
        self.analysisOutputImage = ""

        self.output_dir = Path()

        self.compare_file_used = ""

        self.diffDir = ""
        self.diffOpts = ""

        self.addToCompileString = ""
        self.ignoreGlobalMakeAdditions = 0

        self.runtime_params = ""

        self.reClean = 0  # set automatically, not by users

        self.wall_time = 0  # set automatically, not by users
        self.build_time = 0  # set automatically, not by users
        self.config_time = 0  # set automatically, not by users

        self.nlevels = None  # set but running fboxinfo on the output

        self.config_cmd = None  # set automatically

        self.comp_string = None  # set automatically
        self.run_command = None  # set automatically

        self.job_info_field1 = ""
        self.job_info_field2 = ""
        self.job_info_field3 = ""

        self.has_jobinfo = 0  # filled automatically

        self.backtrace: List[str] = []  # filled automatically

        self.has_stderr = False  # filled automatically

        self.post_successful = False  # filled automatically
        self.compare_successful = False  # filled automatically
        self.analysis_successful = False  # filled automatically

        self.customRunCmd: Optional[List[str]] = None

        self._check_performance = 0
        self._performance_threshold = 1.2
        self._runs_to_average = 5
        self.past_average = None

        self.keywords: List[str] = []

    def bdir(self, suite):
        return suite.testDir / self.buildDir

    def sdir(self, suite):
        return suite.sourceDir / self.buildDir

    def __lt__(self, other):
        return self.value() < other.value()

    def value(self):
        return self.name

    def find_backtrace(self):
        """find any backtrace files produced"""
        return [
            ft
            for ft in os.listdir(self.output_dir)
            if os.path.isfile(ft) and ft.startswith("Backtrace.")
        ]

    def get_compare_file(self, output_dir=None):
        """Find the last plotfile written.  Note: we give an error if the
        last plotfile is 0.  If output_dir is specified, then we use
        that instead of the default
        """

        if output_dir is None:
            output_dir = self.output_dir  # not yet implemented

        plts = [
            d
            for d in os.listdir(output_dir)
            if (os.path.isdir(d) and d.startswith("{}_plt".format(self.name)) and d[-1].isdigit())
            or (
                os.path.isfile(d)
                and d.startswith("{}_plt".format(self.name))
                and d.endswith(".tgz")
            )
        ]

        if len(plts) == 0:
            self.log.warn("test did not produce any output")
            return ""

        plts.sort()
        last_plot = plts.pop()

        if last_plot.endswith("00000"):
            self.log.warn("only plotfile 0 was output -- skipping comparison")
            return ""

        return last_plot

    def measure_performance(self):
        """returns performance relative to past average, as a tuple of:
        meets threshold, percentage slower/faster, whether slower/faster"""

        try:
            ratio = self.wall_time / self.past_average
        except (ZeroDivisionError, TypeError):
            return None, 0.0, "error computing ratio"

        meets_threshold = ratio < self.performance_threshold
        percentage = 100 * (1 - ratio)

        if percentage < 0:
            compare_str = "slower"
        else:
            compare_str = "faster"

        return meets_threshold, abs(percentage), compare_str

    def command(self, suite, base_command):
        cmd_s = " ".join(base_command)
        outfile = self.output_dir / self.outfile
        return (
            cmd_s
            if not self.useMPI
            else (
                suite.MPIcommand.replace("@host@", suite.MPIhost)
                .replace("@output@", str(outfile))
                .replace("@error@", str(outfile.with_suffix(".err")))
                .replace("@nprocs@", str(self.numprocs))
                .replace("@command@", cmd_s)
            )
        )

    def run_actual_test(self, suite, _runtimes, args):
        # ----------------------------------------------------------------------
        # run the test
        # ----------------------------------------------------------------------
        suite.log.log("running the test...")

        self.output_dir.mkdir(exist_ok=True)
        os.chdir(self.output_dir)

        self.wall_time = time.time()
        suite.run_test(self, self.base_command(suite, args), cwd=self.output_dir)
        self.wall_time = time.time() - self.wall_time

    def base_command(self, suite, args):

        if self.customRunCmd is not None:
            return self.customRunCmd

        base_cmd = [str(self.executable(suite)), str(self.inputFile)]

        if suite.plot_file_name != "":
            base_cmd.append(f"{suite.plot_file_name}={self.name}_plt")

        if suite.check_file_name != "none":
            base_cmd.append(f"{suite.check_file_name}={self.name}_chk")

        # keep around the checkpoint files only for the restart runs
        if self.restartTest:
            if suite.check_file_name != "none":
                base_cmd.extend(
                    ["mfix.checkpoint_files_output=1", f"mfix.check_int={self.restartFileNum}"]
                )
        else:
            base_cmd.append("mfix.checkpoint_files_output=0")

        base_cmd.extend(
            [str(suite.globalAddToExecString.strip()), str(self.runtime_params.strip())]
        )

        if args.with_valgrind:
            base_cmd = ["valgrind", str(args.valgrind_options), str(base_cmd)]

        return base_cmd

    #######################################################
    #           Static members and properties             #
    #######################################################

    @property
    def passed(self):
        """Whether the test passed or not"""

        if self.compileTest or not self.post_successful:
            return self.post_successful

        compare = self.compare_successful
        analysis = self.analysisRoutine == "" or self.analysis_successful
        return compare and analysis

    @property
    def crashed(self) -> bool:
        """Whether the test crashed or not"""
        return len(self.backtrace) > 0 or self.return_code

    @property
    def outfile(self):
        """The basename of this run's output file"""

        return "{}.run.out".format(self.name)

    @property
    def errfile(self) -> str:
        """The basename of this run's error file"""

        return "{}.err.out".format(self.name)

    @property
    def comparison_outfile(self):
        """The basename of this run's comparison output file"""

        return "{}.compare.out".format(self.name)

    def executable(self, suite):
        return (self.bdir(suite) / self.name).with_suffix(".ex")

    def record_runtime(self, suite):

        test = self.passed and not self.compileTest
        suite = not suite.args.do_temp_run
        return test and suite

    def set_compile_test(self, value):
        """Sets whether this test is compile-only"""

        self._compileTest = value

    def get_compile_test(self):
        """Returns True if the global --compile_only flag was set or
        this test is compile-only, False otherwise
        """

        return self._compileTest or Test.compile_only

    def get_tolerance(self):
        """Returns the global tolerance if one was set,
        and the test-specific one otherwise.
        """

        if Test.global_tolerance is None:
            return self._tolerance
        return Test.global_tolerance

    def set_tolerance(self, value):
        """Sets the test-specific tolerance to the specified value."""

        self._tolerance = value

    def get_check_performance(self):
        """Returns whether to check performance for this test."""

        return self._check_performance or Test.performance_params

    def set_check_performance(self, value):
        """Setter for check_performance."""

        self._check_performance = value

    def get_performance_threshold(self):
        """Returns the threshold at which to warn of a performance drop."""

        if Test.performance_params:
            return float(Test.performance_params[0])
        if self._check_performance:
            return self._performance_threshold
        return None

    def set_performance_threshold(self, value):
        """Setter for performance_threshold."""

        self._performance_threshold = value

    def get_runs_to_average(self):
        """Returns the number of past runs to include in the running runtime average."""

        if Test.performance_params:
            return int(Test.performance_params[1])
        if self._check_performance:
            return self._runs_to_average
        return None

    def set_runs_to_average(self, value):
        """Setter for runs_to_average."""

        self._runs_to_average = value

    # Static member variables, set explicitly in apply_args in Suite class
    compile_only = False
    global_tolerance = None
    performance_params: List[str] = []

    # Properties - allow for direct access as an attribute
    # (e.g. test.compileTest) while still utilizing getters and setters
    compileTest = property(get_compile_test, set_compile_test)
    tolerance = property(get_tolerance, set_tolerance)
    check_performance = property(get_check_performance, set_check_performance)
    performance_threshold = property(get_performance_threshold, set_performance_threshold)
    runs_to_average = property(get_runs_to_average, set_runs_to_average)


class Suite:
    def __init__(self, args):

        self.args = args
        self.apply_args()

        # this will hold all of the Repo() objects for the AMReX, source,
        # and build directories
        self.repos = {}

        self.test_file_path = os.getcwd() + "/" + self.args.input_file[0]

        self.suiteName = "testDefault"
        self.sub_title = ""

        self.testTopDir: Optional[Path] = None
        self.webTopDir = ""
        self._noWebDir = False
        self.wallclockFile = "wallclock_history"

        self.launch_dir = os.getcwd()

        self.reportCoverage = args.with_coverage

        # set automatically
        self.covered_frac = None
        self.total = None
        self.covered_nonspecific_frac = None
        self.total_nonspecific = None

        self.sourceDir = ""
        self.source_build_dir = ""  # Cmake build dir
        self.workDir: Optional[Path] = None
        self.workTopDir: Optional[Path] = None
        self._full_test_dir: Optional[Path] = None
        self.cmakeSetupOpts = ""

        self.refdataDir: Path = Path()
        self.Label = None

        self.MPIcommand = ""
        self.MPIhost = ""

        self.job_manager = ""
        self.partition = ""
        self.cpus_per_task = 1
        self.slurm_command = ""
        self.ntasks_per_node = None
        self.ntasks_per_socket = None

        self.reportActiveTestsOnly = 0
        self.goUpLink = 0
        self.lenTestName = 0

        self.sendEmailWhenFail = 0
        self.emailFrom = ""
        self.emailTo = []
        self.emailSubject = ""
        self.emailBody = ""

        self.post_only = None

        self.plot_file_name = "mfix.plot_file"
        self.check_file_name = "mfix.check_file"

        self.globalAddToExecString = ""

        # this will be automatically filled
        self.extra_src_comp_string = ""

        # delete all plot/checkfiles but the plotfile used for comparison upon
        # completion
        self.purge_output = 0

        self.log = None

        self.do_timings_plots = DO_TIMINGS_PLOTS

        # default branch -- we use this only for display purposes --
        # if the test was run on a branch other than the default, then
        # an asterisk will appear next to the date in the main page
        self.default_branch = "develop"

    @property
    def buildDir(self) -> Path:
        if self.workDir is None:
            raise RuntimeError("invalid workDir")
        return self.workDir / "build"

    @property
    def testDir(self) -> Path:
        if self.workDir is None:
            raise RuntimeError("invalid workDir")
        return self.workDir / "run"

    @property
    def timing_default(self):
        """Determines the format of the wallclock history JSON file"""

        return {"runtimes": [], "dates": []}

    def refdata(self, test: "Test") -> Path:
        prefix = f"{self.Label}." if self.Label else ""
        return self.refdataDir / test.name / f"{prefix}runningave.dat"

    def get_test_dir(self, dir_name: str) -> Path:
        """given a string representing a directory, check if it points to
        a valid directory.  If so, return the directory name"""

        path = Path(dir_name)
        if path.is_absolute():
            return path

        if self.testTopDir is None:
            return Path.cwd() / dir_name
        return Path(self.testTopDir) / dir_name

    def check_test_dir(self, path):
        """given a string representing a directory, check if it points to
        a valid directory.  If so, return the directory name"""

        if path.is_dir():
            return path

        # we failed :(
        self.log.fail("ERROR: {} is not a valid directory".format(path))
        return None

    def init_build_dir(self, dir_name):
        """
        Sets the suite build directory to dir_name if dir_name is neither null
        nor whitespace, and initializes it to a temporary directory under the
        test directory otherwise.
        """

    def delete_tempdirs(self):
        """
        Removes any temporary directories that were created during the
        current test run.
        """

        if self._noWebDir:
            shutil.rmtree(self.webTopDir)

    def get_tests_to_run(self, test_list_old):
        """perform various tests based on the runtime options to determine
        which of the tests in the input file we run"""

        # if we only want to run the tests that failed previously,
        # remove the others
        if self.args.redo_failed:
            last_run = self.get_last_run()
            failed = self.get_test_failures(last_run)

            test_list = [t for t in test_list_old if t.name in failed]
        else:
            test_list = test_list_old[:]

        # if we only want to run tests of a certain dimensionality, remove
        # the others
        if self.args.d in [1, 2, 3]:
            test_list = [t for t in test_list_old if t.dim == self.args.d]

        # if we specified any keywords, only run those
        if self.args.keyword is not None:
            test_list = [t for t in test_list_old if self.args.keyword in t.keywords]

        # if we are doing a single test, remove all other tests; if we
        # specified a list of tests, check each one; if we did both
        # --single_test and --tests, complain
        if not self.args.single_test == "" and not self.args.tests == "":
            self.log.fail("ERROR: specify tests either by --single_test or --tests, not both")

        if not self.args.single_test == "":
            tests_find = [self.args.single_test]
        elif not self.args.tests == "":
            tests_find = self.args.tests.split()
        else:
            tests_find = []

        if len(tests_find) > 0:
            new_test_list = []
            for test in tests_find:
                _tmp = [o for o in test_list if o.name == test]
                if len(_tmp) == 1:
                    new_test_list += _tmp
                else:
                    self.log.fail("ERROR: {} is not a valid test".format(test))

            test_list = new_test_list

        if len(test_list) == 0:
            self.log.fail("No valid tests defined")

        return test_list

    def get_bench_dir(self):
        bench_dir = self.testTopDir / f"{self.suiteName}-benchmarks"
        if not bench_dir.is_dir():
            bench_dir.mkdir()
        return bench_dir

    def get_wallclock_file(self):
        """returns the path to the json file storing past runtimes for each test"""

        return self.get_bench_dir() / f"{self.wallclockFile}.json"

    @property
    def full_test_dir(self) -> Path:
        if self.workDir is None:
            raise RuntimeError("invalid workDir")
        return self.workDir

    def _getFresh(self, base: Path) -> Path:
        today = str(datetime.date.today())

        # figure out what the current output directory should be
        maxRuns = 100  # maximum number of tests in a given day

        if self.args.do_temp_run:
            return base / "TEMP_RUN"

        for i in range(1, maxRuns):
            full_dir = base / "{}-{:03d}".format(today, i)
            if not full_dir.is_dir():
                return full_dir

        raise RuntimeError("NO TEST DIR")

    def make_test_dirs(self):
        os.chdir(self.testTopDir)

        if self.args.do_temp_run and self.full_test_dir.is_dir():
            shutil.rmtree(self.full_test_dir)

        workTop = self.workTopDir if self.workTopDir is not None else self.testTopDir / "work"
        self.workDir = self._getFresh(workTop)
        if self.post_only is not None:
            self.workDir = workTop / self.post_only
            if not self.workDir.is_dir():
                self.log.fail(f"Missing --post_only directory {self.workDir}")

        webTop = self.webTopDir if self.webTopDir is not None else self.testTopDir / "web"
        self.webDir = webTop / self.workDir.name

        self.log.skip()
        self.log.bold(f"testing directory is: {self.workDir}")
        if not self.buildDir.is_dir():
            self.buildDir.mkdir(parents=True)

        # make the web directory -- this is where all the output and HTML will be
        # put, so it is easy to move the entire test website to a different disk
        self.full_web_dir = Path(self.webTopDir) / self.workDir.name

        if self.args.do_temp_run and self.full_web_dir.is_dir():
            shutil.rmtree(self.full_web_dir)

        self.full_web_dir.mkdir(parents=True, exist_ok=True)

        # copy the test file into the web output directory
        shutil.copy(self.test_file_path, self.full_web_dir)

    def get_run_history(self, active_test_list=None, check_activity=True):
        """return the list of output directories run over the
        history of the suite and a separate list of the tests
        run (unique names)"""

        valid_dirs = []
        all_tests = []

        # start by finding the list of valid test directories
        for f in self.webTopDir.iterdir():
            # look for a directory of the form 20* (this will work up until 2099
            if f.name.startswith("20") and f.is_dir():
                if f.glob("*.status"):
                    valid_dirs.append(f)

        valid_dirs.sort()
        valid_dirs.reverse()

        # now find all of the unique problems in the test directories
        for valid_dir in valid_dirs:

            for f in valid_dir.glob("*.status"):
                test_name = f.stem
                if test_name == "branch":
                    continue
                if test_name not in all_tests:
                    is_active = active_test_list and (test_name in active_test_list)
                    if (not (self.reportActiveTestsOnly and check_activity)) or is_active:
                        all_tests.append(test_name)

        all_tests.sort()

        return valid_dirs, all_tests

    def get_wallclock_history(self):
        """returns the wallclock time history for all the valid tests as a dictionary
        of NumPy arrays. Set filter_times to False to return 0.0 as a placeholder
        when there was no available execution time."""

        def extract_time(file):
            """Helper function for getting runtimes"""

            for line in file:

                if "Execution time" in line:
                    # this is of the form: <li>Execution time: 412.930 s
                    return float(line.split(":")[1].strip().split(" ")[0])

                if "(seconds)" in line:
                    # this is the older form -- split on "="
                    # form: <p><b>Execution Time</b> (seconds) = 399.414828
                    return float(line.split("=")[1])

            raise RuntimeError()

        json_file = self.get_wallclock_file()
        if json_file.is_file():
            with open(json_file, "r") as f:
                try:
                    timings = json.load(f)
                    if not isinstance(timings, dict):
                        self.log.warn(f"Incorrect JSON: {json_file}")
                    return timings
                except json.JSONDecodeError:
                    self.log.warn(f"Unable to parse {json_file}")

        valid_dirs, all_tests = self.get_run_history(check_activity=False)

        # store the timings in a dictionary
        timings = {}

        for dir_ in valid_dirs:

            # Get status files
            dir_path = os.path.join(self.webTopDir, dir_)
            sfiles = glob.glob("{}/*.status".format(dir_path))
            sfiles = list(filter(os.path.isfile, sfiles))

            # Tests that should be counted
            passed = set()

            for i, file in enumerate(map(open, sfiles)):

                contents = file.read()

                if "PASSED" in contents:
                    filename = os.path.basename(sfiles[i])
                    passed.add(filename.split(".")[0])

                file.close()

            for test in all_tests:
                if test not in passed:
                    continue

                fname = "{}/{}.html".format(dir_path, test)
                with open(fname) as file_:
                    extract_time(file_)
                    test_dict = timings.setdefault(test, self.timing_default)
                    test_dict["runtimes"].append(time)
                    test_dict["dates"].append(dir_.name)

        return timings

    def make_timing_plots(self, active_test_list=None, all_tests=None):
        """plot the wallclock time history for all the valid tests"""

        if active_test_list is not None:
            _, all_tests = self.get_run_history(active_test_list)
        timings = self.get_wallclock_history()

        try:
            bokeh
        except NameError:

            convf = dates.datestr2num
            using_mpl = True
            self.plot_ext = "png"

        else:

            convf = lambda s: dt.strptime(s, "%Y-%m-%d")
            using_mpl = False
            self.plot_ext = "html"

        def convert_date(date):
            """Convert to a matplotlib readable date"""

            if len(date) > 10:
                date = date[: date.rfind("-")]
            return convf(date)

        def hover_tool():
            """
            Encapsulates hover tool creation to prevent errors when generating
            multiple documents.
            """

            return HoverTool(
                tooltips=[("date", "@date{%F}"), ("runtime", "@runtime{0.00}")],
                formatters={"date": "datetime"},
            )

        # make the plots
        for t in all_tests:

            try:
                test_dict = timings[t]
            except KeyError:
                continue

            days = list(map(convert_date, test_dict["dates"]))
            times = test_dict["runtimes"]

            if len(times) == 0:
                continue

            if using_mpl:

                plt.clf()
                plt.plot_date(days, times, "o", xdate=True)

                years = dates.YearLocator()  # every year
                months = dates.MonthLocator()
                years_fmt = dates.DateFormatter("%Y")

                ax = plt.gca()
                ax.xaxis.set_major_locator(years)
                ax.xaxis.set_major_formatter(years_fmt)
                ax.xaxis.set_minor_locator(months)

                plt.ylabel("time (seconds)")
                plt.title(t)

                if max(times) / min(times) > 10.0:
                    ax.set_yscale("log")

                fig = plt.gcf()
                fig.autofmt_xdate()

                plt.savefig("{}/{}-timings.{}".format(self.webTopDir, t, self.plot_ext))

            else:

                source = ColumnDataSource(dict(date=days, runtime=times))

                settings = dict(x_axis_type="datetime")
                if max(times) / min(times) > 10.0:
                    settings["y_axis_type"] = "log"
                plot = figure(**settings)
                plot.add_tools(hover_tool())

                plot.circle(  # pylint: disable=too-many-function-args
                    "date",
                    "runtime",
                    source=source,
                )
                plot.xaxis.axis_label = "Date"
                plot.yaxis.axis_label = "Runtime (s)"

                save(
                    plot,
                    resources=CDN,
                    filename="{}/{}-timings.{}".format(self.webTopDir, t, self.plot_ext),
                    title=f"{t} Runtime History",
                )

    def get_last_run(self):
        """return the name of the directory corresponding to the previous
        run of the test suite"""

        outdir = self.testDir

        # this will work through 2099
        if os.path.isdir(outdir):
            dirs = [
                d for d in os.listdir(outdir) if (os.path.isdir(outdir + d) and d.startswith("20"))
            ]
            dirs.sort()

            return dirs[-1]
        return None

    def get_test_failures(self, test_dir):
        """look at the test run in test_dir and return the list of tests that
        failed"""

        cwd = os.getcwd()

        outdir = self.testDir

        os.chdir(outdir + test_dir)

        failed = []

        for test in os.listdir("."):
            if not os.path.isdir(test):
                continue

            # the status files are in the web dir
            status_file = "{}/{}/{}.status".format(self.webTopDir, test_dir, test)
            with open(status_file, "r") as sf:
                for line in sf:
                    if line.find("FAILED") >= 0 or line.find("CRASHED") >= 0:
                        failed.append(test)

        os.chdir(cwd)
        return failed

    def command(self, test, base_command):
        return (
            (
                self.MPIcommand.replace("@host@", self.MPIhost)
                .replace("@output@", str(test.output_dir / test.outfile))
                .replace("@nprocs@", str(test.numprocs))
                .replace("@command@", base_command)
            )
            if test.useMPI
            else base_command
        )

    def run_test(self, test, base_command, cwd):
        test_run_command = test.command(self, base_command)

        self.log.log(test_run_command)
        _, _, ierr = test_util.run(
            test_run_command,
            test.log,
            stdin=True,
            outfile=cwd / test.outfile,
            errfile=test.errfile,
            cwd=cwd,
        )
        test.run_command = test_run_command
        test.return_code = ierr

    def copy_backtrace(self, test):
        """
        if any backtrace files were output (because the run crashed), find them
        and copy them to the web directory
        """
        backtrace = test.find_backtrace()

        for btf in backtrace:
            ofile = f"{self.full_web_dir}/{test.name}.{btf}"
            shutil.copy(btf, ofile)
            test.backtrace.append("{}.{}".format(test.name, btf))

    def apply_args(self):
        """
        makes any necessary adjustments to module settings based on the
        command line arguments supplied to the main module
        """

        args = self.args

        Test.compile_only = args.compile_only
        Test.global_tolerance = args.tolerance
        Test.performance_params = args.check_performance

    #######################################################
    #        CMake utilities                              #
    #######################################################
    def cmake_config(self):
        "Generate Cmake configuration"
        name = self.suiteName
        path = self.sourceDir
        configOpts = self.cmakeSetupOpts
        if name == "AMReX":
            return (None, None)

        self.log.outdent()
        self.log.skip()
        self.log.bold("configuring " + name + " build...")
        self.log.indent()

        # Setup dir names
        builddir = self.buildDir

        # Logfile
        coutfile = self.workDir / f"{name}.cmake.log"

        # Run cmake
        cmd = f"cmake {configOpts} -S{path} -B{builddir} "
        assert name != "AMReX"

        self.log.log(cmd)
        assert isinstance(cmd, str)
        _, _, rc = test_util.run(cmd, self.log, outfile=coutfile)
        self.config_cmd = cmd

        # Check exit condition
        if rc:
            errstr = (
                f"\n \nERROR! Cmake configuration failed for {name} \n"
                f"Check {coutfile} for more information:"
            )
            with open(coutfile) as output:
                errstr += output.read()
            self.log.fail(errstr)
            sys.exit(errstr)

        self.source_build_dir = builddir
        shutil.copy(coutfile, self.full_web_dir / f"{name}.cmake.log")

    def cmake_clean(self, name, path):
        "Clean Cmake build and install directories"

        self.log.outdent()
        self.log.skip()
        self.log.bold("cleaning " + name + " Cmake directories...")
        self.log.indent()

        # Setup dir names
        builddir = path / "builddir"
        installdir = path / "installdir"

        # remove build and installation directories if present
        if os.path.isdir(builddir):
            shutil.rmtree(builddir)

        if os.path.isdir(installdir):
            shutil.rmtree(installdir)

    def cmake_build(self, path, outfile):
        "Build target for a repo configured via cmake"

        self.log.outdent()
        self.log.skip()
        self.log.bold(f"cmake_build:")
        self.log.indent()

        cmd = f"cmake --build {path} -j4"
        self.log.log(cmd)
        _, _, rc = test_util.run(cmd, self.log, outfile=outfile)

        # make returns 0 if everything was good
        if rc:
            errstr = f"Failed to build targets: {target}\n Check {outfile} for more information."
            self.log.fail(errstr)

        comp_string = cmd

        return rc, comp_string

    def build_test_cmake(self, test, outfile):
        """build an executable with CMake build system"""

        rc, comp_string = self.cmake_build(
            path=self.source_build_dir,
            outfile=outfile,
        )

        # make returns 0 if everything was good
        if rc:
            self.log.fail(f"Failed to build test {test.name}")
            return comp_string, rc

        for path_to_exe in self.source_build_dir.glob(f"**/{test.target}"):
            # Symlink executable to test dir
            if test.executable(self).exists():
                test.executable(self).unlink()
            test.executable(self).parent.mkdir(exist_ok=True, parents=True)
            os.link(path_to_exe, test.executable(self))
            break
        else:
            self.log.fail(f"Unable to find {test.target} under {self.source_build_dir}")

        return comp_string, rc


def f_flag(opt, test_not=False):
    """convert a test parameter into t if true for the Fortran build system"""
    if test_not:
        return " " if opt else "t"
    return "t" if opt else " "


def c_flag(opt, test_not=False):
    """convert a test parameter into t if true for the Fortran build system"""
    if test_not:
        return "FALSE" if opt else "TRUE"
    return "TRUE" if opt else "FALSE"
