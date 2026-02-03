import configparser

import getpass
import os
import re
import socket
from tempfile import mkdtemp
from pathlib import Path

from rtl2 import repo
from rtl2 import suite
from rtl2 import test_util


def convert_type(string):
    """return an integer, float, or string from the input string"""
    if string is None:
        return None

    try:
        int(string)
    except:
        pass
    else:
        return int(string)

    try:
        float(string)
    except:
        pass
    else:
        return float(string)

    return string.strip()


def safe_get(cp, sec, opt, default=None):
    try:
        v = cp.get(sec, opt)
    except:
        v = default
    return v


def load_params(args):
    """
    reads the parameter file and creates as list of test objects as well as
    the suite object
    """

    test_list = []

    try:
        cp = configparser.ConfigParser(strict=False)
    except:
        cp = configparser.ConfigParser()

    cp.optionxform = str

    log = test_util.Log(output_file=args.log_file)

    log.bold("loading " + args.input_file[0])

    if not os.path.exists(args.input_file[0]):
        raise OSError("Parameter file {} does not exist".format(args.input_file[0]))

    try:
        cp.read(args.input_file[0])
    except:
        log.fail("ERROR: unable to read parameter file {}".format(args.input_file[0]))

    # "main" is a special section containing the global suite parameters.
    mysuite = suite.Suite(args)
    log.suite = mysuite
    mysuite.log = log

    valid_options = list(mysuite.__dict__.keys())
    for opt in cp.options("main"):

        # get the value of the current option
        value = convert_type(cp.get("main", opt))

        if opt in valid_options or "_" + opt in valid_options:

            if opt == "testTopDir":
                mysuite.testTopDir = mysuite.get_test_dir(value)
            elif opt == "workTopDir":
                mysuite.workTopDir = (
                    Path.cwd() / value if value.strip() else Path(mkdtemp(dir=mysuite.testTopDir))
                )
            elif opt == "webTopDir":
                mysuite.webTopDir = (
                    Path.cwd() / value if value.strip() else Path(mkdtemp(dir=mysuite.testTopDir))
                )
                mysuite._noWebDir = not value.strip()
            elif opt == "sourceDir":
                src = Path(value)
                mysuite.sourceDir = src if src.is_absolute() else Path.cwd() / src
            elif opt == "refdataDir":
                ref = Path(value)
                mysuite.refdataDir = ref if ref.is_absolute() else Path.cwd() / ref
            elif opt == "Label":
                mysuite.Label = value
            elif opt == "reportCoverage":
                mysuite.reportCoverage = mysuite.reportCoverage or value
            elif opt == "emailTo":
                mysuite.emailTo = value.split(",")

            else:
                # generic setting of the object attribute
                setattr(mysuite, opt, value)

        else:
            mysuite.log.warn("suite parameter {} not valid".format(opt))

    mysuite.repos["source"] = repo.Repo(
        mysuite,
        mysuite.sourceDir,
        mysuite.sourceDir.name,
        branch_wanted="develop",
        hash_wanted="HEAD",
        build=0,
        comp_string="",
    )

    # did we override the branch on the commandline?
    if args.source_branch is not None and args.source_pr is not None:
        mysuite.log.fail("ERROR: cannot specify both source_branch and source_pr")

    if args.source_branch is not None:
        mysuite.repos["source"].branch_wanted = args.source_branch

    if args.source_pr is not None:
        mysuite.repos["source"].pr_wanted = args.source_pr

    # now flesh out the compile strings -- they may refer to either themselves
    # or the source dir
    for r_repo in mysuite.repos.values():
        s = r_repo.comp_string
        if s is not None:
            r_repo.comp_string = s.replace("@self@", r_repo.dir.as_posix()).replace(
                "@source@", mysuite.repos["source"].dir.as_posix()
            )

    # the suite needs to know any ext_src_comp_string
    for r_repo in mysuite.repos.values():
        if not r_repo.build == 1:
            if r_repo.comp_string is not None:
                mysuite.extra_src_comp_string += " {} ".format(r_repo.comp_string)

    # checks
    if args.send_no_email:
        mysuite.sendEmailWhenFail = 0

    if args.post_only:
        mysuite.post_only = args.post_only

    suite.SUITE = mysuite

    if mysuite.sendEmailWhenFail:
        if mysuite.emailTo == [] or mysuite.emailBody == "":
            mysuite.log.fail(
                "ERROR: when sendEmailWhenFail = 1, you must specify emailTo and emailBody\n"
            )

        if mysuite.emailFrom == "":
            mysuite.emailFrom = "@".join((getpass.getuser(), socket.getfqdn()))

        if mysuite.emailSubject == "":
            mysuite.emailSubject = mysuite.suiteName + " Regression Test Failed"

    if mysuite.sourceDir == "" or mysuite.testTopDir == "":
        mysuite.log.fail(
            "ERROR: required suite-wide directory not specified\n" "(sourceDir, testTopDir)"
        )

    # all other sections are tests
    mysuite.log.skip()
    mysuite.log.bold("finding tests and checking parameters...")

    for sec in cp.sections():

        if sec in ["main", "AMReX", "source"] or sec.startswith("extra-"):
            continue

        # maximum test name length -- used for HTML formatting
        mysuite.lenTestName = max(mysuite.lenTestName, len(sec))

        # create the test object for this test
        mytest = suite.Test(sec)
        mytest.log = log
        invalid = 0

        # set the test object data by looking at all the options in
        # the current section of the parameter file
        valid_options = list(mytest.__dict__.keys())
        aux_pat = re.compile(r"aux\d+File")
        link_pat = re.compile(r"link\d+File")

        for opt in cp.options(sec):

            # get the value of the current option
            value = convert_type(cp.get(sec, opt))

            if opt in valid_options or "_" + opt in valid_options:

                if opt == "keyword":
                    mytest.keywords = [k.strip() for k in value.split(",")]

                else:
                    # generic setting of the object attribute
                    setattr(mytest, opt, value)

            elif aux_pat.match(opt):
                mytest.auxFiles.append(value)

            elif link_pat.match(opt):

                mytest.linkFiles.append(value)

            else:
                mysuite.log.warn(f"unrecognized parameter {opt} for test {sec}")

        mytest.auxFiles = [mytest.sdir(mysuite) / aux for aux in mytest.auxFiles]

        # make sure all the require parameters are present
        if mytest.compileTest:
            if mytest.buildDir == "":
                mysuite.log.warn("mandatory parameters for test {} not set".format(sec))
                invalid = 1

        else:

            input_file_invalid = mytest.inputFile == ""
            if input_file_invalid or mytest.dim == -1:
                warn_msg = [
                    "required params for test {} not set".format(sec),
                    "buildDir = {}".format(mytest.buildDir),
                    "inputFile = {}".format(mytest.inputFile),
                ]
                warn_msg += ["dim = {}".format(mytest.dim)]
                mysuite.log.warn(warn_msg)

                invalid = 1

        if mytest.inputFile is None:
            mysuite.log.warn("inputsFile not specified")
            invalid = 1
        else:
            mytest.inputFile = mytest.sdir(mysuite) / mytest.inputFile

        # check the optional parameters
        if mytest.restartTest and mytest.restartFileNum == -1:
            mysuite.log.warn("restart-test {} needs a restartFileNum".format(sec))
            invalid = 1

        if mytest.stSuccessString == "":
            mysuite.log.warn("self-test {} needs a stSuccessString".format(sec))
            invalid = 1

        if mytest.useMPI and mytest.numprocs == -1:
            mysuite.log.warn("MPI parallel test {} needs numprocs".format(sec))
            invalid = 1

        if mytest.useOMP and mytest.numthreads == -1:
            mysuite.log.warn("OpenMP parallel test {} needs numthreads".format(sec))
            invalid = 1

        # add the current test object to the master list
        if not invalid:
            test_list.append(mytest)
        else:
            mysuite.log.warn("test {} will be skipped".format(sec))

    # if any runs are parallel, make sure that the MPIcommand is defined
    any_MPI = any(t.useMPI for t in test_list)

    if not mysuite.job_manager:
        mysuite.log.fail(
            f"ERROR: In {args.input_file[0]},"
            " must set 'job_manager' to either 'local' or 'slurm'"
        )

    if mysuite.job_manager == "slurm" and not mysuite.partition:
        mysuite.log.fail(
            f"ERROR: In {args.input_file[0]},"
            " when 'job_manager' is 'slurm',"
            " then 'partition' must be set to the SLURM partition/queue."
        )

    if mysuite.job_manager == "slurm" and not mysuite.slurm_command:
        mysuite.log.fail(
            f"ERROR: In {args.input_file[0]},"
            " when 'job_manager' is 'slurm',"
            " then 'slurm_command' must be set to the 'srun' or 'salloc'."
        )

    if mysuite.slurm_command == "salloc" and not mysuite.ntasks_per_node:
        mysuite.log.fail(
            f"ERROR: In {args.input_file[0]},"
            " when 'slurm_command' is 'salloc',"
            " then 'ntasks_per_node' must be set."
        )

    if mysuite.slurm_command == "salloc" and not mysuite.ntasks_per_socket:
        mysuite.log.fail(
            f"ERROR: In {args.input_file[0]},"
            " when 'slurm_command' is 'salloc',"
            " then 'ntasks_per_socket' must be set."
        )

    if any_MPI and mysuite.MPIcommand == "":
        mysuite.log.fail("ERROR: some tests are MPI parallel, but MPIcommand not defined")

    test_list.sort()

    return mysuite, test_list
