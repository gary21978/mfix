"""This module is used to handle all of the git operations for the
test suite"""

import os
import shutil


from rtl2 import test_util


class Repo:
    """a simple class to manage our git operations"""

    def __init__(
        self,
        suite,
        directory,
        name,
        branch_wanted=None,
        pr_wanted=None,
        hash_wanted=None,
        build=0,
        comp_string=None,
    ):

        self.suite = suite
        self.dir = directory
        self.name = name
        self.branch_wanted = branch_wanted
        self.pr_wanted = pr_wanted
        self.hash_wanted = hash_wanted
        self.log = None

        self.build = build  # does this repo contain build directories?
        self.comp_string = comp_string  # environment vars needed to build

        # for storage
        self.branch_orig = None
        self.hash_current = None

        self.update = True
        if hash_wanted:
            self.update = False

    def get_branch_name(self):
        """for descriptive purposes, return the name of the branch we will
        use.  This could be a PR branch that was fetched"""
        if self.pr_wanted is not None:
            return "pr-{}".format(self.pr_wanted)
        if self.branch_wanted is not None:
            return self.branch_wanted.strip('"')

        return None

    def save_head(self):
        """Save the current head of the repo"""

        os.chdir(self.dir)

        self.suite.log.log("saving git HEAD for {}/".format(self.name))

        stdout, _, _ = test_util.run(
            "git rev-parse HEAD", self.log, outfile=f"git.{self.name}.HEAD"
        )

        self.hash_current = stdout
        shutil.copy("git.{}.HEAD".format(self.name), self.suite.full_web_dir)

    def make_changelog(self):
        """generate a ChangeLog git repository, and copy it to the
        web directory"""

        os.chdir(self.dir)

        self.suite.log.log("generating ChangeLog for {}/".format(self.name))

        test_util.run(
            "git log --name-only",
            self.log,
            outfile="ChangeLog.{}".format(self.name),
            outfile_mode="w",
        )
        shutil.copy("ChangeLog.{}".format(self.name), self.suite.full_web_dir)
