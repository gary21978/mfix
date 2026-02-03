# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

CSG_REPO = "https://mfix.netl.doe.gov/gitlab/exa/csg-eb.git"
CSG_BRANCH = "main"


class CsgEb(CMakePackage):
    """Define Embedded Boundaries (for AMReX) using CSG"""

    homepage = "https://mfix.netl.doe.gov/gitlab/exa/csg-eb"

    version("main", git=CSG_REPO, branch=CSG_BRANCH)

    variant("cgal", default=True, description="Build with CGAL Support")

    depends_on("pegtl")
    depends_on("catch2")
    depends_on("cgal header_only=True", when="+cgal")

    def cmake_args(self):
        return [
            "-DCSG_CGAL_ENABLED={}".format("True" if "+cgal" in self.spec else "False")
        ]

    def setup_run_environment(self, env):
        env.prepend_path("CPATH", self.prefix.include)

    @property
    def headers(self):
        headers = find_all_headers(self.prefix.include)
        headers.directories = [self.prefix.include]
        return headers
