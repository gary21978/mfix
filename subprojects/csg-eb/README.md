# Table of Contents
1. [Dependencies](#dependencies)
2. [Build](#build)
3. [Install with CMake](#install-with-cmake)
4. [Install with Spack](#install-with-spack)
5. [Using CSG-EB in your project](#using-csg-eb-in-your-project)

This library is for parsing [Constructive Solid Geometry
(CSG)](https://github.com/openscad/openscad/wiki/CSG-File-Format) for use in
[AMReX Embedded Boundaries
(EB)](https://amrex-codes.github.io/amrex/docs_html/EB.html).

CSG is a file format used by the open source CAD modelling tool
[OpenSCAD](https://www.openscad.org). CSG is a subset of the SCAD language,
equivalent in the kinds of geometries that can be expressed, but with fewer
primitives (OpenSCAD can export `.scad` files to `.csg`). SCAD code is easier to
read and write by humans; CSG files are easier to parse and process by software.


## Dependencies

If installing with [Spack](#spack), dependencies will be installed automatically.

Install the following dependencies before building `csg-eb`

 - C++17 compiler (GCC >=7.x, Clang >=5.x)
 - [CMake](https://cmake.org) >=3.14
 - [PEGTL](https://github.com/taocpp/PEGTL) for parsing
 - [Catch2](https://github.com/catchorg/Catch2) for testing framework
 - [CGAL](https://www.cgal.org) for geometry computation (CGAL 5.x). Note we currently do not support newer CGAL versions >= 6.0


## Build

```shell
> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
> cmake --build build
```

### Run tests

```shell
> cd build
> ctest
```

## Install with CMake

To install `csg-eb` to a location `$CSGEB_HOME`:
``` shell
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$CSGEB_HOME
cmake --build build --target install
```


## Install with Spack

Another way to install csg-eb is with [Spack ](https://spack.readthedocs.io).

To install with spack:
``` shell
spack repo add .spack/repo  # add repo with csg-eb recipe
spack spec csg-eb           # preview install
spack install csg-eb        # install csg-eb (with CGAL support)
spack install csg-eb ~cgal  # install csg-eb (without CGAL support)
spack load csg-eb           # Adds libcsg-eb.a to LD_LIBRARY_PATH and csg.hpp to CPATH
```


## Using CSG-EB in your project

 - Include the header file `csg.hpp`
 - Call `csg::get_csgif` with the path to a CSG file to create an [Implicit Function object](https://amrex-codes.github.io/amrex/docs_html/EB.html#implicit-function).

```cpp
#include <csg.hpp>

// for EB2::makeShop
#include <AMReX_EB2_GeometryShop.H>

...


auto csg_file = "~/my_geometry.csg";
auto is_internal_flow = true;

auto csg_if = csg::get_csgif(csg_file, is_internal_flow);

auto gshop = EB2::makeShop(*csg_if);

```


### Note on the return value

If `is_internal_flow==true`, the CSG geometry defines a hollow space for the domain:

   - csg_if(x, y, z) < 0 ⇒ (x, y, z) is **INSIDE**
   - csg_if(x, y, z) > 0 ⇒ (x, y, z) is **OUTSIDE**

If `is_internal_flow==false`, the CSG geometry defines a solid boundary to the domain:

   - csg_if(x, y, z) < 0 ⇒ (x, y, z) is **OUTSIDE**
   - csg_if(x, y, z) > 0 ⇒ (x, y, z) is **INSIDE**


### Add to the CMake build

* Add this repo as a CMake subdirectory in your `CMakeLists.txt`:

```cmake
add_subdirectory(/path/to/repo/csg-eb)
```

* Add the `csg` target to the CMake target that uses code like in the the above example:

```cmake
target_link_libraries(<target> PRIVATE csg)
```
