# Parameters for app

# MFIX-Exa plot file to process
plotfile = /nfs/home/jmusser/projects/ALCC/alpha001_ar1.8/plt53750

# Sets the max index for looping through multiple plot files. The
# name (and path) are defined by the plotfile keyword. If you
# want to loop over multiple files, set the following to the largest
# plot file index.
#
#plotfile.max_index = int [undefined]

#==============================================================================
# Defines the mean particle diameter which may be needed for
# various calculations.
#
#  * If the following is left undefined, then it is computed
#    using particle radii values in the plot file.
#  * A provided value overrides a computed mean diameter.
#
# particle.mean_diameter = Real [computed / undefined]

#==============================================================================
# Define constant fluid properties used to computing derived quantities
# like the drag coefficient.

fluid.density =   1.0     # Real [undefined]
fluid.viscosity = 1.8e-5  # Real [undefined]

#==============================================================================
# Compute filtered Eulerian quantities
#
# It is assumed that the fluid velocity components (u_g, v_g, w_g) and particle
# radius (radius) and velocity components (velx, vely, velz) are saved in the
# plot file.
#
# Flag to enable Eulerian drag filtering.
#
# filter.Eulerian = bool [false]
#
# Size of filter (number of grid cells in one direction). The whole domain
# should be uniformly divisible by this quantity.
#
# filter.Eulerian.size = int [undefined]
#
# [Optional] minimum filter size. The filter is refined by a factor of
# two until the minimum filter size is reached. This allows the run to
# process multiple filter sizes without recreating Eulerian solids.
#
# filter.Eulerian.min = int [undefined]
#
# Example: Assume the fluid grid is 512 x 128 x 128, and the following
# input settings:
#
# filter.Eulerian = true
# filter.Eulerian.size = 64
# filter.Eulerian.min = 16
#
# The app will initially generatei 32 (8x2x2) filter regions, each containing
# 64^3 cells. Once the filtered values are computed, the filter size is
# reduced to 32, then 16. Each filter size will be stored in its own csv file.
#
filter.Eulerian = true

filter.Eulerian.size = 64

# Options for smoothing deposited particle data:
# ------------------------------------------------------------//
# Type of smoothing operation:
# none     : (default) no smoothing
# constant : constant filter size
# variable : local filter size computed from alpha_p
#
#smoothing.type = none # string [none], constant, variable
#
# Constant filter size .........................................
#
# The filter width can be explicitly defined. The provided value
# specifies the diffusion coefficient to be used as-is. This
# input overrides a defined width (filter.width)
#
# filter.diff_coeff = Real [undefined]
#
# The filter size can also be defined as a function of the mean
# particle diameter. (See particle.mean_diameter)
#
# filter.width = 7.0 # Real [7.0]
#
# Variable filter size .........................................
#
# A local filter size is computed from the solids volume fraction
#
# filter.sample_size = 10.0 # Real [10.0]
#
#
##==============================================================================
## Settings for diffusion solver
## ------------------------------------------------------------//
#diffusion.verbose        = 0 # int [0]
#diffusion.bottom_verbose = 0 # int [0]
#
#diffusion.max_iter        = 500 # int [100]
#diffusion.bottom_max_iter = 100 # int [100]
#
#diffusion.max_fmg_iter   = 0 # int [0]
#diffusion.linop_maxorder = 2 # int [2]
#
#diffusion.agglomeration  = true  # bool [true]
#diffusion.agg_grid_size  = -1
#
#diffusion.consolidation  = true  # bool [true]
#diffusion.semicoarsening = false # bool [false]
#
#diffusion.max_coarsening_level     = 32 # int [32]
#diffusion.max_semicoarsening_level = 32 # int [32]
#
#diffusion.rtol = 1.0e-4  # Real [1.0e-11]
#diffusion.atol = 1.0e-9  # Real [1.0e-14]
#
#diffusion.bottom_solver = "bicgstab" # string ["bicgstab"]
