AMREX_HOME = ../../subprojects/amrex
AMREX_HYDRO_HOME = ../../subprojects/AMReX-Hydro

HYPRE_DIR = /nfs/home/jmusser/packages/mws015/gcc/9.5.0/hypre/release/v2.29.0/openmpi/4.1.5/cuda/NONE
HYPRE_LIBRARIES=$(HYPRE_DIR)/lib
HYPRE_INCLUDE_DIRS=$(HYPRE_DIR)/include

HDF5_HOME ?= /nfs/home/2/jmusser/packages/joule2/gnu/9.3.0/openmpi/4.0.4/hdf5/1.14.2
ZFP_HOME ?=  /nfs/home/2/jmusser/packages/joule2/gnu/9.3.0/openmpi/4.0.4/zfp/1.0.0
H5Z_HOME  ?= /nfs/home/2/jmusser/packages/joule2/gnu/9.3.0/openmpi/4.0.4/H5Z-ZFP/v1.1.1

TOP = .

EBASE := post-mfix

BL_NO_FORT = TRUE

DEBUG    = TRUE
DEBUG    = FALSE

DIM = 3

COMP = gcc

USE_PARTICLES = TRUE

PRECISION = DOUBLE

USE_MPI  = TRUE
USE_OMP  = FALSE
USE_CUDA = FALSE

USE_HDF5  = FALSE
USE_HDF5_ZFP = FALSE

USE_HYPRE = FALSE

include ./post-mfix/Make.post_mfix
