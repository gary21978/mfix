MFIX = ../../
AMREX_HOME ?= $(MFIX)/subprojects/amrex/

TOP = .

EBASE := mfix_fluid.unit-test

BL_NO_FORT = TRUE

DEBUG    = FALSE
DEBUG    = TRUE

DIM = 3

COMP = gcc

USE_PARTICLES = TRUE
USE_EB        = TRUE
USE_MG        = FALSE

PRECISION = DOUBLE

USE_MPI  = FALSE
USE_OMP  = FALSE
USE_CUDA = FALSE

CPPFLAGS += -DMFIX_NO_ABORT -DMFIX_SILENT_REPORTER

include ./source/Make.mfix_fluid
