MFIX = ../../
AMREX_HOME ?= $(MFIX)/subprojects/amrex/

TOP = .

EBASE := mfix.unit-test.chem

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

USE_ASSERTION = TRUE
TEST = TRUE

#CPPFLAGS += -DCHEM_TEST_EULERIAN01

include ./source/Make.chem
