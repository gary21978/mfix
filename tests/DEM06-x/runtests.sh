#!/bin/bash -lx

set -euo pipefail

RUN_NAME="DEM06"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FJOIN_PAR=$2/fjoin_par
fi

INPUTS=inputs_single
if [ -n "$3" ]; then
    INPUTS=$3
fi
echo "Using INPUTS file ${INPUTS}"

if [ "$MFIX_MPI" -eq "1" ]; then
    if [ "$MFIX_OMP" -eq "1" ]; then
  MPIRUN="mpirun -np 2"
    else
  MPIRUN="mpirun -np 4"
    fi
else
    MPIRUN=""
fi

rm -rf ${RUN_NAME}* const_plt* POST_* &> /dev/null

time -p ${MPIRUN} ${MFIX} "${INPUTS}" "amr.check_input=0" "mfix.stop_time=0.150"
time -p ${MPIRUN} ${MFIX} "${INPUTS}" "amr.check_input=0" "mfix.restart=DEM06_chk00150"

${FJOIN_PAR} -f DEM06_par --end 350 --var  1 --format 4 --dt 0.001 -j POST_POS.NEW
${FJOIN_PAR} -f DEM06_par --end 350 --var  6 --format 4 --dt 0.001 -j POST_VEL.NEW

post_dats=POST*.NEW
for result in ${post_dats}; do
    diff -b -u -I '#.*'  -w -B "../DEM06-y/AUTOTEST/${result}" "${result}"
done
