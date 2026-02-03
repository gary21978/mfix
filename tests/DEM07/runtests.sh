#!/bin/bash -lx

set -euo pipefail

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
    MPIRUN="mpirun -np 4"
else
    MPIRUN=""
fi

rm -rf const_plt* POST_* &> /dev/null
time -p ${MPIRUN} ${MFIX} ${INPUTS}

${FJOIN_PAR} -f DEM07_par --end 25 --var 100 --format 6 --join POST_GRAN_TEMP.NEW

post_dats=POST*.NEW
for result in ${post_dats}; do
    diff -w -B "AUTOTEST/${result}" "${result}"
done
