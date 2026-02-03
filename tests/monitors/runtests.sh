#!/bin/bash -exl

set -euo pipefail

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
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

rm -rf const_plt* POST_* &> /dev/null
time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

post_csvs=*.csv
for result in ${post_csvs}; do
    diff -b -u -I '#.*' "AUTOTEST/${result}" "${result}"
done
