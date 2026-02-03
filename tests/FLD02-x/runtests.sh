#!/bin/bash -lex

set -euo pipefail

RUN_NAME="FLD02"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FEXTRACT=$2/amrex_fextract
fi
if [ -z "${FEXTRACT}" ]; then
    echo "FEXTRACT is not set. Aborting test."
    exit 1
fi

INPUTS=inputs_single
if [ -n "$3" ]; then
    INPUTS=$3
fi
echo "Using INPUTS file ${INPUTS}"

if [ "$MFIX_MPI" -eq "1" ]; then
    MPIRUN="mpirun -np 2"
else
    MPIRUN=""
fi

rm -rf POST_* const_plt* ${RUN_NAME}* &> /dev/null
time -p ${MPIRUN} "${MFIX}" "${INPUTS}"

${FEXTRACT} -d 1 -v u_g -e -p 8 -t 1.0e-10 -s POST_VG.dat FLD0200001

post_dats=POST*.dat
for result in ${post_dats}; do
    diff -u -I '#.*' "../FLD02-y/AUTOTEST/${result}" "${result}"
done
