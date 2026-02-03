#!/bin/bash -lex

RUN_NAME="DEM03"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FJOIN_PAR=$2/fjoin_par
fi

rm -rf ${RUN_NAME}* const_plt* POST_* &> /dev/null

time -p ${MFIX} inputs

${FJOIN_PAR} -f DEM03_par --end 10 --var 2 --format 8 --dt 0.0001 --id 1 -j POST_POS1.dat
${FJOIN_PAR} -f DEM03_par --end 10 --var 2 --format 8 --dt 0.0001 --id 2 -j POST_POS2.dat

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    diff -b -u -I '#.*'  ${test_post_file} $(basename ${test_post_file})
done
