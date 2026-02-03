#!/bin/bash -lex

RUN_NAME="DEM05"

write_data() {
    echo "   " >> $1
    echo "   Particle $2" >> $1
    echo "   " >> $1

    ${FJOIN_PAR} -f DEM05_par --end 25 --id $2 --var $3 --var $4 --var $5 --format 6 --dt 0.0001 -j tmp.out
    cat tmp.out >> $1

}

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FJOIN_PAR=$2/fjoin_par
fi

rm -f ${RUN_NAME}* const_plt* POST_* &> /dev/null

time -p ${MFIX} inputs

for pID in {1..62}; do
    write_data POST_VEL.dat ${pID}  7 6 11
done

post_dats=../DEM05-y/AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    diff -b -u -I '#.*'  ${test_post_file} $(basename ${test_post_file})
done
