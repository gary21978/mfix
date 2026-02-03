#!/bin/bash -lex

RUN_NAME="DEM04"

cleanup() {
    rm -rf ${RUN_NAME}* const_plt* tmp.out &> /dev/null
}

write_data() {
    echo "   " >> $1
    echo "   Friction coefficient $2  (N/m)" >> $1
    echo "   " >> $1

    ${FJOIN_PAR} -f DEM04_par --end 50 --var $3 --var $4 --format 6 --dt 0.0004 -j tmp.out
    cat tmp.out >> $1

}

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FJOIN_PAR=$2/fjoin_par
fi

rm -rf POST_* &> /dev/null

for DES_MEW in 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do

    cleanup

    time -p ${MFIX} inputs \
       "dem.friction_coeff.pp=${DES_MEW}" \
       "dem.friction_coeff.pw=${DES_MEW}"

    write_data POST_VEL.dat ${DES_MEW}  6 11

done

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    diff -b -u -I '#.*'  ${test_post_file} $(basename ${test_post_file})
done
