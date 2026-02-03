#!/bin/bash -lex

cleanup() {
    rm -rf ${RUN_NAME}* const_plt* tmp.out &> /dev/null
}

write_data() {
    echo "   " >> ${1}
    echo "   Normal collision spring constant. $2  (N/m)" >> $1
    echo "   Restitution coefficient.          $3  ( - )" >> $1
    echo "   " >> ${1}
    ${FJOIN_PAR} -f DEM01_par --end 100 --var $4 --format 8 --dt 0.005 -j tmp.out
    cat tmp.out >> $1
}

RUN_NAME="DEM01"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FJOIN_PAR=$2/fjoin_par
fi

cleanup

echo "DEM01: FREELY FALLING PARTICLE W/WALL COLLISION" > POST_POS.dat
echo "DEM01: FREELY FALLING PARTICLE W/WALL COLLISION" > POST_VEL.dat

DES_KN=10000
for DES_ETA in 0.9 0.8 0.7 0.6; do

    cleanup

    time -p ${MFIX} inputs \
       "dem.spring_const.pp=${DES_KN}        dem.spring_const.pw=${DES_KN}" \
       "dem.restitution_coeff.solid0.solid0=${DES_ETA}  dem.restitution_coeff.solid0.wall=${DES_ETA}"

    write_data "POST_POS.dat" ${DES_KN} ${DES_ETA}  2
    write_data "POST_VEL.dat" ${DES_KN} ${DES_ETA}  7

done

for DES_KN in 25000 50000 100000; do
    for DES_ETA in 1.0 0.9 0.8 0.7 0.6; do

        cleanup

        time -p ${MFIX} inputs \
           "dem.spring_const.pp=${DES_KN}        dem.spring_const.pw=${DES_KN}" \
           "dem.restitution_coeff.solid0.solid0=${DES_ETA}  dem.restitution_coeff.solid0.wall=${DES_ETA}"

        write_data "POST_POS.dat" ${DES_KN} ${DES_ETA}  2
        write_data "POST_VEL.dat" ${DES_KN} ${DES_ETA}  7

    done
done

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    diff -b -u -I '#.*'  ${test_post_file} $(basename ${test_post_file})
done
