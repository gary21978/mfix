#!/bin/bash
#BSUB -P CFD122
#BSUB -W 0:20
#BSUB -nnodes 1
#BSUB -J MFIX
#BSUB -o MFIXo.%J
#BSUB -e MFIXe.%J
module load gcc
module load cuda/9.1.85
module list
set -x
omp=1
export OMP_NUM_THREADS=${omp}
EXE="./mfix3d.gnu.MPI.CUDA.ex"
jsrun -n 1 -a 1 -g 1 -c 1 --bind=packed:${omp} ${EXE} inputs > output.txt
