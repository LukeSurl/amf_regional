#!/bin/sh
#BSUB -q lotus-smp
#BSUB -J "amf_calc"
#BSUB -W 4:00
#BSUB -n 8
#BSUB -u nope@nope.nope
#BSUB -o output.out
#BSUB -e error.err
#BSUB -R "span[ptile=8]"

ulimit -s 999999
export OMP_NUM_THREADS=8
./test_scia.sh 12 04 20 20 1 > amf_calc.log
