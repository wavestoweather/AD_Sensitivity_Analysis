#!/bin/bash

KEYS=(Misc Cloud)
OUT_PARAMS=(qc qr qv)

INPUT_FILENAME=("/data/project/wcb/wcb_traj_flag_deriv_mult_min/parquet_concat",10 "anotherfile",4)
AMOUNT=(5 10)

OUTPUT_PATH="pics/"
NODE_TASKS="4"
my_parallel="parallel -j $NODE_TASKS"

$my_parallel "echo 'plot_mogon.py {1} {2}' && sleep 5" ::: ${KEYS[@]} ::: ${OUT_PARAMS[@]} &

KEYS=(Graupel Hail)
OUT_PARAMS=(qc qr qv)

$my_parallel "echo 'plot_mogon.py {1} {2}' && sleep 5" ::: ${KEYS[@]} ::: ${OUT_PARAMS[@]} &
wait

# for i in $INPUT_FILENAME; do IFS=","; set -- $i;
#     parallel -j2 "echo 'running $1 {1}' && sleep 5" ::: {0..$2}
# done
# There are 217 combinations
# $my_parallel "python3 plot_mogon.py ${INPUT_FILENAME} {1} {2} ${OUTPUT_PATH}" ::: Misc Cloud ::: qc qr qv
