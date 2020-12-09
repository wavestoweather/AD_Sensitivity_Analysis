#!/bin/bash

make clean
make release SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE'
rm execute_id*
rm checkpoint_id*

# Environmental conditions
SCALING_FACTOR="1.0"
AUTO_TYPE="3"
INPUT_FILENAME="/data/project/wcb/netcdf/vladiana_met_stats/no_exclusions_conv_400_median.nc_wcb"
TIMESTEP="20"
# Do not use progressbar for ensembles. It clutters the output.
PROGRESSBAR="0"
START_TIME="-1000"
WRITE_INDEX="1000"
SNAPSHOT_INDEX="1"
TARGET_TIME_AFTER_START="26000"
FILENAME="compressed_big"

SUFF="_outSat_sbShape_sbConv"

# OUTPUT_PATH="/data/project/wcb/sim_results/ens_test/${FILENAME}${SUFF}/"
OUTPUT_PATH="/data/project/wcb/netcdf/ens_test_sim_result/${FILENAME}${SUFF}/"

ENSEMBLE_CONFIG="config_ens.json"
if [ ! -d "$OUTPUT_PATH" ]
then
    mkdir -p "$OUTPUT_PATH"
else
    rm "${OUTPUT_PATH}"*.nc_wcb
fi
START_OVER="0"
START_OVER_ENVIRONMENT="1"
FIXED_ITERATION="0"

echo "Running... "
echo "build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
-t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
-i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
-o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj0_t000000_p001" \
-e ${START_OVER_ENVIRONMENT} \
-p ${PROGRESSBAR} \
-n ${START_TIME} \
-l ${INPUT_FILENAME} -r 0 -g 0 \
-m ${ENSEMBLE_CONFIG}"
# valgrind --undef-value-errors=no
build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
-t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
-i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
-o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj0_MAP_t000000_p001" \
-e ${START_OVER_ENVIRONMENT} \
-p ${PROGRESSBAR} \
-n ${START_TIME} \
-l ${INPUT_FILENAME} -r 0 -g 0 \
-m ${ENSEMBLE_CONFIG}

declare -a arr_files
declare -a old_files

array_not_contains () {
    local array="$1[@]"
    local seeking=$2
    local in=0
    for element in "${!array}"; do
        if [[ $element == "$seeking" ]]; then
            in=1
            break
        fi
    done
    return $in
}

while true
do
    for entry in "./execute_id*"
    do
        array_not_contains old_files ${entry} && arr_files=("${arr_files[@]}" "$entry")
    done
    if [ ${#arr_files[@]} == 0 ]
    then
        break
    fi
    for entry in ${arr_files[@]}
    do
        ./${entry}
        old_files=("${old_files[@]}" "$entry")
    done
    unset arr_files
done

# python scripts/merge.py ${OUTPUT_PATH} ${OUTPUT_PATH}