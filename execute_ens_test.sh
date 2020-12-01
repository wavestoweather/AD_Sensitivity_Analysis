#!/bin/bash

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
FILENAME="no_exclusions_conv_400_median"

SUFF="_outSat_sbShape_sbConv"
# make release SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE'
OUTPUT_PATH="/data/project/wcb/sim_results/ens_test/${FILENAME}${SUFF}/"
ENSEMBLE_CONFIG="config_ens_minimal.json"
if [ ! -d "$OUTPUT_PATH" ]
then
    mkdir -p "$OUTPUT_PATH"
else
    rm "${OUTPUT_PATH}"*.txt
fi
START_OVER="0"
START_OVER_ENVIRONMENT="1"
FIXED_ITERATION="0"

echo "Running... "
echo "build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
-t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
-i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
-o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj0_MAP_t000000_p001" \
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
