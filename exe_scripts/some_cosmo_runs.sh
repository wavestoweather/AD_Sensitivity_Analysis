#!/bin/bash


TRAJ=0

# Set to the number of threads or CPUs in case you want to run ensemble simulations
NTASKS=4

# Environmental conditions
SCALING_FACTOR="1.0"
AUTO_TYPE="3"
TIMESTEP="20"

# Output variables
PROGRESSBAR="5"
START_TIME="-2800"
WRITE_INDEX="1000"
SNAPSHOT_INDEX="1"
TARGET_TIME_AFTER_START="26000"

# PROGRESSBAR="0"
# TARGET_TIME_AFTER_START="6400"

# Operational variables
START_OVER="0"
START_OVER_ENVIRONMENT="1"
FIXED_ITERATION="0"

OUTPUT_PATH="/data/project/wcb/netcdf/sim_comp/"
IN_PATH="/data/project/wcb/netcdf/sim_comp_input/"

cd ..
AD_SIM_HOME=$(pwd)
NEW_CHECKPOINTS=${AD_SIM_HOME}/tmp

for FILENAME in "conv_400_0_traj_t000000_p001_met3d" "conv_400_10_traj_t000600_p001_met3d" "conv_600_20_traj_t001200_p001_met3d" "conv_600_40_traj_t002400_p001_met3d"
do

    if [ ! -d ${OUTPUT_PATH}${FILENAME} ]
    then
        mkdir -p ${OUTPUT_PATH}${FILENAME}
    fi
    INPUT_FILENAME=${IN_PATH}${FILENAME}

    ${AD_SIM_HOME}/build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
        -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
        -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
        -o ${OUTPUT_PATH}${FILENAME}"/wcb${TARGET_TIME_AFTER_START}_traj0_MAP_t000000_p001" \
        -e ${START_OVER_ENVIRONMENT} \
        -p ${PROGRESSBAR} \
        -n ${START_TIME} \
        -l ${INPUT_FILENAME}.nc_wcb -r ${TRAJ} -g 0
    mv ${OUTPUT_PATH}${FILENAME}"/wcb${TARGET_TIME_AFTER_START}_traj0_MAP_t000000_p001.nc_wcb" ${OUTPUT_PATH}${FILENAME}".nc_wcb"
    rm -r ${OUTPUT_PATH}${FILENAME}/
done
