#!/bin/bash

# The number of the trajectory in the file for thermodynamics
TRAJ=0

# Set to the number of threads or CPUs in case you want to run ensemble simulations
NTASKS=4

# scaling factor used in a one-moment scheme. It has no effect on the two-moment scheme
SCALING_FACTOR="1.0"
# Type for autoconversion. 1 has not been tested for correctness, 2 skips autoconversion, 3 is the Seifert and Beheng version
AUTO_TYPE="3"
# Time step in seconds for the microphysics
TIMESTEP="20"

# Show a progressbar every 5 steps. Set to 0 to deactivate the progressbar
PROGRESSBAR="5"
# Start time in seconds relative to the start of the ascend
START_TIME="-2800"
# Write every 1000 iteration the results to disk. Increasing it needs more RAM but decreases IO
WRITE_INDEX="1000"
# Create a snapshot every iteration.
SNAPSHOT_INDEX="1"
# How many seconds shall be simulated
TARGET_TIME_AFTER_START="26000"

### Operational variables
# If set to 1: Get mass densities and particle densities every iteration from the trajectory input file
START_OVER="0"
# If set to 1: Get the environment variables (pressure, temperature, vertical ascend) every iteration from the trajectory input file
START_OVER_ENVIRONMENT="1"
# If set to 1: temperature, pressure and vertical ascend are fixed during microphysics
# If set to 0: t, p and w are influenced by microphysics
FIXED_ITERATION="0"


cd ..
AD_SIM_HOME=$(pwd)

# Output path to store the result
OUTPUT_PATH="${AD_SIM_HOME}/data/sim_comp_test/"
# Path to the files that we iterate over in the following loop
IN_PATH="${AD_SIM_HOME}/data/sim_comp_input/"

### Compile the code. You may comment this part if you have already done it
make release SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE -DNPROCS='${NTASKS}
### End compiling

for FILENAME in "conv_400_0_traj_t000000_p001_met3d" "conv_400_10_traj_t000600_p001_met3d" "conv_600_20_traj_t001200_p001_met3d" "conv_600_40_traj_t002400_p001_met3d"
do
    # Create the output folder if it doesn't exist
    if [ ! -d ${OUTPUT_PATH} ]
    then
        mkdir -p ${OUTPUT_PATH}
    fi
    INPUT_FILENAME=${IN_PATH}${FILENAME}
    OUPTUT_FILENAME=${FILENAME/"_traj"}

    ${AD_SIM_HOME}/build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
        -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
        -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
        -o ${OUTPUT_PATH}"_wcb${TARGET_TIME_AFTER_START}_traj0_MAP_t000000_p001_"${OUPTUT_FILENAME} \
        -e ${START_OVER_ENVIRONMENT} \
        -p ${PROGRESSBAR} \
        -n ${START_TIME} \
        -l ${INPUT_FILENAME}.nc_wcb -r ${TRAJ} -g 0
done
