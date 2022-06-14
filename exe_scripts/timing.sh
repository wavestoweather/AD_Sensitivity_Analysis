#!/bin/bash
cd ..
AD_SIM_HOME=$(pwd)

# Set to the number of threads or CPUs in case you want to run ensemble simulations
NTASKS=4

# The simulation mode determines how multiple processes are used
# ensembles at different time steps and sensitvity analysis 0
# sensitivity analysis 1
# ensembles at different time steps 2
# sensitivity analysis on grids 3
# ensembles at fixed intervals and sensitivity analysis 4
SIMULATION_MODE="1"
AUTO_TYPE="3"
# Start time relative to ascend in seconds
START_TIME="-2800"
# Write to disk every WRITE_INDEX steps. If you experience I/O issues,
# then set it to a multiple of SNAPSHOT_INDEX
WRITE_INDEX="20000"
# Time step in seconds
TIMESTEP="20"
# Get results every SNAPSHOT_INDEX steps, in this example every $TIMESTEP (20) seconds
SNAPSHOT_INDEX="1"
# Fix pressure, temperature and ascent during microphysics (=1)
# or take any changes of them into account (=0)
FIXED_ITERATION="0"
# Warm up time in seconds for the model before any data is tracked
WARM_UP="1800"
# Wether to take the data from the netcdf-file every 20 seconds (=1)
# or just use the initial pressure, temperature and ascent
# from the file and simulate microphysics
# until the target time is reached (=0 not recommended)
START_OVER_ENVIRONMENT="1"
FILENAME="no_exclusions_conv_400_median"
# A dummy path
OUTPUT_PATH="${AD_SIM_HOME}/data/dummy/"
if [ ! -d "$OUTPUT_PATH" ]
then
    mkdir -p "$OUTPUT_PATH"
else
    rm "${OUTPUT_PATH}"*.nc_wcb
fi
# A progressbar would be a nuisance if want to measure time
PROGRESSBAR="0"
INPUT_FILENAME="${AD_SIM_HOME}/data/vladiana_trajectories/${FILENAME}.nc_wcb"

TARGET_TIME_AFTER_START="26000"

mpirun -n ${NTASKS} build/bin/./timing \
-w ${WRITE_INDEX} \
-a ${AUTO_TYPE} \
-t ${FIXED_ITERATION} \
-f ${TARGET_TIME_AFTER_START} \
-d ${TIMESTEP} \
-i ${SNAPSHOT_INDEX} \
-b ${SIMULATION_MODE} \
-o ${OUTPUT_PATH}".nc_wcb" \
-e ${START_OVER_ENVIRONMENT} \
-p ${PROGRESSBAR} \
-g "0" \
-l ${INPUT_FILENAME} \
-r 0 \
-u ${WARM_UP}
