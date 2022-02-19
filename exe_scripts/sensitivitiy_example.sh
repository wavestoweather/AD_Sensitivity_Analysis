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
# Write to disk every WRITE_INDEX steps. If you experience I/O issues,
# then set it to a multiple of SNAPSHOT_INDEX
WRITE_INDEX="20000"
# Time step in seconds
TIMESTEP="20"
# Get results every SNAPSHOT_INDEX steps, in this example every $TIMESTEP (20) seconds
SNAPSHOT_INDEX="1"
# Wether to take the data from the netcdf-file every 20 seconds (=1)
# or just use the initial pressure, temperature and ascent
# from the file and simulate microphysics
# until the target time is reached (=0 not recommended)
START_OVER_ENVIRONMENT="1"
START_TIME="-2800"
# The number indicates how many iterations are done between updates of
# the progressbar. On an Intel i7 there are between 1000 and 4000 steps per
# second.
PROGRESSBAR="500"
# Configuration file that defines a subset of parameters that shall be tracked
# We use this here to track only few parameters which speeds up the simulation.
TRACK_FILE="configs/qr_config_all.json"

for FILENAME in "conv_400_0_traj_t000000_p001_met3d" "conv_400_10_traj_t000600_p001_met3d" "conv_600_20_traj_t001200_p001_met3d" "conv_600_40_traj_t002400_p001_met3d"
do

    INPUT_FILENAME="${AD_SIM_HOME}/data/vladiana_input/${FILENAME}.nc_wcb"
    TARGET_TIME_AFTER_START=$(ncdump -h $INPUT_FILENAME | grep -m 1 "time = " | sed 's/[^0-9]//g' )
    TARGET_TIME_AFTER_START=$(($TARGET_TIME_AFTER_START * 20 - 20))
    # TARGET_TIME_AFTER_START=25000

    echo "###################################"
    echo "Running for ${INPUT_FILENAME} until ${TARGET_TIME_AFTER_START}"
    echo ""

    OUTPUT_PATH="data/comparison/"
    if [ ! -d "$OUTPUT_PATH" ]
    then
        mkdir -p "$OUTPUT_PATH"
    fi

    mpirun -n ${NTASKS} build/bin/./trajectories \
    -w ${WRITE_INDEX} \
    -a ${AUTO_TYPE} \
    -f ${TARGET_TIME_AFTER_START} \
    -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} \
    -b ${SIMULATION_MODE} \
    -o ${OUTPUT_PATH}${FILENAME}.nc_wcb \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -n ${START_TIME} \
    -r "0" \
    -l ${INPUT_FILENAME} \
    -s ${TRACK_FILE}
    break
done
