#!/bin/bash
# Run a sensitivity analysis of several ensembles with the top variables for qr
cd ..
# Set to the number of threads or CPUs in case you want to run ensemble simulations
NTASKS=4

# trajectory_sensitvity_perturbance 0
# trajectory_sensitivity 1
# trajectory_perturbance 2
# grid_sensitivity 3
SIMULATION_MODE="1"
AUTO_TYPE="3"

# Write to disk every WRITE_INDEX steps. If you experience I/O issues,
# then set it to a multiple of SNAPSHOT_INDEX
WRITE_INDEX="15000"

TIMESTEP="20"
# Get results every SNAPSHOT_INDEX steps, in this example every $TIMESTEP (20) seconds
SNAPSHOT_INDEX="1"
# Fix pressure, temperature and ascent during microphysics (=1)
# or take any changes of them into account (=0)
FIXED_ITERATION="0"
#
TRACK_FILE="configs/qv_qr_lat_config.json"
# Wether to take the data from the netcdf-file every 20 seconds (=1)
# or just use the initial pressure, temperature and ascent
# from the file and simulate microphysics
# until the target time is reached (=0 not recommended)
START_OVER_ENVIRONMENT="1"

# Warm up time in seconds for the model before any data is tracked
WARM_UP="1800"
# Track sensitivities starting from the initial values.
TRACK_INIT="1"

SUFF="_qv_qr_lat_ic"
INPUT_PATH="/data/project/wcb/netcdf/vladiana_met_updated/"
# INPUT_PATH_TWO="/data/project/wcb/netcdf/vladiana_met/"

OUTPUT_PATH="/data/project/wcb/b8/init_test/"
if [ ! -d "$OUTPUT_PATH" ]
then
    mkdir -p "$OUTPUT_PATH"
fi

for INPUT_FILENAME in ${INPUT_PATH}conv_400_10_traj_t000600_p001_met3d.nc_wcb # conv_400_0_traj_t000000_p001_met3d.nc_wcb
do
    FILENAME=${INPUT_FILENAME##*/}
    FILENAME=${FILENAME%.*}
    PROGRESSBAR="100"

    TARGET_TIME_AFTER_START=$(ncdump -h $INPUT_FILENAME | grep -m 1 "time = " | sed 's/[^0-9]//g' )
    TARGET_TIME_AFTER_START=$(($TARGET_TIME_AFTER_START * 20 - 20))
    echo "Running for ${FILENAME} until ${TARGET_TIME_AFTER_START} with ${NTASKS} tasks"

    mpirun -n ${NTASKS} build/bin/./trajectories \
    -w ${WRITE_INDEX} \
    -a ${AUTO_TYPE} \
    -t ${FIXED_ITERATION} \
    -f ${TARGET_TIME_AFTER_START} \
    -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} \
    -b ${SIMULATION_MODE} \
    -o ${OUTPUT_PATH}${FILENAME}${SUFF}.nc \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -g "0" \
    -l ${INPUT_FILENAME} \
    -s ${TRACK_FILE} \
    -u ${WARM_UP} \
    -x ${TRACK_INIT}
done
