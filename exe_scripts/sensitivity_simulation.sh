#!/bin/bash
# run with
# ./sensitivity_simulation.sh <path-to-folder-or-file>
# where the path starts from data/. To run a sensitivity simulation for
# all files in data/, use
# ./sensitivity_simulation.sh all

cd ..
AD_SIM_HOME=$(pwd)

# Set to the number of threads or CPUs in case you want to run ensemble simulations
NTASKS="6"

# The simulation mode determines how multiple processes are used
# ensembles at different time steps and sensitivity analysis 0
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
TIMESTEP="30"
# Get results every SNAPSHOT_INDEX steps, in this example every $TIMESTEP (20) seconds
SNAPSHOT_INDEX="1"
# Wether to take the data from the netcdf-file every time step if available (=1)
# or just use the initial pressure, temperature and ascent
# from the file and simulate microphysics
# until the target time is reached (=0 not recommended)
START_OVER_ENVIRONMENT="1"
# Fix pressure, temperature and ascent during microphysics (=1)
# or take any changes of them into account (=0)
FIXED_ITERATION="0"
# Warm up time in seconds for the model before any data is tracked
WARM_UP="180"
# Start index along the time dimension at which one should start the
# simulation
# An alternative way would be the option '-n START_TIME', where
# START_TIME is the start time in seconds relative to the start of the
# (convective or slantwise) ascend.
START_TIME_IDX="0"
# The number indicates how many iterations are done between updates of
# the progressbar. On an Intel i7 there are between 80 and 200 steps per
# second. Use "0" to deactivate the progressbar.
PROGRESSBAR="500"
# Configuration file that defines a subset of parameters that shall be tracked
# We use this here to track only few parameters which speeds up the simulation.
TRACK_FILE="configs/b8_config_many.json"

OUTPUT_PATH="data/simulation/"
if [ ! -d "${OUTPUT_PATH}" ]
then
    mkdir -p "${OUTPUT_PATH}"
fi
INPUT_FILEPATH="data/"

function run_sim {
    FILENAME=$1
    INPUT_FILENAME="${INPUT_FILEPATH}${FILENAME}.nc"
    TARGET_TIME_AFTER_START=$(ncdump -h "$INPUT_FILENAME" | grep -m 1 "time = " | sed 's/[^0-9]//g' )
    TARGET_TIME_AFTER_START=$(($TARGET_TIME_AFTER_START * ${TIMESTEP} - 10 * ${TIMESTEP}))

    echo "###################################"
    echo "Running for ${INPUT_FILENAME} until ${TARGET_TIME_AFTER_START}"
    echo ""

    mpirun -n ${NTASKS} ${AD_SIM_HOME}/build/bin/./trajectories \
    -w ${WRITE_INDEX} \
    -a ${AUTO_TYPE} \
    -t ${FIXED_ITERATION} \
    -f ${TARGET_TIME_AFTER_START} \
    -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} \
    -b ${SIMULATION_MODE} \
    -o ${OUTPUT_PATH}"${FILENAME}".nc \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -g ${START_TIME_IDX} \
    -l "${INPUT_FILENAME}" \
    -s ${TRACK_FILE} \
    -u ${WARM_UP}

    echo "Finished run for ${FILENAME}"
}

for var in "$@"
do
  if [ "$var" == "all" ]
  then
      echo "Running for all files in ${INPUT_FILEPATH}"
      for INPUT in ${INPUT_FILEPATH}*.nc
      do
          FILENAME=${INPUT##*/}
          FILENAME=${FILENAME%.*}
          run_sim "$FILENAME"
      done
  else
      # Do it for a specific file
      run_sim "$var"
  fi
done
wait
exit 0
