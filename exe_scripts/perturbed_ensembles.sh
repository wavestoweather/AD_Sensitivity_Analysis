#!/bin/bash
cd ..
AD_SIM_HOME=$(pwd)

# Set to the number of threads or CPUs in case you want to run ensemble simulations
NTASKS=2

# The simulation mode determines how multiple processes are used
# ensembles at different time steps and sensitvity analysis 0
# sensitivity analysis 1
# ensembles at different time steps 2
# sensitivity analysis on grids 3
# ensembles at fixed intervals and sensitivity analysis 4
SIMULATION_MODE="4"
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
# If you want to see a progressbar, set the value to 500 or above.
# The sumulation is rather fast, hence we set it off (=0).
# The number indicates how many iterations are done between updates of
# the progressbar. On an Intel i7 there are between 1000 and 4000 steps per
# second.
PROGRESSBAR="0"

for FILENAME in "no_exclusions_conv_400_quan25" "no_exclusions_conv_400_median" "no_exclusions_conv_400_quan75" "no_exclusions_conv_600_quan25" "no_exclusions_conv_600_median"  "no_exclusions_conv_600_quan75"
do
    OUTPUT_PATH="${AD_SIM_HOME}/data/vladiana_ensembles/${FILENAME}/"
    if [ ! -d "$OUTPUT_PATH" ]
    then
        mkdir -p "$OUTPUT_PATH"
    else
        rm "${OUTPUT_PATH}"*.nc_wcb
    fi
    for ENSEMBLE_CONFIG in ${AD_SIM_HOME}/configs/all/*
    do
        SUFF=${ENSEMBLE_CONFIG##*/}
        SUFF=${SUFF%.*}

        INPUT_FILENAME="${AD_SIM_HOME}/data/vladiana_trajectories/${FILENAME}.nc_wcb"
        TARGET_TIME_AFTER_START="27800"

        echo "###################################"
        echo "Running for ${INPUT_FILENAME} while perturbing ${SUFF}"
        echo ""
        # optional in case of problems occur: mpirun --mca osc pt2pt
        mpirun --mca osc pt2pt -n ${NTASKS} build/bin/./trajectories \
        -w ${WRITE_INDEX} \
        -a ${AUTO_TYPE} \
        -t ${FIXED_ITERATION} \
        -f ${TARGET_TIME_AFTER_START} \
        -d ${TIMESTEP} \
        -i ${SNAPSHOT_INDEX} \
        -b ${SIMULATION_MODE} \
        -o ${OUTPUT_PATH}${SUFF}".nc_wcb" \
        -e ${START_OVER_ENVIRONMENT} \
        -p ${PROGRESSBAR} \
        -l ${INPUT_FILENAME} \
        -n ${START_TIME} \
        -r 0 \
        -m ${ENSEMBLE_CONFIG} \
        -u ${WARM_UP}
        break
    done
    break
done
