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

# Show a progressbar every x steps. Set to 0 to deactivate the progressbar
PROGRESSBAR="0"
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

### Compile the code. You may comment this part if you have already done it
make release SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE -DNPROCS='${NTASKS}
### End compiling

INPUT_FILENAME=$1
FOLDER_NAME=${INPUT_FILENAME##*/}
FOLDER_NAME=${FOLDER_NAME%.*}
# Output path to store the result
OUTPUT_PATH="${AD_SIM_HOME}/data/vladiana_keyparams/${FOLDER_NAME}/"

echo "Using {$INPUT_FILENAME} and storing to {$OUTPUT_PATH}"

cleanup(){
    if [ ${1} == "a_HET" ]
    then
        mv ${OUTPUT_PATH}${1}/id0_* ${OUTPUT_PATH}_notPerturbed.nc_wcb
    fi
    find ${OUTPUT_PATH}${1}/ -type f -name "*.nc_wcb" -exec rm {} \;
}


for ENSEMBLE_CONFIG in ${AD_SIM_HOME}/configs/all/*
do
    SUFF=${ENSEMBLE_CONFIG##*/}
    SUFF=${SUFF%.*}
    echo "Running for ${SUFF}"
    # Path where checkpoints for starting ensembles are saved
    NEW_CHECKPOINTS=${AD_SIM_HOME}/tmp/${SUFF}
    if [ ! -d "$NEW_CHECKPOINTS" ]
    then
        mkdir -p "$NEW_CHECKPOINTS"
    fi
    # Path where temporary results aka the files for a single ensemble member are stored
    if [ ! -d ${OUTPUT_PATH}${SUFF} ]
    then
        mkdir -p ${OUTPUT_PATH}${SUFF}
    fi


    ${AD_SIM_HOME}/build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    -o ${OUTPUT_PATH}${SUFF}"/wcb${TARGET_TIME_AFTER_START}_traj0_MAP_t000000_p001" \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -n ${START_TIME} \
    -l ${INPUT_FILENAME} -r ${TRAJ} -g 0 \
    -m ${ENSEMBLE_CONFIG} \
    -h ${NEW_CHECKPOINTS}

    # Run all the ensembles now
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
        for entry in "${NEW_CHECKPOINTS}/execute_id*"
        do
            array_not_contains old_files ${entry} && arr_files=("${arr_files[@]}" "$entry")
        done
        if [ ${#arr_files[@]} == 0 ]
        then
            break
        fi
        for entry in ${arr_files[@]}
        do
            /bin/bash ${entry}
            old_files=("${old_files[@]}" "$entry")
        done
        unset arr_files
    done

    echo "Finished ensembles for ${SUFF}"

    # Merge the output
    python ${AD_SIM_HOME}/scripts/merge_all.py --data_path ${OUTPUT_PATH}${SUFF}/ --store_path ${OUTPUT_PATH}${SUFF}.nc_wcb
    cleanup $SUFF
done
