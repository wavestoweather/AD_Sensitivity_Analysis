#!/bin/bash
cd ..
# Set to the number of threads or CPUs in case you want to run ensemble simulations
NTASKS=2

# Environmental conditions
SCALING_FACTOR="1.0"
AUTO_TYPE="3"
# Relative to ascend
START_TIME="-1000"
# Write to disk every WRITE_INDEX steps. If you experience I/O issues,
# then set it to a multiple of SNAPSHOT_INDEX
WRITE_INDEX="1000"
TIMESTEP="20"
# Get results every SNAPSHOT_INDEX steps, in this example every $TIMESTEP (20) seconds
SNAPSHOT_INDEX="1"
# Fix pressure, temperature and ascent during microphysics (=1)
# or take any changes of them into account (=0)
FIXED_ITERATION="0"
# Wether to take the data from the netcdf-file every 20 seconds (=1)
# or just use the initial mixing ratios and particles numbers
# from the file and simulate microphysics
# until the target time is reached (=0)
START_OVER="0"
# Wether to take the data from the netcdf-file every 20 seconds (=1)
# or just use the initial pressure, temperature and ascent
# from the file and simulate microphysics
# until the target time is reached (=0 not recommended)
START_OVER_ENVIRONMENT="1"

#  "no_exclusions_conv_400_quan25" "no_exclusions_conv_400_quan75" \
# "no_exclusions_conv_600_median" "no_exclusions_conv_600_quan25" "no_exclusions_conv_600_quan75" \
# "no_exclusions_slan_400_median" "no_exclusions_slan_400_quan25" "no_exclusions_slan_400_quan75" \
# "no_exclusions_slan_600_median" "no_exclusions_slan_600_quan25" "no_exclusions_slan_600_quan75"
ENSEMBLE_CONFIG=configs/all/a_ccn_1.json
SUFF="_mpi_version_ensemble"
NEW_CHECKPOINTS=tmp/${SUFF}
if [ ! -d "$NEW_CHECKPOINTS" ]
    then
        mkdir -p "$NEW_CHECKPOINTS"
    else
        rm "${NEW_CHECKPOINTS}"*.json
        rm "${NEW_CHECKPOINTS}"*.sh
    fi

for FILENAME in "no_exclusions_conv_400_median"
do
    if [[ "$FILENAME" == *"conv_"* ]]
    then
        # Update the progressbar after that many simulation steps
        PROGRESSBAR="20"
        # End time in seconds for the simulation. Should be lower than
        # the maximum time of the netcdf-file
        TARGET_TIME_AFTER_START="26000" # Minus Start_time
        TARGET_TIME_AFTER_START="5400"
    else
        # Update the progressbar after that many simulation steps
        PROGRESSBAR="5"
        # End time in seconds for the simulation. Should be lower than
        # the maximum time of the netcdf-file
        TARGET_TIME_AFTER_START="83000" # Minus Start_time
    fi

    INPUT_FILENAME="/data/project/wcb/netcdf/vladiana_met_stats/${FILENAME}.nc_wcb"

    echo "###################################"
    echo "Running for ${INPUT_FILENAME}"
    echo ""

    # Where to write output files. Keep the naming convention of the files, ie
    # wcb403220_traj0_MAP...
    OUTPUT_PATH="/data/project/wcb/netcdf/sim_result/${FILENAME}${SUFF}/"
    if [ ! -d "$OUTPUT_PATH" ]
    then
        mkdir -p "$OUTPUT_PATH"
    else
        rm "${OUTPUT_PATH}"*.nc_wcb
    fi

    mpirun -n ${NTASKS} build/bin/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj0_MAP_t000000_p001" \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -n ${START_TIME} \
    -l ${INPUT_FILENAME} -r 0 -g 0 \
    -m ${ENSEMBLE_CONFIG} \
    -h ${NEW_CHECKPOINTS}
done
