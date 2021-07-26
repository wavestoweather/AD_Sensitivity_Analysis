#!/bin/bash
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
# Values that are too large or too small may decrease performance.
WRITE_INDEX="1000"
TIMESTEP="20"
# Get results every SNAPSHOT_INDEX steps, in this example every $TIMESTEP (20) seconds
SNAPSHOT_INDEX="1"
# Fix pressure, temperature and ascent during microphysics (=1)
# or take any changes of them into account (=0)
FIXED_ITERATION="0"
#
TRACK_FILE="configs/qr_config.json"
# Wether to take the data from the netcdf-file every 20 seconds (=1)
# or just use the initial pressure, temperature and ascent
# from the file and simulate microphysics
# until the target time is reached (=0 not recommended)
START_OVER_ENVIRONMENT="1"

#  "no_exclusions_conv_400_quan25" "no_exclusions_conv_400_quan75" \
# "no_exclusions_conv_600_median" "no_exclusions_conv_600_quan25" "no_exclusions_conv_600_quan75" \
# "no_exclusions_slan_400_median" "no_exclusions_slan_400_quan25" "no_exclusions_slan_400_quan75" \
# "no_exclusions_slan_600_median" "no_exclusions_slan_600_quan25" "no_exclusions_slan_600_quan75"
SUFF="_mpi_version_ensemble"

# for FILENAME in "no_exclusions_conv_400_median"
for FILENAME in "conv_400_0_traj_t000000_p001_met3d"
do
    PROGRESSBAR="500"
    INPUT_FILENAME="/data/project/wcb/netcdf/vladiana_met_updated/${FILENAME}.nc_wcb"

    TARGET_TIME_AFTER_START=$(ncdump -h $INPUT_FILENAME | grep -m 1 "time = " | sed 's/[^0-9]//g' )
    let TARGET_TIME_AFTER_START=$TARGET_TIME_AFTER_START*20-20
    # TARGET_TIME_AFTER_START="403180"

    echo "###################################"
    echo "Running for ${INPUT_FILENAME} until ${TARGET_TIME_AFTER_START}"
    echo ""

    # Where to write output files. Keep the naming convention of the files, ie
    # wcb403220_traj0_MAP...
    OUTPUT_PATH="/data/project/wcb/netcdf/sim_result_met3d_compatible/${FILENAME}${SUFF}/"
    if [ ! -d "$OUTPUT_PATH" ]
    then
        mkdir -p "$OUTPUT_PATH"
    else
        rm "${OUTPUT_PATH}"*.nc_wcb
    fi

    mpirun -n ${NTASKS} build/bin/./trajectories \
    -w ${WRITE_INDEX} \
    -a ${AUTO_TYPE} \
    -t ${FIXED_ITERATION} \
    -s ${TRACK_FILE} \
    -f ${TARGET_TIME_AFTER_START} \
    -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} \
    -b ${SIMULATION_MODE} \
    -o ${OUTPUT_PATH}"qr_wcb${TARGET_TIME_AFTER_START}_traj0_MAP_t000000_p001" \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -g "0" \
    -l ${INPUT_FILENAME} \
    -r 0
done
