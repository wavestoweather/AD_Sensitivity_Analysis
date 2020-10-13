#!/bin/bash

# Set to the number of threads or CPUs
NTASKS=4

# Environmental conditions
SCALING_FACTOR="1.0"
AUTO_TYPE="3"
# Define here your input data
INPUT_FILENAME="data/conv_400_0_traj_t000000_p001.nc_wcb"
# Time step size in seconds for the microphysics
TIMESTEP="0.01"
# Update the progressbar after that many simulation steps
# Usually there are about 87 to 280 snapshots/second
# Only useful with -u for parallel or 1 process
PROGRESSBAR="300"
# Get results every SNAPSHOT_INDEX steps, in this example
# every $TIMESTEP*SNAPSHOT = 2 seconds
SNAPSHOT_INDEX="200"
# Write to disk every WRITE_INDEX steps. If you experience I/O issues,
# then set it to a multiple of SNAPSHOT_INDEX
WRITE_INDEX="1000"
# Start time in seconds before the ascend starts
START_TIME="-100"
# End time after START_TIME in seconds for the simulation.
TARGET_TIME_AFTER_START="250"

# You can use the loop to run different versions with hard coded parameters
# or with debugging statements. See the cases below for different versions
for SUFF in "_outSat_sbShape_sbConv"
do
    make clean
    echo "###################################"
    echo "Running for ${SUFF}"
    echo ""
    case $SUFF in
        _outSat_sbConv)
            # Saturation adjustment every 20s after microphysics
            # Seifert & Beheng convection parameters
            make release SOURCE='-DMET3D -DSB_CONV'
            ;;
        _inSat_sbConv)
            # Saturation adjustment within microphysics (not recommended, just for educative purpose)
            # Seifert & Beheng convection parameters
            make release SOURCE='-DMET3D -DIN_SAT_ADJ -DSB_CONV'
            ;;
        _outSat)
            # Saturation adjustment every 20s after microphysics
            make release SOURCE='-DMET3D'
            ;;
        _inSat)
            # Saturation adjustment within microphysics (not recommended, just for educative purpose)
            make release SOURCE='-DMET3D -DIN_SAT_ADJ'
            ;;
        _outSat_sbShape_sbConv)
            # Saturation adjustment every 20s after microphysics
            # Seifert & Beheng shape parameters (width and skew)
            # Seifert & Beheng convection parameters
            make release SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE'
            ;;
        _outSat_sbShape_sbConv_debug)
            # Saturation adjustment every 20s after microphysics
            # Seifert & Beheng shape parameters (width and skew)
            # Seifert & Beheng convection parameters
            make release SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE -DTRACE_SEDI -DTRACE_TIME -DTRACE_SAT -DTRACE_QV -DTRACE_QC -DTRACE_QR -DTRACE_ENV -DTRACE_QS -DTRACE_QI -DTRACE_QG -DTRACE_QH'
            ;;
        _inSat_sbShape_sbConv)
            # Saturation adjustment within microphysics (not recommended, just for educative purpose)
            # Seifert & Beheng shape parameters (width and skew)
            # Seifert & Beheng convection parameters
            make release SOURCE='-DMET3D -DIN_SAT_ADJ -DSB_CONV -DSB_SHAPE'
            ;;
        _outSat_sbShape)
            # Saturation adjustment every 20s after microphysics
            # Seifert & Beheng shape parameters (width and skew)
            make release SOURCE='-DMET3D -DSB_SHAPE'
            ;;
        _inSat_sbShape)
            # Saturation adjustment within microphysics (not recommended, just for educative purpose)
            # Seifert & Beheng shape parameters (width and skew)
            make release SOURCE='-DMET3D -DIN_SAT_ADJ -DSB_SHAPE'
            ;;
    esac

    # Where to write output files.
    OUTPUT_PATH="data/sim_results/conv_400_0_t000000_p001_start_over_mult${SUFF}/"
    if [ ! -d "$OUTPUT_PATH" ]
    then
        mkdir -p "$OUTPUT_PATH"
    else
        rm "${OUTPUT_PATH}"*.txt
    fi
    # Wether to take the data from the netcdf-file every 20 seconds (=1)
    # or just use the initial mixing ratios and particles numbers
    # from the file and simulate microphysics
    # until the target time is reached (=0 not recommended)
    START_OVER="1"

    # Wether to take other data from the netcdf-file every 20 seconds (=1)
    # or just use the initial pressure, temperature and ascent
    # from the file and simulate microphysics
    # until the target time is reached (=0 not recommended)
    START_OVER_ENVIRONMENT="1"

    # Fix pressure, temperature and ascent during microphysics (=1)
    # or take any changes of them into account (=0)
    FIXED_ITERATION="0"

    # Only used if you do not want to use parallel
    # TRAJECTORY="0"

    # # Run the stuff (in parallel with GNU parallel or not)
    # build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    # -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    # -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    # -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj${TRAJECTORY}_MAP_t000000_p001" \
    # -e ${START_OVER_ENVIRONMENT} \
    # -p ${PROGRESSBAR} \
    # -n ${START_TIME} \
    # -l ${INPUT_FILENAME} -r ${TRAJECTORY}

    parallel -u -j ${NTASKS} --no-notice --delay .2 build/apps/src/microphysics/./trajectories \
    -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj{1}_MAP_t000000_p001" \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -n ${START_TIME} \
    -l ${INPUT_FILENAME} -r {1} ::: {0..0} # You can run multiple trajectories at once here

    OUTPUT_PATH="data/sim_results/conv_400_0_t000000_p001_mult${SUFF}/"
    START_OVER="0"
    if [ ! -d "$OUTPUT_PATH" ]
    then
        mkdir -p "$OUTPUT_PATH"
    else
        rm "${OUTPUT_PATH}"*.txt
    fi

    # Run without gnu parallel
    # build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    # -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    # -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    # -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj${TRAJECTORY}_MAP_t000000_p001" \
    # -e ${START_OVER_ENVIRONMENT} \
    # -p ${PROGRESSBAR} \
    # -n ${START_TIME} \
    # -l ${INPUT_FILENAME} -r ${TRAJECTORY}

    parallel -u -j ${NTASKS} --no-notice --delay .2 build/apps/src/microphysics/./trajectories \
    -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj{1}_MAP_t000000_p001" \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -n ${START_TIME} \
    -l ${INPUT_FILENAME} -r {1} ::: {0..0} # You can run multiple trajectories at once here
    cd scripts
    echo ""
    # Define file type. We recommend netcdf over parquet
    FILE_TYPE="netcdf"
    # MET3D is a visualization tool. Our most up-to-date version supports it
    INPUT_TYPE="MET3D"
    INPUT_PATH="../data/sim_results/conv_400_0_t000000_p001_start_over_mult${SUFF}"
    STORE_PATH="../data/sim_processed/conv_400_0_t000000_p001_start_over_mult${SUFF}"
    if [ ! -d "$STORE_PATH" ]
    then
        mkdir -p "$STORE_PATH"
    else
        rm "${STORE_PATH}/"*.nc_wcb
    fi
    python Create_processed_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH} ${INPUT_TYPE}

    INPUT_PATH="../data/sim_results/conv_400_0_t000000_p001_mult${SUFF}"
    STORE_PATH="../data/sim_processed/conv_400_0_t000000_p001_mult${SUFF}"
    if [ ! -d "$STORE_PATH" ]
    then
        mkdir -p "$STORE_PATH"
    else
        rm "${STORE_PATH}/"*.nc_wcb
    fi
    python Create_processed_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH} ${INPUT_TYPE}

    # Plot results
    # Define store path and which plots to generate
    # See the file to alter it to your desire.
    python plot_outcomes.py ../pics/ ${SUFF}

    cd ..
done

