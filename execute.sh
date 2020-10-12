#!/bin/bash

# Set to the number of threads
NTASKS=4

# Environmental conditions
SCALING_FACTOR="1.0"
AUTO_TYPE="3"
# Timestep, general, input data _testset/conv_400_median_ansatz4.nc_wcb
INPUT_FILENAME="/data/project/wcb/netcdf/traj_stats2_testset/conv_400_median_long.nc_wcb"
# INPUT_FILENAME="/data/project/wcb/netcdf/vladiana_1_test/no_exclusions_conv_600_median.nc_wcb"
INPUT_FILENAME="/data/project/wcb/netcdf/vladiana_met/conv_400_0_traj_t000000_p001.nc_wcb"
TIMESTEP="0.01"
# Update the progressbar after that many simulation steps
# Usually there are about 87 to 280 snapshots/second
PROGRESSBAR="300" # Only useful with -u for parallel and 1 process
# TIMESTEP="19.99"
# Get results every SNAPSHOT_INDEX steps, in this example every 2 seconds
# SNAPSHOT_INDEX="200"

# Write to disk every WRITE_INDEX steps. If you experience I/O issues,
# then set it to a multiple of SNAPSHOT_INDEX
# WRITE_INDEX="200"

# End time in seconds for the simulation. Should be lower than
# the maximum time of the netcdf-file
# WRITE_INDEX="1000"
# SNAPSHOT_INDEX="100"
# TARGET_TIME_AFTER_START="5000"
# START_TIME="-1000"

# TARGET_TIME_AFTER_START="600"
START_TIME="-1000"
WRITE_INDEX="1000"
SNAPSHOT_INDEX="100"
TARGET_TIME_AFTER_START="7500" # Minus Start_time
# TARGET_TIME_AFTER_START="400"
# for SUFF in "_outSat_sbConv" "_outSat" "_outSat_sbShape_sbConv" "_outSat_sbShape"
# dt 20 s
# TIMESTEP="20"
# SNAPSHOT_INDEX="1"
# TARGET_TIME_AFTER_START="3000"
# for SUFF in "_outSat_sbConv" "_inSat_sbConv" "_outSat" "_inSat" "_outSat_sbShape_sbConv" "_inSat_sbShape_sbConv" "_outSat_sbShape" "_inSat_sbShape"
#  -DTRACE_TIME -DTRACE_SAT -DTRACE_QV -DTRACE_QC -DTRACE_QR -DTRACE_ENV -DTRACE_QS -DTRACE_QI -DTRACE_QG -DTRACE_QH -DTRACE_SEDI
for SUFF in "_outSat_sbShape_sbConv"
do
    make clean
    echo "###################################"
    echo "Running for ${SUFF}"
    echo ""
    case $SUFF in
        _outSat_sbConv)
            make release SOURCE='-DMET3D -DSB_CONV'
            ;;
        _inSat_sbConv)
            make release SOURCE='-DMET3D -DIN_SAT_ADJ -DSB_CONV'
            ;;
        _outSat)
            make release SOURCE='-DMET3D'
            ;;
        _inSat)
            make release SOURCE='-DMET3D -DIN_SAT_ADJ'
            ;;
        _outSat_sbShape_sbConv)
            make release SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE'
            ;;
        _outSat_sbShape_sbConv_debug)
            make release SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE -DTRACE_SEDI -DTRACE_TIME -DTRACE_SAT -DTRACE_QV -DTRACE_QC -DTRACE_QR -DTRACE_ENV -DTRACE_QS -DTRACE_QI -DTRACE_QG -DTRACE_QH'
            ;;
        _inSat_sbShape_sbConv)
            make release SOURCE='-DMET3D -DIN_SAT_ADJ -DSB_CONV -DSB_SHAPE'
            ;;
        _outSat_sbShape)
            make release SOURCE='-DMET3D -DSB_SHAPE'
            ;;
        _inSat_sbShape)
            make release SOURCE='-DMET3D -DIN_SAT_ADJ -DSB_SHAPE'
            ;;
    esac

# SUFF="_inSat"
    # Where to write output files. Keep the naming convention of the files, ie
    # wcb403220_traj0_MAP...
    OUTPUT_PATH="/data/project/wcb/sim_results/traj_stats2_testset/conv_400_0_t000000_p001_start_over_mult${SUFF}/"
    if [ ! -d "$OUTPUT_PATH" ]
    then
        mkdir "$OUTPUT_PATH"
    else
        rm "${OUTPUT_PATH}"*.txt
    fi
    # Wether to take the data from the netcdf-file every 20 seconds (=1)
    # or just use the initial mixing ratios and particles numbers
    # from the file and simulate microphysics
    # until the target time is reached (=0 not recommended)
    START_OVER="1"

    # Wether to take the data from the netcdf-file every 20 seconds (=1)
    # or just use the initial pressure, temperature and ascent
    # from the file and simulate microphysics
    # until the target time is reached (=0 not recommended)
    START_OVER_ENVIRONMENT="1"

    # Fix pressure, temperature and ascent during microphysics (=1)
    # or take any changes of them into account (=0)
    FIXED_ITERATION="0"
    # Only used if you do not want to use parallel
    # TRAJECTORY="3"

    # # Run the stuff (in parallel with GNU parallel or not)
    # build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    # -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    # -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    # -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj${TRAJECTORY}_MAP_t000000_p001" \
    # -e ${START_OVER_ENVIRONMENT} \
    # -p ${PROGRESSBAR} \
    # -n ${START_TIME} \
    # -l ${INPUT_FILENAME} -r ${TRAJECTORY}

    # parallel -u -j ${NTASKS} --no-notice --delay .2 build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    # -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    # -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    # -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj{1}_MAP_t000000_p001" \
    # -e ${START_OVER_ENVIRONMENT} \
    # -p ${PROGRESSBAR} \
    # -n ${START_TIME} \
    # -l ${INPUT_FILENAME} -r {1} ::: {3..3}

    OUTPUT_PATH="/data/project/wcb/sim_results/traj_stats2_testset/conv_400_0_t000000_p001_mult${SUFF}/"
    START_OVER="0"
    if [ ! -d "$OUTPUT_PATH" ]
    then
        mkdir "$OUTPUT_PATH"
    else
        rm "${OUTPUT_PATH}"*.txt
    fi

    # build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    # -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    # -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    # -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj${TRAJECTORY}_MAP_t000000_p001" \
    # -e ${START_OVER_ENVIRONMENT} \
    # -p ${PROGRESSBAR} \
    # -n ${START_TIME} \
    # -l ${INPUT_FILENAME} -r ${TRAJECTORY}


    parallel -u -j ${NTASKS} --no-notice --delay .2 build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
    -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME_AFTER_START} -d ${TIMESTEP} \
    -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
    -o ${OUTPUT_PATH}"wcb${TARGET_TIME_AFTER_START}_traj{1}_MAP_t000000_p001" \
    -e ${START_OVER_ENVIRONMENT} \
    -p ${PROGRESSBAR} \
    -n ${START_TIME} \
    -l ${INPUT_FILENAME} -r {1} ::: {3..5}

    echo ""
    FILE_TYPE="netcdf"
    INPUT_TYPE="MET3D"
    INPUT_PATH="/data/project/wcb/sim_results/traj_stats2_testset/conv_400_0_t000000_p001_start_over_mult${SUFF}"
    STORE_PATH="/data/project/wcb/netcdf/sim_output_testset/conv_400_0_t000000_p001_start_over_mult${SUFF}"
    if [ ! -d "$STORE_PATH" ]
    then
        mkdir "$STORE_PATH"
    else
        rm "${STORE_PATH}/"*.nc_wcb
    fi
    cd scripts
    # python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH} ${INPUT_TYPE}

    INPUT_PATH="/data/project/wcb/sim_results/traj_stats2_testset/conv_400_0_t000000_p001_mult${SUFF}"
    STORE_PATH="/data/project/wcb/netcdf/sim_output_testset/conv_400_0_t000000_p001_mult${SUFF}"
    if [ ! -d "$STORE_PATH" ]
    then
        mkdir "$STORE_PATH"
    else
        rm "${STORE_PATH}/"*.nc_wcb
    fi
    python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH} ${INPUT_TYPE}
    cd ..
done

# FILE_TYPE="parquet"
# INPUT_PATH="/data/project/wcb/sim_results/traj_stats2_testset/conv_400_0_t000000_p001_start_over_ansatz4${SUFF}"
# STORE_PATH="/data/project/wcb/parquet/traj_stats2_testset/conv_400_0_t000000_p001_start_over_ansatz4${SUFF}"
# INPUT_TYPE="MET3D"
# cd scripll -h ts
# python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH} ${INPUT_TYPE}

# INPUT_PATH="/data/project/wcb/sim_results/traj_stats2_testset/conv_400_0_t000000_p001_ansatz4${SUFF}"
# STORE_PATH="/data/project/wcb/parquet/traj_stats2_testset/conv_400_0_t000000_p001_ansatz4${SUFF}"
# python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH} ${INPUT_TYPE}
# pics/ inSat_sbShape
python plot_outcomes.py pics/ outSat_sbShape_sbConv