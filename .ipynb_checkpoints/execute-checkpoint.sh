#!/bin/bash
for execute in trajectories_sb.out # trajectories_sb_noice.out trajectories_sb.out
do
    #
    # Warm cloud parcel
    #

    # Environmental conditions
    SCALING_FACTOR="1.0"
    AUTO_TYPE="3"
    # Timestep, general
    OUTPUT_FILENAME="OUTPUT"
    INPUT_FILENAME="O_WCB_all_20160922_00.nc"
    TIMESTEP="0.01"
    SNAPSHOT_INDEX="1"
    TRAJ_IDX="485"

    # Case 2:
    # Downdraft with existing cloud

    TARGET_TIME="61.0"

    # Take trajectory input and simulate everything
    # START_OVER="0"
    # FIXED_ITERATION="0"
    # ./${execute} -a ${AUTO_TYPE} -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME} -d ${TIMESTEP} -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} -o ${OUTPUT_FILENAME} -l ${INPUT_FILENAME}
    #####
    # Take trajectory input every 20 seconds and simulate everything
    OUTPUT_FILENAME="wcb61_traj${TRAJ_IDX}_start_over_20160922_00"
    START_OVER="1"
    FIXED_ITERATION="0"
    ./${execute} -a ${AUTO_TYPE} -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME} -d ${TIMESTEP} -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} -o ${OUTPUT_FILENAME} -l ${INPUT_FILENAME} -r ${TRAJ_IDX}
    # likwid-perfctr -g DATA -C S0:0-3 -m ./${execute} -a ${AUTO_TYPE} -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME} -d ${TIMESTEP} -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} -o ${OUTPUT_FILENAME} -l ${INPUT_FILENAME}

    # # Fix trajectory input every 20 seconds
    # START_OVER="1"
    # FIXED_ITERATION="1"
    # OUTPUT_FILENAME="OUTPUT_start_over_fixed"
    # ./${execute} -a ${AUTO_TYPE} -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME} -d ${TIMESTEP} -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} -o ${OUTPUT_FILENAME} -l ${INPUT_FILENAME}

    # Set input of trajectory fixed and simulate cloud
    # START_OVER="0"
    # FIXED_ITERATION="1"
    # OUTPUT_FILENAME="tmp/OUTPUT_fixed"
    # ./${execute} -a ${AUTO_TYPE} -t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME} -d ${TIMESTEP} -i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} -o ${OUTPUT_FILENAME} -l ${INPUT_FILENAME}

done
