#!/bin/bash

# Set to the number of threads
NTASKS=6

# Environmental conditions
SCALING_FACTOR="1.0"
AUTO_TYPE="3"
# Timestep, general, input data
INPUT_FILENAME="/data/project/wcb/netcdf/traj_new_appended/traj_t000000_p001.nc_wcb"
TIMESTEP="0.01"
# Get results every SNAPSHOT_INDEX steps, in this example every 2 seconds
SNAPSHOT_INDEX="200"

# Write to disk every WRITE_INDEX steps. If you experience I/O issues,
# then set it to a multiple of SNAPSHOT_INDEX
WRITE_INDEX="200"

# End time in seconds for the simulation. Should be lower than
# the maximum time of the netcdf-file
TARGET_TIME="60"

# Where to write output files. Keep the naming convention of the files, ie
# wcb403220_traj0_MAP...
OUTPUT_PATH="/data/project/wcb/sim_results/sample_vladiana_test/"

# Wether to take the data from the netcdf-file every 20 seconds (=1)
# or just use the initial conditions from the file and simulate microphysics
# until the target time is reached (=0 not recommended)
START_OVER="1"

# Fix pressure, temperature and ascent during microphysics (=1)
# or take any changes of them into account (=0)
FIXED_ITERATION="0"

# Run the stuff (in parallel with GNU parallel)
parallel -j ${NTASKS} --no-notice --delay .2 build/apps/src/microphysics/./trajectories -w ${WRITE_INDEX} -a ${AUTO_TYPE} \
-t ${FIXED_ITERATION} -s ${START_OVER} -f ${TARGET_TIME} -d ${TIMESTEP} \
-i ${SNAPSHOT_INDEX} -b ${SCALING_FACTOR} \
-o ${OUTPUT_PATH}"wcb${TARGET_TIME}_traj{1}_MAP_t000000_p001" \
-l ${INPUT_FILENAME} -r {1} ::: {0..3}

# Do some post processing; may use a lot of memory with default settings.
# If your machine runs out of memory now (or during visualization),
# try "netcdf"
# DISCLAIMER: netcdf is not yet supported!
FILE_TYPE="parquet"
INPUT_PATH="/data/project/wcb/sim_results/sample_vladiana_test"
STORE_PATH="/data/project/wcb/parquet/sample_vladiana_test"
cd scripts
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}