#!/bin/bash
cd ..
AD_SIM_HOME=$(pwd)

# Number of time steps to load at once into RAM for 600 trajectories.
# Higher values make regridding faster at the cost of additional RAM usage.
BUFFER_SIZE="8600"
# Number of bins for longitude and latitude each. Increases RAM usage dramatically where a value of 100 uses about
# 40 GB RAM.
N_BINS="100"
# Lower bound for the time dimension. Other datapoints are dismissed.
MIN_TIME="0"
# Upper bound for the time dimension. Other datapoints are dismissed.
MAX_TIME="259200"
# Lower bound for longitude. Other datapoints are dismissed.
MIN_LON="-68"
# Upper bound for longitude. Other datapoints are dismissed.
MAX_LON="70"
# Lower bound for latitude. Other datapoints are dismissed.
MIN_LAT="17"
# Upper bound for latitude. Other datapoints are dismissed.
MAX_LAT="85"
# The time in seconds that fall within one time bin.
DELTA_TIME="3600"
# Define the number of time steps before the fastest (600 hPa) ascent starts.
  #Datapoints before that are dismissed.
INFLOW_TIME="7200"
# Define the number of time steps after the fastest (600 hPa) ascent ends.
  #Datapoints after that are dismissed.
OUTFLOW_TIME="7200"

OUTPUT_PATH="data/grid/"
if [ ! -d "${OUTPUT_PATH}" ]
then
    mkdir -p "${OUTPUT_PATH}"
fi

# Create a grid with the "normal" time dimension
"${AD_SIM_HOME}"/build/bin/regrid \
  --input_path "${AD_SIM_HOME}"/data/simulation/ \
  --output_path "${AD_SIM_HOME}"/${OUTPUT_PATH}grid.nc \
  --buffer_size ${BUFFER_SIZE} \
  --n_bins ${N_BINS} \
  --min_time ${MIN_TIME} \
  --max_time ${MAX_TIME} \
  --min_lon ${MIN_LON} \
  --max_lon ${MAX_LON} \
  --min_lat ${MIN_LAT} \
  --max_lat ${MAX_LAT} \
  --delta_t ${DELTA_TIME} \
  --inflow_time ${INFLOW_TIME} \
  --outflow_time ${OUTFLOW_TIME}

# Create a grid with dimensions relative to the
# start of the fastest (600 hPa) ascent.
# This also triggers longitude and latitude to be relative to the start of the ascent.
MIN_LON="-43"
MAX_LON="115"
MIN_LAT="-28"
MAX_LAT="46"
"${AD_SIM_HOME}"/build/bin/regrid \
  --input_path "${AD_SIM_HOME}"/data/simulation/ \
  --output_path "${AD_SIM_HOME}"/${OUTPUT_PATH}grid_relative.nc \
  --buffer_size ${BUFFER_SIZE} \
  --n_bins ${N_BINS} \
  --min_time ${MIN_TIME} \
  --max_time ${MAX_TIME} \
  --min_lon ${MIN_LON} \
  --max_lon ${MAX_LON} \
  --min_lat ${MIN_LAT} \
  --max_lat ${MAX_LAT} \
  --delta_t ${DELTA_TIME} \
  --inflow_time ${INFLOW_TIME} \
  --outflow_time ${OUTFLOW_TIME} \
  --relative_time
