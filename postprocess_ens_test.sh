SUFF="_outSat_sbShape_sbConv"
FILENAME="no_exclusions_conv_400_median"
OUTPUT_PATH="/data/project/wcb/sim_results/ens_test/${FILENAME}${SUFF}/"
INPUT_PATH="/data/project/wcb/sim_results/ens_test/${FILENAME}${SUFF}"
STORE_PATH="/data/project/wcb/netcdf/ens_test/${FILENAME}${SUFF}"
FILE_TYPE="netcdf"
INPUT_TYPE="MET3D"
INPUT_FILENAME="/data/project/wcb/netcdf/vladiana_met_stats/no_exclusions_conv_400_median.nc_wcb"
CHANGE_IDS="True"

if [ ! -d "$STORE_PATH" ]
then
    mkdir -p "$STORE_PATH"
else
    rm "${STORE_PATH}/"*.nc_wcb
fi

cd scripts
python Create_processed_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH} ${INPUT_TYPE} ${CHANGE_IDS}

# Plot results
# Define store path and which plots to generate
# See the file to alter it to your desire.
python plot_outcomes.py ${STORE_PATH} ../pics/ ${SUFF} ${INPUT_FILENAME}
cd ..