make clean
make scratch SOURCE='-DMET3D -DSB_CONV -DSB_SHAPE'
build/apps/src/scratch/load_test \
/data/project/wcb/netcdf/vladiana_met_stats/no_exclusions_conv_400_median.nc_wcb \
config.xml