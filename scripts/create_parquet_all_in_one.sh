FILE_TYPE="netcdf"
INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/slan_600_quan25"
STORE_PATH="/data/project/wcb/netcdf/sim_output"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/slan_600_median"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/slan_600_quan75"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/slan_400_quan25"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/slan_400_median"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/slan_400_quan75"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/conv_400_quan25"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/conv_400_median"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/conv_400_quan75"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/conv_600_quan25"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/conv_600_median"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}

INPUT_PATH="/data/project/wcb/sim_results/traj_stats2/conv_600_quan75"
python Create_parquet_local.py ${FILE_TYPE} ${INPUT_PATH} ${STORE_PATH}