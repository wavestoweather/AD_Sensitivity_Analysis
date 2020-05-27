import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
import progressbar
from loader import *
import os

direc = sys.argv[1]

for filename in os.listdir(direc):
    # print("Loading from {}".format(direc + filename))
    ds = xr.open_dataset(direc + filename)
    n_traj = ds.dims["id"]
    target_time = int(ds.dims["time"])*20
    shortened_name = filename[5:-7]
    print('TARGET_TIME="{}"'.format(target_time))
    print('INPUT_FILENAME="/lustre/project/m2_jgu-tapt/cosmo_output/vladiana/traj/{}"'.format(filename))
    print('$my_parallel "$my_srun build/apps/src/microphysics/./trajectories -w ${{WRITE_INDEX}} -a ${{AUTO_TYPE}} -t ${{FIXED_ITERATION}} -s ${{START_OVER}} -f ${{TARGET_TIME}} -d ${{TIMESTEP}} -i ${{SNAPSHOT_INDEX}} -b ${{SCALING_FACTOR}} -o /lustre/project/m2_zdvresearch/mahieron/wcb_traj_flag_deriv_mult_min2_p002/wcb${{TARGET_TIME}}_traj{{}}_MAP_{} -l ${{INPUT_FILENAME}} -r {{}}" ::: {{0..{}}}'.format(shortened_name, n_traj))
    print("\n")
#     text =(
# """
# TARGET_TIME="{}"
# INPUT_FILENAME="/lustre/project/m2_jgu-tapt/cosmo_output/vladiana/traj/{}"
# $my_parallel "$my_srun build/apps/src/microphysics/./trajectories -w $\{WRITE_INDEX\} -a $\{AUTO_TYPE\} -t $\{FIXED_ITERATION\} -s $\{START_OVER\} -f $\{TARGET_TIME\} -d $\{TIMESTEP\} -i $\{SNAPSHOT_INDEX\} -b $\{SCALING_FACTOR\} -o /lustre/project/m2_zdvresearch/mahieron/wcb_traj_flag_deriv_mult_min2_p002/wcb$\{TARGET_TIME\}_traj\{\}_MAP_{} -l $\{INPUT_FILENAME\} -r \{\}" ::: \{0..{}\}


# """.format(target_time, filename, shortened_name, n_traj))
#     print(text)
    # exit()