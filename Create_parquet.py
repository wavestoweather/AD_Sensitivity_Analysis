#!/bin/env python

#SBATCH -p parallel
#SBATCH -A m2_zdvresearch
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --time=120:00:00
#SBATCH --job-name=parquet_ratio
#SBATCH --output=logs/parquet_ratio_%A.out
#SBATCH --error=logs/parquet_ratio_%A.err
#SBATCH --mem-per-cpu=55G

import os
import sys 
sys.path.append(os.getcwd())

from scripts.loader import load_mult_derivates_directory
from scripts.loader import load_nc, rotate_df, norm_deriv, ratio_deriv
import scripts.loader as loader
import numpy as np
from scripts.Deriv import Deriv
from scripts.Sim import Sim
from timeit import default_timer as timer

direc_path = "/lustre/project/m2_zdvresearch/mahieron/wcb_conv_slan"
store_path = "/lustre/project/m2_zdvresearch/mahieron/wcb_conv_slan_parquet"
filt = False
EPSILON = 0.0
ncpus = 2
c_this = False
trajectories = np.arange(63) 
file_list = []
for f in os.listdir(direc_path):
    file_list.append(os.path.join(direc_path, f))
file_list = np.sort(file_list)

for traj in trajectories:
    idx = np.argwhere(["traj{}_".format(traj) in f for f in file_list]).flatten()
    files = file_list[idx]
    if not files:
        break
    t = timer()
    load_f = []
    ref = ""
    sim = " "
    for f in files:
        if "diff" in f:
            load_f.append(f)
        else:
            if "reference" in f:
                ref = f
            else:
                sim = f
                suffix = f[:-4]
                suffix = suffix.split("_")
                suffix = suffix[-2] + "_" + suffix[-1]
    df_dic_mapped = Deriv(direc=direc_path,
                          filt=filt,
                          EPSILON=EPSILON,
                          trajectories=[traj],
                          file_list=load_f,
                          suffix=suffix,
                          threads=ncpus)
    
    df_sim_mapped = Sim()
    df_sim_mapped.load_file(
        filename=sim,
        sep=",",
        change_ref=True,
        refs=ref)
    t2 = timer()
    print("Loading done in {} s".format(t2-t))
    print("Get ratio of data")
    t = timer()
    df_dic_mapped.calculate_ratios()
    t2 = timer()
    print("ratio done in {} s".format(t2-t), flush=True)
    print("Adding columns for output Parameter results")
    t = timer()
    df_dic_mapped.add_param_values(df_sim_mapped.data)
    t2 = timer()
    print("Adding finished in {} s".format(t2-t), flush=True)
    print("Saving as parquet")
    t = timer()
    df_dic_mapped.to_parquet(store_path, compression="snappy")
    t2 = timer()
    print("Saving done in {} s".format(t2-t), flush=True)
    
print("Finished")