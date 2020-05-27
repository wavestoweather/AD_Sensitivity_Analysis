#!/bin/env python

#SBATCH -p parallel
#SBATCH -A m2_zdvresearch
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --job-name=Ascent
#SBATCH --output=logs/Ascent_%A.out
#SBATCH --error=logs/Ascent_%A.err
#SBATCH --mem=100GB

import sys
import os
sys.path.append(os.getcwd())
import multiprocessing
import numpy as np
from pylab import rcParams
from scripts.Deriv_dask import Deriv_dask
from timeit import default_timer as timer
import dask.dataframe as dd

direc_path = "/lustre/project/m2_zdvresearch/mahieron/parquet_complete/wcb_traj_ratio_concat"
store_path = "/lustre/project/m2_zdvresearch/mahieron/parquet_ascent"

ncpus = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
print("Using {} cpus".format(ncpus))

t = timer()
data = Deriv_dask(
    direc=direc_path,
    parquet=True,
    columns=columns,
    backend="matplotlib"
)
t2 = timer()
print("Loading done in {}s".format(t2-t))

df = data.data
t = timer()
trajectories = df.trajectory.unique().compute()
t2 = timer()
print("Computing done in {} s".format(t2-t))
import pandas as pd

window_conv_400 = 1 * 60 * 60
window_conv_600 = 3 * 60 * 60
window_slan_400 = 35 * 6 * 60 
window_slan_400_min = 15 * 6 * 60 
window_slan_600 = 22  * 60 * 60 
window_slan_600_min = 65 * 6 * 60 

differ_400 = lambda x: ((x.max()-x.min()) >= 4000)
differ_600 = lambda x: ((x.max()-x.min()) >= 6000)
t_start = timer()
i = 0

for t in trajectories:
    data_st_1 = df[df.trajectory == t]
    
    for f in data_st_1.netcdf.unique().compute():
        t = timer()
        data_st = data_st_1[data_st_1.netcdf == f]
        i += 1

        # Each entry is 0.01 s
        # We are interested in values of 1.5-3.5 hours
        # and values between 6.5 and 22 hours
        # Calculate the minimum and maximum pressure in those windows
        # Idea use rolling window of maximum time to see, if pressure had been reached
        
        # This is needed for the rolling window stuff
        data_st = data_st.repartition(npartitions=data_st.npartitions // 100)

        conv_400 = data_st["p"].rolling(window_conv_400, min_periods=1).apply(differ_400, raw=True).astype(dtype=bool)
        conv_600 = data_st["p"].rolling(window_conv_600, min_periods=1).apply(differ_600, raw=True).astype(dtype=bool)

        def differ_slan_400(x):
            delta = differ_400(x)
            if not delta:
                return False
            smaller = pd.Series(x).rolling(window_slan_400_min, min_periods=1).apply(differ_400, raw=True).astype(dtype=bool).any()
            if smaller:
                return False
            return True

        def differ_slan_600(x):
            delta = differ_600(x)
            if not delta:
                return False
            smaller = pd.Series(x).rolling(window_slan_600_min, min_periods=1).apply(differ_600, raw=True).astype(dtype=bool).any() 
            if smaller:
                return False
            return True
        
        slan_400 = data_st["p"].rolling(window_slan_400, min_periods=1).apply(differ_slan_400, raw=True).astype(dtype=bool)#.compute()
        slan_600 = data_st["p"].rolling(window_slan_600, min_periods=1).apply(differ_slan_600, raw=True).astype(dtype=bool)#.compute()

        data_st = data_st.assign(conv_400 = conv_400)
        data_st = data_st.assign(conv_600 = conv_600)
        data_st = data_st.assign(slan_400 = slan_400)
        data_st = data_st.assign(slan_600 = slan_600)
        data_st = data_st.assign(traj_index = i)
        
        if i > 0:
            dd.to_parquet(data_st, store_path, append=True, ignore_divisions=True, compression="snappy")
        else:
            dd.to_parquet(data_st, store_path, compression="snappy")
            
        t2 = timer()
        print("Traj {} from {} - {} done after {} s".format(i, t, f, t2-t))

t_end = timer()
print("Finished after {} s".format(t_end-t_start))