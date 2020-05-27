#!/bin/env python

#SBATCH -p parallel
#SBATCH -A m2_zdvresearch
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --job-name=Concat_Flags
#SBATCH --output=logs/Concat_Flags_%A.out
#SBATCH --error=logs/Concat_Flags_%A.err
#SBATCH --mem=100G

import sys
import os
sys.path.append(os.getcwd())
import multiprocessing

import numpy as np
# from progressbar import progressbar as pb
import xarray as xr
import dask.dataframe as dd
import pandas as pd
from timeit import default_timer as timer


direc_path = "/lustre/project/m2_zdvresearch/mahieron/parquet_ratio_complete"
store_path = "/lustre/project/m2_zdvresearch/mahieron/parquet_ascent"

file_list = [os.path.join(direc_path, f) + "/*.parquet" for f in os.listdir(direc_path)]

written = False
get_new = True

window_conv_400 = 1 * 60 * 60
window_conv_600 = 3 * 60 * 60
window_slan_400 = 35 * 6 * 60 
window_slan_400_min = 15 * 6 * 60 
window_slan_600 = 22  * 60 * 60 
window_slan_600_min = 65 * 6 * 60 

differ_400 = lambda x: ((x.max()-x.min()) >= 4000)
differ_600 = lambda x: ((x.max()-x.min()) >= 6000)

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

# for i, f in pb(enumerate(file_list)):
for i, f in enumerate(file_list):
    if get_new:
        t = timer()
        get_new = False
        df = dd.read_parquet(f)
        name = f.split("_ratio.parq")
        name = name[0]
        name = name.split("_t")
        name = name[-1] + "t"
        df["netcdf"] = name
        df["traj_index"] = i
        
        conv_400 = df["p"].rolling(window_conv_400, min_periods=1).apply(differ_400, raw=True).astype(dtype=bool)
        conv_600 = df["p"].rolling(window_conv_600, min_periods=1).apply(differ_600, raw=True).astype(dtype=bool)
        slan_400 = df["p"].rolling(window_slan_400, min_periods=1).apply(differ_slan_400, raw=True).astype(dtype=bool)
        slan_600 = df["p"].rolling(window_slan_600, min_periods=1).apply(differ_slan_600, raw=True).astype(dtype=bool)
        df = df.assign(conv_400 = conv_400)
        df = df.assign(conv_600 = conv_600)
        df = df.assign(slan_400 = slan_400)
        df = df.assign(slan_600 = slan_600)
        
        
    else:
        small_set = dd.read_parquet(f)
        name = f.split("_ratio.parq")
        name = name[0]
        name = name.split("_t")
        name = name[-1] + "t"
        small_set["netcdf"] = name
        small_set["traj_index"] = i
        
        conv_400 = small_set["p"].rolling(window_conv_400, min_periods=1).apply(differ_400, raw=True).astype(dtype=bool)
        conv_600 = small_set["p"].rolling(window_conv_600, min_periods=1).apply(differ_600, raw=True).astype(dtype=bool)
        slan_400 = small_set["p"].rolling(window_slan_400, min_periods=1).apply(differ_slan_400, raw=True).astype(dtype=bool)
        slan_600 = small_set["p"].rolling(window_slan_600, min_periods=1).apply(differ_slan_600, raw=True).astype(dtype=bool)
        small_set = small_set.assign(conv_400 = conv_400)
        small_set = small_set.assign(conv_600 = conv_600)
        small_set = small_set.assign(slan_400 = slan_400)
        small_set = small_set.assign(slan_600 = slan_600)
        
        df = dd.concat([df, small_set])

    if (i%10==0 and i>0) or i == len(file_list)-1:
        if written:
            dd.to_parquet(df, store_path, append=True,  compression="snappy") # ignore_divisions=True,
            t2 = timer()
            print("Writing step {} done in {} s".format(i, t2-t))
            quit()
        else:
            dd.to_parquet(df, store_path, compression="snappy")
            written = True
        t2 = timer()
        print("Writing step {} done in {} s".format(i, t2-t))
        get_new = True
    
    
