import os

from scripts.loader import load_mult_derivates_directory
from scripts.loader import load_nc, rotate_df, norm_deriv, ratio_deriv
import scripts.loader as loader
import scripts.plot_mapped as plotter
import scripts.latexify as latexify
import numpy as np
from pylab import rcParams
import xarray as xr
import dask.dataframe as pd
import sys


direc_path = ".."
store_path = direc_path + "/wcb_traj_ratio_concat"
paths = []
# 54 takes about 575 GBytes of RAM


load_path = direc_path + "/traj_ratio.parquet/*.parquet"
paths.append(load_path)
df = pd.read_parquet(paths)
# print("Types of traj parquet:")

# df2 = pd.read_parquet(store_path)
# print("Types of concat parquet:")
# for stuff, stuff2 in zip(df.dtypes, df2.dtypes):
#     print("{}  --  {}".format(stuff, stuff2))

# print("n_cols: {} vs {}".format(len(df.dtypes), len(df2.dtypes)))

pd.to_parquet(df, store_path, append=True, ignore_divisions=True)



# direc_path = "/data/project/wcb/wcb_traj_flag_deriv_mult_min"
# store_path = direc_path + "/wcb_traj_ratio_concat"
# paths = []
# for t in [0, 1, 2]:
#     load_path = direc_path + "/wcb_traj{}_ratio.parquet/*.parquet".format(t)
#     paths.append(load_path)
# df = pd.read_parquet(paths)

# pd.to_parquet(df, store_path)

