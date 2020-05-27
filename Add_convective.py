# Script to add convective flag to trajectories
# WCB trajectories are convective if
# 400 hPa had been ascended within 1 hour
# or 600 hPa had been ascended within 3 hours
#
# WCB trajectories are slantwise if
# 400 hPa ascent times are between 1.5 and 3.5 hours
# or 600 hPa ascent times are between 6.5 and 22 hours

import dask.dataframe as dd
import numpy as np
import os
import sys


direc_path = "/lustre/project/m2_zdvresearch/mahieron/parquet_complete/wcb_traj_ratio_concat"

data = dd.read_parquet(direc_path)
# Get the different original sources
sources = data["netcdf"].unique().compute()
for s in sources:
    data_s = data[data["netcdf"] == s]
#     trajectories = np.arange(3000) # There are definitely less than that
    trajectories = data_s.trajectory.unique().compute()
    for t in trajectories:
        data_st = data_s[data_s.trajectory == t]
        
        # Each entry is 0.01 s
        # We are interested in values of 1.5-3.5 hours
        # and values between 6.5 and 22 hours
        # Calculate the minimum and maximum pressure in those windows
        # Idea use rolling window of maximum time to see, if pressure had been reached
        # ToDo: How to ensure minimum time?
        window_conv_400 = 1.0 * 60 * 60 * 100
        window_conv_600 = 3.0 * 60 * 60 * 100
        window_slan_400 = 3.5 * 60 * 60 * 100
        window_slan_400_min = 1.5 * 60 * 60 * 100
        window_slan_600 = 22  * 60 * 60 * 100
        window_slan_600_min = 6.5 * 60 * 60 * 100
       
        differ_400 = lambda x: ((x.max()-x.min()) >= 40000)
        differ_600 = lambda x: ((x.max()-x.min()) >= 60000)
            
        conv_400 = data_st.rolling(window_conv_400, min_periods=1).apply(differ_400).astype(dtype=bool)
        conv_600 = data_st.rolling(window_conv_600, min_periods=1).apply(differ_600).astype(dtype=bool)
        
        def differ_slan_400(x):
            delta = differ_400(x)
            if not delta:
                return False
            smaller = x.rolling(window_slan_400_min, min_periods=1).apply(differ_400).astype(dtype=bool).any()
            if smaller:
                return False
            return True
        
        def differ_slan_600(x):
            delta = differ_600(x)
            if not delta:
                return False
            smaller = x.rolling(window_slan_600_min, min_periods=1).apply(differ_600).astype(dtype=bool).any()
            if smaller:
                return False
            return True
        slan_400 = data_st.rolling(window_slan_400, min_periods=1).apply(differ_slan_400).astype(dtype=bool)
        slan_600 = data_st.rolling(window_slan_600, min_periods=1).apply(differ_slan_600).astype(dtype=bool)
        
        
            
        # TODO: Add plots of pressure differences such that I can see if trajecotries always ascend.