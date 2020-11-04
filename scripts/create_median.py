import numpy as np
import os
import pandas as pd
import sys
from timeit import default_timer as timer
import xarray as xr

try:
    import Deriv_dask
    import latexify
    import loader
except:
    import scripts.Deriv_dask as Deriv_dask
    import scripts.latexify as latexify
    import scripts.loader as loader

def norm_time(df, norm_col, group, columns=None, flag=None):
    '''
    Return a view that consists only of entries that are flagged.
    Those are normed along norm_col such that every entry for every
    trajectory starts at norm_col==0. columns is a list of
    columns that the returned view shall have.
    if columns is None, take all columns. If flag is None, take all trajectories.
    '''
    if columns is None:
        df_flagged = df.copy()
    else:
        df_flagged = df[columns + [flag] + [norm_col] + [group]]

    def reducer(x, col):
        mini = x.loc[x[flag] == True][col].min()
        x[col] = x[col] - mini
        return x

    return df_flagged.groupby([group]).apply(reducer, norm_col)

def get_statistics_pandas(df, group, flag=None):
    '''
    Create a median, 25, and 75 percentile trajectory. If flag is set, Use only areas where flag is true.
    '''
    if flag is not None:
        df_tmp = df.loc[df[flag] == True]
    else:
        df_tmp = df.copy()
    drop_vars = ["WCB_flag", "dp2h", "slan_400", "slan_600", "conv_400", "conv_600", "type"]
    groupy = df_tmp.drop(drop_vars, axis=1).groupby(group)

    return groupy.median(), groupy.quantile(.25), groupy.quantile(.75)

def get_statistics_dask(df, group, flag=None):
    '''
    Create a median, 25, and 75 percentile trajectory. If flag is set, Use only areas where flag is true.
    '''
    if flag is not None:
        df_tmp = df.loc[df[flag] == True]
    else:
        df_tmp = df.copy()
    drop_vars = ["WCB_flag", "dp2h", "slan_400", "slan_600", "conv_400", "conv_600", "type"]
    groupy = df_tmp.drop_vars(drop_vars, errors="ignore").groupby(group)
    return groupy.median(keep_attr=True), groupy.quantile(0.25, keep_attrs=True), groupy.quantile(0.75, keep_attrs=True)

np.set_printoptions(threshold=sys.maxsize)

# 400 hPa and 600 hPa ascent
window_conv_400 = 1 * 3 * 60
window_conv_600 = 3 * 3 * 60
window_slan_400 = 35 * 6 * 3
window_slan_400_min = 15 * 6 * 3
window_slan_600 = 22  * 3 * 60
window_slan_600_min = 65 * 6 * 3

def find_runs(x):
    """Find runs of consecutive items in an array."""

    n = x.shape[0]

    # handle empty array
    if n == 0:
        return np.array([]), np.array([]), np.array([])

    else:
        # find run starts
        loc_run_start = np.empty(n, dtype=bool)
        loc_run_start[0] = True
        np.not_equal(x[:-1], x[1:], out=loc_run_start[1:])
        run_starts = np.nonzero(loc_run_start)[0]

        # find run values
        run_values = x[loc_run_start]

        # find run lengths
        run_lengths = np.diff(np.append(run_starts, n))

        return run_values, run_starts.astype(np.int64), run_lengths.astype(np.int64)

def differ(x, axis, hPa, debug=False):
    if debug:
        print("x")
        print(np.shape(x))
    window_size = len(x[0][0][0])
    ascent = np.argmax(x, axis=axis) < np.argmin(x, axis=3)
    amount = np.max(x, axis=axis) - np.min(x, axis=axis) >= hPa*100
    # ascent = np.nanargmax(x, axis=axis) < np.nanargmin(x, axis=3)
    # amount = np.nanmax(x, axis=axis) - np.nanmin(x, axis=axis) >= hPa*100
    both = np.logical_and(ascent, amount)

    # Calculate the differences within every window
    differences = np.diff(x, axis=3) # Ignore nan
    if debug:
        print("diffs")
        print(np.shape(differences))
        print("ascent")
        print(np.shape(ascent))
        print("amount")
        print(np.shape(amount))
        print("both")
        print(np.shape(both))
    # Get minimum length of window with value >= hPa
    min_lengths = np.full(np.shape(both), np.inf)
    counter = 0

    for ens in range(len(differences)):
        if not both[ens].any():
            continue
        for traj in range(len(differences[ens])):
            if not both[ens][traj].any():
                # No timestep in this trajectory satisfies
                # the constraint
                continue
            for timestep in range(len(differences[ens][traj])):
                if not both[ens][traj][timestep].any():
                    continue
                window = differences[ens][traj][timestep]
                start = 0
                end = 0
                curr_sum = 0
                min_len = np.inf
                    # 17926
                while(end < len(window)):
                    # Find a window where the ascend is done by pushing the end further
                    # Reset the start if the end is suddenly a strong descend
                    while( (curr_sum > -hPa*100 and end < len(window)) or (np.isnan(curr_sum)) ):
                        if (curr_sum >= 0 and window[end] < 0 or np.isnan(curr_sum) ):
                            start = end;
                            curr_sum = 0;

                        curr_sum += window[end]
                        end += 1
                    # Check, if a smaller window exists where the ascend is done by pushing the start
                    while(curr_sum <= -hPa*100 and start < len(window)):
                        if(end-start < min_len):
                            min_len = end-start
                        curr_sum -= window[start]
                        start += 1
                min_lengths[ens][traj][timestep] = min_len

    # Take the minimum overall and that's whenever True shall stand
    # Those are minimum window sizes for every trajectory
    return_bools = []
    min_len_traj_all = []
    for ens in range(len(x)):
        min_len_traj = np.nanmin(min_lengths[ens], axis=1)
        min_len_traj_all.append(min_len_traj)
        min_len_traj[min_len_traj == np.inf] = -1
        return_bools_tmp = np.full(np.shape(both[ens]), 0)#, dtype=bool)
        return_bools_tmp = np.transpose(return_bools_tmp)
        min_lengths_trans = min_lengths[ens].transpose()
        if debug:
            print("min_lengths_trans")
            print(np.shape(min_lengths_trans))
            print("min_len_traj")
            print(np.shape(min_len_traj))
            print("return_bools_tmp")
            print(np.shape(return_bools_tmp))
        for timestep in range(len(return_bools_tmp)):
            return_bools_tmp[timestep] = (min_lengths_trans[timestep] == min_len_traj)
        return_bools.append(np.transpose(return_bools_tmp))

    # Shift everything such that the beginning starts at the actual start and ends accordingly
    for ens in range(len(return_bools)):
        for traj in range(len(return_bools[ens])):
            if min_len_traj_all[ens][traj] == -1:
                continue
            min_len_traj = min_len_traj_all[ens]
            vals, start, length = find_runs(return_bools[ens][traj])

            for i in range(len(vals)):
                if vals[i] > 0:
                    set_start = int(start[i] - min_len_traj[traj])
                    set_end = set_start + min_len_traj[traj] + 1

                    if length[i] > min_len_traj[traj]:
                        set_end = set_start + length[i] + 1

                    return_bools[ens][traj][start[i]:length[i]+start[i]] = False
                    return_bools[ens][traj][set_start:int(set_end)] = True

    return return_bools

def differ_slan(x, axis, hPa, min_window):
    window_size = len(x[0][0][0])
    ascent = np.argmax(x, axis=axis) < np.argmin(x, axis=3)
    amount = np.max(x, axis=axis) - np.min(x, axis=axis) >= hPa*100
    both = np.logical_and(ascent, amount)

    # Calculate the differences within every window
    differences = np.diff(x, axis=3)
    # Get minimum length of window with value >= hPa
    min_lengths = np.full(np.shape(both), np.inf)


    for ens in range(len(differences)):
        if not both[ens].any():
            continue
        for traj in range(len(differences[ens])):
            if not both[ens][traj].any():
                # No timestep in this trajectory satisfies
                # the constraint
                continue
            for timestep in range(len(differences[ens][traj])):
                if not both[ens][traj][timestep].any():
                    continue
                window = differences[ens][traj][timestep]
                start = 0
                end = 0
                curr_sum = 0
                min_len = np.inf
                while(end < len(window)):
                    while( (curr_sum > -hPa*100 and end < len(window)) or (np.isnan(curr_sum)) ):
                        if (curr_sum >= 0 and window[end] < 0 or np.isnan(curr_sum) ):
                            start = end
                            curr_sum = 0
                        curr_sum += window[end]
                        end += 1
                    while(curr_sum <= -hPa*100 and start < len(window)):
                        if(end-start < min_len and end-start >= min_window):
                            min_len = end-start
                        curr_sum -= window[start]
                        start += 1
                min_lengths[ens][traj][timestep] = min_len
    # Take the minimum overall and that's whenever True shall stand
    # Those are minimum window sizes for every trajectory
    return_bools = []
    min_len_traj_all = []
    for ens in range(len(x)):
        min_len_traj = np.nanmin(min_lengths[ens], axis=1)
        min_len_traj_all.append(min_len_traj)
        min_len_traj[min_len_traj == np.inf] = -1
        return_bools_tmp = np.full(np.shape(both[ens]), 0)#, dtype=bool)
        return_bools_tmp = np.transpose(return_bools_tmp)
        min_lengths_trans = min_lengths[ens].transpose()
        for timestep in range(len(return_bools_tmp)):
            return_bools_tmp[timestep] = (min_lengths_trans[timestep] == min_len_traj)
        return_bools.append(np.transpose(return_bools_tmp))

    # Shift everything such that the beginning starts at the actual start and ends accordingly
    for ens in range(len(return_bools)):
        for traj in range(len(return_bools[ens])):
            if min_len_traj_all[ens][traj] == -1:
                continue
            vals, start, length = find_runs(return_bools[ens][traj])
            for i in range(len(vals)):
                if vals[i] > 0:
                    set_start = int(start[i] - min_len_traj[traj])
                    set_end = set_start + min_len_traj[traj] + 1

                    if length[i] > min_len_traj[traj]:
                        set_end = set_start + length[i] + 1

                    return_bools[ens][traj][start[i]:length[i]+start[i]] = False
                    return_bools[ens][traj][set_start:int(set_end)] = True
    return return_bools


def norm_time_xarray(df, norm_col, group, columns=None, flag=None):
    '''
    Return a view that consists only of entries that are flagged.
    Those are normed along norm_col such that every entry for every
    trajectory starts at norm_col==0. columns is a list of
    columns that the returned view shall have.
    if columns is None, take all columns. If flag is None, take all trajectories.
    '''
    if columns is None:
        df_flagged = df.copy()
    else:
        df_flagged = df[columns + [flag] + [norm_col] + [group]]

    def reducer(x, col):
        mini = x.where(x[flag])[col].min()
        x[col] = x[col] - mini
        return x
    return df_flagged.groupby(group).apply(reducer, **{"col": norm_col})


def get_statistics_xarray(df, group, flag=None):
    '''
    Create a median, 25, and 75 percentile trajectory. If flag is set, Use only areas where flag is true.
    '''
    if flag is not None:
        df_tmp = df.where(df[flag])
    else:
        df_tmp = df.copy()

    groupy = df_tmp.groupby(group)

    return groupy.median(), groupy.quantile(.25), groupy.quantile(.75)

def add_attrs(ds, ref_ds_path=None, attrs=None):
    if ref_ds_path is not None:
        ds_2 = xr.open_dataset(ref_ds_path, decode_times=False)
        duration = ds_2.attrs["duration_in_sec"]
        pollon = ds_2.attrs["pollon"]
        pollat = ds_2.attrs["pollat"]
        output_timestep_in_sec = ds_2.attrs["output_timestep_in_sec"]

        ds.attrs = {
            "duration_in_sec": duration,
            "pollon": pollon,
            "pollat": pollat,
            "output_timestep_in_sec": output_timestep_in_sec,
            "cloud_type": 2723}

        ds["time"].attrs = ds_2["time"].attrs

    if attrs is not None:
        ds.attrs = attrs["ds"]
        ds["time"].attrs = attrs["time"]

    ds["time_after_ascent"].attrs = {
        "standard_name": "time_after_ascent",
        "long_name": "time after rapid ascent started",
        "units": "seconds since start of convective/slantwise ascent"}
    ds["lon"].attrs = {
        "standard_name": "longitude",
        "long_name": "rotated longitude",
        "units": "degrees"}
    ds["lat"].attrs = {
        "standard_name": "latitude",
        "long_name": "rotated latitude",
        "units": "degrees"}
    ds["pressure"].attrs = {
        "standard_name": "air_pressure",
        "long_name": "pressure",
        "units": "Pa",
        "positive": "down",
        "axis": "z"}
    ds["z"].attrs = {
        "standard_name": "height",
        "long_name": "height above mean sea level",
        "auxiliary_data": "yes",
        "units": "m AMSL"}
    ds["T"].attrs = {
        "standard_name": "air_temperature",
        "long_name": "temperature",
        "auxiliary_data": "yes",
        "units": "K"}
    ds["S"].attrs = {
        "standard_name": "saturation",
        "long_name": "saturation",
        "auxiliary_data": "yes",
        "units": "percentage"}
    ds["conv_400"].attrs = {
        "standard_name": "convective_400hPa_ascent",
        "long_name": "convective 400hPa ascent",
        "auxiliary_data": "yes"}
    ds["conv_600"].attrs = {
        "standard_name": "convective_600hPa_ascent",
        "long_name": "convective 600hPa ascent",
        "auxiliary_data": "yes"}
    ds["slan_400"].attrs = {
        "standard_name": "slantwise_400hPa_ascent",
        "long_name": "slantwise 400hPa ascent",
        "auxiliary_data": "yes"}
    ds["slan_600"].attrs = {
        "standard_name": "slantwise_600hPa_ascent",
        "long_name": "slantwise 600hPa ascent",
        "auxiliary_data": "yes"}
    ds["w"].attrs = {
        "standard_name": "ascend_velocity",
        "long_name": "ascend velocity",
        "auxiliary_data": "yes",
        "units": "m s^-1"}

    ds["QV"].attrs = {
        "standard_name": "specific_humidity",
        "long_name": "specific humidity",
        "auxiliary_data": "yes",
        "units": "kg kg^-1"}
    ds["QC"].attrs = {
        "standard_name": "mass_fraction_of_cloud_liquid_water_in_air",
        "long_name": "specific cloud liquid water content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1"}
    ds["QR"].attrs = {
        "standard_name": "mass_fraction_of_rain_in_air",
        "long_name": "specific rain content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1"}
    ds["QS"].attrs = {
        "standard_name": "mass_fraction_of_snow_in_air",
        "long_name": "specific snow content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1"}
    ds["QI"].attrs = {
        "standard_name": "mass_fraction_of_cloud_ice_in_air",
        "long_name": "specific cloud ice content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1"}
    ds["QG"].attrs = {
        "standard_name": "mass_fraction_of_graupel_in_air",
        "long_name": "specific graupel content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1"}

    ds["QR_IN"].attrs = {
        "standard_name": "sedi_influx_of_rain",
        "long_name": "sedimentation (from above) of rain droplet mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1"}
    ds["QS_IN"].attrs = {
        "standard_name": "sedi_influx_of_snow",
        "long_name": "sedimentation (from above) of snow crystal mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1"}
    ds["QI_IN"].attrs = {
        "standard_name": "sedi_influx_of_cloud_ice",
        "long_name": "sedimentation (from above) of ice crystal mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1"}
    ds["QG_IN"].attrs = {
        "standard_name": "sedi_influx_of_graupel",
        "long_name": "sedimentation (from above) of graupel mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1"}

    ds["QR_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_rain",
        "long_name": "sedimentation of rain droplet mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1"}
    ds["QS_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_snow",
        "long_name": "sedimentation of snow crystal mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1"}
    ds["QI_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_cloud_ice",
        "long_name": "sedimentation of ice crystal mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1"}
    ds["QG_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_graupel",
        "long_name": "sedimentation of graupel mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1"}

    ds["NCCLOUD"].attrs = {
        "standard_name": "specif_number_of_cloud_droplets_in_air",
        "long_name": "specific cloud droplet number",
        "auxiliary_data": "yes",
        "units": "kg^-1"}
    ds["NCRAIN"].attrs = {
        "standard_name": "specif_number_of_rain_drops_in_air",
        "long_name": "specific rain drop number",
        "auxiliary_data": "yes",
        "units": "kg^-1"}
    ds["NCSNOW"].attrs = {
        "standard_name": "specif_number_of_snow_flakes_in_air",
        "long_name": "specific snow flake number",
        "auxiliary_data": "yes",
        "units": "kg^-1"}
    ds["NCICE"].attrs = {
        "standard_name": "specif_number_of_cloud_ice_in_air",
        "long_name": "specific cloud ice number",
        "auxiliary_data": "yes",
        "units": "kg^-1"}
    ds["NCGRAUPEL"].attrs = {
        "standard_name": "specif_number_of_graupel_in_air",
        "long_name": "specific graupel number",
        "auxiliary_data": "yes",
        "units": "kg^-1"}

    ds["NR_IN"].attrs = {
        "standard_name": "sedi_influx_of_rain_number",
        "long_name": "sedimentation (from above) of specific rain drop number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}
    ds["NS_IN"].attrs = {
        "standard_name": "sedi_influx_of_snow_number",
        "long_name": "sedimentation (from above) of specific snow flake number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}
    ds["NI_IN"].attrs = {
        "standard_name": "sedi_influx_of_ics_number",
        "long_name": "sedimentation (from above) of specific cloud ice number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}
    ds["NG_IN"].attrs = {
        "standard_name": "sedi_influx_of_graupel_number",
        "long_name": "sedimentation (from above) of specific graupel number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}

    ds["NR_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_rain_number",
        "long_name": "sedimentation of rain droplet number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}
    ds["NS_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_snow_number",
        "long_name": "sedimentation of snow crystal number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}
    ds["NI_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_ice_number",
        "long_name": "sedimentation of ice crystal number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}
    ds["NG_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_graupel_number",
        "long_name": "sedimentation of graupel number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}

    ds["Q_TURBULENCE"].attrs = {
        "standard_name": "turbulence_flux",
        "long_name": "flux from turbulence",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1"}
    ds["type"].attrs = {
        "standard_name": "trajectory_type",
        "long_name": "trajectory type",
        "auxiliary_data": "yes"}

    return ds

def do_the_stuff_more(store_path, version="no exclusions", fls=["conv_400", "conv_600", "slan_400", "slan_600"], file_list=None):
    """
    Versions
    --------
    no exclusions: Get all trajectories with the corresponding flag
    excl other: Get all trajectories with the corresponding flag without the other type (ie flag is conv_400, excl slan_600)
    excl same: Get all trajectories with the corresponding flag without the same type (ie flag is conv_400, excl conv_600)
    excl all: Get all trajectories with the corresponding flag without all others (ie flag is conv_400, excl conv_600 and slan_600)
    conv_X and slan_X are mutually exclusive by definition already!
    """
    store_path = store_path + version.replace(" ", "_") + "_"
    n = 0
    # datasets = []
    ds = None
    ref_ds_path = None
    # iteri = 0
    for fl in fls:
        print(f"################### Running for {fl} ###################")
        t_c = timer()
        for f in file_list:
            if fl not in f:
                continue
            print(f"Loading {f}")
            if ref_ds_path is None:
                ref_ds_path = f
            # iteri += 1
            # if iteri >= 5:
            #     continue
            # ds_tmp = xr.open_dataset(f, decode_times=False)
            ds_tmp = xr.open_dataset(f, decode_times=False).to_dataframe().reset_index()
            ds_tmp = ds_tmp[(ds_tmp.time_after_ascent >= -10000)]

            if version == "excl other":
                fl_no = "slan_400"
                if fl == "conv_400":
                    fl_no = "slan_600"
                elif fl == "slan_400":
                    fl_no = "conv_600"
                elif fl == "slan_600":
                    fl_no = "conv_400"
                ids = np.unique(ds_tmp.where(ds_tmp[fl] == True)["trajectory"])
                non_ids = np.unique(ds_tmp.where(ds_tmp[fl_no] == True)["trajectory"])
                ids = np.setdiff1d(ids, non_ids)
                ds_tmp = ds_tmp.where(ds_tmp["trajectory"].isin(ids))
            elif version == "excl same":
                fl_no = "conv_400"
                if fl == "conv_400":
                    fl_no = "conv_600"
                elif fl == "slan_400":
                    fl_no = "slan_600"
                elif fl == "slan_600":
                    fl_no = "slan_400"
                ids = np.unique(ds_tmp.where(ds_tmp[fl] == True)["trajectory"])
                non_ids = np.unique(ds_tmp.where(ds_tmp[fl_no] == True)["trajectory"])
                ids = np.setdiff1d(ids, non_ids)
                ds_tmp = ds_tmp.where(ds_tmp["trajectory"].isin(ids))
            elif version == "excl all":
                for fl_no in ["conv_400", "conv_600", "slan_400", "slan_600"]:
                    if fl_no == fl:
                        continue
                    ids = np.unique(ds_tmp.where(ds_tmp[fl] == True)["trajectory"])
                    non_ids = np.unique(ds_tmp.where(ds_tmp[fl_no] == True)["trajectory"])
                    ids = np.setdiff1d(ids, non_ids)
                    ds_tmp = ds_tmp.where(ds_tmp["trajectory"].isin(ids))
            elif version != "no exclusions":
                print(f"version {version} unknown!")
                return

            if ds is not None:
                ds = ds.append(ds_tmp)
            else:
                ds = ds_tmp
            # datasets.append(ds_tmp)

        # ds = xr.merge(datasets, join="outer")
        t_c2 = timer()
        print("loading done in {} s".format(t_c2-t_c), flush=True)
        ds = ds.dropna()
        print(ds.describe())
        t_c = timer()
        # medi, quan25, quan75 = get_statistics_dask(ds, group="time_after_ascent")
        medi, quan25, quan75 = get_statistics_pandas(ds, group="time_after_ascent")
        medi["trajectory"] = n
        n += 1
        quan25["trajectory"] = n
        n += 1
        quan75["trajectory"] = n
        n += 1
        t_c2 = timer()
        print("statistics done in {} s".format(t_c2-t_c), flush=True)

        t_c = timer()
        # Set flags for slantwise or convective parts
        for flag in ["conv_600", "conv_400", "slan_400", "slan_600"]:
            medi[flag] = False
            quan25[flag] = False
            quan75[flag] = False

        # The time needs to be adjusted
        def adjust(this_df):
            this_df = this_df.reset_index()
            this_df["ensemble"] = n%3
            start_time = this_df["time"][0]
            end_time = start_time + 20*len(this_df.index)
            this_df["time"] = np.arange(start_time, end_time, 20)
            return this_df
        medi = adjust(medi)
        quan25 = adjust(quan25)
        quan75 = adjust(quan75)

        medi = xr.Dataset.from_dataframe(medi.set_index(["ensemble", "trajectory", "time"]))
        quan25 = xr.Dataset.from_dataframe(quan25.set_index(["ensemble", "trajectory", "time"]).dropna())
        quan75 = xr.Dataset.from_dataframe(quan75.set_index(["ensemble", "trajectory", "time"]).dropna())

        # Set a new column name to the corresponding trajectories
        type_name = ""
        if fl == "conv_600":
            type_name = "Convective 600hPa"
        elif fl == "conv_400":
            type_name = "Convective 400hPa"
        elif fl == "slan_400":
            type_name = "Slantwise 400hPa"
        else:
            type_name = "Slantwise 600hPa"

        other_name = ""
        if version == "vanilla":
            other_name = "no exclusions"
        medi["type"] = type_name + " 50. Quantile"
        quan25["type"] = type_name + " 25. Quantile"
        quan75["type"] = type_name + " 75. Quantile"

        if fl == "conv_600":
            t_c = timer()
            conv_600 = medi["pressure"].rolling(dim={"time": window_conv_600}, min_periods=1).reduce(
                differ, **{"hPa": 600}).fillna(False).astype(dtype=bool)
            medi = medi.assign(conv_600=conv_600)
            t_c2 = timer()
            print("Got conv_600 in {} s".format(t_c2-t_c), flush=True)

            t_c = timer()
            conv_600 = quan25["pressure"].rolling(dim={"time": window_conv_600}, min_periods=1).reduce(
                differ, **{"hPa": 600, "debug": False}).fillna(False).astype(dtype=bool)
            quan25 = quan25.assign(conv_600=conv_600)
            t_c2 = timer()
            print("Got conv_600 in {} s".format(t_c2-t_c), flush=True)

            t_c = timer()
            conv_600 = quan75["pressure"].rolling(dim={"time": window_conv_600}, min_periods=1).reduce(
                differ, **{"hPa": 600}).fillna(False).astype(dtype=bool)
            quan75 = quan75.assign(conv_600=conv_600)
            t_c2 = timer()
            print("Got conv_600 in {} s".format(t_c2-t_c), flush=True)
        elif fl == "conv_400":
            t_c = timer()
            conv_400 = medi["pressure"].rolling(dim={"time": window_conv_400}, min_periods=1).reduce(
                differ, **{"hPa": 400}).fillna(False).astype(dtype=bool)
            medi = medi.assign(conv_400=conv_400)
            t_c2 = timer()
            print("Got conv_400 in {} s".format(t_c2-t_c), flush=True)

            t_c = timer()
            conv_400 = quan25["pressure"].rolling(dim={"time": window_conv_400}, min_periods=1).reduce(
                differ, **{"hPa": 400}).fillna(False).astype(dtype=bool)
            quan25 = quan25.assign(conv_400=conv_400)
            t_c2 = timer()
            print("Got conv_400 in {} s".format(t_c2-t_c), flush=True)

            t_c = timer()
            conv_400 = quan75["pressure"].rolling(dim={"time": window_conv_400}, min_periods=1).reduce(
                differ, **{"hPa": 400}).fillna(False).astype(dtype=bool)
            quan75 = quan75.assign(conv_400=conv_400)
            t_c2 = timer()
            print("Got conv_400 in {} s".format(t_c2-t_c), flush=True)
        elif fl == "slan_600":
            t_c = timer()
            slan_600 = medi["pressure"].rolling(dim={"time": window_slan_600}, min_periods=1).reduce(
                differ_slan, **{"hPa": 600, "min_window": window_slan_600_min}).fillna(False).astype(dtype=bool)
            medi = medi.assign(slan_600=slan_600)
            t_c2 = timer()
            print("Got slan_600 in {} s".format(t_c2-t_c), flush=True)

            t_c = timer()
            slan_600 = quan25["pressure"].rolling(dim={"time": window_slan_600}, min_periods=1).reduce(
                differ_slan, **{"hPa": 600, "min_window": window_slan_600_min}).fillna(False).astype(dtype=bool)
            quan25 = quan25.assign(slan_600=slan_600)
            t_c2 = timer()
            print("Got slan_600 in {} s".format(t_c2-t_c), flush=True)

            t_c = timer()
            slan_600 = quan75["pressure"].rolling(dim={"time": window_slan_600}, min_periods=1).reduce(
                differ_slan, **{"hPa": 600, "min_window": window_slan_600_min}).fillna(False).astype(dtype=bool)
            quan75 = quan75.assign(slan_600=slan_600)
            t_c2 = timer()
            print("Got slan_600 in {} s".format(t_c2-t_c), flush=True)
        elif fl == "slan_400":
            t_c = timer()
            slan_400 = medi["pressure"].rolling(dim={"time": window_slan_400}, min_periods=1).reduce(
                differ_slan, **{"hPa": 400, "min_window": window_slan_400_min}).fillna(False).astype(dtype=bool)
            medi = medi.assign(slan_400=slan_400)
            t_c2 = timer()
            print("Got slan_400 in {} s".format(t_c2-t_c), flush=True)

            t_c = timer()
            slan_400 = quan25["pressure"].rolling(dim={"time": window_slan_400}, min_periods=1).reduce(
                differ_slan, **{"hPa": 400, "min_window": window_slan_400_min}).fillna(False).astype(dtype=bool)
            quan25 = quan25.assign(slan_400=slan_400)
            t_c2 = timer()
            print("Got slan_400 in {} s".format(t_c2-t_c), flush=True)

            t_c = timer()
            slan_400 = quan75["pressure"].rolling(dim={"time": window_slan_400}, min_periods=1).reduce(
                differ_slan, **{"hPa": 400, "min_window": window_slan_400_min}).fillna(False).astype(dtype=bool)
            quan75 = quan75.assign(slan_400=slan_400)
            t_c2 = timer()
            print("Got slan_400 in {} s".format(t_c2-t_c), flush=True)
        print("Rolling done")

        t_c = timer()
        # Set time_after_ascent to zero where flag first occurs
        def adjust_ascent_time(this_ds):
            start_idx = np.where(~np.isnan(this_ds.where(this_ds[fl] == True)["time_after_ascent"]))[2][0]
            start_time = this_ds.isel(time=start_idx)["time_after_ascent"].values
            this_ds["time_after_ascent"] -= start_time
            return this_ds

        medi = adjust_ascent_time(medi)
        quan25 = adjust_ascent_time(quan25)
        quan75 = adjust_ascent_time(quan75)

        medi = add_attrs(medi, ref_ds_path)
        quan25 = add_attrs(quan25, ref_ds_path)
        quan75 = add_attrs(quan75, ref_ds_path)
        t_c2 = timer()
        print("Set time after ascent in {} s".format(t_c2-t_c), flush=True)

        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in medi.data_vars}
        medi.to_netcdf(store_path + fl + "_median.nc_wcb",
            encoding=encoding, compute=True, engine="netcdf4", format="NETCDF4", mode="w")
        quan25.to_netcdf(store_path + fl + "_quan25.nc_wcb",
            encoding=encoding, compute=True, engine="netcdf4", format="NETCDF4", mode="w")
        quan75.to_netcdf(store_path + fl + "_quan75.nc_wcb",
            encoding=encoding, compute=True, engine="netcdf4", format="NETCDF4", mode="w")
        print("storing done")

def add_flags():
    file_list = []
    for f in os.listdir("/data/project/wcb/netcdf/traj_stats/"):
        file_list.append(os.path.join("/data/project/wcb/netcdf/traj_stats/", f))

    file_list = np.sort(np.asarray(file_list))
    for f in file_list:
        with xr.open_dataset(f) as ds:
            if "conv_400" in f:
                conv_400 = (ds["P"].rolling(dim={"time": window_conv_400}, min_periods=1).reduce(
                    differ, **{"hPa": 400}).fillna(False).astype(dtype=bool))
                ds = ds.assign(conv_400=conv_400)
            elif "conv_600" in f:
                conv_600 = (ds["P"].rolling(dim={"time": window_conv_600}, min_periods=1).reduce(
                    differ, **{"hPa": 600}).fillna(False).astype(dtype=bool))
                ds = ds.assign(conv_600=conv_600)
            elif "slan_400" in f:
                slan_400 = (ds["P"].rolling(dim={"time": window_slan_400}, min_periods=1).reduce(
                    differ_slan, **{"hPa": 400, "min_window": window_slan_400_min}).fillna(False).astype(dtype=bool))
                ds = ds.assign(slan_400=slan_400)
            elif "slan_600" in f:
                slan_600 = (ds["P"].rolling(dim={"time": window_slan_600}, min_periods=1).reduce(
                    differ_slan, **{"hPa": 600, "min_window": window_slan_600_min}).fillna(False).astype(dtype=bool))
                ds = ds.assign(slan_600=slan_600)
            ds.to_netcdf(f + "2")

def do_the_stuff_more_slan_600():
    file_list = []
    for f in os.listdir(netcdf_path):
        if os.path.isfile(os.path.join(netcdf_path, f)):
            file_list.append(os.path.join(netcdf_path, f))

    file_list = np.sort(np.asarray(file_list))

    store_path_tmp = "/data/project/wcb/netcdf/tmp/"
    n = 0
    fl = "slan_600"
    non_fl = ["slan_400", "conv_400"]
    # First we norm the timesteps for every other file
    for i in range(0, len(file_list), 3):
        t_c = timer()
        ds = None
        for f in file_list[i:i+3]:
            ds_tmp = xr.open_dataset(f).to_dataframe().reset_index()
            ids = ds_tmp.loc[ds_tmp[fl] == True]["id"]
            ds_tmp = ds_tmp.loc[ds_tmp["id"].isin(ids)]
            if ds is None:
                ds = ds_tmp
            else:
                ds = ds.append(ds_tmp)

        t_c2 = timer()
        print("loading done in {} s".format(t_c2-t_c), flush=True)
        # Make id and time columns instead of MultiIndex
        t_c = timer()
        normed = norm_time(ds, "time", "id", None, fl)
        t_c2 = timer()
        print("norming done in {} s".format(t_c2-t_c), flush=True)
        for flag in ["conv_600", "conv_400", "slan_400", "slan_600"]:
            normed[flag] = False
        t_c = timer()
        xr.Dataset.from_dataframe(normed.reset_index().set_index(["time", "id"])).to_netcdf(store_path_tmp + fl + "_normed_{}".format(i//3))
        t_c2 = timer()
        print("storing done in {} s".format(t_c2-t_c), flush=True)

    ds = None
    file_list = []
    for f in os.listdir(store_path_tmp):
        if os.path.isfile(os.path.join(store_path_tmp, f)):
            file_list.append(os.path.join(store_path_tmp, f))

    file_list = np.sort(np.asarray(file_list))

    # Load n rows of the trajectories and do statistics on the normed variants
    finished = [False for _ in file_list]
    n = 200000
    min_time = -385400.0
    for i in range(n):
        ds = None
        for j, f in enumerate(file_list):
            if finished[j]:
                continue
            ds_tmp = xr.open_dataset(f).to_dataframe().reset_index()
            # Load only certain timesteps
            ds_tmp = ds_tmp.loc[(ds_tmp.time > (min_time + n*i)) & (ds_tmp.time < (min_time + n*(i+1)))]

            if ds_tmp.empty:
                finished[j] = True
                print("Finished loading for {}".format(f))
                continue
            if ds is None:
                ds = ds_tmp
            else:
                ds = ds.append(ds_tmp)
        if np.sum(finished) == len(file_list):
            break

        t_c = timer()
        medi, quan25, quan75 = get_statistics_pandas(ds, group="time")
        ds = None
        medi["id"] = 10
        quan25["id"] = 11
        quan75["id"] = 12
        t_c2 = timer()
        print("statistics done in {} s".format(t_c2-t_c), flush=True)
        t_c = timer()
        # Set flags for slantwise or convective parts
        for flag in ["conv_600", "conv_400", "slan_400", "slan_600"]:
            medi[flag] = False
            quan25[flag] = False
            quan75[flag] = False

        medi = xr.Dataset.from_dataframe(medi.reset_index().set_index(["time", "id"]).dropna())
        quan25 = xr.Dataset.from_dataframe(quan25.reset_index().set_index(["time", "id"]).dropna())
        quan75 = xr.Dataset.from_dataframe(quan75.reset_index().set_index(["time", "id"]).dropna())

        type_name = "Slantwise 600hPa"
        medi["type"] = type_name + " 50. Quantile"
        quan25["type"] = type_name + " 25. Quantile"
        quan75["type"] = type_name + " 75. Quantile"
        t_c2 = timer()
        print("Resetting done in {} s".format(t_c2-t_c), flush=True)


        t_c = timer()
        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in ds.data_vars}
        medi.to_netcdf(store_path_tmp + "medi/" + fl + "_median_{}.nc_wcb".format(i),
            encoding=encoding, compute=True, engine="netcdf4", format="NETCDF4", mode="w")
        quan25.to_netcdf(store_path_tmp + "quan25/" + fl + "_quan25_{}.nc_wcb".format(i),
            encoding=encoding, compute=True, engine="netcdf4", format="NETCDF4", mode="w")
        quan75.to_netcdf(store_path_tmp + "quan75/" + fl + "_quan75_{}.nc_wcb".format(i),
            encoding=encoding, compute=True, engine="netcdf4", format="NETCDF4", mode="w")
        t_c2 = timer()
        print("storing done in {} s".format(t_c2 - t_c), flush=True)


def merge_stuff(input_path, output,
    ensembles=["conv_400", "conv_600", "slan_400", "slan_600"], dropna=False):
    """
    Merge statistical trajectories to a single file.
    Each category given in ensembles will be merged to a single ensemble.
    The trajectories will be enumerated accordingly.

    Parmaters
    ---------
    input_path: string
        Path to trajectory files.
    output: string
        Path + name where new file shall be stored.
    ensembles: list of string
        Categories to merge into a single ensembles (must be present in filename).

    Returns
    -------
    xarray.DataSet
        The merged dataset.
    """
    file_list = []
    for f in os.listdir(input_path):
        file_list.append(os.path.join(input_path, f))
    file_list = np.sort(np.asarray(file_list))
    ds_list = []
    for ens in ensembles:
        ds_ens = []
        for f in file_list:
            if ens in f:
                ds_tmp = xr.open_dataset(f, decode_times=False)
                # Unfortunately, we need to make sure the time coordinates
                # are somewhat similar ie modulo 20
                # Since we are merging "representatives" the exact time
                # should not matter
                ds_tmp["time"] = ds_tmp["time"] - ds_tmp["time"]%20
                ds_ens.append(ds_tmp)
        ds_list.append(xr.combine_nested(ds_ens, concat_dim="trajectory"))
    merged_ds = xr.combine_nested(ds_list, concat_dim="ensemble")
    min_t = merged_ds["time"].min()
    max_t = merged_ds["time"].max()

    attrs = {
        "time": {
            "standard_name": "time",
            "long_name": "time",
            "units": "seconds since 2016-09-20 00:00:00",
            "trajectory_starttime": "2016-09-20 00:00:00",
            "forecast_inittime": "2016-09-20 00:00:00"
        },
        "ds": {
            "duration_in_sec": (max_t-min_t).item(),
            "pollon": 160.,
            "pollat": 51.,
            "output_timestep_in_sec": 20,
            "cloud_type": 2723
        }
    }
    merged_ds = add_attrs(merged_ds, attrs=attrs)
    if dropna:
        merged_ds = merged_ds.dropna(dim="time", how="any")
        print(merged_ds)

    comp = dict(zlib=True, complevel=9)
    encoding = {var: comp for var in merged_ds.data_vars}

    merged_ds.to_netcdf(output, compute=True, engine="netcdf4",
        encoding=encoding, format="NETCDF4", mode="w")
    return merged_ds


if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1] # /data/project/wcb/netcdf/vladiana_met/
    else:
        path = "../data/sim_processed/conv_400_0_t000000_p001_mult_outSat_sbShape_sbConv/"
    if len(sys.argv) > 2:
        store_path = sys.argv[2] # /data/project/wcb/netcdf/vladiana_met_stats/
    else:
        store_path = "../data/median/"
    if len(sys.argv) > 3:
        flags = [sys.argv[3]]
    else:
        flags = ["slan_400", "slan_600"]

    file_list = []
    for f in os.listdir(path):
        file_list.append(os.path.join(path, f))

    file_list = np.sort(np.asarray(file_list))
    do_the_stuff_more(store_path, "no exclusions", flags, file_list)


