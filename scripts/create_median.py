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
    import convert_to_met3d.find_runs as find_runs
    import convert_to_met3d.differ as differ
    import convert_to_met3d.differ_slan as differ_slan

    import convert_to_met3d.window_conv_400 as window_conv_400
    import convert_to_met3d.window_conv_600 as window_conv_600
    import convert_to_met3d.window_slan_400 as window_slan_400
    import convert_to_met3d.window_slan_400_min as window_slan_400_min
    import convert_to_met3d.window_slan_600 as window_slan_600
    import convert_to_met3d.window_slan_600_min as window_slan_600_min
except:
    import scripts.Deriv_dask as Deriv_dask
    import scripts.latexify as latexify
    import scripts.loader as loader
    import scripts.convert_to_met3d.find_runs as find_runs
    import scripts.convert_to_met3d.differ as differ
    import scripts.convert_to_met3d.differ_slan as differ_slan

    import scripts.convert_to_met3d.window_conv_400 as window_conv_400
    import scripts.convert_to_met3d.window_conv_600 as window_conv_600
    import scripts.convert_to_met3d.window_slan_400 as window_slan_400
    import scripts.convert_to_met3d.window_slan_400_min as window_slan_400_min
    import scripts.convert_to_met3d.window_slan_600 as window_slan_600
    import scripts.convert_to_met3d.window_slan_600_min as window_slan_600_min

pandas.options.mode.chained_assignment = None
np.set_printoptions(threshold=sys.maxsize)


def norm_time(df, norm_col, group, columns=None, flag=None):
    """
    Return a view that consists only of entries that are flagged.
    Those are normed along norm_col such that every entry for every
    trajectory starts at norm_col==0. columns is a list of
    columns that the returned view shall have.
    if columns is None, take all columns. If flag is None, take all trajectories.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns norm_col, group and flag at minimum.
    norm_col : string
        Column where the minimum value where flag is true is substracted from.
        Usually it is "time"
    group : string
        Column for grouping operation, i.e. "trajectory".
    columns : list of string
        List of columns that the returned view shall have. If None is given,
        return all columns.
    flag : string
        Column with bools, i.e. start of an ascend.

    Returns
    -------
    pandas.Dataarray shifted norm_col.
    """
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
    """
    Create a median, 25, and 75 percentile trajectory.
    If flag is set, Use only areas where flag is true.
    Drops flag and type columns since those might not be valid anymore and
    need to be recalculated.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns "WCB_flag", "dp2h", "slan_400", "slan_600",
        "conv_400", "conv_600", "type" and group.
    group : string
        Column to group by, e.g. "time_after_ascent" to calculate the percentiles
        along.
    flag : string
        Column with bools.
    Returns
    -------
    Median, 25, 75 percentile trajectory
    """
    if flag is not None:
        df_tmp = df.loc[df[flag] == True]
    else:
        df_tmp = df.copy()
    drop_vars = [
        "WCB_flag",
        "dp2h",
        "slan_400",
        "slan_600",
        "conv_400",
        "conv_600",
        "type",
    ]
    groupy = df_tmp.drop(drop_vars, axis=1).groupby(group)

    return groupy.median(), groupy.quantile(0.25), groupy.quantile(0.75)


def get_statistics_dask(df, group, flag=None):
    """
    Create a median, 25, and 75 percentile trajectory.
    If flag is set, Use only areas where flag is true.
    Drops flag and type columns since those might not be valid anymore and
    need to be recalculated.

    Parameters
    ----------
    df : dask.Dataframe
        Dataframe with columns "WCB_flag", "dp2h", "slan_400", "slan_600",
        "conv_400", "conv_600", "type" and group.
    group : string
        Column to group by, e.g. "time_after_ascent" to calculate the percentiles
        along.
    flag : string
        Column with bools.
    Returns
    -------
    Median, 25, 75 percentile trajectory
    """
    if flag is not None:
        df_tmp = df.loc[df[flag] == True]
    else:
        df_tmp = df.copy()
    drop_vars = [
        "WCB_flag",
        "dp2h",
        "slan_400",
        "slan_600",
        "conv_400",
        "conv_600",
        "type",
    ]
    groupy = df_tmp.drop_vars(drop_vars, errors="ignore").groupby(group)
    return (
        groupy.median(keep_attr=True),
        groupy.quantile(0.25, keep_attrs=True),
        groupy.quantile(0.75, keep_attrs=True),
    )


def norm_time_xarray(df, norm_col, group, columns=None, flag=None):
    """
    Return a view that consists only of entries that are flagged.
    Those are normed along norm_col such that every entry for every
    trajectory starts at norm_col==0. columns is a list of
    columns that the returned view shall have.
    if columns is None, take all columns. If flag is None, take all trajectories.

    Parameters
    ----------
    df : xarray.Dataframe
        Dataframe with columns norm_col, group and flag at minimum.
    norm_col : string
        Column where the minimum value where flag is true is substracted from.
        Usually it is "time"
    group : string
        Column for grouping operation, i.e. "trajectory".
    columns : list of string
        List of columns that the returned view shall have. If None is given,
        return all columns.
    flag : string
        Column with bools, i.e. start of an ascend.

    Returns
    -------
    xarray.Dataarray shifted norm_col.
    """
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
    """
    Create a median, 25, and 75 percentile trajectory.
    If flag is set, Use only areas where flag is true.
    Drops flag and type columns since those might not be valid anymore and
    need to be recalculated.

    Parameters
    ----------
    df : xarray.Dataframe
        Dataframe with column group.
    group : string
        Column to group by, e.g. "time_after_ascent" to calculate the percentiles
        along.
    flag : string
        Column with bools.
    Returns
    -------
    Median, 25, 75 percentile trajectory
    """
    if flag is not None:
        df_tmp = df.where(df[flag])
    else:
        df_tmp = df.copy()

    groupy = df_tmp.groupby(group)

    return groupy.median(), groupy.quantile(0.25), groupy.quantile(0.75)


def add_attrs(ds, ref_ds_path=None, attrs=None):
    """
    Add attributes to most columns with default values and values
    from a file in ref_ds_path if given.

    Parameters
    ----------
    ds : xarray.Dataset

    ref_ds_path : Path
        Path to NetCDF-file with attributes to read from.
    attrs : dict
        Dictionary with attributes to add. Must have the following keys:
        "ds": Generic attributes for the dataset
        "time": Attributes for the column "time"

    Returns
    -------
    Dataset with added attributes.
    """
    if ref_ds_path is not None:
        ds_2 = xr.open_dataset(ref_ds_path, decode_times=False, engine="netcdf4")
        duration = ds_2.attrs["duration_in_sec"]
        pollon = ds_2.attrs["pollon"]
        pollat = ds_2.attrs["pollat"]
        output_timestep_in_sec = ds_2.attrs["output_timestep_in_sec"]

        ds.attrs = {
            "duration_in_sec": duration,
            "pollon": pollon,
            "pollat": pollat,
            "output_timestep_in_sec": output_timestep_in_sec,
            "cloud_type": 2723,
        }

        ds["time"].attrs = ds_2["time"].attrs

    if attrs is not None:
        ds.attrs = attrs["ds"]
        ds["time"].attrs = attrs["time"]

    ds["time_after_ascent"].attrs = {
        "standard_name": "time_after_ascent",
        "long_name": "time after rapid ascent started",
        "units": "seconds since start of convective/slantwise ascent",
    }
    ds["lon"].attrs = {
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees_east",
    }
    ds["lat"].attrs = {
        "standard_name": "latitude",
        "long_name": "latitude",
        "units": "degrees_north",
    }
    ds["pressure"].attrs = {
        "standard_name": "air_pressure",
        "long_name": "pressure",
        "units": "Pa",
        "positive": "down",
        "axis": "Z",
    }
    ds["z"].attrs = {
        "standard_name": "height",
        "long_name": "height above mean sea level",
        "auxiliary_data": "yes",
        "units": "m AMSL",
    }
    ds["T"].attrs = {
        "standard_name": "air_temperature",
        "long_name": "temperature",
        "auxiliary_data": "yes",
        "units": "K",
    }
    ds["S"].attrs = {
        "standard_name": "saturation",
        "long_name": "saturation",
        "auxiliary_data": "yes",
        "units": "percentage",
    }
    ds["conv_400"].attrs = {
        "standard_name": "convective_400hPa_ascent",
        "long_name": "convective 400hPa ascent",
        "auxiliary_data": "yes",
    }
    ds["conv_600"].attrs = {
        "standard_name": "convective_600hPa_ascent",
        "long_name": "convective 600hPa ascent",
        "auxiliary_data": "yes",
    }
    ds["slan_400"].attrs = {
        "standard_name": "slantwise_400hPa_ascent",
        "long_name": "slantwise 400hPa ascent",
        "auxiliary_data": "yes",
    }
    ds["slan_600"].attrs = {
        "standard_name": "slantwise_600hPa_ascent",
        "long_name": "slantwise 600hPa ascent",
        "auxiliary_data": "yes",
    }
    ds["w"].attrs = {
        "standard_name": "ascend_velocity",
        "long_name": "ascend velocity",
        "auxiliary_data": "yes",
        "units": "m s^-1",
    }

    ds["QV"].attrs = {
        "standard_name": "specific_humidity",
        "long_name": "specific humidity",
        "auxiliary_data": "yes",
        "units": "kg kg^-1",
    }
    ds["QC"].attrs = {
        "standard_name": "mass_fraction_of_cloud_liquid_water_in_air",
        "long_name": "specific cloud liquid water content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1",
    }
    ds["QR"].attrs = {
        "standard_name": "mass_fraction_of_rain_in_air",
        "long_name": "specific rain content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1",
    }
    ds["QS"].attrs = {
        "standard_name": "mass_fraction_of_snow_in_air",
        "long_name": "specific snow content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1",
    }
    ds["QI"].attrs = {
        "standard_name": "mass_fraction_of_cloud_ice_in_air",
        "long_name": "specific cloud ice content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1",
    }
    ds["QG"].attrs = {
        "standard_name": "mass_fraction_of_graupel_in_air",
        "long_name": "specific graupel content",
        "auxiliary_data": "yes",
        "units": "kg kg^-1",
    }

    ds["QR_IN"].attrs = {
        "standard_name": "sedi_influx_of_rain",
        "long_name": "sedimentation (from above) of rain droplet mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1",
    }
    ds["QS_IN"].attrs = {
        "standard_name": "sedi_influx_of_snow",
        "long_name": "sedimentation (from above) of snow crystal mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1",
    }
    ds["QI_IN"].attrs = {
        "standard_name": "sedi_influx_of_cloud_ice",
        "long_name": "sedimentation (from above) of ice crystal mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1",
    }
    ds["QG_IN"].attrs = {
        "standard_name": "sedi_influx_of_graupel",
        "long_name": "sedimentation (from above) of graupel mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1",
    }

    ds["QR_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_rain",
        "long_name": "sedimentation of rain droplet mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1",
    }
    ds["QS_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_snow",
        "long_name": "sedimentation of snow crystal mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1",
    }
    ds["QI_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_cloud_ice",
        "long_name": "sedimentation of ice crystal mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1",
    }
    ds["QG_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_graupel",
        "long_name": "sedimentation of graupel mixing ratio",
        "auxiliary_data": "yes",
        "units": "kg kg^-1 s^-1",
    }

    ds["NCCLOUD"].attrs = {
        "standard_name": "specif_number_of_cloud_droplets_in_air",
        "long_name": "specific cloud droplet number",
        "auxiliary_data": "yes",
        "units": "kg^-1",
    }
    ds["NCRAIN"].attrs = {
        "standard_name": "specif_number_of_rain_drops_in_air",
        "long_name": "specific rain drop number",
        "auxiliary_data": "yes",
        "units": "kg^-1",
    }
    ds["NCSNOW"].attrs = {
        "standard_name": "specif_number_of_snow_flakes_in_air",
        "long_name": "specific snow flake number",
        "auxiliary_data": "yes",
        "units": "kg^-1",
    }
    ds["NCICE"].attrs = {
        "standard_name": "specif_number_of_cloud_ice_in_air",
        "long_name": "specific cloud ice number",
        "auxiliary_data": "yes",
        "units": "kg^-1",
    }
    ds["NCGRAUPEL"].attrs = {
        "standard_name": "specif_number_of_graupel_in_air",
        "long_name": "specific graupel number",
        "auxiliary_data": "yes",
        "units": "kg^-1",
    }

    ds["NR_IN"].attrs = {
        "standard_name": "sedi_influx_of_rain_number",
        "long_name": "sedimentation (from above) of specific rain drop number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }
    ds["NS_IN"].attrs = {
        "standard_name": "sedi_influx_of_snow_number",
        "long_name": "sedimentation (from above) of specific snow flake number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }
    ds["NI_IN"].attrs = {
        "standard_name": "sedi_influx_of_ics_number",
        "long_name": "sedimentation (from above) of specific cloud ice number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }
    ds["NG_IN"].attrs = {
        "standard_name": "sedi_influx_of_graupel_number",
        "long_name": "sedimentation (from above) of specific graupel number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }

    ds["NR_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_rain_number",
        "long_name": "sedimentation of rain droplet number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }
    ds["NS_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_snow_number",
        "long_name": "sedimentation of snow crystal number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }
    ds["NI_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_ice_number",
        "long_name": "sedimentation of ice crystal number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }
    ds["NG_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_graupel_number",
        "long_name": "sedimentation of graupel number",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }

    ds["Q_TURBULENCE"].attrs = {
        "standard_name": "turbulence_flux",
        "long_name": "flux from turbulence",
        "auxiliary_data": "yes",
        "units": "kg^-1 s^-1",
    }
    ds["type"].attrs = {
        "standard_name": "trajectory_type",
        "long_name": "trajectory type",
        "auxiliary_data": "yes",
    }

    return ds


def get_percentiles(
    store_path,
    version="no exclusions",
    fls=["conv_400", "conv_600", "slan_400", "slan_600"],
    file_list=None,
):
    """
    Load NetCDF-files and calculate the median, 25 and 75 percentile
    trajectories. Different versions include different trajectories.
    Stores those trajectories on disk.

    Parameters
    ----------
    store_path : path
        Path where to store the statistic trajectories.
    version : string
        Define which trajectories are included for the statistic trajectories.
        Options are:
        no exclusions: Get all trajectories with the corresponding flag
        excl other: Get all trajectories with the corresponding flag without
            the other type (ie flag is conv_400, excl slan_600)
        excl same: Get all trajectories with the corresponding flag without
            the same type (ie flag is conv_400, excl conv_600)
        excl all: Get all trajectories with the corresponding flag without
            all others (ie flag is conv_400, excl conv_600 and slan_600)
        conv_X and slan_X are mutually exclusive by definition already!
    fls : list of string
        List of flags to calculate the statistic trajectories for. Options are
        "conv_400", "conv_600", "slan_400", "slan_600"
    file_list : list of paths
        List of NetCDF-files to load for calculating the statistic trajectories
        for.
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
            ds_tmp = (
                xr.open_dataset(f, decode_times=False, engine="netcdf4")
                .to_dataframe()
                .reset_index()
            )
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
                    non_ids = np.unique(
                        ds_tmp.where(ds_tmp[fl_no] == True)["trajectory"]
                    )
                    ids = np.setdiff1d(ids, non_ids)
                    ds_tmp = ds_tmp.where(ds_tmp["trajectory"].isin(ids))
            elif version != "no exclusions":
                print(f"version {version} unknown!")
                return

            if ds is not None:
                ds = ds.append(ds_tmp)
            else:
                ds = ds_tmp

        t_c2 = timer()
        print(f"loading done in {t_c2-t_c} s")
        ds = ds.dropna()
        print(ds.describe())
        t_c = timer()

        medi, quan25, quan75 = get_statistics_pandas(ds, group="time_after_ascent")
        medi["trajectory"] = n
        n += 1
        quan25["trajectory"] = n
        n += 1
        quan75["trajectory"] = n
        n += 1
        t_c2 = timer()
        print(f"statistics done in {t_c2-t_c} s")

        t_c = timer()
        # Set flags for slantwise or convective parts
        for flag in ["conv_600", "conv_400", "slan_400", "slan_600"]:
            medi[flag] = False
            quan25[flag] = False
            quan75[flag] = False

        # The time needs to be adjusted
        def adjust(this_df):
            this_df = this_df.reset_index()
            this_df["ensemble"] = n % 3
            start_time = this_df["time"][0]
            end_time = start_time + 20 * len(this_df.index)
            this_df["time"] = np.arange(start_time, end_time, 20)
            return this_df

        medi = adjust(medi)
        quan25 = adjust(quan25)
        quan75 = adjust(quan75)

        medi = xr.Dataset.from_dataframe(
            medi.set_index(["ensemble", "trajectory", "time"])
        )
        quan25 = xr.Dataset.from_dataframe(
            quan25.set_index(["ensemble", "trajectory", "time"]).dropna()
        )
        quan75 = xr.Dataset.from_dataframe(
            quan75.set_index(["ensemble", "trajectory", "time"]).dropna()
        )

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
            conv_600 = (
                medi["pressure"]
                .rolling(dim={"time": window_conv_600}, min_periods=1)
                .reduce(differ, **{"hPa": 600})
                .fillna(False)
                .astype(dtype=bool)
            )
            medi = medi.assign(conv_600=conv_600)
            t_c2 = timer()
            print(f"Got conv_600 in {t_c2-t_c} s")

            t_c = timer()
            conv_600 = (
                quan25["pressure"]
                .rolling(dim={"time": window_conv_600}, min_periods=1)
                .reduce(differ, **{"hPa": 600, "debug": False})
                .fillna(False)
                .astype(dtype=bool)
            )
            quan25 = quan25.assign(conv_600=conv_600)
            t_c2 = timer()
            print(f"Got conv_600 in {t_c2-t_c} s")

            t_c = timer()
            conv_600 = (
                quan75["pressure"]
                .rolling(dim={"time": window_conv_600}, min_periods=1)
                .reduce(differ, **{"hPa": 600})
                .fillna(False)
                .astype(dtype=bool)
            )
            quan75 = quan75.assign(conv_600=conv_600)
            t_c2 = timer()
            print(f"Got conv_600 in {t_c2-t_c} s")
        elif fl == "conv_400":
            t_c = timer()
            conv_400 = (
                medi["pressure"]
                .rolling(dim={"time": window_conv_400}, min_periods=1)
                .reduce(differ, **{"hPa": 400})
                .fillna(False)
                .astype(dtype=bool)
            )
            medi = medi.assign(conv_400=conv_400)
            t_c2 = timer()
            print(f"Got conv_400 in {t_c2-t_c} s")

            t_c = timer()
            conv_400 = (
                quan25["pressure"]
                .rolling(dim={"time": window_conv_400}, min_periods=1)
                .reduce(differ, **{"hPa": 400})
                .fillna(False)
                .astype(dtype=bool)
            )
            quan25 = quan25.assign(conv_400=conv_400)
            t_c2 = timer()
            print(f"Got conv_400 in {t_c2-t_c} s")

            t_c = timer()
            conv_400 = (
                quan75["pressure"]
                .rolling(dim={"time": window_conv_400}, min_periods=1)
                .reduce(differ, **{"hPa": 400})
                .fillna(False)
                .astype(dtype=bool)
            )
            quan75 = quan75.assign(conv_400=conv_400)
            t_c2 = timer()
            print(f"Got conv_400 in {t_c2-t_c} s")
        elif fl == "slan_600":
            t_c = timer()
            slan_600 = (
                medi["pressure"]
                .rolling(dim={"time": window_slan_600}, min_periods=1)
                .reduce(differ_slan, **{"hPa": 600, "min_window": window_slan_600_min})
                .fillna(False)
                .astype(dtype=bool)
            )
            medi = medi.assign(slan_600=slan_600)
            t_c2 = timer()
            print(f"Got slan_600 in {t_c2-t_c} s")

            t_c = timer()
            slan_600 = (
                quan25["pressure"]
                .rolling(dim={"time": window_slan_600}, min_periods=1)
                .reduce(differ_slan, **{"hPa": 600, "min_window": window_slan_600_min})
                .fillna(False)
                .astype(dtype=bool)
            )
            quan25 = quan25.assign(slan_600=slan_600)
            t_c2 = timer()
            print(f"Got slan_600 in {t_c2-t_c} s")

            t_c = timer()
            slan_600 = (
                quan75["pressure"]
                .rolling(dim={"time": window_slan_600}, min_periods=1)
                .reduce(differ_slan, **{"hPa": 600, "min_window": window_slan_600_min})
                .fillna(False)
                .astype(dtype=bool)
            )
            quan75 = quan75.assign(slan_600=slan_600)
            t_c2 = timer()
            print(f"Got slan_600 in {t_c2-t_c} s")
        elif fl == "slan_400":
            t_c = timer()
            slan_400 = (
                medi["pressure"]
                .rolling(dim={"time": window_slan_400}, min_periods=1)
                .reduce(differ_slan, **{"hPa": 400, "min_window": window_slan_400_min})
                .fillna(False)
                .astype(dtype=bool)
            )
            medi = medi.assign(slan_400=slan_400)
            t_c2 = timer()
            print(f"Got slan_400 in {t_c2-t_c} s")

            t_c = timer()
            slan_400 = (
                quan25["pressure"]
                .rolling(dim={"time": window_slan_400}, min_periods=1)
                .reduce(differ_slan, **{"hPa": 400, "min_window": window_slan_400_min})
                .fillna(False)
                .astype(dtype=bool)
            )
            quan25 = quan25.assign(slan_400=slan_400)
            t_c2 = timer()
            print(f"Got slan_400 in {t_c2-t_c} s")

            t_c = timer()
            slan_400 = (
                quan75["pressure"]
                .rolling(dim={"time": window_slan_400}, min_periods=1)
                .reduce(differ_slan, **{"hPa": 400, "min_window": window_slan_400_min})
                .fillna(False)
                .astype(dtype=bool)
            )
            quan75 = quan75.assign(slan_400=slan_400)
            t_c2 = timer()
            print(f"Got slan_400 in {t_c2-t_c} s")
        print("Rolling done")

        t_c = timer()
        # Set time_after_ascent to zero where flag first occurs
        def adjust_ascent_time(this_ds):
            start_idx = np.where(
                ~np.isnan(this_ds.where(this_ds[fl] == True)["time_after_ascent"])
            )[2][0]
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
        print(f"Set time after ascent in {t_c2-t_c} s")

        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in medi.data_vars}
        medi.to_netcdf(
            store_path + fl + "_median.nc_wcb",
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        quan25.to_netcdf(
            store_path + fl + "_quan25.nc_wcb",
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        quan75.to_netcdf(
            store_path + fl + "_quan75.nc_wcb",
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        print("storing done")


def add_flags(path="/data/project/wcb/netcdf/traj_stats/"):
    """
    Load all files in path, calculate where the ascend starts based on the
    filename and store it on disk with "2" added to the name.

    path : path
        Path to folder with NetCDF-files.
    """
    file_list = []
    for f in os.listdir(path):
        file_list.append(os.path.join(path, f))

    file_list = np.sort(np.asarray(file_list))
    for f in file_list:
        with xr.open_dataset(f, engine="netcdf4") as ds:
            if "conv_400" in f:
                conv_400 = (
                    ds["P"]
                    .rolling(dim={"time": window_conv_400}, min_periods=1)
                    .reduce(differ, **{"hPa": 400})
                    .fillna(False)
                    .astype(dtype=bool)
                )
                ds = ds.assign(conv_400=conv_400)
            elif "conv_600" in f:
                conv_600 = (
                    ds["P"]
                    .rolling(dim={"time": window_conv_600}, min_periods=1)
                    .reduce(differ, **{"hPa": 600})
                    .fillna(False)
                    .astype(dtype=bool)
                )
                ds = ds.assign(conv_600=conv_600)
            elif "slan_400" in f:
                slan_400 = (
                    ds["P"]
                    .rolling(dim={"time": window_slan_400}, min_periods=1)
                    .reduce(
                        differ_slan, **{"hPa": 400, "min_window": window_slan_400_min}
                    )
                    .fillna(False)
                    .astype(dtype=bool)
                )
                ds = ds.assign(slan_400=slan_400)
            elif "slan_600" in f:
                slan_600 = (
                    ds["P"]
                    .rolling(dim={"time": window_slan_600}, min_periods=1)
                    .reduce(
                        differ_slan, **{"hPa": 600, "min_window": window_slan_600_min}
                    )
                    .fillna(False)
                    .astype(dtype=bool)
                )
                ds = ds.assign(slan_600=slan_600)
            ds.to_netcdf(f + "2")


def get_percentiles_slan_600():
    """
    Debug method. TODO: Delete this?

    """
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
        for f in file_list[i : i + 3]:
            ds_tmp = xr.open_dataset(f, engine="netcdf4").to_dataframe().reset_index()
            ids = ds_tmp.loc[ds_tmp[fl] == True]["id"]
            ds_tmp = ds_tmp.loc[ds_tmp["id"].isin(ids)]
            if ds is None:
                ds = ds_tmp
            else:
                ds = ds.append(ds_tmp)

        t_c2 = timer()
        print(f"loading done in {t_c2-t_c} s")
        # Make id and time columns instead of MultiIndex
        t_c = timer()
        normed = norm_time(ds, "time", "id", None, fl)
        t_c2 = timer()
        print(f"norming done in {t_c2-t_c} s")
        for flag in ["conv_600", "conv_400", "slan_400", "slan_600"]:
            normed[flag] = False
        t_c = timer()
        xr.Dataset.from_dataframe(
            normed.reset_index().set_index(["time", "id"])
        ).to_netcdf(store_path_tmp + fl + "_normed_{}".format(i // 3))
        t_c2 = timer()
        print(f"storing done in {t_c2-t_c} s")

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
            ds_tmp = xr.open_dataset(f, engine="netcdf4").to_dataframe().reset_index()
            # Load only certain timesteps
            ds_tmp = ds_tmp.loc[
                (ds_tmp.time > (min_time + n * i))
                & (ds_tmp.time < (min_time + n * (i + 1)))
            ]

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
        print(f"statistics done in {t_c2-t_c} s")
        t_c = timer()
        # Set flags for slantwise or convective parts
        for flag in ["conv_600", "conv_400", "slan_400", "slan_600"]:
            medi[flag] = False
            quan25[flag] = False
            quan75[flag] = False

        medi = xr.Dataset.from_dataframe(
            medi.reset_index().set_index(["time", "id"]).dropna()
        )
        quan25 = xr.Dataset.from_dataframe(
            quan25.reset_index().set_index(["time", "id"]).dropna()
        )
        quan75 = xr.Dataset.from_dataframe(
            quan75.reset_index().set_index(["time", "id"]).dropna()
        )

        type_name = "Slantwise 600hPa"
        medi["type"] = type_name + " 50. Quantile"
        quan25["type"] = type_name + " 25. Quantile"
        quan75["type"] = type_name + " 75. Quantile"
        t_c2 = timer()
        print(f"Resetting done in {t_c2-t_c} s")

        t_c = timer()
        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in ds.data_vars}
        medi.to_netcdf(
            store_path_tmp + "medi/" + fl + "_median_{}.nc_wcb".format(i),
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        quan25.to_netcdf(
            store_path_tmp + "quan25/" + fl + "_quan25_{}.nc_wcb".format(i),
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        quan75.to_netcdf(
            store_path_tmp + "quan75/" + fl + "_quan75_{}.nc_wcb".format(i),
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        t_c2 = timer()
        print(f"storing done in {t_c2-t_c} s")


def merge_stuff(
    input_path,
    output,
    ensembles=["conv_400", "conv_600", "slan_400", "slan_600"],
    dropna=None,
    max_diff=None,
    complevel=9,
):
    """
    Merge statistical trajectories to a single file.
    Each category given in ensembles will be merged to a single ensemble.
    The time coordinates are set %20 to align them.
    Longitude and Latitude are set, such that they do not deviate too much
    for plotting in 3D.


    Parmaters
    ---------
    input_path: string
        Path to trajectory files.
    output: string
        Path + name where new file shall be stored.
    ensembles: list of string
        Categories to merge into a single ensembles (must be present in filename).
    dropna: string
        "drop": Remove all NaNs such that only coordinates with values are left.
        "fill": Fill values with first and then last valid numbers
        "999": Fill values with -999.99
        None: Do nothing
    max_diff: float
        Maximum allowed difference for consecutive timesteps of longitude and
        latitude. Use 0.012 for roughly 180 km/h max velocity. If None is given
        the coordinates will not be touched.
    complevel: int
        Compression level to use (0-9)

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
                ds_tmp = xr.open_dataset(f, decode_times=False, engine="netcdf4")
                # Unfortunately, we need to make sure the time coordinates
                # are somewhat similar ie modulo 20
                # Since we are merging "representatives" the exact time
                # should not matter
                ds_tmp["time"] = ds_tmp["time"] - ds_tmp["time"] % 20
                max_t = 25000
                if "slan" in f:
                    max_t = 82000
                ds_tmp = ds_tmp.where(ds_tmp["time_after_ascent"] <= max_t)
                if max_diff is not None:

                    def adjust(ds_tmp, key):
                        ll = np.asarray(ds_tmp[key])
                        before = ll[0, 0, 0]
                        for i, l in enumerate(ll[0, 0, :]):
                            if np.abs(before - l) > max_diff:
                                ll[0, 0, i] = before + np.sign(l - before) * max_diff
                                before = ll[0, 0, i]
                            else:
                                before = l
                        ds_tmp[key] = (ds_tmp[key].dims, ll)
                        return ds_tmp

                    ds_tmp = adjust(ds_tmp, "lon")
                    ds_tmp = adjust(ds_tmp, "lat")

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
            "forecast_inittime": "2016-09-20 00:00:00",
        },
        "ds": {
            "duration_in_sec": (max_t - min_t).item(),
            "pollon": 160.0,
            "pollat": 51.0,
            "output_timestep_in_sec": 20,
            "cloud_type": 2723,
        },
    }
    merged_ds = add_attrs(merged_ds, attrs=attrs)
    if dropna == "drop":
        merged_ds = merged_ds.dropna(dim="time", how="any")
        print(merged_ds)
    elif dropna == "fill":
        params = ["pressure", "lat", "lon", "QV", "QC", "QR"]
        for par in params:
            merged_ds[par] = merged_ds[par].ffill("time").bfill("time")
    elif dropna == "999":
        merged_ds = merged_ds.fillna(-999.99)

    comp = dict(zlib=True, complevel=complevel)
    encoding = {var: comp for var in merged_ds.data_vars}

    merged_ds.to_netcdf(
        output,
        compute=True,
        engine="netcdf4",
        encoding=encoding,
        format="NETCDF4",
        mode="w",
    )
    return merged_ds


# def get_all_stats(mean_dfs, out_params, others_path, sub_folders,
#     ratio_type=["per_timestep", "window"], in_params=None):
#     """
#     TODO: Delete this ?

#     Parameters
#     ----------
#     mean_dfs :
#     out_params :
#     others_path :
#     sub_folders :
#     ratio_type :
#     in_params :
#     """
#     df_dic = {"MSE": np.asarray([]),
#             "Max Error": np.asarray([]),
#             "MSE (no zero)": np.asarray([]),
#             "Mean Error": np.asarray([]),
#             "Mean Absolute Sensitivity": np.asarray([]),
#             "Max Sensitivity": np.asarray([]),
#             "Median Absolute Sensitivity": np.asarray([]),
#             "Mean Absolute Sensitivity (no zero)": np.asarray([]),
#             "Median Absolute Sensitivity (no zero)": np.asarray([]),
#             "Mean Absolute Sensitivity (no zero sens)": np.asarray([]),
#             "Median Absolute Sensitivity (no zero sens)": np.asarray([]),
#             "Output Parameter": np.asarray([]),
#             "Perturbed Parameter": np.asarray([]),
#             "Ratio Type": np.asarray([])}

#     folder_mean_dic = {}
#     folder_sens_dic = {}
#     folder_mean_df = {}
#     # Pre calculate the sensitivities from the not perturbed version
#     for folder in sub_folders:

#         # Load mean for this

#         mean = deriv_dask.cache
#         sens_dic = {}
#         mean_dic = {}
#         for out_param in out_params:
#             mean_dic[out_param] = np.reshape(np.asarray(
#                 mean.loc[mean["Output Parameter"] == out_param][out_param]),
#                 (len(mean.loc[mean["Output Parameter"] == out_param][out_param].index), 1, 1))
#             # print(mean["time"])
#             if isinstance(ratio_type, list):
#                 sens_dic[out_param] = {}
#                 for rt in ratio_type:
#                     sens_dic[out_param][rt] = self._recalc_ratios(mean.loc[mean["Output Parameter"] == out_param],
#                         rt, ratio_window, in_params)
#             else:
#                 sens_dic[out_param] = self._recalc_ratios(mean.loc[mean["Output Parameter"] == out_param],
#                     ratio_type, ratio_window, in_params)
#         folder_mean_dic[folder] = mean_dic
#         folder_sens_dic[folder] = sens_dic
#         folder_mean_df[folder] = mean
#     # Load the perturbed versions
#     for perturbed_param in in_params:
#         # Load the dfs with the perturbed ensembles
#         for folder in sub_folders:
#             min_x = np.min(folder_mean_df[folder]["time"])
#             max_x = np.max(folder_mean_df[folder]["time"])
#             ds = xr.open_dataset(
#                 others_path + "/" + folder + "/d" + perturbed_param + ".nc_wcb",
#                 decode_times=False)[out_params]
#             ds = ds.sel(time=np.arange(min_x, max_x+20, 20)).compute()
#             for out_param in out_params:
#                 def add_to_df(sens_df, rt):
#                     #

#                 if isinstance(ratio_type, list):
#                     for rt in ratio_type:
#                         add_to_df(sens_dic[out_param][rt], rt)
#                 else:
#                     add_to_df(sens_dic[out_param], ratio_type)

#     return pandas.DataFrame.from_dict(df_dic)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Calculate statistic trajectories from "path" for given types via
        "flags" and store them to "store_path".
        """
    )
    parser.add_argument(
        "--path",
        default="../data/sim_processed/conv_400_0_t000000_p001_mult_outSat_sbShape_sbConv/",
        # my path /data/project/wcb/netcdf/vladiana_met/
        help="""
        Path to NetCDF-files stored in Met3D style.
        """,
    )
    parser.add_argument(
        "--store_path",
        default="../data/median/",
        # my path "/data/project/wcb/netcdf/vladiana_met_stats/"
        help="""
        Path where statistic trajectories shall be stored.
        """,
    )
    parser.add_argument(
        "--flags",
        default=["conv_400", "conv_600", "slan_400", "slan_600"],
        type=str,
        nargs="+",
        help="""
        Flags for the type of trajectories to look for.
        """,
    )
    parser.add_argument(
        "--version",
        default="no exclusions",
        type=str,
        help="""
        Define which trajectories are included for the statistic trajectories.
        Options are:
        no exclusions: Get all trajectories with the corresponding flag
        excl other: Get all trajectories with the corresponding flag without
            the other type (ie flag is conv_400, excl slan_600)
        excl same: Get all trajectories with the corresponding flag without
            the same type (ie flag is conv_400, excl conv_600)
        excl all: Get all trajectories with the corresponding flag without
            all others (ie flag is conv_400, excl conv_600 and slan_600)
        conv_X and slan_X are mutually exclusive by definition already!
        """,
    )
    args = parser.parse_args()

    file_list = []
    for f in os.listdir(args.path):
        file_list.append(os.path.join(args.path, f))

    file_list = np.sort(np.asarray(file_list))
    do_the_stuff_more(args.store_path, args.version, args.flags, file_list)
