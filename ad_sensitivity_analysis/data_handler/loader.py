"""Functions to load different datasets.

"""
import os

import numpy as np
from tqdm.auto import tqdm
import xarray as xr

from ad_sensitivity_analysis.data_handler.filter_data import filter_trajectories
from ad_sensitivity_analysis.data_handler.transform import (
    get_average,
    add_liquid_content,
    add_cold_content,
    add_sepcific_humidity,
)


def load_dataset(in_path):
    """
    Load a dataset.

    Parameters
    ----------
    in_path : string
        Path to the NetCDF-file.

    Returns
    -------
    xr.Dataset with trajectory data.
    """
    return xr.open_dataset(in_path, decode_times=False, engine="netcdf4")


def load_dataset_part(
    f,
    only_asc600=False,
    only_phase=None,
    inoutflow_time=-1,
    load_params=None,
    lat_bug=False,
):
    """
    Load an xarray.Dataset and filter it for certain columns or time steps or phases if needed.

    Parameters
    ----------
    f : string
        Path and name of the file to load
    only_asc600 : bool
        Consider only time steps during the ascend.
    only_phase : string
        Consider only time steps with the given phase. Can be combined with only_asc600 or inoutflow_time.
        Possible values are "warm phase", "mixed phase", "ice phase", "neutral phase".
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.
    load_params : list of strings
        If given, load only the provided columns.

    Returns
    -------
    Loaded and, if needed, filtered xarray.Dataset.
    """

    if load_params is not None:
        additional_vars = []
        if only_phase is not None and "phase" not in load_params:
            additional_vars.append("phase")
        if inoutflow_time > 0 and "asc600" not in load_params:
            additional_vars.append("asc600")
        if len(additional_vars) > 0:
            ds = xr.open_dataset(f, decode_times=False, engine="netcdf4")[
                load_params + additional_vars
            ]
        else:
            ds = xr.open_dataset(f, decode_times=False, engine="netcdf4")[load_params]
    else:
        ds = xr.open_dataset(f, decode_times=False, engine="netcdf4")

    return filter_trajectories(
        ds=ds,
        only_asc600=only_asc600,
        only_phase=only_phase,
        lat_bug=lat_bug,
    )


def load_data(
    file_path,
    x,
    only_asc600=False,
    inoutflow_time=-1,
    phase=None,
    averages=False,
    verbose=False,
):
    """
    Load multiple NetCDF-files as xarray.Dataset, filter them if necessary and concatenate all
    data useful for other steps.

    Parameters
    ----------
    file_path : string
        Path to trajectories from a sensitivity analysis.
    x : string or list of strings
        Model state variable(s) or model parameter(s) for calculating the clusters.
    only_asc600 : bool
        Consider only time steps during the fastest 600 hPa ascent.
    inoutflow_time : int
        Consider only time steps during the fastest 600 hPa ascent and this many timesteps before
        and after the ascent (in- and outflow).
    averages : bool
        If true, calculate the averages over time when each file is loaded. Otherwise tries to concatenate
        all data in the end which may need lots of memory if inoutflow_time > 0 or only_asc600 == True
        or phase is not None.
    phase : int
        Only use the given phase.
        0: warm phase
        1: mixed phase
        2: ice phase
        3: neutral phase
    verbose : bool
        If true, get more output.

    Returns
    -------
    xarray.Dataset with concatenated values for all trajectories. Includes a new dimension 'file' to
    backtrack any results to the correct file and trajectory.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    list_of_arrays = []
    if isinstance(x, str):
        variables = [x, "asc600", "phase"]
    else:
        variables = list(x) + ["asc600", "phase"]

    for f in tqdm(files) if verbose else files:
        ds = load_dataset_part(
            file_path + f,
            load_params=variables,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            only_phase=phase,
            lat_bug=True,
        )
        ds = ds[x].dropna("time", how="all")
        ds = ds.expand_dims({"file": [f]})

        if averages:
            ds = get_average(ds)
            list_of_arrays.append(ds)
        elif len(ds["time"]) > 0:
            list_of_arrays.append(ds)
    return xr.concat(list_of_arrays, dim="file", join="outer", combine_attrs="override")


def _simulation_and_orig_load_aux(
    f, cols_final_sim, traj, limits_x, get_limits, verbosity
):
    """

    Parameters
    ----------
    f
    cols_final_sim
    traj
    limits_x
    get_limits
    verbosity

    Returns
    -------

    """
    df = None
    ds = xr.open_dataset(f, decode_times=False, engine="netcdf4")[cols_final_sim]
    if isinstance(traj, int):
        ds = ds.loc[{"trajectory": ds["trajectory"][traj]}]
    if df is not None:
        df = df.append(ds.to_dataframe())
    else:
        if get_limits:
            limits_x[0] = np.min(ds["time_after_ascent"]).item()
            limits_x[1] = np.max(ds["time_after_ascent"]).item()
        else:
            ds = ds.where(
                (
                    (ds["time_after_ascent"] >= limits_x[0])
                    & (ds["time_after_ascent"] <= limits_x[1])
                ),
                drop=True,
            )
        if verbosity > 2:
            print(f"Got min time = {limits_x[0]} s, max time = {limits_x[1]} s.")
        df = ds.to_dataframe()
    return df, limits_x


def load_simulation_and_orig(
    data_cosmo_path,
    data_sim_path,
    verbosity=0,
    traj=None,
):
    """

    Parameters
    ----------
    data_cosmo_path
    data_sim_path
    verbosity
    traj

    Returns
    -------

    """
    file_list_cosmo = [
        os.path.join(data_cosmo_path, f) for f in os.listdir(data_cosmo_path)
    ]
    file_list_cosmo.sort()
    cols_final = []
    for col in list(
        xr.open_dataset(file_list_cosmo[0], decode_times=False, engine="netcdf4").keys()
    ):
        if (
            "_IN" not in col
            and "WCB_flag" not in col
            and "dp2h" not in col
            and "Q_TURBULENCE" not in col
            and "type" not in col
        ):
            cols_final.append(col)

    cols_final_sim = cols_final.copy()
    cols_final_sim.append("QH")
    cols_final_sim.append("QH_OUT")
    file_list = [os.path.join(data_sim_path, f) for f in os.listdir(data_sim_path)]
    file_list.sort()

    if verbosity > 2:
        print(file_list)
        print("original:")
        print(file_list_cosmo)
    for f in tqdm(file_list) if verbosity > 0 else file_list:
        df_sim, limits_x = _simulation_and_orig_load_aux(
            f, cols_final_sim, traj, [], True, verbosity
        )
    df_sim["Simulation"] = "Our Sim."
    df_cosmo = None
    for f in tqdm(file_list_cosmo) if verbosity > 0 else file_list_cosmo:
        df_cosmo, _ = _simulation_and_orig_load_aux(
            f, cols_final_sim, traj, limits_x, False, verbosity
        )
    df_cosmo["Simulation"] = "Original"
    df_cosmo["QH"] = 0
    df_cosmo["QH_OUT"] = 0

    # drop all the NaN entries now that coordinates are collapsed to columns
    df_cosmo = df_cosmo.loc[
        (df_cosmo["time_after_ascent"] >= limits_x[0])
        & (df_cosmo["time_after_ascent"] < limits_x[1])
    ]
    df = df_cosmo.append(df_sim)

    df["time_after_ascent_h"] = df["time_after_ascent"] / 3600
    df["pressure_hPa"] = df["pressure"] / 100
    if verbosity > 2:
        print("Adding liquid and cold contents and specific humidity")
    df = add_liquid_content(df)
    df = add_cold_content(df)
    df = add_sepcific_humidity(df)

    df["Q_total"] = df["Q_cold"] + df["Q_liquid"] + df["Specific_Humidity"]
    if verbosity > 2:
        print("Finished loading data")
    return df
