import warnings

warnings.simplefilter(action="ignore", category=RuntimeWarning)

import hvplot.xarray  # noqa
import itertools
import holoviews as hv
from holoviews import opts
from matplotlib.colors import SymLogNorm, NoNorm
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import os
import panel as pn
import seaborn as sns
from tqdm.auto import tqdm
import xarray as xr

try:
    from scripts.latexify import param_id_map
except:
    from latexify import param_id_map

hv.extension("matplotlib")
pn.extension()


def filter_trajectories(
    ds,
    only_asc600=False,
    sens_model_state_ids=[],
    inoutflow_time=-1,
    min_pressure=None,
    max_pressure=None,
):
    """

    Parameters
    ----------
    ds
    only_asc600
    sens_model_state_ids
    inoutflow_time
    min_pressure
    max_pressure

    Returns
    -------

    """
    if inoutflow_time > 0 and "asc600" in ds:
        ds_flow = ds.where(ds["asc600"] == 1)["asc600"]
        ds_flow = ds_flow.rolling(
            time=inoutflow_time * 2,  # once for inflow, another for outflow
            min_periods=1,
            center=True,
        ).reduce(np.nanmax)
        ds = ds.where(ds_flow == 1)
    elif only_asc600 and "asc600" in ds:
        ds = ds.where(ds["asc600"] == 1)
    if min_pressure is not None:
        ds = ds.where(ds["pressure"] >= min_pressure)
    if max_pressure is not None:
        ds = ds.where(ds["pressure"] <= max_pressure)
    if len(sens_model_state_ids) > 0 and "Output_Parameter_ID" in ds:
        ds = ds.sel({"Output_Parameter_ID": sens_model_state_ids})
    # There seems to be a bug with some output data where the longitude
    # and latitude happen to be zero when the trajectory is already finished.
    # Since this is far away from our domain, we can savely delete that-
    ds["lon"] = ds["lon"].where(ds["lon"] != 0)
    ds["lat"] = ds["lat"].where(ds["lat"] != 0)
    return ds


def load_lon_lat_time(
    file_path,
    sens_model_states=[],
    only_asc600=False,
    inoutflow_time=-1,
    min_pressure=None,
    max_pressure=None,
    n_lons=100,
    n_lats=100,
    delta_time=None,
    relative_lon_lat=False,
    verbose=False,
):
    """

    Parameters
    ----------
    file_path
    sens_model_states
    only_asc600
    inoutflow_time
    min_pressure
    max_pressure
    n_lons
    n_lats
    delta_time
    relative_lon_lat : bool
        Use longitude and latitude degrees relative to the start of the ascent.
    verbose

    Returns
    -------

    """
    if verbose:
        print("Get the longitude, latitude and time values to create gridded data.")
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    lon_limits = []
    lat_limits = []
    time_limits = []
    sens_model_state_ids = []
    if sens_model_states is not None:
        for state in sens_model_states:
            sens_model_state_ids.append(
                np.argwhere(np.asarray(param_id_map) == state).item()
            )

    additional_vars = ["lon", "lat", "asc600", "pressure"]
    if relative_lon_lat:
        additional_vars.extend(["relative_lon", "relative_lat"])
        lon_name = "relative_lon"
        lat_name = "relative_lat"
    else:
        lon_name = "lon"
        lat_name = "lat"
    # Get the limits on lon and lat
    for f in tqdm(files) if verbose else files:
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[
            additional_vars
        ]
        ds = filter_trajectories(
            ds=ds,
            sens_model_state_ids=sens_model_state_ids,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
        )
        if len(lon_limits) > 0:
            v = ds[lon_name].min().values.item()
            if v < lon_limits[0]:
                lon_limits[0] = v
            v = ds[lon_name].max().values.item()
            if v > lon_limits[1]:
                lon_limits[1] = v
            v = ds[lat_name].min().values.item()
            if v < lat_limits[0]:
                lat_limits[0] = v
            v = ds[lat_name].max().values.item()
            if v > lat_limits[1]:
                lat_limits[1] = v
            v = ds["time"].min().values.item()
            if v < time_limits[0]:
                time_limits[0] = v
            v = ds["time"].max().values.item()
            if v > time_limits[1]:
                time_limits[1] = v

        else:
            lon_limits.append(ds[lon_name].min().values.item())
            lon_limits.append(ds[lon_name].max().values.item())
            lat_limits.append(ds[lat_name].min().values.item())
            lat_limits.append(ds[lat_name].max().values.item())
            time_limits.append(ds["time"].min().values.item())
            time_limits.append(ds["time"].max().values.item())

    delta_lon = (lon_limits[1] - lon_limits[0]) / (n_lons - 2)
    lons = np.arange(
        lon_limits[0] - delta_lon, lon_limits[1] + 1.8 * delta_lon, delta_lon
    )
    delta_lat = (lat_limits[1] - lat_limits[0]) / (n_lats - 2)
    lats = np.arange(
        lat_limits[0] - delta_lat, lat_limits[1] + 1.8 * delta_lat, delta_lat
    )
    if delta_time is not None:
        times = np.arange(time_limits[0], time_limits[1] + 0.8 * delta_time, delta_time)
    else:
        times = time_limits
    return lons, delta_lon, lats, delta_lat, times


def load_counts_means(
    file_path,
    variables,
    additional_vars,
    lons,
    delta_lon,
    lats,
    delta_lat,
    sens_model_states=[],
    only_asc600=False,
    inoutflow_time=-1,
    min_pressure=None,
    max_pressure=None,
    pressure_levels=None,
    time_levels=None,
    verbose=False,
    process_file=None,
    relative_lon_lat=False,
):
    """

    Parameters
    ----------
    file_path
    variables
    additional_vars
    lons
    delta_lon
    lats
    delta_lat
    sens_model_states
    only_asc600
    inoutflow_time
    min_pressure
    max_pressure
    pressure_levels
    time_levels
    verbose
    process_file
    relative_lon_lat

    Returns
    -------

    """
    if verbose:
        print(
            "Calculating the mean and number of occurrences at each grid point in this function."
        )
    lon_name = "lon"
    lat_name = "lat"
    if relative_lon_lat:
        lon_name = "relative_lon"
        lat_name = "relative_lat"
    n_lons = len(lons) - 1
    n_lats = len(lats) - 1
    if pressure_levels is not None:
        n_press = len(pressure_levels) - 1
    else:
        n_press = 0
    if time_levels is not None:
        n_times = len(time_levels) - 1
    else:
        n_times = 1
    sens_model_state_ids = []
    if sens_model_states is not None:
        for state in sens_model_states:
            sens_model_state_ids.append(
                np.argwhere(np.asarray(param_id_map) == state).item()
            )

    if pressure_levels is not None and n_times > 1:
        var_shapes = (n_times, n_press, n_lats, n_lons)
        sens_shapes = (len(sens_model_state_ids), n_times, n_press, n_lats, n_lons)
    elif pressure_levels is not None:
        var_shapes = (n_times, n_press, n_lats, n_lons)
        sens_shapes = (len(sens_model_state_ids), n_times, n_press, n_lats, n_lons)
    elif time_levels is not None and n_times > 1:
        var_shapes = (n_times, n_lats, n_lons)
        sens_shapes = (len(sens_model_state_ids), n_times, n_lats, n_lons)
    else:
        var_shapes = (n_lats, n_lons)
        sens_shapes = (len(sens_model_state_ids), n_lats, n_lons)
    if len(sens_model_state_ids) <= 1:
        sums = {v: np.zeros(var_shapes) for v in variables}
    else:
        sums = {}
        for v in variables:
            if v[0] != "d" or v == "deposition":
                sums[v] = np.zeros(var_shapes)
            else:
                sums[v] = np.zeros(sens_shapes)

    if verbose:
        print(
            "Calculating the sums for means and the number of occurrences at each grid point."
        )
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    if process_file is not None:
        files = [files[process_file]]
        if verbose:
            print(f"Processing {files[0]}")
    counts = np.zeros(var_shapes)
    for f in tqdm(files) if verbose else files:
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[
            additional_vars + variables
        ]
        ds = filter_trajectories(
            ds=ds,
            sens_model_state_ids=sens_model_state_ids,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
        )
        if pressure_levels is not None and n_times > 1:
            raise NotImplementedError(
                "Integrating for specific time levels is not yet supported"
            )
        elif pressure_levels is not None:
            for i in range(n_press):
                counts_tmp, _, _ = np.histogram2d(
                    y=ds[lon_name]
                    .where(
                        (ds["pressure"] >= pressure_levels[i])
                        & (ds["pressure"] < pressure_levels[i + 1])
                    )
                    .values.flatten(),
                    x=ds[lat_name]
                    .where(
                        (ds["pressure"] >= pressure_levels[i])
                        & (ds["pressure"] < pressure_levels[i + 1])
                    )
                    .values.flatten(),
                    bins=(lats, lons),
                )

                counts[0, i] += counts_tmp
        elif time_levels is not None and n_times > 1:
            raise NotImplementedError(
                "Integrating for specific time levels is not yet supported"
            )
        else:
            counts_tmp, _, _ = np.histogram2d(
                y=ds[lon_name].values.flatten(),
                x=ds[lat_name].values.flatten(),
                bins=(lats, lons),
            )
            counts += counts_tmp

        # Get the actual values and calculate the mean at the end
        lo_idx = 0
        la_idx = 0
        for la, lo in (
            tqdm(
                itertools.product(lats[:-1], lons[:-1]),
                total=n_lons * n_lats,
                leave=False,
            )
            if verbose
            else itertools.product(lats[:-1], lons[:-1])
        ):
            if pressure_levels is not None and n_times > 1:
                raise NotImplementedError(
                    "Integrating for specific time levels is not yet supported"
                )
            elif pressure_levels is not None:
                for p_i in range(n_press):
                    idx_where = (
                        (ds[lon_name] >= lo)
                        & (ds[lon_name] < lo + delta_lon)
                        & (ds[lat_name] >= la)
                        & (ds[lat_name] < la + delta_lat)
                        & (ds["pressure"] >= pressure_levels[p_i])
                        & (ds["pressure"] < pressure_levels[p_i + 1])
                    )
                    for v in variables:
                        if (
                            v[0] != "d"
                            or v == "deposition"
                            or len(sens_model_states) == 1
                        ):
                            sums[v][0, p_i, la_idx, lo_idx] += np.nansum(
                                ds[v].where(idx_where)
                            )
                        else:
                            sens_idx = 0
                            for id, s in zip(sens_model_state_ids, sens_model_states):
                                ds_tmp = ds.sel({"Output_Parameter_ID": id})
                                idx_where = (
                                    (ds_tmp[lon_name] >= lo)
                                    & (ds_tmp[lon_name] < lo + delta_lon)
                                    & (ds_tmp[lat_name] >= la)
                                    & (ds_tmp[lat_name] < la + delta_lat)
                                    & (ds_tmp["pressure"] >= pressure_levels[p_i])
                                    & (ds_tmp["pressure"] < pressure_levels[p_i + 1])
                                )
                                sums[v][sens_idx, 0, p_i, la_idx, lo_idx] += np.nansum(
                                    ds_tmp[v].where(idx_where)
                                )
                                sens_idx += 1

            elif time_levels is not None and n_times > 1:
                raise NotImplementedError(
                    "Integrating for specific time levels is not yet supported"
                )
            else:
                idx_where = (
                    (ds[lon_name] >= lo)
                    & (ds[lon_name] < lo + delta_lon)
                    & (ds[lat_name] >= la)
                    & (ds[lat_name] < la + delta_lat)
                )
                for v in variables:
                    if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
                        sums[v][la_idx, lo_idx] += np.nansum(ds[v].where(idx_where))
                    else:
                        sens_idx = 0
                        for id, s in zip(sens_model_state_ids, sens_model_states):
                            ds_tmp = ds.sel({"Output_Parameter_ID": id})
                            idx_where = (
                                (ds_tmp[lon_name] >= lo)
                                & (ds_tmp[lon_name] < lo + delta_lon)
                                & (ds_tmp[lat_name] >= la)
                                & (ds_tmp[lat_name] < la + delta_lat)
                            )
                            sums[v][sens_idx, la_idx, lo_idx] += np.nansum(
                                ds_tmp[v].where(idx_where)
                            )
                            sens_idx += 1
            lo_idx += 1
            if lo_idx == n_lons:
                la_idx += 1
                lo_idx = 0
    if verbose:
        print("Divide the sum to get the expected value (here: the mean)")
    lo_idx = 0
    la_idx = 0
    for _, _ in (
        tqdm(itertools.product(lats[:-1], lons[:-1]), total=n_lons * n_lats)
        if verbose
        else itertools.product(lats[:-1], lons[:-1])
    ):
        for v in sums:
            if pressure_levels is not None and n_times > 1:
                raise NotImplementedError(
                    "Integrating for specific time levels is not yet supported"
                )
            elif pressure_levels is not None:
                for p_i in range(n_press):
                    if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
                        sums[v][0, p_i, la_idx, lo_idx] /= counts[
                            0, p_i, la_idx, lo_idx
                        ]
                    else:
                        sums[v][:, 0, p_i, la_idx, lo_idx] /= counts[
                            0, p_i, la_idx, lo_idx
                        ]
            elif time_levels is not None and n_times > 1:
                raise NotImplementedError(
                    "Integrating for specific time levels is not yet supported"
                )
            else:
                if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
                    sums[v][la_idx, lo_idx] /= counts[la_idx, lo_idx]
                else:
                    sums[v][:, la_idx, lo_idx] /= counts[la_idx, lo_idx]
        lo_idx += 1
        if lo_idx == n_lons:
            la_idx += 1
            lo_idx = 0
    return counts, sums


def load_min_max_variance(
    file_path,
    variables,
    additional_vars,
    lons,
    delta_lon,
    lats,
    delta_lat,
    means,
    counts,
    sens_model_states=[],
    only_asc600=False,
    inoutflow_time=-1,
    min_pressure=None,
    max_pressure=None,
    pressure_levels=None,
    time_levels=None,
    verbose=False,
    process_file=None,
    relative_lon_lat=False,
):
    """

    Parameters
    ----------
    file_path
    variables
    additional_vars
    lons
    delta_lon
    lats
    delta_lat
    means
    counts
    sens_model_states
    only_asc600
    inoutflow_time
    min_pressure
    max_pressure
    pressure_levels
    time_levels
    verbose
    process_file
    relative_lon_lat

    Returns
    -------

    """
    if verbose:
        print(
            "Calculating minimum, maximum and variance of the variables in this function."
        )
    lon_name = "lon"
    lat_name = "lat"
    if relative_lon_lat:
        lon_name = "relative_lon"
        lat_name = "relative_lat"
    n_lons = len(lons) - 1
    n_lats = len(lats) - 1
    if pressure_levels is not None:
        n_press = len(pressure_levels) - 1
    else:
        n_press = 0
    if time_levels is not None:
        n_times = len(time_levels) - 1
    else:
        n_times = 1
    sens_model_state_ids = []
    if sens_model_states is not None:
        for state in sens_model_states:
            sens_model_state_ids.append(
                np.argwhere(np.asarray(param_id_map) == state).item()
            )
    # Create arrays
    if pressure_levels is not None and n_times > 1:
        n_var = n_times * n_press * n_lons * n_lats
        n_sens = len(sens_model_state_ids) * n_var
        var_shapes = (n_times, n_press, n_lats, n_lons)
        sens_shapes = (len(sens_model_state_ids), n_times, n_press, n_lats, n_lons)
    elif pressure_levels is not None:
        n_var = n_press * n_lons * n_lats
        n_sens = len(sens_model_state_ids) * n_var
        var_shapes = (n_times, n_press, n_lats, n_lons)
        sens_shapes = (len(sens_model_state_ids), n_times, n_press, n_lats, n_lons)
    elif time_levels is not None and n_times > 1:
        n_var = n_times * n_lons * n_lats
        n_sens = len(sens_model_state_ids) * n_var
        var_shapes = (n_times, n_lats, n_lons)
        sens_shapes = (len(sens_model_state_ids), n_times, n_lats, n_lons)
    else:
        n_var = n_lons * n_lats
        n_sens = len(sens_model_state_ids) * n_var
        var_shapes = (n_lats, n_lons)
        sens_shapes = (len(sens_model_state_ids), n_lats, n_lons)

    if len(sens_model_state_ids) <= 1:
        mins = {v: np.reshape(np.repeat(np.Inf, n_var), var_shapes) for v in variables}
        maxs = {v: np.reshape(np.repeat(np.NINF, n_var), var_shapes) for v in variables}
        variance = {v: np.zeros(var_shapes) for v in variables}
    else:
        variance = {}
        mins = {}
        maxs = {}
        for v in variables:
            if v[0] != "d" or v == "deposition":
                mins[v] = np.reshape(np.repeat(np.Inf, n_var), var_shapes)
                maxs[v] = np.reshape(np.repeat(np.NINF, n_var), var_shapes)
                variance[v] = np.zeros(var_shapes)
            else:
                mins[v] = np.reshape(
                    np.repeat(np.Inf, n_sens),
                    sens_shapes,
                )
                maxs[v] = np.reshape(
                    np.repeat(np.NINF, n_sens),
                    sens_shapes,
                )
                variance[v] = np.zeros(sens_shapes)
    # Get the sums for the variances and minimum and maximum values
    if verbose:
        print("Calculating the sums for variances. Getting minimum and maximum values.")
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    if process_file is not None:
        files = [files[process_file]]
    for f in tqdm(files) if verbose else files:
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[
            additional_vars + variables
        ]
        ds = filter_trajectories(
            ds=ds,
            sens_model_state_ids=sens_model_state_ids,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
        )

        lo_idx = 0
        la_idx = 0
        for la, lo in (
            tqdm(
                itertools.product(lats[:-1], lons[:-1]),
                total=n_lons * n_lats,
                leave=False,
            )
            if verbose
            else itertools.product(lats[:-1], lons[:-1])
        ):

            if pressure_levels is not None and n_times > 1:
                raise NotImplementedError(
                    "Integrating for specific time levels is not yet supported"
                )
            elif pressure_levels is not None:
                for p_i in range(n_press):
                    idx_where = (
                        (ds[lon_name] >= lo)
                        & (ds[lon_name] < lo + delta_lon)
                        & (ds[lat_name] >= la)
                        & (ds[lat_name] < la + delta_lat)
                        & (ds["pressure"] >= pressure_levels[p_i])
                        & (ds["pressure"] < pressure_levels[p_i + 1])
                    )
                    for v in variables:
                        if (
                            v[0] != "d"
                            or v == "deposition"
                            or len(sens_model_states) == 1
                        ):
                            variance[v][0, p_i, la_idx, lo_idx] += np.nansum(
                                (
                                    ds[v].where(idx_where)
                                    - means[v][0, p_i, la_idx, lo_idx]
                                )
                                ** 2
                            )
                            mins[v][0, p_i, la_idx, lo_idx] = np.nanmin(
                                [
                                    mins[v][0, p_i, la_idx, lo_idx],
                                    np.nanmin(ds[v].where(idx_where).values),
                                ]
                            )
                            maxs[v][0, p_i, la_idx, lo_idx] = np.nanmax(
                                [
                                    maxs[v][0, p_i, la_idx, lo_idx],
                                    np.nanmax(ds[v].where(idx_where).values),
                                ]
                            )
                        else:
                            sens_idx = 0
                            for id, s in zip(sens_model_state_ids, sens_model_states):
                                ds_tmp = ds.sel({"Output_Parameter_ID": id})
                                idx_where = (
                                    (ds_tmp[lon_name] >= lo)
                                    & (ds_tmp[lon_name] < lo + delta_lon)
                                    & (ds_tmp[lat_name] >= la)
                                    & (ds_tmp[lat_name] < la + delta_lat)
                                    & (ds_tmp["pressure"] >= pressure_levels[p_i])
                                    & (ds_tmp["pressure"] < pressure_levels[p_i + 1])
                                )
                                variance[v][
                                    sens_idx, 0, p_i, la_idx, lo_idx
                                ] += np.nansum(
                                    (
                                        ds_tmp[v].where(idx_where)
                                        - means[v][sens_idx, 0, p_i, la_idx, lo_idx]
                                    )
                                    ** 2
                                )
                                mins[v][sens_idx, 0, p_i, la_idx, lo_idx] = np.nanmin(
                                    [
                                        mins[v][sens_idx, 0, p_i, la_idx, lo_idx],
                                        np.nanmin(ds_tmp[v].where(idx_where).values),
                                    ]
                                )
                                maxs[v][sens_idx, 0, p_i, la_idx, lo_idx] = np.nanmax(
                                    [
                                        maxs[v][sens_idx, 0, p_i, la_idx, lo_idx],
                                        np.nanmax(ds_tmp[v].where(idx_where).values),
                                    ]
                                )
                                sens_idx += 1
            elif time_levels is not None and n_times > 1:
                raise NotImplementedError(
                    "Integrating for specific time levels is not yet supported"
                )
            else:
                idx_where = (
                    (ds[lon_name] >= lo)
                    & (ds[lon_name] < lo + delta_lon)
                    & (ds[lat_name] >= la)
                    & (ds[lat_name] < la + delta_lat)
                )
                for v in variables:
                    if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
                        variance[v][la_idx, lo_idx] += np.nansum(
                            (ds[v].where(idx_where) - means[v][la_idx, lo_idx]) ** 2
                        )
                        mins[v][la_idx, lo_idx] = np.nanmin(
                            [
                                mins[v][la_idx, lo_idx],
                                np.nanmin(ds[v].where(idx_where).values),
                            ]
                        )
                        maxs[v][la_idx, lo_idx] = np.nanmax(
                            [
                                maxs[v][la_idx, lo_idx],
                                np.nanmax(ds[v].where(idx_where).values),
                            ]
                        )
                    else:
                        sens_idx = 0
                        for id, s in zip(sens_model_state_ids, sens_model_states):
                            ds_tmp = ds.sel({"Output_Parameter_ID": id})
                            idx_where = (
                                (ds_tmp[lon_name] >= lo)
                                & (ds_tmp[lon_name] < lo + delta_lon)
                                & (ds_tmp[lat_name] >= la)
                                & (ds_tmp[lat_name] < la + delta_lat)
                            )
                            variance[v][sens_idx, la_idx, lo_idx] += np.nansum(
                                (
                                    ds_tmp[v].where(idx_where)
                                    - means[v][sens_idx, la_idx, lo_idx]
                                )
                                ** 2
                            )
                            mins[v][sens_idx, la_idx, lo_idx] = np.nanmin(
                                [
                                    mins[v][sens_idx, la_idx, lo_idx],
                                    np.nanmin(ds_tmp[v].where(idx_where).values),
                                ]
                            )
                            maxs[v][sens_idx, la_idx, lo_idx] = np.nanmax(
                                [
                                    maxs[v][sens_idx, la_idx, lo_idx],
                                    np.nanmax(ds_tmp[v].where(idx_where).values),
                                ]
                            )
                            sens_idx += 1
            lo_idx += 1
            if lo_idx == n_lons:
                la_idx += 1
                lo_idx = 0
    if verbose:
        print("Divide the sum to get the expected value (here: the variance)")
    lo_idx = 0
    la_idx = 0
    for _, _ in (
        tqdm(itertools.product(lons[:-1], lats[:-1]), total=n_lons * n_lats)
        if verbose
        else itertools.product(lons[:-1], lats[:-1])
    ):
        for v in variance:
            if pressure_levels is not None and n_times > 1:
                raise NotImplementedError(
                    "Integrating for specific time levels is not yet supported"
                )
            elif pressure_levels is not None:
                for p_i in range(n_press):
                    if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
                        variance[v][0, p_i, la_idx, lo_idx] /= counts[
                            0, p_i, la_idx, lo_idx
                        ]
                    else:
                        variance[v][:, 0, p_i, la_idx, lo_idx] /= counts[
                            0, p_i, la_idx, lo_idx
                        ]
            elif time_levels is not None and n_times > 1:
                raise NotImplementedError(
                    "Integrating for specific time levels is not yet supported"
                )
            else:
                if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
                    variance[v][la_idx, lo_idx] /= counts[la_idx, lo_idx]
                else:
                    variance[v][:, la_idx, lo_idx] /= counts[la_idx, lo_idx]
        lo_idx += 1
        if lo_idx == n_lons:
            la_idx += 1
            lo_idx = 0
    return mins, maxs, variance


def to_Dataset(
    file_path,
    lons,
    lats,
    pressure_levels,
    time_levels,
    delta_time,
    sens_model_states,
    counts,
    means,
    mins,
    maxs,
    variance,
    relative_lon_lat,
    verbose,
):
    """

    Parameters
    ----------
    file_path
    lons
    lats
    pressure_levels
    time_levels
    delta_time
    sens_model_states
    counts
    means
    mins
    maxs
    variance
    verbose

    Returns
    -------

    """
    if verbose:
        print("Setting up dimensions and coordinates for xarray.Dataset")
    dataset_dic = {}
    sens_values = False
    if pressure_levels is not None and delta_time is not None:
        raise NotImplementedError(
            "Integrating for specific time levels is not yet supported"
        )
    elif pressure_levels is not None:
        var_dims = ["time", "pressure", "lon", "lat"]
        sens_dims = ["Output Parameter", "time", "pressure", "lon", "lat"]
        coords = {
            "lon": lons[:-1],
            "lat": lats[:-1],
            "pressure": pressure_levels[:-1],
            "time": [time_levels[-1]],
        }
    elif delta_time is not None:
        raise NotImplementedError(
            "Integrating for specific time levels is not yet supported"
        )
    else:
        var_dims = ["lon", "lat"]
        sens_dims = ["Output Parameter", "lon", "lat"]
        coords = {"lon": lons[:-1], "lat": lats[:-1]}
    for v in means:
        if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
            dataset_dic[f"Mean {v}"] = (var_dims, means[v])
        else:
            dataset_dic[f"Mean {v}"] = (
                sens_dims,
                means[v],
            )
            sens_values = True
    for v in variance:
        if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
            dataset_dic[f"Var {v}"] = (var_dims, variance[v])
        else:
            dataset_dic[f"Var {v}"] = (
                sens_dims,
                variance[v],
            )
            sens_values = True
    for v in mins:
        if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
            dataset_dic[f"Min {v}"] = (var_dims, mins[v])
        else:
            dataset_dic[f"Min {v}"] = (
                sens_dims,
                mins[v],
            )
            sens_values = True
    for v in maxs:
        if v[0] != "d" or v == "deposition" or len(sens_model_states) == 1:
            dataset_dic[f"Max {v}"] = (var_dims, maxs[v])
        else:
            dataset_dic[f"Max {v}"] = (
                sens_dims,
                maxs[v],
            )
            sens_values = True
    dataset_dic["counts"] = (var_dims, counts)

    if sens_values:
        coords["Output Parameter"] = sens_model_states
    ds = xr.Dataset(
        data_vars=dataset_dic,
        coords=coords,
    )
    # set attributes and set NaN where possible
    if "pressure" in ds and not ds["pressure"].dtype == np.float64:
        ds["pressure"] = ds["pressure"].astype(np.float64)
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds_tmp = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    for key in ds:
        if key in ds_tmp:
            ds[key].attrs = ds_tmp[key].attrs
    if relative_lon_lat:
        ds["lon"].attrs = ds_tmp["relative_lon"].attrs
        ds["lat"].attrs = ds_tmp["relative_lat"].attrs
    ds["counts"].attrs = {
        "long_name": "Count",
        "standard_name": "count",
        "auxiliary_data": "yes",
        "_FillValue": np.NaN,
    }
    ds["counts"] = ds["counts"].where(ds["counts"] > 0)
    for v in mins:
        ds[f"Min {v}"] = ds[f"Min {v}"].where(ds["counts"] > 0)
        ds[f"Max {v}"] = ds[f"Max {v}"].where(ds["counts"] > 0)
        ds[f"Var {v}"] = ds[f"Var {v}"].where(ds["counts"] > 0)
        ds[f"Mean {v}"] = ds[f"Mean {v}"].where(ds["counts"] > 0)
        ds[f"Min {v}"].attrs = ds_tmp[v].attrs
        ds[f"Max {v}"].attrs = ds_tmp[v].attrs
        ds[f"Var {v}"].attrs = ds_tmp[v].attrs
        ds[f"Mean {v}"].attrs = ds_tmp[v].attrs

    # Zero values in latitude and longitudes is considered garbage values, so we delete that.
    ds["lat"] = ds["lat"].where(ds["lat"] != 0)
    ds["lon"] = ds["lon"].where(ds["lon"] != 0)
    ds["lon"].attrs = ds_tmp["lon"].attrs
    ds["lat"].attrs = ds_tmp["lat"].attrs
    if "pressure" in ds:
        ds["pressure"].attrs = ds_tmp["pressure"].attrs
    if "time" in ds:
        ds["time"].attrs = ds_tmp["time"].attrs
    if "Output Parameter" in ds:
        ds["Output Parameter"].attrs = {
            "long_name": "gradients are calculated w.r.t. this output parameter",
            "standard_name": "Output Parameter",
            "auxiliary_data": "yes",
        }
    return ds


def load_data(
    file_path,
    variables,
    sens_model_states=[],
    only_asc600=False,
    inoutflow_time=-1,
    min_pressure=None,
    max_pressure=None,
    pressure_levels=None,
    delta_time=None,
    n_lons=100,
    n_lats=100,
    relative_lon_lat=False,
    verbose=False,
):
    """

    Parameters
    ----------
    file_path
    variables
    sens_model_states
    only_asc600
    inoutflow_time
    min_pressure
    max_pressure
    pressure_levels
    delta_time
    n_lons
    n_lats
    verbose

    Returns
    -------

    """
    lons, delta_lon, lats, delta_lat, time_levels = load_lon_lat_time(
        file_path=file_path,
        sens_model_states=sens_model_states,
        only_asc600=only_asc600,
        inoutflow_time=inoutflow_time,
        min_pressure=min_pressure,
        max_pressure=max_pressure,
        n_lons=n_lons,
        n_lats=n_lats,
        delta_time=delta_time,
        relative_lon_lat=relative_lon_lat,
        verbose=verbose,
    )
    if relative_lon_lat:
        additional_vars = [
            "lon",
            "lat",
            "relative_lon",
            "relative_lat",
            "asc600",
            "phase",
        ]
    else:
        additional_vars = ["lon", "lat", "asc600", "phase"]
    if "pressure" not in variables:
        additional_vars.append("pressure")
    counts, means = load_counts_means(
        file_path=file_path,
        variables=variables,
        additional_vars=additional_vars,
        lons=lons,
        delta_lon=delta_lon,
        lats=lats,
        delta_lat=delta_lat,
        sens_model_states=sens_model_states,
        only_asc600=only_asc600,
        inoutflow_time=inoutflow_time,
        min_pressure=min_pressure,
        max_pressure=max_pressure,
        pressure_levels=pressure_levels,
        time_levels=time_levels,
        relative_lon_lat=relative_lon_lat,
        verbose=verbose,
    )
    mins, maxs, variance = load_min_max_variance(
        file_path=file_path,
        variables=variables,
        additional_vars=additional_vars,
        lons=lons,
        delta_lon=delta_lon,
        lats=lats,
        delta_lat=delta_lat,
        means=means,
        counts=counts,
        sens_model_states=sens_model_states,
        only_asc600=only_asc600,
        inoutflow_time=inoutflow_time,
        min_pressure=min_pressure,
        max_pressure=max_pressure,
        pressure_levels=pressure_levels,
        time_levels=time_levels,
        relative_lon_lat=relative_lon_lat,
        verbose=verbose,
    )
    return to_Dataset(
        file_path=file_path,
        lons=lons,
        lats=lats,
        pressure_levels=pressure_levels,
        time_levels=time_levels,
        delta_time=delta_time,
        sens_model_states=sens_model_states,
        counts=counts,
        means=means,
        mins=mins,
        maxs=maxs,
        variance=variance,
        relative_lon_lat=relative_lon_lat,
        verbose=verbose,
    )


def plot_heatmap(ds, col, fig_size=250, aspect=1, cmap="viridis"):
    """
    An easy and fast way to plot a given column along latitude and
    longitude.

    Parameters
    ----------
    ds
    col
    fig_size
    aspect
    cmap

    Returns
    -------
    holoviews.Image
    """
    x = "longitude"
    y = "latitude"
    return hv.Image((ds[x], ds[y], ds[col]), datatype=["grid"]).opts(
        fig_size=fig_size,
        aspect=aspect,
        cmap=cmap,
        colorbar=True,
        xlabel=x,
        ylabel=y,
        clabel=col,
    )


def plot_2dmap_interactive(ds):
    """
    Calling this function from a Jupyter notebook allows to visualize the
    2d maps at different heights for all available columns in ds.

    Parameters
    ----------
    ds : xarray.Dataset

    Returns
    -------

    """
    sns.set_style("darkgrid")

    pressure = pn.widgets.FloatSlider(
        name="pressure",
        start=ds["pressure"].min().values.item(),
        end=ds["pressure"].max().values.item(),
        step=(ds["pressure"].values[1] - ds["pressure"].values[0]),
        value=ds["pressure"].values[-4],
    )
    out_param = pn.widgets.RadioButtonGroup(
        name="Output Parameter",
        value=ds["Output Parameter"].values[0],
        options=list(ds["Output Parameter"].values),
        button_type="primary",
    )
    kind_param = pn.widgets.RadioButtonGroup(
        name="Kind",
        value="Mean",
        options=["Mean", "Min", "Max", "Var"],
        button_type="primary",
    )
    in_params = []
    for col in ds:
        if "Mean " in col:
            in_params.append(col[5:])
    if "counts" in ds:
        in_params.append("counts")
    in_params.sort()
    in_param = pn.widgets.Select(
        name="Color according to",
        value=in_params[-1],
        options=in_params,
    )
    color_map = pn.widgets.Select(
        name="Colormap",
        value="RdBu",
        options=["RdBu", "viridis", "blues"],
    )
    fix = pn.widgets.Toggle(
        name="Fix colorbar over all levels",
        button_type="success",
    )
    static = pn.widgets.StaticText(name="Min, Max Values", value="")
    width_slider = pn.widgets.IntSlider(
        name="Width in inches",
        start=3,
        end=15,
        step=1,
        value=9,
    )
    height_slider = pn.widgets.IntSlider(
        name="Height in inches",
        start=3,
        end=15,
        step=1,
        value=6,
    )
    log_plot = pn.widgets.Toggle(
        name="Use log colorbar",
        value=True,
        button_type="success",
    )
    log_threshold_slider = pn.widgets.FloatSlider(
        name="Log Threshold (set to zero for automatic estimation)",
        start=-25,
        end=0,
        value=0,
        step=1,
    )
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=2,
        step=0.1,
        value=0.7,
    )

    def plot_me(
        o_p,
        i_p,
        k_p,
        p,
        c,
        fix,
        log_plot,
        width,
        height,
        lthresh,
        title,
        font_scale,
        save,
        save_path,
    ):
        if title == "":
            title = None
        if i_p == "counts":
            ds_tmp = ds.isel({"time": 0})[i_p]
        elif i_p[0] == "d" and i_p != "deposition":
            ds_tmp = ds.isel({"time": 0}).sel({"Output Parameter": o_p})[f"{k_p} {i_p}"]
        else:
            ds_tmp = ds.isel({"time": 0})[f"{k_p} {i_p}"]
        mini = ds_tmp.min().values.item()
        maxi = ds_tmp.max().values.item()
        if lthresh == 0:
            linthresh = np.nanmin(np.abs(ds_tmp.where(np.abs(ds_tmp) > 0)))
        else:
            linthresh = 10 ** lthresh
        if i_p == "counts":
            ds_tmp = ds[i_p].isel({"time": 0}).sel({"pressure": p})
        elif i_p[0] == "d" and i_p != "deposition":
            ds_tmp = (
                ds[f"{k_p} {i_p}"]
                .isel({"time": 0})
                .sel({"Output Parameter": o_p, "pressure": p})
            )
        else:
            ds_tmp = ds[f"{k_p} {i_p}"].isel({"time": 0}).sel({"pressure": p})
        min_local = np.nanmin(ds_tmp)
        max_local = np.nanmax(ds_tmp)
        static.value = f"({mini:.2e}, {maxi:.2e}); at {p/100} hPa: ({min_local:.2e}, {max_local:.2e})"
        fig = Figure(
            figsize=(width, height),
        )
        ax = fig.subplots()
        if fix:
            if np.abs(mini) > maxi:
                maxi = np.abs(mini)
            else:
                mini = -maxi
            if log_plot:
                ds_tmp.plot(
                    x="lat",
                    y="lon",
                    cmap=c,
                    norm=SymLogNorm(
                        linthresh=linthresh,
                        vmin=mini,
                        vmax=maxi,
                        base=10,
                    ),
                    ax=ax,
                )
            else:
                ds_tmp.plot(
                    x="lat",
                    y="lon",
                    cmap=c,
                    vmin=mini,
                    vmax=maxi,
                    ax=ax,
                )
            _ = ax.set_title(title, fontsize=int(12 * font_scale))
            ax.tick_params(
                axis="both",
                which="major",
                labelsize=int(10 * font_scale),
            )
            ax.xaxis.get_label().set_fontsize(int(11 * font_scale))
            ax.yaxis.get_label().set_fontsize(int(11 * font_scale))
            ax.yaxis.grid(True, which="major")
            ax.xaxis.grid(True, which="major")
            cbar = ax.collections[-1].colorbar
            cbarax = cbar.ax
            cbarax.tick_params(labelsize=int(10 * font_scale))
            cbar.set_label(
                label=f"{k_p} {i_p}",
                fontsize=int(11 * font_scale),
            )
        else:
            if log_plot:
                if lthresh == 0:
                    linthresh = np.nanmin(np.abs(ds_tmp.where(np.abs(ds_tmp) > 0)))
                else:
                    linthresh = 10 ** lthresh

                ds_tmp.plot(
                    x="lat",
                    y="lon",
                    cmap=c,
                    norm=SymLogNorm(
                        linthresh=linthresh,
                        base=10,
                    ),
                    ax=ax,
                )
            else:
                ds_tmp.plot(
                    x="lat",
                    y="lon",
                    cmap=c,
                    ax=ax,
                )
            _ = ax.set_title(title, fontsize=int(12 * font_scale))
            ax.tick_params(
                axis="both",
                which="major",
                labelsize=int(10 * font_scale),
            )
            ax.xaxis.get_label().set_fontsize(int(11 * font_scale))
            ax.yaxis.get_label().set_fontsize(int(11 * font_scale))
            ax.yaxis.grid(True, which="major")
            ax.xaxis.grid(True, which="major")
            cbar = ax.collections[-1].colorbar
            cbarax = cbar.ax
            cbarax.tick_params(labelsize=int(10 * font_scale))
            cbar.set_label(
                label=f"{k_p} {i_p}",
                fontsize=int(11 * font_scale),
            )
        if save:
            try:
                ax.figure.savefig(save_path, bbox_inches="tight", dpi=300)
            except:
                save_to_field.value = (
                    f"Could not save to {save_path}. Did you forget the filetype?"
                )
                pass
            save_path = None
            save = False
        return fig

    plot_pane = pn.panel(
        pn.bind(
            plot_me,
            o_p=out_param,
            i_p=in_param,
            k_p=kind_param,
            p=pressure,
            c=color_map,
            fix=fix,
            log_plot=log_plot,
            width=width_slider,
            height=height_slider,
            lthresh=log_threshold_slider,
            title=title_widget,
            font_scale=font_slider,
            save=save_button,
            save_path=save_to_field,
        ),
    ).servable()

    return pn.Column(
        "# Plot Grid",
        static,
        pressure,
        kind_param,
        out_param,
        pn.Row(
            in_param,
            color_map,
        ),
        pn.Row(
            fix,
            log_plot,
        ),
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
        ),
        log_threshold_slider,
        pn.Row(
            save_to_field,
            save_button,
        ),
        title_widget,
        plot_pane,
    )


if __name__ == "__main__":
    import argparse
    import pickle
    import textwrap

    from latexify import *

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Load data from a sensitivity simulation and calculate column-wise mean, maximum, variance of
            any given parameter or sensitivity.
            Can optionally set minimum and maximum height via --min_pressure or 
            --max_pressure to exclude certain heights.
            Can include a density map of trajectories, phases and filter by phases. 
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file_path",
        default="../data/vladiana_ensembles/",
        help=textwrap.dedent(
            """\
            Path to a folder with many files from a sensitivity analysis simulation.
            """
        ),
    )
    parser.add_argument(
        "--out_file",
        default="../pics/plots.png",
        help=textwrap.dedent(
            """\
            Path and name to store plots.
            """
        ),
    )
    parser.add_argument(
        "--width",
        default=24,
        type=float,
        help=textwrap.dedent(
            """\
            Width in inches for histogram plots.
            """
        ),
    )
    parser.add_argument(
        "--height",
        default=12,
        type=float,
        help=textwrap.dedent(
            """\
            Height in inches for histogram plots.
            """
        ),
    )
    parser.add_argument(
        "--delta_time",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            Integrate over the given time interval in seconds. If None given, then integrate over all time steps.
            Is not implemented yet
            """
        ),
    )
    parser.add_argument(
        "--n_lons",
        default=100,
        type=int,
        help=textwrap.dedent(
            """\
            The number of bins along longitude.
            """
        ),
    )
    parser.add_argument(
        "--n_lats",
        default=100,
        type=int,
        help=textwrap.dedent(
            """\
            The number of bins along latitude.
            """
        ),
    )
    parser.add_argument(
        "--pressure_levels",
        default=None,
        nargs="+",
        type=float,
        help=textwrap.dedent(
            """\
            A list of pressure levels in Pa to use. Supersedes --min_pressure and --max_pressure
            """
        ),
    )
    parser.add_argument(
        "--only_asc600",
        action="store_true",
        help=textwrap.dedent(
            """\
            Consider only time steps during the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--relative_lon_lat",
        action="store_true",
        help=textwrap.dedent(
            """\
            Use longitude and latitude degrees relative to the start of the ascent.
            """
        ),
    )
    parser.add_argument(
        "--inoutflow_time",
        default=-1,
        type=int,
        help=textwrap.dedent(
            """\
            Consider only time steps during the fastest ascent and within the given range before (inflow) 
            and after (outflow) of the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--select_phase",
        default=None,
        type=str,
        help=textwrap.dedent(
            """\
            Consider only the given phase. Options are:
            'warm', 'mixed', 'ice', 'all'. The latter creates a different plot for each phase.
            Not yet implemented.
            """
        ),
    )
    parser.add_argument(
        "--plot_type",
        default=["density"],
        type=str,
        nargs="+",
        choices=["density", "mean", "min", "max", "sd"],
        help=textwrap.dedent(
            """\
            Define which kind of plot shall be created. Multiple options can be selected to create multiple plots.
            Options are: 
            'density': Plot the amount of datapoints in each grid point as heatmap. 
            'mean': Plot the mean in each grid point as a heatmap.
            'min': Plot the minimum value in each grid point as a heatmap.
            'max': Plot the maximum value in each grid point as a heatmap.
            'sd': Plot the standard deviation in each grid point as a heatmap.
            """
        ),
    )
    parser.add_argument(
        "--var",
        default=None,
        type=str,
        nargs="+",
        help=textwrap.dedent(
            """\
            Plot the given variable on a map. If the variable is a sensitivity, then set --sens_model_states 
            to a model state variable for which the sensitivity is for. 
            Otherwise a sensitivity plot for each model state variable will be made.
            May include multiple variables which generates plots for each of them.
            
            """
        ),
    )
    parser.add_argument(
        "--sens_model_states",
        default=None,
        type=str,
        nargs="+",
        help=textwrap.dedent(
            """\
            If --var is a sensitivity (model parameter), then you may define the model state variable here 
            to calculate and plot only sensitivities regarding this model state variable.
            """
        ),
    )
    parser.add_argument(
        "--min_pressure",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            Filter such that only values at this pressure in Pa or higher (= at this height or lower) are considered.
            """
        ),
    )
    parser.add_argument(
        "--max_pressure",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            Filter such that only values at this pressure in Pa or lower (= at this height or higher) are considered.
            """
        ),
    )
    parser.add_argument(
        "--store_calculated",
        default=None,
        type=str,
        help=textwrap.dedent(
            """\
            Store the processed data such that it can be loaded and plotted faster again or to continue calculation. 
            """
        ),
    )
    parser.add_argument(
        "--load_calculated",
        default=None,
        type=str,
        help=textwrap.dedent(
            """\
            Path to load previously processed data.
            """
        ),
    )
    parser.add_argument(
        "--calculate_only",
        default="all",
        type=str,
        choices=[
            "all",
            "only_lon_lat_time_levels",
            "only_counts_means",
            "only_min_max_variance",
            "merge_counts_means",
            "merge_min_max_variance",
        ],
        help=textwrap.dedent(
            """\
            Calculate only the given step. Previous steps must be loaded using --load_calculated!
            All options except 'all' and 'merge_counts_means' store only *.pkl files for further processing.
            The 'merge_*' options merge all results calculated 
            using the respective steps for different files with --process_file. 
            """
        ),
    )
    parser.add_argument(
        "--process_file",
        default=None,
        type=int,
        help=textwrap.dedent(
            """\
            If multiple files in --file_path are present, you can define which file shall be processed using an integer 
            from 0  to n-1, whith n the number of files. This is useful if the script is run in parallel 
            for multiple files. Makes only sense with --calculate_only 
            and the options 'only_counts_means' or 'only_min_max_variance'.
            """
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=textwrap.dedent(
            """\
            More output regarding the calculations.
            """
        ),
    )
    args = parser.parse_args()

    if args.calculate_only == "only_lon_lat_time_levels":
        if args.pressure_levels is None:
            min_pressure = args.min_pressure
            max_pressure = args.max_pressure
        else:
            min_pressure = args.pressure_levels[0]
            max_pressure = args.pressure_levels[-1]
        lons, delta_lon, lats, delta_lat, time_levels = load_lon_lat_time(
            file_path=args.file_path,
            sens_model_states=args.sens_model_states,
            only_asc600=args.only_asc600,
            inoutflow_time=args.inoutflow_time,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
            n_lons=args.n_lons,
            n_lats=args.n_lats,
            delta_time=args.delta_time,
            relative_lon_lat=args.relative_lon_lat,
            verbose=args.verbose,
        )
        data = {
            "lons": lons,
            "delta_lon": delta_lon,
            "n_lons": args.n_lons,
            "lats": lats,
            "delta_lat": delta_lat,
            "n_lats": args.n_lats,
            "time_levels": time_levels,
            "delta_time": args.delta_time,
            "only_asc600": args.only_asc600,
            "sens_model_states": args.sens_model_states,
            "min_pressure": min_pressure,
            "max_pressure": max_pressure,
            "inoutflow_time": args.inoutflow_time,
            "relative_lon_lat": args.relative_lon_lat,
        }
        if args.verbose:
            print("########### Store longitudes, latitudes and time levels ###########")
        with open(args.store_calculated + "lon_lat_time.pkl", "wb") as f:
            pickle.dump(data, f)
    elif args.calculate_only == "only_counts_means":
        if args.load_calculated is None:
            raise ValueError(
                "Can not use '--calculate_only' for counts and means without a path to lon_lat_time.pkl."
            )
        with open(args.load_calculated + "lon_lat_time.pkl", "rb") as f:
            data = pickle.load(f)
            lons = data["lons"]
            delta_lon = data["delta_lon"]
            n_lons = data["n_lons"]
            lats = data["lats"]
            delta_lat = data["delta_lat"]
            n_lats = data["n_lats"]
            time_levels = data["time_levels"]
            delta_time = data["delta_time"]
            only_asc600 = data["only_asc600"]
            sens_model_states = data["sens_model_states"]
            min_pressure = data["min_pressure"]
            max_pressure = data["max_pressure"]
            inoutflow_time = data["inoutflow_time"]
            relative_lon_lat = data["relative_lon_lat"]
        if relative_lon_lat:
            additional_vars = [
                "lon",
                "lat",
                "relative_lon",
                "relative_lat",
                "asc600",
                "phase",
            ]
        else:
            additional_vars = ["lon", "lat", "asc600", "phase"]

        if "pressure" not in args.var:
            additional_vars.append("pressure")
        counts, means = load_counts_means(
            file_path=args.file_path,
            variables=args.var,
            additional_vars=additional_vars,
            lons=lons,
            delta_lon=delta_lon,
            lats=lats,
            delta_lat=delta_lat,
            sens_model_states=sens_model_states,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
            pressure_levels=args.pressure_levels,
            time_levels=time_levels,
            verbose=args.verbose,
            process_file=args.process_file,
            relative_lon_lat=relative_lon_lat,
        )
        data = {
            "counts": counts,
            "means": means,
            "pressure_levels": args.pressure_levels,
            "variables": args.var,
        }
        if args.verbose:
            print("########### Store counts and means ###########")
        if args.process_file is not None:
            filename = f"counts_means_{args.process_file}.pkl"
        else:
            filename = "counts_means.pkl"
        with open(args.store_calculated + filename, "wb") as f:
            pickle.dump(data, f)
    elif args.calculate_only == "merge_counts_means":
        # There are two ways to do that.
        # a: all files are in a single folder and have counts and means
        # for the same variables
        # b: all files are distributed in different folders where each folder
        # is the name of a variable dimension and within each folder
        # are multiple files just as in a
        files = [
            f
            for f in os.listdir(args.load_calculated)
            if os.path.isfile(args.load_calculated + f)
            and f.startswith("counts_means_")
        ]

        def do_the_merge(files, filepath, leave=False):
            files = np.sort(files)
            means = None
            counts = None
            variables = None
            pressure_levels = None
            for filename in tqdm(files, leave=leave) if args.verbose else files:
                with open(filepath + filename, "rb") as f:
                    data = pickle.load(f)
                    if counts is not None:
                        old_count = counts.copy()
                        loaded_counts = data["counts"]
                        loaded_means = data["means"]
                        loaded_counts[np.isnan(loaded_counts)] = 0
                        counts += loaded_counts
                        for key in means:
                            loaded_means[key][np.isnan(loaded_means[key])] = 0
                            means[key] = (
                                means[key] * old_count / counts
                                + loaded_means[key] * loaded_counts / counts
                            )
                    else:
                        counts = data["counts"]
                        means = data["means"]
                        variables = data["variables"]
                        pressure_levels = data["pressure_levels"]
                    counts[np.isnan(counts)] = 0
                    for key in means:
                        means[key][np.isnan(means[key])] = 0
            return means, variables, counts, pressure_levels

        if len(files) == 0:
            # This is approach b with multiple folders
            folders = [
                f
                for f in os.listdir(args.load_calculated)
                if os.path.isdir(args.load_calculated + f)
            ]
            folders = np.sort(folders)
            counts = None
            pressure_levels = None
            variables = []
            means = None
            for dir in tqdm(folders) if args.verbose else folders:
                filepath = args.load_calculated + dir + "/"
                files = [
                    f
                    for f in os.listdir(filepath)
                    if os.path.isfile(filepath + f) and f.startswith("counts_means_")
                ]
                if counts is None:
                    means, tmp_variables, counts, pressure_levels = do_the_merge(
                        files, filepath
                    )
                else:
                    tmp_means, tmp_variables, _, _ = do_the_merge(files, filepath)
                    for key in tmp_means:
                        means[key] = tmp_means[key]
                variables.extend(tmp_variables)
        else:
            means, variables, counts, pressure_levels = do_the_merge(
                files, args.load_calculated, True
            )

        if args.verbose:
            print("########### Store counts and means ###########")
        data = {
            "counts": counts,
            "means": means,
            "pressure_levels": pressure_levels,
            "variables": variables,
        }
        with open(args.store_calculated + "counts_means.pkl", "wb") as f:
            pickle.dump(data, f)
    elif args.calculate_only == "only_min_max_variance":
        if args.load_calculated is None:
            raise ValueError(
                "Can not use '--calculate_only' for counts and means "
                + "without a path to lon_lat_time.pkl and without a path to counts_means.pkl."
            )
        with open(args.load_calculated + "lon_lat_time.pkl", "rb") as f:
            data = pickle.load(f)
            lons = data["lons"]
            delta_lon = data["delta_lon"]
            n_lons = data["n_lons"]
            lats = data["lats"]
            delta_lat = data["delta_lat"]
            n_lats = data["n_lats"]
            time_levels = data["time_levels"]
            delta_time = data["delta_time"]
            only_asc600 = data["only_asc600"]
            sens_model_states = data["sens_model_states"]
            min_pressure = data["min_pressure"]
            max_pressure = data["max_pressure"]
            inoutflow_time = data["inoutflow_time"]
            relative_lon_lat = data["relative_lon_lat"]

        with open(args.load_calculated + "counts_means.pkl", "rb") as f:
            data = pickle.load(f)
            counts = data["counts"]
            means = data["means"]
            pressure_levels = data["pressure_levels"]
            if isinstance(args.var, str):
                if args.var not in data["variables"]:
                    raise ValueError(
                        f"You asked for variable {args.var} which is not "
                        + f"present in {args.load_calculated}counts_means.pkl"
                    )
                variables = [args.var]
            elif args.var is None:
                variables = data["variables"]
            else:
                variables = args.var
                for var in variables:
                    if var not in data["variables"]:
                        raise ValueError(
                            f"You asked for variable {var} which is not "
                            + f"present in {args.load_calculated}counts_means.pkl"
                        )
        if relative_lon_lat:
            additional_vars = [
                "lon",
                "lat",
                "relative_lon",
                "relative_lat",
                "asc600",
                "phase",
            ]
        else:
            additional_vars = ["lon", "lat", "asc600", "phase"]
        if "pressure" not in variables:
            additional_vars.append("pressure")
        mins, maxs, variance = load_min_max_variance(
            file_path=args.file_path,
            variables=variables,
            additional_vars=additional_vars,
            lons=lons,
            delta_lon=delta_lon,
            lats=lats,
            delta_lat=delta_lat,
            means=means,
            counts=counts,
            sens_model_states=sens_model_states,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
            pressure_levels=pressure_levels,
            time_levels=time_levels,
            verbose=args.verbose,
            process_file=args.process_file,
            relative_lon_lat=relative_lon_lat,
        )
        if args.verbose:
            print("########### Store minimums, maximums and variances ###########")
        data = {
            "mins": mins,
            "maxs": maxs,
            "variance": variance,
        }
        if args.process_file is not None:
            filename = f"min_max_variances_{args.process_file}.pkl"
        else:
            filename = "min_max_variances.pkl"
        with open(args.store_calculated + filename, "wb") as f:
            pickle.dump(data, f)
    elif args.calculate_only == "merge_min_max_variance":
        # There are two ways to do that.
        # a: all files are in a single folder and have counts and means
        # for the same variables
        # b: all files are distributed in different folders where each folder
        # is the name of a variable dimension and within each folder
        # are multiple files just as in a
        files = [
            f
            for f in os.listdir(args.load_calculated)
            if os.path.isfile(args.load_calculated + f)
            and f.startswith("min_max_variances_")
        ]

        def do_the_merge(files, filepath, leave=False):
            files = np.sort(files)
            mins = None
            maxs = None
            variance = None
            for filename in tqdm(files, leave=leave) if args.verbose else files:
                with open(filepath + filename, "rb") as f:
                    data = pickle.load(f)
                    if mins is not None:
                        for key in mins:
                            loaded_var = data["variance"][key]
                            loaded_var[np.isnan(loaded_var)] = 0
                            mins[key] = np.fmin(mins[key], data["mins"][key])
                            maxs[key] = np.fmax(maxs[key], data["maxs"][key])
                            variance[key] += loaded_var
                    else:
                        mins = data["mins"]
                        maxs = data["maxs"]
                        variance = data["variance"]
                        for key in variance:
                            variance[key][np.isnan(variance[key])] = 0
            return mins, maxs, variance

        if len(files) == 0:
            # This is approach b with multiple folders
            folders = [
                f
                for f in os.listdir(args.load_calculated)
                if os.path.isdir(args.load_calculated + f)
            ]
            folders = np.sort(folders)
            mins = None
            maxs = None
            variance = None
            for dir in tqdm(folders) if args.verbose else folders:
                filepath = args.load_calculated + dir + "/"
                files = [
                    f
                    for f in os.listdir(filepath)
                    if os.path.isfile(filepath + f)
                    and f.startswith("min_max_variances_")
                ]
                if mins is None:
                    mins, maxs, variance = do_the_merge(files, filepath)
                else:
                    tmp_mins, tmp_maxs, tmp_variance = do_the_merge(files, filepath)
                    for key in tmp_mins:
                        mins[key] = tmp_mins[key]
                        maxs[key] = tmp_maxs[key]
                        variance[key] = tmp_variance[key]
        else:
            mins, maxs, variance = do_the_merge(files, args.load_calculated, True)
        # During the merge, any NaNs are replaced with for the variance.
        # We need to plug in the NaNs where no data is available.
        with open(args.load_calculated + "counts_means.pkl", "rb") as f:
            data_cm = pickle.load(f)
        for key in variance:
            if key[0] != "d" or key == "deposition":
                variance[key][data_cm["counts"] == 0] = np.NaN
                # While not *necessary* it is consistent to do that for min and max too
                mins[key][data_cm["counts"] == 0] = np.NaN
                maxs[key][data_cm["counts"] == 0] = np.NaN
            else:
                variance[key][:, data_cm["counts"] == 0] = np.NaN
                # While not *necessary* it is consistent to do that for min and max too
                mins[key][:, data_cm["counts"] == 0] = np.NaN
                maxs[key][:, data_cm["counts"] == 0] = np.NaN
        if args.verbose:
            print("########### Store minimums, maximums and variances ###########")
        data = {
            "mins": mins,
            "maxs": maxs,
            "variance": variance,
        }
        with open(args.store_calculated + "min_max_variances.pkl", "wb") as f:
            pickle.dump(data, f)
        if args.file_path is not None:
            pressure_levels = data_cm["pressure_levels"]
            counts = data_cm["counts"]
            means = data_cm["means"]
            with open(args.load_calculated + "lon_lat_time.pkl", "rb") as f:
                data = pickle.load(f)
                lons = data["lons"]
                n_lons = data["n_lons"]
                lats = data["lats"]
                n_lats = data["n_lats"]
                time_levels = data["time_levels"]
                delta_time = data["delta_time"]
                sens_model_states = data["sens_model_states"]
                relative_lon_lat = data["relative_lon_lat"]

            ds = to_Dataset(
                file_path=args.file_path,
                lons=lons,
                lats=lats,
                pressure_levels=pressure_levels,
                time_levels=time_levels,
                delta_time=delta_time,
                sens_model_states=sens_model_states,
                counts=counts,
                means=means,
                mins=mins,
                maxs=maxs,
                variance=variance,
                relative_lon_lat=relative_lon_lat,
                verbose=args.verbose,
            )
            comp = dict(zlib=True, complevel=9)
            encoding = {var: comp for var in ds.data_vars}
            ds.to_netcdf(
                path=f"{args.store_calculated}grid_{n_lons}_{n_lats}.nc",
                encoding=encoding,
                compute=True,
                engine="netcdf4",
                format="NETCDF4",
                mode="w",
            )

    elif args.calculate_only == "all":
        if args.pressure_levels is None:
            min_pressure = args.min_pressure
            max_pressure = args.max_pressure
        else:
            min_pressure = args.pressure_levels[0]
            max_pressure = args.pressure_levels[-1]
        ds = load_data(
            file_path=args.file_path,
            variables=args.var,
            sens_model_states=args.sens_model_states,
            only_asc600=args.only_asc600,
            inoutflow_time=args.inoutflow_time,
            min_pressure=min_pressure,
            max_pressure=max_pressure,
            pressure_levels=args.pressure_levels,
            delta_time=args.delta_time,
            n_lons=args.n_lons,
            n_lats=args.n_lats,
            relative_lon_lat=args.relative_lon_lat,
            verbose=args.verbose,
        )
        if args.store_calculated is not None:
            comp = dict(zlib=True, complevel=9)
            encoding = {var: comp for var in ds.data_vars}
            ds.to_netcdf(
                path=f"{args.store_calculated}grid_{args.n_lons}_{args.n_lats}.nc",
                encoding=encoding,
                compute=True,
                engine="netcdf4",
                format="NETCDF4",
                mode="w",
            )

        def handle_plot(ds, var, args):
            file_tmp = args.out_file.split(".")
            if len(file_tmp) > 2:
                filename = ""
                for i in file_tmp[:-1]:
                    filename += i
            else:
                filename = file_tmp[0]
            filetype = file_tmp[-1]
            if "pressure" in ds[var].coords:
                for p in ds["pressure"]:
                    if "Output Parameter" in ds[var].coords:
                        if args.sens_model_states is not None:
                            plot = plot_heatmap(
                                ds=ds.sel(
                                    {
                                        "Output Parameter": args.sens_model_states,
                                        "pressure": p,
                                    }
                                ),
                                col=var,
                                fig_size=args.width / 10,
                                aspect=args.width / args.height,
                            )
                            hv.save(
                                plot,
                                f"{filename}_d{args.sens_model_states}_{var}_{p//100}.{filetype}",
                            )
                        else:
                            for s in ds["Output Parameter"]:
                                plot = plot_heatmap(
                                    ds=ds.sel({"Output Parameter": s, "pressure": p}),
                                    col=var,
                                    fig_size=args.width / 10,
                                    aspect=args.width / args.height,
                                )
                                hv.save(
                                    plot, f"{filename}_d{s}_{var}_{p//100}.{filetype}"
                                )
                    else:
                        plot = plot_heatmap(
                            ds=ds.sel({"pressure": p}),
                            col=var,
                            fig_size=args.width / 10,
                            aspect=args.width / args.height,
                        )
                        hv.save(plot, f"{filename}_{var}_{p//100}.{filetype}")
            else:
                if "Output Parameter" in ds[var].coords:
                    if args.sens_model_states is not None:
                        plot = plot_heatmap(
                            ds=ds.sel({"Output Parameter": args.sens_model_states}),
                            col=var,
                            fig_size=args.width / 10,
                            aspect=args.width / args.height,
                        )
                        hv.save(
                            plot,
                            f"{filename}_d{args.sens_model_states}_{var}.{filetype}",
                        )
                    else:
                        for s in ds["Output Parameter"]:
                            plot = plot_heatmap(
                                ds=ds.sel({"Output Parameter": s}),
                                col=var,
                                fig_size=args.width / 10,
                                aspect=args.width / args.height,
                            )
                            hv.save(plot, f"{filename}_d{s}_{var}.{filetype}")
                else:
                    plot = plot_heatmap(
                        ds=ds,
                        col=var,
                        fig_size=args.width / 10,
                        aspect=args.width / args.height,
                    )
                    hv.save(plot, f"{filename}_{var}.{filetype}")

        if isinstance(args.var, str):
            handle_plot(ds, args.var, args)
        else:
            for var in args.var:
                handle_plot(ds, var, args)
