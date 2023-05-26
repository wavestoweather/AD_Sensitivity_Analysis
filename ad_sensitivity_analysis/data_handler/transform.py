"""Transform trajectory datasets.

"""
import os

import numpy as np
from tqdm.auto import tqdm
import xarray as xr


def fill_time_coords(ds, min_time, max_time, out_coord, traj_enum):
    """
    Fill time coordinate and add NaNs where necessary

    Parameters
    ----------
    ds
    min_time
    max_time
    out_coord
    traj_enum

    Returns
    -------

    """
    tmp_min = ds["time"].min()
    tmp_max = ds["time"].max()
    if tmp_min > min_time:
        time_vals = np.arange(min_time, tmp_min, 30)
        times = int(tmp_min + 20 - min_time) // 30
    elif tmp_max < max_time:
        time_vals = np.arange(tmp_max + 30, max_time + 20, 30)
        times = int(max_time + 20 - tmp_max) // 30
    if tmp_min > min_time or tmp_max < max_time:
        variables = {}
        if out_coord in ds:
            sens_shape = (len(ds[out_coord]), 1, 1, times)
            add_sens = np.reshape(
                np.repeat(np.NaN, times * len(ds[out_coord])),
                sens_shape,
            )
        val_shape = (1, 1, times)
        add_vals = np.reshape(np.repeat(np.NaN, times), val_shape)
        for var in ds:
            if var[0] != "d" or var == "deposition" or var == "dp2h":
                variables[var] = (["ensemble", "trajectory", "time"], add_vals)
            else:
                variables[var] = (
                    [out_coord, "ensemble", "trajectory", "time"],
                    add_sens,
                )
                if var == "asc600":
                    variables[var] = variables[var].where(
                        ~np.isnan(variables[var]), other=0
                    )
        return xr.concat(
            [
                ds,
                xr.Dataset(
                    data_vars=variables,
                    coords={
                        "ensemble": [0],
                        "trajectory": [traj_enum],
                        "time": time_vals,
                    },
                ),
            ],
            join="outer",
            dim="time",
            combine_attrs="override",
        )
    return ds


def get_average(ds):
    """
    Calculate the average over all time steps.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset loaded with 'load_data()' although any dataset with a coordinate 'time' works as well.

    Returns
    -------
    xarray.Dataset where the 'time' coordinate is ditched and averages are calculated.
    """
    return ds.mean(dim="time", skipna=True, keep_attrs=True)


def set_col_types_fp32(ds):
    """

    Parameters
    ----------
    ds

    Returns
    -------

    """
    if ds["time"].dtype == np.float64:
        ds["time"] = ds["time"].astype(np.float32)
    for col in ds:
        if ds[col].dtype == np.float64:
            ds[col] = ds[col].astype(np.float32)
    return ds


def fix_coords(ds):
    """
    Fix coordinates that should not be zero. Replaces those with interpolations between the last and the next
     non-zero coordinate. Only valid for some simulations. Useful when plotting clustered data.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset with "lon", "lat", and "time".

    Returns
    -------
    xarray.Dataset with fixed zero coordinates
    """
    weird_lon = np.count_nonzero(ds["lon"] == 0)
    weird_lat = np.count_nonzero(ds["lat"] == 0)
    # A rare bug can happen where a single lat or lon is set to zero inbetween valid values
    if weird_lon > 0:
        idxs = np.argwhere(ds["lon"].where(ds["lon"] == 0).values == 0)[:, 2]
        # We assume the idxs are consecutive
        n_fills = len(idxs)
        lon_tmp = ds["lon"].isel(
            {
                "time": [idxs[0] - 1, idxs[-1] + 1],
            }
        )
        val_delta = (
            (lon_tmp.isel({"time": 1}) - lon_tmp.isel({"time": 0})) / n_fills
        ).values.item()
        val_start = lon_tmp.isel({"time": 1}).values.item()
        for i_idx, idx in enumerate(idxs):
            ds["lon"][0, 0, idx] = val_start + (i_idx + 1) * val_delta
    if weird_lat > 0:
        idxs = np.argwhere(ds["lat"].where(ds["lat"] == 0).values == 0)[:, 2]
        # We assume the idxs are consecutive
        n_fills = len(idxs)
        lat_tmp = ds["lat"].isel(
            {
                "time": [idxs[0] - 1, idxs[-1] + 1],
            }
        )
        val_delta = (
            (lat_tmp.isel({"time": 1}) - lat_tmp.isel({"time": 0})) / n_fills
        ).values.item()
        val_start = lat_tmp.isel({"time": 1}).values.item()
        for i_idx, idx in enumerate(idxs):
            ds["lat"][0, 0, idx] = val_start + (i_idx + 1) * val_delta
    return ds


def transform_coord(ds, ds_orig=None, zero_condition="ascent_start"):
    """
    Transform latitude and longitude such that
    zero is at the given condition.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset that contains 'lon', 'lat', and 'time_after_ascent'.
    ds_orig : xarray.Dataset
        You may provide the original dataset without a sensitivity analysis to
        get the relative coordinates. This is necessary if the sensitivity analysis
        used a warm-up time such that the start of the ascent is not in the output.
    zero_condition : string
        Possible condition is 'ascent_start' where the coordinates are realative
        to the position to the start of the ascent of the trajectories.

    Returns
    -------
    xarray.Dataset with transformed coordinates
    """
    if zero_condition == "ascent_start" and ds_orig is None:
        lon_start = ds["lon"].where(ds["time_after_ascent"] == 0)
        lat_start = ds["lat"].where(ds["time_after_ascent"] == 0)
    elif zero_condition == "ascent_start":
        lon_start = ds_orig["lon"].where(ds_orig["time_after_asc_start"] == 0)
        lat_start = ds_orig["lat"].where(ds_orig["time_after_asc_start"] == 0)
    shape = np.shape(ds["lon"])[0:-1] + (1,)
    lon_shift = np.zeros(shape)
    lat_shift = np.zeros(shape)
    if "ensemble" in lon_start:
        for ens_idx, ens in enumerate(lon_start["ensemble"]):
            for traj_idx, traj in enumerate(lon_start["trajectory"]):
                lon_shift[ens_idx, traj_idx] = np.nansum(
                    lon_start.sel({"ensemble": ens, "trajectory": traj})
                )
                lat_shift[ens_idx, traj_idx] = np.nansum(
                    lat_start.sel({"ensemble": ens, "trajectory": traj})
                )
    else:
        for traj_idx, traj in enumerate(lon_start["trajectory"]):
            lon_shift[0, traj_idx] = np.nansum(lon_start.sel({"trajectory": traj}))
            lat_shift[0, traj_idx] = np.nansum(lat_start.sel({"trajectory": traj}))
    ds["relative_lon"] = ds["lon"] - lon_shift
    ds["relative_lat"] = ds["lat"] - lat_shift
    ds["relative_lon"].attrs = {
        "long_name": "longitude relative to start of ascent",
        "standard_name": "relative_longitude",
        "units": "degrees",
    }
    ds["relative_lat"].attrs = {
        "long_name": "latitude relative to start of ascent",
        "standard_name": "relative_latitude",
        "units": "degrees",
    }
    return ds


def correct_coordinates(ds):
    """
    An old version of the sensitivity analysis had wrong time after ascent values. This
    function corrects it.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset that contains 'time_after_ascent'.

    Returns
    -------
    xarray.Dataset with corrected coordinates
    """
    ds["time_after_ascent"] = (ds["time_after_ascent"] - 29) * 30
    return ds


def main(parsed_args):
    """

    Parameters
    ----------
    args

    Returns
    -------

    """
    files = [
        f
        for f in os.listdir(parsed_args.file_path)
        if os.path.isfile(parsed_args.file_path + f)
    ]
    comp = {"zlib": True, "complevel": 9}
    for f in tqdm(files) if parsed_args.verbose else files:
        ds = xr.open_dataset(
            parsed_args.file_path + f,
            decode_times=False,
            engine="netcdf4",
        )
        if parsed_args.correct:
            ds = correct_coordinates(ds)
        ds_orig = None
        if parsed_args.file_path_orig is not None:
            ds_orig = xr.open_dataset(
                parsed_args.file_path_orig + f,
                decode_times=False,
                engine="netcdf4",
            )
        ds = transform_coord(ds, ds_orig)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(
            path=parsed_args.store_path + f,
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )


def add_liquid_content(data):
    """

    Parameters
    ----------
    data

    Returns
    -------

    """
    q_total = None
    for col in data.columns:
        if "QR" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
        elif "QC" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
    data["Q_liquid"] = q_total
    return data


# pylint: disable=too-many-branches
def add_cold_content(data):
    """

    Parameters
    ----------
    data

    Returns
    -------

    """
    q_total = None
    for col in data.columns:
        if "QI" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
        elif "QS" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
        elif "QG" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
        elif "QH" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
    data["Q_cold"] = q_total
    return data


def add_sepcific_humidity(data):
    """

    Parameters
    ----------
    data

    Returns
    -------

    """
    data["Specific_Humidity"] = data["QV"] / (data["QV"] + 1)
    return data


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Load data from a sensitivity simulation and transform the coordinate values.
            Currently supports only longitude and latitude transformation.
            Can be used to correct the time after ascent from old sensitivity analysis.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file_path",
        required=True,
        help=textwrap.dedent(
            """\
            Path to a folder with many files from a sensitivity analysis simulation.
            """
        ),
    )
    parser.add_argument(
        "--file_path_orig",
        required=True,
        help=textwrap.dedent(
            """\
            Path to a folder with many files with trajectories from an ICON simulation.
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        default=None,
        help=textwrap.dedent(
            """\
            Path to a folder to store the new files.
            """
        ),
    )
    parser.add_argument(
        "--correct",
        action="store_true",
        help=textwrap.dedent(
            """\
            Correct the time after ascent from old sensitivity analysis.
            """
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=textwrap.dedent(
            """\
            More output, i.e., a progressbar.
            """
        ),
    )
    args = parser.parse_args()
    main(args)
