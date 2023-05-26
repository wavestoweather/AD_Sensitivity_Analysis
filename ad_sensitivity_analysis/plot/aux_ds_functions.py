"""Auxiliary functions to extract information from datasets or data arrays.

"""
import numpy as np


def get_log_threshold(da, lthresh):
    """

    Parameters
    ----------
    da
    lthresh

    Returns
    -------

    """
    if lthresh == 0:
        linthresh = np.nanmin(np.abs(da.where(np.abs(da) > 0)))
    else:
        linthresh = 10**lthresh
    return linthresh


def get_extreme_vals(da):
    """

    Parameters
    ----------
    da
    lthresh

    Returns
    -------

    """
    mini = da.min().values.item()
    maxi = da.max().values.item()
    if mini == 0:
        mini2 = da.where(da > 0).min().values.item() / 10
    else:
        mini2 = mini
    if maxi == 0:
        maxi2 = da.where(da < 0).max().values.item() / 10
    else:
        maxi2 = maxi
    return mini, maxi, mini2, maxi2


def get_regrid_slice(
    ds,
    in_param,
    out_param,
    kind_param,
    time,
    coord,
    pressure,
    fix=False,
    fix_time=False,
):
    """

    Parameters
    ----------
    ds
    in_param
    out_param
    kind_param
    time
    coord
    pressure
    fix
    fix_time

    Returns
    -------

    """
    if in_param == "Top_Parameter":
        return (
            ds.isel({"time": time}).sel({coord: out_param, "pressure": pressure})[
                "Top_Parameter"
            ],
            None,
        )

    if in_param == "counts":
        ds_tmp = ds[in_param]
    elif in_param[0] == "d" and in_param != "deposition":
        ds_tmp = ds.sel({coord: out_param})[f"{kind_param} {in_param}"]
    else:
        ds_tmp = ds[f"{kind_param} {in_param}"]

    if fix_time and not fix:
        ds_fix = ds_tmp.sel({"pressure": pressure})
    elif fix_time and fix:
        ds_fix = ds_tmp
    elif not fix_time and fix:
        ds_fix = ds_tmp.isel({"time": time})
    else:
        ds_fix = ds_tmp
    return ds_tmp.sel({"pressure": pressure}).isel({"time": time}), ds_fix


def set_regrid_missing_kind(ds, kind_param, in_param):
    """

    Parameters
    ----------
    ds
    kind_param
    in_param

    Returns
    -------

    """
    if kind_param == "Std" and f"Std {in_param}" not in ds:
        ds[f"Std {in_param}"] = np.sqrt(ds[f"Var {in_param}"])
    if kind_param == "Var" and f"Var {in_param}" not in ds:
        ds[f"Var {in_param}"] = ds[f"Std {in_param}"] * ds[f"Std {in_param}"]
    return ds
