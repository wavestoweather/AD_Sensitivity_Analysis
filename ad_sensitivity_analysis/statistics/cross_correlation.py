"""Get cross-correlation between two (discrete) time signals

"""
import os

import numpy as np

from tqdm.auto import tqdm
import scipy.signal as scisig
import xarray as xr

from ad_sensitivity_analysis.data_handler.loader import load_dataset_part
from ad_sensitivity_analysis.plot.latexify import param_id_map


def _get_coords(
    ds=None,
    file_path=None,
    phases=False,
    columns=None,
    inoutflow_time=-1,
    phases_arr=np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"]),
):
    """

    Parameters
    ----------
    ds
    file_path
    phases
    columns
    inoutflow_time
    phases_arr

    Returns
    -------

    """
    coords = {
        "Parameter": [],
        "trajectory": 0,
        "ensemble": 0,
    }
    if ds is None:
        files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
        for f in files:
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[
                ["trajectory", "ensemble"]
            ]
            coords["trajectory"] = (
                len(ds["trajectory"])
                if len(ds["trajectory"]) > coords["trajectory"]
                else coords["trajectory"]
            )
            coords["ensemble"] = (
                len(ds["ensemble"])
                if len(ds["ensemble"]) > coords["ensemble"]
                else coords["ensemble"]
            )
        coords["trajectory"] = np.arange(coords["trajectory"])
        coords["ensemble"] = np.arange(coords["ensemble"])
        ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    else:
        files = None
        coords["trajectory"] = np.arange(len(ds["trajectory"]))
        coords["ensemble"] = np.arange(len(ds["ensemble"]))

    if "Output_Parameter_ID" in ds:
        out_param_coord = "Output_Parameter_ID"
        out_params = list(ds[out_param_coord].values)
        coords["Output Parameter"] = [param_id_map[out_p] for out_p in out_params]
    else:
        out_param_coord = "Output Parameter"
        out_params = list(ds[out_param_coord].values)
        coords["Output Parameter"] = out_params
    if columns is None:
        columns = [
            col
            for col in ds
            if col not in ["phase", "step", "asc600", "time_after_ascent"]
        ]
    load_vars = list(columns)
    phases_idx = None
    if phases:
        load_vars.append("phase")
        if ds["phase"].dtype == str:
            phases_idx = phases_arr
        else:
            phases_idx = np.arange(4)

    for p in columns:
        if (p != "deposition") and (p[0] == "d"):
            for out_name in coords["Output Parameter"]:
                coords["Parameter"].append(f"d{out_name}/{p}")
        else:
            coords["Parameter"].append(p)

    if inoutflow_time > 0:
        load_vars.append("asc600")

    return (
        coords,
        out_param_coord,
        phases_idx,
        files,
        load_vars,
        out_params,
    )


def _get_corr(ds_tmp1, ds_tmp2, col1, col2):
    vals1 = ds_tmp1[col1].values
    vals2 = ds_tmp2[col2].values
    # NaNs *should* always be equally distributed
    # in both vals1 and vals2
    mask = ~(np.isnan(vals1) + np.isnan(vals2))
    # We don't care if less than 10 datapoints (=5 minutes) are available
    if np.sum(mask) < 10:
        return np.nan, np.nan

    vals1 = vals1[mask]
    vals2 = vals2[mask]
    vals1 -= np.nanmean(vals1)
    vals2 -= np.nanmean(vals2)

    std1 = np.nanstd(vals1)
    std2 = np.nanstd(vals2)
    if std1 == 0 or std2 == 0:
        # Unfortunately, only zero values might be possible, leading to a variance of zero.
        return np.nan, np.nan
    vals1 /= std1
    vals2 /= std2
    corr = scisig.correlate(
        vals1,
        vals2,
        mode="full",
        method="auto",
    )

    corr_best = np.nanmax(corr)
    offset = np.nanargmax(corr)
    if np.abs(np.nanmax(corr)) > corr_best:
        corr_best = np.nanmin(corr)
        offset = np.nanargmin(corr)
    corr_best /= len(vals1)
    offset = offset - len(vals1) + 1
    return corr_best, offset


# pylint: disable=too-many-locals
def _get_corr_ds_param_col_sens(
    ds, data_inp, out_param_coord, phase_i, f_i, col, param, param_idx, out_params
):
    """

    Parameters
    ----------
    ds
    data_inp
    out_param_coord
    phase_i
    f_i
    col
    param
    param_idx
    out_params

    Returns
    -------

    """
    out_param2 = param.split("/")[0][1::]
    if out_param_coord == "Output_Parameter_ID":
        out_param_ds_idx2 = np.argwhere(np.asarray(param_id_map) == out_param2).item()
    else:
        out_param_ds_idx2 = out_param2
    ds_col = param.split("/")[1]
    for ens_idx, ens in enumerate(ds["ensemble"]):
        for traj_idx, traj in enumerate(ds["trajectory"]):
            for out_param_idx, out_param in enumerate(out_params):
                ds_tmp = ds.sel(
                    {
                        "ensemble": ens,
                        "trajectory": traj,
                        out_param_coord: out_param,
                    }
                )
                ds_tmp2 = ds.sel(
                    {
                        "ensemble": ens,
                        "trajectory": traj,
                        out_param_coord: out_param_ds_idx2,
                    }
                )
                corr_best, offset = _get_corr(ds_tmp, ds_tmp2, col, ds_col)

                if phase_i is None and f_i is None:
                    data_inp[col][1][
                        out_param_idx, ens_idx, traj_idx, param_idx
                    ] = corr_best
                    data_inp[f"Offset {col}"][1][
                        out_param_idx, ens_idx, traj_idx, param_idx
                    ] = offset
                elif phase_i is None:
                    data_inp[col][1][
                        f_i, out_param_idx, ens_idx, traj_idx, param_idx
                    ] = corr_best
                    data_inp[f"Offset {col}"][1][
                        f_i, out_param_idx, ens_idx, traj_idx, param_idx
                    ] = offset
                elif f_i is None:
                    data_inp[col][1][
                        phase_i,
                        out_param_idx,
                        ens_idx,
                        traj_idx,
                        param_idx,
                    ] = corr_best
                    data_inp[f"Offset {col}"][1][
                        phase_i,
                        out_param_idx,
                        ens_idx,
                        traj_idx,
                        param_idx,
                    ] = offset
                else:
                    data_inp[col][1][
                        phase_i,
                        f_i,
                        out_param_idx,
                        ens_idx,
                        traj_idx,
                        param_idx,
                    ] = corr_best
                    data_inp[f"Offset {col}"][1][
                        phase_i,
                        f_i,
                        out_param_idx,
                        ens_idx,
                        traj_idx,
                        param_idx,
                    ] = offset


# pylint: disable=too-many-locals
def _get_corr_ds_param_sens(
    ds, data_outp, out_param_coord, phase_i, f_i, col, param, param_idx
):
    """

    Parameters
    ----------
    ds
    data_outp
    out_param_coord
    phase_i
    f_i
    col
    param
    param_idx

    Returns
    -------

    """
    out_param = param.split("/")[0][1::]
    if out_param_coord == "Output_Parameter_ID":
        out_param_ds_idx = np.argwhere(np.asarray(param_id_map) == out_param).item()
    else:
        out_param_ds_idx = out_param
    ds_col = param.split("/")[1]
    for ens_idx, ens in enumerate(ds["ensemble"]):
        for traj_idx, traj in enumerate(ds["trajectory"]):
            ds_tmp = ds.sel(
                {
                    "ensemble": ens,
                    "trajectory": traj,
                    out_param_coord: out_param_ds_idx,
                }
            )
            corr_best, offset = _get_corr(ds_tmp, ds_tmp, col, ds_col)

            if phase_i is None and f_i is None:
                data_outp[col][1][ens_idx, traj_idx, param_idx] = corr_best
                data_outp[f"Offset {col}"][1][ens_idx, traj_idx, param_idx] = offset
            elif phase_i is None:
                data_outp[col][1][f_i, ens_idx, traj_idx, param_idx] = corr_best
                data_outp[f"Offset {col}"][1][
                    f_i, ens_idx, traj_idx, param_idx
                ] = offset
            elif f_i is None:
                data_outp[col][1][phase_i, ens_idx, traj_idx, param_idx] = corr_best
                data_outp[f"Offset {col}"][1][
                    phase_i, ens_idx, traj_idx, param_idx
                ] = offset
            else:
                data_outp[col][1][
                    phase_i, f_i, ens_idx, traj_idx, param_idx
                ] = corr_best
                data_outp[f"Offset {col}"][1][
                    phase_i, f_i, ens_idx, traj_idx, param_idx
                ] = offset


# pylint: disable=too-many-locals
def _get_corr_ds_col_sens(
    ds, out_param_coord, data_inp, phase_i, f_i, col, param, param_idx, out_params
):
    """

    Parameters
    ----------
    ds
    out_param_coord
    data_inp
    phase_i
    f_i
    col
    param
    param_idx
    out_params

    Returns
    -------

    """
    for ens_idx, ens in enumerate(ds["ensemble"]):
        for traj_idx, traj in enumerate(ds["trajectory"]):
            for out_param_idx, out_param in enumerate(out_params):
                ds_tmp = ds.sel(
                    {
                        "ensemble": ens,
                        "trajectory": traj,
                        out_param_coord: out_param,
                    }
                )
                corr_best, offset = _get_corr(ds_tmp, ds_tmp, col, param)

                if phase_i is None and f_i is None:
                    data_inp[col][1][
                        out_param_idx, ens_idx, traj_idx, param_idx
                    ] = corr_best
                    data_inp[f"Offset {col}"][1][
                        out_param_idx, ens_idx, traj_idx, param_idx
                    ] = offset
                elif phase_i is None:
                    data_inp[col][1][
                        f_i, out_param_idx, ens_idx, traj_idx, param_idx
                    ] = corr_best
                    data_inp[f"Offset {col}"][1][
                        f_i, out_param_idx, ens_idx, traj_idx, param_idx
                    ] = offset
                elif f_i is None:
                    data_inp[col][1][
                        phase_i,
                        out_param_idx,
                        ens_idx,
                        traj_idx,
                        param_idx,
                    ] = corr_best
                    data_inp[f"Offset {col}"][1][
                        phase_i,
                        out_param_idx,
                        ens_idx,
                        traj_idx,
                        param_idx,
                    ] = offset
                else:
                    data_inp[col][1][
                        phase_i,
                        f_i,
                        out_param_idx,
                        ens_idx,
                        traj_idx,
                        param_idx,
                    ] = corr_best
                    data_inp[f"Offset {col}"][1][
                        phase_i,
                        f_i,
                        out_param_idx,
                        ens_idx,
                        traj_idx,
                        param_idx,
                    ] = offset


def _get_corr_ds_no_sens(ds, phase_i, f_i, col, data_outp, param, param_idx):
    """

    Parameters
    ----------
    ds
    phase_i
    f_i
    col
    data_outp
    param
    param_idx

    Returns
    -------

    """
    for ens_idx, ens in enumerate(ds["ensemble"]):
        for traj_idx, traj in enumerate(ds["trajectory"]):
            ds_tmp = ds.sel({"ensemble": ens, "trajectory": traj})
            corr_best, offset = _get_corr(ds_tmp, ds_tmp, col, param)

            if phase_i is None and f_i is None:
                data_outp[col][1][ens_idx, traj_idx, param_idx] = corr_best
                data_outp[f"Offset {col}"][1][ens_idx, traj_idx, param_idx] = offset
            elif phase_i is None:
                data_outp[col][1][f_i, ens_idx, traj_idx, param_idx] = corr_best
                data_outp[f"Offset {col}"][1][
                    f_i, ens_idx, traj_idx, param_idx
                ] = offset
            elif f_i is None:
                data_outp[col][1][phase_i, ens_idx, traj_idx, param_idx] = corr_best
                data_outp[f"Offset {col}"][1][
                    phase_i, ens_idx, traj_idx, param_idx
                ] = offset
            else:
                data_outp[col][1][
                    phase_i, f_i, ens_idx, traj_idx, param_idx
                ] = corr_best
                data_outp[f"Offset {col}"][1][
                    phase_i, f_i, ens_idx, traj_idx, param_idx
                ] = offset


# pylint: disable=too-many-arguments
def _get_corr_ds(
    data_outp,
    data_inp,
    ds,
    params,
    columns,
    out_param_coord,
    out_params,
    phase_i=None,
    f_i=None,
    verbose=False,
    leave=False,
):
    """

    Parameters
    ----------
    data_outp
    data_inp
    ds
    params
    columns
    out_param_coord
    out_params
    phase_i
    f_i
    verbose
    leave

    Returns
    -------

    """
    for col in tqdm(columns, leave=leave) if verbose else columns:
        col_sens = (col[0] == "d") and (col != "deposition")
        for param_idx, param in enumerate(
            tqdm(params, leave=False) if verbose else params
        ):
            param_sens = (param[0] == "d") and (param != "deposition")
            if param_sens and col_sens:
                _get_corr_ds_param_col_sens(
                    ds=ds,
                    out_param_coord=out_param_coord,
                    data_inp=data_inp,
                    phase_i=phase_i,
                    f_i=f_i,
                    col=col,
                    param=param,
                    param_idx=param_idx,
                    out_params=out_params,
                )
            elif param_sens:
                _get_corr_ds_param_sens(
                    ds=ds,
                    data_outp=data_outp,
                    out_param_coord=out_param_coord,
                    phase_i=phase_i,
                    f_i=f_i,
                    col=col,
                    param=param,
                    param_idx=param_idx,
                )
            elif col_sens:
                _get_corr_ds_col_sens(
                    ds=ds,
                    out_param_coord=out_param_coord,
                    data_inp=data_inp,
                    phase_i=phase_i,
                    f_i=f_i,
                    col=col,
                    param=param,
                    param_idx=param_idx,
                    out_params=out_params,
                )
            else:
                _get_corr_ds_no_sens(
                    ds=ds,
                    phase_i=phase_i,
                    f_i=f_i,
                    col=col,
                    data_outp=data_outp,
                    param=param,
                    param_idx=param_idx,
                )


def _get_new_coord(coords, files, columns, phases_arr, phases):
    """

    Parameters
    ----------
    coords
    files
    columns
    phases_arr
    phases

    Returns
    -------

    """
    x_corrs_outp = {}
    x_corrs_inp = {}
    if phases and files is None:
        coords["phase"] = phases_arr
        out_coords = ["phase", "ensemble", "trajectory", "Parameter"]
        in_coords = ["phase", "Output Parameter", "ensemble", "trajectory", "Parameter"]
        out_shape = (
            len(coords["phase"]),
            len(coords["ensemble"]),
            len(coords["trajectory"]),
            len(coords["Parameter"]),
        )
        in_shape = (
            len(coords["phase"]),
            len(coords["Output Parameter"]),
            len(coords["ensemble"]),
            len(coords["trajectory"]),
            len(coords["Parameter"]),
        )
    elif files is None:
        out_coords = ["ensemble", "trajectory", "Parameter"]
        in_coords = ["Output Parameter", "ensemble", "trajectory", "Parameter"]
        out_shape = (
            len(coords["ensemble"]),
            len(coords["trajectory"]),
            len(coords["Parameter"]),
        )
        in_shape = (
            len(coords["Output Parameter"]),
            len(coords["ensemble"]),
            len(coords["trajectory"]),
            len(coords["Parameter"]),
        )
    elif phases:
        coords["file"] = files
        coords["phase"] = phases_arr
        out_coords = ["phase", "file", "ensemble", "trajectory", "Parameter"]
        in_coords = [
            "phase",
            "file",
            "Output Parameter",
            "ensemble",
            "trajectory",
            "Parameter",
        ]
        out_shape = (
            len(coords["phase"]),
            len(coords["file"]),
            len(coords["ensemble"]),
            len(coords["trajectory"]),
            len(coords["Parameter"]),
        )
        in_shape = (
            len(coords["phase"]),
            len(coords["file"]),
            len(coords["Output Parameter"]),
            len(coords["ensemble"]),
            len(coords["trajectory"]),
            len(coords["Parameter"]),
        )
    else:
        coords["file"] = files
        out_coords = ["file", "ensemble", "trajectory", "Parameter"]
        in_coords = ["file", "Output Parameter", "ensemble", "trajectory", "Parameter"]
        out_shape = (
            len(coords["file"]),
            len(coords["ensemble"]),
            len(coords["trajectory"]),
            len(coords["Parameter"]),
        )
        in_shape = (
            len(coords["file"]),
            len(coords["Output Parameter"]),
            len(coords["ensemble"]),
            len(coords["trajectory"]),
            len(coords["Parameter"]),
        )

    for p in columns:
        if (p != "deposition") and (p[0] == "d"):
            x_corrs_inp[p] = (in_coords, np.empty(in_shape))
            x_corrs_inp[p][1][:] = np.nan
            x_corrs_inp[f"Offset {p}"] = (in_coords, np.empty(in_shape))
            x_corrs_inp[f"Offset {p}"][1][:] = np.nan
        else:
            x_corrs_outp[p] = (out_coords, np.empty(out_shape))
            x_corrs_outp[p][1][:] = np.nan
            x_corrs_outp[f"Offset {p}"] = (out_coords, np.empty(out_shape))
            x_corrs_outp[f"Offset {p}"][1][:] = np.nan
    return x_corrs_outp, x_corrs_inp, coords


def cross_correlation(
    ds=None,
    file_path=None,
    phases=False,
    columns=None,
    only_asc600=False,
    inoutflow_time=-1,
    verbose=False,
):
    """
    Calculate the cross-correlation using two (discrete) time signals.
    Non-linear correlations may not be picked up and lead to low cross-correlation values.
    scipy.signal.correlate() uses either fft or direct estimation.
    Calculates the best time for negative or positive correlation and their correlation values.
    The larger absolute value determines if negative or positive correlation is considered.
    Columns are set to zero mean and variance.

    Parameters
    ----------
    ds : xarray.Dataset
        Result of a sensitivity simulation.
    file_path : string
        Path to files with sensitivity simulations.
    phases : bool
        If true, calculate the auto-correlation for each phase separately.
    columns : list of strings
        Columns to load from the dataset. If none is given, all model state variables and model parameters
        will be loaded. Using fewer columns uses significantly less time.
    only_asc600 : bool
        Consider only time steps during the ascend.
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.
    verbose : bool
        Print when a phase is being evaluated and print progressbars.

    Returns
    -------
    xarray.Dataset with coordinates 'Output Parameter' (for sensitivities), 'Parameter' (the first parameter
    for the cross correlation), 'trajectory', and 'ensemble'. If multiple files have been loaded, an additional
    coordinate 'file' with the name of the file is given. If phases is true, then 'phases' is available as
    coordinate.
    Columns are the cross correlation with the respective parameter and 'Offset Parametername' which gives the offset
    best correlation fit.
    """
    phases_arr = np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
    coords, out_param_coord, phases_idx, files, load_vars, out_params = _get_coords(
        ds=ds,
        file_path=file_path,
        phases=phases,
        columns=columns,
        inoutflow_time=inoutflow_time,
        phases_arr=phases_arr,
    )

    x_corrs_outp, x_corrs_inp, coords = _get_new_coord(
        coords, files, columns, phases_arr, phases
    )

    if files is None:
        if phases:
            for phase_i, phase in enumerate(phases_idx):
                if verbose:
                    print(f"Phase: {phase}")
                _get_corr_ds(
                    data_outp=x_corrs_outp,
                    data_inp=x_corrs_inp,
                    ds=ds.where(ds["phase"] == phase),
                    columns=columns,
                    out_param_coord=out_param_coord,
                    params=coords["Parameter"],
                    out_params=out_params,
                    phase_i=phase_i,
                    verbose=verbose,
                    leave=True,
                )
        else:
            _get_corr_ds(
                data_outp=x_corrs_outp,
                data_inp=x_corrs_inp,
                ds=ds,
                columns=columns,
                out_param_coord=out_param_coord,
                params=coords["Parameter"],
                out_params=out_params,
                verbose=verbose,
                leave=True,
            )
    else:
        for f_i, f in enumerate(tqdm(files) if verbose else files):
            ds = load_dataset_part(
                f=file_path + f,
                only_asc600=only_asc600,
                inoutflow_time=inoutflow_time,
                load_params=load_vars,
            )
            if phases:
                for phase_i, phase in enumerate(
                    tqdm(phases_idx, leave=False) if verbose else phases_idx
                ):
                    _get_corr_ds(
                        data_outp=x_corrs_outp,
                        data_inp=x_corrs_inp,
                        ds=ds.where(ds["phase"] == phase),
                        columns=columns,
                        out_param_coord=out_param_coord,
                        params=coords["Parameter"],
                        f_i=f_i,
                        phase_i=phase_i,
                        out_params=out_params,
                        verbose=verbose,
                        leave=False,
                    )
            else:
                _get_corr_ds(
                    data_outp=x_corrs_outp,
                    data_inp=x_corrs_inp,
                    ds=ds,
                    columns=columns,
                    out_param_coord=out_param_coord,
                    params=coords["Parameter"],
                    f_i=f_i,
                    out_params=out_params,
                    verbose=verbose,
                    leave=False,
                )

    data_vars = x_corrs_outp
    for key, value in x_corrs_inp.items():
        data_vars[key] = value
    return xr.Dataset(data_vars=data_vars, coords=coords)
