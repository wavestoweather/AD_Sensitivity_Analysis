import copy
import itertools
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_col
from matplotlib.figure import Figure
import numpy as np
import os
import pandas as pd
import panel as pn
import pickle

try:
    from tqdm.auto import tqdm
except:
    from progressbar import progressbar as tqdm
import scipy.signal as scisig
import seaborn as sns
import sys
import xarray as xr

try:
    from Deriv_dask import Deriv_dask
    from latexify import in_params_dic, physical_params, in_params_notation_mapping
    from plot_mse import load_and_append, reduce_df
    from segment_identifier import d_unnamed
    from create_mse import load_and_append
    import latexify as latexify
except:
    from scripts.Deriv_dask import Deriv_dask
    from scripts.latexify import (
        in_params_dic,
        physical_params,
        in_params_notation_mapping,
    )
    from scripts.plot_mse import load_and_append, reduce_df
    from scripts.segment_identifier import d_unnamed
    from scripts.create_mse import load_and_append
    import scripts.latexify as latexify


def load_ds(f, only_asc600=False, only_phase=None, inoutflow_time=-1, load_params=None):
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

    phases = np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
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

    if inoutflow_time > 0:
        ds_flow = ds.where(ds["asc600"] == 1)["asc600"]
        ds_flow = ds_flow.rolling(
            time=inoutflow_time * 2,  # once for inflow, another for outflow
            min_periods=1,
            center=True,
        ).reduce(np.nanmax)
        ds = ds.where(ds_flow == 1)
    elif only_asc600:
        ds = ds.where(ds["asc600"] == 1)
    if only_phase is not None:
        if ds["phase"].dtype != str and ds["phase"].dtype != np.uint64:
            ds["phase"] = ds["phase"].astype(np.uint64)
            phase_idx = np.argwhere(phases == only_phase)[0].item()
            ds = ds.where(ds["phase"] == phase_idx)
        elif ds["phase"].dtype == str:
            ds = ds.where(ds["phase"] == only_phase)
        else:
            phase_idx = np.argwhere(phases == only_phase)[0].item()
            ds = ds.where(ds["phase"] == phase_idx)
    return ds


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
    Calculate the cross-correlation using two (discrete-)time signals.
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
    if ds is None:
        files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
        n_ensembles = 0
        n_trajectories = 0
        for f in files:
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[
                ["trajectory", "ensemble"]
            ]
            if len(ds["ensemble"]) > n_ensembles:
                n_ensembles = len(ds["ensemble"])
            if len(ds["trajectory"]) > n_trajectories:
                n_trajectories = len(ds["trajectory"])
        ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    else:
        files = None
        n_ensembles = len(ds["ensemble"])
        n_trajectories = len(ds["trajectory"])
    if "Output_Parameter_ID" in ds:
        out_param_coord = "Output_Parameter_ID"
        out_params = list(ds[out_param_coord].values)
        param_names = []
        for out_p in out_params:
            param_names.append(latexify.param_id_map[out_p])
    else:
        out_param_coord = "Output Parameter"
        out_params = list(ds[out_param_coord].values)
        param_names = out_params
    if phases:
        if ds["phase"].dtype == str:
            phases_idx = phases_arr
        else:
            phases_idx = np.arange(4)

    param_coords = []
    ignore_list = ["phase", "step", "asc600", "time_after_ascent"]
    if columns is None:
        columns = []
        for p in list(ds.keys()):
            if p in ignore_list:
                continue
            columns.append(p)
    for p in columns:
        if (p != "deposition") and (p[0] == "d"):
            for out_name in param_names:
                param_coords.append(f"d{out_name}/{p}")
        else:
            param_coords.append(p)

    load_vars = list(columns)
    if phases:
        load_vars.append("phase")
    if inoutflow_time > 0:
        load_vars.append("asc600")

    coords = {
        "Output Parameter": param_names,
        "Parameter": param_coords,
        "trajectory": np.arange(n_trajectories),
        "ensemble": np.arange(n_ensembles),
    }

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

    def get_corr(ds_tmp1, ds_tmp2, col1, col2):
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

    def get_corr_ds(
        data_outp,
        data_inp,
        ds,
        params,
        columns,
        phase_i=None,
        f_i=None,
        verbose=False,
        leave=False,
    ):
        for col_idx, col in enumerate(
            tqdm(columns, leave=leave) if verbose else columns
        ):
            col_sens = (col[0] == "d") and (col != "deposition")
            for param_idx, param in enumerate(
                tqdm(params, leave=False) if verbose else params
            ):
                param_sens = (param[0] == "d") and (param != "deposition")
                if param_sens and col_sens:
                    tmp_param = param.split("/")
                    out_param2 = tmp_param[0][1::]
                    if out_param_coord == "Output_Parameter_ID":
                        out_param_ds_idx2 = np.argwhere(
                            np.asarray(latexify.param_id_map) == out_param2
                        ).item()
                    else:
                        out_param_ds_idx2 = out_param2
                    ds_col = tmp_param[1]
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
                                corr_best, offset = get_corr(
                                    ds_tmp, ds_tmp2, col, ds_col
                                )

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
                elif param_sens:
                    tmp_param = param.split("/")
                    out_param = tmp_param[0][1::]
                    if out_param_coord == "Output_Parameter_ID":
                        out_param_ds_idx = np.argwhere(
                            np.asarray(latexify.param_id_map) == out_param
                        ).item()
                    else:
                        out_param_ds_idx = out_param
                    ds_col = tmp_param[1]
                    for ens_idx, ens in enumerate(ds["ensemble"]):
                        for traj_idx, traj in enumerate(ds["trajectory"]):
                            ds_tmp = ds.sel(
                                {
                                    "ensemble": ens,
                                    "trajectory": traj,
                                    out_param_coord: out_param_ds_idx,
                                }
                            )
                            corr_best, offset = get_corr(ds_tmp, ds_tmp, col, ds_col)

                            if phase_i is None and f_i is None:
                                data_outp[col][1][
                                    ens_idx, traj_idx, param_idx
                                ] = corr_best
                                data_outp[f"Offset {col}"][1][
                                    ens_idx, traj_idx, param_idx
                                ] = offset
                            elif phase_i is None:
                                data_outp[col][1][
                                    f_i, ens_idx, traj_idx, param_idx
                                ] = corr_best
                                data_outp[f"Offset {col}"][1][
                                    f_i, ens_idx, traj_idx, param_idx
                                ] = offset
                            elif f_i is None:
                                data_outp[col][1][
                                    phase_i, ens_idx, traj_idx, param_idx
                                ] = corr_best
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
                elif col_sens:
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
                                corr_best, offset = get_corr(ds_tmp, ds_tmp, col, param)

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
                else:
                    for ens_idx, ens in enumerate(ds["ensemble"]):
                        for traj_idx, traj in enumerate(ds["trajectory"]):
                            ds_tmp = ds.sel({"ensemble": ens, "trajectory": traj})
                            corr_best, offset = get_corr(ds_tmp, ds_tmp, col, param)

                            if phase_i is None and f_i is None:
                                data_outp[col][1][
                                    ens_idx, traj_idx, param_idx
                                ] = corr_best
                                data_outp[f"Offset {col}"][1][
                                    ens_idx, traj_idx, param_idx
                                ] = offset
                            elif phase_i is None:
                                data_outp[col][1][
                                    f_i, ens_idx, traj_idx, param_idx
                                ] = corr_best
                                data_outp[f"Offset {col}"][1][
                                    f_i, ens_idx, traj_idx, param_idx
                                ] = offset
                            elif f_i is None:
                                data_outp[col][1][
                                    phase_i, ens_idx, traj_idx, param_idx
                                ] = corr_best
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

    if files is None:
        if phases:
            for phase_i, phase in enumerate(phases_idx):
                if verbose:
                    print(f"Phase: {phase}")
                get_corr_ds(
                    data_outp=x_corrs_outp,
                    data_inp=x_corrs_inp,
                    ds=ds.where(ds["phase"] == phase),
                    columns=columns,
                    params=param_coords,
                    phase_i=phase_i,
                    verbose=verbose,
                    leave=True,
                )
        else:
            get_corr_ds(
                data_outp=x_corrs_outp,
                data_inp=x_corrs_inp,
                ds=ds,
                columns=columns,
                params=param_coords,
                verbose=verbose,
                leave=True,
            )
    else:
        for f_i, f in enumerate(tqdm(files) if verbose else files):
            ds = load_ds(
                f=file_path + f,
                only_asc600=only_asc600,
                inoutflow_time=inoutflow_time,
                load_params=load_vars,
            )
            if phases:
                for phase_i, phase in enumerate(
                    tqdm(phases_idx, leave=False) if verbose else phases_idx
                ):
                    get_corr_ds(
                        data_outp=x_corrs_outp,
                        data_inp=x_corrs_inp,
                        ds=ds.where(ds["phase"] == phase),
                        columns=columns,
                        params=param_coords,
                        f_i=f_i,
                        phase_i=phase_i,
                        verbose=verbose,
                        leave=False,
                    )
            else:
                get_corr_ds(
                    data_outp=x_corrs_outp,
                    data_inp=x_corrs_inp,
                    ds=ds,
                    columns=columns,
                    params=param_coords,
                    f_i=f_i,
                    verbose=verbose,
                    leave=False,
                )

    data_vars = x_corrs_outp
    for key in x_corrs_inp:
        data_vars[key] = x_corrs_inp[key]
    return xr.Dataset(data_vars=data_vars, coords=coords)


def auto_correlation(
    ds=None,
    file_path=None,
    delay=10,
    phases=False,
    columns=None,
    only_asc600=False,
    inoutflow_time=-1,
    verbose=False,
):
    """
    Estimate auto-correlation via

    {\displaystyle {\hat {R}}(k)={\frac {1}{(n-k)\sigma_{t} \cdot \sigma_{t+k}}}\sum _{t=1}^{n-k}(X_{t}-\mu_{t} )(X_{t+k}-\mu_{t+k} )}

    for each variable and trajectory.
    Mean and standard deviation are estimated for every trajectory and time slice X_{t} and X_{t+k}.
    The mean is set to zero by substracting the mean for every datapoint.
    The variance is set to one by dividing by the variance.

    Parameters
    ds : xarray.Dataset
        Result of a sensitivity simulation.
    file_path : string
        Path to files with sensitivity simulations.
    delay : int or list-like of int
        Define the delay (in time steps) to calculate the auto-correlation or a range of delays.
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
    xarray.Dataset with the model state variables and parameters as column names and their auto-correlation as values.
    The indices are 'phase' (if phases is true), 'delay', 'file' (if file_path is given), 'trajectory', 'ensemble', and
    'Output Parameter' (for sensitivities to model parameters).
    """
    phases_arr = np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
    out_param_coord = ""
    if ds is None:
        files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
        n_ensembles = 0
        n_trajectories = 0
        for f in files:
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[
                ["trajectory", "ensemble"]
            ]
            if len(ds["ensemble"]) > n_ensembles:
                n_ensembles = len(ds["ensemble"])
            if len(ds["trajectory"]) > n_trajectories:
                n_trajectories = len(ds["trajectory"])
        ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    else:
        files = None
        n_ensembles = len(ds["ensemble"])
        n_trajectories = len(ds["trajectory"])
    if "Output_Parameter_ID" in ds:
        out_params = list(ds["Output_Parameter_ID"].values)
        param_names = []
        for out_p in out_params:
            param_names.append(latexify.param_id_map[out_p])
        out_param_coord = "Output_Parameter_ID"
    else:
        out_params = list(ds["Output Parameter"].values)
        param_names = out_params
        out_param_coord = "Output Parameter"

    auto_corr_coords = {
        "Output Parameter": param_names,
        "trajectory": np.arange(n_trajectories),
        "ensemble": np.arange(n_ensembles),
        "delay": [],
    }

    if columns is None:
        columns = list(ds.keys())
        if "phase" in columns:
            columns.remove("phase")
        if "step" in columns:
            columns.remove("step")
        if "asc600" in columns:
            columns.remove("asc600")
        if "time_after_ascent" in columns:
            columns.remove("time_after_ascent")
    load_vars = columns
    if phases:
        auto_corr_coords["phase"] = phases_arr
        load_vars.append("phase")
        if ds["phase"].dtype == str:
            phases_idx = phases_arr
        else:
            phases_idx = np.arange(4)
    if inoutflow_time > 0:
        load_vars.append("asc600")

    if isinstance(delay, int):
        auto_corr_coords["delay"].append(delay)
    else:
        auto_corr_coords["delay"].extend(delay)

    def get_corr(ds_tmp, d, n, data, d_i):
        ds_tmp1 = ds_tmp.isel({"time": np.arange(n - d)})
        ds_tmp2 = ds_tmp.isel({"time": np.arange(d, n)})
        mean1 = ds_tmp1.mean(dim="time")
        mean2 = ds_tmp2.mean(dim="time")
        sigma1 = ds_tmp1.std(dim="time")
        sigma2 = ds_tmp1.std(dim="time")
        no_delay = ds_tmp1 - mean1
        delayed = ds_tmp2 - mean2
        corr = np.asarray(
            (no_delay * delayed).sum(dim="time") / ((n - d) * sigma1 * sigma2)
        )
        corr_shape = np.shape(corr)
        data_shape = np.shape(data[:, :, d_i])
        if corr_shape != data_shape:
            dim_diff = data_shape[1] - corr_shape[1]
            if dim_diff > 0:
                append_arr = np.empty((corr_shape[0], dim_diff))
                append_arr[:] = np.nan
                corr = np.append(corr, append_arr, axis=1)
                corr = np.reshape(corr, (corr_shape[0], data_shape[1]))
            dim_diff = data_shape[0] - corr_shape[0]
            if dim_diff > 0:
                append_arr = np.empty((dim_diff, data_shape[1]))
                append_arr[:] = np.nan
                corr = np.append(corr, append_arr)
                corr = np.reshape(corr, data_shape)
        data[:, :, d_i] = corr

    def get_corr_ds(
        data_outp,
        data_inp,
        ds,
        params,
        phase_i=None,
        f_i=None,
        verbose=False,
        leave=False,
    ):
        n = len(ds["time"])
        for param_i, col in enumerate(tqdm(params, leave=leave) if verbose else params):
            if col[0] != "d" or col == "deposition":
                if isinstance(delay, int):
                    if phase_i is None and f_i is None:
                        get_corr(ds[col], delay, n, data_outp[col][1], 0)
                    elif phase_i is None:
                        get_corr(ds[col], delay, n, data_outp[col][1][f_i], 0)
                    elif f_i is None:
                        get_corr(
                            ds[col],
                            delay,
                            n,
                            data_outp[col][1][phase_i],
                            0,
                        )
                    else:
                        get_corr(
                            ds[col],
                            delay,
                            n,
                            data_outp[col][1][f_i, phase_i],
                            0,
                        )
                else:
                    for d_i, d in enumerate(delay):
                        if phase_i is None and f_i is None:
                            get_corr(ds[col], d, n, data_outp[col][1], d_i)
                        elif phase_i is None:
                            get_corr(ds[col], d, n, data_outp[col][1][f_i], d_i)
                        elif f_i is None:
                            get_corr(
                                ds[col],
                                d,
                                n,
                                data_outp[col][1][phase_i],
                                d_i,
                            )
                        else:
                            get_corr(
                                ds[col],
                                d,
                                n,
                                data_outp[col][1][f_i, phase_i],
                                d_i,
                            )
            else:
                out_p_i = 0
                for out_p, out_name in zip(out_params, param_names):
                    if isinstance(delay, int):
                        if phase_i is None and f_i is None:
                            get_corr(
                                ds[col].sel({out_param_coord: out_p}),
                                delay,
                                n,
                                data_inp[col][1][out_p_i],
                                0,
                            )
                        elif phase_i is None:
                            get_corr(
                                ds[col].sel({out_param_coord: out_p}),
                                delay,
                                n,
                                data_inp[col][1][f_i, out_p_i],
                                0,
                            )
                        elif f_i is None:
                            get_corr(
                                ds[col].sel({out_param_coord: out_p}),
                                delay,
                                n,
                                data_inp[col][1][phase_i, out_p_i],
                                0,
                            )
                        else:
                            get_corr(
                                ds[col].sel({out_param_coord: out_p}),
                                delay,
                                n,
                                data_inp[col][1][f_i, phase_i, out_p_i],
                                0,
                            )
                    else:
                        for d_i, d in enumerate(delay):
                            if phase_i is None and f_i is None:
                                get_corr(
                                    ds[col].sel({out_param_coord: out_p}),
                                    d,
                                    n,
                                    data_inp[col][1][out_p_i],
                                    d_i,
                                )
                            elif phase_i is None:
                                get_corr(
                                    ds[col].sel({out_param_coord: out_p}),
                                    d,
                                    n,
                                    data_inp[col][1][f_i, out_p_i],
                                    d_i,
                                )
                            elif f_i is None:
                                get_corr(
                                    ds[col].sel({out_param_coord: out_p}),
                                    d,
                                    n,
                                    data_inp[col][1][phase_i, out_p_i],
                                    d_i,
                                )
                            else:
                                get_corr(
                                    ds[col].sel({out_param_coord: out_p}),
                                    d,
                                    n,
                                    data_inp[col][1][f_i, phase_i, out_p_i],
                                    d_i,
                                )
                    out_p_i += 1

    if files is None:
        if phases:
            auto_corrs_outp = {}
            auto_corrs_inp = {}
            for param in columns:
                if param[0] != "d" or param == "deposition":
                    auto_corrs_outp[param] = (
                        ["phase", "ensemble", "trajectory", "delay"],
                        np.empty(
                            (
                                len(auto_corr_coords["phase"]),
                                len(ds["ensemble"]),
                                len(ds["trajectory"]),
                                len(auto_corr_coords["delay"]),
                            ),
                        ),
                    )
                    auto_corrs_outp[param][1][:] = np.nan
                else:
                    auto_corrs_inp[param] = (
                        [
                            "phase",
                            "Output Parameter",
                            "ensemble",
                            "trajectory",
                            "delay",
                        ],
                        np.empty(
                            (
                                len(auto_corr_coords["phase"]),
                                len(auto_corr_coords["Output Parameter"]),
                                len(ds["ensemble"]),
                                len(ds["trajectory"]),
                                len(auto_corr_coords["delay"]),
                            ),
                        ),
                    )
                    auto_corrs_inp[param][1][:] = np.nan

            for phase_i, phase in enumerate(auto_corr_coords["phase"]):
                if verbose:
                    print(f"{phase}")
                if ds["phase"].dtype != str:
                    phase_val = phase_i
                else:
                    phase_val = phase
                get_corr_ds(
                    data_outp=auto_corrs_outp,
                    data_inp=auto_corrs_inp,
                    ds=ds,
                    params=columns,
                    phase_i=phase_i,
                    leave=True,
                    verbose=verbose,
                )
        else:
            auto_corrs_outp = {}
            auto_corrs_inp = {}
            for param in columns:
                if param[0] != "d" or param == "deposition":
                    auto_corrs_outp[param] = (
                        ["ensemble", "trajectory", "delay"],
                        np.empty(
                            (
                                len(ds["ensemble"]),
                                len(ds["trajectory"]),
                                len(auto_corr_coords["delay"]),
                            ),
                        ),
                    )
                    auto_corrs_outp[param][1][:] = np.nan
                else:
                    auto_corrs_inp[param] = (
                        ["Output Parameter", "ensemble", "trajectory", "delay"],
                        np.empty(
                            (
                                len(auto_corr_coords["Output Parameter"]),
                                len(ds["ensemble"]),
                                len(ds["trajectory"]),
                                len(auto_corr_coords["delay"]),
                            ),
                        ),
                    )
                    auto_corrs_inp[param][1][:] = np.nan
            get_corr_ds(
                data_outp=auto_corrs_outp,
                data_inp=auto_corrs_inp,
                ds=ds,
                params=columns,
                leave=True,
                verbose=verbose,
            )
    else:
        auto_corr_coords["file"] = files
        if phases:
            auto_corrs_outp = {}
            auto_corrs_inp = {}
            for param in columns:
                if param[0] != "d" or param == "deposition":
                    auto_corrs_outp[param] = (
                        ["file", "phase", "ensemble", "trajectory", "delay"],
                        np.empty(
                            (
                                len(auto_corr_coords["file"]),
                                len(auto_corr_coords["phase"]),
                                len(auto_corr_coords["ensemble"]),
                                len(auto_corr_coords["trajectory"]),
                                len(auto_corr_coords["delay"]),
                            ),
                        ),
                    )
                    auto_corrs_outp[param][1][:] = np.nan
                else:
                    auto_corrs_inp[param] = (
                        [
                            "file",
                            "phase",
                            "Output Parameter",
                            "ensemble",
                            "trajectory",
                            "delay",
                        ],
                        np.empty(
                            (
                                len(auto_corr_coords["file"]),
                                len(auto_corr_coords["phase"]),
                                len(auto_corr_coords["Output Parameter"]),
                                len(auto_corr_coords["ensemble"]),
                                len(auto_corr_coords["trajectory"]),
                                len(auto_corr_coords["delay"]),
                            ),
                        ),
                    )
                    auto_corrs_inp[param][1][:] = np.nan

            for f_i, f in enumerate(tqdm(files) if verbose else files):
                for phase_i, phase in enumerate(
                    tqdm(auto_corr_coords["phase"], leave=False)
                    if verbose
                    else auto_corr_coords["phase"]
                ):
                    ds = load_ds(
                        f=file_path + f,
                        only_asc600=only_asc600,
                        inoutflow_time=inoutflow_time,
                        load_params=load_vars,
                    )
                    if ds["phase"].dtype != str:
                        phase_val = phase_i
                    else:
                        phase_val = phase
                    get_corr_ds(
                        data_outp=auto_corrs_outp,
                        data_inp=auto_corrs_inp,
                        ds=ds.where(ds["phase"] == phase_val),
                        params=columns,
                        phase_i=phase_i,
                        f_i=f_i,
                        verbose=verbose,
                        leave=False,
                    )
        else:
            auto_corrs_outp = {}
            auto_corrs_inp = {}
            for param in columns:
                if param[0] != "d" or param == "deposition":
                    auto_corrs_outp[param] = (
                        ["file", "ensemble", "trajectory", "delay"],
                        np.empty(
                            (
                                len(auto_corr_coords["file"]),
                                len(auto_corr_coords["ensemble"]),
                                len(auto_corr_coords["trajectory"]),
                                len(auto_corr_coords["delay"]),
                            ),
                        ),
                    )
                    auto_corrs_outp[param][1][:] = np.nan
                else:
                    auto_corrs_inp[param] = (
                        ["file", "Output Parameter", "ensemble", "trajectory", "delay"],
                        np.empty(
                            (
                                len(auto_corr_coords["file"]),
                                len(auto_corr_coords["Output Parameter"]),
                                len(auto_corr_coords["ensemble"]),
                                len(auto_corr_coords["trajectory"]),
                                len(auto_corr_coords["delay"]),
                            ),
                        ),
                    )
                    auto_corrs_inp[param][1][:] = np.nan
            for f_i, f in enumerate(tqdm(files) if verbose else files):
                ds = load_ds(
                    f=file_path + f,
                    only_asc600=only_asc600,
                    inoutflow_time=inoutflow_time,
                    load_params=load_vars,
                )
                get_corr_ds(
                    data_outp=auto_corrs_outp,
                    data_inp=auto_corrs_inp,
                    ds=ds,
                    params=columns,
                    phase_i=None,
                    f_i=f_i,
                    verbose=verbose,
                    leave=False,
                )
    data_vars = auto_corrs_outp
    for key in auto_corrs_inp:
        data_vars[key] = auto_corrs_inp[key]
    return xr.Dataset(data_vars=data_vars, coords=auto_corr_coords)


def get_top_list(ds, print_out=True, verbose=True):
    """

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    print_out : bool
        Print the number of parameters and the dataset.
    verbose : bool
        Print the top parameters for each output variable.

    Returns
    -------
    list of output parameters, list of top20 parameters, list of top10 parameters, dict of top20 parameters,
    dict of top10 parameters.
    """
    out_params = [
        "QV",
        "QC",
        "QR",
        "QG",
        "QH",
        "QI",
        "QS",
        "NCCLOUD",
        "NCRAIN",
        "NCGRAUPEL",
        "NCHAIL",
        "NCICE",
        "NCSNOW",
        "QR_OUT",
        "QG_OUT",
        "QH_OUT",
        "QI_OUT",
        "QS_OUT",
        "NR_OUT",
        "NG_OUT",
        "NH_OUT",
        "NI_OUT",
        "NS_OUT",
    ]
    top20_sens_dic = {}
    top10_sens_dic = {}
    tmp20 = []
    tmp10 = []
    if print_out:
        print("\nGet the top parameters for each output variable\n")
    for out_p in out_params:
        if out_p not in ds["Output Parameter"]:
            continue
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        # Skip those that did not appear
        if np.max(df["Predicted Squared Error"]) == 0:
            continue
        top20_sens_dic[out_p] = list(
            np.unique(df.nlargest(20, "Predicted Squared Error")["Input Parameter"])
        )
        top10_sens_dic[out_p] = list(
            np.unique(df.nlargest(10, "Predicted Squared Error")["Input Parameter"])
        )
        tmp20.extend(top20_sens_dic[out_p])
        tmp10.extend(top10_sens_dic[out_p])
        if verbose:
            print(f"###################{out_p}")
            print(f"Top 10: \n{top10_sens_dic[out_p]}")
            print(f"Top 20: \n{top20_sens_dic[out_p]}")
    top20_list = list(set(tmp20))
    top10_list = list(set(tmp10))
    if print_out or verbose:
        print(
            f"Number of distinct parameters by taking the top 20 for everything: {len(top20_list)}"
        )
        print(
            f"Number of distinct parameters by taking the top 10 for everything: {len(top10_list)}"
        )
        print(ds)
    return (
        list(ds["Output Parameter"].values),
        top20_list,
        top10_list,
        top20_sens_dic,
        top10_sens_dic,
    )


def get_magnitude_list(ds, out_params, print_out=True, verbose=True):
    """
    Get the top parameters within one, two, and three orders of magnitude.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    out_params : list-like of strings
        The model state variables for which sensitivities have been calculated for.
    print_out : bool
        Print the top parameters within one order of magnitude and the number of parameters in each magnitude list.
    verbose : bool
        Print the top parameters for each output variable for each magnitude.

    Returns
    -------
    list of parameters within one order of magnitude, list within two orders of magnitudes, list within three orders of magnitudes
    """
    if print_out:
        print(
            "\nInstead of taking the top 10, take the parameters in the same order of magnitude\n"
        )
    top_one_order = []
    top_two_orders = []
    top_three_orders = []

    for out_p in out_params:
        if out_p not in ds["Output Parameter"]:
            continue
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        # Skip those that did not appear
        if np.max(df["Predicted Squared Error"]) == 0:
            continue
        max_order = np.max(df["Predicted Squared Error"])
        top_one_order_tmp = list(
            np.unique(
                df[df["Predicted Squared Error"] >= max_order / 10]["Input Parameter"]
            )
        )
        top_two_orders_tmp = list(
            np.unique(
                df[df["Predicted Squared Error"] >= max_order / 100]["Input Parameter"]
            )
        )
        top_three_orders_tmp = list(
            np.unique(
                df[df["Predicted Squared Error"] >= max_order / 1000]["Input Parameter"]
            )
        )
        top_one_order.extend(top_one_order_tmp)
        top_two_orders.extend(top_two_orders_tmp)
        top_three_orders.extend(top_three_orders_tmp)
        if verbose:
            print(f"###################{out_p}")
            print(f"Top order: \n{top_one_order_tmp}")
            print(f"Top 2 orders: \n{top_two_orders_tmp}")
            print(f"Top 3 orders: \n{top_three_orders_tmp}")
    top_three_orders_list = list(set(top_three_orders))
    top_two_orders_list = list(set(top_two_orders))
    top_one_order_list = list(set(top_one_order))
    if print_out or verbose:
        print(
            f"Number of distinct parameters by taking the top order of magnitude: {len(top_one_order_list)}"
        )
        print(
            f"Number of distinct parameters by taking the top 2 orders of magnitude: {len(top_two_orders_list)}"
        )
        print(
            f"Number of distinct parameters by taking the top 3 orders of magnitude: {len(top_three_orders_list)}"
        )
        print("Parameters within the top order of magnitude:")
        print(top_one_order_list)
    return top_one_order_list, top_two_orders_list, top_three_orders_list


def print_unique_params(top_sens_dic):
    """
    Print the parameters that appear only for a single output variable.

    Parameters
    ----------
    top_sens_dic : dict of list of strings
        The result of get_top_list() or get_magnitude_list()

    Returns
    -------
    String of all printed statements
    """
    text = "\nParameters that appear only for a single output variable\n"
    print(text)
    unique_list = []
    unique_pairing = []
    not_unique_list = []
    for out_p in top_sens_dic:
        tmp = top_sens_dic[out_p]
        for param in tmp:
            if param not in unique_list and param not in not_unique_list:
                unique_list.append(param)
                unique_pairing.append((out_p, param))
            elif param in unique_list:
                # find the index
                idx = np.argwhere(np.asarray(unique_list) == param)
                not_unique_list.append(param)
                del unique_list[idx[0][0]]
                del unique_pairing[idx[0][0]]
    for pair in unique_pairing:
        print(pair)
        text += f"{pair}\n"
    return text


def print_correlation_broad(ds, out_params):
    """
    Print correlation coefficients (Spearman, Pearson, and Kendall) between predicted and actual error using each
    time step individually with all data for each type of state variable (first moment, second moment,
    sedimentation), and for each output variable.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    out_params : list-like of strings
        The model state variables for which sensitivities have been calculated for.

    Returns
    -------
    String of all printed statements
    """
    text = f"\nCorrelation with all data\n"
    print(text)

    def get_corr(df, kind):
        return df[["Predicted Squared Error", "Mean Squared Error"]].corr(kind)[
            "Predicted Squared Error"
        ][1]

    def print_corr(df):
        global text
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        text += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}\n"
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
        df = df.loc[df["Predicted Squared Error"] != 0]
        n = len(np.unique(df["Input Parameter"]))
        text += f"Correlation without zero parameters; total of {n} parameters"
        print(f"Correlation without zero parameters; total of {n} parameters")
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        text += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}"
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")

    df = ds.to_dataframe().reset_index()
    print_corr(df)

    text += "\nFor each output variable type individually\n"
    print("\nFor each output variable type individually\n")

    second_moment = ["QV", "QC", "QR", "QS", "QG", "QH", "QI"]
    first_moment = ["NCCLOUD", "NCRAIN", "NCGRAUPEL", "NCHAIL", "NCICE", "NCSNOW"]
    second_sed = ["QR_OUT", "QG_OUT", "QH_OUT", "QI_OUT", "QS_OUT"]
    first_sed = ["NR_OUT", "NG_OUT", "NH_OUT", "NI_OUT", "NS_OUT"]

    text += "##################First Moment (Number Count)"
    print("##################First Moment (Number Count)")
    df = ds.sel({"Output Parameter": first_moment}).to_dataframe().reset_index()
    print_corr(df)

    text += "##################Second Moment (Mixing Ratio)"
    print("##################Second Moment (Mixing Ratio)")
    df = ds.sel({"Output Parameter": second_moment}).to_dataframe().reset_index()
    print_corr(df)

    text += "##################First Moment Sedimentation (Number Count)"
    print("##################First Moment Sedimentation (Number Count)")
    df = ds.sel({"Output Parameter": first_sed}).to_dataframe().reset_index()
    print_corr(df)

    text += "##################Second Moment Sedimentation (Mixing Ratio)"
    print("##################Second Moment Sedimentation (Mixing Ratio)")
    df = ds.sel({"Output Parameter": second_sed}).to_dataframe().reset_index()
    print_corr(df)

    text += "\nFor each output variable individually\n"
    print("\nFor each output variable individually\n")

    for out_p in out_params:
        out_p = out_p
        text += f"##################{out_p}"
        print(f"##################{out_p}")
        df = ds.sel({"Output Parameter": [out_p]}).to_dataframe().reset_index()
        print_corr(df)
    return text


def print_correlation_mean(ds, out_params):
    """
    Print correlation coefficients (Spearman, Pearson, and Kendall) between predicted and actual error using the mean
    over time after ascent and trajectory for each type of state variable (first moment, second moment,
    sedimentation), and for each output variable.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    out_params : list-like of strings
        The model state variables for which sensitivities have been calculated for.

    Returns
    -------
    String of all printed statements
    """
    text = "\nInstead of looking at correlations within time steps and trajectories individually, "
    text += "take the mean of those and look at the correlation\n"
    text += f"\nCorrelation with all data\n"
    print(text)

    def get_corr(df, kind):
        return df[["Predicted Squared Error", "Mean Squared Error"]].corr(kind)[
            "Predicted Squared Error"
        ][1]

    def print_corr(df):
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        global text
        text += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}"
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
        df = df.loc[df["Predicted Squared Error"] != 0]
        n = len(np.unique(df["Input Parameter"]))
        text += f"Correlation without zero parameters; total of {n} parameters"
        print(f"Correlation without zero parameters; total of {n} parameters")
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        text += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}"
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")

    df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    print_corr(df)

    text += "\nFor each output variable individually\n"
    print("\nFor each output variable individually\n")
    for out_p in out_params:
        out_p = out_p
        text += f"##################{out_p}"
        print(f"##################{out_p}")
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        print_corr(df)

    text += "\nCorrelation taking different number of top parameters\n"
    print("\nCorrelation taking different number of top parameters\n")
    df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    tuples = []
    for n in [10, 20, 50, 100, 150, 200]:
        text += f"##############n={n}#################"
        print(f"##############n={n}#################")
        param_list = []
        for out_p in out_params:
            df_tmp = df.loc[df["Output Parameter"] == out_p]
            if np.max(df_tmp["Predicted Squared Error"]) == 0:
                continue
            param_list.extend(
                df_tmp.nlargest(n, "Predicted Squared Error")["Input Parameter"].values
            )
        param_set = set(param_list)
        df_tmp = df.loc[df["Input Parameter"].isin(param_set)]
        i = len(np.unique(df_tmp["Input Parameter"]))
        text += f"Correlation with zero parameters; total of {i} parameters\n"
        df_tmp2 = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "spearman"
        )["Predicted Squared Error"][1]
        text += f"{df_tmp2}"
        print(f"Correlation with zero parameters; total of {i} parameters")
        print(df_tmp2)
        pearson = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "pearson"
        )["Predicted Squared Error"][1]
        kendall = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "kendall"
        )["Predicted Squared Error"][1]
        text += f"Pearson: {pearson}, Kendall: {kendall}"
        print(f"Pearson: {pearson}, Kendall: {kendall}")
        tuples.append(
            (
                n,
                df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
                    "spearman"
                )["Predicted Squared Error"][1],
            )
        )
    for t in tuples:
        text += f"n={t[0]}, r={t[1]:1.3f}"
        print(f"n={t[0]}, r={t[1]:1.3f}")

    text += "\nCorrelation only for the 75th percentile\n"
    print("\nCorrelation only for the 75th percentile\n")
    params = []
    for out_p in out_params:
        out_p = out_p
        if "NH_OUT" == out_p or "QH_OUT" == out_p:
            continue
        text += f"##################{out_p}\n"
        print(f"##################{out_p}")
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        quan = df["Predicted Squared Error"].quantile(q=0.75)
        df = df.loc[df["Predicted Squared Error"] >= quan]
        n = len(np.unique(df["Input Parameter"]))
        params.extend(np.unique(df["Input Parameter"]))
        text += f"Correlation only with 75th percentile; total of {n} parameters\n"
        print(f"Correlation only with 75th percentile; total of {n} parameters")
        r_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
            "Predicted Squared Error"
        ][1]
        text += f"{r_2}"
        print(r_2)
        pearson = df[["Predicted Squared Error", "Mean Squared Error"]].corr("pearson")[
            "Predicted Squared Error"
        ][1]
        kendall = df[["Predicted Squared Error", "Mean Squared Error"]].corr("kendall")[
            "Predicted Squared Error"
        ][1]
        text += f"Pearson: {pearson}, Kendall: {kendall}"
        print(f"Pearson: {pearson}, Kendall: {kendall}")
    params = list(set(params))
    text += f"Number of parameters: {len(params)}"
    print(f"Number of parameters: {len(params)}")
    df_tmp = (
        ds.sel({"Input Parameter": params})
        .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    r = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
        "Predicted Squared Error"
    ][1]
    text += f"Correlation for all predicted errors with ensemble errors: {r}"
    print(f"Correlation for all predicted errors with ensemble errors: {r}")
    text += f"\nCreate Latex table with Spearman correlations\n"
    print(f"\nCreate Latex table with Spearman correlations\n")
    table_corr = r"""
\begin{table}[hbt]
    \centering
    \begin{tabular}{l|c|c}
        Model State Variable $y_s$                 & $r(y_s)$ without zero sensitivities & $r(y_s)$ \\ \hline
"""
    for out_p in out_params:
        out_p = out_p
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        r = df[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
            "Predicted Squared Error"
        ][1]
        df = df.loc[df["Predicted Squared Error"] != 0]
        r_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
            "Predicted Squared Error"
        ][1]
        table_corr += (
            "\t\t"
            + latexify.parse_word(out_p).title().replace(" Of ", " of ")
            + "\t\t& $ "
            + f"{r_2:1.3f}"
            + " $ & $ "
            + f"{r:1.3f} $ "
            + r"\\"
            + "\n"
        )
    table_corr += r"""    \end{tabular}
    \caption{}
    \label{tab:validate:correlation}
\end{table}"""
    print(table_corr)
    text += table_corr

    text += f"\nCreate Latex table with Pearson correlations\n"
    print(f"\nCreate Latex table with Pearson correlations\n")
    table_corr = r"""
\begin{table}[hbt]
    \centering
    \begin{tabular}{l|c|c}
        Model State Variable $y_s$                 & $r_{\text{pearson}}(y_s)$ without zero sensitivities & $r_{\text{pearson}}(y_s)$ \\ \hline
"""
    for out_p in out_params:
        out_p = out_p
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        r = df[["Predicted Squared Error", "Mean Squared Error"]].corr("pearson")[
            "Predicted Squared Error"
        ][1]
        df = df.loc[df["Predicted Squared Error"] != 0]
        r_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr("pearson")[
            "Predicted Squared Error"
        ][1]
        table_corr += (
            "\t\t"
            + latexify.parse_word(out_p).title().replace(" Of ", " of ")
            + "\t\t& $ "
            + f"{r_2:1.3f}"
            + " $ & $ "
            + f"{r:1.3f} $ "
            + r"\\"
            + "\n"
        )
    table_corr += r"""    \end{tabular}
    \caption{}
    \label{tab:correlation_pearson}
\end{table}"""
    print(table_corr)
    text += table_corr
    return text


def print_latex_tables(ds, top=10, verbose=True):
    """
    Create a string to use in latex for the top n (='top') parameters for each model state variable.
    A model parameter is listed only for the model state variable where the sensitivity is the highest.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    out_params : list-like of strings
        The model state variables for which sensitivities have been calculated for.
    top : int
        The number of top parameters to print a table for.
    verbose : Bool
        If True: print the parameters while building the table.

    Returns
    -------
    sort_key_list: a sorted list of (predicted squared error, model parameter, model state variable,
    string of a row for model parameters in latex) which is sorted by the name of the model state variable.
    table_dic: Dictionary with keys = model parameters where the value is a string of a row of the latex table.
    printed statements as string.
    """
    text = "\nBuild Latex tables\n"
    print(text)
    tmp_df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )

    table_dic = {}
    sort_key_list = []

    latexify_state = {
        "QV": r"\frac{\partial Q_\vapor}{",
        "QC": r"\frac{\partial Q_\cloud}{",
        "QR": r"\frac{\partial Q_\rain}{",
        "QG": r"\frac{\partial Q_\graupel}{",
        "QH": r"\frac{\partial Q_\hail}{",
        "QI": r"\frac{\partial Q_\ice}{",
        "QS": r"\frac{\partial Q_\snow}{",
        "NCCLOUD": r"\frac{\partial N_\cloud}{",
        "NCRAIN": r"\frac{\partial N_\rain}{",
        "NCGRAUPEL": r"\frac{\partial N_\graupel}{",
        "NCHAIL": r"\frac{\partial N_\hail}{",
        "NCICE": r"\frac{\partial N_\ice}{",
        "NCSNOW": r"\frac{\partial N_\snow}{",
        "QR_OUT": r"\frac{\partial Q_{\rain, \text{out}}}{",
        "QG_OUT": r"\frac{\partial Q_{\graupel, \text{out}}}{",
        "QH_OUT": r"\frac{\partial Q_{\hail, \text{out}}}{",
        "QI_OUT": r"\frac{\partial Q_{\ice, \text{out}}}{",
        "QS_OUT": r"\frac{\partial Q_{\snow, \text{out}}}{",
        "NR_OUT": r"\frac{\partial N_{\rain, \text{out}}}{",
        "NG_OUT": r"\frac{\partial N_{\graupel, \text{out}}}{",
        "NH_OUT": r"\frac{\partial N_{\hail, \text{out}}}{",
        "NI_OUT": r"\frac{\partial N_{\ice, \text{out}}}{",
        "NS_OUT": r"\frac{\partial N_{\snow, \text{out}}}{",
    }

    top_10_table = "\\begin{table}[hbt] \n \t\\centering \n \t\\begin{tabular}{ll}"
    top_10_table += "\n \t\t\\textbf{Model State Variable} \t& \\textbf{Top 10 Parameters} \\\\ \\hline \n"
    sedi_latex = ""
    sedi_started = False
    long_table_dic = {}

    for out_p in latexify_state.keys():
        if "OUT" in out_p:
            if sedi_started:
                sedi_latex = sedi_latex[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
            else:
                sedi_latex = "\t\t Sedimentation \t& $ "
                sedi_started = True
        else:
            top_10_table += "\t\t" + latexify.parse_word(out_p) + "\t& $ "
        if verbose:
            print(f"########################### {out_p} ########################")
        df = tmp_df.loc[tmp_df["Output Parameter"] == out_p]
        # Ignore parameters that never appeared in unperturbed versions
        if np.max(df["Predicted Squared Error"]) == 0:
            continue
        if verbose:
            print("sort by sensitivity")
            print(
                df.nlargest(top, "Predicted Squared Error")[
                    ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
                ]
            )
        tmp = df.nlargest(top, "Predicted Squared Error")[
            ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
        ]
        i = 0
        for idx, row in tmp.iterrows():
            if i == 5:
                if "OUT" in out_p:
                    sedi_latex = sedi_latex[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
                else:
                    top_10_table = top_10_table[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
            i += 1
            if "OUT" in out_p:
                sedi_latex += (
                    latexify.parse_word(row["Input Parameter"])
                    .replace("$", "")
                    .replace("\partial", "")
                    + ", "
                )
            else:
                top_10_table += (
                    latexify.parse_word(row["Input Parameter"])
                    .replace("$", "")
                    .replace("\partial", "")
                    + ", "
                )
            found = False
            for val, param, state_var, l_string in sort_key_list:
                if param == row["Input Parameter"] and (
                    val < row["Predicted Squared Error"]
                    or ("N" in state_var and "N" not in out_p)
                ):
                    if "N" not in state_var and "N" in out_p:
                        break

                    found = True
                    if verbose:
                        print(f"Replace ({val}, {param}, {state_var})")
                        print(f"With (")
                        print(row["Predicted Squared Error"], end=", ")
                        print(row["Input Parameter"], end=", ")
                        print(out_p, end=")\n")

                    sort_key_list.remove((val, param, state_var, l_string))
                    break

            if row["Input Parameter"] not in table_dic or found:

                group = None
                for g in latexify.in_params_grouping:
                    if row["Input Parameter"] in latexify.in_params_grouping[g]:
                        group = g

                def latex_my_number(x):
                    if x == 0:
                        return "$ 0.00 $"
                    if x >= 100:
                        exponent = int(np.log10(x))
                        var = x / 10 ** exponent
                        return f"$ {var:2.2f} \\times 10^{ {exponent} } $"
                    elif x < 0.01:
                        exponent = math.floor(np.log10(x))
                        var = x * 10 ** (-exponent)
                        return f"$ {var:2.2f} \\times 10^{ {exponent} } $"
                    else:
                        err = row["Predicted Squared Error"]
                        return f"$ {err:2.2f} $"

                long_string = (
                    latexify.parse_word(row["Input Parameter"]).replace("\partial", "")
                    + " & "
                    + latex_my_number(row["Mean Squared Error"])
                    + " & "
                    + latex_my_number(row["Predicted Squared Error"])
                    + " & "
                    + "\\textbf{"
                    + group.title()
                    + "}: "
                    + latexify.in_params_descr_dic[row["Input Parameter"]]
                    + " \\\\ "
                )
                if out_p not in long_table_dic:
                    long_table_dic[out_p] = [long_string]
                else:
                    long_table_dic[out_p].append(long_string)
                sort_key_list.append(
                    (
                        row["Predicted Squared Error"],
                        row["Input Parameter"],
                        out_p,
                        long_string,
                    )
                )

                table_dic[row["Input Parameter"]] = (
                    "$ \displaystyle "
                    + latexify_state[out_p]
                    + latexify.parse_word(row["Input Parameter"]).replace("$", "")
                    + r"} $ & "
                    + latex_my_number(row["Mean Squared Error"])
                    + " & "
                    + latex_my_number(row["Predicted Squared Error"])
                    + " & "
                    + "\\textbf{"
                    + group.title()
                    + "}: "
                    + latexify.in_params_descr_dic[row["Input Parameter"]]
                    + " \\\\ "
                )
        if "OUT" not in out_p:
            top_10_table = top_10_table[:-2] + " $ \\\\ \n"
        if verbose:
            print(f"sort by Predicted Squared Error")
            print(
                df.nlargest(top, "Predicted Squared Error")[
                    ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
                ]
            )

    top_10_table += sedi_latex[:-2] + " $ \\\\"
    top_10_table += "\n\t\\end{tabular} \n \t\\caption{} \n"
    top_10_table += "\t\\label{tab:} \n \\end{table} \n"
    text += "\nThe table of top 10 parameters for each state variable:\n"
    text += top_10_table
    print("\nThe table of top 10 parameters for each state variable:\n")
    print(top_10_table)

    if verbose:
        print(f"There are {len(table_dic)} different input parameters")

    text += "\nAppendix table of top parameters\n"
    print("\nAppendix table of top parameters\n")
    tmp_sort = sorted(sort_key_list, key=lambda x: (x[2], x[0]), reverse=True)
    sort_dic_long_table = {}
    sort_dic_short_table = {}
    for sens, key, state_variable, l_string in tmp_sort:
        if "_OUT" in state_variable:
            state_variable2 = "Sedimentation"
        else:
            state_variable2 = state_variable
        if state_variable not in sort_dic_long_table:
            sort_dic_long_table[state_variable] = [(sens, key, l_string)]
        else:
            sort_dic_long_table[state_variable].append((sens, key, l_string))

        if state_variable2 not in sort_dic_short_table:
            sort_dic_short_table[state_variable2] = [(sens, key, l_string)]
        else:
            sort_dic_short_table[state_variable2].append((sens, key, l_string))
    table_text = "\\bgroup\n"
    table_text += "\\def\\arraystretch{1.2} %  1 is the default, we want it slightly larger such that exponents are easier to read\n"
    table_text += "\\begin{tabularx}{\\linewidth}{@{}lccX@{}}\n"
    table_text += "\t\\textbf{Model Param.}  & \\textbf{MSD} & \\textbf{Predicted MSD} & \\textbf{Parameter Description}\n"
    table_text += "\t\\endhead\n"
    i = 0
    for state_variable in latexify_state.keys():
        if state_variable not in sort_dic_long_table:
            continue
        table_text += (
            "\t\t\\hline \\multicolumn{4}{c}{"
            + latexify.parse_word(state_variable).title()
            + "}\t \\\\ \\hline\n"
        )
        for _, _, s in sort_dic_long_table[state_variable]:
            table_text += "\t\t"
            table_text += f"{s}"
            i += 1
    table_text += "\t"
    top_str = "ten"
    if top != 10:
        top_str = str(top)
    table_text += (
        r"\caption{This is the set of parameters if we gather the "
        + top_str
        + r" most important ones for each model state variable. The predicted MSD is defined in Equation~\ref{eq:identification:msd_predict}, where we only show the highest predicted MSD among all mass densities unless the parameter did not have an impact on mass densities. In that case, the predicted deviation on number density and precipitation is considered. There are $ "
        + str(i)
        + r" $ different parameters in total.}"
    )
    table_text += "\t\\label{tab:important_params}\n"
    table_text += "\\end{tabularx}\n"
    table_text += "\\egroup\n"
    text += table_text
    print(table_text)
    return sort_key_list, table_dic, text


def print_variable_with_important_params(sort_key_list):
    """
    Print model state variables and their number of model parameters. The number of parameters is determined
    by associating a model parameter x with a model state variable y if the sensitivity dx/dy is highest for y.

    Parameters
    ----------
    sort_key_list
        A sorted list of (predicted squared error, model parameter, model state variable,
        string of a row for model parameters in latex) which is sorted by the name of the model state variable.
        You may generate this using print_late_tables().

    Returns
    -------
    String of all printed statements
    """
    text = (
        "\nWhich model state variable has the most parameters with a high influence?\n"
    )
    print(text)
    state_counts = {}
    for sens, key, state_variable, desc in sort_key_list:
        if state_variable not in state_counts:
            state_counts[state_variable] = 1
        else:
            state_counts[state_variable] += 1
    for state_variable in state_counts:
        print(f"{state_variable}: {state_counts[state_variable]}")
        text += f"{state_variable}: {state_counts[state_variable]}"
    return text


def print_param_types(ds, table_dic):
    """
    Print the type of parameters that are available in the dataset. Types are 'physical',
    'physical (high variability)', 'artificial', and 'artificial (threshold)'. The categories had been determined
    using the feedback of several meteorologists and are merely a guidance, not a definitive categorization.
    In addition, a categorization into geometric, velocity related, exponents, coefficients, and miscellaneous
    parameters is made too.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    table_dic : dict of strings
        Dictionary with keys = model parameters where the value is a string of a row of the latex table.
        This can be generated using print_latex_tables().

    Returns
    -------
    String of all printed statements
    """
    text = "\nHow many type of parameters are there?\n"
    print(text)
    geo_count = 0
    vel_count = 0
    misc_count = 0
    exp_count = 0
    coeff_count = 0
    else_count = 0
    phys_count = 0
    phys_var_count = 0
    art_count = 0
    art_thresh_count = 0
    else_group_count = 0
    for key in table_dic:
        if "geo" in key:
            geo_count += 1
        elif "vel" in key:
            vel_count += 1
        else:
            misc_count += 1

        if "exponent" in table_dic[key].lower():
            exp_count += 1
        elif "coefficient" in table_dic[key].lower():
            coeff_count += 1
        else:
            else_count += 1

        if (
            "physical" in table_dic[key].lower()
            and not "high variability" in table_dic[key].lower()
        ):
            phys_count += 1
        elif (
            "artificial" in table_dic[key].lower()
            and not "threshold" in table_dic[key].lower()
        ):
            art_count += 1
        elif "physical" in table_dic[key].lower():
            phys_var_count += 1
        elif "threshold" in table_dic[key].lower():
            art_thresh_count += 1
        else:
            else_group_count += 1

    total_phys_count = 0
    total_phys_var_count = 0
    total_art_count = 0
    total_art_thresh_count = 0
    total_parameters = len(np.unique(ds["Input Parameter"]))
    for param in np.unique(ds["Input Parameter"]):
        if param in latexify.in_params_grouping["physical"]:
            total_phys_count += 1
        if param in latexify.in_params_grouping["physical (high variability)"]:
            total_phys_var_count += 1
        if param in latexify.in_params_grouping["artificial"]:
            total_art_count += 1
        if param in latexify.in_params_grouping["artificial (threshold)"]:
            total_art_thresh_count += 1

    print(f"There are {total_parameters} many parameters")
    print(f"There are {geo_count} geometric parameters")
    print(f"There are {vel_count} velocity parameters")
    print(f"There are {misc_count} misc. parameters")
    print(f"There are {exp_count} exponents")
    print(f"There are {coeff_count} coefficients")
    print(
        f"There are {else_count} not determined parameters in terms of coefficient or exponent"
    )
    print(f"There are {phys_count} physical parameters (total: {total_phys_count})")
    print(
        f"There are {phys_var_count} physical parameters with a high variability (total: {total_phys_var_count})"
    )
    print(f"There are {art_count} artificial parameters (total: {total_art_count})")
    print(
        f"There are {art_thresh_count} artificial threshold parameters (total: {total_art_thresh_count})"
    )
    print(f"There are {else_group_count} not determined parameters in terms of groups")
    text += f"There are {total_parameters} many parameters\n"
    text += f"There are {geo_count} geometric parameters\n"
    text += f"There are {vel_count} velocity parameters\n"
    text += f"There are {misc_count} misc. parameters\n"
    text += f"There are {exp_count} exponents\n"
    text += f"There are {coeff_count} coefficients\n"
    text += f"There are {else_count} not determined parameters in terms of coefficient or exponent\n"
    text += f"There are {phys_count} physical parameters (total: {total_phys_count})\n"
    text += f"There are {phys_var_count} physical parameters with a high variability (total: {total_phys_var_count})\n"
    text += f"There are {art_count} artificial parameters (total: {total_art_count})\n"
    text += f"There are {art_thresh_count} artificial threshold parameters (total: {total_art_thresh_count})\n"
    text += (
        f"There are {else_group_count} not determined parameters in terms of groups\n"
    )
    return text


def print_large_impact_no_sens(ds, top=50):
    """
    Print every model parameter that shows a large error when perturbed despite a sensitivity value of zero.
    Large is defined as an error within the top n (='top') parameters.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    top : int
        The number of top parameters to consider having a large impact.

    Returns
    -------
    String of all printed statements
    """
    text = f"\nWhich parameters had a large influence (>{top}) despite showing no sensitivity?\n"
    print(text)
    tmp_df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    tmp_df = tmp_df.loc[tmp_df["Predicted Squared Error"] == 0]
    tmp_df2 = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    tmp_df2 = tmp_df2.loc[tmp_df2["Predicted Squared Error"] != 0]
    error_key = "Mean Squared Error"
    out_params = tmp_df["Output Parameter"]

    for out_p in out_params:
        out_p = out_p.item()
        print(f"########################### {out_p} ########################")
        text += f"########################### {out_p} ########################"
        df = tmp_df.loc[tmp_df["Output Parameter"] == out_p]
        df2 = tmp_df2.loc[tmp_df2["Output Parameter"] == out_p]
        nlargest50 = df2.nlargest(top, error_key)[
            ["Input Parameter", error_key, "Predicted Squared Error"]
        ]
        min_sens = nlargest50[error_key].min()
        df_3 = df.loc[df[error_key] >= min_sens]
        text += f"sort by {error_key} with errors > {min_sens}"
        df_tmp_2 = df_3.nlargest(20, error_key)[
            ["Input Parameter", error_key, "Predicted Squared Error"]
        ]
        text += f"{df_tmp_2}"
        print(f"sort by {error_key} with errors > {min_sens}")
        print(df_tmp_2)
    return text


def print_table_top_lists(top_n_lists, top_orders_lists):
    """
    Print a (pandas) table with all parameters for each top variant.

    Parameters
    ----------
    top_n_lists : list of lists of strings
        List of lists of top n parameters generated using get_top_list().
    top_orders_lists : list of lists of strings
        List of lists of parameters within a given magnitude range generated using get_top_list().

    Returns
    -------
    String of all printed statements
    """
    print("\nTable with all parameters for each top variant\n")
    text = "\nTable with all parameters for each top variant\n"
    dict_tops = {}
    max_length = 0
    for i, t in enumerate(top_n_lists):
        t.sort()
        dict_tops[f"top_n_{i}"] = t
        if len(t) > max_length:
            max_length = len(t)
    for i, t in enumerate(top_orders_lists):
        t.sort()
        dict_tops[f"top_order_{i}"] = t
        if len(t) > max_length:
            max_length = len(t)
    for key in dict_tops:
        if len(dict_tops[key]) < max_length:
            dict_tops[key].extend(
                ["-" for _ in range(max_length - len(dict_tops[key]))]
            )
    df = pd.DataFrame(dict_tops)
    print(df)
    text += f"{df}"
    return text


def print_top_parameters(top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic):
    """
    Print the parameters and number of parameters with a sensitivity within a magnitude
    or within the top 10 for each output parameter.

    Parameters
    ----------
    top_magn_set : Set
        Set of parameters within one order of magnitude.
    top10_set : Set
        Set of top 10 parameters for each model state variable.
    top_magn_sens_dic : Dictionary
        Dictionary of model state variables with a list of parameters within one order of magnitude.
    top_sens_dic : Dictionary
        Dictionary of model state variables with a list of top 10 parameters for each.

    Returns
    -------
    String of all printed statements
    """
    text = f"No. of parameters within magnitude of 10**1: {len(top_magn_set)}\n"
    text += f"{top_magn_set}\n"
    text += "The parameters within a magnitude for each output Parameter:\n"
    for out_p in top_magn_sens_dic.keys():
        text += f"~*~*~*~*~*~* {out_p} ~*~*~*~*~*~*\n"
        for param in top_magn_sens_dic[out_p]:
            try:
                text += f"{param}: {in_params_notation_mapping[param][0]}\n"
            except:
                text += f"{param}\n"
        text += "\n"
    text += f"No. of parameters within the top 10: {len(top10_set)}\n"
    text += f"{top10_set}\n"
    text += "The top parameters 10 for each output Parameter:\n"
    for out_p in top_sens_dic.keys():
        text += f"~*~*~*~*~*~* {out_p} ~*~*~*~*~*~*\n"
        text += f"{top_sens_dic[out_p]}\n"
        for param in top_sens_dic[out_p]:
            try:
                text += f"{param}: {in_params_notation_mapping[param][0]}\n"
            except:
                text += f"{param}\n"
        text += "\n"
    print(text)
    return text


def traj_get_sum_derivatives(file_path):
    """
    Calculate the sum of absolute values of sensitivities for each model state variable and model parameter.

    Parameters
    ----------
    file_path : String
        Path to NetCDF-files with sensitivities that used simulation_mode 1 (sensitivity analysis for trajectories)
        in the AD-based C++ simulation.

    Returns
    -------
    Dictionary with tracked output parameters as keys and a pandas.Dataframe with the sum of
    absolute values of the gradients, and a list of strings with the output parameter names.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    out_param_coord = "Output_Parameter_ID"
    if out_param_coord not in ds:
        out_param_coord = "Output Parameter"
        out_params = list(ds[out_param_coord].values)
        param_name = out_params
    else:
        param_name = []
        out_params = list(ds[out_param_coord].values)
        for idx in out_params:
            param_name.append(latexify.param_id_map[idx.values])

    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]

    sums = {}
    for f in tqdm(files):
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
        for out_p, out_name in zip(out_params, param_name):
            ds[in_params] = np.abs(ds[in_params])
            df = (
                ds[in_params]
                .sel({out_param_coord: out_p})
                .sum(dim=["trajectory", "time"], skipna=True)
                .to_dataframe()
                .reset_index()
            )
            df = df[in_params]

            if out_name in sums.keys():
                sums[out_name] += df
            else:
                sums[out_name] = df
    return sums, param_name


def traj_get_top_params(dict_of_df, param_name, n, orders):
    """
    Given the results of get_sums() or get_sums_phase(), get the parameters with the highest sensitivities,
    within a given order of magnitude or within the top n parameters.

    Parameters
    ----------
    dict_of_df : dict of pandas.DataFrame
        Dictionary with (phase +) model state variables as keys and the sums for each gradients as pandas.DataFrame.
    param_name : list-like of strings
        Keys of dict_of_df.
    n : int
        Get the top n parameters for each tracked model state variable
    orders : int or float
        Get the parameters within orders many orders of magnitude for each tracked model state variable

    Returns
    -------
    Set of parameters within the given order of magnitude, set of top n parameters,

    """
    top_sens_dic = {}
    top_magn_sens_dic = {}

    for out_name in param_name:
        tmp_df = dict_of_df[out_name].T
        max_order = np.max(tmp_df[0])
        if max_order == 0:
            continue
        top_sens_dic[out_name] = dict_of_df[out_name].T.nlargest(n, 0).T.columns.values
        top_magn_sens_dic[out_name] = tmp_df[
            tmp_df[0] >= max_order / (10 ** orders)
        ].T.columns.values

    tmp = []
    tmp2 = []
    for out_name in param_name:
        if out_name in top_magn_sens_dic:
            tmp.extend(top_magn_sens_dic[out_name])
            tmp2.extend(top_sens_dic[out_name])
    top_magn_set = set(tmp)
    top10_set = set(tmp2)
    return top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic


def get_edges(
    min_max,
    min_max_in_params,
    in_params,
    param_name,
    additional_params,
    n_bins,
    verbose=False,
):
    """
    Create edges for 1D or 2D histograms.

    Parameters
    ----------
    min_max : Dict of lists with two floats
        Keys are model state variables, values are [minimum, maximum] floats.
    min_max_in_params : Dict of dict of lists with two floats
        First keys are model state variablees, second keys are model parameters. The values are [minimum, maximum]
        floats.
    in_params : List of strings
        Names of model parameters.
    param_name : List of strings
        Name of model state variables for which sensitivities are available.
    additional_params : List of strings
        List of additional model state variables to get edges for.
    n_bins : int
        Number of bins.
    verbose : bool
        More progressbars.

    Returns
    -------
    Dict of lists of floats (edges for model state variables), Dict of dict of lists of floats (edges
    for model parameters for each model state variable defined in param_name).
    """

    edges = {}
    edges_in_params = {}
    for out_p in param_name:
        # -2 because we add one bin each to the left and right to catch any
        # floats just on the edge.
        delta = (min_max[out_p][1] - min_max[out_p][0]) / (n_bins - 2)
        edges[out_p] = np.arange(
            min_max[out_p][0] - delta, min_max[out_p][1] + 0.5 * delta, delta
        )
        edges_in_params[out_p] = {}
        for in_p in in_params:
            if min_max_in_params[out_p][in_p][0] == min_max_in_params[out_p][in_p][
                1
            ] or np.isnan(min_max_in_params[out_p][in_p][0]):
                continue
            delta = (
                min_max_in_params[out_p][in_p][1] - min_max_in_params[out_p][in_p][0]
            ) / (n_bins - 2)
            edges_in_params[out_p][in_p] = np.arange(
                min_max_in_params[out_p][in_p][0] - delta,
                min_max_in_params[out_p][in_p][1] + 0.5 * delta,
                delta,
            )
    if additional_params is not None:
        for out_p in (
            tqdm(additional_params, leave=False) if verbose else additional_params
        ):
            if min_max[out_p][0] == min_max[out_p][1] or np.isnan(min_max[out_p][0]):
                continue
            delta = (min_max[out_p][1] - min_max[out_p][0]) / (n_bins - 2)
            edges[out_p] = np.arange(
                min_max[out_p][0] - delta, min_max[out_p][1] + 0.5 * delta, delta
            )
    return edges, edges_in_params


def get_histogram(
    file_path,
    in_params=None,
    out_params=None,
    n_bins=100,
    additional_params=None,
    only_asc600=False,
    only_phase=None,
    inoutflow_time=-1,
    filter_mag=None,
    means=None,
    verbose=False,
):
    """
    Calculate 1D histograms.

    Parameters
    ----------
    file_path : string
        Path to a folder with many files from a sensitivity analysis simulation.
    in_params : list-like of strings
        List of model parameters.
    out_params : list of strings or ints
        Define for which output parameters the sensitivities should be considered for.
        Can be either ids ("Output_Parameter_ID") or the names of the parameters, e.g., "QV".
    n_bins : int
        Number of bins.
    additional_params : list-like of strings
        Additional parameters that are not model parameters and for which no sensitivity is available but you
        wish to have a histogram for.
    only_asc600 : bool
        Consider only time steps during the ascend.
    only_phase : string
        Consider only time steps with the given phase. Can be combined with only_asc600 or inoutflow_time.
        Possible values are "warm phase", "mixed phase", "ice phase", "neutral phase".
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.
    filter_mag : float
        Filter all values that are more than the given magnitude larger or smaller
        than the mean.
    means : Dictionary of floats
        Mean values for the parameters. Model parameters are dictionaries (model state for which the sensitivitity is for) of dictionaries (the model parameter).
    verbose : bool
        Additional progressbars.

    Returns
    -------
    Dictionary with the following keys and values:
    hist_out_params: Dictionary (keys = model state variable) of arrays with values of the histogram for the given key.
    hist_in_params: Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        values of the histogram for the given key.
    edges_out_params: Dictionary (keys = model state variable) of arrays with the bin edges for the given key.
    edges_in_params: Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        the bin edges for the given keys.
    """
    phases = np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
    if only_phase is not None and only_phase not in phases:
        raise (
            f"You asked for phase {only_phase}, which does not exist."
            f"Possible phases are {phases}"
        )
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = None
    if in_params is None:
        ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
        in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]

    if out_params is None:
        if ds is None:
            ds = xr.open_dataset(
                file_path + files[0], decode_times=False, engine="netcdf4"
            )
        out_params = ds["Output_Parameter_ID"]

        param_name = []
        for idx in out_params:
            param_name.append(latexify.param_id_map[idx.values])
    else:
        param_name = []
        tmp = []
        for idx in out_params:
            if isinstance(idx, int):
                param_name.append(latexify.param_id_map[idx])
            else:
                param_name.append(idx)
                tmp.append(np.argwhere(np.asarray(latexify.param_id_map) == idx).item())
        for i in range(len(tmp)):
            out_params[i] = tmp[i]

    min_max = {}
    min_max_in_params = {}
    for out_p in param_name:
        min_max_in_params[out_p] = {}
    load_params = param_name.copy()
    if inoutflow_time > 0 or only_asc600:
        load_params.append("asc600")
    if only_phase is not None:
        load_params.append("phase")
    if in_params is not None:
        load_params.extend(in_params)
    if additional_params is not None:
        load_params.extend(additional_params)

    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            only_phase=only_phase,
            inoutflow_time=inoutflow_time,
            load_params=load_params,
        )
        for out_p, out_name in (
            tqdm(zip(out_params, param_name), leave=False, total=len(param_name))
            if verbose
            else zip(out_params, param_name)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})

            min_p = np.nanmin(ds_tmp[out_name])
            max_p = np.nanmax(ds_tmp[out_name])
            if out_name in min_max.keys():
                if min_p < min_max[out_name][0]:
                    min_max[out_name][0] = min_p
                if max_p > min_max[out_name][1]:
                    min_max[out_name][1] = max_p
            else:
                min_max[out_name] = [min_p, max_p]
            for in_p in tqdm(in_params, leave=False) if verbose else in_params:
                if filter_mag is not None and means is not None:
                    mean = means[out_name][in_p]
                    filtered_ds = ds_tmp[in_p]
                    filtered_ds = filtered_ds.where(
                        np.abs(np.log10(np.abs(filtered_ds / mean))) <= filter_mag
                    )
                    min_p = np.nanmin(filtered_ds)
                    max_p = np.nanmax(filtered_ds)
                else:
                    min_p = np.nanmin(ds_tmp[in_p])
                    max_p = np.nanmax(ds_tmp[in_p])
                if in_p in min_max_in_params[out_name]:
                    if min_p < min_max_in_params[out_name][in_p][0]:
                        min_max_in_params[out_name][in_p][0] = min_p
                    if max_p > min_max_in_params[out_name][in_p][1]:
                        min_max_in_params[out_name][in_p][1] = max_p
                else:
                    min_max_in_params[out_name][in_p] = [min_p, max_p]
        if additional_params is not None:
            for out_p in (
                tqdm(additional_params, leave=False) if verbose else additional_params
            ):
                min_p = np.nanmin(ds[out_p])
                max_p = np.nanmax(ds[out_p])
                if out_p in min_max.keys():
                    if min_p < min_max[out_p][0]:
                        min_max[out_p][0] = min_p
                    if max_p > min_max[out_p][1]:
                        min_max[out_p][1] = max_p
                else:
                    min_max[out_p] = [min_p, max_p]

    edges, edges_in_params = get_edges(
        min_max=min_max,
        min_max_in_params=min_max_in_params,
        in_params=in_params,
        param_name=param_name,
        additional_params=additional_params,
        n_bins=n_bins,
        verbose=verbose,
    )

    hist = {}
    hist_in_params = {}
    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            only_phase=only_phase,
            inoutflow_time=inoutflow_time,
            load_params=load_params,
        )
        for out_p, out_name in (
            tqdm(zip(out_params, param_name), leave=False, total=len(param_name))
            if verbose
            else zip(out_params, param_name)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            hist_tmp, _ = np.histogram(ds[out_name], edges[out_name])
            if out_name in hist:
                hist[out_name] += hist_tmp
            else:
                hist[out_name] = hist_tmp
                hist_in_params[out_name] = {}
            for in_p in tqdm(in_params, leave=False) if verbose else in_params:
                if in_p not in edges_in_params[out_name]:
                    continue
                hist_tmp, _ = np.histogram(
                    ds_tmp[in_p], edges_in_params[out_name][in_p]
                )
                if in_p in hist_in_params[out_name]:
                    hist_in_params[out_name][in_p] += hist_tmp
                else:
                    hist_in_params[out_name][in_p] = hist_tmp
        if additional_params is not None:
            for out_p in (
                tqdm(additional_params, leave=False) if verbose else additional_params
            ):
                if out_p not in edges:
                    # In case no values for the additional parameters are available, i.e.,
                    # because only_phase is warm phase but QI is given in additonal_params.
                    continue
                hist_tmp, _ = np.histogram(ds[out_p], edges[out_p])
                if out_p in hist:
                    hist[out_p] += hist_tmp
                else:
                    hist[out_p] = hist_tmp

    return {
        "hist_out_params": hist,
        "hist_in_params": hist_in_params,
        "edges_out_params": edges,
        "edges_in_params": edges_in_params,
    }


def get_histogram_cond(
    file_path,
    in_params=None,
    out_params=None,
    conditional_hist=[],
    n_bins=100,
    additional_params=[],
    only_asc600=False,
    only_phase=None,
    inoutflow_time=-1,
    filter_mag=None,
    means=None,
    verbose=False,
):
    """
    Get a 2D histogram where 'cond' is the parameter for which the edges are calculated. The final 2D histogram
    is the count for values that are in the bin defined by 'cond' and in the bin calculated
    from in_params or out_params.

    Parameters
    ----------
    file_path : string
        Path to a folder with many files from a sensitivity analysis simulation.
    in_params : list-like of strings
        List of model parameters.
    out_params : list of strings or ints
        Define for which output parameters the sensitivities should be considered for.
        Can be either ids ("Output_Parameter_ID") or the names of the parameters, e.g., "QV".
    conditional_hist : string
        The model state variable (aka output parameter) for which additional edges shall be calculated.
    n_bins : int
        Number of bins.
    additional_params : list-like of strings
        Additional parameters that are not model parameters and for which no sensitivity is available but you
        wish to have a histogram for.
    only_asc600 : bool
        Consider only time steps during the ascend.
    only_phase : string
        Consider only time steps with the given phase. Can be combined with only_asc600 or inoutflow_time.
        Possible values are "warm phase", "mixed phase", "ice phase", "neutral phase".
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.
    filter_mag : float
        Filter all values that are more than the given magnitude larger or smaller
        than the mean.
    means : Dictionary of floats
        Mean values for the parameters. Model parameters are dictionaries (model state for which the sensitivitity is for) of dictionaries (the model parameter).
    verbose : bool
        Additional progressbars.

    Returns
    -------
    Dictionary with the following keys:
    'edges_out_params': Dictionary where the keys are model state variables and the values are arrays of bin edges.
    'edges_in_params': Dictionary where the keys are model state variables for which sensitivities are available
        and the values are dictionaries of model parameters with arrays of bin edges.
    model state variables: Each model state variable has a dictionary for 'hist_out_params' and 'hist_in_params'.
    'hist_out_params' is a dictionary of model state variables with arrays of bin counts.
    'hist_in_params' is a dictionary of model state variables for which sensitivities are available
        and the values are dictionaries of model parameters with arrays of bin counts.
    """
    phases = np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
    if only_phase is not None and only_phase not in phases:
        raise (
            f"You asked for phase {only_phase}, which does not exist."
            f"Possible phases are {phases}"
        )
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = None
    if in_params is None:
        ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
        in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]

    if out_params is None:
        if ds is None:
            ds = xr.open_dataset(
                file_path + files[0], decode_times=False, engine="netcdf4"
            )
        out_params = list(ds["Output_Parameter_ID"].values)
        param_name = []
        for idx in out_params:
            param_name.append(latexify.param_id_map[idx])
    else:
        param_name = []
        tmp = []
        for idx in out_params:
            if isinstance(idx, int):
                param_name.append(latexify.param_id_map[idx])
            else:
                param_name.append(idx)
                tmp.append(np.argwhere(np.asarray(latexify.param_id_map) == idx).item())
        for i in range(len(tmp)):
            out_params[i] = tmp[i]

    for cond in conditional_hist:
        if cond not in additional_params:
            additional_params.append(cond)

    min_max = {}
    min_max_in_params = {}
    for out_p in param_name:
        min_max_in_params[out_p] = {}

    load_params = param_name.copy()
    if inoutflow_time > 0 or only_asc600:
        load_params.append("asc600")
    if only_phase is not None:
        load_params.append("phase")
    if in_params is not None:
        load_params.extend(in_params)
    if additional_params is not None:
        load_params.extend(additional_params)

    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            only_phase=only_phase,
            inoutflow_time=inoutflow_time,
            load_params=load_params,
        )
        for out_p, out_name in (
            tqdm(zip(out_params, param_name), leave=False, total=len(param_name))
            if verbose
            else zip(out_params, param_name)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            min_p = np.nanmin(ds_tmp[out_name])
            max_p = np.nanmax(ds_tmp[out_name])
            if out_name in min_max.keys():
                if min_p < min_max[out_name][0]:
                    min_max[out_name][0] = min_p
                if max_p > min_max[out_name][1]:
                    min_max[out_name][1] = max_p
            else:
                min_max[out_name] = [min_p, max_p]
            for in_p in tqdm(in_params, leave=False) if verbose else in_params:
                if filter_mag is not None and means is not None:
                    mean = means[out_name][in_p]
                    filtered_ds = ds_tmp[in_p]
                    filtered_ds = filtered_ds.where(
                        np.abs(np.log10(np.abs(filtered_ds / mean))) <= filter_mag
                    )
                    min_p = np.nanmin(filtered_ds)
                    max_p = np.nanmax(filtered_ds)
                else:
                    min_p = np.nanmin(ds_tmp[in_p])
                    max_p = np.nanmax(ds_tmp[in_p])
                if in_p in min_max_in_params[out_name]:
                    if min_p < min_max_in_params[out_name][in_p][0]:
                        min_max_in_params[out_name][in_p][0] = min_p
                    if max_p > min_max_in_params[out_name][in_p][1]:
                        min_max_in_params[out_name][in_p][1] = max_p
                else:
                    min_max_in_params[out_name][in_p] = [min_p, max_p]
        for out_p in (
            tqdm(additional_params, leave=False) if verbose else additional_params
        ):
            min_p = np.nanmin(ds[out_p])
            max_p = np.nanmax(ds[out_p])
            if out_p in min_max.keys():
                if min_p < min_max[out_p][0]:
                    min_max[out_p][0] = min_p
                if max_p > min_max[out_p][1]:
                    min_max[out_p][1] = max_p
            else:
                min_max[out_p] = [min_p, max_p]

    edges, edges_in_params = get_edges(
        min_max=min_max,
        min_max_in_params=min_max_in_params,
        in_params=in_params,
        param_name=param_name,
        additional_params=additional_params,
        n_bins=n_bins,
        verbose=verbose,
    )

    hist_conditional = {
        "edges_out_params": edges,
        "edges_in_params": edges_in_params,
    }

    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            only_phase=only_phase,
            inoutflow_time=inoutflow_time,
            load_params=load_params,
        )
        for cond in (
            tqdm(conditional_hist, leave=False) if verbose else conditional_hist
        ):
            if cond not in edges:
                continue
            if cond not in hist_conditional:
                hist_conditional[cond] = {
                    "hist_out_params": {},
                    "hist_in_params": {},
                }
            for out_p, out_name in (
                tqdm(zip(out_params, param_name), leave=False, total=len(param_name))
                if verbose
                else zip(out_params, param_name)
            ):
                if out_name not in edges:
                    continue
                ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
                hist_tmp, _, _ = np.histogram2d(
                    ds[cond].values.flatten(),
                    ds[out_name].values.flatten(),
                    [edges[cond], edges[out_name]],
                )
                if out_name in hist_conditional[cond]["hist_out_params"]:
                    hist_conditional[cond]["hist_out_params"][out_name] += hist_tmp
                else:
                    hist_conditional[cond]["hist_out_params"][out_name] = hist_tmp
                    hist_conditional[cond]["hist_in_params"][out_name] = {}
                for in_p in tqdm(in_params, leave=False) if verbose else in_params:
                    if in_p not in edges_in_params[out_name]:
                        continue
                    hist_tmp, _, _ = np.histogram2d(
                        ds_tmp[cond].values.flatten(),
                        ds_tmp[in_p].values.flatten(),
                        [edges[cond], edges_in_params[out_name][in_p]],
                    )
                    if in_p in hist_conditional[cond]["hist_in_params"][out_name]:
                        hist_conditional[cond]["hist_in_params"][out_name][
                            in_p
                        ] += hist_tmp
                    else:
                        hist_conditional[cond]["hist_in_params"][out_name][
                            in_p
                        ] = hist_tmp
            for out_p in (
                tqdm(additional_params, leave=False) if verbose else additional_params
            ):
                if out_p == cond:
                    continue
                if out_p not in edges:
                    # In case no values for the additional parameters are available, i.e.,
                    # because only_phase is warm phase but QI is given in additonal_params.
                    continue
                hist_tmp, _, _ = np.histogram2d(
                    ds[cond].values.flatten(),
                    ds[out_p].values.flatten(),
                    [edges[cond], edges[out_p]],
                )
                if out_p in hist_conditional[cond]["hist_out_params"]:
                    hist_conditional[cond]["hist_out_params"][out_p] += hist_tmp
                else:
                    hist_conditional[cond]["hist_out_params"][out_p] = hist_tmp

    return hist_conditional


def traj_plot_histogram_out(
    out_params,
    filename,
    edges,
    hist,
    log=False,
    width=24,
    height=12,
    font_scale=None,
    title=None,
    save=True,
    interactive=False,
    latex=False,
    ticks_offset=None,
    verbose=False,
):
    """
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot the histogram of an output parameter,
    i.e., QV, latent_heat, latent_cool, etc.

    Parameters
    ----------
    out_params : string or list-like of strings
        Output parameter name or multiple output parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges : Dictionary of list-like of float
        Edges for the histogram. Keys must be in out_params.
    hist : Dictionary of list-like of int
        Number of entries for each bin. Keys must be in out_params.
    log : bool
        Plot the y-axis using log-scale.
    width : float
        Width in inches
    height : float
        Height in inches
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    interactive : bool
        Create a figure for interactive plotting.
    title : string
        Title of the histogram. If none is given, a title will be generated.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.
    ticks_offset : int
        If none, the number of ticks is calculated based on the width of the plot. Otherwise,
        every other "ticks_offset" is used.
    verbose : bool
        Print some additional information.

    Returns
    -------
    If interactive == False, returns None. If interactive == True, returns the matplotlib.figure.Figure with the plot drawn onto it.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    if interactive:
        sns.set()
        fig = Figure()
        ax = fig.subplots()
    else:
        fig = None
        ax = None

    def plot_hist(out_p, title=None, ax=None):
        if title is None:
            title = f"Histogram for {latexify.parse_word(out_p)}"
        if verbose:
            print(f"Plotting histogram for {out_p}")

        if ticks_offset is None:
            if 24 > width > 5:
                local_ticks_offset = 24 // (width - 5)
            elif width <= 5:
                local_ticks_offset = 6
        else:
            local_ticks_offset = ticks_offset
        if ticks_offset >= len(edges[out_p]) - 1:
            local_ticks_offset = len(edges[out_p]) - 2
        if ax is None:
            ax = sns.barplot(
                x=edges[out_p][:-1], y=hist[out_p], log=log, color="seagreen"
            )
        else:
            sns.barplot(
                x=edges[out_p][:-1], y=hist[out_p], log=log, color="seagreen", ax=ax
            )
        x_labels = [f"{tick:1.1e}" for tick in edges[out_p][:-1:local_ticks_offset]]
        x_ticks = np.arange(0, len(edges[out_p]) - 1, local_ticks_offset)
        ax.set(xticks=x_ticks)
        _ = ax.set_xticklabels(
            x_labels, rotation=45, ha="right", rotation_mode="anchor"
        )
        if font_scale is None:
            _ = ax.set_title(title)
            ax.set_ylabel(out_p)
        else:
            ax.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
            _ = ax.set_title(title, fontsize=int(12 * font_scale))
            ax.set_ylabel(f"# Entries", fontsize=int(11 * font_scale))
            xlabel = f"Bins of {latexify.parse_word(out_p)}"
            ax.set_xlabel(xlabel, fontsize=int(11 * font_scale))

        if filename is not None and save:
            plt.tight_layout()
            try:
                i = 0
                store_type = filename.split(".")[-1]
                store_path = filename[0 : -len(store_type) - 1]
                save_name = store_path + "_" + out_p + "_{:03d}.".format(i) + store_type
                while os.path.isfile(save_name):
                    i = i + 1
                    save_name = (
                        store_path + "_" + out_p + "_{:03d}.".format(i) + store_type
                    )
                if not interactive:
                    fig = ax.get_figure()
                fig.savefig(save_name, bbox_inches="tight", dpi=300)
            except:
                print(f"Storing to {save_name} failed.", file=sys.stderr)
            if not interactive:
                plt.clf()

    if isinstance(out_params, list):
        for out_p in tqdm(out_params):
            plot_hist(out_p, title, None)
    else:
        plot_hist(out_params, title, ax)
    if verbose:
        print("All plots finished!")
    if interactive:
        # Just a failsafe to avoid repeating error messages with erroneous filepaths.
        save = False
        filename = None
    return fig


def traj_plot_histogram_inp(
    in_params,
    filename,
    edges_in_params,
    hist_in_params,
    log=False,
    font_scale=None,
    width=24,
    height=12,
    title=None,
    save=True,
    interactive=False,
    latex=False,
    ticks_offset=None,
    verbose=False,
):
    """
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot three histograms per image with
    \partial output / \partial model parameter
    where output is QV, latent_heat and latent_cool. Plot one image per model_parameter.

    Parameters
    ----------
    in_params : string or list-like of strings
        Model parameter name or multiple input parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges_in_params : Dictionary of dictionary of list-like of float
        Edges for the histogram. First keys are output parameters, second keys must be in in_params.
    hist_in_params : Dictionary of dictionary of list-like of int
        Number of entries for each bin. First keys are output parameters, second keys must be in in_params.
    log : bool
        Plot the y-axis using log-scale.
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.
    ticks_offset : int
        If none, the number of ticks is calculated based on the width of the plot. Otherwise,
        every other "ticks_offset" is used.
    interactive : bool
        Create a figure for interactive plotting.

    Returns
    -------
    If interactive == False, returns None. If interactive == True, returns the matplotlib.figure.Figure with the plot drawn onto it.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    out_params = list(edges_in_params.keys())
    if interactive:
        sns.set()
        fig = Figure()
        axs = fig.subplots(
            nrows=len(out_params),
            ncols=1,
        )
    else:
        fig = None
        axs = None

    def plot_hist(out_params, in_p, title=None, axs=None):

        if len(out_params) != 3:
            print("The number of output params should be three.")
            print("Future versions will support varying numbers.")
        if axs is None:
            ax1 = plt.subplot(311)
            ax2 = plt.subplot(312)
            ax3 = plt.subplot(313)

        else:
            ax1 = axs[0]
            ax2 = axs[1]
            ax3 = axs[2]
        if verbose:
            print(f"Plotting histogram w.r.t. {in_p}")

        def create_fig(ax, out_p, in_p, title=None):
            if in_p not in edges_in_params[out_p]:
                return
            if title is None:
                in_p_latex = latexify.parse_word(in_p).replace(r"\partial", "")
                title = f"{latexify.parse_word(out_p)} w.r.t. {in_p_latex}"
            ax_t = sns.barplot(
                x=edges_in_params[out_p][in_p][:-1],
                y=hist_in_params[out_p][in_p],
                color="seagreen",
                ax=ax,
            )
            if log:
                ax_t.set_yscale("log")
            if ticks_offset is None:
                if 24 > width > 5:
                    local_ticks_offset = 24 // (width - 5)
                elif width <= 5:
                    local_ticks_offset = 6
            else:
                local_ticks_offset = ticks_offset
            if local_ticks_offset >= len(edges_in_params[out_p][in_p][:-1]) - 1:
                local_ticks_offset = len(edges_in_params[out_p][in_p][:-1]) - 2
            x_labels = [
                f"{tick:1.1e}"
                for tick in edges_in_params[out_p][in_p][:-1:local_ticks_offset]
            ]
            x_ticks = np.arange(
                0, len(edges_in_params[out_p][in_p][:-1]) - 1, local_ticks_offset
            )
            ax_t.set(xticks=x_ticks)
            _ = ax_t.set_xticklabels(
                x_labels, rotation=45, ha="right", rotation_mode="anchor"
            )
            if font_scale is None:
                _ = ax_t.set_title(title)
            else:
                ax_t.tick_params(
                    axis="both", which="major", labelsize=int(10 * font_scale)
                )
                _ = ax_t.set_title(title, fontsize=int(12 * font_scale))

        create_fig(ax1, out_params[0], in_p, title)
        create_fig(ax2, out_params[1], in_p, title)
        create_fig(ax3, out_params[2], in_p, title)

        if interactive:
            # The view in jupyter notebooks may be different from the stored one without
            # the given padding.
            fig.tight_layout(h_pad=1)

        if filename is not None and save:
            try:
                i = 0
                store_type = filename.split(".")[-1]
                store_path = filename[0 : -len(store_type) - 1]
                save_name = store_path + "_" + in_p + "_{:03d}.".format(i) + store_type
                while os.path.isfile(save_name):
                    i = i + 1
                    save_name = (
                        store_path + "_" + in_p + "_{:03d}.".format(i) + store_type
                    )
                if interactive:
                    fig.savefig(save_name, bbox_inches="tight", dpi=300)
                else:
                    plt.savefig(save_name, bbox_inches="tight", dpi=300)
            except:
                print(f"Storing to {save_name} failed.", file=sys.stderr)
            if not interactive:
                plt.clf()

    if isinstance(in_params, list):
        for in_p in tqdm(in_params):
            plot_hist(out_params, in_p, title, axs)
    else:
        plot_hist(out_params, in_params, title, axs)

    if verbose:
        print("All plots finished!")
    if interactive:
        # Just a failsafe to avoid repeating error messages with erroneous filepaths.
        save = False
        filename = None
    return fig


def plot_heatmap_traj(
    in_params,
    filename,
    edges_in_params,
    hist_in_params,
    width=24,
    height=24,
    title=None,
    verbose=False,
):
    """
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot a heatmap "model parameters" over
    "bin number" where each row is another histogram for a certain model parameter.

    Parameters
    ----------
    in_params : string or list-like of strings
        Model parameter name or multiple input parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges_in_params : Dictionary of dictionary list-like of float
        Edges for the histogram. First keys are output parameters, second keys must be in in_params.
    hist_in_params : Dictionary of dictionary list-like of int
        Number of entries for each bin. First keys are output parameters, second keys must be in in_params.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    verbose : bool
        Additional output.
    """
    sns.set(rc={"figure.figsize": (width, height)})
    # sort the histogram by a simple similarity metric. It is not perfect but better than random.
    out_params = list(edges_in_params.keys())
    for out_p in tqdm(out_params):
        in_params_tmp = []
        for in_p in in_params:
            if in_p in hist_in_params[out_p]:
                in_params_tmp.append(in_p)

        corr_matrix = np.zeros((len(in_params_tmp), len(in_params)))
        if verbose:
            print(f"Create similarity matrix for {out_p}")
        for i in range(len(in_params_tmp)):
            for j in range(len(in_params_tmp)):
                if i == j:
                    corr_matrix[i, j] = 1
                if i >= j:
                    continue
                corr_matrix[i, j] = np.corrcoef(
                    hist_in_params[out_p][in_params_tmp[i]],
                    hist_in_params[out_p][in_params_tmp[j]],
                )[0][1]
                corr_matrix[j, i] = corr_matrix[i, j]
        if verbose:
            print("Sort parameters according to correlation matrix")
        best_row = 0
        best_value = -1
        for i in range(len(in_params_tmp)):
            for j in range(len(in_params_tmp)):
                if i >= j:
                    continue
                if best_value == -1:
                    best_value = corr_matrix[i, j]
                    best_row = i
                elif best_value < corr_matrix[i, j]:
                    best_value = corr_matrix[i, j]
                    best_row = i
        in_params_sorted = [in_params_tmp[best_row]]
        previous_rows = [best_row]
        if verbose:
            print("Generate matrix for plotting")
        for i in range(len(in_params_tmp) - 1):
            next_row = 0
            next_best_value = -1
            for j in range(len(in_params_tmp)):
                if j in previous_rows:
                    continue
                if next_best_value == -1:
                    next_best_value = corr_matrix[i, j]
                    next_row = j
                elif next_best_value < corr_matrix[i, j]:
                    next_best_value = corr_matrix[i, j]
                    next_row = j
            in_params_sorted.append(in_params_tmp[next_row])
            previous_rows.append(next_row)
        # old version with similarity calculated as difference in each time step.
        # similarity_matrix = np.zeros((len(in_params_tmp), len(in_params)))
        # if verbose:
        #     print(f"Create similarity matrix for {out_p}")
        # for i in range(len(in_params_tmp)):
        #     for j in range(len(in_params_tmp)):
        #         if i >= j:
        #             continue
        #         similarity_matrix[i, j] = np.sum(
        #             np.abs( hist_in_params[out_p][in_params_tmp[i]] - hist_in_params[out_p][in_params_tmp[j]])
        #         )
        #         similarity_matrix[j, i] = similarity_matrix[i, j]
        # # get the row with the best/smallest value and let the next row be the most similar one to the previous row
        # best_row = 0
        # best_value = -1
        # if verbose:
        #     print("Sort parameters according to similarity matrix")
        # for i in range(len(in_params_tmp)):
        #     for j in range(len(in_params_tmp)):
        #         if i >= j:
        #             continue
        #         if best_value == -1:
        #             best_value = similarity_matrix[i, j]
        #             best_row = i
        #         elif best_value > similarity_matrix[i, j]:
        #             best_value = similarity_matrix[i, j]
        #             best_row = i
        # in_params_sorted = [in_params_tmp[best_row]]
        # previous_rows = [best_row]
        # if verbose:
        #     print("Generate matrix for plotting")
        # for i in range(len(in_params_tmp)-1):
        #     next_row = 0
        #     next_best_value = -1
        #     for j in range(len(in_params_tmp)):
        #         if j in previous_rows:
        #             continue
        #         if next_best_value == -1:
        #             next_best_value = similarity_matrix[i, j]
        #             next_row = j
        #         elif next_best_value > similarity_matrix[i, j]:
        #             next_best_value = similarity_matrix[i, j]
        #             next_row = j
        #     in_params_sorted.append(in_params_tmp[next_row])
        #     previous_rows.append(next_row)

        # found the correct order. Let's create a corresponding 2D-array and plot it
        n_edges = len(edges_in_params[out_p][in_params_tmp[0]][:-1])
        plot_matrix = np.zeros((len(in_params_tmp), n_edges))
        for i, in_p in enumerate(in_params_sorted):
            plot_matrix[i, :] = hist_in_params[out_p][in_p]
        plot_matrix[plot_matrix == 0] = np.nan
        ax = sns.heatmap(
            plot_matrix,
            cbar=True,
            cmap="viridis",
            norm=mpl_col.LogNorm(),
            yticklabels=in_params_sorted,
        )

        if title is None:
            title2 = f"Histogram for {out_p}"
        else:
            title2 = title
        _ = ax.set_title(title2)
        plt.tight_layout()
        fig = ax.get_figure()
        i = 0
        store_type = filename.split(".")[-1]
        store_path = filename[0 : -len(store_type) - 1]
        save = store_path + "_" + out_p + "_{:03d}.".format(i) + store_type
        while os.path.isfile(save):
            i = i + 1
            save = store_path + "_" + out_p + "_{:03d}.".format(i) + store_type
        fig.savefig(save, dpi=300)
        plt.clf()


def get_sums(
    file_path,
    store_path=None,
    only_asc600=False,
    inoutflow_time=-1,
):
    """
    Calculate the sum for all gradients over all time steps and trajectories for multiple files.

    Parameters
    ----------
    file_path : string
        Path to NetCDF-files with trajectories from a sensitivity analysis simulation.
    store_path : string
        If a store path is given then dumps the sums to 'f{store_path}_sums.pkl'.
    only_asc600 : bool
        Consider only time steps during the ascend.
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.

    Returns
    -------
    Dictionary with model state variables as keys and the sums for each gradients as pandas.DataFrame
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    out_param_coord = "Output_Parameter_ID"
    if "Output_Parameter_ID" not in ds:
        out_param_coord = "Output Parameter"
        out_params = list(ds[out_param_coord].values)
        param_name = out_params
    else:
        out_params = list(ds[out_param_coord].values)
        param_name = []
        for out_p in out_params:
            param_name.append(latexify.param_id_map[out_p])
    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    sums = {}
    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
        )
        ds[in_params] = np.abs(ds[in_params])
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            df = (
                ds[in_params]
                .sel({out_param_coord: out_p})
                .sum(dim=["trajectory", "time"], skipna=True)
                .to_dataframe()
                .reset_index()
            )
            df = df[in_params]

            if out_name in sums.keys():
                sums[out_name] += df
            else:
                sums[out_name] = df
    if store_path is not None and store_path != "no":
        with open(store_path + "_sums.pkl", "wb") as f:
            pickle.dump(sums, f)
    return sums


def get_sums_phase(
    file_path,
    store_path=None,
    only_asc600=False,
    inoutflow_time=-1,
):
    """
    Calculate the sum for all gradients over all time steps and trajectories for each phase.

    Parameters
    ----------
    file_path : string
        Path to NetCDF-files with trajectories from a sensitivity analysis simulation.
    store_path : string
        If a store path is given then dumps the sums to 'f{store_path}_sums_phase.pkl'.
    only_asc600 : bool
        Consider only time steps during the ascend.
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.

    Returns
    -------
    Dictionary with phase + model state variables as keys and the sums for each gradients as pandas.DataFrame
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    out_param_coord = "Output_Parameter_ID"
    if out_param_coord not in ds:
        out_param_coord = "Output Parameter"
        out_params = list(ds[out_param_coord].values)
        param_name = out_params
    else:
        param_name = []
        out_params = list(ds[out_param_coord].values)
        for idx in out_params:
            param_name.append(latexify.param_id_map[idx])

    phases = ["warm phase", "mixed phase", "ice phase", "neutral phase"]
    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    sums = {}
    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
        )
        ds[in_params] = np.abs(ds[in_params])
        if ds["phase"].dtype != str and ds["phase"].dtype != np.uint64:
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({out_param_coord: out_p})
            for phase_i, phase in enumerate(phases):
                if ds_tmp["phase"].dtype == str:
                    idx = np.where(ds_tmp["phase"] == phase)
                else:
                    idx = np.where(ds_tmp["phase"] == phase_i)
                dic = {in_p: None for in_p in in_params}
                for in_p in in_params:
                    dic[in_p] = [np.nansum(ds_tmp[in_p].values[idx])]
                if phase + " " + out_name in sums.keys():
                    sums[phase + " " + out_name] += pd.DataFrame(dic)
                else:
                    sums[phase + " " + out_name] = pd.DataFrame(dic)
    if store_path is not None and store_path != "no":
        with open(store_path + "_sums_phase.pkl", "wb") as f:
            pickle.dump(sums, f)
    return sums


def get_cov_matrix(
    file_path,
    in_params=None,
    store_path=None,
    only_asc600=False,
    inoutflow_time=-1,
):
    """
    Calculate the means and covariance matrices for model state variables and sensitivities. For each model state
    variable with sensitivities there will be a different covariance matrix.

    Parameters
    ----------
    file_path : string
        Path to NetCDF-files with trajectories from a sensitivity analysis simulation.
    in_params : list-like of strings
        List of model parameters.
    store_path : string
        If a store path is given then dumps the means to 'f{store_path}_means.pkl'
        and the covariance matrix to f'{store_path}_covariance_matrix.pkl'.
    only_asc600 : bool
        Consider only time steps during the ascend.
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.

    Returns
    -------
    means, cov, where means is a dictionary of model state variables with available sensitivities and values are
    dictionaries with model parameters/state variables as keys and the mean values as values, i.e., the mean value
    of dQV/db_v is means['QV']['db_v'], the mean value of QV is means[any key]['QV'].
    cov is a dictionary of covariance matrices. The keys are model state variables with available sensitivities.
    The values are matrices where covariances are given for sensitivities for the respective model state variable.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    if in_params is None:
        in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    out_param_coord = "Output_Parameter_ID"
    if out_param_coord not in ds:
        out_param_coord = "Output Parameter"
        out_params = list(ds[out_param_coord].values)
        param_name = out_params
    else:
        param_name = []
        out_params = list(ds[out_param_coord].values)
        for idx in out_params:
            param_name.append(latexify.param_id_map[idx.values])
    more_params = []
    if only_asc600 or inoutflow_time > 0:
        more_params.append("asc600")
    all_params = param_name + in_params
    n = len(all_params)
    means = {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
    cov = {out_p: np.zeros((n, n), dtype=np.float64) for out_p in param_name}
    n_total = {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
        )
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({out_param_coord: out_p})
            means_tmp = ds_tmp.mean(skipna=True)
            count_tmp = (~np.isnan(ds_tmp)).sum()
            for p in means[out_name]:
                if n_total[out_name][p] > 0:
                    n_new = count_tmp[p].values.item() + n_total[out_name][p]
                    if n_new > 0:
                        means[out_name][p] = (
                            means[out_name][p] * n_total[out_name][p] / n_new
                            + means_tmp[p].values.item()
                            * count_tmp[p].values.item()
                            / n_new
                        )
                        n_total[out_name][p] = n_new
                else:
                    n_total[out_name][p] = count_tmp[p].values.item()
                    means[out_name][p] = means_tmp[p].values.item()

    n_total = {out_p: {p: 0.0 for p in all_params} for out_p in param_name}
    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
        )
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({out_param_coord: out_p})
            count_tmp = (~np.isnan(ds_tmp)).sum()
            for i, p in enumerate(means[out_name]):
                if n_total[out_name][p] > 0:
                    n_add = count_tmp[p].values.item()
                    n_new = n_add + n_total[out_name][p]
                    if n_new > 0:
                        for i2, p2 in enumerate(means[out_name]):
                            cov_tmp = np.nanmean(
                                (ds_tmp[p].values - means[out_name][p])
                                * (ds_tmp[p2].values - means[out_name][p2])
                            )
                            cov[out_name][i, i2] = (
                                cov[out_name][i, i2] * n_total[out_name][p] / n_new
                                + cov_tmp * n_add / n_new
                            )
                        n_total[out_name][p] = n_new
                else:
                    n_new = count_tmp[p].values.item() + n_total[out_name][p]
                    if n_new > 0:
                        for i2, p2 in enumerate(means[out_name]):
                            cov[out_name][i, i2] = np.nanmean(
                                (ds_tmp[p].values - means[out_name][p])
                                * (ds_tmp[p2].values - means[out_name][p2])
                            )
                        n_total[out_name][p] = n_new
    if store_path is not None and store_path != "no":
        with open(store_path + "_means.pkl", "wb") as f:
            pickle.dump(means, f)
        with open(store_path + "_covariance_matrix.pkl", "wb") as f:
            pickle.dump(cov, f)
    return means, cov


def get_cov_matrix_phase(
    file_path,
    in_params=None,
    store_path=None,
    only_asc600=False,
    inoutflow_time=-1,
):
    """
    Calculate the means and covariance matrices for model state variables and sensitivities. For each model state
    variable with sensitivities and each phase, there will be a different covariance matrix.

    Parameters
    ----------
    file_path : string
        Path to NetCDF-files with trajectories from a sensitivity analysis simulation.
    in_params : list-like of strings
        List of model parameters.
    store_path : string
        If a store path is given then dumps the means to 'f{store_path}_means_phases.pkl'
        and the covariance matrix to f'{store_path}_covariance_matrix_phases.pkl'.
    only_asc600 : bool
        Consider only time steps during the ascend.
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.

    Returns
    -------
    means, cov, where means is a dictionary of phases and then of model state variables with available sensitivities
    and values are dictionaries with model parameters/state variables as keys and the mean values as values, i.e.,
    the mean value of dQV/db_v during the warm phase is means['warm phase']['QV']['db_v'],
    the mean value of QV during the warm phase is means['warm phase'][any key]['QV'].
    cov is a dictionary of covariance matrices. The first keys are phases, the second keys are model state variables
    with available sensitivities.
    The values are matrices where covariances are given for sensitivities for the respective model state variable.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    if in_params is None:
        in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    out_param_coord = "Output_Parameter_ID"
    if out_param_coord not in ds:
        out_param_coord = "Output Parameter"
        out_params = list(ds[out_param_coord].values)
        param_name = out_params
    else:
        param_name = []
        out_params = list(ds[out_param_coord].values)
        for idx in out_params:
            param_name.append(latexify.param_id_map[idx.values])
    more_params = ["phase"]
    if only_asc600 or inoutflow_time > 0:
        more_params.append("asc600")
    phases = ["warm phase", "mixed phase", "ice phase", "neutral phase"]
    all_params = param_name + in_params
    n = len(all_params)
    means = {
        phase: {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
        for phase in phases
    }
    cov = {
        phase: {out_p: np.zeros((n, n), dtype=np.float64) for out_p in param_name}
        for phase in phases
    }
    n_total = {
        phase: {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
        for phase in phases
    }
    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            load_params=all_params + more_params,
        )
        if ds["phase"].dtype != str and ds["phase"].dtype != np.uint64:
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({out_param_coord: out_p})

            for phase_i, phase in enumerate(tqdm(phases, leave=False)):
                if ds_tmp["phase"].dtype == str:
                    idx = np.where(ds_tmp["phase"] == phase)
                else:
                    idx = np.where(ds_tmp["phase"] == phase_i)
                for p in tqdm(means[phase][out_name], leave=False):
                    count_tmp = np.sum((~np.isnan(ds_tmp[p].values[idx])))
                    if count_tmp == 0:
                        continue
                    mean_tmp = np.nanmean(ds_tmp[p].values[idx])

                    if n_total[phase][out_name][p] > 0:
                        n_new = count_tmp + n_total[phase][out_name][p]
                        if n_new > 0:
                            means[phase][out_name][p] = (
                                means[phase][out_name][p]
                                * n_total[phase][out_name][p]
                                / n_new
                                + mean_tmp * count_tmp / n_new
                            )
                            n_total[phase][out_name][p] = n_new
                    else:
                        n_total[phase][out_name][p] = count_tmp
                        means[phase][out_name][p] = mean_tmp

    n_total = {
        phase: {out_p: {p: 0.0 for p in all_params} for out_p in param_name}
        for phase in phases
    }
    for f in tqdm(files):
        ds = load_ds(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            load_params=all_params + more_params,
        )
        if ds["phase"].dtype != str and ds["phase"].dtype != np.uint64:
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({out_param_coord: out_p})
            for phase_i, phase in enumerate(tqdm(phases, leave=False)):
                if ds_tmp["phase"].dtype == str:
                    idx = np.where(ds_tmp["phase"] == phase)
                else:
                    idx = np.where(ds_tmp["phase"] == phase_i)
                for i, p in enumerate(tqdm(means[phase][out_name], leave=False)):
                    count_tmp = np.sum((~np.isnan(ds_tmp[p].values[idx])))
                    if count_tmp == 0:
                        continue
                    p_variance = ds_tmp[p].values[idx] - means[phase][out_name][p]
                    if n_total[phase][out_name][p] > 0:
                        n_add = count_tmp
                        n_new = n_add + n_total[phase][out_name][p]
                        if n_new > 0 and n_add > 0:
                            for i2, p2 in enumerate(means[phase][out_name]):
                                cov_tmp = np.nanmean(
                                    p_variance
                                    * (
                                        ds_tmp[p2].values[idx]
                                        - means[phase][out_name][p2]
                                    )
                                )
                                cov[phase][out_name][i, i2] = (
                                    cov[phase][out_name][i, i2]
                                    * n_total[phase][out_name][p]
                                    / n_new
                                    + cov_tmp * n_add / n_new
                                )
                            n_total[phase][out_name][p] = n_new
                    else:
                        n_new = count_tmp + n_total[phase][out_name][p]
                        if n_new > 0:
                            for i2, p2 in enumerate(means[phase][out_name]):
                                cov[phase][out_name][i, i2] = np.nanmean(
                                    p_variance
                                    * (
                                        ds_tmp[p2].values[idx]
                                        - means[phase][out_name][p2]
                                    )
                                )
                            n_total[phase][out_name][p] = n_new
    if store_path is not None and store_path != "no":
        with open(store_path + "_means_phases.pkl", "wb") as f:
            pickle.dump(means, f)
        with open(store_path + "_covariance_matrix_phases.pkl", "wb") as f:
            pickle.dump(cov, f)
    return means, cov


def plot_heatmap_cov(
    data_in,
    out_param,
    names,
    in_params=None,
    plot_type="all",
    norm=None,
    title=None,
    filename=None,
    width=30,
    height=16,
):
    """
    Plot a heatmap of a covariance matrix.

    Parameters
    ----------
    data_in : dictionary of matrices
        The keys are model state variables with available sensitivities.
        The values are matrices where covariances are given for sensitivities for the respective model state variable.
        The data can be generated using get_cov_matrix().
    out_param : string
        The model state variable for which the sensitivities shall be plotted for.
    names : list of strings
        The names of each column/row in the covariance matrix. The index of names is the row/column of the
        covariance matrix.
    in_params : list-like of strings
        List of model parameters or model states to plot in the covariance matrix.
    plot_type : string
        Define if only negative (='negative'), positive (='positive'), or all ('all') values shall be plotted.
        If 'all' is chosen and a norm other than SymLogNorm is given, the absolute values will be plotted.
    norm : matplotlib.colors normalization instance
        A normalization for the colormap, such as matplotlib.colors.SymLogNorm()
    title : string
        Optional title for the plot. Otherwise a standard name will be used.
    filename : string
        Path and filename to save the plot on disk. The filename will be numerated.
        If a file with the same name already exists, the number will be incremented.
    width : float
        Width in inches
    height : float
        Height in inches

    Returns
    -------
    If successful, returns matplotlib.axes. Otherwise returns None.
    """
    sns.set(rc={"figure.figsize": (width, height)})
    data = copy.deepcopy(data_in[out_param])
    if in_params is not None and len(in_params) < np.shape(data_in[out_param])[0]:
        if len(in_params) == 0:
            return
        idx = []
        names = np.asarray(list(names))
        for in_p in in_params:
            idx.append(np.where(names == in_p)[0][0])
        idx = np.asarray(idx)
        data = data[idx]
        data = data[:, idx]
        names = names[idx]
    if plot_type == "negative":
        data[np.where(data >= 0)] = np.nan
    elif plot_type == "positive":
        data[np.where(data < 0)] = np.nan
    if (
        norm is not None
        and plot_type != "positive"
        and not isinstance(norm, mpl_col.SymLogNorm)
    ):
        data = np.abs(data)
    try:
        g = sns.heatmap(
            data=data,
            cmap="viridis",
            norm=norm,
            cbar=True,
            yticklabels=names,
            xticklabels=names,
            annot=((len(names) <= 20) and (width > 19)),
            fmt="1.2e",
        )
        if title is None:
            title_ = f"Heatmap of Covariance (Gradients for {out_param})"
        else:
            title_ = title
        _ = g.set_title(title_)
        plt.tight_layout()
        fig = g.get_figure()
        if filename is not None:
            i = 0
            store_type = filename.split(".")[-1]
            store_path = filename[0 : -len(store_type) - 1]
            save = store_path + "_{:03d}.".format(i) + store_type
            while os.path.isfile(save):
                i = i + 1
                save = store_path + "_{:03d}.".format(i) + store_type
            fig.savefig(save, dpi=300)
        return g
    except:
        if filename is not None:
            print(f"Cannot create plot for {filename} ({out_param})")
        else:
            print(f"Cannot create plot for {out_param}")
        return None


def plot_heatmap_histogram(
    hist_conditional,
    in_params,
    out_params,
    conditions,
    title=None,
    filename=None,
    width=17,
    height=16,
    log=True,
    font_scale=None,
    save=True,
    interactive=False,
    latex=False,
    ticks_offset=None,
    verbose=False,
):
    """
    Plot 2D histogram.

    Parameters
    ----------
    hist_conditional : Dictionary with edges and bin counts
        The dictionary generated using get_histogram_cond(). It has the following keys:
        'edges_out_params': Dictionary where the keys are model state variables and the values are arrays of bin edges.
        'edges_in_params': Dictionary where the keys are model state variables for which sensitivities are available
            and the values are dictionaries of model parameters with arrays of bin edges.
        model state variables: Each model state variable has a dictionary for 'hist_out_params' and 'hist_in_params'.
        'hist_out_params' is a dictionary of model state variables with arrays of bin counts.
        'hist_in_params' is a dictionary of model state variables for which sensitivities are available
            and the values are dictionaries of model parameters with arrays of bin counts.
    in_params : list of strings
        A list of model parameters.
    out_params : list of strings
        A list of model state variables with available sensitivities.
    conditions : list of strings
        A list of model state variables for the x-axis.
    title : string
        Title of the histogram. If none is given, a title will be generated.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    width : float
        Width in inches
    height : float
        Height in inches
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    log : bool
        Plot the histograms using a log-scale
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    interactive : bool
        Create a figure for interactive plotting.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.
    ticks_offset : int
        If none, the number of ticks is calculated based on the width of the plot. Otherwise,
        every other "ticks_offset" is used.
    verbose : bool
        Print some additional information.

    Returns
    -------
    If filename is given, returns None. If filename is None, returns the matplotlib.figure.Figure with the plot drawn
     onto it.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    if interactive:
        sns.set()
        fig = Figure()
        ax = fig.subplots()
    else:
        fig = None
        ax = None
    norm = None
    if log:
        norm = mpl_col.LogNorm()

    def plot_hist(
        hist2d, x_name, y_name, x_ticks, y_ticks, title=None, p=None, ax=None
    ):
        if title is None:
            if p is None:
                title = f"Histogram for {latexify.parse_word(y_name)} over {latexify.parse_word(x_name)}"
            else:
                p_name = r"$ \partial$" + latexify.parse_word(p)
                title = f"Histogram for {p_name}/{latexify.parse_word(y_name)} over {latexify.parse_word(x_name)}"

        if verbose:
            if p is None:
                print(f"Plotting histogram for {x_name}, {y_name}")
            else:
                print(f"Plotting histogram for {x_name}, d {p}/{y_name}")
        try:
            if ax is None:
                ax = sns.heatmap(
                    np.transpose(
                        hist2d
                    ),  # np.histogram2d uses the first dimension for the x-axis
                    cmap="viridis",
                    cbar=True,
                    square=True,
                    norm=norm,
                )
            else:
                _ = sns.heatmap(
                    np.transpose(
                        hist2d
                    ),  # np.histogram2d uses the first dimension for the x-axis
                    cmap="viridis",
                    cbar=True,
                    square=True,
                    norm=norm,
                    ax=ax,
                )
        except:
            print(f"Plotting for {x_name}, {y_name} failed")
            return
        if ticks_offset is None:
            if 24 > width > 5:
                local_ticks_offset = 24 // (width - 5)
            elif width <= 5:
                local_ticks_offset = 6
        else:
            local_ticks_offset = ticks_offset

        if local_ticks_offset >= np.shape(hist2d)[0]:
            local_ticks_offset_x = np.shape(hist2d)[0] - 1
        else:
            local_ticks_offset_x = local_ticks_offset
        x_ticks_location = np.arange(0, np.shape(hist2d)[0], local_ticks_offset_x)
        ax.set_xticks(x_ticks_location)
        if local_ticks_offset >= np.shape(hist2d)[1]:
            local_ticks_offset_y = np.shape(hist2d)[1] - 1
        else:
            local_ticks_offset_y = local_ticks_offset
        y_ticks_location = np.arange(0, np.shape(hist2d)[1], local_ticks_offset_y)
        ax.set_yticks(y_ticks_location)
        x_labels = [f"{x_ticks[i]:1.1e}" for i in x_ticks_location]
        _ = ax.set_xticklabels(
            x_labels, rotation=45, ha="right", rotation_mode="anchor"
        )
        y_labels = [f"{y_ticks[i]:1.1e}" for i in y_ticks_location]
        _ = ax.set_yticklabels(y_labels, rotation=0)
        if font_scale is None:
            _ = ax.set_title(title)
        else:
            ax.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
            _ = ax.set_title(title, fontsize=int(12 * font_scale))
            cbar = ax.collections[-1].colorbar
            cbarax = cbar.ax
            cbarax.tick_params(labelsize=int(10 * font_scale))
        ax.set_xlabel(latexify.parse_word(x_name), fontsize=int(11 * font_scale))
        ax.set_ylabel(latexify.parse_word(y_name), fontsize=int(11 * font_scale))

        plt.tight_layout()
        if filename is not None and save:
            fig = ax.get_figure()
            try:
                i = 0
                store_type = filename.split(".")[-1]
                store_path = filename[0 : -len(store_type) - 1]
                if p is None:
                    save_name = (
                        store_path
                        + "_"
                        + x_name
                        + "_"
                        + y_name
                        + "_{:03d}.".format(i)
                        + store_type
                    )
                else:
                    save_name = (
                        store_path
                        + "_"
                        + x_name
                        + "_"
                        + p
                        + "_wrt_"
                        + y_name
                        + "_{:03d}.".format(i)
                        + store_type
                    )
                while os.path.isfile(save_name):
                    i = i + 1
                    if p is None:
                        save_name = (
                            store_path
                            + "_"
                            + x_name
                            + "_"
                            + y_name
                            + "_{:03d}.".format(i)
                            + store_type
                        )
                    else:
                        save_name = (
                            store_path
                            + "_"
                            + x_name
                            + "_"
                            + p
                            + "_wrt_"
                            + y_name
                            + "_{:03d}.".format(i)
                            + store_type
                        )
                fig.savefig(save_name, bbox_inches="tight", dpi=300)
            except:
                print(f"Storing to {save_name} failed.", file=sys.stderr)
            if not interactive:
                plt.clf()

    if interactive:
        c = conditions
        if in_params is None or in_params == "None":
            out_p = out_params
            plot_hist(
                hist_conditional[c]["hist_out_params"][out_p],
                c,
                out_p,
                hist_conditional["edges_out_params"][c][:-1],
                hist_conditional["edges_out_params"][out_p][:-1],
                title,
                ax=ax,
            )
        else:
            p = out_params
            in_p = in_params
            plot_hist(
                hist_conditional[c]["hist_in_params"][p][in_p],
                c,
                in_p,
                hist_conditional["edges_out_params"][c][:-1],
                hist_conditional["edges_in_params"][p][in_p][:-1],
                title,
                p,
                ax=ax,
            )
        save = False
        filename = None
        return fig

    for c in tqdm(conditions):
        for out_p in tqdm(out_params, leave=False):
            if c == out_p:
                continue
            plot_hist(
                hist_conditional[c]["hist_out_params"][out_p],
                c,
                out_p,
                hist_conditional["edges_out_params"][c][:-1],
                hist_conditional["edges_out_params"][out_p][:-1],
                title,
            )
        # Sensitivities are p w.r.t. in_p
        for p in tqdm(list(hist_conditional["edges_in_params"]), leave=False):
            for in_p in tqdm(in_params, leave=False):
                if in_p not in hist_conditional[c]["hist_in_params"][p]:
                    continue
                plot_hist(
                    hist_conditional[c]["hist_in_params"][p][in_p],
                    c,
                    in_p,
                    hist_conditional["edges_out_params"][c][:-1],
                    hist_conditional["edges_in_params"][p][in_p][:-1],
                    title,
                    p,
                )

    if verbose:
        print("All plots finished!")
    return None


def plot_traj_histogram_out_interactive(edges, hist):
    """
    Calling this function from a Jupyter notebook allows to visualize the
    traj_plot_histogram_out interactively. Plots the histogram of an output parameter,
    i.e., QV, latent_heat, latent_cool, etc.

    Parameters
    ----------
    edges : Dictionary of list-like of float
        Edges for the histogram. Keys must be in out_params.
    hist : Dictionary of list-like of int
        Number of entries for each bin. Keys must be in edges.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=list(edges.keys())[0],
        options=list(edges.keys()),
    )
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
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    log_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=True,
        button_type="success",
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    tick_slider = pn.widgets.IntSlider(
        name="Plot every n ticks on the x-axis:",
        start=1,
        end=20,
        step=1,
        value=6,
    )
    plot_pane = pn.panel(
        pn.bind(
            traj_plot_histogram_out,
            out_params=out_param,
            filename=save_to_field,
            edges=edges,
            hist=hist,
            log=log_plot,
            width=width_slider,
            height=height_slider,
            title=title_widget,
            save=save_button,
            interactive=True,
            font_scale=font_slider,
            latex=latex_button,
            ticks_offset=tick_slider,
            verbose=False,
        ),
    ).servable()

    return pn.Column(
        out_param,
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
            tick_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            log_plot,
            latex_button,
        ),
        title_widget,
        plot_pane,
    )


def plot_traj_histogram_inp_interactive(
    edges_in_params,
    hist_in_params,
):
    """
    Can be used in jupyter notebooks to interactively plot traj_plot_histogram_inp(). From traj_plot_histogram_inp():
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot three histograms per image with
    \partial output / \partial model parameter
    where output is QV, latent_heat and latent_cool. Plot one image per model_parameter.

    Parameters
    ----------
    edges_in_params : Dictionary of dictionary list-like of float
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        the bin edges for the given keys. Optional: The first level can be another dictionary of phases.
    hist_in_params : Dictionary of dictionary list-like of int
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        values of the histogram for the given key. Optional: The first level can be another dictionary of phases.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    tmp_params = list(edges_in_params[list(edges_in_params.keys())[0]].keys())
    in_params = []
    for param in tmp_params:
        if param[0] == "d" and param != "deposition":
            in_params.append(param)

    in_param = pn.widgets.Select(
        name="Model Parameter",
        value=in_params[0],
        options=in_params,
    )
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
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    log_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=True,
        button_type="success",
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    tick_slider = pn.widgets.IntSlider(
        name="Plot every n ticks on the x-axis:",
        start=1,
        end=20,
        step=1,
        value=6,
    )

    plot_pane = pn.panel(
        pn.bind(
            traj_plot_histogram_inp,
            filename=save_to_field,
            in_params=in_param,
            edges_in_params=edges_in_params,
            hist_in_params=hist_in_params,
            log=log_plot,
            width=width_slider,
            height=height_slider,
            title=None,
            font_scale=font_slider,
            save=save_button,
            interactive=True,
            latex=latex_button,
            ticks_offset=tick_slider,
            verbose=False,
        ),
    ).servable()

    return pn.Column(
        in_param,
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
            tick_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            log_plot,
            latex_button,
        ),
        plot_pane,
    )


def plot_heatmap_histogram_interactive(hist_conditional):
    """
    Can be used in jupyter notebooks to interactively plot plot_heatmap_histogram() (2D histograms).

    Parameters
    ----------
    hist_conditional : Dictionary of dictionaries with edges and bin counts for 2D histograms.
        Result of get_histogram_cond().
        Dictionary with the following keys:
        'edges_out_params': Dictionary where the keys are model state variables and the values are arrays of bin edges.
        'edges_in_params': Dictionary where the keys are model state variables for which sensitivities are available
            and the values are dictionaries of model parameters with arrays of bin edges.
        model state variables: Each model state variable has a dictionary for 'hist_out_params' and 'hist_in_params'.
        'hist_out_params' is a dictionary of model state variables with arrays of bin counts.
        'hist_in_params' is a dictionary of model state variables for which sensitivities are available
            and the values are dictionaries of model parameters with arrays of bin counts.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    conds = list(hist_conditional.keys())
    conditions = []
    for c in conds:
        if c != "edges_in_params" and c != "edges_out_params":
            conditions.append(c)
    wrt_params = list(hist_conditional["edges_in_params"].keys())
    out_param = pn.widgets.RadioButtonGroup(
        name="Output Parameter",
        value=wrt_params[0],
        options=wrt_params,
        button_type="primary",
    )
    condition = pn.widgets.Select(
        name="X-Axis",
        value=conditions[0],
        options=conditions,
    )
    in_params = ["None"]
    tmp_params = list(hist_conditional["edges_in_params"][wrt_params[0]].keys())
    for param in tmp_params:
        if param[0] == "d" and param != "deposition":
            in_params.append(param)
    in_param = pn.widgets.Select(
        name="Model Parameter (Y-Axis)",
        value=in_params[-1],
        options=in_params,
    )
    width_slider = pn.widgets.IntSlider(
        name="Width in inches",
        start=3,
        end=20,
        step=1,
        value=8,
    )
    height_slider = pn.widgets.IntSlider(
        name="Height in inches",
        start=3,
        end=20,
        step=1,
        value=8,
    )
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    log_plot = pn.widgets.Toggle(
        name="Use log colorbar",
        value=True,
        button_type="success",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    tick_slider = pn.widgets.IntSlider(
        name="Plot every n ticks:",
        start=1,
        end=20,
        step=1,
        value=6,
    )

    plot_pane = pn.panel(
        pn.bind(
            plot_heatmap_histogram,
            hist_conditional=hist_conditional,
            filename=save_to_field,
            in_params=in_param,
            out_params=out_param,
            conditions=condition,
            log=log_plot,
            width=width_slider,
            height=height_slider,
            font_scale=font_slider,
            title=title_widget,
            save=save_button,
            latex=latex_button,
            interactive=True,
            ticks_offset=tick_slider,
            verbose=False,
        ),
    ).servable()

    return pn.Column(
        in_param,
        "w.r.t. Model State (Y-Axis)",
        out_param,
        condition,
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
            tick_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            log_plot,
            latex_button,
        ),
        title_widget,
        plot_pane,
    )


def traj_plot_kde_inp(
    in_params,
    out_param,
    filename,
    edges_in_params,
    hist_in_params,
    linewidth=2,
    bw_adjust=0.1,
    log=False,
    font_scale=1,
    width=24,
    height=12,
    title=None,
    save=False,
    latex=False,
):
    """
    Similar to traj_plot_histogram_inp() but using a kde estimation.

    Parameters
    ----------
    in_params : list-like of strings
        Input parameter names to plot the histogram for.
    out_param : string
        Name od model state variable to plot the sensitivities for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges_in_params : Dictionary of dictionary list-like of float
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        the bin edges for the given keys. Optional: The first level can be another dictionary of phases.
    hist_in_params : Dictionary of dictionary list-like of int
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        values of the histogram for the given key. Optional: The first level can be another dictionary of phases.
    linewidth : float
        Line width of each kde.
    bw_adjust : float
        Used to calculate the kde. Typically, 0.3 is used but here we set the default lower since the true
        distributions tend to be less smooth. Larger values generate smoother estimations.
    log : bool
        Plot the y-axis using log-scale.
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.

    Returns
    -------
    Returns the matplotlib.figure.Figure with the plot drawn onto it.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    common_norm = False
    data_dic = {
        "Parameter": [],
        "weight": [],
        "Impact": [],
    }

    if len(in_params) == 0:
        return None

    for in_p in in_params:
        if log:
            zeros = np.argwhere(hist_in_params[out_param][in_p] == 0)
            weights = np.log10(
                hist_in_params[out_param][in_p],
                where=(hist_in_params[out_param][in_p] >= 0.5),
            )
            weights += 0.1
            weights[zeros] = 0
        else:
            weights = np.ones(np.shape(hist_in_params[out_param][in_p]))
        data_dic["weight"].extend(weights)
        data_dic["Impact"].extend(
            edges_in_params[out_param][in_p][:-1]
            + (
                edges_in_params[out_param][in_p][1]
                - edges_in_params[out_param][in_p][0]
            )
            / 2
        )
        data_dic["Parameter"].extend(np.repeat(in_p, len(weights)))

    df = pd.DataFrame(data_dic)

    fig = Figure()
    ax = fig.subplots()

    _ = sns.kdeplot(
        data=df,
        x="Impact",
        hue="Parameter",
        weights="weight",
        common_norm=common_norm,
        ax=ax,
        bw_adjust=bw_adjust,
        hue_order=np.sort(in_params),
        linewidth=linewidth,
    )

    ax.xaxis.get_label().set_fontsize(int(11 * font_scale))
    ax.yaxis.get_label().set_fontsize(int(11 * font_scale))
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=int(10 * font_scale),
    )
    if title is not None:
        _ = ax.set_title(title, fontsize=int(12 * font_scale))

    plt.setp(ax.get_legend().get_texts(), fontsize=int(9 * font_scale))
    plt.setp(ax.get_legend().get_title(), fontsize=int(11 * font_scale))
    ax.xaxis.get_offset_text().set_size(int(9 * font_scale))
    ax.yaxis.get_offset_text().set_size(int(9 * font_scale))

    if filename is not None and save:
        fig = ax.get_figure()
        try:
            i = 0
            store_type = filename.split(".")[-1]
            store_path = filename[0 : -len(store_type) - 1]
            save_name = store_path + "_{:03d}.".format(i) + store_type

            while os.path.isfile(save_name):
                i = i + 1
                save_name = store_path + "_{:03d}.".format(i) + store_type
            fig.savefig(save_name, bbox_inches="tight", dpi=300)
        except:
            print(f"Storing to {save_name} failed.", file=sys.stderr)

    return fig


def plot_traj_kde_inp_interactive(
    edges_in_params,
    hist_in_params,
):
    """
    Can be used in jupyter notebooks to interactively plot traj_plot_kde_inp().

    Parameters
    ----------
    edges_in_params : Dictionary of dictionary list-like of float
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        the bin edges for the given keys. Optional: The first level can be another dictionary of phases.
    hist_in_params : Dictionary of dictionary list-like of int
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        values of the histogram for the given key. Optional: The first level can be another dictionary of phases.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    tmp_params = list(edges_in_params[list(edges_in_params.keys())[0]].keys())
    in_params = []
    for param in tmp_params:
        if param[0] == "d" and param != "deposition":
            in_params.append(param)

    in_param = pn.widgets.CrossSelector(
        name="Parameter",
        value=in_params[0:2],
        options=in_params,
    )
    out_param_list = list(edges_in_params.keys())
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_param_list[0],
        options=out_param_list,
    )
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
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    line_slider = pn.widgets.FloatSlider(
        name="Change the line width",
        start=1,
        end=10,
        step=0.5,
        value=2,
    )
    bw_slider = pn.widgets.FloatSlider(
        name="Change the bandwidth for the kde calculation",
        start=0.05,
        end=1.5,
        step=0.05,
        value=0.1,
    )
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )

    plot_pane = pn.panel(
        pn.bind(
            traj_plot_kde_inp,
            filename=save_to_field,
            in_params=in_param,
            out_param=out_param,
            linewidth=line_slider,
            bw_adjust=bw_slider,
            edges_in_params=edges_in_params,
            hist_in_params=hist_in_params,
            log=True,
            width=width_slider,
            height=height_slider,
            title=title_widget,
            font_scale=font_slider,
            save=save_button,
            latex=latex_button,
        ),
    ).servable()

    return pn.Column(
        pn.Row(
            in_param,
            out_param,
        ),
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            latex_button,
        ),
        pn.Row(
            line_slider,
            bw_slider,
            title_widget,
        ),
        plot_pane,
    )


def main(args):
    """
    Handle the plotting routines for different data inputs and plots.

    Parameters
    ----------
    args : arparse.ArgumentParser or any class with the same members
        A number of arguments to handle loading and plotting.
    """
    if not args.from_processed:
        if (
            args.plot_type == "all"
            or args.plot_type == "hist_out"
            or args.plot_type == "hist_in"
            or args.plot_type == "heat"
            or args.plot_type == "2D_hist"
            or args.save_histogram != "no"
        ):
            if args.load_histogram != "no":
                print("########### Load histograms ###########")
                file_path = args.load_histogram
                if len(args.conditional_hist) == 0:
                    with open(file_path + "_all_hist.pkl", "rb") as f:
                        all_hist = pickle.load(f)
                    hist = all_hist["hist_out_params"]
                    hist_in_params = all_hist["hist_in_params"]
                    edges = all_hist["edges_out_params"]
                    edges_in_params = all_hist["edges_in_params"]
                else:
                    with open(file_path + "_hist_conditional.pkl", "rb") as f:
                        hist_conditional = pickle.load(f)
            else:
                print("########### Calculate histograms ###########")
                if args.filter > 0 and args.load_covariance == "no":
                    print(
                        "If you want to filter sensitivities for the histogram "
                        "then you have to provide a path to the means using "
                        "'--load_covariance'"
                    )
                    print(
                        "You can generate this file using --plot_type cov_heat."
                        "This generates a covariance matrix for the parameters with "
                        "the largest magnitudes."
                    )
                    return
                elif args.filter:
                    with open(args.load_covariance + "_means.pkl", "rb") as f:
                        means = pickle.load(f)
                else:
                    means = None
                if len(args.conditional_hist) == 0:
                    all_hist = get_histogram(
                        args.file,
                        additional_params=args.additional_hist_params,
                        only_asc600=args.only_asc600,
                        inoutflow_time=args.inoutflow_time,
                        means=means,
                        filter_mag=args.filter,
                        verbose=args.verbose,
                    )
                    hist = all_hist["hist_out_params"]
                    hist_in_params = all_hist["hist_in_params"]
                    edges = all_hist["edges_out_params"]
                    edges_in_params = all_hist["edges_in_params"]
                else:
                    hist_conditional = get_histogram_cond(
                        args.file,
                        cond=args.conditional_hist,
                        additional_params=args.additional_hist_params,
                        only_asc600=args.only_asc600,
                        inoutflow_time=args.inoutflow_time,
                        means=means,
                        filter_mag=args.filter,
                        verbose=args.verbose,
                    )
            if args.save_histogram != "no":
                print("########### Store histograms ###########")
                file_path = args.save_histogram
                if len(args.conditional_hist) == 0:
                    with open(file_path + "_all_hist.pkl", "wb") as f:
                        pickle.dump(all_hist, f)
                else:
                    with open(file_path + "_hist_conditional.pkl", "wb") as f:
                        pickle.dump(hist_conditional, f)

            if args.plot_type == "all" or args.plot_type == "hist_out":
                print(
                    "########### Plot histograms for model state variables ###########"
                )
                traj_plot_histogram_out(
                    out_params=list(edges.keys()),
                    filename=args.out_file,
                    edges=edges,
                    hist=hist,
                    width=args.width,
                    height=args.height,
                    title=None,
                    verbose=args.verbose,
                )
            if args.plot_type == "all" or args.plot_type == "hist_in":
                print("########### Plot histograms for model parameters ###########")
                if len(args.in_params) == 0:
                    in_params_plot = list(edges_in_params[list(edges.keys())[0]].keys())
                else:
                    in_params_plot = args.in_params
                traj_plot_histogram_inp(
                    in_params=in_params_plot,
                    filename=args.out_file,
                    edges_in_params=edges_in_params,
                    hist_in_params=hist_in_params,
                    width=args.width,
                    height=args.height,
                    title=None,
                    verbose=args.verbose,
                )
            if args.plot_type == "all" or args.plot_type == "heat":
                print("########### Plot heatmaps ###########")
                if len(args.in_params) == 0:
                    in_params_plot = list(edges_in_params[list(edges.keys())[0]].keys())
                else:
                    in_params_plot = args.in_params
                plot_heatmap_traj(
                    in_params=in_params_plot,
                    filename=args.out_file,
                    edges_in_params=edges_in_params,
                    hist_in_params=hist_in_params,
                    width=args.width,
                    height=args.height,
                    title=None,
                    verbose=args.verbose,
                )

        print("########### Some statistics ###########")
        if args.load_statistics != "no":
            print("########### Load sums ###########")
            file_path = args.load_statistics
            with open(file_path + "_sums.pkl", "rb") as f:
                sums = pickle.load(f)
            with open(file_path + "_sums_phase.pkl", "rb") as f:
                sums_phase = pickle.load(f)
        else:
            sums = get_sums(
                args.file,
                args.save_statistics,
                only_asc600=args.only_asc600,
                inoutflow_time=args.inoutflow_time,
            )
            sums_phase = get_sums_phase(
                args.file,
                args.save_statistics,
                only_asc600=args.only_asc600,
                inoutflow_time=args.inoutflow_time,
            )

        if args.load_covariance != "no":
            print("########### Load covariance matrix ###########")
            file_path = args.load_covariance
            if args.plot_type == "cov_heat" or args.plot_type == "all":
                with open(file_path + "_means.pkl", "rb") as f:
                    means = pickle.load(f)
                with open(file_path + "_covariance_matrix.pkl", "rb") as f:
                    cov = pickle.load(f)
            if args.plot_type == "cov_heat_phases" or args.plot_type == "all":
                with open(file_path + "_means_phases.pkl", "rb") as f:
                    means_phase = pickle.load(f)
                with open(file_path + "_covariance_matrix_phases.pkl", "rb") as f:
                    cov_phase = pickle.load(f)

        top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic = traj_get_top_params(
            sums, sums.keys(), 10, 1
        )
        (
            top_magn_set_phase,
            top10_set_phase,
            top_magn_sens_dic_phase,
            top_sens_dic_phase,
        ) = traj_get_top_params(sums_phase, sums_phase.keys(), 10, 1)
        text = print_top_parameters(
            top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic
        )
        text += "\n########### With phases ###########\n"
        text += print_top_parameters(
            top_magn_set_phase,
            top10_set_phase,
            top_magn_sens_dic_phase,
            top_sens_dic_phase,
        )
        if args.out_file != "none":
            filename = args.out_file
            ending = filename.split(".")[-1]
            filename = filename[0 : -len(ending) - 1] + ".txt"
            with open(filename, "w+") as f:
                f.write(text)

        if (
            args.plot_type == "all" or args.plot_type == "cov_heat"
        ) and args.load_covariance == "no":
            print("########### Calculate covariance matrix ###########")
            means, cov = get_cov_matrix(
                args.file,
                in_params=list(top_magn_set),
                filepath=args.save_statistics,
                only_asc600=args.only_asc600,
                inoutflow_time=args.inoutflow_time,
            )
        if (
            args.plot_type == "all" or args.plot_type == "cov_heat_phases"
        ) and args.load_covariance == "no":
            print("########### Calculate covariance matrix with phases ###########")
            means_phase, cov_phase = get_cov_matrix_phase(
                args.file,
                in_params=list(top_magn_set_phase),
                filepath=args.save_statistics,
                only_asc600=args.only_asc600,
                inoutflow_time=args.inoutflow_time,
            )

        if len(args.conditional_hist) != 0 and (
            args.plot_type == "all" or args.plot_type == "2D_hist"
        ):
            print("########### Plot 2D histogram ###########")
            if args.conditional_hist[0] == "all":
                conds = list(hist_conditional.keys())
                conditions = []
                for c in conds:
                    if c != "edges_in_params" and c != "edges_out_params":
                        conditions.append(c)
            else:
                conditions = args.conditional_hist
            wrt_params = list(hist_conditional["edges_in_params"].keys())
            plot_heatmap_histogram(
                hist_conditional,
                in_params=list(
                    hist_conditional["edges_in_params"][wrt_params[0]].keys()
                ),
                out_params=list(hist_conditional["edges_out_params"].keys()),
                conditions=conditions,
                title=None,
                filename=args.out_file,
                verbose=args.verbose,
            )

        store_type = args.out_file.split(".")[-1]
        store_path = args.out_file[0 : -len(store_type) - 1]
        if args.plot_type == "all" or args.plot_type == "cov_heat":
            print("########### Plot covariance matrix ###########")
            for out_p in tqdm(cov):
                _ = plot_heatmap_cov(
                    data_in=cov,
                    out_param=out_p,
                    names=list(means[out_p].keys()),
                    title=f"Heatmap of Covariance (Gradients for {out_p})",
                    plot_type="all",
                    norm=mpl_col.LogNorm(),
                    filename=f"{store_path}_{out_p}_all_log.{store_type}",
                    width=args.width,
                    height=args.height,
                )
                plt.clf()
                _ = plot_heatmap_cov(
                    data_in=cov,
                    out_param=out_p,
                    names=list(means[out_p].keys()),
                    title=f"Heatmap of Covariance (Gradients for {out_p})",
                    plot_type="all",
                    filename=f"{store_path}_{out_p}_all.{store_type}",
                    width=args.width,
                    height=args.height,
                )
                plt.clf()
                _ = plot_heatmap_cov(
                    data_in=cov,
                    out_param=out_p,
                    names=list(means[out_p].keys()),
                    in_params=top_magn_sens_dic[out_p],
                    title=f"Heatmap of Covariance (Gradients for {out_p})",
                    plot_type="all",
                    norm=mpl_col.LogNorm(),
                    filename=f"{store_path}_{out_p}_some_log.{store_type}",
                    width=args.width,
                    height=args.height,
                )
                plt.clf()
                _ = plot_heatmap_cov(
                    data_in=cov,
                    out_param=out_p,
                    names=list(means[out_p].keys()),
                    in_params=top_magn_sens_dic[out_p],
                    title=f"Heatmap of Covariance (Gradients for {out_p})",
                    plot_type="all",
                    filename=f"{store_path}_{out_p}_some.{store_type}",
                    width=args.width,
                    height=args.height,
                )
                plt.clf()
        if args.plot_type == "all" or args.plot_type == "cov_heat_phases":
            print("########### Plot covariance matrix with phases ###########")
            for phase in tqdm(cov_phase):
                if "neutral" in phase:
                    continue
                phase_strip = phase.strip()
                for out_p in tqdm(cov_phase[phase], leave=False):
                    _ = plot_heatmap_cov(
                        data_in=cov_phase[phase],
                        out_param=out_p,
                        names=list(means_phase[phase][out_p].keys()),
                        title=f"({phase_strip}) Heatmap of Covariance (Gradients for {out_p})",
                        plot_type="all",
                        norm=mpl_col.LogNorm(),
                        filename=f"{store_path}_{phase_strip}_{out_p}_all_log.{store_type}",
                        width=args.width,
                        height=args.height,
                    )
                    plt.clf()
                    _ = plot_heatmap_cov(
                        data_in=cov_phase[phase],
                        out_param=out_p,
                        names=list(means_phase[phase][out_p].keys()),
                        title=f"({phase_strip}) Heatmap of Covariance (Gradients for {out_p})",
                        plot_type="all",
                        filename=f"{store_path}_{phase_strip}_{out_p}_all.{store_type}",
                        width=args.width,
                        height=args.height,
                    )
                    plt.clf()
                    _ = plot_heatmap_cov(
                        data_in=cov_phase[phase],
                        out_param=out_p,
                        names=means_phase[phase][out_p].keys(),
                        in_params=top_magn_sens_dic_phase[phase + " " + out_p],
                        title=f"({phase_strip}) Heatmap of Covariance (Gradients for {out_p})",
                        plot_type="all",
                        norm=mpl_col.LogNorm(),
                        filename=f"{store_path}_{phase_strip}_{out_p}_some_log.{store_type}",
                        width=args.width,
                        height=args.height,
                    )
                    plt.clf()
                    _ = plot_heatmap_cov(
                        data_in=cov_phase[phase],
                        out_param=out_p,
                        names=means_phase[phase][out_p].keys(),
                        in_params=top_magn_sens_dic_phase[phase + " " + out_p],
                        title=f"({phase_strip}) Heatmap of Covariance (Gradients for {out_p})",
                        plot_type="all",
                        filename=f"{store_path}_{phase_strip}_{out_p}_some.{store_type}",
                        width=args.width,
                        height=args.height,
                    )
                    plt.clf()
    else:
        ds = xr.open_dataset(args.file, decode_times=False, engine="netcdf4")
        (
            out_params,
            top20_list,
            top10_list,
            top20_sens_dic,
            top10_sens_dic,
        ) = get_top_list(ds, True, args.verbose)
        (
            top_one_order_list,
            top_two_orders_list,
            top_three_orders_list,
        ) = get_magnitude_list(ds, out_params, True, args.verbose)

        pd.set_option("display.max_rows", 100)
        pd.set_option("display.max_columns", 10)

        with pd.option_context(
            "display.max_rows",
            100,
            "display.max_columns",
            10,
            "display.expand_frame_repr",
            False,
        ):
            text = print_table_top_lists(
                [top10_list, top20_list],
                [top_one_order_list, top_two_orders_list, top_three_orders_list],
            )

        text += print_unique_params(top10_sens_dic)
        text += print_correlation_broad(ds, out_params)
        text += print_correlation_mean(ds, out_params)
        sort_key_list, table_dic, text_tmp = print_latex_tables(ds, 10, args.verbose)
        text += text_tmp
        text += print_variable_with_important_params(sort_key_list)
        text += print_param_types(ds, table_dic)
        text += print_large_impact_no_sens(ds)
        if args.out_file != "none":
            filename = args.out_file
            ending = filename.split(".")[-1]
            filename = filename[0 : -len(ending) - 1] + ".txt"
            with open(filename, "w+") as f:
                f.write(text)


if __name__ == "__main__":
    import argparse
    import textwrap

    from latexify import *

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Get statistics of a final, post-processed dataset with mean squared deviation and 
            predicted mean squared deviation.
            Or get statistics and plot histograms for files from a sensitivity analysis simulation along
            trajectories, e.g., by using
            python get_stats.py --file /project/meteo/w2w/Z2/Z2_data_gradients/ --out_file /path/to/pics/histogram.png 
            The name of the plots will be changed automatically to store multiple plots.
            Beware that creating the histogram may take a while. You can use 
            --save_histogram /path/to/folder/
            to store the histogram and edges to disk. 
            Some statistics are done after plotting which may take a while as well.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file",
        default="../data/vladiana_ensembles_postprocess/merged_independent.nc",
        help=textwrap.dedent(
            """\
            Path to post-processed file or to a folder with many files from a sensitivity analysis simulation.
            """
        ),
    )
    parser.add_argument(
        "--from_processed",
        action="store_true",
        help=textwrap.dedent(
            """\
            If true, then --file points to a post-processed file, otherwise --file points either to 
            a folder with statistics or with results from a sensitivity analysis.
            """
        ),
    )
    parser.add_argument(
        "--out_file",
        default="../pics/histogram.png",
        help=textwrap.dedent(
            """\
            Path and name to store histogram plots if the input is a set of trajectories with a sensitivity analysis
            simulation. Exchanges the ending with 'txt' and stores the results of any statistics in there unless it is 'none'.
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
        "--only_asc600",
        action="store_true",
        help=textwrap.dedent(
            """\
            Consider only time steps during the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--inoutflow_time",
        default=-1,
        type=int,
        help=textwrap.dedent(
            """\
            Consider only time steps during the fastest ascent and within the given range before (inflow) and after (outflow) of the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--plot_type",
        default="all",
        help=textwrap.dedent(
            """\
            Choose which plots to create. Options are
            all: All plots.
            hist_out: Histogram for output parameters.
            hist_in: Histogram for all input parameters.
            heat: Heatmap for all parameters.
            cov_heat: Heatmap of covariance matrix.
            cov_heat_phases: Heatmap of covariance matrix for each phase.
            2D_hist : 2D histogram (heatmap) where the x-axis is determined by conditional_hist.
            none: No plots.
            """
        ),
    )
    parser.add_argument(
        "--in_params",
        default=[],
        nargs="+",
        type=str,
        help=textwrap.dedent(
            """\
            During plotting of 'hist_in' and 'heat': plot only the given model parameters. If none are given, plot using all model parameters available.
            """
        ),
    )
    parser.add_argument(
        "--conditional_hist",
        type=str,
        nargs="+",
        default=[],
        help=textwrap.dedent(
            """\
            Calculate 2D histograms, where the "x-axis" is given by "conditional_hist".
            If "load_histogram" is set, you may use "all" to use all available parameters for the x-axis given by the dataset.
            """
        ),
    )
    parser.add_argument(
        "--additional_hist_params",
        type=str,
        nargs="+",
        default=[],
        help=textwrap.dedent(
            """\
            Additional parameters to create a histogram for, such as pressure or asc600. 
            """
        ),
    )
    parser.add_argument(
        "--load_histogram",
        default="no",
        help=textwrap.dedent(
            """\
            Load the histogram and edges with pickle from this path.
            """
        ),
    )
    parser.add_argument(
        "--save_histogram",
        default="no",
        help=textwrap.dedent(
            """\
            Store the histogram and edges with pickle to this path.
            """
        ),
    )
    parser.add_argument(
        "--save_statistics",
        default="no",
        help=textwrap.dedent(
            """\
            Store the covariance matrix and means to this path. This includes the matrix with and without
            phases. Also store the sums for all gradients to this path.
            """
        ),
    )
    parser.add_argument(
        "--load_statistics",
        default="no",
        help=textwrap.dedent(
            """\
            Load the sums for all gradients from this path.
            """
        ),
    )
    parser.add_argument(
        "--load_covariance",
        default="no",
        help=textwrap.dedent(
            """\
            Load the covariance matrix and means from this path. This includes the matrix with and without
            phases. 
            """
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=textwrap.dedent(
            """\
            More output, i.e., in each intermediate step for building tables.
            """
        ),
    )

    args = parser.parse_args()
    main(args)
