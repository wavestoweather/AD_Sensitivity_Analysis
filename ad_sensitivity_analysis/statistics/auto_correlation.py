"""Get auto-correlation for variables and trajectories.

"""
import os

import numpy as np

from tqdm.auto import tqdm
import xarray as xr

from ad_sensitivity_analysis.data_handler.loader import load_dataset_part
from ad_sensitivity_analysis.plot.latexify import param_id_map


def _get_corr(ds_tmp, local_delay, n, data, d_i):
    """

    Parameters
    ----------
    ds_tmp
    local_delay
    n
    data
    d_i

    Returns
    -------

    """
    ds_tmp1 = ds_tmp.isel({"time": np.arange(n - local_delay)})
    ds_tmp2 = ds_tmp.isel({"time": np.arange(local_delay, n)})
    no_delay = ds_tmp1 - ds_tmp1.mean(dim="time")
    delayed = ds_tmp2 - ds_tmp2.mean(dim="time")
    corr = np.asarray(
        (no_delay * delayed).sum(dim="time")
        / ((n - local_delay) * ds_tmp1.std(dim="time") * ds_tmp1.std(dim="time"))
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


def _get_corr_ds_out_param(
    data_outp,
    ds,
    col,
    n,
    phase_i=None,
    f_i=None,
    delay=10,
):
    if isinstance(delay, int):
        if phase_i is None and f_i is None:
            _get_corr(ds[col], delay, n, data_outp[col][1], 0)
        elif phase_i is None:
            _get_corr(ds[col], delay, n, data_outp[col][1][f_i], 0)
        elif f_i is None:
            _get_corr(
                ds[col],
                delay,
                n,
                data_outp[col][1][phase_i],
                0,
            )
        else:
            _get_corr(
                ds[col],
                delay,
                n,
                data_outp[col][1][f_i, phase_i],
                0,
            )
    else:
        for delay_idx, delay_value in enumerate(delay):
            if phase_i is None and f_i is None:
                _get_corr(ds[col], delay_value, n, data_outp[col][1], delay_idx)
            elif phase_i is None:
                _get_corr(
                    ds[col],
                    delay_value,
                    n,
                    data_outp[col][1][f_i],
                    delay_idx,
                )
            elif f_i is None:
                _get_corr(
                    ds[col],
                    delay_value,
                    n,
                    data_outp[col][1][phase_i],
                    delay_idx,
                )
            else:
                _get_corr(
                    ds[col],
                    delay_value,
                    n,
                    data_outp[col][1][f_i, phase_i],
                    delay_idx,
                )


def _get_corr_ds_in_param(
    data_inp,
    ds,
    out_params,
    col,
    n,
    out_param_coord,
    phase_i=None,
    f_i=None,
    delay=10,
):
    out_p_i = 0
    for out_p in out_params:
        if isinstance(delay, int):
            if phase_i is None and f_i is None:
                _get_corr(
                    ds[col].sel({out_param_coord: out_p}),
                    delay,
                    n,
                    data_inp[col][1][out_p_i],
                    0,
                )
            elif phase_i is None:
                _get_corr(
                    ds[col].sel({out_param_coord: out_p}),
                    delay,
                    n,
                    data_inp[col][1][f_i, out_p_i],
                    0,
                )
            elif f_i is None:
                _get_corr(
                    ds[col].sel({out_param_coord: out_p}),
                    delay,
                    n,
                    data_inp[col][1][phase_i, out_p_i],
                    0,
                )
            else:
                _get_corr(
                    ds[col].sel({out_param_coord: out_p}),
                    delay,
                    n,
                    data_inp[col][1][f_i, phase_i, out_p_i],
                    0,
                )
        else:
            for delay_idx, delay_val in enumerate(delay):
                if phase_i is None and f_i is None:
                    _get_corr(
                        ds[col].sel({out_param_coord: out_p}),
                        delay_val,
                        n,
                        data_inp[col][1][out_p_i],
                        delay_idx,
                    )
                elif phase_i is None:
                    _get_corr(
                        ds[col].sel({out_param_coord: out_p}),
                        delay_val,
                        n,
                        data_inp[col][1][f_i, out_p_i],
                        delay_idx,
                    )
                elif f_i is None:
                    _get_corr(
                        ds[col].sel({out_param_coord: out_p}),
                        delay_val,
                        n,
                        data_inp[col][1][phase_i, out_p_i],
                        delay_idx,
                    )
                else:
                    _get_corr(
                        ds[col].sel({out_param_coord: out_p}),
                        delay_val,
                        n,
                        data_inp[col][1][f_i, phase_i, out_p_i],
                        delay_idx,
                    )
        out_p_i += 1


# pylint: disable=too-many-arguments, too-many-locals
def _get_corr_ds(
    data_outp,
    data_inp,
    ds,
    params,
    out_params,
    out_param_coord,
    phase_i=None,
    f_i=None,
    verbose=False,
    leave=False,
    delay=10,
):
    """

    Parameters
    ----------
    data_outp
    data_inp
    ds
    params
    out_params
    out_param_coord
    phase_i
    f_i
    verbose
    leave
    delay

    Returns
    -------

    """
    n = len(ds["time"])
    for col in tqdm(params, leave=leave) if verbose else params:
        if col[0] != "d" or col == "deposition":
            _get_corr_ds_out_param(
                data_outp=data_outp,
                ds=ds,
                col=col,
                n=n,
                phase_i=phase_i,
                f_i=f_i,
                delay=delay,
            )
        else:
            _get_corr_ds_in_param(
                data_inp=data_inp,
                ds=ds,
                out_params=out_params,
                col=col,
                n=n,
                out_param_coord=out_param_coord,
                phase_i=phase_i,
                f_i=f_i,
                delay=delay,
            )


def _get_coords(
    ds=None,
    file_path=None,
    delay=10,
    phases=False,
    columns=None,
    inoutflow_time=-1,
):
    """

    Parameters
    ----------
    ds
    file_path
    delay
    phases
    columns
    inoutflow_time

    Returns
    -------

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
            n_ensembles = (
                len(ds["ensemble"])
                if len(ds["ensemble"]) > n_ensembles
                else n_ensembles
            )
            n_trajectories = (
                len(ds["trajectory"])
                if len(ds["trajectory"]) > n_trajectories
                else n_trajectories
            )
        ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    else:
        files = None
        n_ensembles = len(ds["ensemble"])
        n_trajectories = len(ds["trajectory"])
    if "Output_Parameter_ID" in ds:
        out_params = list(ds["Output_Parameter_ID"].values)
        param_names = [param_id_map[out_p] for out_p in out_params]
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
        columns = [
            col
            for col in ds
            if col not in ["phase", "step", "asc600", "time_after_ascent"]
        ]
    load_vars = columns
    if phases:
        auto_corr_coords["phase"] = phases_arr
        load_vars.append("phase")
    if inoutflow_time > 0:
        load_vars.append("asc600")

    if isinstance(delay, int):
        auto_corr_coords["delay"].append(delay)
    else:
        auto_corr_coords["delay"].extend(delay)
    return auto_corr_coords, out_param_coord, files, load_vars, out_params


def _auto_corr_single_file_phases(
    ds, auto_corr_coords, columns, out_params, out_param_coord, delay, verbose
):
    """

    Parameters
    ----------
    ds
    auto_corr_coords
    columns
    out_params
    out_param_coord
    delay
    verbose

    Returns
    -------

    """
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
        _get_corr_ds(
            data_outp=auto_corrs_outp,
            data_inp=auto_corrs_inp,
            ds=ds,
            params=columns,
            out_params=out_params,
            out_param_coord=out_param_coord,
            phase_i=phase_val,
            leave=True,
            verbose=verbose,
            delay=delay,
        )
    return auto_corrs_outp, auto_corrs_inp


def _auto_corr_single_file(
    ds, auto_corr_coords, columns, out_params, out_param_coord, delay, verbose
):
    """

    Parameters
    ----------
    ds
    auto_corr_coords
    columns
    out_params
    out_param_coord
    delay
    verbose

    Returns
    -------

    """
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
    _get_corr_ds(
        data_outp=auto_corrs_outp,
        data_inp=auto_corrs_inp,
        ds=ds,
        params=columns,
        out_params=out_params,
        out_param_coord=out_param_coord,
        leave=True,
        verbose=verbose,
        delay=delay,
    )


def _auto_corr_files_phases(
    files,
    file_path,
    only_asc600,
    inoutflow_time,
    load_vars,
    auto_corr_coords,
    columns,
    out_params,
    out_param_coord,
    delay,
    verbose,
):
    """

    Parameters
    ----------
    files
    file_path
    only_asc600
    inoutflow_time
    load_vars
    auto_corr_coords
    columns
    out_params
    out_param_coord
    delay
    verbose

    Returns
    -------

    """
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
            ds = load_dataset_part(
                f=file_path + f,
                only_asc600=only_asc600,
                inoutflow_time=inoutflow_time,
                load_params=load_vars,
            )
            if ds["phase"].dtype != str:
                phase_val = phase_i
            else:
                phase_val = phase
            _get_corr_ds(
                data_outp=auto_corrs_outp,
                data_inp=auto_corrs_inp,
                ds=ds.where(ds["phase"] == phase_val),
                params=columns,
                out_params=out_params,
                out_param_coord=out_param_coord,
                phase_i=phase_i,
                f_i=f_i,
                verbose=verbose,
                leave=False,
                delay=delay,
            )
    return auto_corrs_outp, auto_corrs_inp


def _auto_corr_files(
    files,
    file_path,
    only_asc600,
    inoutflow_time,
    load_vars,
    auto_corr_coords,
    columns,
    out_params,
    out_param_coord,
    delay,
    verbose,
):
    """

    Parameters
    ----------
    files
    file_path
    only_asc600
    inoutflow_time
    load_vars
    auto_corr_coords
    columns
    out_params
    out_param_coord
    delay
    verbose

    Returns
    -------

    """
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
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            load_params=load_vars,
        )
        _get_corr_ds(
            data_outp=auto_corrs_outp,
            data_inp=auto_corrs_inp,
            ds=ds,
            params=columns,
            out_params=out_params,
            out_param_coord=out_param_coord,
            phase_i=None,
            f_i=f_i,
            verbose=verbose,
            leave=False,
            delay=delay,
        )
    return auto_corrs_outp, auto_corrs_inp


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
    r"""
    Estimate auto-correlation via

    {\displaystyle {\hat {R}}(k)={\frac {1}{(n-k)\sigma_{t} \cdot \sigma_{t+k}}}
    \sum _{t=1}^{n-k}(X_{t}-\mu_{t} )(X_{t+k}-\mu_{t+k} )}

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
    auto_corr_coords, out_param_coord, files, load_vars, out_params = _get_coords(
        ds=ds,
        file_path=file_path,
        delay=delay,
        phases=phases,
        columns=columns,
        inoutflow_time=inoutflow_time,
    )
    if files is None:
        if phases:
            auto_corrs_outp, auto_corrs_inp = _auto_corr_single_file_phases(
                ds=ds,
                auto_corr_coords=auto_corr_coords,
                columns=columns,
                out_params=out_params,
                out_param_coord=out_param_coord,
                delay=delay,
                verbose=verbose,
            )
        else:
            _auto_corr_single_file(
                ds=ds,
                auto_corr_coords=auto_corr_coords,
                columns=columns,
                out_params=out_params,
                out_param_coord=out_param_coord,
                delay=delay,
                verbose=verbose,
            )
    else:
        auto_corr_coords["file"] = files
        if phases:
            auto_corrs_outp, auto_corrs_inp = _auto_corr_files_phases(
                files=files,
                file_path=file_path,
                only_asc600=only_asc600,
                inoutflow_time=inoutflow_time,
                load_vars=load_vars,
                auto_corr_coords=auto_corr_coords,
                columns=columns,
                out_params=out_params,
                out_param_coord=out_param_coord,
                delay=delay,
                verbose=verbose,
            )
        else:
            auto_corrs_outp, auto_corrs_inp = _auto_corr_files(
                files=files,
                file_path=file_path,
                only_asc600=only_asc600,
                inoutflow_time=inoutflow_time,
                load_vars=load_vars,
                auto_corr_coords=auto_corr_coords,
                columns=columns,
                out_params=out_params,
                out_param_coord=out_param_coord,
                delay=delay,
                verbose=verbose,
            )
    data_vars = auto_corrs_outp
    for key, value in auto_corrs_inp.items():
        data_vars[key] = value
    return xr.Dataset(data_vars=data_vars, coords=auto_corr_coords)
