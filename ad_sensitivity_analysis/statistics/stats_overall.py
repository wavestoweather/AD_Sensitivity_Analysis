"""Get statistics over all trajectories.

"""
import pickle
import os

import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import xarray as xr

from ad_sensitivity_analysis.data_handler.loader import load_dataset_part
from ad_sensitivity_analysis.plot.latexify import param_id_map


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
            param_name.append(param_id_map[out_p])
    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    sums = {}
    for f in tqdm(files):
        ds = load_dataset_part(
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

            if out_name in sums:
                sums[out_name] += df
            else:
                sums[out_name] = df
    if store_path is not None and store_path != "no":
        with open(store_path + "_sums.pkl", "wb") as f:
            pickle.dump(sums, f)
    return sums


# pylint: disable=too-many-locals
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
            param_name.append(param_id_map[idx])

    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    sums = {}
    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
        )
        ds[in_params] = np.abs(ds[in_params])
        if ds["phase"].dtype not in (str, np.uint64):
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({out_param_coord: out_p})
            for phase_i, phase in enumerate(
                ["warm phase", "mixed phase", "ice phase", "neutral phase"]
            ):
                if ds_tmp["phase"].dtype == str:
                    idx = np.where(ds_tmp["phase"] == phase)
                else:
                    idx = np.where(ds_tmp["phase"] == phase_i)
                dic = {in_p: None for in_p in in_params}
                for in_p in in_params:
                    dic[in_p] = [np.nansum(ds_tmp[in_p].values[idx])]
                if phase + " " + out_name in sums:
                    sums[phase + " " + out_name] += pd.DataFrame(dic)
                else:
                    sums[phase + " " + out_name] = pd.DataFrame(dic)
    if store_path is not None and store_path != "no":
        with open(store_path + "_sums_phase.pkl", "wb") as f:
            pickle.dump(sums, f)
    return sums


def _get_param_names(file_path, files, in_params):
    """

    Parameters
    ----------
    file_path
    files
    in_params
    only_asc600
    inoutflow_time

    Returns
    -------

    """
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
            param_name.append(param_id_map[idx])
    all_params = param_name + in_params
    return out_params, param_name, all_params, out_param_coord


def _set_means_counts(ds, out_param_coord, out_p, means, n_total, out_name):
    """

    Parameters
    ----------
    ds
    out_param_coord
    out_p
    means
    n_total
    out_name

    Returns
    -------

    """
    ds_tmp = ds.sel({out_param_coord: out_p})
    means_tmp = ds_tmp.mean(skipna=True)
    count_tmp = (~np.isnan(ds_tmp)).sum()
    for p in means[out_name]:
        if n_total[out_name][p] > 0:
            n_new = count_tmp[p].values.item() + n_total[out_name][p]
            if n_new > 0:
                means[out_name][p] = (
                    means[out_name][p] * n_total[out_name][p] / n_new
                    + means_tmp[p].values.item() * count_tmp[p].values.item() / n_new
                )
                n_total[out_name][p] = n_new
        else:
            n_total[out_name][p] = count_tmp[p].values.item()
            means[out_name][p] = means_tmp[p].values.item()


def _fill_cov_matrix(ds, cov, out_param_coord, out_p, means, n_total, out_name):
    """

    Parameters
    ----------
    ds
    cov
    out_param_coord
    out_p
    means
    n_total
    out_name

    Returns
    -------

    """
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
    out_params, param_name, all_params, out_param_coord = _get_param_names(
        file_path, files, in_params
    )

    means = {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
    cov = {
        out_p: np.zeros((len(all_params), len(all_params)), dtype=np.float64)
        for out_p in param_name
    }
    n_total = {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            load_params=all_params,
        )
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            _set_means_counts(ds, out_param_coord, out_p, means, n_total, out_name)

    n_total = {out_p: {p: 0.0 for p in all_params} for out_p in param_name}
    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            load_params=all_params,
        )
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            _fill_cov_matrix(ds, cov, out_param_coord, out_p, means, n_total, out_name)

    if store_path is not None and store_path != "no":
        with open(store_path + "_means.pkl", "wb") as f:
            pickle.dump(means, f)
        with open(store_path + "_covariance_matrix.pkl", "wb") as f:
            pickle.dump(cov, f)
    return means, cov


def _set_means_counts_phase(
    ds, out_param_coord, out_p, means, n_total, out_name, phases
):
    """

    Parameters
    ----------
    ds
    out_param_coord
    out_p
    means
    n_total
    out_name
    phases

    Returns
    -------

    """
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
                        means[phase][out_name][p] * n_total[phase][out_name][p] / n_new
                        + mean_tmp * count_tmp / n_new
                    )
                    n_total[phase][out_name][p] = n_new
            else:
                n_total[phase][out_name][p] = count_tmp
                means[phase][out_name][p] = mean_tmp


def _fill_cov_matrix_phase(
    ds, cov, out_param_coord, out_p, means, n_total, out_name, phases
):
    """

    Parameters
    ----------
    ds
    cov
    out_param_coord
    out_p
    means
    n_total
    out_name
    phases

    Returns
    -------

    """
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
                            * (ds_tmp[p2].values[idx] - means[phase][out_name][p2])
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
                            * (ds_tmp[p2].values[idx] - means[phase][out_name][p2])
                        )
                    n_total[phase][out_name][p] = n_new


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
    out_params, param_name, all_params, out_param_coord = _get_param_names(
        file_path, files, in_params
    )
    more_params = ["phase"]
    if only_asc600 or inoutflow_time > 0:
        more_params.append("asc600")
    phases = ["warm phase", "mixed phase", "ice phase", "neutral phase"]
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
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            load_params=all_params + more_params,
        )
        if ds["phase"].dtype not in (str, np.uint64):
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            _set_means_counts_phase(
                ds, out_param_coord, out_p, means, n_total, out_name, phases
            )

    n_total = {
        phase: {out_p: {p: 0.0 for p in all_params} for out_p in param_name}
        for phase in phases
    }
    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            inoutflow_time=inoutflow_time,
            load_params=all_params + more_params,
        )
        if ds["phase"].dtype not in (str, np.uint64):
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            _fill_cov_matrix_phase(
                ds, cov, out_param_coord, out_p, means, n_total, out_name, phases
            )

    if store_path is not None and store_path != "no":
        with open(store_path + "_means_phases.pkl", "wb") as f:
            pickle.dump(means, f)
        with open(store_path + "_covariance_matrix_phases.pkl", "wb") as f:
            pickle.dump(cov, f)
    return means, cov
