"""Calculate statistics over all trajectories.

"""
import os
import pickle

import numpy as np
from tqdm.auto import tqdm
import xarray as xr

from ad_sensitivity_analysis.data_handler import filter_data
from ad_sensitivity_analysis.data_handler import loader
from ad_sensitivity_analysis.plot import latexify


def get_coords(file_path, files, model_params=None):
    """

    Parameters
    ----------
    file_path
    files
    model_params

    Returns
    -------

    """
    ds = loader.load_dataset_part(file_path + files[0])
    if model_params is None:
        model_params = []
        for col in ds:
            if col[0] == "d" and col != "deposition":
                model_params.append(col)
    if "Output_Parameter_ID" in ds:
        out_coord = "Output_Parameter_ID"
        out_params = []
        for out_p in ds["Output_Parameter_ID"]:
            out_params.append(latexify.param_id_map[out_p.item()])
        out_params_coord_id = ds["Output_Parameter_ID"].values
    else:
        out_coord = "Output Parameter"
        out_params = ds["Output Parameter"].values
        out_params_coord_id = ds["Output Parameter"].values
    max_traj = np.max(ds["trajectory"]).item()
    for f in files:
        ds = loader.load_dataset_part(file_path + f, load_params="trajectory")
        if np.max(ds["trajectory"]).item() > max_traj:
            max_traj = np.max(ds["trajectory"]).item()
    max_traj += 1  # Trajectory enumeration starts with zero...
    return (
        {
            "Output Parameter": out_params,
            "file": files,
            "trajectory": np.arange(max_traj),
            "phase": ["warm phase", "mixed phase", "ice phase", "neutral phase", "any"],
            "flow": ["inflow", "ascent", "outflow", "any"],
        },
        ds["phase"].dtype,
        out_coord,
        out_params_coord_id,
    )


def get_phase(ds_tmp, phase, phase_type):
    """

    Parameters
    ----------
    ds_tmp
    phase
    phase_type

    Returns
    -------

    """
    if phase == "any":
        return ds_tmp
    if phase_type in (str, object):
        phase_idx = phase
    else:
        phase_idx = np.argwhere(
            np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
            == phase
        )[0].item()
    return ds_tmp.where(ds_tmp["phase"] == phase_idx)


def get_ranking(ds, model_params):
    """

    Parameters
    ----------
    ds
    model_params

    Returns
    -------

    """
    ds_avg = ds[model_params].mean(dim="time", skipna=True)
    param_avg = [[] for _ in ds_avg["trajectory"].values]
    for traj_idx in ds_avg["trajectory"].values:
        ds_avg2 = ds_avg.sel({"trajectory": traj_idx})
        for param in ds_avg:
            param_avg[traj_idx].append((param, ds_avg2[param].values.item()))
    return param_avg


# pylint: disable=too-many-locals
def set_ranking_from_dataset(
    ds,
    model_params,
    data_vars,
    f_i,
    coords,
    phase_type,
    inoutflow_time,
    out_params_coord_id,
    out_coord,
    worst_rank,
):
    """

    Parameters
    ----------
    ds
    model_params
    data_vars
    f_i
    coords
    phase_type
    inoutflow_time
    out_params_coord_id
    out_coord
    worst_rank

    Returns
    -------

    """
    ds = ds.isel({"ensemble": 0})
    ds[model_params] = np.abs(ds[model_params])
    ds_lons_lats = ds[["lon", "lat", "asc600"]].where(ds["asc600"])[["lon", "lat"]]
    for traj_idx in range(len(ds_lons_lats["lon"].values)):
        data_vars["lon"][1][f_i, traj_idx] = ds_lons_lats["lon"][traj_idx][
            ~np.isnan(ds_lons_lats["lon"][traj_idx])
        ][0]
        data_vars["lat"][1][f_i, traj_idx] = ds_lons_lats["lat"][traj_idx][
            ~np.isnan(ds_lons_lats["lat"][traj_idx])
        ][0]

    for phase_i, phase in enumerate(coords["phase"]):
        ds_phase = get_phase(ds, phase, phase_type)
        avg_ascent_tmp = ds_phase["asc600"].count(dim="time").values
        asc600_steps_tmp = ds_phase["w"].mean(dim="time")
        for traj_idx, traj_ascent_val in enumerate(avg_ascent_tmp):
            data_vars["avg ascent"][1][f_i, traj_idx, phase_i] = traj_ascent_val
            data_vars["asc600 step count"][1][
                f_i, traj_idx, phase_i
            ] = asc600_steps_tmp[traj_idx]
        for flow_i, flow in enumerate(coords["flow"]):
            if inoutflow_time < 0 and flow != "any" and flow != "ascent":
                continue
            ds_flow = filter_data.filter_by_flow(ds_phase, flow, inoutflow_time)
            for out_p_i, out_p in enumerate(out_params_coord_id):
                ds_out = ds_flow.sel({out_coord: out_p})
                param_values = get_ranking(ds_out, model_params)
                for traj_idx, traj_param_val in enumerate(param_values):
                    sort_params_set_rank(
                        data_vars=data_vars,
                        param_values=param_values,
                        traj_idx=traj_idx,
                        worst_rank=worst_rank,
                        out_p_i=out_p_i,
                        f_i=f_i,
                        phase_i=phase_i,
                        flow_i=flow_i,
                        traj_param_val=traj_param_val,
                    )


def sort_params_set_rank(
    data_vars,
    param_values,
    traj_idx,
    worst_rank,
    out_p_i,
    f_i,
    phase_i,
    flow_i,
    traj_param_val,
):
    """

    Parameters
    ----------
    data_vars
    param_values
    traj_idx
    worst_rank
    out_p_i
    f_i
    phase_i
    flow_i
    traj_param_val

    Returns
    -------

    """
    param_values[traj_idx].sort(key=lambda x: x[1])
    param_values[traj_idx] = traj_param_val[::-1]
    last_val = 0
    rank_offset = 1  # zero is the "invalid number" in our case
    for rank, param_pair in enumerate(param_values[traj_idx]):
        if rank > 0:
            if last_val == param_pair[1]:
                rank_offset -= 1
        else:
            last_val = param_pair[1]
        rank += rank_offset
        if param_pair[1] == 0:
            data_vars[f"{param_pair[0]} rank"][1][
                out_p_i, f_i, traj_idx, phase_i, flow_i
            ] = worst_rank
        elif np.isnan(param_pair[1]):
            data_vars[f"{param_pair[0]} rank"][1][
                out_p_i, f_i, traj_idx, phase_i, flow_i
            ] = 0
        else:
            data_vars[f"{param_pair[0]} rank"][1][
                out_p_i, f_i, traj_idx, phase_i, flow_i
            ] = rank
        data_vars[f"{param_pair[0]} avg"][1][
            out_p_i, f_i, traj_idx, phase_i, flow_i
        ] = param_pair[1]


def create_rank_traj_dataset(file_path, inoutflow_time=240, model_params=None):
    """
    Ranked index of parameters for each trajectory.
    Count model parameter / process index occurrence over all trajectories.
    Additional information that might be useful:
    location of ascent. Average ascent. Number of time steps with ascent.

    Parameters
    ----------
    file_path : string
        Path to a set of NetCDF-files from a sensitivity analysis.
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.
    model_params : List of strings or None
        List of model parameters to create a rank for. If None is given, all model parameters are used.

    Returns
    -------
    xarray.Dataset:
    dims:   Output Parameter (QV, latent_heat, latent_cool)
            trajectory       (0, 1, 2, ...)
            file             (traj20161003.nc, ...)
            # index            (0, 1, 2, ...)
            phase            (warm phase, mixed phase, ice phase, neutral phase, any)
            flow             (inflow, ascent, outflow, any)

    vars:   Model Parameter Rank (Output Parameter, file, trajectory, phase, flow)
            Model Parameter Avg (Output Parameter, file, trajectory, phase, flow)
            lon             (file, trajectory)
            lat             (file, trajectory)
            avg ascent      (file, trajectory, phase)
            asc600 steps    (file, trajectory, phase)
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    coords, phase_type, out_coord, out_params_coord_id = get_coords(
        file_path=file_path,
        files=files,
        model_params=model_params,
    )
    worst_rank = len(model_params) + 1
    data_vars = {
        f"{param} rank": (
            ["Output Parameter", "file", "trajectory", "phase", "flow"],
            np.zeros(
                (
                    len(coords["Output Parameter"]),
                    len(coords["file"]),
                    len(coords["trajectory"]),
                    len(coords["phase"]),
                    len(coords["flow"]),
                )
            ),
        )
        for param in model_params
    }
    for param in model_params:
        data_vars[f"{param} avg"] = (
            ["Output Parameter", "file", "trajectory", "phase", "flow"],
            np.zeros(
                (
                    len(coords["Output Parameter"]),
                    len(coords["file"]),
                    len(coords["trajectory"]),
                    len(coords["phase"]),
                    len(coords["flow"]),
                ),
            ),
        )
    data_vars["lon"] = (
        ["file", "trajectory"],
        np.zeros((len(coords["file"]), len(coords["trajectory"]))),
    )
    data_vars["lat"] = (
        ["file", "trajectory"],
        np.zeros((len(coords["file"]), len(coords["trajectory"]))),
    )
    data_vars["avg ascent"] = (
        ["file", "trajectory", "phase"],
        np.zeros(
            (len(coords["file"]), len(coords["trajectory"]), len(coords["phase"]))
        ),
    )
    data_vars["asc600 step count"] = (
        ["file", "trajectory", "phase"],
        np.zeros(
            (len(coords["file"]), len(coords["trajectory"]), len(coords["phase"]))
        ),
    )

    for f_i, f in enumerate(tqdm(files)):
        ds = loader.load_dataset_part(
            file_path + f,
            inoutflow_time=inoutflow_time,
            load_params=model_params + ["lon", "lat", "w", "asc600", "phase"],
        )
        set_ranking_from_dataset(
            ds=ds,
            model_params=model_params,
            data_vars=data_vars,
            f_i=f_i,
            coords=coords,
            phase_type=phase_type,
            inoutflow_time=inoutflow_time,
            out_params_coord_id=out_params_coord_id,
            out_coord=out_coord,
            worst_rank=worst_rank,
        )

    ds = xr.Dataset(data_vars=data_vars, coords=coords)
    for param in model_params:
        ds[f"{param} rank"].attrs = {
            "no data": 0,
            "rank for zero gradients": worst_rank,
        }
    return ds


# pylint: disable=too-many-nested-blocks,too-many-branches,too-many-locals
def get_matrix(
    ds,
    store_path=None,
    corr=True,
    verbose=False,
):
    """

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with ranks per trajectory created by create_rank_traj_dataset().
    store_path : string
        If a store path is given then dumps the correlation matrix to f'{store_path}_correlation_matrix_per_traj.pkl'.
    corr : bool
        If true, calculate the Pearson correlation. Otherwise, calculate the covariance.
    verbose : bool
        Print progressbars.

    Returns
    -------
    Dictionary with output parameter, flow, and phase as keys and values are correlation/covariance matrices.
    Also a list of model parameter names for each column/row.
    """
    n = len(ds)

    corr_matrix = {
        out_p.item(): {
            flow.item(): {phase.item(): np.zeros((n, n)) for phase in ds["phase"]}
            for flow in ds["flow"]
        }
        for out_p in ds["Output Parameter"]
    }
    col_names = list(ds.keys())
    for out_p in tqdm(ds["Output Parameter"]) if verbose else ds["Output Parameter"]:
        ds_tmp = ds.sel({"Output Parameter": out_p})
        for flow in tqdm(ds["flow"], leave=False) if verbose else ds["flow"]:
            ds_tmp2 = ds_tmp.sel({"flow": flow})
            for phase in tqdm(ds["phase"], leave=False) if verbose else ds["phase"]:
                ds_tmp3 = ds_tmp2.sel({"phase": phase})
                for i, col in enumerate(
                    tqdm(col_names, leave=False) if verbose else col_names
                ):
                    if "rank" in col:
                        ds_tmp4 = ds_tmp3.where(ds_tmp3[col] > 0)
                    else:
                        ds_tmp4 = ds_tmp3
                    for j, col2 in enumerate(col_names):
                        if "rank" in col2:
                            ds_tmp5 = ds_tmp4.where(ds_tmp4[col2] > 0)
                        else:
                            ds_tmp5 = ds_tmp4
                        if corr:
                            val = xr.corr(ds_tmp5[col], ds_tmp5[col2]).item()
                        else:
                            val = xr.cov(ds_tmp5[col], ds_tmp5[col2]).item()
                        corr_matrix[out_p.item()][flow.item()][phase.item()][i, j] = val
    if store_path is not None:
        corr_and_names = {"correlation": corr_matrix, "column names": col_names}
        if corr:
            with open(store_path + "_correlation_matrix_per_traj.pkl", "wb") as f:
                pickle.dump(corr_and_names, f)
        else:
            with open(store_path + "_covariance_matrix_per_traj.pkl", "wb") as f:
                pickle.dump(corr_and_names, f)
    return corr_matrix, col_names
