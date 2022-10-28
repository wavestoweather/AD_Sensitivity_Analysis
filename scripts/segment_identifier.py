from typing import Union

import holoviews as hv
import itertools
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import os.path
import pandas as pd

try:
    from tqdm.auto import tqdm
except:
    from progressbar import progressbar as tqdm
from timeit import default_timer as timer
import warnings
import xarray as xr

from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.multiclass import OneVsRestClassifier
from skmultilearn.model_selection import iterative_train_test_split
from itertools import product

try:
    from latexify import (
        in_params_numeric_value_dic,
        parse_word,
        in_params_dic,
        physical_params,
        in_params_grouping,
        param_id_map,
        in_params_notation_mapping,
    )
except:
    from scripts.latexify import (
        in_params_numeric_value_dic,
        parse_word,
        in_params_dic,
        physical_params,
        in_params_grouping,
        param_id_map,
        in_params_notation_mapping,
    )


def d_unnamed(df):
    """
    Remove unnamed column from dataframe.

    Parameters
    ----------
    df : pandas.Dataframe
        A dataframe with one or more unnamed columns

    Returns
    -------
    Dataframe without unnamed columns.
    """
    return df.loc[:, ~df.columns.str.contains("^Unnamed")]


def rolling_window(a, window):
    """
    Create a rolling window for a given array of arbitrary dimensions.

    Parameters
    ----------
    a : np.ndarray
        Numpy array to create rolling window for.
    window : int
        Size of the window.

    Returns
    -------
    View of np.ndarray with rolling window.
    """
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def load_sensitivity(
    path,
    out_params,
    param_names,
    in_params=["da_1"],
    par_dim_name="Output Parameter",
):
    """
    Load the sensitivities from the unperturbed trajectory.

    Parameters
    ----------
    path: String
        path and name of file to load. Should be the unperturbed
        trajectory with 'Output Parameter' as column.
    out_params: List of string
        List of output parameters to get sensitivities for.
    time_integrated : bool
        If true, calculate the time integrated sensitivities and model states.
        Otherwise get the difference at each time step.
    in_params: List of string
        List of input parameters for the sensitivies.
    par_dim_name: String
        Name of the dimension for the output parameter.
    Returns
    -------
    xarray.DataFrame
        Dataframe with model state parameters, sensitivities
        and time after ascent.
    """
    ds = xr.open_dataset(path, decode_times=False, engine="netcdf4")
    ds = ds.loc[{par_dim_name: out_params}][
        param_names + in_params + ["time_after_ascent"]
    ]
    for in_p in in_params:
        ds[in_p] = ds[in_p]  # * (in_params_numeric_value_dic[in_p] * 0.1)
    return ds


def load_unperturbed(path, out_params, par_dim_name="Output Parameter"):
    """
    Load the unperturbed trajectory output parameters.
    Should not be necessary since load_sensitivity() already
    stores those values.

    Parameters
    ----------
    path: String
        path and name of file to load. Should be the unperturbed
        trajectory with 'Output Parameter' as column.
    out_params: List of string
        List of output parameters.

    Returns
    -------
    xarray.DataFrame
        Dataframe with model state parameters and time after ascent.
    """
    ds = load_sensitivity(path, out_params, par_dim_name)
    return ds[out_params + ["time_after_ascent"]]


def load_dataset(
    path,
    out_params,
    in_params,
    traj_list,
    time_integrated=False,
    traj_offset=0,
    verbosity=0,
    par_dim_name="Output Parameter",
):
    """
    Load all trajectory data and store the MSE for each ensemble
    and parameter as a result.


    Parameters
    ----------
    path: String
        Path to files to load. Within this folder, iterate
        over trajy..trajx where y and x are values from traj_list.
        Loads path/trajx.nc_wcb and path/trajx/in_param.nc_wcb
        where in_param is an input parameter from in_params.
        If traj_list is a list of files, no subsequent folder with trajx
        is assumed and instead only a single trajectory as source for
        this ensemble is assumed. This is usually the case for quantile
        trajectories.
    out_params: List of string or int
        List of output parameters to get MSE for.
    in_params: List of string
        List of input parameters to get MSE for.
    traj_list: List of int
        List of trajectory numbers to load or list of files.
    time_integrated : bool
        If true, calculate the time integrated sensitivities and model states.
        Otherwise get the difference at each time step.
    traj_offset : int
        Offset for the new index of the trajectories.
    verbosity : int
        Set verbosity level.

    Returns
    -------
    xarray.DataFrame
        DataFrame MSE and predicted error for each output and input parameter and trajectory.
        Also has output parameter values for not perturbed trajectory
    """
    trajectories = []
    if time_integrated:
        dim_order = (
            "Output Parameter",
            "Input Parameter",
            "trajectory",
        )
    else:
        dim_order = (
            "Output Parameter",
            "Input Parameter",
            "trajectory",
            "time_after_ascent",
        )
    idx_offset = 0
    n_timesteps = 0
    time_after_ascent = None
    one_iteration = False
    param_names = []
    for i in range(len(out_params)):
        param_names.append(param_id_map[out_params[i]])

    for traj_idx, traj in enumerate(traj_list):
        if os.path.isfile(path + traj):
            trajectories.append(traj_idx + traj_offset)
            one_iteration = True
        elif os.path.isfile(
            path + "traj" + str(traj) + "/" + in_params[0][1::] + ".nc_wcb"
        ):
            trajectories.append(traj + traj_offset)
        else:
            continue

        if n_timesteps == 0 and not time_integrated:
            not_perturbed_path = path + "traj" + str(traj) + "_notPerturbed.nc_wcb"
            # I wasn't consistent with the way not perturbed parameters are stored
            if not os.path.isfile(not_perturbed_path):
                not_perturbed_path = path + "traj" + str(traj) + "/_notPerturbed.nc_wcb"
            if not os.path.isfile(not_perturbed_path):
                not_perturbed_path = path + "_notPerturbed.nc_wcb"
            if not os.path.isfile(not_perturbed_path):
                not_perturbed_path = path + in_params[0][1::] + ".nc_wcb"
            if not os.path.isfile(not_perturbed_path):
                not_perturbed_path = path + in_params[0][1::] + ".nc"
            if not os.path.isfile(not_perturbed_path):
                not_perturbed_path = (
                    path + in_params[0][1::] + "_" + path.split("/")[-2] + ".nc"
                )
            val_df = load_sensitivity(
                path=not_perturbed_path,
                out_params=out_params,
                param_names=param_names,
                in_params=in_params,
                par_dim_name=par_dim_name,
            )
            val_array = val_df.loc[{par_dim_name: out_params[0], "trajectory": 0}][
                param_names[0]
            ]
            n_timesteps = len(val_array["time"])
            time_after_ascent = np.asarray(
                val_df.loc[{par_dim_name: out_params[0], "trajectory": 0}][
                    "time_after_ascent"
                ]
            ).flatten()
        if one_iteration:
            break
    if time_integrated:
        dims = (len(out_params), len(in_params), len(trajectories))
    else:
        dims = (len(out_params), len(in_params), len(trajectories), n_timesteps)
    mse = np.zeros(dims)
    predicted = np.zeros(dims)
    if not time_integrated:
        not_perturbed = np.zeros((len(out_params), len(trajectories), n_timesteps))

    for traj_idx, traj in enumerate(traj_list):
        if not os.path.isfile(path + traj) and not os.path.isfile(
            path + "traj" + str(traj) + "/" + in_params[0][1::] + ".nc_wcb"
        ):
            idx_offset += 1
            continue

        not_perturbed_path = path + "traj" + str(traj) + "_notPerturbed.nc_wcb"
        # I wasn't consistent with the way not perturbed parameters are stored
        if not os.path.isfile(not_perturbed_path):
            not_perturbed_path = path + "traj" + str(traj) + "/_notPerturbed.nc_wcb"
        if not os.path.isfile(not_perturbed_path):
            not_perturbed_path = path + "_notPerturbed.nc_wcb"
        if not os.path.isfile(not_perturbed_path):
            not_perturbed_path = path + in_params[0][1::] + ".nc_wcb"
        if not os.path.isfile(not_perturbed_path):
            not_perturbed_path = path + in_params[0][1::] + ".nc"
        if not os.path.isfile(not_perturbed_path):
            not_perturbed_path = (
                path + in_params[0][1::] + "_" + path.split("/")[-2] + ".nc"
            )

        if verbosity > 1:
            print(f"Loading from {not_perturbed_path} for index {traj_idx}")
        val_df = load_sensitivity(
            not_perturbed_path, out_params, param_names, in_params, par_dim_name
        )
        val_only_df = val_df.loc[{par_dim_name: out_params[0], "trajectory": 0}][
            param_names
        ]
        for out_idx, out_p in enumerate(out_params):
            out_p_name = param_names[out_idx]
            if not time_integrated:
                val_array = val_only_df[out_p_name]
                not_perturbed[out_idx, traj_idx - idx_offset, :] = np.asarray(
                    val_array
                ).flatten()
            sens_df = val_df.loc[{par_dim_name: out_p}][in_params]
            if verbosity > 2:
                for in_idx in tqdm(range(len(in_params))):
                    in_p = in_params[in_idx]
                    load_path = path + "traj" + str(traj) + "/" + in_p[1::] + ".nc_wcb"
                    if not os.path.isfile(load_path):
                        load_path = path + in_p[1::] + ".nc_wcb"
                    if not os.path.isfile(load_path):
                        load_path = path + in_p[1::] + ".nc"
                    if not os.path.isfile(load_path):
                        load_path = path + in_p[1::] + "_" + path.split("/")[-2] + ".nc"
                    ds = xr.open_dataset(
                        load_path,
                        decode_times=False,
                        engine="netcdf4",
                    )
                    ds = ds.loc[{"trajectory": ds["trajectory"][1::]}]
                    tmp1 = np.asarray(ds[out_p_name].load())
                    tmp2 = np.asarray(val_only_df[out_p_name])
                    if time_integrated:
                        tmp1_sum = np.nanmean(np.sum(tmp1, axis=1))
                        tmp2_sum = np.nanmean(np.sum(tmp2, axis=1))
                        mse[out_idx, in_idx, traj_idx - idx_offset] = (
                            tmp1_sum - tmp2_sum
                        ) ** 2
                        predicted[out_idx, in_idx, traj_idx - idx_offset] = np.sum(
                            np.asarray(sens_df[in_p])
                        )
                    else:
                        mse[out_idx, in_idx, traj_idx - idx_offset, :] = np.nanmean(
                            (tmp1 - tmp2) ** 2, axis=1
                        ).flatten()
                        predicted[
                            out_idx, in_idx, traj_idx - idx_offset, :
                        ] = np.asarray(sens_df[in_p]).flatten()
            else:
                for in_idx, in_p in enumerate(in_params):
                    load_path = path + "traj" + str(traj) + "/" + in_p[1::] + ".nc_wcb"
                    if not os.path.isfile(load_path):
                        load_path = path + in_p[1::] + ".nc_wcb"
                    if not os.path.isfile(load_path):
                        load_path = path + in_p[1::] + ".nc"
                    if not os.path.isfile(load_path):
                        load_path = path + in_p[1::] + "_" + path.split("/")[-2] + ".nc"
                    ds = xr.open_dataset(
                        load_path,
                        decode_times=False,
                        engine="netcdf4",
                    )
                    ds = ds.loc[{"trajectory": ds["trajectory"][1::]}]
                    tmp1 = np.asarray(ds[out_p_name].load())
                    tmp2 = np.asarray(val_only_df[out_p_name])
                    if time_integrated:
                        tmp1_sum = np.nanmean(np.sum(tmp1, axis=1))
                        tmp2_sum = np.nanmean(np.sum(tmp2, axis=1))
                        mse[out_idx, in_idx, traj_idx - idx_offset] = (
                            tmp1_sum - tmp2_sum
                        ) ** 2
                        predicted[out_idx, in_idx, traj_idx - idx_offset] = np.sum(
                            np.asarray(sens_df[in_p])
                        )
                    else:
                        mse[out_idx, in_idx, traj_idx - idx_offset, :] = np.nanmean(
                            (tmp1 - tmp2) ** 2, axis=1
                        ).flatten()
                        predicted[
                            out_idx, in_idx, traj_idx - idx_offset, :
                        ] = np.asarray(sens_df[in_p]).flatten()
        if one_iteration:
            break
    if time_integrated:
        return xr.Dataset(
            data_vars={
                "Predicted Squared Error": (list(dim_order), predicted ** 2),
                "Predicted Error": (list(dim_order), predicted),
                "Mean Squared Error": (list(dim_order), mse),
            },
            coords={
                "Output Parameter": param_names,
                "Input Parameter": in_params,
                "trajectory": trajectories,
            },
        )
    else:
        return xr.Dataset(
            data_vars={
                "Predicted Squared Error": (list(dim_order), predicted ** 2),
                "Predicted Error": (list(dim_order), predicted),
                "Mean Squared Error": (list(dim_order), mse),
                "Not Perturbed Value": (
                    ["Output Parameter", "trajectory", "time_after_ascent"],
                    not_perturbed,
                ),
            },
            coords={
                "Output Parameter": param_names,
                "Input Parameter": in_params,
                "trajectory": trajectories,
                "time_after_ascent": time_after_ascent,
            },
        )


def parse_load(
    data_path,
    out_params,
    all_params_list,
    time_integrated=False,
    store_many_appended_data=None,
    load_on_the_fly=False,
    min_time=-1000,
    verbosity=0,
):
    """
    Parse the args.data_path and corresponding arguments.

    Parameters
    ----------
    data_path : string
        Path to folders with ensemble datasets or to single NetCDF file
        with all data concatenated along 'trajectory' axis.
        If a path to numpy arrays is given, it is assumed to be a training
        or test set.
    out_params : list of string
        List of output parameters.
    all_params_list : list of string
        List of all input params to get predicted errors for.
    time_integrated : bool
        If true, calculate the time integrated sensitivities and model states.
        Otherwise get the difference at each time step.
    store_many_appended_data : string
        Store the appended input data to this path as NetCDF file for each appended
        version. Used mainly for debugging.
    load_on_the_fly : boolean
        Load data and find the segments on the fly for predicting a dataset,
        or load precalculated training or test set.
    verbosity : int
        Set verbosity level.
        0: No output except for exceptions
        1: Print datasets
        2: Print loading statements

    Returns
    -------
    Either array of datapaths for loading on the fly or xarray.Dataset with
    already loaded data.
    """
    traj_offset = 0
    if ".nc" in data_path:
        # ie data2_327.nc
        if verbosity > 1:
            print(f"Loading {data_path}")
        data = xr.open_dataset(data_path, decode_times=False, engine="netcdf4")
        if not time_integrated:
            data = data.where(data["time_after_ascent"] >= min_time, drop=True)
    else:
        if verbosity > 1:
            print(f"Checking {data_path}")

        paths = list(os.listdir(data_path))
        if not load_on_the_fly:
            data = None
            if verbosity > 1:
                print(f"Loading from {paths}")
            for p in range(len(paths)):
                path = data_path + paths[p] + "/"
                traj_dirs = list(os.listdir(path))
                if "traj" in traj_dirs[0]:
                    n_trajs = len(traj_dirs)
                    traj_list = np.arange(len(traj_dirs))
                else:
                    n_trajs = 1
                    traj_list = traj_dirs
                if verbosity > 1:
                    print(f"Loading from {path} with {n_trajs} trajectories")

                tmp = load_dataset(
                    path=path,
                    out_params=out_params,
                    in_params=all_params_list,
                    traj_list=traj_list,
                    time_integrated=time_integrated,
                    traj_offset=traj_offset,
                    verbosity=verbosity,
                    par_dim_name="Output_Parameter_ID",
                )

                if tmp is None:
                    continue

                if data is None:
                    data = tmp
                else:
                    data = xr.concat([data, tmp], dim="trajectory", join="outer")
                if store_many_appended_data is not None:
                    comp = dict(zlib=True, complevel=9)
                    encoding = {var: comp for var in data.data_vars}
                    data.to_netcdf(
                        path=f"{store_many_appended_data}data_{traj_offset}.nc",
                        encoding=encoding,
                        compute=True,
                        engine="netcdf4",
                        format="NETCDF4",
                        mode="w",
                    )
                traj_offset += n_trajs
            if not time_integrated:
                data = data.where(data["time_after_ascent"] >= min_time, drop=True)
        else:  # numpy arrays with training data
            for i in range(len(paths)):
                paths[i] = data_path + paths[i]
            data = np.sort(paths)
    return data


def find_segments(df, error_threshold=0, cooldown=0):
    """
    Iterate over time steps to mark the start of a segment with large errors, where large
    is defined via error_threshold. The dataframe needs to have a column (or index)
    "Output Parameter" and should be ordered by time or time_after_ascent.

    Parameters
    ----------
    df : xarray.DataFrame
        Dataframe with MSE created by load_dataset() for every model state
        parameter and all trajectories.
    error_threshold : float
        Threshold for errors to identify true segment starts.
    cooldown : int
        Minimum number of timesteps between last time the error threshold has
        been met and the next time step. This is useful if the original data
        is based on ensembles that start at every few time step which leads to
        an error of zero until the divergence of the ensemble is high enough
        again, resulting in a new segment start although it is just the
        start of the ensemble.

    Returns
    -------
    xarray.DataFrame with additional column "segment_start" (1 for start, else 0)
    """

    if cooldown == 0:

        def start(x, axis):
            return np.logical_and(
                x[:, :, :, :, 0] <= error_threshold, x[:, :, :, :, 1] > error_threshold
            )

        rolled = (
            df["Mean Squared Error"]
            .rolling(time_after_ascent=2, min_periods=2)
            .reduce(start)
            .fillna(False)
            .astype(dtype=bool)
        )

    else:

        def start(x, axis):
            logic = np.logical_and(
                x[:, :, :, :, -2] <= error_threshold,
                x[:, :, :, :, -1] > error_threshold,
            )
            max_val = min(np.shape(x)[-1], cooldown)
            if max_val > np.shape(x)[-1] - 2:
                max_val = np.shape(x)[-1] - 2
            for i in range(max_val):
                logic = np.logical_and(
                    logic,
                    np.logical_or(
                        np.isnan(x[:, :, :, :, i]), x[:, :, :, :, i] <= error_threshold
                    ),
                )
            return logic

        if cooldown + 2 > len(df["time_after_ascent"]):
            cooldown = len(df["time_after_ascent"]) - 2
        rolled = (
            df["Mean Squared Error"]
            .rolling(time_after_ascent=2 + cooldown, min_periods=2)
            .reduce(start)
            .fillna(False)
            .astype(dtype=bool)
        )
    df = df.assign(segment_start=rolled)
    df["segment_start"].attrs = {
        "standard_name": "segment_start",
        "long_name": "start of a segment",
        "auxiliary_data": "yes",
        "threshold": error_threshold,
    }
    return df


def predict_segments(
    df, threshold=0, threshold_jump=None, threshold_acc=None, repeats=1
):
    """
    Given sensitivities in a dataframe, predict when a segment of interest starts.

    Parameters
    ----------
    df: DataFrame
        Dataframe from find_segments() with sensitivities for every model state parameter and all trajectories
        and a column segment_start with boolean wether a segment starts.

    threshold: float or DataArray
        Threshold a gradient needs to overcome to be considered a segment.
        If DataArray, then dimensions need to have at least one dimension with df
        in common for defining thresholds along those dimensions.
    threshold_jump: float
        Threshold a gradient needs to jump in two consecutive
        time steps for a segment of interest (basically second order derivative).
        If DataArray, then dimensions need to have at least one dimension with df
        in common for defining thresholds along those dimensions.
    threshold_acc: float
        Threshold a gradient needs to accelerate (basically third order derivative).
        If DataArray, then dimensions need to have at least one dimension with df
        in common for defining thresholds along those dimensions.
    repeats: int
        Optional number of time steps at which any of the criteria must be
        fullfilled within window_size. Defaults to 1.

    Returns
    -------
    Dataframe with predicted segment starts.
    """
    if threshold is None and threshold_jump is None and threshold_acc is None:
        print("Need to have at least one threshold as criterium to find segments")
        return df

    def get_dims(threshold):
        dims = ()
        if "Output Parameter" in threshold.coords:
            dims = dims + (len(threshold.coords["Output Parameter"]),)
        else:
            dims = dims + (1,)
        if "Input Parameter" in threshold.coords:
            dims = dims + (len(threshold.coords["Input Parameter"]),)
        else:
            dims = dims + (1,)
        if "trajectory" in threshold.coords:
            dims = dims + (len(threshold.coords["trajectory"]),)
        else:
            dims = dims + (1,)
        return dims + (1,)

    def detect_segment_threshold(x, axis, threshold, repeats):
        # check if it is over the threshold when it wasn't a step before
        if isinstance(threshold, float) or isinstance(threshold, int):
            return np.logical_and(
                x[:, :, :, :, 0] <= threshold, x[:, :, :, :, 1] > threshold
            )
        dims = get_dims(threshold)
        return np.logical_and(
            x[:, :, :, :, 0] <= np.asarray(threshold).reshape(dims),
            x[:, :, :, :, 1] > np.asarray(threshold).reshape(dims),
        )

    def detect_segment_jump(x, axis, threshold_jump, repeats):
        if isinstance(threshold_jump, float) or isinstance(threshold_jump, int):
            return np.logical_and(
                x[:, :, :, :, 1] - x[:, :, :, :, 0] <= threshold_jump,
                x[:, :, :, :, 2] - x[:, :, :, :, 1] > threshold_jump,
            )
        dims = get_dims(threshold_jump)
        return np.logical_and(
            x[:, :, :, :, 1] - x[:, :, :, :, 0]
            <= np.asarray(threshold_jump).reshape(dims),
            x[:, :, :, :, 2] - x[:, :, :, :, 1]
            > np.asarray(threshold_jump).reshape(dims),
        )

    def detect_segment_acc(x, axis, threshold_acc, repeats):
        if isinstance(threshold_acc, float) or isinstance(threshold_acc, int):
            return np.logical_and(
                (x[:, :, :, :, 2] - x[:, :, :, :, 1])
                - (x[:, :, :, :, 1] - x[:, :, :, :, 0])
                <= threshold_acc,
                (x[:, :, :, :, 3] - x[:, :, :, :, 2])
                - (x[:, :, :, :, 2] - x[:, :, :, :, 1])
                > threshold_acc,
            )
        dims = get_dims(threshold_acc)
        return np.logical_and(
            (x[:, :, :, :, 2] - x[:, :, :, :, 1])
            - (x[:, :, :, :, 1] - x[:, :, :, :, 0])
            <= np.asarray(threshold_acc).reshape(dims),
            (x[:, :, :, :, 3] - x[:, :, :, :, 2])
            - (x[:, :, :, :, 2] - x[:, :, :, :, 1])
            > np.asarray(threshold_acc).reshape(dims),
        )

    if threshold is not None:
        window_size = 2 * repeats
        detected_segments = (
            df["Predicted Squared Error"]
            .rolling(time_after_ascent=window_size)
            .reduce(
                detect_segment_threshold, **{"threshold": threshold, "repeats": repeats}
            )
            .fillna(False)
            .astype(dtype=bool)
        )
        df = df.assign(detected_segment_thresh=detected_segments)
        df["detected_segment_thresh"].attrs = {
            "standard_name": "detected_segment_thresh",
            "long_name": "detected segment using threshold",
            "auxiliary_data": "yes",
            "threshold": threshold,
        }
    if threshold_jump is not None:
        window_size = 3 * repeats
        detected_segments = (
            df["Predicted Squared Error"]
            .rolling(time_after_ascent=window_size)
            .reduce(
                detect_segment_jump,
                **{"threshold_jump": threshold_jump, "repeats": repeats},
            )
            .fillna(False)
            .astype(dtype=bool)
        )
        df = df.assign(detected_segment_jump=detected_segments)
        df["detected_segment_jump"].attrs = {
            "standard_name": "detected_segment_jump",
            "long_name": "detected segment using threshold jump",
            "auxiliary_data": "yes",
            "threshold": threshold_jump,
        }
    if threshold_acc is not None:
        window_size = 4 * repeats
        detected_segments = (
            df["Predicted Squared Error"]
            .rolling(time_after_ascent=window_size)
            .reduce(
                detect_segment_acc,
                **{"threshold_acc": threshold_acc, "repeats": repeats},
            )
            .fillna(False)
            .astype(dtype=bool)
        )
        df = df.assign(detected_segment_acc=detected_segments)
        df["detected_segment_acc"].attrs = {
            "standard_name": "detected_segment_acc",
            "long_name": "detected segment using threshold acceleration",
            "auxiliary_data": "yes",
            "threshold": threshold_acc,
        }
    return df


def combine_predictions(df, def_ver=True, jum_ver=False, acc_ver=False, how=False):
    """
    Combine different predictions for a segment start.

    Parameters
    ----------
    df : Dataframe
        Dataframe with prediction columns such as
        "detected_segment_thresh", "detected_segment_jump"
        or "detected_segment_acc".
    def_ver : bool
        Use detection by reaching default threshold aka first derivative.
    jum_ver : bool
        Use detection by raching "jump" threshold aka second derivative.
    acc_ver : bool
        Use detection by reaching "acceleration" threshold aka third
        derivative
    how : bool
        If true, combine predictions using "and", otherwise "or".

    Returns
    -------
    Dataarray with combined predictions.
    """
    if how:
        final_pred = None
        if def_ver:
            final_pred = df["detected_segment_thresh"]

        if jum_ver:
            if final_pred is None:
                final_pred = df["detected_segment_jump"]
            else:
                final_pred = final_pred & df["detected_segment_jump"]

        if acc_ver:
            if final_pred is None:
                final_pred = df["detected_segment_acc"]
            else:
                final_pred = final_pred & df["detected_segment_acc"]

    else:
        final_pred = None
        if def_ver:
            final_pred = df["detected_segment_thresh"]

        if jum_ver:
            if final_pred is None:
                final_pred = df["detected_segment_jump"]
            else:
                final_pred = Union[final_pred, df["detected_segment_jump"]]

        if acc_ver:
            if final_pred is None:
                final_pred = df["detected_segment_acc"]
            else:
                final_pred = Union[final_pred, df["detected_segment_acc"]]
    return final_pred


def plot_histogram(df, column, bins=50, by=None, hist_log=True, drop_zero=True):
    """
    Plot a histogram over column for a dataframe with usually but not limited
    to either "segment_start" or "Predicted Segment Start".

    Parameters
    ----------
    df : Dataframe
        A dataframe with the column specified by column.
    bins : int
        Number of bins.
    by : string
        Column name to group by in different colors.
    hist_log : bool
        If true, use logarithmic x-axis.
    drop_zero : bool
        If true, drop rows where the value in column is zero.

    Returns
    -------
    Holoviews histogram plot
    """
    if column == "segment_start":
        hist_df = df.where(df[column], drop=True)
        title = f"Distribution of {np.sum(hist_df[column]).values} True Segment Starts over Predicted Squared Errors"
        column = "Predicted Squared Error"
    elif column == "Predicted Segment Start":
        hist_df = df.where(df[column], drop=True)
        title = f"Distribution of {np.sum(hist_df[column]).values} Predicted Segment Starts over Mean Squared Errors"
        column = "Mean Squared Error"
    else:
        title = f"Distribution of {column}"
        if drop_zero:
            hist_df = df.where(df[column] > 0.0, drop=True)
        else:
            hist_df = df

    if hist_log:
        min_x = df[column].values.flatten()[
            np.ma.masked_invalid(np.log(df[column])).argmin()
        ]
        return hist_df.hvplot.hist(
            column,
            by=by,
            bins=bins,
            alpha=0.3,
            width=1900,
            height=600,
            title=title,
            logx=hist_log,
            ylim=(0, 100),
            xlim=(min_x, df[column].max().values),
            stacked=False,
        )  # stacked would be nice but seems buggy
    else:
        return hist_df.hvplot.hist(
            column,
            by=by,
            bins=bins,
            alpha=0.3,
            width=1900,
            height=600,
            title=title,
            stacked=False,
        )


def fill_stats(
    def_val,
    jum_val,
    acc_val,
    index,
    df_input,
    how,
    step_tol,
    stats,
    def_ver,
    jum_ver,
    acc_ver,
):
    """
    Fill the row index in a (confusion + more info) matrix.
    The values are in the following order:
    True Positive: Number of true positive segment start predictions
    False Negative: Number of false negative predictions
    False Positive: Number of false positive predictions
    True Negative: Number of true negative predictions
    Early Positive: Number of true positive predictions that are earlier
        but within the time window or at exact the right time step
    Late Positive: Number of true positive predictions that are later
        but within the time window
    Positive Actual: Actual number of segment starts
    Positive Actual (Windows): Number of windows in which a segment starts
    Precision: (True Positive Predictions / Positive Predictions)
    Recall: (True Positive Predictions / Positive Actual)
    F1 Score: Balanced F-Score: (2 * (Precision*Recall) / (Precision+Recall))
    False Positive Rate: (False Positive / Negative Predictions)

    Parameters
    ----------
    def_val : float
        Threshold value to predict the start of a segment.
    jum_val : float
        Threshold value for second derivative to predict the start of a segment.
    acc_val : float
        Threshold value for third derivative to predict the start of a segment.
    index : int
        Row number to fill in stats.
    df_input : Dataframe
        Dataframe used to predict segments in and with a column
        "segment_start" for actual segment starts.
    how : bool
        If true, combine predictions using "and", otherwise "or".
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    stats : np.ndarray or list
        Matrix to fill a row in
    def_ver : bool
        Use detection by reaching default threshold aka first derivative.
    jum_ver : bool
        Use detection by raching "jump" threshold aka second derivative.
    acc_ver : bool
        Use detection by reaching "acceleration" threshold aka third
        derivative
    """
    df = predict_segments(df_input, def_val, jum_val, acc_val)

    final_pred = combine_predictions(df, def_ver, jum_ver, acc_ver, how)

    if step_tol == 0:
        tp = np.sum((df["segment_start"] == 1) & (final_pred == 1)).values
        tn = np.sum((df["segment_start"] == 0) & (final_pred == 0)).values
        fp = np.sum((df["segment_start"] == 0) & (final_pred == 1)).values
        fn = np.sum((df["segment_start"] == 1) & (final_pred == 0)).values
        ep = tp
        lp = 0
        p_act = np.sum(df["segment_start"] == 1)
        p_act_win = p_act
    else:
        idx_p_pred = np.argwhere(final_pred.values)
        idx_p = np.argwhere(df["segment_start"].values)
        idx_n_pred = np.argwhere(final_pred.values == 0)

        # Difference as in true start - predicted start
        idx_p_diff = []
        for pp in idx_p_pred:
            p = idx_p[(idx_p[:, 0:3] == pp[0:3]).all(axis=1)]
            if len(p) == 0:
                continue
            idx_p_diff.extend((p - pp)[:, -1])
        idx_p_diff = np.asarray(idx_p_diff)

        fn = 0
        p_act_win = 0
        last_idx = -np.inf
        last_inp = -1
        last_traj = -1
        last_outp = -1
        ep = 0
        lp = 0
        possible_window_size = 0
        for p in idx_p:
            if (
                p[3] - last_idx > step_tol
                or last_traj != p[2]
                or last_inp != p[1]
                or last_outp != p[0]
            ):
                possible_window_size = 2 * step_tol + 1
                p_act_win += possible_window_size
            else:
                possible_window_size = 2 * step_tol + 1 - (p[3] - last_idx - step_tol)
                p_act_win += possible_window_size
            last_idx = p[3]
            last_traj = p[2]
            last_inp = p[1]
            last_outp = p[0]
            # Check if a positive prediction does not exist
            pp = idx_p_pred[(idx_p_pred[:, 0:3] == p[0:3]).all(axis=1)]
            if len(pp) == 0:
                fn += possible_window_size  # There is no positive prediction
                continue
            # Check if positive prediction is within range
            differences = (p - pp)[:, -1]
            not_in_range = not ((np.abs(differences) <= step_tol)).any()
            if not_in_range > 0:
                fn += not_in_range * possible_window_size
            else:
                if ((differences >= -step_tol) & (differences <= 0)).any():
                    ep += 1
                else:
                    lp += 1
        # The version below is somewhat valid. If one segment is detected early and late,
        # it will be counted as two positive detections
        # The new version just takes the earlier detection and dimisses any late ones
        #         ep = np.sum( (idx_p_diff >= 0) & (idx_p_diff <= step_tol) )
        #         lp = np.sum( (idx_p_diff < 0) & (idx_p_diff >= -step_tol) )
        tp = ep + lp
        fp = np.sum(final_pred == 1) - tp
        # All negative windows = All possible windows - all positive windows
        tn = (
            (len(df["time_after_ascent"]) - (2 * step_tol + 1))
            * len(df["Output Parameter"])
            * len(df["Input Parameter"])
            * len(df["trajectory"])
            - np.sum(df["segment_start"]) * (1 + 2 * step_tol)
            - fn
        )  # All negative predictions - false negative
        p_act = np.sum(df["segment_start"] == 1)

    pr = tp / np.sum(final_pred)
    re = tp / np.sum(df["segment_start"])
    fpr = fp / np.sum((df["segment_start"] == 0))

    #  balanced F-score (from sklearn.metrics.f1_score)
    # The F1 score can be interpreted as a weighted average of the precision and recall,
    # where an F1 score reaches its best value at 1 and worst score at 0.
    # The relative contribution of precision and recall to the F1 score are equal.
    # F1 = 2 * (precision * recall) / (precision + recall)
    f1 = 2 * (pr * re) / (pr + re)
    stats[index] = np.array([tp, fn, fp, tn, ep, lp, p_act, p_act_win, pr, re, f1, fpr])


def get_stats(
    df_input,
    def_ver,
    jum_ver,
    acc_ver,
    how,
    def_thresh=None,
    jum_thresh=None,
    acc_thresh=None,
    n=100,
    step_tol=0,
):
    """
    Get statistics in a (confusion + more info) matrix. Each column consists of
    values in the following order:
    True Positive: Number of true positive segment start predictions
    False Negative: Number of false negative predictions
    False Positive: Number of false positive predictions
    True Negative: Number of true negative predictions
    Early Positive: Number of true positive predictions that are earlier
        but within the time window or at exact the right time step
    Late Positive: Number of true positive predictions that are later
        but within the time window
    Positive Actual: Actual number of segment starts
    Positive Actual (Windows): Number of windows in which a segment starts
    Precision: (True Positive Predictions / Positive Predictions)
    Recall: (True Positive Predictions / Positive Actual)
    F1 Score: Balanced F-Score: (2 * (Precision*Recall) / (Precision+Recall))
    False Positive Rate: (False Positive / Negative Predictions)

    Each row is used for a different combination of thresholds for predicting
    segment starts. If *all* combinations are needed,
    use get_stats_combinations(..).

    Parameters
    ----------
    df_input : Dataframe
        Dataframe used to predict segments in and with a column
        "segment_start" for actual segment starts and "time_after_ascent".
    def_ver : bool
        Use detection by reaching default threshold aka first derivative.
    jum_ver : bool
        Use detection by raching "jump" threshold aka second derivative.
    acc_ver : bool
        Use detection by reaching "acceleration" threshold aka third
        derivative
    how : bool
        If true, combine predictions using "and", otherwise "or".
    def_thresh : float or int
        Maximum threshold value to predict the start of a segment.
    jum_thresh : float or int
        Maximum threshold value for second derivative to predict the start of a segment.
    acc_thresh : float or int
        Maximum threshold value for third derivative to predict the start of a segment.
    n : int
        Number of rows in the return matrix aka the number of different
        thresholds used in predicting segments.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.

    Returns
    -------
    2D np.array with rows for each threshold and columns as described above.
    """
    stats = np.zeros((n, 12))
    max_time = np.shape(df_input["time_after_ascent"])[-1]

    for i in range(n):
        if def_ver:
            if isinstance(def_thresh, int) or isinstance(def_thresh, float):
                if def_thresh > 0:
                    def_val = i / def_thresh
                else:
                    def_val = 0
            else:
                def_val = (i / def_thresh).fillna(0)
        else:
            def_val = None

        if jum_ver:
            if isinstance(jum_thresh, int) or isinstance(jum_thresh, float):
                if jum_thresh > 0:
                    jum_val = i / jum_thresh
                else:
                    jum_val = 0
            else:
                jum_val = (i / jum_thresh).fillna(0)
        else:
            jum_val = None

        if acc_ver:
            if isinstance(acc_thresh, int) or isinstance(acc_thresh, float):
                if acc_thresh > 0:
                    acc_val = i / acc_thresh
                else:
                    acc_val = 0
            else:
                acc_val = (i / acc_thresh).fillna(0)
        else:
            acc_val = None

        fill_stats(
            def_val,
            jum_val,
            acc_val,
            i,
            df_input,
            how,
            step_tol,
            stats,
            def_ver,
            jum_ver,
            acc_ver,
        )
    return stats


def confusion_matrix(
    df,
    def_ver,
    jum_ver,
    acc_ver,
    how,
    def_thresh=None,
    jum_thresh=None,
    acc_thresh=None,
    n=100,
    step_tol=0,
    confus=None,
):
    """
    Plot a confusion matrix with F1 score and step tolerance in the title.
    The best result for true positive predictions is used for plotting.

    Parameters
    ----------
    df : Dataframe
        Dataframe used to predict segments in and with a column
        "segment_start" for actual segment starts and "time_after_ascent".
    def_ver : bool
        Use detection by reaching default threshold aka first derivative.
    jum_ver : bool
        Use detection by raching "jump" threshold aka second derivative.
    acc_ver : bool
        Use detection by reaching "acceleration" threshold aka third
        derivative
    how : bool
        If true, combine predictions using "and", otherwise "or".
    def_thresh : float or int
        Maximum threshold value to predict the start of a segment.
    jum_thresh : float or int
        Maximum threshold value for second derivative to predict the start of
        a segment.
    acc_thresh : float or int
        Maximum threshold value for third derivative to predict the start of
        a segment.
    n : int
        Number of rows in the return matrix aka the number of different
        thresholds used in predicting segments.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    confus : 2D np.array
        Confusion matrix from get_stats(..). If None is given, calculate
        the confusion matrix defined using the other parameters.

    Returns
    -------
    Holoviews heatmap plot, confusion matrix (np.ndarray)
    """
    if confus is None:
        confus = get_stats(
            df,
            def_ver,
            jum_ver,
            acc_ver,
            how,
            def_thresh,
            jum_thresh,
            acc_thresh,
            n,
            step_tol,
        )
    pos = [(0, 1), (1, 1), (0, 0), (1, 0), (0, -1), (1, -1), (0, -2), (1, -2)]
    labels = [
        "True Positive",
        "False Negative",
        "False Positive",
        "True Negative",
        "Early Positive",
        "Late Positive",
        "Total Positive (Actual)",
        "Total Positive (Windows)",
    ]
    best_data = confus[np.argmax(confus[:, 0])]
    heatmap_data = list(map(lambda tup, i: tup + (i,), pos, best_data[0:8]))
    labels = list(
        map(lambda tup, lab, i: tup + (f"{lab}\n{i}",), pos, labels, best_data[0:8])
    )
    plot = hv.HeatMap(heatmap_data).opts(
        cmap=cm.Blues,
        title=f"Best Result for TP; F1 score: {best_data[-2]:.2f}, Time Tolerance: +/-{step_tol*20}s",
        colorbar=False,
        xaxis="bare",
        yaxis="bare",
        height=640,
        width=640,
        logz=True,
    )
    return plot * hv.Labels(labels), confus


def AUC(con_mat, extra_title=None):
    """
    Plot false positive over true positive.
    As a rule of thumb, if the curve is very steep
    very early on, the data is rather easy to classify.

    Parameters
    ----------
    con_mat : 2D np.array
        Confusion matrix created by get_stats(..).
    extra_title : string
        Add this to the title.

    Returns
    -------
    holoviews.Curve with AUC
    """
    title = "Area Under The Curve"
    if extra_title is not None:
        title += "\n" + extra_title
    return hv.Curve(con_mat[:, [2, 0]]).opts(
        xlabel="False Positive", ylabel="True Positive", width=1000, title=title
    )


def PRC(con_mat, extra_title=None):
    """
    Plot the Precision-Recall Curve (re over pr)

    Parameters
    ----------
    con_mat : 2D np.array
        Confusion matrix created by get_stats(..).
    extra_title : string
        Add this to the title.

    Returns
    -------
    holoviews.Curve with PRC
    """
    title = "Precision-Recall Curve"
    if extra_title is not None:
        title += "\n" + extra_title
    return hv.Curve(con_mat[:, [7, 6]]).opts(
        xlabel="Recall", ylabel="Precision", width=1000, title=title
    )


def ROC(con_mat, extra_title=None):
    """
    Plot receiver operating characteristic curve, that is
    true positive rate (=recall) against false positive rate with
    tpr = tp/p
    fpr = tn/n

    Parameters
    ----------
    con_mat : 2D np.array
        Confusion matrix created by get_stats(..).
    extra_title : string
        Add this to the title.

    Returns
    -------
    holoviews.Curve with ROC
    """
    title = "Receiver Operating Characteristics"
    if extra_title is not None:
        title += "\n" + extra_title
    return (
        hv.Curve(con_mat[:, [9, 7]]) * hv.Curve([[0, 0], [1, 1]]).opts(color="#D3D3D3")
    ).opts(
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        width=1000,
        xlim=(0, 1),
        ylim=(0, 1),
        title=title,
    )


def get_stats_combinations(
    df_input,
    def_ver,
    jum_ver,
    acc_ver,
    how,
    n=100,
    step_tol=0,
    limits=[(-80, 0), (-80, 0), (-80, 0)],
):
    """
    Get statistics in a (confusion + more info) matrix. Each column consists of
    values in the following order:
    True Positive: Number of true positive segment start predictions
    False Negative: Number of false negative predictions
    False Positive: Number of false positive predictions
    True Negative: Number of true negative predictions
    Early Positive: Number of true positive predictions that are earlier
        but within the time window or at exact the right time step
    Late Positive: Number of true positive predictions that are later
        but within the time window
    Positive Actual: Actual number of segment starts
    Positive Actual (Windows): Number of windows in which a segment starts
    Precision: (True Positive Predictions / Positive Predictions)
    Recall: (True Positive Predictions / Positive Actual)
    F1 Score: Balanced F-Score: (2 * (Precision*Recall) / (Precision+Recall))
    False Positive Rate: (False Positive / Negative Predictions)

    Each row is used for *all* combinations of thresholds for predicting
    segment starts. The rows are ordered by default threshold, jump threshold,
    acceleration threshold (fastest index). There are n*n*n many rows
    if all threshold types are used, n*n many rows if two types are used and
    n if only one is used.

    Parameters
    ----------
    df_input : Dataframe
        Dataframe used to predict segments in and with a column
        "segment_start" for actual segment starts and "time_after_ascent".
    def_ver : bool
        Use detection by reaching default threshold aka first derivative.
    jum_ver : bool
        Use detection by raching "jump" threshold aka second derivative.
    acc_ver : bool
        Use detection by reaching "acceleration" threshold aka third
        derivative
    how : bool
        If true, combine predictions using "and", otherwise "or".
    n : int
        Number of rows in the return matrix aka the number of different
        thresholds used in predicting segments.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    limits : List of tuples
        Upper and lower limits for each threshold type. Threshold using
        def_ver is at index 0, jum_ver at index 1 and acc_ver at index 2.

    Returns
    -------
    2D np.array with rows for each threshold combination
    and columns as described above.
    """
    total_n = 1
    delta_def = None
    if def_ver:
        total_n *= n
        delta_def = (limits[0][1] - limits[0][0]) / (n - 1)
    delta_jum = None
    if jum_ver:
        total_n *= n
        delta_jum = (limits[1][1] - limits[1][0]) / (n - 1)
    delta_acc = None
    if acc_ver:
        total_n *= n
        delta_acc = (limits[2][1] - limits[2][0]) / (n - 1)
    stats = np.zeros((total_n, 12))
    max_time = np.shape(df_input["time_after_ascent"])[-1]

    for i in range(n):
        if def_ver:
            def_val = 10 ** (i * delta_def + limits[0][0])
            if not jum_ver and not acc_ver:
                fill_stats(
                    def_val,
                    None,
                    None,
                    i,
                    df_input,
                    how,
                    step_tol,
                    stats,
                    def_ver,
                    jum_ver,
                    acc_ver,
                )
            elif jum_ver:
                for j in range(n):
                    jum_val = 10 ** (j * delta_jum + limits[1][0])
                    if not acc_ver:
                        fill_stats(
                            def_val,
                            jum_val,
                            None,
                            i * n + j,
                            df_input,
                            how,
                            step_tol,
                            stats,
                            def_ver,
                            jum_ver,
                            acc_ver,
                        )
                    else:
                        for k in range(n):
                            acc_val = 10 ** (k * delta_acc + limits[2][0])
                            fill_stats(
                                def_val,
                                jum_val,
                                acc_val,
                                i * n * n + j * n + k,
                                df_input,
                                how,
                                step_tol,
                                stats,
                                def_ver,
                                jum_ver,
                                acc_ver,
                            )
            elif acc_ver:
                for j in range(n):
                    acc_val = 10 ** (j * delta_acc + limits[2][0])
                    fill_stats(
                        def_val,
                        None,
                        acc_val,
                        i * n + j,
                        df_input,
                        how,
                        step_tol,
                        stats,
                        def_ver,
                        jum_ver,
                        acc_ver,
                    )
        elif jum_ver:
            jum_val = 10 ** (i * delta_jum + limits[1][0])
            if not acc_ver:
                fill_stats(
                    None,
                    jum_val,
                    None,
                    i,
                    df_input,
                    how,
                    step_tol,
                    stats,
                    def_ver,
                    jum_ver,
                    acc_ver,
                )
            else:
                for j in range(n):
                    acc_val = 10 ** (j * delta_acc + limits[2][0])
                    fill_stats(
                        None,
                        jum_val,
                        acc_val,
                        i * n + j,
                        df_input,
                        how,
                        step_tol,
                        stats,
                        def_ver,
                        jum_ver,
                        acc_ver,
                    )
        else:
            acc_val = 10 ** (i * delta_acc + limits[2][0])
            fill_stats(
                None,
                None,
                acc_val,
                i,
                df_input,
                how,
                step_tol,
                stats,
                def_ver,
                jum_ver,
                acc_ver,
            )

    return stats


def create_chached_matrix_dic(
    segment_data,
    segment_threshold=10 ** (-10.3),
    steps=21,
    min_def_thresh=-80,
    max_def_thresh=0,
    min_jum_thresh=-80,
    max_jum_thresh=0,
    min_acc_thresh=-80,
    max_acc_thres=0,
    step_tol=2,
    cooldown=0,
):
    """

    Parameters
    ----------
    segment_data : Dataframe
        Dataframe with MSE created by load_dataset() for every model state
        parameter and all trajectories.
    segment_threshold : float
        Threshold for errors to identify true segment starts.
    steps : int
        Number of steps from minimum to maximum thresholds to calculate
        entries in the confusion matrix for.
    min_def_thresh : float
        Minimum threshold value to predict the start of a segment.
    max_def_thresh : float
        Maximum threshold value to predict the start of a segment.
    min_jum_thresh : float
        Minimum threshold value for second derivative to predict the start of a segment.
    max_jum_thresh : float
        Maximum threshold value for second derivative to predict the start of a segment.
    min_acc_thresh : float
        Minimum threshold value for third derivative to predict the start of a segment.
    max_acc_thres : float
        Maximum threshold value for third derivative to predict the start of a segment.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    cooldown : int
        Minimum number of timesteps between last time the error threshold has
        been met and the next time step. This is useful if the original data
        is based on ensembles that start at every few time step which leads to
        an error of zero until the divergence of the ensemble is high enough
        again, resulting in a new segment start although it is just the
        start of the ensemble.

    Returns
    -------
    Dictionary with confusion matrices for every output parameter.
    """
    cached_data = find_segments(segment_data, segment_threshold, cooldown)
    delta_def = (max_def_thresh - min_def_thresh) / (steps - 1)
    delta_jum = (max_jum_thresh - min_jum_thresh) / (steps - 1)
    cached_matrix_dic = {}
    for out_p in np.unique(cached_data["Output Parameter"]):
        cached_matrix_dic[out_p] = get_stats_combinations(
            df_input=cached_data.sel({"Output Parameter": [out_p]}),
            def_ver=True,
            jum_ver=True,
            acc_ver=False,
            how=False,
            n=steps,
            step_tol=step_tol,
            limits=[
                (min_def_thresh, max_def_thresh),
                (min_jum_thresh, max_jum_thresh),
                (-80, 0),
            ],
        )
    return cached_matrix_dic


def confusion_matrix_faster(index, confus):
    """
    Plot a confusion matrix with F1 score and step tolerance in the title.
    The best result for true positive predictions is used for plotting.
    This is supposedly faster than confusion_matrix(..) where
    the confusion matrix is already calculated.

    Parameters
    ----------
    index : int
        Row index of the confusion matrix to plot.
    confus : 2D np.array
        Confusion matrix created by get_stats(..).

    Returns
    -------
    Holoviews heatmap, np.array with row given by index.
    """
    pos = [(0, 1), (1, 1), (0, 0), (1, 0), (0, -1), (1, -1), (0, -2), (1, -2)]
    labels = [
        "True Positive",
        "False Negative",
        "False Positive",
        "True Negative",
        "Early Positive",
        "Late Positive",
        "Total Positive (Actual)",
        "Total Positive (Windows)",
    ]
    best_data = confus[index]
    heatmap_data = list(map(lambda tup, i: tup + (i,), pos, best_data[0:8]))
    labels = list(
        map(lambda tup, lab, i: tup + (f"{lab}\n{i}",), pos, labels, best_data[0:8])
    )
    min_c = 0
    max_c = np.max([np.max(best_data[0:1]), np.max(best_data[4:8])])
    plot = hv.HeatMap(heatmap_data).opts(
        cmap=cm.Blues,
        title=f"F1 score: {best_data[-2]:.2f}, Time Tolerance: +/-{step_tol*20}s",
        colorbar=False,
        clim=(min_c, max_c),
        xaxis="bare",
        yaxis="bare",
        height=640,
        width=640,
        logz=False,
    )
    return plot * hv.Labels(labels), best_data


def create_df_confusion(
    df,
    segment_threshold=10 ** (-10.3),
    steps=21,
    min_def_thresh=-80,
    max_def_thresh=0,
    min_jum_thresh=-80,
    max_jum_thresh=0,
    min_acc_thresh=-80,
    max_acc_thres=0,
    step_tol=2,
    how=False,
    cooldown=0,
):
    """
    Create an xr.Dataset using get_stats_combinations(..) for
    different output parameters, input parameters and thresholds.

    Parameters
    ----------
    df : Dataframe
        Dataframe used to predict segments in and with a column
        "segment_start" for actual segment starts and "time_after_ascent".
    segment_threshold : float
        Threshold for errors to identify true segment starts.
    steps : int
        Number of steps from minimum to maximum thresholds to calculate
        entries in the confusion matrix for.
    min_def_thresh : float
        Minimum threshold value to predict the start of a segment.
    max_def_thresh : float
        Maximum threshold value to predict the start of a segment.
    min_jum_thresh : float
        Minimum threshold value for second derivative to predict the start of a segment.
    max_jum_thresh : float
        Maximum threshold value for second derivative to predict the start of a segment.
    min_acc_thresh : float
        Minimum threshold value for third derivative to predict the start of a segment.
    max_acc_thres : float
        Maximum threshold value for third derivative to predict the start of a segment.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    how : bool
        If true, combine predictions using "and", otherwise "or".
    cooldown : int
        Minimum number of timesteps between last time the error threshold has
        been met and the next time step. This is useful if the original data
        is based on ensembles that start at every few time step which leads to
        an error of zero until the divergence of the ensemble is high enough
        again, resulting in a new segment start although it is just the
        start of the ensemble.

    Returns
    -------
    xr.Dataset for statistics of predictions.
    Coordinates are "Output Parameter", "Input Parameter", "Default Threshold",
    "Jump Threshold" ("Acceleration Threshold" is not supported).
    Columns are "TP over P", "FP over P", "FN over P_windows", "TP", "P", "FP",
    "FN", "P_windows".
    """
    cached_data = find_segments(df, segment_threshold, cooldown)
    delta_def = (max_def_thresh - min_def_thresh) / (steps - 1)
    delta_jum = (max_jum_thresh - min_jum_thresh) / (steps - 1)

    out_params = list(np.unique(cached_data["Output Parameter"]))
    in_params = list(np.unique(cached_data["Input Parameter"]))
    # These are the thresholds calculated in get_stats_combinations
    def_thresh = np.arange(min_def_thresh, max_def_thresh + delta_def / 2, delta_def)
    jum_thresh = np.arange(min_jum_thresh, max_jum_thresh + delta_jum / 2, delta_jum)
    coords_dic = {
        "Output Parameter": out_params,
        "Input Parameter": in_params,
        "Default Threshold": def_thresh,
        "Jump Threshold": jum_thresh,
    }
    tpp = np.zeros((len(out_params), len(in_params), len(def_thresh), len(jum_thresh)))
    tp = np.zeros((len(out_params), len(in_params), len(def_thresh), len(jum_thresh)))
    p = np.zeros((len(out_params), len(in_params), len(def_thresh), len(jum_thresh)))
    fp = np.zeros((len(out_params), len(in_params), len(def_thresh), len(jum_thresh)))
    fn = np.zeros((len(out_params), len(in_params), len(def_thresh), len(jum_thresh)))
    pw = np.zeros((len(out_params), len(in_params), len(def_thresh), len(jum_thresh)))
    for i, out_p in enumerate(out_params):
        df_out_p = cached_data.sel({"Output Parameter": [out_p]})
        for j, in_p in enumerate(in_params):
            matrix = get_stats_combinations(
                df_input=df_out_p.sel({"Input Parameter": [in_p]}),
                def_ver=True,
                jum_ver=True,
                acc_ver=False,
                how=how,
                n=steps,
                step_tol=step_tol,
                limits=[
                    (min_def_thresh, max_def_thresh),
                    (min_jum_thresh, max_jum_thresh),
                    (-80, 0),
                ],
            )
            matrix = np.reshape(
                matrix,
                (
                    int(np.sqrt(np.shape(matrix)[0])),
                    int(np.sqrt(np.shape(matrix)[0])),
                    12,
                ),
            )
            tp[i, j, :, :] = matrix[:, :, 0]
            fn[i, j, :, :] = matrix[:, :, 1]
            fp[i, j, :, :] = matrix[:, :, 2]
            p[i, j, :, :] = matrix[:, :, 6]
            pw[i, j, :, :] = matrix[:, :, 7]
            tpp[i, j, :, :] = matrix[:, :, 9]

    dims = [
        "Output Parameter",
        "Input Parameter",
        "Default Threshold",
        "Jump Threshold",
    ]
    data_dic = {
        "TP over P": (dims, tpp),
        "FP over P": (dims, fp / p),
        "FN over P_windows": (dims, fn / pw),
        "TP": (dims, tp),
        "P": (dims, p),
        "FP": (dims, fp),
        "FN": (dims, fn),
        "P_windows": (dims, pw),
    }
    return xr.Dataset(data_vars=data_dic, coords=coords_dic)


def create_bar_plot(df, use_percentage, width=1000, height=600, extra_title=""):
    """
    Create a bar plot for different columns of the confusion matrix.

    Parameters
    ----------
    df : xarray.Dataset
        Confusion matrix created by create_df_confusion(..).
    use_percentage : bool
        If true, plot percentages for the bar plot. Otherwise plot
        absolute values.
    width : int
        Width in pixels of the plot.
    height : int
        Height in pixels of the plot.
    extra_title : string
        Addition to the title.

    Returns
    -------
    Holoviews barplot.
    """
    if use_percentage:
        if "P_windows" in df:
            by = ["TP over P", "FP over P", "FN over P_windows"]
        else:
            by = ["TP over P", "FP over P", "FN over P"]
        return df.hvplot.bar(
            x="Output Parameter",
            y=by,
            stacked=False,
            rot=90,
            width=width,
            height=height,
            ylim=(0, 1),
            cmap={
                "TP over P": "seagreen",
                "FP over P": "crimson",
                "FN over P_windows": "royalblue",
            },
        ).opts(title="Prediction Count (%)" + extra_title)
    else:
        if "P_windows" in df:
            by = ["TP", "P", "FP", "FN", "P_windows"]
        else:
            by = ["TP", "P", "FP", "FN"]
        return df.hvplot.bar(
            x="Output Parameter",
            y=by,
            stacked=False,
            rot=90,
            width=width,
            height=height,
            cmap={
                "P_windows": "mediumorchid",
                "TP": "seagreen",
                "FP": "crimson",
                "FN": "royalblue",
                "P": "orange",
            },
        ).opts(title="Prediction Count" + extra_title)


def create_many_bar_plots(df, cols=2, width=500, height=400):
    """
    Plot multiple bar plots (for each input parameter one).

    Parameters
    ----------
    df : Dataframe
        Dataframe used to predict segments in and with a column
        "segment_start" for actual segment starts and "time_after_ascent".
    cols : int
        Number of plots per row.
    width : int
        Width in pixels of the plot.
    height : int
        Height in pixels of the plot.

    Returns
    -------
    Multiple Holoviews barplots.
    """
    plots = None
    for in_p in np.unique(df["Input Parameter"]):
        if plots is not None:
            plots = plots + create_bar_plot(
                df.sel({"Input Parameter": in_p}),
                False,
                width=width,
                height=height,
                extra_title=" for " + in_p,
            )
        else:
            plots = create_bar_plot(
                df.sel({"Input Parameter": in_p}),
                False,
                width=width,
                height=height,
                extra_title=" for " + in_p,
            )
        plots = plots + create_bar_plot(
            df.sel({"Input Parameter": in_p}), True, width=width, height=height
        )
    return plots.opts(shared_axes=False).cols(cols)


def create_input_labels(ds, step_tol, distinct_outparams, verbosity=0):
    """
    Create windows and labels for creating training and testing sets.

    Parameters
    ----------
    ds : Dataset
        Dataset created by find_segments(..) where segments are identified and
        predicted errors are stored.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    distinct_outparams : bool
        If true, try to predict perturbing an input parameter for a segment
        start for each output parameter independently. This may not be useful,
        since one input parameter can have a high impact on multiple output
        parameters.
    verbosity : int
        Set verbosity level.

    Returns
    -------
    X as np.array of shape (n_samples, features),
    labels as np.array of shape (n_samples, n_classes),
    feature_names, class_names
    """
    # what I want:
    # (n_timewindows, window_size*output_param*input_param), where the latter consists of [grad_q1_d1_t-2, t-1, ..., t+2, d2_t-2, ...]
    # Reshape input data
    # The input set has a feature dimension that consists of:
    # window size many gradients concatenated for every gradient
    # and optional
    # gradients for every output parameter
    # [dqv/dx_1|_1, dqv/dx_1|_2, dqv/dx_1|_3.
    #  dqc/dx_1|_1, dqc/dx_1|_2, dqc/dx_1|_3
    #  dqv/dx_2|_1, dqv/dx_2|_2, dqv/dx_2|_3.
    #  dqc/dx_2|_1, dqc/dx_2|_2, dqc/dx_2|_3]
    # old:
    # [dqv/dx_1|_1, dqv/dx_1|_2, dqv/dx_1|_3.
    #  dqv/dx_2|_1, dqv/dx_2|_2, dqv/dx_2|_3
    #  dqc/dx_1|_1, dqc/dx_1|_2, dqc/dx_1|_3.
    #  dqc/dx_2|_1, dqc/dx_2|_2, dqc/dx_2|_3]
    # The labels consist of input parameters many
    # outputs. When perturbing the parameter
    # it doesn't matter which output parameter
    # has the highest real error.
    n_input_params = len(ds["Input Parameter"])
    n_timesteps = len(ds["time_after_ascent"])
    n_out_params = len(ds["Output Parameter"])
    n_trajs = len(ds["trajectory"])
    # (due to padding)
    n_windows = n_timesteps - step_tol
    X = ds["Predicted Squared Error"].values
    X = np.moveaxis(X, [0], [1])
    # pad left zero values
    X_shape = np.shape(X)[0:3] + (step_tol,)
    X = np.insert(X, [0], np.zeros(X_shape), axis=3)
    y = ds["segment_start"].values
    y = np.moveaxis(y, [0], [1])
    # pad right zero values
    y_shape = np.shape(y)[0:3] + (step_tol,)
    y = np.insert(y, [-1], np.zeros(y_shape), axis=3)

    if distinct_outparams:
        feature_names = product(
            np.unique(ds["Input Parameter"].values),
            np.unique(ds["Output Parameter"].values),
        )
        class_names = np.unique(ds["Input Parameter"].values)
        X_final = np.concatenate(
            [
                *(
                    np.vstack(
                        X[:, out_p, traj, j : (step_tol * 2 + 1) + j]
                        for j in range(n_windows)
                    ).reshape((n_windows, -1))
                    for traj in range(n_trajs)
                    for out_p in range(n_out_params)
                )
            ],
            axis=0,
        )
        y_final = np.concatenate(
            [
                *(
                    np.vstack(
                        y[:, out_p, traj, j : (step_tol * 2 + 1) + j]
                        for j in range(n_windows)
                    ).reshape((n_windows, -1))
                    for traj in range(n_trajs)
                    for out_p in range(n_out_params)
                )
            ],
            axis=0,
        )
        y_final = np.reshape(y_final, (-1, n_input_params, (2 * step_tol + 1)))
        class_names = [""]
    else:
        feature_names = np.unique(ds["Input Parameter"].values)
        class_names = np.unique(ds["Input Parameter"].values)
        X_final = np.concatenate(
            [
                *(
                    np.vstack(
                        X[:, :, traj, j : (step_tol * 2 + 1) + j].reshape(
                            (-1, step_tol * 2 + 1)
                        )
                        for j in range(n_windows)
                    ).reshape((n_windows, -1))
                    for traj in range(n_trajs)
                )
            ],
            axis=0,
        )
        y_final = np.concatenate(
            [
                *(
                    np.vstack(
                        y[:, :, traj, j : (step_tol * 2 + 1) + j].reshape(
                            (-1, step_tol * 2 + 1)
                        )
                        for j in range(n_windows)
                    ).reshape((n_windows, -1))
                    for traj in range(n_trajs)
                )
            ],
            axis=0,
        )
        y_final = np.reshape(
            y_final, (-1, n_input_params, (2 * step_tol + 1) * n_out_params)
        )

    y_final = np.sum(y_final, axis=2) > 0
    # Binary classification still needs two dimension
    # where it is either it or not
    if np.shape(y_final)[-1] == 1:
        y_final = np.hstack(y_final, (not y_final))
        # print(f"Shape of y: {np.shape(y_final)}")

    if verbosity > 5:
        sumsies = np.sum(y_final, axis=0)
        print(
            f"Number of true per class (59 classes; shape: {np.shape(sumsies)}): {sumsies}"
        )
        for c, s in zip(class_names, sumsies):
            print(f"{c}: {s}")
    return X_final, y_final, feature_names, class_names


def create_forest(
    ds,
    step_tol,
    test_size=None,
    distinct_outparams=False,
    n_estimators=100,
    max_features="auto",
    save_memory=False,
    no_split=False,
    verbosity=False,
    max_depth=36,
    max_leaf_nodes=1000,
    classifier="random_forest",
    learning_rate=1.0,
):
    """
    Create a random forest and fit it. Returns training and testing set and
    feature names and class names.

    Parameters
    ----------
    ds : Dataset or list of label and data.
        Dataset created by find_segments(..) where segments are identified and
        predicted errors are stored.
        If list, then a list of precalculated labels and their data.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    test_size : float
        Percentage of number of windows to use for testing.
    distinct_outparams : bool
        If true, try to predict perturbing an input parameter for a segment
        start for each output parameter independently. This may not be useful,
        since one input parameter can have a high impact on multiple output
        parameters.
    n_estimators : int
        Number of trees.
    max_features : string or float or int
        From sklearn.ensemble.RandomForestClassifier:
        The number of features to consider when looking for the best split:

        If int, then consider max_features features at each split.
        If float, then max_features is a fraction and
            round(max_features * n_features) features are considered at each split.
        If “auto”, then max_features=sqrt(n_features).
        If “sqrt”, then max_features=sqrt(n_features) (same as “auto”).
        If “log2”, then max_features=log2(n_features).
        If None, then max_features=n_features.

        Note: the search for a split does not stop until at least one valid
        partition of the node samples is found, even if it requires to
        effectively inspect more than max_features features.
    save_memory : bool
        If true and not "no_split", use the first "test_size" many rows for
        training. This means, that the data is not stratified! Do not use this
        unless the classes to predict is roughly uniformly distributed.
    no_split : bool
        If true, do not split data into training and test set and use everything
        for training.
    verbosity : bool
        If true, set verbosity of the random forest to 2.
    max_depth : int
        From sklearn.ensemble.RandomForestClassifier:
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split (here: 2) samples.
    max_leaf_nodes : int
        From sklearn.ensemble.RandomForestClassifier:
        Grow trees with max_leaf_nodes in best-first fashion. Best nodes are
        defined as relative reduction in impurity. If None then unlimited
        number of leaf nodes.
    classifier : string
        Choose a classifier for prediction. Options are "random_forest"
        and "adaboost". The latter ignores "max_leaf_nodes" and "max_depth" and
        rather uses "learning_rate".
    learning_rate : float
        From sklean.ensemble.AdaBoostClassifier:
        Learning rate shrinks the contribution of each classifier by
        "learning_rate". There is a trade-off between "learning_rate" and
        "n_estimators".

    Returns
    -------
    model, train set, test set, train labels, test labels, feature names,
    class names
    """
    # Distinct outparams *should* lead to worse
    # predictions since the forest cannot differentiate
    # if the current window shows gradients for QV or QC
    # but if it is robust in this version, it might be a
    # hint for a generalized formulation
    verbose = 0
    if verbosity:
        verbose = 2

    if classifier == "random_forest":
        model = OneVsRestClassifier(
            RandomForestClassifier(
                n_estimators=n_estimators,
                max_features=max_features,
                n_jobs=-1,
                verbose=verbose,
                max_depth=max_depth,
                max_leaf_nodes=max_leaf_nodes,
            ),
            n_jobs=1,
        )
    else:
        model = OneVsRestClassifier(
            AdaBoostClassifier(
                n_estimators=n_estimators,
                random_state=42,
                learning_rate=learning_rate,
            ),
            n_jobs=1,  # Using more jobs here requires a lot of memory
        )

    if isinstance(ds, list):
        y_final, X_final = ds
        model = model.fit(X_final, y_final)
        return model, X_final, None, y_final, None, None, None
    else:
        X_final, y_final, feature_names, class_names = create_input_labels(
            ds, step_tol, distinct_outparams
        )

    if no_split:
        model = model.fit(X_final, y_final)
        return model, X_final, None, y_final, None, feature_names, class_names

    if save_memory:
        model = model.fit(
            X_final[: int(np.shape(X_final)[0] * (1 - test_size))],
            y_final[: int(np.shape(X_final)[0] * (1 - test_size))],
        )
        return (
            model,
            X_final[0 : int(np.shape(X_final)[0] * (1 - test_size))],
            X_final[int(np.shape(X_final)[0] * (1 - test_size)) : :],
            y_final[0 : int(np.shape(X_final)[0] * (1 - test_size))],
            y_final[int(np.shape(X_final)[0] * (1 - test_size)) : :],
            feature_names,
            class_names,
        )
    else:
        train, train_labels, test, test_labels = iterative_train_test_split(
            X_final, y_final, test_size=test_size
        )

    model = model.fit(train, train_labels)
    return model, train, test, train_labels, test_labels, feature_names, class_names


def get_tree_matrix(trained_model, X, y, only_idx=None, step_tol=None):
    """
    Create a confusion matrix with a single row given a random forest.
    Late and early positive is every segment start, that got detected within
    the given tolerance incl. exact true predictions.
    If "step_tol" is given, then "True Positive" means, for any segment
    start there is at least one window (early or late but within tolerance)
    which detected that. "step_tol" can be any value, since the value
    itself is given by the window size.
    If "step_tol" is not given, then "True Positive" is counted for each
    window individually. A segment start "exists" as many times as the
    window size such that "True Positive" can be higher than actual segment
    starts exist.
    "Early Positive" and "Late Positive" are set to zero.

    Parameters
    ----------
    trained_model : model
        A trained model with a method predict(..).
    X : Array or list
        Array of shape (n_samples, n_features) for prediction.
    y : Array or list
        Array of shape (n_samples, n_classes) with labels associated with X.
    only_idx : int
        Get values in confusion matrix only for the input parameter at the
        given index in n_features.
    step_tol : int or None
        If int is given: Number of steps to tolerate for a prediction
        to be true. Otherwise only look for exact matches.
    Returns
    -------
    list of list with values as in the confusion matrix.
    """
    pred = trained_model.predict(X)

    if only_idx is None:
        p_act_win = np.sum(y)
        n = np.sum((y == 0))
        if step_tol is None:
            p_pred = np.sum(pred)
            tp = np.sum((pred == True) & (y == 1))
            ep = tp
            lp = 0
            fn = p_act_win - tp
            fp = np.sum((pred == True) & (y == 0))
        else:
            y_idx = np.argwhere(y)
            y_idx = sorted(y_idx, key=lambda x: x[1])
            tp = 0
            fn = 0

            got_last = False
            # Go through all segment starts
            for i, y_i in enumerate(y_idx):
                if i > 0:
                    if y_i[1] == y_idx[i - 1][1]:
                        # same input parameter
                        if got_last and y_i[0] - y_idx[i - 1][0] == 1:
                            # same segment start, just another index
                            # in the other window
                            continue
                        elif not got_last and y_i[0] - y_idx[i - 1][0] > 1:
                            # new segment start and the last was not detected
                            fn += 1
                            if pred[y_i[0], y_i[1]]:
                                tp += 1
                                got_last = True
                        elif pred[y_i[0], y_i[1]]:
                            tp += 1
                            got_last = True
                        elif y_i[0] - y_idx[i - 1][0] > 1:
                            got_last = False

                    else:
                        # New input parameter
                        if not got_last:
                            fn += 1
                        got_last = False
                        if pred[y_i[0], y_i[1]]:
                            tp += 1
                            got_last = True

                else:
                    if pred[y_i[0], y_i[1]]:
                        tp += 1
                        got_last = True
            if not got_last:
                fn += 1
            # This would be fp_windows
            # fp = np.sum(pred) - not_fp - tp
            idx = np.argwhere(pred)
            idx = np.asarray(sorted(idx, key=lambda x: x[1]))
            p_pred = 0
            # index offset from start of segment to current index
            got_last = 0
            for i, pred_i in enumerate(idx):
                if i > 0:
                    if pred_i[1] == idx[i - 1][1]:
                        # same input parameter
                        if pred_i[0] - idx[i - got_last][0] <= 2 * step_tol:
                            # within tolerance => same segment
                            got_last += 1
                        else:
                            p_pred += 1
                            got_last = 1
                    else:
                        # new input parameter and therefore segment
                        got_last = 1
                        p_pred += 1
                else:
                    got_last = 1
                    p_pred += 1

            fp = p_pred - tp
            ep = 0
            lp = 0
    else:
        p_act_win = np.sum(y[:, only_idx])
        n = np.sum((y[:, only_idx] == 0))

        if step_tol == None:
            p_pred = np.sum(pred[:, only_idx])
            fp = np.sum((pred[:, only_idx] == True) & (y[:, only_idx] == 0))
            tp = np.sum((pred[:, only_idx] == True) & (y[:, only_idx] == 1))
            ep = tp
            lp = 0
            fn = p_act_win - tp
        else:
            y_idx = np.argwhere(y[:, only_idx])

            ep = 0
            lp = 0
            tp = 0
            fn = 0
            got_last = False
            # Go through all segment starts
            for i, y_i in enumerate(y_idx):
                if i > 0:
                    if got_last and y_i[0] - y_idx[i - 1][0] == 1:
                        # same segment start, just another index
                        # in the other window
                        continue
                    elif not got_last and y_i[0] - y_idx[i - 1][0] > 1:
                        # new segment start and the last was not detected
                        fn += 1
                        if pred[y_i[0], only_idx]:
                            tp += 1
                            got_last = True
                    elif pred[y_i[0], only_idx]:
                        tp += 1
                        got_last = True
                    elif y_i[0] - y_idx[i - 1][0] > 1:
                        got_last = False

                else:
                    if pred[y_i[0], only_idx]:
                        tp += 1
                        got_last = True
            if not got_last:
                fn += 1
            # This would be fp_windows
            # fp = np.sum(pred[:, only_idx]) - not_fp - tp
            idx = np.argwhere(pred[:, only_idx])
            p_pred = 0
            got_last = 0
            for i, pred_i in enumerate(idx):
                if i > 0:
                    if pred_i[0] - idx[i - got_last][0] <= 2 * step_tol:
                        # within tolerance => same segment
                        got_last += 1
                    else:
                        p_pred += 1
                        got_last = 1
                else:
                    got_last = 1
                    p_pred += 1

            fp = p_pred - tp

    if step_tol is None:
        if p_pred == 0:
            if tp == 0:
                pr = 1
            else:
                pr = 0
        else:
            pr = tp / p_pred
    else:
        if tp + fp == 0:
            pr = 0
        else:
            pr = tp / (fp + tp)
    if step_tol is None:
        if p_act_win == 0:
            if tp == 0:
                re = 1
            else:
                re = 0
        else:
            re = tp / p_act_win
    else:
        if tp + fn == 0:
            re = 0
        else:
            re = tp / (tp + fn)

    tn = n - fp
    if pr + re == 0:
        f1 = 0
    else:
        f1 = 2 * (pr * re) / (pr + re)
    fpr = fp / n
    confusion_matrix = [
        tp,
        fn,
        fp,
        tn,
        ep,  # early positive - got a segment although early (or exact)
        lp,  # late positive - positive within a tolerance after the fact
        tp + fn,
        p_act_win,
        pr,
        re,
        f1,
        fpr,
        p_pred,
    ]
    return confusion_matrix


def plot_tree_matrix(
    model,
    train,
    test,
    train_labels,
    test_labels,
    feature_names,
    class_names,
    step_tol=8,
):
    """
    TODO: Plot the single row confusion matrix. Currently this method only
    returns a single row confusion matrix.
    """
    confusion_matrix = get_tree_matrix(model, test, test_labels, step_tol=None)

    # TODO Plot results?
    return confusion_matrix


def show_tree_stats(
    model,
    train,
    test,
    train_labels,
    test_labels,
    feature_names,
    class_names,
    step_tol=8,
):
    """
    TODO: Plot the single row confusion matrix. Currently this method
    only returns a single row confusion matrix and prints some info about
    the model used.
    """
    n_nodes = []
    max_depth = []
    for i in model.estimators_:
        n_nodes.append(i.tree_.node_count)
        max_depth.append(i.tree_.max_depth)
    print(f"Mean number of nodes: {np.mean(n_nodes):1.2e}")
    print(f"Mean max depth: {np.mean(max_depth):1.2e}")

    confusion_matrix = get_tree_matrix(model, test, test_labels, step_tol=None)

    # TODO Plot results?
    return confusion_matrix


def train_many_models(
    data,
    max_features_list,
    min_threshold,
    max_threshold,
    threshold_step,
    precalc=True,
    step_tol=2,
    distinct_outparams=False,
    n_estimators=100,
    save_memory=False,
    no_split=True,
    verbosity=0,
    max_depth=36,
    max_leaf_nodes=1000,
    classifier="random_forest",
    learning_rate=1.0,
    cooldown=0,
):
    """
    Train different models for different max feature splits and thresholds
    for true segment starts.

    Parameters
    ----------
    data : list or np.ndarray of paths or dataset
        Dataset to find the segments for. If list or np.ndarray is given,
        load a pickled numpy array as testset from disk. Must have "test"
        in its name.
    max_features_list : list of string or int or float
        A list of possible number of features for the best split
        to train a model for.
        From sklearn.ensemble.RandomForestClassifier:
        The number of features to consider when looking for the best split:

        If int, then consider max_features features at each split.
        If float, then max_features is a fraction and
            round(max_features * n_features) features are considered at each split.
        If “auto”, then max_features=sqrt(n_features).
        If “sqrt”, then max_features=sqrt(n_features) (same as “auto”).
        If “log2”, then max_features=log2(n_features).
        If None, then max_features=n_features.

        Note: the search for a split does not stop until at least one valid
        partition of the node samples is found, even if it requires to
        effectively inspect more than max_features features.
    min_threshold : float
        Minimum threshold for errors to identify true segment starts.
    max_threshold : float
        Maximum threshold for errors to identify true segment starts.
    threshold_step : float
        Step size between minimum and maximum threshold.
    precalc : bool
        If true, calculate segment starts for different thresholds
        before everything else, otherwise calculate those a new
        for every feature.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    distinct_outparams : bool
        If true, try to predict perturbing an input parameter for a segment
        start for each output parameter independently. This may not be useful,
        since one input parameter can have a high impact on multiple output
        parameters.
    n_estimators : int
        Number of trees.
    save_memory : bool
        If true and not "no_split", use the first "test_size" many rows for
        training. This means, that the data is not stratified! Do not use this
        unless the classes to predict is roughly uniformly distributed.
    no_split : bool
        If true, do not split data into training and test set and use everything
        for training.
    verbosity : int
        Set verbosity level.
    max_depth : int
        From sklearn.ensemble.RandomForestClassifier:
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split (here: 2) samples.
    max_leaf_nodes : int
        From sklearn.ensemble.RandomForestClassifier:
        Grow trees with max_leaf_nodes in best-first fashion. Best nodes are
        defined as relative reduction in impurity. If None then unlimited
        number of leaf nodes.
    classifier : string
        Choose a classifier for prediction. Options are "random_forest"
        and "adaboost". The latter ignores "max_leaf_nodes" and "max_depth" and
        rather uses "learning_rate".
    learning_rate : float
        From sklean.ensemble.AdaBoostClassifier:
        Learning rate shrinks the contribution of each classifier by
        "learning_rate". There is a trade-off between "learning_rate" and
        "n_estimators".
    cooldown : int
        Minimum number of timesteps between last time the error threshold has
        been met and the next time step. This is useful if the original data
        is based on ensembles that start at every few time step which leads to
        an error of zero until the divergence of the ensemble is high enough
        again, resulting in a new segment start although it is just the
        start of the ensemble.

    Returns
    -------
    List of list of trained models, where the first dimension corresponds
    to the max feature split used and the second dimension corresponds
    to the different segment thresholds.
    Also returns np.array of segment thresholds.
    """
    seg_thresholds = np.arange(min_threshold, max_threshold, threshold_step)

    # Pre calculate the segment positions?
    if not isinstance(data, list) and not isinstance(data, np.ndarray):
        load_data = False
        ds_cache = {}
        if precalc:
            for seg_thresh in seg_thresholds:
                ds_cache[seg_thresh] = find_segments(data, 10.0 ** seg_thresh, cooldown)
    else:
        load_data = True

    models = []

    for feat_idx, max_features in enumerate(max_features_list):
        models_inner = []
        for seg_idx, seg_thresh in enumerate(seg_thresholds):
            if load_data:
                ds = []
                for path in data:
                    if "test" in path:
                        continue
                    if str(seg_thresh) in path or f"{seg_thresh:.1e}" in path:
                        ds.append(np.load(path, allow_pickle=True, fix_imports=False))
                    if len(ds) == 2:
                        ds[0] = ds[0].astype(int, copy=False)
                        break
            elif precalc:
                ds = ds_cache[seg_thresh]
            else:
                ds = find_segments(data, 10.0 ** seg_thresh, cooldown)
            if verbosity > 2:
                print(f"Training for {feat_idx}, {seg_idx}")

            model, _, _, _, _, _, _ = create_forest(
                ds,
                step_tol=step_tol,
                distinct_outparams=distinct_outparams,
                n_estimators=n_estimators,
                max_features=max_features,
                save_memory=save_memory,
                no_split=no_split,
                verbosity=verbosity > 3,
                max_depth=max_depth,
                max_leaf_nodes=max_leaf_nodes,
                classifier=classifier,
                learning_rate=learning_rate,
            )
            models_inner.append(model)
        models.append(models_inner)
    return models, seg_thresholds


def create_dataset_pretrained(
    data,
    models,
    seg_thresholds,
    max_features_list,
    step_tol,
    distinct_outparams,
    precalc=True,
    verbosity=0,
    in_params=None,
    predict_train=False,
    independent_windows=True,
    cooldown=0,
):
    """
    Create a dataset with dimensions "Max Features", "Output Parameter",
    "Input Parameter", "Segment Threshold" and columns "TP over P",
    "FP over P", "FN over P_windows", "TP", "P", "FP", "FN", "P_window" by
    predicting segment starts using given trained models.

    Parameters
    ----------
    data : list or np.ndarray of paths or dataset
        Dataset to find the segments for. If list or np.ndarray is given,
        load a pickled numpy array as testset from disk. Must have "test"
        in its name.
    models : List of list of trained models
        These are models trained by train_many_models(..)
    seg_thresholds : list or np.array of float
        Array of segment thresholds for true segment starts.
    max_features_list : list of string or int or float
        A list of possible number of features for the best split
        to train a model for.
        From sklearn.ensemble.RandomForestClassifier:
        The number of features to consider when looking for the best split:

        If int, then consider max_features features at each split.
        If float, then max_features is a fraction and
            round(max_features * n_features) features are considered at each split.
        If “auto”, then max_features=sqrt(n_features).
        If “sqrt”, then max_features=sqrt(n_features) (same as “auto”).
        If “log2”, then max_features=log2(n_features).
        If None, then max_features=n_features.

        Note: the search for a split does not stop until at least one valid
        partition of the node samples is found, even if it requires to
        effectively inspect more than max_features features.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    distinct_outparams : bool
        If true, try to predict perturbing an input parameter for a segment
        start for each output parameter independently. This may not be useful,
        since one input parameter can have a high impact on multiple output
        parameters.
    precalc : bool
        If true, calculate segment starts for different thresholds
        before everything else, otherwise calculate those a new
        for every feature.
    verbosity : int
        Set verbosity level.
    in_params : list
        List of input parameters to create a separate prediction for.
    predict_train : bool
        Only useful if "data" is a list or array of paths to load data from.
        If true, use train data for the confusion matrix, otherwise use
        test data.
    cooldown : int
        Minimum number of timesteps between last time the error threshold has
        been met and the next time step. This is useful if the original data
        is based on ensembles that start at every few time step which leads to
        an error of zero until the divergence of the ensemble is high enough
        again, resulting in a new segment start although it is just the
        start of the ensemble.
    independent_windows : bool or None
        If true: count true positive predictions for each window independently.
        If false: count detected segments within tolerance as true positive
        predictions, ignoring multiple true labels for the same segment
        and multiple true predictions.

    Returns
    -------
    xarray.Dataset of confusion matrix as described above.
    """
    if independent_windows:
        pass_step_tol = None
    else:
        pass_step_tol = step_tol

    if isinstance(data, list) or isinstance(data, np.ndarray):
        precalc = False
        load_data = True
    else:
        load_data = False

    dims_forest = (len(max_features_list), 1, 1 + len(in_params), len(seg_thresholds))
    if predict_train:
        input_param_coords = ["All Input Parameters Train Set"] + in_params
    else:
        input_param_coords = ["All Input Parameters Test Set"] + in_params
    coords_forest_dic = {
        "Max Features": max_features_list,
        "Output Parameter": ["All Output Parameters"],
        "Input Parameter": input_param_coords,
        "Segment Threshold": seg_thresholds,
    }

    tpp_forest = np.zeros(dims_forest)
    fpp_forest = np.zeros(dims_forest)
    fnpw_forest = np.zeros(dims_forest)
    tp_forest = np.zeros(dims_forest, dtype=int)
    p_forest = np.zeros(dims_forest, dtype=int)
    fp_forest = np.zeros(dims_forest, dtype=int)
    fn_forest = np.zeros(dims_forest, dtype=int)
    p_w_forest = np.zeros(dims_forest, dtype=int)

    # Pre calculate the segment positions?
    test_cache = {}
    if precalc:
        for seg_thresh in seg_thresholds:
            ds = find_segments(data, 10.0 ** seg_thresh, cooldown)
            test, test_labels, _, _ = create_input_labels(
                ds, step_tol, distinct_outparams
            )
            test_cache[seg_thresh] = (test, test_labels)

    for seg_idx, seg_thresh in enumerate(seg_thresholds):
        if load_data:
            test = None
            test_labels = None
            for path in data:
                if not "test" in path and not predict_train:
                    continue
                if "test" in path and predict_train:
                    continue
                if not (str(seg_thresh) in path or f"{seg_thresh:.1e}" in path):
                    continue
                if "labels" in path:
                    test_labels = np.load(path, allow_pickle=True, fix_imports=False)
                    test_labels = test_labels.astype(int, copy=False)
                else:
                    test = np.load(path, allow_pickle=True, fix_imports=False)
                if test is not None and test_labels is not None:
                    break
        elif precalc:
            test, test_labels = test_cache[seg_thresh]
        else:
            ds = find_segments(data, 10.0 ** seg_thresh, cooldown)
            test, test_labels, _, _ = create_input_labels(
                ds, step_tol, distinct_outparams
            )
        for feat_idx, max_features in enumerate(max_features_list):

            if verbosity > 2:
                print(
                    f"Running for {feat_idx}, {seg_idx} ({max_features}, {seg_thresh})"
                )
                t1 = timer()
            model = models[feat_idx][seg_idx]
            matrix = get_tree_matrix(model, test, test_labels, step_tol=pass_step_tol)
            # For paranoia reasons
            np.nan_to_num(matrix, copy=False)

            tpp_forest[feat_idx, 0, 0, seg_idx] = matrix[9]
            if matrix[6] != 0:
                fpp_forest[feat_idx, 0, 0, seg_idx] = matrix[2] / matrix[6]
            if matrix[7] != 0:
                fnpw_forest[feat_idx, 0, 0, seg_idx] = matrix[1] / matrix[7]
            tp_forest[feat_idx, 0, 0, seg_idx] = matrix[0]
            p_forest[feat_idx, 0, 0, seg_idx] = matrix[6]
            fp_forest[feat_idx, 0, 0, seg_idx] = matrix[2]
            fn_forest[feat_idx, 0, 0, seg_idx] = matrix[1]
            p_w_forest[feat_idx, 0, 0, seg_idx] = matrix[7]
            if verbosity > 4:
                print(
                    f"\nTP: {matrix[0]}; FP: {matrix[2]}; FN: {matrix[1]}; P_w: {matrix[7]}; P: {matrix[6]}"
                )

            for in_idx, in_p in enumerate(in_params):
                matrix = get_tree_matrix(
                    model, test, test_labels, only_idx=in_idx, step_tol=pass_step_tol
                )
                np.nan_to_num(matrix, copy=False)

                tpp_forest[feat_idx, 0, 1 + in_idx, seg_idx] = matrix[9]
                if matrix[6] != 0:
                    fpp_forest[feat_idx, 0, 1 + in_idx, seg_idx] = matrix[2] / matrix[6]
                if matrix[7] != 0:
                    fnpw_forest[feat_idx, 0, 1 + in_idx, seg_idx] = (
                        matrix[1] / matrix[7]
                    )
                tp_forest[feat_idx, 0, 1 + in_idx, seg_idx] = matrix[0]
                p_forest[feat_idx, 0, 1 + in_idx, seg_idx] = matrix[6]
                fp_forest[feat_idx, 0, 1 + in_idx, seg_idx] = matrix[2]
                fn_forest[feat_idx, 0, 1 + in_idx, seg_idx] = matrix[1]
                p_w_forest[feat_idx, 0, 1 + in_idx, seg_idx] = matrix[7]
                if verbosity > 4:
                    print(f"{in_idx} ({in_p})")
                    print(
                        f"TP: {matrix[0]}; FP: {matrix[2]}; FN: {matrix[1]}; P_w: {matrix[7]}; P: {matrix[6]}"
                    )
            if verbosity > 2:
                t2 = timer()
                print(f" ... done in {t2-t1} s")

    dim_names = [
        "Max Features",
        "Output Parameter",
        "Input Parameter",
        "Segment Threshold",
    ]

    return xr.Dataset(
        data_vars={
            "TP over P": (dim_names, tpp_forest),
            "FP over P": (dim_names, fpp_forest),
            "FN over P_windows": (dim_names, fnpw_forest),
            "TP": (dim_names, tp_forest),
            "P": (dim_names, p_forest),
            "FP": (dim_names, fp_forest),
            "FN": (dim_names, fn_forest),
            "P_windows": (dim_names, p_w_forest),
        },
        coords=coords_forest_dic,
    )


def create_dataset_forest(
    data,
    max_features_list,
    min_threshold,
    max_threshold,
    threshold_step,
    precalc=True,
    step_tol=2,
    test_size=0.25,
    distinct_outparams=False,
    n_estimators=100,
    save_memory=False,
    no_split=False,
    verbosity=0,
    max_depth=36,
    max_leaf_nodes=1000,
    classifier="random_forest",
    learning_rate=1.0,
    cooldown=0,
):
    """
    Create a dataset with dimensions "Max Features", "Output Parameter",
    "Input Parameter", "Segment Threshold" and columns "TP over P",
    "FP over P", "FN over P_windows", "TP", "P", "FP", "FN", "P_window" by
    predicting segment starts. Models are being trained used create_forest(..)
    in this method as well.

    Parameters
    ----------
    data : dataset
        Dataframe with MSE created by load_dataset() for every model state
        parameter and all trajectories to find segments for.
    max_features_list : list of string or int or float
        A list of possible number of features for the best split
        to train a model for.
        From sklearn.ensemble.RandomForestClassifier:
        The number of features to consider when looking for the best split:

        If int, then consider max_features features at each split.
        If float, then max_features is a fraction and
            round(max_features * n_features) features are considered at each split.
        If “auto”, then max_features=sqrt(n_features).
        If “sqrt”, then max_features=sqrt(n_features) (same as “auto”).
        If “log2”, then max_features=log2(n_features).
        If None, then max_features=n_features.

        Note: the search for a split does not stop until at least one valid
        partition of the node samples is found, even if it requires to
        effectively inspect more than max_features features.
    min_threshold : float
        Minimum threshold for errors to identify true segment starts.
    max_threshold : float
        Maximum threshold for errors to identify true segment starts.
    threshold_step : float
        Step size between minimum and maximum threshold.
    precalc : bool
        If true, calculate segment starts for different thresholds
        before everything else, otherwise calculate those a new
        for every feature.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    test_size : float
        Percentage of number of windows to use for testing.
    distinct_outparams : bool
        If true, try to predict perturbing an input parameter for a segment
        start for each output parameter independently. This may not be useful,
        since one input parameter can have a high impact on multiple output
        parameters.
    n_estimators : int
        Number of trees.
    save_memory : bool
        If true and not "no_split", use the first "test_size" many rows for
        training. This means, that the data is not stratified! Do not use this
        unless the classes to predict is roughly uniformly distributed.
    no_split : bool
        If true, do not split data into training and test set and use everything
        for training.
    verbosity : int
        Set verbosity level.
    max_depth : int
        From sklearn.ensemble.RandomForestClassifier:
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split (here: 2) samples.
    max_leaf_nodes : int
        From sklearn.ensemble.RandomForestClassifier:
        Grow trees with max_leaf_nodes in best-first fashion. Best nodes are
        defined as relative reduction in impurity. If None then unlimited
        number of leaf nodes.
    classifier : string
        Choose a classifier for prediction. Options are "random_forest"
        and "adaboost". The latter ignores "max_leaf_nodes" and "max_depth" and
        rather uses "learning_rate".
    learning_rate : float
        From sklean.ensemble.AdaBoostClassifier:
        Learning rate shrinks the contribution of each classifier by
        "learning_rate". There is a trade-off between "learning_rate" and
        "n_estimators".
    cooldown : int
        Minimum number of timesteps between last time the error threshold has
        been met and the next time step. This is useful if the original data
        is based on ensembles that start at every few time step which leads to
        an error of zero until the divergence of the ensemble is high enough
        again, resulting in a new segment start although it is just the
        start of the ensemble.

    Returns
    -------
    xarray.Dataset of confusion matrix as described above.
    """
    # We create a dataset of different thresholds (good for ROC and so on)
    # Fraction for testing and numbers of estimators don't really matter, I guess
    # Test and trainingssets go into different dataset since we cannot differentiate
    # between different output parameters, just input parameters
    seg_thresholds = np.arange(min_threshold, max_threshold, threshold_step)
    in_params = list(np.unique(data["Input Parameter"]))
    dims_forest = (len(max_features_list), 1, 2 + len(in_params), len(seg_thresholds))
    coords_forest_dic = {
        "Max Features": max_features_list,
        "Output Parameter": ["All Output Parameters"],
        "Input Parameter": [
            "All Input Parameters Trainingset",
            "All Input Parameters Testset",
        ]
        + in_params,
        "Segment Threshold": seg_thresholds,
    }

    tpp_forest = np.zeros(dims_forest)
    fpp_forest = np.zeros(dims_forest)
    fnpw_forest = np.zeros(dims_forest)
    tp_forest = np.zeros(dims_forest, dtype=int)
    p_forest = np.zeros(dims_forest, dtype=int)
    fp_forest = np.zeros(dims_forest, dtype=int)
    fn_forest = np.zeros(dims_forest, dtype=int)
    p_w_forest = np.zeros(dims_forest, dtype=int)

    # Pre calculate the segment positions?
    ds_cache = {}
    if precalc:
        for seg_thresh in seg_thresholds:
            ds_cache[seg_thresh] = find_segments(data, 10.0 ** seg_thresh, cooldown)

    for feat_idx, max_features in enumerate(max_features_list):
        for seg_idx, seg_thresh in enumerate(seg_thresholds):
            if precalc:
                ds = ds_cache[seg_thresh]
            else:
                ds = find_segments(data, 10.0 ** seg_thresh, cooldown)
            if verbosity > 2:
                print(f"Running for {feat_idx}, {seg_idx}")

            model, train, test, train_labels, test_labels, _, _ = create_forest(
                ds,
                step_tol=step_tol,
                test_size=test_size,
                distinct_outparams=distinct_outparams,
                n_estimators=n_estimators,
                max_features=max_features,
                save_memory=save_memory,  # Old version set to False
                no_split=no_split,
                verbosity=verbosity > 3,
                max_depth=max_depth,
                max_leaf_nodes=max_leaf_nodes,
                classifier=classifier,
                learning_rate=learning_rate,
            )
            if verbosity > 2:
                print("Trained model")
            test_matrix = get_tree_matrix(model, test, test_labels, step_tol=None)
            train_matrix = get_tree_matrix(model, train, train_labels, step_tol=None)
            if verbosity > 2:
                print("Got confuse matrix")
            tpp_forest[feat_idx, 0, 0, seg_idx] = train_matrix[9]
            fpp_forest[feat_idx, 0, 0, seg_idx] = train_matrix[2] / train_matrix[6]
            fnpw_forest[feat_idx, 0, 0, seg_idx] = train_matrix[1] / train_matrix[7]
            tp_forest[feat_idx, 0, 0, seg_idx] = train_matrix[0]
            p_forest[feat_idx, 0, 0, seg_idx] = train_matrix[6]
            fp_forest[feat_idx, 0, 0, seg_idx] = train_matrix[2]
            fn_forest[feat_idx, 0, 0, seg_idx] = train_matrix[1]
            p_w_forest[feat_idx, 0, 0, seg_idx] = train_matrix[7]

            tpp_forest[feat_idx, 0, 1, seg_idx] = test_matrix[9]
            fpp_forest[feat_idx, 0, 1, seg_idx] = test_matrix[2] / test_matrix[6]
            fnpw_forest[feat_idx, 0, 1, seg_idx] = test_matrix[1] / test_matrix[7]
            tp_forest[feat_idx, 0, 1, seg_idx] = test_matrix[0]
            p_forest[feat_idx, 0, 1, seg_idx] = test_matrix[6]
            fp_forest[feat_idx, 0, 1, seg_idx] = test_matrix[2]
            fn_forest[feat_idx, 0, 1, seg_idx] = test_matrix[1]
            p_w_forest[feat_idx, 0, 1, seg_idx] = test_matrix[7]

            for in_idx, in_p in enumerate(in_params):
                matrix = get_tree_matrix(
                    model, test, test_labels, only_idx=in_idx, step_tol=None
                )

                tpp_forest[feat_idx, 0, 2 + in_idx, seg_idx] = matrix[9]
                fpp_forest[feat_idx, 0, 2 + in_idx, seg_idx] = matrix[2] / matrix[6]
                fnpw_forest[feat_idx, 0, 2 + in_idx, seg_idx] = matrix[1] / matrix[7]
                tp_forest[feat_idx, 0, 2 + in_idx, seg_idx] = matrix[0]
                p_forest[feat_idx, 0, 2 + in_idx, seg_idx] = matrix[6]
                fp_forest[feat_idx, 0, 2 + in_idx, seg_idx] = matrix[2]
                fn_forest[feat_idx, 0, 2 + in_idx, seg_idx] = matrix[1]
                p_w_forest[feat_idx, 0, 2 + in_idx, seg_idx] = matrix[7]

    dim_names = [
        "Max Features",
        "Output Parameter",
        "Input Parameter",
        "Segment Threshold",
    ]
    return xr.Dataset(
        data_vars={
            "TP over P": (dim_names, tpp_forest),
            "FP over P": (dim_names, fpp_forest),
            "FN over P_windows": (dim_names, fnpw_forest),
            "TP": (dim_names, tp_forest),
            "P": (dim_names, p_forest),
            "FP": (dim_names, fp_forest),
            "FN": (dim_names, fn_forest),
            "P_windows": (dim_names, p_w_forest),
        },
        coords=coords_forest_dic,
    )


def get_stratified_sets(
    ds, step_tol, distinct_outparams=False, train_size=0.75, verbosity=0
):
    """
    Create a stratified train set from the given dataset.
    The idea here: Call this function multiple times for subsets
    of the df and concatenate the results. This way, we can create
    a bigger training set without blowing up the needed RAM.

    Parameters
    ----------
    ds : Dataset
        Dataset created by find_segments(..) where segments are identified and
        predicted errors are stored.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    distinct_outparams : bool
        If true, try to predict perturbing an input parameter for a segment
        start for each output parameter independently. This may not be useful,
        since one input parameter can have a high impact on multiple output
        parameters.
    train_size : float
        Size in percentage for the training set.
    verbosity : int
        Set verbosity level.

    Returns
    -------
    Four np.arrays with train set, train labels, test set and test labels.
    """

    X_final, y_final, _, _ = create_input_labels(
        ds, step_tol, distinct_outparams, verbosity
    )

    train, train_labels, test, test_labels = iterative_train_test_split(
        X_final, y_final, test_size=1 - train_size
    )
    return train, train_labels, test, test_labels


def create_big_stratified_set(
    data_path,
    step_tol,
    all_params_list,
    n_trajs_iter,
    out_params,
    threshold,
    distinct_outparams=False,
    train_size=0.75,
    cooldown=0,
    min_time=-1000,
    verbosity=0,
):
    """
    Load part of one file or multiple files iteratively and create
    stratified subsets for training and concatenate them. Store the result
    on disk for training models on this.

    Parameters
    ----------
    data_path : path
        If string ends with ".nc", load a single file.
        Otherwise load and append multiple files to test and training sets.
        Only processes "n_trajs_iter" many trajectories per step.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    all_params_list : list of str
        List of all input params to get predicted errors for.
    n_trajs_iter : int
        Number of trajectories to process per iteration. Using 20
        trajectories may consume around 24 GB of RAM.
    out_params : list of string
        List of output parameters.
    threshold : float
        Exponent with base 10 for theshold of true segment start.
    distinct_outparams : bool
        If true, try to predict perturbing an input parameter for a segment
        start for each output parameter independently. This may not be useful,
        since one input parameter can have a high impact on multiple output
        parameters.
    train_size : float
        Size in percentage for the training set.
    cooldown : int
        Minimum number of timesteps between last time the error threshold has
        been met and the next time step. This is useful if the original data
        is based on ensembles that start at every few time step which leads to
        an error of zero until the divergence of the ensemble is high enough
        again, resulting in a new segment start although it is just the
        start of the ensemble.
    min_time : float or int
        Minimum value of "time_after_ascent". This is useful if ensembles
        started later than the baseline trajectory. If None is given, use
        all time steps.
    verbosity : int
        Set verbosity level.

    Returns
    -------
    Four np.arrays with train set, train labels, test set and test labels.
    """
    all_train = None
    all_labels = None
    all_test = None
    all_test_labels = None

    def append_data(ds):
        nonlocal all_train
        nonlocal all_labels
        nonlocal all_test
        nonlocal all_test_labels

        if min_time is not None:
            ds = ds.where(ds["time_after_ascent"] >= min_time, drop=True)

        total_trajs = len(ds.trajectory)
        if total_trajs < n_trajs_iter:
            traj_step = total_trajs
        else:
            traj_step = n_trajs_iter
        n_steps = int((total_trajs + traj_step - 1) / traj_step)
        for i in range(n_steps):
            if verbosity > 1:
                print(f"{i}: {i*traj_step} - {(i+1)*traj_step}, {total_trajs}")

            if (i + 1) * traj_step >= total_trajs:
                ds_tmp = find_segments(
                    ds.sel(
                        {"trajectory": ds.trajectory[i * traj_step : total_trajs - 1]}
                    ),
                    10.0 ** threshold,
                    cooldown,
                )
            else:
                ds_tmp = find_segments(
                    ds.sel(
                        {
                            "trajectory": ds.trajectory[
                                i * traj_step : (i + 1) * traj_step
                            ]
                        }
                    ),
                    10.0 ** threshold,
                    cooldown,
                )

            train, train_labels, test, test_labels = get_stratified_sets(
                ds=ds_tmp,
                step_tol=step_tol,
                distinct_outparams=distinct_outparams,
                train_size=train_size,
                verbosity=verbosity,
            )
            if all_train is not None:
                all_train = np.append(all_train, train, axis=0)
                all_labels = np.append(all_labels, train_labels, axis=0)
                all_test = np.append(all_test, test, axis=0)
                all_test_labels = np.append(all_test_labels, test_labels, axis=0)
            else:
                all_train = train
                all_labels = train_labels
                all_test = test
                all_test_labels = test_labels

    if ".nc" in data_path:
        # ie data2_327.nc
        if verbosity > 1:
            print(f"Loading {data_path}")
        data = xr.open_dataset(data_path, decode_times=False, engine="netcdf4")
        append_data(data)

    else:
        if verbosity > 1:
            print(f"Checking {data_path}")
        paths = list(os.listdir(data_path))
        if verbosity > 1:
            print(f"Loading from {paths}")
        traj_offset = 0
        for p in range(len(paths)):
            if "quan" in paths[p] or "median" in paths[p] or "perturb" in paths[p]:
                continue
            path = data_path + paths[p] + "/"
            n_trajs = len(list(os.listdir(path)))
            if verbosity > 1:
                print(f"Loading from {path} with {n_trajs} trajectories")
            try:
                data = load_dataset(
                    path=path,
                    out_params=out_params,
                    in_params=all_params_list,
                    traj_list=np.arange(n_trajs),
                    traj_offset=traj_offset,
                    verbosity=verbosity,
                )
            except:
                print(f"Loading      {path} failed.")
                continue
            append_data(data)
            traj_offset += n_trajs
    return all_train, all_labels, all_test, all_test_labels


if __name__ == "__main__":
    import argparse
    from joblib import dump, load
    import os
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Using perturbed ensembles: Create random forests to predict
            segments (start of deviations of ensembles from unperturbed trajectory)
            based on predicted errors from algorithmic differentiation (AD).
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--load_models",
        dest="model_path",
        default=None,
        help=textwrap.dedent(
            """\
            Path to random forest models that had been created earlier.
            They are named as rand_forest_{max_features}_{seg_thresh:.1e}.joblid.
            """
        ),
    )
    parser.add_argument(
        "--store_models",
        default=None,
        help=textwrap.dedent(
            """\
            If a path is given, store models on disk
            named as rand_forest_{max_features}_{seg_thresh:.1e}.joblid.
            """
        ),
    )
    parser.add_argument(
        "--train_subset",
        action="store_true",
        help=textwrap.dedent(
            """\
            Train the models on a small subset of the data. This may be
            needed since stratifying the complete dataset can consume too much memory.
            Using a random sample of the complete dataset leads to poor performance
            since segment starts are rare and negative examples might be too
            dominant without startification.
            If 'no_split' is used, no stratification will be done and all
            data defined by 'n_trajs' is used for training. Otherwise only
            a quarter of 'n_trajs' is used for training but with stratification.
            """
        ),
    )
    parser.add_argument(
        "--data_path",
        default="/data/project/wcb/netcdf/perturbed_ensembles/",
        help=textwrap.dedent(
            """\
            Path to folders with ensemble datasets or to single NetCDF file
            with all data concatenated along 'trajectory' axis.
            If a path to numpy arrays is given, it is assumed to be a training
            or test set.
            """
        ),
    )
    parser.add_argument(
        "--min_threshold",
        type=float,
        default=-40,
        help=textwrap.dedent(
            """\
            Minimum threshold (power of 10) for classifying a time step in a
            perturbed ensemble as a segment start.
            """
        ),
    )
    parser.add_argument(
        "--max_threshold",
        type=float,
        default=0,
        help=textwrap.dedent(
            """\
            Maximum threshold (power of 10) for classifying a time step in a
            perturbed ensemble as a segment start.
            """
        ),
    )
    parser.add_argument(
        "--threshold_step",
        type=float,
        default=2,
        help=textwrap.dedent(
            """\
            Step size from minimum threshold to maximum threshold for
            classifying a time step in a perturbed ensemble as a segment start.
            """
        ),
    )
    parser.add_argument(
        "--step_tol",
        type=int,
        default=8,
        help=textwrap.dedent(
            """\
            Tolerance for predicting a segment start in time steps. The
            size of the resulting sliding window is step_tol*2+1.
            """
        ),
    )
    parser.add_argument(
        "--n_estimators",
        type=int,
        default=100,
        help=textwrap.dedent(
            """\
            Number of estimators for the random forest.
            """
        ),
    )
    parser.add_argument(
        "--n_trajs",
        type=int,
        default=10,
        help=textwrap.dedent(
            """\
            Number of trajectories that will be predicted using the trained models
            at once. If 'train_subset' is set, then this is the number of
            trajectories used for training the models.
            """
        ),
    )
    parser.add_argument(
        "--store_name",
        default="prediction.nc",
        help=textwrap.dedent(
            """\
            Name of NetCDF file where results (confusion matrix) will be stored.
            Must end with '.nc'!
            """
        ),
    )
    parser.add_argument(
        "--save_memory",
        action="store_true",
        help=textwrap.dedent(
            """\
            Save memory by splitting the dataset into a train and test set by
            using the first trajectories as training set without stratification.
            This can lead to poor performance of the models.
            """
        ),
    )
    parser.add_argument(
        "--no_split",
        action="store_true",
        help=textwrap.dedent(
            """\
            Do not split the dataset into a training and test set.
            If 'train_subset' is used, then the subset defined by 'n_trajs' is
            used for training. Otherwise the complete dataset is used for training.
            Overrides 'save_memory'.
            """
        ),
    )
    parser.add_argument(
        "--store_many_appended_data",
        default=None,
        help=textwrap.dedent(
            """\
            Store the appended input data to this path as NetCDF file for each appended
            version. Used mainly for debugging.
            """
        ),
    )
    parser.add_argument(
        "--store_appended_data",
        default=None,
        help=textwrap.dedent(
            """\
            Store the final appended data to this path and name. Must end with
            '.nc' for datasets. For training sets, '.nc' will be stripped away.
            """
        ),
    )
    parser.add_argument(
        "--only_append",
        action="store_true",
        help=textwrap.dedent(
            """\
            Only appending of data. Use 'store_appended_data' to define
            a path and name.
            """
        ),
    )
    parser.add_argument(
        "--only_training",
        action="store_true",
        help=textwrap.dedent(
            """\
            Only model training. Use 'store_models' to define a path to save the
            models.
            """
        ),
    )
    parser.add_argument(
        "--create_trainset",
        action="store_true",
        help=textwrap.dedent(
            """\
            Create a trainingset with labels and store it at 'store_appended_data'
            for training later.
            """
        ),
    )
    parser.add_argument(
        "--predict_trainset",
        action="store_true",
        help=textwrap.dedent(
            """\
            Predict the training set created by 'create_trainset'.
            """
        ),
    )
    parser.add_argument(
        "--feature_split",
        default=["log2", "sqrt", 1.0],
        nargs="+",
        help=textwrap.dedent(
            """\
            A list of possible number of features for the best split
            to train a model for.
            From sklearn.ensemble.RandomForestClassifier:
            The number of features to consider when looking for the best split:
    
            If int, then consider max_features features at each split.
            If float, then max_features is a fraction and
                round(max_features * n_features) features are considered at each split.
            If “auto”, then max_features=sqrt(n_features).
            If “sqrt”, then max_features=sqrt(n_features) (same as “auto”).
            If “log2”, then max_features=log2(n_features).
            If None, then max_features=n_features.
    
            Note: the search for a split does not stop until at least one valid
            partition of the node samples is found, even if it requires to
            effectively inspect more than max_features features.
            """
        ),
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help=textwrap.dedent(
            """\
            Set verbosity level.
            0: No output except for exceptions
            1: Print datasets
            2: Print loading statements
            3: Print building/training statements
            4: Set verbosity of random forest to 2
            5: Get predicted results for each entry
            6: Get count for segment starts for every input parameter
            """
        ),
    )
    parser.add_argument(
        "--max_depth",
        type=int,
        default=36,
        help=textwrap.dedent(
            """\
            From sklearn.ensemble.RandomForestClassifier:
            The maximum depth of the tree. If None, then nodes are expanded until
            all leaves are pure or until all leaves contain less than
            min_samples_split (here: 2) samples.
            """
        ),
    )
    parser.add_argument(
        "--max_leaf_nodes",
        type=int,
        default=1000,
        help=textwrap.dedent(
            """\
            From sklearn.ensemble.RandomForestClassifier:
            Grow trees with max_leaf_nodes in best-first fashion. Best nodes are
            defined as relative reduction in impurity. If None then unlimited
            number of leaf nodes.
            """
        ),
    )
    parser.add_argument(
        "--load_on_the_fly",
        action="store_true",
        help=textwrap.dedent(
            """\
            Load data and find the segments on the fly for predicting a dataset,
            or load precalculated training or test set.
            """
        ),
    )
    parser.add_argument(
        "--classifier",
        type=str,
        default="random_forest",
        help=textwrap.dedent(
            """\
            Choose a classifier for prediction. Options are "random_forest"
            and "adaboost". The latter ignores "max_leaf_nodes" and "max_depth" and
            rather uses "learning_rate".
            """
        ),
    )
    parser.add_argument(
        "--learning_rate",
        type=float,
        default=1.0,
        help=textwrap.dedent(
            """\
            From sklean.ensemble.AdaBoostClassifier:
            Learning rate shrinks the contribution of each classifier by
            "learning_rate". There is a trade-off between "learning_rate" and
            "n_estimators".
            """
        ),
    )
    parser.add_argument(
        "--cooldown",
        type=int,
        default=0,
        help=textwrap.dedent(
            """\
            Minimum number of timesteps between last time the error threshold has
            been met and the next time step. This is useful if the original data
            is based on ensembles that start at every few time step which leads to
            an error of zero until the divergence of the ensemble is high enough
            again, resulting in a new segment start although it is just the
            start of the ensemble.
            """
        ),
    )
    parser.add_argument(
        "--independent_windows",
        action="store_true",
        help=textwrap.dedent(
            """\
            If true: count true positive predictions for each window independently.
            This *should* be done when using stratified datasets.
            If false: count detected segments within tolerance as true positive
            predictions, ignoring multiple true labels for the same segment
            and multiple true predictions. This is arguably closer to real life
            usage of the model.
            """
        ),
    )
    parser.add_argument(
        "--in_params",
        nargs="+",
        default=[],
        help=textwrap.dedent(
            """\
            You can define a subset of input parameters that shall be used.
            This is useful if the ensemble simulation was not run with
            perturbing all parameters but only a few.
            If none are given, i.e. this parameter is not used, all ensemble
            simulations will be loaded.
            """
        ),
    )
    parser.add_argument(
        "--out_params",
        nargs="+",
        default=[
            "QC",
            "QR",
            "QV",
            "NCCLOUD",
            "NCRAIN",
            "QI",
            "NCICE",
            "QS",
            "NCSNOW",
            "QG",
            "NCGRAUPEL",
            "QH",
            "NCHAIL",
            "QI_OUT",
            "QS_OUT",
            "QR_OUT",
            "QG_OUT",
            "QH_OUT",
            "NI_OUT",
            "NS_OUT",
            "NR_OUT",
            "NG_OUT",
            "NH_OUT",
        ],
        type=str,
        help=textwrap.dedent(
            """\
            You can define a subset of output parameters
            (aka model state variables or prognostic variables) that shall be used.
            """
        ),
    )
    parser.add_argument(
        "--time_integrated",
        action="store_true",
        help=textwrap.dedent(
            """\
            If true, calculate the time integrated sensitivities and model states.
            Otherwise get the difference at each time step.
            """
        ),
    )

    args = parser.parse_args()
    if args.classifier != "random_forest" and args.classifier != "adaboost":
        print(
            f"Classifier ({args.classifier}): Classifier must be 'random_forest' or 'adaboost'."
        )

    if args.store_name is not None:
        if ".nc" not in args.store_name:
            print(f"You must add '.nc' to {args.store_name}.")
            print(f"Result will be saved as {args.store_name}.nc")
            store_name = args.store_name + ".nc"
        else:
            store_name = args.store_name
    out_params = args.out_params
    if len(out_params) == 0:
        out_params = np.arange(len(param_id_map))
    elif not out_params[0].isdigit():
        for i in range(len(out_params)):
            out_params[i] = param_id_map.index(out_params[i])
    else:
        for i in range(len(out_params)):
            out_params[i] = int(out_params[i])

    all_params_list = args.in_params
    if len(all_params_list) == 0:
        in_params = np.sort(list(in_params_notation_mapping.keys()))
        for in_p in in_params:
            if (
                "Not tracked" in in_params_notation_mapping[in_p][0]
                or "Not used" in in_params_notation_mapping[in_p][0]
                or "one-moment warm physics" in in_params_notation_mapping[in_p][0]
                or "dependent" == in_params_notation_mapping[in_p][3]
                or "formulation by Miltenberger et al."
                in in_params_notation_mapping[in_p][0]
                or "_ecoll_c" in in_p
            ):
                # or in_p in physical_params):
                continue
            all_params_list.append(in_p)

    if args.create_trainset:
        store_appended_data = "./"
        if args.store_appended_data is not None:
            if ".nc" in args.store_appended_data:
                store_appended_data = args.store_appended_data[:-3]
            else:
                store_appended_data = args.store_appended_data

        seg_thresholds = np.arange(
            args.min_threshold, args.max_threshold, args.threshold_step
        )

        for threshold in seg_thresholds:
            train, labels, test, test_labels = create_big_stratified_set(
                data_path=args.data_path,
                step_tol=args.step_tol,
                all_params_list=all_params_list,
                n_trajs_iter=args.n_trajs,
                out_params=out_params,
                threshold=threshold,
                distinct_outparams=False,
                cooldown=args.cooldown,
                min_time=-1000,  # Ensembles actually start here
                verbosity=args.verbosity,
            )

            np.save(
                file=f"{store_appended_data}_thresh{threshold:.1e}_cool{args.cooldown}_train",
                arr=train,
                allow_pickle=True,
                fix_imports=False,
            )  # I refuse to support Python2 in this year and age
            np.save(
                file=f"{store_appended_data}_thresh{threshold:.1e}_cool{args.cooldown}_labels",
                arr=labels,
                allow_pickle=True,
                fix_imports=False,
            )
            np.save(
                file=f"{store_appended_data}_thresh{threshold:.1e}_cool{args.cooldown}_test",
                arr=test,
                allow_pickle=True,
                fix_imports=False,
            )
            np.save(
                file=f"{store_appended_data}_thresh{threshold:.1e}_cool{args.cooldown}_test_labels",
                arr=test_labels,
                allow_pickle=True,
                fix_imports=False,
            )
        exit()

    data = parse_load(
        data_path=args.data_path,
        out_params=out_params,
        all_params_list=all_params_list,
        time_integrated=args.time_integrated,
        store_many_appended_data=args.store_many_appended_data,
        load_on_the_fly=args.load_on_the_fly,
        min_time=-1000,  # Ensembles actually start here
        verbosity=args.verbosity,
    )

    if args.verbosity > 0:
        print("The loaded dataset")
        print(data)

    if args.store_appended_data is not None:
        if ".nc" not in args.store_appended_data:
            print(f"You must add '.nc' to {args.store_appended_data}.")
            print(f"Result will be saved as {args.store_appended_data}.nc")
            store_appended_data = args.store_appended_data + ".nc"
        else:
            store_appended_data = args.store_appended_data

    if args.store_appended_data is not None:
        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in data.data_vars}
        index = store_appended_data.rfind("/")
        store_path = store_appended_data[:index]
        if not os.path.isdir(store_path):
            os.mkdir(store_path)
        data.to_netcdf(
            path=f"{store_appended_data}",
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
    if args.only_append:
        exit()

    n_trajs = args.n_trajs

    if args.classifier == "adaboost":
        max_features_list = [""]
    else:
        max_features_list = args.feature_split

    distinct_outparams = False
    step_tol = args.step_tol

    # Train the models first or load them
    seg_thresholds = np.arange(
        args.min_threshold, args.max_threshold, args.threshold_step
    )

    if args.model_path is not None:
        models = []
        for feat_idx, max_features in enumerate(max_features_list):
            models_inner = []
            forest_verbosity = 0
            if args.verbosity > 3:
                forest_verbosity = 2
            for seg_idx, seg_thresh in enumerate(seg_thresholds):
                model = load(
                    f"{args.model_path}_{args.classifier}_{max_features}_{seg_thresh:.1e}_{args.cooldown}.joblid"
                )
                if args.classifier == "random_forest":
                    model = model.set_params(**{"estimator__verbose": forest_verbosity})
                models_inner.append(model)
            models.append(models_inner)
    else:
        for i in range(len(max_features_list)):
            if args.classifier != "adaboost":
                try:
                    f = float(max_features_list[i])
                    max_features_list[i] = f
                except:
                    pass
        if args.train_subset:
            models, seg_thresholds = train_many_models(
                data=data.sel({"trajectory": data["trajectory"][0:n_trajs]}),
                max_features_list=max_features_list,
                min_threshold=args.min_threshold,
                max_threshold=args.max_threshold,
                threshold_step=args.threshold_step,
                precalc=True,
                step_tol=step_tol,
                distinct_outparams=distinct_outparams,
                n_estimators=args.n_estimators,
                no_split=args.no_split,
                save_memory=args.save_memory,
                verbosity=args.verbosity,
                max_depth=args.max_depth,
                max_leaf_nodes=args.max_leaf_nodes,
                classifier=args.classifier,
                learning_rate=args.learning_rate,
                cooldown=args.cooldown,
            )
        else:
            # This branch can make use of load-on-the-fly.
            models, seg_thresholds = train_many_models(
                data=data,
                max_features_list=max_features_list,
                min_threshold=args.min_threshold,
                max_threshold=args.max_threshold,
                threshold_step=args.threshold_step,
                precalc=True,
                step_tol=step_tol,
                distinct_outparams=distinct_outparams,
                n_estimators=args.n_estimators,
                no_split=args.no_split,
                save_memory=args.save_memory,
                verbosity=args.verbosity,
                max_depth=args.max_depth,
                max_leaf_nodes=args.max_leaf_nodes,
                classifier=args.classifier,
                learning_rate=args.learning_rate,
                cooldown=args.cooldown,
            )
        if args.store_models is not None:
            for feat_idx, max_features in enumerate(max_features_list):
                for seg_idx, seg_thresh in enumerate(seg_thresholds):
                    if max_features is None:
                        max_features = "None"
                    dump(
                        models[feat_idx][seg_idx],
                        f"{args.store_models}_{args.classifier}_{max_features}_{seg_thresh:.1e}_{args.cooldown}.joblid",
                    )
    if args.only_training:
        exit()

    traj_idx = 0
    # Get the confusion matrix for the training set if possible
    if args.train_subset:
        confus_train = create_dataset_pretrained(
            data=data.sel({"trajectory": data["trajectory"][0:n_trajs]}),
            models=models,
            seg_thresholds=seg_thresholds,
            max_features_list=max_features_list,
            step_tol=step_tol,
            distinct_outparams=distinct_outparams,
            precalc=True,
            cooldown=args.cooldown,
            verbosity=args.verbosity,
            independent_windows=args.independent_windows,
            in_params=all_params_list,
        )
        traj_idx = n_trajs

    # create confusion matrices (without the training set if possible)
    if not isinstance(data, list) and not isinstance(data, np.ndarray):
        max_traj_idx = len(data["trajectory"])
        confus_matrix = None
        if args.verbosity > 2:
            t_start = timer()

        while traj_idx < max_traj_idx:
            traj_idx_end = traj_idx + n_trajs
            if traj_idx_end > max_traj_idx:
                traj_idx_end = max_traj_idx
            if args.verbosity > 2:
                print(
                    f"Predicting trajectories {traj_idx} to {traj_idx_end} of {max_traj_idx}"
                )
                t1 = timer()
            confus_tmp = create_dataset_pretrained(
                data=data.sel(
                    {"trajectory": data["trajectory"][traj_idx:traj_idx_end]}
                ),
                models=models,
                seg_thresholds=seg_thresholds,
                max_features_list=max_features_list,
                step_tol=step_tol,
                distinct_outparams=distinct_outparams,
                precalc=True,
                cooldown=args.cooldown,
                verbosity=args.verbosity,
                independent_windows=args.independent_windows,
                in_params=all_params_list,
            )

            if confus_matrix is None:
                confus_matrix = confus_tmp
            else:
                confus_matrix += confus_tmp
            if args.verbosity > 2:
                t2 = timer()
                t_est = (t2 - t1) * (max_traj_idx / n_trajs - traj_idx_end / n_trajs)
                print(f"Done in {t2-t1} s (total {t2-t_start} s; est end in {t_est} s)")
            traj_idx = traj_idx_end
    else:
        # This path can make use of args.load_on_the_fly
        if args.verbosity > 2:
            print("Predicting starts")
        confus_matrix = create_dataset_pretrained(
            data=data,
            models=models,
            seg_thresholds=seg_thresholds,
            max_features_list=max_features_list,
            step_tol=step_tol,
            distinct_outparams=distinct_outparams,
            verbosity=args.verbosity,
            in_params=all_params_list,
            predict_train=args.predict_trainset,
            cooldown=args.cooldown,
            independent_windows=args.independent_windows,
        )

    confus_matrix["TP over P"] = confus_matrix["TP"] / confus_matrix["P"]
    confus_matrix["FP over P"] = confus_matrix["FP"] / confus_matrix["P"]
    confus_matrix["FN over P_windows"] = (
        confus_matrix["FN"] / confus_matrix["P_windows"]
    )
    if args.verbosity > 0:
        print("Confus matrix done")
        print(confus_matrix)
    comp = dict(zlib=True, complevel=9)
    encoding = {var: comp for var in confus_matrix.data_vars}
    try:
        confus_matrix.to_netcdf(
            path=args.store_name,
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
    except:
        print("Unable to write file. Trying without encoding")
        confus_matrix.to_netcdf(
            path=args.store_name,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
