import holoviews as hv
import itertools
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import os.path
import pandas as pd
from timeit import default_timer as timer
import warnings
import xarray as xr

from sklearn.ensemble import RandomForestClassifier
from skmultilearn.model_selection import iterative_train_test_split
from itertools import product

try:
    from latexify import in_params_numeric_value_dic, parse_word
except:
    from scripts.latexify import in_params_numeric_value_dic, parse_word


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


def load_sensitivity(path, out_params, in_params=["da_1"]):
    """
    Load the sensitivities from the unperturbed trajectory.

    Parameters
    ----------
    path: String
        path and name of file to load. Should be the unperturbed
        trajectory with 'Output Parameter' as column.
    out_params: List of string
        List of output parameters to get sensitivities for.
    in_params: List of string
        List of input parameters for the sensitivies.
    Returns
    -------
    xarray.DataFrame
        Dataframe with model state parameters, sensitivities
        and time after ascent.
    """
    ds = xr.open_dataset(path, decode_times=False)
    ds = ds.loc[{"Output Parameter": out_params}][
        out_params + in_params + ["time_after_ascent"]
    ]
    for in_p in in_params:
        ds[in_p] = ds[in_p] * (in_params_numeric_value_dic[in_p] * 0.1)
    return ds


def load_unperturbed(path, out_params):
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
    ds = load_sensitivity(path, out_params)
    return ds[out_params + ["time_after_ascent"]]


def load_dataset(path, out_params, in_params, traj_list, traj_offset=0, verbosity=0):
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
    out_params: List of string
        List of output parameters to get MSE for.
    in_params: List of string
        List of input parameters to get MSE for.
    traj_list: List of int
        List of trajectory numbers to load.
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
    dim_order = (
        "Output Parameter",
        "Input Parameter",
        "trajectory",
        "time_after_ascent",
    )
    idx_offset = 0
    n_timesteps = 0
    time_after_ascent = None

    for traj in traj_list:
        if not os.path.isfile(
            path + "traj" + str(traj) + "/" + in_params[0][1::] + ".nc_wcb"
        ):
            continue
        trajectories.append(traj + traj_offset)
        if n_timesteps == 0:
            not_perturbed_path = path + "traj" + str(traj) + "_notPerturbed.nc_wcb"
            # I wasn't consistent with the way not perturbed parameters are stored
            if not os.path.isfile(not_perturbed_path):
                not_perturbed_path = path + "traj" + str(traj) + "/_notPerturbed.nc_wcb"

            val_df = load_sensitivity(not_perturbed_path, out_params, in_params)
            val_array = val_df.loc[{"Output Parameter": out_params[0]}][out_params[0]]
            n_timesteps = len(val_array["time"])
            time_after_ascent = np.asarray(
                val_df.loc[{"Output Parameter": out_params[0]}]["time_after_ascent"]
            ).flatten()

    mse = np.zeros((len(out_params), len(in_params), len(trajectories), n_timesteps))
    predicted = np.zeros(
        (len(out_params), len(in_params), len(trajectories), n_timesteps)
    )
    not_perturbed = np.zeros((len(out_params), len(trajectories), n_timesteps))

    for traj_idx, traj in enumerate(traj_list):
        if not os.path.isfile(
            path + "traj" + str(traj) + "/" + in_params[0][1::] + ".nc_wcb"
        ):
            idx_offset += 1
            continue

        not_perturbed_path = path + "traj" + str(traj) + "_notPerturbed.nc_wcb"
        # I wasn't consistent with the way not perturbed parameters are stored
        if not os.path.isfile(not_perturbed_path):
            not_perturbed_path = path + "traj" + str(traj) + "/_notPerturbed.nc_wcb"
        if verbosity > 1:
            print(f"Loading from {not_perturbed_path}")
        val_df = load_sensitivity(not_perturbed_path, out_params, in_params)
        val_only_df = val_df.loc[{"Output Parameter": out_params[0]}][out_params]
        for out_idx, out_p in enumerate(out_params):
            val_array = val_only_df[out_p]
            sens_df = val_df.loc[{"Output Parameter": out_p}][in_params]
            not_perturbed[out_idx, traj_idx - idx_offset, :] = np.asarray(
                val_array
            ).flatten()
            for in_idx, in_p in enumerate(in_params):
                ds = xr.open_dataset(
                    path + "traj" + str(traj) + "/" + in_p[1::] + ".nc_wcb",
                    decode_times=False,
                )
                tmp1 = np.asarray(ds[out_p].load())
                tmp2 = np.asarray(val_only_df[out_p])
                mse[out_idx, in_idx, traj_idx - idx_offset, :] = np.nanmean(
                    (tmp1 - tmp2) ** 2, axis=1
                ).flatten()
                predicted[out_idx, in_idx, traj_idx - idx_offset, :] = np.asarray(
                    sens_df[in_p]
                ).flatten()

    return xr.Dataset(
        data_vars={
            "Predicted Squared Error": (list(dim_order), predicted ** 2),
            "Mean Squared Error": (list(dim_order), mse),
            "Not Perturbed Value": (
                ["Output Parameter", "trajectory", "time_after_ascent"],
                not_perturbed,
            ),
        },
        coords={
            "Output Parameter": out_params,
            "Input Parameter": in_params,
            "trajectory": trajectories,
            "time_after_ascent": time_after_ascent,
        },
    )


def find_segments(df, error_threshold=0):
    """
    Iterate over time steps to mark the start of a segment with large errors, where large
    is defined via error_threshold. The dataframe needs to have a column (or index)
    "Output Parameter" and should be ordered by time or time_after_ascent.

    Parameters
    ----------
    df : DataFrame
        Dataframe with MSE created by load_dataset() for every model state
        parameter and all trajectories.
    error_threshold : float
        Threshold for errors to identify true segment starts.

    Returns
    -------
    Dataframe with additional column "segment_start" (1 for start, else 0)
    """

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
                final_pred = final_pred | df["detected_segment_jump"]

        if acc_ver:
            if final_pred is None:
                final_pred = df["detected_segment_acc"]
            else:
                final_pred = final_pred | df["detected_segment_acc"]
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

    Returns
    -------
    Dictionary with confusion matrices for every output parameter.
    """
    cached_data = find_segments(segment_data, segment_threshold)
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

    Returns
    -------
    xr.Dataset for statistics of predictions.
    Coordinates are "Output Parameter", "Input Parameter", "Default Threshold",
    "Jump Threshold" ("Acceleration Threshold" is not supported).
    Columns are "TP over P", "FP over P", "FN over P_windows", "TP", "P", "FP",
    "FN", "P_windows".
    """
    cached_data = find_segments(df, segment_threshold)
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
    n_windows = n_timesteps - (2 * step_tol + 1)
    X = ds["Predicted Squared Error"].values
    y = ds["segment_start"].values
    if distinct_outparams:
        feature_names = product(
            np.unique(ds["Output Parameter"].values),
            np.unique(ds["Input Parameter"].values),
        )
        class_names = np.unique(ds["Input Parameter"].values)
        X_final = np.concatenate(
            [
                *(
                    np.vstack(
                        X[out_p, :, traj, j : (step_tol * 2 + 1) + j]
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
                        y[out_p, :, traj, j : (step_tol * 2 + 1) + j]
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
    model = RandomForestClassifier(
        n_estimators=n_estimators,
        max_features=max_features,
        n_jobs=-1,
        verbose=verbose,
        max_depth=36,
        max_leaf_nodes=1000,
    )

    if isinstance(ds, list):
        y_final, X_final = ds
        model.fit(X_final, y_final)
        return model, X_final, None, y_final, None, None, None
    else:
        X_final, y_final, feature_names, class_names = create_input_labels(
            ds, step_tol, distinct_outparams
        )

    if no_split:
        model.fit(X_final, y_final)
        return model, X_final, None, y_final, None, feature_names, class_names

    if save_memory:
        model.fit(
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

    model.fit(train, train_labels)
    return model, train, test, train_labels, test_labels, feature_names, class_names


def get_tree_matrix(trained_model, X, y, only_idx=None):
    """
    Create a confusion matrix with a single row.

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

    Returns
    -------
    list of list with values as in the confusion matrix.
    """
    pred = trained_model.predict(X)

    # create the confusion matrix
    # The confusion matrix has the dimensions
    # tp, fn, fp, tn, ep, lp, p_act, p_act_win, pr, re, f1, fpr
    # Using the test set
    if only_idx is None:
        p_act_win = np.sum(y)
        n = np.sum((y == 0))  # np.product(np.shape(test_labels))
        tp = np.sum((pred == True) & (y == 1))
        fp = np.sum((pred == True) & (y == 0))
        p_pred = np.sum(pred)
        pr = tp / p_pred
        re = tp / p_act_win
    else:
        p_act_win = np.sum(y[:, only_idx])
        n = np.sum((y[:, only_idx] == 0))  # np.product(np.shape(test_labels))
        tp = np.sum((pred[:, only_idx] == True) & (y[:, only_idx] == 1))
        fp = np.sum((pred[:, only_idx] == True) & (y[:, only_idx] == 0))
        p_pred = np.sum(pred[:, only_idx])
        pr = tp / p_pred
        re = tp / p_act_win
    fn = p_act_win - tp
    tn = n - fp
    f1 = 2 * (pr * re) / (pr + re)
    fpr = fp / n
    confusion_matrix = [
        tp,
        fn,
        fp,
        tn,
        0,
        0,
        p_act_win,
        p_act_win,
        pr,
        re,
        f1,
        fpr,
        p_pred,
    ]
    return confusion_matrix


def plot_tree_matrix(
    model, train, test, train_labels, test_labels, feature_names, class_names
):
    """
    TODO: Plot the single row confusion matrix. Currently this method only
    returns a single row confusion matrix.
    """
    confusion_matrix = get_tree_matrix(model, test, test_labels)

    # TODO Plot results?
    return confusion_matrix


def show_tree_stats(
    model, train, test, train_labels, test_labels, feature_names, class_names
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

    confusion_matrix = get_tree_matrix(model, test, test_labels)

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
                ds_cache[seg_thresh] = find_segments(data, 10.0 ** seg_thresh)
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
                ds = find_segments(data, 10.0 ** seg_thresh)
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

    Returns
    -------
    xarray.Dataset of confusion matrix as described above.
    """

    if isinstance(data, list) or isinstance(data, np.ndarray):
        precalc = False
        load_data = True
    else:
        load_data = False
        in_params = list(np.unique(data["Input Parameter"]))
    dims_forest = (len(max_features_list), 1, 1 + len(in_params), len(seg_thresholds))
    coords_forest_dic = {
        "Max Features": max_features_list,
        "Output Parameter": ["All Output Parameters"],
        "Input Parameter": ["All Input Parameters Test Set"] + in_params,
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
    # ds_cache = {}
    test_cache = {}
    if precalc:
        for seg_thresh in seg_thresholds:
            ds = find_segments(data, 10.0 ** seg_thresh)
            test, test_labels, _, _ = create_input_labels(
                ds, step_tol, distinct_outparams
            )
            test_cache[seg_thresh] = (test, test_labels)
            # ds_cache[seg_thresh] = find_segments(data, 10.0**seg_thresh)

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
                    # print(path)
                    # print(test_labels)
                    test_labels = test_labels.astype(int, copy=False)
                else:
                    test = np.load(path, allow_pickle=True, fix_imports=False)
                    # print(path)
                    # print(test)
                if test is not None and test_labels is not None:
                    break
        elif precalc:
            test, test_labels = test_cache[seg_thresh]
            # ds = ds_cache[seg_thresh]
        else:
            ds = find_segments(data, 10.0 ** seg_thresh)
            test, test_labels, _, _ = create_input_labels(
                ds, step_tol, distinct_outparams
            )
        for feat_idx, max_features in enumerate(max_features_list):

            if verbosity > 2:
                print(
                    f"Running for {feat_idx}, {seg_idx} ({max_features}, {seg_thresh})",
                    end="",
                )
                t1 = timer()
            model = models[feat_idx][seg_idx]
            matrix = get_tree_matrix(model, test, test_labels)
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
                matrix = get_tree_matrix(model, test, test_labels, only_idx=in_idx)
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
            if verbosity > 2:
                if verbosity > 4:
                    print(f"{in_idx} ({in_p})")
                    print(
                        f"TP: {matrix[0]}; FP: {matrix[2]}; FN: {matrix[1]}; P_w: {matrix[7]}; P: {matrix[6]}"
                    )
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
            ds_cache[seg_thresh] = find_segments(data, 10.0 ** seg_thresh)

    for feat_idx, max_features in enumerate(max_features_list):
        for seg_idx, seg_thresh in enumerate(seg_thresholds):
            if precalc:
                ds = ds_cache[seg_thresh]
            else:
                ds = find_segments(data, 10.0 ** seg_thresh)
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
            )
            if verbosity > 2:
                print("Trained model")
            test_matrix = get_tree_matrix(model, test, test_labels)
            train_matrix = get_tree_matrix(model, train, train_labels)
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
                matrix = get_tree_matrix(model, test, test_labels, only_idx=in_idx)

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
        data = xr.open_dataset(data_path, decode_times=False)
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
            append_data(
                data
            )  # TODO: Why did I have that here? find_segments(data, 10.0**threshold))
            traj_offset += n_trajs
    return all_train, all_labels, all_test, all_test_labels


if __name__ == "__main__":
    import argparse
    from joblib import dump, load
    import os
    import progressbar as pb
    import sys

    parser = argparse.ArgumentParser(
        description="""
        Using perturbed ensembles: Create random forests to predict
        segments (start of deviations of ensembles from unperturbed trajectory)
        based on predicted errors from algorithmic differentiation (AD).
        """
    )
    parser.add_argument(
        "--load_models",
        dest="model_path",
        default=None,
        help="""
        Path to random forest models that had been created earlier.
        They are named as rand_forest_{max_features}_{seg_thresh:.1e}.joblid.
        """,
    )
    parser.add_argument(
        "--store_models",
        default=None,
        help="""
        If a path is given, store models on disk
        named as rand_forest_{max_features}_{seg_thresh:.1e}.joblid.
        """,
    )
    parser.add_argument(
        "--train_subset",
        action="store_true",
        help="""
        Train the models on a small subset of the data. This may be
        needed since stratifying the complete dataset can consume too much memory.
        Using a random sample of the complete dataset leads to poor performance
        since segment starts are rare and negative examples might be too
        dominant without startification.
        If 'no_split' is used, no stratification will be done and all
        data defined by 'n_trajs' is used for training. Otherwise only
        a quarter of 'n_trajs' is used for training but with stratification.
        """,
    )
    parser.add_argument(
        "--data_path",
        default="/data/project/wcb/netcdf/perturbed_ensembles/",
        help="""
        Path to folders with ensemble datasets or to single NetCDF file
        with all data concatenated along 'trajectory' axis.
        If a path to numpy arrays is given, it is assumed to be a training
        set for training only.
        """,
    )
    parser.add_argument(
        "--min_threshold",
        type=float,
        default=-40,
        help="""
        Minimum threshold (power of 10) for classifying a time step in a
        perturbed ensemble as a segment start.
        """,
    )
    parser.add_argument(
        "--max_threshold",
        type=float,
        default=0,
        help="""
        Maximum threshold (power of 10) for classifying a time step in a
        perturbed ensemble as a segment start.
        """,
    )
    parser.add_argument(
        "--threshold_step",
        type=float,
        default=2,
        help="""
        Step size from minimum threshold to maximum threshold for
        classifying a time step in a perturbed ensemble as a segment start.
        """,
    )
    parser.add_argument(
        "--step_tol",
        type=int,
        default=8,
        help="""
        Tolerance for predicting a segment start in time steps. The
        size of the resulting sliding window is step_tol*2+1.
        """,
    )
    parser.add_argument(
        "--n_estimators",
        type=int,
        default=100,
        help="""
        Number of estimators for the random forest.
        """,
    )
    parser.add_argument(
        "--n_trajs",
        type=int,
        default=10,
        help="""
        Number of trajectories that will be predicted using the trained models
        at once. If 'train_subset' is set, then this is the number of
        trajectories used for training the models.
        """,
    )
    parser.add_argument(
        "--store_name",
        default="prediction.nc",
        help="""
        Name of NetCDF file where results (confusion matrix) will be stored.
        Must end with '.nc'!
        """,
    )
    parser.add_argument(
        "--save_memory",
        action="store_true",
        help="""
        Save memory by splitting the dataset into a train and test set by
        using the first trajectories as training set without stratification.
        This can lead to poor performance of the models.
        """,
    )
    parser.add_argument(
        "--no_split",
        action="store_true",
        help="""
        Do not split the dataset into a training and test set.
        If 'train_subset' is used, then the subset defined by 'n_trajs' is
        used for training. Otherwise the complete dataset is used for training.
        Overrides 'save_memory'.
        """,
    )
    parser.add_argument(
        "--store_many_appended_data",
        default=None,
        help="""
        Store the appended input data to this path as NetCDF file for each appended
        version. Used mainly for debugging.
        """,
    )
    parser.add_argument(
        "--store_appended_data",
        default=None,
        help="""
        Store the final appended data to this path and name. Must end with
        '.nc' for datasets. For training sets, '.nc' will be stripped away.
        """,
    )
    parser.add_argument(
        "--only_append",
        action="store_true",
        help="""
        Only appending of data. Use 'store_appended_data' to define
        a path and name.
        """,
    )
    parser.add_argument(
        "--only_training",
        action="store_true",
        help="""
        Only model training. Use 'store_models' to define a path to save the
        models.
        """,
    )
    parser.add_argument(
        "--create_trainset",
        action="store_true",
        help="""
        Create a trainingset with labels and store it at 'store_appended_data'
        for training later.
        """,
    )
    parser.add_argument(
        "--predict_trainset",
        action="store_true",
        help="""
        Predict the training set created by 'create_trainset'.
        """,
    )
    parser.add_argument(
        "--feature_split",
        default=["log2", "sqrt", 1.0],
        nargs="+",
        help="""
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
        """,
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help="""
        Set verbosity level.
        0: No output except for exceptions
        1: Print datasets
        2: Print loading statements
        3: Print building/training statements
        4: Set verbosity of random forest to 2
        5: Get predicted results for each entry
        6: Get count for segment starts for every input parameter
        """,
    )

    args = parser.parse_args()
    if args.store_name is not None:
        if ".nc" not in args.store_name:
            print(f"You must add '.nc' to {args.store_name}.")
            print(f"Result will be saved as {args.store_name}.nc")
            store_name = args.store_name + ".nc"
        else:
            store_name = args.store_name

    out_params = ["QV", "QC", "QR", "QG", "QH", "QI", "QS"]
    # the following lines can give you the most important parameters
    # reduced_df = d_unnamed(pd.read_csv("../stats_full/conv_adjusted_mse_errorMean_sensMean.csv"))
    # reduced_df = reduced_df.loc[reduced_df["Ratio Type"] == "adjusted"]
    # reduced_df = reduced_df.loc[reduced_df["Sensitivity"] != 0]
    # top20_sens_dic = {}
    # all_params_list = []
    # for out_p in out_params:
    #     df = reduced_df.loc[reduced_df["Output Parameter"] == out_p]
    #     top20_sens_dic[out_p] = list(np.unique( df.nlargest(20, "Sensitivity")["Perturbed Parameter"] ))
    #     all_params_list.extend(top20_sens_dic[out_p])
    # del(reduced_df)
    # all_params_list = list(set(all_params_list))
    all_params_list = [
        "dc_ccn_4",
        "dice_b_geo",
        "da_prime",
        "dT_mult_max",
        "dsnow_b_f",
        "dsnow_c_s",
        "dsnow_b_vel",
        "dsnow_a_f",
        "dgraupel_vsedi_min",
        "dT_mult_min",
        "dgraupel_ecoll_c",
        "db_prime",
        "dice_c_s",
        "db_v",
        "dp_sat_melt",
        "drain_b_geo",
        "dhail_b_vel",
        "dsnow_b_geo",
        "dD_rainfrz_gh",
        "db_HET",
        "dgraupel_b_f",
        "dgraupel_b_vel",
        "drain_beta",
        "dice_b_f",
        "drain_mu",
        "dice_b_vel",
        "db_ccn_4",
        "dcloud_b_vel",
        "dsnow_a_geo",
        "dp_sat_ice_const_b",
        "da_HET",
        "drain_c_z",
        "drain_b_vel",
        "drain_a_geo",
        "dhail_b_geo",
        "dgraupel_c_s",
        "dgraupel_b_geo",
        "dice_a_geo",
        "drain_gamma",
        "dinv_z",
        "drain_a_vel",
        "dhail_a_vel",
        "dcloud_b_geo",
        "drain_alpha",
        "drain_cmu2",
        "drho_vel",
        "dhail_a_geo",
        "dgraupel_a_geo",
        "dsnow_a_vel",
        "drain_cmu4",
        "dkin_visc_air",
        "dk_r",
        "drain_cmu3",
        "dice_a_f",
        "drain_nu",
        "dgraupel_max_x",
        "drain_g2",
        "dice_a_vel",
        "dgraupel_a_vel",
    ]

    traj_offset = 0

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
                verbosity=args.verbosity,
            )

            np.save(
                file=f"{store_appended_data}_thresh{threshold:.1e}_train",
                arr=train,
                allow_pickle=True,
                fix_imports=False,
            )  # I refuse to support Python2 in this year and age
            np.save(
                file=f"{store_appended_data}_thresh{threshold:.1e}_labels",
                arr=labels,
                allow_pickle=True,
                fix_imports=False,
            )
            np.save(
                file=f"{store_appended_data}_thresh{threshold:.1e}_test",
                arr=test,
                allow_pickle=True,
                fix_imports=False,
            )
            np.save(
                file=f"{store_appended_data}_thresh{threshold:.1e}_test_labels",
                arr=test_labels,
                allow_pickle=True,
                fix_imports=False,
            )
        exit()

    if ".nc" in args.data_path:
        # ie data2_327.nc
        if args.verbosity > 1:
            print(f"Loading {args.data_path}")
        data = xr.open_dataset(args.data_path, decode_times=False)
    else:
        if args.verbosity > 1:
            print(f"Checking {args.data_path}")

        paths = list(os.listdir(args.data_path))
        if ".nc" in paths[0]:
            data = None
            if args.verbosity > 1:
                print(f"Loading from {paths}")
            for p in range(len(paths)):
                if "quan" in paths[p] or "median" in paths[p] or "perturb" in paths[p]:
                    continue
                path = args.data_path + paths[p] + "/"
                n_trajs = len(list(os.listdir(path)))
                if args.verbosity > 1:
                    print(f"Loading from {path} with {n_trajs} trajectories")
                try:
                    tmp = load_dataset(
                        path=path,
                        out_params=out_params,
                        in_params=all_params_list,
                        traj_list=np.arange(n_trajs),
                        traj_offset=traj_offset,
                        verbosity=args.verbosity,
                    )
                except:
                    print(f"Loading      {path} failed.")
                    tmp = None
                if tmp is None:
                    continue

                if data is None:
                    data = tmp
                else:
                    data = xr.concat([data, tmp], dim="trajectory", join="outer")
                if args.store_many_appended_data is not None:
                    comp = dict(zlib=True, complevel=9)
                    encoding = {var: comp for var in data.data_vars}
                    data.to_netcdf(
                        path=f"{args.store_many_appended_data}data_{traj_offset}.nc",
                        encoding=encoding,
                        compute=True,
                        engine="netcdf4",
                        format="NETCDF4",
                        mode="w",
                    )
                traj_offset += n_trajs
        else:  # numpy arrays with training data
            # create a dictionary for every threshold that holds a list
            # of X and y
            for i in range(len(paths)):
                paths[i] = args.data_path + paths[i]
            data = np.sort(paths)

    if args.verbosity > 0:
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

    max_features_list = (
        args.feature_split
    )  # [1.0]#["log2", "sqrt"]#, 0.25] # 0.25, 0.75
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
                    f"{args.model_path}rand_forest_{max_features}_{seg_thresh:.1e}.joblid"
                )
                model = model.set_params(**{"verbose": forest_verbosity})
                models_inner.append(model)
            models.append(models_inner)
    else:
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
            )
        else:
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
            )
        if args.store_models is not None:
            for feat_idx, max_features in enumerate(max_features_list):
                for seg_idx, seg_thresh in enumerate(seg_thresholds):
                    if max_features is None:
                        max_features = "None"
                    dump(
                        models[feat_idx][seg_idx],
                        f"{args.store_models}rand_forest_{max_features}_{seg_thresh:.1e}.joblid",
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
            verbosity=args.verbosity,
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
                verbosity=args.verbosity,
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
        if args.verbosity > 2:
            print("Predicting starts")
        # TODO: all_params_list might have wrong order
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
    # old version that does all at once without saving models
    # ds_forest = create_dataset_forest(
    #     data=data,
    #     max_features_list=["auto", "log2", "sqrt", None],
    #     min_threshold=-40,
    #     max_threshold=0,
    #     threshold_step=2,
    #     precalc=True,
    #     step_tol=8,
    #     test_size=0.25,
    #     distinct_outparams=False,
    #     n_estimators=100)
    # print("create_dataset_forest done")
    # comp = dict(zlib=True, complevel=9)
    # encoding = {var: comp for var in ds_forest.data_vars}
    # ds_forest.to_netcdf(
    #     path="cached_ds_forest_win8_many.nc",
    #     encoding=encoding,
    #     compute=True,
    #     engine="netcdf4",
    #     format="NETCDF4",
    #     mode="w")
