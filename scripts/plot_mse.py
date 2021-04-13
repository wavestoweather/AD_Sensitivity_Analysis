from glob import glob
from multiprocessing import Pool
import numpy as np
import pandas as pd
from progressbar import progressbar as pb
import xarray as xr

try:
    from Deriv_dask import Deriv_dask
    from latexify import (
        in_params_dic,
        physical_params,
        in_params_grouping,
        in_params_descr_dic,
    )
    from segment_identifier import d_unnamed
    from create_mse import load_and_append
except:
    from scripts.Deriv_dask import Deriv_dask
    from scripts.latexify import (
        in_params_dic,
        physical_params,
        in_params_grouping,
        in_params_descr_dic,
    )
    from scripts.segment_identifier import d_unnamed
    from scripts.create_mse import load_and_append


def reduce_df(df, error_key, sens_kind="mean", error_kind="mean"):
    """
    Given an appended dataframe, reduce the (predicted) errors such that a
    single datapoint for each output parameter, perturbed parameter
    and ratio type is left.

    Parameters
    ----------
    df : pandas.Dataframe

    error_key : string
        Name of the column for predicted errors.
    sens_kind : string
        How to reduce the predicted errors over all trajectories. Options are:
        sum: Take the sum of all errors.
        max: Take the maximum of all errors.
        mean: Take the mean of all errors.
    error_kind : string
        How to reduce the true errors from perturbed ensembles. Options are:
        sum: Take the sum of all errors.
        max: Take the maximum of all errors.
        mean: Take the mean of all errors.

    Returns
    -------
    Reduced pandas.Dataframe
    """
    # get max or mean sensitivity and max or mean error
    if sens_kind == error_kind and error_kind == "sum":
        return (
            df.groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .sum()
            .reset_index()
        )
    elif sens_kind == error_kind and error_kind == "max":
        tmp_df = df.copy()
        tmp_df["Sensitivity"] = np.abs(tmp_df["Sensitivity"])
        return (
            tmp_df.groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .max()
            .reset_index()
        )
    elif sens_kind == error_kind and error_kind == "mean":
        return (
            df.groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .mean()
            .reset_index()
        )

    if sens_kind == "sum":
        sens_df = (
            df[["Output Parameter", "Perturbed Parameter", "Ratio Type", "Sensitivity"]]
            .groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .sum()
        )
    elif sens_kind == "max":
        tmp_df = df.copy()
        tmp_df["Sensitivity"] = np.abs(tmp_df["Sensitivity"])
        sens_df = (
            tmp_df[
                ["Output Parameter", "Perturbed Parameter", "Ratio Type", "Sensitivity"]
            ]
            .groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .max()
        )
    elif sens_kind == "mean":
        sens_df = (
            df[["Output Parameter", "Perturbed Parameter", "Ratio Type", "Sensitivity"]]
            .groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .mean()
        )

    if error_kind == "sum":
        err_df = (
            df[["Output Parameter", "Perturbed Parameter", "Ratio Type", error_key]]
            .groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .sum()
        )
    elif error_kind == "max":
        err_df = (
            df[["Output Parameter", "Perturbed Parameter", "Ratio Type", error_key]]
            .groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .max()
        )
    elif error_kind == "mean":
        err_df = (
            df[["Output Parameter", "Perturbed Parameter", "Ratio Type", error_key]]
            .groupby(["Output Parameter", "Perturbed Parameter", "Ratio Type"])
            .mean()
        )

    return pd.merge(
        sens_df, err_df, how="left", left_index=True, right_index=True
    ).reset_index()


def load_and_plot(
    kind,
    sens_kind,
    error_kind,
    plot_types,
    out_params,
    ratio_type,
    backend,
    filetype,
    files,
    store_path,
):
    """
    Load several files from intermediate results from create_mse.py and
    reduce them and finally plot the results.

    Parameters
    ----------
    kind : string
        Define the operation that had been used to calculate the error of
        the ensembles. Possible options:
        mse: Mean squared error
        maxse: Maximum error
        nozeromse: Mean squared error ignoring time steps with zero error
        sum: Cumulative squared error
        me: Mean error
        mae: Mean absolute error
    sens_kind : string
        Operation to summarize all sensitivities (aka predicted errors) across
        multiple ensembles. Possible options:
        sum: Sum of all errors
        max: Maximum error across all ensembles
        mean: Mean error
    error_kind : string
        Operation to summarize all deviations of ensembles from an
        unperturbed trajectory across multiple ensembles. Possible options:
        sum: Sum of all deviations
        max: Maximum deviation across all ensembles
        mean: Mean deviation
    plot_types : bool
        Wether to plot the errors grouped by the type of error.
    out_params : list of string
        List of output parameters to plot the datapoints for.
    ratio_type : string
        "vanilla": Use the derivative ratio in the file.
        "adjusted": Can be added to any to any type below where the sensitivity
        is adjusted to the parameter value such that perturbing this value by a
        certain percentage gives an approximation of the resulting error/difference
        for any given hydrometeor.
        "per_timestep": Use the highest derivative per timestep as denominator.
        "window": Use the highest derivative in the given window by min_x and max_x.
        "per_xaxis": Use the highest derivative per x_axis value. If x_axis is "timestep"
        it is the same as "per_timestep".
        "x_per_out_param": Replace 'x' with any other option than "vanilla". Use the highest
        derivative but per output parameter. (that *should* be the vanilla version)
        "x_weighted": Append this to any ratio type to use the inverse of the
        output parameter value as weight. Only works with "x_per_out_param".
    backend : string
        Either "matplotlib" or "bokeh" for plotting.
    filetype : string
        Define which type to load. Options are "csv" and "netcdf".
    """

    if kind == "mse":
        error_key = "MSE"
    elif kind == "maxse":
        error_key = "Max Error"
    elif kind == "nozeromse":
        error_key = "MSE (no zero)"
    elif kind == "sum":
        error_key = "Cumulative Squared Error"
    elif kind == "me":
        error_key = "Mean Error"
    elif kind == "mae":
        error_key = "Mean Absolute Error"
    all_df = load_and_append(files, filetype)

    in_params = np.unique(all_df["Perturbed Parameter"])
    # clean up the data. Some parameters might have slipped through that
    # are not used or not tracked with AD.
    del_list = []
    for in_p in in_params:
        if (
            "Not used" == in_params_descr_dic[in_p]
            or "Is not tracked with AD" in in_params_descr_dic[in_p]
        ):
            del_list.append(in_p)
    all_df = all_df[~all_df["Perturbed Parameter"].isin(del_list)]

    reduced_df = reduce_df(
        all_df, error_key, sens_kind=sens_kind, error_kind=error_kind
    )

    # out_params = ["QV"]
    if plot_types:
        # We need to add a column 'Group' to plot it correctly
        tmp_dic = {}
        for in_p in in_params:
            for g in in_params_grouping:
                if in_p in in_params_grouping[g]:
                    tmp_dic[in_p] = g
                    break

        reduced_df["Group"] = reduced_df.apply(
            lambda row: tmp_dic[row["Perturbed Parameter"]], axis=1
        )

    datashade = False
    alpha = 0.5
    s = 6
    f_limits = (-2, 2)
    confidence = 0.90
    hist = True
    plot_kind = "grid_plot"

    # Dummy for plotting
    mean_traj = Deriv_dask(
        direc="",
        parquet=False,
        netcdf=True,
        columns=None,
        backend=backend,
        file_ending="",
    )

    tmp_df = reduced_df.loc[reduced_df["Sensitivity"] != 0]
    if ratio_type is not None:
        tmp_df = reduced_df.loc[reduced_df["Ratio Type"] == ratio_type]

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=hist,
        confidence=None,
        kind=plot_kind,
        error_key=error_key,
        prefix=store_path + kind + "_s_" + sens_kind + "_e_" + error_kind + "_",
    )

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=hist,
        log_x=False,
        log_y=True,
        kind=plot_kind,
        error_key=error_key,
        prefix=store_path + kind + "_s_" + sens_kind + "_e_" + error_kind + "_ly",
    )

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=False,
        confidence=confidence,
        abs_x=False,
        log_x=True,
        log_y=True,
        kind=plot_kind,
        error_key=error_key,
        prefix=store_path + kind + "_s_" + sens_kind + "_e_" + error_kind + "_lxly",
    )

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=False,
        confidence=confidence,
        abs_x=False,
        log_x=True,
        log_y=False,
        kind=plot_kind,
        error_key=error_key,
        prefix=store_path + kind + "_s_" + sens_kind + "_e_" + error_kind + "_lx",
    )

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=False,
        confidence=confidence,
        abs_x=True,
        log_x=True,
        log_y=True,
        kind=plot_kind,
        error_key=error_key,
        prefix=store_path + kind + "_s_" + sens_kind + "_e_" + error_kind + "_lxlyabs",
    )

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=True,
        confidence=confidence,
        abs_x=True,
        log_x=True,
        log_y=True,
        kind=plot_kind,
        error_key=error_key,
        plot_singles=True,
        prefix=store_path
        + kind
        + "_s_"
        + sens_kind
        + "_e_"
        + error_kind
        + "_lxlyabshist",
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Plot mean squared error over maximal sensitivity calculated
        via the given ratio type (i.e. sensitivity per timestep) and given
        the trajectory to use for sensitivities and mean (i.e. the original trajectory
        from which model parameters have been perturbed).
        """
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="""
        Plot all kinds of error calculations and error predictions and the
        operation used to get a synopsis. This can take a long time!
        """,
    )
    parser.add_argument(
        "--kind",
        default="mse",
        help="""
        Define the operation that had been used to calculate the error of
        the ensembles. Possible options:
        mse: Mean squared error
        maxse: Maximum error
        nozeromse: Mean squared error ignoring time steps with zero error
        sum: Cumulative squared error
        me: Mean error
        mae: Mean absolute error
        """,
    )
    parser.add_argument(
        "--sens_kind",
        default="mean",
        help="""
        Operation to summarize all sensitivities (aka predicted errors) across
        multiple ensembles. Possible options:
        sum: Sum of all errors
        max: Maximum error across all ensembles
        mean: Mean error
        """,
    )
    parser.add_argument(
        "--error_kind",
        default="mean",
        help="""
        Operation to summarize all deviations of ensembles from an
        unperturbed trajectory across multiple ensembles. Possible options:
        sum: Sum of all deviations
        max: Maximum deviation across all ensembles
        mean: Mean deviation
        """,
    )
    parser.add_argument(
        "--plot_types",
        action="store_true",
        help="""
        If true: Plot input parameter types in different colors.
        Types means here wether a parameter is physical or rather artificial.
        """,
    )
    parser.add_argument(
        "--out_parameter",
        type=str,
        nargs="+",
        required=True,
        help="""
        Output parameter to plot for. Multiple options are possible.
        Possible options are
        QV, QC, QR, QG, QH, QS, QI.
        """,
    )
    parser.add_argument(
        "--ratio_type",
        default=None,
        help="""
        Plot only the given ratio type. If none is given, plot all ratio types.
        We recommend 'adjusted'.
        """,
    )
    parser.add_argument(
        "--backend",
        default="matplotlib",
        help="""
        Choose a backend for plotting. Options are:
        matplotlib: Most plots should be fine with it.
        bokeh: Can be used in case colors are wrong.
        """,
    )
    parser.add_argument(
        "--filetype",
        default="csv",
        type=str,
        help="""
        Define which type to load. Options are "csv" and "netcdf".
        """,
    )
    parser.add_argument(
        "--files",
        nargs="+",
        default=[
            "stats_full/mse_adjusted_conv_400_0.csv",
            "stats_full/mse_adjusted_conv_600_0.csv",
            "stats_full/mse_adjusted_conv_600_2.csv",
            "stats_full/mse_adjusted_conv_600_3.csv",
        ],
        help="""
        Define the files to load.
        """,
    )
    parser.add_argument(
        "--store_path",
        default="../pics/mse_",
        type=str,
        help="""
        Path to store the generated images.
        """,
    )

    args = parser.parse_args()

    if args.all:
        from itertools import product
        from timeit import default_timer as timer

        kind = ["mse", "maxse", "nozeromse", "sum", "me", "mae"]
        sens_kind = ["max", "mean"]
        error_kind = ["max", "mean"]
        n = len(kind) * len(sens_kind) * len(error_kind)
        i = 0
        t_start = timer()
        for k, s, e in product(kind, sens_kind, error_kind):
            t = timer()
            load_and_plot(
                k,
                s,
                e,
                args.plot_types,
                args.out_parameter,
                args.ratio_type,
                args.backend,
                args.filetype,
                args.files,
                args.store_path,
            )
            t2 = timer()
            ti = t2 - t
            i += 1
            t_all = (t2 - t_start) / i * n
            print(
                f"\n{i:3d}/{n} in {ti:3.3f}s, estimated end total {t_all:3.3f}s, total {t2-t_start:3.3f}s\n"
            )
    else:
        load_and_plot(
            args.kind,
            args.sens_kind,
            args.error_kind,
            args.plot_types,
            args.out_parameter,
            args.ratio_type,
            args.backend,
            args.filetype,
            args.files,
            args.store_path,
        )
