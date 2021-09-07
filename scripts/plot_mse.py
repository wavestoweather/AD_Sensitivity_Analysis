from bokeh.models import Range1d, LinearAxis, GlyphRenderer
from glob import glob
import holoviews as hv
from holoviews.operation.datashader import datashade as dsshade
from holoviews import opts
import hvplot.pandas
import matplotlib
from multiprocessing import Pool
import numpy as np
import os
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
        parse_word,
        get_unit,
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
        parse_word,
        get_unit,
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


def plot_histogram(
    df,
    out_params,
    store_path="pics/correlation",
    backend="bokeh",
    plot_types=True,
    add_zero_sens=False,
    title=None,
    xlabel=None,
    xlabel2=None,
    width=900,
    height=900,
):
    """
    Plot the dataframe which should hold parameters with their sensitivity
    to one model state parameter and the actual deviation when perturbing
    on a log histogram plot.

    Parameters
    ----------
    df : pandas.Dataframe
        A dataframe with columns "Output Parameter" for the model state
        variables, which must have only one value,
        "Input Parameter" for the perturbed model parameter,
        "Predicted Squared Error" for the sensitivity calculated to deviations
        in the next timestep, "Mean Squared Error" for the actual deviations
    out_params : list of string
        *Should* be only a list of a single value, namely the output parameter
        in the given dataframe.
    store_path : string
        Path to folder where to store the images.
    backend : string
        Either "matplotlib" or "bokeh" for plotting.
    plot_types : bool
        Wether to plot the errors grouped by the type of the model parameters.
    add_zero_sens : bool
        Add sensitivities of value zero to the far left and mark it with
        negative infinity.
    title : string
        Title for the plot.
    xlabel : string
        Alternative label for x-axis of upper plot.
    xlabel2 : string
        Alternative label for x-axis of bottom plot.
    width : int
        Width of plot in pixels.
    height : int
        Height of plot in pixels.
    """
    in_params = np.unique(df["Input Parameter"])

    datashade = False
    alpha = 0.5
    f_limits = (-2, 2)
    plot_kind = "paper"

    error_key = "Mean Squared Error"
    sens_key = "Predicted Squared Error"

    fontscale = width / 350
    matplotlib.rcParams["axes.formatter.limits"] = f_limits

    mse_df = df.copy()
    if plot_types:
        # We need to add a column 'Group' to plot it correctly
        tmp_dic = {}
        for in_p in in_params:
            for g in in_params_grouping:
                if in_p in in_params_grouping[g]:
                    tmp_dic[in_p] = g
                    break

        mse_df["Group"] = mse_df.apply(
            lambda row: tmp_dic[row["Input Parameter"]], axis=1
        )

    if add_zero_sens:
        inf_val = np.min(mse_df.loc[mse_df[sens_key] != 0][sens_key]) / 10
        mse_df[sens_key] = mse_df[sens_key].replace({0: inf_val})
    else:
        mse_df = mse_df.loc[df[sens_key] != 0]
        inf_val = None

    hv.extension(backend)

    # log_x
    mse_df[sens_key] = np.log10(mse_df[sens_key])
    if inf_val is not None:
        inf_val = np.log10(inf_val)
    # log_y
    mse_df[error_key] = np.abs(mse_df[error_key])
    mse_df[error_key] = np.log10(mse_df[error_key])

    if xlabel is None or xlabel == "":
        xlabel = "Log MSD"
    if xlabel2 is None or xlabel2 == "":
        xlabel2 = r"Log Predicted MSD"

    if title is None:
        title = "Deviation by Perturbed Parameter"
    # Dummy for colormappings
    dummy = Deriv_dask(
        direc="",
        parquet=False,
        netcdf=True,
        columns=None,
        backend=backend,
        file_ending="",
    )
    cmap_values = []
    groups = np.unique(mse_df["Group"])
    for group in groups:
        cmap_values.append(matplotlib.colors.to_hex(dummy.cmap_types[group]))
    by_col = "Group"
    if backend == "bokeh":
        image = (
            mse_df.hvplot.hist(
                y=error_key,
                by=by_col,
                alpha=alpha,
                legend=False,
                color=cmap_values,
                width=width,
                height=int(height / 2),
                xlabel=xlabel,
            ).opts(fontscale=fontscale)
            + mse_df.hvplot.hist(
                y=sens_key,
                by=by_col,
                alpha=alpha,
                legend="top_left",
                color=cmap_values,
                width=width,
                height=int(height / 2),
                xlabel=xlabel2,
            ).opts(fontscale=fontscale)
        )
        image = image.cols(1)
    else:
        plot1 = mse_df.hvplot.hist(
            y=error_key,
            by=by_col,
            alpha=alpha,
            legend=False,
            color=cmap_values,
            xlabel=xlabel,
        ).opts(fontscale=fontscale)

        plot2 = mse_df.hvplot.hist(
            y=sens_key,
            by=by_col,
            alpha=alpha,
            legend=False,
            color=cmap_values,
            xlabel=xlabel2,
        ).opts(fontscale=fontscale)

        legend_overlay = hv.NdOverlay(
            {
                groups[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                    opts.Scatter(
                        color=cmap_values[i],
                    )
                )
                for i in range(len(groups))
            }
        )
        image = plot1 + plot2 * legend_overlay

    if backend == "matplotlib":
        image = image.opts(
            fig_inches=((width / 300.0, height / 300.0)),
            fontscale=fontscale,
        )
    else:
        image = image.opts(
            width=width,
            height=height,
            title=title,
            fontscale=fontscale,
        )

    renderer = hv.Store.renderers[backend].instance(fig="png", dpi=300)
    i = 0
    save = store_path + "_" + "{:03d}".format(i)
    while os.path.isfile(save + ".png"):
        i = i + 1
        save = store_path + "_" + "{:03d}".format(i)
    renderer.save(image, save)


def plot_mse(
    df,
    out_params,
    store_path="pics/correlation",
    backend="bokeh",
    plot_types=True,
    add_zero_sens=False,
    confidence=0.90,
    title=None,
    xlabel="Sensitivity",
    ylabel=None,
    width=900,
    height=900,
    hist=True,
):
    """
    Plot the dataframe which should hold parameters with their sensitivity
    to one model state parameter and the actual deviation when perturbing
    on a log-log plot.

    Parameters
    ----------
    df : pandas.Dataframe
        A dataframe with columns "Output Parameter" for the model state
        variables, which must have only one value,
        "Input Parameter" for the perturbed model parameter,
        "Predicted Squared Error" for the sensitivity calculated to deviations
        in the next timestep, "Mean Squared Error" for the actual deviations
    out_params : list of string
        *Should* be only a list of a single value, namely the output parameter
        in the given dataframe.
    store_path : string
        Path to folder where to store the images.
    backend : string
        Either "matplotlib" or "bokeh" for plotting.
    plot_types : bool
        Wether to plot the errors grouped by the type of the model parameters.
    add_zero_sens : bool
        Add sensitivities of value zero to the far left and mark it with
        negative infinity.
    confidence : float or None
        Plot a confidence ellipse around each sample with confidence
        between 0 and 1. If none is given, no ellipse will be plotted.
    title : string
        Title for the plot.
    xlabel : string
        Alternative label for x-axis.
    ylabel : string
        Alternative label for y-axis.
    width : int
        Width of plot in pixels.
    height : int
        Height of plot in pixels.
    hist : bool
        Plot histogram around the plot.
    """

    in_params = np.unique(df["Input Parameter"])

    datashade = False
    alpha = 1  # 0.5
    s = 12
    f_limits = (-2, 2)
    plot_kind = "paper"

    error_key = "Mean Squared Error"
    sens_key = "Predicted Squared Error"

    df_tmp = df.copy()
    if plot_types:
        # We need to add a column 'Group' to plot it correctly
        tmp_dic = {}
        for in_p in in_params:
            for g in in_params_grouping:
                if in_p in in_params_grouping[g]:
                    tmp_dic[in_p] = g
                    break

        df_tmp["Group"] = df_tmp.apply(
            lambda row: tmp_dic[row["Input Parameter"]], axis=1
        )

    # Dummy for plotting
    mean_traj = Deriv_dask(
        direc="",
        parquet=False,
        netcdf=True,
        columns=None,
        backend=backend,
        file_ending="",
    )
    if add_zero_sens:
        inf_val = np.min(df_tmp.loc[df_tmp[sens_key] != 0][sens_key]) / 10
        df_tmp[sens_key] = df_tmp[sens_key].replace({0: inf_val})
    else:
        df_tmp = df_tmp.loc[df[sens_key] != 0]
        inf_val = None

    if "NC" in title or "_OUT" in title:
        title = title.replace("Estimation", "Est.")
    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=df_tmp,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        height=height,
        width=width,
        hist=hist,
        confidence=confidence,
        log_func=np.log10,
        abs_x=True,
        log_x=True,
        log_y=True,
        kind=plot_kind,
        error_key=error_key,
        sens_key=sens_key,
        prefix="_s_e_lxlyabshist",
        title=title,
        linewidth=3,
        xlabel=xlabel,
        ylabel=ylabel,
        plot_path=store_path,
        inf_val=inf_val,
    )


def plot_time_evolution(
    df,
    store_path="../pics/time_evo",
    backend="matplotlib",
    title=None,
    xlabel=None,
    ylabel=None,
    twinlabel=None,
    min_x=None,
    max_x=None,
    plot_deviation=False,
    perturb_delta=None,
    alpha=1,
    width=900,
    height=900,
    logy=False,
    logtwin=False,
):
    """
    Plot over time after ascent the model state (left y-axis) and
    sensitivity (right y-axis).
    We cut off the highest predicted value, if the difference to the second
    biggest value is an order of magnitude or more,
    since this can be an outlier which
    makes the rest of the values hard to see.

    Parameters
    ----------
    df : pandas.Dataframe
        A dataframe with columns "Output Parameter" for the model state
        variables, which must have only one value,
        "Input Parameter" for the perturbed model parameter,
        "Predicted Squared Error" for the sensitivity calculated to deviations
        in the next timestep, "Mean Squared Error" for the actual deviations
    store_path : string
        Path to folder where to store the images.
    backend : string
        Either "matplotlib" or "bokeh" for plotting. Currently, only matplotlib
        gives correct plots.
    title : string
        Title for the plot.
    xlabel : string
        Alternative label for x-axis.
    ylabel : string
        Alternative label for y-axis.
    twinlabel : string
        Label for twin axis.
    min_x : float
        Set the start time relative to the start of the ascent. In case
        of plotting the deviation via perturbance it is a good idea to set
        min_x to the first step when it is actually perturbed.
    max_x : float
        Set the last time relative to the start of the ascent.
    plot_deviation : bool
        If true, plot the deviation from perturbing parameters. If false,
        plot the actual model state variable.
    perturb_delta : float
        Time difference between different ensembles where perturbing started.
        This will plot vertical lines at each of those time steps.
    alpha : float
        Alpha value for the dots.
    width : int
        Width of plot in pixels.
    height : int
        Height of plot in pixels.
    logy : bool
        If true, use log10 on y-axis
    logtwin : bool
        If true, use log10 on twin-axis
    """
    s = 14
    import matplotlib.ticker as tick

    fontscale = height / 350
    if backend == "matplotlib":
        fontscale = height / 900 * 0.8
        s = 2
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = "Helvetica"
        title_font_size = 13 * fontscale
        axis_label_text_font_size = matplotlib.rcParams["font.size"] * fontscale
        major_label_text_font_size = matplotlib.rcParams["font.size"] * fontscale
        aspect = width / height
    else:
        axis_label_text_font_size = int(13 * fontscale * 0.8)
        major_label_text_font_size = int(11 * fontscale * 0.8)

    hv.extension(backend)

    if min_x is not None:
        df = df.loc[df["time_after_ascent"] >= min_x]
    if max_x is not None:
        df = df.loc[df["time_after_ascent"] <= max_x]

    # Make time after ascent to minuts
    df["time_after_ascent"] /= 60
    y = "Not Perturbed Value"
    if plot_deviation:
        y = "Mean Squared Error"

    min_log_y = np.NaN
    if logy:
        df[y] = np.log10(df[y])
        min_log_y = np.nanmin(df.loc[df[y] != -np.inf][y]) - 1
        df.replace(-np.inf, min_log_y, inplace=True)
    min_log_twin = np.NaN
    if logtwin:
        df["Predicted Error"] = np.log10(df["Predicted Error"])
        min_log_twin = (
            np.nanmin(df.loc[df["Predicted Error"] != -np.inf]["Predicted Error"]) - 1
        )
        df.replace(-np.inf, min_log_twin, inplace=True)

    df = df.loc[df["Predicted Error"] < np.max(df["Predicted Error"])]

    lower_y = np.min(df["Predicted Error"])
    upper_y = np.max(df["Predicted Error"])
    delta = (upper_y - lower_y) / 12
    lower_y -= delta
    upper_y += delta

    lower_y2 = np.min(df[y])
    upper_y2 = np.max(df[y])
    delta2 = (upper_y2 - lower_y2) / 12
    lower_y2 -= delta2
    upper_y2 += delta2

    renderer = hv.Store.renderers[backend].instance(fig="png", dpi=300)

    min_t = np.min(df["time_after_ascent"])
    max_t = np.max(df["time_after_ascent"])

    if perturb_delta is not None:
        perturb_lines = None
        t = min_t
        while t < max_t:
            if perturb_lines is None:
                perturb_lines = hv.VLine(x=t).opts(color="black")
            else:
                perturb_lines *= hv.VLine(x=t).opts(color="black")
            t += perturb_delta

    if len(np.unique(df["Input Parameter"])) < 10:
        cmap = matplotlib.pyplot.get_cmap("tab10")
    else:
        cmap = matplotlib.pyplot.get_cmap("tab20")

    cmap_map = {}
    cmap_values = [matplotlib.colors.to_hex(cmap(0)[0:-1])]
    data_types = [parse_word(np.unique(df["Output Parameter"])[0])]
    data_types_2 = [np.unique(df["Output Parameter"])[0]]
    for i, rt in enumerate(np.unique(df["Input Parameter"])):
        data_types.append(parse_word(rt).replace(r"\text", r"\mathrm"))
        cmap_values.append(matplotlib.colors.to_hex(cmap(i + 1)[0:-1]))
        cmap_map[rt] = matplotlib.colors.to_hex(cmap(i + 1)[0:-1])
        data_types_2.append(rt)

    if backend == "matplotlib":

        def apply_axis_format(plot, element):
            ax = plot.handles["axis"]
            ax.set_title(
                title,
                loc="left",
                # fontdic={},
                **{
                    "fontweight": "bold",
                    "fontsize": title_font_size,
                },
            )
            ax.set_ylabel(
                ylabel,
                **{
                    "rotation": 90,
                    "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size,
                },
            )
            ax.set_xlabel(
                xlabel,
                **{
                    "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size,
                },
            )
            ax.minorticks_on()

            def format_fn(tick_val, tick_pose):
                # Can implement an even fancier format here if needed
                if logy and tick_val == min_log_y:
                    return "-Infinity"
                return np.format_float_scientific(
                    tick_val, precision=0, exp_digits=1
                ).replace(
                    ".", ""
                )  # f"{tick_val:1.0e}"

            ax.yaxis.set_major_formatter(tick.FuncFormatter(format_fn))
            ax.tick_params(
                axis="both", which="major", labelsize=major_label_text_font_size
            )
            plot.handles["axis"] = ax

        left = (
            df.loc[df["Input Parameter"] == data_types_2[1]]
            .hvplot.scatter(
                x="time_after_ascent",
                y=y,
                alpha=alpha,
                label=data_types[0],
                ylabel=ylabel,
                xlabel=xlabel,
                legend=False,
                aspect=aspect,
                color=[cmap_values[0]],
            )
            .opts(
                initial_hooks=[apply_axis_format],
            )
            .opts(
                opts.Scatter(
                    s=s,
                )
            )
        )

        def twinx_per_timestep(plot, element):
            ax = plot.handles["axis"]
            twinax = ax.twinx()
            twinax.set_ylim((lower_y, upper_y))
            twinax.set_ylabel(
                twinlabel,
                labelpad=9.0,
                **{
                    "rotation": 270,
                    "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size,
                },
            )
            if logtwin:
                yticks_tmp = twinax.get_yticks()
                yticks_tmp[0] = min_log_twin
                if yticks_tmp[1] - yticks_tmp[0] < 3.5:
                    yticks_tmp = np.delete(yticks_tmp, 1)
                twinax.set_yticks(yticks_tmp)
            else:
                twinax.minorticks_on()

            def format_fn(tick_val, tick_pose):
                # Can implement an even fancier format here if needed
                if logtwin and tick_val == min_log_twin:
                    return "-Infinity"
                return np.format_float_scientific(
                    tick_val, precision=0, exp_digits=1
                ).replace(
                    ".", ""
                )  # f"{tick_val:1.0e}"

            twinax.yaxis.set_major_formatter(tick.FuncFormatter(format_fn))
            twinax.tick_params(
                axis="both", which="major", labelsize=major_label_text_font_size
            )
            plot.handles["axis"] = twinax

        def twinx2(plot, element):
            ax = plot.handles["axis"]
            twinax = ax.twinx()
            twinax.set_ylim((lower_y, upper_y))
            twinax.set_yticks([])
            plot.handles["axis"] = twinax

        twin = None
        for i, in_p in enumerate(data_types_2[1::]):
            if twin is None:
                twin = (
                    df.loc[df["Input Parameter"] == in_p]
                    .hvplot.scatter(
                        x="time_after_ascent",
                        y="Predicted Error",
                        alpha=alpha,
                        legend=False,
                        aspect=aspect,
                        color=cmap_values[i + 1],
                    )
                    .opts(
                        initial_hooks=[twinx_per_timestep],
                        apply_ranges=False,
                    )
                    .opts(
                        opts.Scatter(
                            s=s,
                            cmap=cmap_map,
                        ),
                    )
                )
            else:
                twin = twin * df.loc[df["Input Parameter"] == in_p].hvplot.scatter(
                    x="time_after_ascent",
                    y="Predicted Error",
                    alpha=alpha,
                    legend=False,
                    aspect=aspect,
                    color=cmap_values[i + 1],
                ).opts(initial_hooks=[twinx2], apply_ranges=False,).opts(
                    opts.Scatter(
                        s=s,
                        cmap=cmap_map,
                    ),
                )
        legend_overlay = hv.NdOverlay(
            {
                data_types[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                    opts.Scatter(
                        s=s * 8,
                        color=cmap_values[i],
                    )
                )
                for i in range(len(data_types))
            }
        )
        if perturb_delta is not None:
            image = left * twin * perturb_lines * legend_overlay
        else:
            image = left * twin * legend_overlay
    else:

        left = df.hvplot.scatter(
            x="time_after_ascent",
            y=y,
            alpha=alpha,
            label=data_types[0],
            ylabel=ylabel,
            xlabel=xlabel,
            legend="top_right",
            size=s,
            color=[cmap_values[0]],
        )

        def apply_twin_axis(plot, element):
            p = plot.state
            p.extra_y_ranges = {"twiny": Range1d(start=lower_y, end=upper_y)}
            p.add_layout(
                LinearAxis(
                    y_range_name="twiny",
                    axis_label_text_font_size=f"{axis_label_text_font_size}pt",
                    major_label_text_font_size=f"{major_label_text_font_size}pt",
                    axis_label=twinlabel,
                ),
                "right",
            )
            glyph = p.select(dict(type=GlyphRenderer))[0]
            glyph.y_range_name = "twiny"

        df["Input Parameter"] = parse_word(df["Input Parameter"])
        twin = df.hvplot.scatter(
            x="time_after_ascent",
            y="Predicted Error",
            by="Input Parameter",
            alpha=alpha,
            size=s,
            colormap=cmap_values,
        ).opts(
            hooks=[apply_twin_axis],
            apply_ranges=True,
        )

        if perturb_delta is not None:
            image = left * twin * perturb_lines
        else:
            image = left * twin
    if backend == "matplotlib":
        image = image.opts(
            fig_inches=((width / 300.0, height / 300.0)),
            fontscale=fontscale,
        )
    else:
        image = image.opts(
            width=width,
            height=height,
            title=title,
            fontscale=fontscale,
        )

    i = 0
    save = store_path + "_" + "{:03d}".format(i)
    while os.path.isfile(save + ".png"):
        i = i + 1
        save = store_path + "_" + "{:03d}".format(i)
    renderer.save(image, save)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Plot either mean squared deviation from perturbance over mean
        predicted deviation calculated via the sensitivity where the predicted
        axis is at most 1 such that plots for particle numbers are not entirely
        correct. Or plot the model state variable and predicted squared error
        over time.
        """
    )
    parser.add_argument(
        "--data_path",
        type=str,
        required=True,
        help="""
        Path to dataframe as NetCDF-file which had been created via create_mse.py.
        """,
    )
    parser.add_argument(
        "--add_zero_sens",
        action="store_true",
        help="""
        Add sensitivities of value zero to the far left and mark it with
        negative infinity.
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
        default=[],
        help="""
        Output parameter to plot for. Default plots all that are available
        in the dataframe.
        """,
    )
    parser.add_argument(
        "--backend",
        default="matplotlib",
        help="""
        Choose a backend for plotting. Options are:
        matplotlib: Most plots should be fine with it.
        bokeh: Recommended.
        """,
    )
    parser.add_argument(
        "--store_path",
        default="../pics/correlation",
        type=str,
        help="""
        Path to store the generated images.
        """,
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=0.90,
        help="""
        Plot a confidence ellipse around each sample with confidence
        between 0 and 1. If none is given, no ellipse will be plotted.
        """,
    )
    parser.add_argument(
        "--xlabel",
        default="Predicted Log MSD",
        type=str,
        help="""
        Alternative label for x-axis.
        """,
    )
    parser.add_argument(
        "--ylabel",
        default="True Log MSD",
        type=str,
        help="""
        Alternative label for y-axis. If plot_variant is "time_plot", then
        " [out_parameter]" is added.
        """,
    )
    parser.add_argument(
        "--title",
        default="True Deviation vs Prediction",
        type=str,
        help="""
        Title for the plot where " for [out_param]" is added.
        """,
    )
    parser.add_argument(
        "--width",
        default=900,
        type=int,
        help="""
        Width of plot in pixels.
        """,
    )
    parser.add_argument(
        "--height",
        default=900,
        type=int,
        help="""
        Height of plot in pixels.
        """,
    )
    parser.add_argument(
        "--set_zero",
        default=None,
        type=float,
        help="""
        If plot_variant is "correlation".
        Set any predicted squared errors with this value or lower to zero.
        This makes the plots easier to look at, when only a single parameter
        has a predicted error of 1e-200 or less.
        """,
    )
    parser.add_argument(
        "--plot_variant",
        default="correlation",
        type=str,
        help="""
        Plot either correlation plots with true deviation over predicted
        deviation by perturbing a parameter with "correlation",
        "correlation_hist" to add histograms on each axis,
        "histogram" for
        plotting the histogram of true and predicted deviations or
        use "time_plot" to plot (a single or mean) trajectory and
        the predicted deviation over time with the actual model
        state variable.
        """,
    )
    parser.add_argument(
        "--traj",
        type=int,
        default=-1,
        help="""
        If plot_type is "time_plot", the trajectory with this index will
        be plotted. If a value below zero is given, plot the mean of all.
        """,
    )
    parser.add_argument(
        "--in_parameter",
        type=str,
        nargs="+",
        default=[],
        help="""
        If plot_type is "time_plot", then plot the predicted deviation
        for those model parameters. If none are given, plot the top ten
        most influential parameters for each model state parameter.
        This plots all those predictions in one plot.
        """,
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help="""
        If plot_type is "time_plot", plot the y-axis as log10.
        """,
    )
    parser.add_argument(
        "--twinlabel",
        default="Predicted Squared Error",
        type=str,
        help="""
        Only if plot_type is "time_plot". Label for the twin axis.
        """,
    )
    parser.add_argument(
        "--logtwin",
        action="store_true",
        help="""
        If plot_type is "time_plot", plot the twin-axis as log10.
        """,
    )
    parser.add_argument(
        "--n_model_params",
        type=int,
        default=5,
        help="""
        If plot_type is "time_plot", plot this many model parameters.
        """,
    )
    parser.add_argument(
        "--min_time",
        type=float,
        default=None,
        help="""
        If plot_type is "time_plot", use this as start point for the plot
        as in time after ascent.
        """,
    )
    parser.add_argument(
        "--max_time",
        type=float,
        default=None,
        help="""
        If plot_type is "time_plot", use this as last point for the plot
        as in time after ascent.
        """,
    )
    args = parser.parse_args()
    ds = xr.open_dataset(args.data_path, decode_times=False, engine="netcdf4")
    if len(args.out_parameter) == 0:
        out_params = []
        for out_p in ds["Output Parameter"]:
            out_params.append(out_p.item())
    else:
        out_params = args.out_parameter

    param_title_names = {
        "QV": "Water Vapor",
        "QC": "Cloud",
        "QR": "Rain",
        "QS": "Snow",
        "QI": "Ice",
        "QG": "Graupel",
        "QH": "Hail",
        "NCCLOUD": "Cloud",
        "NCRAIN": "Rain",
        "NCSNOW": "Snow",
        "NCICE": "Ice",
        "NCGRAUPEL": "Graupel",
        "NCHAIL": "Hail",
    }

    if "correlation" in args.plot_variant or args.plot_variant == "histogram":
        for out_p in out_params:
            df = (
                ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
                .to_dataframe()
                .reset_index()
            )
            df = df.loc[df["Output Parameter"] == out_p]
            if (
                np.min(df["Predicted Squared Error"]) == 0
                and np.max(df["Predicted Squared Error"]) == 0
            ):
                # nothing to see here; Usually applies to
                # variables with no contents
                continue
            print(f"Plotting for {out_p}")
            if args.set_zero is not None:
                print("Replacing following parameters and values with zero:")
                print(
                    df.loc[
                        (df["Predicted Squared Error"] <= args.set_zero)
                        & (df["Predicted Squared Error"] != 0)
                    ][
                        [
                            "Input Parameter",
                            "Predicted Squared Error",
                            "Mean Squared Error",
                        ]
                    ]
                )
                df["Predicted Squared Error"] = df["Predicted Squared Error"].where(
                    df["Predicted Squared Error"] > args.set_zero, 0.0
                )
            if args.plot_variant == "histogram":
                plot_histogram(
                    df=df,
                    out_params=[out_p],
                    store_path=args.store_path,
                    backend=args.backend,
                    plot_types=args.plot_types,
                    add_zero_sens=args.add_zero_sens,
                    title=args.title + " for " + param_title_names[out_p],
                    xlabel=args.xlabel,
                    xlabel2=args.ylabel,
                    width=args.width,
                    height=args.height,
                )
            else:
                hist = False
                if "correlation_hist" == args.plot_variant:
                    hist = True
                plot_mse(
                    df=df,
                    out_params=[out_p],
                    store_path=args.store_path,
                    backend=args.backend,
                    plot_types=args.plot_types,
                    add_zero_sens=args.add_zero_sens,
                    confidence=args.confidence,
                    title=args.title + " for " + param_title_names[out_p],
                    xlabel=args.xlabel,
                    ylabel=args.ylabel,
                    width=args.width,
                    height=args.height,
                    hist=hist,
                )
    elif args.plot_variant == "time_plot":
        if args.traj < 0:
            df = ds.mean(dim=["trajectory"], skipna=True).to_dataframe().reset_index()
            df_mean = (
                ds.mean(dim=["time_after_ascent", "trajectory"], skipna=True)
                .to_dataframe()
                .reset_index()
            )
        else:
            df_mean = (
                ds.mean(dim=["time_after_ascent"], skipna=True)
                .to_dataframe()
                .reset_index()
            )
            df = ds.to_dataframe().reset_index()
            df = df.loc[df["trajectory"] == args.traj]
            df_mean = df_mean.loc[df_mean["trajectory"] == args.traj]
        for out_p in out_params:
            df_tmp = df.loc[df["Output Parameter"] == out_p]
            df_mean_tmp = df_mean.loc[df_mean["Output Parameter"] == out_p]
            if len(args.in_parameter) == 0:
                in_params = np.unique(
                    df_mean_tmp.nlargest(
                        args.n_model_params, "Predicted Squared Error"
                    )["Input Parameter"]
                )
            else:
                in_params = args.in_parameter
            df_tmp = df_tmp.loc[df_tmp["Input Parameter"].isin(in_params)]
            if (
                np.min(df_tmp["Predicted Squared Error"]) == 0
                and np.max(df_tmp["Predicted Squared Error"]) == 0
            ):
                # nothing to see here
                continue
            print(f"Plotting for {out_p}")
            plot_time_evolution(
                df=df_tmp,
                backend=args.backend,
                store_path=args.store_path,
                title=args.title + " for " + param_title_names[out_p],
                xlabel=args.xlabel,
                ylabel=args.ylabel + " " + parse_word(out_p),
                twinlabel=args.twinlabel,
                logy=args.logy,
                width=args.width,
                height=args.height,
                logtwin=args.logtwin,
                min_x=args.min_time,
                max_x=args.max_time,
            )
    else:
        print(f"plot_variant '{args.plot_variant}': No such plot variant. ABORTING!")
