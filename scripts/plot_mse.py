import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from bokeh.models import Range1d, LinearAxis, GlyphRenderer
import holoviews as hv
from holoviews import opts
import hvplot.pandas
import matplotlib
import numpy as np
import os
import pandas as pd
import panel as pn
import seaborn as sns
import sys
import xarray as xr

try:
    from Deriv_dask import Deriv_dask
    from latexify import (
        param_id_map,
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
        param_id_map,
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


def plot_errors(
    df,
    out_param,
    in_params,
    x_key,
    y_key,
    alpha=0.5,
    plot_types=True,
    add_zero_sens=False,
    n_std=2,
    linewidth=2,
    title=None,
    xlabel="Sensitivity",
    ylabel=None,
    width=12,
    height=12,
    log_x=False,
    log_y=False,
    corr_line=False,
    font_scale=None,
    save=True,
    latex=False,
    filename=None,
    s=2,
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
    out_param : string
       The output parameter in the given dataframe to plot for.
    store_path : string
        Path to folder where to store the image.
    alpha : float
        Alpha value for the dots.
    plot_types : bool
        Wether to plot the errors grouped by the type of the model parameters.
    add_zero_sens : bool
        Add sensitivities of value zero to the far left and mark it with
        negative infinity for the x-axis.
    n_std : float or None
        Plot a confidence ellipse around each sample with confidence
        as number of standard deviations. If none is given, no ellipse will be plotted.
    title : string
        Title for the plot.
    xlabel : string
        Alternative label for x-axis.
    ylabel : string
        Alternative label for y-axis.
    width : int
        Width of plot in inches.
    height : int
        Height of plot in inches.
    plot_kind : string
        "paper" for single plots, "single_plot" for a plot with
        multiple output parameters at once.
    legend_pos : string
        if plot_kind == "paper", then define the legend position here.
    corr_line : bool
        Plot a dashed line to show the 1-to-1 mapping in the plot.
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    latex : bool
        Use latex font.
    """

    def confidence_ellipse(x, y, ax, n_std=3.0, facecolor="none", **kwargs):
        """
        Create a plot of the covariance confidence ellipse of *x* and *y*.
        By https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html

        Parameters
        ----------
        x, y : array-like, shape (n, )
            Input data.

        ax : matplotlib.axes.Axes
            The axes object to draw the ellipse into.

        n_std : float
            The number of standard deviations to determine the ellipse's radiuses.

        **kwargs
            Forwarded to `~matplotlib.patches.Ellipse`

        Returns
        -------
        matplotlib.patches.Ellipse
        """
        if x.size != y.size:
            raise ValueError("x and y must be the same size")

        cov = np.cov(x, y)
        pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
        # Using a special case to obtain the eigenvalues of this
        # two-dimensional dataset.
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse(
            (0, 0),
            width=ell_radius_x * 2,
            height=ell_radius_y * 2,
            facecolor=facecolor,
            **kwargs,
        )

        # Calculating the standard deviation of x from
        # the squareroot of the variance and multiplying
        # with the given number of standard deviations.
        scale_x = np.sqrt(cov[0, 0]) * n_std
        mean_x = np.mean(x)

        # calculating the standard deviation of y ...
        scale_y = np.sqrt(cov[1, 1]) * n_std
        mean_y = np.mean(y)

        transf = (
            transforms.Affine2D()
            .rotate_deg(45)
            .scale(scale_x, scale_y)
            .translate(mean_x, mean_y)
        )

        ellipse.set_transform(transf + ax.transData)
        return ax.add_patch(ellipse)

    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})

    df_tmp = df.copy()
    df_tmp = df_tmp.loc[df_tmp["Output Parameter"] == out_param]
    df_tmp = df_tmp.loc[df_tmp["Input Parameter"].isin(in_params)]
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
        group_colors = {
            "artificial": "tab:orange",
            "artificial (threshold)": "tab:green",
            "physical": "tab:blue",
            "physical (high variability)": "tab:red",
            "1-moment": "k",
        }

    if add_zero_sens:
        inf_val = np.min(df_tmp.loc[df_tmp[x_key] != 0][x_key]) / 10
        df_tmp[x_key] = df_tmp[x_key].replace({0: inf_val})
    else:
        df_tmp = df_tmp.loc[df_tmp[x_key] != 0]
        inf_val = None

    if log_x:
        df_tmp[x_key] = np.log10(np.abs(df_tmp[x_key]))
    if log_y:
        df_tmp[y_key] = np.log10(np.abs(df_tmp[y_key]))

    fig = Figure()
    ax = fig.subplots()
    if plot_types:
        g = sns.scatterplot(
            data=df_tmp,
            x=x_key,
            y=y_key,
            hue="Group",
            ax=ax,
            palette=group_colors,
            alpha=alpha,
            s=s,
        )
        if n_std is not None and n_std > 0:
            for group in np.unique(df_tmp["Group"]):
                df_tmp2 = df_tmp.loc[df_tmp["Group"] == group]
                confidence_ellipse(
                    x=df_tmp2[x_key],
                    y=df_tmp2[y_key],
                    ax=ax,
                    n_std=n_std,
                    edgecolor=group_colors[group],
                    linewidth=linewidth,
                )
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, labels, fontsize=int(9 * font_scale))
    else:
        g = sns.scatterplot(
            data=df_tmp,
            x=x_key,
            y=y_key,
            ax=ax,
            alpha=alpha,
            s=s,
        )
        if n_std is not None and n_std > 0:
            confidence_ellipse(
                x=df_tmp[x_key],
                y=df_tmp[y_key],
                ax=ax,
                n_std=n_std,
                linewidth=linewidth,
                edgecolor="k",
            )
    if corr_line:
        tmp_min = np.nanmin(df_tmp[x_key])
        if tmp_min < np.nanmin(df_tmp[y_key]):
            tmp_min = np.nanmin(df_tmp[y_key])
        tmp_max = np.nanmax(df_tmp[x_key])
        if tmp_max > np.nanmax(df_tmp[y_key]):
            tmp_max = np.nanmax(df_tmp[y_key])
        line = matplotlib.lines.Line2D(
            [tmp_min, tmp_max],
            [tmp_min, tmp_max],
            lw=linewidth,
            color="k",
            linestyle="dashed",
        )
        ax.add_line(line)

    if font_scale is None:
        _ = ax.set_title(title)
    else:
        ax.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
    ax.set_xlabel(xlabel, fontsize=int(11 * font_scale))
    ax.set_ylabel(ylabel, fontsize=int(11 * font_scale))

    plt.tight_layout()
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
    plot_kind="paper",
    legend_pos="top_left",
    corr_line=False,
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
    out_params : list of string or string
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
    plot_kind : string
        "paper" for single plots, "single_plot" for a plot with
        multiple output parameters at once.
    legend_pos : string
        if plot_kind == "paper", then define the legend position here.
    corr_line : bool
        Plot a dashed line to show the 1-to-1 mapping in the plot.
    """

    in_params = np.unique(df["Input Parameter"])

    datashade = False
    if isinstance(out_params, str):
        out_params = [out_params]

    if plot_kind == "paper" or plot_kind == "single_plot":
        alpha = 0.5
    else:
        alpha = 1
    s = 12
    f_limits = (-2, 2)

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
        xlabel=None,
        ylabel=None,
        plot_path=store_path,
        inf_val=inf_val,
        legend_pos=legend_pos,
        corr_line=corr_line,
        plot_types=plot_types,
    )


def plot_errors_interactive(
    df,
):
    """
    Use plot_errors() interactively in a jupyter notebook for perturbation simulations.

    Parameters
    ----------
    df : pandas.DataFrame
        A dataframe with columns "Output Parameter" for the model state
        "Input Parameter" for the perturbed model parameter,
        "Predicted Squared Error" and "Predicted Error" for the sensitivity calculated to deviations
        in the next timestep, "Mean Squared Error" or "Mean Error" for the actual deviations

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    out_param_list = list(np.unique(df["Output Parameter"]))
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_param_list[0],
        options=out_param_list,
    )
    in_params_list = list(np.unique(df["Input Parameter"]))
    in_params = pn.widgets.CrossSelector(
        name="Parameter",
        value=in_params_list[0:2],
        options=in_params_list,
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
    alpha_slider = pn.widgets.FloatSlider(
        name="Alpha",
        start=0.1,
        end=1,
        step=0.1,
        value=0.6,
    )
    dot_slider = pn.widgets.IntSlider(
        name="Change the dot size",
        start=1,
        end=200,
        step=2,
        value=12,
    )
    logx_plot = pn.widgets.Toggle(
        name="Use log x-axis",
        value=False,
        button_type="success",
    )
    logy_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=False,
        button_type="success",
    )
    x_widget = pn.widgets.TextInput(
        name="X-label",
        placeholder="Predicted Error",
        value="Predicted Error",
    )
    y_widget = pn.widgets.TextInput(
        name="Y-label",
        placeholder="Ensemble Error",
        value="Ensemble Error",
    )
    data_variants = []
    for col in df:
        if col != "Output Parameter" and col != "Input Parameter":
            data_variants.append(col)
    x_data = pn.widgets.Select(
        name="X-axis",
        value=data_variants[0],
        options=data_variants,
    )
    y_data = pn.widgets.Select(
        name="Y-axis",
        value=data_variants[2],
        options=data_variants,
    )
    corr_line = pn.widgets.Toggle(
        name="Correlation line",
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
    group_toggle = pn.widgets.Toggle(
        name="Group parameters",
        value=False,
        button_type="success",
    )
    ellipsis_widget = pn.widgets.FloatSlider(
        name="Ellipsis in standard deviations",
        start=0,
        end=5,
        step=0.2,
        value=2,
    )

    plot_pane = pn.panel(
        pn.bind(
            plot_errors,
            df=df,
            out_param=out_param,
            in_params=in_params,
            x_key=x_data,
            y_key=y_data,
            alpha=alpha_slider,
            plot_types=group_toggle,
            n_std=ellipsis_widget,
            linewidth=line_slider,
            title=title_widget,
            xlabel=x_widget,
            ylabel=y_widget,
            width=width_slider,
            height=height_slider,
            log_x=logx_plot,
            log_y=logy_plot,
            corr_line=corr_line,
            font_scale=font_slider,
            save=save_button,
            latex=latex_button,
            filename=save_to_field,
            s=dot_slider,
        ),
    ).servable()
    return pn.Column(
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
        ),
        pn.Row(
            x_data,
            y_data,
            ellipsis_widget,
        ),
        pn.Row(
            save_to_field,
            save_button,
            latex_button,
        ),
        pn.Row(
            logx_plot,
            logy_plot,
            corr_line,
        ),
        pn.Row(
            in_params,
            pn.Column(
                out_param,
                x_widget,
                y_widget,
                group_toggle,
            ),
        ),
        pn.Row(
            dot_slider,
            line_slider,
            alpha_slider,
        ),
        title_widget,
        plot_pane,
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
    x_limits=None,
    plot_deviation=False,
    perturb_delta=None,
    alpha=1,
    width=900,
    height=900,
    logy=False,
    logtwin=False,
    save=False,
    font_scale=None,
    s=2,
    latex=False,
    trajectory=None,
    out_param=None,
    in_params=None,
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
        "Predicted Error" for the sensitivity calculated to deviations
        in the next timestep, (optional; if plot_deviation is true) "Mean Squared Error" for the actual deviations
        or (optional; if plot_deviation is false) "Not Perturbed Value".
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
    x_limits : Tuple of floats
        Set the limits of the x-axis. Values range between zero and one as percentages of the full
        range. Is used for interactive plotting. Is always superceeded by min_x or max_x.
    plot_deviation : bool
        If true, plot the deviation from perturbing parameters. If false,
        plot the actual model state variable.
    perturb_delta : float
        Time difference between different ensembles where perturbing started.
        This will plot vertical lines at each of those time steps.
    alpha : float
        Alpha value for the dots.
    width : int or float
        Width of plot in pixels or inches. Defaults to inches for values lower than 100 or floats.
    height : int or float
        Height of plot in pixels or inches. Defaults to inches for values lower than 100 or floats.
    logy : bool
        If true, use log10 on y-axis
    logtwin : bool
        If true, use log10 on twin-axis
    save : bool
        If true, save the plot to the given path.
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    s : int
        Dot size.
    latex : bool
        Use latex font.
    trajectory : int
        In case df has multiple trajectories, you may select one.
    out_param : string
        Name of the model state variable for which sensitivities shall be plotted in case multiple
        are available in df.
    in_params : list of strings
        Name of the model parameters in df to plot.

    Returns
    -------
    matplotlib.figure.Figure with the plot drawn onto it.
    """
    import matplotlib.ticker as tick

    sns.set(rc={"text.usetex": latex, "axes.grid": True})

    if trajectory is not None:
        df = df.loc[df["trajectory"] == trajectory]
    if out_param is not None:
        df = df.loc[df["Output Parameter"] == out_p]
    if in_params is not None:
        df = df.loc[df["Input Parameter"].isin(in_params)]

    if backend == "matplotlib":
        if font_scale is None:
            if isinstance(height, float) or height < 100:
                font_scale = height / 3 * 0.8
            else:
                font_scale = height / 900 * 0.8
        if s is None:
            s = 2
        # matplotlib.rcParams["font.family"] = "sans-serif"
        # matplotlib.rcParams["font.sans-serif"] = "Helvetica"
        title_font_size = 12  # * font_scale
        axis_label_text_font_size = 10  # * font_scale
        major_label_text_font_size = 10  # * font_scale
        aspect = width / height
    else:
        if font_scale is None:
            if isinstance(height, float) or height < 100:
                font_scale = height * 1.2
            else:
                font_scale = height / 350
        axis_label_text_font_size = int(13 * font_scale * 0.8)
        major_label_text_font_size = int(11 * font_scale * 0.8)
        if s is None:
            s = 14

    hv.extension(backend)

    tmp_max = np.nanmax(df["time_after_ascent"])
    tmp_min = np.nanmin(df["time_after_ascent"])
    if min_x is not None:
        df = df.loc[df["time_after_ascent"] >= min_x]
    elif x_limits is not None:
        min_x = tmp_min + (tmp_max - tmp_min) * x_limits[0]
        df = df.loc[df["time_after_ascent"] >= min_x]
    if max_x is not None:
        df = df.loc[df["time_after_ascent"] <= max_x]
    elif x_limits is not None:
        max_x = tmp_min + (tmp_max - tmp_min) * x_limits[1]
        df = df.loc[df["time_after_ascent"] <= max_x]

    # Make time after ascent to minutes
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
        df["Predicted Error"] = np.log10(np.abs(df["Predicted Error"]))
        min_log_twin = (
            np.nanmin(df.loc[df["Predicted Error"] != -np.inf]["Predicted Error"]) - 1
        )
        df.replace(-np.inf, min_log_twin, inplace=True)
    # df = df.loc[df["Predicted Error"] < np.max(df["Predicted Error"])]
    lower_y = np.min(df["Predicted Error"])
    upper_y = np.max(df["Predicted Error"])

    delta = (upper_y - lower_y) / 12
    # if upper_y == 0:
    #     upper_y += 4 * delta
    #     delta = (upper_y - lower_y) / 12
    # else:
    upper_y += delta
    lower_y -= delta

    lower_y2 = np.min(df[y])
    upper_y2 = np.max(df[y])
    delta2 = (upper_y2 - lower_y2) / 12
    lower_y2 -= delta2
    upper_y2 += delta2

    if "." not in store_path:
        filetype = "png"
        tmp_path = store_path
    else:
        filetype = store_path.split(".")[-1]
        tmp_path = store_path.split(".")[0]
    renderer = hv.Store.renderers[backend].instance(fig=filetype, dpi=300)

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
            ax.set_ylim((lower_y2, upper_y2))
            ax.set_title(
                title,
                loc="left",
                # fontdic={},
                **{
                    "fontweight": "bold",
                    "fontsize": title_font_size * font_scale,
                },
            )
            ax.set_ylabel(
                ylabel,
                **{
                    "rotation": 90,
                    # "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size * font_scale,
                },
            )
            ax.set_xlabel(
                xlabel,
                **{
                    # "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size
                    * font_scale,
                },
            )
            # ax.minorticks_on()
            ax.xaxis.set_ticks_position("none")

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
                axis="both",
                which="major",
                labelsize=major_label_text_font_size * font_scale,
            )
            ax.xaxis.set_ticks_position("none")
            ax.set_yticks(np.linspace(ax.get_ybound()[0], ax.get_ybound()[1], 5))
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
                # yticks=np.arange(lower_y2 + delta, upper_y2 - delta2/2, delta2*2.5),
            )
            .opts(
                opts.Scatter(
                    s=s,
                )
            )
        )

        def twinx_per_timestep(plot, element):
            ax = plot.handles["axis"]
            ax.xaxis.set_ticks_position("none")
            twinax = ax.twinx()

            twinax.set_ylim((lower_y, upper_y))
            twinax.set_ylabel(
                twinlabel,
                labelpad=13.0 * font_scale,
                **{
                    "rotation": 270,
                    # "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size * font_scale,
                },
            )
            if logtwin:
                yticks_tmp = twinax.get_yticks()
                yticks_tmp[0] = min_log_twin
                if yticks_tmp[1] - yticks_tmp[0] < 3.5:
                    yticks_tmp = np.delete(yticks_tmp, 1)
                twinax.set_yticks(yticks_tmp)
            # else:
            # twinax.set_yticks(np.arange(lower_y + delta, upper_y, delta*2.5))
            # else:
            #     twinax.minorticks_on()

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
                axis="both",
                which="major",
                labelsize=major_label_text_font_size * font_scale,
            )
            twinax.xaxis.set_ticks_position("none")
            twinax.set_yticks(
                np.linspace(twinax.get_ybound()[0], twinax.get_ybound()[1], 5)
            )
            plot.handles["axis"] = twinax

        def twinx2(plot, element):
            ax = plot.handles["axis"]
            twinax = ax.twinx()
            twinax.set_ylim((lower_y, upper_y))
            twinax.set_yticks([])
            plot.handles["axis"] = twinax
            twinax.xaxis.set_ticks_position("none")

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
                ).opts(
                    initial_hooks=[twinx2],
                    apply_ranges=False,
                    # yticks=np.arange(lower_y + delta, upper_y - delta/2, delta*2.5),
                ).opts(
                    opts.Scatter(
                        s=s,
                        cmap=cmap_map,
                    ),
                )

        legend_overlay = hv.NdOverlay(
            {
                data_types[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                    opts.Scatter(
                        s=s * 4,
                        color=cmap_values[i],
                    )
                )
                for i in range(len(data_types))
            }
        ).opts(
            fontsize={
                "legend": major_label_text_font_size * font_scale,
                "xticks": major_label_text_font_size * font_scale,
                "yticks": major_label_text_font_size * font_scale,
                "title": title_font_size * font_scale,
                "xlabel": axis_label_text_font_size * font_scale,
                "ylabel": axis_label_text_font_size * font_scale,
            }
        )
        if perturb_delta is not None:
            image = left * twin * perturb_lines * legend_overlay
        else:
            image = legend_overlay * left * twin
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
        if isinstance(width, float) or width < 100:
            image = image.opts(
                fig_inches=(width, height),
                # fontscale=font_scale,
            )
        else:
            image = image.opts(
                fig_inches=(width / 300.0, height / 300.0),
                # fontscale=font_scale,
            )
    else:
        if isinstance(width, float) or width < 100:
            image = image.opts(
                width=int(width * 300),
                height=int(height * 300),
                title=title,
                # fontscale=font_scale,
            )
        else:
            image = image.opts(
                width=width,
                height=height,
                title=title,
                # fontscale=font_scale,
            )

    if save:
        i = 0
        save_path = tmp_path + "_" + "{:03d}".format(i)
        while os.path.isfile(save_path + "." + filetype):
            i = i + 1
            save_path = tmp_path + "_" + "{:03d}".format(i)
        renderer.save(image, save_path)
    return image


def plot_time_evolution_interactive(ds, perturbed=False):
    """
    Use plot_time_evolution() interactively in a jupyter notebook for perturbation or sensitivity simulations.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with trajectories from a perturbed or sensitivity simulation.
    perturbed : bool
        If true, then ds is from a perturbed ensemble simulation. Otherwise a sensitivity simulation is
        assumed.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    if "Output Parameter" in ds:
        out_param_list = ds["Output Parameter"].values.tolist()
    else:
        out_param_list = []
        for out_p in ds["Output_Parameter_ID"]:
            out_param_list.append(param_id_map[out_p.item()])
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_param_list[0],
        options=out_param_list,
    )
    in_params_list = []
    for col in ds:
        if col[0] == "d" and col != "deposition":
            in_params_list.append(col)
    in_params_list = list(np.sort(in_params_list))
    in_params = pn.widgets.CrossSelector(
        name="Parameter",
        value=in_params_list[0:2],
        options=in_params_list,
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
    dot_slider = pn.widgets.IntSlider(
        name="Change the dot size",
        start=1,
        end=200,
        step=2,
        value=12,
    )
    x_slider = pn.widgets.RangeSlider(
        name="X Limits in percent",
        start=0,
        end=1,
        value=(0, 1),
        step=0.001,
    )
    log_twin_plot = pn.widgets.Toggle(
        name="Use log twin y-axis",
        value=False,
        button_type="success",
    )
    log_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=False,
        button_type="success",
    )
    x_widget = pn.widgets.TextInput(
        name="X-label",
        placeholder="Time after ascent [min]",
        value="Time after ascent [min]",
    )
    y_widget = pn.widgets.TextInput(
        name="Y-label",
        placeholder="Specific humidty [kg/kg]",
    )
    y2_widget = pn.widgets.TextInput(
        name="Twin y-label",
        placeholder="Predicted Deviation [kg/kg]",
        value="Predicted Deviation [kg/kg]",
    )
    traj_widget = pn.widgets.IntSlider(
        name="Trajectory",
        start=0,
        end=len(ds["trajectory"]) - 1,
        step=1,
    )

    def prepare_dataset_and_plot(
        ds,
        backend,
        store_path,
        title,
        xlabel,
        ylabel,
        twinlabel,
        logtwin,
        logy,
        width,
        height,
        x_limits,
        save,
        latex,
        font_scale,
        s,
        trajectory,
        out_param,
        in_params,
        plot_deviation,
    ):
        if not perturbed:
            out_p_id = np.argwhere(np.asarray(param_id_map) == out_param).item()
            ds_traj = ds.isel({"ensemble": 0, "trajectory": trajectory}).sel(
                {"Output_Parameter_ID": out_p_id}
            )
            # Sensitivity simulation
            predicted_errors = np.array([])
            in_param_thingy = np.array([])
            unperturbed = np.array([])
            time = np.array([])
            for col in in_params:
                predicted_errors = np.append(predicted_errors, ds_traj[col].values)
                in_param_thingy = np.append(
                    in_param_thingy, np.repeat(col, len(ds_traj[col].values))
                )
                unperturbed = np.append(unperturbed, ds_traj[out_param].values)
                time = np.append(time, ds_traj["time_after_ascent"].values)
            df = pd.DataFrame(
                data={
                    "Predicted Error": predicted_errors,
                    "Predicted Squared Error": predicted_errors ** 2,
                    "Output Parameter": np.repeat(out_param, len(predicted_errors)),
                    "Input Parameter": in_param_thingy,
                    "Not Perturbed Value": unperturbed,
                    "time_after_ascent": time,
                }
            )
            return plot_time_evolution(
                df=df,
                backend=backend,
                store_path=store_path,
                title=title,
                xlabel=xlabel,
                ylabel=ylabel,
                twinlabel=twinlabel,
                logtwin=logtwin,
                logy=logy,
                width=width,
                height=height,
                x_limits=x_limits,
                save=save,
                latex=latex,
                font_scale=font_scale,
                s=s,
                plot_deviation=plot_deviation,
            )
        else:
            # Perturbed ensemble
            df = ds.to_dataframe().reset_index()
            return plot_time_evolution(
                df=df,
                backend=backend,
                store_path=store_path,
                title=title,
                xlabel=xlabel,
                ylabel=ylabel,
                twinlabel=twinlabel,
                logtwin=logtwin,
                logy=logy,
                width=width,
                height=height,
                x_limits=x_limits,
                save=save,
                latex=latex,
                s=s,
                trajectory=trajectory,
                font_scale=font_scale,
                out_param=out_param,
                in_params=in_params,
                plot_deviation=plot_deviation,
            )

    plot_pane = pn.panel(
        pn.bind(
            prepare_dataset_and_plot,
            ds=ds,
            backend="matplotlib",
            store_path=save_to_field,
            title=title_widget,
            xlabel=x_widget,
            ylabel=y_widget,
            twinlabel=y2_widget,
            logtwin=log_twin_plot,
            logy=log_plot,
            width=width_slider,
            height=height_slider,
            x_limits=x_slider,
            save=save_button,
            latex=latex_button,
            s=dot_slider,
            trajectory=traj_widget,
            out_param=out_param,
            in_params=in_params,
            plot_deviation=perturbed,
            font_scale=font_slider,
        ),
    ).servable()

    return pn.Column(
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
            log_plot,
            log_twin_plot,
            traj_widget,
        ),
        pn.Row(
            in_params,
            pn.Column(
                out_param,
                x_widget,
                y_widget,
                y2_widget,
            ),
        ),
        pn.Row(
            x_slider,
            dot_slider,
            title_widget,
        ),
        plot_pane,
    )


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Plot either mean squared deviation/error from perturbation over mean
            predicted deviation/error calculated via the sensitivity where the predicted
            axis is at most 1 such that plots for particle numbers are not entirely
            correct. Or plot the model state variable and predicted squared error
            over time.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        type=str,
        required=True,
        help=textwrap.dedent(
            """\
            Path to dataframe as NetCDF-file which had been created via create_mse.py.
            """
        ),
    )
    parser.add_argument(
        "--add_zero_sens",
        action="store_true",
        help=textwrap.dedent(
            """\
            Add sensitivities of value zero to the far left and mark it with
            negative infinity.
            """
        ),
    )
    parser.add_argument(
        "--plot_types",
        action="store_true",
        help=textwrap.dedent(
            """\
            If true: Plot input parameter types in different colors.
            Types means here if a parameter is physical or rather artificial.
            """
        ),
    )
    parser.add_argument(
        "--out_parameter",
        type=str,
        nargs="+",
        default=[],
        help=textwrap.dedent(
            """\
            Output parameter to plot for. Default plots all that are available
            in the dataframe in a separate plot. If you want a plot with all
            datapoints in one plot, use "all_at_once".
            """
        ),
    )
    parser.add_argument(
        "--backend",
        default="matplotlib",
        help=textwrap.dedent(
            """\
            Choose a backend for plotting. Options are:
            matplotlib: Most plots should be fine with it.
            bokeh: Recommended.
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        default="../pics/correlation",
        type=str,
        help=textwrap.dedent(
            """\
            Path to store the generated images.
            """
        ),
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=None,
        help=textwrap.dedent(
            """\
            Plot a confidence ellipse around each sample with confidence
            between 0 and 1. If none is given, no ellipse will be plotted.
            """
        ),
    )
    parser.add_argument(
        "--xlabel",
        default="Predicted Log MSD",
        type=str,
        help=textwrap.dedent(
            """\
            Alternative label for x-axis.
            """
        ),
    )
    parser.add_argument(
        "--ylabel",
        default="True Log MSD",
        type=str,
        help=textwrap.dedent(
            """\
            Alternative label for y-axis. If plot_variant is "time_plot", then
            " [out_parameter]" is added.
            """
        ),
    )
    parser.add_argument(
        "--title",
        default="True Deviation vs Prediction",
        type=str,
        help=textwrap.dedent(
            """\
            Title for the plot where " for [out_param]" is added.
            """
        ),
    )
    parser.add_argument(
        "--width",
        default=900,
        type=int,
        help=textwrap.dedent(
            """\
            Width of plot in pixels.
            """
        ),
    )
    parser.add_argument(
        "--height",
        default=900,
        type=int,
        help=textwrap.dedent(
            """\
            Height of plot in pixels.
            """
        ),
    )
    parser.add_argument(
        "--set_zero",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            If plot_variant is "correlation".
            Set any predicted squared errors with this value or lower to zero.
            This makes the plots easier to look at, when only a single parameter
            has a predicted error of 1e-200 or less.
            """
        ),
    )
    parser.add_argument(
        "--plot_variant",
        default="correlation",
        type=str,
        help=textwrap.dedent(
            """\
            Plot either correlation plots with true deviation over predicted
            deviation by perturbing a parameter with "correlation",
            "correlation_hist" to add histograms on each axis,
            "histogram" for
            plotting the histogram of true and predicted deviations or
            use "time_plot" to plot (a single or mean) trajectory and
            the predicted deviation over time with the actual model
            state variable.
            """
        ),
    )
    parser.add_argument(
        "--traj",
        type=int,
        default=-1,
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", the trajectory with this index will
            be plotted. If a value below zero is given, plot the mean of all.
            """
        ),
    )
    parser.add_argument(
        "--in_parameter",
        type=str,
        nargs="+",
        default=[],
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", then plot the predicted deviation
            for those model parameters. If none are given, plot the top ten
            most influential parameters for each model state parameter.
            This plots all those predictions in one plot.
            """
        ),
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", plot the y-axis as log10.
            """
        ),
    )
    parser.add_argument(
        "--twinlabel",
        default="Predicted Squared Error",
        type=str,
        help=textwrap.dedent(
            """\
            Only if plot_type is "time_plot". Label for the twin axis.
            """
        ),
    )
    parser.add_argument(
        "--logtwin",
        action="store_true",
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", plot the twin-axis as log10.
            """
        ),
    )
    parser.add_argument(
        "--n_model_params",
        type=int,
        default=5,
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", plot this many model parameters.
            """
        ),
    )
    parser.add_argument(
        "--min_time",
        type=float,
        default=None,
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", use this as start point for the plot
            as in time after ascent.
            """
        ),
    )
    parser.add_argument(
        "--max_time",
        type=float,
        default=None,
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", use this as last point for the plot
            as in time after ascent.
            """
        ),
    )
    parser.add_argument(
        "--legend_pos",
        type=str,
        default="bottom_right",
        help=textwrap.dedent(
            """\
            Define the position of the legend for most plots.
            """
        ),
    )
    parser.add_argument(
        "--corr_line",
        action="store_true",
        help=textwrap.dedent(
            """\
            Add a dashed line for a 1-to-1 map of the data.
            """
        ),
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
        "QC": "Cloud Mass",
        "QR": "Rain Mass",
        "QS": "Snow Mass",
        "QI": "Ice Mass",
        "QG": "Graupel Mass",
        "QH": "Hail Mass",
        "NCCLOUD": "Cloud Number",
        "NCRAIN": "Rain Number",
        "NCSNOW": "Snow Number",
        "NCICE": "Ice Number",
        "NCGRAUPEL": "Graupel Number",
        "NCHAIL": "Hail Number",
        "pressure": "Pressure",
        "QR_OUT": "QR_OUT",
        "QS_OUT": "QS_OUT",
        "QI_OUT": "QI_OUT",
        "QG_OUT": "QG_OUT",
        "QH_OUT": "QH_OUT",
        "NR_OUT": "Precipitation of Rain Droplets",
        "NS_OUT": "NS_OUT",
        "NI_OUT": "NI_OUT",
        "NG_OUT": "NG_OUT",
        "NH_OUT": "NH_OUT",
        "latent_heat": "Latent Heating",
        "latent_cool": "Latent Cooling",
    }

    if "correlation" in args.plot_variant or args.plot_variant == "histogram":
        if out_params[0] == "all_at_once":
            all_df = None
            out_params = []
            for out_p in ds["Output Parameter"]:
                out_params.append(out_p.item())
            for out_p in out_params:
                df = (
                    ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
                    .to_dataframe()
                    .reset_index()
                )
                df = df.loc[df["Output Parameter"] == out_p]

                if len(args.in_parameter) == 0:
                    in_params = list(np.unique(df["Input Parameter"]))
                else:
                    in_params = args.in_parameter
                df = df.loc[df["Input Parameter"].isin(in_params)]

                if (
                    np.min(df["Predicted Squared Error"]) == 0
                    and np.max(df["Predicted Squared Error"]) == 0
                ):
                    # nothing to see here; Usually applies to
                    # variables with no contents
                    continue
                if all_df is None:
                    all_df = df
                else:
                    all_df = all_df.append(df)

            hist = False
            if "correlation_hist" == args.plot_variant:
                hist = True
            plot_mse(
                df=all_df,
                out_params=out_params,
                store_path=args.store_path,
                backend=args.backend,
                plot_types=args.plot_types,
                add_zero_sens=args.add_zero_sens,
                confidence=args.confidence,
                title=args.title,
                xlabel=args.xlabel,
                ylabel=args.ylabel,
                width=args.width,
                height=args.height,
                hist=hist,
                plot_kind="single_plot",
            )
        else:
            for out_p in out_params:
                df = (
                    ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
                    .to_dataframe()
                    .reset_index()
                )
                df = df.loc[df["Output Parameter"] == out_p]

                if len(args.in_parameter) == 0:
                    in_params = list(np.unique(df["Input Parameter"]))
                else:
                    in_params = args.in_parameter
                df = df.loc[df["Input Parameter"].isin(in_params)]

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
                        title=args.title + r" for " + param_title_names[out_p],
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
                        title=args.title + r" for " + param_title_names[out_p],
                        xlabel=args.xlabel,
                        ylabel=args.ylabel,
                        width=args.width,
                        height=args.height,
                        hist=hist,
                        legend_pos=args.legend_pos,
                        corr_line=args.corr_line,
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
            # Result from a perturbation
            if "time_after_ascent" in ds.dims:
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
                        in_params = list(
                            np.unique(
                                df_mean_tmp.nlargest(
                                    args.n_model_params, "Predicted Squared Error"
                                )["Input Parameter"]
                            )
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
                    print(df_tmp.columns)
                    df_tmp["Not Perturbed Value"] = df_tmp["Not Perturbed Value"] * -1
                    df_tmp["Predicted Error"] = df_tmp["Predicted Error"] * -1
                    plot_time_evolution(
                        df=df_tmp,
                        backend=args.backend,
                        store_path=args.store_path,
                        title=args.title + " for " + param_title_names[out_p],
                        xlabel=args.xlabel,
                        ylabel=args.ylabel
                        + " "
                        + parse_word(out_p)
                        + " "
                        + get_unit(out_p, brackets=True),
                        twinlabel=args.twinlabel + " " + get_unit(out_p, brackets=True),
                        logy=args.logy,
                        width=args.width,
                        height=args.height,
                        logtwin=args.logtwin,
                        min_x=args.min_time,
                        max_x=args.max_time,
                        save=True,
                    )
            else:
                # Result from a sensitivity simulation
                for out_p in out_params:
                    out_p_id = np.argwhere(np.asarray(param_id_map) == out_p).item()
                    ds_traj = ds.isel({"ensemble": 0, "trajectory": args.traj}).sel(
                        {"Output_Parameter_ID": out_p_id}
                    )
                    df = ds_traj.to_dataframe().reset_index()
                    if args.min_time is not None and args.max_time is not None:
                        ds_mean = np.abs(
                            ds_traj.where(
                                (ds_traj["time_after_ascent"] >= args.min_time)
                                & (ds_traj["time_after_ascent"] <= args.max_time)
                            )
                        ).mean(dim=["time"], skipna=True)
                    elif args.min_time is not None:
                        ds_mean = np.abs(
                            ds_traj.where(ds_traj["time_after_ascent"] >= args.min_time)
                        ).mean(dim=["time"], skipna=True)
                    elif args.max_time is not None:
                        ds_mean = np.abs(
                            ds_traj.where(ds_traj["time_after_ascent"] <= args.max_time)
                        ).mean(dim=["time"], skipna=True)
                    else:
                        ds_mean = np.abs(ds_traj).mean(dim=["time"], skipna=True)
                    in_params = []
                    max_vals = np.zeros(args.n_model_params)
                    for i in range(args.n_model_params):
                        for key in ds_mean:
                            if key[0] == "d" and key != "deposition":
                                if ds_mean[key].values.item() > max_vals[i]:
                                    if key not in in_params:
                                        if len(in_params) <= i:
                                            in_params.append(key)
                                        else:
                                            in_params[i] = key
                                        max_vals[i] = ds_mean[key].values.item()
                    predicted_errors = np.array([])
                    in_param_thingy = np.array([])
                    unperturbed = np.array([])
                    time = np.array([])
                    for col in in_params:
                        predicted_errors = np.append(
                            predicted_errors, ds_traj[col].values
                        )
                        in_param_thingy = np.append(
                            in_param_thingy, np.repeat(col, len(ds_traj[col].values))
                        )
                        unperturbed = np.append(unperturbed, ds_traj[out_p].values)
                        time = np.append(time, ds_traj["time_after_ascent"].values)
                    df_tmp = pd.DataFrame(
                        data={
                            "Predicted Error": predicted_errors,
                            "Predicted Squared Error": predicted_errors ** 2,
                            "Output Parameter": np.repeat(out_p, len(predicted_errors)),
                            "Input Parameter": in_param_thingy,
                            "Not Perturbed Value": unperturbed,
                            "time_after_ascent": time,
                        }
                    )

                    if (
                        np.min(df_tmp["Predicted Squared Error"]) == 0
                        and np.max(df_tmp["Predicted Squared Error"]) == 0
                    ):
                        # nothing to see here
                        continue
                    print(f"Plotting for {out_p}")
                    print(df_tmp.columns)
                    plot_time_evolution(
                        df=df_tmp,
                        backend=args.backend,
                        store_path=args.store_path,
                        title=args.title + " for " + param_title_names[out_p],
                        xlabel=args.xlabel,
                        ylabel=args.ylabel
                        + " "
                        + parse_word(out_p)
                        + " "
                        + get_unit(out_p, brackets=True),
                        twinlabel=args.twinlabel + " " + get_unit(out_p, brackets=True),
                        logy=args.logy,
                        width=args.width,
                        height=args.height,
                        logtwin=args.logtwin,
                        min_x=args.min_time,
                        max_x=args.max_time,
                    )
    else:
        print(f"plot_variant '{args.plot_variant}': No such plot variant. ABORTING!")
