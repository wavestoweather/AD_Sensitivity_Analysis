"""Plot predicted errors.

"""
import holoviews as hv
from holoviews import Curve as hvCurve
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


from ad_sensitivity_analysis.plot.aux_functions import save_plot, save_plot_renderer
from ad_sensitivity_analysis.plot.aux_functions_mse import (
    _add_histogram,
    _get_color_vals,
    _get_kwargs,
    _get_title_ops,
    _get_x_limits,
    _grid_paper_inf_plot,
    _grid_paper_plot,
    _plot_all_into_one_bokeh,
    _plot_all_into_one_mpl,
    _prepare_errors_plot,
    _prepare_mse_plot,
    _save_multiple_plots,
    _set_abs_log,
    _set_final_opts,
    _set_labels,
)
from ad_sensitivity_analysis.plot.ellipse_funcs import add_ellipse, confidence_ellipse


# pylint: disable=too-many-arguments, too-many-locals
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
    dot_size=2,
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
    dot_size : int
    """
    df_tmp, group_colors = _prepare_errors_plot(
        df,
        out_param,
        in_params,
        latex,
        plot_types,
        width,
        height,
        add_zero_sens,
        log_x,
        log_y,
        x_key,
        y_key,
    )

    fig = Figure()
    # pylint: disable=no-member
    ax = fig.subplots()
    if plot_types:
        sns.scatterplot(
            data=df_tmp,
            x=x_key,
            y=y_key,
            hue="Group",
            ax=ax,
            palette=group_colors,
            alpha=alpha,
            dot_size=dot_size,
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
        ax.legend(handles, labels, fontsize=int(9 * font_scale))
    else:
        sns.scatterplot(
            data=df_tmp,
            x=x_key,
            y=y_key,
            ax=ax,
            alpha=alpha,
            dot_size=dot_size,
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
        tmp_min = np.nanmax(np.nanmin(df_tmp[x_key]), np.nanmin(df_tmp[y_key]))
        tmp_max = np.nanmin(np.nanmax(df_tmp[x_key]), np.nanmax(df_tmp[y_key]))
        line = Line2D(
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
    fig = ax.get_figure()
    if filename is not None and save:
        save_plot(filename, fig)
    return fig


# pylint: disable=too-many-branches, too-many-statements
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
    log_x=False,
    log_y=False,
    log_func=np.log10,
    kind="paper",
    scatter_size=None,
    xticks=10,
    yticks=10,
    plot_singles=False,
    formatter_limits=None,
    linewidth=2,
):
    """
    Plot the dataframe which should hold parameters with their sensitivity
    to one model state parameter and the actual deviation when perturbing
    on a log-log plot.
    Backend "matplotlib" might result in wrong colors.

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
    (
        mse_df,
        inf_val,
        title,
        error_key,
        sens_key,
        legend,
        out_params,
        alpha,
    ) = _prepare_mse_plot(
        df=df,
        out_params=out_params,
        kind=kind,
        plot_kind=plot_kind,
        plot_types=plot_types,
        legend_pos=legend_pos,
        add_zero_sens=add_zero_sens,
        title=title,
        log_func=np.log10,
    )

    layout_kwargs, scatter_kwargs, hist_wh, aspect, fontscale = _get_kwargs(
        width=width,
        height=height,
        backend=backend,
        scatter_size=scatter_size,
        formatter_limits=formatter_limits,
    )
    _set_abs_log(df=mse_df, log=log_x, key=sens_key, x_axis=True, log_func=log_func)
    _set_abs_log(df=mse_df, log=log_y, key=error_key, x_axis=False, log_func=log_func)
    xlabel, ylabel = _set_labels(
        sens_key=sens_key,
        error_key=error_key,
        xlabel=xlabel,
        ylabel=ylabel,
        log_x=log_x,
        log_y=log_y,
    )
    title, opts_dic, opts_dic2 = _get_title_ops(
        title, backend, xlabel, ylabel, fontscale, aspect
    )
    if kind == "single_plot":
        if backend == "matplotlib":
            mse_plot = _plot_all_into_one_mpl(
                mse_df=mse_df,
                out_params=out_params,
                sens_key=sens_key,
                error_key=error_key,
                title=title,
                datashade=False,
                xticks=xticks,
                yticks=yticks,
                xlabel=xlabel,
                ylabel=ylabel,
                fontscale=fontscale,
                aspect=aspect,
                backend=backend,
                alpha=alpha,
                plot_types=plot_types,
            )
        else:
            mse_plot = _plot_all_into_one_bokeh(
                mse_df=mse_df,
                sens_key=sens_key,
                error_key=error_key,
                title=title,
                xticks=xticks,
                xlabel=xlabel,
                ylabel=ylabel,
                fontscale=fontscale,
                backend=backend,
                width=width,
                height=height,
                scatter_kwargs=scatter_kwargs,
                layout_kwargs=layout_kwargs,
                inf_val=inf_val,
                log_x=log_x,
                alpha=alpha,
                plot_types=plot_types,
            )
    if kind in ("grid_plot", "paper"):
        min_x, max_x = _get_x_limits(mse_df, sens_key, log_x, kind)
        cmap_values, colors, by_col = _get_color_vals(mse_df, backend, plot_types)
        all_plots = None
        if kind != "paper":
            opts_dic = {"fontscale": fontscale}
            if backend == "matplotlib":
                opts_dic["aspect"] = aspect
            all_plots = (
                mse_df.hvplot.hist(
                    y=sens_key,
                    by=by_col,
                    alpha=alpha,
                    legend=True,
                    grid=True,
                    title=title,
                    color=cmap_values,
                )
                .opts(**opts_dic)
                .options(xlabel="")
            )
        for out_param in out_params:
            if kind != "paper":
                title = None
            tmp_df = mse_df.loc[mse_df["Output Parameter"] == out_param]
            min_y = tmp_df[error_key].min()
            max_y = tmp_df[error_key].max()
            delta_y = (max_y - min_y) / 20
            min_y -= delta_y
            max_y += delta_y
            if kind == "paper" and inf_val is not None:
                mse_plot = _grid_paper_inf_plot(
                    tmp_df=tmp_df,
                    sens_key=sens_key,
                    error_key=error_key,
                    by_col=by_col,
                    width=width,
                    height=height,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    xticks=xticks,
                    scatter_kwargs=scatter_kwargs,
                    cmap_values=cmap_values,
                    min_x=min_x,
                    max_x=max_x,
                    min_y=min_y,
                    max_y=max_y,
                    opts_dic2=opts_dic2,
                    legend=legend,
                    title=title,
                    inf_val=inf_val,
                    alpha=alpha,
                )

            else:
                mse_plot = _grid_paper_plot(
                    tmp_df=tmp_df,
                    sens_key=sens_key,
                    error_key=error_key,
                    by_col=by_col,
                    width=width,
                    height=height,
                    xlabel=xlabel,
                    ylabel=ylabel,
                    scatter_kwargs=scatter_kwargs,
                    cmap_values=cmap_values,
                    min_x=min_x,
                    max_x=max_x,
                    min_y=min_y,
                    max_y=max_y,
                    opts_dic2=opts_dic2,
                    legend=legend,
                    title=title,
                    alpha=alpha,
                )
            if corr_line:
                tmp_min = min(min_x, min_y)
                tmp_max = min(max_x, max_y)
                line_plot = hvCurve([[tmp_min, tmp_min], [tmp_max, tmp_max]]).opts(
                    line_dash="dashed", color="black"
                )
                mse_plot = mse_plot * line_plot

            if out_param == out_params[-1]:
                mse_plot = mse_plot.options(xlabel=xlabel)

            if confidence is not None:
                mse_plot = add_ellipse(
                    tmp_df,
                    mse_plot,
                    sens_key,
                    error_key,
                    by_col,
                    confidence,
                    colors,
                    linewidth,
                    inf_val,
                )
            if hist:
                xdist = tmp_df.hvplot.hist(
                    y=sens_key,
                    by=by_col,
                    alpha=alpha,
                    legend=False,
                    color=cmap_values,
                ).opts(fontscale=fontscale)
                ydist = tmp_df.hvplot.hist(
                    y=error_key,
                    by=by_col,
                    alpha=alpha,
                    legend=False,
                    color=cmap_values,
                ).opts(fontscale=fontscale)
                mse_plot = _add_histogram(
                    mse_plot,
                    width,
                    height,
                    scatter_size,
                    min_x,
                    max_x,
                    min_y,
                    max_y,
                    xdist,
                    ydist,
                    hist_wh,
                    backend,
                )

            if all_plots is None:
                all_plots = mse_plot
            else:
                all_plots += mse_plot
        mse_plot = _set_final_opts(all_plots, layout_kwargs, kind, backend)

    save_plot_renderer(
        plot_obj=mse_plot,
        store_path=store_path,
        renderer=hv.Store.renderers[backend].instance(fig="png", dpi=300),
    )
    if backend != "bokeh" and plot_singles:
        _save_multiple_plots(
            mse_plot=mse_plot,
            alpha=alpha,
            scatter_kwargs=scatter_kwargs,
            layout_kwargs=layout_kwargs,
            aspect=aspect,
            fontscale=fontscale,
            xlabel=xlabel,
            out_params=out_params,
            store_path=store_path,
            backend=backend,
        )
