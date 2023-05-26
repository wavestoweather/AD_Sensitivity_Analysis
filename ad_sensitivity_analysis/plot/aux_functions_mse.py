"""helper functions to create plots with predicted errors.

"""
import holoviews as hv
from holoviews import opts
from holoviews.operation.datashader import datashade as dsshade
from matplotlib import colors as mpl_colors
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import seaborn as sns

from ad_sensitivity_analysis.plot.colors import get_cmap_types, cmap_particles
from ad_sensitivity_analysis.plot.latexify import in_params_grouping
from ad_sensitivity_analysis.plot.aux_functions import save_plot_renderer


def _set_abs_log(df, log, key, x_axis=False, log_func=np.log10):
    df[key] = np.abs(df[key])
    if "Real Predicted Squared Error" in df and x_axis:
        df["Real Predicted Squared Error"] = np.abs(df["Real Predicted Squared Error"])
    if log:
        df[key] = log_func(df[key])

        if "Real Predicted Squared Error" in df and x_axis:
            df["Real Predicted Squared Error"] = log_func(
                df["Real Predicted Squared Error"]
            )


def _set_labels(
    sens_key, error_key, xlabel=None, ylabel=None, log_x=False, log_y=False
):
    if xlabel is None:
        if log_x:
            if sens_key == "Predicted Squared Error":
                xlabel = r"$$\text{AD-Estimated} \log_{10} \text{MSD}$$"
            else:
                xlabel = "Log " + sens_key
        else:
            xlabel = sens_key
    if ylabel is None:
        if log_y:
            if error_key == "Mean Squared Error":
                ylabel = r"$$\text{Ensemble-Estimated} \log_{10} \text{MSD}$$"
            else:
                ylabel = "Log " + error_key
        else:
            ylabel = error_key
    return xlabel, ylabel


# pylint: disable=too-many-arguments, too-many-locals
def _plot_all_into_one_mpl(
    mse_df,
    out_params,
    sens_key,
    error_key,
    title,
    datashade,
    xticks,
    yticks,
    xlabel,
    ylabel,
    fontscale,
    aspect,
    backend,
    alpha=0.3,
    plot_types=False,
):
    """

    Parameters
    ----------
    mse_df
    out_params
    sens_key
    error_key
    title
    datashade
    xticks
    yticks
    xlabel
    ylabel
    fontscale
    aspect
    backend
    alpha
    plot_types

    Returns
    -------

    """
    if not plot_types:
        cmap_values = []
        for out_param in out_params:
            cmap_values.append(mpl_colors.to_hex(cmap_particles[out_param]))
            by_col = "Output Parameter"
    else:
        cmap_values = []
        cmap_types = get_cmap_types(backend=backend)
        for group in np.unique(mse_df["Group"]):
            cmap_values.append(mpl_colors.to_hex(cmap_types[group]))
        by_col = "Group"

    mse_plot = (
        mse_df.hvplot.scatter(
            x=sens_key,
            y=error_key,
            by=by_col,
            title=title,
            color=cmap_values,
            datashade=datashade,
            alpha=alpha,
            legend=True,
            yticks=yticks,
            xticks=xticks,
        )
        .opts(aspect=aspect, fontscale=fontscale)
        .options(xlabel=xlabel, ylabel=ylabel)
    )
    return mse_plot, by_col


# pylint: disable=too-many-arguments, too-many-locals
def _plot_all_into_one_bokeh(
    mse_df,
    sens_key,
    error_key,
    title,
    xticks,
    xlabel,
    ylabel,
    fontscale,
    backend,
    width,
    height,
    scatter_kwargs,
    layout_kwargs,
    inf_val=None,
    log_x=False,
    alpha=0.3,
    plot_types=False,
):
    """

    Parameters
    ----------
    mse_df
    sens_key
    error_key
    title
    xticks
    xlabel
    ylabel
    fontscale
    backend
    width
    height
    scatter_kwargs
    layout_kwargs
    inf_val
    log_x
    alpha
    plot_types

    Returns
    -------

    """
    if log_x:
        min_x = np.min(mse_df[sens_key])
        max_x = np.max(mse_df[sens_key])

        delta_x = (max_x - min_x) / 20
        max_x += delta_x
        min_x -= delta_x
    else:
        delta_x = np.max(mse_df[sens_key]) - np.min(mse_df[sens_key])
        min_x = np.min(mse_df[sens_key]) - delta_x / 10
        max_x = np.max(mse_df[sens_key]) + delta_x / 10
    if not plot_types:
        cmap = plt.get_cmap("tab10")
        colors = {}
        cmap_values = []
        for i, ratio_type in enumerate(np.unique(mse_df["Ratio Type"])):
            colors[ratio_type] = mpl_colors.to_hex(cmap(i)[0:-1])
            cmap_values.append(colors[ratio_type])
        by_col = "Ratio Type"
    else:
        cmap_values = []
        colors = {}
        cmap_types = get_cmap_types(backend=backend)
        for group in np.unique(mse_df["Group"]):
            colors[group] = cmap_types[group]
            cmap_values.append(cmap_types[group])
        by_col = "Group"

    legend = "top_left"
    opts_dic = {"fontscale": fontscale}

    min_y = mse_df[error_key].min()
    max_y = mse_df[error_key].max()
    delta_y = (max_y - min_y) / 20
    min_y -= delta_y
    max_y += delta_y
    if inf_val is not None:
        xticks_list = [(inf_val, "-Infinity")]
        tick_val = int(np.ceil(inf_val / 10.0)) * 10
        delta_tick = int(-tick_val / (xticks - 1))
        tick_val += delta_tick
        stop_crit = max_x - 2
        if max_x < 0:
            stop_crit = 0
        while tick_val < stop_crit:
            xticks_list.append((tick_val, tick_val))
            tick_val += delta_tick

        mse_plot = (
            mse_df.hvplot.scatter(
                x=sens_key,
                y=error_key,
                by=by_col,
                datashade=False,
                alpha=alpha,
                legend=legend,
                title=title,
                grid=True,
                xlim=(min_x, max_x),
                ylim=(min_y, max_y),
                color=cmap_values,
                width=width,
                height=height,
            )
            .opts(opts.Scatter(**scatter_kwargs))  # pylint: disable=no-member
            .opts(**opts_dic)
            .options(
                ylabel=ylabel,
                xlabel=xlabel,
                xticks=xticks_list,
                width=width,
                height=height,
            )
        )
        # Save space by removing legend title
        mse_plot.get_dimension(by_col).label = ""
    else:
        mse_plot = (
            mse_df.hvplot.scatter(
                x=sens_key,
                y=error_key,
                by=by_col,
                datashade=False,
                alpha=alpha,
                legend=legend,
                title=title,
                grid=True,
                xlim=(min_x, max_x),
                ylim=(min_y, max_y),
                color=cmap_values,
                width=width,
                height=height,
            )
            .opts(opts.Scatter(**scatter_kwargs))  # pylint: disable=no-member
            .opts(**opts_dic)
            .options(
                ylabel=ylabel,
                xlabel=xlabel,
                width=width,
                height=height,
            )
        )
    return mse_plot.opts(opts.Layout(**layout_kwargs))  # pylint: disable=no-member


# pylint: disable=too-many-arguments, too-many-locals
def _prepare_errors_plot(
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
):
    """

    Parameters
    ----------
    df
    out_param
    in_params
    latex
    plot_types
    width
    height
    add_zero_sens
    log_x
    log_y
    x_key
    y_key

    Returns
    -------

    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})

    df_tmp = df.copy()
    df_tmp = df_tmp.loc[df_tmp["Output Parameter"] == out_param]
    df_tmp = df_tmp.loc[df_tmp["Input Parameter"].isin(in_params)]
    if plot_types:
        # We need to add a column 'Group' to plot it correctly
        tmp_dic = {}
        for in_p in in_params:
            for group, params in in_params_grouping.items():
                if in_p in params:
                    tmp_dic[in_p] = group
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

    if log_x:
        df_tmp[x_key] = np.log10(np.abs(df_tmp[x_key]))
    if log_y:
        df_tmp[y_key] = np.log10(np.abs(df_tmp[y_key]))
    return df_tmp, group_colors


def _grid_paper_plot(
    tmp_df,
    sens_key,
    error_key,
    by_col,
    width,
    height,
    xlabel,
    ylabel,
    scatter_kwargs,
    cmap_values,
    min_x,
    max_x,
    min_y,
    max_y,
    opts_dic2,
    legend,
    title,
    alpha=0.3,
):
    """

    Parameters
    ----------
    tmp_df
    sens_key
    error_key
    by_col
    width
    height
    xlabel
    ylabel
    scatter_kwargs
    cmap_values
    min_x
    max_x
    min_y
    max_y
    opts_dic2
    legend
    title
    alpha

    Returns
    -------

    """
    mse_plot = (
        tmp_df.hvplot.scatter(
            x=sens_key,
            y=error_key,
            by=by_col,
            datashade=False,
            alpha=alpha,
            legend=legend,
            title=title,
            grid=True,
            xlim=(min_x, max_x),
            ylim=(min_y, max_y),
            color=cmap_values,
            width=width,
            height=height,
        )
        .opts(opts.Scatter(**scatter_kwargs))  # pylint: disable=no-member
        .opts(**opts_dic2)
        .options(
            ylabel=ylabel,
            xlabel=xlabel,  # ""
            width=width,
            height=height,
        )
    )
    return mse_plot


def _grid_paper_inf_plot(
    tmp_df,
    sens_key,
    error_key,
    by_col,
    width,
    height,
    xlabel,
    ylabel,
    xticks,
    scatter_kwargs,
    cmap_values,
    min_x,
    max_x,
    min_y,
    max_y,
    opts_dic2,
    legend,
    title,
    inf_val=None,
    alpha=0.3,
):
    """

    Parameters
    ----------
    tmp_df
    sens_key
    error_key
    by_col
    width
    height
    xlabel
    ylabel
    xticks
    scatter_kwargs
    cmap_values
    min_x
    max_x
    min_y
    max_y
    opts_dic2
    legend
    title
    inf_val
    alpha

    Returns
    -------

    """
    xticks_list = [(inf_val, "-Infinity")]
    tick_val = int(np.ceil(inf_val / 10.0)) * 10
    delta_tick = int(-tick_val / (xticks - 1))
    tick_val += delta_tick
    stop_crit = max_x - 2
    if max_x < 0:
        stop_crit = 0
    while tick_val < stop_crit:
        xticks_list.append((tick_val, tick_val))
        tick_val += delta_tick
    mse_plot = (
        tmp_df.hvplot.scatter(
            x=sens_key,
            y=error_key,
            by=by_col,
            datashade=False,
            alpha=alpha,
            legend=legend,
            title=title,
            grid=True,
            xlim=(min_x, max_x),
            ylim=(min_y, max_y),
            color=cmap_values,
            width=width,
            height=height,
        )
        .opts(opts.Scatter(**scatter_kwargs))  # pylint: disable=no-member
        .opts(**opts_dic2)
        .options(
            ylabel=ylabel,
            xlabel=xlabel,
            xticks=xticks_list,
            width=width,
            height=height,
        )
    )
    # Save space by removing legend title
    mse_plot.get_dimension(by_col).label = ""
    return mse_plot


def _get_kwargs(
    width,
    height,
    backend,
    scatter_size,
    formatter_limits=None,
):
    """

    Parameters
    ----------
    width
    height
    backend
    scatter_size
    log_func
    formatter_limits

    Returns
    -------

    """
    dsshade.dynamic = False
    fontscale = width / 2200
    if formatter_limits is not None:
        rcParams["axes.formatter.limits"] = formatter_limits
    hv.extension(backend)
    aspect = width / height
    if scatter_size is None:
        scatter_size = int(height / 200)
    hist_wh = int(height / 3)

    if backend == "matplotlib":
        scatter_kwargs = {"s": scatter_size}
    else:
        scatter_kwargs = {"size": scatter_size}
    fig_inches = width / 300
    if width < height:
        fig_inches = height / 300
    layout_kwargs = {}
    if backend == "bokeh":
        layout_kwargs["width"] = width
        layout_kwargs["height"] = height
    else:
        layout_kwargs["fig_inches"] = fig_inches
    return layout_kwargs, scatter_kwargs, hist_wh, aspect, fontscale


def _get_title_ops(title, backend, xlabel, ylabel, fontscale, aspect):
    """

    Parameters
    ----------
    title
    backend
    xlabel
    ylabel
    fontscale
    aspect

    Returns
    -------

    """
    if title is None:
        title = ylabel + " over " + xlabel
    opts_dic = {"fontscale": fontscale}
    opts_dic2 = {"fontscale": fontscale}
    if backend == "matplotlib":
        opts_dic["aspect"] = aspect / 2
        opts_dic2["aspect"] = aspect
    return title, opts_dic, opts_dic2


def _save_multiple_plots(
    mse_plot,
    alpha,
    scatter_kwargs,
    layout_kwargs,
    aspect,
    fontscale,
    xlabel,
    out_params,
    store_path,
    backend="matplotlib",
):
    """

    Parameters
    ----------
    mse_plot
    alpha
    scatter_kwargs
    layout_kwargs
    aspect
    fontscale
    xlabel
    out_params
    store_path
    backend

    Returns
    -------

    """
    for j, plot in enumerate(mse_plot):
        if j > 0:
            if isinstance(plot, hv.core.overlay.Overlay):
                plot.opts(
                    opts.Scatter(  # pylint: disable=no-member
                        alpha=alpha,
                        show_grid=True,
                        show_legend=True,
                        **scatter_kwargs,
                    ),
                    opts.Layout(**layout_kwargs),  # pylint: disable=no-member
                ).opts(
                    aspect=aspect,
                    fontscale=fontscale,
                    xlabel=xlabel,
                    xaxis="bottom",
                    title=title,
                )
            else:
                title = "Validation Sensitivity " + out_params[j - 1]
                # pylint: disable=unnecessary-dunder-call
                plot.__setitem__(
                    "main",
                    plot.main()
                    .opts(
                        opts.Scatter(  # pylint: disable=no-member
                            alpha=alpha,
                            show_grid=True,
                            show_legend=True,
                            **scatter_kwargs,
                        ),
                        opts.Layout(**layout_kwargs),  # pylint: disable=no-member
                    )
                    .opts(
                        aspect=aspect,
                        fontscale=fontscale,
                        xlabel=xlabel,
                        xaxis="bottom",
                        title=title,
                    ),
                )
        else:
            plot.opts(show_legend=False)
        save_plot_renderer(
            plot_obj=mse_plot,
            store_path=store_path,
            renderer=hv.Store.renderers[backend].instance(fig="png", dpi=300),
        )


def _get_x_limits(mse_df, sens_key, log_x, kind):
    """

    Parameters
    ----------
    mse_df
    sens_key
    log_x
    kind

    Returns
    -------

    """
    if log_x:
        if kind == "paper":
            min_x = np.min(mse_df[sens_key])
            max_x = np.max(mse_df[sens_key])

            delta_x = (max_x - min_x) / 20
            max_x += delta_x
            min_x -= delta_x
        else:
            min_x = np.min(mse_df[sens_key]) - np.abs(np.min(mse_df[sens_key])) / 7
            max_x = np.max(mse_df[sens_key]) + np.abs(np.max(mse_df[sens_key])) / 7
    else:
        delta_x = np.max(mse_df[sens_key]) - np.min(mse_df[sens_key])
        min_x = np.min(mse_df[sens_key]) - delta_x / 10
        max_x = np.max(mse_df[sens_key]) + delta_x / 10
    return min_x, max_x


def _get_color_vals(mse_df, backend, plot_types):
    """

    Parameters
    ----------
    mse_df
    backend
    plot_types

    Returns
    -------

    """
    if not plot_types:
        cmap = plt.get_cmap("tab10")
        colors = {}
        cmap_values = []
        if "Ratio Type" in mse_df:
            for i, ratio_type in enumerate(np.unique(mse_df["Ratio Type"])):
                colors[ratio_type] = mpl_colors.to_hex(cmap(i)[0:-1])
                cmap_values.append(colors[ratio_type])
            by_col = "Ratio Type"
        else:
            by_col = None
            cmap_values.append(mpl_colors.to_hex(cmap(0)[0:-1]))
    else:
        cmap_values = []
        colors = {}
        cmap_types = get_cmap_types(backend=backend)
        for group in np.unique(mse_df["Group"]):
            colors[group] = cmap_types[group]
            cmap_values.append(cmap_types[group])
        by_col = "Group"
    return cmap_values, colors, by_col


def _add_histogram(
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
    backend="matplotlib",
):
    """

    Parameters
    ----------
    mse_plot
    width
    height
    scatter_size
    min_x
    max_x
    min_y
    max_y
    xdist
    ydist
    hist_wh
    backend

    Returns
    -------

    """
    if backend == "bokeh":
        mse_plot = (
            mse_plot.opts(width=width, height=height).options(
                # Legend scatter size is very small by default
                **{
                    "glyph_height": scatter_size * 5,
                    "glyph_width": scatter_size * 5,
                }
            )
            << ydist.opts(
                yaxis="bare",
                xaxis="bare",
                xlim=(min_y, max_y),
                width=hist_wh,
                height=height,
            )
            << xdist.opts(
                xaxis="bare",
                yaxis="bare",
                xlim=(min_x, max_x),
                height=hist_wh,
                width=width,
            )
        )
    else:
        mse_plot = (
            mse_plot.options(
                # Legend scatter size is very small by default
                **{
                    "glyph_height": scatter_size * 5,
                    "glyph_width": scatter_size * 5,
                }
            )
            << ydist.opts(yaxis="bare", xaxis="bare", xlim=(min_y, max_y))
            << xdist.opts(xaxis="bare", yaxis="bare", xlim=(min_x, max_x))
        )
    return mse_plot


def _prepare_mse_plot(
    df,
    out_params,
    kind,
    plot_kind,
    plot_types,
    legend_pos,
    add_zero_sens,
    title,
    log_func=np.log10,
):
    """

    Parameters
    ----------
    df
    out_params
    kind
    plot_kind
    plot_types
    legend_pos
    add_zero_sens
    title
    log_func

    Returns
    -------

    """
    in_params = np.unique(df["Input Parameter"])
    legend = False
    if kind == "paper":
        legend = legend_pos
    if isinstance(out_params, str):
        out_params = [out_params]

    if plot_kind in ("paper", "single_plot"):
        alpha = 0.5
    else:
        alpha = 1

    error_key = "Mean Squared Error"
    sens_key = "Predicted Squared Error"

    mse_df = df.copy()
    if plot_types:
        # We need to add a column 'Group' to plot it correctly
        tmp_dic = {}
        for in_p in in_params:
            for group, params in in_params_grouping.items():
                if in_p in params:
                    tmp_dic[in_p] = group
                    break
        mse_df["Group"] = mse_df.apply(
            lambda row: tmp_dic[row["Input Parameter"]], axis=1
        )

    if add_zero_sens:
        inf_val = np.min(mse_df.loc[mse_df[sens_key] != 0][sens_key]) / 10
        mse_df[sens_key] = mse_df[sens_key].replace({0: inf_val})
        inf_val = log_func(inf_val)
    else:
        mse_df = mse_df.loc[mse_df[sens_key] != 0]
        inf_val = None

    if "NC" in title or "_OUT" in title:
        title = title.replace("Estimation", "Est.")
    return mse_df, inf_val, title, error_key, sens_key, legend, out_params, alpha


def _set_final_opts(all_plots, layout_kwargs, kind, backend):
    """

    Parameters
    ----------
    all_plots
    layout_kwargs
    kind
    backend

    Returns
    -------

    """
    if kind == "paper":
        if backend == "matplotlib":
            mse_plot = all_plots.opts(sublabel_format="", tight=True).opts(
                opts.Layout(**layout_kwargs)  # pylint: disable=no-member
            )
        else:
            mse_plot = all_plots.opts(
                opts.Layout(**layout_kwargs)  # pylint: disable=no-member
            )
        return mse_plot
    if backend == "matplotlib":
        mse_plot = (
            all_plots.cols(1)
            .opts(sublabel_format="", tight=True)
            .opts(opts.Layout(**layout_kwargs))  # pylint: disable=no-member
        )
    else:
        mse_plot = all_plots.cols(1).opts(
            opts.Layout(**layout_kwargs)  # pylint: disable=no-member
        )
    return mse_plot
