"""Plot predicted error over time.

"""
from bokeh.models import Range1d, LinearAxis, GlyphRenderer
import holoviews as hv
from holoviews import opts
import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import seaborn as sns

from ad_sensitivity_analysis.plot.latexify import parse_word
from ad_sensitivity_analysis.plot.aux_functions import save_plot_renderer


def _get_sizes(backend, font_scale, height, width, dot_size):
    """

    Parameters
    ----------
    backend
    font_scale
    height
    width
    dot_size

    Returns
    -------

    """
    title_font_size = 12
    if backend == "matplotlib":
        if font_scale is None:
            if isinstance(height, float) or height < 100:
                font_scale = height / 3 * 0.8
            else:
                font_scale = height / 900 * 0.8
        if dot_size is None:
            dot_size = 2

        axis_label_text_font_size = 10
        major_label_text_font_size = 10
        aspect = width / height
    else:
        if font_scale is None:
            if isinstance(height, float) or height < 100:
                font_scale = height * 1.2
            else:
                font_scale = height / 350
        axis_label_text_font_size = int(13 * font_scale * 0.8)
        major_label_text_font_size = int(11 * font_scale * 0.8)
        if dot_size is None:
            dot_size = 14
        aspect = None
    hv.extension(backend)
    fontsize = {
        "legend": major_label_text_font_size * font_scale,
        "xticks": major_label_text_font_size * font_scale,
        "yticks": major_label_text_font_size * font_scale,
        "title": title_font_size * font_scale,
        "xlabel": axis_label_text_font_size * font_scale,
        "ylabel": axis_label_text_font_size * font_scale,
    }
    return fontsize, aspect, dot_size


def _get_df_within_limits(df, min_x, max_x, x_limits, y, logy, logtwin):
    """

    Parameters
    ----------
    df
    min_x
    max_x
    x_limits
    y
    logy
    logtwin

    Returns
    -------

    """

    limits_dict = {}
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
    limits_dict["min_log_y"] = np.NaN
    if logy:
        df[y] = np.log10(df[y])
        limits_dict["min_log_y"] = np.nanmin(df.loc[df[y] != -np.inf][y]) - 1
        df.replace(-np.inf, limits_dict["min_log_y"], inplace=True)
    limits_dict["min_log_twin"] = np.NaN
    if logtwin:
        df["Predicted Error"] = np.log10(np.abs(df["Predicted Error"]))
        limits_dict["min_log_twin"] = (
            np.nanmin(df.loc[df["Predicted Error"] != -np.inf]["Predicted Error"]) - 1
        )
        df.replace(-np.inf, limits_dict["min_log_twin"], inplace=True)
    limits_dict["lower_y"] = np.min(df["Predicted Error"])
    limits_dict["upper_y"] = np.max(df["Predicted Error"])
    delta = (limits_dict["upper_y"] - limits_dict["lower_y"]) / 12
    limits_dict["upper_y"] += delta
    limits_dict["lower_y"] -= delta

    limits_dict["lower_y2"] = np.min(df[y])
    limits_dict["upper_y2"] = np.max(df[y])
    delta2 = (limits_dict["upper_y2"] - limits_dict["lower_y2"]) / 12
    limits_dict["lower_y2"] -= delta2
    limits_dict["upper_y2"] += delta2

    limits_dict["min_t"] = np.min(df["time_after_ascent"])
    limits_dict["max_t"] = np.max(df["time_after_ascent"])

    return df, limits_dict


def _get_colors(df, perturb_delta, limits_dict):
    """

    Parameters
    ----------
    df
    perturb_delta
    limits_dict

    Returns
    -------

    """
    if perturb_delta is not None:
        perturb_lines = None
        t = limits_dict["min_t"]
        while t < limits_dict["max_t"]:
            if perturb_lines is None:
                perturb_lines = hv.VLine(x=t).opts(color="black")
            else:
                perturb_lines *= hv.VLine(x=t).opts(color="black")
            t += perturb_delta

    if len(np.unique(df["Input Parameter"])) < 10:
        cmap = plt.get_cmap("tab10")
    else:
        cmap = plt.get_cmap("tab20")

    cmap_map = {}
    cmap_values = [mpl_colors.to_hex(cmap(0)[0:-1])]
    data_types = [parse_word(np.unique(df["Output Parameter"])[0])]
    data_types_2 = [np.unique(df["Output Parameter"])[0]]
    for i, in_param in enumerate(np.unique(df["Input Parameter"])):
        data_types.append(parse_word(in_param).replace(r"\text", r"\mathrm"))
        cmap_values.append(mpl_colors.to_hex(cmap(i + 1)[0:-1]))
        cmap_map[in_param] = mpl_colors.to_hex(cmap(i + 1)[0:-1])
        data_types_2.append(in_param)
    return data_types, data_types_2, cmap_values, cmap_map, perturb_lines


# pylint: disable=too-many-arguments, too-many-locals
def _plot_mpl(
    df,
    y,
    limits_dict,
    logtwin,
    logy,
    data_types,
    data_types_2,
    aspect,
    cmap_values,
    cmap_map,
    alpha,
    dot_size,
    perturb_delta,
    perturb_lines,
    fontsize,
    twinlabel,
    font_scale,
    xlabel,
    ylabel,
    title,
    precision,
):
    """

    Parameters
    ----------
    df
    y
    limits_dict
    logtwin
    logy
    data_types
    data_types_2
    aspect
    cmap_values
    cmap_map
    alpha
    dot_size
    perturb_delta
    perturb_lines
    fontsize
    twinlabel
    font_scale
    xlabel
    ylabel
    title
    precision

    Returns
    -------

    """

    def apply_axis_format(
        plot,
        element,  # pylint: disable=unused-argument
    ):
        ax = plot.handles["axis"]
        ax.set_ylim((limits_dict["lower_y2"], limits_dict["upper_y2"]))
        ax.set_title(
            title,
            loc="left",
            **{
                "fontweight": "bold",
                "fontsize": fontsize["title"],
            },
        )
        ax.set_ylabel(
            ylabel,
            **{
                "rotation": 90,
                "fontsize": fontsize["ylabel"],
            },
        )
        ax.set_xlabel(
            xlabel,
            **{
                "fontsize": fontsize["xlabel"],
            },
        )
        ax.xaxis.set_ticks_position("none")

        # pylint: disable=unused-argument
        def format_fn(tick_val, tick_pose):
            # Can implement an even fancier format here if needed
            if logy and tick_val == limits_dict["min_log_y"]:
                return "-Infinity"
            scientific = np.format_float_scientific(
                tick_val, precision=precision, exp_digits=1, unique=True
            )
            if precision == 0:
                return scientific.replace(".", "")
            return scientific

        ax.yaxis.set_major_formatter(tick.FuncFormatter(format_fn))
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=fontsize["xticks"],
        )
        ax.xaxis.set_ticks_position("none")
        ax.set_yticks(np.linspace(ax.get_ybound()[0], ax.get_ybound()[1], 5))
        plot.handles["axis"] = ax

    # pylint: disable=no-member
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
                s=dot_size,
            )
        )
    )

    def twinx_per_timestep(
        plot,
        element,  # pylint: disable=unused-argument
    ):
        ax = plot.handles["axis"]
        ax.xaxis.set_ticks_position("none")
        twinax = ax.twinx()

        twinax.set_ylim((limits_dict["lower_y"], limits_dict["upper_y"]))
        twinax.set_ylabel(
            twinlabel,
            labelpad=13.0 * font_scale,
            **{
                "rotation": 270,
                "fontsize": fontsize["ylabel"],
            },
        )
        if logtwin:
            yticks_tmp = twinax.get_yticks()
            yticks_tmp[0] = limits_dict["min_log_twin"]
            if yticks_tmp[1] - yticks_tmp[0] < 3.5:
                yticks_tmp = np.delete(yticks_tmp, 1)
            twinax.set_yticks(yticks_tmp)

        def format_fn(
            tick_val,
            tick_pose,  # pylint: disable=unused-argument
        ):
            # Can implement an even fancier format here if needed
            if logtwin and tick_val == limits_dict["min_log_twin"]:
                return "-Infinity"
            scientific = np.format_float_scientific(
                tick_val, precision=precision, exp_digits=1, unique=True
            )
            if precision == 0:
                return scientific.replace(".", "")
            return scientific

        twinax.yaxis.set_major_formatter(tick.FuncFormatter(format_fn))
        twinax.tick_params(
            axis="both",
            which="major",
            labelsize=fontsize["xticks"],
        )
        twinax.xaxis.set_ticks_position("none")
        twinax.set_yticks(
            np.linspace(twinax.get_ybound()[0], twinax.get_ybound()[1], 5)
        )
        plot.handles["axis"] = twinax

    def twinx2(
        plot,
        element,  # pylint: disable=unused-argument
    ):
        ax = plot.handles["axis"]
        twinax = ax.twinx()
        twinax.set_ylim((limits_dict["lower_y"], limits_dict["upper_y"]))
        twinax.set_yticks([])
        plot.handles["axis"] = twinax
        twinax.xaxis.set_ticks_position("none")

    twin = None
    for i, in_p in enumerate(data_types_2[1::]):
        if twin is None:
            # pylint: disable=no-member
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
                        s=dot_size,
                        cmap=cmap_map,
                    ),
                )
            )
        else:
            # pylint: disable=no-member
            twin = twin * df.loc[df["Input Parameter"] == in_p].hvplot.scatter(
                x="time_after_ascent",
                y="Predicted Error",
                alpha=alpha,
                legend=False,
                aspect=aspect,
                color=cmap_values[i + 1],
            ).opts(initial_hooks=[twinx2], apply_ranges=False,).opts(
                opts.Scatter(
                    s=dot_size,
                    cmap=cmap_map,
                ),
            )
    # pylint: disable=no-member
    legend_overlay = hv.NdOverlay(
        {
            data_types[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                opts.Scatter(
                    s=dot_size * 4,
                    color=cmap_values[i],
                )
            )
            for i in range(len(data_types))
        }
    ).opts(fontsize=fontsize)
    return (
        legend_overlay * left * twin
        if perturb_delta is None
        else left * twin * perturb_lines * legend_overlay
    )


# pylint: disable=too-many-arguments, too-many-locals
def _plot_bokeh(
    df,
    y,
    limits_dict,
    data_types,
    cmap_values,
    alpha,
    dot_size,
    perturb_delta,
    perturb_lines,
    fontsize,
    twinlabel,
    xlabel,
    ylabel,
):
    """

    Parameters
    ----------
    df
    y
    limits_dict
    data_types
    cmap_values
    alpha
    dot_size
    perturb_delta
    perturb_lines
    fontsize
    twinlabel
    xlabel
    ylabel

    Returns
    -------

    """
    left = df.hvplot.scatter(
        x="time_after_ascent",
        y=y,
        alpha=alpha,
        label=data_types[0],
        ylabel=ylabel,
        xlabel=xlabel,
        legend="top_right",
        size=dot_size,
        color=[cmap_values[0]],
    )

    # pylint: disable=unused-argument
    def apply_twin_axis(plot, element):
        p = plot.state
        p.extra_y_ranges = {
            "twiny": Range1d(start=limits_dict["lower_y"], end=limits_dict["upper_y"])
        }
        p.add_layout(
            LinearAxis(
                y_range_name="twiny",
                axis_label_text_font_size=f"{fontsize['ylabel']}pt",
                major_label_text_font_size=f"{fontsize['xticks']}pt",
                axis_label=twinlabel,
            ),
            "right",
        )
        glyph = p.select({"type": GlyphRenderer})[0]
        glyph.y_range_name = "twiny"

    df["Input Parameter"] = parse_word(df["Input Parameter"])
    twin = df.hvplot.scatter(
        x="time_after_ascent",
        y="Predicted Error",
        by="Input Parameter",
        alpha=alpha,
        size=dot_size,
        colormap=cmap_values,
    ).opts(
        hooks=[apply_twin_axis],
        apply_ranges=True,
    )

    if perturb_delta is not None:
        image = left * twin * perturb_lines
    else:
        image = left * twin
    return image


def _set_image_size(image, backend, width, height, title):
    """

    Parameters
    ----------
    image
    backend
    width
    height
    title

    Returns
    -------

    """
    if backend == "matplotlib":
        if isinstance(width, float) or width < 100:
            image = image.opts(
                fig_inches=(width, height),
            )
        else:
            image = image.opts(
                fig_inches=(width / 300.0, height / 300.0),
            )
    else:
        if isinstance(width, float) or width < 100:
            image = image.opts(
                width=int(width * 300),
                height=int(height * 300),
                title=title,
            )
        else:
            image = image.opts(
                width=width,
                height=height,
                title=title,
            )
    return image


# pylint: disable=too-many-arguments, too-many-locals
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
    dot_size=2,
    latex=False,
    trajectory=None,
    out_param=None,
    in_params=None,
    precision=0,
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
    dot_size : int
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
    precision : int
        Maximum number of digits for ticks used in
        numpy.format_float_scientific if backend "matplotlib" is used.

    Returns
    -------
    matplotlib.figure.Figure with the plot drawn onto it.
    """
    sns.set(rc={"text.usetex": latex, "axes.grid": True})

    if trajectory is not None:
        df = df.loc[df["trajectory"] == trajectory]
    if out_param is not None:
        df = df.loc[df["Output Parameter"] == out_param]
    if in_params is not None:
        df = df.loc[df["Input Parameter"].isin(in_params)]
    y = "Not Perturbed Value"
    if plot_deviation:
        y = "Mean Squared Error"

    fontsize, aspect, dot_size = _get_sizes(
        backend, font_scale, height, width, dot_size
    )
    df, limits_dict = _get_df_within_limits(
        df, min_x, max_x, x_limits, y, logy, logtwin
    )
    data_types, data_types_2, cmap_values, cmap_map, perturb_lines = _get_colors(
        df, perturb_delta, limits_dict
    )

    if backend == "matplotlib":

        image = _plot_mpl(
            df,
            y,
            limits_dict,
            logtwin,
            logy,
            data_types,
            data_types_2,
            aspect,
            cmap_values,
            cmap_map,
            alpha,
            dot_size,
            perturb_delta,
            perturb_lines,
            fontsize,
            twinlabel,
            font_scale,
            xlabel,
            ylabel,
            title,
            precision=precision,
        )
    else:
        image = _plot_bokeh(
            df,
            y,
            limits_dict,
            data_types,
            cmap_values,
            alpha,
            dot_size,
            perturb_delta,
            perturb_lines,
            fontsize,
            twinlabel,
            xlabel,
            ylabel,
        )

    image = _set_image_size(image, backend, width, height, title)

    if save:
        if "." not in store_path:
            filetype = "png"
            tmp_path = store_path
        else:
            filetype = store_path.split(".")[-1]
            tmp_path = store_path.split(".")[0]
        save_plot_renderer(
            image, tmp_path, hv.Store.renderers[backend].instance(fig=filetype, dpi=300)
        )
    return image
