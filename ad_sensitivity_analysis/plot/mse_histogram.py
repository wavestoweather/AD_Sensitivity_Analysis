"""Plot predicted errors as histograms.

"""
from bokeh.palettes import Category20c
import holoviews as hv
from holoviews import opts
import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

from ad_sensitivity_analysis.plot.aux_functions import save_plot_renderer
from ad_sensitivity_analysis.plot.latexify import in_params_grouping


def _parse_df(df, plot_types, in_params, add_zero_sens, sens_key, error_key):
    """

    Parameters
    ----------
    df
    plot_types
    in_params
    add_zero_sens
    sens_key
    error_key

    Returns
    -------

    """
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
    else:
        mse_df = mse_df.loc[df[sens_key] != 0]
        inf_val = None

    # log_x
    mse_df[sens_key] = np.log10(mse_df[sens_key])
    if inf_val is not None:
        inf_val = np.log10(inf_val)
    # log_y
    mse_df[error_key] = np.abs(mse_df[error_key])
    mse_df[error_key] = np.log10(mse_df[error_key])
    return mse_df


# pylint: disable=too-many-arguments
def _plot_bokeh(
    mse_df,
    error_key,
    sens_key,
    by_col,
    alpha,
    cmap_values,
    width,
    height,
    xlabel,
    xlabel2,
    fontscale,
):
    """

    Parameters
    ----------
    mse_df
    error_key
    sens_key
    by_col
    alpha
    cmap_values
    width
    height
    xlabel
    xlabel2
    fontscale

    Returns
    -------

    """
    image = mse_df.hvplot.hist(
        y=error_key,
        by=by_col,
        alpha=alpha,
        legend=False,
        color=cmap_values,
        width=width,
        height=int(height / 2),
        xlabel=xlabel,
    ).opts(fontscale=fontscale) + mse_df.hvplot.hist(
        y=sens_key,
        by=by_col,
        alpha=alpha,
        legend="top_left",
        color=cmap_values,
        width=width,
        height=int(height / 2),
        xlabel=xlabel2,
    ).opts(
        fontscale=fontscale
    )
    return image.cols(1)


def _plot_mpl(
    mse_df,
    error_key,
    sens_key,
    by_col,
    alpha,
    cmap_values,
    xlabel,
    xlabel2,
    fontscale,
    groups,
):
    """

    Parameters
    ----------
    mse_df
    error_key
    sens_key
    by_col
    alpha
    cmap_values
    xlabel
    xlabel2
    fontscale
    groups

    Returns
    -------

    """
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
    # pylint: disable=no-member
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
    return plot1 + plot2 * legend_overlay


# pylint: disable=too-many-locals, too-many-arguments
def plot_histogram(
    df,
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

    alpha = 0.5
    f_limits = (-2, 2)

    error_key = "Mean Squared Error"
    sens_key = "Predicted Squared Error"

    fontscale = width / 350
    rcParams["axes.formatter.limits"] = f_limits

    mse_df = _parse_df(df, plot_types, in_params, add_zero_sens, sens_key, error_key)
    hv.extension(backend)

    if xlabel is None or xlabel == "":
        xlabel = "Log MSD"
    if xlabel2 is None or xlabel2 == "":
        xlabel2 = r"Log Predicted MSD"

    if title is None:
        title = "Deviation by Perturbed Parameter"
    if backend == "matplotlib":
        colors = plt.get_cmap("tab20c")
        cmap_types = {
            "artificial": mpl_colors.to_hex(colors(1)[0:-1]),
            "artificial (threshold)": mpl_colors.to_hex(colors(5)[0:-1]),
            "physical": mpl_colors.to_hex(colors(9)[0:-1]),
            "physical (high variability)": mpl_colors.to_hex(colors(13)[0:-1]),
            "1-moment": mpl_colors.to_hex(colors(17)[0:-1]),
        }
    else:
        colors = Category20c[20]
        cmap_types = {
            "artificial": colors[1],
            "artificial (threshold)": colors[5],
            "physical": colors[9],
            "physical (high variability)": colors[13],
            "1-moment": colors[17],
        }
    cmap_values = []
    groups = np.unique(mse_df["Group"])
    for group in groups:
        cmap_values.append(mpl_colors.to_hex(cmap_types[group]))
    by_col = "Group"
    if backend == "bokeh":
        image = _plot_bokeh(
            mse_df,
            error_key,
            sens_key,
            by_col,
            alpha,
            cmap_values,
            width,
            height,
            xlabel,
            xlabel2,
            fontscale,
        )
    else:
        image = _plot_mpl(
            mse_df,
            error_key,
            sens_key,
            by_col,
            alpha,
            cmap_values,
            xlabel,
            xlabel2,
            fontscale,
            groups,
        )

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
    save_plot_renderer(
        image, store_path, hv.Store.renderers[backend].instance(fig="png", dpi=300)
    )
