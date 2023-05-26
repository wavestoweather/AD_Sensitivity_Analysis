"""Helper functions to plot histograms based on statistics.

"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import seaborn as sns

from ad_sensitivity_analysis.plot.aux_functions import save_plot
from ad_sensitivity_analysis.plot.latexify import parse_word


# pylint: disable=too-many-locals, too-many-arguments
def _plot_hist(
    out_params,
    in_p,
    verbose,
    edges_in_params,
    hist_in_params,
    log,
    ticks_offset,
    width,
    fig=None,
    axs=None,
    filename=None,
    font_scale=None,
    title=None,
    save=False,
    interactive=False,
):

    if len(out_params) != 3:
        print("The number of output params should be three.")
        print("Future versions will support varying numbers.")
    if axs is None:
        ax1 = plt.subplot(311)
        ax2 = plt.subplot(312)
        ax3 = plt.subplot(313)

    else:
        ax1 = axs[0]
        ax2 = axs[1]
        ax3 = axs[2]
    if verbose:
        print(f"Plotting histogram w.r.t. {in_p}")

    def create_fig(ax, out_p, in_p, title=None):
        if in_p not in edges_in_params[out_p]:
            return
        if title is None:
            in_p_latex = parse_word(in_p).replace(r"\partial", "")
            title = f"{parse_word(out_p)} w.r.t. {in_p_latex}"
        ax_t = sns.barplot(
            x=edges_in_params[out_p][in_p][:-1],
            y=hist_in_params[out_p][in_p],
            color="seagreen",
            ax=ax,
        )
        if log:
            ax_t.set_yscale("log")
        if ticks_offset is None:
            if 24 > width > 5:
                local_ticks_offset = 24 // (width - 5)
            elif width <= 5:
                local_ticks_offset = 6
        else:
            local_ticks_offset = ticks_offset
        if local_ticks_offset >= len(edges_in_params[out_p][in_p][:-1]) - 1:
            local_ticks_offset = len(edges_in_params[out_p][in_p][:-1]) - 2
        x_labels = [
            f"{tick:1.1e}"
            for tick in edges_in_params[out_p][in_p][:-1:local_ticks_offset]
        ]
        x_ticks = np.arange(
            0, len(edges_in_params[out_p][in_p][:-1]) - 1, local_ticks_offset
        )
        ax_t.set(xticks=x_ticks)
        _ = ax_t.set_xticklabels(
            x_labels, rotation=45, ha="right", rotation_mode="anchor"
        )
        if font_scale is None:
            _ = ax_t.set_title(title)
        else:
            ax_t.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
            _ = ax_t.set_title(title, fontsize=int(12 * font_scale))

    create_fig(ax1, out_params[0], in_p, title)
    create_fig(ax2, out_params[1], in_p, title)
    create_fig(ax3, out_params[2], in_p, title)

    if interactive:
        # The view in jupyter notebooks may be different from the stored one without
        # the given padding.
        fig.tight_layout(h_pad=1)

    if filename is not None and save:
        save_plot(fig, filename)
        if not interactive:
            plt.clf()


def _get_kde_df(
    hist_in_params, edges_in_params, out_param, in_params, log_weights, filter_zero
):
    """

    Parameters
    ----------
    hist_in_params
    edges_in_params
    out_param
    in_params
    log_weights
    filter_zero

    Returns
    -------

    """
    data_dic = {
        "Parameter": [],
        "weight": [],
        "Impact": [],
    }
    for in_p in in_params:
        if log_weights:
            zeros = np.argwhere(hist_in_params[out_param][in_p] == 0)
            weights = np.log10(
                hist_in_params[out_param][in_p],
                where=(hist_in_params[out_param][in_p] >= 0.5),
            )
            weights += 0.1
            weights[zeros] = 0
        else:
            weights = np.ones(np.shape(hist_in_params[out_param][in_p]))

        if filter_zero:
            # We can only approximate zero gradients. We have to
            # find the bin that contains zero and set the
            # corresponding weight to zero
            for i, edge in enumerate(edges_in_params[out_param][in_p][:-1]):
                if edge <= 0 <= edges_in_params[out_param][in_p][i + 1]:
                    weights[i] = 0
                    break
        impact = (
            edges_in_params[out_param][in_p][:-1]
            + (
                edges_in_params[out_param][in_p][1]
                - edges_in_params[out_param][in_p][0]
            )
            / 2
        )
        data_dic["weight"].extend(weights)
        data_dic["Impact"].extend(impact)
        data_dic["Parameter"].extend(np.repeat(parse_word(in_p), len(weights)))

    return pd.DataFrame(data_dic)


def _get_corr_matrix_array(in_params, in_params_tmp, out_p, hist_in_params, verbose):
    """

    Parameters
    ----------
    in_params
    in_params_tmp
    out_p
    hist_in_params
    verbose

    Returns
    -------

    """
    corr_matrix = np.zeros((len(in_params_tmp), len(in_params)))
    if verbose:
        print(f"Create similarity matrix for {out_p}")
    for i, in_params_i in enumerate(in_params_tmp):
        for j, in_params_j in enumerate(in_params_tmp):
            if i == j:
                corr_matrix[i, j] = 1
            if i >= j:
                continue
            corr_matrix[i, j] = np.corrcoef(
                hist_in_params[out_p][in_params_i],
                hist_in_params[out_p][in_params_j],
            )[0][1]
            corr_matrix[j, i] = corr_matrix[i, j]
    return corr_matrix


def _get_params_ordered(corr_matrix, in_params_tmp, verbose):
    """

    Parameters
    ----------
    corr_matrix
    in_params_tmp
    verbose

    Returns
    -------

    """
    best_row = 0
    best_value = -1
    for i in range(len(in_params_tmp)):
        for j in range(len(in_params_tmp)):
            if i >= j:
                continue
            if best_value == -1:
                best_value = corr_matrix[i, j]
                best_row = i
            elif best_value < corr_matrix[i, j]:
                best_value = corr_matrix[i, j]
                best_row = i
    in_params_sorted = [in_params_tmp[best_row]]
    previous_rows = [best_row]
    if verbose:
        print("Generate matrix for plotting")
    for i in range(len(in_params_tmp) - 1):
        next_row = 0
        next_best_value = -1
        for j in range(len(in_params_tmp)):
            if j in previous_rows:
                continue
            if next_best_value == -1:
                next_best_value = corr_matrix[i, j]
                next_row = j
            elif next_best_value < corr_matrix[i, j]:
                next_best_value = corr_matrix[i, j]
                next_row = j
        in_params_sorted.append(in_params_tmp[next_row])
        previous_rows.append(next_row)
    return in_params_sorted


def _set_heat_hist_ticks(hist2d, ax, ticks_offset, x_ticks, y_ticks, width):
    """

    Parameters
    ----------
    hist2d
    ax
    ticks_offset
    x_ticks
    y_ticks
    width

    Returns
    -------

    """
    if ticks_offset is None:
        if 24 > width > 5:
            local_ticks_offset = 24 // (width - 5)
        elif width <= 5:
            local_ticks_offset = 6
    else:
        local_ticks_offset = ticks_offset

    if local_ticks_offset >= np.shape(hist2d)[0]:
        local_ticks_offset_x = np.shape(hist2d)[0] - 1
    else:
        local_ticks_offset_x = local_ticks_offset
    x_ticks_location = np.arange(0, np.shape(hist2d)[0], local_ticks_offset_x)
    ax.set_xticks(x_ticks_location)
    if local_ticks_offset >= np.shape(hist2d)[1]:
        local_ticks_offset_y = np.shape(hist2d)[1] - 1
    else:
        local_ticks_offset_y = local_ticks_offset
    y_ticks_location = np.arange(0, np.shape(hist2d)[1], local_ticks_offset_y)
    ax.set_yticks(y_ticks_location)
    x_labels = [f"{x_ticks[i]:1.1e}" for i in x_ticks_location]
    _ = ax.set_xticklabels(x_labels, rotation=45, ha="right", rotation_mode="anchor")
    y_labels = [f"{y_ticks[i]:1.1e}" for i in y_ticks_location]
    _ = ax.set_yticklabels(y_labels, rotation=0)


# pylint: disable=too-many-arguments
def _plot_heat_hist(
    hist2d,
    x_name,
    y_name,
    x_ticks,
    y_ticks,
    norm,
    width,
    ticks_offset,
    filename,
    interactive=False,
    save=False,
    font_scale=None,
    title=None,
    p=None,
    ax=None,
    verbose=False,
):
    """

    Parameters
    ----------
    hist2d
    x_name
    y_name
    x_ticks
    y_ticks
    norm
    width
    ticks_offset
    filename
    interactive
    save
    font_scale
    title
    p
    ax
    verbose

    Returns
    -------

    """
    if title is None:
        if p is None:
            title = f"Histogram for {parse_word(y_name)} over {parse_word(x_name)}"
        else:
            p_name = r"$ \partial$" + parse_word(p)
            title = (
                f"Histogram for {p_name}/{parse_word(y_name)} over {parse_word(x_name)}"
            )

    if verbose:
        if p is None:
            print(f"Plotting histogram for {x_name}, {y_name}")
        else:
            print(f"Plotting histogram for {x_name}, d {p}/{y_name}")
    if ax is None:
        ax = sns.heatmap(
            np.transpose(
                hist2d
            ),  # np.histogram2d uses the first dimension for the x-axis
            cmap="viridis",
            cbar=True,
            square=True,
            norm=norm,
        )
    else:
        _ = sns.heatmap(
            np.transpose(
                hist2d
            ),  # np.histogram2d uses the first dimension for the x-axis
            cmap="viridis",
            cbar=True,
            square=True,
            norm=norm,
            ax=ax,
        )
    _set_heat_hist_ticks(hist2d, ax, ticks_offset, x_ticks, y_ticks, width)
    if font_scale is None:
        _ = ax.set_title(title)
    else:
        ax.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
        cbar = ax.collections[-1].colorbar
        cbarax = cbar.ax
        cbarax.tick_params(labelsize=int(10 * font_scale))
    ax.set_xlabel(parse_word(x_name), fontsize=int(11 * font_scale))
    ax.set_ylabel(parse_word(y_name), fontsize=int(11 * font_scale))

    plt.tight_layout()
    if filename is not None and save:
        save_plot(filename, ax.get_figure())
        if not interactive:
            plt.clf()
