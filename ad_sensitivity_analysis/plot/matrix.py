"""Plot heatmaps from covariance or correlation matrix.

"""
import copy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import seaborn as sns

from ad_sensitivity_analysis.plot.aux_functions import save_plot


def _get_data_names(data, names, sort):
    if sort:
        pairs = []
        i = 0
        for name, col in zip(names, data):
            pairs.append((np.nansum(np.abs(col)), name, i))
            i += 1
        pairs.sort(key=lambda x: x[0])
        permutation = []
        for p in pairs[::-1]:
            permutation.append(p[2])
        plot_data = data[:, permutation]
        plot_data = plot_data[permutation, :]
        plot_names = np.asarray(names)[permutation]
    else:
        plot_data = data
        plot_names = names
    return plot_data, plot_names


# pylint: disable=too-many-arguments, too-many-locals
def plot_heatmap_matrix(
    data_in,
    out_param,
    names,
    flow=None,
    phase=None,
    latex=False,
    font_scale=1.0,
    save=False,
    in_params=None,
    plot_type="all",
    norm=None,
    title=None,
    filename=None,
    cmap="viridis",
    width=30,
    height=16,
    sort=False,
    annot=False,
):
    """
    Plot a heatmap of a covariance or correlation matrix.

    Parameters
    ----------
    data_in : dictionary of matrices or dictionary of dictionary of dictionary of matrices
        Either: The keys are model state variables with available sensitivities.
        The values are matrices where covariances are given for sensitivities for the respective model state variable.
        The data can be generated using get_cov_matrix().
        Or: The first keys are model state variables, then the flow, then phase used for calculating
        correlation/covariance values. The data can be generated using get_stats_per_traj.get_matrix().
    out_param : string
        The model state variable for which the sensitivities shall be plotted for.
    names : list of strings
        The names of each column/row in the covariance matrix. The index of names is the row/column of the matrix.
    in_params : list-like of strings
        List of model parameters or model states to plot in the covariance matrix.
    plot_type : string
        Define if only negative (='negative'), positive (='positive'), or all ('all') values shall be plotted.
    norm : matplotlib.colors normalization instance
        A normalization for the colormap, such as matplotlib.colors.SymLogNorm()
    title : string
        Optional title for the plot. Otherwise a standard name will be used.
    filename : string
        Path and filename to save the plot on disk. The filename will be numerated.
        If a file with the same name already exists, the number will be incremented.
    cmap : string
        Name of the color palette passed to seaborn.
    width : float
        Width in inches
    height : float
        Height in inches
    sort : bool or list of strings
        If true: sort the rows and columns such that the one with the largest sum is on top/left and the smallest
        sum is on bottom/right.
        If list

    Returns
    -------
    If successful, returns the matplotlib.figure.Figure with the plot drawn onto it. Otherwise returns None.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    fig = Figure()
    if isinstance(data_in[out_param], dict):
        data = copy.deepcopy(data_in[out_param][flow][phase])
    else:
        data = copy.deepcopy(data_in[out_param])
    if in_params is not None and len(in_params) < np.shape(data)[0]:
        if len(in_params) == 0:
            return fig
        idx = []
        names = np.asarray(list(names))
        for in_p in in_params:
            idx.append(np.where(names == in_p)[0][0])
        idx = np.asarray(idx)
        data = data[idx]
        data = data[:, idx]
        names = names[idx]
    if plot_type == "negative":
        data[np.where(data >= 0)] = np.nan
    elif plot_type == "positive":
        data[np.where(data < 0)] = np.nan

    plot_data, plot_names = _get_data_names(data, names, sort)

    # pylint: disable=no-member
    ax = fig.subplots()
    sns.heatmap(
        data=plot_data,
        cmap=cmap,
        norm=norm,
        cbar=True,
        yticklabels=plot_names,
        xticklabels=plot_names,
        annot=annot,
        fmt="1.2e",
        ax=ax,
        annot_kws={"size": int(9 * font_scale)},
    )
    if title is None:
        title = f"Heatmap of Relation (Gradients for {out_param})"
    _ = ax.set_title(title, fontsize=int(12 * font_scale))
    _ = ax.set_yticklabels(
        ax.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor"
    )
    _ = ax.set_xticklabels(
        ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor"
    )
    ax.tick_params(axis="x", which="major", labelsize=int(10 * font_scale))
    ax.tick_params(axis="y", which="major", labelsize=int(10 * font_scale))
    cbar = ax.collections[-1].colorbar
    cbarax = cbar.ax
    cbarax.tick_params(labelsize=int(10 * font_scale))

    plt.tight_layout(h_pad=1)
    if filename is not None and save:
        save_plot(filename, fig)
    return fig
