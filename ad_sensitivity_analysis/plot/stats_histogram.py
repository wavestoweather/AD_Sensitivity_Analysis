"""Histogram plots with statistics from sensitivity simulations.

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_col
from matplotlib.figure import Figure

from tqdm.auto import tqdm
import seaborn as sns

from ad_sensitivity_analysis.plot.aux_functions import save_plot
from ad_sensitivity_analysis.plot.latexify import parse_word
from ad_sensitivity_analysis.plot.stats_histogram_aux import (
    _plot_heat_hist,
    _plot_hist,
    _get_kde_df,
    _get_params_ordered,
    _get_corr_matrix_array,
)


# pylint: disable=too-many-arguments, too-many-locals
def traj_plot_histogram_out(
    out_params,
    filename,
    edges,
    hist,
    log=False,
    width=24,
    height=12,
    font_scale=None,
    title=None,
    save=True,
    interactive=False,
    latex=False,
    ticks_offset=None,
    verbose=False,
):
    """
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot the histogram of an output parameter,
    i.e., QV, latent_heat, latent_cool, etc.

    Parameters
    ----------
    out_params : string or list-like of strings
        Output parameter name or multiple output parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges : Dictionary of list-like of float
        Edges for the histogram. Keys must be in out_params.
    hist : Dictionary of list-like of int
        Number of entries for each bin. Keys must be in out_params.
    log : bool
        Plot the y-axis using log-scale.
    width : float
        Width in inches
    height : float
        Height in inches
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    interactive : bool
        Create a figure for interactive plotting.
    title : string
        Title of the histogram. If none is given, a title will be generated.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.
    ticks_offset : int
        If none, the number of ticks is calculated based on the width of the plot. Otherwise,
        every other "ticks_offset" is used.
    verbose : bool
        Print some additional information.

    Returns
    -------
    If interactive == False, returns None. If interactive == True, returns the matplotlib.figure.Figure
    with the plot drawn onto it.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    if interactive:
        sns.set()
        fig = Figure()
        # pylint: disable=no-member
        ax = fig.subplots()
    else:
        fig = None
        ax = None

    def plot_hist(out_p, title=None, ax=None):
        if title is None:
            title = f"Histogram for {parse_word(out_p)}"
        if verbose:
            print(f"Plotting histogram for {out_p}")

        if ticks_offset is None:
            if 24 > width > 5:
                local_ticks_offset = 24 // (width - 5)
            elif width <= 5:
                local_ticks_offset = 6
        else:
            local_ticks_offset = ticks_offset
        if ticks_offset >= len(edges[out_p]) - 1:
            local_ticks_offset = len(edges[out_p]) - 2
        if ax is None:
            ax = sns.barplot(
                x=edges[out_p][:-1], y=hist[out_p], log=log, color="seagreen"
            )
        else:
            sns.barplot(
                x=edges[out_p][:-1], y=hist[out_p], log=log, color="seagreen", ax=ax
            )
        ax.set(xticks=np.arange(0, len(edges[out_p]) - 1, local_ticks_offset))
        _ = ax.set_xticklabels(
            [f"{tick:1.1e}" for tick in edges[out_p][:-1:local_ticks_offset]],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        if font_scale is None:
            _ = ax.set_title(title)
            ax.set_ylabel(out_p)
        else:
            ax.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
            _ = ax.set_title(title, fontsize=int(12 * font_scale))
            ax.set_ylabel("# Entries", fontsize=int(11 * font_scale))
            ax.set_xlabel(f"Bins of {parse_word(out_p)}", fontsize=int(11 * font_scale))

        if filename is not None and save:
            plt.tight_layout()
            if not interactive:
                save_plot(filename, ax.get_figure())
                plt.clf()
            else:
                save_plot(filename, fig)

    if isinstance(out_params, list):
        for out_p in tqdm(out_params):
            plot_hist(out_p, title, None)
    else:
        plot_hist(out_params, title, ax)
    if interactive:
        # Just a failsafe to avoid repeating error messages with erroneous filepaths.
        save = False
        filename = None
    return fig


def traj_plot_histogram_inp(
    in_params,
    filename,
    edges_in_params,
    hist_in_params,
    log=False,
    font_scale=None,
    width=24,
    height=12,
    title=None,
    save=True,
    interactive=False,
    latex=False,
    ticks_offset=None,
    verbose=False,
):
    r"""
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot three histograms per image with
    \partial output / \partial model parameter
    where output is QV, latent_heat and latent_cool. Plot one image per model_parameter.

    Parameters
    ----------
    in_params : string or list-like of strings
        Model parameter name or multiple input parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges_in_params : Dictionary of dictionary of list-like of float
        Edges for the histogram. First keys are output parameters, second keys must be in in_params.
    hist_in_params : Dictionary of dictionary of list-like of int
        Number of entries for each bin. First keys are output parameters, second keys must be in in_params.
    log : bool
        Plot the y-axis using log-scale.
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.
    ticks_offset : int
        If none, the number of ticks is calculated based on the width of the plot. Otherwise,
        every other "ticks_offset" is used.
    interactive : bool
        Create a figure for interactive plotting.

    Returns
    -------
    If interactive == False, returns None. If interactive == True, returns the matplotlib.figure.Figure with the plot
    drawn onto it.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    out_params = list(edges_in_params.keys())
    if interactive:
        sns.set()
        fig = Figure()
        # pylint: disable=no-member
        axs = fig.subplots(
            nrows=len(out_params),
            ncols=1,
        )
    else:
        fig = None
        axs = None

    if isinstance(in_params, list):
        for in_p in tqdm(in_params):
            _plot_hist(
                out_params=out_params,
                in_p=in_p,
                verbose=verbose,
                edges_in_params=edges_in_params,
                hist_in_params=hist_in_params,
                log=log,
                ticks_offset=ticks_offset,
                width=width,
                fig=fig,
                axs=axs,
                filename=filename,
                font_scale=font_scale,
                title=title,
                save=save,
                interactive=interactive,
            )
    else:
        _plot_hist(
            out_params=out_params,
            in_p=in_params,
            verbose=verbose,
            edges_in_params=edges_in_params,
            hist_in_params=hist_in_params,
            log=log,
            ticks_offset=ticks_offset,
            width=width,
            fig=fig,
            axs=axs,
            filename=filename,
            font_scale=font_scale,
            title=title,
            save=save,
            interactive=interactive,
        )

    if verbose:
        print("All plots finished!")
    if interactive:
        # Just a failsafe to avoid repeating error messages with erroneous filepaths.
        save = False
        filename = None
    return fig


def traj_plot_kde_inp(
    in_params,
    out_param,
    filename,
    edges_in_params,
    hist_in_params,
    filter_zero=False,
    linewidth=2,
    bw_adjust=0.1,
    log_weights=False,
    log_x=False,
    font_scale=1,
    width=24,
    height=12,
    title=None,
    save=False,
    latex=False,
):
    """
    Similar to traj_plot_histogram_inp() but using a kde estimation.

    Parameters
    ----------
    in_params : list-like of strings
        Input parameter names to plot the histogram for.
    out_param : string
        Name od model state variable to plot the sensitivities for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges_in_params : Dictionary of dictionary list-like of float
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        the bin edges for the given keys. Optional: The first level can be another dictionary of phases.
    hist_in_params : Dictionary of dictionary list-like of int
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        values of the histogram for the given key. Optional: The first level can be another dictionary of phases.
    filter_zero : bool
        Filter zero gradients. Does not really help other than there is a dent in the distribution.
    linewidth : float
        Line width of each kde.
    bw_adjust : float
        Used to calculate the kde. Adjusts multiplicatively the bandwidth that is found automatically by
        scipy.stats.gaussian.kde. Larger values generate smoother estimations. Since the underlying data
        here is a set of bins with their counts as weights, the automatically estimated bandwidth is
        smoother than what would be used when using the complete dataset. Therefore, we use a default of 0.1.
    log_weights : bool
        Plot the y-axis using log-scale by using log10 weights.
    log_x : bool
        Plot the x-axis using log-scale.
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.

    Returns
    -------
    Returns the matplotlib.figure.Figure with the plot drawn onto it.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    common_norm = False

    if len(in_params) == 0:
        return None

    df = _get_kde_df(
        hist_in_params=hist_in_params,
        edges_in_params=edges_in_params,
        out_param=out_param,
        in_params=in_params,
        log_weights=log_weights,
        filter_zero=filter_zero,
    )
    fig = Figure()
    # pylint: disable=no-member
    ax = fig.subplots()
    if log_x:
        signs = np.sign(df["Impact"])
        zeros = df["Impact"] == 0
        df["Impact"] = np.log10(np.abs(df["Impact"]))
        min_imp = np.nanmin(df["Impact"])  # * 1.1
        df["Impact"] -= min_imp
        df["Impact"] *= signs
        df["Impact"] = df["Impact"].where(~zeros, 0)

    _ = sns.kdeplot(
        data=df,
        x="Impact",
        hue="Parameter",
        weights="weight",
        common_norm=common_norm,
        ax=ax,
        bw_adjust=bw_adjust,
        hue_order=[parse_word(p) for p in np.sort(in_params)],
        linewidth=linewidth,
    )

    if title is not None:
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
    if log_x:
        old_x_ticks = ax.get_xticks()
        signs = np.sign(old_x_ticks)
        if latex:
            x_labels = []
            for i, x_tick in enumerate(old_x_ticks):
                l = "$"
                if signs[i] < 0:
                    l += "-"
                l += r"10^{" + f"{(abs(x_tick) + min_imp):1.2f}" + r"}$"
                x_labels.append(l)
        else:
            x_labels = [
                f"{signs[i] * 10**(abs(x_tick) + min_imp):1.1e}"
                for i, x_tick in enumerate(old_x_ticks)
            ]
        for i, tick_val in enumerate(old_x_ticks):
            if tick_val == 0:
                if latex:
                    x_labels[i] = r"$\leq 10^{" + f"{min_imp:1.2f}" + r"}$"
                else:
                    x_labels[i] = r"$\leq$" + f"{10**min_imp:1.1e}"
                break
        # For some reason, set_xticklabels() resets the font.
        _ = ax.set_xticklabels(x_labels, fontdict={"usetex": latex})

    ax.xaxis.get_label().set_fontsize(int(11 * font_scale))
    ax.yaxis.get_label().set_fontsize(int(11 * font_scale))
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=int(10 * font_scale),
    )
    plt.setp(ax.get_legend().get_texts(), fontsize=int(9 * font_scale))
    plt.setp(ax.get_legend().get_title(), fontsize=int(11 * font_scale))
    ax.xaxis.get_offset_text().set_size(int(9 * font_scale))
    ax.yaxis.get_offset_text().set_size(int(9 * font_scale))

    if filename is not None and save:
        save_plot(filename, ax.get_figure())
    return fig


def plot_heatmap_traj(
    in_params,
    filename,
    edges_in_params,
    hist_in_params,
    width=24,
    height=24,
    title=None,
    verbose=False,
):
    """
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot a heatmap "model parameters" over
    "bin number" where each row is another histogram for a certain model parameter.

    Parameters
    ----------
    in_params : string or list-like of strings
        Model parameter name or multiple input parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges_in_params : Dictionary of dictionary list-like of float
        Edges for the histogram. First keys are output parameters, second keys must be in in_params.
    hist_in_params : Dictionary of dictionary list-like of int
        Number of entries for each bin. First keys are output parameters, second keys must be in in_params.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    verbose : bool
        Additional output.
    """
    sns.set(rc={"figure.figsize": (width, height)})
    # sort the histogram by a simple similarity metric. It is not perfect but better than random.
    out_params = list(edges_in_params.keys())
    for out_p in tqdm(out_params):
        in_params_tmp = []
        for in_p in in_params:
            if in_p in hist_in_params[out_p]:
                in_params_tmp.append(in_p)

        corr_matrix = _get_corr_matrix_array(
            in_params=in_params,
            in_params_tmp=in_params_tmp,
            out_p=out_p,
            hist_in_params=hist_in_params,
            verbose=verbose,
        )
        if verbose:
            print("Sort parameters according to correlation matrix")
        in_params_sorted = _get_params_ordered(corr_matrix, in_params_tmp, verbose)
        # found the correct order. Let's create a corresponding 2D-array and plot it
        n_edges = len(edges_in_params[out_p][in_params_tmp[0]][:-1])
        plot_matrix = np.zeros((len(in_params_tmp), n_edges))
        for i, in_p in enumerate(in_params_sorted):
            plot_matrix[i, :] = hist_in_params[out_p][in_p]
        plot_matrix[plot_matrix == 0] = np.nan
        ax = sns.heatmap(
            plot_matrix,
            cbar=True,
            cmap="viridis",
            norm=mpl_col.LogNorm(),
            yticklabels=in_params_sorted,
        )

        if title is None:
            title2 = f"Histogram for {out_p}"
        else:
            title2 = title
        _ = ax.set_title(title2)
        plt.tight_layout()
        save_plot(filename, ax.get_figure())
        plt.clf()


# pylint: disable=too-many-arguments
def plot_heatmap_histogram(
    hist_conditional,
    in_params,
    out_params,
    conditions,
    title=None,
    filename=None,
    width=17,
    height=16,
    log=True,
    font_scale=None,
    save=True,
    interactive=False,
    latex=False,
    ticks_offset=None,
    verbose=False,
):
    """
    Plot 2D histogram.

    Parameters
    ----------
    hist_conditional : Dictionary with edges and bin counts
        The dictionary generated using get_histogram_cond(). It has the following keys:
        'edges_out_params': Dictionary where the keys are model state variables and the values are arrays of bin edges.
        'edges_in_params': Dictionary where the keys are model state variables for which sensitivities are available
            and the values are dictionaries of model parameters with arrays of bin edges.
        model state variables: Each model state variable has a dictionary for 'hist_out_params' and 'hist_in_params'.
        'hist_out_params' is a dictionary of model state variables with arrays of bin counts.
        'hist_in_params' is a dictionary of model state variables for which sensitivities are available
            and the values are dictionaries of model parameters with arrays of bin counts.
    in_params : list of strings
        A list of model parameters.
    out_params : list of strings
        A list of model state variables with available sensitivities.
    conditions : list of strings
        A list of model state variables for the x-axis.
    title : string
        Title of the histogram. If none is given, a title will be generated.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    width : float
        Width in inches
    height : float
        Height in inches
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    log : bool
        Plot the histograms using a log-scale
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    interactive : bool
        Create a figure for interactive plotting.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.
    ticks_offset : int
        If none, the number of ticks is calculated based on the width of the plot. Otherwise,
        every other "ticks_offset" is used.
    verbose : bool
        Print some additional information.

    Returns
    -------
    If filename is given, returns None. If filename is None, returns the matplotlib.figure.Figure with the plot drawn
     onto it.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    if interactive:
        sns.set()
        fig = Figure()
        # pylint: disable=no-member
        ax = fig.subplots()
    else:
        fig = None
        ax = None
    norm = None
    if log:
        norm = mpl_col.LogNorm()

    if interactive:
        c = conditions
        if in_params is None or in_params == "None":
            out_p = out_params
            _plot_heat_hist(
                hist2d=hist_conditional[c]["hist_out_params"][out_p],
                x_name=c,
                y_name=out_p,
                x_ticks=hist_conditional["edges_out_params"][c][:-1],
                y_ticks=hist_conditional["edges_out_params"][out_p][:-1],
                title=title,
                ax=ax,
                norm=norm,
                width=width,
                ticks_offset=ticks_offset,
                filename=filename,
                interactive=interactive,
                save=save,
                font_scale=font_scale,
                verbose=verbose,
            )
        else:
            p = out_params
            in_p = in_params
            _plot_heat_hist(
                hist2d=hist_conditional[c]["hist_in_params"][p][in_p],
                x_name=c,
                y_name=in_p,
                x_ticks=hist_conditional["edges_out_params"][c][:-1],
                y_ticks=hist_conditional["edges_in_params"][p][in_p][:-1],
                title=title,
                ax=ax,
                norm=norm,
                width=width,
                ticks_offset=ticks_offset,
                filename=filename,
                interactive=interactive,
                save=save,
                font_scale=font_scale,
                p=p,
                verbose=verbose,
            )
        save = False
        filename = None
        return fig

    for c in tqdm(conditions):
        for out_p in tqdm(out_params, leave=False):
            if c == out_p:
                continue
            _plot_heat_hist(
                hist2d=hist_conditional[c]["hist_out_params"][out_p],
                x_name=c,
                y_name=out_p,
                x_ticks=hist_conditional["edges_out_params"][c][:-1],
                y_ticks=hist_conditional["edges_out_params"][out_p][:-1],
                title=title,
                norm=norm,
                width=width,
                ticks_offset=ticks_offset,
                filename=filename,
                interactive=interactive,
                save=save,
                font_scale=font_scale,
                verbose=verbose,
            )
        # Sensitivities are p w.r.t. in_p
        for p in tqdm(list(hist_conditional["edges_in_params"]), leave=False):
            for in_p in tqdm(in_params, leave=False):
                if in_p not in hist_conditional[c]["hist_in_params"][p]:
                    continue
                _plot_heat_hist(
                    hist2d=hist_conditional[c]["hist_in_params"][p][in_p],
                    x_name=c,
                    y_name=in_p,
                    x_ticks=hist_conditional["edges_out_params"][c][:-1],
                    y_ticks=hist_conditional["edges_in_params"][p][in_p][:-1],
                    title=title,
                    norm=norm,
                    width=width,
                    ticks_offset=ticks_offset,
                    filename=filename,
                    interactive=interactive,
                    save=save,
                    font_scale=font_scale,
                    p=p,
                    verbose=verbose,
                )

    if verbose:
        print("All plots finished!")
    return None
