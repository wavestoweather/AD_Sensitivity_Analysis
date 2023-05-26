"""Plot statistics over all trajectories.

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import seaborn as sns


from ad_sensitivity_analysis.data_handler import filter_data
from ad_sensitivity_analysis.plot import latexify
from ad_sensitivity_analysis.plot import aux_functions


def prepare_plot(width, height, latex, out_param, title=None):
    """

    Parameters
    ----------
    width
    height
    latex
    out_param
    title

    Returns
    -------

    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    fig = Figure()
    # pylint: disable=no-member
    ax = fig.subplots()
    if title is None:
        title = f"Histogram for impacts on {latexify.parse_word(out_param)}"
    return fig, ax, title


def set_kde_plot_labels(ax, font_scale, title, new_labels):
    """

    Parameters
    ----------
    ax
    font_scale
    title
    new_labels

    Returns
    -------

    """
    if new_labels is not None:
        ax.get_legend().remove()
        handles, _ = ax.get_legend_handles_labels()
        ax.legend(handles, new_labels, fontsize=int(9 * font_scale))
    else:
        plt.setp(ax.get_legend().get_texts(), fontsize=int(9 * font_scale))
        plt.setp(ax.get_legend().get_title(), fontsize=int(10 * font_scale))
    if font_scale is None:
        _ = ax.set_title(title)
    else:
        ax.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
    ax.set_xlabel("Rank", fontsize=int(11 * font_scale))
    ax.set_ylabel("Density", fontsize=int(11 * font_scale))
    plt.tight_layout()


def plot_kde_phase_flow(
    df, ax, common_norm, phase_colors, linewidth, bw_adjust, worst_rank
):
    """

    Parameters
    ----------
    df
    ax
    common_norm
    phase_colors
    linewidth
    bw_adjust
    worst_rank

    Returns
    -------

    """
    new_labels = []
    phase_hue_order = np.sort(np.unique(df["phase"]))
    df_tmp = df.loc[df["flow"] == "inflow"]
    _ = sns.kdeplot(
        data=df_tmp,
        x="Rank",
        hue="phase",
        hue_order=phase_hue_order,
        common_norm=common_norm,
        linestyle="dotted",
        palette=phase_colors,
        ax=ax,
        label="inflow",
        linewidth=linewidth,
        bw_adjust=bw_adjust,
        clip=(1, worst_rank),
    )
    for p in phase_hue_order[::-1]:
        if (
            p in df_tmp["phase"].values
            and len(np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])) >= 1
        ):
            new_labels.append(f"inflow {p}")
        # in case there is anything with zero variance, we plot a single bar.
        vals = np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])
        if len(vals) == 1:
            ax.axvline(
                x=vals[0],
                color=phase_colors[p],
                linestyle="dotted",
                label=f"inflow {p}",
            )
    df_tmp = df.loc[df["flow"] == "ascent"]
    _ = sns.kdeplot(
        data=df_tmp,
        x="Rank",
        hue="phase",
        hue_order=phase_hue_order,
        common_norm=common_norm,
        linestyle="solid",
        palette=phase_colors,
        ax=ax,
        label="ascent",
        linewidth=linewidth,
        bw_adjust=bw_adjust,
        clip=(1, worst_rank),
    )
    for p in phase_hue_order[::-1]:
        if (
            p in df_tmp["phase"].values
            and len(np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])) >= 1
        ):
            new_labels.append(f"ascent {p}")
        # in case there is anything with zero variance, we plot a single bar.
        vals = np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])
        if len(vals) == 1:
            ax.axvline(
                x=vals[0],
                color=phase_colors[p],
                linestyle="solid",
                label=f"ascent {p}",
            )
    df_tmp = df.loc[df["flow"] == "outflow"]
    _ = sns.kdeplot(
        data=df_tmp,
        x="Rank",
        hue="phase",
        hue_order=phase_hue_order,
        common_norm=common_norm,
        linestyle="dashed",
        palette=phase_colors,
        ax=ax,
        label="outflow",
        linewidth=linewidth,
        bw_adjust=bw_adjust,
        clip=(1, worst_rank),
    )
    for p in phase_hue_order[::-1]:
        if (
            p in df_tmp["phase"].values
            and len(np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])) >= 1
        ):
            new_labels.append(f"outflow {p}")
        # in case there is anything with zero variance, we plot a single bar.
        vals = np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])
        if len(vals) == 1:
            ax.axvline(
                x=vals[0],
                color=phase_colors[p],
                linestyle="dashed",
                label=f"outflow {p}",
            )
    return new_labels


def plot_kde_phase(df, ax, common_norm, phase_colors, linewidth, bw_adjust, worst_rank):
    """

    Parameters
    ----------
    df
    ax
    common_norm
    phase_colors
    linewidth
    bw_adjust
    worst_rank

    Returns
    -------

    """
    phase_hue_order = np.sort(np.unique(df["phase"]))
    _ = sns.kdeplot(
        data=df,
        x="Rank",
        hue="phase",
        hue_order=phase_hue_order,
        common_norm=common_norm,
        linestyle="solid",
        palette=phase_colors,
        ax=ax,
        linewidth=linewidth,
        bw_adjust=bw_adjust,
        clip=(1, worst_rank),
    )
    # in case there is anything with zero variance, we plot a single bar.
    for phase in np.unique(df["phase"]):
        vals = np.unique(df.loc[df["phase"] == phase]["Rank"])
        if len(vals) == 1:
            ax.axvline(x=vals[0], color=phase_colors[phase])


def plot_kde_flow(df, ax, common_norm, linewidth, bw_adjust, worst_rank):
    """

    Parameters
    ----------
    df
    ax
    common_norm
    linewidth
    bw_adjust
    worst_rank

    Returns
    -------

    """
    param_hue_order = np.sort(np.unique(df["Parameter"]))
    new_labels = []
    df_tmp = df.loc[df["flow"] == "inflow"]
    _ = sns.kdeplot(
        data=df_tmp,
        x="Rank",
        hue="Parameter",
        hue_order=param_hue_order,
        common_norm=common_norm,
        linestyle="dotted",
        ax=ax,
        label="inflow",
        linewidth=linewidth,
        bw_adjust=bw_adjust,
        clip=(1, worst_rank),
    )
    for param in param_hue_order[::-1]:
        if (
            param in df_tmp["Parameter"].values
            and len(np.unique(df_tmp.loc[df_tmp["Parameter"] == param]["Rank"])) > 1
        ):
            new_labels.append(f"inflow {latexify.parse_word(param[:-5])} rank")
    df_tmp = df.loc[df["flow"] == "ascent"]
    _ = sns.kdeplot(
        data=df_tmp,
        x="Rank",
        hue="Parameter",
        hue_order=param_hue_order,
        common_norm=common_norm,
        linestyle="solid",
        ax=ax,
        label="ascent",
        linewidth=linewidth,
        bw_adjust=bw_adjust,
        clip=(1, worst_rank),
    )
    for param in param_hue_order[::-1]:
        if (
            param in df_tmp["Parameter"].values
            and len(np.unique(df_tmp.loc[df_tmp["Parameter"] == param]["Rank"])) > 1
        ):
            new_labels.append(f"ascent {latexify.parse_word(param[:-5])} rank")
    df_tmp = df.loc[df["flow"] == "outflow"]
    _ = sns.kdeplot(
        data=df_tmp,
        x="Rank",
        hue="Parameter",
        hue_order=param_hue_order,
        common_norm=common_norm,
        linestyle="dashed",
        ax=ax,
        label="outflow",
        linewidth=linewidth,
        bw_adjust=bw_adjust,
        clip=(1, worst_rank),
    )
    for param in param_hue_order[::-1]:
        if (
            param in df_tmp["Parameter"].values
            and len(np.unique(df_tmp.loc[df_tmp["Parameter"] == param]["Rank"])) > 1
        ):
            new_labels.append(f"outflow {latexify.parse_word(param[:-5])} rank")
    return new_labels


def plot_kde_params(df, ax, common_norm, linewidth, bw_adjust, worst_rank):
    """

    Parameters
    ----------
    df
    ax
    common_norm
    linewidth
    bw_adjust
    worst_rank

    Returns
    -------

    """
    param_hue_order_tmp = np.sort(np.unique(df["Parameter"]))
    param_hue_order = []
    for p in param_hue_order_tmp:
        vals = np.unique(df.loc[df["Parameter"] == p]["Rank"])
        if len(vals) == 1 and (0 in vals or worst_rank in vals):
            continue
        if len(vals) == 2 and (0 in vals and worst_rank in vals):
            continue
        param_hue_order.append(p)
    df = df.loc[df["Parameter"].isin(param_hue_order)]
    _ = sns.kdeplot(
        data=df,
        x="Rank",
        hue="Parameter",
        hue_order=param_hue_order,
        common_norm=common_norm,
        linestyle="solid",
        ax=ax,
        linewidth=linewidth,
        bw_adjust=bw_adjust,
        clip=(1, worst_rank),
    )
    new_labels = []
    for l in param_hue_order[::-1]:
        new_labels.append(latexify.parse_word(l[:-5]))
        new_labels[-1] += " rank"
    return new_labels


# pylint: disable=too-many-arguments,too-many-locals
def plot_kde_histogram(
    ds,
    in_params,
    out_param,
    flow,
    phase,
    ignore_zero_gradients,
    linewidth=2,
    bw_adjust=1.0,
    title=None,
    filename=None,
    width=17,
    height=16,
    common_norm=False,
    font_scale=None,
    save=True,
    latex=False,
):
    """
    Plot kde histogram of ranks.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset created via create_rank_traj_dataset().
    in_params : list of strings
        A list of model parameters.
    out_param : string
        Model state variable which the sensitivities are for.
    flow : bool
        Plot different flow phases, i.e., inflow, ascent, and outflow.
    phase : bool
        Plot different phases, i.e., warm phase, mixed phase, and ice phase.
    ignore_zero_gradients : bool
        Do not show the rank of gradients that are always zero.
    linewidth : float
        Line width of each kde.
    bw_adjust : float
        Used to calculate the kde. Adjusts multiplicatively the bandwidth that is found automatically by
        scipy.stats.gaussian.kde. Larger values generate smoother estimations.
    title : string
        Title of the histogram. If none is given, a title will be generated.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    width : float
        Width in inches
    height : float
        Height in inches
    common_norm : bool
        All densities sum up to one if true
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    latex : bool
        Use latex font.

    Returns
    -------
    matplotlib.figure.Figure with the plot drawn onto it.
    """
    phase_colors = {
        "warm phase": "tab:orange",
        "mixed phase": "tab:green",
        "ice phase": "tab:blue",
        "any": "k",
    }

    df, worst_rank = filter_data.filter_rank_data(
        ds=ds,
        out_param=out_param,
        in_params=in_params,
        ignore_zero_gradients=ignore_zero_gradients,
        phase=phase,
        flow=flow,
    )

    fig, ax, title = prepare_plot(
        width=width,
        height=height,
        latex=latex,
        out_param=out_param,
        title=title,
    )
    new_labels = None
    if phase and flow:
        new_labels = plot_kde_phase_flow(
            df=df,
            ax=ax,
            common_norm=common_norm,
            phase_colors=phase_colors,
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            worst_rank=worst_rank,
        )
    elif phase:
        if df.empty:
            return fig
        plot_kde_phase(
            df=df,
            ax=ax,
            common_norm=common_norm,
            phase_colors=phase_colors,
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            worst_rank=worst_rank,
        )
    elif flow:
        new_labels = plot_kde_flow(
            df=df,
            ax=ax,
            common_norm=common_norm,
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            worst_rank=worst_rank,
        )
    else:
        plot_kde_params(
            df=df,
            ax=ax,
            common_norm=common_norm,
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            worst_rank=worst_rank,
        )

    set_kde_plot_labels(
        ax=ax,
        font_scale=font_scale,
        title=title,
        new_labels=new_labels,
    )

    if filename is not None and save:
        aux_functions.save_plot(filename, fig)
    return fig
