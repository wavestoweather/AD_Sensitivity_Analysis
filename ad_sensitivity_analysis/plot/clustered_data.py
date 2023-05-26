"""Plot clustered data.

"""
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from ad_sensitivity_analysis.plot.aux_functions import get_save_name
from ad_sensitivity_analysis.plot.latexify import parse_word
from ad_sensitivity_analysis.plot.latexify import mappings as latex_mappings


def parse_axis_name(in_p, out_p, reduce_name):
    """


    Parameters
    ----------
    in_p
    out_p
    reduce_name

    Returns
    -------
    Two strings, one for the
    """
    if in_p[0] == "d" and in_p != "deposition":
        if " " in in_p:
            in_p_parse = in_p.split()[0]
            reduce_name = in_p.split()[1] + " "
        else:
            in_p_parse = in_p
        data_name = f"{reduce_name}d{out_p}/{in_p}"
        label = (
            reduce_name
            + r"$\partial$"
            + parse_word(out_p)
            + f"/{parse_word(in_p_parse)}"
        )
    elif in_p in latex_mappings:
        data_name = f"{reduce_name}{in_p}"
        label = f"{reduce_name}{parse_word(in_p)}"
    else:
        data_name = in_p
        label = reduce_name + in_p
    return data_name, label


# pylint: disable=too-many-arguments, too-many-locals
def plot_cluster_data(
    data,
    in_p_x,
    out_p_x,
    in_p_y,
    out_p_y,
    reduce_name,
    logx,
    logy,
    width,
    height,
    font_scale,
    title,
    save_path,
    latex,
    save,
    dot_size,
):
    """
    Plot clustered data as either histogram or scatterplot.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame generated using get_cluster().
    in_p_x
    out_p_x
    in_p_y
    out_p_y
    reduce_name : string
        Name prepended to the columns. Should relate to the reduction
        applied to the dataset, such as "avg" or "rank" in get_cluster().
    logx
    logy
    width
    height
    font_scale
    title
    save_path
    latex
    save
    dot_size

    Returns
    -------
    matplotlib.axes.Axes created using seaborn plot function.
    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    fig = Figure()
    # pylint: disable=no-member
    ax = fig.subplots()
    x, x_label = parse_axis_name(in_p_x, out_p_x, reduce_name)
    y, y_label = parse_axis_name(in_p_y, out_p_y, reduce_name)
    histogram = (in_p_x == in_p_y and out_p_x == out_p_y) or (
        in_p_x == in_p_y and "/" not in x
    )
    if histogram:
        palette = "tab10"
        if len(set(data["cluster"])) > 10:
            palette = "tab20"
        sns.histplot(
            data=data,
            x=x,
            hue="cluster",
            palette=palette,
            multiple="stack",
            bins=100,
            log_scale=(logx, logy),
            ax=ax,
        )
    else:
        sns.scatterplot(
            data=data,
            x=x,
            y=y,
            hue="cluster",
            palette="tab10",
            s=dot_size,
            ax=ax,
        )
    if logx and not histogram:
        if np.nanmin(data[x]) < 0:
            linthresh = np.nanmin(np.abs(data[x].where(data[x] != 0)))
            ax.set_xscale("symlog", linthresh=linthresh)
        else:
            ax.set_xscale("log")
    if logy and not histogram:
        if np.nanmin(data[y]) < 0:
            linthresh = np.nanmin(np.abs(data[y].where(data[y] != 0)))
            ax.set_yscale("symlog", linthresh=linthresh)
        else:
            ax.set_yscale("log")
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=int(10 * font_scale),
    )
    _ = ax.set_title(title, fontsize=int(12 * font_scale))
    ax.set_xlabel(x_label, fontsize=int(11 * font_scale))
    ax.set_ylabel(y_label, fontsize=int(11 * font_scale))
    legend = ax.get_legend()
    legend.set_title("cluster", prop={"size": int(11 * font_scale)})
    plt.setp(legend.get_texts(), fontsize=int(10 * font_scale))
    ax.yaxis.get_offset_text().set_fontsize(int(11 * font_scale))
    ax.xaxis.get_offset_text().set_fontsize(int(11 * font_scale))
    # You may use the following line to remove the offset label if needed.
    # ax.xaxis.get_offset_text().set(alpha=0)
    if save:
        try:
            save_name = get_save_name(save_path)
            ax.figure.savefig(save_name, bbox_inches="tight", dpi=300)
        except IOError:
            save_path = f"Could not save to {save_path}. Did you forget the filetype?"
        save_path = None
        save = False
    return fig
