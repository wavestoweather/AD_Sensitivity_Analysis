"""Module to plot 2D maps with longitude and latitude.

The functions may be used with the interactive routines.
"""
import holoviews as hv
from matplotlib.colors import SymLogNorm, LogNorm
from matplotlib.figure import Figure
import numpy as np
import seaborn as sns

from ad_sensitivity_analysis.plot import colors as color_funcs
from ad_sensitivity_analysis.plot import aux_ds_functions, aux_functions, latexify


# pylint: disable=too-many-arguments,too-many-locals
def plot_2dmap(
    ds,
    out_param,
    in_param,
    kind_param,
    pressure,
    time,
    cmap,
    fix,
    fix_time,
    log_plot,
    width,
    height,
    lthresh,
    title,
    font_scale,
    save,
    latex,
    save_path,
    static,
    colorblind=True,
):
    """

    Parameters
    ----------
    ds
    coord
    out_param
    in_param
    kind_param
    pressure
    time
    cmap
    fix
    fix_time
    log_plot
    width
    height
    lthresh
    title
    font_scale
    save
    latex
    save_path
    static
    colorblind

    Returns
    -------

    """
    colors, _ = color_funcs.get_b8_colors(colorblind=colorblind)

    if "Output Parameter" in ds:
        coord = "Output Parameter"
    else:
        coord = "Output_Parameter"

    # pylint: disable=no-member
    fig, ax, ds_tmp, mini, maxi, mini2, maxi2, linthresh = prepare_2dplot(
        ds=ds,
        coord=coord,
        out_param=out_param,
        in_param=in_param,
        kind_param=kind_param,
        pressure=pressure,
        time=time,
        fix=fix,
        fix_time=fix_time,
        width=width,
        height=height,
        lthresh=lthresh,
        latex=latex,
        static=static,
    )
    if title == "":
        title = None

    if in_param == "Top_Parameter":
        ds_tmp.plot.pcolormesh(
            y="lat",
            x="lon",
            levels=np.arange(len(colors)),
            colors=colors,
            ax=ax,
        )
    else:
        if log_plot and mini != 0 and maxi != 0:
            if np.abs(mini) > maxi:
                maxi = np.abs(mini)
            else:
                mini = -maxi
            # pylint: disable=unexpected-keyword-arg
            ds_tmp.plot(
                y="lat",
                x="lon",
                cmap=cmap,
                norm=SymLogNorm(
                    linthresh=linthresh,
                    vmin=mini,
                    vmax=maxi,
                    base=10,
                ),
                ax=ax,
            )
        elif log_plot:
            if linthresh != 0:
                mini2 = 10**lthresh
            ds_tmp.plot(
                y="lat",
                x="lon",
                cmap=cmap,
                norm=LogNorm(
                    vmin=mini2,
                    vmax=maxi2,
                ),
                ax=ax,
            )
        else:
            ds_tmp.plot(
                y="lat",
                x="lon",
                cmap=cmap,
                vmin=mini,
                vmax=maxi,
                ax=ax,
            )

    set_2dmap_labels(
        ax=ax,
        font_scale=font_scale,
        in_param=in_param,
        kind_param=kind_param,
        title=title,
    )
    if in_param == "Top_Parameter":
        color_funcs.set_top_param_cbar(
            fig=fig,
            ax=ax,
            font_scale=font_scale,
            colorblind=colorblind,
        )
    if save:
        # The local variables here have to be set in case
        # of interactive plotting.
        save_path, save = aux_functions.save_plot(save_path, ax)
    return fig


def prepare_2dplot(
    ds,
    coord,
    out_param,
    in_param,
    kind_param,
    pressure,
    time,
    fix,
    fix_time,
    width,
    height,
    lthresh,
    latex,
    static,
):
    """

    Parameters
    ----------
    ds
    coord
    out_param
    in_param
    kind_param
    pressure
    time
    fix
    fix_time
    width
    height
    lthresh
    latex
    static

    Returns
    -------

    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})

    ds = aux_ds_functions.set_regrid_missing_kind(
        ds=ds,
        kind_param=kind_param,
        in_param=in_param,
    )
    ds_local, ds_fix = aux_ds_functions.get_regrid_slice(
        ds=ds,
        in_param=in_param,
        out_param=out_param,
        kind_param=kind_param,
        time=time,
        coord=coord,
        pressure=pressure,
        fix=fix,
        fix_time=fix_time,
    )

    fig = Figure()
    # pylint: disable=no-member
    ax = fig.subplots()
    if in_param == "Top_Parameter":
        return fig, ax, ds_local, 0, 0, 0, 0, 0
    mini, maxi, mini2, maxi2 = aux_ds_functions.get_extreme_vals(ds_fix)
    linthresh = aux_ds_functions.get_log_threshold(da=ds_fix, lthresh=lthresh)
    min_local, max_local, min_local2, max_local2 = aux_ds_functions.get_extreme_vals(
        ds_local
    )
    static.value = (
        f"({mini:.2e}, {maxi:.2e}); "
        f"at {pressure/100} hPa: ({min_local:.2e}, {max_local:.2e}) -- "
        f"({min_local2:.2e}, {max_local2:.2e}) {linthresh:.2e}, {lthresh:.2e}"
    )
    return fig, ax, ds_local, mini, maxi, mini2, maxi2, linthresh


def set_2dmap_labels(ax, font_scale, in_param, kind_param, title):
    """

    Parameters
    ----------
    ax
    font_scale
    in_param
    kind_param
    title

    Returns
    -------

    """
    _ = ax.set_title(title, fontsize=int(12 * font_scale))
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=int(10 * font_scale),
    )
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.xaxis.get_label().set_fontsize(int(11 * font_scale))
    ax.yaxis.get_label().set_fontsize(int(11 * font_scale))
    ax.yaxis.grid(True, which="major")
    ax.xaxis.grid(True, which="major")
    cbar = ax.collections[-1].colorbar
    cbarax = cbar.ax
    if in_param == "Top_Parameter":
        cbarlabel = "Index of top parameter"
        cbar_ticks = [4, 14, 23, 30, 36]
        cbar_labels = [
            "evaporation",
            "CCN activation",
            "geometry",
            "fall velocity",
            "miscellaneous",
        ]
        cbar.set_ticks(cbar_ticks)
        cbar.set_ticklabels(cbar_labels)
    else:
        cbarlabel = f"{kind_param} " + latexify.parse_word(in_param).replace(
            r"\partial", ""
        )
    cbarax.tick_params(labelsize=int(10 * font_scale))
    cbar.set_label(
        label=cbarlabel,
        fontsize=int(11 * font_scale),
    )


def plot_heatmap(ds, col, fig_size=250, aspect=1, cmap="viridis"):
    """
    An easy and fast way to plot a given column along latitude and
    longitude.

    Parameters
    ----------
    ds
    col
    fig_size
    aspect
    cmap

    Returns
    -------
    holoviews.Image
    """
    x = "longitude"
    y = "latitude"
    return hv.Image((ds[x], ds[y], ds[col]), datatype=["grid"]).opts(
        fig_size=fig_size,
        aspect=aspect,
        cmap=cmap,
        colorbar=True,
        xlabel=x,
        ylabel=y,
        clabel=col,
    )
