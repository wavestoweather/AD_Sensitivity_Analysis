try:
    from Deriv_dask import Deriv_dask
except:
    from scripts.Deriv_dask import Deriv_dask
try:
    import latexify
except:
    import scripts.latexify as latexify

from cartopy import crs
from holoviews import opts
from holoviews.operation.datashader import datashade as dsshade
import holoviews as hv
import hvplot.pandas
import matplotlib
import numpy as np
import os
import pandas as pd
import xarray as xr


def add_liquid_content(data):
    q_total = None
    for col in data.columns:
        if "QR" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
        elif "QC" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
    data["Q_liquid"] = q_total
    return data


def add_cold_content(data):
    q_total = None
    for col in data.columns:
        if "QI" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
        elif "QS" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
        elif "QG" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
        elif "QH" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
    data["Q_cold"] = q_total
    return data


def add_sepcific_humidity(data):
    data["Specific_Humidity"] = data["QV"] / (data["QV"] + 1)
    return data


def load_data(data_path, traj_type, verbosity=0):
    """
    Load data, get rotation for pole and cut everything of that is not
    between -2800 s and 26000 s of the ascent.
    Change pressure from Pa to hPa and time_after_ascent s to h.
    Adds total liquid and cold content.

    Parameters
    ----------
    data_path : string
        Path to multiple data files in NetCDF-format.
    traj_type : string
        Type of trajectory to load, e.g. "conv" for all convective trajectories,
        "conv_400" for convective 400 hPa ascend trajectories.
    verbosity : int
        Set verbosity level.
        0: No output except for exceptions
        1: Print loading statements
        2: Print printing statements (not in this function)
        3: Print additional statements

    Returns
    -------
    pandas.Dataframe, pole longitude and pole latitude rotation
    """
    file_list = []
    for f in os.listdir(data_path):
        if traj_type in f:
            file_list.append(os.path.join(data_path, f))
    file_list = np.sort(np.asarray(file_list))
    df = None
    pollon = None
    pollat = None
    for f in file_list:
        if verbosity > 0:
            print(f"Loading {f}")
        ds = xr.open_dataset(f, decode_times=False, engine="netcdf4")
        ds = ds.where(
            ((ds["time_after_ascent"] >= -2800) & (ds["time_after_ascent"] <= 26000)),
            drop=True,
        )
        if df is not None:
            df = df.append(ds.to_dataframe())
        else:
            pollon = ds.attrs["pollon"]
            pollat = ds.attrs["pollat"]
            if verbosity > 2:
                print(f"Got lon={pollon}, lat={pollat} for rotated pole.")
            df = ds.to_dataframe()
    # drop all the NaN entries now that coordinates are collapsed to columns
    df = df.loc[(df["time_after_ascent"] >= -2800) & (df["time_after_ascent"] < 26000)]
    df["time_after_ascent_h"] = df["time_after_ascent"] / 3600
    df["pressure_hPa"] = df["pressure"] / 100
    if verbosity > 2:
        print("Adding liquid contents")
    df = add_liquid_content(df)
    if verbosity > 2:
        print("Adding cold contents")
    df = add_cold_content(df)
    if verbosity > 2:
        print("Adding specific humidity")
    df = add_sepcific_humidity(df)
    if verbosity > 2:
        print("Finished loading data")
    return df, pollon, pollat


def plot_line(
    df,
    y_axis,
    store_path,
    width=1959,
    height=1224,
    formatter_limits=None,
    datashade=True,
    s=None,
    title="",
    show_legend=False,
    backend="bokeh",
    verbosity=0,
):
    """
    Create a line plot with a shaded area for min and max values
    and save it with "time_after_ascent"
    on the x-axis.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe loaded with load_data(..).
    y_axis : List of string
        The name of the column for the y-axis for each plot.
    store_path : string
        Path (and name) where to save images. y_axis will be appended.
    formatter_limits : tuple of ints
            Lower and upper limits for formatting x- and y-axis.
    width : int
        The width of the plot (bokeh) or the relative width of the plot
        (matplotlib).
    height : int
        The height of the plot (bokeh) or the relative height of the plot
        (matplotlib).
    datashade : boolean
        Use datashade for the dataset except for those with type "Quantile".
    s : int
        Size of the dots.
    title : string
        Title for this plot.
    show_legend : boolean
        Wether to plot the legend or not.
    backend : string
        Choose a backend for plotting. Options are:
        matplotlib: Most plots should be fine with it.
        bokeh: Recommended.
    verbosity : int
        Set verbosity level.
        0: No output except for exceptions
        1: Print loading statements (not in this function)
        2: Print printing statements
        3: Print additional statements

    Returns
    -------
    Holoviews image.
    """

    x_axis = "time_after_ascent_h"
    plot_obj = Deriv_dask(
        direc="",
        parquet=False,
        netcdf=True,
        columns=None,
        backend=backend,
        file_ending="",
    )
    by = "type"
    lim_multiplier = 0.1
    yticks = 10
    xticks = 10
    dsshade.dynamic = False
    if height > width:
        fontscale = height / 400
    else:
        fontscale = width / 400
    if formatter_limits is not None:
        matplotlib.rcParams["axes.formatter.limits"] = formatter_limits
    hv.extension(backend)

    if backend == "matplotlib":
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = "Helvetica"
        title_font_size = 13 * fontscale
        axis_label_text_font_size = matplotlib.rcParams["font.size"] * fontscale
        major_label_text_font_size = matplotlib.rcParams["font.size"] * fontscale

    aspect = width / height
    # set linewidt (?)

    fig_inches = width / 60
    if width < height:
        fig_inches = height / 60

    y = y_axis
    if y == x_axis:
        df_group = df[[x_axis, by]]
    else:
        df_group = df[[x_axis, y, by]]

    types = df_group[by].unique()
    types = np.sort(types[::-1])

    min_x = df_group[x_axis].min()
    max_x = df_group[x_axis].max()
    min_y = df_group[y].min()
    max_y = df_group[y].max()

    del_x = (max_x - min_x) * 0.1
    del_y = max_y - min_y

    if min_y == max_y:
        min_y = max_y - 1
        max_y = max_y + 1
        del_y = 2

    cmap = {}
    dashmap = {}
    i = 2
    j = 2
    tmp_types = []
    for ty in types:
        if "Quant" not in ty and "Percent" not in ty:
            continue
        tmp_types.append(ty)
        if "Convective 400" in ty:
            cmap[ty] = plot_obj.cmap["Convective 400hPa 25. Quantile"]
        elif "Convective 600" in ty:
            cmap[ty] = plot_obj.cmap["Convective 600hPa 25. Quantile"]
        elif "Slantwise 400" in ty:
            cmap[ty] = plot_obj.cmap["Slantwise 400hPa 25. Quantile"]
        elif "Slantwise 600" in ty:
            cmap[ty] = plot_obj.cmap["Slantwise 600hPa 25. Quantile"]
        if backend == "matplotlib":
            if "25. Quantile" in ty or "25th Percentile" in ty:
                dashmap[ty] = "dashed"
            elif "75. Quantile" in ty or "75th Percentile" in ty:
                dashmap[ty] = "dashdot"
            else:
                dashmap[ty] = "solid"
        else:
            dashmap[ty] = [j, i]
            i += 2
            if i == 6:
                i = 2
                j += 2

    types = tmp_types
    if show_legend and backend == "matplotlib":
        overlay = hv.NdOverlay(
            {
                types[i]: hv.Curve((np.NaN, np.NaN)).opts(
                    opts.Curve(linestyle=dashmap[types[i]], color=cmap[types[i]])
                )
                for i in range(len(types))
            }
        )
    elif show_legend and backend == "bokeh":
        overlay = hv.NdOverlay(
            {
                types[i]: hv.Curve((np.NaN, np.NaN)).opts(
                    opts.Curve(line_dash=dashmap[types[i]], color=cmap[types[i]])
                )
                for i in range(len(types))
            }
        )

    if backend == "matplotlib":

        def format_axis(plot, element):
            ax = plot.handles["axis"]
            ax.set_title(
                title,
                loc="left",
                **{
                    "fontweight": "bold",
                    "fontsize": title_font_size,
                },
            )
            ax.set_ylabel(
                ax.get_ylabel(),
                **{
                    "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size,
                },
            )
            ax.set_xlabel(
                ax.get_xlabel(),
                **{
                    "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size,
                },
            )
            ax.minorticks_on()
            plot.handles["axis"] = ax

        # plot the area
        area_types = [
            "Slantwise 600hPa",
            "Slantwise 400hPa",
            "Convective 600hPa",
            "Convective 400hPa",
        ]
        cmap_area = {}
        for ty in np.unique(df[by]):
            if "Quantile" not in ty and "Percent" not in ty:
                cmap_area[ty] = plot_obj.cmap[ty]
        df_tmp = df_group.loc[df_group[by].isin(area_types)]
        min_df = df_tmp.groupby([x_axis, by])[y_axis].min()
        min_df.name = "min_y"
        max_df = df_tmp.groupby([x_axis, by])[y_axis].max()
        max_df.name = "max_y"
        df_area = pd.concat([min_df, max_df], axis=1)

        area_plot = None
        for ty in cmap_area:
            tmp = df_area.xs(ty, level="type").hvplot.area(
                x_axis,
                "min_y",
                "max_y",
                alpha=0.5,
                color=cmap_area[ty],
                width=width,
                height=height,
                dynamic=False,
                ylim=(
                    min_y - del_y * lim_multiplier,
                    max_y + del_y * lim_multiplier,
                ),
                xlim=(min_x, max_x),
                yticks=yticks,
                xticks=xticks,
                aspect=aspect,
                grid=True,
                legend=False,
            )
            if area_plot is None:
                area_plot = tmp
            else:
                area_plot = area_plot * tmp

        # plot each type
        line_plots = None
        for ty in types:
            if "Quant" not in ty and "Perc" not in ty:
                continue
            tmp = (
                df_group.loc[df_group["type"] == ty]
                .hvplot.line(
                    x=x_axis,
                    y=y,
                    color=cmap[ty],
                    ylabel=latexify.parse_word(y).upper()
                    + " "
                    + latexify.get_unit(y, True),
                    xlabel=latexify.parse_word(x_axis).upper(),
                    width=width,
                    height=height,
                    dynamic=False,
                    ylim=(
                        min_y - del_y * lim_multiplier,
                        max_y + del_y * lim_multiplier,
                    ),
                    xlim=(min_x, max_x),
                    yticks=yticks,
                    xticks=xticks,
                    aspect=aspect,
                    grid=True,
                    legend=False,
                )
                .opts(
                    initial_hooks=[format_axis],
                    linestyle=dashmap[ty],
                )
            )
            if line_plots is None:
                line_plots = tmp
            else:
                line_plots = line_plots * tmp
        # add all together
        if show_legend:
            img = area_plot * line_plots * overlay
        else:
            img = area_plot * line_plots

    if backend == "matplotlib":
        img = img.opts(
            show_grid=True,
            fontscale=fontscale,
            fig_inches=fig_inches,
            # title=title,
        )
    else:
        img = img.opts(
            show_grid=True,
            fontscale=fontscale,
            title=title,
        )

    renderer = hv.Store.renderers[backend].instance(fig="png", dpi=300)
    filetype = ".png"
    plot_path = store_path + "_" + y_axis
    i = 0
    save = plot_path + "_" + "{:03d}".format(i)
    while os.path.isfile(save + filetype):
        i = i + 1
        save = plot_path + "_" + "{:03d}".format(i)

    renderer.save(img, save)

    return img


def plot_scatter(
    df,
    y_axis,
    store_path,
    width=1959,
    height=1224,
    formatter_limits=None,
    datashade=True,
    s=None,
    title="",
    show_legend=False,
    backend="bokeh",
    verbosity=0,
):
    """
    Create a scatter plot and save it with "time_after_ascent"
    on the x-axis.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe loaded with load_data(..).
    y_axis : List of string
        The name of the column for the y-axis for each plot.
    store_path : string
        Path (and name) where to save images. y_axis will be appended.
    formatter_limits : tuple of ints
            Lower and upper limits for formatting x- and y-axis.
    width : int
        The width of the plot (bokeh) or the relative width of the plot
        (matplotlib).
    height : int
        The height of the plot (bokeh) or the relative height of the plot
        (matplotlib).
    datashade : boolean
        Use datashade for the dataset except for those with type "Quantile".
    s : int
        Size of the dots.
    title : string
        Title for this plot.
    show_legend : boolean
        Wether to plot the legend or not.
    backend : string
        Choose a backend for plotting. Options are:
        matplotlib: Most plots should be fine with it.
        bokeh: Recommended.
    verbosity : int
        Set verbosity level.
        0: No output except for exceptions
        1: Print loading statements (not in this function)
        2: Print printing statements
        3: Print additional statements

    Returns
    -------
    Holoviews image.
    """
    x_axis = "time_after_ascent_h"
    plot_obj = Deriv_dask(
        direc="",
        parquet=False,
        netcdf=True,
        columns=None,
        backend=backend,
        file_ending="",
    )
    by = "type"
    lim_multiplier = 0.1
    yticks = 10
    xticks = 10
    dsshade.dynamic = False
    if height > width:
        fontscale = height / 400
    else:
        fontscale = width / 400
    if formatter_limits is not None:
        matplotlib.rcParams["axes.formatter.limits"] = formatter_limits
    hv.extension(backend)

    if backend == "matplotlib":
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = "Helvetica"
        title_font_size = 13 * fontscale
        axis_label_text_font_size = matplotlib.rcParams["font.size"] * fontscale
        major_label_text_font_size = matplotlib.rcParams["font.size"] * fontscale

    aspect = width / height
    if s is None:
        if height > width:
            scatter_size = int(height / 200)
        else:
            scatter_size = int(width / 200)
    else:
        scatter_size = s

    fig_inches = width / 60
    if width < height:
        fig_inches = height / 60

    y = y_axis
    if y == x_axis:
        df_group = df[[x_axis, by]]
    else:
        df_group = df[[x_axis, y, by]]

    types = df_group[by].unique()
    types = np.sort(types[::-1])

    min_x = df_group[x_axis].min()
    max_x = df_group[x_axis].max()
    min_y = df_group[y].min()
    max_y = df_group[y].max()

    del_x = (max_x - min_x) * 0.1
    del_y = max_y - min_y

    # datashade cannot handle only zero valued plots
    if min_y == max_y and datashade:
        return
    if datashade and df_group.empty:
        return

    if min_y == max_y:
        min_y = max_y - 1
        max_y = max_y + 1
        del_y = 2

    cmap = {}
    for ty in types:
        cmap[ty] = plot_obj.cmap[ty]

    if show_legend and backend == "matplotlib":
        overlay = hv.NdOverlay(
            {
                types[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                    opts.Scatter(s=scatter_size * 8, color=cmap[types[i]])
                )
                for i in range(len(types))
            }
        )
    elif show_legend and backend == "bokeh":
        overlay = hv.NdOverlay(
            {
                types[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                    # types[i]: hv.Scatter((0, 0)).opts(
                    opts.Scatter(size=scatter_size * 8, color=cmap[types[i]])
                )
                for i in range(len(types))
            }
        )

    cmap_values = []
    other_types = []
    for ty in types:
        if "Quantile" in ty:
            other_types.append(ty)
            cmap_values.append(plot_obj.cmap[ty])

    if backend == "matplotlib":

        def format_axis(plot, element):
            ax = plot.handles["axis"]
            ax.set_title(
                title,
                loc="left",
                **{
                    "fontweight": "bold",
                    "fontsize": title_font_size,
                },
            )
            ax.set_ylabel(
                ax.get_ylabel(),
                **{
                    "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size,
                },
            )
            ax.set_xlabel(
                ax.get_xlabel(),
                **{
                    "fontstyle": "italic",
                    "fontsize": axis_label_text_font_size,
                },
            )
            ax.minorticks_on()
            plot.handles["axis"] = ax

        img_other = df_group.loc[df_group["type"].isin(other_types)].hvplot.scatter(
            x=x_axis,
            y=y,
            by=by,
            color=cmap_values,
            datashade=False,
            s=scatter_size,
            legend=False,
        )

        img = (
            df_group.hvplot.scatter(
                x=x_axis,
                y=y,
                by=by,
                cmap=cmap,
                label=None,
                ylabel=latexify.parse_word(y).upper()
                + " "
                + latexify.get_unit(y, True),
                xlabel=latexify.parse_word(x_axis).upper(),
                datashade=datashade,
                width=width,
                height=height,
                dynamic=False,
                ylim=(
                    min_y - del_y * lim_multiplier,
                    max_y + del_y * lim_multiplier,
                ),
                xlim=(min_x, max_x),
                dynspread=datashade,
                yticks=yticks,
                xticks=xticks,
                aspect=aspect,
                grid=True,
            ).opts(
                initial_hooks=[format_axis],
            )
            * img_other
        )
        img = img.opts(
            fontscale=fontscale,
            fig_inches=fig_inches,
        )
    else:
        img_other = df_group.loc[df_group["type"].isin(other_types)].hvplot.scatter(
            x=x_axis,
            y=y,
            by=by,
            color=cmap_values,
            datashade=False,
            size=scatter_size,
            legend=False,
            width=width,
            height=height,
        )

        img = (
            df_group.hvplot.scatter(
                x=x_axis,
                y=y,
                by=by,
                cmap=cmap,
                label=None,
                ylabel=latexify.parse_word(y).upper()
                + " "
                + latexify.get_unit(y, True),
                xlabel=latexify.parse_word(x_axis).upper(),
                datashade=datashade,
                width=width,
                height=height,
                dynamic=False,
                ylim=(
                    min_y - del_y * lim_multiplier,
                    max_y + del_y * lim_multiplier,
                ),
                xlim=(min_x, max_x),
                dynspread=datashade,
                yticks=yticks,
                xticks=xticks,
                grid=True,
                legend=False,
            )
            * img_other
        )
        img = img * overlay * img_other

    if show_legend:
        img = img * overlay

    if backend == "matplotlib":
        img = img.opts(
            show_grid=True,
            fontscale=fontscale,
            fig_inches=fig_inches,
        )
    else:
        img = img.opts(
            show_grid=True,
            fontscale=fontscale,
            title=title,
        )

    renderer = hv.Store.renderers[backend].instance(fig="png", dpi=300)
    filetype = ".png"
    plot_path = store_path + "_" + y_axis
    i = 0
    save = plot_path + "_" + "{:03d}".format(i)
    while os.path.isfile(save + filetype):
        i = i + 1
        save = plot_path + "_" + "{:03d}".format(i)

    renderer.save(img, save)

    return img


def plot_map(
    df,
    store_path,
    pollon,
    pollat,
    color_code=None,
    color_label=None,
    cmap="viridis",
    width=1959,
    height=1224,
    datashade=True,
    alpha=None,
    s=None,
    title="",
    tiles="ESRI",
    backend="matplotlib",
):
    """
    Plot trajectories over a map.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe loaded with load_data(..).
    store_path : string
        Path (and name) where to save images. y_axis will be appended.
    pollon : float
        Longitude for the rotated pole.
    pollat : float
        Latitude for the rotated pole.
    color_code : String
        Name of the column for color coding.
    color_label : String
        Label to plot for the colorbar.
    cmap : string
        Name of the colormap for the colorcode.
    width : int
        The width of the plot (bokeh) or the relative width of the plot
        (matplotlib).
    height : int
        The height of the plot (bokeh) or the relative height of the plot
        (matplotlib).
    datashade : boolean
        Use datashade to plot large amounts of data.
    alpha : float
        Alpha value for the dots.
    s : int
        Size of the scatter dots (if datashade is not used).
    title : string
        Title for this plot.
    tiles : string
        Define the tiles to use for the map, such as ESRI, CartoLight,
        EsriNatGeo. A full list of predefined sources is available via
        the dictionary holoviews.element.tiles.tile_sources. Using new
        sources starts a download of those.
    backend : string
        Choose a backend for plotting. Options are:
        matplotlib: Most plots should be fine with it.
        bokeh: Recommended.

    Returns
    -------
    Holoviews image
    """
    hv.extension(backend)
    rotated = crs.RotatedPole(
        pole_longitude=pollon,
        pole_latitude=pollat,
    )
    dsshade.dynamic = False

    fontscale = height / 350
    if backend == "matplotlib":
        fontscale = height / 190
        title_font_size = 13 * fontscale

    fig_inches = width / 60
    if width < height:
        fig_inches = height / 60
    if s is None:
        if height > width:
            scatter_size = int(height / 200)
        else:
            scatter_size = int(width / 200)
    else:
        scatter_size = s

    if datashade:
        img = df.hvplot.points(
            x="lon",
            y="lat",
            geo=True,
            dynamic=False,
            crs=rotated,
            datashade=datashade,
            coastline="10m",
            # tiles=tiles,
            hover=False,
            xlabel="Longitude",
            ylabel="Latitude",
            colorbar=False,
            xticks=10,
            yticks=10,
        )
    else:
        img = df.hvplot.points(
            x="lon",
            y="lat",
            geo=True,
            c=color_code,
            clabel=color_label,
            dynamic=False,
            crs=rotated,
            cmap=cmap,
            datashade=datashade,
            coastline="10m",
            # tiles=tiles,
            s=scatter_size,
            hover=False,
            alpha=alpha,
            xlabel="Longitude",
            ylabel="Latitude",
            xticks=10,
            yticks=10,
        )
    if backend == "matplotlib":

        def format_axis(plot, element):
            ax = plot.handles["axis"]
            ax.set_title(
                title,
                loc="left",
                **{
                    "fontweight": "bold",
                    "fontsize": title_font_size,
                },
            )
            plot.handles["axis"] = ax

        img = img.opts(
            initial_hooks=[format_axis],
            fontscale=fontscale,
            fig_inches=fig_inches,
        )
    else:
        img = img.opts(
            fontscale=fontscale,
            width=width,
            height=height,
            title=title,
        )

    filetype = ".png"
    if datashade:
        plot_path = store_path + "_map_datashade"
    else:
        plot_path = store_path + "_map"
    i = 0
    save = plot_path + "_" + "{:03d}".format(i)
    while os.path.isfile(save + filetype):
        i = i + 1
        save = plot_path + "_" + "{:03d}".format(i)
    renderer = hv.Store.renderers[backend].instance(fig="png", dpi=300)
    renderer.save(img, save)
    return img


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Create scatter plots with temperature, specific humidity,
            liquid hydrometeor content (qc+qr), cold hydrometeor content (qi+qs+qg+qh).
            Also plot 2D plot with pressure color coded on a map.
            The map uses ESRI as default. An overview how to cite this is given at
            https://support.esri.com/en/technical-article/000012040.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        required=True,
        type=str,
        help=textwrap.dedent(
            """\
            Path to folder with NetCDF-files from a COSMO simulation.
            If the path ends with '.h5', a HDF5-file will be loaded which should
            contain a pandas.Dataframe.
            """
        ),
    )
    parser.add_argument(
        "--traj_type",
        type=str,
        default="conv",
        help=textwrap.dedent(
            """\
            Choose which type of trajectories shall be plotted. Options are
            "conv" for all convective trajectories, "conv_400" for convective
            400 hPa trajectories and "conv_600" for 600 hPa trajectories.
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        type=str,
        default="../pics/cosmo_",
        help=textwrap.dedent(
            """\
            Path (and name) where to save images.
            """
        ),
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1200,
        help=textwrap.dedent(
            """\
            Width in pixels for the plot.
            """
        ),
    )
    parser.add_argument(
        "--height",
        type=int,
        default=800,
        help=textwrap.dedent(
            """\
            Height in pixels for the plot.
            """
        ),
    )
    parser.add_argument(
        "--store_data",
        default=None,
        help=textwrap.dedent(
            """\
            Store the data as HDF5 file so that it is faster when loading it again.
            """
        ),
    )
    parser.add_argument(
        "--pollon",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            If a HDF5 file is loaded, you can define here the longitude for the
            rotated pole.
            """
        ),
    )
    parser.add_argument(
        "--pollat",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            If a HDF5 file is loaded, you can define here the latitude for the
            rotated pole.
            """
        ),
    )
    parser.add_argument(
        "--backend",
        default="bokeh",
        help=textwrap.dedent(
            """\
            Choose a backend for plotting. Options are:
            matplotlib: Most plots should be fine with it.
            bokeh: Recommended.
            """
        ),
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help=textwrap.dedent(
            """\
            Set verbosity level.
            0: No output except for exceptions
            1: Print loading statements
            2: Print printing statements
            3: Print additional statements
            """
        ),
    )

    args = parser.parse_args()
    if args.data_path[-3:] == ".h5" and (args.pollon is None or args.pollat is None):
        print("You try to load a HDF5-file but did not define pollon and pollat")
        print("ABORTING")
        quit()
    if args.data_path[-3:] == ".h5":
        if args.verbosity > 0:
            print(f"Loading HDF5 file from {args.data_path}")
        df = pd.read_hdf(args.data_path)
    else:
        df, pollon, pollat = load_data(args.data_path, args.traj_type, args.verbosity)
        if args.store_data is not None:
            store = pd.HDFStore(args.store_data)
            store["df"] = df

    if args.pollon is not None:
        pollon = args.pollon
    if args.pollat is not None:
        pollat = args.pollat

    df["type"] = df["type"].replace(
        {
            "Convective 400hPa 25. Quantile": "Convective 400hPa 25th Percentile",
            "Convective 400hPa 50. Quantile": "Convective 400hPa 50th Percentile",
            "Convective 400hPa 75. Quantile": "Convective 400hPa 75th Percentile",
            "Convective 600hPa 25. Quantile": "Convective 600hPa 25th Percentile",
            "Convective 600hPa 50. Quantile": "Convective 600hPa 50th Percentile",
            "Convective 600hPa 75. Quantile": "Convective 600hPa 75th Percentile",
        },
    )

    if args.verbosity > 1:
        print("Plotting pressure")
    _ = plot_line(
        df=df,
        y_axis="pressure_hPa",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=30,
        title="(b) Pressure of All Trajectories",
        datashade=True,
        show_legend=True,
        backend="matplotlib",
    )
    if args.verbosity > 1:
        print("Plotting liquid water content")
    _ = plot_line(
        df=df,
        y_axis="Q_liquid",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=30,
        title="(b) Liquid Water Content of All Trajectories",
        datashade=True,
        show_legend=False,
        backend="matplotlib",
    )
    if args.verbosity > 1:
        print("Plotting frozen water content")
    _ = plot_line(
        df=df,
        y_axis="Q_cold",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=30,
        title="(c) Frozen Water Content of All Trajectories",
        datashade=True,
        show_legend=False,
        backend="matplotlib",
    )
    if args.verbosity > 1:
        print("Plotting temperature")
    _ = plot_line(
        df=df,
        y_axis="T",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=30,
        title="(d) Temperature of All Trajectories",
        datashade=True,
        show_legend=False,
        backend="matplotlib",
    )
    if args.verbosity > 1:
        print("Plotting specific humidity")
    _ = plot_line(
        df=df,
        y_axis="Specific_Humidity",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=30,
        title="(a) Specific Humidity of All Trajectories",
        datashade=True,
        show_legend=True,
        backend="matplotlib",
    )
    if args.verbosity > 1:
        print("Plotting all trajectories on a map")
    _ = plot_map(
        df=df,
        pollon=pollon,
        pollat=pollat,
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        title="(a) All Trajectories",
        tiles="EsriNatGeo",
        datashade=True,
    )

    if args.verbosity > 1:
        print("Plotting color coded pressure on a map with limited time")
    _ = plot_map(
        df=df.loc[
            (df["time_after_ascent"] >= -2800) & (df["time_after_ascent"] <= 26000)
        ],
        color_code="pressure_hPa",
        pollon=pollon,
        pollat=pollat,
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        title="Convective Trajectories and Pressure",
        tiles="EsriNatGeo",
        s=1,
        alpha=1,
        datashade=False,
        color_label="Pressure [hPa]",
    )
