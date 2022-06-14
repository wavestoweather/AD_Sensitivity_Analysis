try:
    import latexify
    from plot_cosmo import add_liquid_content, add_cold_content, add_sepcific_humidity
except:
    import scripts.latexify as latexify
    from scripts.plot_cosmo import (
        add_liquid_content,
        add_cold_content,
        add_sepcific_humidity,
    )

from holoviews.operation.datashader import datashade as dsshade
import holoviews as hv
from holoviews import opts
import hvplot.pandas
import matplotlib
import numpy as np
import os
import xarray as xr


def plot_scatter(
    df,
    y_axis,
    store_path,
    width=1959,
    height=1224,
    formatter_limits=None,
    s=None,
    title="",
    show_legend=False,
    backend="bokeh",
    verbosity=0,
):
    """
    Create a scatter plot and save it using bokeh with "time_after_ascent"
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
    verbosity : int
        Set verbosity level.
        0: No output except for exceptions
        1: Print loading statements (not in this function)
        2: Print printing statements
        3: Print additional statements

    Returns
    -------
    holoviews image
    """
    x_axis = "time_after_ascent_h"

    by = "Simulation"
    hv.extension(backend)
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

    cmap = matplotlib.pyplot.get_cmap("tab10")
    cmap_map = {}
    cmap_values = []
    legend_names = []
    for i, rt in enumerate(np.unique(df["Simulation"])):
        cmap_values.append(matplotlib.colors.to_hex(cmap(i)[0:-1]))
        cmap_map[rt] = matplotlib.colors.to_hex(cmap(i)[0:-1])
        legend_names.append(rt)

    if backend == "bokeh":
        img = df_group.hvplot.scatter(
            x=x_axis,
            y=y,
            by=by,
            colormap=cmap_values,
            ylabel=latexify.parse_word(y).upper() + " " + latexify.get_unit(y, True),
            xlabel=latexify.parse_word(x_axis).upper(),
            yticks=yticks,
            xticks=xticks,
            width=width,
            height=height,
            alpha=0.5,
            grid=True,
            size=scatter_size,
            legend="top_right",
        )
        img = img.opts(
            opts.Scatter(
                fontscale=fontscale,
                title=title,
            )
        )
    else:
        img = df_group.hvplot.scatter(
            x=x_axis,
            y=y,
            by=by,
            colormap=cmap_values,
            ylabel=latexify.parse_word(y).upper() + " " + latexify.get_unit(y, True),
            xlabel=latexify.parse_word(x_axis).upper(),
            yticks=yticks,
            xticks=xticks,
            width=width,
            height=height,
            alpha=0.5,
            grid=True,
            s=scatter_size,
            legend="top_right",
        )
        img = img.opts(
            opts.Scatter(
                aspect=aspect,
                s=scatter_size,
                show_grid=True,
                fontscale=fontscale,
                fig_inches=fig_inches,
                title=title,
            )
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


def load_data(
    data_cosmo_path,
    data_sim_path,
    verbosity=0,
    traj=0,
):
    file_list_cosmo = []
    for f in os.listdir(data_cosmo_path):
        file_list_cosmo.append(os.path.join(data_cosmo_path, f))
    file_list_cosmo = np.sort(np.asarray(file_list_cosmo))

    ds = xr.open_dataset(file_list_cosmo[0], decode_times=False, engine="netcdf4")
    cols = list(ds.keys())

    cols_final = []
    for col in cols:
        if (
            "_IN" not in col
            and "WCB_flag" not in col
            and "dp2h" not in col
            and "Q_TURBULENCE" not in col
            and "type" not in col
        ):
            cols_final.append(col)

    cols_final_sim = cols_final.copy()
    cols_final_sim.append("QH")
    cols_final_sim.append("QH_OUT")
    print(cols_final)
    file_list = []
    for f in os.listdir(data_sim_path):
        file_list.append(os.path.join(data_sim_path, f))
    file_list = np.sort(np.asarray(file_list))
    df_sim = None
    min_x = None
    max_x = None
    print(file_list)
    print("cosmo:")
    print(file_list_cosmo)
    for f in file_list:
        if verbosity > 0:
            print(f"Loading {f}")
        ds = xr.open_dataset(f, decode_times=False, engine="netcdf4")[cols_final_sim]

        ds = ds.loc[{"trajectory": ds["trajectory"][0]}]
        if df_sim is not None:
            df_sim = df_sim.append(ds.to_dataframe())
        else:
            min_x = np.min(ds["time_after_ascent"]).item()
            max_x = np.max(ds["time_after_ascent"]).item()
            if verbosity > 2:
                print(f"Got min time = {min_x} s, max time = {max_x} s.")
            df_sim = ds.to_dataframe()
    df_sim["Simulation"] = "Our Sim."

    df_cosmo = None
    for f in file_list_cosmo:
        if verbosity > 0:
            print(f"Loading {f}")
        ds = xr.open_dataset(f, decode_times=False, engine="netcdf4")[cols_final]
        ds = ds.loc[{"trajectory": ds["trajectory"][0]}]
        # ds = ds.where(ds["trajectory"] == ds["trajectory"][traj])
        ds = ds.where(
            ((ds["time_after_ascent"] >= min_x) & (ds["time_after_ascent"] <= max_x)),
            drop=True,
        )
        if df_cosmo is not None:
            df_cosmo = df_cosmo.append(ds.to_dataframe())
        else:
            df_cosmo = ds.to_dataframe()
    df_cosmo["Simulation"] = "COSMO"
    df_cosmo["QH"] = 0
    df_cosmo["QH_OUT"] = 0

    # drop all the NaN entries now that coordinates are collapsed to columns
    df_cosmo = df_cosmo.loc[
        (df_cosmo["time_after_ascent"] >= min_x)
        & (df_cosmo["time_after_ascent"] < max_x)
    ]
    df = df_cosmo.append(df_sim)

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

    df["Q_total"] = df["Q_cold"] + df["Q_liquid"] + df["Specific_Humidity"]
    if verbosity > 2:
        print("Finished loading data")
    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        """
        For comparison of cosmo simulation and our simulation.
        Create scatter plots with temperature, specific humidity,
        liquid hydrometeor content (qc+qr), cold hydrometeor content (qi+qs+qg+qh).
        """,
    )
    parser.add_argument(
        "--data_cosmo_path",
        required=True,
        type=str,
        help="""
        Path to folder with NetCDF-files from a COSMO simulation.
        """,
    )
    parser.add_argument(
        "--data_sim_path",
        required=True,
        type=str,
        help="""
        Path to folder with NetCDF-files from our simulation.
        """,
    )
    parser.add_argument(
        "--store_path",
        type=str,
        default="../pics/cosmo_comparison_",
        help="""
        Path (and name) where to save images.
        """,
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1200,
        help="""
        Width in pixels for the plot.
        """,
    )
    parser.add_argument(
        "--height",
        type=int,
        default=800,
        help="""
        Height in pixels for the plot.
        """,
    )
    parser.add_argument(
        "--traj",
        type=int,
        default=0,
        help="""
        Define the trajectory index of the files in data_cosmo_path.
        """,
    )
    parser.add_argument(
        "--backend",
        default="matplotlib",
        help="""
        Choose a backend for plotting. Options are:
        matplotlib: Most plots should be fine with it.
        bokeh: Recommended.
        """,
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help="""
        Set verbosity level.
        0: No output except for exceptions
        1: Print loading statements
        2: Print printing statements
        3: Print additional statements
        """,
    )

    args = parser.parse_args()

    df = load_data(
        data_cosmo_path=args.data_cosmo_path,
        data_sim_path=args.data_sim_path,
        verbosity=args.verbosity,
        traj=args.traj,
    )

    if args.verbosity > 1:
        print("Plotting pressure")
    _ = plot_scatter(
        df=df,
        y_axis="pressure_hPa",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Pressure Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting liquid water content")
    _ = plot_scatter(
        df=df,
        y_axis="Q_liquid",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="(b) Liquid Water Content Comparison",
        show_legend=False,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting cold hydrometeor content")
    _ = plot_scatter(
        df=df,
        y_axis="Q_cold",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="(c) Frozen Water Content Comparison",
        show_legend=False,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting temperature")
    _ = plot_scatter(
        df=df,
        y_axis="T",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Temperature Comparison",
        show_legend=False,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting specific humidity")
    _ = plot_scatter(
        df=df,
        y_axis="Specific_Humidity",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="(a) Specific Humidity Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting graupel mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QG",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Graupel Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting snow mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QS",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="(d) Snow Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting ice mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QI",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Ice Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )

    if args.verbosity > 1:
        print("Plotting hail mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QH",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Hail Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting vapor mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QV",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Water Vapor Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting cloud droplet mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QC",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Cloud Droplet Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting rain droplet mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QR",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Rain Droplet Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting total amount mass density")
    _ = plot_scatter(
        df=df,
        y_axis="Q_total",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Total Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    # Cosmo uses positive values, we use negative values for sedimentation
    df["QR_OUT"] = np.abs(df["QR_OUT"])
    df["QI_OUT"] = np.abs(df["QI_OUT"])
    df["QS_OUT"] = np.abs(df["QS_OUT"])
    df["QG_OUT"] = np.abs(df["QG_OUT"])
    if args.verbosity > 1:
        print("Plotting rain out mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QR_OUT",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Rain Sed. Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting ice out mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QI_OUT",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Ice Sed. Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting snow out mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QS_OUT",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Snow Sed. Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting graupel out mass density")
    _ = plot_scatter(
        df=df,
        y_axis="QG_OUT",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Graupel Sed. Mass Density Comparison",
        show_legend=True,
        backend=args.backend,
    )

    if args.verbosity > 1:
        print("Plotting cloud droplets")
    _ = plot_scatter(
        df=df,
        y_axis="NCCLOUD",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Cloud Number Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting rain number")
    _ = plot_scatter(
        df=df,
        y_axis="NCRAIN",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Rain Number Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting ice number")
    _ = plot_scatter(
        df=df,
        y_axis="NCICE",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Ice Number Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting snow number")
    _ = plot_scatter(
        df=df,
        y_axis="NCSNOW",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Snow Number Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
    if args.verbosity > 1:
        print("Plotting graupel number")
    _ = plot_scatter(
        df=df,
        y_axis="NCGRAUPEL",
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        s=14,
        title="Graupel Number Density Comparison",
        show_legend=True,
        backend=args.backend,
    )
