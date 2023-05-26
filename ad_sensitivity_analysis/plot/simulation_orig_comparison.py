"""Plot our simulation vs original data.

"""
from holoviews.operation.datashader import datashade as dsshade
import holoviews as hv
from holoviews import opts
import matplotlib
import numpy as np

from ad_sensitivity_analysis.data_handler.loader import load_simulation_and_orig
from ad_sensitivity_analysis.plot import latexify
from ad_sensitivity_analysis.plot.aux_functions import save_plot_renderer


# pylint: disable=too-many-locals
def plot_scatter(
    df,
    y_axis,
    store_path,
    width=1959,
    height=1224,
    formatter_limits=None,
    dot_size=None,
    title="",
    backend="bokeh",
):
    """
    Create a scatter plot and save it using bokeh with "time_after_ascent"
    on the x-axis.

    Parameters
    ----------
    df : pandas.DataFrame
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
    dot_size : int
        Passed to hvplot.scatter for the dot size.
    title : string
        Title of the plot
    backend : string
        "bokeh" or "matplotlib" as backend to plot the image.

    Returns
    -------
    holoviews image
    """
    x_axis = "time_after_ascent_h"

    split_by = "Simulation"
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
    if dot_size is None:
        if height > width:
            dot_size = int(height / 200)
        else:
            dot_size = int(width / 200)

    fig_inches = width / 60
    if width < height:
        fig_inches = height / 60

    y = y_axis
    if y == x_axis:
        df_group = df[[x_axis, split_by]]
    else:
        df_group = df[[x_axis, y, split_by]]

    cmap = matplotlib.pyplot.get_cmap("tab10")
    cmap_map = {}
    cmap_values = []
    legend_names = []
    for i, sim in enumerate(np.unique(df["Simulation"])):
        cmap_values.append(matplotlib.colors.to_hex(cmap(i)[0:-1]))
        cmap_map[sim] = matplotlib.colors.to_hex(cmap(i)[0:-1])
        legend_names.append(sim)

    if backend == "bokeh":
        img = df_group.hvplot.scatter(
            x=x_axis,
            y=y,
            by=split_by,
            colormap=cmap_values,
            ylabel=latexify.parse_word(y).upper() + " " + latexify.get_unit(y, True),
            xlabel=latexify.parse_word(x_axis).upper(),
            yticks=yticks,
            xticks=xticks,
            width=width,
            height=height,
            alpha=0.5,
            grid=True,
            size=dot_size,
            legend="top_right",
        )
        # pylint: disable=no-member
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
            by=split_by,
            colormap=cmap_values,
            ylabel=latexify.parse_word(y).upper() + " " + latexify.get_unit(y, True),
            xlabel=latexify.parse_word(x_axis).upper(),
            yticks=yticks,
            xticks=xticks,
            width=width,
            height=height,
            alpha=0.5,
            grid=True,
            s=dot_size,
            legend="top_right",
        )
        # pylint: disable=no-member
        img = img.opts(
            opts.Scatter(
                aspect=aspect,
                s=dot_size,
                show_grid=True,
                fontscale=fontscale,
                fig_inches=fig_inches,
                title=title,
            )
        )
    save_plot_renderer(
        img,
        store_path + "_" + y_axis,
        hv.Store.renderers[backend].instance(fig="png", dpi=300),
    )
    return img


def main(args):
    """
    Plot some examples.

    Parameters
    ----------
    args :
        Arguments parsed via argparse.

    """
    df = load_simulation_and_orig(
        data_cosmo_path=args.data_cosmo_path,
        data_sim_path=args.data_sim_path,
        verbosity=args.verbosity,
        traj=args.traj,
    )
    yaxis_title_pairs = [
        ["pressure_hPa", "Pressure Comparison"],
        ["Q_liquid", "(b) Liquid Water Content Comparison"],
        ["Q_cold", "(c) Frozen Water Content Comparison"],
        ["T", "Temperature Comparison"],
        ["Specific_Humidity", "(a) Specific Humidity Comparison"],
        ["QG", "Graupel Mass Density Comparison"],
        ["QS", "(d) Snow Mass Density Comparison"],
        ["QI", "Ice Mass Density Comparison"],
        ["QH", "Hail Mass Density Comparison"],
        ["QV", "Water Vapor Mass Density Comparison"],
        ["QC", "Cloud Droplet Mass Density Comparison"],
        ["QR", "Rain Droplet Mass Density Comparison"],
        ["Q_total", "Total Mass Density Comparison"],
        ["QR_OUT", "Rain Sed. Mass Density Comparison"],
        ["QI_OUT", "Ice Sed. Mass Density Comparison"],
        ["QS_OUT", "Snow Sed. Mass Density Comparison"],
        ["QG_OUT", "Graupel Sed. Mass Density Comparison"],
        ["NCCLOUD", "Cloud Number Density Comparison"],
        ["NCRAIN", "Rain Number Density Comparison"],
        ["NCICE", "Ice Number Density Comparison"],
        ["NCSNOW", "Snow Number Density Comparison"],
        ["NCGRAUPEL", "Graupel Number Density Comparison"],
    ]
    # Cosmo uses positive values, we use negative values for sedimentation
    df["QR_OUT"] = np.abs(df["QR_OUT"])
    df["QI_OUT"] = np.abs(df["QI_OUT"])
    df["QS_OUT"] = np.abs(df["QS_OUT"])
    df["QG_OUT"] = np.abs(df["QG_OUT"])
    for y_axis, title in yaxis_title_pairs:
        if args.verbosity > 1:
            print(f"Plotting {title}")
        _ = plot_scatter(
            df=df,
            y_axis=y_axis,
            store_path=args.store_path,
            width=args.width,
            height=args.height,
            dot_size=14,
            title=title,
            backend=args.backend,
        )


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            For comparison of cosmo simulation and our simulation.
            Create scatter plots with temperature, specific humidity,
            liquid hydrometeor content (qc+qr), cold hydrometeor content (qi+qs+qg+qh).
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--data_cosmo_path",
        required=True,
        type=str,
        help=textwrap.dedent(
            """\
            Path to folder with NetCDF-files from a COSMO simulation.
            """
        ),
    )
    parser.add_argument(
        "--data_sim_path",
        required=True,
        type=str,
        help=textwrap.dedent(
            """\
            Path to folder with NetCDF-files from our simulation.
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        type=str,
        default="../pics/cosmo_comparison_",
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
        "--traj",
        type=int,
        default=0,
        help=textwrap.dedent(
            """\
            Define the trajectory index of the files in data_cosmo_path.
            """
        ),
    )
    parser.add_argument(
        "--backend",
        default="matplotlib",
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
    main(parser.parse_args())
