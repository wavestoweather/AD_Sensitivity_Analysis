from argparse import ArgumentParser
from iris.analysis.cartography import rotate_pole
from mpl_toolkits.basemap import Basemap
from PIL import Image

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import numpy as np
import pandas as pd
from scipy import stats
import sys
import os


def plot_ratio_deriv_line(df_dict, out_params=None):
    """
    Plot derivative ratios of one trajectory. The x-axis is the timesteps,
    the y-axis are the derivative ratios. Add a bar for flagged timesteps.

    Parameters
    ----------
    df_dict : dic of pandas.Datafram
        A dictionary of pandas.Dataframe with key the output parameter.
        Dataframes have columns trajectory, timestep, out_param, in_param,
        deriv, ratio_deriv and optionally MAP.
    """
    if out_params is None:
        out_params = df_dict.keys()
    for out_param in out_params:
        df = df_dict[out_param]
        if df.empty:
            continue

        min_time = df.timestep.unique().min()
        max_time = df.timestep.unique().max()
        dt = (max_time - min_time + 19) / 20
        x_ticks = np.arange(min_time, max_time + 19, dt)

        _, ax = plt.subplots()

        ax = sns.lineplot(x="timestep", y="ratio_deriv",
                          data=df, hue="in_param", ax=ax)
        ax.set_title("Deriv. ratio of {}".format(out_param))
        ax.set_xticks(x_ticks)

        # Plot the area that had been flagged
        df_mapped = df[df.MAP == True]
        if not df_mapped.empty:
            min_y = df["ratio_deriv"].min()
            max_y = df["ratio_deriv"].max()
            ax.fill_between(df_mapped["timestep"], 0, 1,
                facecolor="khaki", alpha=0.5)

        i = 0
        save = ("pics/line_" + out_param
                + "_" + "{:03d}".format(i) + ".png")
        while os.path.isfile(save):
            i = i+1
            save = ("pics/_line_" + out_param
                    + "_" + "{:03d}".format(i) + ".png")

        print("Saving to " + save)
        plt.show()
        plt.savefig(save, dpi=300)
        plt.close()


def plot_res_line(df, out_param, dots=False):
    """
    Plot results for an out_param of one trajectory. The x-axis is the timesteps,
    the y-axis is the parameter. Add a bar for flagged timesteps.

    Parameters
    ----------
    df : pandas.Datafram
        Dataframe with columns trajectory, timestep, MAP and the out_param.
    out_param : String
        The out_param to plot
    dots : Bool
        Plot dots every 20 seconds
    """
    min_time = df.timestep.unique().min()
    max_time = df.timestep.unique().max()
    dt = (max_time - min_time + 19) / 20
    x_ticks = np.arange(min_time, max_time + 19, dt)

    _, ax = plt.subplots()

    ax = sns.lineplot(x="timestep", y=out_param,
                        data=df, hue="MAP", ax=ax)
    ax.set_title("Simulation of {}".format(out_param))
    ax.set_xticks(x_ticks)

    # Plot the area that had been flagged
    df_mapped = df[df.MAP == True]
    if not df_mapped.empty:
        min_y = df["ratio_deriv"].min()
        max_y = df["ratio_deriv"].max()
        ax.fill_between(df_mapped["timestep"], 0, 1,
            facecolor="khaki", alpha=0.5)

    # Plot dots every 20 seconds
    if dots:
        x_dots = []
        y_dots = []
        for time in df["timestep"]:
            if time%20 == 0:
                x_dots.append(time)
                y_dots.append( df[df.timestep == time][out_param].values[0] )
        plt.scatter(x_dots, y_dots, marker='x', c="black")

    i = 0
    save = ("pics/line_" + out_param
            + "_" + "{:03d}".format(i) + ".png")
    while os.path.isfile(save):
        i = i+1
        save = ("pics/_line_" + out_param
                + "_" + "{:03d}".format(i) + ".png")

    print("Saving to " + save)
    plt.show()
    plt.savefig(save, dpi=300)
    plt.close()

if __name__ == "__main__":
    parser = ArgumentParser(description=
            '''Load a mapped wcb netCDF file or the derivatives from a
            simulation and plot it. In case of a netCDF file, the flagged
            regions are colored differently. In case of the derivative
            file, one can plot either only flagged or not-flagged areas
            or both. Each area (flagged and not flagged) is plotted with
            normalized derivative ratios to make it easier to compare
            if needed.''')
    parser.add_argument("-i", "--input", type=str, default=None,
            help='''Path to the data folder (!) which shall be read.
            The datafiles should look like
            "wcb1000.0_traj1_start_over_MAP_20160922_02_diff_0.txt", where
            "traj" is important.''')
    parser.add_argument("-n", "--netcdf", type=str, default=None,
            help='''Path to a wcb netCDF file with flags "MAP".''')
    parser.add_argument("-t", "--trajectory", nargs='+', type=int, default=1,
            help='''A list of trajectories to consider here. Using a selected
            number of trajectories helps reducing the amount of used RAM
            and makes it easier to see anything in the image.''')
    parser.add_argument("-r", "--rotate", action="store_true", default=False,
            help='''Rotate longitude and latitude coordinates. For
            the derivatives, we assume pol_lat=90.0 and pollon=0.0 if
            no netCDF file is given.''')
    parser.add_argument("--input_param", nargs='+', type=str, default="dzeta",
            help='''A list of input parameters of the model to plot such
            as "dzeta" or "dbeta_r".''')
    parser.add_argument("--output_param", nargs='+', type=str, default="T",
            help='''A list of output parameters of the model to plot such
            as "T", "qr" or "p".''')
    parser.add_argument("--max_time", type=float, default=None,
            help='''Maximum time in seconds to plot.''')
    parser.add_argument("-e", "--epsilon", type=float, default=0.0,
            help='''Value to filter values with. Default is 0.0.''')
    parser.add_argument("-f", "--filter", action="store_true", default=False,
            help='''Filter values smaller than epsilon if used.''')
    parser.add_argument("--norm", action="store_true", default=True,
            help='''Normalize the derivatives for an area (consecutive
            timings with or without flag).''')

    args = parser.parse_args()

    try:
        from loader import load_mult_derivates_directory, load_nc
        from loader import rotate_df, norm_deriv, ratio_deriv
    except:
        from scripts.loader import load_mult_derivates_directory
        from scripts.loader import load_nc, rotate_df, norm_deriv, ratio_deriv
    from pylab import rcParams


    rcParams['figure.figsize'] = (16,10)

    df = None
    net_df = None
    pollon = 0.0
    pollat = 90.0

    if args.netcdf is not None:
        if args.rotate:
            net_df, pollat, pollon = load_nc(inp=args.netcdf,
                                             get_pol=True)
            rotate_df(net_df, pollon, pollat, "lon", "lat")
        else:
            net_df = load_nc(inp=args.netcdf,
                             get_pol=False)
        net_df = net_df[net_df.trajectory in args.trajectory]

        # Remove any nonsense
        net_df = net_df[net_df.time >= 0]
        if args.max_time is not None:
            net_df = net_df[net_df.time <= args.max_time]

    if args.input is not None:
        df = load_mult_derivates_directory(direc=args.input,
                                            filt=args.filter,
                                            EPSILON=args.epsilon,
                                            trajectories=args.trajectory,
                                            suffix=None)

        if args.max_time is not None:
            df = df[df.timestep <= args.max_time]
        if args.rotate:
            rotate_df(df, pollon, pollat, "LONGITUDE", "LATITUDE")

        ratio_dfs = {}
        for out_par in args.output_param:
            ratio_dfs[out_par] = ratio_deriv(df, out_par)
        if args.norm:
            norm_deriv(df)

    # Plot data
    if net_df is not None:
        for in_par in args.input_param:
            for out_par in args.output_param:
                plot_netcdf(df=net_df,
                            in_par=in_par,
                            out_par=out_par)

    if df is not None:
        for in_par in args.input_param:
            for out_par in args.output_param:
                plot_derif(df=df,
                           in_par=in_par,
                           out_par=out_par)



