from argparse import ArgumentParser
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

try:
    import latexify
except:
    import scripts.latexify as latexify


def plot_ratio_deriv_line(df_dict, out_params=None, mapped=True,
                          in_params=None, **kwargs):
    """
    Plot derivative ratios of one trajectory. The x-axis is the timesteps,
    the y-axis are the derivative ratios. Add a bar for flagged timesteps.

    Parameters
    ----------
    df_dict : dic of pandas.Datafram
        A dictionary of pandas.Dataframe with key the output parameter.
        Dataframes have columns trajectory, timestep, out_param, in_param,
        deriv, ratio_deriv and optionally MAP.
    out_params : list of string
        List of keys to plot the derivatives for.
    mapped : boolean
        If true: plot the region, where "MAP" is true, ie where the wcb
        criterion is fullfilled.
    in_params : list of string
        Plot only the derivatives with respect to those in this list.
    kwargs : dict
        Keyword arguments are passed down to matplotlib.axes.Axes.plot() for
        the derivative plots.
    """
    if out_params is None:
        out_params = df_dict.keys()
    for out_param in out_params:
        df = df_dict[out_param]
        if df is None:
            continue
        if df.empty:
            continue

        def plot_helper(df):
            min_time = df.timestep.unique().min()
            max_time = df.timestep.unique().max()
            dt = (max_time - min_time + 19) / 20
            x_ticks = np.arange(min_time, max_time + 19, dt)

            _, ax = plt.subplots()

            ax = sns.lineplot(x="timestep", y="ratio_deriv",
                            data=df, hue="in_param", ax=ax, **kwargs)
            ax.set_title("Deriv. Ratio of {}".format(latexify.parse_word(out_param)))
            ax.set_xticks(x_ticks)

            # Change labels to latex versions
            legend = ax.get_legend()
            _, labels = ax.get_legend_handles_labels()
            for t, old in zip(legend.texts, labels):
                t.set_text(latexify.parse_word(old))
            ax.set_ylabel("Derivative ratio")

            # Plot the area that had been flagged
            if mapped:
                df_mapped = df[df.MAP == True]
                if not df_mapped.empty:
                    min_y = df["ratio_deriv"].min()
                    max_y = df["ratio_deriv"].max()
                    ax.fill_between(df_mapped["timestep"], min_y, max_y,
                        facecolor="khaki", alpha=0.3)

            i = 0
            save = ("pics/line_" + out_param
                    + "_" + "{:03d}".format(i) + ".png")
            while os.path.isfile(save):
                i = i+1
                save = ("pics/_line_" + out_param
                        + "_" + "{:03d}".format(i) + ".png")

            print("Saving to " + save)
            plt.savefig(save, dpi=300)
            plt.show()
            plt.close()

        if in_params is None:
            plot_helper(df)
        else:
            df_this = df.loc[df["in_param"].isin(in_params)]
            plot_helper(df_this)


def plot_res_line(df, out_param, dots=False, mapped=True, **kwargs):
    """
    Plot results for an out_param of one or more trajectories. The x-axis is the timesteps,
    the y-axis is the parameter. Add a bar for flagged timesteps.

    Parameters
    ----------
    df : pandas.Datafram
        Dataframe with columns trajectory, timestep, MAP and the out_param.
    out_param : String
        The out_param to plot
    dots : Bool
        Plot dots every 20 seconds
    mapped : boolean
        If true: plot the region, where "MAP" is true, ie where the wcb
        criterion is fullfilled.
    kwargs : dict
        Keyword arguments are passed down to matplotlib.axes.Axes.plot().
    """
    min_time = df.timestep.unique().min()
    max_time = df.timestep.unique().max()
    dt = (max_time - min_time + 19) / 20
    x_ticks = np.arange(min_time, max_time + 19, dt)

    _, ax = plt.subplots()

    if len(df.trajectory.unique()) > 1:
        ax = sns.lineplot(x="timestep", y=out_param,
                          data=df, hue="trajectory", ax=ax,
                          palette=sns.color_palette("husl", len(df.trajectory.unique())),
                          **kwargs)
    else:
        ax = sns.lineplot(x="timestep", y=out_param,
                          data=df, ax=ax, **kwargs)
        ax.legend(df.trajectory.unique(), title="trajectory")
    # Plot dots every 20 seconds
    if dots:
        x_dots = []
        y_dots = []
        for time in df["timestep"]:
            if time%20 == 0:
                y_values = df[df.timestep == time][out_param].values
                for y in y_values:
                    x_dots.append(time)
                    y_dots.append(y)
        plt.scatter(x_dots, y_dots, marker='x', c="black")


    ax.set_title("Simulation of {}".format(latexify.parse_word(out_param)))
    ax.set_xticks(x_ticks)

    # Plot the area that had been flagged
    if mapped:
        df_mapped = df[df.MAP == True]
        if not df_mapped.empty:
            min_y = df[out_param].min()
            max_y = df[out_param].max()
            ax.fill_between(df_mapped["timestep"], min_y, max_y,
                facecolor="khaki", alpha=0.3)

    i = 0
    save = ("pics/line_" + out_param
            + "_" + "{:03d}".format(i) + ".png")
    while os.path.isfile(save):
        i = i+1
        save = ("pics/line_" + out_param
                + "_" + "{:03d}".format(i) + ".png")

    print("Saving to " + save)
    plt.savefig(save, dpi=300)
    plt.show()
    plt.close()


def plot_same_orders(df_dict, out_params=None, mapped=True, in_params=None, **kwargs):
    """
    For each out_param, plot multiple figures with derivatives, where
    each figure shows derivatives of the same order.
    Parameters
    ----------
    df_dict : dic of pandas.Datafram
        A dictionary of pandas.Dataframe with key the output parameter.
        Dataframes have columns trajectory, timestep, out_param, in_param,
        deriv, ratio_deriv and optionally MAP.
    out_params : list of string
        List of keys to plot the derivatives for.
    mapped : boolean
        If true: plot the region, where "MAP" is true, ie where the wcb
        criterion is fullfilled.
    in_params : list of string
        Plot only the derivatives with respect to those in this list.
    kwargs : dict
        Keyword arguments are passed down to matplotlib.axes.Axes.plot() for
        the derivative plots.
    """
    if out_params is None:
        out_params = df_dict.keys()
    for out_param in out_params:
        df = df_dict[out_param]
        if df is None:
            continue
        if df.empty:
            continue

        if in_params is None:
            in_params_1 = df["in_param"].unique()

        sorted_tuples = []
        for in_p in in_params_1:
            this_df = df.loc[df["in_param"] == in_p]
            value = abs(this_df["ratio_deriv"].min())
            if abs(this_df["ratio_deriv"].max()) > value:
                value = abs(this_df["ratio_deriv"].max())
            sorted_tuples.append((in_p, value))
        sorted_tuples.sort(key=lambda tup: tup[1])
        while len(sorted_tuples) > 0:
            p, v = sorted_tuples.pop()
            in_params_2 = [p]
            while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                   and abs(v/sorted_tuples[-1][1]) < 10):
                p, v = sorted_tuples.pop()
                in_params_2.append(p)
            plot_ratio_deriv_line(
                df_dict,
                out_params=[out_param],
                in_params=in_params_2,
                mapped=mapped,
                **kwargs)


def plot_ratio_deriv_line_x(df_dict, out_params=None, mapped=True,
                          in_params=None, x_axis="timestep", **kwargs):
    """
    Plot derivative ratios of one trajectory. The x-axis is the timesteps
    by default,
    the y-axis are the derivative ratios. Add a bar for flagged timesteps.

    Parameters
    ----------
    df_dict : dic of pandas.Datafram
        A dictionary of pandas.Dataframe with key the output parameter.
        Dataframes have columns trajectory, timestep, out_param, in_param,
        deriv, ratio_deriv, out_param_value and optionally MAP.
    out_params : list of string
        List of keys to plot the derivatives for.
    mapped : boolean
        If true: plot the region, where "MAP" is true, ie where the wcb
        criterion is fullfilled.
    in_params : list of string
        Plot only the derivatives with respect to those in this list.
    x_axis : string
        The column to use as x-axis. Can be either "timestep" or
        "out_param" or an output parameter.
    kwargs : dict
        Keyword arguments are passed down to matplotlib.axes.Axes.plot() for
        the derivative plots.
    """
    if out_params is None:
        out_params = df_dict.keys()
    for out_param in out_params:
        df = df_dict[out_param]
        if df is None:
            continue
        if df.empty:
            continue

        def plot_helper(df):
            if x_axis not in df.keys():
                # Should be an output parameter
                df_tmp = df_dict[x_axis]
                df[x_axis] = df_dict[x_axis][out_param_value]

            min_time = df[x_axis].unique().min()
            max_time = df[x_axis].unique().max()
            dt = (max_time - min_time + 19) / 20
            x_ticks = np.arange(min_time, max_time + 19, dt)

            _, ax = plt.subplots()

            ax = sns.lineplot(x=x_axis, y="ratio_deriv",
                            data=df, hue="in_param", ax=ax, **kwargs)
            ax.set_title("Deriv. Ratio of {}".format(latexify.parse_word(out_param)))
            ax.set_xticks(x_ticks)

            # Change labels to latex versions
            legend = ax.get_legend()
            _, labels = ax.get_legend_handles_labels()
            for t, old in zip(legend.texts, labels):
                t.set_text(latexify.parse_word(old))
            ax.set_ylabel("Derivative ratio")

            # Plot the area that had been flagged
            if mapped:
                df_mapped = df[df.MAP == True]
                if not df_mapped.empty:
                    min_y = df["ratio_deriv"].min()
                    max_y = df["ratio_deriv"].max()
                    ax.fill_between(df_mapped[x_axis], min_y, max_y,
                        facecolor="khaki", alpha=0.3)

            i = 0
            save = ("pics/line_" + out_param
                    + "_" + "{:03d}".format(i) + ".png")
            while os.path.isfile(save):
                i = i+1
                save = ("pics/_line_" + out_param
                        + "_" + "{:03d}".format(i) + ".png")

            print("Saving to " + save)
            plt.savefig(save, dpi=300)
            plt.show()
            plt.close()

        if in_params is None:
            plot_helper(df)
        else:
            df_this = df.loc[df["in_param"].isin(in_params)]
            plot_helper(df_this)


def plot_same_orders_x(df_dict, out_params=None, mapped=True, in_params=None, x_axis="timestep", **kwargs):
    """
    For each out_param, plot multiple figures with derivatives, where
    each figure shows derivatives of the same order.
    Parameters
    ----------
    df_dict : dic of pandas.Datafram
        A dictionary of pandas.Dataframe with key the output parameter.
        Dataframes have columns trajectory, timestep, out_param, in_param,
        deriv, ratio_deriv and optionally MAP.
    out_params : list of string
        List of keys to plot the derivatives for.
    mapped : boolean
        If true: plot the region, where "MAP" is true, ie where the wcb
        criterion is fullfilled.
    in_params : list of string
        Plot only the derivatives with respect to those in this list.
    x_axis : string
        The column to use as x-axis. Can be either "timestep" or
        "out_param" or an output parameter.
    kwargs : dict
        Keyword arguments are passed down to matplotlib.axes.Axes.plot() for
        the derivative plots.
    """
    if out_params is None:
        out_params = df_dict.keys()
    for out_param in out_params:
        df = df_dict[out_param]
        if df is None:
            continue
        if df.empty:
            continue

        if in_params is None:
            in_params_1 = df["in_param"].unique()

        sorted_tuples = []
        for in_p in in_params_1:
            this_df = df.loc[df["in_param"] == in_p]
            value = abs(this_df["ratio_deriv"].min())
            if abs(this_df["ratio_deriv"].max()) > value:
                value = abs(this_df["ratio_deriv"].max())
            sorted_tuples.append((in_p, value))
        sorted_tuples.sort(key=lambda tup: tup[1])
        while len(sorted_tuples) > 0:
            p, v = sorted_tuples.pop()
            in_params_2 = [p]
            while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                   and abs(v/sorted_tuples[-1][1]) < 10):
                p, v = sorted_tuples.pop()
                in_params_2.append(p)
            plot_ratio_deriv_line(
                df_dict,
                out_params=[out_param],
                in_params=in_params_2,
                mapped=mapped,
                x_axis=x_axis,
                **kwargs)


def plot_orders_x(traj, x_axis_list):
    direc_path = "data/sim_fixed/tr{}".format(traj)
    max_time = None
    filt = True
    EPSILON = 0.0
    trajectory = [traj]
    df_mapped = load_mult_derivates_directory(direc=direc_path,
                                        filt=True,
                                        EPSILON=EPSILON,
                                        trajectories=trajectory,
                                        suffix=None)
    df_mapped = df_mapped[df_mapped.MAP == True]
    out_params = df_mapped.out_param.unique()
    lon = df_mapped[df_mapped.in_param == "LONGITUDE"]["deriv"].unique()
    lat = df_mapped[df_mapped.in_param == "LATITUDE"]["deriv"].unique()
    # Remove entries with coordinates. We just stored them in lat and lon
    df_mapped = df_mapped[df_mapped.in_param != "LONGITUDE"]
    df_mapped = df_mapped[df_mapped.in_param != "LATITUDE"]
    # Load output parameter values and store them accordingly
    df_sim_mapped = loader.load_output("data/sim_fixed/tr{}/wcb187920.0_traj{}_start_over_MAP_Flux_2_t000000_p001.txt".format(traj, traj),
                        sep=",",
                        refs="data/sim_fixed/tr{}/wcb187920.0_traj{}_start_over_MAP_Flux_2_t000000_p001_reference_values.txt".format(traj, traj))
    df_sim_mapped = df_sim_mapped[df_sim_mapped.MAP == True]
    def add_column(df, out_par, timestep):
        return df.loc[df["timestep"] == timestep][out_par]
    df_mapped["out_param_value"] = df_mapped.apply(lambda x: add_column(df_sim_mapped,
                                                                        df_mapped["out_param"],
                                                                        df_mapped["timestep"]), axis=1)

    ratio_dfs_mapped = {}
    for out_par in out_params:
        ratio_dfs_mapped[out_par] = ratio_deriv(df_mapped, out_par)
    from pylab import rcParams
    rcParams['figure.figsize'] = (16,10)

    kwargs = {"alpha": 0.4}
    for x_axis in x_axis_list:
        plotter.plot_same_orders(ratio_dfs_mapped, mapped=False, x_axis=x_axis, **kwargs)


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
    parser.add_argument("--input_param", nargs='+', type=str, default=None,
            help='''A list of input parameters of the model to plot such
            as "dzeta" or "dbeta_r".''')
    parser.add_argument("--output_param", nargs='+', type=str, default=None,
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
    latexify.set_size(beamer=True)

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

        lon = df[df.in_param == "LONGITUDE"]["deriv"].unique()
        lat = df[df.in_param == "LATITUDE"]["deriv"].unique()
        # Remove entries with coordinates. We just stored them in lat and lon
        df = df[df.in_param != "LONGITUDE"]
        df = df[df.in_param != "LATITUDE"]

        ratio_dfs = {}
        for out_par in args.output_param:
            ratio_dfs[out_par] = ratio_deriv(df, out_par)
        if args.norm:
            norm_deriv(df)

    # Plot data
    if net_df is not None:
        for out_par in args.output_param:
            plot_res_line(df=net_df,
                          out_param=out_par,
                          dots=True,
                          mapped=True)

    if df is not None:
        kwargs = {"alpha": 0.4}
        plot_same_orders(
            ratio_dfs,
            out_params=args.output_param,
            in_params=args.input_param,
            mapped=True,
            **kwargs)
        # Use the following call if you want everything in one plot:
        # plot_ratio_deriv_line(
        #     ratio_dfs,
        #     out_params=args.output_param,
        #     in_params=args.input_param,
        #     mapped=True)



