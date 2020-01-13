from argparse import ArgumentParser
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import numpy as np
import pandas as pd
import sys
import os
try:
    from loader import *
except:
    from scripts.loader import *


def plot_heat(df_dict, prefix):
    """
    Plot a heatmap where the x-axis consists of timesteps and the y-axis
    are the input parameters of the model. The colormap corresponds to the
    derivative. Plots several heatmaps for each dataframe in df_dict.

    Parameters
    ----------
    df_dict : pandas.Dataframe
        A dictionary of pandas.Dataframe with key the output parameter
        and the columns the input parameters and values are the
        derivatives.
    prefix : string
        Prefix of the filename where the plot shall be stored.

    """
    for out_param in df_dict.keys():
        df = df_dict[out_param]
        fig, ax = plt.subplots()
        # yticklabels = np.arange(0, len(df.columns))
        ax = sns.heatmap(df.transpose(), cmap="viridis", cbar=True, ax=ax)
        i = 0
        save = ("pics/" + prefix + "_heat_" + out_param
                + "_" + "{:03d}".format(i) + ".png")
        while os.path.isfile(save):
            i = i+1
        save = ("pics/" + prefix + "_heat_" + out_param
                + "_" + "{:03d}".format(i) + ".png")

        print("Saving to " + save)
        plt.savefig(save, dpi=300)
        plt.close()


def plot_line(df_dict, prefix):
    """
    Plot a lineplot where the x-axis consists of timesteps and the y-axis
    are the derivatives of the input parameters of the model.
    The colormap corresponds to the derivative, where brighter colors are
    overall more significant. Plots several lineplots for each dataframe in
    df_dict.

    Parameters
    ----------
    df_dict : pandas.Dataframe
        A dictionary of pandas.Dataframe with key the output parameter
        and the columns the input parameters and values are the
        derivatives.
    prefix : string
        Prefix of the filename where the plot shall be stored.

    """
    for out_param in df_dict.keys():
        df = df_dict[out_param]
        if df.empty:
            continue
        x_ticks = np.arange(0, df.timestep.unique().max() + 19, 20)
        fig, ax = plt.subplots()

        ax = sns.lineplot(x="timestep", y="deriv",
                          data=df, hue="param", ax=ax)
        ax.set_title("Deriv. of {}".format(out_param))
        ax.set_xticks(x_ticks)
        i = 0
        save = ("pics/" + prefix + "_line_" + out_param
                + "_" + "{:03d}".format(i) + ".png")
        while os.path.isfile(save):
            i = i+1
            save = ("pics/" + prefix + "_line_" + out_param
                    + "_" + "{:03d}".format(i) + ".png")

        print("Saving to " + save)
        plt.savefig(save, dpi=300)
        plt.close()


def plot_line_res(df, prefix):
    """
    Plot a lineplot where the x-axis consists of timesteps and the y-axis
    are the absolute values of the output parameters.

    Parameters
    ----------
    df : pandas.Dataframe
        A pandas.Dataframe with the columns the output parameters.
    prefix : string
        Prefix of the filename where the plot shall be stored.

    """
    x_ticks = np.arange(0, df.timestep.unique().max() + 19, 20)
    for out_param in df:
        fig, ax = plt.subplots()
        ax = sns.lineplot(x="timestep", y=out_param, data=df, ax=ax)
        ax.set_title("Results of {}".format(out_param))
        ax.set_xticks(x_ticks)
        i = 0
        save = ("pics/" + prefix + "_line_res_" + out_param
                + "_" + "{:03d}".format(i) + ".png")
        while os.path.isfile(save):
            i = i+1
            save = ("pics/" + prefix + "_line_res_" + out_param
                    + "_" + "{:03d}".format(i) + ".png")

        print("Saving to " + save)
        plt.savefig(save, dpi=300)
        plt.close()


if __name__ == "__main__":
    parser = ArgumentParser(description=
            '''Load data given a prefix and suffix of the filenames and plot various figures, mostly heatmaps.''')
    parser.add_argument("-p", "--prefix", type=str, required=True,
            help='''Prefix of the files to load such as "data/sb_ice_traj".''')
    parser.add_argument("-s", "--suffix", type=str, required=True,
            help='''Suffix of the files to load without datatype such as "_start_over_20160922_00"
            or "_start_over" for derivatives.''')
    parser.add_argument("-i", "--image", type=str, default="results",
            help='''Which plot to generate. Can be either "results" or "derivatives" to plot lineplots
            of either the output parameters or the derivatives.''')
    parser.add_argument("-e", "--epsilon", type=float, default=1e-30,
            help='''If "filter" is set to true, use this as filter threshold.''')
    parser.add_argument("-f", "--filter", action="store_true", default=True,
            help='''If you load derivatives, use this to filter derivatives smaller than "epsilon".''')
    args = parser.parse_args()

    sns.set_style("darkgrid")

    if args.image == "results":
        df_results = load_output(prefix=args.prefix, suffix=args.suffix)
        plot_line_res(df=df_results, prefix=args.prefix)
    elif args.image == "derivatives":
        df_deriv = load_derivatives(prefix=args.prefix, suffix=args.suffix,
                                    filt=args.filter, EPSILON=args.epsilon)
        plot_line(df_deriv, prefix=args.prefix)
