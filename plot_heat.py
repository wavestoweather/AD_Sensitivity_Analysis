import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import numpy as np
import pandas as pd
import sys
import os
from loader import *


prefix = sys.argv[1]
prefix_title = sys.argv[1]
sns.set_style("darkgrid")


def plot_heat(df_dict, prefix):
    """
    Plot a heatmap where the x-axis consists of timesteps and the y-axis
    are the input parameters of the model. The colormap corresponds to the
    derivative. Plots several heatmaps for each dataframe in df_dict.

    Parameters
    ----------
    df_dict     A dictionary of pandas.Dataframe with key the output parameter
                and the columns the input parameters and values are the
                derivatives.
    prefix      Prefix of the filename where the plot shall be stored.

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
    df_dict     A dictionary of pandas.Dataframe with key the output parameter
                and the columns the input parameters and values are the
                derivatives.
    prefix      Prefix of the filename where the plot shall be stored.

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
    df          A pandas.Dataframe with the columns the output parameters.
    prefix      Prefix of the filename where the plot shall be stored.

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

# def plot_map(xar, prefix, suffix):
#     outpath = out + name


df_results = load_output(prefix=prefix, suffix="_start_over_20160922_00")
plot_line_res(df=df_results, prefix=prefix)
# df_deriv = load_derivatives(prefix=prefix, suffix="_start_over",
                            # filt=True, EPSILON=1e-8)
# plot_line(df_deriv, prefix=prefix)
