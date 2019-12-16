import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import numpy as np
import pandas as pd
from scipy import stats
import sys
import os
from loader import *


# prefix = sys.argv[1]
# prefix_title = sys.argv[1]
sns.set_style("darkgrid")


def correlation(data, out_param):
    """
    Calculate the correlation for a given out_param along all (important)
    trajectories for each in_param.

    Parameters
    ----------
    data        A pandas dataframe with columns
                trajectory, timestep, out_param, in_param, deriv

    Return
    ------
    A dictionary for each in_param as keys and a correlation matrix.
    """
    corr = {}
    tmp_data = data.loc[data["out_param"] == out_param]
    trajectories = data["trajectory"].unique()
    for in_param in data["in_param"].unique():
        in_data = tmp_data.loc[tmp_data["in_param"] == in_param]
        corr_matrix = np.zeros(shape=(len(trajectories), len(trajectories)))
        corr_matrix.fill(np.nan)
        for i, traj1 in enumerate(trajectories):
            in_data1 = in_data.loc[in_data["trajectory"] == traj1]
            if in_data1.empty:
                continue

            for j, traj2 in enumerate(trajectories):
                in_data2 = in_data.loc[in_data["trajectory"] == traj2]
                if in_data2.empty:
                    continue
                # Remove NaNs
                if(np.isnan(in_data1["deriv"]).any()
                   or np.isnan(in_data2["deriv"]).any()):

                    mask = np.logical_and(~np.isnan(in_data1["deriv"]),
                                         ~np.isnan(in_data2["deriv"]))
                    d1 = in_data1["deriv"].tolist()
                    d2 = in_data2["deriv"].tolist()
                    mask = np.logical_and(~np.isnan(d1),
                                         ~np.isnan(d2))
                    # print(np.shape(in_data1["deriv"]))
                    # print(np.shape(in_data2["deriv"]))
                    # print("mask: \n{}".format(mask))
                    # print("another:\n{}".format(np.isnan(in_data1["deriv"])))
                    # print(np.shape(mask))
                    in1 = in_data1["deriv"][mask]
                    in2 = in_data2["deriv"][mask]
                    r, p = stats.pearsonr(in1, in2)
                    corr_matrix[i, j] = r
                else:
                    r, p = stats.pearsonr(in_data1["deriv"], in_data2["deriv"])
                    corr_matrix[i, j] = r
        corr[in_param] = corr_matrix
    return corr


def chi_squared(data, out_param):
    """
    Do a chi-squared test for a given out_param along all (important)
    trajectories for each in_param.

    Parameters
    ----------
    data        A pandas dataframe with columns
                trajectory, timestep, out_param, in_param, deriv

    Return
    ------
    A dictionary for each in_param as keys and a tupel with a matrix of
    chi_squared values and the p-values as result.
    """
    chi = {}
    tmp_data = data.loc[data["out_param"] == out_param]
    trajectories = data["trajectory"].unique()
    for in_param in data["in_param"].unique():
        in_data = tmp_data.loc[tmp_data["in_param"] == in_param]
        chi_matrix = np.zeros(shape=(len(trajectories), len(trajectories)))
        p_matrix = np.zeros(shape=(len(trajectories), len(trajectories)))
        chi_matrix.fill(np.nan)
        p_matrix.fill(np.nan)
        for i, traj1 in enumerate(trajectories):
            in_data1 = in_data.loc[in_data["trajectory"] == traj1]
            if in_data1.empty:
                continue
            for j, traj2 in enumerate(trajectories):
                in_data2 = in_data.loc[in_data["trajectory"] == traj2]
                if in_data2.empty:
                    continue
                # Bin the data
                f_obs, edges = np.histogram(in_data1["deriv"], bins="auto", density=False)
                f_exp, edges = np.histogram(in_data2["deriv"], bins=len(edges)-1, density=False)
                chisq, p = stats.chisquare(f_obs=f_obs, f_exp=f_exp)
                # print("{}\n{} \nvs\n{}\n\n".format(in_param, in_data1["deriv"], in_data2["deriv"]))
                chi_matrix[i, j] = chisq
                p_matrix[i, j] = p
        chi[in_param] = (chi_matrix, p_matrix)
    return chi


def plot_corr(dic):
    # print(dic)
    for out_param in dic["out_param"].unique():
        out_dic = dic.loc[dic["out_param"] == out_param]
        for in_param in out_dic["in_param"].unique():

            in_dic = out_dic.loc[out_dic["in_param"] == in_param]
            # print("out {}, in {}, \n{}".format(out_param, in_param, in_dic))
            if in_dic.empty:
                continue
            piv = in_dic.pivot("trajectory", "timestep", "deriv")
            ax = sns.heatmap(piv, cmap="viridis", cbar=True)
            plt.savefig("pics/corr/corr_" + out_param + "_" + in_param + ".png",
                        dpi=300)
            plt.close()


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


# res = load_mult_derivates("dataMogon/sb_ice_wcb7200_traj",
#                           "start_over_20160922_00", filt=True, lo=1, hi=865,
#                           EPSILON=1e-30, high=None)
# res.to_pickle("sb_ice_wcb7200_traj_start_over_20160922_00.gzip")
# res.to_csv("sb_ice_wcb7200_traj_no_filter_start_over_20160922_00.csv")
res = pd.read_csv("sb_ice_wcb272280_filt_zero_30_start_over_20160922_00.csv")
# res = pd.read_csv("sb_ice_wcb7200_trajAll_start_over_20160922_00.csv")
# with pd.option_context('display.max_rows', None, 'display.max_columns', None):
#     print(res[0:40])

def plot_avg_df(df):
    """
    Plot a heatmap with summed and normalized derivatives where the
    x-axis is the trajectories,
    the y-axis are the input parameters and the color indicates the
    size of the derivatives. Plots for each output variable another heatmap.


    """
    out_params = df["out_param"].unique()
    for param in out_params:
        df_param = df.loc[df["out_param"] == param]
        dic_plot = {"trajectory": [], "param": [], "deriv": []}
        in_params = df_param["in_param"].unique()
        trajecs = df_param["trajectory"].unique()
        for traj in trajecs:
            df_traj = df_param.loc[df_param["trajectory"] == traj]
            max_val = 0.0
            tmp_param = ""
            # Search max value
            for in_param in in_params:
                dic_plot["param"].append(in_param)
                dic_plot["trajectory"].append(traj)
                if in_param in df_traj["in_param"].values:
                    s = df_traj.loc[df_traj["in_param"] == in_param]["deriv"].sum()
                    if max_val == 0.0 or max_val < abs(s):
                        max_val = abs(s)
                        tmp_param = in_param
            # print("out {}, in {}, max {}".format(param, tmp_param, max_val))
            # Normalize the data and add it
            for in_param in in_params:
                if in_param in df_traj["in_param"].values:
                    s = df_traj.loc[df_traj["in_param"] == in_param]["deriv"].sum()
                    dic_plot["deriv"].append(s/max_val)
                else:
                    dic_plot["deriv"].append(np.nan)

        df_plot = pd.DataFrame(dic_plot)
        # with pd.option_context('display.max_rows', None, 'display.max_columns', None):

        #     print(df_plot)
        df_plot = df_plot.pivot("param", "trajectory", "deriv")
        ax = sns.heatmap(df_plot, cmap="viridis")
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(8)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(8)
        plt.title("WCB: " + param)
        plt.tight_layout()
        print("Save to pics/deriv_filter_zero_30_{}.png".format(param))
        plt.savefig("pics/deriv_filter_zero_30_" + param + ".png", dpi=300)
        plt.close()

    return True

# Remove all entries with "min"
print(res)
res = res[~res.in_param.str.contains("min")]
print(res)
plot_avg_df(res)
# print(res)

### Idea
### Get the sum of the derivatives for each trajectory (perhaps abs())
### Compare that between all trajectories
# deriv1                        | abs sum of derivatives
# deriv2                        |
# deriv3                        |
#           traj1 traj2 traj3
### Another plot: For each derivative
# traj1                 | derivative value
# traj2                 |
# traj3                 |
#           t0 t1 t2 t3
### Problem: Are those values comparable? Perhaps ratio against (abs) highest
### derivative of each trajectory is a better proxy

# print(res.loc[res["out_param"] == "p"])
# res = load_mult_derivates("/data/project/wcb/sb_ice_warming_traj",
#                           "start_over", filt=True, lo=0, hi=3)
# chi = chi_squared(res, "p")
# print(chi)
# for out_param in res["out_param"].unique():
#     corr = correlation(res, out_param)
#     # Remove all matrices with NaN
#     to_remove = []
#     for in_param in corr.keys():
#         if np.isnan(corr[in_param][0, 0]):
#             to_remove.append(in_param)
#     for i in to_remove:
#         del corr[i]
#     print("\nCorrelation matrices for {}".format(out_param))
#     print(corr)
# plot_corr(res)
