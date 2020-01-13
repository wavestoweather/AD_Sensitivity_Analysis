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


def correlation(data, out_param):
    """
    Calculate the correlation for a given out_param along all (important)
    trajectories for each in_param.

    Parameters
    ----------
    data : pandas.Dataframe
        A pandas dataframe with columns
        trajectory, timestep, out_param, in_param, deriv
    out_param : string
        Output parameter to get the correlation matrix for such as "qr".

    Returns
    -------
    dic of 2D np.ndarrays
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
    data : pandas.Dataframe
        A pandas dataframe with columns
        trajectory, timestep, out_param, in_param, deriv
    out_param : string
        Output parameter to get the chi-squared matrix such as "qr".

    Returns
    ------
    dic of tuples of 2D np.ndarrays
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
    """
    Plot the correlation matrix.

    Parameters
    ----------
    dic : dic
        Dictionary from the function "correlation".

    """
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
    df_dict : dic of pandas.Dataframe
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
    df_dict : dic of pandas.Datafram
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
        _, ax = plt.subplots()

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


def plot_avg_df(df):
    """
    Plot a heatmap with summed and normalized derivatives where the
    x-axis is the trajectories,
    the y-axis are the input parameters and the color indicates the
    size of the derivatives. Plots for each output variable another heatmap.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns "out_param", "in_param", "trajectory", "deriv".

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


def plot_weather_deriv(df, in_param, out_param, filename=None):
    """
    Plot the derivatives of the trajectories.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns "out_param", "in_param", "trajectory", "deriv",
        "LONGITUDE" and "LATITUDE".
    in_param : String
        Input parameter of the simulation to plot, such as "dzeta" or "dbeta_r".
    out_param : String
        Output parameter of the simulation to plot, such as "S", "T", "qr".
    filename : string
        Where to store the image.

    Returns
    -------
    pandas.Dataframe
        df with additional columns "xm" and "ym" for projected coordinates.
    """
    # Get the parts which belong to the in_param and out_param
    df_tmp = df.loc[df["in_param"] == in_param]
    df_tmp = df_tmp.loc[df_tmp["out_param"] == out_param]

    if df_tmp.empty:
        print("Empty dataframe")
        return None
    n_traj = df_tmp["trajectory"].nunique()
    llon = df_tmp["LONGITUDE"].min()# - 5
    ulon = df_tmp["LONGITUDE"].max() #+ 5
    llat = df_tmp["LATITUDE"].min() #- 5
    ulat = df_tmp["LATITUDE"].max() #+ 5
    # You may change resolution to "h" or "f" for higher resolutions
    # but you need to install basemap-data-hires first
    my_map = Basemap(projection="merc",
                    resolution="f", area_thresh=1000.0,
                    llcrnrlon=llon, llcrnrlat=llat, # min longitude and latitude
                    urcrnrlon=ulon, urcrnrlat=ulat) # max longitude and latitude

    my_map.bluemarble()
    if "xm" not in df_tmp:
        xs, ys = my_map(np.asarray(df_tmp.LONGITUDE), np.asarray(df_tmp.LATITUDE))
        df_tmp["xm"] = xs.tolist()
        df_tmp["ym"] = ys.tolist()

    color = []
    cmap = matplotlib.cm.get_cmap('Spectral')
    min_vals = []
    max_vals = []
    for i in range(n_traj):
        df_tmp2 = df_tmp.loc[df_tmp["trajectory"] == i]
        mi = df_tmp2["deriv"].min()
        ma = df_tmp2["deriv"].max()
        if ma-mi <= 1e-29:
            print("Traj {} has zero with ma {}, mi {}".format(i, ma, mi))
            ma = ma + 1e-20
        min_vals.append(mi)
        max_vals.append(ma)


    for i, c in enumerate(df_tmp["deriv"]):
        color.append(cmap((c-min_vals[i%n_traj]) / (max_vals[i%n_traj]-min_vals[i%n_traj])))

    my_map.scatter(df_tmp["xm"], df_tmp["ym"], color=color, s=6)
    if filename is not None:
        plt.savefig(filename, dpi=300)

    return df_tmp


def plot_traj_nc(df, out_param, filename=None, trajs=None, max_time=None):
    """
    Plot the trajectories of a netCDF file.

    Parameters
    ----------
    df : pandas.Dataframe
         Dataframe with columns "lon" and "lat" and a column for each out_param.
    out_param : String
        Parameter to plot such as "T" or "qr".
    filename : string
        Where to store the image.
    trajs : list of int
        List of trajectory indices to plot.
    max_time : float
        Maximum time to plot.
    """
    df_tmp = df.copy()
    n_traj = 903
    n_rows = len(df_tmp.index)
    if not trajs is None:
        idx = []
        for i in trajs:
            idx.extend(np.arange(i, n_rows, n_traj))
        df_tmp = df.iloc[idx]
#     df_tmp = df.dropna()
    if not max_time is None:
        df_tmp = df_tmp.loc[df_tmp["time"] <= max_time]
    n_rows = len(df_tmp.index)

    llon = df_tmp["lon"].min() -5
    ulon = df_tmp["lon"].max() + 5
    llat = df_tmp["lat"].min() - 5
    ulat = df_tmp["lat"].max() + 5
    # You may change resolution to "h" or "f" for higher resolutions
    # but you need to install basemap-data-hires first
    my_map = Basemap(projection="merc",
                    resolution="f", area_thresh=1000.0,
                    llcrnrlon=llon, llcrnrlat=llat, # min longitude and latitude
                    urcrnrlon=ulon, urcrnrlat=ulat) # max longitude and latitude

    my_map.bluemarble()
    xs, ys = my_map(np.asarray(df_tmp.lon), np.asarray(df_tmp.lat))
    xs = xs.tolist()
    ys = ys.tolist()

    color = []
    cmap = matplotlib.cm.get_cmap('Spectral')
    min_vals = []
    max_vals = []
    if not trajs is None:
        n_traj = len(trajs)
    for i in range(n_traj):
        df_tmp2 = df_tmp.iloc[np.arange(i, n_rows, n_traj)]
        min_vals.append(df_tmp2[out_param].min())
        max_vals.append(df_tmp2[out_param].max())

    for i, c in enumerate(df_tmp[out_param]):
        color.append(cmap((c-min_vals[i%n_traj]) / (max_vals[i%n_traj]-min_vals[i%n_traj])))

    my_map.scatter(xs, ys, color=color, s=6)
    if filename is not None:
        plt.savefig(filename, dpi=300)


if __name__ == "__main__":
    parser = ArgumentParser(description=
            '''Load derivatives of a simulation or a netCDF file and plot
            trajectories over satellite images.''')
    parser.add_argument("-c", "--csv", type=str, default=None,
            help='''Load an already processed csv file with derivatives.''')
    parser.add_argument("-n", "--netcdf", type=str, default=None,
            help='''Path to a netCDF file. If you want to plot derivatives from
            a csv file that lacks "LONGITUDE" and "LATITUDE", you need to
            provide a netCDF file with that information.''')
    parser.add_argument("-o", "--output", type=str, default="pics/traj",
            help='''Directory and name to save the image.''')
    parser.add_argument("-t", "--trajectory", nargs='+', type=int, default=1,
            help='''A list of trajectories to consider here. Using a selected
            number of trajectories helps reducing the amount of used RAM
            and makes it easier to see anything in the image.''')
    parser.add_argument("-r", "--rotate", action="store_true", default=False,
            help='''Rotate longitude and latitude coordinates.''')
    parser.add_argument("--out_transformed", type=str, default="transformed.csv",
            help='''In case the csv file had to be transformed because of missing
            coordinates, store this as a new csv file. Set to "None" if you
            don't want to store it.''')
    parser.add_argument("--input_param", nargs='+', type=str, default="dzeta",
            help='''A list of input parameters of the model to plot such
            as "dzeta" or "dbeta_r".''')
    parser.add_argument("--output_param", nargs='+', type=str, default="T",
            help='''A list of output parameters of the model to plot such
            as "T", "qr" or "p".''')
    parser.add_argument("--max_time", type=float, default=None,
            help='''Maximum time in seconds to plot.''')
    args = parser.parse_args()

    try:
        from loader import *
    except:
        from scripts.loader import *
    from pylab import rcParams
    rcParams['figure.figsize'] = (14,10)

    if args.csv is not None:
        res = load_output(args.csv)
        if "LONGITUDE" not in res:
            if args.netcdf is None:
                print("csv file has no coordinates. Please specify a netcdf"
                    + " file with the necessary information. Aborting.")
                exit()
            if args.rotate:
                net_df, pollat, pollon = load_nc(
                                        inp=args.netcdf,
                                        get_pol=True)
                rotate_df(net_df, pollon, pollat, "lon", "lat")
            else:
                net_df = load_nc(
                                 inp=args.netcdf,
                                 get_pol=False)
            res = transform_df2(res, net_df)
            if args.out_transformed != "None":
                res.to_csv(args.out_transformed)
        for input_param in args.input_param:
            for output_param in args.output_param:
                filename = (args.output + "_deriv_" + out_param
                            + "_" + input_param + ".png")
                plot_weather_deriv(df=res,
                                   in_param=input_param,
                                   out_param=output_param,
                                   filename=filename)

    if net_df is None and args.netcdf is not None:
        # plot only the netCDF file
        if args.rotate:
            net_df, pollat, pollon = load_nc(
                                    inp=args.netcdf,
                                    get_pol=True)
            rotate_df(net_df, pollon, pollat, "lon", "lat")
        else:
            net_df = load_nc(
                                inp=args.netcdf,
                                get_pol=False)
    if net_df is not None:
        for out_param in args.output_param:
            filename = args.output + "_netcdf_" + out_param + ".png"
            plot_traj_nc(df=net_df,
                         out_param=out_param,
                         filename=filename,
                         trajs=args.trajectory,
                         max_time=args.max_time)

    # Example:
    # res = load_output("sb_ice_wcb272280_filt_zero_30_start_over_20160922_00_traj1_all.csv")
    # net_df, pollat, pollon = load_nc(
    #                                  inp="O_WCB_all_20160922_00.nc",
    #                                  get_pol=True)
    # rotate_df(res, pollon, pollat)
    # rotate_df(net_df, pollon, pollat, "lon", "lat")
    # res = transform_df2(res, net_df)
    # res.to_csv("transformed.csv") # Transformation takes long so let's store it
    # df_dzeta_S = plot_weather_deriv(res, "dzeta", "S")
    # plot_traj_nc(net_df, "T")
    # plot_traj_nc(net_df, "T", np.arange(0,20))


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
