from argparse import ArgumentParser
from iris.analysis.cartography import rotate_pole
import json
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset as netset
import numpy as np
import pandas as pd
from PIL import Image
import xarray as xr


class Configuration:
    """
    Configuration for plotting several images. May be extended with
    information such as dpi.

    Attributes
    ----------
    net_params : list of string
        List of output parameters to plot with information from a netCDF
        file from a Cosmo or ICON simulation.
    deriv_params: list of tuples of string
        Each tuple is a combination of input and output parameters to
        plot the corresponding derivative.

    """
    net_params = []
    deriv_params = []

    def __init__(self, path):
        """
        Initiliazation loads a json file as configuration.

        Parameters
        ----------
        path : string
            Json file with dictionary for all attributes.
        """
        self.load_config(path)

    def load_config(self, path):
        """
        Load configuration file and set members. If path is None, a
        default configuration is loaded.

        Parameters
        ----------
         path : string
            Json file with dictionary for all attributes.
        """
        if f is None:
            self.net_params = ["T", "S"]
            self.deriv_params = [("dzeta", "T"), ("dzeta", "S")]
        else:
            with open(path) as f:
                data = json.load(f)
                self.net_params = data["net_params"]
                self.deriv_params = data["deriv_params"]


def plot_weather_deriv(df, in_param, out_param, res, path="pics/"):
    """
    Plot derivatives along a trajectory given an input and output parameter.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns "in_param", "out_param", "LONGITUDE",
        "LATITUDE", "trajectory", "deriv".
    in_param : string
        Input parameter to plot such as "dzeta".
    out_param : string
        Output parameter to plot such as "qr".
    res : string
        Resolution for background map. Needs basemap-data-hires
        for "f" and "h".

    Returns
    -------
    pandas.Dataframe
        Dataframe with added "xm" and "ym" columns used for plotting
        with mercator projection.
    """
    # Get the parts which belong to the in_param and out_param
    df_tmp = df.loc[df["in_param"] == in_param]
    df_tmp = df_tmp.loc[df_tmp["out_param"] == out_param]
    n_traj = df_tmp["trajectory"].max()
    if df_tmp.empty:
        print("Empty dataframe")
        return None
    llon = df_tmp["LONGITUDE"].min()
    ulon = df_tmp["LONGITUDE"].max()
    llat = df_tmp["LATITUDE"].min()
    ulat = df_tmp["LATITUDE"].max()
    # You may change resolution to "h" or "f" for higher resolutions
    # but you need to install basemap-data-hires first
    my_map = Basemap(projection="merc",
                    resolution=res, area_thresh=1000.0,
                    llcrnrlon=llon, llcrnrlat=llat, # min longitude, latitude
                    urcrnrlon=ulon, urcrnrlat=ulat) # max longitude, latitude

    my_map.bluemarble()
    if "xm" not in df_tmp:
        xs, ys = my_map(np.asarray(df_tmp.LONGITUDE),
                                   np.asarray(df_tmp.LATITUDE))
        df_tmp["xm"] = xs.tolist()
        df_tmp["ym"] = ys.tolist()

    color = []
    cmap = matplotlib.cm.get_cmap('Spectral')
    min_vals = []
    max_vals = []
    n_traj2 = 0
    for i in range(n_traj):
        df_tmp2 = df_tmp.loc[df_tmp["trajectory"] == i]
        if df_tmp2.empty:
            continue
        n_traj2 += 1
        mi = df_tmp2["deriv"].min()
        ma = df_tmp2["deriv"].max()
        if ma-mi <= 1e-29:
            print("Traj {} has zero with ma {}, mi {}".format(i, ma, mi))
            ma = ma + 1e-20
        min_vals.append(mi)
        max_vals.append(ma)
    for i, c in enumerate(df_tmp["deriv"]):
        color.append(cmap((c-min_vals[i%n_traj2])
                     / (max_vals[i%n_traj2]-min_vals[i%n_traj2])))
    my_map.scatter(df_tmp["xm"], df_tmp["ym"], color=color, s=6)
    plt.savefig(path + "deriv_" + in_param + "_" + out_param + ".png")
    return df_tmp


def transform_df(df, net_df):
    """
    Add columns LONGITUDE and LATITUDE to df coming from net_df.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe from a csv file with information about derivatives.
    net_df : pandas.Dataframe
        Dataframe from a netCDF file with results from a Cosmo or ICON
        simulation.

    Returns
    -------
    pandas.Dataframe
        Dataframe with columns "deriv", "in_param", "out_param",
        "timestep", "trajectory", "LONGITUDE", "LATITUDE".
    """
    new_dic = {"deriv": [], "in_param": [], "out_param": [],
               "timestep": [], "trajectory": [],
               "LONGITUDE": [], "LATITUDE": []}
    for traj in df.trajectory.unique():
        df_traj = df.loc[df["trajectory"] == traj]
        net_df_traj = net_df.iloc[np.arange(traj, n_rows, n_traj)]
        for out_param in df_traj.out_param.unique():
            df_out = df_traj.loc[df_traj["out_param"] == out_param]
            for in_param in df_out.in_param.unique():
                df_in = df_out.loc[df_out["in_param"] == in_param]
                max_time = df_in["timestep"].max()
                for t in np.arange(20, max_time+1, 20):
                    net_df_time = net_df_traj.loc[net_df_traj["time"] == t]
                    if net_df_time.empty:
                        continue
                    new_dic["in_param"].append(in_param)
                    new_dic["out_param"].append(out_param)
                    new_dic["timestep"].append(t)
                    summed = df_in["deriv"].sum()/20.0
                    new_dic["deriv"].append(summed)
                    new_dic["LATITUDE"].append(net_df_time["lat"][0])
                    new_dic["LONGITUDE"].append(net_df_time["lon"][0])
                    new_dic["trajectory"].append(traj)
    return pd.DataFrame.from_dict(new_dic)


def plot_traj_nc(df, out_param, res, n_traj=903,
                 trajs=None, max_time=None, path="pics/"):
    """
    Plot trajectories from an ICON or Cosmo simulation.

    Parameters
    ----------
    df : pandas.Dataframe
        A dataframe with columns "lon", "lat" and various out parameters
        as columns.
    out_param : string
        Output parameter to use plot as "QR".
    res : string
        Resolution for background map. Needs basemap-data-hires
        for "f" and "h".
    n_traj : int
        Number of trajectories in df.
    trajs : list of int
        List of trajectories to plot.
    max_time : float
        Maximum time to plot in seconds.

    """
    df_tmp = df.copy()
    n_rows = len(df_tmp.index)
    if not trajs is None:
        idx = []
        for i in trajs:
            idx.extend(np.arange(i, n_rows, n_traj))
        df_tmp = df.iloc[idx]
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
                    resolution=res, area_thresh=1000.0,
                    llcrnrlon=llon, llcrnrlat=llat, # min longitude, latitude
                    urcrnrlon=ulon, urcrnrlat=ulat) # max longitude, latitude

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
        color.append(cmap((c-min_vals[i%n_traj])
        / (max_vals[i%n_traj]-min_vals[i%n_traj])))

    my_map.scatter(xs, ys, color=color, s=6)
    plt.savefig(path + "nc_" + out_param + ".png")


if __name__ == "__main__":
    parser = ArgumentParser(description=
            '''Plot data from a filtered csv dataset
            and from a netcdf file.''')
    parser.add_argument("-c", "--csv", type=str, default=None,
            help='''File with derivatives.''')
    parser.add_argument("-n", "--netcdf", type=str, default=None,
            help='''File with netCDF info from a Cosmo or Icon simulation.''')
    parser.add_argument("-o", "--outdir", type=str, required=True,
            help='''Path to store plots and transformed data if available.''')
    parser.add_argument("-r", "--rotate", action="store_true", default=False,
            help='''Use this if longitudes and latitudes have to be rotated
            on the csv file.''')
    parser.add_argument("-i", "--input_config", type=str, default=None,
            help='''Config file that tells which plots to generate.
            Default plots are S and T with dzeta each.''')
    parser.add_argument("-t", "--trajectory", nargs='+', type=int, default=1,
            help='''A list of trajectories to consider here. Using a selected
            number of trajectories makes plots easier to read. If -1 is used,
            all trajectories will be plotted.''')
    parser.add_argument("--resolution", type=str, default="f",
            help='''Resolution for the background maps. Default is f.
            You may need to install basemap-data-hires first for this or
            choose a lower resolution.''')
    parser.add_argument("-m", "--max_time", type=float, default=None,
            help='''Maximum time to plot in seconds.''')
    args = parser.parse_args()

    matplotlib.rcParams['figure.figsize'] = (14, 10)
    config = Configuration(args.input_config)

    if args.netcdf is not None:
        ds = xr.open_dataset(args.netcdf)
        pollat = ds.attrs["pollat"]
        pollon = ds.attrs["pollon"]
        n_traj = len(ds.ntra)
        n_timesteps = len(ds.ntim)
        net_df = ds.to_dataframe()
        if args.rotate:
            lat, lon = rotate_pole(
                np.asarray(net_df["lon"].tolist()),
                np.asarray(net_df["lat"].tolist()),
                pole_lon=pollon,
                pole_lat=pollat)
            net_df["lon"] = lon
            net_df["lat"] = lat
            # TODO: Save to csv if rotating takes too long.
            # res.to_csv(args.outdir + "net_cdf_with_coords_rotated.csv")
        for out_param in config.net_params:
            plot_traj_nc(net_df, out_param, res=args.resolution,
                         n_traj=n_traj, trajs=args.trajectory,
                         max_time=args.max_time, path=args.outdir)

    if args.csv is not None:
        res = pd.read_csv(args.csv, low_memory=True)

        if args.netcdf is not None and not "LONGITUDE" in res:
            res = transform_df(res, net_df)
            res.to_csv(args.outdir + "_with_coords.csv")

        if args.rotate and args.netcdf is not None:
            lat, lon = rotate_pole(
                np.asarray(res["LONGITUDE"].tolist()),
                np.asarray(res["LATITUDE"].tolist()),
                pole_lon=pollon,
                pole_lat=pollat)
            res["LONGITUDE"] = lon
            res["LATITUDE"] = lat
            res.to_csv(args.outdir + "_with_coords_rotated.csv")
        for in_param, out_param in config.deriv_params:
            plot_weather_deriv(res, in_param, out_param, res=args.resolution,
                               path=args.outdir)