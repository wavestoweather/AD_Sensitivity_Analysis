from argparse import ArgumentParser
import numpy as np
import pandas as pd


def norm_derivatives(df):
    """
    Normalize the derivatives for every timestep individually such that
    values are between :math:`[-1, 1]` for big negative and positive impacts.
    .. math:: 
    \text{norm}(x, t) = \text{sign}(x(t)) \cdot \frac{|x(t)| - \min(x(t))}{\max(x(t)) - \min(x(t))})
    
    This is done inplace.

    Parameters
    ----------
    df : pandas.DataFrame
        Columns are "timestep", "trajectory", "out_param", "in_param" and "deriv".

    """
    idx = []

    for traj in df["trajectory"].unique():
        df_traj = df.loc[df.trajectory == traj]
        for out_param in df["out_param"].unique():
            df_param = df_traj.loc[df_traj.out_param == out_param]
            for t in df["timestep"].unique():
                df_t = df_param.loc[df_param.timestep == t]
                deriv_min = df_t["deriv"].min()
                deriv_max = df_t["deriv"].max()
                df_t["deriv"] = ( np.sign(df_t["deriv"])
                                  * ( np.abs(df_t["deriv"]) - deriv_min )
                                  / ( deriv_max - deriv_min ) )
                df.update(df_t)


def rotate(lon, lat, pollat=90, pollon=-180):
    """
    Calculate the rotated coordinates to non-rotated.
    non-rotated lat-lon grid set uses pollat = 90 and pollon = -180.

    Parameters
    ----------
    lon : list of float
        A list of longitudes.
    lat : list of float 
        A list of latitudes.
    pollat : float
        Rotated south pole coordinates.
    pollon : float
        Rotated south pole coordiantes.

    Returns
    -------
    list, list
        List of rotated latitudes and list of rotated longitudes.
    """
    lon1 = lon*np.pi/180.0
    lat1 = lat*np.pi/180.0

    x = np.cos(lon1)*np.cos(lat1)
    y = np.sin(lon1)*np.cos(lat1)
    z = np.sin(lat1)

    theta = -(90-pollat)*np.pi/180.0
    phi = -pollon*np.pi/180.0

    xx =  np.cos(theta)*np.cos(phi)*x + np.sin(phi)*y + np.sin(theta)*np.cos(phi)*z
    yy = -np.cos(theta)*np.sin(phi)*x + np.cos(phi)*y - np.sin(theta)*np.sin(phi)*z
    zz = -np.sin(theta)*x + np.cos(theta)*z

    lon2 = np.arctan2(yy, xx)
    lat2 = np.arcsin(zz)
    lon2 = lon2*180.0/np.pi
    lat2 = lat2*180.0/np.pi
    return lat2, lon2


def transform_df(res, net_path=None, rotate=False):
    """
    Return a new pandas dataframe that adds longitude and latitude to res
    given information in net_df. 

    Parameters
    ----------
    res : pandas.Dataframe
        A pandas dataframe with "deriv", "in_param", 
        "out_parm", "timestep" and "trajectory".
    net_path : string
        Path to a netCDF file with information
        about the longitude and latitude.
    rotate : bool
        Rotate the latitude and longitude.

    Returns
    -------
    pandas.Dataframe
        res with latitude and longitude added.
    """
    new_dic = {"deriv": [], "in_param": [], 
               "out_param": [], "timestep": [], 
               "trajectory": [], "LONGITUDE": [], 
               "LATITUDE": []}
    ds = xr.open_dataset(net_path)
    pollat = ds.attrs["pollat"]
    pollon = ds.attrs["pollon"]
    net_df = ds.to_dataframe()

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
    new_df = pd.DataFrame.from_dict(new_dic)
    if rotate:
        lat, lon = rotate(new_df["LONGITUDE"], new_df["LATITUDE"],
                          pollon=pollon, pollat=pollat)
        new_df["LONGITUDE"] = lon
        new_df["LATITUDE"] = lat
    return new_df


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=str, default=None,
            help='''Path to folder with all files to load.''')
    parser.add_argument("-o", "--output", type=str, default="file.csv",
            help='''Directory and name to save the normed csv.''')
    parser.add_argument("-e", "--epsilon", type=float, default=0.0,
            help='''Value to filter values with. Default is 0.0.''')
    parser.add_argument("-f", "--filter", action="store_true", default=False,
            help='''Filter values smaller than epsilon if used.''')
    parser.add_argument("-t", "--trajectory", nargs='+', type=int, default=1,
            help='''A list of trajectories to consider here. Using a selected
            number of trajectories helps reducing the amount of used RAM.''')
    parser.add_argument("-c", "--csv", type=str, default=None,
            help='''Load an already filtered csv file and just norm the derivatives.
            Deactivates -t, -f and -e.''')
    parser.add_argument("-n", "--netcdf", type=str, default=None,
            help='''Path to a netCDF file that is used to include longitude and latitude.''')
    parser.add_argument("-r", "--rotate", action="store_true", default=False,
            help='''Rotate longitude and latitude coordinates.''')
    args = parser.parse_args()

    if args.csv is None:
        df_data = load_mult_derivates_directory(
                direc=args.input,
                filt=args.filter,
                EPSILON=args.epsilon,
                trajectories=args.trajectory)
    else:
        df_data = pd.read_csv(args.csv)
        # df_data = pd.read_csv(args.csv, nrows=1000000, names=["deriv", "in_param", "out_param", "timestep", "trajectory"], skiprows=100000)
    print("Max and min deriv: {}, {}".format(df_data["deriv"].max(), df_data["deriv"].min()))
    
    norm_derivatives(df_data)
    df_data = df_data.loc[:, ~df_data.columns.str.contains('^Unnamed')]
    print(df_data)
    print("Max and min deriv: {}, {}".format(df_data["deriv"].max(), df_data["deriv"].min()))
    if not args.netcdf is None:
        df_data = transform_df(df_data, args.netcdf, args.rotate)
    df_data.to_csv(args.output)


