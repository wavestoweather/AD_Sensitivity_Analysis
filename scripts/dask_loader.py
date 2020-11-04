# from iris.analysis.cartography import rotate_pole
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import dask.dataframe as pd
from progressbar import progressbar as pb
import sys
import os
import xarray as xr
from multiprocessing import Pool
from itertools import repeat
from glob import glob

params_dict = {"p": "_diff_0.txt", "T": "_diff_1.txt",
               "w": "_diff_2.txt", "S": "_diff_3.txt", "qc": "_diff_4.txt",
               "qr": "_diff_5.txt", "qv": "_diff_6.txt", "Nc": "_diff_7.txt",
               "Nr": "_diff_8.txt", "Nv": "_diff_9.txt",
               "qi": "_diff_10.txt", "Ni": "_diff_11.txt",
               "vi": "_diff_12.txt", "qs": "_diff_13.txt",
               "Ns": "_diff_14.txt", "qg": "_diff_15.txt",
               "Ng": "_diff_16.txt", "qh": "_diff_17.txt",
               "Nh": "_diff_18.txt", "qiout": "_diff_19.txt",
               "qsout": "_diff_20.txt", "qrout": "_diff_21.txt",
               "qgout": "_diff_22.txt",
               "qhout": "_diff_23.txt", "latent_heat": "_diff_24.txt",
               "latent_cool": "_diff_25.txt"}
params_dict2 = {"_diff_0.txt": "p",
               "_diff_1.txt": "T",
               "_diff_2.txt": "w",
               "_diff_3.txt": "S",
               "_diff_4.txt": "qc",
               "_diff_5.txt": "qr",
               "_diff_6.txt": "qv",
               "_diff_7.txt": "Nc",
               "_diff_8.txt": "Nr",
               "_diff_9.txt": "Nv",
               "_diff_10.txt": "qi",
               "_diff_11.txt": "Ni",
               "_diff_12.txt": "vi",
               "_diff_13.txt": "qs",
               "_diff_14.txt": "Ns",
               "_diff_15.txt": "qg",
               "_diff_16.txt": "Ng",
               "_diff_17.txt": "qh",
               "_diff_18.txt": "Nh",
               "_diff_19.txt": "qiout",
               "_diff_20.txt": "qsout",
               "_diff_21.txt": "qrout",
               "_diff_22.txt": "qgout",
               "_diff_23.txt": "qhout",
               "_diff_24.txt": "latent_heat",
               "_diff_25.txt": "latent_cool",
               "_diff_26.txt": "Niout",
               "_diff_27.txt": "Nsout",
               "_diff_28.txt": "Nrout",
               "_diff_29.txt": "Ngout",
               "_diff_30.txt": "Nhout",
               "_diff_31.txt": "z",
               "_diff_32.txt": "n_inact",
               "_diff_33.txt": "depo",
               "_diff_34.txt": "sub"}


def filter_zeros(df_dict, EPSILON=1e-31):
    """
    Drop all columns that have zero impact

    Parameters
    ----------
    df_dict : Dictionary of pandas.Dataframe
        A dictionary of pandas.Dataframe with key the output parameter
        and the columns the input parameters and values are the
        derivatives.
    EPSILON : float
        If all values in a column are lower than EPSILON,
        drop that one.

    Returns
    -------
    Dictionary of pandas.Dataframe
        Modified df_dict with removed columns in each dataframe where only (near)
        zeros existed before.
    """
    key_drop = []
    ignore_keys = ["timestep", "trajectory", "LONGITUDE", "LATITUDE", "MAP"]
    for key in df_dict:
        to_drop = []
        for column in df_dict[key]:
            if column in ignore_keys:
                continue
            if not (abs(df_dict[key][column]) > abs(EPSILON)).any():
                to_drop.append(column)
        if len(to_drop) > 0:
            df_dict[key] = df_dict[key].drop(columns=to_drop)
            print("Dropped {} columns for {}. Shape: {}".format(
                len(to_drop), key, np.shape(df_dict[key])))
            if (df_dict[key].empty or
                (len(df_dict[key].columns) == 1
                 and df_dict[key].columns[0] == "timestep")):
                key_drop.append(key)

    for key in key_drop:
        print("Dropping {} entirely.".format(key))
        del df_dict[key]
    return df_dict


def filter_high(df_dict, high=1e1):
    """
    Drop all columns that have a suspiciously high impact.

    Parameters
    ----------
    df_dict : Dictionary of pandas.Dataframe
        A dictionary of pandas.Dataframe with key the output parameter
        and the columns the input parameters and values are the
        derivatives.
    high : float
        If some values in a column are higher than high,
        drop that one.

    Returns
    -------
    Dictionary of pandas.Dataframe
        Modified df_dict with removed columns in each dataframe where some values
        were too high before.
    """
    key_drop = []
    for key in df_dict:
        to_drop = []
        for column in df_dict[key]:
            if column == "timestep" or column == "trajectory":
                continue
            if (abs(df_dict[key][column]) >= abs(high)).any():
                to_drop.append(column)
        if len(to_drop) > 0:
            df_dict[key] = df_dict[key].drop(columns=to_drop)
            print("Dropped {} columns for {} (too high values). Shape: {}".format(
                len(to_drop), key, np.shape(df_dict[key])))
            if (df_dict[key].empty or
                (len(df_dict[key].columns) == 1
                 and df_dict[key].columns[0] == "timestep")):

                print("Dropping {} entirely (too high values).".format(key))
                key_drop.append(key)

    for key in key_drop:
        del df_dict[key]
    return df_dict


def transform_df(df):
    """
    Create a new pandas.DataFrame with column "param", "timestep", "deriv"
    that can be used for plotting with seaborn.lineplot.
    Optionally adds columns LONGITUDE, LATITUDE, MAP if available in df.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe where columns are the names of the parameters.

    Returns
    -------
    pandas.Dataframe
        Transformed Dataframe.
    """
    if "MAP" in df:
        dicti = {"param": [], "timestep": [], "deriv": [], "MAP": [],
                 "LONGITUDE": [], "LATITUDE": []}
    else:
        dicti = {"param": [], "timestep": [], "deriv": []}

    key_list = ["timestep", "trajectory", "LONGITUDE", "LATITUDE", "MAP"]
    for key in df:
        if key in key_list:
            continue
        dicti["timestep"].extend(df["timestep"].tolist())
        dicti["deriv"].extend(df[key].tolist())
        dicti["param"].extend([key for i in range(len(df["timestep"]))])
        if "MAP" in df:
            dicti["MAP"].extend(df["MAP"].tolist())
            dicti["LONGITUDE"].extend(df["LONGITUDE"].tolist())
            dicti["LATITUDE"].extend(df["LATITUDE"].tolist())
    return pd.DataFrame(dicti)


def transform_df2(df, net_df, n_traj=903, traj_timestep=20):
    """
    Create a new pandas.DataFrame with column "deriv", "in_param", "out_param",
    "timestep", "trajectory", "LONGITUDE" and "LATITUDE".
    that can be used for plotting with plot_many_traj.plot_weather_deriv(..).
    It is mainly used to sum the derivatives for traj_timestep and to add
    coordinates to the data from df.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of derivatives.
    net_df : pandas.DataFrame
        Dataframe of a netCDF file, where the coordinates had been rotated.
    n_traj : int
        Number of trajectories in the netCDF dataframe.
    traj_timestep : int or float
        Timestep size in seconds of the trajectories in the netCDF dataframe.

    Returns
    -------
    pandas.Dataframe
        Transformed Dataframe.
    """
    n_rows = len(net_df.index)
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
                for t in np.arange(traj_timestep, max_time+1, traj_timestep):
                    net_df_time = net_df_traj.loc[net_df_traj["time"] == t]
                    if net_df_time.empty:
                        continue
                    new_dic["in_param"].append(in_param)
                    new_dic["out_param"].append(out_param)
                    new_dic["timestep"].append(t)
                    summed = df_in["deriv"].sum()/traj_timestep
                    new_dic["deriv"].append(summed)
                    new_dic["LATITUDE"].append(net_df_time["lat"][0])
                    new_dic["LONGITUDE"].append(net_df_time["lon"][0])
                    new_dic["trajectory"].append(traj)
    return pd.DataFrame.from_dict(new_dic)


def load_mult_derivates_direc_dic(direc="", parquet=True, netcdf=False,
    columns=None, file_ending="*.nc_wcb"):
    """
    Create a dictionary with out parameters as keys and dictionaries with columns:
    trajectory, timestep, MAP, LATITUDE, LONGITUDE
    and a column for each in parameter such as "da_1", "da_2", "dsnow_alfa_q", ...
    Out parameters is a string such as 'p', 'T', 'w' etc
    MAP is an optionally
    available flag for interesting timesteps.

    Parameters
    ----------
    direc : string
        A path to a directory with files to read.
    parquet : boolean
        If true: Load a series of preprocessed parquet files, else try loading *.txt files.
        If this fails, netcdf-files are loaded.
    columns : list of strings
        Specify the columns to load.
    file_ending : string
        In case of netcdf-files, specify the file ending here to load
        multiple files (takes long) or a single file.

    Returns
    -------
    Dask.Dataframe
        A delayed dataframe.
    """
    if parquet:
        df = pd.read_parquet(direc + "/", columns=columns)
    elif netcdf:
        if "*" in file_ending:
            files = sorted(glob(direc + "/" + file_ending))
            df = None
            for f in files:
                ds = xr.open_dataset(f, decode_times=False)
                if "Output Parameter" in ds:
                    if df is not None:
                        df = df.append(ds.to_dask_dataframe(
                            dim_order=["Output Parameter", "ensemble", "trajectory", "time"]))
                    else:
                        df = ds.to_dask_dataframe(
                            dim_order=["Output Parameter", "ensemble", "trajectory", "time"])
                else:
                    if df is not None:
                        df = df.append(ds.to_dask_dataframe(
                            dim_order=["ensemble", "trajectory", "time"]))
                    else:
                        df = ds.to_dask_dataframe(
                            dim_order=["ensemble", "trajectory", "time"])
        else:
            df = xr.open_dataset(direc + "/" + file_ending, decode_times=False)
            df = df.to_dask_dataframe(dim_order=["Output Parameter", "ensemble", "trajectory", "time"])
        # The performance of the following command is not the best.

        # df = xr.open_mfdataset(
        #     direc + "/" + file_ending,
        #     parallel=True,
        #     combine="by_coords",
        #     decode_times=False).to_dask_dataframe(dim_order=["Output Parameter", "ensemble", "trajectory", "time"])
    else:
        try:
            df = pd.read_csv(direc + "/*diff*", assume_missing=True) # , blocksize="8GB"
        except:
            try:
                # Try reading netcdf files in bulk
                # Works only if no multiindex is given
                df = xr.open_mfdataset(direc + "/" + file_ending, parallel=True).to_dask_dataframe()
            except:
                # Last fallback which works for most netcdf files
                df = None
                file_list = []
                for f in os.listdir(direc):
                    file_list.append(os.path.join(direc, f))
                file_list = np.sort(np.asarray(file_list))
                for f in file_list:
                    if df is None:
                        df = xr.open_dataset(f).to_dataframe()
                    else:
                        df = df.append(xr.open_dataset(f).to_dataframe())#[columns])
                # Make id and time columns instead of MultiIndex
                df.reset_index(inplace=True)
    return df


def rotate_df(df, pollon, pollat, lon="LONGITUDE", lat="LATITUDE"):
    """
    Rotate the longitude and latitude with the given pole coordinates.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with longitude and latitude
    pollon : float
        Longitude of pole.
    pollat : float
        Latitude of pole.
    lon : String
        "LONGITUDE" for derivative dataframe, "lon" for netCDF dataframe.
    lat : String
        "LATITUDE" for derivative dataframe, "lat" for netCDF dataframe.
    """
    # lat_v, lon_v = rotate_pole(
    #                        np.asarray(df[lon].tolist()),
    #                        np.asarray(df[lat].tolist()),
    #                        pole_lon=pollon,
    #                        pole_lat=pollat)
    # df[lon] = lon_v
    # df[lat] = lat_v
    print("Not implemented")


def norm_deriv(df):
    """
    Given a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv, ratio_deriv, MAP
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param.

    Normalize the ratio of the derivatives for every timestep and add that
    as another column "norm_deriv". The column "MAP" is used to normalize only
        within an area (consecutive timesteps with equal flag).

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns trajectory, timestep, out_param, in_param,
        deriv, ratio_deriv and MAP. On out: Holds another
        column "norm_deriv"
    """
    print("TODO")


def ratio_deriv(df, out_param):
    """
    Given a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param.

    Calculate the ratio of the derivatives for every timestep and add that
    as another column "ratio_deriv".

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns trajectory, timestep, out_param, in_param,
        deriv and optionally MAP.
    out_param : String
        Output parameter to calculate the ratio for

    Returns
    -------
    pandas.Dataframe
        Dataframe with columns trajectory, timestep, out_param, in_param,
        deriv and optionally MAP. Also additional column ratio_deriv
    """
    df_out = df[df.out_param == out_param]
    if df_out.empty:
        print("No such output parameter: {}".format(out_param))
        return None

    denominator = np.abs(df["deriv"]).max()
    df_out["ratio_deriv"] = df_out["deriv"]/denominator
    return df_out


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print("Please give a path with prefix of data and a suffix and an identifier for the csv.")
        print("Falling back to default.")
        prefix = "/data/project/wcb/sb_ice_wcb272280_traj"
        suffix = "start_over_20160922_00"
        ident  = "filt_zero_30"
    else:
        prefix = sys.argv[1]
        suffix = sys.argv[2]
        ident  = sys.argv[3]

    try:
        refs = np.loadtxt('reference_values.txt')
        Tref = refs[0]
        pref = refs[1]
        qref = refs[2]
        Nref = refs[3]
        wref = refs[4]
        tref = refs[5]
    except:
        print("No file with reference values found. Using default values.")
        Tref = 273.15
        pref = 100000
        qref = 1e-06
        Nref = 1
        wref = 1
        tref = 1
        zref = 1

    filt = True
    lo = 1
    hi = 3
    EPSILON = 1e-30
    high = None
    res = load_mult_derivates_big(prefix, suffix,
                                  filt=filt, lo=lo, hi=hi,
                                  EPSILON=EPSILON, high=high)
    print("Trying to store\n")
    saveword = "sb_ice_wcb272280" + "_" + ident + "_" + suffix + ".csv"
    res.to_csv(saveword)

