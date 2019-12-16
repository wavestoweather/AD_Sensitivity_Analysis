import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from progressbar import progressbar as pb
import sys
import os
import xarray as xr


refs = np.loadtxt('reference_values.txt')
Tref = refs[0]
pref = refs[1]
qref = refs[2]
Nref = refs[3]
wref = refs[4]
tref = refs[5]

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
params_dict2 = {"_diff_0.txt": "p", "_diff_1.txt": "T",
               "_diff_2.txt": "w", "_diff_3.txt": "S", "_diff_4.txt": "qc",
               "_diff_5.txt": "qr", "_diff_6.txt": "qv", "_diff_7.txt": "Nc",
               "_diff_8.txt": "Nr", "_diff_9.txt": "Nv",
               "_diff_10.txt": "qi", "_diff_11.txt": "Ni",
               "_diff_12.txt": "vi", "_diff_13.txt": "qs",
               "_diff_14.txt": "Ns", "_diff_15.txt": "qg",
               "_diff_16.txt": "Ng", "_diff_17.txt": "qh",
               "_diff_18.txt": "Nh", "_diff_19.txt": "qiout",
               "_diff_20.txt": "qsout", "_diff_21.txt": "qrout",
               "_diff_22.txt": "qgout",
               "_diff_23.txt": "qhout", "_diff_24.txt": "latent_heat",
               "_diff_25.txt": "latent_cool"}

def filter_zeros(df_dict, EPSILON=1e-31):
    """
    Drop all columns that have zero inpact

    Parameters
    ----------
    df_dict     A dictionary of pandas.Dataframe with key the output parameter
                and the columns the input parameters and values are the
                derivatives.
    EPSILON     If all values in a column are lower than EPSILON,
                drop that one.

    Return
    ------
    Modified df_dict with removed columns in each dataframe where only (near)
    zeros existed before.
    """
    key_drop = []
    for key in df_dict:
        to_drop = []
        for column in df_dict[key]:
            if column == "timestep" or column == "trajectory":
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

                print("Dropping {} entirely.".format(key))
                key_drop.append(key)

    for key in key_drop:
        del df_dict[key]
    return df_dict


def filter_high(df_dict, high=1e1):
    """
    Drop all columns that have a suspiciously high impact.

    Parameters
    ----------
    df_dict     A dictionary of pandas.Dataframe with key the output parameter
                and the columns the input parameters and values are the
                derivatives.
    high        If some values in a column are higher than high,
                drop that one.

    Return
    ------
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


def load_derivatives(prefix="", suffix="", filt=False, EPSILON=1e-31):
    """
    Load several dataframes and store them in a dictionary where
    the key is the output parameter of a model and the columns are the
    input parameters and the values are the derivatives.

    Parameters
    ----------
    prefix      Prefix of the datafile such as "sb_ice".
    suffix      Suffix of the datafile such as "start_over".

    Return
    ------
    Dictionary with dataframes as described above.

    """
    df_dict = {}
    for key in params_dict.keys():
        print("Loading {}".format(key))
        tmp = pd.read_csv(
            "data/" + prefix + suffix + params_dict[key], sep=",")

        df_dict[key] = tmp
    if filt:
        df_dict = filter_zeros(df_dict, EPSILON)
    for key in df_dict.keys():

        if key == "qr":
        #     print("Renaming stuff in")
        #     print(df_dict[key].keys())
            df_dict[key] = df_dict[key].rename(mapper={"dinv_z": r"$z_{inv}$",
                                 "drain_a_geo": r"$a_{rain, geo}$",
                                 "drain_b_geo": r"$b_{rain, geo}$",
                                 "drain_alpha": r"$\alpha_{rain}$",
                                 "drain_nu": r"$\nu_{rain}$"}, axis="columns")
                # df_dict[key] = df_dict[key].rename({"dinv_z": "inv_z",
                #                  "drain_a_geo": "rain_a_geo",
                #                  "drain_b_geo": "rain_b_geo",
                #                  "drain_alpha": "rain_alpha",
                #                  "drain_nu": "rain_nu"})
            print(df_dict[key])
        df_dict[key] = transform_df(df_dict[key])
    return df_dict


def load_nc(inp="/mnt/localscratch/data/project/m2_jgu-tapt/o"
                + "nline_trajectories/foehn201305_case/" +
                "foehn201305_warming1.nc"):
    ds = xr.open_dataset(inp)
    df = ds.to_dataframe()
    return df


def load_output(prefix="sb_ice", suffix="_start_over"):
    data = pd.read_csv("data/" + prefix + suffix + ".txt", sep=",")
    data["p"] = data["p"]*pref/100  # We want hPa
    data["T"] = data["T"]*Tref
    data["w"] = data["w"]*wref
    data["qc"] = data["qc"]*qref
    data["qr"] = data["qr"]*qref
    data["qs"] = data["qs"]*qref
    data["qg"] = data["qg"]*qref
    data["qh"] = data["qh"]*qref
    data["qi"] = data["qi"]*qref
    data["qv"] = data["qv"]*qref
    data["qiout"] = data["qiout"]*qref
    data["qsout"] = data["qsout"]*qref
    data["qrout"] = data["qrout"]*qref
    data["qgout"] = data["qgout"]*qref
    data["qhout"] = data["qhout"]*qref
    data["latent_heat"] = data["latent_heat"]*Tref
    data["latent_cool"] = data["latent_cool"]*Tref
    return data


def transform_df(df):
    """
    Create a new pandas.DataFrame with column "param", "timestep", "deriv"
    that can be used for plotting with seaborn.lineplot.

    Parameters
    ----------
    df      pandas.DataFrame where columns are the names of the parameters.

    Returns
    -------
    Transformed Dataframe.
    """
    dicti = {"param": [], "timestep": [], "deriv": []}
    for key in df:
        if key == "timestep" or key == "trajectory":
            continue
        dicti["timestep"].extend(df["timestep"].tolist())
        dicti["deriv"].extend(df[key].tolist())
        dicti["param"].extend([key for i in range(len(df["timestep"]))])
    return pd.DataFrame(dicti)


def load_mult_derivates(prefix="", suffix="", filt=False, EPSILON=1e-31,
                        lo=0, hi=0, high=None):
    """
    Create a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param

    Parameters
    ----------
    prefix      Prefix of the datafile such as "sb_ice".
    suffix      Suffix of the datafile such as "start_over".

    Return
    ------
    Pandas dataframe as described above.
    """
    df_dict = {"timestep": [], "trajectory": [], "out_param": [],
               "in_param": [], "deriv": []}

    for i in pb(range(lo, hi+1), redirect_stdout=True):
        tmp_dict = {}
        try:
            for out_param in params_dict.keys():
                tmp = pd.read_csv(prefix + str(i) + "_"
                                + suffix + params_dict[out_param], sep=",")
                tmp_dict[out_param] = tmp
            if filt:
                tmp_dict = filter_zeros(tmp_dict, EPSILON)
            if not high is None:
                tmp_dict = filter_high(tmp_dict, high)
            for out_param in tmp_dict.keys():
                tmp_dict[out_param] = transform_df(tmp_dict[out_param])
            for out_param in tmp_dict.keys():
                n_entries = len(tmp_dict[out_param].index)
                dic = tmp_dict[out_param]
                df_dict["trajectory"].extend([i for j in range(n_entries)])
                df_dict["out_param"].extend([out_param for j in range(n_entries)])
                df_dict["timestep"].extend(dic["timestep"])
                df_dict["in_param"].extend(dic["param"])
                df_dict["deriv"].extend(dic["deriv"])

        except:
            pass
    return pd.DataFrame(df_dict)

def load_mult_derivates_big(prefix="", suffix="", filt=False, EPSILON=1e-31,
                        lo=0, hi=0, high=None):
    """
    Create a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param

    Parameters
    ----------
    prefix      Prefix of the datafile such as "sb_ice".
    suffix      Suffix of the datafile such as "start_over".

    Return
    ------
    Pandas dataframe as described above.
    """
    df = pd.DataFrame(data={"timestep": [], "trajectory": [], "out_param": [],
                            "in_param": [], "deriv": []})
    for i in pb(range(lo, hi+1), redirect_stdout=True):
        tmp_dict = {}
        try:
            for out_param in params_dict.keys():
                tmp = pd.read_csv(prefix + str(i) + "_"
                                + suffix + params_dict[out_param], sep=",")
                tmp_dict[out_param] = tmp
            if filt:
                tmp_dict = filter_zeros(tmp_dict, EPSILON)
            if not high is None:
                tmp_dict = filter_high(tmp_dict, high)
            for out_param in tmp_dict.keys():
                tmp_dict[out_param] = transform_df(tmp_dict[out_param])
            for out_param in tmp_dict.keys():
                n_entries = len(tmp_dict[out_param].index)
                # df_tmp = pd.DataFrame(
                #     data={"timestep": tmp_dict[out_param]["timestep"],
                #           "trajectory": [i for j in range(n_entries)],
                #           "out_param": [out_param for j in range(n_entries)],
                #           "in_param": tmp_dict[out_param]["param"],
                #           "deriv": tmp_dict[out_param]["deriv"]})
                df = df.append(pd.DataFrame(
                    data={"timestep": tmp_dict[out_param]["timestep"],
                          "trajectory": [i for j in range(n_entries)],
                          "out_param": [out_param for j in range(n_entries)],
                          "in_param": tmp_dict[out_param]["param"],
                          "deriv": tmp_dict[out_param]["deriv"]}), ignore_index=True)

        except:
            pass
    return df


def load_mult_derivates_directory(direc="", filt=True, 
                                  EPSILON=0.0, trajectories=[1], suffix="20160922_00"):
    """
    Create a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param

    Parameters
    ----------
    prefix      Prefix of the datafile such as "sb_ice".
    suffix      Suffix of the datafile such as "start_over".

    Return
    ------
    Pandas dataframe as described above.
    """
    df = pd.DataFrame(data={"timestep": [], "trajectory": [], "out_param": [],
                            "in_param": [], "deriv": []})
    file_list = [f for f in listdir(direc= if isfile(join(direc, f))]
    file_list2 = []
    for f in file_list:
        s = f.split("traj")
        s = s[1].split("_")
        if int(s[0]) in trajectories:
            file_list2.append(f)
    for f in pb(file_list2, redirect_stdout=True):       
        tmp_dict = {}
        try:
            tmp = pd.read_csv(f, sep=",")
            out_param = params_dict2[f.split(suffix)[1]]
            tmp_dict[out_param] = tmp
            if filt:
                tmp_dict = filter_zeros(tmp_dict, EPSILON)
            for out_param in tmp_dict.keys():
                tmp_dict[out_param] = transform_df(tmp_dict[out_param])
            for out_param in tmp_dict.keys():
                n_entries = len(tmp_dict[out_param].index)
                df = df.append(pd.DataFrame(
                    data={"timestep": tmp_dict[out_param]["timestep"],
                          "trajectory": [i for j in range(n_entries)],
                          "out_param": [out_param for j in range(n_entries)],
                          "in_param": tmp_dict[out_param]["param"],
                          "deriv": tmp_dict[out_param]["deriv"]}), ignore_index=True)
        except:
            pass
    return df
            
    
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

