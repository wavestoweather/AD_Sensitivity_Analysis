from glob import glob
from multiprocessing import Pool
import numpy as np
import pandas as pd
import xarray as xr
from progressbar import progressbar as pb
from scripts.Deriv_dask import Deriv_dask
from scripts.latexify import in_params_dic, physical_params
import sys
# python plot_mse.py [error_kind] [base_path] [n_trajs] [df_name] [cum_type]
# ie python plot_mse.py mse  /data/project/wcb/netcdf/perturbed_ensembles/conv_400_0_traj_t000000_p001/ 27 conv_400.csv max

error_kind = sys.argv[1]
if error_kind == "mse":
    error_key = "MSE"
elif error_kind == "maxse":
    error_key = "Max Error"
elif error_kind == "sum":
    error_key = "Cumulative Squared Error"
elif error_kind == "nozeromse":
    error_key = "MSE (no zero)"
elif error_kind == "me":
    error_key = "Mean Error"
elif error_kind == "mae":
    error_key = "Mean Absolute Error"

out_params = ["QV", "QC", "QR", "QG", "QH", "QI", "QS"]
in_params = []
for key in in_params_dic:
    in_params.extend(in_params_dic[key])
for e in physical_params:
    in_params.remove(e)

datashade = False
alpha = 0.5
s = 13
f_limits = (-2,2)
store_path = "pics/"

confidence = 0.90
hist = True
kind = "grid_plot"
ratio_type = ["per_timestep", "window"]

base_path = "/data/project/wcb/netcdf/perturbed_ensembles/conv_600_0_traj_t000000_p001/"
base_path = sys.argv[2]
n_trajs = int(sys.argv[3])
df_name = sys.argv[4]
mse_df = None
cum_type = sys.argv[5]
if cum_type == "max":
    cum_f = np.max
elif cum_type == "mean":
    cum_f = np.mean
elif cum_type == "min":
    cum_f = np.min

def load_mse(i):
    try:
        mean_file = "traj" + str(i) + "_notPerturbed.nc_wcb"
        others_path = base_path + "traj" + str(i) + "/"
        mean_traj = Deriv_dask(
            direc=base_path,
            parquet=False,
            netcdf=True,
            columns=None,
            backend="matplotlib",
            file_ending=mean_file)
        mean_traj.cache_data(
            in_params=in_params,
            out_params=out_params,
            x_axis="time_after_ascent",
            compute=True)
        df = mean_traj.get_mse(
            out_params=out_params,
            others_path=others_path,
            in_params=in_params,
            ratio_type=ratio_type,
            kind=error_kind)
        return df
    except:
        return None

try:
    mse_df = pd.read_csv("stats_full/" + error_kind + "_" + df_name)
except:
    with Pool(6) as p:
        mse_df = pd.concat(p.map(load_mse, range(n_trajs)))
    mse_df.to_csv("stats_full/" + error_kind + "_" + df_name)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("Finished parsing the data")
exit()
output_params = np.unique(mse_df["Output Parameter"])
ratio_types = np.unique(mse_df["Ratio Type"])
model_params = np.unique(mse_df["Perturbed Parameter"])
df_dic = {
    "MSE": [],
    "Sensitivity": [],
    "Output Parameter": [],
    "Perturbed Parameter": [],
    "Ratio Type": []
}
try:
    mean_df.read_csv(error_kind + "_" + cum_type + "_" + df_name)
except:
    for rt in ratio_types:
        for op in output_params:
            for i in pb(range(len(model_params))):
                mp = model_params[i]
                tmp_df = mse_df.loc[(mse_df["Ratio Type"] == rt) & (mse_df["Output Parameter"] == op) & (mse_df["Perturbed Parameter"] == mp)]
                df_dic["MSE"].append(cum_f(tmp_df["MSE"]))
                df_dic["Sensitivity"].append(cum_f(tmp_df["Sensitivity"]))
                df_dic["Output Parameter"].append(op)
                df_dic["Perturbed Parameter"].append(mp)
                df_dic["Ratio Type"].append(rt)
    mean_df = pd.DataFrame.from_dict(df_dic)
    mean_df.to_csv(error_kind + "_" + cum_type + "_" + df_name)
exit()
mean_traj = Deriv_dask(
    direc=base_path,
    parquet=False,
    netcdf=True,
    columns=None,
    backend="matplotlib",
    file_ending="traj0_notPerturbed.nc_wcb")
mean_traj.cache_data(
    in_params=in_params,
    out_params=out_params,
    x_axis="time_after_ascent",
    compute=True)

for i in [0, 1, 2, 3]:
    if i == 0:
        tmp_df = mse_df
    if i == 1:
        tmp_df = mse_df.loc[mse_df["Sensitivity"] != 0]
    if i == 2:
        mse_df = mean_df
    if i == 3:
        tmp_df = mean_df.loc[mean_df["Sensitivity"] != 0]

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=hist,
        confidence=None,
        kind=kind,
        error_key=error_key)

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=hist,
        confidence=confidence,
        split_ellipse=True,
        kind=kind,
        error_key=error_key)

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=hist,
        log_x=False,
        log_y=True,
        kind=kind,
        error_key=error_key)

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=False,
        confidence=confidence,
        abs_x=False,
        log_x=True,
        log_y=True,
        kind=kind,
        error_key=error_key)

    mean_traj.plot_mse(
        out_params=out_params,
        mse_df_=tmp_df,
        in_params=in_params,
        datashade=datashade,
        alpha=alpha,
        formatter_limits=f_limits,
        s=s,
        hist=False,
        confidence=confidence,
        abs_x=False,
        log_x=True,
        log_y=False,
        kind=kind,
        error_key=error_key)