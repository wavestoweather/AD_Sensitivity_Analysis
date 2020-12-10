import numpy as np
import xarray as xr
from timeit import default_timer as timer
import sys
from glob import glob
import progressbar as pb


input_path = sys.argv[1]
output_path = sys.argv[2]
files = "*.nc_wcb"

dims = ["Output Parameter", "ensemble", "trajectory", "time"]
columns = ["QI", "NCICE", "QC", "NCCLOUD", "QR", "NCRAIN",
            "QS", "NCSNOW", "QH", "NCHAIL", "QG", "NCGRAUPEL",
            "QV", "T", "pressure", "w", "z", "S", "da_1",
            "type", "step", "time_after_ascent"]
t = timer()
ds = None
file_list = sorted(glob(input_path + files))
for i in pb.progressbar(range(len(file_list))):
    f = file_list[i]
    ds_tmp = xr.open_dataset(f, decode_times=False)#.reset_index("trajectory")
    ds_tmp = ds_tmp[columns]
    ds_tmp = ds_tmp.loc[{"Output Parameter": "QV"}]
    ds_tmp["tmp"] = ("trajectory", [i])
    ds_tmp["trajectory_history"] = ("trajectory", ds_tmp["trajectory"])
    ds_tmp = ds_tmp.set_index(trajectory="tmp")
    ds_tmp["tmp"] = ("ensemble", [0])
    ds_tmp["ensemble_history"] = ("trajectory", ds_tmp["ensemble"])
    ds_tmp = ds_tmp.set_index(ensemble="tmp")
    if ds is not None:
        ds = ds.merge(ds_tmp, join="outer")
    else:
        ds = ds_tmp

ds["trajectory"].attrs = {
    "standard_name": "trajectory",
    "long_name": "trajectory"}
ds["trajectory_history"].attrs = {
    "standard_name": "trajectory_history",
    "long_name": "trajectory history",
    "description": "last number is the id of the instance that ran the simulation and the numbers before are the history",
    "auxiliary_data": "yes"}
ds["ensemble"].attrs = {
    "standard_name": "ensemble",
    "long_name": "ensemble"}
ds["ensemble_history"].attrs = {
    "standard_name": "ensemble_history",
    "long_name": "ensemble history",
    "description": "id that is only consistent within a given history via trajectory_history id",
    "auxiliary_data": "yes"}
# ds["type"] = np.unique(ds["type"])[-1]
t = timer()
comp = dict(zlib=True, complevel=9)
encoding = {var: comp for var in ds.data_vars}
ds.to_netcdf(
    path=output_path,
    encoding=encoding,
    compute=True,
    engine="netcdf4",
    format="NETCDF4",
    mode="w")
t2 = timer()
print(f"Writing done in {t2-t} s")