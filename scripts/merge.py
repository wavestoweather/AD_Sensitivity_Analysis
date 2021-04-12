import numpy as np
import xarray as xr
from timeit import default_timer as timer
import sys
from glob import glob
import progressbar as pb


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='''
        Merge multiple outputs from simulations into one NetCDF-file.
        Those simulations are typically single trajectories per file
        with all their sensitivities. Stores the merged dataset with zlib
        compression of level 9.
        ''')
    parser.add_argument('--input_path', required=True,
        help='''
        Path to files ending with ".nc_wcb".
        ''')
    parser.add_argument('--output_path', required=True,
        help='''
        Path with name to store the merged dataset.
        ''')
    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    files = "*.nc_wcb"

    dims = ["Output Parameter", "ensemble", "trajectory", "time"]
    columns = ["QI", "NCICE", "QC", "NCCLOUD", "QR", "NCRAIN",
                "QS", "NCSNOW", "QH", "NCHAIL", "QG", "NCGRAUPEL",
                "QV", "T", "pressure", "w", "z", "S", "da_1",
                "type", "step", "time_after_ascent",
                "QI_OUT", "NI_OUT","QR_OUT", "NR_OUT",
                "QS_OUT", "NS_OUT", "QH_OUT", "NH_OUT",
                "QG_OUT", "NG_OUT", "Inactive", "deposition", "sublimination"]
    t = timer()
    file_list = sorted(glob(input_path + files))
    sets = []
    time_index = None
    min_time = None
    max_time = None

    min_time2 = None
    max_time2 = None
    t = None

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
        min_time_tmp = ds_tmp.coords["time"].min()
        max_time_tmp = ds_tmp.coords["time"].max()

        if time_index is None or (min_time_tmp <= min_time and max_time_tmp > max_time):
            min_time = min_time_tmp
            max_time = max_time_tmp
            time_index = ds_tmp.coords["time"]
            t = np.unique(ds_tmp["type"])[-1]
            min_time2 = min_time
            max_time2 = max_time
        if min_time_tmp < min_time2:
            min_time2 = min_time_tmp
        if max_time_tmp > max_time2:
            max_time2 = max_time_tmp
        sets.append(ds_tmp)

    time_index = np.arange(min_time2, max_time2+19, 20)
    for i in pb.progressbar(range(len(sets))):
        sets[i] = sets[i].reindex({"time": time_index})

    ds = xr.concat(sets, dim="trajectory")

    ds["type"] = t
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