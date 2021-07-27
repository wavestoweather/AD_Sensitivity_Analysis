from glob import glob
import numpy as np
from timeit import default_timer as timer
import xarray as xr


def load_input_datasets(data_path, smooth_coords=True, verbosity=0):
    file_list = sorted(glob(data_path + "*.nc_wcb"))

    delta_coord = 0.000

    traj_index = 0
    sets = []
    time_index = None
    min_time = None
    max_time = None

    min_time2 = None
    max_time2 = None
    for f in file_list:
        if args.verbosity > 1:
            t = timer()
        ds_tmp = xr.open_dataset(f, decode_times=False)
        n_trajs = len(np.unique(ds_tmp["trajectory"]))
        ds_tmp["tmp"] = ("trajectory", np.arange(traj_index, traj_index + n_trajs))
        ds_tmp = ds_tmp.set_index(trajectory="tmp")
        if smooth_coords:

            def smooth_them(coords):
                for i in range(len(coords) - 1):
                    if coords[i + 1] - coords[i] > delta_coord:
                        coords[i + 1] = coords[i] + delta_coord
                    elif coords[i] - coords[i + 1] > delta_coord:
                        coords[i + 1] = coords[i] - delta_coord
                return coords

            ds_tmp["lat"] = smooth_them(ds_tmp["lat"])
            ds_tmp["lon"] = smooth_them(ds_tmp["lon"])

        min_time_tmp = ds_tmp.coords["time"].min()
        max_time_tmp = ds_tmp.coords["time"].max()
        traj_index += n_trajs
        if time_index is None or (min_time_tmp <= min_time and max_time_tmp > max_time):
            min_time = min_time_tmp
            max_time = max_time_tmp
            time_index = ds_tmp.coords["time"]
            min_time2 = min_time
            max_time2 = max_time
        if min_time_tmp < min_time2:
            min_time2 = min_time_tmp
        if max_time_tmp > max_time2:
            max_time2 = max_time_tmp
        sets.append(ds_tmp)
        if args.verbosity > 1:
            t2 = timer()
            print(f"Loading {f} done in {t2-t} s")

    time_index = np.arange(min_time2, max_time2 + 19, 20)
    for i in range(len(sets)):
        sets[i] = sets[i].reindex({"time": time_index})

    if args.verbosity > 1:
        t = timer()
    ds = xr.concat(sets, dim="trajectory", join="inner")
    if args.verbosity > 1:
        t2 = timer()
        print(f"Concatenating all files done in {t2-t} s")

    ds.attrs["duration_in_sec"] = (
        ds.coords["time"].max() - ds.coords["time"].min()
    ).item()
    ds["trajectory"].attrs = {"standard_name": "trajectory", "long_name": "trajectory"}
    ds["ensemble"].attrs = {"standard_name": "ensemble", "long_name": "ensemble"}
    if "Output Parameter" in ds:
        ds = ds.rename({"Output Parameter": "Output_Parameter"})

    for col in ds:
        if col != "pressure" and col != "lon" and col != "lat":
            if "auxiliary_data" not in ds[col].attrs:
                ds[col].attrs["auxiliary_data"] = "yes"
    return ds


def load_perturbed_datasets(data_path, verbosity=0):
    file_list = sorted(glob(data_path + "*.nc_wcb"))

    columns = [
        "QI",
        "NCICE",
        "QC",
        "NCCLOUD",
        "QR",
        "NCRAIN",
        "QS",
        "NCSNOW",
        "QH",
        "NCHAIL",
        "QG",
        "NCGRAUPEL",
        "QV",
        "T",
        "pressure",
        "w",
        "z",
        "S",
        "da_1",
        "type",
        "step",
        "time_after_ascent",
        "QI_OUT",
        "NI_OUT",
        "QR_OUT",
        "NR_OUT",
        "QS_OUT",
        "NS_OUT",
        "QH_OUT",
        "NH_OUT",
        "QG_OUT",
        "NG_OUT",
        "Inactive",
        "deposition",
        "sublimination",
    ]

    t = None
    time_index = None
    min_time = None
    max_time = None

    min_time2 = None
    max_time2 = None
    sets = []
    for i in range(len(file_list)):
        if args.verbosity > 1:
            t = timer()
        f = file_list[i]
        ds_tmp = xr.open_dataset(f, decode_times=False)
        # We cut off the sensitivities and take solely the results
        ds_tmp = ds_tmp[columns]
        ds_tmp = ds_tmp.loc[{"Output Parameter": "QV"}]

        ds_tmp["tmp"] = ("trajectory", [i])
        ds_tmp = ds_tmp.set_index(trajectory="tmp")
        ds_tmp["tmp"] = ("ensemble", [0])
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
        if args.verbosity > 1:
            t2 = timer()
            print(f"Loading {f} done in {t2-t} s")

    time_index = np.arange(min_time2, max_time2 + 19, 20)
    for i in range(len(sets)):
        sets[i] = sets[i].reindex({"time": time_index})
    if args.verbosity > 1:
        t = timer()
    ds = xr.concat(sets, dim="trajectory")
    if args.verbosity > 1:
        t2 = timer()
        print(f"Concatenating all files done in {t2-t} s")
    ds["type"] = t
    ds["trajectory"].attrs = {"standard_name": "trajectory", "long_name": "trajectory"}
    ds["ensemble"].attrs = {"standard_name": "ensemble", "long_name": "ensemble"}

    if "Output Parameter" in ds:
        ds = ds.rename({"Output Parameter": "Output_Parameter"})

    for col in ds:
        if col != "pressure" and col != "lon" and col != "lat":
            if "auxiliary_data" not in ds[col].attrs:
                ds[col].attrs["auxiliary_data"] = "yes"
    return ds


def load_ensemble_datasets(data_path, verbosity=0):
    file_list = sorted(glob(data_path + "*.nc_wcb"))
    file_list.remove(data_path + "_notPerturbed.nc_wcb")

    def load(f):
        perturbed_param = f.split("/")[-1].split(".nc")[0]
        # print(f"Loading {perturbed_param}")
        if perturbed_param == "_notPerturbed":
            return None
        ds_tmp = xr.open_dataset(f, decode_times=False)
        ds_tmp["Perturbed_Parameter"] = perturbed_param
        return ds_tmp

    ds_list = []
    ds_list2 = []
    step = 10
    for i in range(0, len(file_list), step):
        if args.verbosity > 1:
            t = timer()
        ds_list.append(
            xr.concat(
                [load(f) for f in file_list[i : i + step]], dim="Perturbed_Parameter"
            )
        )
        if args.verbosity > 1:
            t2 = timer()
            print(f"Loading files with index {i} to {i+step} done in {t2-t} s")
        if i % (step * step) == 0 and i > 0:
            if args.verbosity > 1:
                t = timer()
            ds_list2.append(xr.concat(ds_list, dim="Perturbed_Parameter"))
            ds_list.clear()
            if args.verbosity > 1:
                t2 = timer()
                print(f"Concatenating {step} files done in {t2-t} s")
    if len(ds_list) > 0:
        if args.verbosity > 1:
            t = timer()
        ds_list2.append(xr.concat(ds_list, dim="Perturbed_Parameter"))
        ds_list.clear()
        if args.verbosity > 1:
            t2 = timer()
            print(f"Concatenating {step} files done in {t2-t} s")
    if args.verbosity > 1:
        t = timer()
    ds = xr.concat(ds_list2, dim="Perturbed_Parameter")
    ds_list2.clear()
    if args.verbosity > 1:
        t2 = timer()
        print(f"Concatenating all files done in {t2-t} s")

    ds["Perturbed_Parameter"].attrs = {
        "standad_name": "Perturbed_Parameter",
        "long_name": "Perturbed Parameter",
        "description": "This Parameter had been perturbed for this ensemble",
        "auxiliary_data": "yes",
    }
    if "Output Parameter" in ds:
        ds = ds.rename({"Output Parameter": "Output_Parameter"})
    for col in ds:
        if col != "pressure" and col != "lon" and col != "lat":
            if "auxiliary_data" not in ds[col].attrs:
                ds[col].attrs["auxiliary_data"] = "yes"

    return ds


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Merge output of an ensemble to fewer NetCDF-files.
        """
    )
    parser.add_argument(
        "--data_path",
        default="/data/project/wcb/netcdf/perturbed_ensembles/",
        help="""
        Path to folders with ensemble datasets.
        """,
    )
    parser.add_argument(
        "--store_path",
        default="/data/project/wcb/netcdf/perturbed_ensembles/",
        help="""
        Path where to store the merged files.
        """,
    )
    parser.add_argument(
        "--merge_multiple_ensembles",
        action="store_true",
        help="""
        if true: Merge multiple ensembles, each in a single NetCDF-file
        of the form "d[param_name].nc_wcb" to a single file.
        Otherwise merge outputs of a simulation to a single NetCDF-file of
        the form "d[param_name].nc_wcb"
        """,
    )
    parser.add_argument(
        "--merge_input_files",
        action="store_true",
        help="""
        Merge multiple input NetCDF-files that are not results of our own
        simulation. This is mainly used to merge the statistic datasets.
        """,
    )
    parser.add_argument(
        "--complevel",
        default=9,
        type=int,
        help="""
        Compression level for zlib when storing files. Set to 0 if something
        is off with other programs, i.e. Met3D can have problems with compression.
        """,
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help="""
        Set verbosity.
        0: No output
        1: Get coarse time information
        2: Get detailed time information
        """,
    )
    args = parser.parse_args()

    if args.merge_multiple_ensembles:
        if args.verbosity > 0:
            t = timer()
        ds = load_ensemble_datasets(args.data_path, verbosity=args.verbosity)
        if args.verbosity > 0:
            t2 = timer()
            print(f"Loading done in {t2-t} s")

        if args.verbosity > 0:
            t = timer()
        store_path = args.store_path
        if store_path[-1] == "/":
            store_path = store_path[0:-1] + ".nc_wcb"
        if store_path[-7::] != ".nc_wcb":
            store_path += ".nc_wcb"
        comp = dict(zlib=True, complevel=args.complevel)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(
            path=store_path,
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        if args.verbosity > 0:
            t2 = timer()
            print(f"Storing done in {t2-t} s")
    elif args.merge_input_files:
        if args.verbosity > 0:
            t = timer()
        ds = load_input_datasets(args.data_path, verbosity=args.verbosity)
        if args.verbosity > 0:
            t2 = timer()
            print(f"Loading done in {t2-t} s")

        if args.verbosity > 0:
            t = timer()
        store_path = args.store_path
        if store_path[-1] == "/":
            store_path = store_path[0:-1] + ".nc_wcb"
        if store_path[-7::] != ".nc_wcb":
            store_path += ".nc_wcb"
        comp = dict(zlib=True, complevel=args.complevel)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(
            path=store_path,
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        if args.verbosity > 0:
            t2 = timer()
            print(f"Storing done in {t2-t} s")
    else:
        if args.verbosity > 0:
            t = timer()
        ds = load_perturbed_datasets(args.data_path, verbosity=args.verbosity)
        if args.verbosity > 0:
            t2 = timer()
            print(f"Loading done in {t2-t} s")

        if args.verbosity > 0:
            t = timer()
        store_path = args.store_path
        comp = dict(zlib=True, complevel=args.complevel)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(
            path=store_path,
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
        if args.verbosity > 0:
            t2 = timer()
            print(f"Storing done in {t2-t} s")
