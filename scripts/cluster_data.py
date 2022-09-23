import os
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from tqdm.auto import tqdm
import xarray as xr


def load_data(file_path, x, only_asc600=False, inoutflow_time=-1):
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    param_ids = None
    list_of_arrays = []
    for f in tqdm(files):
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[
            [x, "asc600", "phase"]
        ]
        if inoutflow_time > 0:
            ds_flow = ds.where(ds["asc600"] == 1)["asc600"]
            ds_flow = ds_flow.rolling(
                time=inoutflow_time * 2,  # once for inflow, another for outflow
                min_periods=1,
                center=True,
            ).reduce(np.nanmax)
            ds = ds.where(ds_flow == 1)
        elif only_asc600:
            ds = ds.where(ds["asc600"] == 1)
        if param_ids is None and x[0] == "d":
            param_ids = ds["Output_parameter_ID"]
        # fill values where
        data = ds[x].dropna("time", how="all")
        data = data.expand_dims({"file": [f]})
        list_of_arrays.append(data)
    data = xr.concat(list_of_arrays, dim="file", join="outer", combine_attrs="override")
    return data


def get_average(data):
    return data.mean(dim="time", skipna=True, keep_attrs=True)


def get_cluster(data, k, tol=1e-4, x=None, param_names=None):
    """

    Parameters
    ----------
    data : xarray.Dataset or xarray.DataArray
        A dataset of trajectories from a sensitivity simulation or one column of the set.
    tol : float
    k : int
    x : string
    out_params : list of strings

    Returns
    -------
    sklearn.cluster.Kmeans of the array or a dictionary of output parameter names and the clusters for each.

    """
    if x is not None and isinstance(data, xr.Dataset):
        data_array = data[x]
        avg_name = "avg " + x
    else:
        data_array = data
        avg_name = "avg " + data.name

    shape = (len(data_array["file"]) * len(data_array["trajectory"]), 1)

    if param_names is not None:
        kmeans_dic = {
            "clusters": np.asarray([]),
            avg_name: np.asarray([]),
            "param": np.asarray([]),
            "trajectory": np.asarray([]),
            "file": np.asarray([]),
        }
        traj_offset = 0
        for i, out_p in enumerate(param_names):
            fit_data = data_array.sel({"Output_Parameter_ID": i}).values
            fit_data = np.reshape(fit_data, shape)
            kmeans_dic[out_p] = [
                fit_data[~np.isnan(fit_data[:, 0])],
                KMeans(n_clusters=k, tol=tol, random_state=42).fit_predict(
                    fit_data[~np.isnan(fit_data[:, 0])]
                ),
            ]
            clusters = KMeans(n_clusters=k, tol=tol, random_state=42).fit_predict(
                fit_data[~np.isnan(fit_data[:, 0])]
            )
            kmeans_dic["param"] = np.append(kmeans_dic["param"], out_p)
            trajectory = data_array["trajectory"]
            for _ in range(len(data_array["file"]) - 1):
                trajectory = np.append(trajectory, data_array["trajectory"])
            filenames = np.repeat(
                data_array["file"].values, len(data_array["trajectory"])
            )
            kmeans_dic["file"] = np.append(
                kmeans_dic["file"], filenames[~np.isnan(fit_data[:, 0])]
            )
            kmeans_dic["trajectory"] = np.append(
                kmeans_dic["trajectory"], trajectory[~np.isnan(fit_data[:, 0])]
            )
            kmeans_dic["clusters"] = np.append(kmeans_dic["clusters"], clusters)
            kmeans_dic[avg_name] = np.append(
                kmeans_dic[avg_name], fit_data[~np.isnan(fit_data[:, 0])].flatten()
            )
        return pd.DataFrame.from_dict(kmeans_dic)
    else:
        fit_data = np.reshape(data_array.values, shape)
        clusters = KMeans(n_clusters=k, tol=tol, random_state=42).fit_predict(
            fit_data[~np.isnan(fit_data[:, 0])]
        )
        trajectory = data_array["trajectory"]
        for _ in range(len(data_array["file"]) - 1):
            trajectory = np.append(trajectory, data_array["trajectory"])
        filenames = np.repeat(data_array["file"].values, len(data_array["trajectory"]))
        return pd.DataFrame.from_dict(
            {
                "cluster": clusters,
                "trajectory": trajectory[~np.isnan(fit_data[:, 0])],
                avg_name: fit_data[~np.isnan(fit_data[:, 0])].flatten(),
                "file": filenames[~np.isnan(fit_data[:, 0])],
            }
        )


def plot_cluster_data(
    data,
    x,
    y=None,
    plot_type="histplot",
    width=16,
    height=12,
    log_scale=(False, False),
    palette="tab10",
    **kwargs
):

    sns.set(rc={"figure.figsize": (width, height)})

    if plot_type == "histplot":
        return sns.histplot(
            data=data,
            x=x,
            hue="cluster",
            palette=palette,
            log_scale=log_scale,
            **kwargs
        )
    elif plot_type == "scatter":
        g = sns.scatterplot(data=data, x=x, y=y, hue="cluster", palette=palette)
        if log_scale[0]:
            g.set_xscale("log")
        if log_scale[1]:
            g.set_yscale("log")
        return g


def get_traj_near_center(data, x, n=1):
    centers = []
    for cluster in np.unique(data["cluster"]):
        centers.append(np.mean(data.loc[data["cluster"] == cluster][x]))
    traj_centers = []
    for center in centers:
        idx = (np.abs(data[x] - center)).argsort()[:n]
        traj_centers.append(data.iloc[idx])
    return traj_centers


def extract_trajs(centers, file_path, store_path, only_asc600=False, inoutflow_time=-1):
    ds_list = []
    traj_enum = 0
    cluster = 0
    min_time = None
    max_time = None
    for center_df in tqdm(centers):
        for idx, row in tqdm(center_df.iterrows(), leave=False, total=len(center_df)):
            f = row["file"]
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
            tmp_min = ds["time"].min()
            if min_time is None or min_time > tmp_min:
                min_time = tmp_min
            tmp_max = ds["time"].max()
            if max_time is None or max_time < tmp_max:
                max_time = tmp_max

    for center_df in tqdm(centers):
        data_arrays = []
        for idx, row in tqdm(center_df.iterrows(), leave=False, total=len(center_df)):
            traj = row["trajectory"]
            f = row["file"]
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
            ds = ds.sel({"trajectory": [traj]})
            if inoutflow_time > 0:
                ds_flow = ds.where(ds["asc600"] == 1)["asc600"]
                ds_flow = ds_flow.rolling(
                    time=inoutflow_time * 2,  # once for inflow, another for outflow
                    min_periods=1,
                    center=True,
                ).reduce(np.nanmax)
                ds = ds.where(ds_flow == 1)
            elif only_asc600:
                ds = ds.where(ds["asc600"] == 1)
            ds["asc600"] = ds["asc600"].where(ds["asc600"] >= 0)
            traj_attrs = ds["trajectory"].attrs
            ens_attrs = ds["ensemble"].attrs
            time_attrs = ds["time"].attrs
            ds["trajectory"] = [traj_enum]
            ds["cluster"] = xr.DataArray(
                data=[[np.repeat(np.float32(cluster), len(ds["time"].values))]],
                coords=dict(
                    ensemble=[0],
                    trajectory=[traj_enum],
                    time=ds["time"].values,
                ),
                dims=["ensemble", "trajectory", "time"],
                name="cluster",
                attrs={
                    "auxiliary_data": "yes",
                    "_FillValue": np.NaN,
                    "standard_name": "cluster",
                    "long_name": "cluster",
                },
            )
            ds["trajectory"].attrs = traj_attrs
            ds["ensemble"].attrs = ens_attrs
            ds["time"].attrs = time_attrs
            weird_lon = np.count_nonzero(ds["lon"] == 0)
            weird_lat = np.count_nonzero(ds["lat"] == 0)
            # A rare bug can happen where a single lat or lon is set to zero inbetween valid values
            if weird_lon > 0:
                idxs = np.argwhere(ds["lon"].where(ds["lon"] == 0).values == 0)[:, 2]
                # We assume the idxs are consecutive
                n_fills = len(idxs)
                lon_tmp = ds["lon"].isel(
                    {
                        "time": [idxs[0] - 1, idxs[-1] + 1],
                    }
                )
                val_delta = (
                    (lon_tmp.isel({"time": 1}) - lon_tmp.isel({"time": 0})) / n_fills
                ).values.item()
                val_start = lon_tmp.isel({"time": 1}).values.item()
                for i_idx, idx in enumerate(idxs):
                    ds["lon"][0, 0, idx] = val_start + (i_idx + 1) * val_delta
            if weird_lat > 0:
                idxs = np.argwhere(ds["lat"].where(ds["lat"] == 0).values == 0)[:, 2]
                # We assume the idxs are consecutive
                n_fills = len(idxs)
                lat_tmp = ds["lat"].isel(
                    {
                        "time": [idxs[0] - 1, idxs[-1] + 1],
                    }
                )
                val_delta = (
                    (lat_tmp.isel({"time": 1}) - lat_tmp.isel({"time": 0})) / n_fills
                ).values.item()
                val_start = lat_tmp.isel({"time": 1}).values.item()
                for i_idx, idx in enumerate(idxs):
                    ds["lat"][0, 0, idx] = val_start + (i_idx + 1) * val_delta
            # Fill time coordinate and add NaNs where necessary
            tmp_min = ds["time"].min()
            tmp_max = ds["time"].max()
            if tmp_min > min_time:
                variables = {}
                times = int(tmp_min + 20 - min_time) // 30
                sens_shape = (len(ds["Output_Parameter_ID"]), 1, 1, times)
                val_shape = (1, 1, times)
                add_sens = np.reshape(
                    np.repeat(np.NaN, times * len(ds["Output_Parameter_ID"])),
                    sens_shape,
                )
                add_vals = np.reshape(np.repeat(np.NaN, times), val_shape)
                for var in ds:
                    if var[0] != "d" or var == "deposition":
                        variables[var] = (["ensemble", "trajectory", "time"], add_vals)
                    else:
                        variables[var] = (
                            ["Output_Parameter_ID", "ensemble", "trajectory", "time"],
                            add_sens,
                        )
                tmp_set = xr.Dataset(
                    data_vars=variables,
                    coords=dict(
                        ensemble=[0],
                        trajectory=[traj_enum],
                        time=np.arange(min_time, tmp_min, 30),
                    ),
                )
                ds = xr.concat(
                    [ds, tmp_set], join="outer", dim="time", combine_attrs="override"
                )
                # ds["cluster"] = ds["cluster"].where(~np.isnan(ds["pressure"]))
            elif tmp_max < max_time:
                variables = {}
                times = int(max_time + 20 - tmp_max) // 30
                sens_shape = (len(ds["Output_Parameter_ID"]), 1, 1, times)
                val_shape = (1, 1, times)
                add_sens = np.reshape(
                    np.repeat(np.NaN, times * len(ds["Output_Parameter_ID"])),
                    sens_shape,
                )
                add_vals = np.reshape(np.repeat(np.NaN, times), val_shape)
                for var in ds:
                    if var[0] != "d" or var == "deposition":
                        variables[var] = (["ensemble", "trajectory", "time"], add_vals)
                    else:
                        variables[var] = (
                            ["Output_Parameter_ID", "ensemble", "trajectory", "time"],
                            add_sens,
                        )
                        if var == "asc600":
                            variables[var] = variables[var].where(
                                ~np.isnan(variables[var]), other=0
                            )
                tmp_set = xr.Dataset(
                    data_vars=variables,
                    coords=dict(
                        ensemble=[0],
                        trajectory=[traj_enum],
                        time=np.arange(tmp_max + 30, max_time + 20, 30),
                    ),
                )
                ds = xr.concat(
                    [ds, tmp_set], join="outer", dim="time", combine_attrs="override"
                )
                # ds["cluster"] = ds["cluster"].where(~np.isnan(ds["pressure"]))

            if ds["time"].dtype == np.float64:
                ds["time"] = ds["time"].astype(np.float32)
            for col in ds:
                if ds[col].dtype == np.float64:
                    ds[col] = ds[col].astype(np.float32)
            traj_enum += 1
            data_arrays.append(ds)
        ds_list.append(
            xr.concat(
                data_arrays, dim="trajectory", join="outer", combine_attrs="override"
            )
        )
        cluster += 1
    ds = xr.concat(ds_list, dim="trajectory", join="outer", combine_attrs="override")
    comp = dict(zlib=True, complevel=9)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(
        path=store_path + "cluster_merged.nc",
        encoding=encoding,
        compute=True,
        engine="netcdf4",
        format="NETCDF4",
        mode="w",
    )


if __name__ == "__main__":
    import argparse
    import textwrap
    from latexify import *

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Cluster files with k-means from a sensitivity analysis simulation along
            trajectories, e.g., by using
            python cluster_data.py --file /project/meteo/w2w/Z2/Z2_data_gradients/ --out_file /path/to/pics/cluster.png 
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file",
        default="../data/vladiana_ensembles/",
        help=textwrap.dedent(
            """\
            Path to a folder with many files from a sensitivity analysis simulation.
            """
        ),
    )
    parser.add_argument(
        "--out_file",
        default="../pics/plots.png",
        help=textwrap.dedent(
            """\
            Path and name to store plots.
            """
        ),
    )
    parser.add_argument(
        "--width",
        default=24,
        type=float,
        help=textwrap.dedent(
            """\
            Width in inches for histogram plots.
            """
        ),
    )
    parser.add_argument(
        "--height",
        default=12,
        type=float,
        help=textwrap.dedent(
            """\
            Height in inches for histogram plots.
            """
        ),
    )
    parser.add_argument(
        "--only_asc600",
        action="store_true",
        help=textwrap.dedent(
            """\
            Consider only time steps during the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--inoutflow_time",
        default=-1,
        type=int,
        help=textwrap.dedent(
            """\
            Consider only time steps during the fastest ascent and within the given range before (inflow) and after (outflow) of the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--cluster_var",
        default="w",
        help=textwrap.dedent(
            """\
            Cluster w.r.t. this variable.
            If the variable is a model parameter such as drain_a_geo, then clusters for all available model outputs with 
            a sensitivity to the model parameter are generated seperately.
            """
        ),
    )
    parser.add_argument(
        "--k",
        type=int,
        default=5,
        help=textwrap.dedent(
            """\
            Number of clusters.
            """
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=textwrap.dedent(
            """\
            More output regarding the calculations.
            """
        ),
    )

    args = parser.parse_args()
