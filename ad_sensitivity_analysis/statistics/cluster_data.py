"""Create clusters for trajectories.

"""
import warnings

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from tqdm.auto import tqdm
import xarray as xr

from ad_sensitivity_analysis.data_handler.loader import load_data
from ad_sensitivity_analysis.data_handler.filter_data import filter_trajectories
from ad_sensitivity_analysis.data_handler.transform import (
    fix_coords,
    set_col_types_fp32,
    get_average,
    fill_time_coords,
)
from ad_sensitivity_analysis.plot.latexify import param_id_map
from ad_sensitivity_analysis.statistics.cluster_parse_aux import (
    parse_ds_single_model_state,
    parse_ds_single_param,
    parse_ds_model_states,
    parse_ds_model_states_params,
)

warnings.simplefilter(action="ignore", category=RuntimeWarning)


def get_time_limits_coord(center_dfs, file_path, verbose):
    """

    Parameters
    ----------
    center_dfs
    file_path
    verbose

    Returns
    -------

    """
    min_time = None
    max_time = None
    out_coord = "Output_Parameter_ID"
    if verbose:
        print("Load the files to check for the time ranges for the final file")
    for center_df in tqdm(center_dfs) if verbose else center_dfs:
        for _, row in (
            tqdm(center_df.iterrows(), leave=False, total=len(center_df))
            if verbose
            else center_df.iterrows()
        ):
            f = row["file"]
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
            if out_coord not in ds:
                out_coord = "Output Parameter"
            tmp_min = ds["time"].min()
            if min_time is None or min_time > tmp_min:
                min_time = tmp_min
            tmp_max = ds["time"].max()
            if max_time is None or max_time < tmp_max:
                max_time = tmp_max
    return min_time, max_time, out_coord


def cluster_to_xarray(cluster, times, traj_enum):
    """

    Parameters
    ----------
    cluster
    times
    traj_enum

    Returns
    -------

    """
    return xr.DataArray(
        data=[[np.repeat(np.float32(cluster), len(times))]],
        coords={
            "ensemble": [0],
            "trajectory": [traj_enum],
            "time": times,
        },
        dims=["ensemble", "trajectory", "time"],
        name="cluster",
        attrs={
            "auxiliary_data": "yes",
            "_FillValue": np.NaN,
            "standard_name": "cluster",
            "long_name": "cluster",
        },
    )


def get_fit_data(data, n_samples, n_features, avg_name_list, reduce_name, out_coord):
    """

    Parameters
    ----------
    data
    n_samples
    n_features
    avg_name_list
    reduce_name
    out_coord

    Returns
    -------

    """
    fit_data = np.zeros((n_samples, n_features))
    for i in range(n_features):
        if "/" in avg_name_list[i]:
            # Model parameter
            out_p = avg_name_list[i].split("/")[0][len(reduce_name) + 1 :]
            in_p = avg_name_list[i].split("/")[1]
            if out_coord == "Output_Parameter_ID":
                param_i = np.argwhere(np.asarray(param_id_map) == out_p).item()
            else:
                param_i = out_p
            fit_data[:, i] = np.reshape(
                data[in_p].sel({out_coord: param_i}).values, n_samples
            )
        else:
            # Model state
            state_name = avg_name_list[i][len(reduce_name) :]
            fit_data[:, i] = np.reshape(data[state_name].values, n_samples)
    fit_mask = None
    for i in range(n_features):
        if fit_mask is not None:
            fit_mask = ~np.isnan(fit_data[:, i]) & fit_mask
        else:
            fit_mask = ~np.isnan(fit_data[:, i])
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    fit_data_normed = scaler.fit_transform(fit_data)
    return fit_data, fit_data_normed, fit_mask


def ndim_cluster(
    data,
    avg_name_list,
    fit_data,
    out_coord,
    non_features_list,
    fit_data_normed,
    fit_mask,
    k,
    tol,
    include_all_data,
):
    """

    Parameters
    ----------
    data
    avg_name_list
    fit_data
    out_coord
    non_features_list
    fit_data_normed
    fit_mask
    k
    tol
    include_all_data

    Returns
    -------

    """
    kmeans_dic = {
        "cluster": np.asarray([]),
        "trajectory": np.asarray([]),
        "file": np.asarray([]),
    }
    trajectory = data["trajectory"]
    for _ in range(len(data["file"]) - 1):
        trajectory = np.append(trajectory, data["trajectory"])
    filenames = np.repeat(data["file"].values, len(data["trajectory"]))
    kmeans_dic["file"] = np.append(kmeans_dic["file"], filenames[fit_mask])
    kmeans_dic["trajectory"] = np.append(kmeans_dic["trajectory"], trajectory[fit_mask])
    kmeans_dic["cluster"] = np.append(
        kmeans_dic["cluster"],
        KMeans(n_clusters=k, tol=tol, init="k-means++", random_state=42).fit_predict(
            fit_data_normed[fit_mask, :]
        ),
    )
    for i in range(fit_data.shape()[1]):
        kmeans_dic[avg_name_list[i]] = fit_data[fit_mask, i].flatten()
    if include_all_data:
        for non_feat in non_features_list:
            if non_feat[1] is None:
                kmeans_dic[non_feat[0]] = np.asarray(
                    data[non_feat[2]].values
                ).flatten()[fit_mask]
            else:
                kmeans_dic[non_feat[0]] = np.asarray(
                    data.sel({out_coord: non_feat[1]})[non_feat[2]].values
                ).flatten()[fit_mask]
    return pd.DataFrame.from_dict(kmeans_dic)


def dataframe_from_single_feature(
    data_array,
    data,
    out_coord,
    shape,
    k,
    tol,
    avg_name,
    include_all_data,
    non_features_list,
    verbose,
):
    """

    Parameters
    ----------
    data_array
    data
    out_coord
    shape
    k
    tol
    avg_name
    include_all_data
    non_features_list
    verbose

    Returns
    -------

    """
    fit_data = np.reshape(data_array.values, shape)
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    trajectory = data_array["trajectory"]
    for _ in (
        tqdm(range(len(data_array["file"]) - 1))
        if verbose
        else range(len(data_array["file"]) - 1)
    ):
        trajectory = np.append(trajectory, data_array["trajectory"])
    kmeans_dic = {
        "cluster": KMeans(n_clusters=k, tol=tol, random_state=42).fit_predict(
            scaler.fit_transform(fit_data)[~np.isnan(fit_data[:, 0])]
        ),
        "trajectory": trajectory[~np.isnan(fit_data[:, 0])],
        avg_name: fit_data[~np.isnan(fit_data[:, 0])].flatten(),
        "file": np.repeat(data_array["file"].values, len(data_array["trajectory"]))[
            ~np.isnan(fit_data[:, 0])
        ],
    }
    if include_all_data:
        for non_feat in non_features_list:
            if non_feat[1] is None:
                kmeans_dic[non_feat[0]] = np.asarray(
                    data[non_feat[2]].values
                ).flatten()[~np.isnan(fit_data[:, 0])]
            else:
                kmeans_dic[non_feat[0]] = np.asarray(
                    data.sel({out_coord: non_feat[1]})[non_feat[2]].values
                ).flatten()[~np.isnan(fit_data[:, 0])]
    return pd.DataFrame.from_dict(kmeans_dic)


# pylint: disable=too-many-locals
def get_cluster(
    data,
    k,
    tol=1e-4,
    x=None,
    param_names=None,
    include_all_data=False,
    reduce_name="",
    verbose=False,
):
    """
    Calculate clusters using k-means for the variable 'x'.

    Parameters
    ----------
    data : xarray.Dataset or xarray.DataArray
        A dataset of trajectories from a sensitivity simulation or one column of the set
        where the time dimension has been reduced using an average metric.
        Must be a dataset with coordinates 'trajectory' and 'file', i.e., loaded using 'load_data()'.
    k : int
        Number of clusters for k-means.
    tol : float
        From scikit-learn.org: Relative tolerance with regards to Frobenius norm of the difference
        in the cluster centers of two consecutive iterations to declare convergence.
    x : string (needed if data is xarray.Dataset) or list of string
        The model state variable or model parameter to create clusters for. If multiple variables or parameters shall
        be considered, then x must be a list of strings.
    param_names : list of strings
        If 'x' is a model parameter, then calculate clusters for each sensitivity of a model state variable towards 'x'.
    include_all_data : bool
        Include all columns in data in the returned dataframe. Otherwise only include columns used for clustering.
    reduce_name : string
        Name to be prepended to the columns. Should relate to the reduction
        applied to the dataset, such as "avg" or "rank".
    verbose : bool
        If true, get more output.

    Returns
    -------
    A pandas.DataFrame with 'cluster' of type KMeans, 'trajectory', and the column defined using 'x'
    or data.name if 'data' is a xarray.DataArray.
    If the cluster is done for sensitivities (model parameters), then the coordinate 'Output Parameter' is added.

    """
    avg_name_list = None
    n_features = 1
    out_coord = "Output_Parameter_ID"
    if out_coord not in data:
        out_coord = "Output Parameter"
    if x is not None and isinstance(data, xr.Dataset):
        if (isinstance(x, str) or len(x) == 1) and param_names is None:
            data_array, non_features_list, avg_name = parse_ds_single_model_state(
                data_trajs=data,
                x=x,
                out_coord=out_coord,
                reduce_name=reduce_name,
                include_all_data=include_all_data,
            )
        elif isinstance(x, str) or len(x) == 1:
            n_features, non_features_list, avg_name_list = parse_ds_single_param(
                data_trajs=data,
                x=x,
                out_coord=out_coord,
                reduce_name=reduce_name,
                include_all_data=include_all_data,
                param_names=param_names,
            )
        elif param_names is None:
            n_features, non_features_list, avg_name_list = parse_ds_model_states(
                data_trajs=data,
                x=x,
                out_coord=out_coord,
                reduce_name=reduce_name,
                include_all_data=include_all_data,
            )

        else:
            n_features, non_features_list, avg_name_list = parse_ds_model_states_params(
                data_trajs=data,
                x=x,
                out_coord=out_coord,
                reduce_name=reduce_name,
                include_all_data=include_all_data,
                param_names=param_names,
            )

    else:
        data_array = data
        avg_name = f"{reduce_name}" + data.name
    n_samples = len(data["file"]) * len(data["trajectory"])
    if verbose:
        print("Calculate the cluster")

    if param_names is not None:
        fit_data, fit_data_normed, fit_mask = get_fit_data(
            data=data,
            n_samples=n_samples,
            n_features=n_features,
            avg_name_list=avg_name_list,
            reduce_name=reduce_name,
            out_coord=out_coord,
        )
        return ndim_cluster(
            data=data,
            avg_name_list=avg_name_list,
            fit_data=fit_data,
            out_coord=out_coord,
            non_features_list=non_features_list,
            fit_data_normed=fit_data_normed,
            fit_mask=fit_mask,
            k=k,
            tol=tol,
            include_all_data=include_all_data,
        )

    if n_features == 1:
        return dataframe_from_single_feature(
            data_array=data_array,
            data=data,
            out_coord=out_coord,
            shape=(n_samples, n_features),
            k=k,
            tol=tol,
            avg_name=avg_name,
            include_all_data=include_all_data,
            non_features_list=non_features_list,
            verbose=verbose,
        )

    fit_data, fit_data_normed, fit_mask = get_fit_data(
        data=data,
        n_samples=n_samples,
        n_features=n_features,
        avg_name_list=avg_name_list,
        reduce_name=reduce_name,
        out_coord=out_coord,
    )
    return ndim_cluster(
        data=data,
        avg_name_list=avg_name_list,
        fit_data=fit_data,
        out_coord=out_coord,
        non_features_list=non_features_list,
        fit_data_normed=fit_data_normed,
        fit_mask=fit_mask,
        k=k,
        tol=tol,
        include_all_data=include_all_data,
    )


def get_traj_near_center(data, x, n=1):
    """
    Extract the number and filenames for trajectories closest to the cluster centers.
    Does not normalize data for distance calculation.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame with clusters from 'get_cluster()'.
    x : string or list of strings
        Variable(s) that has been used for calculating the clusters. Usually starts with 'avg '.
    n : int
        The number of trajectories to extract for each cluster center.

    Returns
    -------
    np.array of pandas.DataFrame with the closest trajectories where each DataFrame is for a different cluster center.
    """
    centers = []
    for cluster in np.unique(data["cluster"]):
        if isinstance(x, str):
            centers.append(np.mean(data.loc[data["cluster"] == cluster][x]))
        else:
            centers_tmp = []
            for col in x:
                centers_tmp.append(np.mean(data.loc[data["cluster"] == cluster][col]))
            centers.append(centers_tmp)
    traj_centers = []
    for center in centers:
        if isinstance(x, str):
            idx = (np.abs(data[x] - center)).argsort()[:n]
            traj_centers.append(data.iloc[idx])
        else:
            distance = np.power(data[x[0]] - center[0], 2)
            for i in range(1, len(center)):
                distance += np.power(data[x[i]] - center[i], 2)
            distance = np.sqrt(distance)
            idx = distance.argsort()[:n]
            traj_centers.append(data.iloc[idx])
    return traj_centers


def save_extracted_trajs(ds, store_path, split_states, out_coord, verbose):
    """

    Parameters
    ----------
    ds
    store_path
    split_states
    out_coord
    verbose

    Returns
    -------

    """
    if store_path[-1] != "/" and "." not in store_path:
        # We assume store_path is just a path without a proper filename.
        filename = store_path + "/" + "cluster_extracted_trajectories.nc"
    elif store_path[-1] == "/":
        filename = store_path + "cluster_extracted_trajectories.nc"
    else:
        filename = store_path
    if verbose:
        print(f"Storing the extracted trajectories to {filename}")
    comp = {"zlib": True, "complevel": 9}
    if split_states and out_coord in ds.dims:
        for out_param in ds[out_coord]:
            ds_tmp = ds.sel({out_coord: out_param})
            encoding = {var: comp for var in ds_tmp.data_vars}
            ds_tmp.to_netcdf(
                path=filename.split(".")[0] + f"_{out_param.values.item()}.nc",
                encoding=encoding,
                compute=True,
                engine="netcdf4",
                format="NETCDF4",
                mode="w",
            )
    else:
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(
            path=filename,
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )


def set_coords(ds, keep_coords):
    """

    Parameters
    ----------
    ds
    keep_coords

    Returns
    -------

    """
    if keep_coords == "relative":
        cols = []
        for col in ds:
            if col not in ("lon", "lat"):
                cols.append(col)
        ds = ds[cols]
    elif keep_coords == "normal":
        cols = []
        for col in ds:
            if col not in ("relative_lon", "relative_lat"):
                cols.append(col)
        ds = ds[cols]
    return ds


def extract_trajs(
    centers,
    file_path,
    store_path,
    only_asc600=False,
    inoutflow_time=-1,
    traj_abs=False,
    split_states=False,
    keep_coords="both",
    fix_weird_coords=False,
    verbose=False,
):
    """
    Given trajectories from 'get_traj_near_center()', extract the original data from the complete dataset and merge them
    to a new dataset. The new dataset includes 'cluster' to identify the affiliation of each trajectory.

    Parameters
    ----------
    centers : np.array of pandas.DataFrame
        For each cluster one pandas.DataFrame with the trajectories closest to the cluster center.
    file_path : string
        Path to trajectories from a sensitivity analysis.
    store_path : string
        Filepath or filepath and name to store the merged dataset as NetCDF-file.
    only_asc600 : bool
        Consider only time steps during the fastest 600 hPa ascent.
    inoutflow_time : int
        Consider only time steps during the fastest 600 hPa ascent and this many timesteps before
        and after the ascent (in- and outflow).
    traj_abs : bool
        Take the absolute values for sensitivities.
    split_states : bool
        Removes the dimension 'Output_Parameter_ID' or 'Output Parameter' if present and
        stores separate NetCDF-files for each dimension value.
    keep_coords : string
        "both": Keep longitude, latitude, and coordinates relative to the start of the ascent.
        "relative": Keep only coordinates relative to the start of the ascent.
        "normal": Keep only longitude and latitude.
    fix_weird_coords : bool
        Some data may have zero lon and lat where no data should be there.
    verbose : bool
        If true, get more output.
    """
    ds_list = []
    traj_enum = 0
    min_time, max_time, out_coord = get_time_limits_coord(centers, file_path, verbose)
    if verbose:
        print(
            "For each center, get all the trajectories and merge them into one dataset."
        )
    for center_df in tqdm(centers) if verbose else centers:
        data_arrays = []
        cluster = np.unique(center_df["cluster"])[0]
        for _, row in (
            tqdm(center_df.iterrows(), leave=False, total=len(center_df))
            if verbose
            else center_df.iterrows()
        ):
            traj = row["trajectory"]
            f = row["file"]
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
            ds = ds.sel({"trajectory": [traj]})
            ds = filter_trajectories(
                ds=ds,
                only_asc600=only_asc600,
                inoutflow_time=inoutflow_time,
            )
            if "asc600" in ds:
                ds["asc600"] = ds["asc600"].where(ds["asc600"] >= 0)
            traj_attrs = ds["trajectory"].attrs
            ens_attrs = ds["ensemble"].attrs
            time_attrs = ds["time"].attrs
            ds["trajectory"] = [traj_enum]
            ds["cluster"] = cluster_to_xarray(cluster, ds["time"].values, traj_enum)
            ds["trajectory"].attrs = traj_attrs
            if "asc600" not in ds:
                ds["ensemble"] = [0]
            ds["ensemble"].attrs = ens_attrs
            ds["time"].attrs = time_attrs
            if fix_weird_coords:
                ds = fix_coords(ds)
            ds = fill_time_coords(ds, min_time, max_time, out_coord, traj_enum)
            ds = set_col_types_fp32(ds)
            traj_enum += 1
            data_arrays.append(ds)
        ds_list.append(
            xr.concat(
                data_arrays, dim="trajectory", join="outer", combine_attrs="override"
            )
        )
    if verbose:
        print("Concatenating the final list")
    ds = xr.concat(ds_list, dim="trajectory", join="outer", combine_attrs="override")
    ds["cluster"] = ds["cluster"].where(~np.isnan(ds["pressure"]))
    if traj_abs:
        for col in ds:
            if col[0] == "d" and col != "deposition":
                ds[col] = np.abs(ds[col])
    ds = set_coords(ds, keep_coords)
    save_extracted_trajs(
        ds=ds,
        store_path=store_path,
        split_states=split_states,
        out_coord=out_coord,
        verbose=verbose,
    )


def main(arguments):
    """

    Parameters
    ----------
    arguments

    Returns
    -------

    """
    data = load_data(
        file_path=arguments.file_path,
        x=arguments.cluster_var,
        only_asc600=arguments.only_asc600,
        inoutflow_time=arguments.inoutflow_time,
        verbose=arguments.verbose,
    )
    data = get_average(
        data
    )  # We don't need the original dataset with all individual time steps anymore.
    clusters = get_cluster(
        data=data,
        k=arguments.k,
        x=arguments.cluster_var,
        param_names=arguments.sens_model_states,
        reduce_name="avg ",
        verbose=arguments.verbose,
    )
    # It would be nice to see the trajectories near the center in a separate file
    if arguments.extract_n_trajs > 0 and arguments.extract_store_path is not None:
        if arguments.verbose:
            print("Get the trajectory indices near the center")
        centers = get_traj_near_center(
            data=clusters,
            x=f"avg {arguments.cluster_var}",
            n=arguments.extract_n_trajs,
        )
        if arguments.extracted_trajs_complete:
            extract_trajs(
                centers=centers,
                file_path=arguments.file_path,
                store_path=arguments.extract_store_path,
                traj_abs=arguments.extracted_trajs_abs,
                verbose=arguments.verbose,
            )
        else:
            extract_trajs(
                centers=centers,
                file_path=arguments.file_path,
                store_path=arguments.extract_store_path,
                only_asc600=arguments.only_asc600,
                inoutflow_time=arguments.inoutflow_time,
                traj_abs=arguments.extracted_trajs_abs,
                verbose=arguments.verbose,
            )


if __name__ == "__main__":
    import argparse
    import textwrap

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
        "--file_path",
        default="../data/vladiana_ensembles/",
        help=textwrap.dedent(
            """\
            Path to a folder with many files from a sensitivity analysis simulation.
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
            Consider only time steps during the fastest ascent and 
            within the given range before (inflow) and after (outflow) of the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--cluster_var",
        default=None,
        type=str,
        nargs="+",
        required=True,
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
        "--extract_n_trajs",
        type=int,
        default=0,
        help=textwrap.dedent(
            """\
            Extract the given number of trajectories for each cluster center from the complete dataset 
            that are close to the centers.
            """
        ),
    )
    parser.add_argument(
        "--extract_store_path",
        type=str,
        default=None,
        help=textwrap.dedent(
            """\
            The path (optional: and the name) to store the extracted trajectories. 
            """
        ),
    )
    parser.add_argument(
        "--extracted_trajs_complete",
        action="store_true",
        help=textwrap.dedent(
            """\
            If this option is set, then 'only_asc600' and 'inoutflow_time' are ignored for the extracted trajectories.
            This does not change the clusters.
            """
        ),
    )
    parser.add_argument(
        "--extracted_trajs_abs",
        action="store_true",
        help=textwrap.dedent(
            """\
          If this option is set, store sensitivities as absolute values.
          This is useful for plotting trajectories in Met3D with a log-scale.
          """
        ),
    )
    parser.add_argument(
        "--sens_model_states",
        default=None,
        type=str,
        nargs="+",
        help=textwrap.dedent(
            """\
            If --cluster_var is a sensitivity (model parameter), then you may define the model state variable here 
            to calculate and plot only averages and their clusters regarding these model state variables.
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
    main(args)
