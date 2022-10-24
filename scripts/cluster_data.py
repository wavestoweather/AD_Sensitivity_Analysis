import warnings

warnings.simplefilter(action="ignore", category=RuntimeWarning)

from matplotlib.figure import Figure
import numpy as np
import os
import pandas as pd
import panel as pn
import seaborn as sns
from sklearn.cluster import KMeans
from tqdm.auto import tqdm
import xarray as xr

try:
    from latexify import param_id_map
except:
    from scripts.latexify import param_id_map


def load_data(file_path, x, only_asc600=False, inoutflow_time=-1, verbose=False):
    """
    Load multiple NetCDF-files as xarray.Dataset, filter them if necessary and concatenate all
    data useful for other steps.

    Parameters
    ----------
    file_path : string
        Path to trajectories from a sensitivity analysis.
    x : string or list of strings
        Model state variable(s) or model parameter(s) for calculating the clusters.
    only_asc600 : bool
        Consider only time steps during the fastest 600 hPa ascent.
    inoutflow_time : int
        Consider only time steps during the fastest 600 hPa ascent and this many timesteps before
        and after the ascent (in- and outflow).
    verbose : bool
        If true, get more output.

    Returns
    -------
    xarray.Dataset with concatenated values for all trajectories. Includes a new dimension 'file' to
    backtrack any results to the correct file and trajectory.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    param_ids = None
    list_of_arrays = []
    if isinstance(x, str):
        vars = [x, "asc600", "phase"]
    else:
        vars = [v for v in x] + ["asc600", "phase"]

    for f in tqdm(files) if verbose else files:
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[vars]
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
            param_ids = ds["Output_Parameter_ID"]
        # fill values where
        data = ds[x].dropna("time", how="all")
        data = data.expand_dims({"file": [f]})
        list_of_arrays.append(data)
    data = xr.concat(list_of_arrays, dim="file", join="outer", combine_attrs="override")
    return data


def get_average(data):
    """
    Calculate the average over all time steps.

    Parameters
    ----------
    data : xarray.Dataset
        A dataset loaded with 'load_data()' although any dataset with a coordinate 'time' works as well.

    Returns
    -------
    xarray.Dataset where the 'time' coordinate is ditched and averages are calculated.
    """
    return data.mean(dim="time", skipna=True, keep_attrs=True)


def get_cluster(data, k, tol=1e-4, x=None, param_names=None, verbose=False):
    """
    Calculate clusters using k-means for the variable 'x'.

    Parameters
    ----------
    data : xarray.Dataset or xarray.DataArray
        A dataset of trajectories from a sensitivity simulation or one column of the set.
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
    verbose : bool
        If true, get more output.

    Returns
    -------
    A pandas.DataFrame with 'cluster' of type KMeans, 'trajectory', and the column defined using 'x'
    or data.name if 'data' is a xarray.DataArray.
    If the cluster is done for sensitivities (model parameters), then the coordinate 'Output Parameter' is added.

    """
    avg_name_list = None
    data_array = None
    n_features = 1
    if x is not None and isinstance(data, xr.Dataset):
        if isinstance(x, str) and param_names is None:
            # A single model state
            data_array = data[x]
            avg_name = f"avg {x}"
        elif isinstance(x, str):
            # A single parameter but for multiple model states
            avg_name_list = []
            n_features = len(param_names)
            for p in param_names:
                avg_name_list.append(f"avg d{p}/{x}")
        elif param_names is None:
            # Multiple model states and no model parameter
            n_features = len(x)
            avg_name_list = [f"avg {v}" for v in x]
        else:
            # Multiple model states and at least one model parameter
            avg_name_list = []
            for v in x:
                if v[0] == "d" and v != "deposition":
                    for p in param_names:
                        avg_name_list.append(f"avg d{p}/{v}")
                else:
                    avg_name_list.append(f"avg {v}")
            n_features = len(avg_name_list)
    else:
        data_array = data
        avg_name = "avg " + data.name
    # n_samples, n_features
    n_samples = len(data["file"]) * len(data["trajectory"])
    shape = (n_samples, n_features)
    if verbose:
        print("Calculate the cluster")

    def get_fit_data():
        fit_data = np.zeros(shape)
        for i in range(n_features):
            name = avg_name_list[i]
            if "/" in name:
                # Model parameter
                out_p = name.split("/")[0][5:]
                in_p = name.split("/")[1]
                param_i = np.argwhere(np.asarray(param_id_map) == out_p).item()
                fit_data[:, i] = np.reshape(
                    data[in_p].sel({"Output_Parameter_ID": param_i}).values, n_samples
                )
            else:
                # Model state
                state_name = name[4:]
                fit_data[:, i] = np.reshape(data[state_name].values, n_samples)
        fit_mask = None
        for i in range(n_features):
            if fit_mask is not None:
                fit_mask = ~np.isnan(fit_data[:, i]) & fit_mask
            else:
                fit_mask = ~np.isnan(fit_data[:, i])
        fit_data_normed = fit_data.copy()
        for i in range(n_features):
            fit_data_normed[:, i] = fit_data_normed[:, i] / np.linalg.norm(
                fit_data_normed[fit_mask, i]
            )
        return fit_data, fit_data_normed, fit_mask

    def ndim_cluster(fit_data, fit_data_normed, fit_mask):
        kmeans_dic = {
            "cluster": np.asarray([]),
            "trajectory": np.asarray([]),
            "file": np.asarray([]),
        }
        for name in avg_name_list:
            kmeans_dic[name] = np.asarray([])
        clusters = KMeans(n_clusters=k, tol=tol, random_state=42).fit_predict(
            fit_data_normed[fit_mask, :]
        )
        trajectory = data["trajectory"]
        for _ in range(len(data["file"]) - 1):
            trajectory = np.append(trajectory, data["trajectory"])
        filenames = np.repeat(data["file"].values, len(data["trajectory"]))
        kmeans_dic["file"] = np.append(kmeans_dic["file"], filenames[fit_mask])
        kmeans_dic["trajectory"] = np.append(
            kmeans_dic["trajectory"], trajectory[fit_mask]
        )
        kmeans_dic["cluster"] = np.append(kmeans_dic["cluster"], clusters)
        for i in range(n_features):
            name = avg_name_list[i]
            kmeans_dic[name] = np.append(
                kmeans_dic[name], fit_data[fit_mask, i].flatten()
            )
        return pd.DataFrame.from_dict(kmeans_dic)

    if param_names is not None:
        fit_data, fit_data_normed, fit_mask = get_fit_data()
        return ndim_cluster(fit_data, fit_data_normed, fit_mask)

    elif n_features == 1:
        fit_data = np.reshape(data_array.values, shape)
        fit_data_normed = fit_data / np.linalg.norm(fit_data)
        clusters = KMeans(n_clusters=k, tol=tol, random_state=42).fit_predict(
            fit_data_normed[~np.isnan(fit_data[:, 0])]
        )
        trajectory = data_array["trajectory"]
        for _ in (
            tqdm(range(len(data_array["file"]) - 1))
            if verbose
            else range(len(data_array["file"]) - 1)
        ):
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
    else:
        fit_data, fit_data_normed, fit_mask = get_fit_data()
        return ndim_cluster(fit_data, fit_data_normed, fit_mask)


def plot_cluster_data(
    data,
    x,
    y=None,
    plot_type="histplot",
    width=16,
    height=12,
    log_scale=(False, False),
    palette="tab10",
    **kwargs,
):
    """
    Use seaborn to create either a scatter plot or a histogram of all datapoints and their cluster in a different color.

    Parameters
    ----------
    data : pandas.DataFrame
        A pandas.DataFrame with 'cluster' of type KMeans, 'trajectory', and the column defined in 'x' and 'y'.
    x : string
        The column of 'data' for the x-axis for the plots.
    y : string or None
        If plot_type == 'scatter': The column of 'data' for the y-axis for the scatter plot.
    plot_type : string
        Use 'histplot' for a histogram and 'scatter' for a scatter plot.
    width : int
        Width of the plot in pixels.
    height : int
        Height of the plot in pixels.
    log_scale : tuple of bool
        If true, the x- and/or y-axis are plotted using a log-scale.
    palette : string
        Palette for cluster colors.
    kwargs : dictionary
        Arguments passed down to matplotlib plotting functions.

    Returns
    -------
    matplotlib.axes.Axes created using seaborn plot function.
    """
    sns.set(rc={"figure.figsize": (width, height)})

    if plot_type == "histplot":
        return sns.histplot(
            data=data,
            x=x,
            hue="cluster",
            palette=palette,
            log_scale=log_scale,
            **kwargs,
        )
    elif plot_type == "scatter":
        g = sns.scatterplot(
            data=data, x=x, y=y, hue="cluster", palette=palette, **kwargs
        )
        if log_scale[0]:
            g.set_xscale("log")
        if log_scale[1]:
            g.set_yscale("log")
        return g


def plot_cluster_data_interactive(data):
    """
    Calling this function from a Jupyter notebook allows to visualize the cluster association with different
    dimensions. Make sure to call pn.extension() from your notebook first.

    Parameters
    ----------
    data

    Returns
    -------

    """
    out_params = []
    in_params = []
    for col in data:
        if "/" in col:
            out_params.append(col.split("/")[0][5:])
            in_params.append(col.split("/")[1])
        elif "avg " in col:
            in_params.append(col[4:])
        else:
            in_params.append(col)
    in_params = list(set(in_params))
    out_params = list(set(out_params))
    if len(out_params) == 0:
        out_params = ["Not available"]
    out_param_x = pn.widgets.RadioButtonGroup(
        name="Output Parameter (if any) for the x-axis",
        value=out_params[0],
        options=out_params,
        button_type="primary",
    )
    in_param_x = pn.widgets.Select(
        name="Model parameter or model state for the x-axis",
        value=in_params[0],
        options=in_params,
        button_type="default",
    )
    out_param_y = pn.widgets.RadioButtonGroup(
        name="Output Parameter (if any) for the y-axis",
        value=out_params[0],
        options=out_params,
        button_type="primary",
    )
    in_param_y = pn.widgets.Select(
        name="Model parameter or model state for the y-axis",
        value=in_params[1],
        options=in_params,
        button_type="default",
    )
    width_slider = pn.widgets.IntSlider(
        name="Width in inches",
        start=3,
        end=15,
        step=1,
        value=9,
    )
    height_slider = pn.widgets.IntSlider(
        name="Height in inches",
        start=3,
        end=15,
        step=1,
        value=6,
    )
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    logx_plot = pn.widgets.Toggle(
        name="Use log-scale for the x-axis",
        value=False,
        button_type="success",
    )
    logy_plot = pn.widgets.Toggle(
        name="Use log-scale for the y-axis",
        value=False,
        button_type="success",
    )
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=2,
        step=0.1,
        value=0.7,
    )

    def get_plot(
        data,
        in_p_x,
        out_p_x,
        in_p_y,
        out_p_y,
        logx,
        logy,
        width,
        height,
        font_scale,
        title,
    ):
        if in_p_x == in_p_y and out_p_x == out_p_y:
            return
        sns.set(rc={"figure.figsize": (width, height)})
        fig = Figure()
        ax = fig.subplots()
        if in_p_x[0] == "d" and in_p_x != "deposition":
            x = f"avg d{out_p_x}/{in_p_x}"
        elif in_p_x != "cluster" and in_p_x != "trajectory" and in_p_x != "file":
            x = f"avg {in_p_x}"
        else:
            x = in_p_x
        if in_p_y[0] == "d" and in_p_y != "deposition":
            y = f"avg d{out_p_y}/{in_p_y}"
        elif in_p_y != "cluster" and in_p_y != "trajectory" and in_p_y != "file":
            y = f"avg {in_p_y}"
        else:
            y = in_p_y
        sns.scatterplot(
            data=data,
            x=x,
            y=y,
            hue="cluster",
            palette="tab10",
            ax=ax,
        )
        if logx:
            ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=int(10 * font_scale),
        )
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
        return fig

    plot_pane = pn.panel(
        pn.bind(
            get_plot,
            data=data,
            in_p_x=in_param_x,
            out_p_x=out_param_x,
            in_p_y=in_param_y,
            out_p_y=out_param_y,
            logx=logx_plot,
            logy=logy_plot,
            width=width_slider,
            height=height_slider,
            font_scale=font_slider,
            title=title_widget,
        ),
    ).servable()

    return pn.Column(
        pn.Row(
            out_param_x,
            in_param_x,
            logx_plot,
        ),
        pn.Row(
            out_param_y,
            in_param_y,
            logy_plot,
        ),
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
        ),
        plot_pane,
    )


def get_traj_near_center(data, x, n=1):
    """
    Extract the number and filenames for trajectories closest to the cluster centers.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame with clusters from 'get_cluster()'.
    x : string
        Variable that has been used for calculating the clusters. Usually starts with 'avg '.
    n : int
        The number of trajectories to extract for each cluster center.

    Returns
    -------
    np.array of pandas.DataFrame with the closest trajectories where each DataFrame is for a different cluster center.
    """
    centers = []
    for cluster in np.unique(data["cluster"]):
        centers.append(np.mean(data.loc[data["cluster"] == cluster][x]))
    traj_centers = []
    for center in centers:
        idx = (np.abs(data[x] - center)).argsort()[:n]
        traj_centers.append(data.iloc[idx])
    return traj_centers


def extract_trajs(
    centers,
    file_path,
    store_path,
    only_asc600=False,
    inoutflow_time=-1,
    traj_abs=False,
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
    verbose : bool
        If true, get more output.
    """
    ds_list = []
    traj_enum = 0
    cluster = 0
    min_time = None
    max_time = None
    if verbose:
        print("Load the files to check for the time ranges for the final file")
    for center_df in tqdm(centers) if verbose else centers:
        for idx, row in (
            tqdm(center_df.iterrows(), leave=False, total=len(center_df))
            if verbose
            else tqdm(center_df.iterrows(), leave=False, total=len(center_df))
        ):
            f = row["file"]
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
            tmp_min = ds["time"].min()
            if min_time is None or min_time > tmp_min:
                min_time = tmp_min
            tmp_max = ds["time"].max()
            if max_time is None or max_time < tmp_max:
                max_time = tmp_max

    for center_df in tqdm(centers) if verbose else centers:
        data_arrays = []
        for idx, row in (
            tqdm(center_df.iterrows(), leave=False, total=len(center_df))
            if verbose
            else center_df.iterrows()
        ):
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
    if verbose:
        print("Concatenating the final list")
    ds = xr.concat(ds_list, dim="trajectory", join="outer", combine_attrs="override")
    ds["cluster"] = ds["cluster"].where(~np.isnan(ds["pressure"]))
    if traj_abs:
        for col in ds:
            ds[col] = np.abs(ds[col])
    if store_path[-1] != "/" and "." not in store_path:
        # We assume store_path is just a path without a proper filename.
        filename = store_path + "/" + "cluster_extracted_trajectories.nc"
    elif store_path[-1] == "/":
        filename = store_path + "cluster_extracted_trajectories.nc"
    else:
        filename = store_path
    if verbose:
        print(f"Storing the extracted trajectories to {filename}")
    comp = dict(zlib=True, complevel=9)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(
        path=filename,
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
        "--file_path",
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
        "--logx",
        action="store_true",
        help=textwrap.dedent(
            """\
            Plot the x-axis using log-scale.
            """
        ),
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help=textwrap.dedent(
            """\
            Plot the y-axis using log-scale.
            """
        ),
    )
    parser.add_argument(
        "--plot_type",
        type=str,
        default="histplot",
        choices=["histplot", "scatter", "both", "none"],
        help=textwrap.dedent(
            """\
            You may choose different kind of plots.
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
    data = load_data(
        file_path=args.file_path,
        x=args.cluster_var,
        only_asc600=args.only_asc600,
        inoutflow_time=args.inoutflow_time,
        verbose=args.verbose,
    )
    data = get_average(
        data
    )  # We don't need the original dataset with all individual time steps anymore.
    clusters = get_cluster(
        data=data,
        k=args.k,
        x=args.cluster_var,
        param_names=args.sens_model_states,
        verbose=args.verbose,
    )
    if args.plot_type == "histplot" or args.plot_type == "both":
        for c in args.cluster_var:
            plot = plot_cluster_data(
                data=clusters,
                x=f"avg {c}",
                plot_type="histplot",
                log_scale=(args.logx, args.logy),
                width=args.width,
                height=args.height,
                **{
                    "multiple": "stack"
                },  # For demonstration purpose, we add another keyword like that.
            )
            if args.plot_type == "both":
                out_file = args.out_file.split(".pn")[0] + f"_{c}_histogr.png"
            else:
                out_file = args.out_file.split(".pn")[0] + f"_{c}.png"

            plot.get_figure().savefig(out_file, dpi=300, bbox_inches="tight")
            plot.clear()

    if args.plot_type == "scatter" or args.plot_type == "both":
        for c in args.cluster_var:
            plot = plot_cluster_data(
                data=clusters,
                x="trajectory",
                y=f"avg {c}",
                plot_type="scatter",
                log_scale=(args.logx, args.logy),
                width=args.width,
                height=args.height,
            )
            if args.plot_type == "both":
                out_file = args.out_file.split(".pn")[0] + f"_{c}_scatter.png"
            else:
                out_file = args.args.out_file.split(".pn")[0] + f"_{c}.png"
            plot.get_figure().savefig(out_file, dpi=300, bbox_inches="tight")

    # It would be nice to see the trajectories near the center in a separate file
    if args.extract_n_trajs > 0 and args.extract_store_path is not None:
        if args.verbose:
            print("Get the trajectory indices near the center")
        centers = get_traj_near_center(
            data=clusters,
            x=f"avg {args.cluster_var}",
            n=args.extract_n_trajs,
        )
        if args.extracted_trajs_complete:
            extract_trajs(
                centers=centers,
                file_path=args.file_path,
                store_path=args.extract_store_path,
                traj_abs=args.extracted_trajs_abs,
                verbose=args.verbose,
            )
        else:
            extract_trajs(
                centers=centers,
                file_path=args.file_path,
                store_path=args.extract_store_path,
                only_asc600=args.only_asc600,
                inoutflow_time=args.inoutflow_time,
                traj_abs=args.extracted_trajs_abs,
                verbose=args.verbose,
            )
