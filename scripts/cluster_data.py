import warnings

warnings.simplefilter(action="ignore", category=RuntimeWarning)

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import panel as pn
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from tqdm.auto import tqdm
import xarray as xr

try:
    from latexify import param_id_map, parse_word
    from latexify import mappings as latex_mappings
except:
    from scripts.latexify import param_id_map, parse_word
    from scripts.latexify import mappings as latex_mappings


def load_data(
    file_path,
    x,
    only_asc600=False,
    inoutflow_time=-1,
    phase=None,
    averages=False,
    verbose=False,
):
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
    averages : bool
        If true, calculate the averages over time when each file is loaded. Otherwise tries to concatenate
        all data in the end which may need lots of memory if inoutflow_time > 0 or only_asc600 == True
        or phase is not None.
    phase : int
        Only use the given phase.
        0: warm phase
        1: mixed phase
        2: ice phase
        3: neutral phase
    verbose : bool
        If true, get more output.

    Returns
    -------
    xarray.Dataset with concatenated values for all trajectories. Includes a new dimension 'file' to
    backtrack any results to the correct file and trajectory.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    list_of_arrays = []
    if isinstance(x, str):
        variables = [x, "asc600", "phase"]
    else:
        variables = [v for v in x] + ["asc600", "phase"]

    for f in tqdm(files) if verbose else files:
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")[
            variables
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

        if phase is not None:
            ds = ds.where(ds["phase"] == phase)
        if (
            "lat" in variables
            or "lon" in variables
            or "relative_lat" in variables
            or "relative_lon" in variables
        ):
            ds = ds.where(ds["lat"] != 0)
        ds = ds[x].dropna("time", how="all")
        ds = ds.expand_dims({"file": [f]})

        if averages:
            # for key in variables:
            #     if "d" != key[0] or "deposition" == key:
            #         break
            # if "d" == key[0] and "deposition" != key:
            #     if "Output Parameter" in ds.dims:
            #         param_coord = "Output Parameter"
            #     else:
            #         param_coord = "Output_Parameter_ID"
            #     weights = (~np.isnan(ds.isel({param_coord: 0})[key])).sum(dim="time")
            # else:
            #     weights = (~np.isnan(ds[key])).sum(dim="time")
            # ds["weights"] = weights
            ds = get_average(ds)
            list_of_arrays.append(ds)
        elif len(ds["time"]) > 0:
            list_of_arrays.append(ds)
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
    non_features_list = []
    data_array = None
    n_features = 1
    out_coord = "Output_Parameter_ID"
    if out_coord not in data:
        out_coord = "Output Parameter"
    if x is not None and isinstance(data, xr.Dataset):
        rank_set = False
        for col in data:
            if "rank" in col or "avg" in col:
                rank_set = True
        if (isinstance(x, str) or len(x) == 1) and param_names is None:
            # A single model state
            if not isinstance(x, str):
                x = x[0]
            data_array = data[x]
            avg_name = f"{reduce_name}{x}"
            if include_all_data:
                if out_coord in data:
                    out_params = data[out_coord].values
                for col in data:
                    if col != x:
                        if (
                            col[0] == "d"
                            and col != "deposition"
                            and col != "deposition rank"
                        ):
                            for out_p in out_params:
                                if out_coord == "Output_Parameter_ID":
                                    non_features_list.append(
                                        [
                                            f"{reduce_name}d{param_id_map[out_p]}/{col}",
                                            out_p,
                                            col,
                                        ]
                                    )
                                else:
                                    non_features_list.append(
                                        [f"{reduce_name}d{out_p}/{col}", out_p, col]
                                    )
                        else:
                            non_features_list.append([f"{reduce_name}{col}", None, col])
        elif isinstance(x, str) or len(x) == 1:
            if not isinstance(x, str):
                x = x[0]
            # A single parameter but for multiple model states
            avg_name_list = []
            n_features = len(param_names)
            for p in param_names:
                avg_name_list.append(f"{reduce_name}d{p}/{x}")
            if include_all_data:
                if out_coord in data:
                    out_params = data[out_coord].values
                for col in data:
                    if col != x:
                        if (
                            col[0] == "d"
                            and col != "deposition"
                            and col != "deposition rank"
                        ):
                            for out_p in out_params:
                                if out_coord == "Output_Parameter_ID":
                                    non_features_list.append(
                                        [
                                            f"{reduce_name}d{param_id_map[out_p]}/{col}",
                                            out_p,
                                            col,
                                        ]
                                    )
                                else:
                                    non_features_list.append(
                                        [f"{reduce_name}d{out_p}/{col}", out_p, col]
                                    )
                        else:
                            non_features_list.append([f"{reduce_name}{col}", None, col])
        elif param_names is None:
            # Multiple model states and no model parameter
            n_features = len(x)
            avg_name_list = [f"{reduce_name}{v}" for v in x]
            if include_all_data:
                if out_coord in data:
                    out_params = data[out_coord].values
                else:
                    out_params = ["model_state"]
                for col in data:
                    if col not in x:
                        if (
                            col[0] == "d"
                            and col != "deposition"
                            and col != "deposition rank"
                        ):
                            for out_p in out_params:
                                if out_coord == "Output_Parameter_ID":
                                    non_features_list.append(
                                        [
                                            f"{reduce_name}d{param_id_map[out_p]}/{col}",
                                            out_p,
                                            col,
                                        ]
                                    )
                                else:
                                    non_features_list.append(
                                        [f"{reduce_name}d{out_p}/{col}", out_p, col]
                                    )
                        else:
                            non_features_list.append([f"{reduce_name}{col}", None, col])
        else:
            # Multiple model states and at least one model parameter
            avg_name_list = []
            for v in x:
                if v[0] == "d" and v != "deposition":
                    for p in param_names:
                        avg_name_list.append(f"{reduce_name}d{p}/{v}")
                else:
                    avg_name_list.append(f"{reduce_name}{v}")
            n_features = len(avg_name_list)
            if include_all_data:
                if out_coord in data:
                    out_params = data[out_coord].values
                for col in data:
                    if col not in x:
                        if (
                            col[0] == "d"
                            and col != "deposition"
                            and col != "deposition rank"
                            and col != "deposition avg"
                        ):
                            for out_p in out_params:
                                if out_coord == "Output_Parameter_ID":
                                    non_features_list.append(
                                        [
                                            f"{reduce_name}d{param_id_map[out_p]}/{col}",
                                            out_p,
                                            col,
                                        ]
                                    )
                                else:
                                    non_features_list.append(
                                        [f"{reduce_name}d{out_p}/{col}", out_p, col]
                                    )
                        else:
                            non_features_list.append([f"{reduce_name}{col}", None, col])
    else:
        rank_set = False
        if "rank" in data.name or "avg" in data.name:
            rank_set = True
        data_array = data
        avg_name = f"{reduce_name}" + data.name
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
                out_p = name.split("/")[0][len(reduce_name) + 1 :]
                in_p = name.split("/")[1]
                if out_coord == "Output_Parameter_ID":
                    param_i = np.argwhere(np.asarray(param_id_map) == out_p).item()
                else:
                    param_i = out_p
                fit_data[:, i] = np.reshape(
                    data[in_p].sel({out_coord: param_i}).values, n_samples
                )
            else:
                # Model state
                state_name = name[len(reduce_name) :]
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

    def ndim_cluster(fit_data, fit_data_normed, fit_mask):
        kmeans_dic = {
            "cluster": np.asarray([]),
            "trajectory": np.asarray([]),
            "file": np.asarray([]),
        }
        clusters = KMeans(
            n_clusters=k, tol=tol, init="k-means++", random_state=42
        ).fit_predict(fit_data_normed[fit_mask, :])
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
            kmeans_dic[name] = fit_data[fit_mask, i].flatten()
        if include_all_data:
            for non_feat in non_features_list:
                if non_feat[1] is None:
                    tmp_arr = np.asarray(data[non_feat[2]].values).flatten()
                else:
                    tmp_arr = np.asarray(
                        data.sel({out_coord: non_feat[1]})[non_feat[2]].values
                    ).flatten()

                kmeans_dic[non_feat[0]] = tmp_arr[fit_mask]
        return pd.DataFrame.from_dict(kmeans_dic)

    if param_names is not None:
        fit_data, fit_data_normed, fit_mask = get_fit_data()
        return ndim_cluster(fit_data, fit_data_normed, fit_mask)

    elif n_features == 1:
        fit_data = np.reshape(data_array.values, shape)
        scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
        fit_data_normed = scaler.fit_transform(fit_data)
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
        kmeans_dic = {
            "cluster": clusters,
            "trajectory": trajectory[~np.isnan(fit_data[:, 0])],
            avg_name: fit_data[~np.isnan(fit_data[:, 0])].flatten(),
            "file": filenames[~np.isnan(fit_data[:, 0])],
        }

        if include_all_data:
            for non_feat in non_features_list:
                if non_feat[1] is None:
                    tmp_arr = np.asarray(data[non_feat[2]].values).flatten()
                else:
                    tmp_arr = np.asarray(
                        data.sel({out_coord: non_feat[1]})[non_feat[2]].values
                    ).flatten()
                kmeans_dic[non_feat[0]] = tmp_arr[~np.isnan(fit_data[:, 0])]
        return pd.DataFrame.from_dict(kmeans_dic)
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


def plot_cluster_data_interactive(data, reduce_name=""):
    """
    Calling this function from a Jupyter notebook allows to visualize the cluster association with different
    dimensions. Make sure to call pn.extension() from your notebook first.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame generated using get_cluster().
    reduce_name : string
        Name prepended to the columns. Should relate to the reduction
        applied to the dataset, such as "avg" or "rank" in get_cluster().
    Returns
    -------

    """
    out_params = []
    in_params = []
    for col in data:
        if "/" in col:
            out_params.append(col.split("/")[0][len(reduce_name) + 1 :])
            in_params.append(col.split("/")[1])
        elif reduce_name in col and len(reduce_name) > 0:
            in_params.append(col[len(reduce_name) :])
        else:
            in_params.append(col)
    in_params = list(set(in_params))
    in_params.sort()
    out_params = list(set(out_params))
    out_params.sort()
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
        end=5,
        step=0.1,
        value=0.7,
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    dot_slider = pn.widgets.IntSlider(
        name="Change the dot size",
        start=1,
        end=200,
        step=2,
        value=12,
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
        save_path,
        latex,
        save,
        s,
    ):

        sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
        fig = Figure()
        ax = fig.subplots()

        if in_p_x[0] == "d" and in_p_x != "deposition":
            if " " in in_p_x:
                in_p_x_parse = in_p_x.split()[0]
                reduce_name_x = in_p_x.split()[1] + " "
            else:
                reduce_name_x = reduce_name
                in_p_x_parse = in_p_x
            x = f"{reduce_name}d{out_p_x}/{in_p_x}"
            xlabel = (
                reduce_name_x + r"$\partial$" + out_p_x + f"/{parse_word(in_p_x_parse)}"
            )
        elif in_p_x in latex_mappings:
            x = f"{reduce_name}{in_p_x}"
            xlabel = f"{reduce_name}{parse_word(in_p_x)}"
        else:
            x = in_p_x
            xlabel = reduce_name + in_p_x

        if in_p_y[0] == "d" and in_p_y != "deposition":
            if " " in in_p_y:
                in_p_y_parse = in_p_y.split()[0]
                reduce_name_y = in_p_y.split()[1] + " "
            else:
                reduce_name_y = reduce_name
                in_p_y_parse = in_p_y
            y = f"{reduce_name}d{out_p_y}/{in_p_y}"
            ylabel = (
                reduce_name_y
                + r"$\partial$"
                + parse_word(out_p_y)
                + f"/{parse_word(in_p_y_parse)}"
            )
        elif in_p_y in latex_mappings:
            y = f"{reduce_name}{in_p_y}"
            ylabel = f"{reduce_name}{parse_word(in_p_y)}"
        else:
            y = in_p_y
            ylabel = reduce_name + in_p_y

        histogram = (in_p_x == in_p_y and out_p_x == out_p_y) or (
            in_p_x == in_p_y and "/" not in x
        )
        if histogram:
            palette = "tab10"
            if len(set(data["cluster"])) > 10:
                palette = "tab20"
            sns.histplot(
                data=data,
                x=x,
                hue="cluster",
                palette=palette,
                multiple="stack",
                bins=100,
                log_scale=(logx, logy),
                ax=ax,
            )
        else:
            sns.scatterplot(
                data=data,
                x=x,
                y=y,
                hue="cluster",
                palette="tab10",
                s=s,
                ax=ax,
            )
        if logx and not histogram:
            if np.nanmin(data[x]) < 0:
                linthresh = np.nanmin(np.abs(data[x].where(data[x] != 0)))
                ax.set_xscale("symlog", linthresh=linthresh)
            else:
                ax.set_xscale("log")
        if logy and not histogram:
            if np.nanmin(data[y]) < 0:
                linthresh = np.nanmin(np.abs(data[y].where(data[y] != 0)))
                ax.set_yscale("symlog", linthresh=linthresh)
            else:
                ax.set_yscale("log")
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=int(10 * font_scale),
        )
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
        ax.set_xlabel(xlabel, fontsize=int(11 * font_scale))
        ax.set_ylabel(ylabel, fontsize=int(11 * font_scale))
        legend = ax.get_legend()
        legend.set_title("cluster", prop={"size": int(11 * font_scale)})
        plt.setp(legend.get_texts(), fontsize=int(10 * font_scale))
        ax.yaxis.get_offset_text().set_fontsize(int(11 * font_scale))
        ax.xaxis.get_offset_text().set_fontsize(int(11 * font_scale))
        # You may use the following line to remove the offset label if needed.
        # ax.xaxis.get_offset_text().set(alpha=0)
        if save:
            try:
                i = 0
                store_type = save_path.split(".")[-1]
                store_path = save_path[0 : -len(store_type) - 1]
                save_name = store_path + "_{:03d}.".format(i) + store_type
                while os.path.isfile(save_name):
                    i = i + 1
                    save_name = store_path + "_{:03d}.".format(i) + store_type
                ax.figure.savefig(save_name, bbox_inches="tight", dpi=300)
            except:
                save_to_field.value = (
                    f"Could not save to {save_path}. Did you forget the filetype?"
                )
                pass
            save_path = None
            save = False
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
            save_path=save_to_field,
            latex=latex_button,
            save=save_button,
            s=dot_slider,
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
        pn.Row(
            dot_slider,
            latex_button,
        ),
        pn.Row(
            save_to_field,
            save_button,
        ),
        title_widget,
        plot_pane,
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
    min_time = None
    max_time = None
    out_coord = "Output_Parameter_ID"
    if verbose:
        print("Load the files to check for the time ranges for the final file")
    for center_df in tqdm(centers) if verbose else centers:
        for idx, row in (
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
    if verbose:
        print(
            "For each center, get all the trajectories and merge them into one dataset."
        )
    for center_df in tqdm(centers) if verbose else centers:
        data_arrays = []
        cluster = np.unique(center_df["cluster"])[0]
        for idx, row in (
            tqdm(center_df.iterrows(), leave=False, total=len(center_df))
            if verbose
            else center_df.iterrows()
        ):
            traj = row["trajectory"]
            f = row["file"]
            ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
            # print(ds["time_after_ascent"].attrs)
            # print(ds["ensemble"].attrs)
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
            if "asc600" in ds:
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
            if "asc600" not in ds:
                ds["ensemble"] = [0]
            ds["ensemble"].attrs = ens_attrs
            ds["time"].attrs = time_attrs
            if fix_weird_coords:
                weird_lon = np.count_nonzero(ds["lon"] == 0)
                weird_lat = np.count_nonzero(ds["lat"] == 0)
                # A rare bug can happen where a single lat or lon is set to zero inbetween valid values
                if weird_lon > 0:
                    idxs = np.argwhere(ds["lon"].where(ds["lon"] == 0).values == 0)[
                        :, 2
                    ]
                    # We assume the idxs are consecutive
                    n_fills = len(idxs)
                    lon_tmp = ds["lon"].isel(
                        {
                            "time": [idxs[0] - 1, idxs[-1] + 1],
                        }
                    )
                    val_delta = (
                        (lon_tmp.isel({"time": 1}) - lon_tmp.isel({"time": 0}))
                        / n_fills
                    ).values.item()
                    val_start = lon_tmp.isel({"time": 1}).values.item()
                    for i_idx, idx in enumerate(idxs):
                        ds["lon"][0, 0, idx] = val_start + (i_idx + 1) * val_delta
                if weird_lat > 0:
                    idxs = np.argwhere(ds["lat"].where(ds["lat"] == 0).values == 0)[
                        :, 2
                    ]
                    # We assume the idxs are consecutive
                    n_fills = len(idxs)
                    lat_tmp = ds["lat"].isel(
                        {
                            "time": [idxs[0] - 1, idxs[-1] + 1],
                        }
                    )
                    val_delta = (
                        (lat_tmp.isel({"time": 1}) - lat_tmp.isel({"time": 0}))
                        / n_fills
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
                if out_coord in ds:
                    sens_shape = (len(ds[out_coord]), 1, 1, times)
                    add_sens = np.reshape(
                        np.repeat(np.NaN, times * len(ds[out_coord])),
                        sens_shape,
                    )
                val_shape = (1, 1, times)
                add_vals = np.reshape(np.repeat(np.NaN, times), val_shape)
                for var in ds:
                    if var[0] != "d" or var == "deposition" or var == "dp2h":
                        variables[var] = (["ensemble", "trajectory", "time"], add_vals)
                    else:
                        variables[var] = (
                            [out_coord, "ensemble", "trajectory", "time"],
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
            elif tmp_max < max_time:
                variables = {}
                times = int(max_time + 20 - tmp_max) // 30
                if out_coord in ds:
                    sens_shape = (len(ds[out_coord]), 1, 1, times)
                    add_sens = np.reshape(
                        np.repeat(np.NaN, times * len(ds[out_coord])),
                        sens_shape,
                    )
                val_shape = (1, 1, times)
                add_vals = np.reshape(np.repeat(np.NaN, times), val_shape)
                for var in ds:
                    if var[0] != "d" or var == "deposition" or var == "dp2h":
                        variables[var] = (["ensemble", "trajectory", "time"], add_vals)
                    else:
                        variables[var] = (
                            [out_coord, "ensemble", "trajectory", "time"],
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
    if keep_coords == "relative":
        cols = []
        for col in ds:
            if (col != "lon") and (col != "lat"):
                cols.append(col)
        ds = ds[cols]
    elif keep_coords == "normal":
        cols = []
        for col in ds:
            if (col != "relative_lon") and (col != "relative_lat"):
                cols.append(col)
        ds = ds[cols]
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
        reduce_name="avg ",
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
