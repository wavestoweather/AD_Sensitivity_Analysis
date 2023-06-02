"""Get top parameters w.r.t. different attributes.

"""
import os

import numpy as np
from tqdm.auto import tqdm
import xarray as xr

from ad_sensitivity_analysis.data_handler.loader import load_dataset_part
from ad_sensitivity_analysis.plot.latexify import param_id_map


def get_top_list(ds, print_out=True, verbose=True):
    """

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    print_out : bool
        Print the number of parameters and the dataset.
    verbose : bool
        Print the top parameters for each output variable.

    Returns
    -------
    list of output parameters, list of top20 parameters, list of top10 parameters, dict of top20 parameters,
    dict of top10 parameters.
    """
    out_params = [
        "QV",
        "QC",
        "QR",
        "QG",
        "QH",
        "QI",
        "QS",
        "NCCLOUD",
        "NCRAIN",
        "NCGRAUPEL",
        "NCHAIL",
        "NCICE",
        "NCSNOW",
        "QR_OUT",
        "QG_OUT",
        "QH_OUT",
        "QI_OUT",
        "QS_OUT",
        "NR_OUT",
        "NG_OUT",
        "NH_OUT",
        "NI_OUT",
        "NS_OUT",
    ]
    top20_sens_dic = {}
    top10_sens_dic = {}
    tmp20 = []
    tmp10 = []
    if print_out:
        print("\nGet the top parameters for each output variable\n")
    for out_p in out_params:
        if out_p not in ds["Output Parameter"]:
            continue
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        # Skip those that did not appear
        if np.max(df["Predicted Squared Error"]) == 0:
            continue
        top20_sens_dic[out_p] = list(
            np.unique(df.nlargest(20, "Predicted Squared Error")["Input Parameter"])
        )
        top10_sens_dic[out_p] = list(
            np.unique(df.nlargest(10, "Predicted Squared Error")["Input Parameter"])
        )
        tmp20.extend(top20_sens_dic[out_p])
        tmp10.extend(top10_sens_dic[out_p])
        if verbose:
            print(f"###################{out_p}")
            print(f"Top 10: \n{top10_sens_dic[out_p]}")
            print(f"Top 20: \n{top20_sens_dic[out_p]}")
    top20_list = list(set(tmp20))
    top10_list = list(set(tmp10))
    if print_out or verbose:
        print(
            f"Number of distinct parameters by taking the top 20 for everything: {len(top20_list)}"
        )
        print(
            f"Number of distinct parameters by taking the top 10 for everything: {len(top10_list)}"
        )
        print(ds)
    return (
        list(ds["Output Parameter"].values),
        top20_list,
        top10_list,
        top20_sens_dic,
        top10_sens_dic,
    )


def get_magnitude_list(ds, out_params, print_out=True, verbose=True):
    """
    Get the top parameters within one, two, and three orders of magnitude.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    out_params : list-like of strings
        The model state variables for which sensitivities have been calculated for.
    print_out : bool
        Print the top parameters within one order of magnitude and the number of parameters in each magnitude list.
    verbose : bool
        Print the top parameters for each output variable for each magnitude.

    Returns
    -------
    list of parameters within one order of magnitude, list within two orders of magnitudes,
    list within three orders of magnitudes
    """
    if print_out:
        print(
            "\nInstead of taking the top 10, take the parameters in the same order of magnitude\n"
        )
    top_one_order = []
    top_two_orders = []
    top_three_orders = []

    for out_p in out_params:
        if out_p not in ds["Output Parameter"]:
            continue
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        # Skip those that did not appear
        if np.max(df["Predicted Squared Error"]) == 0:
            continue
        max_order = np.max(df["Predicted Squared Error"])
        top_one_order_tmp = list(
            np.unique(
                df[df["Predicted Squared Error"] >= max_order / 10]["Input Parameter"]
            )
        )
        top_two_orders_tmp = list(
            np.unique(
                df[df["Predicted Squared Error"] >= max_order / 100]["Input Parameter"]
            )
        )
        top_three_orders_tmp = list(
            np.unique(
                df[df["Predicted Squared Error"] >= max_order / 1000]["Input Parameter"]
            )
        )
        top_one_order.extend(top_one_order_tmp)
        top_two_orders.extend(top_two_orders_tmp)
        top_three_orders.extend(top_three_orders_tmp)
        if verbose:
            print(f"###################{out_p}")
            print(f"Top order: \n{top_one_order_tmp}")
            print(f"Top 2 orders: \n{top_two_orders_tmp}")
            print(f"Top 3 orders: \n{top_three_orders_tmp}")
    top_three_orders_list = list(set(top_three_orders))
    top_two_orders_list = list(set(top_two_orders))
    top_one_order_list = list(set(top_one_order))
    if print_out or verbose:
        print(
            f"Number of distinct parameters by taking the top order of magnitude: {len(top_one_order_list)}"
        )
        print(
            f"Number of distinct parameters by taking the top 2 orders of magnitude: {len(top_two_orders_list)}"
        )
        print(
            f"Number of distinct parameters by taking the top 3 orders of magnitude: {len(top_three_orders_list)}"
        )
        print("Parameters within the top order of magnitude:")
        print(top_one_order_list)
    return top_one_order_list, top_two_orders_list, top_three_orders_list


def traj_get_sum_derivatives(file_path):
    """
    Calculate the sum of absolute values of sensitivities for each model state variable and model parameter.

    Parameters
    ----------
    file_path : String
        Path to NetCDF-files with sensitivities that used simulation_mode 1 (sensitivity analysis for trajectories)
        in the AD-based C++ simulation.

    Returns
    -------
    Dictionary with tracked output parameters as keys and a pandas.Dataframe with the sum of
    absolute values of the gradients, and a list of strings with the output parameter names.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    out_param_coord = "Output_Parameter_ID"
    if out_param_coord not in ds:
        out_param_coord = "Output Parameter"
        out_params = list(ds[out_param_coord].values)
        param_name = out_params
    else:
        param_name = []
        out_params = list(ds[out_param_coord].values)
        for idx in out_params:
            param_name.append(param_id_map[idx.values])

    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]

    sums = {}
    for f in tqdm(files):
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
        for out_p, out_name in zip(out_params, param_name):
            ds[in_params] = np.abs(ds[in_params])
            df = (
                ds[in_params]
                .sel({out_param_coord: out_p})
                .sum(dim=["trajectory", "time"], skipna=True)
                .to_dataframe()
                .reset_index()
            )
            df = df[in_params]

            if out_name in sums:
                sums[out_name] += df
            else:
                sums[out_name] = df
    return sums, param_name


def traj_get_top_params(dict_of_df, param_name, n, orders):
    """
    Given the results of get_sums() or get_sums_phase(), get the parameters with the highest sensitivities,
    within a given order of magnitude or within the top n parameters.

    Parameters
    ----------
    dict_of_df : dict of pandas.DataFrame
        Dictionary with (phase +) model state variables as keys and the sums for each gradients as pandas.DataFrame.
    param_name : list-like of strings
        Keys of dict_of_df.
    n : int
        Get the top n parameters for each tracked model state variable
    orders : int or float
        Get the parameters within orders many orders of magnitude for each tracked model state variable

    Returns
    -------
    Set of parameters within the given order of magnitude, set of top n parameters, dictionaries for magnitude and top
    n parameters for each output parameter.

    """
    top_sens_dic = {}
    top_magn_sens_dic = {}

    for out_name in param_name:
        tmp_df = dict_of_df[out_name].T
        max_order = np.max(tmp_df[0])
        if max_order == 0:
            continue
        top_sens_dic[out_name] = dict_of_df[out_name].T.nlargest(n, 0).T.columns.values
        top_magn_sens_dic[out_name] = tmp_df[
            tmp_df[0] >= max_order / (10**orders)
        ].T.columns.values

    tmp = []
    tmp2 = []
    for out_name in param_name:
        if out_name in top_magn_sens_dic:
            tmp.extend(top_magn_sens_dic[out_name])
            tmp2.extend(top_sens_dic[out_name])
    top_magn_set = set(tmp)
    top10_set = set(tmp2)
    return top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic


def get_edges(
    min_max,
    min_max_in_params,
    in_params,
    param_name,
    additional_params,
    n_bins,
    log=False,
    verbose=False,
):
    """
    Create edges for 1D or 2D histograms.

    Parameters
    ----------
    min_max : Dict of lists with two floats
        Keys are model state variables, values are [minimum, maximum] floats.
    min_max_in_params : Dict of dict of lists with two floats
        First keys are model state variablees, second keys are model parameters. The values are [minimum, maximum]
        floats.
    in_params : List of strings
        Names of model parameters.
    param_name : List of strings
        Name of model state variables for which sensitivities are available.
    additional_params : List of strings
        List of additional model state variables to get edges for.
    n_bins : int
        Number of bins.
    log : bool
        Use log-scale for the edges.
    verbose : bool
        More progressbars.

    Returns
    -------
    Dict of lists of floats (edges for model state variables), Dict of dict of lists of floats (edges
    for model parameters for each model state variable defined in param_name).
    """

    edges = {}
    edges_in_params = {}
    for out_p in param_name:
        # -2 because we add one bin each to the left and right to catch any
        # floats just on the edge.
        delta = (min_max[out_p][1] - min_max[out_p][0]) / (n_bins - 2)
        edges[out_p] = np.arange(
            min_max[out_p][0] - delta, min_max[out_p][1] + 0.5 * delta, delta
        )
        edges_in_params[out_p] = {}
        for in_p in in_params:
            if min_max_in_params[out_p][in_p][0] == min_max_in_params[out_p][in_p][
                1
            ] or np.isnan(min_max_in_params[out_p][in_p][0]):
                continue
            delta = (
                min_max_in_params[out_p][in_p][1] - min_max_in_params[out_p][in_p][0]
            ) / (n_bins - 2)
            if log:
                start = min_max_in_params[out_p][in_p][0] - delta
                stop = min_max_in_params[out_p][in_p][1]
                offset = 0 if start >= 0 else start * 2

                edges_in_params[out_p][in_p] = np.logspace(
                    np.log10(np.abs(start)),
                    np.log10(np.abs(stop - offset)),
                    n_bins,
                )
                if start < 0:
                    edges_in_params[out_p][in_p] += offset
            else:
                edges_in_params[out_p][in_p] = np.arange(
                    min_max_in_params[out_p][in_p][0] - delta,
                    min_max_in_params[out_p][in_p][1] + 0.5 * delta,
                    delta,
                )
    if additional_params is not None:
        for out_p in (
            tqdm(additional_params, leave=False) if verbose else additional_params
        ):
            if min_max[out_p][0] == min_max[out_p][1] or np.isnan(min_max[out_p][0]):
                continue
            delta = (min_max[out_p][1] - min_max[out_p][0]) / (n_bins - 2)
            if log:
                start = min_max[out_p][0] - delta
                stop = min_max[out_p][1]
                offset = 0 if start >= 0 else start * 2
                edges[out_p] = np.logspace(
                    np.log10(np.abs(start)),
                    np.log10(np.abs(stop - offset)),
                    n_bins,
                )
                if start < 0:
                    edges[out_p] += offset
            else:
                edges[out_p] = np.arange(
                    min_max[out_p][0] - delta, min_max[out_p][1] + 0.5 * delta, delta
                )
    return edges, edges_in_params


def _get_param_names(files, file_path, out_params=None, in_params=None):
    ds = None
    if in_params is None:
        ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
        in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]

    if out_params is None:
        if ds is None:
            ds = xr.open_dataset(
                file_path + files[0], decode_times=False, engine="netcdf4"
            )
        out_params = ds["Output_Parameter_ID"]

        param_name = []
        for idx in out_params:
            param_name.append(param_id_map[idx.values])
    else:
        param_name = []
        tmp = []
        for idx in out_params:
            if isinstance(idx, int):
                param_name.append(param_id_map[idx])
            else:
                param_name.append(idx)
                tmp.append(np.argwhere(np.asarray(param_id_map) == idx).item())
        for i, val in enumerate(tmp):
            out_params[i] = val
    return param_name, out_params, in_params


def _set_min_max(
    ds_tmp,
    out_name,
    min_max,
    min_max_in_params,
    in_params,
    filter_mag=None,
    means=None,
    verbose=False,
):
    """

    Parameters
    ----------
    ds_tmp
    out_name
    min_max
    min_max_in_params
    in_params
    filter_mag
    means
    verbose

    Returns
    -------

    """
    min_p = np.nanmin(ds_tmp[out_name])
    max_p = np.nanmax(ds_tmp[out_name])
    if out_name in min_max:
        if min_p < min_max[out_name][0]:
            min_max[out_name][0] = min_p
        if max_p > min_max[out_name][1]:
            min_max[out_name][1] = max_p
    else:
        min_max[out_name] = [min_p, max_p]
    for in_p in tqdm(in_params, leave=False) if verbose else in_params:
        if filter_mag is not None and means is not None:
            mean = means[out_name][in_p]
            filtered_ds = ds_tmp[in_p]
            filtered_ds = filtered_ds.where(
                np.abs(np.log10(np.abs(filtered_ds / mean))) <= filter_mag
            )
            min_p = np.nanmin(filtered_ds)
            max_p = np.nanmax(filtered_ds)
        else:
            min_p = np.nanmin(ds_tmp[in_p])
            max_p = np.nanmax(ds_tmp[in_p])
        if in_p in min_max_in_params[out_name]:
            if min_p < min_max_in_params[out_name][in_p][0]:
                min_max_in_params[out_name][in_p][0] = min_p
            if max_p > min_max_in_params[out_name][in_p][1]:
                min_max_in_params[out_name][in_p][1] = max_p
        else:
            min_max_in_params[out_name][in_p] = [min_p, max_p]


def _set_hist_in_params(
    ds_tmp,
    out_name,
    hist,
    hist_in_params,
    in_params,
    edges,
    edges_in_params,
    verbose=False,
):
    """

    Parameters
    ----------
    ds_tmp
    out_name
    hist
    hist_in_params
    in_params
    edges
    edges_in_params
    verbose

    Returns
    -------

    """
    hist_tmp, _ = np.histogram(ds_tmp[out_name], edges[out_name])
    if out_name in hist:
        hist[out_name] += hist_tmp
    else:
        hist[out_name] = hist_tmp
        hist_in_params[out_name] = {}
    for in_p in tqdm(in_params, leave=False) if verbose else in_params:
        if in_p not in edges_in_params[out_name]:
            continue
        hist_tmp, _ = np.histogram(ds_tmp[in_p], edges_in_params[out_name][in_p])
        if in_p in hist_in_params[out_name]:
            hist_in_params[out_name][in_p] += hist_tmp
        else:
            hist_in_params[out_name][in_p] = hist_tmp


# pylint: disable=too-many-arguments, too-many-locals, too-many-branches
def get_histogram(
    file_path,
    in_params=None,
    out_params=None,
    n_bins=100,
    additional_params=None,
    only_asc600=False,
    only_phase=None,
    inoutflow_time=-1,
    filter_mag=None,
    means=None,
    log=False,
    verbose=False,
):
    """
    Calculate 1D histograms.

    Parameters
    ----------
    file_path : string
        Path to a folder with many files from a sensitivity analysis simulation.
    in_params : list-like of strings
        List of model parameters.
    out_params : list of strings or ints
        Define for which output parameters the sensitivities should be considered for.
        Can be either ids ("Output_Parameter_ID") or the names of the parameters, e.g., "QV".
    n_bins : int
        Number of bins.
    additional_params : list-like of strings
        Additional parameters that are not model parameters and for which no sensitivity is available but you
        wish to have a histogram for.
    only_asc600 : bool
        Consider only time steps during the ascend.
    only_phase : string
        Consider only time steps with the given phase. Can be combined with only_asc600 or inoutflow_time.
        Possible values are "warm phase", "mixed phase", "ice phase", "neutral phase".
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.
    filter_mag : float
        Filter all values that are more than the given magnitude larger or smaller
        than the mean.
    means : Dictionary of floats
        Mean values for the parameters. Model parameters are dictionaries
        (model state for which the sensitivitity is for) of dictionaries (the model parameter).
        Only needed for filtering data.
    log : bool
        Use log-scale for the edges.
    verbose : bool
        Additional progressbars.

    Returns
    -------
    Dictionary with the following keys and values:
    hist_out_params: Dictionary (keys = model state variable) of arrays with values of the histogram for the given key.
    hist_in_params: Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        values of the histogram for the given key.
    edges_out_params: Dictionary (keys = model state variable) of arrays with the bin edges for the given key.
    edges_in_params: Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        the bin edges for the given keys.
    """
    phases = np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
    try:
        if only_phase is not None and only_phase not in phases:
            raise KeyError
    except KeyError:
        print(
            f"You asked for phase {only_phase}, which does not exist."
            f"Possible phases are {phases}"
        )
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    param_name, out_params, in_params = _get_param_names(
        files=files,
        file_path=file_path,
        out_params=out_params,
        in_params=in_params,
    )

    min_max = {}
    min_max_in_params = {}
    for out_p in param_name:
        min_max_in_params[out_p] = {}
    load_params = param_name.copy()
    if in_params is not None:
        load_params.extend(in_params)
    if additional_params is not None:
        load_params.extend(additional_params)

    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            only_phase=only_phase,
            inoutflow_time=inoutflow_time,
            load_params=load_params,
        )
        for out_p, out_name in (
            tqdm(zip(out_params, param_name), leave=False, total=len(param_name))
            if verbose
            else zip(out_params, param_name)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            _set_min_max(
                ds_tmp=ds_tmp,
                out_name=out_name,
                min_max=min_max,
                min_max_in_params=min_max_in_params,
                in_params=in_params,
                filter_mag=filter_mag,
                means=means,
                verbose=verbose,
            )
        if additional_params is not None:
            for out_p in (
                tqdm(additional_params, leave=False) if verbose else additional_params
            ):
                min_p = np.nanmin(ds[out_p])
                max_p = np.nanmax(ds[out_p])
                if out_p in min_max:
                    if min_p < min_max[out_p][0]:
                        min_max[out_p][0] = min_p
                    if max_p > min_max[out_p][1]:
                        min_max[out_p][1] = max_p
                else:
                    min_max[out_p] = [min_p, max_p]

    edges, edges_in_params = get_edges(
        min_max=min_max,
        min_max_in_params=min_max_in_params,
        in_params=in_params,
        param_name=param_name,
        additional_params=additional_params,
        n_bins=n_bins,
        log=log,
        verbose=verbose,
    )

    hist = {}
    hist_in_params = {}
    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            only_phase=only_phase,
            inoutflow_time=inoutflow_time,
            load_params=load_params,
        )
        for out_p, out_name in (
            tqdm(zip(out_params, param_name), leave=False, total=len(param_name))
            if verbose
            else zip(out_params, param_name)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            _set_hist_in_params(
                ds_tmp=ds_tmp,
                out_name=out_name,
                hist=hist,
                hist_in_params=hist_in_params,
                in_params=in_params,
                edges=edges,
                edges_in_params=edges_in_params,
                verbose=verbose,
            )
        if additional_params is not None:
            for out_p in (
                tqdm(additional_params, leave=False) if verbose else additional_params
            ):
                if out_p not in edges:
                    # In case no values for the additional parameters are available, i.e.,
                    # because only_phase is warm phase but QI is given in additonal_params.
                    continue
                hist_tmp, _ = np.histogram(ds[out_p], edges[out_p])
                if out_p in hist:
                    hist[out_p] += hist_tmp
                else:
                    hist[out_p] = hist_tmp

    return {
        "hist_out_params": hist,
        "hist_in_params": hist_in_params,
        "edges_out_params": edges,
        "edges_in_params": edges_in_params,
    }


def _set_hist_in_params_cond(
    ds_tmp,
    out_name,
    hist_conditional,
    cond,
    in_params,
    edges,
    edges_in_params,
    verbose=False,
):
    """

    Parameters
    ----------
    ds_tmp
    out_name
    hist_conditional
    cond
    in_params
    edges
    edges_in_params
    verbose

    Returns
    -------

    """
    hist_tmp, _, _ = np.histogram2d(
        ds_tmp[cond].values.flatten(),
        ds_tmp[out_name].values.flatten(),
        [edges[cond], edges[out_name]],
    )
    if out_name in hist_conditional[cond]["hist_out_params"]:
        hist_conditional[cond]["hist_out_params"][out_name] += hist_tmp
    else:
        hist_conditional[cond]["hist_out_params"][out_name] = hist_tmp
        hist_conditional[cond]["hist_in_params"][out_name] = {}
    for in_p in tqdm(in_params, leave=False) if verbose else in_params:
        if in_p not in edges_in_params[out_name]:
            continue
        hist_tmp, _, _ = np.histogram2d(
            ds_tmp[cond].values.flatten(),
            ds_tmp[in_p].values.flatten(),
            [edges[cond], edges_in_params[out_name][in_p]],
        )
        if in_p in hist_conditional[cond]["hist_in_params"][out_name]:
            hist_conditional[cond]["hist_in_params"][out_name][in_p] += hist_tmp
        else:
            hist_conditional[cond]["hist_in_params"][out_name][in_p] = hist_tmp


def _set_hist_out_params_cond(
    ds, additional_params, cond, edges, hist_conditional, verbose=False
):
    """

    Parameters
    ----------
    ds
    additional_params
    cond
    edges
    hist_conditional
    verbose

    Returns
    -------

    """
    for out_p in tqdm(additional_params, leave=False) if verbose else additional_params:
        if out_p == cond:
            continue
        if out_p not in edges:
            # In case no values for the additional parameters are available, i.e.,
            # because only_phase is warm phase but QI is given in additonal_params.
            continue
        hist_tmp, _, _ = np.histogram2d(
            ds[cond].values.flatten(),
            ds[out_p].values.flatten(),
            [edges[cond], edges[out_p]],
        )
        if out_p in hist_conditional[cond]["hist_out_params"]:
            hist_conditional[cond]["hist_out_params"][out_p] += hist_tmp
        else:
            hist_conditional[cond]["hist_out_params"][out_p] = hist_tmp


# pylint: disable=too-complex, too-many-arguments, too-many-locals, too-many-branches
def get_histogram_cond(
    file_path,
    in_params=None,
    out_params=None,
    conditional_hist=None,
    n_bins=100,
    additional_params=None,
    only_asc600=False,
    only_phase=None,
    inoutflow_time=-1,
    filter_mag=None,
    means=None,
    log=False,
    verbose=False,
):
    """
    Get a 2D histogram where 'conditional_hist' is the lsit of parameters for which the edges are calculated.
    The final 2D histogram is the count for values that are in the bin defined by 'conditional_hist'
    and in the bin calculated from in_params or out_params.

    Parameters
    ----------
    file_path : string
        Path to a folder with many files from a sensitivity analysis simulation.
    in_params : list-like of strings
        List of model parameters.
    out_params : list of strings or ints
        Define for which output parameters the sensitivities should be considered for.
        Can be either ids ("Output_Parameter_ID") or the names of the parameters, e.g., "QV".
    conditional_hist : string
        The model state variable (aka output parameter) for which additional edges shall be calculated.
    n_bins : int
        Number of bins.
    additional_params : list-like of strings
        Additional parameters that are not model parameters and for which no sensitivity is available but you
        wish to have a histogram for.
    only_asc600 : bool
        Consider only time steps during the ascend.
    only_phase : string
        Consider only time steps with the given phase. Can be combined with only_asc600 or inoutflow_time.
        Possible values are "warm phase", "mixed phase", "ice phase", "neutral phase".
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.
    filter_mag : float
        Filter all values that are more than the given magnitude larger or smaller
        than the mean.
    means : Dictionary of floats
        Mean values for the parameters. Model parameters are dictionaries
        (model state for which the sensitivitity is for) of dictionaries (the model parameter).
        Only needed for filtering data.
    log : bool
        Use log-scale for the edges.
    verbose : bool
        Additional progressbars.

    Returns
    -------
    Dictionary with the following keys:
    'edges_out_params': Dictionary where the keys are model state variables and the values are arrays of bin edges.
    'edges_in_params': Dictionary where the keys are model state variables for which sensitivities are available
        and the values are dictionaries of model parameters with arrays of bin edges.
    model state variables: Each model state variable has a dictionary for 'hist_out_params' and 'hist_in_params'.
    'hist_out_params' is a dictionary of model state variables with arrays of bin counts.
    'hist_in_params' is a dictionary of model state variables for which sensitivities are available
        and the values are dictionaries of model parameters with arrays of bin counts.
    """
    if conditional_hist is None:
        conditional_hist = []
    if additional_params is None:
        additional_params = []
    phases = np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
    if only_phase is not None and only_phase not in phases:
        raise (
            f"You asked for phase {only_phase}, which does not exist."
            f"Possible phases are {phases}"
        )
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    param_name, out_params, in_params = _get_param_names(
        files=files,
        file_path=file_path,
        out_params=out_params,
        in_params=in_params,
    )

    for cond in conditional_hist:
        if cond not in additional_params:
            additional_params.append(cond)

    min_max = {}
    min_max_in_params = {}
    for out_p in param_name:
        min_max_in_params[out_p] = {}
    load_params = param_name.copy()
    if in_params is not None:
        load_params.extend(in_params)
    if additional_params is not None:
        load_params.extend(additional_params)

    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            only_phase=only_phase,
            inoutflow_time=inoutflow_time,
            load_params=load_params,
        )
        for out_p, out_name in (
            tqdm(zip(out_params, param_name), leave=False, total=len(param_name))
            if verbose
            else zip(out_params, param_name)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            _set_min_max(
                ds_tmp=ds_tmp,
                out_name=out_name,
                min_max=min_max,
                min_max_in_params=min_max_in_params,
                in_params=in_params,
                filter_mag=filter_mag,
                means=means,
                verbose=verbose,
            )
        for out_p, min_max_vals in (
            tqdm(additional_params.items(), leave=False)
            if verbose
            else additional_params.items()
        ):
            min_p = np.nanmin(ds[out_p])
            max_p = np.nanmax(ds[out_p])
            if out_p in min_max:
                if min_p < min_max_vals[0]:
                    additional_params[out_p][0] = min_p
                if max_p > min_max_vals[1]:
                    additional_params[out_p][1] = max_p
            else:
                additional_params[out_p] = [min_p, max_p]

    edges, edges_in_params = get_edges(
        min_max=min_max,
        min_max_in_params=min_max_in_params,
        in_params=in_params,
        param_name=param_name,
        additional_params=additional_params,
        n_bins=n_bins,
        log=log,
        verbose=verbose,
    )

    hist_conditional = {
        "edges_out_params": edges,
        "edges_in_params": edges_in_params,
    }

    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            only_asc600=only_asc600,
            only_phase=only_phase,
            inoutflow_time=inoutflow_time,
            load_params=load_params,
        )
        for cond in (
            tqdm(conditional_hist, leave=False) if verbose else conditional_hist
        ):
            if cond not in edges:
                continue
            if cond not in hist_conditional:
                hist_conditional[cond] = {
                    "hist_out_params": {},
                    "hist_in_params": {},
                }
            for out_p, out_name in (
                tqdm(zip(out_params, param_name), leave=False, total=len(param_name))
                if verbose
                else zip(out_params, param_name)
            ):
                if out_name not in edges:
                    continue
                ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
                _set_hist_in_params_cond(
                    ds_tmp=ds_tmp,
                    out_name=out_name,
                    hist_conditional=hist_conditional,
                    cond=cond,
                    in_params=in_params,
                    edges=edges,
                    edges_in_params=edges_in_params,
                    verbose=verbose,
                )
            _set_hist_out_params_cond(
                ds=ds,
                additional_params=additional_params,
                cond=cond,
                edges=edges,
                hist_conditional=hist_conditional,
                verbose=verbose,
            )
    return hist_conditional
