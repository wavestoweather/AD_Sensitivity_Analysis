"""Helper functions to calculate clusters on trajectories.

"""
from ad_sensitivity_analysis.plot.latexify import param_id_map


def get_non_features_list(data_trajs, out_coord, out_params, reduce_name, x):
    """

    Parameters
    ----------
    data_trajs
    out_coord
    out_params
    reduce_name
    x

    Returns
    -------

    """
    non_features_list = []
    for col in data_trajs:
        if col not in x:
            if (
                col[0] == "d"
                and col != "deposition"
                and col != "deposition rank"
                and col != "deposition avg"
            ):
                for out_p in out_params:
                    if out_coord == "Output_Parameter_ID":
                        # pylint: disable=invalid-sequence-index
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
    return non_features_list


def parse_ds_single_model_state(
    data_trajs, x, out_coord, reduce_name, include_all_data
):
    """

    Parameters
    ----------
    data_trajs : xarray.Dataset or xarray.DataArray
        A dataset of trajectories from a sensitivity simulation or one column of the set
        where the time dimension has been reduced using an average metric.
        Must be a dataset with coordinates 'trajectory' and 'file', i.e., loaded using 'load_data()'.
    x : string (needed if data_trajs is xarray.Dataset) or list of string
        The model state variable or model parameter to create clusters for. If multiple variables or parameters shall
        be considered, then x must be a list of strings.
    out_coord
    reduce_name : string
        Name to be prepended to the columns. Should relate to the reduction
        applied to the dataset, such as "avg" or "rank".
    include_all_data : bool
        Include all columns in data_trajs in the returned dataframe. Otherwise only include columns used for clustering.

    Returns
    -------

    """
    non_features_list = []
    if not isinstance(x, str):
        x = x[0]
    data_array = data_trajs[x]
    avg_name = f"{reduce_name}{x}"
    if not include_all_data:
        return data_array, non_features_list, avg_name
    if out_coord in data_trajs:
        out_params = data_trajs[out_coord].values
    non_features_list = get_non_features_list(
        data_trajs=data_trajs,
        out_coord=out_coord,
        out_params=out_params,
        reduce_name=reduce_name,
        x=[x],
    )
    return data_array, non_features_list, avg_name


def parse_ds_single_param(
    data_trajs, x, out_coord, reduce_name, include_all_data, param_names
):
    """

    Parameters
    ----------
    data_trajs
    x
    out_coord
    reduce_name
    include_all_data
    param_names

    Returns
    -------

    """
    non_features_list = []
    if not isinstance(x, str):
        x = x[0]
        # A single parameter but for multiple model states
    avg_name_list = []
    n_features = len(param_names)
    for p in param_names:
        avg_name_list.append(f"{reduce_name}d{p}/{x}")
    if not include_all_data:
        return n_features, non_features_list, avg_name_list
    if out_coord in data_trajs:
        out_params = data_trajs[out_coord].values
    non_features_list = get_non_features_list(
        data_trajs=data_trajs,
        out_coord=out_coord,
        out_params=out_params,
        reduce_name=reduce_name,
        x=[x],
    )
    return n_features, non_features_list, avg_name_list


def parse_ds_model_states(data_trajs, x, out_coord, reduce_name, include_all_data):
    """
    Multiple model states and no model parameter.

    Parameters
    ----------
    data_trajs
    x
    out_coord
    reduce_name
    include_all_data

    Returns
    -------

    """
    non_features_list = []
    n_features = len(x)
    avg_name_list = [f"{reduce_name}{v}" for v in x]
    if not include_all_data:
        return n_features, non_features_list, avg_name_list
    if out_coord in data_trajs:
        out_params = data_trajs[out_coord].values
    else:
        out_params = ["model_state"]
    non_features_list = get_non_features_list(
        data_trajs=data_trajs,
        out_coord=out_coord,
        out_params=out_params,
        reduce_name=reduce_name,
        x=x,
    )
    return n_features, non_features_list, avg_name_list


def parse_ds_model_states_params(
    data_trajs, x, out_coord, reduce_name, include_all_data, param_names
):
    """
    Multiple model states and at least one model parameter
    Parameters
    ----------
    data_trajs
    x
    out_coord
    reduce_name
    include_all_data
    param_names

    Returns
    -------

    """
    avg_name_list = []
    for variable in x:
        if variable[0] == "d" and variable != "deposition":
            for p in param_names:
                avg_name_list.append(f"{reduce_name}d{p}/{variable}")
        else:
            avg_name_list.append(f"{reduce_name}{variable}")
    n_features = len(avg_name_list)
    if not include_all_data:
        return n_features, [], avg_name_list
    if out_coord in data_trajs:
        out_params = data_trajs[out_coord].values
    non_features_list = get_non_features_list(
        data_trajs=data_trajs,
        out_coord=out_coord,
        out_params=out_params,
        reduce_name=reduce_name,
        x=x,
    )
    return n_features, non_features_list, avg_name_list
