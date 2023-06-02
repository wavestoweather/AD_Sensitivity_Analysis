"""Get ordered parameters to create nicer tables.

"""
import numpy as np

from ad_sensitivity_analysis.plot.latexify import param_id_map


def _get_order(out_param, top_sens_dic, sums):
    """
    Get the parameters sorted by their impact for each output parameter.

    Parameters
    ----------
    out_param : string
        Output parameter where the parameters shall be sorted.
    top_sens_dic : Dictionary of strings
        Created using ad_sensitivity_analysis.statistics.top_params.traj_get_top_params(..)
    sums : Dictionary with pandas.DataFrame values.
        Created using get_sums_phase(..) or get_sums(..). A dictionary with phase + model state variables as keys and
        the sums for each gradients as pandas.DataFrame or a dictionary with model state variables as keys
        and the sums for each gradients as pandas.DataFrame.
    Returns
    -------

    """
    tuples = []
    for key in top_sens_dic:
        if out_param not in key:
            continue
        for param in top_sens_dic[key]:
            tuples.append([sums[key][param].values[0], param])
    tuples.sort()
    tuples = np.asarray(tuples[::-1])
    order = tuples[:, 1]
    return order


def _get_top_dic(top_sens_dic):
    """
    Get an ordered list of model parameters for each output parameter.

    Parameters
    ----------
    top_sens_dic : Dictionary of strings
        Created using ad_sensitivity_analysis.statistics.top_params.traj_get_top_params(..)

    Returns
    -------
    Dictionary with output parameters as keys and sorted lists of parameters.
    """
    top_dic = {}
    for param in param_id_map:
        for key in top_sens_dic:
            if param in key and param not in top_dic:
                top_dic[param] = []
    for key in top_sens_dic:
        for val in top_sens_dic[key]:
            for out_param, vals in top_dic.items():
                if out_param in key:
                    if val not in vals:
                        vals.append(val)
                    break
    return top_dic


def _get_top_list_ordered(top_list, order):
    ordered_list = []
    for param in order:
        if param in top_list and param not in ordered_list:
            ordered_list.append(param)
    return ordered_list


def get_ordered_top_dic(top_sens_dic, sums):
    """

    Parameters
    ----------
    top_sens_dic : Dictionary of strings
        Created using ad_sensitivity_analysis.statistics.top_params.traj_get_top_params(..)
    sums : Dictionary with pandas.DataFrame values.
        Created using get_sums_phase(..) or get_sums(..). A dictionary with phase + model state variables as keys and
        the sums for each gradients as pandas.DataFrame or a dictionary with model state variables as keys
        and the sums for each gradients as pandas.DataFrame.

    Returns
    -------
    An ordered dictionary that can be used to generate sorted tables.
    """
    top_dic = _get_top_dic(top_sens_dic)
    return {
        key: _get_top_list_ordered(top_params, _get_order(key, top_sens_dic, sums))
        for key, top_params in top_dic.items()
    }
