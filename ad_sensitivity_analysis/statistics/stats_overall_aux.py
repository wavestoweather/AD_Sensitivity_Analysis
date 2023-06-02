"""Helper functions to use with stats_overall.py.

"""
import os

import numpy as np
from tqdm.auto import tqdm

from ad_sensitivity_analysis.data_handler.loader import load_dataset_part


def get_phase_flow_combined_counts(
    file_path,
    inoutflow_time=240,
):
    """

    Parameters
    ----------
    file_path
    inoutflow_time

    Returns
    -------

    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    phase_flows = {
        "warm phase inflow": 0,
        "warm phase ascent": 0,
        "warm phase outflow": 0,
        "mixed phase inflow": 0,
        "mixed phase ascent": 0,
        "mixed phase outflow": 0,
        "ice phase inflow": 0,
        "ice phase ascent": 0,
        "ice phase outflow": 0,
        "neutral phase inflow": 0,
        "neutral phase ascent": 0,
        "neutral phase outflow": 0,
    }

    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            inoutflow_time=inoutflow_time,
            load_params=[["phase", "time_after_ascent", "asc600", "T"]],
            lat_bug=True,
        )
        in_times = np.nansum(
            ~np.isnan(
                ds.where((ds["phase"] == 0) & (ds["time_after_ascent"] < 0))[
                    "phase"
                ].values
            )
        )
        phase_flows["warm phase inflow"] += in_times
        asc_times = np.nansum(
            ~np.isnan(
                ds.where((ds["phase"] == 0) & (ds["asc600"] == 1))["phase"].values
            )
        )
        phase_flows["warm phase ascent"] += asc_times
        all_times = np.nansum(
            ~np.isnan(ds.where(ds["phase"] == 0)["time_after_ascent"].values)
        )
        phase_flows["warm phase outflow"] += all_times - asc_times - in_times

        in_times = np.nansum(
            ~np.isnan(
                ds.where((ds["phase"] == 1) & (ds["time_after_ascent"] < 0))[
                    "phase"
                ].values
            )
        )
        phase_flows["mixed phase inflow"] += in_times
        asc_times = np.nansum(
            ~np.isnan(
                ds.where((ds["phase"] == 1) & (ds["asc600"] == 1))["phase"].values
            )
        )
        phase_flows["mixed phase ascent"] += asc_times
        all_times = np.nansum(
            ~np.isnan(ds.where(ds["phase"] == 1)["time_after_ascent"].values)
        )
        phase_flows["mixed phase outflow"] += all_times - asc_times - in_times

        in_times = np.nansum(
            ~np.isnan(
                ds.where((ds["phase"] == 2) & (ds["time_after_ascent"] < 0))[
                    "phase"
                ].values
            )
        )
        if in_times > 0:
            print(f)
        phase_flows["ice phase inflow"] += in_times
        asc_times = np.nansum(
            ~np.isnan(
                ds.where((ds["phase"] == 2) & (ds["asc600"] == 1))["phase"].values
            )
        )
        phase_flows["ice phase ascent"] += asc_times
        all_times = np.nansum(
            ~np.isnan(ds.where(ds["phase"] == 2)["time_after_ascent"].values)
        )
        phase_flows["ice phase outflow"] += all_times - asc_times - in_times

        in_times = np.nansum(
            ~np.isnan(
                ds.where((ds["phase"] == 3) & (ds["time_after_ascent"] < 0))[
                    "phase"
                ].values
            )
        )
        phase_flows["neutral phase inflow"] += in_times
        asc_times = np.nansum(
            ~np.isnan(
                ds.where((ds["phase"] == 3) & (ds["asc600"] == 1))["phase"].values
            )
        )
        phase_flows["neutral phase ascent"] += asc_times
        all_times = np.nansum(
            ~np.isnan(ds.where(ds["phase"] == 3)["time_after_ascent"].values)
        )
        phase_flows["neutral phase outflow"] += all_times - asc_times - in_times
    return phase_flows


def get_phase_flow_counts(
    file_path,
    inoutflow_time=240,
):
    """

    Parameters
    ----------
    file_path
    inoutflow_time

    Returns
    -------

    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    n_phases = {
        "warm phase": 0,
        "mixed phase": 0,
        "ice phase": 0,
        "neutral phase": 0,
    }
    n_flows = {
        "inflow": 0,
        "ascent": 0,
        "outflow": 0,
    }
    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            inoutflow_time=inoutflow_time,
            load_params=[["phase", "time_after_ascent", "asc600", "lat", "T"]],
            lat_bug=True,
        )
        n_phases["warm phase"] += np.nansum(
            ~np.isnan(ds.where(ds["phase"] == 0)["phase"].values)
        )
        n_phases["mixed phase"] += np.nansum(
            ~np.isnan(ds.where(ds["phase"] == 1)["phase"].values)
        )
        n_phases["ice phase"] += np.nansum(
            ~np.isnan(ds.where(ds["phase"] == 2)["phase"].values)
        )
        n_phases["neutral phase"] += np.nansum(
            ~np.isnan(ds.where(ds["phase"] == 3)["phase"].values)
        )

        in_times = np.nansum(
            ~np.isnan(ds.where(ds["time_after_ascent"] < 0)["time_after_ascent"].values)
        )
        n_flows["inflow"] += in_times
        asc_times = np.nansum(~np.isnan(ds.where(ds["asc600"] == 1)["asc600"]))
        n_flows["ascent"] += asc_times
        all_times = np.nansum(~np.isnan(ds["time_after_ascent"]))
        n_flows["outflow"] += all_times - asc_times - in_times
    return n_phases, n_flows


# pylint: disable=too-many-locals
def get_phase_flow_combined_stats(
    file_path,
    n_phase_flows,
    inoutflow_time=240,
):
    """

    Parameters
    ----------
    file_path
    n_phase_flows
    inoutflow_time

    Returns
    -------

    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    # min, max, mean, std
    temp_dict = {
        "warm phase inflow": [np.Inf, np.NINF, 0, 0],
        "mixed phase inflow": [np.Inf, np.NINF, 0, 0],
        "ice phase inflow": [np.Inf, np.NINF, 0, 0],
        "neutral phase inflow": [np.Inf, np.NINF, 0, 0],
        "warm phase ascent": [np.Inf, np.NINF, 0, 0],
        "mixed phase ascent": [np.Inf, np.NINF, 0, 0],
        "ice phase ascent": [np.Inf, np.NINF, 0, 0],
        "neutral phase ascent": [np.Inf, np.NINF, 0, 0],
        "warm phase outflow": [np.Inf, np.NINF, 0, 0],
        "mixed phase outflow": [np.Inf, np.NINF, 0, 0],
        "ice phase outflow": [np.Inf, np.NINF, 0, 0],
        "neutral phase outflow": [np.Inf, np.NINF, 0, 0],
    }
    phases = [
        "warm phase",
        "mixed phase",
        "ice phase",
        "neutral phase",
    ]
    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            inoutflow_time=inoutflow_time,
            load_params=[["phase", "time_after_ascent", "asc600", "lat", "T"]],
            lat_bug=True,
        )

        def _get_min(ds_tmp, phase, key_idx, time_cond, flow_key):
            t_new = (
                ds_tmp.where((ds_tmp["phase"] == key_idx) & (time_cond))["T"]
                .min()
                .item()
            )
            if temp_dict[f"{phase} {flow_key}"][0] > t_new:
                temp_dict[f"{phase} {flow_key}"][0] = t_new

        def _get_max(ds_tmp, phase, key_idx, time_cond, flow_key):
            t_new = (
                ds_tmp.where((ds_tmp["phase"] == key_idx) & (time_cond))["T"]
                .max()
                .item()
            )
            if temp_dict[f"{phase} {flow_key}"][1] < t_new:
                temp_dict[f"{phase} {flow_key}"][1] = t_new

        def _get_mean(ds_tmp, phase, key_idx, n, time_cond, flow_key):
            if n != 0:
                temp_dict[f"{phase} {flow_key}"][2] += (
                    n
                    * ds_tmp.where((ds_tmp["phase"] == key_idx) & (time_cond))["T"]
                    .mean(skipna=True)
                    .item()
                    / n_phase_flows[f"{phase} {flow_key}"]
                )

        def _get_var(ds_tmp, phase, key_idx, n, time_cond, flow_key):
            if n != 0:
                temp_dict[f"{phase} {flow_key}"][3] += (
                    n
                    * (
                        ds_tmp.where((ds_tmp["phase"] == key_idx) & (time_cond))["T"]
                        ** 2
                    )
                    .mean(skipna=True)
                    .item()
                    / n_phase_flows[f"{phase} {flow_key}"]
                )

        inflow_cond = ds["time_after_ascent"] < 0
        ascent_cond = ds["asc600"] == 1
        outflow_cond = ~(inflow_cond | ascent_cond)

        for flow_key, cond in zip(
            ["inflow", "ascent", "outflow"], [inflow_cond, ascent_cond, outflow_cond]
        ):
            for phase_key, p_i in zip(phases, [0, 1, 2, 3]):
                _get_min(ds, phase_key, p_i, cond, flow_key)
                _get_max(ds, phase_key, p_i, cond, flow_key)
                n = np.nansum(
                    ~np.isnan(ds.where((ds["phase"] == p_i) & cond)["phase"].values)
                )
                if n != 0:
                    _get_mean(ds, phase_key, p_i, n, cond, flow_key)
                    _get_var(ds, phase_key, p_i, n, cond, flow_key)

    for key, stats in temp_dict.items():
        stats[3] = stats[3] - stats[2] ** 2
    return temp_dict


def _get_flow_stats(
    ds,
    temp_flows,
    n_flows,
):
    """

    Parameters
    ----------
    ds
    temp_flows
    n_flows

    Returns
    -------

    """
    t_new = ds.where(ds["time_after_ascent"] < 0)["T"].min().item()
    if temp_flows["inflow"][0] > t_new:
        temp_flows["inflow"][0] = t_new
    t_new = ds.where(ds["asc600"] == 1)["T"].min().item()
    if temp_flows["ascent"][0] > t_new:
        temp_flows["ascent"][0] = t_new
    time_cond = ds["time_after_ascent"] < 0
    asc_cond = ds["asc600"] == 1
    cond = ~(time_cond | asc_cond)
    t_new = ds.where(cond)["T"].min().item()
    if temp_flows["outflow"][0] > t_new:
        temp_flows["outflow"][0] = t_new

    t_new = ds.where(ds["time_after_ascent"] < 0)["T"].max().item()
    if temp_flows["inflow"][1] < t_new:
        temp_flows["inflow"][1] = t_new
    t_new = ds.where(ds["asc600"] == 1)["T"].max().item()
    if temp_flows["ascent"][1] < t_new:
        temp_flows["ascent"][1] = t_new
    t_new = ds.where(cond)["T"].max().item()
    if temp_flows["outflow"][1] < t_new:
        temp_flows["outflow"][1] = t_new

    in_times = np.nansum(
        ~np.isnan(ds.where(ds["time_after_ascent"] < 0)["time_after_ascent"].values)
    )
    asc_times = np.nansum(~np.isnan(ds.where(ds["asc600"] == 1)["asc600"]))
    out_times = np.nansum(~np.isnan(ds["time_after_ascent"])) - asc_times - in_times
    if in_times != 0:
        t_new = ds.where(ds["time_after_ascent"] < 0)["T"].mean(skipna=True).item()
        temp_flows["inflow"][2] += t_new * in_times / n_flows["inflow"]

        t_new = (
            (ds.where(ds["time_after_ascent"] < 0)["T"] ** 2).mean(skipna=True).item()
        )
        temp_flows["inflow"][3] += t_new * in_times / n_flows["inflow"]

    if asc_times != 0:
        t_new = ds.where(ds["asc600"] == 1)["T"].mean(skipna=True).item()
        temp_flows["ascent"][2] += t_new * asc_times / n_flows["ascent"]

        t_new = (ds.where(ds["asc600"] == 1)["T"] ** 2).mean(skipna=True).item()
        temp_flows["ascent"][3] += t_new * asc_times / n_flows["ascent"]

    if out_times != 0:
        t_new = ds.where(cond)["T"].mean(skipna=True).item()
        temp_flows["outflow"][2] += t_new * out_times / n_flows["outflow"]

        t_new = (ds.where(cond)["T"] ** 2).mean(skipna=True).item()
        temp_flows["outflow"][3] += t_new * out_times / n_flows["outflow"]
    return temp_flows


def _get_phase_stats(
    ds,
    temp_phases,
    n_phases,
):
    """

    Parameters
    ----------
    ds
    temp_phases
    n_phases

    Returns
    -------

    """

    def _get_min_phase(phase, key_idx):
        t_tmp = ds.where(ds["phase"] == key_idx)["T"].min().item()
        if temp_phases[phase][0] > t_tmp:
            temp_phases[phase][0] = t_tmp

    def _get_max_phase(phase, key_idx):
        t_tmp = ds.where(ds["phase"] == key_idx)["T"].max().item()
        if temp_phases[phase][1] < t_tmp:
            temp_phases[phase][1] = t_tmp

    def _get_mean_phase(phase, key_idx, n):
        if n != 0:
            temp_phases[phase][2] += (
                n
                * ds.where(ds["phase"] == key_idx)["T"].mean(skipna=True).item()
                / n_phases[phase]
            )

    def _get_var_phase(phase, key_idx, n):
        if n != 0:
            temp_phases[phase][3] += (
                n
                * (ds.where(ds["phase"] == key_idx)["T"] ** 2).mean(skipna=True).item()
                / n_phases[phase]
            )

    _get_min_phase("warm phase", 0)
    _get_min_phase("mixed phase", 1)
    _get_min_phase("ice phase", 2)
    _get_min_phase("neutral phase", 3)

    _get_max_phase("warm phase", 0)
    _get_max_phase("mixed phase", 1)
    _get_max_phase("ice phase", 2)
    _get_max_phase("neutral phase", 3)

    n_warm = np.nansum(~np.isnan(ds.where(ds["phase"] == 0)["phase"].values))
    n_mixed = np.nansum(~np.isnan(ds.where(ds["phase"] == 1)["phase"].values))
    n_ice = np.nansum(~np.isnan(ds.where(ds["phase"] == 2)["phase"].values))
    n_neut = np.nansum(~np.isnan(ds.where(ds["phase"] == 3)["phase"].values))

    _get_mean_phase("warm phase", 0, n_warm)
    _get_mean_phase("mixed phase", 1, n_mixed)
    _get_mean_phase("ice phase", 2, n_ice)
    _get_mean_phase("neutral phase", 3, n_neut)

    _get_var_phase("warm phase", 0, n_warm)
    _get_var_phase("mixed phase", 1, n_mixed)
    _get_var_phase("ice phase", 2, n_ice)
    _get_var_phase("neutral phase", 3, n_neut)
    return temp_phases


def get_phase_flow_stats(
    file_path,
    n_phases,
    n_flows,
    inoutflow_time=240,
):
    """

    Parameters
    ----------
    file_path
    n_phases
    n_flows
    inoutflow_time

    Returns
    -------

    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    files = np.sort(files)
    temp_phases = {
        "warm phase": [np.Inf, np.NINF, 0, 0],
        "mixed phase": [np.Inf, np.NINF, 0, 0],
        "ice phase": [np.Inf, np.NINF, 0, 0],
        "neutral phase": [np.Inf, np.NINF, 0, 0],
    }
    temp_flows = {
        "inflow": [np.Inf, np.NINF, 0, 0],
        "ascent": [np.Inf, np.NINF, 0, 0],
        "outflow": [np.Inf, np.NINF, 0, 0],
    }
    for f in tqdm(files):
        ds = load_dataset_part(
            f=file_path + f,
            inoutflow_time=inoutflow_time,
            load_params=[["phase", "time_after_ascent", "asc600", "lat", "T"]],
            lat_bug=True,
        )

        temp_phases = _get_phase_stats(
            ds=ds,
            temp_phases=temp_phases,
            n_phases=n_phases,
        )
        temp_flows = _get_flow_stats(
            ds=ds,
            temp_flows=temp_flows,
            n_flows=n_flows,
        )

    for key, stats in temp_phases.items():
        stats[3] = stats[3] - stats[2] ** 2
    for key, stats in temp_flows.items():
        stats[3] = stats[3] - stats[2] ** 2
    return temp_phases, temp_flows
