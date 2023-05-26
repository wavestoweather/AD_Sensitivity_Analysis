"""Test if phases are properly set

"""
import os

import numpy as np
from tqdm.auto import tqdm

from ad_sensitivity_analysis.tests.color_scheme import COLOR_RESET, SUCCESS, ERROR


def calculate_phase(ds, thresholds):
    """
    Calculate the phase based on given thresholds for each hydrometeor type.

    Parameters
    ----------
    ds
    thresholds

    Returns
    -------

    """
    phase_col = np.full(
        (len(ds["ensemble"]), len(ds["trajectory"]), len(ds["time"])),
        "             ",
    )

    def warm(ds):
        return (
            (ds["QC"] > thresholds["warm_q_phase_threshold"])
            | (ds["QR"] > thresholds["warm_q_phase_threshold"])
            | (ds["NCCLOUD"] > thresholds["warm_n_phase_threshold"])
            | (ds["NCRAIN"] > thresholds["warm_n_phase_threshold"])
        )

    def cold(ds):
        return (
            (ds["QG"] > thresholds["ice_q_phase_threshold"])
            | (ds["QH"] > thresholds["ice_q_phase_threshold"])
            | (ds["NCGRAUPEL"] > thresholds["ice_n_phase_threshold"])
            | (ds["NCHAIL"] > thresholds["ice_n_phase_threshold"])
            | (ds["QI"] > thresholds["ice_q_phase_threshold"])
            | (ds["QS"] > thresholds["ice_q_phase_threshold"])
            | (ds["NCICE"] > thresholds["ice_n_phase_threshold"])
            | (ds["NCSNOW"] > thresholds["ice_n_phase_threshold"])
        )

    def warm_phase(ds):
        return (warm(ds)) & (~cold(ds))

    def cold_phase(ds):
        return (~warm(ds)) & (cold(ds))

    def mixed_phase(ds):
        return (warm(ds)) & (cold(ds))

    def neutral_phase(ds):
        return (~warm(ds)) & (~cold(ds))

    def nan_phase(ds):
        return np.isnan(ds["QV"])

    phase_col[np.where(warm_phase(ds))] = "warm phase   "
    phase_col[np.where(cold_phase(ds))] = "ice phase    "
    phase_col[np.where(mixed_phase(ds))] = "mixed phase  "
    phase_col[np.where(neutral_phase(ds))] = "neutral phase"
    phase_col[np.where(nan_phase(ds))] = "nan"
    return phase_col


def add_phase(ds, col_name="phase"):
    """

    Parameters
    ----------
    ds
    col_name

    Returns
    -------

    """

    thresholds = {
        "ice_q_phase_threshold": 0.0,
        "ice_n_phase_threshold": 0.0,
        "warm_q_phase_threshold": 0.0,
        "warm_n_phase_threshold": 0.0,
    }
    # Extract used thresholds from C++-code
    current_dir = os.getcwd()
    if "/end_to_end" == current_dir[-11::]:
        constants_path = current_dir + "/../../../include/microphysics/constants.h"
    else:
        constants_path = current_dir + "/include/microphysics/constants.h"
    with open(constants_path, "r", encoding="utf-8") as constants_f:
        for line in constants_f:
            for thresh in thresholds:
                if thresh in line:
                    thresholds[thresh] = float(line.split(" = ")[-1].replace(";", ""))

    ds[col_name] = (("ensemble", "trajectory", "time"), calculate_phase(ds, thresholds))
    return ds


def rename_phase(ds):
    """

    Parameters
    ----------
    ds

    Returns
    -------

    """
    phase_col = np.full(
        (len(ds["ensemble"]), len(ds["trajectory"]), len(ds["time"])),
        "             ",
    )
    phase_col[np.where(ds["phase"] == 0)] = "warm phase   "
    phase_col[np.where(ds["phase"] == 2)] = "ice phase    "
    phase_col[np.where(ds["phase"] == 1)] = "mixed phase  "
    phase_col[np.where(ds["phase"] == 3)] = "neutral phase"
    phase_col[
        np.where(
            (ds["phase"] != 3)
            & (ds["phase"] != 2)
            & (ds["phase"] != 1)
            & (ds["phase"] != 0)
            & (ds["phase"] != 3)
        )
    ] = "nan"
    ds["phase"] = (("ensemble", "trajectory", "time"), phase_col)
    return ds


def test_phases(ds, recalc):
    """
    For each phase, print the amount of trajectories that have at least one time step with this phase.
    ALso prints the percentage and number of time steps for each phase over all trajectories.
    Calculates the phases if necessary.
    If recalc is false, tests if the phases are correct.
    This test is only correct if the simulation did not use sedimentation from above.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    recalc : bool
        Recalculate the phases.

    Returns
    -------
    Number of trajectories without any warm phase.
    """
    print("~*~*~*~Testing amount of trajectories with different phases~*~*~*~")
    err = 0
    n_total_timesteps = len(ds["trajectory"]) * len(ds["time"])
    n_trajectories = len(ds["trajectory"])

    if "phase" not in ds or recalc:
        ds = add_phase(ds)
    else:
        ds = rename_phase(ds)

    if not recalc:
        ds = add_phase(ds, "phase_reference")
        # pylint: disable=no-member
        err += np.sum(ds["phase"] != ds["phase_reference"]).values.item()
        if err == 0:
            print(
                f"{SUCCESS}Phases are correct for all {n_total_timesteps} time steps{COLOR_RESET}\n"
            )
        else:
            perc = err / n_total_timesteps
            print(
                f"{ERROR}Phases are not correct in {err} / {n_total_timesteps} time steps "
                f"({perc*100:2.2f}%){COLOR_RESET}\n"
            )

    n_data = {phase.item(): 0 for phase in np.unique(ds["phase"])}
    n_trajs = {phase.item(): 0 for phase in np.unique(ds["phase"])}

    for phase in tqdm(np.unique(ds["phase"])):
        # pylint: disable=no-member
        n = np.sum(ds["phase"] == phase)
        n_data[phase] += n.values.item()

        n_per_traj = (ds["phase"] == phase).sum(axis=2)
        zero_times = (n_per_traj == 0).sum()
        n_trajs[phase] += (n_trajectories - zero_times).values.item()

    for phase in tqdm(n_data.keys()):
        if phase in ("warm phase", "warm phase   ", 0):
            perc = n_trajs[phase] / n_trajectories
            if n_trajs[phase] < n_trajectories:
                print(
                    f"{ERROR}Failed: {phase} occurs only in {n_trajs[phase]} / {n_trajectories} trajectories "
                    f"({perc*100:2.2f}%){COLOR_RESET}\n"
                )
                err += 1
            else:
                print(
                    f"{SUCCESS}{phase} occurs in {n_trajs[phase]} / {n_trajectories} trajectories "
                    f"({perc*100:2.2f}%){COLOR_RESET}\n"
                )
        else:
            perc = n_trajs[phase] / n_trajectories
            print(
                f"{phase} occurs in {n_trajs[phase]} / {n_trajectories} trajectories ({perc*100:2.2f}%)\n"
            )
        perc = n_data[phase] / n_total_timesteps
        print(
            f"{phase} occurs in {n_data[phase]} / {n_total_timesteps} time steps ({perc*100:2.2f}%)\n"
        )
    return err
