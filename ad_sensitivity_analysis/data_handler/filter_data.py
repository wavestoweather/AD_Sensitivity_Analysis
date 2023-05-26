"""Functions to filter datasets.

"""
import numpy as np


def filter_by_flow(ds, flow, inoutflow_time=-1):
    """

    Parameters
    ----------
    ds
    flow
    inoutflow_time

    Returns
    -------

    """
    if flow == "ascent":
        return ds.where(ds["asc600"] == 1)
    if flow == "any":
        return ds
    if flow == "inflow":
        ds_flow = ds.where(ds["asc600"] == 1)["asc600"]
        ds_flow = (
            ds_flow.rolling(
                time=inoutflow_time * 2,  # once for inflow, another for outflow
                min_periods=1,
                center=True,
            )
            .reduce(lambda x, axis: np.nansum(x, axis=2))
            .diff(dim="time", label="lower")
        )
        return ds.where(ds_flow)
    # outflow
    ds_flow = ds.where(ds["asc600"] == 1)["asc600"]
    ds_flow = (
        ds_flow.rolling(
            time=inoutflow_time * 2,  # once for inflow, another for outflow
            min_periods=1,
            center=True,
        )
        .reduce(lambda x, axis: np.nansum(x, axis=2))
        .diff(dim="time", label="upper")
    )
    return ds.where(ds_flow == -1)


def filter_trajectories(
    ds,
    only_asc600=False,
    sens_model_state_ids=None,
    inoutflow_time=-1,
    min_pressure=None,
    max_pressure=None,
    only_phase=None,
    lat_bug=False,
):
    """
    Filter trajectories based on pressure, inflow- and outflow-time, or ascent.

    Parameters
    ----------
    ds
    only_asc600
    sens_model_state_ids : list of
    inoutflow_time
    min_pressure
    max_pressure

    Returns
    -------

    """
    if inoutflow_time > 0 and "asc600" in ds:
        ds_flow = ds.where(ds["asc600"] == 1)["asc600"]
        ds_flow = ds_flow.rolling(
            time=inoutflow_time * 2,  # once for inflow, another for outflow
            min_periods=1,
            center=True,
        ).reduce(np.nanmax)
        ds = ds.where(ds_flow == 1)
    elif only_asc600 and "asc600" in ds:
        ds = ds.where(ds["asc600"] == 1)
    if min_pressure is not None:
        ds = ds.where(ds["pressure"] >= min_pressure)
    if max_pressure is not None:
        ds = ds.where(ds["pressure"] <= max_pressure)
    if sens_model_state_ids is not None and "Output_Parameter_ID" in ds:
        ds = ds.sel({"Output_Parameter_ID": sens_model_state_ids})
    if only_phase is not None:
        if ds["phase"].dtype not in (str, np.uint64):
            ds["phase"] = ds["phase"].astype(np.uint64)
            phase_idx = np.argwhere(
                np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])
                == only_phase
            )[0].item()
            ds = ds.where(ds["phase"] == phase_idx)
        elif ds["phase"].dtype == str:
            ds = ds.where(ds["phase"] == only_phase)
        else:
            phases = np.asarray(
                ["warm phase", "mixed phase", "ice phase", "neutral phase"]
            )
            phase_idx = np.argwhere(phases == only_phase)[0].item()
            ds = ds.where(ds["phase"] == phase_idx)
    if lat_bug:
        # There seems to be a bug with some output data where the longitude
        # and latitude happen to be zero when the trajectory is already finished.
        # Since this is far away from our domain, we can savely delete that.
        # It is save to assume that only latitude with zero values are invalid
        # ds["lon"] = ds["lon"].where(ds["lon"] != 0)
        ds = ds.where(ds["lat"] != 0)
    return ds


def filter_rank_data(ds, out_param, in_params, ignore_zero_gradients, phase, flow):
    """

    Parameters
    ----------
    ds
    out_param
    in_params
    ignore_zero_gradients
    phase
    flow

    Returns
    -------

    """
    ds_tmp = ds.sel({"Output Parameter": out_param}).drop("Output Parameter")
    df = ds_tmp[in_params].to_dataframe().stack().reset_index()
    df = df.rename(columns={"level_4": "Parameter", 0: "Rank"})
    df = df.loc[df["Rank"] > 0]
    worst_rank = ds[in_params[0]].attrs["rank for zero gradients"]
    if ignore_zero_gradients:
        df = df.loc[df["Rank"] < worst_rank]

    if phase and flow:
        df = df.loc[(df["phase"] != "any") & (df["phase"] != "neutral phase")]
    elif phase:
        df = df.loc[
            (df["phase"] != "any")
            & (df["phase"] != "neutral phase")
            & (df["flow"] == "any")
        ]
    elif flow:
        df = df.loc[df["phase"] == "any"]
    else:
        df = df.loc[df["phase"] == "any"]
        df = df.loc[df["flow"] == "any"]
    return df, worst_rank


def d_unnamed(df):
    """
    Remove unnamed column from dataframe.

    Parameters
    ----------
    df : pandas.Dataframe
        A dataframe with one or more unnamed columns

    Returns
    -------
    Dataframe without unnamed columns.
    """
    return df.loc[:, ~df.columns.str.contains("^Unnamed")]
