import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np
import os
import panel as pn
import pickle
import seaborn as sns
import sys
from tqdm.auto import tqdm
import xarray as xr

try:
    from get_stats import load_ds
    import latexify
except:
    from scripts.get_stats import load_ds
    import scripts.latexify as latexify

pn.extension()


def create_rank_traj_dataset(file_path, inoutflow_time=240, model_params=None):
    """
    Ranked index of parameters for each trajectory.
    Count model parameter / process index occurrence over all trajectories.
    Additional information that might be useful:
    location of ascent. Average ascent. Number of time steps with ascent.

    Parameters
    ----------
    file_path : string
        Path to a set of NetCDF-files from a sensitivity analysis.
    inoutflow_time : int
        Number of time steps before and after the ascent that shall be used additionally.
    model_params : List of strings or None
        List of model parameters to create a rank for. If None is given, all model parameters are used.

    Returns
    -------
    xarray.Dataset:
    dims:   Output Parameter (QV, latent_heat, latent_cool)
            trajectory       (0, 1, 2, ...)
            file             (traj20161003.nc, ...)
            # index            (0, 1, 2, ...)
            phase            (warm phase, mixed phase, ice phase, neutral phase, any)
            flow             (inflow, ascent, outflow, any)

    vars:   Model Parameter Rank (Output Parameter, file, trajectory, phase, flow)
            Model Parameter Avg (Output Parameter, file, trajectory, phase, flow)
            lon             (file, trajectory)
            lat             (file, trajectory)
            avg ascent      (file, trajectory, phase)
            asc600 steps    (file, trajectory, phase)
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    lon_lat_dims = ["file", "trajectory"]
    addit_dims = ["file", "trajectory", "phase"]
    param_dims = ["Output Parameter", "file", "trajectory", "phase", "flow"]
    ds = load_ds(file_path + files[0])
    phase_type = ds["phase"].dtype
    if model_params is None:
        model_params = []
        for col in ds:
            if col[0] == "d" and col != "deposition":
                model_params.append(col)
    if "Output_Parameter_ID" in ds:
        out_coord = "Output_Parameter_ID"
        out_params = []
        for out_p in ds["Output_Parameter_ID"]:
            out_params.append(latexify.param_id_map[out_p.item()])
        out_params_coord_id = ds["Output_Parameter_ID"].values
    else:
        out_coord = "Output Parameter"
        out_params = ds["Output Parameter"].values
        out_params_coord_id = ds["Output Parameter"].values
    max_traj = np.max(ds["trajectory"]).item()
    for f in files:
        ds = load_ds(file_path + f, load_params="trajectory")
        if np.max(ds["trajectory"]).item() > max_traj:
            max_traj = np.max(ds["trajectory"]).item()
    max_traj += 1  # Trajectory enumeration starts with zero...
    coords = {
        "Output Parameter": out_params,
        "file": files,
        "trajectory": np.arange(max_traj),
        "phase": ["warm phase", "mixed phase", "ice phase", "neutral phase", "any"],
        "flow": ["inflow", "ascent", "outflow", "any"],
    }
    model_parameter_index = {
        param: np.zeros(
            (
                len(out_params),
                len(files),
                max_traj,
                len(coords["phase"]),
                len(coords["flow"]),
            )
        )
        for param in model_params
    }
    model_parameter_avg = {
        param: np.zeros(
            (
                len(out_params),
                len(files),
                max_traj,
                len(coords["phase"]),
                len(coords["flow"]),
            )
        )
        for param in model_params
    }
    worst_rank = len(model_params) + 1
    lon = np.zeros((len(files), max_traj))
    lat = np.zeros((len(files), max_traj))
    avg_ascent = np.zeros((len(files), max_traj, len(coords["phase"])))
    asc600_steps = np.zeros((len(files), max_traj, len(coords["phase"])))
    phases = np.asarray(["warm phase", "mixed phase", "ice phase", "neutral phase"])

    def get_phase(ds_tmp, phase):
        if phase == "any":
            return ds_tmp
        if phase_type == str:
            phase_idx = phase
        else:
            phase_idx = np.argwhere(phases == phase)[0].item()
        return ds_tmp.where(ds_tmp["phase"] == phase_idx)

    def get_flow(ds_tmp, flow):
        if flow == "ascent":
            return ds_tmp.where(ds_tmp["asc600"] == 1)
        elif flow == "any":
            return ds_tmp
        if flow == "inflow":
            ds_flow = ds_tmp.where(ds_tmp["asc600"] == 1)["asc600"]
            ds_flow = (
                ds_flow.rolling(
                    time=inoutflow_time * 2,  # once for inflow, another for outflow
                    min_periods=1,
                    center=True,
                )
                .reduce(lambda x, axis: np.nansum(x, axis=2))
                .diff(dim="time", label="lower")
            )
            return ds_tmp.where(ds_flow)
        # outflow
        ds_flow = ds_tmp.where(ds_tmp["asc600"] == 1)["asc600"]
        ds_flow = (
            ds_flow.rolling(
                time=inoutflow_time * 2,  # once for inflow, another for outflow
                min_periods=1,
                center=True,
            )
            .reduce(lambda x, axis: np.nansum(x, axis=2))
            .diff(dim="time", label="upper")
        )
        return ds_tmp.where(ds_flow == -1)

    def get_ranking(ds_tmp):
        ds_avg = ds_tmp[model_params].mean(dim="time", skipna=True)
        param_avg = [[] for _ in ds_avg["trajectory"].values]
        for traj_idx in ds_avg["trajectory"].values:
            ds_avg2 = ds_avg.sel({"trajectory": traj_idx})
            for param in ds_avg:
                param_avg[traj_idx].append((param, ds_avg2[param].values.item()))
        return param_avg

    for f_i, f in enumerate(tqdm(files)):
        ds = load_ds(
            file_path + f,
            inoutflow_time=inoutflow_time,
            load_params=model_params + ["lon", "lat", "w", "asc600", "phase"],
        )
        ds = ds.isel({"ensemble": 0})
        ds[model_params] = np.abs(ds[model_params])
        ds_lons_lats = ds[["lon", "lat", "asc600"]].where(ds["asc600"])[["lon", "lat"]]
        lons = ds_lons_lats["lon"]
        lats = ds_lons_lats["lat"]
        for traj_idx in range(len(lons.values)):
            lon[f_i, traj_idx] = lons[traj_idx][~np.isnan(lons[traj_idx])][0]
            lat[f_i, traj_idx] = lats[traj_idx][~np.isnan(lats[traj_idx])][0]

        for phase_i, phase in enumerate(coords["phase"]):
            ds_phase = get_phase(ds, phase)
            avg_ascent_tmp = ds_phase["asc600"].count(dim="time").values
            asc600_steps_tmp = ds_phase["w"].mean(dim="time")
            for traj_idx in range(len(avg_ascent_tmp)):
                avg_ascent[f_i, traj_idx, phase_i] = avg_ascent_tmp[traj_idx]
                asc600_steps[f_i, traj_idx, phase_i] = asc600_steps_tmp[traj_idx]
            for flow_i, flow in enumerate(coords["flow"]):
                ds_flow = get_flow(ds_phase, flow)
                for out_p_i, out_p in enumerate(out_params_coord_id):
                    ds_out = ds_flow.sel({out_coord: out_p})
                    param_values = get_ranking(ds_out)
                    for traj_idx in range(len(param_values)):
                        param_values[traj_idx].sort(key=lambda x: x[1])
                        param_values[traj_idx] = param_values[traj_idx][::-1]
                        last_val = 0
                        rank_offset = 1  # zero is the "invalid number" in our case
                        for rank, param_pair in enumerate(param_values[traj_idx]):
                            if rank > 0:
                                if last_val == param_pair[1]:
                                    rank_offset -= 1
                            else:
                                last_val = param_pair[1]
                            if param_pair[1] == 0:
                                model_parameter_index[param_pair[0]][
                                    out_p_i, f_i, traj_idx, phase_i, flow_i
                                ] = worst_rank
                            else:
                                model_parameter_index[param_pair[0]][
                                    out_p_i, f_i, traj_idx, phase_i, flow_i
                                ] = (rank + rank_offset)
                            model_parameter_avg[param_pair[0]][
                                out_p_i, f_i, traj_idx, phase_i, flow_i
                            ] = param_pair[1]

    data_vars = {
        f"{param} rank": (param_dims, model_parameter_index[param])
        for param in model_parameter_index
    }
    for param in model_parameter_avg:
        data_vars[f"{param} avg"] = (param_dims, model_parameter_avg[param])
    data_vars["lon"] = (lon_lat_dims, lon)
    data_vars["lat"] = (lon_lat_dims, lat)
    data_vars["avg ascent"] = (addit_dims, avg_ascent)
    data_vars["asc600 step count"] = (addit_dims, asc600_steps)
    ds = xr.Dataset(data_vars=data_vars, coords=coords)
    for param in model_parameter_index:
        ds[f"{param} rank"].attrs = {
            "no data": 0,
            "rank for zero gradients": worst_rank,
        }
    return ds


def create_rank_latex_table(ds, phase, flow):
    """
    Create multiple latex tables with median rank, median absolute deviation, mean rank, and standard deviation
    for a given phase and flow.

    Parameters
    ----------
    ds : xarray.Dataset
        Ranked index of each trajectory created by create_rank_traj_dataset().
    phase : string
        Name of the phase to use. Options are those in ds["phase"], which typically amounts to
        "warm phase", "mixed phase", "ice phase", "any".
    flow : string
        Name of the 'flow' to use. Options are those in ds["flow"], which typically amounts to
        "inflow", "ascent", "outflow", "any".

    Returns
    -------
    Dictionary of raw strings (latex tables). The keys are the model state variables.
    """
    table_text = r"""
    \begin{table}[ht]
        \centering
        \begin{tabular}{l|l|c|c|c|c}
            \textbf{Model State Variable} & \textbf{{Parameter}} & \textbf{Median Rank} & \textbf{MAD} & \textbf{Mean Rank} & \textbf{STD} \\ \hline 
    """
    ds2 = ds.sel({"phase": phase, "flow": flow})
    for o_p in ds["Output Parameter"]:
        ds_op = ds2.sel({"Output Parameter": o_p})
        table_text += f"        {latexify.parse_word(o_p.item()):<25}& "
        p_i = 0
        for param in ds:
            if not "rank" in param:
                continue

            ds_op2 = ds_op[param]
            ds_op2 = ds_op2.where(ds_op2 > 0)
            med = ds_op2.median().item()
            if med > 10 and ds_op2.mean().item() > 10:
                continue
            if np.isnan(med):
                continue
            mad = ds_op2 - med
            mad = mad.median().item()
            if p_i > 0:
                table_text += "                                 & "
            p_i += 1
            word = latexify.parse_word(param[:-5]).replace(r"\partial ", "")
            table_text += f"{word:<40}& "
            table_text += (
                f"{int(med):2d} & {int(mad):2d} & {ds_op2.mean().item():2.2f} & {ds_op2.std().item():2.2f}"
                + r" \\"
                + "\n"
            )

    table_text += r"""
        \end{tabular}
        \caption{"""
    caption = "The median and median absolute deviation of the rank of each model parameter if we calculate the rank over all trajectories. Only those with median or mean lower than 11 are used."
    if phase == "any":
        caption += "Using any phase "
    else:
        caption += f"Using {phase} "
    if flow == "any":
        caption += "and inflow, outflow, and ascent."
    else:
        caption += f"and {flow}."
    caption += f"Using {phase} phase and {flow} flow."
    table_text += f"{caption}"
    table_text += r"""}
        \label{tab:analysis:rank_count}
    \end{table}
    """
    return table_text


def create_rank_latex_table_per_outp(ds, col, only_n=None):
    """
    Create multiple latex tables with the statistic defined in col. The table shows the statistic for any combination
    of flow and phase. Each table is the impact on a different model state variable.

    Parameters
    ----------
    ds : xarray.Dataset
        Ranked index of each trajectory created by create_rank_traj_dataset().
    col : string
        Define which statistic shall be shown. Options are "Median Rank", "Mean Rank", "MAD", "STD". If the
        string is not recognized, it defaults to "STD".
    only_n : int
        Show only model parameters that have a rank lower than 'only_n'.

    Returns
    -------
    Dictionary of raw strings (latex tables). The keys are the model state variables.
    """
    tables = {}
    phases = ["warm phase", "mixed phase", "ice phase"]
    flows = ["inflow", "ascent", "outflow"]
    for o_p in ds["Output Parameter"]:
        ds_op = ds.sel({"Output Parameter": o_p})

        table_text = r"""
\begin{table}[ht]
    \centering
    \begin{tabular}{l|c|c|c}
        \multirow{3}*{\textbf{{Parameter}}}     & \multicolumn{3}{c}{\textbf{"""
        table_text += f"{col}"
        table_text += r"""}} \\ \cline{2-4} 
                                                & inflow & ascent & outflow \\ \cline{2-4} 
                                                & wp, mp, ip & wp, mp, ip & wp, mp, ip \\ \hline
"""
        for param in ds:
            if not "rank" in param:
                continue
            row_vals = "        "
            word = latexify.parse_word(param[:-5]).replace(r"\partial ", "")
            row_vals += f"{word:<40}& "
            got_valid = False
            for flow in flows:
                ds_flow = ds_op.sel({"flow": flow})
                for phase_i, phase in enumerate(phases):
                    ds2 = ds_flow.sel({"phase": phase})
                    ds_op2 = ds2[param]
                    ds_op2 = ds_op2.where(ds_op2 > 0)
                    if col == "Median Rank":
                        val = ds_op2.median().item()
                        if only_n is None:
                            got_valid = True
                        elif only_n >= val:
                            got_valid = True
                    elif col == "Mean Rank":
                        val = ds_op2.mean().item()
                        if only_n is None:
                            got_valid = True
                        elif only_n >= val:
                            got_valid = True
                    elif col == "MAD":
                        got_valid = True
                        med = ds_op2.median().item()
                        diff = ds_op2 - med
                        diff = np.sort(diff.values.flatten())
                        diff = diff[~np.isnan(diff)]
                        if len(diff) == 0:
                            val = np.nan
                        elif len(diff) % 2 == 0:
                            val = (diff[len(diff) // 2] + diff[len(diff) // 2 - 1]) / 2
                        else:
                            val = diff[len(diff) // 2]
                    else:
                        val = ds_op2.std().item()

                    if np.isnan(val):
                        row_vals += " -"
                    else:
                        if col == "Median Rank":
                            row_vals += f"{int(val):02d}"
                        else:
                            row_vals += f"{val:2.2f}"

                    if phase_i < len(phases) - 1:
                        row_vals += ", "
                    elif flow != "outflow":
                        row_vals += " & "
            if got_valid:
                table_text += row_vals + r" \\" + "\n"

        table_text += r"""
    \end{tabular}
    \caption{"""
        caption = "The median and median absolute deviation of the rank of each model parameter if we calculate the rank over all trajectories. "
        caption += "Only those with median or mean lower than 11 are used. Distinction with different phases. "
        caption += "wp = warm phase, mp = mixed phase, ip = ice phase. "
        caption += f"Impact on {latexify.parse_word(o_p.item())} is used here."
        table_text += f"{caption}"
        table_text += r"""}
            \label{tab:analysis:rank_count_variations_"""
        table_text += f"{o_p.item()}"
        table_text += """}
\end{table}
        """
        tables[o_p.item()] = table_text
    return tables


def plot_kde_histogram(
    ds,
    in_params,
    out_param,
    flow,
    phase,
    ignore_zero_gradients,
    linewidth=2,
    bw_adjust=1.0,
    title=None,
    filename=None,
    width=17,
    height=16,
    font_scale=None,
    save=True,
    latex=False,
):
    """
    Plot kde histogram of ranks.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset created via create_rank_traj_dataset().
    in_params : list of strings
        A list of model parameters.
    out_param : string
        Model state variable which the sensitivities are for.
    flow : bool
        Plot different flow phases, i.e., inflow, ascent, and outflow.
    phase : bool
        Plot different phases, i.e., warm phase, mixed phase, and ice phase.
    ignore_zero_gradients : bool
        Do not show the rank of gradients that are always zero.
    linewidth : float
        Line width of each kde.
    bw_adjust : float
        Used to calculate the kde. Adjusts multiplicatively the bandwidth that is found automatically by
        scipy.stats.gaussian.kde. Larger values generate smoother estimations.
    title : string
        Title of the histogram. If none is given, a title will be generated.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    width : float
        Width in inches
    height : float
        Height in inches
    font_scale : float
        Scale the fontsize for the title, labels and ticks.
    save : bool
        Used for interactive plotting. If the save button is pressed (=True) then store to the given file path.
    latex : bool
        Use latex names for any title or axis. Otherwise use the
        code names of variables and such.

    Returns
    -------
    matplotlib.figure.Figure with the plot drawn onto it.
    """
    phase_colors = {
        "warm phase": "tab:orange",
        "mixed phase": "tab:green",
        "ice phase": "tab:blue",
        "any": "k",
    }

    ds_tmp = ds.sel({"Output Parameter": out_param}).drop("Output Parameter")
    common_norm = False  # All densities sum up to one if true
    df = ds_tmp[in_params].to_dataframe().stack().reset_index()
    df = df.rename(columns={"level_4": "Parameter", 0: "Rank"})
    df = df.loc[df["Rank"] > 0]
    worst_rank = ds[in_params[0]].attrs["rank for zero gradients"]
    if ignore_zero_gradients:
        df = df.loc[df["Rank"] < worst_rank]

    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    fig = Figure()
    ax = fig.subplots()
    if title is None:
        title = f"Histogram for impacts on {latexify.parse_word(out_param)}"

    if phase and flow:
        df = df.loc[(df["phase"] != "any") & (df["phase"] != "neutral phase")]
        new_labels = []
        phase_hue_order = np.sort(np.unique(df["phase"]))
        df_tmp = df.loc[df["flow"] == "inflow"]
        _ = sns.kdeplot(
            data=df_tmp,
            x="Rank",
            hue="phase",
            hue_order=phase_hue_order,
            common_norm=common_norm,
            linestyle="dotted",
            palette=phase_colors,
            ax=ax,
            label="inflow",
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            clip=(1, worst_rank),
        )
        for p in phase_hue_order[::-1]:
            if (
                p in df_tmp["phase"].values
                and len(np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])) >= 1
            ):
                new_labels.append(f"inflow {p}")
            # in case there is anything with zero variance, we plot a single bar.
            vals = np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])
            if len(vals) == 1:
                ax.axvline(
                    x=vals[0],
                    color=phase_colors[p],
                    linestyle="dotted",
                    label=f"inflow {p}",
                )
        df_tmp = df.loc[df["flow"] == "ascent"]
        _ = sns.kdeplot(
            data=df_tmp,
            x="Rank",
            hue="phase",
            hue_order=phase_hue_order,
            common_norm=common_norm,
            linestyle="solid",
            palette=phase_colors,
            ax=ax,
            label="ascent",
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            clip=(1, worst_rank),
        )
        for p in phase_hue_order[::-1]:
            if (
                p in df_tmp["phase"].values
                and len(np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])) >= 1
            ):
                new_labels.append(f"ascent {p}")
            # in case there is anything with zero variance, we plot a single bar.
            vals = np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])
            if len(vals) == 1:
                ax.axvline(
                    x=vals[0],
                    color=phase_colors[p],
                    linestyle="solid",
                    label=f"ascent {p}",
                )
        df_tmp = df.loc[df["flow"] == "outflow"]
        _ = sns.kdeplot(
            data=df_tmp,
            x="Rank",
            hue="phase",
            hue_order=phase_hue_order,
            common_norm=common_norm,
            linestyle="dashed",
            palette=phase_colors,
            ax=ax,
            label="outflow",
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            clip=(1, worst_rank),
        )
        for p in phase_hue_order[::-1]:
            if (
                p in df_tmp["phase"].values
                and len(np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])) >= 1
            ):
                new_labels.append(f"outflow {p}")
            # in case there is anything with zero variance, we plot a single bar.
            vals = np.unique(df_tmp.loc[df_tmp["phase"] == p]["Rank"])
            if len(vals) == 1:
                ax.axvline(
                    x=vals[0],
                    color=phase_colors[p],
                    linestyle="dashed",
                    label=f"outflow {p}",
                )
        ax.get_legend().remove()
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, new_labels, fontsize=int(10 * font_scale))
    elif phase:
        df = df.loc[
            (df["phase"] != "any")
            & (df["phase"] != "neutral phase")
            & (df["flow"] == "any")
        ]
        phase_hue_order = np.sort(np.unique(df["phase"]))
        _ = sns.kdeplot(
            data=df,
            x="Rank",
            hue="phase",
            hue_order=phase_hue_order,
            common_norm=common_norm,
            linestyle="solid",
            palette=phase_colors,
            ax=ax,
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            clip=(1, worst_rank),
        )
        # in case there is anything with zero variance, we plot a single bar.
        for phase in np.unique(df["phase"]):
            vals = np.unique(df.loc[df["phase"] == phase]["Rank"])
            if len(vals) == 1:
                ax.axvline(x=vals[0], color=phase_colors[phase])

    elif flow:
        df = df.loc[df["phase"] == "any"]
        param_hue_order = np.sort(np.unique(df["Parameter"]))
        new_labels = []
        df_tmp = df.loc[df["flow"] == "inflow"]
        _ = sns.kdeplot(
            data=df_tmp,
            x="Rank",
            hue="Parameter",
            hue_order=param_hue_order,
            common_norm=common_norm,
            linestyle="dotted",
            ax=ax,
            label="inflow",
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            clip=(1, worst_rank),
        )
        for param in param_hue_order[::-1]:
            if (
                param in df_tmp["Parameter"].values
                and len(np.unique(df_tmp.loc[df_tmp["Parameter"] == param]["Rank"])) > 1
            ):
                new_labels.append(f"inflow {param}")
        df_tmp = df.loc[df["flow"] == "ascent"]
        _ = sns.kdeplot(
            data=df_tmp,
            x="Rank",
            hue="Parameter",
            hue_order=param_hue_order,
            common_norm=common_norm,
            linestyle="solid",
            ax=ax,
            label="ascent",
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            clip=(1, worst_rank),
        )
        for param in param_hue_order[::-1]:
            if (
                param in df_tmp["Parameter"].values
                and len(np.unique(df_tmp.loc[df_tmp["Parameter"] == param]["Rank"])) > 1
            ):
                new_labels.append(f"ascent {param}")
        df_tmp = df.loc[df["flow"] == "outflow"]
        _ = sns.kdeplot(
            data=df_tmp,
            x="Rank",
            hue="Parameter",
            hue_order=param_hue_order,
            common_norm=common_norm,
            linestyle="dashed",
            ax=ax,
            label="outflow",
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            clip=(1, worst_rank),
        )
        for param in param_hue_order[::-1]:
            if (
                param in df_tmp["Parameter"].values
                and len(np.unique(df_tmp.loc[df_tmp["Parameter"] == param]["Rank"])) > 1
            ):
                new_labels.append(f"outflow {param}")
        ax.get_legend().remove()
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles, new_labels, fontsize=int(10 * font_scale))

    else:
        param_hue_order = np.sort(np.unique(df["Parameter"]))
        df = df.loc[df["phase"] == "any"]
        _ = sns.kdeplot(
            data=df.loc[df["flow"] == "any"],
            x="Rank",
            hue="Parameter",
            hue_order=param_hue_order,
            common_norm=common_norm,
            linestyle="solid",
            ax=ax,
            linewidth=linewidth,
            bw_adjust=bw_adjust,
            clip=(1, worst_rank),
        )

    if font_scale is None:
        _ = ax.set_title(title)
    else:
        ax.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
    ax.set_xlabel("Rank", fontsize=int(11 * font_scale))
    ax.set_ylabel("Density", fontsize=int(11 * font_scale))

    plt.tight_layout()
    if filename is not None and save:
        fig = ax.get_figure()
        try:
            i = 0
            store_type = filename.split(".")[-1]
            store_path = filename[0 : -len(store_type) - 1]
            save_name = store_path + "_{:03d}.".format(i) + store_type

            while os.path.isfile(save_name):
                i = i + 1
                save_name = store_path + "_{:03d}.".format(i) + store_type
            fig.savefig(save_name, bbox_inches="tight", dpi=300)
        except:
            print(f"Storing to {save_name} failed.", file=sys.stderr)
    return fig


def get_corr_matrix(
    ds,
    store_path=None,
    verbose=False,
):
    """

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with ranks per trajectory created by create_rank_traj_dataset().
    store_path : string
        If a store path is given then dumps the correlation matrix to f'{store_path}_correlation_matrix_per_traj.pkl'.
    verbose : bool
        Print progressbars.

    Returns
    -------
    Dictionary with output parameter, flow, and phase as keys and values are correlation matrices. Also a list of
    model parameter names for each column/row.
    """
    n = len(ds)

    corr_matrix = {
        out_p.item(): {
            flow.item(): {phase.item(): np.zeros((n, n)) for phase in ds["phase"]}
            for flow in ds["flow"]
        }
        for out_p in ds["Output Parameter"]
    }
    col_names = list(ds.keys())
    for out_p in tqdm(ds["Output Parameter"]) if verbose else ds["Output Parameter"]:
        ds_tmp = ds.sel({"Output Parameter": out_p})
        for flow in tqdm(ds["flow"], leave=False) if verbose else ds["flow"]:
            ds_tmp2 = ds_tmp.sel({"flow": flow})
            for phase in tqdm(ds["phase"], leave=False) if verbose else ds["phase"]:
                ds_tmp3 = ds_tmp2.sel({"phase": phase})
                for i, col in enumerate(
                    tqdm(col_names, leave=False) if verbose else col_names
                ):
                    for j, col2 in enumerate(col_names):
                        corr_matrix[out_p.item()][flow.item()][phase.item()][
                            i, j
                        ] = xr.cov(ds_tmp3[col], ds_tmp3[col2]).item()
    if store_path is not None:
        corr_and_names = {"correlation": corr_matrix, "column names": col_names}
        with open(store_path + "_correlation_matrix_per_traj.pkl", "wb") as f:
            pickle.dump(corr_and_names, f)
    return corr_matrix, col_names


def plot_kde_histogram_interactive(ds):
    """
    Use plot_kde_histogram() interactively in a jupyter notebook.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset created via create_rank_traj_dataset().

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=ds["Output Parameter"].values[0],
        options=ds["Output Parameter"].values.tolist(),
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
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    flow_button = pn.widgets.Toggle(
        name="Show flow",
        button_type="success",
    )
    phase_button = pn.widgets.Toggle(
        name="Show phase",
        button_type="success",
    )
    in_param_values = []
    for col in ds:
        if "rank" in col:
            in_param_values.append(col)
    in_param_values = list(np.sort(in_param_values))
    in_params = pn.widgets.CrossSelector(
        name="Parameter",
        value=in_param_values[0:2],
        options=in_param_values,
    )
    line_slider = pn.widgets.FloatSlider(
        name="Change the line width",
        start=1,
        end=10,
        step=0.5,
        value=2,
    )
    bw_slider = pn.widgets.FloatSlider(
        name="Change the bandwidth for the kde calculation",
        start=0.05,
        end=1.5,
        step=0.05,
        value=1.0,
    )
    ignore_zero_button = pn.widgets.Toggle(
        name="Ignore zero gradients",
        button_type="success",
    )

    plot_pane = pn.panel(
        pn.bind(
            plot_kde_histogram,
            ds=ds,
            in_params=in_params,
            out_param=out_param,
            linewidth=line_slider,
            bw_adjust=bw_slider,
            flow=flow_button,
            phase=phase_button,
            ignore_zero_gradients=ignore_zero_button,
            title=title_widget,
            filename=save_to_field,
            width=width_slider,
            height=height_slider,
            font_scale=font_slider,
            save=save_button,
            latex=latex_button,
        ),
    ).servable()

    return pn.Column(
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            latex_button,
        ),
        pn.Row(
            flow_button,
            phase_button,
            ignore_zero_button,
        ),
        pn.Row(
            in_params,
            out_param,
        ),
        pn.Row(
            line_slider,
            bw_slider,
            title_widget,
        ),
        plot_pane,
    )
