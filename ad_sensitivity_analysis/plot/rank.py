"""Visualize rank and impact.

"""
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd
import seaborn as sns

from ad_sensitivity_analysis.plot.aux_functions import save_plot
from ad_sensitivity_analysis.plot.colors import set_top_param_cbar, get_b8_colors


def _create_rank_impact_dataframe(
    ds, ds_imp, ds_rank, phase, flow, out_param, rank_name
):
    """

    Parameters
    ----------
    ds
    ds_imp
    ds_rank
    phase
    flow
    out_param
    rank_name

    Returns
    -------

    """
    if not phase and not flow:
        if out_param != "all":
            imp_vals = ds_imp.to_array().values.flatten()
            rank_vals = ds_rank.to_array().values.flatten()
            df = pd.DataFrame.from_dict(
                {
                    rank_name: rank_vals,
                    "Mean Impact": imp_vals,
                }
            )
        else:
            data_dic = {
                rank_name: np.asarray([]),
                "Mean Impact": np.asarray([]),
                "Output Parameter": np.asarray([]),
            }
            for out_p in ds["Output Parameter"]:
                imp_vals = (
                    ds_imp.sel({"Output Parameter": out_p}).to_array().values.flatten()
                )
                rank_vals = (
                    ds_rank.sel({"Output Parameter": out_p}).to_array().values.flatten()
                )
                data_dic[rank_name] = np.append(data_dic[rank_name], rank_vals)
                data_dic["Mean Impact"] = np.append(data_dic["Mean Impact"], imp_vals)
                data_dic["Output Parameter"] = np.append(
                    data_dic["Output Parameter"],
                    np.repeat(out_p.item(), len(rank_vals)),
                )
            df = pd.DataFrame.from_dict(data_dic)
    elif phase and flow:
        data_dic = {
            rank_name: np.asarray([]),
            "Mean Impact": np.asarray([]),
            "Flow": np.asarray([]),
        }
        for flow_v, phase_v in itertools.product(ds["flow"], ds["phase"]):
            imp_vals = (
                ds_imp.sel({"flow": flow_v, "phase": phase_v})
                .to_array()
                .values.flatten()
            )
            rank_vals = (
                ds_rank.sel({"flow": flow_v, "phase": phase_v})
                .to_array()
                .values.flatten()
            )
            data_dic[rank_name] = np.append(data_dic[rank_name], rank_vals)
            data_dic["Mean Impact"] = np.append(data_dic["Mean Impact"], imp_vals)
            data_dic["Flow"] = np.append(
                data_dic["Flow"], np.repeat(flow_v.item(), len(imp_vals))
            )
            data_dic["Phase"] = np.append(
                data_dic["Phase"], np.repeat(phase_v.item(), len(imp_vals))
            )
        df = pd.DataFrame.from_dict(data_dic)
    elif phase:
        data_dic = {
            rank_name: np.asarray([]),
            "Mean Impact": np.asarray([]),
            "Phase": np.asarray([]),
        }
        if out_param != "all":
            for phase_v in ds["phase"]:
                imp_vals = ds_imp.sel({"phase": phase_v}).to_array().values.flatten()
                rank_vals = ds_rank.sel({"phase": phase_v}).to_array().values.flatten()
                data_dic[rank_name] = np.append(data_dic[rank_name], rank_vals)
                data_dic["Mean Impact"] = np.append(data_dic["Mean Impact"], imp_vals)
                data_dic["Phase"] = np.append(
                    data_dic["Phase"], np.repeat(phase_v.item(), len(imp_vals))
                )
            df = pd.DataFrame.from_dict(data_dic)
        else:
            data_dic["Output Parameter"] = np.asarray([])
            for phase_v, out_p in itertools.product(
                ds["phase"], ds["Output Parameter"]
            ):
                imp_vals = (
                    ds_imp.sel({"Output Parameter": out_p, "phase": phase_v})
                    .to_array()
                    .values.flatten()
                )
                rank_vals = (
                    ds_rank.sel({"Output Parameter": out_p, "phase": phase_v})
                    .to_array()
                    .values.flatten()
                )
                data_dic[rank_name] = np.append(data_dic[rank_name], rank_vals)
                data_dic["Mean Impact"] = np.append(data_dic["Mean Impact"], imp_vals)
                data_dic["Phase"] = np.append(
                    data_dic["Phase"], np.repeat(phase_v.item(), len(imp_vals))
                )
                data_dic["Output Parameter"] = np.append(
                    data_dic["Output Parameter"],
                    np.repeat(out_p.item(), len(rank_vals)),
                )
            df = pd.DataFrame.from_dict(data_dic)
    elif flow:
        data_dic = {
            rank_name: np.asarray([]),
            "Mean Impact": np.asarray([]),
            "Flow": np.asarray([]),
        }
        for flow_v in ds["flow"]:
            imp_vals = ds_imp.sel({"flow": flow_v}).to_array().values.flatten()
            rank_vals = ds_rank.sel({"flow": flow_v}).to_array().values.flatten()
            data_dic[rank_name] = np.append(data_dic[rank_name], rank_vals)
            data_dic["Mean Impact"] = np.append(data_dic["Mean Impact"], imp_vals)
            data_dic["Flow"] = np.append(
                data_dic["Flow"], np.repeat(flow_v.item(), len(imp_vals))
            )
        df = pd.DataFrame.from_dict(data_dic)

    return df


def _scatter_rank_over_impact(
    df, rank_name, out_param, ax, dot_size, out_markers, colorblind=True
):
    """

    Parameters
    ----------
    df
    rank_name
    out_param
    ax
    dot_size
    out_markers
    colorblind

    Returns
    -------

    """
    param_color_order, color_shades = get_b8_colors(colorblind=colorblind)
    if out_param != "all":
        sns.scatterplot(
            data=df,
            x=rank_name,
            y="Impact",
            ax=ax,
            s=dot_size,
            hue="Parameter",
            hue_order=param_color_order,
            palette=color_shades,
            linewidth=0.7,
            legend=False,
        )
    else:
        sns.scatterplot(
            data=df,
            x=rank_name,
            y="Mean Impact",
            style="Output Parameter",
            ax=ax,
            s=dot_size,
            markers=out_markers,
            hue="Parameter",
            hue_order=param_color_order,
            palette=color_shades,
            linewidth=0.7,
            legend="full",
        )
        handles, labels = ax.get_legend_handles_labels()
        for handle in handles:
            handle._sizes = [dot_size]  # pylint: disable=protected-access
        ax.legend(handles=handles[-4:], labels=labels[-4:])


def _scatter_rank_over_impact_phase(
    df, rank_name, out_param, ax, dot_size, out_markers, phases, phase_colors
):
    """

    Parameters
    ----------
    df
    rank_name
    out_param
    ax
    dot_size
    out_markers
    phases
    phase_colors

    Returns
    -------

    """
    if out_param != "all":
        sns.scatterplot(
            data=df,
            x=rank_name,
            y="Mean Impact",
            style="Output Parameter",
            markers=out_markers,
            ax=ax,
            s=dot_size,
            hue="Phase",
            hue_order=phases,
            palette=phase_colors,
        )
    else:
        sns.scatterplot(
            data=df,
            x=rank_name,
            y="Mean Impact",
            style="Output Parameter",
            markers=out_markers,
            ax=ax,
            s=dot_size,
            hue="Phase",
            hue_order=phases,
            palette=phase_colors,
        )


def _scatter_rank_over_impact_flow(df, rank_name, ax, dot_size):
    """

    Parameters
    ----------
    df
    rank_name
    ax
    dot_size

    Returns
    -------

    """
    sns.scatterplot(
        data=df,
        x=rank_name,
        y="Mean Impact",
        ax=ax,
        s=dot_size,
        markers={
            "Inflow": "hexagon",
            "Ascent": "star",
            "Outflow": "plus (filled)",
        },
    )


def _scatter_rank_over_impact_phase_flow(
    df, rank_name, ax, dot_size, phases, phase_colors
):
    """

    Parameters
    ----------
    df
    rank_name
    ax
    dot_size
    phases
    phase_colors

    Returns
    -------

    """
    sns.scatterplot(
        data=df,
        x=rank_name,
        y="Mean Impact",
        ax=ax,
        s=dot_size,
        markers={
            "Inflow": "hexagon",
            "Ascent": "star",
            "Outflow": "plus (filled)",
        },
        hue="Phase",
        hue_order=phases,
        palette=phase_colors,
    )


def _get_ds_rank_imp_name(
    out_param, phase=False, flow=False, phases=None, flows=None, avg=False
):
    """

    Parameters
    ----------
    out_param
    phase
    flow
    phases
    flows
    avg

    Returns
    -------

    """
    if out_param != "all":
        ds = ds.sel({"Output Parameter": out_param})

    if not phase and not flow:
        ds = ds.sel({"phase": "any", "flow": "any"})
    elif phase:
        ds = ds.sel({"phase": phases, "flow": "any"})
    elif flow:
        ds = ds.sel({"phase": "any", "flow": flows})
    else:
        ds = ds.sel({"phase": phases, "flow": flows})
    ds_rank = ds[[p for p in ds if "rank" in p]]
    ds_imp = ds[[p for p in ds if "avg" in p and p != "avg ascent"]]
    ds_imp = ds_imp.mean(dim=["trajectory", "file"])
    if avg:
        ds_rank = ds_rank.mean(dim=["trajectory", "file"])
        rank_name = "Mean Rank"
    else:
        ds_rank = ds_rank.median(dim=["trajectory", "file"])
        rank_name = "Median Rank"
    return ds_imp, ds_rank, rank_name


# pylint: disable=too-many-arguments, too-many-locals
def plot_rank_over_impact(
    ds,
    out_param,
    width=24,
    height=12,
    font_scale=None,
    filename=None,
    title=None,
    save=False,
    latex=False,
    dot_size=6,
    log=False,
    phase=False,
    flow=False,
    avg=False,
    colorblind=False,
):
    """

    Parameters
    ----------
    ds
    out_param
    width
    height
    font_scale
    filename
    title
    save
    latex
    dot_size
    log
    phase
    flow
    avg
    colorblind

    Returns
    -------

    """
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    fig = Figure()
    # pylint: disable=no-member
    ax = fig.subplots()
    phases = ["warm phase", "mixed phase", "ice phase"]
    flows = ["inflow", "ascent", "outflow"]
    phase_colors = {
        "warm phase": "tab:orange",
        "mixed phase": "tab:green",
        "ice phase": "tab:blue",
        "any": "k",
    }
    if out_param != "all":
        out_markers = None
    else:
        out_markers = {
            "QV": "D",
            "latent_heat": "o",
            "latent_cool": "^",
        }
    ds_imp, ds_rank, rank_name = _get_ds_rank_imp_name(
        out_param=out_param,
        phase=phase,
        flow=flow,
        phases=phases,
        flows=flows,
        avg=avg,
    )

    df = _create_rank_impact_dataframe(
        ds=ds,
        ds_imp=ds_imp,
        ds_rank=ds_rank,
        phase=phase,
        flow=flow,
        out_param=out_param,
        rank_name=rank_name,
    )
    if not phase and not flow:
        _scatter_rank_over_impact(
            df=df,
            rank_name=rank_name,
            out_param=out_param,
            ax=ax,
            dot_size=dot_size,
            out_markers=out_markers,
            colorblind=colorblind,
        )
        set_top_param_cbar(
            fig=fig,
            ax=ax,
            font_scale=font_scale,
            colorblind=colorblind,
        )
    elif phase:
        _scatter_rank_over_impact_phase(
            df=df,
            rank_name=rank_name,
            out_param=out_param,
            ax=ax,
            dot_size=dot_size,
            out_markers=out_markers,
            phases=phases,
            phase_colors=phase_colors,
        )
    elif flow:
        _scatter_rank_over_impact_flow(
            df=df,
            rank_name=rank_name,
            ax=ax,
            dot_size=dot_size,
        )
    elif phase and flow:
        _scatter_rank_over_impact_phase_flow(
            df=df,
            rank_name=rank_name,
            ax=ax,
            dot_size=dot_size,
            phases=phases,
            phase_colors=phase_colors,
        )

    plt.setp(ax.get_legend().get_texts(), fontsize=int(9 * font_scale))
    plt.setp(ax.get_legend().get_title(), fontsize=int(11 * font_scale))
    if log:
        ax.set_yscale("log")
    if title is not None:
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
    ax.xaxis.get_label().set_fontsize(int(11 * font_scale))
    ax.yaxis.get_label().set_fontsize(int(11 * font_scale))
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=int(10 * font_scale),
    )
    ax.xaxis.get_offset_text().set_size(int(9 * font_scale))
    ax.yaxis.get_offset_text().set_size(int(9 * font_scale))

    if filename is not None and save:
        save_plot(filename, ax.get_figure())
    return fig


def _get_probs_top_order_vals(ds, rank=False):
    """

    Parameters
    ----------
    ds
    rank

    Returns
    -------

    """
    vals = {"QV": {}, "latent_heat": {}, "latent_cool": {}}
    median_vals = {"QV": {}, "latent_heat": {}, "latent_cool": {}}
    mean_vals = {"QV": {}, "latent_heat": {}, "latent_cool": {}}
    order = []
    for out_p in ds["Output Parameter"]:
        ds_tmp = ds.sel({"Output Parameter": out_p})
        for col in ds_tmp:
            if "rank" not in col:
                continue
            ds_tmp2 = ds_tmp.where(ds_tmp[col] > 0)
            vals[out_p.item()][col[:-5]] = (
                np.nansum(~np.isnan(ds_tmp2.where(ds_tmp2[col] <= 10)[col]))
                / np.nansum(~np.isnan(ds_tmp2[col]))
                * 100
            )
            if out_p.item() == ds["Output Parameter"].values[0]:
                order.append(col[:-5])
            if not rank:
                median_vals[out_p.item()][col[:-5]] = (
                    ds_tmp2[f"{col[:-5]} avg"].median().item()
                )
                mean_vals[out_p.item()][col[:-5]] = (
                    ds_tmp2[f"{col[:-5]} avg"].mean().item()
                )
            else:
                median_vals[out_p.item()][col[:-5]] = ds_tmp2[col].median().item()
                mean_vals[out_p.item()][col[:-5]] = ds_tmp2[col].mean().item()
    vals2 = {"QV": [], "latent_heat": [], "latent_cool": []}
    median = {"QV": [], "latent_heat": [], "latent_cool": []}
    mean = {"QV": [], "latent_heat": [], "latent_cool": []}

    for t in ds["Output Parameter"]:
        t = t.item()
        for o in order:
            vals2[t].append(vals[t][o])
            median[t].append(median_vals[t][o])
            mean[t].append(mean_vals[t][o])
    return order, vals2, median, mean


def _get_probs_top_params(order, vals, mark_top_n):
    """

    Parameters
    ----------
    order
    vals
    mark_top_n

    Returns
    -------

    """
    o_p = "QV"
    df = pd.DataFrame.from_dict(
        {
            "Probability for high ranking": vals[o_p],
            "Parameter": order,
        }
    )
    top_params1 = np.sort(
        np.unique(
            df.nlargest(n=mark_top_n, columns=["Probability for high ranking"])[
                "Parameter"
            ]
        )
    )

    o_p = "latent_heat"
    df = pd.DataFrame.from_dict(
        {
            "Probability for high ranking": vals[o_p],
            "Parameter": order,
        }
    )
    top_params2 = np.sort(
        np.unique(
            df.nlargest(n=mark_top_n, columns=["Probability for high ranking"])[
                "Parameter"
            ]
        )
    )

    o_p = "latent_cool"
    df = pd.DataFrame.from_dict(
        {
            "Probability for high ranking": vals[o_p],
            "Parameter": order,
        }
    )
    top_params3 = np.sort(
        np.unique(
            df.nlargest(n=mark_top_n, columns=["Probability for high ranking"])[
                "Parameter"
            ]
        )
    )
    top_params = np.sort(
        np.unique(list(top_params1) + list(top_params2) + list(top_params3))
    )
    top_params_order = {
        "QV": top_params1,
        "latent_heat": top_params2,
        "latent_cool": top_params3,
    }
    return top_params, top_params_order


def _plot_rank_probs_single(
    ds, ax, out_param, y, y2, rank, dot_size, mark_top_n, colorblind=True
):
    """

    Parameters
    ----------
    ds
    ax
    out_param
    y
    y2
    rank
    dot_size
    mark_top_n
    colorblind

    Returns
    -------

    """
    order, vals, median, mean = _get_probs_top_order_vals(ds=ds, rank=rank)
    param_color_order, color_shades = get_b8_colors(colorblind=colorblind)

    df = pd.DataFrame.from_dict(
        {
            "Probability for high ranking": vals[out_param],
            "Parameter": order,
            f"Median {y2}": median[out_param],
            f"Mean {y2}": mean[out_param],
        }
    )
    sns.scatterplot(
        data=df,
        x="Probability for high ranking",
        y=f"{y} {y2}",
        facecolor="k",
        ax=ax,
        s=dot_size,
        linewidth=0.7,
        legend=False,
    )
    sns.scatterplot(
        data=df.nlargest(n=mark_top_n, columns=["Probability for high ranking"]),
        x="Probability for high ranking",
        y=f"{y} {y2}",
        ax=ax,
        s=dot_size,
        hue="Parameter",
        palette=color_shades,
        hue_order=param_color_order,
        linewidth=0.7,
        legend=False,
    )


def _plot_rank_probs_multiple(
    ds, ax, y, y2, rank, dot_size, mark_top_n, font_scale, colorblind=True
):
    """

    Parameters
    ----------
    ds
    ax
    y
    y2
    rank
    dot_size
    mark_top_n
    font_scale
    colorblind

    Returns
    -------

    """
    order, vals, median, mean = _get_probs_top_order_vals(ds=ds, rank=rank)
    param_color_order, color_shades = get_b8_colors(colorblind=colorblind)
    out_markers = {
        "QV": "D",
        "latent_heat": "o",
        "latent_cool": "^",
    }
    params = list(out_markers.keys())
    probs = np.array([vals[p] for p in params])
    avgs = np.array([])
    if y == "Median":
        for p in params:
            avgs = np.append(avgs, median[p])
    else:
        for p in params:
            avgs = np.append(avgs, mean[p])
    out_names = np.array([])
    for p in params:
        out_names = np.append(out_names, np.repeat(p, len(mean[p])))
    df = pd.DataFrame.from_dict(
        {
            "Probability for high ranking": probs,
            "Parameter": np.tile(order, len(params)),
            f"{y} {y2}": avgs,
            "Output Parameter": out_names,
        }
    )
    sns.scatterplot(
        data=df,
        x="Probability for high ranking",
        y=f"{y} {y2}",
        facecolor="k",
        ax=ax,
        s=dot_size,
        style="Output Parameter",
        markers=out_markers,
    )
    df_tmp = None
    for p in params:
        df_tmp2 = df.loc[df["Output Parameter"] == p].nlargest(
            n=mark_top_n, columns=["Probability for high ranking"]
        )
        if df_tmp is not None:
            df_tmp = df_tmp.append(df_tmp2)
        else:
            df_tmp = df_tmp2
    sns.scatterplot(
        data=df_tmp,
        x="Probability for high ranking",
        y=f"{y} {y2}",
        ax=ax,
        s=dot_size,
        style="Output Parameter",
        markers=out_markers,
        palette=color_shades,
        hue_order=param_color_order,
        linewidth=0.7,
        hue="Parameter",
        legend="full",
    )
    handles, labels = ax.get_legend_handles_labels()
    for handle in handles:
        handle._sizes = [dot_size]  # pylint: disable=protected-access
    ax.legend(handles=handles[-4:], labels=labels[-4:])
    plt.setp(ax.get_legend().get_texts(), fontsize=int(9 * font_scale))
    plt.setp(ax.get_legend().get_title(), fontsize=int(11 * font_scale))


def plot_rank_probs(
    ds,
    flow,
    phase,
    out_param,
    median=False,
    rank=False,
    mark_top_n=10,
    dot_size=12,
    title=None,
    filename=None,
    width=17,
    height=16,
    font_scale=None,
    save=False,
    latex=False,
    logy=False,
    colorblind=True,
    show_cbar=False,
):
    """

    Parameters
    ----------
    ds
    flow
    phase
    out_param
    median
    rank
    mark_top_n
    dot_size
    title
    filename
    width
    height
    font_scale
    save
    latex
    logy
    colorblind
    show_cbar

    Returns
    -------

    """
    ds = ds.sel({"flow": flow, "phase": phase})
    y = "Median" if median else "Mean"
    y2 = "Rank" if rank else "Impact"
    sns.set(rc={"figure.figsize": (width, height), "text.usetex": latex})
    fig = Figure()
    # pylint: disable=no-member
    ax = fig.subplots()
    if out_param != "all":
        _plot_rank_probs_single(
            ds=ds,
            ax=ax,
            out_param=out_param,
            y=y,
            y2=y2,
            rank=rank,
            dot_size=dot_size,
            mark_top_n=mark_top_n,
            colorblind=colorblind,
        )
    else:
        _plot_rank_probs_multiple(
            ds=ds,
            ax=ax,
            y=y,
            y2=y2,
            rank=rank,
            dot_size=dot_size,
            mark_top_n=mark_top_n,
            font_scale=font_scale,
            colorblind=colorblind,
        )
    if show_cbar:
        set_top_param_cbar(
            fig=fig,
            ax=ax,
            font_scale=font_scale,
            colorblind=colorblind,
        )

    if title is not None:
        _ = ax.set_title(title, fontsize=int(12 * font_scale))
    if logy:
        ax.set_yscale("log")
    ax.tick_params(axis="both", which="major", labelsize=int(10 * font_scale))
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=int(10 * font_scale), markerscale=font_scale)
    ax.set_xlabel("Probability for high ranking", fontsize=int(11 * font_scale))
    ax.yaxis.get_offset_text().set_fontsize(int(11 * font_scale))
    ax.set_ylabel(f"{y} {y2}", fontsize=int(11 * font_scale))
    plt.tight_layout()
    if filename is not None and save:
        save_plot(filename, fig)
    return fig
