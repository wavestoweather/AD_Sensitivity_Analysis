"""Create latex tables with ranks of parameters.

"""
import numpy as np

from ad_sensitivity_analysis.plot.latexify import parse_word


# pylint: disable=too-many-branches, too-many-lines, too-many-statements, too-many-locals
def __rank_latex_per_output(
    ds,
    phase,
    flow,
    multiple_tables,
    table_text,
    param_order,
    all_params,
    avg,
    caption,
):
    """

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
    multiple_tables : bool
        If False, create a single large latex table for all model state variables. Otherwise, create a
        dictionary with model state variables as keys and latex tables as values.
    table_text : string
        Header of latex table.
    param_order : list of strings
        List of output parameters that defines the order of the latex table.
    all_params : bool
        If True, return all model parameters. Otherwise, dismiss those that have a median and mean rank
        over 10.
    avg : bool
        If True, use mean as column, otherwise use median.
    caption : string
        Caption for the table.

    Returns
    -------
    String for a latex table or dictionary of model state variables with latex tables as values.
    """
    ds2 = ds.sel({"phase": phase, "flow": flow})
    table_texts = {}
    for o_p in ds["Output Parameter"]:
        ds_op = ds2.sel({"Output Parameter": o_p})
        if multiple_tables:
            table_texts[o_p.item()] = table_text
        else:
            table_text += f"        {parse_word(o_p.item()):<25}& "
        p_i = 0
        if param_order is not None:
            rows = {}
        for param in ds:
            if "rank" not in param:
                continue
            ds_op2 = ds_op[param]
            ds_op2 = ds_op2.where(ds_op2 > 0)
            med = ds_op2.median().item()
            if med > 10 and ds_op2.mean().item() > 10 and not all_params:
                continue
            if np.isnan(med):
                continue
            table_row = ""
            if multiple_tables:
                table_row += "               "
            else:
                if p_i > 0:
                    table_row += "                                 & "
                p_i += 1
            word = parse_word(param[:-5]).replace(r"\partial ", "")
            table_row += f"{word:<40}& "
            if avg:
                ds_op3 = ds_op[param[:-4] + "avg"]
                table_row += (
                    f"{int(med):2d} & {ds_op3.mean().item():1.2e} "
                    f"& {ds_op2.mean().item():2.2f} & {ds_op2.std().item():2.2f}"
                    + r" \\"
                    + "\n"
                )
            else:
                table_row += (
                    f"{int(med):2d} & {ds_op2.mean().item():2.2f} & {ds_op2.std().item():2.2f}"
                    + r" \\"
                    + "\n"
                )

            if param_order is not None:
                rows[param] = table_row
            elif multiple_tables:
                table_texts[o_p.item()] += table_row
            else:
                table_text += table_row
        if param_order is not None:
            for param in param_order:
                if f"{param} rank" not in rows:
                    continue
                if multiple_tables:
                    table_texts[o_p.item()] += rows[f"{param} rank"]
                else:
                    table_text += rows[f"{param} rank"]
        if multiple_tables:
            table_texts[
                o_p.item()
            ] += r"""
        \end{tabular}
        \caption{"""
            table_texts[o_p.item()] += f"{caption}"
            table_texts[
                o_p.item()
            ] += r"""}
        \label{tab:analysis:rank_count}
    \end{table}
    """
    return table_texts


def create_rank_latex_table(
    ds,
    phase,
    flow,
    avg=False,
    all_params=False,
    multiple_tables=False,
    param_order=None,
):
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
    avg : bool
        If True, use mean as column, otherwise use median.
    all_params : bool
        If True, return all model parameters. Otherwise, dismiss those that have a median and mean rank
        over 10.
    multiple_tables : bool
        If False, create a single large latex table for all model state variables. Otherwise, create a
        dictionary with model state variables as keys and latex tables as values.
    param_order : list of strings
        List of output parameters that defines the order of the latex table.

    Returns
    -------
    String for a latex table or dictionary of model state variables with latex tables as values.
    """
    if multiple_tables:
        if avg:
            table_text = r"""
        \begin{table}[ht]
            \centering
            \begin{tabular}{l|c|c|c|c}
                \textbf{{Parameter}} & \textbf{Median Rank} & \textbf{Mean Impact} & \textbf{Mean Rank} & \textbf{STD} \\ \hline 
            """
        else:
            table_text = r"""
        \begin{table}[ht]
            \centering
            \begin{tabular}{l|c|c|c}
                \textbf{{Parameter}} & \textbf{Median Rank} & \textbf{Mean Rank} & \textbf{STD} \\ \hline 
            """
    else:
        if avg:
            table_text = r"""
            \begin{table}[ht]
                \centering
                \begin{tabular}{l|l|c|c|c|c}
                    
"""
            table_text += (
                r"\textbf{Model State Variable} & \textbf{{Parameter}} "
                r"& \textbf{Median Rank} & \textbf{Mean Impact} "
                r"& \textbf{Mean Rank} & \textbf{STD} \\ \hline "
            )
            table_text += "\n                    "
        else:
            table_text = r"""
            \begin{table}[ht]
                \centering
                \begin{tabular}{l|l|c|c|c}
                    """
            table_text += (
                r"\textbf{Model State Variable} & \textbf{{Parameter}} "
                r"& \textbf{Median Rank} & \textbf{Mean Rank} & \textbf{STD} \\ \hline "
            )
            table_text += "\n            "
    if all_params:
        caption = "The median and median absolute deviation of the rank if we calculate the rank over all trajectories."
    else:
        caption = (
            "The median and median absolute deviation of the rank "
            "if we calculate the rank over all trajectories. "
            "Only those with median or mean lower than 11 are used."
        )
    if phase == "any":
        caption += "Using any phase "
    else:
        caption += f"Using {phase} "
    if flow == "any":
        caption += "and inflow, outflow, and ascent."
    else:
        caption += f"and {flow}."
    caption += f"Using {phase} phase and {flow} flow."

    table_texts = __rank_latex_per_output(
        ds=ds,
        phase=phase,
        flow=flow,
        multiple_tables=multiple_tables,
        table_text=table_text,
        param_order=param_order,
        all_params=all_params,
        avg=avg,
        caption=caption,
    )
    if not multiple_tables:
        table_text += r"""
            \end{tabular}
            \caption{"""
        table_text += f"{caption}"
        table_text += r"""}
            \label{tab:analysis:rank_count}
        \end{table}
        """
        return table_text
    return table_texts


# pylint: disable=too-many-nested-blocks, too-many-statements, too-many-locals
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
            word = parse_word(param[:-5]).replace(r"\partial ", "")
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
        caption = (
            "The median and median absolute deviation of the rank "
            "of each model parameter if we calculate the rank over all trajectories. "
        )
        caption += "Only those with median or mean lower than 11 are used. Distinction with different phases. "
        caption += "wp = warm phase, mp = mixed phase, ip = ice phase. "
        caption += f"Impact on {parse_word(o_p.item())} is used here."
        table_text += f"{caption}"
        table_text += r"""}
            \label{tab:analysis:rank_count_variations_"""
        table_text += f"{o_p.item()}"
        table_text += """}
\\end{table}
        """
        tables[o_p.item()] = table_text
    return tables
