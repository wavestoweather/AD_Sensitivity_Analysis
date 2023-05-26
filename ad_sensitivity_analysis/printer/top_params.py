"""Print top parameters in different ways.

"""
import numpy as np
import pandas as pd

from ad_sensitivity_analysis.plot.latexify import (
    in_params_notation_mapping,
    in_params_grouping,
    parse_word,
)


def _top_list_to_table(top, caption, label, parse=True):
    """
    Create a latex table with model parameters and a description of those.

    Parameters
    ----------
    top : set or list
        A set of model parameter names for a table with descriptions.
    caption : string
        Caption to put in the table.
    label : string
        Label for the table.
    parse : bool
        Parse parameters for latex.

    Returns
    -------
        String that can be printed for a latex table.
    """
    n = 0
    for p in top:
        if parse:
            parsed = parse_word(p).replace(r"\partial ", "")
        else:
            parsed = p
        if len(parsed) > n:
            n = len(parsed)
    table = """\
\\bgroup
\\def\\arraystretch{1.3} 
\\begin{tabularx}{\\linewidth}{@{}l|X@{}}
"""
    table += f"\\caption{{{caption}}}\\\\ \\hline\n"
    tmp = "\\textbf{Model Parameter}"
    table += f"{tmp} & \\textbf{{Description}} \\\\[6pt] \\hline \n"
    for p in top:
        if parse:
            parsed = parse_word(p).replace(r"\partial ", "")
        else:
            parsed = p
        method = "".join([f"{mapped}, " for mapped in in_params_notation_mapping[p][4]])
        method = method.replace("_", "\\_")
        table += f"        {parsed:<{n}} & {in_params_notation_mapping[p][0]} ({method}) \\\\"
        table += "\n"
    table += f"    \\label{{{label}}}\n"
    table += "\\end{tabularx}\n"
    table += "\\egroup\n"
    return table


# pylint: disable=too-many-locals
def _top_dict_phase_to_table(top, caption, label, parse=True):
    """
    Create a latex table with model parameters for different model state variables and phases but also
    without description.

    Parameters
    ----------
    top : Dict
        A dictionary keys with phases and model state
        variable name for a distinction between different phases.
    caption : string
        Caption to put in the table.
    label : string
        Label for the table.
    parse : bool
        Parse parameters for latex.

    Returns
    -------
        String that can be printed for a latex table.
    """
    table = """\
\\begin{table}[ht]
    \\centering
    \\begin{tabular}{l|l|l}
        """
    tmp = "\\textbf{Model State Variable}"
    tmp2 = "\\textbf{Phase}"
    tmp3 = "\\textbf{{Top Parameters}}"
    n = len(tmp)
    n2 = len(tmp2)
    n3 = len(tmp3)
    for var in top:
        n = len(var.split("phase ")[1]) if len(var.split("phase ")[1]) > n else n
        phase = var.split("phase ")[0] + "phase"
        if len(phase) > n3:
            n3 = len(phase)
    table += f"{tmp:<{n}} & {tmp2:<{n2}} & {tmp3:<{n3}} \\\\ \n"
    for var in top:
        var1 = var.split("phase ")[1]
        if parse:
            table += f"        {parse_word(var1):<{n}} & "
        else:
            table += f"        {var1:<{n}} & "
        phase = var.split("phase ")[0] + "phase"
        table += f"{phase:<{n2}} & "
        line = "$"
        empty = " "
        n_breaks = 1
        for param in top[var]:
            if parse:
                line += f"{parse_word(param)}".replace(" $", "").replace("$", "")
            else:
                line += f"{param}"
            if len(line) > 100 * n_breaks and param != top[var][-1]:
                line += f"$ \\\\ \n{empty:<{n+8}} & {empty:<{n2}} & $ "
                n_breaks += 1.5
            else:
                line += ", "
        table += line + "$ \\\\ \n"
    table += "    \\end{tabular}\n"
    table += f"    \\caption{{{caption}}}\n"
    table += f"    \\label{{{label}}}\n"
    table += "\\end{table}"
    return table


def _top_dict_to_table(top, caption, label, parse=True):
    """
    Create either a latex table with model parameters and
    differentiate between sensitivities for different model state variables and
    print the top parameters without description.


    Parameters
    ----------
    top : Dict
        A dictionary with keys the model state variable for a table with top
        parameters for each model state variable.
    caption : string
        Caption to put in the table.
    label : string
        Label for the table.
    parse : bool
        Parse parameters for latex.

    Returns
    -------
        String that can be printed for a latex table.
    """
    table = """\
\\begin{table}[ht]
    \\centering
    \\begin{tabular}{l|l}
        """
    tmp = "\\textbf{Model State Variable}"
    tmp2 = "\\textbf{{Top Parameters}}"
    n = len(tmp)
    n2 = len(tmp2)
    for var in top:
        if parse:
            if len(parse_word(var)) > n:
                n = len(parse_word(var))
        else:
            if len(var) > n:
                n = len(var)

    table += f"{tmp:<{n}} & {tmp2:<{n2}} \\\\ \\hline \n"
    for var in top:
        if parse:
            table += f"        {parse_word(var):<{n}} & "
        else:
            table += f"        {var:<{n}} & "
        line = "$"
        empty = " "
        n_breaks = 1
        for param in top[var]:
            if parse:
                line += f"{parse_word(param)}".replace(" $", "").replace("$", "")
            else:
                line += f"{param}"
            if len(line) > 100 * n_breaks and param != top[var][-1]:
                line += f"$ \\\\ \n{empty:<{n+8}} &  $ "
                n_breaks += 1.5
            else:
                line += ", "
        table += line + "$ \\\\ \n"
    table += "    \\end{tabular}\n"
    table += f"    \\caption{{{caption}}}\n"
    table += f"    \\label{{{label}}}\n"
    table += "\\end{table}"
    return table


def top_to_table(top, caption, label, parse=True):
    """
    Create either a latex table with model parameters and a description of those.
    Or differentiate between sensitivities for different model state variables and
    print the top parameters without description.
    Or differentiate between different model state variables and phases but also
    without description.

    Parameters
    ----------
    top : Dict or set or list
        Either a set of model parameter names for a table with descriptions
        or a dictionary with keys the model state variable for a table with top
        parameters for each model state variable or a dictionary keys with phases and model state
        variable name for a distinction between different phases.
    caption : string
        Caption to put in the table.
    label : string
        Label for the table.
    parse : bool
        Parse parameters for latex.

    Returns
    -------
        String that can be printed for a latex table.
    """
    if isinstance(top, (list, set)):
        return _top_list_to_table(top=top, caption=caption, label=label, parse=parse)
    if isinstance(top, dict):
        if "phase" in list(top.keys())[0]:
            return _top_dict_phase_to_table(
                top=top, caption=caption, label=label, parse=parse
            )
        return _top_dict_to_table(top=top, caption=caption, label=label, parse=parse)
    return ""


def print_unique_params(top_sens_dic):
    """
    Print the parameters that appear only for a single output variable.

    Parameters
    ----------
    top_sens_dic : dict of list of strings
        The result of get_top_list() or get_magnitude_list()

    Returns
    -------
    String of all printed statements
    """
    local_text = "\nParameters that appear only for a single output variable\n"
    print(local_text)
    unique_list = []
    unique_pairing = []
    not_unique_list = []
    for out_p in top_sens_dic:
        tmp = top_sens_dic[out_p]
        for param in tmp:
            if param not in unique_list and param not in not_unique_list:
                unique_list.append(param)
                unique_pairing.append((out_p, param))
            elif param in unique_list:
                # find the index
                idx = np.argwhere(np.asarray(unique_list) == param)
                not_unique_list.append(param)
                del unique_list[idx[0][0]]
                del unique_pairing[idx[0][0]]
    for pair in unique_pairing:
        print(pair)
        local_text += f"{pair}\n"
    return local_text


def get_corr(df_tmp, kind):
    """

    Parameters
    ----------
    df_tmp
    kind

    Returns
    -------

    """
    return df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(kind)[
        "Predicted Squared Error"
    ][1]


def print_corr(df_tmp, text_tmp):
    """

    Parameters
    ----------
    df_tmp
    text_tmp

    Returns
    -------

    """
    spearman = get_corr(df_tmp, "spearman")
    pearson = get_corr(df_tmp, "pearson")
    kendall = get_corr(df_tmp, "kendall")
    text_tmp += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}\n"
    print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
    df_tmp = df_tmp.loc[df_tmp["Predicted Squared Error"] != 0]
    n = len(np.unique(df_tmp["Input Parameter"]))
    text_tmp += f"Correlation without zero parameters; total of {n} parameters"
    print(f"Correlation without zero parameters; total of {n} parameters")
    spearman = get_corr(df_tmp, "spearman")
    pearson = get_corr(df_tmp, "pearson")
    kendall = get_corr(df_tmp, "kendall")
    text_tmp += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}"
    print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
    return text_tmp


def print_correlation_broad(ds, out_params):
    """
    Print correlation coefficients (Spearman, Pearson, and Kendall) between predicted and actual error using each
    time step individually with all data for each type of state variable (first moment, second moment,
    sedimentation), and for each output variable.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    out_params : list-like of strings
        The model state variables for which sensitivities have been calculated for.

    Returns
    -------
    String of all printed statements
    """
    local_text = "\nCorrelation with all data\n"
    print(local_text)

    df = ds.to_dataframe().reset_index()
    local_text = print_corr(df, local_text)

    local_text += "\nFor each output variable type individually\n"
    print("\nFor each output variable type individually\n")

    second_moment = ["QV", "QC", "QR", "QS", "QG", "QH", "QI"]
    first_moment = ["NCCLOUD", "NCRAIN", "NCGRAUPEL", "NCHAIL", "NCICE", "NCSNOW"]
    second_sed = ["QR_OUT", "QG_OUT", "QH_OUT", "QI_OUT", "QS_OUT"]
    first_sed = ["NR_OUT", "NG_OUT", "NH_OUT", "NI_OUT", "NS_OUT"]

    local_text += "##################First Moment (Number Count)"
    print("##################First Moment (Number Count)")
    df = ds.sel({"Output Parameter": first_moment}).to_dataframe().reset_index()
    local_text = print_corr(df, local_text)

    local_text += "##################Second Moment (Mixing Ratio)"
    print("##################Second Moment (Mixing Ratio)")
    df = ds.sel({"Output Parameter": second_moment}).to_dataframe().reset_index()
    local_text = print_corr(df, local_text)

    local_text += "##################First Moment Sedimentation (Number Count)"
    print("##################First Moment Sedimentation (Number Count)")
    df = ds.sel({"Output Parameter": first_sed}).to_dataframe().reset_index()
    local_text = print_corr(df, local_text)

    local_text += "##################Second Moment Sedimentation (Mixing Ratio)"
    print("##################Second Moment Sedimentation (Mixing Ratio)")
    df = ds.sel({"Output Parameter": second_sed}).to_dataframe().reset_index()
    local_text = print_corr(df, local_text)

    local_text += "\nFor each output variable individually\n"
    print("\nFor each output variable individually\n")

    for out_p in out_params:
        local_text += f"##################{out_p}"
        print(f"##################{out_p}")
        df = ds.sel({"Output Parameter": [out_p]}).to_dataframe().reset_index()
        local_text = print_corr(df, local_text)
    return local_text


def _print_correlation_mean_top_params(text, out_params, df):
    """
    Correlation taking different number of top parameters

    Parameters
    ----------
    text
    out_params
    df

    Returns
    -------

    """
    tuples = []
    for n in [10, 20, 50, 100, 150, 200]:
        text += f"##############n={n}#################"
        print(f"##############n={n}#################")
        param_list = []
        for out_p in out_params:
            df_tmp = df.loc[df["Output Parameter"] == out_p]
            if np.max(df_tmp["Predicted Squared Error"]) == 0:
                continue
            param_list.extend(
                df_tmp.nlargest(n, "Predicted Squared Error")["Input Parameter"].values
            )
        param_set = set(param_list)
        df_tmp = df.loc[df["Input Parameter"].isin(param_set)]
        i = len(np.unique(df_tmp["Input Parameter"]))
        text += f"Correlation with zero parameters; total of {i} parameters\n"
        df_tmp2 = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "spearman"
        )["Predicted Squared Error"][1]
        text += f"{df_tmp2}"
        print(f"Correlation with zero parameters; total of {i} parameters")
        print(df_tmp2)
        pearson = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "pearson"
        )["Predicted Squared Error"][1]
        kendall = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "kendall"
        )["Predicted Squared Error"][1]
        text += f"Pearson: {pearson}, Kendall: {kendall}"
        print(f"Pearson: {pearson}, Kendall: {kendall}")
        tuples.append(
            (
                n,
                df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
                    "spearman"
                )["Predicted Squared Error"][1],
            )
        )
    for t in tuples:
        text += f"n={t[0]}, r={t[1]:1.3f}"
        print(f"n={t[0]}, r={t[1]:1.3f}")
    return text


def _print_correlation_mean_quantile(text, out_params, ds):
    """
    Correlation only for the 75th percentile

    Parameters
    ----------
    text
    out_params
    ds

    Returns
    -------

    """
    params = []
    for out_p in out_params:
        if out_p in ("NH_OUT", "QH_OUT"):
            continue
        text += f"##################{out_p}\n"
        print(f"##################{out_p}")
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        quan = df["Predicted Squared Error"].quantile(q=0.75)
        df = df.loc[df["Predicted Squared Error"] >= quan]
        n = len(np.unique(df["Input Parameter"]))
        params.extend(np.unique(df["Input Parameter"]))
        text += f"Correlation only with 75th percentile; total of {n} parameters\n"
        print(f"Correlation only with 75th percentile; total of {n} parameters")
        r_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
            "Predicted Squared Error"
        ][1]
        text += f"{r_2}"
        print(r_2)
        pearson = df[["Predicted Squared Error", "Mean Squared Error"]].corr("pearson")[
            "Predicted Squared Error"
        ][1]
        kendall = df[["Predicted Squared Error", "Mean Squared Error"]].corr("kendall")[
            "Predicted Squared Error"
        ][1]
        text += f"Pearson: {pearson}, Kendall: {kendall}"
        print(f"Pearson: {pearson}, Kendall: {kendall}")
    return text, list(set(params))


def _print_correlation_mean_all_spearman(text, out_params, ds, params):
    """

    Parameters
    ----------
    text
    out_params
    ds
    params

    Returns
    -------

    """
    df_tmp = (
        ds.sel({"Input Parameter": params})
        .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    r_spearman = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
        "spearman"
    )["Predicted Squared Error"][1]
    text += f"Correlation for all predicted errors with ensemble errors: {r_spearman}"
    print(f"Correlation for all predicted errors with ensemble errors: {r_spearman}")
    text += "\nCreate Latex table with Spearman correlations\n"
    print("\nCreate Latex table with Spearman correlations\n")
    table_corr = r"""
\begin{table}[hbt]
    \centering
    \begin{tabular}{l|c|c}
        Model State Variable $y_s$                 & $r(y_s)$ without zero sensitivities & $r(y_s)$ \\ \hline
"""
    for out_p in out_params:
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        r_spear = df[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "spearman"
        )["Predicted Squared Error"][1]
        df = df.loc[df["Predicted Squared Error"] != 0]
        r_spear_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "spearman"
        )["Predicted Squared Error"][1]
        table_corr += (
            "\t\t"
            + parse_word(out_p).title().replace(" Of ", " of ")
            + "\t\t& $ "
            + f"{r_spear_2:1.3f}"
            + " $ & $ "
            + f"{r_spear:1.3f} $ "
            + r"\\"
            + "\n"
        )
    table_corr += r"""    \end{tabular}
    \caption{}
    \label{tab:validate:correlation}
\end{table}"""
    print(table_corr)
    text += table_corr
    return text


def _print_correlation_mean_all_pearson(text, out_params, ds):
    """

    Parameters
    ----------
    text
    out_params
    ds

    Returns
    -------

    """
    text += "\nCreate Latex table with Pearson correlations\n"
    print("\nCreate Latex table with Pearson correlations\n")
    table_corr = r"""
\begin{table}[hbt]
    \centering
    \begin{tabular}{l|c|c}
        Model State Variable $y_s$                 
"""
    table_corr += r"& $r_{\text{pearson}}(y_s)$ without zero sensitivities & $r_{\text{pearson}}(y_s)$ \\ \hline"
    for out_p in out_params:
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        r_pearson = df[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "pearson"
        )["Predicted Squared Error"][1]
        df = df.loc[df["Predicted Squared Error"] != 0]
        r_pearson_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "pearson"
        )["Predicted Squared Error"][1]
        table_corr += (
            "\t\t"
            + parse_word(out_p).title().replace(" Of ", " of ")
            + "\t\t& $ "
            + f"{r_pearson_2:1.3f}"
            + " $ & $ "
            + f"{r_pearson:1.3f} $ "
            + r"\\"
            + "\n"
        )
    table_corr += r"""    \end{tabular}
    \caption{}
    \label{tab:correlation_pearson}
\end{table}"""
    print(table_corr)
    text += table_corr
    return text


def print_correlation_mean(ds, out_params):
    """
    Print correlation coefficients (Spearman, Pearson, and Kendall) between predicted and actual error using the mean
    over time after ascent and trajectory for each type of state variable (first moment, second moment,
    sedimentation), and for each output variable.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    out_params : list-like of strings
        The model state variables for which sensitivities have been calculated for.

    Returns
    -------
    String of all printed statements
    """
    text = "\nInstead of looking at correlations within time steps and trajectories individually, "
    text += "take the mean of those and look at the correlation\n"
    text += "\nCorrelation with all data\n"
    print(text)

    df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    text = print_corr(df, text)

    text += "\nFor each output variable individually\n"
    print("\nFor each output variable individually\n")
    for out_p in out_params:
        text += f"##################{out_p}"
        print(f"##################{out_p}")
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        text = print_corr(df, text)

    text += "\nCorrelation taking different number of top parameters\n"
    print("\nCorrelation taking different number of top parameters\n")
    df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    text = _print_correlation_mean_top_params(text, out_params, df)

    text += "\nCorrelation only for the 75th percentile\n"
    print("\nCorrelation only for the 75th percentile\n")
    text, params = _print_correlation_mean_quantile(text, out_params, ds)
    text += f"Number of parameters: {len(params)}"
    print(f"Number of parameters: {len(params)}")
    text = _print_correlation_mean_all_spearman(text, out_params, ds, params)
    text = _print_correlation_mean_all_pearson(text, out_params, ds)
    return text


def print_variable_with_important_params(sort_key_list):
    """
    Print model state variables and their number of model parameters. The number of parameters is determined
    by associating a model parameter x with a model state variable y if the sensitivity dx/dy is highest for y.

    Parameters
    ----------
    sort_key_list
        A sorted list of (predicted squared error, model parameter, model state variable,
        string of a row for model parameters in latex) which is sorted by the name of the model state variable.
        You may generate this using print_late_tables().

    Returns
    -------
    String of all printed statements
    """
    text = (
        "\nWhich model state variable has the most parameters with a high influence?\n"
    )
    print(text)
    state_counts = {}
    for _, _, state_variable, _ in sort_key_list:
        if state_variable not in state_counts:
            state_counts[state_variable] = 1
        else:
            state_counts[state_variable] += 1
    for state_variable, counts in state_counts.items():
        print(f"{state_variable}: {counts}")
        text += f"{state_variable}: {counts}"
    return text


# pylint: disable=too-many-branches
def print_param_types(ds, table_dic):
    """
    Print the type of parameters that are available in the dataset. Types are 'physical',
    'physical (high variability)', 'artificial', and 'artificial (threshold)'. The categories had been determined
    using the feedback of several meteorologists and are merely a guidance, not a definitive categorization.
    In addition, a categorization into geometric, velocity related, exponents, coefficients, and miscellaneous
    parameters is made too.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    table_dic : dict of strings
        Dictionary with keys = model parameters where the value is a string of a row of the latex table.
        This can be generated using print_latex_tables().

    Returns
    -------
    String of all printed statements
    """
    text = "\nHow many type of parameters are there?\n"
    print(text)
    counts = {
        "geo_count": 0,
        "vel_count": 0,
        "misc_count": 0,
        "exp_count": 0,
        "coeff_count": 0,
        "else_count": 0,
        "phys_count": 0,
        "phys_var_count": 0,
        "art_count": 0,
        "art_thresh_count": 0,
        "else_group_count": 0,
        "total_phys_count": 0,
        "total_phys_var_count": 0,
        "total_art_count": 0,
        "total_art_thresh_count": 0,
        "total_parameters": len(np.unique(ds["Input Parameter"])),
    }
    for key in table_dic:
        if "geo" in key:
            counts["geo_count"] += 1
        elif "vel" in key:
            counts["vel_count"] += 1
        else:
            counts["misc_count"] += 1

        if "exponent" in table_dic[key].lower():
            counts["exp_count"] += 1
        elif "coefficient" in table_dic[key].lower():
            counts["coeff_count"] += 1
        else:
            counts["else_count"] += 1

        if (
            "physical" in table_dic[key].lower()
            and not "high variability" in table_dic[key].lower()
        ):
            counts["phys_count"] += 1
        elif (
            "artificial" in table_dic[key].lower()
            and not "threshold" in table_dic[key].lower()
        ):
            counts["art_count"] += 1
        elif "physical" in table_dic[key].lower():
            counts["phys_var_count"] += 1
        elif "threshold" in table_dic[key].lower():
            counts["art_thresh_count"] += 1
        else:
            counts["else_group_count"] += 1

    for param in np.unique(ds["Input Parameter"]):
        if param in in_params_grouping["physical"]:
            counts["total_phys_count"] += 1
        if param in in_params_grouping["physical (high variability)"]:
            counts["total_phys_var_count"] += 1
        if param in in_params_grouping["artificial"]:
            counts["total_art_count"] += 1
        if param in in_params_grouping["artificial (threshold)"]:
            counts["total_art_thresh_count"] += 1

    local_text = f"There are {counts['total_parameters']} many parameters\n"
    local_text += f"There are {counts['geo_count']} geometric parameters\n"
    local_text += f"There are {counts['vel_count']} velocity parameters\n"
    local_text += f"There are {counts['misc_count']} misc. parameters\n"
    local_text += f"There are {counts['exp_count']} exponents\n"
    local_text += f"There are {counts['coeff_count']} coefficients\n"
    local_text += f"There are {counts['else_count']} not determined parameters in terms of coefficient or exponent\n"
    local_text += f"There are {counts['phys_count']} physical parameters (total: {counts['total_phys_count']})\n"
    local_text += (
        f"There are {counts['phys_var_count']} physical parameters with a high variability "
        f"(total: {counts['total_phys_var_count']})\n"
    )
    local_text += f"There are {counts['art_count']} artificial parameters (total: {counts['total_art_count']})\n"
    local_text += (
        f"There are {counts['art_thresh_count']} artificial threshold parameters "
        f"(total: {counts['total_art_thresh_count']})\n"
    )
    local_text += f"There are {counts['else_group_count']} not determined parameters in terms of groups\n"
    print(local_text)
    return text + local_text


def print_large_impact_no_sens(ds, top=50):
    """
    Print every model parameter that shows a large error when perturbed despite a sensitivity value of zero.
    Large is defined as an error within the top n (='top') parameters.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    top : int
        The number of top parameters to consider having a large impact.

    Returns
    -------
    String of all printed statements
    """
    text = f"\nWhich parameters had a large influence (>{top}) despite showing no sensitivity?\n"
    print(text)
    tmp_df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    tmp_df = tmp_df.loc[tmp_df["Predicted Squared Error"] == 0]
    tmp_df2 = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    tmp_df2 = tmp_df2.loc[tmp_df2["Predicted Squared Error"] != 0]
    error_key = "Mean Squared Error"
    out_params = tmp_df["Output Parameter"]

    for out_p in out_params:
        out_p = out_p.item()
        print(f"########################### {out_p} ########################")
        text += f"########################### {out_p} ########################"
        df = tmp_df.loc[tmp_df["Output Parameter"] == out_p]
        df2 = tmp_df2.loc[tmp_df2["Output Parameter"] == out_p]
        nlargest50 = df2.nlargest(top, error_key)[
            ["Input Parameter", error_key, "Predicted Squared Error"]
        ]
        min_sens = nlargest50[error_key].min()
        df_3 = df.loc[df[error_key] >= min_sens]
        text += f"sort by {error_key} with errors > {min_sens}"
        df_tmp_2 = df_3.nlargest(20, error_key)[
            ["Input Parameter", error_key, "Predicted Squared Error"]
        ]
        text += f"{df_tmp_2}"
        print(f"sort by {error_key} with errors > {min_sens}")
        print(df_tmp_2)
    return text


def print_table_top_lists(top_n_lists, top_orders_lists):
    """
    Print a (pandas) table with all parameters for each top variant.

    Parameters
    ----------
    top_n_lists : list of lists of strings
        List of lists of top n parameters generated using get_top_list().
    top_orders_lists : list of lists of strings
        List of lists of parameters within a given magnitude range generated using get_top_list().

    Returns
    -------
    String of all printed statements
    """
    print("\nTable with all parameters for each top variant\n")
    text = "\nTable with all parameters for each top variant\n"
    dict_tops = {}
    max_length = 0
    for i, t in enumerate(top_n_lists):
        t.sort()
        dict_tops[f"top_n_{i}"] = t
        if len(t) > max_length:
            max_length = len(t)
    for i, t in enumerate(top_orders_lists):
        t.sort()
        dict_tops[f"top_order_{i}"] = t
        if len(t) > max_length:
            max_length = len(t)
    for _, value in dict_tops.items():
        if len(value) < max_length:
            value.extend(["-" for _ in range(max_length - len(value))])
    df = pd.DataFrame(dict_tops)
    print(df)
    text += f"{df}"
    return text


def print_top_parameters(top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic):
    """
    Print the parameters and number of parameters with a sensitivity within a magnitude
    or within the top 10 for each output parameter.

    Parameters
    ----------
    top_magn_set : Set
        Set of parameters within one order of magnitude.
    top10_set : Set
        Set of top 10 parameters for each model state variable.
    top_magn_sens_dic : Dictionary
        Dictionary of model state variables with a list of parameters within one order of magnitude.
    top_sens_dic : Dictionary
        Dictionary of model state variables with a list of top 10 parameters for each.

    Returns
    -------
    String of all printed statements
    """
    text = f"No. of parameters within magnitude of 10**1: {len(top_magn_set)}\n"
    text += f"{top_magn_set}\n"
    text += "The parameters within a magnitude for each output Parameter:\n"
    for out_p in top_magn_sens_dic.keys():
        text += f"~*~*~*~*~*~* {out_p} ~*~*~*~*~*~*\n"
        for param in top_magn_sens_dic[out_p]:
            if param in in_params_notation_mapping:
                text += f"{param}: {in_params_notation_mapping[param][0]}\n"
            else:
                text += f"{param}\n"
        text += "\n"
    text += f"No. of parameters within the top 10: {len(top10_set)}\n"
    text += f"{top10_set}\n"
    text += "The top parameters 10 for each output Parameter:\n"
    for out_p in top_sens_dic.keys():
        text += f"~*~*~*~*~*~* {out_p} ~*~*~*~*~*~*\n"
        text += f"{top_sens_dic[out_p]}\n"
        for param in top_sens_dic[out_p]:
            if param in in_params_notation_mapping:
                text += f"{param}: {in_params_notation_mapping[param][0]}\n"
            else:
                text += f"{param}\n"
        text += "\n"
    print(text)
    return text
