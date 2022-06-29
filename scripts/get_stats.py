import copy
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_col
import numpy as np
import os
import pandas as pd
import pickle

try:
    from tqdm import tqdm
except:
    from progressbar import progressbar as tqdm
import seaborn as sns
import xarray as xr

try:
    from Deriv_dask import Deriv_dask
    from latexify import in_params_dic, physical_params
    from plot_mse import load_and_append, reduce_df
    from segment_identifier import d_unnamed
    from create_mse import load_and_append
    import latexify as latexify
except:
    from scripts.Deriv_dask import Deriv_dask
    from scripts.latexify import in_params_dic, physical_params
    from scripts.plot_mse import load_and_append, reduce_df
    from scripts.segment_identifier import d_unnamed
    from scripts.create_mse import load_and_append
    import scripts.latexify as latexify


def get_top_list(ds, print_out=True, verbose=True):
    """

    Parameters
    ----------
    ds : xarray.Dataset
    print_out : bool
        Print the number of parameters and the dataset.
    verbose : bool
        Print the top parameters for each output variable.
    Returns
    -------
    list of output parameters, list of top20 parameters, list of top10 parameters, dict of top20 parameters, dict of top10 parameters
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
        out_p = out_p
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
    if print_out:
        print(
            f"Number of distinct parameters by taking the top 20 for everything: {len(top20_list)}"
        )
        print(
            f"Number of distinct parameters by taking the top 10 for everything: {len(top10_list)}"
        )
        print(ds)
    return out_params, top20_list, top10_list, top20_sens_dic, top10_sens_dic


def get_magnitude_list(ds, out_params, print_out=True, verbose=True):
    """
    Get the top parameters within one, two, and three orders of magnitude.

    Parameters
    ----------
    ds
    out_params
    print_out
    verbose : bool
        Print the top parameters for each output variable.

    Returns
    -------
    list of parameters withing one order of magnitude, list within two orders of magnitudes, list within three orders of magnitudes
    """
    if print_out:
        print(
            "\nInstead of taking the top 10, take the parameters in the same order of magnitude\n"
        )
    top_one_order = []
    top_two_orders = []
    top_three_orders = []

    for out_p in out_params:
        out_p = out_p
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
    if print_out:
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


def print_unique_params(top_sens_dic):
    """
    Print the parameters that appear only for a single output variable.
    Parameters
    ----------
    top_sens_dic

    Returns
    -------
    String of all printed statements
    """
    text = "\nParameters that appear only for a single output variable\n"
    print(text)
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
        text += f"{pair}\n"
    return text


def print_correlation_broad(ds, out_params):
    """
    Print correlation coefficients (Spearman, Pearson, and Kendall) using each time step individually with all data,
    for each type of state variable (first moment, second moment, sedimentation), and for each output variable.
    Parameters
    ----------
    ds :
    out_params :

    Returns
    -------
    String of all printed statements
    """
    text = f"\nCorrelation with all data\n"
    print(text)

    def get_corr(df, kind):
        return df[["Predicted Squared Error", "Mean Squared Error"]].corr(kind)[
            "Predicted Squared Error"
        ][1]

    def print_corr(df):
        global text
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        text += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}\n"
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
        df = df.loc[df["Predicted Squared Error"] != 0]
        n = len(np.unique(df["Input Parameter"]))
        text += f"Correlation without zero parameters; total of {n} parameters"
        print(f"Correlation without zero parameters; total of {n} parameters")
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        text += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}"
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")

    df = ds.to_dataframe().reset_index()
    print_corr(df)

    text += "\nFor each output variable type individually\n"
    print("\nFor each output variable type individually\n")

    second_moment = ["QV", "QC", "QR", "QS", "QG", "QH", "QI"]
    first_moment = ["NCCLOUD", "NCRAIN", "NCGRAUPEL", "NCHAIL", "NCICE", "NCSNOW"]
    second_sed = ["QR_OUT", "QG_OUT", "QH_OUT", "QI_OUT", "QS_OUT"]
    first_sed = ["NR_OUT", "NG_OUT", "NH_OUT", "NI_OUT", "NS_OUT"]

    text += "##################First Moment (Number Count)"
    print("##################First Moment (Number Count)")
    df = ds.sel({"Output Parameter": first_moment}).to_dataframe().reset_index()
    print_corr(df)

    text += "##################Second Moment (Mixing Ratio)"
    print("##################Second Moment (Mixing Ratio)")
    df = ds.sel({"Output Parameter": second_moment}).to_dataframe().reset_index()
    print_corr(df)

    text += "##################First Moment Sedimentation (Number Count)"
    print("##################First Moment Sedimentation (Number Count)")
    df = ds.sel({"Output Parameter": first_sed}).to_dataframe().reset_index()
    print_corr(df)

    text += "##################Second Moment Sedimentation (Mixing Ratio)"
    print("##################Second Moment Sedimentation (Mixing Ratio)")
    df = ds.sel({"Output Parameter": second_sed}).to_dataframe().reset_index()
    print_corr(df)

    text += "\nFor each output variable individually\n"
    print("\nFor each output variable individually\n")

    for out_p in out_params:
        out_p = out_p
        text += f"##################{out_p}"
        print(f"##################{out_p}")
        df = ds.sel({"Output Parameter": [out_p]}).to_dataframe().reset_index()
        print_corr(df)
    return text


def print_correlation_mean(ds, out_params):
    """

    Parameters
    ----------
    ds
    out_params

    Returns
    -------
    String of all printed statements
    """
    text = "\nInstead of looking at correlations within time steps and trajectories individually, "
    text += "take the mean of those and look at the correlation\n"
    text += f"\nCorrelation with all data\n"
    print(text)

    def get_corr(df, kind):
        return df[["Predicted Squared Error", "Mean Squared Error"]].corr(kind)[
            "Predicted Squared Error"
        ][1]

    def print_corr(df):
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        global text
        text += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}"
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
        df = df.loc[df["Predicted Squared Error"] != 0]
        n = len(np.unique(df["Input Parameter"]))
        text += f"Correlation without zero parameters; total of {n} parameters"
        print(f"Correlation without zero parameters; total of {n} parameters")
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        text += f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}"
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")

    df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    print_corr(df)

    text += "\nFor each output variable individually\n"
    print("\nFor each output variable individually\n")
    for out_p in out_params:
        out_p = out_p
        text += f"##################{out_p}"
        print(f"##################{out_p}")
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        print_corr(df)

    text += "\nCorrelation taking different number of top parameters\n"
    print("\nCorrelation taking different number of top parameters\n")
    df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
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

    text += "\nCorrelation only for the 75th percentile\n"
    print("\nCorrelation only for the 75th percentile\n")
    params = []
    for out_p in out_params:
        out_p = out_p
        if "NH_OUT" == out_p or "QH_OUT" == out_p:
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
    params = list(set(params))
    text += f"Number of parameters: {len(params)}"
    print(f"Number of parameters: {len(params)}")
    df_tmp = (
        ds.sel({"Input Parameter": params})
        .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    r = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
        "Predicted Squared Error"
    ][1]
    text += f"Correlation for all predicted errors with ensemble errors: {r}"
    print(f"Correlation for all predicted errors with ensemble errors: {r}")
    text += f"\nCreate Latex table with Spearman correlations\n"
    print(f"\nCreate Latex table with Spearman correlations\n")
    table_corr = r"""
\begin{table}[hbt]
    \centering
    \begin{tabular}{l|c|c}
        Model State Variable $y_s$                 & $r(y_s)$ without zero sensitivities & $r(y_s)$ \\ \hline
"""
    for out_p in out_params:
        out_p = out_p
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        r = df[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
            "Predicted Squared Error"
        ][1]
        df = df.loc[df["Predicted Squared Error"] != 0]
        r_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
            "Predicted Squared Error"
        ][1]
        table_corr += (
            "\t\t"
            + latexify.parse_word(out_p).title().replace(" Of ", " of ")
            + "\t\t& $ "
            + f"{r_2:1.3f}"
            + " $ & $ "
            + f"{r:1.3f} $ "
            + r"\\"
            + "\n"
        )
    table_corr += r"""    \end{tabular}
    \caption{}
    \label{tab:validate:correlation}
\end{table}"""
    print(table_corr)
    text += table_corr

    text += f"\nCreate Latex table with Pearson correlations\n"
    print(f"\nCreate Latex table with Pearson correlations\n")
    table_corr = r"""
\begin{table}[hbt]
    \centering
    \begin{tabular}{l|c|c}
        Model State Variable $y_s$                 & $r_{\text{pearson}}(y_s)$ without zero sensitivities & $r_{\text{pearson}}(y_s)$ \\ \hline
"""
    for out_p in out_params:
        out_p = out_p
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        r = df[["Predicted Squared Error", "Mean Squared Error"]].corr("pearson")[
            "Predicted Squared Error"
        ][1]
        df = df.loc[df["Predicted Squared Error"] != 0]
        r_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr("pearson")[
            "Predicted Squared Error"
        ][1]
        table_corr += (
            "\t\t"
            + latexify.parse_word(out_p).title().replace(" Of ", " of ")
            + "\t\t& $ "
            + f"{r_2:1.3f}"
            + " $ & $ "
            + f"{r:1.3f} $ "
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


def print_latex_tables(ds, top=10, verbose=True):
    """

    Parameters
    ----------
    ds
    out_params
    top
    verbose : Bool
        If True: print the parameters while building the table.

    Returns
    -------
    sort_key_list, table_dic, the printed statements
    """
    text = "\nBuild Latex tables\n"
    print(text)
    tmp_df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )

    table_dic = {}
    sort_key_list = []

    latexify_state = {
        "QV": r"\frac{\partial Q_\vapor}{",
        "QC": r"\frac{\partial Q_\cloud}{",
        "QR": r"\frac{\partial Q_\rain}{",
        "QG": r"\frac{\partial Q_\graupel}{",
        "QH": r"\frac{\partial Q_\hail}{",
        "QI": r"\frac{\partial Q_\ice}{",
        "QS": r"\frac{\partial Q_\snow}{",
        "NCCLOUD": r"\frac{\partial N_\cloud}{",
        "NCRAIN": r"\frac{\partial N_\rain}{",
        "NCGRAUPEL": r"\frac{\partial N_\graupel}{",
        "NCHAIL": r"\frac{\partial N_\hail}{",
        "NCICE": r"\frac{\partial N_\ice}{",
        "NCSNOW": r"\frac{\partial N_\snow}{",
        "QR_OUT": r"\frac{\partial Q_{\rain, \text{out}}}{",
        "QG_OUT": r"\frac{\partial Q_{\graupel, \text{out}}}{",
        "QH_OUT": r"\frac{\partial Q_{\hail, \text{out}}}{",
        "QI_OUT": r"\frac{\partial Q_{\ice, \text{out}}}{",
        "QS_OUT": r"\frac{\partial Q_{\snow, \text{out}}}{",
        "NR_OUT": r"\frac{\partial N_{\rain, \text{out}}}{",
        "NG_OUT": r"\frac{\partial N_{\graupel, \text{out}}}{",
        "NH_OUT": r"\frac{\partial N_{\hail, \text{out}}}{",
        "NI_OUT": r"\frac{\partial N_{\ice, \text{out}}}{",
        "NS_OUT": r"\frac{\partial N_{\snow, \text{out}}}{",
    }

    top_10_table = "\\begin{table}[hbt] \n \t\\centering \n \t\\begin{tabular}{ll}"
    top_10_table += "\n \t\t\\textbf{Model State Variable} \t& \\textbf{Top 10 Parameters} \\\\ \\hline \n"
    sedi_latex = ""
    sedi_started = False
    long_table_dic = {}

    for out_p in latexify_state.keys():
        if "OUT" in out_p:
            if sedi_started:
                sedi_latex = sedi_latex[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
            else:
                sedi_latex = "\t\t Sedimentation \t& $ "
                sedi_started = True
        else:
            top_10_table += "\t\t" + latexify.parse_word(out_p) + "\t& $ "
        if verbose:
            print(f"########################### {out_p} ########################")
        df = tmp_df.loc[tmp_df["Output Parameter"] == out_p]
        # Ignore parameters that never appeared in unperturbed versions
        if np.max(df["Predicted Squared Error"]) == 0:
            continue
        if verbose:
            print("sort by sensitivity")
            print(
                df.nlargest(top, "Predicted Squared Error")[
                    ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
                ]
            )
        tmp = df.nlargest(top, "Predicted Squared Error")[
            ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
        ]
        i = 0
        for idx, row in tmp.iterrows():
            if i == 5:
                if "OUT" in out_p:
                    sedi_latex = sedi_latex[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
                else:
                    top_10_table = top_10_table[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
            i += 1
            if "OUT" in out_p:
                sedi_latex += (
                    latexify.parse_word(row["Input Parameter"])
                    .replace("$", "")
                    .replace("\partial", "")
                    + ", "
                )
            else:
                top_10_table += (
                    latexify.parse_word(row["Input Parameter"])
                    .replace("$", "")
                    .replace("\partial", "")
                    + ", "
                )
            found = False
            for val, param, state_var, l_string in sort_key_list:
                if param == row["Input Parameter"] and (
                    val < row["Predicted Squared Error"]
                    or ("N" in state_var and "N" not in out_p)
                ):
                    if "N" not in state_var and "N" in out_p:
                        break

                    found = True
                    if verbose:
                        print(f"Replace ({val}, {param}, {state_var})")
                        print(f"With (")
                        print(row["Predicted Squared Error"], end=", ")
                        print(row["Input Parameter"], end=", ")
                        print(out_p, end=")\n")

                    sort_key_list.remove((val, param, state_var, l_string))
                    break

            if row["Input Parameter"] not in table_dic or found:

                group = None
                for g in latexify.in_params_grouping:
                    if row["Input Parameter"] in latexify.in_params_grouping[g]:
                        group = g

                def latex_my_number(x):
                    if x == 0:
                        return "$ 0.00 $"
                    if x >= 100:
                        exponent = int(np.log10(x))
                        var = x / 10 ** exponent
                        return f"$ {var:2.2f} \\times 10^{ {exponent} } $"
                    elif x < 0.01:
                        exponent = math.floor(np.log10(x))
                        var = x * 10 ** (-exponent)
                        return f"$ {var:2.2f} \\times 10^{ {exponent} } $"
                    else:
                        err = row["Predicted Squared Error"]
                        return f"$ {err:2.2f} $"

                long_string = (
                    latexify.parse_word(row["Input Parameter"]).replace("\partial", "")
                    + " & "
                    + latex_my_number(row["Mean Squared Error"])
                    + " & "
                    + latex_my_number(row["Predicted Squared Error"])
                    + " & "
                    + "\\textbf{"
                    + group.title()
                    + "}: "
                    + latexify.in_params_descr_dic[row["Input Parameter"]]
                    + " \\\\ "
                )
                if out_p not in long_table_dic:
                    long_table_dic[out_p] = [long_string]
                else:
                    long_table_dic[out_p].append(long_string)
                sort_key_list.append(
                    (
                        row["Predicted Squared Error"],
                        row["Input Parameter"],
                        out_p,
                        long_string,
                    )
                )

                table_dic[row["Input Parameter"]] = (
                    "$ \displaystyle "
                    + latexify_state[out_p]
                    + latexify.parse_word(row["Input Parameter"]).replace("$", "")
                    + r"} $ & "
                    + latex_my_number(row["Mean Squared Error"])
                    #                 + f"$ {row.MSE:1.2e} $"
                    + " & "
                    + latex_my_number(row["Predicted Squared Error"])
                    #                 + f"$ {row.Sensitivity:1.2e} $"
                    + " & "
                    + "\\textbf{"
                    + group.title()
                    + "}: "
                    + latexify.in_params_descr_dic[row["Input Parameter"]]
                    + " \\\\ "
                )
        if "OUT" not in out_p:
            top_10_table = top_10_table[:-2] + " $ \\\\ \n"
        if verbose:
            print(f"sort by Predicted Squared Error")
            print(
                df.nlargest(top, "Predicted Squared Error")[
                    ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
                ]
            )

    top_10_table += sedi_latex[:-2] + " $ \\\\"
    top_10_table += "\n\t\\end{tabular} \n \t\\caption{} \n"
    top_10_table += "\t\\label{tab:} \n \\end{table} \n"
    text += "\nThe table of top 10 parameters for each state variable:\n"
    text += top_10_table
    print("\nThe table of top 10 parameters for each state variable:\n")
    print(top_10_table)

    # print("\nThe table of top 10 parameters for each state variable sorted by sensitivity:\n")
    # sort_key_list = sorted(sort_key_list, key=lambda x: x[0], reverse=True)
    # for sens, key, state_variable, desc in sort_key_list:
    #     print(table_dic[key])
    if verbose:
        print(f"There are {len(table_dic)} different input parameters")

    text += "\nAppendix table of top parameters\n"
    print("\nAppendix table of top parameters\n")
    tmp_sort = sorted(sort_key_list, key=lambda x: (x[2], x[0]), reverse=True)
    sort_dic_long_table = {}
    sort_dic_short_table = {}
    for sens, key, state_variable, l_string in tmp_sort:
        if "_OUT" in state_variable:
            state_variable2 = "Sedimentation"
        else:
            state_variable2 = state_variable
        if state_variable not in sort_dic_long_table:
            sort_dic_long_table[state_variable] = [(sens, key, l_string)]
        else:
            sort_dic_long_table[state_variable].append((sens, key, l_string))

        if state_variable2 not in sort_dic_short_table:
            sort_dic_short_table[state_variable2] = [(sens, key, l_string)]
        else:
            sort_dic_short_table[state_variable2].append((sens, key, l_string))
    table_text = "\\bgroup\n"
    table_text += "\\def\\arraystretch{1.2} %  1 is the default, we want it slightly larger such that exponents are easier to read\n"
    table_text += "\\begin{tabularx}{\\linewidth}{@{}lccX@{}}\n"
    table_text += "\t\\textbf{Model Param.}  & \\textbf{MSD} & \\textbf{Predicted MSD} & \\textbf{Parameter Description}\n"
    table_text += "\t\\endhead\n"
    i = 0
    for state_variable in latexify_state.keys():
        if state_variable not in sort_dic_long_table:
            continue
        table_text += (
            "\t\t\\hline \\multicolumn{4}{c}{"
            + latexify.parse_word(state_variable).title()
            + "}\t \\\\ \\hline\n"
        )
        for _, _, s in sort_dic_long_table[state_variable]:
            table_text += "\t\t"
            table_text += f"{s}"
            i += 1
    table_text += "\t"
    top_str = "ten"
    if top != 10:
        top_str = str(top)
    table_text += (
        r"\caption{This is the set of parameters if we gather the "
        + top_str
        + r" most important ones for each model state variable. The predicted MSD is defined in Equation~\ref{eq:identification:msd_predict}, where we only show the highest predicted MSD among all mass densities unless the parameter did not have an impact on mass densities. In that case, the predicted deviation on number density and precipitation is considered. There are $ "
        + str(i)
        + r" $ different parameters in total.}"
    )
    table_text += "\t\\label{tab:important_params}\n"
    table_text += "\\end{tabularx}\n"
    table_text += "\\egroup\n"
    text += table_text
    print(table_text)
    return sort_key_list, table_dic, text


def print_variable_with_important_params(sort_key_list):
    """

    Parameters
    ----------
    sort_key_list

    Returns
    -------
    String of all printed statements
    """
    text = (
        "\nWhich model state variable has the most parameters with a high influence?\n"
    )
    print(text)
    state_counts = {}
    for sens, key, state_variable, desc in sort_key_list:
        if state_variable not in state_counts:
            state_counts[state_variable] = 1
        else:
            state_counts[state_variable] += 1
    for state_variable in state_counts:
        print(f"{state_variable}: {state_counts[state_variable]}")
        text += f"{state_variable}: {state_counts[state_variable]}"
    return text


def print_param_types(ds, table_dic):
    """

    Parameters
    ----------
    ds
    table_dic

    Returns
    -------
    String of all printed statements
    """
    text = "\nHow many type of parameters are there?\n"
    print(text)
    geo_count = 0
    vel_count = 0
    misc_count = 0
    exp_count = 0
    coeff_count = 0
    else_count = 0
    phys_count = 0
    phys_var_count = 0
    art_count = 0
    art_thresh_count = 0
    else_group_count = 0
    for key in table_dic:
        if "geo" in key:
            geo_count += 1
        elif "vel" in key:
            vel_count += 1
        else:
            misc_count += 1

        if "exponent" in table_dic[key].lower():
            exp_count += 1
        elif "coefficient" in table_dic[key].lower():
            coeff_count += 1
        else:
            else_count += 1

        if (
            "physical" in table_dic[key].lower()
            and not "high variability" in table_dic[key].lower()
        ):
            phys_count += 1
        elif (
            "artificial" in table_dic[key].lower()
            and not "threshold" in table_dic[key].lower()
        ):
            art_count += 1
        elif "physical" in table_dic[key].lower():
            phys_var_count += 1
        elif "threshold" in table_dic[key].lower():
            art_thresh_count += 1
        else:
            else_group_count += 1

    total_phys_count = 0
    total_phys_var_count = 0
    total_art_count = 0
    total_art_thresh_count = 0
    total_parameters = len(np.unique(ds["Input Parameter"]))
    for param in np.unique(ds["Input Parameter"]):
        if param in latexify.in_params_grouping["physical"]:
            total_phys_count += 1
        if param in latexify.in_params_grouping["physical (high variability)"]:
            total_phys_var_count += 1
        if param in latexify.in_params_grouping["artificial"]:
            total_art_count += 1
        if param in latexify.in_params_grouping["artificial (threshold)"]:
            total_art_thresh_count += 1

    print(f"There are {total_parameters} many parameters")
    print(f"There are {geo_count} geometric parameters")
    print(f"There are {vel_count} velocity parameters")
    print(f"There are {misc_count} misc. parameters")
    print(f"There are {exp_count} exponents")
    print(f"There are {coeff_count} coefficients")
    print(
        f"There are {else_count} not determined parameters in terms of coefficient or exponent"
    )
    print(f"There are {phys_count} physical parameters (total: {total_phys_count})")
    print(
        f"There are {phys_var_count} physical parameters with a high variability (total: {total_phys_var_count})"
    )
    print(f"There are {art_count} artificial parameters (total: {total_art_count})")
    print(
        f"There are {art_thresh_count} artificial threshold parameters (total: {total_art_thresh_count})"
    )
    print(f"There are {else_group_count} not determined parameters in terms of groups")
    text += f"There are {total_parameters} many parameters\n"
    text += f"There are {geo_count} geometric parameters\n"
    text += f"There are {vel_count} velocity parameters\n"
    text += f"There are {misc_count} misc. parameters\n"
    text += f"There are {exp_count} exponents\n"
    text += f"There are {coeff_count} coefficients\n"
    text += f"There are {else_count} not determined parameters in terms of coefficient or exponent\n"
    text += f"There are {phys_count} physical parameters (total: {total_phys_count})\n"
    text += f"There are {phys_var_count} physical parameters with a high variability (total: {total_phys_var_count})\n"
    text += f"There are {art_count} artificial parameters (total: {total_art_count})\n"
    text += f"There are {art_thresh_count} artificial threshold parameters (total: {total_art_thresh_count})\n"
    text += (
        f"There are {else_group_count} not determined parameters in terms of groups\n"
    )
    return text


def print_large_impact_no_sens(ds, top=50):
    """

    Parameters
    ----------
    ds :
    top :

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

    for out_p in out_params:
        out_p = out_p
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

    Parameters
    ----------
    top_n_lists
    top_orders_lists

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
    for key in dict_tops:
        if len(dict_tops[key]) < max_length:
            dict_tops[key].extend(
                ["-" for _ in range(max_length - len(dict_tops[key]))]
            )
    df = pd.DataFrame(dict_tops)
    print(df)
    text += f"{df}"
    return text


def print_top_parameters(top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic):
    """

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
            try:
                text += f"{param}: {in_params_notation_mapping[param][0]}\n"
            except:
                text += f"{param}\n"
        text += "\n"
    text += f"No. of parameters within the top 10: {len(top10_set)}\n"
    text += f"{top10_set}\n"
    text += "The top parameters 10 for each output Parameter:\n"
    for out_p in top_sens_dic.keys():
        text += f"~*~*~*~*~*~* {out_p} ~*~*~*~*~*~*\n"
        text += f"{top_sens_dic[out_p]}\n"
        for param in top_sens_dic[out_p]:
            try:
                text += f"{param}: {in_params_notation_mapping[param][0]}\n"
            except:
                text += f"{param}\n"
        text += "\n"
    print(text)
    return text


def traj_get_sum_derivatives(file_path):
    """

    Parameters
    ----------
    file_path : String
        Path to NetCDF-files with sensitivities that used simulation_mode 1 (sensitivity analysis for trajectories)
        in the AD-based C++ simulation.

    Returns
    -------
    Dictionary with tracked output parameters as keys and a pandas.Dataframe with the sum of absolute values
    of the gradients.
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    out_params = ds["Output_Parameter_ID"]
    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    param_name = []
    for idx in out_params:
        param_name.append(latexify.param_id_map[idx.value])

    sums = {}
    for f in tqdm(files):
        ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
        for out_p, out_name in zip(out_params, param_name):
            ds[in_params] = np.abs(ds[in_params])
            df = (
                ds[in_params]
                .sel({"Output_Parameter_ID": out_p})
                .sum(dim=["trajectory", "time"], skipna=True)
                .to_dataframe()
                .reset_index()
            )
            df = df[in_params]

            if out_name in sums.keys():
                sums[out_name] += df
            else:
                sums[out_name] = df
    return sums, param_name


def traj_get_top_params(dict_of_df, param_name, n, orders):
    """

    Parameters
    ----------
    dict_of_df
    param_name
    n : int
        Get the top n parameters for each tracked model state variable
    orders : int or float
        Get the parameters within orders many orders of magnitude for each tracked model state variable
    Returns
    -------
    Set of parameters within the given order of magnitude, set of top n parameters,

    """
    top_sens_dic = {}
    top_magn_sens_dic = {}

    for out_name in param_name:
        top_sens_dic[out_name] = dict_of_df[out_name].T.nlargest(n, 0).T.columns.values
        tmp_df = dict_of_df[out_name].T
        max_order = np.max(tmp_df[0])
        top_magn_sens_dic[out_name] = tmp_df[
            tmp_df[0] >= max_order / (10 ** orders)
        ].T.columns.values

    tmp = []
    tmp2 = []
    for out_name in param_name:
        tmp.extend(top_magn_sens_dic[out_name])
        tmp2.extend(top_sens_dic[out_name])
    top_magn_set = set(tmp)
    top10_set = set(tmp2)
    return top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic


def get_sum_wrt(
    file_path, wrt, in_params=None, out_params=None, n_bins=100, additional_params=None
):
    """

    Parameters
    ----------
    file_path
    wrt
    in_params
    out_params
    n_bins
    additional_params

    Returns
    -------

    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
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
        param_name.append(latexify.param_id_map[idx.values])

    min_max = None
    for f in tqdm(files):
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
        min_p = np.min(ds[wrt]).values
        max_p = np.max(ds[wrt]).values
        if min_max is not None:
            if min_p < min_max[0]:
                min_max[0] = min_p
            if max_p > min_max[1]:
                min_max[1] = max_p
        else:
            min_max = [min_p, max_p]


def get_histogram(
    file_path, in_params=None, out_params=None, n_bins=100, additional_params=None
):
    """

    Parameters
    ----------
    file_path
    in_params
    out_params
    n_bins

    Returns
    -------

    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
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
        param_name.append(latexify.param_id_map[idx.values])

    min_max = {}
    min_max_in_params = {}
    for out_p in param_name:
        min_max_in_params[out_p] = {}

    for f in tqdm(files):
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")

        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(param_name)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            min_p = np.min(ds_tmp[out_name]).values
            max_p = np.max(ds_tmp[out_name]).values
            if out_name in min_max.keys():
                if min_p < min_max[out_name][0]:
                    min_max[out_name][0] = min_p
                if max_p > min_max[out_name][1]:
                    min_max[out_name][1] = max_p
            else:
                min_max[out_name] = [min_p, max_p]
            for in_p in tqdm(in_params, leave=False):
                min_p = np.min(ds_tmp[in_p]).values
                max_p = np.max(ds_tmp[in_p]).values
                if in_p in min_max_in_params[out_name]:
                    if min_p < min_max_in_params[out_name][in_p][0]:
                        min_max_in_params[out_name][in_p][0] = min_p
                    if max_p > min_max_in_params[out_name][in_p][1]:
                        min_max_in_params[out_name][in_p][1] = max_p
                else:
                    min_max_in_params[out_name][in_p] = [min_p, max_p]
        if additional_params is not None:
            for out_p in tqdm(additional_params, leave=False):
                min_p = np.min(ds[out_p]).values
                max_p = np.max(ds[out_p]).values
                if out_p in min_max.keys():
                    if min_p < min_max[out_p][0]:
                        min_max[out_p][0] = min_p
                    if max_p > min_max[out_p][1]:
                        min_max[out_p][1] = max_p
                else:
                    min_max[out_p] = [min_p, max_p]

    edges = {}
    edges_in_params = {}
    for out_p in param_name:
        delta = (min_max[out_p][1] - min_max[out_p][0]) / n_bins
        edges[out_p] = np.arange(
            min_max[out_p][0], min_max[out_p][1] + delta / 2, delta
        )
        edges_in_params[out_p] = {}
        for in_p in in_params:
            delta = (
                min_max_in_params[out_p][in_p][1] - min_max_in_params[out_p][in_p][0]
            ) / n_bins
            if min_max_in_params[out_p][in_p][0] == min_max_in_params[out_p][in_p][1]:
                continue
            else:
                edges_in_params[out_p][in_p] = np.arange(
                    min_max_in_params[out_p][in_p][0],
                    min_max_in_params[out_p][in_p][1] + delta / 2,
                    delta,
                )
    if additional_params is not None:
        for out_p in tqdm(additional_params, leave=False):
            delta = (min_max[out_p][1] - min_max[out_p][0]) / n_bins
            edges[out_p] = np.arange(
                min_max[out_p][0], min_max[out_p][1] + delta / 2, delta
            )

    hist = {}
    hist_in_params = {}
    for f in tqdm(files):
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(param_name)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            hist_tmp, _ = np.histogram(ds[out_name], edges[out_name])
            if out_name in hist:
                hist[out_name] += hist_tmp
            else:
                hist[out_name] = hist_tmp
                hist_in_params[out_name] = {}
            for in_p in tqdm(in_params, leave=False):
                if in_p not in edges_in_params[out_name]:
                    continue
                hist_tmp, _ = np.histogram(
                    ds_tmp[in_p], edges_in_params[out_name][in_p]
                )
                if in_p in hist_in_params[out_name]:
                    hist_in_params[out_name][in_p] += hist_tmp
                else:
                    hist_in_params[out_name][in_p] = hist_tmp
        if additional_params is not None:
            for out_p in tqdm(additional_params, leave=False):
                hist_tmp, _ = np.histogram(ds[out_p], edges[out_p])
                if out_p in hist:
                    hist[out_p] += hist_tmp
                else:
                    hist[out_p] = hist_tmp

    return hist, hist_in_params, edges, edges_in_params


def traj_plot_histogram_out(
    out_params,
    filename,
    edges,
    hist,
    width=24,
    height=12,
    title=None,
    verbose=False,
):
    """
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot the histogram of an output parameter,
    i.e., QV, latent_heat, latent_cool, etc.

    Parameters
    ----------
    out_params : string or list-like of strings
        Output parameter name or multiple output parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges : Dictionary of list-like of float
        Edges for the histogram. Keys must be in out_params.
    hist : Dictionary of list-like of int
        Number of entries for each bin. Keys must be in out_params.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    """
    sns.set(rc={"figure.figsize": (width, height)})

    def plot_hist(out_p, title=None):
        if title is None:
            title = f"Histogram for {latexify.parse_word(out_p)}"
        if verbose:
            print(f"Plotting histogram for {out_p}")
        ax = sns.barplot(x=edges[out_p][:-1], y=hist[out_p], color="seagreen")
        x_labels = [f"{tick:1.1e}" for tick in edges[out_p][:-1]]
        _ = ax.set_xticklabels(x_labels, rotation=45, ha="right")
        _ = ax.set_title(title)
        fig = ax.get_figure()
        # fig = ax.fig
        i = 0
        store_type = filename.split(".")[-1]
        store_path = filename[0 : -len(store_type) - 1]
        save = store_path + "_" + out_p + "_{:03d}.".format(i) + store_type
        while os.path.isfile(save):
            i = i + 1
            save = store_path + "_" + out_p + "_{:03d}.".format(i) + store_type
        plt.tight_layout()
        fig.savefig(save, dpi=300)
        plt.clf()

    if isinstance(out_params, list):
        for out_p in tqdm(out_params):
            plot_hist(out_p, title)
    else:
        plot_hist(out_params, title)
    if verbose:
        print("All plots finished!")


def traj_plot_histogram_inp(
    in_params,
    filename,
    edges_in_params,
    hist_in_params,
    width=24,
    height=12,
    title=None,
    verbose=False,
):
    """
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot three histograms per image with
    \partial output / \partial model parameter
    where output is QV, latent_heat and latent_cool. Plot one image per model_parameter.

    Parameters
    ----------
    in_params : string or list-like of strings
        Model parameter name or multiple input parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges_in_params : Dictionary of dictionary list-like of float
        Edges for the histogram. First keys are output parameters, second keys must be in in_params.
    hist_in_params : Dictionary of dictionary list-like of int
        Number of entries for each bin. First keys are output parameters, second keys must be in in_params.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    """
    sns.set(rc={"figure.figsize": (width, height)})

    def plot_hist(out_params, in_p, title=None):

        if len(out_params) != 3:
            print("The number of output params should be three.")
            print("Future versions will support varying numbers.")
        ax1 = plt.subplot(311)
        ax2 = plt.subplot(312)
        ax3 = plt.subplot(313)
        if verbose:
            print(f"Plotting histogram w.r.t. {in_p}")

        def create_fig(ax, out_p, in_p, title=None):
            if in_p not in edges_in_params[out_p]:
                return
            if title is None:
                title = f"Histogram for {out_p} w.r.t. {in_p}"
            ax_t = sns.barplot(
                x=edges_in_params[out_p][in_p][:-1],
                y=hist_in_params[out_p][in_p],
                color="seagreen",
                ax=ax,
            )
            ax_t.set_yscale("log")
            x_labels = [f"{tick:1.1e}" for tick in edges_in_params[out_p][in_p][:-1]]
            _ = ax_t.set_xticklabels(x_labels, rotation=45, ha="right")
            _ = ax_t.set_title(title)

        create_fig(ax1, out_params[0], in_p, title)
        create_fig(ax2, out_params[1], in_p, title)
        create_fig(ax3, out_params[2], in_p, title)

        plt.tight_layout()

        i = 0
        store_type = filename.split(".")[-1]
        store_path = filename[0 : -len(store_type) - 1]
        save = store_path + "_" + in_p + "_{:03d}.".format(i) + store_type
        while os.path.isfile(save):
            i = i + 1
            save = store_path + "_" + in_p + "_{:03d}.".format(i) + store_type
        plt.savefig(save, dpi=300)

        plt.clf()

    out_params = list(edges_in_params.keys())
    if isinstance(in_params, list):
        for in_p in tqdm(in_params):
            plot_hist(out_params, in_p, title)
    else:
        plot_hist(out_params, in_params, title)
    if verbose:
        print("All plots finished!")


def plot_heatmap_traj(
    in_params,
    filename,
    edges_in_params,
    hist_in_params,
    width=24,
    height=24,
    title=None,
    verbose=False,
):
    """
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot a heatmap "model parameters" over
    "bin number" where each row is another histogram for a certain model parameter.

    Parameters
    ----------
    in_params : string or list-like of strings
        Model parameter name or multiple input parameter names to plot the histogram for.
    filename : string
        Path and name of the output file. If the file already exists, a number will be appended.
    edges_in_params : Dictionary of dictionary list-like of float
        Edges for the histogram. First keys are output parameters, second keys must be in in_params.
    hist_in_params : Dictionary of dictionary list-like of int
        Number of entries for each bin. First keys are output parameters, second keys must be in in_params.
    width : float
        Width in inches
    height : float
        Height in inches
    title : string
        Title of the histogram. If none is given, a title will be generated.
    """
    sns.set(rc={"figure.figsize": (width, height)})
    # sort the histogram by a simple similarity metric. It is not perfect but better than random.
    out_params = list(edges_in_params.keys())
    for out_p in tqdm(out_params):
        in_params_tmp = []
        for in_p in in_params:
            if in_p in hist_in_params[out_p]:
                in_params_tmp.append(in_p)

        corr_matrix = np.zeros((len(in_params_tmp), len(in_params)))
        if verbose:
            print(f"Create similarity matrix for {out_p}")
        for i in range(len(in_params_tmp)):
            for j in range(len(in_params_tmp)):
                if i == j:
                    corr_matrix[i, j] = 1
                if i >= j:
                    continue
                corr_matrix[i, j] = np.corrcoef(
                    hist_in_params[out_p][in_params_tmp[i]],
                    hist_in_params[out_p][in_params_tmp[j]],
                )[0][1]
                corr_matrix[j, i] = corr_matrix[i, j]
        if verbose:
            print("Sort parameters according to correlation matrix")
        best_row = 0
        best_value = -1
        for i in range(len(in_params_tmp)):
            for j in range(len(in_params_tmp)):
                if i >= j:
                    continue
                if best_value == -1:
                    best_value = corr_matrix[i, j]
                    best_row = i
                elif best_value < corr_matrix[i, j]:
                    best_value = corr_matrix[i, j]
                    best_row = i
        in_params_sorted = [in_params_tmp[best_row]]
        previous_rows = [best_row]
        if verbose:
            print("Generate matrix for plotting")
        for i in range(len(in_params_tmp) - 1):
            next_row = 0
            next_best_value = -1
            for j in range(len(in_params_tmp)):
                if j in previous_rows:
                    continue
                if next_best_value == -1:
                    next_best_value = corr_matrix[i, j]
                    next_row = j
                elif next_best_value < corr_matrix[i, j]:
                    next_best_value = corr_matrix[i, j]
                    next_row = j
            in_params_sorted.append(in_params_tmp[next_row])
            previous_rows.append(next_row)
        # old version with similarity calculated as difference in each time step.
        # similarity_matrix = np.zeros((len(in_params_tmp), len(in_params)))
        # if verbose:
        #     print(f"Create similarity matrix for {out_p}")
        # for i in range(len(in_params_tmp)):
        #     for j in range(len(in_params_tmp)):
        #         if i >= j:
        #             continue
        #         similarity_matrix[i, j] = np.sum(
        #             np.abs( hist_in_params[out_p][in_params_tmp[i]] - hist_in_params[out_p][in_params_tmp[j]])
        #         )
        #         similarity_matrix[j, i] = similarity_matrix[i, j]
        # # get the row with the best/smallest value and let the next row be the most similar one to the previous row
        # best_row = 0
        # best_value = -1
        # if verbose:
        #     print("Sort parameters according to similarity matrix")
        # for i in range(len(in_params_tmp)):
        #     for j in range(len(in_params_tmp)):
        #         if i >= j:
        #             continue
        #         if best_value == -1:
        #             best_value = similarity_matrix[i, j]
        #             best_row = i
        #         elif best_value > similarity_matrix[i, j]:
        #             best_value = similarity_matrix[i, j]
        #             best_row = i
        # in_params_sorted = [in_params_tmp[best_row]]
        # previous_rows = [best_row]
        # if verbose:
        #     print("Generate matrix for plotting")
        # for i in range(len(in_params_tmp)-1):
        #     next_row = 0
        #     next_best_value = -1
        #     for j in range(len(in_params_tmp)):
        #         if j in previous_rows:
        #             continue
        #         if next_best_value == -1:
        #             next_best_value = similarity_matrix[i, j]
        #             next_row = j
        #         elif next_best_value > similarity_matrix[i, j]:
        #             next_best_value = similarity_matrix[i, j]
        #             next_row = j
        #     in_params_sorted.append(in_params_tmp[next_row])
        #     previous_rows.append(next_row)

        # found the correct order. Let's create a corresponding 2D-array and plot it
        n_edges = len(edges_in_params[out_p][in_params_tmp[0]][:-1])
        plot_matrix = np.zeros((len(in_params_tmp), n_edges))
        for i, in_p in enumerate(in_params_sorted):
            plot_matrix[i, :] = hist_in_params[out_p][in_p]
        plot_matrix[plot_matrix == 0] = np.nan
        ax = sns.heatmap(
            plot_matrix,
            cbar=True,
            cmap="viridis",
            norm=mpl_col.LogNorm(),
            yticklabels=in_params_sorted,
        )

        if title is None:
            title2 = f"Histogram for {out_p}"
        else:
            title2 = title
        _ = ax.set_title(title2)
        plt.tight_layout()
        fig = ax.get_figure()
        i = 0
        store_type = filename.split(".")[-1]
        store_path = filename[0 : -len(store_type) - 1]
        save = store_path + "_" + out_p + "_{:03d}.".format(i) + store_type
        while os.path.isfile(save):
            i = i + 1
            save = store_path + "_" + out_p + "_{:03d}.".format(i) + store_type
        fig.savefig(save, dpi=300)
        plt.clf()


def get_sums(file_path, store_path=None):
    """
    Calculate the sum for all gradients over all time steps and trajectories.

    Parameters
    ----------
    file_path : string
        Path to NetCDF-files with trajectories from a sensitivity analysis simulation.

    Returns
    -------
    Dictionary with model state variables as keys and the sums for each gradients as pandas.DataFrame
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    out_params = ds["Output_Parameter_ID"]
    param_name = []
    for out_p in out_params:
        param_name.append(latexify.param_id_map[out_p.values.item()])
    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    sums = {}
    for f in tqdm(files):
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
        ds[in_params] = np.abs(ds[in_params])
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            df = (
                ds[in_params]
                .sel({"Output_Parameter_ID": out_p})
                .sum(dim=["trajectory", "time"], skipna=True)
                .to_dataframe()
                .reset_index()
            )
            df = df[in_params]

            if out_name in sums.keys():
                sums[out_name] += df
            else:
                sums[out_name] = df
    if store_path is not None and store_path != "no":
        with open(store_path + "_sums.pkl", "wb") as f:
            pickle.dump(sums, f)
    return sums


def get_sums_phase(file_path, store_path=None):
    """
    Calculate the sum for all gradients over all time steps and trajectories for each phase.

    Parameters
    ----------
    file_path : string
        Path to NetCDF-files with trajectories from a sensitivity analysis simulation.

    Returns
    -------
    Dictionary with phase + model state variables as keys and the sums for each gradients as pandas.DataFrame
    """
    files = [f for f in os.listdir(file_path) if os.path.isfile(file_path + f)]
    ds = xr.open_dataset(file_path + files[0], decode_times=False, engine="netcdf4")
    out_params = ds["Output_Parameter_ID"]
    param_name = []
    for out_p in out_params:
        param_name.append(latexify.param_id_map[out_p.values.item()])
    phases = ["warm phase", "mixed phase", "ice phase", "neutral phase"]
    in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    sums = {}
    for f in tqdm(files):
        ds = xr.open_dataset(file_path + f, decode_times=False, engine="netcdf4")
        ds[in_params] = np.abs(ds[in_params])
        if ds["phase"].dtype != str and ds["phase"].dtype != np.uint64:
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            for phase_i, phase in enumerate(phases):
                if ds_tmp["phase"].dtype == str:
                    idx = np.where(ds_tmp["phase"] == phase)
                else:
                    idx = np.where(ds_tmp["phase"] == phase_i)
                dic = {in_p: None for in_p in in_params}
                for in_p in in_params:
                    dic[in_p] = [np.nansum(ds_tmp[in_p].values[idx])]
                if phase + " " + out_name in sums.keys():
                    sums[phase + " " + out_name] += pd.DataFrame(dic)
                else:
                    sums[phase + " " + out_name] = pd.DataFrame(dic)
    if store_path is not None and store_path != "no":
        with open(store_path + "_sums_phase.pkl", "wb") as f:
            pickle.dump(sums, f)
    return sums


def get_cov_matrix(input_filepath, in_params=None, filepath=None):
    """

    Parameters
    ----------
    input_filepath
    in_params
    filepath

    Returns
    -------

    """
    files = [
        f for f in os.listdir(input_filepath) if os.path.isfile(input_filepath + f)
    ]
    ds = xr.open_dataset(
        input_filepath + files[0], decode_times=False, engine="netcdf4"
    )
    if in_params is None:
        in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    out_params = ds["Output_Parameter_ID"]
    param_name = []
    for out_p in out_params:
        param_name.append(latexify.param_id_map[out_p.values.item()])
    all_params = param_name + in_params
    n = len(all_params)
    means = {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
    cov = {out_p: np.zeros((n, n), dtype=np.float64) for out_p in param_name}
    n_total = {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
    for f in tqdm(files):
        ds = xr.open_dataset(input_filepath + f, decode_times=False, engine="netcdf4")[
            all_params
        ]
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            means_tmp = ds_tmp.mean(skipna=True)
            count_tmp = (~np.isnan(ds_tmp)).sum()
            for p in means[out_name]:
                if n_total[out_name][p] > 0:
                    n_new = count_tmp[p].values.item() + n_total[out_name][p]
                    if n_new > 0:
                        means[out_name][p] = (
                            means[out_name][p] * n_total[out_name][p] / n_new
                            + means_tmp[p].values.item()
                            * count_tmp[p].values.item()
                            / n_new
                        )
                        n_total[out_name][p] = n_new
                else:
                    n_total[out_name][p] = count_tmp[p].values.item()
                    means[out_name][p] = means_tmp[p].values.item()

    n_total = {out_p: {p: 0.0 for p in all_params} for out_p in param_name}
    for f in tqdm(files):
        ds = xr.open_dataset(input_filepath + f, decode_times=False, engine="netcdf4")[
            all_params
        ]
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            count_tmp = (~np.isnan(ds_tmp)).sum()
            for i, p in enumerate(means[out_name]):
                if n_total[out_name][p] > 0:
                    n_add = count_tmp[p].values.item()
                    n_new = n_add + n_total[out_name][p]
                    if n_new > 0:
                        for i2, p2 in enumerate(means[out_name]):
                            cov_tmp = np.nanmean(
                                (ds_tmp[p].values - means[out_name][p])
                                * (ds_tmp[p2].values - means[out_name][p2])
                            )
                            cov[out_name][i, i2] = (
                                cov[out_name][i, i2] * n_total[out_name][p] / n_new
                                + cov_tmp * n_add / n_new
                            )
                        n_total[out_name][p] = n_new
                else:
                    n_new = count_tmp[p].values.item() + n_total[out_name][p]
                    if n_new > 0:
                        for i2, p2 in enumerate(means[out_name]):
                            cov[out_name][i, i2] = np.nanmean(
                                (ds_tmp[p].values - means[out_name][p])
                                * (ds_tmp[p2].values - means[out_name][p2])
                            )
                        n_total[out_name][p] = n_new
    if filepath is not None and filepath != "no":
        with open(filepath + "_means.pkl", "wb") as f:
            pickle.dump(means, f)
        with open(filepath + "_covariance_matrix.pkl", "wb") as f:
            pickle.dump(cov, f)
    return means, cov


def get_cov_matrix_phase(input_filepath, in_params=None, filepath=None):
    """

    Parameters
    ----------
    input_filepath
    in_params
    filepath

    Returns
    -------

    """
    files = [
        f for f in os.listdir(input_filepath) if os.path.isfile(input_filepath + f)
    ]
    ds = xr.open_dataset(
        input_filepath + files[0], decode_times=False, engine="netcdf4"
    )
    if in_params is None:
        in_params = [d for d in ds if (d[0] == "d" and d != "deposition")]
    out_params = ds["Output_Parameter_ID"]
    param_name = []
    for out_p in out_params:
        param_name.append(latexify.param_id_map[out_p.values.item()])
    phases = ["warm phase", "mixed phase", "ice phase", "neutral phase"]
    all_params = param_name + in_params
    n = len(all_params)
    means = {
        phase: {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
        for phase in phases
    }
    cov = {
        phase: {out_p: np.zeros((n, n), dtype=np.float64) for out_p in param_name}
        for phase in phases
    }
    n_total = {
        phase: {out_p: {in_p: 0.0 for in_p in all_params} for out_p in param_name}
        for phase in phases
    }
    for f in tqdm(files):
        ds = xr.open_dataset(input_filepath + f, decode_times=False, engine="netcdf4")[
            all_params + ["phase"]
        ]
        if ds["phase"].dtype != str and ds["phase"].dtype != np.uint64:
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})

            for phase_i, phase in enumerate(tqdm(phases, leave=False)):
                if ds_tmp["phase"].dtype == str:
                    idx = np.where(ds_tmp["phase"] == phase)
                else:
                    idx = np.where(ds_tmp["phase"] == phase_i)
                for p in tqdm(means[phase][out_name], leave=False):
                    count_tmp = np.sum((~np.isnan(ds_tmp[p].values[idx])))
                    if count_tmp == 0:
                        continue
                    mean_tmp = np.nanmean(ds_tmp[p].values[idx])

                    if n_total[phase][out_name][p] > 0:
                        n_new = count_tmp + n_total[phase][out_name][p]
                        if n_new > 0:
                            means[phase][out_name][p] = (
                                means[phase][out_name][p]
                                * n_total[phase][out_name][p]
                                / n_new
                                + mean_tmp * count_tmp / n_new
                            )
                            n_total[phase][out_name][p] = n_new
                    else:
                        n_total[phase][out_name][p] = count_tmp
                        means[phase][out_name][p] = mean_tmp

    n_total = {
        phase: {out_p: {p: 0.0 for p in all_params} for out_p in param_name}
        for phase in phases
    }
    for f in tqdm(files):
        ds = xr.open_dataset(input_filepath + f, decode_times=False, engine="netcdf4")[
            all_params + ["phase"]
        ]
        if ds["phase"].dtype != str and ds["phase"].dtype != np.uint64:
            ds["phase"] = ds["phase"].astype(np.uint64)
        for out_p, out_name in tqdm(
            zip(out_params, param_name), leave=False, total=len(out_params)
        ):
            ds_tmp = ds.sel({"Output_Parameter_ID": out_p})
            for phase_i, phase in enumerate(tqdm(phases, leave=False)):
                if ds_tmp["phase"].dtype == str:
                    idx = np.where(ds_tmp["phase"] == phase)
                else:
                    idx = np.where(ds_tmp["phase"] == phase_i)
                for i, p in enumerate(tqdm(means[phase][out_name], leave=False)):
                    count_tmp = np.sum((~np.isnan(ds_tmp[p].values[idx])))
                    if count_tmp == 0:
                        continue
                    p_variance = ds_tmp[p].values[idx] - means[phase][out_name][p]
                    if n_total[phase][out_name][p] > 0:
                        n_add = count_tmp
                        n_new = n_add + n_total[phase][out_name][p]
                        if n_new > 0 and n_add > 0:
                            for i2, p2 in enumerate(means[phase][out_name]):
                                cov_tmp = np.nanmean(
                                    p_variance
                                    * (
                                        ds_tmp[p2].values[idx]
                                        - means[phase][out_name][p2]
                                    )
                                )
                                cov[phase][out_name][i, i2] = (
                                    cov[phase][out_name][i, i2]
                                    * n_total[phase][out_name][p]
                                    / n_new
                                    + cov_tmp * n_add / n_new
                                )
                            n_total[phase][out_name][p] = n_new
                    else:
                        n_new = count_tmp + n_total[phase][out_name][p]
                        if n_new > 0:
                            for i2, p2 in enumerate(means[phase][out_name]):
                                cov[phase][out_name][i, i2] = np.nanmean(
                                    p_variance
                                    * (
                                        ds_tmp[p2].values[idx]
                                        - means[phase][out_name][p2]
                                    )
                                )
                            n_total[phase][out_name][p] = n_new
    if filepath is not None and filepath != "no":
        with open(filepath + "_means_phases.pkl", "wb") as f:
            pickle.dump(means, f)
        with open(filepath + "_covariance_matrix_phases.pkl", "wb") as f:
            pickle.dump(cov, f)
    return means, cov


def plot_heatmap(
    data_in,
    out_param,
    names,
    in_params=None,
    plot_type="all",
    norm=None,
    title=None,
    filename=None,
    width=30,
    height=16,
):
    """

    Parameters
    ----------
    data_in
    out_param
    names
    in_params
    plot_type
    norm
    title
    filename
    width
    height

    Returns
    -------

    """
    sns.set(rc={"figure.figsize": (width, height)})
    data = copy.deepcopy(data_in[out_param])
    if in_params is not None and len(in_params) < np.shape(data_in[out_param])[0]:
        if len(in_params) == 0:
            return
        idx = []
        names = np.asarray(list(names))
        for in_p in in_params:
            idx.append(np.where(names == in_p)[0][0])
        idx = np.asarray(idx)
        data = data[idx]
        data = data[:, idx]
        names = names[idx]
    if plot_type == "negative":
        data[np.where(data >= 0)] = np.nan
    elif plot_type == "positive":
        data[np.where(data < 0)] = np.nan
    if norm is not None and plot_type != "positive":
        data = np.abs(data)
    try:
        g = sns.heatmap(
            data=data,
            cmap="viridis",
            norm=norm,
            cbar=True,
            yticklabels=names,
            xticklabels=names,
            annot=(len(names) <= 20),
            fmt="1.2e",
        )
        if title is None:
            title_ = f"Heatmap of Covariance (Gradients for {out_param})"
        else:
            title_ = title
        _ = g.set_title(title_)
        plt.tight_layout()
        fig = g.get_figure()
        if filename is not None:
            i = 0
            store_type = filename.split(".")[-1]
            store_path = filename[0 : -len(store_type) - 1]
            save = store_path + "_{:03d}.".format(i) + store_type
            while os.path.isfile(save):
                i = i + 1
                save = store_path + "_{:03d}.".format(i) + store_type
            fig.savefig(save, dpi=300)
        return g
    except:
        if filename is not None:
            print(f"Cannot create plot for {filename} ({out_param})")
        else:
            print(f"Cannot create plot for {out_param}")
        return None


if __name__ == "__main__":
    import argparse
    from latexify import *

    parser = argparse.ArgumentParser(
        description="""
        Get statistics of a final, post-processed dataset with mean squared deviation and 
        predicted mean squared deviation.
        Or get statistics and plot histograms for files from a sensitivity analysis simulation along
        trajectories, e.g., by using
        python get_stats.py --file /project/meteo/w2w/Z2/Z2_data_gradients/ --out_file /path/to/pics/histogram.png 
        The name of the plots will be changed automatically to store multiple plots.
        Beware that creating the histogram may take a while. You can use 
        --save_histogram /path/to/folder/
        to store the histogram and edges to disk. 
        Some statistics are done after plotting which may take a while as well.
        """
    )
    parser.add_argument(
        "--file",
        default="../data/vladiana_ensembles_postprocess/merged_independent.nc",
        help="""
        Path to post-processed file or to a folder with many files from a sensitivity analysis simulation.
        """,
    )
    parser.add_argument(
        "--out_file",
        default="../pics/histogram.png",
        help="""
        Path and name to store histogram plots if the input is a set of trajectories with a sensitivity analysis
        simulation. Exchanges the ending with 'txt' and stores the results of any statistics in there unless it is 'none'.
        """,
    )
    parser.add_argument(
        "--width",
        default=24,
        type=float,
        help="""
        Width in inches for histogram plots.
        """,
    )
    parser.add_argument(
        "--height",
        default=12,
        type=float,
        help="""
        Height in inches for histogram plots.
        """,
    )
    parser.add_argument(
        "--plot_type",
        default="all",
        help="""
        Choose which plots to create. Options are
        all: All plots.
        hist_out: Histogram for output parameters.
        hist_in: Histogram for all input parameters.
        heat: Heatmap for all parameters.
        cov_heat: Heatmap of covariance matrix.
        cov_heat_phases: Heatmap of covariance matrix for each phase.
        none: No plots.
        """,
    )
    parser.add_argument(
        "--additional_hist_params",
        type=str,
        nargs="+",
        default=[],
        help="""
        Additional parameters to create a histogram for, such as pressure or asc600. 
        """,
    )
    parser.add_argument(
        "--load_histogram",
        default="no",
        help="""
        Load the histogram and edges with pickle from this path.
        """,
    )
    parser.add_argument(
        "--save_histogram",
        default="no",
        help="""
        Store the histogram and edges with pickle to this path.
        """,
    )
    parser.add_argument(
        "--save_statistics",
        default="no",
        help="""
        Store the covariance matrix and means to this path. This includes the matrix with and without
        phases. Also store the sums for all gradients to this path.
        """,
    )
    parser.add_argument(
        "--load_statistics",
        default="no",
        help="""
        Load the sums for all gradients from this path.
        """,
    )
    parser.add_argument(
        "--load_covariance",
        default="no",
        help="""
        Load the covariance matrix and means from this path. This includes the matrix with and without
        phases. 
        """,
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="""
        More output, i.e., in each intermediate step for building tables.
        """,
    )

    args = parser.parse_args()
    if args.file.endswith("/"):
        if (
            args.plot_type == "all"
            or args.plot_type == "hist_out"
            or args.plot_type == "hist_in"
            or args.plot_type == "heat"
            or args.save_histogram != "no"
        ):
            if args.load_histogram != "no":
                print("########### Load histograms ###########")
                file_path = args.load_histogram
                if not file_path.endswith("/"):
                    file_path += "/"
                with open(file_path + "hist.pkl", "rb") as f:
                    hist = pickle.load(f)
                with open(file_path + "edges.pkl", "rb") as f:
                    edges = pickle.load(f)
                with open(file_path + "hist_in_params.pkl", "rb") as f:
                    hist_in_params = pickle.load(f)
                with open(file_path + "edges_in_params.pkl", "rb") as f:
                    edges_in_params = pickle.load(f)
            else:
                print("########### Calculate histograms ###########")
                hist, hist_in_params, edges, edges_in_params = get_histogram(
                    args.file,
                    additional_params=args.additional_hist_params,
                )
            if args.save_histogram != "no":
                print("########### Store histograms ###########")
                file_path = args.save_histogram
                with open(file_path + "hist.pkl", "wb") as f:
                    pickle.dump(hist, f)
                with open(file_path + "edges.pkl", "wb") as f:
                    pickle.dump(edges, f)
                with open(file_path + "hist_in_params.pkl", "wb") as f:
                    pickle.dump(hist_in_params, f)
                with open(file_path + "edges_in_params.pkl", "wb") as f:
                    pickle.dump(edges_in_params, f)

            if args.plot_type == "all" or args.plot_type == "hist_out":
                print(
                    "########### Plot histograms for model state variables ###########"
                )
                traj_plot_histogram_out(
                    out_params=list(edges.keys()),
                    filename=args.out_file,
                    edges=edges,
                    hist=hist,
                    width=args.width,
                    height=args.height,
                    title=None,
                    verbose=args.verbose,
                )
            if args.plot_type == "all" or args.plot_type == "hist_in":
                print("########### Plot histograms for model parameters ###########")
                traj_plot_histogram_inp(
                    in_params=list(edges_in_params[list(edges.keys())[0]].keys()),
                    filename=args.out_file,
                    edges_in_params=edges_in_params,
                    hist_in_params=hist_in_params,
                    width=args.width,
                    height=args.height,
                    title=None,
                    verbose=args.verbose,
                )
            if args.plot_type == "all" or args.plot_type == "heat":
                print("########### Plot heatmaps ###########")
                plot_heatmap_traj(
                    in_params=list(edges_in_params[list(edges.keys())[0]].keys()),
                    filename=args.out_file,
                    edges_in_params=edges_in_params,
                    hist_in_params=hist_in_params,
                    width=args.width,
                    height=args.height,
                    title=None,
                    verbose=args.verbose,
                )

        print("########### Some statistics ###########")
        if args.load_statistics != "no":
            print("########### Load sums ###########")
            file_path = args.load_statistics
            with open(file_path + "_sums.pkl", "rb") as f:
                sums = pickle.load(f)
            with open(file_path + "_sums_phase.pkl", "rb") as f:
                sums_phase = pickle.load(f)
        else:
            sums = get_sums(args.file, args.save_statistics)
            sums_phase = get_sums_phase(args.file, args.save_statistics)

        if args.load_covariance != "no":
            print("########### Load covariance matrix ###########")
            file_path = args.load_covariance
            with open(file_path + "_means.pkl", "rb") as f:
                means = pickle.load(f)
            with open(file_path + "_covariance_matrix.pkl", "rb") as f:
                cov = pickle.load(f)
            with open(file_path + "_means_phases.pkl", "rb") as f:
                means_phase = pickle.load(f)
            with open(file_path + "_covariance_matrix_phases.pkl", "rb") as f:
                cov_phase = pickle.load(f)

        top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic = traj_get_top_params(
            sums, sums.keys(), 10, 1
        )
        (
            top_magn_set_phase,
            top10_set_phase,
            top_magn_sens_dic_phase,
            top_sens_dic_phase,
        ) = traj_get_top_params(sums_phase, sums_phase.keys(), 10, 1)
        text = print_top_parameters(
            top_magn_set, top10_set, top_magn_sens_dic, top_sens_dic
        )
        text += "\n########### With phases ###########\n"
        text += print_top_parameters(
            top_magn_set_phase,
            top10_set_phase,
            top_magn_sens_dic_phase,
            top_sens_dic_phase,
        )
        if args.out_file != "none":
            filename = args.out_file
            ending = filename.split(".")[-1]
            filename = filename[0 : -len(ending) - 1] + ".txt"
            with open(filename, "w+") as f:
                f.write(text)

        if (
            args.plot_type == "all" or args.plot_type == "cov_heat"
        ) and args.load_covariance == "no":
            print("########### Calculate covariance matrix ###########")
            means, cov = get_cov_matrix(
                args.file,
                in_params=list(top_magn_set),
                filepath=args.save_statistics,
            )
        if (
            args.plot_type == "all" or args.plot_type == "cov_heat_phases"
        ) and args.load_covariance == "no":
            print("########### Calculate covariance matrix with phases ###########")
            means_phase, cov_phase = get_cov_matrix_phase(
                args.file,
                in_params=list(top_magn_set_phase),
                filepath=args.save_statistics,
            )

        store_type = args.out_file.split(".")[-1]
        store_path = args.out_file[0 : -len(store_type) - 1]
        if args.plot_type == "all" or args.plot_type == "cov_heat":
            print("########### Plot covariance matrix ###########")
            for out_p in tqdm(cov):
                _ = plot_heatmap(
                    data_in=cov,
                    out_param=out_p,
                    names=list(means[out_p].keys()),
                    title=f"Heatmap of Covariance (Gradients for {out_p})",
                    plot_type="all",
                    norm=mpl_col.LogNorm(),
                    filename=f"{store_path}_{out_p}_all_log.{store_type}",
                    width=args.width,
                    height=args.height,
                )
                plt.clf()
                _ = plot_heatmap(
                    data_in=cov,
                    out_param=out_p,
                    names=list(means[out_p].keys()),
                    title=f"Heatmap of Covariance (Gradients for {out_p})",
                    plot_type="all",
                    filename=f"{store_path}_{out_p}_all.{store_type}",
                    width=args.width,
                    height=args.height,
                )
                plt.clf()
                _ = plot_heatmap(
                    data_in=cov,
                    out_param=out_p,
                    names=list(means[out_p].keys()),
                    in_params=top_magn_sens_dic[out_p],
                    title=f"Heatmap of Covariance (Gradients for {out_p})",
                    plot_type="all",
                    norm=mpl_col.LogNorm(),
                    filename=f"{store_path}_{out_p}_some_log.{store_type}",
                    width=args.width,
                    height=args.height,
                )
                plt.clf()
                _ = plot_heatmap(
                    data_in=cov,
                    out_param=out_p,
                    names=list(means[out_p].keys()),
                    in_params=top_magn_sens_dic[out_p],
                    title=f"Heatmap of Covariance (Gradients for {out_p})",
                    plot_type="all",
                    filename=f"{store_path}_{out_p}_some.{store_type}",
                    width=args.width,
                    height=args.height,
                )
                plt.clf()
        if args.plot_type == "all" or args.plot_type == "cov_heat_phases":
            print("########### Plot covariance matrix with phases ###########")
            for phase in tqdm(cov_phase):
                if "neutral" in phase:
                    continue
                phase_strip = phase.strip()
                for out_p in tqdm(cov_phase[phase], leave=False):
                    _ = plot_heatmap(
                        data_in=cov_phase[phase],
                        out_param=out_p,
                        names=list(means_phase[phase][out_p].keys()),
                        title=f"({phase_strip}) Heatmap of Covariance (Gradients for {out_p})",
                        plot_type="all",
                        norm=mpl_col.LogNorm(),
                        filename=f"{store_path}_{phase_strip}_{out_p}_all_log.{store_type}",
                        width=args.width,
                        height=args.height,
                    )
                    plt.clf()
                    _ = plot_heatmap(
                        data_in=cov_phase[phase],
                        out_param=out_p,
                        names=list(means_phase[phase][out_p].keys()),
                        title=f"({phase_strip}) Heatmap of Covariance (Gradients for {out_p})",
                        plot_type="all",
                        filename=f"{store_path}_{phase_strip}_{out_p}_all.{store_type}",
                        width=args.width,
                        height=args.height,
                    )
                    plt.clf()
                    _ = plot_heatmap(
                        data_in=cov_phase[phase],
                        out_param=out_p,
                        names=means_phase[phase][out_p].keys(),
                        in_params=top_magn_sens_dic_phase[phase + " " + out_p],
                        title=f"({phase_strip}) Heatmap of Covariance (Gradients for {out_p})",
                        plot_type="all",
                        norm=mpl_col.LogNorm(),
                        filename=f"{store_path}_{phase_strip}_{out_p}_some_log.{store_type}",
                        width=args.width,
                        height=args.height,
                    )
                    plt.clf()
                    _ = plot_heatmap(
                        data_in=cov_phase[phase],
                        out_param=out_p,
                        names=means_phase[phase][out_p].keys(),
                        in_params=top_magn_sens_dic_phase[phase + " " + out_p],
                        title=f"({phase_strip}) Heatmap of Covariance (Gradients for {out_p})",
                        plot_type="all",
                        filename=f"{store_path}_{phase_strip}_{out_p}_some.{store_type}",
                        width=args.width,
                        height=args.height,
                    )
                    plt.clf()
    else:
        ds = xr.open_dataset(args.file, decode_times=False)
        (
            out_params,
            top20_list,
            top10_list,
            top20_sens_dic,
            top10_sens_dic,
        ) = get_top_list(ds, True, args.verbose)
        (
            top_one_order_list,
            top_two_orders_list,
            top_three_orders_list,
        ) = get_magnitude_list(ds, out_params, True, args.verbose)

        pd.set_option("display.max_rows", 100)
        pd.set_option("display.max_columns", 10)

        with pd.option_context(
            "display.max_rows",
            100,
            "display.max_columns",
            10,
            "display.expand_frame_repr",
            False,
        ):
            text = print_table_top_lists(
                [top10_list, top20_list],
                [top_one_order_list, top_two_orders_list, top_three_orders_list],
            )

        text += print_unique_params(top10_sens_dic)
        text += print_correlation_broad(ds, out_params)
        text += print_correlation_mean(ds, out_params)
        sort_key_list, table_dic, text_tmp = print_latex_tables(ds, 10, args.verbose)
        text += text_tmp
        text += print_variable_with_important_params(sort_key_list)
        text += print_param_types(ds, table_dic)
        text += print_large_impact_no_sens(ds)
        if args.out_file != "none":
            filename = args.out_file
            ending = filename.split(".")[-1]
            filename = filename[0 : -len(ending) - 1] + ".txt"
            with open(filename, "w+") as f:
                f.write(text)
