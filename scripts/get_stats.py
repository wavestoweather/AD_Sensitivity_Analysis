import math
import numpy as np
import pandas as pd
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
        print(top_two_orders_list)
    return top_one_order_list, top_two_orders_list, top_three_orders_list


def print_unique_params(top_sens_dic):
    """
    Print the parameters that appear only for a single output variable.
    Parameters
    ----------
    top_sens_dic

    """
    print("\nParameters that appear only for a single output variable\n")
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


def print_correlation_broad(ds, out_params):
    """
    Print correlation coefficients (Spearman, Pearson, and Kendall) using each time step individually with all data,
    for each type of state variable (first moment, second moment, sedimentation), and for each output variable.
    Parameters
    ----------
    ds :
    out_params :

    """
    print(f"\nCorrelation with all data\n")

    def get_corr(df, kind):
        return (df[["Predicted Squared Error", "Mean Squared Error"]].corr(kind)[
            "Predicted Squared Error"][1]
        )
    def print_corr(df):
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
        df = df.loc[df["Predicted Squared Error"] != 0]
        n = len(np.unique(df["Input Parameter"]))
        print(f"Correlation without zero parameters; total of {n} parameters")
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
    df = ds.to_dataframe().reset_index()
    print_corr(df)

    print("\nFor each output variable type individually\n")

    second_moment = ["QV", "QC", "QR", "QS", "QG", "QH", "QI"]
    first_moment = ["NCCLOUD", "NCRAIN", "NCGRAUPEL", "NCHAIL", "NCICE", "NCSNOW"]
    second_sed = ["QR_OUT", "QG_OUT", "QH_OUT", "QI_OUT", "QS_OUT"]
    first_sed = ["NR_OUT", "NG_OUT", "NH_OUT", "NI_OUT", "NS_OUT"]

    print("##################First Moment (Number Count)")
    df = ds.sel({"Output Parameter": first_moment}).to_dataframe().reset_index()
    print_corr(df)

    print("##################Second Moment (Mixing Ratio)")
    df = ds.sel({"Output Parameter": second_moment}).to_dataframe().reset_index()
    print_corr(df)

    print("##################First Moment Sedimentation (Number Count)")
    df = ds.sel({"Output Parameter": first_sed}).to_dataframe().reset_index()
    print_corr(df)

    print("##################Second Moment Sedimentation (Mixing Ratio)")
    df = ds.sel({"Output Parameter": second_sed}).to_dataframe().reset_index()
    print_corr(df)

    print("\nFor each output variable individually\n")

    for out_p in out_params:
        out_p = out_p
        print(f"##################{out_p}")
        df = ds.sel({"Output Parameter": [out_p]}).to_dataframe().reset_index()
        print_corr(df)


def print_correlation_mean(ds, out_params):
    """

    Parameters
    ----------
    ds
    out_params

    """
    print(
        "\nInstead of looking at correlations within time steps and trajectories individually, "
    )
    print("take the mean of those and look at the correlation\n")
    print(f"\nCorrelation with all data\n")

    def get_corr(df, kind):
        return (df[["Predicted Squared Error", "Mean Squared Error"]].corr(kind)[
            "Predicted Squared Error"][1]
        )
    def print_corr(df):
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
        df = df.loc[df["Predicted Squared Error"] != 0]
        n = len(np.unique(df["Input Parameter"]))
        print(f"Correlation without zero parameters; total of {n} parameters")
        spearman = get_corr(df, "spearman")
        pearson = get_corr(df, "pearson")
        kendall = get_corr(df, "kendall")
        print(f"Spearman: {spearman}, Kendall: {kendall}, Pearson: {pearson}")
    df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
    )
    print_corr(df)

    print("\nFor each output variable individually\n")
    for out_p in out_params:
        out_p = out_p
        print(f"##################{out_p}")
        df = (
            ds.sel({"Output Parameter": [out_p]})
            .mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        print_corr(df)

    print("\nCorrelation taking different number of top parameters\n")
    df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )
    tuples = []
    for n in [10, 20, 50, 100, 150, 200]:
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
        print(f"Correlation with zero parameters; total of {i} parameters")
        print(
            df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
                "Predicted Squared Error"
            ][1]
        )
        pearson = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "pearson"
        )["Predicted Squared Error"][1]
        kendall = df_tmp[["Predicted Squared Error", "Mean Squared Error"]].corr(
            "kendall"
        )["Predicted Squared Error"][1]
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
        print(f"n={t[0]}, r={t[1]:1.3f}")

    print("\nCorrelation only for the 75th percentile\n")
    params = []
    for out_p in out_params:
        out_p = out_p
        if "NH_OUT" == out_p or "QH_OUT" == out_p:
            continue
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
        print(f"Correlation only with 75th percentile; total of {n} parameters")
        r_2 = df[["Predicted Squared Error", "Mean Squared Error"]].corr("spearman")[
            "Predicted Squared Error"
        ][1]
        print(r_2)
        pearson = df[["Predicted Squared Error", "Mean Squared Error"]].corr("pearson")[
            "Predicted Squared Error"
        ][1]
        kendall = df[["Predicted Squared Error", "Mean Squared Error"]].corr("kendall")[
            "Predicted Squared Error"
        ][1]
        print(f"Pearson: {pearson}, Kendall: {kendall}")
    params = list(set(params))
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
    print(f"Correlation for all predicted errors with ensemble errors: {r}")

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
    sort_key_list, table_dic
    """
    print("\nBuild Latex tables\n")
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
    print("\nThe table of top 10 parameters for each state variable:\n")
    print(top_10_table)

    # print("\nThe table of top 10 parameters for each state variable sorted by sensitivity:\n")
    # sort_key_list = sorted(sort_key_list, key=lambda x: x[0], reverse=True)
    # for sens, key, state_variable, desc in sort_key_list:
    #     print(table_dic[key])
    if verbose:
        print(f"There are {len(table_dic)} different input parameters")

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
    print("\\bgroup")
    print(
        "\\def\\arraystretch{1.2} %  1 is the default, we want it slightly larger such that exponents are easier to read"
    )
    print("\\begin{tabularx}{\\linewidth}{@{}lccX@{}}")
    print(
        "\t\\textbf{Model Param.}  & \\textbf{MSD} & \\textbf{Predicted MSD} & \\textbf{Parameter Description}"
    )
    print("\t\\endhead")
    i = 0
    for state_variable in latexify_state.keys():
        if state_variable not in sort_dic_long_table:
            continue
        print(
            "\t\t\\hline \\multicolumn{4}{c}{"
            + latexify.parse_word(state_variable).title()
            + "}\t \\\\ \\hline"
        )
        for _, _, s in sort_dic_long_table[state_variable]:
            print("\t\t", end="")
            print(s)
            i += 1

    print("\t", end="")
    top_str = "ten"
    if top != 10:
        top_str = str(top)
    print(
        r"\caption{This is the set of parameters if we gather the "
        + top_str
        + r" most important ones for each model state variable. The predicted MSD is defined in Equation~\ref{eq:identification:msd_predict}, where we only show the highest predicted MSD among all mass densities unless the parameter did not have an impact on mass densities. In that case, the predicted deviation on number density and precipitation is considered. There are $ "
        + str(i)
        + r" $ different parameters in total.}"
    )
    print("\t\\label{tab:important_params}")
    print("\\end{tabularx}")
    print("\\egroup")
    return sort_key_list, table_dic


def print_variable_with_important_params(sort_key_list):
    """

    Parameters
    ----------
    sort_key_list


    """
    print(
        "\nWhich model state variable has the most parameters with a high influence?\n"
    )
    state_counts = {}
    for sens, key, state_variable, desc in sort_key_list:
        if state_variable not in state_counts:
            state_counts[state_variable] = 1
        else:
            state_counts[state_variable] += 1
    for state_variable in state_counts:
        print(f"{state_variable}: {state_counts[state_variable]}")


def print_param_types(ds, table_dic):
    """

    Parameters
    ----------
    ds
    table_dic


    """
    print("\nHow many type of parameters are there?\n")
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


def print_large_impact_no_sens(ds, top=50):
    """

    Parameters
    ----------
    ds :
    top :

    """
    print(
        f"\nWhich parameters had a large influence (>{top}) despite showing no sensitivity?\n"
    )
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
        df = tmp_df.loc[tmp_df["Output Parameter"] == out_p]
        df2 = tmp_df2.loc[tmp_df2["Output Parameter"] == out_p]
        nlargest50 = df2.nlargest(top, error_key)[
            ["Input Parameter", error_key, "Predicted Squared Error"]
        ]
        min_sens = nlargest50[error_key].min()
        df_3 = df.loc[df[error_key] >= min_sens]
        print(f"sort by {error_key} with errors > {min_sens}")
        print(
            df_3.nlargest(20, error_key)[
                ["Input Parameter", error_key, "Predicted Squared Error"]
            ]
        )


def print_table_top_lists(top_n_lists, top_orders_lists):
    """

    Parameters
    ----------
    top_n_lists
    top_orders_lists

    """
    print("\nTable with all parameters for each top variant\n")
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


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Get statistics of the final, post-processed dataset.
        """
    )
    parser.add_argument(
        "--file",
        default="../data/vladiana_ensembles_postprocess/merged_independent.nc",
        help="""
        Path to post-processed file.
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
    ds = xr.open_dataset(args.file, decode_times=False)
    out_params, top20_list, top10_list, top20_sens_dic, top10_sens_dic = get_top_list(
        ds, True, args.verbose
    )
    top_one_order_list, top_two_orders_list, top_three_orders_list = get_magnitude_list(
        ds, out_params, True, args.verbose
    )

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
        print_table_top_lists(
            [top10_list, top20_list],
            [top_one_order_list, top_two_orders_list, top_three_orders_list],
        )

    print_unique_params(top10_sens_dic)
    print_correlation_broad(ds, out_params)
    print_correlation_mean(ds, out_params)
    sort_key_list, table_dic = print_latex_tables(ds, 10, args.verbose)
    print_variable_with_important_params(sort_key_list)
    print_param_types(ds, table_dic)
    print_large_impact_no_sens(ds)
