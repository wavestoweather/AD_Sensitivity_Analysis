"""Generate tables for descriptions of processes.

"""
import os

from ad_sensitivity_analysis.plot.latexify import in_params_notation_mapping


def get_procs_comments(procs):
    """
    Go through "include/microphysics/user_functions.h" to extract the description of the given processes.

    Parameters
    ----------
    procs : List-like of strings.
        Name of the processes (functions) for which a description is needed.

    Returns
    -------
    Dictionary with processes as keys and their descriptions as values.
    """
    is_comment = False
    f = "include/microphysics/user_functions.h"
    if not os.path.isfile(f):
        f = "../" + f
    procs_desc = {}
    with open(f, "r", encoding="utf-8") as file_c:
        desc = ""
        for line in file_c:
            if "/**" == line[0:3]:
                is_comment = True
            if is_comment and line[0:3] == " */":
                is_comment = False
            if (
                is_comment
                and line[0:3] != "/**"
                and "@params" not in line
                and r"\@return" not in line
                and len(line) > 2
            ):
                desc += line[3::]
            if not is_comment and len(desc) > 0 and "(" in line:
                proc = line.split()[1]  # remove the return type
                proc = proc[: proc.find("(")]  # remove the opening bracket and beyond
                if proc in procs:
                    procs_desc[proc] = desc
                desc = ""
    return procs_desc


def get_process_desc_table(params, caption, label, breakup=False):
    """
    Given a list of parameters, extract all processes in which those parameters
    are used. Then create a table with those processes and a description.

    Parameters
    ----------
    params : List-like of strings
        Name of model parameters (starting with "d").
    caption : string
        Caption to put in the table.
    label : string
        Label for the table.
    breakup : bool
        If true, use ltablex package to create a table that is automatically split for multiple pages.
        Otherwise, puts the table in a "table" environment.

    Returns
    -------
        String that can be printed for a latex table.
    """
    params = list(set(params))
    procs = []
    for p in params:
        for proc in in_params_notation_mapping[p][4]:
            procs.append(proc)
    procs = list(set(procs))
    procs.sort()
    procs_desc = get_procs_comments(procs)
    empty_string = ""
    if breakup:
        table = "\\begin{{tabularx}}{{\\textwidth}}{{l|X}}\n"
        table += f"{empty_string:<4}\\toprule\n"
    else:
        table = "\\begin{{table}}[ht]\n"
        table += f"{empty_string:<4}\\centering\n"
        table += f"{empty_string:<4}\\begin{{tabularx}}{{\\textwidth}}{{l|X}}\n"

    table += f"{empty_string:<8}\\textbf{{Process}} & \\textbf{{Description}} \\\\ \n"
    if breakup:
        table += f"{empty_string:<4}\\midrule\n"
        table += f"{empty_string:<4}\\endhead\n"
    for proc in procs:
        desc = procs_desc[proc].replace("\n", f"\n{empty_string:<28}   ")
        desc = desc.replace("_", r"\_").replace("&", r"\&")
        proc_ = proc.replace("_", r"\_")
        table += f"{empty_string:<8}{proc_:<20} & {desc} \\\\ \n"
    if breakup:
        table += f"{empty_string:<4}\\bottomrule\n"
        table += f"{empty_string:<4}\\caption{{{caption}}}\n"
        table += f"{empty_string:<4}\\label{{{label}}}\n"
        table += "\\end{{tabularx}}"
    else:
        table += f"{empty_string:<4}\\end{{tabularx}}\n"
        table += f"{empty_string:<4}\\caption{{{caption}}}\n"
        table += f"{empty_string:<4}\\label{{{label}}}\n"
        table += "\\end{{table}}\n"
    return table
