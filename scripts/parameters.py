from latexify import *


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Create a string that is a table of Parameter, Description, Value, Unit
        and print it to standard output.
        """
    )
    parser.add_argument(
        "--output_style",
        default="simple",
        help="""
        Output style. Options are
        latex: Latex tabular.
        csv: comma separated
        simple: intended for output on a console
        """,
    )
    args = parser.parse_args()

    output_style = args.output_style

    if output_style == "simple":
        print("Model Parameters Investigated via AD")
        t1 = "Parameter"
        t2 = "Description"
        t3 = "Value"
        t4 = "Unit"
        print(f"{t1:28}\t{t2:107}\t{t3} \t{t4}")
        for param in in_params_dic:
            for deriv in in_params_dic[param]:
                if "Not used" in in_params_descr_dic[deriv]:
                    continue

                print(
                    f"{deriv:28}\t{in_params_descr_dic[deriv]:107}\t{get_value(deriv)}\t(deriv)"
                )
    elif output_style == "latex":
        table_string = """
\\begin{tabularx}{\\linewidth}{@{}cXcc@{}}
\\caption{Model Parameters Investigated via AD}\\\\
\t\\textbf{Parameter} & \\textbf{Description} & \\textbf{Value} & \\textbf{Unit} \\\\[6pt]
\\endhead
"""
        for param in in_params_dic:
            for deriv in in_params_dic[param]:
                if (
                    "Not used" in get_value(deriv)
                    or "Not used" in in_params_descr_dic[deriv]
                ):
                    continue
                table_string += (
                    "\n\t"
                    + parse_word(deriv).replace("\partial", "")
                    + " & "
                    + in_params_descr_dic[deriv]
                    + " & "
                    + get_value(deriv, latex=True)
                    + " & "
                    + get_unit(deriv)
                    + "\\\\"
                )
        table_string += "\n\\end{tabularx}"
        print(table_string)
    elif output_style == "csv":
        print("Model Parameters Investigated via AD")
        print("Parameter,Description,Value,Unit")
        for param in in_params_dic:
            for deriv in in_params_dic[param]:
                if (
                    "Not used" in get_value(deriv)
                    or "Not used" in in_params_descr_dic[deriv]
                ):
                    continue
                print(
                    f"{deriv},{in_params_descr_dic[deriv]},{get_value(deriv)},{get_unit(deriv)}"
                )
    else:
        print(f"Error: {output_style} - No such output style")
