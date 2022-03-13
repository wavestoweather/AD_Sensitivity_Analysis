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
\\bgroup
\\def\\arraystretch{1.3} 
\\begin{tabularx}{\\linewidth}{@{}l|X@{}|>{\hsize=.8\hsize}X|X}
\\caption{Model Parameters Investigated via AD}\\\\ \\hline
\t\\textbf{Parameter} & \\textbf{Description} & \\textbf{Value} & \\textbf{from} \\\\[6pt] \\hline
\\endhead
"""
        counter = 0
        counter2 = 0
        for param in in_params_dic:
            for deriv in in_params_dic[param]:
                counter2 += 1
                if (
                    "Not used" in get_value(deriv)
                    or "Not used" in in_params_notation_mapping[deriv][0]
                    or "Not tracked" in in_params_notation_mapping[deriv][0]
                    or "one-moment warm physics" in in_params_notation_mapping[deriv][0]
                    or "drain_cmu5" == deriv
                    or "dependent" == in_params_notation_mapping[deriv][3]
                ):
                    continue
                # print(deriv)
                # if get_value(deriv) == "$ 0e+00 $":
                #     print(deriv)
                if in_params_notation_mapping[deriv][1] == "-":
                    notation = in_params_notation_mapping[deriv][2]
                else:
                    notation = (
                        in_params_notation_mapping[deriv][1]
                        + ", "
                        + in_params_notation_mapping[deriv][2]
                    )
                table_string += (
                    "\n\t"
                    + parse_word(deriv)
                    .replace("\partial", "")
                    .replace("\mathrm", "\\text")
                    + " & "
                    + replace_cites(
                        in_params_notation_mapping[deriv][0]
                    )  # in_params_descr_dic[deriv]
                    + " & "
                    + get_value(deriv)
                    + " "
                    + get_unit(deriv)
                    + " & "
                    + replace_cites(notation)
                    + "\\\\ \\hline"
                )
                counter += 1
        table_string += "\n\\end{tabularx}\n\\egroup"
        print(f"Found {counter} parameters; tested {counter2}")
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
