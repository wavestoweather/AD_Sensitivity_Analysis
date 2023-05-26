"""Uses plot.mse.py to plot some examples.

"""
import itertools
import numpy as np
import pandas as pd
import xarray as xr

from ad_sensitivity_analysis.plot.latexify import parse_word, get_unit, param_id_map
from ad_sensitivity_analysis.plot.mse import plot_mse
from ad_sensitivity_analysis.plot.mse_histogram import plot_histogram
from ad_sensitivity_analysis.plot.mse_time_evolution import plot_time_evolution


def _parse_sensitivity_ds(ds_traj, args):
    """

    Parameters
    ----------
    ds_traj
    args

    Returns
    -------

    """
    if args.min_time is not None and args.max_time is not None:
        ds_mean = np.abs(
            ds_traj.where(
                (ds_traj["time_after_ascent"] >= args.min_time)
                & (ds_traj["time_after_ascent"] <= args.max_time)
            )
        ).mean(dim=["time"], skipna=True)
    elif args.min_time is not None:
        ds_mean = np.abs(
            ds_traj.where(ds_traj["time_after_ascent"] >= args.min_time)
        ).mean(dim=["time"], skipna=True)
    elif args.max_time is not None:
        ds_mean = np.abs(
            ds_traj.where(ds_traj["time_after_ascent"] <= args.max_time)
        ).mean(dim=["time"], skipna=True)
    else:
        ds_mean = np.abs(ds_traj).mean(dim=["time"], skipna=True)
    return ds_mean


def parse_sensitivity_plot(ds, out_params, param_title_names, args):
    """

    Parameters
    ----------
    ds
    out_params
    param_title_names
    args

    Returns
    -------

    """
    for out_p in out_params:
        out_p_id = np.argwhere(np.asarray(param_id_map) == out_p).item()
        ds_traj = ds.isel({"ensemble": 0, "trajectory": args.traj}).sel(
            {"Output_Parameter_ID": out_p_id}
        )
        ds_mean = _parse_sensitivity_ds(ds_traj, args)
        in_params = []
        max_vals = np.zeros(args.n_model_params)
        for i, key in itertools.product(range(args.n_model_params), ds_mean):
            if key[0] == "d" and key != "deposition":
                if ds_mean[key].values.item() > max_vals[i] and key not in in_params:
                    if len(in_params) <= i:
                        in_params.append(key)
                    else:
                        in_params[i] = key
                    max_vals[i] = ds_mean[key].values.item()
        data_dict = {
            "Predicted Error": np.array([]),
            "Input Parameter": np.array([]),
            "Not Perturbed Value": np.array([]),
            "time_after_ascent": np.array([]),
        }
        for col in in_params:
            data_dict["Predicted Error"] = np.append(
                data_dict["Predicted Error"], ds_traj[col].values
            )
            data_dict["Input Parameter"] = np.append(
                data_dict["Input Parameter"], np.repeat(col, len(ds_traj[col].values))
            )
            data_dict["Not Perturbed Value"] = np.append(
                data_dict["Not Perturbed Value"], ds_traj[out_p].values
            )
            data_dict["time_after_ascent"] = np.append(
                data_dict["time_after_ascent"], ds_traj["time_after_ascent"].values
            )
        data_dict["Predicted Squared Error"] = data_dict["Predicted Error"] ** 2
        data_dict["Output Parameter"] = np.repeat(
            out_p, len(data_dict["Predicted Error"])
        )
        df_tmp = pd.DataFrame(data=data_dict)

        if (
            np.min(df_tmp["Predicted Squared Error"]) == 0
            and np.max(df_tmp["Predicted Squared Error"]) == 0
        ):
            # nothing to see here
            continue
        print(f"Plotting for {out_p}")
        print(df_tmp.columns)
        plot_time_evolution(
            df=df_tmp,
            backend=args.backend,
            store_path=args.store_path,
            title=args.title + " for " + param_title_names[out_p],
            xlabel=args.xlabel,
            ylabel=args.ylabel
            + " "
            + parse_word(out_p)
            + " "
            + get_unit(out_p, brackets=True),
            twinlabel=args.twinlabel + " " + get_unit(out_p, brackets=True),
            logy=args.logy,
            width=args.width,
            height=args.height,
            logtwin=args.logtwin,
            min_x=args.min_time,
            max_x=args.max_time,
        )


def parse_perturbation_plot(ds, out_params, param_title_names, args):
    """

    Parameters
    ----------
    ds
    out_params
    param_title_names
    args

    Returns
    -------

    """
    df_mean = (
        ds.mean(dim=["time_after_ascent"], skipna=True).to_dataframe().reset_index()
    )
    df = ds.to_dataframe().reset_index()
    df = df.loc[df["trajectory"] == args.traj]
    df_mean = df_mean.loc[df_mean["trajectory"] == args.traj]
    for out_p in out_params:
        df_tmp = df.loc[df["Output Parameter"] == out_p]
        df_mean_tmp = df_mean.loc[df_mean["Output Parameter"] == out_p]
        if len(args.in_parameter) == 0:
            in_params = list(
                np.unique(
                    df_mean_tmp.nlargest(
                        args.n_model_params, "Predicted Squared Error"
                    )["Input Parameter"]
                )
            )
        else:
            in_params = args.in_parameter

        df_tmp = df_tmp.loc[df_tmp["Input Parameter"].isin(in_params)]
        if (
            np.min(df_tmp["Predicted Squared Error"]) == 0
            and np.max(df_tmp["Predicted Squared Error"]) == 0
        ):
            # nothing to see here
            continue
        print(f"Plotting for {out_p}")
        print(df_tmp.columns)
        df_tmp["Not Perturbed Value"] = df_tmp["Not Perturbed Value"] * -1
        df_tmp["Predicted Error"] = df_tmp["Predicted Error"] * -1
        plot_time_evolution(
            df=df_tmp,
            backend=args.backend,
            store_path=args.store_path,
            title=args.title + " for " + param_title_names[out_p],
            xlabel=args.xlabel,
            ylabel=args.ylabel
            + " "
            + parse_word(out_p)
            + " "
            + get_unit(out_p, brackets=True),
            twinlabel=args.twinlabel + " " + get_unit(out_p, brackets=True),
            logy=args.logy,
            width=args.width,
            height=args.height,
            logtwin=args.logtwin,
            min_x=args.min_time,
            max_x=args.max_time,
            save=True,
        )


def parse_correlation_plot(ds, out_params, param_title_names, args):
    """

    Parameters
    ----------
    ds
    out_params
    param_title_names
    args

    Returns
    -------

    """
    for out_p in out_params:
        df = (
            ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        df = df.loc[df["Output Parameter"] == out_p]

        if len(args.in_parameter) == 0:
            in_params = list(np.unique(df["Input Parameter"]))
        else:
            in_params = args.in_parameter
        df = df.loc[df["Input Parameter"].isin(in_params)]

        if (
            np.min(df["Predicted Squared Error"]) == 0
            and np.max(df["Predicted Squared Error"]) == 0
        ):
            # nothing to see here; Usually applies to
            # variables with no contents
            continue
        print(f"Plotting for {out_p}")
        if args.set_zero is not None:
            print("Replacing following parameters and values with zero:")
            print(
                df.loc[
                    (df["Predicted Squared Error"] <= args.set_zero)
                    & (df["Predicted Squared Error"] != 0)
                ][
                    [
                        "Input Parameter",
                        "Predicted Squared Error",
                        "Mean Squared Error",
                    ]
                ]
            )
            df["Predicted Squared Error"] = df["Predicted Squared Error"].where(
                df["Predicted Squared Error"] > args.set_zero, 0.0
            )
        if args.plot_variant == "histogram":
            plot_histogram(
                df=df.loc[df["Output Parameter"] == out_p],
                store_path=args.store_path,
                backend=args.backend,
                plot_types=args.plot_types,
                add_zero_sens=args.add_zero_sens,
                title=args.title + r" for " + param_title_names[out_p],
                xlabel=args.xlabel,
                xlabel2=args.ylabel,
                width=args.width,
                height=args.height,
            )
        else:
            hist = False
            if "correlation_hist" == args.plot_variant:
                hist = True
            plot_mse(
                df=df,
                out_params=[out_p],
                store_path=args.store_path,
                backend=args.backend,
                plot_types=args.plot_types,
                add_zero_sens=args.add_zero_sens,
                confidence=args.confidence,
                title=args.title + r" for " + param_title_names[out_p],
                xlabel=args.xlabel,
                ylabel=args.ylabel,
                width=args.width,
                height=args.height,
                hist=hist,
                legend_pos=args.legend_pos,
                corr_line=args.corr_line,
            )


def parse_correlation_all_plot(ds, args):
    """

    Parameters
    ----------
    ds
    args

    Returns
    -------

    """
    all_df = None
    out_params = []
    for out_p in ds["Output Parameter"]:
        out_params.append(out_p.item())
    for out_p in out_params:
        df = (
            ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
            .to_dataframe()
            .reset_index()
        )
        df = df.loc[df["Output Parameter"] == out_p]

        if len(args.in_parameter) == 0:
            in_params = list(np.unique(df["Input Parameter"]))
        else:
            in_params = args.in_parameter
        df = df.loc[df["Input Parameter"].isin(in_params)]

        if (
            np.min(df["Predicted Squared Error"]) == 0
            and np.max(df["Predicted Squared Error"]) == 0
        ):
            # nothing to see here; Usually applies to
            # variables with no contents
            continue
        if all_df is None:
            all_df = df
        else:
            all_df = all_df.append(df)

    hist = False
    if "correlation_hist" == args.plot_variant:
        hist = True
    plot_mse(
        df=all_df,
        out_params=out_params,
        store_path=args.store_path,
        backend=args.backend,
        plot_types=args.plot_types,
        add_zero_sens=args.add_zero_sens,
        confidence=args.confidence,
        title=args.title,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        width=args.width,
        height=args.height,
        hist=hist,
        plot_kind="single_plot",
    )


def main(args):
    """

    Parameters
    ----------
    args

    Returns
    -------

    """
    ds = xr.open_dataset(args.data_path, decode_times=False, engine="netcdf4")
    if len(args.out_parameter) == 0:
        out_params = []
        for out_p in ds["Output Parameter"]:
            out_params.append(out_p.item())
    else:
        out_params = args.out_parameter

    param_title_names = {
        "QV": "Water Vapor",
        "QC": "Cloud Mass",
        "QR": "Rain Mass",
        "QS": "Snow Mass",
        "QI": "Ice Mass",
        "QG": "Graupel Mass",
        "QH": "Hail Mass",
        "NCCLOUD": "Cloud Number",
        "NCRAIN": "Rain Number",
        "NCSNOW": "Snow Number",
        "NCICE": "Ice Number",
        "NCGRAUPEL": "Graupel Number",
        "NCHAIL": "Hail Number",
        "pressure": "Pressure",
        "QR_OUT": "QR_OUT",
        "QS_OUT": "QS_OUT",
        "QI_OUT": "QI_OUT",
        "QG_OUT": "QG_OUT",
        "QH_OUT": "QH_OUT",
        "NR_OUT": "Precipitation of Rain Droplets",
        "NS_OUT": "NS_OUT",
        "NI_OUT": "NI_OUT",
        "NG_OUT": "NG_OUT",
        "NH_OUT": "NH_OUT",
        "latent_heat": "Latent Heating",
        "latent_cool": "Latent Cooling",
    }

    if "correlation" in args.plot_variant or args.plot_variant == "histogram":
        if out_params[0] == "all_at_once":
            parse_correlation_all_plot(ds, args)
        else:
            parse_correlation_plot(ds, out_params, param_title_names, args)
    elif args.plot_variant == "time_plot":
        # Result from a perturbation
        if "time_after_ascent" in ds.dims:
            parse_perturbation_plot(ds, out_params, param_title_names, args)
        else:
            # Result from a sensitivity simulation
            parse_sensitivity_plot(ds, out_params, param_title_names, args)
    else:
        print(f"plot_variant '{args.plot_variant}': No such plot variant. ABORTING!")


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Plot either mean squared deviation/error from perturbation over mean
            predicted deviation/error calculated via the sensitivity where the predicted
            axis is at most 1 such that plots for particle numbers are not entirely
            correct. Or plot the model state variable and predicted squared error
            over time.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        type=str,
        required=True,
        help=textwrap.dedent(
            """\
            Path to dataframe as NetCDF-file which had been created via create_mse.py.
            """
        ),
    )
    parser.add_argument(
        "--add_zero_sens",
        action="store_true",
        help=textwrap.dedent(
            """\
            Add sensitivities of value zero to the far left and mark it with
            negative infinity.
            """
        ),
    )
    parser.add_argument(
        "--plot_types",
        action="store_true",
        help=textwrap.dedent(
            """\
            If true: Plot input parameter types in different colors.
            Types means here if a parameter is physical or rather artificial.
            """
        ),
    )
    parser.add_argument(
        "--out_parameter",
        type=str,
        nargs="+",
        default=[],
        help=textwrap.dedent(
            """\
            Output parameter to plot for. Default plots all that are available
            in the dataframe in a separate plot. If you want a plot with all
            datapoints in one plot, use "all_at_once".
            """
        ),
    )
    parser.add_argument(
        "--backend",
        default="matplotlib",
        help=textwrap.dedent(
            """\
            Choose a backend for plotting. Options are:
            matplotlib: Most plots should be fine with it.
            bokeh: Recommended.
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        default="../pics/correlation",
        type=str,
        help=textwrap.dedent(
            """\
            Path to store the generated images.
            """
        ),
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=None,
        help=textwrap.dedent(
            """\
            Plot a confidence ellipse around each sample with confidence
            between 0 and 1. If none is given, no ellipse will be plotted.
            """
        ),
    )
    parser.add_argument(
        "--xlabel",
        default="Predicted Log MSD",
        type=str,
        help=textwrap.dedent(
            """\
            Alternative label for x-axis.
            """
        ),
    )
    parser.add_argument(
        "--ylabel",
        default="True Log MSD",
        type=str,
        help=textwrap.dedent(
            """\
            Alternative label for y-axis. If plot_variant is "time_plot", then
            " [out_parameter]" is added.
            """
        ),
    )
    parser.add_argument(
        "--title",
        default="True Deviation vs Prediction",
        type=str,
        help=textwrap.dedent(
            """\
            Title for the plot where " for [out_param]" is added.
            """
        ),
    )
    parser.add_argument(
        "--width",
        default=900,
        type=int,
        help=textwrap.dedent(
            """\
            Width of plot in pixels.
            """
        ),
    )
    parser.add_argument(
        "--height",
        default=900,
        type=int,
        help=textwrap.dedent(
            """\
            Height of plot in pixels.
            """
        ),
    )
    parser.add_argument(
        "--set_zero",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            If plot_variant is "correlation".
            Set any predicted squared errors with this value or lower to zero.
            This makes the plots easier to look at, when only a single parameter
            has a predicted error of 1e-200 or less.
            """
        ),
    )
    parser.add_argument(
        "--plot_variant",
        default="correlation",
        type=str,
        help=textwrap.dedent(
            """\
            Plot either correlation plots with true deviation over predicted
            deviation by perturbing a parameter with "correlation",
            "correlation_hist" to add histograms on each axis,
            "histogram" for
            plotting the histogram of true and predicted deviations or
            use "time_plot" to plot (a single or mean) trajectory and
            the predicted deviation over time with the actual model
            state variable.
            """
        ),
    )
    parser.add_argument(
        "--traj",
        type=int,
        default=-1,
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", the trajectory with this index will
            be plotted. If a value below zero is given, plot the mean of all.
            """
        ),
    )
    parser.add_argument(
        "--in_parameter",
        type=str,
        nargs="+",
        default=[],
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", then plot the predicted deviation
            for those model parameters. If none are given, plot the top ten
            most influential parameters for each model state parameter.
            This plots all those predictions in one plot.
            """
        ),
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", plot the y-axis as log10.
            """
        ),
    )
    parser.add_argument(
        "--twinlabel",
        default="Predicted Squared Error",
        type=str,
        help=textwrap.dedent(
            """\
            Only if plot_type is "time_plot". Label for the twin axis.
            """
        ),
    )
    parser.add_argument(
        "--logtwin",
        action="store_true",
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", plot the twin-axis as log10.
            """
        ),
    )
    parser.add_argument(
        "--n_model_params",
        type=int,
        default=5,
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", plot this many model parameters.
            """
        ),
    )
    parser.add_argument(
        "--min_time",
        type=float,
        default=None,
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", use this as start point for the plot
            as in time after ascent.
            """
        ),
    )
    parser.add_argument(
        "--max_time",
        type=float,
        default=None,
        help=textwrap.dedent(
            """\
            If plot_type is "time_plot", use this as last point for the plot
            as in time after ascent.
            """
        ),
    )
    parser.add_argument(
        "--legend_pos",
        type=str,
        default="bottom_right",
        help=textwrap.dedent(
            """\
            Define the position of the legend for most plots.
            """
        ),
    )
    parser.add_argument(
        "--corr_line",
        action="store_true",
        help=textwrap.dedent(
            """\
            Add a dashed line for a 1-to-1 map of the data.
            """
        ),
    )
    main(parser.parse_args())
