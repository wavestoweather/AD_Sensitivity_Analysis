from joblib import load
import pandas as pd
import xarray as xr
import hvplot.xarray  # noqa
import hvplot.pandas  # noqa

try:
    from segment_identifier import *
except:
    from scripts.segment_identifier import *


def clean_predictions(data, step_tol):
    """
    We iterate over the prediction. The first predicted segment
    makes the next predictions within the step tolerance invalid.
    This is just an approximation to clean the data.
    This is useful to get an idea of the distribution of
    true and false positive for comparison.

    Parameters
    ----------
    data : np.ndarray
        Array of shape (n_samples, n_classes)
    step_tol : int
        Number of steps to tolerate for a prediction to be true.

    Returns
    -------
    np.ndarray from input but cleaned from duplicate predictions.
    """
    window_idx = -1
    y_idx = np.argwhere(data)
    y_idx = sorted(y_idx, key=lambda x: x[1])
    replace_idx = None
    for i, y_i in enumerate(y_idx):
        if i > 0:
            if y_i[1] == y_idx[i - 1][1]:
                if y_i[0] - window_idx >= step_tol * 2 + 1:
                    # new window with same parameter
                    window_idx = y_i[0]
                    y_idx[i]
                elif replace_idx is None:
                    replace_idx = y_i
                    break
            else:
                # new parameter
                window_idx = y_i[0]
        else:
            window_idx = y_i[0]

    if replace_idx is None:
        # Nothing to clean up
        return data

    window_idx = -1
    # Get the index, where the first segment is detected
    # delete it and leave the other in that window
    # Remove those starts given the left overs
    y_idx2 = y_idx.copy()
    for i, y_i in enumerate(y_idx):
        if i > 0:
            if y_i[1] == y_idx[i - 1][1]:
                if y_i[0] - window_idx >= step_tol * 2 + 1:
                    # new window with same parameter
                    window_idx = y_i[0]
                    y_idx[i]
                    y_idx2[i] = replace_idx
            else:
                # new parameter
                window_idx = y_i[0]
                y_idx2[i] = replace_idx
        else:
            window_idx = y_i[0]
            y_idx2[i] = replace_idx
    data[y_idx2] = False
    return data


def plot_segment_hist(
    ds,
    y,
    store_path=None,
    width=1680,
    height=1050,
    bins=30,
    group=None,
    stacked=False,
    alpha=1.0,
):
    """
    Plot a histogram of segment starts.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset created by create_hist_seg_dataset(..).
    y : string
        Column for histogram where a segment start occurs.
    store_path : string
        Path and name where to save images.
    width : int
        Width in pixels for the plot.
    height : int
        Height in pixels for the plot.
    bins : int
        Number of bins.
    group : string
        Column to group by.
    stacked : boolean
        Plot a stacked histogram. Only useful with "group".
    alpha : float
        Alpha value for the histogram.

    Returns
    -------
    Histogram plot.
    """
    # Since "Input Parameter" is not a numerical value, we create that
    # using bar plot; hvplot.hist does not support stacking.
    if y == "Input Parameter" or stacked:
        if group is None:
            tmp_dic = {
                y: [],
                "segment_starts": [],
            }
            for i in np.unique(ds[y]):
                tmp_dic[y].append(i)
                tmp_dic["segment_starts"].append(np.sum((ds[y] == i)))
            df = pd.DataFrame.from_dict(tmp_dic)
            hist = df.hvplot.bar(
                x=y,
                y="segment_starts",
                width=width,
                height=height,
                alpha=alpha,
            )
        else:
            tmp_dic = {
                y: [],
                "segment_starts": [],
                group: [],
            }
            for i in np.unique(ds[y]):
                for j in np.unique(ds[group]):
                    tmp_dic[y].append(i)
                    tmp_dic["segment_starts"].append(
                        np.sum((ds[y] == i) & (ds[group] == j))
                    )
                    tmp_dic[group].append(j)
            df = pd.DataFrame.from_dict(tmp_dic)
            hist = df.hvplot.bar(
                x=y,
                y="segment_starts",
                by=group,
                width=width,
                height=height,
                alpha=alpha,
                stacked=stacked,
            )
    else:
        hist = ds.hvplot.hist(
            y=y,
            bins=bins,
            by=group,
            width=width,
            height=height,
            alpha=alpha,
        )
    if store_path is not None:
        renderer = hv.Store.renderers["bokeh"].instance(fig="png", dpi=300)
        renderer.save(hist, store_path)
    return hist


def create_hist_seg_dataset(
    data,
    all_params_list,
    model_path=None,
    step_tol=0,
    remove_pred_duplicates=True,
    verbosity=0,
):
    """
    Create a xarray.Dataset useful for plotting histograms. Indexes where
    segment starts are, are stored here.

    Parameters
    ----------
    data : list or np.ndarray or xarray.Dataset
        Either a precalculated set of data that had been stored as *.npy
        or a dataset where windows (if any) need to be calculated.
    all_params_list : list of string
        List of input parameters.
    model_path : string
        Path to model saved as *.joblid.
    step_tol : int
        Number of steps to tolerate for a prediction to be true.
    remove_pred_duplicates : boolean
        If true: Remove segment starts that are predicted more than once.
        Only useful if a model is given and "step_tol" > 1.
    verbosity : int
        Set verbosity level.
        0: No output except for exceptions
        1: Print datasets
        2: Print loading statements
        3: Print building/training statements
        4: Set verbosity of random forest to 2
        5: Get predicted results for each entry
        6: Get count for segment starts for every input parameter

    Returns
    -------
    xarray.Dataset either with columns "trajectory", "time_after_ascent"
    and "Input Parameter" or "datapoint" and "Input Parameter" if data
    is a list or np.ndarray.
    """
    if model_path is not None:
        model = load(model_path)
        if verbosity > 3 and "adaboost" not in model_path:
            model = model.set_params(**{"estimator__verbose": 2})
        elif "adaboost" not in model_path:
            model = model.set_params(**{"estimator__verbose": 0})

        if isinstance(data, list) or isinstance(data, np.ndarray):
            # Step tolerance is calculated from the dataset
            # step_tol = (np.shape(data)[1]//(len(out_params)*len(all_params_list)) - 1 ) // 2
            # data is either already precalculated set of windows
            # new_shape = (np.shape(data)[0],) + (len(all_params_list), len(out_params), step_tol*2+1)
            data = model.predict(data)
            # data = np.reshape(data, new_shape)
            # create a dataset such that we can plot it later
            data = xr.Dataset.from_dict(
                {
                    "index": {"dims": ("index"), "data": np.arange(np.sum(data))},
                    "datapoint": {"dims": ("index"), "data": np.argwhere(data)[:, 0]},
                    "Input Parameter": {
                        "dims": ("index"),
                        "data": all_params_list[np.argwhere(data)[:, 1]],
                    },
                    # "Output Parameter": {"dims": ("Output Parameter"), "data": out_params[np.argwhere(data)[:, 2]]},
                    # "window_idx": {"dims": ("window_idx"), "data": np.argwhere(data)[:, 3]},
                }
            )
        else:
            # or windows need to be calculated and predictions need to be made
            new_shape = (
                (
                    len(data["trajectory"]),
                    len(data["time_after_ascent"]) - step_tol,
                    len(all_params_list),
                )
                # len(out_params),
                # args.step_tol*2+1,)
            )
            time_after_ascent = data["time_after_ascent"]
            if "segment_start" not in data:
                data = find_segments(data, 10.0 ** -10.0)
            data, _, _, _ = create_input_labels(
                ds=data,
                step_tol=step_tol,
                distinct_outparams=False,
                verbosity=verbosity,
            )
            data = model.predict(data)
            if remove_pred_duplicates:
                data = clean_predictions(data, step_tol)
            data = np.reshape(data, new_shape)

            data = xr.Dataset.from_dict(
                {
                    "index": {"dims": ("index"), "data": np.arange(np.sum(data))},
                    "trajectory": {"dims": ("index"), "data": np.argwhere(data)[:, 0]},
                    "time_after_ascent": {
                        "dims": ("index"),
                        "data": time_after_ascent[np.argwhere(data)[:, 1]],
                    },
                    "Input Parameter": {
                        "dims": ("index"),
                        "data": all_params_list[np.argwhere(data)[:, 2]],
                    },
                    # "Output Parameter": {"dims": ("Output Parameter"), "data": all_params_list[np.argwhere(data)[:, 3]]},
                    # "window_idx": {"dims": ("window_idx"), "data": np.argwhere(data)[:, 4]},
                }
            )
    else:
        if isinstance(data, list) or isinstance(data, np.ndarray):
            # just labels where we have no information about the order of the
            # datapoints
            data = xr.Dataset.from_dict(
                {
                    "index": {"dims": ("index"), "data": np.arange(np.sum(data))},
                    "datapoint": {"dims": ("index"), "data": np.argwhere(data)[:, 0]},
                    "Input Parameter": {
                        "dims": ("index"),
                        "data": all_params_list[np.argwhere(data)[:, 1]],
                    },
                }
            )
        else:
            # windows need to be calculated
            new_shape = (
                (
                    len(data["trajectory"]),
                    len(data["time_after_ascent"]) - step_tol,
                    len(all_params_list),
                )
                # len(out_params),
                # args.step_tol*2+1,)
            )
            time_after_ascent = data["time_after_ascent"]
            if "segment_start" not in data:
                data = find_segments(data, 10.0 ** -10.0)
            _, data, _, _ = create_input_labels(
                ds=data,
                step_tol=step_tol,
                distinct_outparams=False,
                verbosity=verbosity,
            )
            if remove_pred_duplicates:
                data = clean_predictions(data, step_tol)
            data = np.reshape(data, new_shape)

            data = xr.Dataset.from_dict(
                {
                    "index": {"dims": ("index"), "data": np.arange(np.sum(data))},
                    "trajectory": {"dims": ("index"), "data": np.argwhere(data)[:, 0]},
                    "time_after_ascent": {
                        "dims": ("index"),
                        "data": time_after_ascent[np.argwhere(data)[:, 1]],
                    },
                    "Input Parameter": {
                        "dims": ("index"),
                        "data": all_params_list[np.argwhere(data)[:, 2]],
                    },
                    # "Output Parameter": {"dims": ("Output Parameter"), "data": all_params_list[np.argwhere(data)[:, 3]]},
                    # "window_idx": {"dims": ("window_idx"), "data": np.argwhere(data)[:, 4]},
                }
            )
    return data


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        """
        Plot different histograms of the training and test data for
        segment prediction. Plots are segment start which can be determined
        and grouped by trajectory,
        time_after_ascent, input_parameter and output_parameter.
        If a model path is given, histograms for the prediction are plotted.
        This always uses the bokeh backend of holoviews.
        """,
    )
    parser.add_argument(
        "--data_path",
        required=True,
        type=str,
        help="""
        Path to folders with ensemble datasets or to single NetCDF file
        with all data concatenated along 'trajectory' axis.
        If a path to numpy arrays is given, it is assumed to be a
        stratified training or test set without information regarding
        trajectory or time step.
        """,
    )
    parser.add_argument(
        "--count_column",
        required=True,
        type=str,
        help="""
        Column for histogram where a segment start occurs. Options are
        trajectory: Only if "data_path" does not point to *.npy.
        Input Parameter: Always available
        trajectory: Only if "data_path" does not point to *.npy.
        time_after_ascent: Only if "data_path" does not point to *.npy.
        """,
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=30,
        help="""
        Number of bins for the histogram.
        """,
    )
    parser.add_argument(
        "--group",
        type=str,
        default=None,
        help="""
        Group histogram by this keyword. Possible options are:
        "input_parameter"
        "output_parameter"
        "time_after_ascent": Not available for *.npy files.
        "trajectory": Not available for *.npy files.
        """,
    )
    parser.add_argument(
        "--step_tol",
        type=int,
        default=0,
        help="""
        Plot histograms using sliding window data which increases the total
        amount of (possible) segment starts by taking a tolerance into account.
        This is required if a model is being loaded and "data_path" is not
        pointing to a *.npy file.
        """,
    )
    parser.add_argument(
        "--model_path",
        type=str,
        default=None,
        help="""
        Plot the predicted segment starts given this model saved as *.joblid.
        """,
    )
    parser.add_argument(
        "--store_path",
        type=str,
        default="../pics/segment_histogram.png",
        help="""
        Path and name where to save images.
        """,
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1680,
        help="""
        Width in pixels for the plot.
        """,
    )
    parser.add_argument(
        "--height",
        type=int,
        default=1080,
        help="""
        Height in pixels for the plot.
        """,
    )
    parser.add_argument(
        "--stacked",
        action="store_true",
        help="""
        Plot a stacked histogram. Only useful with "--group".
        """,
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=1.0,
        help="""
        Alpha value for the histogram.
        """,
    )
    parser.add_argument(
        "--remove_pred_dupl",
        action="store_true",
        help="""
        If a model is loaded for predicting the labels and the step tolerance
        is bigger than 1, then a single segment start is counted as many times
        as step tolerance allowes. If "remove_pred_dupl" is set, count the
        amount of segments predicted by removing those which are within
        step tolerance and predicted multiple times.
        """,
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help="""
        Set verbosity level.
        0: No output except for exceptions
        1: Print datasets
        2: Print loading statements
        3: Print building/training statements
        4: Set verbosity of random forest to 2
        5: Get predicted results for each entry
        6: Get count for segment starts for every input parameter
        """,
    )

    args = parser.parse_args()
    # "trajectory" is only available for not precalculated data
    if (
        args.group == "time_after_ascent" or args.group == "trajectory"
    ) and ".npy" in args.path:
        print(
            f"Group by operation for {args.group} is not available for precalculated datasets."
        )
        print("ABORTING")
        exit()

    all_params_list = np.asarray(
        [
            "dc_ccn_4",
            "dice_b_geo",
            "da_prime",
            "dT_mult_max",
            "dsnow_b_f",
            "dsnow_c_s",
            "dsnow_b_vel",
            "dsnow_a_f",
            "dgraupel_vsedi_min",
            "dT_mult_min",
            "dgraupel_ecoll_c",
            "db_prime",
            "dice_c_s",
            "db_v",
            "dp_sat_melt",
            "drain_b_geo",
            "dhail_b_vel",
            "dsnow_b_geo",
            "dD_rainfrz_gh",
            "db_HET",
            "dgraupel_b_f",
            "dgraupel_b_vel",
            "drain_beta",
            "dice_b_f",
            "drain_mu",
            "dice_b_vel",
            "db_ccn_4",
            "dcloud_b_vel",
            "dsnow_a_geo",
            "dp_sat_ice_const_b",
            "da_HET",
            "drain_c_z",
            "drain_b_vel",
            "drain_a_geo",
            "dhail_b_geo",
            "dgraupel_c_s",
            "dgraupel_b_geo",
            "dice_a_geo",
            "drain_gamma",
            "dinv_z",
            "drain_a_vel",
            "dhail_a_vel",
            "dcloud_b_geo",
            "drain_alpha",
            "drain_cmu2",
            "drho_vel",
            "dhail_a_geo",
            "dgraupel_a_geo",
            "dsnow_a_vel",
            "drain_cmu4",
            "dkin_visc_air",
            "dk_r",
            "drain_cmu3",
            "dice_a_f",
            "drain_nu",
            "dgraupel_max_x",
            "drain_g2",
            "dice_a_vel",
            "dgraupel_a_vel",
        ]
    )
    out_params = np.asarray(["QV", "QC", "QR", "QG", "QH", "QI", "QS"])

    if not "labels" in args.data_path and args.model_path is None:
        print("'data_path' must contain labels if no model for prediction is given")
        print("ABORTING")
        exit()
    if "labels" in args.data_path and args.model_path is not None:
        print("'data_path' must contain input data and not labels for prediction")
        print("ABORTING")
        exit()

    if ".npy" in args.data_path:
        data = np.load(args.data_path, allow_pickle=True, fix_imports=False)
    else:
        data = parse_load(
            data_path=args.data_path,
            out_params=out_params,
            all_params_list=all_params_list,
            store_many_appended_data=False,
            load_on_the_fly=False,
            verbosity=args.verbosity,
        )

    data = create_hist_seg_dataset(
        data=data,
        all_params_list=all_params_list,
        model_path=args.model_path,
        step_tol=args.step_tol,
        remove_pred_duplicates=args.remove_pred_dupl,
        verbosity=args.verbosity,
    )
    # Plot the histogram
    _ = plot_segment_hist(
        ds=data,
        y=args.count_column,
        group=args.group,
        store_path=args.store_path,
        width=args.width,
        height=args.height,
        bins=args.bins,
        stacked=args.stacked,
        alpha=np.min(1.0, np.max(0.0, args.alpha)),
    )
