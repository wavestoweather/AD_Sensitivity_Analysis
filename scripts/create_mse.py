from multiprocessing import Pool
import numpy as np
import os
import pandas as pd

try:
    from tqdm import tqdm
except:
    from progressbar import progressbar as tqdm
import sys
import xarray as xr

try:
    from Deriv_dask import Deriv_dask
    from latexify import in_params_dic, physical_params, in_params_grouping
    from segment_identifier import d_unnamed
except:
    from scripts.Deriv_dask import Deriv_dask
    from scripts.latexify import in_params_dic, physical_params, in_params_grouping
    from scripts.segment_identifier import d_unnamed


def load_and_append(name_list, filetype="csv"):
    """
    Load and append multiple files that are temporary results created
    by create_mse.py.

    Parameters
    ----------
    name_list : list of string
        List of paths for files to load and append together.
    filetype : string
        Define which type to load. Options are "csv" and "netcdf".

    Returns
    -------
    pandas.Dataframe with all data appended together.
    """
    all_df = None
    for name in name_list:
        try:
            if all_df is None:
                if filetype == "csv":
                    all_df = d_unnamed(pd.read_csv(name))
                else:
                    all_df = xr.open_dataset(
                        name, decode_times=False, engine="netcdf4"
                    ).to_dataframe()
            else:
                if filetype == "csv":
                    all_df = all_df.append(d_unnamed(pd.read_csv(name)))
                else:
                    all_df = all_df.append(
                        xr.open_dataset(
                            name, decode_times=False, engine="netcdf4"
                        ).to_dataframe()
                    )
        except:
            pass
    return all_df


def get_errors_by_file(
    base_path,
    n_trajs,
    error_kind,
    df_name,
    store_path,
    filetype,
    n_cpus,
    ratio_type,
    sens_kind,
    save_intermediate,
):
    """
    Load multiple files, for each trajectory one, with perturbed ensembles
    and get the error and predicted error. All trajectories originate from
    a single source file that had been used as input for the simulation in
    this package.

    Parameters
    ----------
    base_path : path
        Path to perturbed ensembles. Somewhere must be a file with
        "_notPerturbed.nc_wcb" in its name.
    n_trajs : int
        Number of trajectories in the given path for which ensembles have been
        started.
    error_kind : string
        How to reduce the error from ensembles.
        mse: Calculate mean squared error
        maxse: Calculate maximum squared error
        nozeromse: Calculate mean squared error ignoring time steps
            where all simulations have zero values.
        sum : Cumulative squared error
        me : Mean error
        mae : Mean absolute error
    df_name : string
        Name of the original file that had been used in the original
        simulation and from which all trajectories to process here originate
        from.
    store_path : path
        If "save_intermediate" is true, save the results in this path.
    filetype : string
        If "save_intermediate" is true, save the results from this method
        with the given filetype "csv" or "netcdf".
    n_cpus : int
        Number of processes that process data in parallel.
    ratio_type : string
        "vanilla": Use the derivative ratio in the file.
        "adjusted": Can be added to any to any type below where the sensitivity
        is adjusted to the parameter value such that perturbing this value by a
        certain percentage gives an approximation of the resulting error/difference
        for any given hydrometeor.
        "per_timestep": Use the highest derivative per timestep as denominator.
        "window": Use the highest derivative in the given window by min_x and max_x.
        "per_xaxis": Use the highest derivative per x_axis value. If x_axis is "timestep"
        it is the same as "per_timestep".
        "x_per_out_param": Replace 'x' with any other option than "vanilla". Use the highest
        derivative but per output parameter. (that *should* be the vanilla version)
        "x_weighted": Append this to any ratio type to use the inverse of the
        output parameter value as weight. Only works with "x_per_out_param".
    sens_kind : string

    save_intermediate : bool
        If true, save the results here into "store_path"
        in a given filetype.

    Returns
    -------
    pandas.Dataframe with all trajectories and their errors merged together.
    """
    # Ugly workaround to pickle load_mse(..)
    # I could move this function outside of this method though
    global load_mse
    out_params = [
        "QV",
        "QC",
        "QR",
        "QG",
        "QH",
        "QI",
        "QS",
        "NC",
        "NR",
        "NG",
        "NH",
        "NI",
        "NS",
    ]
    in_params = []
    for key in in_params_grouping:
        if key == "1-moment":
            continue
        in_params.extend(in_params_grouping[key])
    for e in physical_params:
        if e in in_params:
            in_params.remove(e)

    def load_mse(i):
        try:
            mean_file = "traj" + str(i) + "_notPerturbed.nc_wcb"
            if not os.path.isfile(base_path + mean_file):
                mean_file = "traj" + str(i) + "/_notPerturbed.nc_wcb"
            others_path = base_path + "traj" + str(i) + "/"
            if not os.path.isfile(base_path + mean_file):
                mean_file = "_notPerturbed.nc_wcb"
                others_path = base_path
            mean_traj = Deriv_dask(
                direc=base_path,
                parquet=False,
                netcdf=True,
                columns=None,
                backend="matplotlib",
                file_ending=mean_file,
            )
            mean_traj.cache_data(
                in_params=in_params,
                out_params=out_params,
                x_axis="time_after_ascent",
                compute=True,
            )
            df = mean_traj.get_mse(
                out_params=out_params,
                others_path=others_path,
                in_params=in_params,
                ratio_type=ratio_type,
                kind=error_kind,
                sens_kind=sens_kind,
            )
            return df
        except:
            return None

    with Pool(n_cpus) as p:
        mse_df = pd.concat(p.map(load_mse, range(n_trajs)))
    if not save_intermediate:
        return mse_df
    if filetype == "csv":
        mse_df.to_csv(store_path + error_kind + "_" + df_name)
    else:
        mse_ds = mse_df.to_xarray()
        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in mse_ds.data_vars}
        try:
            mse_ds.to_netcdf(
                path=store_path + error_kind + "_" + df_name + ".nc",
                encoding=encoding,
                compute=True,
                engine="netcdf4",
                format="NETCDF4",
                mode="w",
            )
        except:
            print("Unable to write file. Trying without encoding")
            mse_ds.to_netcdf(
                path=store_path + error_kind + "_" + df_name + ".nc",
                compute=True,
                engine="netcdf4",
                format="NETCDF4",
                mode="w",
            )
    return mse_df


def reduce_errors(
    mse_df,
    error_kind,
    df_name,
    reduce_func,
    store_path,
    filetype,
):
    """
    Given a dataframe with errors and predicted errors from multiple
    trajectories, reduce those errors along the trajectory dimension such
    that only one datapoint for each input and output parameter is left.

    Parameters
    ----------
    mse_df : pandas.Dataframe
        Dataframe created by get_errors_by_file(..).
    error_kind : string
        How to reduce the error from ensembles.
        mse: Calculate mean squared error
        maxse: Calculate maximum squared error
        nozeromse: Calculate mean squared error ignoring time steps
            where all simulations have zero values.
        sum : Cumulative squared error
        me : Mean error
        mae : Mean absolute error
    df_name : string
        Name of the original file that had been used in the original
        simulation and from which all trajectories to process here originate
        from.
    reduce_func : string
        Function for reducing each error from all trajectories. Possible
        options are:
        max: Max error
        mean: Mean error
        min: Min error
    store_path : path
        Path where to store the final dataframe.
    filetype : string
        If "save_intermediate" is true, save the results from this method
        with the given filetype "csv" or "netcdf".
    """
    output_params = np.unique(mse_df["Output Parameter"])
    ratio_types = np.unique(mse_df["Ratio Type"])
    model_params = np.unique(mse_df["Perturbed Parameter"])
    df_dic = {
        "MSE": [],
        "Sensitivity": [],
        "Output Parameter": [],
        "Perturbed Parameter": [],
        "Ratio Type": [],
    }
    if reduce_func == "max":
        cum_f = np.max
    elif reduce_func == "mean":
        cum_f = np.mean
    elif reduce_func == "min":
        cum_f = np.min

    for rt in ratio_types:
        for op in output_params:
            for i in tqdm(range(len(model_params))):
                mp = model_params[i]
                tmp_df = mse_df.loc[
                    (mse_df["Ratio Type"] == rt)
                    & (mse_df["Output Parameter"] == op)
                    & (mse_df["Perturbed Parameter"] == mp)
                ]
                df_dic["MSE"].append(cum_f(tmp_df["MSE"]))
                df_dic["Sensitivity"].append(cum_f(tmp_df["Sensitivity"]))
                df_dic["Output Parameter"].append(op)
                df_dic["Perturbed Parameter"].append(mp)
                df_dic["Ratio Type"].append(rt)
    mean_df = pd.DataFrame.from_dict(df_dic)
    if filetype == "csv":
        mean_df.to_csv(store_path + error_kind + "_" + reduce_func + "_" + df_name)
    else:
        mean_df = mean_df.to_xarray()
        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in mean_df.data_vars}
        try:
            mean_df.to_netcdf(
                path=store_path
                + error_kind
                + "_"
                + reduce_func
                + "_"
                + df_name
                + ".nc",
                encoding=encoding,
                compute=True,
                engine="netcdf4",
                format="NETCDF4",
                mode="w",
            )
        except:
            print("Unable to write file. Trying without encoding")
            mean_df.to_netcdf(
                path=store_path
                + error_kind
                + "_"
                + reduce_func
                + "_"
                + df_name
                + ".nc",
                compute=True,
                engine="netcdf4",
                format="NETCDF4",
                mode="w",
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Load ensemble NetCDF-files and calculate the errors
        and reduced predicted errors and store them to disk. This is used to
        prepare data to find any correlations and for plotting.
        In order to find the correct errors, a file with "_notPerturbed.nc_wcb"
        must be present in the given path.
        """,
    )
    parser.add_argument(
        "--error_kind",
        default="mse",
        help="""
        How to reduce the error from ensembles.
        mse: Calculate mean squared error
        maxse: Calculate maximum squared error
        nozeromse: Calculate mean squared error ignoring time steps
            where all simulations have zero values.
        sum : Cumulative squared error
        me : Mean error
        mae : Mean absolute error
        """,
    )
    parser.add_argument(
        "--sens_kind",
        default="squared_mean",
        help="""
        How to reduce the predicted error from sensitivities over all time
        steps along a trajectory.
        max : Max error
        mean : Mean error
        squared_mean : Squared mean error
        """,
    )
    parser.add_argument(
        "--reduce_func",
        default="mean",
        help="""
        Function for reducing each error from all trajectories. Possible
        options are:
        max: Max error
        mean: Mean error
        min: Min error
        """,
    )
    parser.add_argument(
        "--filetype",
        default="netcdf",
        help="""
        Store the result either as comma separated file with "csv" or
        as NetCDF file with "netcdf".
        """,
    )
    parser.add_argument(
        "--path",
        default="/data/project/wcb/netcdf/perturbed_ensembles/conv_600_0_traj_t000000_p001/",
        help="""
        Path to perturbed ensembles. Somewhere must be a file with
        "_notPerturbed.nc_wcb" in its name.
        """,
    )
    parser.add_argument(
        "--store_path",
        default="../data/stats/",
        help="""
        Path to store the resulting statistics.
        """,
    )
    parser.add_argument(
        "--n_trajs",
        type=int,
        required=True,
        help="""
        Number of trajectories in the given path for which ensembles have been
        started.
        """,
    )
    parser.add_argument(
        "--filename",
        type=str,
        required=True,
        help="""
        Name of the file where trajectories and subsequent ensembles are
        stored.
        """,
    )
    parser.add_argument(
        "--save_intermediate",
        action="store_true",
        help="""
        Save intermediate results in a subdirectory "tmp" within the path
        specified by "store_path" and do not remove them.
        """,
    )
    parser.add_argument(
        "--n_cpus",
        type=int,
        default=6,
        help="""
        Number of processes that process data in parallel.
        """,
    )
    parser.add_argument(
        "--ratio_type",
        type=str,
        default="adjusted",
        help="""
        "vanilla": Use the derivative ratio in the file.
        "adjusted": Can be added to any to any type below where the sensitivity
        is adjusted to the parameter value such that perturbing this value by a
        certain percentage gives an approximation of the resulting error/difference
        for any given hydrometeor.
        "per_timestep": Use the highest derivative per timestep as denominator.
        "window": Use the highest derivative in the given window by min_x and max_x.
        "per_xaxis": Use the highest derivative per x_axis value. If x_axis is "timestep"
        it is the same as "per_timestep".
        "x_per_out_param": Replace 'x' with any other option than "vanilla". Use the highest
        derivative but per output parameter. (that *should* be the vanilla version)
        "x_weighted": Append this to any ratio type to use the inverse of the
        output parameter value as weight. Only works with "x_per_out_param".
        """,
    )
    parser.add_argument(
        "--no_reduction",
        action="store_true",
        help="""
        Only save intermediate results and exit without further reduction.
        """,
    )
    parser.add_argument(
        "--load_intermediate_files",
        type=str,
        nargs="+",
        default=[],
        help="""
        In case intermediate files have been created before, you can define
        each separately and load these instead of calculating the errors again.
        """,
    )
    args = parser.parse_args()

    if args.save_intermediate or args.no_reduction:
        save_intermediate = True
        try:
            os.mkdir(args.store_path + "tmp/")
        except:
            pass
    else:
        save_intermediate = False

    if len(args.load_intermediate_files) > 0:
        mse_df = load_and_append(args.load_intermediate_files, args.filetype)
    else:
        mse_df = get_errors_by_file(
            base_path=args.path,
            n_trajs=args.n_trajs,
            error_kind=args.error_kind,
            df_name=args.filename,
            store_path=args.store_path + "tmp/",
            filetype=args.filetype,
            n_cpus=args.n_cpus,
            ratio_type=args.ratio_type,
            sens_kind=args.sens_kind,
            save_intermediate=save_intermediate,
        )

    if not args.no_reduction:
        reduce_errors(
            mse_df=mse_df,
            error_kind=args.error_kind,
            df_name=args.filename,
            reduce_func=args.reduce_func,
            store_path=args.store_path,
            filetype=args.filetype,
        )
