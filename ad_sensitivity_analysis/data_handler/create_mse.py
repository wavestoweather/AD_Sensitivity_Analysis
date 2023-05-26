"""Transform perturbed ensemble simulations for analysis.

"""
import os

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from ad_sensitivity_analysis.data_handler.loader import load_dataset
from ad_sensitivity_analysis.plot.latexify import param_id_map


def _get_n_nans(ds, perturb_names, coord_name, out_ids, out_params):
    """

    Parameters
    ----------
    ds
    perturb_names
    coord_name
    out_ids
    out_params

    Returns
    -------

    """
    n_nans = [-1 for _ in range(len(ds["trajectory"]))]
    for ens_id, param_name in enumerate(perturb_names):
        # In case the number of perturbations is larger than what had been tracked.
        if param_name not in ds:
            continue
        for traj_id in ds["trajectory"].values:
            if coord_name in ds.dims:
                for out_id in out_ids:
                    tmp_n = np.sum(
                        np.isnan(
                            ds.sel(
                                {
                                    "ensemble": ens_id + 1,
                                    "trajectory": traj_id,
                                    coord_name: out_id,
                                }
                            )[out_params[0]]
                        )
                    ).item()
                    if tmp_n < n_nans[traj_id] or n_nans[traj_id] == -1:
                        n_nans[traj_id] = tmp_n
            else:
                tmp_n = np.sum(
                    np.isnan(
                        ds.sel({"ensemble": ens_id + 1, "trajectory": traj_id})[
                            out_params[0]
                        ]
                    )
                ).item()
                if tmp_n < n_nans[traj_id] or n_nans[traj_id] == -1:
                    n_nans[traj_id] = tmp_n
    return n_nans


def _get_ds_intermediate(ds, out_params, perturb_names_avail):
    """

    Parameters
    ----------
    ds
    out_params
    perturb_names_avail

    Returns
    -------

    """
    ds_tmp = (
        ds.isel({"ensemble": np.arange(1, len(ds["ensemble"]))})[out_params]
        - ds.isel({"ensemble": 0})[out_params]
    )

    ds_tmp_mean = ds_tmp.reduce(np.nanmean, dim="time")
    ds_tmp = ds_tmp * ds_tmp
    ds_tmp = ds_tmp.reduce(np.nanmean, dim="time")
    ds_pred = ds.isel({"ensemble": 0})[perturb_names_avail]
    ds_pred_mean = ds_pred.reduce(np.nanmean, dim="time")
    ds_pred = ds_pred * ds_pred
    ds_pred = ds_pred.reduce(np.nanmean, dim="time")
    return ds_tmp, ds_tmp_mean, ds_pred, ds_pred_mean


# pylint: disable=too-many-locals
def _get_error_df(
    ds,
    perturb_names,
    coord_name,
    out_ids,
    out_params,
    n_nans,
    perturb_names_avail,
):
    """

    Parameters
    ----------
    ds
    perturb_names
    coord_name
    out_ids
    out_params
    n_nans
    perturb_names_avail

    Returns
    -------

    """
    ds_tmp, ds_tmp_mean, ds_pred, ds_pred_mean = _get_ds_intermediate(
        ds, out_params, perturb_names_avail
    )
    data_dic = {
        "Mean Squared Error": [],
        "Mean Error": [],
        "Output Parameter": [],
        "Input Parameter": [],
        "Predicted Error": [],
        "Predicted Squared Error": [],
    }

    for out_name, out_id in zip(out_params, out_ids):
        for ens_id, param_name in enumerate(perturb_names):
            # In case the number of perturbations is larger than what had been tracked.
            if param_name not in ds:
                continue
            for traj_id in ds_tmp["trajectory"]:
                if coord_name in ds.dims:
                    n_nans_tmp = np.sum(
                        np.isnan(
                            ds.sel(
                                {
                                    "ensemble": ens_id + 1,
                                    "trajectory": traj_id,
                                    coord_name: out_id,
                                }
                            )[out_params[0]]
                        )
                    )
                else:
                    n_nans_tmp = np.sum(
                        np.isnan(
                            ds.sel({"ensemble": ens_id + 1, "trajectory": traj_id})[
                                out_params[0]
                            ]
                        )
                    )
                if n_nans[traj_id.item()] < n_nans_tmp.item():
                    continue
                data_dic["Mean Squared Error"].append(
                    ds_tmp.sel(
                        {
                            "trajectory": traj_id,
                            "ensemble": ens_id + 1,
                        }
                    )[out_name].values.item()
                )
                data_dic["Mean Error"].append(
                    ds_tmp_mean.sel(
                        {
                            "trajectory": traj_id,
                            "ensemble": ens_id + 1,
                        }
                    )[out_name].values.item()
                )
                if coord_name in ds_pred_mean.dims:
                    if "trajectory" in ds_pred_mean:
                        data_dic["Predicted Error"].append(
                            ds_pred_mean.sel(
                                {"trajectory": traj_id, coord_name: out_id}
                            )[param_name].values.item()
                        )
                        data_dic["Predicted Squared Error"].append(
                            ds_pred.sel({"trajectory": traj_id, coord_name: out_id})[
                                param_name
                            ].values.item()
                        )
                    else:
                        data_dic["Predicted Error"].append(
                            ds_pred_mean.sel({coord_name: out_id})[
                                param_name
                            ].values.item()
                        )
                        data_dic["Predicted Squared Error"].append(
                            ds_pred.sel({coord_name: out_id})[param_name].values.item()
                        )
                else:
                    if "trajectory" in ds_pred_mean:
                        data_dic["Predicted Error"].append(
                            ds_pred_mean.sel({"trajectory": traj_id})[
                                param_name
                            ].values.item()
                        )
                        data_dic["Predicted Squared Error"].append(
                            ds_pred.sel({"trajectory": traj_id})[
                                param_name
                            ].values.item()
                        )
                    else:
                        data_dic["Predicted Error"].append(
                            ds_pred_mean[param_name].values.item()
                        )
                        data_dic["Predicted Squared Error"].append(
                            ds_pred[param_name].values.item()
                        )
                data_dic["Output Parameter"].append(out_name)
                data_dic["Input Parameter"].append(param_name)
    return pd.DataFrame.from_dict(data_dic)


def get_errors(
    ds,
    perturb_names,
    store_path=None,
):
    """
    Given a dataset with multiple ensembles, where ensemble 0 is the unperturbed trajectory and
    all the others consist of trajectories with perturbed members,
    calculate the difference in each timestep. Then reduce it over time to get the mean squared error and mean error.
    This function assumes that every trajectory within one ensemble is a different trajectory.

    Parameters
    ----------
    ds : xarray.Dataset
        Simulation created with the mode "limited_time_ensembles".
    perturb_names : list-like of strings
        Names of the parameters perturbed in each ensemble.
    store_path : path (optional)
        Path and name where to store the final dataframe.

    Returns
    -------
    pandas.DataFrame with columns "Mean Squared Error", "Mean Error", "Output Parameter", "Input Parameter",
    "Predicted Error", and "Predicted Squared Error".
    """
    if "Output_Parameter_ID" in ds:
        coord_name = "Output_Parameter_ID"
        out_params = []
        if coord_name not in ds.dims:
            out_params.append(param_id_map[ds[coord_name].values.item()])
        else:
            for param_id in ds[coord_name].values:
                out_params.append(param_id_map[param_id])
    else:
        coord_name = "Output Parameter"
        out_params = list(ds[coord_name].values)
    perturb_names_avail = []
    for p in np.unique(perturb_names):
        if p in ds:
            perturb_names_avail.append(p)

    out_ids = ds[coord_name].values
    if np.shape(out_ids) == ():
        out_ids = [out_ids]
    min_time_idx = np.argwhere(
        (~np.isnan(ds.isel({"ensemble": 1, "trajectory": 0})[out_params[0]])).values
    ).min()
    ds = ds.isel({"time": np.arange(min_time_idx, len(ds["time"]))})
    n_nans = _get_n_nans(ds, perturb_names, coord_name, out_ids, out_params)

    df = _get_error_df(
        ds,
        perturb_names,
        coord_name,
        out_ids,
        out_params,
        n_nans,
        perturb_names_avail,
    )

    if store_path is not None:
        ds = df.to_xarray()
        comp = {"zlib": True, "complevel": 9}
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(
            path=store_path,
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )
    return df


def main(args):
    """

    Parameters
    ----------
    args :
        Arguments parsed via argparse.
    """

    if args.path[-1] == "/":
        get_errors(
            ds=load_dataset(args.path),
            perturb_names=args.perturb_names,
            store_path=args.store_path + args.path.split("/")[-1],
        )
    else:
        files = [f for f in os.listdir(args.path) if os.path.isfile(args.path + f)]
        files.sort()
        for f in tqdm(files):
            get_errors(
                ds=load_dataset(args.path + f),
                perturb_names=args.perturb_names,
                store_path=args.store_path + f,
            )


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Load ensemble NetCDF-files and calculate the errors
            and reduced predicted errors and store them to disk. This is used to
            prepare data to find any correlations and for plotting.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--path",
        default="data/perturbed_ensembles/",
        help=textwrap.dedent(
            """\
            Path to perturbed ensembles. Should point either to a file or to a directory to process multiple
            files.  
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        default="data/stats/",
        help=textwrap.dedent(
            """\
            Path to store the resulting statistics.
            """
        ),
    )

    parser.add_argument(
        "--filename",
        type=str,
        required=True,
        help=textwrap.dedent(
            """\
            Name of the file where trajectories and subsequent ensembles are
            stored.
            """
        ),
    )
    parser.add_argument(
        "--perturb_names",
        type=str,
        nargs="+",
        default=[],
        help=textwrap.dedent(
            """\
            Name of model parameters perturbed in each ensemble. 
            """
        ),
    )
    main(parser.parse_args())
