import numpy as np
import os
from tqdm.auto import tqdm
import xarray as xr


def transform_coord(ds, ds_orig=None, zero_condition="ascent_start"):
    """
    Transform latitude and longitude such that
    zero is at the given condition.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset that contains 'lon', 'lat', and 'time_after_ascent'.
    ds_orig : xarray.Dataset
        You may provide the original dataset without a sensitivity analysis to
        get the relative coordinates. This is necessary if the sensitivity analysis
        used a warm-up time such that the start of the ascent is not in the output.
    zero_condition : string
        Possible condition is 'ascent_start' where the coordinates are realative
        to the position to the start of the ascent of the trajectories.

    Returns
    -------
    xarray.Dataset with transformed coordinates
    """
    if zero_condition == "ascent_start" and ds_orig is None:
        lon_start = ds["lon"].where(ds["time_after_ascent"] == 0)
        lat_start = ds["lat"].where(ds["time_after_ascent"] == 0)
    elif zero_condition == "ascent_start":
        lon_start = ds_orig["lon"].where(ds_orig["time_after_asc_start"] == 0)
        lat_start = ds_orig["lat"].where(ds_orig["time_after_asc_start"] == 0)
    shape = np.shape(ds["lon"])[0:-1] + (1,)
    lon_shift = np.zeros(shape)
    lat_shift = np.zeros(shape)
    if "ensemble" in lon_start:
        for ens_idx, ens in enumerate(lon_start["ensemble"]):
            for traj_idx, traj in enumerate(lon_start["trajectory"]):
                lon_shift[ens_idx, traj_idx] = np.nansum(
                    lon_start.sel({"ensemble": ens, "trajectory": traj})
                )
                lat_shift[ens_idx, traj_idx] = np.nansum(
                    lat_start.sel({"ensemble": ens, "trajectory": traj})
                )
    else:
        for traj_idx, traj in enumerate(lon_start["trajectory"]):
            lon_shift[0, traj_idx] = np.nansum(lon_start.sel({"trajectory": traj}))
            lat_shift[0, traj_idx] = np.nansum(lat_start.sel({"trajectory": traj}))
    ds["relative_lon"] = ds["lon"] - lon_shift
    ds["relative_lat"] = ds["lat"] - lat_shift
    ds["relative_lon"].attrs = {
        "long_name": "longitude relative to start of ascent",
        "standard_name": "relative_longitude",
        "units": "degrees",
    }
    ds["relative_lat"].attrs = {
        "long_name": "latitude relative to start of ascent",
        "standard_name": "relative_latitude",
        "units": "degrees",
    }
    return ds


def correct_coordinates(ds):
    """
    An old version of the sensitivity analysis had wrong time after ascent values. This
    function corrects it.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset that contains 'time_after_ascent'.

    Returns
    -------
    xarray.Dataset with corrected coordinates
    """
    ds["time_after_ascent"] = (ds["time_after_ascent"] - 29) * 30
    return ds


def main(args):
    files = [
        f for f in os.listdir(args.file_path) if os.path.isfile(args.file_path + f)
    ]
    comp = dict(zlib=True, complevel=9)
    for f in tqdm(files) if args.verbose else files:
        ds = xr.open_dataset(
            args.file_path + f,
            decode_times=False,
            engine="netcdf4",
        )
        if args.correct:
            ds = correct_coordinates(ds)
        ds_orig = None
        if args.file_path_orig is not None:
            ds_orig = xr.open_dataset(
                args.file_path_orig + f,
                decode_times=False,
                engine="netcdf4",
            )
        ds = transform_coord(ds, ds_orig)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(
            path=args.store_path + f,
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w",
        )


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Load data from a sensitivity simulation and transform the coordinate values.
            Currently supports only longitude and latitude transformation.
            Can be used to correct the time after ascent from old sensitivity analysis.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file_path",
        required=True,
        help=textwrap.dedent(
            """\
            Path to a folder with many files from a sensitivity analysis simulation.
            """
        ),
    )
    parser.add_argument(
        "--file_path_orig",
        required=True,
        help=textwrap.dedent(
            """\
            Path to a folder with many files with trajectories from an ICON simulation.
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        default=None,
        help=textwrap.dedent(
            """\
            Path to a folder to store the new files.
            """
        ),
    )
    parser.add_argument(
        "--correct",
        action="store_true",
        help=textwrap.dedent(
            """\
            Correct the time after ascent from old sensitivity analysis.
            """
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=textwrap.dedent(
            """\
            More output, i.e., a progressbar.
            """
        ),
    )
    args = parser.parse_args()
    main(args)
