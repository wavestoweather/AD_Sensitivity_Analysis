import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

ColourReset = "\033[0m"
Warning = "\033[93m"
Success = "\033[92m"
Status = "\033[94m"
Error = "\033[91m"


def load_dataset(input):
    """
    Load a dataset.

    Parameters
    ----------
    input : string
        Path to the NetCDF-file.

    Returns
    -------
    xr.Dataset with trajectory data.
    """
    return xr.open_dataset(input, decode_times=False, engine="netcdf4")


def plot_distributions(ds, path, verbose):
    """
    Plot a histogram for every variable (column). The arrays are flattened.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    path : string
        Path to folder where plots will be stored as png.
    verbose : bool
        If true: print which variable is being plotted.
    """
    print("~*~*~*~Plotting distributions~*~*~*~")
    for col in ds:
        if verbose:
            print(f"Plotting {col}")
        ds[col].plot.hist(size=10, bins=100)
        plt.tight_layout()
        plt.savefig(f"{path}hist_{col}.png", dpi=300)
        plt.close()
        ds[col].plot.hist(size=10, bins=100, log=True)
        plt.tight_layout()
        plt.savefig(f"{path}hist_{col}_log.png", dpi=300)
        plt.close()
    print("done\n")


def plot_columns(ds, path, verbose):
    """
    Plot the columns, where ensemble and Output_Parameter_ID is squeezed
    (if present and possible).
    Otherwise an additional plot for each ensemble and ID

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    path : string
        Path where plots shall be stored.
    verbose : bool
        If true: print which variable is being plotted.
    """
    print("~*~*~*~Plotting columns~*~*~*~")
    ds_tmp = ds.copy()
    if "ensemble" in ds_tmp:
        if len(ds_tmp["ensemble"]) == 1:
            ds_tmp = ds_tmp.squeeze(dim="ensemble", drop=True)
    if "Output_Parameter_ID" in ds_tmp:
        if len(ds_tmp["Output_Parameter_ID"]) == 1:
            ds_tmp = ds_tmp.squeeze(dim="Output_Parameter_ID", drop=True)
    for col in ds_tmp:
        if len(ds_tmp[col]) > 0:
            if verbose:
                print(f"Plotting {col}")
            if "ensemble" in ds_tmp[col].dims:
                if "Output_Parameter_ID" in ds_tmp[col].dims:
                    for e in ds_tmp["ensemble"].values:
                        for i in ds_tmp["Output_Parameter_ID"].values:
                            ds_tmp[col].sel(
                                {"ensemble": e, "Output_Parameter_ID": i}
                            ).plot(size=10)
                            plt.tight_layout()
                            plt.savefig(f"{path}{col}_ens{e}_id{i}.png", dpi=300)
                            plt.close()
                else:
                    for e in ds_tmp["ensemble"].values:
                        ds_tmp[col].sel({"ensemble": e}).plot(size=10)
                        plt.tight_layout()
                        plt.savefig(f"{path}{col}_ens{e}.png", dpi=300)
                        plt.close()
            elif "Output_Parameter_ID" in ds_tmp[col].dims:
                for i in ds_tmp["Output_Parameter_ID"].values:
                    ds_tmp[col].sel({"Output_Parameter_ID": i}).plot(size=10)
                    plt.tight_layout()
                    plt.savefig(f"{path}{col}_id{i}.png", dpi=300)
                    plt.close()
            else:
                ds_tmp[col].plot(size=10)
                plt.tight_layout()
                plt.savefig(f"{path}{col}.png", dpi=300)
                plt.close()
    print("done\n")


def test_sizes(ds, threshold, verbosity):
    """
    Test if any hydrometeor sizes are rather unphysical. We assume standard
    ICON sizes for those.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    threshold : float
        The amount of time steps where we allow unphysical hydrometeors before
        we count it as an error.
    verbosity : int
        0: print only results.
        <3: print each error in the test
        >=3: print the amount of non-physical sizes for each hydrometeor even
        if those are within the given threshold.

    Returns
    -------
    Number of hydrometeors that have at least threshold many unphysical time
    steps.
    """
    print("~*~*~*~Testing hydrometeor sizes~*~*~*~")
    # All possible pairs of mass and numbers
    pairs = [
        ("QC", "NCCLOUD"),
        ("QR", "NCRAIN"),
        ("QI", "NCICE"),
        ("QS", "NCSNOW"),
        ("QG", "NCGRAUPEL"),
        ("QR_IN", "NR_IN"),
        ("QI_IN", "NI_IN"),
        ("QS_IN", "NS_IN"),
        ("QG_IN", "NG_IN"),
        ("QR_OUT", "NR_OUT"),
        ("QI_OUT", "NI_OUT"),
        ("QS_OUT", "NS_OUT"),
        ("QG_OUT", "NG_OUT"),
        ("QC", "QNC"),
        ("QR", "QNR"),
        ("QI", "QNI"),
        ("QS", "QNS"),
        ("QG", "QNG"),
        ("QH", "QNH"),
        ("QR_in", "QNR_in"),
        ("QI_in", "QNI_in"),
        ("QS_in", "QNS_in"),
        ("QG_in", "QNG_in"),
        ("QH_in", "QNH_in"),
    ]
    model_constants = {
        "cloud": {"min": 4.2e-15, "max": 2.6e-10},
        "rain": {"min": 2.6e-10, "max": 3.0e-6},
        "snow": {"min": 1.0e-10, "max": 2.0e-5},
        "ice": {"min": 1.0e-12, "max": 1.0e-5},
        "graupel": {"min": 4.19e-9, "max": 5.3e-4},
        "hail": {"min": 2.6e-9, "max": 5.4e-4},
    }
    err = 0
    n_tests = 0
    n_ens = 1
    if "ensemble" in ds:
        n_ens = len(ds["ensemble"])
    n_steps = len(ds["trajectory"]) * len(ds["time"]) * n_ens

    def check_valid(ty, q, n):
        avg = ds[q] / ds[n]
        if "_in" in q:
            avg = np.abs(avg)
        masked_greater = np.ma.masked_greater(
            avg, model_constants[ty]["max"], copy=False
        )
        masked_less = np.ma.masked_less(avg, model_constants[ty]["min"], copy=False)
        error_larger = np.sum(masked_greater.mask)
        if error_larger > 0 and verbosity >= 3:
            print(
                f"{Warning}{ty}: too big: {error_larger} or {error_larger*100/n_steps:2.2f}% of the time{ColourReset}"
            )
        error_smaller = np.sum(masked_less.mask)
        if error_smaller > 0 and verbosity >= 3:
            print(
                f"{Warning}{ty}: too small {error_smaller} or {error_smaller*100/n_steps:2.2f}% of the time{ColourReset}"
            )
        if (error_larger + error_smaller) / n_steps * 100 >= threshold:
            if verbosity > 0:
                print(
                    f"{Warning}{ty}: too big: {error_larger} or {error_larger*100/n_steps:2.2f}% of the time{ColourReset}"
                )
                print(
                    f"{Warning}{ty}: too small {error_smaller} or {error_smaller*100/n_steps:2.2f}% of the time{ColourReset}"
                )
            return 1
        return 0

    for q, n in pairs:
        if q in ds and n in ds:
            if "QC" in q:
                ty = "cloud"
            elif "QR" in q:
                ty = "rain"
            elif "QI" in q:
                ty = "ice"
            elif "QS" in q:
                ty = "snow"
            elif "QG" in q:
                ty = "graupel"
            else:
                ty = "hail"
            err += check_valid(ty, q, n)
            n_tests += 1

    if err > 0:
        print(
            f"{Error}Failed: Unphysical sizes for at least {threshold:2.2f}% time steps detected for {err} of {n_tests} hydrometeors.{ColourReset}\n"
        )
    else:
        print(
            f"{Success}No more than {threshold:2.2f}% of time steps with unphysical sizes detected for all {n_tests} hydrometeors.{ColourReset}\n"
        )
    return err


def test_nan_col(ds, verbose):
    """
    Test if any NaNs occur in any column (variable) of the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    verbose : bool
        If true: print each error in the NaN test for columns (variables).
        Otherwise print only the results of the test.

    Returns
    -------
    Number of columns with NaNs.
    """
    print("~*~*~*~Testing NaNs in data variables~*~*~*~")
    err = 0
    n_vars = 0
    for col in ds:
        n_vars += 1
        n_nans = np.sum(np.isnan(ds[col]))
        if n_nans > 0:
            shape = np.shape(ds[col])
            n = 1
            for s in shape:
                n *= s
            if verbose:
                print(
                    f"{Warning}There are {n_nans.item()} ({(n_nans/n).item():2.2f}%) NaNs in variable {col}{ColourReset}"
                )
            err += 1
    if err > 0:
        print(
            f"{Error}Failed: NaNs detected for {err} of {n_vars} variables.{ColourReset}\n"
        )
    else:
        print(f"{Success}No NaNs for all {n_vars} variables detected.{ColourReset}\n")
    return err


def test_nan_dim(ds, verbose):
    """
    Test if any NaNs occur in any dimension of the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    verbose : bool
        If true: print all errors for each dimension. Otherwise print only
        the results.

    Returns
    -------
    Number of dimensions with NaNs.
    """
    print("~*~*~*~Testing NaNs in dimensions~*~*~*~")
    err = 0
    n_dims = 0
    for d in ds.dims:
        n_dims += 1
        n_nans = np.sum(np.isnan(ds[d]))
        if n_nans > 0:
            shape = np.shape(ds[d])
            n = 1
            for s in shape:
                n *= s
            if verbose:
                print(
                    f"{Warning}There are {n_nans.item()} ({(n_nans/n).item():2.2f}%) NaNs in dimension {d}{ColourReset}"
                )
            err += 1
    if err > 0:
        print(
            f"{Error}Failed: NaNs detected for {err} of {n_dims} dimensions.{ColourReset}\n"
        )
    else:
        print(f"{Success}No NaNs for all {n_dims} dimensions detected.{ColourReset}\n")
    return err


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Test a given NetCDF-file for possible errors, such as hydrometeors that
        are too large or too small, get the distribution of values, check for NaNs
        in the coordinates (dimensions). Works only for trajectories. Trajectories
        are made via interpolation of gridpoints, hence some unphysical
        hydrometeors can occur.
        """
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="""
        Path to the NetCDF-file. If the given path is a folder of files, all
        files will be used sequentially.
        """,
    )
    parser.add_argument(
        "--output_path",
        type=str,
        default="../pics/tests/",
        help="""
        Path for plotting distributions and different variables if 'plot' is true.
        """,
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="""
        Plot distributions and different variables to 'output_folder'.
        """,
    )
    parser.add_argument(
        "--allow_failure",
        action="store_true",
        help="""
        The program shall not return an errorcode if the dataset is suspicious.
        An errorcode is helpful for gitlab CI.
        """,
    )
    parser.add_argument(
        "--allowed_perc",
        type=float,
        default=13.0,
        help=""""
        The amount of hydrometeors that are allowed to be outside of normal,
        physical values. Interpolation can lead to that. Usually either cloud droplets
        or ice crystals have values of 8 or 12 percent. Other hydrometeors may be
        as low as 2 percent.
        """,
    )
    parser.add_argument(
        "--test_nan_dims",
        action="store_true",
        help="""
        Test NaNs in the dataset for all dimensions. This should not happen and
        usually indicates a problem with the data.
        """,
    )
    parser.add_argument(
        "--test_nan_vars",
        action="store_true",
        help="""
        Test NaNs in the dataset for all variables (not dimensions). This can
        happen, when a trajectory stops early, i.e., because it "crashed" into a
        mountain.
        """,
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help="""
        Control the amount of outputs. Verbosity 0 prints only results of
        tests. Verbosity 1 prints in addition each error in the tests for non-physical
        hydrometeor sizes and which variables are being plotted (if any).
        Verbosity 2 prints in addition each error in the NaN test for
        dimensions. Verbosity 3 prints in addition the amount of non-physical
        hydrometeor sizes even if those are below the given threshold. Verbosity 4
        prints in addition each error in the NaN test for variables (columns).
        """,
    )
    args = parser.parse_args()

    n_files = 0

    def parse_ds(ds, output_path):
        errors = 0
        if args.plot:
            plot_distributions(ds, output_path, (args.verbosity > 0))
            plot_columns(ds, output_path, (args.verbosity > 0))
        if args.test_nan_dims:
            errors += test_nan_dim(ds, (args.verbosity > 1))
        if args.test_nan_vars:
            errors += test_nan_col(ds, (args.verbosity > 3))
        errors += test_sizes(ds, args.allowed_perc, args.verbosity)
        return errors

    err = 0
    if args.input[-1] == "/":
        n = len(os.listdir(args.input))
        for i, f in enumerate(os.listdir(args.input)):
            n_files += 1
            print(f"{Status}~*~*~*~Parsing {f} {i+1}/{n}~*~*~*~{ColourReset}")
            ds = load_dataset(args.input + f)
            err += parse_ds(ds, args.output_path + f.split(".")[0] + "_")
    else:
        n_files += 1
        ds = load_dataset(args.input)
        output_path = args.output_path
        if output_path[-1] != "/":
            output_path += "/"
        err += parse_ds(ds, output_path)
    if err == 0:
        print(f"\n{Success}No errors occured for {n_files} files.{ColourReset}")
    if not args.allow_failure and err > 0:
        raise Exception(
            f"{Error}Failed: {err} errors occured during testing for {n_files} files. Check the output!{ColourReset}"
        )
    elif err > 0:
        print(
            f"{Error}Failed: {err} errors occured during testing for {n_files} files. Check the output!{ColourReset}"
        )
