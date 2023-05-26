"""Test results from a sensitivity simulation.

"""
import os

import numpy as np

from ad_sensitivity_analysis.data_handler.loader import load_dataset
from ad_sensitivity_analysis.exception.test_exception import TestException
from ad_sensitivity_analysis.tests.end_to_end.test_phases import test_phases
from ad_sensitivity_analysis.interface.physics_t import PhysicsT
from ad_sensitivity_analysis.plot.latexify import param_id_map
from ad_sensitivity_analysis.plot.simple_sim_plot import (
    plot_columns,
    plot_distributions,
)
from ad_sensitivity_analysis.tests.color_scheme import (
    COLOR_RESET,
    WARNING,
    SUCCESS,
    STATUS,
    ERROR,
)
from ad_sensitivity_analysis.tests.unit.test_physics import (
    test_ccn_act_akm,
    test_graupel_melting,
)


def test_amounts(ds, verbosity):
    """
    Test if the amount of water vapor, cloud droplets, rain droplets, ice, snow, graupel, and hail is more than
    reasonable.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    verbosity : int
        0: print only results.
        <3: print each trajectory error in the test
        >3: print each hydrometeor error

    Returns
    -------
    Number of trajectories and time steps with unphysical amounts.
    """
    print("~*~*~*~Testing hydrometeor amounts~*~*~*~")
    threshold = 1e-1
    err = 0
    err_traj = 0
    n_traj = len(ds["trajectory"])
    hydrometeors = [
        "QV",
        "QC",
        "QR",
        "QS",
        "QI",
        "QG",
        "QH",
        "QR_OUT",
        "QS_OUT",
        "QI_OUT",
        "QG_OUT",
        "QH_OUT",
    ]
    for i in range(n_traj):
        ds_t = ds.isel({"trajectory": i})
        err_tmp = 0
        for param in hydrometeors:
            err_tmp += np.sum(ds_t[param] >= threshold).item()
            if verbosity > 3:
                n = np.sum(~np.isnan(ds_t[param]))
                print(
                    f"{WARNING}There are {err_tmp} ({(err_tmp/n.item()):2.2f}%) steps with unphysical amounts "
                    f"in trajectory {i} and param {param}{COLOR_RESET}"
                )
        if err_tmp > 0:
            err_traj += 1
            err += err_tmp
            if verbosity > 0:
                n = np.sum(~np.isnan(ds_t[param])) * len(hydrometeors)
                print(
                    f"{WARNING}There are {err_tmp} ({(err_tmp/n.item()):2.2f}%) steps with unphysical amounts "
                    f"in trajectory {i}{COLOR_RESET}"
                )
    if err > 0:
        print(
            f"{ERROR}Failed: Unphysical amounts detected for {err_traj} of {n_traj} trajectories.{COLOR_RESET}\n"
        )
    else:
        print(
            f"{SUCCESS}No unphysical amounts for all {n_traj} trajectories detected.{COLOR_RESET}\n"
        )
    return err, err_traj


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

    def check_valid(hydro_type, hydro_name, n):
        avg = ds[hydro_name] / ds[n]
        if "_in" in hydro_name:
            avg = np.abs(avg)
        # pylint: disable=no-member
        masked_greater = np.ma.masked_greater(
            avg, model_constants[hydro_type]["max"], copy=False
        )
        masked_less = np.ma.masked_less(
            avg, model_constants[hydro_type]["min"], copy=False
        )
        error_larger = np.sum(masked_greater.mask)
        if error_larger > 0 and verbosity >= 3:
            print(
                f"{WARNING}{hydro_type}: too big: {error_larger} "
                f"or {error_larger*100/n_steps:2.2f}% of the time{COLOR_RESET}"
            )
        error_smaller = np.sum(masked_less.mask)
        if error_smaller > 0 and verbosity >= 3:
            print(
                f"{WARNING}{hydro_type}: too small {error_smaller} "
                f"or {error_smaller*100/n_steps:2.2f}% of the time{COLOR_RESET}"
            )
        if (error_larger + error_smaller) / n_steps * 100 >= threshold:
            if verbosity > 0:
                print(
                    f"{WARNING}{hydro_type}: too big: {error_larger} "
                    f"or {error_larger*100/n_steps:2.2f}% of the time{COLOR_RESET}"
                )
                print(
                    f"{WARNING}{hydro_type}: too small {error_smaller} "
                    f"or {error_smaller*100/n_steps:2.2f}% of the time{COLOR_RESET}"
                )
            return 1
        return 0

    for hydro_name, n in pairs:
        if hydro_name in ds and n in ds:
            if "QC" in hydro_name:
                hydro_type = "cloud"
            elif "QR" in hydro_name:
                hydro_type = "rain"
            elif "QI" in hydro_name:
                hydro_type = "ice"
            elif "QS" in hydro_name:
                hydro_type = "snow"
            elif "QG" in hydro_name:
                hydro_type = "graupel"
            else:
                hydro_type = "hail"
            err += check_valid(hydro_type, hydro_name, n)
            n_tests += 1

    if err > 0:
        print(
            f"{ERROR}Failed: Unphysical sizes for at least {threshold:2.2f}% time steps detected "
            f"for {err} of {n_tests} hydrometeors.{COLOR_RESET}\n"
        )
    else:
        print(
            f"{SUCCESS}No more than {threshold:2.2f}% of time steps with unphysical sizes detected "
            f"for all {n_tests} hydrometeors.{COLOR_RESET}\n"
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
        if "efficiency" in col:
            continue
        n_vars += 1
        n_nans = np.sum(np.isnan(ds[col]))
        if n_nans > 0:
            shape = np.shape(ds[col])
            n = 1
            for s_number in shape:
                n *= s_number
            if verbose:
                print(
                    f"{WARNING}There are {n_nans.item()} "
                    f"({(n_nans*100/n).item():2.2f}%) NaNs in variable {col}{COLOR_RESET}"
                )
            err += 1
    if err > 0:
        print(
            f"{ERROR}Failed: NaNs detected for {err} of {n_vars} variables.{COLOR_RESET}\n"
        )
    else:
        print(f"{SUCCESS}No NaNs for all {n_vars} variables detected.{COLOR_RESET}\n")
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
    for dim_name in ds.dims:
        n_dims += 1
        n_nans = np.sum(np.isnan(ds[dim_name]))
        if n_nans > 0:
            shape = np.shape(ds[dim_name])
            n = 1
            for s_number in shape:
                n *= s_number
            if verbose:
                print(
                    f"{WARNING}There are {n_nans.item()} ({(n_nans/n).item():2.2f}%) NaNs "
                    f"in dimension {dim_name}{COLOR_RESET}"
                )
            err += 1
    if err > 0:
        print(
            f"{ERROR}Failed: NaNs detected for {err} of {n_dims} dimensions.{COLOR_RESET}\n"
        )
    else:
        print(f"{SUCCESS}No NaNs for all {n_dims} dimensions detected.{COLOR_RESET}\n")
    return err


def test_saturation(ds, recalc, verbose):
    """
    Test if there are time steps which are over saturated for each trajectory.
    Calculates the saturation if necessary.
    Careful! Some datafiles may use the unit "percentage" when in fact, it is not percentage for "S".
    You may have to recalculate the saturation in those cases.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    recalc : bool
        Recalculate the saturation.
    verbose : bool
        If true: print all over saturations for each trajectory. Otherwise print only
        the results.

    Returns
    -------
    Number of trajectories and time steps with over saturation.
    """
    print("~*~*~*~Testing over saturation in each trajectory~*~*~*~")
    err = 0
    err_traj = 0

    def sat_p_water(T):
        p_sat_low_temp = 610.78
        p_sat_const_a = 17.2693882
        p_sat_const_b = 35.86
        t_sat_low_temp = 273.15
        return p_sat_low_temp * np.exp(
            p_sat_const_a * (T - t_sat_low_temp) / (T - p_sat_const_b)
        )

    def calc_saturation(p, qv, T):
        r_universal = 8.3144598
        m_w = 0.018015265
        m_a = 0.02896546
        r_a = r_universal / m_a
        r_v = r_universal / m_w
        epsilon = r_a / r_v
        return (p * qv) / ((epsilon + qv) * sat_p_water(T))

    def is_oversat(p, qv, T):
        return calc_saturation(p, qv, T) > 1

    n_traj = len(ds["trajectory"])
    if "S" not in ds:
        recalc = True
    si_unit = 1
    sat_unit = 1
    if "S" in ds and not recalc:
        if ds["S"].attrs["units"] == "percentage":
            sat_unit = 100
    if ds["pressure"].attrs["units"] == "hPa":
        si_unit = 100
    for i in range(n_traj):
        ds_t = ds.isel({"trajectory": i})
        if recalc:
            err_tmp = np.sum(
                is_oversat(si_unit * ds_t["pressure"], ds_t["QV"], ds_t["T"])
            ).item()
        else:
            err_tmp = np.sum((np.asarray(ds_t["S"]) > sat_unit))
        if err_tmp > 0:
            err_traj += 1
            if verbose:
                n = len(ds_t["time"]) * len(ds_t["ensemble"])
                print(
                    f"{WARNING}There are {err_tmp} ({(err_tmp/n):2.2f}%) over saturated steps "
                    f"in trajectory {i}{COLOR_RESET}"
                )
        err += err_tmp
    if err > 0:
        print(
            f"{ERROR}Failed: Over saturation detected for {err_traj} of {n_traj} trajectories.{COLOR_RESET}\n"
        )
    else:
        print(
            f"{SUCCESS}No over saturation for all {n_traj} trajectories detected.{COLOR_RESET}\n"
        )
    return err, err_traj


def check_sensitivity_size(
    errors, ds_tmp, col, verbosity, threshold, n_steps, out_name, ntrajs
):
    """
    Check if sensitivities are too large for a given model parameter and threshold.

    Parameters
    ----------
    errors
    ds_tmp
    col
    verbosity
    threshold
    n_steps
    out_name
    ntrajs

    Returns
    -------

    """
    if verbosity >= 5:
        print(
            f"{STATUS}d{out_name}/{col}: {np.nanmin(ds_tmp[col])}, {np.nanmax(ds_tmp[col])}, "
            f"{np.nanmean(ds_tmp[col])}, {np.nanmedian(ds_tmp[col])}{COLOR_RESET}"
        )
    if verbosity >= 4:
        already_error = False
        for traj in np.arange(ntrajs):
            # pylint: disable=no-member
            masked_greater = np.ma.masked_greater(
                np.abs(ds_tmp.isel({"trajectory": traj})[col]), threshold, copy=False
            )
            errors["error_local"] = np.sum(masked_greater.mask)
            errors["err_total"] += errors["error_local"]
            if errors["error_local"] > 0:
                print(
                    f"{WARNING}d{out_name}/{col}, traj {traj}: too big sensitivities: "
                    f"{errors['error_local']} or {errors['error_local']*100/n_steps:2.2f}% "
                    f"of the time{COLOR_RESET}"
                )
            if errors["error_local"] > 0 and not already_error:
                errors["err"] += 1
                already_error = True
    else:
        masked_greater = np.ma.masked_greater(
            np.abs(ds_tmp[col]), threshold, copy=False
        )
        # pylint: disable=no-member
        errors["error_local"] = np.sum(masked_greater.mask)
        errors["err_total"] += errors["error_local"]
        if errors["error_local"] > 0 and verbosity >= 3:
            print(
                f"{WARNING}d{out_name}/{col}: too big sensitivities: "
                f"{errors['error_local']} or {errors['error_local']*100/n_steps:2.2f}% of the time{COLOR_RESET}"
            )
        if errors["error_local"] > 0:
            errors["err"] += 1
    return errors


def test_sensitivities(ds, threshold, verbosity):
    """
    Test if there are time steps and model parameters with too large sensitivities. Uses the absolute value
    to check for size.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data from a sensitivity simulation.
    threshold : float
        The threshold at which a sensitivity is considered too large.
    verbosity : verbosity
        0: print only results.
        >=3: print errors for each model state variable w.r.t. each model parameter.
        >=4: print the trajectory at which the large sensitivities occurred.
        5: print minimum, maximum, mean, and median sensitivities.

    Returns
    -------
    Number of model parameters for each model state variables with large sensitivities.
    """
    if verbosity >= 4:
        print("~*~*~*~Testing sensitivities in each trajectory~*~*~*~")
    else:
        print("~*~*~*~Testing sensitivities~*~*~*~")
    errors = {
        "err": 0,
        "err_total": 0,
        "error_local": 0,
    }
    n_steps = len(ds["trajectory"]) * len(ds["time"]) * len(ds["ensemble"])
    n_params = 0
    out_coord = (
        "Output_Parameter_ID" if "Output_Parameter_ID" in ds else "Output Parameter"
    )
    for out_p in ds[out_coord]:
        out_name = (
            out_p.item()
            if out_coord == "Output Parameter"
            else param_id_map[out_p.item()]
        )
        if out_name[0] == "N":
            continue
        ds_tmp = ds.sel({out_coord: out_p})
        for col in ds:
            if col[0] != "d" or col == "deposition":
                continue
            n_params += 1
            errors = check_sensitivity_size(
                errors=errors,
                ds_tmp=ds_tmp,
                col=col,
                verbosity=verbosity,
                threshold=threshold,
                n_steps=n_steps,
                out_name=out_name,
                ntrajs=len(ds["trajectory"]),
            )
    if errors["err"] > 0:
        print(
            f"{ERROR}Failed: Too large sensitivities detected for {errors['err']} model parameters. "
            f"Error rate of {n_steps*errors['err_total'] /n_params*100:2.2f}.{COLOR_RESET}\n"
        )
    else:
        print(
            f"{SUCCESS}No large sensitivities for all model parameters detected.{COLOR_RESET}\n"
        )
    return errors["err"]


def test_physics(arguments):
    """

    Returns
    -------

    """
    err = 0
    # Pylint assumes that 'n_processes' is a constant.
    # pylint: disable=invalid-name
    n_processes = 0
    physics = PhysicsT(lib_path=arguments.test_physics)
    physics.setup_model_constants(arguments.timestep, 0.1)
    print(f"{STATUS}~*~*~*~Testing graupel melting~*~*~*~{COLOR_RESET}")
    TMP = test_graupel_melting(
        physics=physics,
        T=285,
        qg=2,
        Ng=8,
        verbose=(arguments.verbosity > 4),
    )
    if TMP > 0:
        n_processes += 1
    err += TMP

    print(f"{STATUS}~*~*~*~Testing CCN activation (akm)~*~*~*~{COLOR_RESET}")
    TMP = test_ccn_act_akm(
        physics=physics,
        w=3,
        T=270,
        qv=1e2,
        qc=-6,
        verbose=(arguments.verbosity > 4),
    )
    if TMP > 0:
        n_processes += 1
    err += TMP

    if err == 0:
        print(f"\n{SUCCESS}No errors occured for all physics tests.{COLOR_RESET}")
    if not arguments.allow_failure and err > 0:
        raise TestException(
            f"{ERROR}Failed: {err} errors occured during testing for {n_processes} processes.{COLOR_RESET}"
        )
    if err > 0:
        print(
            f"{ERROR}Failed: {err} errors occured during testing for {n_processes} processes.{COLOR_RESET}"
        )


def test_simulation_output(ds, output_path, arguments):
    """

    Parameters
    ----------
    ds
    output_path
    arguments

    Returns
    -------

    """
    errors = 0
    if arguments.plot:
        plot_distributions(ds, output_path, (arguments.verbosity > 0))
        plot_columns(ds, output_path, (arguments.verbosity > 0))
    if arguments.test_nan_dims:
        errors += test_nan_dim(ds, (arguments.verbosity > 1))
    if arguments.test_nan_vars:
        errors += test_nan_col(ds, (arguments.verbosity > 3))
    errors += test_sizes(ds, arguments.allowed_perc, arguments.verbosity)
    if arguments.test_saturation:
        _, err_traj = test_saturation(
            ds, arguments.calc_saturation, arguments.verbosity
        )
        errors += err_traj
    if arguments.test_phases:
        errors += test_phases(ds, arguments.calc_phases)
    if arguments.test_sensitivities:
        errors += test_sensitivities(ds, 10.0, arguments.verbosity)
    if arguments.test_amounts:
        _, err_traj = test_amounts(ds, arguments.verbosity)
    return errors


def main(arguments):
    """

    Parameters
    ----------
    arguments

    Returns
    -------

    """
    n_files = 0
    if arguments.input is not None:
        err = 0
        if arguments.input[-1] == "/":
            files = [
                f
                for f in os.listdir(arguments.input)
                if os.path.isfile(arguments.input + f)
            ]
            n = len(files)
            for i, f in enumerate(files):
                n_files += 1
                print(f"{STATUS}~*~*~*~Parsing {f} {i+1}/{n}~*~*~*~{COLOR_RESET}")
                ds = load_dataset(arguments.input + f)
                err += test_simulation_output(
                    ds, arguments.output_path + f.split(".")[0] + "_", arguments
                )
        else:
            n_files += 1
            ds = load_dataset(arguments.input)
            output_path = arguments.output_path
            if output_path[-1] != "/":
                output_path += "/"
            err += test_simulation_output(ds, output_path, arguments)

        if err == 0:
            print(f"\n{SUCCESS}No errors occured for {n_files} file(s).{COLOR_RESET}")
        if not arguments.allow_failure and err > 0:
            raise TestException(
                f"{ERROR}Failed: {err} errors occured during testing for {n_files} file(s). "
                f"Check the output!{COLOR_RESET}"
            )
        if err > 0:
            print(
                f"{ERROR}Failed: {err} errors occured during testing for {n_files} file(s). "
                f"Check the output!{COLOR_RESET}"
            )

    if arguments.test_physics != "no":
        test_physics(arguments)


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Test a given NetCDF-file for possible errors, such as hydrometeors that
            are too large or too small, get the distribution of values, check for NaNs
            in the coordinates (dimensions). Works only for trajectories. Trajectories
            are made via interpolation of gridpoints, hence some unphysical
            hydrometeors can occur.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help=textwrap.dedent(
            """\
            Path to the NetCDF-file. If the given path is a folder of files, all
            files will be used sequentially.
            """
        ),
    )
    parser.add_argument(
        "--output_path",
        type=str,
        default="../pics/tests/",
        help=textwrap.dedent(
            """\
            Path for plotting distributions and different variables if 'plot' is true.
            """
        ),
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help=textwrap.dedent(
            """\
            Plot distributions and different variables to 'output_folder'.
            """
        ),
    )
    parser.add_argument(
        "--allow_failure",
        action="store_true",
        help=textwrap.dedent(
            """\
            The program shall not return an errorcode if the dataset is suspicious.
            An errorcode is helpful for gitlab CI.
            """
        ),
    )
    parser.add_argument(
        "--allowed_perc",
        type=float,
        default=13.0,
        help=textwrap.dedent(
            """\
            The amount of hydrometeors that are allowed to be outside of normal,
            physical values. Interpolation can lead to that. Usually either cloud droplets
            or ice crystals have values of 8 or 12 percent. Other hydrometeors may be
            as low as 2 percent.
            """
        ),
    )
    parser.add_argument(
        "--test_nan_dims",
        action="store_true",
        help=textwrap.dedent(
            """\
            Test NaNs in the dataset for all dimensions. This should not happen and
            usually indicates a problem with the data.
            """
        ),
    )
    parser.add_argument(
        "--test_nan_vars",
        action="store_true",
        help=textwrap.dedent(
            """\
            Test NaNs in the dataset for all variables (not dimensions). This can
            happen, when a trajectory stops early, i.e., because it "crashed" into a
            mountain.
            """
        ),
    )
    parser.add_argument(
        "--test_saturation",
        action="store_true",
        help=textwrap.dedent(
            """\
            Test over saturation in the dataset for all trajectories (not dimensions). 
            """
        ),
    )
    parser.add_argument(
        "--calc_saturation",
        action="store_true",
        help=textwrap.dedent(
            """\
            During test for over saturation in the dataset for all trajectories (not dimensions):
            Recalculate the saturation.
            """
        ),
    )
    parser.add_argument(
        "--test_amounts",
        action="store_true",
        help=textwrap.dedent(
            """\
            Test if any trajectory has too large (read: physically impossible) amounts of hydrometeors.
            """
        ),
    )
    parser.add_argument(
        "--test_phases",
        action="store_true",
        help=textwrap.dedent(
            """\
            Test if all trajectories have at least once a warm phase. 
            Prints the percentages for all other phases as well.  Also tests if the phases are correct.
            """
        ),
    )
    parser.add_argument(
        "--calc_phases",
        action="store_true",
        help=textwrap.dedent(
            """\
            Recalculate the phases of each trajectory.
            """
        ),
    )
    parser.add_argument(
        "--test_physics",
        type=str,
        default="no",
        help=textwrap.dedent(
            """\
            The path of the Python interface library.
            Run (some) tests using the Python interface. Makes rudimentary
            checks on different microphysical processes with a timestep of
            '--timestep' seconds.
            """
        ),
    )
    parser.add_argument(
        "--test_sensitivities",
        action="store_true",
        help=textwrap.dedent(
            """\
            Test if any sensitivity has extremely large values excluding any number concentration sensitivities. 
            Prints the minimum, maximum, and mean sensitivities, if verbosity is set to 5. 
            """
        ),
    )
    parser.add_argument(
        "--timestep",
        type=float,
        default=30.0,
        help=textwrap.dedent(
            """\
            Timestep size in seconds for tests using the Python interface. 
            """
        ),
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=0,
        help=textwrap.dedent(
            """\
            Control the amount of outputs. Verbosity 0 prints only results of
            tests. Verbosity 1 prints in addition each error in the tests for non-physical
            hydrometeor sizes and which variables are being plotted (if any).
            Verbosity 2 prints in addition each error in the NaN test for
            dimensions. Verbosity 3 prints in addition the amount of non-physical
            hydrometeor sizes even if those are below the given threshold. 
            Prints the amount of time steps with too big sensitivities in sensitivity tests. 
            Verbosity 4 prints in addition each error in the NaN test for variables (columns).
             Prints the trajectory where large sensitivities occurr in sensitivity tests.
            Verbosity 5 prints every error for the physics tests and minimum, maximum, mean, and median sensitivities 
            in the sensitivity test.
            """
        ),
    )
    main(parser.parse_args())
