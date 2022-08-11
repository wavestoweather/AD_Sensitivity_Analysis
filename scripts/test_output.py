import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

from physics_t import physics_t

from tqdm import tqdm

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
        T_sat_low_temp = 273.15
        return p_sat_low_temp * np.exp(
            p_sat_const_a * (T - T_sat_low_temp) / (T - p_sat_const_b)
        )

    def calc_saturation(p, qv, T):
        R_universal = 8.3144598
        M_w = 0.018015265
        M_a = 0.02896546
        R_a = R_universal / M_a
        R_v = R_universal / M_w
        Epsilon = R_a / R_v
        return (p * qv) / ((Epsilon + qv) * sat_p_water(T))

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
            if "pressure" in ds_t:
                err_tmp = np.sum(
                    is_oversat(si_unit * ds_t["pressure"], ds_t["QV"], ds_t["T"])
                ).item()
            else:
                err_tmp = np.sum(
                    is_oversat(si_unit * ds_t["p"], ds_t["QV"], ds_t["T"])
                ).item()
        else:
            err_tmp = np.sum((np.asarray(ds_t["S"]) > sat_unit))
        if err_tmp > 0:
            err_traj += 1
            if verbose:
                n = len(ds_t["time"]) * len(ds_t["ensemble"])
                print(
                    f"{Warning}There are {err_tmp} ({(err_tmp/n):2.2f}%) over saturated steps in trajectory {i}{ColourReset}"
                )
        err += err_tmp
    if err > 0:
        print(
            f"{Error}Failed: Over saturation detected for {err_traj} of {n_traj} trajectories.{ColourReset}\n"
        )
    else:
        print(
            f"{Success}No over saturation for all {n_traj} trajectories detected.{ColourReset}\n"
        )
    return err, err_traj


def test_phases(ds, recalc):
    """
    For each phase, print the amount of trajectories that have at least one time step with this phase.
    ALso prints the percentage and number of time steps for each phase over all trajectories.
    Calculates the phases if necessary.
    If recalc is false, tests if the phases are correct.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    recalc : bool
        Recalculate the phases.

    Returns
    -------
    Number of trajectories without any warm phase.
    """
    err = 0
    n_total_timesteps = len(ds["trajectory"]) * len(ds["time"])
    n_trajectories = len(ds["trajectory"])

    ice_q_phase_threshold = 0.0
    ice_n_phase_threshold = 0.0
    warm_q_phase_threshold = 0.0
    warm_n_phase_threshold = 0.0
    # Extract used thresholds from C++-code
    current_dir = os.getcwd()
    if "/scripts" == current_dir[-8::]:
        constants_path = current_dir + "/../include/microphysics/constants.h"
    else:
        constants_path = current_dir + "/include/microphysics/constants.h"
    with open(constants_path, "r") as constants_f:
        for line in constants_f:
            if "ice_q_phase_threshold" in line:
                ice_q_phase_threshold = float(line.split(" = ")[-1].replace(";", ""))
            if "ice_n_phase_threshold" in line:
                ice_n_phase_threshold = float(line.split(" = ")[-1].replace(";", ""))
            if "warm_q_phase_threshold" in line:
                warm_q_phase_threshold = float(line.split(" = ")[-1].replace(";", ""))
            if "warm_n_phase_threshold" in line:
                warm_n_phase_threshold = float(line.split(" = ")[-1].replace(";", ""))

    def add_phase(ds, col_name="phase"):
        phase_col = np.full(
            (len(ds["ensemble"]), len(ds["trajectory"]), len(ds["time"])),
            "             ",
        )

        def warm(ds):
            return (
                (ds["QC"] > warm_q_phase_threshold)
                | (ds["QR"] > warm_q_phase_threshold)
                | (ds["NCCLOUD"] > warm_n_phase_threshold)
                | (ds["NCRAIN"] > warm_n_phase_threshold)
            )

        def cold(ds):
            return (
                (ds["QG"] > ice_q_phase_threshold)
                | (ds["QH"] > ice_q_phase_threshold)
                | (ds["NCGRAUPEL"] > ice_n_phase_threshold)
                | (ds["NCHAIL"] > ice_n_phase_threshold)
                | (ds["QI"] > ice_q_phase_threshold)
                | (ds["QS"] > ice_q_phase_threshold)
                | (ds["NCICE"] > ice_n_phase_threshold)
                | (ds["NCSNOW"] > ice_n_phase_threshold)
            )

        def warm_phase(ds):
            return (warm(ds)) & (~cold(ds))

        def cold_phase(ds):
            return (~warm(ds)) & (cold(ds))

        def mixed_phase(ds):
            return (warm(ds)) & (cold(ds))

        def neutral_phase(ds):
            return (~warm(ds)) & (~cold(ds))

        def nan_phase(ds):
            return np.isnan(ds["QV"])

        phase_col[np.where(warm_phase(ds))] = "warm phase   "
        phase_col[np.where(cold_phase(ds))] = "ice phase    "
        phase_col[np.where(mixed_phase(ds))] = "mixed phase  "
        phase_col[np.where(neutral_phase(ds))] = "neutral phase"
        phase_col[np.where(nan_phase(ds))] = "nan"
        ds[col_name] = (("ensemble", "trajectory", "time"), phase_col)
        return ds

    def rename_phase(ds):
        phase_col = np.full(
            (len(ds["ensemble"]), len(ds["trajectory"]), len(ds["time"])),
            "             ",
        )
        phase_col[np.where(ds["phase"] == 0)] = "warm phase   "
        phase_col[np.where(ds["phase"] == 2)] = "ice phase    "
        phase_col[np.where(ds["phase"] == 1)] = "mixed phase  "
        phase_col[np.where(ds["phase"] == 3)] = "neutral phase"
        phase_col[
            np.where(
                (ds["phase"] != 3)
                & (ds["phase"] != 2)
                & (ds["phase"] != 1)
                & (ds["phase"] != 0)
                & (ds["phase"] != 3)
            )
        ] = "nan"
        ds["phase"] = (("ensemble", "trajectory", "time"), phase_col)
        return ds

    if "phase" not in ds or recalc:
        ds = add_phase(ds)
    else:
        ds = rename_phase(ds)

    if not recalc:
        ds = add_phase(ds, "phase_reference")
        err += np.sum(ds["phase"] != ds["phase_reference"]).values.item()
        if err == 0:
            print(
                f"{Success}Phases are correct for all {n_total_timesteps} time steps{ColourReset}\n"
            )
        else:
            perc = err / n_total_timesteps
            print(
                f"{Error}Phases are not correct in {err} / {n_total_timesteps} time steps ({perc*100:2.2f}%){ColourReset}\n"
            )
            df = ds.to_dataframe()
            df2 = df.loc[df["phase"] != df["phase_reference"]]

    n_data = {phase.item(): 0 for phase in np.unique(ds["phase"])}
    n_trajs = {phase.item(): 0 for phase in np.unique(ds["phase"])}

    for phase in tqdm(np.unique(ds["phase"])):
        n = np.sum(ds["phase"] == phase)
        n_data[phase] += n.values.item()

        n_per_traj = (ds["phase"] == phase).sum(axis=2)
        zero_times = (n_per_traj == 0).sum()
        n_trajs[phase] += (n_trajectories - zero_times).values.item()

    for phase in tqdm(n_data.keys()):
        if phase == "warm phase" or phase == "warm phase   " or phase == 0:
            perc = n_trajs[phase] / n_trajectories
            if n_trajs[phase] < n_trajectories:
                print(
                    f"{Error}Failed: {phase} occurs only in {n_trajs[phase]} / {n_trajectories} trajectories ({perc*100:2.2f}%){ColourReset}\n"
                )
                err += 1
            else:
                print(
                    f"{Success}{phase} occurs in {n_trajs[phase]} / {n_trajectories} trajectories ({perc*100:2.2f}%){ColourReset}\n"
                )
        else:
            perc = n_trajs[phase] / n_trajectories
            print(
                f"{phase} occurs in {n_trajs[phase]} / {n_trajectories} trajectories ({perc*100:2.2f}%)\n"
            )
        perc = n_data[phase] / n_total_timesteps
        print(
            f"{phase} occurs in {n_data[phase]} / {n_total_timesteps} time steps ({perc*100:2.2f}%)\n"
        )
    return err


def test_graupel_melting(physics, T, qg, Ng, verbose=False):
    """
    Test graupel melting for different inputs. Graupel should be melting completely in these cases.

    Parameters
    ----------
    physics : physics_t
    T : float
        The maximum temperature for testing.
    qg : float
        The maximum power of 10 to use for graupel mass.
    Ng : float
        The maximum power of 10 to use for graupel number.
    verbose : bool
        If true: print whenever a test case fails. Otherwise print only the overall amount of fails.

    Returns
    -------
    The amount of failed test cases.
    """
    err = 0
    n_tests = 0
    err_mass = 0
    res = np.zeros(physics.get_num_comp(), dtype=np.float64)
    gradients = np.zeros(
        physics.get_num_comp() * physics.get_num_par(), dtype=np.float64
    )
    threshold = 1e-20
    threshold_n = 1
    max_error = 0.0
    max_error_n = 0.0
    mass_threshold = 5e-20
    mass_error = 0.0
    for N_i in np.arange(4, Ng, 1):
        Ng_i = 10 ** N_i
        for q_i in np.arange(-10.0, qg, 0.5):
            qg_i = 10 ** q_i
            for T_i in np.arange(280.0, T, 0.5):
                physics.graupel_melting(T_i, qg_i, Ng_i, res, gradients)
                n_tests += 1
                mass_diff = (
                    qg_i * 1e-6
                    - res[physics.index_dic["qg"]]
                    - res[physics.index_dic["qr"]]
                )
                if np.abs(mass_diff) > np.abs(mass_error):
                    mass_error = mass_diff
                if max_error_n < res[physics.index_dic["Ng"]]:
                    max_error_n = res[physics.index_dic["Ng"]]
                if max_error < res[physics.index_dic["qg"]]:
                    max_error = res[physics.index_dic["qg"]]
                if (
                    res[physics.index_dic["qg"]] >= threshold
                    or res[physics.index_dic["Ng"]] >= threshold_n
                ):
                    qg_left = res[physics.index_dic["qg"]]
                    Ng_left = res[physics.index_dic["Ng"]]
                    err += 1
                    if verbose:
                        print(
                            f"{Error}Failed: With input T={T_i}, qg={qg_i*1e-6}, Ng={Ng_i}, there is  "
                            f"qg={qg_left} and Ng={Ng_left} left over.{ColourReset}\n"
                        )
                if np.abs(mass_diff) >= mass_threshold:
                    err_mass += 1
                    if verbose:
                        print(
                            f"{Error}Failed: With input T={T_i}, qg={qg_i*1e-6}, Ng={Ng_i}, there is  "
                            f"a mass difference of {mass_diff}.{ColourReset}\n"
                        )

    if np.abs(max_error) > 0 and err == 0:
        print(
            f"{Warning}Graupel mass of {max_error} did not melt. The threshold for this process is {threshold} {ColourReset}"
        )
    if np.abs(max_error_n) > 0 and err == 0:
        print(
            f"{Warning}Graupel number of {max_error_n} did not melt. The threshold for this process is {threshold_n} {ColourReset}"
        )
    if np.abs(mass_error) > 0 and err_mass == 0:
        print(
            f"{Warning}A mass difference of {mass_error} occurred. The threshold for this process is {mass_threshold} {ColourReset}"
        )
    if err == 0 and err_mass == 0:
        print(f"{Success}Graupel melting looks good for all {n_tests} tests.\n")
    elif err_mass == 0:
        print(
            f"{Error}Failed: Graupel not completely melted for {err} of {n_tests} tests.{ColourReset}\n"
        )
    elif err == 0:
        print(
            f"{Error}Failed: A mass difference occurred for {err_mass} of {n_tests} tests.{ColourReset}\n"
        )
    else:
        print(
            f"{Error}Failed: Graupel not completely melted for {err} of {n_tests} tests and "
            f"a mass difference occured for {err_mass} of {n_tests} tests.{ColourReset}\n"
        )
    return err + err_mass


def test_ccn_act_akm(physics, w, T, qv, qc, verbose=False):
    """
    Test if the mass is preserved and if the variation for different ascend velocities is large enough.

    Parameters
    ----------
    physics : physics_t
    w : float
        Maximum ascend velocity in meters per second.
    T : float
        The maximum temperature for testing.
    qv : float
        Amount of water vapor (times 1e6).
    qc : float
        Initial amount of cloud mass in 10**qc*1e-6
    verbose : bool
        If true: print whenever a test case fails. Otherwise print only the overall amount of fails.

    Returns
    -------
    The amount of failed test cases.
    """
    Nc = 1
    err = 0
    n_tests = 0
    res = np.zeros(physics.get_num_comp(), dtype=np.float64)
    gradients = np.zeros(
        physics.get_num_comp() * physics.get_num_par(), dtype=np.float64
    )
    delta_w = (w - 0.05) / 11
    max_difference = 0.0
    diff_error = False
    threshold = 5e-20
    for i in np.arange(0.01, 1, 0.01):
        qc_list = []
        for w_i in np.arange(0.05, w, delta_w):
            physics.ccn_act_hande_akm(i, w_i, T, qv, 10 ** qc, Nc, res, gradients)
            qc_list.append(res[physics.index_dic["qc"]])
            n_tests += 1
            difference = (
                qv * 1e-6
                + 10 ** qc * 1e-6
                - res[physics.index_dic["qv"]]
                - res[physics.index_dic["qc"]]
            )
            if difference >= threshold:
                err += 1
                diff_error = True
                if verbose:
                    print(
                        f"{Error}Failed: With input T={T}, w={w_i}, qv={qv*1e-6}, p={i*1e5}. "
                        f"There is a difference of mass of {difference}.{ColourReset}\n"
                    )
            if np.abs(difference) > np.abs(max_difference):
                max_difference = difference
        n_tests += 1
        for qc_l in qc_list:
            if qc_list.count(qc_l) > 1:
                err += 1
                if verbose:
                    print(
                        f"{Error}Failed: With input T={T}, qv={qv*1e-6}, p={i*1e5}. "
                        f"There is not enough variation of qc w.r.t. to w.{ColourReset}\n"
                    )
                break
    if not diff_error and np.abs(max_difference) > 0.0:
        print(
            f"{Warning}A mass difference of {max_difference} occurred. The threshold for this process is {threshold} {ColourReset}"
        )
    if err == 0:
        print(f"{Success}CCN activation (akm) looks good for all {n_tests} tests.\n")
    else:
        print(
            f"{Error}Failed: CCN activation (akm) failed {err} of {n_tests} tests.{ColourReset}\n"
        )
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
        default=None,
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
        "--test_saturation",
        action="store_true",
        help="""
        Test over saturation in the dataset for all trajectories (not dimensions). 
        """,
    )
    parser.add_argument(
        "--calc_saturation",
        action="store_true",
        help="""
        During test for over saturation in the dataset for all trajectories (not dimensions):
        Recalculate the saturation.
        """,
    )
    parser.add_argument(
        "--test_phases",
        action="store_true",
        help="""
        Test if all trajectories have at least once a warm phase. Prints the percentages for all other phases as well.  Also tests if the phases are correct.
        """,
    )
    parser.add_argument(
        "--calc_phases",
        action="store_true",
        help="""
        Recalculate the phases of each trajectory.
        """,
    )
    parser.add_argument(
        "--test_physics",
        type=str,
        default="no",
        help="""
        The path of the Python interface library.
        Run (some) tests using the Python interface. Makes rudimentary
        checks on different microphysical processes with a timestep of
        30 seconds.
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
        Verbosity 5 prints every error for the physics tests.
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
        if args.test_saturation:
            _, err_traj = test_saturation(ds, args.calc_saturation, args.verbosity)
            errors += err_traj
        if args.test_phases:
            errors += test_phases(ds, args.calc_phases)
        return errors

    if args.input is not None:
        err = 0
        if args.input[-1] == "/":
            files = [
                f for f in os.listdir(args.input) if os.path.isfile(args.input + f)
            ]
            n = len(files)
            for i, f in enumerate(files):
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
            print(f"\n{Success}No errors occured for {n_files} file(s).{ColourReset}")
        if not args.allow_failure and err > 0:
            raise Exception(
                f"{Error}Failed: {err} errors occured during testing for {n_files} file(s). Check the output!{ColourReset}"
            )
        elif err > 0:
            print(
                f"{Error}Failed: {err} errors occured during testing for {n_files} file(s). Check the output!{ColourReset}"
            )

    if args.test_physics != "no":
        err = 0
        n_processes = 0
        physics = physics_t(lib_path=args.test_physics)
        physics.setup_model_constants(30.0, 0.1)
        print(f"{Status}~*~*~*~Testing graupel melting~*~*~*~{ColourReset}")
        tmp = test_graupel_melting(
            physics=physics,
            T=285,
            qg=2,
            Ng=8,
            verbose=(args.verbosity > 4),
        )
        if tmp > 0:
            n_processes += 1
        err += tmp

        print(f"{Status}~*~*~*~Testing CCN activation (akm)~*~*~*~{ColourReset}")
        tmp = test_ccn_act_akm(
            physics=physics,
            w=3,
            T=270,
            qv=1e2,
            qc=-6,
            verbose=(args.verbosity > 4),
        )
        if tmp > 0:
            n_processes += 1
        err += tmp

        if err == 0:
            print(f"\n{Success}No errors occured for all physics tests.{ColourReset}")
        if not args.allow_failure and err > 0:
            raise Exception(
                f"{Error}Failed: {err} errors occured during testing for {n_processes} processes.{ColourReset}"
            )
        elif err > 0:
            print(
                f"{Error}Failed: {err} errors occured during testing for {n_processes} processes.{ColourReset}"
            )
