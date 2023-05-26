"""Test C++ physics processes.

"""
import numpy as np

from ad_sensitivity_analysis.tests.color_scheme import (
    COLOR_RESET,
    WARNING,
    SUCCESS,
    ERROR,
)


def print_graupel_melting_results(error_counts, thresholds, n_tests):
    """

    Parameters
    ----------
    error_counts
    thresholds
    n_tests

    Returns
    -------

    """
    if np.abs(error_counts["max_error"]) > 0 and error_counts["err"] == 0:
        print(
            f"{WARNING}Graupel mass of {error_counts['max_error']} did not melt. "
            f"The threshold for this process is {thresholds['thresh_q']} {COLOR_RESET}"
        )
    if np.abs(error_counts["max_error_n"]) > 0 and error_counts["err"] == 0:
        print(
            f"{WARNING}Graupel number of {error_counts['max_error_n']} did not melt. "
            f"The threshold for this process is {thresholds['thresh_n']} {COLOR_RESET}"
        )
    if np.abs(error_counts["mass_error"]) > 0 and error_counts["err_mass"] == 0:
        print(
            f"{WARNING}A mass difference of {error_counts['mass_error']} occurred. "
            f"The threshold for this process is {thresholds['thresh_mass']} {COLOR_RESET}"
        )
    if error_counts["err"] == 0 and error_counts["err_mass"] == 0:
        print(f"{SUCCESS}Graupel melting looks good for all {n_tests} tests.\n")
    elif error_counts["err_mass"] == 0:
        print(
            f"{ERROR}Failed: Graupel not completely melted for {error_counts['err']} of {n_tests} tests.{COLOR_RESET}\n"
        )
    elif error_counts["err"] == 0:
        print(
            f"{ERROR}Failed: A mass difference occurred for {error_counts['err_mass']} "
            f"of {n_tests} tests.{COLOR_RESET}\n"
        )
    else:
        print(
            f"{ERROR}Failed: Graupel not completely melted for {error_counts['err']} of {n_tests} tests and "
            f"a mass difference occured for {error_counts['err_mass']} of {n_tests} tests.{COLOR_RESET}\n"
        )


def test_graupel_melting(physics, T, qg, Ng, verbose=False):
    """
    Test graupel melting for different inputs. Graupel should be melting completely in these cases.

    Parameters
    ----------
    physics : PhysicsT
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
    n_tests = 0
    physics_output = {
        "res": np.zeros(physics.get_num_comp(), dtype=np.float64),
        "gradients": np.zeros(
            physics.get_num_comp() * physics.get_num_par(), dtype=np.float64
        ),
    }
    thresholds = {
        "thresh_q": 1e-20,
        "thresh_n": 1,
        "thresh_mass": 5e-20,
    }
    error_counts = {
        "err": 0,
        "err_mass": 0,
        "max_error": 0.0,
        "max_error_n": 0.0,
        "mass_error": 0.0,
    }
    for n_i_exponent in np.arange(4, Ng, 1):
        Ng_i = 10**n_i_exponent
        for q_i in np.arange(-10.0, qg, 0.5):
            qg_i = 10**q_i
            for T_i in np.arange(280.0, T, 0.5):
                physics.graupel_melting(
                    {
                        "T": T_i,
                        "qg": qg_i,
                        "Ng": Ng_i,
                    },
                    **physics_output,
                )
                n_tests += 1
                mass_diff = (
                    qg_i * 1e-6
                    - physics_output["res"][physics.index_dic["qg"]]
                    - physics_output["res"][physics.index_dic["qr"]]
                )
                if np.abs(mass_diff) > np.abs(error_counts["mass_error"]):
                    error_counts["mass_error"] = mass_diff
                if (
                    error_counts["max_error_n"]
                    < physics_output["res"][physics.index_dic["Ng"]]
                ):
                    error_counts["max_error_n"] = physics_output["res"][
                        physics.index_dic["Ng"]
                    ]
                if (
                    error_counts["max_error"]
                    < physics_output["res"][physics.index_dic["qg"]]
                ):
                    error_counts["max_error"] = physics_output["res"][
                        physics.index_dic["qg"]
                    ]
                if (
                    physics_output["res"][physics.index_dic["qg"]]
                    >= thresholds["thresh_q"]
                    or physics_output["res"][physics.index_dic["Ng"]]
                    >= thresholds["thresh_n"]
                ):
                    error_counts["err"] += 1
                    if verbose:
                        print(
                            f"{ERROR}Failed: With input T={T_i}, qg={qg_i*1e-6}, Ng={Ng_i}, there is  "
                            f"qg={physics_output['res'][physics.index_dic['qg']]} and "
                            f"Ng={physics_output['res'][physics.index_dic['Ng']]} left over.{COLOR_RESET}\n"
                        )
                if np.abs(mass_diff) >= thresholds["thresh_mass"]:
                    error_counts["err_mass"] += 1
                    if verbose:
                        print(
                            f"{ERROR}Failed: With input T={T_i}, qg={qg_i*1e-6}, Ng={Ng_i}, there is  "
                            f"a mass difference of {mass_diff}.{COLOR_RESET}\n"
                        )

    print_graupel_melting_results(error_counts, thresholds, n_tests)
    return error_counts["err"] + error_counts["err_mass"]


def print_ccn_act_akm_results(thresholds_errors, n_tests):
    """

    Parameters
    ----------
    thresholds_errors
    n_tests

    Returns
    -------

    """
    if thresholds_errors["not_available"]:
        print(
            "'ccn_act_hande_akm' is not available. Skipping tests. "
            "Please compile the interface with 'CCN_AKM' for this test."
        )
        return
    if (
        not thresholds_errors["diff_error"]
        and np.abs(thresholds_errors["max_difference"]) > 0.0
    ):
        print(
            f"{WARNING}A mass difference of {thresholds_errors['max_difference']} occurred. "
            f"The threshold for this process is {thresholds_errors['threshold']} {COLOR_RESET}"
        )
    if thresholds_errors["err"] == 0:
        print(f"{SUCCESS}CCN activation (akm) looks good for all {n_tests} tests.\n")
    else:
        print(
            f"{ERROR}Failed: CCN activation (akm) failed {thresholds_errors['err']} of {n_tests} tests.{COLOR_RESET}\n"
        )


def test_ccn_act_akm(physics, w, T, qv, qc, verbose=False):
    """
    Test if the mass is preserved and if the variation for different ascend velocities is large enough.

    Parameters
    ----------
    physics : PhysicsT
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
    n_tests = 0
    physics_output = {
        "res": np.zeros(physics.get_num_comp(), dtype=np.float64),
        "gradients": np.zeros(
            physics.get_num_comp() * physics.get_num_par(), dtype=np.float64
        ),
    }
    delta_w = (w - 0.05) / 11
    thresholds_errors = {
        "threshold": 5e-20,
        "not_available": False,
        "diff_error": False,
        "max_difference": 0.0,
        "err": 0,
    }
    for i in np.arange(0.01, 1, 0.01):
        qc_list = []
        for w_i in np.arange(0.05, w, delta_w):
            if not physics.ccn_act_hande_akm(
                {
                    "qv": qv,
                    "qc": 10**qc,
                    "Nc": 1,
                    "T": T,
                    "p": i * 1e5,
                    "w": w_i,
                },
                **physics_output,
            ):
                thresholds_errors["not_available"] = True
                break
            qc_list.append(physics_output["res"][physics.index_dic["qc"]])
            n_tests += 1
            difference = (
                qv * 1e-6
                + 10**qc * 1e-6
                - physics_output["res"][physics.index_dic["qv"]]
                - physics_output["res"][physics.index_dic["qc"]]
            )
            if difference >= thresholds_errors["threshold"]:
                thresholds_errors["err"] += 1
                thresholds_errors["diff_error"] = True
                if verbose:
                    print(
                        f"{ERROR}Failed: With input T={T}, w={w_i}, qv={qv*1e-6}, p={i*1e5}. "
                        f"There is a difference of mass of {difference}.{COLOR_RESET}\n"
                    )
            if np.abs(difference) > np.abs(thresholds_errors["max_difference"]):
                thresholds_errors["max_difference"] = difference
        if thresholds_errors["not_available"]:
            break
        n_tests += 1
        for qc_l in qc_list:
            if qc_list.count(qc_l) > 1:
                thresholds_errors["err"] += 1
                if verbose:
                    print(
                        f"{ERROR}Failed: With input T={T}, qv={qv*1e-6}, p={i*1e5}. "
                        f"There is not enough variation of qc w.r.t. to w.{COLOR_RESET}\n"
                    )
                break
    print_ccn_act_akm_results(thresholds_errors, n_tests)
    return thresholds_errors["err"]
