"""Interface to C++ physics.

"""
import ctypes
import numpy as np

import holoviews as hv
import pandas as pd
import seaborn as sns

from ad_sensitivity_analysis.interface import aux_interface


class PhysicsT:
    """
    Interface to C++ physics.
    """

    c_double_arr = None
    clib = None
    index_dic = aux_interface.get_index_dic()
    index_in_dic = None

    def __init__(
        self,
        lib_path="../cmake-build-python_interface/lib/libpython_interface.so",
        table_path="../dmin_wetgrowth_lookup.dat",
        b_eight=False,
    ):
        hv.extension("bokeh")
        sns.set(rc={"figure.figsize": (11.7, 8.27)})
        self.clib = ctypes.cdll.LoadLibrary(lib_path)

        self.c_double_arr = ctypes.POINTER(ctypes.c_double)
        self.clib.physics_t_new.argtypes = None
        self.clib.physics_t_new.restype = ctypes.c_void_p

        self.clib.physics_t_get_num_comp.argtypes = [ctypes.c_void_p]
        self.clib.physics_t_get_num_comp.restyp = ctypes.c_int

        self.clib.physics_t_get_num_par.argtypes = [ctypes.c_void_p]
        self.clib.physics_t_get_num_par.restyp = ctypes.c_int

        self.clib.physics_t_set_ref_quants.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
        ]
        self.clib.physics_t_set_ref_quants.restyp = ctypes.c_void_p

        self.clib.physics_t_setup_model_constants.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_double,
        ]
        self.clib.physics_t_setup_model_constants.restyp = ctypes.c_void_p

        self.clib.physics_t_setup_model_constants_uncert.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
        ]
        self.clib.physics_t_setup_model_constants_uncert.restyp = ctypes.c_void_p

        if hasattr(self.clib, "physics_t_py_ccn_act_hande_akm"):
            self.clib.physics_t_py_ccn_act_hande_akm.argtypes = [
                ctypes.c_void_p,
                ctypes.c_double,
                ctypes.c_double,
                ctypes.c_double,
                ctypes.c_double,
                ctypes.c_double,
                ctypes.c_double,
                self.c_double_arr,
                self.c_double_arr,
            ]
            self.clib.physics_t_py_ccn_act_hande_akm.restyp = ctypes.c_void_p

        self.clib.physics_t_py_graupel_melting.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            self.c_double_arr,
            self.c_double_arr,
        ]
        self.clib.physics_t_py_graupel_melting.restyp = ctypes.c_void_p

        self.clib.physics_t_py_riming_ice.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            self.c_double_arr,
            self.c_double_arr,
        ]
        self.clib.physics_t_py_riming_ice.restyp = ctypes.c_void_p

        self.clib.physics_t_py_riming_snow.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            self.c_double_arr,
            self.c_double_arr,
        ]
        self.clib.physics_t_py_riming_snow.restyp = ctypes.c_void_p

        self.clib.physics_t_py_riming_with_depo.argtypes = [
            ctypes.c_void_p,
            ctypes.c_double,  # qv
            ctypes.c_double,  # qc
            ctypes.c_double,  # Nc
            ctypes.c_double,  # qi
            ctypes.c_double,  # Ni
            ctypes.c_double,  # qs
            ctypes.c_double,  # Ns
            ctypes.c_double,  # qr
            ctypes.c_double,  # Nr
            ctypes.c_double,  # T
            ctypes.c_double,  # p
            self.c_double_arr,
            self.c_double_arr,
        ]
        self.clib.physics_t_py_riming_with_depo.restyp = ctypes.c_void_p

        self.obj = self.clib.physics_t_new(table_path.encode("utf-8"))
        offset = 0
        if b_eight:
            offset += 9
        self.index_in_dic = aux_interface.get_index_in_dic(offset)

    @staticmethod
    def _create_plot(data, x, y, hue=None, logy=False):
        """

        Parameters
        ----------
        data : Dict or pandas.DataFrame
        x
        y
        hue
        logy

        Returns
        -------

        """
        if not isinstance(data, pd.DataFrame):
            data = pd.DataFrame.from_dict(data)
        plot_obj = sns.scatterplot(data=data, x=x, y=y, hue=hue)
        if logy:
            plot_obj.axes.set_yscale("log")
        return plot_obj, data

    @staticmethod
    def _set_vals(param, n_bins, param_dict):
        """

        Parameters
        ----------
        param : string

        n_bins
        param_dict

        Returns
        -------

        """
        return np.arange(
            param_dict[param][0],
            param_dict[param][1],
            (param_dict[param][1] - param_dict[param][0]) / (n_bins + 1),
        )

    def _get_gradient(self, out_param, in_param, gradients):
        """

        Parameters
        ----------
        out_param
        in_param
        gradients

        Returns
        -------

        """
        return gradients[
            self.index_dic[out_param] * self.get_num_par() + self.index_in_dic[in_param]
        ]

    def get_num_comp(self):
        """

        Returns
        -------

        """
        return self.clib.physics_t_get_num_comp(self.obj)

    def get_num_par(self):
        """

        Returns
        -------

        """
        return self.clib.physics_t_get_num_par(self.obj)

    def set_ref_quants(self, qref, pref, wref, tref, zref, n_ref, timeref):
        """

        Parameters
        ----------
        qref
        pref
        wref
        tref
        zref
        n_ref
        timeref

        Returns
        -------

        """
        self.clib.physics_t_set_ref_quants(qref, pref, wref, tref, zref, n_ref, timeref)

    def setup_model_constants(self, dt_prime, uncertainty_perc):
        """

        Parameters
        ----------
        dt_prime
        uncertainty_perc

        Returns
        -------

        """
        self.clib.physics_t_setup_model_constants(self.obj, dt_prime, uncertainty_perc)

    def setup_model_constants_uncert(self, uncertainty_perc):
        """

        Parameters
        ----------
        uncertainty_perc

        Returns
        -------

        """
        self.clib.physics_t_setup_model_constants_uncert(self.obj, uncertainty_perc)

    def ccn_act_hande_akm(self, param_dict, res, gradients):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "qv", "qc", "Nc", "T", "p", "w". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        res
        gradients

        Returns
        -------
        True if the function exists and has been used, false otherwise.
        """
        if hasattr(self.clib, "physics_t_py_ccn_act_hande_akm"):
            self.clib.physics_t_py_ccn_act_hande_akm(
                self.obj,
                param_dict["p"],
                param_dict["w"],
                param_dict["T"],
                param_dict["qv"],
                param_dict["qc"],
                param_dict["Nc"],
                res.ctypes.data_as(self.c_double_arr),
                gradients.ctypes.data_as(self.c_double_arr),
            )
            return True
        return False

    def graupel_melting(self, param_dict, res, gradients):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "qg", "Ng", "T". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        res
        gradients

        Returns
        -------

        """
        self.clib.physics_t_py_graupel_melting(
            self.obj,
            param_dict["T"],
            param_dict["qg"],
            param_dict["Ng"],
            res.ctypes.data_as(self.c_double_arr),
            gradients.ctypes.data_as(self.c_double_arr),
        )

    def riming_ice(self, param_dict, res, gradients):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "qc", "Nc", "qi", "Ni", "qr", "Nr", "T". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        res
        gradients

        Returns
        -------

        """
        self.clib.physics_t_py_riming_ice(
            self.obj,
            param_dict["qc"],
            param_dict["Nc"],
            param_dict["qi"],
            param_dict["Ni"],
            param_dict["qr"],
            param_dict["Nr"],
            param_dict["T"],
            res.ctypes.data_as(self.c_double_arr),
            gradients.ctypes.data_as(self.c_double_arr),
        )

    def riming_snow(self, param_dict, res, gradients):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "qc", "Nc", "qs", "Ns", "qr", "Nr", "T". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        res
        gradients

        Returns
        -------

        """
        self.clib.physics_t_py_riming_snow(
            self.obj,
            param_dict["qc"],
            param_dict["Nc"],
            param_dict["qs"],
            param_dict["Ns"],
            param_dict["qr"],
            param_dict["Nr"],
            param_dict["T"],
            res.ctypes.data_as(self.c_double_arr),
            gradients.ctypes.data_as(self.c_double_arr),
            **param_dict,
        )

    def riming_with_depo(self, param_dict, res, gradients):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "qv", "qc", "Nc", "qi", "Ni",
            "qs", "Ns", "qr", "Nr", "T", "p". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        res
        gradients

        Returns
        -------

        """
        self.clib.physics_t_py_riming_with_depo(
            self.obj,
            param_dict["qv"],
            param_dict["qc"],
            param_dict["Nc"],
            param_dict["qi"],
            param_dict["Ni"],
            param_dict["qs"],
            param_dict["Ns"],
            param_dict["qr"],
            param_dict["Nr"],
            param_dict["T"],
            param_dict["p"],
            res.ctypes.data_as(self.c_double_arr),
            gradients.ctypes.data_as(self.c_double_arr),
        )

    def graupel_melting_plot(self, param_dict, x, y, hue, logy):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "T", "qg", "Ng". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        x
        y
        hue
        logy

        Returns
        -------

        """
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        data = {
            "qg": [],
            "Ng": [],
            "T": [],
            "final Ng": [],
            "final qg": [],
            "final Nr": [],
            "final qr": [],
        }
        for n_i in np.arange(0, param_dict["Ng"], 1):
            param_dict["Ng"] = 10**n_i
            for q_i in np.arange(-10.0, param_dict["qg"], 0.5):
                param_dict["qg"] = 10**q_i
                for T_i in np.arange(270.0, param_dict["T"], 0.5):
                    param_dict["T"] = T_i
                    self.graupel_melting(param_dict, res, gradients)
                    data["T"].append(param_dict["T"])
                    data["qg"].append(param_dict["qg"] * 1.0e-6)
                    data["Ng"].append(param_dict["Ng"])
                    data["final qg"].append(res[self.index_dic["qg"]])
                    data["final Ng"].append(res[self.index_dic["Ng"]])
                    data["final qr"].append(res[self.index_dic["qr"]])
                    data["final Nr"].append(res[self.index_dic["Nr"]])
        return self._create_plot(data=data, x=x, y=y, hue=hue, logy=logy)

    def ccn_act_hande_akm_plot(self, param_dict, y, logy):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "w", "T", "qv", "qc". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        y
        logy

        Returns
        -------

        """
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        data = {"p": [], y: [], "w": []}
        param_dict["qc"] = 10 ** param_dict["qc"]
        param_dict["Nc"] = 1
        for i in np.arange(0.01, 1, 0.01):
            param_dict["p"] = i * 1e5
            for w_i in np.arange(0.05, param_dict["w"], (param_dict["w"] - 0.05) / 11):
                param_dict["w"] = w_i
                self.ccn_act_hande_akm(
                    param_dict,
                    res,
                    gradients,
                )
                data["p"].append(i * 1e5)
                data[y].append(res[self.index_dic[y]])
                data["w"].append(w_i)
        return self._create_plot(data=data, x="p", y=y, hue="w", logy=logy)

    def riming_ice_plot(
        self, param_dict, x, y, hue, logy, nx=10, nhue=10, out_param=""
    ):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "qc", "Nc", "qi", "Ni", "qr", "Nr", "T". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        x
        y
        hue
        logy
        nx
        nhue
        out_param

        Returns
        -------

        """
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        data = {x: [], y: [], hue: []}
        hue_vals = self._set_vals(hue, nhue, param_dict)
        for x_v in self._set_vals(x, nx, param_dict):
            param_dict[x] = x_v
            for h_v in hue_vals:
                param_dict[hue] = h_v
                self.riming_ice(
                    param_dict=param_dict,
                    res=res,
                    gradients=gradients,
                )
                data[x].append(x_v)
                if "final " in y:
                    data[y].append(res[self.index_dic[y[6::]]])
                else:
                    data[y].append(self._get_gradient(out_param, y, gradients))
                data[hue].append(h_v)
        return self._create_plot(data=data, x=x, y=y, hue=hue, logy=logy)

    def riming_snow_plot(
        self, param_dict, x, y, hue, logy, nx=10, nhue=10, out_param=""
    ):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "qc", "Nc", "qs", "Ns", "qr", "Nr", "T". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        x
        y
        hue
        logy
        nx
        nhue
        out_param

        Returns
        -------

        """
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        data = {x: [], y: [], hue: []}
        hue_vals = self._set_vals(hue, nhue, param_dict)
        for x_v in self._set_vals(x, nx, param_dict):
            param_dict[x] = x_v
            for h_v in hue_vals:
                param_dict[hue] = h_v
                self.riming_snow(
                    param_dict=param_dict,
                    res=res,
                    gradients=gradients,
                )
                data[x].append(x_v)
                if "final " in y:
                    data[y].append(res[self.index_dic[y[6::]]])
                else:
                    data[y].append(self._get_gradient(out_param, y, gradients))
                data[hue].append(h_v)
        return self._create_plot(data=data, x=x, y=y, hue=hue, logy=logy)

    def riming_with_depo_plot(
        self,
        param_dict,
        x,
        y,
        hue,
        logy,
        nx=10,
        nhue=10,
        out_param="",
    ):
        """

        Parameters
        ----------
        param_dict : Dictionary of np.arrays and floats
            Dictionary with keys "qv", "qc", "Nc", "qi", "Ni",
            "qs", "Ns", "qr", "Nr", "T", "p". Values
            are either lower and upper bounds for the binning or one float with fixed
            values.
        x
        y
        hue
        logy
        nx
        nhue
        out_param

        Returns
        -------

        """
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        data = {x: [], y: [], hue: []}
        hue_vals = self._set_vals(hue, nhue, param_dict)
        for x_v in self._set_vals(x, nx, param_dict):
            param_dict[x] = x_v
            for h_v in hue_vals:
                param_dict[hue] = h_v
                self.riming_with_depo(param_dict, res, gradients)
                data[x].append(x_v)
                if "final " in y:
                    data[y].append(res[self.index_dic[y[6::]]])
                else:
                    data[y].append(self._get_gradient(out_param, y, gradients))
                data[hue].append(h_v)
        return self._create_plot(data=data, x=x, y=y, hue=hue, logy=logy)
