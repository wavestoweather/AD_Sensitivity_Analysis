import ctypes
import numpy as np

import holoviews as hv
import pandas as pd
import seaborn as sns


class physics_t(object):
    c_double_arr = None
    clib = None
    index_dic = {
        "p": 0,
        "T": 1,
        "w": 2,
        "S": 3,
        "qc": 4,
        "qr": 5,
        "qv": 6,
        "Nc": 7,
        "Nr": 8,
        "qi": 9,
        "Ni": 10,
        "qs": 11,
        "Ns": 12,
        "qg": 13,
        "Ng": 14,
        "qh": 15,
        "Nh": 16,
        "qi_out": 17,
        "qs_out": 18,
        "qr_out": 19,
        "qg_out": 20,
        "qh_out": 21,
        "lat_heat": 22,
        "lat_cool": 23,
        "Ni_out": 24,
        "Ns_out": 25,
        "Nr_out": 26,
        "Ng_out": 27,
        "Nh_out": 28,
        "z": 29,
        "n_inact": 30,
        "depo": 31,
        "sub": 32,
    }

    def __init__(
        self,
        lib_path="../cmake-build-python_interface/lib/libpython_interface.so",
        table_path="../dmin_wetgrowth_lookup.dat",
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

        self.obj = self.clib.physics_t_new(table_path.encode("utf-8"))

    def get_num_comp(self):
        return self.clib.physics_t_get_num_comp(self.obj)

    def get_num_par(self):
        return self.clib.physics_t_get_num_par(self.obj)

    def set_ref_quants(self, qref, pref, wref, tref, zref, Nref, timeref):
        self.clib.physics_t_set_ref_quants(qref, pref, wref, tref, zref, Nref, timeref)

    def setup_model_constants(self, dt_prime, uncertainty_perc):
        self.clib.physics_t_setup_model_constants(self.obj, dt_prime, uncertainty_perc)

    def setup_model_constants_uncert(self, uncertainty_perc):
        self.clib.physics_t_setup_model_constants_uncert(self.obj, uncertainty_perc)

    def ccn_act_hande_akm(self, p, w, T, qv, qc, Nc, res, gradients):
        self.clib.physics_t_py_ccn_act_hande_akm(
            self.obj,
            p,
            w,
            T,
            qv,
            qc,
            Nc,
            res.ctypes.data_as(self.c_double_arr),
            gradients.ctypes.data_as(self.c_double_arr),
        )

    def graupel_melting(self, T, qg, Ng, res, gradients):
        self.clib.physics_t_py_graupel_melting(
            self.obj,
            T,
            qg,
            Ng,
            res.ctypes.data_as(self.c_double_arr),
            gradients.ctypes.data_as(self.c_double_arr),
        )

    def graupel_melting_plot(self, T, qg, Ng, x, y, hue, logy):
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
        for N_i in np.arange(0, Ng, 1):
            Ng_i = 10 ** N_i
            for q_i in np.arange(-10.0, qg, 0.5):
                qg_i = 10 ** q_i
                for T_i in np.arange(270.0, T, 0.5):
                    self.graupel_melting(T_i, qg_i, Ng_i, res, gradients)
                    data["T"].append(T_i)
                    data["qg"].append(qg_i * 1.0e-6)
                    data["Ng"].append(Ng_i)
                    data["final qg"].append(res[self.index_dic["qg"]])
                    data["final Ng"].append(res[self.index_dic["Ng"]])
                    data["final qr"].append(res[self.index_dic["qr"]])
                    data["final Nr"].append(res[self.index_dic["Nr"]])
        df = pd.DataFrame.from_dict(data)
        g = sns.scatterplot(data=df, x=x, y=y, hue=hue)
        if logy:
            g.axes.set_yscale("log")
        return g, df

    def ccn_act_hande_akm_plot(self, w, T, qv, qc, y, logy):
        Nc = 1
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        data = {"p": [], y: [], "w": []}
        delta_w = (w - 0.05) / 11
        for i in np.arange(0.01, 1, 0.01):
            for w_i in np.arange(0.05, w, delta_w):
                self.ccn_act_hande_akm(i, w_i, T, qv, 10 ** qc, Nc, res, gradients)
                data["p"].append(i * 1e3)
                data[y].append(res[self.index_dic[y]])
                data["w"].append(w_i)
        df = pd.DataFrame.from_dict(data)
        g = sns.scatterplot(data=df, x="p", y=y, hue="w")
        if logy:
            g.axes.set_yscale("log")
        return g, df


if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(
        description="""
        Load the C-library which provides an interface to microphysics processes. Plot some examples of these 
        processes.
        """,
    )
    parser.add_argument(
        "--library_path",
        type=str,
        required=True,
        help="""
        Path to the C-library. You may find it at cmake-build-python_interface/lib/libpython_interface.so.
        """,
    )
    parser.add_argument(
        "--timestep",
        type=float,
        default=30.0,
        help="""
        The timestep size in seconds.
        """,
    )
    parser.add_argument(
        "--uncertainty",
        type=float,
        default=0.01,
        help="""
        The uncertainty that is used for all model parameters.
        """,
    )
    parser.add_argument(
        "--store_path",
        default="../pics/microphysics",
        type=str,
        help="""
        Path to store the generated images.
        """,
    )
    args = parser.parse_args()
    phys = physics_t(lib_path=args.library_path)
    phys.setup_model_constants(30.0, 0.1)
    renderer = hv.Store.renderers["matplotlib"].instance(fig="png", dpi=300)

    plot, _ = phys.graupel_melting_plot(
        T=280,
        qg=2,
        Ng=8,
        x="qg",
        y="final qr",
        hue="T",
    )

    def save_plot(plot, store_path):
        i = 0
        save = store_path + "_" + "{:03d}".format(i)
        while os.path.isfile(save + ".png"):
            i = i + 1
            save = store_path + "_" + "{:03d}".format(i)
        renderer.save(plot, save)

    save_plot(plot, args.store_path)
    plot, _ = phys.ccn_act_hande_akm_plot(
        w=10,
        T=270,
        qv=1e2,
        qc=-6,
        y="qc",
        logy=True,
    )
    save_plot(plot, args.store_path)
