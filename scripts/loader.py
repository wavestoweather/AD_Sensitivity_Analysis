# from iris.analysis.cartography import rotate_pole
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from progressbar import progressbar as pb
import sys
import os
import xarray as xr
from multiprocessing import Pool
from itertools import repeat


params_dict = {"p": "_diff_0.txt", "T": "_diff_1.txt",
               "w": "_diff_2.txt", "S": "_diff_3.txt", "qc": "_diff_4.txt",
               "qr": "_diff_5.txt", "qv": "_diff_6.txt", "Nc": "_diff_7.txt",
               "Nr": "_diff_8.txt", "qi": "_diff_9.txt",
               "Ni": "_diff_10.txt", "qs": "_diff_11.txt",
               "Ns": "_diff_12.txt", "qg": "_diff_13.txt",
               "Ng": "_diff_14.txt", "qh": "_diff_15.txt",
               "Nh": "_diff_16.txt", "qiout": "_diff_17.txt",
               "qsout": "_diff_18.txt", "qrout": "_diff_19.txt",
               "qgout": "_diff_20.txt", "qhout": "_diff_21.txt",
               "latent_heat": "_diff_22.txt",
               "latent_cool": "_diff_23.txt"}
params_dict2 = {"_diff_0.txt": "p",
               "_diff_1.txt": "T",
               "_diff_2.txt": "w",
               "_diff_3.txt": "S",
               "_diff_4.txt": "qc",
               "_diff_5.txt": "qr",
               "_diff_6.txt": "qv",
               "_diff_7.txt": "Nc",
               "_diff_8.txt": "Nr",
               "_diff_9.txt": "qi",
               "_diff_10.txt": "Ni",
               "_diff_11.txt": "qs",
               "_diff_12.txt": "Ns",
               "_diff_13.txt": "qg",
               "_diff_14.txt": "Ng",
               "_diff_15.txt": "qh",
               "_diff_16.txt": "Nh",
               "_diff_17.txt": "qiout",
               "_diff_18.txt": "qsout",
               "_diff_19.txt": "qrout",
               "_diff_20.txt": "qgout",
               "_diff_21.txt": "qhout",
               "_diff_22.txt": "latent_heat",
               "_diff_23.txt": "latent_cool",
               "_diff_24.txt": "Niout",
               "_diff_25.txt": "Nsout",
               "_diff_26.txt": "Nrout",
               "_diff_27.txt": "Ngout",
               "_diff_28.txt": "Nhout",
               "_diff_29.txt": "z",
               "_diff_30.txt": "n_inact",
               "_diff_31.txt": "depo",
               "_diff_32.txt": "sub"}
params_dict2_met3d = {"_diff_0.txt": "pressure",
               "_diff_1.txt": "T",
               "_diff_2.txt": "w",
               "_diff_3.txt": "S",
               "_diff_4.txt": "QC",
               "_diff_5.txt": "QR",
               "_diff_6.txt": "QV",
               "_diff_7.txt": "NCCLOUD",
               "_diff_8.txt": "NCRAIN",
               "_diff_9.txt": "QI",
               "_diff_10.txt": "NCICE",
               "_diff_11.txt": "QS",
               "_diff_12.txt": "NCSNOW",
               "_diff_13.txt": "QG",
               "_diff_14.txt": "NCGRAUPEL",
               "_diff_15.txt": "QH",
               "_diff_16.txt": "NCHAIL",
               "_diff_17.txt": "QI_OUT",
               "_diff_18.txt": "QS_OUT",
               "_diff_19.txt": "QR_OUT",
               "_diff_20.txt": "QG_OUT",
               "_diff_21.txt": "QH_OUT",
               "_diff_22.txt": "latent_heat",
               "_diff_23.txt": "latent_cool",
               "_diff_24.txt": "NI_OUT",
               "_diff_25.txt": "NS_OUT",
               "_diff_26.txt": "NR_OUT",
               "_diff_27.txt": "NG_OUT",
               "_diff_28.txt": "NH_OUT",
               "_diff_29.txt": "z",
               "_diff_30.txt": "n_inact",
               "_diff_31.txt": "depo",
               "_diff_32.txt": "sub"}
met3d_rename_dic = {"p"         : "pressure",
                    "LONGITUDE" : "lon",
                    "LATITUDE"  : "lat",
                    "MAP"       : "WCB_flag",
                    "qv"        : "QV",
                    "qc"        : "QC",
                    "qr"        : "QR",
                    "qi"        : "QI",
                    "qs"        : "QS",
                    "qg"        : "QG",
                    "qh"        : "QH",
                    "Nrout"     : "NR_OUT",
                    "Niout"     : "NI_OUT",
                    "Nsout"     : "NS_OUT",
                    "Ngout"     : "NG_OUT",
                    "Nhout"     : "NH_OUT",
                    "Nc"        : "NCCLOUD",
                    "Ng"        : "NCGRAUPEL",
                    "Nr"        : "NCRAIN",
                    "Nh"        : "NCHAIL",
                    "Ns"        : "NCSNOW",
                    "Ni"        : "NCICE",
                    "qrout"     : "QR_OUT",
                    "qiout"     : "QI_OUT",
                    "qsout"     : "QS_OUT",
                    "qgout"     : "QG_OUT",
                    "qhout"     : "QH_OUT"}
deriv_type_dic = {
    "step": np.uint64,
    "time": np.float64,
    "timestep": np.float64,
     "trajectory": np.uint64,
     "Output Parameter": "category",
     "lon": np.float64,
     "lat": np.float64,
     "time_after_ascent": np.float64,
     "WCB_flag": np.bool_,
     "dp2h": np.bool_,
     "conv_400": np.bool_,
     "conv_600": np.bool_,
     "slan_400": np.bool_,
     "slan_600": np.bool_,
     "da_1": np.float64,
     "da_2": np.float64,
     "de_1": np.float64,
     "de_2": np.float64,
     "dd": np.float64,
     "dN_c": np.float64,
     "dgamma": np.float64,
     "dbeta_c": np.float64,
     "dbeta_r": np.float64,
     "ddelta1": np.float64,
     "ddelta2": np.float64,
     "dzeta": np.float64,
     "drain_gfak": np.float64,
     "dcloud_k_au": np.float64,
     "dcloud_k_sc": np.float64,
     "dkc_autocon": np.float64,
     "dinv_z": np.float64,
     "drain_a_geo": np.float64,
     "drain_b_geo": np.float64,
     "drain_min_x": np.float64,
     "drain_min_x_act": np.float64,
     "drain_min_x_nuc_homo": np.float64,
     "drain_min_x_nuc_hetero": np.float64,
     "drain_min_x_melt": np.float64,
     "drain_min_x_evap": np.float64,
     "drain_min_x_freezing": np.float64,
     "drain_min_x_depo": np.float64,
     "drain_min_x_collision": np.float64,
     "drain_min_x_collection": np.float64,
     "drain_min_x_conversion": np.float64,
     "drain_min_x_sedimentation": np.float64,
     "drain_min_x_riming": np.float64,
     "drain_max_x": np.float64,
     "drain_sc_theta_q": np.float64,
     "drain_sc_delta_q": np.float64,
     "drain_sc_theta_n": np.float64,
     "drain_sc_delta_n": np.float64,
     "drain_s_vel": np.float64,
     "drain_a_vel": np.float64,
     "drain_b_vel": np.float64,
     "drain_rho_v": np.float64,
     "drain_c_z": np.float64,
     "drain_sc_coll_n": np.float64,
     "drain_cmu0": np.float64,
     "drain_cmu1": np.float64,
     "drain_cmu2": np.float64,
     "drain_cmu3": np.float64,
     "drain_cmu4": np.float64,
     "drain_cmu5": np.float64,
     "drain_alpha": np.float64,
     "drain_beta": np.float64,
     "drain_gamma": np.float64,
     "drain_nu": np.float64,
     "drain_g1": np.float64,
     "drain_g2": np.float64,
     "drain_mu": np.float64,
     "drain_nm1": np.float64,
     "drain_nm2": np.float64,
     "drain_nm3": np.float64,
     "drain_q_crit_c": np.float64,
     "drain_d_crit_c": np.float64,
     "drain_ecoll_c": np.float64,
     "drain_cap": np.float64,
     "drain_a_ven": np.float64,
     "drain_b_ven": np.float64,
     "drain_c_s": np.float64,
     "drain_a_f": np.float64,
     "drain_b_f": np.float64,
     "drain_alfa_n": np.float64,
     "drain_alfa_q": np.float64,
     "drain_lambda": np.float64,
     "drain_vsedi_min": np.float64,
     "drain_vsedi_max": np.float64,
     "dcloud_a_geo": np.float64,
     "dcloud_b_geo": np.float64,
     "dcloud_min_x": np.float64,
     "dcloud_min_x_act": np.float64,
     "dcloud_min_x_nuc_homo": np.float64,
     "dcloud_min_x_nuc_hetero": np.float64,
     "dcloud_min_x_melt": np.float64,
     "dcloud_min_x_evap": np.float64,
     "dcloud_min_x_freezing": np.float64,
     "dcloud_min_x_depo": np.float64,
     "dcloud_min_x_collision": np.float64,
     "dcloud_min_x_collection": np.float64,
     "dcloud_min_x_conversion": np.float64,
     "dcloud_min_x_sedimentation": np.float64,
     "dcloud_min_x_riming": np.float64,
     "dcloud_max_x": np.float64,
     "dcloud_sc_theta_q": np.float64,
     "dcloud_sc_delta_q": np.float64,
     "dcloud_sc_theta_n": np.float64,
     "dcloud_sc_delta_n": np.float64,
     "dcloud_s_vel": np.float64,
     "dcloud_a_vel": np.float64,
     "dcloud_b_vel": np.float64,
     "dcloud_rho_v": np.float64,
     "dcloud_c_z": np.float64,
     "dcloud_sc_coll_n": np.float64,
     "dcloud_cmu0": np.float64,
     "dcloud_cmu1": np.float64,
     "dcloud_cmu2": np.float64,
     "dcloud_cmu3": np.float64,
     "dcloud_cmu4": np.float64,
     "dcloud_cmu5": np.float64,
     "dcloud_alpha": np.float64,
     "dcloud_beta": np.float64,
     "dcloud_gamma": np.float64,
     "dcloud_nu": np.float64,
     "dcloud_g1": np.float64,
     "dcloud_g2": np.float64,
     "dcloud_mu": np.float64,
     "dcloud_nm1": np.float64,
     "dcloud_nm2": np.float64,
     "dcloud_nm3": np.float64,
     "dcloud_q_crit_c": np.float64,
     "dcloud_d_crit_c": np.float64,
     "dcloud_ecoll_c": np.float64,
     "dcloud_cap": np.float64,
     "dcloud_a_ven": np.float64,
     "dcloud_b_ven": np.float64,
     "dcloud_c_s": np.float64,
     "dcloud_a_f": np.float64,
     "dcloud_b_f": np.float64,
     "dcloud_alfa_n": np.float64,
     "dcloud_alfa_q": np.float64,
     "dcloud_lambda": np.float64,
     "dcloud_vsedi_min": np.float64,
     "dcloud_vsedi_max": np.float64,
     "dgraupel_a_geo": np.float64,
     "dgraupel_b_geo": np.float64,
     "dgraupel_min_x": np.float64,
     "dgraupel_min_x_act": np.float64,
     "dgraupel_min_x_nuc_homo": np.float64,
     "dgraupel_min_x_nuc_hetero": np.float64,
     "dgraupel_min_x_melt": np.float64,
     "dgraupel_min_x_evap": np.float64,
     "dgraupel_min_x_freezing": np.float64,
     "dgraupel_min_x_depo": np.float64,
     "dgraupel_min_x_collision": np.float64,
     "dgraupel_min_x_collection": np.float64,
     "dgraupel_min_x_conversion": np.float64,
     "dgraupel_min_x_sedimentation": np.float64,
     "dgraupel_min_x_riming": np.float64,
     "dgraupel_max_x": np.float64,
     "dgraupel_sc_theta_q": np.float64,
     "dgraupel_sc_delta_q": np.float64,
     "dgraupel_sc_theta_n": np.float64,
     "dgraupel_sc_delta_n": np.float64,
     "dgraupel_s_vel": np.float64,
     "dgraupel_a_vel": np.float64,
     "dgraupel_b_vel": np.float64,
     "dgraupel_rho_v": np.float64,
     "dgraupel_c_z": np.float64,
     "dgraupel_sc_coll_n": np.float64,
     "dgraupel_cmu0": np.float64,
     "dgraupel_cmu1": np.float64,
     "dgraupel_cmu2": np.float64,
     "dgraupel_cmu3": np.float64,
     "dgraupel_cmu4": np.float64,
     "dgraupel_cmu5": np.float64,
     "dgraupel_alpha": np.float64,
     "dgraupel_beta": np.float64,
     "dgraupel_gamma": np.float64,
     "dgraupel_nu": np.float64,
     "dgraupel_g1": np.float64,
     "dgraupel_g2": np.float64,
     "dgraupel_mu": np.float64,
     "dgraupel_nm1": np.float64,
     "dgraupel_nm2": np.float64,
     "dgraupel_nm3": np.float64,
     "dgraupel_q_crit_c": np.float64,
     "dgraupel_d_crit_c": np.float64,
     "dgraupel_ecoll_c": np.float64,
     "dgraupel_cap": np.float64,
     "dgraupel_a_ven": np.float64,
     "dgraupel_b_ven": np.float64,
     "dgraupel_c_s": np.float64,
     "dgraupel_a_f": np.float64,
     "dgraupel_b_f": np.float64,
     "dgraupel_alfa_n": np.float64,
     "dgraupel_alfa_q": np.float64,
     "dgraupel_lambda": np.float64,
     "dgraupel_vsedi_min": np.float64,
     "dgraupel_vsedi_max": np.float64,
     "dhail_a_geo": np.float64,
     "dhail_b_geo": np.float64,
     "dhail_min_x": np.float64,
     "dhail_min_x_act": np.float64,
     "dhail_min_x_nuc_homo": np.float64,
     "dhail_min_x_nuc_hetero": np.float64,
     "dhail_min_x_melt": np.float64,
     "dhail_min_x_evap": np.float64,
     "dhail_min_x_freezing": np.float64,
     "dhail_min_x_depo": np.float64,
     "dhail_min_x_collision": np.float64,
     "dhail_min_x_collection": np.float64,
     "dhail_min_x_conversion": np.float64,
     "dhail_min_x_sedimentation": np.float64,
     "dhail_min_x_riming": np.float64,
     "dhail_max_x": np.float64,
     "dhail_sc_theta_q": np.float64,
     "dhail_sc_delta_q": np.float64,
     "dhail_sc_theta_n": np.float64,
     "dhail_sc_delta_n": np.float64,
     "dhail_s_vel": np.float64,
     "dhail_a_vel": np.float64,
     "dhail_b_vel": np.float64,
     "dhail_rho_v": np.float64,
     "dhail_c_z": np.float64,
     "dhail_sc_coll_n": np.float64,
     "dhail_cmu0": np.float64,
     "dhail_cmu1": np.float64,
     "dhail_cmu2": np.float64,
     "dhail_cmu3": np.float64,
     "dhail_cmu4": np.float64,
     "dhail_cmu5": np.float64,
     "dhail_alpha": np.float64,
     "dhail_beta": np.float64,
     "dhail_gamma": np.float64,
     "dhail_nu": np.float64,
     "dhail_g1": np.float64,
     "dhail_g2": np.float64,
     "dhail_mu": np.float64,
     "dhail_nm1": np.float64,
     "dhail_nm2": np.float64,
     "dhail_nm3": np.float64,
     "dhail_q_crit_c": np.float64,
     "dhail_d_crit_c": np.float64,
     "dhail_ecoll_c": np.float64,
     "dhail_cap": np.float64,
     "dhail_a_ven": np.float64,
     "dhail_b_ven": np.float64,
     "dhail_c_s": np.float64,
     "dhail_a_f": np.float64,
     "dhail_b_f": np.float64,
     "dhail_alfa_n": np.float64,
     "dhail_alfa_q": np.float64,
     "dhail_lambda": np.float64,
     "dhail_vsedi_min": np.float64,
     "dhail_vsedi_max": np.float64,
     "dice_a_geo": np.float64,
     "dice_b_geo": np.float64,
     "dice_min_x": np.float64,
     "dice_min_x_act": np.float64,
     "dice_min_x_nuc_homo": np.float64,
     "dice_min_x_nuc_hetero": np.float64,
     "dice_min_x_melt": np.float64,
     "dice_min_x_evap": np.float64,
     "dice_min_x_freezing": np.float64,
     "dice_min_x_depo": np.float64,
     "dice_min_x_collision": np.float64,
     "dice_min_x_collection": np.float64,
     "dice_min_x_conversion": np.float64,
     "dice_min_x_sedimentation": np.float64,
     "dice_min_x_riming": np.float64,
     "dice_max_x": np.float64,
     "dice_sc_theta_q": np.float64,
     "dice_sc_delta_q": np.float64,
     "dice_sc_theta_n": np.float64,
     "dice_sc_delta_n": np.float64,
     "dice_s_vel": np.float64,
     "dice_a_vel": np.float64,
     "dice_b_vel": np.float64,
     "dice_rho_v": np.float64,
     "dice_c_z": np.float64,
     "dice_sc_coll_n": np.float64,
     "dice_cmu0": np.float64,
     "dice_cmu1": np.float64,
     "dice_cmu2": np.float64,
     "dice_cmu3": np.float64,
     "dice_cmu4": np.float64,
     "dice_cmu5": np.float64,
     "dice_alpha": np.float64,
     "dice_beta": np.float64,
     "dice_gamma": np.float64,
     "dice_nu": np.float64,
     "dice_g1": np.float64,
     "dice_g2": np.float64,
     "dice_mu": np.float64,
     "dice_nm1": np.float64,
     "dice_nm2": np.float64,
     "dice_nm3": np.float64,
     "dice_q_crit_c": np.float64,
     "dice_d_crit_c": np.float64,
     "dice_ecoll_c": np.float64,
     "dice_cap": np.float64,
     "dice_a_ven": np.float64,
     "dice_b_ven": np.float64,
     "dice_c_s": np.float64,
     "dice_a_f": np.float64,
     "dice_b_f": np.float64,
     "dice_alfa_n": np.float64,
     "dice_alfa_q": np.float64,
     "dice_lambda": np.float64,
     "dice_vsedi_min": np.float64,
     "dice_vsedi_max": np.float64,
     "dsnow_a_geo": np.float64,
     "dsnow_b_geo": np.float64,
     "dsnow_min_x": np.float64,
     "dsnow_min_x_act": np.float64,
     "dsnow_min_x_nuc_homo": np.float64,
     "dsnow_min_x_nuc_hetero": np.float64,
     "dsnow_min_x_melt": np.float64,
     "dsnow_min_x_evap": np.float64,
     "dsnow_min_x_freezing": np.float64,
     "dsnow_min_x_depo": np.float64,
     "dsnow_min_x_collision": np.float64,
     "dsnow_min_x_collection": np.float64,
     "dsnow_min_x_conversion": np.float64,
     "dsnow_min_x_sedimentation": np.float64,
     "dsnow_min_x_riming": np.float64,
     "dsnow_max_x": np.float64,
     "dsnow_sc_theta_q": np.float64,
     "dsnow_sc_delta_q": np.float64,
     "dsnow_sc_theta_n": np.float64,
     "dsnow_sc_delta_n": np.float64,
     "dsnow_s_vel": np.float64,
     "dsnow_a_vel": np.float64,
     "dsnow_b_vel": np.float64,
     "dsnow_rho_v": np.float64,
     "dsnow_c_z": np.float64,
     "dsnow_sc_coll_n": np.float64,
     "dsnow_cmu0": np.float64,
     "dsnow_cmu1": np.float64,
     "dsnow_cmu2": np.float64,
     "dsnow_cmu3": np.float64,
     "dsnow_cmu4": np.float64,
     "dsnow_cmu5": np.float64,
     "dsnow_alpha": np.float64,
     "dsnow_beta": np.float64,
     "dsnow_gamma": np.float64,
     "dsnow_nu": np.float64,
     "dsnow_g1": np.float64,
     "dsnow_g2": np.float64,
     "dsnow_mu": np.float64,
     "dsnow_nm1": np.float64,
     "dsnow_nm2": np.float64,
     "dsnow_nm3": np.float64,
     "dsnow_q_crit_c": np.float64,
     "dsnow_d_crit_c": np.float64,
     "dsnow_ecoll_c": np.float64,
     "dsnow_cap": np.float64,
     "dsnow_a_ven": np.float64,
     "dsnow_b_ven": np.float64,
     "dsnow_c_s": np.float64,
     "dsnow_a_f": np.float64,
     "dsnow_b_f": np.float64,
     "dsnow_alfa_n": np.float64,
     "dsnow_alfa_q": np.float64,
     "dsnow_lambda": np.float64,
     "dsnow_vsedi_min": np.float64,
     "dsnow_vsedi_max": np.float64
}


def parse_attr(attr):
    """
    Parse a file to global attributes dictionary and a dictionary for every
    column

    Parameters
    ----------
    attr : Path to file in ini format.

    Returns
    -------
    Dictionary of dictionaries
        {"Global attributes": {}, "Non global attributes": {"Column name": {}}}
    """
    attributes = {}
    key = None
    column_key = None
    with open(attr, "r") as opened:
        for line in opened:
            if "[" in line and "]" in line:
                key = line.replace("[", "").replace("]", "").rstrip()
                attributes[key] = {}
            elif key == "Global attributes":
                if "name" in line:
                    gl_name = line.split("=")[1].rstrip()
                elif "values" in line:
                    attributes[key][gl_name] = line.split("=")[1].rstrip()
            else:
                if "column" in line:
                    column_key = line.split("=")[1].rstrip()
                    if "_IN" in column_key or "Q_TURBULENCE" == column_key or "dtype" == column_key:
                        continue
                    attributes[key][column_key] = {}
                else:
                    if "_IN" in column_key or "Q_TURBULENCE" == column_key or "dtype" == column_key:
                        continue
                    val = line.split("=")[1].rstrip()
                    if val == "bool":
                        continue
                    attributes[key][column_key][line.split("=")[0].rstrip()] = val
    return attributes


def filter_zeros(df_dict, EPSILON=1e-31):
    """
    Drop all columns that have zero impact

    Parameters
    ----------
    df_dict : Dictionary of pandas.Dataframe
        A dictionary of pandas.Dataframe with key the output parameter
        and the columns the input parameters and values are the
        derivatives.
    EPSILON : float
        If all values in a column are lower than EPSILON,
        drop that one.

    Returns
    -------
    Dictionary of pandas.Dataframe
        Modified df_dict with removed columns in each dataframe where only (near)
        zeros existed before.
    """
    key_drop = []
    ignore_keys = ["timestep", "trajectory", "LONGITUDE", "LATITUDE", "MAP"]
    for key in df_dict:
        to_drop = []
        for column in df_dict[key]:
            if column in ignore_keys:
                continue
            if not (abs(df_dict[key][column]) > abs(EPSILON)).any():
                to_drop.append(column)
        if len(to_drop) > 0:
            df_dict[key] = df_dict[key].drop(columns=to_drop)
            print("Dropped {} columns for {}. Shape: {}".format(
                len(to_drop), key, np.shape(df_dict[key])))
            if (df_dict[key].empty or
                (len(df_dict[key].columns) == 1
                 and df_dict[key].columns[0] == "timestep")):
                key_drop.append(key)

    for key in key_drop:
        print("Dropping {} entirely.".format(key))
        del df_dict[key]
    return df_dict


def filter_high(df_dict, high=1e1):
    """
    Drop all columns that have a suspiciously high impact.

    Parameters
    ----------
    df_dict : Dictionary of pandas.Dataframe
        A dictionary of pandas.Dataframe with key the output parameter
        and the columns the input parameters and values are the
        derivatives.
    high : float
        If some values in a column are higher than high,
        drop that one.

    Returns
    -------
    Dictionary of pandas.Dataframe
        Modified df_dict with removed columns in each dataframe where some values
        were too high before.
    """
    key_drop = []
    for key in df_dict:
        to_drop = []
        for column in df_dict[key]:
            if column == "timestep" or column == "trajectory":
                continue
            if (abs(df_dict[key][column]) >= abs(high)).any():
                to_drop.append(column)
        if len(to_drop) > 0:
            df_dict[key] = df_dict[key].drop(columns=to_drop)
            print("Dropped {} columns for {} (too high values). Shape: {}".format(
                len(to_drop), key, np.shape(df_dict[key])))
            if (df_dict[key].empty or
                (len(df_dict[key].columns) == 1
                 and df_dict[key].columns[0] == "timestep")):

                print("Dropping {} entirely (too high values).".format(key))
                key_drop.append(key)

    for key in key_drop:
        del df_dict[key]
    return df_dict


def load_derivatives(prefix="", suffix="", filt=False, EPSILON=1e-31):
    """
    Load several dataframes and store them in a dictionary where
    the key is the output parameter of a model and the columns are the
    input parameters and the values are the derivatives.

    Parameters
    ----------
    prefix : string                                                                                                                           Prefix of the datafile such as "sb_ice".                                                                                          suffix : string                                                                                                                           Suffix of the datafile such as "start_over".                                                                                      filt : bool                                                                                                                               Filter the data with values smaller than EPSILON if true.                                                                         EPSILON : float                                                                                                                           If filt is true, filter values smaller than EPSILON out.

    Returns
    -------
    Dictionary with pandas.Dataframe
        Dictionary with dataframes as described above.

    """
    df_dict = {}
    for key in params_dict.keys():
        print("Loading {} from {}".format(key,
            "data/" + prefix + suffix + params_dict[key]))
        tmp = pd.read_csv(
            "data/" + prefix + suffix + params_dict[key], sep=",")

        df_dict[key] = tmp
    if filt:
        df_dict = filter_zeros(df_dict, EPSILON)
    for key in df_dict.keys():

        if key == "qr":
        #     print("Renaming stuff in")
        #     print(df_dict[key].keys())
            df_dict[key] = df_dict[key].rename(mapper={"dinv_z": r"$z_{inv}$",
                                 "drain_a_geo": r"$a_{rain, geo}$",
                                 "drain_b_geo": r"$b_{rain, geo}$",
                                 "drain_alpha": r"$\alpha_{rain}$",
                                 "drain_nu": r"$\nu_{rain}$"}, axis="columns")
                # df_dict[key] = df_dict[key].rename({"dinv_z": "inv_z",
                #                  "drain_a_geo": "rain_a_geo",
                #                  "drain_b_geo": "rain_b_geo",
                #                  "drain_alpha": "rain_alpha",
                #                  "drain_nu": "rain_nu"})
            print(df_dict[key])
        df_dict[key] = transform_df(df_dict[key])
    return df_dict


def load_nc(inp="/mnt/localscratch/data/project/m2_jgu-tapt/o"
                + "nline_trajectories/foehn201305_case/" +
                "foehn201305_warming1.nc", get_pol=False):
    """
    Read a netCDF file and create a pandas.Dataframe.

    Parameters
    ----------
    inp : string
        Path to a nc file.
    get_pol : bool
        If true: return latitude and longitude of the pole and mark NaN rows
        with time -1.

    Returns
    -------
    pandas.Dataframe, [float, float]
        The loaded file as pandas.Dataframe. If get_pol=True: Add lat and lon
        of pole.
    """
    ds = xr.open_dataset(inp)
    df = ds.to_dataframe()
    # Rotate coordinates
    if get_pol:
        pollat = ds.attrs["pollat"]
        pollon = ds.attrs["pollon"]
        df["ntim"] = df.index.get_level_values(0)
        df = df.loc[(df["ntim"] == 0) | (df["time"] > 0.0)]
        df.loc[df["T"] != df["T"], "time"] = -1
        return df, pollat, pollon
    return df


def load_output(filename="sb_ice.txt", sep=None, nrows=None, change_ref=True,
        refs=None, attr=None):
    """
    Read a csv file and return a pandas.Dataframe with
    physical (not normalized) entries.

    Parameters
    ----------
    filename : string
        Filepath and filename of the datafile.
    sep : string
        Separator in the file.
    nrows : int, optional
        Number of rows to read from the datafile.
    change_ref : bool
        If true: Multiply all entries with reference values.
    refs : String
        Path to a file of references for transformation
    attr : String
        Path to a file with attributes of every output column. If given, a
        MET3D styled version is assumed and triggered.

    Returns
    -------
    pandas.Dataframe
        Dataframe with certain columns if change_ref is True. Otherwise
        any dataframe from the given csv file.
        If attr is given, column names is slightly different.
    """
    type_dic = {
        "step": np.float64,
        "trajectory": np.uint64,
        "LONGITUDE": np.float64,
        "LATITUDE": np.float64,
        "MAP": np.bool_,
        "dp2h": np.bool_,
        "conv_400": np.bool_,
        "conv_600": np.bool_,
        "slan_400": np.bool_,
        "slan_600": np.bool_,
        "p": np.float64,
        "T": np.float64,
        "w": np.float64,
        "S": np.float64,
        "qc": np.float64,
        "qr": np.float64,
        "qv": np.float64,
        "Nc": np.float64,
        "Nr": np.float64,
        "qi": np.float64,
        "Ni": np.float64,
        "vi": np.float64,
        "qs": np.float64,
        "Ns": np.float64,
        "qg": np.float64,
        "Ng": np.float64,
        "qh": np.float64,
        "Nh": np.float64,
        "qiout": np.float64,
        "qsout": np.float64,
        "qrout": np.float64,
        "qgout": np.float64,
        "qhout": np.float64,
        "latent_heat": np.float64,
        "latent_cool": np.float64,
        "Niout": np.float64,
        "Nsout": np.float64,
        "Nrout": np.float64,
        "Ngout": np.float64,
        "Nhout": np.float64,
        "z": np.float64,
        "Inactive": np.float64,
        "deposition": np.float64,
        "sublimination": np.float64
    }
    if attr is not None:
        type_dic["time"] = np.float64
        type_dic["time_after_ascent"] = np.float64
        type_dic["type"] = str

    if sep is None:
        try:
            data = pd.read_csv(filename, nrows=nrows, dtype=type_dic)
        except:
            print("Reading data with pre specified data types failed.")
            print("Restarting with automatic type detection")
            data = pd.read_csv(filename, nrows=nrows)
    else:
        try:
            data = pd.read_csv(filename, sep=sep, nrows=nrows, dtype=type_dic)
        except:
            print("Reading data with pre specified data types failed.")
            print("Restarting with automatic type detection")
            data = pd.read_csv(filename, sep=sep, nrows=nrows)

    if refs is not None:
        change_ref = True
    if change_ref:
        Tref = 273.15
        pref = 100000
        qref = 1e-06
        Nref = 1
        wref = 1
        tref = 1
        zref = 1
        if refs is not None:
            refs = np.genfromtxt(refs)
            Tref = refs[0]
            pref = refs[1]
            qref = refs[2]
            Nref = refs[3]
            wref = refs[4]
            tref = refs[5]
            zref = refs[6]
        if attr is not None:
            data["time"]            = data["time"]*tref
            data["pressure"]        = data["pressure"]*pref
            data["T"]               = data["T"]*Tref
            data["w"]               = data["w"]*wref
            data["QC"]              = data["QC"]*qref
            data["QR"]              = data["QR"]*qref
            data["QS"]              = data["QS"]*qref
            data["QG"]              = data["QG"]*qref
            data["QH"]              = data["QH"]*qref
            data["QI"]              = data["QI"]*qref
            data["QV"]              = data["QV"]*qref
            data["QI_OUT"]           = data["QI_OUT"]*qref
            data["QS_OUT"]           = data["QS_OUT"]*qref
            data["QR_OUT"]           = data["QR_OUT"]*qref
            data["QG_OUT"]           = data["QG_OUT"]*qref
            data["QH_OUT"]           = data["QH_OUT"]*qref
            data["NCCLOUD"]         *= Nref
            data["NCRAIN"]          *= Nref
            data["NCSNOW"]          *= Nref
            data["NCGRAUPEL"]       *= Nref
            data["NCHAIL"]          *= Nref
            data["NCICE"]           *= Nref
            data["NI_OUT"]          *= Nref
            data["NS_OUT"]          *= Nref
            data["NR_OUT"]          *= Nref
            data["NG_OUT"]          *= Nref
            data["NH_OUT"]          *= Nref
            data["z"]               *= zref
        else:
            if "time" in data:
                data["time"]            = data["time"]*tref
            data["p"]               = data["p"]*pref
            data["T"]               = data["T"]*Tref
            data["w"]               = data["w"]*wref
            data["qc"]              = data["qc"]*qref
            data["qr"]              = data["qr"]*qref
            if "qs" in data:
                data["qs"]          = data["qs"]*qref
            if "qg" in data:
                data["qg"]          = data["qg"]*qref
            if "qh" in data:
                data["qh"]          = data["qh"]*qref
            if "qi" in data:
                data["qi"]          = data["qi"]*qref
            data["qv"]              = data["qv"]*qref
            if "qiout" in data:
                data["qiout"]       = data["qiout"]*qref
            if "qsout" in data:
                data["qsout"]       = data["qsout"]*qref
            if "qrout" in data:
                data["qrout"]       = data["qrout"]*qref
            if "qgout" in data:
                data["qgout"]       = data["qgout"]*qref
            if "qhout" in data:
                data["qhout"]       = data["qhout"]*qref
            data["Nc"]              *= Nref
            data["Nr"]              *= Nref
            if "Ns" in data:
                data["Ns"]          *= Nref
            if "Ng" in data:
                data["Ng"]          *= Nref
            if "Nh" in data:
                data["Nh"]          *= Nref
            if "Ni" in data:
                data["Ni"]          *= Nref
            if "Niout" in data:
                data["Niout"]       *= Nref
            if "Nsout" in data:
                data["Nsout"]       *= Nref
            if "Nrout" in data:
                data["Nrout"]       *= Nref
            if "Ngout" in data:
                data["Ngout"]       *= Nref
            if "Nhout" in data:
                data["Nhout"]       *= Nref
            if "z" in data:
                data["z"]           *= zref
    if attr is not None:
        # data = data.rename(columns=met3d_rename_dic)
        attributes = parse_attr(attr)
        data["QH"].attrs = {
            "standard_name": "mass_fraction_of_hail_in_air",
            "long_name": "specific hail content",
            "units": "kg kg^-1"}
        data["QH_OUT"].attrs = {
            "standard_name": "sedi_outflux_of_hail",
            "long_name": "sedimentation of hail mixing ratio",
            "units": "kg kg^-1 s^-1"}
        data["NCHAIL"].attrs = {
            "standard_name": "specif_number_of_hail_in_air",
            "long_name": "specific hail number",
            "units": "kg^-1"}
        data["NH_OUT"].attrs = {
            "standard_name": "sedi_outflux_of_hail_number",
            "long_name": "sedimentation of hail number",
            "units": "kg^-1 s^-1"}
        for key in attributes:
            if key == "Global attributes":
                data.attrs = attributes[key]
            else:
                for col in attributes[key]:
                    if col in data:
                        data[col].attrs = attributes[key][col]
    return data


def transform_df(df):
    """
    Create a new pandas.DataFrame with column "param", "timestep", "deriv"
    that can be used for plotting with seaborn.lineplot.
    Optionally adds columns LONGITUDE, LATITUDE, MAP if available in df.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe where columns are the names of the parameters.

    Returns
    -------
    pandas.Dataframe
        Transformed Dataframe.
    """
    if "MAP" in df:
        dicti = {"param": [], "timestep": [], "deriv": [], "MAP": [],
                 "LONGITUDE": [], "LATITUDE": []}
    else:
        dicti = {"param": [], "timestep": [], "deriv": []}

    key_list = ["timestep", "trajectory", "LONGITUDE", "LATITUDE", "MAP"]
    for key in df:
        if key in key_list:
            continue
        dicti["timestep"].extend(df["timestep"].tolist())
        dicti["deriv"].extend(df[key].tolist())
        dicti["param"].extend([key for i in range(len(df["timestep"]))])
        if "MAP" in df:
            dicti["MAP"].extend(df["MAP"].tolist())
            dicti["LONGITUDE"].extend(df["LONGITUDE"].tolist())
            dicti["LATITUDE"].extend(df["LATITUDE"].tolist())
    return pd.DataFrame(dicti)


def transform_df2(df, net_df, n_traj=903, traj_timestep=20):
    """
    Create a new pandas.DataFrame with column "deriv", "in_param", "out_param",
    "timestep", "trajectory", "LONGITUDE" and "LATITUDE".
    that can be used for plotting with plot_many_traj.plot_weather_deriv(..).
    It is mainly used to sum the derivatives for traj_timestep and to add
    coordinates to the data from df.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of derivatives.
    net_df : pandas.DataFrame
        Dataframe of a netCDF file, where the coordinates had been rotated.
    n_traj : int
        Number of trajectories in the netCDF dataframe.
    traj_timestep : int or float
        Timestep size in seconds of the trajectories in the netCDF dataframe.

    Returns
    -------
    pandas.Dataframe
        Transformed Dataframe.
    """
    n_rows = len(net_df.index)
    new_dic = {"deriv": [], "in_param": [], "out_param": [],
               "timestep": [], "trajectory": [],
               "LONGITUDE": [], "LATITUDE": []}
    for traj in df.trajectory.unique():
        df_traj = df.loc[df["trajectory"] == traj]
        net_df_traj = net_df.iloc[np.arange(traj, n_rows, n_traj)]
        for out_param in df_traj.out_param.unique():
            df_out = df_traj.loc[df_traj["out_param"] == out_param]
            for in_param in df_out.in_param.unique():
                df_in = df_out.loc[df_out["in_param"] == in_param]
                max_time = df_in["timestep"].max()
                for t in np.arange(traj_timestep, max_time+1, traj_timestep):
                    net_df_time = net_df_traj.loc[net_df_traj["time"] == t]
                    if net_df_time.empty:
                        continue
                    new_dic["in_param"].append(in_param)
                    new_dic["out_param"].append(out_param)
                    new_dic["timestep"].append(t)
                    summed = df_in["deriv"].sum()/traj_timestep
                    new_dic["deriv"].append(summed)
                    new_dic["LATITUDE"].append(net_df_time["lat"][0])
                    new_dic["LONGITUDE"].append(net_df_time["lon"][0])
                    new_dic["trajectory"].append(traj)
    return pd.DataFrame.from_dict(new_dic)


def load_mult_derivates(prefix="", suffix="", filt=False, EPSILON=1e-31,
                        lo=0, hi=0, high=None):
    """
    Create a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param

    Parameters
    ----------
    prefix : string
        Prefix of the datafile such as "sb_ice".
    suffix : string
        Suffix of the datafile such as "start_over".
    filt : bool
        Filter the data with values smaller than EPSILON if true.
    EPSILON : float
        If filt is true, filter values smaller than EPSILON out.
    lo : int
        Lowest trajectory number to read in.
    hi : int
        Highest trajectory number to read in.
    high : float
        If not None: Filter those derivatives out that have
        higher values.

    Returns
    -------
    pandas.Dataframe
        Pandas dataframe as described above.
    """
    df_dict = {"timestep": [], "trajectory": [], "out_param": [],
               "in_param": [], "deriv": []}

    for i in pb(range(lo, hi+1), redirect_stdout=True):
        tmp_dict = {}
        try:
            for out_param in params_dict.keys():
                tmp = pd.read_csv(prefix + str(i) + "_"
                                + suffix + params_dict[out_param], sep=",")
                tmp_dict[out_param] = tmp
            if filt:
                tmp_dict = filter_zeros(tmp_dict, EPSILON)
            if not high is None:
                tmp_dict = filter_high(tmp_dict, high)
            for out_param in tmp_dict.keys():
                tmp_dict[out_param] = transform_df(tmp_dict[out_param])
            for out_param in tmp_dict.keys():
                n_entries = len(tmp_dict[out_param].index)
                dic = tmp_dict[out_param]
                df_dict["trajectory"].extend([i for j in range(n_entries)])
                df_dict["out_param"].extend([out_param for j in range(n_entries)])
                df_dict["timestep"].extend(dic["timestep"])
                df_dict["in_param"].extend(dic["param"])
                df_dict["deriv"].extend(dic["deriv"])

        except:
            pass
    return pd.DataFrame(df_dict)


def load_mult_derivates_big(prefix="", suffix="", filt=False, EPSILON=1e-31,
                        lo=0, hi=0, high=None, trajectories=None):
    """
    Create a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param

    Parameters
    ----------
    prefix : string
        Prefix of the datafile such as "sb_ice".
    suffix : string
        Suffix of the datafile such as "start_over".
    filt : bool                                                                                                                               Filter the data with values smaller than EPSILON if true.                                                                         EPSILON : float                                                                                                                           If filt is true, filter values smaller than EPSILON out.
    lo : int
        Lowest trajectory number to read in.
    hi : int
        Highest trajectory number to read in.
    high : float
        If not None: Filter those derivatives out that have
        higher values.
    trajectories : list of int
        If given: Load only the trajectories from this list.

    Returns
    -------
    pandas.Dataframe
        Pandas dataframe as described above.
    """
    df = pd.DataFrame(data={"timestep": [], "trajectory": [], "out_param": [],
                            "in_param": [], "deriv": []})
    if trajectories is None:
        trajectories = range(lo, hi+1)
    for i in pb(trajectories, redirect_stdout=True):
        tmp_dict = {}
        try:
            for out_param in params_dict.keys():
                tmp = pd.read_csv(prefix + str(i) + "_"
                                + suffix + params_dict[out_param], sep=",")
                tmp_dict[out_param] = tmp
            if filt:
                tmp_dict = filter_zeros(tmp_dict, EPSILON)
            if not high is None:
                tmp_dict = filter_high(tmp_dict, high)
            for out_param in tmp_dict.keys():
                tmp_dict[out_param] = transform_df(tmp_dict[out_param])
            for out_param in tmp_dict.keys():
                n_entries = len(tmp_dict[out_param].index)
                # df_tmp = pd.DataFrame(
                #     data={"timestep": tmp_dict[out_param]["timestep"],
                #           "trajectory": [i for j in range(n_entries)],
                #           "out_param": [out_param for j in range(n_entries)],
                #           "in_param": tmp_dict[out_param]["param"],
                #           "deriv": tmp_dict[out_param]["deriv"]})
                df = df.append(pd.DataFrame(
                    data={"timestep": tmp_dict[out_param]["timestep"],
                          "trajectory": [i for j in range(n_entries)],
                          "out_param": [out_param for j in range(n_entries)],
                          "in_param": tmp_dict[out_param]["param"],
                          "deriv": tmp_dict[out_param]["deriv"]}), ignore_index=True)

        except:
            pass
    return df


def load_parallel(f, suffix):
    """
    A helper for load_mult_derivates_direc_dic(..) to load using multiple
    processes.
    """
    out_param = params_dict2[f.split("/")[-1].split(suffix)[1]]
#     print("Loading from {}".format(f))
    return (out_param, pd.read_csv(f, sep=",", index_col=False, dtype=deriv_type_dic))
    # return (out_param, pd.read_csv(f, sep=",", index_col=False,
    #     dtype={"timestep": "double", "trajectory": "int64",
    #            "MAP": "bool", "LONGITUDE": "double", "LATITUDE": "double"}))


def load_mult_derivates_direc_dic(direc="", filt=True, file_list2=None,
                                  EPSILON=0.0, trajectories=[1], suffix="20160922_00",
                                  pool=None, met3d=True):
    """
    Create a dictionary with out parameters as keys and dictionaries with columns:
    trajectory, timestep, MAP, LATITUDE, LONGITUDE
    and a column for each in parameter such as "da_1", "da_2", "dsnow_alfa_q", ...
    Out parameters is a string such as 'p', 'T', 'w' etc
    MAP is an optionally
    available flag for interesting timesteps.

    Parameters
    ----------
    direc : string
        A path to a directory wit a list of files to read.
    filt : bool
        Filter the data with values smaller than EPSILON if true.
    EPSILON : float
        If filt is true, filter values smaller than EPSILON out.
    trajectories : list of int
        A list of trajectories to read in.
    suffix : string
        The suffix of the filenames before '_diff_xx.txt'. If none is
        given, the method tries to automatically detect it.
    pool : multiprocessing.Pool
        Used to work on data in parallel

    Returns
    -------
    dic of pandas.Dataframe
        Pandas dataframe as described above.
    """
    if file_list2 is None:
        file_list = [os.path.join(direc, f) for f in os.listdir(direc)
                     if os.path.isfile(os.path.join(direc, f))]
        file_list2 = []
        for f in file_list:
            if "diff" not in f:
                continue
            s = f.split("/")[-1]
            s = s.split("traj")
            s = s[-1].split("_")
            if int(s[0]) in trajectories:
                file_list2.append(f)


    if suffix is None:
        example = file_list[0]
        i = 1
        while "diff" in example or "reference" in example:
            example = file_list[i]
            i += 1
        # The last 4 chars should be ".txt" now.
        example = example[:-4]
        example = example.split("/")[-1]
        example = example.split("_")
        suffix = example[-2] + "_" + example[-1]
        print("Found suffix: {}".format(suffix))
    tmp_dict = {}

    def booler(x):
        return np.bool_(x)
    def inter(x):
        return np.uint32(x)
    if pool is None:
        for f in file_list2:
            if met3d:
                out_param = params_dict2_met3d[f.split("/")[-1].split(suffix)[1]]
            else:
                out_param = params_dict2[f.split("/")[-1].split(suffix)[1]]
            df = pd.read_csv(f, sep=",", index_col=False, dtype=deriv_type_dic)
            if met3d:
                df = df.rename(columns=met3d_rename_dic)
            if out_param in tmp_dict:
                tmp_dict[out_param] = tmp_dict[out_param].append(df)
            else:
                tmp_dict[out_param] = df
        if filt:
            tmp_dict = filter_zeros(tmp_dict, EPSILON)

        if "type" not in list(tmp_dict.values())[0]:
            typename = ""
            if "conv" in file_list2[0]:
                typename = "Convective "
            else:
                typename = "Slantwise "
            if "_400_" in file_list2[0]:
                typename = typename + "400hPa "
            else:
                typename = typename + "600hPa "
            if "quan25" in file_list2[0]:
                typename = typename + "25. Quantile"
            elif "quan50" in file_list2[0] or "median" in file_list2[0]:
                typename = typename + "50. Quantile"
            elif "quan75" in file_list2[0]:
                typename = typename + "75. Quantile"
            for key in tmp_dict:
                tmp_dict[key]["type"] = typename

        return tmp_dict

    else:

        for df_tuple in pb(pool.starmap(load_parallel, zip(file_list2, repeat(suffix))), redirect_stdout=True):
            out_param, df = df_tuple
            if out_param in tmp_dict:
                tmp_dict[out_param] = tmp_dict[out_param].append(df)
            else:
                tmp_dict[out_param] = df
        if filt:
            tmp_dict = filter_zeros(tmp_dict, EPSILON)

        return tmp_dict


def load_mult_derivates_directory(direc="", filt=True,
                                  EPSILON=0.0, trajectories=[1], suffix="20160922_00"):
    """
    Create a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv, MAP, LATITUDE, LONGITUDE
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param. MAP is an optionally
    available flag for interesting timesteps.

    Parameters
    ----------
    direc : string
        A path to a directory wit a list of files to read.
    filt : bool
        Filter the data with values smaller than EPSILON if true.
    EPSILON : float
        If filt is true, filter values smaller than EPSILON out.
    trajectories : list of int
        A list of trajectories to read in.
    suffix : string
        The suffix of the filenames before '_diff_xx.txt'.

    Returns
    -------
    pandas.Dataframe
        Pandas dataframe as described above.
    """
    df = pd.DataFrame(data={"timestep": [], "trajectory": [], "out_param": [],
                            "in_param": [], "deriv": [], "MAP": [],
                            "LATITUDE": [], "LONGITUDE": []})

    file_list = [os.path.join(direc, f) for f in os.listdir(direc)
                 if os.path.isfile(os.path.join(direc, f))]
    file_list2 = []
    for f in file_list:
        s = f.split("traj")
        s = s[1].split("_")
        if int(s[0]) in trajectories:
            file_list2.append(f)

    if suffix is None:
        example = file_list2[0]
        i = 1
        while "diff" in example:
            example = file_list2[i]
            i += 1
        # The last 4 chars should be ".txt" now.
        example = example[:-4]
        example = example.split("_")
        suffix = example[-2] + "_" + example[-1]

    for f in pb(file_list2, redirect_stdout=True):
        tmp_dict = {}
        try:
            out_param = params_dict2[f.split(suffix)[1]]
            tmp_dict = {out_param: pd.read_csv(f, sep=",", index_col=False)}

            s = f.split("traj")
            s = s[1].split("_")
            traj = int(s[0])
            if filt:
                tmp_dict = filter_zeros(tmp_dict, EPSILON)
            tmp_dict[out_param] = transform_df(tmp_dict[out_param])
            n_entries = len(tmp_dict[out_param].index)
            df = df.append(pd.DataFrame(
                data={"timestep": tmp_dict[out_param]["timestep"],
                        "trajectory": [traj for j in range(n_entries)],
                        "out_param": [out_param for j in range(n_entries)],
                        "in_param": tmp_dict[out_param]["param"],
                        "deriv": tmp_dict[out_param]["deriv"],
                        "MAP": tmp_dict[out_param]["MAP"],
                        "LATITUDE": tmp_dict[out_param]["LATITUDE"],
                        "LONGITUDE": tmp_dict[out_param]["LONGITUDE"]
                        }), ignore_index=True)
        except:
            pass
    return df


def rotate_df(df, pollon, pollat, lon="LONGITUDE", lat="LATITUDE"):
    """
    Rotate the longitude and latitude with the given pole coordinates.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with longitude and latitude
    pollon : float
        Longitude of pole.
    pollat : float
        Latitude of pole.
    lon : String
        "LONGITUDE" for derivative dataframe, "lon" for netCDF dataframe.
    lat : String
        "LATITUDE" for derivative dataframe, "lat" for netCDF dataframe.
    """
    # lat_v, lon_v = rotate_pole(
    #                        np.asarray(df[lon].tolist()),
    #                        np.asarray(df[lat].tolist()),
    #                        pole_lon=pollon,
    #                        pole_lat=pollat)
    # df[lon] = lon_v
    # df[lat] = lat_v
    print("Not available")


def norm_deriv(df):
    """
    Given a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv, ratio_deriv, MAP
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param.

    Normalize the ratio of the derivatives for every timestep and add that
    as another column "norm_deriv". The column "MAP" is used to normalize only
        within an area (consecutive timesteps with equal flag).

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns trajectory, timestep, out_param, in_param,
        deriv, ratio_deriv and MAP. On out: Holds another
        column "norm_deriv"
    """
    print("TODO")


def ratio_deriv(df, out_param):
    """
    Given a dataframe with columns:
    trajectory, timestep, out_param, in_param, deriv
    where out_param is a string such as 'p', 'T', 'w' etc
    in_param is a string such as "da_1", "da_2", "dsnow_alfa_q", ...
    deriv is the float derivative of the given in_param.

    Calculate the ratio of the derivatives for every timestep and add that
    as another column "ratio_deriv".

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with columns trajectory, timestep, out_param, in_param,
        deriv and optionally MAP.
    out_param : String
        Output parameter to calculate the ratio for

    Returns
    -------
    pandas.Dataframe
        Dataframe with columns trajectory, timestep, out_param, in_param,
        deriv and optionally MAP. Also additional column ratio_deriv
    """
    df_out = df[df.out_param == out_param]
    if df_out.empty:
        print("No such output parameter: {}".format(out_param))
        return None

    denominator = np.abs(df["deriv"]).max()
    df_out["ratio_deriv"] = df_out["deriv"]/denominator
    return df_out


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print("Please give a path with prefix of data and a suffix and an identifier for the csv.")
        print("Falling back to default.")
        prefix = "/data/project/wcb/sb_ice_wcb272280_traj"
        suffix = "start_over_20160922_00"
        ident  = "filt_zero_30"
    else:
        prefix = sys.argv[1]
        suffix = sys.argv[2]
        ident  = sys.argv[3]

    try:
        refs = np.loadtxt('reference_values.txt')
        Tref = refs[0]
        pref = refs[1]
        qref = refs[2]
        Nref = refs[3]
        wref = refs[4]
        tref = refs[5]
    except:
        print("No file with reference values found. Using default values.")
        Tref = 273.15
        pref = 100000
        qref = 1e-06
        Nref = 1
        wref = 1
        tref = 1
        zref = 1

    filt = True
    lo = 1
    hi = 3
    EPSILON = 1e-30
    high = None
    res = load_mult_derivates_big(prefix, suffix,
                                  filt=filt, lo=lo, hi=hi,
                                  EPSILON=EPSILON, high=high)
    print("Trying to store\n")
    saveword = "sb_ice_wcb272280" + "_" + ident + "_" + suffix + ".csv"
    res.to_csv(saveword)

