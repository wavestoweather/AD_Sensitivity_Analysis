"""Helpers to generate Latex code.

Here we store various methods to create strings to use with latex.
"""
# pylint: disable=too-many-lines
import math

import matplotlib as mpl
import numpy as np

param_id_map = [
    "pressure",
    "T",
    "w",
    "S",
    "QC",
    "QR",
    "QV",
    "NCCLOUD",
    "NCRAIN",
    "QI",
    "NCICE",
    "QS",
    "NCSNOW",
    "QG",
    "NCGRAUPEL",
    "QH",
    "NCHAIL",
    "QI_OUT",
    "QS_OUT",
    "QR_OUT",
    "QG_OUT",
    "QH_OUT",
    "latent_heat",
    "latent_cool",
    "NI_OUT",
    "NS_OUT",
    "NR_OUT",
    "NG_OUT",
    "NH_OUT",
    "z",
    "Inactive",
    "deposition",
    "sublimination",
]

physical_params = [
    "dw",
    "dEpsilon",
    "drho_ice",
    "dT_freeze",
    "dT_f",
    "drho_w",
    "dM_w",
    "dM_a",
    "dgravity_acc",
    "drho_0",
    "dK_T",
    "dR_a",
    "dT_sat_low_temp",
    "dp_sat_const_a",
    "dp_sat_const_b",
    "dp_sat_ice_const_a",
    "dp_sat_ice_const_b",
    "drain_gfak",
    "dp_sat_low_temp",
    "dcloud_nm1",
    "dcloud_nm2",
    "dcloud_nm3",
    "drain_nm1",
    "drain_nm2",
    "drain_nm3",
    "dice_nm1",
    "dice_nm2",
    "dice_nm3",
    "dsnow_nm1",
    "dsnow_nm2",
    "dsnow_nm3",
    "dhail_nm1",
    "dhail_nm2",
    "dhail_nm3",
    "dgraupel_nm1",
    "dgraupel_nm2",
    "dgraupel_nm3",
    "dL_wd",
    "dcp",
    "dD_v",
    "dL_ed",
    "dR_v",
]

json_particle_cons = [
    "a_geo",
    "b_geo",
    "min_x",
    "min_x_act",
    "min_x_nuc_homo",
    "min_x_nuc_hetero",
    "min_x_melt",
    "min_x_evap",
    "min_x_freezing",
    "min_x_depo",
    "min_x_collision",
    "min_x_collection",
    "min_x_conversion",
    "min_x_sedimentation",
    "min_x_riming",
    "max_x",
    "sc_theta_q",
    "sc_delta_q",
    "sc_theta_n",
    "sc_delta_n",
    "s_vel",
    "a_vel",
    "b_vel",
    "rho_v",
    "c_z",
    "sc_coll_n",
    "cmu0",
    "cmu1",
    "cmu2",
    "cmu3",
    "cmu4",
    "cmu5",
    "alpha",
    "beta",
    "gamma",
    "g1",
    "g2",
    "mu",
    "nu",
    "nm1",
    "nm2",
    "nm3",
    "q_crit_c",
    "d_crit_c",
    "ecoll_c",
    "cap",
    "a_ven",
    "b_ven",
    "c_s",
    "a_f",
    "b_f",
    "alfa_n",
    "alfa_q",
    "lambda",
    "vsedi_min",
    "vsedi_max",
]
json_particles = ["rain", "cloud", "graupel", "hail", "ice", "snow"]
json_model_cons = [
    "rain_gfak",
    "cloud_k_au",
    "cloud_k_sc",
    "kc_autocon",
    "inv_z",
    "dw",
    "q_crit_i",
    "D_crit_i",
    "D_conv_i",
    "q_crit_r",
    "D_crit_r",
    "q_crit_fr",
    "D_coll_c",
    "q_crit",
    "D_conv_sg",
    "D_conv_ig",
    "x_conv",
    "parcel_height",
    "alpha_spacefilling",
    "T_nuc",
    "T_freeze",
    "T_f",
    "D_eq",
    "rho_w",
    "rho_0",
    "rho_vel",
    "rho_vel_c",
    "rho_ice",
    "M_w",
    "M_a",
    "R_universal",
    "Epsilon",
    "gravity_acc",
    "R_a",
    "R_v",
    "a_v",
    "b_v",
    "a_prime",
    "b_prime",
    "c_prime",
    "K_T",
    "L_wd",
    "L_ed",
    "D_v",
    "ecoll_min",
    "ecoll_gg",
    "ecoll_gg_wet",
    "kin_visc_air",
    "C_mult",
    "T_mult_min",
    "T_mult_max",
    "T_mult_opt",
    "const0",
    "const3",
    "const4",
    "const5",
    "D_rainfrz_gh",
    "D_rainfrz_ig",
    "dv0",
    "p_sat_melt",
    "cp",
    "k_b",
    "a_HET",
    "b_HET",
    "N_sc",
    "n_f",
    "N_avo",
    "na_dust",
    "na_soot",
    "na_orga",
    "ni_het_max",
    "ni_hom_max",
    "a_dep",
    "b_dep",
    "c_dep",
    "d_dep",
    "nim_imm",
    "nin_dep",
    "alf_imm",
    "bet_dep",
    "bet_imm",
    "r_const",
    "r1_const",
    "cv",
    "p_sat_const_a",
    "p_sat_ice_const_a",
    "p_sat_const_b",
    "p_sat_ice_const_b",
    "p_sat_low_temp",
    "T_sat_low_temp",
    "alpha_depo",
    "r_0",
    "k_1_conv",
    "k_2_conv",
    "k_1_accr",
    "k_r",
    "a_ccn_1",
    "a_ccn_2",
    "a_ccn_3",
    "a_ccn_4",
    "b_ccn_1",
    "b_ccn_2",
    "b_ccn_3",
    "b_ccn_4",
    "c_ccn_1",
    "c_ccn_2",
    "c_ccn_3",
    "c_ccn_4",
    "d_ccn_1",
    "d_ccn_2",
    "d_ccn_3",
    "d_ccn_4",
]

mappings = {
    "lat_heat": "Latent Heating",
    "lat_cool": "Latent Cooling",
    "latent_heat": "Latent Heating",
    "latent_cool": "Latent Cooling",
    "dinv_z": r"$ \partial z^{-1} $",
    "ratio_deriv": "Derivative Ratio",
    "in_param": "Input Parameter",
    "p": "Pressure",
    "T": "Temperature",
    "S": "Saturation",
    "Si": "Saturation w.r.t. Ice",
    "Q_TURBULENCE": "Mass Mixing from Turbulence",
    "qv": "Water Vapor Mass Density",
    "qc": "Cloud Mass Density",
    "qr": "Rain Mass Density",
    "qs": "Snow Mass Density",
    "qi": "Ice Mass Density",
    "qg": "Graupel Mass Density",
    "qh": "Hail Mass Density",
    "Nv": "Water Vapor Particle Density",
    "Nc": "Cloud Droplet Particle Density",
    "Nr": "Rain Droplet Particle Density",
    "Ns": "Snow Particle Density",
    "Ni": "Ice Particle Density",
    "Ng": "Graupel Particle Density",
    "Nh": "Hail Particle Density",
    "qvout": "Precipitation of Water Vapor Mass Density",
    "qcout": "Precipitation of Cloud Droplet Mass Density",
    "qrout": "Precipitation of Rain Droplet Mass Density",
    "qsout": "Precipitation of Snow Mass Density",
    "qiout": "Precipitation of Ice Mass Density",
    "qgout": "Precipitation of Graupel Mass Density",
    "qhout": "Precipitation of Hail Mass Density",
    "Nrout": "Precipitation of Rain Droplets",
    "Nsout": "Precipitation of Snow Crystals",
    "Niout": "Precipitation of Ice Crystals",
    "Ngout": "Precipitation of Graupel Particles",
    "Nhout": "Precipitation of Hail Particles",
    "LATITUDE": "Latitude",
    "LONGITDUE": "Longitude",
    "pressure": "Pressure",
    "pressure_hPa": "Pressure",
    "QV": "Water Vapor Mass Density",
    "QC": "Cloud Mass Density",
    "QR": "Rain Mass Density",
    "QS": "Snow Mass Density",
    "QI": "Ice Mass Density",
    "QG": "Graupel Mass Density",
    "QH": "Hail Mass Density",
    "NCCLOUD": "Cloud Droplet Particle Density",
    "NCRAIN": "Rain Droplet Particle Density",
    "NCSNOW": "Snow Particle Density",
    "NCICE": "Ice Particle Density",
    "NCGRAUPEL": "Graupel Particle Density",
    "NCHAIL": "Hail Particle Density",
    "QR_IN": "sedimentation (from above) of rain droplet mass density",
    "QS_IN": "sedimentation (from above) of snow crystal mass density",
    "QI_IN": "sedimentation (from above) of ice crystal mass density",
    "QG_IN": "sedimentation (from above) of graupel mass density",
    "QR_OUT": "sedimentation of rain droplet mass density",
    "QS_OUT": "sedimentation of snow crystal mass density",
    "QI_OUT": "sedimentation of ice crystal mass density",
    "QG_OUT": "sedimentation of graupel mass density",
    "QH_OUT": "sedimentation of hail mass density",
    "NR_OUT": "Precipitation",
    "NS_OUT": "sedimentation of snow crystals",
    "NI_OUT": "sedimentation of ice crystals",
    "NG_OUT": "sedimentation of graupel particles",
    "NH_OUT": "sedimentation of hail particles",
    "Q_liquid": "Liquid Water Content",
    "Q_cold": "Frozen Water Content",
    "Q_total": "Total Water Content",
    "lat": "Latitude",
    "lon": "longitude",
    "relative_lon": "Relative Longitude",
    "relative_lat": "Relative Latitude",
    "z": "Height [m]",
    "w": "Ascent [m/s]",
    "MAP": "Flag for WCB-criterion",
    "Derivatives": "Derivatives",
    "timestep": "Time [s] after ascent begins",
    "time": "Time [s] after COSMO simulation begins",
    "time_after_ascent": "Time [s] after ascent begins",
    "time_after_ascent_h": "Time [h] after ascent begins",
    "step": "Simulation step",
    "dmin_x_nuc_hetero": r"$ \partial x_{\mathrm{min},\mathrm{nuc},\mathrm{het}} $",
    "dmin_x_nuc_homo": r"$ \partial x_{\mathrm{min},\mathrm{nuc},\mathrm{hom}} $",
    "dmin_x_melt": r"$ \partial x_{\mathrm{min},\mathrm{melt}} $",
    "dmin_x_evap": r"$ \partial x_{\mathrm{min},\mathrm{evap}} $",
    "dmin_x_freezing": r"$ \partial x_{\mathrm{min},\mathrm{frz}} $",
    "dmin_x_depo": r"$ \partial x_{\mathrm{min},\mathrm{dep}} $",
    "dmin_x_collision": r"$ \partial x_{\mathrm{min},\mathrm{coli}} $",
    "dmin_x_collection": r"$ \partial x_{\mathrm{min},\mathrm{coll}} $",
    "dmin_x_conversion": r"$ \partial x_{\mathrm{min},\mathrm{con}} $",
    "dmin_x_sedimentation": r"$ \partial x_{\mathrm{min},\mathrm{sed}} $",
    "dmin_x_riming": r"$ \partial x_{\mathrm{min},\mathrm{rim}} $",
    "dEpsilon": r"$ \partial \varepsilon$",
    "dkin_visc_air": r"$ \partial \nu_{\mathrm{kin}, \mathrm{air}} $",
    "dK_T": r"$K_{\mathrm{air}} $",
    "dp_sat_const_a": r"$ \partial p_{\mathrm{sat}, a} $",
    "dp_sat_ice_const_a": r"$ \partial p_{\mathrm{sat}, \mathrm{ice}, a} $",
    "dp_sat_const_b": r"$ \partial p_{\mathrm{sat}, b} $",
    "dp_sat_ice_const_b": r"$ \partial p_{\mathrm{sat}, \mathrm{ice}, b} $",
    "dp_sat_low_temp": r"$ \partial p_{\mathrm{sat}, \mathrm{low}, T} $",
    "dT_sat_low_temp": r"$ \partial T_{\mathrm{sat}, \mathrm{low}, T} $",
    "ddv0": r"$\partial D_{0, \mathrm{vapor}} $",
    "instance_id": "Instance ID",
    "dT_mult_min": r"$ \partial T_{\mathrm{mult}, \mathrm{min}} $",
    "dT_mult_max": r"$ \partial T_{\mathrm{mult}, \mathrm{max}} $",
    "dp_sat_melt": r"$ \partial p_{\mathrm{sat}, \mathrm{melt}} $",
    "drain_cmu1": r"$ \partial \mu_{\mathrm{rain}, c, 1} $",
    "drain_cmu2": r"$ \partial \mu_{\mathrm{rain}, c, 2} $",
    "drain_cmu3": r"$ \partial \mu_{\mathrm{rain}, c, 3} $",
    "drain_cmu4": r"$ \partial \mu_{\mathrm{rain}, c, 4} $",
    "Specific_Humidity": "specific humidity",
    "dD_rainfrz_ig": r"$ \partial D_{\mathrm{rainfrz}, ig} $",
    "dD_rainfrz_gh": r"$ \partial D_{\mathrm{rainfrz}, gh} $",
    "drain_g1": r"$ \partial g_{1, \mathrm{rain}} $",
    "drain_g2": r"$ \partial g_{2, \mathrm{rain}} $",
    "dgraupel_vsedi_max": r"$ \partial v_{\mathrm{graupel}, \mathrm{sedi}, \mathrm{max}} $",
    "dice_vsedi_max": r"$ \partial v_{\mathrm{ice}, \mathrm{sedi}, \mathrm{max}} $",
    "dsnow_vsedi_max": r"$ \partial v_{\mathrm{snow}, \mathrm{sedi}, \mathrm{max}} $",
    "dhail_vsedi_max": r"$ \partial v_{\mathrm{hail}, \mathrm{sedi}, \mathrm{max}} $",
    "dk_1_accr": r"$ \partial k_{\mathrm{accr}, 1} $",
    "dk_1_conv": r"$ \partial k_{\mathrm{conv}, 1} $",
    "dk_2_conv": r"$ \partial k_{\mathrm{conv}, 2} $",
    "dkc_autocon": r"$ \partial k_{c, \mathrm{autocon}} $",
    "dna_dust": r"$ \partial n_{\mathrm{dust}} $",
    "dna_orga": r"$ \partial n_{\mathrm{orga}} $",
    "dna_soot": r"$ \partial n_{\mathrm{soot}} $",
    "dni_het_max": r"$ \partial n_{\mathrm{ice}, \mathrm{het}, \mathrm{max}} $",
    "dni_hom_max": r"$ \partial n_{\mathrm{ice}, \mathrm{hom}, \mathrm{max}} $",
    "dq_crit_i": r"$ \partial Q_{\mathrm{ice}, \mathrm{crit}} $",
    "dq_crit_r": r"$ \partial Q_{\mathrm{rain}, \mathrm{crit}} $",
    "x_conv": r"$ \partial x_{\mathrm{conv}} $",
    "dcloud_crit_d": r"$ \partial d_{\mathrm{cloud}, \mathrm{crit}} $",
    "drain_crit_d": r"$ \partial d_{\mathrm{rain}, \mathrm{crit}} $",
    "dgraupel_crit_d": r"$ \partial d_{\mathrm{graupel}, \mathrm{crit}} $",
    "dhail_crit_d": r"$ \partial d_{\mathrm{hail}, \mathrm{crit}} $",
    "dsnow_crit_d": r"$ \partial d_{\mathrm{snow}, \mathrm{crit}} $",
    "dice_crit_d": r"$ \partial d_{\mathrm{ice}, \mathrm{crit}} $",
    "drain_cmu0": r"$ \partial \mu_{\mathrm{rain}, c, 0} $",
    "drain_cmu5": r"$ \partial \mu_{\mathrm{rain}, c, 5} $",
    "drain_gfak": r"$ \partial g_{\mathrm{rain}, \mathrm{fak}} $",
    "dD_crit_i": r"$ \partial D_{\mathrm{crit}, \mathrm{ice}} $",
    "dD_conv_i": r"$ \partial D_{\mathrm{conv}, \mathrm{ice}} $",
    "dD_crit_r": r"$ \partial D_{\mathrm{crit}, \mathrm{rain}} $",
    "dq_crit_fr": r"$ \partial Q_{\mathrm{crit}, \mathrm{freeze}} $",
    "dD_coll_c": r"$ \partial D_{\mathrm{crit}, \mathrm{coll}} $",
    "dq_crit": r"$ \partial Q_{\mathrm{crit}} $",
    "drho_ice": r"$ \partial \rho_{\mathrm{ice}} $",
    "drho_vel_c": r"$ \partial \rho_{\mathrm{vel}, c} $",
    "dx_conv": r"$ \partial x_{\mathrm{conv}} $",
    "dD_conv_sg": r"$ \partial D_{\mathrm{conv}, \mathrm{snow}, \mathrm{graupel}} $",
    "dD_conv_ig": r"$ \partial D_{\mathrm{conv}, \mathrm{ice}, \mathrm{graupel}} $",
    "dcloud_k_au": r"$ \partial k_{\mathrm{cloud}, \mathrm{au}} $",
    "dcloud_k_sc": r"$ \partial k_{\mathrm{cloud}, \mathrm{sc}} $",
    "dcp": r"$ \partial c_p $",
    "dconst0": r"$ \partial c_0 $",
    "dconst1": r"$ \partial c_1 $",
    "dconst2": r"$ \partial c_2 $",
    "dconst3": r"$ \partial c_3 $",
    "dconst4": r"$ \partial c_4 $",
    "dconst5": r"$ \partial c_5 $",
    "dgravity_acc": r"$ \partial g $",
    "decoll_min": r"$ \partial e_{\mathrm{coll}, \mathrm{min}} $",
    "decoll_gg": r"$ \partial e_{\mathrm{gg}} $",
    "decoll_gg_wet": r"$ \partial e_{\mathrm{gg}, \mathrm{wet}} $",
    "dhande_ccn_fac": r"$ \partial \mathrm{fac}_{\mathrm{ccn}} $",
}


# A dictionary of all the derivatives where the key refers to the
# particle for which the parameter is defined. However,
# cloud parameters can affect rain droplets!
in_params_dic = {
    "Misc": [
        "da_1",
        "da_2",
        "de_1",
        "de_2",
        "dd",
        "dN_c",
        "dgamma",
        "dbeta_c",
        "dbeta_r",
        "ddelta1",
        "ddelta2",
        "dzeta",
        "drain_gfak",
        "dcloud_k_au",
        "dcloud_k_sc",
        "dkc_autocon",
        "dinv_z",
        "dw",
        "dq_crit_i",
        "dD_crit_i",
        "dD_conv_i",
        "dq_crit_r",
        "dD_crit_r",
        "dq_crit_fr",
        "dD_coll_c",
        "dq_crit",
        "dD_conv_sg",
        "dD_conv_ig",
        "dx_conv",
        "dparcel_height",
        "dalpha_spacefilling",
        "dT_nuc",
        "dT_freeze",
        "dT_f",
        "dD_eq",
        "drho_w",
        "drho_0",
        "drho_vel",
        "drho_vel_c",
        "drho_ice",
        "dM_w",
        "dM_a",
        "dR_universal",
        "dEpsilon",
        "dgravity_acc",
        "dR_a",
        "dR_v",
        "da_v",
        "db_v",
        "da_prime",
        "db_prime",
        "dc_prime",
        "dK_T",
        "dL_wd",
        "dL_ed",
        "dD_v",
        "decoll_min",
        "decoll_gg",
        "decoll_gg_wet",
        "dkin_visc_air",
        "dC_mult",
        "dT_mult_min",
        "dT_mult_max",
        "dT_mult_opt",
        "dconst0",
        "dconst3",
        "dconst4",
        "dconst5",
        "dD_rainfrz_gh",
        "dD_rainfrz_ig",
        "ddv0",
        "dp_sat_melt",
        "dcp",
        "dk_b",
        "da_HET",
        "db_HET",
        "dN_sc",
        "dn_f",
        "dN_avo",
        "dna_dust",
        "dna_soot",
        "dna_orga",
        "dni_het_max",
        "dni_hom_max",
        "da_dep",
        "db_dep",
        "dc_dep",
        "dd_dep",
        "dnim_imm",
        "dnin_dep",
        "dalf_imm",
        "dbet_dep",
        "dbet_imm",
        "dr_const",
        "dr1_const",
        "dcv",
        "dp_sat_const_a",
        "dp_sat_ice_const_a",
        "dp_sat_const_b",
        "dp_sat_ice_const_b",
        "dp_sat_low_temp",
        "dT_sat_low_temp",
        "dalpha_depo",
        "dr_0",
        "dk_1_conv",
        "dk_2_conv",
        "dk_1_accr",
        "dk_r",
        "da_ccn_1",
        "da_ccn_2",
        "da_ccn_3",
        "da_ccn_4",
        "db_ccn_1",
        "db_ccn_2",
        "db_ccn_3",
        "db_ccn_4",
        "dc_ccn_1",
        "dc_ccn_2",
        "dc_ccn_3",
        "dc_ccn_4",
        "dd_ccn_1",
        "dd_ccn_2",
        "dd_ccn_3",
        "dd_ccn_4",
    ],
    "Rain": [
        "drain_a_geo",
        "drain_b_geo",
        "drain_min_x",
        "drain_min_x_act",
        "drain_min_x_nuc_homo",
        "drain_min_x_nuc_hetero",
        "drain_min_x_melt",
        "drain_min_x_evap",
        "drain_min_x_freezing",
        "drain_min_x_depo",
        "drain_min_x_collision",
        "drain_min_x_collection",
        "drain_min_x_conversion",
        "drain_min_x_sedimentation",
        "drain_min_x_riming",
        "drain_max_x",
        "drain_sc_theta_q",
        "drain_sc_delta_q",
        "drain_sc_theta_n",
        "drain_sc_delta_n",
        "drain_s_vel",
        "drain_a_vel",
        "drain_b_vel",
        "drain_rho_v",
        "drain_c_z",
        "drain_sc_coll_n",
        "drain_cmu0",
        "drain_cmu1",
        "drain_cmu2",
        "drain_cmu3",
        "drain_cmu4",
        "drain_cmu5",
        "drain_alpha",
        "drain_beta",
        "drain_gamma",
        "drain_nu",
        "drain_g1",
        "drain_g2",
        "drain_mu",
        "drain_nm1",
        "drain_nm2",
        "drain_nm3",
        "drain_q_crit_c",
        "drain_d_crit_c",
        "drain_ecoll_c",
        "drain_cap",
        "drain_a_ven",
        "drain_b_ven",
        "drain_c_s",
        "drain_a_f",
        "drain_b_f",
        "drain_alfa_n",
        "drain_alfa_q",
        "drain_lambda",
        "drain_vsedi_min",
        "drain_vsedi_max",
    ],
    "Cloud": [
        "dcloud_a_geo",
        "dcloud_b_geo",
        "dcloud_min_x",
        "dcloud_min_x_act",
        "dcloud_min_x_nuc_homo",
        "dcloud_min_x_nuc_hetero",
        "dcloud_min_x_melt",
        "dcloud_min_x_evap",
        "dcloud_min_x_freezing",
        "dcloud_min_x_depo",
        "dcloud_min_x_collision",
        "dcloud_min_x_collection",
        "dcloud_min_x_conversion",
        "dcloud_min_x_sedimentation",
        "dcloud_min_x_riming",
        "dcloud_max_x",
        "dcloud_sc_theta_q",
        "dcloud_sc_delta_q",
        "dcloud_sc_theta_n",
        "dcloud_sc_delta_n",
        "dcloud_s_vel",
        "dcloud_a_vel",
        "dcloud_b_vel",
        "dcloud_rho_v",
        "dcloud_c_z",
        "dcloud_sc_coll_n",
        "dcloud_cmu0",
        "dcloud_cmu1",
        "dcloud_cmu2",
        "dcloud_cmu3",
        "dcloud_cmu4",
        "dcloud_cmu5",
        "dcloud_alpha",
        "dcloud_beta",
        "dcloud_gamma",
        "dcloud_nu",
        "dcloud_g1",
        "dcloud_g2",
        "dcloud_mu",
        "dcloud_nm1",
        "dcloud_nm2",
        "dcloud_nm3",
        "dcloud_q_crit_c",
        "dcloud_d_crit_c",
        "dcloud_ecoll_c",
        "dcloud_cap",
        "dcloud_a_ven",
        "dcloud_b_ven",
        "dcloud_c_s",
        "dcloud_a_f",
        "dcloud_b_f",
        "dcloud_alfa_n",
        "dcloud_alfa_q",
        "dcloud_lambda",
        "dcloud_vsedi_min",
        "dcloud_vsedi_max",
    ],
    "Graupel": [
        "dgraupel_a_geo",
        "dgraupel_b_geo",
        "dgraupel_min_x",
        "dgraupel_min_x_act",
        "dgraupel_min_x_nuc_homo",
        "dgraupel_min_x_nuc_hetero",
        "dgraupel_min_x_melt",
        "dgraupel_min_x_evap",
        "dgraupel_min_x_freezing",
        "dgraupel_min_x_depo",
        "dgraupel_min_x_collision",
        "dgraupel_min_x_collection",
        "dgraupel_min_x_conversion",
        "dgraupel_min_x_sedimentation",
        "dgraupel_min_x_riming",
        "dgraupel_max_x",
        "dgraupel_sc_theta_q",
        "dgraupel_sc_delta_q",
        "dgraupel_sc_theta_n",
        "dgraupel_sc_delta_n",
        "dgraupel_s_vel",
        "dgraupel_a_vel",
        "dgraupel_b_vel",
        "dgraupel_rho_v",
        "dgraupel_c_z",
        "dgraupel_sc_coll_n",
        "dgraupel_cmu0",
        "dgraupel_cmu1",
        "dgraupel_cmu2",
        "dgraupel_cmu3",
        "dgraupel_cmu4",
        "dgraupel_cmu5",
        "dgraupel_alpha",
        "dgraupel_beta",
        "dgraupel_gamma",
        "dgraupel_nu",
        "dgraupel_g1",
        "dgraupel_g2",
        "dgraupel_mu",
        "dgraupel_nm1",
        "dgraupel_nm2",
        "dgraupel_nm3",
        "dgraupel_q_crit_c",
        "dgraupel_d_crit_c",
        "dgraupel_ecoll_c",
        "dgraupel_cap",
        "dgraupel_a_ven",
        "dgraupel_b_ven",
        "dgraupel_c_s",
        "dgraupel_a_f",
        "dgraupel_b_f",
        "dgraupel_alfa_n",
        "dgraupel_alfa_q",
        "dgraupel_lambda",
        "dgraupel_vsedi_min",
        "dgraupel_vsedi_max",
    ],
    "Hail": [
        "dhail_a_geo",
        "dhail_b_geo",
        "dhail_min_x",
        "dhail_min_x_act",
        "dhail_min_x_nuc_homo",
        "dhail_min_x_nuc_hetero",
        "dhail_min_x_melt",
        "dhail_min_x_evap",
        "dhail_min_x_freezing",
        "dhail_min_x_depo",
        "dhail_min_x_collision",
        "dhail_min_x_collection",
        "dhail_min_x_conversion",
        "dhail_min_x_sedimentation",
        "dhail_min_x_riming",
        "dhail_max_x",
        "dhail_sc_theta_q",
        "dhail_sc_delta_q",
        "dhail_sc_theta_n",
        "dhail_sc_delta_n",
        "dhail_s_vel",
        "dhail_a_vel",
        "dhail_b_vel",
        "dhail_rho_v",
        "dhail_c_z",
        "dhail_sc_coll_n",
        "dhail_cmu0",
        "dhail_cmu1",
        "dhail_cmu2",
        "dhail_cmu3",
        "dhail_cmu4",
        "dhail_cmu5",
        "dhail_alpha",
        "dhail_beta",
        "dhail_gamma",
        "dhail_nu",
        "dhail_g1",
        "dhail_g2",
        "dhail_mu",
        "dhail_nm1",
        "dhail_nm2",
        "dhail_nm3",
        "dhail_q_crit_c",
        "dhail_d_crit_c",
        "dhail_ecoll_c",
        "dhail_cap",
        "dhail_a_ven",
        "dhail_b_ven",
        "dhail_c_s",
        "dhail_a_f",
        "dhail_b_f",
        "dhail_alfa_n",
        "dhail_alfa_q",
        "dhail_lambda",
        "dhail_vsedi_min",
        "dhail_vsedi_max",
    ],
    "Ice": [
        "dice_a_geo",
        "dice_b_geo",
        "dice_min_x",
        "dice_min_x_act",
        "dice_min_x_nuc_homo",
        "dice_min_x_nuc_hetero",
        "dice_min_x_melt",
        "dice_min_x_evap",
        "dice_min_x_freezing",
        "dice_min_x_depo",
        "dice_min_x_collision",
        "dice_min_x_collection",
        "dice_min_x_conversion",
        "dice_min_x_sedimentation",
        "dice_min_x_riming",
        "dice_max_x",
        "dice_sc_theta_q",
        "dice_sc_delta_q",
        "dice_sc_theta_n",
        "dice_sc_delta_n",
        "dice_s_vel",
        "dice_a_vel",
        "dice_b_vel",
        "dice_rho_v",
        "dice_c_z",
        "dice_sc_coll_n",
        "dice_cmu0",
        "dice_cmu1",
        "dice_cmu2",
        "dice_cmu3",
        "dice_cmu4",
        "dice_cmu5",
        "dice_alpha",
        "dice_beta",
        "dice_gamma",
        "dice_nu",
        "dice_g1",
        "dice_g2",
        "dice_mu",
        "dice_nm1",
        "dice_nm2",
        "dice_nm3",
        "dice_q_crit_c",
        "dice_d_crit_c",
        "dice_ecoll_c",
        "dice_cap",
        "dice_a_ven",
        "dice_b_ven",
        "dice_c_s",
        "dice_a_f",
        "dice_b_f",
        "dice_alfa_n",
        "dice_alfa_q",
        "dice_lambda",
        "dice_vsedi_min",
        "dice_vsedi_max",
    ],
    "Snow": [
        "dsnow_a_geo",
        "dsnow_b_geo",
        "dsnow_min_x",
        "dsnow_min_x_act",
        "dsnow_min_x_nuc_homo",
        "dsnow_min_x_nuc_hetero",
        "dsnow_min_x_melt",
        "dsnow_min_x_evap",
        "dsnow_min_x_freezing",
        "dsnow_min_x_depo",
        "dsnow_min_x_collision",
        "dsnow_min_x_collection",
        "dsnow_min_x_conversion",
        "dsnow_min_x_sedimentation",
        "dsnow_min_x_riming",
        "dsnow_max_x",
        "dsnow_sc_theta_q",
        "dsnow_sc_delta_q",
        "dsnow_sc_theta_n",
        "dsnow_sc_delta_n",
        "dsnow_s_vel",
        "dsnow_a_vel",
        "dsnow_b_vel",
        "dsnow_rho_v",
        "dsnow_c_z",
        "dsnow_sc_coll_n",
        "dsnow_cmu0",
        "dsnow_cmu1",
        "dsnow_cmu2",
        "dsnow_cmu3",
        "dsnow_cmu4",
        "dsnow_cmu5",
        "dsnow_alpha",
        "dsnow_beta",
        "dsnow_gamma",
        "dsnow_nu",
        "dsnow_g1",
        "dsnow_g2",
        "dsnow_mu",
        "dsnow_nm1",
        "dsnow_nm2",
        "dsnow_nm3",
        "dsnow_q_crit_c",
        "dsnow_d_crit_c",
        "dsnow_ecoll_c",
        "dsnow_cap",
        "dsnow_a_ven",
        "dsnow_b_ven",
        "dsnow_c_s",
        "dsnow_a_f",
        "dsnow_b_f",
        "dsnow_alfa_n",
        "dsnow_alfa_q",
        "dsnow_lambda",
        "dsnow_vsedi_min",
        "dsnow_vsedi_max",
    ],
}

# A dictionary of physical vs non-physical parameters
in_params_phys_dic = {
    "Physical": {
        "Misc": [],
        "Rain": [],
        "Cloud": [],
        "Graupel": [],
        "Hail": [],
        "Ice": [],
        "Snow": [],
    },
    "Non-Physical": {
        "Misc": [],
        "Rain": [],
        "Cloud": [],
        "Graupel": [],
        "Hail": [],
        "Ice": [],
        "Snow": [],
    },
}

# A dictionary of descriptions for each parameter
in_params_descr_dic = {
    "da_1": "Dimensional coefficient used in one-moment warm physics for qc and qr calculation",
    "da_2": "Dimensional coefficient used in one-moment warm physics for qc and qr calculation",
    "de_1": "Dimensional coefficients used in one-moment warm physics for temperature calculation",
    "de_2": "Dimensional coefficients used in one-moment warm physics for temperature calculation",
    "dd": "Dimensional coefficient used in one-moment warm physics qr calculation for sedimentation",
    "dN_c": "Number concentration of cloud droplets needed for one-moment scheme",
    "dgamma": "Exponent used in one-moment warm physics for qc and qr calculation",
    "dbeta_c": "Exponent used in one-moment warm physics for qc and qr calculation",
    "dbeta_r": "Exponent used in one-moment warm physics for qc and qr calculation",
    "ddelta1": "Exponent used in one-moment warm physics for qv, qr, saturation and temperature calculation",
    "ddelta2": "Exponent used in one-moment warm physics for qv, qr, saturation and temperature calculation",
    "dzeta": "Exponents used in one-moment warm physics for qr calculation",
    "drain_gfak": "Coefficient for gamma evaluation in rain evaporation",
    "dcloud_k_au": "Coefficient for autoconversion of cloud to rain",
    "dcloud_k_sc": "Coefficient for autoconversion of cloud to rain",
    "dkc_autocon": "Kernel for autoconversion",
    "dinv_z": "Inverse of air parcel size (height) used in explicit sedimentation",
    "dw": "Change in buoancy",
    "dq_crit_i": "Threshold (ratio mass) for ice selfcollection",
    "dD_crit_i": "Threshold (diameter) for ice selfcollection",
    "dD_conv_i": "Threshold (diameter) for ice conversion in selfcollection",
    "dq_crit_r": "Threshold (ratio mass) for ice rain riming and snow rain riming",
    "dD_crit_r": "Threshold (diameter) for ice rain riming and snow rain riming",
    "dq_crit_fr": "Threshold (ratio mass) for rain freeze and cloud water",
    "dD_coll_c": "Upper bound for diameter in collision efficiency",
    "dq_crit": "Default threshold (ratio mass)",
    "dD_conv_sg": "Threshold (diameter) for conversion snow to graupel, ice to graupel",
    "dD_conv_ig": "Threshold (diameter) for conversion snow to graupel, ice to graupel",
    "dx_conv": "Minimum mass of conversion due to riming",
    "dparcel_height": "Height of the trajectory package",
    "dalpha_spacefilling": "Coefficient used in riming processes for "
    "enhanced melting or conversion of snow or ice to graupel",
    "dT_nuc": "Lower temperature threshold for ice nucleation",
    "dT_freeze": "Lower temperature threshold for raindrop freezing",
    "dT_f": "Lower temperature threshold for (instantaneous) raindrop freezing",
    "dD_eq": "Equilibrium diameter for Seifert & Beheng (2008), ie Eq. 20.",
    "drho_w": r"Density of liquid water in \f$ \mathrm{kg}/\mathrm{m}^3\f $",
    "drho_0": "Norm air density",
    "drho_vel": "Exponent for density correction",
    "drho_vel_c": "Exponent for density correction of cloud droplets",
    "drho_ice": r"Density of ice in $ \mathrm{kg}/\mathrm{m}^3 $",
    "dM_w": "Molar mass of water",
    "dM_a": "Molar mass of dry air",
    "dR_universal": "Universal gas constant",
    "dEpsilon": "Quotient of the individual gas constants",
    "dgravity_acc": "Gravitational acceleration",
    "dR_a": "Gas constant for dry air",
    "dR_v": "Gas constant for water vapor",
    "da_v": "Constant used in rain evaporation after Seifert (2008)",
    "db_v": "Coefficient used in rain evaporation after Seifert (2008)",
    "da_prime": "Constant used to calculate the terminal fall velocity of raindrops during rain evaporation",
    "db_prime": "Coefficient used to calculate the terminal fall velocity of raindrops during rain evaporation",
    "dc_prime": "Exponent used to calculate the terminal fall velocity of raindrops during rain evaporation",
    "dK_T": "Heat conductivity of air",
    "dL_wd": "Latent heat of evaporation of water",
    "dL_ed": " Heat of sublimination ice -> vapor",
    "dD_v": r"Diffusivity of water vapor in air at $0^\circ\mathrm{C} $",
    "decoll_min": "Min. efficiency for collisions graupel - cloud, ice - cloud, snow - cloud",
    "decoll_gg": "Collision efficiency for graupel selfcollection",
    "decoll_gg_wet": "Collision efficiency for wet graupel",
    "dkin_visc_air": "Kinematic viscosity of dry air",
    "dC_mult": "Coefficient for splintering during Hallet-Mossop ice multiplication",
    "dT_mult_min": "Coefficient used in Hallet-Mossop ice multiplication",
    "dT_mult_max": "Coefficient used in Hallet-Mossop ice multiplication",
    "dT_mult_opt": "Coefficient used in Hallet-Mossop ice multiplication",
    "dconst0": "Coefficient used in riming processes",
    "dconst3": "Coefficient used in riming processes for breaking up particles",
    "dconst4": "Coefficient used in riming processes for breaking up particles",
    "dconst5": "Coefficient used in riming processes for enhanced melting or conversion of snow or ice to graupel",
    "dD_rainfrz_gh": " Size thresholds for partitioning of freezing rain in the hail scheme",
    "dD_rainfrz_ig": "Size thresholds for partitioning of freezing rain in the hail scheme",
    "ddv0": r"Diffusivity of water vapor in air at $0^\circ\mathrm{C} $",
    "dp_sat_melt": r"Saturation pressure at $ \mathrm{T}=\mathrm{T}_\mathrm{freeze} $",
    "dcp": r"Specific heat capacity of air at constant pressure in $ \mathrm{J}/\mathrm{K}/\mathrm{kg} $",
    "dk_b": r"Boltzmann constant in $ \mathrm{J}/\mathrm{K} $",
    "da_HET": "Exponent for rain freeze with data of Barklie and Gokhale",
    "db_HET": "Coefficient for rain freeze with data of Barklie and Gokhale",
    "dN_sc": "Schmidt number",
    "dn_f": " Exponent of N_sc in the vent-coeff",
    "dN_avo": r"Avogadro number in $ \mathrm{mol}^{-1} $",
    "dna_dust": r"Initial number density of dust [$m^{-3} $]",
    "dna_soot": r"Initial number density of soot [$m^{-3} $]",
    "dna_orga": r"Initial number density of organics [$m^{-3} $]",
    "dni_het_max": r"Constants for Phillips et al. ice nucleation scheme max number "
    r"of IN between $ 1^{-10} $ per liter, i.e. 1d3-10d3",
    "dni_hom_max": r"Constants for Phillips et al. ice nucleation scheme number "
    r"of liquid aerosols between $ 100-5000 $ per liter",
    "da_dep": "Cons_idxeters for deposition formula (2) of Hande et al.",
    "db_dep": "Cons_idxeters for deposition formula (2) of Hande et al.",
    "dc_dep": "Parameters for deposition formula (2) of Hande et al.",
    "dd_dep": "Cons_idxeters for deposition formula (2) of Hande et al.",
    "dnim_imm": "Parameter for Hande et al. nucleation",
    "dnin_dep": "Parameter for Hande et al. nucleation",
    "dalf_imm": "Parameter for Hande et al. nucleation",
    "dbet_dep": "Parameter for Hande et al. nucleation",
    "dbet_imm": "Parameter for Hande et al. nucleation",
    "da_ccn_1": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "da_ccn_2": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "da_ccn_3": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "da_ccn_4": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "db_ccn_1": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "db_ccn_2": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "db_ccn_3": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "db_ccn_4": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dc_ccn_1": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dc_ccn_2": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dc_ccn_3": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dc_ccn_4": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dd_ccn_1": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dd_ccn_2": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dd_ccn_3": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dd_ccn_4": r"Parameter for calculating CCN concentration during CCN activation \citeA{hande_parameterizing_2016}",
    "dr_const": "Parameter for saturation adjustment",
    "dr1_const": "Parameter for saturation adjustment",
    "dcv": r"Specific heat capacity of water vapor at constant pressure in $ \mathrm{J}/\mathrm{K}/\mathrm{kg} $",
    "dp_sat_const_a": "Parameter for saturation adjustment. Constant saturated water vapor pressure",
    "dp_sat_ice_const_a": "Parameter for saturation adjustment. Constant saturated ice pressure",
    "dp_sat_const_b": "Parameter for saturation adjustment. Constant saturated water vapor pressure",
    "dp_sat_ice_const_b": "Parameter for saturation adjustment. Constant saturated ice pressure",
    "dp_sat_low_temp": r"Parameter for saturation adjustment. Saturated water vapor pressure at $T = 233 $K",
    "dT_sat_low_temp": "Parameter for saturation adjustment",
    "dalpha_depo": "Depostion coefficient for homogeneous ice nucleation",
    "dr_0": "Aerosol particle radius prior to freezing used in homogeneous nucleation",
    "dk_1_conv": "Exponent for autoconversion of qc to qr",
    "dk_2_conv": "Exponent for autoconversion of qc to qr",
    "dk_1_accr": "Coefficient for accretion of qc to qr",
    "dk_r": "Coefficient for accretion of qc to qr",
    # Rain
    "drain_a_geo": "Coefficient for diameter size calculation",
    "drain_b_geo": "Exponent for diameter size calculation",
    "drain_min_x": "Minimum size of the particle used in one-moment schemes",
    "drain_min_x_act": "Not used",
    "drain_min_x_nuc_homo": "Not used",
    "drain_min_x_nuc_hetero": "Not used",
    "drain_min_x_melt": "Not used",
    "drain_min_x_evap": "Minimum size of particle for evaporation",
    "drain_min_x_freezing": "Minimum size of particle for freezing",
    "drain_min_x_depo": "Not used",
    "drain_min_x_collision": "Not used",
    "drain_min_x_collection": "Minimum size of particle for different collision processes",
    "drain_min_x_conversion": "Not used",
    "drain_min_x_sedimentation": "Minimum size of particle for sedimentation",
    "drain_min_x_riming": "Minimum size of particle for riming",
    "drain_max_x": "Maximum size of particle",
    "drain_sc_theta_q": "Not used",
    "drain_sc_delta_q": "Not used",
    "drain_sc_theta_n": "Not used",
    "drain_sc_delta_n": "Not used",
    "drain_s_vel": "Not used",
    "drain_a_vel": "Coefficient for particle velocity",
    "drain_b_vel": "Exponent for particle velocity",
    "drain_rho_v": "Coefficient used in density correction for the "
    "increased terminal fall velocity with decreasing air density",
    "drain_c_z": "Coefficient for 2nd mass moment used in freezing processes",
    "drain_sc_coll_n": "Not used",
    "drain_cmu0": r"Coefficient for calculating the shape parameter $ \mu $ during rain evaporation",
    "drain_cmu1": r"Coefficient for calculating the shape parameter $ \mu $ during rain evaporation",
    "drain_cmu2": r"Coefficient for calculating the shape parameter $ \mu $ during rain evaporation",
    "drain_cmu3": r"Constant for calculating the shape parameter $ \mu $ during rain evaporation",
    "drain_cmu4": r"Constant for calculating the shape parameter $ \mu $ during rain evaporation",
    "drain_cmu5": r"Exponent for calculating the shape parameter $ \mu $ during rain evaporation",
    "drain_alpha": "Constant in rain sedimentation",
    "drain_beta": "Coefficient for rain sedimentation",
    "drain_gamma": "Exponent for rain sedimentation",
    "drain_nu": r"Parameter to calculate the shape of the generalized $ \Gamma $-distribution",
    "drain_g1": r"Right edge of incomplete gamma function,  which had been initialized with $ \mathrm{nm}_1 $",
    "drain_g2": r"Right edge of incomplete gamma function,  which had been initialized with $ \mathrm{nm}_2 $",
    "drain_mu": r"Shape parameter of the generalized $ \Gamma $-distribution",
    "drain_nm1": "Number of bins of the incomplete gamma function lookup table 1",
    "drain_nm2": "Number of bins of the incomplete gamma function lookup table 2",
    "drain_nm3": "Number of bins of the incomplete gamma function lookup table 3",
    "drain_q_crit_c": "Not used",
    "drain_d_crit_c": "Not used",
    "drain_ecoll_c": "Not used",
    "drain_cap": "Coefficient for capacity of particle",
    "drain_a_ven": "Not used",
    "drain_b_ven": "Not used",
    "drain_c_s": "Inverse of capacity. Coefficient in evaporation and vapor deposition",
    "drain_a_f": "Constant for average ventilation. Used in melting and ice-vapor processes.",
    "drain_b_f": "Coefficient for average ventilation",
    "drain_alfa_n": "Not used",
    "drain_alfa_q": "Not used",
    "drain_lambda": "Not used",
    "drain_vsedi_min": "Not used",
    "drain_vsedi_max": "Not used",
    # Cloud
    "dcloud_a_geo": "Coefficient for diameter size calculation",
    "dcloud_b_geo": "Exponent for diameter size calculation",
    "dcloud_min_x": "Minimum size of the particle",
    "dcloud_min_x_act": "Minimum size of particle for CCN activation",
    "dcloud_min_x_nuc_homo": "Not used",
    "dcloud_min_x_nuc_hetero": "Not used",
    "dcloud_min_x_melt": "Not used",
    "dcloud_min_x_evap": "Not used",
    "dcloud_min_x_freezing": "Minimum size of particle for freezing",
    "dcloud_min_x_depo": "Not used",
    "dcloud_min_x_collision": "Not used",
    "dcloud_min_x_collection": "Not used",
    "dcloud_min_x_conversion": "Minimum size of particle for conversion processes",
    "dcloud_min_x_sedimentation": "Not used",
    "dcloud_min_x_riming": "Minimum size of particle for riming",
    "dcloud_max_x": "Maximum size of particle",
    "dcloud_sc_theta_q": "Not used",
    "dcloud_sc_delta_q": "Not used",
    "dcloud_sc_theta_n": "Not used",
    "dcloud_sc_delta_n": "Not used",
    "dcloud_s_vel": "Not used",
    "dcloud_a_vel": "Coefficient for particle velocity",
    "dcloud_b_vel": "Exponent for particle velocity",
    "dcloud_rho_v": "Coefficient used in density correction for the "
    "increased terminal fall velocity with decreasing air density",
    "dcloud_c_z": "Coefficient for 2nd mass moment",
    "dcloud_sc_coll_n": "Not used",
    "dcloud_cmu0": "Not used",
    "dcloud_cmu1": "Not used",
    "dcloud_cmu2": "Not used",
    "dcloud_cmu3": "Not used",
    "dcloud_cmu4": "Not used",
    "dcloud_cmu5": "Not used",
    "dcloud_alpha": "Not used",
    "dcloud_beta": "Not used",
    "dcloud_gamma": "Not used",
    "dcloud_nu": r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution",
    "dcloud_g1": "Not used",
    "dcloud_g2": "Not used",
    "dcloud_mu": r"Shape parameter of the generalized $ \Gamma $-distribution",
    "dcloud_nm1": "Not used",
    "dcloud_nm2": "Not used",
    "dcloud_nm3": "Not used",
    "dcloud_q_crit_c": "Not used",
    "dcloud_d_crit_c": "Not used",
    "dcloud_ecoll_c": "Not used",
    "dcloud_cap": "Coefficient for capacity of particle",
    "dcloud_a_ven": "Used to calculate the constant for average ventilation. Is not tracked with AD.",
    "dcloud_b_ven": "Used to calculate the constant for average ventilation. Is not tracked with AD.",
    "dcloud_c_s": "Inverse of capacity. Coefficient in evaporation and vapor deposition",
    "dcloud_a_f": "Constant for average ventilation. Used in melting and ice-vapor processes.",
    "dcloud_b_f": "Coefficient for average ventilation",
    "dcloud_alfa_n": "Not used",
    "dcloud_alfa_q": "Not used",
    "dcloud_lambda": "Not used",
    "dcloud_vsedi_min": "Not used",
    "dcloud_vsedi_max": "Not used",
    # Graupel
    "dgraupel_a_geo": "Coefficient for diameter size calculation",
    "dgraupel_b_geo": "Exponent for diameter size calculation",
    "dgraupel_min_x": "Minimum size of the particle used in one-moment schemes",
    "dgraupel_min_x_act": "Not used",
    "dgraupel_min_x_nuc_homo": "Not used",
    "dgraupel_min_x_nuc_hetero": "Not used",
    "dgraupel_min_x_melt": "Minimum size of particle for melting",
    "dgraupel_min_x_evap": "Minimum size of particle for evaporation",
    "dgraupel_min_x_freezing": "Not used",
    "dgraupel_min_x_depo": "Minimum size of particle for vapor deposition",
    "dgraupel_min_x_collision": "Not used",
    "dgraupel_min_x_collection": "Minimum size of particle for different collision processes",
    "dgraupel_min_x_conversion": "Minimum size of particle for conversion processes",
    "dgraupel_min_x_sedimentation": "Minimum size of particle for sedimentation",
    "dgraupel_min_x_riming": "Minimum size of particle for riming",
    "dgraupel_max_x": "Maximum size of particle",
    "dgraupel_sc_theta_q": "Not used",
    "dgraupel_sc_delta_q": "Not used",
    "dgraupel_sc_theta_n": "Not used",
    "dgraupel_sc_delta_n": "Not used",
    "dgraupel_s_vel": "Variance for the assumed Gaussian velocity distributions "
    "used in collection and riming processes",
    "dgraupel_a_vel": "Coefficient for particle velocity",
    "dgraupel_b_vel": "Exponent for particle velocity",
    "dgraupel_rho_v": "Coefficient used in density correction for the "
    "increased terminal fall velocity with decreasing air density",
    "dgraupel_c_z": "Coefficient for 2nd mass moment",
    "dgraupel_sc_coll_n": "Coefficient in graupel self collection and cloud riming",
    "dgraupel_cmu0": "Not used",
    "dgraupel_cmu1": "Not used",
    "dgraupel_cmu2": "Not used",
    "dgraupel_cmu3": "Not used",
    "dgraupel_cmu4": "Not used",
    "dgraupel_cmu5": "Not used",
    "dgraupel_alpha": "Not used",
    "dgraupel_beta": "Not used",
    "dgraupel_gamma": "Not used",
    "dgraupel_nu": r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution",
    "dgraupel_g1": r"Right edge of incomplete gamma function,  which had been initialized with $ \mathrm{nm}_1 $",
    "dgraupel_g2": r"Right edge of incomplete gamma function,  which had been initialized with $ \mathrm{nm}_2 $",
    "dgraupel_mu": r"Shape parameter of the generalized $ \Gamma $-distribution",
    "dgraupel_nm1": "Number of bins of the incomplete gamma function lookup table 1",
    "dgraupel_nm2": "Number of bins of the incomplete gamma function lookup table 2",
    "dgraupel_nm3": "Not used",
    "dgraupel_q_crit_c": "Riming parameter",
    "dgraupel_d_crit_c": "Riming parameter",
    "dgraupel_ecoll_c": "Riming coefficient. Maximum collision efficiency with cloud droplets",
    "dgraupel_cap": "Coefficient for capacity of particle",
    "dgraupel_a_ven": "Vapor deposition coefficient",
    "dgraupel_b_ven": "Not used",
    "dgraupel_c_s": "Inverse of capacity. Coefficient in evaporation and vapor deposition",
    "dgraupel_a_f": "Constant for average ventilation. Used in melting and ice-vapor processes.",
    "dgraupel_b_f": "Coefficient for average ventilation",
    "dgraupel_alfa_n": "Sedimentation velocity coefficient",
    "dgraupel_alfa_q": "Sedimentation velocity coefficient",
    "dgraupel_lambda": "Sedimentation velocity coefficient",
    "dgraupel_vsedi_min": "Minimum sedimentation velocity parameter",
    "dgraupel_vsedi_max": "Maximum sedimentation velocity parameter",
    # Hail
    "dhail_a_geo": "Coefficient for diameter size calculation",
    "dhail_b_geo": "Exponent for diameter size calculation",
    "dhail_min_x": "Minimum size of the particle used in one-moment schemes",
    "dhail_min_x_act": "Not used",
    "dhail_min_x_nuc_homo": "Not used",
    "dhail_min_x_nuc_hetero": "Not used",
    "dhail_min_x_melt": "Minimum size of particle for melting",
    "dhail_min_x_evap": "Not used",
    "dhail_min_x_freezing": "Not used",
    "dhail_min_x_depo": "Minimum size of particle for vapor deposition",
    "dhail_min_x_collision": "Not used",
    "dhail_min_x_collection": "Not used",
    "dhail_min_x_conversion": "Not used",
    "dhail_min_x_sedimentation": "Minimum size of particle for sedimentation",
    "dhail_min_x_riming": "Minimum size of particle for riming",
    "dhail_max_x": "Maximum size of particle",
    "dhail_sc_theta_q": "Not used",
    "dhail_sc_delta_q": "Not used",
    "dhail_sc_theta_n": "Not used",
    "dhail_sc_delta_n": "Not used",
    "dhail_s_vel": "Variance for the assumed Gaussian velocity distributions used in collection and riming processes",
    "dhail_a_vel": "Coefficient for particle velocity",
    "dhail_b_vel": "Exponent for particle velocity",
    "dhail_rho_v": "Coefficient used in density correction for the increased "
    "terminal fall velocity with decreasing air density",
    "dhail_c_z": "Coefficient for 2nd mass moment",
    "dhail_sc_coll_n": "Coefficient in graupel self collection and cloud riming",
    "dhail_cmu0": "Not used",
    "dhail_cmu1": "Not used",
    "dhail_cmu2": "Not used",
    "dhail_cmu3": "Not used",
    "dhail_cmu4": "Not used",
    "dhail_cmu5": "Not used",
    "dhail_alpha": "Not used",
    "dhail_beta": "Not used",
    "dhail_gamma": "Not used",
    "dhail_nu": r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution",
    "dhail_g1": "Not used",
    "dhail_g2": "Not used",
    "dhail_mu": r"Shape parameter of the generalized $ \Gamma $-distribution",
    "dhail_nm1": "Not used",
    "dhail_nm2": "Not used",
    "dhail_nm3": "Not used",
    "dhail_q_crit_c": "Riming parameter",
    "dhail_d_crit_c": "Riming parameter",
    "dhail_ecoll_c": "Riming coefficient. Maximum collision efficiency with cloud droplets",
    "dhail_cap": "Coefficient for capacity of particle",
    "dhail_a_ven": "Vapor deposition coefficient",
    "dhail_b_ven": "Not used",
    "dhail_c_s": "Inverse of capacity. Coefficient in evaporation and vapor deposition",
    "dhail_a_f": "Constant for average ventilation. Used in melting and ice-vapor processes.",
    "dhail_b_f": "Coefficient for average ventilation",
    "dhail_alfa_n": "Sedimentation velocity coefficient",
    "dhail_alfa_q": "Sedimentation velocity coefficient",
    "dhail_lambda": "Sedimentation velocity coefficient",
    "dhail_vsedi_min": "Minimum sedimentation velocity parameter",
    "dhail_vsedi_max": "Maximum sedimentation velocity parameter",
    # Ice
    "dice_a_geo": "Coefficient for diameter size calculation",
    "dice_b_geo": "Exponent for diameter size calculation",
    "dice_min_x": "Minimum size of the particle used in one-moment schemes",
    "dice_min_x_act": "Minimum size of particle for ice activation",
    "dice_min_x_nuc_homo": "Minimum size of particle for homogenous nucleation",
    "dice_min_x_nuc_hetero": "Minimum size of particle for heterogeneous nucleation",
    "dice_min_x_melt": "Minimum size of particle for melting",
    "dice_min_x_evap": "Minimum size of particle for evaporation",
    "dice_min_x_freezing": "Not used",
    "dice_min_x_depo": "Minimum size of particle for vapor deposition",
    "dice_min_x_collision": "Minimum size of particle for ice-ice collision",
    "dice_min_x_collection": "Minimum size of particle for different collision processes",
    "dice_min_x_conversion": "Minimum size of particle for conversion processes",
    "dice_min_x_sedimentation": "Minimum size of particle for sedimentation",
    "dice_min_x_riming": "Minimum size of particle for riming",
    "dice_max_x": "Maximum size of particle",
    "dice_sc_theta_q": "Not used",
    "dice_sc_delta_q": "Not used",
    "dice_sc_theta_n": "Coefficient for collision particle density",
    "dice_sc_delta_n": "Coefficient for collision particle density",
    "dice_s_vel": "Variance for the assumed Gaussian velocity distributions used in collection and riming processes",
    "dice_a_vel": "Coefficient for particle velocity",
    "dice_b_vel": "Exponent for particle velocity",
    "dice_rho_v": "Coefficient used in density correction for the increased "
    "terminal fall velocity with decreasing air density",
    "dice_c_z": "Coefficient for 2nd mass moment",
    "dice_sc_coll_n": "Coefficient in graupel self collection and cloud riming",
    "dice_cmu0": "Not used",
    "dice_cmu1": "Not used",
    "dice_cmu2": "Not used",
    "dice_cmu3": "Not used",
    "dice_cmu4": "Not used",
    "dice_cmu5": "Not used",
    "dice_alpha": "Not used",
    "dice_beta": "Not used",
    "dice_gamma": "Not used",
    "dice_nu": r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution",
    "dice_g1": "Not used",
    "dice_g2": "Not used",
    "dice_mu": r"Shape parameter of the generalized $ \Gamma $-distribution",
    "dice_nm1": "Not used",
    "dice_nm2": "Not used",
    "dice_nm3": "Not used",
    "dice_q_crit_c": "Riming parameter",
    "dice_d_crit_c": "Riming parameter",
    "dice_ecoll_c": "Riming coefficient. Maximum collision efficiency with cloud droplets",
    "dice_cap": "Coefficient for capacity of particle",
    "dice_a_ven": "Vapor deposition coefficient",
    "dice_b_ven": "Not used",
    "dice_c_s": "Inverse of capacity. Coefficient in evaporation and vapor deposition",
    "dice_a_f": "Constant for average ventilation. Used in melting and ice-vapor processes.",
    "dice_b_f": "Coefficient for average ventilation",
    "dice_alfa_n": "Sedimentation velocity coefficient",
    "dice_alfa_q": "Sedimentation velocity coefficient",
    "dice_lambda": "Sedimentation velocity coefficient",
    "dice_vsedi_min": "Minimum sedimentation velocity parameter",
    "dice_vsedi_max": "Maximum sedimentation velocity parameter",
    # Snow
    "dsnow_a_geo": "Coefficient for diameter size calculation",
    "dsnow_b_geo": "Exponent for diameter size calculation",
    "dsnow_min_x": "Minimum size of the particle used in one-moment schemes",
    "dsnow_min_x_act": "Not used",
    "dsnow_min_x_nuc_homo": "Not used",
    "dsnow_min_x_nuc_hetero": "Not used",
    "dsnow_min_x_melt": "Minimum size of particle for melting",
    "dsnow_min_x_evap": "Minimum size of particle for evaporation",
    "dsnow_min_x_freezing": "Not used",
    "dsnow_min_x_depo": "Minimum size of particle for vapor deposition",
    "dsnow_min_x_collision": "Not used",
    "dsnow_min_x_collection": "Minimum size of particle for different collision processes",
    "dsnow_min_x_conversion": "Not used",
    "dsnow_min_x_sedimentation": "Minimum size of particle for sedimentation",
    "dsnow_min_x_riming": "Minimum size of particle for riming",
    "dsnow_max_x": "Maximum size of particle",
    "dsnow_sc_theta_q": "Not used",
    "dsnow_sc_delta_q": "Not used",
    "dsnow_sc_theta_n": "Coefficient for collision particle density",
    "dsnow_sc_delta_n": "Coefficient for collision particle density",
    "dsnow_s_vel": "Variance for the assumed Gaussian velocity distributions used in collection and riming processes",
    "dsnow_a_vel": "Coefficient for particle velocity",
    "dsnow_b_vel": "Exponent for particle velocity",
    "dsnow_rho_v": "Coefficient used in density correction for the increased "
    "terminal fall velocity with decreasing air density",
    "dsnow_c_z": "Coefficient for 2nd mass moment",
    "dsnow_sc_coll_n": "Coefficient in graupel self collection and cloud riming",
    "dsnow_cmu0": "Not used",
    "dsnow_cmu1": "Not used",
    "dsnow_cmu2": "Not used",
    "dsnow_cmu3": "Not used",
    "dsnow_cmu4": "Not used",
    "dsnow_cmu5": "Not used",
    "dsnow_alpha": "Not used",
    "dsnow_beta": "Not used",
    "dsnow_gamma": "Not used",
    "dsnow_nu": r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution",
    "dsnow_g1": "Not used",
    "dsnow_g2": "Not used",
    "dsnow_mu": r"Shape parameter of the generalized $ \Gamma $-distribution",
    "dsnow_nm1": "Not used",
    "dsnow_nm2": "Not used",
    "dsnow_nm3": "Not used",
    "dsnow_q_crit_c": "Riming parameter",
    "dsnow_d_crit_c": "Riming parameter",
    "dsnow_ecoll_c": "Riming coefficient. Maximum collision efficiency with cloud droplets",
    "dsnow_cap": "Coefficient for capacity of particle",
    "dsnow_a_ven": "Vapor deposition coefficient",
    "dsnow_b_ven": "Not used",
    "dsnow_c_s": "Inverse of capacity. Coefficient in evaporation and vapor deposition",
    "dsnow_a_f": "Constant for average ventilation. Used in melting and ice-vapor processes.",
    "dsnow_b_f": "Coefficient for average ventilation",
    "dsnow_alfa_n": "Sedimentation velocity coefficient",
    "dsnow_alfa_q": "Sedimentation velocity coefficient",
    "dsnow_lambda": "Sedimentation velocity coefficient",
    "dsnow_vsedi_min": "Minimum sedimentation velocity parameter",
    "dsnow_vsedi_max": "Maximum sedimentation velocity parameter",
}

# grouping input parameters (not a true grouping but rather for convenience)
in_params_grouping = {
    "artificial": [
        "drain_gfak",
        "dcloud_a_geo",
        "dcloud_b_geo",
        "dcloud_k_au",
        "dcloud_k_sc",
        "dkc_autocon",
        "dinv_z",
        "dalpha_spacefilling",
        "dD_eq",
        "drho_vel",
        "drho_vel_c",
        "da_v",
        "db_v",
        "da_prime",
        "db_prime",
        "dc_prime",
        "decoll_min",
        "dC_mult",
        "dT_mult_min",
        "dT_mult_max",
        "dT_mult_opt",
        "dconst0",
        "dconst3",
        "dconst4",
        "dconst5",
        "da_HET",
        "db_HET",
        "da_dep",
        "db_dep",
        "dc_dep",
        "dd_dep",
        "dnim_imm",
        "dnin_dep",
        "dalf_imm",
        "dbet_dep",
        "dbet_imm",
        "dr_const",
        "dr1_const",
        "dp_sat_const_a",
        "dp_sat_ice_const_a",
        "dp_sat_const_b",
        "dp_sat_ice_const_b",
        "dp_sat_low_temp",
        "dT_sat_low_temp",
        "dalpha_depo",
        "dk_1_conv",
        "dk_2_conv",
        "dk_1_accr",
        "dk_r",
        "da_ccn_1",
        "da_ccn_2",
        "da_ccn_3",
        "da_ccn_4",
        "db_ccn_1",
        "db_ccn_2",
        "db_ccn_3",
        "db_ccn_4",
        "dc_ccn_1",
        "dc_ccn_2",
        "dc_ccn_3",
        "dc_ccn_4",
        "dd_ccn_1",
        "dd_ccn_2",
        "dd_ccn_3",
        "dd_ccn_4",
        "drain_a_geo",
        "drain_b_geo",
        "drain_c_z",
        "drain_cmu0",
        "drain_cmu1",
        "drain_cmu2",
        "drain_cmu3",
        "drain_cmu4",
        "drain_cmu5",
        "drain_g1",
        "drain_g2",
        "drain_nm1",
        "drain_nm2",
        "drain_nm3",
        "drain_a_f",
        "drain_b_f",
        "dcloud_c_z",
        "dcloud_a_f",
        "dcloud_b_f",
        "dgraupel_a_geo",
        "dgraupel_b_geo",
        "dgraupel_s_vel",
        "dgraupel_c_z",
        "dgraupel_sc_coll_n",
        "dgraupel_g1",
        "dgraupel_g2",
        "dgraupel_nm1",
        "dgraupel_nm2",
        "dgraupel_ecoll_c",
        "dgraupel_a_f",
        "dgraupel_b_f",
        "dhail_a_geo",
        "dhail_b_geo",
        "dhail_s_vel",
        "dhail_c_z",
        "dhail_sc_coll_n",
        "dhail_ecoll_c",
        "dhail_a_f",
        "dhail_b_f",
        "dice_a_geo",
        "dice_b_geo",
        "dice_sc_theta_n",
        "dice_sc_delta_n",
        "dice_s_vel",
        "dice_rho_v",
        "dice_c_z",
        "dice_sc_coll_n",
        "dice_a_f",
        "dice_b_f",
        "dsnow_a_geo",
        "dsnow_b_geo",
        "dsnow_sc_theta_n",
        "dsnow_sc_delta_n",
        "dsnow_s_vel",
        "dsnow_c_z",
        "dsnow_sc_coll_n",
        "dsnow_ecoll_c",
        "dsnow_a_f",
        "dsnow_b_f",
        ####
        "dice_sc_delta_q",
        "dice_sc_delta_n",
        "dice_sc_theta_n",
        "dice_sc_theta_q",
        "dsnow_sc_delta_q",
        "dsnow_sc_delta_n",
        "dsnow_sc_theta_n",
        "dsnow_sc_theta_q",
        "drain_sc_delta_q",
        "drain_sc_delta_n",
        "drain_sc_theta_n",
        "drain_sc_theta_q",
        "dhail_sc_delta_q",
        "dhail_sc_delta_n",
        "dhail_sc_theta_n",
        "dhail_sc_theta_q",
        "dgraupel_sc_delta_q",
        "dgraupel_sc_delta_n",
        "dgraupel_sc_theta_n",
        "dgraupel_sc_theta_q",
        "dcloud_sc_delta_q",
        "dcloud_sc_delta_n",
        "dcloud_sc_theta_n",
        "dcloud_sc_theta_q",
    ],
    "artificial (threshold)": [
        "dq_crit_i",
        "dD_crit_i",
        "dD_conv_i",
        "dq_crit_r",
        "dD_crit_r",
        "dq_crit_fr",
        "dD_coll_c",
        "dq_crit",
        "dD_conv_sg",
        "dD_conv_ig",
        "dx_conv",
        "dD_rainfrz_gh",
        "dD_rainfrz_ig",
        "dni_het_max",
        "dni_hom_max",
        "drain_min_x",
        "drain_min_x_evap",
        "drain_min_x_freezing",
        "drain_min_x_collection",
        "drain_min_x_sedimentation",
        "drain_min_x_riming",
        "drain_max_x",
        "dcloud_min_x",
        "dcloud_min_x_act",
        "dcloud_min_x_freezing",
        "dcloud_min_x_conversion",
        "dcloud_min_x_riming",
        "dcloud_max_x",
        "dcloud_d_crit_c",
        "dcloud_q_crit_c",
        "dgraupel_min_x",
        "dgraupel_min_x_melt",
        "dgraupel_min_x_evap",
        "dgraupel_min_x_depo",
        "dgraupel_min_x_collection",
        "dgraupel_min_x_conversion",
        "dgraupel_min_x_sedimentation",
        "dgraupel_min_x_riming",
        "dgraupel_max_x",
        "dgraupel_q_crit_c",
        "dgraupel_d_crit_c",
        "dgraupel_vsedi_min",
        "dgraupel_vsedi_max",
        "dhail_min_x",
        "dhail_min_x_melt",
        "dhail_min_x_depo",
        "dhail_min_x_sedimentation",
        "dhail_min_x_riming",
        "dhail_max_x",
        "dhail_q_crit_c",
        "dhail_d_crit_c",
        "dhail_vsedi_min",
        "dhail_vsedi_max",
        "dice_min_x",
        "dice_min_x_act",
        "dice_min_x_nuc_homo",
        "dice_min_x_nuc_hetero",
        "dice_min_x_melt",
        "dice_min_x_evap",
        "dice_min_x_depo",
        "dice_min_x_collision",
        "dice_min_x_collection",
        "dice_min_x_conversion",
        "dice_min_x_sedimentation",
        "dice_min_x_riming",
        "dice_max_x",
        "dice_q_crit_c",
        "dice_d_crit_c",
        "dice_vsedi_min",
        "dice_vsedi_max",
        "dsnow_min_x",
        "dsnow_min_x_melt",
        "dsnow_min_x_evap",
        "dsnow_min_x_depo",
        "dsnow_min_x_collection",
        "dsnow_min_x_sedimentation",
        "dsnow_min_x_riming",
        "dsnow_max_x",
        "dsnow_q_crit_c",
        "dsnow_d_crit_c",
        "dsnow_vsedi_min",
        "dsnow_vsedi_max",
    ],
    "physical": [
        "dw",
        "dparcel_height",
        "dT_nuc",
        "dT_freeze",
        "dT_f",
        "drho_w",
        "drho_0",
        "drho_ice",
        "dM_w",
        "dM_a",
        "dR_universal",
        "dEpsilon",
        "dgravity_acc",
        "dR_a",
        "dR_v",
        "dK_T",
        "dL_wd",
        "dL_ed",
        "dD_v",
        "decoll_gg",
        "decoll_gg_wet",
        "dkin_visc_air",
        "ddv0",
        "dp_sat_melt",
        "dcp",
        "dk_b",
        "dN_sc",
        "dn_f",
        "dN_avo",
        "dcv",
        "drain_rho_v",
        "dhail_rho_v",
        "dsnow_rho_v",
    ],
    "physical (high variability)": [
        "dna_dust",
        "dna_soot",
        "dna_orga",
        "dr_0",
        "drain_a_vel",
        "drain_b_vel",
        "drain_alpha",
        "drain_beta",
        "drain_gamma",
        "drain_nu",
        "drain_mu",
        "drain_cap",
        "drain_c_s",
        "dcloud_a_vel",
        "dcloud_b_vel",
        "dcloud_rho_v",
        "dcloud_nu",
        "dcloud_mu",
        "dcloud_cap",
        "dcloud_c_s",
        "dgraupel_a_vel",
        "dgraupel_b_vel",
        "dgraupel_rho_v",
        "dgraupel_nu",
        "dgraupel_mu",
        "dgraupel_cap",
        "dgraupel_a_ven",
        "dgraupel_c_s",
        "dgraupel_alfa_n",
        "dgraupel_alfa_q",
        "dgraupel_lambda",
        "dhail_a_vel",
        "dhail_b_vel",
        "dhail_nu",
        "dhail_mu",
        "dhail_cap",
        "dhail_a_ven",
        "dhail_c_s",
        "dhail_alfa_n",
        "dhail_alfa_q",
        "dhail_lambda",
        "dice_a_vel",
        "dice_b_vel",
        "dice_nu",
        "dice_mu",
        "dice_ecoll_c",
        "dice_cap",
        "dice_a_ven",
        "dice_c_s",
        "dice_alfa_n",
        "dice_alfa_q",
        "dice_lambda",
        "dsnow_a_vel",
        "dsnow_b_vel",
        "dsnow_nu",
        "dsnow_mu",
        "dsnow_cap",
        "dsnow_a_ven",
        "dsnow_c_s",
        "dsnow_alfa_n",
        "dsnow_alfa_q",
        "dsnow_lambda",
    ],
    "1-moment": [
        "da_1",
        "da_2",
        "de_1",
        "de_2",
        "dd",
        "dN_c",
        "dgamma",
        "dbeta_c",
        "dbeta_r",
        "ddelta1",
        "ddelta2",
        "dzeta",
    ],
}


def __parse_known_suffixes(word):
    """
    Parse a single name of a derivative and return it in a latex conform type if it has
    a known suffix. Otherwise, returns the word unchanged.

    Parameters
    ----------
    word : string
        Word to parse.

    Returns
    -------
    Changed string and True if suffix is known, otherwise unchanged string and False.
    """
    changed = True
    if "_prime" in word:
        word = r"$ \partial " + word[1] + " $"
    elif "dnm" in word:
        word = r"$ \partial n_{m, " + word[-1] + r"} $"
    elif "_nm1" in word:
        word = r"$ \partial n_{\mathrm{" + word.split("_")[0][1::] + r"}, 1} $"
    elif "_nm2" in word:
        word = r"$ \partial n_{\mathrm{" + word.split("_")[0][1::] + r"}, 2} $"
    elif "_nm3" in word:
        word = r"$ \partial n_{\mathrm{" + word.split("_")[0][1::] + r"}, 3} $"
    elif "_c_s" in word[-4::]:
        word = (
            r"$ \partial \kappa_{s, " + r"\mathrm{" + word.split("_")[0][1::] + r"}} $"
        )
    elif "HET" in word:
        word = r"$ \partial " + word[1] + r"_{\mathrm{HET}} $"

    else:
        changed = False
    return word, changed


def __parse_known_words(word):
    """
    Parse a single name of a CCN or alfa derivative and return it in a latex conform type.

    Parameters
    ----------
    word : string
        Word to parse.

    Returns
    -------
    string, bool
        String to use with latex. '$' will be added if it is a CCN or alfa parameter. Bool is True if
        the string had been changed.
    """
    if "_ccn_" in word:
        return r"$ \partial " + word[1] + r"_{\mathrm{ccn}, " + word[-1] + r"} $", True
    if "_alfa_" in word:
        return (
            r"$ \partial "
            + r"\alpha_{\mathrm{"
            + word.split("_")[0][1::]
            + r"}, "
            + word[-1]
            + r"} $",
            True,
        )
    return word, False


def __replace_long_words(word):
    """
    Parse a single name of a derivative and return it in a latex conform type if it contains a long single word.

    Parameters
    ----------
    word : string
        Word to parse.

    Returns
    -------
    string
        String to use with latex with \\partial added.
    """
    long_words = {
        "freezing": "frz",
        "collection": "coll",
        "collision": "coli",
        "sedimentation": "sed",
        "conversion": "con",
        "deposition": "dep",
        "homogenous": "hom",
        "heterogenous": "het",
    }
    for l, replacement in long_words.items():
        if l in word:
            word = word.replace(l, replacement)
    # The first "d" shall be a derivative symbol
    word = r"\partial " + word[1::]

    word = word.replace("delta_", "delta ")
    # Check for this typo
    word = word.replace("alfa", "alpha")
    word = word.replace("vsedi", "sedi_v")
    return word


def __set_subscript(word):
    """
    Replace any subscripts based on hydrometeors.

    Parameters
    ----------
    word : string
        Word to parse.

    Returns
    -------
    string
        String to use with latex.
    """
    # if any of words is in there, make it to subscript
    subscript_no_math = ["snow", "graupel", "rain", "ice", "cloud", "hail", "vapor"]
    no_subscript = True
    for w in subscript_no_math:
        if w in word:
            parts = word.split(" ")
            start = parts[0]
            parts = parts[1].split("_")
            if len(parts) == 5:
                word = (
                    start
                    + " "
                    + parts[2]
                    + r"_{"
                    + r"\mathrm{"
                    + parts[0]
                    + r"},"
                    + r"\mathrm{"
                    + parts[1]
                    + r"},\mathrm{"
                    + parts[3]
                    + r"},\mathrm{"
                    + parts[4]
                    + r"}}"
                )
            elif len(parts) == 4:
                if len(parts[1]) == 1:
                    word = (
                        start
                        + " "
                        + parts[1].capitalize()
                        + r"_{"
                        + r"\mathrm{"
                        + parts[0]
                        + r"},"
                        + r"\mathrm{"
                        + parts[2]
                        + r"},"
                        + parts[3]
                        + r"}"
                    )
                elif parts[2] == "x" or parts[2] == "v":
                    word = (
                        start
                        + " "
                        + parts[2]
                        + r"_{"
                        + r"\mathrm{"
                        + parts[0]
                        + r"},"
                        + r"\mathrm{"
                        + parts[1]
                        + r"},\mathrm{"
                        + parts[3]
                        + r"}}"
                    )
                else:
                    word = (
                        start
                        + " "
                        + parts[2]
                        + r"_{"
                        + parts[1]
                        + r", \mathrm{"
                        + parts[0]
                        + r"},"
                        + parts[3]
                        + r"}"
                    )

            elif len(parts) == 3:
                word = (
                    start
                    + " "
                    + parts[2]
                    + r"_{"
                    + parts[1]
                    + r", \mathrm{"
                    + parts[0]
                    + r"}}"
                )
            else:
                word = start + " " + parts[1] + r"_{" + r"\mathrm{" + parts[0] + r"}}"
            no_subscript = False
            break
    if no_subscript:
        parts = word.split("_")
        if len(parts) == 2:
            word = parts[0] + r"_{" + parts[1] + r"}"
    return word


def __parse_word(word):
    """
    Parse a single name of a derivative and return it in a latex conform type.

    Parameters
    ----------
    word : string
        Word to parse.

    Returns
    -------
    string
        String to use with latex. If a formula has been detected, '$' will
        be added.
    """
    no_math = ["geo", "min", "max", "ven", "vel", "freeze"]
    math_keys = [
        "alpha",
        "gamma",
        "beta",
        "delta",
        "zeta",
        "rho",
        "nu",
        "mu",
        "lambda",
        "theta",
    ]

    if word in mappings:
        return mappings[word]

    word, changed = __parse_known_suffixes(word)
    if changed:
        return word
    word, changed = __parse_known_words(word)
    if changed:
        return word
    word = __replace_long_words(word)
    word = __set_subscript(word)
    word = r"$" + word + r" $"
    for w in no_math:
        word = word.replace(w, r"\mathrm{" + w + r"}")
    if r"ri\mathrm{min]g" in word:
        word = word.replace(r"ri\mathrm{min]g", r"\mathrm{riming}")
    for w in math_keys:
        word = word.replace(w, "\\" + w)
    if "nuc" in word:
        word = word.replace("\\nuc", "nuc")
    if "\\mult" in word:
        word = word.replace("\\mult", r"\mathrm{mult}")
    if r"}\text" in word:
        word = word.replace(r"}\text", r"}, \text")
    return word


def parse_word(word):
    """
    Parse a name or multiple names of a derivative and return it in a latex conform type.

    Parameters
    ----------
    word : string or list like
        Word or iterable to parse.

    Returns
    -------
    string
        String (or list of strings) to use with latex. If a formula has been detected, '$' will
        be added.
    """
    if isinstance(word, str):
        return __parse_word(word)
    return [__parse_word(w) for w in word]


in_params_numeric_value_dic = {
    "da1_scale": 0.001,
    "da2_scale": 0.0121568,
    "de1_scale": 0.0590233,
    "de2_scale": 0.185865,
    "dd_scale": 0.0147626,
    "dc_br": 0,
    "dk_br": 0,
    "dc_ccn_2": 3.22012,
    "db_ccn_4": 0.000198616,
    "db_ccn_3": 0.000184323,
    "db_ccn_2": 4.47319e-05,
    "da_ccn_3": -0.29224,
    "dk_2_conv": 0.7,
    "dk_1_conv": 400,
    "dc_ccn_1": 16.242,
    "db_ccn_1": 0.000198405,
    "dT_sat_low_temp": 273.15,
    "dp_sat_low_temp": 610.78,
    "dp_sat_ice_const_b": 7.66,
    "dD_br": 0,
    "dcv": 718.66,
    "dr1_const": 461.5,
    "dr_const": 287.04,
    "dr_0": 2.5e-07,
    "dbet_imm": 1.2044,
    "dbet_dep": 1.4705,
    "dc_ccn_3": 13.8499,
    "dalf_imm": 0.2622,
    "dnin_dep": 77167,
    "dd_dep": 0.26789,
    "da_dep": 0.27626,
    "dni_hom_max": 5e06,
    "dni_het_max": 500000,
    "dna_orga": 3e07,
    "dN_sc": 0.71,
    "db_HET": 200,
    "dcp": 1004.64,
    "dp_sat_const_b": 35.86,
    "ddv0": 2.22e-05,
    "dna_dust": 1.6e06,
    "dconst5": 0.0109087,
    "dconst4": -0.5,
    "dconst3": 0.333333,
    "dp_sat_ice_const_a": 21.8746,
    "dp_sat_const_a": 17.2694,
    "dconst0": 33333.3,
    "da_ccn_4": 2.2919e08,
    "dT_mult_opt": 268,
    "dn_f": 0.333,
    "dT_mult_max": 270,
    "dkin_visc_air": 1.5e-05,
    "dD_rainfrz_gh": 0.00125,
    "dc_prime": 600,
    "db_prime": 9.8,
    "da_prime": 9.65,
    "dd_ccn_3": 0.890749,
    "dT_mult_min": 265,
    "dD_v": 2.22e-05,
    "db_v": 0.308,
    "dd_ccn_1": 2.87736e08,
    "da_v": 0.78,
    "dR_v": 461.523,
    "dN_avo": 6.02214e23,
    "dgravity_acc": 9.80665,
    "decoll_gg_wet": 0.4,
    "dR_universal": 8.31446,
    "drho_ice": 916.7,
    "drho_vel_c": 1,
    "db_dep": 6.21,
    "da_HET": 0.65,
    "dL_ed": 2.8345e06,
    "drho_vel": 0.4,
    "drho_0": 1.225,
    "dnim_imm": 49920,
    "dc_dep": -1.3107,
    "dD_eq": 0.0011,
    "dk_r": 5.78,
    "dT_f": 233,
    "dT_freeze": 273.15,
    "dT_nuc": 268.15,
    "dalpha_spacefilling": 0.01,
    "dalpha_depo": 0.5,
    "dna_soot": 2.5e07,
    "dparcel_height": 250,
    "dM_a": 0.0289655,
    "dD_conv_sg": 0.0002,
    "dc_ccn_4": 16.2462,
    "dq_crit": 1e-09,
    "dq_crit_fr": 1e-06,
    "dD_crit_r": 0.0001,
    "decoll_gg": 0.1,
    "dK_T": 0.024,
    "dq_crit_r": 1e-05,
    "dp_sat_melt": 610.78,
    "dD_conv_i": 7.5e-05,
    "da_ccn_1": 1.83231e08,
    "dq_crit_i": 1e-06,
    "dD_crit_i": 0.0001,
    "ddw": 0,
    "dR_a": 287.047,
    "dinv_z": 0.004,
    "dC_mult": 3.5e08,
    "dkc_autocon": 9.44e09,
    "dD_rainfrz_ig": 0.0005,
    "dcloud_k_sc": 1.416e10,
    "drain_gfak": 1,
    "drho_w": 1000,
    "dzeta": 1.125,
    "ddelta2": 0.6875,
    "ddelta1": 0.5,
    "dd_ccn_2": 0.625881,
    "da_ccn_2": 0.101474,
    "dbeta_r": 0.875,
    "dx_conv": 1e-10,
    "dbeta_c": 1,
    "dd_ccn_4": 3.60849e08,
    "dL_wd": 2.5008e06,
    "dD_conv_ig": 0.0002,
    "dgamma": 1,
    "dk_1_accr": 0.0005,
    "dk_b": 1.38065e-23,
    "dcloud_k_au": 6.80769e18,
    "dN_c": 0,
    "dd": 0,
    "dM_w": 0.0180153,
    "de_2": 0,
    "dD_br_threshold": 0,
    "dEpsilon": 0.621957,
    "de_1": 0,
    "dD_coll_c": 4e-05,
    "da_2": 0,
    "decoll_min": 0.01,
    "da_1": 0,
    "dcloud_vsedi_max": 1,
    "dcloud_vsedi_min": 0,
    "dcloud_lambda": 0.5,
    "dcloud_alfa_n": 564216,
    "dcloud_a_ven": 0.78,
    "dcloud_d_crit_c": 1e-05,
    "dcloud_q_crit_c": 1e-06,
    "dcloud_b_ven": 0.308,
    "dcloud_nm3": 0,
    "dcloud_nm2": 0,
    "dcloud_nm1": 0,
    "dcloud_ecoll_c": 0,
    "dcloud_g2": 0,
    "dcloud_c_s": 0.5,
    "dcloud_g1": 0,
    "dcloud_gamma": 0,
    "dcloud_beta": 0,
    "dcloud_alpha": 0,
    "dcloud_cmu5": 0,
    "dcloud_nu": 1,
    "dcloud_cmu3": 0,
    "dcloud_cmu2": 0,
    "dcloud_cmu1": 0,
    "dcloud_b_f": 68.6733,
    "dcloud_c_z": 1.5,
    "dcloud_rho_v": 0,
    "dcloud_sc_coll_n": 0,
    "dcloud_b_vel": 0.666667,
    "dcloud_cmu0": 0,
    "dcloud_s_vel": 0,
    "dcloud_sc_delta_n": 0,
    "dcloud_sc_delta_q": 0,
    "dcloud_mu": 1,
    "dcloud_sc_theta_q": 0,
    "dcloud_alfa_q": 752288,
    "dcloud_max_x": 2.6e-10,
    "dcloud_min_x_sedimentation": 4.2e-15,
    "dcloud_a_f": 0.737109,
    "dcloud_min_x_conversion": 4.2e-15,
    "dcloud_min_x_collection": 4.2e-15,
    "dcloud_min_x_collision": 4.2e-15,
    "dcloud_min_x_depo": 4.2e-15,
    "dcloud_min_x_riming": 4.2e-15,
    "dcloud_min_x_freezing": 4.2e-15,
    "dcloud_cap": 2,
    "dcloud_cmu4": 0,
    "dcloud_min_x_evap": 4.2e-15,
    "dcloud_min_x_melt": 4.2e-15,
    "dcloud_min_x_nuc_hetero": 4.2e-15,
    "dcloud_min_x": 4.2e-15,
    "dcloud_min_x_nuc_homo": 4.2e-15,
    "dcloud_b_geo": 0.333333,
    "dcloud_a_vel": 375000,
    "dcloud_sc_theta_n": 0,
    "dcloud_min_x_act": 4.2e-15,
    "dcloud_a_geo": 0.124,
    "drain_vsedi_max": 20,
    "drain_vsedi_min": 0.1,
    "drain_lambda": 0.166667,
    "drain_alfa_n": 103.665,
    "drain_a_ven": 0.78,
    "drain_d_crit_c": 0,
    "drain_q_crit_c": 0,
    "drain_b_ven": 0.308,
    "drain_nm3": 7,
    "drain_nm2": 4,
    "drain_nm1": 1,
    "drain_ecoll_c": 0,
    "drain_g2": 6,
    "drain_c_s": 0.5,
    "drain_g1": 1,
    "drain_gamma": 622.2,
    "drain_beta": 9.623,
    "drain_alpha": 9.292,
    "drain_cmu5": 2,
    "drain_nu": -0.666667,
    "drain_cmu3": 0.0011,
    "drain_cmu2": 1000,
    "drain_cmu1": 30,
    "drain_b_f": 41.1318,
    "drain_c_z": 20,
    "drain_rho_v": 0,
    "drain_sc_coll_n": 0,
    "drain_b_vel": 0.23437,
    "drain_cmu0": 6,
    "drain_s_vel": 0,
    "drain_sc_delta_n": 0,
    "drain_sc_delta_q": 0,
    "drain_mu": 0.333333,
    "drain_sc_theta_q": 0,
    "drain_alfa_q": 294.546,
    "drain_max_x": 3e-06,
    "drain_min_x_sedimentation": 2.6e-10,
    "drain_a_f": 0.429251,
    "drain_min_x_conversion": 2.6e-10,
    "drain_min_x_collection": 2.6e-10,
    "drain_min_x_collision": 2.6e-10,
    "drain_min_x_depo": 2.6e-10,
    "drain_min_x_riming": 2.6e-10,
    "drain_min_x_freezing": 2.6e-10,
    "drain_cap": 2,
    "drain_cmu4": 1,
    "drain_min_x_evap": 2.6e-10,
    "drain_min_x_melt": 2.6e-10,
    "drain_min_x_nuc_hetero": 2.6e-10,
    "drain_min_x": 2.6e-10,
    "drain_min_x_nuc_homo": 2.6e-10,
    "drain_b_geo": 0.333333,
    "drain_a_vel": 114.014,
    "drain_sc_theta_n": 0,
    "drain_min_x_act": 2.6e-10,
    "drain_a_geo": 0.124,
    "dice_vsedi_max": 3,
    "dice_vsedi_min": 0,
    "dice_lambda": 0.00297619,
    "dice_alfa_n": 86.7053,
    "dice_a_ven": 0.78,
    "dice_d_crit_c": 0.00015,
    "dice_q_crit_c": 1e-05,
    "dice_b_ven": 0.308,
    "dice_nm3": 0,
    "dice_nm2": 0,
    "dice_nm1": 0,
    "dice_ecoll_c": 0.8,
    "dice_g2": 0,
    "dice_c_s": 0.5,
    "dice_g1": 0,
    "dice_gamma": 0,
    "dice_beta": 0,
    "dice_alpha": 0,
    "dice_cmu5": 0,
    "dice_nu": 1,
    "dice_cmu3": 0,
    "dice_cmu2": 0,
    "dice_cmu1": 0,
    "dice_b_f": 62.0574,
    "dice_c_z": 2.94643,
    "dice_rho_v": 0,
    "dice_sc_coll_n": 0,
    "dice_b_vel": 0.21579,
    "dice_cmu0": 0,
    "dice_s_vel": 0.05,
    "dice_sc_delta_n": 3.26042,
    "dice_sc_delta_q": 5.39197,
    "dice_mu": 0.333333,
    "dice_sc_theta_q": 0.193013,
    "dice_alfa_q": 113.436,
    "dice_max_x": 1e-05,
    "dice_min_x_sedimentation": 1e-12,
    "dice_a_f": 0.667109,
    "dice_min_x_conversion": 1e-12,
    "dice_min_x_collection": 1e-12,
    "dice_min_x_collision": 1e-12,
    "dice_min_x_depo": 1e-12,
    "dice_min_x_riming": 1e-12,
    "dice_min_x_freezing": 1e-12,
    "dice_cap": 2,
    "dice_cmu4": 0,
    "dice_min_x_evap": 1e-12,
    "dice_min_x_melt": 1e-12,
    "dice_min_x_nuc_hetero": 1e-12,
    "dice_min_x": 1e-12,
    "dice_min_x_nuc_homo": 1e-12,
    "dice_b_geo": 0.39,
    "dice_a_vel": 27.7,
    "dice_sc_theta_n": 0.124677,
    "dice_min_x_act": 1e-12,
    "dice_a_geo": 0.835,
    "dsnow_vsedi_max": 3,
    "dsnow_vsedi_min": 0.1,
    "dsnow_lambda": 0.00297619,
    "dsnow_alfa_n": 19.3054,
    "dsnow_a_ven": 0.78,
    "dsnow_d_crit_c": 0.00015,
    "dsnow_q_crit_c": 1e-05,
    "dsnow_b_ven": 0.308,
    "dsnow_nm3": 0,
    "dsnow_nm2": 0,
    "dsnow_nm1": 0,
    "dsnow_ecoll_c": 0.8,
    "dsnow_g2": 0,
    "dsnow_c_s": 0.5,
    "dsnow_g1": 0,
    "dsnow_gamma": 0,
    "dsnow_beta": 0,
    "dsnow_alpha": 0,
    "dsnow_cmu5": 0,
    "dsnow_nu": 1,
    "dsnow_cmu3": 0,
    "dsnow_cmu2": 0,
    "dsnow_cmu1": 0,
    "dsnow_b_f": 63.268,
    "dsnow_c_z": 2.94643,
    "dsnow_rho_v": 0,
    "dsnow_sc_coll_n": 0,
    "dsnow_b_vel": 0.15,
    "dsnow_cmu0": 0,
    "dsnow_s_vel": 0.25,
    "dsnow_sc_delta_n": 3.35045,
    "dsnow_sc_delta_q": 0,
    "dsnow_mu": 0.333333,
    "dsnow_sc_theta_q": 0,
    "dsnow_alfa_q": 23.3299,
    "dsnow_max_x": 2e-05,
    "dsnow_min_x_sedimentation": 1e-10,
    "dsnow_a_f": 0.663731,
    "dsnow_min_x_conversion": 1e-10,
    "dsnow_min_x_collection": 1e-10,
    "dsnow_min_x_collision": 1e-10,
    "dsnow_min_x_depo": 1e-10,
    "dsnow_min_x_riming": 1e-10,
    "dsnow_min_x_freezing": 1e-10,
    "dsnow_cap": 2,
    "dsnow_cmu4": 0,
    "dsnow_min_x_evap": 1e-10,
    "dsnow_min_x_melt": 1e-10,
    "dsnow_min_x_nuc_hetero": 1e-10,
    "dsnow_min_x": 1e-10,
    "dsnow_min_x_nuc_homo": 1e-10,
    "dsnow_b_geo": 0.455,
    "dsnow_a_vel": 8.8,
    "dsnow_sc_theta_n": 0.0562866,
    "dsnow_min_x_act": 1e-10,
    "dsnow_a_geo": 2.4,
    "dgraupel_vsedi_max": 30,
    "dgraupel_vsedi_min": 0.1,
    "dgraupel_lambda": 0.00297619,
    "dgraupel_alfa_n": 362.91,
    "dgraupel_a_ven": 0.78,
    "dgraupel_d_crit_c": 0.0001,
    "dgraupel_q_crit_c": 1e-06,
    "dgraupel_b_ven": 0.308,
    "dgraupel_nm3": 0,
    "dgraupel_nm2": 9,
    "dgraupel_nm1": 6,
    "dgraupel_ecoll_c": 1,
    "dgraupel_g2": 40320,
    "dgraupel_c_s": 0.5,
    "dgraupel_g1": 120,
    "dgraupel_gamma": 0,
    "dgraupel_beta": 0,
    "dgraupel_alpha": 0,
    "dgraupel_cmu5": 0,
    "dgraupel_nu": 1,
    "dgraupel_cmu3": 0,
    "dgraupel_cmu2": 0,
    "dgraupel_cmu1": 0,
    "dgraupel_b_f": 60.9328,
    "dgraupel_c_z": 2.94643,
    "dgraupel_rho_v": 0,
    "dgraupel_sc_coll_n": 0.563476,
    "dgraupel_b_vel": 0.268325,
    "dgraupel_cmu0": 0,
    "dgraupel_s_vel": 0,
    "dgraupel_sc_delta_n": 0,
    "dgraupel_sc_delta_q": 0,
    "dgraupel_mu": 0.333333,
    "dgraupel_sc_theta_q": 0,
    "dgraupel_alfa_q": 505.11,
    "dgraupel_max_x": 0.00053,
    "dgraupel_min_x_sedimentation": 4.19e-09,
    "dgraupel_a_f": 0.675949,
    "dgraupel_min_x_conversion": 4.19e-09,
    "dgraupel_min_x_collection": 4.19e-09,
    "dgraupel_min_x_collision": 4.19e-09,
    "dgraupel_min_x_depo": 4.19e-09,
    "dgraupel_min_x_riming": 4.19e-09,
    "dgraupel_min_x_freezing": 4.19e-09,
    "dgraupel_cap": 2,
    "dgraupel_cmu4": 0,
    "dgraupel_min_x_evap": 4.19e-09,
    "dgraupel_min_x_melt": 4.19e-09,
    "dgraupel_min_x_nuc_hetero": 4.19e-09,
    "dgraupel_min_x": 4.19e-09,
    "dgraupel_min_x_nuc_homo": 4.19e-09,
    "dgraupel_b_geo": 0.314,
    "dgraupel_a_vel": 86.8937,
    "dgraupel_sc_theta_n": 0,
    "dgraupel_min_x_act": 4.19e-09,
    "dgraupel_a_geo": 0.142,
    "dhail_vsedi_max": 30,
    "dhail_vsedi_min": 0.1,
    "dhail_lambda": 0.00297619,
    "dhail_alfa_n": 94.2824,
    "dhail_a_ven": 0.78,
    "dhail_d_crit_c": 0.0001,
    "dhail_q_crit_c": 1e-06,
    "dhail_b_ven": 0.308,
    "dhail_nm3": 0,
    "dhail_nm2": 0,
    "dhail_nm1": 0,
    "dhail_ecoll_c": 1,
    "dhail_g2": 0,
    "dhail_c_s": 0.5,
    "dhail_g1": 0,
    "dhail_gamma": 0,
    "dhail_beta": 0,
    "dhail_alpha": 0,
    "dhail_cmu5": 0,
    "dhail_nu": 1,
    "dhail_cmu3": 0,
    "dhail_cmu2": 0,
    "dhail_cmu1": 0,
    "dhail_b_f": 60.7445,
    "dhail_c_z": 2.94643,
    "dhail_rho_v": 0,
    "dhail_sc_coll_n": 0,
    "dhail_b_vel": 0.166667,
    "dhail_cmu0": 0,
    "dhail_s_vel": 0,
    "dhail_sc_delta_n": 0,
    "dhail_sc_delta_q": 0,
    "dhail_mu": 0.333333,
    "dhail_sc_theta_q": 0,
    "dhail_alfa_q": 116.275,
    "dhail_max_x": 0.00054,
    "dhail_min_x_sedimentation": 2.6e-09,
    "dhail_a_f": 0.673182,
    "dhail_min_x_conversion": 2.6e-09,
    "dhail_min_x_collection": 2.6e-09,
    "dhail_min_x_collision": 2.6e-09,
    "dhail_min_x_depo": 2.6e-09,
    "dhail_min_x_riming": 2.6e-09,
    "dhail_min_x_freezing": 2.6e-09,
    "dhail_cap": 2,
    "dhail_cmu4": 0,
    "dhail_min_x_evap": 2.6e-09,
    "dhail_min_x_melt": 2.6e-09,
    "dhail_min_x_nuc_hetero": 2.6e-09,
    "dhail_min_x": 2.6e-09,
    "dhail_min_x_nuc_homo": 2.6e-09,
    "dhail_b_geo": 0.333333,
    "dhail_a_vel": 39.3,
    "dhail_sc_theta_n": 0,
    "dhail_min_x_act": 2.6e-09,
    "dhail_a_geo": 0.1366,
}


# r"\citeA{seifert_two-moment_2006}"
# Mapping from input parameter to description, notation in different paper, different paper
# Entry 4 shall be dependent or independent parameter
# Entry 5 is a list of the involved processes
# What about spherical parameters?
in_params_notation_mapping = {
    "da_1": [
        "Dimensional coefficient used in one-moment warm physics for qc and qr calculation",
        "",
        "",
        "independent",
        [],
    ],
    "da_2": [
        "Dimensional coefficient used in one-moment warm physics for qc and qr calculation",
        "",
        "",
        "independent",
        [],
    ],
    "de_1": [
        "Dimensional coefficients used in one-moment warm physics for temperature calculation",
        "",
        "",
        "independent",
        [],
    ],
    "de_2": [
        "Dimensional coefficients used in one-moment warm physics for temperature calculation",
        "",
        "",
        "independent",
        [],
    ],
    "dd": [
        "Dimensional coefficient used in one-moment warm physics qr calculation for sedimentation",
        "",
        "",
        "independent",
        ["ccn_act_hande", "snow_melting", "graupel_melting", "hail_melting"],
    ],
    "dN_c": [
        "Number concentration of cloud droplets needed for one-moment warm physics",
        "",
        "",
        "independent",
        [],
    ],
    "dgamma": [
        "Exponent used in one-moment warm physics for qc and qr calculation",
        "",
        "",
        "independent",
        [],
    ],
    "dbeta_c": [
        "Exponent used in one-moment warm physics for qc and qr calculation",
        "",
        "",
        "independent",
        [],
    ],
    "dbeta_r": [
        "Exponent used in one-moment warm physics for qc and qr calculation",
        "",
        "",
        "independent",
        [],
    ],
    "ddelta1": [
        "Exponent used in one-moment warm physics for qv, qr, saturation and temperature calculation",
        "",
        "",
        "independent",
        [],
    ],
    "ddelta2": [
        "Exponent used in one-moment warm physics for qv, qr, saturation and temperature calculation",
        "",
        "",
        "independent",
        [],
    ],
    "dzeta": [
        "Exponents used in one-moment warm physics for qr calculation",
        "",
        "",
        "independent",
        [],
    ],
    "drain_gfak": [
        "Coefficient for gamma evaluation in rain evaporation. Not tracked.",
        "",
        "",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "dcloud_k_au": [
        "Coefficient for autoconversion of cloud to rain. Not tracked.",
        r"$ \frac{k_{cc}}{20x^\ast} \frac{(\nu_c+2)(\nu_c+4)}{(\nu_c+1)^2} $",
        r"\citeA{seifert_two-moment_2006}",
        "dependent",
        ["auto_conversion_sb"],
    ],
    "dcloud_k_sc": [
        "Coefficient for autoconversion of cloud to rain",
        r"$ k_{cc} \frac{\nu+2}{\nu+1} $",
        r"\citeA{seifert_two-moment_2006}",
        "dependent",
        ["auto_conversion_sb"],
    ],
    "dkc_autocon": [
        "Kernel for autoconversion",
        "$ k_{cc} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [],
    ],
    "dinv_z": [
        "Inverse of air parcel size (height) used in explicit sedimentation",
        "-",
        "Necessary for our box model simulation",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dw": [
        "Change in buoancy. Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dq_crit_i": [
        "Threshold (mass density) for ice selfcollection",
        r"Similar to $ \overline{D}_{i, 0} $ but for mass density",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["ice_self_collection", "auto_conversion_kb"],
    ],
    "dD_crit_i": [
        "Threshold (diameter) for ice selfcollection",
        r"$ \overline{D}_{i, 0} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["ice_self_collection"],
    ],
    "dD_conv_i": [
        "Minimum mean size of ice particles to determine the number of final ice particles in ice selfcollection",
        r"Similar to $ \overline{x}_{i, \mathrm{min}} $ but only for ice self collection",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["ice_self_collection", "ice_riming"],
    ],
    "dq_crit_r": [
        "Threshold (mass density) for ice rain riming and snow rain riming",
        r"Similar to $ \overline{D}_{c, 0} $ but for rain mass density",
        r"\citeA{seifert_two-moment_2006}, Eq. (65)",
        "independent",
        ["riming_rain_core"],
    ],
    "dD_crit_r": [
        "Threshold (diameter) for ice rain riming and snow rain riming",
        r"Similar to $ \overline{D}_{c, 0} $ but for rain",
        r"\citeA{seifert_two-moment_2006}, Eq. (65)",
        "independent",
        ["riming_rain_core"],
    ],
    "dq_crit_fr": [
        "Threshold (mass density) for instantaneous rain freezing to ice, graupel and hail.",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        ["rain_freeze"],
    ],
    "dD_coll_c": [
        "Upper bound for diameter in collision efficiency. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dq_crit": [
        "Default threshold (ratio mass). Not tracked",
        "",
        "",
        "independent",
        [
            "ice_self_collection",
            "snow_self_collection",
            "auto_conversion_kb",
            "auto_conversion_sb",
            "rain_evaporation_sb",
            "sedimentation_explicit",
            "particle_particle_collection",
            "hail_collision",
            "riming_rain_core",
            "particle_rain_riming",
            "rain_freeze",
        ],
    ],
    "dD_conv_sg": [
        "Threshold (diameter) for conversion snow to graupel",
        r"$ \overline{D}_s $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["snow_riming"],
    ],
    "dD_conv_ig": [
        "Threshold (diameter) for conversion ice to graupel",
        r"$ \overline{D}_i $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["ice_riming"],
    ],
    "dx_conv": [
        "Minimum mean mass for conversion of ice or snow to graupel",
        "-",
        "Introduced to avoid the conversion of too many particles due to too small mean diameter",
        "independent",
        ["ice_riming", "snow_riming"],
    ],
    "dparcel_height": [
        "Height of the trajectory package. Not tracked",
        "-",
        "-",
        "independent",
        [],
    ],
    "dalpha_spacefilling": [
        "Coefficient used in riming processes for enhanced melting or conversion of snow or ice to graupel. "
        "Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dT_nuc": [
        "Upper temperature threshold for heterogeneous ice nucleation",
        r"Hidden in $\xi$",
        r"\citeA{phillips_empirical_2008}, page 2763",
        "independent",
        ["ice_activation_phillips"],
    ],
    "dT_freeze": [
        "Threshold for freezing of water. Not tracked.",
        "$ T_3 $",
        r"\cite{seifert_two-moment_2006}",
        "independent",
        [
            "saturation_adjust",
            "ice_nuc_hom",
            "ice_activation_phillips",
            "cloud_freeze_hom",
            "snow_self_collection",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_particle_collection",
            "ice_riming",
            "snow_riming",
            "particle_cloud_riming",
            "particle_rain_riming",
            "rain_freeze",
            "ice_melting",
        ],
    ],
    "dT_f": [
        "Lower temperature threshold for (instantaneous) raindrop freezing. Not tracked.",
        "-",
        "-",
        "independent",
        [
            "saturation_adjust",
            "ice_nuc_hom",
            "ice_activation_phillips",
            "cloud_freeze_hom",
            "snow_self_collection",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_particle_collection",
            "ice_riming",
            "snow_riming",
            "particle_cloud_riming",
            "particle_rain_riming",
            "rain_freeze",
            "ice_melting",
        ],
    ],
    "dD_eq": [
        "Equilibrium diameter for Seifert & Beheng (2008), ie Eq. (20). Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "drho_w": [
        r"Density of liquid water in $ \mathrm{kg}/\mathrm{m}^3 $. Not tracked.",
        r"$ \rho_w $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["cloud_freeze_hom", "rain_evaporation_sb", "rain_freeze"],
    ],
    "drho_0": [
        "Air density at surface conditions. Not tracked.",
        r"$ \rho_0 $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [],
    ],
    "drho_vel": [
        "Exponent for density correction",
        r"$ \gamma $",
        r"\citeA{seifert_two-moment_2006}, Eq. (33)",
        "independent",
        [],
    ],
    "drho_vel_c": [
        "Exponent for density correction of cloud droplets",
        r"$ \gamma $",
        r"\citeA{seifert_two-moment_2006}, Eq. (33)",
        "independent",
        [],
    ],
    "drho_ice": [
        r"Density of ice in $ \mathrm{kg}/\mathrm{m}^3 $. Not tracked.",
        r"$ \rho_\epsilon $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["ice_nuc_hom", "ice_riming", "snow_riming"],
    ],
    "dM_w": [
        "Molar mass of water used for homogeneous nuclation of ice. Not tracked.",
        "$ M_w $",
        r"\citeA{karcher_physically_2006}",
        "independent",
        [
            "ice_nuc_hom",
            "ice_activation_phillips",
            "cloud_freeze_hom",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "evaporation",
            "vapor_dep_relaxation",
            "ice_riming",
            "snow_riming",
            "particle_cloud_riming",
            "particle_rain_riming",
            "rain_freeze",
            "ice_melting",
        ],
    ],
    "dM_a": [
        "Molar mass of dry air. Not tracked.",
        "$ M $",
        r"\citeA{karcher_physically_2006}",
        "independent",
        ["ice_nuc_hom"],
    ],
    "dR_universal": [
        "Universal gas constant. Not tracked",
        "$ R $",
        r"\citeA{karcher_physically_2006}",
        "independent",
        [],
    ],
    "dEpsilon": [
        "Quotient of the gas constants for dry air and for water vapor used in saturation adjustment "
        "and for updating temperature and pressure. Not tracked.",
        r"$ \epsilon $",
        r"\citeA{baumgartner_algorithmic_2019}",
        "independent",
        ["saturation_adjust", "saturation_adjust_legacy", "vapor_dep_relaxation"],
    ],
    "dgravity_acc": [
        "Gravitational acceleration used in homogeneous ice nucleation and for updating temperature and pressure. "
        "Not tracked.",
        "$ g $",
        r"\citeA{karcher_physically_2006}",
        "independent",
        ["ice_nuc_hom"],
    ],
    "dR_a": [
        "Gas constant for dry air used in homogeneous ice nucleation and for updating temperature and pressure. "
        "Not tracked.",
        "$ R_a $",
        r"\citeA{baumgartner_algorithmic_2019}",
        "independent",
        ["saturation_adjust", "ice_nuc_hom"],
    ],
    "dR_v": [
        "Gas constant for water vapor used in evaporation and melting processes. Not tracked.",
        "$ R_v $",
        r"\citeA{baumgartner_algorithmic_2019}",
        "independent",
        [
            "saturation_adjust",
            "ice_nuc_hom",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "rain_evaporation_sb",
            "evaporation",
            "vapor_dep_relaxation",
        ],
    ],
    "da_v": [
        "Constant used in rain evaporation to calculate the ventilation factor",
        "$ a_v $",
        r"\citeA{seifert_parameterization_2008}, Eq. (8)",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "db_v": [
        "Coefficient used in rain evaporation to calculate the ventilation factor",
        "$ b_v $",
        r"\citeA{seifert_parameterization_2008}, Eq. (8)",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "da_prime": [
        "Constant used to calculate the terminal fall velocity of raindrops during rain evaporation "
        "more accurately than with a power law",
        "$ a $",
        r"\citeA{seifert_parameterization_2008}, Eq. (A4)",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "db_prime": [
        "Coefficient used to calculate the terminal fall velocity of raindrops during rain evaporation "
        "more accurately than with a power law",
        "$ b $",
        r"\citeA{seifert_parameterization_2008}, Eq. (A4)",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "dc_prime": [
        "Exponent used to calculate the terminal fall velocity of raindrops during rain evaporation "
        "more accurately than with a power law",
        "$ c $",
        r"\citeA{seifert_parameterization_2008}, Eq. (A4)",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "dK_T": [
        "Heat conductivity of air used for melting and evaporation. Not tracked.",
        "$ K_T $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "rain_evaporation_sb",
            "evaporation",
            "vapor_dep_relaxation",
        ],
    ],
    "dL_wd": [
        "Latent heat of evaporation of water used in saturation ajustment, melting and evaporation. Not tracked.",
        "$ L_{lv} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [
            "saturation_adjust",
            "saturation_adjust_legacy",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "rain_evaporation_sb",
            "evaporation",
        ],
    ],
    "dL_ed": [
        " Heat of sublimation ice to vapor used in sublimation and homogeneous ice nucleation. Not tracked.",
        "$ L_{iv} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["ice_nuc_hom", "vapor_dep_relaxation"],
    ],
    "dD_v": [
        r"Diffusivity of water vapor in air at $0^\circ\mathrm{C} $ used in evaporation of ice particles. Not tracked",
        "$ D_v $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["evaporation"],
    ],
    "decoll_min": [
        "Min. efficiency for riming of cloud with graupel, ice or snow",
        r"similar to $ \overline{E}_{g, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "decoll_gg": [
        "Collision efficiency for graupel selfcollection",
        r"$ \overline{E}_{\mathrm{coll}, gg} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (64)",
        "independent",
        ["particle_particle_collection"],
    ],
    "decoll_gg_wet": [
        "Collision efficiency for wet graupel",
        r"similar to $ \overline{E}_{\mathrm{coll}, gg} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (64)",
        "independent",
        ["particle_particle_collection"],
    ],
    "dkin_visc_air": [
        "Kinematic viscosity of dry air used to calculate the terminal fall velocity "
        "of raindrops during rain evaporation more accurately than with a power law",
        r"$ \nu_{\mathrm{air}} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "dC_mult": [
        "Coefficient for splintering during Hallet-Mossop ice multiplication.",
        r"$ F_\mathrm{splint} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["ice_riming", "snow_riming", "particle_cloud_riming", "particle_rain_riming"],
    ],
    "dT_mult_min": [
        "Coefficient used in Hallet-Mossop ice multiplication",
        r"$ T_{\mathrm{splint}, \mathrm{min}} $",
        r"Based on \citeA{hallett_production_1974}; following notation from "
        r"\citeA{seifert_parametrisierung_2002} Eq. (4.127)",
        "independent",
        ["ice_riming", "snow_riming", "particle_cloud_riming", "particle_rain_riming"],
    ],
    "dT_mult_max": [
        "Coefficient used in Hallet-Mossop ice multiplication",
        r"$ T_{\mathrm{splint}, \mathrm{max}} $",
        r"Based on \citeA{hallett_production_1974}; following notation from "
        r"\citeA{seifert_parametrisierung_2002} Eq. (4.127)",
        "independent",
        ["ice_riming", "snow_riming", "particle_cloud_riming", "particle_rain_riming"],
    ],
    "dT_mult_opt": [
        "Coefficient used in Hallet-Mossop ice multiplication. Not tracked",
        r"$ T_{\mathrm{splint}, \mathrm{opt}} $",
        r"Based on \citeA{hallett_production_1974}; following notation from "
        r"\citeA{seifert_parametrisierung_2002} Eq. (4.127)",
        "independent",
        [],
    ],
    "dconst0": [
        "Coefficient used in riming processes",
        r"$ (\overline{D}_{c,1} - \overline{D}_{c, 0})^{-1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (65)",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "dconst3": [
        "Coefficient used in riming processes for breaking up particles",
        r"$( T_{\mathrm{splint}, \mathrm{opt}} - T_{\mathrm{splint}, \mathrm{min}})^{-1}$",
        r"Based on \citeA{hallett_production_1974}; following notation from "
        r"\citeA{seifert_parametrisierung_2002} Eq. (4.127)",
        "independent",
        ["ice_riming", "snow_riming", "particle_cloud_riming", "particle_rain_riming"],
    ],
    "dconst4": [
        "Coefficient used in riming processes for breaking up particles",
        r"$( T_{\mathrm{splint}, \mathrm{opt}} - T_{\mathrm{splint}, \mathrm{max}})^{-1}$",
        r"Based on \citeA{hallett_production_1974}; following notation from "
        r"\citeA{seifert_parametrisierung_2002} Eq. (4.127)",
        "independent",
        ["ice_riming", "snow_riming", "particle_cloud_riming", "particle_rain_riming"],
    ],
    "dconst5": [
        "Coefficient used for conversion of snow or ice to graupel during riming",
        r"$ \alpha_\circ \frac{\rho_q}{\rho_\epsilon} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["ice_riming", "snow_riming"],
    ],
    "dD_rainfrz_gh": [
        " Size threshold for partitioning of freezing rain in the hail scheme",
        r"Similar to $r_\ast$",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["rain_freeze"],
    ],
    "dD_rainfrz_ig": [
        "Size threshold for partitioning of freezing rain in the hail scheme",
        r"Similar to $r_\ast$",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["rain_freeze"],
    ],
    "ddv0": [
        r"Diffusivity of water vapor in air at $0^\circ\mathrm{C} $",
        "$ D_v $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["snow_melting", "graupel_melting", "hail_melting"],
    ],
    "dp_sat_melt": [
        r"Saturation pressure at $ \mathrm{T}=\mathrm{T}_\mathrm{freeze} $",
        "$ P_{lv}(T_3) $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["snow_melting", "graupel_melting", "hail_melting"],
    ],
    "dcp": [
        r"Specific heat capacity of air at constant pressure in $ \mathrm{J}/\mathrm{K}/\mathrm{kg} $. Not tracked.",
        "$ c_p $",
        r"\citeA{karcher_physically_2006}",
        "independent",
        ["saturation_adjust", "saturation_adjust_legacy", "ice_nuc_hom"],
    ],
    "dk_b": [
        r"Boltzmann constant to calculate the water vapor number density at ice saturation and "
        r"the thermal speed of water molecules. See Equations~\ref{eq:n_sat}~and~\ref{eq:v_th}.",
        "In $b_1$ and $b_2$",
        r"\citeA{karcher_physically_2006}, page 2",
        "independent",
        ["ice_nuc_hom", "rain_self_collection_sb"],
    ],
    "da_HET": [
        "Exponent for heterogeneous rain freezing with data of Barklie and Gokhale",
        r"$ B_{\mathrm{HET}} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["rain_freeze"],
    ],
    "db_HET": [
        "Coefficient for heterogeneous rain freezing with data of Barklie and Gokhale",
        r"$ A_{\mathrm{HET}} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["rain_freeze"],
    ],
    "dN_sc": [
        "Schmidt number used in rain evaporation",
        r"$ N_{\mathrm{sc}} $",
        r"\citeA{seifert_parameterization_2008}, Eqs. (8), (A3), (A6), and (A7)",
        "independent",
        [],
    ],
    "dn_f": [
        " Exponent of b_f in the vent-coeff. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dN_avo": [
        r"Avogadro number used in homogeneous ice nucleation to calculate the mass of a water molecule "
        r"$M_a$ for the thermal speed of water molecules. See Equation~\ref{eq:v_th}.",
        "In $b_1$ and $b_2$",
        r"\citeA{karcher_physically_2006}, page 2",
        "independent",
        ["ice_nuc_hom"],
    ],
    "dna_dust": [
        "Initial number density of dust",
        r"Similar to $ n_{\mathrm{IN}, \mathrm{DM}} $",
        r"\citeA{phillips_empirical_2008}",
        "independent",
        ["ice_activation_phillips"],
    ],
    "dna_soot": [
        "Initial number density of soot",
        r"$ n_{\mathrm{IN}, \mathrm{BC}} $",
        r"\citeA{phillips_empirical_2008}",
        "independent",
        ["ice_activation_phillips"],
    ],
    "dna_orga": [
        "Initial number density of organics",
        r"$ n_{\mathrm{IN}, \mathrm{O}} $",
        r"\citeA{phillips_empirical_2008}",
        "independent",
        ["ice_activation_phillips"],
    ],
    "dni_het_max": [
        "Maximum particle number for heterogeneous ice nucleation per liter",
        r"$ n_{\mathrm{IN}, X} $",
        r"\citeA{phillips_empirical_2008}, Eq. (13)",
        "independent",
        ["ice_activation_phillips"],
    ],
    "dni_hom_max": [
        "Maximum particle number for homogeneous ice nucleation per liter",
        "$ n $",
        r"\citeA{karcher_physically_2006}, Eq. (10)",
        "independent",
        ["ice_nuc_hom"],
    ],
    "da_dep": [
        "Cons_idxeters for deposition formula (2) of Hande et al. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "db_dep": [
        "Cons_idxeters for deposition formula (2) of Hande et al. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dc_dep": [
        "Parameters for deposition formula (2) of Hande et al. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dd_dep": [
        "Cons_idxeters for deposition formula (2) of Hande et al. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dnim_imm": [
        "Parameter for Hande et al. nucleation. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dnin_dep": [
        "Parameter for Hande et al. nucleation. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dalf_imm": [
        "Parameter for Hande et al. nucleation. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dbet_dep": [
        "Parameter for Hande et al. nucleation. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dbet_imm": [
        "Parameter for Hande et al. nucleation. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dhande_ccn_fac": [
        "Parameter to scale the effect of CCN activation.",
        r"$f_{\mathrm{fac}}",
        "Annette's code",
        "independent",
        ["ccn_act_hande_akm"],
    ],
    "di_ccn_1": [
        "Parameter to calculate the CCN concentration.",
        "$i_1",
        "Annette's code",
        "independent",
        ["ccn_act_hande_akm"],
    ],
    "di_ccn_2": [
        "Parameter to calculate the lower limit of CCN concentration during CCN activation.",
        "$i_2",
        "Annette's code",
        "independent",
        ["ccn_act_hande_akm"],
    ],
    "dh_ccn_1": [
        "Parameter to calculate the CCN concentration.",
        "$h_1",
        "Annette's code",
        "independent",
        ["ccn_act_hande_akm"],
    ],
    "dh_ccn_2": [
        "Parameter to calculate the CCN concentration during CCN activation.",
        "$h_2",
        "Annette's code",
        "independent",
        ["ccn_act_hande_akm"],
    ],
    "dg_ccn_1": [
        "Parameter for calculating CCN concentration during CCN activation using a formulation by Miltenberger et al.",
        "$g_1",
        "Annette's code",
        "independent",
        ["ccn_act_hande_akm"],
    ],
    "dg_ccn_2": [
        "Parameter to calculate the CCN concentration during CCN activation.",
        "$g_2",
        "Annette's code",
        "independent",
        ["ccn_act_hande_akm"],
    ],
    "dg_ccn_3": [
        "Parameter to calculate the CCN concentration during CCN activation.",
        "$g_3",
        "Annette's code",
        "independent",
        ["ccn_act_hande_akm"],
    ],
    "da_ccn_1": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ a_1 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "da_ccn_2": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ a_2 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "da_ccn_3": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ a_3 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "da_ccn_4": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ a_4 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "db_ccn_1": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ b_1 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande_akm", "ccn_act_hande"],
    ],
    "db_ccn_2": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ b_2 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande_akm", "ccn_act_hande"],
    ],
    "db_ccn_3": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ b_3 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande_akm", "ccn_act_hande"],
    ],
    "db_ccn_4": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ b_4 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "dc_ccn_1": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ c_1 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande_akm", "ccn_act_hande"],
    ],
    "dc_ccn_2": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ c_2 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande_akm", "ccn_act_hande"],
    ],
    "dc_ccn_3": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ c_3 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande_akm", "ccn_act_hande"],
    ],
    "dc_ccn_4": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ c_4 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "dd_ccn_1": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ d_1 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "dd_ccn_2": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ d_2 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "dd_ccn_3": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ d_3 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "dd_ccn_4": [
        "Parameter for calculating CCN concentration during CCN activation",
        "$ d_4 $",
        r"\citeA{hande_parameterizing_2016}",
        "independent",
        ["ccn_act_hande"],
    ],
    "dr_const": [
        "Parameter for saturation adjustment. Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dr1_const": [
        "Parameter for saturation adjustment. Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcv": [
        r"Specific heat capacity of water vapor at constant pressure in $ \mathrm{J}/\mathrm{K}/\mathrm{kg} $. "
        r"Not used",
        "",
        "",
        "independent",
        ["saturation_adjust"],
    ],
    "dp_sat_const_a": [
        "Parameter for saturation $ S $ adjustment. Not tracked.",
        "",
        "",
        "independent",
        [
            "saturation_adjust",
            "saturation_adjust_legacy",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "vapor_dep_relaxation",
        ],
    ],
    "dp_sat_ice_const_a": [
        "Parameter for saturation adjustment. Constant saturated ice pressure. Not tracked.",
        "",
        "-",
        "independent",
        ["vapor_dep_relaxation"],
    ],
    "dp_sat_const_b": [
        "Parameter for saturation adjustment. Constant saturated water vapor pressure. Not tracked.",
        "",
        "-",
        "independent",
        [
            "saturation_adjust",
            "saturation_adjust_legacy",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "vapor_dep_relaxation",
        ],
    ],
    "dp_sat_ice_const_b": [
        "Parameter for saturation adjustment. Constant saturated ice pressure. Not tracked.",
        "",
        "-",
        "independent",
        ["vapor_dep_relaxation"],
    ],
    "dp_sat_low_temp": [
        "Parameter for saturation adjustment. Saturated water vapor pressure at $T = 233 $K. Not tracked.",
        "",
        "-",
        "independent",
        [
            "saturation_adjust",
            "saturation_adjust_legacy",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "vapor_dep_relaxation",
        ],
    ],
    "dT_sat_low_temp": [
        "Parameter for saturation adjustment. Not tracked.",
        "",
        "",
        "independent",
        [
            "saturation_adjust",
            "saturation_adjust_legacy",
            "snow_melting",
            "graupel_melting",
            "hail_melting",
            "vapor_dep_relaxation",
        ],
    ],
    "dalpha_depo": [
        "Deposition coefficient of water molecules touching ice for homogeneous ice nucleation",
        r"$ \alpha $ (below Eq. (2))",
        r"\citeA{karcher_physically_2006}",
        "independent",
        ["ice_nuc_hom"],
    ],
    "dr_0": [
        "Aerosol particle radius prior to freezing used in homogeneous nucleation",
        "$ r_0 $",
        r"\citeA{karcher_physically_2006}",
        "independent",
        ["ice_nuc_hom"],
    ],
    "dk_1_conv": [
        "Exponent for autoconversion of qc to qr",
        "-",
        r"\citeA{seifert_two-moment_2006}, Eq. (6)",
        "independent",
        ["auto_conversion_sb"],
    ],
    "dk_2_conv": [
        "Exponent for autoconversion of qc to qr",
        "-",
        r"\citeA{seifert_two-moment_2006}, Eq. (6)",
        "independent",
        ["auto_conversion_sb"],
    ],
    "dk_1_accr": [
        "Coefficient for accretion of qc to qr",
        "-",
        r"\citeA{seifert_two-moment_2006}, Eq. (8)",
        "independent",
        ["auto_conversion_sb"],
    ],
    "dk_r": [
        "Coefficient for accretion of qc to qr",
        "$ k_{cr} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["auto_conversion_sb"],
    ],
    "drain_a_geo": [
        "Coefficient for diameter size calculation",
        "$ a $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [
            "rain_self_collection_sb",
            "rain_evaporation_sb",
            "sedimentation_explicit",
            "riming_rain_core",
            "particle_rain_riming",
            "rain_freeze",
        ],
    ],
    "drain_b_geo": [
        "Exponent for diameter size calculation",
        "$ b $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [
            "rain_self_collection_sb",
            "rain_evaporation_sb",
            "sedimentation_explicit",
            "riming_rain_core",
            "particle_rain_riming",
            "rain_freeze",
        ],
    ],
    "drain_min_x": [
        "Minimum size of the particle used after the microphysics",
        r"$ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "rain_self_collection_sb",
            "rain_evaporation_sb",
            "sedimentation_explicit",
            "riming_rain_core",
            "ice_riming",
            "snow_riming",
            "particle_rain_riming",
            "rain_freeze",
        ],
    ],
    "drain_min_x_act": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_min_x_nuc_homo": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_min_x_nuc_hetero": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_min_x_melt": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_min_x_evap": [
        "Minimum size of particle for evaporation",
        r"$ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "drain_min_x_freezing": [
        "Minimum size of particle for freezing",
        r"$ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["rain_freeze"],
    ],
    "drain_min_x_depo": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_min_x_collision": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_min_x_collection": [
        "Minimum size of particle for different collision processes",
        r"$ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["rain_self_collection_sb"],
    ],
    "drain_min_x_conversion": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_min_x_sedimentation": [
        "Minimum size of particle for sedimentation",
        r"$ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "drain_min_x_riming": [
        "Minimum size of particle for riming",
        r"$ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["riming_rain_core", "ice_riming", "snow_riming", "particle_rain_riming"],
    ],
    "drain_max_x": [
        "Maximum size of particle",
        r"$ \overline{x}_{r, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "rain_self_collection_sb",
            "rain_evaporation_sb",
            "sedimentation_explicit",
            "riming_rain_core",
            "ice_riming",
            "snow_riming",
            "particle_rain_riming",
            "rain_freeze",
        ],
    ],
    "drain_sc_theta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_sc_delta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_sc_theta_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_sc_delta_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_s_vel": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_a_vel": [
        "Coefficient for particle velocity",
        r"$ \alpha $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        ["riming_rain_core", "particle_rain_riming"],
    ],
    "drain_b_vel": [
        "Exponent for particle velocity",
        r"$ \beta $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        ["riming_rain_core", "particle_rain_riming"],
    ],
    "drain_rho_v": [
        "Coefficient used in density correction for the increased terminal fall velocity with decreasing air density",
        r"$ \Big( \frac{\rho_0}{\rho} \Big)^\gamma $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "dependent",
        [
            "rain_self_collection_sb",
            "rain_evaporation_sb",
            "sedimentation_explicit",
            "riming_rain_core",
            "particle_rain_riming",
        ],
    ],
    "drain_c_z": [
        "Coefficient for 2nd mass moment used in homogeneous freezing",
        r"Similar to $ \frac{\nu_c + 2}{\nu_c + 1} $",
        r"\citeA{seifert_two-moment_2006} Eq. (50)",
        "dependent",
        ["rain_freeze"],
    ],
    "drain_sc_coll_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_cmu0": [
        r"Coefficient for calculating the shape parameter $ \mu $ during rain evaporation",
        "-",
        r"\citeA{seifert_parameterization_2008}, Eq. (20)",
        "independent",
        ["rain_evaporation_sb", "sedimentation_explicit"],
    ],
    "drain_cmu1": [
        r"Coefficient for calculating the shape parameter $ \mu $ during rain evaporation",
        "-",
        r"\citeA{seifert_parameterization_2008}, Eq. (20)",
        "independent",
        ["rain_evaporation_sb", "sedimentation_explicit"],
    ],
    "drain_cmu2": [
        r"Coefficient for calculating the shape parameter $ \mu $ during rain evaporation",
        "$ c_2 $",
        r"\citeA{seifert_parameterization_2008}, Eq. (20)",
        "independent",
        ["rain_evaporation_sb", "sedimentation_explicit"],
    ],
    "drain_cmu3": [
        r"Constant for calculating the shape parameter $ \mu $ during rain evaporation",
        r"$ D_{\mathrm{eq}} $",
        r"\citeA{seifert_parameterization_2008}, Eq. (20)",
        "independent",
        ["rain_evaporation_sb", "sedimentation_explicit"],
    ],
    "drain_cmu4": [
        r"Constant for calculating the shape parameter $ \mu $ during rain evaporation",
        "-",
        r"\citeA{seifert_parameterization_2008}, Eq. (20)",
        "independent",
        ["rain_evaporation_sb", "sedimentation_explicit"],
    ],
    "drain_cmu5": [
        r"Exponent for calculating the shape parameter $ \mu $ during rain evaporation",
        "-",
        r"\citeA{seifert_parameterization_2008}, Eq. (20)",
        "independent",
        ["rain_evaporation_sb"],
    ],
    "drain_alpha": [
        "Constant in rain sedimentation",
        "$ a $",
        r"\citeA{seifert_parameterization_2008}, Eq. (A10)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "drain_beta": [
        "Coefficient for rain sedimentation",
        "$ b $",
        r"\citeA{seifert_parameterization_2008}, Eq. (A10)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "drain_gamma": [
        "Exponent for rain sedimentation",
        "$ c $",
        r"\citeA{seifert_parameterization_2008}, Eq. (A10)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "drain_nu": [
        r"Parameter to calculate the shape of the generalized $ \Gamma $-distribution",
        r"$ \mu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "drain_g1": [
        "Right edge of incomplete gamma function",
        r"$ \Gamma\Big( \frac{\nu+1}{\mu} \Big) $",
        r"\citeA{seifert_two-moment_2006}, Eq. (80)",
        "dependent",
        ["rain_freeze"],
    ],
    "drain_g2": [
        "Right edge of incomplete gamma function",
        r"$ \Gamma\Big( \frac{\nu+2}{\mu} \Big) $",
        r"\citeA{seifert_two-moment_2006}, Eq. (80)",
        "dependent",
        ["rain_freeze"],
    ],
    "drain_mu": [
        r"Shape parameter of the generalized $ \Gamma $-distribution",
        r"$ \nu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        ["rain_freeze"],
    ],
    "drain_nm1": [
        r"Argument of incomplete $ \Gamma $-distribution for $g_{1, \mathrm{rain}}$",
        r"$ \frac{\nu+1}{\mu} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (80)",
        "dependent",
        ["rain_freeze"],
    ],
    "drain_nm2": [
        r"Argument of incomplete $ \Gamma $-distribution for $g_{2, \mathrm{rain}}$",
        r"$ \frac{\nu+2}{\mu} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (80)",
        "dependent",
        ["rain_freeze"],
    ],
    "drain_nm3": [
        r"Argument of incomplete $ \Gamma $-distribution",
        r"$ \frac{\nu+3}{\mu} $",
        r"\citeA{seifert_two-moment_2006}, see $Z_r$ in Eq. (47)",
        "dependent",
        ["rain_freeze"],
    ],
    "drain_q_crit_c": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_d_crit_c": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_ecoll_c": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_cap": [
        "Coefficient for capacity of particle. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "drain_a_ven": [
        "Used to calculate the constant for average ventilation. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "drain_b_ven": [
        "Used to calculate the constant for average ventilation. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "drain_c_s": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "drain_a_f": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "drain_b_f": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "drain_alfa_n": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "drain_alfa_q": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "drain_lambda": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "drain_vsedi_min": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "drain_vsedi_max": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_a_geo": [
        "Coefficient for diameter size calculation",
        "$ a $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "dcloud_b_geo": [
        "Exponent for diameter size calculation",
        "$ b $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "dcloud_min_x": [
        "Minimum size of the particle used after the microphysics",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "ccn_act_hande_akm",
            "ccn_act_hande",
            "cloud_freeze_hom",
            "auto_conversion_kb",
            "auto_conversion_sb",
            "riming_cloud_core",
            "particle_cloud_riming",
        ],
    ],
    "dcloud_min_x_act": [
        "Minimum size of particle for CCN activation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["ccn_act_hande_akm", "ccn_act_hande"],
    ],
    "dcloud_min_x_nuc_homo": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_min_x_nuc_hetero": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_min_x_melt": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_min_x_evap": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_min_x_freezing": [
        "Minimum size of particle for freezing",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["cloud_freeze_hom"],
    ],
    "dcloud_min_x_depo": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_min_x_collision": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_min_x_collection": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_min_x_conversion": [
        "Minimum size of particle for conversion processes",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["auto_conversion_kb", "auto_conversion_sb"],
    ],
    "dcloud_min_x_sedimentation": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_min_x_riming": [
        "Minimum size of particle for riming",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "dcloud_max_x": [
        "Maximum size of particle",
        r"similar to $ \overline{x}_{r, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "cloud_freeze_hom",
            "auto_conversion_kb",
            "auto_conversion_sb",
            "riming_cloud_core",
            "particle_cloud_riming",
            "ice_melting",
        ],
    ],
    "dcloud_sc_theta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_sc_delta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_sc_theta_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_sc_delta_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_s_vel": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_a_vel": [
        "Coefficient for particle velocity",
        r"$ \alpha $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "dcloud_b_vel": [
        "Exponent for particle velocity",
        r"$ \beta $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "dcloud_rho_v": [
        "Coefficient used in density correction for the increased terminal fall velocity with decreasing air density",
        r"$ \Big( \frac{\rho_0}{\rho} \Big)^\gamma $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "dependent",
        ["auto_conversion_sb", "riming_cloud_core", "particle_cloud_riming"],
    ],
    "dcloud_c_z": [
        "Coefficient for 2nd mass moment used in homogeneous freezing. Not tracked.",
        r"Similar to $ \frac{\nu_c + 2}{\nu_c + 1} $",
        r"\citeA{seifert_two-moment_2006} Eq. (50)",
        "dependent",
        ["cloud_freeze_hom"],
    ],
    "dcloud_sc_coll_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_cmu0": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_cmu1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_cmu2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_cmu3": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_cmu4": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_cmu5": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_alpha": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_beta": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_gamma": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_nu": [
        r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \mu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        [],
    ],
    "dcloud_g1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_g2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_mu": [
        r"Shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \nu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        [],
    ],
    "dcloud_nm1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_nm2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_nm3": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_q_crit_c": [
        "Threshold (mass density) in riming",
        r"Similar to $ \overline{D}_{c, 0} $ but for mass density",
        r"\citeA{seifert_two-moment_2006}, Eq. (65)",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "dcloud_d_crit_c": [
        "Threshold (diameter) in riming",
        r"$ \overline{D}_{c, 0} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (65)",
        "independent",
        ["riming_cloud_core", "particle_cloud_riming"],
    ],
    "dcloud_ecoll_c": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_cap": [
        "Coefficient for capacity of particle. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_a_ven": [
        "Used to calculate the constant for average ventilation. Not tracked",
        "$ a_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "independent",
        [],
    ],
    "dcloud_b_ven": [
        "Used to calculate the constant for average ventilation. Not tracked",
        "$ b_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "independent",
        [],
    ],
    "dcloud_c_s": [
        "Inverse of capacity. Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dcloud_a_f": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dcloud_b_f": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dcloud_alfa_n": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dcloud_alfa_q": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dcloud_lambda": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dcloud_vsedi_min": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dcloud_vsedi_max": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_a_geo": [
        "Coefficient for diameter size calculation",
        "$ a $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [
            "graupel_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "particle_particle_collection",
            "graupel_hail_conv",
            "particle_cloud_riming",
            "particle_rain_riming",
        ],
    ],
    "dgraupel_b_geo": [
        "Exponent for diameter size calculation",
        "$ b $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [
            "graupel_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "particle_particle_collection",
            "graupel_hail_conv",
            "particle_cloud_riming",
            "particle_rain_riming",
        ],
    ],
    "dgraupel_min_x": [
        "Minimum size of the particle used after the microphysics",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "graupel_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "particle_particle_collection",
            "graupel_hail_conv",
            "particle_cloud_riming",
            "particle_rain_riming",
        ],
    ],
    "dgraupel_min_x_act": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_min_x_nuc_homo": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_min_x_nuc_hetero": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_min_x_melt": [
        "Minimum size of particle for melting",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["graupel_melting"],
    ],
    "dgraupel_min_x_evap": [
        "Minimum size of particle for evaporation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["evaporation"],
    ],
    "dgraupel_min_x_freezing": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_min_x_depo": [
        "Minimum size of particle for vapor deposition",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["vapor_dep_relaxation"],
    ],
    "dgraupel_min_x_collision": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_min_x_collection": [
        "Minimum size of particle for different collision processes",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["particle_collection", "particle_particle_collection"],
    ],
    "dgraupel_min_x_conversion": [
        "Minimum size of particle for conversion processes",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["graupel_hail_conv"],
    ],
    "dgraupel_min_x_sedimentation": [
        "Minimum size of particle for sedimentation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_min_x_riming": [
        "Minimum size of particle for riming",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["particle_cloud_riming", "particle_rain_riming"],
    ],
    "dgraupel_max_x": [
        "Maximum size of particle",
        r"similar to $ \overline{x}_{r, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "graupel_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "particle_particle_collection",
            "graupel_hail_conv",
            "particle_cloud_riming",
            "particle_rain_riming",
        ],
    ],
    "dgraupel_sc_theta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_sc_delta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_sc_theta_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_sc_delta_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_s_vel": [
        "Variance for the assumed Gaussian velocity distributions used in collection and riming processes. Not used.",
        r"similar to $ \sigma_i $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["particle_collection"],
    ],
    "dgraupel_a_vel": [
        "Coefficient for particle velocity",
        r"$ \alpha $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        [
            "graupel_melting",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "particle_particle_collection",
            "particle_cloud_riming",
            "particle_rain_riming",
        ],
    ],
    "dgraupel_b_vel": [
        "Exponent for particle velocity",
        r"$ \beta $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        [
            "graupel_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "particle_particle_collection",
            "particle_cloud_riming",
            "particle_rain_riming",
        ],
    ],
    "dgraupel_rho_v": [
        "Coefficient used in density correction for the increased terminal fall velocity with decreasing air density",
        r"$ \Big( \frac{\rho_0}{\rho} \Big)^\gamma $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "dependent",
        [
            "graupel_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "particle_particle_collection",
            "particle_cloud_riming",
            "particle_rain_riming",
        ],
    ],
    "dgraupel_c_z": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dgraupel_sc_coll_n": [
        "Coefficient in graupel self collection and cloud riming",
        r"combination of $ \vartheta $ and $ \delta $",
        r"\citeA{seifert_two-moment_2006}, Eq. (62)",
        "dependent",
        ["particle_particle_collection"],
    ],
    "dgraupel_cmu0": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_cmu1": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_cmu2": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_cmu3": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_cmu4": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_cmu5": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_alpha": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_beta": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_gamma": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_nu": [
        r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \mu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_g1": [
        "Right edge of incomplete gamma function",
        r"$ \Gamma\Big( \frac{\nu+1}{\mu} \Big) $",
        r"\citeA{seifert_two-moment_2006}, Eq. (80)",
        "dependent",
        ["graupel_hail_conv"],
    ],
    "dgraupel_g2": [
        "Right edge of incomplete gamma function",
        r"$ \Gamma\Big( \frac{\nu+2}{\mu} \Big) $",
        r"\citeA{seifert_two-moment_2006}, Eq. (80)",
        "dependent",
        ["graupel_hail_conv"],
    ],
    "dgraupel_mu": [
        r"Shape parameter of the generalized $ \Gamma $-distribution",
        r"$ \nu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        ["graupel_hail_conv"],
    ],
    "dgraupel_nm1": [
        r"Argument of incomplete $ \Gamma $-distribution for $g_{1, \mathrm{rain}}$",
        r"$ \frac{\nu+1}{\mu} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (80)",
        "dependent",
        ["graupel_hail_conv"],
    ],
    "dgraupel_nm2": [
        r"Argument of incomplete $ \Gamma $-distribution for $g_{2, \mathrm{rain}}$",
        r"$ \frac{\nu+2}{\mu} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (80)",
        "dependent",
        ["graupel_hail_conv"],
    ],
    "dgraupel_nm3": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dgraupel_q_crit_c": [
        "Threshold (mass density) in riming",
        r"Similar to $ \overline{D}_{g, 0} $ but for mass density",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["graupel_hail_conv", "particle_cloud_riming"],
    ],
    "dgraupel_d_crit_c": [
        "Threshold (diameter) in riming",
        r"$ \overline{D}_{g, 0} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["particle_cloud_riming"],
    ],
    "dgraupel_ecoll_c": [
        "Riming coefficient. Maximum collision efficiency with cloud droplets",
        r"$ \overline{E}_{g, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["particle_cloud_riming"],
    ],
    "dgraupel_cap": [
        "Coefficient for capacity of particle. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dgraupel_a_ven": [
        "Used to calculate the constant for average ventilation",
        "$ a_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "independent",
        ["vapor_dep_relaxation"],
    ],
    "dgraupel_b_ven": [
        "Used to calculate the constant for average ventilation. Not tracked",
        "$ b_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "independent",
        [],
    ],
    "dgraupel_c_s": [
        "Coefficient in evaporation and vapor deposition for capacity",
        "$ c_g^{-1} $",
        r"\citeA{seifert_two-moment_2006}",
        "dependent",
        ["evaporation", "vapor_dep_relaxation"],
    ],
    "dgraupel_a_f": [
        "Constant for average ventilation. Used in melting and ice-vapor processes",
        r"$ \overline{a}_{\mathrm{vent}, 1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "dependent",
        ["graupel_melting", "evaporation", "vapor_dep_relaxation"],
    ],
    "dgraupel_b_f": [
        "Coefficient for average ventilation",
        r"$ \overline{b}_{\mathrm{vent}, 1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "dependent",
        ["graupel_melting", "evaporation", "vapor_dep_relaxation"],
    ],
    "dgraupel_alfa_n": [
        "Sedimentation velocity coefficient",
        r"Similar to $ \alpha_g \Gamma( \frac{k+\nu_g+\beta_g+1}{\mu_g}) "
        r"\Big( \Gamma( \frac{\nu_g+2}{\mu_g}) \Big)^{-1} $ but for particle number",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_alfa_q": [
        "Sedimentation velocity coefficient",
        r"$ \alpha_g \Gamma( \frac{k+\nu_g+\beta_g+1}{\mu_g}) \Big( \Gamma( \frac{\nu_g+2}{\mu_g}) \Big)^{-1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_lambda": [
        "Sedimentation velocity coefficient",
        r"$ \Gamma \frac{\nu_g+1}{\mu_g} \Big( \frac{\nu_g+2}{\mu_g} \Big)^{-1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_vsedi_min": [
        "Minimum sedimentation velocity parameter",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dgraupel_vsedi_max": [
        "Maximum sedimentation velocity parameter",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dhail_a_geo": [
        "Coefficient for diameter size calculation",
        "$ a $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [],
    ],
    "dhail_b_geo": [
        "Exponent for diameter size calculation",
        "$ b $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [],
    ],
    "dhail_min_x": [
        "Minimum size of the particle used after the microphysics",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [],
    ],
    "dhail_min_x_act": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_min_x_nuc_homo": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_min_x_nuc_hetero": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_min_x_melt": [
        "Minimum size of particle for melting",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [],
    ],
    "dhail_min_x_evap": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_min_x_freezing": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_min_x_depo": [
        "Minimum size of particle for vapor deposition",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [],
    ],
    "dhail_min_x_collision": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_min_x_collection": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_min_x_conversion": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_min_x_sedimentation": [
        "Minimum size of particle for sedimentation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [],
    ],
    "dhail_min_x_riming": [
        "Minimum size of particle for riming",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [],
    ],
    "dhail_max_x": [
        "Maximum size of particle",
        r"similar to $ \overline{x}_{r, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [],
    ],
    "dhail_sc_theta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_sc_delta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_sc_theta_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_sc_delta_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_s_vel": [
        "Variance for the assumed Gaussian velocity distributions used in collection and riming processes. Not used.",
        r"similar to $ \sigma_i $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [],
    ],
    "dhail_a_vel": [
        "Coefficient for particle velocity",
        r"$ \alpha $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        [],
    ],
    "dhail_b_vel": [
        "Exponent for particle velocity",
        r"$ \beta $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        [],
    ],
    "dhail_rho_v": [
        "Coefficient used in density correction for the increased terminal fall velocity with decreasing air density",
        r"$ \Big( \frac{\rho_0}{\rho} \Big)^\gamma $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "dependent",
        [],
    ],
    "dhail_c_z": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dhail_sc_coll_n": [
        "Coefficient in graupel self collection and cloud riming. Not tracked.",
        r"combination of $ \vartheta $ and $ \delta $",
        r"\citeA{seifert_two-moment_2006}, Eq. (62)",
        "independent",
        [],
    ],
    "dhail_cmu0": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_cmu1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_cmu2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_cmu3": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_cmu4": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_cmu5": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_alpha": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_beta": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_gamma": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_nu": [
        r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \mu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        [],
    ],
    "dhail_g1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_g2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_mu": [
        r"Shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \nu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        [],
    ],
    "dhail_nm1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_nm2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_nm3": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_q_crit_c": [
        "Threshold (mass density) in riming",
        r"Similar to $ \overline{D}_{g, 0} $ but for hail mass density",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [],
    ],
    "dhail_d_crit_c": [
        "Threshold (diameter) in riming",
        r"Similar to $ \overline{D}_{g, 0} $ but for hail",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [],
    ],
    "dhail_ecoll_c": [
        "Riming coefficient. Maximum collision efficiency with cloud droplets",
        r"Similar to $ \overline{E}_{g, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [],
    ],
    "dhail_cap": [
        "Coefficient for capacity of particle. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dhail_a_ven": [
        "Used to calculate the constant for average ventilation",
        "$ a_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "independent",
        [],
    ],
    "dhail_b_ven": [
        "Used to calculate the constant for average ventilation. Not tracked",
        "$ b_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "independent",
        [],
    ],
    "dhail_c_s": [
        "Coefficient in evaporation and vapor deposition for capacity",
        "similar to $ c_g^{-1} $",
        r"\citeA{seifert_two-moment_2006}",
        "dependent",
        [],
    ],
    "dhail_a_f": [
        "Constant for average ventilation. Used in melting and ice-vapor processes",
        r"$ \overline{a}_{\mathrm{vent}, 1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "dependent",
        [],
    ],
    "dhail_b_f": [
        "Coefficient for average ventilation",
        r"$ \overline{b}_{\mathrm{vent}, 1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "dependent",
        [],
    ],
    "dhail_alfa_n": [
        "Sedimentation velocity coefficient",
        r"Similar to $ \alpha_e \Gamma( \frac{k+\nu_e+\beta_e+1}{\mu_e}) "
        r"\Big( \Gamma( \frac{\nu_e+2}{\mu_e}) \Big)^{-1} $ but for hail particle number",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        [],
    ],
    "dhail_alfa_q": [
        "Sedimentation velocity coefficient",
        r"Similar to $ \alpha_e \Gamma( \frac{k+\nu_e+\beta_e+1}{\mu_e}) "
        r"\Big( \Gamma( \frac{\nu_e+2}{\mu_e}) \Big)^{-1} $ but for hail",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        [],
    ],
    "dhail_lambda": [
        "Sedimentation velocity coefficient",
        r"Similar to $ \Gamma \frac{\nu_e+1}{\mu_e} \Big( \frac{\nu_e+2}{\mu_e} \Big)^{-1} $ but for hail",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        [],
    ],
    "dhail_vsedi_min": [
        "Minimum sedimentation velocity parameter",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        [],
    ],
    "dhail_vsedi_max": [
        "Maximum sedimentation velocity parameter",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        [],
    ],
    "dice_a_geo": [
        "Coefficient for diameter size calculation",
        "$ a $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [
            "ice_self_collection",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
            "ice_riming",
        ],
    ],
    "dice_b_geo": [
        "Exponent for diameter size calculation",
        "$ b $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [
            "ice_self_collection",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
            "ice_riming",
        ],
    ],
    "dice_min_x": [
        "Minimum size of the particle used after the microphysics",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "ice_nuc_hom",
            "ice_activation_phillips",
            "ice_self_collection",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
            "ice_riming",
            "snow_riming",
            "particle_cloud_riming",
            "particle_rain_riming",
            "ice_melting",
        ],
    ],
    "dice_min_x_act": [
        "Minimum size of particle for ice activation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["ice_activation_phillips"],
    ],
    "dice_min_x_nuc_homo": [
        "Minimum size of particle for homogenous nucleation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["ice_nuc_hom"],
    ],
    "dice_min_x_nuc_hetero": [
        "Minimum size of particle for heterogeneous nucleation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [],
    ],
    "dice_min_x_melt": [
        "Minimum size of particle for melting",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["ice_melting"],
    ],
    "dice_min_x_evap": [
        "Minimum size of particle for evaporation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["evaporation"],
    ],
    "dice_min_x_freezing": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dice_min_x_depo": [
        "Minimum size of particle for vapor deposition",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["vapor_dep_relaxation"],
    ],
    "dice_min_x_collision": [
        "Minimum size of particle for ice-ice collision",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["ice_self_collection"],
    ],
    "dice_min_x_collection": [
        "Minimum size of particle for different collision processes",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["particle_collection"],
    ],
    "dice_min_x_conversion": [
        "Minimum size of particle for conversion processes",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["ice_riming"],
    ],
    "dice_min_x_sedimentation": [
        "Minimum size of particle for sedimentation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_min_x_riming": [
        "Minimum size of particle for riming",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "riming_cloud_core",
            "riming_rain_core",
            "ice_riming",
            "snow_riming",
            "particle_cloud_riming",
            "particle_rain_riming",
        ],
    ],
    "dice_max_x": [
        "Maximum size of particle",
        r"similar to $ \overline{x}_{r, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "ice_nuc_hom",
            "ice_self_collection",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
            "ice_riming",
            "ice_melting",
        ],
    ],
    "dice_sc_theta_q": [
        "Coefficient for collision process ice to snow",
        r"Special case of $ \vartheta $",
        r"\citeA{seifert_two-moment_2006}, Eq. (61)",
        "dependent",
        ["ice_self_collection"],
    ],
    "dice_sc_delta_q": [
        "Coefficient for collision process ice to snow",
        r"Special case of $ \delta $",
        r"\citeA{seifert_two-moment_2006}, Eq. (61)",
        "dependent",
        ["ice_self_collection"],
    ],
    "dice_sc_theta_n": [
        "Coefficient for collision process ice to snow",
        r"Special case of $ \vartheta $",
        r"\citeA{seifert_two-moment_2006}, Eq. (62)",
        "dependent",
        ["ice_self_collection"],
    ],
    "dice_sc_delta_n": [
        "Coefficient for collision process ice to snow",
        r"Special case of $ \delta $",
        r"\citeA{seifert_two-moment_2006}, Eq. (62)",
        "dependent",
        ["ice_self_collection"],
    ],
    "dice_s_vel": [
        "Variance for the assumed Gaussian velocity distributions used in collection and riming processes",
        r"$ \sigma_i $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [
            "ice_self_collection",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
        ],
    ],
    "dice_a_vel": [
        "Coefficient for particle velocity",
        r"$ \alpha $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        [
            "ice_self_collection",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
        ],
    ],
    "dice_b_vel": [
        "Exponent for particle velocity",
        r"$ \beta $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        [
            "ice_self_collection",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
        ],
    ],
    "dice_rho_v": [
        "Coefficient used in density correction for the increased terminal fall velocity with decreasing air density",
        r"$ \Big( \frac{\rho_0}{\rho} \Big)^\gamma $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "dependent",
        [
            "ice_self_collection",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
        ],
    ],
    "dice_c_z": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dice_sc_coll_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dice_cmu0": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_cmu1": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_cmu2": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_cmu3": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_cmu4": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_cmu5": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dice_alpha": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_beta": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_gamma": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_nu": [
        r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \mu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_g1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dice_g2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dice_mu": [
        r"Shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \nu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        [],
    ],
    "dice_nm1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dice_nm2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dice_nm3": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dice_q_crit_c": [
        "Threshold (mass density) in riming",
        r"Similar to $ \overline{D}_{i, 0} $ but for mass density",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["riming_cloud_core"],
    ],
    "dice_d_crit_c": [
        "Threshold (diameter) in riming",
        r"$ \overline{D}_{i, 0} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["riming_cloud_core"],
    ],
    "dice_ecoll_c": [
        "Riming coefficient. Maximum collision efficiency with cloud droplets",
        r"$ \overline{E}_{i, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["riming_cloud_core"],
    ],
    "dice_cap": [
        "Coefficient for capacity of particle. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dice_a_ven": [
        "Used to calculate the constant for average ventilation",
        "$ a_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "independent",
        ["vapor_dep_relaxation"],
    ],
    "dice_b_ven": [
        "Used to calculate the constant for average ventilation. Not tracked",
        "$ b_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "independent",
        [],
    ],
    "dice_c_s": [
        "Coefficient in evaporation and vapor deposition for capacity",
        "$ c_g^{-1} $",
        r"\citeA{seifert_two-moment_2006}",
        "dependent",
        ["evaporation", "vapor_dep_relaxation"],
    ],
    "dice_a_f": [
        "Constant for average ventilation. Used in melting and ice-vapor processes",
        r"$ \overline{a}_{\mathrm{vent}, 1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "dependent",
        ["evaporation", "vapor_dep_relaxation"],
    ],
    "dice_b_f": [
        "Coefficient for average ventilation",
        r"$ \overline{b}_{\mathrm{vent}, 1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "dependent",
        ["evaporation", "vapor_dep_relaxation"],
    ],
    "dice_alfa_n": [
        "Sedimentation velocity coefficient",
        r"Similar to $ \alpha_i \Gamma( \frac{k+\nu_i+\beta_i+1}{\mu_i}) "
        r"\Big( \Gamma( \frac{\nu_i+2}{\mu_i}) \Big)^{-1} $ but for particle number",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dice_alfa_q": [
        "Sedimentation velocity coefficient",
        r"$ \alpha_i \Gamma( \frac{k+\nu_i+\beta_i+1}{\mu_i}) \Big( \Gamma( \frac{\nu_i+2}{\mu_i}) \Big)^{-1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dice_lambda": [
        "Sedimentation velocity coefficient",
        r"$ \Gamma \frac{\nu_i+1}{\mu_i} \Big( \frac{\nu_i+2}{\mu_i} \Big)^{-1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dice_vsedi_min": [
        "Minimum sedimentation velocity parameter. Not tracked.",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dice_vsedi_max": [
        "Maximum sedimentation velocity parameter",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_a_geo": [
        "Coefficient for diameter size calculation",
        "$ a $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [
            "ice_self_collection",
            "snow_self_collection",
            "snow_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
            "snow_riming",
        ],
    ],
    "dsnow_b_geo": [
        "Exponent for diameter size calculation",
        "$ b $",
        r"\citeA{seifert_two-moment_2006} Eq. (32)",
        "independent",
        [
            "ice_self_collection",
            "snow_self_collection",
            "snow_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
            "snow_riming",
        ],
    ],
    "dsnow_min_x": [
        "Minimum size of the particle used after the microphysics",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "snow_self_collection",
            "snow_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
            "snow_riming",
        ],
    ],
    "dsnow_min_x_act": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_min_x_nuc_homo": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_min_x_nuc_hetero": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_min_x_melt": [
        "Minimum size of particle for melting",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["snow_melting"],
    ],
    "dsnow_min_x_evap": [
        "Minimum size of particle for evaporation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["evaporation"],
    ],
    "dsnow_min_x_freezing": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_min_x_depo": [
        "Minimum size of particle for vapor deposition",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["vapor_dep_relaxation"],
    ],
    "dsnow_min_x_collision": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_min_x_collection": [
        "Minimum size of particle for different collision processes",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["snow_self_collection", "particle_collection"],
    ],
    "dsnow_min_x_conversion": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_min_x_sedimentation": [
        "Minimum size of particle for sedimentation",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_min_x_riming": [
        "Minimum size of particle for riming",
        r"similar to $ \overline{x}_{r, \mathrm{min}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        ["riming_cloud_core", "riming_rain_core", "snow_riming"],
    ],
    "dsnow_max_x": [
        "Maximum size of particle",
        r"similar to $ \overline{x}_{r, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}, Eqs. (94), (97)",
        "independent",
        [
            "snow_self_collection",
            "snow_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
            "snow_riming",
        ],
    ],
    "dsnow_sc_theta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_sc_delta_q": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_sc_theta_n": [
        "Coefficient for collision process snow with snow",
        r"Special case of $ \vartheta $",
        r"\citeA{seifert_two-moment_2006}, Eq. (62)",
        "dependent",
        ["snow_self_collection"],
    ],
    "dsnow_sc_delta_n": [
        "Coefficient for collision process snow with snow",
        r"Special case of $ \delta $",
        r"\citeA{seifert_two-moment_2006}, Eq. (62)",
        "dependent",
        ["snow_self_collection"],
    ],
    "dsnow_s_vel": [
        "Variance for the assumed Gaussian velocity distributions used in collection and riming processes",
        r"$ \sigma_s $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        [
            "snow_self_collection",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
        ],
    ],
    "dsnow_a_vel": [
        "Coefficient for particle velocity",
        r"$ \alpha $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        [
            "snow_self_collection",
            "snow_melting",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
        ],
    ],
    "dsnow_b_vel": [
        "Exponent for particle velocity",
        r"$ \beta $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "independent",
        [
            "snow_self_collection",
            "snow_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
        ],
    ],
    "dsnow_rho_v": [
        "Coefficient used in density correction for the increased terminal fall velocity with decreasing air density",
        r"$ \Big( \frac{\rho_0}{\rho} \Big)^\gamma $",
        r"\citeA{seifert_two-moment_2006} Eq. (33)",
        "dependent",
        [
            "snow_self_collection",
            "snow_melting",
            "sedimentation_explicit",
            "evaporation",
            "vapor_dep_relaxation",
            "particle_collection",
            "riming_cloud_core",
            "riming_rain_core",
        ],
    ],
    "dsnow_c_z": [
        "Not used",
        "",
        "",
        "dependent",
        [],
    ],
    "dsnow_sc_coll_n": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_cmu0": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_cmu1": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_cmu2": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_cmu3": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_cmu4": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_cmu5": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_alpha": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_beta": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_gamma": [
        "Not used",
        "",
        "",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_nu": [
        r"Parameter to calculate the shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \mu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_g1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_g2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_mu": [
        r"Shape parameter of the generalized $ \Gamma $-distribution. Not tracked",
        r"$ \nu $",
        r"\citeA{seifert_two-moment_2006}, Eq. (79)",
        "independent",
        [],
    ],
    "dsnow_nm1": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_nm2": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_nm3": [
        "Not used",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_q_crit_c": [
        "Threshold (mass density) in riming",
        r"Similar to $ \overline{D}_{s, 0} $ but for mass density",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["riming_cloud_core"],
    ],
    "dsnow_d_crit_c": [
        "Threshold (diameter) in riming",
        r"$ \overline{D}_{s, 0} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["riming_cloud_core"],
    ],
    "dsnow_ecoll_c": [
        "Riming coefficient. Maximum collision efficiency with cloud droplets",
        r"$ \overline{E}_{s, \mathrm{max}} $",
        r"\citeA{seifert_two-moment_2006}",
        "independent",
        ["riming_cloud_core"],
    ],
    "dsnow_cap": [
        "Coefficient for capacity of particle. Not tracked",
        "",
        "",
        "independent",
        [],
    ],
    "dsnow_a_ven": [
        "Used to calculate the constant for average ventilation",
        "$ a_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "independent",
        ["vapor_dep_relaxation"],
    ],
    "dsnow_b_ven": [
        "Used to calculate the constant for average ventilation. Not tracked",
        "$ b_v $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "independent",
        [],
    ],
    "dsnow_c_s": [
        "Coefficient in evaporation and vapor deposition for capacity",
        "$ c_g^{-1} $",
        r"\citeA{seifert_two-moment_2006}",
        "dependent",
        ["evaporation", "vapor_dep_relaxation"],
    ],
    "dsnow_a_f": [
        "Constant for average ventilation. Used in melting and ice-vapor processes",
        r"$ \overline{a}_{\mathrm{vent}, 1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (88)",
        "dependent",
        ["snow_melting", "evaporation", "vapor_dep_relaxation"],
    ],
    "dsnow_b_f": [
        "Coefficient for average ventilation",
        r"$ \overline{b}_{\mathrm{vent}, 1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (89)",
        "dependent",
        ["snow_melting", "evaporation", "vapor_dep_relaxation"],
    ],
    "dsnow_alfa_n": [
        "Sedimentation velocity coefficient",
        r"Similar to $ \alpha_s \Gamma( \frac{k+\nu_s+\beta_s+1}{\mu_s}) "
        r"\Big( \Gamma( \frac{\nu_s+2}{\mu_s}) \Big)^{-1} $ but for particle number",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dsnow_alfa_q": [
        "Sedimentation velocity coefficient",
        r"$ \alpha_s \Gamma( \frac{k+\nu_s+\beta_s+1}{\mu_s}) "
        r"\Big( \Gamma( \frac{\nu_s+2}{\mu_s}) \Big)^{-1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dsnow_lambda": [
        "Sedimentation velocity coefficient",
        r"$ \Gamma \frac{\nu_s+1}{\mu_s} \Big( \frac{\nu_s+2}{\mu_s} \Big)^{-1} $",
        r"\citeA{seifert_two-moment_2006}, Eq. (78)",
        "dependent",
        ["sedimentation_explicit"],
    ],
    "dsnow_vsedi_min": [
        "Minimum sedimentation velocity parameter",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        ["sedimentation_explicit"],
    ],
    "dsnow_vsedi_max": [
        "Maximum sedimentation velocity parameter",
        "-",
        "Introduced to avoid numerical issues",
        "independent",
        ["sedimentation_explicit"],
    ],
}

# A dictionary of used values for each parameter
in_params_value_dic = {
    "da_1": 1.0e-3,
    "da_2": r"$ 1.72 / (R_a^{7/8}) $",
    "de_1": r"$ 1 / \sqrt{R_a} $",
    "de_2": r"$ 9.1 / R_a^{11/16} $",
    "dd": "see Equation ..",
    "dN_c": 50,
    "dgamma": 1,
    "dbeta_c": 1,
    "dbeta_r": 7.0 / 8.0,
    "ddelta1": 0.5,
    "ddelta2": 11.0 / 16.0,
    "dzeta": 9.0 / 8.0,
    "drain_gfak": 1,
    "dcloud_k_au": " see Equation .. ",
    "dcloud_k_sc": " see Equation .. ",
    "dkc_autocon": 9.44e9,
    "dinv_z": 1.0 / 250.0,
    "dw": "",
    "dq_crit_i": "",
    "dD_crit_i": "",
    "dD_conv_i": "",
    "dq_crit_r": "",
    "dD_crit_r": "",
    "dq_crit_fr": "",
    "dq_crit_c": "",
    "dq_crit": "",
    "dD_conv_sg": "",
    "dD_conv_ig": "",
    "dx_conv": "",
    "dD_crit_c": "",
    "dD_coll_c": "",
    "dparcel_height": "",
    "dT_nuc": "",
    "dT_freeze": 273.15,
    "dT_f": 233.0,
    "dD_eq": "",
    "drho_w": 1000,
    "drho_0": 1.225,
    "drho_vel": 0.4,
    "drho_vel_c": 1,
    "drho_ice": 916.7,
    "dM_w": "",
    "dM_a": "",
    "dR_universal": "",
    "dEpsilon": "",
    "dgravity_acc": "",
    "dR_a": "",
    "dR_v": "",
    "da_v": 0.78,
    "db_v": 0.308,
    "da_prime": 9.65,
    "db_prime": 9.8,
    "dc_prime": 600,
    "dK_T": 2.4e-2,
    "dL_wd": 2.5008e6,
    "dL_ed": 2.8345e6,
    "dD_v": 2.22e-5,
    "decoll_min": 0.01,
    "decoll_gg": 0.10,
    "decoll_gg_wet": 0.4,
    "dkin_visc_air": 1.5e-5,
    "dalpha_spacefilling": 0.01,
    "dC_mult": 3.5e8,
    "dT_mult_min": "",
    "dT_mult_max": "",
    "dT_mult_opt": "",
    "dconst0": "",
    "dconst3": "",
    "dconst4": "",
    "dconst5": "",
    "dD_rainfrz_gh": 1.25e-3,
    "dD_rainfrz_ig": "0.5e-3",
    "ddv0": "",
    "dp_sat_melt": 6.1078e2,
    "dcp": 1004.64,
    "dk_b": "",
    "da_HET": 0.65,
    "db_HET": 200,
    "dN_sc": "",
    "dn_f": "",
    "dN_avo": "",
    "dgrav": "",
    "damd": "",
    "damw": "",
    "dna_dust": "160e4 for nuc_type == 7 or 6 or 5, else 70e4",
    "dna_soot": "25e6 for nuc_type == 7 or 5, 30e6 for nuc_type == 6, else 0",
    "dna_orga": "30e6 for nuc_type == 7 or 5, else 0",
    "dni_het_max": "",
    "dni_hom_max": "",
    "da_dep": "",
    "db_dep": "",
    "dc_dep": "",
    "dd_dep": "",
    "dnim_imm": "",
    "dnin_dep": "",
    "dalf_imm": "",
    "dbet_dep": "",
    "dbet_imm": "",
    "da_ccn_1": "",
    "da_ccn_2": "",
    "da_ccn_3": "",
    "da_ccn_4": "",
    "db_ccn_1": "",
    "db_ccn_2": "",
    "db_ccn_3": "",
    "db_ccn_4": "",
    "dc_ccn_1": "",
    "dc_ccn_2": "",
    "dc_ccn_3": "",
    "dc_ccn_4": "",
    "dd_ccn_1": 287736034.13,
    "dd_ccn_2": 0.6258809883,
    "dd_ccn_3": 0.8907491812,
    "dd_ccn_4": 360848977.55,
    "dr_const": "",
    "dr1_const": "",
    "dcv": "",
    "dp_sat_const_a": "",
    "dp_sat_ice_const_a": "",
    "dp_sat_const_b": "",
    "dp_sat_ice_const_b": "",
    "dp_sat_low_temp": "",
    "dT_sat_low_temp": "",
    "dalpha_depo": "",
    "dr_0": "",
    "dk_1_conv": 400,
    "dk_2_conv": 0.7,
    "dk_1_accr": 5.0e-4,
    "dk_r": 5.78,
    # Cloud
    "dcloud_nu": 1,
    "dcloud_mu": 1,
    "dcloud_max_x": 2.6e-10,
    "dcloud_min_x": 4.2e-15,
    "dcloud_min_x_act": 4.2e-15,
    "dcloud_min_x_nuc_homo": 4.2e-15,
    "dcloud_min_x_nuc_hetero": 4.2e-15,
    "dcloud_min_x_melt": 4.2e-15,
    "dcloud_min_x_evap": 4.2e-15,
    "dcloud_min_x_freezing": 4.2e-15,
    "dcloud_min_x_depo": 4.2e-15,
    "dcloud_min_x_collision": 4.2e-15,
    "dcloud_min_x_collection": 4.2e-15,
    "dcloud_min_x_conversion": 4.2e-15,
    "dcloud_min_x_sedimentation": 4.2e-15,
    "dcloud_min_x_riming": 4.2e-15,
    "dcloud_a_geo": 1.24e-1,
    "dcloud_b_geo": 0.333333,
    "dcloud_a_vel": 3.75e5,
    "dcloud_b_vel": 0.666667,
    "dcloud_a_ven": 0.78,
    "dcloud_b_ven": 0.308,
    "dcloud_cap": 2.0,
    "dcloud_vsedi_max": 1.0,
    "dcloud_vsedi_min": 0.0,
    "dcloud_c_s": r"$1/" + parse_word("dcloud_cap")[1::],
    "dcloud_a_f": "see Equation ..",
    "dcloud_b_f": "see Equation ..",
    "dcloud_c_z": "see Equation ..",
    "dcloud_rho_v": "Density correction for terminal fall velocity",
    # Rain
    "drain_nu": -2.0 / 3.0,  # SB: -2/3 COSMO: 0.0
    "drain_mu": 1.0 / 3.0,  # SB: 1/3 COMSO: 1.0/3.0
    "drain_max_x": 3.0e-6,
    "drain_min_x": 2.6e-10,
    "drain_min_x_act": 2.6e-10,
    "drain_min_x_nuc_homo": 2.6e-10,
    "drain_min_x_nuc_hetero": 2.6e-10,
    "drain_min_x_melt": 2.6e-10,
    "drain_min_x_evap": 2.6e-10,
    "drain_min_x_freezing": 2.6e-10,
    "drain_min_x_depo": 2.6e-10,
    "drain_min_x_collision": 2.6e-10,
    "drain_min_x_collection": 2.6e-10,
    "drain_min_x_conversion": 2.6e-10,
    "drain_min_x_sedimentation": 2.6e-10,
    "drain_min_x_riming": 2.6e-10,
    "drain_a_geo": 1.24e-1,
    "drain_b_geo": 0.333333,
    "drain_a_vel": 114.0137,
    "drain_b_vel": 0.234370,
    "drain_cap": 2.0,
    "drain_alpha": 9.292,
    "drain_beta": 9.623,
    "drain_gamma": 6.222e2,
    "drain_cmu0": 6.0,
    "drain_cmu1": 3.0e1,
    "drain_cmu2": 1.0e3,
    "drain_cmu3": 1.1e-3,
    "drain_cmu4": 1.0,
    "drain_cmu5": 2.0,
    "drain_vsedi_max": 20.0,
    "drain_vsedi_min": 0.1,
    "drain_c_s": r"$1/" + parse_word("drain_cap")[1::],
    "drain_a_f": "see Equation ..",
    "drain_b_f": "see Equation ..",
    "drain_c_z": "see Equation ..",
    "drain_nm1": "see Equation ..",
    "drain_nm2": "see Equation ..",
    "drain_nm3": "see Equation ..",
    "drain_g1": 1,
    "drain_g2": 6,
    "drain_rho_v": "Density correction for terminal fall velocity",
    # Graupel
    "dgraupel_nu": 1.0,  # SB
    "dgraupel_mu": 1.0 / 3.0,  # SB
    "dgraupel_max_x": 5.0e-4,
    "dgraupel_min_x": 1.0e-9,
    "dgraupel_min_x_act": 1.0e-9,
    "dgraupel_min_x_nuc_homo": 1.0e-9,
    "dgraupel_min_x_nuc_hetero": 1.0e-9,
    "dgraupel_min_x_melt": 1.0e-9,
    "dgraupel_min_x_evap": 1.0e-9,
    "dgraupel_min_x_freezing": 1.0e-9,
    "dgraupel_min_x_depo": 1.0e-9,
    "dgraupel_min_x_collision": 1.0e-9,
    "dgraupel_min_x_collection": 1.0e-9,
    "dgraupel_min_x_conversion": 1.0e-9,
    "dgraupel_min_x_sedimentation": 1.0e-9,
    "dgraupel_min_x_riming": 1.0e-9,
    "dgraupel_a_geo": 1.42e-1,
    "dgraupel_b_geo": 0.314,
    "dgraupel_a_vel": 86.89371,
    "dgraupel_b_vel": 0.268325,
    "dgraupel_a_ven": 0.78,
    "dgraupel_b_ven": 0.308,
    "dgraupel_cap": 2.0,
    "dgraupel_vsedi_max": 30.0,
    "dgraupel_vsedi_min": 0.1,
    "dgraupel_d_crit_c": 100.0e-6,
    "dgraupel_q_crit_c": 1.0e-6,
    "dgraupel_s_vel": 0.0,
    "dgraupel_ecoll_c": 1.0,
    "dgraupel_c_s": r"$1/" + parse_word("dgraupel_cap")[1::],
    "dgraupel_nm1": "see Equation ..",
    "dgraupel_nm2": "see Equation ..",
    "dgraupel_a_f": "see Equation ..",
    "dgraupel_b_f": "see Equation ..",
    "dgraupel_c_z": "see Equation ..",
    "dgraupel_sc_coll_n": "see Equation ..",
    "dgraupel_g1": 120,
    "dgraupel_g2": 40320,
    "dgraupel_lambda": "see Equation ..",
    "dgraupel_alfa_n": "see Equation ..",
    "dgraupel_alfa_q": "see Equation ..",
    "dgraupel_rho_v": "Density correction for terminal fall velocity",
    # Hail
    "dhail_nu": 1.0,
    "dhail_mu": 1.0 / 3.0,
    "dhail_max_x": 5.0e-4,
    "dhail_min_x": 2.6e-9,
    "dhail_min_x_act": 2.6e-9,
    "dhail_min_x_nuc_homo": 2.6e-9,
    "dhail_min_x_nuc_hetero": 2.6e-9,
    "dhail_min_x_melt": 2.6e-9,
    "dhail_min_x_evap": 2.6e-9,
    "dhail_min_x_freezing": 2.6e-9,
    "dhail_min_x_depo": 2.6e-9,
    "dhail_min_x_collision": 2.6e-9,
    "dhail_min_x_collection": 2.6e-9,
    "dhail_min_x_conversion": 2.6e-9,
    "dhail_min_x_sedimentation": 2.6e-9,
    "dhail_min_x_riming": 2.6e-9,
    "dhail_a_geo": 0.1366,
    "dhail_b_geo": 1.0 / 3.0,
    "dhail_a_vel": 39.3,
    "dhail_b_vel": 0.166667,
    "dhail_a_ven": 0.78,
    "dhail_b_ven": 0.308,
    "dhail_cap": 2.0,
    "dhail_vsedi_max": 30.0,
    "dhail_vsedi_min": 0.1,
    "dhail_sc_coll_n": 1.0,
    "dhail_d_crit_c": 100.0e-6,
    "dhail_q_crit_c": 1.0e-6,
    "dhail_s_vel": 0.0,
    "dhail_ecoll_c": 1.0,
    "dhail_c_s": r"$1/" + parse_word("dhail_cap")[1::],
    "dhail_a_f": "see Equation ..",
    "dhail_b_f": "see Equation ..",
    "dhail_c_z": "see Equation ..",
    "dhail_lambda": "see Equation ..",
    "dhail_alfa_n": "see Equation ..",
    "dhail_alfa_q": "see Equation ..",
    "dhail_rho_v": "Density correction for terminal fall velocity",
    # Ice
    "dice_nu": 1.0,
    "dice_mu": 1.0 / 3.0,
    "dice_max_x": 1.0e-5,
    "dice_min_x": 1.0e-12,
    "dice_min_x_act": 1.0e-12,
    "dice_min_x_nuc_homo": 1.0e-12,
    "dice_min_x_nuc_hetero": 1.0e-12,
    "dice_min_x_melt": 1.0e-12,
    "dice_min_x_evap": 1.0e-12,
    "dice_min_x_freezing": 1.0e-12,
    "dice_min_x_depo": 1.0e-12,
    "dice_min_x_collision": 1.0e-12,
    "dice_min_x_collection": 1.0e-12,
    "dice_min_x_conversion": 1.0e-12,
    "dice_min_x_sedimentation": 1.0e-12,
    "dice_min_x_riming": 1.0e-12,
    "dice_a_geo": 0.835,
    "dice_b_geo": 0.39,
    "dice_a_vel": 2.77e1,
    "dice_b_vel": 0.21579,
    "dice_a_ven": 0.78,
    "dice_b_ven": 0.308,
    "dice_cap": 2.0,
    "dice_vsedi_max": 3.0,
    "dice_vsedi_min": 0.0,
    "dice_sc_coll_n": 0.8,
    "dice_d_crit_c": 150.0e-6,
    "dice_q_crit_c": 1.0e-5,
    "dice_s_vel": 0.05,
    "dice_ecoll_c": 0.80,
    "dice_c_s": r"$1/" + parse_word("dice_cap")[1::],
    "dice_a_f": "see Equation ..",
    "dice_b_f": "see Equation ..",
    "dice_c_z": "see Equation ..",
    "dice_sc_delta_n": "see Equation ..",
    "dice_sc_delta_q": "see Equation ..",
    "dice_sc_theta_n": "see Equation ..",
    "dice_sc_theta_q": "see Equation ..",
    "dice_lambda": "see Equation ..",
    "dice_alfa_n": "see Equation ..",
    "dice_alfa_q": "see Equation ..",
    "dice_rho_v": "Density correction for terminal fall velocity",
    # Snow
    "dsnow_nu": 1.0,  # COSMO: 0.0, SB 1.0
    "dsnow_mu": 1.0 / 3.0,  # COSMO 0.5, SB: 1.0/3.0
    "dsnow_max_x": 2.0e-5,
    "dsnow_min_x": 1.0e-10,
    "dsnow_min_x_act": 1.0e-10,
    "dsnow_min_x_nuc_homo": 1.0e-10,
    "dsnow_min_x_nuc_hetero": 1.0e-10,
    "dsnow_min_x_melt": 1.0e-10,
    "dsnow_min_x_evap": 1.0e-10,
    "dsnow_min_x_freezing": 1.0e-10,
    "dsnow_min_x_depo": 1.0e-10,
    "dsnow_min_x_collision": 1.0e-10,
    "dsnow_min_x_collection": 1.0e-10,
    "dsnow_min_x_conversion": 1.0e-10,
    "dsnow_min_x_sedimentation": 1.0e-10,
    "dsnow_min_x_riming": 1.0e-10,
    "dsnow_a_geo": 2.4,
    "dsnow_b_geo": 0.455,
    "dsnow_a_vel": 8.8,
    "dsnow_b_vel": 0.15,
    "dsnow_a_ven": 0.78,
    "dsnow_b_ven": 0.308,
    "dsnow_cap": 2.0,
    "dsnow_vsedi_max": 3.0,
    "dsnow_vsedi_min": 0.1,
    "dsnow_sc_coll_n": 0.8,
    "dsnow_d_crit_c": 150.0e-6,
    "dsnow_q_crit_c": 1.0e-5,
    "dsnow_s_vel": 0.25,
    "dsnow_ecoll_c": 0.80,
    "dsnow_c_s": r"$1/" + parse_word("dsnow_cap")[1::],
    "dsnow_a_f": "see Equation ..",
    "dsnow_b_f": "see Equation ..",
    "dsnow_c_z": "see Equation ..",
    "dsnow_sc_delta_n": "see Equation ..",
    "dsnow_sc_theta_n": "see Equation ..",
    "dsnow_lambda": "see Equation ..",
    "dsnow_alfa_n": "see Equation ..",
    "dsnow_alfa_q": "see Equation ..",
    "dsnow_rho_v": "Density correction for terminal fall velocity",
}


def set_size(beamer=True, scale=None):
    """
    Set some options to use latex.

    Parameters
    ----------
    beamer : bool
        Beamer is used for bigger texts.
    scale : float
        Can be used to scale the fontsize.
    """
    if beamer:
        mpl.rcParams.update(
            {
                "text.usetex": False,
                "font.family": "serif",
                "axes.labelsize": int(20 * scale),
                "font.size": int(20 * scale),
                "legend.fontsize": int(16 * scale),
                "xtick.labelsize": int(12 * scale),
                "ytick.labelsize": int(16 * scale),
            }
        )
    else:
        mpl.rcParams.update(
            {
                "text.usetex": False,
                "font.family": "serif",
                "axes.labelsize": int(10 * scale),
                "font.size": int(10 * scale),
                "legend.fontsize": int(10 * scale),
                "xtick.labelsize": int(8 * scale),
                "ytick.labelsize": int(8 * scale),
            }
        )


def replace_cites(str_to_replace):
    r"""

    Parameters
    ----------
    str_to_replace : string
        String with occurrences of the form of \citeA{...}
    Returns
    -------
    Replaced the \citeA with the name of the citation.
    """
    cite_dic = {
        r"\citeA{seifert_two-moment_2006}": "Seifert and Beheng (2006a)",
        r"\citeA{seifert_parameterization_2008}": "Seifert (2008)",
        r"\citeA{hande_parameterizing_2016}": "Hande et al. (2016)",
        r"\citeA{phillips_empirical_2008}": "Phillips et al. (2008)",
        r"\citeA{hallett_production_1974}": "Hallet and Mossop (1974)",
        r"\citeA{seifert_parametrisierung_2002}": "Seifert (2002)",
        r"\citeA{seifert_double-moment_2001}": "Seifert and Beheng (2001)",
        r"\citeA{karcher_physically_2006}": r"Kärcher et al. (2006)",
    }
    for key, replacement in cite_dic.items():
        if key in str_to_replace:
            str_to_replace = str_to_replace.replace(key, replacement)
    return str_to_replace


def cites_to_citep(str_to_replace):
    """
    Replace some strings that might be from an autogenerated table such as from get_process_desc_table(..).

    Parameters
    ----------
    str_to_replace : string
        String with occurrences of the form of https://doi.org/...

    Returns
    -------
    String with citep instead of doi or similar strings.
    """
    cite_dic = {
        "Hande et al 2016": r"\citep{Hande2016}",
        r"https://doi.org/10.5194/acp-16-12059-2016": "",
        r"Seifert \& Beheng (2006)": r"\citep{Seifert2006b}",
        r"Seifert \& Beheng; 2006": r"\citep{Seifert2006b}",
        r"Seifert and Beheng (2006)": r"\citep{Seifert2006b}",
        r"10.1007/s00703-005-0112-4": "",
        "Rutledge, Steven A., and Peter V. Hobbs": r"\cite{Rutledge1984}",
        '"The mesoscale and microscale structure and organization of clouds '
        "and precipitation in midlatitude cyclones.": "",
        'XII: A diagnostic modeling study of precipitation development in narrow cold-frontal rainbands."': "",
        "Journal of Atmospheric Sciences 41.20 (1984): 2949-2972.": "",
        r"Ulrich Blahak": r"\citep{Blahak2008}",
        '"Towards a better representation of high density': "",
        "ice particles in a state-of-the-art two-moment": "",
        'bulk microphysical scheme."': "",
        "Proc. 15th Int. Conf. Clouds and Precip.,": "",
        "Cancun, Mexico. Vol. 20208. 2008": "",
        "Beheng, K. D.": r"\citep{Beheng1982}",
        '"Numerical study on the combined action of droplet coagulation,': "",
        'ice particle riming and the splintering process concerning maritime cumuli."': "",
        "Contrib. Atmos. Phys.;(Germany, Federal Republic of) 55.3 (1982).": "",
        "Seifert and Beheng (2008)": r"\citep{Seifert2008}",
        r"https://www.imk-tro.kit.edu/4437\_1388.php": "",
        "H. Morrison, J.A.Curry, V.I. Khvorostyanov": r"\citep{Morrison2005}",
        '"A New Double-Moment Microphysics Parameterization for Application in Cloud and': "",
        'Climate Models. Part 1: Description" by ': "",
        r"10.1175/2008JAS2586.1": "",
        "Seifert (2008)": r"\citep{Seifert2008}",
    }
    for key, replacement in cite_dic.items():
        str_to_replace = str_to_replace.replace(key, replacement)
    return str_to_replace


def get_value(param):
    """
    Given an input parameter, get the value used in the model.

    Parameters
    ----------
    param : string
        Key from in_params_value_dic (a model input parameter)

    Returns
    -------
    Formatted string with value of the model input parameter.
    """
    bracket = "}"
    if param in in_params_numeric_value_dic:
        if 1000 > in_params_numeric_value_dic[param] >= 1:
            if float(in_params_numeric_value_dic[param]).is_integer():
                return f"$ {in_params_numeric_value_dic[param]} $"
            new_str = f"$ {in_params_numeric_value_dic[param]:1.3E}{bracket} $".replace(
                "E", r"\mathrm{E}{"
            )
            new_str = (
                new_str.replace(r"0\mathrm{E}", r"\mathrm{E}")
                .replace(r"0\mathrm{E}", r"\mathrm{E}")
                .replace(r"0\mathrm{E}", r"\mathrm{E}")
            )
            if r".\mathrm{E}" in new_str:
                new_str = new_str.replace(r".\mathrm{E}", r"\mathrm{E}")
            if r"\mathrm{E}+00" in new_str:
                return new_str.replace(r"\mathrm{E}+00", "")
            return new_str
        if float(in_params_numeric_value_dic[param]).is_integer():
            return f"$ {in_params_numeric_value_dic[param]:1.0E}{bracket} $".replace(
                "E", r"\mathrm{E}{"
            )
        new_str = f"$ {in_params_numeric_value_dic[param]:1.3E}{bracket} $".replace(
            "E", r"\mathrm{E}{"
        )
        if in_params_numeric_value_dic[param] < 1:
            new_str = (
                new_str.replace(r"0\mathrm{E}", r"\mathrm{E}")
                .replace(r"0\mathrm{E}", r"\mathrm{E}")
                .replace(r"0\mathrm{E}", r"\mathrm{E}")
            )
            if r".\mathrm{E}" in new_str:
                new_str = new_str.replace(r".\mathrm{E}", r"\mathrm{E}")
        return new_str
    return ""


def get_unit(param, brackets=False):
    """
    Get the unit for a given parameter.

    Parameters
    ----------
    param : string
        Name of the parameter
    brackets : bool
        If true, add "[]" to output

    Returns
    -------
    string
        Unit string.
    """
    unit_dic = {
        "T": "K",
        "p": "Pa",
        "pressure": "Pa",
        "pressure_hPa": "hPa",
        "time_after_ascent": "s",
        "time_after_ascent_h": "h",
        "time": "s",
        "timestep": "s",
        "QR_OUT": r"kg/m³",
        "QS_OUT": r"kg/m³",
        "QI_OUT": r"kg/m³",
        "QG_OUT": r"kg/m³",
        "Q_total": r"kg/m³",
        "NR_OUT": r"$ \mathrm{m}^{-3} $",
        "QV": r"kg/m³",
        "QC": r"kg/m³",
        "QR": r"kg/m³",
        "QG": r"kg/m³",
        "QH": r"kg/m³",
        "QI": r"kg/m³",
        "QS": r"kg/m³",
        "Q_cold": r"kg/m³",
        "Q_liquid": r"kg/m³",
        "dinv_z": r"$ \mathrm{m}^{-1} $",
        "drain_vsedi_min": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "drain_vsedi_max": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dice_vsedi_min": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dice_vsedi_max": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dsnow_vsedi_min": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dsnow_vsedi_max": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dgraupel_vsedi_min": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dgraupel_vsedi_max": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dhail_vsedi_min": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dhail_vsedi_max": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dkc_autocon": r"$ \mathrm{m}\,\mathrm{kg}^{-2}\mathrm{s}^{-1} $",
        "dq_crit_i": r"$ \mathrm{kg}\,\mathrm{m}^{-3} $",
        "dD_crit_i": r"m",
        "dD_conv_i": r"m",
        "dq_crit_r": r"$ \mathrm{kg}\,\mathrm{m}^{-3} $",
        "dD_crit_r": r"m",
        "dq_crit_fr": r"$ \mathrm{kg}\,\mathrm{m}^{-3} $",
        "dD_conv_sg": r"m",
        "dD_conv_ig": r"m",
        "dx_conv": "kg",
        "dT_nuc": "K",
        "da_prime": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "db_prime": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "dc_prime": r"$ \mathrm{m}^{-1} $",
        "dkin_visc_air": r"$ \mathrm{m}^2\mathrm{s}^{-1} $",
        "dC_mult": r"$ \mathrm{kg}^{-1} $",
        "dT_mult_min": "K",
        "dT_mult_max": "K",
        "dconst0": r"$ \mathrm{m}^{-1} $",
        "dconst3": r"$ \mathrm{K}^{-1} $",
        "dconst4": r"$ \mathrm{K}^{-1} $",
        "dD_rainfrz_gh": "m",
        "dD_rainfrz_ig": "m",
        "ddv0": r"$ \mathrm{m}^2\mathrm{s}^{-1} $",
        "dp_sat_melt": "Pa",
        "dk_b": r"$ \mathrm{J}\,\mathrm{K}^{-1} $",
        "da_HET": r"$ \mathrm{K}^{-1} $",
        "db_HET": r"$ \mathrm{kg}^{-1}\mathrm{s}^{-1} $",
        "dN_avo": r"$ \mathrm{mol}^{-1} $",
        "dna_dust": r"$ \mathrm{m}^{-3} $",
        "dna_soot": r"$ \mathrm{m}^{-3} $",
        "dna_orga": r"$ \mathrm{m}^{-3} $",
        "dni_het_max": r"$ \mathrm{l}^{-1} $",
        "dni_hom_max": r"$ \mathrm{l}^{-1} $",
        "dr_0": "m",
        "dk_r": r"$ \mathrm{m}^3\mathrm{kg}^{-1}\mathrm{s}^{-1} $",
        "da_ccn_1": r"$ \mathrm{m}^{-3} $",
        "da_ccn_2": r"$ \log(\mathrm{m}\,\mathrm{s}^{-1})^{-1} $",
        "da_ccn_4": r"$ \mathrm{m}^{-3} $",
        "db_ccn_1": r"$ \mathrm{Pa}^{-1} $",
        "db_ccn_2": r"$ \mathrm{Pa}^{-1} $",
        "db_ccn_4": r"$ \mathrm{Pa}^{-1} $",
        "dd_ccn_1": r"$ \mathrm{m}^{-3} $",
        "dd_ccn_2": r"$ \log(\mathrm{m}\,\mathrm{s}^{-1})^{-1} $",
        "dd_ccn_4": r"$ \mathrm{m}^{-3} $",
        "drain_cmu2": r"$ \mathrm{m}^{-1} $",
        "drain_cmu3": r"$ \mathrm{m} $",
        "drain_alpha": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "drain_beta": r"$ \mathrm{m}\,\mathrm{s}^{-1} $",
        "drain_gamma": r"$ \mathrm{m}^{-1} $",
    }
    word = ""
    if param in unit_dic:
        if brackets:
            word = "[" + unit_dic[param] + "]"
        else:
            word = unit_dic[param]
    elif param.endswith("_min_x") or "_min_x_" in param:
        word = r"kg"
    elif "_a_geo" in param:
        particle = param[1:].replace("_a_geo", "")
        word = (
            r"$ \mathrm{m}\,\mathrm{kg}^{-\mathrm{geo}_{b, \mathrm{"
            + particle
            + r"}}} $"
        )
    elif "_a_vel" in param:
        particle = param[1:].replace("_a_vel", "")
        word = (
            r"$ \mathrm{m}\,\mathrm{s}^{-1}\mathrm{kg}^{-\mathrm{vel}_{b, \mathrm{"
            + particle
            + r"}}} $"
        )
    elif param.endswith("_max_x"):
        word = "kg"
    elif param.endswith("_q_crit_c"):
        word = r"$ \mathrm{kg}\,\mathrm{m}^{-3} $"
    elif param.endswith("_d_crit_c"):
        word = r"$ \mathrm{m} $"
    elif param.endswith("_s_vel"):
        word = r"$ \mathrm{m}\,\mathrm{s}^{-1} $"
    return word


# pylint: disable=too-many-locals
def top_phase_to_table(top, caption, label, parse=True):
    """
    Create either a latex table with model parameters and a description of those.
    Or differentiate between sensitivities for different model state variables and
    print the top parameters without description.
    Or differentiate between different model state variables and phases but also
    without description.

    Parameters
    ----------
    top : Dict
        A dictionary with keys the model state variable for a table with top
        parameters for each model state variable or a keys with phases and model state
        variable name for a distinction between different phases.
    caption : string
        Cation to put in the table.
    label : string
        Label for the table.
    parse : bool
        Parse parameters for latex.

    Returns
    -------
        String that can be printed for a latex table.
    """
    table = """\
\\begin{table}[ht]
    \\centering
    \\begin{tabular}{l|l|l}
        """
    tmp = "\\textbf{Model State Variable}"
    tmp2 = "\\textbf{Phase}"
    tmp3 = "\\textbf{{Top Parameters}}"
    n = len(tmp)
    n2 = len(tmp2)
    n3 = len(tmp3)
    for var in top:
        n = len(var.split("phase ")[1]) if (len(var.split("phase ")[1]) > n) else n
        phase = var.split("phase ")[0] + "phase"
        n3 = len(phase) if (len(phase) > n3) else n3
    table += f"{tmp:<{n}} & {tmp2:<{n2}} & {tmp3:<{n3}} \\\\ \\hline \n"
    for var in top:
        var1 = var.split("phase ")[1]
        table += (
            f"        {parse_word(var1):<{n}} & "
            if parse
            else f"        {var1:<{n}} & "
        )
        phase = var.split("phase ")[0] + "phase"
        table += f"{phase:<{n2}} & "
        line = "$"
        empty = " "
        n_breaks = 1
        for param in top[var]:
            if parse:
                line += f"{parse_word(param)}".replace(" $", "").replace("$", "")
            else:
                line += f"{param}"
            if len(line) > 100 * n_breaks and param != top[var][-1]:
                line += f"$ \\\\ \n{empty:<{n+8}} & {empty:<{n2}} & $ "
                n_breaks += 1.5
            else:
                line += ", "
        table += line + "$ \\\\ \n"
    table += "    \\end{tabular}\n"
    table += f"    \\caption{{{caption}}}\n"
    table += f"    \\label{{{label}}}\n"
    table += "\\end{table}"
    return table


def top_dict_to_table(top, caption, label, parse):
    """
    Parse a dictionary without phases to a altex table.

    Parameters
    ----------
    top : Dict
        A dictionary with keys the model state variable for a table with top
        parameters for each model state variable.
    parse : bool
        Parse parameters to latex names if True, otherwise use names used in the code.
    caption : string
        Caption for the table.
    label : string
        Label for the table.

    Returns
    -------
    String for a latex table or dictionary of model state variables with latex tables as values.
    """
    table = """\
\\begin{table}[ht]
\\centering
\\begin{tabular}{l|l}
    """
    tmp = "\\textbf{Model State Variable}"
    tmp2 = "\\textbf{{Top Parameters}}"
    n = len(tmp)
    n2 = len(tmp2)
    for var in top:
        if parse:
            if len(parse_word(var)) > n:
                n = len(parse_word(var))
        else:
            if len(var) > n:
                n = len(var)

    table += f"{tmp:<{n}} & {tmp2:<{n2}} "
    table += r"\\ \hline"
    table += "\n"
    for var in top:
        if parse:
            table += f"        {parse_word(var):<{n}} & "
        else:
            table += f"        {var:<{n}} & "
        line = "$"
        empty = " "
        n_breaks = 1
        for param in top[var]:
            if parse:
                line += f"{parse_word(param)}".replace(" $", "").replace("$", "")
            else:
                line += f"{param}"
            if len(line) > 100 * n_breaks and param != top[var][-1]:
                line += f"$ \\\\ \n{empty:<{n+8}} &  $ "
                n_breaks += 1.5
            else:
                line += ", "
        table += line + "$ \\\\ \n"
    table += "    \\end{tabular}\n"
    table += f"    \\caption{{{caption}}}\n"
    table += f"    \\label{{{label}}}\n"
    table += "\\end{table}"
    return table


# pylint: disable=too-complex, too-many-branches, too-many-statements
def print_latex_tables(ds, top=10, verbose=True):
    """
    Create a string to use in latex for the top n (='top') parameters for each model state variable.
    A model parameter is listed only for the model state variable where the sensitivity is the highest.

    Parameters
    ----------
    ds : xarray.Dataset
        Final, post-processed dataset with mean squared deviation and  predicted mean squared deviation.
    out_params : list-like of strings
        The model state variables for which sensitivities have been calculated for.
    top : int
        The number of top parameters to print a table for.
    verbose : Bool
        If True: print the parameters while building the table.

    Returns
    -------
    sort_key_list: a sorted list of (predicted squared error, model parameter, model state variable,
    string of a row for model parameters in latex) which is sorted by the name of the model state variable.
    table_dic: Dictionary with keys = model parameters where the value is a string of a row of the latex table.
    printed statements as string.
    """
    text = "\nBuild Latex tables\n"
    print(text)
    tmp_df = (
        ds.mean(dim=["trajectory", "time_after_ascent"], skipna=True)
        .to_dataframe()
        .reset_index()
    )

    table_dic = {}
    sort_key_list = []

    latexify_state = {
        "QV": r"\frac{\partial Q_\vapor}{",
        "QC": r"\frac{\partial Q_\cloud}{",
        "QR": r"\frac{\partial Q_\rain}{",
        "QG": r"\frac{\partial Q_\graupel}{",
        "QH": r"\frac{\partial Q_\hail}{",
        "QI": r"\frac{\partial Q_\ice}{",
        "QS": r"\frac{\partial Q_\snow}{",
        "NCCLOUD": r"\frac{\partial N_\cloud}{",
        "NCRAIN": r"\frac{\partial N_\rain}{",
        "NCGRAUPEL": r"\frac{\partial N_\graupel}{",
        "NCHAIL": r"\frac{\partial N_\hail}{",
        "NCICE": r"\frac{\partial N_\ice}{",
        "NCSNOW": r"\frac{\partial N_\snow}{",
        "QR_OUT": r"\frac{\partial Q_{\rain, \text{out}}}{",
        "QG_OUT": r"\frac{\partial Q_{\graupel, \text{out}}}{",
        "QH_OUT": r"\frac{\partial Q_{\hail, \text{out}}}{",
        "QI_OUT": r"\frac{\partial Q_{\ice, \text{out}}}{",
        "QS_OUT": r"\frac{\partial Q_{\snow, \text{out}}}{",
        "NR_OUT": r"\frac{\partial N_{\rain, \text{out}}}{",
        "NG_OUT": r"\frac{\partial N_{\graupel, \text{out}}}{",
        "NH_OUT": r"\frac{\partial N_{\hail, \text{out}}}{",
        "NI_OUT": r"\frac{\partial N_{\ice, \text{out}}}{",
        "NS_OUT": r"\frac{\partial N_{\snow, \text{out}}}{",
    }

    top_10_table = "\\begin{table}[hbt] \n \t\\centering \n \t\\begin{tabular}{ll}"
    top_10_table += "\n \t\t\\textbf{Model State Variable} \t& \\textbf{Top 10 Parameters} \\\\ \\hline \n"
    sedi_latex = ""
    sedi_started = False
    long_table_dic = {}

    for out_p, out_p_latexified in latexify_state.items():
        if "OUT" in out_p:
            if sedi_started:
                sedi_latex = sedi_latex[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
            else:
                sedi_latex = "\t\t Sedimentation \t& $ "
                sedi_started = True
        else:
            top_10_table += "\t\t" + parse_word(out_p) + "\t& $ "
        if verbose:
            print(f"########################### {out_p} ########################")
        df = tmp_df.loc[tmp_df["Output Parameter"] == out_p]
        # Ignore parameters that never appeared in unperturbed versions
        if np.max(df["Predicted Squared Error"]) == 0:
            continue
        if verbose:
            print("sort by sensitivity")
            print(
                df.nlargest(top, "Predicted Squared Error")[
                    ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
                ]
            )
        tmp = df.nlargest(top, "Predicted Squared Error")[
            ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
        ]
        i = 0
        for _, row in tmp.iterrows():
            if i == 5:
                if "OUT" in out_p:
                    sedi_latex = sedi_latex[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
                else:
                    top_10_table = top_10_table[:-2] + "$ \\\\ \n\t\t\t\t\t\t & $ "
            i += 1
            if "OUT" in out_p:
                sedi_latex += (
                    parse_word(row["Input Parameter"])
                    .replace("$", "")
                    .replace(r"\partial", "")
                    + ", "
                )
            else:
                top_10_table += (
                    parse_word(row["Input Parameter"])
                    .replace("$", "")
                    .replace(r"\partial", "")
                    + ", "
                )
            found = False
            val = ""
            param = ""
            state_var = ""
            l_string = ""
            for val, param, state_var, l_string in sort_key_list:
                if param == row["Input Parameter"] and (
                    val < row["Predicted Squared Error"]
                    or ("N" in state_var and "N" not in out_p)
                ):
                    if "N" not in state_var and "N" in out_p:
                        break

                    found = True
                    if verbose:
                        print(f"Replace ({val}, {param}, {state_var})")
                        print("With (")
                        print(row["Predicted Squared Error"], end=", ")
                        print(row["Input Parameter"], end=", ")
                        print(out_p, end=")\n")
                    break
            if found:
                sort_key_list.remove((val, param, state_var, l_string))

            if row["Input Parameter"] not in table_dic or found:

                group = None
                for key, grouped in in_params_grouping.items():
                    if row["Input Parameter"] in grouped:
                        group = key

                def latex_my_number(x, pred_err):
                    if x == 0:
                        return "$ 0.00 $"
                    if x >= 100:
                        exponent = int(np.log10(x))
                        var = x / 10**exponent
                        return f"$ {var:2.2f} \\times 10^{ {exponent} } $"
                    if x < 0.01:
                        exponent = math.floor(np.log10(x))
                        var = x * 10 ** (-exponent)
                        return f"$ {var:2.2f} \\times 10^{ {exponent} } $"
                    return f"$ {pred_err:2.2f} $"

                long_string = (
                    parse_word(row["Input Parameter"]).replace(r"\partial", "")
                    + " & "
                    + latex_my_number(
                        row["Mean Squared Error"], row["Predicted Squared Error"]
                    )
                    + " & "
                    + latex_my_number(
                        row["Predicted Squared Error"], row["Predicted Squared Error"]
                    )
                    + " & "
                    + r"\textbf{"
                    + group.title()
                    + "}: "
                    + in_params_descr_dic[row["Input Parameter"]]
                    + " \\\\ "
                )
                if out_p not in long_table_dic:
                    long_table_dic[out_p] = [long_string]
                else:
                    long_table_dic[out_p].append(long_string)
                sort_key_list.append(
                    (
                        row["Predicted Squared Error"],
                        row["Input Parameter"],
                        out_p,
                        long_string,
                    )
                )

                table_dic[row["Input Parameter"]] = (
                    r"$ \displaystyle "
                    + out_p_latexified
                    + parse_word(row["Input Parameter"]).replace("$", "")
                    + r"} $ & "
                    + latex_my_number(
                        row["Mean Squared Error"], row["Predicted Squared Error"]
                    )
                    + " & "
                    + latex_my_number(
                        row["Predicted Squared Error"], row["Predicted Squared Error"]
                    )
                    + " & "
                    + r"\textbf{"
                    + group.title()
                    + "}: "
                    + in_params_descr_dic[row["Input Parameter"]]
                    + " \\\\ "
                )
        if "OUT" not in out_p:
            top_10_table = top_10_table[:-2] + " $ \\\\ \n"
        if verbose:
            print("sort by Predicted Squared Error")
            print(
                df.nlargest(top, "Predicted Squared Error")[
                    ["Input Parameter", "Predicted Squared Error", "Mean Squared Error"]
                ]
            )

    top_10_table += sedi_latex[:-2] + " $ \\\\"
    top_10_table += "\n\t\\end{tabular} \n \t\\caption{} \n"
    top_10_table += "\t\\label{tab:} \n \\end{table} \n"
    text += "\nThe table of top 10 parameters for each state variable:\n"
    text += top_10_table
    print("\nThe table of top 10 parameters for each state variable:\n")
    print(top_10_table)

    if verbose:
        print(f"There are {len(table_dic)} different input parameters")

    text += "\nAppendix table of top parameters\n"
    print("\nAppendix table of top parameters\n")
    tmp_sort = sorted(sort_key_list, key=lambda x: (x[2], x[0]), reverse=True)
    sort_dic_long_table = {}
    sort_dic_short_table = {}
    for sens, key, state_variable, l_string in tmp_sort:
        if "_OUT" in state_variable:
            state_variable2 = "Sedimentation"
        else:
            state_variable2 = state_variable
        if state_variable not in sort_dic_long_table:
            sort_dic_long_table[state_variable] = [(sens, key, l_string)]
        else:
            sort_dic_long_table[state_variable].append((sens, key, l_string))

        if state_variable2 not in sort_dic_short_table:
            sort_dic_short_table[state_variable2] = [(sens, key, l_string)]
        else:
            sort_dic_short_table[state_variable2].append((sens, key, l_string))
    table_text = "\\bgroup\n"
    table_text += (
        "\\def\\arraystretch{1.2} %  1 is the default, "
        "we want it slightly larger such that exponents are easier to read\n"
    )
    table_text += "\\begin{tabularx}{\\linewidth}{@{}lccX@{}}\n"
    table_text += (
        "\t\\textbf{Model Param.}  & \\textbf{MSD} & "
        "\\textbf{Predicted MSD} & \\textbf{Parameter Description}\n"
    )
    table_text += "\t\\endhead\n"
    i = 0
    for state_variable in latexify_state:
        if state_variable not in sort_dic_long_table:
            continue
        table_text += (
            "\t\t\\hline \\multicolumn{4}{c}{"
            + parse_word(state_variable).title()
            + "}\t \\\\ \\hline\n"
        )
        for _, _, iter_string in sort_dic_long_table[state_variable]:
            table_text += "\t\t"
            table_text += f"{iter_string}"
            i += 1
    table_text += "\t"
    top_str = "ten"
    if top != 10:
        top_str = str(top)
    table_text += (
        r"\caption{This is the set of parameters if we gather the "
        + top_str
        + r" most important ones for each model state variable. "
        r"The predicted MSD is defined in Equation~\ref{eq:identification:msd_predict}, "
        r"where we only show the highest predicted MSD among all mass densities unless "
        r"the parameter did not have an impact on mass densities. "
        r"In that case, the predicted deviation on number density and precipitation is considered. There are $ "
        + str(i)
        + r" $ different parameters in total.}"
    )
    table_text += "\t\\label{tab:important_params}\n"
    table_text += "\\end{tabularx}\n"
    table_text += "\\egroup\n"
    text += table_text
    print(table_text)
    return sort_key_list, table_dic, text