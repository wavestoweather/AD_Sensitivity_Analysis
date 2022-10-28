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
        self.index_in_dic = {
            "da_1": 0,
            "da_2": 1,
            "de_1": 2,
            "de_2": 3,
            "dd": 4,
            "dN_c": 5,
            "dgamma": 6,
            "dbeta_c": 7,
            "dbeta_r": 8,
            "ddelta1": 9,
            "ddelta2": 10,
            "dzeta": 11,
            "drain_gfak": 12,
            "dcloud_k_au": 13,
            "dcloud_k_sc": 14,
            "dkc_autocon": 15,
            "dinv_z": 16,
            "dw": 17,
            "dq_crit_i": 18,
            "dD_crit_i": 19,
            "dD_conv_i": 20,
            "dq_crit_r": 21,
            "dD_crit_r": 22,
            "dq_crit_fr": 23,
            "dD_coll_c": 24,
            "dq_crit": 25,
            "dD_conv_sg": 26,
            "dD_conv_ig": 27,
            "dx_conv": 28,
            "dparcel_height": 29,
            "dalpha_spacefilling": 30,
            "dT_nuc": 31,
            "dT_freeze": 32,
            "dT_f": 33,
            "dD_eq": 34,
            "drho_w": 35,
            "drho_0": 36,
            "drho_vel": 37,
            "drho_vel_c": 38,
            "drho_ice": 39,
            "dM_w": 40,
            "dM_a": 41,
            "dR_universal": 42,
            "dEpsilon": 43,
            "dgravity_acc": 44,
            "dR_a": 45,
            "dR_v": 46,
            "da_v": 47,
            "db_v": 48,
            "da_prime": 49,
            "db_prime": 50,
            "dc_prime": 51,
            "dK_T": 52,
            "dL_wd": 53,
            "dL_ed": 54,
            "dD_v": 55,
            "decoll_min": 56,
            "decoll_gg": 57,
            "decoll_gg_wet": 58,
            "dkin_visc_air": 59,
            "dC_mult": 60,
            "dT_mult_min": 61,
            "dT_mult_max": 62,
            "dT_mult_opt": 63,
            "dconst0": 64,
            "dconst3": 65,
            "dconst4": 66,
            "dconst5": 67,
            "dD_rainfrz_gh": 68,
            "dD_rainfrz_ig": 69,
            "ddv0": 70,
            "dp_sat_melt": 71,
            "dcp": 72,
            "dk_b": 73,
            "da_HET": 74,
            "db_HET": 75,
            "dN_sc": 76,
            "dn_f": 77,
            "dN_avo": 78,
            "dna_dust": 79,
            "dna_soot": 80,
            "dna_orga": 81,
            "dni_het_max": 82,
            "dni_hom_max": 83,
            "da_dep": 84,
            "db_dep": 85,
            "dc_dep": 86,
            "dd_dep": 87,
            "dnim_imm": 88,
            "dnin_dep": 89,
            "dalf_imm": 90,
            "dbet_dep": 91,
            "dbet_imm": 92,
            "dr_const": 93,
            "dr1_const": 94,
            "dcv": 95,
            "dp_sat_const_a": 96,
            "dp_sat_ice_const_a": 97,
            "dp_sat_const_b": 98,
            "dp_sat_ice_const_b": 99,
            "dp_sat_low_temp": 100,
            "dT_sat_low_temp": 101,
            "dalpha_depo": 102,
            "dr_0": 103,
            "dk_1_conv": 104,
            "dk_2_conv": 105,
            "dk_1_accr": 106,
            "dk_r": 107,
            "da_ccn_1": 108,
            "da_ccn_2": 109,
            "da_ccn_3": 110,
            "da_ccn_4": 111,
            "db_ccn_1": 112,
            "db_ccn_2": 113,
            "db_ccn_3": 114,
            "db_ccn_4": 115,
            "dc_ccn_1": 116,
            "dc_ccn_2": 117,
            "dc_ccn_3": 118,
            "dc_ccn_4": 119,
            "dd_ccn_1": 120,
            "dd_ccn_2": 121,
            "dd_ccn_3": 122,
            "dd_ccn_4": 123,
            "dp_ccn": 124,
            "dh_ccn_1": 125,
            "dh_ccn_2": 126,
            "dg_ccn_1": 127,
            "dg_ccn_2": 128,
            "dg_ccn_3": 129,
            "di_ccn_1": 130,
            "di_ccn_2": 131,
            "dhande_ccn_fac": 132,
            "dD_br_threshold": 133 + offset,
            "dk_br": 134 + offset,
            "dD_br": 135 + offset,
            "dc_br": 136 + offset,
            "drain_a_geo": 137 + offset,
            "drain_b_geo": 138 + offset,
            "drain_min_x": 139 + offset,
            "drain_min_x_act": 140 + offset,
            "drain_min_x_nuc_homo": 141 + offset,
            "drain_min_x_nuc_hetero": 142 + offset,
            "drain_min_x_melt": 143 + offset,
            "drain_min_x_evap": 144 + offset,
            "drain_min_x_freezing": 145 + offset,
            "drain_min_x_depo": 146 + offset,
            "drain_min_x_collision": 147 + offset,
            "drain_min_x_collection": 148 + offset,
            "drain_min_x_conversion": 149 + offset,
            "drain_min_x_sedimentation": 150 + offset,
            "drain_min_x_riming": 151 + offset,
            "drain_max_x": 152 + offset,
            "drain_sc_theta_q": 153 + offset,
            "drain_sc_delta_q": 154 + offset,
            "drain_sc_theta_n": 155 + offset,
            "drain_sc_delta_n": 156 + offset,
            "drain_s_vel": 157 + offset,
            "drain_a_vel": 158 + offset,
            "drain_b_vel": 159 + offset,
            "drain_rho_v": 160 + offset,
            "drain_c_z": 161 + offset,
            "drain_sc_coll_n": 162 + offset,
            "drain_cmu0": 163 + offset,
            "drain_cmu1": 164 + offset,
            "drain_cmu2": 165 + offset,
            "drain_cmu3": 166 + offset,
            "drain_cmu4": 167 + offset,
            "drain_cmu5": 168 + offset,
            "drain_alpha": 169 + offset,
            "drain_beta": 170 + offset,
            "drain_gamma": 171 + offset,
            "drain_nu": 172 + offset,
            "drain_g1": 173 + offset,
            "drain_g2": 174 + offset,
            "drain_mu": 175 + offset,
            "drain_nm1": 176 + offset,
            "drain_nm2": 177 + offset,
            "drain_nm3": 178 + offset,
            "drain_q_crit_c": 179 + offset,
            "drain_d_crit_c": 180 + offset,
            "drain_ecoll_c": 181 + offset,
            "drain_cap": 182 + offset,
            "drain_a_ven": 183 + offset,
            "drain_b_ven": 184 + offset,
            "drain_c_s": 185 + offset,
            "drain_a_f": 186 + offset,
            "drain_b_f": 187 + offset,
            "drain_alfa_n": 188 + offset,
            "drain_alfa_q": 189 + offset,
            "drain_lambda": 190 + offset,
            "drain_vsedi_min": 191 + offset,
            "drain_vsedi_max": 192 + offset,
            "dcloud_a_geo": 193 + offset,
            "dcloud_b_geo": 194 + offset,
            "dcloud_min_x": 195 + offset,
            "dcloud_min_x_act": 196 + offset,
            "dcloud_min_x_nuc_homo": 197 + offset,
            "dcloud_min_x_nuc_hetero": 198 + offset,
            "dcloud_min_x_melt": 199 + offset,
            "dcloud_min_x_evap": 200 + offset,
            "dcloud_min_x_freezing": 201 + offset,
            "dcloud_min_x_depo": 202 + offset,
            "dcloud_min_x_collision": 203 + offset,
            "dcloud_min_x_collection": 204 + offset,
            "dcloud_min_x_conversion": 205 + offset,
            "dcloud_min_x_sedimentation": 206 + offset,
            "dcloud_min_x_riming": 207 + offset,
            "dcloud_max_x": 208 + offset,
            "dcloud_sc_theta_q": 209 + offset,
            "dcloud_sc_delta_q": 210 + offset,
            "dcloud_sc_theta_n": 211 + offset,
            "dcloud_sc_delta_n": 212 + offset,
            "dcloud_s_vel": 213 + offset,
            "dcloud_a_vel": 214 + offset,
            "dcloud_b_vel": 215 + offset,
            "dcloud_rho_v": 216 + offset,
            "dcloud_c_z": 217 + offset,
            "dcloud_sc_coll_n": 218 + offset,
            "dcloud_cmu0": 219 + offset,
            "dcloud_cmu1": 220 + offset,
            "dcloud_cmu2": 221 + offset,
            "dcloud_cmu3": 222 + offset,
            "dcloud_cmu4": 223 + offset,
            "dcloud_cmu5": 224 + offset,
            "dcloud_alpha": 225 + offset,
            "dcloud_beta": 226 + offset,
            "dcloud_gamma": 227 + offset,
            "dcloud_nu": 228 + offset,
            "dcloud_g1": 229 + offset,
            "dcloud_g2": 230 + offset,
            "dcloud_mu": 231 + offset,
            "dcloud_nm1": 232 + offset,
            "dcloud_nm2": 233 + offset,
            "dcloud_nm3": 234 + offset,
            "dcloud_q_crit_c": 235 + offset,
            "dcloud_d_crit_c": 236 + offset,
            "dcloud_ecoll_c": 237 + offset,
            "dcloud_cap": 238 + offset,
            "dcloud_a_ven": 239 + offset,
            "dcloud_b_ven": 240 + offset,
            "dcloud_c_s": 241 + offset,
            "dcloud_a_f": 242 + offset,
            "dcloud_b_f": 243 + offset,
            "dcloud_alfa_n": 244 + offset,
            "dcloud_alfa_q": 245 + offset,
            "dcloud_lambda": 246 + offset,
            "dcloud_vsedi_min": 247 + offset,
            "dcloud_vsedi_max": 248 + offset,
            "dgraupel_a_geo": 249 + offset,
            "dgraupel_b_geo": 250 + offset,
            "dgraupel_min_x": 251 + offset,
            "dgraupel_min_x_act": 252 + offset,
            "dgraupel_min_x_nuc_homo": 253 + offset,
            "dgraupel_min_x_nuc_hetero": 254 + offset,
            "dgraupel_min_x_melt": 255 + offset,
            "dgraupel_min_x_evap": 256 + offset,
            "dgraupel_min_x_freezing": 257 + offset,
            "dgraupel_min_x_depo": 258 + offset,
            "dgraupel_min_x_collision": 259 + offset,
            "dgraupel_min_x_collection": 260 + offset,
            "dgraupel_min_x_conversion": 261 + offset,
            "dgraupel_min_x_sedimentation": 262 + offset,
            "dgraupel_min_x_riming": 263 + offset,
            "dgraupel_max_x": 264 + offset,
            "dgraupel_sc_theta_q": 265 + offset,
            "dgraupel_sc_delta_q": 266 + offset,
            "dgraupel_sc_theta_n": 267 + offset,
            "dgraupel_sc_delta_n": 268 + offset,
            "dgraupel_s_vel": 269 + offset,
            "dgraupel_a_vel": 270 + offset,
            "dgraupel_b_vel": 271 + offset,
            "dgraupel_rho_v": 272 + offset,
            "dgraupel_c_z": 273 + offset,
            "dgraupel_sc_coll_n": 274 + offset,
            "dgraupel_cmu0": 275 + offset,
            "dgraupel_cmu1": 276 + offset,
            "dgraupel_cmu2": 277 + offset,
            "dgraupel_cmu3": 278 + offset,
            "dgraupel_cmu4": 279 + offset,
            "dgraupel_cmu5": 280 + offset,
            "dgraupel_alpha": 281 + offset,
            "dgraupel_beta": 282 + offset,
            "dgraupel_gamma": 283 + offset,
            "dgraupel_nu": 284 + offset,
            "dgraupel_g1": 285 + offset,
            "dgraupel_g2": 286 + offset,
            "dgraupel_mu": 287 + offset,
            "dgraupel_nm1": 288 + offset,
            "dgraupel_nm2": 289 + offset,
            "dgraupel_nm3": 290 + offset,
            "dgraupel_q_crit_c": 291 + offset,
            "dgraupel_d_crit_c": 292 + offset,
            "dgraupel_ecoll_c": 293 + offset,
            "dgraupel_cap": 294 + offset,
            "dgraupel_a_ven": 295 + offset,
            "dgraupel_b_ven": 296 + offset,
            "dgraupel_c_s": 297 + offset,
            "dgraupel_a_f": 298 + offset,
            "dgraupel_b_f": 299 + offset,
            "dgraupel_alfa_n": 300 + offset,
            "dgraupel_alfa_q": 301 + offset,
            "dgraupel_lambda": 302 + offset,
            "dgraupel_vsedi_min": 303 + offset,
            "dgraupel_vsedi_max": 304 + offset,
            "dhail_a_geo": 305 + offset,
            "dhail_b_geo": 306 + offset,
            "dhail_min_x": 307 + offset,
            "dhail_min_x_act": 308 + offset,
            "dhail_min_x_nuc_homo": 309 + offset,
            "dhail_min_x_nuc_hetero": 310 + offset,
            "dhail_min_x_melt": 311 + offset,
            "dhail_min_x_evap": 312 + offset,
            "dhail_min_x_freezing": 313 + offset,
            "dhail_min_x_depo": 314 + offset,
            "dhail_min_x_collision": 315 + offset,
            "dhail_min_x_collection": 316 + offset,
            "dhail_min_x_conversion": 317 + offset,
            "dhail_min_x_sedimentation": 318 + offset,
            "dhail_min_x_riming": 319 + offset,
            "dhail_max_x": 320 + offset,
            "dhail_sc_theta_q": 321 + offset,
            "dhail_sc_delta_q": 322 + offset,
            "dhail_sc_theta_n": 323 + offset,
            "dhail_sc_delta_n": 324 + offset,
            "dhail_s_vel": 325 + offset,
            "dhail_a_vel": 326 + offset,
            "dhail_b_vel": 327 + offset,
            "dhail_rho_v": 328 + offset,
            "dhail_c_z": 329 + offset,
            "dhail_sc_coll_n": 330 + offset,
            "dhail_cmu0": 331 + offset,
            "dhail_cmu1": 332 + offset,
            "dhail_cmu2": 333 + offset,
            "dhail_cmu3": 334 + offset,
            "dhail_cmu4": 335 + offset,
            "dhail_cmu5": 336 + offset,
            "dhail_alpha": 337 + offset,
            "dhail_beta": 338 + offset,
            "dhail_gamma": 339 + offset,
            "dhail_nu": 340 + offset,
            "dhail_g1": 341 + offset,
            "dhail_g2": 342 + offset,
            "dhail_mu": 343 + offset,
            "dhail_nm1": 344 + offset,
            "dhail_nm2": 345 + offset,
            "dhail_nm3": 346 + offset,
            "dhail_q_crit_c": 347 + offset,
            "dhail_d_crit_c": 348 + offset,
            "dhail_ecoll_c": 349 + offset,
            "dhail_cap": 350 + offset,
            "dhail_a_ven": 351 + offset,
            "dhail_b_ven": 352 + offset,
            "dhail_c_s": 353 + offset,
            "dhail_a_f": 354 + offset,
            "dhail_b_f": 355 + offset,
            "dhail_alfa_n": 356 + offset,
            "dhail_alfa_q": 357 + offset,
            "dhail_lambda": 358 + offset,
            "dhail_vsedi_min": 359 + offset,
            "dhail_vsedi_max": 360 + offset,
            "dice_a_geo": 361 + offset,
            "dice_b_geo": 362 + offset,
            "dice_min_x": 363 + offset,
            "dice_min_x_act": 364 + offset,
            "dice_min_x_nuc_homo": 365 + offset,
            "dice_min_x_nuc_hetero": 366 + offset,
            "dice_min_x_melt": 367 + offset,
            "dice_min_x_evap": 368 + offset,
            "dice_min_x_freezing": 369 + offset,
            "dice_min_x_depo": 370 + offset,
            "dice_min_x_collision": 371 + offset,
            "dice_min_x_collection": 372 + offset,
            "dice_min_x_conversion": 373 + offset,
            "dice_min_x_sedimentation": 374 + offset,
            "dice_min_x_riming": 375 + offset,
            "dice_max_x": 376 + offset,
            "dice_sc_theta_q": 377 + offset,
            "dice_sc_delta_q": 378 + offset,
            "dice_sc_theta_n": 379 + offset,
            "dice_sc_delta_n": 380 + offset,
            "dice_s_vel": 381 + offset,
            "dice_a_vel": 382 + offset,
            "dice_b_vel": 383 + offset,
            "dice_rho_v": 384 + offset,
            "dice_c_z": 385 + offset,
            "dice_sc_coll_n": 386 + offset,
            "dice_cmu0": 387 + offset,
            "dice_cmu1": 388 + offset,
            "dice_cmu2": 389 + offset,
            "dice_cmu3": 390 + offset,
            "dice_cmu4": 391 + offset,
            "dice_cmu5": 392 + offset,
            "dice_alpha": 393 + offset,
            "dice_beta": 394 + offset,
            "dice_gamma": 395 + offset,
            "dice_nu": 396 + offset,
            "dice_g1": 397 + offset,
            "dice_g2": 398 + offset,
            "dice_mu": 399 + offset,
            "dice_nm1": 400 + offset,
            "dice_nm2": 401 + offset,
            "dice_nm3": 402 + offset,
            "dice_q_crit_c": 403 + offset,
            "dice_d_crit_c": 404 + offset,
            "dice_ecoll_c": 405 + offset,
            "dice_cap": 406 + offset,
            "dice_a_ven": 407 + offset,
            "dice_b_ven": 408 + offset,
            "dice_c_s": 409 + offset,
            "dice_a_f": 410 + offset,
            "dice_b_f": 411 + offset,
            "dice_alfa_n": 412 + offset,
            "dice_alfa_q": 413 + offset,
            "dice_lambda": 414 + offset,
            "dice_vsedi_min": 415 + offset,
            "dice_vsedi_max": 416 + offset,
            "dsnow_a_geo": 417 + offset,
            "dsnow_b_geo": 418 + offset,
            "dsnow_min_x": 419 + offset,
            "dsnow_min_x_act": 420 + offset,
            "dsnow_min_x_nuc_homo": 421 + offset,
            "dsnow_min_x_nuc_hetero": 422 + offset,
            "dsnow_min_x_melt": 423 + offset,
            "dsnow_min_x_evap": 424 + offset,
            "dsnow_min_x_freezing": 425 + offset,
            "dsnow_min_x_depo": 426 + offset,
            "dsnow_min_x_collision": 427 + offset,
            "dsnow_min_x_collection": 428 + offset,
            "dsnow_min_x_conversion": 429 + offset,
            "dsnow_min_x_sedimentation": 430 + offset,
            "dsnow_min_x_riming": 431 + offset,
            "dsnow_max_x": 432 + offset,
            "dsnow_sc_theta_q": 433 + offset,
            "dsnow_sc_delta_q": 434 + offset,
            "dsnow_sc_theta_n": 435 + offset,
            "dsnow_sc_delta_n": 436 + offset,
            "dsnow_s_vel": 437 + offset,
            "dsnow_a_vel": 438 + offset,
            "dsnow_b_vel": 439 + offset,
            "dsnow_rho_v": 440 + offset,
            "dsnow_c_z": 441 + offset,
            "dsnow_sc_coll_n": 442 + offset,
            "dsnow_cmu0": 443 + offset,
            "dsnow_cmu1": 444 + offset,
            "dsnow_cmu2": 445 + offset,
            "dsnow_cmu3": 446 + offset,
            "dsnow_cmu4": 447 + offset,
            "dsnow_cmu5": 448 + offset,
            "dsnow_alpha": 449 + offset,
            "dsnow_beta": 450 + offset,
            "dsnow_gamma": 451 + offset,
            "dsnow_nu": 452 + offset,
            "dsnow_g1": 453 + offset,
            "dsnow_g2": 454 + offset,
            "dsnow_mu": 455 + offset,
            "dsnow_nm1": 456 + offset,
            "dsnow_nm2": 457 + offset,
            "dsnow_nm3": 458 + offset,
            "dsnow_q_crit_c": 459 + offset,
            "dsnow_d_crit_c": 460 + offset,
            "dsnow_ecoll_c": 461 + offset,
            "dsnow_cap": 462 + offset,
            "dsnow_a_ven": 463 + offset,
            "dsnow_b_ven": 464 + offset,
            "dsnow_c_s": 465 + offset,
            "dsnow_a_f": 466 + offset,
            "dsnow_b_f": 467 + offset,
            "dsnow_alfa_n": 468 + offset,
            "dsnow_alfa_q": 469 + offset,
            "dsnow_lambda": 470 + offset,
            "dsnow_vsedi_min": 471 + offset,
            "dsnow_vsedi_max": 472 + offset,
        }

    def _get_gradient(self, out_param, in_param, gradients):
        return gradients[
            self.index_dic[out_param] * self.get_num_par() + self.index_in_dic[in_param]
        ]

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

    def riming_ice(self, qc, Nc, qi, Ni, qr, Nr, T, res, gradients):
        self.clib.physics_t_py_riming_ice(
            self.obj,
            qc,
            Nc,
            qi,
            Ni,
            qr,
            Nr,
            T,
            res.ctypes.data_as(self.c_double_arr),
            gradients.ctypes.data_as(self.c_double_arr),
        )

    def riming_snow(self, qc, Nc, qs, Ns, qr, Nr, T, res, gradients):
        self.clib.physics_t_py_riming_snow(
            self.obj,
            qc,
            Nc,
            qs,
            Ns,
            qr,
            Nr,
            T,
            res.ctypes.data_as(self.c_double_arr),
            gradients.ctypes.data_as(self.c_double_arr),
        )

    def riming_with_depo(
        self, qv, qc, Nc, qi, Ni, qs, Ns, qr, Nr, T, p, res, gradients
    ):
        self.clib.physics_t_py_riming_with_depo(
            self.obj,
            qv,
            qc,
            Nc,
            qi,
            Ni,
            qs,
            Ns,
            qr,
            Nr,
            T,
            p,
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

    def riming_ice_plot(
        self, qc, Nc, qi, Ni, qr, Nr, T, x, y, hue, logy, nx=10, nhue=10, out_param=""
    ):
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        if "final " in y:
            y_name = y[6::]
        else:
            y_name = y
        data = {x: [], y: [], hue: []}

        def set_vals(name, n):
            if name == "qc":
                return np.arange(qc[0], qc[1], (qc[1] - qc[0]) / (n + 1))
            elif name == "Nc":
                return np.arange(Nc[0], Nc[1], (Nc[1] - Nc[0]) / (n + 1))
            elif name == "qi":
                return np.arange(qi[0], qi[1], (qi[1] - qi[0]) / (n + 1))
            elif name == "Ni":
                return np.arange(Ni[0], Ni[1], (Ni[1] - Ni[0]) / (n + 1))
            elif name == "qr":
                return np.arange(qr[0], qr[1], (qr[1] - qr[0]) / (n + 1))
            elif name == "Nr":
                return np.arange(Nr[0], Nr[1], (Nr[1] - Nr[0]) / (n + 1))
            elif name == "T":
                return np.arange(T[0], T[1], (T[1] - T[0]) / (n + 1))

        x_vals = set_vals(x, nx)
        hue_vals = set_vals(hue, nhue)

        for x_v in x_vals:
            if x == "qc":
                qc = x_v
            elif x == "Nc":
                Nc = x_v
            elif x == "qi":
                qi = x_v
            elif x == "Ni":
                Ni = x_v
            elif x == "qr":
                qr = x_v
            elif x == "Nr":
                Nr = x_v
            elif x == "T":
                T = x_v

            for h_v in hue_vals:
                if hue == "qc":
                    qc = h_v
                elif hue == "Nc":
                    Nc = h_v
                elif hue == "qi":
                    qi = h_v
                elif hue == "Ni":
                    Ni = h_v
                elif hue == "qr":
                    qr = h_v
                elif hue == "Nr":
                    Nr = h_v
                elif hue == "T":
                    T = h_v
                self.riming_ice(qc, Nc, qi, Ni, qr, Nr, T, res, gradients)
                data[x].append(x_v)
                if "final " in y:
                    data[y].append(res[self.index_dic[y_name]])
                else:
                    data[y].append(self._get_gradient(out_param, y, gradients))
                data[hue].append(h_v)
        df = pd.DataFrame.from_dict(data)
        g = sns.scatterplot(data=df, x=x, y=y, hue=hue)
        if logy:
            g.axes.set_yscale("log")
        return g, df

    def riming_snow_plot(
        self, qc, Nc, qs, Ns, qr, Nr, T, x, y, hue, logy, nx=10, nhue=10, out_param=""
    ):
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        if "final " in y:
            y_name = y[6::]
        else:
            y_name = y
        data = {x: [], y: [], hue: []}

        def set_vals(name, n):
            if name == "qc":
                return np.arange(qc[0], qc[1], (qc[1] - qc[0]) / (n + 1))
            elif name == "Nc":
                return np.arange(Nc[0], Nc[1], (Nc[1] - Nc[0]) / (n + 1))
            elif name == "qs":
                return np.arange(qs[0], qs[1], (qs[1] - qs[0]) / (n + 1))
            elif name == "Ns":
                return np.arange(Ns[0], Ns[1], (Ns[1] - Ns[0]) / (n + 1))
            elif name == "qr":
                return np.arange(qr[0], qr[1], (qr[1] - qr[0]) / (n + 1))
            elif name == "Nr":
                return np.arange(Nr[0], Nr[1], (Nr[1] - Nr[0]) / (n + 1))
            elif name == "T":
                return np.arange(T[0], T[1], (T[1] - T[0]) / (n + 1))

        x_vals = set_vals(x, nx)
        hue_vals = set_vals(hue, nhue)

        for x_v in x_vals:
            if x == "qc":
                qc = x_v
            elif x == "Nc":
                Nc = x_v
            elif x == "qs":
                qs = x_v
            elif x == "Ns":
                Ns = x_v
            elif x == "qr":
                qr = x_v
            elif x == "Nr":
                Nr = x_v
            elif x == "T":
                T = x_v

            for h_v in hue_vals:
                if hue == "qc":
                    qc = h_v
                elif hue == "Nc":
                    Nc = h_v
                elif hue == "qs":
                    qs = h_v
                elif hue == "Ns":
                    Ns = h_v
                elif hue == "qr":
                    qr = h_v
                elif hue == "Nr":
                    Nr = h_v
                elif hue == "T":
                    T = h_v
                self.riming_snow(qc, Nc, qs, Ns, qr, Nr, T, res, gradients)
                data[x].append(x_v)
                if "final " in y:
                    data[y].append(res[self.index_dic[y_name]])
                else:
                    data[y].append(self._get_gradient(out_param, y, gradients))
                data[hue].append(h_v)
        df = pd.DataFrame.from_dict(data)
        g = sns.scatterplot(data=df, x=x, y=y, hue=hue)
        if logy:
            g.axes.set_yscale("log")
        return g, df

    def riming_with_depo_plot(
        self,
        qv,
        qc,
        Nc,
        qi,
        Ni,
        qs,
        Ns,
        qr,
        Nr,
        T,
        p,
        x,
        y,
        hue,
        logy,
        nx=10,
        nhue=10,
        out_param="",
    ):
        res = np.zeros(self.get_num_comp(), dtype=np.float64)
        gradients = np.zeros(self.get_num_comp() * self.get_num_par(), dtype=np.float64)
        if "final " in y:
            y_name = y[6::]
        else:
            y_name = y
        data = {x: [], y: [], hue: []}

        def set_vals(name, n):
            if name == "qv":
                return np.arange(qv[0], qv[1], (qv[1] - qv[0]) / (n + 1))
            elif name == "qc":
                return np.arange(qc[0], qc[1], (qc[1] - qc[0]) / (n + 1))
            elif name == "Nc":
                return np.arange(Nc[0], Nc[1], (Nc[1] - Nc[0]) / (n + 1))
            elif name == "qi":
                return np.arange(qi[0], qi[1], (qi[1] - qi[0]) / (n + 1))
            elif name == "Ni":
                return np.arange(Ni[0], Ni[1], (Ni[1] - Ni[0]) / (n + 1))
            elif name == "qs":
                return np.arange(qs[0], qs[1], (qs[1] - qs[0]) / (n + 1))
            elif name == "Ns":
                return np.arange(Ns[0], Ns[1], (Ns[1] - Ns[0]) / (n + 1))
            elif name == "qr":
                return np.arange(qr[0], qr[1], (qr[1] - qr[0]) / (n + 1))
            elif name == "Nr":
                return np.arange(Nr[0], Nr[1], (Nr[1] - Nr[0]) / (n + 1))
            elif name == "T":
                return np.arange(T[0], T[1], (T[1] - T[0]) / (n + 1))
            elif name == "p":
                return np.arange(p[0], p[1], (p[1] - p[0]) / (n + 1))

        x_vals = set_vals(x, nx)
        hue_vals = set_vals(hue, nhue)

        for x_v in x_vals:
            if x == "qv":
                qv = x_v
            elif x == "qc":
                qc = x_v
            elif x == "Nc":
                Nc = x_v
            elif x == "qi":
                qi = x_v
            elif x == "Ni":
                Ni = x_v
            elif x == "qs":
                qs = x_v
            elif x == "Ns":
                Ns = x_v
            elif x == "qr":
                qr = x_v
            elif x == "Nr":
                Nr = x_v
            elif x == "T":
                T = x_v
            elif x == "p":
                p = x_v

            for h_v in hue_vals:
                if hue == "qv":
                    qv = h_v
                elif hue == "qc":
                    qc = h_v
                elif hue == "Nc":
                    Nc = h_v
                elif hue == "qi":
                    qi = h_v
                elif hue == "Ni":
                    Ni = h_v
                elif hue == "qs":
                    qs = h_v
                elif hue == "Ns":
                    Ns = h_v
                elif hue == "qr":
                    qr = h_v
                elif hue == "Nr":
                    Nr = h_v
                elif hue == "T":
                    T = h_v
                elif hue == "p":
                    p = h_v
                self.riming_with_depo(
                    qv, qc, Nc, qi, Ni, qs, Ns, qr, Nr, T, p, res, gradients
                )
                data[x].append(x_v)
                if "final " in y:
                    data[y].append(res[self.index_dic[y_name]])
                else:
                    data[y].append(self._get_gradient(out_param, y, gradients))
                data[hue].append(h_v)
        df = pd.DataFrame.from_dict(data)
        g = sns.scatterplot(data=df, x=x, y=y, hue=hue)
        if logy:
            g.axes.set_yscale("log")
        return g, df


if __name__ == "__main__":
    import argparse
    import os
    import textwrap

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Load the C-library which provides an interface to microphysics processes. Plot some examples of these 
            processes.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--library_path",
        type=str,
        required=True,
        help=textwrap.dedent(
            """\
            Path to the C-library. You may find it at cmake-build-python_interface/lib/libpython_interface.so.
            """
        ),
    )
    parser.add_argument(
        "--timestep",
        type=float,
        default=30.0,
        help=textwrap.dedent(
            """\
            The timestep size in seconds.
            """
        ),
    )
    parser.add_argument(
        "--uncertainty",
        type=float,
        default=0.01,
        help=textwrap.dedent(
            """\
            The uncertainty that is used for all model parameters.
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        default="../pics/microphysics",
        type=str,
        help=textwrap.dedent(
            """\
            Path to store the generated images.
            """
        ),
    )
    args = parser.parse_args()
    phys = physics_t(lib_path=args.library_path)
    phys.setup_model_constants(30.0, 0.1)
    renderer = hv.Store.renderers["matplotlib"].instance(fig="png", dpi=300)

    def save_plot(plot, store_path):
        i = 0
        save = store_path + "_" + "{:03d}".format(i)
        while os.path.isfile(save + ".png"):
            i = i + 1
            save = store_path + "_" + "{:03d}".format(i)
        renderer.save(plot, save)

    plot, _ = phys.graupel_melting_plot(
        T=280,
        qg=2,
        Ng=8,
        x="qg",
        y="final qr",
        hue="T",
    )
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

    plot, _ = phys.riming_ice_plot(
        qc=[0, 1e-3],
        Nc=1e4,
        qi=[0, 1e-3],
        Ni=0,
        qr=1e-6,
        Nr=1e3,
        T=[260, 290],
        x="qc",
        y="qi",
        hue="T",
        logy=False,
        out_param="QV",
    )
    save_plot(plot, args.store_path)
