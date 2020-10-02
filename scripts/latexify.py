import matplotlib as mpl

mappings = {"lat_heat": "Latent Heating",
            "lat_cool": "Latent Cooling",
            "latent_heat": "Latent Heating",
            "latent_cool": "Latent Cooling",
            "dinv_z": r"$\partial z^{-1}$",
            "ratio_deriv": "Derivative Ratio",
            "in_param": "Input Parameter",
            "p": "Pressure",
            "T": "Temperature",
            "S": "Saturation",
            "qv": "Water Vapor Mixing Ratio",
            "qc": "Cloud Droplet Mixing Ratio",
            "qr": "Rain Droplet Mixing Ratio",
            "qs": "Snow Mixing Ratio",
            "qi": "Ice Mixing Ratio",
            "qg": "Graupel Mixing Ratio",
            "qh": "Hail Mixing Ratio",
            "Nv": "Water Vapor Particle Number",
            "Nc": "Cloud Droplet Particle Number",
            "Nr": "Rain Droplet Particle Number",
            "Ns": "Snow Particle Number",
            "Ni": "Ice Particle Number",
            "Ng": "Graupel Particle Number",
            "Nh": "Hail Particle Number",
            "qvout": "Precipitation of Water Vapor Mixing Ratio",
            "qcout": "Precipitation of Cloud Droplet Mixing Ratio",
            "qrout": "Precipitation of Rain Droplet Mixing Ratio",
            "qsout": "Precipitation of Snow Mixing Ratio",
            "qiout": "Precipitation of Ice Mixing Ratio",
            "qgout": "Precipitation of Graupel Mixing Ratio",
            "qhout": "Precipitation of Hail Mixing Ratio",
            "Nrout": "Precipitation of Rain Droplets",
            "Nsout": "Precipitation of Snow Crystals",
            "Niout": "Precipitation of Ice Crystals",
            "Ngout": "Precipitation of Graupel Particles",
            "Nhout": "Precipitation of Hail Particles",
            "LATITUDE": "Latitude",
            "LONGITDUE": "longitude",

            "pressure": "Pressure",
            "QV": "Water Vapor Mixing Ratio",
            "QC": "Cloud Droplet Mixing Ratio",
            "QR": "Rain Droplet Mixing Ratio",
            "QS": "Snow Mixing Ratio",
            "QI": "Ice Mixing Ratio",
            "QG": "Graupel Mixing Ratio",
            "QH": "Hail Mixing Ratio",
            "NCCLOUD": "Cloud Droplet Particle Number",
            "NCRAIN": "Rain Droplet Particle Number",
            "NCSNOW": "Snow Particle Number",
            "NCICE": "Ice Particle Number",
            "NCGRAUPEL": "Graupel Particle Number",
            "NCHAIL": "Hail Particle Number",
            "QR_IN": "sedimentation (from above) of rain droplet mixing ratio",
            "QS_IN": "sedimentation (from above) of snow crystal mixing ratio",
            "QI_IN": "sedimentation (from above) of ice crystal mixing ratio",
            "QG_IN": "sedimentation (from above) of graupel mixing ratio",
            "QR_OUT": "sedimentation of rain droplet mixing ratio",
            "QS_OUT": "sedimentation of snow crystal mixing ratio",
            "QI_OUT": "sedimentation of ice crystal mixing ratio",
            "QG_OUT": "sedimentation of graupel mixing ratio",
            "QH_OUT": "sedimentation of hail mixing ratio",

            "lat": "Latitude",
            "lon": "longitude",

            "z": "Height [m]",
            "w": "Ascend [m/s]",
            # "w": r"Ascend $[\frac{\mathrm{m}}{\mathrm{s}}]$",
            "MAP": "Flag for WCB-criterion",
            "Derivatives": "Derivatives",
            "timestep": "Time [s] after ascend begins",
            "time": "Time [s] after COSMO simulation begins",
            "time_after_ascent": "Time [s] after ascend begins",
            "step": "Simulation step",
            # Misc
#             "da_1", "da_2", "de_1", "de_2", "dd", "dN_c", "dgamma",
#             "dbeta_c", "dbeta_r", "ddelta1", "ddelta2", "dzeta",
#             "drain_gfak", "dcloud_k_au", "dcloud_k_sc", "dkc_autocon",
#             "dinv_z",
#             # Rain
#             "drain_a_geo", "drain_b_geo", "drain_min_x", "drain_max_x",
#             "drain_sc_theta_q", "drain_sc_delta_q", "drain_sc_theta_n",
#             "drain_sc_delta_n", "drain_s_vel", "drain_a_vel", "drain_b_vel",
#             "drain_rho_v", "drain_c_z", "drain_sc_coll_n", "drain_cmu0",
#             "drain_cmu1", "drain_cmu2", "drain_cmu3", "drain_cmu4",
#             "drain_cmu5", "drain_alpha", "drain_beta", "drain_gamma",
#             "drain_nu", "drain_g1", "drain_g2", "drain_mu", "drain_nm1",
#             "drain_nm2", "drain_nm3", "drain_q_crit_c", "drain_d_crit_c",
#             "drain_ecoll_c", "drain_cap", "drain_a_ven", "drain_b_ven",
#             "drain_c_s", "drain_a_f", "drain_b_f", "drain_alfa_n",
#             "drain_alfa_q", "drain_lambda", "drain_vsedi_min",
#             "drain_vsedi_max"
            # Cloud

            # Graupel

            # Hail

            # Ice

            # Snow


            "dmin_x_nuc_hetero": r"$\partial x_{\mathrm{min},\mathrm{nuc},\mathrm{hetero}}",
            "dmin_x_nuc_homo": r"$\partial x_{\mathrm{min},\mathrm{nuc},\mathrm{homo}}",
            "dmin_x_melt": r"$\partial x_{\mathrm{min},\mathrm{melt}}",
            "dmin_x_evap": r"$\partial x_{\mathrm{min},\mathrm{evap}}",
            "dmin_x_freezing": r"$\partial x_{\mathrm{min},\mathrm{freezing}}",
            "dmin_x_depo": r"$\partial x_{\mathrm{min},\mathrm{depo}}",
            "dmin_x_collision": r"$\partial x_{\mathrm{min},\mathrm{collision}}",
            "dmin_x_collection": r"$\partial x_{\mathrm{min},\mathrm{collection}}",
            "dmin_x_conversion": r"$\partial x_{\mathrm{min},\mathrm{conversion}}",
            "dmin_x_sedimentation": r"$\partial x_{\mathrm{min},\mathrm{sedimentation}}",
            "dmin_x_riming": r"$\partial x_{\mathrm{min},\mathrm{riming}}"}

def set_size(beamer=True):
    """
    Set some options to use latex.

    Parameters
    ----------
    beamer : bool
        Beamer is used for bigger texts.
    """
    if beamer:
        mpl.rcParams.update({
            "text.usetex": False,
            "font.family": "serif",
            "axes.labelsize": 20,
            "font.size": 20,
            "legend.fontsize": 16,
            "xtick.labelsize": 12,
            "ytick.labelsize": 16
        })
    else:
        mpl.rcParams.update({
            "text.usetex": False,
            "font.family": "serif",
            "axes.labelsize": 10,
            "font.size": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8
        })


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
        "timestep": "s"}
    if param in unit_dic:
        if brackets:
            return "[" + unit_dic[param] + "]"
        return unit_dic[param]
    else:
        return ""

def parse_word(word):
    """
    Parse a name of a derivative and return it in a latex conform type.

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
    subscript_no_math = ["snow", "graupel", "rain", "ice", "cloud", "hail", "vapor"]
    no_math = ["geo", "min", "max", "ven", "vel"]
    math_keys = ["alpha", "gamma", "beta", "delta", "zeta",
                 "rho", "nu", "mu", "lambda", "theta"]

    maps = mappings.keys()
    for w in maps:
        if word == w:
            return mappings[word]
    # The first "d" shall be a derivative symbol
    word = r"\partial " + word[1::]

    word = word.replace("delta_", "delta ")
    # Check for this typo
    word = word.replace("alfa", "alpha")
    # if any of words is in there, make it to subscript
    for w in subscript_no_math:
        if w in word:
            parts = word.split(" ")
            start = parts[0]
            parts = parts[1].split("_")
            if len(parts) == 4:
                word = (start + " "
                        + parts[2] + r"_{"
                        + parts[1] + r", \mathrm{"
                        + parts[0] + r","
                        + parts[3] + r"}}")

            elif len(parts) == 3:
                word = (start + " "
                        + parts[2] + r"_{"
                        + parts[1] + r", \mathrm{"
                        + parts[0] + r"}}")
            else:
                word = (start + " "
                        + parts[1] + r"_{"
                        + r"\mathrm{"
                        + parts[0] + r"}}")

            break
    word = r"$" + word + r"$"
    for w in no_math:
        word = word.replace(w, r"\mathrm{" + w + r"}")
    for w in math_keys:
        word = word.replace(w, "\\" + w )
    if "nuc" in word:
        word = word.replace("\\nuc", "nuc")
    return word