import matplotlib as mpl


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
            "text.usetex": True,
            "font.family": "serif",
            "axes.labelsize": 20,
            "font.size": 20,
            "legend.fontsize": 16,
            "xtick.labelsize": 12,
            "ytick.labelsize": 16
        })
    else:
        mpl.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "axes.labelsize": 10,
            "font.size": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8
        })


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
                 "rho", "nu", "mu", "lambda"]
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
                "qvout": "Sedimentation of Water Vapor Mixing Ratio",
                "qcout": "Sedimentation of Cloud Droplet Mixing Ratio",
                "qrout": "Sedimentation of Rain Droplet Mixing Ratio",
                "qsout": "Sedimentation of Snow Mixing Ratio",
                "qiout": "Sedimentation of Ice Mixing Ratio",
                "qgout": "Sedimentation of Graupel Mixing Ratio",
                "qhout": "Sedimentation of Hail Mixing Ratio",
                "LATITUDE": "Latitude",
                "LONGITDUE": "longitude",
                "z": "Height z",
                "w": "Ascend w",
                "MAP": "Flag for WCB-criterion"}
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
            try:
                word = (start + " "
                            + parts[2] + r"_{"
                            + parts[1] + r", \mathrm{"
                            + parts[0] + r"}}")
            except:
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
    return word