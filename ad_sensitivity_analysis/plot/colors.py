"""Definitions of various colors.

"""
from bokeh.palettes import Category20c
from matplotlib import colors as mpl_colors
from matplotlib import cm as mpl_cm
import matplotlib.pyplot as plt
import numpy as np

cmap_particles = {
    "QV": mpl_colors.to_hex([153 / 255, 204 / 255, 255 / 255]),
    "QC": mpl_colors.to_hex([102 / 255, 178 / 255, 255 / 255]),
    "QR": mpl_colors.to_hex([51 / 255, 153 / 255, 255 / 255]),
    "QI": mpl_colors.to_hex([215 / 255, 255 / 255, 255 / 255]),
    "QS": mpl_colors.to_hex([170 / 255, 255 / 255, 255 / 255]),
    "QG": mpl_colors.to_hex([0 / 255, 255 / 255, 255 / 255]),
    "QH": mpl_colors.to_hex([0 / 255, 204 / 255, 204 / 255]),
    "NC": mpl_colors.to_hex([0 / 255, 102 / 255, 255 / 255]),
    "NR": mpl_colors.to_hex([0 / 255, 0 / 255, 255 / 255]),
    "NI": mpl_colors.to_hex([102 / 255, 204 / 255, 255 / 255]),
    "NS": mpl_colors.to_hex([51 / 255, 204 / 255, 255 / 255]),
    "NG": mpl_colors.to_hex([0 / 255, 204 / 255, 255 / 255]),
    "NH": mpl_colors.to_hex([0 / 255, 153 / 255, 204 / 255]),
    "S": mpl_colors.to_hex([0 / 255, 204 / 255, 0 / 255]),
    "latent_heat": mpl_colors.to_hex([255 / 255, 51 / 255, 0 / 255]),
    "latent_cool": mpl_colors.to_hex([255 / 255, 153 / 255, 0 / 255]),
}

param_process_order = {
    "evaporation": [
        "dkin_visc_air",
        "drain_alpha",
        "drain_beta",
        "db_v",
        "drain_cmu3",
        "drain_nu",
    ],
    "geometry": [
        "drain_a_geo",
        "drain_b_geo",
        "dsnow_b_geo",
        "dice_a_geo",
        "dice_b_geo",
        "dgraupel_a_geo",
        "dgraupel_b_geo",
    ],
    "ccn": [
        "db_ccn_1",
        "db_ccn_2",
        "db_ccn_3",
        "dc_ccn_1",
        "dc_ccn_2",
        "dc_ccn_3",
        "dg_ccn_1",
        "dh_ccn_1",
        "di_ccn_1",
        "dhande_ccn_fac",
    ],
    "fall velocity": [
        "drain_a_vel",
        "drain_b_vel",
        "dice_b_vel",
        "dgraupel_b_vel",
    ],
    "misc": [
        "dp_sat_melt",
        "drain_mu",
    ],
}


# The order in the dataset
regrid_idx = [
    "db_ccn_1",
    "db_ccn_2",
    "db_ccn_3",
    "db_v",
    "dc_ccn_1",
    "dc_ccn_2",
    "dc_ccn_3",
    "dg_ccn_1",
    "dgraupel_a_geo",
    "dgraupel_b_geo",
    "dgraupel_b_vel",
    "dh_ccn_1",
    "dhande_ccn_fac",
    "di_ccn_1",
    "dice_a_geo",
    "dice_b_geo",
    "dice_b_vel",
    "dkin_visc_air",
    "dp_sat_melt",
    "drain_a_geo",
    "drain_a_vel",
    "drain_alpha",
    "drain_b_geo",
    "drain_b_vel",
    "drain_beta",
    "drain_cmu3",
    "drain_mu",
    "drain_nu",
    "dsnow_b_geo",
]


def get_cmap_types(backend="matplotlib"):
    """

    Parameters
    ----------
    backend

    Returns
    -------

    """
    if backend == "matplotlib":
        colors = plt.get_cmap("tab20c")
        return {
            "artificial": mpl_colors.to_hex(colors(1)[0:-1]),
            "artificial (threshold)": mpl_colors.to_hex(colors(5)[0:-1]),
            "physical": mpl_colors.to_hex(colors(9)[0:-1]),
            "physical (high variability)": mpl_colors.to_hex(colors(13)[0:-1]),
            "1-moment": mpl_colors.to_hex(colors(17)[0:-1]),
        }
    colors = Category20c[20]
    return {
        "artificial": colors[1],
        "artificial (threshold)": colors[5],
        "physical": colors[9],
        "physical (high variability)": colors[13],
        "1-moment": colors[17],
    }


def get_b8_cbar_colors(colorblind=True):
    """
    Colors to get the correct colorbar.

    Parameters
    ----------
    colorblind

    Returns
    -------

    """
    _, color_shades, _ = get_b8_colors(colorblind=colorblind)
    cbar_colors = []
    for _, params in param_process_order.items():
        for param in params:
            cbar_colors.append(color_shades[param])
    cbar_colors.append([0, 0, 0, 0])
    return cbar_colors


def get_b8_colors(colorblind=True):
    """
    Colors that match parameters.

    Parameters
    ----------
    colorblind

    Returns
    -------

    """

    color_shades = {}
    for key, params in param_process_order.items():
        n = len(params)
        if key == "evaporation":
            color = plt.cm.Blues(np.linspace(0, 1, n))
        elif key == "geometry":
            color = plt.cm.Greens(np.linspace(0.2, 1, n))
        elif key == "misc":
            color = plt.cm.Purples(np.linspace(0.2, 1, n))
        elif key == "ccn":
            color = plt.cm.Reds(np.linspace(0.2, 1, n))
        elif key == "fall velocity":
            if not colorblind:
                color = plt.cm.Oranges(np.linspace(0.2, 1, n))
            else:
                color = plt.cm.Greys(np.linspace(0.2, 1, n))

        # color_shades[key] = []
        for i, param in enumerate(params):
            color_shades[param] = color[i]

    param_color_order = []
    param_order = []
    for p in regrid_idx:
        param_order.append(p)
        param_color_order.append(color_shades[p])
    param_color_order.append([0, 0, 0, 0])
    return param_color_order, color_shades, param_order


def get_b8_process_from_param(param):
    """

    Parameters
    ----------
    param

    Returns
    -------

    """
    for key, params in param_process_order.items():
        if param in params:
            return key
    return None


def get_b8_cbar_ticks_labels(cbar_colors):
    """

    Parameters
    ----------
    cbar_colors

    Returns
    -------

    """
    cbar_labels = list(param_process_order.keys())
    cbar_ticks = np.asarray(
        [len(params) / 2 for key, params in param_process_order.items()]
    )
    n = 0
    for i, _ in enumerate(cbar_ticks):
        if i > 0:
            cbar_ticks[i] += n
            n += len(param_process_order[cbar_labels[i]])
        else:
            n = len(param_process_order[cbar_labels[i]])
    cbar_ticks /= len(cbar_colors) - 1
    return cbar_ticks, cbar_labels


def set_top_param_cbar(
    fig, ax, font_scale=1.0, cbarlabel="Index of Top Parameter", colorblind=True
):
    """

    Parameters
    ----------
    fig
    ax
    font_scale
    cbarlabel
    colorblind

    Returns
    -------

    """
    cbar_colors = get_b8_cbar_colors(colorblind=colorblind)
    cbar_ticks, cbar_labels = get_b8_cbar_ticks_labels(cbar_colors=cbar_colors)
    cmap = mpl_colors.ListedColormap(cbar_colors, N=len(cbar_colors) - 1)

    cbar = fig.colorbar(
        mpl_cm.ScalarMappable(cmap=cmap, norm=None),
        ax=ax,
        orientation="vertical",
        label=cbarlabel,
        ticks=cbar_ticks,
    )
    cbar.set_ticklabels(
        cbar_labels,
        fontsize=int(11 * font_scale),
    )
    cbar.set_label(
        label=cbarlabel,
        fontsize=int(11 * font_scale),
    )
