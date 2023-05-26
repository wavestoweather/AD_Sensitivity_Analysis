"""Functions to save plots and similar things.

"""
import os
import sys

import matplotlib.pyplot as plt


def get_save_name(save_path):
    """

    Parameters
    ----------
    save_path

    Returns
    -------

    """
    i = 0
    store_type = save_path.split(".")[-1]
    store_path = save_path[0 : -len(store_type) - 1]
    save_name = f"{store_path}_{i:03d}.{store_type}"
    while os.path.isfile(save_name):
        i = i + 1
        save_name = f"{store_path}_{i:03d}.{store_type}"
    return save_name


def save_plot_renderer(plot_obj, store_path, renderer):
    """
    Save a plot using a renderer.

    Parameters
    ----------
    plot_obj
    store_path
    renderer

    Returns
    -------

    """
    save = get_save_name(f"{store_path}.png")
    renderer.save(plot_obj, save)


def save_plot(save_path, ax_or_fig):
    """

    Parameters
    ----------
    save_path
    ax

    Returns
    -------
    New values for 'save_path' (None) and 'save' (False)
    """
    try:
        save_name = get_save_name(save_path)
        if isinstance(ax_or_fig, plt.Figure):
            ax_or_fig.savefig(save_name, bbox_inches="tight", dpi=300)
        else:
            ax_or_fig.figure.savefig(save_name, bbox_inches="tight", dpi=300)
    except OSError as exception:
        print(
            f"Could not save to {save_path}. Did you forget the filetype?",
            file=sys.stderr,
        )
        print(exception.errno, file=sys.stderr)
    return None, False


# pylint: disable=invalid-name
def format_tick(v, pos, max_tick):
    """

    Parameters
    ----------
    v
    pos
    max_tick

    Returns
    -------

    """
    if pos:
        upper = v - max_tick
        return f"1e{upper:2.1f}"
    upper = max_tick - v
    return f"-1e{upper:2.1f}"
