"""Interactive plots of histograms.

"""
import numpy as np
import panel as pn

from ad_sensitivity_analysis.plot.stats_histogram import (
    plot_heatmap_histogram,
    traj_plot_histogram_inp,
    traj_plot_kde_inp,
    traj_plot_histogram_out,
)


def plot_traj_histogram_out_interactive(edges, hist):
    """
    Calling this function from a Jupyter notebook allows to visualize the
    traj_plot_histogram_out interactively. Plots the histogram of an output parameter,
    i.e., QV, latent_heat, latent_cool, etc.

    Parameters
    ----------
    edges : Dictionary of list-like of float
        Edges for the histogram. Keys must be in out_params.
    hist : Dictionary of list-like of int
        Number of entries for each bin. Keys must be in edges.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=list(edges.keys())[0],
        options=list(edges.keys()),
    )
    width_slider = pn.widgets.IntSlider(
        name="Width in inches",
        start=3,
        end=15,
        step=1,
        value=9,
    )
    height_slider = pn.widgets.IntSlider(
        name="Height in inches",
        start=3,
        end=15,
        step=1,
        value=6,
    )
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    log_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=True,
        button_type="success",
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    tick_slider = pn.widgets.IntSlider(
        name="Plot every n ticks on the x-axis:",
        start=1,
        end=20,
        step=1,
        value=6,
    )
    plot_pane = pn.panel(
        pn.bind(
            traj_plot_histogram_out,
            out_params=out_param,
            filename=save_to_field,
            edges=edges,
            hist=hist,
            log=log_plot,
            width=width_slider,
            height=height_slider,
            title=title_widget,
            save=save_button,
            interactive=True,
            font_scale=font_slider,
            latex=latex_button,
            ticks_offset=tick_slider,
            verbose=False,
        ),
    ).servable()

    return pn.Column(
        out_param,
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
            tick_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            log_plot,
            latex_button,
        ),
        title_widget,
        plot_pane,
    )


def plot_traj_histogram_inp_interactive(
    edges_in_params,
    hist_in_params,
):
    r"""
    Can be used in jupyter notebooks to interactively plot traj_plot_histogram_inp(). From traj_plot_histogram_inp():
    Giuen histograms from a sensitivity analysis with multiple trajectories, plot three histograms per image with
    \partial output / \partial model parameter
    where output is QV, latent_heat and latent_cool. Plot one image per model_parameter.

    Parameters
    ----------
    edges_in_params : Dictionary of dictionary list-like of float
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        the bin edges for the given keys. Optional: The first level can be another dictionary of phases.
    hist_in_params : Dictionary of dictionary list-like of int
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        values of the histogram for the given key. Optional: The first level can be another dictionary of phases.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    tmp_params = list(edges_in_params[list(edges_in_params.keys())[0]].keys())
    in_params = []
    for param in tmp_params:
        if param[0] == "d" and param != "deposition":
            in_params.append(param)
    in_params = list(np.sort(in_params))
    in_param = pn.widgets.Select(
        name="Model Parameter",
        value=in_params[0],
        options=in_params,
    )
    width_slider = pn.widgets.IntSlider(
        name="Width in inches",
        start=3,
        end=15,
        step=1,
        value=9,
    )
    height_slider = pn.widgets.IntSlider(
        name="Height in inches",
        start=3,
        end=15,
        step=1,
        value=6,
    )
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    log_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=True,
        button_type="success",
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    tick_slider = pn.widgets.IntSlider(
        name="Plot every n ticks on the x-axis:",
        start=1,
        end=20,
        step=1,
        value=6,
    )

    plot_pane = pn.panel(
        pn.bind(
            traj_plot_histogram_inp,
            filename=save_to_field,
            in_params=in_param,
            edges_in_params=edges_in_params,
            hist_in_params=hist_in_params,
            log=log_plot,
            width=width_slider,
            height=height_slider,
            title=None,
            font_scale=font_slider,
            save=save_button,
            interactive=True,
            latex=latex_button,
            ticks_offset=tick_slider,
            verbose=False,
        ),
    ).servable()

    return pn.Column(
        in_param,
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
            tick_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            log_plot,
            latex_button,
        ),
        plot_pane,
    )


# pylint: disable=too-many-locals
def plot_traj_kde_inp_interactive(
    edges_in_params,
    hist_in_params,
):
    """
    Can be used in jupyter notebooks to interactively plot traj_plot_kde_inp().

    Parameters
    ----------
    edges_in_params : Dictionary of dictionary list-like of float
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        the bin edges for the given keys. Optional: The first level can be another dictionary of phases.
    hist_in_params : Dictionary of dictionary list-like of int
        Dictionary (keys = model state variable) of dictionaries (keys = model parameters) of arrays with
        values of the histogram for the given key. Optional: The first level can be another dictionary of phases.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    tmp_params = list(edges_in_params[list(edges_in_params.keys())[0]].keys())
    in_params = []
    for param in tmp_params:
        if param[0] == "d" and param != "deposition":
            in_params.append(param)
    in_params.sort()
    in_param = pn.widgets.CrossSelector(
        name="Parameter",
        value=in_params[0:2],
        options=in_params,
    )
    out_param_list = list(edges_in_params.keys())
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_param_list[0],
        options=out_param_list,
    )
    width_slider = pn.widgets.IntSlider(
        name="Width in inches",
        start=3,
        end=15,
        step=1,
        value=9,
    )
    height_slider = pn.widgets.IntSlider(
        name="Height in inches",
        start=3,
        end=15,
        step=1,
        value=6,
    )
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    line_slider = pn.widgets.FloatSlider(
        name="Change the line width",
        start=1,
        end=10,
        step=0.5,
        value=2,
    )
    bw_slider = pn.widgets.FloatSlider(
        name="Change the bandwidth for the kde calculation",
        start=0.05,
        end=1.5,
        step=0.05,
        value=0.1,
    )
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    log_x = pn.widgets.Toggle(
        name="Log x-axis",
        value=False,
        button_type="success",
    )
    log_weights = pn.widgets.Toggle(
        name="Log weights",
        value=True,
        button_type="success",
    )

    plot_pane = pn.panel(
        pn.bind(
            traj_plot_kde_inp,
            filename=save_to_field,
            in_params=in_param,
            out_param=out_param,
            linewidth=line_slider,
            bw_adjust=bw_slider,
            edges_in_params=edges_in_params,
            hist_in_params=hist_in_params,
            log_weights=log_weights,
            log_x=log_x,
            width=width_slider,
            height=height_slider,
            title=title_widget,
            font_scale=font_slider,
            filter_zero=False,
            save=save_button,
            latex=latex_button,
        ),
    ).servable()

    return pn.Column(
        pn.Row(
            in_param,
            pn.Column(
                out_param,
                log_x,
                log_weights,
            ),
        ),
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            latex_button,
        ),
        pn.Row(
            line_slider,
            bw_slider,
            title_widget,
        ),
        plot_pane,
    )


# pylint: disable=too-many-locals
def plot_heatmap_histogram_interactive(hist_conditional):
    """
    Can be used in jupyter notebooks to interactively plot plot_heatmap_histogram() (2D histograms).

    Parameters
    ----------
    hist_conditional : Dictionary of dictionaries with edges and bin counts for 2D histograms.
        Result of get_histogram_cond().
        Dictionary with the following keys:
        'edges_out_params': Dictionary where the keys are model state variables and the values are arrays of bin edges.
        'edges_in_params': Dictionary where the keys are model state variables for which sensitivities are available
            and the values are dictionaries of model parameters with arrays of bin edges.
        model state variables: Each model state variable has a dictionary for 'hist_out_params' and 'hist_in_params'.
        'hist_out_params' is a dictionary of model state variables with arrays of bin counts.
        'hist_in_params' is a dictionary of model state variables for which sensitivities are available
            and the values are dictionaries of model parameters with arrays of bin counts.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    conds = list(hist_conditional.keys())
    conditions = []
    for c in conds:
        if c not in ("edges_in_params", "edges_out_params"):
            conditions.append(c)
    conditions.sort()
    wrt_params = list(hist_conditional["edges_in_params"].keys())
    out_param = pn.widgets.RadioButtonGroup(
        name="Output Parameter",
        value=wrt_params[0],
        options=wrt_params,
        button_type="primary",
    )
    condition = pn.widgets.Select(
        name="X-Axis",
        value=conditions[0],
        options=conditions,
    )
    in_params = ["None"]
    tmp_params = list(hist_conditional["edges_in_params"][wrt_params[0]].keys())
    for param in tmp_params:
        if param[0] == "d" and param != "deposition":
            in_params.append(param)
    in_params.sort()
    in_param = pn.widgets.Select(
        name="Model Parameter (Y-Axis)",
        value=in_params[-1],
        options=in_params,
    )
    width_slider = pn.widgets.IntSlider(
        name="Width in inches",
        start=3,
        end=20,
        step=1,
        value=8,
    )
    height_slider = pn.widgets.IntSlider(
        name="Height in inches",
        start=3,
        end=20,
        step=1,
        value=8,
    )
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    save_to_field = pn.widgets.TextInput(
        value="Path/to/store/plot.png",
    )
    save_button = pn.widgets.Button(
        name="Save Plot",
        button_type="primary",
    )
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    log_plot = pn.widgets.Toggle(
        name="Use log colorbar",
        value=True,
        button_type="success",
    )
    latex_button = pn.widgets.Toggle(
        name="Latexify",
        value=False,
        button_type="success",
    )
    tick_slider = pn.widgets.IntSlider(
        name="Plot every n ticks:",
        start=1,
        end=20,
        step=1,
        value=6,
    )

    plot_pane = pn.panel(
        pn.bind(
            plot_heatmap_histogram,
            hist_conditional=hist_conditional,
            filename=save_to_field,
            in_params=in_param,
            out_params=out_param,
            conditions=condition,
            log=log_plot,
            width=width_slider,
            height=height_slider,
            font_scale=font_slider,
            title=title_widget,
            save=save_button,
            latex=latex_button,
            interactive=True,
            ticks_offset=tick_slider,
            verbose=False,
        ),
    ).servable()

    return pn.Column(
        in_param,
        "w.r.t. Model State (Y-Axis)",
        out_param,
        condition,
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
            tick_slider,
        ),
        pn.Row(
            save_to_field,
            save_button,
            log_plot,
            latex_button,
        ),
        title_widget,
        plot_pane,
    )
