"""Interactive plots of covariance or correlation matrix.

"""
import pickle

import numpy as np
import matplotlib.colors as mpl_col
import panel as pn

from ad_sensitivity_analysis.plot.matrix import plot_heatmap_matrix


# pylint: disable=too-many-locals, too-many-branches
def plot_heatmap_interactive(
    filepath=None,
    matrix_dic=None,
    names=None,
):
    """
    Interactive plot for plot_heatmap_matrix() that takes correlation and covariance matrices.

    Parameters
    ----------
    filepath : string (optional)
        Path to pickle file with a correlation matrix to load.
    matrix_dic : dictionary of matrices or dictionary of dictionary of dictionary of matrices
        Either: The keys are model state variables with available sensitivities.
        The values are matrices where covariances are given for sensitivities for the respective model state variable.
        The data can be generated using get_cov_matrix().
        Or: The first keys are model state variables, then the flow, then phase used for calculating
        correlation/covariance values.
        The data can be generated using get_stats_per_traj.get_matrix().
    names : list of strings
        The names of each column/row in the covariance matrix. The index of names is the row/column of the matrix.
    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    if filepath is not None:
        with open(filepath, "rb") as f:
            data = pickle.load(f)
            matrix_dic = data["correlation"]
            names = data["column names"]
    elif matrix_dic is None or names is None:
        print(
            "You either have to specify a filepath to a pickle file or "
            "you need to provide the matrix dictionary with the names."
        )
        return None
    out_param_list = list(matrix_dic.keys())
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_param_list[0],
        options=out_param_list,
    )

    width_slider = pn.widgets.IntSlider(
        name="Width in inches",
        start=3,
        end=45,
        step=1,
        value=9,
    )
    height_slider = pn.widgets.IntSlider(
        name="Height in inches",
        start=3,
        end=25,
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
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    cmaps = [
        "viridis",
        "mako",
        "vlag",
        "icefire",
        "Spectral",
        "coolwarm",
        "Blues",
        "Reds",
    ]
    color_widget = pn.widgets.Select(
        name="Colormap",
        value=cmaps[2],
        options=cmaps,
    )
    sort_button = pn.widgets.Toggle(
        name="Sort by correlation sum",
        button_type="success",
    )
    rank_names = False
    for n in names:
        if "rank" in n:
            rank_names = True
            break
    if rank_names:
        names_list = []
        for n in np.sort(names):
            if "rank" in n:
                names_list.append(n)
        for n in np.sort(names):
            if "rank" not in n:
                names_list.append(n)
    else:
        names_list = list(np.sort(names))
    param_select = pn.widgets.CrossSelector(
        name="Parameter",
        value=names_list[0:5],
        options=names_list,
    )
    plot_type = pn.widgets.Select(
        name="Plot values:",
        value="all",
        options=["all", "negative", "positive"],
    )
    norms = [
        mpl_col.Normalize(-1, 1),
        mpl_col.Normalize(1, 0),
        mpl_col.Normalize(0, 1),
        mpl_col.SymLogNorm(1e-2),
        mpl_col.SymLogNorm(1e-10),
        mpl_col.SymLogNorm(1e-15),
        mpl_col.SymLogNorm(1e-20),
        mpl_col.LogNorm(),
        None,
    ]
    norm = pn.widgets.Select(
        name="Norm",
        value=norms[0],
        options=norms,
    )
    annot_button = pn.widgets.Toggle(
        name="Annotate entries",
        button_type="success",
    )
    if isinstance(matrix_dic[out_param_list[0]], dict):
        flow_list = list(matrix_dic[out_param_list[0]].keys())
        flow = pn.widgets.Select(
            name="Flow",
            value=flow_list[0],
            options=flow_list,
        )
        phase_list = list(matrix_dic[out_param_list[0]][flow_list[0]].keys())
        phase = pn.widgets.Select(
            name="Phase",
            value=phase_list[0],
            options=phase_list,
        )
    else:
        phase = None
        flow = None

    plot_pane = pn.panel(
        pn.bind(
            plot_heatmap_matrix,
            data_in=matrix_dic,
            out_param=out_param,
            names=names,
            in_params=param_select,
            plot_type=plot_type,
            norm=norm,
            title=title_widget,
            filename=save_to_field,
            cmap=color_widget,
            width=width_slider,
            height=height_slider,
            sort=sort_button,
            latex=latex_button,
            save=save_button,
            phase=phase,
            flow=flow,
            font_scale=font_slider,
            annot=annot_button,
        )
    ).servable()

    if isinstance(matrix_dic[out_param_list[0]], dict):
        return pn.Column(
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
                param_select,
                pn.Column(
                    sort_button,
                    flow,
                    phase,
                ),
            ),
            pn.Row(
                pn.pane.Markdown(
                    "Available norms: Normalized (-1, 1), Normalized (1, 0), Normalizes (0, -1), "
                    "SymLogNorm(thresholds: 1e-2, 1e-10, 1e-15, or 1e-20), and LogNorm."
                ),
            ),
            pn.Row(
                out_param,
                plot_type,
                norm,
            ),
            pn.Row(
                title_widget,
                color_widget,
                annot_button,
            ),
            plot_pane,
        )
    return pn.Column(
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
            param_select,
            sort_button,
        ),
        pn.Row(
            out_param,
            plot_type,
            pn.pane.Markdown(
                "Available norms: Normalized (-1, 1), SymLogNorm(1e-2, 1e-10, 1e-15, or 1e-20), and LogNorm."
            ),
        ),
        pn.Row(
            title_widget,
            color_widget,
            norm,
        ),
        plot_pane,
    )
