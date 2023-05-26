"""Interactive functions for clustered datasets and notebooks.

"""

import panel as pn

from ad_sensitivity_analysis.plot.clustered_data import plot_cluster_data


# pylint: disable=too-many-locals
def plot_cluster_data_interactive(data, reduce_name=""):
    """
    Calling this function from a Jupyter notebook allows to visualize the cluster association with different
    dimensions. Make sure to call pn.extension() from your notebook first.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame generated using get_cluster().
    reduce_name : string
        Name prepended to the columns. Should relate to the reduction
        applied to the dataset, such as "avg" or "rank" in get_cluster().
    Returns
    -------

    """
    out_params = []
    in_params = []
    for col in data:
        if "/" in col:
            out_params.append(col.split("/")[0][len(reduce_name) + 1 :])
            in_params.append(col.split("/")[1])
        elif reduce_name in col and len(reduce_name) > 0:
            in_params.append(col[len(reduce_name) :])
        else:
            in_params.append(col)
    in_params = list(set(in_params))
    in_params.sort()
    out_params = list(set(out_params))
    out_params.sort()
    if len(out_params) == 0:
        out_params = ["Not available"]
    out_param_x = pn.widgets.RadioButtonGroup(
        name="Output Parameter (if any) for the x-axis",
        value=out_params[0],
        options=out_params,
        button_type="primary",
    )
    in_param_x = pn.widgets.Select(
        name="Model parameter or model state for the x-axis",
        value=in_params[0],
        options=in_params,
    )
    out_param_y = pn.widgets.RadioButtonGroup(
        name="Output Parameter (if any) for the y-axis",
        value=out_params[0],
        options=out_params,
        button_type="primary",
    )
    in_param_y = pn.widgets.Select(
        name="Model parameter or model state for the y-axis",
        value=in_params[1],
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
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    logx_plot = pn.widgets.Toggle(
        name="Use log-scale for the x-axis",
        value=False,
        button_type="success",
    )
    logy_plot = pn.widgets.Toggle(
        name="Use log-scale for the y-axis",
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
    dot_slider = pn.widgets.IntSlider(
        name="Change the dot size",
        start=1,
        end=200,
        step=2,
        value=12,
    )

    plot_pane = pn.panel(
        pn.bind(
            plot_cluster_data,
            data=data,
            in_p_x=in_param_x,
            out_p_x=out_param_x,
            in_p_y=in_param_y,
            out_p_y=out_param_y,
            logx=logx_plot,
            logy=logy_plot,
            width=width_slider,
            height=height_slider,
            font_scale=font_slider,
            title=title_widget,
            save_path=save_to_field,
            latex=latex_button,
            save=save_button,
            s=dot_slider,
        ),
    ).servable()

    return pn.Column(
        pn.Row(
            out_param_x,
            in_param_x,
            logx_plot,
        ),
        pn.Row(
            out_param_y,
            in_param_y,
            logy_plot,
        ),
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
        ),
        pn.Row(
            dot_slider,
            latex_button,
        ),
        pn.Row(
            save_to_field,
            save_button,
        ),
        title_widget,
        plot_pane,
    )
