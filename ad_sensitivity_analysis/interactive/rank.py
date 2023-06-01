"""Visualize rank and impact interactively.

"""
import panel as pn

from ad_sensitivity_analysis.plot.rank import plot_rank_over_impact, plot_rank_probs


# pylint: disable=too-many-locals
def plot_rank_over_impact_interactive(ds):
    """

    Parameters
    ----------
    ds

    Returns
    -------

    """
    out_params = []
    for o_p in ds["Output Parameter"]:
        out_params.append(o_p.item())
    out_params.append("all")
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_params[0],
        options=out_params,
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
    dot_slider = pn.widgets.IntSlider(
        name="Dot size",
        start=1,
        end=500,
        step=2,
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
    log_y = pn.widgets.Toggle(
        name="Log y-axis",
        value=False,
        button_type="success",
    )
    trajectory_toggle = pn.widgets.Toggle(
        name="Average over all trajectories",
        value=True,
        button_type="success",
    )
    phase_toggle = pn.widgets.Toggle(
        name="Show phases",
        value=False,
        button_type="success",
    )
    flow_toggle = pn.widgets.Toggle(
        name="Show flow",
        value=False,
        button_type="success",
    )
    color_blind_toggle = pn.widgets.Toggle(
        name="Adjust for Colorblind",
        value=True,
        button_type="success",
    )
    plot_pane = pn.panel(
        pn.bind(
            plot_rank_over_impact,
            ds=ds,
            out_param=out_param,
            width=width_slider,
            height=height_slider,
            font_scale=font_slider,
            filename=save_to_field,
            title=title_widget,
            save=save_button,
            latex=latex_button,
            dot_size=dot_slider,
            log=log_y,
            phase=phase_toggle,
            flow=flow_toggle,
            avg=trajectory_toggle,
            colorblind=color_blind_toggle,
        ),
    ).servable()
    return pn.Column(
        pn.Row(
            width_slider,
            height_slider,
            dot_slider,
        ),
        pn.Row(
            phase_toggle,
            flow_toggle,
            color_blind_toggle,
        ),
        pn.Row(
            font_slider,
            log_y,
            trajectory_toggle,
        ),
        pn.Row(
            save_to_field,
            save_button,
            latex_button,
        ),
        pn.Row(
            out_param,
            title_widget,
        ),
        plot_pane,
    )


# pylint: disable=too-many-locals
def plot_rank_probs_interactive(ds):
    """

    Parameters
    ----------
    ds

    Returns
    -------

    """
    out_params = []
    for o_p in ds["Output Parameter"]:
        out_params.append(o_p.item())
    out_params.append("all")
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_params[0],
        options=out_params,
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
    flow_values = []
    for flow in ds["flow"].values:
        if "neutral" not in flow:
            flow_values.append(flow)
    flow_select = pn.widgets.Select(
        name="Flow",
        value=flow_values[0],
        options=flow_values,
    )
    phase_select = pn.widgets.Select(
        name="Phase",
        value=ds["phase"].values[0],
        options=ds["phase"].values.tolist(),
    )
    dot_slider = pn.widgets.IntSlider(
        name="Change the dot size",
        start=1,
        end=500,
        step=2,
        value=12,
    )
    median_toggle = pn.widgets.Toggle(
        name="Show median instead of mean on y-axis",
        button_type="success",
        value=True,
    )
    rank_toggle = pn.widgets.Toggle(
        name="Show rank instead of impact on y-axis",
        button_type="success",
    )
    top_slider = pn.widgets.IntSlider(
        name="Color the top n parameters",
        start=1,
        end=38,
        step=1,
        value=10,
    )
    log_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=False,
        button_type="success",
    )
    color_blind_toggle = pn.widgets.Toggle(
        name="Adjust for colorblind",
        value=True,
        button_type="success",
    )
    cbar_toggle = pn.widgets.Toggle(
        name="Show colorbar",
        value=True,
        button_type="success",
    )

    plot_pane = pn.panel(
        pn.bind(
            plot_rank_probs,
            ds=ds,
            flow=flow_select,
            phase=phase_select,
            out_param=out_param,
            median=median_toggle,
            rank=rank_toggle,
            mark_top_n=top_slider,
            dot_size=dot_slider,
            title=title_widget,
            filename=save_to_field,
            width=width_slider,
            height=height_slider,
            font_scale=font_slider,
            save=save_button,
            latex=latex_button,
            logy=log_plot,
            colorblind=color_blind_toggle,
            show_cbar=cbar_toggle,
        ),
    ).servable()

    return pn.Column(
        pn.Row(
            width_slider,
            height_slider,
        ),
        pn.Row(
            color_blind_toggle,
            cbar_toggle,
        ),
        pn.Row(
            save_to_field,
            save_button,
            latex_button,
        ),
        pn.Row(
            flow_select,
            phase_select,
            out_param,
        ),
        pn.Row(
            median_toggle,
            rank_toggle,
            log_plot,
        ),
        pn.Row(
            top_slider,
            font_slider,
            dot_slider,
        ),
        plot_pane,
    )
