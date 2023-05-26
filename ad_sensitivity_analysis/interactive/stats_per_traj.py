"""Interactive plots for statistics over trajectories.

"""
import warnings

import panel as pn

from ad_sensitivity_analysis.interactive import create_widgets
from ad_sensitivity_analysis.plot import stats_per_traj_static

warnings.simplefilter(action="ignore", category=RuntimeWarning)


def plot_kde_histogram_interactive(ds):
    """
    Use plot_kde_histogram() interactively in a jupyter notebook.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset created via create_rank_traj_dataset().

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    out_param = create_widgets.create_out_params_select(ds)
    width_height_font_row = create_widgets.image_font_size_row(pixel=False)

    title_widget = create_widgets.create_title_widget()
    save_field_button_latex_row = create_widgets.create_save_latex_row()

    in_params = create_widgets.create_in_params_selector(ds)
    line_slider = create_widgets.create_line_slider()
    bw_slider = pn.widgets.FloatSlider(
        name="Change the bandwidth for the kde calculation",
        start=0.05,
        end=1.5,
        step=0.05,
        value=1.0,
    )
    flow_phase_ignore_row = create_widgets.create_flow_phase_ignore_row()
    common_toggle = pn.widgets.Toggle(
        name="Common norm",
        button_type="success",
    )

    plot_pane = pn.panel(
        pn.bind(
            stats_per_traj_static.plot_kde_histogram,
            ds=ds,
            in_params=in_params,
            out_param=out_param,
            linewidth=line_slider,
            bw_adjust=bw_slider,
            flow=flow_phase_ignore_row[0],
            phase=flow_phase_ignore_row[1],
            ignore_zero_gradients=flow_phase_ignore_row[2],
            title=title_widget,
            filename=save_field_button_latex_row[0],
            width=width_height_font_row[0],
            height=width_height_font_row[1],
            common_norm=common_toggle,
            font_scale=width_height_font_row[2],
            save=save_field_button_latex_row[1],
            latex=save_field_button_latex_row[2],
        ),
    ).servable()

    return pn.Column(
        width_height_font_row,
        save_field_button_latex_row,
        flow_phase_ignore_row,
        pn.Row(
            in_params,
            pn.Column(
                out_param,
                common_toggle,
            ),
        ),
        pn.Row(
            line_slider,
            bw_slider,
            title_widget,
        ),
        plot_pane,
    )
