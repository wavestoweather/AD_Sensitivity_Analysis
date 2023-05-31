"""Functions to create interactive plots with Panel.

To be used with Jupyter notebooks.
"""
import warnings

import panel as pn
import seaborn as sns

from ad_sensitivity_analysis.interactive import create_widgets
from ad_sensitivity_analysis.plot import project_map

warnings.simplefilter(action="ignore", category=RuntimeWarning)


# pylint: disable=too-many-locals
def plot_2dmap_interactive(ds):
    """
    Calling this function from a Jupyter notebook allows to visualize the
    2d maps at different heights for all available columns in ds.

    Parameters
    ----------
    ds : xarray.Dataset

    Returns
    -------

    """
    sns.set_style("darkgrid")

    pressure = create_widgets.create_float_slider(ds, "pressure")
    time_slider = create_widgets.create_time_slider_select(ds)
    out_param = create_widgets.create_out_params_buttons(ds)
    kind_param, in_param = create_widgets.create_kind_in_params_selector(ds)
    color_map = create_widgets.create_color_select()

    fix = pn.widgets.Toggle(
        name="Fix colorbar over all levels",
        button_type="success",
    )
    fix_time = pn.widgets.Toggle(
        name="Fix colorbar over all time steps",
        button_type="success",
    )
    width_height_font_row = create_widgets.image_font_size_row(pixel=False)

    log_plot = pn.widgets.Toggle(
        name="Use log colorbar",
        value=False,
        button_type="success",
    )
    log_threshold_slider = pn.widgets.FloatSlider(
        name="Log Threshold (set to zero for automatic estimation)",
        start=-25,
        end=0,
        value=0,
        step=1,
    )
    title_widget = pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )
    color_blind_toggle = pn.widgets.Toggle(
        name="Adjust for Colorblind",
        value=False,
        button_type="success",
    )

    save_field_button_latex_widget = create_widgets.create_save_latex_row()

    plot_pane = pn.panel(
        pn.bind(
            project_map.plot_2dmap,
            ds=ds,
            out_param=out_param,
            in_param=in_param,
            kind_param=kind_param,
            pressure=pressure,
            time=time_slider,
            cmap=color_map,
            fix=fix,
            fix_time=fix_time,
            log_plot=log_plot,
            width=width_height_font_row[0],
            height=width_height_font_row[1],
            lthresh=log_threshold_slider,
            title=title_widget,
            font_scale=width_height_font_row[2],
            save=save_field_button_latex_widget[1],
            latex=save_field_button_latex_widget[2],
            save_path=save_field_button_latex_widget[0],
            colorblind=color_blind_toggle,
        ),
    ).servable()

    return pn.Column(
        "# Plot Grid",
        kind_param,
        out_param,
        pn.Row(
            in_param,
            color_map,
        ),
        pn.Row(
            fix,
            log_plot,
            color_blind_toggle,
        ),
        width_height_font_row,
        log_threshold_slider,
        save_field_button_latex_widget,
        title_widget,
        pn.Row(
            pressure,
            time_slider,
        ),
        plot_pane,
    )
