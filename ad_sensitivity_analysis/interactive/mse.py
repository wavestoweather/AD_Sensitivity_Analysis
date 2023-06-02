"""Interactive plots for plotting predicted errors.

"""
import numpy as np
import pandas as pd
import panel as pn

from ad_sensitivity_analysis.plot.latexify import param_id_map
from ad_sensitivity_analysis.plot.mse_time_evolution import plot_time_evolution
from ad_sensitivity_analysis.plot.mse import plot_errors

# pylint: disable=too-many-locals
def plot_errors_interactive(
    df,
):
    """
    Use plot_errors() interactively in a jupyter notebook for perturbation simulations.

    Parameters
    ----------
    df : pandas.DataFrame
        A dataframe with columns "Output Parameter" for the model state
        "Input Parameter" for the perturbed model parameter,
        "Predicted Squared Error" and "Predicted Error" for the sensitivity calculated to deviations
        in the next timestep, "Mean Squared Error" or "Mean Error" for the actual deviations

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    out_param_list = list(np.unique(df["Output Parameter"]))
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_param_list[0],
        options=out_param_list,
    )
    in_params_list = list(np.unique(df["Input Parameter"]))
    in_params = pn.widgets.CrossSelector(
        name="Parameter",
        value=in_params_list[0:2],
        options=in_params_list,
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
    alpha_slider = pn.widgets.FloatSlider(
        name="Alpha",
        start=0.1,
        end=1,
        step=0.1,
        value=0.6,
    )
    dot_slider = pn.widgets.IntSlider(
        name="Change the dot size",
        start=1,
        end=200,
        step=2,
        value=12,
    )
    logx_plot = pn.widgets.Toggle(
        name="Use log x-axis",
        value=False,
        button_type="success",
    )
    logy_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=False,
        button_type="success",
    )
    x_widget = pn.widgets.TextInput(
        name="X-label",
        placeholder="Predicted Error",
        value="Predicted Error",
    )
    y_widget = pn.widgets.TextInput(
        name="Y-label",
        placeholder="Ensemble Error",
        value="Ensemble Error",
    )
    data_variants = []
    for col in df:
        if col not in ("Output Parameter", "Input Parameter"):
            data_variants.append(col)
    x_data = pn.widgets.Select(
        name="X-axis",
        value=data_variants[0],
        options=data_variants,
    )
    y_data = pn.widgets.Select(
        name="Y-axis",
        value=data_variants[2],
        options=data_variants,
    )
    corr_line = pn.widgets.Toggle(
        name="Correlation line",
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
    group_toggle = pn.widgets.Toggle(
        name="Group parameters",
        value=False,
        button_type="success",
    )
    ellipsis_widget = pn.widgets.FloatSlider(
        name="Ellipsis in standard deviations",
        start=0,
        end=5,
        step=0.2,
        value=2,
    )

    plot_pane = pn.panel(
        pn.bind(
            plot_errors,
            df=df,
            out_param=out_param,
            in_params=in_params,
            x_key=x_data,
            y_key=y_data,
            alpha=alpha_slider,
            plot_types=group_toggle,
            n_std=ellipsis_widget,
            linewidth=line_slider,
            title=title_widget,
            xlabel=x_widget,
            ylabel=y_widget,
            width=width_slider,
            height=height_slider,
            log_x=logx_plot,
            log_y=logy_plot,
            corr_line=corr_line,
            font_scale=font_slider,
            save=save_button,
            latex=latex_button,
            filename=save_to_field,
            s=dot_slider,
        ),
    ).servable()
    return pn.Column(
        pn.Row(
            width_slider,
            height_slider,
            font_slider,
        ),
        pn.Row(
            x_data,
            y_data,
            ellipsis_widget,
        ),
        pn.Row(
            save_to_field,
            save_button,
            latex_button,
        ),
        pn.Row(
            logx_plot,
            logy_plot,
            corr_line,
        ),
        pn.Row(
            in_params,
            pn.Column(
                out_param,
                x_widget,
                y_widget,
                group_toggle,
            ),
        ),
        pn.Row(
            dot_slider,
            line_slider,
            alpha_slider,
        ),
        title_widget,
        plot_pane,
    )


# pylint: disable=too-many-arguments, too-many-locals
def _prepare_dataset_and_plot(
    ds,
    backend,
    store_path,
    title,
    xlabel,
    ylabel,
    twinlabel,
    logtwin,
    logy,
    width,
    height,
    x_limits,
    save,
    latex,
    font_scale,
    dot_size,
    trajectory,
    out_param,
    in_params,
    plot_deviation,
    perturbed,
    precision,
):
    """

    Parameters
    ----------
    ds
    backend
    store_path
    title
    xlabel
    ylabel
    twinlabel
    logtwin
    logy
    width
    height
    x_limits
    save
    latex
    font_scale
    dot_size
    trajectory
    out_param
    in_params
    plot_deviation
    perturbed
    precision

    Returns
    -------

    """
    if not perturbed:
        out_p_id = np.argwhere(np.asarray(param_id_map) == out_param).item()
        ds_traj = ds.isel({"ensemble": 0, "trajectory": trajectory}).sel(
            {"Output_Parameter_ID": out_p_id}
        )
        # Sensitivity simulation
        predicted_errors = np.array([])
        in_param_thingy = np.array([])
        unperturbed = np.array([])
        time = np.array([])
        for col in in_params:
            predicted_errors = np.append(predicted_errors, ds_traj[col].values)
            in_param_thingy = np.append(
                in_param_thingy, np.repeat(col, len(ds_traj[col].values))
            )
            unperturbed = np.append(unperturbed, ds_traj[out_param].values)
            time = np.append(time, ds_traj["time_after_ascent"].values)
        df = pd.DataFrame(
            data={
                "Predicted Error": predicted_errors,
                "Predicted Squared Error": predicted_errors**2,
                "Output Parameter": np.repeat(out_param, len(predicted_errors)),
                "Input Parameter": in_param_thingy,
                "Not Perturbed Value": unperturbed,
                "time_after_ascent": time,
            }
        )
        return plot_time_evolution(
            df=df,
            backend=backend,
            store_path=store_path,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            twinlabel=twinlabel,
            logtwin=logtwin,
            logy=logy,
            width=width,
            height=height,
            x_limits=x_limits,
            save=save,
            latex=latex,
            font_scale=font_scale,
            dot_size=dot_size,
            plot_deviation=plot_deviation,
            precision=precision,
        )
    # Perturbed ensemble
    df = ds.to_dataframe().reset_index()
    return plot_time_evolution(
        df=df,
        backend=backend,
        store_path=store_path,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        twinlabel=twinlabel,
        logtwin=logtwin,
        logy=logy,
        width=width,
        height=height,
        x_limits=x_limits,
        save=save,
        latex=latex,
        dot_size=dot_size,
        trajectory=trajectory,
        font_scale=font_scale,
        out_param=out_param,
        in_params=in_params,
        plot_deviation=plot_deviation,
        precision=precision,
    )


# pylint: disable=too-many-locals
def plot_time_evolution_interactive(ds, perturbed=False):
    """
    Use plot_time_evolution() interactively in a jupyter notebook for perturbation or sensitivity simulations.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with trajectories from a perturbed or sensitivity simulation.
    perturbed : bool
        If true, then ds is from a perturbed ensemble simulation. Otherwise a sensitivity simulation is
        assumed.

    Returns
    -------
    panel.layout that can be used in a jupyter notebook.
    """
    if "Output Parameter" in ds:
        out_param_list = ds["Output Parameter"].values.tolist()
    else:
        out_param_list = []
        for out_p in ds["Output_Parameter_ID"]:
            out_param_list.append(param_id_map[out_p.item()])
    out_param = pn.widgets.Select(
        name="Output Parameter",
        value=out_param_list[0],
        options=out_param_list,
    )
    in_params_list = []
    for col in ds:
        if col[0] == "d" and col != "deposition":
            in_params_list.append(col)
    in_params_list = list(np.sort(in_params_list))
    in_params = pn.widgets.CrossSelector(
        name="Parameter",
        value=in_params_list[0:2],
        options=in_params_list,
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
    dot_slider = pn.widgets.IntSlider(
        name="Change the dot size",
        start=1,
        end=200,
        step=2,
        value=12,
    )
    x_slider = pn.widgets.RangeSlider(
        name="X Limits in percent",
        start=0,
        end=1,
        value=(0, 1),
        step=0.001,
    )
    log_twin_plot = pn.widgets.Toggle(
        name="Use log twin y-axis",
        value=False,
        button_type="success",
    )
    log_plot = pn.widgets.Toggle(
        name="Use log y-axis",
        value=False,
        button_type="success",
    )
    x_widget = pn.widgets.TextInput(
        name="X-label",
        placeholder="Time after ascent [min]",
        value="Time after ascent [min]",
    )
    y_widget = pn.widgets.TextInput(
        name="Y-label",
        placeholder="Specific humidty [kg/kg]",
    )
    y2_widget = pn.widgets.TextInput(
        name="Twin y-label",
        placeholder="Predicted Deviation [kg/kg]",
        value="Predicted Deviation [kg/kg]",
    )
    traj_widget = pn.widgets.IntSlider(
        name="Trajectory",
        start=0,
        end=len(ds["trajectory"]) - 1,
        step=1,
    )
    prec_widget = pn.widgets.IntSlider(
        name="Precision",
        start=0,
        end=5,
        step=1,
    )

    plot_pane = pn.panel(
        pn.bind(
            _prepare_dataset_and_plot,
            ds=ds,
            backend="matplotlib",
            store_path=save_to_field,
            title=title_widget,
            xlabel=x_widget,
            ylabel=y_widget,
            twinlabel=y2_widget,
            logtwin=log_twin_plot,
            logy=log_plot,
            width=width_slider,
            height=height_slider,
            x_limits=x_slider,
            save=save_button,
            latex=latex_button,
            dot_size=dot_slider,
            trajectory=traj_widget,
            out_param=out_param,
            in_params=in_params,
            plot_deviation=perturbed,
            font_scale=font_slider,
            perturbed=perturbed,
            precision=prec_widget,
        ),
    ).servable()

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
            log_plot,
            log_twin_plot,
            traj_widget,
        ),
        pn.Row(
            in_params,
            pn.Column(
                out_param,
                x_widget,
                y_widget,
                y2_widget,
                prec_widget,
            ),
        ),
        pn.Row(
            x_slider,
            dot_slider,
            title_widget,
        ),
        plot_pane,
    )
