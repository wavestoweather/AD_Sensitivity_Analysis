"""Create panel widgets for interactive plotting.

"""
import panel as pn


def create_image_size_slider(width=True, pixel=True):
    """

    Parameters
    ----------
    width
    pixel

    Returns
    -------

    """
    if width:
        name = "Width"
    else:
        name = "Height"
    if pixel:
        unit = "pixel"
    else:
        unit = "inches"
    return pn.widgets.IntSlider(
        name=f"{name} in {unit}",
        start=3 - 3 * pixel + 300 * pixel,
        end=15 - 15 * pixel + 3000 * pixel,
        step=1 - pixel + 150 * pixel,
        value=9 - 9 * pixel + 800 * pixel,
    )


def create_float_slider(ds, col):
    """
    Create a float slider based on a column of ds.

    Parameters
    ----------
    ds : xarra.Dataset
        A dataset with at least the column given by 'col'.
    col : string
        Name of the column.

    Returns
    -------
    panel.widgets.FloatSlider with range [min(ds[col]), max(ds[col])] and
    step size is the size between the first and second value.
    """
    return pn.widgets.FloatSlider(
        name=col,
        start=ds[col].min().values.item(),
        end=ds[col].max().values.item(),
        step=(ds[col].values[1] - ds[col].values[0]),
        value=ds[col].values[-4],
    )


def create_time_slider_select(ds):
    """
    Create a slider or a select widget for the index of the time dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with (dimension) 'time'.

    Returns
    -------
        panel.widgets.IntSlider if more than one time value is available.
        panel.widgets.Select otherwise.
    """
    if len(ds["time"]) > 1:
        return pn.widgets.IntSlider(
            name="timestep",
            start=0,
            end=len(ds["time"]) - 1,
            step=1,
            value=0,
        )
    return pn.widgets.Select(
        name="timestep",
        value=0,
        options=[0],
    )


def create_title_widget():
    """

    Returns
    -------

    """
    return pn.widgets.TextInput(
        name="Title",
        placeholder="",
    )


def image_font_size_row(pixel=True):
    """
    Creates a row with widgets to adjust width, height, and font size of a plot.

    Parameters
    ----------
    unit : bool
        If true, create width and height widgets based in pixels. Otherwise, the
        size is given in inches.

    Returns
    -------
    panel.Row with a width_slider, height_slider, and a font_slider in this order.
    """
    width_slider = create_image_size_slider(width=True, pixel=pixel)
    height_slider = create_image_size_slider(width=False, pixel=pixel)
    font_slider = pn.widgets.FloatSlider(
        name="Scale fontsize",
        start=0.2,
        end=5,
        step=0.1,
        value=0.7,
    )
    return pn.Row(
        width_slider,
        height_slider,
        font_slider,
    )


def create_save_latex_row():
    """
    Create a panel.Row with three widgets to save a plot and use Latex to plot the
    text.

    Returns
    -------
    panel.Row with a widget to specify a path to save the image, a button to save it,
    and a toggle to use Latex in this order.
    """
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
    return pn.Row(
        save_to_field,
        save_button,
        latex_button,
    )


def create_in_params_selector(ds):
    """

    Parameters
    ----------
    ds

    Returns
    -------

    """
    in_param_values = []
    for col in ds:
        if "rank" in col:
            in_param_values.append(col)
    in_param_values.sort()
    return pn.widgets.CrossSelector(
        name="Parameter",
        value=in_param_values[0:2],
        options=in_param_values,
    )


def create_kind_in_params_selector(ds):
    """

    Parameters
    ----------
    ds

    Returns
    -------

    """
    in_params = []
    kind_params = []
    for col in ds:
        if "Mean " in col:
            in_params.append(col[5:])
            if "Mean" not in kind_params:
                kind_params.append("Mean")
        if "Var " in col or "Std " in col:
            if "Std" not in kind_params:
                kind_params.append("Std")
            if "Var" not in kind_params:
                kind_params.append("Var")
        if "Min " in col:
            if "Min" not in kind_params:
                kind_params.append("Min")
        if "Max " in col:
            if "Max" not in kind_params:
                kind_params.append("Max")
    if "counts" in ds:
        in_params.append("counts")
    if "Top_Parameter" in ds:
        in_params.append("Top_Parameter")
    in_params.sort()
    in_param = pn.widgets.Select(
        name="Color according to",
        value=in_params[-1],
        options=in_params,
    )
    kind_param = pn.widgets.RadioButtonGroup(
        name="Kind",
        value=kind_params[0],
        options=kind_params,
        button_type="primary",
    )
    return kind_param, in_param


def create_out_params_buttons(ds):
    """

    Parameters
    ----------
    ds
    coord

    Returns
    -------

    """
    if "Output Parameter" in ds:
        coord = "Output Parameter"
    else:
        coord = "Output_Parameter"
    return pn.widgets.RadioButtonGroup(
        name="Output Parameter",
        value=ds[coord].values[0],
        options=list(ds[coord].values),
        button_type="primary",
    )


def create_line_slider():
    """

    Returns
    -------

    """
    return pn.widgets.FloatSlider(
        name="Change the line width",
        start=1,
        end=10,
        step=0.5,
        value=2,
    )


def create_flow_phase_ignore_row():
    """

    Returns
    -------

    """
    flow_button = pn.widgets.Toggle(
        name="Show flow",
        button_type="success",
    )
    phase_button = pn.widgets.Toggle(
        name="Show phase",
        button_type="success",
    )
    ignore_zero_button = pn.widgets.Toggle(
        name="Ignore zero gradients",
        button_type="success",
    )
    return pn.Row(
        flow_button,
        phase_button,
        ignore_zero_button,
    )


def create_out_params_select(ds):
    """

    Parameters
    ----------
    ds

    Returns
    -------

    """
    if "Output Parameter" in ds:
        coord = "Output Parameter"
    else:
        coord = "Output_Parameter"
    return pn.widgets.Select(
        name="Output Parameter",
        value=ds[coord].values[0],
        options=ds[coord].values.tolist(),
    )


def create_color_select():
    """

    Returns
    -------

    """
    return pn.widgets.Select(
        name="Colormap",
        value="RdBu",
        options=[
            "RdBu",
            "viridis",
            "Blues",
            "Reds",
            "PiYG",
            "PRGn",
            "BrBG",
            "PuOr",
            "plasma",
            "cividis",
        ],
    )
