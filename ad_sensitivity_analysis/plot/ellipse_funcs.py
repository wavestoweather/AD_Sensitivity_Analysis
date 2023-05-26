"""Helper functions to add ellipses to a plot.

"""
from holoviews import Ellipse as hvEllipse
from matplotlib.patches import Ellipse
from matplotlib import transforms
import numpy as np
from scipy.stats import chi2


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor="none", **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.
    By https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius = [np.sqrt(1 + pearson), np.sqrt(1 - pearson)]
    ellipse = Ellipse(
        (0, 0),
        width=ell_radius[0],
        height=ell_radius[1],
        facecolor=facecolor,
        **kwargs,
    )

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = (
        transforms.Affine2D()
        .rotate_deg(45)
        .scale(scale_x, scale_y)
        .translate(mean_x, mean_y)
    )

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def create_ellipse(x, y, chisq, **kwargs):
    """

    Parameters
    ----------
    x
    y
    chisq
    kwargs

    Returns
    -------

    """
    cov = np.cov(x, y)
    eigen, eigenv = np.linalg.eig(cov)
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    if eigen[0] > eigen[1]:
        bigger_idx = 0
    else:
        bigger_idx = 1
    rot = np.arctan2(eigenv[bigger_idx][1], eigenv[bigger_idx][0])
    if rot < 0:
        rot += np.pi
    major_scale = chisq * np.sqrt(eigen[bigger_idx])
    minor_scale = chisq * np.sqrt(eigen[(1 + bigger_idx) % 2])
    if bigger_idx == 1:
        scale = (major_scale, minor_scale)
    else:
        scale = (minor_scale, major_scale)
    # Major axis if needed
    # Seems a bit off at least on log log plots
    # slope = eigenv[1][0]/eigenv[0][0]
    # intercept = mean_y - slope*mean_x
    # max_x = np.max(x)
    # curve = hv.Curve(
    #     [(i, intercept + slope*i) for i in np.arange(mean_x, max_x, (max_x-mean_x)/10)]
    # ).opts(color="green", alpha=0.5)
    # return (hv.Ellipse(mean_x, mean_y, scale, orientation=-rot).opts(**kwargs)
    #     * curve
    #      * hv.Text(
    #         x.iloc[1], (intercept + slope*x.iloc[1])*1.1, f"Slope: {slope:.2f}" ) )
    return hvEllipse(mean_x, mean_y, scale, orientation=-rot).opts(**kwargs)


# pylint: disable=too-many-branches
def correl_conf_ell(
    df,
    x,
    y,
    linewidth,
    inf_val=None,
    backend="matplotlib",
    by_var=None,
    confidence=0.95,
    color=None,
    **kwargs
):
    """

    Parameters
    ----------
    df
    x
    y
    linewidth
    inf_val
    backend
    by_var
    confidence
    color
    kwargs

    Returns
    -------

    """
    new_kwargs = kwargs
    if color is not None:
        new_kwargs["color"] = color
        if linewidth is not None:
            if backend == "matplotlib":
                new_kwargs["linewidth"] = linewidth
            else:
                new_kwargs["line_width"] = linewidth
    if by_var is None:
        if inf_val is not None:
            df_tmp = df.loc[df[x] != inf_val]
        else:
            df_tmp = df
        x = df_tmp[x]
        y = df_tmp[y]
        if len(df_tmp[x]) == 1:
            return None
        return create_ellipse(x, y, chi2.ppf(confidence, 2), **new_kwargs)
    by_list = np.unique(df[by_var])
    ells = None
    for by_value in by_list:
        new_kwargs = kwargs
        if color is not None:
            new_kwargs["color"] = color[by_value]
        df_tmp = df.loc[df[by_var] == by_value]
        if inf_val is not None:
            df_tmp = df_tmp.loc[df_tmp[x] != inf_val]
        if len(df_tmp[x]) == 1:
            continue
        if ells is None:
            ells = create_ellipse(
                df_tmp[x], df_tmp[y], chi2.ppf(confidence, 2), **new_kwargs
            )
        else:
            ells = ells * create_ellipse(
                df_tmp[x], df_tmp[y], chi2.ppf(confidence, 2), **new_kwargs
            )
    return ells


def add_ellipse(
    tmp_df,
    mse_plot,
    sens_key,
    error_key,
    by_col,
    confidence,
    colors,
    kind,
    linewidth=2,
    inf_val=None,
):
    """

    Parameters
    ----------
    tmp_df
    mse_plot
    sens_key
    error_key
    by_col
    confidence
    colors
    kind
    linewidth
    inf_val

    Returns
    -------

    """
    if kind == "paper":
        tmp_df_ell = tmp_df.loc[tmp_df[sens_key] > inf_val]
        ells = correl_conf_ell(
            df=tmp_df_ell,
            x=sens_key,
            y=error_key,
            linewidth=linewidth,
            by=by_col,
            confidence=confidence,
            color=colors,
        )
        if ells is not None:
            mse_plot = mse_plot * ells
    else:
        ells = correl_conf_ell(
            df=tmp_df,
            x=sens_key,
            y=error_key,
            linewidth=linewidth,
            by=by_col,
            confidence=confidence,
            color=colors,
        )
        if ells is not None:
            mse_plot = mse_plot * ells
    return mse_plot
