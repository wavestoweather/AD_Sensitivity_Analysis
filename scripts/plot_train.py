import numpy as np
import xarray as xr
import hvplot.xarray  # noqa
from tqdm import tqdm


def load_data(f_name):
    return xr.open_dataset(f_name, decode_times=False, engine="netcdf4")


def get_info(ds):
    nparams = len(ds["perturbed"])
    n = len(ds["ensemble"]) - 1
    ny = len(ds["Output_Parameter"])
    print(f"# of parameters np: {nparams}")
    print(f"# of perturbations n: {n}")
    print(f"# of output parameters ny: {ny}")


def get_training_data(ds):
    n = len(ds["ensemble"]) - 1
    ny = len(ds["Output_Parameter"])
    n_time = len(ds["time"])
    n_traj = len(ds["trajectory"])

    impact_matrix = np.zeros((n_traj, n_time, n, ny))
    prediction_matrix = np.zeros((n_traj, n_time, n, ny))

    for ens in tqdm(ds["ensemble"].values):
        if ens == 0:
            continue
        for p_i, out_p in enumerate(ds["Output_Parameter"].values):
            ens_i = int(ens - 1)
            impact_matrix[:, :, ens_i, p_i] = ds[out_p].sel({"ensemble": 0}) - ds[
                out_p
            ].sel({"ensemble": ens})

    for ens in tqdm(ds["ensemble"].values):
        if ens == 0:
            continue
        ens_i = int(ens - 1)
        for pert in ds["perturbed"].values:
            orig_val = ds["perturbation_value"].sel({"ensemble": 0, "perturbed": pert})
            ds_pert = (
                ds["perturbation_value"].sel({"ensemble": ens, "perturbed": pert})
                - orig_val
            )

            for p_i, out_p in enumerate(ds["Output_Parameter"].values):
                for traj_i, traj in enumerate(ds["trajectory"]):
                    prediction_matrix[traj_i, :, ens_i, p_i] += (
                        ds["d" + pert].sel(
                            {"Output_Parameter": out_p, "trajectory": traj}
                        )
                        * ds_pert
                    )

    ds_training = xr.Dataset(
        coords=dict(
            trajectory=ds["trajectory"],
            time=ds["time"],
            ensemble=np.arange(n),
            Output_Parameter=ds["Output_Parameter"].values,
        ),
        data_vars=dict(
            impact_matrix=(
                ["trajectory", "time", "ensemble", "Output_Parameter"],
                impact_matrix,
            ),
            prediction_matrix=(
                ["trajectory", "time", "ensemble", "Output_Parameter"],
                prediction_matrix,
            ),
        ),
        attrs=dict(description="Training data for nets based on Vladiana."),
    )
    return ds_training


def plot_heatmap(ds, plot="both", width=1600, height=600):
    """

    Parameters
    ----------
    ds : xarray.Dataset
        The training data
    plot: string
        'prediction': Plot the prediction matrix.
        'impact': Plot the impact matrix.
        'both': Plot both.
    width
    height

    Returns
    -------

    """
    if plot == "impact":
        c = "impact_matrix"
        ylabel = "Output Parameter Deviation"
    elif plot == "prediction":
        c = "prediction_matrix"
        ylabel = "Predicted Deviation"
    else:
        imp_plot = np.abs(ds).hvplot.heatmap(
            x="ensemble",
            y="Output_Parameter",
            C="impact_matrix",
            height=height,
            width=width,
            ylabel="Output Parameter Deviation",
            logz=True,
            clim=(1e-14, 1e10),
        )
        pred_plot = np.abs(ds).hvplot.heatmap(
            x="ensemble",
            y="Output_Parameter",
            C="prediction_matrix",
            height=height,
            width=width,
            ylabel="Predicted Deviation",
            logz=True,
            clim=(1e-14, 1e10),
        )
        return (imp_plot + pred_plot).cols(1)
    return np.abs(ds).hvplot.heatmap(
        x="ensemble",
        y="Output_Parameter",
        C=c,
        height=height,
        width=width,
        ylabel=ylabel,
        logz=True,
        clim=(1e-14, 1e10),
    )


def plot_training_time(
    ds,
    min_x=None,
    max_x=None,
    width=1600,
    height=1200,
    logx=False,
    logy=False,
    logy2=False,
):
    """

    Parameters
    ----------
    ds
    min_x
    max_x
    width
    height
    logx
    logy
    logy2

    Returns
    -------

    """
    if min_x is not None:
        ds = ds.where(ds["time"] >= min_x, drop=True)
    if max_x is not None:
        ds = ds.where(ds["time"] <= max_x, drop=True)
    if logy:
        ds["impact_matrix"] = np.abs(ds["impact_matrix"])
    if logy2:
        ds["prediction_matrix"] = np.abs(ds["prediction_matrix"])
        plots = ds.hvplot(
            x="time",
            y="prediction_matrix",
            by="ensemble",
            ylabel="Predicted Deviation",
            logx=logx,
            logy=logy2,
            width=width,
            height=height,
            ylim=(1e-14, 1e3),
        )
    else:
        plots = ds.hvplot(
            x="time",
            y="prediction_matrix",
            by="ensemble",
            ylabel="Predicted Deviation",
            logx=logx,
            logy=logy2,
            width=width,
            height=height,
        )

    return (
        ds.hvplot(
            x="time",
            y="impact_matrix",
            ylabel="True Deviation",
            logx=logx,
            logy=logy,
            width=width,
            height=height,
            by="ensemble",
        )
        + plots
    ).cols(1)


def plot_time(
    ds,
    state_var,
    model_par=None,
    min_x=None,
    max_x=None,
    width=1600,
    height=1200,
    logx=False,
    logy=False,
    logy2=False,
):
    """
    Plot variables over time.

    Parameters
    ----------
    ds : xarray.Dataset
        The output of the simulation
    state_var
    model_par
    min_x
    max_x
    width
    height
    logx
    logy
    logy2

    Returns
    -------

    """
    if min_x is not None:
        ds = ds.where(ds["time"] >= min_x, drop=True)
    if max_x is not None:
        ds = ds.where(ds["time"] <= max_x, drop=True)
    if model_par is None:
        model_par = ds["perturbed"].values
    ylim = None
    if isinstance(model_par, list) or isinstance(model_par, np.ndarray):
        y2_label = "Gradients"
        if model_par[0][0] != "d":
            for i in range(len(model_par)):
                model_par[i] = "d" + model_par[i]
        plots = None
        for p in model_par:
            if logy2:
                ds[p] = np.abs(ds[p])
                y2_label = "log10 |" + y2_label + "|"
                ylim = (1e-14, 1e3)
            if plots is not None:
                plots *= ds.drop_dims(["perturbed", "ensemble"]).hvplot(
                    x="time",
                    y=p,
                    ylabel=y2_label,
                    logx=logx,
                    logy=logy2,
                    width=width,
                    height=height,
                    label=p,
                    ylim=ylim,
                )
            else:
                plots = ds.drop_dims(["perturbed", "ensemble"]).hvplot(
                    x="time",
                    y=p,
                    ylabel=y2_label,
                    logx=logx,
                    logy=logy2,
                    width=width,
                    height=height,
                    label=p,
                    ylim=ylim,
                )
    else:
        pert_value = ds["perturbation_value"].sel({"perturbed": model_par})
        if "d" != model_par[0]:
            model_par = "d" + model_par
        y2_label = model_par
        ds[model_par] *= pert_value
        if logy2:
            ds[model_par] = np.abs(ds[model_par])
            y2_label = "log10 |" + y2_label + "|"
            ylim = (1e-14, 1e3)
        plots = ds.drop_dims(["perturbed", "ensemble"]).hvplot(
            x="time",
            y=model_par,
            ylabel=y2_label,
            logx=logx,
            logy=logy2,
            width=width,
            height=height,
            ylim=ylim,
        )

    return (
        ds.drop_dims("perturbed").hvplot(
            x="time",
            y=state_var,
            logx=logx,
            logy=logy,
            width=width,
            height=height,
            by="ensemble",
        )
        + plots
    ).cols(1)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        """
        Plot the training data either as scatter plot or histogram.
        Plot results from trained networks.
        """,
    )
    ds = load_data("/media/mahieron/Austausch/AD_Vis/train_data_test.nc")
    get_info(ds)
