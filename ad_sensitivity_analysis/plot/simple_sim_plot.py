"""Simple plots for checking simulation results.

"""
import matplotlib.pyplot as plt


def plot_distributions(ds, path, verbose):
    """
    Plot a histogram for every variable (column). The arrays are flattened.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    path : string
        Path to folder where plots will be stored as png.
    verbose : bool
        If true: print which variable is being plotted.
    """
    print("~*~*~*~Plotting distributions~*~*~*~")
    for col in ds:
        if verbose:
            print(f"Plotting {col}")
        ds[col].plot.hist(size=10, bins=100)
        plt.tight_layout()
        plt.savefig(f"{path}hist_{col}.png", dpi=300)
        plt.close()
        ds[col].plot.hist(size=10, bins=100, log=True)
        plt.tight_layout()
        plt.savefig(f"{path}hist_{col}_log.png", dpi=300)
        plt.close()
    print("done\n")


def plot_columns(ds, path, verbose):
    """
    Plot the columns, where ensemble and Output_Parameter_ID is squeezed
    (if present and possible).
    Otherwise an additional plot for each ensemble and ID

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with trajectory data.
    path : string
        Path where plots shall be stored.
    verbose : bool
        If true: print which variable is being plotted.
    """
    print("~*~*~*~Plotting columns~*~*~*~")
    ds_tmp = ds.copy()
    if "ensemble" in ds_tmp and len(ds_tmp["ensemble"]) == 1:
        ds_tmp = ds_tmp.squeeze(dim="ensemble", drop=True)
    if "Output_Parameter_ID" in ds_tmp and len(ds_tmp["Output_Parameter_ID"]) == 1:
        ds_tmp = ds_tmp.squeeze(dim="Output_Parameter_ID", drop=True)
    for col in ds_tmp:
        if len(ds_tmp[col]) == 0:
            continue
        if verbose:
            print(f"Plotting {col}")
        if "ensemble" in ds_tmp[col].dims:
            for ens in ds_tmp["ensemble"].values:
                if "Output_Parameter_ID" in ds_tmp[col].dims:
                    for i in ds_tmp["Output_Parameter_ID"].values:
                        ds_tmp[col].sel(
                            {"ensemble": ens, "Output_Parameter_ID": i}
                        ).plot(size=10)
                        plt.tight_layout()
                        plt.savefig(f"{path}{col}_ens{ens}_id{i}.png", dpi=300)
                        plt.close()
                else:
                    ds_tmp[col].sel({"ensemble": ens}).plot(size=10)
                    plt.tight_layout()
                    plt.savefig(f"{path}{col}_ens{ens}.png", dpi=300)
                    plt.close()
        else:
            ds_tmp[col].plot(size=10)
            plt.tight_layout()
            plt.savefig(f"{path}{col}.png", dpi=300)
            plt.close()
    print("done\n")
