try:
    from Deriv_dask import Deriv_dask
except:
    from scripts.Deriv_dask import Deriv_dask
try:
    from latexify import in_params_dic
except:
    from scripts.latexify import in_params_dic
import numpy as np
import xarray as xr
from glob import glob
import sys



def add_q_total(data):
    q_total = None
    for col in data.columns:
        if "Q" in col and not "IN" in col and not "OUT" in col:
            if q_total is not None:
                q_total += data[col]
            else:
                q_total = data[col].copy(deep=True)
    data["Q_total"] = q_total
    return data

def p_sat_icon(T):
    return 610.78 * np.exp(17.2693882*(T-273.15)/(T-35.86))

def p_sat_ice(T):
    return 610.78 * np.exp(21.8745584*(T-273.15)/(T-7.66))

def add_Si(data):
    data["Si"] = data["S"]*p_sat_icon(data["T"])/p_sat_ice(data["T"])
    return data

def add_p_sat(data):
    Rv = 8.3144598/0.018015265
    Ra = 8.3144598/0.02896546
    eps = Ra/Rv
    def p_sat(T):
        Tinv = 1.0/T
        logT = np.log(T)
        return np.exp( 54.842763 - 6763.22*Tinv - 4.21*logT + 0.000367*T + np.tanh(0.0415*(T-218.8))*(53.878 - 1331.22*Tinv - 9.44523*logT + 0.014025*T) )
    data["p_sat"] = p_sat(data["T"])
    data["p_sat_icon"] = p_sat_icon(data["T"])
    return data


# Helper function to load the data and plot a grid of output parameters
def plot_my_thing(path, netcdf_path, prefix, min_x=None, max_x=None, time_offset=0,
                  trajectories=None, f="*.nc_wcb",
                  x_axis="time_after_ascent", y_axis="pressure",
                  twin_axis=None, vertical_mark=None,
                  cross_mark_bool=False, store_path="pics/", in_params=None,
                  out_params=None,
                  plot_singles=True):

    deriv_4 = Deriv_dask(
        direc=path,
        parquet=False,
        netcdf=True,
        columns=None,
        backend="matplotlib",
        file_ending=f)

    if trajectories is not None:
        for traj in trajectories:
            print(f"Loading {traj}")
            deriv_4.cache_data(
                in_params=in_params,
                out_params=out_params,
                x_axis=x_axis,
                y_axis=y_axis,
                compute=True,
                trajectories=[traj],
                min_x=min_x,
                max_x=max_x)
            for col in out_params:
                if "OUT" in col:
                    deriv_4.cache[col] = np.abs(deriv_4.cache[col])
            # Add total amount of mixing ratios
            deriv_4.cache = add_q_total(deriv_4.cache)
            deriv_4.cache = add_Si(deriv_4.cache)

            cross_mark = None
            if cross_mark_bool:
                if min_x is not None:
                    start = min_x
                else:
                    start = deriv_4.cache[x_axis].min()
                    start = start - start%20
                if max_x is not None:
                    end = max_x
                else:
                    end = deriv_4.cache[x_axis].max()
                cross_mark = {x_axis: np.arange(start, end, 20)}

                ds = xr.open_dataset(netcdf_path, decode_times=False).to_dask_dataframe(dim_order=["ensemble", "trajectory", "time"])
                ds = ds.loc[ds["trajectory"] == traj]
                ds = ds.loc[ds[x_axis].isin(cross_mark[x_axis])].compute()
                ds["Output Parameter"] = out_params[0]
                ds["S"] /= 100
                ds = add_q_total(ds)
                ds = add_Si(ds)
                ds[x_axis] += 0.00001
                cross_mark[x_axis] += 0.00001

                deriv_4.cache = deriv_4.cache.loc[deriv_4.cache["Output Parameter"] == out_params[0]]
                deriv_4.cache = deriv_4.cache[~deriv_4.cache[x_axis].isin(cross_mark[x_axis])]
                deriv_4.cache = deriv_4.cache.append(ds, ignore_index=True)

            out_params.remove("S")
            grid = deriv_4.plot_grid_one_param(
                by="type",
                x_axis=x_axis,
                prefix=prefix,
                twin_axis=twin_axis,
                out_param=out_params[0],
                y_axis=out_params + ["Q_total", "Si", "S"],
                col_wrap=2, # Amount of plots per row
                rolling_avg=None,
                vertical_mark=vertical_mark,
                cross_mark=cross_mark,
                plot_path=store_path,
                alpha=1,
                plot_singles=plot_singles,
                formatter_limits=(-2,2), # Limits for x- and y-axis format
                s=20, # Size of scatter plots
                kind="scatter")
            out_params.append("S")
        return grid
    else:
        deriv_4.cache_data(
            in_params=in_params,
            out_params=out_params,
            x_axis=x_axis,
            y_axis=y_axis,
            compute=True,
            trajectories=trajectories,
            min_x=min_x,
            max_x=max_x)
        for col in out_params:
            if "OUT" in col:
                deriv_4.cache[col] = np.abs(deriv_4.cache[col])
        deriv_4.cache = add_q_total(deriv_4.cache)
        deriv_4.cache = add_Si(deriv_4.cache)

        cross_mark = None
        if cross_mark_bool:
            start = deriv_4.cache[x_axis].min()
            end = deriv_4.cache[x_axis].max()
            cross_mark = {x_axis: np.arange(start, end, 20)}

            ds = xr.open_dataset(netcdf_path, decode_times=False).to_dask_dataframe(dim_order=["ensemble", "trajectory", "time"])
            traj = deriv_4.cache["trajectory"].unique()
            ds = ds.loc[ds["trajectory"] == traj[0]]

            ds = ds.loc[ds[x_axis].isin(cross_mark[x_axis])].compute()

            ds["Output Parameter"] = out_params[0]
            ds["S"] /= 100
            ds = add_q_total(ds)
            ds = add_Si(ds)
            ds[x_axis] += 0.00001
            cross_mark[x_axis] += 0.00001

            deriv_4.cache = deriv_4.cache.loc[deriv_4.cache["Output Parameter"] == out_params[0]]
            deriv_4.cache = deriv_4.cache[~deriv_4.cache[x_axis].isin(cross_mark[x_axis])]
            deriv_4.cache = deriv_4.cache.append(ds, ignore_index=True)

        out_params.remove("S")
        grid = deriv_4.plot_grid_one_param(
            by="type",
            x_axis=x_axis,
            prefix=prefix,
            twin_axis=twin_axis,
            out_param=out_params[0],
            y_axis=out_params + ["Q_total", "Si", "S"],
            col_wrap=2,
            rolling_avg=None,
            vertical_mark=vertical_mark,
            cross_mark=cross_mark,
            plot_path=store_path,
            plot_singles=plot_singles,
            alpha=1,
            formatter_limits=(-2,2),
            s=20,
            kind="scatter")
        out_params.append("S")
        return grid

def plot_all(path, netcdf_path, store_path, suffixes, in_params, out_params,
    plot_singles=False, twin_axis=None, vertical_mark={"T": [273, 235]}, cross_mark=True):
    for suff in suffixes:
        files = sorted(glob(path + "/*"))
        for f in files:
            # print(f"Plotting for {f}")
            just_f = f.split("/")[-1]
            number = just_f.split("derivs_")[-1].split(".nc_wcb")[0]
            plot_my_thing(path,
                        netcdf_path,
                        suff + "_" + number + "_",
                        y_axis="pressure",
                        x_axis="time_after_ascent",
                        f=just_f,
                        twin_axis=twin_axis,
                        vertical_mark=vertical_mark,
            #               min_x=0, # You can use a start time for plotting
            #               max_x=100, # You can use an end time for plotting
                        cross_mark_bool=cross_mark,
                        store_path=store_path,
                        in_params=in_params,
                        plot_singles=plot_singles,
                        out_params=out_params)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = "../data/sim_processed/conv_400_0_t000000_p001_mult_outSat_sbShape_sbConv"
    if len(sys.argv) > 2:
        store_path = sys.argv[2]
    else:
        store_path = "../pics/"
    if len(sys.argv) > 3:
        plot_only = sys.argv[3]
        plot_only = plot_only[1::] # Remove first underscore
    else:
        plot_only = None
    if len(sys.argv) > 4:
        # The path to your input file of the simulation.
        netcdf_path = sys.argv[4]
    else:
        netcdf_path = "../data/conv_400_0_traj_t000000_p001.nc_wcb"

    plot_singles = False

    # Define which output parameters you would like to plot
    out_params = ["QI", "NCICE", "QC", "NCCLOUD", "QR", "NCRAIN",
                "QS", "NCSNOW", "QH", "NCHAIL", "QG", "NCGRAUPEL",
                "QV", "T", "pressure", "w", "z", "S",
                "QR_OUT", "QS_OUT", "QI_OUT", "QG_OUT", "QH_OUT"]

    # Define which derivatives you would like to load by adding their respective keys here
    keys = ["Misc"]
    # and either add all input parameters using this loop
    in_params = []
    for key in keys:
        in_params.extend(in_params_dic[key])
    # or define them directly
    in_params = ["da_1"]

    suffixes = ["outSat_sbConv", "inSat_sbConv", "outSat", "inSat",
                "outSat_sbShape_sbConv", "inSat_sbShape_sbConv", "outSat_sbShape", "inSat_sbShape"]
    if plot_only is not None:
        suffixes = [plot_only]

    # plot_all(path, netcdf_path, store_path, suffixes, in_params,
    #     twin_axis="pressure", cross_mark=True)
    # plot_all(path, netcdf_path, store_path, suffixes, in_params,
    #     twin_axis="pressure", cross_mark=False)
    plot_all(path, netcdf_path, store_path, suffixes, in_params, out_params,
        plot_singles=plot_singles,
        twin_axis=None, cross_mark=True)
    plot_all(path, netcdf_path, store_path, suffixes, in_params, out_params,
        plot_singles=plot_singles,
        twin_axis=None, cross_mark=False)
