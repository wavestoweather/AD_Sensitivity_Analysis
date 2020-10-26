from Deriv_dask import Deriv_dask
import numpy as np
import xarray as xr
from glob import glob
import sys

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

plot_singles = True
# A dictionary of all the derivatives where the key refers to the
# particle for which the parameter is defined. However,
# cloud parameters can affect rain droplets!
in_params_dic = {"Misc":
            ["da_1", "da_2", "de_1", "de_2", "dd", "dN_c", "dgamma",
            "dbeta_c", "dbeta_r", "ddelta1", "ddelta2", "dzeta",
            "drain_gfak", "dcloud_k_au", "dcloud_k_sc", "dkc_autocon",
            "dinv_z"],
            "Rain":
            ["drain_a_geo", "drain_b_geo", "drain_min_x", "drain_min_x_act",
            "drain_min_x_nuc_homo", "drain_min_x_nuc_hetero", "drain_min_x_melt",
            "drain_min_x_evap", "drain_min_x_freezing", "drain_min_x_depo",
            "drain_min_x_collision", "drain_min_x_collection",
            "drain_min_x_conversion", "drain_min_x_sedimentation",
            "drain_min_x_riming", "drain_max_x", "drain_sc_theta_q",
            "drain_sc_delta_q", "drain_sc_theta_n", "drain_sc_delta_n",
            "drain_s_vel", "drain_a_vel", "drain_b_vel", "drain_rho_v", "drain_c_z",
            "drain_sc_coll_n", "drain_cmu0", "drain_cmu1", "drain_cmu2", "drain_cmu3",
            "drain_cmu4", "drain_cmu5", "drain_alpha", "drain_beta", "drain_gamma",
            "drain_nu", "drain_g1", "drain_g2", "drain_mu", "drain_nm1", "drain_nm2",
            "drain_nm3", "drain_q_crit_c", "drain_d_crit_c", "drain_ecoll_c",
            "drain_cap", "drain_a_ven", "drain_b_ven", "drain_c_s", "drain_a_f",
            "drain_b_f", "drain_alfa_n", "drain_alfa_q", "drain_lambda",
            "drain_vsedi_min", "drain_vsedi_max"],
            "Cloud":
            ["dcloud_a_geo", "dcloud_b_geo", "dcloud_min_x", "dcloud_min_x_act",
            "dcloud_min_x_nuc_homo", "dcloud_min_x_nuc_hetero", "dcloud_min_x_melt",
            "dcloud_min_x_evap", "dcloud_min_x_freezing", "dcloud_min_x_depo",
            "dcloud_min_x_collision", "dcloud_min_x_collection",
            "dcloud_min_x_conversion", "dcloud_min_x_sedimentation",
            "dcloud_min_x_riming", "dcloud_max_x", "dcloud_sc_theta_q",
            "dcloud_sc_delta_q", "dcloud_sc_theta_n", "dcloud_sc_delta_n",
            "dcloud_s_vel", "dcloud_a_vel", "dcloud_b_vel", "dcloud_rho_v",
            "dcloud_c_z", "dcloud_sc_coll_n", "dcloud_cmu0", "dcloud_cmu1",
            "dcloud_cmu2", "dcloud_cmu3", "dcloud_cmu4", "dcloud_cmu5",
            "dcloud_alpha", "dcloud_beta", "dcloud_gamma", "dcloud_nu", "dcloud_g1",
            "dcloud_g2", "dcloud_mu", "dcloud_nm1", "dcloud_nm2", "dcloud_nm3",
            "dcloud_q_crit_c", "dcloud_d_crit_c", "dcloud_ecoll_c", "dcloud_cap",
            "dcloud_a_ven", "dcloud_b_ven", "dcloud_c_s", "dcloud_a_f", "dcloud_b_f",
            "dcloud_alfa_n", "dcloud_alfa_q", "dcloud_lambda", "dcloud_vsedi_min",
            "dcloud_vsedi_max"],
            "Graupel":
            ["dgraupel_a_geo", "dgraupel_b_geo", "dgraupel_min_x",
            "dgraupel_min_x_act", "dgraupel_min_x_nuc_homo",
            "dgraupel_min_x_nuc_hetero", "dgraupel_min_x_melt", "dgraupel_min_x_evap",
            "dgraupel_min_x_freezing", "dgraupel_min_x_depo",
            "dgraupel_min_x_collision", "dgraupel_min_x_collection",
            "dgraupel_min_x_conversion", "dgraupel_min_x_sedimentation",
            "dgraupel_min_x_riming", "dgraupel_max_x", "dgraupel_sc_theta_q",
            "dgraupel_sc_delta_q", "dgraupel_sc_theta_n", "dgraupel_sc_delta_n",
            "dgraupel_s_vel", "dgraupel_a_vel", "dgraupel_b_vel", "dgraupel_rho_v",
            "dgraupel_c_z", "dgraupel_sc_coll_n", "dgraupel_cmu0", "dgraupel_cmu1",
            "dgraupel_cmu2", "dgraupel_cmu3", "dgraupel_cmu4", "dgraupel_cmu5",
            "dgraupel_alpha", "dgraupel_beta", "dgraupel_gamma", "dgraupel_nu",
            "dgraupel_g1", "dgraupel_g2", "dgraupel_mu", "dgraupel_nm1",
            "dgraupel_nm2", "dgraupel_nm3", "dgraupel_q_crit_c", "dgraupel_d_crit_c",
            "dgraupel_ecoll_c", "dgraupel_cap", "dgraupel_a_ven", "dgraupel_b_ven",
            "dgraupel_c_s", "dgraupel_a_f", "dgraupel_b_f", "dgraupel_alfa_n",
            "dgraupel_alfa_q", "dgraupel_lambda", "dgraupel_vsedi_min",
            "dgraupel_vsedi_max"],
            "Hail":
            ["dhail_a_geo", "dhail_b_geo", "dhail_min_x", "dhail_min_x_act",
            "dhail_min_x_nuc_homo", "dhail_min_x_nuc_hetero", "dhail_min_x_melt",
            "dhail_min_x_evap", "dhail_min_x_freezing", "dhail_min_x_depo",
            "dhail_min_x_collision", "dhail_min_x_collection",
            "dhail_min_x_conversion", "dhail_min_x_sedimentation",
            "dhail_min_x_riming", "dhail_max_x", "dhail_sc_theta_q",
            "dhail_sc_delta_q", "dhail_sc_theta_n", "dhail_sc_delta_n", "dhail_s_vel",
            "dhail_a_vel", "dhail_b_vel", "dhail_rho_v", "dhail_c_z",
            "dhail_sc_coll_n", "dhail_cmu0", "dhail_cmu1", "dhail_cmu2", "dhail_cmu3",
            "dhail_cmu4", "dhail_cmu5", "dhail_alpha", "dhail_beta", "dhail_gamma",
            "dhail_nu", "dhail_g1", "dhail_g2", "dhail_mu", "dhail_nm1", "dhail_nm2",
            "dhail_nm3", "dhail_q_crit_c", "dhail_d_crit_c", "dhail_ecoll_c",
            "dhail_cap", "dhail_a_ven", "dhail_b_ven", "dhail_c_s", "dhail_a_f",
            "dhail_b_f", "dhail_alfa_n", "dhail_alfa_q", "dhail_lambda",
            "dhail_vsedi_min", "dhail_vsedi_max"],
            "Ice":
            ["dice_a_geo", "dice_b_geo", "dice_min_x", "dice_min_x_act",
            "dice_min_x_nuc_homo", "dice_min_x_nuc_hetero", "dice_min_x_melt",
            "dice_min_x_evap", "dice_min_x_freezing", "dice_min_x_depo",
            "dice_min_x_collision", "dice_min_x_collection", "dice_min_x_conversion",
            "dice_min_x_sedimentation", "dice_min_x_riming", "dice_max_x",
            "dice_sc_theta_q", "dice_sc_delta_q", "dice_sc_theta_n",
            "dice_sc_delta_n", "dice_s_vel", "dice_a_vel", "dice_b_vel", "dice_rho_v",
            "dice_c_z", "dice_sc_coll_n", "dice_cmu0", "dice_cmu1", "dice_cmu2",
            "dice_cmu3", "dice_cmu4", "dice_cmu5", "dice_alpha", "dice_beta",
            "dice_gamma", "dice_nu", "dice_g1", "dice_g2", "dice_mu", "dice_nm1",
            "dice_nm2", "dice_nm3", "dice_q_crit_c", "dice_d_crit_c", "dice_ecoll_c",
            "dice_cap", "dice_a_ven", "dice_b_ven", "dice_c_s", "dice_a_f",
            "dice_b_f", "dice_alfa_n", "dice_alfa_q", "dice_lambda", "dice_vsedi_min",
            "dice_vsedi_max"],
            "Snow":
            ["dsnow_a_geo", "dsnow_b_geo", "dsnow_min_x", "dsnow_min_x_act",
            "dsnow_min_x_nuc_homo", "dsnow_min_x_nuc_hetero", "dsnow_min_x_melt",
            "dsnow_min_x_evap", "dsnow_min_x_freezing", "dsnow_min_x_depo",
            "dsnow_min_x_collision", "dsnow_min_x_collection",
            "dsnow_min_x_conversion", "dsnow_min_x_sedimentation",
            "dsnow_min_x_riming", "dsnow_max_x", "dsnow_sc_theta_q",
            "dsnow_sc_delta_q", "dsnow_sc_theta_n", "dsnow_sc_delta_n", "dsnow_s_vel",
            "dsnow_a_vel", "dsnow_b_vel", "dsnow_rho_v", "dsnow_c_z",
            "dsnow_sc_coll_n", "dsnow_cmu0", "dsnow_cmu1", "dsnow_cmu2", "dsnow_cmu3",
            "dsnow_cmu4", "dsnow_cmu5", "dsnow_alpha", "dsnow_beta", "dsnow_gamma",
            "dsnow_nu", "dsnow_g1", "dsnow_g2", "dsnow_mu", "dsnow_nm1", "dsnow_nm2",
            "dsnow_nm3", "dsnow_q_crit_c", "dsnow_d_crit_c", "dsnow_ecoll_c",
            "dsnow_cap", "dsnow_a_ven", "dsnow_b_ven", "dsnow_c_s", "dsnow_a_f",
            "dsnow_b_f", "dsnow_alfa_n", "dsnow_alfa_q", "dsnow_lambda",
            "dsnow_vsedi_min", "dsnow_vsedi_max"]}

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
                  x_axis="time_after_ascent", y_axis="pressure", traj="Convective 400hPa",
                  twin_axis=None, vertical_mark=None,
                  cross_mark_bool=False):

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

# Plot with "pressure" as twin axis
# Plot with crosses for the input simulation data
# (overrides according datapoints of the output simulation)
for suff in suffixes:
    files = sorted(glob(path + "/*"))
    for f in files:
        just_f = f.split("/")[-1]
        number = just_f.split("derivs_")[-1].split(".nc_wcb")[0]
        plot_my_thing(path,
                    netcdf_path,
                    suff + "_" + number + "_",
                    y_axis="pressure",
                    x_axis="time_after_ascent",
                    traj="Convective 400hPa",
                    f=just_f,
                    twin_axis="pressure",
                    vertical_mark={"T": [273, 235]},
        #               min_x=0, # You can use a start time for plotting
        #               max_x=100, # You can use an end time for plotting
                    cross_mark_bool=True)

# Plot with "pressure" as twin axis
# Plot no crosses
# for suff in suffixes:
#     files = sorted(glob(path + "/*"))
#     for f in files:
#         just_f = f.split("/")[-1]
#         number = just_f.split("derivs_")[-1].split(".nc_wcb")[0]
#         plot_my_thing(path,
#                     netcdf_path,
#                     suff + "_noCross_" + number + "_",
#                     y_axis="pressure",
#                     x_axis="time_after_ascent",
#                     traj="Convective 400hPa",
#                     f=just_f,
#                     twin_axis="pressure",
#                     vertical_mark={"T": [273, 235]},
#         #               min_x=0,
#         #               max_x=100,
#                     cross_mark_bool=False)

# No twin axis
# Plot with crosses for the input simulation data
# (overrides according datapoints of the output simulation)
for suff in suffixes:
    files = sorted(glob(path + "/*"))
    for f in files:
        just_f = f.split("/")[-1]
        number = just_f.split("derivs_")[-1].split(".nc_wcb")[0]
        plot_my_thing(path,
                    netcdf_path,
                    suff + "_no_twin_" + number + "_",
                    y_axis="pressure",
                    x_axis="time_after_ascent",
                    traj="Convective 400hPa",
                    f=just_f,
                    twin_axis=None,
                    vertical_mark={"T": [273, 235]},#, "S": [0.97]},
        #               min_x=0,
        #               max_x=100,
                    cross_mark_bool=True)

# No twin axis
# No crosses
# for suff in suffixes:
#     files = sorted(glob(path + "/*"))
#     for f in files:
#         just_f = f.split("/")[-1]
#         number = just_f.split("derivs_")[-1].split(".nc_wcb")[0]
#         plot_my_thing(path,
#                     netcdf_path,
#                     suff + "_noCross_no_twin_" + number + "_",
#                     y_axis="pressure",
#                     x_axis="time_after_ascent",
#                     traj="Convective 400hPa",
#                     f=just_f,
#                     twin_axis=None,
#                     vertical_mark={"T": [273, 235]},#, "S": [0.97]},
#         #               min_x=0,
#         #               max_x=100,
#                     cross_mark_bool=False)

