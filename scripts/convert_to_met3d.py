# Convert WCB2 netcdf files and outputs of the simulation to comply with Met.3D
# See https://met3d.readthedocs.io/en/latest/dataformats.html#trajectory-data-in-netcdf-format
import os
import pandas as pd
import xarray as xr
import numpy as np
from timeit import default_timer as timer
import progressbar as pb
# import sys

# np.set_printoptions(threshold=sys.maxsize)

met_dims = ["ensemble", "trajectory", "time", "start_lon",
            "start_lat", "start_isobaric"]

path = "/data/project/wcb/netcdf/vladiana"
# path = "/lustre/project/m2_jgu-tapt/cosmo_output/vladiana/traj"
store_path = "/data/project/wcb/netcdf/vladiana_met/"
# /lustre/project/m2_zdvresearch/mahieron/netcdf_vladiana

# 400 hPa and 600 hPa ascent
window_conv_400 = 1 * 3 * 60
window_conv_600 = 3 * 3 * 60
window_slan_400 = 35 * 6 * 3
window_slan_400_min = 15 * 6 * 3
window_slan_600 = 22  * 3 * 60
window_slan_600_min = 65 * 6 * 3

def find_runs(x):
    """Find runs of consecutive items in an array."""

    n = x.shape[0]

    # handle empty array
    if n == 0:
        return np.array([]), np.array([]), np.array([])

    else:
        # find run starts
        loc_run_start = np.empty(n, dtype=bool)
        loc_run_start[0] = True
        np.not_equal(x[:-1], x[1:], out=loc_run_start[1:])
        run_starts = np.nonzero(loc_run_start)[0]

        # find run values
        run_values = x[loc_run_start]

        # find run lengths
        run_lengths = np.diff(np.append(run_starts, n))

        return run_values, run_starts.astype(np.int64), run_lengths.astype(np.int64)

def differ(x, axis, hPa, debug=False):
    if debug:
        print("x")
        print(np.shape(x))
    window_size = len(x[0][0][0])
    ascent = np.argmax(x, axis=axis) < np.argmin(x, axis=3)
    amount = np.max(x, axis=axis) - np.min(x, axis=axis) >= hPa*100
    # ascent = np.nanargmax(x, axis=axis) < np.nanargmin(x, axis=3)
    # amount = np.nanmax(x, axis=axis) - np.nanmin(x, axis=axis) >= hPa*100
    both = np.logical_and(ascent, amount)

    # Calculate the differences within every window
    differences = np.diff(x, axis=3) # Ignore nan
    if debug:
        print("diffs")
        print(np.shape(differences))
        print("ascent")
        print(np.shape(ascent))
        print("amount")
        print(np.shape(amount))
        print("both")
        print(np.shape(both))
    # Get minimum length of window with value >= hPa
    min_lengths = np.full(np.shape(both), np.inf)
    counter = 0

    for ens in range(len(differences)):
        if not both[ens].any():
            continue
        for traj in range(len(differences[ens])):
            if not both[ens][traj].any():
                # No timestep in this trajectory satisfies
                # the constraint
                continue
            for timestep in range(len(differences[ens][traj])):
                if not both[ens][traj][timestep].any():
                    continue
                window = differences[ens][traj][timestep]
                start = 0
                end = 0
                curr_sum = 0
                min_len = np.inf
                    # 17926
                while(end < len(window)):
                    # Find a window where the ascend is done by pushing the end further
                    # Reset the start if the end is suddenly a strong descend
                    while( (curr_sum > -hPa*100 and end < len(window)) or (np.isnan(curr_sum)) ):
                        if (curr_sum >= 0 and window[end] < 0 or np.isnan(curr_sum) ):
                            start = end;
                            curr_sum = 0;

                        curr_sum += window[end]
                        end += 1
                    # Check, if a smaller window exists where the ascend is done by pushing the start
                    while(curr_sum <= -hPa*100 and start < len(window)):
                        if(end-start < min_len):
                            min_len = end-start
                        curr_sum -= window[start]
                        start += 1
                min_lengths[ens][traj][timestep] = min_len

    # Take the minimum overall and that's whenever True shall stand
    # Those are minimum window sizes for every trajectory
    return_bools = []
    min_len_traj_all = []
    for ens in range(len(x)):
        min_len_traj = np.nanmin(min_lengths[ens], axis=1)
        min_len_traj_all.append(min_len_traj)
        min_len_traj[min_len_traj == np.inf] = -1
        return_bools_tmp = np.full(np.shape(both[ens]), 0)#, dtype=bool)
        return_bools_tmp = np.transpose(return_bools_tmp)
        min_lengths_trans = min_lengths[ens].transpose()
        if debug:
            print("min_lengths_trans")
            print(np.shape(min_lengths_trans))
            print("min_len_traj")
            print(np.shape(min_len_traj))
            print("return_bools_tmp")
            print(np.shape(return_bools_tmp))
        for timestep in range(len(return_bools_tmp)):
            return_bools_tmp[timestep] = (min_lengths_trans[timestep] == min_len_traj)
        return_bools.append(np.transpose(return_bools_tmp))

    # Shift everything such that the beginning starts at the actual start and ends accordingly
    for ens in range(len(return_bools)):
        for traj in range(len(return_bools[ens])):
            if min_len_traj_all[ens][traj] == -1:
                continue
            min_len_traj = min_len_traj_all[ens]
            vals, start, length = find_runs(return_bools[ens][traj])

            for i in range(len(vals)):
                if vals[i] > 0:
                    set_start = int(start[i] - min_len_traj[traj])
                    set_end = set_start + min_len_traj[traj] + 1

                    if length[i] > min_len_traj[traj]:
                        set_end = set_start + length[i] + 1

                    return_bools[ens][traj][start[i]:length[i]+start[i]] = False
                    return_bools[ens][traj][set_start:int(set_end)] = True

    return return_bools

def differ_slan(x, axis, hPa, min_window):
    window_size = len(x[0][0][0])
    ascent = np.argmax(x, axis=axis) < np.argmin(x, axis=3)
    amount = np.max(x, axis=axis) - np.min(x, axis=axis) >= hPa*100
    both = np.logical_and(ascent, amount)

    # Calculate the differences within every window
    differences = np.diff(x, axis=3)
    # Get minimum length of window with value >= hPa
    min_lengths = np.full(np.shape(both), np.inf)


    for ens in range(len(differences)):
        if not both[ens].any():
            continue
        for traj in range(len(differences[ens])):
            if not both[ens][traj].any():
                # No timestep in this trajectory satisfies
                # the constraint
                continue
            for timestep in range(len(differences[ens][traj])):
                if not both[ens][traj][timestep].any():
                    continue
                window = differences[ens][traj][timestep]
                start = 0
                end = 0
                curr_sum = 0
                min_len = np.inf
                while(end < len(window)):
                    while( (curr_sum > -hPa*100 and end < len(window)) or (np.isnan(curr_sum)) ):
                        if (curr_sum >= 0 and window[end] < 0 or np.isnan(curr_sum) ):
                            start = end
                            curr_sum = 0
                        curr_sum += window[end]
                        end += 1
                    while(curr_sum <= -hPa*100 and start < len(window)):
                        if(end-start < min_len and end-start >= min_window):
                            min_len = end-start
                        curr_sum -= window[start]
                        start += 1
                min_lengths[ens][traj][timestep] = min_len
    # Take the minimum overall and that's whenever True shall stand
    # Those are minimum window sizes for every trajectory
    return_bools = []
    min_len_traj_all = []
    for ens in range(len(x)):
        min_len_traj = np.nanmin(min_lengths[ens], axis=1)
        min_len_traj_all.append(min_len_traj)
        min_len_traj[min_len_traj == np.inf] = -1
        return_bools_tmp = np.full(np.shape(both[ens]), 0)#, dtype=bool)
        return_bools_tmp = np.transpose(return_bools_tmp)
        min_lengths_trans = min_lengths[ens].transpose()
        for timestep in range(len(return_bools_tmp)):
            return_bools_tmp[timestep] = (min_lengths_trans[timestep] == min_len_traj)
        return_bools.append(np.transpose(return_bools_tmp))

    # Shift everything such that the beginning starts at the actual start and ends accordingly
    for ens in range(len(return_bools)):
        for traj in range(len(return_bools[ens])):
            if min_len_traj_all[ens][traj] == -1:
                continue
            vals, start, length = find_runs(return_bools[ens][traj])
            for i in range(len(vals)):
                if vals[i] > 0:
                    set_start = int(start[i] - min_len_traj[traj])
                    set_end = set_start + min_len_traj[traj] + 1

                    if length[i] > min_len_traj[traj]:
                        set_end = set_start + length[i] + 1

                    return_bools[ens][traj][start[i]:length[i]+start[i]] = False
                    return_bools[ens][traj][set_start:int(set_end)] = True
    return return_bools

def add_norm_time(df, norm_col, group, columns=None, flag=None):
    '''
    Return a view that consists only of entries that are flagged.
    Those are normed along norm_col such that every entry for every
    trajectory starts at norm_col==0. columns is a list of
    columns that the returned view shall have.
    if columns is None, take all columns. If flag is None, take all trajectories.
    '''
    if columns is None:
        df_flagged = df.copy()
    else:
        df_flagged = df[columns + [flag] + [norm_col] + [group]]
    def reducer(x, col):
        mini = x.loc[x[flag] == True][col].min()
        x[col] = x[col] - mini
        return x
    normed = df_flagged.groupby([group]).apply(reducer, norm_col)
    df["time_after_ascent"] = normed["time"]
    return df


def add_sat(ds):
    def p_sat(T):
        return 610.78 * np.exp(17.2693882 * (T-273.16)/(T-35.86))
    def convert_qv_to_S(p, T, qv):
        eps = (8.3144598/0.02896546)/(8.3144598/0.018015265)
        return (p*qv)/((eps+qv)*p_sat(T))*100
    ds["S"] = convert_qv_to_S(ds["pressure"], ds["T"], ds["QV"])
    return ds

def add_ascend_velocity(ds):
    ds["w"] = ds["z"].diff(dim="time", label="lower")/20
    return ds

def add_turb_flux(ds):
    q_total = ds["QV"] + ds["QC"] + ds["QR"] + ds["QS"] + ds["QI"] + ds["QG"]
    q_in_total = ds["QR_IN"] + ds["QS_IN"] + ds["QI_IN"] + ds["QG_IN"]
    q_out_total = ds["QR_OUT"] + ds["QS_OUT"] + ds["QI_OUT"] + ds["QG_OUT"]
    ds["Q_TURBULENCE"] = q_total.diff(dim="time", label="upper")
    ds["Q_TURBULENCE"] = ds["Q_TURBULENCE"] - q_in_total + q_out_total
    return ds

def convert_wcb2(f, store_path, fl, ensemble):
    ds = xr.open_dataset(f).to_dataframe().reset_index()

    ds.rename(columns={
        "id": "trajectory",
        "longitude": "lon",
        "latitude": "lat",
        "P": "pressure"}, inplace=True)
    ds["ensemble"] = ensemble # Number of ensembles
    ds_2 = xr.open_dataset(f[:-4])

    duration = ds_2.attrs["duration_in_sec"]
    pollon = ds_2.attrs["pollon"]
    pollat = ds_2.attrs["pollat"]
    output_timestep_in_sec = ds_2.attrs["output_timestep_in_sec"]
    ref_year = ds_2.attrs["ref_year"]
    ref_month = ds_2.attrs["ref_month"]
    ref_day = ds_2.attrs["ref_day"]
    ref_hour = ds_2.attrs["ref_hour"]
    ref_min = ds_2.attrs["ref_min"]
    ref_sec = ds_2.attrs["ref_sec"]
    ref_time = f"{ref_year}-{ref_month:02}-{ref_day:02} {ref_hour:02}:{ref_min:02}:{ref_sec:02}"

    ds = ds[ds.lon != -999]
    ds = ds[ds.lat != -999]
    ds = ds[ds.z != -999]
    ds = ds.dropna(how="all")

    ds = xr.Dataset.from_dataframe(ds.set_index(["ensemble", "trajectory", "time"])) # ,
    # Set flag
    for flag in ["conv_600", "conv_400", "slan_400", "slan_600"]:
            ds[flag] = False

    for col in ds:
        ds[col] = ds[col].astype("float64", casting="safe", copy=False)

    if fl == "conv_600":
        t_c = timer()
        conv_600 = ds["pressure"].rolling(dim={"time": window_conv_600}, min_periods=1).reduce(
            differ, **{"hPa": 600}).fillna(False).astype(dtype=bool)
        ds = ds.assign(conv_600=conv_600)
        t_c2 = timer()
        print("Got conv_600 in {} s".format(t_c2-t_c), flush=True)
    elif fl == "conv_400":
        t_c = timer()
        conv_400 = ds["pressure"].rolling(dim={"time": window_conv_400}, min_periods=1).reduce(
            differ, **{"hPa": 400, "debug": False}).fillna(False).astype(dtype=bool)
        ds = ds.assign(conv_400=conv_400)
        t_c2 = timer()
        print("Got conv_400 in {} s".format(t_c2-t_c), flush=True)
    elif fl == "slan_600":
        t_c = timer()
        slan_600 = ds["pressure"].rolling(dim={"time": window_slan_600}, min_periods=1).reduce(
            differ_slan, **{"hPa": 600, "min_window": window_slan_600_min}).fillna(False).astype(dtype=bool)
        ds = ds.assign(slan_600=slan_600)
        t_c2 = timer()
        print("Got slan_600 in {} s".format(t_c2-t_c), flush=True)
    elif fl == "slan_400":
        t_c = timer()
        slan_400 = ds["pressure"].rolling(dim={"time": window_slan_400}, min_periods=1).reduce(
            differ_slan, **{"hPa": 400, "min_window": window_slan_400_min}).fillna(False).astype(dtype=bool)
        ds = ds.assign(slan_400=slan_400)
        t_c2 = timer()
        print("Got slan_400 in {} s".format(t_c2-t_c), flush=True)

    # If there is none, return
    if not ds[fl].any():
        print(f"Got no {fl}")
        return
    ds  = ds.to_dataframe().reset_index()
    ids = ds.loc[ds[fl]]["trajectory"].unique()
    ds = ds.loc[ds["trajectory"].isin(ids)]

    # Get either convective or slantwise ascent start time
    ds = add_norm_time(ds, "time", "trajectory", None, fl)
    ds = xr.Dataset.from_dataframe(ds.set_index(["ensemble", "trajectory", "time"]))

    ds = add_sat(ds)

    for col in ds:
        if col in ["QR_IN", "QS_IN", "QI_IN", "QG_IN",
                   "QR_OUT", "QS_OUT", "QI_OUT", "QG_OUT",
                   "NR_IN", "NS_IN", "NI_IN", "NG_IN",
                   "NR_OUT", "NS_OUT", "NI_OUT", "NG_OUT"]:
            ds[col] = np.fabs(ds[col])

    # Calculate turbulence flux
    # q_tot(0) = q_tot(1) + flux_out(1) - flux_in(1) - flux_turb(1)
    ds = add_turb_flux(ds)

    # Calculate ascend velocity
    ds = add_ascend_velocity(ds)

    type_name = ""
    if fl == "conv_600":
        type_name = "Convective 600hPa"
    elif fl == "conv_400":
        type_name = "Convective 400hPa"
    elif fl == "slan_400":
        type_name = "Slantwise 400hPa"
    else:
        type_name = "Slantwise 600hPa"
    ds["type"] = type_name

    for col in ds:
        if col == "type":
            continue
        if col in ["conv_400", "conv_600", "slan_400",
                   "slan_600", "type", "WCB_flag", "dp2h"]:
            ds[col] = ds[col].fillna(False)
            ds[col] = ds[col].astype("bool", casting="unsafe", copy=False)
            continue

    ds.attrs = {
        "duration_in_sec": duration,
        "pollon": pollon,
        "pollat": pollat,
        "output_timestep_in_sec": output_timestep_in_sec,
        "cloud_type": 2723}

    ds["time"].attrs = {
        "standard_name": "time",
        "long_name":     "time",
        "units":  "seconds since " + ref_time,
        "trajectory_starttime": ref_time,
        "forecast_inittime": ref_time}
    ds["time_after_ascent"].attrs = {
        "standard_name": "time_after_ascent",
        "long_name": "time after rapid ascent started",
        "units": "seconds since start of convective/slantwise ascent"}
    ds["lon"].attrs = {
        "standard_name": "longitude",
        "long_name": "rotated longitude",
        "units": "degrees"}
    ds["lat"].attrs = {
        "standard_name": "latitude",
        "long_name": "rotated latitude",
        "units": "degrees"}
    ds["pressure"].attrs = {
        "standard_name": "air_pressure",
        "long_name": "pressure",
        "units": "Pa",
        "positive": "down",
        "axis": "z"}
    ds["z"].attrs = {
        "standard_name": "height",
        "long_name": "height above mean sea level",
        "units": "m AMSL"}
    ds["T"].attrs = {
        "standard_name": "air_temperature",
        "long_name": "temperature",
        "units": "K"}
    ds["S"].attrs = {
        "standard_name": "saturation",
        "long_name": "saturation",
        "units": "percentage"}
    ds["conv_400"].attrs = {
        "standard_name": "convective_400hPa_ascent",
        "long_name": "convective 400hPa ascent"}
    ds["conv_600"].attrs = {
        "standard_name": "convective_600hPa_ascent",
        "long_name": "convective 600hPa ascent"}
    ds["slan_400"].attrs = {
        "standard_name": "slantwise_400hPa_ascent",
        "long_name": "slantwise 400hPa ascent"}
    ds["slan_600"].attrs = {
        "standard_name": "slantwise_600hPa_ascent",
        "long_name": "slantwise 600hPa ascent"}
    ds["WCB_flag"].attrs = {
        "standard_name": "wcb_ascent",
        "long_name": "WCB ascent"}
    ds["dp2h"].attrs = {
        "standard_name": "2h_ascent_rate",
        "long_name": "2h ascent rate"}
    ds["w"].attrs = {
        "standard_name": "ascend_velocity",
        "long_name": "ascend velocity",
        "units": "m s^-1"}

    ds["QV"].attrs = {
        "standard_name": "specific_humidity",
        "long_name": "specific humidity",
        "units": "kg kg^-1"}
    ds["QC"].attrs = {
        "standard_name": "mass_fraction_of_cloud_liquid_water_in_air",
        "long_name": "specific cloud liquid water content",
        "units": "kg kg^-1"}
    ds["QR"].attrs = {
        "standard_name": "mass_fraction_of_rain_in_air",
        "long_name": "specific rain content",
        "units": "kg kg^-1"}
    ds["QS"].attrs = {
        "standard_name": "mass_fraction_of_snow_in_air",
        "long_name": "specific snow content",
        "units": "kg kg^-1"}
    ds["QI"].attrs = {
        "standard_name": "mass_fraction_of_cloud_ice_in_air",
        "long_name": "specific cloud ice content",
        "units": "kg kg^-1"}
    ds["QG"].attrs = {
        "standard_name": "mass_fraction_of_graupel_in_air",
        "long_name": "specific graupel content",
        "units": "kg kg^-1"}

    ds["QR_IN"].attrs = {
        "standard_name": "sedi_influx_of_rain",
        "long_name": "sedimentation (from above) of rain droplet mixing ratio",
        "units": "kg kg^-1 s^-1"}
    ds["QS_IN"].attrs = {
        "standard_name": "sedi_influx_of_snow",
        "long_name": "sedimentation (from above) of snow crystal mixing ratio",
        "units": "kg kg^-1 s^-1"}
    ds["QI_IN"].attrs = {
        "standard_name": "sedi_influx_of_cloud_ice",
        "long_name": "sedimentation (from above) of ice crystal mixing ratio",
        "units": "kg kg^-1 s^-1"}
    ds["QG_IN"].attrs = {
        "standard_name": "sedi_influx_of_graupel",
        "long_name": "sedimentation (from above) of graupel mixing ratio",
        "units": "kg kg^-1 s^-1"}

    ds["QR_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_rain",
        "long_name": "sedimentation of rain droplet mixing ratio",
        "units": "kg kg^-1 s^-1"}
    ds["QS_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_snow",
        "long_name": "sedimentation of snow crystal mixing ratio",
        "units": "kg kg^-1 s^-1"}
    ds["QI_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_cloud_ice",
        "long_name": "sedimentation of ice crystal mixing ratio",
        "units": "kg kg^-1 s^-1"}
    ds["QG_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_graupel",
        "long_name": "sedimentation of graupel mixing ratio",
        "units": "kg kg^-1 s^-1"}

    ds["NCCLOUD"].attrs = {
        "standard_name": "specif_number_of_cloud_droplets_in_air",
        "long_name": "specific cloud droplet number",
        "units": "kg^-1"}
    ds["NCRAIN"].attrs = {
        "standard_name": "specif_number_of_rain_drops_in_air",
        "long_name": "specific rain drop number",
        "units": "kg^-1"}
    ds["NCSNOW"].attrs = {
        "standard_name": "specif_number_of_snow_flakes_in_air",
        "long_name": "specific snow flake number",
        "units": "kg^-1"}
    ds["NCICE"].attrs = {
        "standard_name": "specif_number_of_cloud_ice_in_air",
        "long_name": "specific cloud ice number",
        "units": "kg^-1"}
    ds["NCGRAUPEL"].attrs = {
        "standard_name": "specif_number_of_graupel_in_air",
        "long_name": "specific graupel number",
        "units": "kg^-1"}

    ds["NR_IN"].attrs = {
        "standard_name": "sedi_influx_of_rain_number",
        "long_name": "sedimentation (from above) of specific rain drop number",
        "units": "kg^-1 s^-1"}
    ds["NS_IN"].attrs = {
        "standard_name": "sedi_influx_of_snow_number",
        "long_name": "sedimentation (from above) of specific snow flake number",
        "units": "kg^-1 s^-1"}
    ds["NI_IN"].attrs = {
        "standard_name": "sedi_influx_of_ics_number",
        "long_name": "sedimentation (from above) of specific cloud ice number",
        "units": "kg^-1 s^-1"}
    ds["NG_IN"].attrs = {
        "standard_name": "sedi_influx_of_graupel_number",
        "long_name": "sedimentation (from above) of specific graupel number",
        "units": "kg^-1 s^-1"}

    ds["NR_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_rain_number",
        "long_name": "sedimentation of rain droplet number",
        "units": "kg^-1 s^-1"}
    ds["NS_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_snow_number",
        "long_name": "sedimentation of snow crystal number",
        "units": "kg^-1 s^-1"}
    ds["NI_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_ice_number",
        "long_name": "sedimentation of ice crystal number",
        "units": "kg^-1 s^-1"}
    ds["NG_OUT"].attrs = {
        "standard_name": "sedi_outflux_of_graupel_number",
        "long_name": "sedimentation of graupel number",
        "units": "kg^-1 s^-1"}

    ds["Q_TURBULENCE"].attrs = {
        "standard_name": "turbulence_flux",
        "long_name": "flux from turbulence",
        "units": "kg^-1 s^-1"}
    ds["type"].attrs = {
        "standard_name": "trajectory_type",
        "long_name": "trajectory type"}

    comp = dict(zlib=True, complevel=9)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(
        path=store_path,
        encoding=encoding,
        compute=True,
        engine="netcdf4",
        format="NETCDF4",
        mode="w")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        path = sys.argv[1]
    if len(sys.argv) > 2:
        store_path = sys.argv[2]
        if store_path[-1] != "/":
            store_path += "/"
    if len(sys.argv) > 3:
        flags = [sys.argv[3]]
    else:
        flags = ["conv_400", "conv_600", "slan_400", "slan_600"]
    file_list = []
    for f in os.listdir(path):
        if ".nc_wcb" in f:
            file_list.append(os.path.join(path, f))
    file_list = np.sort(file_list)
    for flag in flags:
        print(f"#################### {flag} ######################")
        for i in pb.progressbar(range(len(file_list)), redirect_stdout=True):
            convert_wcb2(
                f = file_list[i],
                store_path=store_path + flag + "_" + str(i) + "_" + file_list[i].split("/")[-1],
                fl=flag,
                ensemble=i)

