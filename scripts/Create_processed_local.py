import os
import sys
sys.path.append(os.getcwd())

from loader import load_mult_derivates_directory
from loader import load_nc, rotate_df, norm_deriv, ratio_deriv
import loader as loader
import numpy as np
from Deriv import Deriv
from Sim import Sim
from timeit import default_timer as timer
import xarray as xr

file_type = sys.argv[1]
direc_path = sys.argv[2]
store_path = sys.argv[3]
input_type = sys.argv[4]
filt = False
EPSILON = 0.0
ncpus = None
c_this = False
file_list = []
for f in os.listdir(direc_path):
    file_list.append(os.path.join(direc_path, f))
file_list = np.sort(file_list)
print("Running for {} trajectories".format(len(file_list)//35))

# wcb#####_traj#_MAP_t#####_p###
# where # is a number.
processed_trajectories = []
failed_trajectories = []
datasets = {}
def parse_attr(f):
    print(f)
    attributes = {}
    key = None
    column_key = None
    opened = open(f, "r")
    for line in opened:
        print(line)
        if "[" in line and "]" in line:
            key = line.replace("[", "").replace("]", "")
            attributes[key] = []
        elif key == "Global attributes":
            attributes[key].append((line.split("=")[0], line.split("=")[1]))
        else:
            if "column" in line:
                column_key = line.split("=")[1]
                attributes[key][column_key] = []
            else:
                attributes[key][column_key].append((line.split("=")[0], line.split("=")[1]))
    return attributes

for f_this in file_list:
    if "diff" in f_this or "reference" in f_this or "attributes" in f_this:
        continue
    prefix = f_this.split(".tx")[0]
    if prefix in processed_trajectories:
        print("{} already processed. Continue".format(prefix))
        continue

    processed_trajectories.append(prefix)
    idx = np.argwhere([prefix in f for f in file_list]).flatten()
    if len(idx) == 0:
        print("Did not find {}".format(prefix))
    files = file_list[idx]
    t = timer()
    load_f = []
    ref = ""
    sim = ""
    att = None
    traj = -1

    for f in files:
        if "attributes" in f:
            att = f
        elif "diff" in f:
            load_f.append(f)
        elif "reference" in f:
            ref = f
            traj = f.split("_traj")[1]
            traj = traj.split("_MAP")[0]
            traj = int(traj)
        else:
            sim = f
            suffix = f[:-4]
            suffix = suffix.split("/")[-1]
            suffix = suffix.split("_")
            suffix = suffix[-2] + "_" + suffix[-1]
    print("Found reference: {}".format(ref))
    print("Found sim: {}".format(sim))
    if input_type == "MET3D":
        print(f"Found attributes: {att}")
    print

    df_dic_mapped = Deriv(direc=direc_path,
                          filt=filt,
                          EPSILON=EPSILON,
                          trajectories=[traj],
                          file_list=load_f,
                          suffix=suffix,
                          threads=ncpus)

    df_sim_mapped = Sim()
    df_sim_mapped.load_file(
        filename=sim,
        sep=",",
        change_ref=True,
        refs=ref,
        attr=att)
    t2 = timer()
    print("Loading done in {} s".format(t2-t))

    print("Checking for minimum number of rows")
    min_rows = len(df_sim_mapped.data.index)
    crop_data = False
    for key in df_dic_mapped.data:
        this_rows = len(df_dic_mapped.data[key].index)
        if this_rows < min_rows:
            min_rows = this_rows
            crop_data = True
    if crop_data:
        print(f"Some sensitivities lack some rows. We delete the last rows such that {min_rows} are left.")
        print("If anything fails from here on, please check your data!")
        df_sim_mapped.data = df_sim_mapped.data.head(min_rows)
        for key in df_dic_mapped.data:
            df_dic_mapped.data[key] = df_dic_mapped.data[key].head(min_rows)


    print("Get ratio of data")
    t = timer()
    df_dic_mapped.calculate_ratios()
    t2 = timer()
    print("ratio done in {} s".format(t2-t))

    print("Adding columns for output Parameter results")
    t = timer()
    df_dic_mapped.add_param_values(df_sim_mapped.data)
    t2 = timer()
    print("Adding finished in {} s".format(t2-t))

    if input_type != "MET3D":
        print("Shift the timesteps such that t=0 is the start of the ascent")
        t = timer()
        # Get the currently used flag from the filename
        flag = "conv_400"
        if "conv_600" in direc_path:
            flag = "conv_600"
        elif "slan_400" in direc_path:
            flag = "slan_400"
        elif "slan_600" in direc_path:
            flag = "slan_600"
        df_dic_mapped.shift_time(flag=flag)
        t2 = timer()
        print("Shifting done in {} s.".format(t2-t))
    if not input_type == "MET3D":
        print("Saving as {}".format(file_type))
    t = timer()
    try:
        if file_type == "parquet":
            df_dic_mapped.to_parquet(
                store_path, compression="snappy", low_memory=True, attr=att)
        elif file_type == "netcdf":
            f_name = store_path + "/" + suffix
            if input_type == "MET3D":
                if suffix in datasets:
                    datasets[suffix].append(df_dic_mapped.get_netcdf_ready_data(attr=att))
                else:
                    datasets[suffix] = [df_dic_mapped.get_netcdf_ready_data(attr=att)]
            else:
                df_dic_mapped.to_netcdf(f_name)
        else:
            print("No such file format: {}".format(file_type))
            failed_trajectories.append(prefix)
    except Exception as e:
        print("FAILED: {}".format(prefix))
        print(str(e))
        failed_trajectories.append(prefix)
    t2 = timer()
    if not input_type == "MET3D":
        print("Saving done in {} s".format(t2-t))
if input_type == "MET3D":
    t = timer()
    for suffix in datasets:
        for i, data in enumerate(datasets[suffix]):
            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in data.data_vars}
            f_name = store_path + "/" + suffix
            data.to_netcdf(
                    f_name + f"_derivs_{i:03}.nc_wcb",
                    encoding=encoding,
                    compute=True,
                    engine="netcdf4",
                    format="NETCDF4",
                    mode="w")
    t2 = timer()
    print("Saving done in {} s".format(t2-t))


print("Done with following trajectories:\n{}".format(processed_trajectories))
if len(failed_trajectories) > 0:
    print("Failed the following trajectories:\n{}".format(failed_trajectories))
print("Finished")
