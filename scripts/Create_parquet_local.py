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

file_type = sys.argv[1]
direc_path = sys.argv[2]
store_path = sys.argv[3]
filt = False
EPSILON = 0.0
ncpus = None
c_this = False
file_list = []
for f in os.listdir(direc_path):
    file_list.append(os.path.join(direc_path, f))
file_list = np.sort(file_list)
print("Running for {} trajectories".format(len(file_list)//37))

# wcb#####_traj#_MAP_t#####_p###
# where # is a number.
processed_trajectories = []
failed_trajectories = []

for f_this in file_list:
    if "diff" in f_this or "reference" in f_this:
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
    sim = " "
    traj = -1
    for f in files:
        if "diff" in f:
            load_f.append(f)
        else:
            if "reference" in f:
                ref = f
                traj = f.split("_traj")[1]
                traj = traj.split("_MAP")[0]
                traj = int(traj)
            else:
                sim = f
                suffix = f[:-4]
                suffix = suffix.split("_")
                suffix = suffix[-2] + "_" + suffix[-1]

    df_dic_mapped = Deriv(direc=direc_path,
                          filt=filt,
                          EPSILON=EPSILON,
                          trajectories=[traj],
                          file_list=load_f,
                          suffix=suffix,
                          threads=ncpus)
    print("Found reference: {}".format(ref))
    print("Found sim: {}".format(sim))
    df_sim_mapped = Sim()
    df_sim_mapped.load_file(
        filename=sim,
        sep=",",
        change_ref=True,
        refs=ref)
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
        
    print("Saving as {}".format(file_type))
    t = timer()
    try:
        if file_type == "parquet":
            df_dic_mapped.to_parquet(
                store_path, compression="snappy", low_memory=True)
        elif file_type == "netcdf":
            df_dic_mapped.to_netcdf(store_path)
        else:
            print("No such file format: {}".format(file_type))
            failed_trajectories.append(prefix)
    except Exception as e:
        print("FAILED: {}".format(prefix))
        print(str(e))
        failed_trajectories.append(prefix)
    t2 = timer()
    print("Saving done in {} s".format(t2-t))

print("Done with following trajectories:\n{}".format(processed_trajectories))
print("Failed the following trajectories:\n{}".format(failed_trajectories))
print("Finished")
