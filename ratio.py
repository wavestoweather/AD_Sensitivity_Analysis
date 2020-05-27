import os

from scripts.loader import load_mult_derivates_directory
from scripts.loader import load_nc, rotate_df, norm_deriv, ratio_deriv
import scripts.loader as loader
import scripts.plot_mapped as plotter
import scripts.latexify as latexify
import numpy as np
from pylab import rcParams
from scripts.Deriv import Deriv
from scripts.Sim import Sim

rcParams['figure.figsize'] = (16,10)
# latexify.set_size(beamer=True)
kwargs = {"alpha": 0.4}

# direc_path = "/data/project/wcb/sim_min2_fixed"
direc_path = "/data/project/wcb/wcb_traj_flag_deriv_mult_min"
max_time = None
filt = True
EPSILON = 0.0

for t in [2]:
    c_this = False
    trajectory = [t]
    ################################## LOADING AND PREPROCESSING DATA
    print("Load derivatives for {}".format(t))
    df_dic_mapped = Deriv(direc=direc_path,
                        filt=False,
                        EPSILON=EPSILON,
                        trajectories=trajectory,
                        suffix=None,
                        threads=4)
    print("Loading derivatives done")
    print("Loading output parameter values")
    df_sim_mapped = Sim()
    df_sim_mapped.load_path(direc_path, sep=",", trajectories=trajectory)
    print("Loading parameters done")
    quit()
    # That was for broken datasets...
    # for k in df_dic_mapped.data.keys():

    #     df = df_dic_mapped.data[k].loc[np.isnan(df_dic_mapped.data[k]["MAP"])]
    #     if not df.empty:
    #         print(k)
    #         print(df)
    #         print(df_dic_mapped.data[k]["MAP"].unique())
    #         c_this = True
    #         break

    # if c_this:
    #     continue
    #### End broken datasets
    print("Get ratio of data")
    df_dic_mapped.calculate_ratios()
    print("ratio done")
    print("Adding columns for output Parameter results")
    df_dic_mapped.add_param_values(df_sim_mapped.data)
    print("Adding finished")
    print("Saving as parquet")
    df_dic_mapped.to_parquet(direc_path + "/wcb_traj{}_ratio".format(t))
    print("Done")


