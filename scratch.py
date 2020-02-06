from scripts.loader import load_mult_derivates_directory
from scripts.loader import load_nc, rotate_df, norm_deriv, ratio_deriv
import scripts.loader as loader
import scripts.plot_mapped as plotter
import scripts.latexify as latexify
import numpy as np
from pylab import rcParams
from scripts.Deriv import Deriv

rcParams['figure.figsize'] = (16,10)
latexify.set_size(beamer=True)


traj = 0
direc_path = "/data/project/wcb/sim_min2/"
max_time = None
filt = True
EPSILON = 0.0
trajectory = [traj]
print("Load derivatives")
df_dic_mapped = Deriv(direc=direc_path,
                      filt=False,
                      EPSILON=EPSILON,
                      trajectories=trajectory,
                      suffix=None)
print("Loading derivatives done")
print("Loading output parameter values")
df_sim_mapped = loader.load_output(
    "/data/project/wcb/sim_min2/wcb187920.0_traj{}_start_over_MAP_Flux_2_t000000_p001.txt".format(traj, traj),
    sep=",",
    refs="/data/project/wcb/sim_min2/wcb187920.0_traj{}_start_over_MAP_Flux_2_t000000_p001_reference_values.txt".format(traj, traj))
print("Loading parameters done")
print("Look only at mapped data")
df_sim_mapped = df_sim_mapped[df_sim_mapped.MAP == True]
df_dic_mapped.delete_not_mapped()
print("Cutting data done")
####### DEBUG
n = np.Inf
for key in df_dic_mapped.get_out_params():
    n_tmp = len(df_dic_mapped.data[key].index)
    if n_tmp < n:
        n = n_tmp
if len(df_sim_mapped.index) < n:
    n = len(df_sim_mapped.index)
for key in df_dic_mapped.get_out_params():
    df_dic_mapped.data[key] = df_dic_mapped.data[key].iloc[:n]
df_sim_mapped = df_sim_mapped.iloc[:n]
df_dic_mapped.n_timesteps = n
####### END DEBUG
print("Get ratio of data")
df_dic_mapped.calculate_ratios()
print("ratio done")
kwargs = {"alpha": 0.4}
x_axis_list = ["p", "T", "qr", "w", "z"]
df_dic_mapped.plot_same_orders(**kwargs)
for x in x_axis_list:
    df_dic_mapped.add_param_values(x, df_sim_mapped[x].tolist())
    df_dic_mapped.plot_same_orders(x_axis=x, **kwargs)
