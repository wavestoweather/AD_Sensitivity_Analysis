import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
import progressbar
from loader import *

n_traj = 903
prefix = sys.argv[1] # rain_wcb180_traj0
suffix = sys.argv[2] # start_over_200

def show(df, start=0, step=500, n_steps=2000):
    print(np.shape(df))
    stop = n_steps*step + start
    with pd.option_context('display.max_rows', None,
                           'display.max_columns', None):
        print(df[start:stop:step])
        print("Delta z")

        dt = df.diff(periods=n_traj)
        if "z" in dt.keys():
            print(dt.z[start:stop:step])
        else:
            print("T")
            print(np.sum(dt["T"][1:2001]))
            print(np.sum(dt["T"][2001:4001]))
            print(dt["T"][start:stop:step])

        if "latent_heat" in dt.keys():
            print("latent_heat")
            print(np.sum(dt["latent_heat"][1:2001]))
            print(np.sum(dt["latent_heat"][2001:4001]))
            print("latent_cool")
            print(np.sum(dt["latent_cool"][1:2001]))
            print(np.sum(dt["latent_cool"][2001:4001]))



nc = load_nc(inp="/data/project/m2_jgu-tapt/"
             + "online_trajectories/wcb201609_vladiana/"
             + "O_WCB_all_20160922_00.nc")
print("########Trajectory#############")
# print(np.shape(nc))
# print(len(nc.ntra))
n_substeps = 200
show(nc, step=n_traj, n_steps=9)
# exit()
print("########start_over#############")
start_over = load_output(prefix=prefix, suffix=suffix)
show(start_over, step=1000)

exit()
# print("########diff#############")
# diff = pd.read_csv("diff_pd.csv")
# show(diff, step=1000)

# show(start_over, 0, 20, 1, 1)
print("#########start_over parts:")
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(start_over[0:5])
    print(start_over[n_substeps-5:n_substeps+5])
exit()

print("#########diff parts:")
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(diff[0:5])
    print(diff[1995:2005])

print("#########derivatives:")

deriv_data = load_derivatives(prefix=prefix, suffix="_start_over")
deriv_data = filter_zeros(deriv_data, EPSILON=0.0)
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    for key in deriv_data.keys():
        print("*+*+*+*+*+*+ {}:".format(key))
        print(deriv_data[key][0:5])
        print(deriv_data[key][1995:2005])
