import numpy as np
from pylab import rcParams
from scripts.Deriv_dask import Deriv_dask
import scripts.latexify as latexify
from timeit import default_timer as timer
import sys

rcParams['figure.figsize'] = (16,10)
latexify.set_size(beamer=True)
kwargs = {"alpha": 0.4}
trajectories = None # [0, 1, 2, 3, 4, 5]
# out_params = ["qi", "qs", "qg", "qh", "qv", "qc", "qr", "p", "T", "S"]
# in_params = None
# S variant
out_params = ["qc", "S", "qr", "p"]
in_params = ["dcloud_min_x_act",
    "dice_min_x_nuc_hetero",
    "dgraupel_min_x_melt", "drain_min_x_sedimentation", "dsnow_min_x_melt",
    "dice_min_x_depo", "dsnow_min_x_depo",
    "dgraupel_min_x_evap", "dgraupel_min_x_sedimentation", "de_2"]
nplots = 3

direc_path = "/data/project/wcb/wcb_traj_flag_deriv_mult_min2/parquet_concat"
# direc_path = "/data/project/wcb/testcase"

columns = ["timestep", "Output Parameter", "trajectory"] + in_params + out_params

t = timer()
data = Deriv_dask(
    direc=direc_path,
    parquet=True,
    columns=columns
)
t2 = timer()
print("Loading done in {}s".format(t2-t))
print("Got trajectories: {}".format(data.data["trajectory"].unique().compute()))
print(data.data.head())
# print(data.data.loc[data.data["trajectory"] == 0].head())
# print(data.data.head())

# t = timer()
# # Standard
# data.plot_same_orders(
#     out_params=out_params,
#     mapped=False,
#     scatter=False,
#     in_params=in_params,
#     x_axis="timestep",
#     n_plots=nplots,
#     verbose=False,
#     frac=0.001,
#     **kwargs
# )
# t2 = timer()
# print("Plotting done in {}s".format(t2-t))
# out_params = ["qc"]
# print("Sorted tuples (MAPPED) for out_params: {}".format(out_params))
# print(data.get_sorting(out_params=out_params, in_params=None, mapped=True))

# for p in ["qv"]:
#     print("\nSorted tuples (MAPPED) for out_param: {}".format(p), flush=True)
#     print(data.get_sorting(out_params=[p], in_params=None, mapped=True), flush=True)


# # With scatter plot
# print("Scatter plot")
# t = timer()
# data.plot_same_orders(
#     out_params=out_params,
#     mapped=False,
#     scatter=True,
#     in_params=in_params,
#     x_axis="timestep",
#     n_plots=nplots,
#     verbose=False,
#     frac=0.01,
#     **kwargs
# )
# t2 = timer()
# print("Scatter plot done in {}s".format(t2-t))

# print("Plot with p")
# t = timer()
# # # Same but with pressure as x-axis
# data.plot_same_orders(
#     out_params=out_params,
#     mapped=True,
#     scatter=False,
#     in_params=in_params,
#     x_axis="p",
#     n_plots=nplots,
#     frac=0.01,
#     **kwargs
# )
# t2 = timer()
# print("Plot with p done in {}s".format(t2-t))

# print("Scatter plot with p")
# t = timer()
# # # With scatter plot
# data.plot_same_orders(
#     out_params=out_params,
#     mapped=True,
#     scatter=True,
#     in_params=in_params,
#     x_axis="p",
#     n_plots=nplots,
#     frac=0.01,
#     **kwargs
# )
# t2 = timer()
# print("Scatter plot with p done in {}s".format(t2-t))

# Plot violin plots
# data.plot_mapped(
#     out_params=out_params,
#     in_params=in_params,
#     kind="violin",
#     n_plots=nplots
# )
# Takes way too long
# data.plot_mapped(
#     out_params=out_params,
#     in_params=in_params,
#     kind="swarm",
#     n_plots=nplots
# )
# Looks weird
# data.plot_mapped(
#     out_params=out_params,
#     in_params=in_params,
#     kind="bar",
#     n_plots=nplots
# )

# print("Plot two")
# t = timer()
# # Works so far
# data.plot_two(
#     in_params=in_params,
#     out_params=out_params,
#     x_axis="timestep",
#     mapped=False,
#     trajectories=None,
#     scatter=False,
#     frac=None,
#     percentile=None,
#     # min_x=0,
#     # max_x=187940,
#     # nth=10000,
#     line_deriv=True,
#     **kwargs
# )
# t2 = timer()
# print("Plot two done in {}s\n".format(t2-t))
# ################
# t = timer()
# # Works so far
# if sys.argv[1] == "2":
#     data.plot_two_ds2(
#         in_params=in_params,
#         out_params=out_params,
#         x_axis="p",
#         mapped=False,
#         trajectories=trajectories,
#         scatter=False,
#         frac=None,
#         percentile=None,
#         # min_x=20000,
#         # max_x=40000,
#         # nth=10000,
#         line_deriv=True,
#         **kwargs
#     )
#     t2 = timer()
#     print("Plot two (p) done in {}s\n".format(t2-t))
# elif sys.argv[1] == "3":
#     data.plot_two_ds3(
#         in_params=in_params,
#         out_params=out_params,
#         x_axis="p",
#         mapped=False,
#         trajectories=trajectories,
#         scatter=False,
#         frac=None,
#         percentile=None,
#         # min_x=20000,
#         # max_x=40000,
#         # nth=10000,
#         line_deriv=True,
#         **kwargs
#     )
#     t2 = timer()
#     print("Plot two (p) done in {}s\n".format(t2-t))
# else:
#     data.plot_two_ds(
#         in_params=in_params,
#         out_params=out_params,
#         x_axis="p",
#         mapped=False,
#         trajectories=trajectories,
#         scatter=False,
#         frac=None,
#         percentile=None,
#         # min_x=20000,
#         # max_x=40000,
#         # nth=10000,
#         line_deriv=True,
#         **kwargs
#     )
#     t2 = timer()
#     print("Plot two (p) done in {}s\n".format(t2-t))

t = timer()
# Works so far
data.plot_histo(
    in_params=in_params,
    out_params=out_params,
    x_axis="p",
    trajectories=trajectories,
    hex=True,
    logz=True
)
t2 = timer()
print("Plot hex (p) done in {}s\n".format(t2-t))
###### With percentile
# print("Plot two with percentile")
# t = timer()
# # Works so far
# data.plot_two(
#     in_params=in_params,
#     out_params=out_params,
#     x_axis="timestep",
#     mapped=False,
#     trajectories=None,
#     scatter=False,
#     frac=None,
#     percentile=68,
#     # min_x=0,
#     # max_x=187940,
#     # nth=10000,
#     line_deriv=True,
#     **kwargs
# )
# t2 = timer()
# print("Plot two done in {}s\n".format(t2-t))

# t = timer()
# # Works so far
# data.plot_two(
#     in_params=in_params,
#     out_params=out_params,
#     x_axis="p",
#     mapped=False,
#     trajectories=None,
#     scatter=False,
#     frac=None,
#     percentile=68,
#     # min_x=20000,
#     # max_x=40000,
#     # nth=10000,
#     line_deriv=True,
#     **kwargs
# )
# t2 = timer()
# print("Plot two (p) done in {}s\n".format(t2-t))

# print("Plot two with p")
# t = timer()
# data.plot_two(
#     in_params=in_params,
#     out_params=["p", "qc", "S"],
#     x_axis="p",
#     mapped=False,
#     trajectories=None,
#     frac=None,
#     scatter=False,
#     nth=100,
#     **kwargs
# )
# t2 = timer()
# print("Plot two with p done in {}s".format(t2-t))
