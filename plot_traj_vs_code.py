import numpy as np
import xarray as xr
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd
import sys
import time as timer


matplotlib.rcParams['agg.path.chunksize'] = 10000
get_time = True
if get_time:
    start = timer.time()
data = None
data_fixed = None
data_restart = None
data_restart_fixed = None
prefix = sys.argv[1]
prefix_title = sys.argv[1]

sns.set_style("darkgrid")

# Load the references
refs = np.loadtxt('reference_values.txt')
Tref = refs[0]
pref = refs[1]
qref = refs[2]
Nref = refs[3]
wref = refs[4]
tref = refs[5]


def ref(data):
    data["p"] = pref*data["p"]
    data["T"] = Tref*data["T"]
    data["qr"] = qref*data["qr"]
    data["qc"] = qref*data["qc"]
    data["qv"] = qref*data["qv"]
    data["qi"] = qref*data["qi"]
    data["qs"] = qref*data["qs"]
    data["qg"] = qref*data["qg"]
    data["w"] = wref*data["w"]
    data["timestep"] = tref*data["timestep"]


try:
    data = pd.read_csv("data/" + prefix + ".txt", sep=",")
    data = data[data["trajectory"] == 0]
    ref(data)
except:
    pass

try:
    data_restart = pd.read_csv("data/" + prefix + "_start_over.txt", sep=",")
    data_restart = data_restart[data_restart["trajectory"] == 0]
    ref(data_restart)
except:
    pass

try:
    data_restart_fixed = pd.read_csv("data/" + prefix
                                     + "_start_over_fixed.txt", sep=",")
    data_restart_fixed = (data_restart_fixed[data_restart_fixed["trajectory"]
                                             == 0]
    ref(data_restart_fixed)
except:
    pass

try:
    data_fixed = pd.read_csv("data/" + prefix + "_fixed.txt", sep=",")
    data_fixed = data_fixed[data_fixed["trajectory"] == 0]
    ref(data_fixed)
except:
    pass

inp = ("/mnt/localscratch/data/project/m2_jgu-tapt/" +
       "online_trajectories/foehn201305_case/foehn201305_warming1.nc")
ds = xr.open_dataset(inp)
data_nc = ds.to_dataframe()
data_nc.reset_index(drop=True, inplace=True)

data_nc = data_nc.drop(index=[0])
n_rows = len(data_nc.index)
n_times = len(ds["time"])
n_trajectories = len(ds["id"])
time = [20.0*(i%n_times) for i in range(n_rows)]
data_nc["timestep"] = time
numpy_traj = ds["id"]
trajectory = np.asarray([numpy_traj[i//n_times] for i in range(n_rows)])
data_nc["trajectory"] = trajectory
data_nc = data_nc[data_nc["trajectory"] == 0]
data_nc["p"] = data_nc["p"]*100 # Trajectories are in hPa

if data is not None:
    end_time = data["timestep"].max()
elif data_fixed is not None:
    end_time = data_fixed["timestep"].max()
elif data_restart_fixed is not None:
    end_time = data_restart_fixed["timestep"].max()
else:
    end_time = data_nc["timestep"].max()
data_nc = data_nc[data_nc["timestep"] <= end_time]
data_nc = data_nc.rename(columns={"t": "T"})

max_time = data_nc["timestep"].max()
try:
    data = data[data["timestep"] <= max_time]
except:
    pass
try:
    data_restart = data_restart[data_restart["timestep"] <= max_time]
except:
    pass
try:
    tmp = data_restart_fixed[data_restart_fixed["trajectory"] == 0]
    data_restart_fixed = data_restart_fixed[data_restart_fixed["timestep"]
                                            <= max_time]
    tmp = data_restart_fixed[data_restart_fixed["trajectory"] == 0]
except:
    pass
if get_time:
    end = timer.time()
    print("Loading data in {} seconds".format(end-start))
# The first row is *often times* just zero
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    if data is not None:
        print(data[0:25])
    if data_restart is not None:
        print("######restart#############")
        print(data_restart[0:25])
    if data_restart_fixed is not None:
        print("#######restartfixed############")
        print(data_restart_fixed[-25:-1])
        tmp = data_restart_fixed[data_restart_fixed["trajectory"] == 0]
        # print(tmp[-25:-1])
        # print(tmp["timestep"])
        # print(data_restart_fixed["timestep"])
    print("#######nc##########")
    print(data_nc[0:25])
# print(data_restart["timestep"])
def plot(y, data, data_restart=None, data_nc=None, data_fixed=None,
    data_restart_fixed=None):

    fig, ax = plt.subplots()
    if not data is None:
        sns.lineplot("timestep", y, data=data, ax=ax, label="w/o restart")
    if not data_restart is None:
        sns.lineplot("timestep", y, data=data_restart, ax=ax,
                     label="with restart")
    if not data_fixed is None:
        sns.lineplot("timestep", y, data=data_fixed, ax=ax,
                     label="w/o restart, fixed")
    if not data_restart_fixed is None:
        sns.lineplot("timestep", y, data=data_restart_fixed, ax=ax,
                     label="with restart, fixed")
    if not data_nc is None:
        sns.scatterplot("timestep", y, data=data_nc, ax=ax,
                        label="trajectory")
        x_ticks = np.arange(0, data_nc["timestep"].max() + 19, 20)
        ax.set_xticks(x_ticks)
    ax.lines[1].set_linestyle("--")
    ax.lines[2].set_linestyle("-.")
    y_ticks = np.arange(-100, 100, 10)
    ax.set_yticks(y_ticks)
    # ax.lines[3].set_linestyle(":")
    plt.tight_layout()
    plt.legend()
    i = 0
    save = ("pics/" + prefix_title + "_diff_traj_" + y
            + "_" + "{:03d}".format(i) + ".png")
    while os.path.isfile(save):
        i = i+1
        save = ("pics/" + prefix_title + "_diff_traj_" + y
                + "_" + "{:03d}".format(i) + ".png")
    print("Saving to " + save)
    plt.savefig(save, dpi=300)
    plt.close()

def plot_differences(y, data):
    fig, ax = plt.subplots()
    if get_time:
        start = timer.time()
    sns.lineplot("time [s]", y, data=data, ax=ax, hue="method", alpha=0.5)
    if get_time:
        end = timer.time()
        print("Plot within {} seconds".format(end-start))
        start = end
    ax.set_ylabel(y)
    x_ticks = np.arange(0, data["time [s]"].max() + 19, 20)
    ax.set_xticks(x_ticks)
    y_ticks = np.arange(-100, 100, 10)
    ax.set_yticks(y_ticks)
    # Calculate the total difference
    total_diff = 0.0
    diff_restart = data[data["method"] == "with restart"]

    times = diff_restart["time [s]"].to_list()
    diff = diff_restart[y].to_list()
    if get_time:
        end = timer.time()
        print("Stuff done in {} seconds".format(end-start))
        start = end
    for i in range(len(diff_restart["time [s]"])-1):
        total_diff = (total_diff
            + (times[i+1] - times[i])
            * ((diff[i] + diff[i+1])/2.0))
    total_diff = total_diff/(data["time [s]"].max()-data["time [s]"].min())
    if get_time:
        end = timer.time()
        print("Got difference within {} seconds".format(end-start))
        start = end
    title = (prefix_title + ": avg difference " + y
             + " (with restart) {:.2f}".format(total_diff))
    plt.title(title)
    plt.tight_layout()
    plt.legend()
    i = 0
    save = ("pics/" + prefix_title + "_avg_diff_"
            + y.replace('$', '').replace("%", '').replace('[', '').replace(']', '').replace("\\", '').replace(' ', '_') + "_" + "{:03d}".format(i) + ".png")
    while os.path.isfile(save):
        i = i+1
        save = ("pics/" + prefix_title + "_avg_diff_" + y.replace('$', '').replace("%", '').replace('[', '').replace(']', '').replace("\\", '').replace(' ', '_') + "_" + "{:03d}".format(i) + ".png")

    print("Saving to " + save)
    plt.savefig(save, dpi=300)
    if get_time:
        end = timer.time()
        print("Finsihed saving in {} seconds".format(end-start))
    d = data[y]
    print(d[2000*7-10:2000*7+10])
    plt.close()

# plot("p", data, data_restart, data_nc, data_fixed, data_restart_fixed)
# plot("T", data, data_restart, data_nc, data_fixed, data_restart_fixed)
# plot("qr", data, data_restart, data_nc, data_fixed, data_restart_fixed)
# plot("qc", data, data_restart, data_nc, data_fixed, data_restart_fixed)
# plot("qv", data, data_restart, data_nc, data_fixed, data_restart_fixed)
# plot("w", data, data_restart, data_nc, data_fixed, data_restart_fixed)

def interpolate(data, key, n):
    steps = data[key].to_list()
    new_steps = np.linspace(steps[0], steps[-1], (steps[-1]-steps[0])/n)
    ret = np.interp(new_steps, steps, data)

    new_vals = np.asarray([])
    for i in range(len(steps)-1):
        start = steps[i]
        end = steps[i+1]
        m = (end-start)/n
        for j in range(n):
            new_vals = np.append(new_vals, start + m*j)
        if i == len(steps)-2:
            new_vals = np.append(new_vals, end)
    return new_vals


# Create dataframe that holds the differences to the trajectory
# The data of the trajectory is interpolated between each timestep
data_diff_dict = {"time [s]": np.asarray([]),
                  r"$\Delta p [\%]$": np.asarray([]),
                  r"$\Delta T [\%]$": np.asarray([]),
                  r"$\Delta q_r [\%]$": np.asarray([]),
                  r"$\Delta q_c [\%]$": np.asarray([]),
                  r"$\Delta q_v [\%]$": np.asarray([]),
                  r"$\Delta q_i [\%]$": np.asarray([]),
                  r"$\Delta q_s [\%]$": np.asarray([]),
                  r"$\Delta q_g [\%]$": np.asarray([]),
                  r"$\Delta w [\%]$": np.asarray([]),
                  "method": np.asarray([])}
if data is None and data_restart_fixed is None:
    data_diff_dict["time [s]"] = data_restart["timestep"]
    data_diff_dict["method"] = (
        ["with restart" for i in range(len(data_restart["timestep"]))])
elif data is None:
    data_diff_dict["time [s]"] = np.append(data_restart_fixed["timestep"],
                                           data_restart["timestep"])
    data_diff_dict["method"] = np.append(["with restart, fixed"
        for i in range(len(data_restart_fixed["timestep"]))],
        ["with restart" for i in range(len(data_restart["timestep"]))])
else:
    data_diff_dict["time [s]"] = np.append(data["timestep"],
        np.append(data_restart_fixed["timestep"],
        data_restart["timestep"]))
    data_diff_dict["method"] = ["w/o restart"
                                for i in range(len(data["timestep"]))]
    data_diff_dict["method"] = np.append(data_diff_dict["method"],
        np.append(["with restart, fixed"
                   for i in range(len(data_restart_fixed["timestep"]))],
        ["with restart" for i in range(len(data_restart["timestep"]))]))

data_nc_dict = {"p": np.asarray([]),
                "T": np.asarray([]),
                "qr": np.asarray([]),
                "qv": np.asarray([]),
                "qc": np.asarray([]),
                "qi": np.asarray([]),
                "qs": np.asarray([]),
                "qg": np.asarray([]),
                "w": np.asarray([]),
                "timestep": np.asarray([])}
if get_time:
    end2 = timer.time()
    print("Creating data_diff_dict in {} seconds".format(end2-end))
if data is None and data_restart_fixed is None:
    print("iiiiiiiiiiiiiiiiiiii")
    print("interpolate from {} to {} with timestep {}, {}".format(
        data_nc["p"].iloc[-2], data_nc["p"].iloc[-1],
        data_restart["timestep"].iloc[-2], data_restart["timestep"].iloc[-1]
    ))
    print(data_restart["timestep"].iloc[-25::])

    data_nc_dict["p"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["p"])
    data_nc_dict["T"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["T"])
    data_nc_dict["qv"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qv"])
    data_nc_dict["qr"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qr"])
    data_nc_dict["qc"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qc"])
    data_nc_dict["qi"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qi"])
    data_nc_dict["qs"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qs"])
    data_nc_dict["qg"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qg"])
    data_nc_dict["w"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["w"])
    data_nc_dict["timestep"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["timestep"])
elif data is None:

    data_nc_dict["p"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["p"])
    data_nc_dict["T"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["T"])
    data_nc_dict["qv"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qv"])
    data_nc_dict["qr"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qr"])
    data_nc_dict["qc"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qc"])
    data_nc_dict["qi"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qi"])
    data_nc_dict["qs"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qs"])
    data_nc_dict["qg"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["qg"])
    data_nc_dict["w"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["w"])
    data_nc_dict["timestep"] = np.interp( data_restart["timestep"],
        data_nc["timestep"], data_nc["timestep"])
else:
    data_nc_dict["p"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["p"])
    data_nc_dict["T"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["T"])
    data_nc_dict["qv"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["qv"])
    data_nc_dict["qr"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["qr"])
    data_nc_dict["qc"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["qc"])
    data_nc_dict["qi"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["qi"])
    data_nc_dict["qs"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["qs"])
    data_nc_dict["qg"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["qg"])
    data_nc_dict["w"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["w"])
    data_nc_dict["timestep"] = np.interp( data["timestep"],
        data_nc["timestep"] , data_nc["timestep"])

if get_time:
    end = timer.time()
    print("Interpolating done in {} seconds".format(end-end2))
data_nc = pd.DataFrame.from_dict(data_nc_dict)
if data is None and data_restart_fixed is None:
    data_diff_dict[r"$\Delta p [\%]$"] = (
       (np.asarray(data_restart["p"])
        - np.asarray(data_nc["p"]))*100 / np.asarray(data_nc["p"]))
    data_diff_dict[r"$\Delta T [\%]$"] = (
        (np.asarray(data_restart["T"])
         - np.asarray(data_nc["T"]))*100 / np.asarray(data_nc["T"]))
    data_diff_dict[r"$\Delta q_r [\%]$"] = (
        (np.asarray(data_restart["qr"])
         - np.asarray(data_nc["qr"]))*100 / np.asarray(data_nc["qr"]))
    data_diff_dict[r"$\Delta q_v [\%]$"] = (
        (np.asarray(data_restart["qv"])
         - np.asarray(data_nc["qv"]))*100 / np.asarray(data_nc["qv"]))
    data_diff_dict[r"$\Delta q_c [\%]$"] = (
        (np.asarray(data_restart["qc"])
         - np.asarray(data_nc["qc"]))*100 / np.asarray(data_nc["qc"]))
    data_diff_dict[r"$\Delta q_i [\%]$"] = (
        (np.asarray(data_restart["qi"])
         - np.asarray(data_nc["qi"]))*100 / np.asarray(data_nc["qi"]))
    data_diff_dict[r"$\Delta q_s [\%]$"] = (
        (np.asarray(data_restart["qs"])
         - np.asarray(data_nc["qs"]))*100 / np.asarray(data_nc["qs"]))
    data_diff_dict[r"$\Delta q_g [\%]$"] = (
        (np.asarray(data_restart["qg"])
         - np.asarray(data_nc["qg"]))*100 / np.asarray(data_nc["qg"]))
    data_diff_dict[r"$\Delta w [\%]$"] = (
        (np.asarray(data_restart["w"])
         - np.asarray(data_nc["w"]))*100 / np.asarray(data_nc["w"]))
elif data is None:
    tmp = (data_restart_fixed["p"] - data_nc["p"])*100 / data_nc["p"]
    data_diff_dict[r"$\Delta p [\%]$"] = np.append(
       (np.asarray(data_restart_fixed["p"])
        - np.asarray(data_nc["p"]))*100 / np.asarray(data_nc["p"]),
       (np.asarray(data_restart["p"])
        - np.asarray(data_nc["p"]))*100 / np.asarray(data_nc["p"]))
    data_diff_dict[r"$\Delta T [\%]$"] = np.append(
        (np.asarray(data_restart_fixed["T"])
         - np.asarray(data_nc["T"]))*100 / np.asarray(data_nc["T"]),
        (np.asarray(data_restart["T"])
         - np.asarray(data_nc["T"]))*100 / np.asarray(data_nc["T"]))
    data_diff_dict[r"$\Delta q_r [\%]$"] = np.append(
        (np.asarray(data_restart_fixed["qr"])
         - np.asarray(data_nc["qr"]))*100 / np.asarray(data_nc["qr"]),
        (np.asarray(data_restart["qr"])
         - np.asarray(data_nc["qr"]))*100 / np.asarray(data_nc["qr"]))
    data_diff_dict[r"$\Delta q_v [\%]$"] = np.append(
        (np.asarray(data_restart_fixed["qv"])
         - np.asarray(data_nc["qv"]))*100 / np.asarray(data_nc["qv"]),
        (np.asarray(data_restart["qv"])
         - np.asarray(data_nc["qv"]))*100 / np.asarray(data_nc["qv"]))
    data_diff_dict[r"$\Delta q_c [\%]$"] = np.append(
        (np.asarray(data_restart_fixed["qc"])
         - np.asarray(data_nc["qc"]))*100 / np.asarray(data_nc["qc"]),
        (np.asarray(data_restart["qc"])
         - np.asarray(data_nc["qc"]))*100 / np.asarray(data_nc["qc"]))
    data_diff_dict[r"$\Delta q_i [\%]$"] = np.append(
        (np.asarray(data_restart_fixed["qi"])
         - np.asarray(data_nc["qi"]))*100 / np.asarray(data_nc["qi"]),
        (np.asarray(data_restart["qi"]) - np.asarray(data_nc["qi"]))*100 / np.asarray(data_nc["qi"]))
    data_diff_dict[r"$\Delta q_s [\%]$"] = np.append(
        (np.asarray(data_restart_fixed["qs"])
         - np.asarray(data_nc["qs"]))*100 / np.asarray(data_nc["qs"]),
        (np.asarray(data_restart["qs"])
         - np.asarray(data_nc["qs"]))*100 / np.asarray(data_nc["qs"]))
    data_diff_dict[r"$\Delta q_g [\%]$"] = np.append(
        (np.asarray(data_restart_fixed["qg"])
         - np.asarray(data_nc["qg"]))*100 / np.asarray(data_nc["qg"]),
        (np.asarray(data_restart["qg"])
         - np.asarray(data_nc["qg"]))*100 / np.asarray(data_nc["qg"]))
    data_diff_dict[r"$\Delta w [\%]$"] = np.append(
        (np.asarray(data_restart_fixed["w"])
         - np.asarray(data_nc["w"]))*100 / np.asarray(data_nc["w"]),
        (np.asarray(data_restart["w"])
         - np.asarray(data_nc["w"]))*100 / np.asarray(data_nc["w"]))
else:
    data_diff_dict[r"$\Delta p [\%]$"] = np.append(
        (data["p"] - data_nc["p"]) / data_nc["p"],
        np.append((data_restart_fixed["p"] - data_nc["p"]) / data_nc["p"],
        (data_restart["p"] - data_nc["p"]) / data_nc["p"]))
    data_diff_dict[r"$\Delta T [\%]$"] = np.append(
        (data["T"] - data_nc["T"]) / data_nc["T"],
        np.append((data_restart_fixed["T"] - data_nc["T"]) / data_nc["T"],
        (data_restart["T"] - data_nc["T"]) / data_nc["T"]))
    data_diff_dict[r"$\Delta q_r [\%]$"] = np.append(
        (data["qr"] - data_nc["qr"]) / data_nc["qr"],
        np.append((data_restart_fixed["qr"] - data_nc["qr"]) / data_nc["qr"],
        (data_restart["qr"] - data_nc["qr"]) / data_nc["qr"]))
    data_diff_dict[r"$\Delta q_v [\%]$"] = np.append(
        (data["qv"] - data_nc["qv"]) / data_nc["qv"],
        np.append((data_restart_fixed["qv"] - data_nc["qv"]) / data_nc["qv"],
        (data_restart["qv"] - data_nc["qv"]) / data_nc["qv"]))
    data_diff_dict[r"$\Delta q_c [\%]$"] = np.append(
        (data["qc"] - data_nc["qc"]) / data_nc["qc"],
        np.append((data_restart_fixed["qc"] - data_nc["qc"]) / data_nc["qc"],
        (data_restart["qc"] - data_nc["qc"]) / data_nc["qc"]))
    data_diff_dict[r"$\Delta q_i [\%]$"] = np.append(
        (data["qi"] - data_nc["qi"]) / data_nc["qi"],
        np.append((data_restart_fixed["qi"] - data_nc["qi"]) / data_nc["qi"],
        (data_restart["qi"] - data_nc["qi"]) / data_nc["qi"]))
    data_diff_dict[r"$\Delta q_s [\%]$"] = np.append(
        (data["qs"] - data_nc["qs"]) / data_nc["qs"],
        np.append((data_restart_fixed["qs"] - data_nc["qs"]) / data_nc["qs"],
        (data_restart["qs"] - data_nc["qs"]) / data_nc["qs"]))
    data_diff_dict[r"$\Delta q_g [\%]$"] = np.append(
        (data["qg"] - data_nc["qg"]) / data_nc["qg"],
        np.append((data_restart_fixed["qg"] - data_nc["qg"]) / data_nc["qg"],
        (data_restart["qg"] - data_nc["qg"]) / data_nc["qg"]))
    data_diff_dict[r"$\Delta w [\%]$"] = np.append(
        (data["w"] - data_nc["w"]) / data_nc["w"],
        np.append((data_restart_fixed["w"] - data_nc["w"]) / data_nc["w"],
        (data_restart["w"] - data_nc["w"]) / data_nc["w"]))

try:
    print("\np\nLength of restart: {}\n".format(len(data_restart["p"]))
          + "Length of restart fixed: {}\n".format(len(data_restart_fixed["p"]))
          + "Length of data_diff_dict: {}\n".format(len(data_diff_dict[r"$\Delta p [\%]$"]))
          + "Length of data_nc: {}".format(len(data_nc["p"])))
except:
    pass

try:
    print("Length of diff fixed: {}\nLength of diff restart: {}".format(
        len((data_restart_fixed["p"] - data_nc["p"])*100 / data_nc["p"]),
        len((data_restart["p"] - data_nc["p"])*100 / data_nc["p"])
    ))
except:
    pass

try:
    tmp = data_restart_fixed[data_restart_fixed["trajectory"] == 0]
except:
    tmp = data_restart[data_restart["trajectory"] == 0]

for key in data_diff_dict.keys():
    print("{}: {}".format(key, np.shape(data_diff_dict[key])))
data_diff = pd.DataFrame.from_dict(data_diff_dict)
if get_time:
    end2 = timer.time()
    print("Creating data_diff in {} seconds".format(end2-end))
print("~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+")
print(data_nc["timestep"])
print(np.shape(data_nc["timestep"]))
print(data_restart["timestep"])
print(np.shape(data_restart["timestep"]))
d = data_nc["qg"]
print(d[0:2])
print(d[2000-10:2000+10])
print(d[2000*7-10:2000*7+10])
d = data_restart["qg"]
print(d[0:2])
print(d[2000-10:2000+10])
print(d[2000*7-10:2000*7+10])
data_diff.to_csv("diff_pd.csv")
plot_differences(r"$\Delta p [\%]$", data_diff)
if get_time:
    end = timer.time()
    print("One plot done in {} seconds".format(end-end2))
plot_differences(r"$\Delta T [\%]$", data_diff)
plot_differences(r"$\Delta q_r [\%]$", data_diff)
plot_differences(r"$\Delta q_v [\%]$", data_diff)
plot_differences(r"$\Delta q_c [\%]$", data_diff)
plot_differences(r"$\Delta q_i [\%]$", data_diff)
plot_differences(r"$\Delta q_s [\%]$", data_diff)
plot_differences(r"$\Delta q_g [\%]$", data_diff)
plot_differences(r"$\Delta w [\%]$", data_diff)
