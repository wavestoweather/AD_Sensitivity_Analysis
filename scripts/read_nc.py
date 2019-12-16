import cartopy.crs as ccrs
import cartopy
from matplotlib import pyplot as plt
import numpy as np
from netCDF4 import Dataset as netset
import pandas as pd
import progressbar as pb
import seaborn as sns
import sys
import xarray as xr
import progressbar


def plot_wcb(inp, outp, name):

    ds = xr.open_dataset(inp)
    df = ds.to_dataframe()

    n_rows = len(df["time"])
    n_trajectories = len(df["time"][0])
    dic = {
        "time": df["time"].values,
        "T": df["T"].values,
        "z": df["z"].values,
        "trajectory": [i % n_trajectories for i in range(n_rows)]}

    df2 = pd.DataFrame(dic)
    outpath = outp + name
    sns.set(style="darkgrid")

    n_timings = len(df2["time"].unique())
    n_rows = n_timings * n_trajectories
    df2 = df2.drop(df2.index[n_rows::])

    fig = plt.figure()
    slice_df = df2.pivot(index="time", columns="trajectory", values="T")
    s = sns.heatmap(slice_df, cmap="viridis", cbar_kws={"label": "T"})
    plt.tight_layout()
    plt.savefig(outpath + "_traject_T.png", dpi=300)
    plt.close()

    fig = plt.figure()
    slice_df = df2.pivot(index="time", columns="trajectory", values="z")
    s = sns.heatmap(slice_df, cmap="viridis", cbar_kws={"label": "z"})
    plt.tight_layout()
    plt.savefig(outpath + "_traject_z.png", dpi=300)
    plt.close()


def plot_foehn(inp, outp, name):

    outpath = outp + name
    ds = xr.open_dataset(inp)
    df = ds.to_dataframe()
    n_rows = len(df.index)
    n_times = len(ds["time"])
    n_trajectories = len(ds["id"])

    time = [i % n_times for i in range(n_rows)]
    df["time"] = time
    numpy_traj = ds["id"]
    trajectory = np.asarray([numpy_traj[i//n_times] for i in range(n_rows)])
    df["trajectory"] = trajectory

    df = df.dropna()
    plt.figure()
    slice_df = df.pivot(index="time", columns="trajectory", values="t")
    sns.heatmap(slice_df, cmap="viridis", cbar_kws={"label": "t"})
    plt.tight_layout()
    plt.savefig(outpath + "_traject_t.png", dpi=300)
    plt.close()

    plt.figure()
    slice_df = df.pivot(index="time", columns="trajectory", values="z")
    s = sns.heatmap(slice_df, cmap="viridis", cbar_kws={"label": "z"})
    plt.tight_layout()
    plt.savefig(outpath + "_traject_z.png", dpi=300)
    plt.close()

    plt.figure()
    slice_df = df.loc[df["trajectory"] == 0]
    sns.scatterplot("time", "qv", data=df)
    plt.tight_layout()
    plt.savefig(outpath + "_traject_qv.png", dpi=300)
    plt.close()

    plt.figure()
    slice_df = df.pivot(index="time", columns="trajectory", values="qv")
    s = sns.heatmap(slice_df, cmap="viridis", cbar_kws={"label": "qv"})
    plt.tight_layout()
    plt.savefig(outpath + "_traject_heat_qv.png", dpi=300)
    plt.close()

# def plot_map(xar, outpath):
def plot_map(data, outpath, var='qv'):
    extent = [data.variables['lon'][0], data.variables['lon'][-1],
              data.variables['lat'][0], data.variables['lat'][-1]]
    lon = data.variables['lon'][:]
    lat = data.variables['lat'][:]

    # lon = xar.lon.values
    # lat = xar.lat.values
    # extent = [lon[0], lon[len(xar.lon)-1],
    #           lat[0], lat[len(xar.lat)-1],]
    central_lon = (extent[1] - extent[0]) / 2.0 + extent[0]
    central_lat = (extent[3] - extent[2]) / 2.0 + extent[2]
    print("Plotting {}".format(var))
    # print("extent\n {}, \n\nlon\n {}, \n\nlat\n {}".format(extent, central_lon, central_lat))
    # for h in range(len(xar.height)):

    for h in pb.progressbar(range(len(data.variables['height']))):
        ax = plt.axes(projection=ccrs.Orthographic(central_lon, central_lat))
        ax.set_extent(extent)

        ax.coastlines(resolution="50m")
        # color for land and ocean
        ax.add_feature(cartopy.feature.OCEAN, zorder=2)
        ax.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='black')
        # ax.add_feature(cartopy.feature.LAKES, zorder=2, edgecolor='black')
        # ax.add_feature(cartopy.feature.RIVERS, zorder=2)

        # Add some data
        var_data = data.variables[var][0, h, :, :]
        # print(var_data[200:220, 200:220])
        # print(data[0:10, 0:10])
        # sys.exit()
        ax.contourf(lon, lat, var_data, alpha=0.9, zorder=3, transform=ccrs.PlateCarree())
        ax.gridlines(zorder=4)
        plt.savefig(outpath + "_" + var + "_" + "{:02d}".format(h) + "_map.png", dpi=300)
        ax.clear()
        plt.close()

# inp = ("/mnt/localscratch/data/project/m2_jgu-tapt/"
#        + "online_trajectories/foehn201305_case/foehn201305_cooling1.nc")
if len(sys.argv) == 1:
    inp = ("/data/project/m2_jgu-tapt/"
        + "online_trajectories/foehn201305_case/foehn201305_warming.nc")
    inp = ("/data/project/m2_jgu-tapt/online_trajectories/"
        + "wcb201609_vladiana/O_WCB_all_2016092" + "2_00.nc")
else:
    inp = sys.argv[1]

print("READING FROM {}\n".format(inp))
outp = "pics/"
name = "foehn"
# ds = xr.open_dataset(inp)
# print(ds)
if "NWP" in inp:
    # print("Keys:")
    # print(ds.keys())
    # print("qv shape")
    # qv = ds.qv
    # print(np.shape(qv))
    # print(type(qv[0, 0, 0, 0]))
    # print(qv[0, 0, 0, 0])
    # print(ds.qc[0, :, 0, 0])
    # qc = ds.qc[0, :, 0, 0]
    # print(qc)
    # print(qc.to_series())
    # print("Non zeros")
    # print(ds.where(ds.qnc > 0.0))
    # plot_map(ds, "pics/test")

    dataset = netset(inp)
    # d = dataset.variables["qv_2m"]
    # slice_d = d[0, 0, 100, 100]
    # print(slice_d)
    # exit()

    plot_map(dataset, "pics/test", "qr")
    plot_map(dataset, "pics/test", "temp")
    plot_map(dataset, "pics/test", "pres")
    plot_map(dataset, "pics/test", "qv")
    plot_map(dataset, "pics/test", "qc")
    plot_map(dataset, "pics/test", "w")
    # plot_map(dataset, "pics/test", "rain_gsp")

    # qv = xr.DataArray(qv, dims=[], coords={},)
    # print("\nqc:")
    # ds1 = ds
    # print(ds["qc"])
    # print(ds["qc"]["height"])
    # print(ds.sel(space="qc"))

else:
    ds = xr.open_dataset(inp)
    df = ds.to_dataframe()
    ntraj = 903
    steps = 7922
    print(df)

    print(df["QC"][480:490])
    for i in range(10):
        print(df["QC"][ntraj*i + 480 : ntraj*i + 490])
    exit()
    try:
        print(df["p"][0:40])
        print(df["p"][7920:7925])
        print(df["p"][steps*4:steps*4+5])
        print(df["p"][steps*5:steps*5+5])
    except:
        print(df["P"][0:40])
        print(df["P"][7920:7925])
        print(df["P"][steps*4:steps*4+5])
        print(df["P"][steps*5:steps*5+5])
    print(df["z"][steps*4:steps*4+5])
    print(df["z"][steps*5:steps*5+5])
    print("Find a specific value")
    try:
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(df.loc[np.abs(df["p"] - 488.4627) <= 0.001])
        print(df["qc"][0:40])
        print(df["qc"][7920:7925])
        print(df["qc"][steps*4:steps*4+5])
        print(df["qc"][steps*5:steps*5+5])
    except:
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(df.loc[np.abs(df["P"] - 488.4627) <= 0.001])
        print(df["QC"][0:40])
        print(df["QC"][7920:7925])
        print(df["QC"][steps*4:steps*4+5])
        print(df["QC"][steps*5:steps*5+5])
    # plot_foehn(inp, outp, name)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df[0:40])
    print("qv")
    try:
        print(df.loc[df["qv"] < 1e35])
        n_qv_small = len(df.loc[df["qv"] > 1e35].index)
    except:
        print(df.loc[df["QV"] < 1e35])
        n_qv_small = len(df.loc[df["QV"] > 1e35].index)
    n_total = len(df.index)
    print("{} / {}".format(n_qv_small, n_total))
    print("################")

    name = "wcb"
    inp = ("/mnt/localscratch/data/project/m2_jgu-tapt/"
        + "online_trajectories/wcb201609_vladiana/O_WCB_all_20160922_00.nc")
    # plot_wcb(inp, outp, name)
