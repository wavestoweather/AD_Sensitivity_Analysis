try:
    import latexify
except:
    import scripts.latexify as latexify

import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

# Plot large datasets
import datashader as ds

# see https://github.com/holoviz/datashader/issues/282
import holoviews.plotting.mpl
from holoviews.operation.datashader import datashade

from pylab import rcParams
import os

from itertools import repeat
from itertools import product
from progressbar import progressbar as pb

try:
    import dask_loader
except:
    import scripts.dask_loader as dask_loader

import dask.dataframe as pd
from bokeh.io import curdoc
from bokeh.palettes import Category20c
from timeit import default_timer as timer
import pandas
import xarray as xr
from glob import glob
from scipy.stats import chi2


class Deriv_dask:
    """
    Class that holds a dictionary of output parameters with its derivatives
    for every timestep as pandas.Dataframe. It can operate on them, such
    as deleting non-mapped regions and plotting the data.

    Parameters
    ----------
    data : Dic of pd.Dataframe
        Keys are output parameters, values are its derivatives
        for every timestep as pandas.Dataframe.
    clusternames : Dic of lsit of string
        Keys are output parameters where the value is a list of column names
        that hold a clustering assignment.
    threads : int
            Number of threads to use.
    font_dic : Dic
        Dictionary of labels and fontsize for plotting.
    backend : string
        Either "matplotlib" or "bokeh" for plotting.
    widget :
    colors : Colormap
        Colormap used for types of trajectories.
    cmap : Dic
        Dic of typenames of trajectories and colors.
    cache : pandas.Dataframe
        A dataframe after compute() was used that can be used as cache.
    cmap_particles : Dict of colors
        Mapping of output parameters (mixing ratios, saturation,
        particle numbers, latent heating and cooling) to colors.
    plots : list of plots
        Store plots in the list for further manipulation.
    """

    data = {}
    cluster_names = {}
    threads = 0
    font_dic = {}
    backend = "matplotlib"
    widget = None
    colors = None
    cmap = None
    cache = None
    cmap_particles = None
    plots = None

    def __init__(
        self,
        direc,
        parquet=True,
        netcdf=False,
        columns=None,
        backend="matplotlib",
        file_ending="*.nc_wcb",
    ):
        """
        Init class by loading the data from the given path.

        Parameters
        ----------
        direc : string
            A path to a directory wit a list of files to read.
        parquet : boolean
            If true: Load a series of preprocessed parquet files.
            This is no longer supported and will be deleted from future versions.
        netcdf : boolean
            If true: Load a series of preprocessed netcdf files.
        columns : list of strings
            Specify the columns to load.
        backend : String
            Either matplotlib or bokeh for plotting.
        file_ending : String
            File ending to use when loading netcdf files.
        """
        if parquet == True:
            print("Parquet is deprecated. Please use NetCDF.")
        if direc == "" and file_ending == "":
            # For creating a dummy/empty object
            self.data = None
        else:
            self.data = dask_loader.load_mult_derivates_direc_dic(
                direc, parquet, netcdf, columns, file_ending
            )

        self.cluster_names = {}
        # self.font_dic = {
        #     "title": 20,
        #     "labels": 20,
        #     "xticks": 12,
        #     "yticks": 16,
        #     "legend": 16
        # }
        self.font_dic = {
            "title": 10,
            "labels": 10,
            "xticks": 8,
            "yticks": 8,
            "legend": 10,
        }
        self.vertical_mark_fontsize = 10
        self.backend = backend
        self.plots = []
        if backend == "matplotlib":
            colors = plt.get_cmap("tab20c")
            self.cmap = {
                "Slantwise 600hPa 25. Quantile": matplotlib.colors.to_hex(
                    colors(0)[0:-1]
                ),
                "Slantwise 600hPa 50. Quantile": matplotlib.colors.to_hex(
                    colors(1)[0:-1]
                ),
                "Slantwise 600hPa 75. Quantile": matplotlib.colors.to_hex(
                    colors(2)[0:-1]
                ),
                "Slantwise 600hPa": matplotlib.colors.to_hex(colors(3)[0:-1]),
                "Slantwise 400hPa 25. Quantile": matplotlib.colors.to_hex(
                    colors(4)[0:-1]
                ),
                "Slantwise 400hPa 50. Quantile": matplotlib.colors.to_hex(
                    colors(5)[0:-1]
                ),
                "Slantwise 400hPa 75. Quantile": matplotlib.colors.to_hex(
                    colors(6)[0:-1]
                ),
                "Slantwise 400hPa": matplotlib.colors.to_hex(colors(7)[0:-1]),
                "Convective 600hPa 25. Quantile": matplotlib.colors.to_hex(
                    colors(8)[0:-1]
                ),
                "Convective 600hPa 50. Quantile": matplotlib.colors.to_hex(
                    colors(9)[0:-1]
                ),
                "Convective 600hPa 75. Quantile": matplotlib.colors.to_hex(
                    colors(10)[0:-1]
                ),
                "Convective 600hPa": matplotlib.colors.to_hex(colors(11)[0:-1]),
                "Convective 400hPa 25. Quantile": matplotlib.colors.to_hex(
                    colors(12)[0:-1]
                ),
                "Convective 400hPa 50. Quantile": matplotlib.colors.to_hex(
                    colors(13)[0:-1]
                ),
                "Convective 400hPa 75. Quantile": matplotlib.colors.to_hex(
                    colors(14)[0:-1]
                ),
                "Convective 400hPa": matplotlib.colors.to_hex(colors(15)[0:-1]),
            }
            self.colors = colors
            self.cmap_particles = {
                "QV": matplotlib.colors.to_hex([153 / 255, 204 / 255, 255 / 255]),
                "QC": matplotlib.colors.to_hex([102 / 255, 178 / 255, 255 / 255]),
                "QR": matplotlib.colors.to_hex([51 / 255, 153 / 255, 255 / 255]),
                "QI": matplotlib.colors.to_hex([215 / 255, 255 / 255, 255 / 255]),
                "QS": matplotlib.colors.to_hex([170 / 255, 255 / 255, 255 / 255]),
                "QG": matplotlib.colors.to_hex([0 / 255, 255 / 255, 255 / 255]),
                "QH": matplotlib.colors.to_hex([0 / 255, 204 / 255, 204 / 255]),
                "NC": matplotlib.colors.to_hex([0 / 255, 102 / 255, 255 / 255]),
                "NR": matplotlib.colors.to_hex([0 / 255, 0 / 255, 255 / 255]),
                "NI": matplotlib.colors.to_hex([102 / 255, 204 / 255, 255 / 255]),
                "NS": matplotlib.colors.to_hex([51 / 255, 204 / 255, 255 / 255]),
                "NG": matplotlib.colors.to_hex([0 / 255, 204 / 255, 255 / 255]),
                "NH": matplotlib.colors.to_hex([0 / 255, 153 / 255, 204 / 255]),
                "S": matplotlib.colors.to_hex([0 / 255, 204 / 255, 0 / 255]),
                "latent_heat": matplotlib.colors.to_hex([255 / 255, 51 / 255, 0 / 255]),
                "latent_cool": matplotlib.colors.to_hex(
                    [255 / 255, 153 / 255, 0 / 255]
                ),
            }
            self.cmap_types = {
                "artificial": matplotlib.colors.to_hex(colors(1)[0:-1]),
                "artificial (threshold)": matplotlib.colors.to_hex(colors(5)[0:-1]),
                "physical": matplotlib.colors.to_hex(colors(9)[0:-1]),
                "physical (high variability)": matplotlib.colors.to_hex(
                    colors(13)[0:-1]
                ),
                "1-moment": matplotlib.colors.to_hex(colors(17)[0:-1]),
            }
        else:
            colors = Category20c[20]
            self.cmap = {
                "Slantwise 600hPa 25. Quantile": colors[0],
                "Slantwise 600hPa 50. Quantile": colors[1],
                "Slantwise 600hPa 75. Quantile": colors[2],
                "Slantwise 600hPa": colors[3],
                "Slantwise 400hPa 25. Quantile": colors[4],
                "Slantwise 400hPa 50. Quantile": colors[5],
                "Slantwise 400hPa 75. Quantile": colors[6],
                "Slantwise 400hPa": colors[7],
                "Convective 600hPa 25. Quantile": colors[8],
                "Convective 600hPa 50. Quantile": colors[9],
                "Convective 600hPa 75. Quantile": colors[10],
                "Convective 600hPa": colors[11],
                "Convective 400hPa 25. Quantile": colors[12],
                "Convective 400hPa 50. Quantile": colors[13],
                "Convective 400hPa 75. Quantile": colors[14],
                "Convective 400hPa": colors[15],
            }
            self.colors = colors
            self.cmap_types = {
                "artificial": colors[1],
                "artificial (threshold)": colors[5],
                "physical": colors[9],
                "physical (high variability)": colors[13],
                "1-moment": colors[17],
            }

    def to_parquet(self, f_name, compression="snappy"):
        """
        Store the data as parquet files.
        This is deprecated and will be deleted in future version.

        Parameters
        ----------
        f_name : string
            Path to folder where to store the data.
        compression : string or dict
            From dask docs: Either a string like "snappy" or a dictionary mapping
            column names to compressors
            like {"name": "gzip", "values": "snappy"}.
            The default is "default", which uses the default compression for
            whichever engine is selected.
        """
        append = not os.listdir(f_name)
        append = not append
        pd.to_parquet(
            self.data,
            f_name,
            append=append,
            ignore_divisions=append,
            compression=compression,
        )

    def delete_not_mapped(self):
        """
        Delete all entries that are not within a mapped region, where
        mapped usually refers to timesteps where the WCB-criterion is
        satisfied.
        """
        self.data = self.data[self.data["MAP"] == True]

    def get_out_params(self):
        """
        Get all output parameters for which a dictionary of derivatives exists.

        Returns
        -------
        List of string
            Output parameters.
        """
        return self.data.keys()

    def cache_data(
        self,
        in_params,
        out_params,
        x_axis="step",
        y_axis=None,
        mapped=None,
        trajectories=None,
        frac=None,
        min_x=None,
        max_x=None,
        nth=None,
        compute=True,
    ):
        """
        Load some data into memory and store it as a pandas dataframe.
        May speed up some routines if the data fits into memory.

        Parameters
        ----------
        in_params : list of string
            Load the derivatives with respect to those in this list..
        out_params : list of string
            List of keys to plot the derivatives for.
        x_axis : string
            The column to use as x-axis. Can be either "step" or
            an output parameter or a derivative.
        y_axis : string
            y-axis for the upper plot. If none is given, use output parameter.
        mapped : string
            Column name which has to be true such as conv_400, slan_400,
            conv_600, slan_600.
        trajectories : List of int
            The index of the trajectories to plot. If None is given, all
            trajectories will be plotted.
        frac : float
            Sample a given fraction of rows. Deactivates "nth".
        nth : int
            Sample every nth entry. Works fast with "step" as x-axis and
            a given min_x and max_x value. If x_axis is any other than
            "step", an errorband triggered by "percentile" may not
            be plotted.
        min_x : float
            Minimum value for the x-axis.
        max_x : float
            Maximum value for the x-axis.
        compute : bool
            If true, store a pandas dataframe in self.cache. Otherwise dask
            dataframe (not entirely loaded) is stored.
        """
        t = timer()
        if frac is not None:
            df = self.data.sample(frac=frac, replace=False, random_state=42)
        elif nth is not None:
            if min_x is not None and max_x is not None and x_axis == "step":
                steps = np.arange(min_x, max_x, nth * 2.5)
            elif x_axis == "step":
                df_tmp = self.data.step.unique().compute()
                min_val = df_tmp.min()
                max_val = df_tmp.max()
                steps = np.arange(min_val, max_val, nth * 2.5)
            else:
                steps = self.data[x_axis].unique().compute().to_numpy()[::nth]

            df = self.data.loc[self.data[x_axis].isin(steps)]
        else:
            df = self.data
        if min_x is not None:
            df = df.loc[df[x_axis] >= min_x]
        if max_x is not None:
            df = df.loc[df[x_axis] <= max_x]

        if trajectories is not None:
            if isinstance(trajectories[0], int):
                df = df.loc[df["trajectory"].isin(trajectories)]
            else:
                df = df.loc[df["type"].isin(trajectories)]

        if mapped is not None:
            df = df.loc[df[mapped]]

        all_params = list(
            set(["trajectory", "type", "time"] + in_params + [x_axis] + out_params)
        )
        if "instance_id" in df:
            all_params.append("instance_id")
        if "Output Parameter" in df:
            df = df.loc[df["Output Parameter"].isin(out_params)]
            all_params.append("Output Parameter")

        if y_axis is not None and not y_axis in all_params:
            all_params.append(y_axis)
        if compute:
            self.cache = df[all_params].compute()
        else:
            self.cache = df[all_params]
        t2 = timer()
        print("Loading done in {} s".format(t2 - t))
        t = timer()

    def _recalc_ratios(
        self,
        ratio_df,
        ratio_type,
        ratio_window,
        in_params,
        x_axis="pressure",
        verbose=False,
    ):
        """
        Calculate the ratio of the sensitivity and return a dataframe
        with the new sensitivities.

        Parameters
        ----------
        ratio_df : Datframe
            Basically the cache of self or a subset of its cache
        ratio_type : String
            "vanilla": Use the derivative ratio in the file.
            "adjusted": Can be added to any to any type below where the sensitivity
            is adjusted to the parameter value such that perturbing this value by a
            certain percentage gives an approximation of the resulting error/difference
            for any given hydrometeor.
            "per_timestep": Use the highest derivative per timestep as denominator.
            "window": Use the highest derivative in the given window by min_x and max_x.
            "per_xaxis": Use the highest derivative per x_axis value. If x_axis is "timestep"
            it is the same as "per_timestep".
            "x_per_out_param": Replace 'x' with any other option than "vanilla". Use the highest
            derivative but per output parameter. (that *should* be the vanilla version)
            "x_weighted": Append this to any ratio type to use the inverse of the
            output parameter value as weight. Only works with "x_per_out_param".
        ratio_window : dic of list
            Overides ratio_type and switches to a derivative ratio calculation
            per output param if "ratio_type" has "per_out_param" in it.
            Calculate the derivative ratio in windows
            where the key is a column name and the list consists of values
            at which a new window starts, i.e. {"T": [235, 273]} results
            in three windows with T < 235K, 235 <= T < 273K and 237K <= T.
        in_params : list of string
            List of model parameters to consider.
        x_axis : string
            Use this column of the dataframe for ratio_type "per_xaxis".
        verbose : bool
            If true, print timing info.

        Returns
        -------
        Dataframe
            Same as ratio_df but with recalculated sensitivities.
        """

        if "weighted" in ratio_type and "per_out_param" in ratio_type:
            out_par = np.unique(ratio_df["Output Parameter"])[0]
            ratio_df.loc[:, in_params] = ratio_df[in_params].div(
                ratio_df[out_par], axis=0
            )
        if "adjusted" in ratio_type:
            for in_p in in_params:
                ratio_df[in_p] = ratio_df[in_p].mul(
                    latexify.in_params_numeric_value_dic[in_p] * 0.1
                )
        if ratio_window is not None:
            col = list(ratio_window.keys())[0]
            edges = np.sort(ratio_window[col])
            for i in range(len(edges) + 1):
                if i == 0:
                    df_edge = ratio_df.loc[ratio_df[col] < edges[i]]
                    max_val = (
                        df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                    )
                    if max_val == 0:
                        continue
                    ratio_df.loc[:, in_params] = ratio_df.apply(
                        lambda x: x[in_params] / max_val
                        if x[col] < edges[i]
                        else x[in_params],
                        axis=1,
                    )
                elif i == len(edges):
                    df_edge = ratio_df.loc[ratio_df[col] >= edges[i - 1]]
                    max_val = (
                        df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                    )
                    if max_val == 0:
                        continue
                    ratio_df.loc[:, in_params] = ratio_df.apply(
                        lambda x: x[in_params] / max_val
                        if x[col] >= edges[i - 1]
                        else x[in_params],
                        axis=1,
                    )
                else:
                    df_edge = ratio_df.loc[
                        (ratio_df[col] < edges[i]) & (ratio_df[col] >= edges[i - 1])
                    ]
                    max_val = (
                        df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                    )
                    if max_val == 0:
                        continue
                    ratio_df.loc[:, in_params] = ratio_df.apply(
                        lambda x: x[in_params] / max_val
                        if ((x[col] < edges[i]) and (x[col] >= edges[i - 1]))
                        else x[in_params],
                        axis=1,
                    )
        elif "per_timestep" in ratio_type:
            # Get series of max values over all timesteps (equals index)
            max_vals = ratio_df[in_params].apply(lambda x: np.max(np.abs(x)), axis=1)
            ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_vals, axis=0)
        elif "window" in ratio_type:
            max_val = ratio_df[in_params].apply(lambda x: np.max(np.abs(x))).max()
            ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_val)
        elif "per_xaxis" in ratio_type:
            ratio_df = ratio_df.set_index(x_axis)
            max_vals = (
                ratio_df.groupby(x_axis)[in_params]
                .apply(lambda x: np.max(np.abs(x)))
                .max(axis=1)
            )
            if (
                x_axis == "timestep"
                or x_axis == "step"
                or x_axis == "time_after_ascent"
            ):
                ratio_df.loc[:, in_params] = ratio_df[in_params].div(
                    max_vals, axis="index"
                )
            else:
                ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_vals)
        ratio_df.loc[:, in_params] = ratio_df[in_params].fillna(0)
        if verbose:
            t2 = timer()
            print("Recalculating ratios done in {} s".format(t2 - t))
        return ratio_df

    def get_important_sensitivities(
        self,
        in_params,
        out_params=None,
        n=None,
        ratio_type="vanilla",
        ratio_window=None,
        vertical_mark=None,
        verbose=False,
        x_axis="timestep",
    ):
        """
        Parameters
        ----------
        in_params : list of string
            List of all input parameters to be considered.
        out_params : list of string
            List of output parameters to get the most sensitive parameters for.
        n : int
            Get the n most important parameters.
        ratio_type : String
            "vanilla": Use the derivative ratio in the file that *should* use the
            highest derivative over all times for each output parameter as denominator.
            "per_timestep": Use the highest derivative per timestep as denominator.
            "window": Use the highest derivative in the given window by min_x and max_x.
            "per_xaxis": Use the highest derivative per x_axis value. If x_axis is "timestep"
            it is the same as "per_timestep".
            "x_per_out_param": Replace 'x' with any other option than "vanilla". Use the highest
            derivative but per output parameter. (that *should* be the vanilla version)
            "x_weighted": Append this to any ratio type to use the inverse of the
            output parameter value as weight. Only works with "x_per_out_param".
        ratio_window : dic of list
            Overides ratio_type and switches to a derivative ratio calculation
            per output param if "ratio_type" has "per_out_param" in it.
            Calculate the derivative ratio in windows
            where the key is a column name and the list consists of values
            at which a new window starts, i.e. {"T": [235, 273]} results
            in three windows with T < 235K, 235 <= T < 273K and 237K <= T.
        vertical_mark : dic of list
            A dictionary containing column names and values where a horizontal
            line should be created whenever the x_axis value intersects
            with the given value, i.e. {"T": [273, 235]} with x_axis in time
            marks all times, where a trajectory reached that temperature.
            Recommended to use with a single trajectory.

        Returns
        -------
        pandas.DataFrame
            A dataframe with the parameters and their respective,
            highest sensitivities. Column names are "Output Parameter",
            "description", "parameter", "value", "latex", "gradient_value"
        """
        if self.cache is None:
            print("You have to cache the data first using cache_data(..)")
            return None
        df = self.cache.copy()
        if verbose:
            t = timer()

        def recalc_ratios(ratio_df):
            if "weighted" in ratio_type and "per_out_param" in ratio_type:
                out_par = np.unique(ratio_df["Output Parameter"])[0]
                ratio_df.loc[:, in_params] = ratio_df[in_params].div(
                    ratio_df[out_par], axis=0
                )
            if ratio_window is not None:
                col = list(ratio_window.keys())[0]
                edges = np.sort(ratio_window[col])
                for i in range(len(edges) + 1):
                    if i == 0:
                        df_edge = ratio_df.loc[ratio_df[col] < edges[i]]
                        max_val = (
                            df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        )
                        if max_val == 0:
                            continue
                        ratio_df.loc[:, in_params] = ratio_df.apply(
                            lambda x: x[in_params] / max_val
                            if x[col] < edges[i]
                            else x[in_params],
                            axis=1,
                        )
                    elif i == len(edges):
                        df_edge = ratio_df.loc[ratio_df[col] >= edges[i - 1]]
                        max_val = (
                            df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        )
                        if max_val == 0:
                            continue
                        ratio_df.loc[:, in_params] = ratio_df.apply(
                            lambda x: x[in_params] / max_val
                            if x[col] >= edges[i - 1]
                            else x[in_params],
                            axis=1,
                        )
                    else:
                        df_edge = ratio_df.loc[
                            (ratio_df[col] < edges[i]) & (ratio_df[col] >= edges[i - 1])
                        ]
                        max_val = (
                            df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        )
                        if max_val == 0:
                            continue
                        ratio_df.loc[:, in_params] = ratio_df.apply(
                            lambda x: x[in_params] / max_val
                            if ((x[col] < edges[i]) and (x[col] >= edges[i - 1]))
                            else x[in_params],
                            axis=1,
                        )
            elif "per_timestep" in ratio_type:
                # Get series of max values over all timesteps (equals index)
                max_vals = ratio_df[in_params].apply(
                    lambda x: np.max(np.abs(x)), axis=1
                )
                ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_vals, axis=0)
            elif "window" in ratio_type:
                max_val = ratio_df[in_params].apply(lambda x: np.max(np.abs(x))).max()
                ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_val)
            elif "per_xaxis" in ratio_type:
                ratio_df = ratio_df.set_index(x_axis)
                max_vals = (
                    ratio_df.groupby(x_axis)[in_params]
                    .apply(lambda x: np.max(np.abs(x)))
                    .max(axis=1)
                )
                if (
                    x_axis == "timestep"
                    or x_axis == "step"
                    or x_axis == "time_after_ascent"
                ):
                    ratio_df.loc[:, in_params] = ratio_df[in_params].div(
                        max_vals, axis="index"
                    )
                else:
                    ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_vals)
            ratio_df.loc[:, in_params] = ratio_df[in_params].fillna(0)
            if verbose:
                t2 = timer()
                print("Recalculating ratios done in {} s".format(t2 - t))
            return ratio_df

        if out_params is None:
            out_params = df["Output Parameter"].unique()
        if "per_out_param" not in ratio_type:
            df = recalc_ratios(df)
            if verbose:
                t = timer()

        dic_sorted = {
            "Output Parameter": [],
            "description": [],
            "parameter": [],
            "value": [],
            "latex": [],
            "gradient_value": [],
        }
        for out_par in out_params:
            df_tmp_out = df.loc[df["Output Parameter"] == out_par]
            if "per_out_param" in ratio_type:
                df_tmp_out = recalc_ratios(df_tmp_out)
                if verbose:
                    t = timer()

            sorted_tuples = []
            for in_p in in_params:
                value = np.abs(df_tmp_out[in_p].min())
                max_val = np.abs(df_tmp_out[in_p].max())
                if max_val > value:
                    value = max_val
                if value != 0 and not np.isnan(value):
                    sorted_tuples.append((in_p, value))
            sorted_tuples.sort(key=lambda tup: tup[1])
            if verbose:
                t2 = timer()
                print("Sorting done in {} s".format(t2 - t), flush=True)
            for i, tup in enumerate(sorted_tuples[::-1]):
                if n is not None and i == n:
                    break
                p, g = tup
                dic_sorted["Output Parameter"].append(out_par)
                dic_sorted["parameter"].append(p)
                dic_sorted["gradient_value"].append(g)
                try:
                    dic_sorted["value"].append(latexify.in_params_value_dic[p])
                except:
                    dic_sorted["value"].append("no value provided")
                try:
                    dic_sorted["description"].append(latexify.in_params_descr_dic[p])
                except:
                    dic_sorted["description"].append("no description provided")
                dic_sorted["latex"].append(latexify.parse_word(p))
        df_sorted = pandas.DataFrame.from_dict(dic_sorted)
        return df_sorted

    def plot_two_ds(
        self,
        in_params,
        out_params,
        x_axis="timestep",
        y_axis=None,
        twin_axis=None,
        decimals=-3,
        mapped=None,
        trajectories=None,
        scatter=False,
        n_plots=None,
        percentile=None,
        frac=None,
        min_x=None,
        max_x=None,
        nth=None,
        hist=[False, False],
        hexbin=[False, False],
        bins=50,
        log=[False, False],
        sort=True,
        scatter_deriv=False,
        line_deriv=False,
        compute=False,
        use_cache=False,
        errorband=False,
        plot_path="pics/",
        prefix=None,
        fig_type="svg",
        datashade=True,
        by=None,
        alpha=[1, 1],
        rolling_avg=20,
        rolling_avg_par=20,
        max_per_deriv=10,
        width=1959,
        height=1224,
        ratio_type="vanilla",
        ratio_window=None,
        vertical_mark=None,
        plot_singles=False,
        xticks=10,
        param_method="",
        deriv_method="",
        **kwargs,
    ):
        """
        Plot two plots in two rows. At the top: Output parameter.
        At the bottom: Derivative with respect to that output parameter.
        Another plot is created for every gradient of different order.
        If multiple parameters are given, another plot is created for each
        output parameter. This function is mostly intended for multiple
        trajectories.
        For small fractions being plotted, initialization takes long
        (ie 85 seconds) vs plotting (roughly 0.45 seconds per plot).

        Parameters
        ----------
        in_params : list of string
            Plot the derivatives with respect to those in this list.
        out_params : list of string
            List of keys to plot the derivatives for.
        x_axis : string
            The column to use as x-axis. Can be either "timestep" or
            an output parameter or a derivative.
        y_axis : string
            y-axis for the upper plot. If none is given, use output parameter.
        twin_axis : string
            Plot a second y-axis for both plots using the column given
            by "twin_axis". Works only with matplotlib atm. Does not
            take "by" into account and might look cluttered with multiple
            elements from "by".
        decimals : int
            For rounding the right y-axis, if provided by "twin_axis".
            From numpy: Number of decimal places to round to (default: 0).
            If decimals is negative, it specifies the number of
            positions to the left of the decimal point.
        mapped : string
            Column name which has to be true such as conv_400, slan_400,
            conv_600, slan_600.
        trajectories : List of int or list of string
            If int: the index of the trajectories to plot from column "trajectory".
            If string: the type name of the trajectories to plot from column "type".
            If None is given, all trajectories will be plotted.
        scatter : boolean
            Plot a scatter plot or a line plot.
        n_plots : int
            Plot only that many plots. If None is given, plot every possible
            plot.
        percentile : list of int
            Plot the given percentiles along with the spread if
            errorband is not given. Not implemented yet.
        frac : float
            Sample a given fraction of rows. Deactivates "nth".
        min_x : float
            Minimum value for the x-axis.
        max_x : float
            Maximum value for the x-axis.
        nth : int
            Sample every nth entry. Works fast with "timestep" as x-axis and
            a given min_x and max_x value. If x_axis is any other than
            "timestep", an errorband triggered by "percentile" may not
            be plotted.
        hist : list of bools
            Plot a histogram on the side of the top [0] or bottom [1]
            (not implemented) graph.
        hexbin : list of bools
            Plot values as hexbin for the top [0] or bottom [1] graph.
        bins : int
            Number of bins to use on hexbin and histogram plots.
        log : list of bools
            Plot y-axis as log-plot for the top [0] or bottom [1] graph.
        sort : Bool
            If True, sort the derivatives and plot only those within the same
            magnitude in one plot. If False, plot every derivative.
        scatter_deriv : boolean
            Plot the derivative plot with scatterplot. If scatter_deriv and
            line_deriv are set to false, use the same as provided via scatter
            and percentile.
        line_deriv : boolean
            Plot the derivative plot with lineplot. If scatter_deriv and
            line_deriv are set to false, use the same as provided via scatter
            and percentile.
        compute : bool
            Load the data into memory and create a pandas.DataFrame that will
            vbe stored to cache.
        use_cache : bool
            Copy the pandas.DataFrame from the cache and work with that.
            Overrides "compute".
        errorband : Bool
            Plot an errorband (spread) using the standard deviation and min/max values.
            This method is not properly tested!
            TODO: Debugging.
        plot_path : string
            Path to the folder to save the image.
        prefix : string
            Prefix to add to the filename.
        fig_type : string
            Figure type for saving the image (only for bokeh).
        datashade : bool
            Wether to use datashade for plotting the data.
        by : String
            String for groupby.
            Can be either "trajectory" or "type" (recommended).
        alpha : List of floats
            Alpha values for top [0] and bottom [1] graph.
        rolling_avg : int
            Number of timesteps to use for a rolling average of the derivatives.
            If None, no rolling average is being calculated.
            If "by" is set, the rolling average will be calculated
            per instance in the column from "by".
        rolling_avg_par : int
            Same as rollling_avg but for output parameters.
        max_per_deriv : int
            Maximum number of derivatives per plot.
        width : int
            The width of the plot (bokeh) or the relative width of the plot
            (matplotlib).
        height : int
            The height of the plot (bokeh) or the relative height of the plot
            (matplotlib).
        ratio_type : String
            "vanilla": Use the derivative ratio in the file that *should* use the
            highest derivative over all times for each output parameter as denominator.
            "per_timestep": Use the highest derivative per timestep as denominator.
            "window": Use the highest derivative in the given window by min_x and max_x.
            "per_xaxis": Use the highest derivative per x_axis value. If x_axis is "timestep"
            it is the same as "per_timestep".
            "x_per_out_param": Replace 'x' with any other option than "vanilla". Use the highest
            derivative but per output parameter. (that *should* be the vanilla version)
            "x_weighted": Append this to any ratio type to use the inverse of the
            output parameter value as weight. Only works with "x_per_out_param".
        ratio_window : dic of list
            Overides ratio_type and switches to a derivative ratio calculation
            per output param if "ratio_type" has "per_out_param" in it.
            Calculate the derivative ratio in windows
            where the key is a column name and the list consists of values
            at which a new window starts, i.e. {"T": [235, 273]} results
            in three windows with T < 235K, 235 <= T < 273K and 237K <= T.
        vertical_mark : dic of list
            A dictionary containing column names and values where a horizontal
            line should be created whenever the x_axis value intersects
            with the given value, i.e. {"T": [273, 235]} with x_axis in time
            marks all times, where a trajectory reached that temperature.
            Recommended to use with a single trajectory.
        plot_singles : bool
            Plot every plot individually as well.
        xticks : int
            Number of ticks on all x-axis.
        param_method : string
            Alternative way to define which plot to use for the param plot.
            Options are "scatter", "hexbin", "errorband", "percentile".
        deriv_method : string
            Alternative way to define which plot to use for the derivative plot.
            Options are "scatter", "hexbin", "heatmap", "heatmap_error"
        kwargs : dict
            Keyword arguments are passed down matplotlib.
        """
        import hvplot.dask  # adds hvplot method to dask objects
        import hvplot.pandas
        import hvplot
        from holoviews import opts
        import holoviews as hv
        from holoviews.operation import histogram as hv_histo
        import pandas
        import dask.array as da

        hv.extension(self.backend)
        deriv_col_name = "Derivative Ratio"
        matplotlib.rcParams["axes.formatter.limits"] = (-2, 2)

        # Adjust scatter size according to height
        scatter_size = int(height / 200)
        fig_inches = width / 300
        if width < height:
            fig_inches = height / 300
        aspect = width / height

        if not use_cache:
            self.cache_data(
                in_params,
                out_params,
                x_axis,
                y_axis,
                mapped,
                trajectories,
                frac,
                min_x,
                max_x,
                nth,
                compute,
            )
            if compute:
                df = self.cache.copy()
            else:
                df = self.cache
        else:
            df = self.cache.copy()
        t = timer()

        if trajectories is not None:
            if isinstance(trajectories[0], int):
                df = df.loc[df["trajectory"].isin(trajectories)]
            else:
                df = df.loc[df["type"].isin(trajectories)]

        if rolling_avg is not None:
            # For every output parameter, calculate the rolling average
            # for all input parameters at once
            for out_par in out_params:
                df_tmp = df.loc[df["Output Parameter"] == out_par]
                if by is None:
                    # TODO: is that really rolling over consecutive timesteps? What if derivative gets 0?
                    series = (
                        df_tmp[in_params].rolling(rolling_avg, min_periods=1).mean()
                    )
                    series.index = df.loc[df["Output Parameter"] == out_par].index
                    df.update(series)
                else:
                    types = df_tmp[by].unique()
                    for ty in types:
                        df_tmp2 = df_tmp.loc[df_tmp[by] == ty]
                        tmp_par = ["dcloud_a_vel", x_axis]
                        df_tmp2 = df_tmp.loc[df_tmp[by] == ty]
                        series = (
                            df_tmp2[in_params]
                            .rolling(rolling_avg, min_periods=1)
                            .mean()
                        )
                        series.index = df.loc[
                            (df["Output Parameter"] == out_par) & (df[by] == ty)
                        ].index
                        df.update(series)
            t2 = timer()
            print(
                "Calculating rolling average for derivatives done in {} s".format(
                    t2 - t
                )
            )
            t = timer()

        if rolling_avg_par is not None:
            for out_par in out_params:
                df_tmp = df.loc[df["Output Parameter"] == out_par]
                if by is None:
                    series = (
                        df_tmp[out_par].rolling(rolling_avg_par, min_periods=1).mean()
                    )
                    series.index = df.loc[df["Output Parameter"] == out_par].index
                    df.update(series)
                else:
                    types = df_tmp[by].unique()
                    for ty in types:
                        df_tmp2 = df_tmp.loc[df_tmp[by] == ty]
                        series = (
                            df_tmp2[out_par]
                            .rolling(rolling_avg_par, min_periods=1)
                            .mean()
                        )
                        series.index = df.loc[
                            (df["Output Parameter"] == out_par) & (df[by] == ty)
                        ].index
                        df.update(series)
            t2 = timer()
            print("Calculating rolling average for outputs done in {} s".format(t2 - t))
            t = timer()

        def recalc_ratios(ratio_df):
            if "weighted" in ratio_type and "per_out_param" in ratio_type:
                out_par = np.unique(ratio_df["Output Parameter"])[0]
                ratio_df.loc[:, in_params] = ratio_df[in_params].div(
                    ratio_df[out_par], axis=0
                )
            if ratio_window is not None:
                col = list(ratio_window.keys())[0]
                edges = np.sort(ratio_window[col])
                for i in range(len(edges) + 1):
                    if i == 0:
                        df_edge = ratio_df.loc[ratio_df[col] < edges[i]]
                        max_val = (
                            df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        )
                        if max_val == 0:
                            continue
                        ratio_df.loc[:, in_params] = ratio_df.apply(
                            lambda x: x[in_params] / max_val
                            if x[col] < edges[i]
                            else x[in_params],
                            axis=1,
                        )
                    elif i == len(edges):
                        df_edge = ratio_df.loc[ratio_df[col] >= edges[i - 1]]
                        max_val = (
                            df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        )
                        if max_val == 0:
                            continue
                        ratio_df.loc[:, in_params] = ratio_df.apply(
                            lambda x: x[in_params] / max_val
                            if x[col] >= edges[i - 1]
                            else x[in_params],
                            axis=1,
                        )
                    else:
                        df_edge = ratio_df.loc[
                            (ratio_df[col] < edges[i]) & (ratio_df[col] >= edges[i - 1])
                        ]
                        max_val = (
                            df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        )
                        if max_val == 0:
                            continue
                        ratio_df.loc[:, in_params] = ratio_df.apply(
                            lambda x: x[in_params] / max_val
                            if ((x[col] < edges[i]) and (x[col] >= edges[i - 1]))
                            else x[in_params],
                            axis=1,
                        )
            elif "per_timestep" in ratio_type:
                # Get series of max values over all timesteps (equals index)
                max_vals = ratio_df[in_params].apply(
                    lambda x: np.max(np.abs(x)), axis=1
                )
                ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_vals, axis=0)
            elif "window" in ratio_type:
                max_val = ratio_df[in_params].apply(lambda x: np.max(np.abs(x))).max()
                ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_val)
            elif "per_xaxis" in ratio_type:
                ratio_df = ratio_df.set_index(x_axis)
                max_vals = (
                    ratio_df.groupby(x_axis)[in_params]
                    .apply(lambda x: np.max(np.abs(x)))
                    .max(axis=1)
                )
                if (
                    x_axis == "timestep"
                    or x_axis == "step"
                    or x_axis == "time_after_ascent"
                ):
                    ratio_df.loc[:, in_params] = ratio_df[in_params].div(
                        max_vals, axis="index"
                    )
                else:
                    ratio_df.loc[:, in_params] = ratio_df[in_params].div(max_vals)
            ratio_df.loc[:, in_params] = ratio_df[in_params].fillna(0)
            t2 = timer()
            print("Recalculating ratios done in {} s".format(t2 - t))
            return ratio_df

        if not "per_out_param" in ratio_type:
            if by is not None:
                for ty in np.unique(df["type"]):
                    df.loc[df["type"] == ty] = recalc_ratios(df.loc[df["type"] == ty])
            else:
                df = recalc_ratios(df)
            t = timer()

        deriv_title = "Deriv. Ratio"
        if ratio_window is not None:
            deriv_title += " within Vertical Lines"
        elif "per_timestep" in ratio_type:
            deriv_title += " per Time Step"
        elif "per_xaxis" in ratio_type:
            deriv_title += " per X-Axis"

        set_yaxis = False
        if y_axis is None:
            set_yaxis = True
        for out_par in out_params:
            t = timer()
            if set_yaxis:
                y_axis = out_par
            df_tmp_out = df.loc[df["Output Parameter"] == out_par]
            if "per_out_param" in ratio_type:
                if by is not None:
                    for ty in np.unique(df_tmp_out["type"]):
                        df_tmp_out.loc[df_tmp_out["type"] == ty] = recalc_ratios(
                            df_tmp_out.loc[df_tmp_out["type"] == ty]
                        )
                else:
                    df_tmp_out = recalc_ratios(df_tmp_out)
                t = timer()

            # This is for the matplotlib backend
            # For bokeh, see https://github.com/holoviz/holoviews/issues/396
            twin = None
            if twin_axis is not None and self.backend == "matplotlib":
                lower_y = np.around(df_tmp_out[twin_axis].min(), decimals=decimals)
                upper_y = np.around(df_tmp_out[twin_axis].max(), decimals=decimals)
                tick_size = (upper_y - lower_y) / 10

                def twinx(plot, element):
                    ax = plot.handles["axis"]
                    twinax = ax.twinx()
                    twinax.set_ylim(
                        (lower_y - tick_size * 0.75, upper_y + tick_size * 0.75)
                    )
                    twinax.set_yticks(
                        np.arange(lower_y, upper_y + tick_size * 0.75, tick_size)
                    )
                    twinax.set_ylabel(
                        latexify.parse_word(twin_axis)
                        + " "
                        + latexify.get_unit(twin_axis, brackets=True)
                    )
                    plot.handles["axis"] = twinax

                if log[1]:
                    df_tmp_out[twin_axis] = da.log(da.fabs(df_tmp_out[twin_axis]))
                twin = (
                    df_tmp_out.hvplot.scatter(
                        x=x_axis,
                        y=twin_axis,
                        color="gray",
                        datashade=datashade,
                        legend=False,
                    )
                    .opts(initial_hooks=[twinx], apply_ranges=False)
                    .opts(opts.Scatter(s=scatter_size))
                )

            # Sort the derivatives
            if sort:
                sorted_tuples = []

                for in_p in in_params:
                    if log[1]:
                        value = np.log(np.abs(df_tmp_out[in_p].min()))
                        max_val = np.log(np.abs(df_tmp_out[in_p].max()))
                        if max_val > value:
                            value = max_val
                        if not np.isnan(value) and not np.isinf(value):
                            sorted_tuples.append((in_p, value))
                    else:
                        value = np.abs(df_tmp_out[in_p].min())
                        max_val = np.abs(df_tmp_out[in_p].max())
                        if max_val > value:
                            value = max_val
                        if value != 0 and not np.isnan(value):
                            sorted_tuples.append((in_p, value))
                sorted_tuples.sort(key=lambda tup: tup[1])
                t2 = timer()
                print("Sorting done in {} s".format(t2 - t), flush=True)

            def plot_helper(df, in_params, prefix, **kwargs):
                deriv_col_name = "Derivative Ratio"
                # following https://holoviz.org/tutorial/Composing_Plots.html
                t = timer()
                add_cols = []
                if vertical_mark is not None:
                    add_cols = list(vertical_mark.keys())

                if by is not None:
                    df_tmp = df[in_params + [x_axis, by] + add_cols]
                    df_tmp = df_tmp.melt(
                        [x_axis, by] + add_cols,
                        var_name="Derivatives",
                        value_name=deriv_col_name,
                    )
                else:
                    df_tmp = df[in_params + [x_axis] + add_cols]
                    df_tmp = df_tmp.melt(
                        [x_axis] + add_cols,
                        var_name="Derivatives",
                        value_name=deriv_col_name,
                    )
                t2 = timer()
                print("Melting done in {} s".format(t2 - t), flush=True)

                if log[1]:
                    df_tmp[deriv_col_name] = da.log(da.fabs(df_tmp[deriv_col_name]))
                    # Remove zero entries (-inf)
                    df_tmp = df_tmp[~da.isinf(df_tmp[deriv_col_name])]
                    df_tmp.rename(
                        columns={deriv_col_name: "Log Derivative Ratio"}, inplace=True
                    )
                    deriv_col_name = "Log Derivative Ratio"
                    t2 = timer()
                    print("Log done in {} s".format(t2 - t))
                    t = timer()
                df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)

                if "percentile" in globals() or param_method == "percentile":
                    if percentile not in globals():
                        percentile = [0.25, 0.5, 0.75]
                    if y_axis == x_axis:
                        print(
                            "x-axis and y-axis are the same. Cannot plot that with percentiles!"
                        )
                        return
                    else:
                        df_group = df[[x_axis, y_axis] + add_cols]
                    if log[0]:
                        # Apply should be more expensive
                        df_group[y_axis] = da.log(da.fabs(df_group[y_axis]))
                        # Remove zero entries (-inf)
                        df_group = df_group[~da.isinf(df_group[y_axis])]

                    # Group for min, max and percentiles
                    funcs = [np.min, np.max] + [
                        lambda x, perc=perc: np.percentile(x, perc, axis=0)
                        for perc in percentile
                    ]
                    df_min_max = df_group.groupby(x_axis).agg(funcs)[y_axis]

                    # Rename the columns
                    p_list = []
                    p_dic = {}
                    for i, perc in enumerate(percentile):
                        p_list.append("{}. Percentile".format(perc))
                        p_dic["<lambda_{}>".format(i)] = p_list[-1]
                    p_dic["amin"] = "Min"
                    p_dic["amax"] = "Max"
                    df_min_max = df_min_max.rename(columns=p_dic)

                    # Plot
                    param_plot = df_min_max.hvplot.area(
                        x=x_axis,
                        y="Min",
                        y2="Max",
                        alpha=0.5,
                        ylabel=latexify.parse_word(y_axis),
                        title="Values of {}".format(latexify.parse_word(y_axis)),
                        label="Spread",
                        color="grey",
                    ) * df_min_max.hvplot.line(
                        x=x_axis, y=p_list, ylabel=latexify.parse_word(y_axis), **kwargs
                    )
                elif param_method == "errorband" or errorband:
                    df_group = df[[x_axis, y_axis, "trajectory"] + add_cols]
                    if log[0]:
                        # Apply should be more expensive
                        df_group[y_axis] = da.log(da.fabs(df_group[y_axis]))
                    df_min_max = (
                        df_group.groupby([x_axis, "trajectory"])[y_axis]
                        .mean()
                        .groupby(x_axis)
                        .agg([np.min, np.max])
                    )
                    df_std = (
                        df_group.groupby([x_axis, "trajectory"])[y_axis]
                        .mean()
                        .groupby(x_axis)
                        .agg(
                            [
                                lambda x: -1 * np.std(x) + np.mean(x),
                                lambda x: np.std(x) + np.mean(x),
                            ]
                        )
                    )
                    df_mean = df_group.groupby(x_axis)[y_axis].mean()

                    param_plot = (
                        df_min_max.hvplot.area(
                            x=x_axis,
                            y="amin",
                            y2="amax",
                            alpha=0.5,
                            value_label=latexify.parse_word(y_axis),
                            label="Spread",
                        )
                        * df_mean.hvplot()
                        * df_std.hvplot.area(
                            x=x_axis,
                            y="<lambda_0>",
                            y2="<lambda_1>",
                            alpha=0.3,
                            value_label=latexify.parse_word(y_axis),
                            label="sd",
                        )
                    )
                elif hexbin[0] or param_method == "hexbin":
                    if y_axis == x_axis:
                        df_group = df[[x_axis, "trajectory"] + add_cols]
                    else:
                        df_group = df[[x_axis, y_axis, "trajectory"] + add_cols]
                    if log[0]:
                        # Apply should be more expensive
                        df_group[y_axis] = da.log(da.fabs(df_group[y_axis]))
                    param_plot = df_group.hvplot.hexbin(
                        x=x_axis,
                        y=y_axis,
                        title="Values of {}".format(latexify.parse_word(y_axis)),
                        label=None,
                        clabel="Count",
                        cmap="viridis",
                        logz=True,
                        gridsize=100,
                    ).options(ylabel=latexify.parse_word(y_axis))
                elif by is not None:
                    if y_axis == x_axis:
                        df_group = df[[x_axis, by] + add_cols]
                    else:
                        df_group = df[[x_axis, y_axis, by] + add_cols]
                    if log[0]:
                        # Apply should be more expensive
                        df_group[y_axis] = da.log(da.fabs(df_group[y_axis]))
                    types = df_group[df_group[y_axis] != 0][by].unique()
                    types = np.sort(types[::-1])
                    if not datashade:
                        cmap_values = []
                        for ty in types:
                            cmap_values.append(self.cmap[ty])
                    else:
                        cmap = {}
                        for ty in types:
                            cmap[ty] = self.cmap[ty]
                    if twin is not None:
                        cmap_values.append("gray")
                        types = np.append(types, latexify.parse_word(twin_axis))
                    if not datashade:
                        if self.backend == "matplotlib":
                            overlay = hv.NdOverlay(
                                {
                                    types[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                                        opts.Scatter(
                                            s=scatter_size * 8, color=cmap_values[i]
                                        )
                                    )
                                    for i in range(len(types))
                                }
                            )

                            param_plot = (
                                df_group.hvplot.scatter(
                                    x=x_axis,
                                    y=y_axis,
                                    by=by,
                                    title="Values of {}".format(
                                        latexify.parse_word(y_axis)
                                    ),
                                    color=cmap_values,
                                    label=None,
                                    datashade=datashade,
                                    alpha=alpha[0],
                                    legend=False,
                                )
                                .opts(opts.Scatter(s=scatter_size))
                                .opts(aspect=aspect)
                                .options(ylabel=latexify.parse_word(y_axis))
                                * overlay
                            )

                            if hist[0]:
                                xhist = df_group.hvplot.hist(
                                    y=x_axis, bins=bins, height=125
                                )
                                yhist = df_group.hvplot.hist(
                                    y=y_axis, bins=bins, width=125
                                )
                                param_plot = param_plot << yhist << xhist
                        else:
                            param_plot = (
                                df_group.hvplot.scatter(
                                    x=x_axis,
                                    y=y_axis,
                                    by=by,
                                    title="Values of {}".format(
                                        latexify.parse_word(y_axis)
                                    ),
                                    color=cmap_values,
                                    label=None,
                                    datashade=datashade,
                                    alpha=alpha[0],
                                )
                                .opts(opts.Scatter(size=scatter_size))
                                .options(ylabel=latexify.parse_word(y_axis))
                            )
                    else:
                        param_plot = (
                            df_group.hvplot.scatter(
                                x=x_axis,
                                y=y_axis,
                                by=by,
                                title="Values of {}".format(
                                    latexify.parse_word(y_axis)
                                ),
                                cmap=cmap,
                                label=None,
                                datashade=datashade,
                            )
                            .opts(aspect=width / height)
                            .options(ylabel=latexify.parse_word(y_axis))
                        )
                        if self.backend == "bokeh":
                            points = hv.Points(
                                (
                                    [np.NaN for i in range(len(list(cmap.keys())))],
                                    [np.NaN for i in range(len(list(cmap.keys())))],
                                    list(cmap.keys()),
                                ),
                                vdims=by,
                            )
                            legend = points.options(
                                color=by, cmap=cmap, show_legend=True
                            )
                            param_plot = (legend * param_plot).opts(aspect=aspect)
                else:
                    if y_axis == x_axis:
                        df_group = df[[x_axis, "trajectory"] + add_cols]
                    else:
                        df_group = df[[x_axis, y_axis, "trajectory"] + add_cols]
                    if log[0]:
                        # Apply should be more expensive
                        df_group[y_axis] = da.log(da.fabs(df_group[y_axis]))
                    param_plot = df_group.hvplot.scatter(
                        x=x_axis,
                        y=y_axis,
                        title="Values of {}".format(latexify.parse_word(y_axis)),
                        label=None,
                        datashade=datashade,
                    ).options(ylabel=latexify.parse_word(y_axis))

                t2 = timer()
                print("Setting up upper plot done in {} s".format(t2 - t))
                t = timer()

                if hexbin[1] or deriv_method == "hexbin":
                    deriv_plot = df_tmp.hvplot.hexbin(
                        x=x_axis,
                        y=deriv_col_name,
                        by="Derivatives",
                        title=deriv_title + " for {}".format(out_par),
                        label=None,
                        logz=True,
                        gridsize=100,
                        clabel="Count",
                        cmap="viridis",
                        colorbar=True,
                    )
                elif "heatmap" in deriv_method:
                    if self.backend == "matplotlib":
                        timely_data = (
                            df_tmp.drop("type", axis=1)
                            .groupby(x_axis)
                            .sum()
                            .reset_index()
                        )
                        plot_data = pandas.melt(
                            timely_data,
                            id_vars=[x_axis],
                            var_name="Derivatives",
                            value_name=deriv_col_name,
                        )
                        heat_plot = hv.HeatMap(
                            plot_data, label=deriv_title + " for {}".format(out_par)
                        )
                        # Add an errorbar
                        if "error" in deriv_method:
                            aggr = hv.Dataset(heat_plot).aggregate(
                                x_axis, np.mean, np.std
                            )
                            aggr_plot = hv.ErrorBars(aggr) * hv.Curve(aggr)
                            deriv_plot = (heat_plot + aggr_plot).cols(1)
                        else:
                            deriv_plot = heat_plot
                    else:
                        print(f"heatmap not tested for backend {self.backend}")
                else:
                    colors = plt.get_cmap("tab10")
                    all_derives = df_tmp["Derivatives"].unique()
                    if len(all_derives) > 10 and len(all_derives) < 21:
                        colors = plt.get_cmap("tab20")
                    elif len(all_derives) > 20:
                        colors = plt.get_cmap("gist_ncar")
                    all_derives = np.sort(all_derives[::-1])
                    if not datashade:
                        cmap_values = []
                        for i in range(len(all_derives)):
                            cmap_values.append(
                                matplotlib.colors.to_hex(colors(i)[0:-1])
                            )
                        if twin is not None:
                            cmap_values.append("gray")
                            all_derives = np.append(
                                all_derives, latexify.parse_word(twin_axis)
                            )
                        if self.backend == "matplotlib":
                            overlay = hv.NdOverlay(
                                {
                                    derivative: hv.Scatter((np.NaN, np.NaN)).opts(
                                        opts.Scatter(
                                            s=scatter_size * 8, color=cmap_values[i]
                                        )
                                    )
                                    for i, derivative in enumerate(all_derives)
                                }
                            )
                            deriv_plot = (
                                df_tmp.hvplot.scatter(
                                    x=x_axis,
                                    y=deriv_col_name,
                                    by="Derivatives",
                                    title=deriv_title + " for {}".format(out_par),
                                    label=None,
                                    datashade=datashade,
                                    alpha=alpha[1],
                                    legend=False,
                                    cmap=cmap_values,
                                )
                                .opts(opts.Scatter(s=scatter_size))
                                .opts(aspect=aspect)
                            )
                            deriv_plot = deriv_plot * overlay
                        else:
                            deriv_plot = df_tmp.hvplot.scatter(
                                x=x_axis,
                                y=deriv_col_name,
                                by="Derivatives",
                                title=deriv_title + " for {}".format(out_par),
                                label=None,
                                datashade=datashade,
                                alpha=alpha[1],
                                cmap=cmap_values,
                            ).opts(opts.Scatter(size=scatter_size))

                    else:
                        cmap = {}
                        for i, d in enumerate(all_derives):
                            cmap[d] = matplotlib.colors.to_hex(colors(i)[0:-1])
                        deriv_plot = df_tmp.hvplot.scatter(
                            x=x_axis,
                            y=deriv_col_name,
                            by="Derivatives",
                            title=deriv_title + " for {}".format(out_par),
                            label=None,
                            datashade=datashade,
                            cmap=cmap,
                        ).opts(aspect=width / height)

                        if self.backend == "bokeh":
                            points_deriv = hv.Points(
                                (
                                    [np.NaN for i in range(len(all_derives))],
                                    [np.NaN for i in range(len(all_derives))],
                                    list(cmap.keys()),
                                ),
                                vdims="Derivatives",
                            )
                            legend_deriv = points_deriv.options(
                                color="Derivatives", cmap=cmap, show_legend=True
                            )
                            deriv_plot = (legend_deriv * deriv_plot).opts(
                                aspect=width / height
                            )

                t2 = timer()
                print("Setting up lower plot done in {} s".format(t2 - t))
                t = timer()

                # Create vertical lines if needed
                if vertical_mark is not None:
                    marks = None
                    if min_x is None:
                        this_min_x = df_group[x_axis].min()
                    else:
                        this_min_x = min_x
                    if max_x is None:
                        this_max_x = df_group[x_axis].max()
                    else:
                        this_max_x = max_x
                    min_y = df_group[y_axis].min()
                    max_y = df_group[y_axis].max()
                    del_x = (this_max_x - this_min_x) * 0.1
                    del_y = max_y - min_y
                    if min_y == max_y:
                        min_y = max_y - 1
                        max_y = max_y + 1
                        del_y = 2

                    min_y_deriv = df_tmp[deriv_col_name].min()
                    max_y_deriv = df_tmp[deriv_col_name].max()
                    del_y_deriv = max_y_deriv - min_y_deriv
                    if min_y_deriv == max_y_deriv:
                        min_y_deriv = max_y_deriv - 1
                        max_y_deriv = max_y_deriv + 1
                        del_y_deriv = 2

                    if by is not None:
                        df_tmp_group = df_group.groupby(by)
                    else:
                        df_tmp_group = df_group

                    for col in vertical_mark:
                        scale = 1
                        for v in vertical_mark[col]:
                            # Works well as long as dataframe is smaller than ~ 30k elements
                            df_sort = df_tmp_group.apply(
                                lambda x: x.iloc[np.argmin(np.abs(x[col] - v))]
                            )
                            col_values = np.sort(df_sort[col].values)
                            if scale == 8:
                                scale = 1
                            for index, row in df_sort.iterrows():
                                if np.abs(row[col] - v) > 1.0:
                                    break
                                # The one for the plot above
                                text = hv.Text(
                                    row[x_axis] + del_x,
                                    max_y - del_y * 0.1 * scale,
                                    f"{col}={v:.2f}{latexify.get_unit(col)}",
                                    fontsize=self.vertical_mark_fontsize,
                                )
                                mark = (
                                    hv.VLine(
                                        x=row[x_axis], label=col + "=" + str(v)
                                    ).opts(
                                        color="black",
                                        ylim=(min_y - del_y * 0.1, max_y + del_y * 0.1),
                                    )
                                    * text
                                )
                                param_plot = param_plot * mark
                                # For the derivatives
                                text = hv.Text(
                                    row[x_axis] + del_x,
                                    max_y_deriv - del_y_deriv * 0.1 * scale,
                                    f"{col}={v:.2f}{latexify.get_unit(col)}",
                                    fontsize=self.vertical_mark_fontsize,
                                )
                                mark = (
                                    hv.VLine(
                                        x=row[x_axis], label=col + "=" + str(v)
                                    ).opts(
                                        color="black",
                                        ylim=(
                                            min_y_deriv - del_y_deriv * 0.1,
                                            max_y_deriv + del_y_deriv * 0.1,
                                        ),
                                    )
                                    * text
                                )
                                deriv_plot = deriv_plot * mark
                                scale += 1
                    t2 = timer()
                    print("Creating marks done in {} s".format(t2 - t))
                    t = timer()

                # Add second y-axis if given
                if twin is not None:
                    param_plot = param_plot * twin
                    if (
                        "heatmap" not in deriv_method
                        or "hexbin" != deriv_method
                        or not hexbin[1]
                    ):
                        deriv_plot = deriv_plot * twin

                layout = param_plot + deriv_plot

                opts_arg = {}  # Currently empty. Maybe useful in further iterations

                scatter_kwargs = {}
                area_kwargs = {}

                for k in kwargs:
                    scatter_kwargs[k] = kwargs[k]

                if self.backend == "matplotlib":
                    area_kwargs["edgecolor"] = None
                    area_kwargs["color"] = "black"
                    for k in kwargs:
                        area_kwargs[k] = kwargs[k]

                layout_kwargs = {}
                if self.backend == "bokeh":
                    layout_kwargs["width"] = width
                    layout_kwargs["height"] = height
                else:
                    layout_kwargs["fig_inches"] = fig_inches

                # Matplotlib uses a horrible colormap as default...
                curve_kwargs = kwargs.copy()
                if "percentile" in globals():
                    if self.backend == "matplotlib":
                        if len(percentile) <= 10:
                            curve_kwargs["color"] = hv.Cycle("tab10")
                        elif len(percentile) <= 20:
                            curve_kwargs["color"] = hv.Cycle("tab20")
                    else:
                        if len(percentile) <= 10:
                            curve_kwargs["color"] = hv.Cycle("Category10")
                        elif len(percentile) <= 20:
                            curve_kwargs["color"] = hv.Cycle("Category20")

                if errorband or param_method == "errorband":
                    both_plots = layout.opts(
                        opts.Area(
                            xticks=xticks,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **area_kwargs,
                        ),
                        opts.Scatter(
                            xticks=20,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **scatter_kwargs,
                        ),
                        opts.Layout(**layout_kwargs),
                    ).cols(1)
                elif "percentile" in globals() or param_method == "percentile":
                    both_plots = layout.opts(
                        opts.Area(
                            xticks=xticks,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **area_kwargs,
                        ),
                        opts.Scatter(
                            xticks=20,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **scatter_kwargs,
                        ),
                        opts.Curve(**curve_kwargs),
                        opts.Layout(**layout_kwargs),
                    ).cols(1)
                else:
                    if (not hexbin[0] and not hexbin[1]) and (
                        not "hexbin" == param_method and not "hexbin" == deriv_method
                    ):
                        if "heatmap" in deriv_method:
                            both_plots = layout.opts(
                                opts.Scatter(
                                    xticks=xticks,
                                    xaxis="bottom",
                                    fontsize=self.font_dic,
                                    show_grid=True,
                                    show_legend=True,
                                    xlabel=latexify.parse_word(x_axis),
                                    **scatter_kwargs,
                                ),
                                opts.HeatMap(
                                    # width=width,
                                    logz=True,
                                    xaxis=None,
                                    colorbar=True,
                                    cmap="viridis",
                                ),
                                opts.Layout(**layout_kwargs),
                            ).cols(1)
                        else:
                            both_plots = layout.opts(
                                opts.Scatter(
                                    xticks=xticks,
                                    xaxis="bottom",
                                    fontsize=self.font_dic,
                                    show_grid=True,
                                    show_legend=True,
                                    xlabel=latexify.parse_word(x_axis),
                                    **scatter_kwargs,
                                ),
                                opts.Layout(**layout_kwargs),
                            ).cols(1)
                    elif (hexbin[0] and hexbin[1]) or (
                        "hexbin" == deriv_method and "hexbin" == param_method
                    ):
                        both_plots = layout.opts(
                            opts.HexTiles(**opts_arg), opts.Layout(**layout_kwargs)
                        ).cols(1)
                    else:
                        both_plots = layout.opts(
                            opts.Scatter(
                                xticks=xticks,
                                xaxis="bottom",
                                fontsize=self.font_dic,
                                show_grid=True,
                                show_legend=True,
                                xlabel=latexify.parse_word(x_axis),
                                **scatter_kwargs,
                            ),
                            opts.HexTiles(**opts_arg),
                            opts.Layout(**layout_kwargs),
                        ).cols(1)

                if self.backend == "matplotlib":
                    both_plots = both_plots.opts(sublabel_format="", tight=True)
                if hist[0]:
                    print("Histogram for the upper plot is not yet supported.")
                elif by is not None:
                    param_opts = param_plot.opts(xaxis="bare")
                else:
                    param_opts = param_plot.opts(xaxis="bare")

                t2 = timer()
                print("Setting options done in {} s".format(t2 - t))
                t = timer()

                renderer = hv.Store.renderers[self.backend].instance(fig="png", dpi=300)

                i = 0
                if prefix is None:
                    prefix = "plt_line_"
                    if scatter:
                        prefix = "plt_scatter_"
                save = (
                    plot_path
                    + prefix
                    + x_axis
                    + "_"
                    + out_par
                    + "_"
                    + "{:03d}".format(i)
                )
                while os.path.isfile(save + ".png"):
                    i = i + 1
                    save = (
                        plot_path
                        + prefix
                        + x_axis
                        + "_"
                        + out_par
                        + "_"
                        + "{:03d}".format(i)
                    )
                t2 = timer()
                print("Everything else done in {} s".format(t2 - t), flush=True)
                if self.backend == "bokeh":
                    print("Plotting", flush=True)
                    hvplot.show(both_plots)
                    self.plots.append(both_plots)
                    t2 = timer()
                    print("Plotting done in {}s".format(t2 - t), flush=True)
                else:
                    filetype = ".png"
                    if datashade:
                        filetype = ".html"
                    self.plots.append(both_plots)
                    print("Saving to " + save + filetype, flush=True)
                    renderer.save(both_plots, save)
                    display_file = save
                    if plot_singles:
                        for j, pl in enumerate(both_plots):
                            which_plot = "_outParam_"
                            if j > 0:
                                which_plot = "_deriv_"
                            i = 0
                            save = (
                                plot_path
                                + prefix
                                + x_axis
                                + "_"
                                + out_par
                                + which_plot
                                + "{:03d}".format(i)
                            )
                            while os.path.isfile(save + ".png"):
                                i = i + 1
                                save = (
                                    plot_path
                                    + prefix
                                    + x_axis
                                    + "_"
                                    + out_par
                                    + which_plot
                                    + "{:03d}".format(i)
                                )
                            pl = pl.opts(
                                opts.Scatter(
                                    xticks=xticks,
                                    xaxis="bottom",
                                    fontsize=self.font_dic,
                                    show_grid=True,
                                    show_legend=True,
                                    **scatter_kwargs,
                                )
                            ).opts(fig_inches=fig_inches, aspect=aspect, xaxis="bottom")
                            renderer.save(pl, save)
                    t2 = timer()
                    try:
                        from IPython.display import Image, display

                        display(Image(display_file + filetype, width=width))
                    except:
                        pass

                    print("Saving done in {}s".format(t2 - t), flush=True)

            if sort:
                i = 0
                if n_plots is None:
                    n_plots = 9999999
                while len(sorted_tuples) > 0 and i < n_plots:
                    p, v = sorted_tuples.pop()
                    in_params_2 = [p]
                    v_max = v
                    if log[1]:
                        while (
                            len(sorted_tuples) > 0
                            and np.abs(v_max) - np.abs(sorted_tuples[-1][1]) < 15
                        ):
                            p, v = sorted_tuples.pop()
                            in_params_2.append(p)
                            if len(in_params_2) == max_per_deriv:
                                break
                    else:
                        while (
                            len(sorted_tuples) > 0
                            and np.abs(sorted_tuples[-1][1] / v_max) > 0.1
                        ):
                            p, v = sorted_tuples.pop()
                            in_params_2.append(p)
                            if len(in_params_2) == max_per_deriv:
                                break
                    print(f"Plotting for {in_params_2}")
                    plot_helper(
                        df_tmp_out, in_params=in_params_2, prefix=prefix, **kwargs
                    )
                    i += 1
            else:
                # Plot for every input parameter
                for in_param in in_params:
                    plot_helper(
                        df_tmp_out,
                        in_params=[in_param],
                        prefix=prefix + "_" + in_param,
                        **kwargs,
                    )

    def plot_grid_one_param(
        self,
        out_param,
        y_axis,
        x_axis="step",
        twin_axis=None,
        by=None,
        hue="type",
        col_wrap=4,
        trajectories=None,
        width=1959,
        height=1224,
        log=[False, False],
        vertical_mark=None,
        cross_mark=None,
        datashade=False,
        prefix=None,
        alpha=1,
        plot_path="pics/",
        yticks=10,
        xticks=10,
        decimals=-3,
        rolling_avg=20,
        kind="scatter",
        plot_singles=False,
        s=None,
        formatter_limits=None,
        **kwargs,
    ):
        """
        Plot a grid for comparing multiple output parameters or
        multiple derivatives across one output parameter.
        Only uses cached and computed data!

        Parameters
        ----------
        out_param : string
            Use this output parameter as reference, i.e. if y-axis consists
            of derivatives, then those are sensitivities towards the output
            parameter. If the y-axis are output parameters, then out_param
            can be any parameter.
        y_axis : List of string
            The name of the column for the y-axis for each plot.
        x_axis : string
            The column for the x-axis.
        twin_axis : string
            Plot a second y-axis for every plot using the column given
            by "twin_axis". Works only with matplotlib atm. Does not
            take "by" into account and might look cluttered with multiple
            elements from "by".
        by : String
            String for groupby.
            Can be either "trajectory" or "type" (recommended).
        hue : string
            Define the column to colorize multiple outputs in one plot.
        col_wrap : int
            Number of plots per column
        trajectories : List of int or list of string
            If int: the index of the trajectories to plot from column "trajectory".
            If string: the type name of the trajectories to plot from column "type".
            If None is given, all trajectories will be plotted.
        width : int
            The width of the plot (bokeh) or the relative width of the plot
            (matplotlib).
        height : int
            The height of the plot (bokeh) or the relative height of the plot
            (matplotlib).
        log : List of bools
            Plot the left y-axis (log[0]) with logarithmic scale or the right
            y-axis (log[1]) with logarithmic scale.
        vertical_mark : dic of list
            A dictionary containing column names and values where a horizontal
            line should be created whenever the x_axis value intersects
            with the given value, i.e. {"T": [273, 235]} with x_axis in time
            marks all times, where a trajectory reached that temperature.
            Recommended to use with a single trajectory.
        cross_mark : dic of list
            A dictionary containing column names and values where a cross
            should be created whenever the x_axis value intersects
            with the given value, i.e. {"T": [273, 235]} with x_axis in time
            marks all times, where a trajectory reached that temperature.
            Recommended to use with a single trajectory.
        datashade : bool
            Use datashade to plot results.
            Not yet implemented/tested.
        prefix : string
            Prefix to add to the filename.
        alpha : float
            Alpha value to use on the plot aligned to the left y-axis.
        plot_path : string
            Path to the folder to save the image.
        yticks : int
            Number of ticks on all y-axis.
        xticks : int
            Number of ticks on all x-axis.
        decimals : int
            For rounding the right y-axis, if provided by "twin_axis".
            From numpy: Number of decimal places to round to (default: 0).
            If decimals is negative, it specifies the number of
            positions to the left of the decimal point.
        rolling_avg : int
            Number of consecutive rows for a rolling average
            for the left y-axis.
        kind : string
            "scatter" for a scatter plot,
            "hexbin" for hexagonal bins with logarithmic count
            and viridis colormap, else defaults to line plot.
        plot_singles : bool
            Plot every plot as an individual plot in addition to the grid of
            plots.
        s : int
            Can be used to adjust the size of the scatter dots.
        formatter_limits : tuple of ints
            Lower and upper limits for formatting x- and y-axis.
        kwargs : dict
            Keyword arguments are passed down matplotlib.
        """
        import hvplot.pandas
        from holoviews import opts
        import holoviews as hv
        import dask.array as da
        import math
        import pandas
        from holoviews.operation.datashader import datashade as dsshade

        dsshade.dynamic = False
        fontscale = width / 2200
        if formatter_limits is not None:
            matplotlib.rcParams["axes.formatter.limits"] = formatter_limits

        df = self.cache
        if df is None:
            print("Cached data is None. Make sure to cache your data using")
            print("self.cache_data()\nABORTING")
            return
        hv.extension(self.backend)

        aspect = width / height
        if s is None:
            scatter_size = int(height / 200)
        else:
            scatter_size = s

        fig_inches = width / 300
        if width < height:
            fig_inches = height / 300
        layout_kwargs = {}
        if self.backend == "bokeh":
            layout_kwargs["width"] = width
            layout_kwargs["height"] = height
        else:
            layout_kwargs["fig_inches"] = fig_inches

        if kind == "hexbin":
            datashade = False

        plot_list = []
        if "Output Parameter" in df:
            df_tmp_out = df.loc[df["Output Parameter"] == out_param]
        else:
            df_tmp_out = df
        # df_tmp_out["qi"] = da.fabs(df_tmp_out["qi"]) # sign error?
        if trajectories is not None:
            if isinstance(trajectories[0], int):
                df_tmp_out = df_tmp_out.loc[df_tmp_out["trajectory"].isin(trajectories)]
            else:
                df_tmp_out = df_tmp_out.loc[df_tmp_out["type"].isin(trajectories)]
        # This is for the matplotlib backend
        # For bokeh, see https://github.com/holoviz/holoviews/issues/396
        twin = None
        if twin_axis is not None and self.backend == "matplotlib":
            lower_y = np.around(df_tmp_out[twin_axis].min(), decimals=decimals)
            upper_y = np.around(df_tmp_out[twin_axis].max(), decimals=decimals)
            tick_size = (upper_y - lower_y) / yticks

            def twinx(plot, element):
                ax = plot.handles["axis"]
                twinax = ax.twinx()
                twinax.set_ylim(
                    (lower_y - tick_size * 0.75, upper_y + tick_size * 0.75)
                )
                twinax.set_yticks(
                    np.arange(lower_y, upper_y + tick_size * 0.75, tick_size)
                )
                twinax.set_ylabel(
                    latexify.parse_word(twin_axis)
                    + " "
                    + latexify.get_unit(twin_axis, brackets=True)
                )
                plot.handles["axis"] = twinax

            if log[1]:
                df_tmp_out[twin_axis] = da.log(da.fabs(df_tmp_out[twin_axis]))
            twin = (
                df_tmp_out.hvplot.scatter(
                    x=x_axis,
                    y=twin_axis,
                    color="gray",
                    datashade=datashade,
                    legend=False,
                )
                .opts(initial_hooks=[twinx], apply_ranges=False)
                .opts(opts.Scatter(s=scatter_size))
            )

        add_cols = []
        if vertical_mark is not None:
            add_cols = list(vertical_mark.keys())
            if x_axis in add_cols:
                add_cols.remove(x_axis)

        for y in y_axis:
            marks = None
            if by is not None:
                add_col_tmp = None
                if y in add_cols:
                    add_cols.remove(y)
                    add_col_tmp = y

                if y == x_axis:
                    df_group = df_tmp_out[[x_axis, by] + add_cols]
                else:
                    df_group = df_tmp_out[[x_axis, y, by] + add_cols]
                # print(df_group)
                # print("############################")
                # print(df_group[by])
                types = df_group[by].unique()
                types = np.sort(types[::-1])
                if rolling_avg is not None:
                    for ty in types:
                        df_roll = df_group.loc[df_group[by] == ty]
                        series = df_roll[y].rolling(rolling_avg, min_periods=1).mean()
                        series.index = df_group.loc[df_group[by] == ty].index
                        df_group.update(series)
                if log[0]:
                    # Apply should be more expensive
                    df_group[y] = da.log(da.fabs(df_group[y]))

                min_x = df_group[x_axis].min()
                max_x = df_group[x_axis].max()
                min_y = df_group[y].min()
                max_y = df_group[y].max()

                del_x = (max_x - min_x) * 0.1
                del_y = max_y - min_y

                # datashade cannot handle only zero valued plots
                if min_y == max_y and datashade:
                    continue
                if datashade and df_group.empty:
                    continue

                if min_y == max_y:
                    min_y = max_y - 1
                    max_y = max_y + 1
                    del_y = 2

                if kind == "scatter":
                    plot_func = df_group.hvplot.scatter
                    if datashade:
                        options = opts.Scatter(show_grid=True)
                    else:
                        options = opts.Scatter(s=scatter_size, show_grid=True)
                elif kind == "hexbin":

                    def agg(x):
                        # print(f"{y}: {np.shape(x)} -- {np.sum(x)}, {np.log(np.sum(x))}")
                        return np.log(np.sum(x))

                    plot_func = df_group.hvplot.hexbin
                    if self.backend == "matplotlib":
                        options = opts.HexTiles(
                            colorbar=False,
                            cmap="viridis",
                            logz=False,
                            show_grid=True,
                            gridsize=int(np.ceil(width / 30)),
                            aggregator=agg,
                        )
                    else:
                        options = opts.HexTiles(
                            colorbar=False,
                            logz=False,
                            show_grid=True,
                            tools=["hover"],
                            gridsize=int(np.ceil(width / 30)),
                            aggregator=agg,
                        )
                else:
                    plot_func = df_group.hvplot.line
                    if self.backend == "bokeh":
                        options = opts.Curve(
                            line_width=self.vertical_mark_fontsize - 2, show_grid=True
                        )
                    else:
                        options = opts.Curve(
                            linewidth=self.vertical_mark_fontsize - 4, show_grid=True
                        )

                if vertical_mark is not None:
                    if by is not None:
                        df_tmp = df_group.groupby(by)
                    else:
                        df_tmp = df_group

                    for col in vertical_mark:
                        for v in vertical_mark[col]:
                            if y == col:
                                df_sort = df_tmp.apply(
                                    lambda x: x.iloc[np.argmin(np.abs(x[col] - v))]
                                )
                            else:
                                # Works well as long as dataframe is smaller than ~ 30k elements
                                df_sort = df_tmp.apply(
                                    lambda x: x.iloc[np.argmin(np.abs(x[col] - v))]
                                )
                            col_values = np.sort(df_sort[col].values)

                            for index, row in df_sort.iterrows():
                                row_col = row[col]

                                if np.abs(row_col - v) > 1.0:
                                    break
                                text = hv.Text(
                                    row[x_axis] + del_x,
                                    max_y - del_y * 0.1,
                                    col + "=" + str(v) + latexify.get_unit(col),
                                    fontsize=self.vertical_mark_fontsize,
                                )
                                if marks is None:
                                    marks = (
                                        hv.VLine(
                                            x=row[x_axis], label=col + "=" + str(v)
                                        ).opts(
                                            color="black",
                                            ylim=(
                                                min_y - del_y * 0.1,
                                                max_y + del_y * 0.1,
                                            ),
                                        )
                                        * text
                                    )
                                else:
                                    marks = (
                                        marks
                                        * hv.VLine(
                                            x=row[x_axis], label=col + "=" + str(v)
                                        ).opts(
                                            color="black",
                                            ylim=(
                                                min_y - del_y * 0.1,
                                                max_y + del_y * 0.1,
                                            ),
                                        )
                                        * text
                                    )
                if cross_mark is not None:

                    def add_marks(df_tmp, marks):
                        for col in cross_mark:
                            df_mark = df_tmp.loc[df_tmp[col].isin(cross_mark[col])]
                            if marks is None:
                                marks = (
                                    df_mark.hvplot.scatter(
                                        x=x_axis,
                                        y=y,
                                        color="black",
                                        marker="x",
                                        datashade=datashade,
                                        legend=False,
                                    )
                                    .opts(apply_ranges=False)
                                    .opts(opts.Scatter(s=scatter_size))
                                )
                            else:
                                marks = (
                                    marks
                                    * df_mark.hvplot.scatter(
                                        x=x_axis,
                                        y=y,
                                        color="black",
                                        marker="x",
                                        datashade=datashade,
                                        legend=False,
                                    )
                                    .opts(apply_ranges=False)
                                    .opts(opts.Scatter(s=scatter_size))
                                )
                        return marks

                    if by is not None:
                        for byby in df_group[by].unique():
                            df_tmp = df_group.loc[df_group[by] == byby]
                            marks = add_marks(df_tmp, marks)
                    else:
                        df_tmp = df_group
                        marks = add_marks(df_tmp, marks)

                if not datashade:
                    cmap_values = []
                    for ty in types:
                        cmap_values.append(self.cmap[ty])
                    if twin is not None:
                        cmap_values.append("gray")
                        types = np.append(types, latexify.parse_word(twin_axis))
                    lim_multiplier = 0.1
                    if twin is not None:
                        lim_multiplier = 0.7
                    if self.backend == "matplotlib":
                        if len(plot_list) == 0 and kind != "hexbin":
                            overlay = hv.NdOverlay(
                                {
                                    types[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                                        opts.Scatter(
                                            s=scatter_size * 8, color=cmap_values[i]
                                        )
                                    )
                                    for i in range(len(types))
                                }
                            )
                            plot_list.append(
                                plot_func(
                                    x=x_axis,
                                    y=y,
                                    by=by,
                                    color=cmap_values,
                                    label=None,
                                    ylabel=latexify.parse_word(y),
                                    xlabel=latexify.parse_word(x_axis),
                                    datashade=datashade,
                                    alpha=alpha,
                                    legend=False,
                                )
                                .opts(options)
                                .opts(
                                    title="Values of of {}".format(
                                        latexify.parse_word(y)
                                    ),
                                    aspect=aspect,
                                    ylim=(
                                        min_y - del_y * lim_multiplier,
                                        max_y + del_y * lim_multiplier,
                                    ),
                                    xlim=(min_x, max_x),
                                    yticks=yticks,
                                    xticks=xticks,
                                    fontscale=fontscale,
                                )
                                * overlay
                            )
                        else:
                            plot_list.append(
                                plot_func(
                                    x=x_axis,
                                    y=y,
                                    by=by,
                                    color=cmap_values,
                                    label=None,
                                    ylabel=latexify.parse_word(y),
                                    xlabel=latexify.parse_word(x_axis),
                                    datashade=datashade,
                                    alpha=alpha,
                                    legend=False,
                                )
                                .opts(options)
                                .opts(
                                    title="Values of of {}".format(
                                        latexify.parse_word(y)
                                    ),
                                    aspect=aspect,
                                    ylim=(
                                        min_y - del_y * lim_multiplier,
                                        max_y + del_y * lim_multiplier,
                                    ),
                                    xlim=(min_x, max_x),
                                    yticks=yticks,
                                    xticks=xticks,
                                    fontscale=fontscale,
                                )
                            )
                        if marks is not None:
                            plot_list[-1] = plot_list[-1] * marks
                        if twin is not None:
                            plot_list[-1] = plot_list[-1] * twin
                    else:
                        plot_list.append(
                            df_group[df_group[y] != 0]
                            .hvplot.scatter(
                                x=x_axis,
                                y=y,
                                by=by,
                                title="Values of of {}".format(latexify.parse_word(y)),
                                label=None,
                                datashade=datashade,
                                alpha=alpha,
                            )
                            .opts(
                                opts.Scatter(size=scatter_size, fontscale=fontscale),
                                # **layout_kwargs,
                            )
                        )

                else:
                    lim_multiplier = 0.1
                    if twin is not None:
                        lim_multiplier = 0.7
                    if kind == "hexbin":
                        cmap = "viridis"
                    if len(plot_list) == 0 and kind != "hexbin":
                        cmap = {}
                        for ty in types:
                            cmap[ty] = self.cmap[ty]

                        overlay = hv.NdOverlay(
                            {
                                types[i]: hv.Scatter((np.NaN, np.NaN)).opts(
                                    opts.Scatter(
                                        s=scatter_size * 8, color=cmap[types[i]]
                                    )
                                )
                                for i in range(len(types))
                            }
                        )

                        plot_list.append(
                            plot_func(
                                x=x_axis,
                                y=y,
                                by=by,
                                cmap=cmap,
                                label=None,
                                ylabel=latexify.parse_word(y),
                                xlabel=latexify.parse_word(x_axis),
                                datashade=datashade,
                                width=width,
                                height=height,
                                dynamic=False,
                                ylim=(
                                    min_y - del_y * lim_multiplier,
                                    max_y + del_y * lim_multiplier,
                                ),
                                xlim=(min_x, max_x),
                                dynspread=True,
                                yticks=yticks,
                                xticks=xticks,
                                aspect=aspect,
                            )
                            .opts(options)
                            .opts(
                                show_grid=True,
                                fontscale=fontscale,
                                title="Values of of {}".format(latexify.parse_word(y)),
                            )
                            * overlay
                        )

                    else:
                        plot_list.append(
                            plot_func(
                                x=x_axis,
                                y=y,
                                by=by,
                                cmap=cmap,
                                label=None,
                                ylabel=latexify.parse_word(y),
                                xlabel=latexify.parse_word(x_axis),
                                datashade=datashade,
                                width=width,
                                height=height,
                                dynamic=False,
                                ylim=(
                                    min_y - del_y * lim_multiplier,
                                    max_y + del_y * lim_multiplier,
                                ),
                                xlim=(min_x, max_x),
                                dynspread=True,
                                yticks=yticks,
                                xticks=xticks,
                                aspect=aspect,
                            )
                            .opts(options)
                            .opts(
                                show_grid=True,
                                fontscale=fontscale,
                                title="Values of of {}".format(latexify.parse_word(y)),
                            )
                        )

                    if marks is not None:
                        plot_list[-1] = plot_list[-1] * marks
                    if twin is not None:
                        plot_list[-1] = plot_list[-1] * twin
                if add_col_tmp is not None:
                    add_cols.append(add_col_tmp)
            else:
                if y == x_axis:
                    df_group = df_tmp_out[[x_axis, "trajectory"]]
                else:
                    df_group = df_tmp_out[[x_axis, y, "trajectory"]]

                if rolling_avg is not None:
                    series = df_group[y].rolling(rolling_avg, min_periods=1).mean()
                    df_group.update(series)
                if kind == "scatter":
                    plot_func = df_group.hvplot.scatter
                    if datashade:
                        options = opts.Scatter(show_grid=True)
                    else:
                        options = opts.Scatter(s=scatter_size, show_grid=True)
                elif kind == "hexbin":

                    def agg(x):
                        # print(f"{y}: {np.shape(x)} -- {np.sum(x)}, {np.log(np.sum(x))}")
                        return np.log(np.sum(x))

                    plot_func = df_group.hvplot.hexbin
                    if self.backend == "matplotlib":
                        options = opts.HexTiles(
                            colorbar=False,
                            cmap="viridis",
                            logz=False,
                            show_grid=True,
                            gridsize=int(np.ceil(width / 30)),
                            aggregator=agg,
                        )
                    else:
                        options = opts.HexTiles(
                            colorbar=False,
                            logz=False,
                            show_grid=True,
                            tools=["hover"],
                            gridsize=int(np.ceil(width / 30)),
                            aggregator=agg,
                        )
                else:
                    plot_func = df_group.hvplot.line
                    if self.backend == "bokeh":
                        options = opts.Curve(
                            line_width=self.vertical_mark_fontsize - 2, show_grid=True
                        )
                    else:
                        options = opts.Curve(
                            linewidth=self.vertical_mark_fontsize - 4, show_grid=True
                        )

                if log:
                    # Apply should be more expensive
                    df_group[y] = da.log(da.fabs(df_group[y]))
                plot_list.append(
                    plot_func(
                        x=x_axis,
                        y=y,
                        title="Values of of {}".format(latexify.parse_word(y)),
                        label=None,
                        datashade=datashade,
                    ).opts(**layout_kwargs)
                )
        if len(plot_list) > 1:
            all_plots = plot_list[0] + plot_list[1]
            for i in range(len(plot_list) - 2):
                all_plots = all_plots + plot_list[i + 2]
        else:
            all_plots = plot_list[0]
            plot_singles = False
        # layout_kwargs = {}
        # if self.backend == "bokeh":
        #     layout_kwargs["width"] = width
        #     layout_kwargs["height"] = height
        # else:
        #     layout_kwargs["fig_size"] = height/2

        scatter_kwargs = {}

        for k in kwargs:
            scatter_kwargs[k] = kwargs[k]
        if len(plot_list) > 1:
            final_plots = all_plots.opts(
                opts.Scatter(
                    xticks=xticks,
                    xaxis="bottom",
                    # fontsize=self.font_dic,
                    show_grid=True,
                    show_legend=True,
                    **scatter_kwargs,
                ),
                opts.Layout(**layout_kwargs),
            ).cols(col_wrap)
        else:
            final_plots = all_plots.opts(
                opts.Scatter(
                    xticks=xticks,
                    xaxis="bottom",
                    # fontsize=self.font_dic,
                    show_grid=True,
                    show_legend=True,
                    **scatter_kwargs,
                ),
                opts.Layout(**layout_kwargs),
            )

        if self.backend == "matplotlib":
            final_plots = final_plots.opts(sublabel_format="")

        if datashade:
            final_plots = final_plots.opts(plot=dict(width=width, height=height))

        renderer = hv.Store.renderers[self.backend].instance(fig="png", dpi=300)

        i = 0
        if prefix is None:
            prefix = "plt_grid_"

        if self.backend == "bokeh":
            hvplot.show(final_plots)
            self.plots.append(final_plots)
        else:
            filetype = ".png"
            self.plots.append(final_plots)

            save = plot_path + prefix + x_axis + "_" + "{:03d}".format(i)
            while os.path.isfile(save + filetype):
                i = i + 1
                save = plot_path + prefix + x_axis + "_" + "{:03d}".format(i)

            renderer.save(final_plots, save)
            display_file = save
            # Store every image as an individual one as well
            if plot_singles:
                for j, pl in enumerate(final_plots):
                    y_name = y_axis[j]
                    i = 0
                    save = (
                        plot_path + prefix + x_axis + "_" + y_name + "_{:03d}".format(i)
                    )
                    while os.path.isfile(save + filetype):
                        i = i + 1
                        save = (
                            plot_path
                            + prefix
                            + x_axis
                            + "_"
                            + y_name
                            + "_{:03d}".format(i)
                        )

                    pl = pl.opts(
                        opts.Scatter(
                            xticks=xticks,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **scatter_kwargs,
                        )
                    ).opts(aspect=aspect, **layout_kwargs)
                    print(f"Plot to {save}")
                    renderer.save(pl, save)
            try:
                from IPython.display import Image, display

                display(Image(display_file + filetype, width=width))
            except:
                pass

    def get_mse(
        self,
        out_params,
        others_df=None,
        others_path=None,
        in_params=None,
        ratio_type="adjusted",
        ratio_window=None,
        kind="mse",
        sens_kind="mean",
    ):
        """
        Create a pandas.Dataframe with columns
        MSE, Sensitivity, Output Parameter

        Parameters
        ----------
        out_params : list of string
            List of model output parameters to get the sensitivities for.
        others_df : Dataframe
            Dataframe with perturbed ensembles.
        others_path : path
            Path to dataframes with perturbed ensembles.
        in_params : list of string
            List of model parameters to consider.
        ratio_type : string
            "vanilla": Use the derivative ratio in the file.
            "adjusted": Can be added to any to any type below where the sensitivity
            is adjusted to the parameter value such that perturbing this value by a
            certain percentage gives an approximation of the resulting error/difference
            for any given hydrometeor.
            "per_timestep": Use the highest derivative per timestep as denominator.
            "window": Use the highest derivative in the given window by min_x and max_x.
            "per_xaxis": Use the highest derivative per x_axis value. If x_axis is "timestep"
            it is the same as "per_timestep".
            "x_per_out_param": Replace 'x' with any other option than "vanilla". Use the highest
            derivative but per output parameter. (that *should* be the vanilla version)
            "x_weighted": Append this to any ratio type to use the inverse of the
            output parameter value as weight. Only works with "x_per_out_param".
        ratio_window : dic of list
            Overides ratio_type and switches to a derivative ratio calculation
            per output param if "ratio_type" has "per_out_param" in it.
            Calculate the derivative ratio in windows
            where the key is a column name and the list consists of values
            at which a new window starts, i.e. {"T": [235, 273]} results
            in three windows with T < 235K, 235 <= T < 273K and 237K <= T.
        kind : string
            How to reduce the error from ensembles.
            mse: Calculate mean squared error
            maxse: Calculate maximum squared error
            nozeromse: Calculate mean squared error ignoring time steps
                where all simulations have zero values.
            sum : Cumulative squared error
            me : Mean error
            mae : Mean absolute error
        sens_kind : string
            How to reduce the predicted error from sensitivities.
            Be aware, that sensitivities can be calculated as squared
            predicted error already.
            max : Max error
            mean : Mean error
            squared_mean : Squared mean error

        Returns
        -------
        pandas.Dataframe with MSE, Sensitivity, Output Parameter and Perturbed
        Parameter
        """
        import pandas

        pandas.options.mode.chained_assignment = None  # default='warn'

        mean = self.cache
        if mean is None:
            print("Cached data is None. Make sure to cache your data using")
            print("self.cache_data()\nABORTING")
            return
        min_x = np.min(mean["time"])
        max_x = np.max(mean["time"])
        if others_df is None and others_path is None:
            print(
                "Either provide a loaded dataframe others_df or a path to the other files"
            )
            return
        if others_df is None and in_params is None:
            print("You need to provide in_params in addition to to others_df")
            return
        if kind == "mse":
            error_key = "MSE"
        elif kind == "maxse":
            error_key = "Max Error"
        elif kind == "nozeromse":
            error_key = "MSE (no zero)"
        elif kind == "sum":
            error_key = "Cumulative Squared Error"
        elif kind == "me":
            error_key = "Mean Error"
        elif kind == "mae":
            error_key = "Mean Absolute Error"

        if others_path is not None:
            files = sorted(glob(others_path + "*.nc_wcb"))
            df_dic = {
                error_key: np.asarray([]),
                "Sensitivity": np.asarray([]),
                "Output Parameter": np.asarray([]),
                "Perturbed Parameter": np.asarray([]),
                "Ratio Type": np.asarray([]),
            }
            # Pre calculate the sensitivities
            sens_dic = {}
            mean_dic = {}
            for i in pb(range(len(out_params))):
                # out_param in out_params:
                out_param = out_params[i]
                mean_dic[out_param] = np.reshape(
                    np.asarray(
                        mean.loc[mean["Output Parameter"] == out_param][out_param]
                    ),
                    (
                        len(
                            mean.loc[mean["Output Parameter"] == out_param][
                                out_param
                            ].index
                        ),
                        1,
                        1,
                    ),
                )
                # print(mean["time"])
                if isinstance(ratio_type, list):
                    sens_dic[out_param] = {}
                    for rt in ratio_type:
                        sens_dic[out_param][rt] = self._recalc_ratios(
                            mean.loc[mean["Output Parameter"] == out_param],
                            rt,
                            ratio_window,
                            in_params,
                        )
                else:
                    sens_dic[out_param] = self._recalc_ratios(
                        mean.loc[mean["Output Parameter"] == out_param],
                        ratio_type,
                        ratio_window,
                        in_params,
                    )

            for i in pb(range(len(files))):
                f = files[i]
                perturbed_param = "d" + f.split("/")[-1][0:-7]
                if perturbed_param not in in_params:
                    continue
                ds = xr.open_dataset(f, decode_times=False)[out_params]  #
                # print(f"\n{i}: {np.min(ds.time)}, {np.max(ds.time)}")
                ds = ds.sel(time=np.arange(min_x, max_x + 20, 20)).compute()
                for out_param in out_params:

                    def add_to_df(sens_df, rt):
                        if kind == "mse":
                            mean_squared_error = [
                                np.mean((ds[out_param] - mean_dic[out_param]) ** 2)
                            ]
                        if kind == "mae":
                            mean_squared_error = [
                                np.mean(np.abs(ds[out_param] - mean_dic[out_param]))
                            ]
                        if kind == "me":
                            mean_squared_error = [
                                np.mean((ds[out_param] - mean_dic[out_param]))
                            ]
                        elif kind == "maxse":
                            mean_squared_error = [
                                np.max((ds[out_param] - mean_dic[out_param]) ** 2)
                            ]
                        elif kind == "nozeromse":
                            mean_squared_error = (
                                ds[out_param] - mean_dic[out_param]
                            ) ** 2
                            shape = np.shape(mean_squared_error)
                            mean_squared_error = mean_squared_error.values.flatten()
                            mean_squared_error[mean_squared_error == 0] = np.nan
                            mean_squared_error = np.reshape(mean_squared_error, shape)
                            mean_squared_error = [np.nanmean(mean_squared_error)]
                        elif kind == "sum":
                            mean_squared_error = [
                                np.sum((ds[out_param] - mean_dic[out_param]) ** 2)
                            ]
                        df_dic[error_key] = np.append(
                            df_dic[error_key], mean_squared_error
                        )
                        if sens_kind == "max":
                            sens_idx = np.argmax(np.abs(sens_df[perturbed_param]))
                            df_dic["Sensitivity"] = np.append(
                                df_dic["Sensitivity"],
                                [sens_df[perturbed_param].iloc[sens_idx]],
                            )
                        elif sens_kind == "mean":
                            df_dic["Sensitivity"] = np.append(
                                df_dic["Sensitivity"],
                                [np.mean(sens_df[perturbed_param])],
                            )
                        elif sens_kind == "squared_mean":
                            df_dic["Sensitivity"] = np.append(
                                df_dic["Sensitivity"],
                                [
                                    np.mean(
                                        sens_df[perturbed_param]
                                        * sens_df[perturbed_param]
                                    )
                                ],
                            )

                        df_dic["Output Parameter"] = np.append(
                            df_dic["Output Parameter"], [out_param]
                        )
                        df_dic["Perturbed Parameter"] = np.append(
                            df_dic["Perturbed Parameter"], [perturbed_param]
                        )
                        df_dic["Ratio Type"] = np.append(df_dic["Ratio Type"], [rt])

                    if isinstance(ratio_type, list):
                        for rt in ratio_type:
                            add_to_df(sens_dic[out_param][rt], rt)
                    else:
                        add_to_df(sens_dic[out_param], ratio_type)
        else:
            if in_params is None:
                in_params = np.unique(others_df["Perturbed Parameter"])

            # Calculate MSE
            df_dic = {
                error_key: np.asarray([]),
                "Sensitivity": np.asarray([]),
                "Output Parameter": np.asarray([]),
                "Perturbed Parameter": np.asarray([]),
                "Ratio Type": np.asarray([]),
            }
            for i in pb(range(len(out_params))):
                out_param = out_params[i]
                tmp_df = mean.loc[mean["Output Parameter"] == out_param]
                mean_values = np.reshape(
                    np.asarray(tmp_df[out_param]), (1, 1, len(tmp_df[out_param].index))
                )

                def add_to_df(sens_df, rt):
                    if kind == "mse":
                        mean_squared_error = [
                            np.mean(
                                (
                                    others_df.sel({"Perturbed Parameter": mp})[
                                        out_param
                                    ]
                                    - mean_values
                                )
                                ** 2
                            )
                            for mp in in_params
                        ]
                    elif kind == "maxse":
                        # unfortunately numpy uses a lot of memory. This works only for small datasets
                        # mean_squared_error = ( (others_df[out_param]-mean_values)**2).mean( axis=(1, 2, 3) ).values
                        mean_squared_error = [
                            np.max(
                                (
                                    others_df.sel({"Perturbed Parameter": mp})[
                                        out_param
                                    ]
                                    - mean_values
                                )
                                ** 2
                            )
                            for mp in in_params
                        ]
                    elif kind == "nozeromse":
                        mean_squared_error = [
                            (
                                others_df.sel({"Perturbed Parameter": mp})[out_param]
                                - mean_values
                            )
                            ** 2
                            for mp in in_params
                        ]
                        mean_squared_error[mean_squared_error == 0] = np.nan
                        mean_squared_error = [
                            np.mean(mean_squared_error[j])
                            for j in range(len(in_params))
                        ]
                    df_dic[error_key] = np.append(df_dic[error_key], mean_squared_error)
                    sens_idx = [np.argmax(np.abs(sens_df[mp])) for mp in in_params]
                    df_dic["Sensitivity"] = np.append(
                        df_dic["Sensitivity"],
                        [
                            sens_df[mp].iloc[sens_idx[i]]
                            for i, mp in enumerate(in_params)
                        ],
                    )
                    # [ np.max(np.abs(sens_df[mp])) for mp in in_params ])
                    df_dic["Output Parameter"] = np.append(
                        df_dic["Output Parameter"], [out_param for _ in in_params]
                    )
                    df_dic["Perturbed Parameter"] = np.append(
                        df_dic["Perturbed Parameter"], in_params
                    )
                    df_dic["Ratio Type"] = np.append(
                        df_dic["Ratio Type"], [rt for _ in in_params]
                    )

                if isinstance(ratio_type, list):
                    for rt in ratio_type:
                        add_to_df(
                            self._recalc_ratios(tmp_df, rt, ratio_window, in_params), rt
                        )
                else:
                    add_to_df(
                        self._recalc_ratios(
                            tmp_df, ratio_type, ratio_window, in_params
                        ),
                        ratio_type,
                    )
        return pandas.DataFrame.from_dict(df_dic)

    def plot_mse(
        self,
        out_params,
        others_df=None,
        mse_df_=None,
        in_params=None,
        width=1959,
        height=1224,
        datashade=False,
        prefix=None,
        alpha=1,
        plot_path="pics/",
        yticks=10,
        xticks=10,
        decimals=3,
        ratio_type="vanilla",
        ratio_window=None,
        s=None,
        formatter_limits=None,
        log_x=False,
        log_y=False,
        confidence=None,
        hist=False,
        kind="grid_plot",
        abs_x=False,
        abs_y=False,
        split_ellipse=False,
        error_key="MSE",
        xlabel=None,
        sens_key="Sensitivity",
        ylabel=None,
        log_func=np.log,
        plot_singles=False,
        plot_types=True,
        title=None,
        linewidth=None,
        inf_val=None,
        **kwargs,
    ):
        """
        Plot mean squared error over maximal sensitivity calculated
        via the given ratio type (i.e. sensitivity per timestep) and given
        the trajectory to use for sensitivities and mean (i.e. the original trajectory
        from which model parameters have been perturbed).
        Supports only matplotlib at the moment.

        Parameters
        ----------
        out_params : list of string
            List of output parameters to plot the datapoints for.
        others_df : xarray.DataSet
            Dataframe with output parameters for each of the perturbed
            parameters.
            Dimensions are "Perturbed Parameter", "time", "trajectory",
            "ensemble". Either this or mse_df need to be given.
        mse_df : pandas.Dataframe
            Dataframe with pre calculated MSE values from get_mse().
        in_params : list of string
            List of model parameters to consider for the plot. If none is
            given, use all.
        width : int
            The width of the plot (bokeh) or the relative width of the plot
            (matplotlib).
        height : int
            The height of the plot (bokeh) or the relative height of the plot
            (matplotlib).
        datashade : bool
            Wether to use datashade for plotting the data.
        prefix : string
            Prefix to add to the filename.
        alpha : List of floats
            Alpha values for top [0] and bottom [1] graph.
        plot_path : string
            Path to the folder to save the image.
        yticks : int
            Number of ticks on all y-axis.
        xticks : int
            Number of ticks on all x-axis.
        decimals : int
            For rounding the right y-axis, if provided by "twin_axis".
            From numpy: Number of decimal places to round to (default: 0).
            If decimals is negative, it specifies the number of
            positions to the left of the decimal point.
        ratio_type : String
            "vanilla": Use the derivative ratio in the file.
            "adjusted": Can be added to any to any type below where the sensitivity
            is adjusted to the parameter value such that perturbing this value by a
            certain percentage gives an approximation of the resulting error/difference
            for any given hydrometeor.
            "per_timestep": Use the highest derivative per timestep as denominator.
            "window": Use the highest derivative in the given window by min_x and max_x.
            "per_xaxis": Use the highest derivative per x_axis value. If x_axis is "timestep"
            it is the same as "per_timestep".
            "x_per_out_param": Replace 'x' with any other option than "vanilla". Use the highest
            derivative but per output parameter. (that *should* be the vanilla version)
            "x_weighted": Append this to any ratio type to use the inverse of the
            output parameter value as weight. Only works with "x_per_out_param".
        ratio_window : dic of list
            Overides ratio_type and switches to a derivative ratio calculation
            per output param if "ratio_type" has "per_out_param" in it.
            Calculate the derivative ratio in windows
            where the key is a column name and the list consists of values
            at which a new window starts, i.e. {"T": [235, 273]} results
            in three windows with T < 235K, 235 <= T < 273K and 237K <= T.
        s : int
            Can be used to adjust the size of the scatter dots. As a rule of
            thumb: 13 is good for matplotlib backend but too large for
            bokeh.
        formatter_limits : tuple of ints
            Lower and upper limits for formatting x- and y-axis.
        log_x : bool
            Log x-axis.
        log_y : bool
            Log y-axis.
        confidence : float
            Plot a confidence ellipse around each sample with confidence
            between 0 and 1.
        hist : bool
            Plot histograms around each sample
        kind : string
            "single_plot": A single plot with all the data. May be
                off if different output parameters are used due to
                different scalings.
            "grid_plot": Multiple plots that can have histograms and
                confidence ellipses as well.
            "paper": Create a plot for the paper which sets figure
                size and removes the first histogram.
        abs_x : bool
            Use absolute value of sensitivities.
        abs_y : bool
            Use absolute value of MSE.
        split_ellipse : bool
            Create two confidence ellipses, on for sensitivity >= 0
            and another for sensitvity <= 0.
        error_key : string
            Name of column with actual errors.
        xlabel : string
            Alternative label for x-axis.
        sens_key : string
            Name of column with predicted errors aka sensitivities.
        ylabel : string
            Alternative label for y-axis.
        log_func: callable
            Function for logarithm (i.e. np.log10 or np.log (default))
        plot_singles : bool
            Plot every plot as an individual plot in addition to the grid of
            plots.
        plot_types : bool
            Wether to plot the color of the dots grouped by the type of
            input parameter.
        title : string
            Title used for the plot. If None is given, a standard
            title based on different labels is used.
        linewidth : int
            Width of ellipse.
        inf_val : float or None
            If float: Mark values at this x-value as negative infinity. The
            confidence ellipse does not take these values into account.
        kwargs : Dict of args
            Arguments are passed to ellipse, i.e. {"color": "red"}
        """
        if mse_df_ is None and others_df is None:
            print("Either mse_df or others_df for calculating mse_df must be given!")
            return
        if mse_df_ is None:
            mse_df = self.get_mse(
                out_params, others_df, in_params, ratio_type, ratio_window
            )
        else:
            mse_df = mse_df_.copy()
        # mse_df = mse_df.loc[mse_df["Sensitivity"] != 0]

        import hvplot.dask  # adds hvplot method to dask objects
        import hvplot.pandas
        from holoviews.operation.datashader import datashade as dsshade
        from holoviews import opts
        import holoviews as hv

        dsshade.dynamic = False
        if kind == "paper":
            fontscale = width / 350  # 500
        else:
            fontscale = width / 2200
        if formatter_limits is not None:
            matplotlib.rcParams["axes.formatter.limits"] = formatter_limits

        hv.extension(self.backend)

        aspect = width / height
        if s is None:
            scatter_size = int(height / 200)
        else:
            scatter_size = s

        hist_wh = int(height / 3)

        if self.backend == "matplotlib":
            scatter_kwargs = {"s": scatter_size}
        else:
            scatter_kwargs = {"size": scatter_size}

        fig_inches = width / 300
        if width < height:
            fig_inches = height / 300
        layout_kwargs = {}
        if self.backend == "bokeh":
            layout_kwargs["width"] = width
            layout_kwargs["height"] = height
        else:
            layout_kwargs["fig_inches"] = fig_inches

        if abs_x:
            mse_df[sens_key] = np.abs(mse_df[sens_key])
            if log_x:
                mse_df[sens_key] = log_func(mse_df[sens_key])
                if inf_val is not None:
                    inf_val = log_func(inf_val)
        if abs_y:
            mse_df[error_key] = np.abs(mse_df[error_key])
            if log_y:
                mse_df[error_key] = log_func(mse_df[error_key])

        if not abs_x and log_x:
            sign = np.sign(mse_df[sens_key])
            mse_df[sens_key] = log_func(np.abs(mse_df[sens_key]))
            mse_df[sens_key] = mse_df[sens_key] - np.min(mse_df[sens_key])
            mse_df[sens_key] *= sign
            sign = np.sign(inf_val)
            if inf_val is not None:
                inf_val = log_func(inf_val) * sign
        if not abs_y and log_y:
            sign = np.sign(mse_df[error_key])
            mse_df[error_key] = sign * log_func(np.abs(mse_df[error_key]))

        if xlabel is None:
            if log_x:
                xlabel = "Log " + sens_key
            else:
                xlabel = sens_key
        if ylabel is None:
            if log_y:
                ylabel = "Log " + error_key
            else:
                ylabel = error_key
        hspace = 0.05
        if title is None:
            title = ylabel + " over " + xlabel

        # Plot all into one plot
        if kind == "single_plot":
            if not datashade:
                if self.backend == "matplotlib":
                    if not plot_types:
                        cmap_values = []
                        for out_param in out_params:
                            cmap_values.append(
                                matplotlib.colors.to_hex(self.cmap_particles[out_param])
                            )
                            by_col = "Output Parameter"
                    else:
                        cmap_values = []
                        for group in np.unique(mse_df["Group"]):
                            cmap_values.append(
                                matplotlib.colors.to_hex(self.cmap_types[group])
                            )
                        by_col = "Group"

                    mse_plot = (
                        mse_df.hvplot.scatter(
                            x=sens_key,
                            y=error_key,
                            by=by_col,
                            title=title,
                            color=cmap_values,
                            datashade=datashade,
                            alpha=alpha,
                            legend=True,
                            yticks=yticks,
                            xticks=xticks,
                        )
                        .opts(aspect=aspect, fontscale=fontscale)
                        .options(xlabel=xlabel, ylabel=ylabel)
                    )

        all_opts = dict(
            title=title, hspace=0.05, shared_axes=False, fontscale=fontscale
        )
        all_opts_notitle = dict(hspace=0.05, shared_axes=False, fontscale=fontscale)
        if self.backend == "matplotlib":
            all_opts["tight"] = False
            all_opts["sublabel_format"] = ""
            all_opts_notitle["tight"] = False
            all_opts_notitle["sublabel_format"] = ""
        # Gridded plot with every output parameter on y
        # and sensitivity calculations on x
        # and different ratio types as "by"
        if kind == "grid_plot" or kind == "paper":
            # Only tested with matplotlib
            def correl_conf_ell(
                df, x, y, by=None, confidence=0.95, color=None, **kwargs
            ):
                new_kwargs = kwargs
                if color is not None:
                    new_kwargs["color"] = color
                    if linewidth is not None:
                        if self.backend == "matplotlib":
                            new_kwargs["linewidth"] = linewidth
                        else:
                            new_kwargs["line_width"] = linewidth
                chisq = chi2.ppf(confidence, 2)

                def create_ellipse(x, y, **kwargs):
                    cov = np.cov(x, y)
                    eigen, eigenv = np.linalg.eig(cov)
                    mean_x = np.mean(x)
                    mean_y = np.mean(y)
                    if eigen[0] > eigen[1]:
                        bigger_idx = 0
                    else:
                        bigger_idx = 1
                    rot = np.arctan2(eigenv[bigger_idx][1], eigenv[bigger_idx][0])
                    if rot < 0:
                        rot += np.pi
                    major_scale = chisq * np.sqrt(eigen[bigger_idx])
                    minor_scale = chisq * np.sqrt(eigen[(1 + bigger_idx) % 2])
                    if bigger_idx == 1:
                        scale = (major_scale, minor_scale)
                    else:
                        scale = (minor_scale, major_scale)
                    scale = (major_scale, minor_scale)
                    # Major axis if needed
                    # Seems a bit off at least on log log plots
                    # slope = eigenv[1][0]/eigenv[0][0]
                    # intercept = mean_y - slope*mean_x
                    # max_x = np.max(x)
                    # curve = hv.Curve(
                    #     [(i, intercept + slope*i) for i in np.arange(mean_x, max_x, (max_x-mean_x)/10)]
                    # ).opts(color="green", alpha=0.5)
                    # return (hv.Ellipse(mean_x, mean_y, scale, orientation=-rot).opts(**kwargs)
                    #     * curve
                    #      * hv.Text(
                    #         x.iloc[1], (intercept + slope*x.iloc[1])*1.1, f"Slope: {slope:.2f}" ) )
                    return hv.Ellipse(mean_x, mean_y, scale, orientation=-rot).opts(
                        **kwargs
                    )

                if by is None:
                    if inf_val is not None:
                        df_tmp = df.loc[df[x] != inf_val]
                    else:
                        df_tmp = df
                    x = df_tmp[x]
                    y = df_tmp[y]
                    if len(df_tmp[x]) == 1:
                        return None
                    return create_ellipse(x, y, **new_kwargs)
                else:
                    by_list = np.unique(df[by])
                    ells = None
                    for b in by_list:
                        new_kwargs = kwargs
                        if color is not None:
                            new_kwargs["color"] = color[b]
                        df_tmp = df.loc[df[by] == b]
                        if inf_val is not None:
                            df_tmp = df_tmp.loc[df_tmp[x] != inf_val]
                        if len(df_tmp[x]) == 1:
                            continue
                        if ells is None:
                            ells = create_ellipse(df_tmp[x], df_tmp[y], **new_kwargs)
                        else:
                            ells = ells * create_ellipse(
                                df_tmp[x], df_tmp[y], **new_kwargs
                            )
                    return ells

            if log_x:
                if kind == "paper":
                    tmp_df = mse_df.loc[mse_df["Output Parameter"] == out_params[0]]
                    min_x = (
                        np.min(tmp_df[sens_key]) - np.abs(np.min(tmp_df[sens_key])) / 10
                    )
                    max_x = 2
                else:
                    min_x = (
                        np.min(mse_df[sens_key]) - np.abs(np.min(mse_df[sens_key])) / 7
                    )
                    max_x = (
                        np.max(mse_df[sens_key]) + np.abs(np.max(mse_df[sens_key])) / 7
                    )
            else:
                delta_x = np.max(mse_df[sens_key]) - np.min(mse_df[sens_key])
                min_x = np.min(mse_df[sens_key]) - delta_x / 10
                max_x = np.max(mse_df[sens_key]) + delta_x / 10

            if not plot_types:
                cmap = plt.get_cmap("tab10")
                colors = {}
                cmap_values = []
                for i, rt in enumerate(np.unique(mse_df["Ratio Type"])):
                    colors[rt] = matplotlib.colors.to_hex(cmap(i)[0:-1])
                    cmap_values.append(colors[rt])
                by_col = "Ratio Type"
            else:
                cmap_values = []
                colors = {}
                for group in np.unique(mse_df["Group"]):
                    colors[group] = self.cmap_types[group]
                    cmap_values.append(self.cmap_types[group])
                by_col = "Group"

            def format_tick(v, pos, max_tick):
                if pos:
                    upper = v - max_tick
                    return f"1e{upper:2.1f}"
                else:
                    upper = max_tick - v
                    return f"-1e{upper:2.1f}"

            this_x_ticks = int(np.floor(xticks / 2))
            all_plots = None
            legend = False
            if kind == "paper":
                legend = "top_left"  # bottom
            else:
                if log_x and not abs_x:
                    pos_plot = mse_df.loc[mse_df[sens_key] >= 0].hvplot.hist(
                        y=sens_key,
                        by=by_col,
                        alpha=0.5,
                        legend=True,
                        title=title,
                        grid=False,
                        color=cmap_values,
                    )

                    max_y = 0
                    for key in pos_plot:
                        max_x_tick = np.max(key[sens_key])
                        min_x_tick = 0
                        delta_tick = max_x_tick / (this_x_ticks - 1)
                        max_y = np.max(key["Sensitivity_count"])

                    pos_x_ticks = [
                        (
                            min_x_tick + i * delta_tick,
                            format_tick(min_x_tick + i * delta_tick, True, max_x_tick),
                        )
                        for i in range(this_x_ticks)
                    ]

                    neg_plot = mse_df.loc[mse_df[sens_key] < 0].hvplot.hist(
                        y=sens_key,
                        by=by_col,
                        alpha=0.5,
                        legend=False,
                        grid=False,
                        color=cmap_values,
                    )

                    for key in neg_plot:
                        if np.max(key["Sensitivity_count"]) > max_y:
                            max_y = np.max(key["Sensitivity_count"])
                        max_x_tick = np.min(key[sens_key])
                        min_x_tick = 0
                        delta_tick = max_x_tick / (this_x_ticks - 1)

                    neg_x_ticks = [
                        (
                            min_x_tick + i * delta_tick,
                            format_tick(min_x_tick + i * delta_tick, False, max_x_tick),
                        )
                        for i in range(this_x_ticks)
                    ]
                    y_lim = (0, max_y)
                    opts_dic = {"fontscale": fontscale}
                    if self.backend == "matplotlib":
                        opts_dic["aspect"] = aspect / 2
                    pos_plot = (
                        pos_plot.opts(**opts_dic)
                        .options(xlabel="", yaxis=False)
                        .opts(xticks=pos_x_ticks, ylim=y_lim)
                    )
                    neg_plot = (
                        neg_plot.opts(**opts_dic)
                        .options(xlabel="")
                        .opts(xticks=neg_x_ticks, ylim=y_lim)
                    )

                    text_pos = hv.Text(0, 0, "/          ")
                    text_neg = hv.Text(0, 0, "          /")

                    all_plots = (neg_plot * text_neg + pos_plot * text_pos).opts(
                        **all_opts
                    )
                else:
                    opts_dic = {"fontscale": fontscale}
                    if self.backend == "matplotlib":
                        opts_dic["aspect"] = aspect
                    all_plots = (
                        mse_df.hvplot.hist(
                            y=sens_key,
                            by=by_col,
                            alpha=alpha,
                            legend=True,
                            grid=True,
                            title=title,
                            color=cmap_values,
                        )
                        .opts(**opts_dic)
                        .options(xlabel="")
                    )

            opts_dic = {"fontscale": fontscale}
            if self.backend == "matplotlib":
                opts_dic["aspect"] = aspect / 2
            opts_dic2 = {"fontscale": fontscale}
            if self.backend == "matplotlib":
                opts_dic2["aspect"] = aspect

            for out_param in out_params:
                if kind != "paper":
                    title = None
                tmp_df = mse_df.loc[mse_df["Output Parameter"] == out_param]
                min_y = tmp_df[error_key].min()
                max_y = tmp_df[error_key].max()
                delta_y = (max_y - min_y) / 20
                min_y -= delta_y
                max_y += delta_y
                if log_x and not abs_x:
                    if hist:
                        hspace = 0

                    split_ellipse = True
                    this_x_ticks = int(np.floor(xticks / 2))
                    # Make two scatter plots with correct ticks
                    pos_plot = tmp_df.loc[tmp_df[sens_key] >= 0].hvplot.scatter(
                        x=sens_key,
                        title=title,
                        y=error_key,
                        by=by_col,
                        datashade=False,
                        alpha=alpha,
                        legend=legend,
                        grid=False,
                        xlim=(-2, max_x),
                        ylim=(min_y, max_y),
                        color=cmap_values,
                    )

                    # Get cosest value to zero for proper ticks
                    max_x_tick = np.max(tmp_df[sens_key])
                    if -1 * np.min(tmp_df[sens_key]) > max_x_tick:
                        max_x_tick = -1 * np.min(tmp_df[sens_key])
                    delta_tick = max_x_tick / (this_x_ticks - 1)
                    pos_x_ticks = [
                        (
                            min_x_tick + i * delta_tick,
                            format_tick(min_x_tick + i * delta_tick, True, max_x_tick),
                        )
                        for i in range(this_x_ticks)
                    ]

                    neg_plot = tmp_df.loc[tmp_df[sens_key] < 0].hvplot.scatter(
                        x=sens_key,
                        y=error_key,
                        by=by_col,
                        datashade=False,
                        alpha=alpha,
                        legend=False,
                        grid=False,
                        xlim=(min_x, 2),
                        ylim=(min_y, max_y),
                        color=cmap_values,
                    )

                    neg_x_ticks = [
                        (
                            np.min(tmp_df[sens_key]) + i * delta_tick,
                            format_tick(
                                np.min(tmp_df[sens_key]) + i * delta_tick,
                                False,
                                -max_x_tick,
                            ),
                        )
                        for i in range(this_x_ticks)
                    ]

                    pos_plot = (
                        pos_plot.opts(opts.Scatter(**scatter_kwargs))  # s=scatter_size
                        .opts(xticks=pos_x_ticks, **opts_dic)
                        .options(
                            ylabel=ylabel,
                            xlabel="",
                            yaxis=False,
                        )
                    )

                    neg_plot = (
                        neg_plot.opts(opts.Scatter(**scatter_kwargs))  # s=scatter_size
                        .opts(xticks=neg_x_ticks, **opts_dic)
                        .options(ylabel=ylabel, xlabel="")
                    )

                    # Make three histograms
                    if hist:
                        xdist_pos = (
                            tmp_df.loc[tmp_df[sens_key] >= 0]
                            .hvplot.hist(
                                y=sens_key,
                                by=by_col,
                                alpha=alpha,
                                legend=False,
                                color=cmap_values,
                            )
                            .opts(fontscale=fontscale)
                        )
                        xdist_neg = (
                            tmp_df.loc[tmp_df[sens_key] < 0]
                            .hvplot.hist(
                                y=sens_key,
                                by=by_col,
                                alpha=alpha,
                                legend=False,
                                color=cmap_values,
                            )
                            .opts(fontscale=fontscale)
                        )
                        ydist = tmp_df.hvplot.hist(
                            y=error_key,
                            by=by_col,
                            alpha=alpha,
                            legend=False,
                            color=cmap_values,
                        ).opts(fontscale=fontscale)

                    # Make four ellipses
                    if confidence is not None:
                        ells_pos = correl_conf_ell(
                            df=tmp_df.loc[tmp_df[sens_key] >= 0],
                            x=sens_key,
                            y=error_key,
                            by=by_col,
                            confidence=confidence,
                            color=colors,
                        )
                        ells_neg = correl_conf_ell(
                            df=tmp_df.loc[tmp_df[sens_key] <= 0],
                            x=sens_key,
                            y=error_key,
                            by=by_col,
                            confidence=confidence,
                            color=colors,
                        )
                        if ells_neg is not None:
                            neg_plot = neg_plot * ells_neg
                        if ells_pos is not None:
                            pos_plot = pos_plot * ells_pos

                    # Add small break lines
                    text_pos = hv.Text(0, min_y, "/   ")
                    text_neg = hv.Text(0, min_y, "   /")

                    # put it all together
                    if hist:
                        if self.backend == "bokeh":
                            all_neg_plot = (
                                neg_plot * text_neg
                                << hv.Empty()
                                << xdist_neg.opts(
                                    xaxis="bare", xlim=(min_x, 2), width=hist_wh
                                )
                            )
                            all_pos_plot = (
                                pos_plot * text_pos
                                << ydist.opts(
                                    yaxis="bare", xlim=(min_y, max_y), height=hist_wh
                                )
                                << xdist_pos.opts(
                                    xaxis="bare", xlim=(-2, max_x), width=hist_wh
                                )
                            )
                        else:
                            all_neg_plot = (
                                neg_plot * text_neg
                                << hv.Empty()
                                << xdist_neg.opts(xaxis="bare", xlim=(min_x, 2))
                            )
                            all_pos_plot = (
                                pos_plot * text_pos
                                << ydist.opts(yaxis="bare", xlim=(min_y, max_y))
                                << xdist_pos.opts(xaxis="bare", xlim=(-2, max_x))
                            )
                        mse_plot = all_neg_plot + all_pos_plot
                    else:
                        mse_plot = neg_plot * text_neg + pos_plot * text_pos

                    if all_plots is None:
                        all_plots = mse_plot
                    else:
                        all_plots += mse_plot
                else:
                    # experimental
                    # Assume 10^-80 (or here -80) for kind==paper is
                    # negative infinity
                    if (
                        kind == "paper" and inf_val is not None
                    ):  # np.min(tmp_df[sens_key]) == -80:
                        xticks_list = [(inf_val, "-Infinity")]
                        tick_val = int(np.ceil(inf_val / 10.0)) * 10
                        delta_tick = int(-tick_val / (xticks - 1))
                        tick_val += delta_tick
                        while tick_val < 0:
                            xticks_list.append((tick_val, tick_val))
                            tick_val += delta_tick
                        xticks_list.append((0, 0))

                        mse_plot = (
                            tmp_df.hvplot.scatter(
                                x=sens_key,
                                y=error_key,
                                by=by_col,
                                datashade=False,
                                alpha=alpha,
                                legend=legend,
                                title=title,
                                grid=True,
                                xlim=(min_x, max_x),
                                ylim=(min_y, max_y),
                                color=cmap_values,
                                width=width,
                                height=height,
                            )
                            .opts(opts.Scatter(**scatter_kwargs))  # s=scatter_size
                            .opts(**opts_dic2)
                            .options(
                                ylabel=ylabel,
                                xlabel="",
                                xticks=xticks_list,
                            )
                        )
                        # Save space by removing legend title
                        mse_plot.get_dimension(by_col).label = ""
                    else:
                        mse_plot = (
                            tmp_df.hvplot.scatter(
                                x=sens_key,
                                y=error_key,
                                by=by_col,
                                datashade=False,
                                alpha=alpha,
                                legend=legend,
                                title=title,
                                grid=True,
                                xlim=(min_x, max_x),
                                ylim=(min_y, max_y),
                                color=cmap_values,
                                width=width,
                                height=height,
                            )
                            .opts(opts.Scatter(**scatter_kwargs))  # s=scatter_size
                            .opts(**opts_dic2)
                            .options(ylabel=ylabel, xlabel="")
                        )

                    # Adding histograms around those
                    if hist:
                        xdist = tmp_df.hvplot.hist(
                            y=sens_key,
                            by=by_col,
                            alpha=alpha,
                            legend=False,
                            color=cmap_values,
                        ).opts(fontscale=fontscale)
                        ydist = tmp_df.hvplot.hist(
                            y=error_key,
                            by=by_col,
                            alpha=alpha,
                            legend=False,
                            color=cmap_values,
                        ).opts(fontscale=fontscale)

                    if out_param == out_params[-1]:
                        mse_plot = mse_plot.options(xlabel=xlabel)

                    if confidence is not None:
                        if split_ellipse:
                            ells_pos = correl_conf_ell(
                                df=tmp_df.loc[tmp_df[sens_key] >= 0],
                                x=sens_key,
                                y=error_key,
                                by=by_col,
                                confidence=confidence,
                                color=colors,
                            )
                            ells_neg = correl_conf_ell(
                                df=tmp_df.loc[tmp_df[sens_key] <= 0],
                                x=sens_key,
                                y=error_key,
                                by=by_col,
                                confidence=confidence,
                                color=colors,
                            )
                            if ells_pos is not None and ells_neg is not None:
                                mse_plot = mse_plot * ells_pos * ells_neg
                            elif ells_pos is None:
                                mse_plot = mse_plot * ells_neg
                            elif ells_neg is None:
                                mse_plot = mse_plot * ells_pos
                        elif kind == "paper":
                            # Assuming -80 are the infinity placements
                            tmp_df_ell = tmp_df.loc[tmp_df[sens_key] > -80]
                            ells = correl_conf_ell(
                                df=tmp_df_ell,
                                x=sens_key,
                                y=error_key,
                                by=by_col,
                                confidence=confidence,
                                color=colors,
                            )
                            if ells is not None:
                                mse_plot = mse_plot * ells
                        else:
                            ells = correl_conf_ell(
                                df=tmp_df,
                                x=sens_key,
                                y=error_key,
                                by=by_col,
                                confidence=confidence,
                                color=colors,
                            )
                            if ells is not None:
                                mse_plot = mse_plot * ells
                    if hist:
                        if self.backend == "bokeh":
                            mse_plot = (
                                mse_plot.opts(width=width, height=height).options(
                                    # Legend scatter size is very small by default
                                    **{
                                        "glyph_height": scatter_size * 5,
                                        "glyph_width": scatter_size * 5,
                                    }
                                )
                                << ydist.opts(
                                    yaxis="bare",
                                    xaxis="bare",
                                    xlim=(min_y, max_y),
                                    width=hist_wh,
                                    height=height,
                                )
                                << xdist.opts(
                                    xaxis="bare",
                                    yaxis="bare",
                                    xlim=(min_x, max_x),
                                    height=hist_wh,
                                    width=width,
                                )
                            )
                        else:
                            mse_plot = (
                                mse_plot.options(
                                    # Legend scatter size is very small by default
                                    **{
                                        "glyph_height": scatter_size * 5,
                                        "glyph_width": scatter_size * 5,
                                    }
                                )
                                << ydist.opts(
                                    yaxis="bare", xaxis="bare", xlim=(min_y, max_y)
                                )
                                << xdist.opts(
                                    xaxis="bare", yaxis="bare", xlim=(min_x, max_x)
                                )
                            )

                    if all_plots is None:
                        all_plots = mse_plot
                    else:
                        all_plots += mse_plot
            if log_x and not abs_x:
                if kind == "paper":
                    mse_plot = all_plots.opts(**all_opts_notitle).opts(
                        opts.Layout(**layout_kwargs)
                    )
                else:
                    mse_plot = (
                        all_plots.opts(**all_opts_notitle)
                        .cols(2)
                        .opts(opts.Layout(**layout_kwargs))
                    )
            else:
                if kind == "paper":
                    if self.backend == "matplotlib":
                        mse_plot = all_plots.opts(sublabel_format="", tight=True).opts(
                            opts.Layout(**layout_kwargs)
                        )
                    else:
                        mse_plot = all_plots.opts(opts.Layout(**layout_kwargs))
                else:
                    if self.backend == "matplotlib":
                        mse_plot = (
                            all_plots.cols(1)
                            .opts(sublabel_format="", tight=True)
                            .opts(opts.Layout(**layout_kwargs))
                        )
                    else:
                        mse_plot = all_plots.cols(1).opts(opts.Layout(**layout_kwargs))

        # Save image to disk
        # The renderer should be able to switch between
        # matplotlib and bokeh but the colors can be wrong with
        # matplotlib
        renderer = hv.Store.renderers[self.backend].instance(fig="png", dpi=300)
        filetype = ".png"

        i = 0
        if prefix is None:
            prefix = error_key

        if self.backend == "bokeh":
            # from bokeh.io import export_svgs
            self.plots.append(mse_plot)
            save = plot_path + prefix + "_" + "{:03d}".format(i)
            while os.path.isfile(save + filetype):
                i = i + 1
                save = plot_path + prefix + "_" + "{:03d}".format(i)
            renderer.save(mse_plot, save)
            hvplot.show(mse_plot)
        else:

            self.plots.append(mse_plot)

            save = plot_path + prefix + "_" + "{:03d}".format(i)
            while os.path.isfile(save + filetype):
                i = i + 1
                save = plot_path + prefix + "_" + "{:03d}".format(i)

            renderer.save(mse_plot, save)
            display_file = save

            if plot_singles:
                for j, pl in enumerate(mse_plot):
                    i = 0
                    save = plot_path + prefix + "_" + "{:03d}".format(i)
                    while os.path.isfile(save + filetype):
                        i = i + 1
                        save = plot_path + prefix + "_" + "{:03d}".format(i)
                    if j > 0:
                        if isinstance(pl, hv.core.overlay.Overlay):
                            pl = pl.opts(
                                opts.Scatter(
                                    # s=scatter_size,
                                    alpha=alpha,
                                    show_grid=True,
                                    show_legend=True,
                                    **scatter_kwargs,
                                ),  # s=scatter_size),
                                opts.Layout(**layout_kwargs),
                            ).opts(
                                aspect=aspect,
                                fontscale=fontscale,
                                xlabel=xlabel,
                                xaxis="bottom",
                                title=title,
                            )
                        else:
                            title = "Validation Sensitivity " + out_params[j - 1]
                            # hv.core.layout.AdjointLayout
                            pl.__setitem__(
                                "main",
                                pl.main()
                                .opts(
                                    opts.Scatter(
                                        # s=scatter_size,
                                        alpha=alpha,
                                        show_grid=True,
                                        show_legend=True,
                                        **scatter_kwargs,
                                    ),
                                    opts.Layout(**layout_kwargs),
                                )
                                .opts(
                                    aspect=aspect,
                                    fontscale=fontscale,
                                    xlabel=xlabel,
                                    xaxis="bottom",
                                    title=title,
                                ),
                            )
                    else:
                        pl = pl.opts(show_legend=False)
                    print(f"Plot to {save}")
                    renderer.save(pl, save)

            try:
                from IPython.display import Image, display

                display(Image(display_file + filetype, width=width))
            except:
                pass
