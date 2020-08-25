try:
    import scripts.loader as loader
except:
    import loader
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
    n_timesteps : int
        Number of timesteps.
    clusternames : Dic of lsit of string
        Keys are output parameters where the value is a list of column names
        that hold a clustering assignment.
    threads : int
            Number of threads to use.
    cmap : Dic
        Dic of typenames of trajectories and colors.
    cache : pandas.Dataframe
        A dataframe after compute() was used that can be used as cache.
    """
    data = {}
    n_timesteps = 0
    cluster_names = {}
    threads = 0
    font_dic = {}
    backend = "matplotlib"
    widget = None
    colors = None
    cmap = None
    cache = None

    def __init__(self, direc, parquet=True, columns=None, backend="matplotlib"):
        """
        Init class by loading the data from the given path.

        Parameters
        ----------
        direc : string
            A path to a directory wit a list of files to read.
        parquet : boolean
            If true: Load a series of preprocessed parquet files,
            else load *.txt files.
        columns : list of strings
            Specify the columns to load.
        backend : String
            Either matplotlib or bokeh
        """
        self.data = dask_loader.load_mult_derivates_direc_dic(
            direc, parquet, columns
        )

        self.n_timesteps = len(self.data["timestep"].unique().compute())
        self.cluster_names = {}
        self.font_dic = {
            "title": 20,
            "labels": 20,
            "xticks": 12,
            "yticks": 16,
            "legend": 16
        }
        self.backend = backend
        self.plots = []
        if backend == "matplotlib":
            colors = plt.get_cmap("tab20c")
            self.cmap = {"Slantwise 600hPa 25. Quantile":  matplotlib.colors.to_hex(colors(0)[0:-1]),
                "Slantwise 600hPa 50. Quantile":  matplotlib.colors.to_hex(colors(1)[0:-1]),
                "Slantwise 600hPa 75. Quantile":  matplotlib.colors.to_hex(colors(2)[0:-1]),
                "Slantwise 400hPa 25. Quantile":  matplotlib.colors.to_hex(colors(4)[0:-1]),
                "Slantwise 400hPa 50. Quantile":  matplotlib.colors.to_hex(colors(5)[0:-1]),
                "Slantwise 400hPa 75. Quantile":  matplotlib.colors.to_hex(colors(6)[0:-1]),
                "Convective 600hPa 25. Quantile": matplotlib.colors.to_hex(colors(8)[0:-1]),
                "Convective 600hPa 50. Quantile": matplotlib.colors.to_hex(colors(9)[0:-1]),
                "Convective 600hPa 75. Quantile": matplotlib.colors.to_hex(colors(10)[0:-1]),
                "Convective 400hPa 25. Quantile": matplotlib.colors.to_hex(colors(12)[0:-1]),
                "Convective 400hPa 50. Quantile": matplotlib.colors.to_hex(colors(13)[0:-1]),
                "Convective 400hPa 75. Quantile": matplotlib.colors.to_hex(colors(14)[0:-1])}
            self.colors = colors
        else:
            colors = Category20c[20]
            self.cmap = {"Slantwise 600hPa 25. Quantile":  colors[0],
                "Slantwise 600hPa 50. Quantile":  colors[1],
                "Slantwise 600hPa 75. Quantile":  colors[2],
                "Slantwise 400hPa 25. Quantile":  colors[4],
                "Slantwise 400hPa 50. Quantile":  colors[5],
                "Slantwise 400hPa 75. Quantile":  colors[6],
                "Convective 600hPa 25. Quantile": colors[8],
                "Convective 600hPa 50. Quantile": colors[9],
                "Convective 600hPa 75. Quantile": colors[10],
                "Convective 400hPa 25. Quantile": colors[12],
                "Convective 400hPa 50. Quantile": colors[13],
                "Convective 400hPa 75. Quantile": colors[14]}
            self.colors = colors


    def to_parquet(self, f_name, compression="snappy"):
        append = not os.listdir(f_name)
        append = not append
        pd.to_parquet(self.data, f_name, append=append, ignore_divisions=append, compression=compression)

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

    def get_n_timesteps(self):
        """
        Get the number of timesteps of the data.

        Returns
        -------
        int
            Number of timesteps
        """
        return self.n_timesteps

    def cache_data(self, in_params, out_params, x_axis="timestep",
        y_axis=None,mapped=None,
        trajectories=None, frac=None, min_x=None, max_x=None, nth=None,
        compute=True):
        """
        Load some data into memory and store it as a pandas dataframe.
        May speed up some routines if the data fits into memory.

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
        mapped : string
            Column name which has to be true such as conv_400, slan_400,
            conv_600, slan_600.
        trajectories : List of int
            The index of the trajectories to plot. If None is given, all
            trajectories will be plotted.
        frac : float
            Sample a given fraction of rows. Deactivates "nth".
        nth : int
            Sample every nth entry. Works fast with "timestep" as x-axis and
            a given min_x and max_x value. If x_axis is any other than
            "timestep", an errorband triggered by "percentile" may not
            be plotted.
        min_x : float
            Minimum value for the x-axis.
        max_x : float
            Maximum value for the x-axis.
        nth : int
            Sample every nth entry. Works fast with "timestep" as x-axis and
            a given min_x and max_x value. If x_axis is any other than
            "timestep", an errorband triggered by "percentile" may not
            be plotted.
        compute : bool
            If true, store a pandas dataframe in self.cache. Otherwise dask
            dataframe (not entirely loaded) is stored.
        """
        t = timer()
        if frac is not None:
            df = self.data.sample(frac=frac, replace=False, random_state=42)
        elif nth is not None:
            if min_x is not None and max_x is not None and x_axis == "timestep":
                steps = np.arange(min_x, max_x, nth*2.5)
            elif x_axis == "timestep":
                df_tmp = self.data.timestep.unique().compute()
                min_val = df_tmp.min()
                max_val = df_tmp.max()
                steps = np.arange(min_val, max_val, nth*2.5)
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
            df = df.loc[df.trajectory.isin(trajectories)]
        if mapped is not None:
            df = df.loc[df[mapped]]

        df = df.loc[df["Output Parameter"].isin(out_params)]
        all_params = list(set(["Output Parameter", "trajectory", "type"] + in_params + [x_axis] + out_params))
        if y_axis is not None and not y_axis in all_params:
            all_params.append(y_axis)
        if compute:
            self.cache = df[all_params].compute()
        else:
            self.cache = df[all_params]
        t2 = timer()
        print("Loading done in {} s".format(t2-t))
        t = timer()

    def plot_two_ds(self, in_params, out_params, x_axis="timestep",
        y_axis=None, twin_axis=None, decimals=-3, mapped=None,
        trajectories=None, scatter=False, n_plots=None, percentile=None,
        frac=None, min_x=None, max_x=None, nth=None, hist=[False, False],
        hexbin=[False, False], bins=50,  log=[False, False], sort=True,
        scatter_deriv=False, line_deriv=False, compute=False,
        use_cache=False, errorband=False, plot_path="pics/", prefix=None,
        fig_type='svg', datashade=True, by=None,  alpha=[1, 1],
        rolling_avg=20, rolling_avg_par=20, max_per_deriv=10,
        width=1280, height=800, ratio_type="vanilla", ratio_window=None,
        vertical_mark=None, **kwargs):
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
        kwargs : dict
            Keyword arguments are passed down matplotlib.
        """
        import hvplot.dask # adds hvplot method to dask objects
        import hvplot.pandas
        import hvplot
        from holoviews import opts
        import holoviews as hv
        from holoviews.operation import histogram as hv_histo
        import pandas
        import dask.array as da

        hv.extension(self.backend)
        deriv_col_name = "Derivative Ratio"

        if not use_cache:
            self.cache_data(in_params, out_params, x_axis, y_axis, mapped,
                trajectories, frac, min_x, max_x, nth, compute)
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
                    series = df_tmp[in_params].rolling(
                        rolling_avg, min_periods=1).mean()
                    series.index = df.loc[df["Output Parameter"] == out_par].index
                    df.update(series)
                else:
                    types = df_tmp[by].unique()
                    for ty in types:
                        df_tmp2 = df_tmp.loc[df_tmp[by] == ty]
                        tmp_par = ["dcloud_a_vel", x_axis]
                        df_tmp2 = df_tmp.loc[df_tmp[by] == ty]
                        series = df_tmp2[in_params].rolling(
                            rolling_avg, min_periods=1).mean()
                        series.index = df.loc[
                            (df["Output Parameter"] == out_par) & (df[by] == ty) ].index
                        df.update(series)
            t2 = timer()
            print("Calculating rolling average for derivatives done in {} s".format(t2-t))
            t = timer()

        if rolling_avg_par is not None:
            for out_par in out_params:
                df_tmp = df.loc[df["Output Parameter"] == out_par]
                if by is None:
                    series = df_tmp[out_par].rolling(
                        rolling_avg_par, min_periods=1).mean()
                    series.index = df.loc[df["Output Parameter"] == out_par].index
                    df.update(series)
                else:
                    types = df_tmp[by].unique()
                    for ty in types:
                        df_tmp2 = df_tmp.loc[df_tmp[by] == ty]
                        series = df_tmp2[out_par].rolling(
                            rolling_avg_par, min_periods=1).mean()
                        series.index = df.loc[
                            (df["Output Parameter"] == out_par) & (df[by] == ty) ].index
                        df.update(series)
            t2 = timer()
            print("Calculating rolling average for outputs done in {} s".format(t2-t))
            t = timer()

        def recalc_ratios(ratio_df):
            # TODO: Add if by is not None: go over all types and do that stuff.
            if ratio_window is not None:
                col = list(ratio_window.keys())[0]
                edges = np.sort(ratio_window[col])
                for i in range(len(edges)+1):
                    if i == 0:
                        df_edge = ratio_df.loc[ratio_df[col] < edges[i]]
                        max_val = df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        if max_val == 0:
                            continue
                        ratio_df[in_params] = ratio_df.apply(
                            lambda x: x[in_params]/max_val if x[col] < edges[i] else x[in_params], axis=1)
                    elif i == len(edges):
                        df_edge = ratio_df.loc[ratio_df[col] >= edges[i-1]]
                        max_val = df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        if max_val == 0:
                            continue
                        ratio_df[in_params] = ratio_df.apply(
                            lambda x: x[in_params]/max_val if x[col] >= edges[i-1] else x[in_params], axis=1)
                    else:
                        df_edge = ratio_df.loc[(ratio_df[col] < edges[i]) & (ratio_df[col] >= edges[i-1])]
                        max_val = df_edge[in_params].apply(lambda x: np.max(np.abs(x))).max()
                        if max_val == 0:
                            continue
                        ratio_df[in_params] = ratio_df.apply(
                            lambda x: x[in_params]/max_val if x[col] < edges[i] and x[col] >= edges[i-1] else x[in_params], axis=1)
            elif "per_timestep" in ratio_type:
                # Get series of max values over all timesteps (equals index)
                max_vals = ratio_df[in_params].apply(lambda x: np.max(np.abs(x)), axis=1)
                ratio_df[in_params] = ratio_df[in_params].div(max_vals, axis=0)
                t2 = timer()
                print("Recalculating ratios done in {} s".format(t2-t))
            elif "window" in ratio_type:
                max_val = ratio_df[in_params].apply(lambda x: np.max(np.abs(x))).max()
                ratio_df[in_params] = ratio_df[in_params].div(max_val)
                t2 = timer()
                print("Recalculating ratios done in {} s".format(t2-t))
            elif "per_xaxis" in ratio_type:
                ratio_df = ratio_df.set_index(x_axis)
                max_vals = ratio_df.groupby(x_axis)[in_params].apply(lambda x: np.max(np.abs(x))).max(axis=1)
                ratio_df[in_params] = ratio_df[in_params].div(max_vals, axis="index")
                t2 = timer()
                print("Recalculating ratios done in {} s".format(t2-t))
            return ratio_df

        if not "per_out_param" in ratio_type:
            df = recalc_ratios(df)
            t = timer()

        set_yaxis = False
        if y_axis is None:
            set_yaxis = True
        for out_par in out_params:
            t = timer()
            if set_yaxis:
                y_axis = out_par
            df_tmp_out = df.loc[df["Output Parameter"] == out_par]
            if "per_out_param" in ratio_type:
                df_tmp_out = recalc_ratios(df_tmp_out)
                t = timer()

            # This is for the matplotlib backend
            # For bokeh, see https://github.com/holoviz/holoviews/issues/396
            twin = None
            if twin_axis is not None and self.backend == "matplotlib":
                lower_y = np.around(df_tmp_out[twin_axis].min(), decimals=decimals)
                upper_y = np.around(df_tmp_out[twin_axis].max(), decimals=decimals)
                tick_size = (upper_y - lower_y)/10
                def twinx(plot, element):
                    ax = plot.handles["axis"]
                    twinax = ax.twinx()
                    twinax.set_ylim((lower_y - tick_size * 0.75,
                                    upper_y + tick_size * 0.75))
                    twinax.set_yticks(np.arange(lower_y, upper_y + tick_size*0.75, tick_size))
                    twinax.set_ylabel(
                        latexify.parse_word(twin_axis)
                        + " " + latexify.get_unit(twin_axis, brackets=True))
                    plot.handles["axis"] = twinax
                if log[1]:
                    df_tmp_out[twin_axis] = da.log(da.fabs(df_tmp_out[twin_axis]))
                twin = df_tmp_out.hvplot.scatter(
                    x=x_axis,
                    y=twin_axis,
                    color="gray",
                    datashade=datashade,
                    legend=False,
                ).opts(initial_hooks=[twinx], apply_ranges=False).opts(opts.Scatter(s=8))

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
            print("Sorting done in {} s".format(t2-t), flush=True)
            def plot_helper(df, in_params, prefix, **kwargs):
                deriv_col_name = "Derivative Ratio"
                # following https://holoviz.org/tutorial/Composing_Plots.html
                t = timer()
                add_cols = []
                if vertical_mark is not None:
                    add_cols = list(vertical_mark.keys())

                if by is not None:
                    df_tmp = df[in_params+[x_axis, by] + add_cols]
                    df_tmp = df_tmp.melt([x_axis, by] + add_cols, var_name="Derivatives",
                                    value_name=deriv_col_name)
                else:
                    df_tmp = df[in_params+[x_axis] + add_cols]
                    df_tmp = df_tmp.melt([x_axis] + add_cols, var_name="Derivatives",
                                         value_name=deriv_col_name)
                t2 = timer()
                print("Melting done in {} s".format(t2-t), flush=True)

                if log[1]:
                    df_tmp[deriv_col_name] = da.log(da.fabs(df_tmp["Derivative Ratio"]))
                    # Remove zero entries (-inf)
                    df_tmp = df_tmp[~da.isinf(df_tmp[deriv_col_name])]
                    df_tmp.rename(columns={deriv_col_name: "Log Derivative Ratio"}, inplace=True)
                    deriv_col_name = "Log Derivative Ratio"
                    t2 = timer()
                    print("Log done in {} s".format(t2-t))
                    t = timer()
                df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)

                if percentile is not None:
                    if y_axis == x_axis:
                        print("x-axis and y-axis are the same. Cannot plot that with percentiles!")
                        return
                    else:
                        df_group = df[[x_axis, y_axis] + add_cols]
                    if log[0]:
                        # Apply should be more expensive
                        df_group[y_axis] = da.log(da.fabs(df_group[y_axis]))
                        # Remove zero entries (-inf)
                        df_group = df_group[~da.isinf(df_group[y_axis])]

                    # Group for min, max and percentiles
                    funcs = [np.min, np.max] + [lambda x, perc=perc: np.percentile(x, perc, axis=0) for perc in percentile]
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
                    param_plot = (
                        df_min_max.hvplot.area(
                            x=x_axis,
                            y="Min",
                            y2="Max",
                            alpha=0.5,
                            ylabel=latexify.parse_word(y_axis),
                            title="Values of of {}".format(latexify.parse_word(y_axis)),
                            label="Spread",
                            color="grey")
                        * df_min_max.hvplot.line(
                            x=x_axis,
                            y=p_list,
                            ylabel=latexify.parse_word(y_axis),
                            **kwargs)
                    )
                elif errorband:
                    df_group = df[[x_axis, y_axis, "trajectory"] + add_cols]
                    if log[0]:
                        # Apply should be more expensive
                        df_group[y_axis] = da.log(da.fabs(df_group[y_axis]))
                    df_min_max = df_group.groupby([x_axis, "trajectory"])[y_axis].mean().groupby(x_axis).agg([np.min, np.max])
                    df_std = df_group.groupby([x_axis, "trajectory"])[y_axis].mean().groupby(x_axis).agg([lambda x: -1*np.std(x)+np.mean(x), lambda x: np.std(x)+np.mean(x)])
                    df_mean = df_group.groupby(x_axis)[y_axis].mean()

                    param_plot = (
                        df_min_max.hvplot.area(
                            x=x_axis,
                            y="amin",
                            y2="amax",
                            alpha=0.5,
                            value_label=latexify.parse_word(y_axis),
                            label="Spread")
                        * df_mean.hvplot()
                        * df_std.hvplot.area(
                            x=x_axis,
                            y="<lambda_0>",
                            y2="<lambda_1>",
                            alpha=0.3,
                            value_label=latexify.parse_word(y_axis),
                            label="sd")
                    )
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
                                {types[i]: hv.Scatter((np.NaN, np.NaN)).opts(opts.Scatter(s=50, color=cmap_values[i]))
                                for i in range(len(types)) }
                            )

                            param_plot = df_group[df_group[y_axis] != 0].hvplot.scatter(
                                x=x_axis,
                                y=y_axis,
                                by=by,
                                title="Values of of {}".format(latexify.parse_word(y_axis)),
                                color=cmap_values,
                                label=None,
                                datashade=datashade,
                                alpha=alpha[0],
                                legend=False
                            ).opts(opts.Scatter(s=8)) * overlay

                            if hist[0]:
                                xhist = df_group[df_group[y_axis] != 0].hvplot.hist(y=x_axis, bins=bins, height=125)
                                yhist = df_group[df_group[y_axis] != 0].hvplot.hist(y=y_axis, bins=bins, width=125)
                                param_plot = param_plot << yhist << xhist
                        else:
                            param_plot = df_group[df_group[y_axis] != 0].hvplot.scatter(
                                x=x_axis,
                                y=y_axis,
                                by=by,
                                title="Values of of {}".format(latexify.parse_word(y_axis)),
                                color=cmap_values,
                                label=None,
                                datashade=datashade,
                                alpha=alpha[0]
                            ).opts(opts.Scatter(size=2))
                    else:
                        param_plot = df_group[df_group[y_axis] != 0].hvplot.scatter(
                            x=x_axis,
                            y=y_axis,
                            by=by,
                            title="Values of of {}".format(latexify.parse_word(y_axis)),
                            cmap=cmap,
                            label=None,
                            datashade=datashade
                        ).opts(aspect=3.2)
                        if self.backend == "bokeh":
                            points = hv.Points(
                                ([np.NaN for i in range(len(list(cmap.keys())))],
                                 [np.NaN for i in range(len(list(cmap.keys())))],
                                 list(cmap.keys())), vdims=by)
                            legend = points.options(
                                color=by,
                                cmap=cmap,
                                show_legend=True)
                            param_plot = (legend * param_plot).opts(aspect=3.2)
                elif hexbin[0]:
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
                        title="Values of of {}".format(latexify.parse_word(y_axis)),
                        label=None,
                        clabel="Count",
                        cmap="viridis",
                        logz=True,
                        gridsize=100
                    )
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
                        title="Values of of {}".format(latexify.parse_word(y_axis)),
                        label=None,
                        datashade=datashade
                    )

                t2 = timer()
                print("Setting up upper plot done in {} s".format(t2-t))
                t = timer()

                layout_kwargs = {}
                if self.backend == "bokeh":
                    layout_kwargs["width"] = width
                    layout_kwargs["height"] = height
                else:
                    if not hist[0]:
                        layout_kwargs["aspect"] = width/height
                    layout_kwargs["fig_size"] = height/2

                if hexbin[1]:
                    deriv_plot = df_tmp.hvplot.hexbin(
                        x=x_axis,
                        y=deriv_col_name,
                        by="Derivatives",
                        title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
                        label=None,
                        logz=True,
                        gridsize=100,
                        clabel="Count",
                        cmap="viridis",
                        colorbar=True
                    )
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
                            cmap_values.append(matplotlib.colors.to_hex(colors(i)[0:-1]))
                        if twin is not None:
                            cmap_values.append("gray")
                            all_derives = np.append(all_derives, latexify.parse_word(twin_axis))
                        if self.backend == "matplotlib":
                            overlay = hv.NdOverlay(
                                {derivative: hv.Scatter((np.NaN, np.NaN)).opts(opts.Scatter(s=50, color=cmap_values[i]))
                                for i, derivative in enumerate(all_derives)}
                            )
                            deriv_plot = df_tmp.hvplot.scatter(
                                x=x_axis,
                                y=deriv_col_name,
                                by="Derivatives",
                                title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
                                label=None,
                                datashade=datashade,
                                alpha=alpha[1],
                                legend=False,
                                cmap=cmap_values
                            ).opts(opts.Scatter(s=8)).opts(aspect=3.2)
                            deriv_plot = (deriv_plot * overlay)
                        else:
                            deriv_plot = df_tmp.hvplot.scatter(
                                x=x_axis,
                                y=deriv_col_name,
                                by="Derivatives",
                                title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
                                label=None,
                                datashade=datashade,
                                alpha=alpha[1],
                                cmap=cmap_values
                            ).opts(opts.Scatter(size=2))

                    else:
                        cmap = {}
                        for i, d in enumerate(all_derives):
                            cmap[d] = matplotlib.colors.to_hex(colors(i)[0:-1])
                        deriv_plot = df_tmp.hvplot.scatter(
                            x=x_axis,
                            y=deriv_col_name,
                            by="Derivatives",
                            title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
                            label=None,
                            datashade=datashade,
                            cmap=cmap
                        ).opts(aspect=3.2)

                        if self.backend == "bokeh":
                            points_deriv = hv.Points(
                                ([np.NaN for i in range(len(all_derives))],
                                 [np.NaN for i in range(len(all_derives))],
                                 list(cmap.keys())), vdims="Derivatives")
                            legend_deriv = points_deriv.options(
                                color="Derivatives",
                                cmap=cmap,
                                show_legend=True)
                            deriv_plot = (legend_deriv * deriv_plot).opts(aspect=3.2)

                t2 = timer()
                print("Setting up lower plot done in {} s".format(t2-t))
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
                    del_x = (this_max_x-this_min_x)*0.1
                    del_y = (max_y-min_y)
                    if min_y == max_y:
                        min_y = max_y - 1
                        max_y = max_y + 1
                        del_y = 2

                    min_y_deriv = df_tmp["Derivative Ratio"].min()
                    max_y_deriv = df_tmp["Derivative Ratio"].max()
                    del_y_deriv = (max_y_deriv-min_y_deriv)
                    if min_y_deriv == max_y_deriv:
                        min_y_deriv = max_y_deriv - 1
                        max_y_deriv = max_y_deriv + 1
                        del_y_deriv = 2

                    if by is not None:
                        df_tmp_group = df_group.groupby(by)
                    else:
                        df_tmp_group = df_group

                    for col in vertical_mark:
                        for v in vertical_mark[col]:
                            # Works well as long as dataframe is smaller than ~ 30k elements
                            df_sort = df_tmp_group.apply(lambda x: x.iloc[np.argmin(np.abs(x[col]-v))] )
                            col_values = np.sort(df_sort[col].values)
                            for index, row in df_sort.iterrows():
                                if np.abs(row[col]-v) > 1.0:
                                    break
                                # The one for the plot above
                                text = hv.Text(
                                    row[x_axis]+del_x,
                                    max_y-del_y*0.1,
                                    col + "=" + str(v) + latexify.get_unit(col),
                                    fontsize=24)
                                mark = hv.VLine(
                                        x=row[x_axis],
                                        label=col + "=" + str(v)
                                    ).opts(
                                        color="black",
                                        ylim=(min_y-del_y*0.1,
                                              max_y+del_y*0.1)) * text
                                param_plot = param_plot * mark
                                # For the derivatives
                                text = hv.Text(
                                    row[x_axis]+del_x,
                                    max_y_deriv-del_y_deriv*0.1,
                                    col + "=" + str(v) + latexify.get_unit(col),
                                    fontsize=24)
                                mark = hv.VLine(
                                        x=row[x_axis],
                                        label=col + "=" + str(v)
                                    ).opts(
                                        color="black",
                                        ylim=(min_y_deriv-del_y_deriv*0.1,
                                              max_y_deriv+del_y_deriv*0.1)) * text
                                deriv_plot = deriv_plot * mark
                    t2 = timer()
                    print("Creating marks done in {} s".format(t2-t))
                    t = timer()

                # Add second y-axis if given
                if twin is not None:
                    param_plot = param_plot * twin
                    deriv_plot = deriv_plot * twin

                if hist[0] and not hist[1]:
                    layout = param_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)
                else:
                    layout = param_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)

                opts_arg = {} # Currently empty. Maybe useful in further iterations

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
                if  self.backend == "bokeh":
                    layout_kwargs["width"] = width
                    layout_kwargs["height"] = height
                else:
                    layout_kwargs["fig_size"] = height/2

               # Matplotlib uses a horrible colormap as default...
                curve_kwargs = kwargs.copy()
                if percentile is not None:
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

                if errorband:
                    both_plots = layout.opts(
                        opts.Area(
                            xticks=20,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **area_kwargs),
                        opts.Scatter(
                            xticks=20,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **scatter_kwargs),
                        opts.Layout(**layout_kwargs)
                    ).cols(1)
                elif percentile is not None:
                    both_plots = layout.opts(
                        opts.Area(
                            xticks=20,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **area_kwargs),
                        opts.Scatter(
                            xticks=20,
                            xaxis="bottom",
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True,
                            **scatter_kwargs),
                        opts.Curve(**curve_kwargs),
                        opts.Layout(**layout_kwargs)
                    ).cols(1)
                else:
                    if not hexbin[0] and not hexbin[1]:
                        both_plots = layout.opts(
                            opts.Scatter(
                                xticks=20,
                                xaxis="bottom",
                                fontsize=self.font_dic,
                                show_grid=True,
                                show_legend=True,
                                **scatter_kwargs),
                            opts.Layout(**layout_kwargs)
                        ).cols(1)
                    elif hexbin[0] and hexbin[1]:
                        both_plots = layout.opts(
                            opts.HexTiles(**opts_arg),
                            opts.Layout(**layout_kwargs)
                        ).cols(1)
                    else:
                        both_plots = layout.opts(
                            opts.Scatter(
                                xticks=20,
                                xaxis="bottom",
                                fontsize=self.font_dic,
                                show_grid=True,
                                show_legend=True,
                                **scatter_kwargs),
                            opts.HexTiles(**opts_arg),
                            opts.Layout(**layout_kwargs)
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
                print("Setting options done in {} s".format(t2-t))
                t = timer()

                renderer = hv.Store.renderers[self.backend].instance(
                    fig='png', dpi=300)
                latexify.set_size(True)

                i = 0
                if prefix is None:
                    prefix = "plt_line_"
                    if scatter:
                        prefix = "plt_scatter_"
                save = (plot_path + prefix + x_axis + "_" + out_par
                        + "_" + "{:03d}".format(i))
                while os.path.isfile(save + ".png"):
                    i = i+1
                    save = (plot_path + prefix + x_axis + "_" + out_par
                            + "_" + "{:03d}".format(i))
                t2 = timer()
                print("Everything else done in {} s".format(t2-t), flush=True)
                if self.backend == "bokeh":
                    print("Plotting", flush=True)
                    hvplot.show(both_plots)
                    self.plots.append(both_plots)
                    t2 = timer()
                    print("Plotting done in {}s".format(t2-t), flush=True)
                else:
                    filetype = ".png"
                    if datashade:
                        filetype = ".html"
                        self.plots.append(both_plots)
                    print("Saving to " + save + filetype, flush=True)
                    renderer.save(both_plots, save)
                    t2 = timer()
                    try:
                        from IPython.display import Image, display
                        display(Image(save + filetype, width=width))
                    except:
                        pass
                    print("Saving done in {}s".format(t2-t), flush=True)
            if sort:
                i = 0
                if n_plots is None:
                    n_plots = 9999999

                while len(sorted_tuples) > 0 and i < n_plots:
                    p, v = sorted_tuples.pop()
                    in_params_2 = [p]
                    v_max = v
                    if log[1]:
                        while (len(sorted_tuples) > 0 and np.abs(v_max) - np.abs(sorted_tuples[-1][1]) < 15):
                            p, v = sorted_tuples.pop()
                            in_params_2.append(p)
                            if len(in_params_2) == max_per_deriv:
                                break
                    else:
                        while (len(sorted_tuples) > 0
                               and np.abs(sorted_tuples[-1][1]/v_max) > 0.1):
                            p, v = sorted_tuples.pop()
                            in_params_2.append(p)
                            if len(in_params_2) == max_per_deriv:
                                break
                    plot_helper(df_tmp_out, in_params=in_params_2,
                                prefix=prefix, **kwargs)
                    i += 1
            else:
                # Plot for every input parameter
                for in_param in in_params:
                    plot_helper(df_tmp_out, in_params=[in_param],
                                prefix=prefix + "_" + in_param, **kwargs)

    def plot_grid_one_param(self, out_param, y_axis, x_axis="timestep",
        twin_axis=None, by=None, hue="type", col_wrap=4, trajectories=None,
        width=1280, height=800, log=[False, False],
        vertical_mark=None, datashade=False, prefix=None, alpha=1,
        plot_path="pics/", yticks=10, decimals=-3, rolling_avg=20, **kwargs):
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
        datashade : bool
            Use datashade to plot results.
            Not yet implemented/tested.
        prefix : string
            Prefix to add to the filename.
        alpha : float
            Alpha value to use on the plot aligned to the left y-axis.
        plot_path : string
            Path to the folder to save the image.
        ytick : int
            Number of ticks on all y-axis.
        decimals : int
            For rounding the right y-axis, if provided by "twin_axis".
            From numpy: Number of decimal places to round to (default: 0).
            If decimals is negative, it specifies the number of
            positions to the left of the decimal point.
        rolling_avg : int
            Number of consecutive rows for a rolling average
            for the left y-axis.
        """
        import hvplot.pandas
        from holoviews import opts
        import holoviews as hv
        import dask.array as da
        import math

        df = self.cache
        hv.extension(self.backend)

        aspect = width/height
        layout_kwargs = {}
        if  self.backend == "bokeh":
            layout_kwargs["width"] = width
            layout_kwargs["height"] = height
        else:
            layout_kwargs["fig_size"] = height/2

        plot_list = []
        df_tmp_out = df.loc[df["Output Parameter"] == out_param]
        df_tmp_out["qi"] = da.fabs(df_tmp_out["qi"]) # sign error?
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
            tick_size = (upper_y - lower_y)/yticks
            def twinx(plot, element):
                ax = plot.handles["axis"]
                twinax = ax.twinx()
                twinax.set_ylim((lower_y - tick_size * 0.75,
                                 upper_y + tick_size * 0.75))
                twinax.set_yticks(np.arange(lower_y, upper_y + tick_size*0.75, tick_size))
                twinax.set_ylabel(
                    latexify.parse_word(twin_axis)
                    + " " + latexify.get_unit(twin_axis, brackets=True))
                plot.handles["axis"] = twinax
            if log[1]:
                df_tmp_out[twin_axisy] = da.log(da.fabs(df_tmp_out[twin_axis]))
            twin = df_tmp_out.hvplot.scatter(
                x=x_axis,
                y=twin_axis,
                color="gray",
                datashade=datashade,
                legend=False,
            ).opts(initial_hooks=[twinx], apply_ranges=False).opts(opts.Scatter(s=8))

        add_cols = []
        if vertical_mark is not None:
            add_cols = list(vertical_mark.keys())
        for y in y_axis:
            marks = None
            if by is not None:
                if y == x_axis:
                    df_group = df_tmp_out[[x_axis, by] + add_cols]
                else:
                    df_group = df_tmp_out[[x_axis, y, by] + add_cols]
                types = df_group[by].unique()
                types = np.sort(types[::-1])
                if rolling_avg is not None:
                    for ty in types:
                        df_roll = df_group.loc[df_group[by] == ty]
                        series = df_roll[y].rolling(
                            rolling_avg, min_periods=1).mean()
                        series.index = df_group.loc[df_group[by] == ty].index
                        df_group.update(series)
                if log[0]:
                    # Apply should be more expensive
                    df_group[y] = da.log(da.fabs(df_group[y]))

                min_x = df_group[x_axis].min()
                max_x = df_group[x_axis].max()
                min_y = df_group[y].min()
                max_y = df_group[y].max()
                del_x = (max_x-min_x)*0.1
                del_y = (max_y-min_y)
                if min_y == max_y:
                    min_y = max_y - 1
                    max_y = max_y + 1
                    del_y = 2

                if vertical_mark is not None:
                    if by is not None:
                        df_tmp = df_group.groupby(by)
                    else:
                        df_tmp = df_group

                    for col in vertical_mark:
                        for v in vertical_mark[col]:
                            # Works well as long as dataframe is smaller than ~ 30k elements
                            df_sort = df_tmp.apply(lambda x: x.iloc[np.argmin(np.abs(x[col]-v))] )
                            col_values = np.sort(df_sort[col].values)
                            for index, row in df_sort.iterrows():
                                if np.abs(row[col]-v) > 1.0:
                                    break
                                text = hv.Text(row[x_axis]+del_x, max_y-del_y*0.1, col + "=" + str(v) + latexify.get_unit(col),
                                    fontsize=24)
                                if marks is None:
                                    marks = hv.VLine(
                                        x=row[x_axis],
                                        label=col + "=" + str(v)
                                    ).opts(
                                        color="black",
                                        ylim=(min_y-del_y*0.1, max_y+del_y*0.1)) * text
                                else:
                                    marks = marks * hv.VLine(
                                        x=row[x_axis],
                                        label=col + "=" + str(v)
                                    ).opts(
                                        color="black",
                                        ylim=(min_y-del_y*0.1, max_y+del_y*0.1)) * text
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
                        if len(plot_list) == 0:
                            overlay = hv.NdOverlay(
                                {types[i]: hv.Scatter((np.NaN, np.NaN)).opts(opts.Scatter(s=50, color=cmap_values[i]))
                                for i in range(len(types)) }
                            )
                            plot_list.append(df_group.hvplot.scatter(
                                x=x_axis,
                                y=y,
                                by=by,
                                title="Values of of {}".format(latexify.parse_word(y)),
                                color=cmap_values,
                                label=None,
                                datashade=datashade,
                                alpha=alpha,
                                legend=False
                            ).opts(opts.Scatter(s=8)).opts(
                                aspect=aspect,
                                ylim=(min_y-del_y*lim_multiplier, max_y+del_y*lim_multiplier),
                                xlim=(min_x, max_x),
                                yticks=yticks) * overlay)
                        else:
                            plot_list.append(df_group.hvplot.scatter(
                                x=x_axis,
                                y=y,
                                by=by,
                                title="Values of of {}".format(latexify.parse_word(y)),
                                color=cmap_values,
                                label=None,
                                datashade=datashade,
                                alpha=alpha,
                                legend=False
                            ).opts(opts.Scatter(s=8)).opts(
                                aspect=aspect,
                                ylim=(min_y-del_y*lim_multiplier, max_y+del_y*lim_multiplier),
                                xlim=(min_x, max_x),
                                yticks=yticks))

                        if marks is not None:
                            plot_list[-1] = (plot_list[-1] * marks)
                        if twin is not None:
                            plot_list[-1] = (plot_list[-1] * twin)
                    else:
                        param_plot = df_group[df_group[y] != 0].hvplot.scatter(
                            x=x_axis,
                            y=y,
                            by=by,
                            title="Values of of {}".format(latexify.parse_word(y)),
                            label=None,
                            datashade=datashade,
                            alpha=alpha
                        ).opts(opts.Scatter(size=2), **layout_kwargs)

                else:
                    cmap = {}
                    for ty in types:
                        cmap[ty] = self.cmap[ty]
                    param_plot = df_group[df_group[y] != 0].hvplot.scatter(
                        x=x_axis,
                        y=y,
                        by=by,
                        title="Values of of {}".format(latexify.parse_word(y)),
                        cmap=cmap,
                        label=None,
                        datashade=datashade
                    ).opts(aspect=aspect)
                    if self.backend == "bokeh":
                        points = hv.Points(
                            ([np.NaN for i in range(len(list(cmap.keys())))],
                             [np.NaN for i in range(len(list(cmap.keys())))],
                             list(cmap.keys())), vdims=by)
                        legend = points.options(
                            color=by,
                            cmap=cmap,
                            show_legend=True)
                        param_plot = (legend * param_plot).opts(aspect=aspect, **layout_kwargs)
            else:
                if y == x_axis:
                    df_group = df_tmp_out[[x_axis, "trajectory"]]
                else:
                    df_group = df_tmp_out[[x_axis, y, "trajectory"]]

                if rolling_avg is not None:
                    series = df_group[y].rolling(
                        rolling_avg, min_periods=1).mean()
                    df_group.update(series)

                if log:
                    # Apply should be more expensive
                    df_group[y] = da.log(da.fabs(df_group[y]))
                param_plot = df_group.hvplot.scatter(
                    x=x_axis,
                    y=y,
                    title="Values of of {}".format(latexify.parse_word(y)),
                    label=None,
                    datashade=datashade
                ).opts(**layout_kwargs)
        all_plots = plot_list[0] + plot_list[1]
        for i in range(len(plot_list)-2):
            all_plots = all_plots + plot_list[i+2]
        layout_kwargs = {}
        if self.backend == "bokeh":
            layout_kwargs["width"] = width
            layout_kwargs["height"] = height
        else:
            layout_kwargs["fig_size"] = height/2

        scatter_kwargs = {}

        for k in kwargs:
            scatter_kwargs[k] = kwargs[k]

        final_plots = all_plots.opts(
            opts.Scatter(
                xticks=20,
                xaxis="bottom",
                fontsize=self.font_dic,
                show_grid=True,
                show_legend=True,
                **scatter_kwargs),
            opts.Layout(**layout_kwargs)
        ).cols(col_wrap)
        if self.backend == "matplotlib":
            final_plots = final_plots.opts(sublabel_format="", tight=False)

        renderer = hv.Store.renderers[self.backend].instance(
                    fig='png', dpi=300)
        latexify.set_size(True)

        i = 0
        if prefix is None:
            prefix = "plt_grid_"

        save = (plot_path + prefix + x_axis + "_" + "{:03d}".format(i))
        while os.path.isfile(save + ".png"):
            i = i+1
            save = (plot_path + prefix + x_axis + "_" + "{:03d}".format(i))

        if self.backend == "bokeh":
            hvplot.show(final_plots)
            self.plots.append(final_plots)
        else:
            filetype = ".png"
            if datashade:
                filetype = ".html"
                self.plots.append(final_plots)
            print("Saving to " + save + filetype, flush=True)
            renderer.save(final_plots, save)
            try:
                from IPython.display import Image, display
                display(Image(save + filetype, width=width))
            except:
                pass
