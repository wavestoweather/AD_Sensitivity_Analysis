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

from sklearn.cluster import MiniBatchKMeans, SpectralClustering, DBSCAN
from sklearn.mixture import BayesianGaussianMixture
from sklearn.metrics import adjusted_rand_score

from dask_ml.cluster import KMeans

from itertools import repeat
from itertools import product
from progressbar import progressbar as pb

try:
    import dask_loader
except:
    import scripts.dask_loader as dask_loader

import dask.dataframe as pd
from bokeh.io import curdoc

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
    """
    data = {}
    n_timesteps = 0
    cluster_names = {}
    threads = 0
    font_dic = {}
    backend = "matplotlib"
    widget = None

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

    def to_parquet(self, f_name):
        pd.to_parquet(self.data, f_name + ".parquet")

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

    def plot_same_orders(self, out_params=None, mapped=True, scatter=False,
                         in_params=None, x_axis="timestep", n_plots=5, frac=None,
                         verbose=False, **kwargs):
        """
        For each out_param, plot multiple figures with derivatives, where
        each figure shows derivatives of the same order.

        Parameters
        ----------
        out_params : list of string
            List of keys to plot the derivatives for. If None is given, all
            available parameters will be plotted.
        mapped : boolean
            If true: plot the region, where "MAP" is true, ie where the wcb
            criterion is fullfilled.
        scatter : boolean
            Plot a scatter plot or a line plot.
        in_params : list of string
            Plot only the derivatives with respect to those in this list.
        x_axis : string
            The column to use as x-axis. Can be either "timestep" or
            "out_param" or an output parameter.
        n_plots : int
            Plot only that many plots. If None is given, plot every possible
            plot.
        kwargs : dict
            Keyword arguments are passed down to seaborn.scatterplott() or
            seaborn.lineplot() for the derivative plots.
        """
        if frac is not None:
            if verbose:
                print("Sample")
            df = self.data.sample(frac=frac, replace=False, random_state=42).compute()
        else:
            df = self.data.compute()
        if verbose:
            print("Sampling done")
        if out_params is not None:
            df = df.loc[df["Output Parameter"].isin(out_params)]


        if in_params is None:
            in_params = []
            head = df.columns
            for h in head:
                if h[0] == 'd':
                    in_params.append(h)
            if verbose:
                print("Got in_params: {}".format(in_params))
        if verbose:
            print("Sorting the derivatives")
        # Sort the derivatives
        sorted_tuples = []
        for in_p in in_params:
            if verbose:
                print("Look at {}".format(in_p))
            value = np.abs(df[in_p].min())
            if verbose:
                print("Got value {}".format(value))
            max_val = np.abs(df[in_p].max())
            if verbose:
                print("Or is it {}?".format(max_val))
            if max_val > value:
                value = max_val
            if value != 0 and not np.isnan(value):
                sorted_tuples.append((in_p, value))
        sorted_tuples.sort(key=lambda tup: tup[1])

        if verbose:
            print("Sorting finished")

        def plot_helper(df, in_params, out_param=None, **kwargs):
            min_time = df[x_axis].min()
            max_time = df[x_axis].max()
            dt1 = (max_time - min_time) / 20
            dt = (max_time - min_time + dt1) / 20

            x_ticks = np.arange(min_time, max_time + dt1, dt)

            if verbose:
                print("Melt dataframe")

            _, ax = plt.subplots()
            df_tmp = df[in_params+[x_axis]]
            df_tmp = df_tmp.melt(x_axis, var_name="Derivatives",
                                    value_name="Derivative Ratio")
            if verbose:
                print("Plot on axis")

            if scatter:
                ax = sns.scatterplot(x=x_axis, y="Derivative Ratio",
                                        data=df_tmp, hue="Derivatives", ax=ax, linewidth=0,
                                        **kwargs)
            else:
                ax = sns.lineplot(x=x_axis, y="Derivative Ratio",
                                    data=df_tmp, hue="Derivatives", ax=ax,
                                    **kwargs)
            if verbose:
                print("Plotting finished")
            if out_param is None:
                out_param = "all Parameters"
            ax.set_title("Deriv. Ratio of {}".format(latexify.parse_word(out_param)))
            ax.set_xticks(x_ticks)

            # Change labels to latex versions
            legend = ax.get_legend()
            _, labels = ax.get_legend_handles_labels()
            for t, old in zip(legend.texts, labels):
                t.set_text(latexify.parse_word(old))
            ax.set_ylabel("Derivative ratio")
            plt.ticklabel_format(style="scientific", axis="y", scilimits=(0,0))

            # Plot the area that had been flagged
            if mapped:
                df_mapped = df[df.MAP == True]
                if not df_mapped.empty:
                    xmin = df_mapped[x_axis].min()
                    xmax = df_mapped[x_axis].max()
                    plt.axvspan(xmin=xmin, xmax=xmax,
                        facecolor="khaki", alpha=0.3)

            # Change the limits for the y-axis because sometimes that
            # can be off and it is hard to see anything.
            min_y = df_tmp["Derivative Ratio"].min()
            max_y = df_tmp["Derivative Ratio"].max()
            plt.ylim(min_y, max_y)
            # The same is true for the x-axis
            plt.xlim(x_ticks[0], x_ticks[-1])
            i = 0
            prefix = "line_"
            if scatter:
                prefix = "scatter_"
            save = ("pics/" + prefix + x_axis + "_" + out_param
                    + "_" + "{:03d}".format(i) + ".png")
            while os.path.isfile(save):
                i = i+1
                save = ("pics/" + prefix + x_axis + "_" + out_param
                        + "_" + "{:03d}".format(i) + ".png")

            print("Saving to " + save)
            success = False
            for i in range(100):
                try:
                    plt.savefig(save, dpi=300)
                    plt.show()
                    plt.close()
                except:
                    continue
                success = True
                break
            if not success:
                print("Saving image {} failed.".format(save))

        # Plot them
        i = 0
        while len(sorted_tuples) > 0 and i < n_plots:
            p, v = sorted_tuples.pop()
            in_params_2 = [p]
            while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                and np.abs(v/sorted_tuples[-1][1]) < 10):
                p, v = sorted_tuples.pop()
                in_params_2.append(p)
            if verbose:
                print("Plot {}".format(in_params_2))
            plot_helper(df, in_params=in_params_2, **kwargs)
            i += 1

    def plot_mapped(self, out_params=None, in_params=None, kind="violin",
                    x_label="WCB criterion", n_plots=None, **kwargs):
        """
        Create a scatter plot that uses the x-axis for the boolean MAP.
        Input parameter values of the same order are plotted on the same
        plot.

        Parameters
        ----------
        out_params : list of string
            List of keys to plot the derivatives for. If None is given, all
            available parameters will be plotted.
        in_params : list of strings
            List of derivatives to plot for. If none is given, a plot
            will be created for every input parameter.
        kind : String
            Can bei either of "swarm", "violin" (both recommended, but swarm
            may take long for many datapoints) or
            "point", "bar", "strip", "box", or "boxen". Is passed
            to seaborn.catplot(..).
        x_label : String
            This will be plotted on the x-axis.
        n_plots : int
            Plot only that many plots. If None is given, plot every possible
            plot.
        kwargs : dict
            Keyword arguments are passed down to seaborn.catplot(..) for
            the derivative plots.
        """
        if out_params is not None:
            df = self.data.loc[self.data["Output Parameter"].isin(out_params)]
        else:
            df = self.data

        def plot_helper(df, in_params, out_param, **kwargs):
            df_tmp = df[in_params+["MAP"]]

            df_tmp = df_tmp.melt("MAP", var_name="Derivatives",
                                value_name="Derivative Ratio").compute()
            df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)

            _, ax = plt.subplots()
            # Just a workaround. One could use sns.catplot(kind=kind) but
            # it does not work with plt.subplots() and it does not always
            # take any preceeding settings for matplotlib into account
            if kind == "violin":
                sns.violinplot(x="MAP", y="Derivative Ratio", data=df_tmp,
                            hue="Derivatives", **kwargs)
            elif kind == "swarm":
                sns.swarmplot(x="MAP", y="Derivative Ratio", data=df_tmp,
                            hue="Derivatives", **kwargs)
            elif kind == "strip":
                sns.stripplot(x="MAP", y="Derivative Ratio", data=df_tmp,
                            hue="Derivatives", **kwargs)
            elif kind == "box":
                sns.boxplot(x="MAP", y="Derivative Ratio", data=df_tmp,
                            hue="Derivatives", **kwargs)
            elif kind == "boxen":
                sns.boxenplot(x="MAP", y="Derivative Ratio", data=df_tmp,
                            hue="Derivatives", **kwargs)
            elif kind == "point":
                sns.pointplot(x="MAP", y="Derivative Ratio", data=df_tmp,
                            hue="Derivatives", **kwargs)
            elif kind == "bar":
                sns.barplot(x="MAP", y="Derivative Ratio", data=df_tmp,
                            hue="Derivatives", **kwargs)
            elif kind == "count":
                sns.countplot(x="MAP", y="Derivative Ratio", data=df_tmp,
                            hue="Derivatives", **kwargs)
            else:
                print("No such kind: {}".format(kind))
                return

            ax.set_title("Deriv. Ratio of {}".format([latexify.parse_word(p) for p in out_param]))
            ax.set_ylabel("Derivative ratio")
            ax.set_xlabel(x_label)
            plt.ticklabel_format(style="scientific", axis="y", scilimits=(0,0))

            # Change the limits for the y-axis because sometimes that
            # can be off and it is hard to see anything.
            min_y = df_tmp["Derivative Ratio"].min()
            max_y = df_tmp["Derivative Ratio"].max()
            plt.ylim(min_y, max_y)

            i = 0
            save = ("pics/" + kind + "_MAP_ " + out_param[0]
                    + "_" + "{:03d}".format(i) + ".png")
            while os.path.isfile(save):
                i = i+1
                save = ("pics/" + kind + "_MAP_ " + out_param[0]
                        + "_" + "{:03d}".format(i) + ".png")

            print("Saving to " + save)
            plt.savefig(save, dpi=300)
            plt.show()
            plt.close()

        if in_params is None:
            in_params = []
            head = df.columns
            for h in head:
                if h[0] == 'd':
                    in_params.append(h)

        # Sort the derivatives
        sorted_tuples = []
        for in_p in in_params:
            value = np.abs(df[in_p].min().compute())
            max_val = np.abs(df[in_p].max().compute())
            if max_val > value:
                value = max_val
            if value != 0 and not np.isnan(value):
                sorted_tuples.append((in_p, value))
        sorted_tuples.sort(key=lambda tup: tup[1])

        # Plot them
        i = 0
        while len(sorted_tuples) > 0 and i < n_plots:
            p, v = sorted_tuples.pop()
            in_params_2 = [p]
            while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                and np.abs(v/sorted_tuples[-1][1]) < 10):
                p, v = sorted_tuples.pop()
                in_params_2.append(p)
            plot_helper(df, in_params=in_params_2, out_param=out_params, **kwargs)
            i += 1

    def cluster(self, k, method, out_params=None, features=None,
                new_col="cluster", truth=None, thresh=0.10):
        """
        Cluster the data where "features" are columns, "truth" can be a column
        that can be used to identify a ground truth or a list of assignments.
        If truth is given, purity and adjusted RAND index are calculated.
        "method" can be used to choose different clustering methods. Discards
        all columns that consist of at least one NaN.

        Parameters
        ----------
        k : int
            Number of clusters to cluster for. Is ignored for certain cluster
            methods.
        method : string
            Choose the clustering method. Options are "kmeans", "spectral",
            "gaussian" or "dbscan".
        out_params : list of strings
            The parameter for which the derivatives shall be considered
            for. If none is given, all output parameters will be considered.
        features : list of strings
            Columns of the data that shall be used. If it is None, all
            derivatives will be used.
        new_col : string
            Name of the column where the cluster assignment shall be stored.
        truth : String or list of int
            Either a column of the data or a list of assignments.
        thresh : float
            Threshold for dropping columns. if more than this amount of values
            is NaN, drop that column.

        Returns
        -------
        Array of float, optionally float
            Cluster centers if "kmeans" is selected,
            adjusted RAND index if truth is given.
        """

        if features is None:
            features_tmp = self.data.columns.values
            features = []
            for f in features_tmp:
                if f[0] == "d":
                    features.append(f)

        cols_to_delete = self.data.columns[self.data.isnull().sum()/len(self.data)  > thresh]
        for col in cols_to_delete:
            try:
                features.remove(col)
            except:
                pass
        if features == []:
            print("No derivatives in this dataset. Returning.")
            if method == "gaussian":
                return None, None, None
            return None, None

        tmp_df = self.data.dropna(subset=features)
        if len(tmp_df.index) == 0:
            print("Empty dataframe")
            return None, None
        print("Clustering starts")
        clustering = KMeans(n_clusters=k).fit(tmp_df[features])
        print("Clustering finished")
        if truth is None:
            return clustering.cluster_centers_, None
        print("Calculating RAND index")
        if isinstance(truth, str):
            adj_index = adjusted_rand_score(self.data[truth], clustering.labels_)
        else:
            adj_index = adjusted_rand_score(truth, clustering.labels_)
        print("Done")
        return clustering.cluster_centers_, adj_index

    def get_sorting(self, out_params=None, in_params=None, mapped=False, verbose=False):
        """
        Get a list of tuples with the derivatives order.

        Parameters
        ----------
        out_params : list of string
            List of keys to get the derivatives for. If None is given, all
            available parameters will be considered.
        in_params : list of string
            Get only the derivatives with respect to those in this list.
        """
        if out_params is not None:
            df = self.data.loc[self.data["Output Parameter"].isin(out_params)]
        else:
            df = self.data

        if in_params is None:
            in_params = []
            head = df.columns
            for h in head:
                if h[0] == 'd':
                    in_params.append(h)
            if verbose:
                print("Got in_params: {}".format(in_params))
        if verbose:
            print("Sorting the derivatives")

        if mapped:
            df = df.loc[df.MAP == True]
            df = df.compute()
            # Sort the derivatives
            sorted_tuples = []
            for in_p in in_params:
                if verbose:
                    print("Look at {}".format(in_p))
                value = np.abs(df[in_p].min())
                if verbose:
                    print("Got value {}".format(value))
                max_val = np.abs(df[in_p].max())
                if verbose:
                    print("Or is it {}?".format(max_val))
                if max_val > value:
                    value = max_val
                if value != 0 and not np.isnan(value):
                    sorted_tuples.append((in_p, value))
            sorted_tuples.sort(key=lambda tup: tup[1])
            return sorted_tuples
        else:
            # Sort the derivatives
            sorted_tuples = []
            for in_p in in_params:
                if verbose:
                    print("Look at {}".format(in_p))
                value = np.abs(df[in_p].min().compute())
                if verbose:
                    print("Got value {}".format(value))
                max_val = np.abs(df[in_p].max().compute())
                if verbose:
                    print("Or is it {}?".format(max_val))
                if max_val > value:
                    value = max_val
                if value != 0 and not np.isnan(value):
                    sorted_tuples.append((in_p, value))
            sorted_tuples.sort(key=lambda tup: tup[1])
            return sorted_tuples

    def plot_two_ds(self, in_params, out_params, x_axis="timestep", mapped=True,
        trajectories=None, scatter=False, n_plots=None, percentile=None,
        frac=None, min_x=None, max_x=None, nth=None, hist=[False, False],
        hexbin=[False, False], log=[False, False], sort=True,
        scatter_deriv=False, line_deriv=False, prefix=None, c=False, compute=False,
        errorband=False, bins=50, plot_path="pics/", fig_type='svg', **kwargs):
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
        mapped : boolean
            If true: plot the region, where "MAP" is true, ie where the wcb
            criterion is fullfilled.
        trajectories : int
            The index of the trajectories to plot. If None is given, all
            trajectories will be plotted.
        scatter : boolean
            Plot a scatter plot or a line plot.
        n_plots : int
            Plot only that many plots. If None is given, plot every possible
            plot.
        percentile : list of int
            Plot the given percentiles along with the spread if
            errorband is not given. Not implemented yet
        errorband : Bool
            Plot an errorband (spread) using the standard deviation and min/max values.
        frac : float
            Sample a given fraction of rows. Deactivates "nth".
        nth : int
            Sample every nth entry. Works fast with "timestep" as x-axis and
            a given min_x and max_x value. If x_axis is any other than
            "timestep", an errorband triggered by "percentile" may not
            be plotted.
        scatter_deriv : boolean
            Plot the derivative plot with scatterplot. If scatter_deriv and
            line_deriv are set to false, use the same as provided via scatter
            and percentile.
        line_deriv : boolean
            Plot the derivative plot with lineplot. If scatter_deriv and
            line_deriv are set to false, use the same as provided via scatter
            and percentile.
        prefix : string
            Prefix to add to the filename.
        sort : Bool
            If True, sort the derivatives and plot only those within the same 
            magnitude in one plot. If False, plot every derivative.
        kwargs : dict
            Keyword arguments are passed down matplotlib.
        """
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
        if mapped:
            df = df.loc[df.MAP == True]
        if compute:
            df = df.loc[df["Output Parameter"].isin(out_params)].compute()
        else:
            df = df.loc[df["Output Parameter"].isin(out_params)]
        import hvplot.dask # adds hvplot method to dask objects
        import hvplot.pandas
        import hvplot
        from holoviews import opts
        import holoviews as hv
        from timeit import default_timer as timer
        from holoviews.operation import histogram as hv_histo
        import pandas
        import dask.array as da


        hv.extension(self.backend)

        for out_par in out_params:
            df_tmp_out = df.loc[df["Output Parameter"] == out_par]
            t = timer()
            # Sort the derivatives
            if sort:
                sorted_tuples = []
                for in_p in in_params:
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
                # following https://holoviz.org/tutorial/Composing_Plots.html
                t = timer()
                df_tmp = df[in_params+[x_axis]]
                df_tmp = df_tmp.melt(x_axis, var_name="Derivatives",
                                    value_name="Derivative Ratio")
                t2 = timer()
                print("Melting done in {} s".format(t2-t), flush=True)
                t = timer()
                if log[1]:
                    df_tmp["Derivative Ratio"] = df_tmp["Derivative Ratio"].apply(lambda x: np.log(np.abs(x)))
                    # Remove zero entries (-inf)
                    df_tmp = df_tmp[~da.isinf(df_tmp["Derivative Ratio"])]
                df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)

                if percentile is not None:

                    if out_par == x_axis:
                        print("x-axis and y-axis are the same. Cannot plot that with percentiles!")
                        return
                    else:
                        df_group = df[[x_axis, out_par]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                        # Remove zero entries (-inf)
                        df_group = df_group[~da.isinf(df_group[out_par])]

                    # Group for min, max and percentiles
                    funcs = [np.min, np.max] + [lambda x, perc=perc: np.percentile(x, perc, axis=0) for perc in percentile]
                    df_min_max = df_group.groupby(x_axis).agg(funcs)[out_par]

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
                            ylabel=latexify.parse_word(out_par),
                            title="Values of of {}".format(latexify.parse_word(out_par)),
                            label="Spread",
                            color="grey")
                        * df_min_max.hvplot.line(
                            x=x_axis,
                            y=p_list,
                            ylabel=latexify.parse_word(out_par),
                            **kwargs)
                    )
                elif errorband:
                    df_group = df[[x_axis, out_par, "trajectory"]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                    df_min_max = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([np.min, np.max])
                    df_std = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([lambda x: -1*np.std(x)+np.mean(x), lambda x: np.std(x)+np.mean(x)])
                    df_mean = df_group.groupby(x_axis)[out_par].mean()

                    param_plot = (
                        df_min_max.hvplot.area(
                            x=x_axis,
                            y="amin",
                            y2="amax",
                            alpha=0.5,
                            value_label=latexify.parse_word(out_par),
                            label="Spread")
                            # color="grey")
                        * df_mean.hvplot()
                        * df_std.hvplot.area(
                            x=x_axis,
                            y="<lambda_0>",
                            y2="<lambda_1>",
                            alpha=0.3,
                            value_label=latexify.parse_word(out_par),
                            label="sd")
                            # color="grey")
                    )
                elif c:
                    if out_par == x_axis:
                        df_group = df[[x_axis, "trajectory"]]
                    else:
                        df_group = df[[x_axis, out_par, "trajectory"]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                    param_plot = df_group.hvplot.scatter(
                        x=x_axis,
                        y=out_par,
                        c="trajectory",
                        title="Values of of {}".format(latexify.parse_word(out_par)),
                        label=None
                        )
                elif hexbin[0]:
                    if out_par == x_axis:
                        df_group = df[[x_axis, "trajectory"]]
                    else:
                        df_group = df[[x_axis, out_par, "trajectory"]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                    param_plot = df_group.hvplot.hexbin(
                        x=x_axis,
                        y=out_par,
                        title="Values of of {}".format(latexify.parse_word(out_par)),
                        label=None,
                        clabel="Count",
                        cmap="viridis",
                        logz=True,
                        gridsize=100
                        )
                else:
                    if out_par == x_axis:
                        df_group = df[[x_axis, "trajectory"]]
                    else:
                        df_group = df[[x_axis, out_par, "trajectory"]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                    param_plot = df_group.hvplot.scatter(
                        x=x_axis,
                        y=out_par,
                        title="Values of of {}".format(latexify.parse_word(out_par)),
                        label=None
                    )

                layout_kwargs = {}
                if self.backend == "bokeh":
                    layout_kwargs["width"] = 1600
                    layout_kwargs["height"] = 500

                if hist[0]:
                    param_plot.opts(**layout_kwargs)
                    param_hist_plot = param_plot.hist(dimension=[x_axis, out_par], bins=bins, range=[0, 10000])
                

                if hexbin[1]:
                    deriv_plot = df_tmp.hvplot.hexbin(
                        x=x_axis,
                        y="Derivative Ratio",
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
#                     df_tmp = df_tmp.compute()
                    deriv_plot = df_tmp.hvplot.scatter(
                        x=x_axis,
                        y="Derivative Ratio",
                        by="Derivatives",
                        title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
                        label=None, 
#                         datashade=True # Useful for bokeh images, I guess
                    )
                  # Does not work: Rectangle object has no property dimension  
#                 if hist[1]:
#                     deriv_hist_plot = df_tmp.hist(dimension=[x_axis, "Derivative Ratio"], bins=bins)
                hist[1] = False
    
    
                if hist[0] and not hist[1]:
                    layout = param_hist_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)
                elif hist[0] and hist[1]:
                    layout = param_hist_plot.opts(**layout_kwargs) + deriv_hist_plot.opts(**layout_kwargs)
                elif not hist[0] and hist[1]:
                    layout = param_plot.opts(**layout_kwargs) + deriv_hist_plot.opts(**layout_kwargs)
                else:
                    layout = param_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)

                opts_arg = {}
                if self.backend == "matplotlib":
                    opts_arg["aspect"] = 3.2
                    opts_arg["fig_latex"] = False

                scatter_kwargs = opts_arg.copy()
                area_kwargs = opts_arg.copy()

                for k in kwargs:
                    scatter_kwargs[k] = kwargs[k]

                if self.backend == "bokeh":
                    scatter_kwargs["size"] = 5
                else:
                    scatter_kwargs["s"] = 8

                if self.backend == "matplotlib":
                    area_kwargs["edgecolor"] = None
                    area_kwargs["color"] = "black"
                    for k in kwargs:
                        area_kwargs[k] = kwargs[k]


                layout_kwargs = {}
                if self.backend == "matplotlib":
                    layout_kwargs["fig_size"] = 400

               # Matplotlib uses a horrible colormap as default...
                curve_kwargs = kwargs.copy()
                if percentile is not None:
                    if self.backend == "matplotlib":
                        if len(percentile) <= 10:
                            curve_kwargs["color"] = hv.Cycle("tab10")
                        elif len(percentile) <= 20:
                            curve_kwargs["color"] = hv.Cycle("tab20")
                        # More seems convoluted to me
#                         else:
#                             curve_kwargs["color"] = hv.Cycle("viridis")
                    else:
                        if len(percentile) <= 10:
                            curve_kwargs["color"] = hv.Cycle("Category10")
                        elif len(percentile) <= 20:
                            curve_kwargs["color"] = hv.Cycle("Category20")
#                         else:
#                             curve_kwargs["color"] = "viridis" # "colorcet"

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
                    if hexbin[0] and hexbin[1]:
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
                    param_opts = param_plot.opts(xaxis="bare", alpha=0.1)
                elif c:
                    param_opts = param_plot.opts(xaxis="bare", alpha=1.0)
                else:
                    param_opts = param_plot.opts(xaxis="bare")


                renderer = hv.Store.renderers[self.backend].instance(
                    fig='png', dpi=300)
                latexify.set_size(True)

                i = 0
                if prefix is None:
                    prefix = "plt_1line_"
                    if scatter:
                        prefix = "plt_1scatter_"
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
                    t2 = timer()
                    print("Plotting done in {}s".format(t2-t), flush=True)
                else:
                    print("Saving to " + save + ".png", flush=True)
                    renderer.save(both_plots, save)
                    t2 = timer()
                    try:
                        from IPython.display import Image, display
                        display(Image(save + ".png", width=1600))
                    except:
                        pass
                    print("Saving done in {}s".format(t2-t), flush=True)
            if sort:
                i = 0
                if n_plots is None:
                    n_plots = 9999999
    #             print("Creating {} plots".format(n_plots))
                while len(sorted_tuples) > 0 and i < n_plots:
                    p, v = sorted_tuples.pop()
                    in_params_2 = [p]
                    while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                        and np.abs(v/sorted_tuples[-1][1]) < 10):
                        p, v = sorted_tuples.pop()
                        in_params_2.append(p)
                    plot_helper(df_tmp_out, in_params=in_params_2, prefix=prefix, **kwargs)
                    i += 1
            else:
                # Plot for every input parameter
                for in_param in in_params:
                    plot_helper(df_tmp_out, in_params=[in_param], prefix=prefix + "_" + in_param, **kwargs)

    def plot_two_ds_2(self, in_params, out_params, x_axis="timestep", mapped=True,
        trajectories=None, scatter=False, n_plots=None, percentile=None,
        frac=None, min_x=None, max_x=None, nth=None, hist=False, hexbin=[False, False],
        log=[False, False],
        scatter_deriv=False, line_deriv=False, prefix=None, c=False,
        param_method="scatter",
        errorband=False, bins=50, plot_path="pics/", fig_type='svg', **kwargs):
        """
        Test case with datashader
        Plot two plots in two rows. At the top: Output parameter.
        At the bottom: Derivative with respect to that output parameter.
        Another plot is created for every gradient of different order.
        If multiple parameters are given, another plot is created for each
        output parameter. This function is mostly intended for multiple 
        trajectories.
        Does take long to plot each individual plot but may need less memory.

        Parameters
        ----------
        in_params : list of string
            Plot the derivatives with respect to those in this list.
        out_params : list of string
            List of keys to plot the derivatives for.
        x_axis : string
            The column to use as x-axis. Can be either "timestep" or
            an output parameter or a derivative.
        mapped : boolean
            If true: plot the region, where "MAP" is true, ie where the wcb
            criterion is fullfilled.
        trajectories : int
            The index of the trajectories to plot. If None is given, all
            trajectories will be plotted.
        hexbin : List of boolean
            First entry: Create a hexbin of simulation results.
            Second entry: Create a hexbin of derivatives.
        log : List of boolean
            Use log for y-axis for [0] simulation results, [1] derivatives.
        scatter : boolean
            Plot a scatter plot or a line plot.
        n_plots : int
            Plot only that many plots. If None is given, plot every possible
            plot.
        percentile : list of int
            Plot the given percentiles along with the spread if
            errorband is not given. Not implemented yet
        errorband : Bool
            Plot an errorband (spread) using the standard deviation and min/max values.
        frac : float
            Sample a given fraction of rows. Deactivates "nth".
        nth : int
            Sample every nth entry. Works fast with "timestep" as x-axis and
            a given min_x and max_x value. If x_axis is any other than
            "timestep", an errorband triggered by "percentile" may not
            be plotted.
        scatter_deriv : boolean
            Plot the derivative plot with scatterplot. If scatter_deriv and
            line_deriv are set to false, use the same as provided via scatter
            and percentile.
        line_deriv : boolean
            Plot the derivative plot with lineplot. If scatter_deriv and
            line_deriv are set to false, use the same as provided via scatter
            and percentile.
        prefix : string
            Prefix to add to the filename.
        kwargs : dict
            Keyword arguments are passed down matplotlib.
        """
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
        if mapped:
            df = df.loc[df.MAP == True]

        df = df.loc[df["Output Parameter"].isin(out_params)].compute()

        import hvplot.dask # adds hvplot method to dask objects
        import hvplot.pandas
        import hvplot
        from holoviews import opts
        import holoviews as hv
        from timeit import default_timer as timer
        from holoviews.operation import histogram as hv_histo
#             from holoviews.operation.datashader import aggregate
        import pandas
        from holoviews.operation.datashader import datashade


        hv.extension(self.backend)

        for out_par in out_params:
            df_tmp_out = df.loc[df["Output Parameter"] == out_par]

            # Sort the derivatives
            sorted_tuples = []
            for in_p in in_params:
                value = np.abs(df_tmp_out[in_p].min())
                max_val = np.abs(df_tmp_out[in_p].max())
                if max_val > value:
                    value = max_val
                if value != 0 and not np.isnan(value):
                    sorted_tuples.append((in_p, value))
            sorted_tuples.sort(key=lambda tup: tup[1])

            def plot_helper(df, in_params, prefix, **kwargs):
                # following https://holoviz.org/tutorial/Composing_Plots.html
                t = timer()
                df_tmp = df[in_params+[x_axis]]
                df_tmp = df_tmp.melt(x_axis, var_name="Derivatives",
                                    value_name="Derivative Ratio")
                if log[1]:
                    df_tmp["Derivative Ratio"] = df_tmp["Derivative Ratio"].apply(lambda x: np.log(np.abs(x)))
                df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)
            
                # Output simulation plot
                if percentile is not None:
                    if out_par == x_axis:
                        print("x-axis and y-axis are the same. Cannot plot that with percentiles!")
                        return
                    else:
                        df_group = df[[x_axis, out_par]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))

                    # Group for min, max and percentiles
                    funcs = [np.min, np.max] + [lambda x, perc=perc: np.percentile(x, perc, axis=0) for perc in percentile]
                    df_min_max = df_group.groupby(x_axis).agg(funcs)[out_par]

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
                        datashade(df_min_max.hvplot.area(
                            x=x_axis,
                            y="Min",
                            y2="Max",
                            alpha=0.5,
                            value_label=latexify.parse_word(out_par),
                            title="Values of of {}".format(latexify.parse_word(out_par)),
                            label="Spread",
                            color="grey"))
                        * datashade(df_min_max.hvplot.line(
                            x=x_axis,
                            y=p_list,
                            value_label=latexify.parse_word(out_par),
                            **kwargs))
                    )
                elif errorband:
                    df_group = df[[x_axis, out_par, "trajectory"]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                    df_min_max = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([np.min, np.max])
                    df_std = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([lambda x: -1*np.std(x)+np.mean(x), lambda x: np.std(x)+np.mean(x)])
                    df_mean = df_group.groupby(x_axis)[out_par].mean()

                    param_plot = (
                        df_min_max.hvplot.area(
                            x=x_axis,
                            y="amin",
                            y2="amax",
                            alpha=0.5,
                            value_label=latexify.parse_word(out_par),
                            label="Spread")
                            # color="grey")
                        * df_mean.hvplot()
                        * df_std.hvplot.area(
                            x=x_axis,
                            y="<lambda_0>",
                            y2="<lambda_1>",
                            alpha=0.3,
                            value_label=latexify.parse_word(out_par),
                            label="sd")
                            # color="grey")
                    )
                elif c:
                    if out_par == x_axis:
                        df_group = df[[x_axis, "trajectory"]]
                    else:
                        df_group = df[[x_axis, out_par, "trajectory"]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                    param_plot = df_group.hvplot.scatter(
                        x=x_axis,
                        y=out_par,
                        c="trajectory",
                        title="Values of of {}".format(latexify.parse_word(out_par)),
                        label=None, datashade=True
                        )
                elif hexbin[0]:
                    if out_par == x_axis:
                        df_group = df[[x_axis, "trajectory"]]
                    else:
                        df_group = df[[x_axis, out_par, "trajectory"]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                    param_plot = df_group.hvplot.hexbin(
                        x=x_axis,
                        y=out_par,
                        title="Values of of {}".format(latexify.parse_word(out_par)),
                        label=None,
                        clabel="Count",
                        cmap="viridis",
                        logz=True,
                        gridsize=100,
                        colorbar=True
                        )
                else:
                    if out_par == x_axis:
                        df_group = df[[x_axis, "trajectory"]]
                    else:
                        df_group = df[[x_axis, out_par, "trajectory"]]
                    if log[0]:
                        df_group[out_par] = df_group[out_par].apply(lambda x: np.log(np.abs(x)))
                    param_plot = df_group.hvplot.scatter(
                        x=x_axis,
                        y=out_par,
                        title="Values of of {}".format(latexify.parse_word(out_par)),
                        label=None, datashade=True
                    )

                layout_kwargs = {}
                if self.backend == "bokeh":
                    layout_kwargs["width"] = 1600
                    layout_kwargs["height"] = 500

                if hist:
                    param_plot.opts(**layout_kwargs)
                    param_hist_plot = param_plot.hist(dimension=[x_axis, out_par], bins=bins, datashade=True)
                # Derivatives plot
                if hexbin[1]:
                    deriv_plot = df_tmp.hvplot.hexbin(
                        x=x_axis,
                        y="Derivative Ratio",
                        by="Derivatives",
                        title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
                        label=None, 
                        logz=True,
                        gridsize=100,
                        clabel="Count",
                        cmap="viridis",
                        colorbar=True,
                        width=1600,
                        height=500
                    )
                else:
                    deriv_plot = df_tmp.hvplot.scatter(
                        x=x_axis,
                        y="Derivative Ratio",
                        by="Derivatives",
                        title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
                        label=None, datashade=True
                    )

                if hist:
                    layout = param_hist_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)
                else:
                    layout = param_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)

                opts_arg = {}
                if self.backend == "matplotlib":
                    opts_arg["aspect"] = 3.2
                    opts_arg["fig_latex"] = False

                scatter_kwargs = opts_arg.copy()
                area_kwargs = opts_arg.copy()

                for k in kwargs:
                    scatter_kwargs[k] = kwargs[k]

                if self.backend == "bokeh":
                    scatter_kwargs["size"] = 5
                else:
                    scatter_kwargs["s"] = 8

                if self.backend == "matplotlib":
                    area_kwargs["edgecolor"] = None
                    area_kwargs["color"] = "black"
                    for k in kwargs:
                        area_kwargs[k] = kwargs[k]


                layout_kwargs = {}
                if self.backend == "matplotlib":
                    layout_kwargs["fig_size"] = 400

               # Matplotlib uses a horrible colormap as default...
                curve_kwargs = kwargs.copy()
                if percentile is not None:
                    if self.backend == "matplotlib":
                        if len(percentile) <= 10:
                            curve_kwargs["color"] = hv.Cycle("tab10")
                        elif len(percentile) <= 20:
                            curve_kwargs["color"] = hv.Cycle("tab20")
                        # More seems convoluted to me
#                         else:
#                             curve_kwargs["color"] = hv.Cycle("viridis")
                    else:
                        if len(percentile) <= 10:
                            curve_kwargs["color"] = hv.Cycle("Category10")
                        elif len(percentile) <= 20:
                            curve_kwargs["color"] = hv.Cycle("Category20")
#                         else:
#                             curve_kwargs["color"] = "viridis" # "colorcet"

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
                    if hexbin[0] and hexbin[1]:
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
                if hist:
                    param_opts = param_plot.opts(xaxis="bare", alpha=0.1)
                elif c:
                    param_opts = param_plot.opts(xaxis="bare", alpha=1.0)
                else:
                    param_opts = param_plot.opts(xaxis="bare")


                renderer = hv.Store.renderers[self.backend].instance(
                    fig=fig_type) # , dpi=300)
                latexify.set_size(True)

                i = 0
                if prefix is None:
                    prefix = "plt_1line_"
                    if scatter:
                        prefix = "plt_1scatter_"
                save = (plot_path + prefix + x_axis + "_" + out_par
                        + "_" + "{:03d}".format(i))
                while os.path.isfile(save + "." + fig_type):
                    i = i+1
                    save = (plot_path + prefix + x_axis + "_" + out_par
                            + "_" + "{:03d}".format(i))

                if self.backend == "bokeh":
                    print("Plotting")
                    hvplot.show(both_plots)
                    t2 = timer()
                    print("Plotting done in {}s".format(t2-t))
                else:
                    print("Saving to " + save + "." + fig_type)
                    renderer.save(both_plots, save)
                    t2 = timer()
                    try:
                        from IPython.display import Image, display
                        display(Image(save + "." + fig_type, width=1600))
                    except:
                        pass
                    print("Saving done in {}s".format(t2-t))

            i = 0
            if n_plots is None:
                n_plots = 9999999
            print("Creating {} plots".format(n_plots))
            while len(sorted_tuples) > 0 and i < n_plots:
                p, v = sorted_tuples.pop()
                in_params_2 = [p]
                while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                    and np.abs(v/sorted_tuples[-1][1]) < 10):
                    p, v = sorted_tuples.pop()
                    in_params_2.append(p)
                plot_helper(df_tmp_out, in_params=in_params_2, prefix=prefix, **kwargs)
                i += 1
                    
    def plot_inter(self, in_params, out_params, x_axis="timestep", mapped=True,
        trajectories=None, scatter=False, n_plots=None, percentile=None,
        frac=None, min_x=None, max_x=None, nth=None, hist=False,
        scatter_deriv=False, line_deriv=False, prefix=None, c=False,
        errorband=False, bins=50, plot_path="pics/", **kwargs):
        """
        Test case with interactive usage ie panels
        Plot two plots in two rows. At the top: Output parameter.
        At the bottom: Derivative with respect to that output parameter.
        Another plot is created for every gradient of different order.
        If multiple parameters are given, another plot is created for each
        output parameter. This function is mostly intended for multiple 
        trajectories.
        Does take long to plot each individual plot but may need less memory.

        Parameters
        ----------
        in_params : list of string
            Plot the derivatives with respect to those in this list.
        out_params : list of string
            List of keys to plot the derivatives for.
        x_axis : string
            The column to use as x-axis. Can be either "timestep" or
            an output parameter or a derivative.
        mapped : boolean
            If true: plot the region, where "MAP" is true, ie where the wcb
            criterion is fullfilled.
        trajectories : int
            The index of the trajectories to plot. If None is given, all
            trajectories will be plotted.
        scatter : boolean
            Plot a scatter plot or a line plot.
        n_plots : int
            Plot only that many plots. If None is given, plot every possible
            plot.
        percentile : list of int
            Plot the given percentiles along with the spread if
            errorband is not given. Not implemented yet
        errorband : Bool
            Plot an errorband (spread) using the standard deviation and min/max values.
        frac : float
            Sample a given fraction of rows. Deactivates "nth".
        nth : int
            Sample every nth entry. Works fast with "timestep" as x-axis and
            a given min_x and max_x value. If x_axis is any other than
            "timestep", an errorband triggered by "percentile" may not
            be plotted.
        scatter_deriv : boolean
            Plot the derivative plot with scatterplot. If scatter_deriv and
            line_deriv are set to false, use the same as provided via scatter
            and percentile.
        line_deriv : boolean
            Plot the derivative plot with lineplot. If scatter_deriv and
            line_deriv are set to false, use the same as provided via scatter
            and percentile.
        prefix : string
            Prefix to add to the filename.
        kwargs : dict
            Keyword arguments are passed down matplotlib.
        """
        import hvplot.dask # adds hvplot method to dask objects
        import hvplot.pandas
        import hvplot
        from holoviews import opts
        import holoviews as hv
        from timeit import default_timer as timer
        from holoviews.operation import histogram as hv_histo
#             from holoviews.operation.datashader import aggregate
        import pandas
        from holoviews.operation.datashader import datashade
#         import panel as pn
        from bokeh.layouts import column, row
        from bokeh.models import MultiSelect, RadioButtonGroup, Toggle, CheckboxGroup, Select
#         from bokeh.io import curdoc

#         pn.extension()
        hv.extension(self.backend)

        df = self.data

        def select_out_param(out_param="p"):
            return df.loc[df["Output Parameter"] == out_param]

        def select_traj(traj=0):
            return df.loc[df.trajectory == traj].compute()
        
        traj = df.trajectory.unique().compute().tolist()
#         traj_slider = pn.widgets.IntSlider(name="trajectory", start=int(traj[0]), end=int(traj[len(traj)-1]), value=int(traj[0]))
        
        
        
#         traj_select = pn.widgets.MultiSelect(name="Trajectory", value=traj, size=10, options=traj)
#         progress = pn.widgets.Progress(name="Progress", value=0, width=400)
#         in_selector = pn.widgets.CrossSelector(name="Input Parameters", value=[], options=in_params)
#         out_button = pn.widgets.RadioButtonGroup(name="Output Parameter", options=out_params, button_type="default")
#         mapped_toggle = pn.widgets.Toggle(name="Mapped", button_type="primary")
        
        
        
        
#         @pn.depends(out_button.param.value, traj_select.param.value)
        def title_fn(param, traj):
            return "## Output parameter {} - traj {}".format(param, traj[0])
        
#         @pn.depends(traj_select.param.value, out_button.param.value, in_selector.param.value, mapped_toggle.param.value)
        def plot_helper(trajs, out_param, in_params, mapped):
            # following https://holoviz.org/tutorial/Composing_Plots.html
            t = timer()
            df_tmp = df[in_params + [x_axis] + [out_param]]
            if mapped:
                df_tmp = df_tmp.loc[df_tmp.MAP == True]
            df_tmp = df_tmp.melt(x_axis, var_name="Derivatives",
                                value_name="Derivative Ratio")
            df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)

            if percentile is not None:

                if out_param == x_axis:
                    print("x-axis and y-axis are the same. Cannot plot that with percentiles!")
                    return
                else:
                    df_group = df[[x_axis, out_param]]

                # Group for min, max and percentiles
                funcs = [np.min, np.max] + [lambda x, perc=perc: np.percentile(x, perc, axis=0) for perc in percentile]
                df_min_max = df_group.groupby(x_axis).agg(funcs)[out_param]

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
                    datashade(df_min_max.hvplot.area(
                        x=x_axis,
                        y="Min",
                        y2="Max",
                        alpha=0.5,
                        value_label=latexify.parse_word(out_param),
                        title="Values of of {}".format(latexify.parse_word(out_param)),
                        label="Spread",
                        color="grey"))
                    * datashade(df_min_max.hvplot.line(
                        x=x_axis,
                        y=p_list,
                        value_label=latexify.parse_word(out_param),
                        **kwargs))
                )
            elif errorband:
                df_group = df[[x_axis, out_param, "trajectory"]]
                df_min_max = df_group.groupby([x_axis, "trajectory"])[out_param].mean().groupby(x_axis).agg([np.min, np.max])
                df_std = df_group.groupby([x_axis, "trajectory"])[out_param].mean().groupby(x_axis).agg([lambda x: -1*np.std(x)+np.mean(x), lambda x: np.std(x)+np.mean(x)])
                df_mean = df_group.groupby(x_axis)[out_param].mean()

                param_plot = (
                    df_min_max.hvplot.area(
                        x=x_axis,
                        y="amin",
                        y2="amax",
                        alpha=0.5,
                        value_label=latexify.parse_word(out_param),
                        label="Spread")
                        # color="grey")
                    * df_mean.hvplot()
                    * df_std.hvplot.area(
                        x=x_axis,
                        y="<lambda_0>",
                        y2="<lambda_1>",
                        alpha=0.3,
                        value_label=latexify.parse_word(out_param),
                        label="sd")
                        # color="grey")
                )
            elif c:
                if out_param == x_axis:
                    df_group = df[[x_axis, "trajectory"]]
                else:
                    df_group = df[[x_axis, out_param, "trajectory"]]
                param_plot = df_group.hvplot.scatter(
                    x=x_axis,
                    y=out_param,
                    c="trajectory",
                    title="Values of of {}".format(latexify.parse_word(out_param)),
                    label=None, datashade=True
                    )
            else:
                if out_param == x_axis:
                    df_group = df[[x_axis, "trajectory"]]
                else:
                    df_group = df[[x_axis, out_param, "trajectory"]]
                param_plot = df_group.hvplot.scatter(
                    x=x_axis,
                    y=out_param,
                    title="Values of of {}".format(latexify.parse_word(out_param)),
                    label=None, datashade=True
                )

            layout_kwargs = {}
            if self.backend == "bokeh":
                layout_kwargs["width"] = 1600
                layout_kwargs["height"] = 500

            if hist:
                param_plot.opts(**layout_kwargs)
                param_hist_plot = param_plot.hist(dimension=[x_axis, out_param], bins=bins, datashade=True)

            deriv_plot = df_tmp.hvplot.scatter(
                x=x_axis,
                y="Derivative Ratio",
                by="Derivatives",
                title="Deriv. Ratio of {}".format(latexify.parse_word(out_param)),
                label=None, datashade=True
            )

            if hist:
                layout = param_hist_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)
            else:
                layout = param_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)

            opts_arg = {}
            if self.backend == "matplotlib":
                opts_arg["aspect"] = 3.2
                opts_arg["fig_latex"] = False

            scatter_kwargs = opts_arg.copy()
            area_kwargs = opts_arg.copy()

            for k in kwargs:
                scatter_kwargs[k] = kwargs[k]

            if self.backend == "bokeh":
                scatter_kwargs["size"] = 5
            else:
                scatter_kwargs["s"] = 8

            if self.backend == "matplotlib":
                area_kwargs["edgecolor"] = None
                area_kwargs["color"] = "black"
                for k in kwargs:
                    area_kwargs[k] = kwargs[k]


            layout_kwargs = {}
            if self.backend == "matplotlib":
                layout_kwargs["fig_size"] = 400

           # Matplotlib uses a horrible colormap as default...
            curve_kwargs = kwargs.copy()
            if percentile is not None:
                if self.backend == "matplotlib":
                    if len(percentile) <= 10:
                        curve_kwargs["color"] = hv.Cycle("tab10")
                    elif len(percentile) <= 20:
                        curve_kwargs["color"] = hv.Cycle("tab20")
                    # More seems convoluted to me
#                         else:
#                             curve_kwargs["color"] = hv.Cycle("viridis")
                else:
                    if len(percentile) <= 10:
                        curve_kwargs["color"] = hv.Cycle("Category10")
                    elif len(percentile) <= 20:
                        curve_kwargs["color"] = hv.Cycle("Category20")
#                         else:
#                             curve_kwargs["color"] = "viridis" # "colorcet"

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

            if self.backend == "matplotlib":
                both_plots = both_plots.opts(sublabel_format="", tight=True)
            if hist:
                param_opts = param_plot.opts(xaxis="bare", alpha=0.1)
            elif c:
                param_opts = param_plot.opts(xaxis="bare", alpha=1.0)
            else:
                param_opts = param_plot.opts(xaxis="bare")
            return both_plots
        
        def update_plot(attrname, old, new):
            trajs = traj_select.value
            out_param = out_button.value
            in_params = in_selector.value
            mapped = mapped_toggle.value
            plot_helper(trajs, out_param, in_params, mapped)
            
        plot = plot_helper(trajectories, out_params[0], in_params, mapped)
        
#         traj_select = MultiSelect(title="Trajectory", value=traj[0], options=traj)
#         in_selector = CheckboxGroup(labels=in_params, active=[0])
#         out_button = RadioButtonGroup(labels=out_params, active=0)
#         mapped_toggle = Toggle(label="Mapped", button_type="default")
        
        traj_select = Select(value=traj[0], title="Trajectory", options=traj)
        in_selector = Select(value=in_params[0], title="Input Parameters", options=in_params)
        out_button = Select(value=out_params[0], title="Output Parameter", options=out_params)
        mapped_toggle = Select(value=False, title="Mapped", options=[True, False])
            
        traj_select.on_change('value', update_plot)
        in_selector.on_change('value', update_plot)
        out_button.on_change('value', update_plot)
        mapped_toggle.on_change('value', update_plot)
        
        selectors = row(
            column(mapped_toggle, out_button),
            column(traj_select, in_selector))
        
        self.widget = row(selectors, plot)
        curdoc().add_root(self.widget)
        curdoc().title("Weather")
        
#         header = pn.Row(pn.panel(title_fn, width=400))
#         selectors = pn.Column(
#             pn.Row(
#                 pn.Column(mapped_toggle), pn.Column(out_button)),
#             pn.Row(
#                 pn.Column(traj_select), pn.Column(in_selector)),
#             pn.Row(
#                 pn.Column(progress))
#         )
#         self.widget = pn.Column(header, selectors)
        
#         plots = pn.Row(hv.DynamicMap(plot_helper))
#         self.widget = pn.Column(header, selectors, plots)

        
#         pn_traj = pn.interact(select_traj, traj=(0, df.trajectory.unique()))
#         pn_param = pn.interact(select_out_param, out_param=df["Output Parameter"].unique())
#         title = pn.panel("# Output simulation", width=400)
#         header = pn.Row(title)
#         body = pn.Row(
#             pn.Column("### Choose a Trajectory", traj_slider),
#             pn.Column("### Output parameter", pn_param),
# #             pn.Column("### Some plot", plot)
#         )
#         return pn.Column(header, body)
        
#         return pn_traj

#         if frac is not None:
#             df = self.data.sample(frac=frac, replace=False, random_state=42)
#         elif nth is not None:
#             if min_x is not None and max_x is not None and x_axis == "timestep":
#                 steps = np.arange(min_x, max_x, nth*2.5)
#             elif x_axis == "timestep":
#                 df_tmp = self.data.timestep.unique().compute()
#                 min_val = df_tmp.min()
#                 max_val = df_tmp.max()
#                 steps = np.arange(min_val, max_val, nth*2.5)
#             else:
#                 steps = self.data[x_axis].unique().compute().to_numpy()[::nth]

#             df = self.data.loc[self.data[x_axis].isin(steps)]
#         else:
#             df = self.data
#         if min_x is not None:
#             df = df.loc[df[x_axis] >= min_x]
#         if max_x is not None:
#             df = df.loc[df[x_axis] <= max_x]
#         if trajectories is not None:
#             df = df.loc[df.trajectory.isin(trajectories)]
#         if mapped:
#             df = df.loc[df.MAP == True]



#         df = df.loc[df["Output Parameter"].isin(out_params)]



#         for out_par in out_params:
#             df_tmp_out = df.loc[df["Output Parameter"] == out_par]

#             # Sort the derivatives
#             sorted_tuples = []
#             for in_p in in_params:
#                 value = np.abs(df_tmp_out[in_p].min())
#                 max_val = np.abs(df_tmp_out[in_p].max())
#                 if max_val > value:
#                     value = max_val
#                 if value != 0 and not np.isnan(value):
#                     sorted_tuples.append((in_p, value))
#             sorted_tuples.sort(key=lambda tup: tup[1])

#             def plot_helper(df, in_params, prefix, **kwargs):
#                 # following https://holoviz.org/tutorial/Composing_Plots.html
#                 t = timer()
#                 df_tmp = df[in_params+[x_axis]]
#                 df_tmp = df_tmp.melt(x_axis, var_name="Derivatives",
#                                     value_name="Derivative Ratio")
#                 df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)

#                 if percentile is not None:

#                     if out_par == x_axis:
#                         print("x-axis and y-axis are the same. Cannot plot that with percentiles!")
#                         return
#                     else:
#                         df_group = df[[x_axis, out_par]]

#                     # Group for min, max and percentiles
#                     funcs = [np.min, np.max] + [lambda x, perc=perc: np.percentile(x, perc, axis=0) for perc in percentile]
#                     df_min_max = df_group.groupby(x_axis).agg(funcs)[out_par]

#                     # Rename the columns
#                     p_list = []
#                     p_dic = {}
#                     for i, perc in enumerate(percentile):
#                         p_list.append("{}. Percentile".format(perc))
#                         p_dic["<lambda_{}>".format(i)] = p_list[-1]
#                     p_dic["amin"] = "Min"
#                     p_dic["amax"] = "Max"
#                     df_min_max = df_min_max.rename(columns=p_dic)

#                     # Plot
#                     param_plot = (
#                         datashade(df_min_max.hvplot.area(
#                             x=x_axis,
#                             y="Min",
#                             y2="Max",
#                             alpha=0.5,
#                             value_label=latexify.parse_word(out_par),
#                             title="Values of of {}".format(latexify.parse_word(out_par)),
#                             label="Spread",
#                             color="grey"))
#                         * datashade(df_min_max.hvplot.line(
#                             x=x_axis,
#                             y=p_list,
#                             value_label=latexify.parse_word(out_par),
#                             **kwargs))
#                     )
#                 elif errorband:
#                     df_group = df[[x_axis, out_par, "trajectory"]]
#                     df_min_max = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([np.min, np.max])
#                     df_std = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([lambda x: -1*np.std(x)+np.mean(x), lambda x: np.std(x)+np.mean(x)])
#                     df_mean = df_group.groupby(x_axis)[out_par].mean()

#                     param_plot = (
#                         df_min_max.hvplot.area(
#                             x=x_axis,
#                             y="amin",
#                             y2="amax",
#                             alpha=0.5,
#                             value_label=latexify.parse_word(out_par),
#                             label="Spread")
#                             # color="grey")
#                         * df_mean.hvplot()
#                         * df_std.hvplot.area(
#                             x=x_axis,
#                             y="<lambda_0>",
#                             y2="<lambda_1>",
#                             alpha=0.3,
#                             value_label=latexify.parse_word(out_par),
#                             label="sd")
#                             # color="grey")
#                     )
#                 elif c:
#                     if out_par == x_axis:
#                         df_group = df[[x_axis, "trajectory"]]
#                     else:
#                         df_group = df[[x_axis, out_par, "trajectory"]]
#                     param_plot = df_group.hvplot.scatter(
#                         x=x_axis,
#                         y=out_par,
#                         c="trajectory",
#                         title="Values of of {}".format(latexify.parse_word(out_par)),
#                         label=None, datashade=True
#                         )
#                 else:
#                     if out_par == x_axis:
#                         df_group = df[[x_axis, "trajectory"]]
#                     else:
#                         df_group = df[[x_axis, out_par, "trajectory"]]
#                     param_plot = df_group.hvplot.scatter(
#                         x=x_axis,
#                         y=out_par,
#                         title="Values of of {}".format(latexify.parse_word(out_par)),
#                         label=None, datashade=True
#                     )

#                 layout_kwargs = {}
#                 if self.backend == "bokeh":
#                     layout_kwargs["width"] = 1600
#                     layout_kwargs["height"] = 500

#                 if hist:
#                     param_plot.opts(**layout_kwargs)
#                     param_hist_plot = param_plot.hist(dimension=[x_axis, out_par], bins=bins, datashade=True)

#                 deriv_plot = df_tmp.hvplot.scatter(
#                     x=x_axis,
#                     y="Derivative Ratio",
#                     by="Derivatives",
#                     title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
#                     label=None, datashade=True
#                 )

#                 if hist:
#                     layout = param_hist_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)
#                 else:
#                     layout = param_plot.opts(**layout_kwargs) + deriv_plot.opts(**layout_kwargs)

#                 opts_arg = {}
#                 if self.backend == "matplotlib":
#                     opts_arg["aspect"] = 3.2
#                     opts_arg["fig_latex"] = False

#                 scatter_kwargs = opts_arg.copy()
#                 area_kwargs = opts_arg.copy()

#                 for k in kwargs:
#                     scatter_kwargs[k] = kwargs[k]

#                 if self.backend == "bokeh":
#                     scatter_kwargs["size"] = 5
#                 else:
#                     scatter_kwargs["s"] = 8

#                 if self.backend == "matplotlib":
#                     area_kwargs["edgecolor"] = None
#                     area_kwargs["color"] = "black"
#                     for k in kwargs:
#                         area_kwargs[k] = kwargs[k]


#                 layout_kwargs = {}
#                 if self.backend == "matplotlib":
#                     layout_kwargs["fig_size"] = 400

#                # Matplotlib uses a horrible colormap as default...
#                 curve_kwargs = kwargs.copy()
#                 if percentile is not None:
#                     if self.backend == "matplotlib":
#                         if len(percentile) <= 10:
#                             curve_kwargs["color"] = hv.Cycle("tab10")
#                         elif len(percentile) <= 20:
#                             curve_kwargs["color"] = hv.Cycle("tab20")
#                         # More seems convoluted to me
# #                         else:
# #                             curve_kwargs["color"] = hv.Cycle("viridis")
#                     else:
#                         if len(percentile) <= 10:
#                             curve_kwargs["color"] = hv.Cycle("Category10")
#                         elif len(percentile) <= 20:
#                             curve_kwargs["color"] = hv.Cycle("Category20")
# #                         else:
# #                             curve_kwargs["color"] = "viridis" # "colorcet"

#                 if errorband:
#                     both_plots = layout.opts(
#                         opts.Area(
#                             xticks=20,
#                             xaxis="bottom",
#                             fontsize=self.font_dic,
#                             show_grid=True,
#                             show_legend=True,
#                             **area_kwargs),
#                         opts.Scatter(
#                             xticks=20,
#                             xaxis="bottom",
#                             fontsize=self.font_dic,
#                             show_grid=True,
#                             show_legend=True,
#                             **scatter_kwargs),
#                         opts.Layout(**layout_kwargs)
#                     ).cols(1)
#                 elif percentile is not None:
#                     both_plots = layout.opts(
#                         opts.Area(
#                             xticks=20,
#                             xaxis="bottom",
#                             fontsize=self.font_dic,
#                             show_grid=True,
#                             show_legend=True,
#                             **area_kwargs),
#                         opts.Scatter(
#                             xticks=20,
#                             xaxis="bottom",
#                             fontsize=self.font_dic,
#                             show_grid=True,
#                             show_legend=True,
#                             **scatter_kwargs),
#                         opts.Curve(**curve_kwargs),
#                         opts.Layout(**layout_kwargs)
#                     ).cols(1)
#                 else:
#                     both_plots = layout.opts(
#                         opts.Scatter(
#                             xticks=20,
#                             xaxis="bottom",
#                             fontsize=self.font_dic,
#                             show_grid=True,
#                             show_legend=True,
#                             **scatter_kwargs),
#                         opts.Layout(**layout_kwargs)
#                     ).cols(1)

#                 if self.backend == "matplotlib":
#                     both_plots = both_plots.opts(sublabel_format="", tight=True)
#                 if hist:
#                     param_opts = param_plot.opts(xaxis="bare", alpha=0.1)
#                 elif c:
#                     param_opts = param_plot.opts(xaxis="bare", alpha=1.0)
#                 else:
#                     param_opts = param_plot.opts(xaxis="bare")


#                 renderer = hv.Store.renderers[self.backend].instance(
#                     fig='svg') # , dpi=300)
#                 latexify.set_size(True)

#                 i = 0
#                 if prefix is None:
#                     prefix = "plt_1line_"
#                     if scatter:
#                         prefix = "plt_1scatter_"
#                 save = (plot_path + prefix + x_axis + "_" + out_par
#                         + "_" + "{:03d}".format(i))
#                 while os.path.isfile(save + ".svg"):
#                     i = i+1
#                     save = (plot_path + prefix + x_axis + "_" + out_par
#                             + "_" + "{:03d}".format(i))

#                 if self.backend == "bokeh":
#                     print("Plotting")
#                     hvplot.show(both_plots)
#                     t2 = timer()
#                     print("Plotting done in {}s".format(t2-t))
#                 else:
#                     print("Saving to " + save + ".svg")
#                     renderer.save(both_plots, save)
#                     t2 = timer()
#                     from IPython.display import Image, display
#                     display(Image(save + ".svg", width=1600))
#                     print("Saving done in {}s".format(t2-t))

#             i = 0
#             if n_plots is None:
#                 n_plots = 9999999
#             print("Creating {} plots".format(n_plots))
#             while len(sorted_tuples) > 0 and i < n_plots:
#                 p, v = sorted_tuples.pop()
#                 in_params_2 = [p]
#                 while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
#                     and np.abs(v/sorted_tuples[-1][1]) < 10):
#                     p, v = sorted_tuples.pop()
#                     in_params_2.append(p)
#                 plot_helper(df_tmp_out, in_params=in_params_2, prefix=prefix, **kwargs)
#                 i += 1

## Idea: https://www.holoviews.org/user_guide/Linking_Plots.html
