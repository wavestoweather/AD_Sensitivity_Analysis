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
# matplotlib.use('Agg')
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

    def plot_two(self, in_params, out_params, x_axis="timestep", mapped=True,
        trajectories=None, scatter=False, n_plots=None, percentile=None,
        frac=None, min_x=None, max_x=None, nth=None,
        scatter_deriv=False, line_deriv=False, **kwargs):
        """
        Plot two plots in two rows. At the top: Output parameter.
        At the bottom: Derivative with respect to that output parameter.
        If multiple parameters are given, another plot is created for each
        output parameter and another plot is created for every gradient of
        different order.

        Parameters
        ----------
        in_params : list of string
            Plot the derivatives with respect to those in this list.
        out_param : list of string
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
        percentile : int
            If None: Plot different trajectories as distinct lines.
            If value : Plot a mean trajectory and given percentile
            (kexword "ci" from seaborn).
            If percentile is given, dots and scatter is deactivated.
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
        kwargs : dict
            Keyword arguments are passed down to seaborn.scatterplott() or
            seaborn.lineplot() for the derivative plots.
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

        df = df.loc[df["Output Parameter"].isin(out_params)].compute()
        # Memory usage tool
        print("Calculate memory usage", flush=True)
        memory_usage = df.memory_usage(deep=True)
        # print("Memory usage in MB: {}".format(memory_usage/1e6), flush=True)
        print("Total memory needed: {} GB".format(memory_usage.sum()/1e9), flush=True)

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

            def plot_helper(df, in_params, **kwargs):
                fig, axs = plt.subplots(2, 1, sharex=True, sharey=False)
                min_time = df[x_axis].min()
                max_time = df[x_axis].max()
                dt1 = (max_time - min_time) / 20
                dt = (max_time - min_time + dt1) / 20
                x_ticks = np.arange(min_time, max_time + dt1, dt)

                if len(df.trajectory.unique()) > 1:
                    if percentile is not None:
                        sns.lineplot(x=x_axis, y=out_par, ci=percentile,
                                        data=df, ax=axs[0], err_style="band",
                                        estimator="mean",
                                        **kwargs)
                    elif scatter:
                        sns.scatterplot(x=x_axis, y=out_par,
                                        data=df, hue="trajectory", ax=axs[0],
                                        palette=sns.color_palette("husl", len(df.trajectory.unique())),
                                        **kwargs)
                    else:
                        # n_trajs = len(df.trajectory.unique())
                        # color = sns.color_palette("husl", n_trajs).as_hex()
                        # for traj in range(n_trajs):
                        #     axs[0].plot(x_axis,
                        #             out_par,
                        #             data=df.loc[df.trajectory == traj],
                        #             # fmt=color[traj],
                        #             **kwargs)
                        # axs[0].legend(df.trajectory.unique(), title="trajectory")
                        sns.lineplot(x=x_axis, y=out_par, ci=None,
                                        data=df, hue="trajectory", ax=axs[0],
                                        palette=sns.color_palette("husl", len(df.trajectory.unique())),
                                        **kwargs)
                else:
                    if scatter:
                        sns.scatterplot(x=x_axis, y=out_par,
                                        data=df, ax=axs[0], **kwargs)
                    else:

                        sns.lineplot(x=x_axis, y=out_par, ci=None,
                                        data=df, ax=axs[0], **kwargs)
                    axs[0].legend(df.trajectory.unique(), title="trajectory")

                df_tmp = df[in_params+[x_axis]]
                df_tmp = df_tmp.melt(x_axis, var_name="Derivatives",
                                    value_name="Derivative Ratio")

                if not scatter_deriv and not line_deriv:
                    if len(df.trajectory.unique()) > 1 and percentile is not None:
                        sns.lineplot(x=x_axis, y="Derivative Ratio", ci=percentile,
                                        data=df_tmp, hue="Derivatives",
                                        ax=axs[1], err_style="band",
                                        estimator="mean",
                                        **kwargs)
                    elif scatter:
                        sns.scatterplot(x=x_axis, y="Derivative Ratio",
                                            data=df_tmp, hue="Derivatives",
                                            ax=axs[1], linewidth=0,
                                            **kwargs)
                    else:
                        # derivs = df.Derivatives.unique()
                        # n_derivs = len(derivs)
                        # color = sns.color_palette("husl", n_derivs).as_hex()
                        # for i, deriv in enumerate(derivs):
                        #     axs[1].plot(x_axis,
                        #             "Derivative Ratio",
                        #             data=df.loc[df.Derivatives == deriv],
                        #             # fmt=color[i],
                        #             **kwargs)
                        # axs[1].legend(derivs, title="Derivatives")
                        sns.lineplot(x=x_axis, y="Derivative Ratio",
                                        data=df_tmp, hue="Derivatives",
                                        ax=axs[1],
                                        **kwargs)
                elif scatter_deriv:
                    sns.scatterplot(x=x_axis, y="Derivative Ratio",
                                            data=df_tmp, hue="Derivatives",
                                            ax=axs[1], linewidth=0,
                                            **kwargs)
                else:
                    sns.lineplot(x=x_axis, y="Derivative Ratio",
                                        data=df_tmp, hue="Derivatives",
                                        ax=axs[1],
                                        **kwargs)

                axs[1].set_title("Deriv. Ratio of {}".format(latexify.parse_word(out_par)))
                axs[1].set_xticks(x_ticks)

                # Change labels to latex versions
                legend = axs[1].get_legend()
                _, labels = axs[1].get_legend_handles_labels()
                for t, old in zip(legend.texts, labels):
                    t.set_text(latexify.parse_word(old))
                axs[1].set_ylabel("Derivative ratio")
                plt.ticklabel_format(style="scientific", axis="y", scilimits=(0,0))

                # Plot the area that had been flagged
                if mapped:
                    df_mapped = df[df.MAP == True]
                    if not df_mapped.empty:
                        xmin = df_mapped[x_axis].min()
                        xmax = df_mapped[x_axis].max()
                        axs[0].axvspan(xmin=xmin, xmax=xmax,
                            facecolor="khaki", alpha=0.3)
                        axs[1].axvspan(xmin=xmin, xmax=xmax,
                            facecolor="khaki", alpha=0.3)

                # Change the limits for the y-axis because sometimes that
                # can be off and it is hard to see anything.
                min_y = df_tmp["Derivative Ratio"].min()
                max_y = df_tmp["Derivative Ratio"].max()
                axs[1].set_ylim(min_y, max_y)

                min_y = df[out_par].min()
                max_y = df[out_par].max()
                axs[0].set_ylim(min_y, max_y)
                # The same is true for the x-axis
                plt.xlim(x_ticks[0], x_ticks[-1])

                i = 0
                prefix = "plt_2line_"
                if scatter:
                    prefix = "plt_2scatter_"
                save = ("pics/" + prefix + x_axis + "_" + out_par
                        + "_" + "{:03d}".format(i) + ".png")
                while os.path.isfile(save):
                    i = i+1
                    save = ("pics/" + prefix + x_axis + "_" + out_par
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

            i = 0
            if n_plots is None:
                n_plots = 9999999
            while len(sorted_tuples) > 0 and i < n_plots:
                p, v = sorted_tuples.pop()
                in_params_2 = [p]
                while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                    and np.abs(v/sorted_tuples[-1][1]) < 10):
                    p, v = sorted_tuples.pop()
                    in_params_2.append(p)
                plot_helper(df_tmp_out, in_params=in_params_2, **kwargs)
                i += 1

    def plot_two_ds(self, in_params, out_params, x_axis="timestep", mapped=True,
            trajectories=None, scatter=False, n_plots=None, percentile=None,
            frac=None, min_x=None, max_x=None, nth=None, hist=False,
            scatter_deriv=False, line_deriv=False, prefix=None, c=False,
            errorband=False, bins=50, plot_path="pics/", **kwargs):
            """
            Plot two plots in two rows. At the top: Output parameter.
            At the bottom: Derivative with respect to that output parameter.
            If multiple parameters are given, another plot is created for each
            output parameter and another plot is created for every gradient of
            different order.

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

            df = df.loc[df["Output Parameter"].isin(out_params)].compute()
            # Memory usage tool
            print("Calculate memory usage", flush=True)
            memory_usage = df.memory_usage(deep=True)
            # print("Memory usage in MB: {}".format(memory_usage/1e6), flush=True)
            print("Total memory needed: {} GB".format(memory_usage.sum()/1e9), flush=True)

            import hvplot.dask # adds hvplot method to dask objects
            import hvplot.pandas
            import hvplot
            from holoviews import opts
            import holoviews as hv
            from timeit import default_timer as timer
            from holoviews.operation import histogram as hv_histo
            import pandas
            # Hotfix
            # matplotlib.rcParams['agg.path.chunksize'] = 9e16

            # matplotlib.rcParams['agg.path.chunksize'] = 0
            # matplotlib.rcParams.update( matplotlib.rc_params() )

            hv.extension(self.backend)
            # hv.extension('bokeh')

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
                    df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)

                    if percentile is not None:
                        # print("Select values.")
                        if out_par == x_axis:
                            df_group = df[[x_axis, "trajectory"]]
                        else:
                            df_group = df[[x_axis, out_par, "trajectory"]]
                        x_values = np.sort(df_group[x_axis].unique())
                        # print("Selecting done.")
                        percentiles = []
                        min_vals = []
                        max_vals = []
                        # print("Getting percentiles")
                        for val in x_values:
                            df_tmp2 = df_group.loc[df_group[x_axis] == val]
                            min_vals.append(df_tmp2[out_par].min())
                            max_vals.append(df_tmp2[out_par].max())

                            # percentiles.append( df_tmp2.array.percentile() )

                            percentiles.append(np.percentile(
                                df_tmp2[out_par], percentile, axis=0))
                        # print("Getting percentiles done")
                        pandas_df = {
                            x_axis: x_values,
                            "Min": min_vals,
                            "Max": max_vals}
                        percentiles = np.transpose(percentiles)
                        p_list = []
                        for i, perc in enumerate(percentile):
                            pandas_df["{}. Percentile".format(perc)] = percentiles[i]
                            p_list.append("{}. Percentile".format(perc))

                        # print("Percentiles shape: {}".format(np.shape(percentiles)))
                        # print("Plot percentiles.")
                        # for key in pandas_df:
                        #     print("{}: {}".format(key, len(pandas_df[key])))
                        pandas_df = pandas.DataFrame(pandas_df)
                        # perc_plot = pandas_df.hvplot(x=x_axis, y=p_list,
                        #     value_label=latexify.parse_word(out_par))
                        # print("Plot area.")
                        # param_plot = perc_plot.hvplot.area(
                        #     x=x_axis, y="Min", y2="Max", alpha=0.2)
                        # print("Plots created")
                        print(pandas_df.describe())

                        df_min_max = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([np.min, np.max])
                        # print(df_min_max.describe())
                        # print(p_list)
                        # print(pandas_df.head(100))
                        param_plot = (
                            df_min_max.hvplot.area(
                                x=x_axis,
                                y="amin",
                                y2="amax",
                                alpha=0.5,
                                value_label=latexify.parse_word(out_par),
                                title="Values of of {}".format(latexify.parse_word(out_par)),
                                label="Spread",
                                color="grey")
                            * pandas_df.hvplot.line(
                                x=x_axis,
                                y=p_list,
                                value_label=latexify.parse_word(out_par),
                                **kwargs)
                        )

                        # df_tmp2 = df_group.melt()
                        # percentiles = df.array.percentile(df_group[out_par], percentile, axis=)
                        # for perc in percentile:
                        #     df_group["{}.Percentile".format(perc)] = (
                        #         df.array.percentile(df_group[out_par], perc))


                        # param_plot = df_group.hvplot.scatter(
                        #     x=x_axis,
                        #     y=out_par,
                        #     title="Values of of {}".format(latexify.parse_word(out_par)),
                        #     label=None
                        #     )
                    elif errorband:
                        # print("Current limit {}".format(matplotlib.rcParams['agg.path.chunksize']))
                        df_group = df[[x_axis, out_par, "trajectory"]]
                        df_min_max = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([np.min, np.max])
                        df_std = df_group.groupby([x_axis, "trajectory"])[out_par].mean().groupby(x_axis).agg([lambda x: -1*np.std(x)+np.mean(x), lambda x: np.std(x)+np.mean(x)])
                        df_mean = df_group.groupby(x_axis)[out_par].mean()
                        # print(df_min_max.describe())
                        # print(df_std.describe())
                        # print(df_mean.describe())
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
                        # Alternative with Spread? (maybe faster)
                        # param_plot = (
                        #     df_min_max.hvplot.area(
                        #         x=x_axis,
                        #         y="amin",
                        #         y2="amax",
                        #         alpha=0.3,
                        #         value_label=latexify.parse_word(out_par),
                        #         label="Spread")
                        #         # color="blue")
                        #     * df_mean.hvplot()
                        #     * df_mean.hvplot.Spread(
                        #         x=x_axis,
                        #         y=out_par,
                        #         err=,
                        #         alpha=0.2,
                        #         value_label=latexify.parse_word(out_par),
                        #         label="sd")
                        #         # color="blue")
                        # )
                    elif c:
                        if out_par == x_axis:
                            df_group = df[[x_axis, "trajectory"]]
                        else:
                            df_group = df[[x_axis, out_par, "trajectory"]]
                        param_plot = df_group.hvplot.scatter(
                            x=x_axis,
                            y=out_par,
                            c="trajectory",
                            title="Values of of {}".format(latexify.parse_word(out_par)),
                            label=None
                            )
                    else:
                        if out_par == x_axis:
                            df_group = df[[x_axis, "trajectory"]]
                        else:
                            df_group = df[[x_axis, out_par, "trajectory"]]
                        param_plot = df_group.hvplot.scatter(
                        x=x_axis,
                        y=out_par,
                        title="Values of of {}".format(latexify.parse_word(out_par)),
                        label=None
                        )

                    if hist:
                        param_hist_plot = param_plot.hist(dimension=[x_axis, out_par], bins=bins)

                    deriv_plot = df_tmp.hvplot.scatter(
                        x=x_axis,
                        y="Derivative Ratio",
                        by="Derivatives",
                        title="Deriv. Ratio of {}".format(latexify.parse_word(out_par)),
                        label=None
                        )
                    # print(deriv_plot)
                    if hist:
                    #     # deriv_hist_plot = deriv_plot.hist(dimension=["Derivative Ratio"], bins=bins)
                        layout = param_hist_plot + deriv_plot
                    else:
                        layout = param_plot + deriv_plot

                    if errorband:
                        both_plots = layout.opts(
                            opts.Area(aspect=3.2,
                                edgecolor=None, # or linewidth=0
                                xticks=20,
                                xaxis="bottom",
                                fontsize=self.font_dic,
                                show_grid=True,
                                show_legend=True,
                                fig_latex=False),
                            opts.Scatter(aspect=3.2,
                                xticks=20,
                                xaxis="bottom",
                                s=5,
                                fontsize=self.font_dic,
                                show_grid=True,
                                show_legend=True,
                                fig_latex=False,
                                **kwargs),
                            opts.Layout(fig_size=400)
                        ).cols(1)
                    elif percentile is not None:
                        both_plots = layout.opts(
                            opts.Area(aspect=3.2,
                                edgecolor=None, # or linewidth=0
                                xticks=20,
                                xaxis="bottom",
                                fontsize=self.font_dic,
                                show_grid=True,
                                show_legend=True,
                                fig_latex=False),
                            opts.Scatter(aspect=3.2,
                                xticks=20,
                                xaxis="bottom",
                                s=5,
                                fontsize=self.font_dic,
                                show_grid=True,
                                show_legend=True,
                                fig_latex=False,
                                **kwargs),
                            opts.Layout(fig_size=400)
                        ).cols(1)
                    else:
                        both_plots = layout.opts(
                            opts.Scatter(aspect=3.2,
                                xticks=20,
                                xaxis="bottom",
                                s=5,
                                fontsize=self.font_dic,
                                show_grid=True,
                                show_legend=True,
                                fig_latex=False,
                                **kwargs),
                            opts.Layout(fig_size=400)
                        ).cols(1)

                    both_plots = both_plots.opts(sublabel_format="", tight=True)
                    if hist:
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


                    if self.backend == "bokeh":
                        print("Plotting")
                        hvplot.show(both_plots)
                        t2 = timer()
                        print("Plotting done in {}s".format(t2-t))
                    else:
                        print("Saving to " + save + ".png")
                        renderer.save(both_plots, save)
                        t2 = timer()
                        from IPython.display import Image, display
                        display(Image(save + ".png", width=1600))
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


    def plot_andrew(self, group, out_params, params, plot_path="pics/",
                    prefix=None, sort=False, trajectories=None, **kwargs):
            """
            Plot an Andrew plot where group is used to colorize every
            datapoint, params consists of the column names to be used.
            The dimensions are sorted by their magnitude. If sort is True,
            only dimensions of the same magnitude are used together.

            Parameters
            ----------
            group : string
                Group the datapoints according to this group, i.e. trajectory or p
            out_params : list of string
                List of output parameters for whose derivatives are used.
            params : list of string
                List of columns used for the observations.
            plot_path : string
                # Path to folder for saving the pictures.
            sort : boolean
                Sort the derivatives and plot Andrew curves only for those
                which are in the same order of magnitude.
            trajectories : int
                The index of the trajectories to plot. If None is given, all
                trajectories will be plotted.
            kwargs : dict
                Keyword arguments are passed down matplotlib.
            """
            all_cols = params + [group]
            if "trajectory" not in all_cols:
                all_cols.append("trajectory")
            if "Output Parameter" not in all_cols:
                all_cols.append("Output Parameter")
            df = self.data[all_cols]
            if trajectories is not None:
                df = df.loc[df.trajectory.isin(trajectories)]

            df = df.loc[df["Output Parameter"].isin(out_params)].compute()
            # Memory usage tool
            print("Calculate memory usage", flush=True)
            memory_usage = df.memory_usage(deep=True)
            # print("Memory usage in MB: {}".format(memory_usage/1e6), flush=True)
            print("Total memory needed: {} GB".format(memory_usage.sum()/1e9), flush=True)

            import hvplot.dask # adds hvplot method to dask objects
            import hvplot.pandas
            from holoviews import opts
            import holoviews as hv
            from timeit import default_timer as timer
            hv.extension(self.backend)

            for out_par in out_params:
                df_tmp_out = df.loc[df["Output Parameter"] == out_par]

                if sort:
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

                def plot_helper(df, params, prefix, **kwargs):
                    t = timer()
                    df_tmp = df[params+[group]]
                    alpha = 0.5
                    # if "alpha" in kwargs:
                    #     alpha = kwargs["alpha"]
                    samples = 200
                    if "samples" in kwargs:
                        samples = kwargs["samples"]
                    colormap = None
                    if "colormap" in kwargs:
                        colormap = kwargs["colormap"]

                    print(df_tmp.head())
                    print(df_tmp.describe())

                    # plot = df_tmp.hvplot.andrews_curves(
                    #     class_column=group,
                    #     samples=samples,
                    #     alpha=alpha,
                    #     colormap=colormap
                    # )
                    # plot_layout = plot.opts(
                    #     opts.Curve(aspect=3.2,
                    #         xticks=20,
                    #         fontsize=self.font_dic,
                    #         show_grid=True,
                    #         show_legend=True,
                    #         fig_latex=True),
                    #     opts.Layout(fig_size=400)
                    # )
                    # Works using bokeh and less than 4000 datapoints
                    plot = hvplot.andrews_curves(
                        data=df_tmp,
                        class_column=group,
                        # samples=samples,
                        # alpha=alpha,
                        # colormap=colormap
                    )
                    plot_layout = plot.opts(
                        opts.Curve(aspect=3.2,
                            xticks=20,
                            fontsize=self.font_dic,
                            show_grid=True,
                            show_legend=True),
                            # fig_latex=False),
                        # opts.Layout(fig_size=400)
                    )

                    renderer = hv.Store.renderers[self.backend].instance(
                        fig='png', dpi=300)
                    latexify.set_size(True)

                    i = 0
                    if prefix is None:
                        prefix = "andrew_"
                    save = (plot_path + prefix
                            + "_" + "{:03d}".format(i))
                    while os.path.isfile(save + ".png"):
                        i = i+1
                        save = (plot_path + prefix
                                + "_" + "{:03d}".format(i))

                    print("Saving to " + save + ".png")
                    renderer.save(plot, save)
                    t2 = timer()
                    print("Saving done in {}s".format(t2-t))

                if sort:
                    while len(sorted_tuples) > 0:
                        p, v = sorted_tuples.pop()
                        in_params_2 = [p]
                        while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                            and np.abs(v/sorted_tuples[-1][1]) < 10):
                            p, v = sorted_tuples.pop()
                            in_params_2.append(p)
                        plot_helper(df_tmp_out, params=in_params_2, **kwargs)
                else:
                    plot_helper(df_tmp_out, params=params, prefix=prefix, **kwargs)


## Idea: https://www.holoviews.org/user_guide/Linking_Plots.html
