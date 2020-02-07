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
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm

from iris.analysis.cartography import rotate_pole
from mpl_toolkits.basemap import Basemap
from PIL import Image
from pylab import rcParams
import os

class Deriv:
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
    """
    data = {}
    n_timesteps = 0

    def __init__(self, direc, filt=True,
                 EPSILON=0.0, trajectories=[1], suffix=None):
        """
        Init class by loading the data from the given path.

        Parameters
        ----------
        direc : string
            A path to a directory wit a list of files to read.
        filt : bool
            Filter the data with values smaller than EPSILON if true.
        EPSILON : float
            If filt is true, filter values smaller than EPSILON out.
        trajectories : list of int
            A list of trajectories to read in.
        suffix : string
            The suffix of the filenames before '_diff_xx.txt'. If none is
            given, the method tries to automatically detect it.

        """
        self.data = loader.load_mult_derivates_direc_dic(
            direc, filt, EPSILON, trajectories, suffix
        )
        df = list(self.data.values())[0]
        self.n_timesteps = len(df.index)

    def delete_not_mapped(self):
        """
        Delete all entries that are not within a mapped region, where
        mapped usually refers to timesteps where the WCB-criterion is
        satisfied.
        """
        for key in self.data:
            self.data[key] = self.data[key][self.data[key]["MAP"] == True]
        df = list(self.data.values())[0]
        self.n_timesteps = len(df.index)

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

    def add_param_values(self, header, values, key=None):
        """
        Add to dataframe with key key another column header with the values.

        Parameters
        ----------
        header : string
            The name of the new column
        values : List
            The values to add. Must have the same size as the dataframes.
        key : string or list of strings
            The key of the dataframe. If None is given, the procedure
            is applied on all dataframes.
        """
        if len(values) != self.n_timesteps:
            print("Cannot add column {}.".format(header))
            print("Number of timesteps needed: {} But got {}".format(
                self.n_timesteps, len(values)))
            return
        if isinstance(key, str):
            self.data[key][header] = values
        elif key is None:
            for k in self.data:
                self.data[k][header] = values
        else:
            for k in key:
                self.data[k][header] = values

    def calculate_ratios(self, key=None):
        """
        Calculate the ratio of the derivatives for every timestep.

        Parameters
        ----------
        key : String or list of strings
            Output parameter to calculate the ratio for. If None, then
            it will be done for all available parameters.
        """
        def get_max_denom(df):
            denom = -1.0
            for deriv in df:
                if deriv[0] != 'd':
                    continue
                d = np.abs(df[deriv]).max()
                if d > denom:
                    denom = d
            return denom

        if isinstance(key, str):
            denom = get_max_denom(self.data[key])
            for deriv in self.data[key]:
                if deriv[0] != 'd':
                    continue
                self.data[key][deriv] = self.data[key][deriv]/denom
        elif key is None:
            for k in self.data:
                denom = get_max_denom(self.data[k])
                for deriv in self.data[k]:
                    if deriv[0] != 'd':
                        continue
                    self.data[k][deriv] = self.data[k][deriv]/denom
        else:
            for k in key:
                denom = get_max_denom(self.data[k])
                for deriv in self.data[k]:
                    if deriv[0] != 'd':
                        continue
                    self.data[k][deriv] = self.data[k][deriv]/denom

    def plot_same_orders(self, out_params=None, mapped=True, scatter=False,
                         in_params=None, x_axis="timestep", **kwargs):
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
        kwargs : dict
            Keyword arguments are passed down to matplotlib.axes.Axes.plot() for
            the derivative plots.
        """
        if out_params is None:
            out_params = self.data.keys()
        elif isinstance(out_params, str):
            out_params = [out_params]

        def plot_helper(df, in_params, out_param, **kwargs):
            min_time = df[x_axis].min()
            max_time = df[x_axis].max()
            dt = (max_time - min_time) / 20
            dt = (max_time - min_time + dt) / 20
            x_ticks = np.arange(min_time, max_time + 19, dt)

            _, ax = plt.subplots()
            df_tmp = df[in_params+[x_axis]]
            df_tmp = df_tmp.melt(x_axis, var_name="Derivatives", value_name="Derivative Ratio")
            if scatter:
                ax = sns.scatterplot(x=x_axis, y="Derivative Ratio",
                                     data=df_tmp, hue="Derivatives", ax=ax,
                                     **kwargs)
            else:
                ax = sns.lineplot(x=x_axis, y="Derivative Ratio",
                                 data=df_tmp, hue="Derivatives", ax=ax,
                                 **kwargs)
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
                    min_y = df_tmp["Derivative Ratio"].min()
                    max_y = df_tmp["Derivative Ratio"].max()
                    ax.fill_between(df_mapped[x_axis], min_y, max_y,
                        facecolor="khaki", alpha=0.3)
            # Change the limits for the y-axis because sometimes that
            # can be off and it is hard to see anything.
            min_y = df_tmp["Derivative Ratio"].min()
            max_y = df_tmp["Derivative Ratio"].max()
            plt.ylim(min_y - min_y/10, max_y + max_y/10)
            i = 0
            prefix = "line_"
            if scatter:
                prefix = "scatter_"
            save = ("pics/" + prefix + x_axis + "_ " + out_param
                    + "_" + "{:03d}".format(i) + ".png")
            while os.path.isfile(save):
                i = i+1
                save = ("pics/" + prefix + x_axis + "_ " + out_param
                        + "_" + "{:03d}".format(i) + ".png")

            print("Saving to " + save)
            plt.savefig(save, dpi=300)
            plt.show()
            plt.close()

        for out_param in out_params:
            df = self.data[out_param]
            if df is None:
                continue
            if df.empty:
                continue

            if in_params is None:
                in_params_tmp = list(df)
                in_params = []
                for i in range(len(in_params_tmp)):
                    if in_params_tmp[i][0] == 'd':
                        in_params.append(in_params_tmp[i])

            # Sort the derivatives
            sorted_tuples = []
            for in_p in in_params:
                value = np.abs(df[in_p].min())
                if np.abs(df[in_p].max()) > value:
                    value = np.abs(df[in_p].max())
                if value != 0 and not np.isnan(value):
                    sorted_tuples.append((in_p, value))
            sorted_tuples.sort(key=lambda tup: tup[1])

            # Plot them
            while len(sorted_tuples) > 0:
                p, v = sorted_tuples.pop()
                in_params_2 = [p]
                while (len(sorted_tuples) > 0 and sorted_tuples[-1][1] > 0
                    and np.abs(v/sorted_tuples[-1][1]) < 10):
                    p, v = sorted_tuples.pop()
                    in_params_2.append(p)
                plot_helper(df, in_params=in_params_2, out_param=out_param, **kwargs)