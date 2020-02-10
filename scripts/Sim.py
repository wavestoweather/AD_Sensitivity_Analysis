try:
    import scripts.loader as loader
except:
    import loader
try:
    import latexify
except:
    import scripts.latexify as latexify

import numpy as np
import pandas as pd
from progressbar import progressbar as pb
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm

from iris.analysis.cartography import rotate_pole
from mpl_toolkits.basemap import Basemap
from PIL import Image
from pylab import rcParams
import os

class Sim:
    """
    Class to load the results of a simulation and do various operations
    on it as well as provide functions to plot it. Stores different
    trajectory results in a dictionary with key the trajectory index.

    Parameters
    ----------
    data : pd.Dataframe
        Dataframe with the columns "p", "T", "w",
        "qc", "qr", "qs", "qg", "qh", "qi", "qv",
        "qiout", "qsout", "qrout", "qgout", "qhout",
        "latent_heat", "latent_cool".
    n_timesteps : Dic of int
        Number of timestep for each trajectory
    """

    def __init__(self):
        self.data = pd.DataFrame()
        self.n_timesteps = {}

    def load_file(self, filename, sep=None, nrows=None,
                 change_ref=True, refs=None):
        """
        Read a csv file and init a pandas.Dataframe with
        physical (not normalized) entries.

        Parameters
        ----------
        filename : string
            Filepath and filename of the datafile.
        sep : string
            Separator in the file.
        nrows : int, optional
            Number of rows to read from the datafile.
        change_ref : bool
            If true: Multiply all entries with reference values.
        refs : String
            Path to a file of references for transformation
        """
        df = loader.load_output(filename, sep, nrows, change_ref, refs)
        self.n_timesteps[df.trajectory.unique()[0]] = len(df.index)
        if self.data.empty:
            self.data = df
        else:
            self.data = self.data.append(df, ignore_index=True)

    def load_path(self, direc, sep=None, nrows=None,
                 change_ref=True, trajectories=None):
        """
        Read a csv file and init a pandas.Dataframe with
        physical (not normalized) entries.

        Parameters
        ----------
        direc : string
            Filepath to datafiles.
        sep : string
            Separator in the file.
        nrows : int, optional
            Number of rows to read from the datafile.
        change_ref : bool
            If true: Multiply all entries with reference values.
        trajectories : List of int
            List of trajectories to load. If none is given, it will load all
            available trajectories.
        """
        file_list = [os.path.join(direc, f) for f in os.listdir(direc)
                 if os.path.isfile(os.path.join(direc, f))]
        file_list2 = {}
        ref_list = {}
        for f in file_list:
            if "diff" in f:
                continue
            s = f.split("traj")
            s = s[1].split("_")
            if trajectories is None:
                if "reference" in f:
                    ref_list[int(s[0])] = f
                    continue
                file_list2[int(s[0])] = f
            elif int(s[0]) in trajectories:
                if "reference" in f:
                    ref_list[int(s[0])] = f
                    continue
                file_list2[int(s[0])] = f
        for t in pb(file_list2, redirect_stdout=True):
            refs = None
            if t in ref_list:
                refs = ref_list[t]
            df = loader.load_output(file_list2[t], sep, nrows, change_ref, refs)
            self.n_timesteps[df.trajectory.unique()[0]] = len(df.index)
            if self.data.empty:
                self.data = df
            else:
                self.data = self.data.append(df, ignore_index=True)

    def delete_not_mapped(self):
        """
        Delete all entries that are not within a mapped region, where
        mapped usually refers to timesteps where the WCB-criterion is
        satisfied.
        """
        self.data = self.data[self.data["MAP"] == True]
        trajectories = self.data["trajectory"].unique()
        for t in trajectories:
            tmp = self.data[self.data["trajectory"] == t]
            self.n_timesteps[t] = len(tmp.index)

    def get_n_timesteps(self):
        """
        Get the number of timesteps of the data.

        Returns
        -------
        Dic of int
            Number of timesteps (=value) for each trajectory (=key)
        """
        return self.n_timesteps

    def get_sim(self, out_param):
        """
        Return a list of values from the simulation.

        Parameters
        ----------
        out_param : string
            Key value (an output parameter).
        """
        return self.data[out_param].tolist()

    def plot(self, out_params, trajectories=None, dots=False, mapped=True,
             scatter=False, x_axis="timestep", **kwargs):
        """
        Plot results for an out_param of one or more trajectories.
        The x-axis is the timesteps or any other output parameter,
        the y-axis is the parameter. Add a bar for flagged timesteps.

        Parameters
        ----------
        out_params : list of strings
            The output parameters to plot.
        trajectories : list of int
            List of trajectories to plot.
        dots : Bool
            Plot dots every 20 seconds ie the datapoints from the
            Cosmo or ICON simulation.
        mapped : boolean
            If true: plot the region, where "MAP" is true, ie where the wcb
            criterion is fullfilled.
        scatter : boolean
            Plot a scatter plot or a line plot.
        x_axis : string
            The column to use as x-axis. Can be either "timestep" or
            "out_param" or an output parameter.
        kwargs : dict
            Keyword arguments are passed down to matplotlib.axes.Axes.plot().
        """
        def plot_helper(df, out_param, **kwargs):
            min_time = df[x_axis].unique().min()
            max_time = df[x_axis].unique().max()
            dt1 = (max_time - min_time) / 20
            dt = (max_time - min_time + dt1) / 20
            x_ticks = np.arange(min_time, max_time + dt1, dt)

            _, ax = plt.subplots()

            if len(df.trajectory.unique()) > 1:
                ax = sns.lineplot(x=x_axis, y=out_param,
                                data=df, hue="trajectory", ax=ax,
                                palette=sns.color_palette("husl", len(df.trajectory.unique())),
                                **kwargs)
            else:
                ax = sns.lineplot(x=x_axis, y=out_param,
                                data=df, ax=ax, **kwargs)
                ax.legend(df.trajectory.unique(), title="trajectory")
            # Plot dots every 20 seconds
            if dots:
                x_dots = []
                y_dots = []
                for time in df["timestep"]:
                    if time%20 == 0:
                        y_values = df[df.timestep == time][out_param].values
                        x_values = df[df.timestep == time][x_axis].values
                        x_dots.extend(x_values)
                        y_dots.extend(y_values)
                plt.scatter(x_dots, y_dots, marker='x', c="black")


            ax.set_title("Simulation of {}".format(latexify.parse_word(out_param)))
            ax.set_xticks(x_ticks)

            # Plot the area that had been flagged
            if mapped:
                df_mapped = df[df.MAP == True]
                if not df_mapped.empty:
                    min_y = df[out_param].min()
                    max_y = df[out_param].max()
                    ax.fill_between(df_mapped[x_axis], min_y, max_y,
                        facecolor="khaki", alpha=0.3)

            i = 0
            prefix = "line_"
            if scatter:
                prefix = "scatter_"
            save = ("pics/" + prefix + out_param
                    + "_" + "{:03d}".format(i) + ".png")
            while os.path.isfile(save):
                i = i+1
                save = ("pics/" + prefix + out_param
                        + "_" + "{:03d}".format(i) + ".png")

            print("Saving to " + save)
            plt.savefig(save, dpi=300)
            plt.show()
            plt.close()
        if trajectories is None:
            for out_param in out_params:
                plot_helper(self.data, out_param, **kwargs)
        else:
            tmp_df = self.data[self.data["trajectory"].isin(trajectories)]
            for out_param in out_params:
                plot_helper(tmp_df, out_param, **kwargs)