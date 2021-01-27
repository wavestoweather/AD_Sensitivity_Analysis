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

from pylab import rcParams
import os

from multiprocessing import Pool
from itertools import repeat
from progressbar import progressbar as pb


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
    clusternames : Dic of lsit of string
        Keys are output parameters where the value is a list of column names
        that hold a clustering assignment.
    pool : multiprocessing.Pool
        Used to work on data in parallel
    threads : int
            Number of threads to use.
    """
    data = {}
    n_timesteps = 0
    cluster_names = {}
    pool = None
    threads = 0

    def __init__(self, direc, filt=True, file_list=None,
                 EPSILON=0.0, trajectories=[1], suffix=None, threads=8):
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
        threads : int
            Number of threads to use for plotting.

        """
        if threads is None:
            self.pool = None
        else:
            self.pool = Pool(processes=threads)
        self.threads = threads
        self.data = loader.load_mult_derivates_direc_dic(
            direc=direc,
            filt=filt,
            EPSILON=EPSILON,
            trajectories=trajectories,
            suffix=suffix,
            pool=self.pool,
            file_list2=file_list
        )
        df = list(self.data.values())[0]
        self.n_timesteps = len(df.index)
        self.cluster_names = {}

    def to_parquet(self, f_name, compression, add_columns=None, low_memory=False, attr=None):
        import dask.dataframe as dd
        if attr is not None:
            attributes = loader.parse_attr(attr)
        if low_memory:
            for k in self.data:
                tmp_df = self.data[k]
                if add_columns is not None:
                    for col in add_columns:
                        tmp_df[col] = add_columns[col]
                append = not os.listdir(f_name)
                append = not append
                if attr is not None:
                    for key in attributes:
                        if key == "Global attributes":
                            tmp_df.attrs = attributes[key]
                        else:
                            for col in attributes[key]:
                                tmp_df[col].attrs = attributes[key][col]
                dd.from_pandas(tmp_df, chunksize=3000000).to_parquet(f_name, append=append, ignore_divisions=append, compression=compression)
        else:
            tmp_df = None
            for k in self.data:
                if tmp_df is not None:
                    tmp_df = tmp_df.append(self.data[k], ignore_index=True)
                else:
                    tmp_df = self.data[k]
            if add_columns is not None:
                for col in add_columns:
                    tmp_df[col] = add_columns[col]
            append = not os.listdir(f_name)
            append = not append
            if attr is not None:
                for key in attributes:
                    if key == "Global attributes":
                        tmp_df.attrs = attributes[key]
                    else:
                        for col in attributes[key]:
                            tmp_df[col].attrs = attributes[key][col]
            dd.from_pandas(tmp_df, chunksize=3000000).to_parquet(f_name, append=append, ignore_divisions=append, compression=compression)

    def get_netcdf_ready_data(self, add_columns=None, dropna=False, attr=None):
        import xarray as xr

        df = None
        for k in self.data:
            tmp_df = self.data[k]
            if add_columns is not None:
                for col in add_columns:
                    tmp_df[col] = add_columns[col]
            if df is None:
                df = tmp_df
            else:
                df = df.append(tmp_df)
        df["Output Parameter"] = df["Output Parameter"].astype(str)
        df["Output Parameter"].attrs = {
            "standard_name": "output_parameter",
            "long_name": "Output parameter for sensitivities"}
        df["instance_id"] = df["instance_id"].astype(str)
        df["instance_id"].attrs = {
            "standard_name": "instance_id",
            "long_name": "Instance ID"}

        if dropna:
            ds_complete = xr.Dataset.from_dataframe(df.set_index(
                ["Output Parameter", "ensemble", "trajectory", "time"]).dropna())
        else:
            ds_complete = xr.Dataset.from_dataframe(df.set_index(
                ["Output Parameter", "ensemble", "trajectory", "time"]))
        # Add attributes where applicable
        if attr is not None:
            attributes = loader.parse_attr(attr)
            for key in attributes:
                if key == "Global attributes":
                    ds_complete.attrs = attributes[key]
                else:
                    for col in attributes[key]:
                        if col in ds_complete:
                            ds_complete[col].attrs = attributes[key][col]
        return ds_complete

    def get_dataframe(self, add_columns=None, dropna=False):
        df = None
        for k in self.data:
            tmp_df = self.data[k]
            if add_columns is not None:
                for col in add_columns:
                    tmp_df[col] = add_columns[col]
            if df is None:
                df = tmp_df
            else:
                df = df.append(tmp_df)
        df["Output Parameter"] = df["Output Parameter"].astype(str)
        df["Output Parameter"].attrs = {
            "standard_name": "output_parameter",
            "long_name": "Output parameter for sensitivities"}
        df["instance_id"] = df["instance_id"].astype(str)
        df["instance_id"].attrs = {
            "standard_name": "instance_id",
            "long_name": "Instance ID"}

        if dropna:
            df_complete = df.set_index(
                ["Output Parameter", "ensemble", "trajectory", "time"]).dropna()
        else:
            df_complete = df.set_index(
                ["Output Parameter", "ensemble", "trajectory", "time"])
        return df_complete

    def to_netcdf(self, f_name, add_columns=None, dropna=False, met3d=False, attr=None):
        """
        Store the data to netcdf-4 files where the index consists of
        (timestep, trajectory, LONGITUDE, LATITUDE), and columns are the derivatives.

        Parameters
        ----------
        f_name : string
                 Path and first name of the netcdf-file
        """
        import xarray as xr

        df = None
        for k in self.data:
            tmp_df = self.data[k]
            if add_columns is not None:
                for col in add_columns:
                    tmp_df[col] = add_columns[col]
            if df is None:
                df = tmp_df
            else:
                df = df.append(tmp_df)
        if met3d:
            df["Output Parameter"] = df["Output Parameter"].astype(str)
            df["Output Parameter"].attrs = {
                "standard_name": "output_parameter",
                "long_name": "Output parameter for sensitivities"}
            if dropna:
                ds_complete = xr.Dataset.from_dataframe(df.set_index(
                    ["Output Parameter", "ensemble", "trajectory", "time"]).dropna())
            else:
                ds_complete = xr.Dataset.from_dataframe(df.set_index(
                    ["Output Parameter", "ensemble", "trajectory", "time"]))
        else:
            if dropna:
                ds_complete = xr.Dataset.from_dataframe(df.set_index(
                    ["Output Parameter", "timestep", "trajectory"]).dropna())
            else:
                ds_complete = xr.Dataset.from_dataframe(df.set_index(
                    ["Output Parameter", "timestep", "trajectory"]))

        # Choose encoding
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds_complete.data_vars}

        # Add attributes where applicable
        if attr is not None:
            attributes = loader.parse_attr(attr)
            for key in attributes:
                if key == "Global attributes":
                    ds_complete.attrs = attributes[key]
                else:
                    for col in attributes[key]:
                        ds_complete[col].attrs = attributes[key][col]
        ds_complete.to_netcdf(
            f_name + "_derivs.nc_wcb",
            encoding=encoding,
            compute=True,
            engine="netcdf4",
            format="NETCDF4",
            mode="w")

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

    def add_param_values(self, df):
        """
        Add to dataframe a column header with the values.

        Parameters
        ----------
        df : pandas.Dataframe
            Dataframe with output parameters of simulation

        """
        cols = []
        for col in df:
            if col in ["LONGITUDE", "LATITUDE", "MAP", "dp2h",
                       "conv_400", "conv_600", "slan_400", "slan_600",
                       "lon", "lat", "WCB_flag", "instance_id"]:
                continue
            cols.append(col)
        for k in self.data:
            self.data[k] = self.data[k].merge(df[cols], how='right')

    def shift_time(self, flag, debug=False):
        """
        Shift the time axis such that at t=0 the first occurence of flag appears.

        Parameters
        ----------
        flag : string
               Name of the (boolean type) column, ie conv_400, slan_600, conv_600, slan_400.
        """
        start_time = None
        for k in self.data:
            if start_time is None:
                start_time = self.data[k].loc[self.data[k][flag] == True]["timestep"].min()
                if debug:
                    print(f"Using start time {start_time} s for flag {flag}")
            if np.isnan(start_time):
                print(f"No region with {flag} found")
                return
            self.data[k]["timestep"] = self.data[k]["timestep"] - start_time
            if debug:
                new_min_time = self.data[k].loc[self.data[k][flag] == True]["timestep"].min()
                print(f"New start time is {new_min_time} for {flag} (if this is not 0, something went wrong)")
                assert new_min_time == 0

    @staticmethod
    def parallel_ratio(df, k):
        """
        Used to calculate ratios in parallel using multiprocessing.
        Finds the maximum derivative for a given output parameter.

        Parameters
        ----------
        df : pandas.Dataframe
            The dataframw with derivatives.
        k : The output parameter for which the derivatives are calculated for.

        Returns
        -------
        float, string
            Maximum absolute derivative and output parameter.
        """
        denom = -1.0
        for deriv in df:
            if deriv[0] != 'd':
                continue
            d = np.abs(df[deriv]).max()
            if d > denom:
                denom = d
        return denom, k

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
                if deriv[0] != 'd' or deriv == "dp2h" or deriv == "depo":
                    continue
                d = np.abs(df[deriv]).max()
                if d > denom:
                    denom = d
            return denom

        if isinstance(key, str):
            denom = get_max_denom(self.data[key])
            if denom == 0:
                return
            for deriv in self.data[key]:
                if deriv[0] != 'd' or deriv == "dp2h" or deriv == "depo":
                    continue
                self.data[key][deriv] = self.data[key][deriv]/denom
            return
        elif key is None:
            if self.pool is None:
                for k in self.data.keys():
                    denom = get_max_denom(self.data[k])
                    if denom == 0:
                        continue
                    for deriv in self.data[k]:
                        if deriv[0] != 'd' or deriv == "dp2h" or deriv == "depo":
                            continue
                        self.data[k][deriv] = self.data[k][deriv]/denom
            else:
                for denom, k in pb( self.pool.starmap(self.parallel_ratio,
                    zip([self.data[k] for k in self.data.keys()], self.data.keys())) ):
                    if denom == 0:
                        continue
                    for deriv in self.data[k]:
                        if deriv[0] != 'd' or deriv == "dp2h" or deriv == "depo":
                            continue
                        self.data[k][deriv] = self.data[k][deriv]/denom
        else:
            if self.pool is None:
                for k in key:
                    denom = get_max_denom(self.data[k])
                    if denom == 0:
                        continue
                    for deriv in self.data[k]:
                        if deriv[0] != 'd' or deriv == "dp2h" or deriv == "depo":
                            continue
                        self.data[k][deriv] = self.data[k][deriv]/denom
            else:
                for denom, k in pb( self.pool.starmap(self.parallel_ratio,
                    zip([self.data[k] for k in key], key)) ):
                    if denom == 0:
                        continue
                    for deriv in self.data[k]:
                        if deriv[0] != 'd' or deriv == "dp2h" or deriv == "depo":
                            continue
                        self.data[k][deriv] = self.data[k][deriv]/denom

    @staticmethod
    def parallel_plot_orders(df, out_param, mapped, scatter,
        in_params, x_axis, kwargs):
        """
        Helper method to plot in parallel.
        """
        def plot_helper(df, in_params, out_param, **kwargs):
            min_time = df[x_axis].min()
            max_time = df[x_axis].max()
            dt1 = (max_time - min_time) / 20
            dt = (max_time - min_time + dt1) / 20
            x_ticks = np.arange(min_time, max_time + dt1, dt)

            _, ax = plt.subplots()
            df_tmp = df[in_params+[x_axis]]
            df_tmp = df_tmp.melt(x_axis, var_name="Derivatives",
                                    value_name="Derivative Ratio")
            if scatter:
                ax = sns.scatterplot(x=x_axis, y="Derivative Ratio",
                                        data=df_tmp, hue="Derivatives", ax=ax, linewidth=0,
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
                    xmin = df_mapped[x_axis].min()
                    xmax = df_mapped[x_axis].max()
                    plt.axvspan(xmin=xmin, xmax=xmax,
                        facecolor="khaki", alpha=0.3)

            # Change the limits for the y-axis because sometimes that
            # can be off and it is hard to see anything.
            min_y = df_tmp["Derivative Ratio"].min()
            max_y = df_tmp["Derivative Ratio"].max()
            plt.ylim(min_y - min_y/10, max_y + max_y/10)
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

        if df is None:
            return
        if df.empty:
            return

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
            Keyword arguments are passed down to seaborn.scatterplott() or
            seaborn.lineplot() for the derivative plots.
        """
        if out_params is None:
            out_params = self.data.keys()
        elif isinstance(out_params, str):
            out_params = [out_params]
        self.pool.starmap(self.parallel_plot_orders,
            zip([self.data[p] for p in out_params],
            out_params, repeat(mapped), repeat(scatter), repeat(in_params),
            repeat(x_axis), repeat(kwargs)))

    @staticmethod
    def parallel_plot_mapped(df, out_param, in_params, kind, x_label, kwargs):
        def plot_helper(df, in_params, out_param, **kwargs):
            df_tmp = df[in_params+["MAP"]]
            df_tmp = df_tmp.melt("MAP", var_name="Derivatives",
                                 value_name="Derivative Ratio")
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
            ax.set_title("Deriv. Ratio of {}".format(latexify.parse_word(out_param)))
            ax.set_ylabel("Derivative ratio")
            ax.set_xlabel(x_label)
            plt.ticklabel_format(style="scientific", axis="y", scilimits=(0,0))

            # Change the limits for the y-axis because sometimes that
            # can be off and it is hard to see anything.
            min_y = df_tmp["Derivative Ratio"].min()
            max_y = df_tmp["Derivative Ratio"].max()
            plt.ylim(min_y - min_y/10, max_y + max_y/10)

            i = 0
            save = ("pics/" + kind + "_MAP_" + out_param
                    + "_" + "{:03d}".format(i) + ".png")
            while os.path.isfile(save):
                i = i+1
                save = ("pics/" + kind + "_MAP_" + out_param
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

        if df is None:
            return
        if df.empty:
            return

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


    def plot_mapped(self, out_params=None, in_params=None, kind="violin",
                    x_label="WCB criterion", **kwargs):
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
            “point”, “bar”, “strip”, “box”, or “boxen”. Is passed
            to seaborn.catplot(..).
        x_label : String
            This will be plotted on the x-axis.
        kwargs : dict
            Keyword arguments are passed down to seaborn.catplot(..) for
            the derivative plots.
        """
        if out_params is None:
            out_params = self.data.keys()
        elif isinstance(out_params, str):
            out_params = [out_params]
        if self.pool is not None:
            self.pool.starmap(self.parallel_plot_mapped,
                zip([self.data[out_param] for out_param in out_params],
                out_params, repeat(in_params), repeat(kind), repeat(x_label),
                repeat(kwargs)))
        else:
            def plot_helper(df, in_params, out_param, **kwargs):
                df_tmp = df[in_params+["MAP"]]
                print(out_param)
                print(df_tmp)
                df_tmp = df_tmp.melt("MAP", var_name="Derivatives",
                                    value_name="Derivative Ratio")
                df_tmp["Derivatives"] = df_tmp["Derivatives"].apply(latexify.parse_word)
                print(df_tmp)
                print(df_tmp["Derivative Ratio"].unique())
                print(df_tmp["MAP"].unique())
                print("######################################################")
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

                ax.set_title("Deriv. Ratio of {}".format(latexify.parse_word(out_param)))
                ax.set_ylabel("Derivative ratio")
                ax.set_xlabel(x_label)
                plt.ticklabel_format(style="scientific", axis="y", scilimits=(0,0))

                # Change the limits for the y-axis because sometimes that
                # can be off and it is hard to see anything.
                min_y = df_tmp["Derivative Ratio"].min()
                max_y = df_tmp["Derivative Ratio"].max()
                plt.ylim(min_y - min_y/10, max_y + max_y/10)

                i = 0
                save = ("pics/" + kind + "_MAP_ " + out_param
                        + "_" + "{:03d}".format(i) + ".png")
                while os.path.isfile(save):
                    i = i+1
                    save = ("pics/" + kind + "_MAP_ " + out_param
                            + "_" + "{:03d}".format(i) + ".png")

                print("Saving to " + save)
                plt.savefig(save, dpi=300)
                plt.show()
                plt.close()

            for out_param in out_params:
                df = self.data[out_param]
                if df is None:
                    return
                if df.empty:
                    return

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
