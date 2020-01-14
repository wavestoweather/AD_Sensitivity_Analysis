from argparse import ArgumentParser
from iris.analysis.cartography import rotate_pole
from mpl_toolkits.basemap import Basemap
from PIL import Image

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import numpy as np
import pandas as pd
from scipy import stats
import sys
import os


if __name__ == "__main__":
    parser = ArgumentParser(description=
            '''Load a mapped wcb netCDF file or the derivatives from a
            simulation and plot it. In case of a netCDF file, the flagged
            regions are colored differently. In case of the derivative
            file, one can plot either only flagged or not-flagged areas
            or both. Each area (flagged and not flagged) is plotted with
            normalized derivative ratios to make it easier to compare
            if needed.''')
    parser.add_argument("-i", "--input", type=str, default=None,
            help='''Path to the data folder (!) which shall be read.
            The datafiles should look like
            "wcb1000.0_traj1_start_over_MAP_20160922_02_diff_0.txt", where
            "traj" is important.''')
    parser.add_argument("-n", "--netcdf", type=str, default=None,
            help='''Path to a wcb netCDF file with flags "MAP".''')
    parser.add_argument("-t", "--trajectory", nargs='+', type=int, default=1,
            help='''A list of trajectories to consider here. Using a selected
            number of trajectories helps reducing the amount of used RAM
            and makes it easier to see anything in the image.''')
    parser.add_argument("-r", "--rotate", action="store_true", default=False,
            help='''Rotate longitude and latitude coordinates. For
            the derivatives, we assume pol_lat=90.0 and pollon=0.0 if
            no netCDF file is given.''')
    parser.add_argument("--input_param", nargs='+', type=str, default="dzeta",
            help='''A list of input parameters of the model to plot such
            as "dzeta" or "dbeta_r".''')
    parser.add_argument("--output_param", nargs='+', type=str, default="T",
            help='''A list of output parameters of the model to plot such
            as "T", "qr" or "p".''')
    parser.add_argument("--max_time", type=float, default=None,
            help='''Maximum time in seconds to plot.''')
    parser.add_argument("-e", "--epsilon", type=float, default=0.0,
            help='''Value to filter values with. Default is 0.0.''')
    parser.add_argument("-f", "--filter", action="store_true", default=False,
            help='''Filter values smaller than epsilon if used.''')
    parser.add_argument("--norm", action="store_true", default=True,
            help='''Normalize the derivatives for an area (consecutive
            timings with or without flag).''')

    args = parser.parse_args()

    try:
        from loader import load_mult_derivates_directory, load_nc, rotate_df
    except:
        from scripts.loader import load_mult_derivates_directory, load_nc, rotate_df
    from pylab import rcParams


    rcParams['figure.figsize'] = (16,10)

    df = None
    net_df = None
    pollon = 0.0
    pollat = 90.0

    if args.netcdf is not None:
        if args.rotate:
            net_df, pollat, pollon = load_nc(inp=args.netcdf,
                                             get_pol=True)
            rotate_df(net_df, pollon, pollat, "lon", "lat")
        else:
            net_df = load_nc(inp=args.netcdf,
                             get_pol=False)
        #TODO Remove all trajectories that are not needed


    if args.input is not None:
        df = load_mult_derivates_directory(direc=args.input,
                                            filt=args.filter,
                                            EPSILON=args.epsilon,
                                            trajectories=args.trajectory,
                                            suffix=None)
        if args.rotate:
            rotate_df(df, pollon, pollat, "LONGITUDE", "LATITUDE")

        if args.norm:


    # Plot data
    if net_df is not None:
        for in_par in args.input_param:
            for out_par in args.output_param:
                plot_netcdf(df=net_df,
                            max_time=args.max_time,
                            in_par=in_par,
                            out_par=out_par)

    if df is not None:




