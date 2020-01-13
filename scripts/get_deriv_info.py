from argparse import ArgumentParser
from iris.analysis.cartography import rotate_pole
import numpy as np
import pandas as pd
import xarray as xr
from loader import load_mult_derivates_directory


def norm_derivatives(df):
    """
    Normalize the derivatives for every timestep individually such that
    values are between [-1, 1] for big negative and positive impacts.

    .. math::
        \\text{norm}(x, t) = \\text{sign}(x(t)) \cdot \\frac{|x(t)| - \min(x(t))}{\max(x(t)) - \min(x(t))}

    This is done inplace.

    Parameters
    ----------
    df : pandas.DataFrame
        Columns are "timestep", "trajectory", "out_param", "in_param" and "deriv".

    """

    for traj in df["trajectory"].unique():
        df_traj = df.loc[df.trajectory == traj]
        for out_param in df["out_param"].unique():
            df_param = df_traj.loc[df_traj.out_param == out_param]
            for t in df["timestep"].unique():
                df_t = df_param.loc[df_param.timestep == t]
                deriv_min = df_t["deriv"].min()
                deriv_max = df_t["deriv"].max()
                df_t["deriv"] = ( np.sign(df_t["deriv"])
                                  * ( np.abs(df_t["deriv"]) - deriv_min )
                                  / ( deriv_max - deriv_min ) )
                df.update(df_t)


if __name__ == "__main__":
    parser = ArgumentParser(description=
            '''Load derivatives of a simulation and store that to a csv
            file after being filtered and transformed if needed.''')
    parser.add_argument("-i", "--input", type=str, default=None,
            help='''Path to folder with all files to load.''')
    parser.add_argument("-o", "--output", type=str, default="file.csv",
            help='''Directory and name to save the normed csv.''')
    parser.add_argument("-e", "--epsilon", type=float, default=0.0,
            help='''Value to filter values with. Default is 0.0.''')
    parser.add_argument("-f", "--filter", action="store_true", default=False,
            help='''Filter values smaller than epsilon if used.''')
    parser.add_argument("-t", "--trajectory", nargs='+', type=int, default=1,
            help='''A list of trajectories to consider here. Using a selected
            number of trajectories helps reducing the amount of used RAM.''')
    parser.add_argument("-c", "--csv", type=str, default=None,
            help='''Load an already filtered csv file and just norm the derivatives.
            Deactivates -t, -f and -e.''')
    parser.add_argument("-n", "--netcdf", type=str, default=None,
            help='''Path to a netCDF file that is used to include longitude and latitude.''')
    parser.add_argument("-r", "--rotate", action="store_true", default=False,
            help='''Rotate longitude and latitude coordinates.''')
    args = parser.parse_args()

    if args.csv is None:
        df_data = load_mult_derivates_directory(
                direc=args.input,
                filt=args.filter,
                EPSILON=args.epsilon,
                trajectories=args.trajectory)
    else:
        df_data = pd.read_csv(args.csv)
        # df_data = pd.read_csv(args.csv, nrows=1000000, names=["deriv", "in_param", "out_param", "timestep", "trajectory"], skiprows=100000)
    print("Max and min deriv: {}, {}".format(df_data["deriv"].max(), df_data["deriv"].min()))

    norm_derivatives(df_data)
    df_data = df_data.loc[:, ~df_data.columns.str.contains('^Unnamed')]

    print(df_data)
    print("Max and min deriv: {}, {}".format(df_data["deriv"].max(), df_data["deriv"].min()))
    if not args.netcdf is None:
        df_data = transform_df(df_data, args.netcdf, args.rotate)
    df_data.to_csv(args.output)


