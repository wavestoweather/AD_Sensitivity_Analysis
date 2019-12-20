from argparse import ArgumentParser
import pandas as pd
from progressbar import progressbar as pb


def get_traj(df, traj=1):
    """
    Given a pandas dataframe and a trajectory number or a list of numbers,
    get only those trajectories and discard the others.

    Parameters
    ----------
    df : pandas.Dataframe
        Dataframe with a column "trajectory".
    traj : int or list of int
        Trajectory index to filter for.

    Returns
    -------
    pandas.Dataframe
        Dataframe with only certain trajectories selected.

    """
    if isinstance(traj, list):
        return df.loc[df["trajectory"] in traj]
    else:
        return df.loc[df["trajectory"] == traj]


def filt_highest_in(df, out_param, n=3):
    """
    Given a pandas dataframe and an out_param, get the n most important
    in_params (sum over all timesteps) and discard all others. Deletes
    other out_params as well.

    Paramters
    ---------
    df : pandas.Dataframe
        Dataframe with columns "in_param" and "deriv".
    out_param : string
        The output parameter to consider here such as "qr".
    n : int
        Get the n most important derivatives.

    Returns
    -------
    pandas.Dataframe
        Dataframe with only out_param and the n most important 
        in_param.

    """
    max_deriv = [0.0 for i in range(n)]
    max_params = [" " for i in range(n)]
    df_out = df.loc[df["out_param"] == out_param]
    in_params = df_out["in_param"].unique()

    for in_param in pb(in_params):
        this_df = df_out.loc[df_out["in_param"] == in_param]
        this_sum = this_df["deriv"].sum()
        for i in range(n):
            if this_sum > max_deriv[i]:
                max_deriv[i] = this_sum
                max_params[i] = in_param
    df_filtered = df_out.loc[df_out.in_param == max_params[0]]
    for i in range(1, len(max_params)):
        df_filtered = pd.concat(
		[df_filtered, df_out.loc[df_out.in_param == max_params[i]]], 
		ignore_index=True)
    return df_filtered


if __name__ == "__main__":
    parser = ArgumentParser(description=
            '''Filter a csv file with a given trajectory and get the n most important derivatives
            for a given out_param.''')
    parser.add_argument("-i", "--input", type=str, default="sb_ice_wcb272280_filt_zero_30_start_over_20160922_00.csv",
            help='''Path to csv file.''')
    parser.add_argument("-n", type=int, default=-1, 
            help='''The number of derivatives on output. If n=-1 no filtering is being done.''')
    parser.add_argument("-o", "--output", type=str, default=None, required=True,
            help='''Path and name where to store the filtered file.''')
    parser.add_argument("-p", "--param", type=str, default=None,
            help='''The out parameter to filter for.''')
    parser.add_argument("-t", "--trajectory", type=int, default=1,
            help='''The trajectory to filter for.''')
    args = parser.parser_args()
    if args.trajectory < 0:
        raise ValueError("You used {} for trajectories but it should be positive.".format(args.trajectory))

    res = pd.read_csv(args.input)
    res = get_traj(df=res, traj=args.trajectory)
    if args.p is None:
        res.to_csv(args.output)
    else:
        if args.n < 1:
            raise ValueError("You used {} for the number of derivatives.".format(args.n))
        res = filt_highest_in(df=res, out_param=args.param, n=args.n)
        res.to_csv(args.output)

