import pandas as pd
from progressbar import progressbar as pb


def get_traj(df, traj=1):
    """
    Given a pandas dataframe and a trajectory number or a list of numbers,
    get only those trajectories and discard the others.

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
   # conditions = [df_out.in_param == max_param for max_param in max_params]
   # conditions2 = conditions[0]
   # for cond in conditions:
   #     conditions2 = conditions2 or cond
   # return df_out[conditions2]


if __name__ == "__main__":
    res = pd.read_csv("sb_ice_wcb272280_filt_zero_30_start_over_20160922_00.csv")
    res = get_traj(df=res, traj=1)
    res.to_csv("sb_ice_wcb272280_filt_zero_30_start_over_20160922_00_traj1_all.csv")
    # res = filt_highest_in(df=res, out_param="qr", n=10)
    # res.to_csv("sb_ice_wcb272280_filt_zero_30_start_over_20160922_00_traj1_qr_top10.csv")
