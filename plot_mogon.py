# Usage
# plot_mogon <path to parquet files> <key> <out_param>
# Key can be: 'Misc', 'Rain', 'Cloud', 'Graupel', 'Hail', 'Ice', 'Snow'
# out_param can be:
#   p,T,w,S,qc,qr,qv,Nc,Nr,Nv,qi,Ni,vi,
#   qs,Ns,qg,Ng,qh,Nh,qiout,qsout,qrout,qgout,qhout,
#   latent_heat,latent_cool,Niout,Nsout,Nrout,Ngout,Nhout
#
import numpy as np
from pylab import rcParams
from scripts.Deriv_dask import Deriv_dask
import scripts.latexify as latexify
from timeit import default_timer as timer
import sys

rcParams['figure.figsize'] = (16,10)
latexify.set_size(beamer=True)
kwargs = {"alpha": 0.4}
trajectories = None
in_params_dic = {"Misc":
            ["da_1", "da_2", "de_1", "de_2", "dd", "dN_c", "dgamma",
            "dbeta_c", "dbeta_r", "ddelta1", "ddelta2", "dzeta",
            "drain_gfak", "dcloud_k_au", "dcloud_k_sc", "dkc_autocon",
            "dinv_z"],
            "Rain":
            ["drain_a_geo", "drain_b_geo", "drain_min_x", "drain_max_x",
            "drain_sc_theta_q", "drain_sc_delta_q", "drain_sc_theta_n",
            "drain_sc_delta_n", "drain_s_vel", "drain_a_vel", "drain_b_vel",
            "drain_rho_v", "drain_c_z", "drain_sc_coll_n", "drain_cmu0",
            "drain_cmu1", "drain_cmu2", "drain_cmu3", "drain_cmu4",
            "drain_cmu5", "drain_alpha", "drain_beta", "drain_gamma",
            "drain_nu", "drain_g1", "drain_g2", "drain_mu", "drain_nm1",
            "drain_nm2", "drain_nm3", "drain_q_crit_c", "drain_d_crit_c",
            "drain_ecoll_c", "drain_cap", "drain_a_ven", "drain_b_ven",
            "drain_c_s", "drain_a_f", "drain_b_f", "drain_alfa_n",
            "drain_alfa_q", "drain_lambda", "drain_vsedi_min",
            "drain_vsedi_max"],
            "Cloud":
            ["dcloud_a_geo", "dcloud_b_geo", "dcloud_min_x",
            "dcloud_max_x", "dcloud_sc_theta_q", "dcloud_sc_delta_q",
            "dcloud_sc_theta_n", "dcloud_sc_delta_n", "dcloud_s_vel",
            "dcloud_a_vel", "dcloud_b_vel", "dcloud_rho_v", "dcloud_c_z",
            "dcloud_sc_coll_n", "dcloud_cmu0", "dcloud_cmu1", "dcloud_cmu2",
            "dcloud_cmu3", "dcloud_cmu4", "dcloud_cmu5", "dcloud_alpha",
            "dcloud_beta", "dcloud_gamma", "dcloud_nu", "dcloud_g1",
            "dcloud_g2", "dcloud_mu", "dcloud_nm1", "dcloud_nm2",
            "dcloud_nm3", "dcloud_q_crit_c", "dcloud_d_crit_c",
            "dcloud_ecoll_c", "dcloud_cap", "dcloud_a_ven", "dcloud_b_ven",
            "dcloud_c_s", "dcloud_a_f", "dcloud_b_f", "dcloud_alfa_n",
            "dcloud_alfa_q", "dcloud_lambda", "dcloud_vsedi_min",
            "dcloud_vsedi_max"],
            "Graupel":
            ["dgraupel_a_geo", "dgraupel_b_geo", "dgraupel_min_x",
            "dgraupel_max_x", "dgraupel_sc_theta_q", "dgraupel_sc_delta_q",
            "dgraupel_sc_theta_n", "dgraupel_sc_delta_n", "dgraupel_s_vel",
            "dgraupel_a_vel", "dgraupel_b_vel", "dgraupel_rho_v",
            "dgraupel_c_z", "dgraupel_sc_coll_n", "dgraupel_cmu0",
            "dgraupel_cmu1", "dgraupel_cmu2", "dgraupel_cmu3", "dgraupel_cmu4",
            "dgraupel_cmu5", "dgraupel_alpha", "dgraupel_beta",
            "dgraupel_gamma", "dgraupel_nu", "dgraupel_g1", "dgraupel_g2",
            "dgraupel_mu", "dgraupel_nm1", "dgraupel_nm2", "dgraupel_nm3",
            "dgraupel_q_crit_c", "dgraupel_d_crit_c", "dgraupel_ecoll_c",
            "dgraupel_cap", "dgraupel_a_ven", "dgraupel_b_ven", "dgraupel_c_s",
            "dgraupel_a_f", "dgraupel_b_f", "dgraupel_alfa_n", "dgraupel_alfa_q",
            "dgraupel_lambda", "dgraupel_vsedi_min", "dgraupel_vsedi_max"],
            "Hail":
            ["dhail_a_geo", "dhail_b_geo", "dhail_min_x", "dhail_max_x",
            "dhail_sc_theta_q", "dhail_sc_delta_q", "dhail_sc_theta_n",
            "dhail_sc_delta_n", "dhail_s_vel", "dhail_a_vel", "dhail_b_vel",
            "dhail_rho_v", "dhail_c_z", "dhail_sc_coll_n", "dhail_cmu0",
            "dhail_cmu1", "dhail_cmu2", "dhail_cmu3", "dhail_cmu4",
            "dhail_cmu5", "dhail_alpha", "dhail_beta", "dhail_gamma",
            "dhail_nu", "dhail_g1", "dhail_g2", "dhail_mu", "dhail_nm1",
            "dhail_nm2", "dhail_nm3", "dhail_q_crit_c", "dhail_d_crit_c",
            "dhail_ecoll_c", "dhail_cap", "dhail_a_ven", "dhail_b_ven",
            "dhail_c_s", "dhail_a_f", "dhail_b_f", "dhail_alfa_n",
            "dhail_alfa_q", "dhail_lambda", "dhail_vsedi_min",
            "dhail_vsedi_max"],
            "Ice":
            ["dice_a_geo", "dice_b_geo", "dice_min_x", "dice_max_x",
            "dice_sc_theta_q", "dice_sc_delta_q", "dice_sc_theta_n",
            "dice_sc_delta_n", "dice_s_vel", "dice_a_vel", "dice_b_vel",
            "dice_rho_v", "dice_c_z", "dice_sc_coll_n", "dice_cmu0",
            "dice_cmu1", "dice_cmu2", "dice_cmu3", "dice_cmu4", "dice_cmu5",
            "dice_alpha", "dice_beta", "dice_gamma", "dice_nu", "dice_g1",
            "dice_g2", "dice_mu", "dice_nm1", "dice_nm2", "dice_nm3",
            "dice_q_crit_c", "dice_d_crit_c", "dice_ecoll_c", "dice_cap",
            "dice_a_ven", "dice_b_ven", "dice_c_s", "dice_a_f", "dice_b_f",
            "dice_alfa_n", "dice_alfa_q", "dice_lambda", "dice_vsedi_min",
            "dice_vsedi_max"],
            "Snow":
            ["dsnow_a_geo", "dsnow_b_geo", "dsnow_min_x", "dsnow_max_x",
            "dsnow_sc_theta_q", "dsnow_sc_delta_q", "dsnow_sc_theta_n",
            "dsnow_sc_delta_n", "dsnow_s_vel", "dsnow_a_vel", "dsnow_b_vel",
            "dsnow_rho_v", "dsnow_c_z", "dsnow_sc_coll_n", "dsnow_cmu0",
            "dsnow_cmu1", "dsnow_cmu2", "dsnow_cmu3", "dsnow_cmu4",
            "dsnow_cmu5", "dsnow_alpha", "dsnow_beta", "dsnow_gamma",
            "dsnow_nu", "dsnow_g1", "dsnow_g2", "dsnow_mu", "dsnow_nm1",
            "dsnow_nm2", "dsnow_nm3", "dsnow_q_crit_c", "dsnow_d_crit_c",
            "dsnow_ecoll_c", "dsnow_cap", "dsnow_a_ven", "dsnow_b_ven",
            "dsnow_c_s", "dsnow_a_f", "dsnow_b_f", "dsnow_alfa_n",
            "dsnow_alfa_q", "dsnow_lambda", "dsnow_vsedi_min",
            "dsnow_vsedi_max"]}

direc_path = sys.argv[1]
key =  sys.argv[2]
out_params = [sys.argv[3]]

print("Plotting for {} with input ".format(key, sys.argv[3]))
in_params = in_params_dic[key]

columns = ["timestep", "Output Parameter", "trajectory"] + in_params + out_params
if "p" not in out_params:
    columns += ["p"]

t = timer()
data = Deriv_dask(
    direc=direc_path,
    parquet=True,
    columns=columns
)
t2 = timer()
print("Loading done in {}s".format(t2-t))

print("Got trajectories: {}".format(data.data["trajectory"].unique().compute()))
t = timer()
data.plot_two_ds(
    in_params=in_params,
    out_params=out_params,
    x_axis="p",
    mapped=False,
    trajectories=trajectories,
    scatter=False,
    frac=None,
    percentile=None,
    line_deriv=True,
    prefix=key,
    **kwargs
)
t2 = timer()
print("Plot two (p) done in {}s\n".format(t2-t))
