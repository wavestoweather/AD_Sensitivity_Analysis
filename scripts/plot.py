import numpy as np
import scipy.special as sp
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import sys
import os


prefix = sys.argv[1] + "_OUTPUT"
prefix_title = sys.argv[1]

# Load the references
refs = np.loadtxt('reference_values.txt')
Tref = refs[0]
pref = refs[1]/100 # We want hPa
qref = refs[2]
Nref = refs[3]
wref = refs[4]
tref = refs[5]

SMALL_SIZE = 8
MEDIUM_SIZE = 10
matplotlib.rcParams.update({'font.size': SMALL_SIZE})
matplotlib.rcParams.update({'axes.titlesize': SMALL_SIZE})
matplotlib.rcParams.update({'axes.labelsize': SMALL_SIZE})
matplotlib.rcParams.update({'xtick.labelsize': SMALL_SIZE})
matplotlib.rcParams.update({'ytick.labelsize': SMALL_SIZE})
matplotlib.rcParams.update({'legend.fontsize': SMALL_SIZE})
sns.set_style("darkgrid")

def plot_something(x, y, data, title, xlabel, ylabel, suffix="", one_traj=False,
        y_list1=None, y_list2=None):
    x_ticks = np.arange(0, data[x].max() + 19, 20)
    fig = plt.figure()
    ax = fig.add_subplot(221)
    ax.set_title(title.capitalize(), fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if one_traj:
        ax.plot(x, y, data=data)
    else:
        sns.lineplot(x, y, hue="trajectory", data=data, ax=ax)

    ax.get_yaxis().get_major_formatter().set_powerlimits((-1,1))
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.set_xticks(x_ticks)

    ax = fig.add_subplot(222)
    ax.set_title('Deriv. of {}'.format(title), fontsize=8)
    ax.set_xlabel(xlabel)
    if y_list1 is None:
        y_list1 = [("da_1", r"$\mathrm{d}" + y + r"/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}" + y + r"/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}" + y + r"/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}" + y + r"/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}" + y + r"/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}" + y + r"/\mathrm{d}N_c$")]

    for key, label in y_list1:
        if one_traj:
            ax.plot(x, key, data=data, label=label)
        else:
            sns.lineplot(x, key, hue="trajectory", data=data, label=label)

    ax.get_yaxis().get_major_formatter().set_powerlimits((-1,1))
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.legend(loc="best")
    ax.set_xticks(x_ticks)

    ax = fig.add_subplot(223)
    ax.set_title(prefix_title + ': deriv. of {}'.format(title), fontsize=8)
    ax.set_xlabel(xlabel)

    if y_list2 is None:
        y_list2 = [("dgamma", r"$\mathrm{d}" + y + r"/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}" + y + r"/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}" + y + r"/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}" + y + r"/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}" + y + r"/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}" + y + r"/\mathrm{d}\zeta$")]

    for key, label in y_list2:
        if one_traj:
            ax.plot(x, key, data=data, label=label)
        else:
            sns.lineplot(x, key, hue="trajectory", data=data, label=label)

    ax.get_yaxis().get_major_formatter().set_powerlimits((-1,1))
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.legend(loc="best")
    ax.set_xticks(x_ticks)
    plt.tight_layout()
    i = 0
    filename = ("pics/" + prefix_title + "_derivatives_"
                + title + suffix + "_" + "{:03d}".format(i) + ".png")
    while os.path.isfile(filename):
        i = i+1
        filename = ("pics/" + prefix_title + "_derivatives_"
                    + title + suffix + "_" + "{:03d}".format(i) + ".png")

    plt.savefig(filename, dpi=300)
    plt.close()


def plot_all(prefix="", suffix="", one_traj=False):
    # Load the single solution
    print("Plotting with suffix " + suffix)
    data = pd.read_csv("data/" + prefix + suffix + ".txt", sep=",")
    diffp = pd.read_csv("data/" + prefix + suffix + "_diff_0.txt", sep=",")
    diffp["pressure"] = pref*data["p"]
    diffp["timestep"] = data["timestep"]
    if one_traj:
        diffp = diffp[diffp["trajectory"] == 0]
        data = data[data["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}p/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}p/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}p/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}p/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}p/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}p/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}p/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}p/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}p/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}p/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}p/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}p/\mathrm{d}\zeta$")]
    plot_something("timestep", "pressure", diffp,
        "pressure", "Time (s)", "Pressure (Pa)", suffix,
        one_traj, y_list1, y_list2)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}p/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}p/\mathrm{d}b_r$"),
        ("drain_min_x", r"$\mathrm{d}p/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x", r"$\mathrm{d}p/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}p/\mathrm{d}\rho_v$")]
    y_list2 = [("drain_cmu0", r"$\mathrm{d}p/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}p/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}p/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}p/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}p/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}p/\mathrm{d}c_5$")]
    plot_something("timestep", "pressure", diffp,
        "pressure", "Time (s)", "Pressure (Pa)", suffix,
        one_traj, y_list1, y_list2)

    diffT = pd.read_csv("data/" + prefix + suffix + "_diff_1.txt", sep=",")
    diffT["timestep"] = data["timestep"]
    diffT["temperature"] = Tref*data["T"]
    if one_traj:
        diffT = diffT[diffT["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}T/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}T/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}T/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}T/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}T/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}T/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}T/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}T/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}T/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}T/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}T/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}T/\mathrm{d}\zeta$")]
    plot_something("timestep", "temperature", diffT,
        "temperature", "Time (s)", "Temperature (K)", suffix,
        one_traj, y_list1, y_list2)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}T/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}T/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}T/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}T/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}T/\mathrm{d}\rho_v$")]
    y_list2 = [("drain_cmu0", r"$\mathrm{d}T/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}T/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}T/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}T/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}T/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}T/\mathrm{d}c_5$")]
    plot_something("timestep", "temperature", diffT,
        "temperature", "Time (s)", "Temperature (K)", suffix,
        one_traj, y_list1, y_list2)


    diffw = pd.read_csv("data/" + prefix + suffix + "_diff_2.txt", sep=",")
    diffw["timestep"] = data["timestep"]
    diffw["uplift"] = wref*data["w"]
    if one_traj:
        diffw = diffw[diffw["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}w/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}w/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}w/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}w/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}w/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}w/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}w/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}w/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}w/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}w/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}w/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}w/\mathrm{d}\zeta$")]
    plot_something("timestep", "uplift", diffw,
        "Vertical velocity", "Time (s)", "W (m/s)", suffix,
        one_traj, y_list1, y_list2)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}w/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}w/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}w/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}w/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}w/\mathrm{d}\rho_v$")]
    y_list2 = [("drain_cmu0", r"$\mathrm{d}w/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}w/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}w/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}w/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}w/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}w/\mathrm{d}c_5$")]
    plot_something("timestep", "uplift", diffw,
        "Vertical velocity", "Time (s)", "W (m/s)", suffix,
        one_traj, y_list1, y_list2)

    diffS = pd.read_csv("data/" + prefix + suffix + "_diff_3.txt", sep=",")
    diffS["timestep"] = data["timestep"]
    diffS["saturation"] = data["S"]
    if one_traj:
        diffS = diffS[diffS["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}S/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}S/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}S/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}S/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}S/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}S/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}S/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}S/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}S/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}S/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}S/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}S/\mathrm{d}\zeta$")]
    plot_something("timestep", "saturation", diffS,
        "Saturation ratio", "Time (s)", "Saturation (1)", suffix,
        one_traj, y_list1, y_list2)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}S/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}S/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}S/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}S/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}S/\mathrm{d}\rho_v$")]
    y_list2 = [("drain_cmu0", r"$\mathrm{d}S/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}S/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}S/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}S/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}S/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}S/\mathrm{d}c_5$")]
    plot_something("timestep", "saturation", diffS,
        "Saturation ratio", "Time (s)", "Saturation (1)", suffix,
        one_traj, y_list1, y_list2)

    diffqc = pd.read_csv("data/" + prefix + suffix + "_diff_4.txt", sep=",")
    diffqc["timestep"] = data["timestep"]
    diffqc["qc"] = qref*data["qc"]
    if one_traj:
        diffqc = diffqc[diffqc["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}q_c/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}q_c/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}q_c/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}q_c/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}q_c/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}q_c/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}q_c/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}q_c/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}q_c/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}q_c/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}q_c/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}q_c/\mathrm{d}\zeta$")]
    plot_something("timestep", "qc", diffqc,
        "Cloud droplet mixing-ratio", "Time (s)", "qc (kg/kg)", suffix,
        one_traj, y_list1, y_list2)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}q_c/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}q_c/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}q_c/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}q_c/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}q_c/\mathrm{d}\rho_v$")]

    y_list2 = [("drain_cmu0", r"$\mathrm{d}q_c/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}q_c/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}q_c/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}q_c/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}q_c/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}q_c/\mathrm{d}c_5$")]
    plot_something("timestep", "qc", diffqc,
        "Cloud droplet mixing-ratio", "Time (s)", "qc (kg/kg)", suffix,
        one_traj, y_list1, y_list2)

    diffqr = pd.read_csv("data/" + prefix + suffix + "_diff_5.txt", sep=",")
    diffqr["timestep"] = data["timestep"]
    diffqr["qr"] = qref*data["qr"]
    if one_traj:
        diffqr = diffqr[diffqr["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}q_r/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}q_r/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}q_r/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}q_r/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}q_r/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}q_r/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}q_r/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}q_r/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}q_r/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}q_r/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}q_r/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}q_r/\mathrm{d}\zeta$")]
    plot_something("timestep", "qr", diffqr,
        "Rain drop mixing-ratio", "Time (s)", "qr (kg/kg)", suffix,
        one_traj, y_list1, y_list2)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}q_r/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}q_r/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}q_r/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}q_r/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}q_r/\mathrm{d}\rho_v$")]

    y_list2 = [("drain_cmu0", r"$\mathrm{d}q_r/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}q_r/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}q_r/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}q_r/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}q_r/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}q_r/\mathrm{d}c_5$")]
    plot_something("timestep", "qr", diffqr,
        "Rain drop mixing-ratio", "Time (s)", "qr (kg/kg)", suffix, one_traj,
        y_list1, y_list2)

    diffqv = pd.read_csv("data/" + prefix + suffix + "_diff_6.txt", sep=",")
    diffqv["timestep"] = data["timestep"]
    diffqv["qv"] = qref*data["qv"]
    if one_traj:
        diffqv = diffqv[diffqv["trajectory"] == 0]

    y_list1 = [("da_1", r"$\mathrm{d}q_v/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}q_v/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}q_v/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}q_v/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}q_v/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}q_v/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}q_v/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}q_v/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}q_v/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}q_v/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}q_v/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}q_v/\mathrm{d}\zeta$")]
    plot_something("timestep", "qv", diffqv,
        "Water vapor mixing-ratio", "Time (s)", "qv (kg/kg)", suffix,
        one_traj, y_list1, y_list2)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}q_v/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}q_v/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}q_v/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}q_v/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}q_v/\mathrm{d}\rho_v$")]

    y_list2 = [("drain_cmu0", r"$\mathrm{d}q_v/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}q_v/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}q_v/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}q_v/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}q_v/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}q_v/\mathrm{d}c_5$")]
    plot_something("timestep", "qv", diffqv,
        "Water vapor mixing-ratio", "Time (s)", "qv (kg/kg)", suffix, one_traj,
        y_list1, y_list2)

    diffqi = pd.read_csv("data/" + prefix + suffix + "_diff_10.txt", sep=",")
    diffqi["timestep"] = data["timestep"]
    diffqi["qi"] = qref*data["qi"]
    if one_traj:
        diffqi = diffqi[diffqi["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}q_i/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}q_i/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}q_i/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}q_i/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}q_i/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}q_i/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}q_i/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}q_i/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}q_i/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}q_i/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}q_i/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}q_i/\mathrm{d}\zeta$")]
    plot_something("timestep", "qi", diffqi,
        "Ice mixing-ratio", "Time (s)", "qi (kg/kg)", suffix, one_traj)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}q_i/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}q_i/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}q_i/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}q_i/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}q_i/\mathrm{d}\rho_v$")]

    y_list2 = [("drain_cmu0", r"$\mathrm{d}q_i/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}q_i/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}q_i/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}q_i/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}q_i/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}q_i/\mathrm{d}c_5$")]
    plot_something("timestep", "qi", diffqi,
        "Ice mixing-ratio", "Time (s)", "qi (kg/kg)", suffix, one_traj,
        y_list1, y_list2)
    ###
    diffqs = pd.read_csv("data/" + prefix + suffix + "_diff_15.txt", sep=",")
    diffqs["timestep"] = data["timestep"]
    diffqs["qs"] = qref*data["qs"]
    if one_traj:
        diffqs = diffqs[diffqs["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}q_s/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}q_s/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}q_s/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}q_s/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}q_s/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}q_s/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}q_s/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}q_s/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}q_s/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}q_s/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}q_s/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}q_s/\mathrm{d}\zeta$")]
    plot_something("timestep", "qs", diffqs,
        "Snow mixing-ratio", "Time (s)", "qs (kg/kg)", suffix, one_traj)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}q_s/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}q_s/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}q_s/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}q_s/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}q_s/\mathrm{d}\rho_v$")]

    y_list2 = [("drain_cmu0", r"$\mathrm{d}q_s/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}q_s/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}q_s/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}q_s/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}q_s/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}q_s/\mathrm{d}c_5$")]
    plot_something("timestep", "qs", diffqs,
        "Snow mixing-ratio", "Time (s)", "qs (kg/kg)", suffix, one_traj,
        y_list1, y_list2)
    ###
    diffqh = pd.read_csv("data/" + prefix + suffix + "_diff_17.txt", sep=",")
    diffqh["timestep"] = data["timestep"]
    diffqh["qh"] = qref*data["qh"]
    if one_traj:
        diffqh = diffqh[diffqh["trajectory"] == 0]
    y_list1 = [("da_1", r"$\mathrm{d}q_h/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}q_h/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}q_h/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}q_h/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}q_h/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}q_h/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}q_h/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}q_h/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}q_h/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}q_h/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}q_h/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}q_h/\mathrm{d}\zeta$")]
    plot_something("timestep", "qh", diffqh,
        "Hail mixing-ratio", "Time (s)", "qh (kg/kg)", suffix, one_traj)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}q_h/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}q_h/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}q_h/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}q_h/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}q_h/\mathrm{d}\rho_v$")]

    y_list2 = [("drain_cmu0", r"$\mathrm{d}q_h/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}q_h/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}q_h/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}q_h/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}q_h/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}q_h/\mathrm{d}c_5$")]
    plot_something("timestep", "qh", diffqh,
        "Hail mixing-ratio", "Time (s)", "qh (kg/kg)", suffix, one_traj,
        y_list1, y_list2)

    # Plot the sum of the mixing ratios
    diffsum = diffqh.copy()
    diffsum["qs"] = diffqs["qs"]
    diffsum["qc"] = diffqc["qc"]
    diffsum["qi"] = diffqi["qi"]
    diffsum["qv"] = diffqv["qv"]
    diffsum["qr"] = diffqr["qr"]
    # diffsum = diffqh.add(diffqs.add(diffqi.add(diffqv.add(diffqr.add(diffqc)))))
    diffsum["qsum"] = diffsum[["qh", "qs", "qi", "qv", "qr", "qc"]].sum(axis=1)

    y_list1 = [("da_1", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}a_1$"),
            ("da_2", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}a_2$"),
            ("de_1", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}e_1$"),
            ("de_2", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}e_2$"),
            ("dd", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}d$"),
            ("dN_c", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}N_c$")]
    y_list2 = [("dgamma", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\gamma$"),
            ("dbeta_c", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\beta_c$"),
            ("dbeta_r", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\beta_r$"),
            ("ddelta1", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\delta_1$"),
            ("ddelta2", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\delta_2$"),
            ("dzeta", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\zeta$")]
    plot_something("timestep", "qsum", diffsum,
        "All mixing-ratios", "Time (s)", "q_sum (kg/kg)", suffix, one_traj)

    y_list1 = [("drain_a_geo", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}a_r$"),
        ("drain_b_geo", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}b_r$"),
        ("drain_min_x",
         r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\overline{x}_{\mathrm{min}}$"),
        ("drain_max_x",
         r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\overline{x}_{\mathrm{max}}$"),
        ("drain_rho_v", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}\rho_v$")]

    y_list2 = [("drain_cmu0", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}c_0$"),
        ("drain_cmu1", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}c_1$"),
        ("drain_cmu2", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}c_2$"),
        ("drain_cmu3", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}c_3$"),
        ("drain_cmu4", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}c_4$"),
        ("drain_cmu5", r"$\mathrm{d}q_{\mathrm{sum}}/\mathrm{d}c_5$")]
    plot_something("timestep", "qsum", diffsum,
        "All mixing-ratios", "Time (s)", "q_sum (kg/kg)", suffix, one_traj,
        y_list1, y_list2)

# plot_all(prefix=prefix, one_traj=True)
try:
    plot_all(prefix=prefix, suffix="_start_over", one_traj=True)
except:
    print("Failed to plot start_over")
try:
    plot_all(prefix=prefix, suffix="_start_over_fixed", one_traj=True)
except:
    print("Failed to plot start_over_fixed")
# plot_all("_fixed", True)
# plot_all(prefix=prefix, suffix="_start_over_fixed", one_traj=True)
