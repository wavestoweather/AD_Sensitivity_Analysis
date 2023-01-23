#ifndef INCLUDE_MISC_CONSTANTS_REGRID_H_
#define INCLUDE_MISC_CONSTANTS_REGRID_H_

#include <array>
#include <string>
#include <vector>

const std::array<std::string, 38> order_sens = {
        "dkin_visc_air", "da_prime", "db_prime", "dc_prime", "da_v", "db_v", "drain_cmu3", "drain_nu",
        "db_ccn_1", "db_ccn_2", "db_ccn_3", "dc_ccn_1", "dc_ccn_2", "dc_ccn_3", "dg_ccn_1", "dg_ccn_2", "dh_ccn_1",
        "di_ccn_1", "dhande_ccn_fac",
        "drain_a_geo", "drain_b_geo", "dsnow_a_geo", "dsnow_b_geo", "dice_a_geo", "dice_b_geo", "dgraupel_a_geo",
        "dgraupel_b_geo",
        "dcloud_b_vel", "drain_a_vel", "drain_b_vel", "dsnow_b_vel", "dice_b_vel", "dgraupel_a_vel",
        "dgraupel_b_vel",
        "dT_mult_min", "dp_sat_melt", "da_HET", "drain_mu"
};
const std::array<std::string, 24> order_par = {
        "pressure",
        "T",
        "w",
        "S",
        "QC",
        "QR",
        "QV",
        "NCCLOUD",
        "NCRAIN",
        "QI",
        "NCICE",
        "QS",
        "NCSNOW",
        "QG",
        "NCGRAUPEL",
        "QH",
        "NCHAIL",
        "latent_heat",
        "latent_cool",
        "lat",
        "lon",
        "time_after_ascent",
        "asc600",
        "time"
};

const std::array<std::string, 2> order_misc = {
        "Top_Parameter",
        "counts"
};

enum Par_idx {
    pressure,
    temperature,
    ascent,
    saturation,
    qc,
    qr,
    qv,
    nc,
    nr,
    qi,
    ni,
    qs,
    ns,
    qg,
    ng,
    qh,
    nh,
    heat,
    cool,
    lat,
    lon,
    time_after_ascent,
    asc600,
    time_val,
    n_pars
};

enum Misc_idx {
    ranking,
    counts,
    n_misc
};

enum Sens_par_idx {
    dkin_visc_air,
    da_prime,
    db_prime,
    dc_prime,
    da_v,
    db_v,
    drain_cmu3,
    drain_nu,
    db_ccn_1,
    db_ccn_2,
    db_ccn_3,
    dc_ccn_1,
    dc_ccn_2,
    dc_ccn_3,
    dg_ccn_1,
    dg_ccn_2,
    dh_ccn_1,
    di_ccn_1,
    dhande_ccn_fac,
    drain_a_geo,
    drain_b_geo,
    dsnow_a_geo,
    dsnow_b_geo,
    dice_a_geo,
    dice_b_geo,
    dgraupel_a_geo,
    dgraupel_b_geo,
    dcloud_b_vel,
    drain_a_vel,
    drain_b_vel,
    dsnow_b_vel,
    dice_b_vel,
    dgraupel_a_vel,
    dgraupel_b_vel,
    dT_mult_min,
    dp_sat_melt,
    da_HET,
    drain_mu,
    n_sens_pars
};
enum Dim_idx {
    time_dim_idx,
    trajectory_dim_idx,
    ensemble_dim_idx,
    output_para_idx,
    n_dims
};

enum Grid_dim_idx {
    outp_dim_idx,
    time_griddim_idx,
    p_dim_idx,
    lat_dim_idx,
    lon_dim_idx,
    n_griddims
};


#endif  // INCLUDE_MISC_CONSTANTS_REGRID_H_
