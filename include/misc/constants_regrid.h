#ifndef INCLUDE_MISC_CONSTANTS_REGRID_H_
#define INCLUDE_MISC_CONSTANTS_REGRID_H_

#include <array>
#include <string>
#include <vector>

#define DELTA_TIMESTEP 30
/**
 * Define the compression level. Level 9 is the best compression but level 6
 * is much faster and uses something between 2 and 5 times as much space.
 */
#define COMPRESSION_LEVEL 6

const std::array<std::string, 29> order_sens = {
        "db_ccn_1", "db_ccn_2", "db_ccn_3", "db_v",
        "dc_ccn_1", "dc_ccn_2", "dc_ccn_3", "dg_ccn_1",
        "dgraupel_a_geo", "dgraupel_b_geo", "dgraupel_b_vel", "dh_ccn_1",
        "dhande_ccn_fac", "di_ccn_1", "dice_a_geo", "dice_b_geo",
        "dice_b_vel", "dkin_visc_air", "dp_sat_melt", "drain_a_geo",
        "drain_a_vel", "drain_alpha", "drain_b_geo", "drain_b_vel",
        "drain_beta", "drain_cmu3", "drain_mu", "drain_nu",
        "dsnow_b_geo"
};
const std::array<std::string, 26> order_par = {
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
        "relative_lat",
        "relative_lon",
        "time_after_ascent",
        "asc600",
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
    rel_lat,
    rel_lon,
    time_after_ascent,
    asc600,
    n_pars
};

enum Misc_idx {
    ranking,
    counts,
    n_misc
};

enum Sens_par_idx {
    db_ccn_1,
    db_ccn_2,
    db_ccn_3,
    db_v,
    dc_ccn_1,
    dc_ccn_2,
    dc_ccn_3,
    dg_ccn_1,
    dgraupel_a_geo,
    dgraupel_b_geo,
    dgraupel_b_vel,
    dh_ccn_1,
    dhande_ccn_fac,
    di_ccn_1,
    dice_a_geo,
    dice_b_geo,
    dice_b_vel,
    dkin_visc_air,
    dp_sat_melt,
    drain_a_geo,
    drain_a_vel,
    drain_alpha,
    drain_b_geo,
    drain_b_vel,
    drain_beta,
    drain_cmu3,
    drain_mu,
    drain_nu,
    dsnow_b_geo,
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
