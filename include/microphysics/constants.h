#pragma once

#include <cmath>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "codi.hpp"

/** @defgroup constants Constants
 * Various constants for accessing data in the right order and model constants
 * for which most are rather uninteresting (i.e. no AD available).
 * @{
 */

typedef bool(*track_func)(const int&, const bool&);

////////////////////////////////////////////////////////////////////////////////
// Indices of the output parameters
////////////////////////////////////////////////////////////////////////////////
#define p_idx 0             /*!< Pressure index */
#define T_idx 1             /*!< Temperature index */
#define w_idx 2             /*!< Vertical acceleration index */
#define S_idx 3             /*!< Satruation index */
#define qc_idx 4            /*!< Cloud droplet mixing ratio index */
#define qr_idx 5            /*!< Rain droplet mixing ratio index */
#define qv_idx 6            /*!< Water vapor mixing ratio index */
#define Nc_idx 7            /*!< Number of cloud droplets index */
#define Nr_idx 8            /*!< Number of rain droplets index */
#define qi_idx 9           /*!< Ice mixing ratio index */
#define Ni_idx 10           /*!< Number of ice crystals index */
#define qs_idx 11           /*!< Snow mixing ratio index */
#define Ns_idx 12           /*!< Number of snow particles index */
#define qg_idx 13           /*!< Graupel mixing ratio index */
#define Ng_idx 14           /*!< Number of graupel particles index */
#define qh_idx 15           /*!< Hail mixing ratio index */
#define Nh_idx 16           /*!< Number of hail particles index */
#define qi_out_idx 17       /*!< Ice mixing ratio precipitation index */
#define qs_out_idx 18       /*!< Snow mixing ratio precipitation index */
#define qr_out_idx 19       /*!< Rain mixing ratio precipitation index */
#define qg_out_idx 20       /*!< Graupel mixing ratio precipitation index */
#define qh_out_idx 21       /*!< Hail mixing ratio precipitation index */
#define lat_heat_idx 22     /*!< Latent heating index */
#define lat_cool_idx 23     /*!< Latent cooling index */
#define Ni_out_idx 24       /*!< Ice particles precipitation index */
#define Ns_out_idx 25       /*!< Snow particles ratio precipitation index */
#define Nr_out_idx 26       /*!< Rain droplets ratio precipitation index */
#define Ng_out_idx 27       /*!< Graupel particles ratio precipitation index */
#define Nh_out_idx 28       /*!< Hail particles ratio precipitation index */
#define z_idx 29            /*!< Altitude */
#define n_inact_idx 30      /*!< Number of inactive nuclei (ie due to being activated before) */
#define depo_idx 31         /*!< Number of deposited nuclei */
#define sub_idx 32          /*!< Sublimination number */

#if defined(RK4_ONE_MOMENT)
#define num_comp 9          /*!< Number of output elements of a model */
#define num_par 12          /*!< Number of gradients */

#elif defined(RK4ICE) || defined(RK4NOICE)
#define num_comp 33         /*!< Number of output elements of a model */
#if defined(B_EIGHT)
#define num_par (56*6+124+18)
#else
#define num_par (56*6+134+18)  /*!< Number of gradients. 56 for each particle + model constants + initial conditions */
#endif

#endif

#define num_par_init 18  /*!< Number of gradients w.r.t. initial conditions. */

// Those are for an inflow vector
#define qi_in_idx 0         /*!< Ice input index for another vector */
#define qs_in_idx 1         /*!< Snow input index for another vector */
#define qr_in_idx 2         /*!< Rain input index for another vector */
#define qg_in_idx 3         /*!< Graupel input index for another vector */
#define Ni_in_idx 4         /*!< Ice input index for another vector */
#define Ns_in_idx 5         /*!< Snow input index for another vector */
#define Nr_in_idx 6         /*!< Rain input index for another vector */
#define Ng_in_idx 7         /*!< Graupel input index for another vector */
#if defined B_EIGHT && !defined(TURBULENCE)
#define qh_in_idx 8         /*!< Hail input index for another vector */
#define Nh_in_idx 9         /*!< Hail input index for another vector */
#define num_inflows 10       /*!< Number of parameters for inflowing stuff */
#elif defined(MET3D) && defined(TURBULENCE) && !defined(B_EIGHT)
#define qv_in_idx 8         /*!< Vapor input index for another vector */
#define num_inflows 9       /*!< Number of parameters for inflowing stuff */
#elif defined(MET3D) && defined(TURBULENCE) && defined(B_EIGHT)
#define qh_in_idx 8         /*!< Hail input index for another vector */
#define Nh_in_idx 9         /*!< Hail input index for another vector */
#define qv_in_idx 10         /*!< Vapor input index for another vector */
#define num_inflows 11       /*!< Number of parameters for inflowing stuff */
#else
#define num_inflows 8       /*!< Number of parameters for inflowing stuff */
#endif

#define CHECKPOINT_MESSAGE 1
#define SIMULATION_MESSAGE 2
#define SIMULATION_MESSAGE_FLAGS 3
#define SIMULATION_MESSAGE_STR 4
#define SIMULATION_MESSAGE_INT 5
#define SIMULATION_MESSAGE_NSNAPS 6

////////////////////////////////////////////////////////////////////////////////
// Simulation modes
////////////////////////////////////////////////////////////////////////////////
/**
 * Input data are trajectories. Perturbed ensembles are possible.
 * The output contains sensitivities.
 */
#define trajectory_sensitvity_perturbance 0
/**
 * Input data are trajectories. Perturbed ensembles are not possible.
 * The output contains sensitivities.
 */
#define trajectory_sensitivity 1
/**
 * Input data are trajectories. Perturbed ensembles are possible.
 * The output does not contain sensitivities.
 */
#define trajectory_perturbance 2
/**
 * Input data is a grid. Perturbed ensembles are not possible.
 * The output contains sensitivities.
 */
#define grid_sensitivity 3
/**
 * Input data is a single trajectory. Perturbed ensembles are started according
 * to the ensemble configuration file.
 * The output contains sensitivities.
 */
#define limited_time_ensembles 4


////////////////////////////////////////////////////////////////////////////////
// Random generators
////////////////////////////////////////////////////////////////////////////////

extern std::random_device rand_device;
extern std::mt19937 rand_generator;

////////////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////////////
extern bool trace;
#if defined(TRACE_TIME)
// Relative to ascent time
extern double trace_time;
const double trace_start = 0;
const double trace_end = 180;
#endif

#if defined(RK4_ONE_MOMENT)
/**
 * Used for header files of output parameters.
 */
const std::vector<std::string> output_par_idx = {
    "p", "T", "w", "S", "qc", "qr", "qv", "Nc", "Nr"};

/**
 * Used for header files of gradients.
 */
const std::vector<std::string> output_grad_idx = {
    "da_1", "da_2", "de_1", "de_2", "dd", "dN_c", "dgamma", "dbeta_c",
    "dbeta_r", "ddelta1", "ddelta2", "dzeta"};
#elif defined(RK4ICE) || defined(RK4NOICE)
#ifdef MET3D
const std::vector<std::string> output_par_idx = {
    "pressure", "T", "w", "S", "QC", "QR", "QV", "NCCLOUD", "NCRAIN",
    "QI", "NCICE", "QS", "NCSNOW", "QG", "NCGRAUPEL", "QH", "NCHAIL",
    "QI_OUT", "QS_OUT", "QR_OUT", "QG_OUT", "QH_OUT",
    "latent_heat", "latent_cool", "NI_OUT", "NS_OUT", "NR_OUT",
    "NG_OUT", "NH_OUT", "z", "Inactive", "deposition", "sublimination"};
#else
/**
 * Used for header files of output parameters.
 */
const std::vector<std::string> output_par_idx = {
    "p", "T", "w", "S", "qc", "qr", "qv", "Nc", "Nr",
    "qi", "Ni", "qs", "Ns", "qg", "Ng", "qh", "Nh",
    "qiout", "qsout", "qrout", "qgout", "qhout",
    "latent_heat", "latent_cool", "Niout", "Nsout", "Nrout",
    "Ngout", "Nhout", "z", "Inactive", "deposition", "sublimination"};
#endif
/**
 * Used for header files of gradients.
 */
const std::vector<std::string> output_grad_idx = {
    // model parameters
    "da_1", "da_2", "de_1", "de_2", "dd", "dN_c", "dgamma", "dbeta_c",
    "dbeta_r", "ddelta1", "ddelta2", "dzeta",
    "drain_gfak", "dcloud_k_au",
    "dcloud_k_sc", "dkc_autocon", "dinv_z", "dw",
    "dq_crit_i", "dD_crit_i", "dD_conv_i", "dq_crit_r",
    "dD_crit_r", "dq_crit_fr", "dD_coll_c", "dq_crit",
    "dD_conv_sg", "dD_conv_ig", "dx_conv", "dparcel_height",
    "dalpha_spacefilling", "dT_nuc", "dT_freeze", "dT_f",
    "dD_eq", "drho_w", "drho_0", "drho_vel",
    "drho_vel_c", "drho_ice", "dM_w", "dM_a",
    "dR_universal", "dEpsilon", "dgravity_acc", "dR_a",
    "dR_v", "da_v", "db_v", "da_prime",
    "db_prime", "dc_prime", "dK_T", "dL_wd",
    "dL_ed", "dD_v", "decoll_min", "decoll_gg",
    "decoll_gg_wet", "dkin_visc_air", "dC_mult", "dT_mult_min",
    "dT_mult_max", "dT_mult_opt", "dconst0", "dconst3",
    "dconst4", "dconst5", "dD_rainfrz_gh", "dD_rainfrz_ig", "ddv0",
    "dp_sat_melt", "dcp", "dk_b", "da_HET",
    "db_HET", "dN_sc", "dn_f", "dN_avo",
    "dna_dust", "dna_soot", "dna_orga", "dni_het_max",
    "dni_hom_max", "da_dep", "db_dep", "dc_dep",
    "dd_dep", "dnim_imm", "dnin_dep", "dalf_imm",
    "dbet_dep", "dbet_imm", "dr_const", "dr1_const",
    "dcv", "dp_sat_const_a", "dp_sat_ice_const_a", "dp_sat_const_b",
    "dp_sat_ice_const_b", "dp_sat_low_temp", "dT_sat_low_temp", "dalpha_depo",
    "dr_0", "dk_1_conv", "dk_2_conv", "dk_1_accr",
    "dk_r", "da_ccn_1", "da_ccn_2", "da_ccn_3", "da_ccn_4",
    "db_ccn_1", "db_ccn_2", "db_ccn_3", "db_ccn_4",
    "dc_ccn_1", "dc_ccn_2", "dc_ccn_3", "dc_ccn_4",
    "dd_ccn_1", "dd_ccn_2", "dd_ccn_3", "dd_ccn_4",
#if defined(B_EIGHT)
    "dp_ccn", "dh_ccn_1", "dh_ccn_2", "dh_ccn_3",
    "dg_ccn_1", "dg_ccn_2", "dg_ccn_3",
    "di_ccn_1", "di_ccn_2", "dhande_ccn_fac",
#endif
    // Rain
    "drain_a_geo", "drain_b_geo", "drain_min_x", "drain_min_x_act",
    "drain_min_x_nuc_homo", "drain_min_x_nuc_hetero", "drain_min_x_melt",
    "drain_min_x_evap", "drain_min_x_freezing", "drain_min_x_depo",
    "drain_min_x_collision", "drain_min_x_collection",
    "drain_min_x_conversion", "drain_min_x_sedimentation",
    "drain_min_x_riming", "drain_max_x", "drain_sc_theta_q",
    "drain_sc_delta_q", "drain_sc_theta_n", "drain_sc_delta_n",
    "drain_s_vel", "drain_a_vel", "drain_b_vel", "drain_rho_v", "drain_c_z",
    "drain_sc_coll_n", "drain_cmu0", "drain_cmu1", "drain_cmu2", "drain_cmu3",
    "drain_cmu4", "drain_cmu5", "drain_alpha", "drain_beta", "drain_gamma",
    "drain_nu", "drain_g1", "drain_g2", "drain_mu", "drain_nm1", "drain_nm2",
    "drain_nm3", "drain_q_crit_c", "drain_d_crit_c", "drain_ecoll_c",
    "drain_cap", "drain_a_ven", "drain_b_ven", "drain_c_s", "drain_a_f",
    "drain_b_f", "drain_alfa_n", "drain_alfa_q", "drain_lambda",
    "drain_vsedi_min", "drain_vsedi_max",
    // Cloud
    "dcloud_a_geo", "dcloud_b_geo", "dcloud_min_x", "dcloud_min_x_act",
    "dcloud_min_x_nuc_homo", "dcloud_min_x_nuc_hetero", "dcloud_min_x_melt",
    "dcloud_min_x_evap", "dcloud_min_x_freezing", "dcloud_min_x_depo",
    "dcloud_min_x_collision", "dcloud_min_x_collection",
    "dcloud_min_x_conversion", "dcloud_min_x_sedimentation",
    "dcloud_min_x_riming", "dcloud_max_x", "dcloud_sc_theta_q",
    "dcloud_sc_delta_q", "dcloud_sc_theta_n", "dcloud_sc_delta_n",
    "dcloud_s_vel", "dcloud_a_vel", "dcloud_b_vel", "dcloud_rho_v",
    "dcloud_c_z", "dcloud_sc_coll_n", "dcloud_cmu0", "dcloud_cmu1",
    "dcloud_cmu2", "dcloud_cmu3", "dcloud_cmu4", "dcloud_cmu5",
    "dcloud_alpha", "dcloud_beta", "dcloud_gamma", "dcloud_nu", "dcloud_g1",
    "dcloud_g2", "dcloud_mu", "dcloud_nm1", "dcloud_nm2", "dcloud_nm3",
    "dcloud_q_crit_c", "dcloud_d_crit_c", "dcloud_ecoll_c", "dcloud_cap",
    "dcloud_a_ven", "dcloud_b_ven", "dcloud_c_s", "dcloud_a_f", "dcloud_b_f",
    "dcloud_alfa_n", "dcloud_alfa_q", "dcloud_lambda", "dcloud_vsedi_min",
    "dcloud_vsedi_max",
    // Graupel
    "dgraupel_a_geo", "dgraupel_b_geo", "dgraupel_min_x",
    "dgraupel_min_x_act", "dgraupel_min_x_nuc_homo",
    "dgraupel_min_x_nuc_hetero", "dgraupel_min_x_melt", "dgraupel_min_x_evap",
    "dgraupel_min_x_freezing", "dgraupel_min_x_depo",
    "dgraupel_min_x_collision", "dgraupel_min_x_collection",
    "dgraupel_min_x_conversion", "dgraupel_min_x_sedimentation",
    "dgraupel_min_x_riming", "dgraupel_max_x", "dgraupel_sc_theta_q",
    "dgraupel_sc_delta_q", "dgraupel_sc_theta_n", "dgraupel_sc_delta_n",
    "dgraupel_s_vel", "dgraupel_a_vel", "dgraupel_b_vel", "dgraupel_rho_v",
    "dgraupel_c_z", "dgraupel_sc_coll_n", "dgraupel_cmu0", "dgraupel_cmu1",
    "dgraupel_cmu2", "dgraupel_cmu3", "dgraupel_cmu4", "dgraupel_cmu5",
    "dgraupel_alpha", "dgraupel_beta", "dgraupel_gamma", "dgraupel_nu",
    "dgraupel_g1", "dgraupel_g2", "dgraupel_mu", "dgraupel_nm1",
    "dgraupel_nm2", "dgraupel_nm3", "dgraupel_q_crit_c", "dgraupel_d_crit_c",
    "dgraupel_ecoll_c", "dgraupel_cap", "dgraupel_a_ven", "dgraupel_b_ven",
    "dgraupel_c_s", "dgraupel_a_f", "dgraupel_b_f", "dgraupel_alfa_n",
    "dgraupel_alfa_q", "dgraupel_lambda", "dgraupel_vsedi_min",
    "dgraupel_vsedi_max",
    // Hail
    "dhail_a_geo", "dhail_b_geo", "dhail_min_x", "dhail_min_x_act",
    "dhail_min_x_nuc_homo", "dhail_min_x_nuc_hetero", "dhail_min_x_melt",
    "dhail_min_x_evap", "dhail_min_x_freezing", "dhail_min_x_depo",
    "dhail_min_x_collision", "dhail_min_x_collection",
    "dhail_min_x_conversion", "dhail_min_x_sedimentation",
    "dhail_min_x_riming", "dhail_max_x", "dhail_sc_theta_q",
    "dhail_sc_delta_q", "dhail_sc_theta_n", "dhail_sc_delta_n", "dhail_s_vel",
    "dhail_a_vel", "dhail_b_vel", "dhail_rho_v", "dhail_c_z",
    "dhail_sc_coll_n", "dhail_cmu0", "dhail_cmu1", "dhail_cmu2", "dhail_cmu3",
    "dhail_cmu4", "dhail_cmu5", "dhail_alpha", "dhail_beta", "dhail_gamma",
    "dhail_nu", "dhail_g1", "dhail_g2", "dhail_mu", "dhail_nm1", "dhail_nm2",
    "dhail_nm3", "dhail_q_crit_c", "dhail_d_crit_c", "dhail_ecoll_c",
    "dhail_cap", "dhail_a_ven", "dhail_b_ven", "dhail_c_s", "dhail_a_f",
    "dhail_b_f", "dhail_alfa_n", "dhail_alfa_q", "dhail_lambda",
    "dhail_vsedi_min", "dhail_vsedi_max",
    // Ice
    "dice_a_geo", "dice_b_geo", "dice_min_x", "dice_min_x_act",
    "dice_min_x_nuc_homo", "dice_min_x_nuc_hetero", "dice_min_x_melt",
    "dice_min_x_evap", "dice_min_x_freezing", "dice_min_x_depo",
    "dice_min_x_collision", "dice_min_x_collection", "dice_min_x_conversion",
    "dice_min_x_sedimentation", "dice_min_x_riming", "dice_max_x",
    "dice_sc_theta_q", "dice_sc_delta_q", "dice_sc_theta_n",
    "dice_sc_delta_n", "dice_s_vel", "dice_a_vel", "dice_b_vel", "dice_rho_v",
    "dice_c_z", "dice_sc_coll_n", "dice_cmu0", "dice_cmu1", "dice_cmu2",
    "dice_cmu3", "dice_cmu4", "dice_cmu5", "dice_alpha", "dice_beta",
    "dice_gamma", "dice_nu", "dice_g1", "dice_g2", "dice_mu", "dice_nm1",
    "dice_nm2", "dice_nm3", "dice_q_crit_c", "dice_d_crit_c", "dice_ecoll_c",
    "dice_cap", "dice_a_ven", "dice_b_ven", "dice_c_s", "dice_a_f",
    "dice_b_f", "dice_alfa_n", "dice_alfa_q", "dice_lambda", "dice_vsedi_min",
    "dice_vsedi_max",
    // Snow
    "dsnow_a_geo", "dsnow_b_geo", "dsnow_min_x", "dsnow_min_x_act",
    "dsnow_min_x_nuc_homo", "dsnow_min_x_nuc_hetero", "dsnow_min_x_melt",
    "dsnow_min_x_evap", "dsnow_min_x_freezing", "dsnow_min_x_depo",
    "dsnow_min_x_collision", "dsnow_min_x_collection",
    "dsnow_min_x_conversion", "dsnow_min_x_sedimentation",
    "dsnow_min_x_riming", "dsnow_max_x", "dsnow_sc_theta_q",
    "dsnow_sc_delta_q", "dsnow_sc_theta_n", "dsnow_sc_delta_n", "dsnow_s_vel",
    "dsnow_a_vel", "dsnow_b_vel", "dsnow_rho_v", "dsnow_c_z",
    "dsnow_sc_coll_n", "dsnow_cmu0", "dsnow_cmu1", "dsnow_cmu2", "dsnow_cmu3",
    "dsnow_cmu4", "dsnow_cmu5", "dsnow_alpha", "dsnow_beta", "dsnow_gamma",
    "dsnow_nu", "dsnow_g1", "dsnow_g2", "dsnow_mu", "dsnow_nm1", "dsnow_nm2",
    "dsnow_nm3", "dsnow_q_crit_c", "dsnow_d_crit_c", "dsnow_ecoll_c",
    "dsnow_cap", "dsnow_a_ven", "dsnow_b_ven", "dsnow_c_s", "dsnow_a_f",
    "dsnow_b_f", "dsnow_alfa_n", "dsnow_alfa_q", "dsnow_lambda",
    "dsnow_vsedi_min", "dsnow_vsedi_max",
    // initial conditions
    "dpressure_init", "dT_init", "dw_init", "dS_init", "dQC_init", "dQR_init",
    "dQV_init", "dNCCLOUD_init",  "dNCRAIN_init", "dQI_init", "dNCICE_init",
    "dQS_init", "dNCSNOW_init", "dQG_init", "dNCGRAUPEL_init", "dQH_init",
    "dNCHAIL_init", "dz_init"};
#endif

const std::vector<std::string> init_grad_idx = {
    "dpressure_init", "dT_init", "dw_init", "dS_init", "dQC_init",
    "dQR_init", "dQV_init", "dNCCLOUD_init", "dNCRAIN_init", "dQI_init",
    "dNCICE_init", "dQS_init", "dNCSNOW_init", "dQG_init", "dNCGRAUPEL_init",
    "dQH_init", "dNCHAIL_init", "dz_init"};

enum class Init_cons_idx: uint32_t{
    pressure,
    T,
    w,
    S,
    qc,
    qr,
    qv,
    Nc,
    Nr,
    qi,
    Ni,
    qs,
    Ns,
    qg,
    Ng,
    qh,
    Nh,
    z,
    n_items
};

/**
 * Mapping of json configuration names to parameters in a model_constants_t.
 * Only for initial conditions.
 */
std::unordered_map<std::string, Init_cons_idx> const table_init_param = {
    {"pressure_init", Init_cons_idx::pressure}, {"T_init", Init_cons_idx::T},
    {"w_init", Init_cons_idx::w}, {"S_init", Init_cons_idx::S},
    {"QC_init", Init_cons_idx::qc}, {"QR_init", Init_cons_idx::qr},
    {"QV_init", Init_cons_idx::qv}, {"NCCLOUD_init", Init_cons_idx::Nc},
    {"NCRAIN_init", Init_cons_idx::Nr}, {"QI_init", Init_cons_idx::qi},
    {"NCICE_init", Init_cons_idx::Ni}, {"QS_init", Init_cons_idx::qs},
    {"NCSNOW_init", Init_cons_idx::Ns}, {"QG_init", Init_cons_idx::qg},
    {"NCGRAUPEL_init", Init_cons_idx::Ng}, {"QH_init", Init_cons_idx::qh},
    {"NCHAIL_init", Init_cons_idx::Nh}, {"z_init", Init_cons_idx::z}
};

enum class Cons_idx: uint32_t{
    a1_prime,   /*!< Dimensional coefficient used in one-moment warm physics for qc and qr calculation */
    a2_prime,   /*!< Dimensional coefficient used in one-moment warm physics for qc and qr calculation */
    e1_prime,   /*!< Dimensional coefficients used in one-moment warm physics for temperature calculation */
    e2_prime,   /*!< Dimensional coefficients used in one-moment warm physics for temperature calculation */
    d_prime,    /*!< Dimensional coefficient used in one-moment warm physics qr calculation for sedimentation*/
    Nc_prime,   /*!< Number concentration of cloud droplets needed for one-moment scheme */

    gamma,      /*!< Exponent used in one-moment warm physics for qc and qr calculation */
    betac,      /*!< Exponent used in one-moment warm physics for qc and qr calculation */
    betar,      /*!< Exponent used in one-moment warm physics for qc and qr calculation */
    delta1,     /*!< Exponent used in one-moment warm physics for qv, qr, saturation and temperature calculation */
    delta2,     /*!< Exponent used in one-moment warm physics for qv, qr, saturation and temperature calculation */
    zeta,       /*!< Exponents used in one-moment warm physics for qr calculation */
#if defined(RK4ICE) || defined(RK4NOICE)
    rain_gfak,  /*!< Coefficient for gamma evaluation in rain evaporation */
    cloud_k_au, /*!< Coefficient for autoconversion of cloud to rain */
    cloud_k_sc, /*!< Coefficient for autoconversion of cloud to rain */
    kc_autocon, /*!< Kernel for autoconversion */
    /**
     * Inverse layer thickness. Used for sedimentation.
     * In Miltenberger (2016) the trajectories start every \f$100 \text{m}\f$
     * between the surface and \f$4 \text{km}\f$ altitude using COSMO-2, which
     * uses a mean spacing of \f$388 \text{m}\f$
     * with \f$13 \text{m}\f$ close to the surface and \f$1190 \text{m}\f$
     * at \f$23 \text{km}\f$.
     */
    inv_z,

    dw, /*!< Change in buoancy */
    q_crit_i,
    D_crit_i,
    D_conv_i,
    q_crit_r,
    D_crit_r,
    q_crit_fr,
    D_coll_c,
    q_crit,
    D_conv_sg,
    D_conv_ig,
    x_conv,
    parcel_height,
    alpha_spacefilling,
    T_nuc,
    T_freeze,
    T_f,
    D_eq,
    rho_w,
    rho_0,
    rho_vel,
    rho_vel_c,
    rho_ice,
    M_w,
    M_a,
    R_universal,
    Epsilon,
    gravity_acc,
    R_a,
    R_v,
    a_v,
    b_v,
    a_prime,
    b_prime,
    c_prime,
    K_T,
    L_wd,
    L_ed,
    D_v,
    ecoll_min,
    ecoll_gg,
    ecoll_gg_wet,
    kin_visc_air,
    C_mult,
    T_mult_min,
    T_mult_max,
    T_mult_opt,

    /**
     * Constant used in cloud riming.
     */
    const0,
    /**
     * Hallet-Mossop ice multiplication.
     * Constant used in ice - x and snow - x riming.
     */
    const3,
    /**
     * Hallet-Mossop ice multiplication.
     * Constant used in ice - x and snow - x riming.
     */
    const4,
    /**
     * Constant for conversions ice -> graupel, snow -> graupel,
     * melting (used in riming).
     */
    const5,
    D_rainfrz_ig,
    D_rainfrz_gh,
    dv0,
    p_sat_melt,
    cp,
    k_b,
    a_HET,
    b_HET,
    N_sc,
    n_f,
    N_avo,
    na_dust,
    na_soot,
    na_orga,
    ni_het_max,
    ni_hom_max,
    a_dep,
    b_dep,
    c_dep,
    d_dep,
    nim_imm,
    nin_dep,
    alf_imm,
    bet_dep,
    bet_imm,
    r_const,
    r1_const,
    cv,
    p_sat_const_a,
    p_sat_ice_const_a,
    p_sat_const_b,
    p_sat_ice_const_b,
    p_sat_low_temp,
    T_sat_low_temp,
    alpha_depo,
    r_0,
    k_1_conv,
    k_2_conv,
    k_1_accr,
    k_r,
    a_ccn_1, a_ccn_2, a_ccn_3, a_ccn_4,
    b_ccn_1, b_ccn_2, b_ccn_3, b_ccn_4,
    c_ccn_1, c_ccn_2, c_ccn_3, c_ccn_4,
    d_ccn_1, d_ccn_2, d_ccn_3, d_ccn_4,
#endif
#if defined(B_EIGHT)
    p_ccn,
    h_ccn_1, h_ccn_2, h_ccn_3,
    g_ccn_1, g_ccn_2, g_ccn_3,
    i_ccn_1, i_ccn_2,
    hande_ccn_fac,
#endif
    D_br_threshold, k_br, D_br, c_br,
    n_items
};

/**
 * Mapping of json configuration names to parameters in a model_constants_t
 */
std::unordered_map<std::string, Cons_idx> const table_param = {
    {"a_1", Cons_idx::a1_prime}, {"a_2", Cons_idx::a2_prime},
    {"e_1", Cons_idx::e1_prime}, {"e_2", Cons_idx::e2_prime},
    {"d", Cons_idx::d_prime}, {"N_c", Cons_idx::Nc_prime},
    {"gamma", Cons_idx::gamma}, {"beta_c", Cons_idx::betac}, {"beta_r", Cons_idx::betar},
    {"delta1", Cons_idx::delta1}, {"delta2", Cons_idx::delta2}, {"zeta", Cons_idx::zeta},
#if defined(RK4ICE) || defined(RK4NOICE)
    {"rain_gfak", Cons_idx::rain_gfak}, {"cloud_k_au", Cons_idx::cloud_k_au},
    {"cloud_k_sc", Cons_idx::cloud_k_sc}, {"kc_autocon", Cons_idx::kc_autocon},
    {"inv_z", Cons_idx::inv_z}, {"dw", Cons_idx::dw}, {"q_crit_i", Cons_idx::q_crit_i},
    {"D_crit_i", Cons_idx::D_crit_i}, {"D_conv_i", Cons_idx::D_conv_i},
    {"q_crit_r", Cons_idx::q_crit_r}, {"D_crit_r", Cons_idx::D_crit_r},
    {"q_crit_fr", Cons_idx::q_crit_fr}, {"D_coll_c", Cons_idx::D_coll_c},
    {"q_crit", Cons_idx::q_crit}, {"D_conv_sg", Cons_idx::D_conv_sg},
    {"D_conv_ig", Cons_idx::D_conv_ig}, {"x_conv", Cons_idx::x_conv},
    {"parcel_height", Cons_idx::parcel_height}, {"alpha_spacefilling", Cons_idx::alpha_spacefilling},
    {"T_nuc", Cons_idx::T_nuc}, {"T_freeze", Cons_idx::T_freeze}, {"T_f", Cons_idx::T_f},
    {"D_eq", Cons_idx::D_eq}, {"rho_w", Cons_idx::rho_w}, {"rho_0", Cons_idx::rho_0},
    {"rho_vel", Cons_idx::rho_vel}, {"rho_vel_c", Cons_idx::rho_vel_c}, {"rho_ice", Cons_idx::rho_ice},
    {"M_w", Cons_idx::M_w}, {"M_a", Cons_idx::M_a}, {"R_universal", Cons_idx::R_universal},
    {"Epsilon", Cons_idx::Epsilon}, {"gravity_acc", Cons_idx::gravity_acc},
    {"R_a", Cons_idx::R_a}, {"R_v", Cons_idx::R_v}, {"a_v", Cons_idx::a_v},
    {"b_v", Cons_idx::b_v}, {"a_prime", Cons_idx::a_prime}, {"b_prime", Cons_idx::b_prime},
    {"c_prime", Cons_idx::c_prime}, {"K_T", Cons_idx::K_T}, {"L_wd", Cons_idx::L_wd},
    {"L_ed", Cons_idx::L_ed}, {"D_v", Cons_idx::D_v}, {"ecoll_min", Cons_idx::ecoll_min},
    {"ecoll_gg", Cons_idx::ecoll_gg}, {"ecoll_gg_wet", Cons_idx::ecoll_gg_wet},
    {"kin_visc_air", Cons_idx::kin_visc_air}, {"C_mult", Cons_idx::C_mult},
    {"T_mult_min", Cons_idx::T_mult_min}, {"T_mult_max", Cons_idx::T_mult_max},
    {"T_mult_opt", Cons_idx::T_mult_opt}, {"const0", Cons_idx::const0}, {"const3", Cons_idx::const3},
    {"const4", Cons_idx::const4}, {"const5", Cons_idx::const5}, {"D_rainfrz_ig", Cons_idx::D_rainfrz_ig},
    {"D_rainfrz_gh", Cons_idx::D_rainfrz_gh}, {"dv0", Cons_idx::dv0}, {"p_sat_melt", Cons_idx::p_sat_melt},
    {"cp", Cons_idx::cp}, {"k_b", Cons_idx::k_b}, {"a_HET", Cons_idx::a_HET}, {"b_HET", Cons_idx::b_HET},
    {"N_sc", Cons_idx::N_sc}, {"n_f", Cons_idx::n_f}, {"N_avo", Cons_idx::N_avo}, {"na_dust", Cons_idx::na_dust},
    {"na_soot", Cons_idx::na_soot}, {"na_orga", Cons_idx::na_orga}, {"ni_het_max", Cons_idx::ni_het_max},
    {"ni_hom_max", Cons_idx::ni_hom_max}, {"a_dep", Cons_idx::a_dep}, {"b_dep", Cons_idx::b_dep},
    {"c_dep", Cons_idx::c_dep}, {"d_dep", Cons_idx::d_dep}, {"nim_imm", Cons_idx::nim_imm},
    {"nin_dep", Cons_idx::nin_dep}, {"alf_imm", Cons_idx::alf_imm}, {"bet_dep", Cons_idx::bet_dep},
    {"bet_imm", Cons_idx::bet_imm}, {"r_const", Cons_idx::r_const}, {"r1_const", Cons_idx::r1_const},
    {"cv", Cons_idx::cv}, {"p_sat_const_a", Cons_idx::p_sat_const_a},
    {"p_sat_ice_const_a", Cons_idx::p_sat_ice_const_a}, {"p_sat_const_b", Cons_idx::p_sat_const_b},
    {"p_sat_ice_const_b", Cons_idx::p_sat_ice_const_b}, {"p_sat_low_temp", Cons_idx::p_sat_low_temp},
    {"T_sat_low_temp", Cons_idx::T_sat_low_temp}, {"alpha_depo", Cons_idx::alpha_depo},
    {"r_0", Cons_idx::r_0}, {"k_1_conv", Cons_idx::k_1_conv}, {"k_2_conv", Cons_idx::k_2_conv},
    {"k_1_accr", Cons_idx::k_1_accr}, {"k_r", Cons_idx::k_r},
    {"a_ccn_1", Cons_idx::a_ccn_1}, {"a_ccn_2", Cons_idx::a_ccn_2},
    {"a_ccn_3", Cons_idx::a_ccn_3}, {"a_ccn_4", Cons_idx::a_ccn_4},
    {"b_ccn_1", Cons_idx::b_ccn_1}, {"b_ccn_2", Cons_idx::b_ccn_2},
    {"b_ccn_3", Cons_idx::b_ccn_3}, {"b_ccn_4", Cons_idx::b_ccn_4},
    {"c_ccn_1", Cons_idx::c_ccn_1}, {"c_ccn_2", Cons_idx::c_ccn_2},
    {"c_ccn_3", Cons_idx::c_ccn_3}, {"c_ccn_4", Cons_idx::c_ccn_4},
    {"d_ccn_1", Cons_idx::d_ccn_1}, {"d_ccn_2", Cons_idx::d_ccn_2},
    {"d_ccn_3", Cons_idx::d_ccn_3}, {"d_ccn_4", Cons_idx::d_ccn_4},
#endif
#if defined(B_EIGHT)
    {"p_ccn", Cons_idx::p_ccn}, {"h_ccn_1", Cons_idx::h_ccn_1},
    {"h_ccn_2", Cons_idx::h_ccn_2}, {"h_ccn_3", Cons_idx::h_ccn_3},
    {"g_ccn_1", Cons_idx::g_ccn_1}, {"g_ccn_2", Cons_idx::g_ccn_2},
    {"g_ccn_3", Cons_idx::g_ccn_3}, {"i_ccn_1", Cons_idx::i_ccn_1},
    {"i_ccn_2", Cons_idx::i_ccn_2}, {"hande_ccn_fac", Cons_idx::hande_ccn_fac},
#endif
    {"D_br_threshold", Cons_idx::D_br_threshold}, {"k_br", Cons_idx::k_br},
    {"D_br", Cons_idx::D_br}, {"c_br", Cons_idx::c_br}
};

enum class Particle_cons_idx: uint32_t{
    /**
     * Geometry coefficients.
     */
    a_geo, /*!< Coefficient for diameter size calculation */
    b_geo, /*!< Exponent for diameter size calculation */

    /**
     * Minimum size of particle for mean meass calculation.
     */
    min_x,
    /**
     * Minimum size of particle for CCN activation (cloud) and
     * ice activation (ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_act,
    /**
     * Minimum size of particle for homogenous nucleation (ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_nuc_homo,
    /**
     * Minimum size of particle for heterogeneous nucleation (ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_nuc_hetero,
    /**
     * Minimum size of particle for melting (snow, graupel, ice, hail).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_melt,
    /**
     * Minimum size of particle for evaporation (rain, snow, graupel, ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_evap,
    /**
     * Minimum size of particle for freezing (rain, cloud).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_freezing,
    /**
     * Minimum size of particle for vapor deposition (ice, snow, graupel, hail).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_depo,
    /**
     * Minimum size of particle for ice-ice collision.
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_collision,
    /**
     * Minimum size of particle for different collision processes
     * (snow, rain, ice, snow, graupel).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_collection,
    /**
     * Minimum size of particle for conversion processes (cloud, graupel, ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_conversion,
    /**
     * Minimum size of particle for sedimentation (rain, ice, snow, graupel, hail).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_sedimentation,
    /**
     * Minimum size of particle for riming (cloud, rain, ice, snow, hail, graupel).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    min_x_riming,

    max_x, /*!< Maximum size of particle. */
    sc_theta_q, /*!< Coefficient for ice collision ratio mass. */
    sc_delta_q, /*!< Coefficient for ice collision ratio mass. */
    sc_theta_n, /*!< Coefficient for collision particle number (ice, snow). */
    sc_delta_n, /*!< Coefficient for collision particle number (ice, snow). */
    /**
     * Variance for the assumed Gaussian velocity distributions used in collection and riming processes.
     */
    s_vel,
    a_vel,    /*!< Coefficient for particle velocity. */
    b_vel,    /*!< Exponent for particle velocity. */
    /**
     * Coefficient used in density correction for the increased terminal
     * fall velocity with decreasing air density.
     * \f[ \rho_v = (\rho/\rho_0)^{-\rho_{\text{vel}}} \f]
     */
    rho_v,
    c_z, /*!< Coefficient for 2nd mass moment. */
    sc_coll_n,    /*!< Coefficient in graupel self collection and cloud riming. */
    cmu0, /*!< Coefficient for calculating the shape parameter \f$\mu\f$. */
    cmu1, /*!< Coefficient for calculating the shape parameter \f$\mu\f$. */
    cmu2, /*!< Coefficient for calculating the shape parameter \f$\mu\f$. */
    cmu3, /*!< Constant for calculating the shape parameter \f$\mu\f$. */
    cmu4, /*!< Constant for calculating the shape parameter \f$\mu\f$. */
    cmu5, /*!< Exponent for calculating the shape parameter \f$\mu\f$. */
    alpha, /*!< Constant in rain sedimentation. */
    beta, /*!< Coefficient for rain sedimentation. */
    gamma, /*!< Exponent for rain sedimentation. */
    /**
     * Shape parameter of the generalized \f$\Gamma$\f-distribution.
     * i.e. used in rain sedimentation as coefficient.
     */
    nu,
    /**
     * Right edge of incomplete gamma function,
     * which had been initialized with \f[\text{nm}_1\f].
     */
    g1,
    /**
     * Right edge of incomplete gamma function,
     * which had been initialized with \f[\text{nm}_2\f].
     */
    g2,
    /**
     * Shape parameter of the generalized \f$\Gamma$\f-distribution.
     */
    mu,
    /**
     * Used for initializing the incomplete
     * gamma function lookup table 1.
     * Number of bins.
     */
    nm1,
    /**
     * Used for initializing the incomplete
     * gamma function lookup table 2.
     * Number of bins.
     */
    nm2,
    /**
     * Used for initializing the incomplete
     * gamma function lookup table 3.
     * Number of bins.
     */
    nm3,

    q_crit_c, /*!<  Riming parameter. */
    d_crit_c, /*!<  Riming parameter. */
    ecoll_c,  /*!<  Riming coefficient. */
    /**
     * Coefficient for capacity of particle.
     */
    cap,
    a_ven,    /*!< Vapor deposition coefficient. */
    b_ven,    /*!< Currently unused parameter. */

    c_s,  /*!< Inverse of capacity. Coefficient in evaporation and vapor deposition. */
    a_f, /*!< Constant for average ventilation. Used in melting and ice-vapor processes. */
    b_f, /*!< Coefficient for average ventilation. */

    alfa_n,       /*!<  Sedimentation velocity coefficient. */
    alfa_q,       /*!<  Sedimentation velocity coefficient. */
    lambda,       /*!<  Sedimentation velocity coefficient. */
    vsedi_min,    /*!<  Minimum sedimentation velocity parameter. */
    vsedi_max,    /*!<  Maximum sedimentation velocity parameter. */
    n_items
};

template<
    class float_t,
    class enum_t>
inline float_t get_at(std::vector<float_t> vec, enum_t idx) {
    return vec[static_cast<int>(idx)];
}

template<
    class float_t,
    class enum_t>
inline float_t get_at(std::array<float_t, static_cast<int>(Cons_idx::n_items)> arr, enum_t idx) {
    return arr[static_cast<int>(idx)];
}

/**
 * Mapping of json configuration names to parameters in a particle_model_constants_t
 */
std::unordered_map<std::string, Particle_cons_idx> const table_particle_param = {
    {"a_geo", Particle_cons_idx::a_geo}, {"b_geo", Particle_cons_idx::b_geo}, {"min_x", Particle_cons_idx::min_x},
    {"min_x_act", Particle_cons_idx::min_x_act}, {"min_x_nuc_homo", Particle_cons_idx::min_x_nuc_homo},
    {"min_x_nuc_hetero", Particle_cons_idx::min_x_nuc_hetero}, {"min_x_melt", Particle_cons_idx::min_x_melt},
    {"min_x_evap", Particle_cons_idx::min_x_evap}, {"min_x_freezing", Particle_cons_idx::min_x_freezing},
    {"min_x_depo", Particle_cons_idx::min_x_depo}, {"min_x_collision", Particle_cons_idx::min_x_collision},
    {"min_x_collection", Particle_cons_idx::min_x_collection},
    {"min_x_conversion", Particle_cons_idx::min_x_conversion},
    {"min_x_sedimentation", Particle_cons_idx::min_x_sedimentation},
    {"min_x_riming", Particle_cons_idx::min_x_riming}, {"max_x", Particle_cons_idx::max_x},
    {"sc_theta_q", Particle_cons_idx::sc_theta_q}, {"sc_delta_q", Particle_cons_idx::sc_delta_q},
    {"sc_theta_n", Particle_cons_idx::sc_theta_n}, {"sc_delta_n", Particle_cons_idx::sc_delta_n},
    {"s_vel", Particle_cons_idx::s_vel}, {"a_vel", Particle_cons_idx::a_vel}, {"b_vel", Particle_cons_idx::b_vel},
    {"rho_v", Particle_cons_idx::rho_v}, {"c_z", Particle_cons_idx::c_z},
    {"sc_coll_n", Particle_cons_idx::sc_coll_n}, {"cmu0", Particle_cons_idx::cmu0},
    {"cmu1", Particle_cons_idx::cmu1}, {"cmu2", Particle_cons_idx::cmu2}, {"cmu3", Particle_cons_idx::cmu3},
    {"cmu4", Particle_cons_idx::cmu4}, {"cmu5", Particle_cons_idx::cmu5}, {"alpha", Particle_cons_idx::alpha},
    {"beta", Particle_cons_idx::beta}, {"gamma", Particle_cons_idx::gamma}, {"nu", Particle_cons_idx::nu},
    {"g1", Particle_cons_idx::g1}, {"g2", Particle_cons_idx::g2}, {"mu", Particle_cons_idx::mu},
    {"nm1", Particle_cons_idx::nm1}, {"nm2", Particle_cons_idx::nm2}, {"nm3", Particle_cons_idx::nm3},
    {"q_crit_c", Particle_cons_idx::q_crit_c}, {"d_crit_c", Particle_cons_idx::d_crit_c},
    {"ecoll_c", Particle_cons_idx::ecoll_c}, {"cap", Particle_cons_idx::cap}, {"a_ven", Particle_cons_idx::a_ven},
    {"b_ven", Particle_cons_idx::b_ven}, {"c_s", Particle_cons_idx::c_s}, {"a_f", Particle_cons_idx::a_f},
    {"b_f", Particle_cons_idx::b_f}, {"alfa_n", Particle_cons_idx::alfa_n}, {"alfa_q", Particle_cons_idx::alfa_q},
    {"lambda", Particle_cons_idx::lambda}, {"vsedi_min", Particle_cons_idx::vsedi_min},
    {"vsedi_max", Particle_cons_idx::vsedi_max}
};

extern double sediment_q;
extern double sediment_n;
extern double sediment_q_total;
extern double sediment_n_total;

/////// Various parameters for hydrometeors
#ifdef SB_SHAPE
/**
 * Shape parameter of the generalized \f$\Gamma$\f-distribution.
 */
const double cloud_nu = 1;

/**
 * Shape parameter of the generalized \f$\Gamma$\f-distribution.
 */
const double cloud_mu = 1;
#else
/**
 * Shape parameter of the generalized \f$\Gamma$\f-distribution.
 */
const double cloud_nu = 0;

/**
 * Shape parameter of the generalized \f$\Gamma$\f-distribution.
 */
const double cloud_mu = 1.0/3.0;
#endif

/**
 * Maximum size of cloud droplet.
 */
const double cloud_max_x = 2.6e-10;

/**
 * Minimum size of cloud droplet.
 */
const double cloud_min_x = 4.2e-15;

/**
 *
 */
const double cloud_a_geo = 1.24e-1;

/**
 *
 */
const double cloud_b_geo = 1.0/3.0;

/**
 *
 */
const double cloud_a_vel = 3.75e5;

/**
 *
 */
const double cloud_b_vel = 2.0/3.0;

/**
 *
 */
const double cloud_a_ven = 0.78;

/**
 *
 */
const double cloud_b_ven = 0.308;

/**
 *
 */
const double cloud_cap = 2;

/**
 *
 */
const double cloud_vsedi_max = 1;

/**
 *
 */
const double cloud_vsedi_min = 0;

/**
 *
 */
const double cloud_q_crit_c = 1.0e-6;

/**
 *
 */
const double cloud_d_crit_c = 1.0e-5;

#ifdef SB_SHAPE
/**
 *
 */
const double rain_nu = -2.0/3.0;
#else
/**
 *
 */
const double rain_nu = 0;


#endif
/**
 *
 */
const double rain_mu = 1.0/3.0;

/**
 *
 */
const double rain_max_x = 3.0e-6;

/**
 *
 */
const double rain_min_x = 2.6e-10;

/**
 *
 */
const double rain_a_geo = 1.24e-1;

/**
 *
 */
const double rain_b_geo = 1.0/3.0;

/**
 *
 */
const double rain_a_vel = 114.0137;

/**
 *
 */
const double rain_b_vel = 0.23437;

/**
 *
 */
const double rain_a_ven = 0.78;

/**
 *
 */
const double rain_b_ven = 0.308;

/**
 *
 */
const double rain_cap = 2;

/**
 *
 */
const double rain_vsedi_max = 20;

/**
 *
 */
const double rain_vsedi_min = 0.1;

/**
 *
 */
const double rain_alpha = 9.292;

/**
 *
 */
const double rain_beta = 9.623;

/**
 *
 */
const double rain_gamma = 6.222e2;

/**
 * Depends on mue-D-relation of raindrops.
 * Possible values are
 * 0
 * 6 (default)
 * 11 for increased evaporation
 * 19 (Milbrandt & Yau, 2005)
 */
const double rain_cmu0 = 6;

/**
 * Depends on mue-D-relation of raindrops.
 * Possible values are
 * 0
 * 30 (default)
 * 30 for increased evaporation
 * 19 (Milbrandt & Yau, 2005)
 */
const double rain_cmu1 = 30;

/**
 * Depends on mue-D-relation of raindrops.
 * Possible values are
 * 1
 * 1000 (default)
 * 1000 for increased evaporation
 * 600 (Milbrandt & Yau, 2005)
 */
const double rain_cmu2 = 1000;

/**
 * Depends on mue-D-relation of raindrops.
 * Possible values are
 * 1
 * 1.1e-3 (default)
 * 1.1e-3 for increased evaporation
 * 1.8e-3 (Milbrandt & Yau, 2005)
 */
const double rain_cmu3 = 1.1e-3;

/**
 * Depends on mue-D-relation of raindrops.
 * Possible values are
 * (nu+1)/b_geo - 1
 * 1 (default)
 * 4 for increased evaporation
 * 17 (Milbrandt & Yau, 2005)
 */
const double rain_cmu4 = 1;

/**
 * Depends on mue-D-relation of raindrops.
 * Possible values are
 * 1
 * 2 (default)
 * 2 for increased evaporation
 * 1 (Milbrandt & Yau, 2005)
 */
const double rain_cmu5 = 2;

/**
 * Depends on mue-D-relation of raindrops.
 * Possible values are
 * -1
 * 1 (default)
 * 0.5 for increased evaporation
 * -1 (Milbrandt & Yau, 2005)
 */
const double rain_rain_gfak = 1;

/**
 *
 */
const double graupel_nu = 1;

/**
 *
 */
const double graupel_mu = 1.0/3.0;

/**
 *
 */
const double graupel_max_x = 5.3e-4;

/**
 *
 */
const double graupel_min_x = 4.19e-9;

/**
 *
 */
const double graupel_a_geo = 1.42e-1;

/**
 *
 */
const double graupel_b_geo = 0.314;

/**
 *
 */
const double graupel_a_vel = 86.89371;

/**
 *
 */
const double graupel_b_vel = 0.268325;

/**
 *
 */
const double graupel_a_ven = 0.78;

/**
 *
 */
const double graupel_b_ven = 0.308;

/**
 *
 */
const double graupel_cap = 2;

/**
 *
 */
const double graupel_vsedi_max = 30;

/**
 *
 */
const double graupel_vsedi_min = 0.1;

/**
 *
 */
const double graupel_q_crit_c = 1.0e-6;

/**
 *
 */
const double graupel_d_crit_c = 100.0e-6;

/**
 *
 */
const double graupel_s_vel = 0;

/**
 *
 */
const double graupel_ecoll_c = 1;

/**
 *
 */
const double hail_nu = 1;

/**
 *
 */
const double hail_mu = 1.0/3.0;

/**
 *
 */
const double hail_max_x = 5.4e-4;

/**
 *
 */
const double hail_min_x = 2.6e-9;

/**
 *
 */
const double hail_a_geo = 0.1366;

/**
 *
 */
const double hail_b_geo = 1.0/3.0;

/**
 *
 */
const double hail_a_vel = 39.3;

/**
 *
 */
const double hail_b_vel = 1.0/6.0;

/**
 *
 */
const double hail_a_ven = 0.78;

/**
 *
 */
const double hail_b_ven = 0.308;

/**
 *
 */
const double hail_cap = 2;

/**
 *
 */
const double hail_vsedi_max = 30;

/**
 *
 */
const double hail_vsedi_min = 0.1;

/**
 *
 */
const double hail_q_crit_c = 1.0e-6;

/**
 *
 */
const double hail_d_crit_c = 100.0e-6;

/**
 *
 */
const double hail_s_vel = 0;

/**
 *
 */
const double hail_ecoll_c = 1;

#ifdef SB_SHAPE
/**
 *
 */
const double ice_nu = 1;
#else
/**
 *
 */
const double ice_nu = 0;
#endif
/**
 *
 */
const double ice_mu = 1.0/3.0;

/**
 *
 */
const double ice_max_x = 1.0e-5;

/**
 *
 */
const double ice_min_x = 1.0e-12;

/**
 *
 */
const double ice_a_geo = 0.835;

/**
 *
 */
const double ice_b_geo = 0.39;

/**
 *
 */
const double ice_a_vel = 27.7;

/**
 *
 */
const double ice_b_vel = 0.21579;

/**
 *
 */
const double ice_a_ven = 0.78;

/**
 *
 */
const double ice_b_ven = 0.308;

/**
 *
 */
const double ice_cap = 2;

/**
 *
 */
const double ice_vsedi_max = 3;

/**
 *
 */
const double ice_vsedi_min = 0;

/**
 *
 */
const double ice_q_crit_c = 1.0e-5;

/**
 *
 */
const double ice_d_crit_c = 150.0e-6;

/**
 *
 */
const double ice_s_vel = 0.05;

/**
 *
 */
const double ice_ecoll_c = 0.8;

#ifdef SB_SHAPE
/**
 *
 */
const double snow_nu = 1;

/**
 *
 */
const double snow_mu = 1.0/3.0;
#else
/**
 *
 */
const double snow_nu = 0;

/**
 *
 */
const double snow_mu = 0.5;
#endif

/**
 *
 */
const double snow_max_x = 2.0e-5;

/**
 *
 */
const double snow_min_x = 1.0e-10;

/**
 *
 */
const double snow_a_geo = 2.4;

/**
 *
 */
const double snow_b_geo = 0.455;

/**
 *
 */
const double snow_a_vel = 8.8;

/**
 *
 */
const double snow_b_vel = 0.15;

/**
 *
 */
const double snow_a_ven = 0.78;

/**
 *
 */
const double snow_b_ven = 0.308;

/**
 *
 */
const double snow_cap = 2;

/**
 *
 */
const double snow_vsedi_max = 3;

/**
 *
 */
const double snow_vsedi_min = 0.1;

/**
 *
 */
const double snow_q_crit_c = 1.0e-5;

/**
 *
 */
const double snow_d_crit_c = 150.0e-6;

/**
 *
 */
const double snow_s_vel = 0.25;

/**
 *
 */
const double snow_ecoll_c = 0.8;
/////// More parameters

/**
 * Factor for raindrop breakup.
 */
const double k_br = 1000;

/**
 * Threshold at which raindrops break up during self collection.
 */
const double D_br_threshold = 0.3e-3;

/**
 * Factor for raindrop breakup.
 */
const double D_br = 1.1e-3;

/**
 * Factor for raindrop breakup.
 */
const double c_br = 1;

/**
 * Threshold used during CCN activation from A. Miltenberger.
 * Determines how the number of nuclei is calculated.
 *
 */
const double p_ccn = 80000;

/**
 * Coefficient used during CCN activation from A. Miltenberger.
 */
const double h_ccn_1 = 250;

/**
 * Coefficient used during CCN activation from A. Miltenberger.
 */
const double h_ccn_2 = 7;

/**
 * Used during CCN activation from A. Miltenberger.
 * The number of nuclei that can be activated.
 */
const double h_ccn_3 = 257;

/**
 * Exponent used during CCN activation from A. Miltenberger.
 */
const double g_ccn_1 = 800;

/**
 * Exponent used during CCN activation from A. Miltenberger.
 */
const double g_ccn_2 = 150;

/**
 * Exponent used during CCN activation from A. Miltenberger.
 */
const double g_ccn_3 = 400;

/**
 * Coefficient for the number of CCNs used during CCN activation from
 * A. Miltenberger.
 */
const double i_ccn_1 = 1e6;

/**
 * Minimum number of CCNs used during CCN activation from
 * A. Miltenberger.
 */
const double i_ccn_2 = 10.0e-6;

/**
 * Parameter for scaling Hande CCN activation in general.
 */
const double hande_ccn_fac = 1.0;

/**
 * Universal gas constant, unit: J/(mol*K)
 * Source: http://physics.nist.gov/cuu/Constants/
 */
const double R_universal = 8.3144598;

/**
 * Molar mass of water, unit: kg/mol
 * Source: http://www1.lsbu.ac.uk/water/water_properties.html
 */
const double M_w = 0.018015265;

/**
 * Molar mass of dry air, unit: kg/mol
 * Source: Picard et al, 2008: Revised formula for the density of moist air
 */
const double M_a = 0.02896546;

/**
 * Gas constant for water vapor, unit: \f$\text{J}/\text{K}/\text{kg}\f$
 */
const double R_v = R_universal/M_w;  //  461.51 in ICON

/**
 * Gas constant for dry air, unit: \f$\text{J}/\text{K}/\text{kg}\f$
 */
const double R_a = R_universal/M_a;

/**
 * Quotient of the individual gas constants
 */
const double Epsilon = R_a/R_v;

/**
 * Aerosol particle radius prior to freezing used in homogeneous nucleation.
 */
const double r_0 = 0.25e-6;

/**
 * Depostion coefficient for homogeneous ice nucleation.
 * See Spichtinger & Gierens 2009.
 */
const double alpha_depo = 0.5;

/**
 * Gravitational acceleration (m/s^2)
 */
const double gravity_acc = 9.80665;

/**
 * Treshold (ratio mass) for ice selfcollection
 */
const double q_crit_i = 1.0e-6;
/**
 * Treshold (diameter) for ice selfcollection
 */
const double D_crit_i = 1.0e-4;

/**
 * Threshold (diameter) for ice conversion in selfcollection
 */
const double D_conv_i = 75.0e-6;

/**
 * Threshold (ratio mass) for ice rain riming and snow rain riming
 */
const double q_crit_r = 1.0e-5;
/**
 * Threshold (diameter) for ice rain riming and snow rain riming
 */
const double D_crit_r = 1.0e-4;

/**
 * Threshold (ratio mass) for rain freeze and cloud water
 */
const double q_crit_fr = 1.0e-6;

/**
 * Default threshold (ratio mass) is 1e-4 g/m^3
 */
const double q_crit = 1.0e-9;

/**
 * Threshold (diameter) for conversion snow to graupel, ice to graupel
 */
const double D_conv_sg = 2.0e-4;
/**
 * Threshold (diameter) for conversion snow to graupel, ice to graupel
 */
const double D_conv_ig = 2.0e-4;

/**
 * Minimum mass of conversion due to riming
 */
const double x_conv = 1.0e-10;

/**
 * Upper bound for diameter in collision efficiency
 */
const double D_coll_c = 4.0e-5;

/**
 * Kernel for autoconversion
 */
const double kc_autocon = 9.44e9;

/**
 * Height of the trajectory package [m].
 */
const double parcel_height = 250.0;

/**
 * Lower temperature threshold for ice nucleation, -5°C
 */
const double T_nuc = 268.15;

/**
 * Lower temperature threshold for raindrop freezing
 */
const double T_freeze = 273.15;

/**
 * Lower temperature threshold for (instantaneous) raindrop freezing
 */
const double T_f = 233.0;

/**
 * Equilibrium diameter for Seifert & Beheng (2008), ie Eq. 20.
 */
const double D_eq = 1.1e-3;

/**
 * Density of liquid water in \f$\text{kg}/\text{m}^3\f$
 */
const double rho_w = 1000.0;

/**
 * Norm air density.
 */
const double rho_0 = 1.225;

/**
 * Exponent for density correction
 */
const double rho_vel = 0.4;

/**
 * Exponent for density correction of cloud droplets
 */
const double rho_vel_c = 1;  // 0.2;

/**
 * Density of ice in \f$\text{kg}/\text{m}^3\f$
 */
const double rho_ice = 916.7;


/**
 * Constant used in rain evaporation after Seifert (2008)
 */
const double a_v = 0.78;
/**
 * Coefficient used in rain evaporation after Seifert (2008)
 */
const double b_v = 0.308;

/**
 *  Variuous constants from ICON regarding evaporation from melting ice particles
 */
const double N_Sc = 0.71;
/**
 * Various coefficients from ICON regarding evaporation from melting ice particles
 */const double a_prime = 9.65;
/**
 * Various coefficients from ICON regarding evaporation from melting ice particles
 */const double b_prime = 9.80;
/**
 * Various coefficients from ICON regarding evaporation from melting ice particles
 */const double c_prime = 600;

/**
 * Heat conductivity of air in [J/m/s/K].
 */
const double K_T = 2.4e-2;

/**
 * Latent heat of evaporation of water (water->vapor)
 */
const double L_wd = 2.5008e6;
/**
 * Heat of sublimination ice -> vapor
 */
const double L_ed = 2.8345e6;
/**
 * Heat of fusion ice -> water
 */
const double L_ew = L_ed - L_wd;
/**
 * Diffusivity of water vapor in air at \f$0^\circ\text{C}\f$
 */
const double D_v = 2.22e-5;
/**
 * Min. efficiency for collisions graupel - cloud, ice - cloud, snow - cloud
 */
const double ecoll_min = 0.01;
/**
 * Collision efficiency for graupel selfcollection
 */
const double ecoll_gg = 0.10;
/**
 * Collision efficiency for wet graupel
 */
const double ecoll_gg_wet = 0.40;
/**
 * Kinematic viscosity of dry air in
 * \f$\text{m}^2/\text{s}\f$.
 */
const double kin_visc_air = 1.5e-5;
/**
 * Maxmimum is 0.68.
 */
const double alpha_spacefilling = 0.01;

/**
 * Hallet-Mossop ice multiplication: coefficient for splintering
 */
const double C_mult = 3.5e8;
/**
 * Hallet-Mossop ice multiplication
 */
const double T_mult_min = 265.0;
/**
 * Hallet-Mossop ice multiplication
 */
const double T_mult_max = 270.0;
/**
 * Hallet-Mossop ice multiplication
 */
const double T_mult_opt = 268.0;

/**
 * Constant wether to use ice multiplication.
 */
const bool ice_multiplication = true;
/**
 * Constant wether to enhance riming with melting.
 */
const bool enhanced_melting = true;

/**
 * Size thresholds for partitioning of freezing rain in the hail scheme
 */
const double D_rainfrz_ig = 0.5e-3;
/**
 * Size thresholds for partitioning of freezing rain in the hail scheme
 */
const double D_rainfrz_gh = 1.25e-3;

/**
 * Diffusivity of water vapor in air at \f$0^\circ\text{C}\f$
 */
const double dv0 = 2.22e-5;
/**
 * Saturation pressure at \f$\text{T}=\text{T}_\text{freeze}\f$, called e_3 in ICON.
 */
const double p_sat_melt = 6.1078e2;

/**
 * Specific heat capacity of air at constant pressure in
 * \f$\text{J}/\text{K}/\text{kg}\f$
 */
const double cp = 1004.64;  // COSMO: 1005.7

/**
 * Boltzmann constant in \f$\text{J}/\text{K}\f$
 */
const double k_b = 1.3806504e-23;

// A switch for constant drops
#if CONSTANT_DROP
/**
 * Do not use constant drops.
 */
const bool nuc_c_type = true;
#else
/**
 * Use constant drops.
 */
const bool nuc_c_type = false;
#endif

/**
 * Exponent for rain freeze with data of Barklie and Gokhale (PK page 350).
 */
const double a_HET = 0.65;
/**
 * Coefficient for rain freeze with data of Barklie and Gokhale (PK page 350)
 */
const double b_HET = 200.0;

/**
 * Schmidt number (PK, page 541).
 */
const double N_sc  = 0.71;
/**
 * Exponent of N_sc in the vent-coeff. (PK, page 541)
 */
const double n_f = 0.333;

/**
 * Avogadro number in \f$\text{mol}^{-1}\f$
 */
const double N_avo = 6.02214179e23;
/**
 *  Molar weight of dry air in \f$\text{g}\cdot\text{mol}^{-1}\f$
 */
const double amd = 28.97;
/**
 *  Molar weight of water in \f$\text{g}\cdot\text{mol}^{-1}\f$
 */
const double amw = 18.0154;

/**
 * Constants for Phillips et al. ice nucleation scheme
 * initial number density of dust in \f$\text{m}^{-3}\f$, Phillips08
 */
const double na_dust_phillips = 162.0e3;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * initial number density of soot in \f$\text{m}^{-3}\f$ Phillips08
 */
const double na_soot_phillips = 15.0e6;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * initial number density of organics in \f$\text{m}^{-3}\f$, Phillips08
 */
const double na_orga_phillips = 177.0e6;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * max number of IN between \f$1-10\f$ per liter, i.e. 1d3-10d3
*/
const double ni_het_max = 500.0e3;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * number of liquid aerosols between \f$100-5000\f$ per liter
 */
const double ni_hom_max = 5000.0e3;

/**
 * Parameters for deposition formula (2) of Hande et al.
 */
const double a_dep = 0.27626;
/**
 * Parameters for deposition formula (2) of Hande et al.
 */
const double b_dep = 6.21;
/**
 * parameters for deposition formula (2) of Hande et al.
 */
const double c_dep = -1.3107;
/**
 * Parameters for deposition formula (2) of Hande et al.
 */
const double d_dep = 0.26789;

//// more parameters for Hande et al. nucleation for HDCP2 simulations
#if defined(SPRING)
const double nim_imm = 1.5684e5;    /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double nin_dep = 1.7836e5;    /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_imm = 0.2466;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_dep = 0.0075;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_dep = 2.0341;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_imm = 1.2293;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
#endif
#if defined(SUMMER)
const double nim_imm = 2.9694e4;    /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double nin_dep = 2.6543e4;    /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_imm = 0.2813;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_dep = 0.0020;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_dep = 2.5128;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_imm = 1.1778;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
#endif
#if defined(AUTUMN)
const double nim_imm = 4.9920e4;    /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double nin_dep = 7.7167e4;    /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_imm = 0.2622;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_dep = 0.0406;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_dep = 1.4705;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_imm = 1.2044;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
#endif
#if defined(WINTER)
const double nim_imm = 1.0259e5;    /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double nin_dep = 1.1663e4;    /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_imm = 0.2073;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_dep = 0.0194;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_dep = 1.6943;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_imm = 1.2873;      /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
#endif
#if defined(SPRING95)
const double nim_imm = 1.5684e5 * 17.82; /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double nin_dep = 1.7836e5 * 5.87;  /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_imm = 0.2466;           /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double alf_dep = 0.0075;           /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_dep = 2.0341;           /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
const double bet_imm = 1.2293;           /*!< More parameters for Hande et al. nucleation for HDCP2 simulations */
#endif

/**
 * Different autoconversion types. Use 1 for KB and Beheng (1994), 3 for
 * Seifert & Beheng. 2 is currently not supported.
 */
extern uint32_t auto_type;

/**
 * ccn_activation_hdcp2 after Hande et al (2016)
 */
const std::vector<double> a_ccn = {183230691.161, 0.10147358938,
                                -0.2922395814, 229189886.226};
/**
 * ccn_activation_hdcp2 after Hande et al (2016)
 */
const std::vector<double> b_ccn = {0.0001984051994, 4.473190485e-05,
                                0.0001843225275, 0.0001986158191};
/**
 * ccn_activation_hdcp2 after Hande et al (2016)
 */
const std::vector<double> c_ccn = {16.2420263911, 3.22011836758,
                                13.8499423719, 16.2461600644};
/**
 * ccn_activation_hdcp2 after Hande et al (2016)
 */
const std::vector<double> d_ccn = {287736034.13, 0.6258809883,
                                0.8907491812, 360848977.55};

/** Nucleation types
 * 0: force constant cloud drop number, not implemented
 * <6: ccn_activation_hdcp2 (Hande et al)
 * 6: ccn_activation_sk (Segal & Khain), not implemented
 * >6: SB (2006) from Cosmo 5.2 (cloud_nucleation(..))
 */
const int nuc_type = 5;

/**
 * Coefficient for accretion of qc to qr.
 */
const double k_1_accr = 5.0e-4;
/**
 * Coeficcient for accretion of qc to qr.
 */
const double k_r = 5.78;
// const double k_r = 5.28;

#ifdef SB_CONV
/**
 * Exponent for autoconversion of qc to qr.
 */
const double k_1_conv = 400;
/**
 * Exponent for autoconversion of qc to qr.
 */
const double k_2_conv = 0.7;
#else
/**
 * Exponent for autoconversion of qc to qr.
 */
const double k_1_conv = 600;
/**
 * Exponent for autoconversion of qc to qr.
 */
const double k_2_conv = 0.68;
#endif

/**
 * Initial number density of dust [1/m^3] for nuc_type == 6
 */
const double na_dust = 160e4;
/**
 * Initial number density of dust [1/m^3] for nuc_type == 7 or 5
 */
const double na_dust_2 = 160e4;
/**
 * Initial number density of dust [1/m^3] for nuc_type == 8
 */
const double na_dust_3 = 70e4;

/**
 * Initial number density of soot [1/m^3] for nuc_type == 6
 */
const double na_soot = 30e6;
/**
 * Initial number density of soot [1/m^3] for nuc_type == 7 or 5
 */
const double na_soot_2 = 25e6;
/**
 * Initial number density of soot [1/m^3] for nuc_type == 8
 */
const double na_soot_3 = 0;

/**
 * Initial number density of organics [1/m^3] for nuc_type == 6
 */
const double na_orga = 0;
/**
 * Initial number density of organics [1/m^3] for nuc_type == 7 or 5
 */
const double na_orga_2 = 30e6;
/**
 * Initial number density of organics [1/m^3] for nuc_type == 8
 */
const double na_orga_3 = 0;

/** Use nucleation based either on Hande et al. (true)
 * or Phillips et al. (false). This *should* depend
 * on the nucleation type.
 */
const bool use_hdcp2_het = false;

/**
 * Temperature limit for Phillips et al. nucleation look-up table
 */
const uint32_t t_tmax = 30;

/**
 * Supersaturation limit for Phillips et al. nucleation look-up table
 */
const uint32_t s_smax = 60;

/**
 * Increment for temperature for Phillips et al. nucleation look-up table
 */
const uint32_t t_tstep = 2;

/**
 * Increment for ice supersaturation for Phillips et al. nucleation look-up table
 */
const uint32_t s_sstep = 1;

/**
 * Parameter for saturation adjustment
 */
const double r_const = 287.04;

/**
 * Parameter for saturation adjustment
 */
const double r1_const = 461.5;


/**
 * Specific heat capacity of water vapor at constant pressure in
 * \f$\text{J}/\text{K}/\text{kg}\f$
 */
const double cv = 718.66;

/**
 * Always use saturation adjustement (=True) or only adjust if T > 233K (=False)
 */
const bool always_sat_adj = true;

/**
 * Parameter for saturation adjustment. Constant saturated water vapor pressure
 */
const double p_sat_const_a = 17.2693882;

/**
 * Parameter for saturation adjustment. Constant saturated ice pressure
 */
const double p_sat_ice_const_a = 21.8745584;

/**
 * Parameter for saturation adjustment. Constant saturated water vapor pressure
 */
const double p_sat_const_b = 35.86;

/**
 * Parameter for saturation adjustment. Constant saturated ice pressure
 */
const double p_sat_ice_const_b = 7.66;

/**
 * Parameter for saturation adjustment. Saturated water vapor pressure at T = 233K
 */
const double p_sat_low_temp = 610.78;

/**
 * Parameter for saturation adjustment.
 */
const double T_sat_low_temp = 273.15;

const std::vector<std::vector<double> > afrac_dust = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.99e-06},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 1.01e-06, 9.03e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.38e-07, 1.56e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.13e-07, 7.97e-05, 0.000286},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.55e-06, 2.32e-06, 1.16e-06, 3.22e-07, 1.76e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.23e-06, 8.01e-05, 0.000278, 0.000572},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.42e-06, 7.53e-06, 5.44e-06, 3.38e-06, 1.72e-06, 6.01e-07, 5.98e-08, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.01e-06, 9.32e-05, 0.000277, 0.000566, 0.000934},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.42e-06, 1.61e-05, 1.31e-05, 9.9e-06, 6.93e-06, 4.4e-06, 2.86e-06, 1.49e-06,
     2.75e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.37e-05, 0.00012, 0.000305, 0.00058, 0.000937, 0.00134},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.42e-06, 1.95e-05, 2.43e-05, 2.03e-05, 1.61e-05, 1.22e-05, 1.04e-05, 1.04e-05,
   9.89e-06, 7.63e-06, 4.11e-06, 7.05e-06, 2.93e-06, 6.94e-06, 8.2e-06, 2.97e-05, 7.98e-05, 0.000179, 0.000356,
     0.000624, 0.000962, 0.00136, 0.00177},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.95e-05, 3.97e-05, 3.47e-05, 2.96e-05, 2.43e-05, 2.33e-05, 2.66e-05,
   3.14e-05, 3.38e-05, 2.85e-05, 6.91e-05, 0.000105, 0.000113, 0.000138, 0.000194, 0.000296, 0.000461, 0.000702,
     0.00102, 0.0014, 0.00181, 0.00219},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.95e-05, 4.86e-05, 5.33e-05, 4.73e-05, 4.12e-05, 4.16e-05, 5.07e-05,
   6.43e-05, 7.57e-05, 7.21e-05, 0.000197, 0.000338, 0.00036, 0.000399, 0.000485, 0.000627, 0.000837, 0.00113, 0.00148,
     0.00187, 0.00225, 0.00257},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.86e-05, 7.57e-05, 6.95e-05, 6.27e-05, 6.56e-05, 8.21e-05, 0.000108,
   0.000131, 0.000132, 0.000378, 0.000675, 0.000726, 0.00077, 0.000883, 0.00105, 0.00129, 0.00161, 0.00197, 0.00233,
     0.00265, 0.00285},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.86e-05, 9.34e-05, 9.51e-05, 8.85e-05, 9.45e-05, 0.00012, 0.000159,
   0.000197, 0.000204, 0.000603, 0.0011, 0.00118, 0.00123, 0.00136, 0.00154, 0.0018, 0.00211, 0.00244, 0.00274, 0.00296,
     0.00301},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.86e-05, 9.34e-05, 0.000123, 0.000117, 0.000128, 0.000162, 0.000216,
   0.00027, 0.000284, 0.000857, 0.00158, 0.00169, 0.00175, 0.00188, 0.00207, 0.00231, 0.0026, 0.00287, 0.00307, 0.00313,
     0.00313},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.34e-05, 0.0, 0.000148, 0.000162, 0.000208, 0.000276, 0.000343,
    0.000366, 0.00112, 0.00209, 0.00223, 0.00229, 0.00241, 0.00259, 0.0028, 0.00302, 0.0032, 0.00326, 0.00326, 0.00326},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.34e-05, 0.0, 0.000176, 0.000195, 0.000252, 0.000332, 0.000413,
   0.000444, 0.00138, 0.00258, 0.00276, 0.00282, 0.00291, 0.00306, 0.00321, 0.00334, 0.00338, 0.00338, 0.00339,
     0.00339},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000202, 0.000223, 0.000289, 0.00038, 0.000471,
    0.000511, 0.0016, 0.00302, 0.00324, 0.00327, 0.00334, 0.00343, 0.0035, 0.00352, 0.00352, 0.00352, 0.00352, 0.00352},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000208, 0.00025, 0.000319, 0.000413, 0.000509,
   0.000558, 0.00177, 0.00336, 0.00359, 0.00361, 0.00364, 0.00365, 0.00365, 0.00365, 0.00366, 0.00366, 0.00366,
     0.00366},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000208, 0.000281, 0.000351, 0.000443, 0.000534,
     0.000582, 0.00185, 0.00355, 0.00379, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000208, 0.0003, 0.000387, 0.000475, 0.00056,
   0.000605, 0.00193, 0.00369, 0.00394, 0.00394, 0.00395, 0.00395, 0.00395, 0.00395, 0.00395, 0.00395, 0.00395,
     0.00395},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000208, 0.0003, 0.000426, 0.000508, 0.000587,
     0.000629, 0.002, 0.00383, 0.0041, 0.0041, 0.0041, 0.0041, 0.0041, 0.0041, 0.00411, 0.00411, 0.00411, 0.00411},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0003, 0.0, 0.000545, 0.000616, 0.000654,
     0.00208, 0.00398, 0.00426, 0.00426, 0.00426, 0.00426, 0.00426, 0.00427, 0.00427, 0.00427, 0.00427, 0.00427},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0003, 0.0, 0.000584, 0.000645, 0.00068, 0.00217,
     0.00414, 0.00443, 0.00443, 0.00443, 0.00443, 0.00443, 0.00443, 0.00443, 0.00444, 0.00444, 0.00444},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000625, 0.000676, 0.000708, 0.00225,
     0.0043, 0.0046, 0.0046, 0.0046, 0.00461, 0.00461, 0.00461, 0.00461, 0.00461, 0.00461, 0.00461},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000709, 0.000736, 0.00234,
     0.00447, 0.00478, 0.00478, 0.00479, 0.00479, 0.00479, 0.00479, 0.00479, 0.00479, 0.00479, 0.00479},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000743, 0.000765, 0.00243,
     0.00465, 0.00497, 0.00497, 0.00497, 0.00498, 0.00498, 0.00498, 0.00498, 0.00498, 0.00498, 0.00498},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.000795, 0.00253,
     0.00483, 0.00517, 0.00517, 0.00517, 0.00517, 0.00517, 0.00517, 0.00517, 0.00518, 0.00518, 0.00518},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.000827, 0.00263,
     0.00502, 0.00537, 0.00537, 0.00537, 0.00537, 0.00537, 0.00538, 0.00538, 0.00538, 0.00538, 0.00538},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.00086, 0.00273,
     0.00522, 0.00558, 0.00558, 0.00558, 0.00558, 0.00559, 0.00559, 0.00559, 0.00559, 0.00559, 0.00559},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.000867, 0.00284,
     0.00542, 0.0058, 0.0058, 0.0058, 0.0058, 0.0058, 0.00581, 0.00581, 0.00581, 0.00581, 0.00581},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.000867, 0.00295,
     0.00564, 0.00603, 0.00603, 0.00603, 0.00603, 0.00603, 0.00603, 0.00604, 0.00604, 0.00604, 0.00604},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000867, 0.00303, 0.00586,
     0.00626, 0.00626, 0.00627, 0.00627, 0.00627, 0.00627, 0.00627, 0.00627, 0.00627, 0.00628},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000867, 0.00303, 0.00609,
     0.00651, 0.00651, 0.00651, 0.00651, 0.00651, 0.00652, 0.00652, 0.00652, 0.00652, 0.00652},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00303, 0.00633,
     0.00676, 0.00676, 0.00677, 0.00677, 0.00677, 0.00677, 0.00677, 0.00677, 0.00678, 0.00678},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00303, 0.00637,
     0.00703, 0.00703, 0.00703, 0.00703, 0.00703, 0.00704, 0.00704, 0.00704, 0.00704, 0.00704},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00303, 0.00637, 0.0073,
     0.0073, 0.0073, 0.00731, 0.00731, 0.00731, 0.00731, 0.00731, 0.00732, 0.00732},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00637, 0.0075,
     0.00759, 0.00759, 0.00759, 0.00759, 0.0076, 0.0076, 0.0076, 0.0076, 0.0076},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00637, 0.0075,
     0.00788, 0.00789, 0.00789, 0.00789, 0.00789, 0.00789, 0.0079, 0.0079, 0.0079},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0075,
     0.00819, 0.00819, 0.0082, 0.0082, 0.0082, 0.0082, 0.0082, 0.00821, 0.00821},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0075,
     0.00826, 0.00851, 0.00851, 0.00852, 0.00852, 0.00852, 0.00852, 0.00853, 0.00853},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0075,
     0.00826, 0.00884, 0.00885, 0.00885, 0.00885, 0.00885, 0.00886, 0.00886, 0.00886},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00826,
     0.00911, 0.00919, 0.00919, 0.0092, 0.0092, 0.0092, 0.0092, 0.0092},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00826,
     0.00911, 0.00955, 0.00955, 0.00955, 0.00955, 0.00956, 0.00956, 0.00956},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00826,
     0.00911, 0.00992, 0.00992, 0.00992, 0.00993, 0.00993, 0.00993, 0.00993},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.00911, 0.00988, 0.0103, 0.0103, 0.0103, 0.0103, 0.0103, 0.0103},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.00911, 0.00988, 0.0107, 0.0107, 0.0107, 0.0107, 0.0107, 0.0107},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.00988, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.00988, 0.0111, 0.0116, 0.0116, 0.0116, 0.0116, 0.0116},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.00988, 0.0111, 0.012, 0.012, 0.012, 0.012, 0.012},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0111, 0.0122, 0.0125, 0.0125, 0.0125, 0.0125},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0111, 0.0122, 0.0129, 0.0129, 0.013, 0.013},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0111, 0.0122, 0.0134, 0.0134, 0.0135, 0.0135},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0122, 0.0134, 0.014, 0.014, 0.014},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0122, 0.0134, 0.0145, 0.0145, 0.0145},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 7.71e-05, 0.000525, 0.00101, 0.00135, 0.00181, 0.00245, 0.00331, 0.00451, 0.00615,
   0.00841, 0.0107, 0.0118, 0.013, 0.0143, 0.0158, 0.0174, 0.0192, 0.0211, 0.0233, 0.0256, 0.0282, 0.031, 0.034, 0.0373,
     0.0408},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};

const std::vector<std::vector<double> > afrac_orga = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.03e-12, 1.32e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     8.79e-10, 7.91e-07, 1.71e-07, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     7.93e-09, 2.84e-06, 3.18e-06, 1.2e-06, 6.83e-08, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.42e-08, 6.06e-06, 9.23e-06, 5.96e-06, 2.89e-06, 9.1e-07, 4e-08, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     4.79e-08, 1.01e-05, 1.78e-05, 1.38e-05, 9.1e-06, 5.44e-06, 2.76e-06, 8.43e-07},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     7.56e-08, 1.49e-05, 2.82e-05, 2.43e-05, 1.84e-05, 1.33e-05, 8.96e-06, 5.48e-06},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.05e-07, 2.01e-05, 4.01e-05, 3.66e-05, 2.99e-05, 2.38e-05, 1.83e-05, 1.36e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.32e-07, 2.51e-05, 5.22e-05, 5e-05, 4.31e-05, 3.65e-05, 3.02e-05, 2.45e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.56e-07, 2.97e-05, 6.42e-05, 6.36e-05, 5.71e-05, 5.06e-05, 4.39e-05, 3.77e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.72e-07, 3.33e-05, 7.46e-05, 7.62e-05, 7.08e-05, 6.48e-05, 5.85e-05, 5.22e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.8e-07, 3.55e-05, 8.21e-05, 8.67e-05, 8.28e-05, 7.81e-05, 7.29e-05, 6.69e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.89e-07, 3.7e-05, 8.64e-05, 9.37e-05, 9.22e-05, 8.95e-05, 8.56e-05, 8.11e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.98e-07, 3.84e-05, 8.98e-05, 9.74e-05, 9.74e-05, 9.73e-05, 9.57e-05, 9.32e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.07e-07, 4e-05, 9.34e-05, 0.000101, 0.000101, 0.000101, 0.000101, 0.000101},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.16e-07, 4.15e-05, 9.71e-05, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.26e-07, 4.32e-05, 0.000101, 0.000109, 0.00011, 0.00011, 0.00011, 0.00011},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.34e-07, 4.49e-05, 0.000105, 0.000114, 0.000114, 0.000114, 0.000114, 0.000114},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.34e-07, 4.67e-05, 0.000109, 0.000118, 0.000118, 0.000118, 0.000118, 0.000119},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.34e-07, 4.86e-05, 0.000114, 0.000123, 0.000123, 0.000123, 0.000123, 0.000123},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.34e-07, 4.6e-05, 0.000118, 0.000128, 0.000128, 0.000128, 0.000128, 0.000128},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     2.34e-07, 4.6e-05, 0.000123, 0.000133, 0.000133, 0.000133, 0.000133, 0.000133},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     4.6e-05, 0.000127, 0.000138, 0.000138, 0.000138, 0.000139, 0.000139},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     4.6e-05, 0.000127, 0.000144, 0.000144, 0.000144, 0.000144, 0.000144},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     4.6e-05, 0.000127, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.000127, 0.000152, 0.000156, 0.000156, 0.000156, 0.000156},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.000127, 0.000152, 0.000162, 0.000162, 0.000162, 0.000162},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.000127, 0.000152, 0.000167, 0.000168, 0.000168, 0.000168},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.000152, 0.000167, 0.000175, 0.000175, 0.000175},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.000152, 0.000167, 0.000182, 0.000182, 0.000182},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 2.46e-07, 3.05e-06, 3.99e-06, 5.24e-06, 6.92e-06, 9.18e-06, 1.22e-05, 1.64e-05, 2.2e-05, 2.98e-05, 4.04e-05,
   5.51e-05, 7.55e-05, 0.000104, 0.000133, 0.000147, 0.000162, 0.000179, 0.000199, 0.00022, 0.000243, 0.000269,
     0.000299, 0.000331, 0.000366, 0.000405, 0.000448, 0.000495, 0.000545},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};

const std::vector<std::vector<double> > afrac_soot = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     6.81e-10, 1.37e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     3.58e-08, 3.47e-06, 8.39e-07, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.05e-07, 1.22e-05, 1.36e-05, 5.2e-06, 3.99e-07, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.97e-07, 2.59e-05, 3.94e-05, 2.55e-05, 1.24e-05, 3.98e-06, 2.73e-07, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     3.09e-07, 4.32e-05, 7.57e-05, 5.9e-05, 3.88e-05, 2.33e-05, 1.19e-05, 3.7e-06},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     4.31e-07, 6.36e-05, 0.00012, 0.000103, 7.84e-05, 5.67e-05, 3.82e-05, 2.34e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     5.56e-07, 8.55e-05, 0.000171, 0.000156, 0.000127, 0.000101, 7.79e-05, 5.82e-05},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     6.74e-07, 0.000107, 0.000222, 0.000213, 0.000184, 0.000155, 0.000129, 0.000104},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     7.73e-07, 0.000127, 0.000273, 0.000271, 0.000243, 0.000216, 0.000187, 0.00016},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     8.4e-07, 0.000142, 0.000317, 0.000324, 0.000301, 0.000276, 0.000249, 0.000222},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     8.77e-07, 0.000151, 0.000349, 0.000369, 0.000352, 0.000332, 0.00031, 0.000285},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     9.12e-07, 0.000157, 0.000368, 0.000399, 0.000392, 0.000381, 0.000364, 0.000345},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     9.49e-07, 0.000164, 0.000382, 0.000414, 0.000414, 0.000414, 0.000407, 0.000397},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     9.87e-07, 0.00017, 0.000397, 0.000431, 0.000431, 0.000431, 0.000431, 0.000431},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.03e-06, 0.000177, 0.000413, 0.000448, 0.000448, 0.000448, 0.000448, 0.000448},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.07e-06, 0.000184, 0.00043, 0.000466, 0.000466, 0.000466, 0.000466, 0.000466},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.1e-06, 0.000191, 0.000447, 0.000484, 0.000484, 0.000485, 0.000485, 0.000485},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.1e-06, 0.000199, 0.000464, 0.000504, 0.000504, 0.000504, 0.000504, 0.000504},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.1e-06, 0.000207, 0.000483, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.1e-06, 0.000196, 0.000502, 0.000544, 0.000545, 0.000545, 0.000545, 0.000545},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1.1e-06, 0.000196, 0.000522, 0.000566, 0.000566, 0.000566, 0.000567, 0.000567},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.000196, 0.00054, 0.000589, 0.000589, 0.000589, 0.000589, 0.000589},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.000196, 0.00054, 0.000612, 0.000612, 0.000612, 0.000613, 0.000613},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.000196, 0.00054, 0.000636, 0.000637, 0.000637, 0.000637, 0.000637},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.00054, 0.000646, 0.000662, 0.000662, 0.000662, 0.000662},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.00054, 0.000646, 0.000688, 0.000688, 0.000689, 0.000689},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.00054, 0.000646, 0.000712, 0.000716, 0.000716, 0.000716},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.000646, 0.000712, 0.000744, 0.000744, 0.000745},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.000646, 0.000712, 0.000774, 0.000774, 0.000774},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.37e-07, 2.76e-05, 8.66e-05, 0.000127, 0.000172, 0.000235, 0.000321, 0.000442,
   0.000567, 0.000625, 0.000691, 0.000763, 0.000844, 0.000934, 0.00103, 0.00115, 0.00127, 0.0014, 0.00155, 0.00172,
     0.0019, 0.0021, 0.00231},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};

const uint32_t n_lookup = 2000;
const uint32_t n_lookup_highres = 10000;
const uint32_t n_lookup_hr_dummy = 10;


/** @} */  // end of group constants
