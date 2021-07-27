#pragma once

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
#include <boost/property_tree/ptree.hpp>

#include "include/misc/error.h"
#include "include/types/input_parameters_t.h"
#include "include/types/model_constants_t.h"
#include "include/types/param_t.h"
#include "include/types/reference_quantities_t.h"

namespace pt = boost::property_tree;

struct segment_t {
    std::vector<param_t> params;
    double value;
    /**
     * Used to track the value for perturbing in case the tracked value
     * is not exactly equal to value but rather jumps around value.
     */
    double old_value;
    uint32_t n_members;
    int value_name;
    int value_name_sig; /*<! Which parameter had the most significant impact so far >*/
    int out_param; /*<! Which output parameter in case of sensitivity methods >*/
    uint32_t n_segments;
    uint32_t err;
    uint32_t old_sign;
    double duration; /*<! Maximum duration integration time for this ensemble >*/
    std::unordered_map<std::string, std::string> tree_strings;
    bool activated;
    /**
     * Allowed tolerance for value_method. If the value jumps, i.e. due to
     * large time steps, it may jump by tol*value from below/above value
     * to above/below it.
     */
    const double tol = 1e-5;

    enum Method {
        impact_change, sign_flip, value_method, repeated_time
    };
    std::unordered_map<std::string, Method> const table_method = {
        {"impact_change", Method::impact_change}, {"sign_flip", Method::sign_flip},
        {"value", Method::value_method}, {"repeated_time", Method::repeated_time}
    };
    int method;
    enum Param {
        pressure, T, w, S, QC, QR, QV, NCCLOUD, NCRAIN,
        QI, NCICE, QS, NCSNOW, QG, NCGRAUPEL, QH, NCHAIL,
        QI_OUT, QS_OUT, QR_OUT, QG_OUT, QH_OUT, latent_heat,
        latent_cool, NI_OUT, NS_OUT, NR_OUT, NG_OUT,
        NH_OUT, z, Inactive, deposition, sublimination,
        a_1, a_2, e_1, e_2, d, N_c, gamma, beta_c,
        beta_r, delta1, delta2, zeta, rain_gfak, cloud_k_au,
        cloud_k_sc, kc_autocon, inv_z, dw, q_crit_i,
        D_crit_i, D_conv_i, q_crit_r, D_crit_r, q_crit_fr,
        D_coll_c, q_crit, D_conv_sg, D_conv_ig, x_conv,
        parcel_height, alpha_spacefilling, T_nuc, T_freeze,
        T_f, D_eq, rho_w, rho_0, rho_vel, rho_vel_c,
        rho_ice, M_w, M_a, R_universal, Epsilon, gravity_acc,
        R_a, R_v, a_v, b_v, a_prime, b_prime, c_prime,
        K_T, L_wd, L_ed, D_v, ecoll_min, ecoll_gg,
        ecoll_gg_wet, kin_visc_air, C_mult, T_mult_min,
        T_mult_max, T_mult_opt, const0, const3, const4,
        const5, D_rainfrz_gh, D_rainfrz_ig, dv0, p_sat_melt, cp,
        k_b, a_HET, b_HET, N_sc, n_f, N_avo, na_dust,
        na_soot, na_orga, ni_het_max, ni_hom_max, a_dep,
        b_dep, c_dep, d_dep, nim_imm, nin_dep, alf_imm,
        bet_dep, bet_imm, r_const, r1_const, cv, p_sat_const_a,
        p_sat_ice_const_a, p_sat_const_b, p_sat_ice_const_b,
        p_sat_low_temp, T_sat_low_temp, alpha_depo,
        r_0, k_1_conv, k_2_conv, k_1_accr, k_r, a_ccn_1,
        a_ccn_2, a_ccn_3, a_ccn_4, b_ccn_1, b_ccn_2,
        b_ccn_3, b_ccn_4, c_ccn_1, c_ccn_2, c_ccn_3,
        c_ccn_4, d_ccn_1, d_ccn_2, d_ccn_3, d_ccn_4,
        rain_a_geo, rain_b_geo, rain_min_x, rain_min_x_act,
        rain_min_x_nuc_homo, rain_min_x_nuc_hetero,
        rain_min_x_melt, rain_min_x_evap, rain_min_x_freezing,
        rain_min_x_depo, rain_min_x_collision, rain_min_x_collection,
        rain_min_x_conversion, rain_min_x_sedimentation,
        rain_min_x_riming, rain_max_x, rain_sc_theta_q,
        rain_sc_delta_q, rain_sc_theta_n, rain_sc_delta_n,
        rain_s_vel, rain_a_vel, rain_b_vel, rain_rho_v,
        rain_c_z, rain_sc_coll_n, rain_cmu0, rain_cmu1,
        rain_cmu2, rain_cmu3, rain_cmu4, rain_cmu5,
        rain_alpha, rain_beta, rain_gamma, rain_nu,
        rain_g1, rain_g2, rain_mu, rain_nm1, rain_nm2,
        rain_nm3, rain_q_crit_c, rain_d_crit_c, rain_ecoll_c,
        rain_cap, rain_a_ven, rain_b_ven, rain_c_s,
        rain_a_f, rain_b_f, rain_alfa_n, rain_alfa_q,
        rain_lambda, rain_vsedi_min, rain_vsedi_max,
        cloud_a_geo, cloud_b_geo, cloud_min_x, cloud_min_x_act,
        cloud_min_x_nuc_homo, cloud_min_x_nuc_hetero,
        cloud_min_x_melt, cloud_min_x_evap, cloud_min_x_freezing,
        cloud_min_x_depo, cloud_min_x_collision, cloud_min_x_collection,
        cloud_min_x_conversion, cloud_min_x_sedimentation,
        cloud_min_x_riming, cloud_max_x, cloud_sc_theta_q,
        cloud_sc_delta_q, cloud_sc_theta_n, cloud_sc_delta_n,
        cloud_s_vel, cloud_a_vel, cloud_b_vel, cloud_rho_v,
        cloud_c_z, cloud_sc_coll_n, cloud_cmu0, cloud_cmu1,
        cloud_cmu2, cloud_cmu3, cloud_cmu4, cloud_cmu5,
        cloud_alpha, cloud_beta, cloud_gamma, cloud_nu,
        cloud_g1, cloud_g2, cloud_mu, cloud_nm1, cloud_nm2,
        cloud_nm3, cloud_q_crit_c, cloud_d_crit_c,
        cloud_ecoll_c, cloud_cap, cloud_a_ven, cloud_b_ven,
        cloud_c_s, cloud_a_f, cloud_b_f, cloud_alfa_n,
        cloud_alfa_q, cloud_lambda, cloud_vsedi_min,
        cloud_vsedi_max, graupel_a_geo, graupel_b_geo,
        graupel_min_x, graupel_min_x_act, graupel_min_x_nuc_homo,
        graupel_min_x_nuc_hetero, graupel_min_x_melt,
        graupel_min_x_evap, graupel_min_x_freezing,
        graupel_min_x_depo, graupel_min_x_collision,
        graupel_min_x_collection, graupel_min_x_conversion,
        graupel_min_x_sedimentation, graupel_min_x_riming,
        graupel_max_x, graupel_sc_theta_q, graupel_sc_delta_q,
        graupel_sc_theta_n, graupel_sc_delta_n, graupel_s_vel,
        graupel_a_vel, graupel_b_vel, graupel_rho_v,
        graupel_c_z, graupel_sc_coll_n, graupel_cmu0,
        graupel_cmu1, graupel_cmu2, graupel_cmu3, graupel_cmu4,
        graupel_cmu5, graupel_alpha, graupel_beta,
        graupel_gamma, graupel_nu, graupel_g1, graupel_g2,
        graupel_mu, graupel_nm1, graupel_nm2, graupel_nm3,
        graupel_q_crit_c, graupel_d_crit_c, graupel_ecoll_c,
        graupel_cap, graupel_a_ven, graupel_b_ven,
        graupel_c_s, graupel_a_f, graupel_b_f, graupel_alfa_n,
        graupel_alfa_q, graupel_lambda, graupel_vsedi_min,
        graupel_vsedi_max, hail_a_geo, hail_b_geo,
        hail_min_x, hail_min_x_act, hail_min_x_nuc_homo,
        hail_min_x_nuc_hetero, hail_min_x_melt, hail_min_x_evap,
        hail_min_x_freezing, hail_min_x_depo, hail_min_x_collision,
        hail_min_x_collection, hail_min_x_conversion,
        hail_min_x_sedimentation, hail_min_x_riming,
        hail_max_x, hail_sc_theta_q, hail_sc_delta_q,
        hail_sc_theta_n, hail_sc_delta_n, hail_s_vel,
        hail_a_vel, hail_b_vel, hail_rho_v, hail_c_z,
        hail_sc_coll_n, hail_cmu0, hail_cmu1, hail_cmu2,
        hail_cmu3, hail_cmu4, hail_cmu5, hail_alpha,
        hail_beta, hail_gamma, hail_nu, hail_g1, hail_g2,
        hail_mu, hail_nm1, hail_nm2, hail_nm3, hail_q_crit_c,
        hail_d_crit_c, hail_ecoll_c, hail_cap, hail_a_ven,
        hail_b_ven, hail_c_s, hail_a_f, hail_b_f, hail_alfa_n,
        hail_alfa_q, hail_lambda, hail_vsedi_min, hail_vsedi_max,
        ice_a_geo, ice_b_geo, ice_min_x, ice_min_x_act,
        ice_min_x_nuc_homo, ice_min_x_nuc_hetero, ice_min_x_melt,
        ice_min_x_evap, ice_min_x_freezing, ice_min_x_depo,
        ice_min_x_collision, ice_min_x_collection,
        ice_min_x_conversion, ice_min_x_sedimentation,
        ice_min_x_riming, ice_max_x, ice_sc_theta_q,
        ice_sc_delta_q, ice_sc_theta_n, ice_sc_delta_n,
        ice_s_vel, ice_a_vel, ice_b_vel, ice_rho_v,
        ice_c_z, ice_sc_coll_n, ice_cmu0, ice_cmu1,
        ice_cmu2, ice_cmu3, ice_cmu4, ice_cmu5, ice_alpha,
        ice_beta, ice_gamma, ice_nu, ice_g1, ice_g2,
        ice_mu, ice_nm1, ice_nm2, ice_nm3, ice_q_crit_c,
        ice_d_crit_c, ice_ecoll_c, ice_cap, ice_a_ven,
        ice_b_ven, ice_c_s, ice_a_f, ice_b_f, ice_alfa_n,
        ice_alfa_q, ice_lambda, ice_vsedi_min, ice_vsedi_max,
        snow_a_geo, snow_b_geo, snow_min_x, snow_min_x_act,
        snow_min_x_nuc_homo, snow_min_x_nuc_hetero,
        snow_min_x_melt, snow_min_x_evap, snow_min_x_freezing,
        snow_min_x_depo, snow_min_x_collision, snow_min_x_collection,
        snow_min_x_conversion, snow_min_x_sedimentation,
        snow_min_x_riming, snow_max_x, snow_sc_theta_q,
        snow_sc_delta_q, snow_sc_theta_n, snow_sc_delta_n,
        snow_s_vel, snow_a_vel, snow_b_vel, snow_rho_v,
        snow_c_z, snow_sc_coll_n, snow_cmu0, snow_cmu1,
        snow_cmu2, snow_cmu3, snow_cmu4, snow_cmu5,
        snow_alpha, snow_beta, snow_gamma, snow_nu,
        snow_g1, snow_g2, snow_mu, snow_nm1, snow_nm2,
        snow_nm3, snow_q_crit_c, snow_d_crit_c, snow_ecoll_c,
        snow_cap, snow_a_ven, snow_b_ven, snow_c_s,
        snow_a_f, snow_b_f, snow_alfa_n, snow_alfa_q,
        snow_lambda, snow_vsedi_min, snow_vsedi_max
    };
    std::unordered_map<std::string, Param> const table_param = {
        {"pressure", Param::pressure}, {"T", Param::T}, {"w", Param::w}, {"S", Param::S}, {"QC", Param::QC},
        {"QR", Param::QR}, {"QV", Param::QV}, {"NCCLOUD", Param::NCCLOUD}, {"NCRAIN", Param::NCRAIN},
        {"QI", Param::QI}, {"NCICE", Param::NCICE}, {"QS", Param::QS}, {"NCSNOW", Param::NCSNOW}, {"QG", Param::QG},
        {"NCGRAUPEL", Param::NCGRAUPEL}, {"QH", Param::QH}, {"NCHAIL", Param::NCHAIL}, {"QI_OUT", Param::QI_OUT},
        {"QS_OUT", Param::QS_OUT}, {"QR_OUT", Param::QR_OUT}, {"QG_OUT", Param::QG_OUT}, {"QH_OUT", Param::QH_OUT},
        {"latent_heat", Param::latent_heat}, {"latent_cool", Param::latent_cool}, {"NI_OUT", Param::NI_OUT},
        {"NS_OUT", Param::NS_OUT}, {"NR_OUT", Param::NR_OUT}, {"NG_OUT", Param::NG_OUT}, {"NH_OUT", Param::NH_OUT},
        {"z", Param::z}, {"Inactive", Param::Inactive}, {"deposition", Param::deposition},
        {"sublimination", Param::sublimination}, {"a_1", Param::a_1}, {"a_2", Param::a_2}, {"e_1", Param::e_1},
        {"e_2", Param::e_2}, {"d", Param::d}, {"N_c", Param::N_c}, {"gamma", Param::gamma}, {"beta_c", Param::beta_c},
        {"beta_r", Param::beta_r}, {"delta1", Param::delta1}, {"delta2", Param::delta2}, {"zeta", Param::zeta},
        {"rain_gfak", Param::rain_gfak}, {"cloud_k_au", Param::cloud_k_au}, {"cloud_k_sc", Param::cloud_k_sc},
        {"kc_autocon", Param::kc_autocon}, {"inv_z", Param::inv_z}, {"dw", Param::dw}, {"q_crit_i", Param::q_crit_i},
        {"D_crit_i", Param::D_crit_i}, {"D_conv_i", Param::D_conv_i}, {"q_crit_r", Param::q_crit_r},
        {"D_crit_r", Param::D_crit_r}, {"q_crit_fr", Param::q_crit_fr}, {"D_coll_c", Param::D_coll_c},
        {"q_crit", Param::q_crit}, {"D_conv_sg", Param::D_conv_sg}, {"D_conv_ig", Param::D_conv_ig},
        {"x_conv", Param::x_conv}, {"parcel_height", Param::parcel_height},
        {"alpha_spacefilling", Param::alpha_spacefilling}, {"T_nuc", Param::T_nuc}, {"T_freeze", Param::T_freeze},
        {"T_f", Param::T_f}, {"D_eq", Param::D_eq}, {"rho_w", Param::rho_w}, {"rho_0", Param::rho_0},
        {"rho_vel", Param::rho_vel}, {"rho_vel_c", Param::rho_vel_c}, {"rho_ice", Param::rho_ice}, {"M_w", Param::M_w},
        {"M_a", Param::M_a}, {"R_universal", Param::R_universal}, {"Epsilon", Param::Epsilon},
        {"gravity_acc", Param::gravity_acc}, {"R_a", Param::R_a}, {"R_v", Param::R_v}, {"a_v", Param::a_v},
        {"b_v", Param::b_v}, {"a_prime", Param::a_prime}, {"b_prime", Param::b_prime}, {"c_prime", Param::c_prime},
        {"K_T", Param::K_T}, {"L_wd", Param::L_wd}, {"L_ed", Param::L_ed}, {"D_v", Param::D_v},
        {"ecoll_min", Param::ecoll_min}, {"ecoll_gg", Param::ecoll_gg}, {"ecoll_gg_wet", Param::ecoll_gg_wet},
        {"kin_visc_air", Param::kin_visc_air}, {"C_mult", Param::C_mult}, {"T_mult_min", Param::T_mult_min},
        {"T_mult_max", Param::T_mult_max}, {"T_mult_opt", Param::T_mult_opt}, {"const0", Param::const0},
        {"const3", Param::const3}, {"const4", Param::const4}, {"const5", Param::const5},
        {"D_rainfrz_gh", Param::D_rainfrz_gh}, {"D_rainfrz_ig", Param::D_rainfrz_ig}, {"dv0", Param::dv0},
        {"p_sat_melt", Param::p_sat_melt}, {"cp", Param::cp}, {"k_b", Param::k_b}, {"a_HET", Param::a_HET},
        {"b_HET", Param::b_HET}, {"N_sc", Param::N_sc}, {"n_f", Param::n_f}, {"N_avo", Param::N_avo},
        {"na_dust", Param::na_dust}, {"na_soot", Param::na_soot}, {"na_orga", Param::na_orga},
        {"ni_het_max", Param::ni_het_max}, {"ni_hom_max", Param::ni_hom_max}, {"a_dep", Param::a_dep},
        {"b_dep", Param::b_dep}, {"c_dep", Param::c_dep}, {"d_dep", Param::d_dep}, {"nim_imm", Param::nim_imm},
        {"nin_dep", Param::nin_dep}, {"alf_imm", Param::alf_imm}, {"bet_dep", Param::bet_dep},
        {"bet_imm", Param::bet_imm}, {"r_const", Param::r_const}, {"r1_const", Param::r1_const}, {"cv", Param::cv},
        {"p_sat_const_a", Param::p_sat_const_a}, {"p_sat_ice_const_a", Param::p_sat_ice_const_a},
        {"p_sat_const_b", Param::p_sat_const_b}, {"p_sat_ice_const_b", Param::p_sat_ice_const_b},
        {"p_sat_low_temp", Param::p_sat_low_temp}, {"T_sat_low_temp", Param::T_sat_low_temp},
        {"alpha_depo", Param::alpha_depo}, {"r_0", Param::r_0}, {"k_1_conv", Param::k_1_conv},
        {"k_2_conv", Param::k_2_conv}, {"k_1_accr", Param::k_1_accr}, {"k_r", Param::k_r}, {"a_ccn_1", Param::a_ccn_1},
        {"a_ccn_2", Param::a_ccn_2}, {"a_ccn_3", Param::a_ccn_3}, {"a_ccn_4", Param::a_ccn_4},
        {"b_ccn_1", Param::b_ccn_1}, {"b_ccn_2", Param::b_ccn_2}, {"b_ccn_3", Param::b_ccn_3},
        {"b_ccn_4", Param::b_ccn_4}, {"c_ccn_1", Param::c_ccn_1}, {"c_ccn_2", Param::c_ccn_2},
        {"c_ccn_3", Param::c_ccn_3}, {"c_ccn_4", Param::c_ccn_4}, {"d_ccn_1", Param::d_ccn_1},
        {"d_ccn_2", Param::d_ccn_2}, {"d_ccn_3", Param::d_ccn_3}, {"d_ccn_4", Param::d_ccn_4},
        {"rain_a_geo", Param::rain_a_geo}, {"rain_b_geo", Param::rain_b_geo}, {"rain_min_x", Param::rain_min_x},
        {"rain_min_x_act", Param::rain_min_x_act}, {"rain_min_x_nuc_homo", Param::rain_min_x_nuc_homo},
        {"rain_min_x_nuc_hetero", Param::rain_min_x_nuc_hetero}, {"rain_min_x_melt", Param::rain_min_x_melt},
        {"rain_min_x_evap", Param::rain_min_x_evap}, {"rain_min_x_freezing", Param::rain_min_x_freezing},
        {"rain_min_x_depo", Param::rain_min_x_depo}, {"rain_min_x_collision", Param::rain_min_x_collision},
        {"rain_min_x_collection", Param::rain_min_x_collection},
        {"rain_min_x_conversion", Param::rain_min_x_conversion},
        {"rain_min_x_sedimentation", Param::rain_min_x_sedimentation}, {"rain_min_x_riming", Param::rain_min_x_riming},
        {"rain_max_x", Param::rain_max_x}, {"rain_sc_theta_q", Param::rain_sc_theta_q},
        {"rain_sc_delta_q", Param::rain_sc_delta_q}, {"rain_sc_theta_n", Param::rain_sc_theta_n},
        {"rain_sc_delta_n", Param::rain_sc_delta_n}, {"rain_s_vel", Param::rain_s_vel},
        {"rain_a_vel", Param::rain_a_vel}, {"rain_b_vel", Param::rain_b_vel}, {"rain_rho_v", Param::rain_rho_v},
        {"rain_c_z", Param::rain_c_z}, {"rain_sc_coll_n", Param::rain_sc_coll_n}, {"rain_cmu0", Param::rain_cmu0},
        {"rain_cmu1", Param::rain_cmu1}, {"rain_cmu2", Param::rain_cmu2}, {"rain_cmu3", Param::rain_cmu3},
        {"rain_cmu4", Param::rain_cmu4}, {"rain_cmu5", Param::rain_cmu5}, {"rain_alpha", Param::rain_alpha},
        {"rain_beta", Param::rain_beta}, {"rain_gamma", Param::rain_gamma}, {"rain_nu", Param::rain_nu},
        {"rain_g1", Param::rain_g1}, {"rain_g2", Param::rain_g2}, {"rain_mu", Param::rain_mu},
        {"rain_nm1", Param::rain_nm1}, {"rain_nm2", Param::rain_nm2}, {"rain_nm3", Param::rain_nm3},
        {"rain_q_crit_c", Param::rain_q_crit_c}, {"rain_d_crit_c", Param::rain_d_crit_c},
        {"rain_ecoll_c", Param::rain_ecoll_c}, {"rain_cap", Param::rain_cap}, {"rain_a_ven", Param::rain_a_ven},
        {"rain_b_ven", Param::rain_b_ven}, {"rain_c_s", Param::rain_c_s}, {"rain_a_f", Param::rain_a_f},
        {"rain_b_f", Param::rain_b_f}, {"rain_alfa_n", Param::rain_alfa_n}, {"rain_alfa_q", Param::rain_alfa_q},
        {"rain_lambda", Param::rain_lambda}, {"rain_vsedi_min", Param::rain_vsedi_min},
        {"rain_vsedi_max", Param::rain_vsedi_max}, {"cloud_a_geo", Param::cloud_a_geo},
        {"cloud_b_geo", Param::cloud_b_geo}, {"cloud_min_x", Param::cloud_min_x},
        {"cloud_min_x_act", Param::cloud_min_x_act}, {"cloud_min_x_nuc_homo", Param::cloud_min_x_nuc_homo},
        {"cloud_min_x_nuc_hetero", Param::cloud_min_x_nuc_hetero}, {"cloud_min_x_melt", Param::cloud_min_x_melt},
        {"cloud_min_x_evap", Param::cloud_min_x_evap}, {"cloud_min_x_freezing", Param::cloud_min_x_freezing},
        {"cloud_min_x_depo", Param::cloud_min_x_depo}, {"cloud_min_x_collision", Param::cloud_min_x_collision},
        {"cloud_min_x_collection", Param::cloud_min_x_collection},
        {"cloud_min_x_conversion", Param::cloud_min_x_conversion},
        {"cloud_min_x_sedimentation", Param::cloud_min_x_sedimentation},
        {"cloud_min_x_riming", Param::cloud_min_x_riming}, {"cloud_max_x", Param::cloud_max_x},
        {"cloud_sc_theta_q", Param::cloud_sc_theta_q}, {"cloud_sc_delta_q", Param::cloud_sc_delta_q},
        {"cloud_sc_theta_n", Param::cloud_sc_theta_n}, {"cloud_sc_delta_n", Param::cloud_sc_delta_n},
        {"cloud_s_vel", Param::cloud_s_vel}, {"cloud_a_vel", Param::cloud_a_vel}, {"cloud_b_vel", Param::cloud_b_vel},
        {"cloud_rho_v", Param::cloud_rho_v}, {"cloud_c_z", Param::cloud_c_z},
        {"cloud_sc_coll_n", Param::cloud_sc_coll_n}, {"cloud_cmu0", Param::cloud_cmu0},
        {"cloud_cmu1", Param::cloud_cmu1}, {"cloud_cmu2", Param::cloud_cmu2}, {"cloud_cmu3", Param::cloud_cmu3},
        {"cloud_cmu4", Param::cloud_cmu4}, {"cloud_cmu5", Param::cloud_cmu5}, {"cloud_alpha", Param::cloud_alpha},
        {"cloud_beta", Param::cloud_beta}, {"cloud_gamma", Param::cloud_gamma}, {"cloud_nu", Param::cloud_nu},
        {"cloud_g1", Param::cloud_g1}, {"cloud_g2", Param::cloud_g2}, {"cloud_mu", Param::cloud_mu},
        {"cloud_nm1", Param::cloud_nm1}, {"cloud_nm2", Param::cloud_nm2}, {"cloud_nm3", Param::cloud_nm3},
        {"cloud_q_crit_c", Param::cloud_q_crit_c}, {"cloud_d_crit_c", Param::cloud_d_crit_c},
        {"cloud_ecoll_c", Param::cloud_ecoll_c}, {"cloud_cap", Param::cloud_cap}, {"cloud_a_ven", Param::cloud_a_ven},
        {"cloud_b_ven", Param::cloud_b_ven}, {"cloud_c_s", Param::cloud_c_s}, {"cloud_a_f", Param::cloud_a_f},
        {"cloud_b_f", Param::cloud_b_f}, {"cloud_alfa_n", Param::cloud_alfa_n}, {"cloud_alfa_q", Param::cloud_alfa_q},
        {"cloud_lambda", Param::cloud_lambda}, {"cloud_vsedi_min", Param::cloud_vsedi_min},
        {"cloud_vsedi_max", Param::cloud_vsedi_max}, {"graupel_a_geo", Param::graupel_a_geo},
        {"graupel_b_geo", Param::graupel_b_geo}, {"graupel_min_x", Param::graupel_min_x},
        {"graupel_min_x_act", Param::graupel_min_x_act}, {"graupel_min_x_nuc_homo", Param::graupel_min_x_nuc_homo},
        {"graupel_min_x_nuc_hetero", Param::graupel_min_x_nuc_hetero},
        {"graupel_min_x_melt", Param::graupel_min_x_melt}, {"graupel_min_x_evap", Param::graupel_min_x_evap},
        {"graupel_min_x_freezing", Param::graupel_min_x_freezing}, {"graupel_min_x_depo", Param::graupel_min_x_depo},
        {"graupel_min_x_collision", Param::graupel_min_x_collision},
        {"graupel_min_x_collection", Param::graupel_min_x_collection},
        {"graupel_min_x_conversion", Param::graupel_min_x_conversion},
        {"graupel_min_x_sedimentation", Param::graupel_min_x_sedimentation},
        {"graupel_min_x_riming", Param::graupel_min_x_riming}, {"graupel_max_x", Param::graupel_max_x},
        {"graupel_sc_theta_q", Param::graupel_sc_theta_q}, {"graupel_sc_delta_q", Param::graupel_sc_delta_q},
        {"graupel_sc_theta_n", Param::graupel_sc_theta_n}, {"graupel_sc_delta_n", Param::graupel_sc_delta_n},
        {"graupel_s_vel", Param::graupel_s_vel}, {"graupel_a_vel", Param::graupel_a_vel},
        {"graupel_b_vel", Param::graupel_b_vel}, {"graupel_rho_v", Param::graupel_rho_v},
        {"graupel_c_z", Param::graupel_c_z}, {"graupel_sc_coll_n", Param::graupel_sc_coll_n},
        {"graupel_cmu0", Param::graupel_cmu0}, {"graupel_cmu1", Param::graupel_cmu1},
        {"graupel_cmu2", Param::graupel_cmu2}, {"graupel_cmu3", Param::graupel_cmu3},
        {"graupel_cmu4", Param::graupel_cmu4}, {"graupel_cmu5", Param::graupel_cmu5},
        {"graupel_alpha", Param::graupel_alpha}, {"graupel_beta", Param::graupel_beta},
        {"graupel_gamma", Param::graupel_gamma}, {"graupel_nu", Param::graupel_nu}, {"graupel_g1", Param::graupel_g1},
        {"graupel_g2", Param::graupel_g2}, {"graupel_mu", Param::graupel_mu}, {"graupel_nm1", Param::graupel_nm1},
        {"graupel_nm2", Param::graupel_nm2}, {"graupel_nm3", Param::graupel_nm3},
        {"graupel_q_crit_c", Param::graupel_q_crit_c}, {"graupel_d_crit_c", Param::graupel_d_crit_c},
        {"graupel_ecoll_c", Param::graupel_ecoll_c}, {"graupel_cap", Param::graupel_cap},
        {"graupel_a_ven", Param::graupel_a_ven}, {"graupel_b_ven", Param::graupel_b_ven},
        {"graupel_c_s", Param::graupel_c_s}, {"graupel_a_f", Param::graupel_a_f}, {"graupel_b_f", Param::graupel_b_f},
        {"graupel_alfa_n", Param::graupel_alfa_n}, {"graupel_alfa_q", Param::graupel_alfa_q},
        {"graupel_lambda", Param::graupel_lambda}, {"graupel_vsedi_min", Param::graupel_vsedi_min},
        {"graupel_vsedi_max", Param::graupel_vsedi_max}, {"hail_a_geo", Param::hail_a_geo},
        {"hail_b_geo", Param::hail_b_geo}, {"hail_min_x", Param::hail_min_x}, {"hail_min_x_act", Param::hail_min_x_act},
        {"hail_min_x_nuc_homo", Param::hail_min_x_nuc_homo}, {"hail_min_x_nuc_hetero", Param::hail_min_x_nuc_hetero},
        {"hail_min_x_melt", Param::hail_min_x_melt}, {"hail_min_x_evap", Param::hail_min_x_evap},
        {"hail_min_x_freezing", Param::hail_min_x_freezing}, {"hail_min_x_depo", Param::hail_min_x_depo},
        {"hail_min_x_collision", Param::hail_min_x_collision}, {"hail_min_x_collection", Param::hail_min_x_collection},
        {"hail_min_x_conversion", Param::hail_min_x_conversion},
        {"hail_min_x_sedimentation", Param::hail_min_x_sedimentation}, {"hail_min_x_riming", Param::hail_min_x_riming},
        {"hail_max_x", Param::hail_max_x}, {"hail_sc_theta_q", Param::hail_sc_theta_q},
        {"hail_sc_delta_q", Param::hail_sc_delta_q}, {"hail_sc_theta_n", Param::hail_sc_theta_n},
        {"hail_sc_delta_n", Param::hail_sc_delta_n}, {"hail_s_vel", Param::hail_s_vel},
        {"hail_a_vel", Param::hail_a_vel}, {"hail_b_vel", Param::hail_b_vel}, {"hail_rho_v", Param::hail_rho_v},
        {"hail_c_z", Param::hail_c_z}, {"hail_sc_coll_n", Param::hail_sc_coll_n}, {"hail_cmu0", Param::hail_cmu0},
        {"hail_cmu1", Param::hail_cmu1}, {"hail_cmu2", Param::hail_cmu2}, {"hail_cmu3", Param::hail_cmu3},
        {"hail_cmu4", Param::hail_cmu4}, {"hail_cmu5", Param::hail_cmu5}, {"hail_alpha", Param::hail_alpha},
        {"hail_beta", Param::hail_beta}, {"hail_gamma", Param::hail_gamma}, {"hail_nu", Param::hail_nu},
        {"hail_g1", Param::hail_g1}, {"hail_g2", Param::hail_g2}, {"hail_mu", Param::hail_mu},
        {"hail_nm1", Param::hail_nm1}, {"hail_nm2", Param::hail_nm2}, {"hail_nm3", Param::hail_nm3},
        {"hail_q_crit_c", Param::hail_q_crit_c}, {"hail_d_crit_c", Param::hail_d_crit_c},
        {"hail_ecoll_c", Param::hail_ecoll_c}, {"hail_cap", Param::hail_cap}, {"hail_a_ven", Param::hail_a_ven},
        {"hail_b_ven", Param::hail_b_ven}, {"hail_c_s", Param::hail_c_s}, {"hail_a_f", Param::hail_a_f},
        {"hail_b_f", Param::hail_b_f}, {"hail_alfa_n", Param::hail_alfa_n}, {"hail_alfa_q", Param::hail_alfa_q},
        {"hail_lambda", Param::hail_lambda}, {"hail_vsedi_min", Param::hail_vsedi_min},
        {"hail_vsedi_max", Param::hail_vsedi_max}, {"ice_a_geo", Param::ice_a_geo}, {"ice_b_geo", Param::ice_b_geo},
        {"ice_min_x", Param::ice_min_x}, {"ice_min_x_act", Param::ice_min_x_act},
        {"ice_min_x_nuc_homo", Param::ice_min_x_nuc_homo}, {"ice_min_x_nuc_hetero", Param::ice_min_x_nuc_hetero},
        {"ice_min_x_melt", Param::ice_min_x_melt}, {"ice_min_x_evap", Param::ice_min_x_evap},
        {"ice_min_x_freezing", Param::ice_min_x_freezing}, {"ice_min_x_depo", Param::ice_min_x_depo},
        {"ice_min_x_collision", Param::ice_min_x_collision}, {"ice_min_x_collection", Param::ice_min_x_collection},
        {"ice_min_x_conversion", Param::ice_min_x_conversion},
        {"ice_min_x_sedimentation", Param::ice_min_x_sedimentation}, {"ice_min_x_riming", Param::ice_min_x_riming},
        {"ice_max_x", Param::ice_max_x}, {"ice_sc_theta_q", Param::ice_sc_theta_q},
        {"ice_sc_delta_q", Param::ice_sc_delta_q}, {"ice_sc_theta_n", Param::ice_sc_theta_n},
        {"ice_sc_delta_n", Param::ice_sc_delta_n}, {"ice_s_vel", Param::ice_s_vel}, {"ice_a_vel", Param::ice_a_vel},
        {"ice_b_vel", Param::ice_b_vel}, {"ice_rho_v", Param::ice_rho_v}, {"ice_c_z", Param::ice_c_z},
        {"ice_sc_coll_n", Param::ice_sc_coll_n}, {"ice_cmu0", Param::ice_cmu0}, {"ice_cmu1", Param::ice_cmu1},
        {"ice_cmu2", Param::ice_cmu2}, {"ice_cmu3", Param::ice_cmu3}, {"ice_cmu4", Param::ice_cmu4},
        {"ice_cmu5", Param::ice_cmu5}, {"ice_alpha", Param::ice_alpha}, {"ice_beta", Param::ice_beta},
        {"ice_gamma", Param::ice_gamma}, {"ice_nu", Param::ice_nu}, {"ice_g1", Param::ice_g1},
        {"ice_g2", Param::ice_g2}, {"ice_mu", Param::ice_mu}, {"ice_nm1", Param::ice_nm1}, {"ice_nm2", Param::ice_nm2},
        {"ice_nm3", Param::ice_nm3}, {"ice_q_crit_c", Param::ice_q_crit_c}, {"ice_d_crit_c", Param::ice_d_crit_c},
        {"ice_ecoll_c", Param::ice_ecoll_c}, {"ice_cap", Param::ice_cap}, {"ice_a_ven", Param::ice_a_ven},
        {"ice_b_ven", Param::ice_b_ven}, {"ice_c_s", Param::ice_c_s}, {"ice_a_f", Param::ice_a_f},
        {"ice_b_f", Param::ice_b_f}, {"ice_alfa_n", Param::ice_alfa_n}, {"ice_alfa_q", Param::ice_alfa_q},
        {"ice_lambda", Param::ice_lambda}, {"ice_vsedi_min", Param::ice_vsedi_min},
        {"ice_vsedi_max", Param::ice_vsedi_max}, {"snow_a_geo", Param::snow_a_geo}, {"snow_b_geo", Param::snow_b_geo},
        {"snow_min_x", Param::snow_min_x}, {"snow_min_x_act", Param::snow_min_x_act},
        {"snow_min_x_nuc_homo", Param::snow_min_x_nuc_homo}, {"snow_min_x_nuc_hetero", Param::snow_min_x_nuc_hetero},
        {"snow_min_x_melt", Param::snow_min_x_melt}, {"snow_min_x_evap", Param::snow_min_x_evap},
        {"snow_min_x_freezing", Param::snow_min_x_freezing}, {"snow_min_x_depo", Param::snow_min_x_depo},
        {"snow_min_x_collision", Param::snow_min_x_collision}, {"snow_min_x_collection", Param::snow_min_x_collection},
        {"snow_min_x_conversion", Param::snow_min_x_conversion},
        {"snow_min_x_sedimentation", Param::snow_min_x_sedimentation}, {"snow_min_x_riming", Param::snow_min_x_riming},
        {"snow_max_x", Param::snow_max_x}, {"snow_sc_theta_q", Param::snow_sc_theta_q},
        {"snow_sc_delta_q", Param::snow_sc_delta_q}, {"snow_sc_theta_n", Param::snow_sc_theta_n},
        {"snow_sc_delta_n", Param::snow_sc_delta_n}, {"snow_s_vel", Param::snow_s_vel},
        {"snow_a_vel", Param::snow_a_vel}, {"snow_b_vel", Param::snow_b_vel}, {"snow_rho_v", Param::snow_rho_v},
        {"snow_c_z", Param::snow_c_z}, {"snow_sc_coll_n", Param::snow_sc_coll_n}, {"snow_cmu0", Param::snow_cmu0},
        {"snow_cmu1", Param::snow_cmu1}, {"snow_cmu2", Param::snow_cmu2}, {"snow_cmu3", Param::snow_cmu3},
        {"snow_cmu4", Param::snow_cmu4}, {"snow_cmu5", Param::snow_cmu5}, {"snow_alpha", Param::snow_alpha},
        {"snow_beta", Param::snow_beta}, {"snow_gamma", Param::snow_gamma}, {"snow_nu", Param::snow_nu},
        {"snow_g1", Param::snow_g1}, {"snow_g2", Param::snow_g2}, {"snow_mu", Param::snow_mu},
        {"snow_nm1", Param::snow_nm1}, {"snow_nm2", Param::snow_nm2}, {"snow_nm3", Param::snow_nm3},
        {"snow_q_crit_c", Param::snow_q_crit_c}, {"snow_d_crit_c", Param::snow_d_crit_c},
        {"snow_ecoll_c", Param::snow_ecoll_c}, {"snow_cap", Param::snow_cap}, {"snow_a_ven", Param::snow_a_ven},
        {"snow_b_ven", Param::snow_b_ven}, {"snow_c_s", Param::snow_c_s}, {"snow_a_f", Param::snow_a_f},
        {"snow_b_f", Param::snow_b_f}, {"snow_alfa_n", Param::snow_alfa_n}, {"snow_alfa_q", Param::snow_alfa_q},
        {"snow_lambda", Param::snow_lambda}, {"snow_vsedi_min", Param::snow_vsedi_min},
        {"snow_vsedi_max", Param::snow_vsedi_max}
    };

    segment_t();
    void add_param(param_t &param);
    void add_value(double v);
    void add_value_name(std::string n);
    void add_method(std::string m);
    void add_counter(uint32_t c);
    void add_duration(double t);
    void add_out_param(std::string p);
    void add_amount(int a);

    int check();

    /**
     * Check if perturbing is necessary.
     *
     */
    bool perturb_check(
        const model_constants_t &cc,
        const std::vector< std::array<double, num_par > > &gradients,
        const std::vector<codi::RealReverse> &y,
        const double timestep);

    /**
     * Deactivates the segment so it is not used anymore/decreases the amount
     * of possible usages of this segment. If this is a repeated_time segment
     * without a fixed duration it will not decrease the possible amounts of
     * usages. The same applies with duration if keep_repeated is true, which
     * can be used to spawn huge ensembles for a limited time.
     */
    void deactivate(const bool keep_repeated = false);

    void perturb(model_constants_t &cc,
        const reference_quantities_t &ref_quant,
        input_parameters_t &input,
        std::string &descr);

    void put(pt::ptree &ptree) const;

    /**
     * Used to read from a checkpoint file where the mean for the gaussians
     * to draw from is given.
     */
    int from_pt(pt::ptree &ptree, model_constants_t &cc);
};
