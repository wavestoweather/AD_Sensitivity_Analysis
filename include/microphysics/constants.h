#ifndef CONSTANTS_PHYSICS_H
#define CONSTANTS_PHYSICS_H

#include "codi.hpp"
#include <cmath>
#include "types.h"

/** @defgroup constants Constants
 * Various constants for accessing data in the right order and model constants
 * for which most are rather uninteresting (i.e. no AD available).
 * @{
 */

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
#define Nv_idx 9            /*!< Number of water vapor droplets index */
#define qi_idx 10           /*!< Ice mixing ratio index */
#define Ni_idx 11           /*!< Number of ice crystals index */
#define vi_idx 12           /*!< Vertical acceleration of ice index */
#define qs_idx 13           /*!< Snow mixing ratio index */
#define Ns_idx 14           /*!< Number of snow particles index */
#define qg_idx 15           /*!< Graupel mixing ratio index */
#define Ng_idx 16           /*!< Number of graupel particles index */
#define qh_idx 17           /*!< Hail mixing ratio index */
#define Nh_idx 18           /*!< Number of hail particles index */
#define qi_out_idx 19       /*!< Ice mixing ratio precipitation index */
#define qs_out_idx 20       /*!< Snow mixing ratio precipitation index */
#define qr_out_idx 21       /*!< Rain mixing ratio precipitation index */
#define qg_out_idx 22       /*!< Graupel mixing ratio precipitation index */
#define qh_out_idx 23       /*!< Hail mixing ratio precipitation index */
#define lat_heat_idx 24     /*!< Latent heating index */
#define lat_cool_idx 25     /*!< Latent cooling index */
#define Ni_out_idx 26       /*!< Ice particles precipitation index */
#define Ns_out_idx 27       /*!< Snow particles ratio precipitation index */
#define Nr_out_idx 28       /*!< Rain droplets ratio precipitation index */
#define Ng_out_idx 29       /*!< Graupel particles ratio precipitation index */
#define Nh_out_idx 30       /*!< Hail particles ratio precipitation index */
#define z_idx 31            /*!< Altitude */
#define n_inact_idx 32      /*!< Number of inactive nuclei (ie due to being activated before) */
#define depo_idx 33         /*!< Number of deposited nuclei */
#define sub_idx 34          /*!< Sublimination number */

#if defined(RK4_ONE_MOMENT)
#define num_comp 10         /*!< Number of output elements of a model */
#define num_par 12          /*!< Number of gradients */

#elif defined(RK4ICE) || defined(RK4NOICE)
#define num_comp 35         /*!< Number of output elements of a model */
#define num_par 56*6+17     /*!< Number of gradients */

#endif

// Those are for an inflow vector
#define qi_in_idx 0         /*!< Ice input index for another vector */
#define qs_in_idx 1         /*!< Snow input index for another vector */
#define qr_in_idx 2         /*!< Rain input index for another vector */
#define qg_in_idx 3         /*!< Graupel input index for another vector */
#define Ni_in_idx 4         /*!< Ice input index for another vector */
#define Ns_in_idx 5         /*!< Snow input index for another vector */
#define Nr_in_idx 6         /*!< Rain input index for another vector */
#define Ng_in_idx 7         /*!< Graupel input index for another vector */
#define num_inflows 8       /*!< Number of parameters for inflowing stuff */

////////////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////////////

#if defined(RK4_ONE_MOMENT)
/**
 * Used for header files of output parameters.
 */
const std::vector<std::string> output_par_idx =
    {"p", "T", "w", "S", "qc", "qr", "qv", "Nc", "Nr", "Nv"};

/**
 * Used for header files of gradients.
 */
const std::vector<std::string> output_grad_idx =
    {"da_1", "da_2", "de_1", "de_2", "dd", "dN_c", "dgamma", "dbeta_c",
    "dbeta_r", "ddelta1", "ddelta2", "dzeta"};
#elif defined(RK4ICE) || defined(RK4NOICE)
/**
 * Used for header files of output parameters.
 */
const std::vector<std::string> output_par_idx =
    {"p", "T", "w", "S", "qc", "qr", "qv", "Nc", "Nr", "Nv",
     "qi", "Ni", "vi", "qs", "Ns", "qg", "Ng", "qh", "Nh",
     "qiout", "qsout", "qrout", "qgout", "qhout",
     "latent_heat", "latent_cool", "Niout", "Nsout", "Nrout",
     "Ngout", "Nhout", "z", "Inactive", "deposition", "sublimination"};

/**
 * Used for header files of gradients.
 */
const std::vector<std::string> output_grad_idx =
    {"da_1", "da_2", "de_1", "de_2", "dd", "dN_c", "dgamma", "dbeta_c",
    "dbeta_r", "ddelta1", "ddelta2", "dzeta", "dcc.rain_gfak", "dcloud_k_au",
    "dcloud_k_sc", "dkc_autocon", "dinv_z",
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
    "dsnow_vsedi_min", "dsnow_vsedi_max"};
#endif


/**
 * Universal gas constant, unit: J/(mol*K)
 * Source: http://physics.nist.gov/cuu/Constants/
 */
const double R_universal = 8.3144598;

/**
 * Molar mass of water, unit: kg/mol
 * Source: http://www1.lsbu.ac.uk/water/water_properties.html
 */
const double Mw = 0.018015265;

/**
 * Molar mass of dry air, unit: kg/mol
 * Source: Picard et al, 2008: Revised formula for the density of moist air
 */
const double Ma = 0.02896546;

/**
 * Gas constant for water vapor, unit: J/(kg*K)
 */
const double Rv = R_universal/Mw;

/**
 * Gas constant for dry air, unit: J/(kg*K)
 */
const double Ra = R_universal/Ma;

/**
 * Quotient of the individual gas constants
 */
const double Epsilon = Ra/Rv;

/**
 * Gravitational acceleration (m/s^2)
 */
const double gravity_acc = 9.81;

/**
 * Melting temperature of ice/snow
 */
const double tmelt = 273.15;

/**
 * Treshold for ice selfcollection
 */
const double q_crit_i = 1.0e-6;
/**
 * Treshold for ice selfcollection
 */
const double D_crit_i = 1.0e-4;

/**
 * Threshold for ice conversion in selfcollection
 */
const double D_conv_i = 75.0e-6;

/**
 * Threshold for ice rain riming and snow rain riming
 */
const double q_crit_r = 1.0e-5;
/**
 * Threshold for ice rain riming and snow rain riming
 */
const double D_crit_r = 1.0e-4;

/**
 * Threshold for rain freeze and cloud water
 */
const double q_crit_fr = 1.0e-6;
/**
 * Threshold for rain freeze and cloud water
 */
const double q_crit_c = 1.0e-6;

/**
 * Default threshold is 1e-4 g/m^3
 */
const double q_crit = 1.0e-7;

/**
 * Threshold for conversion snow to graupel, ice to graupel
 */
const double D_conv_sg = 2.0e-4;
/**
 * Threshold for conversion snow to graupel, ice to graupel
 */
const double D_conv_ig = 2.0e-4;

/**
 * Minimum mass of conversion due to riming
 */
const double x_conv = 1.0e-10;

/**
 * Threshold for cloud drop collection efficiency
 */
const double D_crit_c = 1.0e-5;

/**
 * Upper bound for diameter in collision efficiency
 */
const double D_coll_c = 4.0e-5;

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
const double rho_vel_c = 0.2;

/**
 * Density of ice in \f$\text{kg}/\text{m}^3\f$
 */
const double rho_ice = 916.7;

/**
 * Gas constant of water vapor from ICON mo_physical_constants
 */
const double R_v = 461.51;

/**
 * Various constants from ICON regarding evaporation from melting ice particles
 */
const double a_v = 0.78;
/**
 * Various constants from ICON regarding evaporation from melting ice particles
 */const double b_v = 0.308;
//! Variuous constants from ICON regarding evaporation from melting ice particles
const double N_Sc = 0.71;
/**
 * Various constants from ICON regarding evaporation from melting ice particles
 */const double a_prime = 9.65;
/**
 * Various constants from ICON regarding evaporation from melting ice particles
 */const double b_prime = 9.80;
/**
 * Various constants from ICON regarding evaporation from melting ice particles
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
 * Kinematic viscosity of dry air from mo_physical_constants.f90 in
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
 * Constant used in cloud riming.
 */
const double const0 = 1.0/(D_coll_c - D_crit_c);
/**
 * Hallet-Mossop ice multiplication.
 * Constant used in ice - x and snow - x riming.
 */
const double const3 = 1.0/(T_mult_opt - T_mult_min);
/**
 * Hallet-Mossop ice multiplication.
 * Constant used in ice - x and snow - x riming.
 */
const double const4 = 1.0/(T_mult_opt - T_mult_max);
/**
 * Constant for conversions ice -> graupel, snow -> graupel,
 * melting (used in riming).
 */
const double const5 = alpha_spacefilling * rho_w/rho_ice;

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
 * Saturation pressure at \f$\text{T}=\text{tmelt}\f$, called e_3 in ICON.
 */
const double p_sat_melt = 6.1078e2;

/**
 * Specific heat capacity of air at constant pressure in
 * \f$\text{J}/\text{K}/\text{kg}\f$
 */
const double cp = 1004.64;

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
 * Parameters for rain freeze with data of Barklie and Gokhale (PK page 350).
 */
const double a_HET = 0.65;
/**
 * Parameters for rain freeze with data of Barklie and Gokhale (PK page 350)
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
 * Average gravity
 */
const double grav = 9.80665;
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
const double na_dust = 162.0e3;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * initial number density of soot in \f$\text{m}^{-3}\f$ Phillips08
 */
const double na_soot = 15.0e6;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * initial number density of organics in \f$\text{m}^{-3}\f$, Phillips08
 */
const double na_orga = 177.0e6;
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
uint32_t auto_type = 3;

/**
 * ccn_activation_hdcp2 after Hande et al (2015)
 */
const std::vector<double> a_ccn = {183230691.161, 0.10147358938,
                                -0.2922395814, 229189886.226};
/**
 * ccn_activation_hdcp2 after Hande et al (2015)
 */
const std::vector<double> b_ccn = {0.0001984051994, 4.473190485e-05,
                                0.0001843225275, 0.0001986158191};
/**
 * ccn_activation_hdcp2 after Hande et al (2015)
 */
const std::vector<double> c_ccn = {16.2420263911, 3.22011836758,
                                13.8499423719, 16.2461600644};
/**
 * ccn_activation_hdcp2 after Hande et al (2015)
 */
const std::vector<double> d_ccn = {287736034.13, 0.6258809883,
                                0.8907491812, 360848977.55};

/** Nucleation types
 * 0: force constant cloud drop number, not implemented
 * <6: ccn_activation_hdcp2 (Hande et al)
 * 6: ccn_activation_sk (Segal & Khain), not implemented
 * >6: SB (2006) from Cosmo 5.2 (cloud_nucleation(..))
 */
const int nuc_type = 7;

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

const std::vector<std::vector<double> > afrac_dust = {
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.99e-06},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.01e-06, 9.03e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.38e-07, 1.56e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.13e-07, 7.97e-05, 0.000286},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.55e-06, 2.32e-06, 1.16e-06, 3.22e-07, 1.76e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.23e-06, 8.01e-05, 0.000278, 0.000572},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.42e-06, 7.53e-06, 5.44e-06, 3.38e-06, 1.72e-06, 6.01e-07, 5.98e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.01e-06, 9.32e-05, 0.000277, 0.000566, 0.000934},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.42e-06, 1.61e-05, 1.31e-05, 9.9e-06, 6.93e-06, 4.4e-06, 2.86e-06, 1.49e-06, 2.75e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.37e-05, 0.00012, 0.000305, 0.00058, 0.000937, 0.00134},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.42e-06, 1.95e-05, 2.43e-05, 2.03e-05, 1.61e-05, 1.22e-05, 1.04e-05, 1.04e-05, 9.89e-06, 7.63e-06, 4.11e-06, 7.05e-06, 2.93e-06, 6.94e-06, 8.2e-06, 2.97e-05, 7.98e-05, 0.000179, 0.000356, 0.000624, 0.000962, 0.00136, 0.00177},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.95e-05, 3.97e-05, 3.47e-05, 2.96e-05, 2.43e-05, 2.33e-05, 2.66e-05, 3.14e-05, 3.38e-05, 2.85e-05, 6.91e-05, 0.000105, 0.000113, 0.000138, 0.000194, 0.000296, 0.000461, 0.000702, 0.00102, 0.0014, 0.00181, 0.00219},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.95e-05, 4.86e-05, 5.33e-05, 4.73e-05, 4.12e-05, 4.16e-05, 5.07e-05, 6.43e-05, 7.57e-05, 7.21e-05, 0.000197, 0.000338, 0.00036, 0.000399, 0.000485, 0.000627, 0.000837, 0.00113, 0.00148, 0.00187, 0.00225, 0.00257},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.86e-05, 7.57e-05, 6.95e-05, 6.27e-05, 6.56e-05, 8.21e-05, 0.000108, 0.000131, 0.000132, 0.000378, 0.000675, 0.000726, 0.00077, 0.000883, 0.00105, 0.00129, 0.00161, 0.00197, 0.00233, 0.00265, 0.00285},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.86e-05, 9.34e-05, 9.51e-05, 8.85e-05, 9.45e-05, 0.00012, 0.000159, 0.000197, 0.000204, 0.000603, 0.0011, 0.00118, 0.00123, 0.00136, 0.00154, 0.0018, 0.00211, 0.00244, 0.00274, 0.00296, 0.00301},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.86e-05, 9.34e-05, 0.000123, 0.000117, 0.000128, 0.000162, 0.000216, 0.00027, 0.000284, 0.000857, 0.00158, 0.00169, 0.00175, 0.00188, 0.00207, 0.00231, 0.0026, 0.00287, 0.00307, 0.00313, 0.00313},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.34e-05, 0.0, 0.000148, 0.000162, 0.000208, 0.000276, 0.000343, 0.000366, 0.00112, 0.00209, 0.00223, 0.00229, 0.00241, 0.00259, 0.0028, 0.00302, 0.0032, 0.00326, 0.00326, 0.00326},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.34e-05, 0.0, 0.000176, 0.000195, 0.000252, 0.000332, 0.000413, 0.000444, 0.00138, 0.00258, 0.00276, 0.00282, 0.00291, 0.00306, 0.00321, 0.00334, 0.00338, 0.00338, 0.00339, 0.00339},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000202, 0.000223, 0.000289, 0.00038, 0.000471, 0.000511, 0.0016, 0.00302, 0.00324, 0.00327, 0.00334, 0.00343, 0.0035, 0.00352, 0.00352, 0.00352, 0.00352, 0.00352},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000208, 0.00025, 0.000319, 0.000413, 0.000509, 0.000558, 0.00177, 0.00336, 0.00359, 0.00361, 0.00364, 0.00365, 0.00365, 0.00365, 0.00366, 0.00366, 0.00366, 0.00366},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000208, 0.000281, 0.000351, 0.000443, 0.000534, 0.000582, 0.00185, 0.00355, 0.00379, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038, 0.0038},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000208, 0.0003, 0.000387, 0.000475, 0.00056, 0.000605, 0.00193, 0.00369, 0.00394, 0.00394, 0.00395, 0.00395, 0.00395, 0.00395, 0.00395, 0.00395, 0.00395, 0.00395},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000208, 0.0003, 0.000426, 0.000508, 0.000587, 0.000629, 0.002, 0.00383, 0.0041, 0.0041, 0.0041, 0.0041, 0.0041, 0.0041, 0.00411, 0.00411, 0.00411, 0.00411},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0003, 0.0, 0.000545, 0.000616, 0.000654, 0.00208, 0.00398, 0.00426, 0.00426, 0.00426, 0.00426, 0.00426, 0.00427, 0.00427, 0.00427, 0.00427, 0.00427},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0003, 0.0, 0.000584, 0.000645, 0.00068, 0.00217, 0.00414, 0.00443, 0.00443, 0.00443, 0.00443, 0.00443, 0.00443, 0.00443, 0.00444, 0.00444, 0.00444},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000625, 0.000676, 0.000708, 0.00225, 0.0043, 0.0046, 0.0046, 0.0046, 0.00461, 0.00461, 0.00461, 0.00461, 0.00461, 0.00461, 0.00461},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000709, 0.000736, 0.00234, 0.00447, 0.00478, 0.00478, 0.00479, 0.00479, 0.00479, 0.00479, 0.00479, 0.00479, 0.00479, 0.00479},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000743, 0.000765, 0.00243, 0.00465, 0.00497, 0.00497, 0.00497, 0.00498, 0.00498, 0.00498, 0.00498, 0.00498, 0.00498, 0.00498},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.000795, 0.00253, 0.00483, 0.00517, 0.00517, 0.00517, 0.00517, 0.00517, 0.00517, 0.00517, 0.00518, 0.00518, 0.00518},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.000827, 0.00263, 0.00502, 0.00537, 0.00537, 0.00537, 0.00537, 0.00537, 0.00538, 0.00538, 0.00538, 0.00538, 0.00538},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.00086, 0.00273, 0.00522, 0.00558, 0.00558, 0.00558, 0.00558, 0.00559, 0.00559, 0.00559, 0.00559, 0.00559, 0.00559},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.000867, 0.00284, 0.00542, 0.0058, 0.0058, 0.0058, 0.0058, 0.0058, 0.00581, 0.00581, 0.00581, 0.00581, 0.00581},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000771, 0.000867, 0.00295, 0.00564, 0.00603, 0.00603, 0.00603, 0.00603, 0.00603, 0.00603, 0.00604, 0.00604, 0.00604, 0.00604},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000867, 0.00303, 0.00586, 0.00626, 0.00626, 0.00627, 0.00627, 0.00627, 0.00627, 0.00627, 0.00627, 0.00627, 0.00628},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000867, 0.00303, 0.00609, 0.00651, 0.00651, 0.00651, 0.00651, 0.00651, 0.00652, 0.00652, 0.00652, 0.00652, 0.00652},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00303, 0.00633, 0.00676, 0.00676, 0.00677, 0.00677, 0.00677, 0.00677, 0.00677, 0.00677, 0.00678, 0.00678},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00303, 0.00637, 0.00703, 0.00703, 0.00703, 0.00703, 0.00703, 0.00704, 0.00704, 0.00704, 0.00704, 0.00704},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00303, 0.00637, 0.0073, 0.0073, 0.0073, 0.00731, 0.00731, 0.00731, 0.00731, 0.00731, 0.00732, 0.00732},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00637, 0.0075, 0.00759, 0.00759, 0.00759, 0.00759, 0.0076, 0.0076, 0.0076, 0.0076, 0.0076},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00637, 0.0075, 0.00788, 0.00789, 0.00789, 0.00789, 0.00789, 0.00789, 0.0079, 0.0079, 0.0079},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0075, 0.00819, 0.00819, 0.0082, 0.0082, 0.0082, 0.0082, 0.0082, 0.00821, 0.00821},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0075, 0.00826, 0.00851, 0.00851, 0.00852, 0.00852, 0.00852, 0.00852, 0.00853, 0.00853},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0075, 0.00826, 0.00884, 0.00885, 0.00885, 0.00885, 0.00885, 0.00886, 0.00886, 0.00886},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00826, 0.00911, 0.00919, 0.00919, 0.0092, 0.0092, 0.0092, 0.0092, 0.0092},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00826, 0.00911, 0.00955, 0.00955, 0.00955, 0.00955, 0.00956, 0.00956, 0.00956},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00826, 0.00911, 0.00992, 0.00992, 0.00992, 0.00993, 0.00993, 0.00993, 0.00993},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00911, 0.00988, 0.0103, 0.0103, 0.0103, 0.0103, 0.0103, 0.0103},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00911, 0.00988, 0.0107, 0.0107, 0.0107, 0.0107, 0.0107, 0.0107},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00988, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111, 0.0111},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00988, 0.0111, 0.0116, 0.0116, 0.0116, 0.0116, 0.0116},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00988, 0.0111, 0.012, 0.012, 0.012, 0.012, 0.012},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0111, 0.0122, 0.0125, 0.0125, 0.0125, 0.0125},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0111, 0.0122, 0.0129, 0.0129, 0.013, 0.013},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0111, 0.0122, 0.0134, 0.0134, 0.0135, 0.0135},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0122, 0.0134, 0.014, 0.014, 0.014},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0122, 0.0134, 0.0145, 0.0145, 0.0145},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 7.71e-05, 0.000525, 0.00101, 0.00135, 0.00181, 0.00245, 0.00331, 0.00451, 0.00615, 0.00841, 0.0107, 0.0118, 0.013, 0.0143, 0.0158, 0.0174, 0.0192, 0.0211, 0.0233, 0.0256, 0.0282, 0.031, 0.034, 0.0373, 0.0408},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};

const std::vector<std::vector<double> > afrac_orga = {
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.03e-12, 1.32e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.79e-10, 7.91e-07, 1.71e-07, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.93e-09, 2.84e-06, 3.18e-06, 1.2e-06, 6.83e-08, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.42e-08, 6.06e-06, 9.23e-06, 5.96e-06, 2.89e-06, 9.1e-07, 4e-08, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.79e-08, 1.01e-05, 1.78e-05, 1.38e-05, 9.1e-06, 5.44e-06, 2.76e-06, 8.43e-07},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.56e-08, 1.49e-05, 2.82e-05, 2.43e-05, 1.84e-05, 1.33e-05, 8.96e-06, 5.48e-06},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.05e-07, 2.01e-05, 4.01e-05, 3.66e-05, 2.99e-05, 2.38e-05, 1.83e-05, 1.36e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.32e-07, 2.51e-05, 5.22e-05, 5e-05, 4.31e-05, 3.65e-05, 3.02e-05, 2.45e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.56e-07, 2.97e-05, 6.42e-05, 6.36e-05, 5.71e-05, 5.06e-05, 4.39e-05, 3.77e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.72e-07, 3.33e-05, 7.46e-05, 7.62e-05, 7.08e-05, 6.48e-05, 5.85e-05, 5.22e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.8e-07, 3.55e-05, 8.21e-05, 8.67e-05, 8.28e-05, 7.81e-05, 7.29e-05, 6.69e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.89e-07, 3.7e-05, 8.64e-05, 9.37e-05, 9.22e-05, 8.95e-05, 8.56e-05, 8.11e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.98e-07, 3.84e-05, 8.98e-05, 9.74e-05, 9.74e-05, 9.73e-05, 9.57e-05, 9.32e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.07e-07, 4e-05, 9.34e-05, 0.000101, 0.000101, 0.000101, 0.000101, 0.000101},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.16e-07, 4.15e-05, 9.71e-05, 0.000105, 0.000105, 0.000105, 0.000105, 0.000105},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.26e-07, 4.32e-05, 0.000101, 0.000109, 0.00011, 0.00011, 0.00011, 0.00011},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.34e-07, 4.49e-05, 0.000105, 0.000114, 0.000114, 0.000114, 0.000114, 0.000114},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.34e-07, 4.67e-05, 0.000109, 0.000118, 0.000118, 0.000118, 0.000118, 0.000119},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.34e-07, 4.86e-05, 0.000114, 0.000123, 0.000123, 0.000123, 0.000123, 0.000123},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.34e-07, 4.6e-05, 0.000118, 0.000128, 0.000128, 0.000128, 0.000128, 0.000128},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.34e-07, 4.6e-05, 0.000123, 0.000133, 0.000133, 0.000133, 0.000133, 0.000133},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.6e-05, 0.000127, 0.000138, 0.000138, 0.000138, 0.000139, 0.000139},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.6e-05, 0.000127, 0.000144, 0.000144, 0.000144, 0.000144, 0.000144},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.6e-05, 0.000127, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000127, 0.000152, 0.000156, 0.000156, 0.000156, 0.000156},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000127, 0.000152, 0.000162, 0.000162, 0.000162, 0.000162},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000127, 0.000152, 0.000167, 0.000168, 0.000168, 0.000168},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000152, 0.000167, 0.000175, 0.000175, 0.000175},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000152, 0.000167, 0.000182, 0.000182, 0.000182},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 2.46e-07, 3.05e-06, 3.99e-06, 5.24e-06, 6.92e-06, 9.18e-06, 1.22e-05, 1.64e-05, 2.2e-05, 2.98e-05, 4.04e-05, 5.51e-05, 7.55e-05, 0.000104, 0.000133, 0.000147, 0.000162, 0.000179, 0.000199, 0.00022, 0.000243, 0.000269, 0.000299, 0.000331, 0.000366, 0.000405, 0.000448, 0.000495, 0.000545},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};

const std::vector<std::vector<double> > afrac_soot = {
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.81e-10, 1.37e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.58e-08, 3.47e-06, 8.39e-07, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.05e-07, 1.22e-05, 1.36e-05, 5.2e-06, 3.99e-07, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.97e-07, 2.59e-05, 3.94e-05, 2.55e-05, 1.24e-05, 3.98e-06, 2.73e-07, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.09e-07, 4.32e-05, 7.57e-05, 5.9e-05, 3.88e-05, 2.33e-05, 1.19e-05, 3.7e-06},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.31e-07, 6.36e-05, 0.00012, 0.000103, 7.84e-05, 5.67e-05, 3.82e-05, 2.34e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.56e-07, 8.55e-05, 0.000171, 0.000156, 0.000127, 0.000101, 7.79e-05, 5.82e-05},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.74e-07, 0.000107, 0.000222, 0.000213, 0.000184, 0.000155, 0.000129, 0.000104},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.73e-07, 0.000127, 0.000273, 0.000271, 0.000243, 0.000216, 0.000187, 0.00016},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.4e-07, 0.000142, 0.000317, 0.000324, 0.000301, 0.000276, 0.000249, 0.000222},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.77e-07, 0.000151, 0.000349, 0.000369, 0.000352, 0.000332, 0.00031, 0.000285},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.12e-07, 0.000157, 0.000368, 0.000399, 0.000392, 0.000381, 0.000364, 0.000345},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.49e-07, 0.000164, 0.000382, 0.000414, 0.000414, 0.000414, 0.000407, 0.000397},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.87e-07, 0.00017, 0.000397, 0.000431, 0.000431, 0.000431, 0.000431, 0.000431},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.03e-06, 0.000177, 0.000413, 0.000448, 0.000448, 0.000448, 0.000448, 0.000448},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.07e-06, 0.000184, 0.00043, 0.000466, 0.000466, 0.000466, 0.000466, 0.000466},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1e-06, 0.000191, 0.000447, 0.000484, 0.000484, 0.000485, 0.000485, 0.000485},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1e-06, 0.000199, 0.000464, 0.000504, 0.000504, 0.000504, 0.000504, 0.000504},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1e-06, 0.000207, 0.000483, 0.000524, 0.000524, 0.000524, 0.000524, 0.000524},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1e-06, 0.000196, 0.000502, 0.000544, 0.000545, 0.000545, 0.000545, 0.000545},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1e-06, 0.000196, 0.000522, 0.000566, 0.000566, 0.000566, 0.000567, 0.000567},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000196, 0.00054, 0.000589, 0.000589, 0.000589, 0.000589, 0.000589},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000196, 0.00054, 0.000612, 0.000612, 0.000612, 0.000613, 0.000613},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000196, 0.00054, 0.000636, 0.000637, 0.000637, 0.000637, 0.000637},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00054, 0.000646, 0.000662, 0.000662, 0.000662, 0.000662},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00054, 0.000646, 0.000688, 0.000688, 0.000689, 0.000689},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00054, 0.000646, 0.000712, 0.000716, 0.000716, 0.000716},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000646, 0.000712, 0.000744, 0.000744, 0.000745},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000646, 0.000712, 0.000774, 0.000774, 0.000774},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.37e-07, 2.76e-05, 8.66e-05, 0.000127, 0.000172, 0.000235, 0.000321, 0.000442, 0.000567, 0.000625, 0.000691, 0.000763, 0.000844, 0.000934, 0.00103, 0.00115, 0.00127, 0.0014, 0.00155, 0.00172, 0.0019, 0.0021, 0.00231},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};

/**
 * Structure to hold the new equidistant lookup table for
 * graupel wetgrowth diameter
 */
table_t ltabdminwgg;
gamma_table_t table_g1, table_g2, table_r1, table_r2, table_r3;

const uint64_t n_lookup = 2000;
const uint64_t n_lookup_highres = 10000;
const uint64_t n_lookup_hr_dummy = 10;

/**
 * Used for writing outputs.
 */
std::stringstream out_tmp;
std::ofstream outfile;
std::ofstream out_diff[num_comp];
std::stringstream out_diff_tmp[num_comp];

/** @} */ // end of group constants

#endif
