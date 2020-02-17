#ifndef CONSTANTS_PHYSICS_H
#define CONSTANTS_PHYSICS_H

#include "codi.hpp"
#include <cmath>
#include "physical_parameterizations.h"
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
#define num_comp 32         /*!< Number of output elements of a model */
// Those are for a different vector
#define qi_in_idx 0         /*!< Ice input index for another vector */
#define qs_in_idx 1         /*!< Snow input index for another vector */
#define qr_in_idx 2         /*!< Rain input index for another vector */
#define qg_in_idx 3         /*!< Graupel input index for another vector */
#define Ni_in_idx 4         /*!< Ice input index for another vector */
#define Ns_in_idx 5         /*!< Snow input index for another vector */
#define Nr_in_idx 6         /*!< Rain input index for another vector */
#define Ng_in_idx 7         /*!< Graupel input index for another vector */
#define num_inflows 8       /*!< Number of parameters for inflowing stuff */

#define num_par 56*6+17    /*!< Number of gradients */

////////////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////////////

/**
 * Used for output derivatives.
 */
const std::string output_par_idx[32] =
    {"p", "T", "w", "S", "qc", "qr", "qv", "Nc", "Nr", "Nv",
     "qi", "Ni", "vi", "qs", "Ns", "qg", "Ng", "qh", "Nh",
     "qiout", "qsout", "qrout", "qgout", "qhout",
     "latent_heat", "latent_cool", "Niout", "Nsout", "Nrout",
     "Ngout", "Nhout", "z"};

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

double rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2,
    graupel_nm1, graupel_nm2, graupel_g1, graupel_g2;

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
 * See mo_2mom_mcrph_main.f90 line 830 following of ICON.
 */
codi::RealReverse rain_gfak = -1.0;

codi::RealReverse cloud_k_au; /*!< Parameter for autoconversion Seifert & Beheng. */
codi::RealReverse cloud_k_sc; /*!< Parameter for autoconversion Seifert & Beheng. */

/**
 * Kernel for autoconversion
 */
codi::RealReverse kc_autocon = 9.44e9;

/**
 * Inverse layer thickness. Used for sedimentation.
 * In Miltenberger (2016) the trajectories start every \f$100 \text{m}\f$
 * between the surface and \f$4 \text{km}\f$ altitude using COSMO-2, which
 * uses a mean spacing of \f$388 \text{m}\f$
 * with \f$13 \text{m}\f$ close to the surface and \f$1190 \text{m}\f$
 * at \f$23 \text{km}\f$.
 */
codi::RealReverse inv_z = 1.0/150.0;

/**
 * ccn_activation_hdcp2 after Hande et al (2015)
 */
std::vector<double> a_ccn = {183230691.161, 0.10147358938,
                                -0.2922395814, 229189886.226};
/**
 * ccn_activation_hdcp2 after Hande et al (2015)
 */
std::vector<double> b_ccn = {0.0001984051994, 4.473190485e-05,
                                0.0001843225275, 0.0001986158191};
/**
 * ccn_activation_hdcp2 after Hande et al (2015)
 */
std::vector<double> c_ccn = {16.2420263911, 3.22011836758,
                                13.8499423719, 16.2461600644};
/**
 * ccn_activation_hdcp2 after Hande et al (2015)
 */
std::vector<double> d_ccn = {287736034.13, 0.6258809883,
                                0.8907491812, 360848977.55};

/** Nucleation types
 * 0: force constant cloud drop number, not implemented
 * <6: ccn_activation_hdcp2 (Hande et al)
 * 6: ccn_activation_sk (Segal & Khain), not implemented
 * >6: SB (2006) from Cosmo 5.2 (cloud_nucleation(..))
 */
const int nuc_type = 7;
/** @} */ // end of group constants

#endif
