#ifndef CONSTANTS_PHYSICS_H
#define CONSTANTS_PHYSICS_H

#include "codi.hpp"
#include <cmath>
#include "physical_parameterizations.h"
#include "types.h"

////////////////////////////////////////////////////////////////////////////////
// Indices of the parameters
////////////////////////////////////////////////////////////////////////////////
#define p_idx 0
#define T_idx 1
#define w_idx 2
#define S_idx 3
#define qc_idx 4
#define qr_idx 5
#define qv_idx 6
#define Nc_idx 7
#define Nr_idx 8
#define Nv_idx 9
#define qi_idx 10
#define Ni_idx 11
#define vi_idx 12
#define qs_idx 13
#define Ns_idx 14
#define qg_idx 15
#define Ng_idx 16
#define qh_idx 17
#define Nh_idx 18
#define qi_out_idx 19
#define qs_out_idx 20
#define qr_out_idx 21
#define qg_out_idx 22
#define qh_out_idx 23
#define lat_heat_idx 24
#define lat_cool_idx 25
#define num_comp 26
// Those are for a different vector
#define qi_in_idx 0
#define qs_in_idx 1
#define qr_in_idx 2
#define qg_in_idx 3

////////////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////////////

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
//! Lower temperature threshold for (instantaneous) raindrop freezing
const double T_f = 233.0;


const double D_eq = 1.1e-3;

//! [kg/m^3] density of liquid water
const double rho_w = 1000.0;

//! norm air density
const double rho_0 = 1.225;

//! exponent for density correction
const double rho_vel = 0.4;

//! [kg/m^3] density of ice
const double rho_ice = 916.7;

//! gas constant of water vapor from ICON mo_physical_constants
const double R_v = 461.51;

//! Variuous constants from ICON regarding evaporation from melting ice particles
const double a_v = 0.78;
//! Variuous constants from ICON regarding evaporation from melting ice particles
const double b_v = 0.308;
//! Variuous constants from ICON regarding evaporation from melting ice particles
const double N_Sc = 0.71;
//! Variuous constants from ICON regarding evaporation from melting ice particles
const double a_prime = 9.65;
//! Variuous constants from ICON regarding evaporation from melting ice particles
const double b_prime = 9.80;
//! Variuous constants from ICON regarding evaporation from melting ice particles
const double c_prime = 600;
//! heat conductivity of air
const double K_T = 1; // TODO
//! Latent heat of evaporation of water (wasser->dampf)
const double L_wd = 2.5008e6;
//! heat of sublimination ice -> vapor
const double L_ed = 2.8345e6;
//! heat of fusion ice -> water
const double L_ew = L_ed - L_wd;
//! Diffusivity of water vapor in air at 0°C
const double D_v = 2.22e-5;
//! Min. efficiency for collisions graupel - cloud, ice - cloud, snow - cloud
const double ecoll_min = 0.01;
//! Collision efficiency for graupel selfcollection
const double ecoll_gg = 0.10;
//! Collision efficiency for wet graupel
const double ecoll_gg_wet = 0.40;
//! [m^2/s]  kinematic viscosity of dry air from mo_physical_constants.f90
const double kin_visc_air = 1.5e-5;
//! max 0.68
const double alpha_spacefilling = 0.01;
//! Hallet-Mossop ice multiplication
/*!
    Coefficient for splintering
*/
const double C_mult = 3.5e8;
//! Hallet-Mossop ice multiplication
const double T_mult_min = 265.0;
//! Hallet-Mossop ice multiplication
const double T_mult_max = 270.0;
//! Hallet-Mossop ice multiplication
const double T_mult_opt = 268.0;

const double const0 = 1.0/(D_coll_c - D_crit_c);
const double const3 = 1.0/(T_mult_opt - T_mult_min);
const double const4 = 1.0/(T_mult_opt - T_mult_max);
//! const5 = c_w / L_ew ?
const double const5 = alpha_spacefilling * rho_w/rho_ice;

const bool ice_multiplication = true;
const bool enhanced_melting = true;

//! Size thresholds for partitioning of freezing rain in the hail scheme
const double D_rainfrz_ig = 0.5e-3;
//! Size thresholds for partitioning of freezing rain in the hail scheme
const double D_rainfrz_gh = 1.25e-3;

//! Diffusivity of water vapor in air at 0 C
const double dv0 = 2.22e-5;
//! Saturation pressure at T=tmelt, e_3 in ICON
const double p_sat_melt = 6.1078e2;

//! specific heat capacity of air at constant pressure [J/K/kg]
const double cp = 1004.64;

//! Boltzmann constant [J/K]
const double k_b = 1.3806504e-23;

// A switch for constant drops
#if CONSTANT_DROP
//! A switch for constant drops
const bool nuc_c_type = true;
#else
//! A switch for constant drops
const bool nuc_c_type = false;
#endif

double rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2,
    graupel_nm1, graupel_nm2, graupel_g1, graupel_g2;

//! Parameters for rain freeze with data of Barklie and Gokhale (PK page 350)
const double a_HET = 0.65;
//! Parameters for rain freeze with data of Barklie and Gokhale (PK page 350)
const double b_HET = 200.0;

//! Schmidt number (PK, page 541)
const double N_sc  = 0.71;
/**
 * Exponent of N_sc in the vent-coeff. (PK, page 541)
 */
const double n_f = 0.333;

/**
 * Avogadro number [1/mol]
 */
const double N_avo = 6.02214179e23;
/**
 * average gravity
 */
const double grav = 9.80665;
/**
 *  molar weight of dry air
 */
const double amd = 28.97;
/**
 *  molar weight of water [g/mol]
 */
const double amw = 18.0154;

/**
 * Constants for Phillips et al. ice nucleation scheme
 * initial number density of dust [1/m³], Phillips08 value 162e3
 */
const double na_dust = 162.0e3;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * initial number density of soot [1/m³], Phillips08 value 15e6
 */
const double na_soot = 15.0e6;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * initial number density of organics [1/m3], Phillips08 value 177e6
 */
const double na_orga = 177.0e6;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * max number of IN between 1-10 per liter, i.e. 1d3-10d3
*/
const double ni_het_max = 500.0e3;
/**
 * Constants for Phillips et al. ice nucleation scheme
 * number of liquid aerosols between 100-5000 per liter
 */
const double ni_hom_max = 5000.0e3;

/**
 * parameters for deposition formula (2) of Hande et al.
 */
const double a_dep = 0.27626;
/**
 * parameters for deposition formula (2) of Hande et al.
 */
const double b_dep = 6.21;
/**
 * parameters for deposition formula (2) of Hande et al.
 */
const double c_dep = -1.3107;
/**
 * parameters for deposition formula (2) of Hande et al.
 */
const double d_dep = 0.26789;

//// more parameters for Hande et al. nucleation for HDCP2 simulations
#if defined(SPRING)
const double nim_imm = 1.5684e5;
const double nin_dep = 1.7836e5;
const double alf_imm = 0.2466;
const double alf_dep = 0.0075;
const double bet_dep = 2.0341;
const double bet_imm = 1.2293;
#endif
#if defined(SUMMER)
const double nim_imm = 2.9694e4;
const double nin_dep = 2.6543e4;
const double alf_imm = 0.2813;
const double alf_dep = 0.0020;
const double bet_dep = 2.5128;
const double bet_imm = 1.1778;
#endif
#if defined(AUTUMN)
const double nim_imm = 4.9920e4;
const double nin_dep = 7.7167e4;
const double alf_imm = 0.2622;
const double alf_dep = 0.0406;
const double bet_dep = 1.4705;
const double bet_imm = 1.2044;
#endif
#if defined(WINTER)
const double nim_imm = 1.0259e5;
const double nin_dep = 1.1663e4;
const double alf_imm = 0.2073;
const double alf_dep = 0.0194;
const double bet_dep = 1.6943;
const double bet_imm = 1.2873;
#endif
#if defined(SPRING95)
const double nim_imm = 1.5684e5 * 17.82;
const double nin_dep = 1.7836e5 * 5.87;
const double alf_imm = 0.2466;
const double alf_dep = 0.0075;
const double bet_dep = 2.0341;
const double bet_imm = 1.2293;
#endif

uint32_t auto_type = 3;

/**
 * see mo_2mom_mcrph_main.f90 line 830 following
 */
codi::RealReverse rain_gfak = -1.0;

codi::RealReverse cloud_k_au;
codi::RealReverse cloud_k_sc;

/**
 * Kernel for autoconversion
 */
codi::RealReverse kc_autocon = 9.44e9;

/**
 * Inverse layer thickness. Used for sedimentation.
 * In Miltenberger (2016) the trajectories start every 100 m between the surface
 * and 4 km altitude using COSMO-2, which uses a mean spacing of 388 m
 * with 13 m close to the surface and 1190 m at 23 km.
 */
codi::RealReverse inv_z = 1.0/1000.0;

void print_reference_quantities(reference_quantities_t &ref)
{
  std::cout << "\nReference quantities\n"
	    << "--------------------\n"
        << "Temperature: " << ref.Tref << " Kelvin\n"
        << "Pressure: " << ref.pref << " Pascal\n"
        << "Mixing-ratio: " << ref.qref << "\n"
        << "Vertical velocity: " << ref.wref << " meter per second\n"
        << "Time: " << ref.tref << " Second\n"
        << std::endl;
}

void print_constants(model_constants_t &cc)
{

  std::cout << "\nModel constants:\n"
	    << "----------------\n"
	    << "\n"
	    << "Final integration time: " << cc.T_end_prime << " seconds\n"
        << "Nondimensional final integration time: " << cc.T_end << "\n"
	    << "Timestep: " << cc.dt_prime << " seconds\n"
	    << "Snapshot Index: " << cc.snapshot_index << "\n"
	    << "Nondimensional timestep: " << cc.dt << "\n"
	    << "Number of iterations: " << cc.num_steps << "\n"
        << "Number of substeps: " << cc.num_sub_steps << "\n"
	    << "a1_scale: " << cc.a1_scale << "\n"
        << "a2_scale: " << cc.a2_scale << "\n"
        << "e1_scale: " << cc.e1_scale << "\n"
        << "e2_scale: " << cc.e2_scale << "\n"
        << "d_scale: " << cc.d_scale << "\n"
	    << "Scaling factor: " << cc.scaling_fact << "\n"
	    << std::endl;
}

#endif
