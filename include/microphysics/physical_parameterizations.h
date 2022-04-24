#pragma once

#include <algorithm>
#include <cmath>

#include <boost/math/special_functions/gamma.hpp>
#include "codi.hpp"

#include "include/microphysics/constants.h"
#include "include/types/collection_model_constants_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/particle_model_constants_t.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/table_t.h"


/** @defgroup parametrizations Physical Parametrizations
 * Various functions for physical parametrizations.
 * @{
 */


/**
 * Evaluate a polynomial
 * \f[ a_0 + a_1*x + a_2*x^2 \f]
 * of order 2 using Horner's method.
 *
 * @param a0 Coefficient of the polynomial.
 * @param a1 Coefficient of the polynomial.
 * @param a2 Coefficient of the polynomial.
 * @param x Value at which to evaluate the polynomial.
 * @return Value at x
 */
template <class float_t>
inline float_t polyval2(
    const double a0,
    const double a1,
    const double a2,
    float_t x) {

    return (a0 + x*(a1 + a2*x));
}


/**
 * Evaluate a polynomial
 * \f[ a_0 + a_1*x + a_2*x^2 + a_3*x^3 + a_4*x^4 \f]
 * of order 4 using Horner's method.
 *
 * @param a0 Coefficient of the polynomial.
 * @param a1 Coefficient of the polynomial.
 * @param a2 Coefficient of the polynomial.
 * @param a3 Coefficient of the polynomial.
 * @param a4 Coefficient of the polynomial.
 * @param x Value at which to evaluate the polynomial.
 * @return Value at x
 */
template <class float_t>
inline float_t polyval4(
    const double a0,
    const double a1,
    const double a2,
    const double a3,
    const double a4,
    float_t x) {

    return (a0 + x*(a1 + x*(a2 + x*(a3 + x*a4))));
}


/**
 * Evaluate a polynomial
 * \f[ a_0 + a_1*x + a_2*x^2 + a_3*x^3 + a_4*x^4 + a_5*x^5 \f]
 * of order 5 using Horner's method.
 *
 * @param a0 Coefficient of the polynomial.
 * @param a1 Coefficient of the polynomial.
 * @param a2 Coefficient of the polynomial.
 * @param a3 Coefficient of the polynomial.
 * @param a4 Coefficient of the polynomial.
 * @param a5 Coefficient of the polynomial.
 * @param x Value at which to evaluate the polynomial.
 * @return Value at x
 */
template <class float_t>
inline float_t polyval5(
    const double a0,
    const double a1,
    const double a2,
    const double a3,
    const double a4,
    const double a5,
    float_t x) {

    return (a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*a5)))));
}


/**
 * Evaluate a polynomial
 * \f[ a_0 + a_1*x + a_2*x^2 + a_3*x^3 + a_4*x^4 + a_5*x^5 + a_6*x^6 \f]
 * of order 6 using Horner's method.
 *
 * @param a0 Coefficient of the polynomial.
 * @param a1 Coefficient of the polynomial.
 * @param a2 Coefficient of the polynomial.
 * @param a3 Coefficient of the polynomial.
 * @param a4 Coefficient of the polynomial.
 * @param a5 Coefficient of the polynomial.
 * @param a6 Coefficient of the polynomial.
 * @param x Value at which to evaluate the polynomial.
 * @return Value at x
 */
template <class float_t>
inline float_t polyval6(
    const double a0,
    const double a1,
    const double a2,
    const double a3,
    const double a4,
    const double a5,
    const double a6,
    float_t x) {

    return (a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + a6*x))))));
}


/**
 * Diffusivity of dry air in \f$ \text{m}^2/\text{s} \f$.
 * From Hall and Pruppacher (1976)
 * Validity range: \f$ 193.15 \text{K} <= \text{T} <= 313.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin.
 * @param p Pressure in Pascal.
 * @return Diffusivity
 */
template <class float_t>
inline float_t diffusivity(
    float_t T,
    float_t p) {

    return ((2.11e-5)*(101325.0/p)*pow(T/273.15 , 1.94));
}


/**
 * Thermal conductivity of dry air in \f$ \text{W}/(\text{m} \cdot \text{K}) \f$.
 * From Kannuluik and Carman (1951)
 * Validity range: \f$ 90.15 \text{K} <= \text{T} <= 491.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Thermal conductivity
 */
template <class float_t>
inline float_t thermal_conductivity_dry_air(float_t T) {
    // Change temperature to Celsius-scale
    float_t T_cel = T - 273.15;
    return (418.68*(5.75e-5)*polyval2(1.0, 0.00317, -0.0000021, T_cel));
}


/**
 * Thermal conductivity of moist air in \f$ \text{W}/(\text{m} \cdot \text{K}) \f$.
 * From Pruppacher and Klett (1997), Equation 13-18c
 * Validity range: ?
 *
 * @param T Temperature in Kelvin
 * @param qv Mixing-ratio of water vapor
 * @return Thermal conductivity
 */
template <class float_t>
inline float_t thermal_conductivity_moist_air(
    float_t T,
    float_t qv,
    float_t Epsilon) {

    // Thermal conductivity dry air
    float_t Kt = thermal_conductivity_dry_air(T);

    float_t Kt_tilde = Kt/418.68;
    float_t Kv_tilde = (3.78e-5) + (2.0e-7)*(T - 273.15);

    return (Kt*(1.0 - (1.17 - 1.02*(Kt_tilde/Kv_tilde))*(qv/(qv+Epsilon))));
}


/**
 * Density of liquid water in \f$ \text{kg}/\text{m}^3 \f$.
 * From Hare and Sorensen (1987); Kell (1975)
 * Validity range: \f$ 239.74 \text{K} <= \text{T} <= 423.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Density
 */
template <class float_t>
inline float_t density_water(float_t T) {
    // Change to Celsius scale
    float_t T_cel = T - 273.15;
    float_t F_T;
    float_t denom;
    double a0, a1, a2, a3, a4, a5, a6;

    if (273.15 >= T) {
        // Parameterization from Hare and Sorensen (1987)
        a0 = 999.86;
        a1 = 6.69e-2;
        a2 = -8.486e-3;
        a3 = 1.518e-4;
        a4 = -6.9484e-6;
        a5 = -3.6449e-7;
        a6 = -7.497e-9;

        F_T = polyval6(a0, a1, a2, a3, a4, a5, a6, T_cel);
        denom = 1.0;

        } else {
        // Parameterization from Kell (1975)

        a0 = 999.8395;
        a1 = 16.945176;
        a2 = -7.9870401e-3;
        a3 = -46.170461e-6;
        a4 = 105.56302e-9;
        a5 = -280.54253e-12;

        F_T = polyval5(a0, a1, a2, a3, a4, a5, T_cel);
        denom = 1.0 + (16.87985e-3)*T_cel;
    }
    return (F_T/denom);
}


/**
 * Density of ice in \f$ \text{kg}/\text{m}^3 \f$.
 * From Bogdan (1997)
 * Validity range: \f$ 228 \text{K} <= \text{T} <= 273.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Density
 */
template <class float_t>
inline float_t density_ice(float_t T) {
    float_t T_cel = T - 273.15;
    return polyval2(916.7, -0.175, -0.0005, T_cel);
}


/**
 * Specific heat capacity of dry air in \f$ \text{J}/(\text{kg} \cdot \text{K}) \f$.
 * From Rogers and Yau (1989)
 * For now we use a fixed constant that is valid in constant pressure.
 * Validity range: ?
 *
 * @param T Temperature in Kelvin
 * @return Specific heat capacity
 */
template <class float_t>
inline float_t specific_heat_dry_air(float_t T) {
    return cp;  // 1005.0
}


/**
 * Specific heat capacity of water vapor in \f$ \text{J}/(\text{kg} \cdot \text{K}) \f$.
 * From Pruppacher and Klett (1997)
 * Validity range: ?
 *
 * @param T Temperature in Kelvin
 * @return Specific heat capacity
 */
template <class float_t>
inline float_t specific_heat_water_vapor() {
    return 1884.06;
}


/**
 * Specific heat capacity of liquid water in \f$ \text{J}/(\text{kg} \cdot \text{K}) \f$.
 * From Pruppacher and Klett (1997), equations 3-15 and 3-16
 * Validity range: \f$ 236.15 \text{K} <= \text{T} <= 308.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Specific heat capacity
 */
template <class float_t>
inline float_t specific_heat_water(float_t T) {
    float_t T_cel = T - 273.15;
    double a0, a1, a2, a3, a4;

    if (273.15 > T) {
        a0 = 1.000983;
        a1 = -2.7052e-3;
        a2 = -2.3235e-5;
        a3 = 4.3778e-6;
        a4 = 2.7136e-7;
    } else {
        a0 = 0.9979;
        a1 = 0.0;
        a2 = 3.1e-6;
        a3 = 0.0;
        a4 = 3.8e-9;
    }
    return (4186.8*polyval4(a0, a1, a2, a3, a4, T_cel));
}


/**
 * Specific heat capacity of ice in \f$ \text{J}/(\text{kg} \cdot \text{K}) \f$.
 * From Murphy and Koop (2005)
 * Validity range: \f$ 20 \text{K} <= \text{T} \f$
 *
 * @param T Temperature in Kelvin
 * @return Specific heat capacity
 */
template <class float_t>
inline float_t specific_heat_ice(float_t T, float_t M_w) {
    float_t T_frac = T/125.1;
    return ((-2.0572 + 0.14644*T + 0.06163*T*exp(-T_frac*T_frac))/M_w);
}


/**
 * Latent heat of vaporization of supercooled liquid water in
 * \f$ \text{J}/\text{kg} \f$.
 * From Murphy and Koop (2005)
 * Validity range: \f$ 236 \text{K} <= \text{T} <= 273.16 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Latent heat
 */
template <class float_t>
inline float_t latent_heat_water_supercooled(float_t T, float_t M_w) {
    return ((56579.0 - 42.212*T + exp(0.1149*(281.6-T)))/M_w);
}


/**
 * Latent heat of vaporization of liquid water in
 * \f$ \text{J}/\text{kg} \f$.
 *
 * @param T Temperature in Kelvin
 * @param L_wd Latent heat of evaporation of water (water->vapor)
 * @param cv Specific heat capacity of water vapor at constant pressure in \f$\text{J}/\text{K}/\text{kg}\f$
 * @param cp Specific heat capacity of air at constant pressure in \f$\text{J}/\text{K}/\text{kg}\f$
 * @param T_freeze Lower temperature threshold for raindrop freezing
 * @param R_v  Gas constant for water vapor, unit in \f$\text{J}/\text{K}/\text{kg}\f$
 * @param R_a Gas constant for dry air, unit: \f$\text{J}/\text{K}/\text{kg}\f$
 * @return Latent heat
 */
template <class float_t>
inline float_t latent_heat_water(
    float_t T,
    float_t L_wd,
    float_t cv,
    float_t cp,
    float_t T_freeze,
    float_t R_v,
    float_t R_a) {
    // 4.1733*cp = liquid water heat capacity
    // cp - R_a = specific heat at constant volume
    return (L_wd + (cv - 4.1733*cp) * (T - T_freeze) - R_v * T) / (cp - R_a);
}

/**
 * Latent heat of sublimation of ice in
 * \f$ \text{J}/\text{kg} \f$.
 * From Murphy and Koop (2005)
 * Validity range: \f$ 30 \text{K} <= \text{T} \f$
 *
 * @param T Temperature in Kelvin
 * @return Latent heat
 */
template <class float_t>
inline float_t latent_heat_ice(float_t T, float_t M_w) {
    float_t T_frac = T/123.75;
    return ((polyval2(46782.5, 35.8925, -0.07414, T)
        + 541.5*exp(-T_frac*T_frac))/M_w);
}


/**
 * Latent heat of melting of ice in
 * \f$ \text{J}/\text{kg} \f$.
 * From Rasmussen and Heymsfield -
 * Melting and Shedding of Graupel and Hail. Part I: Model Physics
 * Updated with ICON (instead of cal/g we use J/kg)
 * Validity range: ?
 *
 * @param T Temperature in Kelvin
 * @return Latent heat
 */
template <class float_t>
inline float_t latent_heat_melt(float_t T, float_t T_freeze) {
    // Table A1
    return 4.184e3 * (79.7+0.485*(T-T_freeze)
        - 2.5e-3*(T-T_freeze)*(T-T_freeze));
}


/**
 * Latent heat of evaporation of water in
 * \f$ \text{J}/\text{kg} \f$.
 * From Rasmussen and Heymsfield -
 * Melting and Shedding of Graupel and Hail. Part I: Model Physics
 * Updated with ICON (instead of cal/g we use J/kg)
 * Validity range: ?
 *
 * @param T Temperature in Kelvin
 * @return Latent heat
 */
template <class float_t>
inline float_t latent_heat_evap(float_t T) {
    // Table A1
    float_t lh_e0 = 2.5006e6;
    float_t gam = 0.167 + 3.67e-4 * T;
    return lh_e0 * pow(T_freeze/T, gam);
}


/**
 * Saturation vapor pressure over a flat surface of liquid water in Pa.
 * From ICON.
 * Validity range: ?
 *
 * Vanilla: From Murphy and Koop (2005).
 * Validity range: \f$ 123 \text{K} <= \text{T} <= 332 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Saturation vapor pressure
 */
template <class float_t>
inline float_t saturation_pressure_water(
    float_t T,
    float_t p_sat_low_temp,
    float_t p_sat_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_const_b) {
#ifdef VANILLA_PRESSURE
    float_t Tinv = 1.0/T;
    float_t logT = log(T);
    return (exp(54.842763 - 6763.22*Tinv - 4.21*logT + 0.000367*T
        + tanh(0.0415*(T-218.8))*(53.878 - 1331.22*Tinv - 9.44523*logT + 0.014025*T)));
#else
    return p_sat_low_temp
        * exp(p_sat_const_a * (T-T_sat_low_temp) / (T-p_sat_const_b));
#endif
}


/**
 * Saturation vapor pressure over a flat surface of ice in Pa.
 *
 * For vanilla: From Murphy and Koop (2005).
 * Validity range: \f$ 110 \text{K} <= \text{T} \f$
 *
 * @param T Temperature in Kelvin
 * @return Saturation vapor pressure
 */
template <class float_t>
inline float_t saturation_pressure_ice(
    float_t T,
    float_t p_sat_low_temp,
    float_t p_sat_ice_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_ice_const_b) {
#ifdef VANILLA_PRESSURE
    return (exp(9.550426 - (5723.265/T) + 3.53068*log(T) - 0.00728332*T));
#else
    return (p_sat_low_temp * exp(p_sat_ice_const_a
            * (T-T_sat_low_temp)
        / (T-p_sat_ice_const_b)));
#endif
}


/**
 * Surface tension of liquid water in \f$ \text{N}/\text{m} \f$.
 * From Prof. Borrmann, based on equation 5-12 from Pruppacher and Klett (1997)
 * Validity range: \f$ 233.15 \text{K} <= \text{T} <= 373.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Surface tension
 */
template <class float_t>
inline float_t surface_tension_water(float_t T) {
    float_t T_cel = T - 273.15;
    double a0, a1, a2, a3, a4, a5, a6;

    a0 = 75.7901;
    a1 = -0.139649;
    a2 = -4.62016e-4;
    a3 = -2.92323e-5;
    a4 = 1.29287e-6;
    a5 = -1.70799e-8;
    a6 = 7.25066e-11;

    return (0.001*polyval6(a0, a1, a2, a3, a4, a5, a6, T_cel));
}


/**
 * Mean free path in air in m.
 * From somewhere in Pruppacher (?)
 * Validity range: ?
 *
 * @param p Pressure in Pascal
 * @param T Temperature in Kelvin
 * @return Mean free path
 */
template <class float_t>
inline float_t mean_free_path(
    float_t &p,
    float_t &T) {

    return ((6.6e-8)*(T/293.15)*(101325.0/p));
}


/**
 * Partial pressure of water vapor in Pa.
 *
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 * @return Partial pressure
 */
template <class float_t>
inline float_t compute_pv(
    float_t T,
    float_t S,
    float_t p_sat_low_temp,
    float_t p_sat_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_const_b) {

    return (S*saturation_pressure_water(T, p_sat_low_temp, p_sat_const_a,
        T_sat_low_temp, p_sat_const_b));
}


/**
 * Partial pressure of dry air in Pa.
 *
 * @param p Total pressure in Pascal
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 * @return Partial pressure
 */
template <class float_t>
inline float_t compute_pa(
    float_t p,
    float_t T,
    float_t S,
    float_t p_sat_low_temp,
    float_t p_sat_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_const_b) {

    return (p - compute_pv(T, S, p_sat_low_temp, p_sat_const_a,
        T_sat_low_temp, p_sat_const_b));
}


/**
 * Density of dry air in \f$ \text{kg}/\text{m}^3 \f$.
 *
 * @param p Pressure in Pascal
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 * @return Density
 */
template <class float_t>
inline float_t compute_rhoa(
    float_t p,
    float_t T,
    float_t S,
    float_t p_sat_low_temp,
    float_t p_sat_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_const_b,
    float_t R_a) {

    return (compute_pa(p, T, S, p_sat_low_temp, p_sat_const_a,
        T_sat_low_temp, p_sat_const_b) / (R_a*T));
}


/**
 * Density of moist air in \f$ \text{kg}/\text{m}^3 \f$.
 *
 * @param p Pressure in Pascal
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 *
 * @return Density
 */
template <class float_t>
inline float_t compute_rhoh(
    float_t p,
    float_t T,
    float_t S,
    float_t p_sat_low_temp,
    float_t p_sat_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_const_b,
    float_t R_a,
    float_t R_v) {

    return (compute_pa(p, T, S, p_sat_low_temp, p_sat_const_a,
        T_sat_low_temp, p_sat_const_b)
        / (R_a*T)
        + compute_pv(T, S, p_sat_low_temp, p_sat_const_a,
            T_sat_low_temp, p_sat_const_b)
        / (R_v*T));
}

/**
 * Convert saturation ratio to water vapor mixing-ratio.
 *
 * @param p Total pressure in Pascal
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 * @return Water vapor mixing ratio
 */
template <class float_t>
inline float_t convert_S_to_qv(float_t p,
    float_t T,
    float_t S,
    float_t p_sat_low_temp,
    float_t p_sat_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_const_b,
    float_t Epsilon) {
// #ifdef B_EIGHT
//     float_t r_d = 287.04;
//     float_t r_v = 461.51;
//     float_t qv_sat = r_d/r_v * saturation_pressure_water(
        // T, p_sat_low_temp, p_sat_const_a, T_sat_low_temp, p_sat_const_b)
//         / (p - (1-r_d/r_v)
//             * saturation_pressure_water(
        // T, p_sat_low_temp, p_sat_const_a, T_sat_low_temp, p_sat_const_b));
//     return S*qv_sat*100.0;
// #else
    return (Epsilon*(
        compute_pv(T, S, p_sat_low_temp, p_sat_const_a, T_sat_low_temp, p_sat_const_b)
        / compute_pa(p, T, S, p_sat_low_temp, p_sat_const_a, T_sat_low_temp, p_sat_const_b)));
// #endif
}


/**
 * Convert ice saturation ratio to water vapor mixing-ratio.
 *
 * @param p Total pressure in Pascal
 * @param T Temperature in Kelvin
 * @param Si Ice saturation ratio
 * @return Water vapor mixing ratio
 */
template <class float_t>
inline float_t convert_Si_to_qv(float_t p,
    float_t T,
    float_t Si,
    float_t p_sat_low_temp,
    float_t p_sat_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_const_b,
    float_t Epsilon,
    float_t p_sat_ice_const_a,
    float_t p_sat_ice_const_b) {

    float_t S = Si * (
        saturation_pressure_ice(T, p_sat_low_temp, p_sat_ice_const_a, T_sat_low_temp, p_sat_ice_const_b)
        / saturation_pressure_water(T, p_sat_low_temp, p_sat_const_a, T_sat_low_temp, p_sat_const_b));

  return convert_S_to_qv(p, T, S, p_sat_low_temp, p_sat_const_a,
    T_sat_low_temp, p_sat_const_b, Epsilon);
}


/**
 * Convert water vapor mixing-ratio to saturation ratio.
 *
 * @param p Total pressure in Pascal
 * @param T Temperature in Kelvin
 * @param qv Water vapor mixing ratio
 * @return Saturation ratio
 */
template <class float_t>
inline float_t convert_qv_to_S(
    float_t p,
    float_t T,
    float_t qv,
    float_t p_sat_low_temp,
    float_t p_sat_const_a,
    float_t T_sat_low_temp,
    float_t p_sat_const_b,
    float_t Epsilon) {
// #ifdef B_EIGHT
//     float_t r_d = 287.04;
//     float_t r_v = 461.51;
//     float_t qv_sat = r_d/r_v * saturation_pressure_water(
        // T, p_sat_low_temp, p_sat_const_a, T_sat_low_temp, p_sat_const_b)
//         / (p - (1-r_d/r_v)
//             * saturation_pressure_water(
        // T, p_sat_low_temp, p_sat_const_a, T_sat_low_temp, p_sat_const_b));
//     return qv/qv_sat/100.0;
// #else
    return ((p*qv)/((Epsilon + qv)
        * saturation_pressure_water(T, p_sat_low_temp, p_sat_const_a, T_sat_low_temp, p_sat_const_b)));
// #endif
}


/**
 * Get mean mass of particle assuming a gamma distribution
 * (TODO: Is it though? It assumes a rather high shape parameter \f$k\f$).
 *
 * @param q Particle mixing ratio [kg/kg]
 * @param n Number of particles
 * @param min_x Minimum size of particle
 * @param max_x Maximum size of particle
 * @return Mean mass of particle in kg
 */
template <class float_t>
inline float_t particle_mean_mass(
    const float_t &q,
    const float_t &n,
    const float_t &min_x,
    const float_t &max_x) {

    return min(max(q/(n+DBL_EPSILON), min_x), max_x);
}


/**
 * Get diameter of the particle by calculating
 * \f[ a_{\text{geo}} \cdot x^{b_{\text{geo}}} \f]
 * with \f$x\f4 the particle mass in kg.
 * Seifert & Beheng; 2006 Eq. (32)
 *
 * @param x Particle mass in kg
 * @param a_geo Coefficient from particle model constant
 * @param b_geo Exponent from particle model constant
 * @return Particle diameter in m
 */
template <class float_t>
inline float_t particle_diameter(
    const float_t &x,
    const float_t &a_geo,
    const float_t &b_geo) {

    return a_geo * pow(x, b_geo);
}

/**
 * Get particle velocity of the particle by calculating
 * \f[ a_{\text{vel}} \cdot x^{b_{\text{vel}}} \f]
 * with \f$x\f4 the particle mass in kg.
 * Seifert & Beheng; 2006 Eq. (33)
 *
 * @param x Particle mass in kg
 * @param a_vel Coefficient from particle model constant
 * @param b_vel Exponent from particle model constant
 * @return Particle velocity in \f$ \text{m}/\text{s} \f$
 */
template <class float_t>
inline float_t particle_velocity(
    const float_t &x,
    const float_t &a_vel,
    const float_t &b_vel) {

    return a_vel * pow(x, b_vel);
}

/**
 * Calculate the shape parameter \f$\mu\f$ as in ICON.
 *
 * @param D_m Volume diameter
 * @param cmu0 Parametrization value
 * @param cmu1 Parametrization value
 * @param cmu2 Parametrization value
 * @param cmu3 Parametrization value
 * @param cmu4 Parametrization value
 * @return \f$mu\f$
 */
template <class float_t>
inline float_t rain_mue_dm_relation(
    const float_t &D_m,
    const float_t &cmu0,
    const float_t &cmu1,
    const float_t &cmu2,
    const float_t &cmu3,
    const float_t &cmu4) {

    float_t delta = cmu2*(D_m-cmu3);
    if (D_m <= cmu3)
        return cmu0*tanh(16.0*delta*delta) + cmu4;
    return cmu1*tanh(delta*delta) + cmu4;
}

/**
 * Calculate the shape parameter \f$\mu\f$ as in Seifert & Beheng (2008), Eq. 20.
 *
 * @param D_m Volume diameter
 * @param cmu0 Parametrization value
 * @param cmu1 Parametrization value
 * @param cmu2 Parametrization value
 * @param cmu3 Parametrization value
 * @param cmu4 Parametrization value
 * @return \f$mu\f$
 */
template <class float_t>
inline float_t rain_mue_dm_relation_sb(
    const float_t &D_m,
    const float_t &cmu0,
    const float_t &cmu1,
    const float_t &cmu2,
    const float_t &cmu3,
    const float_t &cmu4) {

    float_t delta = D_m-cmu3;
    if (D_m <= cmu3)
        return cmu0*tanh(4.0e-3*delta)*tanh(4.0*delta) + cmu4;
    return cmu1*tanh(1e-3*delta)*tanh(1e-3*delta) + cmu4;
}


/**
 * Calculate the separation diameter for wet growth of supercooled
 * water to ice.
 *
 * @param p Pressure in Pascal
 * @param T Temperature in Kelvin
 * @param qw Mixing ratio of rain and cloud water
 * @param qi Mixing ratio of ice
 * @param table Lookup-table (usually ltabdminwgg)
 * @return Separation diameter in m
 */
template <class float_t>
float_t wet_growth_diam(
    const float_t &p,
    const float_t &T,
    const float_t &qw,
    const float_t &qi,
    const table_t<float_t> &table);

/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 90.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param n \f$ k \f$ from Eq. 90
 * @return \f$\delta_b^k\f$
 */
template <class float_t = codi::RealReverse>
inline float_t coll_delta(
    particle_model_constants_t<float_t> &pc1,
    const uint64_t n) {

    float_t tmp1 = (2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)
        + get_at(pc1.constants, Particle_cons_idx::nu)+1.0+n)
        / get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp2 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp3 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp4 = (get_at(pc1.constants, Particle_cons_idx::nu)+2.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    return tgamma(tmp1)
        / tgamma(tmp2)
        * pow(tgamma(tmp3), 2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)+n)
        / pow(tgamma(tmp4), 2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)+n);
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 90.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 90
 * @return \f$\delta_b^k\f$
 */
template <class float_t = codi::RealReverse>
inline float_t coll_delta_11(
    particle_model_constants_t<float_t> &pc1,
    const uint64_t n) {

    return coll_delta<float_t>(pc1, n);
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 90.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 90
 * @return \f$\delta_a^k\f$
 */
template <class float_t = codi::RealReverse>
inline float_t coll_delta_22(
    particle_model_constants_t<float_t> &pc,
    const uint64_t n) {

    return coll_delta<float_t>(pc, n);
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 91.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 91
 * @return \f$\delta_{ba}^k\f$
 */
template <class float_t = codi::RealReverse>
inline float_t coll_delta_12(
    particle_model_constants_t<float_t> &pc1,
    particle_model_constants_t<float_t> &pc2,
    const uint64_t n) {

    float_t tmp1 = (get_at(pc1.constants, Particle_cons_idx::b_geo)+get_at(pc1.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp2 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp3 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp4 = (get_at(pc1.constants, Particle_cons_idx::nu)+2.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t res =  2.0 * tgamma(tmp1)
        / tgamma(tmp2)
        * pow(tgamma(tmp3), get_at(pc1.constants, Particle_cons_idx::b_geo))
        / pow(tgamma(tmp4), get_at(pc1.constants, Particle_cons_idx::b_geo));
    tmp1 = (get_at(pc2.constants, Particle_cons_idx::b_geo)+get_at(pc2.constants, Particle_cons_idx::nu)+1.0+n)
        / get_at(pc2.constants, Particle_cons_idx::mu);
    tmp2 = (get_at(pc2.constants, Particle_cons_idx::nu)+1.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    tmp3 = (get_at(pc2.constants, Particle_cons_idx::nu)+1.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    tmp4 = (get_at(pc2.constants, Particle_cons_idx::nu)+2.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    return res * tgamma(tmp1)
        / tgamma(tmp2)
        * pow(tgamma(tmp3), get_at(pc2.constants, Particle_cons_idx::b_geo)+n)
        / pow(tgamma(tmp4), get_at(pc2.constants, Particle_cons_idx::b_geo)+n);
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 92.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param n \f$ k \f$ from Eq. 92
 * @return \f$\vartheta_b^k\f$
 */
template <class float_t = codi::RealReverse>
inline float_t coll_theta(
    particle_model_constants_t<float_t> &pc1,
    const uint64_t n) {

    float_t tmp1 = (2.0*get_at(pc1.constants, Particle_cons_idx::b_vel)
        + 2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)
        + get_at(pc1.constants, Particle_cons_idx::nu)+1.0+n)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp2 = (2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)
        + get_at(pc1.constants, Particle_cons_idx::nu)+1.0+n)
        / get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp3 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp4 = (get_at(pc1.constants, Particle_cons_idx::nu)+2.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    return tgamma(tmp1)
        / tgamma(tmp2)
        * pow(tgamma(tmp3), 2.0*get_at(pc1.constants, Particle_cons_idx::b_vel))
        / pow(tgamma(tmp4), 2.0*get_at(pc1.constants, Particle_cons_idx::b_vel));
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 92.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 92
 * @return \f$\vartheta_b^k\f$
 */
template <class float_t = codi::RealReverse>
inline float_t coll_theta_11(
    particle_model_constants_t<float_t> &pc,
    const uint64_t n) {

    return coll_theta<float_t>(pc, n);
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 92.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 92
 * @return \f$\vartheta_a^k\f$
 */
template <class float_t = codi::RealReverse>
inline float_t coll_theta_22(
    particle_model_constants_t<float_t> &pc,
    const uint64_t n) {

    return coll_theta<float_t>(pc, n);
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 93.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 93
 * @return \f$\vartheta_{ba}^k\f$
 */
template <class float_t = codi::RealReverse>
inline float_t coll_theta_12(
    particle_model_constants_t<float_t> &pc1,
    particle_model_constants_t<float_t> &pc2,
    const uint64_t n) {

    float_t tmp1 = (get_at(pc1.constants, Particle_cons_idx::b_vel)+2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)
        + get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp2 = (2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)
        + get_at(pc1.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp3 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t tmp4 = (get_at(pc1.constants, Particle_cons_idx::nu)+2.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    float_t res = 2.0 * tgamma(tmp1)
        / tgamma(tmp2)
        * pow(tgamma(tmp3), get_at(pc1.constants, Particle_cons_idx::b_vel))
        / pow(tgamma(tmp4), get_at(pc1.constants, Particle_cons_idx::b_vel));
    tmp1 = (get_at(pc2.constants, Particle_cons_idx::b_vel)+2.0*get_at(pc2.constants, Particle_cons_idx::b_geo)
        + get_at(pc2.constants, Particle_cons_idx::nu)+1.0+n)/get_at(pc2.constants, Particle_cons_idx::mu);
    tmp2 = (2.0*get_at(pc2.constants, Particle_cons_idx::b_geo)+get_at(pc2.constants, Particle_cons_idx::nu)+1.0+n)
        / get_at(pc2.constants, Particle_cons_idx::mu);
    tmp3 = (get_at(pc2.constants, Particle_cons_idx::nu)+1.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    tmp4 = (get_at(pc2.constants, Particle_cons_idx::nu)+2.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    return res * tgamma(tmp1)
        / tgamma(tmp2)
        * pow(tgamma(tmp3), get_at(pc2.constants, Particle_cons_idx::b_vel))
        / pow(tgamma(tmp4), get_at(pc2.constants, Particle_cons_idx::b_vel));
}


/**
 * Calculate the model constant \f$a_f = \overline{a}_{\text{vent},n}\f$
 * of a particle model that is used
 * for vaporization and melting processes with formulas such as
 * \f[ a_f + b_f \sqrt{d \cdot v} \f]
 * where \f$d\f$ is the diameter of the particle to melt/vaporize and
 * \f$v\f$ is the particle velocity.
 * From Seifert & Beheng (2006), Eq. 88.
 *
 * @param pc Model constants for a particle type
 * @param n 1 for cosmo5
 * @return \f$a_f\f$
 */
template <class float_t = codi::RealReverse>
inline float_t vent_coeff_a(
    particle_model_constants_t<float_t> &pc,
    uint64_t n) {

    return get_at(pc.constants, Particle_cons_idx::a_ven)
        * tgamma((get_at(pc.constants, Particle_cons_idx::nu)+n
        + get_at(pc.constants, Particle_cons_idx::b_geo))/get_at(pc.constants, Particle_cons_idx::mu))
        / tgamma((get_at(pc.constants, Particle_cons_idx::nu)+1.0)/get_at(pc.constants, Particle_cons_idx::mu))
        * pow(tgamma((get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu)) / tgamma((get_at(pc.constants, Particle_cons_idx::nu)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu)), get_at(pc.constants, Particle_cons_idx::b_geo)+n-1.0);
}


/**
 * Calculate the model constant \f$b_f = \overline{b}_{\text{vent},n}\f$ of a particle model that is used
 * for vaporization and melting processes with formulas such as
 * \f[ a_f + b_f \sqrt{d \cdot v} \f]
 * where \f$d\f$ is the diameter of the particle to melt/vaporize and
 * \f$v\f$ is the particle velocity.
 * From Seifert & Beheng (2006), Eq. 89:
 *
 * @param pc Model constants for a particle type
 * @param n 1 for cosmo5
 * @return \f$b_f\f$
 */
template <class float_t = codi::RealReverse>
inline float_t vent_coeff_b(
    particle_model_constants_t<float_t> &pc,
    uint64_t n) {

    const double m_f = 0.5;  // From PK, page 541
    return get_at(pc.constants, Particle_cons_idx::b_ven)
        * tgamma((get_at(pc.constants, Particle_cons_idx::nu)+n+(m_f+1.0)
        * get_at(pc.constants, Particle_cons_idx::b_geo)+m_f
        * get_at(pc.constants, Particle_cons_idx::b_vel))/get_at(pc.constants, Particle_cons_idx::mu))
        / tgamma((get_at(pc.constants, Particle_cons_idx::nu)+1.0)/get_at(pc.constants, Particle_cons_idx::mu))
        * pow(tgamma((get_at(pc.constants, Particle_cons_idx::nu)+1.0)/get_at(pc.constants, Particle_cons_idx::mu))
        / tgamma((get_at(pc.constants, Particle_cons_idx::nu)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu)), (m_f+1.0) * get_at(pc.constants, Particle_cons_idx::b_geo)
        + m_f*get_at(pc.constants, Particle_cons_idx::b_vel)+n-1.0);
}


/**
 * Calculate the complete mass moment of particle size distribution from
 * Seifert & Beheng (2006), Eq. 82.
 *
 * @param pc Model constants for a particle type
 * @param n power (2 for cosmo5)
 * @return Complete mass moment \f$M^n\f$
 */
template <class float_t = codi::RealReverse>
inline float_t moment_gamma(
    particle_model_constants_t<float_t> &pc,
    uint64_t n) {
    return tgamma((n+get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu))
        / tgamma((get_at(pc.constants, Particle_cons_idx::nu)+1.0)/get_at(pc.constants, Particle_cons_idx::mu))
        * pow(tgamma((get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu))
        / tgamma((get_at(pc.constants, Particle_cons_idx::nu)+2.0)/get_at(pc.constants, Particle_cons_idx::mu)), n);
}


/**
 * Setup for bulk sedimentation velocity.
 *
 * @param pc Model constants for a certain particle type.
 */
template<class float_t = codi::RealReverse>
void setup_bulk_sedi(
    particle_model_constants_t<float_t> &pc);

/**
 * Initialize the constants for the particle collection using
 * Seifert & Beheng (2006) Eq. 90-93. This function is for all collections
 * except snow - rain and ice - rain collection.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param c Model constants for particle collection
 */
template<class float_t = codi::RealReverse>
void init_particle_collection_1(
    particle_model_constants_t<float_t> &pc1,
    particle_model_constants_t<float_t> &pc2,
    collection_model_constants_t<float_t> &c);

/**
 * Initialize the constants for the particle collection using
 * Seifert & Beheng (2006) Eq. 90-93. This function is for snow - rain and
 * ice - rain collection.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param c Model constants for particle collection
 */
template<class float_t = codi::RealReverse>
void init_particle_collection_2(
    particle_model_constants_t<float_t> &pc1,
    particle_model_constants_t<float_t> &pc2,
    collection_model_constants_t<float_t> &c);

/** @} */  // end of group parametrizations
