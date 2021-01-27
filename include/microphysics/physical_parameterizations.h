#ifndef PHYSICAL_PARAMETERIZATIONS_H
#define PHYSICAL_PARAMETERIZATIONS_H

#include "codi.hpp"
#include <cmath>
#include "types.h"
#include "constants.h"

#include <boost/math/special_functions/gamma.hpp>


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
template <class A>
inline A polyval2(const double a0,
		  const double a1,
		  const double a2,
		  A x)
{
  return ( a0 + x*(a1 + a2*x) );
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
template <class A>
inline A polyval4(const double a0,
		  const double a1,
		  const double a2,
		  const double a3,
		  const double a4,
		  A x)
{
  return ( a0 + x*(a1 + x*(a2 + x*(a3 + x*a4))) );
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
template <class A>
inline A polyval5(const double a0,
		  const double a1,
		  const double a2,
		  const double a3,
		  const double a4,
		  const double a5,
		  A x)
{
  return ( a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*a5)))) );
}

////////////////////////////////////////////////////////////////////////////////
// Evaluate a polynomial
// a_0 + a_1*x + a_2*x^2 + a_3*x^3 + a_4*x^4 + a_5*x^5 + a_6*x^6
// of order 6 using Horner's method
//
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
template <class A>
inline A polyval6(const double a0,
		  const double a1,
		  const double a2,
		  const double a3,
		  const double a4,
		  const double a5,
		  const double a6,
		  A x)
{
  return ( a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + a6*x))))) );
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
template <class A>
inline A diffusivity(A T,
		     A p)
{

  return ( (2.11e-5)*(101325.0/p)*pow(T/273.15 , 1.94) );

}


/**
 * Thermal conductivity of dry air in \f$ \text{W}/(\text{m} \cdot \text{K}) \f$.
 * From Kannuluik and Carman (1951)
 * Validity range: \f$ 90.15 \text{K} <= \text{T} <= 491.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Thermal conductivity
 */
template <class A>
inline A thermal_conductivity_dry_air(A T)
{
  // Change temperature to Celsius-scale
  A T_cel = T - 273.15;

  return ( 418.68*(5.75e-5)*polyval2(1.0, 0.00317, -0.0000021, T_cel) );
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
template <class A>
inline A thermal_conductivity_moist_air(A T, A qv, model_constants_t &cc)
{
  // Thermal conductivity dry air
  A Kt = thermal_conductivity_dry_air(T);

  A Kt_tilde = Kt/418.68;
  A Kv_tilde = (3.78e-5) + (2.0e-7)*(T - 273.15);

  return ( Kt*( 1.0 - (1.17 - 1.02*(Kt_tilde/Kv_tilde))*(qv/(qv+get_at(cc.constants, Cons_idx::Epsilon))) ) );
}


/**
 * Density of liquid water in \f$ \text{kg}/\text{m}^3 \f$.
 * From Hare and Sorensen (1987); Kell (1975)
 * Validity range: \f$ 239.74 \text{K} <= \text{T} <= 423.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Density
 */
template <class A>
inline A density_water(A T)
{

  // Change to Celsius scale
  A T_cel = T - 273.15;

  A F_T;
  A denom;
  double a0, a1, a2, a3, a4, a5, a6;

  if( 273.15 >= T ){
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

  }else{
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

  return ( F_T/denom );
}


/**
 * Density of ice in \f$ \text{kg}/\text{m}^3 \f$.
 * From Bogdan (1997)
 * Validity range: \f$ 228 \text{K} <= \text{T} <= 273.15 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Density
 */
template <class A>
inline A density_ice(A T)
{
  A T_cel = T - 273.15;

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
template <class A>
inline A specific_heat_dry_air(A T)
{

  return cp; // 1005.0

}


/**
 * Specific heat capacity of water vapor in \f$ \text{J}/(\text{kg} \cdot \text{K}) \f$.
 * From Pruppacher and Klett (1997)
 * Validity range: ?
 *
 * @param T Temperature in Kelvin
 * @return Specific heat capacity
 */
template <class A>
inline A specific_heat_water_vapor(A T)
{

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
template <class A>
inline A specific_heat_water(A T)
{
  A T_cel = T - 273.15;
  double a0, a1, a2, a3, a4;

  if( 273.15 > T ){
    a0 = 1.000983;
    a1 = -2.7052e-3;
    a2 = -2.3235e-5;
    a3 = 4.3778e-6;
    a4 = 2.7136e-7;
  }else{
    a0 = 0.9979;
    a1 = 0.0;
    a2 = 3.1e-6;
    a3 = 0.0;
    a4 = 3.8e-9;
  }

  return ( 4186.8*polyval4(a0, a1, a2, a3, a4, T_cel) );
}


/**
 * Specific heat capacity of ice in \f$ \text{J}/(\text{kg} \cdot \text{K}) \f$.
 * From Murphy and Koop (2005)
 * Validity range: \f$ 20 \text{K} <= \text{T} \f$
 *
 * @param T Temperature in Kelvin
 * @return Specific heat capacity
 */
template <class A>
inline A specific_heat_ice(A T, model_constants_t &cc)
{

  A T_frac = T/125.1;

  return ( (-2.0572 + 0.14644*T + 0.06163*T*exp( -T_frac*T_frac ))/get_at(cc.constants, Cons_idx::M_w) );

}


/**
 * Latent heat of vaporization of (supercooled) liquid water in
 * \f$ \text{J}/\text{kg} \f$.
 * From Murphy and Koop (2005)
 * Validity range: \f$ 236 \text{K} <= \text{T} <= 273.16 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Latent heat
 */
template <class A>
inline A latent_heat_water(A T, model_constants_t &cc)
{

  return ( ( 56579.0 - 42.212*T + exp( 0.1149*(281.6-T) ) )/get_at(cc.constants, Cons_idx::M_w) );

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
template <class A>
inline A latent_heat_ice(A T, model_constants_t &cc)
{

  A T_frac = T/123.75;

  return ( ( polyval2(46782.5, 35.8925, -0.07414, T) + 541.5*exp(-T_frac*T_frac) )/get_at(cc.constants, Cons_idx::M_w) );

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
template <class A>
inline A latent_heat_melt(A T, model_constants_t &cc)
{
  // Table A1
  return 4.184e3 * (79.7+0.485*(T-get_at(cc.constants, Cons_idx::T_freeze))
    - 2.5e-3*(T-get_at(cc.constants, Cons_idx::T_freeze))*(T-get_at(cc.constants, Cons_idx::T_freeze)));
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
template <class A>
inline A latent_heat_evap(A T)
{
  // Table A1
  A lh_e0 = 2.5006e6;
  A gam = 0.167 + 3.67e-4 * T;
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
template <class A>
inline A saturation_pressure_water(A T, model_constants_t &cc)
{
#ifdef VANILLA_PRESSURE
  A Tinv = 1.0/T;
  A logT = log(T);
  return ( exp( 54.842763 - 6763.22*Tinv - 4.21*logT + 0.000367*T
    + tanh(0.0415*(T-218.8))*(53.878 - 1331.22*Tinv - 9.44523*logT + 0.014025*T) ) );
#else
  return get_at(cc.constants, Cons_idx::p_sat_low_temp)
    * exp(get_at(cc.constants, Cons_idx::p_sat_const_a)
    * (T-get_at(cc.constants, Cons_idx::T_sat_low_temp))
    / (T-get_at(cc.constants, Cons_idx::p_sat_const_b)));
#endif
}


/**
 * Calculate the water vapor mixing ratio at saturation = 1 using
 * Dotzek (4.33).
 *
 * @param T Temperature in Kelvin
 * @param p Pressure in Pa
 * @return Mixing ratio at saturation = 1
 */
template <class A>
inline A water_vapor_sat_ratio_dotzek(
    A p,
    A T,
    model_constants_t &cc)
{
    A p_sat = saturation_pressure_water(T);
    return (get_at(cc.constants, Cons_idx::r_const)
        / get_at(cc.constants, Cons_idx::r1_const))
        / (p/p_sat + (get_at(cc.constants, Cons_idx::r_const)
        / get_at(cc.constants, Cons_idx::r1_const))-1);
}


/**
 * Calculate the water vapor mixing ratio at saturation = 1.
 *
 * @param T Temperature in Kelvin
 * @param p Pressure in Pa
 * @return Mixing ratio at saturation = 1
 */
template <class A>
inline A water_vapor_sat_ratio(
    A p,
    A T,
    model_constants_t &cc)
{
    A p_sat = saturation_pressure_water(T);
    return get_at(cc.constants, Cons_idx::Epsilon)*( p_sat/(p - p_sat) );
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
template <class A>
inline A saturation_pressure_ice(A T, model_constants_t &cc)
{
#ifdef VANILLA_PRESSURE
    return ( exp( 9.550426 - (5723.265/T) + 3.53068*log(T) - 0.00728332*T ) );
#else
    return ( get_at(cc.constants, Cons_idx::p_sat_low_temp)
        * exp(get_at(cc.constants, Cons_idx::p_sat_ice_const_a)
        * (T-get_at(cc.constants, Cons_idx::T_sat_low_temp))
        / (T-get_at(cc.constants, Cons_idx::p_sat_ice_const_b))) );
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
template <class A>
inline A surface_tension_water(A T)
{
  A T_cel = T - 273.15;
  double a0, a1, a2, a3, a4, a5, a6;

  a0 = 75.7901;
  a1 = -0.139649;
  a2 = -4.62016e-4;
  a3 = -2.92323e-5;
  a4 = 1.29287e-6;
  a5 = -1.70799e-8;
  a6 = 7.25066e-11;

  return ( 0.001*polyval6(a0, a1, a2, a3, a4, a5, a6, T_cel) );

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
template <class A>
inline A mean_free_path(A p,
			A T)
{

  return ( ( 6.6e-8 )*( T/293.15 )*( 101325.0/p ) );

}


/**
 * Partial pressure of water vapor in Pa.
 *
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 * @return Partial pressure
 */
template <class A>
inline A compute_pv(A T,
		    A S,
            model_constants_t &cc)
{

  return ( S*saturation_pressure_water(T, cc) );

}


/**
 * Partial pressure of dry air in Pa.
 *
 * @param p Total pressure in Pascal
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 * @return Partial pressure
 */
template <class A>
inline A compute_pa(A p,
		    A T,
		    A S,
            model_constants_t &cc)
{

  return (p - compute_pv(T,S, cc));

}


/**
 * Density of dry air in \f$ \text{kg}/\text{m}^3 \f$.
 *
 * @param p Pressure in Pascal
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 * @return Density
 */
template <class A>
inline A compute_rhoa(A p,
		      A T,
		      A S,
    model_constants_t &cc)
{

  return ( compute_pa(p,T,S, cc)/(get_at(cc.constants, Cons_idx::R_a)*T) );

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
template <class A>
inline A compute_rhoh(
    A p,
    A T,
    A S,
    model_constants_t &cc)
{
    // auto p_sat = p_sat_low_temp*exp(p_sat_const_a*(T-T_sat_low_temp)/(T-p_sat_const_b));
    // auto p_sat = exp( 9.550426 - (5723.265/T) + 3.53068*log(T) - 0.00728332*T );
    // return ( (p-S*p_sat)/(Ra*T) + S*p_sat/(R_v*T) );
    return ( compute_pa(p, T, S, cc)/(get_at(cc.constants, Cons_idx::R_a)*T)
        + compute_pv(T, S, cc)/(get_at(cc.constants, Cons_idx::R_v)*T) );

}


/**
 * Convert saturation ratio to water vapor mixing-ratio.
 *
 * @param p Total pressure in Pascal
 * @param T Temperature in Kelvin
 * @param S Saturation ratio
 * @return Water vapor mixing ratio
 */
template <class A>
inline A convert_S_to_qv(A p,
			 A T,
			 A S,
            model_constants_t &cc)
{

  return ( get_at(cc.constants, Cons_idx::Epsilon)*( compute_pv(T,S, cc)/compute_pa(p,T,S, cc) ) );

}


/**
 * Convert ice saturation ratio to water vapor mixing-ratio.
 *
 * @param p Total pressure in Pascal
 * @param T Temperature in Kelvin
 * @param Si Ice saturation ratio
 * @return Water vapor mixing ratio
 */
template <class A>
inline A convert_Si_to_qv(A p,
			  A T,
			  A Si,
            model_constants_t &cc)
{

  A S = Si * ( saturation_pressure_ice(T, cc)/saturation_pressure_water(T, cc) );

  return convert_S_to_qv(p, T, S, cc);

}


/**
 * Convert water vapor mixing-ratio to saturation ratio.
 *
 * @param p Total pressure in Pascal
 * @param T Temperature in Kelvin
 * @param qv Water vapor mixing ratio
 * @return Saturation ratio
 */
template <class A>
inline A convert_qv_to_S(A p,
			 A T,
			 A qv,
             model_constants_t &cc)
{

  return ( (p*qv)/((get_at(cc.constants, Cons_idx::Epsilon) + qv)
    * saturation_pressure_water(T, cc)) );

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
template <class A>
inline A particle_mean_mass(
    const A &q,
    const A &n,
    const A &min_x,
    const A &max_x)
{
    return min(max(q/(n+DBL_EPSILON), min_x), max_x);
}


/**
 * Get diameter of the particle by calculating
 * \f[ a_{\text{geo}} \cdot x^{b_{\text{geo}}} \f]
 * with \f$x\f4 the particle mass in kg.
 *
 * @param x Particle mass in kg
 * @param a_geo Coefficient from particle model constant
 * @param b_geo Exponent from particle model constant
 * @return Particle diameter in m
 */
template <class A>
inline A particle_diameter(
    const A &x,
    const A &a_geo,
    const A &b_geo)
{
    return a_geo * pow(x, b_geo);
}

/**
 * Get particle velocity of the particle by calculating
 * \f[ a_{\text{vel}} \cdot x^{b_{\text{vel}}} \f]
 * with \f$x\f4 the particle mass in kg.
 *
 * @param x Particle mass in kg
 * @param a_vel Coefficient from particle model constant
 * @param b_vel Exponent from particle model constant
 * @return Particle velocity in \f$ \text{m}/\text{s} \f$
 */
template <class A>
inline A particle_velocity(
    const A &x,
    const A &a_vel,
    const A &b_vel)
{
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
template <class A>
inline A rain_mue_dm_relation(
    const A &D_m,
    const A &cmu0,
    const A &cmu1,
    const A &cmu2,
    const A &cmu3,
    const A &cmu4)
{
    A delta = cmu2*(D_m-cmu3);
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
template <class A>
inline A rain_mue_dm_relation_sb(
    const A &D_m,
    const A &cmu0,
    const A &cmu1,
    const A &cmu2,
    const A &cmu3,
    const A &cmu4)
{
    A delta = D_m-cmu3;
    if (D_m <= cmu3)
        return cmu0*tanh(4.0e-3*delta)*tanh(4.0*delta) + cmu4;
    return cmu1*tanh(1e-3*delta)*tanh(1e-3*delta) + cmu4;
}


//dmin_wg_gr_ltab_equi
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
template <class A>
A wet_growth_diam(
    const A &p,
    const A &T,
    const A &qw,
    const A &qi,
    const table_t &table)
{
    A dmin_loc;
    if(T >= table.x2[table.n2-1])
    {
        dmin_loc = 0.0;
    } else if(T > table.x2[0])
    {
        dmin_loc = 999.99;
    } else
    {
        A p_lok = min( max( p, table.x1[0] ), table.x1[table.n1-1] );
        A tmp = (p_lok - table.x1[0]) * table.odx1;
        uint64_t tmp_d = floor(tmp)+1;
        uint64_t iu =  std::min(tmp_d, table.n1-1 );
        A T_lok = min( max( T, table.x2[0] ), table.x2[table.n2-1] );
        tmp = (T_lok - table.x2[0]) * table.odx2;
        tmp_d = floor(tmp)+1;
        uint64_t ju = std::min(tmp_d, table.n2-1);
        A qw_lok = min( max( qw, table.x3[0] ), table.x3[table.n3-1] );
        tmp = (qw_lok - table.x3[0]) * table.odx3;
        tmp_d = floor(tmp)+1;
        uint64_t ku = std::min(tmp_d, table.n3-1 );
        A qi_lok = min( max( qi, table.x4[0] ), table.x4[table.n4-1] );
        tmp = (qi_lok - table.x4[0]) * table.odx4;
        tmp_d = floor(tmp)+1;
        uint64_t lu = std::min(tmp_d, table.n4-1);

        std::vector<A> h1(16);
        std::vector<A> h2(8);
        std::vector<A> h3(4);
        std::vector<A> h4(2);
        // Tetra linear interpolation by Dmin
        for(uint64_t i=0; i<16; ++i)
            h1[i] = table.get(iu + i/8, ju + (i%8)/4, ku + (i%4)/2, lu+i%2);

        for(uint64_t i=0; i<8; ++i)
        {
            h2[i] = h1[i] + (h1[8+i]-h1[i]) * table.odx1 * (p_lok-table.x1[iu]);
        }
        for(uint64_t i=0; i<4; ++i)
        {
            h3[i] = h2[i] + (h2[4+i]-h2[i]) * table.odx2 * (T_lok-table.x2[iu]);
        }

        h4[0] = h3[0] + (h3[2]-h3[0])   * table.odx3 * (qw_lok-table.x3[ku]);
        h4[1] = h3[1] + (h3[3] - h3[1]) * table.odx3 * (qw_lok-table.x3[ku]);
        dmin_loc = h4[0] + (h4[1]-h4[0]) * table.odx4 * (qi_lok-table.x4[lu]);
#ifdef TRACE_GROWTH
        std::cout << "\nh4_0: " << h4[0]
                  << "\nh4_1: " << h4[1]
                  << "\ntable.odx1: " << table.odx1
                  << "\ntable.odx2: " << table.odx2
                  << "\ntable.odx3: " << table.odx3
                  << "\ntable.odx4: " << table.odx4
                  << "\nku: " << ku
                  << "\nlu: " << lu
                  << "\niu: " << iu
                  << "\nlo: " << lo
                  << "\nju: " << ju
                  << "\njo: " << jo << "\n";
#endif

    }
    return dmin_loc;
}

////////////////////////////////////////////////////////////////////////////////
// Variuous functions for collision integrals of SB2006 (Eq. 90, 91, 92, 93)

/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 90.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param n \f$ k \f$ from Eq. 90
 * @return \f$\delta_b^k\f$
 */
template <class A = codi::RealReverse>
inline A coll_delta(
    particle_model_constants_t &pc1,
    const uint64_t n)
{
    A tmp1 = (2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)+get_at(pc1.constants, Particle_cons_idx::nu)+1.0+n)
        / get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp2 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp3 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp4 = (get_at(pc1.constants, Particle_cons_idx::nu)+2.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    return tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), 2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)+n)
        / pow(tgamma( tmp4 ), 2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)+n);
}

/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 90.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 90
 * @return \f$\delta_b^k\f$
 */
template <class A = codi::RealReverse>
inline A coll_delta_11(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    return coll_delta<A>(pc1, n);
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 90.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 90
 * @return \f$\delta_a^k\f$
 */
template <class A = codi::RealReverse>
inline A coll_delta_22(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    return coll_delta<A>(pc2, n);
}

/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 91.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 91
 * @return \f$\delta_{ba}^k\f$
 */
template <class A = codi::RealReverse>
inline A coll_delta_12(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    A tmp1 = (get_at(pc1.constants, Particle_cons_idx::b_geo)+get_at(pc1.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp2 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp3 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp4 = (get_at(pc1.constants, Particle_cons_idx::nu)+2.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A res =  2.0 * tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), get_at(pc1.constants, Particle_cons_idx::b_geo))
        / pow(tgamma( tmp4 ), get_at(pc1.constants, Particle_cons_idx::b_geo));
    tmp1 = (get_at(pc2.constants, Particle_cons_idx::b_geo)+get_at(pc2.constants, Particle_cons_idx::nu)+1.0+n)
        / get_at(pc2.constants, Particle_cons_idx::mu);
    tmp2 = (get_at(pc2.constants, Particle_cons_idx::nu)+1.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    tmp3 = (get_at(pc2.constants, Particle_cons_idx::nu)+1.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    tmp4 = (get_at(pc2.constants, Particle_cons_idx::nu)+2.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    return res * tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), get_at(pc2.constants, Particle_cons_idx::b_geo)+n)
        / pow(tgamma( tmp4 ), get_at(pc2.constants, Particle_cons_idx::b_geo)+n);
}

/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 92.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param n \f$ k \f$ from Eq. 92
 * @return \f$\vartheta_b^k\f$
 */
template <class A = codi::RealReverse>
inline A coll_theta(
    particle_model_constants_t &pc1,
    const uint64_t n)
{
    A tmp1 = (2.0*get_at(pc1.constants, Particle_cons_idx::b_vel)+2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)
        + get_at(pc1.constants, Particle_cons_idx::nu)+1.0+n)/get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp2 = (2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)+get_at(pc1.constants, Particle_cons_idx::nu)+1.0+n)
        / get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp3 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp4 = (get_at(pc1.constants, Particle_cons_idx::nu)+2.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    return tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), 2.0*get_at(pc1.constants, Particle_cons_idx::b_vel))
        / pow(tgamma( tmp4 ), 2.0*get_at(pc1.constants, Particle_cons_idx::b_vel));
}


/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 92.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 92
 * @return \f$\vartheta_b^k\f$
 */
template <class A = codi::RealReverse>
inline A coll_theta_11(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    return coll_theta<A>(pc1, n);
}

/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 92.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 92
 * @return \f$\vartheta_a^k\f$
 */
template <class A = codi::RealReverse>
inline A coll_theta_22(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    return coll_theta<A>(pc2, n);
}

/**
 * Function for collision integrals of Seifert & Beheng (2006), Eq. 93.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param n \f$ k \f$ from Eq. 93
 * @return \f$\vartheta_{ba}^k\f$
 */
template <class A = codi::RealReverse>
inline A coll_theta_12(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    A tmp1 = (get_at(pc1.constants, Particle_cons_idx::b_vel)+2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)
        + get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp2 = (2.0*get_at(pc1.constants, Particle_cons_idx::b_geo)+get_at(pc1.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp3 = (get_at(pc1.constants, Particle_cons_idx::nu)+1.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A tmp4 = (get_at(pc1.constants, Particle_cons_idx::nu)+2.0)/get_at(pc1.constants, Particle_cons_idx::mu);
    A res = 2.0 * tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), get_at(pc1.constants, Particle_cons_idx::b_vel))
        / pow(tgamma( tmp4 ), get_at(pc1.constants, Particle_cons_idx::b_vel));
    tmp1 = (get_at(pc2.constants, Particle_cons_idx::b_vel)+2.0*get_at(pc2.constants, Particle_cons_idx::b_geo)
        + get_at(pc2.constants, Particle_cons_idx::nu)+1.0+n)/get_at(pc2.constants, Particle_cons_idx::mu);
    tmp2 = (2.0*get_at(pc2.constants, Particle_cons_idx::b_geo)+get_at(pc2.constants, Particle_cons_idx::nu)+1.0+n)
        / get_at(pc2.constants, Particle_cons_idx::mu);
    tmp3 = (get_at(pc2.constants, Particle_cons_idx::nu)+1.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    tmp4 = (get_at(pc2.constants, Particle_cons_idx::nu)+2.0)/get_at(pc2.constants, Particle_cons_idx::mu);
    return res * tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), get_at(pc2.constants, Particle_cons_idx::b_vel))
        / pow(tgamma( tmp4 ), get_at(pc2.constants, Particle_cons_idx::b_vel));
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
inline codi::RealReverse vent_coeff_a(
    particle_model_constants_t &pc,
    uint64_t n)
{
    return get_at(pc.constants, Particle_cons_idx::a_ven)
        * tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+n
        + get_at(pc.constants, Particle_cons_idx::b_geo))/get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)/get_at(pc.constants, Particle_cons_idx::mu) )
        * pow( tgamma((get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu)) / tgamma((get_at(pc.constants, Particle_cons_idx::nu)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu)), get_at(pc.constants, Particle_cons_idx::b_geo)+n-1.0 );
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
inline codi::RealReverse vent_coeff_b(
    particle_model_constants_t &pc,
    uint64_t n)
{
    const double m_f = 0.5; // From PK, page 541
    return get_at(pc.constants, Particle_cons_idx::b_ven)
        * tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+n+(m_f+1.0)
        * get_at(pc.constants, Particle_cons_idx::b_geo)+m_f
        * get_at(pc.constants, Particle_cons_idx::b_vel))/get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)/get_at(pc.constants, Particle_cons_idx::mu) )
        * pow( tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)/get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+2.0)
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
inline codi::RealReverse moment_gamma(
    particle_model_constants_t &pc,
    uint64_t n)
{
    return tgamma( (n+get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)/get_at(pc.constants, Particle_cons_idx::mu) )
        * pow(tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+2.0)/get_at(pc.constants, Particle_cons_idx::mu)), n);
}


/**
 * Set the constants for the cloud model from given environmental conditions.
 *
 * @param y Vector of initial conditions for pressure and temperature
 * @param cc Pointer to constants from the model. On out: modified constants
 * @param ref Pointer to reference quantities to transform between units
 */
void setCoefficients(
    std::vector<codi::RealReverse> & y,
    model_constants_t& cc,
    reference_quantities_t& ref)
{
    codi::RealReverse p_prime = y[p_idx]*ref.pref;
    codi::RealReverse T_prime = y[T_idx]*ref.Tref;

    codi::RealReverse rho_prime = p_prime /( get_at(cc.constants, Cons_idx::R_a) * T_prime );
    codi::RealReverse L_vap_prime = latent_heat_water(T_prime, cc);
    codi::RealReverse Ka_prime = thermal_conductivity_dry_air(T_prime);
    codi::RealReverse psat_prime = saturation_pressure_water(T_prime, cc);
    codi::RealReverse A_pp = (L_vap_prime/(Ka_prime*T_prime))
        * ((L_vap_prime/(get_at(cc.constants, Cons_idx::R_v)*T_prime)) - 1.0);
    codi::RealReverse B_pp = (get_at(cc.constants, Cons_idx::R_v)*T_prime)/((2.21/p_prime)*psat_prime);


    cc.constants[static_cast<int>(Cons_idx::e1_prime)] = cc.e1_scale
        * ( pow(rho_prime, 2.0*cc.alpha_r-2.0)/(A_pp + B_pp) );
    cc.constants[static_cast<int>(Cons_idx::e2_prime)] = cc.e2_scale
        * ( pow(rho_prime, cc.alpha_r*cc.epsilonr - (7.0/4.0))/(A_pp + B_pp) );

    cc.constants[static_cast<int>(Cons_idx::a1_prime)] = cc.a1_scale;	// Constant coefficient
    cc.constants[static_cast<int>(Cons_idx::a2_prime)] = cc.a2_scale;	// Constant coefficient
    cc.constants[static_cast<int>(Cons_idx::d_prime)] = cc.d_scale;	// Constant coefficient
}

/**
 * Set the constants for the cloud model from given environmental conditions.
 *
 * @param p_prime Initial pressure in Pa
 * @param T_prime Initial temperature in Kelvin
 * @param cc Pointer to constants from the model. On out: modified constants
 */
void setCoefficients(
    codi::RealReverse p_prime,
    codi::RealReverse T_prime,
    model_constants_t &cc)
{
    codi::RealReverse rho_prime = p_prime /( get_at(cc.constants, Cons_idx::R_a) * T_prime );
    codi::RealReverse L_vap_prime = latent_heat_water(T_prime, cc);
    codi::RealReverse Ka_prime = thermal_conductivity_dry_air(T_prime);
    codi::RealReverse psat_prime = saturation_pressure_water(T_prime, cc);
    codi::RealReverse A_pp = (L_vap_prime/(Ka_prime*T_prime))*((L_vap_prime/(get_at(cc.constants, Cons_idx::R_v)*T_prime)) - 1.0);
    codi::RealReverse B_pp = (get_at(cc.constants, Cons_idx::R_v)*T_prime)/((2.21/p_prime)*psat_prime);

    cc.constants[static_cast<int>(Cons_idx::a1_prime)] = cc.a1_scale;	// Constant coefficient
    cc.constants[static_cast<int>(Cons_idx::a2_prime)] = cc.a2_scale;	// Constant coefficient

    cc.constants[static_cast<int>(Cons_idx::e1_prime)] = cc.e1_scale * ( pow(rho_prime, 2.0*cc.alpha_r-2.0)/(A_pp + B_pp) );
    cc.constants[static_cast<int>(Cons_idx::e2_prime)] = cc.e2_scale * ( pow(rho_prime, cc.alpha_r*cc.epsilonr - (7.0/4.0))/(A_pp + B_pp) );

    cc.constants[static_cast<int>(Cons_idx::d_prime)] = cc.d_scale;	// Constant coefficient
}



/**
 * Setup the cloud autoconversion parameters.
 *
 * @param pc Model constants for a certain particle type.
 */
void setup_cloud_autoconversion(
    particle_model_constants_t &pc,
    model_constants_t &cc)
{
    auto nu = get_at(pc.constants, Particle_cons_idx::nu) + 1.0;
    auto mu = get_at(pc.constants, Particle_cons_idx::mu);
    if(get_at(pc.constants, Particle_cons_idx::mu) == 1.0)
    {
        cc.constants[static_cast<int>(Cons_idx::cloud_k_au)] =
            get_at(cc.constants, Cons_idx::kc_autocon)
            / get_at(pc.constants, Particle_cons_idx::max_x) * 0.05
            * (nu+1.0)*(nu+3.0) / pow(nu, 2);
        cc.constants[static_cast<int>(Cons_idx::cloud_k_sc)] =
            get_at(cc.constants, Cons_idx::kc_autocon) * (nu+1.0)/(nu);
    } else
    {
        cc.constants[static_cast<int>(Cons_idx::cloud_k_au)] =
            get_at(cc.constants, Cons_idx::kc_autocon)
            / get_at(pc.constants, Particle_cons_idx::max_x) * 0.05
            * (2.0 * tgamma((nu+3.0)/mu)
            * tgamma((nu+1.0)/mu) * pow(tgamma((nu)/mu), 2)
            - 1.0 * pow(tgamma((nu+2.0)/mu), 2) * pow(tgamma((nu)/mu), 2))
            / pow(tgamma((nu+1.0)/mu), 4);
        cc.constants[static_cast<int>(Cons_idx::cloud_k_sc)] =
            get_at(cc.constants, Cons_idx::kc_autocon)
            * get_at(pc.constants, Particle_cons_idx::c_z);
    }
}

/**
 * Setup for bulk sedimentation velocity.
 *
 * @param pc Model constants for a certain particle type.
 */
void setup_bulk_sedi(
    particle_model_constants_t &pc)
{
    get_at(pc.constants, Particle_cons_idx::alfa_n) =
        get_at(pc.constants, Particle_cons_idx::a_vel)
        * tgamma( (get_at(pc.constants, Particle_cons_idx::nu)
        + get_at(pc.constants, Particle_cons_idx::b_vel)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu));
    get_at(pc.constants, Particle_cons_idx::alfa_q) =
        get_at(pc.constants, Particle_cons_idx::a_vel)
        * tgamma( (get_at(pc.constants, Particle_cons_idx::nu)
        + get_at(pc.constants, Particle_cons_idx::b_vel)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu) );
    get_at(pc.constants, Particle_cons_idx::lambda) =
        tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu) );
}

/**
 * Initialize the constants for the particle collection using
 * Seifert & Beheng (2006) Eq. 90-93. This function is for all collections
 * except snow - rain and ice - rain collection.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param c Model constants for particle collection
 */
void init_particle_collection_1(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    collection_model_constants_t &c)
{
    c.delta_n_aa = coll_delta_11(pc1, pc2, 0);
    c.delta_n_ab = coll_delta_12(pc1, pc2, 0);
    c.delta_n_bb = coll_delta_22(pc1, pc2, 0);
    c.delta_q_aa = coll_delta_11(pc1, pc2, 0);
    c.delta_q_ab = coll_delta_12(pc1, pc2, 1);
    c.delta_q_bb = coll_delta_22(pc1, pc2, 1);

    c.theta_n_aa = coll_theta_11(pc1, pc2, 0);
    c.theta_n_ab = coll_theta_12(pc1, pc2, 0);
    c.theta_n_bb = coll_theta_22(pc1, pc2, 0);
    c.theta_q_aa = coll_theta_11(pc1, pc2, 0);
    c.theta_q_ab = coll_theta_12(pc1, pc2, 1);
    c.theta_q_bb = coll_theta_22(pc1, pc2, 1);
}


/**
 * Initialize the constants for the particle collection using
 * Seifert & Beheng (2006) Eq. 90-93. This function is for snow - rain and
 * ice - rain collection.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param c Model constants for particle collection
 */
void init_particle_collection_2(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    collection_model_constants_t &c)
{
    c.delta_n_aa = coll_delta_11(pc1, pc2, 0);
    c.delta_n_ab = coll_delta_12(pc1, pc2, 0);
    c.delta_n_bb = coll_delta_22(pc1, pc2, 0);
    c.delta_q_aa = coll_delta_11(pc1, pc2, 1);
    c.delta_q_ab = coll_delta_12(pc1, pc2, 1);
    c.delta_q_ba = coll_delta_12(pc2, pc1, 1);
    c.delta_q_bb = coll_delta_22(pc1, pc2, 1);

    c.theta_n_aa = coll_theta_11(pc1, pc2, 0);
    c.theta_n_ab = coll_theta_12(pc1, pc2, 0);
    c.theta_n_bb = coll_theta_22(pc1, pc2, 0);
    c.theta_q_aa = coll_theta_11(pc1, pc2, 1);
    c.theta_q_ab = coll_theta_12(pc1, pc2, 1);
    c.theta_q_ba = coll_theta_12(pc2, pc1, 1);
    c.theta_q_bb = coll_theta_22(pc1, pc2, 1);
}


/**
 * Setup model constants and gamma tables.
 */
void setup_model_constants(
      const input_parameters_t &input,
      model_constants_t &cc,
      const reference_quantities_t &ref_quant)
{
    cc.id = std::to_string(input.id);
    // Set constants
    cc.constants[static_cast<int>(Cons_idx::q_crit_i)] = q_crit_i;
    cc.constants[static_cast<int>(Cons_idx::D_crit_i)] = D_crit_i;
    cc.constants[static_cast<int>(Cons_idx::D_conv_i)] = D_conv_i;
    cc.constants[static_cast<int>(Cons_idx::q_crit_r)] = q_crit_r;
    cc.constants[static_cast<int>(Cons_idx::D_crit_r)] = D_crit_r;
    cc.constants[static_cast<int>(Cons_idx::q_crit_fr)] = q_crit_fr;
    cc.constants[static_cast<int>(Cons_idx::D_coll_c)] = D_coll_c;
    cc.constants[static_cast<int>(Cons_idx::q_crit)] = q_crit;
    cc.constants[static_cast<int>(Cons_idx::D_conv_sg)] = D_conv_sg;
    cc.constants[static_cast<int>(Cons_idx::D_conv_ig)] = D_conv_ig;
    cc.constants[static_cast<int>(Cons_idx::x_conv)] = x_conv;
    cc.constants[static_cast<int>(Cons_idx::parcel_height)] = parcel_height;
    cc.constants[static_cast<int>(Cons_idx::inv_z)] = 1.0/parcel_height;
    cc.constants[static_cast<int>(Cons_idx::alpha_spacefilling)] = alpha_spacefilling;
    cc.constants[static_cast<int>(Cons_idx::T_nuc)] = T_nuc;
    cc.constants[static_cast<int>(Cons_idx::T_freeze)] = T_freeze;
    cc.constants[static_cast<int>(Cons_idx::T_f)] = T_f;
    cc.constants[static_cast<int>(Cons_idx::D_eq)] = D_eq;
    cc.constants[static_cast<int>(Cons_idx::rho_w)] = rho_w;
    cc.constants[static_cast<int>(Cons_idx::rho_0)] = rho_0;
    cc.constants[static_cast<int>(Cons_idx::rho_vel)] = rho_vel;
    cc.constants[static_cast<int>(Cons_idx::rho_vel_c)] = rho_vel_c;
    cc.constants[static_cast<int>(Cons_idx::rho_ice)] = rho_ice;
    cc.constants[static_cast<int>(Cons_idx::M_w)] = M_w;
    cc.constants[static_cast<int>(Cons_idx::M_a)] = M_a;
    cc.constants[static_cast<int>(Cons_idx::R_universal)] = R_universal;
    cc.constants[static_cast<int>(Cons_idx::Epsilon)] = Epsilon;
    cc.constants[static_cast<int>(Cons_idx::gravity_acc)] = gravity_acc;
    cc.constants[static_cast<int>(Cons_idx::R_a)] = R_a;
    cc.constants[static_cast<int>(Cons_idx::R_v)] = R_v;
    cc.constants[static_cast<int>(Cons_idx::a_v)] = a_v;
    cc.constants[static_cast<int>(Cons_idx::b_v)] = b_v;
    cc.constants[static_cast<int>(Cons_idx::a_prime)] = a_prime;
    cc.constants[static_cast<int>(Cons_idx::b_prime)] = b_prime;
    cc.constants[static_cast<int>(Cons_idx::c_prime)] = c_prime;
    cc.constants[static_cast<int>(Cons_idx::K_T)] = K_T;
    cc.constants[static_cast<int>(Cons_idx::L_wd)] = L_wd;
    cc.constants[static_cast<int>(Cons_idx::L_ed)] = L_ed;
    cc.constants[static_cast<int>(Cons_idx::D_v)] = D_v;
    cc.constants[static_cast<int>(Cons_idx::ecoll_min)] = ecoll_min;
    cc.constants[static_cast<int>(Cons_idx::ecoll_gg)] = ecoll_gg;
    cc.constants[static_cast<int>(Cons_idx::ecoll_gg_wet)] = ecoll_gg_wet;
    cc.constants[static_cast<int>(Cons_idx::kin_visc_air)] = kin_visc_air;
    cc.constants[static_cast<int>(Cons_idx::C_mult)] = C_mult;
    cc.constants[static_cast<int>(Cons_idx::T_mult_min)] = T_mult_min;
    cc.constants[static_cast<int>(Cons_idx::T_mult_max)] = T_mult_max;
    cc.constants[static_cast<int>(Cons_idx::T_mult_opt)] = T_mult_opt;
    cc.constants[static_cast<int>(Cons_idx::kc_autocon)] = kc_autocon;
    cc.constants[static_cast<int>(Cons_idx::D_rainfrz_gh)] = D_rainfrz_gh;
    cc.constants[static_cast<int>(Cons_idx::D_rainfrz_ig)] = D_rainfrz_ig;
    cc.constants[static_cast<int>(Cons_idx::dv0)] = dv0;
    cc.constants[static_cast<int>(Cons_idx::p_sat_melt)] = p_sat_melt;
    cc.constants[static_cast<int>(Cons_idx::cp)] = cp;
    cc.constants[static_cast<int>(Cons_idx::k_b)] = k_b;
    cc.constants[static_cast<int>(Cons_idx::a_HET)] = a_HET;
    cc.constants[static_cast<int>(Cons_idx::b_HET)] = b_HET;
    cc.constants[static_cast<int>(Cons_idx::N_sc)] = N_sc;
    cc.constants[static_cast<int>(Cons_idx::n_f)] = n_f;
    cc.constants[static_cast<int>(Cons_idx::N_avo)] = N_avo;
    cc.constants[static_cast<int>(Cons_idx::na_dust)] = na_dust;
    cc.constants[static_cast<int>(Cons_idx::na_soot)] = na_soot;
    cc.constants[static_cast<int>(Cons_idx::na_orga)] = na_orga;
    cc.constants[static_cast<int>(Cons_idx::ni_het_max)] = ni_het_max;
    cc.constants[static_cast<int>(Cons_idx::ni_hom_max)] = ni_hom_max;
    cc.constants[static_cast<int>(Cons_idx::a_dep)] = a_dep;
    cc.constants[static_cast<int>(Cons_idx::b_dep)] = b_dep;
    cc.constants[static_cast<int>(Cons_idx::c_dep)] = c_dep;
    cc.constants[static_cast<int>(Cons_idx::d_dep)] = d_dep;
    cc.constants[static_cast<int>(Cons_idx::nim_imm)] = nim_imm;
    cc.constants[static_cast<int>(Cons_idx::nin_dep)] = nin_dep;
    cc.constants[static_cast<int>(Cons_idx::alf_imm)] = alf_imm;
    cc.constants[static_cast<int>(Cons_idx::bet_dep)] = bet_dep;
    cc.constants[static_cast<int>(Cons_idx::bet_imm)] = bet_imm;
    cc.constants[static_cast<int>(Cons_idx::r_const)] = r_const;
    cc.constants[static_cast<int>(Cons_idx::r1_const)] = r1_const;
    cc.constants[static_cast<int>(Cons_idx::cv)] = cv;
    cc.constants[static_cast<int>(Cons_idx::p_sat_const_a)] = p_sat_const_a;
    cc.constants[static_cast<int>(Cons_idx::p_sat_ice_const_a)] = p_sat_ice_const_a;
    cc.constants[static_cast<int>(Cons_idx::p_sat_const_b)] = p_sat_const_b;
    cc.constants[static_cast<int>(Cons_idx::p_sat_ice_const_b)] = p_sat_ice_const_b;
    cc.constants[static_cast<int>(Cons_idx::p_sat_low_temp)] = p_sat_low_temp;
    cc.constants[static_cast<int>(Cons_idx::T_sat_low_temp)] = T_sat_low_temp;
    cc.constants[static_cast<int>(Cons_idx::alpha_depo)] = alpha_depo;
    cc.constants[static_cast<int>(Cons_idx::r_0)] = r_0;

    cc.constants[static_cast<int>(Cons_idx::k_1_conv)] = k_1_conv;
    cc.constants[static_cast<int>(Cons_idx::k_2_conv)] = k_2_conv;
    cc.constants[static_cast<int>(Cons_idx::k_1_accr)] = k_1_accr;
    cc.constants[static_cast<int>(Cons_idx::k_r)] = k_r;


    // Numerics
#ifdef MET3D
    cc.t_end_prime = input.t_end_prime + input.start_time;;
#else
    cc.t_end_prime = input.t_end_prime;
#endif

    cc.t_end = cc.t_end_prime/ref_quant.tref;
    // Time of the substeps
    cc.dt = input.dt_prime/ref_quant.tref;
    cc.dt_prime = input.dt_prime;
    cc.dt_traject_prime = 20.0;
    cc.dt_traject = cc.dt_traject_prime/ref_quant.tref;
    // The trajectories are calculated with 20 s timesteps.
    cc.num_sub_steps = (floor( cc.dt_traject_prime/cc.dt ) < 1) ? 1 : floor( cc.dt_traject_prime/cc.dt );

    // Evaluate the general performance constants
    cc.dt_half = cc.dt*0.5;
    cc.dt_third = cc.dt/3.0;
    cc.dt_sixth = cc.dt/6.0;

    // Scaling factor from input
    cc.scaling_fact = input.scaling_fact;

    // Accomodation coefficient
    cc.alpha_d = 1.0;

    // // Performance constants for warm cloud; COSMO
    cc.a1_scale = 1.0e-3;
    cc.a2_scale = 1.72 / pow(get_at(cc.constants, Cons_idx::R_a) , 7./8.);
    cc.e1_scale = 1.0 / sqrt(get_at(cc.constants, Cons_idx::R_a));
    cc.e2_scale = 9.1 / pow(get_at(cc.constants, Cons_idx::R_a) , 11./16.);
    cc.d_scale = ( 130.0*tgamma(4.5) )/( 6.0*(1.0e3)*pow(M_PI*(8.0e6)*get_at(cc.constants, Cons_idx::R_a) , 1.0/8.0) );

    // Performance constants for warm cloud; IFS
    // The file constants.h also defines some constants as nar, ...
#ifndef MET3D
    const double Nc = 50; 	// 50 over ocean; 300 over land
    const double F_aut = 1.5;
    const double F_acc = 2.0;
    const double lambda_pp = pow(cc.nar * cc.ar * tgamma(cc.br + 1.0) , cc.alpha_r);
#endif
    for(uint32_t i=0; i<4; i++)
    {
        cc.constants[static_cast<int>(Cons_idx::a_ccn_1)+i] = a_ccn[i];
        cc.constants[static_cast<int>(Cons_idx::b_ccn_1)+i] = b_ccn[i];
        cc.constants[static_cast<int>(Cons_idx::c_ccn_1)+i] = c_ccn[i];
        cc.constants[static_cast<int>(Cons_idx::d_ccn_1)+i] = d_ccn[i];
    }

    // cc.a1_scale = (1350. * F_aut)/pow(Nc , 1.79);
    // cc.a2_scale = 67.0 * F_acc;
    // cc.e1_scale = 2.0 * M_PI * cc.nar * ( (0.78 * tgamma(2.0 - cc.nbr))/(lambda_pp*lambda_pp) );
    // cc.e2_scale = cc.scaling_fact * 2.0 * M_PI * cc.nar * 0.31
    //     * pow(cc.cr/get_at(cc.constants, Particle_cons_idx::mu), 0.5) * pow(cc.Sc, 1.0/3.0) * pow(cc.rho0, 0.25)
    //     * (tgamma(cc.epsilonr + cc.nbr)/pow(lambda_pp ,cc.epsilonr));
    // cc.d_scale = 4.0e-3;

    if(nuc_type == 6)
    {
        cc.constants[static_cast<int>(Cons_idx::na_dust)] = na_dust;
        cc.constants[static_cast<int>(Cons_idx::na_soot)] = na_soot;
        cc.constants[static_cast<int>(Cons_idx::na_orga)] = na_orga;
    } else if(nuc_type == 7 || nuc_type == 5)
    {
        // Standard values
        cc.constants[static_cast<int>(Cons_idx::na_dust)] = na_dust_2;
        cc.constants[static_cast<int>(Cons_idx::na_soot)] = na_soot_2;
        cc.constants[static_cast<int>(Cons_idx::na_orga)] = na_orga_2;
    } else if(nuc_type == 8)
    {
        cc.constants[static_cast<int>(Cons_idx::na_dust)] = na_dust_3;
        cc.constants[static_cast<int>(Cons_idx::na_soot)] = na_soot_3;
        cc.constants[static_cast<int>(Cons_idx::na_orga)] = na_orga_3;
    }

    // Inflow from above
    // cc.B_prime = 0.0; //1.0e-7;

    // // Exponents of the cloud model
    // // COSMO
    cc.constants[static_cast<int>(Cons_idx::gamma)] = 1.0;
    cc.constants[static_cast<int>(Cons_idx::betac)] = 1.0;
    cc.constants[static_cast<int>(Cons_idx::betar)] = 7./8.;
    cc.constants[static_cast<int>(Cons_idx::delta1)] = 0.5;
    cc.constants[static_cast<int>(Cons_idx::delta2)] = 11./16.;
    cc.constants[static_cast<int>(Cons_idx::zeta)] = 9./8.;
#if defined(RK4ICE) || defined(RK4NOICE)
    // Exponents of the cloud model
    // IFS
    // get_at(cc.constants, Particle_cons_idx::gamma) = 2.47;
    // get_at(cc.constants, Cons_idx::betac) = 1.15;
    // get_at(cc.constants, Cons_idx::betar) = 1.15;
    // get_at(cc.constants, Cons_idx::delta1) = 2.0/( cc.br + 1.0 - cc.nbr );
    // get_at(cc.constants, Cons_idx::delta2) = ( 0.5*cc.dr + 2.5 - cc.nbr )/( cc.br + 1.0 - cc.nbr );
    // get_at(cc.constants, Cons_idx::zeta) = 1.0;

    // ==================================================
    // Set rain constants
    // See init_2mom_scheme_once in mo_2mom_mcrph_main.f90
    // ==================================================
    // Cosmo5 although nue1nue1 would be okay too, I guess
    //// Cloud
#ifdef SB_SHAPE
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::nu)] = 1;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::mu)] = 1;
#else
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::nu)] = 0;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0;
#endif
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::max_x)] = 2.6e-10;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 4.2e-15;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 1.24e-1;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.333333;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 3.75e5;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.666667;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 1.0;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.0;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-6;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 1.0e-5;
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0
        / get_at(cc.cloud.constants, Particle_cons_idx::cap);
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(cc.cloud, 1);
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(cc.cloud, 1)
        * pow(get_at(cc.constants, Cons_idx::N_sc), get_at(cc.constants, Cons_idx::n_f))
        / sqrt(get_at(cc.constants, Cons_idx::kin_visc_air));
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(cc.cloud, 2);

    setup_cloud_autoconversion(cc.cloud, cc);
    setup_bulk_sedi(cc.cloud);

    //// Rain
#ifdef SB_SHAPE
    cc.rain.constants[static_cast<int>(Particle_cons_idx::nu)] = -2.0/3.0; // SB: -2/3 COSMO: 0.0
#else
     cc.rain.constants[static_cast<int>(Particle_cons_idx::nu)] = 0; // SB: -2/3 COSMO: 0.0
#endif
    cc.rain.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0; // SB: 1/3 COMSO: 1.0/3.0
    cc.rain.constants[static_cast<int>(Particle_cons_idx::max_x)] = 3.0e-6;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 2.6e-10;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 1.24e-1;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.333333;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 114.0137;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.234370;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;

    // From rainSBBcoeffs
    cc.rain.constants[static_cast<int>(Particle_cons_idx::alpha)] = 9.292;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::beta)] = 9.623;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::gamma)] = 6.222e2;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::cmu0)] = 6.0;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::cmu1)] = 3.0e1;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::cmu2)] = 1.0e3;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::cmu3)] = 1.1e-3;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::cmu4)] = 1.0;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::cmu5)] = 2.0;
    cc.constants[static_cast<int>(Cons_idx::rain_gfak)] = 1.0;

    cc.rain.constants[static_cast<int>(Particle_cons_idx::nm1)] =
        (get_at(cc.rain.constants, Particle_cons_idx::nu)+1.0)/get_at(cc.rain.constants, Particle_cons_idx::mu);
    cc.rain.constants[static_cast<int>(Particle_cons_idx::nm2)] =
        (get_at(cc.rain.constants, Particle_cons_idx::nu)+2.0)/get_at(cc.rain.constants, Particle_cons_idx::mu);
    cc.rain.constants[static_cast<int>(Particle_cons_idx::nm3)] =
        (get_at(cc.rain.constants, Particle_cons_idx::nu)+3.0)/get_at(cc.rain.constants, Particle_cons_idx::mu);
    cc.rain.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 20.0;
    cc.rain.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.1;
    cc.table_r1.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(cc.rain.constants, Particle_cons_idx::nm1).getValue());
    cc.table_r2.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(cc.rain.constants, Particle_cons_idx::nm2).getValue());
    cc.table_r3.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(cc.rain.constants, Particle_cons_idx::nm3).getValue());
    cc.rain.constants[static_cast<int>(Particle_cons_idx::g1)] = cc.table_r1.igf[cc.table_r1.n_bins-1];
    cc.rain.constants[static_cast<int>(Particle_cons_idx::g2)] = cc.table_r2.igf[cc.table_r2.n_bins-1];
    cc.rain.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(cc.rain.constants, Particle_cons_idx::cap);
    cc.rain.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(cc.rain, 1);
    cc.rain.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(cc.rain, 1)
        * pow(get_at(cc.constants, Cons_idx::N_sc), get_at(cc.constants, Cons_idx::n_f))
        / sqrt(get_at(cc.constants, Cons_idx::kin_visc_air));
    cc.rain.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(cc.rain, 2);
    setup_bulk_sedi(cc.rain);

    //// Graupel
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::nu)] = 1.0; // SB
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0; // SB
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::max_x)] = 5.0e-4;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 1.0e-9;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 1.42e-1;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.314;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 86.89371;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.268325;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 30.0;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.1;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 100.0e-6;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-6;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::s_vel)] = 0.0;

    cc.graupel.constants[static_cast<int>(Particle_cons_idx::nm1)] =
        (get_at(cc.graupel.constants, Particle_cons_idx::nu)+1.0)/get_at(cc.graupel.constants, Particle_cons_idx::mu);
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::nm2)] =
        (get_at(cc.graupel.constants, Particle_cons_idx::nu)+2.0)/get_at(cc.graupel.constants, Particle_cons_idx::mu);
    codi::RealReverse a = (get_at(cc.graupel.constants, Particle_cons_idx::nu)+1.0)/get_at(cc.graupel.constants, Particle_cons_idx::mu);
    cc.table_g1.init_gamma_table(n_lookup, n_lookup_hr_dummy, get_at(cc.graupel.constants, Particle_cons_idx::nm1).getValue());
    a = (get_at(cc.graupel.constants, Particle_cons_idx::nu)+2.0)/get_at(cc.graupel.constants, Particle_cons_idx::mu);
    cc.table_g2.init_gamma_table(n_lookup, n_lookup_hr_dummy, get_at(cc.graupel.constants, Particle_cons_idx::nm2).getValue());
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::g1)] = cc.table_g1.igf[cc.table_g1.n_bins-1];
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::g2)] = cc.table_g2.igf[cc.table_g2.n_bins-1];
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(cc.graupel.constants, Particle_cons_idx::cap);
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = 1.0;
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(cc.graupel, 1);
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(cc.graupel, 1)
        * pow(get_at(cc.constants, Cons_idx::N_sc), get_at(cc.constants, Cons_idx::n_f))
        / sqrt(get_at(cc.constants, Cons_idx::kin_visc_air));
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(cc.graupel, 2);
    setup_bulk_sedi(cc.graupel);

    //// Hail
    cc.hail.constants[static_cast<int>(Particle_cons_idx::nu)] = 1.0;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::max_x)] = 5.0e-4;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 2.6e-9;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 0.1366;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 1.0/3.0;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 39.3;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.166667;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 30.0;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.1;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = 1.0;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 100.0e-6;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-6;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::s_vel)] = 0.0;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(cc.hail.constants, Particle_cons_idx::cap);
    cc.hail.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = 1.0;
    cc.hail.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(cc.hail, 1);
    cc.hail.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(cc.hail, 1)
        * pow(get_at(cc.constants, Cons_idx::N_sc), get_at(cc.constants, Cons_idx::n_f))
        / sqrt(get_at(cc.constants, Cons_idx::kin_visc_air));
    cc.hail.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(cc.hail, 2);
    setup_bulk_sedi(cc.hail);

    //// Ice
#ifdef SB_SHAPE
    cc.ice.constants[static_cast<int>(Particle_cons_idx::nu)] = 1.0;
#else
    cc.ice.constants[static_cast<int>(Particle_cons_idx::nu)] = 0.0; // COSMO 0.0, SB: 1.0
#endif
    cc.ice.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::max_x)] = 1.0e-5;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 1.0e-12;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 0.835;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.39;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 2.77e1;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.21579;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 3.0;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.0;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = 0.8;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 150.0e-6;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-5;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::s_vel)] = 0.05;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(cc.ice.constants, Particle_cons_idx::cap);
    cc.ice.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = 0.80;
    cc.ice.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(cc.ice, 1);
    cc.ice.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(cc.ice, 1)
        * pow(get_at(cc.constants, Cons_idx::N_sc), get_at(cc.constants, Cons_idx::n_f))
        / sqrt(get_at(cc.constants, Cons_idx::kin_visc_air));
    cc.ice.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(cc.ice, 2);
    setup_bulk_sedi(cc.ice);

    //// Snow
#ifdef SB_SHAPE
    cc.snow.constants[static_cast<int>(Particle_cons_idx::nu)] = 1.0; // COSMO: 0.0, SB 1.0
    cc.snow.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0; // COSMO 0.5, SB: 1.0/3.0
#else
    cc.snow.constants[static_cast<int>(Particle_cons_idx::nu)] = 0.0; // COSMO: 0.0, SB 1.0
    cc.snow.constants[static_cast<int>(Particle_cons_idx::mu)] = 0.5; // COSMO 0.5, SB: 1.0/3.0
#endif
    cc.snow.constants[static_cast<int>(Particle_cons_idx::max_x)] = 2.0e-5;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 1.0e-10;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 2.4;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.455;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 8.8;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.15;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 3.0;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.1;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = 0.8;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 150.0e-6;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-5;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::s_vel)] = 0.25;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(cc.snow.constants, Particle_cons_idx::cap);
    cc.snow.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = 0.80;
    cc.snow.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(cc.snow, 1);
    cc.snow.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(cc.snow, 1)
        * pow(get_at(cc.constants, Cons_idx::N_sc), get_at(cc.constants, Cons_idx::n_f))
        / sqrt(get_at(cc.constants, Cons_idx::kin_visc_air));
    cc.snow.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(cc.snow, 2);
    setup_bulk_sedi(cc.snow);

    cc.constants[static_cast<int>(Cons_idx::const0)] =
        1.0/(get_at(cc.constants, Cons_idx::D_coll_c)
        - get_at(cc.cloud.constants, Particle_cons_idx::d_crit_c));
    cc.constants[static_cast<int>(Cons_idx::const3)] =
        1.0/(get_at(cc.constants, Cons_idx::T_mult_opt) - get_at(cc.constants, Cons_idx::T_mult_min));
    cc.constants[static_cast<int>(Cons_idx::const4)] =
        1.0/(get_at(cc.constants, Cons_idx::T_mult_opt) - get_at(cc.constants, Cons_idx::T_mult_max));
    cc.constants[static_cast<int>(Cons_idx::const5)] =
        get_at(cc.constants, Cons_idx::alpha_spacefilling)
        * get_at(cc.constants, Cons_idx::rho_w)/get_at(cc.constants, Cons_idx::rho_ice);

    init_particle_collection_1(cc.snow, cc.cloud, cc.coeffs_scr);
    init_particle_collection_2(cc.snow, cc.rain, cc.coeffs_srr);
    init_particle_collection_2(cc.ice, cc.rain, cc.coeffs_irr);
    init_particle_collection_1(cc.ice, cc.cloud, cc.coeffs_icr);
    init_particle_collection_1(cc.hail, cc.rain, cc.coeffs_hrr);
    init_particle_collection_1(cc.graupel, cc.rain, cc.coeffs_grr); // Cosmo uses 2, ICON uses 1
    init_particle_collection_1(cc.hail, cc.cloud, cc.coeffs_hcr);
    init_particle_collection_1(cc.graupel, cc.cloud, cc.coeffs_gcr);
    init_particle_collection_1(cc.snow, cc.ice, cc.coeffs_sic);
    init_particle_collection_1(cc.hail, cc.ice, cc.coeffs_hic);
    init_particle_collection_1(cc.graupel, cc.ice, cc.coeffs_gic);
    init_particle_collection_1(cc.hail, cc.snow, cc.coeffs_hsc);
    init_particle_collection_1(cc.graupel, cc.snow, cc.coeffs_gsc);

    // Setup graupel, snow and ice selfcollection

    cc.graupel.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = M_PI/8.0
        * (2.0*coll_delta_11(cc.graupel, cc.graupel, 0)
           + coll_delta_12(cc.graupel, cc.graupel, 0))
        * sqrt((2.0*coll_theta_11(cc.graupel, cc.graupel, 0)
           - coll_theta_12(cc.graupel, cc.graupel, 0)));

    cc.snow.constants[static_cast<int>(Particle_cons_idx::sc_delta_n)] = (2.0*coll_delta_11(cc.snow, cc.snow, 0)
        + coll_delta_12(cc.snow, cc.snow, 0));
    cc.snow.constants[static_cast<int>(Particle_cons_idx::sc_theta_n)] = (2.0*coll_theta_11(cc.snow, cc.snow, 0)
        - coll_theta_12(cc.snow, cc.snow, 0));

    cc.ice.constants[static_cast<int>(Particle_cons_idx::sc_delta_n)] = coll_delta_11(cc.ice, cc.ice, 0)
        + coll_delta_12(cc.ice, cc.ice, 0)
        + coll_delta_22(cc.ice, cc.ice, 0);
    cc.ice.constants[static_cast<int>(Particle_cons_idx::sc_delta_q)] = coll_delta_11(cc.ice, cc.ice, 0)
        + coll_delta_12(cc.ice, cc.ice, 1)
        + coll_delta_22(cc.ice, cc.ice, 1);
    cc.ice.constants[static_cast<int>(Particle_cons_idx::sc_theta_n)] = coll_theta_11(cc.ice, cc.ice, 0)
        - coll_theta_12(cc.ice, cc.ice, 0)
        + coll_theta_22(cc.ice, cc.ice, 0);
    cc.ice.constants[static_cast<int>(Particle_cons_idx::sc_theta_q)] = coll_theta_11(cc.ice, cc.ice, 0)
        - coll_theta_12(cc.ice, cc.ice, 1)
        + coll_theta_22(cc.ice, cc.ice, 1);
#endif
}

/** @} */ // end of group parametrizations


#endif
