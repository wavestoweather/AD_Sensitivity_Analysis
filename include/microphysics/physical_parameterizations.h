#ifndef PHYSICAL_PARAMETERIZATIONS_H
#define PHYSICAL_PARAMETERIZATIONS_H

#include "codi.hpp"
#include <cmath>
#include "types.h"
#include "constants.h"

#include <boost/math/special_functions/gamma.hpp>


/** @defgroup parametrizations Physical Parametrizations
 * Various functions for physical parametrizations, initializing lookup tables
 * and particle collection constants.
 * @{
 */

/**
 * Structure to hold the new equidistant lookup table for
 * graupel wetgrowth diameter
 */
table_t ltabdminwgg;
gamma_table_t table_g1, table_g2, table_r1, table_r2, table_r3;


/** Init lookup table for the incomplete gamma function.
 * From ICON mo_2mom_mcrph_util.f90 incgfct_lower_lookupcreate
 * Create a lookup-table for the lower incomplete gamma function
 * \f[ \text{int}(0)(x) \exp(-t) t^{a-1} \text{d}t \f]
 * with constant a from x=0 to the 99.5% value of the normalized incomplete
 * gamma function. This 99.5 % - value has been fitted
 * with high accuracy as function of a in the range a in [0;20], but can
 * safely be applied also to higher values of a. (Fit created with the
 * matlab-program "gamma_unvoll_lower_lookup.m" by Ulrich Blahak, 2008/11/13).
 *
 * The last value in the table corresponds to x = infinity, so that
 * during the reconstruction of incgfct-values from the table,
 * the x-value can safely be truncated at the maximum table x-value.
 *
 * @param t Gamma table that shall be initialized.
 * @param nl Number of bins in the table.
 * @param nl_highres Optional number of bins for the high resolution table.
 * @param a Value where the incomplete gamma function shall be interpolated around.
 */
void init_gamma_table(
    gamma_table_t &t,
    const uint64_t &nl,
    const uint64_t &nl_highres,
    const double &a)
{

    const double c1 = 36.629433904824623;
    const double c2 = -0.119475603955226;
    const double c3 = 0.339332937820052;
    const double c4 = 1.156369000458310;

    t.n_bins = nl;
    t.n_bins_highres = nl_highres;
    t.a = a;

    t.x.resize(nl);
    t.x_highres.resize(nl_highres);
    t.igf.resize(nl);
    t.igf_highres.resize(nl_highres);

    // Low resolution
    // maximum x-value (99.5%)
    t.x[t.n_bins-2] = c1 * (1.0-exp(c2*pow(a, c3))) + c4*a;
    t.dx = t.x[t.n_bins-2] / (t.n_bins-2.0);
    t.odx = 1.0/t.dx;
    for(uint64_t i=0; i<t.n_bins-2; ++i)
    {
        t.x[i] = (i-1) * t.dx;
        t.igf[i] = boost::math::tgamma_lower(a, t.x[i]);
    }
    // for x -> infinity:
    t.x[t.n_bins-1] = (t.n_bins-1)*t.dx;
    t.igf[t.n_bins-1] = std::tgamma(a);

    // High resolution (lowest 2% of the x-values)
    t.dx_highres = t.x[std::round(0.01*(t.n_bins-1))] / (t.n_bins_highres - 1.0);
    t.odx_highres = 1.0/t.dx_highres;
    for(uint64_t i=0; i<t.n_bins_highres; ++i)
    {
        t.x_highres[i] = (i-1) * t.dx_highres;
        t.igf_highres[i] = boost::math::tgamma_lower(a, t.x_highres[i]);
    }

}


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
inline A thermal_conductivity_moist_air(A T, A qv)
{
  // Thermal conductivity dry air
  A Kt = thermal_conductivity_dry_air(T);

  A Kt_tilde = Kt/418.68;
  A Kv_tilde = (3.78e-5) + (2.0e-7)*(T - 273.15);

  return ( Kt*( 1.0 - (1.17 - 1.02*(Kt_tilde/Kv_tilde))*(qv/(qv+Epsilon)) ) );
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
 * Validity range: ?
 *
 * @param T Temperature in Kelvin
 * @return Specific heat capacity
 */
template <class A>
inline A specific_heat_dry_air(A T)
{

  return 1005.0;

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
inline A specific_heat_ice(A T)
{

  A T_frac = T/125.1;

  return ( (-2.0572 + 0.14644*T + 0.06163*T*exp( -T_frac*T_frac ))/Mw );

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
inline A latent_heat_water(A T)
{

  return ( ( 56579.0 - 42.212*T + exp( 0.1149*(281.6-T) ) )/Mw );

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
inline A latent_heat_ice(A T)
{

  A T_frac = T/123.75;

  return ( ( polyval2(46782.5, 35.8925, -0.07414, T) + 541.5*exp(-T_frac*T_frac) )/Mw );

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
inline A latent_heat_melt(A T)
{
  // Table A1
  return 4.184e3 * (79.7+0.485*(T-tmelt) - 2.5e-3*(T-tmelt)*(T-tmelt));
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
  return lh_e0 * pow(tmelt/T, gam);
}


/**
 * Saturation vapor pressure over a flat surface of liquid water in Pa.
 * From ICON.
 * Validity range: ?
 *
 * @param T Temperature in Kelvin
 * @return Saturation vapor pressure
 */
template <class A>
inline A saturation_pressure_water_icon(A T)
{
  return 610.78*exp(17.269*(T-273.15)/(T-35.86));
}


/**
 * Saturation vapor pressure over a flat surface of liquid water in Pa.
 * From Murphy and Koop (2005).
 * Validity range: \f$ 123 \text{K} <= \text{T} <= 332 \text{K} \f$
 *
 * @param T Temperature in Kelvin
 * @return Saturation vapor pressure
 */
template <class A>
inline A saturation_pressure_water(A T)
{

  A Tinv = 1.0/T;
  A logT = log(T);

  return ( exp( 54.842763 - 6763.22*Tinv - 4.21*logT + 0.000367*T + tanh(0.0415*(T-218.8))*(53.878 - 1331.22*Tinv - 9.44523*logT + 0.014025*T) ) );

}


/**
 * Saturation vapor pressure over a flat surface of ice in Pa.
 * From Murphy and Koop (2005).
 * Validity range: \f$ 110 \text{K} <= \text{T} \f$
 *
 * @param T Temperature in Kelvin
 * @return Saturation vapor pressure
 */
template <class A>
inline A saturation_pressure_ice(A T)
{

  return ( exp( 9.550426 - (5723.265/T) + 3.53068*log(T) - 0.00728332*T ) );

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
		    A S)
{

  return ( S*saturation_pressure_water(T) );

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
		    A S)
{

  return (p - compute_pv(T,S));

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
		      A S)
{

  return ( compute_pa(p,T,S)/(Ra*T) );

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
			 A S)
{

  return ( Epsilon*( compute_pv(T,S)/compute_pa(p,T,S) ) );

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
			  A Si)
{

  A S = Si * ( saturation_pressure_ice(T)/saturation_pressure_water(T) );

  return ( Epsilon*( compute_pv(T,S)/compute_pa(p,T,S) ) );

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
			 A qv)
{

  return ( (p*qv)/((Epsilon + qv)*saturation_pressure_water(T)) );

}


/**
 * Get mean mass of particle assuming a gamma distribution (?).
 *
 * @param q Particle mixing ratio
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

// As seen in ICON
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

// From Seifert & Beheng (2008), Eq. 20
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
        // uint64_t iu =  (uint64_t) min(tmp_d, (codi::RealReverse) (table.n1-1 )).getValue();
        uint64_t io = iu + 1;
        A T_lok = min( max( T, table.x2[0] ), table.x2[table.n2-1] );
        tmp = (T_lok - table.x2[0]) * table.odx2;
        tmp_d = floor(tmp)+1;
        uint64_t ju = std::min(tmp_d, table.n2-1);
        uint64_t jo = ju + 1;
        A qw_lok = min( max( qw, table.x3[0] ), table.x3[table.n3-1] );
        tmp = (qw_lok - table.x3[0]) * table.odx3;
        tmp_d = floor(tmp)+1;
        uint64_t ku = std::min(tmp_d, table.n3-1 );
        uint64_t ko = ku + 1;
        A qi_lok = min( max( qi, table.x4[0] ), table.x4[table.n4-1] );
        tmp = (qi_lok - table.x4[0]) * table.odx4;
        tmp_d = floor(tmp)+1;
        uint64_t lu = std::min(tmp_d, table.n4-1);
        uint64_t lo = lu + 1;

        std::vector<A> h1(16);
        std::vector<A> h2(8);
        std::vector<A> h3(4);
        std::vector<A> h4(2);
        // Tetra linear interpolation by Dmin
        for(uint64_t i=0; i<16; ++i)
            h1[i] = table.get(iu + i/8, ju + (i%8)/4, ku + (i%4)/2, lo+i%2);

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

    }
    return dmin_loc;
}

////////////////////////////////////////////////////////////////////////////////
// Variuous functions for collision integrals of SB2006 (Eq. 90, 91, 92, 93)

// Eq. 90
template <class A = codi::RealReverse>
inline A coll_delta(
    particle_model_constants_t &pc1,
    const uint64_t n)
{
    A tmp1 = (2.0*pc1.b_geo+pc1.nu+1.0+n)/pc1.mu;
    A tmp2 = (pc1.nu+1.0)/pc1.mu;
    A tmp3 = (pc1.nu+1.0)/pc1.mu;
    A tmp4 = (pc1.nu+2.0)/pc1.mu;
    return tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), 2.0*pc1.b_geo+n)
        / pow(tgamma( tmp4 ), 2.0*pc1.b_geo+n);
}

template <class A = codi::RealReverse>
inline A coll_delta_11(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    return coll_delta<A>(pc1, n);
}

template <class A = codi::RealReverse>
inline A coll_delta_22(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    return coll_delta<A>(pc2, n);
}

template <class A = codi::RealReverse>
inline A coll_delta_12(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    A tmp1 = (pc1.b_geo+pc1.nu+1.0)/pc1.mu;
    A tmp2 = (pc1.nu+1.0)/pc1.mu;
    A tmp3 = (pc1.nu+1.0)/pc1.mu;
    A tmp4 = (pc1.nu+2.0)/pc1.mu;
    A res =  2.0 * tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), pc1.b_geo)
        / pow(tgamma( tmp4 ), pc1.b_geo);
    tmp1 = (pc2.b_geo+pc2.nu+1.0+n)/pc2.mu;
    tmp2 = (pc2.nu+1.0)/pc2.mu;
    tmp3 = (pc2.nu+1.0)/pc2.mu;
    tmp4 = (pc2.nu+2.0)/pc2.mu;
    return res * tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), pc2.b_geo+n)
        / pow(tgamma( tmp4 ), pc2.b_geo+n);
}

// Eq. 92
template <class A = codi::RealReverse>
inline A coll_theta(
    particle_model_constants_t &pc1,
    const uint64_t n)
{
    A tmp1 = (2.0*pc1.b_vel+2.0*pc1.b_geo+pc1.nu+1.0)/pc1.mu;
    A tmp2 = (2.0*pc1.b_geo+pc1.nu+1.0+n)/pc1.mu;
    A tmp3 = (pc1.nu+1.0)/pc1.mu;
    A tmp4 = (pc1.nu+2.0)/pc1.mu;
    return tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), 2.0*pc1.b_vel)
        / pow(tgamma( tmp4 ), 2.0*pc1.b_vel);
}

template <class A = codi::RealReverse>
inline A coll_theta_11(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    return coll_theta<A>(pc1, n);
}

template <class A = codi::RealReverse>
inline A coll_theta_22(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    return coll_theta<A>(pc2, n);
}

// Eq. 93
template <class A = codi::RealReverse>
inline A coll_theta_12(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const uint64_t n)
{
    A tmp1 = (pc1.b_vel+2.0*pc1.b_geo+pc1.nu+1.0)/pc1.mu;
    A tmp2 = (2.0*pc1.b_geo+pc1.nu+1.0)/pc1.mu;
    A tmp3 = (pc1.nu+1.0)/pc1.mu;
    A tmp4 = (pc1.nu+2.0)/pc1.mu;
    A res = 2.0 * tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), pc1.b_vel)
        / pow(tgamma( tmp4 ), pc1.b_vel);
    tmp1 = (pc2.b_vel+2.0*pc2.b_geo+pc2.nu+1.0+n)/pc2.mu;
    tmp2 = (2.0*pc2.b_geo+pc2.nu+1.0+n)/pc2.mu;
    tmp3 = (pc2.nu+1.0)/pc2.mu;
    tmp4 = (pc2.nu+2.0)/pc2.mu;
    return res * tgamma( tmp1 )
        / tgamma( tmp2 )
        * pow(tgamma( tmp3 ), pc2.b_vel)
        / pow(tgamma( tmp4 ), pc2.b_vel);
}


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

inline codi::RealReverse vent_coeff_a(
    particle_model_constants_t &pc,
    uint64_t n)
{
    return pc.a_ven * tgamma( (pc.nu+n+pc.b_geo)/pc.mu )
        / tgamma( (pc.nu+1.0)/pc.mu )
        * pow( tgamma((pc.nu+1.0)/pc.mu) / tgamma((pc.nu+2.0)/pc.mu), pc.b_geo+n-1.0 );
}

inline codi::RealReverse vent_coeff_b(
    particle_model_constants_t &pc,
    uint64_t n)
{
    const double m_f = 0.5; // From PK, page 541
    return pc.b_ven * tgamma( (pc.nu+n+(m_f+1.0)*pc.b_geo+m_f*pc.b_vel)/pc.mu )
        / tgamma( (pc.nu+1.0)/pc.mu )
        * pow( tgamma( (pc.nu+1.0)/pc.mu )
            / tgamma( (pc.nu+2.0)/pc.mu), (m_f+1.0) * pc.b_geo + m_f*pc.b_vel+n-1.0);
}

// complete mass moment of particle size distribution, Eq (82) of SB2006
inline codi::RealReverse moment_gamma(
    particle_model_constants_t &pc,
    uint64_t n)
{
    return tgamma( (n+pc.nu+1.0)/pc.mu ) / tgamma( (pc.nu+1.0)/pc.mu )
        * pow(tgamma( (pc.nu+1.0)/pc.mu ) / tgamma( (pc.nu+2.0)/pc.mu), n);
}

/** @} */ // end of group parametrizations


#endif
