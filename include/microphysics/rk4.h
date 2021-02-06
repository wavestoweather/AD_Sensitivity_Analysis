#ifndef RK4_H
#define RK4_H

#include "codi.hpp"
#include <math.h>
#include <vector>

#include "include/microphysics/constants.h"
#include "include/types/model_constants_t.h"
#include "include/types/nc_parameters_t.h"
#include "include/types/reference_quantities_t.h"

#include "user_functions.h"

/** @defgroup rk Runge-Kutta 4 Method
 * This file provides the function for a single step with the
 * Runge-Kutta 4 method.
 * @{
 */

/**
 * Compute a single step using the Runge-Kutta 4 method for the ODE
 * \f[ y' = \text{RHS}(y) \f]
 * with old system state yold and new system state ynew. It uses the
 * 1 moment bulk scheme from 10.1175/JAS3980
 *
 * @param ynew On out: new system state
 * @param yold Old system state
 * @param ref Pointer to reference quantities to transform between units
 * @param cc Pointer to constants from the model
 * @param nc Pointer to parameters from the netCDF file
 * @param fixed If True: Do not change pressure, temperature and ascent (w)
 */
void RK4_step(
    std::vector<codi::RealReverse> &ynew,
    std::vector<codi::RealReverse> &yold,
    const reference_quantities_t& ref,
    model_constants_t& cc,
    nc_parameters_t& nc,
    bool fixed)
{
    uint32_t update_idx = 3;
    // Define temporary variables
    std::vector<codi::RealReverse> k(num_comp);
    std::vector<codi::RealReverse> ytmp(num_comp);

    // Reset the result vector with the current state
    for(int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];

    //
    // Do all computations involving k1
    //
#if defined(RK4_ONE_MOMENT)
    RHS(k, yold, ref, cc, nc); // k = k1
#endif
#if defined(OTHER)
    RHS_other();
#endif

    for(int ii = 0 ; ii < num_comp ; ii++){
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii]; // Add k1-part to the result
#ifdef TRACE_ENV
        if(trace && ii == T_idx)
            std::cout << "RK1, yold[T] = " << yold[ii]
                      << "\nk[T] = " << k[ii]
                      << "\ndT = " << yold[ii] + cc.dt_half*k[ii]
                      << "\ndT Result = " << cc.dt_sixth*k[ii]
                      << "\n";
#endif
    }
    if(fixed)
        for(uint32_t ii=0; ii<update_idx; ++ii)
            ytmp[ii] = yold[ii];

    //
    // Do all computations involving k2
    //
#if defined(RK4_ONE_MOMENT)
    RHS(k, ytmp, ref, cc, nc); // k = k2
#endif

    for(int ii = 0 ; ii < num_comp ; ii++){
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii]; // Add k2-part to the result
#ifdef TRACE_ENV
        if(trace && ii == T_idx)
            std::cout << "RK2, yold[T] = " << yold[ii]
                      << "\nk[T] = " << k[ii]
                      << "\ndT = " << yold[ii] + cc.dt_half*k[ii]
                      << "\ndT Result = " << cc.dt_third*k[ii]
                      << "\n";
#endif
    }
    if(fixed)
        for(uint32_t ii=0; ii<update_idx; ++ii)
            ytmp[ii] = yold[ii];

    //
    // Do all computations involving k3
    //
#if defined(RK4_ONE_MOMENT)
    RHS(k, ytmp, ref, cc, nc); // k = k3
#endif

    for(int ii = 0 ; ii < num_comp ; ii++){
        ytmp[ii] = yold[ii] + cc.dt*k[ii]; // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii]; // Add k3-part to the result
#ifdef TRACE_ENV
        if(trace && ii == T_idx)
            std::cout << "RK3, yold[T] = " << yold[ii]
                      << "\nk[T] = " << k[ii]
                      << "\ndT = " << yold[ii] + cc.dt*k[ii]
                      << "\ndT Result = " << cc.dt_third*k[ii]
                      << "\n";
#endif
    }
    if(fixed)
        for(uint32_t ii=0; ii<update_idx; ++ii)
            ytmp[ii] = yold[ii];

    //
    // Do all computations involving k4
    //
#if defined(RK4_ONE_MOMENT)
    RHS(k, ytmp, ref, cc, nc); // k = k4
#endif

    for(int ii = 0 ; ii < num_comp ; ii++){
        ynew[ii] += cc.dt_sixth*k[ii];
#ifdef TRACE_ENV
        if(trace && ii == T_idx)
            std::cout << "RK4, k[T] = " << k[ii]
                      << "\ndT Result = " << cc.dt_sixth*k[ii]
                      << "\n";
#endif
    }

    if(fixed)
        for(uint32_t ii=0; ii < update_idx; ii++)
            ynew[ii] = yold[ii];

}


/**
 * Compute a single step using the Runge-Kutta 4 method for the ODE
 * \f[ y' = \text{RHS}(y) \f]
 * with old system state yold and new system state ynew.  It uses the
 * 2 moment bulk scheme from 10.1175/JAS3980 without ice.
 *
 * @param ynew On out: new system state
 * @param yold Old system state
 * @param ref Pointer to reference quantities to transform between units
 * @param cc Pointer to constants from the model
 * @param nc Pointer to parameters from the netCDF file
 * @param fixed If True: Do not change pressure, temperature and ascent (w)
 */
void RK4_step_2_no_ice(
    std::vector<codi::RealReverse> &ynew,
    std::vector<codi::RealReverse> &yold,
    const reference_quantities_t& ref,
    model_constants_t& cc,
    nc_parameters_t& nc,
    bool fixed)
{
    // Define temporary variables
    std::vector<codi::RealReverse> k(num_comp);
    std::vector<codi::RealReverse> ytmp(num_comp);

    // Reset the result vector with the current state
    for(int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];

    //
    // Do all computations involving k1
    //
    RHS_SB_no_ice(k, yold, ref, cc, nc, fixed); // k = k1

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii]; // Add k1-part to the result

    }

    //
    // Do all computations involving k2
    //
    RHS_SB_no_ice(k, ytmp, ref, cc, nc, fixed); // k = k2

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii]; // Add k2-part to the result
    }

    //
    // Do all computations involving k3
    //
    RHS_SB_no_ice(k, ytmp, ref, cc, nc, fixed); // k = k3

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        ytmp[ii] = yold[ii] + cc.dt*k[ii]; // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii]; // Add k3-part to the result
    }

    //
    // Do all computations involving k4
    //
    RHS_SB_no_ice(k, ytmp, ref, cc, nc, fixed); // k = k4

    for(int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] += cc.dt_sixth*k[ii];

}


template<class float_t>
void set_limits(
    std::vector<float_t> &y,
    const reference_quantities_t& ref,
    model_constants_t &cc)
{
    if(nuc_type > 0)
        y[Nc_idx] = min(min(max(y[Nc_idx], y[qc_idx]*ref.qref/get_at(cc.cloud.constants, Particle_cons_idx::max_x)),
            y[qc_idx]*ref.qref/get_at(cc.cloud.constants, Particle_cons_idx::min_x)), 5e9);
    y[Nr_idx] = min(max(y[Nr_idx], y[qr_idx]*ref.qref/get_at(cc.rain.constants, Particle_cons_idx::max_x)),
        y[qr_idx]*ref.qref/get_at(cc.rain.constants, Particle_cons_idx::min_x));
    y[Ni_idx] = min(max(y[Ni_idx], y[qi_idx]*ref.qref/get_at(cc.ice.constants, Particle_cons_idx::max_x)),
        y[qi_idx]*ref.qref/get_at(cc.ice.constants, Particle_cons_idx::min_x));
    y[Ns_idx] = min(max(y[Ns_idx], y[qs_idx]*ref.qref/get_at(cc.snow.constants, Particle_cons_idx::max_x)),
        y[qs_idx]*ref.qref/get_at(cc.snow.constants, Particle_cons_idx::min_x));
    y[Ng_idx] = min(max(y[Ng_idx], y[qg_idx]*ref.qref/get_at(cc.graupel.constants, Particle_cons_idx::max_x)),
        y[qg_idx]*ref.qref/get_at(cc.graupel.constants, Particle_cons_idx::min_x));
    y[Nh_idx] = min(max(y[Nh_idx], y[qh_idx]*ref.qref/get_at(cc.hail.constants, Particle_cons_idx::max_x)),
        y[qh_idx]*ref.qref/get_at(cc.hail.constants, Particle_cons_idx::min_x));
    // Set everything negative to zero
    y[qv_idx] = (y[qv_idx] < 0) ? 0 : y[qv_idx];
    y[qc_idx] = (y[qc_idx] < 0) ? 0 : y[qc_idx];
    y[qr_idx] = (y[qr_idx] < 0) ? 0 : y[qr_idx];
    y[qs_idx] = (y[qs_idx] < 0) ? 0 : y[qs_idx];
    y[qi_idx] = (y[qi_idx] < 0) ? 0 : y[qi_idx];
    y[qg_idx] = (y[qg_idx] < 0) ? 0 : y[qg_idx];
    y[qh_idx] = (y[qh_idx] < 0) ? 0 : y[qh_idx];
    y[Nc_idx] = (y[Nc_idx] < 0) ? 0 : y[Nc_idx];
    y[Nr_idx] = (y[Nr_idx] < 0) ? 0 : y[Nr_idx];
    y[Ns_idx] = (y[Ns_idx] < 0) ? 0 : y[Ns_idx];
    y[Ni_idx] = (y[Ni_idx] < 0) ? 0 : y[Ni_idx];
    y[Ng_idx] = (y[Ng_idx] < 0) ? 0 : y[Ng_idx];
    y[Nh_idx] = (y[Nh_idx] < 0) ? 0 : y[Nh_idx];
}

/**
 * Compute a single step using the Runge-Kutta 4 method for the ODE
 * \f[ y' = \text{RHS}(y) \f]
 * with old system state yold and new system state ynew.  It uses the
 * 2 moment bulk scheme from 10.1175/JAS3980 with ice.
 *
 * @param ynew On out: new system state
 * @param yold Old system state
 * @param ref Pointer to reference quantities to transform between units
 * @param cc Pointer to constants from the model
 * @param nc Pointer to parameters from the netCDF file
 * @param fixed If True: Do not change pressure, temperature and ascent (w)
 */
void RK4_step_2_sb_ice(
    std::vector<codi::RealReverse> &ynew,
    std::vector<codi::RealReverse> &yold,
    const reference_quantities_t& ref,
    model_constants_t& cc,
    nc_parameters_t& nc,
    bool fixed)
{

    // Define temporary variables
    std::vector<codi::RealReverse> k(num_comp);
    std::vector<codi::RealReverse> ytmp(num_comp);

    // Reset the result vector with the current state
    for(int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];

    //
    // Do all computations involving k1
    //
    RHS_SB(k, yold, ref, cc, nc, cc.dt_sixth, fixed); // k = k1

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
#ifdef TRACE_ENV
        if(trace && ii == T_idx)
            std::cout << "RK1, yold[T] = " << yold[ii]
                      << "\nk[T] = " << k[ii]
                      << "\nk[T] prime = " << k[ii]*ref.Tref
                      << "\nk[T] prime time = " << k[ii]*ref.Tref*cc.dt_half
                      << "\ndT = " << yold[ii] + cc.dt_half*k[ii]
                      << "\ndT Result = " << cc.dt_sixth*k[ii]
                      << "\ndT prime Result = " << cc.dt_sixth*k[ii]*ref.Tref
                      << "\n";
#endif
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii]; // Add k1-part to the result
    }
    set_limits(ytmp, ref, cc);
    sediment_q_total += cc.dt_sixth*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt_sixth*sediment_n;
    sediment_n = 0;

    //
    // Do all computations involving k2
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_third, fixed); // k = k2

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
#ifdef TRACE_ENV
        if(trace && ii == T_idx)
            std::cout << "RK2, yold[T] = " << yold[ii]
                      << "\nk[T] = " << k[ii]
                      << "\nk[T] prime = " << k[ii]*ref.Tref
                      << "\nk[T] prime time = " << k[ii]*ref.Tref*cc.dt_half
                      << "\ndT = " << yold[ii] + cc.dt_half*k[ii]
                      << "\ndT Result = " << cc.dt_third*k[ii]
                      << "\ndT prime Result = " << cc.dt_third*k[ii]*ref.Tref
                      << "\n";
#endif
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii]; // Add k2-part to the result
    }
    set_limits(ytmp, ref, cc);
    sediment_q_total += cc.dt_third*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt_third*sediment_n;
    sediment_n = 0;

    //
    // Do all computations involving k3
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_third, fixed); // k = k3

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
#ifdef TRACE_ENV
        if(trace && ii == T_idx)
            std::cout << "RK3, yold[T] = " << yold[ii]
                      << "\nk[T] = " << k[ii]
                      << "\nk[T] prime = " << k[ii]*ref.Tref
                      << "\nk[T] prime time = " << k[ii]*ref.Tref*cc.dt
                      << "\ndT = " << yold[ii] + cc.dt*k[ii]
                      << "\ndT Result = " << cc.dt_third*k[ii]
                      << "\ndT prime Result = " << cc.dt_third*k[ii]*ref.Tref
                      << "\n";
#endif
        ytmp[ii] = yold[ii] + cc.dt*k[ii]; // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii]; // Add k3-part to the result
    }
    set_limits(ytmp, ref, cc);
    sediment_q_total += cc.dt_third*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt_third*sediment_n;
    sediment_n = 0;

    //
    // Do all computations involving k4
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_sixth, fixed); // k = k4

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
#ifdef TRACE_ENV
        if(trace && ii == T_idx)
            std::cout << "RK4, k[T] = " << k[ii]
                        << "\nk[T] prime = " << k[ii]*ref.Tref
                      << "\ndT Result = " << cc.dt_sixth*k[ii]
                      << "\ndT prime Result = " << cc.dt_sixth*k[ii]*ref.Tref
                      << "\n";
#endif
        ynew[ii] += cc.dt_sixth*k[ii];
    }
    set_limits(ynew, ref, cc);
    sediment_q_total += cc.dt_sixth*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt_sixth*sediment_n;
    sediment_n = 0;
    // Explicit calculation of saturation
    codi::RealReverse T_prime = ynew[T_idx]*ref.Tref;
    codi::RealReverse p_prime = ynew[p_idx]*ref.pref;
    codi::RealReverse qv_prime = ynew[qv_idx]*ref.qref;
    ynew[S_idx] = convert_qv_to_S(p_prime, T_prime, qv_prime,
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                 get_at(cc.constants, Cons_idx::Epsilon));
}

/** @} */ // end of group rk

#endif
