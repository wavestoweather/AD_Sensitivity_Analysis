#pragma once

#include <math.h>

#include <algorithm>
#include <vector>

#include <codi.hpp>

#include "include/microphysics/constants.h"
#include "include/microphysics/user_functions.h"
#include "include/types/model_constants_t.h"
#include "include/types/reference_quantities_t.h"


/** @defgroup rk Runge-Kutta 4 Method
 * This file provides the function for a single step with the
 * Runge-Kutta 4 method.
 * @{
 */


template<class float_t>
void set_limits(
    std::vector<float_t> &y,
    const reference_quantities_t& ref,
    model_constants_t<float_t> &cc) {

    // Set everything negative and very small to zero.
    // Those should be artifacts from instantaneous processes
    y[qv_idx] = (y[qv_idx]*ref.qref < artifact_thresh) ? 0 : y[qv_idx];
    y[qc_idx] = (y[qc_idx]*ref.qref < artifact_thresh) ? 0 : y[qc_idx];
    y[qr_idx] = (y[qr_idx]*ref.qref < artifact_thresh) ? 0 : y[qr_idx];
    y[qs_idx] = (y[qs_idx]*ref.qref < artifact_thresh) ? 0 : y[qs_idx];
    y[qi_idx] = (y[qi_idx]*ref.qref < artifact_thresh) ? 0 : y[qi_idx];
    y[qg_idx] = (y[qg_idx]*ref.qref < artifact_thresh) ? 0 : y[qg_idx];
    y[qh_idx] = (y[qh_idx]*ref.qref < artifact_thresh) ? 0 : y[qh_idx];
    y[Nc_idx] = (y[Nc_idx] < artifact_thresh) ? 0 : y[Nc_idx];
    y[Nr_idx] = (y[Nr_idx] < artifact_thresh) ? 0 : y[Nr_idx];
    y[Ns_idx] = (y[Ns_idx] < artifact_thresh) ? 0 : y[Ns_idx];
    y[Ni_idx] = (y[Ni_idx] < artifact_thresh) ? 0 : y[Ni_idx];
    y[Ng_idx] = (y[Ng_idx] < artifact_thresh) ? 0 : y[Ng_idx];
    y[Nh_idx] = (y[Nh_idx] < artifact_thresh) ? 0 : y[Nh_idx];

    if (nuc_type > 0)
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
 * @param fixed If True: Do not change pressure, temperature and ascent (w)
 */
template<class float_t>
void RK4_step_2_sb_ice(
    std::vector<float_t> &ynew,
    std::vector<float_t> &yold,
    const reference_quantities_t& ref,
    model_constants_t<float_t>& cc,
    bool fixed) {

    // Define temporary variables
    std::vector<float_t> k(num_comp);
    std::vector<float_t> ytmp(num_comp);

    // Reset the result vector with the current state
    for (int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];
    //
    // Do all computations involving k1
    //
    RHS_SB(k, yold, ref, cc, cc.dt_half, fixed);  // k = k1

    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii];  // Add k1-part to the result
    }
    set_limits(ytmp, ref, cc);
    sediment_q_total += cc.dt_sixth*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt_sixth*sediment_n;
    sediment_n = 0;

    //
    // Do all computations involving k2
    //
    RHS_SB(k, ytmp, ref, cc, cc.dt_half, fixed);  // k = k2

    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii];  // Add k2-part to the result
    }
    set_limits(ytmp, ref, cc);
    sediment_q_total += cc.dt_third*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt_third*sediment_n;
    sediment_n = 0;

    //
    // Do all computations involving k3
    //
    RHS_SB(k, ytmp, ref, cc, cc.dt_half, fixed);  // k = k3

    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt*k[ii];  // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii];  // Add k3-part to the result
    }
    set_limits(ytmp, ref, cc);
    sediment_q_total += cc.dt_third*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt_third*sediment_n;
    sediment_n = 0;

    //
    // Do all computations involving k4
    //
    RHS_SB(k, ytmp, ref, cc, cc.dt_half, fixed);  // k = k4

    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] += cc.dt_sixth*k[ii];
    }
    set_limits(ynew, ref, cc);
    sediment_q_total += cc.dt_sixth*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt_sixth*sediment_n;
    sediment_n = 0;
    // Explicit calculation of saturation
    float_t T_prime = ynew[T_idx]*ref.Tref;
    float_t p_prime = ynew[p_idx]*ref.pref;
    float_t qv_prime = ynew[qv_idx]*ref.qref;
    ynew[S_idx] = convert_qv_to_S(p_prime, T_prime, qv_prime,
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                get_at(cc.constants, Cons_idx::Epsilon));
#ifdef TRACE_SAT
    std::cout << "traj: " << cc.traj_id << ", S (end of RK4): " << ynew[S_idx] << "\n";
#endif
    cc.increment_w_idx();
}

/** @} */  // end of group rk
