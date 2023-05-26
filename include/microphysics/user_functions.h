#pragma once

#include <math.h>
#include <tgmath.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <codi.hpp>

#include "include/microphysics/constants.h"
#include "include/microphysics/physical_parameterizations.h"
#include "include/misc/error.h"
#include "include/types/model_constants_t.h"
#include "include/types/reference_quantities_t.h"


/** @defgroup ufunc User Functions
 * Functions to be defined by the user, specifically the definition of
 * the right-hand side function RHS of the ODE
 * \f[ y' = \text{RHS}(y) \f]
 * which is solved.
 * @{
 */

codi::RealReverse trunc(codi::RealReverse x) {
    int left = x.getValue();
    codi::RealReverse r = left;
    r.setGradient(x.getGradient());
    return r;
}

codi::RealForwardVec<num_par_init> trunc(codi::RealForwardVec<num_par_init> x) {
    int left = x.getValue();
    codi::RealForwardVec<num_par_init> r = left;
    for (auto i=0; i < num_par_init; i++) {
        r.gradient()[i] = x.getGradient()[i];
    }
    return r;
}

template<class float_t>
void precipitation_efficiency(
    std::vector<float_t> &res
) {
    float_t formed_hydro = res[qr_idx];
    float_t prec_hydro = res[qr_out_idx];
#ifndef RK4_ONE_MOMENT
    formed_hydro += res[qi_idx] + res[qs_idx] + res[qg_idx] + res[qh_idx];
    prec_hydro +=  res[qi_out_idx] + res[qs_out_idx] + res[qg_out_idx] + res[qh_out_idx];
#endif
    res[out_eff_idx] += formed_hydro/prec_hydro;
}

/**
 * Temperature and water vapor contents are adjusted. Heat capacity is calculated using
 * dry air (errors should be rather small).
 * Total density is calculated using moist air
 *
 * @params T_prime Current temperature [K]
 * @params p_prime Current pressure [Pa]
 * @params p_sat Pressure at saturation [Pa]
 * @params qv_prime Current water vapor mixing ratio
 * @params qc_prime Current cloud water mixing ratio
 * @params res Vector to store the changes at qv_idx, qc_idx and T_idx
 */
template<class float_t>
void saturation_adjust(
    const float_t &T_prime,
    const float_t &p_prime,
    const float_t &qv_prime,
    const float_t &qc_prime,
    std::vector<float_t> &res,
    const model_constants_t<float_t> &cc) {
    float_t S = convert_qv_to_S(
        p_prime,
        T_prime,
        qv_prime,
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b),
        get_at(cc.constants, Cons_idx::Epsilon));

#ifdef TRACE_ENV
        if (trace)
        std::cout << "\nbefore sat ad S calc "
            << S
            << ", QC: " << qc_prime << "\n";
#endif

    if (S < 1 && qc_prime <= 0) return;  // Nothing to evaporate.
    float_t q_total = qv_prime + qc_prime;

    float_t lat_heat_vapor = latent_heat_water(
        T_prime,
        get_at(cc.constants, Cons_idx::L_wd),
        get_at(cc.constants, Cons_idx::cv),
        get_at(cc.constants, Cons_idx::cp),
        get_at(cc.constants, Cons_idx::T_freeze),
        get_at(cc.constants, Cons_idx::R_v),
        get_at(cc.constants, Cons_idx::R_a));
    float_t T_test = T_prime - lat_heat_vapor * qc_prime;

    float_t rho_total = compute_rhoh(
        p_prime, T_prime, S,
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b),
        get_at(cc.constants, Cons_idx::R_a),
        get_at(cc.constants, Cons_idx::R_v));
    float_t p_vapor_sat = compute_pv(
        T_test,
        float_t(1),
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b));
    float_t q_test = p_vapor_sat
        / (rho_total * get_at(cc.constants, Cons_idx::R_v) * T_test);

    // not saturated, even if all qc evaporates?
    if (q_total <= q_test) {
#ifdef IN_SAT_ADJ
        res[qv_idx] += res[qc_idx] + qc_prime/dt;
        res[qc_idx] -= res[qc_idx] + qc_prime/dt;
        res[T_idx] -= lat_heat_vapor * (res[qc_idx] + qc_prime/dt);
#if defined(TRACE_ENV) || defined(TRACE_QV) || defined(TRACE_QC)
        if (trace)
            std::cout << "traj: " << cc.traj_id << " saturation_adjust (in, new) dT "
                << -lat_heat_vapor * (res[qc_idx] + qc_prime/dt) << ", "
                << "dqv: " << res[qc_idx] + qc_prime/dt << ", "
                << "dqc: " << -res[qc_idx] + qc_prime/dt << ", "
                << "T_prime: " << T_prime << ", "
                << "qv_prime: " << qv_prime << ", "
                << "qc_prime: " << qc_prime
                << "\n";
#endif
#else
        res[qv_idx] += q_total - qv_prime;
        res[qc_idx] += -qc_prime;
        res[T_idx] += T_test - T_prime;

        float_t delta_q = -qc_prime;
        float_t delta_e = latent_heat_evap(T_prime) * delta_q / specific_heat_dry_air(T_prime);
        if (delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        else
            res[lat_heat_idx] += delta_e;
#if defined(TRACE_ENV) || defined(TRACE_QV) || defined(TRACE_QC)
        if (trace)
            std::cout << "traj: " << cc.traj_id << " saturation_adjust (out, new) dT "
                << T_test - T_prime << ", "
                << "dqv: " << q_total - qv_prime << ", "
                << "dqc: " << -qc_prime << ", "
                << "T_prime: " << T_prime << ", "
                << "qv_prime: " << qv_prime << ", "
                << "qc_prime: " << qc_prime
                << "\n";
#endif
#endif

    } else {
        const int maxiter = 10;
        const float_t tolerance = 1e-3;
        // Newton: let (a) be the difference in temperature before and after
        // evaporating cloud droplets and (b) be the difference in temperature
        // from latent heat by evaporating cloud mass. Then f = (a) - (b)
        // shall be zero with no change in the next Newton step.
        float_t T_new = T_prime;
        float_t T_old;
        float_t q_new, dq_new_dT, f, df_dT;
        float_t p_sat_const = get_at(cc.constants, Cons_idx::p_sat_const_a)
            * (get_at(cc.constants, Cons_idx::T_freeze) - get_at(cc.constants, Cons_idx::p_sat_const_b));
        int i = 0;
        do {
            i++;
            T_old = T_new;
            p_vapor_sat = compute_pv(
                T_new,
                float_t(1),
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b));
            rho_total = compute_rhoh(
                p_prime, T_new, S,
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                get_at(cc.constants, Cons_idx::R_a),
                get_at(cc.constants, Cons_idx::R_v));
            q_new = p_vapor_sat / (rho_total * get_at(cc.constants, Cons_idx::R_v) * T_new);
            dq_new_dT =  p_sat_const
                / std::pow(T_new - get_at(cc.constants, Cons_idx::p_sat_const_b), 2)
                - 1/T_new * q_new;
            lat_heat_vapor = latent_heat_water(
                T_new,
                get_at(cc.constants, Cons_idx::L_wd),
                get_at(cc.constants, Cons_idx::cv),
                get_at(cc.constants, Cons_idx::cp),
                get_at(cc.constants, Cons_idx::T_freeze),
                get_at(cc.constants, Cons_idx::R_v),
                get_at(cc.constants, Cons_idx::R_a));
            f = T_new - T_prime + lat_heat_vapor * (q_new - qv_prime);
            df_dT = 1 + lat_heat_vapor * dq_new_dT;
            T_new = T_new - f/df_dT;
        } while ((std::abs(T_new - T_old) > tolerance) && (i < maxiter));
#ifdef IN_SAT_ADJ
        float_t delta_q = p_vapor_sat / (rho_total * get_at(cc.constants, Cons_idx::R_v) * T_new);
        float_t delta_qc = std::max(qc_prime + qv_prime - delta_q, 1.0e-20) - qc_prime;
        if (delta_qc > 0) {
            // cloud droplets evaporate
            delta_qc = std::min(delta_qc, res[qc_idx]+qc_prime/dt);
        } else {
            // water vapor gets cloud droplets
            delta_qc = -std::min(-delta_qc, res[qv_idx]+qv_prime/dt);
        }
        float_t delta_T = lat_heat_vapor * delta_qc - T_prime;
        res[qc_idx] -= delta_qc;
        res[qv_idx] += delta_qc;
        res[T_idx] += delta_T;
#if defined(TRACE_ENV) || defined(TRACE_QV) || defined(TRACE_QC)
        if (trace)
            std::cout << "traj: " << cc.traj_id << " saturation_adjust (in, new, after Newton) dT "
                << delta_T << ", "
                << "dqv: " << delta_qc << ", "
                << "dqc: " << -delta_qc << ", "
                << "T_prime: " << T_prime << ", "
                << "qv_prime: " << qv_prime << ", "
                << "qc_prime: " << qc_prime
                << "\n";
#endif
#else
        res[T_idx] = T_new - T_prime;
        p_vapor_sat = compute_pv(
            T_new,
            float_t(1),
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));
        rho_total = compute_rhoh(
                p_prime, T_new, S,
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                get_at(cc.constants, Cons_idx::R_a),
                get_at(cc.constants, Cons_idx::R_v));
        float_t delta_q = p_vapor_sat / (rho_total * get_at(cc.constants, Cons_idx::R_v) * T_new);
        // I'm paranoid and this case should be in the other branch
        // of the if-clause.
        delta_q = (delta_q > q_total) ? q_total : delta_q;
        float_t delta_qc = std::max(q_total - delta_q, 1.0e-20) - qc_prime;
        res[qc_idx] += delta_qc;
        res[qv_idx] += delta_q - qv_prime;

        float_t delta_e = latent_heat_evap(T_prime) * delta_qc / specific_heat_dry_air(T_prime);
        if (delta_qc < 0.0)
            res[lat_cool_idx] += delta_e;
        else
            res[lat_heat_idx] += delta_e;
#if defined(TRACE_ENV) || defined(TRACE_QV) || defined(TRACE_QC)
        if (trace)
            std::cout << "traj: " << cc.traj_id << " saturation_adjust (out, new, after Newton) dT "
                << res[T_idx] << ", "
                << "dqv: " << res[qv_idx] << ", "
                << "dqc: " << res[qc_idx] << ", "
                << "T_prime: " << T_prime << ", "
                << "qv_prime: " << qv_prime << ", "
                << "qc_prime: " << qc_prime
                << "\n";
#endif
#endif
    }
}


/**
 * Saturation adjustment that also adjusts temperature (no explicit calculation
 * of latent heating/cooling; temperature change is given directly!).
 *
 * @params T_prime Current temperature [K]
 * @params p_prime Current pressure [Pa]
 * @params p_sat Pressure at saturation [Pa]
 * @params qv_prime Current water vapor mixing ratio
 * @params qc_prime Current cloud water mixing ratio
 * @params res Vector to store the changes at qv_idx, qc_idx and T_idx
 */
template<class float_t>
void saturation_adjust_legacy(
    const float_t &T_prime,
    const float_t &p_prime,
    const float_t &qv_prime,
    const float_t &qc_prime,
    std::vector<float_t> &res,
    const model_constants_t<float_t> &cc) {
    float_t qv = qv_prime;
    float_t qc = qc_prime;
    float_t T = T_prime;

    T = T - get_at(cc.constants, Cons_idx::L_wd)*qc_prime/get_at(cc.constants, Cons_idx::cp);
    qv = qc_prime + qv_prime;
    qc = 0;
    float_t delta_q = convert_S_to_qv(p_prime, T, float_t(1),
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b),
        get_at(cc.constants, Cons_idx::Epsilon)) - qv;

    if (delta_q < 0) {
        // Handle over saturation
        // Adjust temperature
        float_t Th = get_at(cc.constants, Cons_idx::cp)*T + get_at(cc.constants, Cons_idx::L_wd)*qv;
        float_t T_qd0 = convert_S_to_qv(
            p_prime, T_prime, float_t(1),
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon));
        float_t T_dt0;
        // Do the Newton
        for (uint32_t i=0; i < 1; ++i) {
            T_dt0 = get_at(cc.constants, Cons_idx::p_sat_const_a)
                * (get_at(cc.constants, Cons_idx::T_sat_low_temp)-get_at(cc.constants, Cons_idx::p_sat_const_b))
                * (1+(get_at(cc.constants, Cons_idx::Epsilon)-1)*T_qd0) * T_qd0
                / ((T_prime-get_at(cc.constants, Cons_idx::p_sat_const_b))
                    * (T_prime-get_at(cc.constants, Cons_idx::p_sat_const_b)));
            T = (Th-get_at(cc.constants, Cons_idx::L_wd)*(T_qd0-T_dt0*T_prime))
                / (get_at(cc.constants, Cons_idx::cp)+get_at(cc.constants, Cons_idx::L_wd)*T_dt0);
            T_qd0 = convert_S_to_qv(p_prime, T, float_t(1),
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                get_at(cc.constants, Cons_idx::Epsilon));
        }

        // Get ratio mixing "back"
        T_dt0 = get_at(cc.constants, Cons_idx::p_sat_const_a)
                * (get_at(cc.constants, Cons_idx::T_sat_low_temp)-get_at(cc.constants, Cons_idx::p_sat_const_b))
                * (1+(get_at(cc.constants, Cons_idx::Epsilon)-1)*T_qd0) * T_qd0
                / ((T-get_at(cc.constants, Cons_idx::p_sat_const_b))
                    * (T-get_at(cc.constants, Cons_idx::p_sat_const_b)));
        float_t T_gn = (Th - get_at(cc.constants, Cons_idx::L_wd)*(T_qd0-T_dt0*T))
            / (get_at(cc.constants, Cons_idx::cp)+get_at(cc.constants, Cons_idx::L_wd)*T_dt0);
        T_qd0 += T_dt0*(T_gn-T);

        qv = T_qd0;
        qc = std::max(qc_prime+qv_prime - T_qd0, 1.0e-20);
        delta_q = qv - qv_prime + res[qv_idx];
#ifdef IN_SAT_ADJ
        if (delta_q > 0)
            delta_q = std::min(delta_q, res[qc_idx] + qc_prime/dt);
        else
            delta_q = -min(abs(delta_q), res[qv_idx] + qv_prime/dt);
        T_gn = T_prime - get_at(cc.constants, Cons_idx::L_wd)*delta_q/get_at(cc.constants, Cons_idx::cp);
#endif
        res[qv_idx] += delta_q;
        res[qc_idx] -= delta_q;
        // T <- T_gn
        res[T_idx] +=  T_gn - T_prime;
#ifdef TRACE_QV
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " saturation_adjust (over saturated for new qv = qv+qc) dqv " << delta_q << "\n";

#endif
#ifdef TRACE_ENV
        if (trace)
            std::cout << "traj: " << cc.traj_id << " saturation_adjust dT " << T_gn - T_prime
                  << "\n";
#endif
#ifdef TRACE_QC
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " saturation_adjust (over saturated for new qv = qv+qc) dqc " << -delta_q << "\n";
#endif
    } else {
#ifdef IN_SAT_ADJ
        res[qv_idx] += res[qc_idx] + qc_prime/dt;
        res[qc_idx] -= res[qc_idx] + qc_prime/dt;
        res[T_idx] += - get_at(cc.constants, Cons_idx::L_wd)
            * (res[qc_idx]+qc_prime)/(get_at(cc.constants, Cons_idx::cp)*dt);
#else
        // Throw everything from qc at qv -> qv will be as saturated as possible
        res[qv_idx] += qc_prime;
        res[qc_idx] -= qc_prime;
        res[T_idx] += T-T_prime;
#endif
#ifdef TRACE_QV
        if (trace)
            std::cout << "traj: " << cc.traj_id << " saturation_adjust (under saturated) dqv " << qc_prime << "\n";
#endif
#ifdef TRACE_QC
        if (trace)
            std::cout << "traj: " << cc.traj_id << " saturation_adjust (under saturated) dqc " << -qc_prime << "\n";
#endif
#ifdef TRACE_ENV
        if (trace)
            std::cout << "traj: " << cc.traj_id << " saturation_adjust (under saturated) dT " << T-T_prime << "\n";
#endif
    }
}


#if defined(CCN_AKM)
/**
 * CCN activation similar to Hande et al 2016.
 * https://doi.org/10.5194/acp-16-12059-2016
 * Adjusted by AKM to fit NAWDEX FAAM PCASP observations
 * using a series of exponential functions.
 * Reimplemented in C++ by MH.
 *
 * @params p_prime Pressure [Pa]
 * @params w_prime Ascend velocity [m s^-1]
 * @params T_prime Temperature [K]
 * @params qv_prime Water vapor mixing ratio
 * @params qc_prime Cloud water mixing ratio
 * @params Nc Number of cloud droplets
 * @params EPSILON Small value that qc needs to exceed
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void ccn_act_hande_akm(
    float_t &p_prime,
    float_t &w_prime,
    float_t &T_prime,
    float_t &qv_prime,
    float_t &qc_prime,
    float_t &Nc,
    const double &EPSILON,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {
    // non maritime case
    if (qc_prime > EPSILON && w_prime > 0.0 && T_prime >= (T_freeze - 38)) {
        float_t bcoeff = get_at(cc.constants, Cons_idx::b_ccn_1)
            * exp(-get_at(cc.constants, Cons_idx::b_ccn_2) * p_prime/100 + get_at(cc.constants, Cons_idx::b_ccn_3));
        float_t ccoeff = get_at(cc.constants, Cons_idx::c_ccn_1)
            * exp(-get_at(cc.constants, Cons_idx::c_ccn_2) * p_prime/100 + get_at(cc.constants, Cons_idx::c_ccn_3));

        float_t Na;
        // The formulation is for hPa, hence dividing by 100.
        if (p_prime < get_at(cc.constants, Cons_idx::p_ccn)) {
            Na = get_at(cc.constants, Cons_idx::h_ccn_1)
                * exp((p_prime/100-get_at(cc.constants, Cons_idx::g_ccn_1))/get_at(cc.constants, Cons_idx::g_ccn_2))
                    + get_at(cc.constants, Cons_idx::h_ccn_2)
                * exp((p_prime/100-get_at(cc.constants, Cons_idx::g_ccn_1))/get_at(cc.constants, Cons_idx::g_ccn_3));
        } else {
            Na = get_at(cc.constants, Cons_idx::h_ccn_1) + get_at(cc.constants, Cons_idx::h_ccn_2);
        }

        // concentration of ccn
        float_t delta_n = get_at(cc.constants, Cons_idx::hande_ccn_fac) * (
            Na * (1/(1+exp(-bcoeff*log(w_prime)-ccoeff))) * get_at(cc.constants, Cons_idx::i_ccn_1));
        delta_n = max(max(delta_n, get_at(cc.constants, Cons_idx::i_ccn_2)) - Nc, 0);
        // Problem here: delta_n shall be the difference between *should be*
        // nuc_n and *is* Nc but this needs to be scaled with our timestep
        // if the timestep is bigger than 1s or else we overshoot
        // The same applies to other cases like this, however this introduces minor numerical errors.
        delta_n /= dt;
        float_t delta_q = min(delta_n * get_at(cc.cloud.constants, Particle_cons_idx::min_x_act),
            res[qv_idx] + qv_prime/dt);
        delta_n = delta_q / get_at(cc.cloud.constants, Particle_cons_idx::min_x_act);

        res[Nc_idx] += delta_n;
        res[qc_idx] += delta_q;
        res[qv_idx] -= delta_q;

#ifdef TRACE_QC
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Ascent dqc " << delta_q << ", dNc " << delta_n
                << ", Nc " << Nc << ", rest " << get_at(cc.constants, Cons_idx::i_ccn_2) << "\n";
#endif
#ifdef TRACE_QV
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Ascent dqv " << -delta_q << "\n";
#endif
        float_t delta_e = latent_heat_evap(T_prime) * delta_q / specific_heat_dry_air(T_prime);
        // Evaporation
        if (delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        else
            res[lat_heat_idx] += delta_e;
    }
}
#else

/**
 * CCN activation after Hande et al 2016.
 * https://doi.org/10.5194/acp-16-12059-2016
 * CCN depends on vertical velocity at different
 * pressure levels.
 * Given the low spatial and the high temporal
 * variability of CCN concentrations, the parameterization
 * is most representative for short time periods.
 * The parameterization consists of a series of arc-tan
 * functions.
 *
 *
 * @params p_prime Pressure [Pa]
 * @params w_prime Ascend velocity [m s^-1]
 * @params T_prime Temperature [K]
 * @params qv_prime Water vapor mixing ratio
 * @params qc_prime Cloud water mixing ratio
 * @params Nc Number of cloud droplets
 * @params EPSILON Small value that qc needs to exceed
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void ccn_act_hande(
    float_t &p_prime,
    float_t &w_prime,
    float_t &T_prime,
    float_t &qv_prime,
    float_t &qc_prime,
    float_t &Nc,
    const double &EPSILON,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {
    // non maritime case
    if (qc_prime > EPSILON && w_prime > 0.0) {
        float_t acoeff = get_at(cc.constants, Cons_idx::a_ccn_1)
            * atan(get_at(cc.constants, Cons_idx::b_ccn_1) * p_prime
                - get_at(cc.constants, Cons_idx::c_ccn_1)) + get_at(cc.constants, Cons_idx::d_ccn_1);
        float_t bcoeff = get_at(cc.constants, Cons_idx::a_ccn_2)
            * atan(get_at(cc.constants, Cons_idx::b_ccn_2) * p_prime
                - get_at(cc.constants, Cons_idx::c_ccn_2)) + get_at(cc.constants, Cons_idx::d_ccn_2);
        float_t ccoeff = get_at(cc.constants, Cons_idx::a_ccn_3)
            * atan(get_at(cc.constants, Cons_idx::b_ccn_3) * p_prime
                - get_at(cc.constants, Cons_idx::c_ccn_3)) + get_at(cc.constants, Cons_idx::d_ccn_3);
        float_t dcoeff = get_at(cc.constants, Cons_idx::a_ccn_4)
            * atan(get_at(cc.constants, Cons_idx::b_ccn_4) * p_prime
                - get_at(cc.constants, Cons_idx::c_ccn_4)) + get_at(cc.constants, Cons_idx::d_ccn_4);
        // concentration of ccn
        float_t nuc_n = acoeff * atan(bcoeff * log(w_prime) + ccoeff) + dcoeff;

        // we need to substract the already "used" ccn in cloud droplets
        // the rest can create new cloud droplets
        // float_t Nc_tmp = qv_prime / get_at(cc.cloud.constants, Particle_cons_idx::min_x);
        float_t delta_n = std::max(nuc_n, float_t(10.0e-6));
        delta_n -= Nc;

        // Problem here: delta_n shall be the difference between *should be*
        // nuc_n and *is* Nc but this needs to be scaled with our timestep
        // if the timestep is bigger than 1s or else we overshoot
        // The same applies to other cases like this, however this introduces minor numerical errors.
        if (dt > 1) delta_n /= dt;
        float_t delta_q = delta_n * get_at(cc.cloud.constants, Particle_cons_idx::min_x_act);
        delta_q = std::min(delta_q, float_t(res[qv_idx] + qv_prime/dt));
        delta_n = std::max(delta_n, float_t(0.0));
        delta_n = delta_q / get_at(cc.cloud.constants, Particle_cons_idx::min_x_act);

        res[Nc_idx] += delta_n;
        res[qc_idx] += delta_q;
        res[qv_idx] -= delta_q;
#ifdef TRACE_QC
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Ascent dqc " << delta_q << ", dNc " << delta_n
                << ", nuc_n " << nuc_n << ", Nc " << Nc << ", rest " << 10.0e-6 << "\n";
#endif
#ifdef TRACE_QV
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Ascent dqv " << -delta_q << "\n";
#endif
        float_t delta_e = latent_heat_evap(T_prime) * delta_q / specific_heat_dry_air(T_prime);
        // Evaporation
        if (delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        else
            res[lat_heat_idx] += delta_e;
    }
}
#endif


/**
 * (optional) homogeneous nucleation using Kaercher et al (2006)
 *
 * @params T_prime Temperature [K]
 * @params w_prime Ascend velocity [m s^-1]
 * @params p_prime Pressure [Pa]
 * @params qv_prime Water vapor mixing ratio
 * @params qi_prime Ice crystal mixing ratio
 * @params Ni Number of ice crystals
 * @params S_i Over saturation regarding to ice
 * @params p_sat_ice Saturation pressure of ice
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void ice_nuc_hom(
    float_t &T_prime,
    float_t &w_prime,
    float_t &p_prime,
    float_t &qv_prime,
    float_t &qi_prime,
    float_t &Ni,
    float_t &S_i,
    float_t &p_sat_ice,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    float_t s_crit = 2.349 - T_prime * (1.0/259.0);
    // mass of water molecule [kg]
    float_t ma_w = get_at(cc.constants, Cons_idx::M_w)/get_at(cc.constants, Cons_idx::N_avo);
    float_t svol = ma_w / get_at(cc.constants, Cons_idx::rho_ice);  // specific volume of a water molecule in ice

    if (S_i > s_crit && T_prime < 235.0 && Ni < get_at(cc.constants, Cons_idx::ni_hom_max)) {
        float_t x_i = particle_mean_mass(qi_prime, Ni,
            get_at(cc.ice.constants, Particle_cons_idx::min_x_nuc_homo),
            get_at(cc.ice.constants, Particle_cons_idx::max_x));
        float_t r_i = pow(x_i/(4.0/3.0*M_PI*get_at(cc.constants, Cons_idx::rho_ice)), 1.0/3.0);
        float_t v_th = sqrt(8.0*get_at(cc.constants, Cons_idx::k_b)*T_prime/(M_PI*ma_w));
        float_t flux = get_at(cc.constants, Cons_idx::alpha_depo) * v_th/4.0;
        float_t n_sat = p_sat_ice/(get_at(cc.constants, Cons_idx::k_b)*T_prime);

        // coeffs of supersaturation equation
        std::vector<float_t> acoeff(3);
        acoeff[0] = (get_at(cc.constants, Cons_idx::L_ed) * get_at(cc.constants, Cons_idx::gravity_acc))
            / (get_at(cc.constants, Cons_idx::cp) * get_at(cc.constants, Cons_idx::R_v) * T_prime*T_prime)
            - get_at(cc.constants, Cons_idx::gravity_acc)/(get_at(cc.constants, Cons_idx::R_a) * T_prime);
        acoeff[1] = 1.0/n_sat;
        acoeff[2] = (get_at(cc.constants, Cons_idx::L_ed)*get_at(cc.constants, Cons_idx::L_ed)
            * get_at(cc.constants, Cons_idx::M_w) * ma_w)
            / (get_at(cc.constants, Cons_idx::cp) * p_prime * T_prime * get_at(cc.constants, Cons_idx::M_a));

        // coeffs of depositional growth equation
        std::vector<float_t> bcoeff(2);
        bcoeff[0] = flux * svol * n_sat * (S_i - 1.0);
        bcoeff[1] = flux/diffusivity(T_prime, p_prime);

        // pre-existing ice crystals included as reduced updraft speed
        float_t ri_dot = bcoeff[0] / (1.0 + bcoeff[1] * r_i);
        float_t R_ik = (4.0*M_PI) / svol * Ni * r_i*r_i * ri_dot;
        float_t w_pre = std::max(float_t(0.0), (acoeff[1] + acoeff[2] * S_i)/(acoeff[0]*S_i)*R_ik);  // Eq. 19

        // homogenous nucleation event
        if (w_prime > w_pre) {
            float_t cool = get_at(cc.constants, Cons_idx::gravity_acc) / get_at(cc.constants, Cons_idx::cp)*w_prime;
            float_t ctau = T_prime * (0.004*T_prime - 2.0) + 304.4;
            float_t tau = 1.0/(ctau*cool);
            float_t delta = bcoeff[1] * get_at(cc.constants, Cons_idx::r_0);
            float_t phi = acoeff[0]*S_i / (acoeff[1]+acoeff[2]*S_i) * (w_prime - w_pre);

            // monodisperse approximation
            float_t kappa = 2.0 * bcoeff[0]*bcoeff[1]*tau/((1.0+delta)*(1.0+delta));
            float_t sqrtkap = sqrt(kappa);
            float_t ren = 3.0*sqrtkap / (2.0 + sqrt(1.0+9.0*kappa/M_PI));
            float_t R_imfc = 4.0 * M_PI * bcoeff[0]/(bcoeff[1]*bcoeff[1]) / svol;
            float_t R_im = R_imfc / (1.0+delta) * (delta*delta - 1.0
                + (1.0+0.5*kappa*(1.0+delta)*(1.0+delta)) * ren/sqrtkap);

            // number concentration and radius of ice particles
            float_t ni_hom = phi/R_im;
            float_t ri_0 = 1.0 + 0.5*sqrtkap * ren;
            float_t ri_hom = (ri_0 * (1.0+delta) - 1.0) / bcoeff[1];
            float_t mi_hom = (4.0/3.0 * M_PI * get_at(cc.constants, Cons_idx::rho_ice))
                * ni_hom * ri_hom*ri_hom*ri_hom;
            mi_hom = std::max(mi_hom, get_at(cc.ice.constants, Particle_cons_idx::min_x_nuc_homo));

            float_t delta_n;
            float_t delta_n_min = get_at(cc.constants, Cons_idx::ni_hom_max)/dt;
            delta_n = std::max(std::min(ni_hom, delta_n_min), float_t(0.0));
            float_t delta_q_min = res[qv_idx] + qv_prime/dt;
            float_t delta_q = std::max(float_t(0), std::min(delta_n * mi_hom, delta_q_min));

            res[Ni_idx] += delta_n;
            res[qi_idx] += delta_q;
            res[qv_idx] -= delta_q;
#ifdef TRACE_QV
            if (trace)
                std::cout << "traj: " << cc.traj_id << " KHL06 dqv " << -delta_q << "\n";
#endif
#ifdef TRACE_QI
            if (trace)
                std::cout << "traj: " << cc.traj_id << " KHL06 dqi " << delta_q << ", dNi " << delta_n << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * delta_q
                / specific_heat_dry_air(T_prime);
            // Sublimation, cooling
            if (delta_e < 0.0)
                res[lat_cool_idx] += delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] += delta_e;
        }
    }
}


/**
 * Heterogeneous nucleation using Phillips et al. (2008)
 * 10.1175/2007JAS2546.1
 * Implementation by Carmen Koehler and AS
 * modified for C++ and Codipack by Maicon Hieronymus
 *
 * @params qc_prime Cloud water mixing ratio
 * @params qv_prime Water vapor mixing ratio
 * @params T_prime Temperature [K]
 * @params S_i Over saturation regarding to ice
 * @params n_inact Number of inactivated nuclei
 * @params use_prog_in Not yet implemented
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void ice_activation_phillips(
    float_t &qc_prime,
    float_t &qv_prime,
    float_t &T_prime,
    float_t &S_i,
    float_t &n_inact,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    const double EPSILON = 1.0e-20;
#ifdef TRACE_QI
    if (trace)
        std::cout << "traj: " << cc.traj_id << " T_prime " << T_prime << "\n"
                  << "T_nuc " << get_at(cc.constants, Cons_idx::T_nuc) << "\n"
                  << "S_i " << S_i << "\n"
                  << "n_inact " << n_inact << "\n"
                  << "get_at(cc.constants, Cons_idx::ni_het_max) "
                  << get_at(cc.constants, Cons_idx::ni_het_max) << "\n";
#endif
    if (T_prime < get_at(cc.constants, Cons_idx::T_nuc) && T_prime > 180.0 && S_i > 1.0
        && n_inact < get_at(cc.constants, Cons_idx::ni_het_max)) {
        float_t x_t = (274.0 - T_prime) / t_tstep;
        x_t = std::min(x_t, float_t(t_tmax-1));
        int tt = static_cast<int>(x_t.getValue()) - 1;

        std::vector<float_t> infrac(3);
        if (qc_prime > EPSILON) {
            // Immersion freezing at water saturation
            infrac[0] = (trunc(x_t)+1.0-x_t) * afrac_dust[98][tt]
                + (x_t-trunc(x_t)) * afrac_dust[98][tt+1];
            infrac[1] = (trunc(x_t)+1.0-x_t) * afrac_soot[98][tt]
                + (x_t-trunc(x_t)) * afrac_soot[98][tt+1];
            infrac[2] = (trunc(x_t)+1.0-x_t) * afrac_orga[98][tt]
                + (x_t-trunc(x_t)) * afrac_orga[98][tt+1];
        } else {
            // deposition nucleation below water saturation
            // Indices for 2D look-up tables
            float_t x_s = 100.0*(S_i-1.0) / s_sstep;
            x_s = std::min(x_s, float_t(s_smax-1));
            int ss = std::max(static_cast<int>(0), static_cast<int>(x_s.getValue()-1));
            float_t S_sr = std::max(float_t(1.0), trunc(x_s));
            infrac[0] = (trunc(x_t)+1.0-x_t) * (S_sr+1.0-x_s)
                * afrac_dust[ss][tt]
                + (x_t-trunc(x_t)) * (S_sr+1.0-x_s)
                * afrac_dust[ss][tt+1]
                + (x_t-trunc(x_t)+1.0-x_t) * (x_s-S_sr)
                * afrac_dust[ss+1][tt]
                + (x_t-trunc(x_t)) * (x_s-S_sr)
                * afrac_dust[ss+1][tt+1];
            infrac[1] = (trunc(x_t)+1.0-x_t) * (S_sr+1.0-x_s)
                * afrac_soot[ss][tt]
                + (x_t-trunc(x_t)) * (S_sr+1.0-x_s)
                * afrac_dust[ss][tt+1]
                + (x_t-trunc(x_t)+1.0-x_t) * (x_s-S_sr)
                * afrac_soot[ss+1][tt]
                + (x_t-trunc(x_t)) * (x_s-S_sr)
                * afrac_soot[ss+1][tt+1];
            infrac[2] = (trunc(x_t)+1.0-x_t) * (S_sr+1.0-x_s)
                * afrac_orga[ss][tt]
                + (x_t-trunc(x_t)) * (S_sr+1.0-x_s)
                * afrac_orga[ss][tt+1]
                + (x_t-trunc(x_t)+1.0-x_t) * (x_s-S_sr)
                * afrac_orga[ss+1][tt]
                + (x_t-trunc(x_t)) * (x_s-S_sr)
                * afrac_dust[ss+1][tt+1];
        }
        float_t ndiag = infrac[0]*get_at(cc.constants, Cons_idx::na_dust)
            + get_at(cc.constants, Cons_idx::na_soot)*infrac[1]
            + get_at(cc.constants, Cons_idx::na_orga)*infrac[2];
#ifdef TRACE_QI
        if (trace)
            std::cout << "traj: " << cc.traj_id << " \n\nqc>EPSILON: " << (qc_prime > EPSILON)
                      << "\ninfac[0]: " << infrac[0]
                      << "\ninfrac[1]: " << infrac[1]
                      << "\ninfrac[2]: " << infrac[2]
                      << "\nna_soot: " << get_at(cc.constants, Cons_idx::na_soot)
                      << "\nna_orga: " << get_at(cc.constants, Cons_idx::na_orga)
                      << "\ntt: " << tt
                      << "\nafrac_dust[98][tt]: " << afrac_dust[98][tt]
                      << "\nqv: " << qv_prime
                      << "\nqc: " << qc_prime
                      << "\nT: " << T_prime
                      << "\nn_inact: " << n_inact
                      << "\nS_i: " << S_i
                      << "\nndiag: " << ndiag
                      << "\nni_het_max: " << get_at(cc.constants, Cons_idx::ni_het_max)
                      << "\nn_inact: " << n_inact
                      << "\nna_dust: " << get_at(cc.constants, Cons_idx::na_dust) << "\n";
#endif

        float_t tmp_ndiag = get_at(cc.constants, Cons_idx::ni_het_max)/dt;
        ndiag = std::min(ndiag, tmp_ndiag);
        float_t delta_n = std::max(ndiag-n_inact, float_t(0.0))/dt;
        float_t delta_q_max = res[qv_idx] + qv_prime/dt;
        float_t delta_q =
            std::max(float_t(0), std::min(delta_n*get_at(cc.ice.constants, Particle_cons_idx::min_x_act), delta_q_max));

        delta_n = delta_q/get_at(cc.ice.constants, Particle_cons_idx::min_x_act);
        res[Ni_idx] += delta_n;
        res[qi_idx] += delta_q;
        res[qv_idx] -= delta_q;
        res[n_inact_idx] += delta_n;
        res[depo_idx] += delta_q;
#ifdef TRACE_QV
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Phillips nucleation dqv: " << -delta_q << "\n";
#endif
#ifdef TRACE_QI
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " Phillips nucleation dqi: " << delta_q << ", dNi: " << delta_n << "\n";
#endif
        // latent heating and cooling
        float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * delta_q
            / specific_heat_dry_air(T_prime);
        // Sublimation, cooling
        if (delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] += delta_e;
    }
}

/**
 * Homogeneous freezing of cloud droplets. Cotton and Field (2002), Seifert & Beheng (2006) for reference
 * Ice nucleation characteristics of an isolated wave cloud
 * Freeze only if cloud droplets to freeze are available
 * and the temperature is low enough.
 *
 * @params qc_prime Cloud water mixing ratio
 * @params Nc Number of cloud droplets
 * @params T_prime Temperature [K]
 * @params T_c Difference between temperature and freezing temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void cloud_freeze_hom(
    float_t &qc_prime,
    float_t &Nc,
    float_t &T_prime,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (qc_prime > 0.0 && T_c < -30.0) {
        float_t x_c = particle_mean_mass(
            qc_prime, Nc, get_at(cc.cloud.constants, Particle_cons_idx::min_x_freezing),
            get_at(cc.cloud.constants, Particle_cons_idx::max_x));

        float_t delta_qi;
        float_t delta_ni;
        // instantaneous freezing for temperatures below -50 °C
        if (T_c < -50.0) {
            delta_qi = res[qc_idx] + qc_prime/dt;
            delta_ni = res[Nc_idx] + Nc/dt;
        } else {
            float_t j_hom;
            if (T_c > -30.0)
                j_hom = 1.0e6 / get_at(cc.constants, Cons_idx::rho_w) * pow(10,
                    -7.63-2.996*(T_c+30.0));
            else
                j_hom = 1.0e6 / get_at(cc.constants, Cons_idx::rho_w) * pow(10,
                    - 243.4
                    - 14.75 * T_c
                    - 0.307 * T_c * T_c
                    - 0.00287 * T_c * T_c * T_c
                    - 0.0000102 * pow(T_c, 4));
            delta_ni = j_hom * (res[qc_idx] + qc_prime/dt);
            delta_qi = j_hom * (res[qc_idx] + qc_prime/dt)
                * x_c * get_at(cc.cloud.constants, Particle_cons_idx::c_z);
#ifdef TRACE_QI
            if (trace)
                std::cout << "rho_w: " << get_at(cc.constants, Cons_idx::rho_w) << "\n"
                    << "x_c: " << x_c << "\n"
                    << "c_z: " << get_at(cc.cloud.constants, Particle_cons_idx::c_z) << "\n"
                    << "j_hom: " << j_hom << "\n"
                    << "res[qc]: " << res[qc_idx] << "\n"
                    << "T_c: " << T_c << "\n"
                    << "T_prime: " << T_prime << "\n"
                    << "qc_prime: " << qc_prime << "\n";
#endif
            float_t delta_min = res[Nc_idx] + Nc/dt;
            delta_ni = std::max(float_t(0), std::min(delta_ni, delta_min));
            delta_min = res[qc_idx] + qc_prime/dt;
            delta_qi = std::max(float_t(0), std::min(delta_qi, delta_min));
        }
        // Remove cloud droplets
        res[qc_idx] -= delta_qi;
        res[Nc_idx] -= delta_ni;
#ifdef TRACE_QC
        if (trace)
            std::cout << "traj: " << cc.traj_id << " cloud freeze dqc " << -delta_qi << ", dNc " << -delta_ni << "\n";
#endif
        // The amount of ice crystals should be capped by the maximum size
        // of cloud droplets since big cloud droplets are rain droplets per definition...
        float_t max_delta_ni = delta_qi/get_at(cc.cloud.constants, Particle_cons_idx::max_x);
        delta_ni = std::max(delta_ni, max_delta_ni);

        res[qi_idx] += delta_qi;
        res[Ni_idx] += delta_ni;
#ifdef TRACE_QI
        if (trace)
            std::cout << "traj: " << cc.traj_id << " cloud freeze dqi " << delta_qi << ", dNi " << delta_ni << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * delta_qi
            / specific_heat_dry_air(T_prime);
        // Melting, cooling
        if (delta_qi < 0.0)
            res[lat_cool_idx] += delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] += delta_e;
    }
}


/**
 * Ice-ice collection. Seifert & Beheng (2006), Eq. 61-67
 *
 * @params qi_prime Ice crystal mixing ratio
 * @params Ni Number of ice crystals
 * @params Difference between temperature and freezing temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void ice_self_collection(
    float_t &qi_prime,
    float_t &Ni,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    float_t x_i = particle_mean_mass(qi_prime, Ni,
        get_at(cc.ice.constants, Particle_cons_idx::min_x_collision),
        get_at(cc.ice.constants, Particle_cons_idx::max_x));
    float_t D_i = particle_diameter(x_i,
        get_at(cc.ice.constants, Particle_cons_idx::a_geo),
        get_at(cc.ice.constants, Particle_cons_idx::b_geo));
    if (Ni > 0.0 && qi_prime > get_at(cc.constants, Cons_idx::q_crit_i)
        && D_i > get_at(cc.constants, Cons_idx::D_crit_i)) {
        float_t x_conv_i = pow(get_at(cc.constants, Cons_idx::D_conv_i)
            / get_at(cc.snow.constants, Particle_cons_idx::a_geo),
            1.0/get_at(cc.snow.constants, Particle_cons_idx::b_geo));
        // efficiency depends on temperature here (Cotton et al 1986)
        // also Straka 1989, page 53
        float_t e_coll = std::min(pow(10, 0.035*T_c-0.7), 0.2);
        float_t vel_i = particle_velocity(x_i,
            get_at(cc.ice.constants, Particle_cons_idx::a_vel),
            get_at(cc.ice.constants, Particle_cons_idx::b_vel))
            * get_at(cc.ice.constants, Particle_cons_idx::rho_v);

        float_t delta_n = M_PI/4.0 * e_coll * get_at(cc.ice.constants, Particle_cons_idx::sc_delta_n)
            * Ni * Ni * D_i * D_i * sqrt(
                get_at(cc.ice.constants, Particle_cons_idx::sc_theta_n) * vel_i * vel_i
                + 2.0 * get_at(cc.ice.constants, Particle_cons_idx::s_vel)
                * get_at(cc.ice.constants, Particle_cons_idx::s_vel));
        float_t delta_q = M_PI/4.0 * e_coll * get_at(cc.ice.constants, Particle_cons_idx::sc_delta_q)
            * Ni * qi_prime * D_i * D_i * sqrt(
                get_at(cc.ice.constants, Particle_cons_idx::sc_theta_q) * vel_i * vel_i
                + 2.0 * get_at(cc.ice.constants, Particle_cons_idx::s_vel)
                * get_at(cc.ice.constants, Particle_cons_idx::s_vel));

        float_t delta_min = res[qi_idx] + qi_prime/dt;
        delta_q = std::max(float_t(0), std::min(delta_q, delta_min));
        delta_min = res[Ni_idx] + Ni/dt;
        float_t delta_min2 = delta_q/x_conv_i;
        delta_n = std::max(float_t(0), std::min(std::min(delta_n, delta_min2), delta_min));

        res[qi_idx] -= delta_q;
        res[qs_idx] += delta_q;
        res[Ni_idx] -= delta_n;
        res[Ns_idx] += delta_n/2.0;
#ifdef TRACE_QI
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " ice self colletion dqi " << -delta_q << ", dNi " << -delta_n << "\n";
#endif
#ifdef TRACE_QS
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " ice self colletion dqs " << delta_q << ", dNs " << delta_n/2.0 << "\n";
#endif
    }
}

/**
 * Snow-snow collection after Seifert & Beheng (2006).
 * 10.1007/s00703-005-0112-4
 *
 * @params qs_prime Snowflakes mixing ratio
 * @params Ns Number of snowflakes
 * @params T_prime Temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 */
template<class float_t>
void snow_self_collection(
    float_t &qs_prime,
    float_t &Ns,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc) {

    if (qs_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        // temperature dependent sticking efficiency Lin (1983)
        float_t e_coll = std::max(0.1, std::min(exp(0.09*(T_prime-get_at(cc.constants, Cons_idx::T_freeze))), 1.0));
        float_t x_s = particle_mean_mass(qs_prime, Ns, get_at(cc.snow.constants, Particle_cons_idx::min_x_collection),
            get_at(cc.snow.constants, Particle_cons_idx::max_x));
        float_t D_s = particle_diameter(x_s, get_at(cc.snow.constants, Particle_cons_idx::a_geo),
            get_at(cc.snow.constants, Particle_cons_idx::b_geo));
        float_t vel_s = particle_velocity(x_s, get_at(cc.snow.constants, Particle_cons_idx::a_vel),
            get_at(cc.snow.constants, Particle_cons_idx::b_vel)) * get_at(cc.snow.constants, Particle_cons_idx::rho_v);

        res[Ns_idx] -= M_PI/8.0 * e_coll * Ns * Ns * get_at(cc.snow.constants, Particle_cons_idx::sc_delta_n)
            * D_s * D_s * sqrt(get_at(cc.snow.constants, Particle_cons_idx::sc_theta_n) * vel_s * vel_s
            + 2.0 * get_at(cc.snow.constants, Particle_cons_idx::s_vel)
            * get_at(cc.snow.constants, Particle_cons_idx::s_vel));
#ifdef TRACE_QS
        if (trace)
            std::cout << "traj: " << cc.traj_id << " snow self collection dNs " << - (M_PI/8.0 * e_coll * Ns * Ns
                * get_at(cc.snow.constants, Particle_cons_idx::sc_delta_n) * D_s * D_s
                * sqrt(get_at(cc.snow.constants, Particle_cons_idx::sc_theta_n) * vel_s * vel_s + 2.0
                * get_at(cc.snow.constants, Particle_cons_idx::s_vel)
                * get_at(cc.snow.constants, Particle_cons_idx::s_vel))) << "\n";
#endif
    }
}


/**
 * Melting of snow for temperatures above freezing temperature after Seifert & Beheng (2006),
 * Eqs. 72-77, and Eqs. 85-89.
 * 10.1007/s00703-005-0112-4
 *
 * @params qs_prime Snowflakes mixing ratio
 * @params Ns Number of snowflakes
 * @params T_prime Temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void snow_melting(
    float_t &qs_prime,
    float_t &Ns,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (T_prime > get_at(cc.constants, Cons_idx::T_freeze) && qs_prime > 0.0) {
        float_t p_sat = saturation_pressure_water(
            T_prime, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));

        float_t x_s = particle_mean_mass(qs_prime, Ns, get_at(cc.snow.constants, Particle_cons_idx::min_x_melt),
            get_at(cc.snow.constants, Particle_cons_idx::max_x));
        float_t d_s = particle_diameter(x_s, get_at(cc.snow.constants, Particle_cons_idx::a_geo),
            get_at(cc.snow.constants, Particle_cons_idx::b_geo));
        float_t v_s = particle_velocity(x_s, get_at(cc.snow.constants, Particle_cons_idx::a_vel),
            get_at(cc.snow.constants, Particle_cons_idx::b_vel)) * get_at(cc.snow.constants, Particle_cons_idx::rho_v);

        float_t fv_q = get_at(cc.snow.constants, Particle_cons_idx::a_f)
            + get_at(cc.snow.constants, Particle_cons_idx::b_f) * sqrt(v_s*d_s);
        // From Rasmussen and Heymsfield (1987)
        float_t fh_q = 1.05 * fv_q;
        float_t melt = 2.0*M_PI/L_ew * d_s * Ns;
        float_t melt_h = melt * get_at(cc.constants, Cons_idx::K_T)
            * (T_prime - get_at(cc.constants, Cons_idx::T_freeze));
        float_t melt_v = melt * get_at(cc.constants, Cons_idx::dv0)
            * get_at(cc.constants, Cons_idx::L_wd)/get_at(cc.constants, Cons_idx::R_v)
            * (p_sat/T_prime - get_at(cc.constants, Cons_idx::p_sat_melt)/get_at(cc.constants, Cons_idx::T_freeze));
        float_t melt_q = (melt_h * fh_q + melt_v * fv_q);
        float_t melt_n = res[Ns_idx] + Ns/dt;
        melt_n = std::min(std::max((melt_q-qs_prime)/x_s + Ns, float_t(0.0)), melt_n);
        float_t tmp_melt = res[qs_idx] + qs_prime/dt;
        melt_q = std::min(tmp_melt, std::max(melt_q, float_t(0.0)));
        tmp_melt = res[Ns_idx] + Ns/dt;
        melt_n = std::min(tmp_melt, std::max(melt_n, float_t(0.0)));
        if (T_prime - get_at(cc.constants, Cons_idx::T_freeze) > 10.0) {
            melt_q = qs_prime/dt;
            melt_n = Ns/dt;
        }

        // Snow
        res[qs_idx] -= melt_q;
        // Snow N
        res[Ns_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " snow melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
#ifdef TRACE_QS
        if (trace)
            std::cout << "traj: " << cc.traj_id << " snow melting dqs " << -melt_q << ", dNs " << -melt_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * melt_q
            / specific_heat_dry_air(T_prime);
        // Melting, cooling
        if (melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 * Melting of snow for temperatures above freezing temperature after Seifert & Beheng (2006),
 * Eqs. 72-77, and Eqs. 85-89.
 * 10.1007/s00703-005-0112-4
 *
 * @params qg_prime Graupel mixing ratio
 * @params Ng Number of graupel particles
 * @params T_prime Temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void graupel_melting(
    float_t &qg_prime,
    float_t &Ng,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (T_prime > get_at(cc.constants, Cons_idx::T_freeze) && qg_prime > 0.0) {
        float_t p_sat = saturation_pressure_water(
            T_prime, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));
        float_t x_g = particle_mean_mass(qg_prime, Ng, get_at(cc.graupel.constants, Particle_cons_idx::min_x_melt),
            get_at(cc.graupel.constants, Particle_cons_idx::max_x));
        float_t d_g = particle_diameter(x_g, get_at(cc.graupel.constants, Particle_cons_idx::a_geo),
            get_at(cc.graupel.constants, Particle_cons_idx::b_geo));
        float_t v_g = particle_velocity(x_g, get_at(cc.graupel.constants, Particle_cons_idx::a_vel),
            get_at(cc.graupel.constants, Particle_cons_idx::b_vel))
            * get_at(cc.graupel.constants, Particle_cons_idx::rho_v);

        float_t fv_q = get_at(cc.graupel.constants, Particle_cons_idx::a_f)
            + get_at(cc.graupel.constants, Particle_cons_idx::b_f) * sqrt(v_g*d_g);
        float_t fh_q = 1.05 * fv_q;
        float_t melt = 2.0*M_PI/L_ew * d_g * Ng;
        float_t melt_h = melt * get_at(cc.constants, Cons_idx::K_T)
            * (T_prime - get_at(cc.constants, Cons_idx::T_freeze));
        float_t melt_v = melt * get_at(cc.constants, Cons_idx::dv0)
            * get_at(cc.constants, Cons_idx::L_wd)/get_at(cc.constants, Cons_idx::R_v)
            * (p_sat/T_prime - get_at(cc.constants, Cons_idx::p_sat_melt)/get_at(cc.constants, Cons_idx::T_freeze));
        float_t melt_q = (melt_h * fh_q + melt_v * fv_q);
        float_t melt_tmp = res[Ng_idx] + Ng/dt;
        float_t melt_n = std::min(std::max((melt_q-qg_prime)/x_g + Ng, float_t(0.0)), melt_tmp);

        melt_q = std::max(float_t(0.0), melt_q);
        melt_n = std::max(float_t(0.0), melt_n);

        melt_tmp = res[qg_idx] + qg_prime/dt;
        melt_q = std::max(float_t(0.0), std::min(melt_q, melt_tmp));
        melt_tmp = res[Ng_idx] + Ng/dt;
        melt_n = std::max(float_t(0.0), std::max(melt_n, melt_tmp));

        // Graupel
        res[qg_idx] -= melt_q;
        // Graupel N
        res[Ng_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " graupel melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
#ifdef TRACE_QG
        if (trace)
            std::cout << "traj: " << cc.traj_id << " graupel melting dqg " << -melt_q << ", dNg " << -melt_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * melt_q
            / specific_heat_dry_air(T_prime);
        // Melting, cooling
        if (melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 * Melting of snow for temperatures above freezing temperature after Seifert & Beheng (2006),
 * Eqs. 72-77, and Eqs. 85-89.
 * 10.1007/s00703-005-0112-4
 *
 * @params qh_prime Hail mixing ratio
 * @params Nh Number of hail particles
 * @params T_prime Temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void hail_melting(
    float_t &qh_prime,
    float_t &Nh,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (T_prime > get_at(cc.constants, Cons_idx::T_freeze) && qh_prime > 0.0) {
        float_t p_sat = saturation_pressure_water(
            T_prime, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));
        float_t x_h = particle_mean_mass(qh_prime, Nh, get_at(cc.hail.constants, Particle_cons_idx::min_x_melt),
            get_at(cc.hail.constants, Particle_cons_idx::max_x));
        float_t d_h = particle_diameter(x_h, get_at(cc.hail.constants, Particle_cons_idx::a_geo),
            get_at(cc.hail.constants, Particle_cons_idx::b_geo));
        float_t v_h = particle_velocity(x_h, get_at(cc.hail.constants, Particle_cons_idx::a_vel),
            get_at(cc.hail.constants, Particle_cons_idx::b_vel)) * get_at(cc.hail.constants, Particle_cons_idx::rho_v);

        float_t fv_q = get_at(cc.hail.constants, Particle_cons_idx::a_f)
            + get_at(cc.hail.constants, Particle_cons_idx::b_f) * sqrt(v_h*d_h);
        float_t fh_q = 1.05 * fv_q;
        float_t melt = 2.0*M_PI/L_ew * d_h * Nh;
        float_t melt_h = melt * get_at(cc.constants, Cons_idx::K_T)
            * (T_prime - get_at(cc.constants, Cons_idx::T_freeze));
        float_t melt_v = melt * get_at(cc.constants, Cons_idx::dv0)
            * get_at(cc.constants, Cons_idx::L_wd)/get_at(cc.constants, Cons_idx::R_v)
            * (p_sat/T_prime - get_at(cc.constants, Cons_idx::p_sat_melt)/get_at(cc.constants, Cons_idx::T_freeze));
        float_t melt_q = (melt_h * fh_q + melt_v * fv_q);
        float_t melt_tmp = res[Nh_idx] + Nh/dt;
        float_t melt_n = std::min(std::max((melt_q-qh_prime)/x_h + Nh, float_t(0.0)), melt_tmp);

        melt_tmp = res[qh_idx] + qh_prime/dt;
#ifdef TRACE_MELT
        if (trace)
            std::cout << "Melt melt_q " << melt_q
                      << "\nqh_prime " << qh_prime
                      << "\nmelt_tmp " << melt_tmp
                      << "\nmin: " << std::min(melt_q, melt_tmp)
                      << "\nmelt_q after " << std::max(float_t(0.0), std::min(melt_q, melt_tmp))
                      << "\n";
#endif
//        melt_q = std::max(float_t(0.0), melt_q);
        melt_q = std::max(float_t(0.0), std::min(melt_q, melt_tmp));
        melt_tmp = res[Nh_idx] + Nh/dt;
        melt_n = std::max(float_t(0.0), std::max(melt_n, melt_tmp));

        // Hail
        res[qh_idx] -= melt_q;
        // Hail N
        res[Nh_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " hail melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
#ifdef TRACE_QH
        if (trace)
            std::cout << "traj: " << cc.traj_id << " hail melting dqh " << -melt_q << ", dNh " << -melt_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * melt_q
            / specific_heat_dry_air(T_prime);
        // Melting, cooling
        if (melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 * Conversion of cloud droplets to rain droplets after Seifert & Beheng (2006).
 * 10.1007/s00703-005-0112-4
 *
 * @params qc_prime Cloud mixing ratio
 * @params Nc Number of cloud droplets
 * @params T_prime Temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void auto_conversion_kb(
    float_t &qc_prime,
    float_t &Nc,
    float_t &qr_prime,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    // autoconversionKB
    float_t k_a = 6.0 + 25 * pow(9.59, -1.7);
    float_t x_s_i = 1.0/get_at(cc.cloud.constants, Particle_cons_idx::max_x);
    float_t x_c = particle_mean_mass(
        qc_prime, Nc, get_at(cc.cloud.constants, Particle_cons_idx::min_x_conversion),
        get_at(cc.cloud.constants, Particle_cons_idx::max_x));
    // Using Beheng 1994
    float_t au = k_a * pow(x_c*1e3, 3.3) * pow(qc_prime*1e3, 1.4) * 1e3;
    float_t tmp = res[qc_idx] + qc_prime/dt;
    au = std::min(tmp, au);
    res[Nr_idx] += au*x_s_i;
    res[qr_idx] += au;
    res[Nc_idx] -= au*x_s_i*2.0;
    res[qc_idx] -= au;
#ifdef TRACE_QC
    if (trace)
        if (abs(au) > 0)
            std::cout << "traj: " << cc.traj_id << " autoconversion dqc " << -au << ", dNc " << -au*x_s_i*2.0 << "\n";
#endif
#ifdef TRACE_QR
    if (trace)
        std::cout << "traj: " << cc.traj_id << " autoconversion dqr " << au << ", dNr " << au*x_s_i << "\n";
#endif
    // accretionKB
    if (qc_prime > get_at(cc.constants, Cons_idx::q_crit_i) && qr_prime > get_at(cc.constants, Cons_idx::q_crit_i)) {
        // k_r = 6.0 from Beheng (1994)
        float_t ac = 6.0 * qc_prime * qr_prime;
        float_t ac_tmp = res[qc_idx] + qc_prime/dt;
        ac = std::min(ac_tmp, ac);
        res[qr_idx] += ac;
        res[qc_idx] -= ac;
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " accretionKB dqr " << ac << "\n";
#endif
#ifdef TRACE_QC
        if (trace)
            if (abs(ac) > 0)
                std::cout << "traj: " << cc.traj_id << " accretionKB dqc " << -ac << "\n";
#endif
    }
}

/**
 * Formation of raindrops by coagulating cloud droplets and growth
 * of raindrops collecting cloud droplets using Seifert and Beheng (2006)
 * section 2.2.1.
 * 10.1007/s00703-005-0112-4
 *
 * @params qc_prime Cloud mixing ratio
 * @params Nc Number of cloud droplets
 * @params qr_prime Rain mixing ratio
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void auto_conversion_sb(
    float_t &qc_prime,
    float_t &Nc,
    float_t &qr_prime,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    const double EPSILON = 1.0e-25;

    if (qc_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        float_t x_c = particle_mean_mass(
            qc_prime, Nc, get_at(cc.cloud.constants, Particle_cons_idx::min_x_conversion),
            get_at(cc.cloud.constants, Particle_cons_idx::max_x));
        float_t au = get_at(cc.constants, Cons_idx::cloud_k_au) * qc_prime*qc_prime
                                * x_c*x_c * get_at(cc.cloud.constants, Particle_cons_idx::rho_v);
        float_t tau = std::min(std::max(1.0-qc_prime/
                                (qc_prime+qr_prime+EPSILON), EPSILON), 0.9);
        float_t phi = get_at(cc.constants, Cons_idx::k_1_conv)
            * pow(tau, get_at(cc.constants, Cons_idx::k_2_conv))
            * pow(1.0-pow(tau, get_at(cc.constants, Cons_idx::k_2_conv)), 3);
        au *= (1.0 + phi/pow(1.0-tau, 2));
        float_t tmp = res[qc_idx] + qc_prime/dt;
        au = std::max(std::min(tmp, au), float_t(0.0));

        float_t sc = get_at(cc.constants, Cons_idx::cloud_k_sc) * qc_prime*qc_prime
            * get_at(cc.cloud.constants, Particle_cons_idx::rho_v);

        res[qr_idx] += au;
        res[Nr_idx] += au / get_at(cc.cloud.constants, Particle_cons_idx::max_x);
        res[qc_idx] -= au;
        tmp = res[Nc_idx] + Nc/dt;
        res[Nc_idx] -= std::min(tmp, sc);
#ifdef TRACE_QC
        if (trace)
            if (abs(au) > 0)
                std::cout << "traj: " << cc.traj_id
                    << " autoconversion dqc " << -au << ", dNc " << -min(Nc, sc) << "\n";
#endif
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " autoconversion dqr " << au
                      << ", dNr " << au / get_at(cc.cloud.constants, Particle_cons_idx::max_x) << "\n";
#endif
    }

    // accretion
    if (qc_prime > 0.0 && qr_prime > 0.0) {
        float_t tau = std::min(std::max(1.0-qc_prime/
                                (qc_prime+qr_prime+EPSILON), EPSILON), 1.0);
        float_t phi = pow(tau/(tau+get_at(cc.constants, Cons_idx::k_1_accr)), 4);
        float_t ac = get_at(cc.constants, Cons_idx::k_r) * qc_prime * qr_prime * phi;
        float_t tmp = res[qc_idx] + qc_prime/dt;
        ac = std::min(tmp, ac);
        float_t x_c = particle_mean_mass(
            qc_prime, Nc, get_at(cc.cloud.constants, Particle_cons_idx::min_x_conversion),
            get_at(cc.cloud.constants, Particle_cons_idx::max_x));
        res[qr_idx] += ac;
        res[qc_idx] -= ac;
        tmp = res[Nc_idx] + Nc/dt;
        res[Nc_idx] -= std::min(tmp, x_c);
#ifdef TRACE_QC
        if (trace)
            std::cout << "traj: " << cc.traj_id << " accretionSB dqc " << -ac << ", dNc " << -min(Nc, x_c) << "\n";
#endif
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " accretionSB dqr " << ac << "\n";
#endif
    }
}


/**
 * Rain self collection after Seifert and Beheng (2008), includes breakup
 * of droplets that are bigger than 0.3e-3.
 * https://www.imk-tro.kit.edu/4437_1388.php
 *
 * @params qr_prime Rain mixing ratio
 * @params Nr Number of rain droplets
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void rain_self_collection_sb(
    float_t &qr_prime,
    float_t &Nr,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (qr_prime > 0) {
        float_t x_r = particle_mean_mass(qr_prime, Nr, get_at(cc.rain.constants, Particle_cons_idx::min_x_collection),
            get_at(cc.rain.constants, Particle_cons_idx::max_x));
        float_t D_r = particle_diameter(x_r, get_at(cc.rain.constants, Particle_cons_idx::a_geo),
            get_at(cc.rain.constants, Particle_cons_idx::b_geo));
        // Parameters based on Seifert (2008)
        float_t sc = 4.33 * Nr * qr_prime * get_at(cc.rain.constants, Particle_cons_idx::rho_v);
        // Breakup Seifert (2008), Eq. A13
        float_t breakup = 0.0;
        if (D_r > get_at(cc.constants, Cons_idx::D_br_threshold))
            breakup = sc * (get_at(cc.constants, Cons_idx::k_br)
                * (D_r - get_at(cc.constants, Cons_idx::D_br)) + get_at(cc.constants, Cons_idx::c_br));
        float_t tmp = res[Nr_idx] + Nr/dt;
        res[Nr_idx] -= std::min(tmp, sc-breakup);
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " self collection dNr " << -min(Nr, sc-breakup) << "\n";
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Nr " << Nr
                      << "\nsc: " << sc
                      << "\nbreakup: " << breakup
                      << "\nrain.constants[Particle_cons_idx::rho_v]: "
                      << get_at(cc.rain.constants, Particle_cons_idx::rho_v) << "\n";
#endif
    }
}


/**
 * Rain evaporation after Seifert (2008)
 * 10.1175/2008JAS2586.1
 *
 * @params qr_prime Rain mixing ratio
 * @params Nr Number of rain droplets
 * @params qv_prime Water vapor mixing ratio
 * @params qc_prime Cloud water mixing ratio
 * @params T_prime Temperature [K]
 * @params p_prime Pressure [Pa]
 * @params s_sw Over saturation of water
 * @params p_sat Saturation pressure of water [Pa]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void rain_evaporation_sb(
    float_t &qr_prime,
    float_t &Nr,
    float_t &qc_prime,
    float_t &T_prime,
    float_t &p_prime,
    float_t &s_sw,
    float_t &p_sat,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (s_sw < 0.0 && qr_prime > 0.0 && qc_prime < get_at(cc.constants, Cons_idx::q_crit)) {
        float_t D_v_local = diffusivity(T_prime, p_prime);
        // Equation (A2) of Seifert (2008)
        float_t g_d = 2*M_PI /
            ((get_at(cc.constants, Cons_idx::R_v)*T_prime)/(D_v_local*p_sat)
            + (get_at(cc.constants, Cons_idx::L_wd)*get_at(cc.constants, Cons_idx::L_wd))
            / (get_at(cc.constants, Cons_idx::K_T)*get_at(cc.constants, Cons_idx::R_v)*T_prime*T_prime));
        float_t x_r = particle_mean_mass(qr_prime, Nr, get_at(cc.rain.constants, Particle_cons_idx::min_x_evap),
            get_at(cc.rain.constants, Particle_cons_idx::max_x));
        float_t D_r = particle_diameter(x_r, get_at(cc.rain.constants, Particle_cons_idx::a_geo),
            get_at(cc.rain.constants, Particle_cons_idx::b_geo));

        float_t mue;
        // Equation 20 of Seifert (2008)
        if (D_r <= get_at(cc.rain.constants, Particle_cons_idx::cmu3))
            mue = get_at(cc.rain.constants, Particle_cons_idx::cmu0)
                * tanh(pow(4.0*get_at(cc.rain.constants, Particle_cons_idx::cmu2)
                * (D_r-get_at(cc.rain.constants, Particle_cons_idx::cmu3)),
                    get_at(cc.rain.constants, Particle_cons_idx::cmu5)))
                    + get_at(cc.rain.constants, Particle_cons_idx::cmu4);
        else
            mue = get_at(cc.rain.constants, Particle_cons_idx::cmu1)
                * tanh(pow(get_at(cc.rain.constants, Particle_cons_idx::cmu2)
                * (D_r-get_at(cc.rain.constants, Particle_cons_idx::cmu3)),
                    get_at(cc.rain.constants, Particle_cons_idx::cmu5)))
                    + get_at(cc.rain.constants, Particle_cons_idx::cmu4);
        // Equation A8
        float_t lambda = pow(
            M_PI/6.0*get_at(cc.constants, Cons_idx::rho_w)*(mue+3.0)*(mue+2.0)*(mue+1.0)/x_r, 1.0/3.0);

        // Approximation of Gamma(mue+5/2) / Gamma(mue+2)
        float_t gamma_approx = 0.1357940435E+01
            + mue * (0.3033273220E+00
            + mue * (-0.1299313363E-01
            + mue * (0.4002257774E-03
            - mue * 0.4856703981E-05)));
        // ventilation factor (A7) with (A5) and (A9)
        // Approximation for terminal fall velocity of raindrops
        float_t f_v = get_at(cc.constants, Cons_idx::a_v) + get_at(cc.constants, Cons_idx::b_v) * pow(N_Sc, 1.0/3.0)
            * gamma_approx
            * sqrt(get_at(cc.rain.constants, Particle_cons_idx::alpha)/get_at(cc.constants, Cons_idx::kin_visc_air)
                * get_at(cc.rain.constants, Particle_cons_idx::rho_v)/lambda)
            * (1.0
                - 0.5 * (get_at(cc.rain.constants, Particle_cons_idx::beta)
                    /get_at(cc.rain.constants, Particle_cons_idx::alpha))
                    * pow(lambda/(get_at(cc.rain.constants, Particle_cons_idx::gamma)+lambda), (mue+2.5))
                - 0.125 * pow(get_at(cc.rain.constants, Particle_cons_idx::beta)
                    /get_at(cc.rain.constants, Particle_cons_idx::alpha), 2.0)
                    * pow(lambda/(2*get_at(cc.rain.constants, Particle_cons_idx::gamma)+lambda), (mue+2.5))
                - 1.0/16.0 * pow(get_at(cc.rain.constants, Particle_cons_idx::beta)
                    /get_at(cc.rain.constants, Particle_cons_idx::alpha), 3.0)
                    * pow(lambda/(3*get_at(cc.rain.constants, Particle_cons_idx::gamma)+lambda), (mue+2.5))
                - 5.0/127.0 * pow(get_at(cc.rain.constants, Particle_cons_idx::beta)
                    /get_at(cc.rain.constants, Particle_cons_idx::alpha), 4.0)
                    * pow(lambda/(4*get_at(cc.rain.constants, Particle_cons_idx::gamma)+lambda), (mue+2.5)));
        float_t gamma_eva;

        if (get_at(cc.constants, Cons_idx::rain_gfak) > 0.0)
            gamma_eva = get_at(cc.constants, Cons_idx::rain_gfak) * (1.1e-3/D_r) * exp(-0.2*mue);
        else
            gamma_eva = 1.0;

        // Equation A5 with A9
        float_t delta_qv = g_d * Nr * (mue+1.0) / lambda * f_v * s_sw;
        float_t delta_nv = gamma_eva * delta_qv/x_r;

        delta_qv = std::max(-delta_qv, float_t(0.0));
        delta_nv = std::max(-delta_nv, float_t(0.0));
        float_t delta_tmp = res[qr_idx] + qr_prime/dt;
        delta_qv = std::min(delta_qv, delta_tmp);
        delta_tmp = res[Nr_idx] + Nr/dt;
        delta_nv = std::min(delta_nv, delta_tmp);

        res[qv_idx] += delta_qv;
        res[qr_idx] -= delta_qv;
        res[Nr_idx] -= delta_nv;

#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " rain evaporation after Seifert (2008) dqr " << -delta_qv
                      << "\ndNr " << -delta_nv
                      << "\n";
#endif
#ifdef TRACE_QV
        if (trace)
            std::cout << "traj: " << cc.traj_id << " rain evaporation after Seifert (2008) dqv " << delta_qv
                << " g_d " << g_d << " Nr " << Nr << " mue " << mue
                << " lambda " << lambda << " f_v " << f_v
                << " s_sw " << s_sw
                << " gamma_eva " << gamma_eva
                << " x_r " << x_r
                << "\n";
#endif
        float_t delta_e = latent_heat_evap(T_prime) * delta_qv
            / specific_heat_dry_air(T_prime);
        // Evaporation, cooling
        if (delta_qv > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 * Sedimentation after Seifert and Beheng (2006). This function
 * processes rain, snow, ice, hail and graupel sedimentation at once.
 * The main difference in sedi_icon_core to ICON and COSMO here is that we do not
 * approximate the fall velocity based on the cell above the current one.
 * We simply don't have that information in a simulation along trajectories.
 *
 * Furthermore the improved scheme in sedi_icon_box_core originally
 * uses the information of each flux in all boxes above the current one.
 * This information isn't available as well.
 *
 * In short: Only sedimentation from the current box is calculated here.
 * Sedimentation that may come from further above is not used here.
 *
 * See http://www.cosmo-model.org/content/model/documentation/core/docu_sedi_twomom.pdf
 * for more details.
 *
 * @params T_prime Temperature [K]
 * @params S Saturation
 * @params qc_prime Cloud water mixing ratio
 * @params qr_prime Rain mixing ratio
 * @params Nr Number of rain droplets
 * @params qs_prime Snowflakes mixing ratio
 * @params Ns Number of snowflakes
 * @params qi_prime Ice crystal mixing ratio
 * @params Ni Number of ice crystals
 * @params qh_prime Hail mixing ratio
 * @params Nh Number of hail particles
 * @params qg_prime Graupel mixing ratio
 * @params Ng Number of graupel particles
 * @params p_prime Pressure [Pa]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void sedimentation_explicit(
    float_t &qc_prime,
    float_t &qr_prime,
    float_t &Nr,
    float_t &qs_prime,
    float_t &Ns,
    float_t &qi_prime,
    float_t &Ni,
    float_t &qh_prime,
    float_t &Nh,
    float_t &qg_prime,
    float_t &Ng,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {
    // Every other than from cloud should be the same
    float_t rhocorr = get_at(cc.rain.constants, Particle_cons_idx::rho_v);
    float_t v_n_sedi = 0.0;
    float_t v_q_sedi = 0.0;

    auto sedi_icon_core = [&](
        float_t &q,
        float_t &N,
        float_t &v_q_sedi,
        float_t &v_n_sedi,
        float_t &resQ,
        float_t &resN,
        float_t &resQOut,
        float_t &resNOut) {
        float_t v_nv = v_n_sedi;
        float_t v_qv = v_q_sedi;  // percentage how much trickles down
        // Assuming v_nv, v_qv is negative
        float_t z = 1/get_at(cc.constants, Cons_idx::inv_z);
        // Courant number; is it going to leave the box?
        float_t c_nv = -v_nv * get_at(cc.constants, Cons_idx::inv_z);
        float_t c_qv = -v_qv * get_at(cc.constants, Cons_idx::inv_z);
        // Can't loose more than there is in the box, hence the min().
        float_t flux_nv = N * std::min(c_nv, z);
        // This would be a good place to get the influx from cells above.
        // In case a grid based version will be implemented.
        // We could also make a guess of the following form
        // int32_t i = 1;
        // while (c_nv > 1) {
        //     c_nv--;
        //     i++;
        //     // we reduce N, assuming there is less above the current cell
        //     s_nv += N * pow(2, -i) * z * std::min(c_nv, 1);
        // }
        // But we already have added sedimentation from above in the main
        // function, so this might be an overestimation.
        float_t flux_qv = q * std::min(c_qv, z);
        // Same as above: Grid based implementations would have to add
        // calculations of layers above here.
        // A grid based implementation would need to get the net
        // of in and output here and store the flux for other
        // cells. We can skip this here and directly add, how much
        // is gone.
        flux_nv = abs(flux_nv);
        flux_qv = abs(flux_qv);
        // Avoid negative values
        float_t s_tmp = resN + N/dt;
        flux_nv = std::min(s_tmp, flux_nv);
        s_tmp = resQ + q/dt;
        flux_qv = std::min(s_tmp, flux_qv);

        // abs is used for paranoia reasons and should never be needed
        resN -= abs(flux_nv);
        resQ -= abs(flux_qv);
        resQOut -= abs(flux_qv);
        resNOut -= abs(flux_nv);
#ifdef TRACE_SEDI
        if (trace)
            std::cout << "traj: " << cc.traj_id << " \nWtihin sedi_icon_core\n"
                    << "v_qv: " << v_q_sedi
                    << "\nc_qv: " << -v_qv * get_at(cc.constants, Cons_idx::inv_z)
                    << "\ns_qv (<=1): " << v_qv*q
                    << "\ns_qv (>1):  " << q*get_at(cc.constants, Cons_idx::inv_z)
                    << "s_qv: " << s_qv
                    << "\nmax_s_qv: " << q/dt
                    << "\nc_qv: " << c_qv
                    << "\nv_qv: " << v_qv
                    << "\nv_q_sedi: " << v_q_sedi
                    << "\ns_nv: " << s_nv
                    << "\nc_nv: " << c_nv
                    << "\nv_nv: " << v_nv
                    << "\nv_n_sedi: " << v_n_sedi
                    << "\nq: " << q
                    << "\ninv_z: " << get_at(cc.constants, Cons_idx::inv_z)
                    << "\ns_nv: " << s_nv
                    << "\nresN: " << resN << "\n";
#endif
    };

    auto sedi_icon_sphere = [&](
        float_t &q,
        float_t &N,
        float_t &resQ,
        float_t &resN,
        float_t &resQOut,
        float_t &resNOut,
        particle_model_constants_t<float_t> &pc) {
        float_t v_n_sedi = 0.0;
        float_t v_q_sedi = 0.0;

        if (q > get_at(cc.constants, Cons_idx::q_crit)) {
            float_t x = particle_mean_mass(
                q, N, get_at(pc.constants, Particle_cons_idx::min_x_sedimentation),
                get_at(pc.constants, Particle_cons_idx::max_x));
            float_t lam = pow(get_at(pc.constants, Particle_cons_idx::lambda)*x,
                get_at(pc.constants, Particle_cons_idx::b_vel));
            float_t tmp = get_at(pc.constants, Particle_cons_idx::vsedi_min)/dt;
            float_t v_n = std::max(get_at(pc.constants, Particle_cons_idx::alfa_n) * lam, tmp);
            tmp = get_at(pc.constants, Particle_cons_idx::vsedi_min)/dt;
            float_t v_q = std::max(get_at(pc.constants, Particle_cons_idx::alfa_q) * lam, tmp);
            tmp = get_at(pc.constants, Particle_cons_idx::vsedi_max)/dt;
            v_n = std::min(v_n, tmp);
            tmp = get_at(pc.constants, Particle_cons_idx::vsedi_max)/dt;
            v_q = std::min(v_q, tmp);
            v_n *= rhocorr;
            v_q *= rhocorr;

            v_n_sedi -= v_n;
            v_q_sedi -= v_q;
#ifdef TRACE_SEDI
            if (trace)
                std::cout << "traj: " << cc.traj_id << " \nWithin sedi_icon_sphere"
                    << "\nalfa_q * lam: " << get_at(pc.constants, Particle_cons_idx::alfa_q) * lam
                    << "\nalfa_n * lam: " << get_at(pc.constants, Particle_cons_idx::alfa_n) * lam
                    << "\nv_q_sedi: " << v_q_sedi << "\nq: " << q
                    << "\nv_n_sedi: " << v_n_sedi << "\nN: " << N << "\nresN: " << resN
                    << "\nrhocorr: " << rhocorr
                    << "\nvsedi_max: " << get_at(pc.constants, Particle_cons_idx::vsedi_max)
                    << "\nalfa_q: " << get_at(pc.constants, Particle_cons_idx::alfa_q)
                    << "\nalfa_n: " << get_at(pc.constants, Particle_cons_idx::alfa_n)
                    << "\nvsedi_min: " << get_at(pc.constants, Particle_cons_idx::vsedi_min)
                    << "\nlambda: " << get_at(pc.constants, Particle_cons_idx::lambda)
                    << "\nb_vel: " << get_at(pc.constants, Particle_cons_idx::b_vel) << "\n";
#endif
        }
        float_t sedi_q = 0;
        float_t sedi_n = 0;
        sedi_icon_core(q, N, v_q_sedi, v_n_sedi, resQ, resN, resQOut, resNOut);
    };

    auto sedi_icon_sphere_lwf = [&](
       ) {
    };

    if (qr_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        float_t x_r = particle_mean_mass(
            qr_prime, Nr, get_at(cc.rain.constants, Particle_cons_idx::min_x_sedimentation),
            get_at(cc.rain.constants, Particle_cons_idx::max_x));
        float_t D_r = particle_diameter(
            x_r, get_at(cc.rain.constants, Particle_cons_idx::a_geo),
            get_at(cc.rain.constants, Particle_cons_idx::b_geo));
        float_t mue = (qc_prime >= get_at(cc.constants, Cons_idx::q_crit))
            ? (get_at(cc.rain.constants, Particle_cons_idx::nu)+1.0)
                / get_at(cc.rain.constants, Particle_cons_idx::b_geo) - 1.0
            : rain_mue_dm_relation(
                D_r, get_at(cc.rain.constants, Particle_cons_idx::cmu0),
                get_at(cc.rain.constants, Particle_cons_idx::cmu1),
                get_at(cc.rain.constants, Particle_cons_idx::cmu2),
                get_at(cc.rain.constants, Particle_cons_idx::cmu3),
                get_at(cc.rain.constants, Particle_cons_idx::cmu4));
        // inverse of lambda in Eq. (A08) SB (2008)
        float_t D_p = D_r * pow((mue+3)*(mue+2)*(mue+1), -1.0/3.0);

        // SB (2008), Eq. (A10)
        float_t v_nr =
            get_at(cc.rain.constants, Particle_cons_idx::alpha) - get_at(cc.rain.constants, Particle_cons_idx::beta)
            * pow(1.0+get_at(cc.rain.constants, Particle_cons_idx::gamma)*D_p, -mue-1.0);
        float_t v_qr =
            get_at(cc.rain.constants, Particle_cons_idx::alpha) - get_at(cc.rain.constants, Particle_cons_idx::beta)
            * pow(1.0+get_at(cc.rain.constants, Particle_cons_idx::gamma)*D_p, -mue-4.0);

        // Seifert (2008), Eq. A10
        // rhocorr: height dependency of particle fall speed
        v_nr *= rhocorr;
        v_qr *= rhocorr;
        v_n_sedi -= v_nr;
        v_q_sedi -= v_qr;
    }
#ifdef TRACE_QR
    auto beforeN = res[Nr_idx];
    auto beforeq = res[qr_idx];
#endif
    sedi_icon_core(qr_prime, Nr, v_q_sedi, v_n_sedi, res[qr_idx], res[Nr_idx], res[qr_out_idx], res[Nr_out_idx]);

#ifdef TRACE_QR
    if (trace)
        std::cout << "traj: " << cc.traj_id
            << " sedi_icon_core dqr " << res[qr_idx] - beforeq << ", dNr " << res[Nr_idx] - beforeN << "\n";
#endif

    // sedi_icon_sphere ice
    if (qi_prime > 0.0) {
#ifdef TRACE_QI
        auto before_q = res[qi_idx];
        auto before_n = res[Ni_idx];
#endif
        sedi_icon_sphere(qi_prime, Ni, res[qi_idx], res[Ni_idx], res[qi_out_idx], res[Ni_out_idx], cc.ice);
#ifdef TRACE_QI
        if (trace) {
            std::cout << "traj: " << cc.traj_id << " N before: " << before_n << " and after " << res[Ni_idx] << "\n";
            std::cout << "traj: " << cc.traj_id << " sedi icon sphere dqi " << res[qi_idx] - before_q
                      << ", dNi " << res[Ni_idx] - before_n << "\n";
        }
#endif
    }

    // sedi_icon_sphere snow
    if (qs_prime > 0.0) {
#ifdef TRACE_QS
        auto before_q = res[qs_idx];
        auto before_n = res[Ns_idx];
#endif
        sedi_icon_sphere(qs_prime, Ns, res[qs_idx], res[Ns_idx], res[qs_out_idx], res[Ns_out_idx], cc.snow);
#ifdef TRACE_QS
        if (trace)
            std::cout << "traj: " << cc.traj_id << " sedi icon sphere dqs " << res[qs_idx] - before_q
                      << ", dNs " << res[Ns_idx] - before_n << "\n";
#endif
    }
    // sedi_icon_sphere graupel
    const bool lprogmelt = false;  // true is not implemented
    if (qg_prime > 0.0) {
        if (lprogmelt) {
            sedi_icon_sphere_lwf();
        } else {
#ifdef TRACE_QG

            auto before_q = res[qg_idx];
            auto before_n = res[Ng_idx];
#endif
            sedi_icon_sphere(qg_prime, Ng, res[qg_idx], res[Ng_idx], res[qg_out_idx], res[Ng_out_idx], cc.graupel);
#ifdef TRACE_QG
            sediment_q = sediment_q + res[qg_idx].getValue() - before_q.getValue();
            sediment_n = sediment_n + res[Ng_idx].getValue() - before_n.getValue();
            if (trace)
                std::cout << "traj: " << cc.traj_id << " sedi icon sphere dqg " << res[qg_idx] - before_q
                      << ", dNg " << res[Ng_idx] - before_n << "\n";
#endif
        }
    }

    // sedi_icon_sphere hail
    if (qh_prime > 0.0) {
        if (lprogmelt) {
            sedi_icon_sphere_lwf();
        } else {
#ifdef TRACE_QH
            auto before_q = res[qh_idx];
            auto before_n = res[Nh_idx];
#endif
            sedi_icon_sphere(qh_prime, Nh, res[qh_idx], res[Nh_idx], res[qh_out_idx], res[Nh_out_idx], cc.hail);
#ifdef TRACE_QH
            if (trace)
                std::cout << "traj: " << cc.traj_id << " sedi icon sphere dqh " << res[qh_idx] - before_q
                          << ", dNh " << res[Nh_idx] - before_n << "\n";
#endif
        }
    }
}


/**
 * Evaporation from melting ice particles after Seifert & Beheng (2006), Section 2.5 and 3.3.
 * 10.1007/s00703-005-0112-4
 * The formulation is similar as for depositional growth of
 * a single hydrometeor (ice, snow, graupel) or evaporation
 * of rain droplets but with a surface temperature of 273.15 K.
 *
 * @params qv_prime Water vapor mixing ratio
 * @params e_d Partial pressure of water vapor
 * @params p_sat Saturation pressure of water [Pa]
 * @params s_sw Super saturation of water
 * @params T_prime Temperature [K]
 * @params q1 Mixing ratio of any ice particle
 * @params N1 Particle number of any ice particle
 * @params resq On out: difference of the ice particle mixing ratio
 * @params pc1 Model constants specific for the given ice particle
 * @params res Vector to store the changes
 * @params dt Time step size [s]
 */
template<class float_t>
void evaporation(
    float_t &p_sat,
    float_t &s_sw,
    float_t &T_prime,
    float_t &q1,
    float_t &N1,
    float_t &resq,
    particle_model_constants_t<float_t> &pc1,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (q1 > 0.0 && T_prime > get_at(cc.constants, Cons_idx::T_freeze)) {
        // Eq. 37 / Eq. 23
        float_t g_d = 4.0*M_PI / (get_at(cc.constants, Cons_idx::L_wd)*get_at(cc.constants, Cons_idx::L_wd)
            / (get_at(cc.constants, Cons_idx::K_T) * get_at(cc.constants, Cons_idx::R_v)
            * get_at(cc.constants, Cons_idx::T_freeze)*get_at(cc.constants, Cons_idx::T_freeze))
            + get_at(cc.constants, Cons_idx::R_v)*get_at(cc.constants, Cons_idx::T_freeze)
            / (get_at(cc.constants, Cons_idx::D_v) * p_sat));

        float_t x_1 = particle_mean_mass(
            q1, N1, get_at(pc1.constants, Particle_cons_idx::min_x_evap),
            get_at(pc1.constants, Particle_cons_idx::max_x));
        float_t d_1 = particle_diameter(
            x_1, get_at(pc1.constants, Particle_cons_idx::a_geo), get_at(pc1.constants, Particle_cons_idx::b_geo));
        // Eq. 28
        float_t v_1 = particle_velocity(
            x_1, get_at(pc1.constants, Particle_cons_idx::a_vel),
            get_at(pc1.constants, Particle_cons_idx::b_vel)) * get_at(pc1.constants, Particle_cons_idx::rho_v);
        // Eq. 40 / Eq. 30
        float_t f_v = get_at(pc1.constants, Particle_cons_idx::a_f)
            + get_at(pc1.constants, Particle_cons_idx::b_f) * sqrt(v_1*d_1);
        // Eq. 42
        float_t delta_q = g_d * N1 * get_at(pc1.constants, Particle_cons_idx::c_s) * d_1 * f_v * s_sw;
        float_t tmp = resq + q1/dt;
#ifdef TRACE_EVAP
        if (trace)
            std::cout << "Evap tmp " << tmp
                << "\ndelta_q " << delta_q
                << "\ng_d " << g_d
                << "\nN1 " << N1
                << "\nc_s " << get_at(pc1.constants, Particle_cons_idx::c_s)
                << "\nd_1 " << d_1
                << "\nf_v " << f_v
                << "\ns_sw " << s_sw
                << "\nv_1 " << v_1 << "\n";
#endif
        delta_q = std::min(tmp, std::max(-delta_q, float_t(0.0)));

        // Vapor
        res[qv_idx] += delta_q;
#ifdef TRACE_QV
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Evaporation from melting ice particles dqv " << delta_q << "\n";
#endif
        resq -= delta_q;

        float_t delta_e =  delta_q * latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze))
                / specific_heat_dry_air(T_prime);
        // Sublimination, cooling
        if (delta_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 * Depositional growth of all ice particles.
 * Deposition and sublimation is calculated, where deposition rates of ice and snow are
 * being stored.
 * Seifert (2006) Section 3.3
 * 10.1007/s00703-005-0112-4
 *  Depositional growth is based on
 * "A New Double-Moment Microphysics Parameterization for Application in Cloud and
 * Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov
 *
 * @params qv_prime Water vapor mixing ratio
 * @params qi_prime Ice crystal mixing ratio
 * @params Ni Number of ice crystals
 * @params qs_prime Snowflakes mixing ratio
 * @params Ns Number of snowflakes
 * @params qg_prime Graupel mixing ratio
 * @params Ng Number of graupel particles
 * @params qh_prime Hail mixing ratio
 * @params Nh Number of hail particles
 * @params s_si Over saturation of ice
 * @params p_sat_ice Saturation pressure of ice [Pa]
 * @params T_prime Temperature [K]
 * @params EPSILON Small value that the evaporation mass needs to exceed
 * @params dep_rate_ice On out: deposition rate of ice mass
 * @params dep_rate_snow On out: deposition rate of snow mass
 * @params D_vtp Diffusivity of dry air [m^2 s^-1]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void vapor_dep_relaxation(
    float_t &qv_prime,
    float_t &qi_prime,
    float_t &Ni,
    float_t &qs_prime,
    float_t &Ns,
    float_t &qg_prime,
    float_t &Ng,
    float_t &qh_prime,
    float_t &Nh,
    float_t &s_si,
    float_t &p_sat_ice,
    float_t &p_prime,
    float_t &T_prime,
    const double &EPSILON,
    float_t &dep_rate_ice,
    float_t &dep_rate_snow,
    float_t &D_vtp,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (T_prime < get_at(cc.constants, Cons_idx::T_freeze)) {
        float_t dep_ice = 0.0;
        float_t dep_snow = 0.0;
        float_t dep_graupel = 0.0;
        float_t dep_hail = 0.0;
        // Eq. 37
        float_t g_i = 4.0*M_PI / (get_at(cc.constants, Cons_idx::L_ed)*get_at(cc.constants, Cons_idx::L_ed)
            / (get_at(cc.constants, Cons_idx::K_T)*get_at(cc.constants, Cons_idx::R_v)*T_prime*T_prime)
            + get_at(cc.constants, Cons_idx::R_v)*T_prime / (D_vtp*p_sat_ice));

        auto vapor_deposition = [&](
            float_t &q,
            float_t &N,
            particle_model_constants_t<float_t> &pc,
            float_t &dep) {
            if (q <= 0.0) {
                dep = 0.0;
            } else {
                float_t x = particle_mean_mass(
                    q, N, get_at(pc.constants, Particle_cons_idx::min_x_depo),
                    get_at(pc.constants, Particle_cons_idx::max_x));
                float_t d = particle_diameter(
                    x, get_at(pc.constants, Particle_cons_idx::a_geo),
                    get_at(pc.constants, Particle_cons_idx::b_geo));
                float_t v = particle_velocity(
                    x, get_at(pc.constants, Particle_cons_idx::a_vel),
                    get_at(pc.constants, Particle_cons_idx::b_vel)) * get_at(pc.constants, Particle_cons_idx::rho_v);
                float_t f_v = get_at(pc.constants, Particle_cons_idx::a_f)
                    + get_at(pc.constants, Particle_cons_idx::b_f) * sqrt(d*v);
                // Eq. 40
                f_v = std::max(f_v,
                    get_at(pc.constants, Particle_cons_idx::a_f)/get_at(pc.constants, Particle_cons_idx::a_ven));
                // Eq. 42
                dep = g_i * N * get_at(pc.constants, Particle_cons_idx::c_s) * d * f_v * s_si;
            }
        };

        vapor_deposition(qi_prime, Ni, cc.ice, dep_ice);
        vapor_deposition(qs_prime, Ns, cc.snow, dep_snow);
        vapor_deposition(qg_prime, Ng, cc.graupel, dep_graupel);
        vapor_deposition(qh_prime, Nh, cc.hail, dep_hail);
        // Saturation after qv has been put towards ice such that
        // saturation of ice is 1.
        float_t qvsidiff = qv_prime - convert_Si_to_qv(p_prime, T_prime, float_t(1),
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_a),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_b));
        if (abs(qvsidiff) > EPSILON) {
            float_t tau_i_i = 1.0/qvsidiff*dep_ice/dt;
            float_t tau_s_i = 1.0/qvsidiff*dep_snow/dt;
            float_t tau_g_i = 1.0/qvsidiff*dep_graupel/dt;
            float_t tau_h_i = 1.0/qvsidiff*dep_hail/dt;

            float_t xi_i = tau_i_i + tau_s_i + tau_g_i + tau_h_i;

            float_t xfac = (xi_i < EPSILON) ?
                (float_t) 0.0 : qvsidiff/xi_i * (1.0-exp(-xi_i));

            dep_ice     = xfac * tau_i_i;
            dep_snow    = xfac * tau_s_i;
            dep_graupel = xfac * tau_g_i;
            dep_hail    = xfac * tau_h_i;

            if (qvsidiff < 0.0) {
                float_t dep_max_tmp = -res[qi_idx]-qi_prime/dt;
                dep_ice     = std::max(dep_ice, dep_max_tmp);
                dep_max_tmp = -res[qs_idx]-qs_prime/dt;
                dep_snow    = std::max(dep_snow, dep_max_tmp);
                dep_max_tmp = -res[qg_idx]-qg_prime/dt;
                dep_graupel = std::max(dep_graupel, dep_max_tmp);
                dep_max_tmp = -res[qh_idx]-qh_prime/dt;
                dep_hail    = std::max(dep_hail, dep_max_tmp);
            } else {
                float_t tmp_sum = dep_ice + dep_graupel + dep_snow + dep_hail;
                if (tmp_sum > res[qv_idx] + qv_prime/dt) {
                    dep_ice = dep_ice/tmp_sum * (res[qv_idx] + qv_prime/dt);
                    dep_snow = dep_snow/tmp_sum * (res[qv_idx] + qv_prime/dt);
                    dep_graupel = dep_graupel/tmp_sum * (res[qv_idx] + qv_prime/dt);
                    dep_hail = dep_hail/tmp_sum * (res[qv_idx] + qv_prime/dt);
                }
            }

            float_t dep_sum = dep_ice + dep_graupel + dep_snow + dep_hail;

            res[qi_idx] += dep_ice;
            res[qs_idx] += dep_snow;
            res[qg_idx] += dep_graupel;
            res[qh_idx] += dep_hail;
            res[qv_idx] -= dep_sum;
#ifdef TRACE_QI
            if (trace) {
                std::cout << "traj: " << cc.traj_id
                    << " Depos qvsi " << p_sat_ice /(get_at(cc.constants, Cons_idx::R_v)*T_prime) << "\n";
                std::cout << "traj: " << cc.traj_id << " Depos qv " << qv_prime << "\n";
                std::cout << "traj: " << cc.traj_id << " Depos growth dqi " << dep_ice << "\n";
            }
#endif
#ifdef TRACE_QS
            if (trace)
                std::cout << "traj: " << cc.traj_id << " Depos growth dqs " << dep_snow << "\n";
#endif
#ifdef TRACE_QG
            if (trace) {
                std::cout << "traj: " << cc.traj_id << " Depos growth dqg " << dep_graupel << "\n";
            }
#endif
#ifdef TRACE_QH
            if (trace)
                std::cout << "traj: " << cc.traj_id << " Depos growth dqh " << dep_hail << "\n";
#endif
#ifdef TRACE_QV
            if (trace)
                std::cout << "traj: " << cc.traj_id << " Depos growth dqv " << -dep_sum << "\n";
#endif

            dep_rate_ice += dep_ice;
            dep_rate_snow += dep_snow;

            if (dep_ice > 0) {
                res[depo_idx] += dep_ice;
            } else {
                res[sub_idx] += dep_ice;
            }
            if (dep_snow > 0) {
                res[depo_idx] += dep_snow;
            } else {
                res[sub_idx] += dep_snow;
            }
            if (dep_graupel > 0) {
                res[depo_idx] += dep_graupel;
            } else {
                res[sub_idx] += dep_graupel;
            }
            if (dep_hail > 0) {
                res[depo_idx] += dep_hail;
            } else {
                res[sub_idx] += dep_hail;
            }
            float_t delta_e = latent_heat_melt(
                T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * dep_sum
                / specific_heat_dry_air(T_prime);
            // Sublimation, cooling
            if (dep_sum > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] -= delta_e;
        }
    }
}


/**
 * Particle collection for ice particles such as:
 * graupel and ice  to graupel
 * graupel and snow to graupel
 * hail and ice     to hail
 * hail and snow    to hail
 * snow and ice     to snow
 * This function processes only one of those at a time.
 * After Seifert & Beheng (2006), Section 3.5.1.
 * 10.1007/s00703-005-0112-4
 * The collision integrals are replaced by a mean
 * efficiency that depends on the mean masses,
 * a characteristic velocity difference and the
 * diameter of the involved particles.
 * The mean velocity is approximated using a variance ansatz.
 *
 * @params q1 Mixing mass ratio of collecting particle
 * @params q2 Mixing mass ratio of particle being collected
 * @params N1 Number of particles of collecting particle
 * @params N2 Number of particles of particle being collected
 * @params T_c Difference between temperature and freezing temperature [K]
 * @params coeffs Model constants for collision processes
 * @params pc1 Model constants for a particle of collecting particle
 * @params pc2 Model constants for a particle of particle being collected
 * @params dt Time step size [s]
 *
 * @return Vector of two float_t with the change in particle number (0)
 * and change in mass mixing ratio (1)
 */
template<class float_t>
std::vector<float_t> particle_collection(
    float_t &q1,
    float_t &q2,
    float_t &N1,
    float_t &N2,
    float_t &T_c,
    collection_model_constants_t<float_t> &coeffs,
    particle_model_constants_t<float_t> &pc1,
    particle_model_constants_t<float_t> &pc2,
    const double &dt) {

    float_t e_coll = std::min(exp(0.09*T_c), 1.0);
    float_t x_1 = particle_mean_mass(
        q1, N1, get_at(pc1.constants, Particle_cons_idx::min_x_collection),
        get_at(pc1.constants, Particle_cons_idx::max_x));
    float_t d_1 = particle_diameter(
        x_1, get_at(pc1.constants, Particle_cons_idx::a_geo),
        get_at(pc1.constants, Particle_cons_idx::b_geo));
    float_t v_1 = particle_velocity(
        x_1, get_at(pc1.constants, Particle_cons_idx::a_vel),
        get_at(pc1.constants, Particle_cons_idx::b_vel)) * get_at(pc1.constants, Particle_cons_idx::rho_v);

    float_t x_2 = particle_mean_mass(
        q2, N2, get_at(pc2.constants, Particle_cons_idx::min_x_collection),
        get_at(pc2.constants, Particle_cons_idx::max_x));
    float_t d_2 = particle_diameter(
        x_2, get_at(pc2.constants, Particle_cons_idx::a_geo),
        get_at(pc2.constants, Particle_cons_idx::b_geo));
    float_t v_2 = particle_velocity(
        x_2, get_at(pc2.constants, Particle_cons_idx::a_vel),
        get_at(pc2.constants, Particle_cons_idx::b_vel)) * get_at(pc2.constants, Particle_cons_idx::rho_v);

    float_t coll_n = M_PI/4 * N1 * N2 * e_coll
        * (coeffs.delta_n_aa * d_2 * d_2
            + coeffs.delta_n_ab * d_2 * d_1
            + coeffs.delta_n_bb * d_1 * d_1)
        * sqrt(coeffs.theta_n_aa * v_2 * v_2
            - coeffs.theta_n_ab * v_2 * v_1
            + coeffs.theta_n_bb * v_1 * v_1
            + get_at(pc1.constants, Particle_cons_idx::s_vel) * get_at(pc1.constants, Particle_cons_idx::s_vel));
    float_t coll_q = M_PI/4 * q1 * N2 * e_coll
        * (coeffs.delta_q_aa * d_2 * d_2
            + coeffs.delta_q_ab * d_2 * d_1
            + coeffs.delta_q_bb * d_1 * d_1)
        * sqrt(coeffs.theta_q_aa * v_2 * v_2
            - coeffs.theta_q_ab * v_2 * v_1
            + coeffs.theta_q_bb * v_1 * v_1
            + get_at(pc1.constants, Particle_cons_idx::s_vel) * get_at(pc1.constants, Particle_cons_idx::s_vel));

    coll_n = std::min(N1/dt, coll_n);
    coll_q = std::min(q1/dt, coll_q);

    std::vector<float_t> r(2);
    r[0] = coll_n;
    r[1] = coll_q;
    return r;
}


/**
 * Ice particle collections for the following processes:
 * graupel and ice      to graupel
 * graupel and snow     to graupel
 * snow and ice         to snow
 * graupel and graupel  to graupel
 * All those processes are done in this function at once.
 * Technically, this could be used for hail instead of graupel as well.
 * After Seifert & Beheng (2006), Sections 3.5.1 and 3.5.2.
 * 10.1007/s00703-005-0112-4
 * The collision integrals are replaced by a mean
 * efficiency that depends on the mean masses,
 * a characteristic velocity difference and the
 * diameter of the involved particles.
 * The mean velocity is approximated using a variance ansatz.
 *
 * @params qi_prime Ice crystal mixing ratio
 * @params Ni Number of ice crystals
 * @params qs_prime Snowflakes mixing ratio
 * @params Ns Number of snowflakes
 * @params qg_prime Graupel mixing ratio
 * @params Ng Number of graupel particles
 * @params T_prime Temperature [K]
 * @params T_c Difference between temperature and freezing temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void particle_particle_collection(
    float_t &qi_prime,
    float_t &Ni,
    float_t &qs_prime,
    float_t &Ns,
    float_t &qg_prime,
    float_t &Ng,
    float_t &T_prime,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    //// particle_collection snow
    if (qi_prime > get_at(cc.constants, Cons_idx::q_crit) && qs_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        std::vector<float_t> delta = particle_collection(
            qi_prime, qs_prime, Ni, Ns, T_c, cc.coeffs_sic, cc.ice, cc.snow, dt);
        res[qs_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
#ifdef TRACE_QI
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " ice-snow collision dqi " << -delta[1] << ", dNi " << -delta[0] << "\n";
#endif
#ifdef TRACE_QS
        if (trace)
            std::cout << "traj: " << cc.traj_id << " ice-snow collision dqs " << delta[1] << "\n";
#endif
    }

    //// graupel self collection
    if (qg_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        float_t x_g = particle_mean_mass(qg_prime, Ng,
            get_at(cc.graupel.constants, Particle_cons_idx::min_x_collection),
            get_at(cc.graupel.constants, Particle_cons_idx::max_x));
        float_t d_g = particle_diameter(x_g,
            get_at(cc.graupel.constants, Particle_cons_idx::a_geo),
            get_at(cc.graupel.constants, Particle_cons_idx::b_geo));
        float_t v_g = particle_velocity(x_g,
            get_at(cc.graupel.constants, Particle_cons_idx::a_vel),
            get_at(cc.graupel.constants, Particle_cons_idx::b_vel))
            * get_at(cc.graupel.constants, Particle_cons_idx::rho_v);
        float_t delta_n = get_at(cc.graupel.constants, Particle_cons_idx::sc_coll_n)
            * Ng * Ng * d_g * d_g * v_g;
        // sticking efficiency does only distinguish dry and wet
        delta_n *= (T_prime > get_at(cc.constants, Cons_idx::T_freeze))
            ? get_at(cc.constants, Cons_idx::ecoll_gg_wet) : get_at(cc.constants, Cons_idx::ecoll_gg);
        float_t delta_n_min = res[Ng_idx] + Ng/dt;
        delta_n = std::min(delta_n, delta_n_min);

        res[Ng_idx] -= delta_n;
#ifdef TRACE_QG
        if (trace)
            std::cout << "traj: " << cc.traj_id << " graupel self collection dNg " << -delta_n << "\n";
#endif
    }

    // particle particle collection
    // ice and graupel collision
    if (qi_prime > get_at(cc.constants, Cons_idx::q_crit) && qg_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        std::vector<float_t> delta = particle_collection(
            qi_prime, qg_prime, Ni, Ng, T_c, cc.coeffs_gic, cc.ice, cc.graupel, dt);
        res[qg_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
#ifdef TRACE_QG
        if (trace)
            std::cout << "traj: " << cc.traj_id << " ice-graupel collision dqg " << delta[1] << "\n";
#endif
#ifdef TRACE_QI
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " ice-graupel collision dqi " << -delta[1] << ", dNi " << -delta[0] << "\n";
#endif
    }

    // particle particle collection
    // snow and graupel collision
    if (qs_prime > get_at(cc.constants, Cons_idx::q_crit) && qg_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        std::vector<float_t> delta = particle_collection(
            qs_prime, qg_prime, Ns, Ng, T_c, cc.coeffs_gsc, cc.snow, cc.graupel, dt);
        res[qg_idx] += delta[1];
        res[qs_idx] -= delta[1];
        res[Ns_idx] -= delta[0];
#ifdef TRACE_QG
        if (trace)
            std::cout << "traj: " << cc.traj_id << " snow-graupel collision dqg " << delta[1] << "\n";
#endif
#ifdef TRACE_QS
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " snow-graupel collision dqs " << -delta[1] << ", dNs " << -delta[0] << "\n";
#endif
    }
}


/**
 * Conversion graupel to hail based on look-up tables
 * where supercooled liquid water is present.
 * Ansatz by Ulrich Blahak,
 * "Towards a better representation of high density
 * ice particles in a state-of-the-art two-moment
 * bulk microphysical scheme."
 * Proc. 15th Int. Conf. Clouds and Precip.,
 * Cancun, Mexico. Vol. 20208. 2008
 *
 *
 * @params qc_prime Cloud water mixing ratio
 * @params qr_prime Rain mixing ratio
 * @params qi_prime Ice crystal mixing ratio
 * @params qg_prime Graupel mixing ratio
 * @params Ng Number of graupel particles
 * @params qh_prime Hail mixing ratio
 * @params Nh Number of hail particles
 * @params p_prime Pressure [Pa]
 * @params T_prime Temperature [K]
 * @params T_c Difference between temperature and freezing temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void graupel_hail_conv(
    float_t &qc_prime,
    float_t &qr_prime,
    float_t &qi_prime,
    float_t &qg_prime,
    float_t &Ng,
    float_t &p_prime,
    float_t &T_prime,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    float_t x_g = particle_mean_mass(qg_prime, Ng,
            get_at(cc.graupel.constants, Particle_cons_idx::min_x_conversion),
            get_at(cc.graupel.constants, Particle_cons_idx::max_x));
    float_t d_g = particle_diameter(x_g,
        get_at(cc.graupel.constants, Particle_cons_idx::a_geo),
        get_at(cc.graupel.constants, Particle_cons_idx::b_geo));
    float_t Ng_tmp = qg_prime/x_g;
    // supercooled liquid water = rain + cloud water
    float_t qwa_prime = qr_prime + qc_prime;

    if (qwa_prime > 1e-3 && T_c < 0 && qg_prime > get_at(cc.graupel.constants, Particle_cons_idx::q_crit_c)) {
        float_t d_sep = wet_growth_diam(p_prime, T_prime, qwa_prime,
            qi_prime, cc.ltabdminwgg);
        if (d_sep > 0.0 && d_sep < 10.0*d_g) {
            float_t xmin = pow(d_sep/get_at(cc.graupel.constants, Particle_cons_idx::a_geo),
                1.0/get_at(cc.graupel.constants, Particle_cons_idx::b_geo));
            float_t lam = pow(get_at(cc.graupel.constants, Particle_cons_idx::g2)
                / (get_at(cc.graupel.constants, Particle_cons_idx::g1)*x_g),
                get_at(cc.graupel.constants, Particle_cons_idx::mu));
            xmin = pow(xmin, get_at(cc.graupel.constants, Particle_cons_idx::mu));
            float_t n_0 = get_at(cc.graupel.constants, Particle_cons_idx::mu) * Ng
                * pow(lam, get_at(cc.graupel.constants, Particle_cons_idx::nm1))
                / get_at(cc.graupel.constants, Particle_cons_idx::g1);
            float_t lam_xmin = lam*xmin;

            float_t conv_n = n_0 / (get_at(cc.graupel.constants, Particle_cons_idx::mu)
                * pow(lam, get_at(cc.graupel.constants, Particle_cons_idx::nm1)))
                * cc.table_g1.look_up(lam_xmin);
            float_t conv_q = n_0 / (get_at(cc.graupel.constants, Particle_cons_idx::mu)
                * pow(lam, get_at(cc.graupel.constants, Particle_cons_idx::nm2)))
                * cc.table_g2.look_up(lam_xmin);
            float_t conv_tmp = res[Ng_idx] + Ng/dt;
            conv_n = std::min(conv_n, conv_tmp);
            conv_tmp = res[qg_idx] + qg_prime/dt;
            conv_q = std::min(conv_q, conv_tmp);

            // Graupel
            res[qg_idx] -= conv_q;
            res[Ng_idx] -= conv_n;
            // Hail
            res[qh_idx] += conv_q;
            res[Nh_idx] += conv_n;

#ifdef TRACE_QG
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " conversion graupel->hail dqg " << - conv_q << ", dNq " << - conv_n << "\n";
#endif
#ifdef TRACE_QH
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " conversion graupel->hail dqh " <<  conv_q << ", dNh " << conv_n << "\n";
#endif
        }
    }
}


/**
 * Hail collision with ice and snow after Seifert & Beheng (2006), Section 3.5.1.
 * 10.1007/s00703-005-0112-4
 * Uses particle_collection().
 *
 * @params qh_prime Hail mixing ratio
 * @params Nh Number of hail particle
 * @params qs_prime Snowflakes mixing ratio
 * @params Ns Number of snowflakes
 * @params qi_prime Ice crystal mixing ratio
 * @params Ni Number of ice crystals
 * @params T_c Difference between temperature and freezing temperature [K]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void hail_collision(
    float_t &qh_prime,
    float_t &Nh,
    float_t &qs_prime,
    float_t &Ns,
    float_t &qi_prime,
    float_t &Ni,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    // ice and hail
    if (qi_prime > get_at(cc.constants, Cons_idx::q_crit) && qh_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        std::vector<float_t> delta = particle_collection(
            qi_prime, qh_prime, Ni, Nh, T_c, cc.coeffs_hic, cc.ice, cc.hail, dt);
        res[qh_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
#ifdef TRACE_QH
        if (trace)
            std::cout << "traj: " << cc.traj_id << " ice-hail collision dqh " << delta[1] << "\n";
#endif
#ifdef TRACE_QI
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " ice-hail collision dqi " << -delta[1] << ", dNi " << -delta[0] << "\n";
#endif
    }

    // snow and hail collision
    if (qs_prime > get_at(cc.constants, Cons_idx::q_crit) && qh_prime > get_at(cc.constants, Cons_idx::q_crit)) {
        std::vector<float_t> delta = particle_collection(
            qs_prime, qh_prime, Ns, Nh, T_c, cc.coeffs_hsc, cc.snow, cc.hail, dt);
        res[qh_idx] += delta[1];
        res[qs_idx] -= delta[1];
        res[Ns_idx] -= delta[0];
#ifdef TRACE_QH
        if (trace)
            std::cout << "traj: " << cc.traj_id << " snow-hail collision dqh " << delta[1] << "\n";
#endif
#ifdef TRACE_QS
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " snow-hail collision dqs " << -delta[1] << ", dNs " << -delta[0] << "\n";
#endif
    }
}


/**
 * Rate of ice or snow collecting cloud droplets after Seifert & Beheng (2006), Eqs. 61-67.
 * 10.1007/s00703-005-0112-4
 * Zero velocity variance is assumed for cloud droplets.
 *
 * @params qc_prime Cloud water mixing ratio
 * @params Nc Number of cloud droplets
 * @params q1 Mixing mass ratio
 * @params N1 Number of particles
 * @params pc1 Model constants specific for the given ice particle
 * @params coeffs Model constants for collision processes
 * @params rime_rate_qb On out: Riming rate for mixing mass
 * @params rime_rate_nb On out: Riming rate for number of particles
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void riming_cloud_core(
    float_t &qc_prime,
    float_t &Nc,
    float_t &q1,
    float_t &N1,
    particle_model_constants_t<float_t> &pc1,
    collection_model_constants_t<float_t> &coeffs,
    float_t &rime_rate_qb,
    float_t &rime_rate_nb,
    model_constants_t<float_t> &cc,
    const double &dt) {

    float_t x_1 = particle_mean_mass(
        q1, N1, get_at(pc1.constants, Particle_cons_idx::min_x_riming),
        get_at(pc1.constants, Particle_cons_idx::max_x));
    float_t d_1 = particle_diameter(
        x_1, get_at(pc1.constants, Particle_cons_idx::a_geo),
        get_at(pc1.constants, Particle_cons_idx::b_geo));
    float_t x_c = particle_mean_mass(
        qc_prime, Nc, get_at(cc.cloud.constants, Particle_cons_idx::min_x_riming),
        get_at(cc.cloud.constants, Particle_cons_idx::max_x));
    float_t d_c = particle_diameter(
        x_c, get_at(cc.cloud.constants, Particle_cons_idx::a_geo),
        get_at(cc.cloud.constants, Particle_cons_idx::b_geo));

    float_t const1 = get_at(cc.constants, Cons_idx::const0) * get_at(pc1.constants, Particle_cons_idx::ecoll_c);

    if (qc_prime > get_at(cc.cloud.constants, Particle_cons_idx::q_crit_c)
        && q1 > get_at(pc1.constants, Particle_cons_idx::q_crit_c)
        && d_c > get_at(cc.cloud.constants, Particle_cons_idx::d_crit_c)
        && d_1 > get_at(pc1.constants, Particle_cons_idx::d_crit_c)) {
        float_t v_1 = particle_velocity(x_1, get_at(pc1.constants, Particle_cons_idx::a_vel),
            get_at(pc1.constants, Particle_cons_idx::b_vel)) * get_at(pc1.constants, Particle_cons_idx::rho_v);
        float_t v_c = particle_velocity(x_c, get_at(cc.cloud.constants, Particle_cons_idx::a_vel),
            get_at(cc.cloud.constants, Particle_cons_idx::b_vel))
            * get_at(cc.cloud.constants, Particle_cons_idx::rho_v);
        float_t tmp = const1*(d_c - get_at(cc.cloud.constants, Particle_cons_idx::d_crit_c));
        float_t coll_tmp1 = get_at(pc1.constants, Particle_cons_idx::ecoll_c)/dt;
        float_t coll_tmp2 = get_at(cc.constants, Cons_idx::ecoll_min)/dt;
        float_t e_coll = std::min(coll_tmp1, std::max(tmp, coll_tmp2));

        rime_rate_qb = M_PI/4.0 * e_coll * N1 * qc_prime
            * (coeffs.delta_q_aa * d_1*d_1
                + coeffs.delta_q_ab * d_1*d_c
                + coeffs.delta_q_bb * d_c*d_c)
            * sqrt(coeffs.theta_q_aa * v_1*v_1
                - coeffs.theta_q_ab * v_1*v_c
                + coeffs.theta_q_bb * v_c*v_c
                + get_at(pc1.constants, Particle_cons_idx::s_vel)*get_at(pc1.constants, Particle_cons_idx::s_vel));

        rime_rate_nb = M_PI/4.0 * e_coll * N1 * Nc
            * (coeffs.delta_n_aa * d_1*d_1
                + coeffs.delta_n_ab * d_1*d_c
                + coeffs.delta_n_bb * d_c*d_c)
            * sqrt(coeffs.theta_n_aa * v_1*v_1
                - coeffs.theta_n_ab * v_1*v_c
                + coeffs.theta_n_bb * v_c*v_c
                + get_at(pc1.constants, Particle_cons_idx::s_vel)*get_at(pc1.constants, Particle_cons_idx::s_vel));
    } else {
        rime_rate_qb = 0.0;
        rime_rate_nb = 0.0;
    }
}


/**
 * Riming of rain droplets with ice or snow particles after Seifert & Beheng (2006), Eqs 61-63.
 * 10.1007/s00703-005-0112-4
 * Zero velocity variance is assumed for rain droplets.
 *
 * @params qr_prime Rain mixing ratio
 * @params Nr Number of rain droplets
 * @params q1 Mixing ratio of any ice particle
 * @params N1 Particle number of any ice particle
 * @params pc1 Model constants specific for the given ice particle
 * @params coeffs Model constants for collision processes
 * @params rime_rate_qa On out: Riming rate for mixing mass
 * @params rime_rate_qb On out: Riming rate for rain mixing mass
 * @params rime_rate_nb On out: Riming rate for number of rain droplets
 */
template<class float_t>
void riming_rain_core(
    float_t &qr_prime,
    float_t &Nr,
    float_t &q1,
    float_t &N1,
    particle_model_constants_t<float_t> &pc1,
    collection_model_constants_t<float_t> &coeffs,
    float_t &rime_rate_qa,
    float_t &rime_rate_qb,
    float_t &rime_rate_nb,
    model_constants_t<float_t> &cc) {

    float_t x_1 = particle_mean_mass(q1, N1, get_at(pc1.constants, Particle_cons_idx::min_x_riming),
        get_at(pc1.constants, Particle_cons_idx::max_x));
    float_t d_1 = particle_diameter(x_1, get_at(pc1.constants, Particle_cons_idx::a_geo),
        get_at(pc1.constants, Particle_cons_idx::b_geo));

    if (qr_prime > get_at(cc.constants, Cons_idx::q_crit)
        && q1 > get_at(cc.constants, Cons_idx::q_crit_r) && d_1 > get_at(cc.constants, Cons_idx::D_crit_r)) {
        float_t x_r = particle_mean_mass(qr_prime, Nr, get_at(cc.rain.constants, Particle_cons_idx::min_x_riming),
            get_at(cc.rain.constants, Particle_cons_idx::max_x));
        float_t d_r = particle_diameter(x_r, get_at(cc.rain.constants, Particle_cons_idx::a_geo),
            get_at(cc.rain.constants, Particle_cons_idx::b_geo));

        float_t v_1 = particle_velocity(x_1, get_at(pc1.constants, Particle_cons_idx::a_vel),
            get_at(pc1.constants, Particle_cons_idx::b_vel)) * get_at(pc1.constants, Particle_cons_idx::rho_v);
        float_t v_r = particle_velocity(x_r, get_at(cc.rain.constants, Particle_cons_idx::a_vel),
            get_at(cc.rain.constants, Particle_cons_idx::b_vel)) * get_at(cc.rain.constants, Particle_cons_idx::rho_v);

        rime_rate_qb = M_PI/4.0 * N1 * qr_prime
            * (coeffs.delta_n_aa * d_1*d_1
                + coeffs.delta_q_ab * d_1*d_r
                + coeffs.delta_q_bb * d_r*d_r)
            * sqrt(coeffs.theta_n_aa * v_1*v_1
                - coeffs.theta_q_ab * v_1*v_r
                + coeffs.theta_q_bb * v_r*v_r
                + get_at(pc1.constants, Particle_cons_idx::s_vel)*get_at(pc1.constants, Particle_cons_idx::s_vel));

        rime_rate_qa = M_PI/4.0 * Nr * q1
            * (coeffs.delta_q_aa * d_1*d_1
                + coeffs.delta_q_ba * d_1*d_r
                + coeffs.delta_n_bb * d_r*d_r)
            * sqrt(coeffs.theta_q_aa * v_1*v_1
                - coeffs.theta_q_ba * v_1*v_r
                + coeffs.theta_n_bb * v_r*v_r
                + get_at(pc1.constants, Particle_cons_idx::s_vel)*get_at(pc1.constants, Particle_cons_idx::s_vel));

        rime_rate_nb = M_PI/4.0 * N1 * Nr
            * (coeffs.delta_n_aa * d_1*d_1
                + coeffs.delta_n_ab * d_1*d_r
                + coeffs.delta_n_bb * d_r*d_r)
            * sqrt(coeffs.theta_n_aa * v_1*v_1
                - coeffs.theta_n_ab * v_1*v_r
                + coeffs.theta_n_bb * v_r*v_r
                + get_at(pc1.constants, Particle_cons_idx::s_vel)*get_at(pc1.constants, Particle_cons_idx::s_vel));
    } else {
        rime_rate_qa = 0.0;
        rime_rate_qb = 0.0;
        rime_rate_nb = 0.0;
    }
}


/**
 * Riming of cloud and rain droplets on ice.
 * 10.1007/s00703-005-0112-4
 * Uses rates calculated via riming_rain_core() and riming_cloud_core() (each based on after Seifert & Beheng; 2006)
 * and breaks up ice if ice gets too big (Hallet-Mossop process) and converts small ice particles to graupel.
 * Multiplication of ice is based on Beheng, K. D. "Numerical study on the combined action of droplet coagulation,
 * ice particle riming and the splintering process concerning maritime cumuli."
 * Contrib. Atmos. Phys.;(Germany, Federal Republic of) 55.3 (1982)., Eq. 7, conversion to graupel Eqs. 70-71.
 *
 *
 * @params qc_prime Cloud water mixing ratio
 * @params Nc Number of cloud droplets
 * @params qr_prime Rain mixing ratio
 * @params Nr Number of rain droplets
 * @params qi_prime Ice crystal mixing ratio
 * @params Ni Number of ice crystals
 * @params dep_rate_ice Deposition rate of ice mass
 * @params rime_rate_qc Riming rate for cloud mixing mass
 * @params rime_rate_nc Riming rate for number of cloud droplets
 * @params rime_rate_qr Riming rate for rain mixing mass
 * @params rime_rate_nr Riming rate for number of rain droplets
 * @params rime_rate_qi Riming rate for ice mixing mass
 * @params T_prime Current temperature [K]
 * @params dt Timestep size [s]
 * @params res Vector to store the changes
 * @params cc Model constants
 */
template<class float_t>
void ice_riming(
    float_t &qc_prime,
    float_t &Nc,
    float_t &qr_prime,
    float_t &Nr,
    float_t &qi_prime,
    float_t &Ni,
    float_t &dep_rate_ice,
    float_t &rime_rate_qc,
    float_t &rime_rate_nc,
    float_t &rime_rate_qr,
    float_t &rime_rate_nr,
    float_t &rime_rate_qi,
    float_t &T_prime,
    const double &dt,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc) {

    if (dep_rate_ice > 0.0 && dep_rate_ice >= rime_rate_qc+rime_rate_qr) {
        // Depositional growth is stronger than riming growth, therefore ice stays ice
        // ice cloud riming
        if (rime_rate_qc > 0.0) {
            float_t rime_tmp = res[qc_idx] + qc_prime/dt;
            float_t rime_q = std::max(float_t(0), std::min(rime_tmp, rime_rate_qc));
            rime_tmp = res[Nc_idx] + Nc/dt;
            float_t rime_n = std::max(float_t(0), std::min(rime_tmp, rime_rate_nc));
            // Ice
            res[qi_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;

            // Cloud N
            res[Nc_idx] -= rime_n;

#ifdef TRACE_QC
            if (trace)
                std::cout << "traj: " << cc.traj_id << " ice riming dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
#ifdef TRACE_QI
            if (trace)
                std::cout << "traj: " << cc.traj_id << " ice riming dqi " << rime_q << "\n";
#endif
            if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
                float_t mult_1 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
                float_t mult_2 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
                mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
                mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
                // Ice N
                res[Ni_idx] += get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q;
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " ice riming with mult dNi "
                              <<  get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q << "\n";
#endif
            }
        }
        // ice rain riming
        if (rime_rate_qr > 0.0) {
            float_t rime_tmp = res[qr_idx] + qr_prime/dt;
            float_t rime_q = std::max(float_t(0), std::min(rime_rate_qr, rime_tmp));
            rime_tmp = res[Nr_idx] + Nr/dt;
            float_t rime_n = std::max(float_t(0), std::min(rime_tmp, rime_rate_nr));

            // Ice
            res[qi_idx] += rime_q;
            // Rain
            res[qr_idx] -= rime_q;
            // Rain N
            res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " ice rain riming dqr " << -rime_q << ", dNr " << -rime_n << "\n";
#endif
#ifdef TRACE_QI
            if (trace)
                std::cout << "traj: " << cc.traj_id << " ice rain riming A dqi " << rime_q << "\n";
#endif
            if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
                float_t mult_1 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
                float_t mult_2 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
                mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
                mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
                // Ice N
                res[Ni_idx] += get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q;
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " ice rain rimingwith mult dNi "
                              << get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q << "\n";
#endif
            }

            float_t delta_e = latent_heat_melt(
                T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_q
                / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;
        }
    } else {
        // Depositional growth negative or smaller than riming growth,
        // therefore ice is allowed to convert to graupel and / or hail
        // ice cloud riming
        if (rime_rate_qc > 0.0) {
            float_t x_i = particle_mean_mass(qi_prime, Ni,
                get_at(cc.ice.constants, Particle_cons_idx::min_x_riming),
                get_at(cc.ice.constants, Particle_cons_idx::max_x));
            float_t d_i = particle_diameter(x_i,
                get_at(cc.ice.constants, Particle_cons_idx::a_geo),
                get_at(cc.ice.constants, Particle_cons_idx::b_geo));
            float_t rime_tmp = res[qc_idx] + qc_prime/dt;
            float_t rime_q = std::max(float_t(0), std::min(rime_rate_qc, rime_tmp));
            rime_tmp = res[Nc_idx] + Nc/dt;
            float_t rime_n = std::max(float_t(0), std::min(rime_rate_nc, rime_tmp));

            // Ice
            res[qi_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;

            // Cloud N
            res[Nc_idx] -= rime_n;
#ifdef TRACE_QC
            if (trace)
                if (abs(rime_q) > 0)
                    std::cout << "traj: " << cc.traj_id
                        << " ice cloud riming dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
#ifdef TRACE_QI
            if (trace)
                std::cout << "traj: " << cc.traj_id << " ice cloud riming dqi " << rime_q << "\n";
#endif
            if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
                float_t mult_1 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
                float_t mult_2 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
                mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
                mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
                // Ice N
                res[Ni_idx] += get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q;
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " ice cloud rimingwith mult dNi "
                              << get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q << "\n";
#endif
            }

            float_t delta_e = latent_heat_melt(
                T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_q
                / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            // Conversion ice -> graupel
            // Technically I had to recalculate x_i given the new qi, Ni
            if (d_i > get_at(cc.constants, Cons_idx::D_conv_ig)) {
                float_t conv_q = rime_q
                    / (get_at(cc.constants, Cons_idx::const5)
                        * (M_PI/6.0 * get_at(cc.constants, Cons_idx::rho_ice) * d_i*d_i*d_i/x_i -1.0));
                conv_q = std::min(qi_prime, conv_q);
                float_t qi_tmp = qi_prime+dt*res[qi_idx];
                x_i = particle_mean_mass(qi_tmp, Ni,
                    get_at(cc.ice.constants, Particle_cons_idx::min_x_conversion),
                    get_at(cc.ice.constants, Particle_cons_idx::max_x));
                float_t tmp = conv_q / std::max(x_i, get_at(cc.constants, Cons_idx::x_conv));
                float_t conv_n = res[Ni_idx] + Ni/dt;
                conv_n = std::min(tmp, conv_n);
                conv_q = conv_n = 0;
                // Ice
                res[qi_idx] -= conv_q;
                // Graupel
                res[qg_idx] += conv_q;
                // Ice N
                res[Ni_idx] -= conv_n;
                // Graupel N
                res[Ng_idx] += conv_n;
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " conv ice->graupel dqi " << -conv_q << ", dNi " << -conv_n << "\n";
#endif
#ifdef TRACE_QG
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " conv ice->graupel dqg " << conv_q << ", dNg " << conv_n << "\n";
#endif
            }
        }

        // ice rain riming
        if (rime_rate_qi > 0.0) {
            float_t rime_tmp = res[qi_idx] + qi_prime/dt;
            float_t rime_qi = std::max(float_t(0), std::min(rime_rate_qi, rime_tmp));
            rime_tmp = res[qr_idx] + qr_prime/dt;
            float_t rime_qr = std::max(float_t(0), std::min(rime_rate_qr, rime_tmp));
            rime_tmp = res[Nr_idx] + Nr/dt;
            float_t rime_tmp2 = res[Ni_idx] + Ni/dt;
            float_t rime_n = std::max(float_t(0), std::min(std::min(rime_rate_nr, rime_tmp), rime_tmp2));

            // Ice
            res[qi_idx] -= rime_qi;
            // Rain
            res[qr_idx] -= rime_qr;
            // Ice N
            res[Ni_idx] -= rime_n;
            // Rain N
            res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " ice rain riming dqr " << -rime_qr << ", dNr " << -rime_n << "\n";
#endif
#ifdef TRACE_QI
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " ice rain riming B dqi " << -rime_qi << ", dNi " << -rime_n << "\n";
#endif
            float_t mult_q = 0.0;
            float_t mult_n = 0.0;
            if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
                float_t mult_1 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
                float_t mult_2 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
                mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
                mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
                mult_n = get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_qr;
                float_t tmp = mult_n*get_at(cc.ice.constants, Particle_cons_idx::min_x_riming);
                mult_q = std::min(rime_qr, tmp);
            }

            float_t delta_e = latent_heat_melt(
                T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_qi
                / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (rime_qi > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] -= delta_e;

            if (T_prime >= get_at(cc.constants, Cons_idx::T_freeze)) {
                float_t qr_tmp = qr_prime+dt*res[qr_idx];
                float_t Nr_tmp = Nr+res[Nr_idx]*dt;
                float_t x_r = particle_mean_mass(
                    qr_tmp, Nr_tmp,
                    get_at(cc.rain.constants, Particle_cons_idx::min_x_riming),
                    get_at(cc.rain.constants, Particle_cons_idx::max_x));
                // Ice
                res[qi_idx] += rime_qi;
                // Rain
                res[qr_idx] += rime_qr;
                // Ice N
                res[Ni_idx] += rime_n;
                // Rain N
                res[Nr_idx] += rime_qr/x_r;
#ifdef TRACE_QR
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " Melting dqr " << rime_qr << ", dNr " << rime_qr/x_r << "\n";
#endif
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " Melting dqi " << rime_qi << ", dNi " << rime_n << "\n";
#endif
                delta_e = latent_heat_melt(
                    T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_qi
                    / specific_heat_dry_air(T_prime);
                // Melting, cooling
                if (rime_qi > 0.0)
                    res[lat_cool_idx] -= delta_e;
                // Freezing, heating
                else
                    res[lat_heat_idx] -= delta_e;
            } else {
                // from multiplication
                // Ice
                res[qi_idx] += mult_q;
                // Ice N
                res[Ni_idx] += mult_n;
                // riming to graupel
                // Graupel
                res[qg_idx] += rime_qi + rime_qr - mult_q;
                // Graupel N
                res[Ng_idx] += rime_n;
#ifdef TRACE_QG
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " Melting with mult dqg " << rime_qi + rime_qr - mult_q << ", dNg " << rime_n << "\n";
#endif
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " Melting with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
            }
        }
    }
}


/**
 * Riming of cloud and rain droplets on snow.
 * 10.1007/s00703-005-0112-4
 * Uses rates calculated via riming_rain_core() and riming_cloud_core() (each based on after Seifert & Beheng; 2006)
 * and breaks up snow if snow gets too big and converts small snowflakes to graupel.
 * Multiplication of ice is based on Beheng, K. D. "Numerical study on the combined action of droplet coagulation,
 * ice particle riming and the splintering process concerning maritime cumuli."
 * Contrib. Atmos. Phys.;(Germany, Federal Republic of) 55.3 (1982)., Eq. 7, conversion to graupel Eqs. 70-71. (2006), Section 3.5.3
 * 10.1007/s00703-005-0112-4
 *
 * @params qc_prime Cloud water mixing ratio
 * @params Nc Number of cloud droplets
 * @params qr_prime Rain mixing ratio
 * @params Nr Number of rain droplets
 * @params qs_prime Snowflakes mixing ratio
 * @params Ns Number of snowflakes
 * @params dep_rate_snow Deposition rate of snow mass
 * @params rime_rate_qc Riming rate for cloud mixing mass
 * @params rime_rate_nc Riming rate for number of cloud droplets
 * @params rime_rate_qr Riming rate for rain mixing mass
 * @params rime_rate_nr Riming rate for number of rain droplets
 * @params rime_rate_qs Riming rate for snow mixing mass
 * @params T_prime Current temperature [K]
 * @params dt Timestep size [s]
 * @params res Vector to store the changes
 * @params cc Model constants
 */
template<class float_t>
void snow_riming(
    float_t &qc_prime,
    float_t &Nc,
    float_t &qr_prime,
    float_t &Nr,
    float_t &qs_prime,
    float_t &Ns,
    float_t &dep_rate_snow,
    float_t &rime_rate_qc,
    float_t &rime_rate_nc,
    float_t &rime_rate_qr,
    float_t &rime_rate_nr,
    float_t &rime_rate_qs,
    float_t &T_prime,
    const double &dt,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc) {

    if (dep_rate_snow > 0.0 && dep_rate_snow >= rime_rate_qc+rime_rate_qr) {
        // Depositional growth is stronger than riming growth, therefore ice stays ice
        // ice cloud riming
        if (rime_rate_qc > 0.0) {
            float_t rime_tmp = res[qc_idx] + qc_prime/dt;
            float_t rime_q = std::max(float_t(0), std::min(rime_tmp, rime_rate_qc));
            rime_tmp = res[Nc_idx] + Nc/dt;
            float_t rime_n = std::max(float_t(0), std::min(rime_tmp, rime_rate_nc));

            // Snow
            res[qs_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;

            // Cloud N
            res[Nc_idx] -= rime_n;
#ifdef TRACE_QC
            if (trace)
                if (abs(rime_q) > 0)
                    std::cout << "traj: " << cc.traj_id
                        << " Snow riming dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
#ifdef TRACE_QS
            if (trace)
                std::cout << "traj: " << cc.traj_id << " Snow riming dqs " << rime_q << "\n";
#endif
            float_t delta_e = latent_heat_melt(
                T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_q
                / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
                float_t mult_1 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
                float_t mult_2 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
                mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
                mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
                float_t mult_n = get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q;
                float_t mult_q = mult_n * get_at(cc.ice.constants, Particle_cons_idx::min_x_riming);
                mult_q = std::min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " Snow riming with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
#ifdef TRACE_QS
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " Snow riming with mult dqs " << -mult_q << "\n";
#endif
            }
        }
        // snow rain riming
        if (rime_rate_qr > 0.0) {
            float_t rime_tmp = res[qr_idx] + qr_prime/dt;
            float_t rime_q = std::max(float_t(0), std::min(rime_rate_qr, rime_tmp));
            rime_tmp = res[Nr_idx] + Nr/dt;
            float_t rime_n = std::max(float_t(0), std::min(rime_tmp, rime_rate_nr));

            // Snow
            res[qs_idx] += rime_q;
            // Rain
            res[qr_idx] -= rime_q;
            // Rain N
            res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " snow rain riming dqr " << -rime_q << ", dNr " << -rime_n << "\n";
#endif
#ifdef TRACE_QS
            if (trace)
                std::cout << "traj: " << cc.traj_id << " Snow rain riming dqs " << rime_q << "\n";
#endif
            float_t delta_e = latent_heat_melt(
                T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_q
                / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
                float_t mult_1 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
                float_t mult_2 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
                mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
                mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
                float_t mult_n = get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q;
                float_t mult_q = mult_n * get_at(cc.ice.constants, Particle_cons_idx::min_x_riming);
                mult_q = std::min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " snow rain riming with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
#ifdef TRACE_QS
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " snow rain riming with mult dqs " << -mult_q << "\n";
#endif
            }
        }
    } else {
        // Depositional growth negative or smaller than riming growth,
        // therefore snow is allowed to convert to graupel and / or hail
        // snow cloud riming
        if (rime_rate_qc > 0.0) {
            float_t x_s = particle_mean_mass(qs_prime, Ns,
                get_at(cc.snow.constants, Particle_cons_idx::min_x_riming),
                get_at(cc.snow.constants, Particle_cons_idx::max_x));
            float_t d_s = particle_diameter(x_s,
                get_at(cc.snow.constants, Particle_cons_idx::a_geo),
                get_at(cc.snow.constants, Particle_cons_idx::b_geo));
            float_t rime_tmp = res[qc_idx] + qc_prime/dt;
            float_t rime_q = std::max(float_t(0), std::min(rime_rate_qc, rime_tmp));
            rime_tmp = res[Nc_idx] + Nc/dt;
            float_t rime_n = std::max(float_t(0), std::min(rime_rate_nc, rime_tmp));

            // Snow
            res[qs_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;

            // Cloud N
            res[Nc_idx] -= rime_n;
#ifdef TRACE_QC
            if (trace)
                if (abs(rime_q) > 0)
                    std::cout << "traj: " << cc.traj_id
                        << " snow depos dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
#ifdef TRACE_QS
            if (trace)
                std::cout << "traj: " << cc.traj_id << " snow depos dqs " << rime_q << "\n";
#endif
            float_t delta_e = latent_heat_melt(
                T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_q
                / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            float_t mult_q = 0.0;
            if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
                float_t mult_1 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
                float_t mult_2 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
                mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
                mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
                float_t mult_n = get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q;
                mult_q = mult_n * get_at(cc.ice.constants, Particle_cons_idx::min_x_riming);
                mult_q = std::min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " snow depos with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
#ifdef TRACE_QS
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " snow depos with mult dqs " << -mult_q << "\n";
#endif
            }

            // Conversion snow -> graupel
            if (d_s > get_at(cc.constants, Cons_idx::D_conv_sg)) {
                float_t conv_q = (rime_q - mult_q)
                    / (get_at(cc.constants, Cons_idx::const5)*(M_PI/6.0
                    * get_at(cc.constants, Cons_idx::rho_ice) * d_s*d_s*d_s/x_s -1.0));
                float_t conv_tmp = res[qs_idx] + qs_prime/dt;
                conv_q = std::max(float_t(0), std::min(conv_tmp, conv_q));

                x_s = particle_mean_mass(qs_prime, Ns,
                    get_at(cc.snow.constants, Particle_cons_idx::min_x_riming),
                    get_at(cc.snow.constants, Particle_cons_idx::max_x));
                float_t tmp = conv_q / std::max(x_s, get_at(cc.constants, Cons_idx::x_conv));
                conv_tmp = res[Ns_idx] + Ns/dt;
                float_t conv_n = std::max(float_t(0), std::min(tmp, conv_tmp));

                // Snow
                res[qs_idx] -= conv_q;
                // Graupel
                res[qg_idx] += conv_q;
                // Snow N
                res[Ns_idx] -= conv_n;
                // Graupel N
                res[Ng_idx] += conv_n;
#ifdef TRACE_QS
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " conversion snow->graupel dqs " << -conv_q << ", dNs " << -conv_n << "\n";
#endif
#ifdef TRACE_QG
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " conversion snow->graupel dqg " << conv_q << ", dNg " << conv_n << "\n";
#endif
            }
        }

        // Snow rain riming
        if (rime_rate_qs > 0.0) {
            float_t rime_tmp = res[qs_idx] + qs_prime/dt;
            float_t rime_qs = std::max(float_t(0), std::min(rime_rate_qs, rime_tmp));
            rime_tmp = res[qr_idx] + qr_prime/dt;
            float_t rime_qr = std::max(float_t(0), std::min(rime_rate_qr, rime_tmp));
            rime_tmp = res[Nr_idx] + Nr/dt;
            float_t rime_tmp2 = res[Ns_idx] + Ns/dt;
            float_t rime_n = std::max(float_t(0), std::min(std::min(rime_rate_nr,
                rime_tmp), rime_tmp2));

            // Snow
            res[qs_idx] -= rime_qs;
            // Rain
            res[qr_idx] -= rime_qr;
            // Snow N
            res[Ns_idx] -= rime_n;
            // Rain N
            res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " snow rain riming 2 dqr " << -rime_qr << ", dNr " << -rime_n << "\n";
#endif
#ifdef TRACE_QS
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " snow rain riming 2 dqs " << -rime_qs << ", dNs " << -rime_n << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_qr
                / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (rime_qr < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            float_t mult_q = 0.0;
            float_t mult_n = 0.0;
            if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
                float_t mult_1 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
                float_t mult_2 =
                    (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
                mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
                mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
                mult_n = get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_qr;
                float_t tmp = mult_n*get_at(cc.ice.constants, Particle_cons_idx::min_x_riming);
                mult_q = std::min(rime_qr, tmp);
            }
            if (T_prime >= get_at(cc.constants, Cons_idx::T_freeze)) {
                float_t qr_tmp = qr_prime+dt*res[qr_idx];
                float_t Nr_tmp = Nr+res[Nr_idx]*dt;
                float_t x_r = particle_mean_mass(
                    qr_tmp, Nr_tmp,
                    get_at(cc.rain.constants, Particle_cons_idx::min_x_riming),
                    get_at(cc.rain.constants, Particle_cons_idx::max_x));

                // Snow
                res[qs_idx] += rime_qs;
                // Rain
                res[qr_idx] += rime_qr;
                // Snow N
                res[Ns_idx] += rime_n;
                // Rain N
                res[Nr_idx] += rime_qr/x_r;
#ifdef TRACE_QR
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " More melting dqr " << rime_qr << ", dNr " << rime_qr/x_r << "\n";
#endif
#ifdef TRACE_QS
                if (trace)
                    std::cout << "traj: " << cc.traj_id
                        << " More melting dqs " << rime_qs << ", dNs " << rime_n << "\n";
#endif
                float_t delta_e_enhance = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_qr
                                            / specific_heat_dry_air(T_prime);
                // Melting, cooling
                if (rime_qr > 0.0)
                    res[lat_cool_idx] -= delta_e_enhance;
                // Freezing, heating
                else
                    res[lat_heat_idx] -= delta_e_enhance;
            } else {
                // from multiplication
                // Ice
                res[qi_idx] += mult_q;
                // Ice N
                res[Ni_idx] += mult_n;
                // riming to graupel
                // Graupel
                res[qg_idx] += rime_qs + rime_qr - mult_q;
                // Graupel N
                res[Ng_idx] += rime_n;
#ifdef TRACE_QI
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " More melting with mult dqi " << mult_q
                              << ", dNi " << mult_n << "\n";
#endif
#ifdef TRACE_QG
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " More melting with mult dqg " << rime_qs + rime_qr - mult_q
                              << ", dNg " << rime_n << "\n";
#endif
            }
        }
    }
}


/**
 * Riming of cloud droplets with graupel or hail after Seifert & Beheng (2006), Eqs. 61-67.
 * 10.1007/s00703-005-0112-4
 * Melting of graupel and hail is calculated with Rutledge, Steven A., and Peter V. Hobbs.
 * "The mesoscale and microscale structure and organization of clouds and precipitation in midlatitude cyclones.
 * XII: A diagnostic modeling study of precipitation development in narrow cold-frontal rainbands."
 * Journal of Atmospheric Sciences 41.20 (1984): 2949-2972., Eq. A22.
 * Applies ice multiplication via Hallet-Mossop process.
 * Multiplication of ice is based on Beheng, K. D. "Numerical study on the combined action of droplet coagulation,
 * ice particle riming and the splintering process concerning maritime cumuli."
 * Contrib. Atmos. Phys.;(Germany, Federal Republic of) 55.3 (1982)., Eq. 7.
 *
 * @params qc_prime Cloud water mixing ratio
 * @params Nc Number of cloud droplets
 * @params T_prime Current temperature [K]
 * @params q1 Mixing ratio of any ice particle
 * @params N1 Particle number of any ice particle
 * @params resq On out: difference of the ice particle mixing ratio
 * @params resn On out: difference of the ice particle number
 * @params coeffs Model constants for collision processes
 * @params pc1 Model constants for a particle
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void particle_cloud_riming(
    float_t &qc_prime,
    float_t &Nc,
    float_t &T_prime,
    float_t &q1,
    float_t &N1,
    float_t &resq,
    float_t &resn,
    collection_model_constants_t<float_t> &coeffs,
    particle_model_constants_t<float_t> &pc1,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    float_t const1 = get_at(cc.constants, Cons_idx::const0)
        * get_at(pc1.constants, Particle_cons_idx::ecoll_c);
    float_t x_c = particle_mean_mass(qc_prime, Nc, get_at(cc.cloud.constants, Particle_cons_idx::min_x_riming),
        get_at(cc.cloud.constants, Particle_cons_idx::max_x));
    float_t d_c = particle_diameter(x_c, get_at(cc.cloud.constants, Particle_cons_idx::a_geo),
        get_at(cc.cloud.constants, Particle_cons_idx::b_geo));
    float_t x_1 = particle_mean_mass(q1, N1, get_at(pc1.constants, Particle_cons_idx::min_x_riming),
        get_at(pc1.constants, Particle_cons_idx::max_x));
    float_t d_1 = particle_diameter(x_1, get_at(pc1.constants, Particle_cons_idx::a_geo),
        get_at(pc1.constants, Particle_cons_idx::b_geo));

    if (qc_prime > get_at(cc.cloud.constants, Particle_cons_idx::q_crit_c)
        && q1 > get_at(pc1.constants, Particle_cons_idx::q_crit_c)
        && d_1 > get_at(pc1.constants, Particle_cons_idx::d_crit_c)
        && d_c > get_at(cc.cloud.constants, Particle_cons_idx::d_crit_c)) {
        float_t v_1 = particle_velocity(x_1, get_at(pc1.constants, Particle_cons_idx::a_vel),
            get_at(pc1.constants, Particle_cons_idx::b_vel)) * get_at(pc1.constants, Particle_cons_idx::rho_v);
        float_t v_c = particle_velocity(
            x_c, get_at(cc.cloud.constants, Particle_cons_idx::a_vel),
            get_at(cc.cloud.constants, Particle_cons_idx::b_vel))
            * get_at(cc.cloud.constants, Particle_cons_idx::rho_v);
        float_t coll_tmp1 = get_at(pc1.constants, Particle_cons_idx::ecoll_c)/dt;
        float_t coll_tmp2 = get_at(cc.constants, Cons_idx::ecoll_min)/dt;
        float_t e_coll_n = std::min(coll_tmp1,
            std::max(const1*(d_c-get_at(cc.cloud.constants, Particle_cons_idx::d_crit_c)),
                coll_tmp2));
        float_t e_coll_q = e_coll_n;

        float_t rime_n = M_PI/4.0 * e_coll_n * N1 * Nc
            * (coeffs.delta_n_aa * d_1 * d_1
                + coeffs.delta_n_ab * d_1 * d_c
                + coeffs.delta_n_bb * d_c * d_c)
            * sqrt(coeffs.theta_n_aa * v_1 * v_1
                - coeffs.theta_n_ab * v_1 * v_c
                + coeffs.theta_n_bb * v_c * v_c);
        float_t rime_q = M_PI/4.0 * e_coll_q * N1 * qc_prime
            * (coeffs.delta_q_aa * d_1 * d_1
                + coeffs.delta_q_ab * d_1 * d_c
                + coeffs.delta_q_bb * d_c * d_c)
            * sqrt(coeffs.theta_q_aa * v_1 * v_1
                - coeffs.theta_q_ab * v_1 * v_c
                + coeffs.theta_q_bb * v_c * v_c);
        float_t rime_tmp = res[qc_idx] + qc_prime/dt;
        rime_q = std::max(float_t(0), std::min(rime_tmp, rime_q));
        rime_tmp = res[Nc_idx] + Nc/dt;
        rime_n = std::max(float_t(0), std::min(rime_tmp, rime_n));
        resq += rime_q;
        // Cloud
        res[qc_idx] -= rime_q;

        // Cloud N
        res[Nc_idx] -= rime_n;
#ifdef TRACE_QC
        if (trace)
            if (abs(rime_q) > 0)
                std::cout << "traj: " << cc.traj_id
                    << " Particle cloud riming dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_q
            / specific_heat_dry_air(T_prime);
        // Sublimination, cooling
        if (rime_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] += delta_e;

        // ice multiplication
        if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
            float_t mult_1 =
                (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
            float_t mult_2 =
                (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
            mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
            mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
            float_t mult_n = get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q;
            float_t mult_q = mult_n * get_at(cc.ice.constants, Particle_cons_idx::min_x_riming);
            mult_q = std::min(rime_q, mult_q);

            // Ice
            res[qi_idx] += mult_q;
            // Ice N
            res[Ni_idx] += mult_n;
            resq -= mult_q;
#ifdef TRACE_QI
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " Particle cloud riming dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
        }

        // Enhancement of melting
        if (T_prime > get_at(cc.constants, Cons_idx::T_freeze) && enhanced_melting) {
            float_t tmp_const = specific_heat_water(T_prime)
                / latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze));
            float_t melt_q = (T_prime-get_at(cc.constants, Cons_idx::T_freeze))
                * tmp_const*rime_q;
            float_t melt_n = melt_q/x_1;
            float_t melt_tmp = resq + q1/dt;
            melt_q = std::max(float_t(0), std::min(melt_tmp, melt_q));
            melt_tmp = resn + N1/dt;
            melt_n = std::max(float_t(0), std::min(melt_tmp, melt_n));

            resq -= melt_q;
            resn -= melt_n;
            // Rain
            res[qr_idx] += melt_q;
            // Rain N
            res[Nr_idx] += melt_n;
#ifdef TRACE_QR
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " enhancement of melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
            float_t delta_e_enhance = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * melt_q
                                        / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (melt_q > 0.0)
                res[lat_cool_idx] -= delta_e_enhance;
            // Freezing, heating
            else
                res[lat_heat_idx] -= delta_e_enhance;
        }
    }
}


/**
 * Riming of rain droplets with graupel or hail after Seifert & Beheng (2006), Eqs. 61-63.
 * 10.1007/s00703-005-0112-4
 * Melting of graupel and hail is calculated with Rutledge, Steven A., and Peter V. Hobbs.
 * "The mesoscale and microscale structure and organization of clouds and precipitation in midlatitude cyclones.
 * XII: A diagnostic modeling study of precipitation development in narrow cold-frontal rainbands."
 * Journal of Atmospheric Sciences 41.20 (1984): 2949-2972., Eq. A21.
 * Applies ice multiplication via Hallet-Mossop process.
 * Multiplication of ice is based on Beheng, K. D. "Numerical study on the combined action of droplet coagulation,
 * ice particle riming and the splintering process concerning maritime cumuli."
 * Contrib. Atmos. Phys.;(Germany, Federal Republic of) 55.3 (1982)., Eq. 7.
 *
 * @params qr_prime Rain mixing ratio
 * @params Nr Number of rain droplets
 * @params T_prime Temperature [K]
 * @params q1 Mixing ratio of any ice particle
 * @params N1 Particle number of any ice particle
 * @params resq On out: difference of the ice particle mixing ratio
 * @params resn On out: difference of the ice particle number
 * @params coeffs Model constants for collision processes
 * @params pc1 Model constants for a particle
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void particle_rain_riming(
    float_t &qr_prime,
    float_t &Nr,
    float_t &T_prime,
    float_t &q1,
    float_t &N1,
    float_t &resq,
    float_t &resn,
    collection_model_constants_t<float_t> &coeffs,
    particle_model_constants_t<float_t> &pc1,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (qr_prime > get_at(cc.constants, Cons_idx::q_crit) && q1 > get_at(cc.constants, Cons_idx::q_crit)) {
        float_t x_r = particle_mean_mass(
            qr_prime, Nr, get_at(cc.rain.constants, Particle_cons_idx::min_x_riming),
            get_at(cc.rain.constants, Particle_cons_idx::max_x));
        float_t d_r = particle_diameter(
            x_r, get_at(cc.rain.constants, Particle_cons_idx::a_geo),
            get_at(cc.rain.constants, Particle_cons_idx::b_geo));
        float_t x_1 = particle_mean_mass(
            q1, N1, get_at(pc1.constants, Particle_cons_idx::min_x_riming),
            get_at(pc1.constants, Particle_cons_idx::max_x));
        float_t d_1 = particle_diameter(
            x_1, get_at(pc1.constants, Particle_cons_idx::a_geo),
            get_at(pc1.constants, Particle_cons_idx::b_geo));

        float_t v_1 = particle_velocity(
            x_1, get_at(pc1.constants, Particle_cons_idx::a_vel),
            get_at(pc1.constants, Particle_cons_idx::b_vel)) * get_at(pc1.constants, Particle_cons_idx::rho_v);
        float_t v_r = particle_velocity(
            x_r, get_at(cc.rain.constants, Particle_cons_idx::a_vel),
            get_at(cc.rain.constants, Particle_cons_idx::b_vel)) * get_at(cc.rain.constants, Particle_cons_idx::rho_v);

        float_t rime_n = M_PI/4.0 * N1 * Nr
            * (coeffs.delta_n_aa * d_1 * d_1
                + coeffs.delta_n_ab * d_1 * d_r
                + coeffs.delta_n_bb * d_r * d_r)
            * sqrt(coeffs.theta_n_aa * v_1 * v_1
                - coeffs.theta_n_ab * v_1 * v_r
                + coeffs.theta_n_bb * v_r * v_r);
        float_t rime_q = M_PI/4.0 * N1 * qr_prime
            * (coeffs.delta_n_aa * d_1 * d_1
                + coeffs.delta_q_ab * d_1 * d_r
                + coeffs.delta_q_bb * d_r * d_r)
            * sqrt(coeffs.theta_n_aa * v_1 * v_1
                - coeffs.theta_q_ab * v_1 * v_r
                + coeffs.theta_q_bb * v_r * v_r);
        float_t rime_tmp = res[qr_idx] + qr_prime/dt;
        rime_q = std::max(float_t(0), std::min(rime_tmp, rime_q));
        rime_tmp = res[Nr_idx] + Nr/dt;
        rime_n = std::max(float_t(0), std::min(rime_tmp, rime_n));
        resq += rime_q;
        // Rain
        res[qr_idx] -= rime_q;
        // Rain N
        res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id
                << " particle rain riming dqr " << -rime_q << ", dNr " << -rime_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * rime_q
            / specific_heat_dry_air(T_prime);
        // Melting, cooling
        if (rime_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] += delta_e;

        // ice multiplication (Hallett and Mossop)
        if (T_prime < get_at(cc.constants, Cons_idx::T_freeze) && ice_multiplication) {
            float_t mult_1 =
                (T_prime - get_at(cc.constants, Cons_idx::T_mult_min))*get_at(cc.constants, Cons_idx::const3);
            float_t mult_2 =
                (T_prime - get_at(cc.constants, Cons_idx::T_mult_max))*get_at(cc.constants, Cons_idx::const4);
            mult_1 = std::max(float_t(0.0), std::min(mult_1, float_t(1.0)));
            mult_2 = std::max(float_t(0.0), std::min(mult_2, float_t(1.0)));
            float_t mult_n = get_at(cc.constants, Cons_idx::C_mult) * mult_1 * mult_2 * rime_q;
            float_t mult_q = mult_n * get_at(cc.ice.constants, Particle_cons_idx::min_x_riming);
            mult_q = std::min(rime_q, mult_q);

            // Ice
            res[qi_idx] += mult_q;
            // Ice N
            res[Ni_idx] += mult_n;
#ifdef TRACE_QI
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " particle rain riming dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
            resq -= mult_q;
        }

        // Enhancement of melting
        if (T_prime > get_at(cc.constants, Cons_idx::T_freeze) && enhanced_melting) {
            float_t tmp_const = specific_heat_water(T_prime)
                / latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze));
            float_t melt_q = (T_prime-get_at(cc.constants, Cons_idx::T_freeze))
                * tmp_const*rime_q;
            float_t melt_n = melt_q/x_1;
            float_t melt_tmp = resq + q1/dt;
            melt_q = std::max(float_t(0), std::min(melt_tmp, melt_q));
            melt_tmp = resn + N1/dt;
            melt_n = std::max(float_t(0), std::min(melt_tmp, melt_n));

            resq -= melt_q;
            resn -= melt_n;
            // Rain
            res[qr_idx] += melt_q;
            // Rain N
            res[Nr_idx] += melt_n;
#ifdef TRACE_QR
            if (trace)
                std::cout << "traj: " << cc.traj_id
                    << " particle rain riming enhancement dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
            float_t delta_e_enhance = latent_heat_melt(
                T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * melt_q
                / specific_heat_dry_air(T_prime);
            // Melting, cooling
            if (melt_q > 0.0)
                res[lat_cool_idx] -= delta_e_enhance;
            // Freezing, heating
            else
                res[lat_heat_idx] -= delta_e_enhance;
        }
    }
}


/**
 * Freezing of rain and conversion to ice, graupel, hail after Seifert & Beheng (2006), Section 3.4.
 * 10.1007/s00703-005-0112-4
 * The heterogeneous freezing follows a stochastic model where the rate follows an exponential form that depends
 * on the temperature.
 * Small amounts of rain freeze instantaneously.
 *
 * @params qr_prime Rain mixing ratio
 * @params Nr Number of rain droplets
 * @params T_prime Temperature [K]
 * @params dt Timestep size [s]
 * @params res Vector to store the changes
 * @params cc Model constants
 */
template<class float_t>
void rain_freeze(
    float_t &qr_prime,
    float_t &Nr,
    float_t &T_prime,
    const double &dt,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc) {

    if (T_prime < get_at(cc.constants, Cons_idx::T_freeze)) {
        float_t xmax_ice = pow(pow(get_at(cc.constants, Cons_idx::D_rainfrz_ig)
            / get_at(cc.rain.constants, Particle_cons_idx::a_geo),
            1.0/get_at(cc.rain.constants, Particle_cons_idx::b_geo)), get_at(cc.rain.constants, Particle_cons_idx::mu));
        float_t xmax_gr  = pow(pow(get_at(cc.constants, Cons_idx::D_rainfrz_gh)
            / get_at(cc.rain.constants, Particle_cons_idx::a_geo),
            1.0/get_at(cc.rain.constants, Particle_cons_idx::b_geo)), get_at(cc.rain.constants, Particle_cons_idx::mu));
        float_t fr_q, fr_n, fr_n_i, fr_q_i, fr_n_g, fr_q_g, fr_n_h,
            fr_q_h, fr_n_tmp, fr_q_tmp;
        fr_q = fr_n = fr_n_i = fr_q_i = fr_n_g = fr_q_g = fr_n_h =
            fr_q_h = fr_n_tmp = fr_q_tmp = 0.0;
        float_t Nr_tmp = Nr;
        if (qr_prime <= get_at(cc.constants, Cons_idx::q_crit_fr)) {
            if (T_prime < get_at(cc.constants, Cons_idx::T_f)) {
                // instantaneous freezing
                fr_q = fr_q_i = qr_prime/dt;
                fr_n = fr_n_i = Nr/dt;
                fr_n_tmp = fr_q_tmp = 1.0;
            }
        } else {
            float_t x_r = particle_mean_mass(qr_prime, Nr,
                get_at(cc.rain.constants, Particle_cons_idx::min_x_freezing),
                get_at(cc.rain.constants, Particle_cons_idx::max_x));
            if (T_prime < get_at(cc.constants, Cons_idx::T_f)) {
                fr_q = qr_prime;
                fr_n = Nr_tmp;
                float_t lam = pow(get_at(cc.rain.constants, Particle_cons_idx::g1)
                    / get_at(cc.rain.constants, Particle_cons_idx::g2)*x_r,
                        - get_at(cc.rain.constants, Particle_cons_idx::mu));
                float_t N0 = get_at(cc.rain.constants, Particle_cons_idx::mu) * Nr_tmp
                    * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm1))
                    / get_at(cc.rain.constants, Particle_cons_idx::g1);
                float_t tmp = lam*xmax_ice;
                fr_n_i = N0/(get_at(cc.rain.constants, Particle_cons_idx::mu)
                    * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm1))) * cc.table_r1.look_lo(tmp);
                fr_q_i = N0/(get_at(cc.rain.constants, Particle_cons_idx::mu)
                    * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm2))) * cc.table_r2.look_lo(tmp);
                tmp = lam*xmax_gr;
                fr_n_g = N0/(get_at(cc.rain.constants, Particle_cons_idx::mu)
                    * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm1))) * cc.table_r1.look_lo(tmp);
                fr_q_g = N0/(get_at(cc.rain.constants, Particle_cons_idx::mu)
                    * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm2))) * cc.table_r2.look_lo(tmp);
                fr_n_h = fr_n - fr_n_g;
                fr_q_h = fr_q - fr_q_g;
                fr_n_g = fr_n_g - fr_n_i;
                fr_q_g = fr_q_g - fr_q_i;
                float_t fr_tmp = (res[Nr_idx] + Nr_tmp)*dt;
                if (Nr_tmp == 0 || (fr_n == fr_tmp == 0)) {
                    fr_n_tmp = 0;
                } else {
                    fr_n_tmp = Nr_tmp / std::max(fr_n, fr_tmp);
                }
                fr_tmp = (res[qr_idx] + qr_prime)*dt;
                fr_q_tmp = qr_prime/std::max(fr_q, fr_tmp);
            } else {
                // heterogeneous freezing
                float_t j_het = std::max(get_at(cc.constants, Cons_idx::b_HET) *
                    (exp(get_at(cc.constants, Cons_idx::a_HET)
                        * (get_at(cc.constants, Cons_idx::T_freeze)-T_prime)) - 1.0),
                    0.0) / get_at(cc.constants, Cons_idx::rho_w);

                if (j_het >= 1.0e-20/dt) {
                    fr_n = j_het * qr_prime;
                    fr_q = j_het * qr_prime * x_r * get_at(cc.rain.constants, Particle_cons_idx::c_z);  // rain_coeffs
                    float_t lam = pow(get_at(cc.rain.constants, Particle_cons_idx::g1)
                        / get_at(cc.rain.constants, Particle_cons_idx::g2) * x_r,
                            -get_at(cc.rain.constants, Particle_cons_idx::mu));
                    float_t N0 = get_at(cc.rain.constants, Particle_cons_idx::mu)
                        * Nr_tmp * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm1))
                        / get_at(cc.rain.constants, Particle_cons_idx::g1);

                    float_t tmp = lam*xmax_ice;
                    fr_n_i = j_het * N0/(get_at(cc.rain.constants, Particle_cons_idx::mu)
                        * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm2))) * cc.table_r2.look_lo(tmp);
                    fr_q_i = j_het * N0/(get_at(cc.rain.constants, Particle_cons_idx::mu)
                        * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm3))) * cc.table_r3.look_lo(tmp);

                    tmp = lam*xmax_gr;
                    fr_n_g = j_het * N0/(get_at(cc.rain.constants, Particle_cons_idx::mu)
                        * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm2))) * cc.table_r2.look_lo(tmp);
                    fr_q_g = j_het * N0/(get_at(cc.rain.constants, Particle_cons_idx::mu)
                        * pow(lam, get_at(cc.rain.constants, Particle_cons_idx::nm3))) * cc.table_r3.look_lo(tmp);

                    fr_n_h = fr_n - fr_n_g;
                    fr_q_h = fr_q - fr_q_g;
                    fr_n_g = fr_n_g - fr_n_i;
                    fr_q_g = fr_q_g - fr_q_i;
                    float_t fr_tmp = (res[Nr_idx] + Nr_tmp)*dt;
                    fr_n_tmp = Nr_tmp/std::max(fr_n, fr_tmp);
                    fr_tmp = (res[qr_idx] + qr_prime)*dt;
                    fr_q_tmp = qr_prime/std::max(fr_q, fr_tmp);
                } else {
                    fr_n = fr_q = fr_n_i = fr_q_i = fr_n_g = fr_q_g
                        = fr_n_h = fr_q_h = fr_n_tmp = fr_q_tmp = 0.0;
                }
#if defined(TRACE_QR) || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
                if (trace)
                    std::cout << "traj: " << cc.traj_id << " Hetero j_het: " << j_het
                              << ", val: " << get_at(cc.constants, Cons_idx::b_HET)
                                * (exp(get_at(cc.constants, Cons_idx::a_HET)
                                    * (get_at(cc.constants, Cons_idx::T_freeze)-T_prime)) - 1.0)
                              << "\n";
#endif
            }
            fr_n = fr_n * fr_n_tmp;
            fr_q = fr_q * fr_q_tmp;
            fr_n_i = fr_n_i * fr_n_tmp;
            fr_n_g = fr_n_g * fr_n_tmp;
            fr_n_h = fr_n_h * fr_n_tmp;
            fr_q_i = fr_q_i * fr_q_tmp;
            fr_q_g = fr_q_g * fr_q_tmp;
            fr_q_h = fr_q_h * fr_q_tmp;
        }
        // Rain
        res[qr_idx] -= fr_q;
        // Rain N
        res[Nr_idx] -= fr_n;
        // res[Nr_idx] = Nr_tmp - fr_n;
#if defined(TRACE_QR) || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
        if (trace)
            std::cout << "traj: " << cc.traj_id
                      << " qr_prime <= q_crit_fr " << (qr_prime <= get_at(cc.constants, Cons_idx::q_crit_fr))
                      << "\nT_prime < get_at(cc.constants, Cons_idx::T_freeze) "
                      << (T_prime < get_at(cc.constants, Cons_idx::T_freeze))
                      << "\nT_prime < get_at(cc.constants, Cons_idx::T_f) "
                      << (T_prime < get_at(cc.constants, Cons_idx::T_f))
                      << "\nqr_prime <= get_at(cc.constants, Cons_idx::q_crit_fr) "
                      << (qr_prime <= get_at(cc.constants, Cons_idx::q_crit_fr))
                      << "\nfr_q " << -fr_q
                      << "\nfr_q_i " << fr_q_i
                      << "\nfr_q_tmp " << fr_q_tmp
                      << "\n";
#endif
#ifdef TRACE_QR
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Freezing dqr " << -fr_q << ", dNr " << -fr_n << "\n";
#endif

        // Snow
        res[qs_idx] += fr_q_i;
        // Snow N
        res[Ns_idx] += fr_n_i;
#ifdef TRACE_QS
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Freezing dqs " << fr_q_i << ", dNs " << fr_n_i << "\n";
#endif

        // Graupel
        res[qg_idx] += fr_q_g;
        // Graupel N
        res[Ng_idx] += fr_n_g;
#ifdef TRACE_QG
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Freezing dqg " << fr_q_g << ", dNg " << fr_n_g << "\n";
#endif

        // Hail
        res[qh_idx] += fr_q_h;
        // Hail N
        res[Nh_idx] += fr_n_h;
#ifdef TRACE_QH
        if (trace)
            std::cout << "traj: " << cc.traj_id << " Freezing dqh " << fr_q_h << ", dNh " << fr_n_h << "\n";
#endif

        float_t delta_e = latent_heat_melt(T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * fr_q
                                    / specific_heat_dry_air(T_prime);
        // Melting, cooling
        if (fr_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] += delta_e;
    }
}


/**
 * Ice melts instantaneously to rain or cloud droplets here.
 *
 * @params qi_prime Ice crystal mixing ratio
 * @params Ni Number of ice crystals
 * @params T_prime Current temperature [K]
 * @params dt Timestep size [s]
 * @params res Vector to store the changes
 * @params cc Model constants
 * @params dt Time step size [s]
 */
template<class float_t>
void ice_melting(
    float_t &qi_prime,
    float_t &Ni,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t<float_t> &cc,
    const double &dt) {

    if (T_prime > get_at(cc.constants, Cons_idx::T_freeze) && qi_prime > 0.0) {
        float_t x_i = particle_mean_mass(
            qi_prime, Ni, get_at(cc.ice.constants, Particle_cons_idx::min_x_melt),
            get_at(cc.ice.constants, Particle_cons_idx::max_x));
        // Complete melting

        float_t melt_q = qi_prime / dt;  // Instantanuous, hence "/dt"
        float_t melt_n = Ni / dt;

        // Ice
        res[qi_idx] = -melt_q;
        // Ice N
        res[Ni_idx] = -melt_n;
#ifdef TRACE_QI
        if (trace) {
            std::cout << "traj: " << cc.traj_id << " ice_melting dqi " << -melt_q << ", dNi " << -melt_n << "\n";
        }
#endif
        // melt into cloud or rain
        if (x_i > get_at(cc.cloud.constants, Particle_cons_idx::max_x)) {
            // Rain
            res[qr_idx] += melt_q;
            // Rain N
            res[Nr_idx] += melt_n;
#ifdef TRACE_QR
            if (trace)
                std::cout << "traj: " << cc.traj_id << " ice melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
        } else {
            // Cloud
            res[qc_idx] += melt_q;

            // Cloud N
            res[Nc_idx] += melt_n;
#ifdef TRACE_QC
            if (trace)
                if (abs(melt_q) > 0)
                    std::cout << "traj: " << cc.traj_id << " ice melting dqc " << melt_q << ", dNc " << melt_n << "\n";
#endif
        }

        float_t delta_e = latent_heat_melt(
            T_prime, get_at(cc.constants, Cons_idx::T_freeze)) * melt_q
            / specific_heat_dry_air(T_prime);
        // Melting, cooling
        if (melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}

/**
 * This function evaluates the RHS function of the ODE. It uses the 2 moment
 * cloud scheme after Seifert and Beheng (2006),
 * see https://doi.org/10.1007/s00703-005-0112-4
 *
 * @param res On out: system state difference
 * @param y Old system state
 * @param ref Pointer to reference quantities to transform between units
 * @param cc Pointer to constants from the model
 * @param dt Timestep size
 * @param fixed If True: Reset change of pressure, temperature and ascent (w)
 *              at the end to zero
 */
template<class float_t>
void RHS_SB(std::vector<float_t> &res,
    std::vector<float_t> &y,
    const reference_quantities_t &ref,
    model_constants_t<float_t> &cc,
    const double &dt,
    bool fixed = false) {
    // // Decrypt the variables
    float_t p = y[p_idx];
    float_t T = y[T_idx];
    float_t w = y[w_idx];
    float_t z = y[z_idx];

    float_t S = y[S_idx];
    float_t qc = y[qc_idx];
    float_t qr = y[qr_idx];
    float_t qv = y[qv_idx];
    float_t Nc = y[Nc_idx];
    float_t Nr = y[Nr_idx];
    float_t qi = y[qi_idx];
    float_t Ni = y[Ni_idx];
    float_t qs = y[qs_idx];
    float_t Ns = y[Ns_idx];
    float_t qg = y[qg_idx];
    float_t Ng = y[Ng_idx];
    float_t qh = y[qh_idx];
    float_t Nh = y[Nh_idx];
    float_t n_inact = y[n_inact_idx];

    for (auto &r : res) r = 0;

    // Safety measure: ensure positiveness
    if (0. > qc)
        qc = 0.0;

    if (0. > qr)
        qr = 0.0;

    if (0. > qv)
        qv = 0.0;

    if (0. > Nc)
        Nc = 0.0;

    if (0. > qi)
        qi = 0.0;

    if (0. > qs)
        qs = 0.0;

    if (0. > qg)
        qg = 0.0;

    if (0. > qh)
        qh = 0.0;

    if (0. > Nr)
        Nr = 0.0;

    if (0. > Ni)
        Ni = 0.0;

    if (0. > Ns)
        Ns = 0.0;

    if (0. > Ng)
        Ng = 0.0;

    if (0. > Nh)
        Nh = 0.0;

    float_t dep_rate_ice = 0.0;
    float_t dep_rate_snow = 0.0;
    // Change to dimensional variables
    float_t p_prime = ref.pref * p;
    float_t T_prime = ref.Tref * T;
    float_t w_prime = ref.wref * w;
    float_t qc_prime = ref.qref * qc;
    float_t qr_prime = ref.qref * qr;
    float_t qv_prime = ref.qref * qv;
    float_t qi_prime = ref.qref * qi;
    float_t qs_prime = ref.qref * qs;
    float_t qg_prime = ref.qref * qg;
    float_t qh_prime = ref.qref * qh;
    float_t z_prime = ref.zref * z;
    // Additional variables such as super saturation
    float_t T_c = T_prime - get_at(cc.constants, Cons_idx::T_freeze);
    float_t p_sat = saturation_pressure_water(
        T_prime, get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b));
    float_t p_sat_ice = saturation_pressure_ice(
        T_prime, get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_ice_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_ice_const_b));
    float_t D_vtp = diffusivity(T_prime, p_prime);
    S = convert_qv_to_S(
        p_prime,
        T_prime,
        qv_prime,
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b),
        get_at(cc.constants, Cons_idx::Epsilon));
    float_t e_d = compute_pv(
        T_prime, S, get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b));

    float_t S_i =
        (T_prime < get_at(cc.constants, Cons_idx::T_freeze)) ? e_d / p_sat_ice : float_t(1);

    float_t s_sw = S - 1.0;   // super saturation over water
    float_t s_si = S_i - 1.0;  // super saturation over ice
#if defined(TRACE_SAT) || defined(TRACE_QI) || defined(TRACE_QS) || defined(TRACE_QV) || defined(TRACE_QR) \
    || defined(TRACE_QC) || defined(TRACE_QG) || defined(TRACE_QH)
    if (trace) {
        float_t x = particle_mean_mass(
            qr_prime, Nr, get_at(cc.rain.constants, Particle_cons_idx::min_x_depo),
            get_at(cc.rain.constants, Particle_cons_idx::max_x));
        float_t q_sat = get_at(cc.constants, Cons_idx::Epsilon)*(p_sat/(p_prime - p_sat));
        float_t q_sat_2 = p_sat/(get_at(cc.constants, Cons_idx::R_v)*T_prime);
        float_t p_sat_cosmo = saturation_pressure_water(
            T_prime, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));
        float_t p_sat_ice_old = saturation_pressure_ice(
            T_prime, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_b));
        float_t q_sat_cosmo = get_at(cc.constants, Cons_idx::Epsilon)*(p_sat_cosmo/(p_prime - p_sat_cosmo));
        float_t q_sat_2_cosmo = p_sat_cosmo/(get_at(cc.constants, Cons_idx::R_v)*T_prime);
        float_t S_i_2 = p_sat/(p_prime-p_sat) * (p_prime - p_sat_ice)/p_sat_ice;
        std::cout << "------------------- Substep with dt " << dt << "\n"
                  << "traj: " << cc.traj_id << " \np_prime: " << p_prime.getValue()
                  << "\np_sat: " << p_sat.getValue()
                  << "\np_sat_ice: " << p_sat_ice.getValue()
                  << "\nS_i: " << S_i.getValue()
                  << "\nS_i_old: "
                  << ((T_prime < get_at(cc.constants, Cons_idx::T_freeze)) ? e_d/p_sat_ice_old : float_t(1))
                  << "\nS_i2: " << S_i_2.getValue()
                  << "\ns_sw: " << s_sw.getValue()
                  << "\ns_si: " << s_si.getValue()
                  << "\nS: " << S.getValue()
                  << "\nT: " << T_prime.getValue()
                  << "\nqv: " << qv_prime.getValue()
                  << "\nqc: " << qc_prime.getValue()
                  << "\nqr: " << qr_prime.getValue()
                  << "\nqi: " << qi_prime.getValue()
                  << "\nqg: " << qg_prime.getValue()
                  << "\nqh: " << qh_prime.getValue()
                  << "\nqs: " << qs_prime.getValue()
                  << "\nNc: " << Nc.getValue()
                  << "\nNr: " << Nr.getValue()
                  << "\nmeanNr: " << qr_prime/x
                  << "\nmaxNr: " << qr_prime/get_at(cc.rain.constants, Particle_cons_idx::min_x_depo)
                  << "\nminNr: " << qr_prime/get_at(cc.rain.constants, Particle_cons_idx::max_x)
                  << "\nNi: " << Ni.getValue()
                  << "\nNg: " << Ng.getValue()
                  << "\nNs: " << Ns.getValue()
                  << "\nNh: " << Nh.getValue()
                  << "\n\n";
    }
#endif

    float_t rime_rate_qc, rime_rate_qr, rime_rate_qi, rime_rate_qs;
    float_t rime_rate_nc, rime_rate_nr;
    float_t rho_inter = log(compute_rhoh(p_prime, T_prime, S,
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b),
        get_at(cc.constants, Cons_idx::R_a),
        get_at(cc.constants, Cons_idx::R_v)) / get_at(cc.constants, Cons_idx::rho_0));
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
        exp(-get_at(cc.constants, Cons_idx::rho_vel_c) * rho_inter);
    cc.rain.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
        exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
        exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    cc.hail.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
        exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    cc.ice.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
        exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    cc.snow.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
        exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
#if defined(TRACE_SAT) || defined(TRACE_QI) || defined(TRACE_QS) || defined(TRACE_QV) || \
    defined(TRACE_QR) || defined(TRACE_QC) || defined(TRACE_QG) || defined(TRACE_QH)
    if (trace)
        std::cout << "traj: " << cc.traj_id << " \ncloud.constants[Particle_cons_idx::rho_v]:    "
                  << get_at(cc.cloud.constants, Particle_cons_idx::rho_v)
                  << "\nrain.constants[Particle_cons_idx::rho_v]:     "
                  << get_at(cc.rain.constants, Particle_cons_idx::rho_v)
                  << "\ngraupel.constants[Particle_cons_idx::rho_v]:  "
                  << get_at(cc.graupel.constants, Particle_cons_idx::rho_v)
                  << "\nhail.constants[Particle_cons_idx::rho_v]:     "
                  << get_at(cc.hail.constants, Particle_cons_idx::rho_v)
                  << "\nice.constants[Particle_cons_idx::rho_v]:      "
                  << get_at(cc.ice.constants, Particle_cons_idx::rho_v)
                  << "\nsnow.constants[Particle_cons_idx::rho_v]:     "
                  << get_at(cc.snow.constants, Particle_cons_idx::rho_v)
                 << "\n";
#endif
    const double EPSILON = 1.0e-20;

#if defined(CCN_AKM)
    ccn_act_hande_akm(p_prime, w_prime, T_prime, qv_prime, qc_prime, Nc,
        EPSILON, res, cc, dt);
#else
    ccn_act_hande(p_prime, w_prime, T_prime, qv_prime, qc_prime, Nc,
        EPSILON, res, cc, dt);
#endif

    ice_activation_phillips(qc_prime, qv_prime, T_prime,
         S_i, n_inact, res, cc, dt);

    // (optional) homogeneous nucleation using KHL06
    ice_nuc_hom(T_prime, w_prime, p_prime, qv_prime,
        qi_prime, Ni, S_i, p_sat_ice, res, cc, dt);

    cloud_freeze_hom(qc_prime, Nc, T_prime, T_c, res, cc, dt);

    vapor_dep_relaxation(qv_prime, qi_prime, Ni,
        qs_prime, Ns, qg_prime, Ng,
        qh_prime, Nh, s_si, p_sat_ice, p_prime,
        T_prime, EPSILON, dep_rate_ice, dep_rate_snow, D_vtp, res, cc, dt);

    ////////////// ice-ice collisions
    ice_self_collection(qi_prime, Ni, T_c, res, cc, dt);

    snow_self_collection(qs_prime, Ns, T_prime, res, cc);

    particle_particle_collection(qi_prime, Ni, qs_prime, Ns,
        qg_prime, Ng, T_prime, T_c, res, cc, dt);

    graupel_hail_conv(qc_prime, qr_prime, qi_prime, qg_prime, Ng,
        p_prime, T_prime, T_c, res, cc, dt);

    hail_collision(qh_prime, Nh, qs_prime, Ns, qi_prime, Ni, T_c, res, cc, dt);

    ////////////// Riming of ice with cloud and rain droplets and conversion to graupel
    riming_cloud_core(qc_prime, Nc, qi_prime, Ni,
        cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc, dt);
    riming_rain_core(qr_prime, Nr, qi_prime, Ni,
        cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(qc_prime, Nc, qr_prime, Nr, qi_prime, Ni,
        dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
        rime_rate_qi, T_prime, dt, res, cc);

    // snow riming
    riming_cloud_core(qc_prime, Nc, qs_prime, Ns,
        cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc, dt);
    riming_rain_core(qr_prime, Nr, qs_prime, Ns,
        cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(qc_prime, Nc, qr_prime, Nr, qs_prime, Ns,
        dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
        rime_rate_qs, T_prime, dt, res, cc);

#ifdef TRACE_QH
    auto qh_before = res[qh_idx];
    auto Nh_before = res[Nh_idx];
#endif
    //// hail cloud riming
    particle_cloud_riming(qc_prime, Nc, T_prime, qh_prime, Nh,
        res[qh_idx], res[Nh_idx], cc.coeffs_hcr, cc.hail, res, cc, dt);

#ifdef TRACE_QH
    if (trace)
        std::cout << "traj: " << cc.traj_id
            << " hail cloud riming dqh " << res[qh_idx]-qh_before << ", dNh " << res[Nh_idx]-Nh_before << "\n";
    qh_before = res[qh_idx];
    Nh_before = res[Nh_idx];
#endif
    //// hail rain riming
    particle_rain_riming(qr_prime, Nr, T_prime, qh_prime, Nh,
        res[qh_idx], res[Nh_idx], cc.coeffs_hrr, cc.hail, res, cc, dt);
#ifdef TRACE_QH
    if (trace)
        std::cout << "traj: " << cc.traj_id
            << " hail rain riming dqh " << res[qh_idx]-qh_before << ", dNh " << res[Nh_idx]-Nh_before << "\n";
#endif
#ifdef TRACE_QG
    auto qg_before = res[qg_idx];
    auto Ng_before = res[Ng_idx];
#endif
    //// graupel cloud riming
    particle_cloud_riming(qc_prime, Nc, T_prime, qg_prime, Ng,
        res[qg_idx], res[Ng_idx], cc.coeffs_gcr, cc.graupel, res, cc, dt);

#ifdef TRACE_QG
    if (trace)
        std::cout << "traj: " << cc.traj_id
            << " graupel cloud riming dqg " << res[qg_idx]-qg_before << ", dNg " << res[Ng_idx]-Ng_before << "\n";
    qg_before = res[qg_idx];
    Ng_before = res[Ng_idx];
#endif
    //// graupel rain riming
    particle_rain_riming(qr_prime, Nr, T_prime, qg_prime, Ng,
        res[qg_idx], res[Ng_idx], cc.coeffs_grr, cc.graupel, res, cc, dt);
#ifdef TRACE_QG
    if (trace)
        std::cout << "traj: " << cc.traj_id
            << " graupel rain riming dqg " << res[qg_idx]-qg_before << ", dNg " << res[Ng_idx]-Ng_before << "\n";
#endif

    rain_freeze(qr_prime, Nr, T_prime, dt, res, cc);

    ice_melting(qi_prime, Ni, T_prime, res, cc, dt);

    snow_melting(qs_prime, Ns, T_prime, res, cc, dt);

    graupel_melting(qg_prime, Ng, T_prime, res, cc, dt);

    hail_melting(qh_prime, Nh, T_prime, res, cc, dt);

#ifdef TRACE_QS
    auto qs_before = res[qs_idx];
#endif
    evaporation(p_sat, s_sw, T_prime,
        qs_prime, Ns, res[qs_idx], cc.snow, res, cc, dt);

#ifdef TRACE_QS
    if (trace)
        std::cout << "traj: " << cc.traj_id << " evaporation dqs " << res[qs_idx]-qs_before << "\n";
#endif

#ifdef TRACE_QG
    qg_before = res[qg_idx];
#endif
    evaporation(p_sat, s_sw, T_prime,
        qg_prime, Ng, res[qg_idx], cc.graupel, res, cc, dt);
#ifdef TRACE_QG
    if (trace)
        std::cout << "traj: " << cc.traj_id << " evaporation dqg " << res[qg_idx]-qg_before << "\n";
#endif

#ifdef TRACE_QI
    auto qi_before = res[qi_idx];
#endif
    evaporation(p_sat, s_sw, T_prime,
        qi_prime, Ni, res[qi_idx], cc.ice, res, cc, dt);
#ifdef TRACE_QI
    if (trace)
        std::cout << "traj: " << cc.traj_id << " evaporation dqi " << res[qi_idx]-qi_before << "\n";
#endif

    if (auto_type == 1) {
        auto_conversion_kb(qc_prime, Nc, qr_prime, res, cc, dt);

    } else if (auto_type == 2) {
        // Not implemented since it appears to be not very interesting
    } else {
        auto_conversion_sb(qc_prime, Nc, qr_prime, res, cc, dt);
    }

    rain_self_collection_sb(qr_prime, Nr, res, cc, dt);

    rain_evaporation_sb(qr_prime, Nr, qc_prime,
        T_prime, p_prime, s_sw, p_sat, res, cc, dt);

    sedimentation_explicit(
        qc_prime, qr_prime, Nr,
        qs_prime, Ns, qi_prime, Ni,
        qh_prime, Nh, qg_prime, Ng,
        res, cc, dt);

    precipitation_efficiency(res);
#if defined TRACE_QC && defined IN_SAT_ADJ
    if (trace)
        std::cout << "traj: " << cc.traj_id
            << " before saturation adj dqc " << res[qc_idx] << ", dNc " << res[Nc_idx] << "\n";
#endif
#if defined TRACE_QV && defined IN_SAT_ADJ
    if (trace)
        std::cout << "traj: " << cc.traj_id << " before saturation adj dqv " << res[qv_idx] << "\n";
#endif

#ifdef IN_SAT_ADJ
    saturation_adjust(T_prime, p_prime,
                      qv_prime, qc_prime, res, cc);
#endif

#ifdef TRACE_QV
    if (trace)
        std::cout << "traj: " << cc.traj_id << " end dqv " << res[qv_idx] << "\n";
#endif

#ifdef TRACE_QI
    if (trace)
        std::cout << "traj: " << cc.traj_id << " end dqi " << res[qi_idx] << ", dNi " << res[Ni_idx] << "\n";
#endif
#ifdef TRACE_QS
    if (trace)
        std::cout << "traj: " << cc.traj_id << " end dqs " << res[qs_idx] << ", dNs " << res[Ns_idx] << "\n";
#endif
#ifdef TRACE_QR
    if (trace)
        std::cout << "traj: " << cc.traj_id << " end dqr " << res[qr_idx] << ", dNr " << res[Nr_idx] << "\n";
#endif
#ifdef TRACE_QC
    if (trace)
        std::cout << "traj: " << cc.traj_id << " end dqc " << res[qc_idx] << ", dNc " << res[Nc_idx] << "\n";
#endif
#ifdef TRACE_QG
    if (trace)
        std::cout << "traj: " << cc.traj_id << " end dqg " << res[qg_idx] << ", dNg " << res[Ng_idx] << "\n";
#endif
#ifdef TRACE_QH
    if (trace)
        std::cout << "traj: " << cc.traj_id << " end dqh " << res[qh_idx] << ", dNh " << res[Nh_idx] << "\n";
#endif

    ////// Get back to non-prime
    res[qc_idx] /= ref.qref;
    res[qr_idx] /= ref.qref;
    res[qv_idx] /= ref.qref;
    res[qi_idx] /= ref.qref;
    res[qs_idx] /= ref.qref;
    res[qg_idx] /= ref.qref;
    res[qh_idx] /= ref.qref;
    res[qr_out_idx] /= ref.qref;
    res[qi_out_idx] /= ref.qref;
    res[qs_out_idx] /= ref.qref;
    res[qg_out_idx] /= ref.qref;
    res[qh_out_idx] /= ref.qref;
    res[T_idx] /= ref.Tref;

    if (fixed) {
        res[p_idx] = 0;
        res[T_idx] = 0;
        res[w_idx] = 0;
        res[z_idx] = 0;
    } else {
        // Compute nondimensional coefficients
        float_t C1 = (ref.tref*get_at(cc.constants, Cons_idx::gravity_acc)*ref.wref)
            / (get_at(cc.constants, Cons_idx::R_a)*ref.Tref);
        float_t C2 = ((1.0-get_at(cc.constants, Cons_idx::Epsilon))*ref.qref)
            / get_at(cc.constants, Cons_idx::Epsilon);
        // Calculate pressure and temperature change in non-prime directly
        // Pressure change from ascent (C1), change in partial pressure from water vapor.
        res[p_idx] = -(C1/(1.0 + C2*(qv/(1.0 + qv_prime))))*((p*w)/T);
        res[T_idx] += (res[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
            + res[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp)
            - get_at(cc.constants, Cons_idx::gravity_acc) * w_prime/get_at(cc.constants, Cons_idx::cp)) / ref.Tref;
        res[w_idx] += cc.get_dw();
//        res[w_idx] += get_at(cc.constants, Cons_idx::dw)/ref.wref;
        res[z_idx] += w_prime/ref.zref;
    }

#ifdef TRACE_ENV
    if (trace) {
        std::cout << "traj: " << cc.traj_id << " \nPressure " << res[p_idx] << "\n"
                  << "Temperat " << res[T_idx] << "\n"
                  << "ascent   " << res[w_idx] << "\n"
                  << "height   " << res[z_idx] << "\n";
        std::cout << "traj: " << cc.traj_id << " \nlat_heat " << res[lat_heat_idx]
                  << "\nlat_cool " << res[lat_cool_idx]
                  << "\nlat_heat_air " << res[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                  << "\nlat_cool " << res[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp)
                  << "\ndT ascent " << - get_at(cc.constants, Cons_idx::gravity_acc) * w_prime
                    / get_at(cc.constants, Cons_idx::cp) << "\n";
        std::cout << "traj: " << cc.traj_id << " w_prime " << w_prime << "\n"
                  << "z_prime " << z_prime << "\n"
                  << "zref " << ref.zref << "\n";
    }
#endif
}
/** @} */  // end of group ufunc
