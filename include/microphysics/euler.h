#pragma once

#include <math.h>

#include <algorithm>
#include <vector>

#include "codi.hpp"

#include "include/microphysics/constants.h"
#include "include/microphysics/user_functions.h"
#include "include/types/model_constants_t.h"
#include "include/types/reference_quantities_t.h"
#include "include/microphysics/euler.h"

// This method computes a single Euler-step with timestep dt for the ODE
// y' = RHS(y)
// taking yold as the initial value and storing the result in ynew.
template<class float_t>
void euler_step(
    std::vector<float_t> &ynew,
    std::vector<float_t> &yold,
    const reference_quantities_t& ref,
    model_constants_t<float_t>& cc,
    bool fixed) {

    std::vector<float_t> k1(num_comp);
    // Evaluate the RHS
    RHS_SB(k1, yold, ref, cc, cc.dt, fixed);

    // Do the step
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] = yold[ii] + cc.dt*k1[ii];
    }

    set_limits(ynew, ref, cc);
    sediment_q_total += cc.dt*sediment_q;
    sediment_q = 0;
    sediment_n_total += cc.dt*sediment_n;
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
}
