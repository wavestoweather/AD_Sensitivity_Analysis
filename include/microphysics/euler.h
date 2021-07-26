#pragma once

#include "codi.hpp"

#include "include/microphysics/constants.h"
#include "include/microphysics/user_functions.h"

////////////////////////////////////////////////////////////////////////////////
// =============================================================================
// This file provides the function for a single step with the
// explicit Euler method
// =============================================================================
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// This method computes a single Euler-step with timestep dt for the ODE
// y' = RHS(y)
// taking yold as the initial value and storing the result in ynew.
void euler_step(
    codi::RealReverse ynew[],
    codi::RealReverse yold[],
    const int num_comp,
    const reference_quantities_t &ref,
    model_constants_t &cc) {

    codi::RealReverse k1[num_comp];
    // Evaluate the RHS
    RHS(k1, yold, num_comp, ref, cc);

    // Do the step
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] = yold[ii] + cc.dt*k1[ii];
    }
}
