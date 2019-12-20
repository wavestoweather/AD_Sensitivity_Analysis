#ifndef RK4_H
#define RK4_H

#include "codi.hpp"

#include "constants.h"
#include "user_functions.h"
#include <vector>
#include <math.h>
#include "types.h"

using v_rev = std::vector<codi::RealReverse>;
////////////////////////////////////////////////////////////////////////////////
// =============================================================================
// This file provides the function for a single step with the
// Runge-Kutta 4 method
// =============================================================================
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// This method computes a single step using the Runge-Kutta 4 method
// for the ODE
// y' = RHS(y)
// with old system state yold and new system state ynew.
//
void RK4_step(
    v_rev &ynew,
    v_rev &yold,
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
    RHS(k, yold, ref, cc, nc); // k = k1

    for(int ii = 0 ; ii < num_comp ; ii++){
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii]; // Add k1-part to the result

    }
    if(fixed)
        for(int ii=0; ii<update_idx; ++ii)
            ytmp[ii] = yold[ii];

    //
    // Do all computations involving k2
    //
    RHS(k, ytmp, ref, cc, nc); // k = k2

    for(int ii = 0 ; ii < num_comp ; ii++){
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii]; // Add k2-part to the result
    }
    if(fixed)
        for(int ii=0; ii<update_idx; ++ii)
            ytmp[ii] = yold[ii];

    //
    // Do all computations involving k3
    //
    RHS(k, ytmp, ref, cc, nc); // k = k3

    for(int ii = 0 ; ii < num_comp ; ii++){
        ytmp[ii] = yold[ii] + cc.dt*k[ii]; // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii]; // Add k3-part to the result
    }
    if(fixed)
        for(int ii=0; ii<update_idx; ++ii)
            ytmp[ii] = yold[ii];

    //
    // Do all computations involving k4
    //
    RHS(k, ytmp, ref, cc, nc); // k = k4

    for(int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] += cc.dt_sixth*k[ii];

    if(fixed)
        for(int ii=0; ii < update_idx; ii++)
            ynew[ii] = yold[ii];
    if(k[qv_idx] != 0.0)
    {
        std::cout << "dqv: " << k[qv_idx] << "\n";
    }
    // std::cout << "\n#####################################################\n";

} // End of method RK4_step
////////////////////////////////////////////////////////////////////////////////

// 2 moment bulk scheme 10.1175/JAS3980 without ice
void RK4_step_2_no_ice(
    v_rev &ynew,
    v_rev &yold,
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

// 2 moment bulk scheme 10.1175/JAS3980 without ice
void RK4_step_2_sb_ice(
    v_rev &ynew,
    v_rev &yold,
    const reference_quantities_t& ref,
    model_constants_t& cc,
    nc_parameters_t& nc,
    bool fixed,
    const int counter)
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
    RHS_SB(k, yold, ref, cc, nc, cc.dt_sixth, fixed, counter); // k = k1

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii]; // Add k1-part to the result

    }

    //
    // Do all computations involving k2
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_third, fixed, counter); // k = k2

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii]; // Add k2-part to the result
    }

    //
    // Do all computations involving k3
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_third, fixed, counter); // k = k3

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        ytmp[ii] = yold[ii] + cc.dt*k[ii]; // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii]; // Add k3-part to the result
    }

    //
    // Do all computations involving k4
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_sixth, fixed, counter); // k = k4

    for(int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] += cc.dt_sixth*k[ii];


}
#endif
