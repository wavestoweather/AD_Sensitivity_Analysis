#ifndef RK4_H
#define RK4_H

#include "codi.hpp"

#include "constants.h"
#include "user_functions.h"
#include <vector>
#include <math.h>
#include "types.h"

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
#if defined(RK4_ONE_MOMENT)
    RHS(k, ytmp, ref, cc, nc); // k = k2
#endif

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
#if defined(RK4_ONE_MOMENT)
    RHS(k, ytmp, ref, cc, nc); // k = k3
#endif

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
#if defined(RK4_ONE_MOMENT)
    RHS(k, ytmp, ref, cc, nc); // k = k4
#endif

    for(int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] += cc.dt_sixth*k[ii];

    if(fixed)
        for(int ii=0; ii < update_idx; ii++)
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
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii]; // Add k1-part to the result

    }

    //
    // Do all computations involving k2
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_third, fixed); // k = k2

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii]; // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii]; // Add k2-part to the result
    }

    //
    // Do all computations involving k3
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_third, fixed); // k = k3

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        ytmp[ii] = yold[ii] + cc.dt*k[ii]; // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii]; // Add k3-part to the result
    }

    //
    // Do all computations involving k4
    //
    RHS_SB(k, ytmp, ref, cc, nc, cc.dt_sixth, fixed); // k = k4

    for(int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] += cc.dt_sixth*k[ii];
}

/** @} */ // end of group rk

#endif
