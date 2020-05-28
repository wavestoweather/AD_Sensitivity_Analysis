#ifndef GRADIENT_HANDLE_H
#define GRADIENT_HANDLE_H

#include <stdlib.h>
#include <cmath>
#include <string>
#include <netcdf>
#include <vector>
#include "constants.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include "codi.hpp"
#include "types.h"


/** @defgroup gradients Handle CODIPACK Gradients
 * Functions to get gradients and handle them.
 * @{
 */


/**
 * Register all model parameters on tape for gradient evaluation.
 */
void register_everything(
    codi::RealReverse::TapeType &tape,
    model_constants_t &cc)
{
#if defined(RK4_ONE_MOMENT)
    // Dimensional coefficients
    tape.registerInput(cc.a1_prime);    // Autoconversion
    tape.registerInput(cc.a2_prime);    // Accretion
    tape.registerInput(cc.e1_prime);    // Evaporation
    tape.registerInput(cc.e2_prime);    // Evaporation
    tape.registerInput(cc.d_prime);     // Sedimentation
    tape.registerInput(cc.Nc_prime);    // Concentration of cloud droplets

    // Exponents
    tape.registerInput(cc.gamma);       // Autoconversion
    tape.registerInput(cc.betac);       // Accretion
    tape.registerInput(cc.betar);       // Accretion
    tape.registerInput(cc.delta1);      // Evaporation
    tape.registerInput(cc.delta2);      // Evaporation
    tape.registerInput(cc.zeta);        // Sedimentation

#elif defined(RK4ICE) || defined(RK4NOICE)
    // Dimensional coefficients
    tape.registerInput(cc.a1_prime);    // Autoconversion
    tape.registerInput(cc.a2_prime);    // Accretion
    tape.registerInput(cc.e1_prime);    // Evaporation
    tape.registerInput(cc.e2_prime);    // Evaporation
    tape.registerInput(cc.d_prime);     // Sedimentation
    tape.registerInput(cc.Nc_prime);    // Concentration of cloud droplets

    // Exponents
    tape.registerInput(cc.gamma);       // Autoconversion
    tape.registerInput(cc.betac);       // Accretion
    tape.registerInput(cc.betar);       // Accretion
    tape.registerInput(cc.delta1);      // Evaporation
    tape.registerInput(cc.delta2);      // Evaporation
    tape.registerInput(cc.zeta);        // Sedimentation

    // ICON parameters
    tape.registerInput(cc.rain_gfak);
    tape.registerInput(cc.cloud_k_au);
    tape.registerInput(cc.cloud_k_sc);
    tape.registerInput(cc.kc_autocon);
    tape.registerInput(cc.inv_z);
    cc.rain.register_input(tape);
    cc.cloud.register_input(tape);
    cc.hail.register_input(tape);
    cc.ice.register_input(tape);
    cc.snow.register_input(tape);
    cc.graupel.register_input(tape);
#endif
}

/**
 *
 */
void get_gradients(
    std::vector<codi::RealReverse> &y_single_new,
    std::vector< std::array<double, num_par > > &y_diff,
    model_constants_t &cc,
    codi::RealReverse::TapeType &tape)
{
    for(int ii = 0 ; ii < num_comp ; ii++)
        tape.registerOutput(y_single_new[ii]);

    tape.setPassive();
    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        y_single_new[ii].setGradient(1.0);
        tape.evaluate();

#if defined(RK4_ONE_MOMENT)
        y_diff[ii][0] = cc.a1_prime.getGradient();
        y_diff[ii][1] = cc.a2_prime.getGradient();
        y_diff[ii][2] = cc.e1_prime.getGradient();
        y_diff[ii][3] = cc.e2_prime.getGradient();
        y_diff[ii][4] = cc.d_prime.getGradient();

        y_diff[ii][5] = cc.gamma.getGradient();
        y_diff[ii][6] = cc.betac.getGradient();
        y_diff[ii][7] = cc.betar.getGradient();
        y_diff[ii][8] = cc.delta1.getGradient();
        y_diff[ii][9] = cc.delta2.getGradient();
        y_diff[ii][10] = cc.zeta.getGradient();

        y_diff[ii][11] = cc.Nc_prime.getGradient();

#elif defined(RK4ICE) || defined(RK4NOICE)
        y_diff[ii][0] = cc.a1_prime.getGradient();
        y_diff[ii][1] = cc.a2_prime.getGradient();
        y_diff[ii][2] = cc.e1_prime.getGradient();
        y_diff[ii][3] = cc.e2_prime.getGradient();
        y_diff[ii][4] = cc.d_prime.getGradient();

        y_diff[ii][5] = cc.gamma.getGradient();
        y_diff[ii][6] = cc.betac.getGradient();
        y_diff[ii][7] = cc.betar.getGradient();
        y_diff[ii][8] = cc.delta1.getGradient();
        y_diff[ii][9] = cc.delta2.getGradient();
        y_diff[ii][10] = cc.zeta.getGradient();

        y_diff[ii][11] = cc.Nc_prime.getGradient();

        y_diff[ii][12] = cc.rain_gfak.getGradient();
        y_diff[ii][13] = cc.cloud_k_au.getGradient();
        y_diff[ii][14] = cc.cloud_k_sc.getGradient();
        y_diff[ii][15] = cc.kc_autocon.getGradient();
        y_diff[ii][16] = cc.inv_z.getGradient();

        uint64_t idx = 17;
        cc.rain.get_gradient(y_diff[ii], idx);
        cc.cloud.get_gradient(y_diff[ii], idx);
        cc.graupel.get_gradient(y_diff[ii], idx);
        cc.hail.get_gradient(y_diff[ii], idx);
        cc.ice.get_gradient(y_diff[ii], idx);
        cc.snow.get_gradient(y_diff[ii], idx);
#endif

        tape.clearAdjoints();
    }
    tape.reset();
}

/** @} */ // end of group io

#endif
