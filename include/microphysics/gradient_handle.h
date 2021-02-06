#ifndef GRADIENT_HANDLE_H
#define GRADIENT_HANDLE_H

// #include <stdlib.h>
// #include <cmath>
// #include <string>
// #include <netcdf>
// #include <vector>
// #include <iostream>
// #include <fstream>
// #include <iterator>
// #include "codi.hpp"

// #include "include/microphysics/constants.h"
// #include "include/types/model_constants_t.h"


// /** @defgroup gradients Handle CODIPACK Gradients
//  * Functions to get gradients and handle them.
//  * @{
//  */


// /**
//  *
//  */
// void get_gradients(
//     std::vector<codi::RealReverse> &y_single_new,
//     std::vector< std::array<double, num_par > > &y_diff,
//     model_constants_t &cc,
//     codi::RealReverse::TapeType &tape)
// {
//     for(uint32_t ii = 0 ; ii < num_comp ; ii++)
//         tape.registerOutput(y_single_new[ii]);

//     tape.setPassive();
//     for(uint32_t ii = 0 ; ii < num_comp ; ii++)
//     {
//         y_single_new[ii].setGradient(1.0);
//         tape.evaluate();

//         cc.get_gradient(y_diff[ii]);
//         tape.clearAdjoints();
//     }
//     tape.reset();
// }

/** @} */ // end of group io

#endif
