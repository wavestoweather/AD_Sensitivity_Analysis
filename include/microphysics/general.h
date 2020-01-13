#ifndef GENERAL_H
#define GENERAL_H

#include <stdio.h>
#include "codi.hpp"
#include "constants.h"


/** @defgroup general General Purpose Functions
 * This file contains several general purpose functions.
 * @{
 */

/**
 * This function prints the contents of a vector y with length
 * num_comp to stdout.
 *
 * @param y Contents to print of size num_comp
 */
void printy(const double y[])
{
  std::cout << "\n";
  for(int ii = 0 ; ii < num_comp ; ii++){
    std::cout << "y[" << std::to_string(ii) << "] = " << y[ii] << "\n";
  }
  std::cout << std::endl;
}

/**
 * This function computes
 * \f[ \vec{y} = \vec{v1} + a*\vec{v2} \f]
 * for vectors \f$\vec{y}, \vec{v1}, \vec{v2}\f$ and a scalar \f$a\f$.
 *
 * @param y On out: the result of \f4\vec{y} = \vec{v1} + a*\vec{v2}\f$
 * @param v1 \f$\vec{v1}\f$
 * @param v2 \f$\vec{v2}\f$
 * @param a \f$a\f$
 */
void v1pav2(double y[],
	    const double v1[],
	    const double v2[],
	    const double a)
{

  for(int ii = 0 ; ii < num_comp ; ii++){
    y[ii] = v1[ii] + a*v2[ii];
  }

}

/**
 * Print reference quantities temperature, pressure,
 * mixing-ratio, vertical velocity, time.
 *
 * @param ref Struct with reference quantities.
 */
void print_reference_quantities(reference_quantities_t &ref)
{
  std::cout << "\nReference quantities\n"
	    << "--------------------\n"
        << "Temperature: " << ref.Tref << " Kelvin\n"
        << "Pressure: " << ref.pref << " Pascal\n"
        << "Mixing-ratio: " << ref.qref << "\n"
        << "Vertical velocity: " << ref.wref << " meter per second\n"
        << "Time: " << ref.tref << " Second\n"
        << std::endl;
}

/**
 * Print constants given a model constants structure, namely integration
 * time, number of steps, scaling factors.
 *
 * @param cc A structure with model constants.
 */
void print_constants(model_constants_t &cc)
{

  std::cout << "\nModel constants:\n"
	    << "----------------\n"
	    << "\n"
	    << "Final integration time: " << cc.t_end_prime << " seconds\n"
        << "Nondimensional final integration time: " << cc.t_end << "\n"
	    << "Timestep: " << cc.dt_prime << " seconds\n"
	    << "Snapshot Index: " << cc.snapshot_index << "\n"
	    << "Nondimensional timestep: " << cc.dt << "\n"
	    << "Number of iterations: " << cc.num_steps << "\n"
        << "Number of substeps: " << cc.num_sub_steps << "\n"
	    << "a1_scale: " << cc.a1_scale << "\n"
        << "a2_scale: " << cc.a2_scale << "\n"
        << "e1_scale: " << cc.e1_scale << "\n"
        << "e2_scale: " << cc.e2_scale << "\n"
        << "d_scale: " << cc.d_scale << "\n"
	    << "Scaling factor: " << cc.scaling_fact << "\n"
	    << std::endl;
}
/** @} */ // end of group general

#endif
