#pragma once

#include <stdio.h>
#include <sys/stat.h>

#include <string>
#include <vector>

#include <codi.hpp>

#include "include/microphysics/constants.h"
#include "include/types/particle_model_constants_t.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/model_constants_t.h"
#include "include/types/segment_t.h"

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
void printy(const double y[]);

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
        const double a);

/**
 * Print parameters for given particle.
 *
 * @param pc Particle constants
 * @param title Name of particle
 */
template<class float_t>
void print_particle_params(
    const particle_model_constants_t<float_t> &pc,
    const std::string title);

/**
 * Print reference quantities temperature, pressure,
 * mixing-ratio, vertical velocity, time.
 *
 * @param ref Struct with reference quantities.
 */
void print_reference_quantities(const reference_quantities_t &ref);

/**
 * Print segments for ensembles.
 *
 * @param segments Vector of segments.
 */
void print_segments(const std::vector<segment_t> &segments);

/**
 * Check if a file exists.
 */
inline bool exists(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

/** @} */  // end of group general
