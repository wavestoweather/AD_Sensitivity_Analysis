#pragma once

#include <vector>
#include <codi.hpp>

#include "include/microphysics/constants.h"

// A 4D lookup table
/**
 * A 4D lookup table on a grid used in ICON.
 */
template<class float_t>
struct table_t{
    /**
     * Number of grid points in direction 1.
     */
    uint64_t n1;
    /**
     * Number of grid points in direction 1.
     */
    uint64_t n2;
    /**
     * Number of grid points in direction 1.
     */
    uint64_t n3;
    /**
     * Number of grid points in direction 1.
     */
    uint64_t n4;

    std::vector<float_t> x1; /*!< Grid vector */
    std::vector<float_t> x2; /*!< Grid vector */
    std::vector<float_t> x3; /*!< Grid vector */
    std::vector<float_t> x4; /*!< Grid vector */
    float_t dx1; /*!< Grid distances for vector x1 */
    float_t dx2; /*!< Grid distances for vector x2 */
    float_t dx3; /*!< Grid distances for vector x3 */
    float_t dx4; /*!< Grid distances for vector x4 */
    float_t odx1; /*!< One over dx1 */
    float_t odx2; /*!< One over dx2 */
    float_t odx3; /*!< One over dx3 */
    float_t odx4; /*!< One over dx4 */
    std::vector<float_t> table; /*!< The table values */

    table_t();

    /**
     * Get the value at a given index.
     *
     * @param i Index along vector 1
     * @param j Index along vector 2
     * @param k Index along vector 3
     * @param l Index along vector 4
     */
    float_t get(
        const uint64_t i,
        const uint64_t j,
        const uint64_t k,
        const uint64_t l) const;
};
