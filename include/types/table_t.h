#pragma once

#include <vector>
#include "codi.hpp"

// A 4D lookup table
/**
 * A 4D lookup table on a grid used in ICON.
 */
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

    std::vector<codi::RealReverse> x1; /*!< Grid vector */
    std::vector<codi::RealReverse> x2; /*!< Grid vector */
    std::vector<codi::RealReverse> x3; /*!< Grid vector */
    std::vector<codi::RealReverse> x4; /*!< Grid vector */
    codi::RealReverse dx1; /*!< Grid distances for vector x1 */
    codi::RealReverse dx2; /*!< Grid distances for vector x2 */
    codi::RealReverse dx3; /*!< Grid distances for vector x3 */
    codi::RealReverse dx4; /*!< Grid distances for vector x4 */
    codi::RealReverse odx1; /*!< One over dx1 */
    codi::RealReverse odx2; /*!< One over dx2 */
    codi::RealReverse odx3; /*!< One over dx3 */
    codi::RealReverse odx4; /*!< One over dx4 */
    std::vector<codi::RealReverse> table; /*!< The table values */

    table_t();

    /**
     * Get the value at a given index.
     *
     * @param i Index along vector 1
     * @param j Index along vector 2
     * @param k Index along vector 3
     * @param l Index along vector 4
     */
    codi::RealReverse get(
        uint64_t i,
        uint64_t j,
        uint64_t k,
        uint64_t l) const;
};
