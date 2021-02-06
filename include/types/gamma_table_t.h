#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include "codi.hpp"
#include <vector>

/**
 * A lookup table for the lower incomplete gamma function from ICON
 * mo_2mom_mcrph_util.f90 at incgfct_lower_lookupcreate:
 * \f[ \text{int}(0)(x) \exp(-t) t^{a-1} \text{d}t \f]
 *
 * with constant a from x=0 to the 99.5% value of the normalized incomplete
 * gamma function. This 99.5 % - value has been fitted
 * with high accuracy as function of a in the range a in [0;20], but can
 * safely be applied also to higher values of a. (Fit created with the
 * matlab-program "gamma_unvoll_lower_lookup.m" by Ulrich Blahak, 2008/11/13).

 * The last value in the table corresponds to x = infinity, so that
 * during the reconstruction of incgfct-values from the table,
 * the x-value can safely be truncated at the maximum table x-value.
 */
struct gamma_table_t{
    uint64_t n_bins;    /*!< Number of bins of the lookup table. */
    /**
     * Number of bins of the lookup table with higher resolution if desired.
     */
    uint64_t n_bins_highres;
    double a; /*!< Value at which the function is being fit. */
    std::vector<double> x; /*!< always starts at 0 and has distant dx. */
    std::vector<double> x_highres; /*!< as x but with higher resolution. */
    double dx; /*!< Distance of each value in x. */
    double dx_highres; /*!< Distance of each value in x_highres. */
    double odx; /*!< 1/dx. */
    double odx_highres; /*!< 1/dx_highres. */
    std::vector<double> igf; /*!< Values of the inc. gamma function at (a,x). */
    /**
     * Values of the inc. gamma function at (a,x) with higher resolution.
     */
    std::vector<double> igf_highres;

    gamma_table_t();

    /**
     * Get the lower incomplete gamma function value for a given x coordinate.
     *
     * @param x Coordinate at which to evaluate the gamma function.
     */
    template <class A = codi::RealReverse>
    A look_lo(A x) const;

    /**
     * Get the upper incomplete gamma function value for a given x coordinate.
     *
     * @param x Coordinate at which to evaluate the gamma function.
     */
    template <class A = codi::RealReverse>
    A look_up(A x) const;

    /** Init lookup table for the incomplete gamma function.
     * From ICON mo_2mom_mcrph_util.f90 incgfct_lower_lookupcreate
     * Create a lookup-table for the lower incomplete gamma function
     * \f[ \text{int}(0)(x) \exp(-t) t^{a-1} \text{d}t \f]
     * with constant a from x=0 to the 99.5% value of the normalized incomplete
     * gamma function. This 99.5 % - value has been fitted
     * with high accuracy as function of a in the range a in [0;20], but can
     * safely be applied also to higher values of a. (Fit created with the
     * matlab-program "gamma_unvoll_lower_lookup.m" by Ulrich Blahak, 2008/11/13).
     *
     * The last value in the table corresponds to x = infinity, so that
     * during the reconstruction of incgfct-values from the table,
     * the x-value can safely be truncated at the maximum table x-value.
     *
     * @param nl Number of bins in the table.
     * @param nl_highres Optional number of bins for the high resolution table.
     * @param a Value where the incomplete gamma function shall be interpolated around.
     */
    void init_gamma_table(
        const uint64_t &nl,
        const uint64_t &nl_highres,
        const double &a);
};