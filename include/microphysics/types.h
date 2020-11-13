#pragma once
#include "codi.hpp"
#include <netcdf>
#include <random>
#include <unordered_map>
#include <boost/math/special_functions/gamma.hpp>
#include "include/misc/error.h"

using namespace netCDF;


/** @defgroup types Commonly Used Types
 * Various types for storing model constants, reading data from netCDF files
 * and lookup tables used in Icon.
 * @{
 */

/**
 *  Struct to hold all reference values.
 */
struct reference_quantities_t{
  double Tref;          /*!< Reference temperature */
  double pref;          /*!< Reference pressure */
  double qref;          /*!< Reference mixing-ratio */
  double wref;          /*!< Reference vertical velocity */
  double tref;          /*!< Reference time */
  double zref;          /*!< Reference distances */

  double Nref;          /*!< DUMMY */
};

/**
 * Struct to hold all model constants regarding particles that we want to
 * know the gradients using AD of.
 * Has getter and register functions for codi::RealReverse for its
 * members.
 */
struct particle_model_constants_t{

    /**
     * Geometry coefficients.
     */
    codi::RealReverse a_geo; /*!< Coefficient for diameter size calculation */
    codi::RealReverse b_geo; /*!< Exponent for diameter size calculation */

    /**
     * Minimum size of particle for mean meass calculation.
     */
    codi::RealReverse min_x;
    /**
     * Minimum size of particle for CCN activation (cloud) and
     * ice activation (ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_act;
    /**
     * Minimum size of particle for homogenous nucleation (ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_nuc_homo;
    /**
     * Minimum size of particle for heterogeneous nucleation (ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_nuc_hetero;
    /**
     * Minimum size of particle for melting (snow, graupel, ice, hail).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_melt;
    /**
     * Minimum size of particle for evaporation (rain, snow, graupel, ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_evap;
    /**
     * Minimum size of particle for freezing (rain, cloud).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_freezing;
    /**
     * Minimum size of particle for vapor deposition (ice, snow, graupel, hail).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_depo;
    /**
     * Minimum size of particle for ice-ice collision.
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_collision;
    /**
     * Minimum size of particle for different collision processes
     * (snow, rain, ice, snow, graupel).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_collection;
    /**
     * Minimum size of particle for conversion processes (cloud, graupel, ice).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_conversion;
    /**
     * Minimum size of particle for sedimentation (rain, ice, snow, graupel, hail).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_sedimentation;
    /**
     * Minimum size of particle for riming (cloud, rain, ice, snow, hail, graupel).
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_riming;

    codi::RealReverse max_x; /*!< Maximum size of particle. */
    codi::RealReverse sc_theta_q; /*!< Coefficient for ice collision ratio mass. */
    codi::RealReverse sc_delta_q; /*!< Coefficient for ice collision ratio mass. */
    codi::RealReverse sc_theta_n; /*!< Coefficient for collision particle number (ice, snow). */
    codi::RealReverse sc_delta_n; /*!< Coefficient for collision particle number (ice, snow). */
    /**
     * Variance for the assumed Gaussian velocity distributions used in collection and riming processes.
     */
    codi::RealReverse s_vel;
    codi::RealReverse a_vel;    /*!< Coefficient for particle velocity. */
    codi::RealReverse b_vel;    /*!< Exponent for particle velocity. */
    /**
     * Coefficient used in density correction for the increased terminal
     * fall velocity with decreasing air density.
     * \f[ \rho_v = (\rho/\rho_0)^{-\rho_{\text{vel}}} \f]
     */
    codi::RealReverse rho_v;
    codi::RealReverse c_z; /*!< Coefficient for 2nd mass moment. */
    codi::RealReverse sc_coll_n;    /*!< Coefficient in graupel self collection and cloud riming. */
    codi::RealReverse cmu0; /*!< Coefficient for calculating the shape parameter \f$\mu\f$. */
    codi::RealReverse cmu1; /*!< Coefficient for calculating the shape parameter \f$\mu\f$. */
    codi::RealReverse cmu2; /*!< Coefficient for calculating the shape parameter \f$\mu\f$. */
    codi::RealReverse cmu3; /*!< Constant for calculating the shape parameter \f$\mu\f$. */
    codi::RealReverse cmu4; /*!< Constant for calculating the shape parameter \f$\mu\f$. */
    codi::RealReverse cmu5; /*!< Exponent for calculating the shape parameter \f$\mu\f$. */
    codi::RealReverse alpha; /*!< Constant in rain sedimentation. */
    codi::RealReverse beta; /*!< Coefficient for rain sedimentation. */
    codi::RealReverse gamma; /*!< Exponent for rain sedimentation. */
    /**
     * Right edge of incomplete gamma function,
     * which had been initialized with \f[\text{nm}_1\f].
     */
    codi::RealReverse g1;
    /**
     * Right edge of incomplete gamma function,
     * which had been initialized with \f[\text{nm}_2\f].
     */
    codi::RealReverse g2;
    /**
     * Shape parameter of the generalized \f$\Gamma$\f-distribution.
     */
    codi::RealReverse mu;
    /**
     * Shape parameter of the generalized \f$\Gamma$\f-distribution.
     * i.e. used in rain sedimentation as coefficient.
     */
    codi::RealReverse nu;
    /**
     * Used for initializing the incomplete
     * gamma function lookup table 1.
     * Number of bins.
     */
    codi::RealReverse nm1;
    /**
     * Used for initializing the incomplete
     * gamma function lookup table 2.
     * Number of bins.
     */
    codi::RealReverse nm2;
    /**
     * Used for initializing the incomplete
     * gamma function lookup table 3.
     * Number of bins.
     */
    codi::RealReverse nm3;

    codi::RealReverse q_crit_c; /*!<  Riming parameter. */
    codi::RealReverse d_crit_c; /*!<  Riming parameter. */
    codi::RealReverse ecoll_c;  /*!<  Riming coefficient. */
    /**
     * Coefficient for capacity of particle.
     */
    codi::RealReverse cap;
    codi::RealReverse a_ven;    /*!< Vapor deposition coefficient. */
    codi::RealReverse b_ven;    /*!< Currently unused parameter. */

    codi::RealReverse c_s;  /*!< Inverse of capacity. Coefficient in evaporation and vapor deposition. */
    codi::RealReverse a_f; /*!< Constant for average ventilation. Used in melting and ice-vapor processes. */
    codi::RealReverse b_f; /*!< Coefficient for average ventilation. */

    codi::RealReverse alfa_n;       /*!<  Sedimentation velocity coefficient. */
    codi::RealReverse alfa_q;       /*!<  Sedimentation velocity coefficient. */
    codi::RealReverse lambda;       /*!<  Sedimentation velocity coefficient. */
    codi::RealReverse vsedi_min;    /*!<  Minimum sedimentation velocity parameter. */
    codi::RealReverse vsedi_max;    /*!<  Maximum sedimentation velocity parameter. */

    /**
     * Register the model parameters on the tape for codi::RealReverse.
     *
     * @param tape Tape where the parameters should be registered on.
     */
    void register_input(codi::RealReverse::TapeType &tape)
    {
        tape.registerInput(this->a_geo);
        tape.registerInput(this->b_geo);
        tape.registerInput(this->min_x);
        tape.registerInput(this->min_x_act);
        tape.registerInput(this->min_x_nuc_homo);
        tape.registerInput(this->min_x_nuc_hetero);
        tape.registerInput(this->min_x_melt);
        tape.registerInput(this->min_x_evap);
        tape.registerInput(this->min_x_freezing);
        tape.registerInput(this->min_x_depo);
        tape.registerInput(this->min_x_collision);
        tape.registerInput(this->min_x_collection);
        tape.registerInput(this->min_x_conversion);
        tape.registerInput(this->min_x_sedimentation);
        tape.registerInput(this->min_x_riming);
        tape.registerInput(this->max_x);
        tape.registerInput(this->sc_theta_q);
        tape.registerInput(this->sc_delta_q);
        tape.registerInput(this->sc_theta_n);
        tape.registerInput(this->sc_delta_n);
        tape.registerInput(this->s_vel);
        tape.registerInput(this->a_vel);
        tape.registerInput(this->b_vel);
        tape.registerInput(this->rho_v);
        tape.registerInput(this->c_z);
        tape.registerInput(this->sc_coll_n);
        tape.registerInput(this->cmu0);
        tape.registerInput(this->cmu1);
        tape.registerInput(this->cmu2);
        tape.registerInput(this->cmu3);
        tape.registerInput(this->cmu4);
        tape.registerInput(this->cmu5);
        tape.registerInput(this->alpha);
        tape.registerInput(this->beta);
        tape.registerInput(this->gamma);
        tape.registerInput(this->nu);
        tape.registerInput(this->g1);
        tape.registerInput(this->g2);
        tape.registerInput(this->mu);
        tape.registerInput(this->nm1);
        tape.registerInput(this->nm2);
        tape.registerInput(this->nm3);
        tape.registerInput(this->q_crit_c);
        tape.registerInput(this->d_crit_c);
        tape.registerInput(this->ecoll_c);
        tape.registerInput(this->cap);
        tape.registerInput(this->a_ven);
        tape.registerInput(this->b_ven);
        tape.registerInput(this->c_s);
        tape.registerInput(this->a_f);
        tape.registerInput(this->b_f);
        tape.registerInput(this->alfa_n);
        tape.registerInput(this->alfa_q);
        tape.registerInput(this->lambda);
        tape.registerInput(this->vsedi_min);
        tape.registerInput(this->vsedi_max);
    }

    /**
     * Get the gradients of all its members. You need to register them on a
     * type before to get meaningful values.
     *
     * @param out_vec On out: Stores all gradients.
     * @param idx Start index of out_vec where the gradients should be stored.
     */
    template<class T>
    void get_gradient(T &out_vec, uint64_t &idx)
    {
        out_vec[idx] = this->a_geo.getGradient();
        idx++;
        out_vec[idx] = this->b_geo.getGradient();
        idx++;
        out_vec[idx] = this->min_x.getGradient();
        idx++;
        out_vec[idx] = this->min_x_act.getGradient();
        idx++;
        out_vec[idx] = this->min_x_nuc_homo.getGradient();
        idx++;
        out_vec[idx] = this->min_x_nuc_hetero.getGradient();
        idx++;
        out_vec[idx] = this->min_x_melt.getGradient();
        idx++;
        out_vec[idx] = this->min_x_evap.getGradient();
        idx++;
        out_vec[idx] = this->min_x_freezing.getGradient();
        idx++;
        out_vec[idx] = this->min_x_depo.getGradient();
        idx++;
        out_vec[idx] = this->min_x_collision.getGradient();
        idx++;
        out_vec[idx] = this->min_x_collection.getGradient();
        idx++;
        out_vec[idx] = this->min_x_conversion.getGradient();
        idx++;
        out_vec[idx] = this->min_x_sedimentation.getGradient();
        idx++;
        out_vec[idx] = this->min_x_riming.getGradient();
        idx++;
        out_vec[idx] = this->max_x.getGradient();
        idx++;
        out_vec[idx] = this->sc_theta_q.getGradient();
        idx++;
        out_vec[idx] = this->sc_delta_q.getGradient();
        idx++;
        out_vec[idx] = this->sc_theta_n.getGradient();
        idx++;
        out_vec[idx] = this->sc_delta_n.getGradient();
        idx++;
        out_vec[idx] = this->s_vel.getGradient();
        idx++;
        out_vec[idx] = this->a_vel.getGradient();
        idx++;
        out_vec[idx] = this->b_vel.getGradient();
        idx++;
        out_vec[idx] = this->rho_v.getGradient();
        idx++;
        out_vec[idx] = this->c_z.getGradient();
        idx++;
        out_vec[idx] = this->sc_coll_n.getGradient();
        idx++;
        out_vec[idx] = this->cmu0.getGradient();
        idx++;
        out_vec[idx] = this->cmu1.getGradient();
        idx++;
        out_vec[idx] = this->cmu2.getGradient();
        idx++;
        out_vec[idx] = this->cmu3.getGradient();
        idx++;
        out_vec[idx] = this->cmu4.getGradient();
        idx++;
        out_vec[idx] = this->cmu5.getGradient();
        idx++;
        out_vec[idx] = this->alpha.getGradient();
        idx++;
        out_vec[idx] = this->beta.getGradient();
        idx++;
        out_vec[idx] = this->gamma.getGradient();
        idx++;
        out_vec[idx] = this->nu.getGradient();
        idx++;
        out_vec[idx] = this->g1.getGradient();
        idx++;
        out_vec[idx] = this->g2.getGradient();
        idx++;
        out_vec[idx] = this->mu.getGradient();
        idx++;
        out_vec[idx] = this->nm1.getGradient();
        idx++;
        out_vec[idx] = this->nm2.getGradient();
        idx++;
        out_vec[idx] = this->nm3.getGradient();
        idx++;
        out_vec[idx] = this->q_crit_c.getGradient();
        idx++;
        out_vec[idx] = this->d_crit_c.getGradient();
        idx++;
        out_vec[idx] = this->ecoll_c.getGradient();
        idx++;
        out_vec[idx] = this->cap.getGradient();
        idx++;
        out_vec[idx] = this->a_ven.getGradient();
        idx++;
        out_vec[idx] = this->b_ven.getGradient();
        idx++;
        out_vec[idx] = this->c_s.getGradient();
        idx++;
        out_vec[idx] = this->a_f.getGradient();
        idx++;
        out_vec[idx] = this->b_f.getGradient();
        idx++;
        out_vec[idx] = this->alfa_n.getGradient();
        idx++;
        out_vec[idx] = this->alfa_q.getGradient();
        idx++;
        out_vec[idx] = this->lambda.getGradient();
        idx++;
        out_vec[idx] = this->vsedi_min.getGradient();
        idx++;
        out_vec[idx] = this->vsedi_max.getGradient();
        idx++;
    }
};


/**
 * Struct to hold model constants for particle collection that we *usually*
 * are not interested in regarding their influence on the overall model.
 * It is stored within a struct model_constants_t for different collisions
 * such as hail and ice collection or ice and rain riming.
 */
struct collection_model_constants_t{
    codi::RealReverse delta_n_aa;
    codi::RealReverse delta_n_ab;
    codi::RealReverse delta_n_bb;
    codi::RealReverse delta_q_aa;
    codi::RealReverse delta_q_ab;
    codi::RealReverse delta_q_bb;
    codi::RealReverse delta_q_ba;

    codi::RealReverse theta_n_aa;
    codi::RealReverse theta_n_ab;
    codi::RealReverse theta_n_bb;
    codi::RealReverse theta_q_aa;
    codi::RealReverse theta_q_ab;
    codi::RealReverse theta_q_bb;
    codi::RealReverse theta_q_ba;
};

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
        uint64_t l) const
    {
        return table[i*n3*n2*n1 + j*n2*n1 + k*n1 + l];
    }


};


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

    /**
     * Get the lower incomplete gamma function value for a given x coordinate.
     *
     * @param x Coordinate at which to evaluate the gamma function.
     */
    template <class A>
    A look_lo(A x)
    {
        A xt = max(min(x, this->x[this->n_bins-1]), 0.0);
        double tmp = xt.getValue()*this->odx;
        tmp = floor(tmp);
        uint64_t iu = std::min((uint64_t) tmp, this->n_bins-2);
        uint64_t io = iu + 1;
        return this->igf[iu] + (this->igf[io] - this->igf[iu]) * this->odx*(xt-this->x[iu]);
    }

    /**
     * Get the upper incomplete gamma function value for a given x coordinate.
     *
     * @param x Coordinate at which to evaluate the gamma function.
     */
    template <class A>
    A look_up(A x)
    {
        A xt = max(min(x, this->x[this->n_bins-1]), 0.0);
        double tmp = xt.getValue()*this->odx;
        tmp = floor(tmp) ;
        uint64_t iu = std::min((uint64_t) tmp, this->n_bins-2);
        uint64_t io = iu + 1;
        A lookup = this->igf[this->n_bins-1] - this->igf[iu]
            - (this->igf[io] - this->igf[iu]) * this->odx*(xt-this->x[iu]);
        return max(lookup, codi::RealReverse(0.0));
    }

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
        const double &a)
    {
        const double c1 = 36.629433904824623;
        const double c2 = -0.119475603955226;
        const double c3 = 0.339332937820052;
        const double c4 = 1.156369000458310;

        n_bins = nl;
        n_bins_highres = nl_highres;

        x.resize(nl);
        x_highres.resize(nl_highres);
        igf.resize(nl);
        igf_highres.resize(nl_highres);

        // Low resolution
        // maximum x-value (99.5%)
        x[n_bins-2] = c1 * (1.0-exp(c2*pow(a, c3))) + c4*a;
        dx = x[n_bins-2] / (n_bins-2.0);
        odx = 1.0/dx;
        for(uint64_t i=0; i<n_bins-2; ++i)
        {
            x[i] = (i-1) * dx;
            igf[i] = boost::math::tgamma_lower(a, x[i]);
        }
        // for x -> infinity:
        x[n_bins-1] = (n_bins-1)*dx;
        igf[n_bins-1] = std::tgamma(a);

        // High resolution (lowest 2% of the x-values)
        dx_highres = x[std::round(0.01*(n_bins-1))] / (n_bins_highres - 1.0);
        odx_highres = 1.0/dx_highres;
        for(uint64_t i=0; i<n_bins_highres; ++i)
        {
            x_highres[i] = (i-1) * dx_highres;
            igf_highres[i] = boost::math::tgamma_lower(a, x_highres[i]);
        }
    }
};


/**
 * Structure for constants of a model. Includes particle constants as well.
 */
struct model_constants_t{
    //
    // Physical constants warm cloud
    //
    double alpha_d; /*!< Accomodation coefficient */

    codi::RealReverse Nc_prime; /*!< Number concentration of cloud droplets needed for one-moment scheme */

    codi::RealReverse a1_prime; /*!< Dimensional coefficient used in one-moment warm physics for qc and qr calculation */
    codi::RealReverse a2_prime; /*!< Dimensional coefficient used in one-moment warm physics for qc and qr calculation */
    codi::RealReverse e1_prime; /*!< Dimensional coefficients used in one-moment warm physics for temperature calculation */
    codi::RealReverse e2_prime; /*!< Dimensional coefficients used in one-moment warm physics for temperature calculation */
    codi::RealReverse B_prime;  /*!< Dimensional coefficient to simulate inflow from above in one-moment warm physics */
    codi::RealReverse d_prime;  /*!< Dimensional coefficient used in one-moment warm physics qr calculation for sedimentation*/

    codi::RealReverse dw; /*!< Change in buoancy */

    codi::RealReverse gamma;  /*!< Exponent used in one-moment warm physics for qc and qr calculation */
    codi::RealReverse betac;  /*!< Exponent used in one-moment warm physics for qc and qr calculation */
    codi::RealReverse betar;  /*!< Exponent used in one-moment warm physics for qc and qr calculation */
    codi::RealReverse delta1; /*!< Exponent used in one-moment warm physics for qv, qr, saturation and temperature calculation */
    codi::RealReverse delta2; /*!< Exponent used in one-moment warm physics for qv, qr, saturation and temperature calculation */
    codi::RealReverse zeta;   /*!< Exponents used in one-moment warm physics for qr calculation */

    /**
     * Model constants for hail.
     */
    particle_model_constants_t hail;
    /**
     * Model constants for ice.
     */
    particle_model_constants_t ice;
    /**
     * Model constants for snow.
     */
    particle_model_constants_t snow;
    /**
     * Model constants for cloud.
     */
    particle_model_constants_t cloud;
    /**
     * Model constants for rain.
     */
    particle_model_constants_t rain;
    /**
     * Model constants for graupel.
     */
    particle_model_constants_t graupel;

    /**
     * Snow - cloud droplet riming.
     */
    collection_model_constants_t coeffs_scr;
    /**
     * Snow - rain droplet riming.
     */
    collection_model_constants_t coeffs_srr;
    /**
     * Ice - rain droplet riming.
     */
    collection_model_constants_t coeffs_irr;
    /**
     * Ice - cloud droplet riming.
     */
    collection_model_constants_t coeffs_icr;
    /**
     * Hail - rain droplet riming.
     */
    collection_model_constants_t coeffs_hrr;
    /**
     * Graupel - rain droplet riming.
     */
    collection_model_constants_t coeffs_grr;
    /**
     * Hail - cloud droplet riming.
     */
    collection_model_constants_t coeffs_hcr;
    /**
     * Graupel - cloud droplet riming.
     */
    collection_model_constants_t coeffs_gcr;
    /**
     * Snow - ice collection.
     */
    collection_model_constants_t coeffs_sic;
    /**
     * Hail - ice collection.
     */
    collection_model_constants_t coeffs_hic;
    /**
     * Graupel - ice collection.
     */
    collection_model_constants_t coeffs_gic;
    /**
     * Hail - snow collection.
     */
    collection_model_constants_t coeffs_hsc;
    /**
     * Graupel - snow collection.
     */
    collection_model_constants_t coeffs_gsc;

    //
    // Technical constants
    //
    double t_end_prime;       /*!< End time in seconds for the simulation. */
    double t_end;             /*!< End time for the simulation. */
    double dt_prime;          /*!< Timestep size in seconds for the simulation. */
    double dt;                /*!< Timestep size for the simulation. */
    double dt_traject_prime;  /*!< Timestep size of the trajectory from the netCDF file. */
    double dt_traject;        /*!< Timestep size of the trajectory from the netCDF file. */
    uint64_t num_steps;       /*!< Number of timesteps to read from the netCDF file. */

    /**
     * Number of timesteps to simulate between each timestep of the netCDF file.
     */
    uint64_t num_sub_steps;

    //
    // General performance constants
    //
    double dt_half;   /*!< dt/2 */
    double dt_sixth;  /*!< dt/6 */
    double dt_third;  /*!< dt/3 */

    codi::RealReverse a1_scale; /*!< Performance constants warm cloud */
    codi::RealReverse a2_scale; /*!< Performance constants warm cloud */
    codi::RealReverse e1_scale; /*!< Performance constants warm cloud */
    codi::RealReverse e2_scale; /*!< Performance constants warm cloud */
    codi::RealReverse d_scale;  /*!< Performance constants warm cloud */
    codi::RealReverse rain_gfak = 1.0; /*!< Coefficient for gamma evaluation in rain evaporation */
    codi::RealReverse cloud_k_au; /*!< Coefficient for autoconversion of cloud to rain */
    codi::RealReverse cloud_k_sc; /*!< Coefficient for autoconversion of cloud to rain */

    /**
     * Kernel for autoconversion
     */
    codi::RealReverse kc_autocon = 9.44e9;

    /**
     * Inverse layer thickness. Used for sedimentation.
     * In Miltenberger (2016) the trajectories start every \f$100 \text{m}\f$
     * between the surface and \f$4 \text{km}\f$ altitude using COSMO-2, which
     * uses a mean spacing of \f$388 \text{m}\f$
     * with \f$13 \text{m}\f$ close to the surface and \f$1190 \text{m}\f$
     * at \f$23 \text{km}\f$.
     */
    codi::RealReverse inv_z;

    /// Parameters used in warm physics
    const double nar = 0.22;      /*!< Constants for the IFS model. */
    const double nbr = 2.2;       /*!< Constants for the IFS model. */
    const double ar = M_PI / 6.0; /*!< Constants for the IFS model. */
    const double br = 3.0;        /*!< Constants for the IFS model. */
    const double cr = 386.8;      /*!< Constants for the IFS model. */
    const double dr = 0.67;       /*!< Constants for the IFS model. */
    const double Sc = 0.6;        /*!< Constants for the IFS model. */
    const double mu = 16.0e-6;    /*!< Constants for the IFS model. */
    const double rho0 = 1.0;      /*!< Constants for the IFS model. Why not 1.225? */

    const double alpha_r = 1.0/(br + 1.0 - nbr);   /*!< Constants for the IFS model. */
    const double epsilonr = 0.5*dr + 2.5 - nbr;   /*!< Constants for the IFS model. */

    double scaling_fact; /*!< Scaling factor. */

    // See constants.h for a descripition of those.
    codi::RealReverse q_crit_i;
    codi::RealReverse D_crit_i;
    codi::RealReverse D_conv_i;
    codi::RealReverse q_crit_r;
    codi::RealReverse D_crit_r;
    codi::RealReverse q_crit_fr;
    codi::RealReverse D_coll_c;
    codi::RealReverse q_crit;
    codi::RealReverse D_conv_sg;
    codi::RealReverse D_conv_ig;
    codi::RealReverse x_conv;
    codi::RealReverse parcel_height;
    codi::RealReverse alpha_spacefilling;
    codi::RealReverse T_nuc;
    codi::RealReverse T_freeze;
    codi::RealReverse T_f;
    codi::RealReverse D_eq;
    codi::RealReverse rho_w;
    codi::RealReverse rho_0;
    codi::RealReverse rho_vel;
    codi::RealReverse rho_vel_c;
    codi::RealReverse rho_ice;
    codi::RealReverse M_w;
    codi::RealReverse M_a;
    codi::RealReverse R_universal;
    codi::RealReverse Epsilon;
    codi::RealReverse gravity_acc;
    codi::RealReverse R_a;
    codi::RealReverse R_v;
    codi::RealReverse a_v;
    codi::RealReverse b_v;
    codi::RealReverse a_prime;
    codi::RealReverse b_prime;
    codi::RealReverse c_prime;
    codi::RealReverse K_T;
    codi::RealReverse L_wd;
    codi::RealReverse L_ed;
    codi::RealReverse D_v;
    codi::RealReverse ecoll_min;
    codi::RealReverse ecoll_gg;
    codi::RealReverse ecoll_gg_wet;
    codi::RealReverse kin_visc_air;
    codi::RealReverse C_mult;
    codi::RealReverse T_mult_min;
    codi::RealReverse T_mult_max;
    codi::RealReverse T_mult_opt;

    /**
     * Constant used in cloud riming.
     */
    codi::RealReverse const0;
    /**
     * Hallet-Mossop ice multiplication.
     * Constant used in ice - x and snow - x riming.
     */
    codi::RealReverse const3;
    /**
     * Hallet-Mossop ice multiplication.
     * Constant used in ice - x and snow - x riming.
     */
    codi::RealReverse const4;
    /**
     * Constant for conversions ice -> graupel, snow -> graupel,
     * melting (used in riming).
     */
    codi::RealReverse const5;
    codi::RealReverse D_rainfrz_gh;
    codi::RealReverse D_rainfrz_ig;
    codi::RealReverse dv0;
    codi::RealReverse p_sat_melt;
    codi::RealReverse cp;
    codi::RealReverse k_b;
    codi::RealReverse a_HET;
    codi::RealReverse b_HET;
    codi::RealReverse N_sc;
    codi::RealReverse n_f;
    codi::RealReverse N_avo;
    codi::RealReverse na_dust;
    codi::RealReverse na_soot;
    codi::RealReverse na_orga;
    codi::RealReverse ni_het_max;
    codi::RealReverse ni_hom_max;
    codi::RealReverse a_dep;
    codi::RealReverse b_dep;
    codi::RealReverse c_dep;
    codi::RealReverse d_dep;
    codi::RealReverse nim_imm;
    codi::RealReverse nin_dep;
    codi::RealReverse alf_imm;
    codi::RealReverse bet_dep;
    codi::RealReverse bet_imm;
    std::vector<codi::RealReverse> a_ccn;
    std::vector<codi::RealReverse> b_ccn;
    std::vector<codi::RealReverse> c_ccn;
    std::vector<codi::RealReverse> d_ccn;
    codi::RealReverse r_const;
    codi::RealReverse r1_const;
    codi::RealReverse cv;
    codi::RealReverse p_sat_const_a;
    codi::RealReverse p_sat_ice_const_a;
    codi::RealReverse p_sat_const_b;
    codi::RealReverse p_sat_ice_const_b;
    codi::RealReverse p_sat_low_temp;
    codi::RealReverse T_sat_low_temp;
    codi::RealReverse alpha_depo;
    codi::RealReverse r_0;

    codi::RealReverse k_1_conv;
    codi::RealReverse k_2_conv;
    codi::RealReverse k_1_accr;
    codi::RealReverse k_r;

    /**
     * Register the model parameters on the tape for codi::RealReverse.
     *
     * @param tape Tape where the parameters should be registered on.
     */
    void register_input(codi::RealReverse::TapeType &tape)
    {
#if defined(RK4_ONE_MOMENT)
        // Dimensional coefficients
        tape.registerInput(this->a1_prime);    // Autoconversion
        tape.registerInput(this->a2_prime);    // Accretion
        tape.registerInput(this->e1_prime);    // Evaporation
        tape.registerInput(this->e2_prime);    // Evaporation
        tape.registerInput(this->d_prime);     // Sedimentation
        tape.registerInput(this->Nc_prime);    // Concentration of cloud droplets

        // Exponents
        tape.registerInput(this->gamma);       // Autoconversion
        tape.registerInput(this->betac);       // Accretion
        tape.registerInput(this->betar);       // Accretion
        tape.registerInput(this->delta1);      // Evaporation
        tape.registerInput(this->delta2);      // Evaporation
        tape.registerInput(this->zeta);        // Sedimentation

#elif defined(RK4ICE) || defined(RK4NOICE)
        // Dimensional coefficients
        tape.registerInput(this->a1_prime);    // Autoconversion
        tape.registerInput(this->a2_prime);    // Accretion
        tape.registerInput(this->e1_prime);    // Evaporation
        tape.registerInput(this->e2_prime);    // Evaporation
        tape.registerInput(this->d_prime);     // Sedimentation
        tape.registerInput(this->Nc_prime);    // Concentration of cloud droplets

        // Exponents
        tape.registerInput(this->gamma);       // Autoconversion
        tape.registerInput(this->betac);       // Accretion
        tape.registerInput(this->betar);       // Accretion
        tape.registerInput(this->delta1);      // Evaporation
        tape.registerInput(this->delta2);      // Evaporation
        tape.registerInput(this->zeta);        // Sedimentation

        // ICON parameters
        tape.registerInput(this->rain_gfak);
        tape.registerInput(this->cloud_k_au);
        tape.registerInput(this->cloud_k_sc);
        tape.registerInput(this->kc_autocon);
        tape.registerInput(this->inv_z);

        // Everything else
        tape.registerInput(this->dw);
        tape.registerInput(this->q_crit_i);
        tape.registerInput(this->D_crit_i);
        tape.registerInput(this->D_conv_i);
        tape.registerInput(this->q_crit_r);
        tape.registerInput(this->D_crit_r);
        tape.registerInput(this->q_crit_fr);
        tape.registerInput(this->D_coll_c);
        tape.registerInput(this->q_crit);
        tape.registerInput(this->D_conv_sg);
        tape.registerInput(this->D_conv_ig);
        tape.registerInput(this->x_conv);
        tape.registerInput(this->parcel_height);
        tape.registerInput(this->alpha_spacefilling);
        tape.registerInput(this->T_nuc);
        tape.registerInput(this->T_freeze);
        tape.registerInput(this->T_f);
        tape.registerInput(this->D_eq);
        tape.registerInput(this->rho_w);
        tape.registerInput(this->rho_0);
        tape.registerInput(this->rho_vel);
        tape.registerInput(this->rho_vel_c);
        tape.registerInput(this->rho_ice);
        tape.registerInput(this->M_w);
        tape.registerInput(this->M_a);
        tape.registerInput(this->R_universal);
        tape.registerInput(this->Epsilon);
        tape.registerInput(this->gravity_acc);
        tape.registerInput(this->R_a);
        tape.registerInput(this->R_v);
        tape.registerInput(this->a_v);
        tape.registerInput(this->b_v);
        tape.registerInput(this->a_prime);
        tape.registerInput(this->b_prime);
        tape.registerInput(this->c_prime);
        tape.registerInput(this->K_T);
        tape.registerInput(this->L_wd);
        tape.registerInput(this->L_ed);
        tape.registerInput(this->D_v);
        tape.registerInput(this->ecoll_min);
        tape.registerInput(this->ecoll_gg);
        tape.registerInput(this->ecoll_gg_wet);
        tape.registerInput(this->kin_visc_air);
        tape.registerInput(this->C_mult);
        tape.registerInput(this->T_mult_min);
        tape.registerInput(this->T_mult_max);
        tape.registerInput(this->T_mult_opt);
        tape.registerInput(this->const0);
        tape.registerInput(this->const3);
        tape.registerInput(this->const4);
        tape.registerInput(this->const5);
        tape.registerInput(this->D_rainfrz_ig);
        tape.registerInput(this->dv0);
        tape.registerInput(this->p_sat_melt);
        tape.registerInput(this->cp);
        tape.registerInput(this->k_b);
        tape.registerInput(this->a_HET);
        tape.registerInput(this->b_HET);
        tape.registerInput(this->N_sc);
        tape.registerInput(this->n_f);
        tape.registerInput(this->N_avo);
        tape.registerInput(this->na_dust);
        tape.registerInput(this->na_soot);
        tape.registerInput(this->na_orga);
        tape.registerInput(this->ni_het_max);
        tape.registerInput(this->ni_hom_max);
        tape.registerInput(this->a_dep);
        tape.registerInput(this->b_dep);
        tape.registerInput(this->c_dep);
        tape.registerInput(this->d_dep);
        tape.registerInput(this->nim_imm);
        tape.registerInput(this->nin_dep);
        tape.registerInput(this->alf_imm);
        tape.registerInput(this->bet_dep);
        tape.registerInput(this->bet_imm);
        tape.registerInput(this->r_const);
        tape.registerInput(this->r1_const);
        tape.registerInput(this->cv);
        tape.registerInput(this->p_sat_const_a);
        tape.registerInput(this->p_sat_ice_const_a);
        tape.registerInput(this->p_sat_const_b);
        tape.registerInput(this->p_sat_ice_const_b);
        tape.registerInput(this->p_sat_low_temp);
        tape.registerInput(this->T_sat_low_temp);
        tape.registerInput(this->alpha_depo);
        tape.registerInput(this->r_0);
        tape.registerInput(this->k_1_conv);
        tape.registerInput(this->k_2_conv);
        tape.registerInput(this->k_1_accr);
        tape.registerInput(this->k_r);
        for(auto &i: this->a_ccn)
            tape.registerInput(i);
        for(auto &i: this->b_ccn)
            tape.registerInput(i);
        for(auto &i: this->c_ccn)
            tape.registerInput(i);
        for(auto &i: this->d_ccn)
            tape.registerInput(i);
#endif
        this->rain.register_input(tape);
        this->cloud.register_input(tape);
        this->hail.register_input(tape);
        this->ice.register_input(tape);
        this->snow.register_input(tape);
        this->graupel.register_input(tape);
    }

    /**
     * Get the gradients of all its members. You need to register them on a
     * type before to get meaningful values.
     *
     * @param out_vec On out: Stores all gradients.
     */
    template<class T>
    void get_gradient(T &out_vec)
    {
#if defined(RK4_ONE_MOMENT)
        out_vec[0] = this->a1_prime.getGradient();
        out_vec[1] = this->a2_prime.getGradient();
        out_vec[2] = this->e1_prime.getGradient();
        out_vec[3] = this->e2_prime.getGradient();
        out_vec[4] = this->d_prime.getGradient();

        out_vec[5] = this->gamma.getGradient();
        out_vec[6] = this->betac.getGradient();
        out_vec[7] = this->betar.getGradient();
        out_vec[8] = this->delta1.getGradient();
        out_vec[9] = this->delta2.getGradient();
        out_vec[10] = this->zeta.getGradient();

        out_vec[11] = this->Nc_prime.getGradient();

#elif defined(RK4ICE) || defined(RK4NOICE)
        out_vec[0] = this->a1_prime.getGradient();
        out_vec[1] = this->a2_prime.getGradient();
        out_vec[2] = this->e1_prime.getGradient();
        out_vec[3] = this->e2_prime.getGradient();
        out_vec[4] = this->d_prime.getGradient();

        out_vec[5] = this->gamma.getGradient();
        out_vec[6] = this->betac.getGradient();
        out_vec[7] = this->betar.getGradient();
        out_vec[8] = this->delta1.getGradient();
        out_vec[9] = this->delta2.getGradient();
        out_vec[10] = this->zeta.getGradient();

        out_vec[11] = this->Nc_prime.getGradient();

        out_vec[12] = this->rain_gfak.getGradient();
        out_vec[13] = this->cloud_k_au.getGradient();
        out_vec[14] = this->cloud_k_sc.getGradient();
        out_vec[15] = this->kc_autocon.getGradient();
        out_vec[16] = this->inv_z.getGradient();

        out_vec[17] = this->dw.getGradient();
        out_vec[18] = this->q_crit_i.getGradient();
        out_vec[19] = this->D_crit_i.getGradient();
        out_vec[20] = this->D_conv_i.getGradient();
        out_vec[21] = this->q_crit_r.getGradient();
        out_vec[22] = this->D_crit_r.getGradient();
        out_vec[23] = this->q_crit_fr.getGradient();
        out_vec[24] = this->D_coll_c.getGradient();
        out_vec[25] = this->q_crit.getGradient();
        out_vec[26] = this->D_conv_sg.getGradient();
        out_vec[27] = this->D_conv_ig.getGradient();
        out_vec[28] = this->x_conv.getGradient();
        out_vec[29] = this->parcel_height.getGradient();
        out_vec[30] = this->alpha_spacefilling.getGradient();
        out_vec[31] = this->T_nuc.getGradient();
        out_vec[32] = this->T_freeze.getGradient();
        out_vec[33] = this->T_f.getGradient();
        out_vec[34] = this->D_eq.getGradient();
        out_vec[35] = this->rho_w.getGradient();
        out_vec[36] = this->rho_0.getGradient();
        out_vec[37] = this->rho_vel.getGradient();
        out_vec[38] = this->rho_vel_c.getGradient();
        out_vec[39] = this->rho_ice.getGradient();
        out_vec[40] = this->M_w.getGradient();
        out_vec[41] = this->M_a.getGradient();
        out_vec[42] = this->R_universal.getGradient();
        out_vec[43] = this->Epsilon.getGradient();
        out_vec[44] = this->gravity_acc.getGradient();
        out_vec[45] = this->R_a.getGradient();
        out_vec[46] = this->R_v.getGradient();
        out_vec[47] = this->a_v.getGradient();
        out_vec[48] = this->b_v.getGradient();
        out_vec[49] = this->a_prime.getGradient();
        out_vec[50] = this->b_prime.getGradient();
        out_vec[51] = this->c_prime.getGradient();
        out_vec[52] = this->K_T.getGradient();
        out_vec[53] = this->L_wd.getGradient();
        out_vec[54] = this->L_ed.getGradient();
        out_vec[55] = this->D_v.getGradient();
        out_vec[56] = this->ecoll_min.getGradient();
        out_vec[57] = this->ecoll_gg.getGradient();
        out_vec[58] = this->ecoll_gg_wet.getGradient();
        out_vec[59] = this->kin_visc_air.getGradient();
        out_vec[60] = this->C_mult.getGradient();
        out_vec[61] = this->T_mult_min.getGradient();
        out_vec[62] = this->T_mult_max.getGradient();
        out_vec[63] = this->T_mult_opt.getGradient();
        out_vec[64] = this->const0.getGradient();
        out_vec[65] = this->const3.getGradient();
        out_vec[66] = this->const4.getGradient();
        out_vec[67] = this->const5.getGradient();
        out_vec[68] = this->D_rainfrz_ig.getGradient();
        out_vec[69] = this->dv0.getGradient();
        out_vec[70] = this->p_sat_melt.getGradient();
        out_vec[71] = this->cp.getGradient();
        out_vec[72] = this->k_b.getGradient();
        out_vec[73] = this->a_HET.getGradient();
        out_vec[74] = this->b_HET.getGradient();
        out_vec[75] = this->N_sc.getGradient();
        out_vec[76] = this->n_f.getGradient();
        out_vec[77] = this->N_avo.getGradient();
        out_vec[78] = this->na_dust.getGradient();
        out_vec[79] = this->na_soot.getGradient();
        out_vec[80] = this->na_orga.getGradient();
        out_vec[81] = this->ni_het_max.getGradient();
        out_vec[82] = this->ni_hom_max.getGradient();
        out_vec[83] = this->a_dep.getGradient();
        out_vec[84] = this->b_dep.getGradient();
        out_vec[85] = this->c_dep.getGradient();
        out_vec[86] = this->d_dep.getGradient();
        out_vec[87] = this->nim_imm.getGradient();
        out_vec[88] = this->nin_dep.getGradient();
        out_vec[89] = this->alf_imm.getGradient();
        out_vec[90] = this->bet_dep.getGradient();
        out_vec[91] = this->bet_imm.getGradient();
        out_vec[92] = this->r_const.getGradient();
        out_vec[93] = this->r1_const.getGradient();
        out_vec[94] = this->cv.getGradient();
        out_vec[95] = this->p_sat_const_a.getGradient();
        out_vec[96] = this->p_sat_ice_const_a.getGradient();
        out_vec[97] = this->p_sat_const_b.getGradient();
        out_vec[98] = this->p_sat_ice_const_b.getGradient();
        out_vec[99] = this->p_sat_low_temp.getGradient();
        out_vec[100] = this->T_sat_low_temp.getGradient();
        out_vec[101] = this->alpha_depo.getGradient();
        out_vec[102] = this->r_0.getGradient();
        out_vec[103] = this->k_1_conv.getGradient();
        out_vec[104] = this->k_2_conv.getGradient();
        out_vec[105] = this->k_1_accr.getGradient();
        out_vec[106] = this->k_r.getGradient();


        uint64_t idx = 107;
        for(auto &i: this->a_ccn)
        {
            out_vec[idx] = i.getGradient();
            idx++;
        }
        for(auto &i: this->b_ccn)
        {
            out_vec[idx] = i.getGradient();
            idx++;
        }
        for(auto &i: this->c_ccn)
        {
            out_vec[idx] = i.getGradient();
            idx++;
        }
        for(auto &i: this->d_ccn)
        {
            out_vec[idx] = i.getGradient();
            idx++;
        }

        this->rain.get_gradient(out_vec, idx);
        this->cloud.get_gradient(out_vec, idx);
        this->graupel.get_gradient(out_vec, idx);
        this->hail.get_gradient(out_vec, idx);
        this->ice.get_gradient(out_vec, idx);
        this->snow.get_gradient(out_vec, idx);
#endif
    }
};


/**
 * Structure to collect all nc parameters.
 */
struct nc_parameters_t{

    uint64_t n_trajectories = 30; /*!< Number of trajectories in the netCDF file. */
    uint64_t n_timesteps = 7922; /*!< Number of timesteps in the netCDF file. */
#ifdef MET3D
    double z;
    std::vector<double> time_abs;
    uint64_t time_idx = 0;
#else
    std::vector<double> z;
#endif
#if defined MET3D && defined TURBULENCE
    NcVar qturb_var;
    double qturb;
#endif
    std::vector<double> w, lat, lon;
    double  t, p, time_rel,
            qc, qr, qi, qs, qg, qv, S, dw, dlat, dlon,
            QIin, QSin, QRin, QGin, QIout, QSout, QRout, QGout,
            NIin, NSin, NRin, NGin, NIout, NSout, NRout, NGout,
            Nc, Nr, Ni, Ns, Ng;
    bool ascent_flag, conv_400, conv_600, slan_400, slan_600, dp2h;
    char* type[1];
    NcVar   lat_var, lon_var, z_var, t_var, p_var, w_var, time_rel_var,
            qc_var, qr_var, qi_var, qs_var, qg_var, qv_var, S_var,
            QIin_var, QSin_var, QRin_var, QGin_var, QIout_var, QSout_var,
            QRout_var, QGout_var, ascent_flag_var,
            NIin_var, NSin_var, NRin_var, NGin_var, NIout_var, NSout_var,
            NRout_var, NGout_var, dp2h_var,
            Nc_var, Nr_var, Ni_var, Ns_var, Ng_var,
#ifdef MET3D
            type_var, time_abs_var,
#endif
            conv_400_var, conv_600_var, slan_400_var, slan_600_var;

};


/**
 * Helper structure to handle command line arguments.
 */
struct global_args_t{

  int final_time_flag; /*!< Using a final simulation time? */
  char* final_time_string;

  int timestep_flag; /*!< Timestep in seconds specified? */
  char* timestep_string;

  int snapshot_index_flag; /*!< Snapshot every x iterations specified? */
  char* snapshot_index_string;

  int output_flag; /*!< Output path specified? */
  char* output_string;

  int input_flag; /*!< Input netCDF file specified? */
  char* input_file;

  int scaling_fact_flag; /*!< Scaling factor specified? */
  char* scaling_fact_string;

  int start_over_flag; /*!< Reload mixing ratios and particle numbers from trajectory every few seconds? */
  char* start_over_string;

  int start_over_env_flag; /*!< Reload pressure, temperature and ascent from trajectory every few seconds? */
  char* start_over_env_string;

  int fixed_iteration_flag; /*!< Fix p, T, w during microphysics? */
  char* fixed_iteration_string;

  int auto_type_flag; /*!< Particle type specified? */
  char* auto_type_string;

  int traj_flag; /*!< Trajectory to use specified? */
  char* traj_string;

  int write_flag; /*!< Snapshot is flushed every x iterations. */
  char* write_string;

  int progress_index_flag; /*!< Progressbar is updated every x iterations. */
  char* progress_index_string;
#ifdef MET3D
  int delay_start_flag; /*!< Simulation starts at this time relative to ascend. */
  char* delay_start_string;
#endif
};


/**
 * Structure to collect all input parameters.
 */
struct input_parameters_t{

  // Numerics
  double t_end_prime; /*!< End simulation time in seconds. */
  double dt_prime; /*!< Timestep size in seconds for the simulation. */
  double dt_traject_prime; /*!< Timestep size in seconds of the trajectory in the netCDF file. */
  double dt_traject; /*!< Timestep size of the trajectory in the netCDF file. */
#ifdef MET3D
  double start_time;
#endif
  int snapshot_index; /*!< Save a snapshot every snapshot_index iteration. */
  /**
   * Number of timesteps for the simulation between two
   * datapoints from the netCDF file.
   */
  uint64_t num_sub_steps;
  std::string OUTPUT_FILENAME; /*!< Filename for output. */

  std::string INPUT_FILENAME; /*!< Filename for input netCDF file. */

  bool start_over; /*!< Start over at new timestep of trajectory? */
  bool start_over_env; /*!< Start over environment variables at new timestep of trajectory? */
  bool fixed_iteration; /*!< Fix temperature and pressure at every iteration? */

  double scaling_fact; /*!< Scaling factor. */

  uint32_t auto_type; /*!< Particle type. */
  uint32_t traj; /*!< Trajectory index to load from the netCDF file. */
  uint32_t write_index; /*!< Write stringstream every x iterations to disk. */
  uint32_t progress_index; /*!< Index for updating progressbar. */
  uint32_t ensemble = 0; /*!< Index of ensemble. */
};

struct param_t{
    double mean = std::nan("");
    double sigma = std::nan("");
    double sigma_perc = std::nan("");
    int err = 0;
    int name = -1;
    int out_name = -1;
    bool particle_param = false;
    // std::function<double()> dis;
    std::normal_distribution<double> dis;

    enum Param {
        a_1, a_2, e_1, e_2, d,
        N_c, model_gamma, beta_c, beta_r,
        delta1, delta2, zeta, rain_gfak,
        cloud_k_au, cloud_k_sc, kc_autocon, inv_z,
        w, q_crit_i, D_crit_i, D_conv_i,
        q_crit_r, D_crit_r, q_crit_fr, D_coll_c,
        q_crit, D_conv_sg, D_conv_ig, x_conv,
        parcel_height, alpha_spacefilling, T_nuc, T_freeze,
        T_f, D_eq, rho_w, rho_0,
        rho_vel, rho_vel_c, rho_ice, M_w,
        M_a, R_universal, Epsilon, gravity_acc,
        R_a, R_v, a_v, b_v,
        a_prime, b_prime, c_prime, K_T,
        L_wd, L_ed, D_v, ecoll_min,
        ecoll_gg, ecoll_gg_wet, kin_visc_air, C_mult,
        T_mult_min, T_mult_max, T_mult_opt, const0,
        const3, const4, const5, D_rainfrz_ig,
        dv0, p_sat_melt, cp, k_b,
        a_HET, b_HET, N_sc, n_f,
        N_avo, na_dust, na_soot, na_orga,
        ni_het_max, ni_hom_max, a_dep, b_dep,
        c_dep, d_dep, nim_imm, nin_dep,
        alf_imm, bet_dep, bet_imm, r_const,
        r1_const, cv, p_sat_const_a, p_sat_ice_const_a,
        p_sat_const_b, p_sat_ice_const_b, p_sat_low_temp, T_sat_low_temp,
        alpha_depo, r_0, k_1_conv, k_2_conv,
        k_1_accr, k_r, a_ccn_1, a_ccn_2,
        a_ccn_3, a_ccn_4, b_ccn_1, b_ccn_2,
        b_ccn_3, b_ccn_4, c_ccn_1, c_ccn_2,
        c_ccn_3, c_ccn_4, d_ccn_1, d_ccn_2,
        d_ccn_3, d_ccn_4
    };
    enum ParticleParam {
        a_geo, b_geo,
        min_x, min_x_act, min_x_nuc_homo, min_x_nuc_hetero,
        min_x_melt, min_x_evap, min_x_freezing, min_x_depo,
        min_x_collision, min_x_collection, min_x_conversion, min_x_sedimentation,
        min_x_riming, max_x, sc_theta_q, sc_delta_q,
        sc_theta_n, sc_delta_n, s_vel, a_vel,
        b_vel, rho_v, c_z, sc_coll_n,
        cmu0, cmu1, cmu2, cmu3,
        cmu4, cmu5, alpha, beta,
        gamma, nu, g1, g2,
        mu, nm1, nm2, nm3,
        q_crit_c, d_crit_c, ecoll_c, cap,
        a_ven, b_ven, c_s, a_f,
        b_f, alfa_n, alfa_q, lambda,
        vsedi_min, vsedi_max
    };
    std::unordered_map<std::string, Param> const table_param = {
        {"a_1", Param::a_1}, {"a_2", Param::a_2},
        {"e_1", Param::e_1}, {"e_2", Param::e_2},
        {"d", Param::d}, {"N_c", Param::N_c},
        {"gamma", Param::model_gamma}, {"beta_c", Param::beta_c},
        {"beta_r", Param::beta_r}, {"delta1", Param::delta1},
        {"delta2", Param::delta2}, {"zeta", Param::zeta},
        {"rain_gfak", Param::rain_gfak}, {"cloud_k_au", Param::cloud_k_au},
        {"cloud_k_sc", Param::cloud_k_sc}, {"kc_autocon", Param::kc_autocon},
        {"inv_z", Param::inv_z}, {"w", Param::w},
        {"q_crit_i", Param::q_crit_i}, {"D_crit_i", Param::D_crit_i},
        {"D_conv_i", Param::D_conv_i}, {"q_crit_r", Param::q_crit_r},
        {"D_crit_r", Param::D_crit_r}, {"q_crit_fr", Param::q_crit_fr},
        {"D_coll_c", Param::D_coll_c}, {"q_crit", Param::q_crit},
        {"D_conv_sg", Param::D_conv_sg}, {"D_conv_ig", Param::D_conv_ig},
        {"x_conv", Param::x_conv}, {"parcel_height", Param::parcel_height},
        {"alpha_spacefilling", Param::alpha_spacefilling}, {"T_nuc", Param::T_nuc},
        {"T_freeze", Param::T_freeze}, {"T_f", Param::T_f},
        {"D_eq", Param::D_eq}, {"rho_w", Param::rho_w},
        {"rho_0", Param::rho_0}, {"rho_vel", Param::rho_vel},
        {"rho_vel_c", Param::rho_vel_c}, {"rho_ice", Param::rho_ice},
        {"M_w", Param::M_w}, {"M_a", Param::M_a},
        {"R_universal", Param::R_universal}, {"Epsilon", Param::Epsilon},
        {"gravity_acc", Param::gravity_acc}, {"R_a", Param::R_a},
        {"R_v", Param::R_v}, {"a_v", Param::a_v},
        {"b_v", Param::b_v}, {"a_prime", Param::a_prime},
        {"b_prime", Param::b_prime}, {"c_prime", Param::c_prime},
        {"K_T", Param::K_T}, {"L_wd", Param::L_wd},
        {"L_ed", Param::L_ed}, {"D_v", Param::D_v},
        {"ecoll_min", Param::ecoll_min}, {"ecoll_gg", Param::ecoll_gg},
        {"ecoll_gg_wet", Param::ecoll_gg_wet}, {"kin_visc_air", Param::kin_visc_air},
        {"C_mult", Param::C_mult}, {"T_mult_min", Param::T_mult_min},
        {"T_mult_max", Param::T_mult_max}, {"T_mult_opt", Param::T_mult_opt},
        {"const0", Param::const0}, {"const3", Param::const3},
        {"const4", Param::const4}, {"const5", Param::const5},
        {"D_rainfrz_ig", Param::D_rainfrz_ig}, {"dv0", Param::dv0},
        {"p_sat_melt", Param::p_sat_melt}, {"cp", Param::cp},
        {"k_b", Param::k_b}, {"a_HET", Param::a_HET},
        {"b_HET", Param::b_HET}, {"N_sc", Param::N_sc},
        {"n_f", Param::n_f}, {"N_avo", Param::N_avo},
        {"na_dust", Param::na_dust}, {"na_soot", Param::na_soot},
        {"na_orga", Param::na_orga}, {"ni_het_max", Param::ni_het_max},
        {"ni_hom_max", Param::ni_hom_max}, {"a_dep", Param::a_dep},
        {"b_dep", Param::b_dep}, {"c_dep", Param::c_dep},
        {"d_dep", Param::d_dep}, {"nim_imm", Param::nim_imm},
        {"nin_dep", Param::nin_dep}, {"alf_imm", Param::alf_imm},
        {"bet_dep", Param::bet_dep}, {"bet_imm", Param::bet_imm},
        {"r_const", Param::r_const}, {"r1_const", Param::r1_const},
        {"cv", Param::cv}, {"p_sat_const_a", Param::p_sat_const_a},
        {"p_sat_ice_const_a", Param::p_sat_ice_const_a}, {"p_sat_const_b", Param::p_sat_const_b},
        {"p_sat_ice_const_b", Param::p_sat_ice_const_b}, {"p_sat_low_temp", Param::p_sat_low_temp},
        {"T_sat_low_temp", Param::T_sat_low_temp}, {"alpha_depo", Param::alpha_depo},
        {"r_0", Param::r_0}, {"k_1_conv", Param::k_1_conv},
        {"k_2_conv", Param::k_2_conv}, {"k_1_accr", Param::k_1_accr},
        {"k_r", Param::k_r}, {"a_ccn_1", Param::a_ccn_1},
        {"a_ccn_2", Param::a_ccn_2}, {"a_ccn_3", Param::a_ccn_3},
        {"a_ccn_4", Param::a_ccn_4}, {"b_ccn_1", Param::b_ccn_1},
        {"b_ccn_2", Param::b_ccn_2}, {"b_ccn_3", Param::b_ccn_3},
        {"b_ccn_4", Param::b_ccn_4}, {"c_ccn_1", Param::c_ccn_1},
        {"c_ccn_2", Param::c_ccn_2}, {"c_ccn_3", Param::c_ccn_3},
        {"c_ccn_4", Param::c_ccn_4}, {"d_ccn_1", Param::d_ccn_1},
        {"d_ccn_2", Param::d_ccn_2}, {"d_ccn_3", Param::d_ccn_3},
        {"d_ccn_4", Param::d_ccn_4}
    };
    std::unordered_map<std::string, ParticleParam> const table_particle_param = {
        {"a_geo", ParticleParam::a_geo},
        {"b_geo", ParticleParam::b_geo}, {"min_x", ParticleParam::min_x},
        {"min_x_act", ParticleParam::min_x_act}, {"min_x_nuc_homo", ParticleParam::min_x_nuc_homo},
        {"min_x_nuc_hetero", ParticleParam::min_x_nuc_hetero}, {"min_x_melt", ParticleParam::min_x_melt},
        {"min_x_evap", ParticleParam::min_x_evap}, {"min_x_freezing", ParticleParam::min_x_freezing},
        {"min_x_depo", ParticleParam::min_x_depo}, {"min_x_collision", ParticleParam::min_x_collision},
        {"min_x_collection", ParticleParam::min_x_collection}, {"min_x_conversion", ParticleParam::min_x_conversion},
        {"min_x_sedimentation", ParticleParam::min_x_sedimentation}, {"min_x_riming", ParticleParam::min_x_riming},
        {"max_x", ParticleParam::max_x}, {"sc_theta_q", ParticleParam::sc_theta_q},
        {"sc_delta_q", ParticleParam::sc_delta_q}, {"sc_theta_n", ParticleParam::sc_theta_n},
        {"sc_delta_n", ParticleParam::sc_delta_n}, {"s_vel", ParticleParam::s_vel},
        {"a_vel", ParticleParam::a_vel}, {"b_vel", ParticleParam::b_vel},
        {"rho_v", ParticleParam::rho_v}, {"c_z", ParticleParam::c_z},
        {"sc_coll_n", ParticleParam::sc_coll_n}, {"cmu0", ParticleParam::cmu0},
        {"cmu1", ParticleParam::cmu1}, {"cmu2", ParticleParam::cmu2},
        {"cmu3", ParticleParam::cmu3}, {"cmu4", ParticleParam::cmu4},
        {"cmu5", ParticleParam::cmu5}, {"alpha", ParticleParam::alpha},
        {"beta", ParticleParam::beta}, {"gamma", ParticleParam::gamma},
        {"nu", ParticleParam::nu}, {"g1", ParticleParam::g1},
        {"g2", ParticleParam::g2}, {"mu", ParticleParam::mu},
        {"nm1", ParticleParam::nm1}, {"nm2", ParticleParam::nm2},
        {"nm3", ParticleParam::nm3}, {"q_crit_c", ParticleParam::q_crit_c},
        {"d_crit_c", ParticleParam::d_crit_c}, {"ecoll_c", ParticleParam::ecoll_c},
        {"cap", ParticleParam::cap}, {"a_ven", ParticleParam::a_ven},
        {"b_ven", ParticleParam::b_ven}, {"c_s", ParticleParam::c_s},
        {"a_f", ParticleParam::a_f}, {"b_f", ParticleParam::b_f},
        {"alfa_n", ParticleParam::alfa_n}, {"alfa_q", ParticleParam::alfa_q},
        {"lambda", ParticleParam::lambda}, {"vsedi_min", ParticleParam::vsedi_min},
        {"vsedi_max", ParticleParam::vsedi_max}
    };

    enum OutParam {
        model, cloud, rain, ice, graupel, hail, snow
    };
    std::unordered_map<std::string, OutParam> const table_out_param = {
        {"model", OutParam::model},
        {"cloud", OutParam::cloud}, {"rain", OutParam::rain},
        {"ice", OutParam::ice}, {"graupel", OutParam::graupel},
        {"hail", OutParam::hail}, {"snow", OutParam::snow}
    };

    param_t(std::string param_type)
    {
        auto it = table_out_param.find(param_type);
        if(it != table_out_param.end())
        {
            out_name = it->second;
            if(out_name != model)
                particle_param = true;
        } else
        {
            err = OUTPARAM_CONFIG_ERR;
        }
    }

    void add_name(std::string n, model_constants_t &cc)
    {
        if(particle_param)
        {
            auto it = table_particle_param.find(n);
            if(it != table_particle_param.end())
            {
                particle_model_constants_t *pt;
                switch(out_name)
                {
                    case cloud:
                        pt = &(cc.cloud);
                        break;
                    case rain:
                        pt = &(cc.rain);
                        break;
                    case snow:
                        pt = &(cc.snow);
                        break;
                    case graupel:
                        pt = &(cc.graupel);
                        break;
                    case hail:
                        pt = &(cc.hail);
                        break;
                    case ice:
                        pt = &(cc.ice);
                        break;
                }
                name = it->second;
                switch(name)
                {
                    case a_geo:
                        mean = pt->a_geo.getValue();
                        break;
                    case b_geo:
                        mean = pt->b_geo.getValue();
                        break;
                    case min_x:
                        mean = pt->min_x.getValue();
                        break;
                    case min_x_act:
                        mean = pt->min_x_act.getValue();
                        break;
                    case min_x_nuc_homo:
                        mean = pt->min_x_nuc_homo.getValue();
                        break;
                    case min_x_nuc_hetero:
                        mean = pt->min_x_nuc_hetero.getValue();
                        break;
                    case min_x_melt:
                        mean = pt->min_x_melt.getValue();
                        break;
                    case min_x_evap:
                        mean = pt->min_x_evap.getValue();
                        break;
                    case min_x_freezing:
                        mean = pt->min_x_freezing.getValue();
                        break;
                    case min_x_depo:
                        mean = pt->min_x_depo.getValue();
                        break;
                    case min_x_collision:
                        mean = pt->min_x_collision.getValue();
                        break;
                    case min_x_collection:
                        mean = pt->min_x_collection.getValue();
                        break;
                    case min_x_conversion:
                        mean = pt->min_x_conversion.getValue();
                        break;
                    case min_x_sedimentation:
                        mean = pt->min_x_sedimentation.getValue();
                        break;
                    case min_x_riming:
                        mean = pt->min_x_riming.getValue();
                        break;
                    case max_x:
                        mean = pt->max_x.getValue();
                        break;
                    case sc_theta_q:
                        mean = pt->sc_theta_q.getValue();
                        break;
                    case sc_delta_q:
                        mean = pt->sc_delta_q.getValue();
                        break;
                    case sc_theta_n:
                        mean = pt->sc_theta_n.getValue();
                        break;
                    case sc_delta_n:
                        mean = pt->sc_delta_n.getValue();
                        break;
                    case s_vel:
                        mean = pt->s_vel.getValue();
                        break;
                    case a_vel:
                        mean = pt->a_vel.getValue();
                        break;
                    case b_vel:
                        mean = pt->b_vel.getValue();
                        break;
                    case rho_v:
                        mean = pt->rho_v.getValue();
                        break;
                    case c_z:
                        mean = pt->c_z.getValue();
                        break;
                    case sc_coll_n:
                        mean = pt->sc_coll_n.getValue();
                        break;
                    case cmu0:
                        mean = pt->cmu0.getValue();
                        break;
                    case cmu1:
                        mean = pt->cmu1.getValue();
                        break;
                    case cmu2:
                        mean = pt->cmu2.getValue();
                        break;
                    case cmu3:
                        mean = pt->cmu3.getValue();
                        break;
                    case cmu4:
                        mean = pt->cmu4.getValue();
                        break;
                    case cmu5:
                        mean = pt->cmu5.getValue();
                        break;
                    case alpha:
                        mean = pt->alpha.getValue();
                        break;
                    case beta:
                        mean = pt->beta.getValue();
                        break;
                    case gamma:
                        mean = pt->gamma.getValue();
                        break;
                    case nu:
                        mean = pt->nu.getValue();
                        break;
                    case g1:
                        mean = pt->g1.getValue();
                        break;
                    case g2:
                        mean = pt->g2.getValue();
                        break;
                    case mu:
                        mean = pt->mu.getValue();
                        break;
                    case nm1:
                        mean = pt->nm1.getValue();
                        break;
                    case nm2:
                        mean = pt->nm2.getValue();
                        break;
                    case nm3:
                        mean = pt->nm3.getValue();
                        break;
                    case q_crit_c:
                        mean = pt->q_crit_c.getValue();
                        break;
                    case d_crit_c:
                        mean = pt->d_crit_c.getValue();
                        break;
                    case ecoll_c:
                        mean = pt->ecoll_c.getValue();
                        break;
                    case cap:
                        mean = pt->cap.getValue();
                        break;
                    case a_ven:
                        mean = pt->a_ven.getValue();
                        break;
                    case b_ven:
                        mean = pt->b_ven.getValue();
                        break;
                    case c_s:
                        mean = pt->c_s.getValue();
                        break;
                    case a_f:
                        mean = pt->a_f.getValue();
                        break;
                    case b_f:
                        mean = pt->b_f.getValue();
                        break;
                    case alfa_n:
                        mean = pt->alfa_n.getValue();
                        break;
                    case alfa_q:
                        mean = pt->alfa_q.getValue();
                        break;
                    case lambda:
                        mean = pt->lambda.getValue();
                        break;
                    case vsedi_min:
                        mean = pt->vsedi_min.getValue();
                        break;
                    case vsedi_max:
                        mean = pt->vsedi_max.getValue();
                        break;
                }
            } else
            {
                err = PARAM_CONFIG_ERR;
            }
        }
        else
        {
            auto it = table_param.find(n);
            if(it != table_param.end())
            {
                name = it->second;
                switch(name)
                {
                    case a_1:
                        mean = cc.a1_prime.getValue();
                        break;
                    case a_2:
                        mean = cc.a2_prime.getValue();
                        break;
                    case e_1:
                        mean = cc.e1_prime.getValue();
                        break;
                    case e_2:
                        mean = cc.e2_prime.getValue();
                        break;
                    case d:
                        mean = cc.d_prime.getValue();
                        break;
                    case N_c:
                        mean = cc.Nc_prime.getValue();
                        break;
                    case model_gamma:
                        mean = cc.gamma.getValue();
                        break;
                    case beta_c:
                        mean = cc.betac.getValue();
                        break;
                    case beta_r:
                        mean = cc.betar.getValue();
                        break;
                    case delta1:
                        mean = cc.delta1.getValue();
                        break;
                    case delta2:
                        mean = cc.delta2.getValue();
                        break;
                    case zeta:
                        mean = cc.zeta.getValue();
                        break;
                    case rain_gfak:
                        mean = cc.rain_gfak.getValue();
                        break;
                    case cloud_k_au:
                        mean = cc.cloud_k_au.getValue();
                        break;
                    case cloud_k_sc:
                        mean = cc.cloud_k_sc.getValue();
                        break;
                    case kc_autocon:
                        mean = cc.kc_autocon.getValue();
                        break;
                    case inv_z:
                        mean = cc.inv_z.getValue();
                        break;
                    case w:
                        mean = cc.dw.getValue();
                        break;
                    case q_crit_i:
                        mean = cc.q_crit_i.getValue();
                        break;
                    case D_crit_i:
                        mean = cc.D_crit_i.getValue();
                        break;
                    case D_conv_i:
                        mean = cc.D_conv_i.getValue();
                        break;
                    case q_crit_r:
                        mean = cc.q_crit_r.getValue();
                        break;
                    case D_crit_r:
                        mean = cc.D_crit_r.getValue();
                        break;
                    case q_crit_fr:
                        mean = cc.q_crit_fr.getValue();
                        break;
                    case D_coll_c:
                        mean = cc.D_coll_c.getValue();
                        break;
                    case q_crit:
                        mean = cc.q_crit.getValue();
                        break;
                    case D_conv_sg:
                        mean = cc.D_conv_sg.getValue();
                        break;
                    case D_conv_ig:
                        mean = cc.D_conv_ig.getValue();
                        break;
                    case x_conv:
                        mean = cc.x_conv.getValue();
                        break;
                    case parcel_height:
                        mean = cc.parcel_height.getValue();
                        break;
                    case alpha_spacefilling:
                        mean = cc.alpha_spacefilling.getValue();
                        break;
                    case T_nuc:
                        mean = cc.T_nuc.getValue();
                        break;
                    case T_freeze:
                        mean = cc.T_freeze.getValue();
                        break;
                    case T_f:
                        mean = cc.T_f.getValue();
                        break;
                    case D_eq:
                        mean = cc.D_eq.getValue();
                        break;
                    case rho_w:
                        mean = cc.rho_w.getValue();
                        break;
                    case rho_0:
                        mean = cc.rho_0.getValue();
                        break;
                    case rho_vel:
                        mean = cc.rho_vel.getValue();
                        break;
                    case rho_vel_c:
                        mean = cc.rho_vel_c.getValue();
                        break;
                    case rho_ice:
                        mean = cc.rho_ice.getValue();
                        break;
                    case M_w:
                        mean = cc.M_w.getValue();
                        break;
                    case M_a:
                        mean = cc.M_a.getValue();
                        break;
                    case R_universal:
                        mean = cc.R_universal.getValue();
                        break;
                    case Epsilon:
                        mean = cc.Epsilon.getValue();
                        break;
                    case gravity_acc:
                        mean = cc.gravity_acc.getValue();
                        break;
                    case R_a:
                        mean = cc.R_a.getValue();
                        break;
                    case R_v:
                        mean = cc.R_v.getValue();
                        break;
                    case a_v:
                        mean = cc.a_v.getValue();
                        break;
                    case b_v:
                        mean = cc.b_v.getValue();
                        break;
                    case a_prime:
                        mean = cc.a_prime.getValue();
                        break;
                    case b_prime:
                        mean = cc.b_prime.getValue();
                        break;
                    case c_prime:
                        mean = cc.c_prime.getValue();
                        break;
                    case K_T:
                        mean = cc.K_T.getValue();
                        break;
                    case L_wd:
                        mean = cc.L_wd.getValue();
                        break;
                    case L_ed:
                        mean = cc.L_ed.getValue();
                        break;
                    case D_v:
                        mean = cc.D_v.getValue();
                        break;
                    case ecoll_min:
                        mean = cc.ecoll_min.getValue();
                        break;
                    case ecoll_gg:
                        mean = cc.ecoll_gg.getValue();
                        break;
                    case ecoll_gg_wet:
                        mean = cc.ecoll_gg_wet.getValue();
                        break;
                    case kin_visc_air:
                        mean = cc.kin_visc_air.getValue();
                        break;
                    case C_mult:
                        mean = cc.C_mult.getValue();
                        break;
                    case T_mult_min:
                        mean = cc.T_mult_min.getValue();
                        break;
                    case T_mult_max:
                        mean = cc.T_mult_max.getValue();
                        break;
                    case T_mult_opt:
                        mean = cc.T_mult_opt.getValue();
                        break;
                    case const0:
                        mean = cc.const0.getValue();
                        break;
                    case const3:
                        mean = cc.const3.getValue();
                        break;
                    case const4:
                        mean = cc.const4.getValue();
                        break;
                    case const5:
                        mean = cc.const5.getValue();
                        break;
                    case D_rainfrz_ig:
                        mean = cc.D_rainfrz_ig.getValue();
                        break;
                    case dv0:
                        mean = cc.dv0.getValue();
                        break;
                    case p_sat_melt:
                        mean = cc.p_sat_melt.getValue();
                        break;
                    case cp:
                        mean = cc.cp.getValue();
                        break;
                    case k_b:
                        mean = cc.k_b.getValue();
                        break;
                    case a_HET:
                        mean = cc.a_HET.getValue();
                        break;
                    case b_HET:
                        mean = cc.b_HET.getValue();
                        break;
                    case N_sc:
                        mean = cc.N_sc.getValue();
                        break;
                    case n_f:
                        mean = cc.n_f.getValue();
                        break;
                    case N_avo:
                        mean = cc.N_avo.getValue();
                        break;
                    case na_dust:
                        mean = cc.na_dust.getValue();
                        break;
                    case na_soot:
                        mean = cc.na_soot.getValue();
                        break;
                    case na_orga:
                        mean = cc.na_orga.getValue();
                        break;
                    case ni_het_max:
                        mean = cc.ni_het_max.getValue();
                        break;
                    case ni_hom_max:
                        mean = cc.ni_hom_max.getValue();
                        break;
                    case a_dep:
                        mean = cc.a_dep.getValue();
                        break;
                    case b_dep:
                        mean = cc.b_dep.getValue();
                        break;
                    case c_dep:
                        mean = cc.c_dep.getValue();
                        break;
                    case d_dep:
                        mean = cc.d_dep.getValue();
                        break;
                    case nim_imm:
                        mean = cc.nim_imm.getValue();
                        break;
                    case nin_dep:
                        mean = cc.nin_dep.getValue();
                        break;
                    case alf_imm:
                        mean = cc.alf_imm.getValue();
                        break;
                    case bet_dep:
                        mean = cc.bet_dep.getValue();
                        break;
                    case bet_imm:
                        mean = cc.bet_imm.getValue();
                        break;
                    case r_const:
                        mean = cc.r_const.getValue();
                        break;
                    case r1_const:
                        mean = cc.r1_const.getValue();
                        break;
                    case cv:
                        mean = cc.cv.getValue();
                        break;
                    case p_sat_const_a:
                        mean = cc.p_sat_const_a.getValue();
                        break;
                    case p_sat_ice_const_a:
                        mean = cc.p_sat_ice_const_a.getValue();
                        break;
                    case p_sat_const_b:
                        mean = cc.p_sat_const_b.getValue();
                        break;
                    case p_sat_ice_const_b:
                        mean = cc.p_sat_ice_const_b.getValue();
                        break;
                    case p_sat_low_temp:
                        mean = cc.p_sat_low_temp.getValue();
                        break;
                    case T_sat_low_temp:
                        mean = cc.T_sat_low_temp.getValue();
                        break;
                    case alpha_depo:
                        mean = cc.alpha_depo.getValue();
                        break;
                    case r_0:
                        mean = cc.r_0.getValue();
                        break;
                    case k_1_conv:
                        mean = cc.k_1_conv.getValue();
                        break;
                    case k_2_conv:
                        mean = cc.k_2_conv.getValue();
                        break;
                    case k_1_accr:
                        mean = cc.k_1_accr.getValue();
                        break;
                    case k_r:
                        mean = cc.k_r.getValue();
                        break;
                    case a_ccn_1:
                        mean = cc.a_ccn[0].getValue();
                        break;
                    case a_ccn_2:
                        mean = cc.a_ccn[1].getValue();
                        break;
                    case a_ccn_3:
                        mean = cc.a_ccn[2].getValue();
                        break;
                    case a_ccn_4:
                        mean = cc.a_ccn[3].getValue();
                        break;
                    case b_ccn_1:
                        mean = cc.b_ccn[0].getValue();
                        break;
                    case b_ccn_2:
                        mean = cc.b_ccn[1].getValue();
                        break;
                    case b_ccn_3:
                        mean = cc.b_ccn[2].getValue();
                        break;
                    case b_ccn_4:
                        mean = cc.b_ccn[3].getValue();
                        break;
                    case c_ccn_1:
                        mean = cc.c_ccn[0].getValue();
                        break;
                    case c_ccn_2:
                        mean = cc.c_ccn[1].getValue();
                        break;
                    case c_ccn_3:
                        mean = cc.c_ccn[2].getValue();
                        break;
                    case c_ccn_4:
                        mean = cc.c_ccn[3].getValue();
                        break;
                    case d_ccn_1:
                        mean = cc.d_ccn[0].getValue();
                        break;
                    case d_ccn_2:
                        mean = cc.d_ccn[1].getValue();
                        break;
                    case d_ccn_3:
                        mean = cc.d_ccn[2].getValue();
                        break;
                    case d_ccn_4:
                        mean = cc.d_ccn[3].getValue();
                        break;
                }
            } else
            {
                err = PARAM_CONFIG_ERR;
            }
        }
        if(!isnan(sigma) && err != PARAM_CONFIG_ERR)
        {
            add_sigma(sigma);
        }
        else if(!isnan(sigma_perc) && err != PARAM_CONFIG_ERR)
        {
            add_sigma_perc(sigma_perc);
        }
    }

    void add_sigma(double s)
    {
        if(!isnan(mean))
        {
            // dis.param(std::normal_distribution<double>::param_type(mean, s));
            dis = std::normal_distribution<double>(mean, s);
            sigma = s;
        }
        else
            sigma = s;
    }

    void add_sigma_perc(double s)
    {
        if(!isnan(mean))
        {
            dis = std::normal_distribution<double>(mean, s*mean);
            sigma_perc = s;
            sigma = s*mean;
        }
        else
            sigma_perc = s;
    }

    int check()
    {
        if(name == -1)
        {
            std::cout << "Error in config file:\n"
                      << "You did not specify the parameter to perturb or "
                      << "you have a typo at <name>typo</name>\n"
                      << "Skipping this parameter.\n";
            err = MISSING_PARAM_CONFIG_ERR;
            return err;
        }
        if(isnan(sigma) && isnan(sigma_perc))
        {
            std::cout << "Error in config file:\n"
                      << "You did not specify the variance for "
                      << "perturbing the parameter.\n";
            err = MISSING_VARIANCE_CONFIG_ERR;
            return err;
        }

        switch(err)
        {
            case PARAM_CONFIG_ERR:
                std::cout << "Error in config file:\n"
                          << "You used a parameter to perturb that does "
                          << "not exist.\n"
                          << "Skipping this parameter.\n";
                return err;

            case OUTPARAM_CONFIG_ERR:
                std::cout << "Error in config file:\n"
                          << "You used an output parameter that does "
                          << "not exist.\n"
                          << "Skipping this parameter.\n";
                return err;

            default:
                return err;

        }
    }
};

struct segment_t
{
    std::vector<param_t> params;
    double value = std::nan("");
    uint32_t n_members = 1;
    int value_name = -1;
    int out_param = -1;
    uint32_t n_segments = 1;
    uint32_t err = 0;
    byte old_sign = 0;

    enum Method {
        impact_change, sign_flip, value_method
    };
    std::unordered_map<std::string, Method> const table_method = {
        {"impact_change", Method::impact_change}, {"sign_flip", Method::sign_flip},
        {"value", Method::value_method}
    };
    int method = value_method;
    enum Param {
        pressure, T, w, S, QC, QR, QV, NCCLOUD, NCRAIN,
        QI, NCICE, QS, NCSNOW, QG, NCGRAUPEL, QH, NCHAIL,
        QI_OUT, QS_OUT, QR_OUT, QG_OUT, QH_OUT, latent_heat,
        latent_cool, NI_OUT, NS_OUT, NR_OUT, NG_OUT,
        NH_OUT, z, Inactive, eposition, sublimination,
        a_1, a_2, e_1, e_2, d, N_c, gamma, beta_c,
        beta_r, delta1, delta2, zeta, rain_gfak, cloud_k_au,
        cloud_k_sc, kc_autocon, inv_z, dw, q_crit_i,
        D_crit_i, D_conv_i, q_crit_r, D_crit_r, q_crit_fr,
        D_coll_c, q_crit, D_conv_sg, D_conv_ig, x_conv,
        parcel_height, alpha_spacefilling, T_nuc, T_freeze,
        T_f, D_eq, rho_w, rho_0, rho_vel, rho_vel_c,
        rho_ice, M_w, M_a, R_universal, Epsilon, gravity_acc,
        R_a, R_v, a_v, b_v, a_prime, b_prime, c_prime,
        K_T, L_wd, L_ed, D_v, ecoll_min, ecoll_gg,
        ecoll_gg_wet, kin_visc_air, C_mult, T_mult_min,
        T_mult_max, T_mult_opt, const0, const3, const4,
        const5, D_rainfrz_ig, dv0, p_sat_melt, cp,
        k_b, a_HET, b_HET, N_sc, n_f, N_avo, na_dust,
        na_soot, na_orga, ni_het_max, ni_hom_max, a_dep,
        b_dep, c_dep, d_dep, nim_imm, nin_dep, alf_imm,
        bet_dep, bet_imm, r_const, r1_const, cv, p_sat_const_a,
        p_sat_ice_const_a, p_sat_const_b, p_sat_ice_const_b,
        p_sat_low_temp, T_sat_low_temp, alpha_depo,
        r_0, k_1_conv, k_2_conv, k_1_accr, k_r, a_ccn_1,
        a_ccn_2, a_ccn_3, a_ccn_4, b_ccn_1, b_ccn_2,
        b_ccn_3, b_ccn_4, c_ccn_1, c_ccn_2, c_ccn_3,
        c_ccn_4, d_ccn_1, d_ccn_2, d_ccn_3, d_ccn_4,
        rain_a_geo, rain_b_geo, rain_min_x, rain_min_x_act,
        rain_min_x_nuc_homo, rain_min_x_nuc_hetero,
        rain_min_x_melt, rain_min_x_evap, rain_min_x_freezing,
        rain_min_x_depo, rain_min_x_collision, rain_min_x_collection,
        rain_min_x_conversion, rain_min_x_sedimentation,
        rain_min_x_riming, rain_max_x, rain_sc_theta_q,
        rain_sc_delta_q, rain_sc_theta_n, rain_sc_delta_n,
        rain_s_vel, rain_a_vel, rain_b_vel, rain_rho_v,
        rain_c_z, rain_sc_coll_n, rain_cmu0, rain_cmu1,
        rain_cmu2, rain_cmu3, rain_cmu4, rain_cmu5,
        rain_alpha, rain_beta, rain_gamma, rain_nu,
        rain_g1, rain_g2, rain_mu, rain_nm1, rain_nm2,
        rain_nm3, rain_q_crit_c, rain_d_crit_c, rain_ecoll_c,
        rain_cap, rain_a_ven, rain_b_ven, rain_c_s,
        rain_a_f, rain_b_f, rain_alfa_n, rain_alfa_q,
        rain_lambda, rain_vsedi_min, rain_vsedi_max,
        cloud_a_geo, cloud_b_geo, cloud_min_x, cloud_min_x_act,
        cloud_min_x_nuc_homo, cloud_min_x_nuc_hetero,
        cloud_min_x_melt, cloud_min_x_evap, cloud_min_x_freezing,
        cloud_min_x_depo, cloud_min_x_collision, cloud_min_x_collection,
        cloud_min_x_conversion, cloud_min_x_sedimentation,
        cloud_min_x_riming, cloud_max_x, cloud_sc_theta_q,
        cloud_sc_delta_q, cloud_sc_theta_n, cloud_sc_delta_n,
        cloud_s_vel, cloud_a_vel, cloud_b_vel, cloud_rho_v,
        cloud_c_z, cloud_sc_coll_n, cloud_cmu0, cloud_cmu1,
        cloud_cmu2, cloud_cmu3, cloud_cmu4, cloud_cmu5,
        cloud_alpha, cloud_beta, cloud_gamma, cloud_nu,
        cloud_g1, cloud_g2, cloud_mu, cloud_nm1, cloud_nm2,
        cloud_nm3, cloud_q_crit_c, cloud_d_crit_c,
        cloud_ecoll_c, cloud_cap, cloud_a_ven, cloud_b_ven,
        cloud_c_s, cloud_a_f, cloud_b_f, cloud_alfa_n,
        cloud_alfa_q, cloud_lambda, cloud_vsedi_min,
        cloud_vsedi_max, graupel_a_geo, graupel_b_geo,
        graupel_min_x, graupel_min_x_act, graupel_min_x_nuc_homo,
        graupel_min_x_nuc_hetero, graupel_min_x_melt,
        graupel_min_x_evap, graupel_min_x_freezing,
        graupel_min_x_depo, graupel_min_x_collision,
        graupel_min_x_collection, graupel_min_x_conversion,
        graupel_min_x_sedimentation, graupel_min_x_riming,
        graupel_max_x, graupel_sc_theta_q, graupel_sc_delta_q,
        graupel_sc_theta_n, graupel_sc_delta_n, graupel_s_vel,
        graupel_a_vel, graupel_b_vel, graupel_rho_v,
        graupel_c_z, graupel_sc_coll_n, graupel_cmu0,
        graupel_cmu1, graupel_cmu2, graupel_cmu3, graupel_cmu4,
        graupel_cmu5, graupel_alpha, graupel_beta,
        graupel_gamma, graupel_nu, graupel_g1, graupel_g2,
        graupel_mu, graupel_nm1, graupel_nm2, graupel_nm3,
        graupel_q_crit_c, graupel_d_crit_c, graupel_ecoll_c,
        graupel_cap, graupel_a_ven, graupel_b_ven,
        graupel_c_s, graupel_a_f, graupel_b_f, graupel_alfa_n,
        graupel_alfa_q, graupel_lambda, graupel_vsedi_min,
        graupel_vsedi_max, hail_a_geo, hail_b_geo,
        hail_min_x, hail_min_x_act, hail_min_x_nuc_homo,
        hail_min_x_nuc_hetero, hail_min_x_melt, hail_min_x_evap,
        hail_min_x_freezing, hail_min_x_depo, hail_min_x_collision,
        hail_min_x_collection, hail_min_x_conversion,
        hail_min_x_sedimentation, hail_min_x_riming,
        hail_max_x, hail_sc_theta_q, hail_sc_delta_q,
        hail_sc_theta_n, hail_sc_delta_n, hail_s_vel,
        hail_a_vel, hail_b_vel, hail_rho_v, hail_c_z,
        hail_sc_coll_n, hail_cmu0, hail_cmu1, hail_cmu2,
        hail_cmu3, hail_cmu4, hail_cmu5, hail_alpha,
        hail_beta, hail_gamma, hail_nu, hail_g1, hail_g2,
        hail_mu, hail_nm1, hail_nm2, hail_nm3, hail_q_crit_c,
        hail_d_crit_c, hail_ecoll_c, hail_cap, hail_a_ven,
        hail_b_ven, hail_c_s, hail_a_f, hail_b_f, hail_alfa_n,
        hail_alfa_q, hail_lambda, hail_vsedi_min, hail_vsedi_max,
        ice_a_geo, ice_b_geo, ice_min_x, ice_min_x_act,
        ice_min_x_nuc_homo, ice_min_x_nuc_hetero, ice_min_x_melt,
        ice_min_x_evap, ice_min_x_freezing, ice_min_x_depo,
        ice_min_x_collision, ice_min_x_collection,
        ice_min_x_conversion, ice_min_x_sedimentation,
        ice_min_x_riming, ice_max_x, ice_sc_theta_q,
        ice_sc_delta_q, ice_sc_theta_n, ice_sc_delta_n,
        ice_s_vel, ice_a_vel, ice_b_vel, ice_rho_v,
        ice_c_z, ice_sc_coll_n, ice_cmu0, ice_cmu1,
        ice_cmu2, ice_cmu3, ice_cmu4, ice_cmu5, ice_alpha,
        ice_beta, ice_gamma, ice_nu, ice_g1, ice_g2,
        ice_mu, ice_nm1, ice_nm2, ice_nm3, ice_q_crit_c,
        ice_d_crit_c, ice_ecoll_c, ice_cap, ice_a_ven,
        ice_b_ven, ice_c_s, ice_a_f, ice_b_f, ice_alfa_n,
        ice_alfa_q, ice_lambda, ice_vsedi_min, ice_vsedi_max,
        snow_a_geo, snow_b_geo, snow_min_x, snow_min_x_act,
        snow_min_x_nuc_homo, snow_min_x_nuc_hetero,
        snow_min_x_melt, snow_min_x_evap, snow_min_x_freezing,
        snow_min_x_depo, snow_min_x_collision, snow_min_x_collection,
        snow_min_x_conversion, snow_min_x_sedimentation,
        snow_min_x_riming, snow_max_x, snow_sc_theta_q,
        snow_sc_delta_q, snow_sc_theta_n, snow_sc_delta_n,
        snow_s_vel, snow_a_vel, snow_b_vel, snow_rho_v,
        snow_c_z, snow_sc_coll_n, snow_cmu0, snow_cmu1,
        snow_cmu2, snow_cmu3, snow_cmu4, snow_cmu5,
        snow_alpha, snow_beta, snow_gamma, snow_nu,
        snow_g1, snow_g2, snow_mu, snow_nm1, snow_nm2,
        snow_nm3, snow_q_crit_c, snow_d_crit_c, snow_ecoll_c,
        snow_cap, snow_a_ven, snow_b_ven, snow_c_s,
        snow_a_f, snow_b_f, snow_alfa_n, snow_alfa_q,
        snow_lambda, snow_vsedi_min, snow_vsedi_max
    };
    std::unordered_map<std::string, Param> const table_param = {
        {"pressure", Param::pressure}, {"T", Param::T},
        {"w", Param::w}, {"S", Param::S},
        {"QC", Param::QC}, {"QR", Param::QR},
        {"QV", Param::QV}, {"NCCLOUD", Param::NCCLOUD},
        {"NCRAIN", Param::NCRAIN}, {"QI", Param::QI},
        {"NCICE", Param::NCICE}, {"QS", Param::QS},
        {"NCSNOW", Param::NCSNOW}, {"QG", Param::QG},
        {"NCGRAUPEL", Param::NCGRAUPEL}, {"QH", Param::QH},
        {"NCHAIL", Param::NCHAIL}, {"QI_OUT", Param::QI_OUT},
        {"QS_OUT", Param::QS_OUT}, {"QR_OUT", Param::QR_OUT},
        {"QG_OUT", Param::QG_OUT}, {"QH_OUT", Param::QH_OUT},
        {"latent_heat", Param::latent_heat}, {"latent_cool", Param::latent_cool},
        {"NI_OUT", Param::NI_OUT}, {"NS_OUT", Param::NS_OUT},
        {"NR_OUT", Param::NR_OUT}, {"NG_OUT", Param::NG_OUT},
        {"NH_OUT", Param::NH_OUT}, {"z", Param::z},
        {"Inactive", Param::Inactive}, {"eposition", Param::eposition},
        {"sublimination", Param::sublimination}, {"a_1", Param::a_1},
        {"a_2", Param::a_2}, {"e_1", Param::e_1},
        {"e_2", Param::e_2}, {"d", Param::d},
        {"N_c", Param::N_c}, {"gamma", Param::gamma},
        {"beta_c", Param::beta_c}, {"beta_r", Param::beta_r},
        {"delta1", Param::delta1}, {"delta2", Param::delta2},
        {"zeta", Param::zeta}, {"rain_gfak", Param::rain_gfak},
        {"cloud_k_au", Param::cloud_k_au}, {"cloud_k_sc", Param::cloud_k_sc},
        {"kc_autocon", Param::kc_autocon}, {"inv_z", Param::inv_z},
        {"dw", Param::dw}, {"q_crit_i", Param::q_crit_i},
        {"D_crit_i", Param::D_crit_i}, {"D_conv_i", Param::D_conv_i},
        {"q_crit_r", Param::q_crit_r}, {"D_crit_r", Param::D_crit_r},
        {"q_crit_fr", Param::q_crit_fr}, {"D_coll_c", Param::D_coll_c},
        {"q_crit", Param::q_crit}, {"D_conv_sg", Param::D_conv_sg},
        {"D_conv_ig", Param::D_conv_ig}, {"x_conv", Param::x_conv},
        {"parcel_height", Param::parcel_height}, {"alpha_spacefilling", Param::alpha_spacefilling},
        {"T_nuc", Param::T_nuc}, {"T_freeze", Param::T_freeze},
        {"T_f", Param::T_f}, {"D_eq", Param::D_eq},
        {"rho_w", Param::rho_w}, {"rho_0", Param::rho_0},
        {"rho_vel", Param::rho_vel}, {"rho_vel_c", Param::rho_vel_c},
        {"rho_ice", Param::rho_ice}, {"M_w", Param::M_w},
        {"M_a", Param::M_a}, {"R_universal", Param::R_universal},
        {"Epsilon", Param::Epsilon}, {"gravity_acc", Param::gravity_acc},
        {"R_a", Param::R_a}, {"R_v", Param::R_v},
        {"a_v", Param::a_v}, {"b_v", Param::b_v},
        {"a_prime", Param::a_prime}, {"b_prime", Param::b_prime},
        {"c_prime", Param::c_prime}, {"K_T", Param::K_T},
        {"L_wd", Param::L_wd}, {"L_ed", Param::L_ed},
        {"D_v", Param::D_v}, {"ecoll_min", Param::ecoll_min},
        {"ecoll_gg", Param::ecoll_gg}, {"ecoll_gg_wet", Param::ecoll_gg_wet},
        {"kin_visc_air", Param::kin_visc_air}, {"C_mult", Param::C_mult},
        {"T_mult_min", Param::T_mult_min}, {"T_mult_max", Param::T_mult_max},
        {"T_mult_opt", Param::T_mult_opt}, {"const0", Param::const0},
        {"const3", Param::const3}, {"const4", Param::const4},
        {"const5", Param::const5}, {"D_rainfrz_ig", Param::D_rainfrz_ig},
        {"dv0", Param::dv0}, {"p_sat_melt", Param::p_sat_melt},
        {"cp", Param::cp}, {"k_b", Param::k_b},
        {"a_HET", Param::a_HET}, {"b_HET", Param::b_HET},
        {"N_sc", Param::N_sc}, {"n_f", Param::n_f},
        {"N_avo", Param::N_avo}, {"na_dust", Param::na_dust},
        {"na_soot", Param::na_soot}, {"na_orga", Param::na_orga},
        {"ni_het_max", Param::ni_het_max}, {"ni_hom_max", Param::ni_hom_max},
        {"a_dep", Param::a_dep}, {"b_dep", Param::b_dep},
        {"c_dep", Param::c_dep}, {"d_dep", Param::d_dep},
        {"nim_imm", Param::nim_imm}, {"nin_dep", Param::nin_dep},
        {"alf_imm", Param::alf_imm}, {"bet_dep", Param::bet_dep},
        {"bet_imm", Param::bet_imm}, {"r_const", Param::r_const},
        {"r1_const", Param::r1_const}, {"cv", Param::cv},
        {"p_sat_const_a", Param::p_sat_const_a}, {"p_sat_ice_const_a", Param::p_sat_ice_const_a},
        {"p_sat_const_b", Param::p_sat_const_b}, {"p_sat_ice_const_b", Param::p_sat_ice_const_b},
        {"p_sat_low_temp", Param::p_sat_low_temp}, {"T_sat_low_temp", Param::T_sat_low_temp},
        {"alpha_depo", Param::alpha_depo}, {"r_0", Param::r_0},
        {"k_1_conv", Param::k_1_conv}, {"k_2_conv", Param::k_2_conv},
        {"k_1_accr", Param::k_1_accr}, {"k_r", Param::k_r},
        {"a_ccn_1", Param::a_ccn_1}, {"a_ccn_2", Param::a_ccn_2},
        {"a_ccn_3", Param::a_ccn_3}, {"a_ccn_4", Param::a_ccn_4},
        {"b_ccn_1", Param::b_ccn_1}, {"b_ccn_2", Param::b_ccn_2},
        {"b_ccn_3", Param::b_ccn_3}, {"b_ccn_4", Param::b_ccn_4},
        {"c_ccn_1", Param::c_ccn_1}, {"c_ccn_2", Param::c_ccn_2},
        {"c_ccn_3", Param::c_ccn_3}, {"c_ccn_4", Param::c_ccn_4},
        {"d_ccn_1", Param::d_ccn_1}, {"d_ccn_2", Param::d_ccn_2},
        {"d_ccn_3", Param::d_ccn_3}, {"d_ccn_4", Param::d_ccn_4},
        {"rain_a_geo", Param::rain_a_geo}, {"rain_b_geo", Param::rain_b_geo},
        {"rain_min_x", Param::rain_min_x}, {"rain_min_x_act", Param::rain_min_x_act},
        {"rain_min_x_nuc_homo", Param::rain_min_x_nuc_homo}, {"rain_min_x_nuc_hetero", Param::rain_min_x_nuc_hetero},
        {"rain_min_x_melt", Param::rain_min_x_melt}, {"rain_min_x_evap", Param::rain_min_x_evap},
        {"rain_min_x_freezing", Param::rain_min_x_freezing}, {"rain_min_x_depo", Param::rain_min_x_depo},
        {"rain_min_x_collision", Param::rain_min_x_collision}, {"rain_min_x_collection", Param::rain_min_x_collection},
        {"rain_min_x_conversion", Param::rain_min_x_conversion}, {"rain_min_x_sedimentation", Param::rain_min_x_sedimentation},
        {"rain_min_x_riming", Param::rain_min_x_riming}, {"rain_max_x", Param::rain_max_x},
        {"rain_sc_theta_q", Param::rain_sc_theta_q}, {"rain_sc_delta_q", Param::rain_sc_delta_q},
        {"rain_sc_theta_n", Param::rain_sc_theta_n}, {"rain_sc_delta_n", Param::rain_sc_delta_n},
        {"rain_s_vel", Param::rain_s_vel}, {"rain_a_vel", Param::rain_a_vel},
        {"rain_b_vel", Param::rain_b_vel}, {"rain_rho_v", Param::rain_rho_v},
        {"rain_c_z", Param::rain_c_z}, {"rain_sc_coll_n", Param::rain_sc_coll_n},
        {"rain_cmu0", Param::rain_cmu0}, {"rain_cmu1", Param::rain_cmu1},
        {"rain_cmu2", Param::rain_cmu2}, {"rain_cmu3", Param::rain_cmu3},
        {"rain_cmu4", Param::rain_cmu4}, {"rain_cmu5", Param::rain_cmu5},
        {"rain_alpha", Param::rain_alpha}, {"rain_beta", Param::rain_beta},
        {"rain_gamma", Param::rain_gamma}, {"rain_nu", Param::rain_nu},
        {"rain_g1", Param::rain_g1}, {"rain_g2", Param::rain_g2},
        {"rain_mu", Param::rain_mu}, {"rain_nm1", Param::rain_nm1},
        {"rain_nm2", Param::rain_nm2}, {"rain_nm3", Param::rain_nm3},
        {"rain_q_crit_c", Param::rain_q_crit_c}, {"rain_d_crit_c", Param::rain_d_crit_c},
        {"rain_ecoll_c", Param::rain_ecoll_c}, {"rain_cap", Param::rain_cap},
        {"rain_a_ven", Param::rain_a_ven}, {"rain_b_ven", Param::rain_b_ven},
        {"rain_c_s", Param::rain_c_s}, {"rain_a_f", Param::rain_a_f},
        {"rain_b_f", Param::rain_b_f}, {"rain_alfa_n", Param::rain_alfa_n},
        {"rain_alfa_q", Param::rain_alfa_q}, {"rain_lambda", Param::rain_lambda},
        {"rain_vsedi_min", Param::rain_vsedi_min}, {"rain_vsedi_max", Param::rain_vsedi_max},
        {"cloud_a_geo", Param::cloud_a_geo}, {"cloud_b_geo", Param::cloud_b_geo},
        {"cloud_min_x", Param::cloud_min_x}, {"cloud_min_x_act", Param::cloud_min_x_act},
        {"cloud_min_x_nuc_homo", Param::cloud_min_x_nuc_homo}, {"cloud_min_x_nuc_hetero", Param::cloud_min_x_nuc_hetero},
        {"cloud_min_x_melt", Param::cloud_min_x_melt}, {"cloud_min_x_evap", Param::cloud_min_x_evap},
        {"cloud_min_x_freezing", Param::cloud_min_x_freezing}, {"cloud_min_x_depo", Param::cloud_min_x_depo},
        {"cloud_min_x_collision", Param::cloud_min_x_collision}, {"cloud_min_x_collection", Param::cloud_min_x_collection},
        {"cloud_min_x_conversion", Param::cloud_min_x_conversion}, {"cloud_min_x_sedimentation", Param::cloud_min_x_sedimentation},
        {"cloud_min_x_riming", Param::cloud_min_x_riming}, {"cloud_max_x", Param::cloud_max_x},
        {"cloud_sc_theta_q", Param::cloud_sc_theta_q}, {"cloud_sc_delta_q", Param::cloud_sc_delta_q},
        {"cloud_sc_theta_n", Param::cloud_sc_theta_n}, {"cloud_sc_delta_n", Param::cloud_sc_delta_n},
        {"cloud_s_vel", Param::cloud_s_vel}, {"cloud_a_vel", Param::cloud_a_vel},
        {"cloud_b_vel", Param::cloud_b_vel}, {"cloud_rho_v", Param::cloud_rho_v},
        {"cloud_c_z", Param::cloud_c_z}, {"cloud_sc_coll_n", Param::cloud_sc_coll_n},
        {"cloud_cmu0", Param::cloud_cmu0}, {"cloud_cmu1", Param::cloud_cmu1},
        {"cloud_cmu2", Param::cloud_cmu2}, {"cloud_cmu3", Param::cloud_cmu3},
        {"cloud_cmu4", Param::cloud_cmu4}, {"cloud_cmu5", Param::cloud_cmu5},
        {"cloud_alpha", Param::cloud_alpha}, {"cloud_beta", Param::cloud_beta},
        {"cloud_gamma", Param::cloud_gamma}, {"cloud_nu", Param::cloud_nu},
        {"cloud_g1", Param::cloud_g1}, {"cloud_g2", Param::cloud_g2},
        {"cloud_mu", Param::cloud_mu}, {"cloud_nm1", Param::cloud_nm1},
        {"cloud_nm2", Param::cloud_nm2}, {"cloud_nm3", Param::cloud_nm3},
        {"cloud_q_crit_c", Param::cloud_q_crit_c}, {"cloud_d_crit_c", Param::cloud_d_crit_c},
        {"cloud_ecoll_c", Param::cloud_ecoll_c}, {"cloud_cap", Param::cloud_cap},
        {"cloud_a_ven", Param::cloud_a_ven}, {"cloud_b_ven", Param::cloud_b_ven},
        {"cloud_c_s", Param::cloud_c_s}, {"cloud_a_f", Param::cloud_a_f},
        {"cloud_b_f", Param::cloud_b_f}, {"cloud_alfa_n", Param::cloud_alfa_n},
        {"cloud_alfa_q", Param::cloud_alfa_q}, {"cloud_lambda", Param::cloud_lambda},
        {"cloud_vsedi_min", Param::cloud_vsedi_min}, {"cloud_vsedi_max", Param::cloud_vsedi_max},
        {"graupel_a_geo", Param::graupel_a_geo}, {"graupel_b_geo", Param::graupel_b_geo},
        {"graupel_min_x", Param::graupel_min_x}, {"graupel_min_x_act", Param::graupel_min_x_act},
        {"graupel_min_x_nuc_homo", Param::graupel_min_x_nuc_homo}, {"graupel_min_x_nuc_hetero", Param::graupel_min_x_nuc_hetero},
        {"graupel_min_x_melt", Param::graupel_min_x_melt}, {"graupel_min_x_evap", Param::graupel_min_x_evap},
        {"graupel_min_x_freezing", Param::graupel_min_x_freezing}, {"graupel_min_x_depo", Param::graupel_min_x_depo},
        {"graupel_min_x_collision", Param::graupel_min_x_collision}, {"graupel_min_x_collection", Param::graupel_min_x_collection},
        {"graupel_min_x_conversion", Param::graupel_min_x_conversion}, {"graupel_min_x_sedimentation", Param::graupel_min_x_sedimentation},
        {"graupel_min_x_riming", Param::graupel_min_x_riming}, {"graupel_max_x", Param::graupel_max_x},
        {"graupel_sc_theta_q", Param::graupel_sc_theta_q}, {"graupel_sc_delta_q", Param::graupel_sc_delta_q},
        {"graupel_sc_theta_n", Param::graupel_sc_theta_n}, {"graupel_sc_delta_n", Param::graupel_sc_delta_n},
        {"graupel_s_vel", Param::graupel_s_vel}, {"graupel_a_vel", Param::graupel_a_vel},
        {"graupel_b_vel", Param::graupel_b_vel}, {"graupel_rho_v", Param::graupel_rho_v},
        {"graupel_c_z", Param::graupel_c_z}, {"graupel_sc_coll_n", Param::graupel_sc_coll_n},
        {"graupel_cmu0", Param::graupel_cmu0}, {"graupel_cmu1", Param::graupel_cmu1},
        {"graupel_cmu2", Param::graupel_cmu2}, {"graupel_cmu3", Param::graupel_cmu3},
        {"graupel_cmu4", Param::graupel_cmu4}, {"graupel_cmu5", Param::graupel_cmu5},
        {"graupel_alpha", Param::graupel_alpha}, {"graupel_beta", Param::graupel_beta},
        {"graupel_gamma", Param::graupel_gamma}, {"graupel_nu", Param::graupel_nu},
        {"graupel_g1", Param::graupel_g1}, {"graupel_g2", Param::graupel_g2},
        {"graupel_mu", Param::graupel_mu}, {"graupel_nm1", Param::graupel_nm1},
        {"graupel_nm2", Param::graupel_nm2}, {"graupel_nm3", Param::graupel_nm3},
        {"graupel_q_crit_c", Param::graupel_q_crit_c}, {"graupel_d_crit_c", Param::graupel_d_crit_c},
        {"graupel_ecoll_c", Param::graupel_ecoll_c}, {"graupel_cap", Param::graupel_cap},
        {"graupel_a_ven", Param::graupel_a_ven}, {"graupel_b_ven", Param::graupel_b_ven},
        {"graupel_c_s", Param::graupel_c_s}, {"graupel_a_f", Param::graupel_a_f},
        {"graupel_b_f", Param::graupel_b_f}, {"graupel_alfa_n", Param::graupel_alfa_n},
        {"graupel_alfa_q", Param::graupel_alfa_q}, {"graupel_lambda", Param::graupel_lambda},
        {"graupel_vsedi_min", Param::graupel_vsedi_min}, {"graupel_vsedi_max", Param::graupel_vsedi_max},
        {"hail_a_geo", Param::hail_a_geo}, {"hail_b_geo", Param::hail_b_geo},
        {"hail_min_x", Param::hail_min_x}, {"hail_min_x_act", Param::hail_min_x_act},
        {"hail_min_x_nuc_homo", Param::hail_min_x_nuc_homo}, {"hail_min_x_nuc_hetero", Param::hail_min_x_nuc_hetero},
        {"hail_min_x_melt", Param::hail_min_x_melt}, {"hail_min_x_evap", Param::hail_min_x_evap},
        {"hail_min_x_freezing", Param::hail_min_x_freezing}, {"hail_min_x_depo", Param::hail_min_x_depo},
        {"hail_min_x_collision", Param::hail_min_x_collision}, {"hail_min_x_collection", Param::hail_min_x_collection},
        {"hail_min_x_conversion", Param::hail_min_x_conversion}, {"hail_min_x_sedimentation", Param::hail_min_x_sedimentation},
        {"hail_min_x_riming", Param::hail_min_x_riming}, {"hail_max_x", Param::hail_max_x},
        {"hail_sc_theta_q", Param::hail_sc_theta_q}, {"hail_sc_delta_q", Param::hail_sc_delta_q},
        {"hail_sc_theta_n", Param::hail_sc_theta_n}, {"hail_sc_delta_n", Param::hail_sc_delta_n},
        {"hail_s_vel", Param::hail_s_vel}, {"hail_a_vel", Param::hail_a_vel},
        {"hail_b_vel", Param::hail_b_vel}, {"hail_rho_v", Param::hail_rho_v},
        {"hail_c_z", Param::hail_c_z}, {"hail_sc_coll_n", Param::hail_sc_coll_n},
        {"hail_cmu0", Param::hail_cmu0}, {"hail_cmu1", Param::hail_cmu1},
        {"hail_cmu2", Param::hail_cmu2}, {"hail_cmu3", Param::hail_cmu3},
        {"hail_cmu4", Param::hail_cmu4}, {"hail_cmu5", Param::hail_cmu5},
        {"hail_alpha", Param::hail_alpha}, {"hail_beta", Param::hail_beta},
        {"hail_gamma", Param::hail_gamma}, {"hail_nu", Param::hail_nu},
        {"hail_g1", Param::hail_g1}, {"hail_g2", Param::hail_g2},
        {"hail_mu", Param::hail_mu}, {"hail_nm1", Param::hail_nm1},
        {"hail_nm2", Param::hail_nm2}, {"hail_nm3", Param::hail_nm3},
        {"hail_q_crit_c", Param::hail_q_crit_c}, {"hail_d_crit_c", Param::hail_d_crit_c},
        {"hail_ecoll_c", Param::hail_ecoll_c}, {"hail_cap", Param::hail_cap},
        {"hail_a_ven", Param::hail_a_ven}, {"hail_b_ven", Param::hail_b_ven},
        {"hail_c_s", Param::hail_c_s}, {"hail_a_f", Param::hail_a_f},
        {"hail_b_f", Param::hail_b_f}, {"hail_alfa_n", Param::hail_alfa_n},
        {"hail_alfa_q", Param::hail_alfa_q}, {"hail_lambda", Param::hail_lambda},
        {"hail_vsedi_min", Param::hail_vsedi_min}, {"hail_vsedi_max", Param::hail_vsedi_max},
        {"ice_a_geo", Param::ice_a_geo}, {"ice_b_geo", Param::ice_b_geo},
        {"ice_min_x", Param::ice_min_x}, {"ice_min_x_act", Param::ice_min_x_act},
        {"ice_min_x_nuc_homo", Param::ice_min_x_nuc_homo}, {"ice_min_x_nuc_hetero", Param::ice_min_x_nuc_hetero},
        {"ice_min_x_melt", Param::ice_min_x_melt}, {"ice_min_x_evap", Param::ice_min_x_evap},
        {"ice_min_x_freezing", Param::ice_min_x_freezing}, {"ice_min_x_depo", Param::ice_min_x_depo},
        {"ice_min_x_collision", Param::ice_min_x_collision}, {"ice_min_x_collection", Param::ice_min_x_collection},
        {"ice_min_x_conversion", Param::ice_min_x_conversion}, {"ice_min_x_sedimentation", Param::ice_min_x_sedimentation},
        {"ice_min_x_riming", Param::ice_min_x_riming}, {"ice_max_x", Param::ice_max_x},
        {"ice_sc_theta_q", Param::ice_sc_theta_q}, {"ice_sc_delta_q", Param::ice_sc_delta_q},
        {"ice_sc_theta_n", Param::ice_sc_theta_n}, {"ice_sc_delta_n", Param::ice_sc_delta_n},
        {"ice_s_vel", Param::ice_s_vel}, {"ice_a_vel", Param::ice_a_vel},
        {"ice_b_vel", Param::ice_b_vel}, {"ice_rho_v", Param::ice_rho_v},
        {"ice_c_z", Param::ice_c_z}, {"ice_sc_coll_n", Param::ice_sc_coll_n},
        {"ice_cmu0", Param::ice_cmu0}, {"ice_cmu1", Param::ice_cmu1},
        {"ice_cmu2", Param::ice_cmu2}, {"ice_cmu3", Param::ice_cmu3},
        {"ice_cmu4", Param::ice_cmu4}, {"ice_cmu5", Param::ice_cmu5},
        {"ice_alpha", Param::ice_alpha}, {"ice_beta", Param::ice_beta},
        {"ice_gamma", Param::ice_gamma}, {"ice_nu", Param::ice_nu},
        {"ice_g1", Param::ice_g1}, {"ice_g2", Param::ice_g2},
        {"ice_mu", Param::ice_mu}, {"ice_nm1", Param::ice_nm1},
        {"ice_nm2", Param::ice_nm2}, {"ice_nm3", Param::ice_nm3},
        {"ice_q_crit_c", Param::ice_q_crit_c}, {"ice_d_crit_c", Param::ice_d_crit_c},
        {"ice_ecoll_c", Param::ice_ecoll_c}, {"ice_cap", Param::ice_cap},
        {"ice_a_ven", Param::ice_a_ven}, {"ice_b_ven", Param::ice_b_ven},
        {"ice_c_s", Param::ice_c_s}, {"ice_a_f", Param::ice_a_f},
        {"ice_b_f", Param::ice_b_f}, {"ice_alfa_n", Param::ice_alfa_n},
        {"ice_alfa_q", Param::ice_alfa_q}, {"ice_lambda", Param::ice_lambda},
        {"ice_vsedi_min", Param::ice_vsedi_min}, {"ice_vsedi_max", Param::ice_vsedi_max},
        {"snow_a_geo", Param::snow_a_geo}, {"snow_b_geo", Param::snow_b_geo},
        {"snow_min_x", Param::snow_min_x}, {"snow_min_x_act", Param::snow_min_x_act},
        {"snow_min_x_nuc_homo", Param::snow_min_x_nuc_homo}, {"snow_min_x_nuc_hetero", Param::snow_min_x_nuc_hetero},
        {"snow_min_x_melt", Param::snow_min_x_melt}, {"snow_min_x_evap", Param::snow_min_x_evap},
        {"snow_min_x_freezing", Param::snow_min_x_freezing}, {"snow_min_x_depo", Param::snow_min_x_depo},
        {"snow_min_x_collision", Param::snow_min_x_collision}, {"snow_min_x_collection", Param::snow_min_x_collection},
        {"snow_min_x_conversion", Param::snow_min_x_conversion}, {"snow_min_x_sedimentation", Param::snow_min_x_sedimentation},
        {"snow_min_x_riming", Param::snow_min_x_riming}, {"snow_max_x", Param::snow_max_x},
        {"snow_sc_theta_q", Param::snow_sc_theta_q}, {"snow_sc_delta_q", Param::snow_sc_delta_q},
        {"snow_sc_theta_n", Param::snow_sc_theta_n}, {"snow_sc_delta_n", Param::snow_sc_delta_n},
        {"snow_s_vel", Param::snow_s_vel}, {"snow_a_vel", Param::snow_a_vel},
        {"snow_b_vel", Param::snow_b_vel}, {"snow_rho_v", Param::snow_rho_v},
        {"snow_c_z", Param::snow_c_z}, {"snow_sc_coll_n", Param::snow_sc_coll_n},
        {"snow_cmu0", Param::snow_cmu0}, {"snow_cmu1", Param::snow_cmu1},
        {"snow_cmu2", Param::snow_cmu2}, {"snow_cmu3", Param::snow_cmu3},
        {"snow_cmu4", Param::snow_cmu4}, {"snow_cmu5", Param::snow_cmu5},
        {"snow_alpha", Param::snow_alpha}, {"snow_beta", Param::snow_beta},
        {"snow_gamma", Param::snow_gamma}, {"snow_nu", Param::snow_nu},
        {"snow_g1", Param::snow_g1}, {"snow_g2", Param::snow_g2},
        {"snow_mu", Param::snow_mu}, {"snow_nm1", Param::snow_nm1},
        {"snow_nm2", Param::snow_nm2}, {"snow_nm3", Param::snow_nm3},
        {"snow_q_crit_c", Param::snow_q_crit_c}, {"snow_d_crit_c", Param::snow_d_crit_c},
        {"snow_ecoll_c", Param::snow_ecoll_c}, {"snow_cap", Param::snow_cap},
        {"snow_a_ven", Param::snow_a_ven}, {"snow_b_ven", Param::snow_b_ven},
        {"snow_c_s", Param::snow_c_s}, {"snow_a_f", Param::snow_a_f},
        {"snow_b_f", Param::snow_b_f}, {"snow_alfa_n", Param::snow_alfa_n},
        {"snow_alfa_q", Param::snow_alfa_q}, {"snow_lambda", Param::snow_lambda},
        {"snow_vsedi_min", Param::snow_vsedi_min}, {"snow_vsedi_max", Param::snow_vsedi_max}
    };

    enum OutParam {
        vapor, cloud, rain, ice, graupel, hail, snow
    };
    std::unordered_map<std::string, OutParam> const table_out_param = {
        {"snow", OutParam::snow}, {"vapor", OutParam::vapor},
        {"cloud", OutParam::cloud}, {"rain", OutParam::rain},
        {"ice", OutParam::ice}, {"graupel", OutParam::graupel},
        {"hail", OutParam::hail}
    };


    void add_param(param_t &param)
    {
        params.push_back(param);
    }

    void add_value(double v)
    {
        value = v;
    }

    void add_value_name(std::string n)
    {
        auto it = table_param.find(n);
        if(it != table_param.end())
        {
            value_name = it->second;
        } else
        {
            err = VALUE_NAME_CONFIG_ERR;
        }
    }

    void add_method(std::string m)
    {
        auto it = table_method.find(m);
        if(it != table_method.end())
        {
            method = it->second;
        } else
        {
            err = METHOD_CONFIG_ERR;
        }
    }

    void add_counter(int c)
    {
        if(c > 0)
            n_segments = c;
        else
            err = N_SEGMENTS_CONFIG_ERR;
    }

    void add_out_param(std::string p)
    {
        auto it = table_out_param.find(p);
        if(it != table_out_param.end())
        {
            out_param = it->second;
        } else
        {
            err = OUTPARAM_CONFIG_ERR;
        }
    }

    void add_amount(int a)
    {
        if(a > 1)
            n_members = a;
        else
            err = N_MEMBERS_CONFIG_ERR;
    }

    int check()
    {
        if(n_members == 1)
            err = N_MEMBERS_CONFIG_ERR;
        else if(n_segments < 1)
            err = N_SEGMENTS_CONFIG_ERR;
        else if(value_name == -1)
            err = VALUE_NAME_CONFIG_ERR;
        else if(method == -1)
            err = METHOD_CONFIG_ERR;
        else if(out_param == -1)
            err = OUTPARAM_CONFIG_ERR;
        switch(err)
        {
            case VALUE_NAME_CONFIG_ERR:
                std::cout << "Error in config file:\n"
                          << "The name of the parameter used to determine "
                          << "when to start an ensemble is not defined or "
                          << "does not exist.\n"
                          << "Skipping this segment\n";
                n_segments = 0;
                return err
            case METHOD_CONFIG_ERR
                std::cout << "Error in config file:\n"
                          << "The name of the method used to determine "
                          << "when to start an ensemble is not defined or "
                          << "does not exist.\n"
                          << "Skipping this segment\n";
                n_segments = 0;
                return err
            case OUTPARAM_CONFIG_ERR:
                std::cout << "Error in config file:\n"
                          << "The name of the output parameter used to determine "
                          << "when to start an ensemble is not defined or "
                          << "does not exist.\n"
                          << "Skipping this segment\n";
                n_segments = 0;
                return err
            case N_MEMBERS_CONFIG_ERR:
                std::cout << "Error in config file:\n"
                          << "The number of members for an ensemble must be "
                          << "2 or higher. One simulation always runs without "
                          << "perturbed parameters for comparison.\n"
                          << "Skipping this segment\n";
                n_segments = 0;
                return err;
            case N_SEGMENTS_CONFIG_ERR:
                std::cout << "Error in config file:\n"
                          << "The number of segments must be at least 1 in "
                          << "order to start at least 1 ensemble.\n"
                          << "Skipping this segment\n";
                n_segments = 0;
                return err;
        }
        if(err != 0)
            n_segments = 0;
        return err;
    }

    template<class float_t>
    void perturb_check(model_constants_t &cc, std::vector<float_t> &gradients)
    {
        // TODO: Add derivatives and current state
        if(n_segments == 0)
            return;
        switch(method)
        {
            case impact_change:
                int idx = value_name - num_comp;


                break;
            case sign_flip:
                int idx = value_name - num_comp;
                if(byte == 0)
                {
                    // set sign
                    if(gradients[idx] == 0)
                        break;
                    byte = (gradients[idx] > 0) ? 1 : 2;
                } else
                {
                    if(byte == 1 && gradients[idx] < 0)
                    {
                        // Perturb parameters
                        // idea: Run a new script that allocates
                        // the necessary stuff using
                        // stdlib.h
                        // system("next_job.job")
                        // Other idea: create a checkpoint that
                        // can be used as start for more processes
                        n_segments--;
                    }else if(byte == 0 && gradients[idx] > 0)
                    {
                        // Perturb parameters
                        n_segments--;
                    }
                }
                break;
            case value_method:
                if(out_param != -1)
                {
                    int idx = value_name - num_comp;
                } else
                {

                }
                break;
        }
    }
};

/** @} */ // end of group types
