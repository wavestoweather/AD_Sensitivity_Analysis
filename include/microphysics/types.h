#pragma once
#include "codi.hpp"
#include <netcdf>
#include <boost/math/special_functions/gamma.hpp>

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

/** @} */ // end of group types
