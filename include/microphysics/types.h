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

    codi::RealReverse a_geo;
    codi::RealReverse b_geo;
    /**
     * Minimum size of particle for mean meass calculation.
     */
    codi::RealReverse min_x;
    /**
     * Minimum size of particle for activation.
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_act;
    /**
     * Minimum size of particle for nucleation.
     * *Should* be the same as min_x but is used to distinguish the
     * influence of those processes.
     */
    codi::RealReverse min_x_nuc_homo;

    codi::RealReverse min_x_nuc_hetero;
    codi::RealReverse min_x_melt;
    codi::RealReverse min_x_evap;
    codi::RealReverse min_x_freezing;
    codi::RealReverse min_x_depo;
    codi::RealReverse min_x_collision;
    codi::RealReverse min_x_collection;
    codi::RealReverse min_x_conversion;
    codi::RealReverse min_x_sedimentation;
    codi::RealReverse min_x_riming;

    codi::RealReverse max_x;
    codi::RealReverse sc_theta_q;   /*!< For snow collision. */
    codi::RealReverse sc_delta_q;
    codi::RealReverse sc_theta_n;
    codi::RealReverse sc_delta_n;
    codi::RealReverse s_vel;        /*!< Also known as sigma_vel in ICON. */
    codi::RealReverse a_vel;
    codi::RealReverse b_vel;
    codi::RealReverse rho_v;
    codi::RealReverse c_z;
    codi::RealReverse sc_coll_n;    /*!< Also known as e_coll in ICON. */
    codi::RealReverse cmu0, cmu1, cmu2, cmu3, cmu4, cmu5, alpha, beta, gamma;
    codi::RealReverse nu;

    codi::RealReverse g1, g2, mu, nm1, nm2, nm3;

    codi::RealReverse q_crit_c; /*!<  Riming parameter. */
    codi::RealReverse d_crit_c; /*!<  Riming parameter. */
    codi::RealReverse ecoll_c;  /*!<  Riming parameter. */
    /**
     * Riming parameter. Coefficient for capacity of particle.
     */
    codi::RealReverse cap;
    codi::RealReverse a_ven;    /*!<  Riming parameter. */
    codi::RealReverse b_ven;    /*!<  Riming parameter. */

    codi::RealReverse c_s;  /*!< Also known as c_i in ICON.*/
    codi::RealReverse a_f;
    codi::RealReverse b_f;

    codi::RealReverse alfa_n;       /*!<  Bulk sedimentation velocity parameter. */
    codi::RealReverse alfa_q;       /*!<  Bulk sedimentation velocity parameter. */
    codi::RealReverse lambda;       /*!<  Bulk sedimentation velocity parameter. */
    codi::RealReverse vsedi_min;    /*!<  Bulk sedimentation velocity parameter. */
    codi::RealReverse vsedi_max;    /*!<  Bulk sedimentation velocity parameter. */

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

  codi::RealReverse rho_a_prime; /*!< Density of dry air */

  codi::RealReverse Nc_prime; /*!< Number concentration of cloud droplets */

  codi::RealReverse a1_prime; /*!< Dimensional coefficients */
  codi::RealReverse a2_prime; /*!< Dimensional coefficients */
  codi::RealReverse e1_prime; /*!< Dimensional coefficients */
  codi::RealReverse e2_prime; /*!< Dimensional coefficients */
  codi::RealReverse B_prime;  /*!< Dimensional coefficients */
  codi::RealReverse d_prime;  /*!< Dimensional coefficients */

  codi::RealReverse dw; /*!< Change in buoancy */

  codi::RealReverse gamma;  /*!< Exponents */
  codi::RealReverse betac;  /*!< Exponents */
  codi::RealReverse betar;  /*!< Exponents */
  codi::RealReverse delta1; /*!< Exponents */
  codi::RealReverse delta2; /*!< Exponents */
  codi::RealReverse zeta;   /*!< Exponents */

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
  double dt_traject_prime;  /*!< Timestep size in seconds of the trajectory from the netCDF file. */
  double dt_traject;        /*!< Timestep size of the trajectory from the netCDF file. */
  uint64_t num_steps;       /*!< Number of timesteps to read from the netCDF file. */
  /**
   * Number of timesteps to simulate between each timestep of the netCDF file.
   */
  uint64_t num_sub_steps;
  /**
   * Number of simulation steps before a snapshot shall be stored on disk.
   */
  uint64_t snapshot_index;

  /**
   * Number of simulation steps before a snapshot shall stored on disk.
   */
  uint64_t write_index;

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
  codi::RealReverse rain_gfak = -1.0; /*!<  See mo_2mom_mcrph_main.f90 line 830 following of ICON. */
  codi::RealReverse cloud_k_au; /*!< Parameter for autoconversion Seifert & Beheng. */
  codi::RealReverse cloud_k_sc; /*!< Parameter for autoconversion Seifert & Beheng. */

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
  codi::RealReverse inv_z = 1.0/150.0;

  const double nar = 0.22;      /*!< Constants for the IFS model. */
  const double nbr = 2.2;       /*!< Constants for the IFS model. */
  const double ar = M_PI / 6.0; /*!< Constants for the IFS model. */
  const double br = 3.0;        /*!< Constants for the IFS model. */
  const double cr = 386.8;      /*!< Constants for the IFS model. */
  const double dr = 0.67;       /*!< Constants for the IFS model. */
  const double Sc = 0.6;        /*!< Constants for the IFS model. */
  const double mu = 16.0e-6;    /*!< Constants for the IFS model. */
  const double rho0 = 1.0;      /*!< Constants for the IFS model. Why not 1.225? */

  const double alphar = 1.0/(br + 1.0 - nbr);   /*!< Constants for the IFS model. */
  const double epsilonr = 0.5*dr + 2.5 - nbr;   /*!< Constants for the IFS model. */

  double scaling_fact; /*!< Scaling factor. */

};


/**
 * Structure to collect all nc parameters.
 */
struct nc_parameters_t{

    uint32_t n_trajectories = 30; /*!< Number of trajectories in the netCDF file. */
    uint32_t n_timesteps = 7922; /*!< Number of timesteps in the netCDF file. */
    std::vector<double> w, z, lat, lon;
    double  t, p, time_rel,
            qc, qr, qi, qs, qg, qv, S, dw, dlat, dlon,
            QIin, QSin, QRin, QGin, QIout, QSout, QRout, QGout,
            NIin, NSin, NRin, NGin, NIout, NSout, NRout, NGout,
            Nc, Nr, Ni, Ns, Ng;
    bool ascent_flag, conv_400, conv_600, slan_400, slan_600, dp2h;
    NcVar   lat_var, lon_var, z_var, t_var, p_var, w_var, time_rel_var,
            qc_var, qr_var, qi_var, qs_var, qg_var, qv_var, S_var,
            QIin_var, QSin_var, QRin_var, QGin_var, QIout_var, QSout_var,
            QRout_var, QGout_var, ascent_flag_var,
            NIin_var, NSin_var, NRin_var, NGin_var, NIout_var, NSout_var,
            NRout_var, NGout_var, dp2h_var,
            Nc_var, Nr_var, Ni_var, Ns_var, Ng_var,
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

  int start_over_flag; /*!< Reload data from trajectory every few seconds? */
  char* start_over_string;

  int fixed_iteration_flag; /*!< Fix p, T, w during simulation? */
  char* fixed_iteration_string;

  int auto_type_flag; /*!< Particle type specified? */
  char* auto_type_string;

  int traj_flag; /*!< Trajectory to use specified? */
  char* traj_string;

  int write_flag; /*!< Snapshot is flushed every x iterations. */
  char* write_string;
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
  int snapshot_index; /*!< Save a snapshot every snapshot_index iteration. */
  /**
   * Number of timesteps for the simulation between two
   * datapoints from the netCDF file.
   */
  uint64_t num_sub_steps;
  std::string OUTPUT_FILENAME; /*!< Filename for output. */

  std::string INPUT_FILENAME; /*!< Filename for input netCDF file. */

  bool start_over; /*!< Start over at new timestep of trajectory? */
  bool fixed_iteration; /*!< Fix temperature, pressure and ascension at every iteration? */

  double scaling_fact; /*!< Scaling factor. */

  uint32_t auto_type; /*!< Particle type. */
  uint32_t traj; /*!< Trajectory index to load from the netCDF file. */
  uint32_t write_index; /*!< Write stringstream every x iterations to disk. */
};

/** @} */ // end of group types
