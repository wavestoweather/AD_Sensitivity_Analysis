#pragma once
#include <boost/property_tree/ptree.hpp>
#include "codi.hpp"
#include <netcdf>
#include "ncType.h"
#include <random>
#include <unordered_map>
#include <stdlib.h>
#include <boost/math/special_functions/gamma.hpp>
#include "include/misc/error.h"
#include "constants.h"

using namespace netCDF;
namespace pt = boost::property_tree;


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

    std::vector<codi::RealReverse> constants;
    std::vector<uint32_t> perturbed_idx;

    particle_model_constants_t()
    {
        constants.resize(static_cast<int>(Particle_cons_idx::n_items));
        // This makes debugging easier, so pleaase leave it.
        std::fill(constants.begin(), constants.end(), 0);
    }

    /**
     * Register the model parameters on the tape for codi::RealReverse.
     *
     * @param tape Tape where the parameters should be registered on.
     */
    void register_input(codi::RealReverse::TapeType &tape)
    {
        for(auto &c: this->constants)
            tape.registerInput(c);
    }

    /**
     * Get the gradients of all its members. You need to register them on a
     * type before to get meaningful values.
     *
     * @param out_vec On out: Stores all gradients.
     * @param idx Start index of out_vec where the gradients should be stored.
     */
    template<class T>
    void get_gradient(T &out_vec, uint32_t &idx) const
    {
        for(auto &c: this->constants)
        {
            out_vec[idx] = c.getGradient();
            idx++;
        }
    }

    void put(pt::ptree &ptree, const std::string &type_name) const
    {
        if(perturbed_idx.empty())
            return;

        pt::ptree perturbed;

        for(uint32_t idx: perturbed_idx)
        {
            perturbed.put(std::to_string(idx), constants[idx]);
        }
        pt::ptree perturbed_vals;
        perturbed_vals.add_child("perturbed", perturbed);
        ptree.add_child(type_name, perturbed_vals);
    }

    /**
     * Set any perturbed parameter from the property tree.
     *
     * @params ptree Property tree with key = idx of constants, values =
     *         the perturbed values and one list 'perturbed' of indices.
     *
     * @returns Errorcode
     */
    int from_pt(pt::ptree &ptree)
    {
        int err = 0;
        for(auto &it: ptree.get_child("perturbed"))
        {
            uint32_t idx = std::stoi(it.first);
            this->constants[idx] = it.second.get_value<double>();
            perturbed_idx.push_back(idx);
        }
        return err;
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
    A look_lo(A x) const
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
    A look_up(A x) const
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
     * with high accuracy as function of a in the range in [0;20], but can
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

    /**
     * Initial id of this simulation. Emerging ensembles from this
     * have other ids which are set by reading a checkpoint file.
     */
    std::string id = "0";
    /**
     *  Running id for the ensembles starting from this instance.
     */
    uint64_t ensemble_id = 0;
    //
    // Physical constants warm cloud
    //
    double alpha_d; /*!< Accomodation coefficient */

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

    // See constants.h for a description of those.
    std::vector<codi::RealReverse> constants;

    /**
     * Store any idx from perturbed parameters.
     */
    std::vector<uint32_t> perturbed_idx;

    /**
     * Structure to hold the new equidistant lookup table for
     * graupel wetgrowth diameter
     */
    table_t ltabdminwgg;
    gamma_table_t table_g1, table_g2, table_r1, table_r2, table_r3;

    model_constants_t()
    {
        constants.resize(static_cast<int>(Cons_idx::n_items));
        std::fill(constants.begin(), constants.end(), 0);
    }

    /**
     * Register the model parameters on the tape for codi::RealReverse.
     *
     * @param tape Tape where the parameters should be registered on.
     */
    void register_input(codi::RealReverse::TapeType &tape)
    {
        for(auto &c: this->constants)
            tape.registerInput(c);
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
    void get_gradient(T &out_vec) const
    {
        for(int i=0; i<static_cast<int>(Cons_idx::n_items); ++i)
            out_vec[i] = this->constants[i].getGradient();
#if defined(RK4ICE) || defined(RK4NOICE)

        uint32_t idx = static_cast<uint32_t>(Cons_idx::n_items);
        this->rain.get_gradient(out_vec, idx);
        this->cloud.get_gradient(out_vec, idx);
        this->graupel.get_gradient(out_vec, idx);
        this->hail.get_gradient(out_vec, idx);
        this->ice.get_gradient(out_vec, idx);
        this->snow.get_gradient(out_vec, idx);
#endif
    }

    /**
     * Put any perturbed parameter to a property tree.
     * This will compare the parameters to the constants
     * available in constants.h, assuming this is only called during
     * checkpoint writing.
     *
     * @params ptree Property tree, where a tree "model_constants" is being added.
     */
    void put(pt::ptree &ptree) const
    {
        pt::ptree model_cons;
        model_cons.put("id", id);

        // technical parameters
        model_cons.put("t_end_prime", t_end_prime);
        model_cons.put("t_end", t_end);
        model_cons.put("dt_prime", dt_prime);
        model_cons.put("dt", dt);
        model_cons.put("dt_traject_prime", dt_traject_prime);
        model_cons.put("dt_traject", dt_traject);
        model_cons.put("num_steps", num_steps);
        model_cons.put("num_sub_steps", num_sub_steps);
        model_cons.put("dt_half", dt_half);
        model_cons.put("dt_sixth", dt_sixth);
        model_cons.put("dt_third", dt_third);
        model_cons.put("ensemble_id", ensemble_id);

        if(!perturbed_idx.empty())
        {
            pt::ptree perturbed;
            for(uint32_t idx: perturbed_idx)
                perturbed.put(std::to_string(idx), constants[idx]);
            model_cons.add_child("perturbed", perturbed);
        }

        // particle_model_constants
        hail.put(model_cons, "hail");
        ice.put(model_cons, "ice");
        snow.put(model_cons, "snow");
        cloud.put(model_cons, "cloud");
        rain.put(model_cons, "rain");
        graupel.put(model_cons, "graupel");

        // collection coefficients depend completely on particle
        // model constants and can be derived from there. So
        // we skip adding that to the checkpoint.

        ptree.add_child("model_constants", model_cons);
    }

    int from_pt(pt::ptree &ptree)
    {
        int err = 0;
        for(auto &it: ptree.get_child("model_constants"))
        {
            auto first = it.first;
            if(first == "id")
            {
                id = it.second.get_value<std::string>() + "-" + id;
            } else if(first == "t_end_prime")
            {
                t_end_prime = it.second.get_value<double>();
            } else if(first == "t_end")
            {
                t_end = it.second.get_value<double>();
            } else if(first == "dt_prime")
            {
                dt_prime = it.second.get_value<double>();
            } else if(first == "dt")
            {
                dt = it.second.get_value<double>();
            } else if(first == "dt_traject_prime")
            {
                dt_traject_prime = it.second.get_value<double>();
            } else if(first == "dt_traject")
            {
                dt_traject = it.second.get_value<double>();
            } else if(first == "num_steps")
            {
                num_steps = it.second.get_value<uint64_t>();
            } else if(first == "num_sub_steps")
            {
                num_sub_steps = it.second.get_value<uint64_t>();
            } else if(first == "dt_half")
            {
                dt_half = it.second.get_value<double>();
            } else if(first == "dt_sixth")
            {
                dt_sixth = it.second.get_value<double>();
            } else if(first == "dt_third")
            {
                dt_third = it.second.get_value<double>();
            } else if(first == "ensemble_id")
            {
                ensemble_id = it.second.get_value<uint64_t>();
            } else if(first == "perturbed")
            {
                for(auto &it2: ptree.get_child("model_constants.perturbed"))
                {
                    uint32_t idx = std::stoi(it2.first);
                    perturbed_idx.push_back(idx);
                    this->constants[idx] = it2.second.get_value<double>();
                }
            // below from here: perturbed particle models
            } else if(first == "hail")
            {
                for(auto &it2: ptree.get_child("model_constants.hail.perturbed"))
                {
                    uint32_t idx = std::stoi(it2.first);
                    this->hail.perturbed_idx.push_back(idx);
                    this->hail.constants[idx] = it2.second.get_value<double>();
                }
            } else if(first == "ice")
            {
                for(auto &it2: ptree.get_child("model_constants.ice.perturbed"))
                {
                    uint32_t idx = std::stoi(it2.first);
                    this->ice.perturbed_idx.push_back(idx);
                    this->ice.constants[idx] = it2.second.get_value<double>();
                }
            } else if(first == "snow")
            {
                for(auto &it2: ptree.get_child("model_constants.snow.perturbed"))
                {
                    uint32_t idx = std::stoi(it2.first);
                    this->snow.perturbed_idx.push_back(idx);
                    this->snow.constants[idx] = it2.second.get_value<double>();
                }
            } else if(first == "cloud")
            {
                for(auto &it2: ptree.get_child("model_constants.cloud.perturbed"))
                {
                    uint32_t idx = std::stoi(it2.first);
                    this->cloud.perturbed_idx.push_back(idx);
                    this->cloud.constants[idx] = it2.second.get_value<double>();
                }
            } else if(first == "rain")
            {
                for(auto &it2: ptree.get_child("model_constants.rain.perturbed"))
                {
                    uint32_t idx = std::stoi(it2.first);
                    this->rain.perturbed_idx.push_back(idx);
                    this->rain.constants[idx] = it2.second.get_value<double>();
                }
            } else if(first == "graupel")
            {
                for(auto &it2: ptree.get_child("model_constants.graupel.perturbed"))
                {
                    uint32_t idx = std::stoi(it2.first);
                    this->graupel.perturbed_idx.push_back(idx);
                    this->graupel.constants[idx] = it2.second.get_value<double>();
                }
            } else
            {
                err = MODEL_CONS_CHECKPOINT_ERR;
            }
        }
        return err;
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

    int ens_config_flag; /*!< Configuration file for ensembles. */
    char* ens_config_string;

    int checkpoint_flag; /*!< Checkpoint file for the simulation. */
    char* checkpoint_string;

    int gnu_id_flag; /*!< ID given for this instance, i.e. thread_id or id by GNU parallel. */
    char* gnu_id_string;

    int folder_name_flag;
    char* folder_name_string;

    global_args_t()
    {
        final_time_flag = 0;
        final_time_string = nullptr;

        timestep_flag = 0;
        timestep_string = nullptr;

        snapshot_index_flag = 0;
        snapshot_index_string = nullptr;

        output_flag = 0;
        output_string = nullptr;

        scaling_fact_flag = 0;
        scaling_fact_string = nullptr;

        input_flag = 0;
        input_file = nullptr;

        start_over_flag = 0;
        start_over_string = nullptr;

        start_over_env_flag = 0;
        start_over_env_string = nullptr;

        fixed_iteration_flag = 0;
        fixed_iteration_string = nullptr;

        auto_type_flag = 0;
        auto_type_string = nullptr;

        traj_flag = 0;
        traj_string = nullptr;

        write_flag = 0;
        write_string = nullptr;

        progress_index_flag = 0;
        progress_index_string = nullptr;

#ifdef MET3D
        delay_start_flag = 0;
        delay_start_string = nullptr;
#endif

        ens_config_flag = 0;
        ens_config_string = nullptr;

        checkpoint_flag = 0;
        checkpoint_string = nullptr;

        gnu_id_flag = 0;
        gnu_id_string = nullptr;

        folder_name_flag = 0;
        folder_name_string = nullptr;
    }
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
    double start_time; /*!< start time in seconds relativ to ascend */
#endif
    int snapshot_index; /*!< Save a snapshot every snapshot_index iteration. */
    /**
     * Number of timesteps for the simulation between two
     * datapoints from the netCDF file.
     */
    uint64_t num_sub_steps;

    std::string OUTPUT_FILENAME; /*!< Filename for output. */
    std::string INPUT_FILENAME; /*!< Filename for input netCDF file. */
    std::string ENS_CONFIG_FILENAME; /*!< Filename for ensemble configuration file. */
    std::string CHECKPOINT_FILENAME; /*!< Filename for checkpoint file. */
    std::string FOLDER_NAME; /*!< Folder name for newly generated checkpoints. */
    uint32_t id; /*!< ID given for this instance, i.e. thread_id or id by GNU parallel. */

    bool start_over; /*!< Start over at new timestep of trajectory? */
    bool start_over_env; /*!< Start over environment variables at new timestep of trajectory? */
    bool fixed_iteration; /*!< Fix temperature and pressure at every iteration? */

    double scaling_fact; /*!< Scaling factor. */

    uint32_t auto_type; /*!< Particle type. */
    uint32_t traj; /*!< Trajectory index to load from the netCDF file. */
    uint32_t write_index; /*!< Write stringstream every x iterations to disk. */
    uint64_t progress_index; /*!< Index for updating progressbar. */
    uint32_t ensemble; /*!< Index of ensemble. */
    double current_time; /*!< Time for and from checkpoint files. */

    input_parameters_t()
    {
        // Numerics
        t_end_prime = 100.0;	// Seconds
        dt_prime = 0.01;		// Seconds
        snapshot_index = 200;
        dt_traject = 20;       // Seconds; fixed from paper
        // Filename for output
#if defined(RK4)
        OUTPUT_FILENAME = "data/id0_rain_OUTPUT.txt";
#endif
#if defined(RK4NOICE)
        OUTPUT_FILENAME = "data/id0_sb_OUTPUT.txt";
#endif
#if defined(RK4ICE)
        OUTPUT_FILENAME = "data/id0_sb_ice_OUTPUT.txt";
#endif
        CHECKPOINT_FILENAME = "";
        ENS_CONFIG_FILENAME = "";
        FOLDER_NAME = "";
        // Filename for input
        INPUT_FILENAME = "/mnt/localscratch/data/project/m2_jgu-tapt/online_trajectories/foehn201305_case/foehn201305_warming.nc";

        // Scaling factor
        scaling_fact = 1.0;	// No scaling
        start_over = true;
        start_over_env = true;
        fixed_iteration = false;
        auto_type = 3;
        traj = 0;
        write_index = 100000;
        progress_index = 1000;
        ensemble = 0;
#ifdef MET3D
        start_time = std::nan("");
#endif
        current_time = std::nan("");;
        id = 0;
    }

    /**
     * Change the output filename such that it starts with "idx-y-z_"
     * where x-y are the ids of preceding trajectories and z is the current id
     *
     * @params all_ids String of form x-y-z with x,y and z positive ints
     *                 and z being the current id.
     */
    void set_outputfile_id(
        const std::string &all_ids,
        uint64_t ensemble_id)
    {
        std::string preceeding_ids = "";
        if(all_ids != "0")
        {
            auto found = all_ids.find_last_of("-");
            preceeding_ids = all_ids.substr(0, found);
        }

        std::string ens_string = "_ensID_";
        if(ensemble_id == 0)
            ens_string += "000";
        else
            for(int64_t i=ensemble_id; i<1000; i*=10)
                ens_string += "0";

        // master branch
        if(preceeding_ids == "")
        {
            std::string id_str = "id" + std::to_string(id) + ens_string + std::to_string(ensemble_id) + "_";
            auto pos = OUTPUT_FILENAME.rfind("/");
            if(pos == std::string::npos)
            {
                OUTPUT_FILENAME.insert(0, id_str);
            } else
            {
                OUTPUT_FILENAME.insert(pos+1, id_str);
            }

        } else // subsequent ensemble members
        {
            std::string id_str = "id" + all_ids + ens_string + std::to_string(ensemble_id) + "_";
            std::string to_replace = "id" + preceeding_ids + "_";
            auto pos = OUTPUT_FILENAME.rfind(to_replace);

            if(pos != std::string::npos)
            {
                // This *should* always be the correct branch
                OUTPUT_FILENAME.replace(pos, to_replace.length()+11, id_str);
            } else
            {
                // If for any weird reason, "id" is not part of OUTPUT_FILENAME:
                // Add to beginning of the filename. Check for any path characters
                auto pos_folder = OUTPUT_FILENAME.rfind("/");
                if(pos == std::string::npos)
                {
                    OUTPUT_FILENAME.insert(0, id_str);
                } else
                {
                    OUTPUT_FILENAME.insert(pos_folder+1, id_str);
                }
            }
        }
    }

    void put(pt::ptree &ptree, const double &time) const
    {
        pt::ptree input_params;
        input_params.put<double>("t_end_prime", t_end_prime);
        input_params.put<double>("dt_prime", dt_prime);
        input_params.put<double>("dt_traject_prime", dt_traject_prime);
        input_params.put<double>("dt_traject", dt_traject);
#ifdef MET3D
        input_params.put<double>("start_time", start_time);
#endif
        input_params.put<int>("snapshot_index", snapshot_index);
        input_params.put<uint64_t>("num_sub_steps", num_sub_steps);
        input_params.put<std::string>("OUTPUT_FILENAME", OUTPUT_FILENAME);
        input_params.put<std::string>("INPUT_FILENAME", INPUT_FILENAME);
        input_params.put<bool>("start_over", start_over);
        input_params.put<bool>("start_over_env", start_over_env);
        input_params.put<bool>("fixed_iteration", fixed_iteration);
        input_params.put<double>("scaling_fact", scaling_fact);
        input_params.put<uint32_t>("auto_type", auto_type);
        input_params.put<uint32_t>("traj", traj);
        input_params.put<uint32_t>("write_index", write_index);
        input_params.put<uint64_t>("progress_index", progress_index);
        input_params.put<uint32_t>("ensemble", ensemble);
        input_params.put<double>("current_time", time);
        input_params.put<std::string>("FOLDER_NAME", FOLDER_NAME);
        ptree.add_child("input_params", input_params);
    }

    void put(pt::ptree &ptree) const
    {
        put(ptree, current_time);
    }

    /**
     * Set values from property tree used in reading checkpoint files.
     */
    int from_pt(pt::ptree &ptree)
    {
        int err = 0;
        for(auto &it: ptree.get_child("input_params"))
        {
            auto first = it.first;
            if(first == "t_end_prime")
            {
                t_end_prime = it.second.get_value<double>();
            } else if(first == "dt_prime")
            {
                dt_prime = it.second.get_value<double>();
            } else if(first == "dt_traject_prime")
            {
                dt_traject_prime = it.second.get_value<double>();
            } else if(first == "dt_traject")
            {
                dt_traject = it.second.get_value<double>();
#ifdef MET3D
            } else if(first == "start_time")
            {
                start_time = it.second.get_value<double>();
#endif
            } else if(first == "snapshot_index")
            {
                snapshot_index = it.second.get_value<int>();
            } else if(first == "num_sub_steps")
            {
                num_sub_steps = it.second.get_value<uint64_t>();
            } else if(first == "OUTPUT_FILENAME")
            {
                OUTPUT_FILENAME = it.second.data();
            } else if(first == "INPUT_FILENAME")
            {
                INPUT_FILENAME = it.second.data();
            } else if(first == "start_over")
            {
                start_over = (it.second.data() == "1"
                    || it.second.data() == "true") ? true : false;
            } else if(first == "start_over_env")
            {
                start_over_env = (it.second.data() == "1"
                    || it.second.data() == "true") ? true : false;
            } else if(first == "fixed_iteration")
            {
                fixed_iteration = (it.second.data() == "1"
                    || it.second.data() == "true") ? true : false;
            } else if(first == "scaling_fact")
            {
                scaling_fact = it.second.get_value<double>();
            } else if(first == "auto_type")
            {
                auto_type = it.second.get_value<uint32_t>();
            } else if(first == "traj")
            {
                traj = it.second.get_value<uint32_t>();
            } else if(first == "write_index")
            {
                write_index = it.second.get_value<uint32_t>();
            } else if(first == "progress_index")
            {
                progress_index = it.second.get_value<uint64_t>();
            } else if(first == "ensemble")
            {
                ensemble = it.second.get_value<uint32_t>();
            } else if(first == "current_time")
            {
                current_time = it.second.get_value<double>();
            } else if(first == "FOLDER_NAME")
            {
                FOLDER_NAME = it.second.data();
            } else
            {
                err = INPUT_NAME_CHECKPOINT_ERR;
            }
        }
        return err;
    }
};

struct param_t{
    double mean;
    double sigma;
    double sigma_perc;
    int err;
    int name;
    int out_name;
    bool particle_param;

    std::function<double()> get_rand; /*!< distribution used for random number generation. */
    std::normal_distribution<double> normal_dis;
    std::uniform_real_distribution<double> uniform_dis;
    std::string param_name;
    std::string outparam_name;
    std::string func_name;

    enum class OutParam: uint32_t {
        model, cloud, rain, ice, graupel, hail, snow
    };
    std::unordered_map<std::string, OutParam> const table_out_param = {
        {"model", OutParam::model},
        {"cloud", OutParam::cloud}, {"rain", OutParam::rain},
        {"ice", OutParam::ice}, {"graupel", OutParam::graupel},
        {"hail", OutParam::hail}, {"snow", OutParam::snow}
    };

    param_t()
    {
        mean            = std::nan("");
        sigma           = std::nan("");
        sigma_perc      = std::nan("");
        err             = 0;
        name            = -1;
        out_name        = -1;
        particle_param  = false;
        func_name       = "";
    }

    param_t(std::string param_type)
    {
       add_type(param_type);
    }

    void add_type(std::string param_type)
    {
        auto it = table_out_param.find(param_type);
        if(it != table_out_param.end())
        {
            out_name = static_cast<int>(it->second);
            if(out_name != static_cast<int>(OutParam::model))
                particle_param = true;
            outparam_name = param_type;
        } else
        {
            err = OUTPARAM_CONFIG_ERR;
        }
    }

    void add_mean(double m)
    {
        mean = m;
        if(!isnan(sigma_perc) && isnan(sigma))
            sigma = mean*sigma_perc/100;
    }

    void add_name(std::string n, model_constants_t &cc)
    {
        param_name = n;
        if(particle_param)
        {
            auto it = table_particle_param.find(n);
            if(it != table_particle_param.end())
            {
                name = static_cast<uint32_t>(it->second);
                if(std::isnan(mean))
                {
                    particle_model_constants_t *pt_model;
                    switch(out_name)
                    {
                        case static_cast<uint32_t>(OutParam::cloud):
                            pt_model = &(cc.cloud);
                            break;
                        case static_cast<uint32_t>(OutParam::rain):
                            pt_model = &(cc.rain);
                            break;
                        case static_cast<uint32_t>(OutParam::snow):
                            pt_model = &(cc.snow);
                            break;
                        case static_cast<uint32_t>(OutParam::graupel):
                            pt_model = &(cc.graupel);
                            break;
                        case static_cast<uint32_t>(OutParam::hail):
                            pt_model = &(cc.hail);
                            break;
                        case static_cast<uint32_t>(OutParam::ice):
                            pt_model = &(cc.ice);
                            break;
                    }
                    mean = pt_model->constants[name].getValue();
                }
            } else
            {
                err = PARAM_CONFIG_ERR;
            }
        } else
        {
            auto it = table_param.find(n);
            if(it != table_param.end())
            {
                name = static_cast<int>(it->second);
                if(std::isnan(mean))
                    mean = get_at(cc.constants, name).getValue();
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
        sigma = s;
        if(!isnan(mean) && func_name != "")
            add_rand_function(func_name);
    }

    void add_sigma_perc(double s)
    {
        sigma_perc = s;
        if(!isnan(mean))
        {
            sigma = sigma_perc*mean/100;
            if(func_name != "")
                add_rand_function(func_name);
        }
    }

    void add_rand_function(std::string name)
    {
        if(func_name == "")
            func_name = name;
        if(!isnan(sigma_perc) && isnan(sigma) && !isnan(mean))
            sigma = mean*sigma_perc/100;
        if(!isnan(mean) && !isnan(sigma))
        {
            if(func_name == "normal")
            {
                normal_dis = std::normal_distribution<double>(mean, sigma);
                get_rand = std::bind(normal_dis, rand_generator);
            } else if(func_name == "uniform")
            {
                uniform_dis = std::uniform_real_distribution<double>(mean-sigma, mean+sigma);
                get_rand = std::bind(uniform_dis, rand_generator);
            } else
            {
                err = DISTRIBUTION_CONFIG_ERR;
            }
        }
    }

    int check()
    {
        if(name == -1)
        {
            std::cerr << "Error in config file:\n"
                      << "You did not specify the parameter to perturb or "
                      << "you have a typo at <name>typo</name>\n";
            err = MISSING_PARAM_CONFIG_ERR;
            return err;
        }
        if(isnan(sigma) && isnan(sigma_perc))
        {
            std::cerr << "Error in config file:\n"
                      << "You did not specify the variance for "
                      << "perturbing the parameter.\n";
            err = MISSING_VARIANCE_CONFIG_ERR;
            return err;
        }
        if(func_name == "")
            err = DISTRIBUTION_CONFIG_ERR;

        switch(err)
        {
            case PARAM_CONFIG_ERR:
                std::cerr << "Error in config file:\n"
                          << "You used a parameter to perturb that does "
                          << "not exist.\n";
                return err;

            case OUTPARAM_CONFIG_ERR:
                std::cerr << "Error in config file:\n"
                          << "You used an output parameter that does "
                          << "not exist.\n";
                return err;
            case DISTRIBUTION_CONFIG_ERR:
                std::cerr << "Error in config file:\n"
                          << "No such function for generating random "
                          << "numbers in perturbing parameters:\n"
                          << "name: " << func_name << "\n"
                          << "Options are:\n"
                          << "normal: normal distribution with mean=mean "
                          << "or mean=default value of parameter and sigma"
                          << "=sigma or sigma=mean*sigma_perc/100\n"
                          << "uniform: uniform distribution from mean-sigma "
                          << "to mean+sigma and mean, sigma as above\n";
                return err;
            default:
                return err;

        }
    }

    void put(pt::ptree &ptree) const
    {
        if(err != 0)
            return;
        pt::ptree param;
        param.put("mean", mean);
        param.put("name", param_name);
        if(!isnan(sigma))
            param.put("sigma", sigma);
        else
            param.put("sigma_perc", sigma_perc);
        param.put("rand_func", func_name);
        param.put("type", outparam_name);
        // ptree.add_child("params.", param);
        ptree.push_back(std::make_pair("", param));
    }

    int from_pt(pt::ptree &ptree, model_constants_t &cc)
    {
        int err = 0;
        add_type(ptree.get<std::string>("type"));
        if(err != 0)
            return err;
        for(auto &it: ptree)
        {
            auto first = it.first;
            if(first == "name")
            {
                add_name(it.second.get_value<std::string>(), cc);
            } else if(first == "sigma_perc")
            {
                add_sigma_perc(it.second.get_value<double>());
            } else if(first == "sigma")
            {
                add_sigma(it.second.get_value<double>());
            } else if(first == "mean")
            {
                add_mean(it.second.get_value<double>());
            } else if(first == "type")
            {
                // Needs to be done before the for-loop.
            } else if(first == "rand_func")
            {
                add_rand_function(it.second.get_value<std::string>());
            }
            else
            {
                err = PARAM_CONFIG_ERR;
            }
        }
        return err;
    }

    void perturb(model_constants_t &cc) const
    {
        if(particle_param)
        {
            particle_model_constants_t *pt_model = nullptr;
            switch(out_name)
            {
                case static_cast<uint32_t>(OutParam::cloud):
                    pt_model = &(cc.cloud);
                    break;
                case static_cast<uint32_t>(OutParam::rain):
                    pt_model = &(cc.rain);
                    break;
                case static_cast<uint32_t>(OutParam::snow):
                    pt_model = &(cc.snow);
                    break;
                case static_cast<uint32_t>(OutParam::graupel):
                    pt_model = &(cc.graupel);
                    break;
                case static_cast<uint32_t>(OutParam::hail):
                    pt_model = &(cc.hail);
                    break;
                case static_cast<uint32_t>(OutParam::ice):
                    pt_model = &(cc.ice);
                    break;
                default:
                    std::cout << "Error in perturbing...\n";
            }
            pt_model->constants[name] = get_rand();
            pt_model->perturbed_idx.push_back(name);
        } else
        {
            cc.constants[name] = get_rand();
            cc.perturbed_idx.push_back(name);
        }
    }
};

struct segment_t
{
    std::vector<param_t> params;
    double value;
    /**
     * Used to track the value for perturbing in case the tracked value
     * is not exactly equal to value but rather jumps around value.
     */
    double old_value;
    uint32_t n_members;
    int value_name;
    int value_name_sig; /*<! Which parameter had the most significant impact so far >*/
    int out_param; /*<! Which output parameter in case of sensitivity methods >*/
    uint32_t n_segments;
    uint32_t err;
    uint32_t old_sign;
    double duration; /*<! Maximum duration integration time for this ensemble >*/
    std::unordered_map<std::string, std::string> tree_strings;
    bool activated;
    /**
     * Allowed tolerance for value_method. If the value jumps, i.e. due to
     * large time steps, it may jump by tol*value from below/above value
     * to above/below it.
     */
    const double tol = 1e-5;

    enum Method {
        impact_change, sign_flip, value_method, repeated_time
    };
    std::unordered_map<std::string, Method> const table_method = {
        {"impact_change", Method::impact_change}, {"sign_flip", Method::sign_flip},
        {"value", Method::value_method}, {"repeated_time", Method::repeated_time}
    };
    int method;
    enum Param {
        pressure, T, w, S, QC, QR, QV, NCCLOUD, NCRAIN,
        QI, NCICE, QS, NCSNOW, QG, NCGRAUPEL, QH, NCHAIL,
        QI_OUT, QS_OUT, QR_OUT, QG_OUT, QH_OUT, latent_heat,
        latent_cool, NI_OUT, NS_OUT, NR_OUT, NG_OUT,
        NH_OUT, z, Inactive, deposition, sublimination,
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
        const5, D_rainfrz_gh, D_rainfrz_ig, dv0, p_sat_melt, cp,
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
        {"Inactive", Param::Inactive}, {"deposition", Param::deposition},
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
        {"const5", Param::const5}, {"D_rainfrz_gh", Param::D_rainfrz_gh},
        {"D_rainfrz_ig", Param::D_rainfrz_ig},
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

    segment_t()
    {
        value           = std::nan("");
        old_value       = std::nan("");
        n_members       = 1;
        value_name      = -1;
        value_name_sig  = -1;
        out_param       = -1;
        n_segments      = 1;
        err             = 0;
        old_sign        = 0;
        duration        = 0;
        activated       = false;
        method          = value_method;
    }

    // enum OutParam {
    //     vapor, cloud, rain, ice, graupel, hail, snow
    // };
    // std::unordered_map<std::string, OutParam> const table_out_param = {
    //     {"snow", OutParam::snow}, {"vapor", OutParam::vapor},
    //     {"cloud", OutParam::cloud}, {"rain", OutParam::rain},
    //     {"ice", OutParam::ice}, {"graupel", OutParam::graupel},
    //     {"hail", OutParam::hail}
    // };


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
        if(method == impact_change)
        {
            value_name = -1;
            return;
        }
        auto it = table_param.find(n);
        if(it != table_param.end())
        {
            value_name = it->second;
            std::string key = "when_name";
            tree_strings.insert(std::make_pair(key, n));
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
            std::string key = "when_method";
            tree_strings.insert(std::make_pair(key, m));
            if(method == impact_change)
                value_name = -1;
        } else
        {
            err = METHOD_CONFIG_ERR;
        }
    }

    void add_counter(uint32_t c)
    {
        n_segments = c;
    }

    void add_duration(double t)
    {
        duration = t;
    }

    void add_out_param(std::string p)
    {
        auto it = std::find(output_par_idx.begin(), output_par_idx.end(), p);
        if(it != output_par_idx.end())
        {
            out_param = std::distance(output_par_idx.begin(), it);
            std::string key = "when_sens";
            tree_strings.insert(std::make_pair(key, p));
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
        else if(n_segments < 0)
            err = N_SEGMENTS_CONFIG_ERR;
        else if( value_name == -1
                && method != Method::impact_change
                && method != Method::repeated_time )
            err = VALUE_NAME_CONFIG_ERR;
        else if(method == -1)
            err = METHOD_CONFIG_ERR;
        else if(out_param == -1 && value_name>=Param::rain_a_geo)
            err = OUTPARAM_CONFIG_ERR;
        switch(err)
        {
            case VALUE_NAME_CONFIG_ERR:
                std::cerr << "Error in config file:\n"
                          << "The name of the parameter used to determine "
                          << "when to start an ensemble is not defined or "
                          << "does not exist.\n";
                n_segments = 0;
                return err;
            case METHOD_CONFIG_ERR:
                std::cerr << "Error in config file:\n"
                          << "The name of the method used to determine "
                          << "when to start an ensemble is not defined or "
                          << "does not exist.\n";
                n_segments = 0;
                return err;
            case OUTPARAM_CONFIG_ERR:
                std::cerr << "Error in config file:\n"
                          << "The name of the output parameter used to determine "
                          << "when to start an ensemble is not defined or "
                          << "does not exist.\n";
                n_segments = 0;
                return err;
            case N_MEMBERS_CONFIG_ERR:
                std::cerr << "Error in config file:\n"
                          << "The number of members for an ensemble must be "
                          << "2 or higher. One simulation always runs without "
                          << "perturbed parameters for comparison.\n";
                n_segments = 0;
                return err;
            case N_SEGMENTS_CONFIG_ERR:
                std::cerr << "Error in config file:\n"
                          << "The number of segments must be at least 1 in "
                          << "order to start at least 1 ensemble.\n";
                n_segments = 0;
                return err;
        }
        if(err != 0)
            n_segments = 0;
        return err;
    }

    /**
     * Check if perturbing is necessary.
     *
     */
    template<class float_t>
    bool perturb_check(
        const model_constants_t &cc,
        const std::vector< std::array<double, num_par > > &gradients,
        const std::vector<float_t> &y,
        const double timestep)
    {
        if(n_segments == 0)
            return false;
        int idx;
        switch(method)
        {
            case impact_change:
            {
                auto it_minmax = std::minmax_element(gradients[out_param].begin(),
                    gradients[out_param].end());
                auto it_max = (*it_minmax.second > abs(*it_minmax.first)) ? it_minmax.second : it_minmax.first;
                idx = std::distance(gradients[out_param].begin(), it_max);
                if(value_name_sig != -1)
                {
                    if(idx != value_name_sig)
                    {
                        value_name_sig = idx;
                        activated = true;
                        return true;
                    }
                    return false;
                } else // First time we check the highest impact
                {
                    value_name_sig = idx;
                    return false;
                }
                return false;
            }
            case sign_flip:
            {
                // the first num_comp many values refer to output parameters
                idx = value_name - num_comp;
                // std::cout << "idx " << idx << ", grad " << gradients[out_param][idx] << "\n"
                //           << "old_sign " << old_sign << ", out_param " << out_param << "\n";
                if(old_sign == 0)
                {
                    // set sign; no perturbing needed.
                    if(gradients[out_param][idx] == 0)
                        return false;
                    old_sign = (gradients[out_param][idx] > 0) ? 1 : 2;
                    return false;
                } else
                {
                    // check if a sign change happend
                    if(old_sign == 1 && gradients[out_param][idx] < 0)
                    {
                        activated = true;
                        old_sign = 2;
                        return true;
                    }else if(old_sign == 2 && gradients[out_param][idx] > 0)
                    {
                        // Perturb parameters
                        activated = true;
                        old_sign = 1;
                        return true;
                    }
                    return false;
                }
                return false;
            }
            case value_method:
            {
                double current_value;
                // If value is an output parameter or a sensitivity
                if(out_param != -1)
                    current_value = gradients[out_param][value_name - num_comp];
                else
                {
                    current_value = y[value_name].getValue();
                    if(value_name == T_idx)
                        current_value *= 273.15;
                    else if(value_name == p_idx)
                        current_value *= 1.0e5;
                }
                if(current_value == value)
                {
                    activated = true;
                    return true;
                }

                if(std::isnan(old_value))
                {
                    old_value = current_value;
                    return false;
                }

                if( signbit(value-current_value) != signbit(value-old_value)
                    || fabs(value-current_value) < fabs(value*tol) )
                {
                    activated = true;
                    return true;
                }
                old_value = current_value;
                return false;
            case repeated_time:
                if(std::fmod(timestep, value) == 0)
                {
                    activated = true;
                    return true;
                }
                return false;
            }
        }
        return false;
    }

    /**
     * Deactivates the segment so it is not used anymore/decreases the amount
     * of possible usages of this segment. If this is a repeated_time segment
     * without a fixed duration it will not decrease the possible amounts of
     * usages. The same applies with duration if keep_repeated is true, which
     * can be used to spawn huge ensembles for a limited time.
     */
    void deactivate(const bool keep_repeated=false)
    {
        if(!activated) return;
        // Perturbing every few timesteps is done indefenitely
        // unless a fixed duration is given at which the ensembles should stop
        if( (method != repeated_time && n_segments > 0)
            || (method == repeated_time && n_segments > 0 && duration != 0 && !keep_repeated) )
            n_segments--;
        activated = false;
    }

    void perturb(model_constants_t &cc, const reference_quantities_t &ref_quant,
        input_parameters_t &input)
    {
        // Sanity check if had been done already
        if(n_segments == 0)
            return;

        // Change the number of time steps if a fixed duration is given
        if(duration != 0)
        {
            if(duration + input.current_time + input.start_time < cc.t_end_prime)
            {
                cc.t_end_prime = duration + input.current_time + input.start_time;
                cc.t_end = cc.t_end_prime/ref_quant.tref;
                input.t_end_prime = duration;
            }
        }
        // Perturb every param
        for(auto &p: params)
            p.perturb(cc);
        // When perturbing is done, deativate
        deactivate();
    }

    void put(pt::ptree &ptree) const
    {
        if( (err != 0 || n_segments < 1) && !activated)
            return;
        pt::ptree segment;
        if(!isnan(value))
        {
            segment.put("when_value", value);
        }
        if(value_name != -1)
        {
            segment.put("when_name", tree_strings.find("when_name")->second);
        }
        if(out_param != -1)
        {
            segment.put("when_sens", tree_strings.find("when_sens")->second);
        }
        if(n_segments != 1)
        {
            segment.put("when_counter", n_segments);
        }
        if(method != value_method)
        {
            segment.put("when_method",  tree_strings.find("when_method")->second);
        }
        if(n_members != 1)
        {
            segment.put("amount", n_members);
        }
        if(activated)
        {
            segment.put("activated", true);
        }
        if(duration > 0)
        {
            segment.put("duration", duration);
        }
        pt::ptree param_tree;
        for(auto &p: params)
            p.put(param_tree);
        segment.add_child("params", param_tree);
        ptree.push_back(std::make_pair("", segment));
    }

    /**
     * Used to read from a checkpoint file where the mean for the gaussians
     * to draw from is given.
     */
    int from_pt(pt::ptree &ptree, model_constants_t &cc)
    {
        int err = 0;
        for(auto &it: ptree)
        {
            auto first = it.first;
            if(first == "when_value")
            {
                add_value(it.second.get_value<double>());
            } else if(first == "when_name")
            {
                add_value_name(it.second.get_value<std::string>());
            } else if(first == "amount")
            {
                add_amount(it.second.get_value<uint32_t>());
            } else if(first == "when_method")
            {
                add_method(it.second.get_value<std::string>());
            } else if(first == "when_counter")
            {
                add_counter(it.second.get_value<uint32_t>());
            } else if(first == "when_sens")
            {
                add_out_param(it.second.get_value<std::string>());
            } else if(first == "activated")
            {
                activated = it.second.get_value<bool>();
            } else if(first == "duration")
            {
                add_duration(it.second.get_value<double>());
            } else if(first == "params")
            {
                for(auto &param_it: ptree.get_child(first))
                {
                    param_t param;
                    err = param.from_pt(param_it.second, cc);
                    add_param(param);
                }
            } else
            {
                err = SEGMENTS_CHECKPOINT_ERR;
            }
        }
        return err;
    }
};


struct IO_handle_t{
    // for txt files
    std::stringstream out_tmp;
    std::ofstream outfile;
    std::ofstream out_diff[num_comp];
    std::stringstream out_diff_tmp[num_comp];
    uint64_t n_snapshots; // number of buffered snapshots
    uint64_t flushed_snapshots;

    std::string filetype;
    // for netCDF files and a vector for each column
    // columns: output parameter + gradients + time_afer_ascent + type + flags
    // fast index: record
    // each array = one column
    // slow index: num_comp
    std::array<std::vector<double>, num_comp+num_par+4 > output_buffer;
    std::array<std::vector<unsigned char>, 4 > output_buffer_flags;
    std::array<std::vector<std::string>, 1 > output_buffer_str;
    std::array<std::vector<uint64_t>, 1 > output_buffer_int;
    NcFile datafile;
    std::vector<NcDim> dim_vector;
    std::vector<NcVar> var_vector;
    uint64_t n_ens;
    uint64_t n_trajs;
    uint64_t total_snapshots;
    std::string filename;

    IO_handle_t(
        const std::string filetype,
        const std::string filename,
        const uint64_t n_trajs,
        const uint64_t n_ens,
        const model_constants_t &cc,
        const reference_quantities_t &ref_quant,
        const std::string in_filename,
        const uint32_t write_index,
        const uint32_t snapshot_index)
    {
        this->n_trajs = n_trajs;
        this->n_ens = n_ens;
        this->filetype = filetype;
        this->filename = filename;
        if(filetype == "netcdf")
        {
            const int deflateLevel = 9; // compression level
            const bool enableShuffleFilter = true; // increases compression with little cost
            const bool enableDeflateFilter = true; // enable compression
            flushed_snapshots = 0;
            n_snapshots = 0;
            // Allocate memory for the buffer
            // maximum number of snapshots we are going to get
            total_snapshots = std::ceil( ((float)write_index)/snapshot_index ) + 1;
            const uint64_t vec_size = n_trajs * n_ens * total_snapshots; // n_snapshots * num_comp;
            const uint64_t vec_size_grad = num_comp * n_trajs * n_ens * total_snapshots;
            for(uint32_t i=0; i<num_comp; i++)
                output_buffer[i].resize(vec_size);
            for(uint32_t i=num_comp; i<num_comp+num_par; i++)
                output_buffer[i].resize(vec_size_grad);

            output_buffer[num_comp+num_par].resize(vec_size);          // time after ascent
            output_buffer[num_comp+num_par+1].resize(total_snapshots); // just time index
            output_buffer[num_comp+num_par+2].resize(vec_size);        // lat
            output_buffer[num_comp+num_par+3].resize(vec_size);        // lon

            for(uint32_t i=0; i<output_buffer_flags.size(); i++)
                output_buffer_flags[i].resize(vec_size);
            for(uint32_t i=0; i<output_buffer_str.size(); i++)
                output_buffer_str[i].resize(vec_size);
            for(uint32_t i=0; i<output_buffer_int.size(); i++)
                output_buffer_int[i].resize(vec_size);
            try
            {
                datafile.open(filename + ".nc_wcb", NcFile::newFile);
            }
            catch(const std::exception& e)
            {
                std::cerr << "Error creating output file:\n"
                          << filename << ".nc_wcb Does the file already exist?"
                          << "\nAborting";
                std::cerr << e.what() << '\n';
                exit(NC_ERR);
            }

             // Add dimensions
            NcDim param_dim = datafile.addDim("Output Parameter", num_comp);
            NcDim ens_dim = datafile.addDim("ensemble", n_ens);
            NcDim traj_dim = datafile.addDim("trajectory", n_trajs);
            NcDim time_dim = datafile.addDim("time"); // unlimited dimension

            NcVar param_var = datafile.addVar("Output Parameter", ncString, param_dim);
            NcVar ens_var = datafile.addVar("ensemble", ncInt64, ens_dim);
            NcVar traj_var = datafile.addVar("trajectory", ncString, traj_dim);

            ens_var.setCompression(
                enableShuffleFilter, enableDeflateFilter, deflateLevel);
            traj_var.setCompression(
                enableShuffleFilter, enableDeflateFilter, deflateLevel);
            param_var.setCompression(
                enableShuffleFilter, enableDeflateFilter, deflateLevel);

            // Add dim data
            // the ensemble number is taken from the filename after idxx_y
            // where y is the ensemble number
            auto sub_str_it = filename.find("_wcb");
            uint32_t ens = std::stoi(filename.substr(sub_str_it-4, sub_str_it-1));
            ens_var.putVar(&ens);
            traj_var.putVar(&cc.id);

            std::vector<uint64_t> startp, countp;
            startp.push_back(0);
            countp.push_back(1);
            for(const auto &p: output_par_idx)
            {
                param_var.putVar(startp, countp, &p);
                startp[0]++;
            }

            dim_vector.push_back(time_dim);
            dim_vector.push_back(traj_dim);
            dim_vector.push_back(ens_dim);

            for(auto &out_par: output_par_idx)
                var_vector.emplace_back(datafile.addVar(out_par, ncDouble, dim_vector));
            std::vector<NcDim> dim_vector_grad;
            dim_vector_grad.push_back(time_dim);
            dim_vector_grad.push_back(traj_dim);
            dim_vector_grad.push_back(ens_dim);
            dim_vector_grad.push_back(param_dim);
            for(auto &out_grad: output_grad_idx)
                var_vector.emplace_back(datafile.addVar(out_grad, ncDouble, dim_vector_grad));

            var_vector.emplace_back(datafile.addVar("time_after_ascent", ncDouble, dim_vector));
            var_vector.emplace_back(datafile.addVar("time", ncDouble, time_dim));
            var_vector.emplace_back(datafile.addVar("conv_400", ncByte, dim_vector));
            var_vector.emplace_back(datafile.addVar("conv_600", ncByte, dim_vector));
            var_vector.emplace_back(datafile.addVar("slan_400", ncByte, dim_vector));
            var_vector.emplace_back(datafile.addVar("slan_600", ncByte, dim_vector));
            var_vector.emplace_back(datafile.addVar("type", ncString, dim_vector));
            var_vector.emplace_back(datafile.addVar("lat", ncDouble, dim_vector));
            var_vector.emplace_back(datafile.addVar("lon", ncDouble, dim_vector));
            var_vector.emplace_back(datafile.addVar("step", ncUint64, dim_vector));

            for(auto &var: var_vector)
                var.setCompression(
                    enableShuffleFilter, enableDeflateFilter, deflateLevel);

            NcFile in_datafile(in_filename, NcFile::read);

            // global attributes
            auto attributes = in_datafile.getAtts();
            for(auto & name_attr: attributes)
            {
                auto attribute = name_attr.second;
                NcType type = attribute.getType();
                if(type.getName() == "double")
                {
                    std::vector<float> values(1);
                    attribute.getValues(values.data());
                    datafile.putAtt(attribute.getName(), type, values[0]);
                } else if(type.getName() == "int64"
                    || type.getName() == "int32" || type.getName() == "int")
                {
                    std::vector<int> values(1);
                    attribute.getValues(values.data());
                    datafile.putAtt(attribute.getName(), type, values[0]);
                } else if(type.getName() == "char")
                {
                    std::string values;
                    attribute.getValues(values);
                    datafile.putAtt(attribute.getName(), values);
                }
            }
            // column attributes
            // Output Parameter is new. Hence we add it like that:
            std::string att_name = "long_name";
            std::string att_val = "gradients are calculated w.r.t. this output parameter";
            param_var.putAtt(
                att_name,
                att_val);
            att_name = "standard_name";
            att_val = "output_parameter";
            param_var.putAtt(
                att_name,
                att_val);
            // Same for step
            uint32_t offset = num_comp+num_par;
            NcVar *var = &var_vector[offset+9];
            att_val = "step";
            var->putAtt(
                att_name,
                att_val);
            att_name = "long_name";
            att_val = "simulation step";
            var->putAtt(
                att_name,
                att_val);
            att_name = "auxiliary_data";
            att_val = "yes";
            var->putAtt(att_name, att_val);
            // and hail N
            att_name = "long_name";
            att_val = "specific hail number";
            var_vector[Nh_idx].putAtt(
                att_name,
                att_val);
            att_name = "standard_name";
            att_val = "specific_number_of_hail_in_air";
            var_vector[Nh_idx].putAtt(
                att_name,
                att_val);
            att_name = "units";
            att_val = "kg^-1";
            var_vector[Nh_idx].putAtt(
                att_name,
                att_val);
            att_name = "auxiliary_data";
            att_val = "yes";
            var_vector[Nh_idx].putAtt(att_name, att_val);
            // hail out
            att_name = "long_name";
            att_val = "sedimentation of specific hail number";
            var_vector[Nh_out_idx].putAtt(
                att_name,
                att_val);
            att_name = "standard_name";
            att_val = "sedi_outflux_of_hail_number";
            var_vector[Nh_out_idx].putAtt(
                att_name,
                att_val);
            att_name = "units";
            att_val = "kg^-1 s^-1";
            var_vector[Nh_out_idx].putAtt(
                att_name,
                att_val);
            att_name = "auxiliary_data";
            att_val = "yes";
            var_vector[Nh_out_idx].putAtt(att_name, att_val);

            // hail q
            att_name = "long_name";
            att_val = "specific hail content";
            var_vector[qh_idx].putAtt(
                att_name,
                att_val);
            att_name = "standard_name";
            att_val = "mass_fraction_of_hail_in_air";
            var_vector[qh_idx].putAtt(
                att_name,
                att_val);
            att_name = "units";
            att_val = "kg kg^-1";
            var_vector[qh_idx].putAtt(
                att_name,
                att_val);
            att_name = "auxiliary_data";
            att_val = "yes";
            var_vector[qh_idx].putAtt(att_name, att_val);

            // hail q out
            att_name = "long_name";
            att_val = "sedimentation of hail mixing ratio";
            var_vector[qh_out_idx].putAtt(
                att_name,
                att_val);
            att_name = "standard_name";
            att_val = "sedi_outflux_of_hail";
            var_vector[qh_out_idx].putAtt(
                att_name,
                att_val);
            att_name = "units";
            att_val = "kg kg^-1 s^-1";
            var_vector[qh_out_idx].putAtt(
                att_name,
                att_val);
            att_name = "auxiliary_data";
            att_val = "yes";
            var_vector[qh_out_idx].putAtt(att_name, att_val);

            auto vars = in_datafile.getVars();
            for(auto & name_var: vars)
            {
                auto column_name = name_var.first;
                auto it = std::find(output_par_idx.begin(), output_par_idx.end(), column_name);
                if(it != output_par_idx.end())
                {
                    uint32_t idx = std::distance(output_par_idx.begin(), it);
                    var = &var_vector[idx];
                }
                else if(column_name == "time_after_ascent")
                    var = &var_vector[offset];
                else if(column_name == "time")
                    var = &var_vector[offset+1];
                else if(column_name == "conv_400")
                    var = &var_vector[offset+2];
                else if(column_name == "conv_600")
                    var = &var_vector[offset+3];
                else if(column_name == "slan_400")
                    var = &var_vector[offset+4];
                else if(column_name == "slan_600")
                    var = &var_vector[offset+5];
                else if(column_name == "type")
                    var = &var_vector[offset+6];
                else if(column_name == "lat")
                    var = &var_vector[offset+7];
                else if(column_name == "lon")
                    var = &var_vector[offset+8];
                else if(column_name == "ensemble")
                    var = &ens_var;
                else if(column_name == "trajectory")
                    var = &traj_var;
                else if(column_name == "Q_TURBULENCE")
#if defined MET3D && defined TURBULENCE
                    var = &var_vector[offset+10];
#else
                    continue;
#endif
                else
                    continue;
                // add auxiliary column to nearly all columns by default
                if(column_name != "time" && column_name != "time_after_ascent"
                    && column_name != "lon" && column_name != "lat"
                    && column_name != "pressure" && column_name != "trajectory"
                    && column_name != "ensemble")
                {
                    att_name = "auxiliary_data";
                    att_val = "yes";
                    var->putAtt(att_name, att_val);
                }
                if(column_name == "ensemble")
                {
                    att_name = "standard_name";
                    att_val = "ensemble";
                    var->putAtt(att_name, att_val);
                    att_name = "long_name";
                    att_val = "ensemble id that is only consistent within a given history via trajectory id";
                    var->putAtt(att_name, att_val);
                    continue;
                }
                if(column_name == "trajectory")
                {
                    att_name = "standard_name";
                    att_val = "trajectory";
                    var->putAtt(att_name, att_val);
                    att_name = "long_name";
                    att_val = "trajectory";
                    var->putAtt(att_name, att_val);
                    att_name = "description";
                    att_val = "last number is the id of the instance that ran "
                        "the simulation and the numbers before are the history";
                    var->putAtt(att_name, att_val);
                    continue;
                }
                auto attributes = name_var.second.getAtts();
                for(auto & name_attr: attributes)
                {
                    // ie long_name, standard_name
                    auto attribute = name_attr.second;
                    if(attribute.getName() == "_FillValue")
                        continue;
                    NcType type = attribute.getType();
                    if(type.getName() == "double" || type.getName() == "float")
                    {
                        std::vector<double> values(1);
                        var->putAtt(
                            attribute.getName(),
                            type,
                            values.size(),
                            values.data());
                    } else if(type.getName() == "int64"
                        || type.getName() == "int32" || type.getName() == "int")
                    {
                        std::vector<int> values(1);
                        var->putAtt(
                            attribute.getName(),
                            type,
                            values.size(),
                            values.data());
                    } else if(type.getName() == "uint64"
                        || type.getName() == "uint32" || type.getName() == "uint")
                    {
                        std::vector<uint64_t> values(1);
                        var->putAtt(
                            attribute.getName(),
                            type,
                            values.size(),
                            values.data());
                    } else if(type.getName() == "char")
                    {
                        std::string values;
                        attribute.getValues(values);
                        var->putAtt(
                            attribute.getName(),
                            values);
                    }
                }
            }
            // all the gradients are auxiliary data
            for(uint64_t j=0; j<num_par; j++)
            {
                att_name = "auxiliary data";
                att_val = "yes";
                var_vector[num_comp+j].putAtt(att_name, att_val);
                // add descriptions to all gradients
                att_name = "standard_name";
                att_val = output_grad_idx[j].substr(1);
                var_vector[num_comp+j].putAtt(att_name, att_val);
                att_name = "long_name";
                att_val = output_grad_descr[j];
                var_vector[num_comp+j].putAtt(att_name, att_val);
            }


            in_datafile.close();
        } else
        {
            // write attributes
#ifdef MET3D
            std::ofstream outfile_att;
            outfile_att.open(filename + "_attributes.txt");
            outfile_att.precision(10);
            if( !outfile_att.is_open() )
            {
                std::cerr << "ERROR while opening the attribute file:\n"
                        << filename << "_attributes.txt\n. Aborting.\n";
                exit(EXIT_FAILURE);
            }
            // Global attributes
            outfile_att << "[Global attributes]\n";
            NcFile datafile(in_filename, NcFile::read);
            auto attributes = datafile.getAtts(); // <string, NcGroupAtt>
            for(auto & name_attr: attributes)
            {
                NcGroupAtt attribute = name_attr.second;
                NcType type = attribute.getType();
                if(type.getName() == "double")
                {
                    std::vector<float> values(1);
                    attribute.getValues(values.data());
                    outfile_att << "name=" << attribute.getName() << "\n"
                                << "type=" << type.getName() << "\n"
                                << "values=" << values[0] << "\n";
                } else if(type.getName() == "int64"
                    || type.getName() == "int32" || type.getName() == "int")
                {
                    std::vector<int> values(1);
                    attribute.getValues(values.data());
                    outfile_att << "name=" << attribute.getName() << "\n"
                                << "type=" << type.getName() << "\n"
                                << "values=" << values[0] << "\n";
                } else if(type.getName() == "char")
                {
                    std::string values;
                    attribute.getValues(values);
                    outfile_att << "name=" << attribute.getName() << "\n"
                            << "type=" << type.getName() << "\n"
                            << "values=" << values << "\n";
                }
            }

            // Column attributes
            outfile_att << "[Non global attributes]\n";
            auto vars = datafile.getVars();
            for(auto & name_var: vars)
            {
                auto var = name_var.second;
                auto name = name_var.first;
                outfile_att << "column=" << name << "\n";
                auto attributes = var.getAtts();
                for(auto & name_attr: attributes)
                {
                    auto attribute = name_attr.second;
                    NcType type = attribute.getType();
                    if(type.getName() == "double" || type.getName() == "float")
                    {
                        std::vector<double> values(1);
                        attribute.getValues(values.data());
                        outfile_att << attribute.getName() << "=" << values[0] << "\n";
                    } else if(type.getName() == "int64"
                        || type.getName() == "int32" || type.getName() == "int")
                    {
                        std::vector<int> values(1);
                        attribute.getValues(values.data());
                        outfile_att << attribute.getName() << "=" << values[0] << "\n";
                    } else if(type.getName() == "char")
                    {
                        std::string values;
                        attribute.getValues(values);
                        outfile_att << attribute.getName() << "=" << values << "\n";
                    }
                }
            }
            outfile_att.close();
#endif
            // write reference quantities
            std::ofstream outfile_refs;
            outfile_refs.open(filename + "_reference_values.txt");
            outfile_refs.precision(10);

            if( !outfile_refs.is_open() )
            {
                std::cerr << "ERROR while opening the reference file. Aborting." << std::endl;
                exit(EXIT_FAILURE);
            }

            // Write the reference quantities
            outfile_refs << ref_quant.Tref << " "
                << ref_quant.pref << " "
                << ref_quant.qref << " "
                << ref_quant.Nref << " "
                << ref_quant.wref << " "
                << ref_quant.tref << " "
                << ref_quant.zref << "\n";

            outfile_refs.close();

            // write headers
            std::string suffix = ".txt";
            std::string full_filename;
            full_filename = filename;
            full_filename += suffix;

            outfile.open(full_filename);
            outfile.precision(10);

            if( !outfile.is_open() )
            {
                std::cerr << "ERROR while opening the outputfile. Aborting." << std::endl;
                exit(EXIT_FAILURE);
            }

            // Append the initial values and write headers
            out_tmp << std::setprecision(10) << "step,trajectory,lon,lat,"
#if defined WCB
                << "MAP,";
#endif
#if defined WCB2
                << "WCB_flag,"
                << "dp2h,"
#endif
#if defined WCB2 || defined MET3D
                << "conv_400,"
                << "conv_600,"
                << "slan_400,"
                << "slan_600,";
#endif
#if defined MET3D
            out_tmp
                << "time,"
                << "time_after_ascent,"
                << "type,"
                << "ensemble,"
                << "instance_id,";
#endif
            for(uint32_t i=0; i<output_par_idx.size(); ++i)
                out_tmp << output_par_idx[i]  <<
                    ((i < output_par_idx.size()-1) ? "," : "\n");

            std::string basename = "_diff_";
            std::string fname;

            for(int ii = 0; ii < num_comp; ii++)
            {
                fname = filename;
                fname += basename;
                fname += std::to_string(ii);
                fname += suffix;

                out_diff[ii].open(fname);

                if( !out_diff[ii].is_open() )
                {
                    std::cerr << "ERROR while opening outputfile. Aborting." << std::endl;
                    exit(EXIT_FAILURE);
                }
                out_diff[ii].precision(10);
                out_diff_tmp[ii]
                    << std::setprecision(10)
                    << "step,"
                    << "trajectory,"
                    << "Output Parameter,"
                    << "lon,"
                    << "lat,"
#if defined WCB
                    << "MAP,";
#endif
#if defined WCB2
                    << "WCB_flag,"
                    << "dp2h,"
#endif
#if defined WCB2 || defined MET3D
                    << "conv_400,"
                    << "conv_600,"
                    << "slan_400,"
                    << "slan_600,";
#endif
#if defined MET3D
                out_diff_tmp[ii]
                    << "time,"
                    << "time_after_ascent,"
                    << "type,"
                    << "ensemble,"
                    << "instance_id,";
#endif
                for(uint32_t i=0; i<output_grad_idx.size(); ++i)
                    out_diff_tmp[ii] << output_grad_idx[i]  <<
                        ((i < output_grad_idx.size()-1) ? "," : "\n");
            } // End loop over all components
        }
    }

    /**
     * Writes data either to a stringstream for txt files or to different vectors
     * for netCDF files.
     */
    void buffer(const model_constants_t &cc,
        const nc_parameters_t &nc_params,
        const std::vector<codi::RealReverse> &y_single_new,
        const std::vector< std::array<double, num_par > >  &y_diff,
        const uint32_t sub,
        const uint32_t t,
        const double time_new,
        const uint32_t traj_id,
        const uint32_t ensemble,
        const reference_quantities_t &ref_quant)
    {
        if(filetype == "netcdf")
        {
            const uint64_t offset = n_trajs*n_ens;

            // output parameters
            for(uint64_t i=0; i<num_comp; i++)
            {
                switch(i)
                {
                    case p_idx:
                        output_buffer[i][n_snapshots*offset] =
                            y_single_new[i].getValue() * ref_quant.pref;
                        break;
                    case T_idx:
                        output_buffer[i][n_snapshots*offset] =
                            y_single_new[i].getValue() * ref_quant.Tref;
                        break;
                    case w_idx:
                        output_buffer[i][n_snapshots*offset] =
                            y_single_new[i].getValue() * ref_quant.wref;
                        break;
                    case z_idx:
                        output_buffer[i][n_snapshots*offset] =
                            y_single_new[i].getValue() * ref_quant.zref;
                        break;
                    case qc_idx:
                    case qr_idx:
                    case qv_idx:
                    case qi_idx:
                    case qs_idx:
                    case qg_idx:
                    case qh_idx:
                    case qi_out_idx:
                    case qs_out_idx:
                    case qr_out_idx:
                    case qg_out_idx:
                    case qh_out_idx:
                        output_buffer[i][n_snapshots*offset] =
                            y_single_new[i].getValue() * ref_quant.qref;
                        break;
                    case Nc_idx:
                    case Nr_idx:
                    case Ni_idx:
                    case Ns_idx:
                    case Ng_idx:
                    case Nh_idx:
                    case Ni_out_idx:
                    case Ns_out_idx:
                    case Nr_out_idx:
                    case Ng_out_idx:
                    case Nh_out_idx:
                        output_buffer[i][n_snapshots*offset] =
                            y_single_new[i].getValue() * ref_quant.Nref;
                        break;
                    default:
                        output_buffer[i][n_snapshots*offset] =
                            y_single_new[i].getValue();
                        break;
                }
            }

            // gradients
            for(uint64_t i=0; i<num_comp; i++) // gradient sensitive to output parameter i
                for(uint64_t j=0; j<num_par; j++) // gradient of input parameter j
                    output_buffer[num_comp+j][i + n_snapshots*offset*num_comp] = y_diff[i][j];

            // time after ascent
            output_buffer[num_comp+num_par][n_snapshots*offset] =
                nc_params.time_rel + sub*cc.dt;

            // time index
            output_buffer[num_comp+num_par+1][n_snapshots] =
                nc_params.time_abs[t + nc_params.time_idx] + sub*cc.dt;

            // lat
            output_buffer[num_comp+num_par+2][n_snapshots*offset] =
                (nc_params.lat[0] + sub*nc_params.dlat);

            // lon
            output_buffer[num_comp+num_par+3][n_snapshots*offset] =
                (nc_params.lon[0] + sub*nc_params.dlon);

            // flags
            output_buffer_flags[0][n_snapshots*offset] = nc_params.conv_400;
            output_buffer_flags[1][n_snapshots*offset] = nc_params.conv_600;
            output_buffer_flags[2][n_snapshots*offset] = nc_params.slan_400;
            output_buffer_flags[3][n_snapshots*offset] = nc_params.slan_600;

            // type
            output_buffer_str[0][n_snapshots*offset] = nc_params.type[0];

            // simulation step
            output_buffer_int[0][n_snapshots*offset] = sub + t*cc.num_sub_steps;
        } else
        {
#if defined WCB || defined WCB2
            out_tmp << time_new << "," << traj_id << ","
                    << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                    << (nc_params.lat[0] + sub*nc_params.dlat) << ","
                    << nc_params.ascent_flag << ",";
#elif defined MET3D
            out_tmp << sub + t*cc.num_sub_steps << "," << traj_id << ","
                    << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                    << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#else
            out_tmp << time_new << "," << traj_id << ","
                    << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                    << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#endif
#if defined WCB2 || defined MET3D
            out_tmp
#if defined WCB2
                    << nc_params.dp2h << ","
#endif
                    << nc_params.conv_400 << ","
                    << nc_params.conv_600 << ","
                    << nc_params.slan_400 << ","
                    << nc_params.slan_600 << ",";
#endif
#ifdef MET3D
            out_tmp << nc_params.time_abs[t + nc_params.time_idx] + sub*cc.dt << ","
                    << nc_params.time_rel + sub*cc.dt << ","
                    << nc_params.type[0] << ","
                    << ensemble << ","
                    << cc.id << ",";
#endif
            for(int ii = 0 ; ii < num_comp; ii++)
                out_tmp << y_single_new[ii]
                    << ((ii == num_comp-1) ? "\n" : ",");

            for(int ii = 0 ; ii < num_comp ; ii++)
            {
#if defined WCB || defined WCB2
                out_diff_tmp[ii] << time_new << "," << traj_id << ","
                                << output_par_idx[ii] << ","
                                << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                                << (nc_params.lat[0] + sub*nc_params.dlat) << ","
                                << nc_params.ascent_flag << ",";
#elif defined MET3D
                out_diff_tmp[ii] << sub + t*cc.num_sub_steps << "," << traj_id << ","
                                << output_par_idx[ii] << ","
                                << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                                << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#else
                out_diff_tmp[ii] << time_new << "," << traj_id << ","
                                << output_par_idx[ii] << ","
                                << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                                << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#endif
#if defined WCB2 || defined MET3D
                out_diff_tmp[ii]
#if defined WCB2
                    << nc_params.dp2h << ","
#endif
                    << nc_params.conv_400 << ","
                    << nc_params.conv_600 << ","
                    << nc_params.slan_400 << ","
                    << nc_params.slan_600 << ",";
#endif
#if defined MET3D
                out_diff_tmp[ii] << nc_params.time_abs[t + nc_params.time_idx] + sub*cc.dt << ","
                            << nc_params.time_rel + sub*cc.dt << ","
                            << nc_params.type[0] << ","
                            << ensemble << ","
                            << cc.id << ",";
#endif
                for(int jj = 0 ; jj < num_par ; jj++)
                    out_diff_tmp[ii] << y_diff[ii][jj]
                        << ((jj==num_par-1) ? "\n" : ",");
            }
        }
        n_snapshots++;
    }

    /**
     * Write the buffered data to disk.
     */
    void flush_buffer()
    {
        if(filetype == "netcdf")
        {
            std::vector<uint64_t> startp, countp;
            startp.push_back(flushed_snapshots);
            countp.push_back(n_snapshots); // number of snapshots so far
            // time index
            var_vector[num_comp+num_par+1].putVar(startp, countp,
                output_buffer[num_comp+num_par+1].data());

            startp.push_back(0);
            startp.push_back(0);
            countp.push_back(n_ens);
            countp.push_back(n_trajs);

            for(uint64_t i=0; i<num_comp; i++)
                var_vector[i].putVar(startp, countp,
                    output_buffer[i].data());

            uint64_t offset = num_comp+num_par;
            // time after ascent
            var_vector[offset].putVar(startp, countp,
                output_buffer[num_comp+num_par].data());

            offset += 2;
            // flags
            for(uint64_t i=0; i<output_buffer_flags.size(); i++)
                var_vector[offset+i].putVar(startp, countp,
                output_buffer_flags[i].data());

            offset += output_buffer_flags.size();
            // type
            std::vector<uint64_t> startp_str, countp_str;
            for(auto &p: startp)
                startp_str.push_back(p);
            for(auto &p: countp)
                countp_str.push_back(p);
            countp_str[0] = 1;

            // std::cout << "\nstartp_str: ";
            // for(auto &p: startp_str) std::cout << p << ", ";
            // std::cout << "\ncountp_str: ";
            // for(auto &p: countp_str) std::cout << p << ", ";
            // std::cout << "\nn_snapshots " << n_snapshots << "\n";



            for(uint64_t i=0; i<output_buffer_str.size(); i++)
            {
                // startp_str[0] = 0;
                // write one string at a time.
                for(const auto &t: output_buffer_str[i])
                {
                    var_vector[offset+i].putVar(
                        startp_str, countp_str, &t);
                    startp_str[0]++;
                    if(startp_str[0]-startp[0] == n_snapshots)
                        break;
                }
            }

            offset += output_buffer_str.size();
            // lat
            var_vector[offset].putVar(startp, countp,
                output_buffer[num_comp+num_par+2].data());

            offset += 1;
            // lon
            var_vector[offset].putVar(startp, countp,
                output_buffer[num_comp+num_par+3].data());

            offset += 1;
            // step
            var_vector[offset].putVar(startp, countp,
                output_buffer_int[0].data());

            offset = num_comp;
            // gradients
            startp.push_back(0);
            countp.push_back(num_comp);

            for(uint64_t j=0; j<num_par; j++)
                var_vector[offset+j].putVar(startp, countp,
                    output_buffer[num_comp+j].data());
        } else
        {
            outfile << out_tmp.rdbuf();
            for(int ii = 0 ; ii < num_comp ; ii++)
            {
                out_diff[ii] << out_diff_tmp[ii].rdbuf();
                out_diff_tmp[ii].str( std::string() );
                out_diff_tmp[ii].clear();
            }
            out_tmp.str( std::string() );
            out_tmp.clear();
        }
        flushed_snapshots += n_snapshots;
        n_snapshots = 0;

    }
};
/** @} */ // end of group types
