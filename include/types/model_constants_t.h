#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <codi.hpp>

#include "include/misc/error.h"
#include "include/types/collection_model_constants_t.h"
#include "include/types/gamma_table_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/particle_model_constants_t.h"
#include "include/microphysics/physical_parameterizations.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/table_t.h"

// namespace pt = boost::property_tree;

/**
 * Structure for constants of a model. Includes particle constants as well.
 */
template<class float_t>
struct model_constants_t {
    /**
     * Initial id of this simulation. Emerging ensembles from this
     * have other ids which are set by reading a checkpoint file.
     */
    std::string id;
    /**
     *  Running id for the ensembles starting from this instance.
     */
    uint64_t ensemble_id;
    int rank;

    /**
     * Running id for the trajectory in this ensemble.
     */
    uint64_t traj_id;

    /**
     * Maximum number of trajectories in every ensemble.
     * Needed to write proper output files.
     */
    uint64_t max_n_trajs;

    /**
     * Number of trajectories in the current ensemble.
     */
    uint64_t n_trajs;

    /**
     * Total number of ensembles.
     */
    uint64_t n_ensembles;

    /**
     * Description of the ensemble i.e.
     * "root - 0,0 perturb ccn_a1 - 1,3 perturb cloud_max_x"
     * where the first number is the ensemble id, the second number is
     * the trajectory it originates from
     */
    std::string ens_desc;

    //
    // Physical constants warm cloud
    //
    double alpha_d; /*!< Accomodation coefficient */

    /**
     * Model constants for hail.
     */
    particle_model_constants_t<float_t> hail;
    /**
     * Model constants for ice.
     */
    particle_model_constants_t<float_t> ice;
    /**
     * Model constants for snow.
     */
    particle_model_constants_t<float_t> snow;
    /**
     * Model constants for cloud.
     */
    particle_model_constants_t<float_t> cloud;
    /**
     * Model constants for rain.
     */
    particle_model_constants_t<float_t> rain;
    /**
     * Model constants for graupel.
     */
    particle_model_constants_t<float_t> graupel;

    /**
     * Snow - cloud droplet riming.
     */
    collection_model_constants_t<float_t> coeffs_scr;
    /**
     * Snow - rain droplet riming.
     */
    collection_model_constants_t<float_t> coeffs_srr;
    /**
     * Ice - rain droplet riming.
     */
    collection_model_constants_t<float_t> coeffs_irr;
    /**
     * Ice - cloud droplet riming.
     */
    collection_model_constants_t<float_t> coeffs_icr;
    /**
     * Hail - rain droplet riming.
     */
    collection_model_constants_t<float_t> coeffs_hrr;
    /**
     * Graupel - rain droplet riming.
     */
    collection_model_constants_t<float_t> coeffs_grr;
    /**
     * Hail - cloud droplet riming.
     */
    collection_model_constants_t<float_t> coeffs_hcr;
    /**
     * Graupel - cloud droplet riming.
     */
    collection_model_constants_t<float_t> coeffs_gcr;
    /**
     * Snow - ice collection.
     */
    collection_model_constants_t<float_t> coeffs_sic;
    /**
     * Hail - ice collection.
     */
    collection_model_constants_t<float_t> coeffs_hic;
    /**
     * Graupel - ice collection.
     */
    collection_model_constants_t<float_t> coeffs_gic;
    /**
     * Hail - snow collection.
     */
    collection_model_constants_t<float_t> coeffs_hsc;
    /**
     * Graupel - snow collection.
     */
    collection_model_constants_t<float_t> coeffs_gsc;

    //
    // Technical constants
    //
    double t_end_prime;       /*!< End time in seconds for the simulation. */
    double t_end;             /*!< End time for the simulation. */
    double dt_prime;          /*!< Timestep size in seconds for the simulation. */
    double dt;                /*!< Timestep size for the simulation. */
    double dt_traject_prime;  /*!< Timestep size of the trajectory from the netCDF file. */
    double dt_traject;        /*!< Timestep size of the trajectory from the netCDF file. */
    uint64_t num_steps;       /*!< Number of time steps to read from the netCDF file. */
    uint64_t done_steps;      /*!< Number of time steps already done before loading a checkpoint. */
    /**
     * Number of already done steps from the checkpoint. Needed to get the correct index for the netcdf reader.
     */
    uint64_t checkpoint_steps;
    double start_time;        /*!< Start time in seconds from the netCDF file. */

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

    float_t a1_scale; /*!< Performance constants warm cloud */
    float_t a2_scale; /*!< Performance constants warm cloud */
    float_t e1_scale; /*!< Performance constants warm cloud */
    float_t e2_scale; /*!< Performance constants warm cloud */
    float_t d_scale;  /*!< Performance constants warm cloud */

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

    // Parameters to calculate the change in buoyancy
    double poly2, poly1, poly0;
    uint32_t current_w_poly_idx;

    // See constants.h for a description of those.
    std::array<float_t, static_cast<int>(Cons_idx::n_items)> constants;

    /*
    * Uncertainty for every parameter.
    * Currently sets everything to 10% of every parameter.
    */
    std::array<double, static_cast<uint32_t>(Cons_idx::n_items)
        + static_cast<uint32_t>(Init_cons_idx::n_items)> uncertainty;
    // std::array<float_t, static_cast<uint32_t>(Init_cons_idx::n_items)> initial_conditions;
    std::array<codi::RealForwardVec<num_par_init>, static_cast<uint32_t>(Init_cons_idx::n_items)> initial_conditions;

    /**
     * Store any idx from perturbed parameters.
     */
    std::vector<uint32_t> perturbed_idx;

    /**
     * Structure to hold the new equidistant lookup table for
     * graupel wetgrowth diameter
     */
    table_t<float_t> ltabdminwgg;
    gamma_table_t table_g1, table_g2, table_r1, table_r2, table_r3;

    int local_num_comp;
    int local_num_par;
    int local_ic_par;

    explicit model_constants_t(
        const std::string &tracking_filename,
        const bool &track_initial_cond);

    /**
     * Register the model parameters for initial conditions using
     * the forward mode.
     *
     */
    void register_input();

    /**
     * Register the model parameters on the tape for codi::RealReverseIndex.
     *
     * @param tape Tape where the parameters should be registered on.
     */
    void register_input(
        codi::RealReverseIndex::Tape &tape);

    /**
     * Evaluate the gradients in forward mode. Used to track initial
     * conditions.
     */
    void get_gradients(
        std::vector<float_t> &y_single_new,
        std::vector< std::array<double, num_par > > &y_diff,
        const reference_quantities_t &ref_quant) const;

    /**
     * Register output on tape, evaluate it, get all gradients and
     * reset the tape.
     */
    void get_gradients(
        std::vector<float_t> &y_single_new,
        std::vector< std::array<double, num_par > > &y_diff,
        codi::RealReverseIndex::Tape &tape,
        const reference_quantities_t &ref_quant);

    /**
     * Get the gradients of all its members. You need to register them on a
     * tape before to get meaningful values  in the backward mode.
     *
     * @param out_vec On out: Stores all gradients.
     */
    void get_gradient(
        std::array<double, num_par> &out_vec,
        std::vector<codi::RealForwardVec<num_par_init> > &y_single_new,
        uint32_t ii,
        const double &ref_value) const;

    void get_gradient(
        std::array<double, num_par> &out_vec,
        const double &ref_value) const;


#if defined(RK4_ONE_MOMENT)
    /**
     * Set the constants for the cloud model from given environmental conditions.
     *
     * @param y Vector of initial conditions for pressure and temperature
     * @param cc Pointer to constants from the model. On out: modified constants
     * @param ref Pointer to reference quantities to transform between units
     */
    void setCoefficients(
        std::vector<float_t> & y,
        reference_quantities_t& ref);

    /**
     * Set the constants for the cloud model from given environmental conditions.
     *
     * @param p_prime Initial pressure in Pa
     * @param T_prime Initial temperature in Kelvin
     * @param cc Pointer to constants from the model. On out: modified constants
     */
    void setCoefficients(
        float_t p_prime,
        float_t T_prime);
#endif

    /**
     * Setup the cloud autoconversion parameters.
     *
     * @param pc Model constants for a certain particle type.
     */
    void setup_cloud_autoconversion(particle_model_constants_t<float_t> &pc);

    /**
     * Setup all model constants and gamma tables including dependent constants.
     * @param input
     * @param ref_quant
     */
    void setup_model_constants(
        const input_parameters_t &input,
        const reference_quantities_t &ref_quant);

    /**
     * Setup dependent variables.
     *
     * @param ref_quant
     */
    void setup_dependent_model_constants();

    /**
     * Set the uncertainty for every parameter with a default value of 10% of each parameter.
     */
    void set_uncertainty();

    /**
     * Set the uncertainty for every parameter with scale * parameter value, i.e., scale = 0.1 results in
     * 10% of each parameter.
     * @param scale
     */
    void set_uncertainty(double scale);

    /**
     * Set time step from input netcdf file.
     */
    void set_dt(const double dt_prime, const reference_quantities_t &ref_quant);

    /**
     * Print all parameters that are being tracked with algorithmic
     * differentiation.
     */
    void print(std::ostream &os = std::cout);

    /**
     * Check if a certain model state variable or model parameter
     * should be written to the output.
     *
     * @param idx Index of the model state variable or model parameter
     *              to check.
     * @param state_param   0: idx is a model parameter
     *                      1: idx is a model state variable
     *                      2: idx is an initial condition
     */
    bool trace_check(const int &idx, const int state_param) const;

    /**
     * Define which data shall be written by loading a configuration file.
     *
     * @param filename Path to json configuration file
     */
    void load_configuration(const std::string &filename);

    /**
     * Return the values of all perturbed parameters.
     *
     * @return Vector with perturbed values of model parameters.
     */
#ifdef OUT_DOUBLE
    void get_perturbed_info(
        std::vector<double> &p_values,
        std::vector<uint64_t> &p_idx) const;
#else
    void get_perturbed_info(
        std::vector<float> &perturbed,
        std::vector<uint64_t> &param_idx) const;
#endif

    /**
     * Given three datapoints for the ascent velocity, fit
     * a polynomial of degree 2. This is needed if the time step of the simulation
     * is smaller than what is provided in the dataset.
     *
     * @param w0 Velocity [m/s] at the current time step.
     * @param w1 Velocity [m/s] at the next valid time step.
     * @param w2 Velocity [m/s] at the next valid time step after w1.
     * @param offset_1 Time [s] between w0 and w1.
     * @param offset_2 Time [s] between w1 and w2.
     */
    void set_dw_polynomial(
        const double &w0,
        const double &w1,
        const double &w2,
        const int &offset_1,
        const int &offset_2);

    /**
     * Get the current change in ascent velocity.
     *
     * @return Change in ascent velocity [m/s^2]
     */
    double get_dw() const;

    /**
     * Get the current ascent velocity.
     *
     * @return Ascent velocity [m/s]
     */
    double get_w() const;

    /**
     * Needed between each sub step of the simulation to get the interpolated ascent velocity
     * between two time steps of the provided dataset.
     * This sets the current index that is used in get_dw() and get_w().
     */
    void increment_w_idx();

    /**
     * Used to switch on or off certain trackings.
     */
    uint64_t track_state;
    uint64_t track_ic;
    std::vector<uint64_t> track_param;

    int from_json(const nlohmann::json &j);
};

template<class float_t>
void to_json(nlohmann::json& j, const model_constants_t<float_t>& cc);
