#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "codi.hpp"

#include "include/misc/error.h"
#include "include/types/collection_model_constants_t.h"
#include "include/types/gamma_table_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/particle_model_constants_t.h"
#include "include/microphysics/physical_parameterizations.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/table_t.h"

namespace pt = boost::property_tree;

/**
 * Structure for constants of a model. Includes particle constants as well.
 */
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
    uint64_t num_steps;       /*!< Number of time steps to read from the netCDF file. */
    uint64_t done_steps;      /*!< Number of time steps already done before loading a checkpoint. */
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

    // See constants.h for a description of those.
    std::vector<codi::RealReverse> constants;

    /*
    * Uncertainty for every parameter.
    * Currently sets everything to 10% of every parameter.
    */
    std::array<double, static_cast<uint32_t>(Cons_idx::n_items)> uncertainty;

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

    int local_num_comp;
    int local_num_par;

    explicit model_constants_t(const std::string &tracking_filename);

    /**
     * Register the model parameters on the tape for codi::RealReverse.
     *
     * @param tape Tape where the parameters should be registered on.
     */
    void register_input(
        codi::RealReverse::TapeType &tape);

    /**
     * Register output on tape, evaluate it, get all gradients and
     * reset the tape.
     */
    void get_gradients(
        std::vector<codi::RealReverse> &y_single_new,
        std::vector< std::array<double, num_par > > &y_diff,
        codi::RealReverse::TapeType &tape) const;

    /**
     * Get the gradients of all its members. You need to register them on a
     * tape before to get meaningful values.
     *
     * @param out_vec On out: Stores all gradients.
     */
    void get_gradient(
        std::array<double, num_par> &out_vec, uint32_t ii) const;

    /**
     * Put any perturbed parameter to a property tree.
     * This will compare the parameters to the constants
     * available in constants.h, assuming this is only called during
     * checkpoint writing.
     *
     * @params ptree Property tree, where a tree "model_constants" is being added.
     */
    void put(pt::ptree &ptree) const;

    int from_pt(pt::ptree &ptree);

    /**
     * Set the constants for the cloud model from given environmental conditions.
     *
     * @param y Vector of initial conditions for pressure and temperature
     * @param cc Pointer to constants from the model. On out: modified constants
     * @param ref Pointer to reference quantities to transform between units
     */
    void setCoefficients(
        std::vector<codi::RealReverse> & y,
        reference_quantities_t& ref);

    /**
     * Set the constants for the cloud model from given environmental conditions.
     *
     * @param p_prime Initial pressure in Pa
     * @param T_prime Initial temperature in Kelvin
     * @param cc Pointer to constants from the model. On out: modified constants
     */
    void setCoefficients(
        codi::RealReverse p_prime,
        codi::RealReverse T_prime);

    /**
     * Setup the cloud autoconversion parameters.
     *
     * @param pc Model constants for a certain particle type.
     */
    void setup_cloud_autoconversion(particle_model_constants_t &pc);

    /**
     * Setup model constants and gamma tables.
     */
    void setup_model_constants(
        const input_parameters_t &input,
        const reference_quantities_t &ref_quant);

    /**
     * Print all parameters that are being tracked with algorithmic
     * differentiation.
     */
    void print();

    /**
     * Check if a certain model state variable or model parameter
     * should be written to the output.
     *
     * @param idx Index of the model state variable or model parameter
     *              to check.
     * @param state_param If true: idx is a model state variable
     */
    bool trace_check(const int &idx, const bool state_param) const;

    /**
     * Define which data shall be written by loading a configuration file.
     *
     * @param filename Path to json configuration file
     */
    void load_configuration(const std::string &filename);

//  private:
    /**
     * Used to switch on or off certain trackings.
     */
    uint64_t track_state;
    std::vector<uint64_t> track_param;
};
