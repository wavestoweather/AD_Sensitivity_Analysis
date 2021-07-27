#pragma once

#include <cmath>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include "include/misc/error.h"
#include "include/types/global_args_t.h"

namespace pt = boost::property_tree;

/**
 * Structure to collect all input parameters.
 */
struct input_parameters_t {
    // Numerics
    double t_end_prime; /*!< End simulation time in seconds. */
    double dt_prime; /*!< Timestep size in seconds for the simulation. */
    double dt_traject_prime; /*!< Timestep size in seconds of the trajectory in the netCDF file. */
    double dt_traject; /*!< Timestep size of the trajectory in the netCDF file. */
    int32_t start_time_idx; /*!< Timestep index to start (multiply by dt_traject_prime to get the seconds). */
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
    std::string tracking_filename; /*!< File to switch on specific variables to track. */
    uint32_t id; /*!< ID given for this instance, i.e. thread_id or id by GNU parallel. */
    uint32_t n_ensembles;

    // bool start_over; /*!< Start over at new timestep of trajectory? */
    bool start_over_env; /*!< Start over environment variables at new timestep of trajectory? */
    bool fixed_iteration; /*!< Fix temperature and pressure at every iteration? */

    uint32_t auto_type; /*!< Particle type. */
    uint32_t traj; /*!< Trajectory index to load from the netCDF file. */
    uint32_t write_index; /*!< Write stringstream every x iterations to disk. */
    uint64_t progress_index; /*!< Index for updating progressbar. */
    uint32_t ensemble; /*!< Index of ensemble. */
    double current_time; /*!< Time for and from checkpoint files. */

    int simulation_mode; /*!< Simulation mode. */

    input_parameters_t();

    void put(pt::ptree &ptree, const double &time) const;

    void put(pt::ptree &ptree) const;

    /**
     * Set values from property tree used in reading checkpoint files.
     */
    int from_pt(pt::ptree &ptree);

    /**
     * Set the input parameters with the data from the global arguments.
     *
     * @param arg Stuct with command line arguments.
     * @param in Struct where the input parameters will be stored.
     */
    void set_input_from_arguments(global_args_t &arg);

    /**
     * Print all given input parameters.
     */
    void print_parameters();
};
