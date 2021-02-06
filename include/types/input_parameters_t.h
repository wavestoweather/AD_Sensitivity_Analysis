#pragma once

#include <boost/property_tree/ptree.hpp>
#include <cmath>
#include <string>

#include "include/misc/error.h"

namespace pt = boost::property_tree;

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

    input_parameters_t();

    /**
     * Change the output filename such that it starts with "idx-y-z_"
     * where x-y are the ids of preceding trajectories and z is the current id
     *
     * @params all_ids String of form x-y-z with x,y and z positive ints
     *                 and z being the current id.
     */
    void set_outputfile_id(
        const std::string &all_ids,
        uint64_t ensemble_id);

    void put(pt::ptree &ptree, const double &time) const;

    void put(pt::ptree &ptree) const;

    /**
     * Set values from property tree used in reading checkpoint files.
     */
    int from_pt(pt::ptree &ptree);
};