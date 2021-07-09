#pragma once

#include <stdlib.h>
#include <iostream>
// #include <sys/stat.h>
#include <getopt.h>
#include "include/misc/error.h"
#include "include/microphysics/constants.h"

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

    /**
     * Use a simulation mode for sensitivity analysis, perturbance only,
     * perturbance with sensitivity analysis for trajectories and grid based data.
     */
    int simulation_mode_flag;
    char* simulation_mode_string;

    int output_flag; /*!< Output path specified? */
    char* output_string;

    int input_flag; /*!< Input netCDF file specified? */
    char* input_file;

    // int start_over_flag; /*!< Reload mixing ratios and particle numbers from trajectory every few seconds? */
    // char* start_over_string;

    int time_start_idx_flag;
    char* time_start_idx_string;

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

    int tracking_file_flag;
    char* tracking_file_string;

    // int gnu_id_flag; /*!< ID given for this instance, i.e. thread_id or id by GNU parallel. */
    // char* gnu_id_string;

    int folder_name_flag;
    char* folder_name_string;

    int n_ens_flag;
    char* n_ens_string;

    global_args_t();

    /**
     * Parse the arguments and store them in global_args.
     *
     * @param argc Number of arguments
     * @param argv Pointer to arguments
     *
     * @return Error code.
     */
    int parse_arguments(
        const int argc,
        char* const * argv,
        const int &rank,
        const int &n_processes);
    /**
     * Display a help message on how to use this program.
     */
    void display_usage();

    /**
     * Display an error message when command line arguments are faulty.
     */
    void display_error_on_command_line();
};