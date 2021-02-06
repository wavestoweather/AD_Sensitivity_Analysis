#pragma once

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

    global_args_t();

};