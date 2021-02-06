#include "include/types/global_args_t.h"

global_args_t::global_args_t()
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