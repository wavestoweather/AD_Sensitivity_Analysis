#include "include/types/input_parameters_t.h"

input_parameters_t::input_parameters_t()
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

void input_parameters_t::set_outputfile_id(
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

void input_parameters_t::put(
    pt::ptree &ptree,
    const double &time) const
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

void input_parameters_t::put(
    pt::ptree &ptree) const
{
    put(ptree, current_time);
}

int input_parameters_t::from_pt(
    pt::ptree &ptree)
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

void input_parameters_t::set_input_from_arguments(
    global_args_t &arg)
{
    // Final time
    if(1 == arg.final_time_flag){
        this->t_end_prime = std::strtod(arg.final_time_string, nullptr);
    }

    // Timestep
    if(1 == arg.timestep_flag){
        this->dt_prime = std::strtod(arg.timestep_string, nullptr);
    }

    // Snapshot index
    if(1 == arg.snapshot_index_flag){
        this->snapshot_index = std::stoi(arg.snapshot_index_string);
    }

    // Output
    if(1 == arg.output_flag){
        this->OUTPUT_FILENAME = arg.output_string;
    }

    // Input
    if(1 == arg.input_flag){
        this->INPUT_FILENAME = arg.input_file;
    }

    // Scaling factor
    if(1 == arg.scaling_fact_flag){
        this->scaling_fact = std::strtod(arg.scaling_fact_string, nullptr);
    }

    // Starting over mixing ratios and particle numbers
    if(1 == arg.start_over_flag){
        this->start_over = (strcmp(arg.start_over_string, "0"));
    }

    // Starting over environment variables (p, T, w)
    if(1 == arg.start_over_env_flag){
        this->start_over_env = (strcmp(arg.start_over_env_string, "0"));
    }

    if(1 == arg.fixed_iteration_flag){
        this->fixed_iteration = (strcmp(arg.fixed_iteration_string, "0"));
    }

    // Auto type
    if(1 == arg.auto_type_flag){
        this->auto_type = std::stoi(arg.auto_type_string);
    }

    // Trajectory
    if(1 == arg.traj_flag){
        this->traj = std::stoi(arg.traj_string);
    }

    // Write index
    if(1 == arg.write_flag){
        this->write_index = std::stoi(arg.write_string);
    }

    // Progressbar index
    if(1 == arg.progress_index_flag){
        this->progress_index = std::stoull(arg.progress_index_string);
    }
#ifdef MET3D
    // Simulation start time
    if(1 == arg.delay_start_flag){
        this->start_time = std::strtod(arg.delay_start_string, nullptr);
    }
#endif

    // Ensemble configuration file
    if(1 == arg.ens_config_flag){
        this->ENS_CONFIG_FILENAME = arg.ens_config_string;
    }

    // Checkpoint file
    if(1 == arg.checkpoint_flag){
        this->CHECKPOINT_FILENAME = arg.checkpoint_string;
    }

    // ID for this process
    if(1 == arg.gnu_id_flag){
        this->id = std::stoi(arg.gnu_id_string);
    }

    // Folder name for new generated checkpoints
    if(1 == arg.folder_name_flag){
        this->FOLDER_NAME = arg.folder_name_string;
    }
}

void input_parameters_t::print_parameters()
{
#ifdef SILENT_MODE
    return;
#endif
    std::cout << "\n"
        << "Technical input parameters:\n"
        << "---------------------------\n"
        << "Time to integrate: " << this->t_end_prime << " Seconds\n"
        << "Timestep: " << this->dt_prime << " Seconds\n"
#ifdef MET3D
        << "Start time (relative to ascend): " << this->start_time << " Seconds" << "\n"
#endif
        << ( (this->CHECKPOINT_FILENAME != "")
        ?   "Start time from checkpoint (relative to ascend): " + std::to_string(this->current_time + this->start_time) + "\n"
        :   "" )
        << "Snapshot index: " << this->snapshot_index << "\n"
        << "Write index: " << this->write_index << "\n"
        << "Progressbar index: " << this->progress_index << "\n"
        << "Name of output file: " << this->OUTPUT_FILENAME << "\n"
        << "Scaling factor: " << this->scaling_fact << "\n"
        << "Name of input file: " << this->INPUT_FILENAME << "\n"
        << "Start over mixing ratios and particle numbers at each timestep of a trajectory?: " << this->start_over << "\n"
        << "Start over pressure, temperature and ascent at each timestep of a trajectory?: " << this->start_over_env << "\n"
        << "Fix temperature and pressure at each substep?: " << this->fixed_iteration << "\n"
        << "Auto type for rain evaporation (1, 2, 3): " << this->auto_type << "\n"
        << "Trajectory used: " << this->traj << "\n"
        << "Instance id: " << this->id << "\n"
        << ( (this->ENS_CONFIG_FILENAME != "")
        ?   "Ensemble configuration file: " + this->ENS_CONFIG_FILENAME + "\n"
        :   "" )
        << ( (this->CHECKPOINT_FILENAME != "")
        ?   "Checkpoint file: " + this->CHECKPOINT_FILENAME + "\n"
        :   "" )
        << "Folder name for newly generated checkpoints: " << this->FOLDER_NAME << "\n"
        << std::endl << std::flush;
}