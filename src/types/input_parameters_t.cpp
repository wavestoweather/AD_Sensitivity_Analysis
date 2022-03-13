#include "include/types/input_parameters_t.h"


input_parameters_t::input_parameters_t() {
    // Numerics
    t_end_prime = 100.0;  // Seconds
    dt_prime = 0.01;  // Seconds
    snapshot_index = 200;
    num_sub_steps = 1;
    // Filename for output
    OUTPUT_FILENAME = "data/output.nc";
    CHECKPOINT_FILENAME = "";
    ENS_CONFIG_FILENAME = "";
    FOLDER_NAME = "";
    // Filename for input
    INPUT_FILENAME = "";
    tracking_filename = "";

    // start_over = true;
    start_over_env = true;
    fixed_iteration = false;
    auto_type = 3;
    traj = 0;
    write_index = 100000;
    progress_index = 1000;
    ensemble = 0;
    start_time_idx = -1;
#ifdef MET3D
    start_time = std::nan("");
#endif
    current_time = std::nan("");
    id = 0;
    n_ensembles = 0;
    simulation_mode = trajectory_sensitvity_perturbance;
    delay_time_store = 0;
    track_initial_cond = false;
}


void input_parameters_t::to_json(
    nlohmann::json& j,
    const double &time) const {
    j["t_end_prime"] = this->t_end_prime;
    j["dt_prime"] = this->dt_prime;
    j["start_time_idx"] = this->start_time_idx;
#ifdef MET3D
    if (start_time_idx == -1)
        j["start_time"] = this->start_time;
#endif
    j["snapshot_index"] = this->snapshot_index;
    j["num_sub_steps"] = this->num_sub_steps;
    j["OUTPUT_FILENAME"] = this->OUTPUT_FILENAME;
    j["INPUT_FILENAME"] = this->INPUT_FILENAME;
    j["start_over_env"] = this->start_over_env;
    j["fixed_iteration"] = this->fixed_iteration;
    j["auto_type"] = this->auto_type;
    j["traj"] = this->traj;
    j["write_index"] = this->write_index;
    j["progress_index"] = this->progress_index;
    j["ensemble"] = this->ensemble;
    j["n_ensembles"] = this->n_ensembles;
    j["current_time"] = time;
    j["FOLDER_NAME"] = this->FOLDER_NAME;
    j["simulation_mode"] = this->simulation_mode;
    j["delay_time_store"] = this->delay_time_store;
}


void to_json(
    nlohmann::json& j,
    const input_parameters_t &input) {
    j["t_end_prime"] = input.t_end_prime;
    j["dt_prime"] = input.dt_prime;
    j["start_time_idx"] = input.start_time_idx;
#ifdef MET3D
    if (input.start_time_idx == -1)
        j["start_time"] = input.start_time;
#endif
    j["snapshot_index"] = input.snapshot_index;
    j["num_sub_steps"] = input.num_sub_steps;
    j["OUTPUT_FILENAME"] = input.OUTPUT_FILENAME;
    j["INPUT_FILENAME"] = input.INPUT_FILENAME;
    j["start_over_env"] = input.start_over_env;
    j["fixed_iteration"] = input.fixed_iteration;
    j["auto_type"] = input.auto_type;
    j["traj"] = input.traj;
    j["write_index"] = input.write_index;
    j["progress_index"] = input.progress_index;
    j["ensemble"] = input.ensemble;
    j["n_ensembles"] = input.n_ensembles;
    j["current_time"] = input.current_time;
    j["FOLDER_NAME"] = input.FOLDER_NAME;
    j["simulation_mode"] = input.simulation_mode;
    j["delay_time_store"] = input.delay_time_store;
}


void from_json(
    const nlohmann::json &j,
    input_parameters_t &input) {
    for (auto& el : j.items()) {
        auto first = el.key();
        if (first == "t_end_prime") {
            j.at(first).get_to(input.t_end_prime);
        } else if (first == "dt_prime") {
            j.at(first).get_to(input.dt_prime);
        } else if (first == "start_time_idx") {
            j.at(first).get_to(input.start_time_idx);
#ifdef MET3D
        } else if (first == "start_time") {
            j.at(first).get_to(input.start_time);
#endif
        } else if (first == "snapshot_index") {
            j.at(first).get_to(input.snapshot_index);
        } else if (first == "num_sub_steps") {
            j.at(first).get_to(input.num_sub_steps);
        } else if (first == "OUTPUT_FILENAME") {
            j.at(first).get_to(input.OUTPUT_FILENAME);
        } else if (first == "INPUT_FILENAME") {
            j.at(first).get_to(input.INPUT_FILENAME);
        } else if (first == "start_over_env") {
            j.at(first).get_to(input.start_over_env);
        } else if (first == "fixed_iteration") {
            j.at(first).get_to(input.fixed_iteration);
        } else if (first == "auto_type") {
            j.at(first).get_to(input.auto_type);
        } else if (first == "traj") {
            j.at(first).get_to(input.traj);
        } else if (first == "write_index") {
            j.at(first).get_to(input.write_index);
        } else if (first == "progress_index") {
            j.at(first).get_to(input.progress_index);
        } else if (first == "ensemble") {
            j.at(first).get_to(input.ensemble);
        } else if (first == "n_ensembles") {
            j.at(first).get_to(input.n_ensembles);
        } else if (first == "current_time") {
            j.at(first).get_to(input.current_time);
        } else if (first == "FOLDER_NAME") {
            j.at(first).get_to(input.FOLDER_NAME);
        } else if (first == "simulation_mode") {
            j.at(first).get_to(input.simulation_mode);
        } else if (first == "delay_time_store") {
            j.at(first).get_to(input.delay_time_store);
        } else {
            SUCCESS_OR_DIE(INPUT_NAME_CHECKPOINT_ERR);
        }
    }
}

// void input_parameters_t::put(
//     pt::ptree &ptree,
//     const double &time) const {
//     pt::ptree input_params;
//     input_params.put<double>("t_end_prime", t_end_prime);
//     input_params.put<double>("dt_prime", dt_prime);
//     input_params.put<int32_t>("start_time_idx", start_time_idx);
// #ifdef MET3D
//     if (start_time_idx == -1)
//         input_params.put<double>("start_time", start_time);
// #endif
//     input_params.put<int>("snapshot_index", snapshot_index);
//     input_params.put<uint64_t>("num_sub_steps", num_sub_steps);
//     input_params.put<std::string>("OUTPUT_FILENAME", OUTPUT_FILENAME);
//     input_params.put<std::string>("INPUT_FILENAME", INPUT_FILENAME);
//     input_params.put<bool>("start_over_env", start_over_env);
//     input_params.put<bool>("fixed_iteration", fixed_iteration);
//     input_params.put<uint32_t>("auto_type", auto_type);
//     input_params.put<uint32_t>("traj", traj);
//     input_params.put<uint32_t>("write_index", write_index);
//     input_params.put<uint64_t>("progress_index", progress_index);
//     input_params.put<uint32_t>("ensemble", ensemble);
//     input_params.put<uint32_t>("n_ensembles", n_ensembles);
//     input_params.put<double>("current_time", time);
//     input_params.put<std::string>("FOLDER_NAME", FOLDER_NAME);
//     input_params.put<int>("simulation_mode", simulation_mode);
//     input_params.put<double>("delay_time_store", delay_time_store);
//     ptree.add_child("input_params", input_params);
// }


// void input_parameters_t::put(
//     pt::ptree &ptree) const {
//     put(ptree, current_time);
// }


// int input_parameters_t::from_pt(
//     pt::ptree &ptree) {
//     int err = 0;
//     for (auto &it : ptree.get_child("input_params")) {
//         auto first = it.first;
//         if (first == "t_end_prime") {
//             t_end_prime = it.second.get_value<double>();
//         } else if (first == "dt_prime") {
//             dt_prime = it.second.get_value<double>();
//         } else if (first == "start_time_idx") {
//             start_time_idx = it.second.get_value<int32_t>();
// #ifdef MET3D
//         } else if (first == "start_time") {
//             start_time = it.second.get_value<double>();
// #endif
//         } else if (first == "snapshot_index") {
//             snapshot_index = it.second.get_value<int>();
//         } else if (first == "num_sub_steps") {
//             num_sub_steps = it.second.get_value<uint64_t>();
//         } else if (first == "OUTPUT_FILENAME") {
//             OUTPUT_FILENAME = it.second.data();
//         } else if (first == "INPUT_FILENAME") {
//             INPUT_FILENAME = it.second.data();
//         } else if (first == "start_over_env") {
//             start_over_env = (it.second.data() == "1"
//                 || it.second.data() == "true") ? true : false;
//         } else if (first == "fixed_iteration") {
//             fixed_iteration = (it.second.data() == "1"
//                 || it.second.data() == "true") ? true : false;
//         } else if (first == "auto_type") {
//             auto_type = it.second.get_value<uint32_t>();
//         } else if (first == "traj") {
//             traj = it.second.get_value<uint32_t>();
//         } else if (first == "write_index") {
//             write_index = it.second.get_value<uint32_t>();
//         } else if (first == "progress_index") {
//             progress_index = it.second.get_value<uint64_t>();
//         } else if (first == "ensemble") {
//             ensemble = it.second.get_value<uint32_t>();
//         } else if (first == "n_ensembles") {
//             n_ensembles = it.second.get_value<uint32_t>();
//         } else if (first == "current_time") {
//             current_time = it.second.get_value<double>();
//         } else if (first == "FOLDER_NAME") {
//             FOLDER_NAME = it.second.data();
//         } else if (first == "simulation_mode") {
//             simulation_mode = it.second.get_value<int>();
//         } else if (first == "delay_time_store") {
//             delay_time_store = it.second.get_value<double>();
//         } else {
//             err = INPUT_NAME_CHECKPOINT_ERR;
//         }
//     }
//     return err;
// }


void input_parameters_t::set_input_from_arguments(
    global_args_t &arg) {
    // Final time
    if (1 == arg.final_time_flag) {
        this->t_end_prime = std::strtod(arg.final_time_string, nullptr);
    }

    // Timestep
    if (1 == arg.timestep_flag) {
        this->dt_prime = std::strtod(arg.timestep_string, nullptr);
    }

    // Snapshot index
    if (1 == arg.snapshot_index_flag) {
        this->snapshot_index = std::stoi(arg.snapshot_index_string);
    }

    // Output
    if (1 == arg.output_flag) {
        this->OUTPUT_FILENAME = arg.output_string;
    }

    // Input
    if (1 == arg.input_flag) {
        this->INPUT_FILENAME = arg.input_file;
    }

    // Simulation mode
    if (1 == arg.simulation_mode_flag) {
        this->simulation_mode = std::stoi(arg.simulation_mode_string);
        switch (this->simulation_mode) {
            case trajectory_sensitvity_perturbance:
            case trajectory_sensitivity:
            case trajectory_perturbance:
            case limited_time_ensembles:
                break;
            case grid_sensitivity:
                std::cout << "Grid based sensitivity analysis is not supported yet\n";
                SUCCESS_OR_DIE(ARGUMENT_ERR);
                break;
            default:
                std::cout << "No such simulation mode. Aborting.\n";
                SUCCESS_OR_DIE(ARGUMENT_ERR);
                break;
        }
    }

    // Check if any tracking configuration is available
    if (1 == arg.tracking_file_flag) {
        this->tracking_filename = arg.tracking_file_string;
    }

    // Starting over environment variables (p, T, w)
    if (1 == arg.start_over_env_flag) {
        this->start_over_env = (strcmp(arg.start_over_env_string, "0"));
    }

    if (1 == arg.fixed_iteration_flag) {
        this->fixed_iteration = (strcmp(arg.fixed_iteration_string, "0"));
    }

    // Auto type
    if (1 == arg.auto_type_flag) {
        this->auto_type = std::stoi(arg.auto_type_string);
    }

    // Trajectory
    if (1 == arg.traj_flag) {
        this->traj = std::stoi(arg.traj_string);
    }

    // Write index
    if (1 == arg.write_flag) {
        this->write_index = std::stoi(arg.write_string);
    }

    // Progressbar index
    if (1 == arg.progress_index_flag) {
        this->progress_index = std::stoull(arg.progress_index_string);
    }

    // Start time index
    if (1 == arg.time_start_idx_flag) {
        this->start_time_idx = std::stol(arg.time_start_idx_string);
    }

#ifdef MET3D
    // Simulation start time
    if ((1 == arg.delay_start_flag) && (0 == arg.time_start_idx_flag)) {
        this->start_time = std::strtod(arg.delay_start_string, nullptr);
    }
#endif

    // Ensemble configuration file
    if (1 == arg.ens_config_flag) {
        this->ENS_CONFIG_FILENAME = arg.ens_config_string;
    }

    // Checkpoint file
    if (1 == arg.checkpoint_flag) {
        this->CHECKPOINT_FILENAME = arg.checkpoint_string;
    }

    // Folder name for new generated checkpoints
    if (1 == arg.folder_name_flag) {
        this->FOLDER_NAME = arg.folder_name_string;
    }

    // Maximum number of ensembles in the output
    if (1 == arg.n_ens_flag) {
        this->n_ensembles = std::stoi(arg.n_ens_string);
    }

    // Delay time in seconds before data is stored
    if (1 == arg.warm_up_flag) {
        this->delay_time_store = std::strtod(arg.warm_up_string, nullptr);
    }

    // Track sensitivity to initial conditions instead of model parameters
    if (1 == arg.track_ic_flag) {
        this->track_initial_cond =  (strcmp(arg.track_ic_string, "0"));
    }
    // If we want to run ensembles, we need a configuration file for that.
    if (1 != arg.ens_config_flag
        && (this->simulation_mode == trajectory_sensitvity_perturbance
         || this->simulation_mode == trajectory_perturbance
         || this->simulation_mode == limited_time_ensembles)) {
        std::cerr << "Simulation mode creates ensembles but no configuration file "
            << "is given.\n";
        SUCCESS_OR_DIE(ARGUMENT_ERR);
    }
}


void input_parameters_t::print_parameters() {
#ifdef SILENT_MODE
    return;
#endif
    std::cout << "\n"
        << "Technical input parameters:\n"
        << "---------------------------\n"
        << "Time to integrate: " << this->t_end_prime << " Seconds\n"
        << "Timestep: " << this->dt_prime << " Seconds\n"
        << (((this->CHECKPOINT_FILENAME == "") && (this->start_time_idx != -1))
        ? "Start time index: " + std::to_string(this->start_time_idx) + "\n"
        : "")
#ifdef MET3D
        << (((this->CHECKPOINT_FILENAME == "") && (this->start_time_idx == -1))
        ? "Start time (relative to ascend): " + std::to_string(this->start_time) + " Seconds\n"
        : "")
#endif
        << ((this->CHECKPOINT_FILENAME != "")
        ?   "Start time from checkpoint (relative to ascend): "
            + std::to_string(this->current_time + this->start_time) + "\n"
        :   "")
        << "Snapshot index: " << this->snapshot_index << "\n"
        << "Write index: " << this->write_index << "\n"
        << "Progressbar index: " << this->progress_index << "\n"
        << "Name of output file: " << this->OUTPUT_FILENAME << "\n"
        << "Name of input file: " << this->INPUT_FILENAME << "\n"
        << ((this->tracking_filename != "")
        ? "Tracking specified by file: " + this->tracking_filename + "\n"
        : "")
        << "Start over pressure, temperature and ascent at each timestep of a trajectory?: "
        << this->start_over_env << "\n"
        << "Fix temperature and pressure at each substep?: " << this->fixed_iteration << "\n"
        << "Auto type for rain evaporation (1, 2, 3): " << this->auto_type << "\n"
        << "Trajectory used: " << this->traj << "\n"
        << "Instance id: " << this->id << "\n"
        << ((this->ENS_CONFIG_FILENAME != "")
        ?   "Ensemble configuration file: " + this->ENS_CONFIG_FILENAME + "\n"
        :   "")
        << ((this->CHECKPOINT_FILENAME != "")
        ?   "Checkpoint file: " + this->CHECKPOINT_FILENAME + "\n"
        :   "")
        << ((this->FOLDER_NAME != "")
        ? "Folder name for newly generated checkpoints: " + this->FOLDER_NAME + "\n"
        : "")
        << (((this->simulation_mode == trajectory_sensitvity_perturbance)
            || (this->simulation_mode == trajectory_perturbance))
        ? "Maximum number of ensembles in the output file: " + std::to_string(this->n_ensembles) + "\n"
        : "")
        << ((this->track_initial_cond)
        ? "Tracking sensitivity to initial conditions instead of model parameters\n"
        : "")
        << "Simulation mode: " << this->simulation_mode << "\n"
        << "Warm-up time: " << this->delay_time_store << "s\n"
        << std::endl << std::flush;
}
