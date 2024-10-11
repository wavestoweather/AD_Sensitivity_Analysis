//
// Created by mahieron on 6/19/23.
// Calculate the histograms for leading order magnitudes
// relative to the time of ascent over all trajectories.
//
#include <mpi.h>
#include <netcdf.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <vector>

#include <chrono>

#include "include/types/netcdf_simulation_reader_t.h"
#include "include/misc/constants_sim_reader.h"
#include "include/misc/error.h"
#include "include/misc/pbar.h"

const double FILLVALUE = std::numeric_limits<double>::quiet_NaN();
const nc_type NC_FLOAT_T = NC_DOUBLE;


void set_indices(
        netcdf_simulation_reader_t &netcdf_reader,
        uint64_t &buffer_offset,
        uint64_t &start_time_idx,
        const double &inflow_time,
        const uint32_t &traj) {
    netcdf_reader.init_netcdf(traj);
    buffer_offset = (traj % netcdf_reader.n_traj_buffer) * netcdf_reader.read_time_buffer;
    for (uint32_t i = 0; i < 10; i++) {
        if (std::isnan(netcdf_reader.buffer[Par_idx::time_after_ascent][i + buffer_offset])) continue;
        if (inflow_time <= netcdf_reader.buffer[Par_idx::time_after_ascent][i + buffer_offset]) {
            start_time_idx = 0;
        } else {
            auto rel_time = netcdf_reader.buffer[Par_idx::time_after_ascent][i + buffer_offset];
            start_time_idx = i + (inflow_time - rel_time) / DELTA_TIMESTEP;
            break;
        }
    }
    buffer_offset += start_time_idx;
}


void get_loom(
        netcdf_simulation_reader_t &netcdf_reader,
        std::vector<double> &loom,
        std::vector<double> &loom_tmp,
        std::array<double, 2> &time_limits,
        std::vector<uint64_t> &counter_arr,
        std::vector<uint64_t> &n_count_avg,
        uint64_t &n_trajectories,
        const double &inflow_time,
        const double &outflow_time,
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const bool &max_loom,
        const bool &per_traj,
        const bool &ignore_phases,
        const uint32_t &n_phases,
        const uint32_t &count_avg) {
    ProgressBar pbar = ProgressBar(netcdf_reader.n_trajectories, 50, "Trajectory", std::cout);
    n_trajectories += netcdf_reader.n_trajectories;
    for (uint32_t traj = 0; traj < netcdf_reader.n_trajectories; traj++) {
        uint64_t buffer_offset;
        uint64_t start_time_idx;
        set_indices(netcdf_reader, buffer_offset, start_time_idx, inflow_time, traj);

        // For every output parameter
        for (uint32_t j = 0; j < n_out_params; j++) {
            bool ascent_started = false;
            uint32_t outflow_counter = 0;
            // For every time step
            for (uint32_t t = 0; t < netcdf_reader.n_timesteps_in - start_time_idx; t++) {
                auto current_phase = netcdf_reader.buffer[Par_idx::phase][t + buffer_offset];
                // Check if still in valid time
                // Finished if outflow limit is reached
                if (ascent_started) {
                    if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset])
                        outflow_counter += 1;
                } else {
                    if (netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset])
                        ascent_started = true;
                }
                if (outflow_counter * DELTA_TIMESTEP > outflow_time) {
                    break;
                }
                if (current_phase == 3) continue;
                if (ignore_phases) current_phase = 0;
                // For every model parameter
                for (uint32_t i = 0; i < n_params; i++) {
                    auto val = netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset];
                    if (std::isnan(val)) continue;
                    // We calculate loom per trajectories later on the fly
                    if (!per_traj) {
                        if (max_loom) {
                            // Get max value
                            loom[j * n_phases + current_phase] =
                                    (loom[j * n_phases + current_phase] < fabs(val)) ? fabs(val)
                                                                                 : loom[j * n_phases + current_phase];
                        } else {
                            if (i == 0) {
                                counter_arr[j * n_phases + current_phase]++;
                                if (counter_arr[j * n_phases + current_phase] % count_avg == 0) {
                                    n_count_avg[j * n_phases + current_phase]++;
                                }
                            }
                            loom_tmp[(i + 1) * n_phases * n_out_params + j * n_phases + current_phase] += fabs(val);
                            if (counter_arr[j * n_phases + current_phase] % count_avg == 0) {
                                loom[(i + 1) * n_phases * n_out_params + j * n_phases + current_phase] +=
                                        loom_tmp[(i + 1) * n_phases * n_out_params + j * n_phases + current_phase] /
                                        count_avg;
                                loom_tmp[(i + 1) * n_phases * n_out_params + j * n_phases + current_phase] = 0;
                            }
                        }
                    }
                    time_limits[0] = (
                            time_limits[0] > netcdf_reader.buffer[Par_idx::time_after_ascent][t + buffer_offset])
                                 ? netcdf_reader.buffer[Par_idx::time_after_ascent][t + buffer_offset]
                                 : time_limits[0];
                    time_limits[1] = (
                            time_limits[1] < netcdf_reader.buffer[Par_idx::time_after_ascent][t + buffer_offset])
                                 ? netcdf_reader.buffer[Par_idx::time_after_ascent][t + buffer_offset]
                                 : time_limits[1];
                }
            }
        }
        pbar.progress();
    }
}

void get_means(
        std::vector<double> &loom,
        std::vector<double> &loom_tmp,
        std::vector<uint64_t> &counter_arr,
        std::vector<uint64_t> &n_count_avg,
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const uint32_t &n_phases,
        const uint32_t &count_avg) {
    for (uint32_t i = 0; i < n_params; i++) {
        for (uint32_t j = 0; j < n_out_params; j++) {
            for (uint32_t phase_idx = 0; phase_idx < n_phases; phase_idx++) {
                auto val = loom[(i + 1) * n_phases * n_out_params + j * n_phases + phase_idx];
                if (counter_arr[j * n_phases + phase_idx] == 0) {
                    continue;
                }
                if (counter_arr[j * n_phases + phase_idx]%count_avg != 0) {
                    auto left_over = counter_arr[j * n_phases + phase_idx]%count_avg;
                    auto n_count_this = n_count_avg[j * n_phases + phase_idx];
                    val = loom_tmp[(i+1) * n_phases * n_out_params + j * n_phases + phase_idx]/count_avg
                            / (static_cast<double>(n_count_this)
                                + (static_cast<double>(left_over)/static_cast<double>(count_avg)))
                           + val / (n_count_this + static_cast<double>(left_over)/static_cast<double>(count_avg));
                } else {
                    val /= n_count_avg[j * n_phases + phase_idx];
                }
                loom[(i + 1) * n_phases * n_out_params + j * n_phases + phase_idx] = val;
                if (val > loom[j * n_phases + phase_idx]) {
                    loom[j * n_phases + phase_idx] = val;
                }
            }
        }
    }
}


void get_loom_per_traj(
        const netcdf_simulation_reader_t &netcdf_reader,
        std::vector<double> &loom,
        const double &outflow_time,
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const bool &max_loom,
        const bool &ignore_phases,
        const uint64_t &start_time_idx,
        const uint64_t &buffer_offset,
        const uint32_t &n_phases,
        std::vector<uint64_t> &tmp_counter) {
    std::fill(tmp_counter.begin(), tmp_counter.end(), 0);
    std::fill(loom.begin(), loom.end(), 0);
    for (uint32_t i = 0; i < n_params; i++) {
        for (uint32_t j = 0; j < n_out_params; j++) {
            bool ascent_started = false;
            uint32_t outflow_counter = 0;
            for (uint32_t t = 0; t < netcdf_reader.n_timesteps_in - start_time_idx; t++) {
                auto current_phase = netcdf_reader.buffer[Par_idx::phase][t + buffer_offset];
                // Check if still in valid time
                // Finished if outflow limit is reached
                if (ascent_started) {
                    if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset])
                        outflow_counter += 1;
                } else {
                    if (netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset])
                        ascent_started = true;
                }
                if (outflow_counter * DELTA_TIMESTEP > outflow_time) {
                    break;
                }
                if (current_phase == 3) continue;
                if (ignore_phases) current_phase = 0;
                auto val = fabs(netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset]);
                if (std::isnan(val)) continue;
                if (max_loom && val > loom[j * n_phases + current_phase]) {
                    loom[j * n_phases + current_phase] = val;
                } else if (!max_loom) {
                    loom[(i + 1) * n_phases * n_out_params + j * n_phases + current_phase] += val;
                    if ((i == 0) && (j == 0)) tmp_counter[current_phase]++;
                }
            }
        }
    }
    if (!max_loom) {
        for (uint32_t i = 0; i < n_params; i++) {
            for (uint32_t j = 0; j < n_out_params; j++) {
                for (uint32_t phase_idx = 0; phase_idx < n_phases; phase_idx++) {
                    if (ignore_phases && phase_idx > 0)
                        break;
                    auto val = loom[(i + 1) * n_phases * n_out_params + j * n_phases + phase_idx];
                    val /= tmp_counter[phase_idx];
                    loom[j * n_phases + phase_idx] = (
                                                             loom[j * n_phases + phase_idx] < val)
                                                     ? val : loom[j * n_phases + phase_idx];
                }
            }
        }
    }
}


void count_loom_in_traj(
        const netcdf_simulation_reader_t &netcdf_reader,
        const std::vector<double> &loom,
        const std::array<double, 2> &time_limits,
        std::vector<uint64_t> &loom_count,
        std::vector<uint64_t> &count,
        std::vector<uint64_t> &traj_counter,
        const double &outflow_time,
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const bool &ignore_phases,
        const uint32_t &corr_time,
        const uint64_t &start_time_idx,
        const uint64_t &buffer_offset,
        const uint32_t &n_phases,
        const uint32_t &traj,
        const uint64_t &n_times,
        uint32_t &traj_old) {
    // For every model parameter
    uint32_t i_old = -1;
    for (uint32_t i = 0; i < n_params; i++) {
        // For every output parameter
        uint32_t j_old = -1;
        for (uint32_t j = 0; j < n_out_params; j++) {
            bool ascent_started = false;
            bool ignore_parameter_target = false;
            // Use a temporary buffer for the counts to filter
            // any sensitivities that are not at least
            // corr_time [%] long during the ascent.
            if (corr_time > 0) {
                double loom_counter_tmp = 0;
                double ascent_time = 0;
                for (uint32_t t = 0; t < netcdf_reader.n_timesteps_in - start_time_idx; t++) {
                    auto current_phase = netcdf_reader.buffer[Par_idx::phase][t + buffer_offset];
                    if (current_phase == 3) continue;
                    if (ignore_phases) current_phase = 0;
                    if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset] && !ascent_started) continue;
                    if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset]) break;
                    ascent_started = true;
                    ascent_time += 1;
                    auto val = netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset];
                    if (std::isnan(val)) continue;
                    if ((fabs(val) >= loom[j * n_phases + current_phase] / 10) &&
                        (loom[j * n_phases + current_phase] > 0))
                        loom_counter_tmp += 1;
                }
                if (loom_counter_tmp / ascent_time * 100 < corr_time) ignore_parameter_target = true;
            }
            ascent_started = false;
            uint32_t outflow_counter = 0;
            std::array<bool, 3> phase_counted = {false, false, false};
            // For every time step
            for (uint32_t t = 0; t < netcdf_reader.n_timesteps_in - start_time_idx; t++) {
                int current_phase = netcdf_reader.buffer[Par_idx::phase][t + buffer_offset];
                // Check if still in valid time
                // Use only asc600 if correlation time is used
                if (corr_time > 0) {
                    if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset] && !ascent_started) continue;
                    if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset]) break;
                    ascent_started = true;
                } else {
                    // Finished if outflow limit is reached
                    if (ascent_started) {
                        if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset])
                            outflow_counter += 1;
                    } else {
                        if (netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset])
                            ascent_started = true;
                    }
                    if (outflow_counter * DELTA_TIMESTEP > outflow_time) {
                        break;
                    }
                }

                if (current_phase == 3) continue;
                if (ignore_phases) current_phase = 0;
                auto val = netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset];
                if (std::isnan(val)) continue;
                // Get index for relative time
                auto current_time = netcdf_reader.buffer[Par_idx::time_after_ascent][t + buffer_offset];
                uint32_t time_idx = (current_time - time_limits[0])/DELTA_TIMESTEP;
                if ((i == j) && (j == 0)) count[time_idx + current_phase * n_times]++;
                uint64_t loom_idx =
                        t + current_phase * n_times
                        + j * n_phases * n_times
                        + i * n_out_params * n_phases * n_times;
                if ((fabs(val) >= loom[j * n_phases + current_phase]/10) && (loom[j * n_phases + current_phase] > 0)
                    && !ignore_parameter_target) {
                    loom_count[loom_idx]++;
                    if (traj != traj_old || j != j_old || i != i_old || !phase_counted[current_phase]) {
                        traj_counter[j + current_phase*n_out_params + i * n_out_params * n_phases ]++;
                        traj_old = traj;
                        j_old = j;
                        i_old = i;
                        phase_counted[current_phase] = true;
                    }
                }
            }
        }
    }
}


void count_loom_per_time(
        const netcdf_simulation_reader_t &netcdf_reader,
        const std::array<double, 2> &time_limits,
        std::vector<uint64_t> &loom_count,
        std::vector<uint64_t> &count,
        std::vector<uint64_t> &traj_counter,
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const uint32_t &corr_time,
        const uint64_t &start_time_idx,
        const uint64_t &buffer_offset,
        const uint32_t &n_phases,
        const uint32_t &traj,
        const uint64_t &n_times,
        uint32_t &traj_old) {
    std::vector<bool> ignore_params(n_params * n_out_params);
    std::fill(ignore_params.begin(), ignore_params.end(), 0);
    if (corr_time > 0) {
        bool ascent_started = false;
        std::vector<double> loom_counter_tmp(n_params * n_out_params);
        std::fill(loom_counter_tmp.begin(), loom_counter_tmp.end(), 0);
        double ascent_time = 0;
        for (uint32_t t = 0; t < netcdf_reader.n_timesteps_in - start_time_idx; t++) {
            if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset] && !ascent_started) continue;
            if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset]) break;
            ascent_started = true;
            ascent_time += 1;
            auto current_phase = netcdf_reader.buffer[Par_idx::phase][t + buffer_offset];
            if (current_phase == 3) continue;
            current_phase = 0;
            for (uint32_t j = 0; j < n_out_params; j++) {
                double max_val = 0;
                // Get the top magnitude first
                for (uint32_t i = 0; i < n_params; i++) {
                    auto val = netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset];
                    if (std::isnan(val)) continue;
                    max_val = (fabs(val) > max_val) ? fabs(val) : max_val;
                }
                for (uint32_t i = 0; i < n_params; i++) {
                    auto val = netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset];
                    if (std::isnan(val)) continue;
                    if ((fabs(val) >= max_val / 10) &&
                        (max_val > 0))
                        loom_counter_tmp[i * 3 + j] += 1;
                }
            }
        }
        for (uint32_t j = 0; j < n_out_params; j++) {
            for (uint32_t i = 0; i < n_params; i++) {
                if (loom_counter_tmp[i * 3 + j] / ascent_time * 100 < corr_time) ignore_params[i * 3 + j] = true;
            }
        }
    }
    // For every time step
    bool ascent_started = false;
    uint32_t outflow_counter = 0;
    for (uint32_t t = 0; t < netcdf_reader.n_timesteps_in - start_time_idx; t++) {
        // This approach does not consider inflow or outflow.
        if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset] && !ascent_started) continue;
        if (!netcdf_reader.buffer[Par_idx::asc600][t + buffer_offset]) break;
        ascent_started = true;
        int current_phase = netcdf_reader.buffer[Par_idx::phase][t + buffer_offset];
        // To stay consistent with other methods, ignore steps with no hydrometeor in it.
        if (current_phase == 3) continue;
        // This approach ignores phases
        current_phase = 0;
        auto current_time = netcdf_reader.buffer[Par_idx::time_after_ascent][t + buffer_offset];
        // For every output parameter
        uint32_t j_old = -1;
        for (uint32_t j = 0; j < n_out_params; j++) {
            double max_val = 0;
            // Get the top magnitude first
            for (uint32_t i = 0; i < n_params; i++) {
                auto val = netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset];
                if (std::isnan(val)) continue;
                max_val = (fabs(val) > max_val) ? fabs(val) : max_val;
            }
            // For every model parameter
            uint32_t i_old = -1;
            for (uint32_t i = 0; i < n_params; i++) {
                if (ignore_params[i * 3 + j]) continue;
                auto val = netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset];
                if (std::isnan(val)) continue;
                uint32_t time_idx = (current_time - time_limits[0])/DELTA_TIMESTEP;
                if ((i == j) && (j == 0)) count[time_idx + current_phase * n_times]++;
                uint64_t loom_idx =
                        t + current_phase * n_times
                        + j * n_phases * n_times
                        + i * n_out_params * n_phases * n_times;
                if ((fabs(val) >= max_val/10) && max_val > 0) {
                    loom_count[loom_idx]++;
                    if (traj != traj_old || j != j_old || i != i_old) {
                        traj_counter[j + current_phase*n_out_params + i * n_out_params * n_phases]++;
                        traj_old = traj;
                        j_old = j;
                        i_old = i;
                    }
                }
            }
        }
    }
}


void get_loom_histogram(
        netcdf_simulation_reader_t &netcdf_reader,
        std::vector<double> &loom,
        std::array<double, 2> &time_limits,
        std::vector<uint64_t> &loom_count,
        std::vector<uint64_t> &count,
        std::vector<uint64_t> &traj_counter,
        const double &inflow_time,
        const double &outflow_time,
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const uint32_t &n_phases,
        const bool &max_loom,
        const bool &per_traj,
        const bool &per_time,
        const bool &ignore_phases,
        const uint32_t &corr_time) {
    const uint64_t n_times = (time_limits[1] - time_limits[0])/DELTA_TIMESTEP + 1;
    ProgressBar pbar = ProgressBar(netcdf_reader.n_trajectories, 50, "Trajectory", std::cout);
    std::vector<uint64_t> tmp_counter(n_phases);
    uint32_t traj_old = -1;
    for (uint32_t traj = 0; traj < netcdf_reader.n_trajectories; traj++) {
        uint64_t buffer_offset;
        uint64_t start_time_idx;
        set_indices(netcdf_reader, buffer_offset, start_time_idx, inflow_time, traj);

        if (per_traj && !per_time) {
            get_loom_per_traj(
                    netcdf_reader,
                    loom,
                    outflow_time,
                    n_params,
                    n_out_params,
                    max_loom,
                    ignore_phases,
                    start_time_idx,
                    buffer_offset,
                    n_phases,
                    tmp_counter);
            count_loom_in_traj(
                    netcdf_reader,
                    loom,
                    time_limits,
                    loom_count,
                    count,
                    traj_counter,
                    outflow_time,
                    n_params,
                    n_out_params,
                    ignore_phases,
                    corr_time,
                    start_time_idx,
                    buffer_offset,
                    n_phases,
                    traj,
                    n_times,
                    traj_old);
        } else {
            count_loom_per_time(
                    netcdf_reader,
                    time_limits,
                    loom_count,
                    count,
                    traj_counter,
                    n_params,
                    n_out_params,
                    corr_time,
                    start_time_idx,
                    buffer_offset,
                    n_phases,
                    traj,
                    n_times,
                    traj_old);
        }
        pbar.progress();
    }
}


void parse(
        int argc,
        char** argv,
        uint32_t &buffer_size,
        double &inflow_time,
        double &outflow_time,
        bool &max_loom,
        bool &per_traj,
        bool &per_time,
        bool &ignore_phases,
        uint32_t &corr_time,
        std::string &store_path,
        std::string &load_path,
        uint32_t &n_traj_buffer) {
    if (argc > 64) {
        throw std::runtime_error("You provided more than 64 input parameters.");
    }
    // some default values
    buffer_size = 8600;
    inflow_time = 7200;
    outflow_time = 7200;
    max_loom = false;
    per_traj = false;
    ignore_phases = false;
    corr_time = 0;
    n_traj_buffer = 1000;
    per_time = false;
    const std::vector<std::string> args(argv + 1, argv + argc);
    for (auto i = 0; i < args.size(); i += 2) {
        auto arg = args[i];
        if (arg == "-b" || arg == "--buffer_size") {
            buffer_size = std::stoi(args[i+1]);
        } else if (arg == "--inflow_time") {
            inflow_time = std::stod(args[i+1], nullptr)  * (-1);
        } else if (arg == "--outflow_time") {
            outflow_time = std::stod(args[i+1], nullptr);
        } else if (arg == "-i" || arg == "--input_path") {
            load_path = args[i+1];
        } else if (arg == "-o" || arg == "--output_path") {
            store_path = args[i+1];
        } else if (arg == "-m" || arg == "--max_loom") {
            max_loom = true;
            i--;
        } else if (arg == "--per_traj") {
            per_traj = true;
            i--;
        } else if (arg == "--per_time") {
            per_time = true;
            per_traj = true;
            ignore_phases = true;
            i--;
        } else if (arg == "--traj_buffer") {
            n_traj_buffer = std::stoi(args[i + 1]);
        } else if (arg == "-c" || arg == "--correlation_time") {
            corr_time = std::stoi(args[i + 1], nullptr, 0);
        } else if (arg == "--ignore_phases") {
            ignore_phases = true;
            i--;
        } else {
            std::cout << "Argument '" << arg << "' not found. It is ignored.\n"
            << "Possible options are:\n"
            << "-b, --buffer_size: Adjust the number of time steps to store in buffer. Should be enough "
            << "to hold all time steps of a trajectory.\n"
            << "--inflow_time: Time in seconds to consider before the ascent starts.\n"
            << "--outflow_time: Time in seconds to consider after the ascent ends.\n"
            << "-i, --input_path: Path to folder with NetCDF-files or to a single NetCDF-file from a "
            << "sensitivity simulation.\n"
            << "-o, --output_path: Path and name (*.nc) to store the final data.\n"
            << "-m, --max_loom: If set, use the maximum leading order of magnitude for each phase as threshold. "
            << "Otherwise, use the mean.\n"
            << "--per_traj: If set, calculate the leading order of magnitude per trajectory. "
            << "Otherwise, calculate the leading order of magnitude over all trajectories.\n"
            << "--per_time: If set, calculate the leading order of magnitude per time step. "
            << "Automatically sets per_traj and ignore_phases. Ignores inflow_time and outflow_time (set to zero).\n"
            << "--traj_buffer: The number of trajectories to store in RAM at once."
            << "-c, --correlation_time: Count only those leading order of magnitude if "
            << "the total amount of datapoints in a trajectory exceeds -c % of "
            << "all datapoints during the ascent. Ignored if --per_traj is not set.\n"
            << "--ignore_phases: If set, ignore cloud phases for leading order calculation.\n";
            i--;
        }
    }
}

int create_dims(
        const std::string& store_path,
        const uint32_t& n_phases,
        const uint32_t& n_times,
        std::vector<int> &dimid,
        const bool &ignore_phases) {
    int ncid;
    SUCCESS_OR_DIE(nc_create(
            store_path.c_str(),
            NC_NETCDF4,
            &ncid));
    // Create dimensions
    SUCCESS_OR_DIE(
        nc_def_dim(
            ncid,                // ncid
            "Output_Parameter",  // name
            3,                   // length
            &dimid[Loom_dim_idx::outp_loomdim_idx]));   // idp
    SUCCESS_OR_DIE(
        nc_def_dim(
            ncid,
            "time",
            n_times,
            &dimid[Loom_dim_idx::time_loomdim_idx]));
    if (!ignore_phases) {
        SUCCESS_OR_DIE(
            nc_def_dim(
                ncid,
                "phase",
                n_phases,
                &dimid[Loom_dim_idx::phase_loomdim_idx]));;
    } else {
        dimid[Loom_dim_idx::phase_loomdim_idx] = dimid[Loom_dim_idx::time_loomdim_idx];
    }

    return ncid;
}

void define_vars(
        const int &ncid,
        std::vector<int>& varid,
        std::vector<int>& dimid,
        std::vector<int>& loom_order_varid,
        const uint32_t &n_params,
        const bool &max_loom,
        const bool &ignore_phases,
        const bool get_loom_floats = false) {
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,                   // ncid
            "Output_Parameter",     // name
            NC_STRING,              // type
            1,                      // ndims (0=scalar, 1=vector, 2=matrix, ...)
            &dimid[Loom_dim_idx::outp_loomdim_idx],     // dimid
            &varid[Loom_dim_idx::outp_loomdim_idx]));   // varid
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,
            "time",
            NC_FLOAT_T,
            1,
            &dimid[Loom_dim_idx::time_loomdim_idx],
            &varid[Loom_dim_idx::time_loomdim_idx]));
    if (!ignore_phases) {
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                "phase",
                NC_STRING,
                1,
                &dimid[Loom_dim_idx::phase_loomdim_idx],
                &varid[Loom_dim_idx::phase_loomdim_idx]));
    }

    int ndims = (ignore_phases) ? 1 : 2;
    // Now all the other things
    // Loom count and trajectory with loom count
    for (uint32_t i = 0; i < n_params; i++) {
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                loom_order_sens[i].c_str(),
                NC_UINT64,
                ndims+1,
                &dimid[Loom_dim_idx::outp_loomdim_idx],
                &loom_order_varid[i]));
        // Trajectory counts where loom was reached
        std::string name = loom_order_sens[i] + " trajectory counts";
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                name.c_str(),
                NC_UINT64,
                ndims,
                &dimid[Loom_dim_idx::outp_loomdim_idx],
                &loom_order_varid[n_params + n_params * get_loom_floats + get_loom_floats + 1 + i]));
    }
    // Loom average/max
    if (get_loom_floats) {
        for (uint32_t i = 0; i < n_params; i++) {
            std::string name = loom_order_sens[i] + " sensitivity";
            SUCCESS_OR_DIE(
                nc_def_var(
                    ncid,
                    name.c_str(),
                    NC_DOUBLE,
                    ndims,
                    &dimid[Loom_dim_idx::outp_loomdim_idx],
                    &loom_order_varid[i + n_params]));
        }
        // Actual thresholds
        std::string name = (max_loom) ? "Maximum" : "Average";
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                name.c_str(),
                NC_DOUBLE,
                ndims,
                &dimid[Loom_dim_idx::outp_loomdim_idx],
                &loom_order_varid[n_params + n_params + 1]));
    }
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,
            "counts",
            NC_UINT64,
            ndims,
            &dimid[Loom_dim_idx::phase_loomdim_idx],
            &loom_order_varid[n_params + n_params * get_loom_floats]));
}

void set_attributes(
        const int &ncid,
        const bool &max_loom,
        const std::vector<int>& loom_order_varid,
        const uint32_t &n_params,
        const bool &get_loom_floats,
        const bool &per_traj,
        const bool &per_time,
        const bool &ignore_phases,
        const double &inflow_time,
        const double &outflow_time,
        const uint32_t &corr_time,
        const uint64_t &n_trajectories) {
    std::string name = "Time before ascent starts (s)";
    SUCCESS_OR_DIE(
        nc_put_att_double(
            ncid,
            NC_GLOBAL,
            name.c_str(),
            NC_DOUBLE,
            1,
            &inflow_time));
    name = "Time after ascent ended (s)";
    SUCCESS_OR_DIE(
        nc_put_att_double(
            ncid,
            NC_GLOBAL,
            name.c_str(),
            NC_DOUBLE,
            1,
            &outflow_time));
    name = "Selection based on threshold per trajectory";
    const int per_traj_int = static_cast<int>(per_traj);
    SUCCESS_OR_DIE(
        nc_put_att_int(
            ncid,
            NC_GLOBAL,
            name.c_str(),
            NC_INT,
            1,
            &per_traj_int));
    name = "Selection based on threshold per time step";
    const int per_time_int = static_cast<int>(per_time);
    SUCCESS_OR_DIE(
        nc_put_att_int(
            ncid,
            NC_GLOBAL,
            name.c_str(),
            NC_INT,
            1,
            &per_time_int));
    name = "Threshold based on";
    std::string descr = (max_loom) ? "Maximum sensitivity" : "Average sensitivity";
    SUCCESS_OR_DIE(
        nc_put_att_text(
            ncid,
            NC_GLOBAL,
            name.c_str(),
            strlen(descr.c_str()),
            descr.c_str()));
    if (corr_time > 0) {
        name = "Minimum fraction of ascent time for selection in percentage";
        const int corr_time_int = static_cast<int>(corr_time);
        SUCCESS_OR_DIE(
            nc_put_att_int(
                ncid,
                NC_GLOBAL,
                name.c_str(),
                NC_INT,
                1,
                &corr_time_int));
    }
    name = "Total number of trajectories.";
    const int n_trajs = static_cast<int>(n_trajectories);
    SUCCESS_OR_DIE(
        nc_put_att_int(
            ncid,
            NC_GLOBAL,
            name.c_str(),
            NC_INT,
            1,
            &n_trajs));
    if (ignore_phases && per_traj) {
        name = "Ignore phases";
        descr = "Phases dimension has no meaning.";
        SUCCESS_OR_DIE(
            nc_put_att_text(
                ncid,
                NC_GLOBAL,
                name.c_str(),
                strlen(descr.c_str()),
                descr.c_str()));
    }
    // Loom average/max
    descr = (max_loom) ? "Maximum sensitivity over all trajectories and time steps"
                                   : "Average sensitivity over all trajectories and time steps";
    std::string unit = "kg/kg/30s";
    if (get_loom_floats) {
        for (uint32_t i = 0; i < n_params; i++) {
            SUCCESS_OR_DIE(
                nc_put_att_text(
                    ncid,
                    loom_order_varid[i + n_params],
                    "long_name",
                    strlen(descr.c_str()),
                    descr.c_str()));
            SUCCESS_OR_DIE(
                nc_put_att_text(
                    ncid,
                    loom_order_varid[i + n_params],
                    "unit",
                    strlen(unit.c_str()),
                    unit.c_str()));
        }
        descr = (max_loom) ? "Maximum sensitivity over all trajectories, time steps, and parameters"
                           : "Maximum over all parameters of their average sensitivity " \
                             "over all trajectories and time steps";
        SUCCESS_OR_DIE(
            nc_put_att_text(
                ncid,
                loom_order_varid[2*n_params + 1],
                "long_name",
                strlen(descr.c_str()),
                descr.c_str()));
        SUCCESS_OR_DIE(
            nc_put_att_text(
                ncid,
                loom_order_varid[2*n_params + 1],
                "unit",
                strlen(unit.c_str()),
                unit.c_str()));
    }
    // Loom count
    descr = "Amount of times the parameter is within a magnitude of the threshold "
            "over all trajectories at the given time step";
    unit = "#";
    for (uint32_t i = 0; i < n_params; i++) {
        SUCCESS_OR_DIE(
            nc_put_att_text(
                ncid,
                loom_order_varid[i],
                "long_name",
                strlen(descr.c_str()),
                descr.c_str()));
        SUCCESS_OR_DIE(
            nc_put_att_text(
                ncid,
                loom_order_varid[i],
                "unit",
                strlen(unit.c_str()),
                unit.c_str()));
    }
    descr = "Amount of trajectories where threshold is met.";
    for (uint32_t i = 0; i < n_params; i++) {
        SUCCESS_OR_DIE(
            nc_put_att_text(
                ncid,
                loom_order_varid[n_params + n_params * get_loom_floats + get_loom_floats + 1 + i],
                "long_name",
                strlen(descr.c_str()),
                descr.c_str()));
    }
    // Total count
    descr = "Total number of datapoints";
    SUCCESS_OR_DIE(
        nc_put_att_text(
            ncid,
            loom_order_varid[n_params + n_params * get_loom_floats],
            "long_name",
            strlen(descr.c_str()),
            descr.c_str()));
    SUCCESS_OR_DIE(
        nc_put_att_text(
            ncid,
            loom_order_varid[n_params + n_params * get_loom_floats],
            "unit",
            strlen(unit.c_str()),
            unit.c_str()));
}

void write_dim_values(
        const int &ncid,
        const std::vector<int> &varid,
        const std::array<double, 2> &time_limits,
        const uint32_t &n_times,
        const bool &ignore_phases) {
    std::vector<size_t> startp, countp;
    startp.push_back(0);
    std::vector<std::string> out_params = {"QV", "latent_heat", "latent_cool"};

    std::vector<const char*> cstrings_outp;
    cstrings_outp.reserve(3);
    cstrings_outp.push_back("QV");
    cstrings_outp.push_back("latent_heat");
    cstrings_outp.push_back("latent_cool");
    for (auto s : cstrings_outp) {
        SUCCESS_OR_DIE(
            nc_put_var1_string(
                ncid,
                varid[Loom_dim_idx::outp_loomdim_idx],
                startp.data(),
                &s));
        startp[0]++;
    }
    startp[0] = 0;

    std::vector<double> times(n_times);
    uint32_t idx = 0;
    for (double t = time_limits[0]; t <= time_limits[1]; t += DELTA_TIMESTEP) {
        times[idx] = t;
        idx++;
    }
    countp.push_back(n_times);
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Loom_dim_idx::time_loomdim_idx],
            startp.data(),
            countp.data(),
            times.data()));
    if (!ignore_phases) {
        cstrings_outp[0] = "warm phase";
        cstrings_outp[1] = "mixed phase";
        cstrings_outp[2] = "ice phase";
        for (auto s : cstrings_outp) {
            SUCCESS_OR_DIE(
                    nc_put_var1_string(
                        ncid,
                    varid[Loom_dim_idx::phase_loomdim_idx],
                    startp.data(),
                    &s));
            startp[0]++;
        }
    }
}

void set_compression(
        const int &ncid,
        const std::vector<int> &varid,
        const bool &ignore_phases) {
    for (int i=0; i < varid.size(); i++) {
        if (!ignore_phases || i != Loom_dim_idx::phase_loomdim_idx) {
            SUCCESS_OR_DIE(
                nc_def_var_deflate(
                    ncid,
                    varid[i],
                    1,  // shuffle
                    1,  // deflate
                    COMPRESSION_LEVEL));  // compression
        }
    }
}


int main(int argc, char** argv) {
    int rank;
    int n_processes;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string load_path_str;
    std::string store_path;
    uint32_t buffer_size, n_traj_buffer, corr_time;
    double inflow_time, outflow_time;
    bool max_loom, per_traj, ignore_phases, per_time;
    parse(
        argc,
        argv,
        buffer_size,
        inflow_time,
        outflow_time,
        max_loom,
        per_traj,
        per_time,
        ignore_phases,
        corr_time,
        store_path,
        load_path_str,
        n_traj_buffer);
    const std::filesystem::path load_path = {load_path_str};
    const uint32_t n_params = loom_order_sens.size();
    netcdf_simulation_reader_t netcdf_reader(buffer_size, false, n_traj_buffer);
    const uint32_t n_out_params = 3;
    const uint32_t n_phases = (ignore_phases) ? 1 : 3;
    uint32_t loom_size = n_out_params * n_phases;
    const uint32_t counter_size = n_phases * n_out_params;
    std::vector<uint64_t> counter_arr, n_count_avg, traj_counter;
    std::vector<double> loom_tmp;
    const uint32_t count_avg = 10000;
    uint64_t n_trajectories = 0;
    bool get_loom_floats = false;
    // We use the first n_out_params * n_phases to store the maximum average sensitivity.
    if (!max_loom) {
        loom_size *= (n_params+1);
        counter_arr.resize(counter_size);
        loom_tmp.resize(loom_size);
        n_count_avg.resize(counter_size);
        if (!per_traj) get_loom_floats = true;
    }
    std::vector<double> loom(loom_size);
    std::array<double, 2> time_limits = {
            std::numeric_limits<double>::infinity(),
            -1*std::numeric_limits<double>::infinity()};
    double mem_usage = sizeof(double) *
                       (loom.size() + time_limits.size() + loom_tmp.size())
                       + sizeof(uint64_t) * (counter_arr.size() + n_count_avg.size());
    if (max_loom) {
        std::cout << "Using leading order of magnitude based on the overall maximum sensitivity.\n";
    } else {
        std::cout << "Using leading order of magnitude based on the maximum average sensitivity.\n";
    }
    if (per_traj) {
        std::cout << "Using leading order of magnitude independently for each trajectory.\n";
    } else {
        std::cout << "Using leading order of magnitude over all trajectories.\n";
    }
    if (per_time) {
        std::cout << "Using leading order of magnitude independently for each time step.\n";
    }
    if (corr_time > 0 && per_traj) {
         std::cout << "Filter leading order of magnitude by requiring at least "
                   << corr_time << " % of the ascent time.\n";
        if (ignore_phases) {
            std::cout << "Ignoring cloud phases.\n";
        }
    }

    std::cout << "For calculating the loom this program allocates about "
              << mem_usage/(1024*1024) << " MByte\n";
    std::cout << "Total usage: "
              << netcdf_reader.mem_usage/(1024*1024) + mem_usage/(1024*1024*1024) << " GByte.\n";

    std::cout << "Step 1/3: Get the leading order magnitudes and time limits.\n";
    auto t_first = std::chrono::system_clock::now();
    uint32_t n_files = 0;
    if (std::filesystem::is_directory(load_path)) {
        for (auto _ : std::filesystem::directory_iterator(load_path)) {
            n_files++;
        }
    }
    uint32_t counter = 0;
    if (std::filesystem::is_directory(load_path)) {
        for (auto const &dir_entry : std::filesystem::directory_iterator(load_path)) {
            counter++;
            // Load the file
            std::cout << "Loading " << dir_entry.path().c_str() << ".\n";
            netcdf_reader.load_vars(dir_entry.path().c_str());
            get_loom(
                netcdf_reader,
                loom,
                loom_tmp,
                time_limits,
                counter_arr,
                n_count_avg,
                n_trajectories,
                inflow_time,
                outflow_time,
                n_params,
                n_out_params,
                max_loom,
                per_traj,
                ignore_phases,
                n_phases,
                count_avg);
            netcdf_reader.close_netcdf();
            auto now = std::chrono::system_clock::now();
            double dt_total = ((std::chrono::duration<double>) (now - t_first)).count();
            std::string time = " in ";
            if (dt_total >= 3600)
                time = time + std::to_string(static_cast<int>(dt_total / 3600)) + "h ";
            if (dt_total >= 60)
                time = time + std::to_string(static_cast<int>(dt_total) % 3600 / 60) + "min ";

            time = time + std::to_string(static_cast<int>(dt_total) % 60) + "s ";
            std::cout << std::fixed << std::setprecision(3)
                      << "Done " << counter << " files in " << time;
            double remaining = dt_total / counter * (n_files - counter);
            time = "";
            if (remaining >= 3600)
                time = time + std::to_string(static_cast<int>(remaining / 3600)) + "h ";
            if (remaining >= 60)
                time = time + std::to_string(static_cast<int>(remaining) % 3600 / 60) + "min ";
            time = time + std::to_string(static_cast<int>(remaining) % 60) + "s";
            std::cout << " remaining " << time << "\n";
        }
    } else {
        // Load the file
        netcdf_reader.load_vars(load_path.c_str());
        get_loom(
                netcdf_reader,
                loom,
                loom_tmp,
                time_limits,
                counter_arr,
                n_count_avg,
                n_trajectories,
                inflow_time,
                outflow_time,
                n_params,
                n_out_params,
                max_loom,
                per_traj,
                ignore_phases,
                n_phases,
                count_avg);
    }
    if (!per_traj && !max_loom) {
        get_means(
            loom,
            loom_tmp,
            counter_arr,
            n_count_avg,
            n_params,
            n_out_params,
            n_phases,
            count_avg);
    }
    std::cout << "Step 2/3: Calculate the loom histogram\n";
    const uint32_t n_times = (time_limits[1] - time_limits[0])/DELTA_TIMESTEP + 1;
    std::vector<uint64_t> loom_count(n_phases * n_times * n_out_params * n_params);
    std::vector<uint64_t> count(n_phases * n_times);
    traj_counter.resize(n_out_params * n_params * n_phases);
    mem_usage = sizeof(uint64_t) * (loom_count.size() + count.size());
    std::cout << "For calculating the histograms this program allocates about "
              << mem_usage/(1024*1024) << " MByte\n";
    mem_usage = sizeof(double) *
                (loom.size() + time_limits.size())
                + sizeof(uint64_t) * (count.size() + loom_count.size() + traj_counter.size());
    std::cout << "Total usage: "
              << netcdf_reader.mem_usage/(1024*1024) + mem_usage/(1024*1024*1024) << " GByte\n";
    counter = 0;
    t_first = std::chrono::system_clock::now();
    if (std::filesystem::is_directory(load_path)) {
        for (auto const &dir_entry : std::filesystem::directory_iterator(load_path)) {
            counter++;
            // Load the file
            std::cout << "Loading " << dir_entry.path().c_str() << ".\n";
            netcdf_reader.load_vars(dir_entry.path().c_str());
            get_loom_histogram(
                    netcdf_reader,
                    loom,
                    time_limits,
                    loom_count,
                    count,
                    traj_counter,
                    inflow_time,
                    outflow_time,
                    n_params,
                    n_out_params,
                    n_phases,
                    max_loom,
                    per_traj,
                    per_time,
                    ignore_phases,
                    corr_time);
            netcdf_reader.close_netcdf();
            auto now = std::chrono::system_clock::now();
            double dt_total = ((std::chrono::duration<double>) (now - t_first)).count();
            std::string time = " in ";
            if (dt_total >= 3600)
                time = time + std::to_string(static_cast<int>(dt_total / 3600)) + "h ";
            if (dt_total >= 60)
                time = time + std::to_string(static_cast<int>(dt_total) % 3600 / 60) + "min ";

            time = time + std::to_string(static_cast<int>(dt_total) % 60) + "s ";
            std::cout << std::fixed << std::setprecision(3)
                      << "Done " << counter << " files in " << time;
            double remaining = dt_total / counter * (n_files - counter);
            time = "";
            if (remaining >= 3600)
                time = time + std::to_string(static_cast<int>(remaining / 3600)) + "h ";
            if (remaining >= 60)
                time = time + std::to_string(static_cast<int>(remaining) % 3600 / 60) + "min ";
            time = time + std::to_string(static_cast<int>(remaining) % 60) + "s";
            std::cout << " remaining " << time << "\n";
        }
    } else {
        // Load the file
        netcdf_reader.load_vars(load_path.c_str());
        get_loom_histogram(
                netcdf_reader,
                loom,
                time_limits,
                loom_count,
                count,
                traj_counter,
                inflow_time,
                outflow_time,
                n_params,
                n_out_params,
                n_phases,
                max_loom,
                per_traj,
                per_time,
                ignore_phases,
                corr_time);
    }

    std::cout << "Step 3/3: Flushing the results. \n" << std::flush;
    auto start_flush = std::chrono::system_clock::now();
    std::vector<int> dimid(Loom_dim_idx::n_loomdims);
    int ncid = create_dims(store_path, n_phases, n_times, dimid, ignore_phases);

    std::vector<int> varid, loom_order_varid;
    varid.resize(Loom_dim_idx::n_loomdims);
    if (get_loom_floats) {
        loom_order_varid.resize(n_params*3 + 2);
    } else {
        loom_order_varid.resize(n_params*2 + 1);
    }
    define_vars(ncid, varid, dimid, loom_order_varid, n_params, max_loom, ignore_phases, get_loom_floats);

    set_attributes(
            ncid,
            max_loom,
            loom_order_varid,
            n_params,
            get_loom_floats,
            per_traj,
            per_time,
            ignore_phases,
            inflow_time,
            outflow_time,
            corr_time,
            n_trajectories);
    set_compression(ncid, loom_order_varid, ignore_phases);
    write_dim_values(ncid, varid, time_limits, n_times, ignore_phases);
    std::vector<size_t> startp, countp;
    startp.push_back(0);
    if (!ignore_phases) {
        startp.push_back(0);
        countp.push_back(n_phases);
    }
    countp.push_back(n_times);
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            loom_order_varid[n_params + n_params * get_loom_floats],
            startp.data(),
            countp.data(),
            count.data()));

    startp.push_back(0);
    countp.insert(countp.begin(), n_out_params);
    const uint64_t offset_params = n_out_params * n_phases * n_times;
    for (uint32_t i = 0; i < n_params; i++) {
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                loom_order_varid[i],
                startp.data(),
                countp.data(),
                &loom_count[offset_params * i]));
    }
    const uint64_t offset_traj_counter =  n_out_params * n_phases;
    for (uint32_t i = 0; i < n_params; i++) {
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                loom_order_varid[n_params + n_params * get_loom_floats + get_loom_floats + 1 + i],
                startp.data(),
                countp.data(),
                &traj_counter[offset_traj_counter * i]));
    }
    if (get_loom_floats) {
        for (uint32_t i = 0; i < n_params; i++) {
            SUCCESS_OR_DIE(
                nc_put_vara(
                    ncid,
                    loom_order_varid[i + n_params],
                    startp.data(),
                    countp.data(),
                    &loom[n_phases * n_out_params * (i + 1)]));
        }
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                loom_order_varid[2*n_params + 1],
                startp.data(),
                countp.data(),
                &loom[0]));
    }

    auto end_flush = std::chrono::system_clock::now();
    double dt_total = ((std::chrono::duration<double>) (end_flush - start_flush)).count();
    std::string time = "Done in ";
    if (dt_total >= 3600)
        time = time + std::to_string(static_cast<int>(dt_total / 3600)) + "h ";
    if (dt_total >= 60)
        time = time + std::to_string(static_cast<int>(dt_total) % 3600 / 60) + "min ";
    time = time + std::to_string(static_cast<int>(dt_total) % 60) + "s.\n";
}
