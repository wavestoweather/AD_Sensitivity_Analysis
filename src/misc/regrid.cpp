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
#include "include/misc/constants_regrid.h"
#include "include/misc/error.h"
#include "include/misc/pbar.h"

const double FILLVALUE = std::numeric_limits<double>::quiet_NaN();
const nc_type NC_FLOAT_T = NC_DOUBLE;

void get_rank(
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const int &t_idx,
        const int &p_idx,
        const int &lat_idx,
        const int &lon_idx,
        const uint32_t &offset_lat,
        const uint32_t &offset_p,
        const uint32_t &offset_t,
        const uint32_t &offset_params,
        const uint32_t &rank_offset,
        std::vector<double> &tmp_avg,
        std::vector<uint64_t> &ranking) {
    std::vector<int> max_idx = {-1, -1, -1};
    std::vector<double> max_sens = {0, 0, 0};
    for (uint32_t i = 0; i < n_params; i++) {
        for (uint32_t j = 0; j < n_out_params; j++) {
            if (tmp_avg[j*n_params + i] > max_sens[j]) {
                max_sens[j] = tmp_avg[j*n_params + i];
                max_idx[j] = i;
            }
        }
    }
    uint32_t idx = t_idx * offset_t + p_idx * offset_p + lat_idx * offset_lat + lon_idx;
    for (uint32_t j = 0; j < n_out_params; j++) {
        if (max_idx[j] < 0) continue;
        ranking[idx + j*offset_params + max_idx[j]*rank_offset] += 1;
    }
    for (auto &v : tmp_avg) v = 0;
}

void process_file(
        netcdf_simulation_reader_t &netcdf_reader,
        std::vector<double> &means_model,
        std::vector<double> &mins_model,
        std::vector<double> &maxs_model,
        std::vector<double> &std_model,
        std::vector<uint64_t> &counts,
        std::vector<uint64_t> &ranking,
        std::vector<double> &means,
        std::vector<double> &mins,
        std::vector<double> &maxs,
        std::vector<double> &stds,
        const uint32_t &offset_lat,
        const uint32_t &offset_p,
        const uint32_t &offset_t,
        const uint32_t &offset_params,
        const uint32_t &offset_sens,
        const uint32_t &rank_offset,
        const double &inflow_time,
        const double &outflow_time,
        const std::vector<double> &times,
        const std::vector<double> &pressure_levels,
        const std::vector<double> &lons,
        const std::vector<double> &lats,
        const uint32_t &counter,
        const uint32_t &n_params,
        const uint32_t &n_model_params,
        const uint32_t &n_out_params,
        const int &n_times,
        const bool &relative_time) {
    ProgressBar pbar = ProgressBar(netcdf_reader.n_trajectories, 50, "Trajectory", std::cout);
    for (uint32_t traj = 0; traj < netcdf_reader.n_trajectories; traj++) {
        bool got_data = false;
        uint64_t n_data = 0;
        netcdf_reader.init_netcdf(traj);
        // Do the calculations
        uint32_t current_idx = 0;
        bool not_finished = true;
        uint64_t buffer_offset = (traj%netcdf_reader.n_traj_buffer)*netcdf_reader.read_time_buffer;
        auto start_time_idx = 0;
        for (uint32_t i = (traj%netcdf_reader.n_traj_buffer)*10; i < ((traj%netcdf_reader.n_traj_buffer)+1)*10; i++) {
            if (std::isnan(netcdf_reader.relative_time_buffer[i])) continue;
            if (inflow_time < netcdf_reader.relative_time_buffer[i]) {
                start_time_idx = 0;
            } else {
                auto rel_time = netcdf_reader.relative_time_buffer[i];
                start_time_idx = i-(traj%netcdf_reader.n_traj_buffer)*10
                                 + (inflow_time - rel_time) / DELTA_TIMESTEP;

                break;
            }
        }
        buffer_offset += start_time_idx;

        int t_idx = -1;
        double current_time = netcdf_reader.init_time + DELTA_TIMESTEP*start_time_idx;
        for (const auto &t : times) {
            if (current_time >= t) {
                t_idx++;
            }
        }
        uint32_t outflow_counter = 0;
        bool ascent_started = false;
        std::vector<double> tmp_avg(n_params*3);
        int old_t_idx = -1;
        int old_p_idx = -1;
        int old_lon_idx = -1;
        int old_lat_idx = -1;
        int old_idx = -1;
        bool got_rank = false;
        while (not_finished) {
            got_rank = false;
            double current_lon;
            double current_lat;
            if (relative_time) {
                current_lon = netcdf_reader.buffer[Par_idx::rel_lon][current_idx + buffer_offset];
                current_lat = netcdf_reader.buffer[Par_idx::rel_lat][current_idx + buffer_offset];
            } else {
                current_lon = netcdf_reader.buffer[Par_idx::lon][current_idx + buffer_offset];
                current_lat = netcdf_reader.buffer[Par_idx::lat][current_idx + buffer_offset];
            }
            double current_p = netcdf_reader.buffer[Par_idx::pressure][current_idx + buffer_offset];
            int lon_idx = -1;
            int lat_idx = -1;
            int p_idx = -1;
            for (const auto &p : pressure_levels) {
                if (current_p >= p) {
                    p_idx++;
                }
            }
            for (const auto &l : lons) {
                if (current_lon >= l) {
                    lon_idx++;
                }
            }
            for (const auto &l : lats) {
                if (current_lat >= l) {
                    lat_idx++;
                }
            }
            if (lon_idx < 0 || lat_idx < 0 || p_idx < 0
                || lon_idx >= lons.size()-1 || lat_idx >= lats.size()-1 || p_idx >= pressure_levels.size()-1) {
                current_time += DELTA_TIMESTEP;
                current_idx++;
                if (current_idx + start_time_idx >= netcdf_reader.n_timesteps_in) {
                    break;
                }
                if (current_time >= times[t_idx+1]) {
                    t_idx++;
                }
                if (t_idx == n_times) {
                    break;
                }
                continue;
            }
            got_data = true;
            n_data++;
            int idx = t_idx * offset_t + p_idx * offset_p + lat_idx * offset_lat + lon_idx;
            if (old_idx > -1 && old_idx != idx) {
                get_rank(
                        n_params,
                        n_out_params,
                        old_t_idx,
                        old_p_idx,
                        old_lat_idx,
                        old_lon_idx,
                        offset_lat,
                        offset_p,
                        offset_t,
                        offset_params,
                        rank_offset,
                        tmp_avg,
                        ranking);
                old_idx = idx;
                old_t_idx = t_idx;
                old_p_idx = p_idx;
                old_lat_idx = lat_idx;
                old_lon_idx = lon_idx;
            } else if (old_idx == -1) {
                old_idx = idx;
                old_t_idx = t_idx;
                old_p_idx = p_idx;
                old_lat_idx = lat_idx;
                old_lon_idx = lon_idx;
            }
            for (uint32_t i = 1; i < n_model_params+1; i++) {
                auto val = netcdf_reader.buffer[i][current_idx + buffer_offset];
                means_model[idx + offset_params*(i-1)] += val;
                std_model[idx + offset_params*(i-1)] += val*val;
                if (counter != 1 || traj != 0) {
                    if (val < mins_model[idx + offset_params*(i-1)])
                        mins_model[idx + offset_params*(i-1)] = val;
                    if (val > maxs_model[idx + offset_params*(i-1)])
                        maxs_model[idx + offset_params*(i-1)] = val;
                } else {
                    mins_model[idx + offset_params*(i-1)] = val;
                    maxs_model[idx + offset_params*(i-1)] = val;
                }
            }
            idx = t_idx * offset_t + p_idx * offset_p + lat_idx * offset_lat + lon_idx;
            for (uint32_t i = 0; i < n_params; i++) {
                for (uint32_t j = 0; j < n_out_params; j++) {
                    auto val = netcdf_reader.buffer_sens[i*3 + j][current_idx + buffer_offset];
                    means[idx + i*offset_sens + j*offset_params] += val;
                    stds[idx + i*offset_sens + j*offset_params] += val*val;
                    tmp_avg[j*n_params + i] += std::fabs(val);
                    if (counter != 1 || traj != 0) {
                        if (val < mins[idx + i*offset_sens + j*offset_params])
                            mins[idx + i*offset_sens + j*offset_params] = val;
                        if (val > maxs[idx + i*offset_sens + j*offset_params])
                            maxs[idx + i*offset_sens + j*offset_params] = val;
                    } else {
                        mins[idx + i*offset_sens + j*offset_params] = val;
                        maxs[idx + i*offset_sens + j*offset_params] = val;
                    }
                }
            }
            idx = t_idx * offset_t + p_idx * offset_p + lat_idx * offset_lat + lon_idx;
            counts[idx]++;
            // check if the next time step would be a new time idx and then
            // get the ranks
            if (DELTA_TIMESTEP + current_time >= times[t_idx+1] && t_idx-1 == n_times) {
                if (!got_rank) {
                    get_rank(
                            n_params,
                            n_out_params,
                            t_idx,
                            p_idx,
                            lat_idx,
                            lon_idx,
                            offset_lat,
                            offset_p,
                            offset_t,
                            offset_params,
                            rank_offset,
                            tmp_avg,
                            ranking);
                }
                break;
            }
            if (current_idx + start_time_idx >= netcdf_reader.n_timesteps_in) {
                if (!got_rank) {
                    get_rank(
                            n_params,
                            n_out_params,
                            t_idx,
                            p_idx,
                            lat_idx,
                            lon_idx,
                            offset_lat,
                            offset_p,
                            offset_t,
                            offset_params,
                            rank_offset,
                            tmp_avg,
                            ranking);
                }
                break;
            }

            // Finished if outflow limit is reached
            if (ascent_started) {
                if (!netcdf_reader.buffer[Par_idx::asc600][current_idx + buffer_offset])
                    outflow_counter += 1;
            } else {
                if (netcdf_reader.buffer[Par_idx::asc600][current_idx + buffer_offset])
                    ascent_started = true;
            }
            if (outflow_counter * DELTA_TIMESTEP >= outflow_time) {
                if (!got_rank) {
                    get_rank(
                            n_params,
                            n_out_params,
                            t_idx,
                            p_idx,
                            lat_idx,
                            lon_idx,
                            offset_lat,
                            offset_p,
                            offset_t,
                            offset_params,
                            rank_offset,
                            tmp_avg,
                            ranking);
                }
                break;
            }
            current_time += DELTA_TIMESTEP;
            current_idx++;
            if (current_time >= times[t_idx+1]) {
                t_idx++;
                if (t_idx == n_times || current_idx + start_time_idx >= netcdf_reader.n_timesteps_in) {
                    if (!got_rank) {
                        get_rank(
                                n_params,
                                n_out_params,
                                t_idx-1,
                                p_idx,
                                lat_idx,
                                lon_idx,
                                offset_lat,
                                offset_p,
                                offset_t,
                                offset_params,
                                rank_offset,
                                tmp_avg,
                                ranking);
                    }
                    break;
                }
            }
        }
        pbar.progress();
    }
}

int create_dims(
        const std::string& store_path,
        const uint32_t& n_plevels,
        const uint32_t& n_bins,
        const uint32_t& n_times,
        std::vector<int> &dimid) {
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
                    &dimid[Grid_dim_idx::outp_dim_idx]));   // idp
    SUCCESS_OR_DIE(
            nc_def_dim(
                    ncid,
                    "time",
                    n_times,
                    &dimid[Grid_dim_idx::time_griddim_idx]));
    SUCCESS_OR_DIE(
            nc_def_dim(
                    ncid,
                    "pressure",
                    n_plevels,
                    &dimid[Grid_dim_idx::p_dim_idx]));
    SUCCESS_OR_DIE(
            nc_def_dim(
                    ncid,
                    "lat",
                    n_bins,
                    &dimid[Grid_dim_idx::lat_dim_idx]));
    SUCCESS_OR_DIE(
            nc_def_dim(
                    ncid,
                    "lon",
                    n_bins,
                    &dimid[Grid_dim_idx::lon_dim_idx]));
    return ncid;
}

void define_vars(
        const int &ncid,
        std::vector<int>& varid,
        std::vector<int>& dimid,
        std::vector<int>& varid_means,
        std::vector<int>& varid_mins,
        std::vector<int>& varid_maxs,
        std::vector<int>& varid_stds,
        std::vector<int>& varid_means_sens,
        std::vector<int>& varid_mins_sens,
        std::vector<int>& varid_maxs_sens,
        std::vector<int>& varid_stds_sens,
        std::vector<int>& varid_misc,
        const uint32_t &n_model_params,
        const uint32_t &n_params) {
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,                   // ncid
            "Output_Parameter",     // name
            NC_STRING,              // type
            1,                      // ndims (0=scalar, 1=vector, 2=matrix, ...)
            &dimid[Grid_dim_idx::outp_dim_idx],     // dimid
            &varid[Grid_dim_idx::outp_dim_idx]));   // varid
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,
            "time",
            NC_FLOAT_T,
            1,
            &dimid[Grid_dim_idx::time_griddim_idx],
            &varid[Grid_dim_idx::time_griddim_idx]));
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,
            "pressure",
            NC_FLOAT_T,
            1,
            &dimid[Grid_dim_idx::p_dim_idx],
            &varid[Grid_dim_idx::p_dim_idx]));
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,
            "lat",
            NC_FLOAT_T,
            1,
            &dimid[Grid_dim_idx::lat_dim_idx],
            &varid[Grid_dim_idx::lat_dim_idx]));
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,
            "lon",
            NC_FLOAT_T,
            1,
            &dimid[Grid_dim_idx::lon_dim_idx],
            &varid[Grid_dim_idx::lon_dim_idx]));
    // Now all the other things
    // Model state
    for (uint32_t i = 1; i < n_model_params+1; i++) {
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Mean " + order_par[i]).c_str(),
                NC_FLOAT_T,
                4,
                &dimid[Grid_dim_idx::time_griddim_idx],
                &varid_means[i-1]));
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Min " + order_par[i]).c_str(),
                NC_FLOAT_T,
                4,
                &dimid[Grid_dim_idx::time_griddim_idx],
                &varid_mins[i-1]));
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Max " + order_par[i]).c_str(),
                NC_FLOAT_T,
                4,
                &dimid[Grid_dim_idx::time_griddim_idx],
                &varid_maxs[i-1]));
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Std " + order_par[i]).c_str(),
                NC_FLOAT_T,
                4,
                &dimid[Grid_dim_idx::time_griddim_idx],
                &varid_stds[i-1]));
    }
    // Sensitivities
    for (uint32_t i = 0; i < n_params; i++) {
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Mean " + order_sens[i]).c_str(),
                NC_FLOAT_T,
                5,
                &dimid[Grid_dim_idx::outp_dim_idx],
                &varid_means_sens[i]));
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Min " + order_sens[i]).c_str(),
                NC_FLOAT_T,
                5,
                &dimid[Grid_dim_idx::outp_dim_idx],
                &varid_mins_sens[i]));
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Max " + order_sens[i]).c_str(),
                NC_FLOAT_T,
                5,
                &dimid[Grid_dim_idx::outp_dim_idx],
                &varid_maxs_sens[i]));
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Std " + order_sens[i]).c_str(),
                NC_FLOAT_T,
                5,
                &dimid[Grid_dim_idx::outp_dim_idx],
                &varid_stds_sens[i]));
    }
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,
            order_misc[Misc_idx::ranking].c_str(),
            NC_UINT64,
            5,
            &dimid[Grid_dim_idx::outp_dim_idx],
            &varid_misc[Misc_idx::ranking]));
    SUCCESS_OR_DIE(
        nc_def_var(
            ncid,
            order_misc[Misc_idx::counts].c_str(),
            NC_UINT64,
            4,
            &dimid[Grid_dim_idx::time_griddim_idx],
            &varid_misc[Misc_idx::counts]));
}

void set_attributes(
        const int &ncid) {
    // Not now, maybe later
}

void write_dim_values(
        const int &ncid,
        const std::vector<int> &varid,
        const std::vector<double> &times,
        const std::vector<double> &pressure_levels,
        const std::vector<double> &lons,
        const std::vector<double> &lats) {
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
                varid[Grid_dim_idx::outp_dim_idx],
                startp.data(),
                &s));
        startp[0]++;
    }
    startp[0] = 0;

    countp.push_back(times.size()-1);
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Grid_dim_idx::time_griddim_idx],
            startp.data(),
            countp.data(),
            times.data()));
    countp[0] = pressure_levels.size()-1;
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Grid_dim_idx::p_dim_idx],
            startp.data(),
            countp.data(),
            pressure_levels.data()));
    countp[0] = lons.size()-1;
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Grid_dim_idx::lon_dim_idx],
            startp.data(),
            countp.data(),
            lons.data()));
    countp[0] = lats.size()-1;
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Grid_dim_idx::lat_dim_idx],
            startp.data(),
            countp.data(),
            lats.data()));
}

void set_compression(
        const int &ncid,
        const std::vector<int> &varid) {
    for (const auto v : varid)
        SUCCESS_OR_DIE(
            nc_def_var_deflate(
                ncid,
                v,
                1,  // shuffle
                1,  // deflate
                COMPRESSION_LEVEL));  // compression
}

void write_other_values(
        const int &ncid,
        const std::vector<double> &means,
        const std::vector<double> &mins,
        const std::vector<double> &maxs,
        const std::vector<double> &stds,
        const std::vector<int> &varid_means,
        const std::vector<int> &varid_mins,
        const std::vector<int> &varid_maxs,
        const std::vector<int> &varid_stds,
        const uint32_t &n_model_params,
        const uint32_t &offset_params,
        const std::vector<size_t> &startp,
        const std::vector<size_t> &countp) {
    for (uint32_t i = 0; i < n_model_params; i++) {
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid_means[i],
                startp.data(),
                countp.data(),
                &means[offset_params * i]));
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid_stds[i],
                startp.data(),
                countp.data(),
                &stds[offset_params * i]));
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid_mins[i],
                startp.data(),
                countp.data(),
                &mins[offset_params * i]));
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid_maxs[i],
                startp.data(),
                countp.data(),
                &maxs[offset_params * i]));
    }
}


void parse(
        int argc,
        char** argv,
        uint32_t &buffer_size,
        uint32_t &n_bins,
        double &min_time,
        double &max_time,
        double &delta_time,
        double &min_lon,
        double &max_lon,
        double &min_lat,
        double &max_lat,
        double &inflow_time,
        double &outflow_time,
        std::string &store_path,
        std::string &load_path,
        bool &relative_time) {
    if (argc > 64) {
        throw std::runtime_error("You provided more than 64 input parameters.");
    }
    // some default values
    buffer_size = 8600;
    n_bins = 100;
    inflow_time = 240;
    outflow_time = 240;
    min_time = 0;  // 270
    max_time = 60*60*24*3;  // alternatively for relative time: 255870
    delta_time = 60*60*24*3;  // alternatively: 3600; 60*60*24*3
    min_lon = -68;
    max_lon = 70;
    min_lat = 17;
    max_lat = 85;
    relative_time = false;

    const std::vector<std::string> args(argv + 1, argv + argc);
    for (auto i = 0; i < args.size(); i += 2) {
        auto arg = args[i];
        if (arg == "-b" || arg == "--buffer_size") {
            buffer_size = std::stoi(args[i+1]);
        } else if (arg == "-n" || arg == "--n_bins") {
            n_bins = std::stoi(args[i+1]);
        } else if (arg == "--min_time") {
            min_time = std::stod(args[i+1], nullptr);
        } else if (arg =="--max_time") {
            max_time = std::stod(args[i+1], nullptr);
        } else if (arg == "-d" || arg == "--delta_t") {
            delta_time = std::stod(args[i+1], nullptr);
        } else if (arg == "--min_lon") {
            min_lon = std::stod(args[i+1], nullptr);
        } else if (arg == "--max_lon") {
            max_lon = std::stod(args[i+1], nullptr);
        } else if (arg == "--min_lat") {
            min_lat = std::stod(args[i+1], nullptr);
        } else if (arg == "--max_lat") {
            max_lat = std::stod(args[i+1], nullptr);
        } else if (arg == "--inflow_time") {
            inflow_time = std::stod(args[i+1], nullptr)  * (-1);
        } else if (arg == "--outflow_time") {
            outflow_time = std::stod(args[i+1], nullptr);
        } else if (arg == "-i" || arg == "--input_path") {
            load_path = args[i+1];
        } else if (arg == "-o" || arg == "--output_path") {
            store_path = args[i+1];
        } else if (arg == "-r" || arg == "--relative_time") {
            relative_time = true;
            i--;
        }
    }
}


/*
 * Regrid trajectory data.
 *
 *
 */
int main(int argc, char** argv) {
    int rank;
    int n_processes;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string load_path_str;
    std::string store_path;
    uint32_t buffer_size;
    double inflow_time, outflow_time;
    uint32_t n_bins;
    double min_time, max_time, delta_time;
    double min_lon, max_lon, min_lat, max_lat;
    bool relative_time;

    parse(
            argc,
            argv,
            buffer_size,
            n_bins,
            min_time,
            max_time,
            delta_time,
            min_lon,
            max_lon,
            min_lat,
            max_lat,
            inflow_time,
            outflow_time,
            store_path,
            load_path_str,
            relative_time);

    const std::filesystem::path load_path = {load_path_str};
    const uint32_t n_params = order_sens.size();
    const uint32_t n_model_params = 18;

    netcdf_simulation_reader_t netcdf_reader(buffer_size);
    std::vector<double> pressure_levels = {
            0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000
    };
    const uint32_t n_plevels = pressure_levels.size()-1;
    std::vector<double> lons(n_bins+1);
    std::vector<double>  lats(n_bins+1);

    double delta_lon = (max_lon - min_lon)/n_bins;
    double delta_lat = (max_lat - min_lat)/n_bins;
    uint32_t counter = 0;
    for (auto &l : lons) {
        l = min_lon + counter*delta_lon;
        counter++;
    }
    counter = 0;
    for (auto &l : lats) {
        l = min_lat + counter*delta_lat;
        counter++;
    }
    int n_times = (max_time - min_time)/delta_time;
    counter = 0;
    std::vector<double> times(n_times+1);
    for (auto &t : times) {
        t = min_time + counter*delta_time;
        counter++;
    }
    const uint32_t n_out_params = 3;
    std::vector<double> means_model(n_bins*n_bins*n_times*n_plevels*n_model_params);
    std::vector<double> mins_model(n_bins*n_bins*n_times*n_plevels*n_model_params);
    std::vector<double> maxs_model(n_bins*n_bins*n_times*n_plevels*n_model_params);
    std::vector<double> std_model(n_bins*n_bins*n_times*n_plevels*n_model_params);
    std::vector<uint64_t> counts(n_bins*n_bins*n_times*n_plevels);

    std::vector<uint64_t> ranking(n_bins*n_bins*n_times*n_out_params*n_plevels*n_params);
    std::vector<uint64_t> ranking_vals(n_bins*n_bins*n_times*n_plevels*n_out_params);
    std::vector<double> means(n_bins*n_bins*n_times*n_out_params*n_params*n_plevels);
    std::vector<double> mins(n_bins*n_bins*n_times*n_out_params*n_params*n_plevels);
    std::vector<double> maxs(n_bins*n_bins*n_times*n_out_params*n_params*n_plevels);
    std::vector<double> stds(n_bins*n_bins*n_times*n_out_params*n_params*n_plevels);

    double mem_usage = sizeof(double) *
                       (means_model.size() + mins_model.size() + maxs_model.size()
                        + std_model.size() + means.size() + mins.size() + maxs.size()
                        + stds.size()) + sizeof(uint64_t) * (counts.size() + ranking.size()
                                                             + ranking_vals.size());
    std::cout << "For calculating the gridded data this program allocates about "
              << mem_usage/(1024*1024*1024) << " GByte\n";
    std::cout << "Total usage: "
              << netcdf_reader.mem_usage/(1024*1024) + mem_usage/(1024*1024*1024) << " GByte\n";
    // Get the counts, means, min. max, std, ranking
    // Get means:
    // Sum everything and then divide by counts in the end.
    // Get min/max:
    // Check for min/max.
    // Get standard deviation:
    // Increment the squares and in the end divide it by counts and subtract the squared means.
    // Get ranking:
    // For each grid cell: Get a list of all parameters and increment whenever a parameter
    // is number one.
    auto t_first = std::chrono::system_clock::now();
    uint32_t n_files = 0;
    if (std::filesystem::is_directory(load_path)) {
        for (auto _ : std::filesystem::directory_iterator(load_path)) {
            n_files++;
        }
    }
    counter = 0;

    const uint32_t offset_lat = n_bins;
    const uint32_t offset_p = n_bins*offset_lat;
    const uint32_t offset_t = n_plevels * offset_p;
    const uint32_t offset_params = offset_t * n_times;
    const uint32_t offset_sens = offset_params * n_out_params;
    const uint32_t rank_offset = offset_params * n_out_params;

    if (std::filesystem::is_directory(load_path)) {
        for (auto const &dir_entry : std::filesystem::directory_iterator(load_path)) {
            counter += 1;
            // Load the file
            std::cout << "Loading " << dir_entry.path().c_str() << "\n";
            netcdf_reader.load_vars(dir_entry.path().c_str());

            process_file(
                    netcdf_reader,
                    means_model,
                    mins_model,
                    maxs_model,
                    std_model,
                    counts,
                    ranking,
                    means,
                    mins,
                    maxs,
                    stds,
                    offset_lat,
                    offset_p,
                    offset_t,
                    offset_params,
                    offset_sens,
                    rank_offset,
                    inflow_time,
                    outflow_time,
                    times,
                    pressure_levels,
                    lons,
                    lats,
                    counter,
                    n_params,
                    n_model_params,
                    n_out_params,
                    n_times,
                    relative_time);
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
                      << "\nDone " << counter << " files in " << time;
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
        counter += 1;
        // Load the file
        netcdf_reader.load_vars(load_path.c_str());

        process_file(
                netcdf_reader,
                means_model,
                mins_model,
                maxs_model,
                std_model,
                counts,
                ranking,
                means,
                mins,
                maxs,
                stds,
                offset_lat,
                offset_p,
                offset_t,
                offset_params,
                offset_sens,
                rank_offset,
                inflow_time,
                outflow_time,
                times,
                pressure_levels,
                lons,
                lats,
                counter,
                n_params,
                n_model_params,
                n_out_params,
                n_times,
                relative_time);
    }

    // Calculate means and std
    std::cout << "Calculate the means and variances. ";
    auto start_means = std::chrono::system_clock::now();
    for (uint64_t i = 0; i < counts.size(); i++) {
        if (counts[i] == 0) continue;
        for (uint32_t j = 0; j < n_model_params; j++) {
            means_model[i + j * counts.size()] /= counts[i];
            std_model[i + j * counts.size()] /= counts[i];
        }
        for (uint32_t j = 0; j < n_params; j++) {
            means[i + j * counts.size() * n_out_params] /= counts[i];
            means[i + j * counts.size() * n_out_params + counts.size()] /= counts[i];
            means[i + j * counts.size() * n_out_params + 2 * counts.size()] /= counts[i];
            stds[i + j * counts.size() * n_out_params] /= counts[i];
            stds[i + j * counts.size() * n_out_params + counts.size()] /= counts[i];
            stds[i + j * counts.size() * n_out_params + 2 * counts.size()] /= counts[i];
        }
    }
    for (uint64_t i = 0; i < counts.size(); i++) {
        if (counts[i] == 0) continue;
        for (uint32_t j = 0; j < n_model_params; j++) {
            std_model[i + j * counts.size()] -= (means_model[i + j * counts.size()] *
                                                 means_model[i + j * counts.size()]);
        }
        for (uint32_t j = 0; j < n_params; j++) {
            stds[i + j * counts.size() * n_out_params] -= (means[i + j * counts.size() * 3] *
                                                           means[i + j * counts.size() * n_out_params]);
            stds[i + j * counts.size() * n_out_params + counts.size()] -= (
                    means[i + j * counts.size() * n_out_params + counts.size()] *
                    means[i + j * counts.size() * n_out_params + counts.size()]);
            stds[i + j * counts.size() * n_out_params + 2 * counts.size()] -= (
                    means[i + j * counts.size() * n_out_params + 2 * counts.size()] *
                    means[i + j * counts.size() * n_out_params + 2 * counts.size()]);
        }
    }
    auto end_means = std::chrono::system_clock::now();
    double dt_total = ((std::chrono::duration<double>) (end_means - start_means)).count();
    std::string time = "Done in ";
    if (dt_total >= 3600)
        time = time + std::to_string(static_cast<int>(dt_total / 3600)) + "h ";
    if (dt_total >= 60)
        time = time + std::to_string(static_cast<int>(dt_total) % 3600 / 60) + "min ";

    time = time + std::to_string(static_cast<int>(dt_total) % 60) + "s.\n";
    std::cout << std::fixed << std::setprecision(3) << time;

    // Get the rank
    std::cout << "Get the top parameters. ";
    auto start_rank = std::chrono::system_clock::now();
    for (uint32_t model_i = 0; model_i < n_out_params; model_i++) {
        for (int t_i = 0; t_i < n_times; t_i++) {
            for (uint32_t p_i = 0; p_i < n_plevels; p_i++) {
                for (uint32_t lat_i = 0; lat_i < n_bins; lat_i++) {
                    for (uint32_t lon_i = 0; lon_i < n_bins; lon_i++) {
                        uint32_t max_idx = order_sens.size();
                        uint32_t max_count = 0;
                        uint32_t this_idx = t_i * offset_t + p_i * offset_p + lat_i * offset_lat + lon_i;
                        for (uint32_t param_idx = 0; param_idx < n_params; param_idx++) {
                            if (ranking[this_idx + model_i * offset_params + param_idx * rank_offset] > max_count) {
                                max_idx = param_idx;
                                max_count = ranking[this_idx + model_i * offset_params + param_idx * rank_offset];
                            }
                        }
                        ranking_vals[this_idx + model_i * offset_params] = max_idx;
                    }
                }
            }
        }
    }
    auto end_rank = std::chrono::system_clock::now();
    dt_total = ((std::chrono::duration<double>) (end_rank - start_rank)).count();
    time = "Done in ";
    if (dt_total >= 3600)
        time = time + std::to_string(static_cast<int>(dt_total / 3600)) + "h ";
    if (dt_total >= 60)
        time = time + std::to_string(static_cast<int>(dt_total) % 3600 / 60) + "min ";

    time = time + std::to_string(static_cast<int>(dt_total) % 60) + "s.\n";
    std::cout << std::fixed << std::setprecision(3) << time;

    // flush results
    std::cout << "Flushing the results. " << std::flush;
    auto start_flush = std::chrono::system_clock::now();
    std::vector<int> dimid(Grid_dim_idx::n_griddims);
    int ncid = create_dims(store_path, n_plevels, n_bins, n_times, dimid);
    std::vector<int> varid, varid_means, varid_mins, varid_maxs, varid_stds;
    std::vector<int> varid_means_sens, varid_mins_sens, varid_maxs_sens, varid_stds_sens;
    varid.resize(Grid_dim_idx::n_griddims);
    std::vector<int> varid_misc(2);
    varid_means.resize(n_model_params);
    varid_mins.resize(n_model_params);
    varid_maxs.resize(n_model_params);
    varid_stds.resize(n_model_params);
    varid_means_sens.resize(n_params);
    varid_mins_sens.resize(n_params);
    varid_maxs_sens.resize(n_params);
    varid_stds_sens.resize(n_params);
    define_vars(ncid, varid, dimid, varid_means, varid_mins, varid_maxs, varid_stds,
                varid_means_sens, varid_mins_sens, varid_maxs_sens, varid_stds_sens,
                varid_misc, n_model_params, n_params);

    set_attributes(ncid);
    set_compression(ncid, varid_means);
    set_compression(ncid, varid_mins);
    set_compression(ncid, varid_maxs);
    set_compression(ncid, varid_stds);
    set_compression(ncid, varid_means_sens);
    set_compression(ncid, varid_mins_sens);
    set_compression(ncid, varid_maxs_sens);
    set_compression(ncid, varid_stds_sens);
    set_compression(ncid, varid_misc);
    write_dim_values(ncid, varid, times, pressure_levels, lons, lats);
    std::vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);

    countp.push_back(n_times);      // time
    countp.push_back(n_plevels);    // pressure
    countp.push_back(n_bins);       // lat
    countp.push_back(n_bins);       // lon

    write_other_values(
            ncid,
            means_model,
            mins_model,
            maxs_model,
            std_model,
            varid_means,
            varid_mins,
            varid_maxs,
            varid_stds,
            n_model_params,
            offset_params,
            startp,
            countp);

    // flush counts
    SUCCESS_OR_DIE(
            nc_put_vara(
                    ncid,
                    varid_misc[Misc_idx::counts],
                    startp.data(),
                    countp.data(),
                    counts.data()));

    std::vector<size_t> startp2, countp2;
    startp2.push_back(0);
    startp2.push_back(0);
    startp2.push_back(0);
    startp2.push_back(0);
    startp2.push_back(0);

    countp2.push_back(3);
    countp2.push_back(n_times);     // time
    countp2.push_back(n_plevels);   // pressure
    countp2.push_back(n_bins);      // lat
    countp2.push_back(n_bins);      // lon

    write_other_values(
            ncid,
            means,
            mins,
            maxs,
            stds,
            varid_means_sens,
            varid_mins_sens,
            varid_maxs_sens,
            varid_stds_sens,
            n_params,
            offset_sens,
            startp2,
            countp2);
    // flush ranking
    SUCCESS_OR_DIE(
            nc_put_vara(
                    ncid,
                    varid_misc[Misc_idx::ranking],
                    startp2.data(),
                    countp2.data(),
                    ranking_vals.data()));
    auto end_flush = std::chrono::system_clock::now();
    dt_total = ((std::chrono::duration<double>) (end_flush - start_flush)).count();
    time = "Done in ";
    if (dt_total >= 3600)
        time = time + std::to_string(static_cast<int>(dt_total / 3600)) + "h ";
    if (dt_total >= 60)
        time = time + std::to_string(static_cast<int>(dt_total) % 3600 / 60) + "min ";

    time = time + std::to_string(static_cast<int>(dt_total) % 60) + "s.\n";
    std::cout << std::fixed << std::setprecision(3) << time;
}