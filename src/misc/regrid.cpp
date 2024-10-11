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


// Calculate LOoM for every phase for a given trajectory
void get_loom_phases(
    const netcdf_simulation_reader_t &netcdf_reader,
    const double &outflow_time,
    std::vector<double> &loom,
    const uint32_t &n_params,
    const uint32_t &n_out_params,
    const uint32_t &n_phases,
    const uint64_t &buffer_offset,
    const uint64_t &start_time_idx) {
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
                // For every phase
                for (uint32_t phase_idx = 0; phase_idx < n_phases; phase_idx++) {
                    // Check if in phase
                    if (current_phase != phase_idx) continue;
                    auto val = fabs(netcdf_reader.buffer_sens[i * 3 + j][t + buffer_offset]);
                    if (val > loom[j * n_phases + phase_idx]) {
                        loom[j * n_phases + phase_idx] = val;
                    }
                }
            }
        }
    }
}

// Given a LOoM threshold, trajectory and a bin, get the LOoM count.
void get_loom_count(
        netcdf_simulation_reader_t &netcdf_reader,
        std::vector<double> &loom,
        std::vector<uint64_t> &loom_count,
        const uint32_t &n_phases,
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const uint64_t &buffer_idx,
        const uint64_t grid_idx,
        const uint64_t &loom_offset
        ) {
    auto current_phase = netcdf_reader.buffer[Par_idx::phase][buffer_idx];
    if (current_phase == 3) return;
    // For every model parameter
    for (uint32_t i = 0; i < n_params; i++) {
        // For every output parameter
        for (uint32_t j = 0; j < n_out_params; j++) {
            // For every phase
            for (uint32_t phase_idx = 0; phase_idx < n_phases; phase_idx++) {
                // Check if in phase
                if (current_phase != phase_idx) continue;

                auto val = netcdf_reader.buffer_sens[i * n_out_params + j][buffer_idx];
                if ((fabs(val) >= loom[j * n_phases + phase_idx]/10)
                        && (loom[j * n_phases + phase_idx] > 0))
                    loom_count[j*loom_offset + i*(loom_offset*n_out_params) + grid_idx]++;
            }
        }
    }
}

void process_file(
        netcdf_simulation_reader_t &netcdf_reader,
        std::vector<uint64_t> &counts,
        std::vector<uint64_t> &loom_count,
        const uint32_t &offset_lat,
        const uint32_t &offset_p,
        const uint32_t &offset_t,
        const uint32_t &n_phases,
        const double &inflow_time,
        const double &outflow_time,
        const std::vector<double> &times,
        const std::vector<double> &pressure_levels,
        const std::vector<double> &lons,
        const std::vector<double> &lats,
        const uint32_t &n_params,
        const uint32_t &n_out_params,
        const int &n_times,
        const bool &relative_time) {
    ProgressBar pbar = ProgressBar(netcdf_reader.n_trajectories, 50, "Trajectory", std::cout);
    const uint64_t loom_offset = n_times*(lons.size()-1)*(lats.size()-1)*(pressure_levels.size()-1);
    for (uint32_t traj = 0; traj < netcdf_reader.n_trajectories; traj++) {
        uint64_t n_data = 0;
        netcdf_reader.init_netcdf(traj);
        // Do the calculations
        uint32_t current_idx = 0;
        bool not_finished = true;
        uint64_t buffer_offset = (traj%netcdf_reader.n_traj_buffer)*netcdf_reader.read_time_buffer;
        auto start_time_idx = 0;
        for (uint32_t i = 0; i < 10; i++) {
            if (std::isnan(netcdf_reader.buffer[Par_idx::time_after_ascent][i + buffer_offset])) continue;
            if (inflow_time < netcdf_reader.buffer[Par_idx::time_after_ascent][i + buffer_offset]) {
                start_time_idx = 0;
            } else {
                auto rel_time = netcdf_reader.buffer[Par_idx::time_after_ascent][i + buffer_offset];
                start_time_idx = i + (inflow_time - rel_time) / DELTA_TIMESTEP;
                break;
            }
        }
        buffer_offset += start_time_idx;
        std::vector<double> loom(n_phases*n_out_params);
        get_loom_phases(
            netcdf_reader,
            outflow_time,
            loom,
            n_params,
            n_out_params,
            n_phases,
            buffer_offset,
            start_time_idx);

        int t_idx = -1;
        double current_time = DELTA_TIMESTEP*start_time_idx;
        for (const auto &t : times) {
            if (current_time >= t) {
                t_idx++;
            }
        }
        uint32_t outflow_counter = 0;
        bool ascent_started = false;
        std::vector<double> tmp_avg(n_params*3);
        while (not_finished) {
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
            n_data++;
            int idx = t_idx * offset_t + p_idx * offset_p + lat_idx * offset_lat + lon_idx;
            counts[idx]++;
            get_loom_count(
                    netcdf_reader,
                    loom,
                    loom_count,
                    n_phases,
                    n_params,
                    n_out_params,
                    current_idx + buffer_offset,
                    idx,
                    loom_offset);

            if (current_idx + start_time_idx >= netcdf_reader.n_timesteps_in) {
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
                break;
            }
            current_time += DELTA_TIMESTEP;
            current_idx++;
            if (current_time >= times[t_idx+1]) {
                t_idx++;
                if (t_idx == n_times || current_idx + start_time_idx >= netcdf_reader.n_timesteps_in) {
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
        std::vector<int>& varid_sens,
        std::vector<int>& varid_misc,
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
    // Sensitivities
    for (uint32_t i = 0; i < n_params; i++) {
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                ("Loom count " + loom_order_sens[i]).c_str(),
                NC_UINT64,
                5,
                &dimid[Grid_dim_idx::outp_dim_idx],
                &varid_sens[i]));
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
        const std::vector<uint64_t> &counts,
        const std::vector<uint64_t> &loom_count,
        const std::vector<uint64_t> &max_count,
        const std::vector<int> &varid_sens,
        const std::vector<int> &varid_misc,
        const uint32_t &n_model_params,
        const uint32_t &n_out_params,
        const uint32_t &n_times,
        const uint32_t &n_plevels,
        const uint32_t &n_bins) {
    std::vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);

    countp.push_back(n_times);      // time
    countp.push_back(n_plevels);    // pressure
    countp.push_back(n_bins);       // lat
    countp.push_back(n_bins);       // lon

    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid_misc[Misc_idx::counts],
            startp.data(),
            countp.data(),
            counts.data()));

    countp.insert(countp.begin(), n_out_params);
    startp.push_back(0);
    const auto offset_params = n_times * n_plevels * n_bins * n_bins * n_out_params;
    for (uint32_t i = 0; i < n_model_params; i++) {
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid_sens[i],
                startp.data(),
                countp.data(),
                &loom_count[offset_params * i]));
    }
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid_misc[Misc_idx::ranking],
            startp.data(),
            countp.data(),
            max_count.data()));
}


void parse(
        int argc,
        char** argv,
        uint32_t &buffer_size,
        uint32_t &traj_buffer_size,
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
    traj_buffer_size = 300;

    const std::vector<std::string> args(argv + 1, argv + argc);
    for (uint32_t i = 0; i < args.size(); i += 2) {
        auto arg = args[i];
        if (arg == "-b" || arg == "--buffer_size") {
            buffer_size = std::stoi(args[i + 1]);
        } else if (arg == "-t" || arg == "--traj_buffer") {
            traj_buffer_size = std::stoi(args[i + 1]);
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
        } else {
            std::cout << "This program regrids trajectories and counts for each parameter and grid point "
                << "how often it had been in the leading order of magnitude for a target parameter and phase  "
                << "per trajectory.\nPossible options are:\n"
                << "-h, --help: Display this message\n"
                << "-b, --buffer_size: Adjust the number of time steps to store in buffer. Should be enough "
                << "--traj_buffer: The number of trajectories to store in RAM at once.\n"
                << "to hold all time steps of a trajectory.\n"
                << "--inflow_time: Time in seconds to consider before the ascent starts.\n"
                << "--outflow_time: Time in seconds to consider after the ascent ends.\n"
                << "--min_time: Start time in seconds to consider.\n"
                << "--max_time: End time in seconds to consider.\n"
                << "-i, --input_path: Path to folder with NetCDF-files or to a single NetCDF-file from a "
                << "sensitivity simulation.\n"
                << "-o, --output_path: Path and name (*.nc) to store the final data.\n"
                << "--min_lat: Minimum degree of latitude for the grid. Values outside of this are ignored.\n"
                << "--max_lat: Maximum degree of latitude for the grid. Values outside of this are ignored.\n"
                << "--min_lon: Minimum degree of longitude for the grid. Values outside of this are ignored.\n"
                << "--max_lon: Maximum degree of longitude for the grid. Values outside of this are ignored.\n"
                << "-r, --relative_time: Create a grid with time relative to the start of the ascent.\n"
                << "-n, --n_bins: The number of bins along longitude and latitude.\n"
                << "-d, --delta_t: The time in seconds for each bin.\n";
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
    uint32_t buffer_size, traj_buffer_size;
    double inflow_time, outflow_time;
    uint32_t n_bins;
    double min_time, max_time, delta_time;
    double min_lon, max_lon, min_lat, max_lat;
    bool relative_time;

    parse(
            argc,
            argv,
            buffer_size,
            traj_buffer_size,
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
    const uint32_t n_params = loom_order_sens.size();

    netcdf_simulation_reader_t netcdf_reader(buffer_size, false, traj_buffer_size);
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
    const uint32_t n_phases = 3;
    std::vector<uint64_t> counts(n_bins*n_bins*n_times*n_plevels);
    std::vector<uint64_t> loom_count(n_bins*n_bins*n_times*n_plevels*n_params*n_out_params);
    std::vector<uint64_t> max_count(n_bins*n_bins*n_times*n_plevels*n_out_params);

    double mem_usage = sizeof(double) * (n_phases * n_out_params  // for intermediate values
            + lons.size() + lats.size() + pressure_levels.size() + times.size())
            + sizeof(uint64_t) * (loom_count.size()+max_count.size());
    std::cout << "For calculating the gridded data this program allocates about "
              << mem_usage/(1024*1024*1024) << " GByte\n";
    std::cout << "Total usage: "
              << netcdf_reader.mem_usage/(1024*1024) + mem_usage/(1024*1024*1024) << " GByte\n";

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

    if (std::filesystem::is_directory(load_path)) {
        for (auto const &dir_entry : std::filesystem::directory_iterator(load_path)) {
            counter += 1;
            // Load the file
            std::cout << "Loading " << dir_entry.path().c_str() << "\n";
            netcdf_reader.load_vars(dir_entry.path().c_str());

            process_file(
                    netcdf_reader,
                    counts,
                    loom_count,
                    offset_lat,
                    offset_p,
                    offset_t,
                    n_phases,
                    inflow_time,
                    outflow_time,
                    times,
                    pressure_levels,
                    lons,
                    lats,
                    n_params,
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
        std::cout << "Loading " << load_path.c_str() << "\n";
        netcdf_reader.load_vars(load_path.c_str());

        process_file(
                netcdf_reader,
                counts,
                loom_count,
                offset_lat,
                offset_p,
                offset_t,
                n_phases,
                inflow_time,
                outflow_time,
                times,
                pressure_levels,
                lons,
                lats,
                n_params,
                n_out_params,
                n_times,
                relative_time);
    }

    // Calculate top params
    std::cout << "Calculate the top parameters\n";
    auto start_means = std::chrono::system_clock::now();
//    const uint64_t grid_size = n_times*(lons.size()-1)*(lats.size()-1)*(pressure_levels.size()-1);
    const uint64_t grid_size = n_times * n_plevels * n_bins * n_bins;
    std::fill(max_count.begin(), max_count.end(), n_params);
    for (uint64_t i = 0; i < counts.size(); i++) {
//        if (counts[i] == 0) continue;
        std::vector<uint64_t> max_val(n_out_params);
        for (uint64_t param = 0; param < n_params; param++) {
            for (uint32_t j = 0; j < n_out_params; j++) {
                auto current_val = loom_count[i + j*grid_size + param*(grid_size*n_out_params)];
                if (current_val > max_val[j]) {
                    max_count[i + j*grid_size] = param;
                    max_val[j] = current_val;
                }
            }
        }
    }

//    for (uint32_t j = 0; j < n_out_params; j++) {
//        for (uint64_t i = 0; i < counts.size(); i++) {
//            uint64_t max_val = 0;
//            uint64_t max_idx = n_params;
//            for (uint64_t param = 0; param < n_params; param++) {
//                auto current_val = loom_count[i + j*grid_size + param*grid_size*n_out_params];
//                if (current_val > max_val) {
//                    max_val = current_val;
//                    max_idx = param;
//                }
//            }
//            if (max_idx < n_params) {
//                max_count.push_back(loom_order_sens[max_idx]);
//            } else {
//                max_count.push_back("None");
//            }
//        }
//    }


    auto end_means = std::chrono::system_clock::now();
    double dt_total = ((std::chrono::duration<double>) (end_means - start_means)).count();
    std::string time = "Done in ";
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
    std::vector<int> varid, varid_sens;
    varid.resize(Grid_dim_idx::n_griddims);
    std::vector<int> varid_misc(2);
    varid_sens.resize(n_params);
    define_vars(ncid, varid, dimid,
                varid_sens,
                varid_misc, n_params);
    set_attributes(ncid);
    set_compression(ncid, varid_sens);
    set_compression(ncid, varid_misc);
    write_dim_values(ncid, varid, times, pressure_levels, lons, lats);
    write_other_values(
            ncid,
            counts,
            loom_count,
            max_count,
            varid_sens,
            varid_misc,
            n_params,
            n_out_params,
            n_times,
            n_plevels,
            n_bins);

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
