#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include "ncType.h"
#include <netcdf>
#include <sstream>
#include <string>
#include <vector>

// #include "include/microphysics/constants.h"
#include "include/misc/error.h"
#include "include/types/model_constants_t.h"
#include "include/types/nc_parameters_t.h"
#include "include/types/reference_quantities_t.h"

using namespace netCDF;

struct output_handle_t{
    // for txt files
    std::stringstream out_tmp;
    std::ofstream outfile;
    std::ofstream out_diff[num_comp];
    std::stringstream out_diff_tmp[num_comp];
    uint64_t n_snapshots; // number of buffered snapshots
    uint64_t flushed_snapshots;

    std::string filetype;
    // for netCDF files and a vector for each column
    // columns: output parameter + gradients + time_afer_ascent + type + flags
    // fast index: record
    // each array = one column
    // slow index: num_comp
    std::array<std::vector<double>, num_comp+num_par+4 > output_buffer;
    std::array<std::vector<unsigned char>, 4 > output_buffer_flags;
    std::array<std::vector<std::string>, 1 > output_buffer_str;
    std::array<std::vector<uint64_t>, 1 > output_buffer_int;
    NcFile datafile;
    std::vector<NcDim> dim_vector;
    std::vector<NcVar> var_vector;
    uint64_t n_ens;
    uint64_t n_trajs;
    uint64_t total_snapshots;
    std::string filename;

    output_handle_t(
        const std::string filetype,
        const std::string filename,
        const uint64_t n_trajs,
        const uint64_t n_ens,
        const model_constants_t &cc,
        const reference_quantities_t &ref_quant,
        const std::string in_filename,
        const uint32_t write_index,
        const uint32_t snapshot_index);

    /**
     * Writes data either to a stringstream for txt files or to different vectors
     * for netCDF files.
     */
    void buffer(const model_constants_t &cc,
        const nc_parameters_t &nc_params,
        const std::vector<codi::RealReverse> &y_single_new,
        const std::vector< std::array<double, num_par > >  &y_diff,
        const uint32_t sub,
        const uint32_t t,
        const double time_new,
        const uint32_t traj_id,
        const uint32_t ensemble,
        const reference_quantities_t &ref_quant);

    /**
     * Write the buffered data to disk.
     */
    void flush_buffer();

    /**
     * Buffer the current data  (model state and gradients) and
     * flush it to disk if necessary (depending on the time step.).
     */
    void process_step(
        const model_constants_t &cc,
        const nc_parameters_t &nc_params,
        const std::vector<codi::RealReverse> &y_single_new,
        const std::vector< std::array<double, num_par > >  &y_diff,
        const uint32_t sub,
        const uint32_t t,
        const double time_new,
        const uint32_t traj_id,
        const uint32_t write_index,
        const uint32_t snapshot_index,
#ifdef MET3D
        const uint32_t ensemble,
#endif
        const bool last_step,
        const reference_quantities_t &ref_quant);

    /**
     * Receive the buffer from a different MPI process.
     */
    void receive_buffer();

    /**
     * Send the buffer to a different MPI process.
     */
    void send_buffer();
};