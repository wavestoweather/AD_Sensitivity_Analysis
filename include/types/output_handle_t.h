#pragma once

#include <array>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <sstream>
#include <string>
#include <vector>

#include "include/misc/error.h"
#include "include/types/model_constants_t.h"
#include "include/types/netcdf_reader_t.h"
#include "include/types/reference_quantities_t.h"

struct output_handle_t{
    // for txt files
    std::stringstream out_tmp;
    std::ofstream outfile;
    std::ofstream out_diff[num_comp];
    std::stringstream out_diff_tmp[num_comp];
    uint64_t n_snapshots; // number of buffered snapshots
    uint64_t flushed_snapshots;
    int ncid; // ID of output file
    int local_num_comp;
    int local_num_par;

    std::string filetype;
    // for netCDF files and a vector for each column
    // columns: output parameter + gradients + time_afer_ascent + type + flags
    // fast index: record
    // each array = one column
    // slow index: num_comp
#ifdef OUT_DOUBLE
    std::array<std::vector<double>, num_comp+num_par+4 > output_buffer;
#else
    std::array<std::vector<float>, num_comp+num_par+4 > output_buffer;
#endif
    std::array<std::vector<unsigned char>, 4 > output_buffer_flags;
    // std::array<std::vector<std::string>, 1 > output_buffer_str;
    std::array<std::vector<uint64_t>, 1 > output_buffer_int;
    /**
     * ID for dimensions of output file.
     */
    std::vector<int> dimid;
    /**
     * ID for variables of output file.
     */
    std::vector<int> varid;

    uint64_t n_ens;
    uint64_t n_trajs;
    uint64_t n_trajs_file;
    uint64_t traj;
    uint64_t ens;
    uint64_t total_snapshots;
    std::string filename;


    /**
     * Number of time steps to write in the output.
     */
    uint64_t num_time;
    /**
     * Number of ensembles in the output.
     */
    uint64_t num_ens;

    /**
     * Simulation mode determines the order of the output file and
     * the mpi access mode.
     */
    int simulation_mode;

    enum Dim_idx
    {
        out_param_dim,
        ensemble_dim,
        trajectory_dim,
        time_dim,
        n_dims
    };

    enum Var_idx
    {
        // this should be the same as the *_idx in constants.h
        pressure,
        temperature,
        ascent,
        sat,
        qc,
        qr,
        qv,
        ncloud,
        nrain,
        qi,
        nice,
        qs,
        nsnow,
        qg,
        ngraupel,
        qh,
        nhail,
        qi_out,
        qs_out,
        qr_out,
        qg_out,
        qh_out,
        lat_heat,
        lat_cool,
        ni_out,
        ns_out,
        nr_out,
        ng_out,
        nh_out,
        height,
        inactive,
        dep,
        sub,

        time,
        trajectory,
        ensemble,
        out_param,

        time_ascent,
        lat,
        lon,
        conv_400,
        conv_600,
        slan_400,
        slan_600,
        step,

        // We do not clutter the gradient values here
        // The index is given by n_vars + i
        n_vars
    };
    /**
     * Indices within the double buffer
     */
    enum Buffer_idx
    {
        // The first entries are already given via constants.h *_idx
        time_ascent_buf = num_comp,
        lat_buf = num_comp+1,
        lon_buf = num_comp+2,
        n_buffer = num_comp+3
        // We do not clutter the gradient values here
        // The index is given by n_buffer + i
    };

    output_handle_t();

    output_handle_t(
        const std::string filetype,
        const std::string filename,
        const model_constants_t &cc,
        const reference_quantities_t &ref_quant,
        const std::string in_filename,
        const uint32_t write_index,
        const uint32_t snapshot_index,
        const int &rank,
        const int &simulation_mode);

    void setup(
        const std::string filetype,
        const std::string filename,
        const model_constants_t &cc,
        const reference_quantities_t &ref_quant,
        const std::string in_filename,
        const uint32_t write_index,
        const uint32_t snapshot_index,
        const int &rank);

    /**
     * Reset the number of flushed and recorded snapshots and update
     * the current trajectory and ensemble id. This function is needed
     * for static scheduling.
     *
     * @param traj_id New trajectory id for the output
     * @param ens_id New ensemble id for the output
     */
    void reset(const uint32_t traj_id, const uint32_t ens_id);

    /**
     * Buffer the gradients. If snapshot_index > 1, then the average gradient
     * is saved.
     *
     * @param cc
     * @param y_diff
     * @param snpashot_index
     */
    void buffer_gradient(
        const model_constants_t &cc,
        const std::vector< std::array<double, num_par > >  &y_diff,
        const uint32_t snapshot_index);

    /**
     * Writes data either to a stringstream for txt files or to different vectors
     * for netCDF files.
     *
     *
     */
    void buffer(const model_constants_t &cc,
        const netcdf_reader_t &netcdf_reader,
        const std::vector<codi::RealReverse> &y_single_new,
        const std::vector< std::array<double, num_par > >  &y_diff,
        const uint32_t sub,
        const uint32_t t,
        const double time_new,
        const uint32_t traj_id,
        const uint32_t ensemble,
        const reference_quantities_t &ref_quant,
        const uint32_t snapshot_index);

    /**
     * Write the buffered data to disk.
     *
     * @param cc
     */
    void flush_buffer(const model_constants_t &cc);

    /**
     * Buffer the current data  (model state and gradients) and
     * flush it to disk if necessary (depending on the time step.).
     *
     * @param cc
     */
    void process_step(
        const model_constants_t &cc,
        const netcdf_reader_t &netcdf_reader,
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
};