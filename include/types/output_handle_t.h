#pragma once

#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "include/misc/error.h"
#include "include/types/model_constants_t.h"
#include "include/types/netcdf_reader_t.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/segment_t.h"

struct output_handle_t{
    uint64_t n_snapshots;  // number of buffered snapshots
    uint64_t flushed_snapshots;
    int ncid;  // ID of output file
    int local_num_comp;
    int local_num_par;
    int local_ic_par;
    // Tracking either initial conditions or model parameters.
    // Both at the same time is not possible (or rather would take too long)
    bool track_ic;
#ifdef COMPRESS_OUTPUT
    int n_processes;
#endif
#ifdef OUT_DOUBLE
    const double FILLVALUE = std::numeric_limits<double>::quiet_NaN();
    const nc_type NC_FLOAT_T = NC_DOUBLE;
#else
    const float FILLVALUE = std::numeric_limits<float>::quiet_NaN();
    const nc_type NC_FLOAT_T = NC_FLOAT;
#endif

    std::string filetype;
    // for netCDF files and a vector for each column
    // columns: output parameter + gradients + time_after_ascent + type + flags
    // fast index: record
    // each array = one column
    // slow index: num_comp
#ifdef OUT_DOUBLE
    std::array<std::vector<double>, num_comp+num_par+5+static_cast<uint32_t>(Init_cons_idx::n_items) > output_buffer;
#else
    std::array<std::vector<float>, num_comp+num_par+5+static_cast<uint32_t>(Init_cons_idx::n_items) > output_buffer;
#endif
#if !defined B_EIGHT
    std::array<std::vector<unsigned char>, 5 > output_buffer_flags;
#else
    std::array<std::vector<unsigned char>, 1 > output_buffer_flags;
#endif
    std::array<std::vector<uint64_t>, 2 > output_buffer_int;
    /**
     * ID for dimensions of output file.
     */
    std::vector<int> dimid;
    /**
     * ID for variables of output file.
     */
    std::vector<int> varid;

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

    /**
     * For simulation mode create_train_set:
     * Idx in the output array for each parameter that can potentially be perturbed.
     * Values of -1 indicate that an output index has not been set.
     */
    std::vector<int> perturbed_idx;
    std::vector<std::string> perturbed_names;
    uint64_t n_perturbed_params;
#ifdef OUT_DOUBLE
    std::vector<double> unperturbed_vals;
#else
    std::vector<float> unperturbed_vals;
#endif
    enum Dim_idx {
        perturb_param_dim,
        out_param_dim,
        ensemble_dim,
        trajectory_dim,
        time_dim,
        n_dims
    };

    enum Var_idx {
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
        eff,

        time,
        trajectory,
        ensemble,
        out_param,

        time_ascent,
        lat,
        lon,
        asc600,
#if !defined(B_EIGHT)
        conv_400,
        conv_600,
        slan_400,
        slan_600,
#endif
        step,
        phase,

        perturbed,  // the dimension
        perturbation_value,
//        perturbed_param,
        // We do not clutter the gradient values here
        // The index is given by n_vars + i
        n_vars
    };
    /**
     * Indices within the double buffer
     */
    enum Buffer_idx {
        // The first entries are already given via constants.h *_idx
        time_ascent_buf = num_comp,
        lat_buf = num_comp+1,
        lon_buf = num_comp+2,
        perturb_buf = num_comp+3,
        n_buffer = num_comp+4
        // We do not clutter the gradient values here
        // The index is given by n_buffer + i
    };

    output_handle_t();

    template<class float_t>
    output_handle_t(
        const std::string filetype,
        const std::string filename,
        const model_constants_t<float_t> &cc,
        const std::string in_filename,
        const uint32_t write_index,
        const uint32_t snapshot_index,
        const int &rank,
        const int &simulation_mode,
        const bool &initial_cond,
#ifdef COMPRESS_OUTPUT
        const int &n_processes,
#endif
        const double delay_out_time = 0);

    template<class float_t>
    output_handle_t(
        const std::string filetype,
        const std::string filename,
        const model_constants_t<float_t> &cc,
        const std::string in_filename,
        const uint32_t write_index,
        const uint32_t snapshot_index,
        const int &rank,
        const int &simulation_mode,
        const bool &initial_cond,
#ifdef COMPRESS_OUTPUT
        const int &n_processes,
#endif
        const double delay_out_time,
        const std::vector<segment_t> &segments);

    ~output_handle_t() {nc_close(ncid);}

    template<class float_t>
    void setup(
        const std::string filetype,
        const std::string filename,
        const model_constants_t<float_t> &cc,
        const std::string in_filename,
        const uint32_t write_index,
        const uint32_t snapshot_index,
        const int &rank,
        const double delay_out_time = 0);

    /**
     * Reset the number of flushed and recorded snapshots and update
     * the current trajectory and ensemble id. This function is needed
     * for static scheduling.
     *
     * @param traj_id New trajectory id for the output
     * @param ens_id New ensemble id for the output
     */
    void reset(
        const uint32_t traj_id,
        const uint32_t ens_id,
        const uint64_t n_flushed = 0);

    /**
     * Buffer the gradients. If snapshot_index > 1, then the average gradient
     * is saved.
     *
     * @param cc
     * @param y_diff
     * @param snpashot_index
     */
    template<class float_t>
    void buffer_gradient(
        const model_constants_t<float_t> &cc,
        const std::vector< std::array<double, num_par > > &y_diff,
        const uint32_t snapshot_index);

    /**
     * Writes data either to a stringstream for txt files or to different vectors
     * for netCDF files.
     *
     *
     */
    template<class float_t>
    void buffer(const model_constants_t<float_t> &cc,
        const netcdf_reader_t &netcdf_reader,
        const std::vector<float_t> &y_single_new,
        const std::vector< std::array<double, num_par > > &y_diff,
        const uint32_t sub,
        const uint32_t t,
        const reference_quantities_t &ref_quant,
        const uint32_t snapshot_index,
        const bool previous_step_with_warm,
        const bool previous_step_with_ice = false);

    /**
     * Write the buffered data to disk or call certain functions for collective
     * access of the data if needed.
     *
     * @param cc
     * @param no_flush If true: do not really flush anything and just call
     *                 the flushing routines, which is needed if compression
     *                 is enabled.
     */
    template<class float_t>
    bool flush_buffer(const model_constants_t<float_t> &cc, bool no_flush = false);

    /**
     * Buffer the current data  (model state and gradients) and
     * flush it to disk if necessary (depending on the time step.).
     *
     * @param cc
     */
    template<class float_t>
    void process_step(
        const model_constants_t<float_t> &cc,
        const netcdf_reader_t &netcdf_reader,
        const std::vector<float_t> &y_single_new,
        const std::vector< std::array<double, num_par > > &y_diff,
        const uint32_t sub,
        const uint32_t t,
        const uint32_t write_index,
        const uint32_t snapshot_index,
        const bool last_step,
        const reference_quantities_t &ref_quant,
        const bool previous_step_with_warm,
        const bool previous_step_with_ice = false);

 private:
    /**
     * Setup and define variables for gradients in the NetCDF-file. This method sets
     * mode create_train_set: only an Out_Parameter_ID, a trajectory, and a time dimension for
     * the gradients.
     * mode limited_time_ens: only an Out_Parameter_ID, a trajectory, and a time dimension for
     * the gradients
     * else: an ensemble, trajectory and time dimension and an Out_Parameter_ID dimension if needed for
     * the gradients. Calls define_var_gradients().
     *
     * @param cc
     */
    template<class float_t>
    void setup_gradients(const model_constants_t<float_t> &cc);

    /**
     * Define the variables for gradients in the NetCDF-file.
     *
     * @param cc
     * @param dim_pointer Pointer to vector with all dim_ids.
     * @param n_dims Number of dimensions (=length of vector dim_pointer is pointing to).
     */
    template<class float_t>
    void define_var_gradients(
        const model_constants_t<float_t> &cc,
        const int *dim_pointer,
        const int &n_dims);

    /**
     * Define all variables for the NetCDF-file.
     *
     * @param cc
     */
    template<class float_t>
    void define_vars(const model_constants_t<float_t> &cc);

    template<class float_t>
    void set_attributes(
        const model_constants_t<float_t> &cc,
        const std::string in_filename);

    /**
     * In theory, one can apply compression
     * but this needs HDF5 >=1.10.2 (and Lustre).
     * Needs netcdf-c >= 4.7.4
     * All variables must be written in collective mode.
     * Perturbed simulations do not work with compression enabled!
     *
     * @param cc
     */
    template<class float_t>
    void set_compression(const model_constants_t<float_t> &cc);

    /**
     * Write the values for the dimensions.
     *
     * @param cc
     * @param delay_out_time
     */
    template<class float_t>
    void write_dimension_values(
        const model_constants_t<float_t> &cc,
        const double delay_out_time);

    /**
     * Open the recently created output file and prepare it for parallel writes.
     *
     * @param cc
     * @param file_string
     */
    template<class float_t>
    void set_parallel_access(
        const model_constants_t<float_t> &cc,
        const std::string file_string);
};
