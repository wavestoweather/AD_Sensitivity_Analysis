#pragma once

#include <netcdf.h>
#include <netcdf_par.h>

#include <array>
#include <vector>

#include "codi.hpp"

#include "include/misc/error.h"
#include "include/types/model_constants_t.h"
#include "include/types/reference_quantities_t.h"


struct netcdf_reader_t {
    uint64_t n_ensembles;
    uint64_t n_trajectories;
    // double dlon, dlat;
    uint64_t start_time_idx; /*!< The initial start time index. */
    uint64_t time_idx; /*!< Current index to read from netcdf file. */

    /**
     *
     * @param buffer_size Number of time steps to read at once
     */
    explicit netcdf_reader_t(const uint32_t &buffer_size);

    /**
     * Read data from the buffer and store in in y_single_old and inflows.
     * If the buffer has been read completely, the next few datapoints are
     * loaded from the netCDF file via buffer_params().
     *
     * @param cc Model constants. On out: Added number of simulation steps
     * @param ref_quant Reference quantities used to change from netcdf units to
     *                  simulation units
     * @param y_single_old Vector with (possibly updated) model state variables
                        after return
     * @param inflows Vector with inflow from above after return
     * @param step Current iteration step
     * @param checkpoint_flag If true: do not update model state
     * @param start_over_env If true: update environment/thermodynamics
     */
    void read_buffer(
        model_constants_t &cc,
        const reference_quantities_t &ref_quant,
        std::vector<codi::RealReverse> &y_single_old,
        std::vector<codi::RealReverse> &inflows,
        const uint32_t &step,
        const bool &checkpoint_flag,
        const bool &start_over_env);

    /**
     * Set dimids of this object and set the maximum amount of trajectories and
     * ensembles in cc.
     *
     * @param input_file Path to file to read
     * @param cc Model constants with number of ensembles and trajectories on return
     * @param simulation_mode The simulation mode defined on the command line
     */
    void set_dims(
        const char *input_file,
        model_constants_t &cc,
        const int &simulation_mode);

    /**
    * Read initial values from the netcdf file and stores them to y_init.
    * Also stores the amount of trajectories in the input file and
    * several quantities to cc such as the number of steps to simulate.
    *
    *
    * @param input_file Path to input netcdf file as char array
    * @param checkpoint_flag Is this already initialized via a checkpoint?
    * @param cc Model constants. On out: Added number of simulation steps
    * @param simulation_mode Needed to define if parallel access to the file is
    *                   possible.
    * @param current_time Time in seconds from a checkpoint
    * @return Errorcode (0=no errors; 1=simulation breaking error)
    */
    void init_netcdf(
#ifdef MET3D
        double &start_time,
#endif
        const char *input_file,
        const bool &checkpoint_flag,
        model_constants_t &cc,
        const int &simulation_mode,
        const double current_time);

    void read_initial_values(
        std::vector<double> &y_init,
        const reference_quantities_t &ref_quant,
        model_constants_t &cc,
        const bool &checkpoint_flag,
        const uint64_t &traj_id,
        const uint64_t &ens_id);
    /**
     * Read values for initializing a simulation.
     *
     * @param y_init Array of num_comp many doubles
     * @param ref_quant Reference quantities used to change from netcdf units to
     *                  simulation units
     * @param checkpoint_flag Is this already initialized via a checkpoint?
     * @param cc Model constants. On out: Added number of simulation steps
     * @param checkpoint_flag Wether the data is already there due to loading a checkpoint
     */
    template<class float_t>
    void read_initial_values(
        std::vector<float_t> &y_init,
        const reference_quantities_t &ref_quant,
        model_constants_t &cc,
        const bool &checkpoint_flag);

    double get_lat(const uint32_t &t) const {return buffer[Par_idx::lat][(t)%n_timesteps_buffer];}
    double get_lon(const uint32_t &t) const {return buffer[Par_idx::lon][(t)%n_timesteps_buffer];}
#ifdef MET3D
    double get_relative_time(const uint32_t &t) const {
        return buffer[Par_idx::time_after_ascent][(t)%n_timesteps_buffer];
    }
    bool get_conv_400(const uint32_t &t) const {return (buffer[Par_idx::conv_400][t] > 0.5);}
    bool get_conv_600(const uint32_t &t) const {return (buffer[Par_idx::conv_600][t] > 0.5);}
    bool get_slan_400(const uint32_t &t) const {return (buffer[Par_idx::slan_400][t] > 0.5);}
    bool get_slan_600(const uint32_t &t) const {return (buffer[Par_idx::slan_600][t] > 0.5);}

#endif

 private:
    std::vector<int> startp, countp;
    int ncid;
    uint64_t n_timesteps_in; /*!< Total number of time steps that can be read from the input file. */
    uint32_t n_timesteps_buffer; /*!< Number of time steps that can be stored in the buffer. */
    uint64_t time_buffer_idx; /*!< Current index to read from the buffer. */

    uint32_t traj_idx; /*!< Index of trajectory to read from. */
    uint32_t ens_idx; /*!< Index of ensemble to read from. */
    bool already_open; /*!< Is the netCDF file already open? */

    /**
     * ID for dimensions of output file.
     */
    std::vector<int> dimid;
    /**
     * ID for variables of output file.
     */
    std::vector<int> varid, varid_once;

    enum Par_idx {
        pressure,
        temperature,
        ascent,
        lat,
        lon,
#if defined(FLUX) && !defined(WCB)
        qr_in,
        nr_in,
    #if defined(RK4ICE)
        qi_in,
        qs_in,
        qg_in,
        // qh_in,  // there is no dataset yet with hail in it
        ni_in,
        ns_in,
        ng_in,
        // nh_in,
    #endif
#endif
        height,
#ifdef MET3D
        time_after_ascent,
        conv_400,  // might make trouble type wise
        conv_600,
        slan_400,
        slan_600,
    #ifdef TURBULENCE
        q_turb,
    #endif
#endif
        n_pars
    };
    enum Dim_idx {
        time_dim_idx,
        trajectory_dim_idx,
        ensemble_dim_idx,
        n_dims
    };
    // Index for variable ids that are read during initialization
    enum Par_once_idx {
        qv,
        qc,
        qr,
        qi,
        qs,
        qg,
        nc,
        nr,
        ni,
        ns,
        ng,
        time,
        n_pars_once
    };
    std::array<std::vector<double>, Par_idx::n_pars > buffer;

    /**
     * Load variables from the netCDF file such that data can be loaded using
     * load_params()
     *
     */
    void load_vars();

    /**
     * Read data from the netCDF file and store it in the buffer.
     *
     * @param ref_quant Reference quantities used to change from netcdf units to
     *                  simulation units
     */
    void buffer_params(const reference_quantities_t &ref_quant);
};
