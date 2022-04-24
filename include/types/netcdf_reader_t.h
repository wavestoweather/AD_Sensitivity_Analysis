#pragma once

#include <netcdf.h>
#include <netcdf_par.h>

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "codi.hpp"

#include "include/misc/error.h"
#include "include/types/model_constants_t.h"
#include "include/types/reference_quantities_t.h"

struct netcdf_reader_t {
    uint64_t n_ensembles;
    uint64_t n_trajectories;
    // double dlon, dlat;
    bool start_time_idx_given; /*!< True if start_time_idx is given by the user. */
    uint64_t start_time_idx; /*!< The start time index. */
    uint64_t start_time_idx_original; /*!< The initial start time index. */
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
     *
     * @returns SUCCESS or with #def TRUSTED_DATA INPUT_NAN_ERR, otherwise
     *          throws an error.
     */
    template<class float_t>
    int read_buffer(
        model_constants_t<float_t> &cc,
        const reference_quantities_t &ref_quant,
        std::vector<float_t> &y_single_old,
        std::vector<float_t> &inflows,
        const uint32_t &step,
        const bool &checkpoint_flag,
        const bool &start_over_env);

    /**
     * Set dimids of this object and set the maximum amount of trajectories and
     * ensembles in cc. Sets pressure in reference quantities. Default is 1e5
     * for Pa, otherwise hPa with 1e3 is possible as well.
     *
     * @param input_file Path to file to read
     * @param cc Model constants with number of ensembles and trajectories on return
     * @param simulation_mode The simulation mode defined on the command line
     */
    template<class float_t>
    void set_dims(
        const char *input_file,
        model_constants_t<float_t> &cc,
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
    template<class float_t>
    void init_netcdf(
#ifdef MET3D
        double &start_time,
#endif
//        const char *input_file,
        const bool &checkpoint_flag,
        model_constants_t<float_t> &cc,
//        const int &simulation_mode,
        const double current_time,
        const reference_quantities_t &ref_quant);

    template<class float_t>
    void read_initial_values(
        std::vector<double> &y_init,
        const reference_quantities_t &ref_quant,
        model_constants_t<float_t> &cc,
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
        std::vector<double> &y_init,
        const reference_quantities_t &ref_quant,
        model_constants_t<float_t> &cc,
        const bool &checkpoint_flag);

    double get_lat(const uint32_t &t, const uint32_t &sub) const;
    double get_lon(const uint32_t &t, const uint32_t &sub) const;
#if defined(MET3D) || defined(B_EIGHT)
    double get_relative_time(const uint32_t &t) const {
        return buffer[Par_idx::time_after_ascent][t%buffer[Par_idx::time_after_ascent].size()];
    }
#if !defined(B_EIGHT)
    bool get_conv_400(const uint32_t &t) const {return (buffer[Par_idx::conv_400][
        t%buffer[Par_idx::conv_400].size()] > 0.5);}
    bool get_conv_600(const uint32_t &t) const {return (buffer[Par_idx::conv_600][
        t%buffer[Par_idx::conv_600].size()] > 0.5);}
    bool get_slan_400(const uint32_t &t) const {return (buffer[Par_idx::slan_400][
        t%buffer[Par_idx::slan_400].size()] > 0.5);}
    bool get_slan_600(const uint32_t &t) const {return (buffer[Par_idx::slan_600][
        t%buffer[Par_idx::slan_600].size()] > 0.5);}
#endif

#endif
    uint32_t get_traj_idx() const {return traj_idx;}

 private:
    std::vector<int> startp, countp;
    int ncid;
    uint64_t n_timesteps_in; /*!< Total number of time steps that can be read from the input file. */
    uint32_t n_timesteps_buffer; /*!< Number of time steps that can be stored in the buffer. */
    uint64_t time_buffer_idx; /*!< Current index to read from the buffer. */

    uint32_t traj_idx; /*!< Index of trajectory to read from. */
    uint32_t ens_idx; /*!< Index of ensemble to read from. */
    bool already_open; /*!< Is the netCDF file already open? */
    uint32_t n_subs; /*!< Time step size of model from file / size of this model . */
    double pascal_conv; /*!< Convert pressure to Pa even if it is stored as hPa. */

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
    #if defined(B_EIGHT)
        qh_in,
    #endif
        ni_in,
        ns_in,
        ng_in,
    #if defined(B_EIGHT)
        nh_in,
    #endif
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
#if !defined(B_EIGHT)
        ensemble_dim_idx,
#endif
        n_dims
    };
    // Index for variable ids that are read during initialization
    enum Par_once_idx {
        qv,
        qc,
        qr,
#if defined(RK4ICE)
        qi,
        qs,
        qg,
    #if defined(B_EIGHT)
        qh,
    #endif
        nc,
        nr,
        ni,
        ns,
        ng,
    #if defined(B_EIGHT)
        nh,
    #endif
#endif
        time,
        n_pars_once
    };

    std::array<std::vector<double>, Par_idx::n_pars > buffer;
#if defined B_EIGHT
    std::unordered_map<std::string, std::string> const reader_names = {
        {"QV", "QV"}, {"QC", "QC"}, {"QR", "QR"},
        {"QI", "QI"}, {"QG", "QG"}, {"QH", "QH"},
        {"QS", "QS"}, {"NCGRAUPEL", "QNG"}, {"NCICE", "QNI"},
        {"NCSNOW", "QNS"}, {"NCHAIL", "QNH"}, {"NCCLOUD", "QNC"},
        {"NCRAIN", "QNR"}, {"time", "time"}, {"lat", "lat"}, {"lon", "lon"},
        {"time_rel", "time_after_asc_start"}, {"pressure", "p"}, {"T", "T"},
        {"QI_IN", "QI_in"}, {"QS_IN", "QS_in"}, {"QR_IN", "QR_in"},
        {"QG_IN", "QG_in"}, {"NI_IN", "QNI_in"}, {"NS_IN", "QNS_in"},
        {"NR_IN", "QNR_in"}, {"NG_IN", "QNG_in"}, {"z", "z"}, {"w", "w"},
        {"QH_IN", "QH_in"}, {"NH_IN", "NH_in"},
        {"conv_400", "conv_400"}, {"conv_600", "conv_600"}, {"slan_400", "slan_400"}, {"slan_600", "slan_600"},
        {"Q_TURBULENCE", "Q_TURBULENCE"}, {"trajectory", "trajectory"}, {"ensemble", "ensemble"}
    };
#elif defined WCB
    std::unordered_map<std::string, std::string> const reader_names = {
        {"QV", "QV"}, {"QC", "QC"}, {"QR", "QR"},
        {"QI", "QI"}, {"QG", "QG"}, {"QH", "QH"},
        {"QS", "QS"}, {"NCGRAUPEL", "NCGRAUPEL"}, {"NCICE", "NCICE"},
        {"NCSNOW", "NCSNOW"}, {"NCHAIL", "NCHAIL"}, {"NCCLOUD", "NCCLOUD"}, {"NCRAIN", "NCRAIN"},
        {"time", "ntim"}, {"lat", "lat"}, {"lon", "lon"},
        {"time_rel", "time_rel"}, {"pressure", "P"}, {"T", "T"},
        {"QI_IN", "QI_IN"}, {"QS_IN", "QS_IN"}, {"QR_IN", "QR_IN"},
        {"QG_IN", "QG_IN"}, {"NI_IN", "NI_IN"}, {"NS_IN", "NS_IN"},
        {"NR_IN", "NR_IN"}, {"NG_IN", "NG_IN"}, {"z", "z"}, {"w", "w"},
        {"QH_IN", "QH_IN"}, {"NH_IN", "NH_IN"},
        {"conv_400", "conv_400"}, {"conv_600", "conv_600"}, {"slan_400", "slan_400"}, {"slan_600", "slan_600"},
        {"Q_TURBULENCE", "Q_TURBULENCE"}, {"trajectory", "ntra"}, {"ensemble", "ensemble"}
#elif defined WCB2
    std::unordered_map<std::string, std::string> const reader_names = {
        {"QV", "QV"}, {"QC", "QC"}, {"QR", "QR"},
        {"QI", "QI"}, {"QG", "QG"}, {"QH", "QH"},
        {"QS", "QS"}, {"NCGRAUPEL", "NCGRAUPEL"}, {"NCICE", "NCICE"},
        {"NCSNOW", "NCSNOW"}, {"NCCLOUD", "NCCLOUD"}, {"NCRAIN", "NCRAIN"},
        {"time", "time"}, {"lat", "latitude"}, {"lon", "longitude"},
        {"time_rel", "time_rel"}, {"pressure", "P"}, {"T", "T"},
        {"QI_IN", "QI_IN"}, {"QS_IN", "QS_IN"}, {"QR_IN", "QR_IN"},
        {"QG_IN", "QG_IN"}, {"NI_IN", "NI_IN"}, {"NS_IN", "NS_IN"},
        {"NR_IN", "NR_IN"}, {"NG_IN", "NG_IN"}, {"z", "z"}, {"w", "w"},
        {"QH_IN", "QH_IN"}, {"NH_IN", "NH_IN"},
        {"conv_400", "conv_400"}, {"conv_600", "conv_600"}, {"slan_400", "slan_400"}, {"slan_600", "slan_600"},
        {"Q_TURBULENCE", "Q_TURBULENCE"}, {"trajectory", "id"}, {"ensemble", "ensemble"}
    };
#elif defined MET3D
    std::unordered_map<std::string, const char*> const reader_names = {
        {"QV", "QV"}, {"QC", "QC"}, {"QR", "QR"},
        {"QI", "QI"}, {"QG", "QG"}, {"QH", "QH"},
        {"QS", "QS"}, {"NCGRAUPEL", "NCGRAUPEL"}, {"NCICE", "NCICE"},
        {"NCSNOW", "NCSNOW"}, {"NCCLOUD", "NCCLOUD"},
        {"NCRAIN", "NCRAIN"}, {"NCHAIL", "NCHAIL"},
        {"time", "time"}, {"lat", "lat"},
        {"lon", "lon"}, {"time_rel", "time_after_ascent"},
        {"pressure", "pressure"}, {"T", "T"},
        {"QI_IN", "QI_IN"}, {"QS_IN", "QS_IN"}, {"QR_IN", "QR_IN"},
        {"QG_IN", "QG_IN"}, {"NI_IN", "NI_IN"}, {"NS_IN", "NS_IN"},
        {"NR_IN", "NR_IN"}, {"NG_IN", "NG_IN"}, {"z", "z"}, {"w", "w"},
        {"QH_IN", "QH_IN"}, {"NH_IN", "NH_IN"},
        {"conv_400", "conv_400"}, {"conv_600", "conv_600"},
        {"slan_400", "slan_400"}, {"slan_600", "slan_600"},
        {"Q_TURBULENCE", "Q_TURBULENCE"}, {"trajectory", "trajectory"}, {"ensemble", "ensemble"}
    };
#else
    std::unordered_map<std::string, std::string> const reader_names = {
        {"QV", "qv"}, {"QC", "qc"}, {"QR", "qr"},
        {"qi", "qi"}, {"qg", "qg"}, {"qs", "qs"},
        {"qh", "qh"}, {"time", "time"}, {"lat", "lat"}, {"lon", "lon"}, {"time_rel", "time_rel"}, {"pressure", "p"},
        {"T", "t"}, {"z", "z"}, {"w", "w"},
        {"conv_400", "conv_400"}, {"conv_600", "conv_600"}, {"slan_400", "slan_400"}, {"slan_600", "slan_600"},
        {"Q_TURBULENCE", "Q_TURBULENCE"}, {"trajectory", "id"}, {"ensemble", "ensemble"}
    };
#endif

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
