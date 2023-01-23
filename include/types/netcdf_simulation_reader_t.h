#ifndef TRAJECTORIES_NETCDF_SIMULATION_READER_T_H
#define TRAJECTORIES_NETCDF_SIMULATION_READER_T_H

#include <cmath>
#include <array>
#include <netcdf.h>
#include <string>
#include <vector>

#include "include/misc/error.h"
#include "include/misc/constants_regrid.h"

struct netcdf_simulation_reader_t {
    uint64_t n_trajectories;
    uint64_t start_time_idx; /*!< The start time index. */
    uint64_t time_idx; /*!< Current index to read from netcdf file. */
    uint64_t n_readable_timesteps; /*!< Amount of timtesteps that have been loaded. */

    /**
     *
     * @param buffer_size Number of time steps to read at once
     */
    explicit netcdf_simulation_reader_t(const uint32_t &buffer_size);
    void init_netcdf(double &start_time);
    void set_ensemble_idx(uint32_t idx) {this->ens_idx = idx;}
    void set_traj_idx(uint32_t idx) {this->traj_idx = idx;}
    uint64_t get_buffer_size() {return n_timesteps_buffer;}

    double get_relative_time(const uint32_t &t) const {
        return buffer[Par_idx::time_after_ascent][t%buffer[Par_idx::time_after_ascent].size()];
    }
    bool get_asc600(const uint32_t &t) const {return (buffer[Par_idx::asc600][
                                                              t%buffer[Par_idx::asc600].size()] > 0.5);}
    uint32_t get_traj_idx() const {return traj_idx;}

    void load_vars(const char* input_file);
    void init_netcdf(
        const double &start_time,
        const uint32_t &traj_idx);

    std::array<std::vector<double>, Par_idx::n_pars > buffer;
    std::array<std::vector<double>, Sens_par_idx::n_sens_pars*3 > buffer_sens;
  private:
    std::vector<int> startp, countp;
    int ncid;
    uint64_t n_timesteps_in; /*!< Total number of time steps that can be read from the input file. */
    uint32_t n_timesteps_buffer; /*!< Number of time steps that can be stored in the buffer. */
    uint64_t time_buffer_idx; /*!< Current index to read from the buffer. */

    uint32_t traj_idx; /*!< Index of trajectory to read from. */
    uint32_t ens_idx; /*!< Index of ensemble to read from. */
    bool already_open; /*!< Is the netCDF file already open? */
    double pascal_conv; /*!< Convert pressure to Pa even if it is stored as hPa. */

    /**
     * ID for dimensions of output file.
     */
    std::vector<int> dimid;
    /**
     * ID for variables of output file.
     */
    std::vector<int> varid, sensvarid;





};

#endif //TRAJECTORIES_NETCDF_SIMULATION_READER_T_H
