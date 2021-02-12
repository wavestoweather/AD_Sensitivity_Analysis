#pragma once

#include <cmath>
#include "ncType.h"
#include <netcdf>
#include <cstring>
#include <vector>

#include "include/types/reference_quantities_t.h"

using namespace netCDF;

/**
 * Structure to collect all nc parameters.
 */
struct nc_parameters_t{

    uint32_t n_trajectories = 30; /*!< Number of trajectories in the netCDF file. */
    uint32_t n_timesteps = 7922; /*!< Number of timesteps in the netCDF file. */
#ifdef MET3D
    double z;
    std::vector<double> time_abs;
    uint64_t time_idx = 0;
#else
    std::vector<double> z;
#endif
#if defined MET3D && defined TURBULENCE
    NcVar qturb_var;
    double qturb;
#endif
    std::vector<double> w, lat, lon;
    double  t, p, time_rel,
            qc, qr, qi, qs, qg, qv, S, dw, dlat, dlon,
            QIin, QSin, QRin, QGin, QIout, QSout, QRout, QGout,
            NIin, NSin, NRin, NGin, NIout, NSout, NRout, NGout,
            Nc, Nr, Ni, Ns, Ng;
    bool ascent_flag, conv_400, conv_600, slan_400, slan_600, dp2h;
    char* type[1];
    NcVar   lat_var, lon_var, z_var, t_var, p_var, w_var, time_rel_var,
            qc_var, qr_var, qi_var, qs_var, qg_var, qv_var, S_var,
            QIin_var, QSin_var, QRin_var, QGin_var, QIout_var, QSout_var,
            QRout_var, QGout_var, ascent_flag_var,
            NIin_var, NSin_var, NRin_var, NGin_var, NIout_var, NSout_var,
            NRout_var, NGout_var, dp2h_var,
            Nc_var, Nr_var, Ni_var, Ns_var, Ng_var,
#ifdef MET3D
            type_var, time_abs_var,
#endif
            conv_400_var, conv_600_var, slan_400_var, slan_600_var;

    /**
     * Load variables from the netCDF file such that data can be loaded using
     * load_params()
     *
     * @param nc Struct where to load the variables.
     * @param datafile The netCDF file.
     */
    void load_vars(NcFile &datafile);
    /**
     * Load the parameters from the netCDF file where load_vars()
     * must have been called beforehand.
     *
     * @param nc s Struct where to store the values.
     * @param startp Must have two values with index from where to load values.
     *               Depending on the NetCDF file, startp[0] may refer to the
     *               trajectory id.
     * @param countp Must have two values for how many values to load.
     *               Depending on the NetCDF file, countp[0] may refer to the
     *               trajectory id. Usually we set the values to 1.
     * @param ref_quant Reference quantities to transform between units.
     */
    void load_params(std::vector<size_t> &startp,
        std::vector<size_t> &countp,
        const reference_quantities_t &ref_quant,
        uint64_t num_sub_steps);

    /**
     * Initialize the nc parameters to default values. and allocate memory.
     *
     * @param nc The struct where the parameters are initialized.
     * @param n The number of trajectories in the netCDF file.
     * @param n_timesteps The maximum number of timesteps in the netCDF file.
     */
    void init_params(
        uint32_t n=32,
        uint32_t n_time=7922);
};
