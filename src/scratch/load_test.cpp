#include "codi.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <tgmath.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>
#include "include/microphysics/physical_parameterizations.h"

#include "include/microphysics/program_io.h"
#include "include/microphysics/constants.h"
#include "include/microphysics/general.h"

#include <netcdf>

int main(int argc, char** argv)
{
    nc_parameters_t nc_params;
    reference_quantities_t ref_quant;

    ref_quant.Tref = 273.15;
    ref_quant.pref = 1.0e5; // 1.0e3
    ref_quant.qref = 1.0e-4; // 1.0e-4
    ref_quant.wref = 1.; // 10.0
    ref_quant.tref = 1.0;

    ref_quant.Nref = 1.0; 	// DUMMY
    const uint32_t traj = 0;
    // char* file = "/mnt/localscratch/data/project/m2_jgu-tapt/online_trajectories/wcb201609_vladiana/O_WCB_all_20160922_00.nc";
#ifdef MET3D
    const char* file = "/data/project/wcb/netcdf/vladiana_met/conv_400_0_traj_t000000_p001.nc_wcb";
#else
    const char* file = "/data/project/m2_jgu-tapt/online_trajectories/traj_new_appended/traj_t000120_p001.nc_wcb";
#endif
    int traj_id;

    try
    {
        int dimid, ncid;
        size_t lenp, n_timesteps;
        // Get the amount of trajectories
        nc_open(file, NC_NOWRITE, &ncid);
#ifdef WCB
        nc_inq_dimid(ncid, "ntra", &dimid);
#elif defined MET3D
        nc_inq_dimid(ncid, "trajectory", &dimid);
#else
        nc_inq_dimid(ncid, "id", &dimid);
#endif
        nc_inq_dimlen(ncid, dimid, &lenp);
        // Get the amount of timesteps
#ifdef WCB
        nc_inq_dimid(ncid, "ntim", &dimid);
#else
        nc_inq_dimid(ncid, "time", &dimid);
#endif
        nc_inq_dimlen(ncid, dimid, &n_timesteps);
        std::cout << "Number of trajectories in netCDF file: " << lenp << "\n";
        std::cout << "Number of timesteps: " << n_timesteps << "\n";
        if(lenp <= 0)
        {
            std::cout << "You asked for trajectory with index " << 0
                      << " which does not exist. ABORTING.\n";
            return 1;
        }

        init_nc_parameters(nc_params, lenp, n_timesteps);

        netCDF::NcFile datafile(file, netCDF::NcFile::read);
        load_nc_parameters_var(nc_params, datafile);
        netCDF::NcVar id_var;
#ifdef MET3D
        id_var = datafile.getVar("trajectory");
#else
        id_var = datafile.getVar("id");
#endif
        std::vector<int> ids(lenp);
        id_var.getVar(ids.data());
        traj_id = ids[0];

        std::vector<size_t> startp, countp;
#ifdef MET3D
        startp.push_back(0);
        countp.push_back(1);
#endif
        startp.push_back(0); // start row, trajectory id or for wcb time
        startp.push_back(0); // start column, time or for wcb trajectory
        countp.push_back(1);
        countp.push_back(1);
        load_nc_parameters(nc_params, startp, countp, ref_quant, 1);
        std::cout << std::setprecision(8) << std::fixed
                  << std::setfill('0') << std::setw(6);
        std::cout << "Trajectory 0 at t=0\n";
        std::cout << "p\t\tT\t\tw\t\tw2\t\tqc\t\tqr\t\tqv\t\tqi\t\tqs\t\tqg\t\tlat\t\tlon\t\ttraj_id\n"
            << nc_params.p << "\t" << nc_params.t << "\t" << nc_params.w[0] << "\t"
            << nc_params.w[1] << "\t" << nc_params.qc << "\t" << nc_params.qr << "\t"
            << nc_params.qv << "\t" << nc_params.qi << "\t" << nc_params.qs << "\t"
            << nc_params.qg << "\t" << nc_params.lat[0] << "\t" << nc_params.lon[0] << "\t"
            << traj_id << "\n";

        load_nc_parameters_var(nc_params, datafile);
        load_nc_parameters(nc_params, startp, countp, ref_quant, 1);
        std::cout << "Trajectory 0 at t=0\n";
        std::cout << "p\t\tT\t\tw\t\tw2\t\tqc\t\tqr\t\tqv\t\tqi\t\tqs\t\tqg\t\tlat\t\tlon\t\ttraj_id\n"
            << nc_params.p << "\t" << nc_params.t << "\t" << nc_params.w[0] << "\t"
            << nc_params.w[1] << "\t" << nc_params.qc << "\t" << nc_params.qr << "\t"
            << nc_params.qv << "\t" << nc_params.qi << "\t" << nc_params.qs << "\t"
            << nc_params.qg << "\t" << nc_params.lat[0] << "\t" << nc_params.lon[0] << "\t"
            << traj_id << "\n";
#ifdef MET3D
        startp[2] += 1;
#else
        startp[0] += 1;
#endif
        load_nc_parameters_var(nc_params, datafile);
        load_nc_parameters(nc_params, startp, countp, ref_quant, 1);
        std::cout << "Trajectory 0 at t=1\n";
        std::cout << "p\t\tT\t\tw\t\tw2\t\tqc\t\tqr\t\tqv\t\tqi\t\tqs\t\tqg\t\tlat\t\tlon\t\ttraj_id\n"
            << nc_params.p << "\t" << nc_params.t << "\t" << nc_params.w[0] << "\t"
            << nc_params.w[1] << "\t" << nc_params.qc << "\t" << nc_params.qr << "\t"
            << nc_params.qv << "\t" << nc_params.qi << "\t" << nc_params.qs << "\t"
            << nc_params.qg << "\t" << nc_params.lat[0] << "\t" << nc_params.lon[0] << "\t"
            << traj_id << "\n";
        traj_id = ids[1];
#ifdef MET3D
        startp[1] = 1;
        startp[2] = 1;
#else
        startp[0] = 1;
        startp[1] = 1;
#endif
        load_nc_parameters_var(nc_params, datafile);
        load_nc_parameters(nc_params, startp, countp, ref_quant, 1);
        std::cout << "Trajectory 1 at t=1\n";
        std::cout << "p\t\tT\t\tw\t\tw2\t\tqc\t\tqr\t\tqv\t\tqi\t\tqs\t\tqg\t\tlat\t\tlon\t\ttraj_id\n"
            << nc_params.p << "\t" << nc_params.t << "\t" << nc_params.w[0] << "\t"
            << nc_params.w[1] << "\t" << nc_params.qc << "\t" << nc_params.qr << "\t"
            << nc_params.qv << "\t" << nc_params.qi << "\t" << nc_params.qs << "\t"
            << nc_params.qg << "\t" << nc_params.lat[0] << "\t" << nc_params.lon[0] << "\t"
            << traj_id << "\n";
#ifdef MET3D
        startp[2] = 2;
#else
        startp[1] = 2;
#endif
        load_nc_parameters_var(nc_params, datafile);
        load_nc_parameters(nc_params, startp, countp, ref_quant, 1);
        std::cout << "Trajectory 1 at t=2\n";
        std::cout << "p\t\tT\t\tw\t\tw2\t\tqc\t\tqr\t\tqv\t\tqi\t\tqs\t\tqg\t\tlat\t\tlon\t\ttraj_id\n"
            << nc_params.p << "\t" << nc_params.t << "\t" << nc_params.w[0] << "\t"
            << nc_params.w[1] << "\t" << nc_params.qc << "\t" << nc_params.qr << "\t"
            << nc_params.qv << "\t" << nc_params.qi << "\t" << nc_params.qs << "\t"
            << nc_params.qg << "\t" << nc_params.lat[0] << "\t" << nc_params.lon[0] << "\t"
            << traj_id << "\n";
    } catch(netCDF::exceptions::NcException& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "ABORTING." << std::endl;
        return 1;
    }
}