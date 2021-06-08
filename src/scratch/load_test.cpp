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

void print_params(
    nc_parameters_t &nc_params,
    const int &traj_id)
{
#ifdef MET3D
    std::cout << "p\t\tT\t\tw\t\tw2\t\tqc\t\tqr\t\tqv\t\tqi\t\tqs\t\tqg\t\tlat\t\tlon"
              << "\t\t\ttraj_id\t\ttime\t\ttime_after_ascent\t\ttype\n"
#else
    std::cout << "p\t\tT\t\tw\t\tw2\t\tqc\t\tqr\t\tqv\t\tqi\t\tqs\t\tqg\t\tlat\t\tlon\t\ttraj_id\n"
#endif
        << nc_params.p << "\t" << nc_params.t << "\t" << nc_params.w[0] << "\t"
        << nc_params.w[1] << "\t" << nc_params.qc << "\t" << nc_params.qr << "\t"
        << nc_params.qv << "\t" << nc_params.qi << "\t" << nc_params.qs << "\t"
        << nc_params.qg << "\t" << nc_params.lat[0] << "\t" << nc_params.lon[0] << "\t\t"
        << traj_id
#ifdef MET3D
        << "\t\t" << nc_params.time_abs[0] << "\t"
        << nc_params.time_rel << "\t\t" << nc_params.type[0] << "\n";
#else
        << "\n";
#endif
}


int main(int argc, char** argv)
{
    std::cout << "~+~+~+~+~Starting Load Tests~+~+~+~+~\n";
    nc_parameters_t nc_params;
    reference_quantities_t ref_quant;

    ref_quant.Tref = 273.15;
    ref_quant.pref = 1.0e5; // 1.0e3
    ref_quant.qref = 1.0e-4; // 1.0e-4
    ref_quant.wref = 1.; // 10.0
    ref_quant.tref = 1.0;

    ref_quant.Nref = 1.0; 	// DUMMY
    const uint32_t traj = 0;
    const char* file = argv[1];
    // Some datafiles I used for testing
    // const char* file = "/data/project/wcb/netcdf/vladiana_met/conv_400_0_traj_t000000_p001.nc_wcb";
    // const char* file = "/data/project/wcb/netcdf/vladiana_met_stats/no_exclusions_conv_400_median.nc_wcb";
    const char* file_config = argv[2];
    model_constants_t cc;
    std::vector<segment_t> segments;
    SUCCESS_OR_DIE(load_ens_config(file_config, cc, segments));
    for(auto &s: segments)
        SUCCESS_OR_DIE(s.check());
    print_segments(segments);

    std::vector<double> y(num_comp);
    // Populate y with some numbers
    for(uint32_t i=0; i<num_comp; ++i)
        y[i] = file[i%std::strlen(file)];

    input_parameters_t input;
    double current_time = 0;
    // Write a checkpoint file
    write_checkpoint("tmp_checkpoint", cc, y, segments, input, current_time);

    // and load it again
    SUCCESS_OR_DIE(load_checkpoint("tmp_checkpoint_0.json",
        cc, y, segments, input, ref_quant));
    for(auto &s: segments)
        SUCCESS_OR_DIE(s.check());

    // Perturb parameters for one segment
    segments[0].perturb(cc);
    // Loading it later and checking it would fail, if the following line
    // were missing. This is only needed for the test case here.
    // During operation, we would load the checkpoint first and then
    // perturb.
    segments[0].n_segments++;

    // Write again as checkpoint file
    write_checkpoint("tmp_checkpoint_perturbed", cc, y, segments, input, current_time);

    // Load it again
    SUCCESS_OR_DIE(load_checkpoint("tmp_checkpoint_perturbed_0.json",
        cc, y, segments, input, ref_quant));

    for(auto &s: segments)
        SUCCESS_OR_DIE(s.check());

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
#ifdef MET3D
        // Read global attributes
        std::cout << "Global attributes:\n";
        auto attributes = datafile.getAtts();
        for(auto & name_attr: attributes)
        {
            auto attribute = name_attr.second;
            netCDF::NcType type = attribute.getType();
            if(type.getName() == "double")
            {
                std::vector<float> values(1);
                attribute.getValues(values.data());
                std::cout << attribute.getName() << "\n\t"
                            << "type: " << type.getName() << "\n\t"
                            << "values: " << values[0] << "\n";
            } else if(type.getName() == "int64" || type.getName() == "int32" || type.getName() == "int")
            {
                std::vector<int> values(1);
                attribute.getValues(values.data());
                std::cout << attribute.getName() << "\n\t"
                          << "type: " << type.getName() << "\n\t"
                          << "values: " << values[0] << "\n";
            }
        }

        // Read attributes from each column
        std::cout << "Non global attributes:\n";
        auto vars = datafile.getVars();
        for(auto & name_var: vars)
        {
            auto var = name_var.second;
            auto name = name_var.first;
            std::cout << name << "\n";
            auto attributes = var.getAtts();
            for(auto & name_attr: attributes)
            {
                std::cout << "\t";
                auto attribute = name_attr.second;
                netCDF::NcType type = attribute.getType();
                if(type.getName() == "double" || type.getName() == "float")
                {
                    std::vector<double> values(1);
                    attribute.getValues(values.data());
                    std::cout << attribute.getName() << "=" << values[0] << "\n";
                } else if(type.getName() == "int64" || type.getName() == "int32" || type.getName() == "int")
                {
                    std::vector<int> values(1);
                    attribute.getValues(values.data());
                    std::cout << attribute.getName() << "=" << values[0] << "\n";
                } else if(type.getName() == "char")
                {
                    std::string values;
                    attribute.getValues(values);
                    std::cout << attribute.getName() << "=" << values << "\n";
                }
            }
        }


#endif
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
                  << std::setfill('0') << std::setw(8);
        std::cout << "Trajectory 0 at t=0\n";

        print_params(nc_params, traj_id);

        load_nc_parameters_var(nc_params, datafile);
        load_nc_parameters(nc_params, startp, countp, ref_quant, 1);
        std::cout << "Trajectory 0 at t=0\n";
        print_params(nc_params, traj_id);

#ifdef MET3D
        startp[2] += 1;
#else
        startp[0] += 1;
#endif
        load_nc_parameters_var(nc_params, datafile);
        load_nc_parameters(nc_params, startp, countp, ref_quant, 1);
        std::cout << "Trajectory 0 at t=1\n";
        print_params(nc_params, traj_id);
        if(lenp > 1)
        {
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
            print_params(nc_params, traj_id);
#ifdef MET3D
            startp[2] = 2;
#else
            startp[1] = 2;
#endif
            load_nc_parameters_var(nc_params, datafile);
            load_nc_parameters(nc_params, startp, countp, ref_quant, 1);
            std::cout << "Trajectory 1 at t=2\n";
            print_params(nc_params, traj_id);
        }
    } catch(netCDF::exceptions::NcException& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "ABORTING." << std::endl;
        return 1;
    }
}