#ifndef PROGRAM_IO_H
#define PROGRAM_IO_H

#include <stdlib.h>
#include <cmath>
#include <string>
#include <netcdf>
#include <vector>
#include "constants.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include "codi.hpp"
#include "types.h"
using namespace netCDF;


/** @defgroup io Input/Output Functions
 * Functions to load data and setup and initialize structs.
 * @{
 */

/**
 * Based on
 * init_dmin_wg_gr_ltab_equi('dmin_wetgrowth_lookup.dat', unitnr, 61, ltabdminwgg)
 * Loads a lookup table for calculating the growth of rain droplets or something.
 *
 * @param table Reference to the dataobject where the lookup table is going to
 *              be stored.
 * @param filename Path to the data file.
 * @param ndT Number of equidistant table vectors.
 */
bool load_lookup_table(
    table_t &table,
    std::string filename = "dmin_wetgrowth_lookup.dat",
    uint32_t ndT = 61)
{
    // T may be nonequidistant
    // p, qw and qi must be equidistant
    std::ifstream data;
    data.open(filename);
    // Read file
    if(data.is_open())
    {
        std::string line;

        std::getline(data, line);
        std::istringstream stream(line);
        std::vector<std::string> line_vec(
            std::istream_iterator<std::string>{stream},
            std::istream_iterator<std::string>());

        table.n1 = std::stoi(line_vec[0]);
        table.n3 = std::stoi(line_vec[2]);
        table.n4 = std::stoi(line_vec[3]);
        int anzT_wg_loc = std::stoi(line_vec[1]);

        std::vector<double> Tvec_wg_g_loc(anzT_wg_loc);
        std::vector<double> dmin_wg_g_loc(anzT_wg_loc*table.n1*table.n3*table.n4);
        table.x1.resize(table.n1);
        table.x3.resize(table.n3);
        table.x4.resize(table.n4);

        uint64_t counter = 0;
        // The first line holds x1 (p)
        std::getline(data, line);
        stream = std::istringstream(line);
        line_vec = std::vector<std::string> (
            std::istream_iterator<std::string>{stream},
            std::istream_iterator<std::string>());
        for(auto &val: line_vec)
        {
            table.x1[counter] = std::stod(val);
            counter++;
        }

        // The second line holds Tvec_wg_g_loc (T)
        counter = 0;
        std::getline(data, line);
        stream = std::istringstream(line);
        line_vec = std::vector<std::string> (
            std::istream_iterator<std::string>{stream},
            std::istream_iterator<std::string>());
        for(auto &val: line_vec)
        {
            Tvec_wg_g_loc[counter] = std::stod(val);
            counter++;
        }

        // The third line holds x3 (qw)
        counter = 0;
        std::getline(data, line);
        stream = std::istringstream(line);
        line_vec = std::vector<std::string> (
            std::istream_iterator<std::string>{stream},
            std::istream_iterator<std::string>());
        for(auto &val: line_vec)
        {
            table.x3[counter] = std::stod(val);
            counter++;
        }

        // The fourth line holds x4 (qi)
        counter = 0;
        std::getline(data, line);
        stream = std::istringstream(line);
        line_vec = std::vector<std::string> (
            std::istream_iterator<std::string>{stream},
            std::istream_iterator<std::string>());
        for(auto &val: line_vec)
        {
            table.x4[counter] = std::stod(val);
            counter++;
        }

        // Read the values, one value per line
        for(uint64_t n4=0; n4<table.n4; ++n4)
            for(uint64_t n3=0; n3<table.n3; ++n3)
                for(uint64_t i=0; i<anzT_wg_loc; ++i)
                    for(uint64_t n1=0; n1<table.n1; ++n1)
                    {
                        std::getline(data, line);
                        dmin_wg_g_loc[
                              n4*table.n1*anzT_wg_loc*table.n3
                            + n3*table.n1*anzT_wg_loc
                            + i*table.n1
                            + n1] = stod(line);
                    }
        data.close();

        // Generate equidistant table vectors by linear oversampling
        table.n2 = ndT;
        table.x2.resize(ndT);
        table.table.resize(table.n1*table.n2*table.n3*table.n4);
        double minT = Tvec_wg_g_loc[0];
        double maxT = Tvec_wg_g_loc[Tvec_wg_g_loc.size()-1];
        table.dx1 = table.x1[1] - table.x1[0];
        table.odx1 = 1.0/table.dx1;
        table.dx2 = (maxT-minT)/(ndT-1.0);
        table.odx2 = 1.0/table.dx2;
        table.dx3 = table.x3[1] - table.x3[0];
        table.odx3 = 1.0/table.dx3;
        table.dx4 = table.x4[1] - table.x4[0];
        table.odx4 = 1.0/table.dx4;

        for(uint64_t i=0; i<table.n2; ++i)
            table.x2[i] = minT + i*table.dx2;

        // Linear interpolation w.r.t. T
        for(uint64_t i=0; i<table.n2; ++i)
        {
            uint64_t iu = 0;
            for(uint64_t j=0; j<anzT_wg_loc-1; ++j)
            {
                if(table.x2[i] >= Tvec_wg_g_loc[j]
                    && table.x2[i] <= Tvec_wg_g_loc[j+1])
                {
                    iu = j;
                    break;
                }
            }
            uint64_t io = iu+1;

            // actual linear interpolation
            for(uint64_t n4=0; n4<table.n4; ++n4)
                for(uint64_t n3=0; n3<table.n3; ++n3)
                    for(uint64_t n1=0; n1<table.n1; ++n1)
                        table.table[
                              n4*table.n1*table.n2*table.n3
                            + n3*table.n1*table.n2
                            + i*table.n1
                            + n1]
                            =  dmin_wg_g_loc[
                                  n4*table.n1*anzT_wg_loc*table.n3
                                + n3*table.n1*anzT_wg_loc
                                + iu*table.n1
                                + n1]
                              + (dmin_wg_g_loc[
                                  n4*table.n1*anzT_wg_loc*table.n3
                                + n3*table.n1*anzT_wg_loc
                                + io*table.n1
                                + n1]
                              - dmin_wg_g_loc[
                                  n4*table.n1*anzT_wg_loc*table.n3
                                + n3*table.n1*anzT_wg_loc
                                + iu*table.n1
                                + n1]) / (Tvec_wg_g_loc[io]-Tvec_wg_g_loc[iu])
                                    * (table.x2[i] - Tvec_wg_g_loc[iu]);
        }

        return true;
    } else
    {
        std::cout << "Error loading " << filename << ". Does the file exist?\n";
        return false;
    }
}


/**
 * Initialize the nc parameters to default values. and allocate memory.
 *
 * @param nc The struct where the parameters are initialized.
 * @param n The number of trajectories in the netCDF file.
 * @param n_timesteps The maximum number of timesteps in the netCDF file.
 */
void init_nc_parameters(
    nc_parameters_t &nc,
    uint32_t n=32,
    uint32_t n_timesteps=7922)
{
    nc.n_trajectories = n;
    nc.n_timesteps = n_timesteps;
    nc.w.resize(2);
    nc.z.resize(4);
    nc.lat.resize(2);
    nc.lon.resize(2);
}


/**
 * Load variables from the netCDF file such that data can be loaded using
 * load_nc_parameters()
 *
 * @param nc Struct where to load the variables.
 * @param datafile The netCDF file.
 */
void load_nc_parameters_var(
    nc_parameters_t &nc,
    netCDF::NcFile &datafile)
{
#ifdef WCB2
    nc.lat_var      = datafile.getVar("latitude");
    nc.lon_var      = datafile.getVar("longitude");
#else
    nc.lat_var      = datafile.getVar("lat");
    nc.lon_var      = datafile.getVar("lon");
#endif
    nc.z_var        = datafile.getVar("z");

#if defined WCB || defined WCB2
    nc.p_var        = datafile.getVar("P");
    nc.t_var        = datafile.getVar("T");
    nc.qc_var       = datafile.getVar("QC");
    nc.qv_var       = datafile.getVar("QV");
    nc.qr_var       = datafile.getVar("QR");
    nc.qi_var       = datafile.getVar("QI");
    nc.qs_var       = datafile.getVar("QS");
    nc.time_rel_var = datafile.getVar("time");
#else
    nc.p_var        = datafile.getVar("p");
    nc.t_var        = datafile.getVar("t");
    nc.w_var        = datafile.getVar("w");
    nc.time_rel_var = datafile.getVar("time_rel");
    nc.qc_var       = datafile.getVar("qc");
    nc.qr_var       = datafile.getVar("qr");
    nc.qi_var       = datafile.getVar("qi");
    nc.qs_var       = datafile.getVar("qs");
    nc.qg_var       = datafile.getVar("qg");
    nc.qv_var       = datafile.getVar("qv");
    nc.QIin_var     = datafile.getVar("QIin");
    nc.QSin_var     = datafile.getVar("QSin");
    nc.QRin_var     = datafile.getVar("QRin");
    nc.QGin_var     = datafile.getVar("QGin");
    nc.QIout_var    = datafile.getVar("QIout");
    nc.QSout_var    = datafile.getVar("QSout");
    nc.QRout_var    = datafile.getVar("QRout");
    nc.QGout_var    = datafile.getVar("QGout");
#endif
#ifdef WCB
    // specific humidity
    nc.S_var        = datafile.getVar("RELHUM");
    // Flag wether an effective ascent region is reached
    nc.ascent_flag_var = datafile.getVar("MAP");
    // Potential vorticity (German: Wirbelstaerke)
    // nc.pot_vortic   = datafile.getVar("POT_VORTIC")
#endif
#ifdef WCB2
    nc.qg_var       = datafile.getVar("QG");
    nc.QIin_var     = datafile.getVar("QI_IN");
    nc.QSin_var     = datafile.getVar("QS_IN");
    nc.QRin_var     = datafile.getVar("QR_IN");
    nc.QGin_var     = datafile.getVar("QG_IN");
    nc.QIout_var    = datafile.getVar("QI_OUT");
    nc.QSout_var    = datafile.getVar("QS_OUT");
    nc.QRout_var    = datafile.getVar("QR_OUT");
    nc.QGout_var    = datafile.getVar("QG_OUT");

    nc.NIin_var     = datafile.getVar("NI_IN");
    nc.NSin_var     = datafile.getVar("NS_IN");
    nc.NRin_var     = datafile.getVar("NR_IN");
    nc.NGin_var     = datafile.getVar("NG_IN");
    nc.NIout_var    = datafile.getVar("NI_OUT");
    nc.NSout_var    = datafile.getVar("NS_OUT");
    nc.NRout_var    = datafile.getVar("NR_OUT");
    nc.NGout_var    = datafile.getVar("NG_OUT");

    nc.Nc_var       = datafile.getVar("NCCLOUD");
    nc.Nr_var       = datafile.getVar("NCRAIN");
    nc.Ni_var       = datafile.getVar("NCICE");
    nc.Ns_var       = datafile.getVar("NCSNOW");
    nc.Ng_var       = datafile.getVar("NCGRAUPEL");
    // Flag wether an effective ascent region is reached
    nc.ascent_flag_var = datafile.getVar("WCB_flag");
    // 2h ascent rate after Oertel et al. (2019)
    nc.dp2h_var     = datafile.getVar("dp2h");
    nc.conv_400_var = datafile.getVar("conv_400");
    nc.conv_600_var = datafile.getVar("conv_600");
    nc.slan_400_var = datafile.getVar("slan_400");
    nc.slan_600_var = datafile.getVar("slan_600");
#endif
}

/**
 * Load the parameters from the netCDF file where load_nc_parameters_var()
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
void load_nc_parameters(
    nc_parameters_t &nc,
    std::vector<size_t> &startp,
    std::vector<size_t> &countp,
    reference_quantities_t &ref_quant,
    uint64_t num_sub_steps)
{
    // startp[0] <- trajectory id
    // startp[1] <- timestep


#if defined WCB || defined WCB2
    countp[0]++;
    countp[0]++;
    nc.z_var.getVar(startp, countp, nc.z.data());
    nc.lat_var.getVar(startp, countp, nc.lat.data());
    nc.lon_var.getVar(startp, countp, nc.lon.data());
    countp[0]--;
    countp[0]--;
#else
    countp[1]++;
    countp[1]++;
    nc.z_var.getVar(startp, countp, nc.z.data());
    nc.lat_var.getVar(startp, countp, nc.lat.data());
    nc.lon_var.getVar(startp, countp, nc.lon.data());
    countp[1]--;
    countp[1]--;
#endif

#if defined WCB || defined WCB2
    int map = 0;
    nc.ascent_flag_var.getVar(startp, countp, &map);
    nc.ascent_flag = (map > 0) ? true : false;
#endif

    nc.t_var.getVar(startp, countp, &nc.t);
    nc.p_var.getVar(startp, countp, &nc.p);
    nc.time_rel_var.getVar(startp, countp, &nc.time_rel);
    nc.qc_var.getVar(startp, countp, &nc.qc);
    nc.qr_var.getVar(startp, countp, &nc.qr);
    nc.qi_var.getVar(startp, countp, &nc.qi);
    nc.qs_var.getVar(startp, countp, &nc.qs);
    nc.qv_var.getVar(startp, countp, &nc.qv);

    double psat_prime = saturation_pressure_water(nc.t);

#if !defined WCB2
    // We are reading in hPa. Convert to Pa
    nc.p        *= 100;
#endif
    nc.p        /= ref_quant.pref;
    nc.t        /= ref_quant.Tref;
    nc.qc       /= ref_quant.qref;
    nc.qr       /= ref_quant.qref;
    nc.qv       /= ref_quant.qref;
    nc.qi       /= ref_quant.qref;
    nc.qs       /= ref_quant.qref;

#if defined WCB2
    nc.NRin_var.getVar(startp, countp, &nc.NRin);
    nc.NIin_var.getVar(startp, countp, &nc.NIin);
    nc.NSin_var.getVar(startp, countp, &nc.NSin);
    nc.NGin_var.getVar(startp, countp, &nc.NGin);

    nc.NRout_var.getVar(startp, countp, &nc.NRout);
    nc.NIout_var.getVar(startp, countp, &nc.NIout);
    nc.NSout_var.getVar(startp, countp, &nc.NSout);
    nc.NGout_var.getVar(startp, countp, &nc.NGout);

    nc.Nc_var.getVar(startp, countp, &nc.Nc);
    nc.Nr_var.getVar(startp, countp, &nc.Nr);
    nc.Ni_var.getVar(startp, countp, &nc.Ni);
    nc.Ns_var.getVar(startp, countp, &nc.Ns);
    nc.Ng_var.getVar(startp, countp, &nc.Ng);

    nc.NRin     /= ref_quant.Nref;
    nc.NIin     /= ref_quant.Nref;
    nc.NSin     /= ref_quant.Nref;
    nc.NGin     /= ref_quant.Nref;

    nc.NRout    /= ref_quant.Nref;
    nc.NIout    /= ref_quant.Nref;
    nc.NSout    /= ref_quant.Nref;
    nc.NGout    /= ref_quant.Nref;

    nc.Nc       /= ref_quant.Nref;
    nc.Nr       /= ref_quant.Nref;
    nc.Ni       /= ref_quant.Nref;
    nc.Ns       /= ref_quant.Nref;
    nc.Ng       /= ref_quant.Nref;

    nc.NRin     = abs(nc.NRin);
    nc.NIin     = abs(nc.NIin);
    nc.NSin     = abs(nc.NSin);
    nc.NGin     = abs(nc.NGin);

    // Some additional flags for convective and
    // slantwise trajectories
    nc.conv_400_var.getVar(startp, countp, &map);
    nc.conv_400 = (map > 0) ? true : false;
    nc.conv_600_var.getVar(startp, countp, &map);
    nc.conv_600 = (map > 0) ? true : false;
    nc.slan_400_var.getVar(startp, countp, &map);
    nc.slan_400 = (map > 0) ? true : false;
    nc.slan_600_var.getVar(startp, countp, &map);
    nc.slan_600 = (map > 0) ? true : false;
    nc.dp2h_var.getVar(startp, countp, &map);
    nc.dp2h = (map > 0) ? true : false;
#endif

#if !defined WCB
    nc.S = 1.0;
    nc.qg_var.getVar(startp, countp, &nc.qg);
    nc.QIin_var.getVar(startp, countp, &nc.QIin);
    nc.QSin_var.getVar(startp, countp, &nc.QSin);
    nc.QRin_var.getVar(startp, countp, &nc.QRin);
    nc.QGin_var.getVar(startp, countp, &nc.QGin);
    nc.QIout_var.getVar(startp, countp, &nc.QIout);
    nc.QSout_var.getVar(startp, countp, &nc.QSout);
    nc.QRout_var.getVar(startp, countp, &nc.QRout);
    nc.QGout_var.getVar(startp, countp, &nc.QGout);

    nc.qg       /= ref_quant.qref;
    nc.QRin     /= ref_quant.qref;
    nc.QRout    /= ref_quant.qref;
    nc.QIin     /= ref_quant.qref;
    nc.QIout    /= ref_quant.qref;
    nc.QSin     /= ref_quant.qref;
    nc.QSout    /= ref_quant.qref;
    nc.QGin     /= ref_quant.qref;
    nc.QGout    /= ref_quant.qref;

    nc.QRin     = abs(nc.QRin);
    nc.QRout    = abs(nc.QRout);
    nc.QIin     = abs(nc.QIin);
    nc.QIout    = abs(nc.QIout);
    nc.QSin     = abs(nc.QSin);
    nc.QSout    = abs(nc.QSout);
    nc.QGin     = abs(nc.QGin);
    nc.QGout    = abs(nc.QGout);
#endif

#if !defined WCB && !defined WCB2
    nc.w_var.getVar(startp, countp, nc.w.data());
    countp[1]++;
    nc.w_var.getVar(startp, countp, nc.w.data());
    countp[1]--;

    nc.w[0]     /= ref_quant.wref;
    nc.w[1]     /= ref_quant.wref;
    nc.dw = nc.w[1] - nc.w[0];

    // I am not sure why, but those values can be negative...
    nc.QRin     = abs(nc.QRin);
    nc.QRout    = abs(nc.QRout);
    nc.QIin     = abs(nc.QIin);
    nc.QIout    = abs(nc.QIout);
    nc.QSin     = abs(nc.QSin);
    nc.QSout    = abs(nc.QSout);
    nc.QGin     = abs(nc.QGin);
    nc.QGout    = abs(nc.QGout);
#elif defined WCB
    nc.qc /= 1.0e6;
    nc.qr /= 1.0e6;
    nc.qv /= 1.0e6;
    nc.qi /= 1.0e6;
    nc.qs /= 1.0e6;
    // Calculate w by getting the z-coordinates
    // and divide it by the amount of substeps
    nc.w[0] = (nc.z[1] - nc.z[0]) / 20.0;
    nc.w[1] = (nc.z[2] - nc.z[1]) / 20.0;
    nc.S_var.getVar(startp, countp, &nc.S);
    nc.S /= 100; // to percentage
    // we do not change w
    nc.dw = 0;
#else
    // WCB 2
    // Calculate w by getting the z-coordinates
    // and divide it by the amount of substeps
    nc.w[0] = (nc.z[1] - nc.z[0]) / 20.0;
    nc.w[1] = (nc.z[2] - nc.z[1]) / 20.0;
    nc.dw = 0; // We do not interpolate w for now
#endif

    nc.dlat = (nc.lat[1] - nc.lat[0]) / num_sub_steps;
    nc.dlon = (nc.lon[1] - nc.lon[0]) / num_sub_steps;
}


/**
 * Initialize the input parameters to default values.
 *
 * @param in Structure where to store the input parameters.
 */
void init_input_parameters(input_parameters_t &in)
{
  // Numerics
  in.t_end_prime = 100.0;	// Seconds
  in.dt_prime = 0.01;		// Seconds
  in.snapshot_index = 200;
  in.dt_traject = 20;       // Seconds; fixed from paper
  // Filename for output
#if defined(RK4)
   in.OUTPUT_FILENAME = "data/rain_OUTPUT.txt";
#endif
#if defined(RK4NOICE)
    in.OUTPUT_FILENAME = "data/sb_OUTPUT.txt";
#endif
#if defined(RK4ICE)
    in.OUTPUT_FILENAME = "data/sb_ice_OUTPUT.txt";
#endif

  // Filename for input
  in.INPUT_FILENAME = "/mnt/localscratch/data/project/m2_jgu-tapt/online_trajectories/foehn201305_case/foehn201305_warming.nc";

  // Scaling factor
  in.scaling_fact = 1.0;	// No scaling
  in.start_over = true;
  in.fixed_iteration = false;
  in.auto_type = 3;
  in.traj = 0;
  in.write_index = 100000;
}

/**
 * String used to parse commandline input.
 */
static const char *optString = "w:f:d:i:b:o:l:s:t:a:r:?";


/**
 * Initialize global args for parsing command line arguments.
 *
 * @param arg Struct that shall be initialized.
 */
void init_global_args(global_args_t &arg)
{

  arg.final_time_flag = 0;
  arg.final_time_string = nullptr;

  arg.timestep_flag = 0;
  arg.timestep_string = nullptr;

  arg.snapshot_index_flag = 0;
  arg.snapshot_index_string = nullptr;

  arg.output_flag = 0;
  arg.output_string = nullptr;

  arg.scaling_fact_flag = 0;
  arg.scaling_fact_string = nullptr;

  arg.input_flag = 0;
  arg.input_file = nullptr;

  arg.start_over_flag = 0;
  arg.start_over_string = nullptr;

  arg.fixed_iteration_flag = 0;
  arg.fixed_iteration_string = nullptr;

  arg.auto_type_flag = 0;
  arg.auto_type_string = nullptr;

  arg.traj_flag = 0;
  arg.traj_string = nullptr;

  arg.write_flag = 0;
  arg.write_string = nullptr;
}

/**
 * Set the input parameters with the data from the global arguments.
 *
 * @param arg Stuct with command line arguments.
 * @param in Struct where the input parameters will be stored.
 */
void set_input_from_arguments(global_args_t &arg ,
			      input_parameters_t &in )
{
  // Final time
  if(1 == arg.final_time_flag){
    in.t_end_prime = std::strtod(arg.final_time_string, nullptr);
  }

  // Timestep
  if(1 == arg.timestep_flag){
    in.dt_prime = std::strtod(arg.timestep_string, nullptr);
  }

  // Snapshot index
  if(1 == arg.snapshot_index_flag){
    in.snapshot_index = std::stoi(arg.snapshot_index_string);
  }

  // Output
  if(1 == arg.output_flag){
    in.OUTPUT_FILENAME = arg.output_string;
  }

  // Input
  if(1 == arg.input_flag){
    in.INPUT_FILENAME = arg.input_file;
  }

  // Scaling factor
  if(1 == arg.scaling_fact_flag){
    in.scaling_fact = std::strtod(arg.scaling_fact_string, nullptr);
  }

  // Starting over
  if(1 == arg.start_over_flag){
    in.start_over = (strcmp(arg.start_over_string, "0"));
  }

  if(1 == arg.fixed_iteration_flag){
      in.fixed_iteration = (strcmp(arg.fixed_iteration_string, "0"));
  }

  // Auto type
  if(1 == arg.auto_type_flag){
    in.auto_type = std::stoi(arg.auto_type_string);
  }

  // Trajectory
  if(1 == arg.traj_flag){
    in.traj = std::stoi(arg.traj_string);
  }

  // Write index
  if(1 == arg.write_flag){
      in.write_index = std::stoi(arg.write_string);
  }
}


/**
 * Write the file with reference values ending with
 * "__reference_values.txt" that can be read in Python with
 * Numpy. The order is:
 * Temperature, pressure, mixing ratio, particle number, ascent velocity,
 * time, height.
 *
 * @param out_filename String with filename.
 * @param ref_quant reference_quantities_t with all the reference values.
 * @return Errorcode (0=no errors; 1=simulation breaking error)
 */
int write_reference_quantities(
    std::string &out_filename,
    reference_quantities_t &ref_quant)
{
    std::ofstream outfile_refs;
    outfile_refs.open(out_filename + "_reference_values.txt");
    outfile_refs.precision(10);

    if( !outfile_refs.is_open() )
    {
        std::cout << "ERROR while opening the outputfile. Aborting." << std::endl;
        return 1;
    }

    // Write the reference quantities
    // Write the reference quantities
    outfile_refs << ref_quant.Tref << " "
	       << ref_quant.pref << " "
	       << ref_quant.qref << " "
	       << ref_quant.Nref << " "
	       << ref_quant.wref << " "
	       << ref_quant.tref << " "
           << ref_quant.zref << "\n";

    outfile_refs.close();
    return 0;
}


/**
 * Write the header for simulation results files and for files with
 * derivatives which have "_diff_" in their name.
 *
 * @param out_filename String with filename.
 * @return Errorcode (0=no errors; 1=simulation breaking error)
 */
int write_headers(
    std::string &out_filename)
{
    std::string suffix = ".txt";
    std::string full_filename;
    full_filename = out_filename;
    full_filename += suffix;

    outfile.open(full_filename);
    outfile.precision(10);

    if( !outfile.is_open() )
    {
        std::cout << "ERROR while opening the outputfile. Aborting." << std::endl;
        return 1;
    }

    // Append the initial values and write headers
    out_tmp << "timestep,trajectory,LONGITUDE,LATITUDE,"
#if defined WCB
        << "MAP,";
#endif
#if defined WCB2
        << "MAP,"
        << "dp2h,"
        << "conv_400,"
        << "conv_600,"
        << "slan_400,"
        << "slan_600,";
#endif
        for(uint32_t i=0; i<output_par_idx.size(); ++i)
            out_tmp << output_par_idx[i]  <<
                ((i < output_par_idx.size()-1) ? "," : "\n");

    std::string basename = "_diff_";
    std::string fname;

    for(int ii = 0; ii < num_comp; ii++)
    {
        fname = out_filename;
        fname += basename;
        fname += std::to_string(ii);
        fname += suffix;

        out_diff[ii].open(fname);
        out_diff[ii].precision(10);
        if( !out_diff[ii].is_open() )
        {
            std::cout << "ERROR while opening outputfile. Aborting." << std::endl;
            return 1;
        }
        out_diff_tmp[ii]
            << "timestep,"
            << "trajectory,"
            << "Output Parameter,"
            << "LONGITUDE,"
            << "LATITUDE,"
#if defined WCB
            << "MAP,";
#endif
#if defined WCB2
            << "MAP,"
            << "dp2h,"
            << "conv_400,"
            << "conv_600,"
            << "slan_400,"
            << "slan_600,";
#endif
        for(uint32_t i=0; i<output_grad_idx.size(); ++i)
            out_diff_tmp[ii] << output_grad_idx[i]  <<
                ((i < output_grad_idx.size()-1) ? "," : "\n");
    } // End loop over all components

    return 0;
}


/**
 * Read initial values from the netcdf file and stores them to y_init.
 * Also stores the amount of trajectories in the input file to lenp and
 * several quantities to cc such as the number of steps to simulate and
 * to nc_params, which is used to read from netcdf files.
 *
 * @param y_init Array of num_comp many doubles
 * @param nc_params Struct used for reading netcdf files
 * @param lenp On out: Number of trajectories available
 * @param ref_quant Reference quantities used to change from netcdf units to
 *                  simulation units
 * @param input_file Path to input netcdf file as char array
 * @param traj ID of input trajectory to read
 * @param cc Model constants. On out: Added number of simulation steps
 * @return Errorcode (0=no errors; 1=simulation breaking error)
 */
int read_init_netcdf(
    std::vector<double> &y_init,
    nc_parameters_t &nc_params,
    size_t &lenp,
    reference_quantities_t &ref_quant,
    const char *input_file,
    const uint32_t traj,
    model_constants_t &cc)
{
    try
    {
        int dimid, ncid;
        size_t n_timesteps;
        // Get the amount of trajectories
        nc_open(input_file, NC_NOWRITE, &ncid);
#ifdef WCB
        nc_inq_dimid(ncid, "ntra", &dimid);
#else
        nc_inq_dimid(ncid, "id", &dimid);
#endif
        nc_inq_dimlen(ncid, dimid, &lenp);
        std::cout << "Number of trajectories in netCDF file: " << lenp << "\n";
        if(lenp <= traj)
        {
            std::cout << "You asked for trajectory with index " << traj
                      << " which does not exist. ABORTING.\n";
            return 1;
        }
        // Get the amount of timesteps
#ifdef WCB
        nc_inq_dimid(ncid, "ntim", &dimid);
#else
        nc_inq_dimid(ncid, "time", &dimid);
#endif
        nc_inq_dimlen(ncid, dimid, &n_timesteps);
        uint64_t n_timesteps_input = ceil(cc.t_end/20.0)-1;

        cc.num_steps = (n_timesteps-1 > n_timesteps_input) ? n_timesteps_input : n_timesteps-1;
        init_nc_parameters(nc_params, lenp, n_timesteps);
        netCDF::NcFile datafile(input_file, netCDF::NcFile::read);
        load_nc_parameters_var(nc_params, datafile);

        std::vector<size_t> startp, countp;
        // wcb files have a different ordering
#if defined WCB || defined WCB2
        startp.push_back(0); // time
        startp.push_back(traj); // trajectory id
#else
        startp.push_back(traj); // trajectory id
        startp.push_back(1); // time (where time == 0 only has zeros)
#endif
        countp.push_back(1);
        countp.push_back(1);
        load_nc_parameters(nc_params, startp, countp,
                           ref_quant, cc.num_sub_steps);

        y_init[p_idx] = nc_params.p;
        y_init[T_idx] = nc_params.t;

        y_init[S_idx] = nc_params.S;
#ifdef SAT_CALC
        y_init[S_idx] = nc_params.qv*ref_quant.qref * Rv * nc_params.t*ref_quant.Tref
            / saturation_pressure_water_icon(nc_params.t*ref_quant.Tref);
#endif
        y_init[qc_idx] = nc_params.qc;
        y_init[qr_idx] = nc_params.qr;
        y_init[qv_idx] = nc_params.qv;
#if defined(RK4ICE)
        y_init[qi_idx] = nc_params.qi;
        y_init[qs_idx] = nc_params.qs;
#endif
#ifdef WCB
        y_init[w_idx] = 0;
  #if defined(RK4ICE)
        y_init[qg_idx] = 0;
  #endif
#else
        y_init[w_idx] = nc_params.w[0];
  #if defined(RK4ICE)
        y_init[qg_idx] = nc_params.qg;
  #endif
#endif
#if defined(RK4ICE)
        y_init[qh_idx] = 0.0; // hail that is not in the trajectory
        y_init[qh_out_idx] = 0.0;
        y_init[Nh_out_idx] = 0.0;
#endif
#ifdef WCB2

        y_init[Nr_idx] = nc_params.Nr;
        y_init[Nc_idx] = nc_params.Nc;
  #if defined(RK4ICE)
        y_init[qr_out_idx] = nc_params.QRout;
        y_init[Nr_out_idx] = nc_params.NRout;
        y_init[qi_out_idx] = nc_params.QIout;
        y_init[qs_out_idx] = nc_params.QSout;
        y_init[qg_out_idx] = nc_params.QGout;
        y_init[Ni_out_idx] = nc_params.NIout;
        y_init[Ns_out_idx] = nc_params.NSout;
        y_init[Ng_out_idx] = nc_params.NGout;
        y_init[Ni_idx] = nc_params.Ni;
        y_init[Ns_idx] = nc_params.Ns;
        y_init[Ng_idx] = nc_params.Ng;
  #endif
#else
        y_init[qr_out_idx] = 0.0;
        y_init[Nr_out_idx] = 0;
        y_init[Nc_idx] = 0;
        y_init[Nr_idx] = 0;
  #if defined(RK4ICE)
        // We initialize the sedimentation with 0 for the stepper
        y_init[qi_out_idx] = 0.0;
        y_init[qs_out_idx] = 0.0;
        y_init[qg_out_idx] = 0.0;
        y_init[Ni_out_idx] = 0;
        y_init[Ns_out_idx] = 0;
        y_init[Ng_out_idx] = 0;
        y_init[Ni_idx] = 0;
        y_init[Ns_idx] = 0;
        y_init[Ng_idx] = 0;
  #endif
#endif
        y_init[Nv_idx] = 0;
#if defined(RK4ICE) || defined(RK4NOICE)
        y_init[z_idx] = nc_params.z[0];
        y_init[n_inact_idx] = 0;
        y_init[depo_idx] = 0;
        y_init[sub_idx] = 0;
#endif
    } catch(netCDF::exceptions::NcException& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "ABORTING." << std::endl;
        return 1;
    }
    return 0;
}



/**
 * Open NETCDF file for reading and store some information in ncid, startp,
 * countp.
 *
 * @param ncid On out: contains id of netcdf file (needed for closing it)
 * @param startp On out: contains dimensions info for reading
 * @params countp On out: contains dimensions info for reading
 * @params input_file Char array of input file
 */
void open_netcdf(
    int &ncid,
    std::vector<size_t> &startp,
    std::vector<size_t> &countp,
    const char *input_file,
    const uint32_t traj)
{
    nc_open(input_file, NC_NOWRITE, &ncid);
#if defined WCB || defined WCB2
    startp.push_back(1);          // time
    startp.push_back(traj); // trajectory
#else
    startp.push_back(traj); // trajectory
    startp.push_back(1);          // time
#endif
    countp.push_back(1);
    countp.push_back(1);
}


/**
 * Read from netcdf file. Alter current values used for simlation if
 * given timestep is encountered. Writes values to output stringstream.
 *
 *
 */
void read_netcdf_write_stream(
    const char *input_file,
    std::vector<size_t> &startp,
    std::vector<size_t> &countp,
    nc_parameters_t &nc_params,
    model_constants_t &cc,
    input_parameters_t &input,
    reference_quantities_t &ref_quant,
    std::vector<codi::RealReverse> &y_single_old,
    std::vector<codi::RealReverse> &inflow,
    std::vector<int> &ids,
    int &traj_id,
    const uint32_t t)
{
#if defined WCB || defined WCB2
    startp[0] = t;
#else
    startp[1] = t+1;
#endif

    netCDF::NcFile datafile(input_file, netCDF::NcFile::read);
    load_nc_parameters_var(nc_params, datafile);
    load_nc_parameters(nc_params, startp, countp,
                        ref_quant, cc.num_sub_steps);

    netCDF::NcVar id_var;
    id_var = datafile.getVar("id");
    id_var.getVar(ids.data());
    traj_id = ids[input.traj];
    // Set values from a given trajectory

    if(t==0 || !input.fixed_iteration)
    {
        y_single_old[p_idx]  = nc_params.p;     // p
        y_single_old[T_idx]  = nc_params.t;     // T
    }
    if(t==0 || input.start_over)
    {
        y_single_old[S_idx]  = nc_params.S;     // S
#ifdef SAT_CALC
        y_single_old[S_idx]  = nc_params.qv*ref_quant.qref * Rv * nc_params.t*ref_quant.Tref
            / saturation_pressure_water_icon(nc_params.t*ref_quant.Tref);
#endif
        y_single_old[qc_idx] = nc_params.qc;    // qc
        y_single_old[qr_idx] = nc_params.qr;    // qr
        y_single_old[qv_idx] = nc_params.qv;    // qv
#if defined(RK4ICE)
        y_single_old[qi_idx] = nc_params.qi;    // qi
        y_single_old[qs_idx] = nc_params.qs;    // qs
#endif
#if !defined(WCB) && defined(RK4ICE)
        y_single_old[qg_idx] = nc_params.qg;    // qg
#elif defined(RK4ICE)
        y_single_old[qg_idx] = 0;
#endif
#if defined(RK4ICE)
        y_single_old[qh_idx] = 0.0; // qh. We don't have hail in the trajectoris
        y_single_old[Nh_idx] = 0.0; // Nh. We don't have hail in the trajectoris
#endif
        codi::RealReverse denom = 0;
#if defined(RK4ICE) && defined(WCB2)
        y_single_old[Ng_idx] = nc_params.Ng;
        y_single_old[Ni_idx] = nc_params.Ni;
        y_single_old[Ns_idx] = nc_params.Ns;
        y_single_old[Ng_out_idx] = nc_params.NGout;
        y_single_old[Ni_out_idx] = nc_params.NIout;
        y_single_old[Ns_out_idx] = nc_params.NSout;
        y_single_old[Nr_out_idx] = nc_params.NRout;
#endif
#ifdef WCB2
        y_single_old[Nc_idx] = nc_params.Nc;
        y_single_old[Nr_idx] = nc_params.Nr;
#else
        denom = (cc.cloud.max_x - cc.cloud.min_x) / 2.0 + cc.cloud.min_x;
        y_single_old[Nc_idx] = y_single_old[qc_idx] * ref_quant.qref / (denom); //*10e2);  // Nc
        denom = (cc.rain.max_x - cc.rain.min_x) / 2 + cc.rain.min_x;
        y_single_old[Nr_idx] = y_single_old[qr_idx] * ref_quant.qref / (denom); //*10e2);  // Nr
        denom = cc.cloud.min_x / 2.0;
        y_single_old[Nv_idx] = y_single_old[qv_idx] * ref_quant.qref / (denom); //*10e2);  // Nv
#if defined(RK4ICE)
        denom = (cc.ice.max_x - cc.ice.min_x) / 2.0 + cc.ice.min_x;
        y_single_old[Ni_idx] = y_single_old[qi_idx] * ref_quant.qref / (denom); //*10e2); // Ni
        denom = (cc.snow.max_x - cc.snow.min_x) / 2.0 + cc.snow.min_x;
        y_single_old[Ns_idx] = y_single_old[qs_idx] * ref_quant.qref / (denom); //*10e2); // Ns
        denom = (cc.graupel.max_x - cc.graupel.min_x) / 2.0 + cc.graupel.min_x;
        y_single_old[Ng_idx] = y_single_old[qg_idx] * ref_quant.qref / (denom); //*10e2); // Ng
#endif
#endif
        cc.Nc_prime = y_single_old[Nc_idx];

        cc.rho_a_prime = compute_rhoa(nc_params.p*ref_quant.pref,//*100,
            nc_params.t*ref_quant.Tref, nc_params.S);
        y_single_old[w_idx]  = nc_params.w[0]; // w
        cc.dw = nc_params.dw / (cc.dt*cc.num_sub_steps);

        denom = cc.cloud.min_x / 2.0;
        y_single_old[Nv_idx] = y_single_old[qv_idx] * ref_quant.qref / (denom); //*10e2);  // Nv
#if defined(RK4ICE) || defined(RK4NOICE)
        y_single_old[z_idx] = nc_params.z[0];
#endif

#if defined(FLUX) && !defined(WCB)
        inflow[qr_in_idx] = nc_params.QRin;
  #if defined(RK4ICE)
        inflow[qi_in_idx] = nc_params.QIin;
        inflow[qs_in_idx] = nc_params.QSin;
        inflow[qg_in_idx] = nc_params.QGin;
  #endif
#else
        inflow[qr_in_idx] = 0;
  #if defined(RK4ICE)
        inflow[qi_in_idx] = 0;
        inflow[qs_in_idx] = 0;
        inflow[qg_in_idx] = 0;
  #endif
#endif
#if defined(FLUX) && defined(WCB2)
        inflow[Nr_in_idx] = nc_params.NRin;
  #if defined(RK4ICE)
        inflow[Ni_in_idx] = nc_params.NIin;
        inflow[Ns_in_idx] = nc_params.NSin;
        inflow[Ng_in_idx] = nc_params.NGin;
  #endif
#else
        inflow[Nr_in_idx] = 0;
  #if defined(RK4ICE)
        inflow[Ni_in_idx] = 0;
        inflow[Ns_in_idx] = 0;
        inflow[Ng_in_idx] = 0;
  #endif
#endif
    }

#if defined WCB || defined WCB2
    out_tmp << (t*cc.num_sub_steps)*cc.dt << "," << traj_id << ","
            << nc_params.lon[0] << "," << nc_params.lat[0] << ","
            << nc_params.ascent_flag << ",";
#else
    out_tmp << (t*cc.num_sub_steps)*cc.dt << "," << traj_id << ","
            << nc_params.lon[0] << "," << nc_params.lat[0] << ",";
#endif
#if defined WCB2
    out_tmp << nc_params.dp2h << "," << nc_params.conv_400 << ","
            << nc_params.conv_600 << "," << nc_params.slan_400 << ","
            << nc_params.slan_600 << ",";
#endif
    for(int ii = 0 ; ii < num_comp; ii++)
        out_tmp << y_single_old[ii] <<
            ((ii == num_comp-1) ? "\n" : ",");

    for(int ii = 0 ; ii < num_comp ; ii++)
    {
#if defined WCB || defined WCB2
        out_diff_tmp[ii] << t*cc.num_sub_steps*cc.dt << ","
                        << traj_id << ","
                        << output_par_idx[ii] << ","
                        << nc_params.lon[0] << ","
                        << nc_params.lat[0] << ","
                        << nc_params.ascent_flag << ",";
#else
        out_diff_tmp[ii] << t*cc.num_sub_steps*cc.dt << ","
                        << traj_id << ","
                        << output_par_idx[ii] << ","
                        << nc_params.lon[0] << ","
                        << nc_params.lat[0] << ",";
#endif
#if defined WCB2
        out_diff_tmp[ii] << nc_params.dp2h << "," << nc_params.conv_400 << ","
                            << nc_params.conv_600 << "," << nc_params.slan_400 << ","
                            << nc_params.slan_600 << ",";
#endif
        for(int jj = 0 ; jj < num_par ; jj++)
            out_diff_tmp[ii] << 0.0
                << ((jj==num_par-1) ? "\n" : ",");
    }
}


/**
 * Store simulation results and gradients to string stream and dump it
 * if timestep is reached to do so.
 */
void write_output(
    const model_constants_t &cc,
    const nc_parameters_t &nc_params,
    const std::vector<codi::RealReverse> &y_single_new,
    const std::vector< std::array<double, num_par > >  &y_diff,
    const uint32_t sub,
    const uint32_t t,
    const double time_new,
    const uint32_t traj_id,
    const uint32_t write_index,
    const uint32_t snapshot_index)
{
    if( (0 == (sub + t*cc.num_sub_steps) % snapshot_index)
        || ( t == cc.num_steps-1 && sub == cc.num_sub_steps-1 ) )
    {
        // Write the results to the output file
#if defined WCB || defined WCB2
        out_tmp << time_new << "," << traj_id << ","
                << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                << (nc_params.lat[0] + sub*nc_params.dlat) << ","
                << nc_params.ascent_flag << ",";
#else
        out_tmp << time_new << "," << traj_id << ","
                << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#endif
#if defined WCB2
        out_tmp << nc_params.dp2h << "," << nc_params.conv_400 << ","
                << nc_params.conv_600 << "," << nc_params.slan_400 << ","
                << nc_params.slan_600 << ",";
#endif
        for(int ii = 0 ; ii < num_comp; ii++)
            out_tmp << y_single_new[ii]
                << ((ii == num_comp-1) ? "\n" : ",");

        // CODIPACK: BEGIN
        for(int ii = 0 ; ii < num_comp ; ii++)
        {
#if defined WCB || defined WCB2
            out_diff_tmp[ii] << time_new << "," << traj_id << ","
                            << output_par_idx[ii] << ","
                            << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                            << (nc_params.lat[0] + sub*nc_params.dlat) << ","
                            << nc_params.ascent_flag << ",";
#else
            out_diff_tmp[ii] << time_new << "," << traj_id << ","
                            << output_par_idx[ii] << ","
                            << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                            << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#endif
#if defined WCB2
            out_diff_tmp[ii] << nc_params.dp2h << "," << nc_params.conv_400 << ","
                            << nc_params.conv_600 << "," << nc_params.slan_400 << ","
                            << nc_params.slan_600 << ",";
#endif
            for(int jj = 0 ; jj < num_par ; jj++)
                out_diff_tmp[ii] << y_diff[ii][jj]
                    << ((jj==num_par-1) ? "\n" : ",");
        }
        // CODIPACK: END
    }
    if( (0 == (sub + t*cc.num_sub_steps) % write_index)
        || ( t == cc.num_steps-1 && sub == cc.num_sub_steps-1 ) )
    {
        outfile << out_tmp.rdbuf();
        for(int ii = 0 ; ii < num_comp ; ii++)
        {
            out_diff[ii] << out_diff_tmp[ii].rdbuf();
            out_diff_tmp[ii].str( std::string() );
            out_diff_tmp[ii].clear();
        }
        out_tmp.str( std::string() );
        out_tmp.clear();
    }
}

/** @} */ // end of group io

#endif
