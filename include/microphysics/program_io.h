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
        table.dx2 = table.x2[1] - table.x2[0];
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
    nc.ascent_flag_var.getVar(startp, countp, &nc.ascent_flag);
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
  in.snapshot_index = 1000;
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
  in.fixed_iteration = true;
  in.auto_type = 1;
  in.traj = 0;
}

/**
 * String used to parse commandline input.
 */
static const char *optString = "f:d:i:b:o:l:s:t:a:r:?";


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

  if(1 == arg.traj_flag){
    in.traj = std::stoi(arg.traj_string);
  }
}

/**
 * Setup the cloud autoconversion parameters.
 *
 * @param pc Model constants for a certain particle type.
 */
void setup_cloud_autoconversion(
    particle_model_constants_t &pc)
{
    auto nu = pc.nu + 1.0;
    auto mu = pc.mu;
    if(pc.mu == 1.0)
    {
        cloud_k_au = kc_autocon / pc.max_x * 0.05
            * (nu+1.0)*(nu+3.0) / pow(nu, 2);
        cloud_k_sc = kc_autocon * (nu+1.0)/(nu);
    } else
    {
        cloud_k_au = kc_autocon / pc.max_x * 0.05
            * (2.0 * tgamma((nu+3.0)/mu)
            * tgamma((nu+1.0)/mu) * pow(tgamma((nu)/mu), 2)
            - 1.0 * pow(tgamma((nu+2.0)/mu), 2) * pow(tgamma((nu)/mu), 2))
            / pow(tgamma((nu+1.0)/mu), 4);
        cloud_k_sc = kc_autocon * pc.c_z;
    }
}

/**
 * Setup for bulk sedimentation velocity.
 *
 * @param pc Model constants for a certain particle type.
 */
void setup_bulk_sedi(
    particle_model_constants_t &pc)
{
    pc.alfa_n = pc.a_vel * tgamma( (pc.nu+pc.b_vel+1.0)/pc.mu )
        / tgamma( (pc.nu+1.0)/pc.mu);
    pc.alfa_q = pc.a_vel * tgamma( (pc.nu+pc.b_vel+2.0)/pc.mu )
        / tgamma( (pc.nu+2.0)/pc.mu );
    pc.lambda = tgamma( (pc.nu+1.0)/pc.mu )
        / tgamma( (pc.nu+2.0)/pc.mu );
}

/** @} */ // end of group io

#endif
