#ifndef PROGRAM_IO_H
#define PROGRAM_IO_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
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

#include "include/misc/general.h"
#include "include/misc/error.h"
#include "include/types/global_args_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/output_handle_t.h"
#include "include/types/model_constants_t.h"
#include "include/types/nc_parameters_t.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/segment_t.h"
#include "include/types/table_t.h"

#include "include/misc/general.h"
using namespace netCDF;

namespace pt = boost::property_tree;

/** @defgroup io Input/Output Functions
 * Functions to load data and setup and initialize structs.
 * @{
 */


/**
 * Parse an ensemble configuration as an XML-file.
 *
 */
int load_ens_config(
    std::string filename,
    model_constants_t &cc,
    std::vector<segment_t> &segments,
    input_parameters_t &input,
    const reference_quantities_t &ref_quant)
{
    int err = 0;
    pt::ptree pt;
    boost::property_tree::read_json(filename, pt);

    for(auto &it: pt.get_child("segments"))
    {
        segment_t segment;
        SUCCESS_OR_DIE(segment.from_pt(it.second, cc));
        segments.push_back(segment);
        if(segment.activated)
            segment.perturb(cc, ref_quant, input);
    }
    return err;
}


/**
 * Reads a checkpoint file including perturbed parameters if any, current
 * time step, new id as string.
 */
// template<class float_t>
// int load_checkpoint(
    // const std::string &filename,
    // model_constants_t &cc,
    // std::vector<float_t> &y,
    // std::vector<segment_t> &segments,
    // input_parameters_t &input,
    // const reference_quantities_t &ref_quant)
// {
//     int err = 0;
//     pt::ptree pt;
//     boost::property_tree::read_json(filename, pt);
//     // Parse the input parameters
//     SUCCESS_OR_DIE(input.from_pt(pt));
//     cc.setup_model_constants(input, ref_quant);
//     // Parse the model constants
//     SUCCESS_OR_DIE(cc.from_pt(pt));
//     // cc.id are the preceeding ids
//     input.set_outputfile_id(cc.id, cc.ensemble_id);

//     // Parse the segments
//     for(auto &it: pt.get_child("segments"))
//     {
//         segment_t segment;
//         SUCCESS_OR_DIE(segment.from_pt(it.second, cc));

//         if(segment.activated)
//         {
//             segment.perturb(cc, ref_quant, input);
//         }

//         segments.push_back(segment);
//     }
//     for(auto &it: pt.get_child("Output Parameters"))
//     {
//         y[std::stoi(it.first)] = std::stod(it.second.data());
//     }

//     return err;
// }




/**
 * Creates a job-script for MOGON II and puts it on queue.
 * If this is not possible (i.e. running on a local machine),
 * the next jobs are started on the local machine.
 * Also creates an XML-file as checkpoint for initialization to read the
 * current status of the simulation for the next job(s).
 *
 * @param filename Filename for the checkpoint file. On out: added _idx-y.json
 */
// template<class float_t>
// void write_checkpoint(
//     std::string &filename,
//     model_constants_t &cc,
//     const std::vector<float_t> &y,
//     std::vector<segment_t> &segments,
//     const input_parameters_t &input,
//     const double &current_time)
// {
//     pt::ptree checkpoint;
//     // First we add the ensemble configuration
//     pt::ptree segment_tree;
//     for(auto &s: segments)
//         s.put(segment_tree);

//     checkpoint.add_child("segments", segment_tree);
//     // configuration from input_parameters_t
//     input.put(checkpoint, current_time);
//     // Model constants
//     cc.put(checkpoint);
//     // Current status of y
//     pt::ptree output_parameters;
//     for(uint32_t i=0; i<num_comp; i++)
//         output_parameters.put(std::to_string(i), y[i]);
//     checkpoint.add_child("Output Parameters", output_parameters);

//     uint64_t i = 0;
//     std::string actual_filename = filename + "/checkpoint_id" + cc.id + "_0000.json";
//     while(exists(actual_filename))
//     {
//         i++;
//         if(i < 10)
//             actual_filename = filename + "/checkpoint_id" + cc.id + "_000" + std::to_string(i) + ".json";
//         else if(i < 100)
//             actual_filename = filename + "/checkpoint_id" + cc.id + "_00" + std::to_string(i) + ".json";
//         else if(i < 1000)
//             actual_filename = filename + "/checkpoint_id" + cc.id + "_0" + std::to_string(i) + ".json";
//         else
//             actual_filename = filename + "/checkpoint_id" + cc.id + "_" + std::to_string(i) + ".json";
//     }
//     std::fstream outstream(actual_filename, std::ios::out);
//     filename = actual_filename;
//     pt::write_json(outstream, checkpoint);
//     outstream.close();
//     // deactivate all segments, so we know, another instance is going
//     // to process this
//     for(auto &s: segments)
//         s.deactivate(true);

//     cc.ensemble_id++;
// }


/**
 * Create a bash script that executes this very same program by loading
 * a checkpoint file
 * (gnuparallel) starting n processes
 * (slurm) submitting it as a slurm job with n processes
 */
void create_run_script(
    const std::string &foldername,
    const std::string &checkpoint_file,
    const model_constants_t &cc,
    const uint32_t n_processes,
    const std::string which="gnuparallel",
    const bool start=false)
{
    if(which == "gnuparallel")
    {
        uint32_t j = n_processes - 1;
        if(j>NPROCS)
            j = NPROCS;
        std::string script = "parallel -u -j " + std::to_string(j)
            + " --no-notice --delay .2 ${AD_SIM_HOME}/build/apps/src/microphysics/./trajectories "
            "-c " + checkpoint_file + " -g {1} ::: {0.."
            + std::to_string(n_processes-2) + "}";

        // Save as a script
        std::string script_name = foldername + "/execute_id" + cc.id + "_0000.sh";
        uint32_t i = 0;
        while(exists(script_name))
        {
            i++;
            if(i < 10)
                script_name = foldername + "/execute_id" + cc.id + "_000" + std::to_string(i) + ".sh";
            else if(i < 100)
                script_name = foldername + "/execute_id" + cc.id + "_00" + std::to_string(i) + ".sh";
            else if(i < 1000)
                script_name = foldername + "/execute_id" + cc.id + "_0" + std::to_string(i) + ".sh";
            else
                script_name = foldername + "/execute_id" + cc.id + "_" + std::to_string(i) + ".sh";
        }
        std::ofstream out(script_name);
        out << "#!/bin/bash\n" << script << ( (start) ? " &\n" : "\n" );
        out.close();

        if(start)
            std::system(("./" + script_name).c_str());

    } else if(which == "slurm")
    {

    }
}



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
                for(int i=0; i<anzT_wg_loc; ++i)
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
            for(int j=0; j<anzT_wg_loc-1; ++j)
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
        std::cerr << "Error loading " << filename << ". Does the file exist?\n";
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
// void init_nc_parameters(
//     nc_parameters_t &nc,
//     uint32_t n=32,
//     uint32_t n_timesteps=7922)
// {
//     nc.n_trajectories = n;
//     nc.n_timesteps = n_timesteps;
// #if defined MET3D
//     nc.w.resize(2);
//     nc.time_abs.resize(n_timesteps);
//     nc.type[0] = (char*) "";
// #else
//     nc.w.resize(2);
//     nc.z.resize(4);
// #endif
//     nc.lat.resize(2);
//     nc.lon.resize(2);
// }


/**
 * Load variables from the netCDF file such that data can be loaded using
 * load_nc_parameters()
 *
 * @param nc Struct where to load the variables.
 * @param datafile The netCDF file.
 */
// void load_nc_parameters_var(
//     nc_parameters_t &nc,
//     NcFile &datafile)
// {
// #ifdef WCB2
//     nc.lat_var      = datafile.getVar("latitude");
//     nc.lon_var      = datafile.getVar("longitude");
// #else
//     nc.lat_var      = datafile.getVar("lat");
//     nc.lon_var      = datafile.getVar("lon");
// #endif
//     nc.z_var        = datafile.getVar("z");
// #if defined WCB || defined WCB2
//     nc.p_var        = datafile.getVar("P");
//     nc.t_var        = datafile.getVar("T");
//     nc.qc_var       = datafile.getVar("QC");
//     nc.qv_var       = datafile.getVar("QV");
//     nc.qr_var       = datafile.getVar("QR");
//     nc.qi_var       = datafile.getVar("QI");
//     nc.qs_var       = datafile.getVar("QS");
//     nc.time_rel_var = datafile.getVar("time");
// #elif defined MET3D
//     nc.p_var        = datafile.getVar("pressure");
//     nc.t_var        = datafile.getVar("T");
//     nc.qc_var       = datafile.getVar("QC");
//     nc.qv_var       = datafile.getVar("QV");
//     nc.qr_var       = datafile.getVar("QR");
//     nc.qi_var       = datafile.getVar("QI");
//     nc.qs_var       = datafile.getVar("QS");
//     nc.time_rel_var = datafile.getVar("time_after_ascent");
//     nc.time_abs_var = datafile.getVar("time");
//     nc.w_var        = datafile.getVar("w");
//     nc.S_var        = datafile.getVar("S");
//     nc.type_var     = datafile.getVar("type");
// #else
//     nc.p_var        = datafile.getVar("p");
//     nc.t_var        = datafile.getVar("t");
//     nc.w_var        = datafile.getVar("w");
//     nc.time_rel_var = datafile.getVar("time_rel");
//     nc.qc_var       = datafile.getVar("qc");
//     nc.qr_var       = datafile.getVar("qr");
//     nc.qi_var       = datafile.getVar("qi");
//     nc.qs_var       = datafile.getVar("qs");
//     nc.qg_var       = datafile.getVar("qg");
//     nc.qv_var       = datafile.getVar("qv");
//     nc.QIin_var     = datafile.getVar("QIin");
//     nc.QSin_var     = datafile.getVar("QSin");
//     nc.QRin_var     = datafile.getVar("QRin");
//     nc.QGin_var     = datafile.getVar("QGin");
//     nc.QIout_var    = datafile.getVar("QIout");
//     nc.QSout_var    = datafile.getVar("QSout");
//     nc.QRout_var    = datafile.getVar("QRout");
//     nc.QGout_var    = datafile.getVar("QGout");
// #endif
// #ifdef WCB
//     // specific humidity
//     nc.S_var        = datafile.getVar("RELHUM");
//     // Flag wether an effective ascent region is reached
//     nc.ascent_flag_var = datafile.getVar("MAP");
//     // Potential vorticity (German: Wirbelstaerke)
//     // nc.pot_vortic   = datafile.getVar("POT_VORTIC")
// #endif

// #if defined WCB2 || defined MET3D
//     nc.qg_var       = datafile.getVar("QG");
//     nc.QIin_var     = datafile.getVar("QI_IN");
//     nc.QSin_var     = datafile.getVar("QS_IN");
//     nc.QRin_var     = datafile.getVar("QR_IN");
//     nc.QGin_var     = datafile.getVar("QG_IN");
//     nc.QIout_var    = datafile.getVar("QI_OUT");
//     nc.QSout_var    = datafile.getVar("QS_OUT");
//     nc.QRout_var    = datafile.getVar("QR_OUT");
//     nc.QGout_var    = datafile.getVar("QG_OUT");

//     nc.NIin_var     = datafile.getVar("NI_IN");
//     nc.NSin_var     = datafile.getVar("NS_IN");
//     nc.NRin_var     = datafile.getVar("NR_IN");
//     nc.NGin_var     = datafile.getVar("NG_IN");
//     nc.NIout_var    = datafile.getVar("NI_OUT");
//     nc.NSout_var    = datafile.getVar("NS_OUT");
//     nc.NRout_var    = datafile.getVar("NR_OUT");
//     nc.NGout_var    = datafile.getVar("NG_OUT");

//     nc.Nc_var       = datafile.getVar("NCCLOUD");
//     nc.Nr_var       = datafile.getVar("NCRAIN");
//     nc.Ni_var       = datafile.getVar("NCICE");
//     nc.Ns_var       = datafile.getVar("NCSNOW");
//     nc.Ng_var       = datafile.getVar("NCGRAUPEL");

//     nc.conv_400_var = datafile.getVar("conv_400");
//     nc.conv_600_var = datafile.getVar("conv_600");
//     nc.slan_400_var = datafile.getVar("slan_400");
//     nc.slan_600_var = datafile.getVar("slan_600");
// #endif
// #if defined WCB2
//     // Flag wether an effective ascent region is reached
//     nc.ascent_flag_var = datafile.getVar("WCB_flag");
//     // 2h ascent rate after Oertel et al. (2019)
//     nc.dp2h_var     = datafile.getVar("dp2h");
// #endif
// #if defined MET3D && defined TURBULENCE
//     nc.qturb_var    = datafile.getVar("Q_TURBULENCE");
// #endif
// }


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
// void load_nc_parameters(
//     nc_parameters_t &nc,
//     std::vector<size_t> &startp,
//     std::vector<size_t> &countp,
//     const reference_quantities_t &ref_quant,
//     uint64_t num_sub_steps)
// {
//     // startp[0] <- ensemble (if available)
//     // startp[1] <- trajectory id
//     // startp[2] <- timestep


// #if defined WCB || defined WCB2
//     countp[0]++;
//     countp[0]++;
//     nc.z_var.getVar(startp, countp, nc.z.data());
//     nc.lat_var.getVar(startp, countp, nc.lat.data());
//     nc.lon_var.getVar(startp, countp, nc.lon.data());
//     countp[0]--;
//     countp[0]--;
// #elif defined MET3D
//     nc.z_var.getVar(startp, countp, &nc.z);
//     countp[2]++;
//     nc.lat_var.getVar(startp, countp, nc.lat.data());
//     nc.lon_var.getVar(startp, countp, nc.lon.data());
//     nc.w_var.getVar(startp, countp, nc.w.data());
//     countp[2]--;
// #else
//     countp[1]++;
//     countp[1]++;
//     nc.z_var.getVar(startp, countp, nc.z.data());
//     nc.lat_var.getVar(startp, countp, nc.lat.data());
//     nc.lon_var.getVar(startp, countp, nc.lon.data());
//     countp[1]--;
//     countp[1]--;
// #endif

// #if defined WCB || defined WCB2
//     int map = 0;
//     nc.ascent_flag_var.getVar(startp, countp, &map);
//     nc.ascent_flag = (map > 0) ? true : false;
// #endif
// #if defined WCB2
//     nc.dp2h_var.getVar(startp, countp, &map);
//     nc.dp2h = (map > 0) ? true : false;
// #endif
// #if defined MET3D
//     int map = 0;
// #endif
// #if defined MET3D || defined WCB2
//     nc.conv_400_var.getVar(startp, countp, &map);
//     nc.conv_400 = (map > 0) ? true : false;
//     nc.conv_600_var.getVar(startp, countp, &map);
//     nc.conv_600 = (map > 0) ? true : false;
//     nc.slan_400_var.getVar(startp, countp, &map);
//     nc.slan_400 = (map > 0) ? true : false;
//     nc.slan_600_var.getVar(startp, countp, &map);
//     nc.slan_600 = (map > 0) ? true : false;
// #endif

//     nc.t_var.getVar(startp, countp, &nc.t);
//     nc.p_var.getVar(startp, countp, &nc.p);
//     nc.time_rel_var.getVar(startp, countp, &nc.time_rel);
//     nc.qc_var.getVar(startp, countp, &nc.qc);
//     nc.qr_var.getVar(startp, countp, &nc.qr);
//     nc.qi_var.getVar(startp, countp, &nc.qi);
//     nc.qs_var.getVar(startp, countp, &nc.qs);
//     nc.qv_var.getVar(startp, countp, &nc.qv);

// #if !defined WCB2 && !defined MET3D
//     // We are reading in hPa. Convert to Pa
//     nc.p        *= 100;
// #endif
//     nc.p        /= ref_quant.pref;
//     nc.t        /= ref_quant.Tref;
//     nc.qc       /= ref_quant.qref;
//     nc.qr       /= ref_quant.qref;
//     nc.qv       /= ref_quant.qref;
//     nc.qi       /= ref_quant.qref;
//     nc.qs       /= ref_quant.qref;

// #if defined WCB2 || defined MET3D
//     nc.NRin_var.getVar(startp, countp, &nc.NRin);
//     nc.NIin_var.getVar(startp, countp, &nc.NIin);
//     nc.NSin_var.getVar(startp, countp, &nc.NSin);
//     nc.NGin_var.getVar(startp, countp, &nc.NGin);

//     nc.NRout_var.getVar(startp, countp, &nc.NRout);
//     nc.NIout_var.getVar(startp, countp, &nc.NIout);
//     nc.NSout_var.getVar(startp, countp, &nc.NSout);
//     nc.NGout_var.getVar(startp, countp, &nc.NGout);

//     nc.Nc_var.getVar(startp, countp, &nc.Nc);
//     nc.Nr_var.getVar(startp, countp, &nc.Nr);
//     nc.Ni_var.getVar(startp, countp, &nc.Ni);
//     nc.Ns_var.getVar(startp, countp, &nc.Ns);
//     nc.Ng_var.getVar(startp, countp, &nc.Ng);

//     nc.NRin     /= ref_quant.Nref;
//     nc.NIin     /= ref_quant.Nref;
//     nc.NSin     /= ref_quant.Nref;
//     nc.NGin     /= ref_quant.Nref;

//     nc.NRout    /= ref_quant.Nref;
//     nc.NIout    /= ref_quant.Nref;
//     nc.NSout    /= ref_quant.Nref;
//     nc.NGout    /= ref_quant.Nref;

//     nc.Nc       /= ref_quant.Nref;
//     nc.Nr       /= ref_quant.Nref;
//     nc.Ni       /= ref_quant.Nref;
//     nc.Ns       /= ref_quant.Nref;
//     nc.Ng       /= ref_quant.Nref;

//     nc.NRin     = abs(nc.NRin);
//     nc.NIin     = abs(nc.NIin);
//     nc.NSin     = abs(nc.NSin);
//     nc.NGin     = abs(nc.NGin);
// #endif

// #if defined MET3D && defined TURBULENCE
//     nc.qturb_var.getVar(startp, countp, &nc.qturb);
//     nc.qturb /= ref_quant.qref;
// #endif
// #if !defined WCB
//     nc.S = 1.0;
//     nc.qg_var.getVar(startp, countp, &nc.qg);
//     nc.QIin_var.getVar(startp, countp, &nc.QIin);
//     nc.QSin_var.getVar(startp, countp, &nc.QSin);
//     nc.QRin_var.getVar(startp, countp, &nc.QRin);
//     nc.QGin_var.getVar(startp, countp, &nc.QGin);
//     nc.QIout_var.getVar(startp, countp, &nc.QIout);
//     nc.QSout_var.getVar(startp, countp, &nc.QSout);
//     nc.QRout_var.getVar(startp, countp, &nc.QRout);
//     nc.QGout_var.getVar(startp, countp, &nc.QGout);

//     nc.qg       /= ref_quant.qref;
//     nc.QRin     /= ref_quant.qref;
//     nc.QRout    /= ref_quant.qref;
//     nc.QIin     /= ref_quant.qref;
//     nc.QIout    /= ref_quant.qref;
//     nc.QSin     /= ref_quant.qref;
//     nc.QSout    /= ref_quant.qref;
//     nc.QGin     /= ref_quant.qref;
//     nc.QGout    /= ref_quant.qref;

//     nc.QRin     = abs(nc.QRin);
//     nc.QRout    = abs(nc.QRout);
//     nc.QIin     = abs(nc.QIin);
//     nc.QIout    = abs(nc.QIout);
//     nc.QSin     = abs(nc.QSin);
//     nc.QSout    = abs(nc.QSout);
//     nc.QGin     = abs(nc.QGin);
//     nc.QGout    = abs(nc.QGout);
// #endif
// #ifdef MET3D
//     nc.S_var.getVar(startp, countp, &nc.S);
//     nc.S /= 100; // from percentage
//     nc.time_rel_var.getVar(startp, countp, &nc.time_rel);
//     // there is only a single value for that since each type
//     // is divided into different files in the input
//     if(std::strcmp(nc.type[0], "") == 0)
//     {
//         nc.type_var.getVar(nc.type);
//     }
//     nc.w[0]     /= ref_quant.wref;
//     nc.w[1]     /= ref_quant.wref;
// #elif !defined WCB && !defined WCB2
//     nc.w_var.getVar(startp, countp, nc.w.data());
//     countp[1]++;
//     nc.w_var.getVar(startp, countp, nc.w.data());
//     countp[1]--;

//     nc.w[0]     /= ref_quant.wref;
//     nc.w[1]     /= ref_quant.wref;
// #elif defined WCB
//     nc.qc /= 1.0e6;
//     nc.qr /= 1.0e6;
//     nc.qv /= 1.0e6;
//     nc.qi /= 1.0e6;
//     nc.qs /= 1.0e6;
//     // Calculate w by getting the z-coordinates
//     // and divide it by the amount of substeps
//     nc.w[0] = (nc.z[1] - nc.z[0]) / 20.0;
//     nc.w[1] = (nc.z[2] - nc.z[1]) / 20.0;
//     nc.S_var.getVar(startp, countp, &nc.S);
//     nc.S /= 100; // from percentage
// #else
//     // WCB 2
//     // Calculate w by getting the z-coordinates
//     // and divide it by the amount of substeps
//     nc.w[0] = (nc.z[1] - nc.z[0]) / 20.0;
//     nc.w[1] = (nc.z[2] - nc.z[1]) / 20.0;
// #endif
//     nc.dlat = (nc.lat[1] - nc.lat[0]) / num_sub_steps;
//     nc.dlon = (nc.lon[1] - nc.lon[0]) / num_sub_steps;
// }


/**
 * Set the input parameters with the data from the global arguments.
 *
//  * @param arg Stuct with command line arguments.
//  * @param in Struct where the input parameters will be stored.
//  */
// void set_input_from_arguments(global_args_t &arg ,
// 			      input_parameters_t &in )
// {
//   // Final time
//   if(1 == arg.final_time_flag){
//     in.t_end_prime = std::strtod(arg.final_time_string, nullptr);
//   }

//   // Timestep
//   if(1 == arg.timestep_flag){
//     in.dt_prime = std::strtod(arg.timestep_string, nullptr);
//   }

//   // Snapshot index
//   if(1 == arg.snapshot_index_flag){
//     in.snapshot_index = std::stoi(arg.snapshot_index_string);
//   }

//   // Output
//   if(1 == arg.output_flag){
//     in.OUTPUT_FILENAME = arg.output_string;
//   }

//   // Input
//   if(1 == arg.input_flag){
//     in.INPUT_FILENAME = arg.input_file;
//   }

//   // Scaling factor
//   if(1 == arg.scaling_fact_flag){
//     in.scaling_fact = std::strtod(arg.scaling_fact_string, nullptr);
//   }

//   // Starting over mixing ratios and particle numbers
//   if(1 == arg.start_over_flag){
//     in.start_over = (strcmp(arg.start_over_string, "0"));
//   }

//   // Starting over environment variables (p, T, w)
//   if(1 == arg.start_over_env_flag){
//     in.start_over_env = (strcmp(arg.start_over_env_string, "0"));
//   }

//   if(1 == arg.fixed_iteration_flag){
//       in.fixed_iteration = (strcmp(arg.fixed_iteration_string, "0"));
//   }

//   // Auto type
//   if(1 == arg.auto_type_flag){
//     in.auto_type = std::stoi(arg.auto_type_string);
//   }

//   // Trajectory
//   if(1 == arg.traj_flag){
//     in.traj = std::stoi(arg.traj_string);
//   }

//   // Write index
//   if(1 == arg.write_flag){
//       in.write_index = std::stoi(arg.write_string);
//   }

//   // Progressbar index
//   if(1 == arg.progress_index_flag){
//       in.progress_index = std::stoull(arg.progress_index_string);
//   }
// #ifdef MET3D
//   // Simulation start time
//   if(1 == arg.delay_start_flag){
//       in.start_time = std::strtod(arg.delay_start_string, nullptr);
//   }
// #endif

//   // Ensemble configuration file
//   if(1 == arg.ens_config_flag){
//     in.ENS_CONFIG_FILENAME = arg.ens_config_string;
//   }

//   // Checkpoint file
//   if(1 == arg.checkpoint_flag){
//     in.CHECKPOINT_FILENAME = arg.checkpoint_string;
//   }

//   // ID for this process
//   if(1 == arg.gnu_id_flag){
//     in.id = std::stoi(arg.gnu_id_string);
//   }

//   // Folder name for new generated checkpoints
//   if(1 == arg.folder_name_flag){
//     in.FOLDER_NAME = arg.folder_name_string;
//   }
// }


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
#ifdef MET3D
    double &start_time,
#endif
    const char *input_file,
    const uint32_t traj,
    const bool checkpoint_flag,
    model_constants_t &cc,
    double current_time)
{
    try
    {
        netCDF::NcFile datafile(input_file, netCDF::NcFile::read);
#ifdef WCB
        lenp = datafile.getDim("ntra").getSize();
#elif defined MET3D
        lenp = datafile.getDim("trajectory").getSize();
#else
        lenp = datafile.getDim("id").getSize();
#endif
#ifndef SILENT_MODE
        std::cout << "Number of trajectories in netCDF file: " << lenp << "\n" << std::flush;
#endif
        if(lenp <= traj)
        {
            std::cerr << "You asked for trajectory with index " << traj
                      << " which does not exist. ABORTING.\n";
            return NC_TRAJ_IDX_ERR;
        }
#ifdef WCB
        uint32_t n_timesteps = datafile.getDim("ntim").getSize();
#else
        uint32_t n_timesteps = datafile.getDim("time").getSize();
#endif

        nc_params.init_params(lenp, n_timesteps);
        nc_params.load_vars(datafile);
#ifdef MET3D
        // Get the time coordinates
        nc_params.time_abs_var.getVar(nc_params.time_abs.data());
#endif
        std::vector<size_t> startp, countp;
        countp.push_back(1);
        countp.push_back(1);
#ifdef MET3D
        countp.push_back(1);
#endif
        // wcb files have a different ordering
#if defined WCB || defined WCB2
        startp.push_back(0); // time
        startp.push_back(traj); // trajectory id
#elif defined MET3D

        startp.push_back(0); // ensemble
        startp.push_back(traj); // trajectory
        startp.push_back(0); // time
        uint64_t start_time_idx = 0;

        if(!std::isnan(start_time) && !checkpoint_flag)
        {
            double rel_start_time;
            // Get relative start time of trajectory
            nc_params.time_rel_var.getVar(startp, countp, &rel_start_time);
            // Calculate the needed index
            start_time_idx = (start_time-rel_start_time)/cc.dt_traject;
            uint64_t n_timesteps_input = ceil((cc.t_end-start_time)/20.0)-1;
            cc.num_steps = (n_timesteps-1 > n_timesteps_input) ? n_timesteps_input : n_timesteps-1;
        } else if(checkpoint_flag && !std::isnan(current_time) && std::isnan(start_time))
        {
            start_time_idx = ceil(current_time/cc.dt_traject);
            uint64_t n_timesteps_input = ceil((cc.t_end-current_time)/20.0)-1;
            cc.num_steps = (n_timesteps-1 > n_timesteps_input) ? n_timesteps_input : n_timesteps-1;
        } else if(checkpoint_flag && !std::isnan(current_time) && !std::isnan(start_time))
        {
            double rel_start_time;
            // Get relative start time of trajectory
            nc_params.time_rel_var.getVar(startp, countp, &rel_start_time);
            // Calculate the needed index
            start_time_idx = ceil(start_time-rel_start_time + current_time)/cc.dt_traject;
            uint64_t n_timesteps_input = ceil((cc.t_end-current_time-start_time)/20.0)-1;
            cc.num_steps = (n_timesteps-1 > n_timesteps_input) ? n_timesteps_input : n_timesteps-1;
        } else
        {
            uint64_t n_timesteps_input = ceil(cc.t_end/20.0)-1 - start_time_idx;
            cc.num_steps = (n_timesteps-1 > n_timesteps_input) ? n_timesteps_input : n_timesteps-1;
        }

        nc_params.time_idx = start_time_idx;
        startp[2] = start_time_idx;
#else
        startp.push_back(traj); // trajectory id
        startp.push_back(1); // time (where time == 0 only has zeros)
#endif
        nc_params.load_params(startp, countp,
                           ref_quant, cc.num_sub_steps);
        if(!checkpoint_flag)
        {
            y_init[p_idx] = nc_params.p;
            y_init[T_idx] = nc_params.t;

            y_init[S_idx] = nc_params.S;
#ifdef SAT_CALC
            y_init[S_idx] = nc_params.qv*ref_quant.qref * Rv * nc_params.t*ref_quant.Tref
                / saturation_pressure_water(nc_params.t*ref_quant.Tref);
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
#if defined WCB2 || defined MET3D
            y_init[Nr_idx] = nc_params.Nr;
            y_init[Nc_idx] = nc_params.Nc;
#if defined(RK4ICE)
            // We can read the sedimentation from the original file like this
            // y_init[qr_out_idx] = nc_params.QRout;
            // y_init[Nr_out_idx] = nc_params.NRout;
            // y_init[qi_out_idx] = nc_params.QIout;
            // y_init[qs_out_idx] = nc_params.QSout;
            // y_init[qg_out_idx] = nc_params.QGout;
            // y_init[Ni_out_idx] = nc_params.NIout;
            // y_init[Ns_out_idx] = nc_params.NSout;
            // y_init[Ng_out_idx] = nc_params.NGout;
            // But we actually want to set them to zero and calculate them
            y_init[qr_out_idx] = 0;
            y_init[Nr_out_idx] = 0;
            y_init[qi_out_idx] = 0;
            y_init[qs_out_idx] = 0;
            y_init[qg_out_idx] = 0;
            y_init[qh_out_idx] = 0;
            y_init[Nh_out_idx] = 0;
            y_init[Ni_out_idx] = 0;
            y_init[Ns_out_idx] = 0;
            y_init[Ng_out_idx] = 0;

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
            y_init[qh_out_idx] = 0;
            y_init[Nh_out_idx] = 0;
            y_init[Ni_out_idx] = 0;
            y_init[Ns_out_idx] = 0;
            y_init[Ng_out_idx] = 0;
            y_init[Ni_idx] = 0;
            y_init[Ns_idx] = 0;
            y_init[Ng_idx] = 0;
#endif
#endif
#if defined(RK4ICE) || defined(RK4NOICE)
    #ifdef MET3D
            y_init[z_idx] = nc_params.z;
    #else
            y_init[z_idx] = nc_params.z[0];
    #endif
#endif
#if defined(RK4ICE) || defined(RK4NOICE) || defined(MET3D)
            nc_params.dw = (nc_params.w[1] - nc_params.w[0]);
            cc.constants[static_cast<int>(Cons_idx::dw)] = nc_params.dw / (cc.dt*cc.num_sub_steps);
            y_init[n_inact_idx] = 0;
            y_init[depo_idx] = 0;
            y_init[sub_idx] = 0;
#endif

#ifdef SAT_CALC
            y_init[S_idx]  = convert_qv_to_S(
                y_init[p_idx]*ref_quant.pref,
                y_init[T_idx]*ref_quant.Tref,
                y_init[qv_idx]*ref_quant.qref,
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                 get_at(cc.constants, Cons_idx::Epsilon));
#endif
        }
    } catch(netCDF::exceptions::NcException& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "In read_init_netcdf()\nABORTING." << std::endl;
        return 1;
    }
    return 0;
}



/**
 * Resize startp and countp.
 *
 * @param startp On out: contains dimensions info for reading
 * @params countp On out: contains dimensions info for reading
 * @params input_file Char array of input file
 */
void resize_counter(
    std::vector<size_t> &startp,
    std::vector<size_t> &countp,
    const uint32_t traj)
{
#if defined WCB || defined WCB2
    startp.push_back(1);          // time
    startp.push_back(traj); // trajectory
#elif defined MET3D
    startp.push_back(0); // ensemble; currently only one is supported
    startp.push_back(traj); // trajectory
    startp.push_back(1); // time
#else
    startp.push_back(traj); // trajectory
    startp.push_back(1);          // time
#endif
    countp.push_back(1);
    countp.push_back(1);
#ifdef MET3D
    countp.push_back(1);
#endif
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
    const reference_quantities_t &ref_quant,
    std::vector<codi::RealReverse> &y_single_old,
    std::vector<codi::RealReverse> &inflow,
    std::vector<int> &ids,
    int &traj_id,
#ifdef MET3D
    uint32_t &ensemble,
#endif
    const uint32_t t,
    const bool checkpoint_flag,
    output_handle_t &io_handler)
{
#if defined WCB || defined WCB2
    startp[0] = t;
#elif defined MET3D
    startp[2] = t + nc_params.time_idx;
#else
    startp[1] = t+1;
#endif
    netCDF::NcFile datafile(input_file, netCDF::NcFile::read);
    nc_params.load_vars(datafile);
    nc_params.load_params(startp, countp,
                        ref_quant, cc.num_sub_steps);
    netCDF::NcVar id_var;
#ifdef MET3D
    id_var = datafile.getVar("trajectory");
#else
    id_var = datafile.getVar("id");
#endif
    id_var.getVar(ids.data());
    traj_id = ids[input.traj];
#ifdef MET3D
    id_var = datafile.getVar("ensemble");
    id_var.getVar(ids.data());
    ensemble = ids[input.ensemble];
#endif
    // Reset outflow
    y_single_old[qi_out_idx] = 0;
    y_single_old[qs_out_idx] = 0;
    y_single_old[qr_out_idx] = 0;
    y_single_old[qg_out_idx] = 0;
    y_single_old[qh_out_idx] = 0;
    y_single_old[Ni_out_idx] = 0;
    y_single_old[Ns_out_idx] = 0;
    y_single_old[Nr_out_idx] = 0;
    y_single_old[Ng_out_idx] = 0;
    y_single_old[Nh_out_idx] = 0;

    // Set values from a given trajectory
    if((t==0 && !checkpoint_flag) || input.start_over_env)
    {
        y_single_old[p_idx]  = nc_params.p;     // p
        y_single_old[T_idx]  = nc_params.t;     // T
#if defined(RK4ICE) || defined(RK4NOICE)
        y_single_old[w_idx]  = nc_params.w[0]; // w
        nc_params.dw = (nc_params.w[1] - nc_params.w[0]);
        cc.constants[static_cast<int>(Cons_idx::dw)] = nc_params.dw / (cc.dt*cc.num_sub_steps);
#endif
#ifdef MET3D
        y_single_old[z_idx] = nc_params.z;
#else
        y_single_old[z_idx] = nc_params.z[0];
#endif

#if defined(FLUX) && (defined(WCB2) || defined(MET3D))
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
#if defined MET3D && defined TURBULENCE
        inflow[qv_in_idx] = nc_params.qturb;
#endif
    }

    if((t==0 && !checkpoint_flag) || input.start_over)
    {
        y_single_old[S_idx]  = nc_params.S;     // S
#ifdef SAT_CALC
        y_single_old[S_idx]  = nc_params.qv*ref_quant.qref * Rv * nc_params.t*ref_quant.Tref
            / saturation_pressure_water(nc_params.t*ref_quant.Tref);
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
#if defined(RK4ICE) && (defined(WCB2) || defined(MET3D))
        y_single_old[Ng_idx] = nc_params.Ng;
        y_single_old[Ni_idx] = nc_params.Ni;
        y_single_old[Ns_idx] = nc_params.Ns;
        // Reading precipation can be done like this
        // y_single_old[Ng_out_idx] = nc_params.NGout;
        // y_single_old[Ni_out_idx] = nc_params.NIout;
        // y_single_old[Ns_out_idx] = nc_params.NSout;
        // y_single_old[Nr_out_idx] = nc_params.NRout;
        // But we calculate them by ourself
        y_single_old[Ng_out_idx] = 0;
        y_single_old[Ni_out_idx] = 0;
        y_single_old[Ns_out_idx] = 0;
        y_single_old[Nr_out_idx] = 0;
        y_single_old[Nh_out_idx] = 0;
        y_single_old[qr_out_idx] = 0;
        y_single_old[qi_out_idx] = 0;
        y_single_old[qs_out_idx] = 0;
        y_single_old[qg_out_idx] = 0;
        y_single_old[qh_out_idx] = 0;
#endif
#if defined WCB2 || defined MET3D
        y_single_old[qr_out_idx] = 0;
        y_single_old[Nr_out_idx] = 0;
        y_single_old[Nc_idx] = nc_params.Nc;
        y_single_old[Nr_idx] = nc_params.Nr;
#else
        denom = (get_at(cc.cloud.max_x - cc.cloud.constants, Particle_cons_idx::min_x)) / 2.0 + get_at(cc.cloud.constants, Particle_cons_idx::min_x);
        y_single_old[Nc_idx] = y_single_old[qc_idx] * ref_quant.qref / (denom); //*10e2);  // Nc
        denom = (get_at(cc.rain.max_x - cc.rain.constants, Particle_cons_idx::min_x)) / 2 + get_at(cc.rain.constants, Particle_cons_idx::min_x);
        y_single_old[Nr_idx] = y_single_old[qr_idx] * ref_quant.qref / (denom); //*10e2);  // Nr
        denom = get_at(cc.cloud.constants, Particle_cons_idx::min_x) / 2.0;
    #if defined(RK4ICE)
        denom = (get_at(cc.ice.max_x - cc.ice.constants, Particle_cons_idx::min_x)) / 2.0 + get_at(cc.ice.constants, Particle_cons_idx::min_x);
        y_single_old[Ni_idx] = y_single_old[qi_idx] * ref_quant.qref / (denom); //*10e2); // Ni
        denom = (get_at(cc.snow.max_x - cc.snow.constants, Particle_cons_idx::min_x)) / 2.0 + get_at(cc.snow.constants, Particle_cons_idx::min_x);
        y_single_old[Ns_idx] = y_single_old[qs_idx] * ref_quant.qref / (denom); //*10e2); // Ns
        denom = (get_at(cc.graupel.max_x - cc.graupel.constants, Particle_cons_idx::min_x)) / 2.0 + get_at(cc.graupel.constants, Particle_cons_idx::min_x);
        y_single_old[Ng_idx] = y_single_old[qg_idx] * ref_quant.qref / (denom); //*10e2); // Ng
    #endif
#endif
        // Actually only needed for single moment schemes
        cc.constants[static_cast<int>(Cons_idx::Nc_prime)] = y_single_old[Nc_idx];

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
#if defined(FLUX) && (defined(WCB2) || defined(MET3D))
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

#if defined MET3D && defined TURBULENCE
        inflow[qv_in_idx] = nc_params.qturb;
#endif

#ifdef SAT_CALC
        y_single_old[S_idx]  = convert_qv_to_S(
            y_single_old[p_idx].getValue()*ref_quant.pref,
            y_single_old[T_idx].getValue()*ref_quant.Tref,
            y_single_old[qv_idx].getValue()*ref_quant.qref,
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon));
#endif
        datafile.close();
        std::vector< std::array<double, num_par > >  y_diff(num_comp);
        for(auto &y_d: y_diff)
            y_d.fill(0);
        io_handler.buffer(cc, nc_params, y_single_old, y_diff, 0, t,
            (t*cc.num_sub_steps)*cc.dt, traj_id, ensemble, ref_quant);
    }
}


/**
 * Store simulation results and gradients to string stream and dump it
 * if timestep is reached to do so.
 */
// void write_output(
//     const model_constants_t &cc,
//     const nc_parameters_t &nc_params,
//     const std::vector<codi::RealReverse> &y_single_new,
//     const std::vector< std::array<double, num_par > >  &y_diff,
//     const uint32_t sub,
//     const uint32_t t,
//     const double time_new,
//     const uint32_t traj_id,
//     const uint32_t write_index,
//     const uint32_t snapshot_index,
// #ifdef MET3D
//     const uint32_t ensemble,
// #endif
//     const bool last_step,
//     output_handle_t &io_handler,
//     const reference_quantities_t &ref_quant)
// {
//     if( (0 == (sub + t*cc.num_sub_steps) % snapshot_index)
//         || ( t == cc.num_steps-1 && last_step ) )
//     {
//         io_handler.buffer(cc, nc_params, y_single_new, y_diff, sub, t,
//             time_new, traj_id, ensemble, ref_quant);

//     }

//     if( (0 == (sub + t*cc.num_sub_steps) % write_index)
//         || ( t == cc.num_steps-1 && last_step ) )
//     {
//         io_handler.flush_buffer();
//     }
// }


// /**
//  * Parse the arguments and store them in global_args.
//  *
//  * @param argc Number of arguments
//  * @param argv Pointer to arguments
//  * @param global_args Structure to store global args.
//  *
//  * @return Error code.
//  */
// int parse_arguments(
//     const int argc,
//     char* const * argv,
//     global_args_t &global_args,
//     const int &rank,
//     const int &n_processes)
// {
//     /**
//      * String used to parse commandline input.
//      */
//     static const char *optString = "w:f:d:e:i:b:o:l:s:t:a:r:p:n:m:c:g:h:?";
//     bool need_to_abort = false;
//     int opt;

//     if(argc < 2)
//     {
//         need_to_abort = true;
//         display_usage();
//     }else
//     {
//         opt = getopt(argc, argv, optString);

//         while(-1 != opt)
//         {
//             switch(opt)
//             {
//                 case 'f':
//                 {
//                     global_args.final_time_flag = 1;
//                     global_args.final_time_string = optarg;
//                     break;
//                 }
//                 case 'd':
//                 {
//                     global_args.timestep_flag = 1;
//                     global_args.timestep_string = optarg;
//                     break;
//                 }
//                 case 'i':
//                 {
//                     global_args.snapshot_index_flag = 1;
//                     global_args.snapshot_index_string = optarg;
//                     break;
//                 }
//                 case 'b':
//                 {
//                     global_args.scaling_fact_flag = 1;
//                     global_args.scaling_fact_string = optarg;
//                     break;
//                 }
//                 case 'o':
//                 {
//                     global_args.output_flag = 1;
//                     global_args.output_string = optarg;
//                     break;
//                 }
//                 case 'l':
//                 {
//                     global_args.input_flag = 1;
//                     global_args.input_file = optarg;
//                     break;
//                 }
//                 case 's':
//                 {
//                     global_args.start_over_flag = 1;
//                     global_args.start_over_string = optarg;
//                     break;
//                 }
//                 case 'e':
//                 {
//                     global_args.start_over_env_flag = 1;
//                     global_args.start_over_env_string = optarg;
//                     break;
//                 }
//                 case 't':
//                 {
//                     global_args.fixed_iteration_flag = 1;
//                     global_args.fixed_iteration_string = optarg;
//                     break;
//                 }
//                 case 'a':
//                 {
//                     global_args.auto_type_flag = 1;
//                     global_args.auto_type_string = optarg;
//                     break;
//                 }
//                 case 'r':
//                 {
//                     global_args.traj_flag = 1;
//                     global_args.traj_string = optarg;
//                     break;
//                 }
//                 case 'w':
//                 {
//                     global_args.write_flag = 1;
//                     global_args.write_string = optarg;
//                     break;
//                 }
//                 case 'p':
//                 {
//                     global_args.progress_index_flag = 1;
//                     global_args.progress_index_string = optarg;
//                     break;
//                 }
// #ifdef MET3D
//                 case 'n':
//                 {
//                     global_args.delay_start_flag = 1;
//                     global_args.delay_start_string = optarg;
//                     break;
//                 }
// #endif
//                 case 'm':
//                 {
//                     global_args.ens_config_flag = 1;
//                     global_args.ens_config_string = optarg;
//                     break;
//                 }
//                 case 'c':
//                 {
//                     global_args.checkpoint_flag = 1;
//                     global_args.checkpoint_string = optarg;
//                     break;
//                 }
//                 case 'g':
//                 {
//                     global_args.gnu_id_flag = 1;
//                     global_args.gnu_id_string = optarg;
//                     break;
//                 }
//                 case 'h':
//                 {
//                     global_args.folder_name_flag = 1;
//                     global_args.folder_name_string = optarg;
//                     break;
//                 }
//                 case '?':
//                 {
//                     need_to_abort = true;
//                     display_usage();
//                     break;
//                 }
//                 default:
//                 {
//                     need_to_abort = true;
//                     display_error_on_command_line();
//                     display_usage();
//                     break;
//                 }
//             }

//             opt = getopt(argc, argv, optString);
//         }
//     }

//     if(need_to_abort){
//         return ARGUMENT_ERR;
//     }
//     return SUCCESS;
// }

/** @} */ // end of group io

#endif
