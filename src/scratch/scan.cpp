#include "codi.hpp"
#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tgmath.h>
#include <vector>

#include "include/microphysics/constants.h"
#include "include/microphysics/general.h"
#include "include/microphysics/physical_parameterizations.h"
#include "include/microphysics/program_io.h"
#include "include/microphysics/types.h"
#include "include/microphysics/user_functions.h"


static struct option long_options[] =
{
    {"pressure", required_argument, NULL, 0},
    {"temperature", required_argument, NULL, 1},
    {"ascent", required_argument, NULL, 2},
    {"qv", required_argument, NULL, 3},
    {"qc", required_argument, NULL, 4},
    {"qr", required_argument, NULL, 5},
    {"qs", required_argument, NULL, 6},
    {"qi", required_argument, NULL, 7},
    {"qh", required_argument, NULL, 8},
    {"qg", required_argument, NULL, 9},
    {"temp_crit", required_argument, NULL, 10},
    {"n_inact", required_argument, NULL, 11},
    {"x_i", required_argument, NULL, 12},
    {"D_i", required_argument, NULL, 13},
    {"dep_snow", required_argument, NULL, 14},
    {"dep_ice", required_argument, NULL, 15},
    {"rime_qc", required_argument, NULL, 16},
    {"rime_qr", required_argument, NULL, 17},
    {"rime_qi", required_argument, NULL, 18},
    {"rime_qs", required_argument, NULL, 19},

    {"p_min", required_argument, NULL, 20},
    {"p_max", required_argument, NULL, 21},
    {"temp_min", required_argument, NULL, 22},
    {"temp_max", required_argument, NULL, 23},
    {"ascent_min", required_argument, NULL, 24},
    {"ascent_max", required_argument, NULL, 25},
    {"qv_min", required_argument, NULL, 26},
    {"qv_max", required_argument, NULL, 27},
    {"qc_min", required_argument, NULL, 28},
    {"qc_max", required_argument, NULL, 29},
    {"qr_min", required_argument, NULL, 30},
    {"qr_max", required_argument, NULL, 31},
    {"qs_min", required_argument, NULL, 32},
    {"qs_max", required_argument, NULL, 33},
    {"qi_min", required_argument, NULL, 34},
    {"qi_max", required_argument, NULL, 35},
    {"qh_min", required_argument, NULL, 36},
    {"qh_max", required_argument, NULL, 37},
    {"qg_min", required_argument, NULL, 38},
    {"qg_max", required_argument, NULL, 39},
    {"n_inact_min", required_argument, NULL, 40},
    {"n_inact_max", required_argument, NULL, 41},
    {"x_i_min", required_argument, NULL, 42},
    {"x_i_max", required_argument, NULL, 43},
    {"D_i_min", required_argument, NULL, 44},
    {"D_i_max", required_argument, NULL, 45},
    {"dep_snow_min", required_argument, NULL, 46},
    {"dep_snow_max", required_argument, NULL, 47},
    {"dep_ice_min", required_argument, NULL, 48},
    {"dep_ice_max", required_argument, NULL, 49},
    {"rime_qc_min", required_argument, NULL, 50},
    {"rime_qr_min", required_argument, NULL, 51},
    {"rime_qi_min", required_argument, NULL, 52},
    {"rime_qs_min", required_argument, NULL, 53},
    {"rime_qc_max", required_argument, NULL, 54},
    {"rime_qr_max", required_argument, NULL, 55},
    {"rime_qi_max", required_argument, NULL, 56},
    {"rime_qs_max", required_argument, NULL, 57},

    {"n1", required_argument, NULL, 58},
    {"n2", required_argument, NULL, 59},
    {"n3", required_argument, NULL, 60},

    {"function", required_argument, NULL, 61},

    {NULL, 0, NULL, 0}
};



int main(int argc, char** argv)
{

    uint32_t n1 = 500;
    uint32_t n2 = 500;
    uint32_t n3 = 1; // Default: No third iteration
    const double NOT_USED = -9999999;
    std::string func_name = "";

    codi::RealReverse p_prime = NOT_USED;
    codi::RealReverse T_prime = NOT_USED;
    codi::RealReverse ascent = NOT_USED;
    codi::RealReverse qv_prime = NOT_USED;
    codi::RealReverse qc_prime = NOT_USED;
    codi::RealReverse qr_prime = NOT_USED;
    codi::RealReverse qs_prime = NOT_USED;
    codi::RealReverse qi_prime = NOT_USED;
    codi::RealReverse qh_prime = NOT_USED;
    codi::RealReverse qg_prime = NOT_USED;
    codi::RealReverse temp_crit = NOT_USED;
    codi::RealReverse n_inact = NOT_USED;
    codi::RealReverse x_i = NOT_USED;
    codi::RealReverse D_i = NOT_USED;
    codi::RealReverse dep_snow = NOT_USED;
    codi::RealReverse dep_ice = NOT_USED;
    codi::RealReverse rime_qc = NOT_USED;
    codi::RealReverse rime_qr = NOT_USED;
    codi::RealReverse rime_qi = NOT_USED;
    codi::RealReverse rime_qs = NOT_USED;
    codi::RealReverse p_min = NOT_USED;
    codi::RealReverse p_max = NOT_USED;
    codi::RealReverse temp_min = NOT_USED;
    codi::RealReverse temp_max = NOT_USED;
    codi::RealReverse ascent_min = NOT_USED;
    codi::RealReverse ascent_max = NOT_USED;
    codi::RealReverse qv_min = NOT_USED;
    codi::RealReverse qv_max = NOT_USED;
    codi::RealReverse qc_min = NOT_USED;
    codi::RealReverse qc_max = NOT_USED;
    codi::RealReverse qr_min = NOT_USED;
    codi::RealReverse qr_max = NOT_USED;
    codi::RealReverse qs_min = NOT_USED;
    codi::RealReverse qs_max = NOT_USED;
    codi::RealReverse qi_min = NOT_USED;
    codi::RealReverse qi_max = NOT_USED;
    codi::RealReverse qh_min = NOT_USED;
    codi::RealReverse qh_max = NOT_USED;
    codi::RealReverse qg_min = NOT_USED;
    codi::RealReverse qg_max = NOT_USED;
    codi::RealReverse n_inact_min = NOT_USED;
    codi::RealReverse n_inact_max = NOT_USED;
    codi::RealReverse x_i_min = NOT_USED;
    codi::RealReverse x_i_max = NOT_USED;
    codi::RealReverse D_i_min = NOT_USED;
    codi::RealReverse D_i_max = NOT_USED;
    codi::RealReverse dep_snow_min = NOT_USED;
    codi::RealReverse dep_snow_max = NOT_USED;
    codi::RealReverse dep_ice_min = NOT_USED;
    codi::RealReverse dep_ice_max = NOT_USED;
    codi::RealReverse rime_qc_min = NOT_USED;
    codi::RealReverse rime_qr_min = NOT_USED;
    codi::RealReverse rime_qi_min = NOT_USED;
    codi::RealReverse rime_qs_min = NOT_USED;
    codi::RealReverse rime_qc_max = NOT_USED;
    codi::RealReverse rime_qr_max = NOT_USED;
    codi::RealReverse rime_qi_max = NOT_USED;
    codi::RealReverse rime_qs_max = NOT_USED;

    int ch;
    // loop over all of the options
    while ((ch = getopt_long(argc, argv, "", long_options, NULL)) != -1)
    {
        switch (ch)
        {
             case 0:
                p_prime = std::stod(optarg);
                break;
            case 1:
                T_prime = std::stod(optarg);
                break;
            case 2:
                ascent = std::stod(optarg);
                break;
            case 3:
                qv_prime = std::stod(optarg);
                break;
            case 4:
                qc_prime = std::stod(optarg);
                break;
            case 5:
                qr_prime = std::stod(optarg);
                break;
            case 6:
                qs_prime = std::stod(optarg);
                break;
            case 7:
                qi_prime = std::stod(optarg);
                break;
            case 8:
                qh_prime = std::stod(optarg);
                break;
            case 9:
                qg_prime = std::stod(optarg);
                break;
            case 10:
                temp_crit = std::stod(optarg);
                break;
            case 11:
                n_inact = std::stod(optarg);
                break;
            case 12:
                x_i = std::stod(optarg);
                break;
            case 13:
                D_i = std::stod(optarg);
                break;
            case 14:
                dep_snow = std::stod(optarg);
                break;
            case 15:
                dep_ice = std::stod(optarg);
                break;
            case 16:
                rime_qc = std::stod(optarg);
                break;
            case 17:
                rime_qr = std::stod(optarg);
                break;
            case 18:
                rime_qi = std::stod(optarg);
                break;
            case 19:
                rime_qs = std::stod(optarg);
                break;
            case 20:
                p_min = std::stod(optarg);
                break;
            case 21:
                p_max = std::stod(optarg);
                break;
            case 22:
                temp_min = std::stod(optarg);
                break;
            case 23:
                temp_max = std::stod(optarg);
                break;
            case 24:
                ascent_min = std::stod(optarg);
                break;
            case 25:
                ascent_max = std::stod(optarg);
                break;
            case 26:
                qv_min = std::stod(optarg);
                break;
            case 27:
                qv_max = std::stod(optarg);
                break;
            case 28:
                qc_min = std::stod(optarg);
                break;
            case 29:
                qc_max = std::stod(optarg);
                break;
            case 30:
                qr_min = std::stod(optarg);
                break;
            case 31:
                qr_max = std::stod(optarg);
                break;
            case 32:
                qs_min = std::stod(optarg);
                break;
            case 33:
                qs_max = std::stod(optarg);
                break;
            case 34:
                qi_min = std::stod(optarg);
                break;
            case 35:
                qi_max = std::stod(optarg);
                break;
            case 36:
                qh_min = std::stod(optarg);
                break;
            case 37:
                qh_max = std::stod(optarg);
                break;
            case 38:
                qg_min = std::stod(optarg);
                break;
            case 39:
                qg_max = std::stod(optarg);
                break;
            case 40:
                n_inact_min = std::stod(optarg);
                break;
            case 41:
                n_inact_max = std::stod(optarg);
                break;
            case 42:
                x_i_min = std::stod(optarg);
                break;
            case 43:
                x_i_max = std::stod(optarg);
                break;
            case 44:
                D_i_min = std::stod(optarg);
                break;
            case 45:
                D_i_max = std::stod(optarg);
                break;
            case 46:
                dep_snow_min = std::stod(optarg);
                break;
            case 47:
                dep_snow_max = std::stod(optarg);
                break;
            case 48:
                dep_ice_min = std::stod(optarg);
                break;
            case 49:
                dep_ice_max = std::stod(optarg);
                break;
            case 50:
                rime_qc_min = std::stod(optarg);
                break;
            case 51:
                rime_qr_min = std::stod(optarg);
                break;
            case 52:
                rime_qi_min = std::stod(optarg);
                break;
            case 53:
                rime_qs_min = std::stod(optarg);
                break;
            case 54:
                rime_qc_max = std::stod(optarg);
                break;
            case 55:
                rime_qr_max = std::stod(optarg);
                break;
            case 56:
                rime_qi_max = std::stod(optarg);
                break;
            case 57:
                rime_qs_max = std::stod(optarg);
                break;
            case 58:
                n1 = std::stoi(optarg);
                break;
            case 59:
                n2 = std::stoi(optarg);
                break;
            case 60:
                n3 = std::stoi(optarg);
                break;
            case 61:
                func_name = optarg;
                break;
            default:
                std::cout << "No such option " << ch << " with " << optarg << "\n";
                break;
        }
    }

    reference_quantities_t ref_quant;

    // Set standard values for input
    input_parameters_t input;
    init_input_parameters(input);
    auto_type = input.auto_type;
    load_lookup_table(ltabdminwgg);

    ref_quant.Tref = 273.15;
#ifdef WCB
    ref_quant.qref = 1.0e-6;
#else
    ref_quant.qref = 1.0e-4;
#endif
    ref_quant.pref = 1.0e5;
    ref_quant.wref = 1.; // 10.0
    ref_quant.tref = 1.0;
    ref_quant.zref = 1.0;

    codi::RealReverse qv = qv_prime/ref_quant.qref;
    codi::RealReverse qc = qc_prime/ref_quant.qref;
    codi::RealReverse qr = qr_prime/ref_quant.qref;
    codi::RealReverse qs = qs_prime/ref_quant.qref;
    codi::RealReverse qi = qi_prime/ref_quant.qref;
    codi::RealReverse qh = qh_prime/ref_quant.qref;
    codi::RealReverse qg = qg_prime/ref_quant.qref;

    model_constants_t cc;
    setup_model_constants(input, cc, ref_quant);

    std::vector<codi::RealReverse> y(num_comp);

    if(func_name.compare("ccn_act_seifert") == 0)
    {

    } else if(func_name.compare("ice_nuc_hom") == 0)
    {

    } else if(func_name.compare("ice_activation_hande") == 0)
    {

    } else if(func_name.compare("ice_activation_phillips") == 0)
    {

    } else if(func_name.compare("cloud_freeze_hom") == 0)
    {

    } else if(func_name.compare("ice_self_collection") == 0)
    {

    } else if(func_name.compare("snow_self_collection") == 0)
    {

    } else if(func_name.compare("graupel_melting") == 0)
    {

    } else if(func_name.compare("hail_melting") == 0)
    {

    } else if(func_name.compare("auto_conversion_kb") == 0)
    {

    } else if(func_name.compare("auto_conversion_sb") == 0)
    {

    } else if(func_name.compare("rain_self_collection_sb") == 0)
    {

    } else if(func_name.compare("rain_evaporation_sb") == 0)
    {

    } else if(func_name.compare("sedimentation_explicit") == 0)
    {

    } else if(func_name.compare("evaporation") == 0)
    {

    } else if(func_name.compare("vapor_dep_relaxation") == 0)
    {

    } else if(func_name.compare("particle_collection") == 0)
    {

    } else if(func_name.compare("particle_particle_collection") == 0)
    {

    } else if(func_name.compare("graupel_hail_conv") == 0)
    {

    } else if(func_name.compare("hail_collision") == 0)
    {

    } else if(func_name.compare("riming_cloud_core") == 0)
    {

    } else if(func_name.compare("riming_rain_core") == 0)
    {

    } else if(func_name.compare("ice_riming") == 0)
    {

    } else if(func_name.compare("snow_riming") == 0)
    {

    } else if(func_name.compare("particle_cloud_riming") == 0)
    {

    } else if(func_name.compare("particle_rain_riming") == 0)
    {

    } else if(func_name.compare("rain_freeze") == 0)
    {

    } else if(func_name.compare("ice_melting") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        std::cout << "qi,Ni,T,dt,delta_qi,delta_Ni,delta_qr,delta_Nr,"
                  << "delta_qc,delta_Nc,delta_lat_cool,delta_lat_heat\n";
        for(uint32_t i=0; i<n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;

            for(uint32_t j=0; j<n2; ++j)
            {
                if(qi_min != NOT_USED && qi_max != NOT_USED)
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;

                codi::RealReverse qi_in = qi_prime_in/ref_quant.qref;
                codi::RealReverse Ni = qi_in / (cc.ice.max_x - cc.ice.min_x_melt);
                std::cout << qi_prime_in.getValue() << "," << Ni.getValue() << ","
                          << T_prime_in.getValue() << "," << cc.dt_prime << ",";
                ice_melting(qi_prime_in, qi_in, Ni, T_prime_in,
                    cc.dt_prime, y, cc);
                std::cout << y[qi_idx].getValue() << "," << y[Ni_idx].getValue() << ","
                          << y[qr_idx].getValue() << "," << y[Nr_idx].getValue() << ","
                          << y[qc_idx].getValue() << "," << y[Nc_idx].getValue() << ","
                          << y[lat_cool_idx].getValue() << ","
                          << y[lat_heat_idx].getValue() << "\n";
            }
        }
    } else
    {
        std::cout << "No such method: " << func_name << "\n";
        return 1;
    }


//     if(std::strcmp("ice_table", argv[1]) == 0)
//     {
//         std::cout << "x,y,dust,soot,orga\n";
//         for(uint32_t i=0; i<n1; i++)
//         {
//             for(uint32_t j=0; j<n2; j++)
//             {
//                 // Not implemented yet
// //                 ice_table_scan(i, j, y);
//                 std::cout << i << "," << j << ","
//                           << y[0] << "," << y[1] << "," << y[2] << "\n";
//             }
//         }
//     }
//     else if(std::strcmp("ice_act_phillips", argv[1]) == 0)
//     {

//         bool use_prog_in = false;
//         std::cout << "inact,Sat_ice,T,qv,delta_ni,delta_qi\n";
//         for(uint32_t i=0; i<n1; i++)
//         {
//             codi::RealReverse qv = min_1 + (max_1-min_1)/n1 * i;
//             for(uint32_t j=0; j<n2; j++)
//             {
//                 codi::RealReverse T = min_2 + (max_2-min_2)/n2 * j;
//                 codi::RealReverse p_sat_ice = saturation_pressure_ice(T);
//                 codi::RealReverse ssi = qv * Rv * T / p_sat_ice;
//                 codi::RealReverse n_inact = 0;
//                 y[qi_idx] = y[Ni_idx] = 0.0;
//                 ice_activation_phillips(qc, qv, T, p_sat_ice, ssi,
//                     n_inact, use_prog_in, y, cc);
//                 std::cout << n_inact << ","
//                           << ssi << ","
//                           << T << ","
//                           << qv << ","
//                           << y[qi_idx].getValue() << ","
//                           << y[Ni_idx].getValue() << "\n";
//             }
//         }
//     }
//     else if(std::strcmp("infrac_scan", argv[1]) == 0)
//     {
//         std::cout << "x,y,T,dust,soot,orga\n";
//         for(uint32_t i=0; i<n1; i++)
//         {
//             codi::RealReverse T = min_1 + (max_1-min_1)/n1 * i;
//             // Not implemented yet
// //             infrac_table_scan(T, y);
//             std::cout << T.getValue() << "," << y[0].getValue()
//                       << "," << y[1].getValue() << "," << y[2].getValue() << "\n";
//         }
//     }
//     else
//     {
//         std::cout << "No such method: " << argv[1] << "\n";
//         return 1;
//     }
}