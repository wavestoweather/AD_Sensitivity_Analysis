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
    {"height", required_argument, NULL, 62},
    {"height_min", required_argument, NULL, 63},
    {"height_max", required_argument, NULL, 64},

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
    codi::RealReverse z_prime = NOT_USED;
    codi::RealReverse z_min = NOT_USED;
    codi::RealReverse z_max = NOT_USED;

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
            case 62:
                z_prime = std::stod(optarg);
                break;
            case 63:
                z_min = std::stod(optarg);
                break;
            case 64:
                z_max = std::stod(optarg);
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

    const double EPSILON = 1.0e-20;

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

    if(func_name.compare("ccn_act_hande") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse w_prime_in = ascent;
        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qc,Nc,T,S,qv,p,w,delta_qc,delta_Nc,delta_qv,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(ascent_min != NOT_USED && ascent_max != NOT_USED)
            {
                w_prime_in = i * (ascent_max-ascent_min) / n1 + ascent_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                {
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && ascent_min != NOT_USED && ascent_max != NOT_USED)
                {
                    w_prime_in = j * (ascent_max-ascent_min) / n2 + ascent_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && ascent_min != NOT_USED && ascent_max != NOT_USED)
                        w_prime_in = k * (ascent_max-ascent_min) / n3 + ascent_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);

                    std::cout << qc_prime_in.getValue() << ","
                              << Nc.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << w_prime_in.getValue() << ",";

                    ccn_act_hande(p_prime_in, w_prime_in, T_prime_in,
                        qv_prime_in, qc_prime_in, Nc, EPSILON, y, cc);

                    std::cout << y[qc_idx].getValue() << ","
                              << y[Nc_idx].getValue() << ","
                              << y[qv_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("ccn_act_seifert") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse w_prime_in = ascent;
        codi::RealReverse z_prime_in = z_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qc,Nc,T,S,qv,p,w,z,qr,Nr,delta_qc,delta_Nc,delta_qv,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(ascent_min != NOT_USED && ascent_max != NOT_USED)
            {
                w_prime_in = i * (ascent_max-ascent_min) / n1 + ascent_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 5;
            }
            else if(z_min != NOT_USED && z_max != NOT_USED)
            {
                z_prime_in = i * (z_max-z_min) / n1 + z_min;
                used_parameter = 6;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                {
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && ascent_min != NOT_USED && ascent_max != NOT_USED)
                {
                    w_prime_in = j * (ascent_max-ascent_min) / n2 + ascent_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }
                else if(used_parameter < 5 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 5;
                }
                else if(used_parameter < 6 && z_min != NOT_USED && z_max != NOT_USED)
                {
                    z_prime_in = j * (z_max-z_min) / n2 + z_min;
                    used_parameter2 = 6;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && ascent_min != NOT_USED && ascent_max != NOT_USED)
                        w_prime_in = k * (ascent_max-ascent_min) / n3 + ascent_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;
                    else if(used_parameter2 < 5 && qr_min != NOT_USED && qr_max != NOT_USED)
                        qr_prime_in = k * (qr_max-qr_min) / n3 + qr_min;
                    else if(used_parameter2 < 6 && z_min != NOT_USED && z_max != NOT_USED)
                        z_prime_in = k * (z_max-z_min) / n3 + z_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse Nr = qr_prime_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse p = p_prime_in / ref_quant.pref;
                    codi::RealReverse T = T_prime_in / ref_quant.Tref;
                    codi::RealReverse qv = qv_prime_in / ref_quant.qref;
                    codi::RealReverse qc = qc_prime_in / ref_quant.qref;
                    codi::RealReverse qr = qr_prime_in / ref_quant.qref;
                    codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime_in);

                    std::cout << qc_prime_in.getValue() << ","
                              << Nc.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << w_prime_in.getValue() << ","
                              << z_prime_in.getValue() << ","
                              << qr_prime_in.getValue() << ","
                              << Nr.getValue() << ",";

                    ccn_act_seifert(p_prime_in, p, T_prime_in, T,
                        qv_prime_in, qv, qc_prime_in, qc, Nc, qr,
                        z_prime_in, cc.dt_prime, w_prime_in, S, p_sat,
                        ref_quant, y, cc);

                    std::cout << y[qc_idx].getValue() << ","
                              << y[Nc_idx].getValue() << ","
                              << y[qv_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("ice_nuc_hom") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse w_prime_in = ascent;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qi,Ni,T,S,S_i,qv,p,w,delta_qi,delta_Ni,delta_qv,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(ascent_min != NOT_USED && ascent_max != NOT_USED)
            {
                w_prime_in = i * (ascent_max-ascent_min) / n1 + ascent_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qi_min != NOT_USED && qi_max != NOT_USED)
                {
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && ascent_min != NOT_USED && ascent_max != NOT_USED)
                {
                    w_prime_in = j * (ascent_max-ascent_min) / n2 + ascent_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && ascent_min != NOT_USED && ascent_max != NOT_USED)
                        w_prime_in = k * (ascent_max-ascent_min) / n3 + ascent_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Ni = qi_prime_in
                        / ( (cc.ice.max_x - cc.ice.min_x)/2 + cc.ice.min_x );

                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);

                    codi::RealReverse p_sat_ice = saturation_pressure_ice(T_prime_in);
                    codi::RealReverse ssi = qv_prime_in * Rv * T_prime_in / p_sat_ice;

                    std::cout << qi_prime_in.getValue() << ","
                              << Ni.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << ssi.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << w_prime_in.getValue() << ",";

                    ice_nuc_hom(T_prime_in, w_prime_in, p_prime_in,
                        qv_prime_in, qi_prime_in, Ni, ssi, p_sat_ice,
                        y, cc);


                    std::cout << y[qi_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[qv_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("ice_activation_hande") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse n_inact_in = n_inact;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qc,Nc,T,S,S_si,qv,delta_qi,delta_Ni,delta_qv,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(n_inact_min != NOT_USED && n_inact_max != NOT_USED)
            {
                n_inact_in = i * (n_inact_max-n_inact_min) / n1 + n_inact_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                {
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && n_inact_min != NOT_USED && n_inact_max != NOT_USED)
                {
                    n_inact_in = j * (n_inact_max-n_inact_min) / n2 + n_inact_min;
                    used_parameter = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && n_inact_min != NOT_USED && n_inact_max != NOT_USED)
                        n_inact_in = k * (n_inact_max-n_inact_min) / n3 + n_inact_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);

                    codi::RealReverse p_sat_ice = saturation_pressure_ice(T_prime_in);
                    codi::RealReverse ssi = qv_prime_in * Rv * T_prime_in / p_sat_ice;
                    bool ndiag_mask = false;

                    codi::RealReverse n_inact_tmp = n_inact_in;

                    std::cout << qc_prime_in.getValue() << ","
                              << Nc.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << ssi.getValue() << ","
                              << qv_prime_in.getValue() << ",";

                    ice_activation_hande(qc_prime_in, qv_prime_in, T_prime_in,
                        ssi, n_inact_in, ndiag_mask, y, cc);

                    std::cout << y[qi_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[qv_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                    n_inact_in = n_inact_tmp;
                }
            }
        }
    } else if(func_name.compare("ice_activation_phillips") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse n_inact_in = n_inact;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qc,Nc,T,S,qv,delta_qi,delta_Ni,delta_qv,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(n_inact_min != NOT_USED && n_inact_max != NOT_USED)
            {
                n_inact_in = i * (n_inact_max-n_inact_min) / n1 + n_inact_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                {
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && n_inact_min != NOT_USED && n_inact_max != NOT_USED)
                {
                    n_inact_in = j * (n_inact_max-n_inact_min) / n2 + n_inact_min;
                    used_parameter = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && n_inact_min != NOT_USED && n_inact_max != NOT_USED)
                        n_inact_in = k * (n_inact_max-n_inact_min) / n3 + n_inact_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);

                    codi::RealReverse p_sat_ice = saturation_pressure_ice(T_prime_in);
                    codi::RealReverse ssi = qv_prime_in * Rv * T_prime_in / p_sat_ice;
                    bool use_prog_in = false;

                    std::cout << qc_prime_in.getValue() << ","
                              << Nc.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << qv_prime_in.getValue() << ",";

                    ice_activation_phillips(qc_prime_in, qv_prime_in, T_prime_in,
                        p_sat_ice, ssi, n_inact_in, use_prog_in, y, cc);

                    std::cout << y[qi_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[qv_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("cloud_freeze_hom") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qc_prime_in = qc_prime;

        uint32_t used_parameter = 0;

        std::cout << "qc,Nc,T,delta_qi,delta_Ni,delta_qc,delta_Nc,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;

                for(auto& val: y)
                    val = 0;

                codi::RealReverse Nc = qc_prime_in
                    / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );

                std::cout << qc_prime_in.getValue() << ","
                          << Nc.getValue() << ","
                          << T_prime_in.getValue() << ",";

                codi::RealReverse T_c = T_prime_in - tmelt;

                cloud_freeze_hom(qc_prime_in, Nc, T_prime_in, T_c, y, cc);

                std::cout << y[qi_idx].getValue() << ","
                            << y[Ni_idx].getValue() << ","
                            << y[qc_idx].getValue() << ","
                            << y[Nc_idx].getValue() << ","
                            << y[lat_cool_idx].getValue() << ","
                            << y[lat_heat_idx].getValue() << "\n";
            }
        }
    } else if(func_name.compare("ice_self_collection") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qi,Ni,T,qv,p,S,delta_qi,delta_Ni,delta_qs,delta_Ns,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qi_min != NOT_USED && qi_max != NOT_USED)
                {
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Ni = qi_prime_in
                        / ( (cc.ice.max_x - cc.ice.min_x)/2 + cc.ice.min_x );

                    codi::RealReverse T_c = T_prime_in - tmelt;
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                            / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in,
                        T_prime_in, S)/rho_0);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);

                    std::cout << qi_prime_in.getValue() << ","
                            << Ni.getValue() << ","
                            << T_prime_in.getValue() << ","
                            << qv_prime_in.getValue() << ","
                            << p_prime_in.getValue() << ","
                            << S.getValue() << ",";

                    ice_self_collection(qi_prime_in, Ni, T_c, y, cc);

                    std::cout << y[qi_idx].getValue() << ","
                                << y[Ni_idx].getValue() << ","
                                << y[qs_idx].getValue() << ","
                                << y[Ns_idx].getValue() << ","
                                << y[lat_cool_idx].getValue() << ","
                                << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("snow_self_collection") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qs_prime_in = qs_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qs,Ns,T,qv,p,S,delta_Ns\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qs_min != NOT_USED && qs_max != NOT_USED)
            {
                qs_prime_in = i * (qs_max-qs_min) / n1 + qs_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qs_min != NOT_USED && qs_max != NOT_USED)
                {
                    qs_prime_in = j * (qs_max-qs_min) / n2 + qs_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Ns = qs_prime_in
                        / ( (cc.snow.max_x - cc.snow.min_x)/2 + cc.snow.min_x );

                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                            / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in,
                        T_prime_in, S)/rho_0);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);

                    std::cout << qs_prime_in.getValue() << ","
                            << Ns.getValue() << ","
                            << T_prime_in.getValue() << ","
                            << qv_prime_in.getValue() << ","
                            << p_prime_in.getValue() << ","
                            << S.getValue() << ",";

                    snow_self_collection(qs_prime_in, Ns, T_prime, y, cc);

                    std::cout << y[Ns_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("graupel_melting") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qg_prime_in = qg_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qg,Ng,T,qv,p,S,delta_qg,delta_Ng,"
                  << "delta_qr,delta_Nr,delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qg_min != NOT_USED && qg_max != NOT_USED)
            {
                qg_prime_in = i * (qg_max-qg_min) / n1 + qg_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qg_min != NOT_USED && qg_max != NOT_USED)
                {
                    qg_prime_in = j * (qg_max-qg_min) / n2 + qg_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Ng = qg_prime_in
                        / ( (cc.graupel.max_x - cc.graupel.min_x)/2 + cc.graupel.min_x );

                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                            / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in,
                        T_prime_in, S)/rho_0);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);

                    std::cout << qg_prime_in.getValue() << ","
                            << Ng.getValue() << ","
                            << T_prime_in.getValue() << ","
                            << qv_prime_in.getValue() << ","
                            << p_prime_in.getValue() << ","
                            << S.getValue() << ",";

                    graupel_melting(qg_prime_in, Ng, T_prime_in, y, cc);

                    std::cout << y[qg_idx].getValue() << ","
                              << y[Ng_idx].getValue() << ","
                              << y[qr_idx].getValue() << ","
                              << y[Nr_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("hail_melting") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qh_prime_in = qh_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qh,Nh,T,qv,p,S,delta_qh,delta_Nh,"
                  << "delta_qr,delta_Nr,delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qh_min != NOT_USED && qh_max != NOT_USED)
            {
                qh_prime_in = i * (qh_max-qh_min) / n1 + qh_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qh_min != NOT_USED && qh_max != NOT_USED)
                {
                    qh_prime_in = j * (qh_max-qh_min) / n2 + qh_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Nh = qh_prime_in
                        / ( (cc.hail.max_x - cc.hail.min_x)/2 + cc.hail.min_x );

                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                            / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in,
                        T_prime_in, S)/rho_0);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);

                    std::cout << qh_prime_in.getValue() << ","
                            << Nh.getValue() << ","
                            << T_prime_in.getValue() << ","
                            << qv_prime_in.getValue() << ","
                            << p_prime_in.getValue() << ","
                            << S.getValue() << ",";

                    hail_melting(qh_prime_in, Nh, T_prime_in, y, cc);

                    std::cout << y[qh_idx].getValue() << ","
                              << y[Nh_idx].getValue() << ","
                              << y[qr_idx].getValue() << ","
                              << y[Nr_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("auto_conversion_kb") == 0)
    {
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        std::cout << "qc,Nc,qr,delta_qr,delta_Nr,delta_qc,delta_Nc\n";

        uint32_t used_parameter = 0;

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(qc_min != NOT_USED && qc_max != NOT_USED)
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
            else if (qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 1;
            }
            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qr_min != NOT_USED && qr_max != NOT_USED)
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;

                for(auto& val: y)
                    val = 0;

                codi::RealReverse Nc = qc_prime_in
                    / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                std::cout << qc_prime_in.getValue() << ","
                          << Nc.getValue() << ","
                          << qr_prime_in.getValue() << ",";

                auto_conversion_kb(qc_prime_in, Nc, qr_prime_in, y, cc);

                std::cout << y[qr_idx].getValue() << ","
                          << y[Nr_idx].getValue() << ","
                          << y[qc_idx].getValue() << ","
                          << y[Nc_idx].getValue() << "\n";
            }
        }
    } else if(func_name.compare("auto_conversion_sb") == 0)
    {
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;

        std::cout << "qc,Nc,qr,T,p,qv,delta_qr,delta_Nr,delta_qc,delta_Nc\n";

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }
            else if (qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 2;
            }
             else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                {
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qr_min != NOT_USED && qr_max != NOT_USED)
                        qr_prime_in = k * (qr_max-qr_min) / n3 + qr_min;
                    else if(used_parameter2 < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                                / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);

                    std::cout << qc_prime_in.getValue() << ","
                            << Nc.getValue() << ","
                            << qr_prime_in.getValue() << ","
                            << T_prime_in.getValue() << ","
                            << p_prime_in.getValue() << ","
                            << qv_prime_in.getValue() << ",";

                    auto_conversion_sb(qc_prime_in, Nc, qr_prime_in, y, cc);

                    std::cout << y[qr_idx].getValue() << ","
                            << y[Nr_idx].getValue() << ","
                            << y[qc_idx].getValue() << ","
                            << y[Nc_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("rain_self_collection_sb") == 0)
    {
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;

        std::cout << "qr,Nr,T,p,qv,delta_Nr\n";

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 1;
            }
            else if (qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Nr = qr_prime_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                                / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in,
                        T_prime_in, S)/rho_0);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);

                    std::cout << qr_prime_in.getValue() << ","
                              << Nr.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << qv_prime_in.getValue() << ",";

                    rain_self_collection_sb(qr_prime_in, Nr, y, cc);

                    std::cout << y[Nr_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("rain_evaporation_sb") == 0)
    {
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse p_prime_in = p_prime;

        std::cout << "qr,Nr,T,p,qv,s_sw,p_sat,delta_qv,delta_qr,delta_Nr,"
                  << "delta_lat_cool,delta_lat_heat\n";

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 1;
            }
            else if (qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Nr = qr_prime_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                                / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse s_sw = S - 1.0;
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in,
                        T_prime_in, S)/rho_0);
                    codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime_in);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);

                    std::cout << qr_prime_in.getValue() << ","
                              << Nr.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << s_sw.getValue() << ","
                              << p_sat.getValue() << ",";
                    // Needs to be smaller than q_crit in order to trigger
                    // the process without any further involvement of
                    // the overall process. Hence we set it to zero:
                    qc_prime = 0;
                    rain_evaporation_sb(qr_prime_in, Nr, qv_prime_in, qc_prime,
                        T_prime_in, p_prime_in, s_sw, p_sat, y, cc);

                    std::cout << y[qv_idx].getValue() << ","
                              << y[qr_idx].getValue() << ","
                              << y[Nr_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("sedimentation_explicit") == 0)
    {
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse qs_prime_in = qs_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        codi::RealReverse qh_prime_in = qh_prime;
        codi::RealReverse qg_prime_in = qg_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse T_prime_in = T_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "T,S,qv,qc,qr,Nr,qs,Ns,qi,Ni,qh,Nh,qg,Ng,p,"
                  << "delta_qr,delta_Nr,delta_qi,delta_Ni,"
                  << "delta_qs,delta_Ns,delta_qg,delta_Ng,"
                  << "delta_qh,delta_Nh,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(qs_min != NOT_USED && qs_max != NOT_USED)
                qs_prime_in = i * (qs_max-qs_min) / n1 + qs_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 2;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 3;
            }
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 4;
            }
            else if(qh_min != NOT_USED && qh_max != NOT_USED)
            {
                qh_prime_in = i * (qh_max-qh_min) / n1 + qh_min;
                used_parameter = 5;
            }
            else if(qg_min != NOT_USED && qg_max != NOT_USED)
            {
                qg_prime_in = i * (qg_max-qg_min) / n1 + qg_min;
                used_parameter = 6;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 7;
            }
            else if(temp_min != NOT_USED && temp_max != NOT_USED)
            {
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
                used_parameter = 8;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                {
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qi_min != NOT_USED && qi_max != NOT_USED)
                {
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 4;
                }
                else if(used_parameter < 5 && qh_min != NOT_USED && qh_max != NOT_USED)
                {
                    qh_prime_in = j * (qh_max-qh_min) / n2 + qh_min;
                    used_parameter2 = 5;
                }
                else if(used_parameter < 6 && qg_min != NOT_USED && qg_max != NOT_USED)
                {
                    qg_prime_in = j * (qg_max-qg_min) / n2 + qg_min;
                    used_parameter2 = 6;
                }
                else if(used_parameter < 7 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 7;
                }
                else if(used_parameter < 8 && temp_min != NOT_USED && temp_max != NOT_USED)
                {
                    T_prime_in = j * (temp_max-temp_min) / n2 + temp_min;
                    used_parameter2 = 8;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qi_min != NOT_USED && qi_max != NOT_USED)
                        qi_prime_in = k * (qi_max-qi_min) / n3 + qi_min;
                    else if(used_parameter2 < 3 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;
                    else if(used_parameter2 < 4 && qr_min != NOT_USED && qr_max != NOT_USED)
                        qr_prime_in = k * (qr_max-qr_min) / n3 + qr_min;
                    else if(used_parameter2 < 5 && qh_min != NOT_USED && qh_max != NOT_USED)
                        qh_prime_in = k * (qh_max-qh_min) / n3 + qh_min;
                    else if(used_parameter2 < 6 && qg_min != NOT_USED && qg_max != NOT_USED)
                        qg_prime_in = k * (qg_max-qg_min) / n3 + qg_min;
                    else if(used_parameter2 < 7 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    if(used_parameter2 < 8 && temp_min != NOT_USED && temp_max != NOT_USED)
                        T_prime_in = k * (temp_max-temp_min) / n3 + temp_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                                / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse Nr = qr_prime_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                    codi::RealReverse Ns = qs_prime_in
                        / ( (cc.snow.max_x - cc.snow.min_x)/2 + cc.snow.min_x );
                    codi::RealReverse Ni = qi_prime_in
                        / ( (cc.ice.max_x - cc.ice.min_x)/2 + cc.ice.min_x );
                    codi::RealReverse Nh = qh_prime_in
                        / ( (cc.hail.max_x - cc.hail.min_x)/2 + cc.hail.min_x );
                    codi::RealReverse Ng = qg_prime_in
                        / ( (cc.graupel.max_x - cc.graupel.min_x)/2 + cc.graupel.min_x );

                    std::cout << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << qc_prime_in.getValue() << ","
                              << qr_prime_in.getValue() << ","
                              << Nr.getValue() << ","
                              << qs_prime_in.getValue() << ","
                              << Ns.getValue() << ","
                              << qi_prime_in.getValue() << ","
                              << Ni.getValue() << ","
                              << qh_prime_in.getValue() << ","
                              << Nh.getValue() << ","
                              << qg_prime_in.getValue() << ","
                              << Ng.getValue() << ","
                              << p_prime_in.getValue() << ",";

                    cc.rho_a_prime = compute_rhoa(p_prime_in,
                        T_prime_in, S);

                    sedimentation_explicit(T_prime_in, S,
                        qc_prime_in, qr_prime_in, Nr,
                        qs_prime_in, Ns, qi_prime_in, Ni, qh_prime_in, Nh,
                        qg_prime_in, Ng, p_prime_in, y, cc);

                    std::cout << y[qr_idx].getValue() << ","
                              << y[Nr_idx].getValue() << ","
                              << y[qi_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[qs_idx].getValue() << ","
                              << y[Ns_idx].getValue() << ","
                              << y[qg_idx].getValue() << ","
                              << y[Ng_idx].getValue() << ","
                              << y[qh_idx].getValue() << ","
                              << y[Nh_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("evaporation") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse q2_prime_in, q2_min, q2_max;
        codi::RealReverse avg_size;
        uint32_t q2_idx, N2_idx;
        char q2 = 'f';
        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;
        collection_model_constants_t *coeffs;
        particle_model_constants_t *pc;
        if(qg_prime != NOT_USED || qg_min != NOT_USED)
        {
            q2 = 'g';
            q2_min = qg_min;
            q2_max = qg_max;
            q2_prime_in = qg_prime;
            avg_size = (cc.graupel.max_x - cc.graupel.min_x) / 2 + cc.graupel.min_x;
            q2_idx = qg_idx;
            N2_idx = Ng_idx;
            coeffs = &cc.coeffs_gcr;
            pc = &cc.graupel;
        }
        else if(qs_prime != NOT_USED || qs_min != NOT_USED)
        {
            q2 = 's';
            q2_min = qs_min;
            q2_max = qs_max;
            q2_prime_in = qs_prime;
            avg_size = (cc.snow.max_x - cc.snow.min_x) / 2 + cc.snow.min_x;
            q2_idx = qs_idx;
            N2_idx = Ns_idx;
            pc = &cc.snow;
        }
        else if(qi_prime != NOT_USED || qi_min != NOT_USED)
        {
            q2 = 'i';
            q2_min = qi_min;
            q2_max = qi_max;
            q2_prime_in = qi_prime;
            avg_size = (cc.ice.max_x - cc.ice.min_x) / 2 + cc.ice.min_x;
            q2_idx = qi_idx;
            N2_idx = Ni_idx;
            pc = &cc.ice;
        }

        if(q2 == 'f')
        {
            std::cout << "evaporation needs either snow "
                     <<  "or graupel or ice particles\nABORTING!\n";
            return 1;
        }
        std::cout << "T,S,p,q" << q2 << ",N" << q2 << ","
                  << "qv,delta_qv,delta_q" << q2 << ","
                  << "delta_lat_cool,delta_lat_heat\n";
        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(q2_min != NOT_USED && q2_max != NOT_USED)
            {
                q2_prime_in = i * (q2_max-q2_min) / n1 + q2_min;
                used_parameter = 1;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 2;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 3;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && q2_min != NOT_USED && q2_max != NOT_USED)
                {
                    q2_prime_in = j * (q2_max-q2_min) / n2 + q2_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 3;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 3 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse N2 = q2_prime_in / avg_size;
                    // Update vertical velocity parameters
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse s_sw = S - 1.0;
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
                    cc.hail.rho_v = exp(-rho_vel * rho_inter);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);
                    codi::RealReverse e_d = qv_prime_in * Rv * T_prime_in;
                    codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime_in);

                    std::cout << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << q2_prime_in.getValue() << ","
                              << N2.getValue() << ","
                              << qv_prime_in.getValue() << ",";

                    evaporation(qv_prime_in, e_d, p_sat, s_sw, T_prime_in,
                        q2_prime_in, N2, y[q2_idx], *pc, y);

                    std::cout << y[qv_idx].getValue() << ","
                              << y[q2_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("vapor_dep_relaxation") == 0)
    {
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse qs_prime_in = qs_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        codi::RealReverse qh_prime_in = qh_prime;
        codi::RealReverse qg_prime_in = qg_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse T_prime_in = T_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qv,qi,Ni,qs,Ns,qg,Ng,qh,Nh,S_si,p_sat_ice,T,p,"
                  << "delta_qi,delta_qs,delta_qg,delta_qh,delta_qv,"
                  << "delta_deposition,delta_sublimination,"
                  << "delta_dep_rate_ice,delta_dep_rate_snow,"
                  << "delta_lat_cool,delta_lat_heat\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(qs_min != NOT_USED && qs_max != NOT_USED)
                qs_prime_in = i * (qs_max-qs_min) / n1 + qs_min;
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 1;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 2;
            }
            else if(qh_min != NOT_USED && qh_max != NOT_USED)
            {
                qh_prime_in = i * (qh_max-qh_min) / n1 + qh_min;
                used_parameter = 3;
            }
            else if(qg_min != NOT_USED && qg_max != NOT_USED)
            {
                qg_prime_in = i * (qg_max-qg_min) / n1 + qg_min;
                used_parameter = 4;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 5;
            }
            else if(temp_min != NOT_USED && temp_max != NOT_USED)
            {
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
                used_parameter = 6;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qi_min != NOT_USED && qi_max != NOT_USED)
                {
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qh_min != NOT_USED && qh_max != NOT_USED)
                {
                    qh_prime_in = j * (qh_max-qh_min) / n2 + qh_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && qg_min != NOT_USED && qg_max != NOT_USED)
                {
                    qg_prime_in = j * (qg_max-qg_min) / n2 + qg_min;
                    used_parameter2 = 4;
                }
                else if(used_parameter < 5 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 5;
                }
                else if(used_parameter < 6 && temp_min != NOT_USED && temp_max != NOT_USED)
                {
                    T_prime_in = j * (temp_max-temp_min) / n2 + temp_min;
                    used_parameter2 = 6;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;
                    else if(used_parameter2 < 3 && qh_min != NOT_USED && qh_max != NOT_USED)
                        qh_prime_in = k * (qh_max-qh_min) / n3 + qh_min;
                    else if(used_parameter2 < 4 && qg_min != NOT_USED && qg_max != NOT_USED)
                        qg_prime_in = k * (qg_max-qg_min) / n3 + qg_min;
                    else if(used_parameter2 < 5 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    if(used_parameter2 < 6 && temp_min != NOT_USED && temp_max != NOT_USED)
                        T_prime_in = k * (temp_max-temp_min) / n3 + temp_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Ns = qs_prime_in
                        / ( (cc.snow.max_x - cc.snow.min_x)/2 + cc.snow.min_x );
                    codi::RealReverse Ni = qi_prime_in
                        / ( (cc.ice.max_x - cc.ice.min_x)/2 + cc.ice.min_x );
                    codi::RealReverse Nh = qh_prime_in
                        / ( (cc.hail.max_x - cc.hail.min_x)/2 + cc.hail.min_x );
                    codi::RealReverse Ng = qg_prime_in
                        / ( (cc.graupel.max_x - cc.graupel.min_x)/2 + cc.graupel.min_x );
                    codi::RealReverse D_vtp = diffusivity(T_prime_in, p_prime_in);

                    codi::RealReverse p_sat_ice = saturation_pressure_ice(T_prime_in);
                    codi::RealReverse e_d = qv_prime_in * Rv * T_prime_in;
                    codi::RealReverse s_si = e_d / p_sat_ice - 1.0;

                    codi::RealReverse dep_rate_snow, dep_rate_ice;
                    dep_rate_snow = dep_rate_ice = 0;

                    std::cout << qv_prime_in.getValue() << ","
                              << qi_prime_in.getValue() << ","
                              << Ni.getValue() << ","
                              << qs_prime_in.getValue() << ","
                              << Ns.getValue() << ","
                              << qg_prime_in.getValue() << ","
                              << Ng.getValue() << ","
                              << qh_prime_in.getValue() << ","
                              << Nh.getValue() << ","
                              << s_si.getValue() << ","
                              << p_sat_ice.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << p_prime_in.getValue() << ",";

                    vapor_dep_relaxation(qv_prime_in,
                        qi_prime_in, Ni, qs_prime_in, Ns,
                        qg_prime_in, Ng, qh_prime_in, Nh,
                        s_si, p_sat_ice, T_prime_in, EPSILON,
                        dep_rate_ice, dep_rate_snow,
                        D_vtp, y, cc);

                    std::cout << y[qi_idx].getValue() << ","
                              << y[qs_idx].getValue() << ","
                              << y[qg_idx].getValue() << ","
                              << y[qh_idx].getValue() << ","
                              << y[qv_idx].getValue() << ","
                              << y[depo_idx].getValue() << ","
                              << y[sub_idx].getValue() << ","
                              << dep_rate_ice.getValue() << ","
                              << dep_rate_snow.getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("particle_collection") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse q2_prime_in, q2_min, q2_max, q1_prime_in, q1_min, q1_max;
        codi::RealReverse avg_size1, avg_size2;
        uint32_t q2_idx, N2_idx, q1_idx, N1_idx;
        char q2 = 'f';
        char q1 = 'f';
        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;
        collection_model_constants_t *coeffs;
        particle_model_constants_t *pc1, *pc2;

        if(qi_prime != NOT_USED || qi_min != NOT_USED || qi_prime != NOT_USED)
        {
            q1 = 'i';
            q1_min = qi_min;
            q1_max = qi_max;
            q1_prime_in = qi_prime;
            avg_size1 = (cc.ice.max_x - cc.ice.min_x) / 2 + cc.ice.min_x;
            q1_idx = qi_idx;
            N1_idx = Ni_idx;
            pc1 = &cc.ice;

            if(qh_prime != NOT_USED || qh_min != NOT_USED || qh_prime != NOT_USED)
            {
                q2 = 'h';
                q2_min = qh_min;
                q2_max = qh_max;
                q2_prime_in = qh_prime;
                avg_size2 = (cc.hail.max_x - cc.hail.min_x) / 2 + cc.hail.min_x;
                q2_idx = qh_idx;
                N2_idx = Nh_idx;
                pc2 = &cc.hail;
                coeffs = &cc.coeffs_hic;
            }
            else if(qs_prime != NOT_USED || qs_min != NOT_USED || qs_prime != NOT_USED)
            {
                q2 = 's';
                q2_min = qs_min;
                q2_max = qs_max;
                q2_prime_in = qs_prime;
                avg_size2 = (cc.snow.max_x - cc.snow.min_x) / 2 + cc.snow.min_x;
                q2_idx = qs_idx;
                N2_idx = Ns_idx;
                pc2 = &cc.snow;
                coeffs = &cc.coeffs_sic;
            }
            else if(qg_prime != NOT_USED || qg_min != NOT_USED || qg_prime != NOT_USED)
            {
                q2 = 'g';
                q2_min = qg_min;
                q2_max = qg_max;
                q2_prime_in = qg_prime;
                avg_size2 = (cc.graupel.max_x - cc.graupel.min_x) / 2 + cc.graupel.min_x;
                q2_idx = qg_idx;
                N2_idx = Ng_idx;
                pc2 = &cc.graupel;
                coeffs = &cc.coeffs_gic;
            }
        } else if(qs_prime != NOT_USED || qs_min != NOT_USED || qs_prime != NOT_USED)
        {
            q1 = 's';
            q1_min = qs_min;
            q1_max = qs_max;
            q1_prime_in = qs_prime;
            avg_size1 = (cc.snow.max_x - cc.snow.min_x) / 2 + cc.snow.min_x;
            q1_idx = qs_idx;
            N1_idx = Ns_idx;
            pc1 = &cc.snow;

            if(qh_prime != NOT_USED || qh_min != NOT_USED || qh_prime != NOT_USED)
            {
                q2 = 'h';
                q2_min = qh_min;
                q2_max = qh_max;
                q2_prime_in = qh_prime;
                avg_size2 = (cc.hail.max_x - cc.hail.min_x) / 2 + cc.hail.min_x;
                q2_idx = qh_idx;
                N2_idx = Nh_idx;
                pc2 = &cc.hail;
                coeffs = &cc.coeffs_hsc;
            }
            else if(qg_prime != NOT_USED || qg_min != NOT_USED || qg_prime != NOT_USED)
            {
                q2 = 'g';
                q2_min = qg_min;
                q2_max = qg_max;
                q2_prime_in = qg_prime;
                avg_size2 = (cc.graupel.max_x - cc.graupel.min_x) / 2 + cc.graupel.min_x;
                q2_idx = qg_idx;
                N2_idx = Ng_idx;
                pc2 = &cc.graupel;
                coeffs = &cc.coeffs_gsc;
            }
        }

        if(q2 == 'f' || q1 == 'f')
        {
            std::cout << "particle_collection needs either ice with hail "
                      << ", snow or graupel or snow woth hail or graupel"
                      <<  "\nABORTING!\n";
            return 1;
        }

        std::cout << "p,T,S,qv,q" << q1 << ",N" << q1 << ",q" << q2 << ",N" << q2 << ","
                  << "delta_q" << q1 << ",delta_q" << q2 << ",delta_N"
                  << q1 << "\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 1;
            }
            else if(q1_min != NOT_USED && q1_max != NOT_USED)
            {
                q1_prime_in = i * (q1_max-q1_min) / n1 + q1_min;
                used_parameter = 2;
            }
            else if(q2_min != NOT_USED && q2_max != NOT_USED)
            {
                q2_prime_in = i * (q2_max-q2_min) / n1 + q2_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && q1_min != NOT_USED && q1_max != NOT_USED)
                {
                    q1_prime_in = j * (q1_max-q1_min) / n2 + q1_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && q2_min != NOT_USED && q2_max != NOT_USED)
                {
                    q2_prime_in = j * (q2_max-q2_min) / n2 + q2_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && q1_min != NOT_USED && q1_max != NOT_USED)
                        q1_prime_in = k * (q1_max-q1_min) / n3 + q1_min;
                    else if(used_parameter2 < 3 && q2_min != NOT_USED && q2_max != NOT_USED)
                        q2_prime_in = k * (q2_max-q2_min) / n3 + q2_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse N1 = q1_prime_in / avg_size1;
                    codi::RealReverse N2 = q2_prime_in / avg_size2;
                    // Update vertical velocity parameters
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
                    cc.hail.rho_v = exp(-rho_vel * rho_inter);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);
                    std::cout << p_prime_in.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << q1_prime_in.getValue() << ","
                              << N1.getValue() << ","
                              << q2_prime_in.getValue() << ","
                              << N2.getValue() << ",";
                    codi::RealReverse T_c = T_prime_in - tmelt;

                    std::vector<codi::RealReverse> res = particle_collection(
                        q1_prime_in, q2_prime_in, N1, N2, T_c,
                        *coeffs, *pc1, *pc2);

                    std::cout << -res[1].getValue() << ","
                              << res[1].getValue() << ","
                              << -res[0].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("particle_particle_collection") == 0)
    {
        codi::RealReverse qs_prime_in = qs_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        codi::RealReverse qg_prime_in = qg_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse T_prime_in = T_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qv,qi,Ni,qs,Ns,qg,Ng,S,T,p,"
                  << "delta_qi,delta_qs,delta_qg,"
                  << "delta_Ni,delta_Ns,delta_Ng\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(qs_min != NOT_USED && qs_max != NOT_USED)
                qs_prime_in = i * (qs_max-qs_min) / n1 + qs_min;
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 1;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 2;
            }
            else if(qg_min != NOT_USED && qg_max != NOT_USED)
            {
                qg_prime_in = i * (qg_max-qg_min) / n1 + qg_min;
                used_parameter = 3;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 4;
            }
            else if(temp_min != NOT_USED && temp_max != NOT_USED)
            {
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
                used_parameter = 5;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qi_min != NOT_USED && qi_max != NOT_USED)
                {
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qg_min != NOT_USED && qg_max != NOT_USED)
                {
                    qg_prime_in = j * (qg_max-qg_min) / n2 + qg_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 4;
                }
                else if(used_parameter < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                {
                    T_prime_in = j * (temp_max-temp_min) / n2 + temp_min;
                    used_parameter2 = 5;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;
                    else if(used_parameter2 < 3 && qg_min != NOT_USED && qg_max != NOT_USED)
                        qg_prime_in = k * (qg_max-qg_min) / n3 + qg_min;
                    else if(used_parameter2 < 4 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    if(used_parameter2 < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                        T_prime_in = k * (temp_max-temp_min) / n3 + temp_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Ns = qs_prime_in
                        / ( (cc.snow.max_x - cc.snow.min_x)/2 + cc.snow.min_x );
                    codi::RealReverse Ni = qi_prime_in
                        / ( (cc.ice.max_x - cc.ice.min_x)/2 + cc.ice.min_x );
                    codi::RealReverse Ng = qg_prime_in
                        / ( (cc.graupel.max_x - cc.graupel.min_x)/2 + cc.graupel.min_x );

                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
                    cc.hail.rho_v = exp(-rho_vel * rho_inter);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);

                    std::cout << qv_prime_in.getValue() << ","
                              << qi_prime_in.getValue() << ","
                              << Ni.getValue() << ","
                              << qs_prime_in.getValue() << ","
                              << Ns.getValue() << ","
                              << qg_prime_in.getValue() << ","
                              << Ng.getValue() << ","
                              << S.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << p_prime_in.getValue() << ",";
                    codi::RealReverse T_c = T_prime_in - tmelt;

                    particle_particle_collection(
                        qi_prime_in, Ni, qs_prime_in, Ns,
                        qg_prime_in, Ng, T_prime_in, T_c,
                        y, cc);

                    std::cout << y[qi_idx].getValue() << ","
                              << y[qs_idx].getValue() << ","
                              << y[qg_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[Ns_idx].getValue() << ","
                              << y[Ng_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("graupel_hail_conv") == 0)
    {
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        codi::RealReverse qg_prime_in = qg_prime;
        codi::RealReverse qh_prime_in = qh_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse T_prime_in = T_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qc,qr,qi,qg,Ng,qh,Nh,p,T,"
                  << "delta_qg,delta_Ng,delta_qh,delta_Nh\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(qc_min != NOT_USED && qc_max != NOT_USED)
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 1;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 2;
            }
            else if(qg_min != NOT_USED && qg_max != NOT_USED)
            {
                qg_prime_in = i * (qg_max-qg_min) / n1 + qg_min;
                used_parameter = 3;
            }
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 4;
            }
            else if(temp_min != NOT_USED && temp_max != NOT_USED)
            {
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
                used_parameter = 5;
            }
            else if(qh_min != NOT_USED && qh_max != NOT_USED)
            {
                qh_prime_in = i * (qh_max-qh_min) / n1 + qh_min;
                used_parameter = 6;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qi_min != NOT_USED && qi_max != NOT_USED)
                {
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qg_min != NOT_USED && qg_max != NOT_USED)
                {
                    qg_prime_in = j * (qg_max-qg_min) / n2 + qg_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 4;
                }
                else if(used_parameter < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                {
                    T_prime_in = j * (temp_max-temp_min) / n2 + temp_min;
                    used_parameter2 = 5;
                }
                else if(used_parameter < 6 && qh_min != NOT_USED && qh_max != NOT_USED)
                {
                    qh_prime_in = j * (qh_max-qh_min) / n2 + qh_min;
                    used_parameter2 = 6;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;
                    else if(used_parameter2 < 3 && qg_min != NOT_USED && qg_max != NOT_USED)
                        qg_prime_in = k * (qg_max-qg_min) / n3 + qg_min;
                    else if(used_parameter2 < 4 && qr_min != NOT_USED && qr_max != NOT_USED)
                        qr_prime_in = k * (qr_max-qr_min) / n3 + qr_min;
                    if(used_parameter2 < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                        T_prime_in = k * (temp_max-temp_min) / n3 + temp_min;
                    else if(used_parameter2 < 6 && qh_min != NOT_USED && qh_max != NOT_USED)
                        qh_prime_in = k * (qh_max-qh_min) / n3 + qh_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Nh = qh_prime_in
                        / ( (cc.hail.max_x - cc.hail.min_x)/2 + cc.hail.min_x );
                    codi::RealReverse Ng = qg_prime_in
                        / ( (cc.graupel.max_x - cc.graupel.min_x)/2 + cc.graupel.min_x );

                    std::cout << qc_prime_in.getValue() << ","
                              << qr_prime_in.getValue() << ","
                              << qi_prime_in.getValue() << ","
                              << qg_prime_in.getValue() << ","
                              << Ng.getValue() << ","
                              << qh_prime_in.getValue() << ","
                              << Nh.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << T_prime_in.getValue() << ",";
                    codi::RealReverse T_c = T_prime_in - tmelt;

                    graupel_hail_conv(qc_prime_in, qr_prime_in,
                        qi_prime_in, qg_prime_in, Ng,
                        qh_prime_in, Nh, p_prime_in, T_prime_in, T_c,
                        y, cc);

                    std::cout << y[qg_idx].getValue() << ","
                              << y[Ng_idx].getValue() << ","
                              << y[qh_idx].getValue() << ","
                              << y[Nh_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("hail_collision") == 0)
    {
        codi::RealReverse qs_prime_in = qs_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        codi::RealReverse qh_prime_in = qh_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse T_prime_in = T_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qv,qi,Ni,qs,Ns,qh,Nh,S,T,p,"
                  << "delta_qi,delta_qs,delta_qh,"
                  << "delta_Ni,delta_Ns,delta_Nh\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(qs_min != NOT_USED && qs_max != NOT_USED)
                qs_prime_in = i * (qs_max-qs_min) / n1 + qs_min;
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 1;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 2;
            }
            else if(qh_min != NOT_USED && qh_max != NOT_USED)
            {
                qh_prime_in = i * (qh_max-qh_min) / n1 + qh_min;
                used_parameter = 3;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 4;
            }
            else if(temp_min != NOT_USED && temp_max != NOT_USED)
            {
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
                used_parameter = 5;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qi_min != NOT_USED && qi_max != NOT_USED)
                {
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qh_min != NOT_USED && qh_max != NOT_USED)
                {
                    qh_prime_in = j * (qh_max-qh_min) / n2 + qh_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 4;
                }
                else if(used_parameter < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                {
                    T_prime_in = j * (temp_max-temp_min) / n2 + temp_min;
                    used_parameter2 = 5;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;
                    else if(used_parameter2 < 3 && qh_min != NOT_USED && qh_max != NOT_USED)
                        qh_prime_in = k * (qh_max-qh_min) / n3 + qh_min;
                    else if(used_parameter2 < 4 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    if(used_parameter2 < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                        T_prime_in = k * (temp_max-temp_min) / n3 + temp_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Ns = qs_prime_in
                        / ( (cc.snow.max_x - cc.snow.min_x)/2 + cc.snow.min_x );
                    codi::RealReverse Ni = qi_prime_in
                        / ( (cc.ice.max_x - cc.ice.min_x)/2 + cc.ice.min_x );
                    codi::RealReverse Nh = qh_prime_in
                        / ( (cc.hail.max_x - cc.hail.min_x)/2 + cc.hail.min_x );

                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
                    cc.hail.rho_v = exp(-rho_vel * rho_inter);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);

                    std::cout << qv_prime_in.getValue() << ","
                              << qi_prime_in.getValue() << ","
                              << Ni.getValue() << ","
                              << qs_prime_in.getValue() << ","
                              << Ns.getValue() << ","
                              << qh_prime_in.getValue() << ","
                              << Nh.getValue() << ","
                              << S.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << p_prime_in.getValue() << ",";
                    codi::RealReverse T_c = T_prime_in - tmelt;

                    hail_collision(qh_prime_in, Nh,
                        qs_prime_in, Ns, qi_prime_in, Ni, T_c, y, cc);

                    std::cout << y[qi_idx].getValue() << ","
                              << y[qs_idx].getValue() << ","
                              << y[qh_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[Ns_idx].getValue() << ","
                              << y[Nh_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("riming_cloud_core") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse q2_prime_in, q2_min, q2_max;
        codi::RealReverse avg_size;
        uint32_t q2_idx, N2_idx;
        char q2 = 'f';
        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;
        collection_model_constants_t *coeffs;
        particle_model_constants_t *pc;
        if(qs_prime != NOT_USED || qs_min != NOT_USED || qs_prime != NOT_USED)
        {
            q2 = 's';
            q2_min = qs_min;
            q2_max = qs_max;
            q2_prime_in = qs_prime;
            avg_size = (cc.snow.max_x - cc.snow.min_x) / 2 + cc.snow.min_x;
            q2_idx = qs_idx;
            N2_idx = Ns_idx;
            coeffs = &cc.coeffs_scr;
            pc = &cc.snow;
        }
        else if(qi_prime != NOT_USED || qi_min != NOT_USED || qi_prime != NOT_USED)
        {
            q2 = 'i';
            q2_min = qi_min;
            q2_max = qi_max;
            q2_prime_in = qi_prime;
            avg_size = (cc.ice.max_x - cc.ice.min_x) / 2 + cc.ice.min_x;
            q2_idx = qi_idx;
            N2_idx = Ni_idx;
            coeffs = &cc.coeffs_icr;
            pc = &cc.ice;
        }

        if(q2 == 'f')
        {
            std::cout << "riming_cloud_core needs either ice "
                     <<  "or snow particles\nABORTING!\n";
            return 1;
        }
        std::cout << "qc,Nc,T,S,p,qv,q" << q2 << ",N" << q2
                  << ",delta_rime_rate_q" << q2
                  << ",delta_rime_rate_N" << q2 << "\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }
            else if(q2_min != NOT_USED && q2_max != NOT_USED)
            {
                q2_prime_in = i * (q2_max-q2_min) / n1 + q2_min;
                used_parameter = 2;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                {
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && q2_min != NOT_USED && q2_max != NOT_USED)
                {
                    q2_prime_in = j * (q2_max-q2_min) / n2 + q2_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && q2_min != NOT_USED && q2_max != NOT_USED)
                        q2_prime_in = k * (q2_max-q2_min) / n3 + q2_min;
                    else if(used_parameter2 < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse N2 = q2_prime_in / avg_size;
                    // Update vertical velocity parameters
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
                    cc.hail.rho_v = exp(-rho_vel * rho_inter);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);
                    std::cout << qc_prime_in.getValue() << ","
                              << Nc.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << q2_prime_in.getValue() << ","
                              << N2.getValue() << ",";
                    codi::RealReverse rime_rate_qb, rime_rate_nb;
                    rime_rate_qb = rime_rate_nb = 0;

                    riming_cloud_core(qc_prime_in, Nc, q2_prime_in, N2,
                        *pc, *coeffs, rime_rate_qb, rime_rate_nb, cc);

                    std::cout << rime_rate_qb.getValue() << ","
                              << rime_rate_nb.getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("riming_rain_core") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse q2_prime_in, q2_min, q2_max;
        codi::RealReverse avg_size;
        uint32_t q2_idx, N2_idx;
        char q2 = 'f';
        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;
        collection_model_constants_t *coeffs;
        particle_model_constants_t *pc;
        if(qs_prime != NOT_USED || qs_min != NOT_USED || qs_prime != NOT_USED)
        {
            q2 = 's';
            q2_min = qs_min;
            q2_max = qs_max;
            q2_prime_in = qs_prime;
            avg_size = (cc.snow.max_x - cc.snow.min_x) / 2 + cc.snow.min_x;
            q2_idx = qs_idx;
            N2_idx = Ns_idx;
            coeffs = &cc.coeffs_scr;
            pc = &cc.snow;
        }
        else if(qi_prime != NOT_USED || qi_min != NOT_USED || qi_prime != NOT_USED)
        {
            q2 = 'i';
            q2_min = qi_min;
            q2_max = qi_max;
            q2_prime_in = qi_prime;
            avg_size = (cc.ice.max_x - cc.ice.min_x) / 2 + cc.ice.min_x;
            q2_idx = qi_idx;
            N2_idx = Ni_idx;
            coeffs = &cc.coeffs_icr;
            pc = &cc.ice;
        }

        if(q2 == 'f')
        {
            std::cout << "riming_rain_core needs either ice "
                     <<  "or snow particles\nABORTING!\n";
            return 1;
        }
        std::cout << "qr,Nr,T,S,p,qv,q" << q2 << ",N" << q2
                  << ",delta_rime_rate_q" << q2
                  << ",delta_rime_rate_qr,delta_rime_rate_Nr\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 1;
            }
            else if(q2_min != NOT_USED && q2_max != NOT_USED)
            {
                q2_prime_in = i * (q2_max-q2_min) / n1 + q2_min;
                used_parameter = 2;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && q2_min != NOT_USED && q2_max != NOT_USED)
                {
                    q2_prime_in = j * (q2_max-q2_min) / n2 + q2_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && q2_min != NOT_USED && q2_max != NOT_USED)
                        q2_prime_in = k * (q2_max-q2_min) / n3 + q2_min;
                    else if(used_parameter2 < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse Nr = qr_prime_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                    codi::RealReverse N2 = q2_prime_in / avg_size;
                    // Update vertical velocity parameters
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
                    cc.hail.rho_v = exp(-rho_vel * rho_inter);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);
                    std::cout << qr_prime_in.getValue() << ","
                              << Nr.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << q2_prime_in.getValue() << ","
                              << N2.getValue() << ",";

                    codi::RealReverse rime_rate_qa, rime_rate_qb, rime_rate_nb;
                    rime_rate_qa = rime_rate_qb = rime_rate_nb = 0;

                    riming_rain_core(qr_prime_in, Nr, q2_prime_in, N2,
                        *pc, *coeffs, rime_rate_qa, rime_rate_qb,
                        rime_rate_nb, cc);

                    std::cout << rime_rate_qa.getValue() << ","
                              << rime_rate_qb.getValue() << ","
                              << rime_rate_nb.getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("ice_riming") == 0)
    {
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        codi::RealReverse dep_rate_ice_in = dep_ice;
        codi::RealReverse rime_rate_qc_in = rime_qc;
        codi::RealReverse rime_rate_qr_in = rime_qr;
        codi::RealReverse rime_rate_qi_in = rime_qi;
        codi::RealReverse T_prime_in = T_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qc,Nc,qr,Nr,qi,Ni,dep_rate_ice,rime_rate_qc,"
                  << "rime_rate_nc,rime_rate_qr,rime_rate_nr,"
                  << "rime_rate_qi,T,dt,"
                  << "delta_qc,delta_Nc,delta_qi,delta_Ni,"
                  << "delta_qr,delta_Nr,delta_qs,delta_qg,delta_Ng,"
                  << "delta_lat_cool_idx,delta_lat_heat_idx\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(qc_min != NOT_USED && qc_max != NOT_USED)
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 1;
            }
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 2;
            }
            else if(dep_ice_min != NOT_USED && dep_ice_max != NOT_USED)
            {
                dep_rate_ice_in = i * (dep_ice_max-dep_ice_min) / n1 + dep_ice_min;
                used_parameter = 3;
            }
            else if(rime_qc_min != NOT_USED && rime_qc_max != NOT_USED)
            {
                rime_rate_qc_in = i * (rime_qc_max-rime_qc_min) / n1 + rime_qc_min;
                used_parameter = 4;
            }
            else if(temp_min != NOT_USED && temp_max != NOT_USED)
            {
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
                used_parameter = 5;
            }
            else if(rime_qr_min != NOT_USED && rime_qr_max != NOT_USED)
            {
                rime_rate_qr_in = i * (rime_qr_max-rime_qr_min) / n1 + rime_qr_min;
                used_parameter = 6;
            }
            else if(rime_qi_min != NOT_USED && rime_qi_max != NOT_USED)
            {
                rime_rate_qi_in = i * (rime_qi_max-rime_qi_min) / n1 + rime_qi_min;
                used_parameter = 7;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qi_min != NOT_USED && qi_max != NOT_USED)
                {
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && dep_ice_min != NOT_USED && dep_ice_max != NOT_USED)
                {
                    dep_rate_ice_in = j * (dep_ice_max-dep_ice_min) / n2 + dep_ice_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && rime_qc_min != NOT_USED && rime_qc_max != NOT_USED)
                {
                    rime_rate_qc_in = j * (rime_qc_max-rime_qc_min) / n2 + rime_qc_min;
                    used_parameter2 = 4;
                }
                else if(used_parameter < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                {
                    T_prime_in = j * (temp_max-temp_min) / n2 + temp_min;
                    used_parameter2 = 5;
                }
                else if(used_parameter < 6 && rime_qr_min != NOT_USED && rime_qr_max != NOT_USED)
                {
                    rime_rate_qr_in = j * (rime_qr_max-rime_qr_min) / n2 + rime_qr_min;
                    used_parameter2 = 6;
                }
                else if(used_parameter < 4 && rime_qi_min != NOT_USED && rime_qi_max != NOT_USED)
                {
                    rime_rate_qi_in = j * (rime_qi_max-rime_qi_min) / n2 + rime_qi_min;
                    used_parameter2 = 7;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qr_min != NOT_USED && qr_max != NOT_USED)
                        qr_prime_in = k * (qr_max-qr_min) / n3 + qr_min;
                    else if(used_parameter2 < 3 && dep_ice_min != NOT_USED && dep_ice_max != NOT_USED)
                        dep_rate_ice_in = k * (dep_ice_max-dep_ice_min) / n3 + dep_ice_min;
                    else if(used_parameter2 < 4 && rime_qc_min != NOT_USED && rime_qc_max != NOT_USED)
                        rime_rate_qc_in = k * (rime_qc_max-rime_qc_min) / n3 + rime_qc_min;
                    if(used_parameter2 < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                        T_prime_in = k * (temp_max-temp_min) / n3 + temp_min;
                    else if(used_parameter2 < 6 && rime_qr_min != NOT_USED && rime_qr_max != NOT_USED)
                        rime_rate_qr_in = k * (rime_qr_max-rime_qr_min) / n3 + rime_qr_min;
                    else if(used_parameter2 < 7 && rime_qi_min != NOT_USED && rime_qi_max != NOT_USED)
                        rime_rate_qi_in = k * (rime_qi_max-rime_qi_min) / n3 + rime_qi_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse Nr = qr_prime_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                    codi::RealReverse Ni = qi_prime_in
                        / ( (cc.ice.max_x - cc.ice.min_x)/2 + cc.ice.min_x );

                    codi::RealReverse rime_rate_nc = rime_rate_qc_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse rime_rate_nr = rime_rate_qr_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );

                    std::cout << qc_prime_in.getValue() << ","
                              << Nc.getValue() << ","
                              << qr_prime_in.getValue() << ","
                              << Nr.getValue() << ","
                              << qi_prime_in.getValue() << ","
                              << Ni.getValue() << ","
                              << dep_rate_ice_in.getValue() << ","
                              << rime_rate_qc_in.getValue() << ","
                              << rime_rate_nc.getValue() << ","
                              << rime_rate_qr_in.getValue() << ","
                              << rime_rate_nr.getValue() << ","
                              << rime_rate_qi_in.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << cc.dt_prime << ",";

                    ice_riming(qc_prime_in, Nc, qr_prime_in, Nr,
                        qi_prime_in, Ni, dep_rate_ice_in, rime_rate_qc_in,
                        rime_rate_nc, rime_rate_qr_in, rime_rate_nr,
                        rime_rate_qi_in, T_prime_in, cc.dt_prime, y, cc);

                    std::cout << y[qc_idx].getValue() << ","
                              << y[Nc_idx].getValue() << ","
                              << y[qi_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[qr_idx].getValue() << ","
                              << y[Nr_idx].getValue() << ","
                              << y[qs_idx].getValue() << ","
                              << y[qg_idx].getValue() << ","
                              << y[Ng_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("snow_riming") == 0)
    {
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse qs_prime_in = qs_prime;
        codi::RealReverse dep_rate_snow_in = dep_snow;
        codi::RealReverse rime_rate_qc_in = rime_qc;
        codi::RealReverse rime_rate_qr_in = rime_qr;
        codi::RealReverse rime_rate_qs_in = rime_qs;
        codi::RealReverse T_prime_in = T_prime;

        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;

        std::cout << "qc,Nc,qr,Nr,qs,Ns,dep_rate_snow,rime_rate_qc,"
                  << "rime_rate_nc,rime_rate_qr,rime_rate_nr,"
                  << "rime_rate_qs,T,dt,"
                  << "delta_qc,delta_Nc,delta_qi,delta_Ni,"
                  << "delta_qr,delta_Nr,delta_qs,delta_Ns,"
                  << "delta_qg,delta_Ng,"
                  << "delta_lat_cool_idx,delta_lat_heat_idx\n";

        for(uint32_t i=0; i<=n1; ++i)
        {
            if(qc_min != NOT_USED && qc_max != NOT_USED)
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
            else if(qs_min != NOT_USED && qs_max != NOT_USED)
            {
                qs_prime_in = i * (qs_max-qs_min) / n1 + qs_min;
                used_parameter = 1;
            }
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 2;
            }
            else if(dep_snow_min != NOT_USED && dep_snow_max != NOT_USED)
            {
                dep_rate_snow_in = i * (dep_snow_max-dep_snow_min) / n1 + dep_snow_min;
                used_parameter = 3;
            }
            else if(rime_qc_min != NOT_USED && rime_qc_max != NOT_USED)
            {
                rime_rate_qc_in = i * (rime_qc_max-rime_qc_min) / n1 + rime_qc_min;
                used_parameter = 4;
            }
            else if(temp_min != NOT_USED && temp_max != NOT_USED)
            {
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
                used_parameter = 5;
            }
            else if(rime_qr_min != NOT_USED && rime_qr_max != NOT_USED)
            {
                rime_rate_qr_in = i * (rime_qr_max-rime_qr_min) / n1 + rime_qr_min;
                used_parameter = 6;
            }
            else if(rime_qs_min != NOT_USED && rime_qs_max != NOT_USED)
            {
                rime_rate_qs_in = i * (rime_qs_max-rime_qs_min) / n1 + rime_qs_min;
                used_parameter = 7;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qs_min != NOT_USED && qs_max != NOT_USED)
                {
                    qs_prime_in = j * (qs_max-qs_min) / n2 + qs_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && dep_snow_min != NOT_USED && dep_snow_max != NOT_USED)
                {
                    dep_rate_snow_in = j * (dep_snow_max-dep_snow_min) / n2 + dep_snow_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && rime_qc_min != NOT_USED && rime_qc_max != NOT_USED)
                {
                    rime_rate_qc_in = j * (rime_qc_max-rime_qc_min) / n2 + rime_qc_min;
                    used_parameter2 = 4;
                }
                else if(used_parameter < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                {
                    T_prime_in = j * (temp_max-temp_min) / n2 + temp_min;
                    used_parameter2 = 5;
                }
                else if(used_parameter < 6 && rime_qr_min != NOT_USED && rime_qr_max != NOT_USED)
                {
                    rime_rate_qr_in = j * (rime_qr_max-rime_qr_min) / n2 + rime_qr_min;
                    used_parameter2 = 6;
                }
                else if(used_parameter < 4 && rime_qs_min != NOT_USED && rime_qs_max != NOT_USED)
                {
                    rime_rate_qs_in = j * (rime_qs_max-rime_qs_min) / n2 + rime_qs_min;
                    used_parameter2 = 7;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && qr_min != NOT_USED && qr_max != NOT_USED)
                        qr_prime_in = k * (qr_max-qr_min) / n3 + qr_min;
                    else if(used_parameter2 < 3 && dep_snow_min != NOT_USED && dep_snow_max != NOT_USED)
                        dep_rate_snow_in = k * (dep_snow_max-dep_snow_min) / n3 + dep_snow_min;
                    else if(used_parameter2 < 4 && rime_qc_min != NOT_USED && rime_qc_max != NOT_USED)
                        rime_rate_qc_in = k * (rime_qc_max-rime_qc_min) / n3 + rime_qc_min;
                    if(used_parameter2 < 5 && temp_min != NOT_USED && temp_max != NOT_USED)
                        T_prime_in = k * (temp_max-temp_min) / n3 + temp_min;
                    else if(used_parameter2 < 6 && rime_qr_min != NOT_USED && rime_qr_max != NOT_USED)
                        rime_rate_qr_in = k * (rime_qr_max-rime_qr_min) / n3 + rime_qr_min;
                    else if(used_parameter2 < 7 && rime_qs_min != NOT_USED && rime_qs_max != NOT_USED)
                        rime_rate_qs_in = k * (rime_qs_max-rime_qs_min) / n3 + rime_qs_min;

                    for(auto& val: y)
                        val = 0;

                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse Nr = qr_prime_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                    codi::RealReverse Ns = qs_prime_in
                        / ( (cc.snow.max_x - cc.snow.min_x)/2 + cc.snow.min_x );

                    codi::RealReverse rime_rate_nc = rime_rate_qc_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse rime_rate_nr = rime_rate_qr_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );

                    std::cout << qc_prime_in.getValue() << ","
                              << Nc.getValue() << ","
                              << qr_prime_in.getValue() << ","
                              << Nr.getValue() << ","
                              << qs_prime_in.getValue() << ","
                              << Ns.getValue() << ","
                              << dep_rate_snow_in.getValue() << ","
                              << rime_rate_qc_in.getValue() << ","
                              << rime_rate_nc.getValue() << ","
                              << rime_rate_qr_in.getValue() << ","
                              << rime_rate_nr.getValue() << ","
                              << rime_rate_qs_in.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << cc.dt_prime << ",";

                    snow_riming(qc_prime_in, Nc, qr_prime_in, Nr,
                        qs_prime_in, Ns, dep_rate_snow_in, rime_rate_qc_in,
                        rime_rate_nc, rime_rate_qr_in, rime_rate_nr,
                        rime_rate_qs_in, T_prime_in, cc.dt_prime, y, cc);

                    std::cout << y[qc_idx].getValue() << ","
                              << y[Nc_idx].getValue() << ","
                              << y[qi_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[qr_idx].getValue() << ","
                              << y[Nr_idx].getValue() << ","
                              << y[qs_idx].getValue() << ","
                              << y[Ns_idx].getValue() << ","
                              << y[qg_idx].getValue() << ","
                              << y[Ng_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("particle_cloud_riming") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qc_prime_in = qc_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse q2_prime_in, q2_min, q2_max;
        codi::RealReverse avg_size;
        uint32_t q2_idx, N2_idx;
        char q2 = 'f';
        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;
        collection_model_constants_t *coeffs;
        particle_model_constants_t *pc;
        if(qg_prime != NOT_USED || qg_min != NOT_USED || qg_prime != NOT_USED)
        {
            q2 = 'g';
            q2_min = qg_min;
            q2_max = qg_max;
            q2_prime_in = qg_prime;
            avg_size = (cc.graupel.max_x - cc.graupel.min_x) / 2 + cc.graupel.min_x;
            q2_idx = qg_idx;
            N2_idx = Ng_idx;
            coeffs = &cc.coeffs_gcr;
            pc = &cc.graupel;
        }
        else if(qh_prime != NOT_USED || qh_min != NOT_USED || qh_prime != NOT_USED)
        {
            q2 = 'h';
            q2_min = qh_min;
            q2_max = qh_max;
            q2_prime_in = qh_prime;
            avg_size = (cc.hail.max_x - cc.hail.min_x) / 2 + cc.hail.min_x;
            q2_idx = qh_idx;
            N2_idx = Nh_idx;
            coeffs = &cc.coeffs_hcr;
            pc = &cc.hail;
        }

        if(q2 == 'f')
        {
            std::cout << "particle_cloud_riming needs either hail "
                     <<  "or graupel particles\nABORTING!\n";
            return 1;
        }
        std::cout << "qc,Nc,T,S,p,qv,q" << q2 << ",N" << q2
                  << ",delta_qc,delta_Nc,delta_q" << q2
                  << ",delta_N" << q2 << ","
                  << "delta_qi,delta_Ni,delta_qr,delta_Nr,delta_lat_cool,delta_lat_heat\n";
        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qc_min != NOT_USED && qc_max != NOT_USED)
            {
                qc_prime_in = i * (qc_max-qc_min) / n1 + qc_min;
                used_parameter = 1;
            }
            else if(q2_min != NOT_USED && q2_max != NOT_USED)
            {
                q2_prime_in = i * (q2_max-q2_min) / n1 + q2_min;
                used_parameter = 2;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qc_min != NOT_USED && qc_max != NOT_USED)
                {
                    qc_prime_in = j * (qc_max-qc_min) / n2 + qc_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && q2_min != NOT_USED && q2_max != NOT_USED)
                {
                    q2_prime_in = j * (q2_max-q2_min) / n2 + q2_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && q2_min != NOT_USED && q2_max != NOT_USED)
                        q2_prime_in = k * (q2_max-q2_min) / n3 + q2_min;
                    else if(used_parameter2 < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse Nc = qc_prime_in
                        / ( (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x );
                    codi::RealReverse N2 = q2_prime_in / avg_size;
                    // Update vertical velocity parameters
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
                    cc.hail.rho_v = exp(-rho_vel * rho_inter);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);
                    std::cout << qc_prime_in.getValue() << ","
                              << Nc.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << p_prime_in.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << q2_prime_in.getValue() << ","
                              << N2.getValue() << ",";

                    particle_cloud_riming(qc_prime_in, Nc, T_prime_in,
                        q2_prime_in, N2, y[q2_idx], y[N2_idx], *coeffs, *pc, y, cc);

                    std::cout << y[qc_idx].getValue() << ","
                              << y[Nc_idx].getValue() << ","
                              << y[q2_idx].getValue() << ","
                              << y[N2_idx].getValue() << ","
                              << y[qi_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[qr_idx].getValue() << ","
                              << y[Nr_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("particle_rain_riming") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        codi::RealReverse qv_prime_in = qv_prime;
        codi::RealReverse p_prime_in = p_prime;
        codi::RealReverse q2_prime_in, q2_min, q2_max;
        codi::RealReverse avg_size;
        uint32_t q2_idx, N2_idx;
        char q2 = 'f';
        uint32_t used_parameter = 0;
        uint32_t used_parameter2 = 0;
        collection_model_constants_t *coeffs;
        particle_model_constants_t *pc;
        if(qg_prime != NOT_USED || qg_min != NOT_USED || qg_prime != NOT_USED)
        {
            q2 = 'g';
            q2_min = qg_min;
            q2_max = qg_max;
            q2_prime_in = qg_prime;
            avg_size = (cc.graupel.max_x - cc.graupel.min_x) / 2 + cc.graupel.min_x;
            q2_idx = qg_idx;
            N2_idx = Ng_idx;
            coeffs = &cc.coeffs_grr;
            pc = &cc.graupel;
        }
        else if(qh_prime != NOT_USED || qh_min != NOT_USED || qh_prime != NOT_USED)
        {
            q2 = 'h';
            q2_min = qh_min;
            q2_max = qh_max;
            q2_prime_in = qh_prime;
            avg_size = (cc.hail.max_x - cc.hail.min_x) / 2 + cc.hail.min_x;
            q2_idx = qh_idx;
            N2_idx = Nh_idx;
            coeffs = &cc.coeffs_hrr;
            pc = &cc.hail;
        }

        if(q2 == 'f')
        {
            std::cout << "particle_rain_riming needs either hail "
                     <<  "or graupel particles\nABORTING!\n";
            return 1;
        }
        std::cout << "qr,Nr,T,S,qv,q" << q2 << ",N" << q2
                  << ",delta_qr,delta_Nr,delta_q" << q2
                  << ",delta_N" << q2 << ","
                  << "delta_qi,delta_Ni,delta_lat_cool,delta_lat_heat\n";
        for(uint32_t i=0; i<=n1; ++i)
        {

            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 1;
            }
            else if(q2_min != NOT_USED && q2_max != NOT_USED)
            {
                q2_prime_in = i * (q2_max-q2_min) / n1 + q2_min;
                used_parameter = 2;
            }
            else if(qv_min != NOT_USED && qv_max != NOT_USED)
            {
                qv_prime_in = i * (qv_max-qv_min) / n1 + qv_min;
                used_parameter = 3;
            }
            else if(p_min != NOT_USED && p_max != NOT_USED)
            {
                p_prime_in = i * (p_max-p_min) / n1 + p_min;
                used_parameter = 4;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qr_min != NOT_USED && qr_max != NOT_USED)
                {
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;
                    used_parameter2 = 1;
                }
                else if(used_parameter < 2 && q2_min != NOT_USED && q2_max != NOT_USED)
                {
                    q2_prime_in = j * (q2_max-q2_min) / n2 + q2_min;
                    used_parameter2 = 2;
                }
                else if(used_parameter < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                {
                    qv_prime_in = j * (qv_max-qv_min) / n2 + qv_min;
                    used_parameter2 = 3;
                }
                else if(used_parameter < 4 && p_min != NOT_USED && p_max != NOT_USED)
                {
                    p_prime_in = j * (p_max-p_min) / n2 + p_min;
                    used_parameter2 = 4;
                }

                for(uint32_t k=0; k<=n3; ++k)
                {
                    if(used_parameter2 < 2 && q2_min != NOT_USED && q2_max != NOT_USED)
                        q2_prime_in = k * (q2_max-q2_min) / n3 + q2_min;
                    else if(used_parameter2 < 3 && qv_min != NOT_USED && qv_max != NOT_USED)
                        qv_prime_in = k * (qv_max-qv_min) / n3 + qv_min;
                    else if(used_parameter2 < 4 && p_min != NOT_USED && p_max != NOT_USED)
                        p_prime_in = k * (p_max-p_min) / n3 + p_min;

                    for(auto& val: y)
                        val = 0;
                    codi::RealReverse Nr = qr_prime_in
                        / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                    codi::RealReverse N2 = q2_prime_in / avg_size;
                    // Update vertical velocity parameters
                    codi::RealReverse S = qv_prime_in * Rv * T_prime_in
                        / saturation_pressure_water_icon(T_prime_in);
                    codi::RealReverse rho_inter = log(compute_rhoh(p_prime_in, T_prime_in, S)/rho_0);
                    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
                    cc.rain.rho_v = exp(-rho_vel * rho_inter);
                    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
                    cc.hail.rho_v = exp(-rho_vel * rho_inter);
                    cc.ice.rho_v = exp(-rho_vel * rho_inter);
                    cc.snow.rho_v = exp(-rho_vel * rho_inter);
                    std::cout << qr_prime_in.getValue() << ","
                              << Nr.getValue() << ","
                              << T_prime_in.getValue() << ","
                              << S.getValue() << ","
                              << qv_prime_in.getValue() << ","
                              << q2_prime_in.getValue() << ","
                              << N2.getValue() << ",";

                    particle_rain_riming(qr_prime_in, Nr, T_prime_in,
                        q2_prime_in, N2, y[q2_idx], y[N2_idx], *coeffs, *pc, y, cc);

                    std::cout << y[qr_idx].getValue() << ","
                              << y[Nr_idx].getValue() << ","
                              << y[q2_idx].getValue() << ","
                              << y[N2_idx].getValue() << ","
                              << y[qi_idx].getValue() << ","
                              << y[Ni_idx].getValue() << ","
                              << y[lat_cool_idx].getValue() << ","
                              << y[lat_heat_idx].getValue() << "\n";
                }
            }
        }
    } else if(func_name.compare("rain_freeze") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qr_prime_in = qr_prime;
        std::cout << "qr,Nr,T,dt,delta_qr,delta_Nr,delta_qs,delta_Ns,"
                  << "delta_qg,delta_Ng,delta_qh,delta_Nh,delta_lat_cool,delta_lat_heat\n";
        uint32_t used_parameter = 0;
        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qr_min != NOT_USED && qr_max != NOT_USED)
            {
                qr_prime_in = i * (qr_max-qr_min) / n1 + qr_min;
                used_parameter = 1;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qr_min != NOT_USED && qr_max != NOT_USED)
                    qr_prime_in = j * (qr_max-qr_min) / n2 + qr_min;

                for(auto& val: y)
                    val = 0;

                codi::RealReverse Nr = qr_prime_in
                    / ( (cc.rain.max_x - cc.rain.min_x)/2 + cc.rain.min_x );
                std::cout << qr_prime_in.getValue() << "," << Nr.getValue() << ","
                          << T_prime_in.getValue() << "," << cc.dt_prime << ",";
                rain_freeze(qr_prime_in, Nr, T_prime_in, cc.dt_prime, y, cc);
                std::cout << y[qr_idx].getValue() << "," << y[Nr_idx].getValue() << ","
                          << y[qs_idx].getValue() << "," << y[Ns_idx].getValue() << ","
                          << y[qg_idx].getValue() << "," << y[Ng_idx].getValue() << ","
                          << y[qh_idx].getValue() << "," << y[Nh_idx].getValue() << ","
                          << y[lat_cool_idx].getValue() << ","
                          << y[lat_heat_idx].getValue() << "\n";
            }
        }
    } else if(func_name.compare("ice_melting") == 0)
    {
        codi::RealReverse T_prime_in = T_prime;
        codi::RealReverse qi_prime_in = qi_prime;
        std::cout << "qi,Ni,T,dt,delta_qi,delta_Ni,delta_qr,delta_Nr,"
                  << "delta_qc,delta_Nc,delta_lat_cool,delta_lat_heat\n";
        uint32_t used_parameter = 0;
        for(uint32_t i=0; i<=n1; ++i)
        {
            if(temp_min != NOT_USED && temp_max != NOT_USED)
                T_prime_in = i * (temp_max-temp_min) / n1 + temp_min;
            else if(qi_min != NOT_USED && qi_max != NOT_USED)
            {
                qi_prime_in = i * (qi_max-qi_min) / n1 + qi_min;
                used_parameter = 1;
            }

            for(uint32_t j=0; j<=n2; ++j)
            {
                if(used_parameter < 1 && qi_min != NOT_USED && qi_max != NOT_USED)
                    qi_prime_in = j * (qi_max-qi_min) / n2 + qi_min;

                codi::RealReverse qi_in = qi_prime_in/ref_quant.qref;
                codi::RealReverse Ni = qi_prime_in
                    / ( (cc.ice.max_x - cc.ice.min_x_melt)/2 + cc.ice.min_x_melt );
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
}