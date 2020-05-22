#include "codi.hpp"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tgmath.h>
#include <vector>

#include "include/microphysics/types.h"
#include "include/microphysics/constants.h"
#include "include/microphysics/general.h"
#include "include/microphysics/user_functions.h"
#include "include/microphysics/physical_parameterizations.h"



int main(int argc, char** argv)
{

    uint32_t n1 = 500;
    uint32_t n2 = 500;
    codi::RealReverse min_1 = 0.0;
    codi::RealReverse min_2 = 180.0;
    codi::RealReverse max_1 = 0.1;
    codi::RealReverse max_2 = T_nuc;
    codi::RealReverse qc = 0.0;

    if(argc > 2)
        n1 = std::stoi(argv[2]);
    if(argc > 3)
        n2 = std::stoi(argv[3]);
    if(argc > 4)
        min_1 = std::stod(argv[4]);
    if(argc > 5)
        max_1 = std::stod(argv[5]);
    if(argc > 6)
        min_2 = std::stod(argv[6]);
    if(argc > 7)
        max_2 = std::stod(argv[7]);
    if(argc > 8)
        qc = std::stod(argv[8]);


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

    model_constants_t cc;
    setup_model_constants(input, cc, ref_quant);

    std::vector<double> y(num_comp);

    if(std::strcmp("ice_table", argv[1]) == 0)
    {
        std::cout << "x,y,dust,soot,orga\n";
        for(uint32_t i=0; i<n1; i++)
        {
            for(uint32_t j=0; j<n2; j++)
            {
                ice_table_scan(i, j, y);
                std::cout << i << "," << j << ","
                          << y[0] << "," << y[1] << "," << y[2] << "\n";
            }
        }
    }
    else if(std::strcmp("ice_act_phillips", argv[1]) == 0)
    {

        bool use_prog_in = false;
        std::cout << "inact,Sat_ice,T,qv,delta_ni,delta_qi\n";
        for(uint32_t i=0; i<n1; i++)
        {
            codi::RealReverse qv = min_1 + (max_1-min_1)/n1 * i;
            for(uint32_t j=0; j<n2; j++)
            {
                codi::RealReverse T = min_2 + (max_2-min_2)/n2 * j;
                codi::RealReverse p_sat_ice = saturation_pressure_ice(T);
                codi::RealReverse ssi = qv * Rv * T / p_sat_ice;
                codi::RealReverse n_inact = 0;
                y[qi_idx] = y[Ni_idx] = 0.0;
                ice_activation_phillips(qc, qv, T, p_sat_ice, ssi,
                    n_inact, use_prog_in, y, cc);
                std::cout << n_inact << ","
                          << ssi << ","
                          << T << ","
                          << qv << ","
                          << y[qi_idx].getValue() << ","
                          << y[Ni_idx].getValue() << "\n";
            }
        }
    }
    else if(std::strcmp("infrac_scan", argv[1]) == 0)
    {
        std::cout << "x,y,T,dust,soot,orga\n";
        for(uint32_t i=0; i<n1; i++)
        {
            codi::RealReverse T = min_1 + (max_1-min_1)/n1 * i;
            infrac_table_scan(T, y);
            std::cout << T.getValue() << "," << y[0].getValue()
                      << "," << y[1].getValue() << "," << y[2].getValue() << "\n";
        }
    }
    else
    {
        std::cout << "No such method: " << argv[1] << "\n";
        return 1;
    }
}