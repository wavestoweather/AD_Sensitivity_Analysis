#include "codi.hpp"
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
#include "include/scratch/user_functions.h"


int main(int argc, char** argv)
{
    std::vector<codi::RealReverse> y(2);
    uint32_t n1 = 500;
    uint32_t n2 = 500;
    codi::RealReverse min_1 = 0.0;
    codi::RealReverse min_2 = 180.0;
    codi::RealReverse max_1 = 0.1;
    codi::RealReverse max_2 = T_nuc;
    codi::RealReverse qc = 0.0;

    if(argc > 1)
        n1 = std::stoi(argv[1]);
    if(argc > 2)
        n2 = std::stoi(argv[2]);
    if(argc > 3)
        min_2 = std::stod(argv[3]);
    if(argc > 4)
        max_2 = std::stod(argv[4]);
    if(argc > 5)
        min_1 = std::stod(argv[5]);
    if(argc > 6)
        max_1 = std::stod(argv[6]);
    if(argc > 7)
        qc = std::stod(argv[7]);


    reference_quantities_t ref_quant;

    ref_quant.Tref = 1.0;
    ref_quant.qref = 1.0;

    model_constants_t cc;
    cc.ice.min_x = 1.0e-12;

    std::cout << "qv,T,delta_n,delta_q\n";

    for(uint32_t i=0; i<n1; i++)
    {
        codi::RealReverse qv = min_1 + max_1/n1 * i;
        for(uint32_t j=0; j<n2; j++)
        {
            codi::RealReverse T = min_2 + (max_2-min_2)/n2 * j;
            y[0] = y[1] = 0.0;
            ice_activation(qc, qv, T, y, ref_quant, cc);
            std::cout << qv << "," << T << "," << y[0].getValue() << "," << y[1].getValue() << "\n";
        }
    }
}