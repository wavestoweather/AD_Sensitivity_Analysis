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
    uint32_t n1 = 50;
    uint32_t n2 = 50;
    uint32_t n3 = 0;
    // if(argc > 1)
    //     n1 = (uint32_t) argv[1];
    // if(argc > 2)
    //     n2 = (uint32_t) argv[2];
    // if(argc > 3)
    //     n3 = (uint32_t) argv[3];

    reference_quantities_t ref_quant;

    ref_quant.Tref = 1.0;
    ref_quant.qref = 1.0;

    model_constants_t cc;
    cc.ice.min_x = 1.0e-12;
    codi::RealReverse qc = n3;
    std::string out_filename = "scan_qv_T_ice_activation.csv";
    std::ofstream outfile;
    outfile.open(out_filename);

    if( !outfile.is_open() )
    {
        std::cout << "ERROR while opening the outputfile. Aborting." << std::endl;
        return 1;
    }
    outfile << "qv,T,delta_n,delta_q\n";

    for(uint32_t i=0; i<n1; i++)
    {
        codi::RealReverse qv = 0.0 + 0.1/n1 * i;
        for(uint32_t j=0; j<n2; j++)
        {
            codi::RealReverse T = 180.0 + (T_nuc-180.0)/n2 * j;
            y[0] = y[1] = 0.0;
            ice_activation(qc, qv, T, y, ref_quant, cc);
            outfile << qv << "," << T << "," << y[0].getValue() << "," << y[1].getValue() << "\n";
        }
    }
    outfile.close();
}