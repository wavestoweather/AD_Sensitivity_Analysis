#pragma once

#include <stdlib.h>
#include <netcdf.h>

#include <cmath>
#include <fstream>
#include <iterator>
#include <iostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>
#include "codi.hpp"

#include "include/misc/error.h"
#include "include/microphysics/constants.h"
#include "include/types/global_args_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/output_handle_t.h"
#include "include/types/model_constants_t.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/segment_t.h"
#include "include/types/table_t.h"

#include "include/misc/general.h"

/** @defgroup io Input/Output Functions
 * Functions to load data and setup and initialize structs.
 * @{
 */


/**
 * Parse an ensemble configuration as an XML-file.
 *
 */
template<class float_t>
void load_ens_config(
    std::string filename,
    model_constants_t<float_t> &cc,
    std::vector<segment_t> &segments,
    input_parameters_t &input,
    const reference_quantities_t &ref_quant) {

    std::ifstream i(filename);
    nlohmann::json config;
    i >> config;
    std::string ens_desc;
    cc.max_n_trajs = 0;
    if (config.find("segments") != config.end()) {
        for (const auto &s_config : config["segments"]) {
            segment_t s;
            s.from_json(s_config, cc);
            segments.push_back(s);
            if (s.n_members > cc.max_n_trajs)
                cc.max_n_trajs = s.n_members;
            if (s.activated)
                s.perturb(cc, ref_quant, input, ens_desc);
        }
    }

    cc.ens_desc += ens_desc;
    if (input.n_ensembles != 0) {
        cc.n_ensembles = input.n_ensembles;
    } else {
        cc.n_ensembles = segments.size() + 1;
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
template<class float_t>
bool load_lookup_table(
    table_t<float_t> &table,
    std::string filename = "dmin_wetgrowth_lookup.dat",
    uint32_t ndT = 61) {
    // T may be nonequidistant
    // p, qw and qi must be equidistant
    std::ifstream data;
    data.open(filename);
    // Read file
    if (data.is_open()) {
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
        for (auto &val : line_vec) {
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
        for (auto &val : line_vec) {
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
        for (auto &val : line_vec) {
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
        for (auto &val : line_vec) {
            table.x4[counter] = std::stod(val);
            counter++;
        }

        // Read the values, one value per line
        for (uint64_t n4=0; n4 < table.n4; ++n4)
            for (uint64_t n3=0; n3 < table.n3; ++n3)
                for (int i=0; i < anzT_wg_loc; ++i)
                    for (uint64_t n1=0; n1 < table.n1; ++n1) {
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

        for (uint64_t i=0; i < table.n2; ++i)
            table.x2[i] = minT + i*table.dx2;

        // Linear interpolation w.r.t. T
        for (uint64_t i=0; i < table.n2; ++i) {
            uint64_t iu = 0;
            for (int j=0; j < anzT_wg_loc-1; ++j) {
                if (table.x2[i] >= Tvec_wg_g_loc[j]
                    && table.x2[i] <= Tvec_wg_g_loc[j+1]) {
                    iu = j;
                    break;
                }
            }
            uint64_t io = iu+1;

            // actual linear interpolation
            for (uint64_t n4=0; n4 < table.n4; ++n4)
                for (uint64_t n3=0; n3 < table.n3; ++n3)
                    for (uint64_t n1=0; n1 < table.n1; ++n1)
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
    } else {
        std::cerr << "Error loading " << filename << ". Does the file exist?\n";
        return false;
    }
}

/** @} */  // end of group io
