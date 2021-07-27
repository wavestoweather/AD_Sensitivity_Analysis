#pragma once

#include <string>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include "codi.hpp"

#include "include/microphysics/constants.h"

namespace pt = boost::property_tree;

/**
 * Struct to hold all model constants regarding particles that we want to
 * know the gradients using AD of.
 * Has getter and register functions for codi::RealReverse for its
 * members.
 */
struct particle_model_constants_t{
    std::vector<codi::RealReverse> constants;
    std::vector<uint32_t> perturbed_idx;

    particle_model_constants_t();

    /**
     * Register the model parameters on the tape for codi::RealReverse.
     *
     * @param tape Tape where the parameters should be registered on.
     */
    void register_input(
        codi::RealReverse::TapeType &tape,
        uint32_t &id);

    /**
     * Get the gradients of all its members. You need to register them on a
     * type before to get meaningful values.
     *
     * @param out_vec On out: Stores all gradients.
     * @param idx Start index of out_vec where the gradients should be stored.
     */
    void get_gradient(
        std::array<double, num_par> &out_vec,
        uint32_t &idx) const;

    void put(pt::ptree &ptree, const std::string &type_name) const;

    /**
     * Set any perturbed parameter from the property tree.
     *
     * @params ptree Property tree with key = idx of constants, values =
     *         the perturbed values and one list 'perturbed' of indices.
     *
     * @returns Errorcode
     */
    int from_pt(pt::ptree &ptree);

    /**
     * Print the constants used for this particle.
     *
     * @param title Name of particle
     */
    void print(const std::string &title);
};
