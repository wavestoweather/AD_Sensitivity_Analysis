#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <codi.hpp>

#include "include/microphysics/constants.h"
#include "include/misc/error.h"


/**
 * Struct to hold all model constants regarding particles that we want to
 * know the gradients using AD of.
 * Has getter and register functions for codi::RealReverseIndex for its
 * members.
 */
template<class float_t>
struct particle_model_constants_t{
    std::vector<float_t> constants;
    std::vector<uint32_t> perturbed_idx;

    /*
    * Uncertainty for every parameter related to cloud water vapor.
    * Currently sets everything to 10% of every parameter.
    */
    std::array<double, static_cast<uint32_t>(Particle_cons_idx::n_items)> uncertainty;

    particle_model_constants_t();

    /**
     * Register the model parameters on the tape for codi::RealReverseIndex.
     *
     * @param tape Tape where the parameters should be registered on.
     */
    void register_input(
        codi::RealReverseIndex::Tape &tape,
        uint32_t &id);

    /**
     * Get the gradients of all its members. You need to register them on a
     * type before to get meaningful values.
     *
     * @param out_vec On out: Stores all gradients.
     * @param idx Start index of out_vec where the gradients should be stored.
     * @param info Can be used for debugging.
     */
    void get_gradient(
        std::array<double, num_par> &out_vec,
        uint32_t &idx,
        const double &ref_value) const;

    // void put(pt::ptree &ptree, const std::string &type_name) const;

    /**
     * Set any perturbed parameter from the property tree.
     *
     * @params ptree Property tree with key = idx of constants, values =
     *         the perturbed values and one list 'perturbed' of indices.
     *
     * @returns Errorcode
     */
    // int from_pt(pt::ptree &ptree);

    /**
     * Print the constants used for this particle.
     *
     * @param title Name of particle
     */
    void print(const std::string &title, std::ostream &os = std::cout);

    int from_json(const nlohmann::json& j);
};

template<class float_t>
void to_json(nlohmann::json& j, const particle_model_constants_t<float_t>& p);
