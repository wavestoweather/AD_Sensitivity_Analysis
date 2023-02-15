#pragma once

#include <codi.hpp>

/**
 * Struct to hold model constants for particle collection that we *usually*
 * are not interested in regarding their influence on the overall model.
 * It is stored within a struct model_constants_t for different collisions
 * such as hail and ice collection or ice and rain riming.
 */
template<class float_t>
struct collection_model_constants_t{
    float_t delta_n_aa;
    float_t delta_n_ab;
    float_t delta_n_bb;
    float_t delta_q_aa;
    float_t delta_q_ab;
    float_t delta_q_bb;
    float_t delta_q_ba;

    float_t theta_n_aa;
    float_t theta_n_ab;
    float_t theta_n_bb;
    float_t theta_q_aa;
    float_t theta_q_ab;
    float_t theta_q_bb;
    float_t theta_q_ba;
};
