#pragma once

#include "codi.hpp"

/**
 * Struct to hold model constants for particle collection that we *usually*
 * are not interested in regarding their influence on the overall model.
 * It is stored within a struct model_constants_t for different collisions
 * such as hail and ice collection or ice and rain riming.
 */
struct collection_model_constants_t{
    codi::RealReverse delta_n_aa;
    codi::RealReverse delta_n_ab;
    codi::RealReverse delta_n_bb;
    codi::RealReverse delta_q_aa;
    codi::RealReverse delta_q_ab;
    codi::RealReverse delta_q_bb;
    codi::RealReverse delta_q_ba;

    codi::RealReverse theta_n_aa;
    codi::RealReverse theta_n_ab;
    codi::RealReverse theta_n_bb;
    codi::RealReverse theta_q_aa;
    codi::RealReverse theta_q_ab;
    codi::RealReverse theta_q_bb;
    codi::RealReverse theta_q_ba;
};
