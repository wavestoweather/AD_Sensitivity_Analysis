#pragma once

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////
// Error codes
////////////////////////////////////////////////////////////////////////////////
enum {
    SUCCESS = 0,
    // The first roughly hundert codes are reserved for MPI
    MISSING_PARAM_CONFIG_ERR = 101,
    MISSING_VARIANCE_CONFIG_ERR = 102,
    PARAM_CONFIG_ERR = 103,
    OUTPARAM_CONFIG_ERR = 104,
    N_MEMBERS_CONFIG_ERR = 105,
    VALUE_NAME_CONFIG_ERR = 106,
    METHOD_CONFIG_ERR = 107,
    N_SEGMENTS_CONFIG_ERR = 108,
    INPUT_NAME_CHECKPOINT_ERR = 109,
    NC_TRAJ_IDX_ERR = 110,
    MODEL_CONS_CHECKPOINT_ERR = 111,
    SEGMENTS_CHECKPOINT_ERR = 112,
    DISTRIBUTION_CONFIG_ERR = 113,
    ARGUMENT_ERR = 114,
    NC_ERR = 115,
    CHECKPOINT_LOAD_ERR = 116,
    PARAM_ADD_NAME_ERR = 117,
    INPUT_NAN_ERR = 118
};

// Credits to http://www.gpi-site.com/ for this function and the name
#define SUCCESS_OR_DIE(f...)                                                \
    do {                                                                    \
        const int err = f;                                                  \
        if (err != SUCCESS) {                                               \
            std::cerr << "Error: " << #f << " [" << __FILE__                \
                << ":" << __LINE__ << "]: " << err << "\n" << std::flush;   \
            MPI_Abort(MPI_COMM_WORLD, err);                                 \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while (0)
