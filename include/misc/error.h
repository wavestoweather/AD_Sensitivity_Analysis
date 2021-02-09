#ifndef ERROR_H
#define ERROR_H

#include <stdio.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////
// Error codes
////////////////////////////////////////////////////////////////////////////////
enum
{
    SUCCESS = 0,
    MISSING_PARAM_CONFIG_ERR = 1,
    MISSING_VARIANCE_CONFIG_ERR = 2,
    PARAM_CONFIG_ERR = 3,
    OUTPARAM_CONFIG_ERR = 4,
    N_MEMBERS_CONFIG_ERR = 5,
    VALUE_NAME_CONFIG_ERR = 6,
    METHOD_CONFIG_ERR = 7,
    N_SEGMENTS_CONFIG_ERR = 8,
    INPUT_NAME_CHECKPOINT_ERR = 9,
    NC_TRAJ_IDX_ERR = 10,
    MODEL_CONS_CHECKPOINT_ERR = 11,
    SEGMENTS_CHECKPOINT_ERR = 12,
    DISTRIBUTION_CONFIG_ERR = 13,
    ARGUMENT_ERR = 14,
    NC_ERR = 15,
    CHECKPOINT_LOAD_ERR = 16
};

// Credits to http://www.gpi-site.com/ for this function and the name
#define SUCCESS_OR_DIE(f...)                                                \
    do                                                                      \
    {                                                                       \
        const int err = f;                                                  \
        if(err != SUCCESS)                                                  \
        {                                                                   \
            std::cout << "Error: " << #f << " [" << __FILE__                \
                << ":" << __LINE__ << "]: " << err << "\n" << std::flush;   \
            exit(EXIT_FAILURE);                                             \
        }                                                                   \
    } while(0)

#endif