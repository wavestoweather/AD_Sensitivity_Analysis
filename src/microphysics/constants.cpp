#include "include/microphysics/constants.h"
/**
 * Yes, I know these aren't constants and more global variables.
 * Sorry.
 */

std::random_device rand_device{};
std::mt19937 rand_generator{rand_device()};

#if defined(TRACE_TIME)
// Relative to ascent time
double trace_time = 0;
bool trace = false;
#else
bool trace = true;
#endif

double sediment_q = 0;
double sediment_n = 0;
double sediment_q_total = 0;
double sediment_n_total = 0;

uint32_t auto_type = 3;
