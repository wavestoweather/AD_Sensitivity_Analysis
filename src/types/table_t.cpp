#include "include/types/table_t.h"

table_t::table_t()
{

}

codi::RealReverse table_t::get(
    uint64_t i,
    uint64_t j,
    uint64_t k,
    uint64_t l) const
{
    return table[i*n3*n2*n1 + j*n2*n1 + k*n1 + l];
}