#include "include/types/table_t.h"


template<class float_t>
table_t<float_t>::table_t() {
}


template<class float_t>
float_t table_t<float_t>::get(
    const uint64_t i,
    const uint64_t j,
    const uint64_t k,
    const uint64_t l) const {
    return table[i*n3*n2*n1 + j*n2*n1 + k*n1 + l];
}


template class table_t<codi::RealReverse>;
template class table_t<codi::RealForwardVec<num_par_init> >;
