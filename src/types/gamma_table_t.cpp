#include "include/types/gamma_table_t.h"


gamma_table_t::gamma_table_t() {
}


template <class A>
A gamma_table_t::look_lo(A x) const {
    A xt = max(min(x, this->x[this->n_bins-1]), 0.0);
    double tmp = xt.getValue()*this->odx;
    tmp = floor(tmp);
    uint64_t iu = std::min((uint64_t) tmp, this->n_bins-2);
    uint64_t io = iu + 1;
    return this->igf[iu] + (this->igf[io] - this->igf[iu]) * this->odx*(xt-this->x[iu]);
}


template <class A>
A gamma_table_t::look_up(
    A x) const {
    A xt = max(min(x, this->x[this->n_bins-1]), 0.0);
    double tmp = xt.getValue()*this->odx;
    tmp = floor(tmp);
    uint64_t iu = std::min((uint64_t) tmp, this->n_bins-2);
    uint64_t io = iu + 1;
    A lookup = this->igf[this->n_bins-1] - this->igf[iu]
        - (this->igf[io] - this->igf[iu]) * this->odx*(xt-this->x[iu]);
    return max(lookup, codi::RealReverseIndex(0.0));
}


void gamma_table_t::init_gamma_table(
    const uint64_t &nl,
    const uint64_t &nl_highres,
    const double &a) {
    const double c1 = 36.629433904824623;
    const double c2 = -0.119475603955226;
    const double c3 = 0.339332937820052;
    const double c4 = 1.156369000458310;

    n_bins = nl;
    n_bins_highres = nl_highres;

    x.resize(nl);
    x_highres.resize(nl_highres);
    igf.resize(nl);
    igf_highres.resize(nl_highres);

    // Low resolution
    // maximum x-value (99.5%)
    x[n_bins-2] = c1 * (1.0-exp(c2*pow(a, c3))) + c4*a;
    dx = x[n_bins-2] / (n_bins-2.0);
    odx = 1.0/dx;
    for (uint64_t i=0; i < n_bins-2; ++i) {
        x[i] = (i-1) * dx;
        igf[i] = boost::math::tgamma_lower(a, x[i]);
    }
    // for x -> infinity:
    x[n_bins-1] = (n_bins-1)*dx;
    igf[n_bins-1] = std::tgamma(a);

    // High resolution (lowest 2% of the x-values)
    dx_highres = x[std::round(0.01*(n_bins-1))] / (n_bins_highres - 1.0);
    odx_highres = 1.0/dx_highres;
    for (uint64_t i=0; i < n_bins_highres; ++i) {
        x_highres[i] = (i-1) * dx_highres;
        igf_highres[i] = boost::math::tgamma_lower(a, x_highres[i]);
    }
}

template codi::RealReverseIndex gamma_table_t::look_lo<codi::RealReverseIndex>(codi::RealReverseIndex) const;
template codi::RealReverseIndex gamma_table_t::look_up<codi::RealReverseIndex>(codi::RealReverseIndex) const;

template codi::RealForwardVec<num_par_init> gamma_table_t::look_lo<codi::RealForwardVec<num_par_init> >(
    codi::RealForwardVec<num_par_init>) const;
template codi::RealForwardVec<num_par_init> gamma_table_t::look_up<codi::RealForwardVec<num_par_init> >(
    codi::RealForwardVec<num_par_init>) const;

