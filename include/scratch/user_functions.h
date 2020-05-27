#ifndef SCRATCH_USER_FUNCTIONS_H
#define SCRATCH_USER_FUNCTIONS_H

#include "codi.hpp"
#include <math.h>
#include <cmath>
#include <tgmath.h>
#include <vector>
#include "include/microphysics/types.h"
#include "include/microphysics/constants.h"
#include "include/microphysics/physical_parameterizations.h"


codi::RealReverse trunc(codi::RealReverse x)
{
    int left = x.getValue();
    codi::RealReverse r = left;
    r.setGradient(x.getGradient());
    return r;
}

template<class float_t>
void ice_activation(
    float_t qc,
    float_t qv,
    float_t T,
    std::vector<float_t> &y,
    const reference_quantities_t &ref,
    model_constants_t &cc)
{
    float_t T_prime  = ref.Tref * T;
    float_t qc_prime = ref.qref * qc;
    float_t qv_prime = ref.qref * qv;
    float_t p_sat_ice = saturation_pressure_ice(T_prime);
    bool use_prog_in = false;

    float_t S_si = qv_prime * Rv * T_prime / p_sat_ice;

    double n_inact = 0;

    double na_dust;
    double na_soot;
    double na_orga;
    if(nuc_type == 6)
    {
        na_dust = 160e4;
        na_soot = 30e6;
        na_orga = 0.0;
    } else if(nuc_type == 7)
    {
        // Standard values
        na_dust = 160e4;
        na_soot = 25e6;
        na_orga = 30e6;
    } else if(nuc_type == 8)
    {
        na_dust = 70e4;
        na_soot = 0.0;
        na_orga = 0.0;
    }

    const double EPSILON = 1.0e-20;

    if(T_prime < T_nuc && T_prime > 180.0 && S_si > 1.0
        && n_inact < ni_het_max)
    {
        float_t x_t = (274.0 - T_prime) / t_tstep;
        x_t = min(x_t, t_tmax-1);
        int tt = (int) x_t.getValue() - 1;

        std::vector<float_t> infrac(3);
        if(qc_prime > EPSILON)
        {
            // Immersion freezing at water saturation
            infrac[0] = ( trunc(x_t)+1.0-x_t ) * afrac_dust[98][tt]
                + ( x_t-trunc(x_t) ) * afrac_dust[98][tt+1];
            infrac[1] = ( trunc(x_t)+1.0-x_t ) * afrac_soot[98][tt]
                + ( x_t-trunc(x_t) ) * afrac_soot[98][tt+1];
            infrac[2] = ( trunc(x_t)+1.0-x_t ) * afrac_orga[98][tt]
                + ( x_t-trunc(x_t) ) * afrac_orga[98][tt+1];
        } else
        {
            // deposition nucleation below water saturation
            // Indices for 2D look-up tables
            float_t x_s = 100.0*(S_si-1.0) / s_sstep;
            x_s = min(x_s, s_smax-1);
            int ss = std::max(0, (int) (x_s.getValue()-1) );
            float_t S_sr = max(1.0, trunc(x_s));
            infrac[0] = ( trunc(x_t)+1.0-x_t ) * ( S_sr+1.0-x_s )
                * afrac_dust[ss][tt]
                + ( x_t-trunc(x_t) ) * ( S_sr+1.0-x_s )
                * afrac_dust[ss][tt+1]
                + ( x_t-trunc(x_t)+1.0-x_t ) * ( x_s-S_sr )
                * afrac_dust[ss+1][tt]
                + ( x_t-trunc(x_t) ) * ( x_s-S_sr )
                * afrac_dust[ss+1][tt+1];
            infrac[1] = ( trunc(x_t)+1.0-x_t ) * ( S_sr+1.0-x_s )
                * afrac_soot[ss][tt]
                + ( x_t-trunc(x_t) ) * ( S_sr+1.0-x_s )
                * afrac_dust[ss][tt+1]
                + ( x_t-trunc(x_t)+1.0-x_t ) * ( x_s-S_sr )
                * afrac_soot[ss+1][tt]
                + ( x_t-trunc(x_t) ) * ( x_s-S_sr )
                * afrac_soot[ss+1][tt+1];
            infrac[2] = ( trunc(x_t)+1.0-x_t ) * ( S_sr+1.0-x_s )
                * afrac_orga[ss][tt]
                + ( x_t-trunc(x_t) ) * ( S_sr+1.0-x_s )
                * afrac_orga[ss][tt+1]
                + ( x_t-trunc(x_t)+1.0-x_t ) * ( x_s-S_sr )
                * afrac_orga[ss+1][tt]
                + ( x_t-trunc(x_t) ) * ( x_s-S_sr )
                * afrac_dust[ss+1][tt+1];
        }
        float_t ndiag = infrac[0] + na_soot*infrac[1] + na_orga*infrac[2];
        float_t ndiag_dust;
        float_t ndiag_all;
        if(use_prog_in)
        {
            // Not implemented
            // n_inpot replaces n_dust
            // ndiag = n_inpot * ndiag;
            // ndiag_dust = n_inpot*infrac[0];
            // ndiag_all = ndiag;
        } else
        {
            ndiag = na_dust * ndiag;
        }
        std::cout << ndiag << "," << ni_het_max << "," << n_inact << "," 
                  << S_si << "," << T_prime << ",";
        ndiag = min(ndiag, ni_het_max);
        float_t delta_n = max(ndiag-n_inact, 0.0);
        // TODO: Check if min_x or min_x_nuc
        float_t delta_q = min(delta_n*cc.ice.min_x, qv_prime);
        delta_n = delta_q/cc.ice.min_x;
        y[0] = delta_n;
        y[1] = delta_q;
    } else 
    {
        std::cout << "0," << ni_het_max << "," << n_inact << "," 
                  << S_si << "," << T_prime << ",";
    }
}


void ice_table_scan(
    uint32_t x,
    uint32_t y,
    std::vector<double> &out)
{
    out[0] = afrac_dust[x][y];
    out[1] = afrac_soot[x][y];
    out[2] = afrac_orga[x][y];
}

template<class float_t>
void infrac_table_scan(
    float_t T_prime,
    std::vector<float_t> &infrac)
{
    float_t x_t = (274.0 - T_prime) / t_tstep;
    x_t = min(x_t, t_tmax-1);
    int tt = (int) x_t.getValue() - 1;
    std::cout << "98" << "," << x_t << ",";

    // Immersion freezing at water saturation
    infrac[0] = ( trunc(x_t)+1.0-x_t ) * afrac_dust[98][tt]
        + ( x_t-trunc(x_t) ) * afrac_dust[98][tt+1];
    infrac[1] = ( trunc(x_t)+1.0-x_t ) * afrac_soot[98][tt]
        + ( x_t-trunc(x_t) ) * afrac_soot[98][tt+1];
    infrac[2] = ( trunc(x_t)+1.0-x_t ) * afrac_orga[98][tt]
        + ( x_t-trunc(x_t) ) * afrac_orga[98][tt+1];
}
#endif