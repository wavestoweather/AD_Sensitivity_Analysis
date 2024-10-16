#pragma once

#include <include/microphysics/constants.h>
#include <include/microphysics/user_functions.h>
#include <include/microphysics/program_io.h>

#include <iostream>
#include <string>
#include <vector>

#include <fstream>


struct physics_t {
 private:
    model_constants_t<codi::RealReverseIndex> cc;
    input_parameters_t input;
    reference_quantities_t ref_quant;
    std::vector<codi::RealReverseIndex> yold;
    std::vector<codi::RealReverseIndex> ynew;
    std::vector<codi::RealReverseIndex> ytmp;
    std::vector<codi::RealReverseIndex> k;
    std::vector< std::array<double, num_par > > gradients_;
    double uncertainty_scale;

 public:
    explicit physics_t(std::string table_path = "dmin_wetgrowth_lookup.dat");

    void set_limits(
            std::vector<codi::RealReverseIndex> &y,
            model_constants_t<codi::RealReverseIndex> &cc);

    int get_num_comp() {return num_comp;}
    int get_num_par() {return num_par;}

    void compute_nondimensional_effects(
            std::vector<codi::RealReverseIndex> &res,
            const codi::RealReverseIndex &p,
            const codi::RealReverseIndex &w,
            const codi::RealReverseIndex &w_prime,
            const codi::RealReverseIndex &T,
            const codi::RealReverseIndex &qv,
            const codi::RealReverseIndex &qv_prime);

    void set_ref_quants(
            const double qref,
            const double pref,
            const double wref,
            const double tref,
            const double zref,
            const double Nref,
            const double timeref);

    void setup_model_constants(
            const double dt_prime,
            const double uncertainty_perc);

    void setup_model_constants(
            const double uncertainty_perc);

    void setup_model_constants();

    codi::RealReverseIndex::Tape& prepare_call();

    void finish_call(
            codi::RealReverseIndex::Tape& tape,
            double *res,
            double *gradients);
#ifdef CCN_AKM
    void py_ccn_act_hande_akm(
        const double p,
        const double w,
        const double T,
        const double qv,
        const double qc,
        const double Nc,
        double *res,
        double *gradients);
#endif
    void py_graupel_melting(
            const double T,
            const double qg,
            const double Ng,
            double *res,
            double *gradients);

    void py_saturation_adjust(
            const double p,
            const double T,
            const double qv,
            const double qc,
            double *res,
            double *gradients);

    void py_riming_ice(
            const double qc,
            const double Nc,
            const double qi,
            const double Ni,
            const double qr,
            const double Nr,
            const double T,
            double *res,
            double *gradients);

    void py_riming_snow(
            const double qc,
            const double Nc,
            const double qs,
            const double Ns,
            const double qr,
            const double Nr,
            const double T,
            double *res,
            double *gradients);

    void py_riming_with_depo(
            const double qv,
            const double qc,
            const double Nc,
            const double qi,
            const double Ni,
            const double qs,
            const double Ns,
            const double qr,
            const double Nr,
            const double T,
            const double p,
            double *res,
            double *gradients);
};

extern "C" {
physics_t* physics_t_new(char* table_path) {return new physics_t(table_path);}
int physics_t_get_num_comp(physics_t* phys) {return phys->get_num_comp();}
int physics_t_get_num_par(physics_t* phys) {return phys->get_num_par();}
void physics_t_set_ref_quants(
        physics_t* phys,
        const double qref,
        const double pref,
        const double wref,
        const double tref,
        const double zref,
        const double Nref,
        const double timeref) {phys->set_ref_quants(qref, pref, wref, tref, zref, Nref, timeref);}
void physics_t_setup_model_constants(
        physics_t* phys,
        double dt_prime,
        double uncertainty_perc) {phys->setup_model_constants(dt_prime, uncertainty_perc);}
void physics_t_setup_model_constants_uncert(
        physics_t* phys,
        double uncertainty_perc) {phys->setup_model_constants(uncertainty_perc);}
#ifdef CCN_AKM
void physics_t_py_ccn_act_hande_akm(
        physics_t* phys,
        double p,
        double w,
        double T,
        double qv,
        double qc,
        double Nc,
        double *res,
        double *gradients) {phys->py_ccn_act_hande_akm(p, w, T, qv, qc, Nc, res, gradients);}
#endif
void physics_t_py_graupel_melting(
        physics_t* phys,
        double T,
        double qg,
        double Ng,
        double *res,
        double *gradients) {phys->py_graupel_melting(T, qg, Ng, res, gradients);}
void physics_t_py_saturation_adjust(
        physics_t* phys,
        double p,
        double T,
        double qv,
        double qc,
        double *res,
        double *gradients) {phys->py_saturation_adjust(p, T, qv, qc, res, gradients);}
void physics_t_py_riming_ice(
        physics_t* phys,
        double qc,
        double Nc,
        double qi,
        double Ni,
        double qr,
        double Nr,
        double T,
        double *res,
        double *gradients) {phys->py_riming_ice(qc, Nc, qi, Ni, qr, Nr, T, res, gradients);}
void physics_t_py_riming_snow(
        physics_t* phys,
        double qc,
        double Nc,
        double qs,
        double Ns,
        double qr,
        double Nr,
        double T,
        double *res,
        double *gradients) {phys->py_riming_snow(qc, Nc, qs, Ns, qr, Nr, T, res, gradients);}
void physics_t_py_riming_with_depo(
        physics_t* phys,
        double qv,
        double qc,
        double Nc,
        double qi,
        double Ni,
        double qs,
        double Ns,
        double qr,
        double Nr,
        double T,
        double p,
        double *res,
        double *gradients) {phys->py_riming_with_depo(qv, qc, Nc, qi, Ni, qs, Ns, qr, Nr, T, p, res, gradients);}
}

