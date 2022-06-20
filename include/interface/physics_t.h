#pragma once

#include <include/microphysics/constants.h>
#include <include/microphysics/user_functions.h>
#include <include/microphysics/program_io.h>

#include <iostream>
#include <vector>

#include <fstream>


struct physics_t {
 private:
    model_constants_t<codi::RealReverse> cc;
    input_parameters_t input;
    reference_quantities_t ref_quant;
    std::vector<codi::RealReverse> yold;
    std::vector<codi::RealReverse> ynew;
    std::vector<codi::RealReverse> ytmp;
    std::vector<codi::RealReverse> k;
    std::vector< std::array<double, num_par > > gradients_;
    double uncertainty_scale;

 public:
    physics_t(std::string table_path="dmin_wetgrowth_lookup.dat");

    void set_limits(
        std::vector<codi::RealReverse> &y,
        model_constants_t<codi::RealReverse> &cc);

    int get_num_comp() {return num_comp;}
    int get_num_par() {return num_par;}

    void compute_nondimensional_effects(
        std::vector<codi::RealReverse> &res,
        const codi::RealReverse &p,
        const codi::RealReverse &w,
        const codi::RealReverse &w_prime,
        const codi::RealReverse &T,
        const codi::RealReverse &qv,
        const codi::RealReverse &qv_prime);

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

    void py_ccn_act_hande_akm(
        const double p,
        const double w,
        const double T,
        const double qv,
        const double qc,
        const double Nc,
        double *res,
        double *gradients);

    void py_graupel_melting(
        const double T,
        const double qg,
        const double Ng,
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
    void physics_t_py_graupel_melting(
        physics_t* phys,
        double T,
        double qg,
        double Ng,
        double *res,
        double *gradients) {phys->py_graupel_melting(T, qg, Ng, res, gradients);}
}

