#pragma once

#include <boost/property_tree/ptree.hpp>
#include <cmath>
#include <functional>
#include <random>
#include <string>
#include <utility>
#include <unordered_map>

// #include "include/microphysics/constants.h"
#include "include/misc/error.h"
#include "include/types/model_constants_t.h"

namespace pt = boost::property_tree;

struct param_t{
    double mean;
    double sigma;
    double sigma_perc;
    int err;
    int name;
    int out_name;
    bool particle_param;

    std::function<double()> get_rand; /*!< distribution used for random number generation. */
    std::normal_distribution<double> normal_dis;
    std::uniform_real_distribution<double> uniform_dis;
    std::string param_name;
    std::string outparam_name;
    std::string func_name;

    enum class OutParam: uint32_t {
        model, cloud, rain, ice, graupel, hail, snow
    };
    std::unordered_map<std::string, OutParam> const table_out_param = {
        {"model", OutParam::model},
        {"cloud", OutParam::cloud}, {"rain", OutParam::rain},
        {"ice", OutParam::ice}, {"graupel", OutParam::graupel},
        {"hail", OutParam::hail}, {"snow", OutParam::snow}
    };

    param_t();

    param_t(std::string param_type);

    void add_type(std::string param_type);

    void add_mean(double m);

    void add_name(std::string n, model_constants_t &cc);

    void add_sigma(double s);

    void add_sigma_perc(double s);

    void add_rand_function(std::string name);

    int check();

    void put(pt::ptree &ptree) const;

    int from_pt(pt::ptree &ptree, model_constants_t &cc);

    void perturb(model_constants_t &cc) const;

    std::string get_name() const;
};