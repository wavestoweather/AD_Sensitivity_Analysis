#include "include/types/param_t.h"


param_t::param_t() {
    mean            = std::nan("");
    sigma           = std::nan("");
    sigma_perc      = std::nan("");
    err             = 0;
    name            = -1;
    out_name        = -1;
    particle_param  = false;
    func_name       = "";
    positive        = true;
}


param_t::param_t(
    std::string param_type) {
    add_type(param_type);
}


void param_t::add_type(
    std::string param_type) {
    auto it = table_out_param.find(param_type);
    if (it != table_out_param.end()) {
        out_name = static_cast<int>(it->second);
        if (out_name != static_cast<int>(OutParam::model))
            particle_param = true;
        outparam_name = param_type;
    } else {
        SUCCESS_OR_DIE(OUTPARAM_CONFIG_ERR);
    }
}


void param_t::add_mean(
    double m) {
    mean = m;
    if (!isnan(sigma_perc) && isnan(sigma))
        sigma = mean*sigma_perc/100;
}


template<class float_t>
int param_t::add_name(
    std::string n,
    model_constants_t<float_t> &cc) {
    param_name = n;
    if (particle_param) {
        auto it = table_particle_param.find(n);
        if (it != table_particle_param.end()) {
            name = static_cast<uint32_t>(it->second);
            if (std::isnan(mean)) {
                particle_model_constants_t<float_t> *pt_model;
                switch (out_name) {
                    case static_cast<uint32_t>(OutParam::cloud):
                        pt_model = &(cc.cloud);
                        break;
                    case static_cast<uint32_t>(OutParam::rain):
                        pt_model = &(cc.rain);
                        break;
                    case static_cast<uint32_t>(OutParam::snow):
                        pt_model = &(cc.snow);
                        break;
                    case static_cast<uint32_t>(OutParam::graupel):
                        pt_model = &(cc.graupel);
                        break;
                    case static_cast<uint32_t>(OutParam::hail):
                        pt_model = &(cc.hail);
                        break;
                    case static_cast<uint32_t>(OutParam::ice):
                        pt_model = &(cc.ice);
                        break;
                    default:
                        return PARAM_ADD_NAME_ERR;
                }
                mean = pt_model->constants[name].getValue();
            }
        } else {
            err = PARAM_CONFIG_ERR;
        }
    } else {
        auto it = table_param.find(n);
        if (it != table_param.end()) {
            name = static_cast<int>(it->second);
            if (std::isnan(mean))
                mean = get_at(cc.constants, name).getValue();
        } else {
            err = PARAM_CONFIG_ERR;
        }
    }

    if (!isnan(sigma) && err != PARAM_CONFIG_ERR) {
        add_sigma(sigma);
    } else if (!isnan(sigma_perc) && err != PARAM_CONFIG_ERR) {
        add_sigma_perc(sigma_perc);
    }
    return SUCCESS;
}


void param_t::add_sigma(
    double s) {
    sigma = s;
    if (!isnan(mean) && func_name != "")
        add_rand_function(func_name);
}


void param_t::add_sigma_perc(
    double s) {
    sigma_perc = s;
    if (!isnan(mean)) {
        sigma = sigma_perc*mean/100;
        if (func_name != "")
            add_rand_function(func_name);
    }
}


void param_t::add_rand_function(
    std::string name) {
    if (func_name == "")
        func_name = name;
    if (!isnan(sigma_perc) && isnan(sigma) && !isnan(mean))
        sigma = mean*sigma_perc/100;
    if (!isnan(mean) && !isnan(sigma)) {
        if (func_name == "normal") {
            normal_dis = std::normal_distribution<double>(mean, sigma);
            get_rand = std::bind(normal_dis, rand_generator);
        } else if (func_name == "uniform") {
            uniform_dis = std::uniform_real_distribution<double>(mean-sigma, mean+sigma);
            get_rand = std::bind(uniform_dis, rand_generator);
        } else if (func_name == "fixed") {
            // Nothing to be done here
        } else {
        err = DISTRIBUTION_CONFIG_ERR;
        }
    }
}


int param_t::check() {
    if (name == -1) {
        std::cerr << "Error in config file:\n"
                    << "You did not specify the parameter to perturb or "
                    << "you have a typo at <name>typo</name>\n";
        err = MISSING_PARAM_CONFIG_ERR;
        return err;
    }
    if (isnan(sigma) && isnan(sigma_perc)) {
        std::cerr << "Error in config file:\n"
                    << "You did not specify the variance for "
                    << "perturbing the parameter.\n";
        err = MISSING_VARIANCE_CONFIG_ERR;
        return err;
    }
    if (func_name == "")
        err = DISTRIBUTION_CONFIG_ERR;

    switch (err) {
        case PARAM_CONFIG_ERR:
            std::cerr << "Error in config file:\n"
                        << "You used a parameter to perturb that does "
                        << "not exist.\n";
            return err;

        case OUTPARAM_CONFIG_ERR:
            std::cerr << "Error in config file:\n"
                        << "You used an output parameter that does "
                        << "not exist.\n";
            return err;
        case DISTRIBUTION_CONFIG_ERR:
            std::cerr << "Error in config file:\n"
                        << "No such function for generating random "
                        << "numbers in perturbing parameters:\n"
                        << "name: " << func_name << "\n"
                        << "Options are:\n"
                        << "normal: normal distribution with mean=mean "
                        << "or mean=default value of parameter and sigma"
                        << "=sigma or sigma=mean*sigma_perc/100\n"
                        << "uniform: uniform distribution from mean-sigma "
                        << "to mean+sigma and mean, sigma as above\n";
            return err;
        default:
            return err;
    }
}


void to_json(
    nlohmann::json& j,
    const param_t& p) {
    if (p.err != 0)
        return;
    j["mean"] = p.mean;
    j["name"] = p.param_name;
    if (!isnan(p.sigma))
        j["sigma"] = p.sigma;
    else
        j["sigma_perc"] = p.sigma_perc;
    j["rand_func"] = p.func_name;
    j["type"] = p.outparam_name;
}


template<class float_t>
int param_t::from_json(
    const nlohmann::json& j,
    model_constants_t<float_t> &cc) {
    std::string o_name;
    j.at("type").get_to(o_name);
    this->add_type(o_name);
    for (auto &it : j.items()) {
        auto first = it.key();
        if (first == "mean") {
            double m;
            j.at(first).get_to(m);
            this->add_mean(m);
        } else if (first == "name") {
            std::string p_name;
            j.at(first).get_to(p_name);
            err = add_name(p_name, cc);
            if (err != SUCCESS) return err;
        } else if (first == "sigma") {
            double s;
            j.at(first).get_to(s);
            this->add_sigma(s);
        } else if (first == "sigma_perc") {
            double s;
            j.at(first).get_to(s);
            this->add_sigma_perc(s);
        } else if (first == "rand_func") {
            std::string f_name;
            j.at(first).get_to(f_name);
            this->add_rand_function(f_name);
        } else if (first == "type") {
            // Needs to be done before the for-loop
        } else {
            return PARAM_CONFIG_ERR;
        }
    }
    return SUCCESS;
}


template<class float_t>
int param_t::check_name(
    model_constants_t<float_t> &cc) {
    return add_name(param_name, cc);
}


template<class float_t>
void param_t::perturb(
    model_constants_t<float_t> &cc) const {
    if (particle_param) {
        particle_model_constants_t<float_t> *pt_model = nullptr;
        switch (out_name) {
            case static_cast<uint32_t>(OutParam::cloud):
                pt_model = &(cc.cloud);
                break;
            case static_cast<uint32_t>(OutParam::rain):
                pt_model = &(cc.rain);
                break;
            case static_cast<uint32_t>(OutParam::snow):
                pt_model = &(cc.snow);
                break;
            case static_cast<uint32_t>(OutParam::graupel):
                pt_model = &(cc.graupel);
                break;
            case static_cast<uint32_t>(OutParam::hail):
                pt_model = &(cc.hail);
                break;
            case static_cast<uint32_t>(OutParam::ice):
                pt_model = &(cc.ice);
                break;
            default:
                std::cout << "Particle id " << out_name << " does not exist\n";
                SUCCESS_OR_DIE(PERTURB_ERR);
        }
        if (name >= static_cast<int>(pt_model->constants.size())
            || name < 0) {
            std::cerr << "Cannot perturb parameter " << name
            << "/" << static_cast<int>(pt_model->constants.size())
            << " of particle id " << out_name << "\n";
            SUCCESS_OR_DIE(PERTURB_ERR);
        }
        if (func_name == "fixed") {
            if (positive) {
                pt_model->constants[name] = mean+sigma;
            } else {
                pt_model->constants[name] = mean-sigma;
            }
        } else {
            pt_model->constants[name] = get_rand();
        }
        pt_model->perturbed_idx.push_back(name);
    } else {
        if (name >= static_cast<int>(cc.constants.size())
            || name < 0) {
            std::cerr << "Cannot perturb parameter " << name << "/" << static_cast<int>(cc.constants.size()) << "\n";
            SUCCESS_OR_DIE(PERTURB_ERR);
        }
        if (func_name == "fixed") {
            if (positive) {
                cc.constants[name] = mean+sigma;
            } else {
                cc.constants[name] = mean-sigma;
            }
        } else {
            cc.constants[name] = get_rand();
        }
#ifdef DEBUG_SEG
        std::cout << "rank " << cc.rank << " traj " << cc.traj_id << " attempt to add " << name <<
        " with size " << cc.perturbed_idx.size() << "\n"
                      << std::flush;
#endif
        cc.perturbed_idx.push_back(name);
#ifdef DEBUG_SEG
        std::cout << "rank " << cc.rank << " traj " << cc.traj_id << " after add " << name <<
                  " with size " << cc.perturbed_idx.size() << "\n"
                  << std::flush;
#endif
    }
}


template<class float_t>
void param_t::reset(
    model_constants_t<float_t> &cc) const {
    if (particle_param) {
        particle_model_constants_t<float_t> *pt_model = nullptr;
        switch (out_name) {
            case static_cast<uint32_t>(OutParam::cloud):
                pt_model = &(cc.cloud);
                break;
            case static_cast<uint32_t>(OutParam::rain):
                pt_model = &(cc.rain);
                break;
            case static_cast<uint32_t>(OutParam::snow):
                pt_model = &(cc.snow);
                break;
            case static_cast<uint32_t>(OutParam::graupel):
                pt_model = &(cc.graupel);
                break;
            case static_cast<uint32_t>(OutParam::hail):
                pt_model = &(cc.hail);
                break;
            case static_cast<uint32_t>(OutParam::ice):
                pt_model = &(cc.ice);
                break;
            default:
                std::cout << "Error in perturbing...\n";
        }
        pt_model->constants[name] = mean;
        pt_model->perturbed_idx.pop_back();
    } else {
#ifdef DEBUG_SEG
        std::cout << "rank " << cc.rank << " before pop of " << name << " -- "
            << cc.perturbed_idx[cc.perturbed_idx.size()-1] << " size " << cc.perturbed_idx.size() << "\n";
#endif
        cc.constants[name] = mean;
        cc.perturbed_idx.pop_back();
#ifdef DEBUG_SEG
        std::cout << "rank " << cc.rank << " after pop " << cc.perturbed_idx.size() << "\n";
#endif
    }
}


std::string param_t::get_name() const {
    if (particle_param) {
        return (outparam_name + "_" + param_name);
    } else {
        return param_name;
    }
}


template int param_t::check_name<codi::RealReverse>(
    model_constants_t<codi::RealReverse>&);

template int param_t::check_name<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&);

template int param_t::add_name<codi::RealReverse>(
    std::string, model_constants_t<codi::RealReverse>&);

template int param_t::add_name<codi::RealForwardVec<num_par_init> >(
    std::string, model_constants_t<codi::RealForwardVec<num_par_init> >&);

template void param_t::perturb<codi::RealReverse>(
    model_constants_t<codi::RealReverse>&) const;

template void param_t::perturb<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&) const;

template void param_t::reset<codi::RealReverse>(
    model_constants_t<codi::RealReverse>&) const;

template void param_t::reset<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&) const;

template int param_t::from_json<codi::RealReverse>(
    const nlohmann::json&,
    model_constants_t<codi::RealReverse>&);
template int param_t::from_json<codi::RealForwardVec<num_par_init> >(
    const nlohmann::json&,
    model_constants_t<codi::RealForwardVec<num_par_init> >&);
