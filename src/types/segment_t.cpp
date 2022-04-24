#include "include/types/segment_t.h"


segment_t::segment_t() {
    value           = std::nan("");
    old_value       = std::nan("");
    n_members       = 1;
    value_name      = -1;
    value_name_sig  = -1;
    out_param       = -1;
    n_segments      = 1;
    err             = 0;
    old_sign        = 0;
    duration        = 0;
    activated       = false;
    method          = value_method;
}


void segment_t::add_param(
    param_t &param) {
    params.push_back(param);
}


void segment_t::add_value(
    double v) {
    value = v;
}

            // this->params = j["params"].get<std::vector<param_t> >();;
            // j.at("params").get_to(this->params);
void segment_t::add_value_name(
    std::string n) {
    if (method == impact_change) {
        value_name = -1;
        return;
    }
    auto it = table_param.find(n);
    if (it != table_param.end()) {
        value_name = it->second;
        std::string key = "when_name";
        tree_strings.insert(std::make_pair(key, n));
    } else {
        err = VALUE_NAME_CONFIG_ERR;
    }
}


void segment_t::add_method(
    std::string m) {
    auto it = table_method.find(m);
    if (it != table_method.end()) {
        method = it->second;
        std::string key = "when_method";
        tree_strings.insert(std::make_pair(key, m));
        if (method == impact_change)
            value_name = -1;
    } else {
        err = METHOD_CONFIG_ERR;
    }
}


void segment_t::add_counter(
    uint32_t c) {
    n_segments = c;
}


void segment_t::add_duration(
    double t) {
    duration = t;
}


void segment_t::add_out_param(
    std::string p) {
    auto it = std::find(output_par_idx.begin(), output_par_idx.end(), p);
    if (it != output_par_idx.end()) {
        out_param = std::distance(output_par_idx.begin(), it);
        std::string key = "when_sens";
        tree_strings.insert(std::make_pair(key, p));
    } else {
        err = OUTPARAM_CONFIG_ERR;
    }
}


void segment_t::add_amount(
    int a) {
    if (a > 1)
        n_members = a;
    else
        err = N_MEMBERS_CONFIG_ERR;
}


int segment_t::check() {
    if (n_members == 1 && method != Method::full_perturbation)
        err = N_MEMBERS_CONFIG_ERR;
    else if (n_segments < 0 && method != Method::full_perturbation)
        err = N_SEGMENTS_CONFIG_ERR;
    else if (value_name == -1
            && method != Method::impact_change
            && method != Method::repeated_time
            && method != Method::full_perturbation)
        err = VALUE_NAME_CONFIG_ERR;
    else if (method == -1)
        err = METHOD_CONFIG_ERR;
    else if (out_param == -1 && value_name >= Param::rain_a_geo)
        err = OUTPARAM_CONFIG_ERR;
    switch (err) {
        case VALUE_NAME_CONFIG_ERR:
            std::cerr << "Error in config file:\n"
                        << "The name of the parameter used to determine "
                        << "when to start an ensemble is not defined or "
                        << "does not exist.\n";
            n_segments = 0;
            return err;
        case METHOD_CONFIG_ERR:
            std::cerr << "Error in config file:\n"
                        << "The name of the method used to determine "
                        << "when to start an ensemble is not defined or "
                        << "does not exist.\n";
            n_segments = 0;
            return err;
        case OUTPARAM_CONFIG_ERR:
            std::cerr << "Error in config file:\n"
                        << "The name of the output parameter used to determine "
                        << "when to start an ensemble is not defined or "
                        << "does not exist.\n";
            n_segments = 0;
            return err;
        case N_MEMBERS_CONFIG_ERR:
            std::cerr << "Error in config file:\n"
                        << "The number of members for an ensemble must be "
                        << "2 or higher. One simulation always runs without "
                        << "perturbed parameters for comparison.\n";
            n_segments = 0;
            return err;
        case N_SEGMENTS_CONFIG_ERR:
            std::cerr << "Error in config file:\n"
                        << "The number of segments must be at least 1 in "
                        << "order to start at least 1 ensemble.\n";
            n_segments = 0;
            return err;
    }
    if (err != 0)
        n_segments = 0;
    return err;
}


template<class float_t>
bool segment_t::perturb_check(
//    const model_constants_t<float_t> &cc,
    const std::vector< std::array<double, num_par > > &gradients,
    const std::vector<float_t> &y,
    const double timestep) {
    if (n_segments == 0 || method == full_perturbation)
        return false;
    int idx;
    switch (method) {
        case impact_change: {
            auto it_minmax = std::minmax_element(gradients[out_param].begin(),
                gradients[out_param].end());
            auto it_max = (*it_minmax.second > abs(*it_minmax.first)) ? it_minmax.second : it_minmax.first;
            idx = std::distance(gradients[out_param].begin(), it_max);
            if (value_name_sig != -1) {
                if (idx != value_name_sig) {
                    value_name_sig = idx;
                    activated = true;
                    return true;
                }
                return false;
            } else {
                // First time we check the highest impact
                value_name_sig = idx;
                return false;
            }
            return false;
        }
        case sign_flip: {
            // the first num_comp many values refer to output parameters
            idx = value_name - num_comp;
            if (old_sign == 0) {
                // set sign; no perturbing needed.
                if (gradients[out_param][idx] == 0)
                    return false;
                old_sign = (gradients[out_param][idx] > 0) ? 1 : 2;
                return false;
            } else {
                // check if a sign change happend
                if (old_sign == 1 && gradients[out_param][idx] < 0) {
                    activated = true;
                    old_sign = 2;
                    return true;
                } else if (old_sign == 2 && gradients[out_param][idx] > 0) {
                    // Perturb parameters
                    activated = true;
                    old_sign = 1;
                    return true;
                }
                return false;
            }
            return false;
        }
        case value_method: {
            double current_value;
            // If value is an output parameter or a sensitivity
            if (out_param != -1) {
                current_value = gradients[out_param][value_name - num_comp];
            } else {
                current_value = y[value_name].getValue();
                if (value_name == T_idx)
                    current_value *= 273.15;
                else if (value_name == p_idx)
                    current_value *= 1.0e5;
            }
            if (current_value == value) {
                activated = true;
                return true;
            }

            if (std::isnan(old_value)) {
                old_value = current_value;
                return false;
            }

            if (signbit(value-current_value) != signbit(value-old_value)
                || fabs(value-current_value) < fabs(value*tol)) {
                activated = true;
                return true;
            }
            old_value = current_value;
            return false;
        case repeated_time:
            if (std::fmod(timestep, value) == 0) {
                activated = true;
                return true;
            }
            return false;
        }
        default: {
            return false;
        }
    }
    return false;
}


void segment_t::deactivate(
    const bool keep_repeated) {
    if (!activated) return;
    // Perturbing every few timesteps is done indefenitely
    // unless a fixed duration is given at which the ensembles should stop
    if ((method != repeated_time && n_segments > 0)
        || (method == repeated_time && n_segments > 0 && duration != 0 && !keep_repeated))
        n_segments--;
    activated = false;
}


template<class float_t>
void segment_t::perturb(
    model_constants_t<float_t> &cc,
    const reference_quantities_t &ref_quant,
    input_parameters_t &input,
    std::string &descr) {
    // Sanity check if had been done already
    if (n_segments == 0 && method != full_perturbation)
        return;
    // Change the number of time steps if a fixed duration is given
    if (duration != 0) {
        if (duration + input.current_time + input.start_time < cc.t_end_prime) {
            cc.t_end_prime = duration + input.current_time + input.start_time;
            cc.t_end = cc.t_end_prime/ref_quant.tref;
            input.t_end_prime = duration;
        }
    }
    // Perturb every param
    for (auto &p : params) {
        p.positive = cc.traj_id%2;
        p.perturb(cc);
        descr += p.get_name() + " ";
    }
    // When perturbing is done, deactivate
    if (method != full_perturbation)
        deactivate();
}


template<class float_t>
void segment_t::reset_variables(
    model_constants_t<float_t> &cc) {
    n_segments++;
    // Change the previously perturbed values to the original ones
    for (auto &p : params) {
        p.reset(cc);
    }
}

void to_json(
    nlohmann::json& j,
    const segment_t& s) {
    if ((s.err != 0 || s.n_segments < 1) && !s.activated)
        return;
    if (!isnan(s.value)) {
        j["when_value"] = s.value;
    }
    if (s.value_name != -1) {
        j["when_name"] = s.tree_strings.find("when_name")->second;
    }
    if (s.out_param != -1) {
        j["when_sens"] = s.tree_strings.find("when_sens")->second;
    }
    if (s.n_segments != 1) {
        j["when_counter"] = s.n_segments;
    }
    if (s.method != s.value_method) {
        j["when_method"] = s.tree_strings.find("when_method")->second;
    }
    if (s.n_members != 1) {
        j["amount"] = s.n_members;
    }
    if (s.activated) {
        j["activated"] = true;
    }
    if (s.duration > 0) {
        j["duration"] = s.duration;
    }
    j["params"] = s.params;
}


template<class float_t>
void segment_t::from_json(
    const nlohmann::json& j,
    model_constants_t<float_t> &cc) {

    for (auto &it : j.items()) {
        auto first = it.key();
        if (first == "when_value") {
            this->add_value(static_cast<double>(it.value()));
        } else if (first == "when_name") {
            this->add_value_name(static_cast<std::string>(it.value()));
        } else if (first == "amount") {
            this->add_amount(static_cast<uint32_t>(it.value()));
        } else if (first == "when_method") {
            this->add_method(static_cast<std::string>(it.value()));
        } else if (first == "when_counter") {
            this->add_counter(static_cast<uint32_t>(it.value()));
        } else if (first == "when_sens") {
            this->add_out_param(static_cast<std::string>(it.value()));
        } else if (first == "activated") {
            this->activated = static_cast<bool>(it.value());
        } else if (first == "duration") {
            this->add_duration(static_cast<double>(it.value()));
        } else if (first == "params") {
            this->params.clear();
            for (const auto &p_config : j["params"]) {
                param_t p;
                SUCCESS_OR_DIE(p.from_json(p_config, cc));
                this->params.push_back(p);
            }
            // for (auto &param : this->params) SUCCESS_OR_DIE(param.check_name(cc));
        } else {
            SUCCESS_OR_DIE(SEGMENTS_CHECKPOINT_ERR);
        }
    }
}

double segment_t::limit_duration() const {
    if (method == repeated_time) {
        return value;
    }
    return 0;
}

template bool segment_t::perturb_check<codi::RealReverse>(
//    const model_constants_t<codi::RealReverse>&,
    const std::vector< std::array<double, num_par > >&,
    const std::vector<codi::RealReverse>&,
    const double);

template bool segment_t::perturb_check<codi::RealForwardVec<num_par_init> >(
//    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector< std::array<double, num_par > >&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    const double);

template void segment_t::perturb<codi::RealReverse>(
    model_constants_t<codi::RealReverse>&,
    const reference_quantities_t&,
    input_parameters_t&,
    std::string&);

template void segment_t::perturb<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const reference_quantities_t&,
    input_parameters_t&,
    std::string&);

template void segment_t::reset_variables<codi::RealReverse>(
    model_constants_t<codi::RealReverse>&);

template void segment_t::reset_variables<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&);

template void segment_t::from_json<codi::RealReverse>(
    const nlohmann::json&,
    model_constants_t<codi::RealReverse>&);

template void segment_t::from_json<codi::RealForwardVec<num_par_init> >(
    const nlohmann::json&,
    model_constants_t<codi::RealForwardVec<num_par_init> >&);
