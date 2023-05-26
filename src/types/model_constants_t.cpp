#include "include/types/model_constants_t.h"

template<class float_t>
model_constants_t<float_t>::model_constants_t(
    const std::string &tracking_filename,
    const bool &track_initial_cond) {
    id = "0";
    ensemble_id = 0;
    traj_id = 0;
    n_trajs = 1;
    n_ensembles = 1;
    ens_desc = "root ";
    done_steps = 0;
    std::fill(constants.begin(), constants.end(), 0);
    std::fill(initial_conditions.begin(), initial_conditions.end(), 0);
    int track_amount = num_par / 64 + (num_par%64 != 0);
    track_param.resize(track_amount);
    track_state = 0;
    for (auto &t : track_param)
        t = 0;
    if (tracking_filename.compare("") != 0) {
        local_num_comp = 0;
        local_num_par = 0;
        local_ic_par = 0;
        this->load_configuration(tracking_filename);
    } else {
        local_num_comp = num_comp;
        local_num_par = num_par;
        track_state = -1;

        if ( track_initial_cond ) {
            track_ic = -1;
            local_ic_par = static_cast<int>(Init_cons_idx::n_items);
            // Tracking initial conditions and model parameters is not
            // possible. Hence we set local_num_par.
            local_num_par = local_ic_par;
        } else {
            track_ic = 0;
            local_ic_par = 0;
        }
        for (auto &t : track_param)
            t = -1;
    }
}


template<>
void model_constants_t<codi::RealForwardVec<num_par_init> >::register_input() {
    for (int i=0; i<static_cast<int>(Init_cons_idx::n_items); i++)
        if (trace_check(i, 2))
            this->initial_conditions[i].gradient()[i] = 1;
}


template<>
void model_constants_t<codi::RealReverse>::register_input() {
    // Nothing to do for the backward mode
}


template<>
void model_constants_t<codi::RealReverse>::register_input(
    codi::RealReverse::Tape &tape) {
#ifdef DEVELOP
    std::cout << "register input of at most " << static_cast<int>(Cons_idx::n_items)
              << " vs size " << this->constants.size() << "\n";
#endif
    for (int i=0; i<static_cast<int>(Cons_idx::n_items); i++) {
        if (trace_check(i, false)) {
#ifdef DEVELOP
            std::cout << "register input " << i << "\n";
#endif
            tape.registerInput(this->constants[i]);
        }
    }
    uint32_t offset = static_cast<uint32_t>(Cons_idx::n_items);
    if (local_num_par == num_par) {
        this->rain.register_input(tape, offset);
        this->cloud.register_input(tape, offset);
#if defined(RK4ICE)
        this->graupel.register_input(tape, offset);
        this->hail.register_input(tape, offset);
        this->ice.register_input(tape, offset);
        this->snow.register_input(tape, offset);
#endif
    } else {
        for (auto &c : this->rain.constants) {
            if (trace_check(offset, false))
                tape.registerInput(c);
            offset++;
        }
        for (auto &c : this->cloud.constants) {
            if (trace_check(offset, false)) {
                tape.registerInput(c);
            }
            offset++;
        }
#if defined(RK4ICE)
        for (auto &c : this->graupel.constants) {
            if (trace_check(offset, false))
                tape.registerInput(c);
            offset++;
        }
        for (auto &c : this->hail.constants) {
            if (trace_check(offset, false))
                tape.registerInput(c);
            offset++;
        }
        for (auto &c : this->ice.constants) {
            if (trace_check(offset, false))
                tape.registerInput(c);
            offset++;
        }
        for (auto &c : this->snow.constants) {
            if (trace_check(offset, false))
                tape.registerInput(c);
            offset++;
        }
#endif
    }
}


template<>
void model_constants_t<codi::RealForwardVec<num_par_init> >::get_gradient(
    std::array<double, num_par> &out_vec,
    std::vector<codi::RealForwardVec<num_par_init> > &y_single_new,
    uint32_t ii,
    const double &ref_value) const {
#ifdef DEVELOP
    std::cout << "get_gradient real_forward\n";
#endif
    uint32_t offset = out_vec.size() - static_cast<int>(Init_cons_idx::n_items);
    for (int i=0; i<static_cast<int>(Init_cons_idx::n_items); ++i) {
        if (trace_check(i, 2)) {
            out_vec[i + offset] = y_single_new[ii].getGradient()[i]
                * uncertainty[i + static_cast<uint32_t>(Cons_idx::n_items)] * ref_value;
        }
    }
}


template<>
void model_constants_t<codi::RealReverse>::get_gradient(
    std::array<double, num_par> &out_vec,
    const double &ref_value) const {
#ifdef DEVELOP
    std::cout << "get_gradient items = " << static_cast<int>(Cons_idx::n_items)
              << " out_vec: " << out_vec.size() << " constants " << this->constants.size()
              << " uncertainty: " << uncertainty.size() << "\n";
#endif
    for (int i=0; i<static_cast<int>(Cons_idx::n_items); ++i)
        if (trace_check(i, false)) {
            out_vec[i] = this->constants[i].getGradient() * uncertainty[i] * ref_value;
        }

    uint32_t idx = static_cast<uint32_t>(Cons_idx::n_items);
#ifdef DEVELOP
    std::cout << "get_gradient idx = " << idx << "\n";
#endif
    if (local_num_par == num_par) {
        this->rain.get_gradient(out_vec, idx, ref_value);
        this->cloud.get_gradient(out_vec, idx, ref_value);
#if defined(RK4ICE)
        this->graupel.get_gradient(out_vec, idx, ref_value);
        this->hail.get_gradient(out_vec, idx, ref_value);
        this->ice.get_gradient(out_vec, idx, ref_value);
        this->snow.get_gradient(out_vec, idx, ref_value);
#endif
#ifdef DEVELOP
        std::cout << "local_num_par == num_par\n";
#endif
    } else {
        uint32_t start_idx = idx;
        for (auto &c : this->rain.constants) {
            if (trace_check(idx, false))
                out_vec[idx] = c.getGradient() * this->rain.uncertainty[idx-start_idx] * ref_value;
            idx++;
        }
        start_idx = idx;
        for (auto &c : this->cloud.constants) {
            if (trace_check(idx, false))
                out_vec[idx] = c.getGradient() * this->cloud.uncertainty[idx-start_idx] * ref_value;
            idx++;
        }
#if defined(RK4ICE)
        start_idx = idx;
        for (auto &c : this->graupel.constants) {
            if (trace_check(idx, false))
                out_vec[idx] = c.getGradient() * this->graupel.uncertainty[idx-start_idx] * ref_value;
            idx++;
        }
        start_idx = idx;
        for (auto &c : this->hail.constants) {
            if (trace_check(idx, false))
                out_vec[idx] = c.getGradient() * this->hail.uncertainty[idx-start_idx] * ref_value;
            idx++;
        }
        start_idx = idx;
        for (auto &c : this->ice.constants) {
            if (trace_check(idx, false))
                out_vec[idx] = c.getGradient() * this->ice.uncertainty[idx-start_idx] * ref_value;
            idx++;
        }
        start_idx = idx;
        for (auto &c : this->snow.constants) {
            if (trace_check(idx, false))
                out_vec[idx] = c.getGradient() * this->snow.uncertainty[idx-start_idx] * ref_value;
            idx++;
        }
#endif
    }
}


template<>
void model_constants_t<codi::RealForwardVec<num_par_init> >::get_gradients(
    std::vector<codi::RealForwardVec<num_par_init> > &y_single_new,
    std::vector< std::array<double, num_par > > &y_diff,
    const reference_quantities_t &ref_quant) const {

    for (uint32_t ii = 0 ; ii < num_comp ; ii++) {
        if (trace_check(ii, true)) {
            double ref_value = 1.0;
            switch (ii) {
                case p_idx:
                    ref_value = ref_quant.pref;
                    break;
                case T_idx:
                    ref_value = ref_quant.Tref;
                    break;
                case w_idx:
                    ref_value = ref_quant.wref;
                    break;
                case qv_idx:
                case qc_idx:
                case qr_idx:
                case qs_idx:
                case qi_idx:
                case qg_idx:
                case qh_idx:
                case qr_out_idx:
                case qi_out_idx:
                case qs_out_idx:
                case qg_out_idx:
                case qh_out_idx:
                    ref_value = ref_quant.qref;
                    break;
                case Nc_idx:
                case Nr_idx:
                case Ns_idx:
                case Ni_idx:
                case Ng_idx:
                case Nh_idx:
                case Nr_out_idx:
                case Ns_out_idx:
                case Ni_out_idx:
                case Ng_out_idx:
                case Nh_out_idx:
                case depo_idx:
                case sub_idx:
                    ref_value = ref_quant.Nref;
                    break;
                case z_idx:
                    ref_value = ref_quant.zref;
                    break;
                default:
                    ref_value = 1.0;
            }
            this->get_gradient(y_diff[ii], y_single_new, ii, ref_value);
        }
    }
}


template<>
void model_constants_t<codi::RealReverse>::get_gradients(
    std::vector<codi::RealReverse> &y_single_new,
    std::vector< std::array<double, num_par > > &y_diff,
    codi::RealReverse::Tape &tape,
    const reference_quantities_t &ref_quant)  {

    for (uint32_t ii = 0 ; ii < num_comp ; ii++) {
        if (trace_check(ii, true)) {
#ifdef DEVELOP
        std::cout << "register Output " << ii << "\n";
#endif
            tape.registerOutput(y_single_new[ii]);
        }
    }
#ifdef DEVELOP
    std::cout << "get_gradients after register, size of y_single_new: " << y_single_new.size()
              <<  ", y_diff: " << y_diff.size() << "\n";
    std::cout << "num_comp " << num_comp << " vs Cons_idx " <<  static_cast<int>(Cons_idx::n_items) << "\n";
#endif
    tape.setPassive();
    for (uint32_t ii = 0 ; ii < num_comp ; ii++) {
#ifdef DEVELOP
        std::cout << "get_gradients ii = " << ii << "\n";
#endif
        if (!trace_check(ii, true))
            continue;
#ifdef DEVELOP
        std::cout << "get_gradients after trace_check ii = " << ii
                  << ", y_diff[ii] " << y_diff[ii].size()
                  << ", single_new " << y_single_new.size() << "\n";
#endif
        y_single_new[ii].setGradient(1.0);
#ifdef DEVELOP
        std::cout << "get_gradients before tape evaluate\n";
#endif
        tape.evaluate();
#ifdef DEVELOP
        std::cout << "get_gradients after tape evaluate\n" << std::flush;
#endif
        double ref_value = 1.0;
        switch (ii) {
            case p_idx:
                ref_value = ref_quant.pref;
                break;
            case T_idx:
                ref_value = ref_quant.Tref;
                break;
            case w_idx:
                ref_value = ref_quant.wref;
                break;
            case qv_idx:
            case qc_idx:
            case qr_idx:
            case qs_idx:
            case qi_idx:
            case qg_idx:
            case qh_idx:
            case qr_out_idx:
            case qi_out_idx:
            case qs_out_idx:
            case qg_out_idx:
            case qh_out_idx:
                ref_value = ref_quant.qref;
                break;
            case Nc_idx:
            case Nr_idx:
            case Ns_idx:
            case Ni_idx:
            case Ng_idx:
            case Nh_idx:
            case Nr_out_idx:
            case Ns_out_idx:
            case Ni_out_idx:
            case Ng_out_idx:
            case Nh_out_idx:
            case depo_idx:
            case sub_idx:
                ref_value = ref_quant.Nref;
                break;
            case z_idx:
                ref_value = ref_quant.zref;
                break;
            default:
                ref_value = 1.0;
        }
        this->get_gradient(y_diff[ii], ref_value);
        tape.clearAdjoints();
    }

#ifdef DEVELOP
    std::cout << "get_gradients end\n";
#endif
    tape.reset();
}


template<class float_t>
void to_json(
    nlohmann::json& j,
    const model_constants_t<float_t>& cc) {
    j["id"] = cc.id;
    j["ensemble_id"] = cc.ensemble_id;
    j["traj_id"] = cc.traj_id;
    j["n_trajs"] = cc.n_trajs;
    j["ens_desc"] = cc.ens_desc;
    j["num_steps"] = cc.num_steps;
    j["n_ensembles"] = cc.n_ensembles;

    // technical parameters
    j["t_end_prime"] = cc.t_end_prime;
    j["t_end"] = cc.t_end;
    j["dt_prime"] = cc.dt_prime;
    j["dt"] = cc.dt;
    j["dt_traject_prime"] = cc.dt_traject_prime;
    j["dt_traject"] = cc.dt_traject;
    j["num_steps"] = cc.num_steps;
    j["num_sub_steps"] = cc.num_sub_steps;
    j["done_steps"] = cc.done_steps;
    j["dt_half"] = cc.dt_half;
    j["dt_sixth"] = cc.dt_sixth;
    j["dt_third"] = cc.dt_third;
    j["local_num_comp"] = cc.local_num_comp;
    j["local_num_par"] = cc.local_num_par;
    j["track_state"] = cc.track_state;
    j["checkpoint_steps"] = cc.checkpoint_steps;
    j["track_param"] = cc.track_param;

    if (!cc.perturbed_idx.empty()) {
        std::map<uint32_t, double> perturbed;
        for (uint32_t idx : cc.perturbed_idx) {
            perturbed[idx] = cc.constants[idx].getValue();
        }
        j["perturbed"] = perturbed;
    }
    j["hail"] = cc.hail;
    j["ice"] = cc.ice;
    j["snow"] = cc.snow;
    j["cloud"] = cc.cloud;
    j["rain"] = cc.rain;
    j["graupel"] = cc.graupel;
    // collection coefficients depend completely on particle
    // model constants and can be derived from there. So
    // we skip adding that to the checkpoint.
}


template<class float_t>
int model_constants_t<float_t>::from_json(
    const nlohmann::json &j) {
    int err = 0;
    for (auto &it : j.items()) {
        auto first = it.key();
        if (first == "id") {
            std::string new_id;
            j.at(first).get_to(new_id);
            this->id = new_id + "-" + id;
        } else if (first == "ensemble_id") {
            j.at(first).get_to(this->ensemble_id);
        } else if (first == "traj_id") {
            j.at(first).get_to(this->traj_id);
        } else if (first == "n_trajs") {
            j.at(first).get_to(this->n_trajs);
        } else if (first == "ens_desc") {
            j.at(first).get_to(this->ens_desc);
        } else if (first == "t_end_prime") {
            j.at(first).get_to(this->t_end_prime);
        } else if (first == "t_end") {
            j.at(first).get_to(this->t_end);
        } else if (first == "dt_prime") {
            j.at(first).get_to(this->dt_prime);
        } else if (first == "dt") {
            j.at(first).get_to(this->dt);
        } else if (first == "dt_traject_prime") {
            j.at(first).get_to(this->dt_traject_prime);
        } else if (first == "dt_traject") {
            j.at(first).get_to(this->dt_traject);
        } else if (first == "num_steps") {
            j.at(first).get_to(this->num_steps);
        } else if (first == "num_sub_steps") {
            j.at(first).get_to(this->num_sub_steps);
        } else if (first == "checkpoint_steps") {
            j.at(first).get_to(this->checkpoint_steps);
        } else if (first == "dt_half") {
            j.at(first).get_to(this->dt_half);
        } else if (first == "dt_sixth") {
            j.at(first).get_to(this->dt_sixth);
        } else if (first == "dt_third") {
            j.at(first).get_to(this->dt_third);
        } else if (first == "n_ensembles") {
            j.at(first).get_to(this->n_ensembles);
        } else if (first == "num_steps") {
            j.at(first).get_to(this->num_steps);
        } else if (first == "done_steps") {
            j.at(first).get_to(this->done_steps);
        } else if (first == "local_num_comp") {
            j.at(first).get_to(this->local_num_comp);
        } else if (first == "local_num_par") {
            j.at(first).get_to(this->local_num_par);
        } else if (first == "track_state") {
            j.at(first).get_to(this->track_state);
        } else if (first == "track_param") {
            this->track_param.clear();
            j.at(first).get_to(this->track_param);
        } else if (first == "perturbed") {
            this->perturbed_idx.clear();
            std::map<uint32_t, double> perturbed;
            j.at(first).get_to(perturbed);
            for (auto const &p : perturbed) {
                if (p.first >= this->constants.size())
                    err = MODEL_CONS_CHECKPOINT_ERR;
                this->perturbed_idx.push_back(p.first);
                this->constants[p.first] = p.second;
            }
        // below from here: perturbed particle models
        } else if (first == "hail") {
            err = this->hail.from_json(j.at(first));
        } else if (first == "ice") {
            err = this->ice.from_json(j.at(first));
        } else if (first == "snow") {
            err = this->snow.from_json(j.at(first));
        } else if (first == "cloud") {
            err = this->cloud.from_json(j.at(first));
        } else if (first == "rain") {
            err = this->rain.from_json(j.at(first));
        } else if (first == "graupel") {
            err = this->graupel.from_json(j.at(first));
        } else {
            std::cout << "Got '" << first << "' from property tree "
            << "which does not exist in model_constants_t.\n";
            err = MODEL_CONS_CHECKPOINT_ERR;
        }
    }
    return err;
}


#if defined(RK4_ONE_MOMENT)
template<class float_t>
void model_constants_t<float_t>::setCoefficients(
    std::vector<codi::RealReverse> & y,
    reference_quantities_t& ref) {
    float_t p_prime = y[p_idx]*ref.pref;
    float_t T_prime = y[T_idx]*ref.Tref;

    float_t rho_prime = p_prime /(get_at(this->constants, Cons_idx::R_a) * T_prime);
    float_t L_vap_prime = latent_heat_water_supercooled(T_prime, get_at(this->constants, Cons_idx::M_w));
    float_t Ka_prime = thermal_conductivity_dry_air(T_prime);
    float_t psat_prime = saturation_pressure_water(
        T_prime,
        get_at(this->constants, Cons_idx::p_sat_low_temp),
        get_at(this->constants, Cons_idx::p_sat_const_a),
        get_at(this->constants, Cons_idx::T_sat_low_temp),
        get_at(this->constants, Cons_idx::p_sat_const_b));
    float_t A_pp = (L_vap_prime/(Ka_prime*T_prime))
        * ((L_vap_prime/(get_at(this->constants, Cons_idx::R_v)*T_prime)) - 1.0);
    float_t B_pp = (get_at(this->constants, Cons_idx::R_v)*T_prime)/((2.21/p_prime)*psat_prime);


    this->constants[static_cast<int>(Cons_idx::e1_prime)] = this->e1_scale
        * (pow(rho_prime, 2.0*this->alpha_r-2.0)/(A_pp + B_pp));
    this->constants[static_cast<int>(Cons_idx::e2_prime)] = this->e2_scale
        * (pow(rho_prime, this->alpha_r*this->epsilonr - (7.0/4.0))/(A_pp + B_pp));

    this->constants[static_cast<int>(Cons_idx::a1_prime)] = this->a1_scale;  // Constant coefficient
    this->constants[static_cast<int>(Cons_idx::a2_prime)] = this->a2_scale;  // Constant coefficient
    this->constants[static_cast<int>(Cons_idx::d_prime)] = this->d_scale;  // Constant coefficient
}


template<class float_t>
void model_constants_t<float_t>::setCoefficients(
    float_t p_prime,
    float_t T_prime) {
    float_t rho_prime = p_prime /(get_at(this->constants, Cons_idx::R_a) * T_prime);
    float_t L_vap_prime = latent_heat_water_supercooled(T_prime, get_at(this->constants, Cons_idx::M_w));
    float_t Ka_prime = thermal_conductivity_dry_air(T_prime);
    float_t psat_prime = saturation_pressure_water(
        T_prime,
        get_at(this->constants, Cons_idx::p_sat_low_temp),
        get_at(this->constants, Cons_idx::p_sat_const_a),
        get_at(this->constants, Cons_idx::T_sat_low_temp),
        get_at(this->constants, Cons_idx::p_sat_const_b));
    float_t A_pp =
        (L_vap_prime/(Ka_prime*T_prime))*((L_vap_prime/(get_at(this->constants, Cons_idx::R_v)*T_prime)) - 1.0);
    float_t B_pp = (get_at(this->constants, Cons_idx::R_v)*T_prime)/((2.21/p_prime)*psat_prime);

    this->constants[static_cast<int>(Cons_idx::a1_prime)] = this->a1_scale;  // Constant coefficient
    this->constants[static_cast<int>(Cons_idx::a2_prime)] = this->a2_scale;  // Constant coefficient
    this->constants[static_cast<int>(Cons_idx::e1_prime)] =
        this->e1_scale * (pow(rho_prime, 2.0*this->alpha_r-2.0)/(A_pp + B_pp));
    this->constants[static_cast<int>(Cons_idx::e2_prime)] =
        this->e2_scale * (pow(rho_prime, this->alpha_r*this->epsilonr - (7.0/4.0))/(A_pp + B_pp));
    this->constants[static_cast<int>(Cons_idx::d_prime)] = this->d_scale;  // Constant coefficient
}
#endif

template<class float_t>
void model_constants_t<float_t>::setup_cloud_autoconversion(
    particle_model_constants_t<float_t> &pc) {
    float_t nu_local = get_at(pc.constants, Particle_cons_idx::nu) + 1.0;
    float_t mu_local = get_at(pc.constants, Particle_cons_idx::mu);
    if (get_at(pc.constants, Particle_cons_idx::mu) == 1.0) {
        this->constants[static_cast<int>(Cons_idx::cloud_k_au)] =
            get_at(this->constants, Cons_idx::kc_autocon)
            / get_at(pc.constants, Particle_cons_idx::max_x) * 0.05
            * (nu_local+1.0)*(nu_local+3.0) / pow(nu_local, 2);
        this->constants[static_cast<int>(Cons_idx::cloud_k_sc)] =
            get_at(this->constants, Cons_idx::kc_autocon) * (nu_local+1.0)/(nu_local);
    } else {
        this->constants[static_cast<int>(Cons_idx::cloud_k_au)] =
            get_at(this->constants, Cons_idx::kc_autocon)
            / get_at(pc.constants, Particle_cons_idx::max_x) * 0.05
            * (2.0 * tgamma((nu_local+3.0)/mu_local)
            * tgamma((nu_local+1.0)/mu_local) * pow(tgamma((nu_local)/mu_local), 2)
            - 1.0 * pow(tgamma((nu_local+2.0)/mu_local), 2) * pow(tgamma((nu_local)/mu_local), 2))
            / pow(tgamma((nu_local+1.0)/mu_local), 4);
        this->constants[static_cast<int>(Cons_idx::cloud_k_sc)] =
            get_at(this->constants, Cons_idx::kc_autocon)
            * get_at(pc.constants, Particle_cons_idx::c_z);
    }
}


template<class float_t>
void model_constants_t<float_t>::setup_model_constants(
    const input_parameters_t &input,
    const reference_quantities_t &ref_quant) {
    this->id = std::to_string(input.id);
    // Set constants
    this->constants[static_cast<int>(Cons_idx::q_crit_i)] = q_crit_i;
    this->constants[static_cast<int>(Cons_idx::D_crit_i)] = D_crit_i;
    this->constants[static_cast<int>(Cons_idx::D_conv_i)] = D_conv_i;
    this->constants[static_cast<int>(Cons_idx::q_crit_r)] = q_crit_r;
    this->constants[static_cast<int>(Cons_idx::D_crit_r)] = D_crit_r;
    this->constants[static_cast<int>(Cons_idx::q_crit_fr)] = q_crit_fr;
    this->constants[static_cast<int>(Cons_idx::D_coll_c)] = D_coll_c;
    this->constants[static_cast<int>(Cons_idx::q_crit)] = q_crit;
    this->constants[static_cast<int>(Cons_idx::D_conv_sg)] = D_conv_sg;
    this->constants[static_cast<int>(Cons_idx::D_conv_ig)] = D_conv_ig;
    this->constants[static_cast<int>(Cons_idx::x_conv)] = x_conv;
    this->constants[static_cast<int>(Cons_idx::parcel_height)] = parcel_height;
    this->constants[static_cast<int>(Cons_idx::alpha_spacefilling)] = alpha_spacefilling;
    this->constants[static_cast<int>(Cons_idx::T_nuc)] = T_nuc;
    this->constants[static_cast<int>(Cons_idx::T_freeze)] = T_freeze;
    this->constants[static_cast<int>(Cons_idx::T_f)] = T_f;
    this->constants[static_cast<int>(Cons_idx::D_eq)] = D_eq;
    this->constants[static_cast<int>(Cons_idx::rho_w)] = rho_w;
    this->constants[static_cast<int>(Cons_idx::rho_0)] = rho_0;
    this->constants[static_cast<int>(Cons_idx::rho_vel)] = rho_vel;
    this->constants[static_cast<int>(Cons_idx::rho_vel_c)] = rho_vel_c;
    this->constants[static_cast<int>(Cons_idx::rho_ice)] = rho_ice;
    this->constants[static_cast<int>(Cons_idx::M_w)] = M_w;
    this->constants[static_cast<int>(Cons_idx::M_a)] = M_a;
    this->constants[static_cast<int>(Cons_idx::R_universal)] = R_universal;
    this->constants[static_cast<int>(Cons_idx::Epsilon)] = Epsilon;
    this->constants[static_cast<int>(Cons_idx::gravity_acc)] = gravity_acc;
    this->constants[static_cast<int>(Cons_idx::R_a)] = R_a;
    this->constants[static_cast<int>(Cons_idx::R_v)] = R_v;
    this->constants[static_cast<int>(Cons_idx::a_v)] = a_v;
    this->constants[static_cast<int>(Cons_idx::b_v)] = b_v;
    this->constants[static_cast<int>(Cons_idx::a_prime)] = a_prime;
    this->constants[static_cast<int>(Cons_idx::b_prime)] = b_prime;
    this->constants[static_cast<int>(Cons_idx::c_prime)] = c_prime;
    this->constants[static_cast<int>(Cons_idx::K_T)] = K_T;
    this->constants[static_cast<int>(Cons_idx::L_wd)] = L_wd;
    this->constants[static_cast<int>(Cons_idx::L_ed)] = L_ed;
    this->constants[static_cast<int>(Cons_idx::D_v)] = D_v;
    this->constants[static_cast<int>(Cons_idx::ecoll_min)] = ecoll_min;
    this->constants[static_cast<int>(Cons_idx::ecoll_gg)] = ecoll_gg;
    this->constants[static_cast<int>(Cons_idx::ecoll_gg_wet)] = ecoll_gg_wet;
    this->constants[static_cast<int>(Cons_idx::kin_visc_air)] = kin_visc_air;
    this->constants[static_cast<int>(Cons_idx::C_mult)] = C_mult;
    this->constants[static_cast<int>(Cons_idx::T_mult_min)] = T_mult_min;
    this->constants[static_cast<int>(Cons_idx::T_mult_max)] = T_mult_max;
    this->constants[static_cast<int>(Cons_idx::T_mult_opt)] = T_mult_opt;
    this->constants[static_cast<int>(Cons_idx::kc_autocon)] = kc_autocon;
    this->constants[static_cast<int>(Cons_idx::D_rainfrz_gh)] = D_rainfrz_gh;
    this->constants[static_cast<int>(Cons_idx::D_rainfrz_ig)] = D_rainfrz_ig;
    this->constants[static_cast<int>(Cons_idx::dv0)] = dv0;
    this->constants[static_cast<int>(Cons_idx::p_sat_melt)] = p_sat_melt;
    this->constants[static_cast<int>(Cons_idx::cp)] = cp;
    this->constants[static_cast<int>(Cons_idx::k_b)] = k_b;
    this->constants[static_cast<int>(Cons_idx::a_HET)] = a_HET;
    this->constants[static_cast<int>(Cons_idx::b_HET)] = b_HET;
    this->constants[static_cast<int>(Cons_idx::N_sc)] = N_sc;
    this->constants[static_cast<int>(Cons_idx::n_f)] = n_f;
    this->constants[static_cast<int>(Cons_idx::N_avo)] = N_avo;
    this->constants[static_cast<int>(Cons_idx::na_dust)] = na_dust;
    this->constants[static_cast<int>(Cons_idx::na_soot)] = na_soot;
    this->constants[static_cast<int>(Cons_idx::na_orga)] = na_orga;
    this->constants[static_cast<int>(Cons_idx::ni_het_max)] = ni_het_max;
    this->constants[static_cast<int>(Cons_idx::ni_hom_max)] = ni_hom_max;
    this->constants[static_cast<int>(Cons_idx::a_dep)] = a_dep;
    this->constants[static_cast<int>(Cons_idx::b_dep)] = b_dep;
    this->constants[static_cast<int>(Cons_idx::c_dep)] = c_dep;
    this->constants[static_cast<int>(Cons_idx::d_dep)] = d_dep;
    this->constants[static_cast<int>(Cons_idx::nim_imm)] = nim_imm;
    this->constants[static_cast<int>(Cons_idx::nin_dep)] = nin_dep;
    this->constants[static_cast<int>(Cons_idx::alf_imm)] = alf_imm;
    this->constants[static_cast<int>(Cons_idx::bet_dep)] = bet_dep;
    this->constants[static_cast<int>(Cons_idx::bet_imm)] = bet_imm;
    this->constants[static_cast<int>(Cons_idx::r_const)] = r_const;
    this->constants[static_cast<int>(Cons_idx::r1_const)] = r1_const;
    this->constants[static_cast<int>(Cons_idx::cv)] = cv;
    this->constants[static_cast<int>(Cons_idx::p_sat_const_a)] = p_sat_const_a;
    this->constants[static_cast<int>(Cons_idx::p_sat_ice_const_a)] = p_sat_ice_const_a;
    this->constants[static_cast<int>(Cons_idx::p_sat_const_b)] = p_sat_const_b;
    this->constants[static_cast<int>(Cons_idx::p_sat_ice_const_b)] = p_sat_ice_const_b;
    this->constants[static_cast<int>(Cons_idx::p_sat_low_temp)] = p_sat_low_temp;
    this->constants[static_cast<int>(Cons_idx::T_sat_low_temp)] = T_sat_low_temp;
    this->constants[static_cast<int>(Cons_idx::alpha_depo)] = alpha_depo;
    this->constants[static_cast<int>(Cons_idx::r_0)] = r_0;

    this->constants[static_cast<int>(Cons_idx::k_1_conv)] = k_1_conv;
    this->constants[static_cast<int>(Cons_idx::k_2_conv)] = k_2_conv;
    this->constants[static_cast<int>(Cons_idx::k_1_accr)] = k_1_accr;
    this->constants[static_cast<int>(Cons_idx::k_r)] = k_r;

#if defined(CCN_AKM)
    this->constants[static_cast<int>(Cons_idx::p_ccn)] = p_ccn;
    this->constants[static_cast<int>(Cons_idx::h_ccn_1)] = h_ccn_1;
    this->constants[static_cast<int>(Cons_idx::h_ccn_2)] = h_ccn_2;
    this->constants[static_cast<int>(Cons_idx::g_ccn_1)] = g_ccn_1;
    this->constants[static_cast<int>(Cons_idx::g_ccn_2)] = g_ccn_2;
    this->constants[static_cast<int>(Cons_idx::g_ccn_3)] = g_ccn_3;
    this->constants[static_cast<int>(Cons_idx::i_ccn_1)] = i_ccn_1;
    this->constants[static_cast<int>(Cons_idx::i_ccn_2)] = i_ccn_2;
    this->constants[static_cast<int>(Cons_idx::hande_ccn_fac)] = hande_ccn_fac;
#endif
    // Numerics
    this->t_end_prime = input.t_end_prime;
    this->t_end = this->t_end_prime/ref_quant.tref;
    // Time of the substeps
    this->dt = input.dt_prime/ref_quant.tref;
    this->dt_prime = input.dt_prime;
    // Evaluate the general performance constants
    this->dt_half = this->dt*0.5;
    this->dt_third = this->dt/3.0;
    this->dt_sixth = this->dt/6.0;
    // Accomodation coefficient
    this->alpha_d = 1.0;

    // Performance constants for warm cloud; COSMO
    this->a1_scale = 1.0e-3;

    // Performance constants for warm cloud; IFS
    // The file constants.h also defines some constants as nar, ...
#ifndef MET3D
    const double Nc = 50;  // 50 over ocean; 300 over land
    const double F_aut = 1.5;
    const double F_acc = 2.0;
    // const double lambda_pp = pow(this->nar * this->ar * tgamma(this->br + 1.0) , this->alpha_r);
#endif
    for (uint32_t i=0; i < 4; i++) {
        this->constants[static_cast<int>(Cons_idx::a_ccn_1)+i] = a_ccn[i];
        this->constants[static_cast<int>(Cons_idx::b_ccn_1)+i] = b_ccn[i];
        this->constants[static_cast<int>(Cons_idx::c_ccn_1)+i] = c_ccn[i];
        this->constants[static_cast<int>(Cons_idx::d_ccn_1)+i] = d_ccn[i];
    }

    if (nuc_type == 6) {
        this->constants[static_cast<int>(Cons_idx::na_dust)] = na_dust;
        this->constants[static_cast<int>(Cons_idx::na_soot)] = na_soot;
        this->constants[static_cast<int>(Cons_idx::na_orga)] = na_orga;
    } else if (nuc_type == 7 || nuc_type == 5) {
        // Standard values
        this->constants[static_cast<int>(Cons_idx::na_dust)] = na_dust_2;
        this->constants[static_cast<int>(Cons_idx::na_soot)] = na_soot_2;
        this->constants[static_cast<int>(Cons_idx::na_orga)] = na_orga_2;
    } else if (nuc_type == 8) {
        this->constants[static_cast<int>(Cons_idx::na_dust)] = na_dust_3;
        this->constants[static_cast<int>(Cons_idx::na_soot)] = na_soot_3;
        this->constants[static_cast<int>(Cons_idx::na_orga)] = na_orga_3;
    }

    //// Exponents of the cloud model
    this->constants[static_cast<int>(Cons_idx::gamma)] = 1.0;
    this->constants[static_cast<int>(Cons_idx::betac)] = 1.0;
    this->constants[static_cast<int>(Cons_idx::betar)] = 7./8.;
    this->constants[static_cast<int>(Cons_idx::delta1)] = 0.5;
    this->constants[static_cast<int>(Cons_idx::delta2)] = 11./16.;
    this->constants[static_cast<int>(Cons_idx::zeta)] = 9./8.;
#if defined(RK4ICE) || defined(RK4NOICE)
    //// Cloud
    this->cloud.constants[static_cast<int>(Particle_cons_idx::nu)] = cloud_nu;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::mu)] = cloud_mu;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::max_x)] = cloud_max_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = cloud_min_x;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::a_geo)] = cloud_a_geo;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::b_geo)] = cloud_b_geo;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::a_vel)] = cloud_a_vel;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::b_vel)] = cloud_b_vel;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::a_ven)] = cloud_a_ven;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::b_ven)] = cloud_b_ven;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::cap)] = cloud_cap;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = cloud_vsedi_max;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = cloud_vsedi_min;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = cloud_q_crit_c;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = cloud_d_crit_c;

    //// Rain
    this->rain.constants[static_cast<int>(Particle_cons_idx::nu)] = rain_nu;
    this->rain.constants[static_cast<int>(Particle_cons_idx::mu)] = rain_mu;
    this->rain.constants[static_cast<int>(Particle_cons_idx::max_x)] = rain_max_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = rain_min_x;
    this->rain.constants[static_cast<int>(Particle_cons_idx::a_geo)] = rain_a_geo;
    this->rain.constants[static_cast<int>(Particle_cons_idx::b_geo)] = rain_b_geo;
    this->rain.constants[static_cast<int>(Particle_cons_idx::a_vel)] = rain_a_vel;
    this->rain.constants[static_cast<int>(Particle_cons_idx::b_vel)] = rain_b_vel;
    this->rain.constants[static_cast<int>(Particle_cons_idx::a_ven)] = rain_a_ven;
    this->rain.constants[static_cast<int>(Particle_cons_idx::b_ven)] = rain_b_ven;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cap)] = rain_cap;

    this->rain.constants[static_cast<int>(Particle_cons_idx::alpha)] = rain_alpha;
    this->rain.constants[static_cast<int>(Particle_cons_idx::beta)] = rain_beta;
    this->rain.constants[static_cast<int>(Particle_cons_idx::gamma)] = rain_gamma;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu0)] = rain_cmu0;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu1)] = rain_cmu1;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu2)] = rain_cmu2;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu3)] = rain_cmu3;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu4)] = rain_cmu4;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu5)] = rain_cmu5;
    this->constants[static_cast<int>(Cons_idx::rain_gfak)] = rain_rain_gfak;

    this->rain.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = rain_vsedi_max;
    this->rain.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = rain_vsedi_min;

    //// Graupel
    this->graupel.constants[static_cast<int>(Particle_cons_idx::nu)] = graupel_nu;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::mu)] = graupel_mu;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::max_x)] = graupel_max_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = graupel_min_x;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::a_geo)] = graupel_a_geo;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::b_geo)] = graupel_b_geo;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::a_vel)] = graupel_a_vel;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::b_vel)] = graupel_b_vel;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::a_ven)] = graupel_a_ven;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::b_ven)] = graupel_b_ven;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::cap)] = graupel_cap;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = graupel_vsedi_max;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = graupel_vsedi_min;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = graupel_d_crit_c;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = graupel_q_crit_c;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::s_vel)] = graupel_s_vel;

    //// Hail
    this->hail.constants[static_cast<int>(Particle_cons_idx::nu)] = hail_nu;
    this->hail.constants[static_cast<int>(Particle_cons_idx::mu)] = hail_mu;
    this->hail.constants[static_cast<int>(Particle_cons_idx::max_x)] = hail_max_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = hail_min_x;
    this->hail.constants[static_cast<int>(Particle_cons_idx::a_geo)] = hail_a_geo;
    this->hail.constants[static_cast<int>(Particle_cons_idx::b_geo)] = hail_b_geo;
    this->hail.constants[static_cast<int>(Particle_cons_idx::a_vel)] = hail_a_vel;
    this->hail.constants[static_cast<int>(Particle_cons_idx::b_vel)] = hail_b_vel;
    this->hail.constants[static_cast<int>(Particle_cons_idx::a_ven)] = hail_a_ven;
    this->hail.constants[static_cast<int>(Particle_cons_idx::b_ven)] = hail_b_ven;
    this->hail.constants[static_cast<int>(Particle_cons_idx::cap)] = hail_cap;
    this->hail.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = hail_vsedi_max;
    this->hail.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = hail_vsedi_min;
    this->hail.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = hail_d_crit_c;
    this->hail.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = hail_q_crit_c;
    this->hail.constants[static_cast<int>(Particle_cons_idx::s_vel)] = hail_s_vel;

    //// Ice
    this->ice.constants[static_cast<int>(Particle_cons_idx::nu)] = ice_nu;
    this->ice.constants[static_cast<int>(Particle_cons_idx::mu)] = ice_mu;
    this->ice.constants[static_cast<int>(Particle_cons_idx::max_x)] = ice_max_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = ice_min_x;
    this->ice.constants[static_cast<int>(Particle_cons_idx::a_geo)] = ice_a_geo;
    this->ice.constants[static_cast<int>(Particle_cons_idx::b_geo)] = ice_b_geo;
    this->ice.constants[static_cast<int>(Particle_cons_idx::a_vel)] = ice_a_vel;
    this->ice.constants[static_cast<int>(Particle_cons_idx::b_vel)] = ice_b_vel;
    this->ice.constants[static_cast<int>(Particle_cons_idx::a_ven)] = ice_a_ven;
    this->ice.constants[static_cast<int>(Particle_cons_idx::b_ven)] = ice_b_ven;
    this->ice.constants[static_cast<int>(Particle_cons_idx::cap)] = ice_cap;
    this->ice.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = ice_vsedi_max;
    this->ice.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = ice_vsedi_min;
    this->ice.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = ice_d_crit_c;
    this->ice.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = ice_q_crit_c;
    this->ice.constants[static_cast<int>(Particle_cons_idx::s_vel)] = ice_s_vel;

    //// Snow
    this->snow.constants[static_cast<int>(Particle_cons_idx::nu)] = snow_nu;
    this->snow.constants[static_cast<int>(Particle_cons_idx::mu)] = snow_mu;
    this->snow.constants[static_cast<int>(Particle_cons_idx::max_x)] = snow_max_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = snow_min_x;
    this->snow.constants[static_cast<int>(Particle_cons_idx::a_geo)] = snow_a_geo;
    this->snow.constants[static_cast<int>(Particle_cons_idx::b_geo)] = snow_b_geo;
    this->snow.constants[static_cast<int>(Particle_cons_idx::a_vel)] = snow_a_vel;
    this->snow.constants[static_cast<int>(Particle_cons_idx::b_vel)] = snow_b_vel;
    this->snow.constants[static_cast<int>(Particle_cons_idx::a_ven)] = snow_a_ven;
    this->snow.constants[static_cast<int>(Particle_cons_idx::b_ven)] = snow_b_ven;
    this->snow.constants[static_cast<int>(Particle_cons_idx::cap)] = snow_cap;
    this->snow.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = snow_vsedi_max;
    this->snow.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = snow_vsedi_min;
    this->snow.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = snow_d_crit_c;
    this->snow.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = snow_q_crit_c;
    this->snow.constants[static_cast<int>(Particle_cons_idx::s_vel)] = snow_s_vel;
#endif
    this->constants[static_cast<int>(Cons_idx::inv_z)] = 1.0/parcel_height;
    setup_dependent_model_constants();
    // Set the uncertainty of every parameter to scale the gradients.
    if (input.simulation_mode == create_train_set) {
        set_uncertainty(1.0);
    } else {
        // Default value of 10%.
        set_uncertainty();
    }
}


template<class float_t>
void model_constants_t<float_t>::setup_dependent_model_constants() {
    // Performance constants for warm cloud; COSMO
    this->a2_scale = 1.72 / pow(get_at(this->constants, Cons_idx::R_a) , 7./8.);
    this->e1_scale = 1.0 / sqrt(get_at(this->constants, Cons_idx::R_a));
    this->e2_scale = 9.1 / pow(get_at(this->constants, Cons_idx::R_a) , 11./16.);
    this->d_scale = (130.0*tgamma(4.5))
        / (6.0*(1.0e3)*pow(M_PI*(8.0e6)*get_at(this->constants, Cons_idx::R_a) , 1.0/8.0));
#ifndef MET3D
    // Performance constants for warm cloud; IFS
    const double lambda_pp = pow(this->nar * this->ar * tgamma(this->br + 1.0) , this->alpha_r);
#endif
#if defined(RK4ICE) || defined(RK4NOICE)
    this->cloud.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0
        / get_at(this->cloud.constants, Particle_cons_idx::cap);
    this->cloud.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->cloud, 1);
    this->cloud.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->cloud, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->cloud.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->cloud, 2);

    this->setup_cloud_autoconversion(this->cloud);
    setup_bulk_sedi(this->cloud);

    this->rain.constants[static_cast<int>(Particle_cons_idx::nm1)] =
        (get_at(this->rain.constants, Particle_cons_idx::nu)+1.0)/get_at(this->rain.constants, Particle_cons_idx::mu);
    this->rain.constants[static_cast<int>(Particle_cons_idx::nm2)] =
        (get_at(this->rain.constants, Particle_cons_idx::nu)+2.0)/get_at(this->rain.constants, Particle_cons_idx::mu);
    this->rain.constants[static_cast<int>(Particle_cons_idx::nm3)] =
        (get_at(this->rain.constants, Particle_cons_idx::nu)+3.0)/get_at(this->rain.constants, Particle_cons_idx::mu);

    this->table_r1.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(this->rain.constants, Particle_cons_idx::nm1).getValue());
    this->table_r2.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(this->rain.constants, Particle_cons_idx::nm2).getValue());
    this->table_r3.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(this->rain.constants, Particle_cons_idx::nm3).getValue());
    this->rain.constants[static_cast<int>(Particle_cons_idx::g1)] = this->table_r1.igf[this->table_r1.n_bins-1];
    this->rain.constants[static_cast<int>(Particle_cons_idx::g2)] = this->table_r2.igf[this->table_r2.n_bins-1];
    this->rain.constants[static_cast<int>(Particle_cons_idx::c_s)] =
        1.0 / get_at(this->rain.constants, Particle_cons_idx::cap);
    this->rain.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->rain, 1);
    this->rain.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->rain, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->rain.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->rain, 2);
    setup_bulk_sedi(this->rain);

    this->graupel.constants[static_cast<int>(Particle_cons_idx::nm1)] =
        (get_at(this->graupel.constants, Particle_cons_idx::nu)+1.0)
        / get_at(this->graupel.constants, Particle_cons_idx::mu);
    this->graupel.constants[static_cast<int>(Particle_cons_idx::nm2)] =
        (get_at(this->graupel.constants, Particle_cons_idx::nu)+2.0)
        / get_at(this->graupel.constants, Particle_cons_idx::mu);
    float_t a =
        (get_at(this->graupel.constants, Particle_cons_idx::nu)+1.0)
        / get_at(this->graupel.constants, Particle_cons_idx::mu);
    this->table_g1.init_gamma_table(
        n_lookup, n_lookup_hr_dummy, get_at(this->graupel.constants, Particle_cons_idx::nm1).getValue());
    a = (get_at(this->graupel.constants, Particle_cons_idx::nu)+2.0)
        / get_at(this->graupel.constants, Particle_cons_idx::mu);
    this->table_g2.init_gamma_table(
        n_lookup, n_lookup_hr_dummy, get_at(this->graupel.constants, Particle_cons_idx::nm2).getValue());
    this->graupel.constants[static_cast<int>(Particle_cons_idx::g1)] = this->table_g1.igf[this->table_g1.n_bins-1];
    this->graupel.constants[static_cast<int>(Particle_cons_idx::g2)] = this->table_g2.igf[this->table_g2.n_bins-1];
    this->graupel.constants[static_cast<int>(Particle_cons_idx::c_s)] =
        1.0 / get_at(this->graupel.constants, Particle_cons_idx::cap);
    this->graupel.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = graupel_ecoll_c;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->graupel, 1);
    this->graupel.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->graupel, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->graupel.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->graupel, 2);
    setup_bulk_sedi(this->graupel);

    this->hail.constants[static_cast<int>(Particle_cons_idx::c_s)] =
        1.0 / get_at(this->hail.constants, Particle_cons_idx::cap);
    this->hail.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = hail_ecoll_c;
    this->hail.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->hail, 1);
    this->hail.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->hail, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->hail.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->hail, 2);
    setup_bulk_sedi(this->hail);

    this->ice.constants[static_cast<int>(Particle_cons_idx::c_s)] =
        1.0 / get_at(this->ice.constants, Particle_cons_idx::cap);
    this->ice.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = ice_ecoll_c;
    this->ice.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->ice, 1);
    this->ice.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->ice, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->ice.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->ice, 2);
    setup_bulk_sedi(this->ice);

    this->snow.constants[static_cast<int>(Particle_cons_idx::c_s)] =
        1.0 / get_at(this->snow.constants, Particle_cons_idx::cap);
    this->snow.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = snow_ecoll_c;
    this->snow.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->snow, 1);
    this->snow.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->snow, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->snow.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->snow, 2);
    setup_bulk_sedi(this->snow);

    this->constants[static_cast<int>(Cons_idx::const0)] =
        1.0/(get_at(this->constants, Cons_idx::D_coll_c)
        - get_at(this->cloud.constants, Particle_cons_idx::d_crit_c));
    this->constants[static_cast<int>(Cons_idx::const3)] =
        1.0/(get_at(this->constants, Cons_idx::T_mult_opt) - get_at(this->constants, Cons_idx::T_mult_min));
    this->constants[static_cast<int>(Cons_idx::const4)] =
        1.0/(get_at(this->constants, Cons_idx::T_mult_opt) - get_at(this->constants, Cons_idx::T_mult_max));
    this->constants[static_cast<int>(Cons_idx::const5)] =
        get_at(this->constants, Cons_idx::alpha_spacefilling)
        * get_at(this->constants, Cons_idx::rho_w)/get_at(this->constants, Cons_idx::rho_ice);

    init_particle_collection_1(this->snow, this->cloud, this->coeffs_scr);
    init_particle_collection_2(this->snow, this->rain, this->coeffs_srr);
    init_particle_collection_2(this->ice, this->rain, this->coeffs_irr);
    init_particle_collection_1(this->ice, this->cloud, this->coeffs_icr);
    init_particle_collection_1(this->hail, this->rain, this->coeffs_hrr);
    init_particle_collection_1(this->graupel, this->rain, this->coeffs_grr);  // Cosmo uses 2, ICON uses 1
    init_particle_collection_1(this->hail, this->cloud, this->coeffs_hcr);
    init_particle_collection_1(this->graupel, this->cloud, this->coeffs_gcr);
    init_particle_collection_1(this->snow, this->ice, this->coeffs_sic);
    init_particle_collection_1(this->hail, this->ice, this->coeffs_hic);
    init_particle_collection_1(this->graupel, this->ice, this->coeffs_gic);
    init_particle_collection_1(this->hail, this->snow, this->coeffs_hsc);
    init_particle_collection_1(this->graupel, this->snow, this->coeffs_gsc);

    // Setup graupel, snow and ice selfcollection
    this->graupel.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = M_PI/8.0
        * (2.0*coll_delta_11(this->graupel,  0)
           + coll_delta_12(this->graupel, this->graupel, 0))
        * sqrt((2.0*coll_theta_11(this->graupel,  0)
           - coll_theta_12(this->graupel, this->graupel, 0)));

    this->snow.constants[static_cast<int>(Particle_cons_idx::sc_delta_n)] = (
        2.0*coll_delta_11(this->snow, 0) + coll_delta_12(this->snow, this->snow, 0));
    this->snow.constants[static_cast<int>(Particle_cons_idx::sc_theta_n)] = (
        2.0*coll_theta_11(this->snow,  0) - coll_theta_12(this->snow, this->snow, 0));

    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_delta_n)] = coll_delta_11(this->ice, 0)
        + coll_delta_12(this->ice, this->ice, 0)
        + coll_delta_22(this->ice, 0);
    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_delta_q)] = coll_delta_11(this->ice, 0)
        + coll_delta_12(this->ice, this->ice, 1)
        + coll_delta_22(this->ice, 1);
    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_theta_n)] = coll_theta_11(this->ice,  0)
        - coll_theta_12(this->ice, this->ice, 0)
        + coll_theta_22(this->ice, 0);
    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_theta_q)] = coll_theta_11(this->ice, 0)
        - coll_theta_12(this->ice, this->ice, 1)
        + coll_theta_22(this->ice, 1);
#endif
}


template<class float_t>
void model_constants_t<float_t>::set_uncertainty() {
    set_uncertainty(0.1);
}


template<class float_t>
void model_constants_t<float_t>::set_uncertainty(double scale) {
    for (uint32_t i=0; i < static_cast<uint32_t>(Cons_idx::n_items); ++i) {
        this->uncertainty[i] = this->constants[i].getValue() * scale;
    }
    for (uint32_t i=0; i < static_cast<uint32_t>(Particle_cons_idx::n_items); ++i) {
        this->rain.uncertainty[i] = this->rain.constants[i].getValue() * scale;
    }
    for (uint32_t i=0; i < static_cast<uint32_t>(Particle_cons_idx::n_items); ++i) {
        this->cloud.uncertainty[i] = this->cloud.constants[i].getValue() * scale;
    }
    for (uint32_t i=0; i < static_cast<uint32_t>(Particle_cons_idx::n_items); ++i) {
        this->graupel.uncertainty[i] = this->graupel.constants[i].getValue() * scale;
    }
    for (uint32_t i=0; i < static_cast<uint32_t>(Particle_cons_idx::n_items); ++i) {
        this->hail.uncertainty[i] = this->hail.constants[i].getValue() * scale;
    }
    for (uint32_t i=0; i < static_cast<uint32_t>(Particle_cons_idx::n_items); ++i) {
        this->ice.uncertainty[i] = this->ice.constants[i].getValue() * scale;
    }
    for (uint32_t i=0; i < static_cast<uint32_t>(Particle_cons_idx::n_items); ++i) {
        this->snow.uncertainty[i] = this->snow.constants[i].getValue() * scale;
    }
    for (uint32_t i=0; i < static_cast<uint32_t>(Init_cons_idx::n_items); ++i) {
        this->uncertainty[i + static_cast<uint32_t>(Cons_idx::n_items)] =
                this->initial_conditions[i].getValue() * scale;
    }
}


template<class float_t>
void model_constants_t<float_t>::set_dt(
    const double dt_prime,
    const reference_quantities_t &ref_quant) {
    this->dt_traject_prime = dt_prime;
    this->dt_traject = this->dt_traject_prime/ref_quant.tref;
    this->num_steps = floor(this->t_end_prime/this->dt_traject_prime);
    // The trajectories from input files are calculated with dt_traject_prime s timesteps.
    this->num_sub_steps = (floor(this->dt_traject_prime/this->dt_prime) < 1) ?
            1 : floor(this->dt_traject_prime/this->dt_prime);
    }

template<class float_t>
void model_constants_t<float_t>::print(std::ostream &os) {
#ifdef SILENT_MODE
    return;
#endif
  os << "\nModel constants:\n"
        << "----------------\n"
        << "Code for tracking model states = " << this->track_state << "\n"
        << "Codes for tracking model parameters = ";
    for (auto const &t : track_param) {
        os << t << ", ";
    }
    os << "\nCode for tracking initial conditions = " << this->track_ic << "\n"
        << "Final integration time = " << this->t_end_prime << " seconds\n"
        << "Nondimensional final integration time = " << this->t_end << "\n"
        << "Timestep = " << this->dt_prime << " seconds\n"
        << "Nondimensional timestep = " << this->dt << "\n"
        << "Number of iterations = " << this->num_steps << "\n"
        << "Number of substeps = " << this->num_sub_steps << "\n";
    for (auto const &t : table_param) {
        os << t.first << " = " << get_at(this->constants, t.second) << "\n";
    }
    os << std::endl << std::flush;
}


template<class float_t>
bool model_constants_t<float_t>::trace_check(
    const int &idx,
    const int state_param) const {
    if (state_param == 1) {
        return (track_state & (((uint64_t) 1) << idx));
    } else if (state_param == 0) {
        uint32_t i = idx/64;
        return (track_param[i] & (((uint64_t) 1) << (idx%64)));
    } else {
        return (track_ic & (((uint64_t) 1) << idx));
    }
}


template<class float_t>
void model_constants_t<float_t>::load_configuration(
    const std::string &filename) {
    nlohmann::json config;
    std::ifstream ifstr(filename);
    ifstr >> config;

    if (config.find("model_state_variable") != config.end()) {
        local_num_comp = 0;
        track_state = 0;
        for (const auto &c : config["model_state_variable"].items()) {
            const uint32_t id = c.value();
            track_state = track_state | (((uint64_t) 1) << id);
            local_num_comp++;
        }
    } else {
        track_state = -1;
        local_num_comp = num_comp;
    }

    if (config.find("out_params") != config.end()) {
        local_num_par = 0;
        for (auto &t : track_param) t = 0;
        for (const auto &c : config["out_params"].items()) {
            const std::string id_name = c.value();
            auto it_tmp = std::find(
                output_grad_idx.begin(),
                output_grad_idx.end(),
                id_name);
            int id = std::distance(output_grad_idx.begin(), it_tmp);
            uint32_t idx = id/64;
            track_param[idx] = track_param[idx] | (((uint64_t) 1) << (id%64));
            local_num_par++;
        }
    } else {
        for (auto &t : track_param)
            t = -1;
        local_num_par = num_par;
    }

    if (config.find("initial_condition") != config.end()) {
        track_ic = 0;
        local_ic_par = 0;
        for (const auto &c : config["initial_condition"].items()) {
            const std::string id_name = c.value();
            auto it_tmp = std::find(
                init_grad_idx.begin(),
                init_grad_idx.end(),
                id_name);
            int id = std::distance(init_grad_idx.begin(), it_tmp);
            track_ic = track_ic | (((uint64_t) 1) << id);
            local_ic_par++;
        }
    } else {
        track_ic = -1;
        local_ic_par = static_cast<uint32_t>(Init_cons_idx::n_items);
    }
}

#ifdef OUT_DOUBLE
template<class float_t>
void model_constants_t<float_t>::get_perturbed_info(
    std::vector<double> &perturbed,
    std::vector<uint64_t> &param_idx) const {
#else
template<class float_t>
void model_constants_t<float_t>::get_perturbed_info(
    std::vector<float> &perturbed,
    std::vector<uint64_t> &param_idx) const {
#endif

    for (const auto idx : perturbed_idx) {
#ifdef DEBUG_SEG
        std::cout << "rank " << rank << " idx " << idx << " size " << constants.size() <<
        " pert_size " << perturbed_idx.size() <<
        " empty? " << perturbed_idx.empty() << "\n";
#endif
        perturbed.push_back(constants[idx].getValue());
        param_idx.push_back(idx);
    }
    uint64_t offset = constants.size();
    for (const auto idx : cloud.perturbed_idx) {
        perturbed.push_back(cloud.constants[idx].getValue());
        param_idx.push_back(idx + offset);
    }
    offset += cloud.constants.size();
    for (const auto idx : rain.perturbed_idx) {
        perturbed.push_back(rain.constants[idx].getValue());
        param_idx.push_back(idx + offset);
    }
    offset += rain.constants.size();
    for (const auto idx : ice.perturbed_idx) {
        perturbed.push_back(ice.constants[idx].getValue());
        param_idx.push_back(idx + offset);
    }
    offset += ice.constants.size();
    for (const auto idx : snow.perturbed_idx) {
        perturbed.push_back(snow.constants[idx].getValue());
        param_idx.push_back(idx + offset);
    }
    offset += snow.constants.size();
    for (const auto idx : graupel.perturbed_idx) {
        perturbed.push_back(graupel.constants[idx].getValue());
        param_idx.push_back(idx + offset);
    }
    offset += graupel.constants.size();
    for (const auto idx : hail.perturbed_idx) {
        perturbed.push_back(hail.constants[idx].getValue());
        param_idx.push_back(idx + offset);
    }
}

template<class float_t>
void model_constants_t<float_t>::set_dw_polynomial(
    const double &w0,
    const double &w1,
    const double &w2,
    const int &offset_1,
    const int &offset_2
    ) {
    poly0 = w0;
    if (offset_1 == 0 && offset_2 == 0) {
        poly1 = 0;
        poly2 = 0;
    } else if (offset_2 == 0) {
        poly1 = (w1-w0) / offset_1;
        poly2 = 0;
    } else {
        // x0 is 0
        double x1 = offset_1;
        double x2 = offset_2 + x1;
        poly2 = (w1 - w2*x1/x2 - w0 * (1-x1/x2)) / (x1*x1 - x2*x1);
        poly1 = w1/x1 - poly2*x1 - w0/x1;
    }
    current_w_poly_idx = 0;
}

template<class float_t>
double model_constants_t<float_t>::get_dw(
) const {
    return 2*poly2*current_w_poly_idx + poly1;
}

template<class float_t>
void model_constants_t<float_t>::increment_w_idx() {
    current_w_poly_idx += dt_prime;
}

template<class float_t>
double model_constants_t<float_t>::get_w(
) const {
    return poly2*current_w_poly_idx*current_w_poly_idx + poly1*current_w_poly_idx + poly0;
}

template class model_constants_t<codi::RealReverse>;
template class model_constants_t<codi::RealForwardVec<num_par_init> >;

template void to_json<codi::RealReverse>(
    nlohmann::json&, const model_constants_t<codi::RealReverse>&);
template void to_json<codi::RealForwardVec<num_par_init> >(
    nlohmann::json&, const model_constants_t<codi::RealForwardVec<num_par_init> >&);
