#include "include/types/model_constants_t.h"

model_constants_t::model_constants_t()
{
    id = "0";
    ensemble_id = 0;
    traj_id = 0;
    n_trajs = 1;
    n_ensembles = 1;
    ens_desc = "root ";
    constants.resize(static_cast<int>(Cons_idx::n_items));
    std::fill(constants.begin(), constants.end(), 0);
}

void model_constants_t::register_input(
    codi::RealReverse::TapeType &tape)
{
    for(auto &c: this->constants)
        tape.registerInput(c);
    this->rain.register_input(tape);
    this->cloud.register_input(tape);
    this->hail.register_input(tape);
    this->ice.register_input(tape);
    this->snow.register_input(tape);
    this->graupel.register_input(tape);
}


void model_constants_t::get_gradients(
    std::vector<codi::RealReverse> &y_single_new,
    std::vector< std::array<double, num_par > > &y_diff,
    codi::RealReverse::TapeType &tape)
{
    for(uint32_t ii = 0 ; ii < num_comp ; ii++)
        tape.registerOutput(y_single_new[ii]);

    tape.setPassive();
    for(uint32_t ii = 0 ; ii < num_comp ; ii++)
    {
        y_single_new[ii].setGradient(1.0);
        tape.evaluate();

        this->get_gradient(y_diff[ii]);
        tape.clearAdjoints();
    }
    tape.reset();
}

void model_constants_t::get_gradient(
    std::array<double, num_par> &out_vec) const
{
    for(int i=0; i<static_cast<int>(Cons_idx::n_items); ++i)
        out_vec[i] = this->constants[i].getGradient();
#if defined(RK4ICE) || defined(RK4NOICE)

    uint32_t idx = static_cast<uint32_t>(Cons_idx::n_items);
    this->rain.get_gradient(out_vec, idx);
    this->cloud.get_gradient(out_vec, idx);
    this->graupel.get_gradient(out_vec, idx);
    this->hail.get_gradient(out_vec, idx);
    this->ice.get_gradient(out_vec, idx);
    this->snow.get_gradient(out_vec, idx);
#endif
}
// template void model_constants_t::get_gradient<std::array<double, num_par> >(std::array<double, num_par>);

void model_constants_t::put(
    pt::ptree &ptree) const
{
    pt::ptree model_cons;
    model_cons.put("id", id);
    model_cons.put("ensemble_id", ensemble_id);
    model_cons.put("traj_id", traj_id);
    model_cons.put("n_trajs", n_trajs);
    model_cons.put("ens_desc", ens_desc);

    // technical parameters
    model_cons.put("t_end_prime", t_end_prime);
    model_cons.put("t_end", t_end);
    model_cons.put("dt_prime", dt_prime);
    model_cons.put("dt", dt);
    model_cons.put("dt_traject_prime", dt_traject_prime);
    model_cons.put("dt_traject", dt_traject);
    model_cons.put("num_steps", num_steps);
    model_cons.put("num_sub_steps", num_sub_steps);
    model_cons.put("dt_half", dt_half);
    model_cons.put("dt_sixth", dt_sixth);
    model_cons.put("dt_third", dt_third);
    model_cons.put("ensemble_id", ensemble_id);

    if(!perturbed_idx.empty())
    {
        pt::ptree perturbed;
        for(uint32_t idx: perturbed_idx)
            perturbed.put(std::to_string(idx), constants[idx]);
        model_cons.add_child("perturbed", perturbed);
    }

    // particle_model_constants
    hail.put(model_cons, "hail");
    ice.put(model_cons, "ice");
    snow.put(model_cons, "snow");
    cloud.put(model_cons, "cloud");
    rain.put(model_cons, "rain");
    graupel.put(model_cons, "graupel");

    // collection coefficients depend completely on particle
    // model constants and can be derived from there. So
    // we skip adding that to the checkpoint.

    ptree.add_child("model_constants", model_cons);
}

int model_constants_t::from_pt(
    pt::ptree &ptree)
{
    int err = 0;
    for(auto &it: ptree.get_child("model_constants"))
    {
        auto first = it.first;
        if(first == "id")
        {
            id = it.second.get_value<std::string>() + "-" + id;
        } else if(first == "ensemble_id")
        {
            ensemble_id = it.second.get_value<std::uint64_t>();
        } else if(first == "traj_id")
        {
            traj_id = it.second.get_value<std::uint64_t>();
        } else if(first == "n_trajs")
        {
            n_trajs = it.second.get_value<std::uint64_t>();
        } else if(first == "ens_desc")
        {
            ens_desc = it.second.get_value<std::string>();
        } else if(first == "t_end_prime")
        {
            t_end_prime = it.second.get_value<double>();
        } else if(first == "t_end")
        {
            t_end = it.second.get_value<double>();
        } else if(first == "dt_prime")
        {
            dt_prime = it.second.get_value<double>();
        } else if(first == "dt")
        {
            dt = it.second.get_value<double>();
        } else if(first == "dt_traject_prime")
        {
            dt_traject_prime = it.second.get_value<double>();
        } else if(first == "dt_traject")
        {
            dt_traject = it.second.get_value<double>();
        } else if(first == "num_steps")
        {
            num_steps = it.second.get_value<uint64_t>();
        } else if(first == "num_sub_steps")
        {
            num_sub_steps = it.second.get_value<uint64_t>();
        } else if(first == "dt_half")
        {
            dt_half = it.second.get_value<double>();
        } else if(first == "dt_sixth")
        {
            dt_sixth = it.second.get_value<double>();
        } else if(first == "dt_third")
        {
            dt_third = it.second.get_value<double>();
        } else if(first == "ensemble_id")
        {
            ensemble_id = it.second.get_value<uint64_t>();
        } else if(first == "perturbed")
        {
            for(auto &it2: ptree.get_child("model_constants.perturbed"))
            {
                uint32_t idx = std::stoi(it2.first);
                perturbed_idx.push_back(idx);
                this->constants[idx] = it2.second.get_value<double>();
            }
        // below from here: perturbed particle models
        } else if(first == "hail")
        {
            for(auto &it2: ptree.get_child("model_constants.hail.perturbed"))
            {
                uint32_t idx = std::stoi(it2.first);
                this->hail.perturbed_idx.push_back(idx);
                this->hail.constants[idx] = it2.second.get_value<double>();
            }
        } else if(first == "ice")
        {
            for(auto &it2: ptree.get_child("model_constants.ice.perturbed"))
            {
                uint32_t idx = std::stoi(it2.first);
                this->ice.perturbed_idx.push_back(idx);
                this->ice.constants[idx] = it2.second.get_value<double>();
            }
        } else if(first == "snow")
        {
            for(auto &it2: ptree.get_child("model_constants.snow.perturbed"))
            {
                uint32_t idx = std::stoi(it2.first);
                this->snow.perturbed_idx.push_back(idx);
                this->snow.constants[idx] = it2.second.get_value<double>();
            }
        } else if(first == "cloud")
        {
            for(auto &it2: ptree.get_child("model_constants.cloud.perturbed"))
            {
                uint32_t idx = std::stoi(it2.first);
                this->cloud.perturbed_idx.push_back(idx);
                this->cloud.constants[idx] = it2.second.get_value<double>();
            }
        } else if(first == "rain")
        {
            for(auto &it2: ptree.get_child("model_constants.rain.perturbed"))
            {
                uint32_t idx = std::stoi(it2.first);
                this->rain.perturbed_idx.push_back(idx);
                this->rain.constants[idx] = it2.second.get_value<double>();
            }
        } else if(first == "graupel")
        {
            for(auto &it2: ptree.get_child("model_constants.graupel.perturbed"))
            {
                uint32_t idx = std::stoi(it2.first);
                this->graupel.perturbed_idx.push_back(idx);
                this->graupel.constants[idx] = it2.second.get_value<double>();
            }
        } else
        {
            err = MODEL_CONS_CHECKPOINT_ERR;
        }
    }
    return err;
}

void model_constants_t::setCoefficients(
    std::vector<codi::RealReverse> & y,
    reference_quantities_t& ref)
{
    codi::RealReverse p_prime = y[p_idx]*ref.pref;
    codi::RealReverse T_prime = y[T_idx]*ref.Tref;

    codi::RealReverse rho_prime = p_prime /( get_at(this->constants, Cons_idx::R_a) * T_prime );
    codi::RealReverse L_vap_prime = latent_heat_water(T_prime, get_at(this->constants, Cons_idx::M_w));
    codi::RealReverse Ka_prime = thermal_conductivity_dry_air(T_prime);
    codi::RealReverse psat_prime = saturation_pressure_water(
        T_prime,
        get_at(this->constants, Cons_idx::p_sat_low_temp),
        get_at(this->constants, Cons_idx::p_sat_const_a),
        get_at(this->constants, Cons_idx::T_sat_low_temp),
        get_at(this->constants, Cons_idx::p_sat_const_b));
    codi::RealReverse A_pp = (L_vap_prime/(Ka_prime*T_prime))
        * ((L_vap_prime/(get_at(this->constants, Cons_idx::R_v)*T_prime)) - 1.0);
    codi::RealReverse B_pp = (get_at(this->constants, Cons_idx::R_v)*T_prime)/((2.21/p_prime)*psat_prime);


    this->constants[static_cast<int>(Cons_idx::e1_prime)] = this->e1_scale
        * ( pow(rho_prime, 2.0*this->alpha_r-2.0)/(A_pp + B_pp) );
    this->constants[static_cast<int>(Cons_idx::e2_prime)] = this->e2_scale
        * ( pow(rho_prime, this->alpha_r*this->epsilonr - (7.0/4.0))/(A_pp + B_pp) );

    this->constants[static_cast<int>(Cons_idx::a1_prime)] = this->a1_scale;	// Constant coefficient
    this->constants[static_cast<int>(Cons_idx::a2_prime)] = this->a2_scale;	// Constant coefficient
    this->constants[static_cast<int>(Cons_idx::d_prime)] = this->d_scale;	// Constant coefficient
}

void model_constants_t::setCoefficients(
    codi::RealReverse p_prime,
    codi::RealReverse T_prime)
{
    codi::RealReverse rho_prime = p_prime /( get_at(this->constants, Cons_idx::R_a) * T_prime );
    codi::RealReverse L_vap_prime = latent_heat_water(T_prime, get_at(this->constants, Cons_idx::M_w));
    codi::RealReverse Ka_prime = thermal_conductivity_dry_air(T_prime);
    codi::RealReverse psat_prime = saturation_pressure_water(
        T_prime,
        get_at(this->constants, Cons_idx::p_sat_low_temp),
        get_at(this->constants, Cons_idx::p_sat_const_a),
        get_at(this->constants, Cons_idx::T_sat_low_temp),
        get_at(this->constants, Cons_idx::p_sat_const_b));
    codi::RealReverse A_pp = (L_vap_prime/(Ka_prime*T_prime))*((L_vap_prime/(get_at(this->constants, Cons_idx::R_v)*T_prime)) - 1.0);
    codi::RealReverse B_pp = (get_at(this->constants, Cons_idx::R_v)*T_prime)/((2.21/p_prime)*psat_prime);

    this->constants[static_cast<int>(Cons_idx::a1_prime)] = this->a1_scale;	// Constant coefficient
    this->constants[static_cast<int>(Cons_idx::a2_prime)] = this->a2_scale;	// Constant coefficient

    this->constants[static_cast<int>(Cons_idx::e1_prime)] = this->e1_scale * ( pow(rho_prime, 2.0*this->alpha_r-2.0)/(A_pp + B_pp) );
    this->constants[static_cast<int>(Cons_idx::e2_prime)] = this->e2_scale * ( pow(rho_prime, this->alpha_r*this->epsilonr - (7.0/4.0))/(A_pp + B_pp) );

    this->constants[static_cast<int>(Cons_idx::d_prime)] = this->d_scale;	// Constant coefficient
}

void model_constants_t::setup_cloud_autoconversion(
    particle_model_constants_t &pc)
{
    auto nu = get_at(pc.constants, Particle_cons_idx::nu) + 1.0;
    auto mu = get_at(pc.constants, Particle_cons_idx::mu);
    if(get_at(pc.constants, Particle_cons_idx::mu) == 1.0)
    {
        this->constants[static_cast<int>(Cons_idx::cloud_k_au)] =
            get_at(this->constants, Cons_idx::kc_autocon)
            / get_at(pc.constants, Particle_cons_idx::max_x) * 0.05
            * (nu+1.0)*(nu+3.0) / pow(nu, 2);
        this->constants[static_cast<int>(Cons_idx::cloud_k_sc)] =
            get_at(this->constants, Cons_idx::kc_autocon) * (nu+1.0)/(nu);
    } else
    {
        this->constants[static_cast<int>(Cons_idx::cloud_k_au)] =
            get_at(this->constants, Cons_idx::kc_autocon)
            / get_at(pc.constants, Particle_cons_idx::max_x) * 0.05
            * (2.0 * tgamma((nu+3.0)/mu)
            * tgamma((nu+1.0)/mu) * pow(tgamma((nu)/mu), 2)
            - 1.0 * pow(tgamma((nu+2.0)/mu), 2) * pow(tgamma((nu)/mu), 2))
            / pow(tgamma((nu+1.0)/mu), 4);
        this->constants[static_cast<int>(Cons_idx::cloud_k_sc)] =
            get_at(this->constants, Cons_idx::kc_autocon)
            * get_at(pc.constants, Particle_cons_idx::c_z);
    }
}

void model_constants_t::setup_model_constants(
    const input_parameters_t &input,
    const reference_quantities_t &ref_quant)
{
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
    this->constants[static_cast<int>(Cons_idx::inv_z)] = 1.0/parcel_height;
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


    // Numerics
#ifdef MET3D
    this->t_end_prime = input.t_end_prime;// + input.start_time;;
#else
    this->t_end_prime = input.t_end_prime;
#endif

    this->t_end = this->t_end_prime/ref_quant.tref;
    // Time of the substeps
    this->dt = input.dt_prime/ref_quant.tref;
    this->dt_prime = input.dt_prime;
    this->dt_traject_prime = 20.0;
    this->dt_traject = this->dt_traject_prime/ref_quant.tref;
    this->num_steps = ceil(this->t_end_prime/this->dt_traject_prime);
    std::cout << "Got num_steps " << this->num_steps << " "
              << this->t_end_prime << "/" << this->dt_traject_prime << "\n";
    // The trajectories from input files are calculated with 20 s timesteps.
    this->num_sub_steps = (floor( this->dt_traject_prime/this->dt ) < 1) ? 1 : floor( this->dt_traject_prime/this->dt );

    // Evaluate the general performance constants
    this->dt_half = this->dt*0.5;
    this->dt_third = this->dt/3.0;
    this->dt_sixth = this->dt/6.0;

    // Scaling factor from input
    this->scaling_fact = input.scaling_fact;

    // Accomodation coefficient
    this->alpha_d = 1.0;

    // // Performance constants for warm cloud; COSMO
    this->a1_scale = 1.0e-3;
    this->a2_scale = 1.72 / pow(get_at(this->constants, Cons_idx::R_a) , 7./8.);
    this->e1_scale = 1.0 / sqrt(get_at(this->constants, Cons_idx::R_a));
    this->e2_scale = 9.1 / pow(get_at(this->constants, Cons_idx::R_a) , 11./16.);
    this->d_scale = ( 130.0*tgamma(4.5) )/( 6.0*(1.0e3)*pow(M_PI*(8.0e6)*get_at(this->constants, Cons_idx::R_a) , 1.0/8.0) );

    // Performance constants for warm cloud; IFS
    // The file constants.h also defines some constants as nar, ...
#ifndef MET3D
    const double Nc = 50; 	// 50 over ocean; 300 over land
    const double F_aut = 1.5;
    const double F_acc = 2.0;
    const double lambda_pp = pow(this->nar * this->ar * tgamma(this->br + 1.0) , this->alpha_r);
#endif
    for(uint32_t i=0; i<4; i++)
    {
        this->constants[static_cast<int>(Cons_idx::a_ccn_1)+i] = a_ccn[i];
        this->constants[static_cast<int>(Cons_idx::b_ccn_1)+i] = b_ccn[i];
        this->constants[static_cast<int>(Cons_idx::c_ccn_1)+i] = c_ccn[i];
        this->constants[static_cast<int>(Cons_idx::d_ccn_1)+i] = d_ccn[i];
    }

    // this->a1_scale = (1350. * F_aut)/pow(Nc , 1.79);
    // this->a2_scale = 67.0 * F_acc;
    // this->e1_scale = 2.0 * M_PI * this->nar * ( (0.78 * tgamma(2.0 - this->nbr))/(lambda_pp*lambda_pp) );
    // this->e2_scale = this->scaling_fact * 2.0 * M_PI * this->nar * 0.31
    //     * pow(this->cr/get_at(this->constants, Particle_cons_idx::mu), 0.5) * pow(this->Sc, 1.0/3.0) * pow(this->rho0, 0.25)
    //     * (tgamma(this->epsilonr + this->nbr)/pow(lambda_pp ,this->epsilonr));
    // this->d_scale = 4.0e-3;

    if(nuc_type == 6)
    {
        this->constants[static_cast<int>(Cons_idx::na_dust)] = na_dust;
        this->constants[static_cast<int>(Cons_idx::na_soot)] = na_soot;
        this->constants[static_cast<int>(Cons_idx::na_orga)] = na_orga;
    } else if(nuc_type == 7 || nuc_type == 5)
    {
        // Standard values
        this->constants[static_cast<int>(Cons_idx::na_dust)] = na_dust_2;
        this->constants[static_cast<int>(Cons_idx::na_soot)] = na_soot_2;
        this->constants[static_cast<int>(Cons_idx::na_orga)] = na_orga_2;
    } else if(nuc_type == 8)
    {
        this->constants[static_cast<int>(Cons_idx::na_dust)] = na_dust_3;
        this->constants[static_cast<int>(Cons_idx::na_soot)] = na_soot_3;
        this->constants[static_cast<int>(Cons_idx::na_orga)] = na_orga_3;
    }

    // Inflow from above
    // this->B_prime = 0.0; //1.0e-7;

    // // Exponents of the cloud model
    // // COSMO
    this->constants[static_cast<int>(Cons_idx::gamma)] = 1.0;
    this->constants[static_cast<int>(Cons_idx::betac)] = 1.0;
    this->constants[static_cast<int>(Cons_idx::betar)] = 7./8.;
    this->constants[static_cast<int>(Cons_idx::delta1)] = 0.5;
    this->constants[static_cast<int>(Cons_idx::delta2)] = 11./16.;
    this->constants[static_cast<int>(Cons_idx::zeta)] = 9./8.;
#if defined(RK4ICE) || defined(RK4NOICE)
    // Exponents of the cloud model
    // IFS
    // get_at(this->constants, Particle_cons_idx::gamma) = 2.47;
    // get_at(this->constants, Cons_idx::betac) = 1.15;
    // get_at(this->constants, Cons_idx::betar) = 1.15;
    // get_at(this->constants, Cons_idx::delta1) = 2.0/( this->br + 1.0 - this->nbr );
    // get_at(this->constants, Cons_idx::delta2) = ( 0.5*this->dr + 2.5 - this->nbr )/( this->br + 1.0 - this->nbr );
    // get_at(this->constants, Cons_idx::zeta) = 1.0;

    // ==================================================
    // Set rain constants
    // See init_2mom_scheme_once in mo_2mom_mcrph_main.f90
    // ==================================================
    // Cosmo5 although nue1nue1 would be okay too, I guess
    //// Cloud
#ifdef SB_SHAPE
    this->cloud.constants[static_cast<int>(Particle_cons_idx::nu)] = 1;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::mu)] = 1;
#else
    this->cloud.constants[static_cast<int>(Particle_cons_idx::nu)] = 0;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0;
#endif
    this->cloud.constants[static_cast<int>(Particle_cons_idx::max_x)] = 2.6e-10;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 4.2e-15;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 1.24e-1;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.333333;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 3.75e5;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.666667;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 1.0;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.0;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-6;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 1.0e-5;
    this->cloud.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0
        / get_at(this->cloud.constants, Particle_cons_idx::cap);
    this->cloud.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->cloud, 1);
    this->cloud.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->cloud, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->cloud.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->cloud, 2);

    this->setup_cloud_autoconversion(this->cloud);
    setup_bulk_sedi(this->cloud);

    //// Rain
#ifdef SB_SHAPE
    this->rain.constants[static_cast<int>(Particle_cons_idx::nu)] = -2.0/3.0; // SB: -2/3 COSMO: 0.0
#else
     this->rain.constants[static_cast<int>(Particle_cons_idx::nu)] = 0; // SB: -2/3 COSMO: 0.0
#endif
    this->rain.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0; // SB: 1/3 COMSO: 1.0/3.0
    this->rain.constants[static_cast<int>(Particle_cons_idx::max_x)] = 3.0e-6;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 2.6e-10;
    this->rain.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 1.24e-1;
    this->rain.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.333333;
    this->rain.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 114.0137;
    this->rain.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.234370;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;

    // From rainSBBcoeffs
    this->rain.constants[static_cast<int>(Particle_cons_idx::alpha)] = 9.292;
    this->rain.constants[static_cast<int>(Particle_cons_idx::beta)] = 9.623;
    this->rain.constants[static_cast<int>(Particle_cons_idx::gamma)] = 6.222e2;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu0)] = 6.0;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu1)] = 3.0e1;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu2)] = 1.0e3;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu3)] = 1.1e-3;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu4)] = 1.0;
    this->rain.constants[static_cast<int>(Particle_cons_idx::cmu5)] = 2.0;
    this->constants[static_cast<int>(Cons_idx::rain_gfak)] = 1.0;

    this->rain.constants[static_cast<int>(Particle_cons_idx::nm1)] =
        (get_at(this->rain.constants, Particle_cons_idx::nu)+1.0)/get_at(this->rain.constants, Particle_cons_idx::mu);
    this->rain.constants[static_cast<int>(Particle_cons_idx::nm2)] =
        (get_at(this->rain.constants, Particle_cons_idx::nu)+2.0)/get_at(this->rain.constants, Particle_cons_idx::mu);
    this->rain.constants[static_cast<int>(Particle_cons_idx::nm3)] =
        (get_at(this->rain.constants, Particle_cons_idx::nu)+3.0)/get_at(this->rain.constants, Particle_cons_idx::mu);
    this->rain.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 20.0;
    this->rain.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.1;
    this->table_r1.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(this->rain.constants, Particle_cons_idx::nm1).getValue());
    this->table_r2.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(this->rain.constants, Particle_cons_idx::nm2).getValue());
    this->table_r3.init_gamma_table(n_lookup, n_lookup_hr_dummy,
        get_at(this->rain.constants, Particle_cons_idx::nm3).getValue());
    this->rain.constants[static_cast<int>(Particle_cons_idx::g1)] = this->table_r1.igf[this->table_r1.n_bins-1];
    this->rain.constants[static_cast<int>(Particle_cons_idx::g2)] = this->table_r2.igf[this->table_r2.n_bins-1];
    this->rain.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(this->rain.constants, Particle_cons_idx::cap);
    this->rain.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->rain, 1);
    this->rain.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->rain, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->rain.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->rain, 2);
    setup_bulk_sedi(this->rain);

    //// Graupel
    this->graupel.constants[static_cast<int>(Particle_cons_idx::nu)] = 1.0; // SB
    this->graupel.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0; // SB
    this->graupel.constants[static_cast<int>(Particle_cons_idx::max_x)] = 5.0e-4;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 1.0e-9;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 1.42e-1;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.314;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 86.89371;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.268325;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 30.0;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.1;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 100.0e-6;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-6;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::s_vel)] = 0.0;

    this->graupel.constants[static_cast<int>(Particle_cons_idx::nm1)] =
        (get_at(this->graupel.constants, Particle_cons_idx::nu)+1.0)/get_at(this->graupel.constants, Particle_cons_idx::mu);
    this->graupel.constants[static_cast<int>(Particle_cons_idx::nm2)] =
        (get_at(this->graupel.constants, Particle_cons_idx::nu)+2.0)/get_at(this->graupel.constants, Particle_cons_idx::mu);
    codi::RealReverse a = (get_at(this->graupel.constants, Particle_cons_idx::nu)+1.0)/get_at(this->graupel.constants, Particle_cons_idx::mu);
    this->table_g1.init_gamma_table(n_lookup, n_lookup_hr_dummy, get_at(this->graupel.constants, Particle_cons_idx::nm1).getValue());
    a = (get_at(this->graupel.constants, Particle_cons_idx::nu)+2.0)/get_at(this->graupel.constants, Particle_cons_idx::mu);
    this->table_g2.init_gamma_table(n_lookup, n_lookup_hr_dummy, get_at(this->graupel.constants, Particle_cons_idx::nm2).getValue());
    this->graupel.constants[static_cast<int>(Particle_cons_idx::g1)] = this->table_g1.igf[this->table_g1.n_bins-1];
    this->graupel.constants[static_cast<int>(Particle_cons_idx::g2)] = this->table_g2.igf[this->table_g2.n_bins-1];
    this->graupel.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(this->graupel.constants, Particle_cons_idx::cap);
    this->graupel.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = 1.0;
    this->graupel.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->graupel, 1);
    this->graupel.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->graupel, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->graupel.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->graupel, 2);
    setup_bulk_sedi(this->graupel);

    //// Hail
    this->hail.constants[static_cast<int>(Particle_cons_idx::nu)] = 1.0;
    this->hail.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0;
    this->hail.constants[static_cast<int>(Particle_cons_idx::max_x)] = 5.0e-4;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 2.6e-9;
    this->hail.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 0.1366;
    this->hail.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 1.0/3.0;
    this->hail.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 39.3;
    this->hail.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.166667;
    this->hail.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    this->hail.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    this->hail.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    this->hail.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 30.0;
    this->hail.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.1;
    this->hail.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = 1.0;
    this->hail.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 100.0e-6;
    this->hail.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-6;
    this->hail.constants[static_cast<int>(Particle_cons_idx::s_vel)] = 0.0;
    this->hail.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(this->hail.constants, Particle_cons_idx::cap);
    this->hail.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = 1.0;
    this->hail.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->hail, 1);
    this->hail.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->hail, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->hail.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->hail, 2);
    setup_bulk_sedi(this->hail);

    //// Ice
#ifdef SB_SHAPE
    this->ice.constants[static_cast<int>(Particle_cons_idx::nu)] = 1.0;
#else
    this->ice.constants[static_cast<int>(Particle_cons_idx::nu)] = 0.0; // COSMO 0.0, SB: 1.0
#endif
    this->ice.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0;
    this->ice.constants[static_cast<int>(Particle_cons_idx::max_x)] = 1.0e-5;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 1.0e-12;
    this->ice.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 0.835;
    this->ice.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.39;
    this->ice.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 2.77e1;
    this->ice.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.21579;
    this->ice.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    this->ice.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    this->ice.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    this->ice.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 3.0;
    this->ice.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.0;
    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = 0.8;
    this->ice.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 150.0e-6;
    this->ice.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-5;
    this->ice.constants[static_cast<int>(Particle_cons_idx::s_vel)] = 0.05;
    this->ice.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(this->ice.constants, Particle_cons_idx::cap);
    this->ice.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = 0.80;
    this->ice.constants[static_cast<int>(Particle_cons_idx::a_f)] = vent_coeff_a(this->ice, 1);
    this->ice.constants[static_cast<int>(Particle_cons_idx::b_f)] = vent_coeff_b(this->ice, 1)
        * pow(get_at(this->constants, Cons_idx::N_sc), get_at(this->constants, Cons_idx::n_f))
        / sqrt(get_at(this->constants, Cons_idx::kin_visc_air));
    this->ice.constants[static_cast<int>(Particle_cons_idx::c_z)] = moment_gamma(this->ice, 2);
    setup_bulk_sedi(this->ice);

    //// Snow
#ifdef SB_SHAPE
    this->snow.constants[static_cast<int>(Particle_cons_idx::nu)] = 1.0; // COSMO: 0.0, SB 1.0
    this->snow.constants[static_cast<int>(Particle_cons_idx::mu)] = 1.0/3.0; // COSMO 0.5, SB: 1.0/3.0
#else
    this->snow.constants[static_cast<int>(Particle_cons_idx::nu)] = 0.0; // COSMO: 0.0, SB 1.0
    this->snow.constants[static_cast<int>(Particle_cons_idx::mu)] = 0.5; // COSMO 0.5, SB: 1.0/3.0
#endif
    this->snow.constants[static_cast<int>(Particle_cons_idx::max_x)] = 2.0e-5;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_act)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_homo)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_nuc_hetero)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_melt)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_evap)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_freezing)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_depo)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_collision)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_collection)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_conversion)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_sedimentation)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::min_x_riming)] = 1.0e-10;
    this->snow.constants[static_cast<int>(Particle_cons_idx::a_geo)] = 2.4;
    this->snow.constants[static_cast<int>(Particle_cons_idx::b_geo)] = 0.455;
    this->snow.constants[static_cast<int>(Particle_cons_idx::a_vel)] = 8.8;
    this->snow.constants[static_cast<int>(Particle_cons_idx::b_vel)] = 0.15;
    this->snow.constants[static_cast<int>(Particle_cons_idx::a_ven)] = 0.78;
    this->snow.constants[static_cast<int>(Particle_cons_idx::b_ven)] = 0.308;
    this->snow.constants[static_cast<int>(Particle_cons_idx::cap)] = 2.0;
    this->snow.constants[static_cast<int>(Particle_cons_idx::vsedi_max)] = 3.0;
    this->snow.constants[static_cast<int>(Particle_cons_idx::vsedi_min)] = 0.1;
    this->snow.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = 0.8;
    this->snow.constants[static_cast<int>(Particle_cons_idx::d_crit_c)] = 150.0e-6;
    this->snow.constants[static_cast<int>(Particle_cons_idx::q_crit_c)] = 1.0e-5;
    this->snow.constants[static_cast<int>(Particle_cons_idx::s_vel)] = 0.25;
    this->snow.constants[static_cast<int>(Particle_cons_idx::c_s)] = 1.0 / get_at(this->snow.constants, Particle_cons_idx::cap);
    this->snow.constants[static_cast<int>(Particle_cons_idx::ecoll_c)] = 0.80;
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
    init_particle_collection_1(this->graupel, this->rain, this->coeffs_grr); // Cosmo uses 2, ICON uses 1
    init_particle_collection_1(this->hail, this->cloud, this->coeffs_hcr);
    init_particle_collection_1(this->graupel, this->cloud, this->coeffs_gcr);
    init_particle_collection_1(this->snow, this->ice, this->coeffs_sic);
    init_particle_collection_1(this->hail, this->ice, this->coeffs_hic);
    init_particle_collection_1(this->graupel, this->ice, this->coeffs_gic);
    init_particle_collection_1(this->hail, this->snow, this->coeffs_hsc);
    init_particle_collection_1(this->graupel, this->snow, this->coeffs_gsc);

    // Setup graupel, snow and ice selfcollection

    this->graupel.constants[static_cast<int>(Particle_cons_idx::sc_coll_n)] = M_PI/8.0
        * (2.0*coll_delta_11(this->graupel, this->graupel, 0)
           + coll_delta_12(this->graupel, this->graupel, 0))
        * sqrt((2.0*coll_theta_11(this->graupel, this->graupel, 0)
           - coll_theta_12(this->graupel, this->graupel, 0)));

    this->snow.constants[static_cast<int>(Particle_cons_idx::sc_delta_n)] = (2.0*coll_delta_11(this->snow, this->snow, 0)
        + coll_delta_12(this->snow, this->snow, 0));
    this->snow.constants[static_cast<int>(Particle_cons_idx::sc_theta_n)] = (2.0*coll_theta_11(this->snow, this->snow, 0)
        - coll_theta_12(this->snow, this->snow, 0));

    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_delta_n)] = coll_delta_11(this->ice, this->ice, 0)
        + coll_delta_12(this->ice, this->ice, 0)
        + coll_delta_22(this->ice, this->ice, 0);
    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_delta_q)] = coll_delta_11(this->ice, this->ice, 0)
        + coll_delta_12(this->ice, this->ice, 1)
        + coll_delta_22(this->ice, this->ice, 1);
    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_theta_n)] = coll_theta_11(this->ice, this->ice, 0)
        - coll_theta_12(this->ice, this->ice, 0)
        + coll_theta_22(this->ice, this->ice, 0);
    this->ice.constants[static_cast<int>(Particle_cons_idx::sc_theta_q)] = coll_theta_11(this->ice, this->ice, 0)
        - coll_theta_12(this->ice, this->ice, 1)
        + coll_theta_22(this->ice, this->ice, 1);
#endif
}

void model_constants_t::print()
{
#ifdef SILENT_MODE
    return;
#endif
  std::cout << "\nModel constants:\n"
        << "----------------\n"
        << "Final integration time = " << this->t_end_prime << " seconds\n"
        << "Nondimensional final integration time = " << this->t_end << "\n"
        << "Timestep = " << this->dt_prime << " seconds\n"
        << "Nondimensional timestep = " << this->dt << "\n"
        << "Number of iterations = " << this->num_steps << "\n"
        << "Number of substeps = " << this->num_sub_steps << "\n"
        << "a1_scale = " << this->a1_scale << "\n"
        << "a2_scale = " << this->a2_scale << "\n"
        << "e1_scale = " << this->e1_scale << "\n"
        << "e2_scale = " << this->e2_scale << "\n"
        << "d_scale = " << this->d_scale << "\n"
        << "Scaling factor = " << this->scaling_fact << "\n";

    for(auto const &t: table_param)
    {
        std::cout << t.first << " = " << get_at(this->constants, t.second) << "\n";
    }
    std::cout << std::endl << std::flush;
}