#include <include/interface/physics_t.h>

physics_t::physics_t(std::string table_path) : cc("", false) {
    yold.resize(num_comp);
    ynew.resize(num_comp);
    ytmp.resize(num_comp);
    k.resize(num_comp);
    gradients_.resize(num_comp + static_cast<uint32_t>(Init_cons_idx::n_items));
    set_ref_quants(1.0e-6, 1.0e5, 1.0, 1.0, 1.0, 1.0, 1.0);
    setup_model_constants(30.0, 0.1);
    load_lookup_table(cc.ltabdminwgg, table_path);
    for (auto &y : yold) y = 0;
}


codi::RealReverse::Tape& physics_t::prepare_call() {
    codi::RealReverse::Tape& tape = codi::RealReverse::getTape();
    tape.setActive();
    cc.register_input(tape);
    cc.setup_dependent_model_constants();

    yold[T_idx] *= ref_quant.Tref;
    yold[p_idx] *= ref_quant.pref;
    yold[w_idx] *= ref_quant.wref;
    yold[qg_idx] *= ref_quant.qref;
    yold[Ng_idx] *= ref_quant.Nref;
    yold[qr_idx] *= ref_quant.qref;
    yold[Nr_idx] *= ref_quant.Nref;
    yold[qc_idx] *= ref_quant.qref;
    yold[Nc_idx] *= ref_quant.Nref;
    yold[qs_idx] *= ref_quant.qref;
    yold[Ns_idx] *= ref_quant.Nref;
    yold[qi_idx] *= ref_quant.qref;
    yold[Ni_idx] *= ref_quant.Nref;
    yold[qh_idx] *= ref_quant.qref;
    yold[Nh_idx] *= ref_quant.Nref;
    yold[qv_idx] *= ref_quant.qref;

    codi::RealReverse p_prime = yold[p_idx];
    codi::RealReverse T_prime = yold[T_idx];
    codi::RealReverse qv_prime = yold[qv_idx];
    yold[S_idx] = convert_qv_to_S(
            p_prime,
            T_prime,
            qv_prime,
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon));

    codi::RealReverse rho_inter = log(compute_rhoh(p_prime, T_prime, yold[S_idx],
                                         get_at(cc.constants, Cons_idx::p_sat_low_temp),
                                         get_at(cc.constants, Cons_idx::p_sat_const_a),
                                         get_at(cc.constants, Cons_idx::T_sat_low_temp),
                                         get_at(cc.constants, Cons_idx::p_sat_const_b),
                                         get_at(cc.constants, Cons_idx::R_a),
                                         get_at(cc.constants, Cons_idx::R_v)) / get_at(cc.constants, Cons_idx::rho_0));
    cc.cloud.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
            exp(-get_at(cc.constants, Cons_idx::rho_vel_c) * rho_inter);
    cc.rain.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
            exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    cc.graupel.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
            exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    cc.hail.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
            exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    cc.ice.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
            exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    cc.snow.constants[static_cast<int>(Particle_cons_idx::rho_v)] =
            exp(-get_at(cc.constants, Cons_idx::rho_vel) * rho_inter);
    return tape;
}

void physics_t::finish_call(
    codi::RealReverse::Tape& tape,
    double *res,
    double *gradients) {

    auto T_prime = ynew[T_idx];
    auto p_prime = ynew[p_idx];
    auto qv_prime = ynew[qv_idx];
    ynew[S_idx] = convert_qv_to_S(p_prime, T_prime, qv_prime,
                                  get_at(cc.constants, Cons_idx::p_sat_low_temp),
                                  get_at(cc.constants, Cons_idx::p_sat_const_a),
                                  get_at(cc.constants, Cons_idx::T_sat_low_temp),
                                  get_at(cc.constants, Cons_idx::p_sat_const_b),
                                  get_at(cc.constants, Cons_idx::Epsilon));

    ynew[T_idx] /= ref_quant.Tref;
    ynew[p_idx] /= ref_quant.pref;
    ynew[w_idx] /= ref_quant.wref;
    ynew[qg_idx] /= ref_quant.qref;
    ynew[Ng_idx] /= ref_quant.Nref;
    ynew[qr_idx] /= ref_quant.qref;
    ynew[Nr_idx] /= ref_quant.Nref;
    ynew[qc_idx] /= ref_quant.qref;
    ynew[Nc_idx] /= ref_quant.Nref;
    ynew[qs_idx] /= ref_quant.qref;
    ynew[Ns_idx] /= ref_quant.Nref;
    ynew[qi_idx] /= ref_quant.qref;
    ynew[Ni_idx] /= ref_quant.Nref;
    ynew[qh_idx] /= ref_quant.qref;
    ynew[Nh_idx] /= ref_quant.Nref;
    ynew[qv_idx] /= ref_quant.qref;

    // Get the gradients
    cc.get_gradients(ynew, gradients_, tape, ref_quant);

    for (auto ii = 0 ; ii < num_comp; ii++) {
        for (auto i = 0; i < num_par; i++) {
            gradients[i + ii*num_par] = gradients_[ii][i];
        }
    }
    // copy results
    res[p_idx] = ynew[p_idx].getValue() * ref_quant.pref;
    res[w_idx] = ynew[w_idx].getValue() * ref_quant.wref;
    res[T_idx] = ynew[T_idx].getValue() * ref_quant.Tref;
    res[qg_idx] = ynew[qg_idx].getValue() * ref_quant.qref;
    res[Ng_idx] = ynew[Ng_idx].getValue() * ref_quant.Nref;
    res[qr_idx] = ynew[qr_idx].getValue() * ref_quant.qref;
    res[Nr_idx] = ynew[Nr_idx].getValue() * ref_quant.Nref;
    res[qv_idx] = ynew[qv_idx].getValue() * ref_quant.qref;
    res[qc_idx] = ynew[qc_idx].getValue() * ref_quant.qref;
    res[Nc_idx] = ynew[Nc_idx].getValue() * ref_quant.Nref;
    res[qs_idx] = ynew[qs_idx].getValue() * ref_quant.qref;
    res[Ns_idx] = ynew[Ns_idx].getValue() * ref_quant.Nref;
    res[qi_idx] = ynew[qi_idx].getValue() * ref_quant.qref;
    res[Ni_idx] = ynew[Ni_idx].getValue() * ref_quant.Nref;
    res[qh_idx] = ynew[qh_idx].getValue() * ref_quant.qref;
    res[Nh_idx] = ynew[Nh_idx].getValue() * ref_quant.Nref;
    res[S_idx] = ynew[S_idx].getValue();
}

void physics_t::set_ref_quants(
        const double qref,
        const double pref,
        const double wref,
        const double Tref,
        const double zref,
        const double Nref,
        const double timeref) {
    ref_quant.qref = qref;
    ref_quant.pref = pref;
    ref_quant.wref = wref;
    ref_quant.Tref = Tref;
    ref_quant.zref = zref;
    ref_quant.Nref = Nref;
    ref_quant.tref = timeref;
}

void physics_t::setup_model_constants(
    const double dt_prime,
    const double uncertainty_perc) {

    input.dt_prime = dt_prime;
    uncertainty_scale = uncertainty_perc;
    setup_model_constants();
}

void physics_t::setup_model_constants(
    const double uncertainty_perc) {

    uncertainty_scale = uncertainty_perc;
    setup_model_constants();
}

void physics_t::setup_model_constants() {
    auto_type = input.auto_type;
    cc.setup_model_constants(input, ref_quant);
    cc.set_uncertainty(uncertainty_scale);
}


void physics_t::set_limits(
    std::vector<codi::RealReverse> &y,
    model_constants_t<codi::RealReverse> &cc) {
    if (nuc_type > 0)
        y[Nc_idx] = min(min(max(y[Nc_idx], y[qc_idx]/get_at(cc.cloud.constants, Particle_cons_idx::max_x)),
                            y[qc_idx]/get_at(cc.cloud.constants, Particle_cons_idx::min_x)), 5e9);
    y[Nr_idx] = min(max(y[Nr_idx], y[qr_idx]/get_at(cc.rain.constants, Particle_cons_idx::max_x)),
                    y[qr_idx]/get_at(cc.rain.constants, Particle_cons_idx::min_x));
    y[Ni_idx] = min(max(y[Ni_idx], y[qi_idx]/get_at(cc.ice.constants, Particle_cons_idx::max_x)),
                    y[qi_idx]/get_at(cc.ice.constants, Particle_cons_idx::min_x));
    y[Ns_idx] = min(max(y[Ns_idx], y[qs_idx]/get_at(cc.snow.constants, Particle_cons_idx::max_x)),
                    y[qs_idx]/get_at(cc.snow.constants, Particle_cons_idx::min_x));
    y[Ng_idx] = min(max(y[Ng_idx], y[qg_idx]/get_at(cc.graupel.constants, Particle_cons_idx::max_x)),
                    y[qg_idx]/get_at(cc.graupel.constants, Particle_cons_idx::min_x));
    y[Nh_idx] = min(max(y[Nh_idx], y[qh_idx]/get_at(cc.hail.constants, Particle_cons_idx::max_x)),
                    y[qh_idx]/get_at(cc.hail.constants, Particle_cons_idx::min_x));
    // Set everything negative to zero
    y[qv_idx] = (y[qv_idx] < 0) ? 0 : y[qv_idx];
    y[qc_idx] = (y[qc_idx] < 0) ? 0 : y[qc_idx];
    y[qr_idx] = (y[qr_idx] < 0) ? 0 : y[qr_idx];
    y[qs_idx] = (y[qs_idx] < 0) ? 0 : y[qs_idx];
    y[qi_idx] = (y[qi_idx] < 0) ? 0 : y[qi_idx];
    y[qg_idx] = (y[qg_idx] < 0) ? 0 : y[qg_idx];
    y[qh_idx] = (y[qh_idx] < 0) ? 0 : y[qh_idx];
    y[Nc_idx] = (y[Nc_idx] < 0) ? 0 : y[Nc_idx];
    y[Nr_idx] = (y[Nr_idx] < 0) ? 0 : y[Nr_idx];
    y[Ns_idx] = (y[Ns_idx] < 0) ? 0 : y[Ns_idx];
    y[Ni_idx] = (y[Ni_idx] < 0) ? 0 : y[Ni_idx];
    y[Ng_idx] = (y[Ng_idx] < 0) ? 0 : y[Ng_idx];
    y[Nh_idx] = (y[Nh_idx] < 0) ? 0 : y[Nh_idx];
}


void physics_t::compute_nondimensional_effects(
    std::vector<codi::RealReverse> &res,
    const codi::RealReverse &p,
    const codi::RealReverse &w,
    const codi::RealReverse &w_prime,
    const codi::RealReverse &T,
    const codi::RealReverse &qv,
    const codi::RealReverse &qv_prime) {
    codi::RealReverse C1 = (ref_quant.tref*get_at(cc.constants, Cons_idx::gravity_acc)*ref_quant.wref)
                 / (get_at(cc.constants, Cons_idx::R_a)*ref_quant.Tref);
    codi::RealReverse C2 = ((1.0-get_at(cc.constants, Cons_idx::Epsilon))*ref_quant.qref)
                 / get_at(cc.constants, Cons_idx::Epsilon);
    // Calculate pressure and temperature change in non-prime directly
    // Pressure change from ascent (C1), change in partial pressure from water vapor.
    res[p_idx] = -(C1/(1.0 + C2*(qv/(1.0 + qv_prime))))*((p*w)/T);
    res[T_idx] += (res[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                   + res[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp)
                   - get_at(cc.constants, Cons_idx::gravity_acc)
                   * w_prime/get_at(cc.constants, Cons_idx::cp)) / ref_quant.Tref;
    res[w_idx] += get_at(cc.constants, Cons_idx::dw)/ref_quant.wref;
    res[z_idx] += w_prime/ref_quant.zref;
}


void physics_t::py_ccn_act_hande_akm(
    const double p,
    const double w,
    const double T,
    const double qv,
    const double qc,
    const double Nc,
    double *res,
    double *gradients) {
#ifdef DEVELOP
    std::ofstream of;
    of.open("physics_t.log");
#endif
    for (auto &y : yold) y = 0;
    yold[p_idx] = p;
    yold[w_idx] = w;
    yold[T_idx] = T;
    yold[qv_idx] = qv;
    yold[qc_idx] = qc;
    yold[Nc_idx] = Nc;
#ifdef DEVELOP
    of << "params:\n"
       << "p: " << p
       << ", w: " << w
       << ", T: " << T
       << ", qv: " << qv
       << ", qc: " << qc
       << ", Nc: " << Nc << "\n";

    cc.print(of);
    cc.cloud.print("cloud", of);
    cc.rain.print("rain", of);
    cc.ice.print("ice", of);
    cc.snow.print("snow", of);
    cc.graupel.print("graupel", of);
    cc.hail.print("hail", of);
#endif
    const double EPSILON = 1.0e-20;

    codi::RealReverse::Tape& tape = prepare_call();
#ifdef DEVELOP
    of << "yold:\n"
       << "p: " << yold[p_idx]
       << ", w: " << yold[w_idx]
       << ", T: " << yold[T_idx]
       << ", qv: " << yold[qv_idx]
       << ", qc: " << yold[qc_idx]
       << ", Nc: " << yold[Nc_idx]
       << ", S: " << yold[p_idx] << "\n";
#endif
    // RK4
    for (auto &kv : k) kv = 0;
    for (int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];

    // k, yold, dt_sixth
    ccn_act_hande_akm(yold[p_idx], yold[w_idx], yold[T_idx], yold[qv_idx],
                      yold[qc_idx], yold[Nc_idx], EPSILON, k, cc);
    compute_nondimensional_effects(
        k, yold[p_idx]/ref_quant.pref, yold[w_idx]/ref_quant.wref, yold[w_idx],
        yold[T_idx]/ref_quant.Tref, yold[qv_idx]/ref_quant.qref, yold[qv_idx]);
#ifdef DEVELOP
    of << "After first call\n";
#endif
    for (int ii = 0 ; ii < num_comp ; ii++) {
#ifdef DEVELOP
        of << ii << ", " << yold[ii] << " * " << cc.dt_half << ", " << k[ii] << " * " << cc.dt_sixth << "\n";
#endif
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii];  // Add k1-part to the result
    }
#ifdef DEVELOP
    of << "\n";
#endif
    set_limits(ytmp, cc);

    for (auto &kv : k) kv = 0;
    // k, ytmp, dt_third
    ccn_act_hande_akm(ytmp[p_idx], ytmp[w_idx], ytmp[T_idx], ytmp[qv_idx],
                      ytmp[qc_idx], ytmp[Nc_idx], EPSILON, k, cc);
    compute_nondimensional_effects(
            k, ytmp[p_idx]/ref_quant.pref, ytmp[w_idx]/ref_quant.wref, ytmp[w_idx],
            ytmp[T_idx]/ref_quant.Tref, ytmp[qv_idx]/ref_quant.qref, ytmp[qv_idx]);
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii];  // Add k2-part to the result
    }
    set_limits(ytmp, cc);

    for (auto &kv : k) kv = 0;
    // k, ytmp, dt_third
    ccn_act_hande_akm(ytmp[p_idx], ytmp[w_idx], ytmp[T_idx], ytmp[qv_idx],
                      ytmp[qc_idx], ytmp[Nc_idx], EPSILON, k, cc);
    compute_nondimensional_effects(
            k, ytmp[p_idx]/ref_quant.pref, ytmp[w_idx]/ref_quant.wref, ytmp[w_idx],
            ytmp[T_idx]/ref_quant.Tref, ytmp[qv_idx]/ref_quant.qref, ytmp[qv_idx]);
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt*k[ii];  // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii];  // Add k3-part to the result
    }
    set_limits(ytmp, cc);

    for (auto &kv : k) kv = 0;
    // k, ytmp, dt_sixth
    ccn_act_hande_akm(ytmp[p_idx], ytmp[w_idx], ytmp[T_idx], ytmp[qv_idx],
                      ytmp[qc_idx], ytmp[Nc_idx], EPSILON, k, cc);
    compute_nondimensional_effects(
            k, ytmp[p_idx]/ref_quant.pref, ytmp[w_idx]/ref_quant.wref, ytmp[w_idx],
            ytmp[T_idx]/ref_quant.Tref, ytmp[qv_idx]/ref_quant.qref, ytmp[qv_idx]);
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] += cc.dt_sixth*k[ii];
    }
    set_limits(ynew, cc);

    finish_call(tape, res, gradients);
#ifdef DEVELOP
    of << "Test end\n";
    of << "p: " << res[p_idx]
        << ", w: " << res[w_idx]
        << ", T: " << res[T_idx]
        << ", qv: " << res[qv_idx]
        << ", qc: " << res[qc_idx]
        << ", Nc: " << res[Nc_idx]
        << ", S: " << res[p_idx] << "\n";
    of.close();
#endif
}

void physics_t::py_graupel_melting(
    const double T,
    const double qg,
    const double Ng,
    double *res,
    double *gradients) {

    for (auto &y : yold) y = 0;
    yold[T_idx] = T;
    yold[qg_idx] = qg;
    yold[Ng_idx] = Ng;

    codi::RealReverse::Tape& tape = prepare_call();

    // RK4
    for (auto &kv : k) kv = 0;
    for (int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];
    graupel_melting(yold[qg_idx], yold[Ng_idx], yold[T_idx], k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                   + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii];  // Add k1-part to the result
    }
    set_limits(ytmp, cc);

    for (auto &kv : k) kv = 0;
    graupel_melting(ytmp[qg_idx], ytmp[Ng_idx], ytmp[T_idx], k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii];  // Add k2-part to the result
    }
    set_limits(ytmp, cc);

    for (auto &kv : k) kv = 0;
    graupel_melting(ytmp[qg_idx], ytmp[Ng_idx], ytmp[T_idx], k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt*k[ii];  // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii];  // Add k3-part to the result
    }
    set_limits(ytmp, cc);

    cc.dt_prime = cc.dt_sixth;
    for (auto &kv : k) kv = 0;
    graupel_melting(ytmp[qg_idx], ytmp[Ng_idx], ytmp[T_idx], k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;

    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] += cc.dt_sixth*k[ii];
    }
    set_limits(ynew, cc);

    finish_call(tape, res, gradients);
}


void physics_t::py_saturation_adjust(
    const double p,
    const double T,
    const double qv,
    const double qc,
    double *res,
    double *gradients) {

    for (auto &y : yold) y = 0;
    yold[p_idx] = p;
    yold[T_idx] = T;
    yold[qv_idx] = qv;
    yold[qc_idx] = qc;

    codi::RealReverse::Tape& tape = prepare_call();

    // RK4
    for (auto &kv : k) kv = 0;
    for (int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];
    saturation_adjust(yold[T_idx], yold[p_idx],
                      yold[qv_idx], yold[qc_idx], k, cc);
    compute_nondimensional_effects(
            k, yold[p_idx]/ref_quant.pref, yold[w_idx]/ref_quant.wref, yold[w_idx],
            yold[T_idx]/ref_quant.Tref, yold[qv_idx]/ref_quant.qref, yold[qv_idx]);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii];  // Add k1-part to the result
    }
    set_limits(ytmp, cc);

    // k, ytmp, dt_third
    for (auto &kv : k) kv = 0;
    saturation_adjust(ytmp[T_idx], ytmp[p_idx],
                      ytmp[qv_idx], ytmp[qc_idx], k, cc);
    compute_nondimensional_effects(
            k, ytmp[p_idx]/ref_quant.pref, ytmp[w_idx]/ref_quant.wref, ytmp[w_idx],
            ytmp[T_idx]/ref_quant.Tref, ytmp[qv_idx]/ref_quant.qref, ytmp[qv_idx]);
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii];  // Add k2-part to the result
    }
    set_limits(ytmp, cc);

    for (auto &kv : k) kv = 0;
    // k, ytmp, dt_third
    saturation_adjust(ytmp[T_idx], ytmp[p_idx],
                      ytmp[qv_idx], ytmp[qc_idx], k, cc);
    compute_nondimensional_effects(
            k, ytmp[p_idx]/ref_quant.pref, ytmp[w_idx]/ref_quant.wref, ytmp[w_idx],
            ytmp[T_idx]/ref_quant.Tref, ytmp[qv_idx]/ref_quant.qref, ytmp[qv_idx]);
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii];  // Add k2-part to the result
    }
    set_limits(ytmp, cc);

    for (auto &kv : k) kv = 0;
    // dt_sixth
    saturation_adjust(ytmp[T_idx], ytmp[p_idx],
                      ytmp[qv_idx], ytmp[qc_idx], k, cc);
    compute_nondimensional_effects(
            k, ytmp[p_idx]/ref_quant.pref, ytmp[w_idx]/ref_quant.wref, ytmp[w_idx],
            ytmp[T_idx]/ref_quant.Tref, ytmp[qv_idx]/ref_quant.qref, ytmp[qv_idx]);
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] += cc.dt_sixth*k[ii];
    }
    set_limits(ynew, cc);

    finish_call(tape, res, gradients);
}


void physics_t::py_riming_ice(
    const double qc,
    const double Nc,
    const double qi,
    const double Ni,
    const double qr,
    const double Nr,
    const double T,
    double *res,
    double *gradients) {

    for (auto &y : yold) y = 0;
    yold[qc_idx] = qc;
    yold[Nc_idx] = Nc;
    yold[qi_idx] = qi;
    yold[Ni_idx] = Ni;
    yold[qr_idx] = qr;
    yold[Nr_idx] = Nr;
    yold[T_idx] = T;

    codi::RealReverse::Tape& tape = prepare_call();

    codi::RealReverse rime_rate_qc, rime_rate_nc, rime_rate_qi, rime_rate_qr, rime_rate_nr;
    codi::RealReverse dep_rate_ice = 0.0;
    // RK4
    for (auto &kv : k) kv = 0;
    for (int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];
    riming_cloud_core(yold[qc_idx], yold[Nc_idx], yold[qi_idx], yold[Ni_idx],
                      cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(yold[qr_idx], yold[Nr_idx], yold[qi_idx], yold[Ni_idx],
                     cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(yold[qc_idx], yold[Nc_idx], yold[qr_idx], yold[Nr_idx], yold[qi_idx], yold[Ni_idx],
               dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qi, yold[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii];  // Add k1-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_ice = 0.0;
    for (auto &kv : k) kv = 0;
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qi_idx], ytmp[Ni_idx],
                      cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
                     cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
               dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qi, ytmp[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii];  // Add k2-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_ice = 0.0;
    for (auto &kv : k) kv = 0;
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qi_idx], ytmp[Ni_idx],
                      cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
                     cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
               dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qi, ytmp[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt*k[ii];  // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii];  // Add k3-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_ice = 0.0;
    cc.dt_prime = cc.dt_sixth;
    for (auto &kv : k) kv = 0;
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qi_idx], ytmp[Ni_idx],
                      cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
                     cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
               dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qi, ytmp[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;

    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] += cc.dt_sixth*k[ii];
    }
    set_limits(ynew, cc);

    finish_call(tape, res, gradients);
}


void physics_t::py_riming_snow(
    const double qc,
    const double Nc,
    const double qs,
    const double Ns,
    const double qr,
    const double Nr,
    const double T,
    double *res,
    double *gradients) {

    for (auto &y : yold) y = 0;
    yold[qc_idx] = qc;
    yold[Nc_idx] = Nc;
    yold[qs_idx] = qs;
    yold[Ns_idx] = Ns;
    yold[qr_idx] = qr;
    yold[Nr_idx] = Nr;
    yold[T_idx] = T;

    codi::RealReverse::Tape& tape = prepare_call();

    codi::RealReverse rime_rate_qc, rime_rate_nc, rime_rate_qs, rime_rate_qr, rime_rate_nr;
    codi::RealReverse dep_rate_snow = 0.0;
    // RK4
    for (auto &kv : k) kv = 0;
    for (int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];
    riming_cloud_core(yold[qc_idx], yold[Nc_idx], yold[qs_idx], yold[Ns_idx],
                      cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(yold[qr_idx], yold[Nr_idx], yold[qs_idx], yold[Ns_idx],
                     cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(yold[qc_idx], yold[Nc_idx], yold[qr_idx], yold[Nr_idx], yold[qs_idx], yold[Ns_idx],
               dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qs, yold[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii];  // Add k1-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_snow = 0.0;
    for (auto &kv : k) kv = 0;
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qs_idx], ytmp[Ns_idx],
                      cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                     cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
               dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qs, ytmp[T_idx], cc.dt, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii];  // Add k2-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_snow = 0.0;
    for (auto &kv : k) kv = 0;
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qs_idx], ytmp[Ns_idx],
                      cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                     cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
               dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qs, ytmp[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt*k[ii];  // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii];  // Add k3-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_snow = 0.0;
    cc.dt_prime = cc.dt_sixth;
    for (auto &kv : k) kv = 0;
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qs_idx], ytmp[Ns_idx],
                      cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                     cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
               dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qs, ytmp[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;

    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] += cc.dt_sixth*k[ii];
    }
    set_limits(ynew, cc);

    finish_call(tape, res, gradients);
}

void physics_t::py_riming_with_depo(
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
    double *gradients) {

    for (auto &y : yold) y = 0;
    yold[qv_idx] = qv;
    yold[qc_idx] = qc;
    yold[Nc_idx] = Nc;
    yold[qi_idx] = qi;
    yold[Ni_idx] = Ni;
    yold[qs_idx] = qs;
    yold[Ns_idx] = Ns;
    yold[qr_idx] = qr;
    yold[Nr_idx] = Nr;
    yold[T_idx] = T;
    yold[p_idx] = p;

    const double EPSILON = 1.0e-20;
    codi::RealReverse::Tape& tape = prepare_call();

    codi::RealReverse rime_rate_qc, rime_rate_nc, rime_rate_qi, rime_rate_qs, rime_rate_qr, rime_rate_nr;

    codi::RealReverse dep_rate_ice = 0.0;
    codi::RealReverse dep_rate_snow = 0.0;
    codi::RealReverse p_sat_ice = saturation_pressure_ice(
            yold[T_idx], get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_b));
    codi::RealReverse D_vtp = diffusivity(yold[T_idx], yold[p_idx]);
    codi::RealReverse S = convert_qv_to_S(
            yold[p_idx],
            yold[T_idx],
            yold[qv_idx],
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon));
    codi::RealReverse e_d = compute_pv(
            yold[T_idx], S, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));
    codi::RealReverse S_i =
        (yold[T_idx] < get_at(cc.constants, Cons_idx::T_freeze)) ? e_d / p_sat_ice : codi::RealReverse(1);
    codi::RealReverse s_si = S_i - 1.0;
    // RK4
    for (auto &kv : k) kv = 0;
    for (int ii = 0 ; ii < num_comp ; ii++)
        ynew[ii] = yold[ii];
    vapor_dep_relaxation(yold[qv_idx], yold[qi_idx], yold[Ni_idx],
                         yold[qs_idx], yold[Ns_idx], yold[qg_idx], yold[Ng_idx],
                         yold[qh_idx], yold[Nh_idx], s_si, p_sat_ice, yold[p_idx],
                         yold[T_idx], EPSILON, dep_rate_ice, dep_rate_snow, D_vtp, k, cc);
    riming_cloud_core(yold[qc_idx], yold[Nc_idx], yold[qi_idx], yold[Ni_idx],
                      cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(yold[qr_idx], yold[Nr_idx], yold[qi_idx], yold[Ni_idx],
                     cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(yold[qc_idx], yold[Nc_idx], yold[qr_idx], yold[Nr_idx], yold[qi_idx], yold[Ni_idx],
               dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qi, yold[T_idx], cc.dt_prime, k, cc);
    riming_cloud_core(yold[qc_idx], yold[Nc_idx], yold[qs_idx], yold[Ns_idx],
                      cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(yold[qr_idx], yold[Nr_idx], yold[qs_idx], yold[Ns_idx],
                     cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(yold[qc_idx], yold[Nc_idx], yold[qr_idx], yold[Nr_idx], yold[qs_idx], yold[Ns_idx],
                dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
                rime_rate_qs, yold[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k1 for k2
        ynew[ii] += cc.dt_sixth*k[ii];  // Add k1-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_ice = 0.0;
    dep_rate_snow = 0.0;
    p_sat_ice = saturation_pressure_ice(
            ytmp[T_idx], get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_b));
    D_vtp = diffusivity(ytmp[T_idx], ytmp[p_idx]);
    S = convert_qv_to_S(
            ytmp[p_idx],
            ytmp[T_idx],
            ytmp[qv_idx],
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon));
    e_d = compute_pv(
            ytmp[T_idx], S, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));
    S_i = (ytmp[T_idx] < get_at(cc.constants, Cons_idx::T_freeze)) ? e_d / p_sat_ice : codi::RealReverse(1);
    s_si = S_i - 1.0;
    for (auto &kv : k) kv = 0;
    vapor_dep_relaxation(ytmp[qv_idx], ytmp[qi_idx], ytmp[Ni_idx],
                         ytmp[qs_idx], ytmp[Ns_idx], ytmp[qg_idx], ytmp[Ng_idx],
                         ytmp[qh_idx], ytmp[Nh_idx], s_si, p_sat_ice, ytmp[p_idx],
                         ytmp[T_idx], EPSILON, dep_rate_ice, dep_rate_snow, D_vtp, k, cc);
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qi_idx], ytmp[Ni_idx],
                      cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
                     cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
               dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qi, ytmp[T_idx], cc.dt_prime, k, cc);
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qs_idx], ytmp[Ns_idx],
                      cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                     cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
                rime_rate_qs, ytmp[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt_half*k[ii];  // y_0 + (dt/2)*k2 for k3
        ynew[ii] += cc.dt_third*k[ii];  // Add k2-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_ice = 0.0;
    dep_rate_snow = 0.0;
    p_sat_ice = saturation_pressure_ice(
            ytmp[T_idx], get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_b));
    D_vtp = diffusivity(ytmp[T_idx], ytmp[p_idx]);
    S = convert_qv_to_S(
            ytmp[p_idx],
            ytmp[T_idx],
            ytmp[qv_idx],
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon));
    e_d = compute_pv(
            ytmp[T_idx], S, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));
    S_i = (ytmp[T_idx] < get_at(cc.constants, Cons_idx::T_freeze)) ? e_d / p_sat_ice : codi::RealReverse(1);
    s_si = S_i - 1.0;
    for (auto &kv : k) kv = 0;
    vapor_dep_relaxation(ytmp[qv_idx], ytmp[qi_idx], ytmp[Ni_idx],
                         ytmp[qs_idx], ytmp[Ns_idx], ytmp[qg_idx], ytmp[Ng_idx],
                         ytmp[qh_idx], ytmp[Nh_idx], s_si, p_sat_ice, ytmp[p_idx],
                         ytmp[T_idx], EPSILON, dep_rate_ice, dep_rate_snow, D_vtp, k, cc);
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qi_idx], ytmp[Ni_idx],
                      cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
                     cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
               dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qi, ytmp[T_idx], cc.dt_prime, k, cc);
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qs_idx], ytmp[Ns_idx],
                      cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                     cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
                rime_rate_qs, ytmp[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;
    for (int ii = 0 ; ii < num_comp ; ii++) {
        ytmp[ii] = yold[ii] + cc.dt*k[ii];  // y_0 + dt*k3 for k4
        ynew[ii] += cc.dt_third*k[ii];  // Add k3-part to the result
    }
    set_limits(ytmp, cc);

    dep_rate_ice = 0.0;
    dep_rate_snow = 0.0;
    p_sat_ice = saturation_pressure_ice(
            ytmp[T_idx], get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_ice_const_b));
    D_vtp = diffusivity(ytmp[T_idx], ytmp[p_idx]);
    S = convert_qv_to_S(
            ytmp[p_idx],
            ytmp[T_idx],
            ytmp[qv_idx],
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon));
    e_d = compute_pv(
            ytmp[T_idx], S, get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b));
    S_i = (ytmp[T_idx] < get_at(cc.constants, Cons_idx::T_freeze)) ? e_d / p_sat_ice : codi::RealReverse(1);
    s_si = S_i - 1.0;
    cc.dt_prime = cc.dt_sixth;
    for (auto &kv : k) kv = 0;
    vapor_dep_relaxation(ytmp[qv_idx], ytmp[qi_idx], ytmp[Ni_idx],
                         ytmp[qs_idx], ytmp[Ns_idx], ytmp[qg_idx], ytmp[Ng_idx],
                         ytmp[qh_idx], ytmp[Nh_idx], s_si, p_sat_ice, ytmp[p_idx],
                         ytmp[T_idx], EPSILON, dep_rate_ice, dep_rate_snow, D_vtp, k, cc);
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qi_idx], ytmp[Ni_idx],
                      cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
                     cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qi_idx], ytmp[Ni_idx],
               dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
               rime_rate_qi, ytmp[T_idx], cc.dt_prime, k, cc);
    riming_cloud_core(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qs_idx], ytmp[Ns_idx],
                      cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                     cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(ytmp[qc_idx], ytmp[Nc_idx], ytmp[qr_idx], ytmp[Nr_idx], ytmp[qs_idx], ytmp[Ns_idx],
                dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
                rime_rate_qs, ytmp[T_idx], cc.dt_prime, k, cc);
    k[T_idx] += (k[lat_heat_idx]/get_at(cc.constants, Cons_idx::cp)
                 + k[lat_cool_idx]/get_at(cc.constants, Cons_idx::cp))/ref_quant.Tref;

    for (int ii = 0 ; ii < num_comp ; ii++) {
        ynew[ii] += cc.dt_sixth*k[ii];
    }
    set_limits(ynew, cc);

    finish_call(tape, res, gradients);
}
