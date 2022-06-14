#include <include/interface/physics_t.h>

physics_t::physics_t() : cc("", false) {
    yold.resize(num_comp);
    ynew.resize(num_comp);
    ytmp.resize(num_comp);
    k.resize(num_comp);
    gradients_.resize(num_comp + static_cast<uint32_t>(Init_cons_idx::n_items));
    set_ref_quants(1.0e-6, 1.0e5, 1.0, 1.0, 1.0, 1.0, 1.0);
    setup_model_constants(30.0, 0.1);
    load_lookup_table(cc.ltabdminwgg);
    for (auto &y : yold) y = 0;
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
    // Does it make any difference if we use that before or after registering the  data?
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

    codi::RealReverse::Tape& tape = codi::RealReverse::getTape();
    tape.setActive();
    cc.register_input(tape);

    yold[p_idx] *= ref_quant.pref;
    yold[w_idx] *= ref_quant.wref;
    yold[T_idx] *= ref_quant.Tref;
    yold[qv_idx] *= ref_quant.qref;
    yold[qc_idx] *= ref_quant.qref;
    yold[Nc_idx] *= ref_quant.Nref;

    cc.setup_dependent_model_constants();

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

    // Explicit calculation of saturation
    T_prime = ynew[T_idx];
    p_prime = ynew[p_idx];
    qv_prime = ynew[qv_idx];
    ynew[S_idx] = convert_qv_to_S(p_prime, T_prime, qv_prime,
                                  get_at(cc.constants, Cons_idx::p_sat_low_temp),
                                  get_at(cc.constants, Cons_idx::p_sat_const_a),
                                  get_at(cc.constants, Cons_idx::T_sat_low_temp),
                                  get_at(cc.constants, Cons_idx::p_sat_const_b),
                                  get_at(cc.constants, Cons_idx::Epsilon));

    yold[p_idx] /= ref_quant.pref;
    yold[w_idx] /= ref_quant.wref;
    yold[T_idx] /= ref_quant.Tref;
    yold[qv_idx] /= ref_quant.qref;
    yold[qc_idx] /= ref_quant.qref;
    yold[Nc_idx] /= ref_quant.Nref;
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
    res[qv_idx] = ynew[qv_idx].getValue() * ref_quant.qref;
    res[qc_idx] = ynew[qc_idx].getValue() * ref_quant.qref;
    res[Nc_idx] = ynew[Nc_idx].getValue() * ref_quant.Nref;
    res[S_idx] = ynew[S_idx].getValue();
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

    codi::RealReverse::Tape& tape = codi::RealReverse::getTape();
    tape.setActive();
    cc.register_input(tape);
    cc.setup_dependent_model_constants();

    yold[T_idx] *= ref_quant.Tref;
    yold[qg_idx] *= ref_quant.qref;
    yold[Ng_idx] *= ref_quant.Nref;

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

    ynew[T_idx] /= ref_quant.Tref;
    ynew[qg_idx] /= ref_quant.qref;
    ynew[Ng_idx] /= ref_quant.Nref;
    ynew[qr_idx] /= ref_quant.qref;
    ynew[Nr_idx] /= ref_quant.Nref;

    // Get the gradients
    cc.get_gradients(ynew, gradients_, tape, ref_quant);

    for (auto ii = 0 ; ii < num_comp; ii++) {
        for (auto i = 0; i < num_par; i++) {
            gradients[i + ii*num_par] = gradients_[ii][i];
        }
    }
    // copy results
    res[T_idx] = ynew[T_idx].getValue() * ref_quant.Tref;
    res[qg_idx] = ynew[qg_idx].getValue() * ref_quant.qref;
    res[Ng_idx] = ynew[Ng_idx].getValue() * ref_quant.Nref;
    res[qr_idx] = ynew[qr_idx].getValue() * ref_quant.qref;
    res[Nr_idx] = ynew[Nr_idx].getValue() * ref_quant.Nref;
}
