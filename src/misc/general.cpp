#include "include/misc/general.h"

void printy(const double y[]) {
#ifdef SILENT_MODE
    return;
#endif
    std::cout << "\n";
    for (int ii = 0 ; ii < num_comp ; ii++) {
        std::cout << "y[" << std::to_string(ii) << "] = " << y[ii] << "\n";
    }
    std::cout << std::endl;
}

/**
 * This function computes
 * \f[ \vec{y} = \vec{v1} + a*\vec{v2} \f]
 * for vectors \f$\vec{y}, \vec{v1}, \vec{v2}\f$ and a scalar \f$a\f$.
 *
 * @param y On out: the result of \f4\vec{y} = \vec{v1} + a*\vec{v2}\f$
 * @param v1 \f$\vec{v1}\f$
 * @param v2 \f$\vec{v2}\f$
 * @param a \f$a\f$
 */
void v1pav2(double y[],
    const double v1[],
    const double v2[],
    const double a) {

    for (int ii = 0 ; ii < num_comp ; ii++) {
        y[ii] = v1[ii] + a*v2[ii];
    }
}

/**
 * Print parameters for given particle.
 *
 * @param pc Particle constants
 * @param title Name of particle
 */
void print_particle_params(
    const particle_model_constants_t &pc,
    const std::string title) {
#ifdef SILENT_MODE
    return;
#endif
    std::cout << title << "\n"
        << "a_geo = " << get_at(pc.constants, Particle_cons_idx::a_geo) << "\n"
        << "b_geo = " << get_at(pc.constants, Particle_cons_idx::b_geo) << "\n"
        << "min_x = " << get_at(pc.constants, Particle_cons_idx::min_x) << "\n"
        << "max_x = " << get_at(pc.constants, Particle_cons_idx::max_x) << "\n"
        << "sc_theta_q = " << get_at(pc.constants, Particle_cons_idx::sc_theta_q) << "\n"
        << "sc_delta_q = " << get_at(pc.constants, Particle_cons_idx::sc_delta_q) << "\n"
        << "sc_theta_n = " << get_at(pc.constants, Particle_cons_idx::sc_theta_n) << "\n"
        << "sc_delta_n = " << get_at(pc.constants, Particle_cons_idx::sc_delta_n) << "\n"
        << "s_vel = " << get_at(pc.constants, Particle_cons_idx::s_vel) << "\n"
        << "a_vel = " << get_at(pc.constants, Particle_cons_idx::a_vel) << "\n"
        << "b_vel = " << get_at(pc.constants, Particle_cons_idx::b_vel) << "\n"
        << "rho_v = " << get_at(pc.constants, Particle_cons_idx::rho_v) << "\n"
        << "c_z = " << get_at(pc.constants, Particle_cons_idx::c_z) << "\n"
        << "sc_coll_n = " << get_at(pc.constants, Particle_cons_idx::sc_coll_n) << "\n"
        << "nu = " << get_at(pc.constants, Particle_cons_idx::nu) << "\n"
        << "mu = " << get_at(pc.constants, Particle_cons_idx::mu) << "\n"
        << "q_crit_c = " << get_at(pc.constants, Particle_cons_idx::q_crit_c) << "\n"
        << "d_crit_c = " << get_at(pc.constants, Particle_cons_idx::d_crit_c) << "\n"
        << "ecoll_c = " << get_at(pc.constants, Particle_cons_idx::ecoll_c) << "\n"
        << "cap = " << get_at(pc.constants, Particle_cons_idx::cap) << "\n"
        << "a_ven = " << get_at(pc.constants, Particle_cons_idx::a_ven) << "\n"
        << "b_ven = " << get_at(pc.constants, Particle_cons_idx::b_ven) << "\n"
        << "c_s = " << get_at(pc.constants, Particle_cons_idx::c_s) << "\n"
        << "a_f = " << get_at(pc.constants, Particle_cons_idx::a_f) << "\n"
        << "b_f = " << get_at(pc.constants, Particle_cons_idx::b_f) << "\n"
        << "b_f (COSMO variant) = "
        << get_at(pc.constants, Particle_cons_idx::b_f) / pow(N_sc, n_f) * sqrt(kin_visc_air) << "\n"
        << "alfa_n = " << get_at(pc.constants, Particle_cons_idx::alfa_n) << "\n"
        << "alfa_q = " << get_at(pc.constants, Particle_cons_idx::alfa_q) << "\n"
        << "lambda = " << get_at(pc.constants, Particle_cons_idx::lambda) << "\n"
        << "vsedi_min = " << get_at(pc.constants, Particle_cons_idx::vsedi_min) << "\n"
        << "vsedi_max = " << get_at(pc.constants, Particle_cons_idx::vsedi_max) << "\n"
        << "g1 = " << get_at(pc.constants, Particle_cons_idx::g1) << "\n"
        << "g2 = " << get_at(pc.constants, Particle_cons_idx::g2) << "\n"
        << "\n";
}

/**
 * Print reference quantities temperature, pressure,
 * mixing-ratio, vertical velocity, time.
 *
 * @param ref Struct with reference quantities.
 */
void print_reference_quantities(const reference_quantities_t &ref) {
#ifdef SILENT_MODE
    return;
#endif
  std::cout << "\nReference quantities\n"
        << "--------------------\n"
        << "Temperature: " << ref.Tref << " Kelvin\n"
        << "Pressure: " << ref.pref << " Pascal\n"
        << "Mixing-ratio: " << ref.qref << "\n"
        << "Vertical velocity: " << ref.wref << " meter per second\n"
        << "Time: " << ref.tref << " Second\n"
        << std::endl << std::flush;
}

/**
 * Print segments for ensembles.
 *
 * @param segments Vector of segments.
 */
void print_segments(const std::vector<segment_t> &segments) {
#ifdef SILENT_MODE
    return;
#endif
    std::cout << "\nSegments for ensembles:\n"
              << "----------------\n";
    uint32_t i = 1;
    uint32_t total_members = 1;
    for (auto &segment : segments) {
        std::cout << "----Segment No. " << i << "\n"
                  << "----Number of members: " << segment.n_members << "\n"
                  << "----Number of perturbed parameters: " << segment.params.size() << "\n"
                  << "----Number of segments defined with this configuration: " << segment.n_segments << "\n";
        if (segment.method != -1)
            std::cout << "----Method to determine segment start: " << segment.method << "\n";
        if (segment.value_name != -1)
            std::cout << "----Parameter that has to reach given value for ensemble start: "
                      << segment.value_name << "\n";
        if (!std::isnan(segment.value))
            std::cout << "----Value for parameter to start an ensemble: " << segment.value << "\n";
        if (segment.out_param != -1)
            std::cout << "----Parameter is a sensitivity to output parameter: " << segment.out_param << "\n";

        uint32_t j = 1;
        for (auto &param : segment.params) {
            std::cout << "--------Param No. " << j << "\n"
                      << "------------Parameter to perturb: " << param.name << "\n"
                      << "------------Output parameter type: " << param.out_name << "\n"
                      << "------------Mean for normal distribution used in perturbing: " << param.mean << "\n"
                      << "------------Is parameter specific for particle: " << param.particle_param << "\n";
            if (!isnan(param.sigma))
                std::cout << "------------Variance for normal distribution: " << param.sigma << "\n";
            if (!isnan(param.sigma_perc))
                std::cout << "------------Variance in percentage: " << param.sigma_perc << "\n";

            j++;
        }
        total_members *= pow(segment.n_members, segment.n_segments);
        i++;
    }
    std::cout << "Total number of trajectories (worst case): " << total_members << "\n\n";
}
