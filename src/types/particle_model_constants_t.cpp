#include "include/types/particle_model_constants_t.h"


particle_model_constants_t::particle_model_constants_t() {
    constants.resize(static_cast<int>(Particle_cons_idx::n_items));
    // This makes debugging easier, so pleaase leave it.
    std::fill(constants.begin(), constants.end(), 0);
}


void particle_model_constants_t::register_input(
    codi::RealReverse::TapeType &tape,
    uint32_t &idx) {
    for (auto &c : this->constants) {
        tape.registerInput(c);
        idx++;
    }
}


void particle_model_constants_t::get_gradient(
    std::array<double, num_par> &out_vec,
    uint32_t &idx,
    const bool info) const {

    const uint32_t start_idx = idx;
    for (auto &c : this->constants) {
        out_vec[idx] = c.getGradient() * uncertainty[idx-start_idx];
        idx++;
    }
}


void particle_model_constants_t::put(
    pt::ptree &ptree,
    const std::string &type_name) const {
    if (perturbed_idx.empty())
        return;

    pt::ptree perturbed;

    for (uint32_t idx : perturbed_idx) {
        perturbed.put(std::to_string(idx), constants[idx]);
    }
    pt::ptree perturbed_vals;
    perturbed_vals.add_child("perturbed", perturbed);
    ptree.add_child(type_name, perturbed_vals);
}


int particle_model_constants_t::from_pt(
    pt::ptree &ptree) {
    int err = 0;
    for (auto &it : ptree.get_child("perturbed")) {
        uint32_t idx = std::stoi(it.first);
        this->constants[idx] = it.second.get_value<double>();
        perturbed_idx.push_back(idx);
    }
    return err;
}


void particle_model_constants_t::print(
    const std::string &title) {
#ifdef SILENT_MODE
    return;
#endif

    std::cout << title << "\n";
    for (auto const &t : table_particle_param) {
        std::cout << t.first << " = " << get_at(this->constants, t.second) << "\n";
    }
    std::cout << std::endl << std::flush;
}
