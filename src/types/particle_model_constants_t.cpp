#include "include/types/particle_model_constants_t.h"


template<class float_t>
particle_model_constants_t<float_t>::particle_model_constants_t() {
    constants.resize(static_cast<int>(Particle_cons_idx::n_items));
    // This makes debugging easier, so pleaase leave it.
    std::fill(constants.begin(), constants.end(), 0);
}


template<>
void particle_model_constants_t<codi::RealReverse>::register_input(
    codi::RealReverse::Tape &tape,
    uint32_t &idx) {
    for (auto &c : this->constants) {
        tape.registerInput(c);
        idx++;
    }
}


template<>
void particle_model_constants_t<codi::RealForwardVec<num_par_init> >::register_input(
    codi::RealReverse::Tape &tape,
    uint32_t &idx) {

    // Nothing to be done for forward mode.
}


template<>
void particle_model_constants_t<codi::RealReverse>::get_gradient(
    std::array<double, num_par> &out_vec,
    uint32_t &idx,
    const bool info) const {

    const uint32_t start_idx = idx;
    for (auto &c : this->constants) {
        out_vec[idx] = c.getGradient() * uncertainty[idx-start_idx];
        idx++;
    }
}


template<>
void particle_model_constants_t<codi::RealForwardVec<num_par_init> >::get_gradient(
    std::array<double, num_par> &out_vec,
    uint32_t &idx,
    const bool info) const {

    // Nothing to be done for forward mode
}




template<class float_t>
int particle_model_constants_t<float_t>::from_json(
    const nlohmann::json& j) {
    int err = 0;
    std::map<uint32_t, double> perturbed;
    j.at("perturbed").get_to(perturbed);
    perturbed_idx.clear();
    for (auto const& pert : perturbed) {
        this->constants[pert.first] = pert.second;
        perturbed_idx.push_back(pert.first);
    }
    return err;
}

// template<class float_t>
// void particle_model_constants_t<float_t>::put(
//     pt::ptree &ptree,
//     const std::string &type_name) const {
//     if (perturbed_idx.empty())
//         return;

//     pt::ptree perturbed;

//     for (uint32_t idx : perturbed_idx) {
//         perturbed.put(std::to_string(idx), constants[idx].getValue());
//     }
//     pt::ptree perturbed_vals;
//     perturbed_vals.add_child("perturbed", perturbed);
//     ptree.add_child(type_name, perturbed_vals);
// }


// template<class float_t>
// int particle_model_constants_t<float_t>::from_pt(
//     pt::ptree &ptree) {
//     int err = 0;
//     for (auto &it : ptree.get_child("perturbed")) {
//         uint32_t idx = std::stoi(it.first);
//         this->constants[idx] = it.second.get_value<double>();
//         perturbed_idx.push_back(idx);
//     }
//     return err;
// }


template<class float_t>
void particle_model_constants_t<float_t>::print(
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


template class particle_model_constants_t<codi::RealReverse>;
template class particle_model_constants_t<codi::RealForwardVec<num_par_init> >;

template<class float_t>
void to_json(
    nlohmann::json& j,
    const particle_model_constants_t<float_t>& p) {
    if (p.perturbed_idx.empty())
        return;
    
    std::map<uint32_t, double> perturbed;
    for (uint32_t idx : p.perturbed_idx) {
        perturbed[idx] = p.constants[idx].getValue();
    }
    j["perturbed"] = perturbed;
}

template void to_json<codi::RealReverse>(
    nlohmann::json&,
    const particle_model_constants_t<codi::RealReverse>& p);

template void to_json<codi::RealForwardVec<num_par_init> >(
    nlohmann::json&,
    const particle_model_constants_t<codi::RealForwardVec<num_par_init> >& p);
