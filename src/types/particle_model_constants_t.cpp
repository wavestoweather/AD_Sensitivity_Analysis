#include "include/types/particle_model_constants_t.h"


template<class float_t>
particle_model_constants_t<float_t>::particle_model_constants_t() {
    constants.resize(static_cast<int>(Particle_cons_idx::n_items));
    // This makes debugging easier, so please leave it.
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
void particle_model_constants_t<codi::RealReverse>::get_gradient(
    std::array<double, num_par> &out_vec,
    uint32_t &idx,
    const double &ref_value) const {

    const uint32_t start_idx = idx;
    for (auto &c : this->constants) {
#ifdef DEVELOP
        std::cout << "get_gradient particle size uncert = " << uncertainty.size()
                  << " idx: " << idx << " start_idx: " << start_idx << " diff " << idx - start_idx << "\n";
#endif
        out_vec[idx] = c.getGradient() * uncertainty[idx-start_idx] * ref_value;
        idx++;
    }
}


template<class float_t>
int particle_model_constants_t<float_t>::from_json(
    const nlohmann::json& j) {
    int err = 0;
    if (j.find("perturbed") != j.end()) {
        std::map<uint32_t, double> perturbed;
        j.at("perturbed").get_to(perturbed);
        perturbed_idx.clear();
        for (auto const& pert : perturbed) {
            if (pert.first >= this->constants.size())
                err = MODEL_CONS_CHECKPOINT_ERR;
            this->constants[pert.first] = pert.second;
            perturbed_idx.push_back(pert.first);
        }
    }
    return err;
}


template<class float_t>
void particle_model_constants_t<float_t>::print(
    const std::string &title,
    std::ostream &os) {
#ifdef SILENT_MODE
    return;
#endif

    os << title << "\n";
    for (auto const &t : table_particle_param) {
        os << t.first << " = " << get_at(this->constants, t.second) << "\n";
    }
    os << std::endl << std::flush;
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
