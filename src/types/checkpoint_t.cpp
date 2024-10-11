#include "include/types/checkpoint_t.h"

checkpoint_t::checkpoint_t() {
}

template<class float_t>
checkpoint_t::checkpoint_t(
    const model_constants_t<float_t> &cc,
    const std::vector<float_t> &y,
    const std::vector<segment_t> &segments,
    const input_parameters_t &input,
    const double &current_time) {
    this->create_checkpoint(cc, y, segments, input, current_time);
}

template<class float_t>
checkpoint_t::checkpoint_t(
    model_constants_t<float_t> &cc,
    const std::vector<float_t> &y,
    const std::vector<segment_t> &segments,
    const input_parameters_t &input,
    const double &current_time,
    const uint32_t &id,
    const uint64_t &ens_id,
    const uint64_t &n_trajs,
    const double duration) {
    auto old_id = cc.ensemble_id;
    cc.ensemble_id = ens_id;
    auto old_tid = cc.traj_id;
    cc.traj_id = id;
    auto old_n = cc.n_trajs;
    cc.n_trajs = n_trajs;
    auto old_desc = cc.ens_desc;
    auto old_t_end_prime = cc.t_end_prime;
    cc.t_end_prime -= input.delay_time_store;
    auto old_num_steps = cc.num_steps;
    if (duration > 0) {
        cc.num_steps = duration;
    } else {
        cc.num_steps -= input.delay_time_store/cc.dt_prime;
    }

    this->create_checkpoint(cc, y, segments, input, current_time);
    cc.ensemble_id = old_id;
    cc.traj_id = old_tid;
    cc.n_trajs = old_n;
    cc.ens_desc = old_desc;
    cc.t_end_prime = old_t_end_prime;
    cc.num_steps = old_num_steps;
}

template<class float_t>
void checkpoint_t::create_checkpoint(
    const model_constants_t<float_t> &cc,
    const std::vector<float_t> &y,
    const std::vector<segment_t> &segments,
    const input_parameters_t &input,
    const double &current_time) {
    // First we add the ensemble configuration
    checkpoint["segments"] = segments;
    // configuration from input_parameters_t
    input.to_json(checkpoint["input"], current_time);
    // Model constants
    checkpoint["model constants"] = cc;
    // Current status of y
    std::map<uint32_t, double> param_map;
    for (uint32_t i=0; i < num_comp; i++) {
        param_map[i] = y[i].getValue();
    }
    checkpoint["Output Parameters"] = param_map;
}

template<class float_t>
int checkpoint_t::load_checkpoint(
    model_constants_t<float_t> &cc,
    std::vector<double> &y,
    std::vector<segment_t> &segments,
    input_parameters_t &input) {
    if (checkpoint.empty()) return CHECKPOINT_LOAD_ERR;
    // Parse the input parameters
    input = checkpoint["input"];
    // Parse the model constants
    cc.from_json(checkpoint["model constants"]);
    cc.setup_dependent_model_constants();
    // Parse the segments and store which parameters had been perturbed
    std::string ens_desc;
    segments.clear();
    for (const auto &s_config : checkpoint["segments"]) {
        segment_t s;
        s.from_json(s_config, cc);
        segments.push_back(s);
    }
    cc.ens_desc += ens_desc;
    std::map<uint32_t, double> param_map = checkpoint["Output Parameters"];
    for (auto const& i : param_map) {
        y[i.first] = i.second;
    }
    return 0;
}

template<class float_t>
int checkpoint_t::load_checkpoint(
    model_constants_t<float_t> &cc,
    std::vector<double> &y,
    std::vector<segment_t> &segments,
    input_parameters_t &input,
    output_handle_t &out_handler) {
    int err = this->load_checkpoint(cc, y, segments, input);
    out_handler.flushed_snapshots = cc.done_steps;
    out_handler.traj = cc.traj_id;
    out_handler.ens = cc.ensemble_id;
    return err;
}

template<class float_t>
int checkpoint_t::load_checkpoint(
    const std::string &filename,
    model_constants_t<float_t> &cc,
    std::vector<double> &y,
    std::vector<segment_t> &segments,
    input_parameters_t &input) {
    std::ifstream i(filename);
    i >> checkpoint;
    return this->load_checkpoint(cc, y, segments, input);
}

void checkpoint_t::print_checkpoint() {
    if (checkpoint.empty()) {
        return;
    }
    std::cout << checkpoint;
}

void checkpoint_t::send_checkpoint(
    const int send_id) {
    s = checkpoint.dump();
    SUCCESS_OR_DIE(
        MPI_Isend(
            s.c_str(),
            s.size(),
            MPI_CHAR,
            send_id,
            CHECKPOINT_MESSAGE,
            MPI_COMM_WORLD,
            &request));
}

bool checkpoint_t::receive_checkpoint() {
    MPI_Status status;
    int count;
    int got_something;
    SUCCESS_OR_DIE(
        MPI_Iprobe(
            MPI_ANY_SOURCE,
            CHECKPOINT_MESSAGE,
            MPI_COMM_WORLD,
            &got_something,
            &status));

    if (!got_something) return got_something;
    SUCCESS_OR_DIE(
        MPI_Get_count(
            &status,
            MPI_CHAR,
            &count));
    char *buff = new char[count];
    SUCCESS_OR_DIE(
        MPI_Irecv(
            buff,
            count,
            MPI_CHAR,
            status.MPI_SOURCE,
            CHECKPOINT_MESSAGE,
            MPI_COMM_WORLD,
            &request));
    MPI_Wait(&request, &status);
    std::istringstream stream(std::string(buff, count));
    try {
        checkpoint.clear();
        stream >> checkpoint;
    } catch (...) {
        std::cout << "receive checkpoint failed: \n"
        << std::string(buff, count) << "\ncount: " << count << "\n"
        << "from " << status.MPI_SOURCE
        << " with tag " << status.MPI_TAG << "\n";
        SUCCESS_OR_DIE(PERTURB_ERR);
    }
    delete [] buff;
    return got_something;
}

bool checkpoint_t::checkpoint_available() const {
    return !checkpoint.empty();
}


template checkpoint_t::checkpoint_t<codi::RealReverseIndex>(
    model_constants_t<codi::RealReverseIndex>&,
    const std::vector<codi::RealReverseIndex>&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&,
    const uint32_t&,
    const uint64_t&,
    const uint64_t&,
    const double);

    template checkpoint_t::checkpoint_t<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&,
    const uint32_t&,
    const uint64_t&,
    const uint64_t&,
    const double);

template checkpoint_t::checkpoint_t<codi::RealReverseIndex>(
    const model_constants_t<codi::RealReverseIndex>&,
    const std::vector<codi::RealReverseIndex>&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template checkpoint_t::checkpoint_t<codi::RealForwardVec<num_par_init> >(
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template void checkpoint_t::create_checkpoint<codi::RealReverseIndex>(
    const model_constants_t<codi::RealReverseIndex>&,
    const std::vector<codi::RealReverseIndex>&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template void checkpoint_t::create_checkpoint<codi::RealForwardVec<num_par_init> >(
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template int checkpoint_t::load_checkpoint<codi::RealReverseIndex>(
    const std::string&,
    model_constants_t<codi::RealReverseIndex>&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&);

template int checkpoint_t::load_checkpoint<codi::RealReverseIndex>(
    model_constants_t<codi::RealReverseIndex>&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&,
    output_handle_t&);

template int checkpoint_t::load_checkpoint<codi::RealReverseIndex>(
    model_constants_t<codi::RealReverseIndex>&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&);

template int checkpoint_t::load_checkpoint<codi::RealForwardVec<num_par_init> >(
    const std::string&,
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&);

template int checkpoint_t::load_checkpoint<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&,
    output_handle_t&);

template int checkpoint_t::load_checkpoint<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&);
