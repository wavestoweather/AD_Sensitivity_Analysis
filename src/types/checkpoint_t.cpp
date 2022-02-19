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
    // std::string f = "tmp";
    // this->write_checkpoint(f, cc, segments);
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
    // pt::ptree segment_tree;
    // for (auto &s : segments)
    //     s.put(segment_tree);


    // checkpoint.add_child("segments", segment_tree);
    // std::vector<nlohmann::json> segments_json;
    // for (auto &s : segments)
    //     segments_json.emplace_back(s);
    checkpoint["segments"] = segments;
    // configuration from input_parameters_t
    // input.put(checkpoint, current_time);
    // auto tmp_time = input.current_time;
    // input.current_time = current_time;
    input.to_json(checkpoint["input"], current_time);
    // checkpoint["input"] = input;
    // input.current_time = tmp_time;
    // Model constants
    // cc.put(checkpoint);
    checkpoint["model constants"] = cc;
    // Current status of y
    // pt::ptree output_parameters;
    // for (uint32_t i=0; i < num_comp; i++) {
    //     output_parameters.put(std::to_string(i), y[i].getValue());
    // }
    // checkpoint.add_child("Output Parameters", output_parameters);

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
    input_parameters_t &input,
    const reference_quantities_t &ref_quant) {
    if (checkpoint.empty()) return CHECKPOINT_LOAD_ERR;
    // Parse the input parameters
    // SUCCESS_OR_DIE(input.from_pt(checkpoint));
    input = checkpoint["input"];

    // Parse the model constants
    // SUCCESS_OR_DIE(cc.from_pt(checkpoint));
    cc.from_json(checkpoint["model constants"]);
    cc.setup_model_constants(input, ref_quant);

    // Parse the segments and store which parameters had been perturbed
    // in ens_descr
    std::string ens_desc;
    segments.clear();
    // for (auto &it : checkpoint.get_child("segments")) {
    //     segment_t segment;
    //     SUCCESS_OR_DIE(segment.from_pt(it.second, cc));

    //     if (segment.activated) {
    //         segment.perturb(cc, ref_quant, input, ens_desc);
    //     }
    //     segments.push_back(segment);
    // }
    for (const auto &s_config : checkpoint["segments"]) {
        segment_t s;
        s.from_json(s_config, cc);
        segments.push_back(s);
    }
    // segments = checkpoint["segments"];
    // for (auto &s : segments) {
    //     for (auto &param : s.params) SUCCESS_OR_DIE(param.check_name(cc));
    // }

    cc.ens_desc += ens_desc;
    // for (auto &it : checkpoint.get_child("Output Parameters")) {
    //     try {
    //         y[std::stoi(it.first)] = std::stod(it.second.data());
    //     } catch (const std::out_of_range& e) {
    //         // out_of_range means underflow with some libraries
    //         y[std::stoi(it.first)] = 0;
    //     }
    // }
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
    const reference_quantities_t &ref_quant,
    output_handle_t &out_handler) {
    int err = this->load_checkpoint(cc, y, segments, input, ref_quant);
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
    input_parameters_t &input,
    const reference_quantities_t &ref_quant) {
    std::ifstream i(filename);
    i >> checkpoint;
    // boost::property_tree::read_json(filename, checkpoint);
    return this->load_checkpoint(cc, y, segments, input, ref_quant);
}

template<class float_t>
void checkpoint_t::write_checkpoint(
    std::string &filename,
    model_constants_t<float_t> &cc,
    const std::vector<float_t> &y,
    std::vector<segment_t> &segments,
    const input_parameters_t &input,
    const double &current_time) {
    this->create_checkpoint(cc, y, segments, input, current_time);
    this->write_checkpoint(filename, cc, segments);
}

template<class float_t>
void checkpoint_t::write_checkpoint(
    std::string &filename,
    const model_constants_t<float_t> &cc,
    const std::vector<segment_t> &segments) {
    if (checkpoint.empty()) {
        return;
    }
    uint64_t i = 0;
    std::string actual_filename = filename + "/checkpoint_id" + cc.id + "_0000.json";
    while (exists(actual_filename)) {
        i++;
        if (i < 10)
            actual_filename = filename + "/checkpoint_id" + cc.id + "_000" + std::to_string(i) + ".json";
        else if (i < 100)
            actual_filename = filename + "/checkpoint_id" + cc.id + "_00" + std::to_string(i) + ".json";
        else if (i < 1000)
            actual_filename = filename + "/checkpoint_id" + cc.id + "_0" + std::to_string(i) + ".json";
        else
            actual_filename = filename + "/checkpoint_id" + cc.id + "_" + std::to_string(i) + ".json";
    }
    std::fstream outstream(actual_filename, std::ios::out);
    filename = actual_filename;
    outstream << checkpoint << std::endl;
    // pt::write_json(outstream, checkpoint);
    outstream.close();
    // deactivate all segments, so we know, another instance is going
    // to process this
    // for (auto &s : segments)
    //     s.deactivate(true);
    // cc.ensemble_id++;
}

void checkpoint_t::print_checkpoint() {
    if (checkpoint.empty()) {
        return;
    }
    // pt::write_json(std::cout, checkpoint);
    std::cout << checkpoint;
}

void checkpoint_t::send_checkpoint(
    const int send_id) {
    MPI_Request request;
    // std::stringstream ss;
    // pt::json_parser::write_json(ss, checkpoint);
    // std::string s = ss.str();
    std::string s = checkpoint.dump();

    // std::cout << "send checkpoint to " << send_id
    //     << " with size " << s.size() << "\n" << s << "\n";
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
    // get source tag
    char *buff = new char[count];
    SUCCESS_OR_DIE(
        MPI_Recv(
            buff,
            count,
            MPI_CHAR,
            status.MPI_SOURCE,
            CHECKPOINT_MESSAGE,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE));
    // std::string s(buff, count);
    // delete [] buff;
    std::istringstream stream(std::string(buff, count));

    // boost::iostreams::stream<boost::iostreams::array_source> stream(
    //     buff, count);
    // std::cout << "size: " << count << "\n" << buff << "\n";
    try {
        checkpoint.clear();
        stream >> checkpoint;
        // std::cout << "put everything in j1\n";
        // std::cout << j1;
        // pt::read_json(stream, checkpoint);
        // std::cout << "Got something\n" << std::flush;
    } catch (...) {
        std::cout << "receive checkpoint failed: \n"
        << buff << "\ncount: " << count << "\n"
        << "from " << status.MPI_SOURCE << "\n";
    }
    delete [] buff;
    return got_something;
}

bool checkpoint_t::checkpoint_available() const {
    return !checkpoint.empty();
}

template checkpoint_t::checkpoint_t<codi::RealReverse>(
    model_constants_t<codi::RealReverse>&,
    const std::vector<codi::RealReverse>&,
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

template checkpoint_t::checkpoint_t<codi::RealReverse>(
    const model_constants_t<codi::RealReverse>&,
    const std::vector<codi::RealReverse>&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template checkpoint_t::checkpoint_t<codi::RealForwardVec<num_par_init> >(
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template void checkpoint_t::create_checkpoint<codi::RealReverse>(
    const model_constants_t<codi::RealReverse>&,
    const std::vector<codi::RealReverse>&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template void checkpoint_t::create_checkpoint<codi::RealForwardVec<num_par_init> >(
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    const std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template void checkpoint_t::write_checkpoint<codi::RealReverse>(
    std::string&,
    model_constants_t<codi::RealReverse>&,
    const std::vector<codi::RealReverse>&,
    std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template void checkpoint_t::write_checkpoint<codi::RealForwardVec<num_par_init> >(
    std::string&,
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    std::vector<segment_t>&,
    const input_parameters_t&,
    const double&);

template void checkpoint_t::write_checkpoint<codi::RealReverse>(
    std::string&,
    const model_constants_t<codi::RealReverse>&,
    const std::vector<segment_t>&);

template void checkpoint_t::write_checkpoint<codi::RealForwardVec<num_par_init> >(
    std::string&,
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector<segment_t>&);

template int checkpoint_t::load_checkpoint<codi::RealReverse>(
    const std::string&,
    model_constants_t<codi::RealReverse>&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&,
    const reference_quantities_t&);

template int checkpoint_t::load_checkpoint<codi::RealReverse>(
    model_constants_t<codi::RealReverse>&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&,
    const reference_quantities_t&,
    output_handle_t&);

template int checkpoint_t::load_checkpoint<codi::RealReverse>(
    model_constants_t<codi::RealReverse>&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&,
    const reference_quantities_t&);

template int checkpoint_t::load_checkpoint<codi::RealForwardVec<num_par_init> >(
    const std::string&,
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&,
    const reference_quantities_t&);

template int checkpoint_t::load_checkpoint<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&,
    const reference_quantities_t&,
    output_handle_t&);

template int checkpoint_t::load_checkpoint<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    std::vector<double>&,
    std::vector<segment_t>&,
    input_parameters_t&,
    const reference_quantities_t&);
