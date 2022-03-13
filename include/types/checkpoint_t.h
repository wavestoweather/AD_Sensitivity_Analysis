#pragma once

#include <mpi.h>

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include <nlohmann/json.hpp>

#include "include/misc/error.h"
#include "include/misc/general.h"
#include "include/types/param_t.h"
#include "include/types/segment_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/model_constants_t.h"
#include "include/types/output_handle_t.h"
#include "include/types/reference_quantities_t.h"

struct checkpoint_t {
 private:
    nlohmann::json checkpoint;

 public:
    checkpoint_t();

    template<class float_t>
    checkpoint_t(
        const model_constants_t<float_t> &cc,
        const std::vector<float_t> &y,
        const std::vector<segment_t> &segments,
        const input_parameters_t &input,
        const double &current_time);

    template<class float_t>
    checkpoint_t(
        model_constants_t<float_t> &cc,
        const std::vector<float_t> &y,
        const std::vector<segment_t> &segments,
        const input_parameters_t &input,
        const double &current_time,
        const uint32_t &id,
        const uint64_t &ens_id,
        const uint64_t &n_trajs,
        const double duration = 0);

    template<class float_t>
    void create_checkpoint(
        const model_constants_t<float_t> &cc,
        const std::vector<float_t> &y,
        const std::vector<segment_t> &segments,
        const input_parameters_t &input,
        const double &current_time);

    /**
     * Reads a checkpoint file including perturbed parameters if any, current
     * time step, new id as string and sets all objects such as
     * ref_quant, model_constants_t etc.
     */
    template<class float_t>
    int load_checkpoint(
        const std::string &filename,
        model_constants_t<float_t> &cc,
        std::vector<double> &y,
        std::vector<segment_t> &segments,
        input_parameters_t &input,
        const reference_quantities_t &ref_quant);
    /**
     * Same as above with an already loaded checkpoint instead of reading a
     * file.
     */
    template<class float_t>
    int load_checkpoint(
        model_constants_t<float_t> &cc,
        std::vector<double> &y,
        std::vector<segment_t> &segments,
        input_parameters_t &input,
        const reference_quantities_t &ref_quant);
    /**
     * Same as above with an already loaded checkpoint instead of reading a
     * file and adjusting time steps for flushing output.
     */
    template<class float_t>
    int load_checkpoint(
        model_constants_t<float_t> &cc,
        std::vector<double> &y,
        std::vector<segment_t> &segments,
        input_parameters_t &input,
        const reference_quantities_t &ref_quant,
        output_handle_t &out_handler);

    /**
     * Store all data in a property tree and write it as a json file to disk.
     */
    template<class float_t>
    void write_checkpoint(
        std::string &filename,
        model_constants_t<float_t> &cc,
        const std::vector<float_t> &y,
        std::vector<segment_t> &segments,
        const input_parameters_t &input,
        const double &current_time);

    /**
     * Write checkpoint to disk as a json file. A property tree must have been
     * created before or nothing will be written.
     */
    template<class float_t>
    void write_checkpoint(
        std::string &filename,
        const model_constants_t<float_t> &cc,
        const std::vector<segment_t> &segments);

    /**
     * Print the checkpoint to std::out. For debugging purpose.
     */
    void print_checkpoint();

    /**
     * Send a checkpoint to the specified process.
     */
    void send_checkpoint(const int send_id);
    /**
     *
     * @return true if checkpoint is received, false otherwise.
     */
    bool receive_checkpoint();

    /**
     * Check if a checkpoint had been stored already.
     */
    bool checkpoint_available() const;
};
