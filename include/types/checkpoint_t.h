#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include "include/misc/error.h"
#include "include/misc/general.h"
#include "include/types/param_t.h"
#include "include/types/segment_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/model_constants_t.h"
#include "include/types/reference_quantities_t.h"

namespace pt = boost::property_tree;

struct checkpoint_t
{
  private:
    pt::ptree checkpoint;

  public:
    checkpoint_t();

    template<class float_t>
    checkpoint_t(model_constants_t &cc,
        const std::vector<float_t> &y,
        std::vector<segment_t> &segments,
        const input_parameters_t &input,
        const double &current_time);

    template<class float_t>
    void create_checkpoint(model_constants_t &cc,
        const std::vector<float_t> &y,
        std::vector<segment_t> &segments,
        const input_parameters_t &input,
        const double &current_time);

    /**
     * Reads a checkpoint file including perturbed parameters if any, current
     * time step, new id as string and sets all objects such as
     * ref_quant, model_constants_t etc.
     */
    template<class float_t>
    int load_checkpoint(const std::string &filename,
        model_constants_t &cc,
        std::vector<float_t> &y,
        std::vector<segment_t> &segments,
        input_parameters_t &input,
        const reference_quantities_t &ref_quant);
    /**
     * Same as above with an already loaded checkpoint instead of reading a
     * file.
     */
    template<class float_t>
    int load_checkpoint(model_constants_t &cc,
        std::vector<float_t> &y,
        std::vector<segment_t> &segments,
        input_parameters_t &input,
        const reference_quantities_t &ref_quant);

    /**
     * Store all data in a property tree and write it as a json file to disk.
     */
    template<class float_t>
    void write_checkpoint(std::string &filename,
        model_constants_t &cc,
        const std::vector<float_t> &y,
        std::vector<segment_t> &segments,
        const input_parameters_t &input,
        const double &current_time);
    /**
     * Write checkpoint to disk as a json file. A property tree must have been
     * created before or nothing will be written.
     */
    void write_checkpoint(std::string &filename,
        model_constants_t &cc,
        std::vector<segment_t> &segments);

    void send_checkpoint();
    void receive_checkpoint();

    // /**
    //  * Getter for the underlying property tree
    //  */
    // pt::ptree get_checkpoint();
    /**
     * Check if a checkpoint had been stored already.
     */
    bool checkpoint_available() const;
};