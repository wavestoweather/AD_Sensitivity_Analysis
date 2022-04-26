#include "include/types/output_handle_t.h"


output_handle_t::output_handle_t() {
}


template<class float_t>
output_handle_t::output_handle_t(
        const std::string filetype,
        const std::string filename,
        const model_constants_t<float_t> &cc,
        const std::string in_filename,
        const uint32_t write_index,
        const uint32_t snapshot_index,
        const int &rank,
        const int &simulation_mode,
        const bool &initial_cond,
#ifdef COMPRESS_OUTPUT
        const int &n_processes,
#endif
        const double delay_out_time) {
    this->simulation_mode = simulation_mode;
    local_num_comp = cc.local_num_comp;
    local_num_par = cc.local_num_par;
    local_ic_par = cc.local_ic_par;
    this->track_ic = initial_cond;
#ifdef COMPRESS_OUTPUT
    this->n_processes = n_processes;
#endif
    this->setup(filetype, filename, cc,
                in_filename, write_index, snapshot_index,
                rank, delay_out_time);
}


template<class float_t>
output_handle_t::output_handle_t(
    const std::string filetype,
    const std::string filename,
    const model_constants_t<float_t> &cc,
    const std::string in_filename,
    const uint32_t write_index,
    const uint32_t snapshot_index,
    const int &rank,
    const int &simulation_mode,
    const bool &initial_cond,
#ifdef COMPRESS_OUTPUT
    const int &n_processes,
#endif
    const double delay_out_time,
    const std::vector<segment_t> &segments) {

    this->simulation_mode = simulation_mode;
    local_num_comp = cc.local_num_comp;
    local_num_par = cc.local_num_par;
    local_ic_par = cc.local_ic_par;
    this->track_ic = initial_cond;
#ifdef COMPRESS_OUTPUT
    this->n_processes = n_processes;
#endif
    if (simulation_mode == create_train_set && segments.size() > 0) {
        uint64_t param_size = cc.constants.size()
            + cc.cloud.constants.size() + cc.rain.constants.size()
            + cc.ice.constants.size() + cc.snow.constants.size()
            + cc.graupel.constants.size() + cc.hail.constants.size();

        perturbed_idx.resize(param_size);
        std::fill(perturbed_idx.begin(), perturbed_idx.end(), -1);
        // Define an order of the perturbed parameters for each ensemble
        // for the case of create_train_set, that can be reused in the
        // output file.
        // Get the idx of the parameter and its first occurrence in the segments.
        // Some parameters might be perturbed in multiple ensembles.
        n_perturbed_params = 0;
        for (const auto &s : segments) {
            // Get all parameters that shall be perturbed and their index.
            // Is it already in perturbed_idx? Then skip it.
            for (const auto &p : s.params) {
                std::string p_name = p.get_name();
                auto it = std::find(perturbed_names.begin(), perturbed_names.end(), p_name);
                if (it == perturbed_names.end()) {
                    perturbed_names.push_back(p_name);
                    uint64_t param_idx = p.get_idx();
                    if (p.outparam_name == "cloud") {
                        unperturbed_vals.push_back(cc.cloud.constants[param_idx].getValue());
                        param_idx += cc.constants.size();
                    } else if (p.outparam_name == "rain") {
                        unperturbed_vals.push_back(cc.rain.constants[param_idx].getValue());
                        param_idx += cc.constants.size() + cc.cloud.constants.size();
                    } else if (p.outparam_name == "ice") {
                        unperturbed_vals.push_back(cc.ice.constants[param_idx].getValue());
                        param_idx += cc.constants.size() + cc.cloud.constants.size() + cc.rain.constants.size();
                    } else if (p.outparam_name == "snow") {
                        unperturbed_vals.push_back(cc.snow.constants[param_idx].getValue());
                        param_idx += cc.constants.size() + cc.cloud.constants.size() + cc.rain.constants.size()
                                     + cc.ice.constants.size();
                    } else if (p.outparam_name == "graupel") {
                        unperturbed_vals.push_back(cc.graupel.constants[param_idx].getValue());
                        param_idx += cc.constants.size() + cc.cloud.constants.size() + cc.rain.constants.size()
                                     + cc.ice.constants.size() + cc.snow.constants.size();
                    } else if (p.outparam_name == "hail") {
                        unperturbed_vals.push_back(cc.hail.constants[param_idx].getValue());
                        param_idx += cc.constants.size() + cc.cloud.constants.size() + cc.rain.constants.size()
                                     + cc.ice.constants.size() + cc.snow.constants.size()
                                     + cc.graupel.constants.size();
                    } else {
                        unperturbed_vals.push_back(cc.constants[param_idx].getValue());
                    }
                        perturbed_idx[param_idx] = n_perturbed_params;
                    n_perturbed_params++;
                }
            }
        }
    }
    this->setup(filetype, filename, cc,
        in_filename, write_index, snapshot_index,
        rank, delay_out_time);
}


template<class float_t>
void output_handle_t::setup_gradients(
    const model_constants_t<float_t> &cc) {

    if (this->simulation_mode == limited_time_ensembles) {
        std::vector<int> dimid_tmp;
        auto dim_pointer = &dimid[Dim_idx::time_dim];
        int n_dims = 1;
        if (local_num_comp > 1) {
            dimid_tmp.push_back(dimid[Dim_idx::out_param_dim]);
            dimid_tmp.push_back(dimid[Dim_idx::time_dim]);
            dim_pointer = &dimid_tmp[0];
            n_dims = 2;
        }
        define_var_gradients(cc, dim_pointer, n_dims);
    } else if (this->simulation_mode == create_train_set) {
        std::vector<int> dimid_tmp;
        if (local_num_comp > 1) {
            dimid_tmp.push_back(dimid[Dim_idx::out_param_dim]);
        }
        dimid_tmp.push_back(dimid[Dim_idx::trajectory_dim]);
        dimid_tmp.push_back(dimid[Dim_idx::time_dim]);
        define_var_gradients(cc, &dimid_tmp[0], dimid_tmp.size());
    } else {
        auto dim_pointer = &dimid[Dim_idx::out_param_dim];
        int n_dims = 4;
        if (local_num_comp == 1) {
            dim_pointer = &dimid[Dim_idx::ensemble_dim];
            n_dims = 3;
        }
        define_var_gradients(cc, dim_pointer, n_dims);
    }
}


template<class float_t>
void output_handle_t::define_var_gradients(
    const model_constants_t<float_t> &cc,
    const int *dim_pointer,
    const int &n_dims) {

    if (track_ic) {
        // initial condition sensitivity
        for (uint32_t i = 0; i < num_par_init; ++i) {
            if (cc.trace_check(i, 2)) {
                SUCCESS_OR_DIE(nc_def_var(
                       ncid,
                       init_grad_idx[i].c_str(),
                       NC_FLOAT_T,
                       n_dims,
                       dim_pointer,
                       &varid[Var_idx::n_vars + i]));
            }
        }
    } else {
        for (uint32_t i = 0; i < num_par - num_par_init; ++i) {
            if (cc.trace_check(i, false)) {
                SUCCESS_OR_DIE(nc_def_var(
                       ncid,
                       output_grad_idx[i].c_str(),
                       NC_FLOAT_T,
                       n_dims,
                       dim_pointer,
                       &varid[Var_idx::n_vars + i]));
            }
        }
    }
}


template<class float_t>
void output_handle_t::define_vars(const model_constants_t<float_t> &cc) {
    // Define variables for dimensions
    if (simulation_mode == create_train_set) {
        SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "Output_Parameter",
                NC_STRING,
                (local_num_comp > 1),
                &dimid[Dim_idx::out_param_dim],
                &varid[Var_idx::out_param]));
    } else {
        SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "Output_Parameter_ID",
                NC_UINT64,  // type
                (local_num_comp > 1),          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::out_param_dim],    // ignored for ndims = 0
                &varid[Var_idx::out_param]));
    }
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "trajectory",
            NC_UINT64,
            1,
            &dimid[Dim_idx::trajectory_dim],
            &varid[Var_idx::trajectory]));
    if (simulation_mode == create_train_set) {
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "perturbed",
            NC_STRING,
            (n_perturbed_params > 1),
            &dimid[Dim_idx::perturb_param_dim],
            &varid[Var_idx::perturbed]));
    }
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "ensemble",
            NC_UINT64,
            1,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::ensemble]));

    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "time",
            NC_FLOAT_T,
           1,
           &dimid[Dim_idx::time_dim],
           &varid[Var_idx::time]));

    // model state
    for (uint32_t i = 0; i < num_comp; ++i)
        SUCCESS_OR_DIE(
            nc_def_var(
                ncid,
                output_par_idx[i].c_str(),
                NC_FLOAT_T,
                3,
                &dimid[Dim_idx::ensemble_dim],
                &varid[i]));

    // gradients
    setup_gradients(cc);

    // the rest
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "time_after_ascent",
            NC_FLOAT_T,
           3,
           &dimid[Dim_idx::ensemble_dim],
           &varid[Var_idx::time_ascent]));
#if !defined B_EIGHT
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "conv_400",
            NC_BYTE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::conv_400]));
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "conv_600",
            NC_BYTE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::conv_600]));
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "slan_400",
            NC_BYTE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::slan_400]));
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "slan_600",
            NC_BYTE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::slan_600]));
#endif
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "lat",
            NC_FLOAT_T,
           3,
           &dimid[Dim_idx::ensemble_dim],
           &varid[Var_idx::lat]));
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "lon",
            NC_FLOAT_T,
           3,
           &dimid[Dim_idx::ensemble_dim],
           &varid[Var_idx::lon]));
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "step",
            NC_UINT64,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::step]));

    // Phase of the trajectory
    SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "phase",
            NC_UINT64,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::phase]));

    if (simulation_mode == create_train_set) {
        std::vector<int> dimid_tmp;
        dimid_tmp.push_back(dimid[Dim_idx::ensemble_dim]);
        dimid_tmp.push_back(dimid[Dim_idx::time_dim]);
        dimid_tmp.push_back(dimid[Dim_idx::perturb_param_dim]);
        SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "perturbation_value",
               NC_FLOAT_T,
                3,
                &dimid_tmp[0],
                &varid[Var_idx::perturbation_value]));
    }
}


template<class float_t>
void output_handle_t::set_attributes(
    const model_constants_t<float_t> &cc,
    const std::string in_filename) {
#ifdef DEVELOP
    std::cout << "attempt to open " << in_filename << "\n";
#endif
    // Add attributes; read attributes from in_filename first
    int in_ncid;
    SUCCESS_OR_DIE(nc_open(in_filename.c_str(), NC_NOWRITE, &in_ncid));
#ifdef DEVELOP
    std::cout << "opened\n" << std::flush;
#endif
    // get amount of attributes
    int n_atts;
    SUCCESS_OR_DIE(nc_inq(in_ncid, NULL, NULL, &n_atts, NULL));

    // for every attribute, get the name, type and length, and the values
    for (int i=0; i < n_atts; i++) {
        char att_name[NC_MAX_NAME];
        SUCCESS_OR_DIE(nc_inq_attname(in_ncid, NC_GLOBAL, i, att_name));
        nc_type att_type;
        size_t att_len;
        SUCCESS_OR_DIE(nc_inq_att(
                in_ncid, NC_GLOBAL, att_name, &att_type, &att_len));
        std::vector<nc_vlen_t> att_val(att_len);
        SUCCESS_OR_DIE(nc_get_att(
                in_ncid, NC_GLOBAL, att_name, att_val.data()));
        SUCCESS_OR_DIE(nc_put_att(
                ncid, NC_GLOBAL, att_name, att_type, att_len, att_val.data()));
    }
    // time
    int in_time_id;
    SUCCESS_OR_DIE(nc_inq_varid(in_ncid, "time", &in_time_id));
    SUCCESS_OR_DIE(nc_inq_varnatts(in_ncid, in_time_id, &n_atts));
    for (int i=0; i < n_atts; i++) {
        char att_name[NC_MAX_NAME];
        SUCCESS_OR_DIE(nc_inq_attname(in_ncid, in_time_id, i, att_name));
        if (!std::strcmp(att_name, "_FillValue")) {
            continue;
        }
        size_t att_len;
        SUCCESS_OR_DIE(nc_inq_att(
                in_ncid, in_time_id, att_name, NULL, &att_len));
        std::vector<char> att_val(att_len);
        SUCCESS_OR_DIE(nc_get_att(
                in_ncid, in_time_id, att_name, att_val.data()));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time],
                att_name,
                att_len,
                att_val.data()));
    }
    SUCCESS_OR_DIE(ncclose(in_ncid));

    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::out_param],
            "long_name",
            strlen("gradients are calculated w.r.t. this output parameter"),
            "gradients are calculated w.r.t. this output parameter"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::out_param],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    if (this->simulation_mode != create_train_set) {
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::out_param],
                "standard_name",
                strlen("output_parameter_id"),
                "output_parameter_id"));
    } else {
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::out_param],
                "standard_name",
                strlen("output_parameter"),
                "output_parameter"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::perturbed],
                "long_name",
                strlen("Name of the perturbed parameter (each ensemble can have different and " \
                       "multiple perturbed parameters)."),
                "Name of the perturbed parameter (each ensemble can have different and " \
                        "multiple perturbed parameters)."));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::perturbed],
                "standard_name",
                strlen("perturbed_parameter"),
                "perturbed_parameter"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::perturbed],
                "auxiliary_data",
                strlen("yes"),
                "yes"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::perturbation_value],
                "long_name",
                strlen("Value of the perturbation."),
                "Value of the perturbation."));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::perturbation_value],
                "standard_name",
                strlen("perturbation_value"),
                "perturbation_value"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::perturbation_value],
                "auxiliary_data",
                strlen("yes"),
                "yes"));
    }

    auto put_att_mass = [&](
            const char *mass_name,
            const char *long_mass_name,
            auto &varid) {
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "units",
                strlen("kg m^-3"),
                "kg m^-3"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "auxiliary_data",
                strlen("yes"),
                "yes"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "long_name",
                strlen(long_mass_name),
                long_mass_name));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "standard_name",
                strlen(mass_name),
                mass_name));
        SUCCESS_OR_DIE(nc_put_att(
                ncid,
                varid,
                _FillValue,
               NC_FLOAT_T,
               1,
               &FILLVALUE));
    };
    put_att_mass("water_vapor_mass_density", "water vapor mass density", varid[Var_idx::qv]);
    put_att_mass("cloud_droplet_mass_density", "cloud droplet mass density", varid[Var_idx::qc]);
    put_att_mass("rain_droplet_mass_density", "rain droplet mass density", varid[Var_idx::qr]);
    put_att_mass("snow_mass_density", "snow mass density", varid[Var_idx::qs]);
    put_att_mass("ice_mass_density", "ice mass density", varid[Var_idx::qi]);
    put_att_mass("graupel_mass_density", "graupel mass density", varid[Var_idx::qg]);
    put_att_mass("hail_mass_density", "hail mass density", varid[Var_idx::qh]);

    auto put_att_nums = [&](
            const char *name,
            const char *long_name,
            auto &varid) {
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "units",
                strlen("m^-3"),
                "m^-3"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "auxiliary_data",
                strlen("yes"),
                "yes"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "long_name",
                strlen(long_name),
                long_name));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "standard_name",
                strlen(name),
                name));
        SUCCESS_OR_DIE(nc_put_att(
                ncid,
                varid,
                _FillValue,
               NC_FLOAT_T,
               1,
               &FILLVALUE));
    };
    put_att_nums("cloud_droplet_number_density", "cloud droplet number density", varid[Var_idx::ncloud]);
    put_att_nums("rain_droplet_number_density", "rain droplet number density", varid[Var_idx::nrain]);
    put_att_nums("snow_number_density", "snow number density", varid[Var_idx::nsnow]);
    put_att_nums("ice_number_density", "ice number density", varid[Var_idx::nice]);
    put_att_nums("graupel_number_density", "graupel number density", varid[Var_idx::ngraupel]);
    put_att_nums("hail_number_density", "hail number density", varid[Var_idx::nhail]);

    auto put_att_mass_sed = [&](
            const char *mass_name,
            const char *long_mass_name,
            auto &varid) {
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "units",
                strlen("kg m^-3 s^-1"),
                "kg m^-3 s^-1"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "auxiliary_data",
                strlen("yes"),
                "yes"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "long_name",
                strlen(long_mass_name),
                long_mass_name));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "standard_name",
                strlen(mass_name),
                mass_name));
        SUCCESS_OR_DIE(nc_put_att(
                ncid,
                varid,
                _FillValue,
               NC_FLOAT_T,
               1,
               &FILLVALUE));
    };
    put_att_mass_sed(
            "sedi_outflux_of_rain_droplet_mass",
            "sedimentation of rain droplet mass",
            varid[Var_idx::qr_out]);
    put_att_mass_sed(
            "sedi_outflux_of_snow_mass",
            "sedimentation of snow mass",
            varid[Var_idx::qs_out]);
    put_att_mass_sed(
            "sedi_outflux_of_ice_mass",
            "sedimentation of ice mass",
            varid[Var_idx::qi_out]);
    put_att_mass_sed(
            "sedi_outflux_of_graupel_mass",
            "sedimentation of graupel mass",
            varid[Var_idx::qg_out]);
    put_att_mass_sed(
            "sedi_outflux_of_hail_mass",
            "sedimentation of hail mass",
            varid[Var_idx::qh_out]);

    auto put_att_nums_sed = [&](
            const char *name,
            const char *long_name,
            auto &varid) {
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "units",
                strlen("m^-3 s^-1"),
                "m^-3 s^-1"));
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "auxiliary_data",
                strlen("yes"),
                "yes"));
        std::string tmp_string = "sedimentation of " + std::string(long_name) + " number";
        const char *att_val = tmp_string.c_str();
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "long_name",
                strlen(att_val),
                att_val));
        std::string tmp_string_2 = "sedi_outflux_of_" + std::string(name) + "_number";
        const char *att_val_2 = tmp_string_2.c_str();
        SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "standard_name",
                strlen(att_val_2),
                att_val_2));
        SUCCESS_OR_DIE(nc_put_att(
                ncid,
                varid,
                _FillValue,
                NC_FLOAT_T,
               1,
               &FILLVALUE));
    };
    put_att_nums_sed("sedi_outflux_of_rain_droplet_number", "sedimentation of rain droplet number",
                     varid[Var_idx::nr_out]);
    put_att_nums_sed("sedi_outflux_of_snow_number", "sedimentation of snow number", varid[Var_idx::ns_out]);
    put_att_nums_sed("sedi_outflux_of_ice_number", "sedimentation of ice number", varid[Var_idx::ni_out]);
    put_att_nums_sed("sedi_outflux_of_graupel_number", "sedimentation of graupel number", varid[Var_idx::ng_out]);
    put_att_nums_sed("sedi_outflux_of_hail_number", "sedimentation of hail number", varid[Var_idx::nh_out]);

    if (!track_ic) {
        // all gradients are auxiliary data
        for (int i=0; i < num_par-num_par_init; i++) {
            if (cc.trace_check(i, false)) {
                SUCCESS_OR_DIE(nc_put_att_text(
                        ncid,
                        varid[Var_idx::n_vars + i],
                        "auxiliary_data",
                        strlen("yes"),
                        "yes"));
                SUCCESS_OR_DIE(nc_put_att(
                        ncid,
                        varid[Var_idx::n_vars + i],
                        _FillValue,
                       NC_FLOAT_T,
                       1,
                       &FILLVALUE));
            }
        }
    } else {
        // initial condition gradients
        for (int i=0; i < num_par_init; i++) {
            if (cc.trace_check(i, 2)) {
                SUCCESS_OR_DIE(nc_put_att_text(
                        ncid,
                        varid[Var_idx::n_vars + i],
                        "auxiliary_data",
                        strlen("yes"),
                        "yes"));
                SUCCESS_OR_DIE(nc_put_att(
                        ncid,
                        varid[Var_idx::n_vars + i],
                        _FillValue,
                       NC_FLOAT_T,
                       1,
                       &FILLVALUE));
            }
        }
    }
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "long_name",
            strlen("pressure"),
            "pressure"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "standard_name",
            strlen("air_pressure"),
            "air_pressure"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "units",
            strlen("Pa"),
            "Pa"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "positive",
            strlen("down"),
            "down"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "axis",
            strlen("Z"),
            "Z"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::pressure],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::temperature],
            "long_name",
            strlen("temperature"),
            "temperature"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::temperature],
            "standard_name",
            strlen("air_temperature"),
            "air_temperature"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::temperature],
            "units",
            strlen("K"),
            "K"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::temperature],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::temperature],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::ascent],
            "long_name",
            strlen("ascend velocity"),
            "ascend velocity"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::ascent],
            "standard_name",
            strlen("ascend_velocity"),
            "ascend_velocity"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::ascent],
            "units",
            strlen("m s^-1"),
            "m s^-1"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::ascent],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::ascent],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sat],
            "long_name",
            strlen("saturation"),
            "saturation"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sat],
            "standard_name",
            strlen("saturation"),
            "saturation"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sat],
            "units",
            strlen("1"),
            "percentage"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sat],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::sat],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::height],
            "long_name",
            strlen("height above mean sea level"),
            "height above mean sea level"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::height],
            "standard_name",
            strlen("height"),
            "height"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::height],
            "units",
            strlen("m AMSL"),
            "m AMSL"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::height],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::height],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::inactive],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::inactive],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::dep],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::dep],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sub],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::sub],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat_heat],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::lat_heat],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat_cool],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::lat_cool],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::time_ascent],
            "long_name",
            strlen("time after rapid ascent started"),
            "time after rapid ascent started"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::time_ascent],
            "standard_name",
            strlen("time_after_ascent"),
            "time_after_ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::time_ascent],
            "units",
            strlen("seconds since start of convective/slantwise ascent"),
            "seconds since start of convective/slantwise ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::time_ascent],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::time_ascent],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat],
            "long_name",
            strlen("rotated latitude"),
            "rotated latitude"));

    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat],
            "standard_name",
            strlen("latitude"),
            "latitude"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat],
            "units",
            strlen("degrees"),
            "degrees"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::lat],
            _FillValue,
            NC_FLOAT_T,
           1,
           &FILLVALUE));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lon],
            "long_name",
            strlen("rotated longitude"),
            "rotated longitude"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lon],
            "standard_name",
            strlen("longitude"),
            "longitude"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lon],
            "units",
            strlen("degrees"),
            "degrees"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::lon],
            _FillValue,
            NC_FLOAT_T,
            1,
            &FILLVALUE));
#if !defined B_EIGHT
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_400],
            "long_name",
            strlen("convective 400hPa ascent"),
            "convective 400hPa ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_400],
            "standard_name",
            strlen("convective_400hPa_ascent"),
            "convective_400hPa_ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_400],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_600],
            "long_name",
            strlen("convective 600hPa ascent"),
            "convective 600hPa ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_600],
            "standard_name",
            strlen("convective_600hPa_ascent"),
            "convective_600hPa_ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_600],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_400],
            "long_name",
            strlen("slantwise 400hPa ascent"),
            "slantwise 400hPa ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_400],
            "standard_name",
            strlen("slantwise_400hPa_ascent"),
            "slantwise_400hPa_ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_400],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_600],
            "long_name",
            strlen("slantwise 600hPa ascent"),
            "slantwise 600hPa ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_600],
            "standard_name",
            strlen("slantwise_600hPa_ascent"),
            "slantwise_600hPa_ascent"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_600],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
#endif
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::step],
            "long_name",
            strlen("simulation step"),
            "simulation step"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::step],
            "standard_name",
            strlen("step"),
            "step"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::step],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    const uint64_t FILLINT = 0;
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::step],
            _FillValue,
            NC_UINT64,
            1,
            &FILLINT));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::phase],
            "long_name",
            strlen("phase"),
            "phase"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::phase],
            "standard_name",
            strlen("phase"),
            "phase"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::phase],
            "description",
            strlen("0: warm phase, 1: mixed phase, 2: ice phase, 3: neutral phase"),
            "0: warm phase, 1: mixed phase, 2: ice phase, 3: neutral phase"));
    SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::phase],
            "auxiliary_data",
            strlen("yes"),
            "yes"));
    SUCCESS_OR_DIE(nc_put_att(
            ncid,
            varid[Var_idx::phase],
            _FillValue,
            NC_UINT64,
            1,
            &FILLINT));
}


template<class float_t>
void output_handle_t::set_compression(const model_constants_t<float_t> &cc) {
#ifdef DEVELOP
    std::cout << "Using compression\n";
#endif
    if (this->simulation_mode == limited_time_ensembles
        || this->simulation_mode == trajectory_sensitivity
        || this->simulation_mode == grid_sensitivity
        || this->simulation_mode == create_train_set) {
        // Compressing this is buggy
        for (uint32_t i=0; i < Var_idx::n_vars; ++i) {
            // This does not work on scalars; we need to filter those out
            if ((local_num_comp > 1 || i != static_cast<int>(Var_idx::out_param)) &&
              (n_perturbed_params > 1 || i != static_cast<int>(Var_idx::perturbed)))
                if (simulation_mode == create_train_set || i != static_cast<int>(Var_idx::perturbation_value) )
                    // This could be a version for szip
    //                         SUCCESS_OR_DIE(
    //                             nc_def_var_szip(
    //                                 ncid,
    //                                 varid[i],
    //                                 NC_SZIP_NN,
    //                                 8));  // pixels per block
                    // zlib version
                    SUCCESS_OR_DIE(
                            nc_def_var_deflate(
                                    ncid,
                                    varid[i],
                                    1,  // shuffle
                                    1,  // deflate
                                    9));  // max compression
        }
        if (!track_ic) {
            // gradients
            for (uint32_t i=0; i < num_par-num_par_init; ++i) {
                if (cc.trace_check(i, false))
//                         SUCCESS_OR_DIE(
//                             nc_def_var_szip(
//                                 ncid,
//                                 varid[Var_idx::n_vars + i],
//                                 NC_SZIP_NN,
//                                 8));
                    SUCCESS_OR_DIE(
                            nc_def_var_deflate(
                                    ncid,
                                    varid[Var_idx::n_vars + i],
                                    1,  // shuffle
                                    1,  // deflate
                                    9));  // max compression
            }
        } else {
            // initial conditions
            for (uint32_t i=0; i < num_par_init; ++i) {
                if (cc.trace_check(i, 2))
//                         SUCCESS_OR_DIE(
//                             nc_def_var_szip(
//                                 ncid,
//                                 varid[Var_idx::n_vars + i],
//                                 NC_SZIP_NN,
//                                 8));
                    SUCCESS_OR_DIE(
                            nc_def_var_deflate(
                                    ncid,
                                    varid[Var_idx::n_vars + i],
                                    1,  // shuffle
                                    1,  // deflate
                                    9));  // max compression
            }
        }
    }
#ifdef DEVELOP
    std::cout << "all done; attempt to close\n" << std::flush;
#endif
}


template<class float_t>
void output_handle_t::write_dimension_values(
    const model_constants_t<float_t> &cc,
    const double delay_out_time) {

    std::vector<size_t> startp, countp;
    startp.push_back(0);
    countp.push_back(num_time);
#ifdef OUT_DOUBLE
    std::vector<double> time_steps(num_time);
#else
    std::vector<float> time_steps(num_time);
#endif
    for (uint32_t i=0; i < num_time; i++)
        time_steps[i] = cc.dt*i + cc.start_time + delay_out_time;
#ifdef DEVELOP
    std::cout << "attempt to write\n" << std::flush;
#endif
    SUCCESS_OR_DIE(
            nc_put_vara(
                    ncid,                               // ncid
                    varid[Var_idx::time],               // varid
                    startp.data(),                      // startp
                    countp.data(),                      // countp
                    time_steps.data()));                // op
#ifdef DEVELOP
    std::cout << "done time_steps\n" << std::flush;
#endif
    countp[0] = n_trajs_file;
    std::vector<uint64_t> data(n_trajs_file);
    for (uint32_t i=0; i < n_trajs_file; i++)
        data[i] = i;
    SUCCESS_OR_DIE(
            nc_put_vara(
                    ncid,
                    varid[Var_idx::trajectory],
                    startp.data(),
                    countp.data(),
                    data.data()));
#ifdef DEVELOP
    std::cout << "done n_trajs; starting with num_ens " << num_ens << "\n" << std::flush;
#endif
    countp[0] = num_ens;
    data.resize(num_ens);
    for (uint32_t i=0; i < num_ens; i++)
        data[i] = i;
    SUCCESS_OR_DIE(
            nc_put_vara(
                    ncid,
                    varid[Var_idx::ensemble],
                    startp.data(),
                    countp.data(),
                    data.data()));
#ifdef DEVELOP
    std::cout << "out_param\n" << std::flush;
#endif
    countp[0] = local_num_comp;
    if (this->simulation_mode != create_train_set) {
        data.resize(local_num_comp);
        uint32_t counter = 0;
        for (uint32_t i = 0; i < num_comp; i++) {
            if (cc.trace_check(i, true)) {
                data[counter] = i;
                counter++;
            }
        }
        SUCCESS_OR_DIE(
                nc_put_vara(
                        ncid,
                        varid[Var_idx::out_param],
                        startp.data(),
                        countp.data(),
                        data.data()));
    }
#ifdef DEVELOP
    std::cout << "out_param done\n" << std::flush;
#endif
#ifdef DEVELOP
    std::cout << "write perturbed_id\n" << std::flush;
#endif
    if (this->simulation_mode == create_train_set) {
        std::vector<const char*> cstrings_outp;
        cstrings_outp.reserve(local_num_comp);
        for (uint32_t i = 0; i < num_comp; i++) {
            if (cc.trace_check(i, true)) {
                cstrings_outp.push_back(&(output_par_idx[i])[0]);
            }
        }
        for (auto s : cstrings_outp) {
            SUCCESS_OR_DIE(
                    nc_put_var1_string(
                            ncid,
                            varid[Var_idx::out_param],
                            startp.data(),
                            &s));
            startp[0]++;
        }
        startp[0] = 0;
        countp[0] = n_perturbed_params;
        std::vector<const char*> cstrings;
        cstrings.reserve(perturbed_names.size());

        for (const auto& s : perturbed_names)
            cstrings.push_back(&s[0]);

        for (auto s : cstrings) {
            SUCCESS_OR_DIE(
                nc_put_var1_string(
                    ncid,
                    varid[Var_idx::perturbed],
                    startp.data(),
                    &s));
            startp[0]++;
        }
    }

#ifdef DEVELOP
    std::cout << "perturbed_id done\n" << std::flush;
#endif
    SUCCESS_OR_DIE(ncclose(ncid));
#ifdef DEVELOP
    std::cout << "closed\n" << std::flush;
#endif
}


template<class float_t>
void output_handle_t::set_parallel_access(
    const model_constants_t<float_t> &cc,
    const std::string file_string) {
#ifdef DEVELOP
    std::cout << "attempt to open " << filename << "\n" << std::flush;
#endif
    SUCCESS_OR_DIE(
        nc_open_par(
            file_string.c_str(),
            NC_WRITE,
            MPI_COMM_WORLD,
            MPI_INFO_NULL,
            &ncid));
#ifdef DEVELOP
    std::cout << "get varids\n" << std::flush;
#endif
    // gather all necessary variable ids
    if (this->simulation_mode == create_train_set) {
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "Output_Parameter",
                &varid[Var_idx::out_param]));
    } else {
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "Output_Parameter_ID",
                &varid[Var_idx::out_param]));
    }
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "ensemble",
            &varid[Var_idx::ensemble]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "trajectory",
            &varid[Var_idx::trajectory]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "time",
            &varid[Var_idx::time]));
    if (this->simulation_mode == create_train_set) {
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "perturbed",
                &varid[Var_idx::perturbed]));
    }
#ifdef DEVELOP
    std::cout << "get model varids\n" << std::flush;
#endif
    // model state
    for (uint32_t i=0; i < num_comp; ++i)
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                output_par_idx[i].c_str(),
                &varid[i]));
#ifdef DEVELOP
    std::cout << "get track varids\n" << std::flush;
#endif
    if (!track_ic) {
        // gradients
        for (uint32_t i=0; i < num_par-num_par_init; ++i) {
            if (cc.trace_check(i, false)) {
                SUCCESS_OR_DIE(
                    nc_inq_varid(
                        ncid,
                        output_grad_idx[i].c_str(),
                        &varid[Var_idx::n_vars + i]));
            }
        }
    } else {
        // initial conditions
        for (uint32_t i=0; i < num_par_init; ++i) {
            if (cc.trace_check(i, 2)) {
                SUCCESS_OR_DIE(
                    nc_inq_varid(
                        ncid,
                        init_grad_idx[i].c_str(),
                        &varid[Var_idx::n_vars + i]));
            }
        }
    }
#ifdef DEVELOP
    std::cout << "get more time\n" << std::flush;
#endif
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "time_after_ascent",
            &varid[Var_idx::time_ascent]));
#if !defined(B_EIGHT)
#ifdef DEVELOP
    std::cout << "get conv\n" << std::flush;
#endif
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "conv_400",
            &varid[Var_idx::conv_400]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "conv_600",
            &varid[Var_idx::conv_600]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "slan_400",
            &varid[Var_idx::slan_400]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "slan_600",
            &varid[Var_idx::slan_600]));
#endif
#ifdef DEVELOP
    std::cout << "get lat\n" << std::flush;
#endif
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "lat",
            &varid[Var_idx::lat]));
#ifdef DEVELOP
    std::cout << "get lon\n" << std::flush;
#endif
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "lon",
            &varid[Var_idx::lon]));
#ifdef DEVELOP
    std::cout << "get step\n" << std::flush;
#endif
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "step",
            &varid[Var_idx::step]));
#ifdef DEVELOP
    std::cout << "get phase\n" << std::flush;
#endif
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "phase",
            &varid[Var_idx::phase]));

    if (simulation_mode == create_train_set) {
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "perturbation_value",
                &varid[Var_idx::perturbation_value]));
    }

    if ((this->simulation_mode == trajectory_sensitvity_perturbance)
        || (this->simulation_mode == trajectory_perturbance)) {
#ifdef DEVELOP
        std::cout << "different accesses\n" << std::flush;
#endif
        // Make the access independent which is a must due to the dynamic
        // work schedule; This can be expensive though.
        for (uint32_t i=0; i < Var_idx::n_vars; i++)
            SUCCESS_OR_DIE(nc_var_par_access(ncid, varid[i], NC_INDEPENDENT));
        if (!track_ic) {
            for (uint32_t i=0; i < num_par-num_par_init; i++)
                if (cc.trace_check(i, false))
                    SUCCESS_OR_DIE(nc_var_par_access(ncid, varid[i+Var_idx::n_vars], NC_INDEPENDENT));
        } else {
            for (uint32_t i=0; i < num_par_init; i++)
                if (cc.trace_check(i, 2))
                    SUCCESS_OR_DIE(nc_var_par_access(ncid, varid[i+Var_idx::n_vars], NC_INDEPENDENT));
        }
    }
#ifdef COMPRESS_OUTPUT
    // We need to explicitly tell that we want to collectively write the data
    if ((this->simulation_mode != trajectory_sensitvity_perturbance)
        && (this->simulation_mode != trajectory_perturbance)) {
#ifdef DEVELOP
        std::cout << "different accesses\n" << std::flush;
#endif
        // Make the access independent which is a must due to the dynamic
        // work schedule; This can be expensive though.
        for (uint32_t i=0; i < Var_idx::n_vars; i++) {
            SUCCESS_OR_DIE(nc_var_par_access(ncid, varid[i], NC_COLLECTIVE));
        }
        if (!track_ic) {
            for (uint32_t i=0; i < num_par-num_par_init; i++)
                if (cc.trace_check(i, false))
                    SUCCESS_OR_DIE(nc_var_par_access(ncid, varid[i+Var_idx::n_vars], NC_COLLECTIVE));
        } else {
            for (uint32_t i=0; i < num_par_init; i++)
                if (cc.trace_check(i, 2))
                    SUCCESS_OR_DIE(nc_var_par_access(ncid, varid[i+Var_idx::n_vars], NC_COLLECTIVE));
        }
    }
#endif
}

template<class float_t>
void output_handle_t::setup(
    const std::string filetype,
    const std::string filename,
    const model_constants_t<float_t> &cc,
    const std::string in_filename,
    const uint32_t write_index,
    const uint32_t snapshot_index,
    const int &rank,
    const double delay_out_time) {

    this->n_trajs_file = cc.max_n_trajs;
    this->traj = cc.traj_id;
    this->ens = cc.ensemble_id;
    this->num_ens = cc.n_ensembles;
    this->num_time = cc.num_steps * cc.num_sub_steps;
    if (delay_out_time > 0) {
        this->num_time -= delay_out_time / (cc.dt_prime * cc.num_sub_steps)-cc.num_sub_steps;
    }

    this->filetype = filetype;
    this->filename = filename;
    dimid.resize(Dim_idx::n_dims);
    if (track_ic) {
        varid.resize(Var_idx::n_vars + num_par_init);
    } else {
        varid.resize(Var_idx::n_vars + num_par-num_par_init);
    }


    const std::string ending = ".nc_wcb";
    const std::string ending2 = ".nc";
    if (!std::equal(ending.rbegin(), ending.rend(), this->filename.rbegin())) {
        if (!std::equal(ending2.rbegin(), ending2.rend(), this->filename.rbegin())) {
            this->filename += ".nc";
        }
    }
    if (rank == 0)
        std::cout << "Creating " << this->filename << " and defining dimensions "
            << "and attributes. This can take a while for some simulation "
            << "modes\n";

    flushed_snapshots = 0;
    n_snapshots = 0;
    // Allocate memory for the buffer
    // maximum number of snapshots we are going to get
    total_snapshots = std::ceil((static_cast<float>(write_index))/snapshot_index);
    const uint64_t vec_size = total_snapshots;
    const uint64_t vec_size_grad = local_num_comp * total_snapshots;
    // model state
    for (uint32_t i=0; i < num_comp; i++)
        output_buffer[i].resize(vec_size);
    output_buffer[Buffer_idx::time_ascent_buf].resize(vec_size);
    output_buffer[Buffer_idx::lat_buf].resize(vec_size);
    output_buffer[Buffer_idx::lon_buf].resize(vec_size);

    if (this->simulation_mode == create_train_set)
        output_buffer[Buffer_idx::perturb_buf].resize(vec_size*n_perturbed_params);
    if (!track_ic) {
        // gradients
        for (uint32_t i=Buffer_idx::n_buffer; i < Buffer_idx::n_buffer+num_par-num_par_init; i++)
            if (cc.trace_check(i-Buffer_idx::n_buffer, false))
                output_buffer[i].resize(vec_size_grad);
    } else {
        // the initial condition sensitivities are stored here
        for (uint32_t i=Buffer_idx::n_buffer;
            i < Buffer_idx::n_buffer+num_par_init; i++)
            if (cc.trace_check(i-Buffer_idx::n_buffer, 2))
                output_buffer[i].resize(vec_size_grad);
    }

    for (uint32_t i=0; i < output_buffer_flags.size(); i++)
        output_buffer_flags[i].resize(vec_size);

    for (uint32_t i=0; i < output_buffer_int.size(); i++)
        output_buffer_int[i].resize(vec_size);

    // Unfortunately, it is likely to run into an HDF Error if
    // many parameters are investigated in combination with the
    // amount of attributes for the time dimension.
    // We have to create the file as a single process, close it and
    // open it in parallel again for writing purpose
    if (rank == 0) {
#ifdef DEVELOP
        std::cout << "Using local_num_comp " << local_num_comp << "\n"
            << "num_ens " << num_ens << "\n"
            << "n_trajs_file " << n_trajs_file << "\n"
            << "num_time " << num_time << "\n" << std::flush;
#endif
        SUCCESS_OR_DIE(nc_create(
            filename.c_str(),
            NC_NETCDF4,
            &ncid));
        // Create dimensions
        // If there is only one model state variable for which we gather
        // sensitivities for, drop that dimension and add that as a scalar.
        // This is useful as long as Met3D does not support this dimension.
        if (local_num_comp > 1) {
            if (this->simulation_mode == create_train_set) {
                SUCCESS_OR_DIE(nc_def_dim(
                    ncid, "Output_Parameter", local_num_comp, &dimid[Dim_idx::out_param_dim]));
            } else {
                SUCCESS_OR_DIE(nc_def_dim(
                    ncid, "Output_Parameter_ID", local_num_comp, &dimid[Dim_idx::out_param_dim]));
            }
        }
        SUCCESS_OR_DIE(nc_def_dim(
                ncid, "trajectory", n_trajs_file, &dimid[Dim_idx::trajectory_dim]));
        SUCCESS_OR_DIE(nc_def_dim(
                ncid, "ensemble", num_ens, &dimid[Dim_idx::ensemble_dim]));
        SUCCESS_OR_DIE(nc_def_dim(
                ncid, "time", num_time, &dimid[Dim_idx::time_dim]));
        // The training set stores the perturbed parameter for each ensemble as additional
        // information
        if (simulation_mode == create_train_set) {
            SUCCESS_OR_DIE(nc_def_dim(
                    ncid, "perturbed", n_perturbed_params, &dimid[Dim_idx::perturb_param_dim]));
        }
#ifdef DEBUG_SEG
        std::cout << "trajs " << n_trajs_file << " ens " << num_ens << " time " << num_time << "\n";
#endif
        define_vars(cc);
        set_attributes(cc, in_filename);
#ifdef COMPRESS_OUTPUT
        set_compression(cc);
#endif
        SUCCESS_OR_DIE(nc_enddef(ncid));
#ifdef DEVELOP
        std::cout << "closed\n" << std::flush;
#endif
        write_dimension_values(cc, delay_out_time);
    }
    // Open it again for writing with all processes
    set_parallel_access(cc, filename);
#ifdef DEVELOP
    std::cout << "done\n" << std::flush;
#endif
}


void output_handle_t::reset(
    const uint32_t traj_id,
    const uint32_t ens_id,
    const uint64_t n_flushed) {

    this->flushed_snapshots = n_flushed;
    this->n_snapshots = 0;
    this->traj = traj_id;
    this->ens = ens_id;
}


template<class float_t>
void output_handle_t::buffer_gradient(
    const model_constants_t<float_t> &cc,
    const std::vector< std::array<double, num_par > >  &y_diff,
    const uint32_t snapshot_index) {

    uint64_t comp_idx = 0;
    for (uint64_t i=0; i < num_comp; i++) {
        // gradient sensitive to output parameter i
        if (!cc.trace_check(i, true))
            continue;
        if (!track_ic) {
            for (uint64_t j=0; j < num_par-num_par_init; j++)  // gradient of input parameter j
                if (cc.trace_check(j, false)) {
                    if (n_snapshots%snapshot_index == 0) {
                        output_buffer[Buffer_idx::n_buffer+j][comp_idx*total_snapshots + n_snapshots] =
                            y_diff[i][j]/snapshot_index;
                    } else {
                        output_buffer[Buffer_idx::n_buffer+j][comp_idx*total_snapshots + n_snapshots] +=
                            y_diff[i][j]/snapshot_index;
                    }
                }
        } else {
            for (uint64_t j=0; j < num_par_init; j++)  // gradient of initial condition j
                if (cc.trace_check(j, 2)) {
                    if (n_snapshots%snapshot_index == 0) {
                        output_buffer[Buffer_idx::n_buffer+j][comp_idx*total_snapshots + n_snapshots] =
                            y_diff[i][j]/snapshot_index;
                    } else {
                        output_buffer[Buffer_idx::n_buffer+j][comp_idx*total_snapshots + n_snapshots] +=
                            y_diff[i][j]/snapshot_index;
                    }
                }
        }
        comp_idx++;
    }
}


template<class float_t>
void output_handle_t::buffer(
    const model_constants_t<float_t> &cc,
    const netcdf_reader_t &netcdf_reader,
    const std::vector<float_t> &y_single_new,
    const std::vector< std::array<double, num_par > >  &y_diff,
    const uint32_t sub,
    const uint32_t t,
    const reference_quantities_t &ref_quant,
    const uint32_t snapshot_index) {

    // output parameters
    for (uint64_t i=0; i < num_comp; i++) {
        switch (i) {
            case p_idx:
                output_buffer[i][n_snapshots] =
                    y_single_new[i].getValue() * ref_quant.pref;
                break;
            case T_idx:
                output_buffer[i][n_snapshots] =
                    y_single_new[i].getValue() * ref_quant.Tref;
                break;
            case w_idx:
                output_buffer[i][n_snapshots] =
                    y_single_new[i].getValue() * ref_quant.wref;
                break;
            case z_idx:
                output_buffer[i][n_snapshots] =
                    y_single_new[i].getValue() * ref_quant.zref;
                break;
            case qc_idx:
            case qr_idx:
            case qv_idx:
            case qi_idx:
            case qs_idx:
            case qg_idx:
            case qh_idx:
            case qi_out_idx:
            case qs_out_idx:
            case qr_out_idx:
            case qg_out_idx:
            case qh_out_idx:
                output_buffer[i][n_snapshots] =
                    y_single_new[i].getValue() * ref_quant.qref;
                break;
            case Nc_idx:
            case Nr_idx:
            case Ni_idx:
            case Ns_idx:
            case Ng_idx:
            case Nh_idx:
            case Ni_out_idx:
            case Ns_out_idx:
            case Nr_out_idx:
            case Ng_out_idx:
            case Nh_out_idx:
                output_buffer[i][n_snapshots] =
                    y_single_new[i].getValue() * ref_quant.Nref;
                break;
            default:
                output_buffer[i][n_snapshots] =
                    y_single_new[i].getValue();
                break;
        }
    }
    if ( (this->simulation_mode == create_train_set && cc.ensemble_id == 0)
        || (this->simulation_mode == limited_time_ensembles && cc.traj_id == 0)
        || (this->simulation_mode != create_train_set && this->simulation_mode != limited_time_ensembles) )
        buffer_gradient(cc, y_diff, snapshot_index);
    // lat
    output_buffer[Buffer_idx::lat_buf][n_snapshots] = netcdf_reader.get_lat(t, sub);

    // lon
    output_buffer[Buffer_idx::lon_buf][n_snapshots] = netcdf_reader.get_lon(t, sub);

#ifdef MET3D
    // time after ascent
    output_buffer[Buffer_idx::time_ascent_buf][n_snapshots] =
        netcdf_reader.get_relative_time(t) + sub*cc.dt;
    // flags
#if !defined(B_EIGHT)
    output_buffer_flags[0][n_snapshots] = netcdf_reader.get_conv_400(t);
    output_buffer_flags[1][n_snapshots] = netcdf_reader.get_conv_600(t);
    output_buffer_flags[2][n_snapshots] = netcdf_reader.get_slan_400(t);
    output_buffer_flags[3][n_snapshots] = netcdf_reader.get_slan_600(t);
#endif
#endif
    // Value of perturbed parameter
    if (this->simulation_mode == create_train_set) {
        std::vector<float> perturbed;
        std::vector<uint64_t> param_idx;
        cc.get_perturbed_info(perturbed, param_idx);
        // This stores the difference to the unperturbed parameters
        std::fill(output_buffer[Buffer_idx::perturb_buf].begin() + n_snapshots*n_perturbed_params,
                  output_buffer[Buffer_idx::perturb_buf].begin() + (n_snapshots+1)*n_perturbed_params, 0);
        for (uint32_t i=0; i < perturbed.size(); ++i) {
            const uint64_t idx = param_idx[i];
            output_buffer[Buffer_idx::perturb_buf][n_snapshots * n_perturbed_params + perturbed_idx[idx]] =
                    perturbed[i] - unperturbed_vals[perturbed_idx[idx]];
        }
        // The following stores the new values
//        std::copy(output_buffer[Buffer_idx::perturb_buf].begin() + n_snapshots*n_perturbed_params,
//                  output_buffer[Buffer_idx::perturb_buf].begin() + (n_snapshots+1)*n_perturbed_params,
//                  unperturbed_vals);
//        for (uint32_t i=0; i < perturbed.size(); ++i) {
//            const uint64_t idx = param_idx[i];
//            output_buffer[Buffer_idx::perturb_buf][n_snapshots * n_perturbed_params + perturbed_idx[idx]] =
//            perturbed[i];
//        }
    }

    // simulation step
    output_buffer_int[0][n_snapshots] = sub + t*cc.num_sub_steps;

    // phase, 0: warm phase, 1: mixed phase, 2: ice phase,
    // 3: only water vapor or nothing at all
    uint64_t current_phase = 3;
    if (y_single_new[qi_idx] > 0 || y_single_new[qs_idx] > 0
        || y_single_new[qh_idx] > 0 || y_single_new[qg_idx] > 0
        || y_single_new[Ni_idx] > 0 || y_single_new[Ns_idx] > 0
        || y_single_new[Nh_idx] > 0 || y_single_new[Ng_idx] > 0) {
        current_phase = 2;
    }
    if (y_single_new[qc_idx] > 0 || y_single_new[qr_idx] > 0
        || y_single_new[Nc_idx] > 0 || y_single_new[Nr_idx] > 0) {
        current_phase = (current_phase == 2) ? 1 : 0;
    }
    output_buffer_int[1][n_snapshots] = current_phase;

    n_snapshots++;
}


template<class float_t>
bool output_handle_t::flush_buffer(
    const model_constants_t<float_t> &cc,
    bool no_flush) {
#ifdef COMPRESS_OUTPUT
    int needed = (no_flush) ? 0 : 1;
    std::vector<int> needed_receive(n_processes, 0);
    SUCCESS_OR_DIE(
            MPI_Allgather(
                    &needed,
                    1,
                    MPI_INT,
                    needed_receive.data(),
                    1,
                    MPI_INT,
                    MPI_COMM_WORLD));
    if (no_flush) {
        bool return_early = true;
        for (auto const &n : needed_receive) {
            if (n == 1) {
                return_early = false;
                break;
            }
        }
        if (return_early) {
            return false;
        }
    }
#endif
    std::vector<size_t> startp, countp;
    if (no_flush) {
        startp.push_back(0);
        startp.push_back(0);
        startp.push_back(0);
        countp.push_back(0);
        countp.push_back(0);
        countp.push_back(0);
    } else {
        startp.push_back(ens);
        startp.push_back(traj);
        startp.push_back(flushed_snapshots);
        countp.push_back(1);
        countp.push_back(1);
        countp.push_back(n_snapshots);
    }
#ifdef DEVELOP
    std::cout << "traj: " << traj << " at " << flushed_snapshots
              << " with n_snapshots: " << n_snapshots << "\n";
#endif
    for (uint64_t i=0; i < num_comp; i++) {
#if defined(DEVELOP)
//        if (i >= 12 && i < 21) continue;
        std::cout << "traj: " << traj << " at " << flushed_snapshots
                  << " i: " << i << "/" << num_comp
                  << " output_size: " << output_buffer[i].size()
                  << " varid: " << varid[i] << "\n";

#endif
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid[i],
                startp.data(),
                countp.data(),
                output_buffer[i].data()));
    }
#ifdef DEVELOP
    std::cout << "traj: " << traj << " at " << flushed_snapshots
              << " flushed the model state variables" << "\n";
#endif
    // time after ascent
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::time_ascent],
            startp.data(),
            countp.data(),
            output_buffer[Buffer_idx::time_ascent_buf].data()));
#ifdef DEVELOP
    std::cout << "traj: " << traj << " at " << flushed_snapshots
              << " flushed time_ascent" << "\n";
#endif
    // flags
    for (uint64_t i=0; i < output_buffer_flags.size(); i++) {
#if !defined B_EIGHT
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid[Var_idx::conv_400+i],
                startp.data(),
                countp.data(),
                output_buffer_flags[i].data()));
#endif
    }
    // lat
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::lat],
            startp.data(),
            countp.data(),
            output_buffer[Buffer_idx::lat_buf].data()));
    // lon
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::lon],
            startp.data(),
            countp.data(),
            output_buffer[Buffer_idx::lon_buf].data()));
    // step
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::step],
            startp.data(),
            countp.data(),
            output_buffer_int[0].data()));

    // phase
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::phase],
            startp.data(),
            countp.data(),
            output_buffer_int[1].data()));

    // Perturbed parameter value
    if (this->simulation_mode == create_train_set) {
        std::vector<size_t> startp2, countp2;
        if (no_flush) {
            startp2.push_back(0);
            startp2.push_back(0);
            startp2.push_back(0);
            countp2.push_back(0);
            countp2.push_back(0);
            countp2.push_back(0);
        } else {
            startp2.push_back(ens);
            startp2.push_back(flushed_snapshots);
            startp2.push_back(0);
            countp2.push_back(1);
            countp2.push_back(n_snapshots);
            countp2.push_back(n_perturbed_params);
        }
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid[Var_idx::perturbation_value],
                startp2.data(),
                countp2.data(),
                output_buffer[Buffer_idx::perturb_buf].data()));
    }

    // gradients
    if (this->simulation_mode == create_train_set) {
        if (local_num_comp <= 1) {
            // remove ensemble dimension since this does not exist in that
            // simulation mode and it cannot be reused for dimension output_parameter
            startp.erase(startp.begin());
            countp.erase(countp.begin());
        } else {
            startp[0] = 0;
            countp[0] = (no_flush) ? 0 : local_num_comp;
        }
        if (cc.ensemble_id != 0) {
            no_flush = true;
            for (auto &s : startp) s = 0;
            for (auto &c : countp) c = 0;
        }
    }

    if (this->simulation_mode != limited_time_ensembles) {
        if (local_num_comp > 1 && this->simulation_mode != create_train_set) {
            // This is done for create_train_set already
            startp.insert(startp.begin(), 0);
            if (no_flush) {
                countp.insert(countp.begin(), 0);
            } else {
                countp.insert(countp.begin(), local_num_comp);
            }
        }
#ifdef DEVELOP
        std::cout << "traj: " << traj << " at " << flushed_snapshots
                  << " not limited time ensemble; flush gradients\n";
#endif
        // Use an offset if the number of snapshots does not fit
        // This is necessary since the slow index is [0] (Output Parameter)
        // and the fast index is [3] (time), which has gaps now
        if (local_num_comp > 1 && countp[3] != total_snapshots) {
#ifdef COMPRESS_OUTPUT
            // Write strided data for the defined simulation modes
            if (this->simulation_mode == trajectory_sensitivity
                || this->simulation_mode == grid_sensitivity
                || this->simulation_mode == create_train_set) {
                // uint64_t n_snapshots_per_comp = n_snapshots/local_num_comp;
                // For compression: loop over output variables for which
                // we gather sensitivities.
                countp[0] = (no_flush) ? 0 : 1;
                // The dimension are now
                // for create_train_set: (output_parameter_id), trajectory, time
                // all other modes: (output_parameter_id), ensemble, trajectory, time
                uint64_t comp_idx = 0;
                if (!track_ic) {
                    for (int i=0; i < local_num_comp; i++) {
                        // gradient sensitive to output parameter i
                        if (!cc.trace_check(i, true))
                            continue;
                        startp[0] = (no_flush) ? 0 : comp_idx;

                        for (uint64_t j=0; j < num_par-num_par_init; j++) {
                            if (cc.trace_check(j, false)) {
#ifdef DEVELOP
                                std::cout << "traj: " << traj << " at " << flushed_snapshots
                                          << " param " << i << "/" << local_num_comp
                                          << " gradient " << j << "/" << num_par-num_par_init << "\n";
#endif
                                SUCCESS_OR_DIE(
                                    nc_put_vara(
                                        ncid,
                                        varid[Var_idx::n_vars + j],
                                        startp.data(),
                                        countp.data(),
                                        output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()
                                        + comp_idx*total_snapshots));
                            }
                        }
                        comp_idx++;
                    }
                } else {
                    // initial conditions
                    for (int i=0; i < local_num_comp; i++) {
                        // gradient sensitive to output parameter i
                        if (!cc.trace_check(i, true))
                            continue;
                        startp[0] = (no_flush) ? 0 : comp_idx;

                        for (uint64_t j=0; j < num_par_init; j++) {
                            if (cc.trace_check(j, 2)) {
#ifdef DEVELOP
                                std::cout << "traj: " << traj << " at " << flushed_snapshots
                                          << " IC param " << i << "/" << local_num_comp
                                          << " gradient " << j << "/" << num_par-num_par_init << "\n";
#endif
                                SUCCESS_OR_DIE(
                                    nc_put_vara(
                                        ncid,
                                        varid[Var_idx::n_vars + j],
                                        startp.data(),
                                        countp.data(),
                                        output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()
                                        + comp_idx*total_snapshots));
                            }
                        }
                        comp_idx++;
                    }
                }
#ifdef DEVELOP
                std::cout << "traj: " << traj << " at " << flushed_snapshots
                  << " not limited time ensemble; flushed gradients\n";
#endif
            } else {
#endif  // End COMPRESS_OUTPUT
                // In case of compressed output: Write strided data for simulation modes
                // other than trajectory_sensitivity, create_train_set and grid_sensitivity
                // For not compressed output: Strided output is the same for all simulation modes
                // if compression is dsiabled

                // nc_put_varm is discouraged. We use it only if no compression is
                // enabled, since it fails sometimes with compression.
                std::vector<std::ptrdiff_t> stridep, imap;
                if (no_flush) {
                    stridep.push_back(1);
                    stridep.push_back(1);
                    stridep.push_back(1);
                    stridep.push_back(1);
                    imap.push_back(1);
                    imap.push_back(1);
                    imap.push_back(1);
                    imap.push_back(1);
                } else {
                    stridep.push_back(1);
                    stridep.push_back(1);
                    stridep.push_back(1);
                    stridep.push_back(1);
                    imap.push_back(total_snapshots);
                    imap.push_back(1);
                    imap.push_back(1);
                    imap.push_back(1);
                }


                if (!track_ic) {
                    for (uint64_t j=0; j < num_par-num_par_init; j++) {
                        if (cc.trace_check(j, false)) {
#ifdef DEVELOP
                            std::cout << "traj: " << traj << " at " << flushed_snapshots
                                      << " gradient " << j << "/" << num_par-num_par_init << "\n";
#endif
                            SUCCESS_OR_DIE(
                                nc_put_varm(
                                    ncid,
                                    varid[Var_idx::n_vars + j],
                                    startp.data(),
                                    countp.data(),
                                    stridep.data(),
                                    imap.data(),
                                    output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                        }
                    }

                } else {
                    // initial conditions sensitivity
                    for (uint64_t j=0; j < num_par_init; j++) {
                        if (cc.trace_check(j, 2)) {
#ifdef DEVELOP
                            std::cout << "traj: " << traj << " at " << flushed_snapshots
                                      << " IC gradient " << j << "/" << num_par-num_par_init << "\n";
#endif
                            SUCCESS_OR_DIE(
                                nc_put_varm(
                                    ncid,
                                    varid[Var_idx::n_vars + j],
                                    startp.data(),
                                    countp.data(),
                                    stridep.data(),
                                    imap.data(),
                                    output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                        }
                    }
                }
#ifdef COMPRESS_OUTPUT
            }
#endif
        } else {
            // Compressed output: No strided output needed
            // No compression: No strided output needed
            if (!track_ic) {
                for (uint64_t j=0; j < num_par-num_par_init; j++) {
                    if (cc.trace_check(j, false)) {
#ifdef DEVELOP
                        std::cout << "traj: " << traj << " at " << flushed_snapshots
                                  << " Two gradient " << j << "/" << num_par-num_par_init
                                  << " varid: " << varid[Var_idx::n_vars + j] << "\n";
#endif
                        SUCCESS_OR_DIE(
                            nc_put_vara(
                                ncid,
                                varid[Var_idx::n_vars + j],
                                startp.data(),
                                countp.data(),
                                output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                    }
                }
            } else {
                // initial conditions
                for (uint64_t j=0; j < num_par_init; j++) {
                    if (cc.trace_check(j, 2)) {
#ifdef DEVELOP
                        std::cout << "traj: " << traj << " at " << flushed_snapshots
                                  << " IC Two gradient " << j << "/" << num_par-num_par_init << "\n";
#endif
                        SUCCESS_OR_DIE(
                            nc_put_vara(
                                ncid,
                                varid[Var_idx::n_vars + j],
                                startp.data(),
                                countp.data(),
                                output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                    }
                }
            }
        }
#ifdef COMPRESS_OUTPUT
    } else {
        // Compression doesn't work in the simulation mode limited_time_ensemble with
        // strided writes. We have to do that manually.
        std::vector<size_t> startp2, countp2;
        startp2.push_back(0);
        if (no_flush || cc.traj_id != 0) {
            startp2.push_back(0);
            countp2.push_back(0);
            countp2.push_back(0);
        } else {
            startp2.push_back(flushed_snapshots);
            countp2.push_back(1);
            countp2.push_back(n_snapshots);
        }
        if (!track_ic) {
            for (uint64_t j=0; j < num_par-num_par_init; j++) {
                if (cc.trace_check(j, false)) {
                    for (int i = 0; i < local_num_comp; i++) {
#ifdef DEVELOP
                        std::cout << "traj: " << traj << " at " << flushed_snapshots
                                  << " Manual param " << i << "/" << local_num_comp
                                  << " gradient " << j << "/" << num_par-num_par_init << "\n";
#endif
                        startp2[0] = i;
                        SUCCESS_OR_DIE(
                            nc_put_vara(
                                ncid,
                                varid[Var_idx::n_vars + j],
                                startp2.data(),
                                countp2.data(),
                                output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                    }
                }
            }
        } else {
            // initial conditions
            for (uint64_t j=0; j < num_par_init; j++) {
                if (cc.trace_check(j, 2)) {
                    for (int i = 0; i < local_num_comp; i++) {
                        startp2[0] = i;
#ifdef DEVELOP
                        std::cout << "traj: " << traj << " at " << flushed_snapshots
                                  << " Manual IC param " << i << "/" << local_num_comp
                                  << " gradient " << j << "/" << num_par-num_par_init << "\n";
#endif
                        SUCCESS_OR_DIE(
                            nc_put_vara(
                                ncid,
                                varid[Var_idx::n_vars + j],
                                startp2.data(),
                                countp2.data(),
                                output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                    }
                }
            }
        }
    }
#else
    } else if (cc.traj_id == 0) {
        // no compression: limited_time_ensembles with and without strided output
        // Only traj_id = 0 writes gradients
        std::vector<size_t> startp2, countp2;
        startp2.push_back(0);
        if (no_flush) {
            startp2.push_back(0);
            countp2.push_back(0);
            countp2.push_back(0);
        } else {
            startp2.push_back(flushed_snapshots);
            countp2.push_back(local_num_comp);
            countp2.push_back(n_snapshots);
        }

        if (countp[1] != total_snapshots) {
            std::vector<std::ptrdiff_t> stridep, imap;
            stridep.push_back(1);
            stridep.push_back(1);
            if (no_flush) {
                imap.push_back(0);
            } else {
                imap.push_back(total_snapshots);
            }
            imap.push_back(1);

            if (!track_ic) {
                for (uint64_t j=0; j < num_par-num_par_init; j++) {
                    if (cc.trace_check(j, false))
                        SUCCESS_OR_DIE(
                            nc_put_varm(
                                ncid,
                                varid[Var_idx::n_vars + j],
                                startp2.data(),
                                countp2.data(),
                                stridep.data(),
                                imap.data(),
                                output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                }
            } else {
                // initial conditions
                for (uint64_t j=0; j < num_par_init; j++) {
                    if (cc.trace_check(j, 2))
                        SUCCESS_OR_DIE(
                            nc_put_varm(
                                ncid,
                                varid[Var_idx::n_vars + j],
                                startp2.data(),
                                countp2.data(),
                                stridep.data(),
                                imap.data(),
                                output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                }
            }
        } else {
            if (!track_ic) {
                for (uint64_t j=0; j < num_par-num_par_init; j++) {
                    if (cc.trace_check(j, false))
                        SUCCESS_OR_DIE(
                            nc_put_vara(
                                ncid,
                                varid[Var_idx::n_vars + j],
                                startp2.data(),
                                countp2.data(),
                                output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                }
            } else {
                // initial conditions
                for (uint64_t j=0; j < num_par_init; j++) {
                    if (cc.trace_check(j, 2))
                        SUCCESS_OR_DIE(
                            nc_put_vara(
                                ncid,
                                varid[Var_idx::n_vars + j],
                                startp2.data(),
                                countp2.data(),
                                output_buffer[static_cast<uint32_t>(Buffer_idx::n_buffer) + j].data()));
                }
            }
        }
    }
#endif
    if (!no_flush) flushed_snapshots += n_snapshots;
    n_snapshots = 0;
    return true;
}


template<class float_t>
void output_handle_t::process_step(
    const model_constants_t<float_t> &cc,
    const netcdf_reader_t &netcdf_reader,
    const std::vector<float_t> &y_single_new,
    const std::vector< std::array<double, num_par > >  &y_diff,
    const uint32_t sub,
    const uint32_t t,
    const uint32_t write_index,
    const uint32_t snapshot_index,
    const bool last_step,
    const reference_quantities_t &ref_quant) {

    if ((0 == (sub + t*cc.num_sub_steps) % snapshot_index)
        || (t == cc.num_steps-1 && last_step)) {
        this->buffer(cc, netcdf_reader, y_single_new, y_diff, sub, t,
            ref_quant, snapshot_index);
    } else if (this->simulation_mode != limited_time_ensembles || cc.traj_id == 0) {
        // Why do we have this case here? This would store gradients at every step.
        this->buffer_gradient(cc, y_diff, snapshot_index);
    }

    if (((0 == (sub + t*cc.num_sub_steps) % write_index)
        && (sub != 0) && (t != 0))
        || (t == cc.num_steps-1 && last_step)) {
        this->flush_buffer(cc);
    }
}

template output_handle_t::output_handle_t<codi::RealReverse>(
    const std::string,
    const std::string,
    const model_constants_t<codi::RealReverse>&,
    const std::string,
    const uint32_t,
    const uint32_t ,
    const int&,
    const int&,
    const bool&,
#ifdef COMPRESS_OUTPUT
    const int&,
#endif
    const double);

template output_handle_t::output_handle_t<codi::RealForwardVec<num_par_init> >(
    const std::string,
    const std::string,
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::string,
    const uint32_t,
    const uint32_t,
    const int&,
    const int&,
    const bool&,
#ifdef COMPRESS_OUTPUT
    const int&,
#endif
    const double);

template output_handle_t::output_handle_t<codi::RealReverse>(
        const std::string,
        const std::string,
        const model_constants_t<codi::RealReverse>&,
        const std::string,
        const uint32_t,
        const uint32_t ,
        const int&,
        const int&,
        const bool&,
#ifdef COMPRESS_OUTPUT
        const int&,
#endif
        const double,
        const std::vector<segment_t>&);

template output_handle_t::output_handle_t<codi::RealForwardVec<num_par_init> >(
        const std::string,
        const std::string,
        const model_constants_t<codi::RealForwardVec<num_par_init> >&,
        const std::string,
        const uint32_t,
        const uint32_t,
        const int&,
        const int&,
        const bool&,
#ifdef COMPRESS_OUTPUT
        const int&,
#endif
        const double,
        const std::vector<segment_t>&);

template void output_handle_t::setup_gradients<codi::RealReverse>(
        const model_constants_t<codi::RealReverse> &);

template void output_handle_t::setup_gradients<codi::RealForwardVec<num_par_init> >(
        const model_constants_t<codi::RealForwardVec<num_par_init> > &);

template void output_handle_t::define_var_gradients<codi::RealReverse>(
        const model_constants_t<codi::RealReverse> &,
        const int*,
        const int&);

template void output_handle_t::define_var_gradients<codi::RealForwardVec<num_par_init> >(
        const model_constants_t<codi::RealForwardVec<num_par_init> > &,
        const int*,
        const int&);

template void output_handle_t::define_vars<codi::RealReverse>(
        const model_constants_t<codi::RealReverse> &);

template void output_handle_t::define_vars<codi::RealForwardVec<num_par_init> >(
        const model_constants_t<codi::RealForwardVec<num_par_init> > &);

template void output_handle_t::set_attributes<codi::RealReverse>(
        const model_constants_t<codi::RealReverse> &,
        const std::string);

template void output_handle_t::set_attributes<codi::RealForwardVec<num_par_init> >(
        const model_constants_t<codi::RealForwardVec<num_par_init> > &,
        const std::string);

template void output_handle_t::set_compression<codi::RealReverse>(
        const model_constants_t<codi::RealReverse> &);

template void output_handle_t::set_compression<codi::RealForwardVec<num_par_init> >(
        const model_constants_t<codi::RealForwardVec<num_par_init> > &);

template void output_handle_t::write_dimension_values<codi::RealReverse>(
        const model_constants_t<codi::RealReverse> &,
        const double);

template void output_handle_t::write_dimension_values<codi::RealForwardVec<num_par_init> >(
        const model_constants_t<codi::RealForwardVec<num_par_init> > &,
        const double);

template void output_handle_t::set_parallel_access<codi::RealReverse>(
        const model_constants_t<codi::RealReverse> &,
        const std::string);

template void output_handle_t::set_parallel_access<codi::RealForwardVec<num_par_init> >(
        const model_constants_t<codi::RealForwardVec<num_par_init> > &,
        const std::string);

template void output_handle_t::setup<codi::RealReverse>(
    const std::string,
    const std::string,
    const model_constants_t<codi::RealReverse>&,
    const std::string,
    const uint32_t,
    const uint32_t,
    const int&,
    const double);

template void output_handle_t::setup<codi::RealForwardVec<num_par_init> >(
    const std::string,
    const std::string,
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::string,
    const uint32_t,
    const uint32_t,
    const int&,
    const double);

template void output_handle_t::buffer_gradient<codi::RealReverse>(
    const model_constants_t<codi::RealReverse>&,
    const std::vector< std::array<double, num_par > >&,
    const uint32_t);

template void output_handle_t::buffer_gradient<codi::RealForwardVec<num_par_init> >(
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const std::vector< std::array<double, num_par > >&,
    const uint32_t);

template void output_handle_t::buffer<codi::RealReverse>(
    const model_constants_t<codi::RealReverse>&,
    const netcdf_reader_t&,
    const std::vector<codi::RealReverse>&,
    const std::vector< std::array<double, num_par > >&,
    const uint32_t,
    const uint32_t,
    const reference_quantities_t&,
    const uint32_t);

template void output_handle_t::buffer<codi::RealForwardVec<num_par_init> >(
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const netcdf_reader_t&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    const std::vector< std::array<double, num_par > >&,
    const uint32_t,
    const uint32_t,
    const reference_quantities_t&,
    const uint32_t);

template bool output_handle_t::flush_buffer<codi::RealReverse>(
    const model_constants_t<codi::RealReverse>&, const bool);

template bool output_handle_t::flush_buffer<codi::RealForwardVec<num_par_init> >(
    const model_constants_t<codi::RealForwardVec<num_par_init> >&, const bool);

template void output_handle_t::process_step<codi::RealReverse>(
    const model_constants_t<codi::RealReverse>&,
    const netcdf_reader_t&,
    const std::vector<codi::RealReverse>&,
    const std::vector< std::array<double, num_par > >&,
    const uint32_t,
    const uint32_t,
    const uint32_t,
    const uint32_t,
    const bool,
    const reference_quantities_t&);

template void output_handle_t::process_step<codi::RealForwardVec<num_par_init> >(
    const model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const netcdf_reader_t&,
    const std::vector<codi::RealForwardVec<num_par_init> >&,
    const std::vector< std::array<double, num_par > >&,
    const uint32_t,
    const uint32_t,
    const uint32_t,
    const uint32_t,
    const bool,
    const reference_quantities_t&);
