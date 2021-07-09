#include "include/types/output_handle_t.h"

output_handle_t::output_handle_t()
{

}

output_handle_t::output_handle_t(
    const std::string filetype,
    const std::string filename,
    const model_constants_t &cc,
    const reference_quantities_t &ref_quant,
    const std::string in_filename,
    const uint32_t write_index,
    const uint32_t snapshot_index,
    const int &rank,
    const int &simulation_mode)
{
    this->simulation_mode = simulation_mode;
    local_num_comp = cc.local_num_comp;
    local_num_par = cc.local_num_par;

    this->setup(filetype, filename, cc, ref_quant,
        in_filename, write_index, snapshot_index, rank);
}

void output_handle_t::setup(
    const std::string filetype,
    const std::string filename,
    const model_constants_t &cc,
    const reference_quantities_t &ref_quant,
    const std::string in_filename,
    const uint32_t write_index,
    const uint32_t snapshot_index,
    const int &rank)
{
    this->n_trajs_file = cc.max_n_trajs;
    this->traj = cc.traj_id;
    this->ens = cc.ensemble_id;
    this->num_ens = cc.n_ensembles;
    this->num_time = cc.num_steps * cc.num_sub_steps;

    this->filetype = filetype;
    this->filename = filename;
    dimid.resize(Dim_idx::n_dims);
    varid.resize(Var_idx::n_vars + num_par);

    flushed_snapshots = 0;
    n_snapshots = 0;
    // Allocate memory for the buffer
    // maximum number of snapshots we are going to get
    total_snapshots = std::ceil( ((float)write_index)/snapshot_index ) + 1;
    const uint64_t vec_size = total_snapshots; // n_snapshots * num_comp;
    const uint64_t vec_size_grad = num_comp * total_snapshots;
    for(uint32_t i=0; i<num_comp; i++)
        output_buffer[i].resize(vec_size);
    for(uint32_t i=Buffer_idx::n_buffer; i<Buffer_idx::n_buffer+num_par; i++)
        output_buffer[i].resize(vec_size_grad);

    output_buffer[Buffer_idx::time_ascent_buf].resize(vec_size);
    // output_buffer[num_comp+num_par+1].resize(total_snapshots); // just time index
    output_buffer[Buffer_idx::lat_buf].resize(vec_size);        // lat
    output_buffer[Buffer_idx::lon_buf].resize(vec_size);        // lon

    for(uint32_t i=0; i<output_buffer_flags.size(); i++)
        output_buffer_flags[i].resize(vec_size);

    for(uint32_t i=0; i<output_buffer_int.size(); i++)
        output_buffer_int[i].resize(vec_size);

    // Unfortunately, it is likely to run into an HDF Error if
    // many parameters are investigated in combination with the
    // amount of attributes for the time dimension.
    // We have to create the file as a single process, close it and
    // open it in parallel again for writing purpose
    if(rank == 0)
    {
        SUCCESS_OR_DIE(nc_create(
            (filename + ".nc_wcb").c_str(), // path
            NC_NETCDF4,                     // creation mode
            &ncid)
        );
        // Create dimensions
        SUCCESS_OR_DIE(nc_def_dim(
            ncid, "Output_Parameter", local_num_comp, &dimid[Dim_idx::out_param_dim])
        );
        SUCCESS_OR_DIE(nc_def_dim(
            ncid, "ensemble", num_ens, &dimid[Dim_idx::ensemble_dim])
        );
        SUCCESS_OR_DIE(nc_def_dim(
            ncid, "trajectory", n_trajs_file, &dimid[Dim_idx::trajectory_dim])
        );
        SUCCESS_OR_DIE(nc_def_dim(
            ncid, "time", num_time, &dimid[Dim_idx::time_dim])
        );

        // Create variables
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "Output_Parameter_ID",
            NC_UINT64,  // type
            1,          // ndims 2: matrix, 1: vector, 0: scalar
            &dimid[Dim_idx::out_param_dim],    // ignored for ndims = 0
            &varid[Var_idx::out_param])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "ensemble",
            NC_UINT64,
            1,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::ensemble])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "trajectory",
            NC_UINT64,
            1,
            &dimid[Dim_idx::trajectory_dim],
            &varid[Var_idx::trajectory])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "time",
            NC_DOUBLE,
            1,
            &dimid[Dim_idx::time_dim],
            &varid[Var_idx::time])
        );

        // model state
        for(uint32_t i=0; i<num_comp; ++i)
            SUCCESS_OR_DIE(
                nc_def_var(
                    ncid,
                    output_par_idx[i].c_str(),
                    NC_DOUBLE,
                    3,
                    &dimid[Dim_idx::ensemble_dim],
                    &varid[i]
                )
            );
        // gradients
        for(uint32_t i=0; i<output_grad_idx.size(); ++i)
        {
            if(cc.trace_check(i, false))
                SUCCESS_OR_DIE(nc_def_var(
                    ncid,
                    output_grad_idx[i].c_str(),
                    NC_DOUBLE,
                    4,
                    &dimid[Dim_idx::out_param_dim],
                    &varid[Var_idx::n_vars + i])
                );
        }

        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "time_after_ascent",
            NC_DOUBLE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::time_ascent])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "conv_400",
            NC_BYTE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::conv_400])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "conv_600",
            NC_BYTE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::conv_600])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "slan_400",
            NC_BYTE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::slan_400])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "slan_600",
            NC_BYTE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::slan_600])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "lat",
            NC_DOUBLE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::lat])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "lon",
            NC_DOUBLE,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::lon])
        );
        SUCCESS_OR_DIE(nc_def_var(
            ncid,
            "step",
            NC_UINT64,
            3,
            &dimid[Dim_idx::ensemble_dim],
            &varid[Var_idx::step])
        );

        // Add attributes; read attributes from in_filename first
        int in_ncid;
        SUCCESS_OR_DIE(nc_open(in_filename.c_str(), NC_NOWRITE, &in_ncid));

        // get amount of attributes
        int n_atts;
        SUCCESS_OR_DIE(nc_inq(in_ncid, NULL, NULL, &n_atts, NULL));

        // for every attribute, get the name, type and length, and the values
        for(int i=0; i<n_atts; i++)
        {
            char att_name[NC_MAX_NAME];
            SUCCESS_OR_DIE(nc_inq_attname(in_ncid, NC_GLOBAL, i, att_name));
            nc_type att_type;
            size_t att_len;
            SUCCESS_OR_DIE(nc_inq_att(
                in_ncid, NC_GLOBAL, att_name, &att_type, &att_len)
            );
            std::vector<nc_vlen_t> att_val(att_len);
            SUCCESS_OR_DIE(nc_get_att(
                in_ncid, NC_GLOBAL, att_name, att_val.data())
            );
            SUCCESS_OR_DIE(nc_put_att(
                ncid, NC_GLOBAL, att_name, att_type, att_len, att_val.data())
            );
        }
        // time
        int in_time_id;
        SUCCESS_OR_DIE(nc_inq_varid(in_ncid, "time", &in_time_id));
        SUCCESS_OR_DIE(nc_inq_varnatts(in_ncid, in_time_id, &n_atts));
        for(int i=0; i<n_atts; i++)
        {
            char att_name[NC_MAX_NAME];
            SUCCESS_OR_DIE(nc_inq_attname(in_ncid, in_time_id, i, att_name));
            if(!std::strcmp(att_name, "_FillValue"))
            {
                continue;
            }
            // nc_type att_type;
            size_t att_len;
            SUCCESS_OR_DIE(nc_inq_att(
                in_ncid, in_time_id, att_name, NULL, &att_len)
            );
            char att_val[att_len];
            SUCCESS_OR_DIE(nc_get_att(
                in_ncid, in_time_id, att_name, att_val)
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time],
                att_name,
                att_len,
                att_val)
            );
        }
        SUCCESS_OR_DIE(ncclose(in_ncid));

        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::out_param],
            "long_name",
            strlen("gradients are calculated w.r.t. this output parameter"),
            "gradients are calculated w.r.t. this output parameter")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::out_param],
            "standard_name",
            strlen("output_parameter_id"),
            "output_parameter_id")
        );

        auto put_att_mass = [&](
            const char *mass_name,
            const char *long_mass_name,
            auto &varid)
        {
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "units",
                strlen("kg m^-3"),
                "kg m^-3")
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "auxiliary_data",
                strlen("yes"),
                "yes")
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "long_name",
                strlen(long_mass_name),
                long_mass_name)
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "standard_name",
                strlen(mass_name),
                mass_name)
            );
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
            auto &varid)
        {
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "units",
                strlen("m^-3"),
                "m^-3")
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "auxiliary_data",
                strlen("yes"),
                "yes")
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "long_name",
                strlen(long_name),
                long_name)
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "standard_name",
                strlen(name),
                name)
            );
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
            auto &varid)
        {
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "units",
                strlen("kg m^-3 s^-1"),
                "kg m^-3 s^-1")
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "auxiliary_data",
                strlen("yes"),
                "yes")
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "long_name",
                strlen(long_mass_name),
                long_mass_name)
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "standard_name",
                strlen(mass_name),
                mass_name)
            );
        };
        put_att_mass_sed(
            "sedi_outflux_of_rain_droplet_mass",
            "sedimentation of rain droplet mass",
            varid[Var_idx::qr_out]
        );
        put_att_mass_sed(
            "sedi_outflux_of_snow_mass",
            "sedimentation of snow mass",
            varid[Var_idx::qs_out]
        );
        put_att_mass_sed(
            "sedi_outflux_of_ice_mass",
            "sedimentation of ice mass",
            varid[Var_idx::qi_out]
        );
        put_att_mass_sed(
            "sedi_outflux_of_graupel_mass",
            "sedimentation of graupel mass",
            varid[Var_idx::qg_out]
        );
        put_att_mass_sed(
            "sedi_outflux_of_hail_mass",
            "sedimentation of hail mass",
            varid[Var_idx::qh_out]
        );

        auto put_att_nums_sed = [&](
            const char *name,
            const char *long_name,
            auto &varid)
        {
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "units",
                strlen("m^-3 s^-1"),
                "m^-3 s^-1")
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "auxiliary_data",
                strlen("yes"),
                "yes")
            );
            std::string tmp_string = "sedimentation of " + std::string(long_name) + " number";
            const char *att_val = tmp_string.c_str();
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "long_name",
                strlen(att_val),
                att_val)
            );
            std::string tmp_string_2 = "sedi_outflux_of_" + std::string(name) + "_number";
            const char *att_val_2 = tmp_string_2.c_str();
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid,
                "standard_name",
                strlen(att_val_2),
                att_val_2)
            );
        };
        put_att_nums_sed("sedi_outflux_of_rain_droplet_number", "sedimentation of rain droplet number", varid[Var_idx::nr_out]);
        put_att_nums_sed("sedi_outflux_of_snow_number", "sedimentation of snow number", varid[Var_idx::ns_out]);
        put_att_nums_sed("sedi_outflux_of_ice_number", "sedimentation of ice number", varid[Var_idx::ni_out]);
        put_att_nums_sed("sedi_outflux_of_graupel_number", "sedimentation of graupel number", varid[Var_idx::ng_out]);
        put_att_nums_sed("sedi_outflux_of_hail_number", "sedimentation of hail number", varid[Var_idx::nh_out]);

        // all gradients are auxiliary data
        for(int i=0; i<num_par; i++)
        {
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::n_vars + i],
                "auxiliary_data",
                strlen("yes"),
                "yes")
            );
        }
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "long_name",
            strlen("pressure"),
            "pressure")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "standard_name",
            strlen("air_pressure"),
            "air_pressure")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "units",
            strlen("Pa"),
            "Pa")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "positive",
            strlen("down"),
            "down")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::pressure],
            "axis",
            strlen("Z"),
            "Z")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::temperature],
            "long_name",
            strlen("temperature"),
            "temperature")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::temperature],
            "standard_name",
            strlen("air_temperature"),
            "air_temperature")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::temperature],
            "units",
            strlen("K"),
            "K")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::temperature],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::ascent],
            "long_name",
            strlen("ascend velocity"),
            "ascend velocity")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::ascent],
            "standard_name",
            strlen("ascend_velocity"),
            "ascend_velocity")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::ascent],
            "units",
            strlen("m s^-1"),
            "m s^-1")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::ascent],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sat],
            "long_name",
            strlen("saturation"),
            "saturation")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sat],
            "standard_name",
            strlen("saturation"),
            "saturation")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sat],
            "units",
            strlen("percentage"),
            "percentage")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sat],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::height],
            "long_name",
            strlen("height above mean sea level"),
            "height above mean sea level")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::height],
            "standard_name",
            strlen("height"),
            "height")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::height],
            "units",
            strlen("m AMSL"),
            "m AMSL")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::height],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::inactive],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::dep],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::sub],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat_heat],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat_cool],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::time_ascent],
            "long_name",
            strlen("time after rapid ascent started"),
            "time after rapid ascent started")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::time_ascent],
            "standard_name",
            strlen("time_after_ascent"),
            "time_after_ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::time_ascent],
            "units",
            strlen("seconds since start of convective/slantwise ascent"),
            "seconds since start of convective/slantwise ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::time_ascent],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat],
            "long_name",
            strlen("rotated latitude"),
            "rotated latitude")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat],
            "standard_name",
            strlen("latitude"),
            "latitude")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lat],
            "units",
            strlen("degrees"),
            "degrees")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lon],
            "long_name",
            strlen("rotated longitude"),
            "rotated longitude")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lon],
            "standard_name",
            strlen("longitude"),
            "longitude")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::lon],
            "units",
            strlen("degrees"),
            "degrees")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_400],
            "long_name",
            strlen("convective 400hPa ascent"),
            "convective 400hPa ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_400],
            "standard_name",
            strlen("convective_400hPa_ascent"),
            "convective_400hPa_ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_400],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_600],
            "long_name",
            strlen("convective 600hPa ascent"),
            "convective 600hPa ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_600],
            "standard_name",
            strlen("convective_600hPa_ascent"),
            "convective_600hPa_ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::conv_600],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_400],
            "long_name",
            strlen("slantwise 400hPa ascent"),
            "slantwise 400hPa ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_400],
            "standard_name",
            strlen("slantwise_400hPa_ascent"),
            "slantwise_400hPa_ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_400],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_600],
            "long_name",
            strlen("slantwise 600hPa ascent"),
            "slantwise 600hPa ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_600],
            "standard_name",
            strlen("slantwise_600hPa_ascent"),
            "slantwise_600hPa_ascent")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::slan_600],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::step],
            "long_name",
            strlen("simulation step"),
            "simulation step")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::step],
            "standard_name",
            strlen("step"),
            "step")
        );
        SUCCESS_OR_DIE(nc_put_att_text(
            ncid,
            varid[Var_idx::step],
            "auxiliary_data",
            strlen("yes"),
            "yes")
        );
        // in theory, one could apply compression here
        // if(!(this->simulation_mode == trajectory_sensitvity_perturbance)
        //     && !(this->simulation_mode == trajectory_perturbance))
        // {
        //     for(auto &id: varid)
        //         SUCCESS_OR_DIE(
        //             nc_def_var_deflate(
        //                 ncid,
        //                 id,
        //                 1, // shuffle
        //                 1, // deflate
        //                 9 // max compression
        //             )
        //         );
        // }

        SUCCESS_OR_DIE(nc_enddef(ncid));
        // Write dimensions here
        std::vector<size_t> startp, countp;
        startp.push_back(0);
        countp.push_back(num_time);
        std::vector<double> time_steps(num_time);
        for(uint32_t i=0; i<num_time; i++)
            time_steps[i] = cc.dt*i + cc.start_time; // start + i*dt
        SUCCESS_OR_DIE(
            nc_put_vara_double(
                ncid,                       // ncid
                varid[Var_idx::time],       // varid
                startp.data(),              // startp
                countp.data(),              // countp
                time_steps.data()           // op
            )
        );

        countp[0] = n_trajs_file;
        std::vector<uint64_t> data(n_trajs_file);
        for(uint32_t i=0; i<n_trajs_file; i++)
            data[i] = i;
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid[Var_idx::trajectory],
                startp.data(),
                countp.data(),
                data.data()
            )
        );

        countp[0] = num_ens;
        data.resize(num_ens);
        for(uint32_t i=0; i<num_ens; i++)
            data[i] = i;
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid[Var_idx::ensemble],
                startp.data(),
                countp.data(),
                data.data()
            )
        );

        countp[0] = local_num_comp;
        data.resize(local_num_comp);
        uint32_t counter = 0;
        for(uint32_t i=0; i<num_comp; i++)
        {
            if(cc.trace_check(i, true))
            {
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
                data.data()
            )
        );
        SUCCESS_OR_DIE(ncclose(ncid));
    }
    // Open it again for writing with all processes
    std::string file_string = filename + ".nc_wcb";
    SUCCESS_OR_DIE(
        nc_open_par(
            file_string.c_str(),
            NC_WRITE,
            MPI_COMM_WORLD,
            MPI_INFO_NULL,
            &ncid
        )
    );

    // gather all necessary variable ids
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "Output_Parameter_ID",
            &varid[Var_idx::out_param]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "ensemble",
            &varid[Var_idx::ensemble]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "trajectory",
            &varid[Var_idx::trajectory]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "time",
            &varid[Var_idx::time]
        )
    );

    // model state
    for(uint32_t i=0; i<num_comp; ++i)
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                output_par_idx[i].c_str(),
                &varid[i]
            )
        );

    // gradients
    for(uint32_t i=0; i<output_grad_idx.size(); ++i)
    {
        if(cc.trace_check(i, false))
        {
            SUCCESS_OR_DIE(
                nc_inq_varid(
                    ncid,
                    output_grad_idx[i].c_str(),
                    &varid[Var_idx::n_vars + i]
                )
            );
        }
    }

    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "time_after_ascent",
            &varid[Var_idx::time_ascent]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "conv_400",
            &varid[Var_idx::conv_400]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "conv_600",
            &varid[Var_idx::conv_600]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "slan_400",
            &varid[Var_idx::slan_400]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "slan_600",
            &varid[Var_idx::slan_600]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "lat",
            &varid[Var_idx::lat]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "lon",
            &varid[Var_idx::lon]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "step",
            &varid[Var_idx::step]
        )
    );

    if((this->simulation_mode == trajectory_sensitvity_perturbance)
        || (this->simulation_mode == trajectory_perturbance))
    {
        // Make the access independent which is a must due to the dynamic
        // work schedule; This can be expensive though.
        for(uint32_t i=0; i<Var_idx::n_vars; i++)
            SUCCESS_OR_DIE(nc_var_par_access(ncid, varid[i], NC_INDEPENDENT));
        for(uint32_t i=0; i<output_grad_idx.size(); i++)
            if(cc.trace_check(i+Var_idx::n_vars, false))
                SUCCESS_OR_DIE(nc_var_par_access(ncid, varid[i+Var_idx::n_vars], NC_INDEPENDENT));
        // for(auto &id: varid)
        //     SUCCESS_OR_DIE(nc_var_par_access(ncid, id, NC_INDEPENDENT));
    }
}


void output_handle_t::reset(
    const uint32_t traj_id,
    const uint32_t ens_id)
{
    flushed_snapshots = 0;
    n_snapshots = 0;
    this->traj = traj_id;
    this->ens = ens_id;
}


void output_handle_t::buffer(
    const model_constants_t &cc,
    const netcdf_reader_t &netcdf_reader,
    // const nc_parameters_t &nc_params,
    const std::vector<codi::RealReverse> &y_single_new,
    const std::vector< std::array<double, num_par > >  &y_diff,
    const uint32_t sub,
    const uint32_t t,
    const double time_new,
    const uint32_t traj_id,
    const uint32_t ensemble,
    const reference_quantities_t &ref_quant)
{
    const uint64_t offset = 1;

    // output parameters
    for(uint64_t i=0; i<num_comp; i++)
    {
        switch(i)
        {
            case p_idx:
                output_buffer[i][n_snapshots*offset] =
                    y_single_new[i].getValue() * ref_quant.pref;
                break;
            case T_idx:
                output_buffer[i][n_snapshots*offset] =
                    y_single_new[i].getValue() * ref_quant.Tref;
                break;
            case w_idx:
                output_buffer[i][n_snapshots*offset] =
                    y_single_new[i].getValue() * ref_quant.wref;
                break;
            case z_idx:
                output_buffer[i][n_snapshots*offset] =
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
                output_buffer[i][n_snapshots*offset] =
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
                output_buffer[i][n_snapshots*offset] =
                    y_single_new[i].getValue() * ref_quant.Nref;
                break;
            default:
                output_buffer[i][n_snapshots*offset] =
                    y_single_new[i].getValue();
                break;
        }
    }

    for(uint64_t i=0; i<num_comp; i++) // gradient sensitive to output parameter i
    {
        if(!cc.trace_check(i, true))
            continue;
        for(uint64_t j=0; j<num_par; j++) // gradient of input parameter j
            if(cc.trace_check(j, false))
                output_buffer[Buffer_idx::n_buffer+j][i + n_snapshots*offset*num_comp] = y_diff[i][j];
    }

    // lat
    output_buffer[Buffer_idx::lat_buf][n_snapshots*offset] =
        (netcdf_reader.get_lat(t) + sub*(netcdf_reader.get_lat(t+1)-netcdf_reader.get_lat(t)));

    // lon
    output_buffer[Buffer_idx::lon_buf][n_snapshots*offset] =
        (netcdf_reader.get_lon(t) + sub*(netcdf_reader.get_lon(t+1)-netcdf_reader.get_lon(t)));

#ifdef MET3D
    // time after ascent
    output_buffer[Buffer_idx::time_ascent_buf][n_snapshots*offset] =
        netcdf_reader.get_relative_time(t) + sub*cc.dt;
    // flags
    output_buffer_flags[0][n_snapshots*offset] = netcdf_reader.get_conv_400(t);
    output_buffer_flags[1][n_snapshots*offset] = netcdf_reader.get_conv_600(t);
    output_buffer_flags[2][n_snapshots*offset] = netcdf_reader.get_slan_400(t);
    output_buffer_flags[3][n_snapshots*offset] = netcdf_reader.get_slan_600(t);
#endif
    // simulation step
    output_buffer_int[0][n_snapshots*offset] = sub + t*cc.num_sub_steps;

    n_snapshots++;
}

void output_handle_t::flush_buffer(
    const model_constants_t &cc)
{
    std::vector<size_t> startp, countp;
    startp.push_back(ens);
    startp.push_back(traj);
    startp.push_back(flushed_snapshots);
    countp.push_back(1);
    countp.push_back(1);
    countp.push_back(n_snapshots);

    for(uint64_t i=0; i<num_comp; i++)
    {
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid[i],
                startp.data(),
                countp.data(),
                output_buffer[i].data()
            )
        );
    }
    // time after ascent
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::time_ascent],
            startp.data(),
            countp.data(),
            output_buffer[Buffer_idx::time_ascent_buf].data()
        )
    );
    // flags
    for(uint64_t i=0; i<output_buffer_flags.size(); i++)
    {
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid[Var_idx::conv_400+i],
                startp.data(),
                countp.data(),
                output_buffer_flags[i].data()
            )
        );

    }
    // lat
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::lat],
            startp.data(),
            countp.data(),
            output_buffer[Buffer_idx::lat_buf].data()
        )
    );
    // lon
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::lon],
            startp.data(),
            countp.data(),
            output_buffer[Buffer_idx::lon_buf].data()
        )
    );
    // step
    SUCCESS_OR_DIE(
        nc_put_vara(
            ncid,
            varid[Var_idx::step],
            startp.data(),
            countp.data(),
            output_buffer_int[0].data()
        )
    );
    // gradients
    startp.insert(startp.begin(), 0);
    countp.insert(countp.begin(), local_num_comp);
    for(uint64_t j=0; j<num_par; j++)
    {
        if(cc.trace_check(j, false))
            SUCCESS_OR_DIE(
                nc_put_vara(
                    ncid,
                    varid[Var_idx::n_vars + j],
                    startp.data(),
                    countp.data(),
                    output_buffer[Buffer_idx::n_buffer + j].data()
                )
            );
    }
    flushed_snapshots += n_snapshots;
    n_snapshots = 0;
}

void output_handle_t::process_step(
    const model_constants_t &cc,
    const netcdf_reader_t &netcdf_reader,
    // const nc_parameters_t &nc_params,
    const std::vector<codi::RealReverse> &y_single_new,
    const std::vector< std::array<double, num_par > >  &y_diff,
    const uint32_t sub,
    const uint32_t t,
    const double time_new,
    const uint32_t traj_id,
    const uint32_t write_index,
    const uint32_t snapshot_index,
#ifdef MET3D
    const uint32_t ensemble,
#endif
    const bool last_step,
    const reference_quantities_t &ref_quant)
{
    if( (0 == (sub + t*cc.num_sub_steps) % snapshot_index)
        || ( t == cc.num_steps-1 && last_step ) )
    {
        this->buffer(cc, netcdf_reader, y_single_new, y_diff, sub, t,
            time_new, traj_id, ensemble, ref_quant);
    }

    if( (0 == (sub + t*cc.num_sub_steps) % write_index)
        || ( t == cc.num_steps-1 && last_step ) )
    {
        this->flush_buffer(cc);
    }
}

