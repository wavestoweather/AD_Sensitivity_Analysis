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
    const int &rank)
{
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
    this->num_time = cc.num_steps * cc.num_sub_steps + 1; // + 1 for initial time

    this->filetype = filetype;
    this->filename = filename;
    dimid.resize(Dim_idx::n_dims);
    varid.resize(output_par_idx.size() + output_grad_idx.size() + Dim_idx::n_dims + 9);

    if(filetype == "netcdf")
    {
        const int deflateLevel = 9; // compression level
        const bool enableShuffleFilter = true; // increases compression with little cost
        const bool enableDeflateFilter = true; // enable compression
        flushed_snapshots = 0;
        n_snapshots = 0;
        // Allocate memory for the buffer
        // maximum number of snapshots we are going to get
        total_snapshots = std::ceil( ((float)write_index)/snapshot_index ) + 1;
        const uint64_t vec_size = total_snapshots; // n_snapshots * num_comp;
        const uint64_t vec_size_grad = num_comp * total_snapshots;
        for(uint32_t i=0; i<num_comp; i++)
            output_buffer[i].resize(vec_size);
        for(uint32_t i=num_comp; i<num_comp+num_par; i++)
            output_buffer[i].resize(vec_size_grad);

        output_buffer[num_comp+num_par].resize(vec_size);          // time after ascent
        output_buffer[num_comp+num_par+1].resize(total_snapshots); // just time index
        output_buffer[num_comp+num_par+2].resize(vec_size);        // lat
        output_buffer[num_comp+num_par+3].resize(vec_size);        // lon

        for(uint32_t i=0; i<output_buffer_flags.size(); i++)
            output_buffer_flags[i].resize(vec_size);
        // for(uint32_t i=0; i<output_buffer_str.size(); i++)
        //     output_buffer_str[i].resize(vec_size);
        for(uint32_t i=0; i<output_buffer_int.size(); i++)
            output_buffer_int[i].resize(vec_size);

        // auto fill_vectors = [&]()
        // {
        //     var_vector.clear();
        //     // File exists. No need to set it up. Load the columns
        //     bool got_it = false;
        //     do{
        //         try
        //         {
        //             datafile.open(filename + ".nc_wcb", NcFile::read);
        //             got_it = true;
        //         } catch(const std::exception& e2)
        //         {
        //             got_it = false;
        //         }
        //     } while(!got_it);

        //     // std::cout << "opened\n";
        //     for(auto &out_par: output_par_idx)
        //     {
        //         var_vector.push_back(datafile.getVar(out_par));
        //     }

        //     for(auto &out_grad: output_grad_idx)
        //     {
        //         var_vector.push_back(datafile.getVar(out_grad));
        //     }
        //     var_vector.push_back(datafile.getVar("time_after_ascent"));
        //     var_vector.push_back(datafile.getVar("time"));
        //     var_vector.push_back(datafile.getVar("conv_400"));
        //     var_vector.push_back(datafile.getVar("conv_600"));
        //     var_vector.push_back(datafile.getVar("slan_400"));
        //     var_vector.push_back(datafile.getVar("slan_600"));
        //     var_vector.push_back(datafile.getVar("type"));
        //     var_vector.push_back(datafile.getVar("lat"));
        //     var_vector.push_back(datafile.getVar("lon"));
        //     var_vector.push_back(datafile.getVar("step"));
        //     datafile.close();
        // };
        // if(traj > 0)
        // {
        //     fill_vectors();
        //     return;
        // }

        // int ncid;
        // Must be called by all processes defined in the communicator
        SUCCESS_OR_DIE(nc_create_par(
            (filename + ".nc_wcb").c_str(), // path
            NC_NETCDF4,           // creation mode
            MPI_COMM_WORLD,                 // MPI communicator
            MPI_INFO_NULL,        // MPI_Info
            &ncid)
        );
        // if(rank == 0)
        // {
            std::cout << "Creating file with " << num_comp << " output parameters\n";
            std::cout << "with " << num_ens << " ensembles\n";
            std::cout << "with " << n_trajs_file << " trajectories\n";
            std::cout << "with " << num_time << " time steps\n";
            std::cout << "with " << total_snapshots << " total snapshots\n";
            // Create dimensions
            SUCCESS_OR_DIE(nc_def_dim(
                ncid, "Output_Parameter", num_comp, &dimid[Dim_idx::out_param_dim])
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
                NC_UINT64,       // type
                1,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::out_param_dim],    // ignored for ndims = 0
                &varid[Var_idx::out_param])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "ensemble",
                NC_UINT64,       // type
                1,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::ensemble_dim],    // ignored for ndims = 0
                &varid[Var_idx::ensemble])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "trajectory",
                NC_UINT64,       // type
                1,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::trajectory_dim],    // ignored for ndims = 0
                &varid[Var_idx::trajectory])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "time",
                NC_DOUBLE,       // type
                1,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::time])
            );
        // }
        // if(rank == -1)
        // {
            // uint64_t offset = Var_idx::out_param + 1;
            // model state
            for(uint32_t i=0; i<num_comp; ++i)
                SUCCESS_OR_DIE(
                    nc_def_var(
                        ncid,
                        output_par_idx[i].c_str(),
                        NC_DOUBLE,       // type
                        3,          // ndims 2: matrix, 1: vector, 0: scalar
                        &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                        &varid[i]
                    )
                );
            // gradients
            uint64_t offset = Var_idx::n_vars;
            for(uint32_t i=0; i<output_grad_idx.size(); ++i)
                SUCCESS_OR_DIE(nc_def_var(
                    ncid,
                    output_grad_idx[i].c_str(),
                    NC_DOUBLE,       // type
                    4,          // ndims 2: matrix, 1: vector, 0: scalar
                    &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                    &varid[offset + i])
                );

            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "time_after_ascent",
                NC_DOUBLE,       // type
                3,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::time_ascent])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "conv_400",
                NC_BYTE,       // type
                3,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::conv_400])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "conv_600",
                NC_BYTE,       // type
                3,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::conv_600])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "slan_400",
                NC_BYTE,       // type
                3,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::slan_400])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "slan_600",
                NC_BYTE,       // type
                3,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::slan_600])
            );
            // SUCCESS_OR_DIE(nc_def_var(
            //     ncid,
            //     "type",
            //     NC_STRING,       // type
            //     3,          // ndims 2: matrix, 1: vector, 0: scalar
            //     &dimid[Dim_idx::ensemble_dim],    // ignored for ndims = 0
            //     &varid[Var_idx::type])
            // );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "lat",
                NC_DOUBLE,       // type
                3,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::lat])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "lon",
                NC_DOUBLE,       // type
                3,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::lon])
            );
            SUCCESS_OR_DIE(nc_def_var(
                ncid,
                "step",
                NC_UINT64,       // type
                3,          // ndims 2: matrix, 1: vector, 0: scalar
                &dimid[Dim_idx::time_dim],    // ignored for ndims = 0
                &varid[Var_idx::step])
            );
        if(rank == -1)
        {
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
            ncclose(in_ncid);
            // column attributes
            std::string att_val = "gradients are calculated w.r.t. this output parameter";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::out_param],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "output_parameter";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::out_param],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            auto put_att_mass = [&](
                const char *mass_name,
                const char *long_mass_name,
                auto &varid)
            {
                att_val = "kg m^-3";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "units",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = "yes";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "auxiliary_data",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = std::string(long_mass_name) + " mass density";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "long_name",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = std::string(mass_name) + "_mass_density";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "standard_name",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );
            };

            put_att_mass("water_vapor", "water vapor", varid[Var_idx::qv]);
            put_att_mass("cloud_droplet", "cloud droplet", varid[Var_idx::qc]);
            put_att_mass("rain_droplet", "rain droplet", varid[Var_idx::qr]);
            put_att_mass("snow", "snow", varid[Var_idx::qs]);
            put_att_mass("ice", "ice", varid[Var_idx::qi]);
            put_att_mass("graupel", "graupel", varid[Var_idx::qg]);
            put_att_mass("hail", "hail", varid[Var_idx::qh]);

            auto put_att_nums = [&](
                const char *name,
                const char *long_name,
                auto &varid)
            {
                att_val = "m^-3";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "units",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = "yes";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "auxiliary_data",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = std::string(long_name) + " number density";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "long_name",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = std::string(name) + "_number_density";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "standard_name",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );
            };

            put_att_nums("cloud_droplet", "cloud droplet", varid[Var_idx::ncloud]);
            put_att_nums("rain_droplet", "rain droplet", varid[Var_idx::nrain]);
            put_att_nums("snow", "snow", varid[Var_idx::nsnow]);
            put_att_nums("ice", "ice", varid[Var_idx::nice]);
            put_att_nums("graupel", "graupel", varid[Var_idx::ngraupel]);
            put_att_nums("hail", "hail", varid[Var_idx::nhail]);

            auto put_att_mass_sed = [&](
                const char *mass_name,
                const char *long_mass_name,
                auto &varid)
            {
                att_val = "kg m^-3 s^-1";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "units",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = "yes";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "auxiliary_data",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = "sedimentation of " + std::string(long_mass_name) + " mass";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "long_name",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = "sedi_outflux_of_" + std::string(mass_name) + "_mass";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "standard_name",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );
            };

            put_att_mass_sed("rain_droplet", "rain droplet", varid[Var_idx::qr_out]);
            put_att_mass_sed("snow", "snow", varid[Var_idx::qs_out]);
            put_att_mass_sed("ice", "ice", varid[Var_idx::qi_out]);
            put_att_mass_sed("graupel", "graupel", varid[Var_idx::qg_out]);
            put_att_mass_sed("hail", "hail", varid[Var_idx::qh_out]);

            auto put_att_nums_sed = [&](
                const char *name,
                const char *long_name,
                auto &varid)
            {
                att_val = "m^-3 s^-1";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "units",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = "yes";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "auxiliary_data",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = "sedimentation of " + std::string(long_name) + " number";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "long_name",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );

                att_val = "sedi_outflux_of_" + std::string(name) + "_number";
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid,
                    "standard_name",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );
            };

            put_att_nums_sed("rain_droplet", "rain droplet", varid[Var_idx::nr_out]);
            put_att_nums_sed("snow", "snow", varid[Var_idx::ns_out]);
            put_att_nums_sed("ice", "ice", varid[Var_idx::ni_out]);
            put_att_nums_sed("graupel", "graupel", varid[Var_idx::ng_out]);
            put_att_nums_sed("hail", "hail", varid[Var_idx::nh_out]);

            // all gradients are auxiliary data
            offset = Var_idx::time + output_par_idx.size();
            att_val = "yes";
            for(int i=1; i<=output_grad_idx.size()+6; i++)
                SUCCESS_OR_DIE(nc_put_att_text(
                    ncid,
                    varid[offset + i],
                    "auxiliary_data",
                    strlen(att_val.c_str()),
                    att_val.c_str())
                );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[offset + output_grad_idx.size() + 9],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // Pressure
            att_val = "pressure";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::pressure],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "air_pressure";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::pressure],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "Pa";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::pressure],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "down";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::pressure],
                "positive",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "Z";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::pressure],
                "axis",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // Temperature
            att_val = "temperature";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::temperature],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "air_temperature";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::temperature],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "K";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::temperature],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::temperature],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // Ascent
            att_val = "ascend velocity";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::ascent],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "ascend_velocity";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::ascent],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "m s^-1";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::ascent],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::ascent],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // Saturation
            att_val = "saturation";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::sat],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::sat],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "percentage";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::sat],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::sat],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // height
            att_val = "height above mean sea level";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::height],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "height";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::height],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "m AMSL";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::height],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::height],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // inactive CCN number
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::inactive],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // depostion
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::dep],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // sublimination
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::sub],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // latent heat
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::lat_heat],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // latent cool
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::lat_cool],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // time
            att_val = "time";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "time";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "seconds since"; // TODO reference time
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = ""; // TODO reference time
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time],
                "forecast_inittime",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = ""; // TODO reference time
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time],
                "trajectory_starttime",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // time after ascent
            att_val = "time after rapid ascent started";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time_ascent],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "time_after_ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time_ascent],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "seconds since start of convective/slantwise ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time_ascent],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::time_ascent],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // lat
            att_val = "rotated latitude";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::lat],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "latitude";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::lat],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "degrees";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::lat],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // lon
            att_val = "rotated longitude";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::lat],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "longitude";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::lat],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "degrees";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::lat],
                "units",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // conv_400
            att_val = "convective 400hPa ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::conv_400],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "convective_400hPa_ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::conv_400],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::conv_400],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // conv_600
            att_val = "convective 600hPa ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::conv_600],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "convective_600hPa_ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::conv_600],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::conv_600],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // slan_400
            att_val = "slantwise 400hPa ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::slan_400],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "slantwise_400hPa_ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::slan_400],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::slan_400],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // slan_600
            att_val = "slantwise 600hPa ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::slan_600],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "slantwise_600hPa_ascent";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::slan_600],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::slan_600],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );

            // type
            // att_val = "trajectory type";
            // SUCCESS_OR_DIE(nc_put_att_text(
            //     ncid,
            //     varid[Var_idx::type],
            //     "long_name",
            //     strlen(att_val.c_str()),
            //     att_val.c_str())
            // );
            // att_val = "trajectory_type";
            // SUCCESS_OR_DIE(nc_put_att_text(
            //     ncid,
            //     varid[Var_idx::type],
            //     "standard_name",
            //     strlen(att_val.c_str()),
            //     att_val.c_str())
            // );
            // att_val = "yes";
            // SUCCESS_OR_DIE(nc_put_att_text(
            //     ncid,
            //     varid[Var_idx::type],
            //     "auxiliary_data",
            //     strlen(att_val.c_str()),
            //     att_val.c_str())
            // );

            // step
            att_val = "simulation step";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::step],
                "long_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "step";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::step],
                "standard_name",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
            att_val = "yes";
            SUCCESS_OR_DIE(nc_put_att_text(
                ncid,
                varid[Var_idx::step],
                "auxiliary_data",
                strlen(att_val.c_str()),
                att_val.c_str())
            );
        }

        // std::cout << "Rank: " << rank << " Attempt to par access\n";
        // SUCCESS_OR_DIE(nc_var_par_access(ncid, NC_GLOBAL, NC_INDEPENDENT));
        for(auto &id: varid)
            SUCCESS_OR_DIE(nc_var_par_access(ncid, id, NC_INDEPENDENT));
        // std::cout << "Rank: " << rank << " Attempt to nc_enddef\n";
        SUCCESS_OR_DIE(nc_enddef(ncid));
        // std::cout << "Rank: " << rank << " after nc_enddef\n";
        if(rank == 0)
        {
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
                    time_steps.data()            // op
                )
            );

            countp[0] = n_trajs_file;
            std::vector<uint64_t> data(n_trajs_file);
            for(uint32_t i=0; i<n_trajs_file; i++)
                data[i] = i;
            SUCCESS_OR_DIE(
                nc_put_vara(
                    ncid,                       // ncid
                    varid[Var_idx::trajectory],       // varid
                    startp.data(),              // startp
                    countp.data(),              // countp
                    data.data()            // op
                )
            );

            countp[0] = num_ens;
            data.resize(num_ens);
            for(uint32_t i=0; i<num_ens; i++)
                data[i] = i;
            SUCCESS_OR_DIE(
                nc_put_vara(
                    ncid,                       // ncid
                    varid[Var_idx::ensemble],       // varid
                    startp.data(),              // startp
                    countp.data(),              // countp
                    data.data()            // op
                )
            );

            countp[0] = num_comp;
            data.resize(num_comp);
            for(uint32_t i=0; i<num_comp; i++)
                data[i] = i;
            SUCCESS_OR_DIE(
                nc_put_vara(
                    ncid,                       // ncid
                    varid[Var_idx::out_param],       // varid
                    startp.data(),              // startp
                    countp.data(),              // countp
                    data.data()            // op
                )
            );
        }
        // std::vector<nc_vlen_t> att_in;
        // SUCCESS_OR_DIE(nc_get_att(in_ncid, NC_GLOBAL, "", att_in));

        // NcFile in_datafile(in_filename, NcFile::read);
        // // global attributes
        // auto attributes = in_datafile.getAtts();
        // for(auto & name_attr: attributes)
        // {
        //     auto attribute = name_attr.second;
        //     NcType type = attribute.getType();
        //     if(type.getName() == "double")
        //     {
        //         std::vector<float> values(1);
        //         attribute.getValues(values.data());
        //         datafile.putAtt(attribute.getName(), type, values[0]);
        //     } else if(type.getName() == "int64"
        //         || type.getName() == "int32" || type.getName() == "int")
        //     {
        //         std::vector<int> values(1);
        //         attribute.getValues(values.data());
        //         datafile.putAtt(attribute.getName(), type, values[0]);
        //     } else if(type.getName() == "char")
        //     {
        //         std::string values;
        //         attribute.getValues(values);
        //         datafile.putAtt(attribute.getName(), values);
        //     }
        // }

        // // column attributes






//         if(rank != 0)
//         {
//             datafile.open(filename + ".nc_wcb", NcFile::write);
//         } else
//         {
//             std::cout << "nc_var_par_access\n";
//             SUCCESS_OR_DIE(nc_var_par_access(ncid, NC_GLOBAL, NC_INDEPENDENT));
//             std::cout << "Return\n";
//             return;
//         }

//         // try
//         // {
//         //     std::cout << "File exists?\n";
//         //     datafile.open(filename + ".nc_wcb", NcFile::newFile);
//         //     std::cout << "No\n";
//         // }
//         // catch(const std::exception& e)
//         // {
//         //     std::cout << "Yes\n";
//         //     SUCCESS_OR_DIE(nc_var_par_access(ncid, NC_GLOBAL, NC_INDEPENDENT));
//         //     std::cout << "Return\n";
//         //     // Apparently the file exists already
//         //     // fill_vectors();
//         //     return;
//         // }
//         // Add dimensions
//         NcDim param_dim = datafile.addDim("Output Parameter", num_comp);
//         NcDim ens_dim = datafile.addDim("ensemble", 1);
//         NcDim traj_dim = datafile.addDim("trajectory", n_trajs_file);
//         NcDim time_dim = datafile.addDim("time"); // unlimited dimension

//         NcVar param_var = datafile.addVar("Output Parameter", ncString, param_dim);
//         NcVar ens_var = datafile.addVar("ensemble", ncInt64, ens_dim);
//         NcVar traj_var = datafile.addVar("trajectory", ncInt64, traj_dim);

//         ens_var.setCompression(
//             enableShuffleFilter, enableDeflateFilter, deflateLevel);
//         traj_var.setCompression(
//             enableShuffleFilter, enableDeflateFilter, deflateLevel);
//         param_var.setCompression(
//             enableShuffleFilter, enableDeflateFilter, deflateLevel);

//         std::vector<uint64_t> startp, countp;
//         startp.push_back(0);
//         countp.push_back(n_trajs_file);
//         std::vector<uint64_t> traj_ids(n_trajs_file);
//         for(uint64_t i=0; i<n_trajs_file; i++)
//         {
//             traj_ids[i] = i;
//         }
//         ens_var.putVar(&cc.ensemble_id);
//         traj_var.putVar(startp, countp, traj_ids.data());


//         startp[0] = 0;
//         countp[0] = 1;
//         for(const auto &p: output_par_idx)
//         {
//             param_var.putVar(startp, countp, &p);
//             startp[0]++;
//         }
//         std::vector<NcDim> dim_vector;
//         dim_vector.push_back(ens_dim);
//         dim_vector.push_back(traj_dim);
//         dim_vector.push_back(time_dim);

//         for(auto &out_par: output_par_idx)
//             var_vector.emplace_back(datafile.addVar(out_par, ncDouble, dim_vector));
//         std::vector<NcDim> dim_vector_grad;

//         dim_vector_grad.push_back(param_dim);
//         dim_vector_grad.push_back(ens_dim);
//         dim_vector_grad.push_back(traj_dim);
//         dim_vector_grad.push_back(time_dim);
//         for(auto &out_grad: output_grad_idx)
//             var_vector.emplace_back(datafile.addVar(out_grad, ncDouble, dim_vector_grad));

//         var_vector.emplace_back(datafile.addVar("time_after_ascent", ncDouble, dim_vector));
//         var_vector.emplace_back(datafile.addVar("time", ncDouble, time_dim));
//         var_vector.emplace_back(datafile.addVar("conv_400", ncByte, dim_vector));
//         var_vector.emplace_back(datafile.addVar("conv_600", ncByte, dim_vector));
//         var_vector.emplace_back(datafile.addVar("slan_400", ncByte, dim_vector));
//         var_vector.emplace_back(datafile.addVar("slan_600", ncByte, dim_vector));
//         var_vector.emplace_back(datafile.addVar("type", ncString, dim_vector));
//         var_vector.emplace_back(datafile.addVar("lat", ncDouble, dim_vector));
//         var_vector.emplace_back(datafile.addVar("lon", ncDouble, dim_vector));
//         var_vector.emplace_back(datafile.addVar("step", ncUint64, dim_vector));

//         for(auto &var: var_vector)
//             var.setCompression(
//                 enableShuffleFilter, enableDeflateFilter, deflateLevel);

//         NcFile in_datafile(in_filename, NcFile::read);

//         // global attributes
//         auto attributes = in_datafile.getAtts();
//         for(auto & name_attr: attributes)
//         {
//             auto attribute = name_attr.second;
//             NcType type = attribute.getType();
//             if(type.getName() == "double")
//             {
//                 std::vector<float> values(1);
//                 attribute.getValues(values.data());
//                 datafile.putAtt(attribute.getName(), type, values[0]);
//             } else if(type.getName() == "int64"
//                 || type.getName() == "int32" || type.getName() == "int")
//             {
//                 std::vector<int> values(1);
//                 attribute.getValues(values.data());
//                 datafile.putAtt(attribute.getName(), type, values[0]);
//             } else if(type.getName() == "char")
//             {
//                 std::string values;
//                 attribute.getValues(values);
//                 datafile.putAtt(attribute.getName(), values);
//             }
//         }
//         // add another global attribute describing where this ensemble
//         // originates from
//         std::string att_name = "Ensemble History";
//         std::string att_val = cc.ens_desc;
//         datafile.putAtt(
//             att_name,
//             att_val);
//         // column attributes
//         // Output Parameter is new. Hence we add it like that:
//         att_name = "long_name";
//         att_val = "gradients are calculated w.r.t. this output parameter";
//         param_var.putAtt(
//             att_name,
//             att_val);
//         att_name = "standard_name";
//         att_val = "output_parameter";
//         param_var.putAtt(
//             att_name,
//             att_val);
//         // Same for step
//         uint32_t offset = num_comp+num_par;
//         NcVar *var = &var_vector[offset+9];
//         att_val = "step";
//         var->putAtt(
//             att_name,
//             att_val);
//         att_name = "long_name";
//         att_val = "simulation step";
//         var->putAtt(
//             att_name,
//             att_val);
//         att_name = "auxiliary_data";
//         att_val = "yes";
//         var->putAtt(att_name, att_val);
//         // and hail N
//         att_name = "long_name";
//         att_val = "specific hail number";
//         var_vector[Nh_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "standard_name";
//         att_val = "specific_number_of_hail_in_air";
//         var_vector[Nh_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "units";
//         att_val = "kg^-1";
//         var_vector[Nh_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "auxiliary_data";
//         att_val = "yes";
//         var_vector[Nh_idx].putAtt(att_name, att_val);
//         // hail out
//         att_name = "long_name";
//         att_val = "sedimentation of specific hail number";
//         var_vector[Nh_out_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "standard_name";
//         att_val = "sedi_outflux_of_hail_number";
//         var_vector[Nh_out_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "units";
//         att_val = "kg^-1 s^-1";
//         var_vector[Nh_out_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "auxiliary_data";
//         att_val = "yes";
//         var_vector[Nh_out_idx].putAtt(att_name, att_val);

//         // hail q
//         att_name = "long_name";
//         att_val = "specific hail content";
//         var_vector[qh_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "standard_name";
//         att_val = "mass_fraction_of_hail_in_air";
//         var_vector[qh_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "units";
//         att_val = "kg kg^-1";
//         var_vector[qh_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "auxiliary_data";
//         att_val = "yes";
//         var_vector[qh_idx].putAtt(att_name, att_val);

//         // hail q out
//         att_name = "long_name";
//         att_val = "sedimentation of hail mixing ratio";
//         var_vector[qh_out_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "standard_name";
//         att_val = "sedi_outflux_of_hail";
//         var_vector[qh_out_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "units";
//         att_val = "kg kg^-1 s^-1";
//         var_vector[qh_out_idx].putAtt(
//             att_name,
//             att_val);
//         att_name = "auxiliary_data";
//         att_val = "yes";
//         var_vector[qh_out_idx].putAtt(att_name, att_val);

//         auto vars = in_datafile.getVars();
//         for(auto & name_var: vars)
//         {
//             auto column_name = name_var.first;
//             auto it = std::find(output_par_idx.begin(), output_par_idx.end(), column_name);
//             if(it != output_par_idx.end())
//             {
//                 uint32_t idx = std::distance(output_par_idx.begin(), it);
//                 var = &var_vector[idx];
//             }
//             else if(column_name == "time_after_ascent")
//                 var = &var_vector[offset];
//             else if(column_name == "time")
//                 var = &var_vector[offset+1];
//             else if(column_name == "conv_400")
//                 var = &var_vector[offset+2];
//             else if(column_name == "conv_600")
//                 var = &var_vector[offset+3];
//             else if(column_name == "slan_400")
//                 var = &var_vector[offset+4];
//             else if(column_name == "slan_600")
//                 var = &var_vector[offset+5];
//             else if(column_name == "type")
//                 var = &var_vector[offset+6];
//             else if(column_name == "lat")
//                 var = &var_vector[offset+7];
//             else if(column_name == "lon")
//                 var = &var_vector[offset+8];
//             else if(column_name == "ensemble")
//                 var = &ens_var;
//             else if(column_name == "trajectory")
//                 var = &traj_var;
//             else if(column_name == "Q_TURBULENCE")
// #if defined MET3D && defined TURBULENCE
//                 var = &var_vector[offset+10];
// #else
//                 continue;
// #endif
//             else
//                 continue;
//             // add auxiliary column to nearly all columns by default
//             if(column_name != "time" && column_name != "time_after_ascent"
//                 && column_name != "lon" && column_name != "lat"
//                 && column_name != "pressure" && column_name != "trajectory"
//                 && column_name != "ensemble")
//             {
//                 att_name = "auxiliary_data";
//                 att_val = "yes";
//                 var->putAtt(att_name, att_val);
//             }
//             if(column_name == "ensemble")
//             {
//                 att_name = "standard_name";
//                 att_val = "ensemble";
//                 var->putAtt(att_name, att_val);
//                 att_name = "long_name";
//                 att_val = "ensemble id that is only consistent within a given history via trajectory id";
//                 var->putAtt(att_name, att_val);
//                 continue;
//             }
//             if(column_name == "trajectory")
//             {
//                 att_name = "standard_name";
//                 att_val = "trajectory";
//                 var->putAtt(att_name, att_val);
//                 att_name = "long_name";
//                 att_val = "trajectory";
//                 var->putAtt(att_name, att_val);
//                 att_name = "description";
//                 att_val = "last number is the id of the instance that ran "
//                     "the simulation and the numbers before are the history";
//                 var->putAtt(att_name, att_val);
//                 continue;
//             }
//             auto attributes = name_var.second.getAtts();
//             for(auto & name_attr: attributes)
//             {
//                 // ie long_name, standard_name
//                 auto attribute = name_attr.second;
//                 if(attribute.getName() == "_FillValue")
//                     continue;
//                 NcType type = attribute.getType();
//                 if(type.getName() == "double" || type.getName() == "float")
//                 {
//                     std::vector<double> values(1);
//                     var->putAtt(
//                         attribute.getName(),
//                         type,
//                         values.size(),
//                         values.data());
//                 } else if(type.getName() == "int64"
//                     || type.getName() == "int32" || type.getName() == "int")
//                 {
//                     std::vector<int> values(1);
//                     var->putAtt(
//                         attribute.getName(),
//                         type,
//                         values.size(),
//                         values.data());
//                 } else if(type.getName() == "uint64"
//                     || type.getName() == "uint32" || type.getName() == "uint")
//                 {
//                     std::vector<uint64_t> values(1);
//                     var->putAtt(
//                         attribute.getName(),
//                         type,
//                         values.size(),
//                         values.data());
//                 } else if(type.getName() == "char")
//                 {
//                     std::string values;
//                     attribute.getValues(values);
//                     var->putAtt(
//                         attribute.getName(),
//                         values);
//                 }
//             }
//         }
//         // all the gradients are auxiliary data
//         for(uint64_t j=0; j<num_par; j++)
//         {
//             att_name = "auxiliary data";
//             att_val = "yes";
//             var_vector[num_comp+j].putAtt(att_name, att_val);
//             // add descriptions to all gradients
//             att_name = "standard_name";
//             att_val = output_grad_idx[j].substr(1);
//             var_vector[num_comp+j].putAtt(att_name, att_val);
//             att_name = "long_name";
//             att_val = output_grad_descr[j];
//             var_vector[num_comp+j].putAtt(att_name, att_val);
//         }
//         in_datafile.close();
//         datafile.close();
        // Make the access independent which is a must due to the dynamic
        // work schedule; This can be expensive though.
        // SUCCESS_OR_DIE(nc_var_par_access(ncid, NC_GLOBAL, NC_INDEPENDENT));
    } else
    {
        // write attributes
#ifdef MET3D
        std::ofstream outfile_att;
        outfile_att.open(filename + "_attributes.txt");
        outfile_att.precision(10);
        if( !outfile_att.is_open() )
        {
            std::cerr << "ERROR while opening the attribute file:\n"
                    << filename << "_attributes.txt\n. Aborting.\n";
            exit(EXIT_FAILURE);
        }
        // Global attributes
        outfile_att << "[Global attributes]\n";
        NcFile datafile(in_filename, NcFile::read);
        auto attributes = datafile.getAtts(); // <string, NcGroupAtt>
        for(auto & name_attr: attributes)
        {
            NcGroupAtt attribute = name_attr.second;
            NcType type = attribute.getType();
            if(type.getName() == "double")
            {
                std::vector<float> values(1);
                attribute.getValues(values.data());
                outfile_att << "name=" << attribute.getName() << "\n"
                            << "type=" << type.getName() << "\n"
                            << "values=" << values[0] << "\n";
            } else if(type.getName() == "int64"
                || type.getName() == "int32" || type.getName() == "int")
            {
                std::vector<int> values(1);
                attribute.getValues(values.data());
                outfile_att << "name=" << attribute.getName() << "\n"
                            << "type=" << type.getName() << "\n"
                            << "values=" << values[0] << "\n";
            } else if(type.getName() == "char")
            {
                std::string values;
                attribute.getValues(values);
                outfile_att << "name=" << attribute.getName() << "\n"
                        << "type=" << type.getName() << "\n"
                        << "values=" << values << "\n";
            }
        }

        // Column attributes
        outfile_att << "[Non global attributes]\n";
        auto vars = datafile.getVars();
        for(auto & name_var: vars)
        {
            auto var = name_var.second;
            auto name = name_var.first;
            outfile_att << "column=" << name << "\n";
            auto attributes = var.getAtts();
            for(auto & name_attr: attributes)
            {
                auto attribute = name_attr.second;
                NcType type = attribute.getType();
                if(type.getName() == "double" || type.getName() == "float")
                {
                    std::vector<double> values(1);
                    attribute.getValues(values.data());
                    outfile_att << attribute.getName() << "=" << values[0] << "\n";
                } else if(type.getName() == "int64"
                    || type.getName() == "int32" || type.getName() == "int")
                {
                    std::vector<int> values(1);
                    attribute.getValues(values.data());
                    outfile_att << attribute.getName() << "=" << values[0] << "\n";
                } else if(type.getName() == "char")
                {
                    std::string values;
                    attribute.getValues(values);
                    outfile_att << attribute.getName() << "=" << values << "\n";
                }
            }
        }
        outfile_att.close();
#endif
        // write reference quantities
        std::ofstream outfile_refs;
        outfile_refs.open(filename + "_reference_values.txt");
        outfile_refs.precision(10);

        if( !outfile_refs.is_open() )
        {
            std::cerr << "ERROR while opening the reference file. Aborting." << std::endl;
            exit(EXIT_FAILURE);
        }

        // Write the reference quantities
        outfile_refs << ref_quant.Tref << " "
            << ref_quant.pref << " "
            << ref_quant.qref << " "
            << ref_quant.Nref << " "
            << ref_quant.wref << " "
            << ref_quant.tref << " "
            << ref_quant.zref << "\n";

        outfile_refs.close();

        // write headers
        std::string suffix = ".txt";
        std::string full_filename;
        full_filename = filename;
        full_filename += suffix;

        outfile.open(full_filename);
        outfile.precision(10);

        if( !outfile.is_open() )
        {
            std::cerr << "ERROR while opening the outputfile. Aborting." << std::endl;
            exit(EXIT_FAILURE);
        }

        // Append the initial values and write headers
        out_tmp << std::setprecision(10) << "step,trajectory,lon,lat,"
#if defined WCB
            << "MAP,";
#endif
#if defined WCB2
            << "WCB_flag,"
            << "dp2h,"
#endif
#if defined WCB2 || defined MET3D
            << "conv_400,"
            << "conv_600,"
            << "slan_400,"
            << "slan_600,";
#endif
#if defined MET3D
        out_tmp
            << "time,"
            << "time_after_ascent,"
            << "type,"
            << "ensemble,"
            << "instance_id,";
#endif
        for(uint32_t i=0; i<output_par_idx.size(); ++i)
            out_tmp << output_par_idx[i]  <<
                ((i < output_par_idx.size()-1) ? "," : "\n");

        std::string basename = "_diff_";
        std::string fname;

        for(int ii = 0; ii < num_comp; ii++)
        {
            fname = filename;
            fname += basename;
            fname += std::to_string(ii);
            fname += suffix;

            out_diff[ii].open(fname);

            if( !out_diff[ii].is_open() )
            {
                std::cerr << "ERROR while opening outputfile. Aborting." << std::endl;
                exit(EXIT_FAILURE);
            }
            out_diff[ii].precision(10);
            out_diff_tmp[ii]
                << std::setprecision(10)
                << "step,"
                << "trajectory,"
                << "Output Parameter,"
                << "lon,"
                << "lat,"
#if defined WCB
                << "MAP,";
#endif
#if defined WCB2
                << "WCB_flag,"
                << "dp2h,"
#endif
#if defined WCB2 || defined MET3D
                << "conv_400,"
                << "conv_600,"
                << "slan_400,"
                << "slan_600,";
#endif
#if defined MET3D
            out_diff_tmp[ii]
                << "time,"
                << "time_after_ascent,"
                << "type,"
                << "ensemble,"
                << "instance_id,";
#endif
            for(uint32_t i=0; i<output_grad_idx.size(); ++i)
                out_diff_tmp[ii] << output_grad_idx[i]  <<
                    ((i < output_grad_idx.size()-1) ? "," : "\n");
        } // End loop over all components
    }
}

void output_handle_t::buffer(const model_constants_t &cc,
    const nc_parameters_t &nc_params,
    const std::vector<codi::RealReverse> &y_single_new,
    const std::vector< std::array<double, num_par > >  &y_diff,
    const uint32_t sub,
    const uint32_t t,
    const double time_new,
    const uint32_t traj_id,
    const uint32_t ensemble,
    const reference_quantities_t &ref_quant)
{
    if(filetype == "netcdf")
    {
        // this->traj = cc.traj_id;
        // this->ens = cc.ensemble_id;

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

        // gradients
        for(uint64_t i=0; i<num_comp; i++) // gradient sensitive to output parameter i
            for(uint64_t j=0; j<num_par; j++) // gradient of input parameter j
                output_buffer[Buffer_idx::n_buffer+j][i + n_snapshots*offset*num_comp] = y_diff[i][j];

        // time after ascent
        output_buffer[Buffer_idx::time_ascent_buf][n_snapshots*offset] =
            nc_params.time_rel + sub*cc.dt;

        // time index
        // output_buffer[num_comp+num_par+1][n_snapshots] =
        //     nc_params.time_abs[t + nc_params.time_idx] + sub*cc.dt;

        // lat
        output_buffer[Buffer_idx::lat_buf][n_snapshots*offset] =
            (nc_params.lat[0] + sub*nc_params.dlat);

        // lon
        output_buffer[Buffer_idx::lon_buf][n_snapshots*offset] =
            (nc_params.lon[0] + sub*nc_params.dlon);

        // flags
        output_buffer_flags[0][n_snapshots*offset] = nc_params.conv_400;
        output_buffer_flags[1][n_snapshots*offset] = nc_params.conv_600;
        output_buffer_flags[2][n_snapshots*offset] = nc_params.slan_400;
        output_buffer_flags[3][n_snapshots*offset] = nc_params.slan_600;

        // type
        // output_buffer_str[0][n_snapshots*offset] = nc_params.type[0];

        // simulation step
        output_buffer_int[0][n_snapshots*offset] = sub + t*cc.num_sub_steps;
    } else
    {
#if defined WCB || defined WCB2
        out_tmp << time_new << "," << traj_id << ","
                << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                << (nc_params.lat[0] + sub*nc_params.dlat) << ","
                << nc_params.ascent_flag << ",";
#elif defined MET3D
        out_tmp << sub + t*cc.num_sub_steps << "," << traj_id << ","
                << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#else
        out_tmp << time_new << "," << traj_id << ","
                << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#endif
#if defined WCB2 || defined MET3D
        out_tmp
#if defined WCB2
                << nc_params.dp2h << ","
#endif
                << nc_params.conv_400 << ","
                << nc_params.conv_600 << ","
                << nc_params.slan_400 << ","
                << nc_params.slan_600 << ",";
#endif
#ifdef MET3D
        out_tmp << nc_params.time_abs[t + nc_params.time_idx] + sub*cc.dt << ","
                << nc_params.time_rel + sub*cc.dt << ","
                << nc_params.type[0] << ","
                << ensemble << ","
                << cc.id << ",";
#endif
        for(int ii = 0 ; ii < num_comp; ii++)
            out_tmp << y_single_new[ii]
                << ((ii == num_comp-1) ? "\n" : ",");

        for(int ii = 0 ; ii < num_comp ; ii++)
        {
#if defined WCB || defined WCB2
            out_diff_tmp[ii] << time_new << "," << traj_id << ","
                            << output_par_idx[ii] << ","
                            << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                            << (nc_params.lat[0] + sub*nc_params.dlat) << ","
                            << nc_params.ascent_flag << ",";
#elif defined MET3D
            out_diff_tmp[ii] << sub + t*cc.num_sub_steps << "," << traj_id << ","
                            << output_par_idx[ii] << ","
                            << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                            << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#else
            out_diff_tmp[ii] << time_new << "," << traj_id << ","
                            << output_par_idx[ii] << ","
                            << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                            << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#endif
#if defined WCB2 || defined MET3D
            out_diff_tmp[ii]
#if defined WCB2
                << nc_params.dp2h << ","
#endif
                << nc_params.conv_400 << ","
                << nc_params.conv_600 << ","
                << nc_params.slan_400 << ","
                << nc_params.slan_600 << ",";
#endif
#if defined MET3D
            out_diff_tmp[ii] << nc_params.time_abs[t + nc_params.time_idx] + sub*cc.dt << ","
                        << nc_params.time_rel + sub*cc.dt << ","
                        << nc_params.type[0] << ","
                        << ensemble << ","
                        << cc.id << ",";
#endif
            for(int jj = 0 ; jj < num_par ; jj++)
                out_diff_tmp[ii] << y_diff[ii][jj]
                    << ((jj==num_par-1) ? "\n" : ",");
        }
    }
    n_snapshots++;
}

void output_handle_t::flush_buffer()
{
    if(filetype == "netcdf")
    {
        std::cout << "Attempt to flush with flushed "
                  << flushed_snapshots
                  << " snaps: " << n_snapshots
                  << " traj: " << traj
                  << " ens: " << ens
                  << " \n";
        // bool got_it = false;
        // do{
        //     try
        //     {
        //         datafile.open(filename + ".nc_wcb", NcFile::write);
        //         got_it = true;
        //     } catch(const std::exception& e2)
        //     {
        //         got_it = false;
        //     }
        // } while(!got_it);

        std::vector<size_t> startp, countp;
        // startp.push_back(flushed_snapshots);
        // countp.push_back(n_snapshots); // number of snapshots so far
        // std::cout << "Attempt to flush time \n";
        // std::cout << "len buffer "
        //           << output_buffer[num_comp+num_par+1].size()
        //           << "\n";
        // std::cout << "start " << startp[0]
        //           << ", count " << countp[0]
        //           << ", in buffer " << output_buffer[num_comp+num_par+1][0]
        //           << ", " << output_buffer[num_comp+num_par+1][1]
        //           << ", " << output_buffer[num_comp+num_par+1][2]
        //           << "\n";
        // time index
        // SUCCESS_OR_DIE(
        //     nc_put_vara_double(
        //         ncid,                       // ncid
        //         varid[Var_idx::time],       // varid
        //         startp.data(),              // startp
        //         countp.data(),              // countp
        //         output_buffer[num_comp+num_par+1].data()    // op
        //     )
        // );


        // var_vector[num_comp+num_par+1].putVar(startp, countp,
        //     output_buffer[num_comp+num_par+1].data());

        // startp.insert(startp.begin(), traj);
        // startp.insert(startp.begin(), ens);
        // countp.insert(countp.begin(), 1);
        // countp.insert(countp.begin(), 1);

        startp.push_back(flushed_snapshots);
        startp.push_back(traj);
        startp.push_back(ens);

        countp.push_back(n_snapshots);
        countp.push_back(1);
        countp.push_back(1);

        // std::cout << "Attempt to flush model state \n";
        // uint64_t offset = Var_idx::out_param+1;
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
        // std::cout << "Attempt to flush time_after_ascent \n";
        // offset += num_comp;
        SUCCESS_OR_DIE(
            nc_put_vara(
                ncid,
                varid[Var_idx::time_ascent],
                startp.data(),
                countp.data(),
                output_buffer[Buffer_idx::time_ascent_buf].data()
            )
        );
        // time after ascent
        // var_vector[offset].putVar(startp, countp,
        //     output_buffer[num_comp+num_par].data());
        // std::cout << "Attempt to flush flags \n";
        // uint64_t offset = Var_idx::conv_400;
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
            // var_vector[offset+i].putVar(startp, countp,
            // output_buffer_flags[i].data());
        // std::cout << "Attempt to flush type \n";
        // offset += output_buffer_flags.size();
        // // type
        // std::vector<uint64_t> startp_str, countp_str;
        // for(auto &p: startp)
        //     startp_str.push_back(p);
        // for(auto &p: countp)
        //     countp_str.push_back(p);
        // countp_str[countp_str.size()-1] = 1;
        // std::cout << "Attempt to flush strings \n";
        // for(uint64_t i=0; i<output_buffer_str.size(); i++)
        // {
        //     // write one string at a time.
        //     for(const auto &t: output_buffer_str[i])
        //     {
        //         std::cout << "Attempt to flush string at " << offset << " + " << i << " \n";
        //         var_vector[offset+i].putVar(
        //             startp_str, countp_str, &t);
        //         startp_str[startp_str.size()-1]++;
        //         if(startp_str[startp_str.size()-1]-startp[startp_str.size()-1] == n_snapshots)
        //             break;
        //     }
        // }

        // offset += output_buffer_str.size();
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
        // var_vector[offset].putVar(startp, countp,
        //     output_buffer[num_comp+num_par+2].data());

        // offset += 1;
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
        // var_vector[offset].putVar(startp, countp,
        //     output_buffer[num_comp+num_par+3].data());

        // offset += 1;
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
        // var_vector[offset].putVar(startp, countp,
        //     output_buffer_int[0].data());

        // offset = Var_idx::n_vars;
        // gradients
        // startp.insert(startp.begin(), 0);
        // countp.insert(countp.begin(), num_comp);
        startp.push_back(0);
        countp.push_back(num_comp);

        std::cout << "Attempt to flush gradients with\n"
                  << "startp: " << startp[0] << ", " << startp[1] << ", " << startp[2] << ", " << startp[3] << "\n"
                  << "countp: " << countp[0] << ", " << countp[1] << ", " << countp[2] << ", " << countp[3] << "\n";

        for(uint64_t j=0; j<num_par; j++)
        {
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
            // var_vector[offset+j].putVar(startp, countp,
            //     output_buffer[num_comp+j].data());
        // datafile.close();
    } else
    {
        outfile << out_tmp.rdbuf();
        for(int ii = 0 ; ii < num_comp ; ii++)
        {
            out_diff[ii] << out_diff_tmp[ii].rdbuf();
            out_diff_tmp[ii].str( std::string() );
            out_diff_tmp[ii].clear();
        }
        out_tmp.str( std::string() );
        out_tmp.clear();
    }
    flushed_snapshots += n_snapshots;
    n_snapshots = 0;
}

void output_handle_t::process_step(
    const model_constants_t &cc,
    const nc_parameters_t &nc_params,
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
        this->buffer(cc, nc_params, y_single_new, y_diff, sub, t,
            time_new, traj_id, ensemble, ref_quant);
    }

    if( (0 == (sub + t*cc.num_sub_steps) % write_index)
        || ( t == cc.num_steps-1 && last_step ) )
    {
        std::cout << "flush buffer\n";
        this->flush_buffer();
        std::cout << "flush buffer done\n";
    }
}