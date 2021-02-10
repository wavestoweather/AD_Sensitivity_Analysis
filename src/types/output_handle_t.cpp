#include "include/types/output_handle_t.h"

output_handle_t::output_handle_t()
{

}

output_handle_t::output_handle_t(
    const std::string filetype,
    const std::string filename,
    const uint64_t n_trajs,
    const uint64_t n_ens,
    const model_constants_t &cc,
    const reference_quantities_t &ref_quant,
    const std::string in_filename,
    const uint32_t write_index,
    const uint32_t snapshot_index)
{
    this->setup(filetype, filename, n_trajs, n_ens, cc, ref_quant,
        in_filename, write_index, snapshot_index);
}

void output_handle_t::setup(
    const std::string filetype,
    const std::string filename,
    const uint64_t n_trajs,
    const uint64_t n_ens,
    const model_constants_t &cc,
    const reference_quantities_t &ref_quant,
    const std::string in_filename,
    const uint32_t write_index,
    const uint32_t snapshot_index)
{
    this->n_trajs = n_trajs;
    this->n_ens = n_ens;
    this->filetype = filetype;
    this->filename = filename;
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
        const uint64_t vec_size = n_trajs * n_ens * total_snapshots; // n_snapshots * num_comp;
        const uint64_t vec_size_grad = num_comp * n_trajs * n_ens * total_snapshots;
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
        for(uint32_t i=0; i<output_buffer_str.size(); i++)
            output_buffer_str[i].resize(vec_size);
        for(uint32_t i=0; i<output_buffer_int.size(); i++)
            output_buffer_int[i].resize(vec_size);
        try
        {
            datafile.open(filename + ".nc_wcb", NcFile::newFile);
        }
        catch(const std::exception& e)
        {
            std::cerr << "Error creating output file:\n"
                        << filename << ".nc_wcb Does the file already exist?"
                        << "\nAborting";
            std::cerr << e.what() << '\n';
            exit(NC_ERR);
        }

            // Add dimensions
        NcDim param_dim = datafile.addDim("Output Parameter", num_comp);
        NcDim ens_dim = datafile.addDim("ensemble", n_ens);
        NcDim traj_dim = datafile.addDim("trajectory", n_trajs);
        NcDim time_dim = datafile.addDim("time"); // unlimited dimension

        NcVar param_var = datafile.addVar("Output Parameter", ncString, param_dim);
        NcVar ens_var = datafile.addVar("ensemble", ncInt64, ens_dim);
        NcVar traj_var = datafile.addVar("trajectory", ncString, traj_dim);

        ens_var.setCompression(
            enableShuffleFilter, enableDeflateFilter, deflateLevel);
        traj_var.setCompression(
            enableShuffleFilter, enableDeflateFilter, deflateLevel);
        param_var.setCompression(
            enableShuffleFilter, enableDeflateFilter, deflateLevel);

        // Add dim data
        // the ensemble number is taken from the filename after idxx_y
        // where y is the ensemble number
        auto sub_str_it = filename.find("_wcb");
        uint32_t ens = std::stoi(filename.substr(sub_str_it-4, sub_str_it-1));
        ens_var.putVar(&ens);
        traj_var.putVar(&cc.id);

        std::vector<uint64_t> startp, countp;
        startp.push_back(0);
        countp.push_back(1);
        for(const auto &p: output_par_idx)
        {
            param_var.putVar(startp, countp, &p);
            startp[0]++;
        }

        dim_vector.push_back(time_dim);
        dim_vector.push_back(traj_dim);
        dim_vector.push_back(ens_dim);

        for(auto &out_par: output_par_idx)
            var_vector.emplace_back(datafile.addVar(out_par, ncDouble, dim_vector));
        std::vector<NcDim> dim_vector_grad;
        dim_vector_grad.push_back(time_dim);
        dim_vector_grad.push_back(traj_dim);
        dim_vector_grad.push_back(ens_dim);
        dim_vector_grad.push_back(param_dim);
        for(auto &out_grad: output_grad_idx)
            var_vector.emplace_back(datafile.addVar(out_grad, ncDouble, dim_vector_grad));

        var_vector.emplace_back(datafile.addVar("time_after_ascent", ncDouble, dim_vector));
        var_vector.emplace_back(datafile.addVar("time", ncDouble, time_dim));
        var_vector.emplace_back(datafile.addVar("conv_400", ncByte, dim_vector));
        var_vector.emplace_back(datafile.addVar("conv_600", ncByte, dim_vector));
        var_vector.emplace_back(datafile.addVar("slan_400", ncByte, dim_vector));
        var_vector.emplace_back(datafile.addVar("slan_600", ncByte, dim_vector));
        var_vector.emplace_back(datafile.addVar("type", ncString, dim_vector));
        var_vector.emplace_back(datafile.addVar("lat", ncDouble, dim_vector));
        var_vector.emplace_back(datafile.addVar("lon", ncDouble, dim_vector));
        var_vector.emplace_back(datafile.addVar("step", ncUint64, dim_vector));

        for(auto &var: var_vector)
            var.setCompression(
                enableShuffleFilter, enableDeflateFilter, deflateLevel);

        NcFile in_datafile(in_filename, NcFile::read);

        // global attributes
        auto attributes = in_datafile.getAtts();
        for(auto & name_attr: attributes)
        {
            auto attribute = name_attr.second;
            NcType type = attribute.getType();
            if(type.getName() == "double")
            {
                std::vector<float> values(1);
                attribute.getValues(values.data());
                datafile.putAtt(attribute.getName(), type, values[0]);
            } else if(type.getName() == "int64"
                || type.getName() == "int32" || type.getName() == "int")
            {
                std::vector<int> values(1);
                attribute.getValues(values.data());
                datafile.putAtt(attribute.getName(), type, values[0]);
            } else if(type.getName() == "char")
            {
                std::string values;
                attribute.getValues(values);
                datafile.putAtt(attribute.getName(), values);
            }
        }
        // column attributes
        // Output Parameter is new. Hence we add it like that:
        std::string att_name = "long_name";
        std::string att_val = "gradients are calculated w.r.t. this output parameter";
        param_var.putAtt(
            att_name,
            att_val);
        att_name = "standard_name";
        att_val = "output_parameter";
        param_var.putAtt(
            att_name,
            att_val);
        // Same for step
        uint32_t offset = num_comp+num_par;
        NcVar *var = &var_vector[offset+9];
        att_val = "step";
        var->putAtt(
            att_name,
            att_val);
        att_name = "long_name";
        att_val = "simulation step";
        var->putAtt(
            att_name,
            att_val);
        att_name = "auxiliary_data";
        att_val = "yes";
        var->putAtt(att_name, att_val);
        // and hail N
        att_name = "long_name";
        att_val = "specific hail number";
        var_vector[Nh_idx].putAtt(
            att_name,
            att_val);
        att_name = "standard_name";
        att_val = "specific_number_of_hail_in_air";
        var_vector[Nh_idx].putAtt(
            att_name,
            att_val);
        att_name = "units";
        att_val = "kg^-1";
        var_vector[Nh_idx].putAtt(
            att_name,
            att_val);
        att_name = "auxiliary_data";
        att_val = "yes";
        var_vector[Nh_idx].putAtt(att_name, att_val);
        // hail out
        att_name = "long_name";
        att_val = "sedimentation of specific hail number";
        var_vector[Nh_out_idx].putAtt(
            att_name,
            att_val);
        att_name = "standard_name";
        att_val = "sedi_outflux_of_hail_number";
        var_vector[Nh_out_idx].putAtt(
            att_name,
            att_val);
        att_name = "units";
        att_val = "kg^-1 s^-1";
        var_vector[Nh_out_idx].putAtt(
            att_name,
            att_val);
        att_name = "auxiliary_data";
        att_val = "yes";
        var_vector[Nh_out_idx].putAtt(att_name, att_val);

        // hail q
        att_name = "long_name";
        att_val = "specific hail content";
        var_vector[qh_idx].putAtt(
            att_name,
            att_val);
        att_name = "standard_name";
        att_val = "mass_fraction_of_hail_in_air";
        var_vector[qh_idx].putAtt(
            att_name,
            att_val);
        att_name = "units";
        att_val = "kg kg^-1";
        var_vector[qh_idx].putAtt(
            att_name,
            att_val);
        att_name = "auxiliary_data";
        att_val = "yes";
        var_vector[qh_idx].putAtt(att_name, att_val);

        // hail q out
        att_name = "long_name";
        att_val = "sedimentation of hail mixing ratio";
        var_vector[qh_out_idx].putAtt(
            att_name,
            att_val);
        att_name = "standard_name";
        att_val = "sedi_outflux_of_hail";
        var_vector[qh_out_idx].putAtt(
            att_name,
            att_val);
        att_name = "units";
        att_val = "kg kg^-1 s^-1";
        var_vector[qh_out_idx].putAtt(
            att_name,
            att_val);
        att_name = "auxiliary_data";
        att_val = "yes";
        var_vector[qh_out_idx].putAtt(att_name, att_val);

        auto vars = in_datafile.getVars();
        for(auto & name_var: vars)
        {
            auto column_name = name_var.first;
            auto it = std::find(output_par_idx.begin(), output_par_idx.end(), column_name);
            if(it != output_par_idx.end())
            {
                uint32_t idx = std::distance(output_par_idx.begin(), it);
                var = &var_vector[idx];
            }
            else if(column_name == "time_after_ascent")
                var = &var_vector[offset];
            else if(column_name == "time")
                var = &var_vector[offset+1];
            else if(column_name == "conv_400")
                var = &var_vector[offset+2];
            else if(column_name == "conv_600")
                var = &var_vector[offset+3];
            else if(column_name == "slan_400")
                var = &var_vector[offset+4];
            else if(column_name == "slan_600")
                var = &var_vector[offset+5];
            else if(column_name == "type")
                var = &var_vector[offset+6];
            else if(column_name == "lat")
                var = &var_vector[offset+7];
            else if(column_name == "lon")
                var = &var_vector[offset+8];
            else if(column_name == "ensemble")
                var = &ens_var;
            else if(column_name == "trajectory")
                var = &traj_var;
            else if(column_name == "Q_TURBULENCE")
#if defined MET3D && defined TURBULENCE
                var = &var_vector[offset+10];
#else
                continue;
#endif
            else
                continue;
            // add auxiliary column to nearly all columns by default
            if(column_name != "time" && column_name != "time_after_ascent"
                && column_name != "lon" && column_name != "lat"
                && column_name != "pressure" && column_name != "trajectory"
                && column_name != "ensemble")
            {
                att_name = "auxiliary_data";
                att_val = "yes";
                var->putAtt(att_name, att_val);
            }
            if(column_name == "ensemble")
            {
                att_name = "standard_name";
                att_val = "ensemble";
                var->putAtt(att_name, att_val);
                att_name = "long_name";
                att_val = "ensemble id that is only consistent within a given history via trajectory id";
                var->putAtt(att_name, att_val);
                continue;
            }
            if(column_name == "trajectory")
            {
                att_name = "standard_name";
                att_val = "trajectory";
                var->putAtt(att_name, att_val);
                att_name = "long_name";
                att_val = "trajectory";
                var->putAtt(att_name, att_val);
                att_name = "description";
                att_val = "last number is the id of the instance that ran "
                    "the simulation and the numbers before are the history";
                var->putAtt(att_name, att_val);
                continue;
            }
            auto attributes = name_var.second.getAtts();
            for(auto & name_attr: attributes)
            {
                // ie long_name, standard_name
                auto attribute = name_attr.second;
                if(attribute.getName() == "_FillValue")
                    continue;
                NcType type = attribute.getType();
                if(type.getName() == "double" || type.getName() == "float")
                {
                    std::vector<double> values(1);
                    var->putAtt(
                        attribute.getName(),
                        type,
                        values.size(),
                        values.data());
                } else if(type.getName() == "int64"
                    || type.getName() == "int32" || type.getName() == "int")
                {
                    std::vector<int> values(1);
                    var->putAtt(
                        attribute.getName(),
                        type,
                        values.size(),
                        values.data());
                } else if(type.getName() == "uint64"
                    || type.getName() == "uint32" || type.getName() == "uint")
                {
                    std::vector<uint64_t> values(1);
                    var->putAtt(
                        attribute.getName(),
                        type,
                        values.size(),
                        values.data());
                } else if(type.getName() == "char")
                {
                    std::string values;
                    attribute.getValues(values);
                    var->putAtt(
                        attribute.getName(),
                        values);
                }
            }
        }
        // all the gradients are auxiliary data
        for(uint64_t j=0; j<num_par; j++)
        {
            att_name = "auxiliary data";
            att_val = "yes";
            var_vector[num_comp+j].putAtt(att_name, att_val);
            // add descriptions to all gradients
            att_name = "standard_name";
            att_val = output_grad_idx[j].substr(1);
            var_vector[num_comp+j].putAtt(att_name, att_val);
            att_name = "long_name";
            att_val = output_grad_descr[j];
            var_vector[num_comp+j].putAtt(att_name, att_val);
        }


        in_datafile.close();
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
        const uint64_t offset = n_trajs*n_ens;

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
                output_buffer[num_comp+j][i + n_snapshots*offset*num_comp] = y_diff[i][j];

        // time after ascent
        output_buffer[num_comp+num_par][n_snapshots*offset] =
            nc_params.time_rel + sub*cc.dt;

        // time index
        output_buffer[num_comp+num_par+1][n_snapshots] =
            nc_params.time_abs[t + nc_params.time_idx] + sub*cc.dt;

        // lat
        output_buffer[num_comp+num_par+2][n_snapshots*offset] =
            (nc_params.lat[0] + sub*nc_params.dlat);

        // lon
        output_buffer[num_comp+num_par+3][n_snapshots*offset] =
            (nc_params.lon[0] + sub*nc_params.dlon);

        // flags
        output_buffer_flags[0][n_snapshots*offset] = nc_params.conv_400;
        output_buffer_flags[1][n_snapshots*offset] = nc_params.conv_600;
        output_buffer_flags[2][n_snapshots*offset] = nc_params.slan_400;
        output_buffer_flags[3][n_snapshots*offset] = nc_params.slan_600;

        // type
        output_buffer_str[0][n_snapshots*offset] = nc_params.type[0];

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
        std::vector<uint64_t> startp, countp;
        startp.push_back(flushed_snapshots);
        countp.push_back(n_snapshots); // number of snapshots so far
        // time index
        var_vector[num_comp+num_par+1].putVar(startp, countp,
            output_buffer[num_comp+num_par+1].data());

        startp.push_back(0);
        startp.push_back(0);
        countp.push_back(n_ens);
        countp.push_back(n_trajs);

        for(uint64_t i=0; i<num_comp; i++)
            var_vector[i].putVar(startp, countp,
                output_buffer[i].data());

        uint64_t offset = num_comp+num_par;
        // time after ascent
        var_vector[offset].putVar(startp, countp,
            output_buffer[num_comp+num_par].data());

        offset += 2;
        // flags
        for(uint64_t i=0; i<output_buffer_flags.size(); i++)
            var_vector[offset+i].putVar(startp, countp,
            output_buffer_flags[i].data());

        offset += output_buffer_flags.size();
        // type
        std::vector<uint64_t> startp_str, countp_str;
        for(auto &p: startp)
            startp_str.push_back(p);
        for(auto &p: countp)
            countp_str.push_back(p);
        countp_str[0] = 1;

        for(uint64_t i=0; i<output_buffer_str.size(); i++)
        {
            // write one string at a time.
            for(const auto &t: output_buffer_str[i])
            {
                var_vector[offset+i].putVar(
                    startp_str, countp_str, &t);
                startp_str[0]++;
                if(startp_str[0]-startp[0] == n_snapshots)
                    break;
            }
        }

        offset += output_buffer_str.size();
        // lat
        var_vector[offset].putVar(startp, countp,
            output_buffer[num_comp+num_par+2].data());

        offset += 1;
        // lon
        var_vector[offset].putVar(startp, countp,
            output_buffer[num_comp+num_par+3].data());

        offset += 1;
        // step
        var_vector[offset].putVar(startp, countp,
            output_buffer_int[0].data());

        offset = num_comp;
        // gradients
        startp.push_back(0);
        countp.push_back(num_comp);

        for(uint64_t j=0; j<num_par; j++)
            var_vector[offset+j].putVar(startp, countp,
                output_buffer[num_comp+j].data());
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
        this->flush_buffer();
    }
}