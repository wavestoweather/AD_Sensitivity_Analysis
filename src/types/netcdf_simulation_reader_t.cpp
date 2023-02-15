#include "include/types/netcdf_simulation_reader_t.h"

netcdf_simulation_reader_t::netcdf_simulation_reader_t(
    const uint32_t &buffer_size) {

    dimid.resize(Dim_idx::n_dims);
    varid.resize(Par_idx::n_pars);
    sensvarid.resize(Sens_par_idx::n_sens_pars);
    this->n_timesteps_buffer = buffer_size;
    this->time_buffer_idx = 0;
    this->n_traj_buffer = 600;
    n_readable_timesteps.resize(n_traj_buffer);
    relative_time_buffer.resize(10 * this->n_traj_buffer);
    mem_usage = 0;
    for (uint32_t i=0; i < buffer.size(); i++) {
        this->buffer[i].resize(this->n_timesteps_buffer * this->n_traj_buffer);
        mem_usage += (this->buffer[i].size() * sizeof(double))/1024;
    }
    for (uint32_t i=0; i < buffer_sens.size(); i++) {
        this->buffer_sens[i].resize(this->n_timesteps_buffer * this->n_traj_buffer);
        mem_usage += (this->buffer_sens[i].size() * sizeof(double))/1024;
    }
    mem_usage += (n_readable_timesteps.size() * sizeof(uint64_t))/1024;
    mem_usage += (relative_time_buffer.size() * sizeof(double))/1024;
    already_open = false;
    std::cout << "Reading data allocated " << mem_usage/1024 << " MBytes as buffer\n";
}


void netcdf_simulation_reader_t::load_vars(
    const char* input_file) {

    SUCCESS_OR_DIE(
        nc_open(
            input_file,
            NC_NOWRITE,
            &ncid));
    SUCCESS_OR_DIE(
        nc_inq_dimid(
            ncid,
            "trajectory",
            &dimid[Dim_idx::trajectory_dim_idx]));
    SUCCESS_OR_DIE(
        nc_inq_dimid(
            ncid,
            "ensemble",
            &dimid[Dim_idx::ensemble_dim_idx]));
    SUCCESS_OR_DIE(
        nc_inq_dimid(
            ncid,
            "time",
            &dimid[Dim_idx::time_dim_idx]));
    SUCCESS_OR_DIE(
        nc_inq_dimid(
            ncid,
            "Output_Parameter_ID",
            &dimid[Dim_idx::output_para_idx]));
    for (uint32_t i = 0; i < Par_idx::n_pars; i++) {
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                order_par[i].c_str(),
                &varid[i]));
    }
    for (uint32_t i = 0; i < Sens_par_idx::n_sens_pars; i++) {
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                order_sens[i].c_str(),
                &sensvarid[i]));
    }
    SUCCESS_OR_DIE(
        nc_inq_dimlen(
            ncid,
            dimid[Dim_idx::time_dim_idx],
            &n_timesteps_in));
    SUCCESS_OR_DIE(
        nc_inq_dimlen(
            ncid,
            dimid[Dim_idx::trajectory_dim_idx],
            &n_trajectories));
    const size_t indexp = 0;
    SUCCESS_OR_DIE(
        nc_get_var1_double(
            ncid,
            dimid[Dim_idx::time_dim_idx],
            &indexp,
            &init_time));
}


void netcdf_simulation_reader_t::init_netcdf(
    const uint32_t &traj_idx) {

    this->traj_idx = traj_idx;
    if (traj_idx%n_traj_buffer != 0) return;

    std::vector<size_t> startp, countp;
    startp.push_back(0);            // ensemble
    startp.push_back(traj_idx);     // trajectory
    startp.push_back(0);            // time
    countp.push_back(1);
    countp.push_back(n_traj_buffer);
    countp.push_back(10);

    std::vector<size_t> startp2, countp2;
    startp2.push_back(0);               // Out index
    startp2.push_back(0);               // ensemble
    startp2.push_back(traj_idx);        // trajectory
    startp2.push_back(0);  // time
    countp2.push_back(1);
    countp2.push_back(1);
    countp2.push_back(n_traj_buffer);
    countp2.push_back(n_timesteps_buffer);

    if (startp[1] + countp[1] >= n_trajectories) {
        countp[1] = n_trajectories - startp[1];
        countp2[2] = countp[1];
    }
    SUCCESS_OR_DIE(
        nc_get_vara_double(
            ncid,
            varid[Par_idx::time_after_ascent],
            startp.data(),
            countp.data(),
            relative_time_buffer.data()));
    countp[2] = n_timesteps_buffer;

    if (startp[2] + countp[2] >= n_timesteps_in) {
        countp[2] = n_timesteps_in - startp[2];
        countp2[3] = countp[2];
        read_time_buffer = countp[2];
    }
//    for (auto traj=traj_idx; traj < traj_idx+n_traj_buffer; traj++) {
//        if (traj >= n_trajectories) break;
        // get the relative time
//        startp[1] = traj;
//        countp[2] = 10;
//        std::vector<double> rel_start_time(10);
//        SUCCESS_OR_DIE(
//                nc_get_vara_double(
//                        ncid,
//                        varid[Par_idx::time_after_ascent],
//                        startp.data(),
//                        countp.data(),
//                        rel_start_time.data()));
//        for (uint32_t i = 0; i < rel_start_time.size()-1; i++) {
//            if (std::isnan(rel_start_time[i])) continue;
//            if (start_time < rel_start_time[i]) {
//                start_time_idx = 0;
//            } else {
//                start_time_idx = i + (start_time - rel_start_time[i]) / DELTA_TIMESTEP;
//            }
//        }
//        startp[2] = start_time_idx;
//        if (start_time_idx + n_timesteps_buffer >= n_timesteps_in) {
//            countp[2] = n_timesteps_in-start_time_idx-1;
//        } else {
//            countp[2] = n_timesteps_buffer;
//        }
//        n_readable_timesteps[traj%n_traj_buffer] = countp[2];
//        if (n_readable_timesteps[traj%n_traj_buffer] > n_timesteps_buffer) {
//            std::cout << "The amount of readable time steps is larger than "
//                      << "the buffer size. You should increase the buffer size.\n";
//            countp[2] = n_timesteps_buffer;
//        }
//        std::cout << "buffer time done in " << dt_total << " s\n";
        // Load all the data into buffer
//        for (uint32_t i = 0; i < Par_idx::n_pars; i++) {
//            SUCCESS_OR_DIE(
//                    nc_get_vara_double(
//                            ncid,
//                            varid[i],
//                            startp.data(),
//                            countp.data(),
//                            buffer[i].data() + (traj%n_traj_buffer)*n_timesteps_buffer));
//        }
//        auto mod_time = std::chrono::system_clock::now();
//        dt_total = ((std::chrono::duration<double>)(mod_time - buffer_time)).count();
//        std::cout << "mod done in " << dt_total << " s\n";
//        startp2[2] = traj;
//        startp2[3] = start_time_idx;
//        countp2[3] = countp[2];
//
//        for (uint32_t o_id = 0; o_id < 3; o_id++) {
//            startp2[0] = o_id;
//            for (uint32_t i = 0; i < Sens_par_idx::n_sens_pars; i++) {
//                SUCCESS_OR_DIE(
//                        nc_get_vara_double(
//                                ncid,
//                                sensvarid[i],
//                                startp2.data(),
//                                countp2.data(),
//                                buffer_sens[i*3 + o_id].data() + (traj%n_traj_buffer)*n_timesteps_buffer));
//            }
//        }

    for (uint32_t i = 0; i < Par_idx::n_pars; i++) {
        SUCCESS_OR_DIE(
            nc_get_vara_double(
                ncid,
                varid[i],
                startp.data(),
                countp.data(),
                buffer[i].data()));
    }
//    startp2[2] = traj_idx;
//    startp2[3] = start_time_idx;
//    countp2[3] = countp[2];

    for (uint32_t o_id = 0; o_id < 3; o_id++) {
        startp2[0] = o_id;
        for (uint32_t i = 0; i < Sens_par_idx::n_sens_pars; i++) {
            SUCCESS_OR_DIE(
                    nc_get_vara_double(
                            ncid,
                            sensvarid[i],
                            startp2.data(),
                            countp2.data(),
                            buffer_sens[i*3 + o_id].data()));
        }
    }
//    }
}

void netcdf_simulation_reader_t::close_netcdf() {
    SUCCESS_OR_DIE(ncclose(ncid));
}
