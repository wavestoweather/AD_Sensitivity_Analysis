#include "include/types/netcdf_reader_t.h"

netcdf_reader_t::netcdf_reader_t(
    const uint32_t &buffer_size)
{
    dimid.resize(Dim_idx::n_dims);
    varid.resize(Par_idx::n_pars);
    this->n_timesteps_buffer = buffer_size;
    std::cout << "Got buffer_size " << buffer_size
              << ", n_t_buffer: " << this->n_timesteps_buffer << "\n";
    this->time_buffer_idx = 0;
    for(auto &b: buffer)
        b.resize(this->n_timesteps_buffer);
}


void netcdf_reader_t::load_vars()
{
    varid.resize(Par_idx::n_pars);
#ifdef WCB2
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "latitude",
            &varid[Par_idx::lat]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "longitude",
            &varid[Par_idx::lon]
        )
    );
#else
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "lat",
            &varid[Par_idx::lat]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "lon",
            &varid[Par_idx::lon]
        )
    );
#endif
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "z",
            &varid[Par_idx::height]
        )
    );
#if defined WCB || defined WCB2
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "z",
            &varid[Par_idx::height]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "P",
            &varid[Par_idx::pressure]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "T",
            &varid[Par_idx::temperature]
        )
    );
#elif defined MET3D
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "time_after_ascent",
            &varid[Par_idx::time_after_ascent]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "pressure",
            &varid[Par_idx::pressure]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "T",
            &varid[Par_idx::temperature]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "w",
            &varid[Par_idx::ascent]
        )
    );
#else
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "time_rel",
            &varid[Par_idx::time_after_ascent]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "p",
            &varid[Par_idx::pressure]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "t",
            &varid[Par_idx::temperature]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "w",
            &varid[Par_idx::ascent]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QIin",
            &varid[Par_idx::qi_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QSin",
            &varid[Par_idx::qs_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QRin",
            &varid[Par_idx::qr_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QGin",
            &varid[Par_idx::qg_in]
        )
    );
#endif
#if defined WCB2 || defined MET3D
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QI_IN",
            &varid[Par_idx::qi_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QS_IN",
            &varid[Par_idx::qs_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QR_IN",
            &varid[Par_idx::qr_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QG_IN",
            &varid[Par_idx::qg_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "NI_IN",
            &varid[Par_idx::qi_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "NS_IN",
            &varid[Par_idx::qs_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "NR_IN",
            &varid[Par_idx::qr_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "NG_IN",
            &varid[Par_idx::qg_in]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "conv_400",
            &varid[Par_idx::conv_400]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "conv_600",
            &varid[Par_idx::conv_600]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "slan_400",
            &varid[Par_idx::slan_400]
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "slan_600",
            &varid[Par_idx::slan_600]
        )
    );
#endif
#if defined MET3D && defined TURBULENCE
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "Q_TURBULENCE",
            &varid[Par_idx::q_turb]
        )
    );
#endif
}

void netcdf_reader_t::buffer_params(
    const reference_quantities_t &ref_quant)
{
    // check amount of timesteps left to read
    if(time_idx + n_timesteps_buffer > n_timesteps_in)
    {
        n_timesteps_buffer = n_timesteps_in - time_idx;
    }
    // Check for n_timesteps_in == 0?
    std::vector<size_t> startp, countp;
#if defined WCB || defined WCB2
    startp.push_back(time_idx);
    startp.push_back(traj_idx);
    countp.push_back(n_timesteps_buffer);
    countp.push_back(1);
#elif defined MET3D
    startp.push_back(ens_idx);
    startp.push_back(traj_idx);
    startp.push_back(time_idx);

    countp.push_back(1);
    countp.push_back(1);
    countp.push_back(n_timesteps_buffer);
#else
    startp.push_back(traj_idx);
    startp.push_back(time_idx);
    countp.push_back(1);
    countp.push_back(n_timesteps_buffer);
#endif
    // Error here
    std::cout << "startp: " << startp[0]
              << ", " << startp[1]
              << ", " << startp[2]
              << "\n";
    std::cout << "countp: " << countp[0]
              << ", " << countp[1]
              << ", " << countp[2]
              << "\n";
    for(int i=0; i<Par_idx::n_pars; i++)
    {
        SUCCESS_OR_DIE(
            nc_get_vara(
                ncid,
                varid[i],
                startp.data(),
                countp.data(),
                buffer[i].data()
            )
        );
    }

    for(auto &v: buffer[Par_idx::pressure])
    {
        v /= ref_quant.pref;
#if !defined WCB2 && !defined MET3D
        v *= 100;
#endif
    }
    for(auto &v: buffer[Par_idx::temperature])
        v /= ref_quant.Tref;
#if defined(FLUX) && !defined(WCB)
    for(auto &v: buffer[Par_idx::qr_in])
        v /= ref_quant.qref;
    for(auto &v: buffer[Par_idx::nr_in])
        v /= ref_quant.Nref;
  #if defined(RK4ICE)
    for(auto &v: buffer[Par_idx::qi_in])
        v /= ref_quant.qref;
    for(auto &v: buffer[Par_idx::qs_in])
        v /= ref_quant.qref;
    for(auto &v: buffer[Par_idx::qg_in])
        v /= ref_quant.qref;
    for(auto &v: buffer[Par_idx::ni_in])
        v /= ref_quant.Nref;
    for(auto &v: buffer[Par_idx::ns_in])
        v /= ref_quant.Nref;
    for(auto &v: buffer[Par_idx::ng_in])
        v /= ref_quant.Nref;
  #endif
#endif
    for(auto &v: buffer[Par_idx::ascent])
        v /= ref_quant.wref;
#if defined MET3D && defined TURBULENCE
    for(auto &v: buffer[Par_idx::q_turb])
        v /= ref_quant.qref;
#endif
    time_buffer_idx = 0;
}


void netcdf_reader_t::read_buffer(
    model_constants_t &cc,
    const reference_quantities_t &ref_quant,
    std::vector<codi::RealReverse> &y_single_old,
    std::vector<codi::RealReverse> &inflows,
    const uint32_t &step,
    const bool &checkpoint_flag,
    const bool &start_over_env)
{
    time_buffer_idx += step+start_time_idx-time_idx;
    time_idx = step+start_time_idx;
    // check if buffer needs to be reloaded
    if(time_buffer_idx >= n_timesteps_buffer-1)
    {
        buffer_params(ref_quant);
    }
    // Reset outflow
    y_single_old[qi_out_idx] = 0;
    y_single_old[qs_out_idx] = 0;
    y_single_old[qr_out_idx] = 0;
    y_single_old[qg_out_idx] = 0;
    y_single_old[qh_out_idx] = 0;
    y_single_old[Ni_out_idx] = 0;
    y_single_old[Ns_out_idx] = 0;
    y_single_old[Nr_out_idx] = 0;
    y_single_old[Ng_out_idx] = 0;
    y_single_old[Nh_out_idx] = 0;

    // Set values from a given trajectory
    if((step==0 && !checkpoint_flag) || start_over_env)
    {
        y_single_old[p_idx] = buffer[Par_idx::pressure][time_buffer_idx];
        y_single_old[T_idx] = buffer[Par_idx::temperature][time_buffer_idx];
#if defined(RK4ICE) || defined(RK4NOICE)
        y_single_old[w_idx] = buffer[Par_idx::ascent][time_buffer_idx];
        if(time_idx == n_timesteps_in)
        {
            cc.constants[static_cast<int>(Cons_idx::dw)] = 0;
        } else
        {
            cc.constants[static_cast<int>(Cons_idx::dw)] = (
                (buffer[Par_idx::ascent][time_buffer_idx+1] - buffer[Par_idx::ascent][time_buffer_idx])
                / (cc.dt*cc.num_sub_steps)
            );
        }
#endif
        y_single_old[z_idx] = buffer[Par_idx::height][0];
        // set inflow
#if defined(FLUX) && !defined(WCB)
        inflows[qr_in_idx] =  buffer[Par_idx::qr_in][time_buffer_idx];
        inflows[Nr_in_idx] =  buffer[Par_idx::nr_in][time_buffer_idx];
    #if defined(RK4ICE)
        inflows[qi_in_idx] =  buffer[Par_idx::qi_in][time_buffer_idx];
        inflows[Ni_in_idx] =  buffer[Par_idx::ni_in][time_buffer_idx];
        inflows[qs_in_idx] =  buffer[Par_idx::qs_in][time_buffer_idx];
        inflows[Ns_in_idx] =  buffer[Par_idx::ns_in][time_buffer_idx];
        inflows[qg_in_idx] =  buffer[Par_idx::qg_in][time_buffer_idx];
        inflows[Ng_in_idx] =  buffer[Par_idx::ng_in][time_buffer_idx];
    #endif
#else
        inflows[qr_in_idx] = 0;
        inflows[nr_in_idx] = 0;
    #if defined(RK4ICE)
        inflows[qi_in_idx] = 0;
        inflows[qs_in_idx] = 0;
        inflows[qg_in_idx] = 0;
        inflows[ni_in_idx] = 0;
        inflows[ns_in_idx] = 0;
        inflows[ng_in_idx] = 0;
    #endif
#endif
#if defined MET3D && defined TURBULENCE
        inflows[qv_in_idx] = buffer[Par_idx::turbulence][time_buffer_idx];
#endif
    }

    if((step==0 && !checkpoint_flag))
    {
        std::vector<size_t> startp;
#if defined WCB || defined WCB2
        startp.push_back(time_idx);
        startp.push_back(traj_idx);
#elif defined MET3D
        startp.push_back(ens_idx);
        startp.push_back(traj_idx);
        startp.push_back(time_idx);
#else
        startp.push_back(traj_idx);
        startp.push_back(time_idx);
#endif
        int tmp_var_id;
        auto load_var1 = [&](
            const char *name,
            auto target)
        {
            SUCCESS_OR_DIE(
                nc_inq_varid(
                    ncid,
                    name,
                    &tmp_var_id
                )
            );
            SUCCESS_OR_DIE(
                nc_get_var1(
                    ncid,
                    tmp_var_id,
                    startp.data(),
                    target
                )
            );
        };
        load_var1(
#if defined WCB || defined WCB2 || defined MET3D
            "QC",
#else
            "qc",
#endif
            &(y_single_old[qc_idx])
        );
#if defined(RK4ICE)
        load_var1(
    #if defined WCB || defined WCB2 || defined MET3D
            "QI",
    #else
            "qi",
    #endif
            &(y_single_old[qi_idx])
        );
        load_var1(
    #if defined WCB || defined WCB2 || defined MET3D
            "QS",
    #else
            "qs",
    #endif
            &(y_single_old[qs_idx])
        );
#endif
#if !defined(WCB) && defined(RK4ICE)
        load_var1(
    #if defined WCB || defined WCB2 || defined MET3D
            "QG",
    #else
            "qg",
    #endif
            &(y_single_old[qg_idx])
        );
#elif defined(RK4ICE)
        y_single_old[qg_idx] = 0;
#endif
#if defined(RK4ICE)
        y_single_old[qh_idx] = 0.0; // qh. We don't have hail in the trajectoris
        y_single_old[Nh_idx] = 0.0; // Nh. We don't have hail in the trajectoris
#endif

#if defined(RK4ICE) && (defined(WCB2) || defined(MET3D))
        load_var1(
            "NCGRAUPEL",
            &(y_single_old[Ng_idx])
        );
        load_var1(
            "NCICE",
            &(y_single_old[Ni_idx])
        );
        load_var1(
            "NCSNOW",
            &(y_single_old[Ns_idx])
        );
#endif
#if defined WCB2 || defined MET3D
        y_single_old[qr_out_idx] = 0;
        y_single_old[Nr_out_idx] = 0;
        load_var1(
            "NCCLOUD",
            &(y_single_old[Nc_idx])
        );
        load_var1(
            "NCRAIN",
            &(y_single_old[Nr_idx])
        );
#else
        codi::RealReverse denom = 0;
        denom = (get_at(cc.cloud.max_x - cc.cloud.constants, Particle_cons_idx::min_x)) / 2.0 + get_at(cc.cloud.constants, Particle_cons_idx::min_x);
        y_single_old[Nc_idx] = y_single_old[qc_idx] * ref_quant.qref / (denom); //*10e2);  // Nc
        denom = (get_at(cc.rain.max_x - cc.rain.constants, Particle_cons_idx::min_x)) / 2 + get_at(cc.rain.constants, Particle_cons_idx::min_x);
        y_single_old[Nr_idx] = y_single_old[qr_idx] * ref_quant.qref / (denom); //*10e2);  // Nr
        denom = get_at(cc.cloud.constants, Particle_cons_idx::min_x) / 2.0;
    #if defined(RK4ICE)
        denom = (get_at(cc.ice.max_x - cc.ice.constants, Particle_cons_idx::min_x)) / 2.0 + get_at(cc.ice.constants, Particle_cons_idx::min_x);
        y_single_old[Ni_idx] = y_single_old[qi_idx] * ref_quant.qref / (denom); //*10e2); // Ni
        denom = (get_at(cc.snow.max_x - cc.snow.constants, Particle_cons_idx::min_x)) / 2.0 + get_at(cc.snow.constants, Particle_cons_idx::min_x);
        y_single_old[Ns_idx] = y_single_old[qs_idx] * ref_quant.qref / (denom); //*10e2); // Ns
        denom = (get_at(cc.graupel.max_x - cc.graupel.constants, Particle_cons_idx::min_x)) / 2.0 + get_at(cc.graupel.constants, Particle_cons_idx::min_x);
        y_single_old[Ng_idx] = y_single_old[qg_idx] * ref_quant.qref / (denom); //*10e2); // Ng
    #endif
#endif
        // Actually only needed for single moment schemes
        cc.constants[static_cast<int>(Cons_idx::Nc_prime)] = y_single_old[Nc_idx];
    }
}


void netcdf_reader_t::init_netcdf(
#ifdef MET3D
    double &start_time,
#endif
    const char *input_file,
    const bool &checkpoint_flag,
    model_constants_t &cc,
    const int &simulation_mode,
    const double current_time)
{
    time_buffer_idx = 0;
    traj_idx = 0;//cc.traj_id; // either the same via input.traj for all or static schedule without perturbance
    ens_idx = 0;//cc.ensemble_id; // check cc name
    // num_sub_steps = cc.num_sub_steps;

    if(!already_open)
    {
        if((simulation_mode == grid_sensitivity) || (simulation_mode == trajectory_sensitivity))
        {
            SUCCESS_OR_DIE(
                nc_open_par(
                    input_file,
                    NC_NOWRITE,
                    MPI_COMM_WORLD,
                    MPI_INFO_NULL,
                    &ncid
                )
            );
        } else
        {
            SUCCESS_OR_DIE(
                nc_open(
                    input_file,
                    NC_NOWRITE,
                    &ncid
                )
            );
        }

        SUCCESS_OR_DIE(
            nc_inq_dimid(
                ncid,
#ifdef WCB
                "ntra",
#elif defined MET3D
                "trajectory",
#else
                "id",
#endif
                &dimid[Dim_idx::trajectory_dim_idx]
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_dimid(
                ncid,
                "ensemble",
                &dimid[Dim_idx::ensemble_dim_idx]
            )
        );
    }
    SUCCESS_OR_DIE(
        nc_inq_dimlen(
            ncid,
            dimid[Dim_idx::trajectory_dim_idx],
            &n_trajectories
        )
    );
    SUCCESS_OR_DIE(
        nc_inq_dimlen(
            ncid,
            dimid[Dim_idx::ensemble_dim_idx],
            &n_ensembles
        )
    );
    if(n_trajectories <= traj_idx)
    {
        std::cerr << "Number of trajectories in netCDF file: " << n_trajectories << "\n";
        std::cerr << "You asked for trajectory with index " << traj_idx
                    << " which does not exist. ABORTING.\n";
        SUCCESS_OR_DIE(NC_TRAJ_IDX_ERR);
    }
    if(!already_open)
    {
        SUCCESS_OR_DIE(
            nc_inq_dimid(
                ncid,
#ifdef WCB
                "ntim",
#else
                "time",
#endif
                &dimid[Dim_idx::time_dim_idx]
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_dimlen(
                ncid,
                dimid[Dim_idx::time_dim_idx],
                &n_timesteps_in
            )
        );

        load_vars();
        already_open = true;
    }
    std::vector<size_t> startp, countp;
#ifdef MET3D
    countp.push_back(1);
    countp.push_back(1);
    countp.push_back(1);

    startp.push_back(ens_idx); // ensemble
    startp.push_back(traj_idx); // trajectory
    startp.push_back(0); // time
    uint64_t start_time_idx = 0;
    int time_ascent_id;
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "time_after_ascent",
            &time_ascent_id
        )
    );
    double rel_start_time;
    SUCCESS_OR_DIE(
        nc_get_var1(
            ncid,
            time_ascent_id,
            startp.data(),
            &rel_start_time
        )
    );
    if(!std::isnan(start_time) && !checkpoint_flag)
    {
        // Calculate the needed index
        start_time_idx = (start_time-rel_start_time)/cc.dt_traject;
    } else if(checkpoint_flag && !std::isnan(start_time))
    {
        // Calculate the needed index
        start_time_idx = ceil(start_time-rel_start_time + current_time)/cc.dt_traject;
    }
    // n_timesteps_in -= start_time_idx;
    time_idx = start_time_idx;
    this->start_time_idx = start_time_idx;
#else
    time_idx = 1;
    this->start_time_idx = 1;
#endif
}


void netcdf_reader_t::read_initial_values(
    std::vector<double> &y_init,
    reference_quantities_t &ref_quant,
    model_constants_t &cc,
    bool const &checkpoint_flag)
{
    traj_idx = cc.traj_id;
    ens_idx = cc.ensemble_id;
    buffer_params(ref_quant);

    std::vector<size_t> startp, countp;

    if(!checkpoint_flag)
    {
#if defined MET3D
        startp.push_back(ens_idx);
        startp.push_back(traj_idx);
        startp.push_back(time_idx);
        countp.push_back(0);
        countp.push_back(0);
        countp.push_back(0);
#elif defined WCB || WCB2
        startp.push_back(time_idx);
        startp.push_back(traj_idx);
        countp.push_back(1);
        countp.push_back(1);
#else
        startp.push_back(traj_idx);
        startp.push_back(time_idx);
        countp.push_back(1);
        countp.push_back(1);
#endif
        y_init[T_idx] = buffer[Par_idx::temperature][0];
        y_init[p_idx] = buffer[Par_idx::pressure][0];

        int tmp_id;
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
#if defined WCB || defined WCB2 || MET3D
                "QV",
#else
                "qv",
#endif
                &tmp_id
            )
        );
        // std::cout << "startp: " << startp[0] << ", " << startp[1] << ", " << startp[2] << "\n";
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[qv_idx])
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
#if defined WCB || defined WCB2 || MET3D
                "QC",
#else
                "qc",
#endif
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[qc_idx])
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
#if defined WCB || defined WCB2 || MET3D
                "QR",
#else
                "qr",
#endif
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[qr_idx])
            )
        );
#if defined(RK4ICE)
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
    #if defined WCB || defined WCB2 || MET3D
                "QI",
    #else
                "qi",
    #endif
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[qi_idx])
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
    #if defined WCB || defined WCB2 || MET3D
                "QS",
    #else
                "qs",
    #endif
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[qs_idx])
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
    #if defined WCB || defined WCB2 || MET3D
                "QG",
    #else
                "qg",
    #endif
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[qg_idx])
            )
        );
        y_init[qh_idx] = 0.0; // hail that is not in the trajectory
        y_init[qh_out_idx] = 0.0;
        y_init[Nh_out_idx] = 0.0;
#endif
#ifdef WCB
        y_init[w_idx] = 0;
#else
        y_init[w_idx] = buffer[Par_idx::ascent][0];
#endif

#if defined WCB2 || defined MET3D
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "NCRAIN",
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[Nr_idx])
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "NCCLOUD",
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[Nc_idx])
            )
        );
    #if defined(RK4ICE)
        // We can read the sedimentation from the original file
        // But we actually want to set them to zero and calculate them
        y_init[qr_out_idx] = 0;
        y_init[Nr_out_idx] = 0;
        y_init[qi_out_idx] = 0;
        y_init[qs_out_idx] = 0;
        y_init[qg_out_idx] = 0;
        y_init[qh_out_idx] = 0;
        y_init[Nh_out_idx] = 0;
        y_init[Ni_out_idx] = 0;
        y_init[Ns_out_idx] = 0;
        y_init[Ng_out_idx] = 0;

        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "NCICE",
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[Ni_idx])
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "NCSNOW",
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[Ns_idx])
            )
        );
        SUCCESS_OR_DIE(
            nc_inq_varid(
                ncid,
                "NCGRAUPEL",
                &tmp_id
            )
        );
        SUCCESS_OR_DIE(
            nc_get_var1(
                ncid,
                tmp_id,
                startp.data(),
                &(y_init[Ng_idx])
            )
        );
    #endif
#else
        y_init[qr_out_idx] = 0.0;
        y_init[Nr_out_idx] = 0;
        y_init[Nc_idx] = 0;
        y_init[Nr_idx] = 0;
    #if defined(RK4ICE)
        // We initialize the sedimentation with 0 for the stepper
        y_init[qi_out_idx] = 0.0;
        y_init[qs_out_idx] = 0.0;
        y_init[qg_out_idx] = 0.0;
        y_init[qh_out_idx] = 0;
        y_init[Nh_out_idx] = 0;
        y_init[Ni_out_idx] = 0;
        y_init[Ns_out_idx] = 0;
        y_init[Ng_out_idx] = 0;
        y_init[Ni_idx] = 0;
        y_init[Ns_idx] = 0;
        y_init[Ng_idx] = 0;
    #endif
#endif
        y_init[qv_idx]      /= ref_quant.qref;
        y_init[qc_idx]      /= ref_quant.qref;
        y_init[qr_idx]      /= ref_quant.qref;
        y_init[qi_idx]      /= ref_quant.qref;
        y_init[qs_idx]      /= ref_quant.qref;
        y_init[qg_idx]      /= ref_quant.qref;
        y_init[qh_idx]      /= ref_quant.qref;

        y_init[Nc_idx]      /= ref_quant.Nref;
        y_init[Nr_idx]      /= ref_quant.Nref;
        y_init[Ni_idx]      /= ref_quant.Nref;
        y_init[Ns_idx]      /= ref_quant.Nref;
        y_init[Ng_idx]      /= ref_quant.Nref;
        y_init[Nh_idx]      /= ref_quant.Nref;
        y_init[z_idx] = buffer[Par_idx::height][0];

#if defined(RK4ICE) || defined(RK4NOICE) || defined(MET3D)
        double dw = buffer[Par_idx::ascent][1] - buffer[Par_idx::ascent][0];
        cc.constants[static_cast<int>(Cons_idx::dw)] = dw / (cc.dt*cc.num_sub_steps);
        y_init[n_inact_idx] = 0;
        y_init[depo_idx] = 0;
        y_init[sub_idx] = 0;
#endif
#ifdef MET3D
    // Get the time coordinates
    startp[3] = time_idx;
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "time",
            &tmp_id
        )
    );
    SUCCESS_OR_DIE(
        nc_get_var1(
            ncid,
            tmp_id,
            startp.data(),
            &(cc.start_time)
        )
    );
#else
    cc.start_time = 0;
#endif
    }
}