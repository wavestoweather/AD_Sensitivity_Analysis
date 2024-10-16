#include "include/types/netcdf_reader_t.h"


netcdf_reader_t::netcdf_reader_t(
    const uint32_t &buffer_size) {
    dimid.resize(Dim_idx::n_dims);
    varid.resize(Par_idx::n_pars);
    this->n_timesteps_buffer = buffer_size;
    this->time_buffer_idx = 0;
    for (uint32_t i=0; i < buffer.size(); i++) {
        if (i == Par_idx::ascent || i == Par_idx::lat || i == Par_idx::lon) {
#if defined B_EIGHT
            buffer[i].resize(this->n_timesteps_buffer+21);
#else
            this->buffer[i].resize(this->n_timesteps_buffer+1);
#endif
        } else {
            this->buffer[i].resize(this->n_timesteps_buffer);
        }
    }
    std::fill(this->buffer[Par_idx::dom_i].begin(), this->buffer[Par_idx::dom_i].end(), 1);
    already_open = false;
    start_time_idx_given = false;
}


void netcdf_reader_t::load_vars() {
    varid.resize(Par_idx::n_pars);
    varid_once.resize(Par_once_idx::n_pars_once);

    int nvars = 0;
    asc600_bool = false;
    SUCCESS_OR_DIE(
        nc_inq(
            ncid,
            NULL,
            &nvars,
            NULL,
            NULL));
    for (int i=0; i < nvars; ++i) {
        char name[NC_MAX_NAME + 1];
        SUCCESS_OR_DIE(
            nc_inq_varname(
                ncid,
                i,
                name));
        if (std::strcmp(name, "asc600") == 0) {
            asc600_bool = true;
            varid[Par_idx::asc600] = i;
            break;
        }
    }
    dom_i_bool = false;
    SUCCESS_OR_DIE(
        nc_inq(
            ncid,
            NULL,
            &nvars,
            NULL,
            NULL));
    for (int i=0; i < nvars; ++i) {
        char name[NC_MAX_NAME + 1];
        SUCCESS_OR_DIE(
            nc_inq_varname(
                ncid,
                i,
                name));
        if (std::strcmp(name, "DOM_i") == 0) {
            dom_i_bool = true;
            varid[Par_idx::dom_i] = i;
            break;
        }
    }
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QV"),
            &varid_once[Par_once_idx::qv]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QC"),
            &varid_once[Par_once_idx::qc]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QR"),
            &varid_once[Par_once_idx::qr]));

#if defined(RK4ICE)
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QI"),
            &varid_once[Par_once_idx::qi]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QS"),
            &varid_once[Par_once_idx::qs]));
#endif
#if !defined(WCB) && defined(RK4ICE)
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QG"),
            &varid_once[Par_once_idx::qg]));
#endif

#if defined(B_EIGHT) && defined(RK4ICE)
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            "QH",
            &varid_once[Par_once_idx::qh]));
#endif

#if defined(RK4ICE) && (defined(WCB2) || defined(MET3D))

    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NCGRAUPEL"),
            &varid_once[Par_once_idx::ng]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NCICE"),
            &varid_once[Par_once_idx::ni]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NCSNOW"),
            &varid_once[Par_once_idx::ns]));
#endif

#if defined(B_EIGHT) && defined(RK4ICE)
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NCHAIL"),
            &varid_once[Par_once_idx::nh]));
#endif

#if defined WCB2 || defined MET3D || defined B_EIGHT
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NCCLOUD"),
            &varid_once[Par_once_idx::nc]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NCRAIN"),
            &varid_once[Par_once_idx::nr]));
#endif

    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("time"),
            &varid_once[Par_once_idx::time]));

    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("z"),
            &varid[Par_idx::height]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("lat"),
            &varid[Par_idx::lat]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("lon"),
            &varid[Par_idx::lon]));

    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("time_rel"),
            &varid[Par_idx::time_after_ascent]));

    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("pressure"),
            &varid[Par_idx::pressure]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("T"),
            &varid[Par_idx::temperature]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("w"),
            &varid[Par_idx::ascent]));

#if defined WCB2 || defined MET3D || defined B_EIGHT
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QI_IN"),
            &varid[Par_idx::qi_in]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QS_IN"),
            &varid[Par_idx::qs_in]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QR_IN"),
            &varid[Par_idx::qr_in]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QG_IN"),
            &varid[Par_idx::qg_in]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NI_IN"),
            &varid[Par_idx::ni_in]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NS_IN"),
            &varid[Par_idx::ns_in]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NR_IN"),
            &varid[Par_idx::nr_in]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NG_IN"),
            &varid[Par_idx::ng_in]));
    #if defined(B_EIGHT)
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("QH_IN"),
            &varid[Par_idx::qh_in]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("NH_IN"),
            &varid[Par_idx::nh_in]));
    #endif
    #if !defined(B_EIGHT)
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("conv_400"),
            &varid[Par_idx::conv_400]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("conv_600"),
            &varid[Par_idx::conv_600]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("slan_400"),
            &varid[Par_idx::slan_400]));
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("slan_600"),
            &varid[Par_idx::slan_600]));
    #endif
#endif
#if defined MET3D && defined TURBULENCE
    SUCCESS_OR_DIE(
        nc_inq_varid(
            ncid,
            reader_names.at("Q_TURBULENCE"),
            &varid[Par_idx::q_turb]));
#endif
}


void netcdf_reader_t::buffer_params(
    const reference_quantities_t &ref_quant) {
    // check amount of timesteps left to read
    // We need to reset n_timesteps_buffer in case it had been reduced
    // for another trajecotry or grid point that the process has
    // finished before
    n_timesteps_buffer = buffer[0].size();
    if (time_idx + n_timesteps_buffer > n_timesteps_in) {
        n_timesteps_buffer = n_timesteps_in - time_idx;
    }

    std::vector<size_t> startp, countp, countp2;
#if defined WCB || defined WCB2
    startp.push_back(time_idx);
    startp.push_back(traj_idx);
    countp.push_back(n_timesteps_buffer);
    countp.push_back(1);
    int additional_buffer = n_timesteps_in - (time_idx + n_timesteps_buffer);
    if (additional_buffer < 0) additional_buffer = 0;
    if (additional_buffer > 10) additional_buffer = 10;

    countp2.push_back(n_timesteps_buffer + additional_buffer);
    countp2.push_back(1);
#elif defined(MET3D) && !defined(B_EIGHT)
    startp.push_back(ens_idx);
    startp.push_back(traj_idx);
    startp.push_back(time_idx);

    countp.push_back(1);
    countp.push_back(1);
    countp.push_back(n_timesteps_buffer);
    int additional_buffer = n_timesteps_in - (time_idx + n_timesteps_buffer);
    if (additional_buffer < 0) additional_buffer = 0;
    if (additional_buffer > 10) additional_buffer = 10;

    countp2.push_back(1);
    countp2.push_back(1);
    countp2.push_back(n_timesteps_buffer + additional_buffer);
#else
    startp.push_back(traj_idx);
    startp.push_back(time_idx);

    countp.push_back(1);
    countp.push_back(n_timesteps_buffer);
    int additional_buffer = n_timesteps_in - (time_idx + n_timesteps_buffer);
    if (additional_buffer < 0) additional_buffer = 0;
    if (additional_buffer > 10) additional_buffer = 10;

    countp2.push_back(1);
    countp2.push_back(n_timesteps_buffer + additional_buffer);
#endif
    for (int i=0; i < Par_idx::n_pars; i++) {
        if (!asc600_bool && i == Par_idx::asc600) continue;
        if (!dom_i_bool && i == Par_idx::dom_i) continue;
        if (i == Par_idx::ascent || i == Par_idx::lat || i == Par_idx::lon) {
            SUCCESS_OR_DIE(
                nc_get_vara_double(
                    ncid,
                    varid[i],
                    startp.data(),
                    countp2.data(),
                    buffer[i].data()));
        } else {
            SUCCESS_OR_DIE(
                nc_get_vara_double(
                    ncid,
                    varid[i],
                    startp.data(),
                    countp.data(),
                    buffer[i].data()));
        }
    }
    for (auto &v : buffer[Par_idx::pressure]) {
        v *= pascal_conv/ref_quant.pref;
    }
    for (auto &v : buffer[Par_idx::temperature])
        v /= ref_quant.Tref;
#if defined(FLUX) && !defined(WCB)

    for (auto &v : buffer[Par_idx::qr_in])
        v /= ref_quant.qref;
    for (auto &v : buffer[Par_idx::nr_in])
        v /= ref_quant.Nref;
  #if defined(RK4ICE)
    for (auto &v : buffer[Par_idx::qi_in])
        v /= ref_quant.qref;
    for (auto &v : buffer[Par_idx::qs_in])
        v /= ref_quant.qref;
    for (auto &v : buffer[Par_idx::qg_in])
        v /= ref_quant.qref;
    for (auto &v : buffer[Par_idx::ni_in])
        v /= ref_quant.Nref;
    for (auto &v : buffer[Par_idx::ns_in])
        v /= ref_quant.Nref;
    for (auto &v : buffer[Par_idx::ng_in])
        v /= ref_quant.Nref;
    #if defined B_EIGHT
    for (auto &v : buffer[Par_idx::nh_in])
        v /= ref_quant.Nref;
    for (auto &v : buffer[Par_idx::qh_in])
        v /= ref_quant.qref;
    #endif
  #endif
#endif
    for (auto &v : buffer[Par_idx::ascent])
        v /= ref_quant.wref;
#if defined MET3D && defined TURBULENCE
    for (auto &v : buffer[Par_idx::q_turb])
        v /= ref_quant.qref;
#endif
    time_buffer_idx = 0;
}


template<class float_t>
int netcdf_reader_t::read_buffer(
    model_constants_t<float_t> &cc,
    const reference_quantities_t &ref_quant,
    std::vector<float_t> &y_single_old,
    std::vector<float_t> &inflows,
    const uint32_t &step,
    const bool &checkpoint_flag,
    const bool &start_over_env) {
    time_buffer_idx += step+start_time_idx-time_idx;
    time_idx = step+start_time_idx;

    // check if buffer needs to be reloaded
#if defined B_EIGHT
    if (time_buffer_idx + 20 >= n_timesteps_buffer) {
#else
    if (time_buffer_idx >= n_timesteps_buffer) {
#endif
        buffer_params(ref_quant);
        time_buffer_idx = 0;
    }
#ifdef TRACE_COMM
    std::cout << "read buffer: set time_buffer_idx: " << time_buffer_idx
              << ", with n_timesteps_buffer: " << n_timesteps_buffer << "\n";
#endif
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
    // Reset precipitation efficiency
    y_single_old[out_eff_idx] = 0;

    int err = SUCCESS;

    // Set values from a given trajectory
    if ((step == 0 && !checkpoint_flag) || start_over_env) {
        uint32_t i = 0;
        while (check_invalid_buffer(buffer[Par_idx::pressure], i)) {
            i++;
            if (i == 11) {
#ifdef TRUSTED_DATA
                err = INPUT_NAN_ERR;
#else
                std::cerr  << "Error at trajectory index " << traj_idx
                    << " and time index " << time_idx << ".\n"
                    << "At leat one trajectory has more than 10 "
                    << "consecutive time steps with NaNs in it. Check "
                    << "your dataset! Aborting now.\n";
                SUCCESS_OR_DIE(INPUT_NAN_ERR);
#endif
            }
        }
        if (i > 0 && std::isnan(buffer[Par_idx::pressure][time_buffer_idx+i])) {
            // The last couple of time steps are only NaN. Instead of aborting, we simply
            // reuse the old data, but now without any influx of course.
            // We also stop any ascent
            inflows[qr_in_idx] = 0;
            inflows[Nr_in_idx] = 0;
#if defined(RK4ICE)
            inflows[qi_in_idx] = 0;
            inflows[qs_in_idx] = 0;
            inflows[qg_in_idx] = 0;
            inflows[Ni_in_idx] = 0;
            inflows[Ns_in_idx] = 0;
            inflows[Ng_in_idx] = 0;
    #if defined(B_EIGHT)
            inflows[qh_in_idx] = 0;
            inflows[Nh_in_idx] = 0;
    #endif
#endif
#if defined MET3D && defined TURBULENCE
            inflows[qv_in_idx] = 0;
#endif
            y_single_old[w_idx] = 0;
            cc.set_dw_polynomial(0, 0, 0, 0, 0);
        } else if (i > 0) {
#ifdef DEVELOP
            std::cout << "read buffer: write at " << time_buffer_idx+i
                      << " with size " << buffer.size() << " and " << buffer[Par_idx::pressure].size() << "\n";
#endif
            y_single_old[p_idx] += (buffer[Par_idx::pressure][time_buffer_idx+i] - y_single_old[p_idx]) / (i+1);
            // If the dataset isn't broken, then I can just reuse the found i.
            // A check at the end of the function will be done in case the
            // given dataset is broken.
            y_single_old[T_idx] += (buffer[Par_idx::temperature][time_buffer_idx+i] - y_single_old[T_idx]) / (i+1);
#if defined(RK4ICE) || defined(RK4NOICE)
            y_single_old[w_idx] += (buffer[Par_idx::ascent][time_buffer_idx+i] - y_single_old[w_idx]) / (i+1);
            // No need to change dw here. It has been calculated for the case i == 0 or during initialization.
#endif
            y_single_old[z_idx] += (buffer[Par_idx::height][time_buffer_idx+i] - y_single_old[z_idx]) / (i+1);
            // set inflow
#if defined(FLUX) && !defined(WCB)
            inflows[qr_in_idx] =  (buffer[Par_idx::qr_in][time_buffer_idx+i] - inflows[qr_in_idx]) / (i+1);
            inflows[Nr_in_idx] =  (buffer[Par_idx::nr_in][time_buffer_idx+i] - inflows[Nr_in_idx]) / (i+1);
    #if defined(RK4ICE)
            inflows[qi_in_idx] =  (buffer[Par_idx::qi_in][time_buffer_idx+i] - inflows[qi_in_idx]) / (i+1);
            inflows[Ni_in_idx] =  (buffer[Par_idx::ni_in][time_buffer_idx+i] - inflows[Ni_in_idx]) / (i+1);
            inflows[qs_in_idx] =  (buffer[Par_idx::qs_in][time_buffer_idx+i] - inflows[qs_in_idx]) / (i+1);
            inflows[Ns_in_idx] =  (buffer[Par_idx::ns_in][time_buffer_idx+i] - inflows[Ns_in_idx]) / (i+1);
            inflows[qg_in_idx] =  (buffer[Par_idx::qg_in][time_buffer_idx+i] - inflows[qg_in_idx]) / (i+1);
            inflows[Ng_in_idx] =  (buffer[Par_idx::ng_in][time_buffer_idx+i] - inflows[Ng_in_idx]) / (i+1);
        #if defined(B_EIGHT)
            inflows[qh_in_idx] =  (buffer[Par_idx::qh_in][time_buffer_idx+i] - inflows[qh_in_idx]) / (i+1);
            inflows[Nh_in_idx] =  (buffer[Par_idx::nh_in][time_buffer_idx+i] - inflows[Nh_in_idx]) / (i+1);
        #endif
    #endif
#else
            inflows[qr_in_idx] = 0;
            inflows[Nr_in_idx] = 0;
    #if defined(RK4ICE)
            inflows[qi_in_idx] = 0;
            inflows[qs_in_idx] = 0;
            inflows[qg_in_idx] = 0;
            inflows[Ni_in_idx] = 0;
            inflows[Ns_in_idx] = 0;
            inflows[Ng_in_idx] = 0;
        #if defined(B_EIGHT)
            inflows[qh_in_idx] = 0;
            inflows[Nh_in_idx] = 0;
        #endif
    #endif
#endif
#if defined MET3D && defined TURBULENCE
            inflows[qv_in_idx] = (buffer[Par_idx::turbulence][time_buffer_idx+i] - inflows[qv_in_idx]) / (i+1);
#endif

        } else {
#ifdef DEVELOP
            std::cout << "read buffer: read at " << time_buffer_idx
                      << " with size " << buffer.size() << " and " << buffer[Par_idx::pressure].size() << "\n";
#endif
            y_single_old[p_idx] = buffer[Par_idx::pressure][time_buffer_idx];
            y_single_old[T_idx] = buffer[Par_idx::temperature][time_buffer_idx];
#if defined(RK4ICE) || defined(RK4NOICE)
            y_single_old[w_idx] = buffer[Par_idx::ascent][time_buffer_idx];
#ifdef DEVELOP
            std::cout << " test \n";
#endif
            if (time_idx == n_timesteps_in) {
                cc.set_dw_polynomial(buffer[Par_idx::ascent][time_buffer_idx], 0, 0, 0, 0);
            } else {
                // Here we need to find the next valid w value again
                // Unfortunately, some datasets may fill w with zeros instead of NaNs
                // therefore we check if pressure is NaN as well.
#ifdef DEVELOP
                std::cout << " test 2\n";
#endif
                uint32_t j = 1;
                while (check_invalid_buffer(buffer[Par_idx::ascent], j)) {
                    j++;
                    if (j == 11) break;
                }
                uint32_t j2 = j+1;
                while (check_invalid_buffer(buffer[Par_idx::ascent], j2)) {
                    j2++;
                    if (j2 == 11) break;
                }
#ifdef DEVELOP
                std::cout << " test 3\n";
#endif
                // Only NaNs
                if (j == 11 || std::isnan(buffer[Par_idx::ascent][time_buffer_idx+j])) {
                    cc.set_dw_polynomial(buffer[Par_idx::ascent][time_buffer_idx], 0, 0, 0, 0);
                } else if (j2 == 1 || std::isnan(buffer[Par_idx::ascent][time_buffer_idx+j2])) {
                    // Only one more valid datapoint
                    cc.set_dw_polynomial(
                            buffer[Par_idx::ascent][time_buffer_idx],
                            buffer[Par_idx::ascent][time_buffer_idx+j],
                            0,
                            j*cc.dt_traject_prime,
                            0);
                } else {
                    cc.set_dw_polynomial(
                        buffer[Par_idx::ascent][time_buffer_idx],
                        buffer[Par_idx::ascent][time_buffer_idx+j],
                        buffer[Par_idx::ascent][time_buffer_idx+j2],
                        j*cc.dt_traject_prime,
                        (j2-j)*cc.dt_traject_prime);
#ifdef DEVELOP
                    std::cout << " test 4\n";
#endif
                }
            }
#endif
            y_single_old[z_idx] = buffer[Par_idx::height][time_buffer_idx];
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
        #if defined(B_EIGHT)
            inflows[qh_in_idx] =  buffer[Par_idx::qh_in][time_buffer_idx];
            inflows[Nh_in_idx] =  buffer[Par_idx::nh_in][time_buffer_idx];
        #endif
    #endif
#else
            inflows[qr_in_idx] = 0;
            inflows[Nr_in_idx] = 0;
    #if defined(RK4ICE)
            inflows[qi_in_idx] = 0;
            inflows[qs_in_idx] = 0;
            inflows[qg_in_idx] = 0;
            inflows[Ni_in_idx] = 0;
            inflows[Ns_in_idx] = 0;
            inflows[Ng_in_idx] = 0;
        #if defined(B_EIGHT)
            inflows[qh_in_idx] = 0;
            inflows[Nh_in_idx] = 0;
        #endif
    #endif
#endif
#if defined MET3D && defined TURBULENCE
            inflows[qv_in_idx] = buffer[Par_idx::turbulence][time_buffer_idx];
#endif
        }
        for (auto &v : inflows) {
            if (std::isnan(v)) {
#ifdef TRUSTED_DATA
                err = INPUT_NAN_ERR;
#else
                std::cerr  << "Error at trajectory index " << traj_idx
                    << " and time index " << time_idx << ".\n"
                    << "At leat one trajectory has more than 10 "
                    << "consecutive time steps with NaNs in inflows. Check "
                    << "your dataset! Aborting now.\n";
                SUCCESS_OR_DIE(INPUT_NAN_ERR);
#endif
            } else {
                // some datasets set inflows as outflows from a box above
                // which can result in a negative value...
                v = std::abs(v);
            }
        }
        for (const auto &v : y_single_old) {
            if (std::isnan(v)) {
#ifdef TRUSTED_DATA
                err = INPUT_NAN_ERR;
#else
                std::cerr  << "Error at trajectory index " << traj_idx
                    << " and time index " << time_idx << ".\n"
                    << "At leat one trajectory has more than 10 "
                    << "consecutive time steps with NaNs in it. Check "
                    << "your dataset! Aborting now.\n";
                SUCCESS_OR_DIE(INPUT_NAN_ERR);
#endif
            }
        }
    }
#ifdef DEVELOP
    std::cout << " test end\n";
#endif
    if (inflows[qr_in_idx] > 0 && inflows[Nr_in_idx] == 0) {
        auto avg_size = get_at(cc.rain.constants, Particle_cons_idx::min_x) + (
            get_at(cc.rain.constants, Particle_cons_idx::max_x)
            - get_at(cc.rain.constants, Particle_cons_idx::min_x)) / 2;
        inflows[Nr_in_idx] = inflows[qr_in_idx] / avg_size;
    }
    if (inflows[qi_in_idx] > 0 && inflows[Ni_in_idx] == 0) {
        auto avg_size = get_at(cc.ice.constants, Particle_cons_idx::min_x) + (
            get_at(cc.ice.constants, Particle_cons_idx::max_x)
            - get_at(cc.ice.constants, Particle_cons_idx::min_x)) / 2;
        inflows[Ni_in_idx] = inflows[qi_in_idx] / avg_size;
    }
    if (inflows[qs_in_idx] > 0 && inflows[Ns_in_idx] == 0) {
        auto avg_size = get_at(cc.snow.constants, Particle_cons_idx::min_x) + (
            get_at(cc.snow.constants, Particle_cons_idx::max_x)
            - get_at(cc.snow.constants, Particle_cons_idx::min_x)) / 2;
        inflows[Ns_in_idx] = inflows[qs_in_idx] / avg_size;
    }
    if (inflows[qg_in_idx] > 0 && inflows[Ng_in_idx] == 0) {
        auto avg_size = get_at(cc.graupel.constants, Particle_cons_idx::min_x) + (
            get_at(cc.graupel.constants, Particle_cons_idx::max_x)
            - get_at(cc.graupel.constants, Particle_cons_idx::min_x)) / 2;
        inflows[Ng_in_idx] = inflows[qg_in_idx] / avg_size;
    }
#if defined(B_EIGHT)
    if (inflows[qh_in_idx] > 0 && inflows[Nh_in_idx] == 0) {
        auto avg_size = get_at(cc.hail.constants, Particle_cons_idx::min_x)
            + (get_at(cc.hail.constants, Particle_cons_idx::max_x)
            - get_at(cc.hail.constants, Particle_cons_idx::min_x)) / 2;
        inflows[Nh_in_idx] = inflows[qh_in_idx] / avg_size;
    }
#endif
    return err;
}


template<class float_t>
void netcdf_reader_t::set_dims(
    const char *input_file,
    model_constants_t<float_t> &cc,
    const int &simulation_mode) {

    set_dims(input_file, simulation_mode);
    if (simulation_mode != limited_time_ensembles) {
        cc.n_ensembles = n_ensembles;
        cc.max_n_trajs = n_trajectories;
    }
}


void netcdf_reader_t::set_dims(
    const char *input_file,
    const int &simulation_mode) {
    time_buffer_idx = 0;
    traj_idx = 0;
    ens_idx = 0;

    if (!already_open) {
        if ((simulation_mode == grid_sensitivity) || (simulation_mode == trajectory_sensitivity)) {
            SUCCESS_OR_DIE(
                nc_open_par(
                    input_file,
                    NC_NOWRITE,
                    MPI_COMM_WORLD,
                    MPI_INFO_NULL,
                    &ncid));
        } else {
            SUCCESS_OR_DIE(
                nc_open(
                    input_file,
                    NC_NOWRITE,
                    &ncid));
        }

        SUCCESS_OR_DIE(
            nc_inq_dimid(
                ncid,
                reader_names.at("trajectory"),
                &dimid[Dim_idx::trajectory_dim_idx]));
#if !defined(B_EIGHT)
        SUCCESS_OR_DIE(
            nc_inq_dimid(
                ncid,
                reader_names.at("ensemble"),
                &dimid[Dim_idx::ensemble_dim_idx]));
#endif
    }
    SUCCESS_OR_DIE(
        nc_inq_dimlen(
            ncid,
            dimid[Dim_idx::trajectory_dim_idx],
            &n_trajectories));
#if !defined(B_EIGHT)
    SUCCESS_OR_DIE(
        nc_inq_dimlen(
            ncid,
            dimid[Dim_idx::ensemble_dim_idx],
            &n_ensembles));
#else
    n_ensembles = 1;
#endif
    if (n_trajectories <= traj_idx) {
        std::cerr << "Number of trajectories in netCDF file: " << n_trajectories << "\n";
        std::cerr << "You asked for trajectory with index " << traj_idx
                    << " which does not exist. ABORTING.\n";
        SUCCESS_OR_DIE(NC_TRAJ_IDX_ERR);
    }

    if (!already_open) {
        SUCCESS_OR_DIE(
            nc_inq_dimid(
                ncid,
                reader_names.at("time"),
                &dimid[Dim_idx::time_dim_idx]));
        SUCCESS_OR_DIE(
            nc_inq_dimlen(
                ncid,
                dimid[Dim_idx::time_dim_idx],
                &n_timesteps_in));

        load_vars();
        already_open = true;
        // Check if the data is stored as hPa or Pa and change pref in
        // reference quantities accordingly.
        size_t att_len;
        if (0 == nc_inq_att(ncid, varid[Par_idx::pressure], "units", NULL, &att_len)) {
            std::vector<char> att_val(att_len+1);
            SUCCESS_OR_DIE(nc_get_att(
                ncid, varid[Par_idx::pressure], "units", att_val.data()));
            att_val[att_len] = '\0';
            if (std::strcmp(att_val.data(), "hPa") == 0) {
                pascal_conv = 100;
            } else if (std::strcmp(att_val.data(), "Pa") == 0) {
                pascal_conv = 1;
            } else {
                pascal_conv = 1;  // default
            }
        } else {
            pascal_conv = 1;  // default
        }
    }
}


template<class float_t>
void netcdf_reader_t::init_netcdf(
#ifdef MET3D
    double &start_time,
#endif
    const bool &checkpoint_flag,
    model_constants_t<float_t> &cc,
    const double current_time,
    const reference_quantities_t &ref_quant) {
    std::vector<size_t> startp, countp;
#if defined(MET3D) && !defined(B_EIGHT)
    countp.push_back(1);
    countp.push_back(1);
    countp.push_back(1);

    startp.push_back(ens_idx);      // ensemble
    startp.push_back(traj_idx);     // trajectory
    startp.push_back(0);            // time
    uint64_t start_time_idx = 0;
    double rel_start_time;

    SUCCESS_OR_DIE(
        nc_get_var1_double(
            ncid,
            varid[Par_idx::time_after_ascent],
            startp.data(),
            &rel_start_time));

    std::vector<double> time(2);
    std::vector<size_t> startp2, countp2;
    countp2.push_back(2);
    startp2.push_back(0);

    SUCCESS_OR_DIE(
        nc_get_vara_double(
            ncid,
            varid_once[Par_once_idx::time],
            startp2.data(),
            countp2.data(),
            time.data()));
    cc.set_dt(time[1]-time[0], ref_quant);

    if (this->start_time_idx_given) {
#ifdef TRACE_COMM
      std::cout << "traj: " << traj_idx << " start_time_idx given with " << this->start_time_idx << "\n";
#endif
        start_time_idx = this->start_time_idx;
    } else if (!std::isnan(start_time) && !checkpoint_flag) {
        // Calculate the needed index
        start_time_idx = (start_time-rel_start_time)/cc.dt_traject;
#ifdef TRACE_COMM
        std::cout << "traj: " << traj_idx << " start_time not nan and not a checkpoint. start_time: " << start_time
                  << ", rel_start_time: " << rel_start_time << ", dt_traject: " << cc.dt_traject
                  << ", start_time_idx: " << start_time_idx << "\n";
#endif
    } else if (checkpoint_flag && !std::isnan(start_time)) {
        // Calculate the needed index
        start_time_idx = ceil(start_time-rel_start_time + current_time)/cc.dt_traject;
#ifdef TRACE_COMM
        std::cout << "traj: " << traj_idx << " start_time not nan and is a checkpoint. start_time: " << start_time
                  << ", rel_start_time: " << rel_start_time << ", dt_traject: " << cc.dt_traject
                  << ", current_time: " << current_time << ", start_time_idx: " << start_time_idx << "\n";
#endif
    }
    time_idx = start_time_idx;
    this->start_time_idx = start_time_idx;
#elif defined(B_EIGHT)
    countp.push_back(1);
    countp.push_back(10);

    startp.push_back(traj_idx);     // trajectory
    startp.push_back(0);            // time

    uint64_t start_time_idx = 0;
    if (this->start_time_idx_given) {
        start_time_idx = this->start_time_idx;
    }
    std::vector<double> rel_start_time(10);
    SUCCESS_OR_DIE(
        nc_get_vara_double(
            ncid,
            varid[Par_idx::time_after_ascent],
            startp.data(),
            countp.data(),
            rel_start_time.data()));

    std::vector<double> time(2);
    std::vector<size_t> startp2, countp2;
    countp2.push_back(2);
    startp2.push_back(0);

    SUCCESS_OR_DIE(
        nc_get_vara_double(
            ncid,
            varid_once[Par_once_idx::time],
            startp2.data(),
            countp2.data(),
            time.data()));

    cc.set_dt(time[1]-time[0], ref_quant);
    if (!this->start_time_idx_given) {
        bool all_nan = true;
        for (uint32_t i = 0; i < rel_start_time.size() - 1; i++) {
            if (std::isnan(rel_start_time[i])) continue;
            if (!std::isnan(start_time) && !checkpoint_flag) {
                // Calculate the needed index
                if (rel_start_time[i + 1] - rel_start_time[i] < cc.dt_traject) {
                    start_time_idx = i + (start_time - rel_start_time[i] * cc.dt_traject) / cc.dt_traject;
                } else {
                    start_time_idx = i + (start_time - rel_start_time[i]) / cc.dt_traject;
                }
                all_nan = false;
                break;
            } else if (checkpoint_flag && !std::isnan(start_time)) {
                // Calculate the needed index
                if (rel_start_time[i + 1] - rel_start_time[i] < cc.dt_traject) {
                    start_time_idx = i + ceil(start_time - rel_start_time[i]
                                                           * cc.dt_traject + current_time) / cc.dt_traject;
                } else {
                    start_time_idx = i + ceil(start_time - rel_start_time[i] + current_time) / cc.dt_traject;
                }
                all_nan = false;
                break;
            }
            all_nan = false;
        }
        if (all_nan) {
            std::cerr << "Error at trajectory index " << traj_idx << "\n"
                      << " and time index " << time_idx << ".\n"
                      << "The first 10 values of the column 'time_after_asc_start' are "
                      << "NaN. Check your dataset. Aborting.\n";
            SUCCESS_OR_DIE(NC_TRAJ_IDX_ERR);
        }
    }
    // Check if the index begins with non-NaN values to get started properly
    startp[1] = start_time_idx;
    SUCCESS_OR_DIE(
        nc_get_vara_double(
            ncid,
            varid[Par_idx::pressure],  // In case the time is always non-NaN
            startp.data(),
            countp.data(),
            rel_start_time.data()));
    for (uint32_t i=0; i < rel_start_time.size(); i++) {
        if (std::isnan(rel_start_time[i])) {
            start_time_idx++;
            cc.num_steps--;
            continue;
        } else {
            break;
        }
    }
    time_idx = start_time_idx;
    this->start_time_idx = start_time_idx;
#else
    cc.set_dt(20);
    time_idx = 1;
    this->start_time_idx = 1;
    cc.start_time = 0;
#endif
    this->start_time_idx_original = this->start_time_idx;
#ifdef TRACE_COMM
    std::cout << "Init_netcdf, start_time_idx: " << this->start_time_idx << "\n";
#endif

#if defined(MET3D) || defined(B_EIGHT)
    std::vector<size_t> startp3;
    startp3.push_back(time_idx);
    SUCCESS_OR_DIE(
        nc_get_var1_double(
            ncid,
            varid_once[Par_once_idx::time],
            startp3.data(),
            &(cc.start_time)));
#endif
}


template<class float_t>
void netcdf_reader_t::read_initial_values(
    std::vector<double> &y_init,
    const reference_quantities_t &ref_quant,
    model_constants_t<float_t> &cc,
    const bool &checkpoint_flag,
    const uint64_t &traj_id,
    const uint64_t &ens_id) {
    traj_idx = traj_id;
    ens_idx = ens_id;
    this->read_initial_values(y_init, ref_quant, cc, checkpoint_flag);
}


template<class float_t>
void netcdf_reader_t::read_initial_values(
    std::vector<double> &y_init,
    const reference_quantities_t &ref_quant,
    model_constants_t<float_t> &cc,
    const bool &checkpoint_flag) {

    buffer_params(ref_quant);
    std::vector<size_t> startp;
#ifdef TRACE_COMM
    std::cout << "read_initial_values checkpoint: " << checkpoint_flag
        << ", traj_idx: " << traj_idx << ", cc.traj_id: " << cc.traj_id
        << ", start_time_idx: " << start_time_idx << "\n";
#endif
    if (!checkpoint_flag) {
#if defined B_EIGHT
        startp.push_back(traj_idx);
        startp.push_back(start_time_idx);
#elif defined MET3D
        startp.push_back(ens_idx);
        startp.push_back(traj_idx);
        startp.push_back(start_time_idx);
#elif defined WCB || WCB2
        startp.push_back(start_time_idx);
        startp.push_back(traj_idx);
#else
        startp.push_back(traj_idx);
        startp.push_back(start_time_idx);
#endif
        y_init[T_idx] = buffer[Par_idx::temperature][0];
        y_init[p_idx] = buffer[Par_idx::pressure][0];

        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qv],
                startp.data(),
                &(y_init[qv_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qc],
                startp.data(),
                &(y_init[qc_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qr],
                startp.data(),
                &(y_init[qr_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qv],
                startp.data(),
                &(y_init[qv_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qc],
                startp.data(),
                &(y_init[qc_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qr],
                startp.data(),
                &(y_init[qr_idx])));
#if defined(RK4ICE)
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qi],
                startp.data(),
                &(y_init[qi_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qs],
                startp.data(),
                &(y_init[qs_idx])));
#endif
#if !defined(WCB) && defined(RK4ICE)
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qg],
                startp.data(),
                &(y_init[qg_idx])));
#elif defined(RK4ICE)
        y_init[qg_idx] = 0;
#endif
#if defined(RK4ICE)
    #if defined(B_EIGHT)
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::qh],
                startp.data(),
                &(y_init[qh_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::nh],
                startp.data(),
                &(y_init[Nh_idx])));
    #else
        y_init[qh_idx] = 0.0;  // qh. We don't have hail in the trajectoris
        y_init[Nh_idx] = 0.0;  // Nh. We don't have hail in the trajectoris
    #endif
#endif
#ifdef DEVELOP
        std::cout << "Got a lot\n" << std::flush;
#endif
#if defined(RK4ICE) && (defined(WCB2) || defined(MET3D))
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::ng],
                startp.data(),
                &(y_init[Ng_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::ni],
                startp.data(),
                &(y_init[Ni_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::ns],
                startp.data(),
                &(y_init[Ns_idx])));

        y_init[qi_out_idx] = 0;
        y_init[qs_out_idx] = 0;
        y_init[qg_out_idx] = 0;
        y_init[qh_out_idx] = 0;
        y_init[Nh_out_idx] = 0;
        y_init[Ni_out_idx] = 0;
        y_init[Ns_out_idx] = 0;
        y_init[Ng_out_idx] = 0;
#endif
        y_init[qr_out_idx] = 0;
        y_init[Nr_out_idx] = 0;
#if defined WCB2 || defined MET3D
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::nc],
                startp.data(),
                &(y_init[Nc_idx])));
        SUCCESS_OR_DIE(
            nc_get_var1_double(
                ncid,
                varid_once[Par_once_idx::nr],
                startp.data(),
                &(y_init[Nr_idx])));
#else
        float_t denom = 0;
        denom = (get_at(cc.cloud.max_x - cc.cloud.constants, Particle_cons_idx::min_x))
            / 2.0 + get_at(cc.cloud.constants, Particle_cons_idx::min_x);
        y_init[Nc_idx] = y_init[qc_idx] / (denom);
        denom = (get_at(cc.rain.max_x - cc.rain.constants, Particle_cons_idx::min_x))
            / 2 + get_at(cc.rain.constants, Particle_cons_idx::min_x);
        y_init[Nr_idx] = y_init[qr_idx] / (denom);
        denom = get_at(cc.cloud.constants, Particle_cons_idx::min_x) / 2.0;
    #if defined(RK4ICE)
        denom = (get_at(cc.ice.max_x - cc.ice.constants, Particle_cons_idx::min_x))
            / 2.0 + get_at(cc.ice.constants, Particle_cons_idx::min_x);
        y_init[Ni_idx] = y_init[qi_idx] / (denom);
        denom = (get_at(cc.snow.max_x - cc.snow.constants, Particle_cons_idx::min_x))
            / 2.0 + get_at(cc.snow.constants, Particle_cons_idx::min_x);
        y_init[Ns_idx] = y_init[qs_idx] / (denom);
        denom = (get_at(cc.graupel.max_x - cc.graupel.constants, Particle_cons_idx::min_x))
            / 2.0 + get_at(cc.graupel.constants, Particle_cons_idx::min_x);
        y_init[Ng_idx] = y_init[qg_idx] / (denom);
    #endif
#endif
        // Actually only needed for single moment schemes
        cc.constants[static_cast<int>(Cons_idx::Nc_prime)] = y_init[Nc_idx];

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
    #if defined(B_EIGHT)
        uint32_t j = 1;
        while (j < 11 && std::isnan(buffer[Par_idx::ascent][j])) {
            j++;
        }
        uint32_t j2 = j+1;
        while (check_invalid_buffer(buffer[Par_idx::ascent], j2)) {
            j2++;
            if (j2 == 11) break;
        }
        // Only NaNs
        if (j == 11) {
#ifdef TRUSTED_DATA
            cc.set_dw_polynomial(buffer[Par_idx::ascent][time_buffer_idx], 0, 0, 0, 0);
#else
            std::cerr  << "Error at trajectory index " << traj_idx
                << " and time index " << time_idx << ".\n"
                << "At leat one trajectory has more than 10 "
                << "consecutive time steps with NaNs in 'w'. Check "
                << "your dataset! Aborting now.\n";
            SUCCESS_OR_DIE(NC_TRAJ_IDX_ERR);
#endif
        } else if (j2 == 11) {
            // Only one more valid datapoint
            cc.set_dw_polynomial(
                    buffer[Par_idx::ascent][0],
                    buffer[Par_idx::ascent][j],
                    0,
                    j * cc.dt_traject_prime,
                    0);
        } else {
            cc.set_dw_polynomial(
                buffer[Par_idx::ascent][0],
                buffer[Par_idx::ascent][j],
                buffer[Par_idx::ascent][j2],
                j * cc.dt_traject_prime,
                (j2-j) * cc.dt_traject_prime);
        }
    #else
        cc.set_dw_polynomial(
                    buffer[Par_idx::ascent][0],
                    buffer[Par_idx::ascent][1],
                    buffer[Par_idx::ascent][2],
                    cc.dt_traject_prime,
                    cc.dt_traject_prime);
    #endif
        y_init[n_inact_idx] = 0;
        y_init[depo_idx] = 0;
        y_init[sub_idx] = 0;
#endif
#ifdef DEVELOP
        std::cout << "Got time coordinate\n" << std::flush;
#endif
#ifdef WCB
        y_init[w_idx] = 0;
#else
        y_init[w_idx] = buffer[Par_idx::ascent][0];
#endif
        n_subs = cc.num_sub_steps;
    }
}


double netcdf_reader_t::get_lat(
    const uint32_t &t,
    const uint32_t &sub) const {

    double lower, upper;
    uint32_t j = 0;
    while (check_invalid_buffer_lonlat(buffer[Par_idx::lat], t, j)) {
        j++;
    }
    lower = buffer[Par_idx::lat][(t%(buffer[Par_idx::lat].size()-1))+j];
    uint32_t difference = j;
    if (j == 11) {
#ifdef TRUSTED_DATA
        return buffer[Par_idx::lat][t%(buffer[Par_idx::lat].size()-1)];
#else
        std::cerr  << "Error at trajectory index " << traj_idx
                    << " and time index " << time_idx + t << ".\n"
                    << "At leat one trajectory has more than 10 "
                    << "consecutive time steps with NaNs in 'latitude'. Check "
                    << "your dataset! Aborting now.\n";
        SUCCESS_OR_DIE(INPUT_NAN_ERR);
#endif
    }
    j++;
    while (check_invalid_buffer_lonlat(buffer[Par_idx::lat], t, j)) {
        j++;
    }
    if (j >= 11) {
#ifdef TRUSTED_DATA
        return buffer[Par_idx::lat][t%(buffer[Par_idx::lat].size()-1)];
#else
        std::cerr  << "Error at trajectory index " << traj_idx
                    << " and time index " << time_idx + t << ".\n"
                    << "At leat one trajectory has more than 10 "
                    << "consecutive time steps with NaNs in 'latitude'. Check "
                    << "your dataset! Aborting now.\n";
        SUCCESS_OR_DIE(INPUT_NAN_ERR);
#endif
    }
    upper = buffer[Par_idx::lat][(t%(buffer[Par_idx::lat].size()-1))+j];
    difference = j - difference;
    return lower + sub * (upper - lower) / (n_subs*difference);
}


double netcdf_reader_t::get_lon(
    const uint32_t &t,
    const uint32_t &sub) const {
    double lower, upper;
    uint32_t j = 0;
    while (check_invalid_buffer_lonlat(buffer[Par_idx::lon], t, j)) {
        j++;
    }
    lower = buffer[Par_idx::lon][(t%(buffer[Par_idx::lon].size()-1))+j];
    uint32_t difference = j;
    if (j == 11) {
#ifdef TRUSTED_DATA
        return buffer[Par_idx::lon][t%(buffer[Par_idx::lon].size()-1)];
#else
        std::cerr  << "Error at trajectory index " << traj_idx
                    << " and time index " << time_idx + t << ".\n"
                    << "At leat one trajectory has more than 10 "
                    << "consecutive time steps with NaNs in 'latitude'. Check "
                    << "your dataset! Aborting now.\n";
        SUCCESS_OR_DIE(INPUT_NAN_ERR);
#endif
    }
    j++;
    while (check_invalid_buffer_lonlat(buffer[Par_idx::lon], t, j)) {
        j++;
    }
    if (j >= 11) {
#ifdef TRUSTED_DATA
        return buffer[Par_idx::lon][t%(buffer[Par_idx::lon].size()-1)];
#else
        std::cerr  << "Error at trajectory index " << traj_idx
                    << " and time index " << time_idx + t << ".\n"
                    << "At leat one trajectory has more than 10 "
                    << "consecutive time steps with NaNs in 'latitude'. Check "
                    << "your dataset! Aborting now.\n";
        SUCCESS_OR_DIE(INPUT_NAN_ERR);
#endif
    }
    upper = buffer[Par_idx::lon][(t%(buffer[Par_idx::lon].size()-1))+j];
    difference = j - difference;
    return lower + sub * (upper - lower) / (n_subs*difference);
}


bool netcdf_reader_t::check_invalid_buffer(
    const std::vector<double> &buffer_var,
    const uint32_t idx) const {

    // if different domains are available, check the domain to skip invalid values as well.
    return time_idx+idx <= n_timesteps_in
        && (std::isnan(buffer_var[time_buffer_idx+idx]) || buffer[Par_idx::dom_i][time_buffer_idx+idx] == 0);
}


bool netcdf_reader_t::check_invalid_buffer_lonlat(
    const std::vector<double> &buffer_var,
    const uint32_t t,
    const uint32_t idx) const {

    return idx < 11
        && (std::isnan(buffer_var[(t%(buffer_var.size()-1))+idx])
            || buffer[Par_idx::dom_i][(t%(buffer_var.size()-1))+idx] == 0);
}

template int netcdf_reader_t::read_buffer<codi::RealReverseIndex>(
    model_constants_t<codi::RealReverseIndex>&,
    const reference_quantities_t&,
    std::vector<codi::RealReverseIndex>&,
    std::vector<codi::RealReverseIndex>&,
    const uint32_t&,
    const bool&,
    const bool&);

template int netcdf_reader_t::read_buffer<codi::RealForwardVec<num_par_init> >(
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const reference_quantities_t&,
    std::vector<codi::RealForwardVec<num_par_init> >&,
    std::vector<codi::RealForwardVec<num_par_init> >&,
    const uint32_t&,
    const bool&,
    const bool&);

template void netcdf_reader_t::set_dims<codi::RealReverseIndex>(
    const char*,
    model_constants_t<codi::RealReverseIndex>&,
    const int&);

template void netcdf_reader_t::set_dims<codi::RealForwardVec<num_par_init> >(
    const char*,
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const int&);

template void netcdf_reader_t::init_netcdf<codi::RealReverseIndex>(
#ifdef MET3D
    double&,
#endif
    const bool&,
    model_constants_t<codi::RealReverseIndex>&,
    const double,
    const reference_quantities_t &);

template void netcdf_reader_t::init_netcdf<codi::RealForwardVec<num_par_init> >(
#ifdef MET3D
    double&,
#endif
    const bool&,
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const double,
    const reference_quantities_t &);

template void netcdf_reader_t::read_initial_values<codi::RealReverseIndex>(
    std::vector<double>&,
    const reference_quantities_t&,
    model_constants_t<codi::RealReverseIndex>&,
    const bool&,
    const uint64_t&,
    const uint64_t&);

template void netcdf_reader_t::read_initial_values<codi::RealForwardVec<num_par_init> >(
    std::vector<double>&,
    const reference_quantities_t&,
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const bool&,
    const uint64_t&,
    const uint64_t&);

template void netcdf_reader_t::read_initial_values<codi::RealReverseIndex>(
    std::vector<double>&,
    const reference_quantities_t&,
    model_constants_t<codi::RealReverseIndex>&,
    const bool&);

template void netcdf_reader_t::read_initial_values<codi::RealForwardVec<num_par_init> >(
    std::vector<double>&,
    const reference_quantities_t&,
    model_constants_t<codi::RealForwardVec<num_par_init> >&,
    const bool&);
