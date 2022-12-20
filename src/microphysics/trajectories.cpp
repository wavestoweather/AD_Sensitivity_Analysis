#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>

#include "codi.hpp"

#include "include/microphysics/constants.h"
#include "include/types/checkpoint_t.h"
#include "include/types/global_args_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/output_handle_t.h"
#include "include/types/model_constants_t.h"
#include "include/types/netcdf_reader_t.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/segment_t.h"
#include "include/types/task_scheduler_t.h"

#include "include/microphysics/physical_parameterizations.h"
#include "include/microphysics/program_io.h"

#include "include/misc/general.h"


#if defined(RK4) || defined(RK4NOICE) || defined(RK4ICE) || defined(RK4_ONE_MOMENT)
#include "include/microphysics/rk4.h"
#endif

#ifdef EXPLICIT_EULER
#include "include/microphysics/euler.h"
#endif

#ifdef IMPLICIT_EULER
#include "include/microphysics/implicit_euler.h"
#endif

#include "include/misc/pbar.h"


void parse_args(
    const int &argc,
    char* const * argv,
    const int &rank,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant) {
    SUCCESS_OR_DIE(global_args.parse_arguments(argc, argv));

    ref_quant.Tref = 273.15;
#ifdef WCB
    ref_quant.qref = 1.0e-6;
#else
    ref_quant.qref = 1.0e-4;
#endif
    ref_quant.pref = 1.0e5;
    ref_quant.wref = 1.;  // 10.0
    ref_quant.tref = 1.0;
    ref_quant.zref = 1.0;
    ref_quant.Nref = 1.0;   // DUMMY
    input.set_input_from_arguments(global_args);
    if (rank == 0) print_reference_quantities(ref_quant);
}

template<class float_t>
model_constants_t<float_t> prepare_constants(
    const int &rank,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    std::vector<segment_t> &segments,
    std::vector<double> &y_init,
    checkpoint_t &checkpoint) {

    model_constants_t<float_t> cc(input.tracking_filename, input.track_initial_cond);

    if (global_args.checkpoint_flag && rank == 0) {
        checkpoint.load_checkpoint(global_args.checkpoint_string, cc,
            y_init, segments, input);
        print_segments(segments);
    } else {
        cc.setup_model_constants(input, ref_quant);
        if (global_args.ens_config_flag && (rank == 0 || input.simulation_mode == create_train_set)) {
            load_ens_config(global_args.ens_config_string, cc,
                segments, input, ref_quant);
            for (auto &s : segments)
                SUCCESS_OR_DIE(s.check());
            if (rank == 0)
                print_segments(segments);
            if (input.simulation_mode == limited_time_ensembles) {
                cc.n_ensembles = 1;
                cc.max_n_trajs = segments[0].n_members;
            }
        }
    }

    // Communicate amount of ensembles and trajectories
    SUCCESS_OR_DIE(MPI_Bcast(
        &cc.n_ensembles, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD));
    SUCCESS_OR_DIE(MPI_Bcast(
        &cc.max_n_trajs, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD));

    auto_type = input.auto_type;
    load_lookup_table(cc.ltabdminwgg);
    return cc;
}

template<class float_t>
void setup_simulation_base(
    const int &rank,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    model_constants_t<float_t> &cc,
    bool &already_loaded,
    netcdf_reader_t &netcdf_reader) {
    if (!already_loaded) {
        netcdf_reader.init_netcdf(
#ifdef MET3D
            input.start_time,
#endif
            (global_args.checkpoint_flag || already_loaded),
            cc,
            input.current_time, ref_quant);
        if (rank != 0)
            already_loaded = true;
    } else {
        netcdf_reader.time_idx = netcdf_reader.start_time_idx;
    }
#if defined(RK4_ONE_MOMENT)
    cc.setCoefficients(y_init[0] , y_init[1]);
#endif
    if (rank == 0 && !already_loaded) {
        cc.print();
        input.print_parameters();
    }

    if (rank == 0 && !already_loaded) {
#ifdef TRACE_QC
        cc.cloud.print("cloud");
#endif

#ifdef TRACE_QR
        cc.rain.print("rain");
#endif

#ifdef TRACE_QI
        cc.ice.print("ice");
#endif

#ifdef TRACE_QS
        cc.snow.print("snow");
#endif

#ifdef TRACE_QG
        cc.graupel.print("graupel");
#endif

#ifdef TRACE_QH
        cc.hail.print("hail");
#endif
        already_loaded = true;
    }
}


template<class float_t>
void setup_simulation(
    const int &rank,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    std::vector<segment_t> &segments,
    model_constants_t<float_t> &cc,
    std::vector<double> &y_init,
    checkpoint_t &checkpoint,
    output_handle_t &out_handler,
    bool &already_loaded,
    netcdf_reader_t &netcdf_reader) {
    checkpoint.load_checkpoint(cc, y_init, segments, input, out_handler);
    if (cc.checkpoint_steps > 0)
        netcdf_reader.start_time_idx = netcdf_reader.start_time_idx_original + cc.checkpoint_steps - 1;
    setup_simulation_base(rank, input,
            global_args, ref_quant, cc, already_loaded, netcdf_reader);
}

#if defined(TRACE_SAT) || defined(TRACE_QR) || defined(TRACE_QV) || defined(TRACE_QC) || defined(TRACE_QI) \
    || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH) || defined(TRACE_TIME)
template<class float_t>
void substep_trace(
    const uint32_t &sub,
    const uint32_t &t,
    const model_constants_t<float_t> &cc,
    const input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const std::vector<float_t> &y_single_old,
    const std::vector<float_t> &inflow) {

#if defined(TRACE_TIME)
#if defined(MET3D)
    trace = (((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) + input.start_time >= trace_start)
        && ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) + input.start_time <= trace_end));
#else
    trace = (((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime)  >= trace_start)
        && ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) <= trace_end));
#endif
#endif
    if (trace)
        std::cout << "\n################################\n"
                  << cc.id << " timestep : " << (sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime)
                  << "\n";
#if defined(TRACE_SAT)
    if (trace)
        std::cout << "Start qv_prime " << y_single_old[qv_idx]*ref_quant.qref << "\n";
#endif
#ifdef TRACE_QR
    if (trace) {
        std::cout << "Adding qr " << inflow[qr_in_idx]/cc.num_sub_steps*ref_quant.qref
            << ", Nr " << inflow[Nr_in_idx]/cc.num_sub_steps << "\n";
        std::cout << "qr_old " << y_single_old[qr_idx]*ref_quant.qref
            << ", Nr_old " << y_single_old[Nr_idx] << "\n";
    }
#endif
#ifdef TRACE_QI
    if (trace) {
        std::cout << "Adding qi " << inflow[qi_in_idx]/cc.num_sub_steps*ref_quant.qref
            << ", Ni " << inflow[Ni_in_idx]/cc.num_sub_steps << "\n";
        std::cout << "qi_old " << y_single_old[qi_idx]*ref_quant.qref
            << ", Ni_old " << y_single_old[Ni_idx] << "\n";
    }
#endif
#ifdef TRACE_QS
    if (trace) {
        std::cout << "Adding qs " << inflow[qs_in_idx]/cc.num_sub_steps*ref_quant.qref
            << ", Ns " << inflow[Ns_in_idx]/cc.num_sub_steps << "\n";
        std::cout << "qs_old " << y_single_old[qs_idx]*ref_quant.qref
            << ", Ns_old " << y_single_old[Ns_idx] << "\n";
    }
#endif
#ifdef TRACE_QG
    if (trace) {
        std::cout << "Adding qg " << inflow[qg_in_idx]/cc.num_sub_steps*ref_quant.qref
            << ", Ng " << inflow[Ng_in_idx]/cc.num_sub_steps << "\n";
        std::cout << "qg_old " << y_single_old[qg_idx]*ref_quant.qref
            << ", Ng_old " << y_single_old[Ng_idx] << "\n";
    }
#endif

#if defined TRACE_QV && defined MET3D && defined TURBULENCE
    if (trace) {
        std::cout << "Adding qv " << inflow[qv_in_idx]/cc.num_sub_steps*ref_quant.qref
            << "\n";
        std::cout << "qv_old " << y_single_old[qv_idx]*ref_quant.qref
            << "\n";
    }
#endif
}
#endif

template<class float_t>
void parameter_check(
    std::vector<segment_t> &segments,
    model_constants_t<float_t> &cc,
    const double &time_old,
    const std::vector< std::array<double, num_par > > &y_diff,
    const std::vector<float_t> &y_single_old,
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    task_scheduler_t &scheduler) {
    for (auto &s : segments) {
        bool perturb = s.perturb_check(
            y_diff, y_single_old, time_old);
        if (perturb) {
            // Perturb this instance
            if (s.n_members == 1) {
                std::string descr;
                s.perturb(cc, ref_quant, input, descr);
            } else {
                // Create a checkpoint for new members
                uint64_t ens_id;
                if (scheduler.my_rank != 0) {
                    SUCCESS_OR_DIE(
                        MPI_Win_lock(MPI_LOCK_SHARED, scheduler.my_rank, 0, scheduler.ens_window));
                    uint64_t result = scheduler.max_ensemble_id;
                    MPI_Win_unlock(scheduler.my_rank, scheduler.ens_window);
                    do {
                        ens_id = result + 1;
                        SUCCESS_OR_DIE(
                            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, scheduler.my_rank, 0, scheduler.ens_window));
                        scheduler.max_ensemble_id = result;
                        MPI_Win_unlock(scheduler.my_rank, scheduler.ens_window);
                        SUCCESS_OR_DIE(
                        MPI_Compare_and_swap(
                            &ens_id,                        // origin
                            &scheduler.max_ensemble_id,   // compare
                            &result,                        // result
                            MPI_UINT64_T,                     // datatype
                            0,                              // target rank
                            0,                              // target displ
                            scheduler.ens_window));
                    } while (result != scheduler.max_ensemble_id);
                } else if (input.simulation_mode != limited_time_ensembles) {
                    SUCCESS_OR_DIE(
                        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, scheduler.my_rank, 0, scheduler.ens_window));
                    ens_id = scheduler.max_ensemble_id + 1;
                    scheduler.max_ensemble_id = ens_id;
                    MPI_Win_unlock(scheduler.my_rank, scheduler.ens_window);
                } else {
                    ens_id = scheduler.max_ensemble_id;
                }
                if (ens_id < cc.n_ensembles) {
                    for (uint32_t i=1; i < s.n_members; ++i) {
                        const uint64_t total_members = s.n_members;
                        uint64_t duration = s.limit_duration() / (cc.dt_prime * cc.num_sub_steps);
                        if (s.limit_duration() + time_old > (cc.num_steps - cc.done_steps)*cc.dt_prime) {
                            duration = (cc.num_steps - cc.done_steps) + 1
                                - (time_old/(cc.dt_prime * cc.num_sub_steps));
                        }
                        cc.checkpoint_steps = abs(time_old)/cc.dt_prime;
                        if (input.simulation_mode == limited_time_ensembles) {
                            std::string throw_away = "";
                            s.activated = true;
                            s.perturb(cc, ref_quant, input, throw_away);
                        }
                        checkpoint_t checkpoint(
                            cc,
                            y_single_old,
                            segments,
                            input,
                            time_old,
                            i,
                            ens_id,
                            total_members,
                            duration);
                        scheduler.send_new_task(checkpoint);
                        if (input.simulation_mode == limited_time_ensembles) {
                            s.reset_variables(cc);
                        }
                    }
                }
            }
        }
    }
}

template<class float_t>
void finish_last_step(
    std::vector<float_t> &y_single_new,
    const reference_quantities_t &ref_quant,
    const model_constants_t<float_t> &cc) {
    float_t T_prime = y_single_new[T_idx]*ref_quant.Tref;
    float_t p_prime = y_single_new[p_idx]*ref_quant.pref;
    float_t qv_prime = y_single_new[qv_idx]*ref_quant.qref;
    float_t qc_prime = y_single_new[qc_idx]*ref_quant.qref;

    // Large gradients can be avoided by using the following line somehow?
    // There has been a single datapoint amongst 90 million where this line is
    // needed
    y_single_new[S_idx] = convert_qv_to_S(
            p_prime,
            T_prime,
            qv_prime,
            get_at(cc.constants, Cons_idx::p_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_a),
            get_at(cc.constants, Cons_idx::T_sat_low_temp),
            get_at(cc.constants, Cons_idx::p_sat_const_b),
            get_at(cc.constants, Cons_idx::Epsilon));
#ifdef TRACE_SAT
    if (trace)
        std::cout << "traj: " << cc.traj_id << ", S (last step before saturation adjustment): "
                  << y_single_new[S_idx] << "\n";
#endif
#ifdef TRACE_ENV
    if (trace)
        std::cout << "traj: " << cc.traj_id << ", before sat ad S " << y_single_new[S_idx];
#endif
    std::vector<float_t> res(num_comp);
    for (auto& r : res) r = 0;
    saturation_adjust(
        T_prime,
        p_prime,
        qv_prime,
        qc_prime,
        res,
        cc);
    y_single_new[qv_idx] += res[qv_idx]/ref_quant.qref;
    y_single_new[qv_idx] = (y_single_new[qv_idx] < 0) ? 0 : y_single_new[qv_idx];
    y_single_new[qc_idx] += res[qc_idx]/ref_quant.qref;
    // values of 1e-50 are likely which can lead to
    // negative qc due to rounding errors with ref_quant.qref.
    y_single_new[qc_idx] = (y_single_new[qc_idx] < 0) ? 0 : y_single_new[qc_idx];
    y_single_new[T_idx] += res[T_idx]/ref_quant.Tref;
    y_single_new[lat_cool_idx] += res[lat_cool_idx];
    y_single_new[lat_heat_idx] += res[lat_heat_idx];
    p_prime = y_single_new[p_idx]*ref_quant.pref;
    T_prime = y_single_new[T_idx]*ref_quant.Tref;
    qv_prime = y_single_new[qv_idx]*ref_quant.qref;
    y_single_new[S_idx] = convert_qv_to_S(
        p_prime,
        T_prime,
        qv_prime,
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b),
        get_at(cc.constants, Cons_idx::Epsilon));
#ifdef TRACE_SAT
    if (trace)
        std::cout << "traj: " << cc.traj_id << ", S (last step after saturation adjustment): "
                  << y_single_new[S_idx] << "\n";
#endif
#ifdef TRACE_ENV
    if (trace)
        std::cout << "traj: " << cc.traj_id << ", sat ad S " << y_single_new[S_idx]
            << "\nsat ad T " << y_single_new[T_idx]*ref_quant.Tref
            << ", QC: " << y_single_new[qc_idx]*ref_quant.qref << "\n";
#endif
}


template<class float_t>
void determine_phase_start(
    const std::vector<float_t> &y_single_old,
    bool &previous_step_with_ice,
    bool &previous_step_with_warm) {

    previous_step_with_ice = false;
#if defined(RK4ICE)
    if (y_single_old[qi_idx] > ice_q_phase_threshold || y_single_old[Ni_idx] > ice_n_phase_threshold
        || y_single_old[qg_idx] > ice_q_phase_threshold || y_single_old[Ng_idx] > ice_n_phase_threshold
        || y_single_old[qs_idx] > ice_q_phase_threshold || y_single_old[Ns_idx] > ice_n_phase_threshold
        || y_single_old[qh_idx] > ice_q_phase_threshold || y_single_old[Nh_idx]  > ice_n_phase_threshold) {
        previous_step_with_ice = true;
    }
#endif
    previous_step_with_warm = false;
    if (y_single_old[qc_idx] > warm_q_phase_threshold || y_single_old[qr_idx] > warm_q_phase_threshold
        || y_single_old[Nc_idx] > warm_n_phase_threshold || y_single_old[Nr_idx] > warm_n_phase_threshold ) {
        previous_step_with_warm = true;
    }
}


void run_substeps(
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const uint32_t &t,
    model_constants_t<codi::RealReverse> &cc,
    std::vector<codi::RealReverse> &y_single_old,
    std::vector<codi::RealReverse> &inflow,
    std::vector<codi::RealReverse> &y_single_new,
    netcdf_reader_t &netcdf_reader,
    std::vector< std::array<double, num_par > > &y_diff,
    output_handle_t &out_handler,
    const uint32_t &sub_start,
    std::vector<segment_t> &segments,
    ProgressBar &pbar,
    task_scheduler_t &scheduler,
    const double delay_out_time = 0) {

    codi::RealReverse::Tape& tape = codi::RealReverse::getTape();
    double time_old, time_new;
    for (uint32_t sub=sub_start; sub <= cc.num_sub_steps; ++sub) {
#if defined(TRACE_SAT) || defined(TRACE_QR) || defined(TRACE_QV) || defined(TRACE_QC) || defined(TRACE_QI) \
    || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH) || defined(TRACE_TIME)
        substep_trace(sub, t, cc, input, ref_quant, y_single_old, inflow);
#endif
        bool last_step = (((sub+1 + t*cc.num_sub_steps) >= ((t+1)*cc.num_sub_steps + 1))
            || (sub == cc.num_sub_steps));
#if defined(RK4_ONE_MOMENT)
        // Set the coefficients from the last timestep and from
        // the input files
        // *Should* only be necessary when parameters from the
        // trajectory are used as start point
        cc.setCoefficients(y_single_old, ref_quant);
#endif
        // CODIPACK
        tape.setActive();
#if defined(FLUX) && !defined(WCB)
        //     Add the inflow
        y_single_old[qr_idx] += inflow[qr_in_idx]/cc.num_sub_steps;
        y_single_old[Nr_idx] += inflow[Nr_in_idx]/cc.num_sub_steps;
    #if defined(RK4ICE)
        y_single_old[qi_idx] += inflow[qi_in_idx]/cc.num_sub_steps;
        y_single_old[qs_idx] += inflow[qs_in_idx]/cc.num_sub_steps;
        y_single_old[qg_idx] += inflow[qg_in_idx]/cc.num_sub_steps;
        y_single_old[Ni_idx] += inflow[Ni_in_idx]/cc.num_sub_steps;
        y_single_old[Ns_idx] += inflow[Ns_in_idx]/cc.num_sub_steps;
        y_single_old[Ng_idx] += inflow[Ng_in_idx]/cc.num_sub_steps;
    #if defined(B_EIGHT)
        y_single_old[qh_idx] += inflow[qh_in_idx]/cc.num_sub_steps;
        y_single_old[Nh_idx] += inflow[Nh_in_idx]/cc.num_sub_steps;
    #endif
    #endif
#endif
#if defined MET3D && defined TURBULENCE
        y_single_old[qv_idx] += inflow[qv_in_idx]/cc.num_sub_steps;
#endif
        if (!input.track_initial_cond)
            cc.register_input(tape);
        // We always have to set the dependent model constants to get
        // the correct impact of the parameters on which those depend on.
        cc.setup_dependent_model_constants();
        if (sub == 1) {
            codi::RealReverse p_prime = y_single_old[p_idx]*ref_quant.pref;
            codi::RealReverse T_prime = y_single_old[T_idx]*ref_quant.Tref;
            codi::RealReverse qv_prime = y_single_old[qv_idx]*ref_quant.qref;
            y_single_old[S_idx] = convert_qv_to_S(
                p_prime,
                T_prime,
                qv_prime,
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                get_at(cc.constants, Cons_idx::Epsilon));
#ifdef TRACE_SAT
            if (trace)
                std::cout << "traj: " << cc.traj_id << ", S (sub==1): " << y_single_old[S_idx] << "\n";
#endif
        }
#ifdef DEVELOP
        std::cout << "run RK4\n";
#endif
//////////////// Add any different scheme and model here
#if defined(RK4) || defined(RK4_ONE_MOMENT) || defined(OTHER)
        // Not implemented
#elif defined(EXPLICIT_EULER)
        euler_step(y_single_new, y_single_old, ref_quant, cc,
            input.fixed_iteration);
#elif defined(RK4ICE)
        RK4_step_2_sb_ice(y_single_new, y_single_old, ref_quant, cc,
            input.fixed_iteration);
#endif
#ifdef DEVELOP
        std::cout << "done with RK4\n";
#endif
#ifndef IN_SAT_ADJ
        if (last_step) {
            finish_last_step(y_single_new, ref_quant, cc);
        }
#endif
#ifdef DEVELOP
        std::cout << "get gradients\n";
#endif
        if ( (input.simulation_mode == create_train_set && cc.ensemble_id == 0)
             || (input.simulation_mode == limited_time_ensembles && cc.traj_id == 0)
             || (input.simulation_mode != create_train_set && input.simulation_mode != limited_time_ensembles) )
            cc.get_gradients(y_single_new, y_diff, tape, ref_quant);
#ifdef DEVELOP
        std::cout << "got gradients\n";
#endif
        // Time update
        time_new = (sub + (t+cc.done_steps)*cc.num_sub_steps)*cc.dt;
        if (time_new >= delay_out_time) {
#ifdef DEVELOP
            std::cout << "process_step\n";
#endif
            bool previous_step_with_ice, previous_step_with_warm;
            determine_phase_start(y_single_old, previous_step_with_ice, previous_step_with_warm);
            // TODO(mahieron): what if delay_out_time is not a multiple of dt_prime?
            out_handler.process_step(cc, netcdf_reader, y_single_new, y_diff,
                sub, t,
                input.write_index,
                input.snapshot_index,
                last_step, ref_quant, previous_step_with_warm, previous_step_with_ice);
        }
        // Interchange old and new for next step
        time_old = time_new;
        y_single_old.swap(y_single_new);
        if (time_new != cc.t_end_prime) {
            // Check if parameter shall be perturbed
#ifdef DEVELOP
            std::cout << "check parameters\n";
#endif
            parameter_check(segments, cc, time_old, y_diff,
                y_single_old, input, ref_quant, scheduler);
        }

        pbar.progress();

        // In case that the sub timestep%timestep is not 0
        if (last_step)
            break;
    }  // End substep
#if defined(TRACE_COMM_DEBUG) || defined(DEVELOP)
    std::cout << "after substep time \n";
#endif
}

void run_substeps(
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const uint32_t &t,
    model_constants_t<codi::RealForwardVec<num_par_init> > &cc,
    std::vector<codi::RealForwardVec<num_par_init> > &y_single_old,
    std::vector<codi::RealForwardVec<num_par_init> > &inflow,
    std::vector<codi::RealForwardVec<num_par_init> > &y_single_new,
    netcdf_reader_t &netcdf_reader,
    std::vector< std::array<double, num_par > > &y_diff,
    output_handle_t &out_handler,
    const uint32_t &sub_start,
    std::vector<segment_t> &segments,
    ProgressBar &pbar,
    task_scheduler_t &scheduler,
    const double delay_out_time = 0) {

    double time_old, time_new;
    for (uint32_t sub=sub_start; sub <= cc.num_sub_steps; ++sub) {
#if defined(TRACE_SAT) || defined(TRACE_QR) || defined(TRACE_QV) || defined(TRACE_QC) || defined(TRACE_QI) \
    || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH) || defined(TRACE_TIME)
        substep_trace(sub, t, cc, input, ref_quant, y_single_old, inflow);
#endif
        bool last_step = (((sub+1 + t*cc.num_sub_steps) >= ((t+1)*cc.num_sub_steps + 1))
            || (sub == cc.num_sub_steps));
#if defined(RK4_ONE_MOMENT)
        // Set the coefficients from the last timestep and from
        // the input files
        // *Should* only be necessary when parameters from the
        // trajectory are used as start point
        cc.setCoefficients(y_single_old, ref_quant);
#endif
#if defined(FLUX) && !defined(WCB)
        //     Add the inflow
        y_single_old[qr_idx] += inflow[qr_in_idx]/cc.num_sub_steps;
        y_single_old[Nr_idx] += inflow[Nr_in_idx]/cc.num_sub_steps;
    #if defined(RK4ICE)
        y_single_old[qi_idx] += inflow[qi_in_idx]/cc.num_sub_steps;
        y_single_old[qs_idx] += inflow[qs_in_idx]/cc.num_sub_steps;
        y_single_old[qg_idx] += inflow[qg_in_idx]/cc.num_sub_steps;
        y_single_old[Ni_idx] += inflow[Ni_in_idx]/cc.num_sub_steps;
        y_single_old[Ns_idx] += inflow[Ns_in_idx]/cc.num_sub_steps;
        y_single_old[Ng_idx] += inflow[Ng_in_idx]/cc.num_sub_steps;
    #if defined(B_EIGHT)
        y_single_old[qh_idx] += inflow[qh_in_idx]/cc.num_sub_steps;
        y_single_old[Nh_idx] += inflow[Nh_in_idx]/cc.num_sub_steps;
    #endif
    #endif
#endif
#if defined MET3D && defined TURBULENCE
        y_single_old[qv_idx] += inflow[qv_in_idx]/cc.num_sub_steps;
#endif
        // We always have to set the dependent model constants to get
        // the correct impact of the parameters on which those depend on.
        cc.setup_dependent_model_constants();
        if (sub == 1) {
            codi::RealForwardVec<num_par_init> p_prime = y_single_old[p_idx]*ref_quant.pref;
            codi::RealForwardVec<num_par_init> T_prime = y_single_old[T_idx]*ref_quant.Tref;
            codi::RealForwardVec<num_par_init> qv_prime = y_single_old[qv_idx]*ref_quant.qref;
            y_single_old[S_idx] = convert_qv_to_S(
                p_prime,
                T_prime,
                qv_prime,
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                get_at(cc.constants, Cons_idx::Epsilon));
#ifdef TRACE_SAT
            if (trace)
                std::cout << "traj: " << cc.traj_id << ", S (sub==1): " << y_single_old[S_idx] << "\n";
#endif
        }
//////////////// Add any different scheme and model here
#if defined(RK4) || defined(RK4_ONE_MOMENT) || defined(OTHER)
        // Not implemented
#elif defined(RK4NOICE)
        // Not implemented
#elif defined(RK4ICE)
        RK4_step_2_sb_ice(y_single_new, y_single_old, ref_quant, cc,
            input.fixed_iteration);
#endif
#ifndef IN_SAT_ADJ
        if (last_step) {
            finish_last_step(y_single_new, ref_quant, cc);
        }
#endif
        if ( (input.simulation_mode == create_train_set && cc.ensemble_id == 0)
             || (input.simulation_mode == limited_time_ensembles && cc.traj_id == 0)
             || (input.simulation_mode != create_train_set && input.simulation_mode != limited_time_ensembles) )
            cc.get_gradients(y_single_new, y_diff, ref_quant);

        // Time update
        time_new = (sub + (t+cc.done_steps)*cc.num_sub_steps)*cc.dt;
        if (time_new >= delay_out_time) {
            // TODO(mahieron): what if delay_out_time is not a multiple of dt_prime?
            bool previous_step_with_ice, previous_step_with_warm;
            determine_phase_start(y_single_old, previous_step_with_ice, previous_step_with_warm);

            out_handler.process_step(cc, netcdf_reader, y_single_new, y_diff,
                sub, t,
                input.write_index,
                input.snapshot_index,
                last_step, ref_quant, previous_step_with_warm, previous_step_with_ice);
        }
        // Interchange old and new for next step
        time_old = time_new;
        y_single_old.swap(y_single_new);
        if (time_new != cc.t_end_prime) {
            // Check if parameter shall be perturbed
            parameter_check(segments, cc, time_old, y_diff,
                y_single_old, input, ref_quant, scheduler);
        }
        pbar.progress();

        // In case that the sub timestep%timestep is not 0
        if (last_step)
            break;
    }  // End substep
}

template<class float_t>
int run_simulation(
    const int &rank,
    model_constants_t<float_t> &cc,
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const global_args_t &global_args,
    std::vector<float_t> &y_single_old,
    std::vector< std::array<double, num_par > > &y_diff,
    std::vector<float_t> &y_single_new,
    std::vector<float_t> &inflow,
    output_handle_t &out_handler,
    std::vector<segment_t> &segments,
    task_scheduler_t &scheduler,
    netcdf_reader_t &netcdf_reader,
    const double delay_out_time = 0,
    const uint32_t start_step = 0) {
    // force any process that is not root to disable pbar
    const uint64_t progress_index = (rank != 0) ? 0 : input.progress_index;
    ProgressBar pbar = ProgressBar((cc.num_sub_steps)*cc.num_steps,
        progress_index, "simulation step", std::cout);
    return run_simulation(
        cc,
        input,
        ref_quant,
        global_args,
        y_single_old,
        y_diff,
        y_single_new,
        inflow,
        out_handler,
        segments,
        scheduler,
        netcdf_reader,
        pbar,
        delay_out_time,
        start_step);
}

template<class float_t>
int run_simulation(
    model_constants_t<float_t> &cc,
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const global_args_t &global_args,
    std::vector<float_t> &y_single_old,
    std::vector< std::array<double, num_par > > &y_diff,
    std::vector<float_t> &y_single_new,
    std::vector<float_t> &inflow,
    output_handle_t &out_handler,
    std::vector<segment_t> &segments,
    task_scheduler_t &scheduler,
    netcdf_reader_t &netcdf_reader,
    ProgressBar &pbar,
    const double delay_out_time = 0,
    const uint32_t start_step = 0,
    const double finish_progress = true) {

    uint32_t sub_start = 1;
    if (global_args.checkpoint_flag && std::fmod(input.current_time, cc.dt_prime) != 0
        && !std::isnan(input.current_time))
        sub_start = std::fmod(input.current_time, cc.dt_prime)
                / (cc.dt_prime/(cc.num_sub_steps));
    if (input.track_initial_cond) cc.register_input();
    sediment_n_total = 0;
    sediment_q_total = 0;
    // Loop over every timestep that is usually fixed to 20 s
    for (uint32_t t=start_step; t < cc.num_steps - cc.done_steps; ++t) {
        if (netcdf_reader.read_buffer(cc, ref_quant, y_single_old,
            inflow, t, global_args.checkpoint_flag, input.start_over_env) != SUCCESS) {
            // If the input file consists of (multiple) NaNs, we do not
            // need a simulation.
#ifdef DEVELOP
            std::cout << cc.rank << " flush buffer\n";
#endif
            out_handler.flush_buffer(cc);
            break;
        }
#if defined(TRACE_COMM_DEBUG) || defined(DEVELOP)
        std::cout << cc.rank << ", sim " << t << " / " << cc.num_steps - cc.done_steps - 1 << "\n" << std::flush;
#endif
        // Iterate over each substep
        run_substeps(input, ref_quant, t, cc, y_single_old,
            inflow, y_single_new, netcdf_reader, y_diff, out_handler,
            sub_start,
            segments, pbar, scheduler,
            delay_out_time);
#if defined(TRACE_COMM_DEBUG) || defined(DEVELOP)
       std::cout << cc.rank << ", sim done " << t << " / " << cc.num_steps - cc.done_steps - 1 << "\n";
#endif
#ifdef TRACE_QG
        if (trace)
            std::cout << "\nSediment total q: " << sediment_q_total
                    << "\nSediment total N: " << sediment_n_total << "\n" << std::flush;
        sediment_n_total = 0;
        sediment_q_total = 0;
#endif
        sub_start = 1;
        checkpoint_t throw_away;
        scheduler.send_task(throw_away, false);
    }
#if defined(TRACE_COMM_DEBUG) || defined(DEVELOP)
    std::cout << cc.rank << ", sim end. pbar.finish()\n";
#endif
    if (finish_progress)
        pbar.finish();

    return 0;
}


template<class float_t>
void busy_flush(
    model_constants_t<float_t> &cc,
    output_handle_t &out_handler) {
#ifdef COMPRESS_OUTPUT
    // With compression enabled we need to flush the data in a collective manner.
    // If the workload is not equally distributed (i.e. number of trajectories
    // is not a multiple of the number of processes), we need to force some
    // processes to access the flush routine without flushing anything.
    // rank 0 should always check if all processes are done.
    do {
    } while (out_handler.flush_buffer(cc, true));
#endif
}


template<class float_t>
void only_sensitivity_simulation(
    const int &rank,
    const int &n_processes,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    task_scheduler_t &scheduler,
    netcdf_reader_t &netcdf_reader) {

    std::vector<double> y_init(num_comp);
    checkpoint_t checkpoint;
    std::vector<segment_t> segments;
    std::vector< std::array<double, num_par > >  y_diff(num_comp + static_cast<uint32_t>(Init_cons_idx::n_items));

    bool already_loaded = false;
    std::vector<float_t> y_single_old(num_comp);
    std::vector<float_t> y_single_new(num_comp);
    std::vector<float_t> inflow(num_inflows);

    model_constants_t<float_t> cc = prepare_constants<float_t>(rank, input, global_args,
        ref_quant, segments, y_init, checkpoint);
    netcdf_reader.set_dims(input.INPUT_FILENAME.c_str(), cc, input.simulation_mode);

    // static scheduling with parallel read and write enabled
    setup_simulation_base(rank, input,
        global_args, ref_quant, cc,
        already_loaded, netcdf_reader);
#ifdef COMPRESS_OUTPUT
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
        input.INPUT_FILENAME, input.write_index,
        input.snapshot_index, rank, input.simulation_mode,
        input.track_initial_cond, n_processes, input.delay_time_store);
#else
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
        input.INPUT_FILENAME, input.write_index,
        input.snapshot_index, rank, input.simulation_mode,
        input.track_initial_cond, input.delay_time_store);
#endif
#ifdef DEVELOP
        std::cout << "out_handler done\n" << std::flush;
#endif
    scheduler.set_n_ensembles(netcdf_reader.n_ensembles);
    scheduler.set_n_trajectories(netcdf_reader.n_trajectories);
#ifdef TRACE_COMM
    std::cout << "rank " << rank << ", n_ens: " << netcdf_reader.n_ensembles
        << ", n_trajs: " << netcdf_reader.n_trajectories << "\n";
#endif
    const uint64_t progress_index = (rank != 0) ? 0 : input.progress_index;
    int sims_for_r0 = (netcdf_reader.n_ensembles*netcdf_reader.n_trajectories + n_processes-1)/n_processes;

    ProgressBar pbar = ProgressBar(cc.num_sub_steps*cc.num_steps*sims_for_r0,
        progress_index, "simulation step", std::cout);

    if (scheduler.receive_task(checkpoint)) {
        cc.traj_id = scheduler.current_traj;
        cc.ensemble_id = scheduler.current_ens;
#ifdef TRACE_COMM
        std::cout << "rank " << rank << ", received traj: " << scheduler.current_traj
            << ", ens: " << scheduler.current_ens << "\n";
#endif
        netcdf_reader.read_initial_values(y_init, ref_quant, cc,
            global_args.checkpoint_flag, cc.traj_id, cc.ensemble_id);
        // Set "old" values as temporary holder of values.
        for (int ii = 0 ; ii < num_comp ; ii++)
            y_single_old[ii] = y_init[ii];
#ifdef TRACE_COMM
        std::cout << "rank " << rank << ", init pressure: "
            << y_init[p_idx]*ref_quant.pref << "\n";
#endif
        out_handler.reset(scheduler.current_traj, scheduler.current_ens);
        int sim_counter = 1;
        cc.rank = rank;
        // run simulation
        SUCCESS_OR_DIE(run_simulation(cc, input, ref_quant,
            global_args, y_single_old, y_diff, y_single_new, inflow,
            out_handler, segments, scheduler, netcdf_reader, pbar, input.delay_time_store, 0, false));

        while (scheduler.receive_task(checkpoint)) {
            cc.traj_id = scheduler.current_traj;
            cc.ensemble_id = scheduler.current_ens;
            pbar.set_current_step(sim_counter*cc.num_sub_steps*cc.num_steps-1);
            pbar.progress();
            sim_counter++;

            setup_simulation_base(rank, input,
                global_args, ref_quant, cc,
                already_loaded, netcdf_reader);
            netcdf_reader.read_initial_values(y_init, ref_quant, cc,
                global_args.checkpoint_flag, scheduler.current_traj, scheduler.current_ens);
            out_handler.reset(scheduler.current_traj, scheduler.current_ens);

#ifdef TRACE_COMM
            std::cout << "rank " << rank << ", received traj: " << scheduler.current_traj
                << ", ens: " << scheduler.current_ens << "\n";
#endif
            // Set "old" values as temporary holder of values.
            for (int ii = 0 ; ii < num_comp ; ii++)
                y_single_old[ii] = y_init[ii];
#ifdef TRACE_COMM
            std::cout << "rank " << rank << ", init pressure: "
                << y_init[p_idx]*ref_quant.pref << "\n";
#endif
            cc.rank = rank;
            // run simulation
            SUCCESS_OR_DIE(run_simulation(cc, input, ref_quant,
                global_args, y_single_old, y_diff, y_single_new, inflow,
                out_handler, segments, scheduler, netcdf_reader, pbar, input.delay_time_store, 0, false));
        }
    }
#ifdef DEVELOP
    std::cout << "\n" << rank << " before busy_flush\n";
#endif
    busy_flush(cc, out_handler);
#ifdef DEVELOP
    std::cout << "\n" << rank << " after busy_flush\n";
#endif
    pbar.finish();
}

template<class float_t>
void limited_time_ensemble_simulation(
    const int &rank,
    const int &n_processes,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    task_scheduler_t &scheduler,
    netcdf_reader_t &netcdf_reader) {

    std::vector<double> y_init(num_comp);
    checkpoint_t checkpoint;
    std::vector<segment_t> segments;
    std::vector< std::array<double, num_par > >  y_diff(num_comp + static_cast<uint32_t>(Init_cons_idx::n_items));

    bool already_loaded = false;
    std::vector<float_t> y_single_old(num_comp);
    std::vector<float_t> y_single_new(num_comp);
    std::vector<float_t> inflow(num_inflows);

    model_constants_t<float_t> cc = prepare_constants<float_t>(rank, input, global_args,
        ref_quant, segments, y_init, checkpoint);
    netcdf_reader.set_dims(input.INPUT_FILENAME.c_str(), cc, input.simulation_mode);

    // static scheduling with parallel read and write enabled
    // The output is based on the ensemble configuration file.
    setup_simulation_base(rank, input,
        global_args, ref_quant, cc,
        already_loaded, netcdf_reader);
#ifdef COMPRESS_OUTPUT
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
        input.INPUT_FILENAME, input.write_index,
        input.snapshot_index, rank, input.simulation_mode,
        input.track_initial_cond, n_processes, input.delay_time_store);
#else
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
        input.INPUT_FILENAME, input.write_index,
        input.snapshot_index, rank, input.simulation_mode,
        input.track_initial_cond, input.delay_time_store);
#endif
    // The progressbar here is just an estimate since the members
    // are distributed dynamically
    const uint64_t progress_index = (rank != 0) ? 0 : input.progress_index;
    uint64_t sims_for_r0 = (netcdf_reader.n_ensembles + n_processes-1)/n_processes;

    uint64_t steps_members = 0;
    for (auto &s : segments) {
        const uint64_t n_repeats = (cc.num_sub_steps * cc.num_steps * cc.dt_prime - s.value+1) / s.value - 1;

        steps_members += n_repeats * (s.n_members + n_processes - 2)/(n_processes/2)
            * s.n_segments * s.duration/cc.dt_prime;
    }
    ProgressBar pbar = ProgressBar(
        cc.num_sub_steps*cc.num_steps * sims_for_r0 + steps_members,
        progress_index, "simulation step", std::cout);
    SUCCESS_OR_DIE(MPI_Win_lock_all(0, scheduler.free_window));
    if (rank == 0) {
        scheduler.set_n_ensembles(1);
        scheduler.set_n_trajectories(segments[0].n_members);
        netcdf_reader.read_initial_values(y_init, ref_quant, cc,
            global_args.checkpoint_flag, input.traj, input.ensemble);
#ifdef DEVELOP
        std::cout << rank << " setting y for " << num_comp << " variables. sizes: "
                  << y_single_old.size() << ", " << y_init.size() << "\n";
#endif
        // Set "old" values as temporary holder of values.
        for (int ii = 0 ; ii < num_comp ; ii++)
            y_single_old[ii] = y_init[ii];
#ifdef DEVELOP
        std::cout << rank << " resetting out_handler\n";
#endif
        out_handler.reset(scheduler.current_traj, scheduler.current_ens);
#ifdef DEVELOP
        std::cout << rank << " run simulation \n";
#endif
        SUCCESS_OR_DIE(run_simulation(cc, input, ref_quant,
            global_args, y_single_old, y_diff, y_single_new, inflow,
            out_handler, segments, scheduler, netcdf_reader, pbar, input.delay_time_store, 0, false));
    }
    uint64_t pbar_counter =  cc.num_sub_steps*cc.num_steps-1;

    while (scheduler.receive_task(checkpoint)) {
#ifdef DEVELOP
        std::cout << " rank " << rank << " received task \n";
#endif
        global_args.checkpoint_flag = true;

        setup_simulation(rank, input,
            global_args, ref_quant, segments, cc, y_init,
            checkpoint, out_handler, already_loaded, netcdf_reader);

        // the scheduler sets already the correct values at out_handler
        scheduler.current_traj = cc.traj_id;
        scheduler.current_ens = cc.ensemble_id;
        netcdf_reader.read_initial_values(y_init, ref_quant, cc,
            global_args.checkpoint_flag, input.traj, input.ensemble);

#ifdef TRACE_COMM
        std::cout << "rank " << rank << ", received traj: " << scheduler.current_traj
            << ", ens: " << scheduler.current_ens << "\n";
#endif

        // Set "old" values as temporary holder of values.
        for (int ii = 0 ; ii < num_comp ; ii++)
            y_single_old[ii] = y_init[ii];

        // Need to find the number of flushed snapshots; already done in load_checkpoint
        const uint64_t out_timestep = (input.current_time - input.delay_time_store)
            / (cc.dt_prime * cc.num_sub_steps);
        out_handler.reset(scheduler.current_traj, scheduler.current_ens, out_timestep);

        // buffer the initial values for this trajectory
        // -> one step is already done this way
        out_handler.buffer(
            cc,
            netcdf_reader,
            y_single_old,
            y_diff,
            0,
            0,
            ref_quant,
            input.snapshot_index,
            false,
            false);

        pbar.set_current_step(pbar_counter);
        pbar.progress();
#ifdef TRACE_COMM
        std::cout << "rank " << rank << ", run simulation\n";
#endif
        // run simulation
        SUCCESS_OR_DIE(run_simulation(cc, input, ref_quant,
            global_args, y_single_old, y_diff, y_single_new, inflow,
            out_handler, segments, scheduler, netcdf_reader, pbar, 0, 1, false));
        pbar_counter += (cc.num_steps - cc.done_steps)*cc.num_sub_steps;
#ifdef TRACE_COMM
        std::cout << "rank " << rank << ", simulation done\n";
#endif
    }
#ifdef TRACE_COMM
        std::cout << "rank " << rank << " busy\n";
#endif
    busy_flush(cc, out_handler);
    pbar.finish();
#ifdef TRACE_COMM
        std::cout << "rank " << rank << " all done\n";
#endif
    SUCCESS_OR_DIE(MPI_Win_unlock_all(scheduler.free_window));
}

template<class float_t>
void create_set_simulation(
        const int &rank,
        const int &n_processes,
        input_parameters_t &input,
        global_args_t &global_args,
        reference_quantities_t &ref_quant,
        task_scheduler_t &scheduler,
        netcdf_reader_t &netcdf_reader) {
    std::vector<double> y_init(num_comp);
    checkpoint_t checkpoint;
    std::vector<segment_t> segments;
    std::vector< std::array<double, num_par > >  y_diff(num_comp + static_cast<uint32_t>(Init_cons_idx::n_items));

    bool already_loaded = false;
    std::vector<float_t> y_single_old(num_comp);
    std::vector<float_t> y_single_new(num_comp);
    std::vector<float_t> inflow(num_inflows);

    model_constants_t<float_t> cc = prepare_constants<float_t>(rank, input, global_args,
                                                               ref_quant, segments, y_init, checkpoint);
    netcdf_reader.set_dims(input.INPUT_FILENAME.c_str(), cc, input.simulation_mode);

    // static scheduling with parallel read and write enabled
    // The output is based on the ensemble configuration file.
    setup_simulation_base(rank, input,
                          global_args, ref_quant, cc,
                          already_loaded, netcdf_reader);
    cc.n_ensembles = segments.size() + 1;
#ifdef COMPRESS_OUTPUT
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
                                input.INPUT_FILENAME, input.write_index,
                                input.snapshot_index, rank, input.simulation_mode,
                                input.track_initial_cond, n_processes, input.delay_time_store, segments);
#else
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
        input.INPUT_FILENAME, input.write_index,
        input.snapshot_index, rank, input.simulation_mode,
        input.track_initial_cond, input.delay_time_store, segments);
#endif
    const uint64_t progress_index = (rank != 0) ? 0 : input.progress_index;
    uint64_t sims_for_r0 = ((segments.size()+1)*netcdf_reader.n_trajectories + n_processes-1)/n_processes;

    ProgressBar pbar = ProgressBar(cc.num_sub_steps*cc.num_steps*sims_for_r0,
                                   progress_index, "simulation step", std::cout);
    scheduler.set_n_ensembles(segments.size() + 1);
    scheduler.set_n_trajectories(netcdf_reader.n_trajectories);

    uint64_t old_ensemble_id = 0;
    if (scheduler.receive_task(checkpoint)) {
        cc.traj_id = scheduler.current_traj;
        cc.ensemble_id = scheduler.current_ens;
#if defined(DEBUG_SEG) || defined(TRACE_COMM)
        std::cout << "rank " << rank << ", received traj: " << scheduler.current_traj
            << ", ens: " << scheduler.current_ens << "\n";
        std::cout << "rank " << rank << " pert " << cc.perturbed_idx.size() << "\n";
#endif
        netcdf_reader.read_initial_values(y_init, ref_quant, cc,
                                          global_args.checkpoint_flag, cc.traj_id, 0);
        // Set "old" values as temporary holder of values.
        for (int ii = 0 ; ii < num_comp ; ii++)
            y_single_old[ii] = y_init[ii];
#ifdef TRACE_COMM
        std::cout << "rank " << rank << ", init pressure: "
            << y_init[p_idx]*ref_quant.pref << "\n";
#endif
        out_handler.reset(scheduler.current_traj, scheduler.current_ens);
        int sim_counter = 1;
        cc.rank = rank;
        if (cc.ensemble_id != 0) {
            std::string descr;
            segments[cc.ensemble_id - 1].perturb(cc, ref_quant, input, descr);
            old_ensemble_id = cc.ensemble_id;
        }
#ifdef DEBUG_SEG
        std::cout << "rank " << rank << " run simulation\n";
        std::cout << "rank " << rank << " pert2 " << cc.perturbed_idx.size() << "\n";
        std::cout << "rank " << rank << " pert3 " << cc.perturbed_idx.size() << "\n";
#endif
        // run simulation
        SUCCESS_OR_DIE(run_simulation(cc, input, ref_quant,
                                      global_args, y_single_old, y_diff, y_single_new, inflow,
                                      out_handler, segments, scheduler, netcdf_reader, pbar,
                                      input.delay_time_store, 0, false));

        while (scheduler.receive_task(checkpoint)) {
            cc.traj_id = scheduler.current_traj;
            cc.ensemble_id = scheduler.current_ens;
            pbar.set_current_step(sim_counter*cc.num_sub_steps*cc.num_steps-1);
            pbar.progress();
            sim_counter++;
#ifdef DEBUG_SEG
            std::cout << "rank " << rank << " before reset " << cc.perturbed_idx.size() << "\n";
#endif
            if (cc.ensemble_id != old_ensemble_id && old_ensemble_id != 0) {
                segments[old_ensemble_id - 1].reset_variables(cc);
            }
#ifdef DEBUG_SEG
            std::cout << "rank " << rank << " after reset " << cc.perturbed_idx.size() << "\n";
#endif
            setup_simulation_base(rank, input,
                                  global_args, ref_quant, cc,
                                  already_loaded, netcdf_reader);
#ifdef DEBUG_SEG
            std::cout << "rank " << rank << " after setup_base " << cc.perturbed_idx.size() << "\n";
#endif
            netcdf_reader.read_initial_values(y_init, ref_quant, cc,
                                              global_args.checkpoint_flag,
                                              scheduler.current_traj, 0);
#ifdef DEBUG_SEG
            std::cout << "rank " << rank << " after read init " << cc.perturbed_idx.size() << "\n";
#endif
            out_handler.reset(scheduler.current_traj, scheduler.current_ens);
            if (cc.ensemble_id != old_ensemble_id && cc.ensemble_id != 0) {
                std::string descr;
#ifdef DEBUG_SEG
                std::cout << "rank " << rank << " before pert " << cc.perturbed_idx.size() << "\n";
#endif
                segments[cc.ensemble_id - 1].perturb(cc, ref_quant, input, descr);
#ifdef DEBUG_SEG
                std::cout << "rank " << rank << " after pert " << cc.perturbed_idx.size() << "\n";
#endif
                old_ensemble_id = cc.ensemble_id;
            }
#if defined(DEBUG_SEG) || defined(TRACE_COMM)
            std::cout << "rank " << rank << ", received traj: " << scheduler.current_traj
                << ", ens: " << scheduler.current_ens << "\n";
#endif
            // Set "old" values as temporary holder of values.
            for (int ii = 0 ; ii < num_comp ; ii++)
                y_single_old[ii] = y_init[ii];
#ifdef TRACE_COMM
            std::cout << "rank " << rank << ", init pressure: "
                << y_init[p_idx]*ref_quant.pref << "\n";
#endif
            cc.rank = rank;
#ifdef DEBUG_SEG
            std::cout << "rank " << rank << " pert4 " << cc.perturbed_idx.size() << "\n";
            std::cout << "rank " << rank << " run simulation\n";
#endif
            // run simulation
            SUCCESS_OR_DIE(run_simulation(cc, input, ref_quant,
                                          global_args, y_single_old, y_diff, y_single_new, inflow,
                                          out_handler, segments, scheduler, netcdf_reader, pbar,
                                          input.delay_time_store, 0, false));
        }
    }
#ifdef DEVELOP
    std::cout << "\n" << rank << " before busy_flush\n";
#endif
    busy_flush(cc, out_handler);
#ifdef DEVELOP
    std::cout << "\n" << rank << " after busy_flush\n";
#endif
    pbar.finish();
}

template<class float_t>
void dynamic_ensemble_simulation(
    const int &rank,
    const int &n_processes,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    task_scheduler_t &scheduler,
    netcdf_reader_t &netcdf_reader) {

    std::vector<double> y_init(num_comp);
    checkpoint_t checkpoint;
    std::vector<segment_t> segments;
    std::vector< std::array<double, num_par > >  y_diff(num_comp + static_cast<uint32_t>(Init_cons_idx::n_items));

    bool already_loaded = false;
    std::vector<float_t> y_single_old(num_comp);
    std::vector<float_t> y_single_new(num_comp);
    std::vector<float_t> inflow(num_inflows);

    model_constants_t<float_t> cc = prepare_constants<float_t>(rank, input, global_args,
        ref_quant, segments, y_init, checkpoint);
#ifdef COMPRESS_OUTPUT
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
        input.INPUT_FILENAME, input.write_index,
        input.snapshot_index, rank, input.simulation_mode,
        input.track_initial_cond, n_processes);
#else
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
        input.INPUT_FILENAME, input.write_index,
        input.snapshot_index, rank, input.simulation_mode,
        input.track_initial_cond);
#endif
    netcdf_reader.set_dims(input.INPUT_FILENAME.c_str(), cc, input.simulation_mode);
    if (rank == 0) {
        setup_simulation(rank, input,
            global_args, ref_quant, segments, cc, y_init,
            checkpoint, out_handler, already_loaded, netcdf_reader);

        netcdf_reader.read_initial_values(y_init, ref_quant, cc,
            global_args.checkpoint_flag, input.traj, input.ensemble);

        // Set "old" values as temporary holder of values.
        for (int ii = 0 ; ii < num_comp ; ii++)
            y_single_old[ii] = y_init[ii];
        SUCCESS_OR_DIE(run_simulation(rank, cc, input, ref_quant,
            global_args, y_single_old, y_diff, y_single_new, inflow,
            out_handler, segments, scheduler, netcdf_reader));
    }

    while (scheduler.receive_task(checkpoint)) {
        global_args.checkpoint_flag = true;
        setup_simulation(rank, input,
            global_args, ref_quant, segments, cc, y_init,
            checkpoint, out_handler, already_loaded, netcdf_reader);

        netcdf_reader.read_initial_values(y_init, ref_quant, cc,
            global_args.checkpoint_flag, input.traj, input.ensemble);

        // Set "old" values as temporary holder of values.
        for (int ii = 0 ; ii < num_comp ; ii++)
            y_single_old[ii] = y_init[ii];
        // run simulation
        SUCCESS_OR_DIE(run_simulation(rank, cc, input, ref_quant,
            global_args, y_single_old, y_diff, y_single_new, inflow,
            out_handler, segments, scheduler, netcdf_reader));
    }
    busy_flush(cc, out_handler);
}


/**
 * Based on
 * Revisiting the latent heating contribution to foehn warming:
 * Lagrangian analysis of two foehn events over the Swiss Alps
 * DOI:10.1002/qj.2816
 *
 * Horizontal grid spacing of 2.2 km
 * Initial conditions from COMSO model simulation with horizontal grid spacing
 * of 7 km
 * 60 levels with mean spacing of 388 m
 *
 * Start at 2100 UTC on 3 March 2013 for the dry foehn for 102 h
 * Start at 0300 UTC on 14 May 2013 for the moist foehn  for 93 h
 *
 * temporal resolution of 20 s
 *
 * 14125 trajectories for the dry foehn event
 * 27491 trajectories for the moist foehn event
 *
 * important variables:
 * - air temperature T
 * - pressure P
 * - mixing ratio of water vapour, cloud droplets, rain, ice, graupel and snow
 * - velocity in three dimensions
 *
 */
int main(int argc, char** argv) {
    int rank;
    int n_processes;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
    n_processes = 1;
#endif
    input_parameters_t input;
    global_args_t global_args;
    reference_quantities_t ref_quant;
    parse_args(argc, argv, rank, input,
        global_args, ref_quant);

    netcdf_reader_t netcdf_reader(input.write_index);
    if (global_args.time_start_idx_flag) {
        netcdf_reader.start_time_idx_given = true;
        netcdf_reader.start_time_idx = input.start_time_idx;
    }
    task_scheduler_t scheduler(rank, n_processes, input.simulation_mode);

    if ((input.simulation_mode == grid_sensitivity)
        || (input.simulation_mode == trajectory_sensitivity)) {
        if (input.track_initial_cond) {
            only_sensitivity_simulation<codi::RealForwardVec<num_par_init> >(
                rank, n_processes, input, global_args,
                ref_quant, scheduler, netcdf_reader);
        } else {
            only_sensitivity_simulation<codi::RealReverse>(
                rank, n_processes, input, global_args,
                ref_quant, scheduler, netcdf_reader);
        }
    } else if (input.simulation_mode == limited_time_ensembles) {
        if (input.track_initial_cond) {
            limited_time_ensemble_simulation<codi::RealForwardVec<num_par_init> >(
                    rank, n_processes, input, global_args,
                    ref_quant, scheduler, netcdf_reader);
        } else {
            limited_time_ensemble_simulation<codi::RealReverse>(
                    rank, n_processes, input, global_args,
                    ref_quant, scheduler, netcdf_reader);
        }
    } else if (input.simulation_mode == create_train_set) {
        if (input.track_initial_cond) {
            create_set_simulation<codi::RealForwardVec<num_par_init> >(
                    rank, n_processes, input, global_args,
                    ref_quant, scheduler, netcdf_reader);
        } else {
            create_set_simulation<codi::RealReverse>(
                    rank, n_processes, input, global_args,
                    ref_quant, scheduler, netcdf_reader);
        }
    } else {   // dynamic scheduling with parallel read and write disabled
        if (input.track_initial_cond) {
            dynamic_ensemble_simulation<codi::RealForwardVec<num_par_init> >(
                rank, n_processes, input, global_args,
                ref_quant, scheduler, netcdf_reader);
        } else {
            dynamic_ensemble_simulation<codi::RealReverse>(
                rank, n_processes, input, global_args,
                ref_quant, scheduler, netcdf_reader);
        }
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
#ifdef SILENT_MODE
    exit(0);
#else
    if (rank == 0)
        std::cout << "-------------------FINISHED-------------------\n\n\n";
    exit(0);
#endif
}
