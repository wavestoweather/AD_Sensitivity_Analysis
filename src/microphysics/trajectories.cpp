#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <unistd.h>
#include <netcdf>
#include <cmath>
#include <fstream>
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


model_constants_t parse_args(
    const int &argc,
    char* const * argv,
    const int &rank,
    const int &n_processes,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    std::vector<segment_t> &segments,
    std::vector<double> &y_init,
    checkpoint_t &checkpoint) {
    SUCCESS_OR_DIE(global_args.parse_arguments(argc, argv, rank, n_processes));

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
    if (rank == 0) {
        print_reference_quantities(ref_quant);
    }

    input.set_input_from_arguments(global_args);

    model_constants_t cc(input.tracking_filename);

    if (global_args.checkpoint_flag && rank == 0) {
        checkpoint.load_checkpoint(global_args.checkpoint_string, cc,
            y_init, segments, input, ref_quant);
        print_segments(segments);
    } else {
        cc.setup_model_constants(input, ref_quant);
        if (global_args.ens_config_flag && rank == 0) {
            load_ens_config(global_args.ens_config_string, cc,
                segments, input, ref_quant);
            for (auto &s : segments)
                SUCCESS_OR_DIE(s.check());
            print_segments(segments);
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


void setup_simulation_base(
    const int &argc,
    char* const * argv,
    const int &rank,
    const int &n_processes,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    std::vector<segment_t> &segments,
    model_constants_t &cc,
    std::vector<double> &y_init,
    std::vector<codi::RealReverse> &y_single_old,
    bool &already_loaded,
    netcdf_reader_t &netcdf_reader) {
    if (!already_loaded) {
        netcdf_reader.init_netcdf(
#ifdef MET3D
            input.start_time,
#endif
            input.INPUT_FILENAME.c_str(),
            (global_args.checkpoint_flag || already_loaded),
            cc, input.simulation_mode, input.current_time);
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

    // // Set "old" values as temporary holder of values.
    // for (int ii = 0 ; ii < num_comp ; ii++)
    //     y_single_old[ii] = y_init[ii];

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
    }
    already_loaded = true;
}

void setup_simulation(
    const int &argc,
    char* const * argv,
    const int &rank,
    const int &n_processes,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    std::vector<segment_t> &segments,
    model_constants_t &cc,
    std::vector<double> &y_init,
    std::vector<codi::RealReverse> &y_single_old,
    checkpoint_t &checkpoint,
    output_handle_t &out_handler,
    bool &already_loaded,
    netcdf_reader_t &netcdf_reader) {
    checkpoint.load_checkpoint(cc, y_init, segments, input, ref_quant, out_handler);
    setup_simulation_base(argc, argv, rank, n_processes, input,
            global_args, ref_quant, segments, cc, y_init, y_single_old,
            already_loaded, netcdf_reader);
}

void substep_trace(
    const uint32_t &sub,
    const uint32_t &t,
    const model_constants_t &cc,
    const input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const std::vector<codi::RealReverse> &y_single_old,
    const std::vector<codi::RealReverse> &inflow) {
#if defined(TRACE_SAT) || defined(TRACE_QR) || defined(TRACE_QV) || defined(TRACE_QC) || defined(TRACE_QI) \
    || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
#if defined(TRACE_TIME)
#if defined(MET3D)
    trace = (((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) + input.start_time >= trace_start)
        && ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) + input.start_time <= trace_end)) ? true : false;
#else
    trace = (((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime)  >= trace_start)
        && ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) <= trace_end)) ? true : false;
#endif
#endif
    if (trace)
        std::cout << cc.id << " timestep : " << (sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) << "\n";

#endif
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

void parameter_check(
    std::vector<segment_t> &segments,
    model_constants_t &cc,
    const double &time_old,
    const std::vector< std::array<double, num_par > > &y_diff,
    const std::vector<codi::RealReverse> &y_single_old,
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    task_scheduler_t &scheduler) {
    for (auto &s : segments) {
        bool perturb = s.perturb_check(
            cc, y_diff, y_single_old, time_old);
        if (perturb) {
            // Perturb this instance
            if (s.n_members == 1) {
                std::string descr;
                s.perturb(cc, ref_quant, input, descr);
            } else {
                // Create a checkpoint for new members
#ifdef USE_MPI
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
                            &ens_id,  // origin
                            &scheduler.max_ensemble_id,  // compare
                            &result,  // result
                            MPI_UINT64_T,  // datatype
                            0,  // target rank
                            0,  // target displ
                            scheduler.ens_window));
                    } while (result != scheduler.max_ensemble_id);
                } else {
                    SUCCESS_OR_DIE(
                        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, scheduler.my_rank, 0, scheduler.ens_window));
                    ens_id = scheduler.max_ensemble_id + 1;
                    scheduler.max_ensemble_id = ens_id;
                    MPI_Win_unlock(scheduler.my_rank, scheduler.ens_window);
                }
                if (ens_id < cc.n_ensembles) {
                    for (uint32_t i=0; i < s.n_members; ++i) {
                        const uint64_t total_members = s.n_members;
                        checkpoint_t checkpoint(
                            cc,
                            y_single_old,
                            segments,
                            input,
                            time_old,
                            i,
                            ens_id,
                            total_members);
                        scheduler.send_new_task(checkpoint);
                    }
                }
#else
                checkpoint_t checkpoint(
                    cc,
                    y_single_old,
                    segments,
                    input,
                    time_old);
                std::string checkpoint_filename = input.FOLDER_NAME;
                checkpoint.write_checkpoint(checkpoint_filename, cc, segments);
                create_run_script(
                    input.FOLDER_NAME,
                    checkpoint_filename,
                    cc,
                    s.n_members,
                    "gnuparallel",
                    false);
#endif
            }
        }
    }
}

void finish_last_step(
    std::vector<codi::RealReverse> &y_single_new,
    const reference_quantities_t &ref_quant,
    const model_constants_t &cc) {
    codi::RealReverse T_prime = y_single_new[T_idx]*ref_quant.Tref;
    codi::RealReverse p_prime = y_single_new[p_idx]*ref_quant.pref;
    codi::RealReverse qv_prime = y_single_new[qv_idx]*ref_quant.qref;
    codi::RealReverse qc_prime = y_single_new[qc_idx]*ref_quant.qref;
    codi::RealReverse p_sat = saturation_pressure_water(
        T_prime,
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b));
#ifdef TRACE_ENV
    if (trace)
        std::cout << "before sat ad S " << y_single_new[S_idx]
            << "\nbefore sat ad S calc "
            << convert_qv_to_S(
                p_prime,
                T_prime,
                qv_prime,
                get_at(cc.constants, Cons_idx::p_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_a),
                get_at(cc.constants, Cons_idx::T_sat_low_temp),
                get_at(cc.constants, Cons_idx::p_sat_const_b),
                get_at(cc.constants, Cons_idx::Epsilon)) << "\n";
#endif
    std::vector<codi::RealReverse> res(7);
    for (auto& r : res) r = 0;
    saturation_adjust(
        T_prime,
        p_prime,
        p_sat,
        qv_prime,
        qc_prime,
        res,
        cc);
    y_single_new[qv_idx] += res[qv_idx]/ref_quant.qref;
    y_single_new[qc_idx] += res[qc_idx]/ref_quant.qref;
    y_single_new[T_idx] += res[T_idx]/ref_quant.Tref;
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
#ifdef TRACE_ENV
    if (trace)
        std::cout << "sat ad S " << y_single_new[S_idx]
            << "\nsat ad T " << y_single_new[T_idx]*ref_quant.Tref << "\n";
#endif
}

void run_substeps(
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const uint32_t &t,
    model_constants_t &cc,
    std::vector<codi::RealReverse> &y_single_old,
    std::vector<codi::RealReverse> &inflow,
    codi::RealReverse::TapeType &tape,
    std::vector<codi::RealReverse> &y_single_new,
    netcdf_reader_t &netcdf_reader,
    std::vector< std::array<double, num_par > > &y_diff,
    output_handle_t &out_handler,
    const uint32_t &sub_start,
    const uint32_t &ensemble,
    std::vector<segment_t> &segments,
    ProgressBar &pbar,
    task_scheduler_t &scheduler,
    const uint64_t &progress_index) {

    double time_old, time_new;
    for (uint32_t sub=sub_start; sub <= cc.num_sub_steps; ++sub) {
        substep_trace(sub, t, cc, input, ref_quant, y_single_old, inflow);

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
    #endif
#endif
#if defined MET3D && defined TURBULENCE
        y_single_old[qv_idx] += inflow[qv_in_idx]/cc.num_sub_steps;
#endif
        cc.register_input(tape);
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
        }

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
        cc.get_gradients(y_single_new, y_diff, tape);

        // Time update
        time_new = (sub + t*cc.num_sub_steps)*cc.dt;
        out_handler.process_step(cc, netcdf_reader, y_single_new, y_diff, sub, t,
            time_new, input.write_index,
            input.snapshot_index,
#ifdef MET3D
            ensemble,
#endif
            last_step, ref_quant);
        // Interchange old and new for next step
        time_old = time_new;
        y_single_old.swap(y_single_new);

        if (time_new != cc.t_end_prime) {
            // Check if parameter shall be perturbed
            parameter_check(segments, cc, time_old, y_diff,
                y_single_old, input, ref_quant, scheduler);
        }

        // While debugging, the bar is not useful.
#if !defined(TRACE_SAT) && !defined(TRACE_ENV) && !defined(TRACE_QV) && !defined(TRACE_QC) && !defined(TRACE_QR)
#if !defined(TRACE_QS) && !defined(TRACE_QI) && !defined(TRACE_QG) && !defined(TRACE_QH)
        if (progress_index > 0)
            pbar.progress(sub + t*cc.num_sub_steps);
#endif
#endif
        // In case that the sub timestep%timestep is not 0
        if (last_step)
            break;
    }  // End substep
}

int run_simulation(
    const int &rank,
    const int &n_processes,
    model_constants_t &cc,
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const global_args_t &global_args,
    std::vector<codi::RealReverse> &y_single_old,
    std::vector< std::array<double, num_par > > &y_diff,
    std::vector<codi::RealReverse> &y_single_new,
    std::vector<codi::RealReverse> &inflow,
    output_handle_t &out_handler,
    std::vector<segment_t> &segments,
    task_scheduler_t &scheduler,
    netcdf_reader_t &netcdf_reader) {

#ifdef MET3D
    uint32_t ensemble;
#endif
    // force any process that is not root to disable pbar
    const uint64_t progress_index = (rank != 0) ? 0 : input.progress_index;
    ProgressBar pbar = ProgressBar((cc.num_sub_steps)*cc.num_steps,
        progress_index, "simulation step", std::cout);
    // Loop for timestepping: BEGIN
    try {
        codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
        uint32_t sub_start = 1;
        if (global_args.checkpoint_flag && std::fmod(input.current_time, cc.dt_prime) != 0)
            sub_start = std::fmod(input.current_time, cc.dt_prime)
                    / (cc.dt_prime/(cc.num_sub_steps));

        // Loop over every timestep that is usually fixed to 20 s
        for (uint32_t t=0; t < cc.num_steps - cc.done_steps; ++t) {
            netcdf_reader.read_buffer(cc, ref_quant, y_single_old,
                inflow, t, global_args.checkpoint_flag, input.start_over_env);
            // Iterate over each substep
            run_substeps(input, ref_quant, t, cc, y_single_old,
                inflow, tape, y_single_new, netcdf_reader, y_diff, out_handler,
                sub_start, ensemble, segments, pbar, scheduler,
                progress_index);
#ifdef TRACE_QG
            if (trace)
                std::cout << "\nSediment total q: " << sediment_q_total
                        << "\nSediment total N: " << sediment_n_total << "\n";
            sediment_n_total = 0;
            sediment_q_total = 0;
#endif
            sub_start = 1;
            checkpoint_t throw_away;
            scheduler.send_task(throw_away, false);
        }
        if (progress_index > 0)
            pbar.finish();
    } catch(netCDF::exceptions::NcException& e) {
        if (progress_index > 0)
            pbar.finish();
        std::cout << e.what() << std::endl;
        std::cout << "rank " << rank << ": ABORTING. (sorry)" << std::endl;
        return 1;
    }
    return 0;
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
    std::vector<segment_t> segments;
    reference_quantities_t ref_quant;

    // Collect the initial values in a separated vector
    // See constants.h for storage order
    std::vector<double> y_init(num_comp);
    checkpoint_t checkpoint;
    bool already_loaded = false;
    // Hold the derivatives of all components
    std::vector< std::array<double, num_par > >  y_diff(num_comp);

    // Allocate vectors for the single solution with unperturbed
    // parameter and assign initial values
    std::vector<codi::RealReverse> y_single_old(num_comp);
    std::vector<codi::RealReverse> y_single_new(num_comp);
    std::vector<codi::RealReverse> inflow(num_inflows);

    model_constants_t cc = parse_args(argc, argv, rank, n_processes, input,
        global_args, ref_quant, segments, y_init, checkpoint);

    netcdf_reader_t netcdf_reader(input.write_index);
    if ((input.simulation_mode == grid_sensitivity) || (input.simulation_mode == trajectory_sensitivity)) {
        netcdf_reader.set_dims(input.INPUT_FILENAME.c_str(), cc, input.simulation_mode);
    }
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, cc,
            ref_quant, input.INPUT_FILENAME, input.write_index,
            input.snapshot_index, rank, input.simulation_mode);
    task_scheduler_t scheduler(rank, n_processes, input.simulation_mode);

    if ((input.simulation_mode == grid_sensitivity) || (input.simulation_mode == trajectory_sensitivity)) {
        // static scheduling with parallel read and write enabled
        setup_simulation_base(argc, argv, rank, n_processes, input,
            global_args, ref_quant, segments, cc, y_init, y_single_old,
            already_loaded, netcdf_reader);

        scheduler.set_n_ensembles(netcdf_reader.n_ensembles);
        scheduler.set_n_trajectories(netcdf_reader.n_trajectories);

        if (scheduler.receive_task(checkpoint)) {
            cc.traj_id = scheduler.current_traj;
            cc.ensemble_id = scheduler.current_ens;

            netcdf_reader.read_initial_values(y_init, ref_quant, cc,
                global_args.checkpoint_flag, cc.traj_id, cc.ensemble_id);
            // Set "old" values as temporary holder of values.
            for (int ii = 0 ; ii < num_comp ; ii++)
                y_single_old[ii] = y_init[ii];

            out_handler.reset(scheduler.current_traj, scheduler.current_ens);
            // run simulation
            SUCCESS_OR_DIE(run_simulation(rank, n_processes, cc, input, ref_quant,
                global_args, y_single_old, y_diff, y_single_new, inflow,
                out_handler, segments, scheduler, netcdf_reader));
            while (scheduler.receive_task(checkpoint)) {
                setup_simulation_base(argc, argv, rank, n_processes, input,
                    global_args, ref_quant, segments, cc, y_init, y_single_old,
                    already_loaded, netcdf_reader);

                netcdf_reader.read_initial_values(y_init, ref_quant, cc,
                    global_args.checkpoint_flag, scheduler.current_traj, scheduler.current_ens);
                out_handler.reset(scheduler.current_traj, scheduler.current_ens);
                cc.traj_id = scheduler.current_traj;
                cc.ensemble_id = scheduler.current_ens;
                // Set "old" values as temporary holder of values.
                for (int ii = 0 ; ii < num_comp ; ii++)
                    y_single_old[ii] = y_init[ii];

                // run simulation
                SUCCESS_OR_DIE(run_simulation(rank, n_processes, cc, input, ref_quant,
                    global_args, y_single_old, y_diff, y_single_new, inflow,
                    out_handler, segments, scheduler, netcdf_reader));
            }
        }
    } else {   // dynamic scheduling with parallel read and write disabled
        netcdf_reader.set_dims(input.INPUT_FILENAME.c_str(), cc, input.simulation_mode);
        if (rank == 0) {
            setup_simulation(argc, argv, rank, n_processes, input,
                global_args, ref_quant, segments, cc, y_init, y_single_old,
                checkpoint, out_handler, already_loaded, netcdf_reader);

            netcdf_reader.read_initial_values(y_init, ref_quant, cc,
                global_args.checkpoint_flag, input.traj, input.ensemble);

            // Set "old" values as temporary holder of values.
            for (int ii = 0 ; ii < num_comp ; ii++)
                y_single_old[ii] = y_init[ii];
            SUCCESS_OR_DIE(run_simulation(rank, n_processes, cc, input, ref_quant,
                global_args, y_single_old, y_diff, y_single_new, inflow,
                out_handler, segments, scheduler, netcdf_reader));
        }

        while (scheduler.receive_task(checkpoint)) {
            setup_simulation(argc, argv, rank, n_processes, input,
                global_args, ref_quant, segments, cc, y_init, y_single_old,
                checkpoint, out_handler, already_loaded, netcdf_reader);

            netcdf_reader.read_initial_values(y_init, ref_quant, cc,
                global_args.checkpoint_flag, input.traj, input.ensemble);
            // Set "old" values as temporary holder of values.
            for (int ii = 0 ; ii < num_comp ; ii++)
                y_single_old[ii] = y_init[ii];

            // run simulation
            SUCCESS_OR_DIE(run_simulation(rank, n_processes, cc, input, ref_quant,
                global_args, y_single_old, y_diff, y_single_new, inflow,
                out_handler, segments, scheduler, netcdf_reader));
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
