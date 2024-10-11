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

#include <codi.hpp>

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

    input.set_input_from_arguments(global_args);

    model_constants_t cc(input.tracking_filename);

    cc.setup_model_constants(input, ref_quant);
    if (global_args.ens_config_flag && rank == 0) {
        load_ens_config(global_args.ens_config_string, cc,
            segments, input, ref_quant);
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
    model_constants_t &cc,
    std::vector<double> &y_init,
    std::vector<codi::RealReverseIndex> &y_single_old,
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
        if (rank != 0)
            already_loaded = true;
    } else {
        netcdf_reader.time_idx = netcdf_reader.start_time_idx;
    }
#if defined(RK4_ONE_MOMENT)
    cc.setCoefficients(y_init[0] , y_init[1]);
#endif

    if (rank == 0 && !already_loaded) {
        already_loaded = true;
    }
}

void finish_last_step(
    std::vector<codi::RealReverseIndex> &y_single_new,
    const reference_quantities_t &ref_quant,
    const model_constants_t &cc) {
    codi::RealReverseIndex T_prime = y_single_new[T_idx]*ref_quant.Tref;
    codi::RealReverseIndex p_prime = y_single_new[p_idx]*ref_quant.pref;
    codi::RealReverseIndex qv_prime = y_single_new[qv_idx]*ref_quant.qref;
    codi::RealReverseIndex qc_prime = y_single_new[qc_idx]*ref_quant.qref;
    codi::RealReverseIndex p_sat = saturation_pressure_water(
        T_prime,
        get_at(cc.constants, Cons_idx::p_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_a),
        get_at(cc.constants, Cons_idx::T_sat_low_temp),
        get_at(cc.constants, Cons_idx::p_sat_const_b));

    std::vector<codi::RealReverseIndex> res(7);
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
}

void run_substeps(
    input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const uint32_t &t,
    model_constants_t &cc,
    std::vector<codi::RealReverseIndex> &y_single_old,
    std::vector<codi::RealReverseIndex> &inflow,
    codi::RealReverseIndex::Tape &tape,
    std::vector<codi::RealReverseIndex> &y_single_new,
    netcdf_reader_t &netcdf_reader,
    std::vector< std::array<double, num_par > > &y_diff,
    const uint32_t &sub_start,
    const uint32_t &ensemble,
    std::vector<segment_t> &segments,
    task_scheduler_t &scheduler,
    const uint64_t &progress_index,
    const double delay_out_time = 0) {

    for (uint32_t sub=sub_start; sub <= cc.num_sub_steps; ++sub) {
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
            codi::RealReverseIndex p_prime = y_single_old[p_idx]*ref_quant.pref;
            codi::RealReverseIndex T_prime = y_single_old[T_idx]*ref_quant.Tref;
            codi::RealReverseIndex qv_prime = y_single_old[qv_idx]*ref_quant.qref;
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

        // Interchange old and new for next step
        y_single_old.swap(y_single_new);
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
    std::vector<codi::RealReverseIndex> &y_single_old,
    std::vector< std::array<double, num_par > > &y_diff,
    std::vector<codi::RealReverseIndex> &y_single_new,
    std::vector<codi::RealReverseIndex> &inflow,
    std::vector<segment_t> &segments,
    task_scheduler_t &scheduler,
    netcdf_reader_t &netcdf_reader,
    const double delay_out_time = 0) {

#ifdef MET3D
    uint32_t ensemble;
#endif
    // force any process that is not root to disable pbar
    const uint64_t progress_index = (rank != 0) ? 0 : input.progress_index;
    // Loop for timestepping: BEGIN
    try {
        codi::RealReverseIndex::Tape& tape = codi::RealReverseIndex::getTape();
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
                inflow, tape, y_single_new, netcdf_reader, y_diff,
                sub_start, ensemble, segments, scheduler,
                progress_index, delay_out_time);
            sub_start = 1;
        }
    } catch(netCDF::exceptions::NcException& e) {
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
    if (rank == 0)
        std::cout << "#timesteps,#input_params,#output_params,time (s)\n";

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
    std::vector<codi::RealReverseIndex> y_single_old(num_comp);
    std::vector<codi::RealReverseIndex> y_single_new(num_comp);
    std::vector<codi::RealReverseIndex> inflow(num_inflows);

    model_constants_t cc = parse_args(argc, argv, rank, n_processes, input,
        global_args, ref_quant, segments, y_init, checkpoint);

    netcdf_reader_t netcdf_reader(input.write_index);
    if ((input.simulation_mode == grid_sensitivity)
        || (input.simulation_mode == trajectory_sensitivity)
        || (input.simulation_mode == limited_time_ensembles)) {
        netcdf_reader.set_dims(input.INPUT_FILENAME.c_str(), cc, input.simulation_mode);
    }
    task_scheduler_t scheduler(rank, n_processes, input.simulation_mode);

    if (rank == 0) {
        for (int iters=0; iters < 20; iters++) {
            for (int n_out=num_comp; n_out >= 0; n_out--) {
                for (int n_in=num_par; n_in >= 0; n_in--) {
                    // state_param from 0 to Cons_idx::n_items
                    cc.track_state = 0;
                    for (auto &tp : cc.track_param) tp = 0;

                    for (int no=0; no < n_out; no++) cc.track_state += pow(2, no);
                    for (int ni=0; ni < n_in; ni++) {
                        uint32_t idx = ni/64;
                        cc.track_param[idx] += pow(2, ni - idx*64);
                    }

                    setup_simulation_base(argc, argv, rank, n_processes, input,
                        global_args, ref_quant, cc, y_init, y_single_old,
                        already_loaded, netcdf_reader);

                    scheduler.set_n_ensembles(netcdf_reader.n_ensembles);
                    scheduler.set_n_trajectories(netcdf_reader.n_trajectories);
                    cc.traj_id = scheduler.current_traj;
                    cc.ensemble_id = scheduler.current_ens;
                    cc.traj_id = 0;
                    cc.ensemble_id = 0;
                    netcdf_reader.read_initial_values(y_init, ref_quant, cc,
                        global_args.checkpoint_flag, cc.traj_id, cc.ensemble_id);
                    // Set "old" values as temporary holder of values.
                    for (int ii = 0 ; ii < num_comp ; ii++)
                        y_single_old[ii] = y_init[ii];

                    auto t_first = std::chrono::system_clock::now();
                    // run simulation
                    SUCCESS_OR_DIE(run_simulation(rank, n_processes, cc, input, ref_quant,
                        global_args, y_single_old, y_diff, y_single_new, inflow,
                        segments, scheduler, netcdf_reader));

                    auto now = std::chrono::system_clock::now();
                    double time = ((std::chrono::duration<double>)(now - t_first)).count();
                    std::cout << cc.num_steps - cc.done_steps << ","
                        << n_in << ","
                        << n_out << ","
                        << time << "\n";
                }
            }
        }
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif
#ifdef SILENT_MODE
    exit(0);
#else
    exit(0);
#endif
}
