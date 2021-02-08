#include "codi.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <fstream>
#include <cmath>
#include <tgmath.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>

#include "include/microphysics/constants.h"
#include "include/types/global_args_t.h"
#include "include/types/input_parameters_t.h"
#include "include/types/output_handle_t.h"
#include "include/types/model_constants_t.h"
#include "include/types/nc_parameters_t.h"
#include "include/types/reference_quantities_t.h"
#include "include/types/segment_t.h"

#include "include/microphysics/physical_parameterizations.h"
#include "include/microphysics/program_io.h"

#include "include/microphysics/general.h"

#include <netcdf>

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

void parse_args_and_checkpoints(
    const int &argc,
    char* const * argv,
    const int &rank,
    const int &n_processes,
    input_parameters_t &input,
    global_args_t &global_args,
    reference_quantities_t &ref_quant,
    std::vector<segment_t> &segments,
    model_constants_t &cc,
    std::vector<double> &y_init)
{
    SUCCESS_OR_DIE(parse_arguments(argc, argv, global_args, rank, n_processes));

    ref_quant.Tref = 273.15;
#ifdef WCB
    ref_quant.qref = 1.0e-6;
#else
    ref_quant.qref = 1.0e-4;
#endif
    ref_quant.pref = 1.0e5;
    ref_quant.wref = 1.; // 10.0
    ref_quant.tref = 1.0;
    ref_quant.zref = 1.0;
    ref_quant.Nref = 1.0; 	// DUMMY
    if(rank == 0)
    {
        print_reference_quantities(ref_quant);
    }
    set_input_from_arguments(global_args, input);

    if(rank == 0)
    {
        if(global_args.checkpoint_flag)
        {
            load_checkpoint(global_args.checkpoint_string, cc, y_init, segments, input, ref_quant);
            print_segments(segments);
        }
        else
        {
            cc.setup_model_constants(input, ref_quant);
            if(global_args.ens_config_flag)
            {
                load_ens_config(global_args.ens_config_string, cc, segments, input, ref_quant);
                for(auto &s: segments)
                    SUCCESS_OR_DIE(s.check());
                print_segments(segments);
            }
            input.set_outputfile_id(cc.id, cc.ensemble_id);
        }
    }

    auto_type = input.auto_type;
    load_lookup_table(cc.ltabdminwgg);
}

void substep_trace(
    const uint32_t &sub,
    const model_constants_t &cc,
    const input_parameters_t &input,
    const reference_quantities_t &ref_quant,
    const std::vector<codi::RealReverse> &y_single_old,
    const std::vector<codi::RealReverse> &inflow)
{
#if defined(TRACE_SAT) || defined(TRACE_QR) || defined(TRACE_QV) || defined(TRACE_QC) || defined(TRACE_QI) || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
#if defined(TRACE_TIME)
#if defined(MET3D)
    trace = ( ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) + input.start_time >= trace_start) && ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) + input.start_time <= trace_end) ) ? true : false;
#else
    trace = ( ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime)  >= trace_start) && ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) <= trace_end) ) ? true : false;
#endif
#endif
    if(trace)
        std::cout << cc.id << " timestep : " << (sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) << "\n";

#endif
#if defined(TRACE_SAT)
    if(trace)
        std::cout << "Start qv_prime " << y_single_old[qv_idx]*ref_quant.qref << "\n";
#endif
#ifdef TRACE_QR
    if(trace)
    {
        std::cout << "Adding qr " << inflow[qr_in_idx]/cc.num_sub_steps*ref_quant.qref
            << ", Nr " << inflow[Nr_in_idx]/cc.num_sub_steps << "\n";
        std::cout << "qr_old " << y_single_old[qr_idx]*ref_quant.qref
            << ", Nr_old " << y_single_old[Nr_idx] << "\n";
    }
#endif
#ifdef TRACE_QI
    if(trace)
    {
        std::cout << "Adding qi " << inflow[qi_in_idx]/cc.num_sub_steps*ref_quant.qref
            << ", Ni " << inflow[Ni_in_idx]/cc.num_sub_steps << "\n";
        std::cout << "qi_old " << y_single_old[qi_idx]*ref_quant.qref
            << ", Ni_old " << y_single_old[Ni_idx] << "\n";
    }
#endif
#ifdef TRACE_QS
    if(trace)
    {
        std::cout << "Adding qs " << inflow[qs_in_idx]/cc.num_sub_steps*ref_quant.qref
            << ", Ns " << inflow[Ns_in_idx]/cc.num_sub_steps << "\n";
        std::cout << "qs_old " << y_single_old[qs_idx]*ref_quant.qref
            << ", Ns_old " << y_single_old[Ns_idx] << "\n";
    }
#endif
#ifdef TRACE_QG
    if(trace)
    {
        std::cout << "Adding qg " << inflow[qg_in_idx]/cc.num_sub_steps*ref_quant.qref
            << ", Ng " << inflow[Ng_in_idx]/cc.num_sub_steps << "\n";
        std::cout << "qg_old " << y_single_old[qg_idx]*ref_quant.qref
            << ", Ng_old " << y_single_old[Ng_idx] << "\n";
    }
#endif

#if defined TRACE_QV && defined MET3D && defined TURBULENCE
    if(trace)
    {
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
    const reference_quantities_t &ref_quant)
{
    for(auto &s: segments)
    {
        bool perturb = s.perturb_check(
            cc, y_diff, y_single_old, time_old);
        if(perturb)
        {
            // Perturb this instance
            if(s.n_members == 1)
            {
                s.perturb(cc, ref_quant, input);
            } else // Create a checkpoint for new members
            {
                std::string checkpoint_filename = input.FOLDER_NAME;
                write_checkpoint(
                    checkpoint_filename,
                    cc,
                    y_single_old,
                    segments,
                    input,
                    time_old);
                create_run_script(
                    input.FOLDER_NAME,
                    checkpoint_filename,
                    cc,
                    s.n_members,
                    "gnuparallel",
                    false);
            }
        }
    }
}

void finish_last_step(
    std::vector<codi::RealReverse> &y_single_new,
    const reference_quantities_t &ref_quant,
    const model_constants_t &cc)
{
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
    if(trace)
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
    for(auto& r: res) r = 0;
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
    if(trace)
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
    nc_parameters_t &nc_params,
    std::vector< std::array<double, num_par > > &y_diff,
    output_handle_t &out_handler,
    const uint32_t &sub_start,
    const int &traj_id,
    const uint32_t &ensemble,
    std::vector<segment_t> &segments,
    ProgressBar &pbar)
{
    double time_old, time_new;
    for(uint32_t sub=sub_start; sub<=cc.num_sub_steps-input.start_over; ++sub) // cc.num_sub_steps
    {
        substep_trace(sub, cc, input, ref_quant, y_single_old, inflow);

        bool last_step = ( ((sub+1 + t*cc.num_sub_steps) >= ((t+1)*cc.num_sub_steps + !input.start_over))
            || (sub == cc.num_sub_steps-input.start_over) );
#if defined(RK4_ONE_MOMENT)
        // Set the coefficients from the last timestep and from
        // the input files
        // *Should* only be necessary when parameters from the
        // trajectory are used as start point
        cc.setCoefficients(y_single_old, ref_quant);
#endif
        // CODIPACK
        tape.setActive();

        //	 Add the inflow
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
#if defined MET3D && defined TURBULENCE
        y_single_old[qv_idx] += inflow[qv_in_idx]/cc.num_sub_steps;
#endif
        cc.register_input(tape);
        if(sub == 1)
        {
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
        RK4_step(y_single_new, y_single_old, ref_quant, cc,
            nc_params, input.fixed_iteration);
#elif defined(RK4NOICE)
        RK4_step_2_no_ice(y_single_new, y_single_old, ref_quant, cc,
            nc_params, input.fixed_iteration);
#elif defined(RK4ICE)
        RK4_step_2_sb_ice(y_single_new, y_single_old, ref_quant, cc,
            nc_params, input.fixed_iteration);
#endif
#ifndef IN_SAT_ADJ
        if(last_step)
        {
            finish_last_step(y_single_new, ref_quant, cc);
        }
#endif
        cc.get_gradients(y_single_new, y_diff, tape);

        // Time update
        time_new = (sub + t*cc.num_sub_steps)*cc.dt;

        // Output *if needed* (checks within function)
        write_output(cc, nc_params, y_single_new, y_diff, sub, t,
            time_new, traj_id, input.write_index,
            input.snapshot_index,
#ifdef MET3D
            ensemble,
#endif
            last_step, out_handler, ref_quant);

        // Interchange old and new for next step
        time_old = time_new;
        y_single_old.swap(y_single_new);

        // Check if parameter shall be perturbed
        parameter_check(segments, cc, time_old, y_diff,
            y_single_old, input, ref_quant);

        // While debugging, the bar is not useful.
#if !defined(TRACE_SAT) && !defined(TRACE_ENV) && !defined(TRACE_QV) && !defined(TRACE_QC) && !defined(TRACE_QR) && !defined(TRACE_QS) && !defined(TRACE_QI) && !defined(TRACE_QG) && !defined(TRACE_QH)
        if(input.progress_index > 0)
            pbar.progress(sub + t*cc.num_sub_steps);
#endif
        // In case that the sub timestep%timestep is not 0
        if( last_step )
            break;
    } // End substep
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
    nc_parameters_t &nc_params,
    output_handle_t &out_handler,
    std::vector<segment_t> &segments,
    const size_t &lenp)
{
    int traj_id = -1;
#ifdef MET3D
    uint32_t ensemble;
#endif

    ProgressBar pbar = ProgressBar((cc.num_sub_steps-input.start_over)*cc.num_steps, input.progress_index, "simulation step", std::cout);
    // Loop for timestepping: BEGIN
    try
    {
        // Needed to read from trajectory file
        std::vector<int> ids(lenp);
        std::vector<size_t> startp, countp;
        resize_counter(startp, countp, input.traj);
        codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
        uint32_t sub_start = 1;
        if(global_args.checkpoint_flag && std::fmod(input.current_time, cc.dt_prime) != 0)
            sub_start = std::fmod(input.current_time, cc.dt_prime)
                    / (cc.dt_prime/(cc.num_sub_steps-input.start_over));
        // Loop over every timestep that is usually fixed to 20 s
        for(uint32_t t=0; t<cc.num_steps; ++t) //
        {
            read_netcdf_write_stream(input.INPUT_FILENAME.c_str(), startp, countp,
                nc_params, cc, input, ref_quant, y_single_old, inflow, ids,
                traj_id,
#ifdef MET3D
                ensemble,
#endif
                t, global_args.checkpoint_flag, out_handler);

            // Iterate over each substep
            run_substeps(input, ref_quant, t, cc, y_single_old,
                inflow, tape, y_single_new, nc_params, y_diff, out_handler,
                sub_start, traj_id, ensemble, segments, pbar);

#ifdef TRACE_QG
            if(trace)
                std::cout << "\nSediment total q: " << sediment_q_total
                        << "\nSediment total N: " << sediment_n_total << "\n";
            sediment_n_total = 0;
            sediment_q_total = 0;
#endif
            sub_start = 1;
        }
        if(input.progress_index > 0)
            pbar.finish();
    } catch(netCDF::exceptions::NcException& e)
    {
        if(input.progress_index > 0)
            pbar.finish();
        std::cout << e.what() << std::endl;
        std::cout << "ABORTING." << std::endl;
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
int main(int argc, char** argv)
{
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
    model_constants_t cc;
    // Collect the initial values in a separated vector
    // See constants.h for storage order
    std::vector<double> y_init(num_comp);

    parse_args_and_checkpoints(argc, argv, rank, n_processes, input, global_args,
        ref_quant, segments, cc, y_init);

    nc_parameters_t nc_params;
    size_t lenp;
    SUCCESS_OR_DIE(read_init_netcdf(y_init, nc_params, lenp, ref_quant,
#ifdef MET3D
        input.start_time,
#endif
        input.INPUT_FILENAME.c_str(), input.traj, global_args.checkpoint_flag,
        cc, input.current_time));

#if defined(RK4_ONE_MOMENT)
    cc.setCoefficients(y_init[0] , y_init[1]);
#endif
    if(rank == 0)
    {
        print_constants(cc);
        print_input_parameters(input);
    }
    // Hold the derivatives of all components
    std::vector< std::array<double, num_par > >  y_diff(num_comp);

    // Allocate vectors for the single solution with unperturbed
    // parameter and assign initial values
    std::vector<codi::RealReverse> y_single_old(num_comp);
    std::vector<codi::RealReverse> y_single_new(num_comp);
    std::vector<codi::RealReverse> inflow(num_inflows);

    // Set "old" values as temporary holder of values.
    for(int ii = 0 ; ii < num_comp ; ii++)
        y_single_old[ii] = y_init[ii];

    // Currently we only load one trajectory per instance
    output_handle_t out_handler("netcdf", input.OUTPUT_FILENAME, 1, 1, cc,
        ref_quant, input.INPUT_FILENAME, input.write_index, input.snapshot_index);
    if(rank == 0)
    {
#ifdef TRACE_QC
        print_particle_params(cc.cloud, "cloud");
#endif
#ifdef TRACE_QR
        print_particle_params(cc.rain, "rain");
#endif

#ifdef TRACE_QI
        print_particle_params(cc.ice, "ice");
#endif

#ifdef TRACE_QS
        print_particle_params(cc.snow, "snow");
#endif

#ifdef TRACE_QG
        print_particle_params(cc.graupel, "graupel");
#endif

#ifdef TRACE_QH
        print_particle_params(cc.hail, "hail");
#endif
    }


    SUCCESS_OR_DIE(run_simulation(rank, n_processes, cc, input, ref_quant,
        global_args, y_single_old, y_diff, y_single_new, inflow, nc_params,
        out_handler, segments, lenp));

    // Better:
    // something like
    // if(scheduler.check_queue()){
    //     run_simulation();
    // }
    // else end this
#ifdef USE_MPI
    MPI_Finalize();
#endif
#ifdef SILENT_MODE
    exit(0);
#else
    std::cout << "-------------------FINISHED-------------------\n\n\n";
    exit(0);
#endif
}
