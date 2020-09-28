#include "codi.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tgmath.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>
#include "include/microphysics/physical_parameterizations.h"

#include "include/microphysics/program_io.h"
#include "include/microphysics/constants.h"
#include "include/microphysics/general.h"
#include "include/microphysics/gradient_handle.h"

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
    input_parameters_t input;
    init_input_parameters(input);
    global_args_t global_args;
    int err;
    if( (err = parse_arguments(argc, argv, global_args)) != 0 ) return err;

    set_input_from_arguments(global_args, input);
    auto_type = input.auto_type;
    load_lookup_table(ltabdminwgg);

    int traj_id;

    // ==================================================
    // Define the reference quantities
    // ==================================================
    reference_quantities_t ref_quant;

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

    // Print the reference quantities
    print_reference_quantities(ref_quant);

    // ==================================================
    // Setup the model constants
    // ==================================================
    model_constants_t cc;
    setup_model_constants(input, cc, ref_quant);

    // ==================================================
    // Allocate the memory
    // ==================================================

    // Current time used in the iterations
    double time_old = 0.0;
    double time_new = 0.0;

    // Collect the initial values in a separated vector
    // Storage:
    // See constants.h for storage order
    std::vector<double> y_init(num_comp);
    nc_parameters_t nc_params;
    size_t lenp;
    if(read_init_netcdf(y_init, nc_params, lenp, ref_quant,
#ifdef MET3D
        input.start_time,
#endif
        global_args.input_file, input.traj, cc) != 0)
    {
        return 1;
    }
#ifdef MET3D
    if(write_attributes(input.INPUT_FILENAME, input.OUTPUT_FILENAME) != 0)
    {
        return 1;
    }
#endif

    setCoefficients(y_init[0] , y_init[1], cc);

    // Print the constants
    print_constants(cc);

    // Print the input parameters
    print_input_parameters(input);

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

    // ==================================================
    // Arrange the Output
    // ==================================================
    if(write_reference_quantities(input.OUTPUT_FILENAME, ref_quant) != 0)
    {
        return 1;
    }
    print_reference_quantities(ref_quant);


    if(write_headers(input.OUTPUT_FILENAME) != 0)
    {
        return 1;
    }

    ProgressBar pbar = ProgressBar((cc.num_sub_steps-2)*cc.num_steps, input.progress_index, "simulation step", std::cout);
    // Loop for timestepping: BEGIN
    try
    {
        // Needed to read from trajectory file
        std::vector<int> ids(lenp);
        int ncid;
        std::vector<size_t> startp, countp;
        open_netcdf(ncid, startp, countp, global_args.input_file, input.traj);
        codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

        // Loop over every timestep that is usually fixed to 20 s
        for(uint32_t t=0; t<cc.num_steps; ++t) //
        {
            read_netcdf_write_stream(global_args.input_file, startp, countp,
                nc_params, cc, input, ref_quant, y_single_old, inflow, ids,
                traj_id, t);

            // Iterate over each substep
            for(uint32_t sub=1; sub<cc.num_sub_steps; ++sub) // cc.num_sub_steps
            {
#if defined(TRACE_SAT) || defined(TRACE_QR) || defined(TRACE_QV) || defined(TRACE_QC) || defined(TRACE_QI) || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
#if defined(TRACE_TIME)
                trace = ( ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) >= trace_start) && ((sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) <= trace_end) ) ? true : false;
#endif
                if(trace)
                    std::cout << "timestep : " << (sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) << "\n";

#endif
                bool last_step = ((sub+1 + t*cc.num_sub_steps)*cc.dt >= ((t+1)*cc.num_sub_steps)*cc.dt) || (sub == cc.num_sub_steps-input.start_over);
                // Set the coefficients from the last timestep and from
                // the input files
                // *Should* only be necessary when parameters from the
                // trajectory are used as start point
                setCoefficients(y_single_old, cc, ref_quant);
#if defined(TRACE_SAT)
            std::cout << "Start qv_prime " << y_single_old[qv_idx]*ref_quant.qref << "\n";
#endif
                // CODIPACK
                tape.setActive();
#ifdef TRACE_QR
                if(trace)
                {
                    std::cout << "Adding qr " << inflow[qr_in_idx]/cc.num_sub_steps
                        << ", Nr " << inflow[Nr_in_idx]/cc.num_sub_steps << "\n";
                    std::cout << "qr_old " << y_single_old[qr_idx]
                        << ", Nr_old " << y_single_old[Nr_idx] << "\n";
                }
#endif
#ifdef TRACE_QI
                if(trace)
                {
                    std::cout << "Adding qi " << inflow[qi_in_idx]/cc.num_sub_steps
                        << ", Ni " << inflow[Ni_in_idx]/cc.num_sub_steps << "\n";
                    std::cout << "qi_old " << y_single_old[qi_idx]
                        << ", Ni_old " << y_single_old[Ni_idx] << "\n";
                }
#endif
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
                register_everything(tape, cc);
                if(sub == 1)
                {
                                    //	 Add the inflow
//                     y_single_old[qr_idx] += inflow[qr_in_idx];
//                     y_single_old[Nr_idx] += inflow[Nr_in_idx];
// #if defined(RK4ICE)
//                     y_single_old[qi_idx] += inflow[qi_in_idx];
//                     y_single_old[qs_idx] += inflow[qs_in_idx];
//                     y_single_old[qg_idx] += inflow[qg_in_idx];
//                     y_single_old[Ni_idx] += inflow[Ni_in_idx];
//                     y_single_old[Ns_idx] += inflow[Ns_in_idx];
//                     y_single_old[Ng_idx] += inflow[Ng_in_idx];
// #endif
                    // codi::RealReverse T_prime = y_single_old[T_idx]*ref_quant.Tref;
                    // codi::RealReverse p_prime = y_single_old[p_idx]*ref_quant.pref;
                    // codi::RealReverse qv_prime = y_single_old[qv_idx]*ref_quant.qref;
                    // codi::RealReverse qc_prime = y_single_old[qc_idx]*ref_quant.qref;
                    // codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime);
                    // std::vector<codi::RealReverse> res(7);
                    // for(auto& r: res) r = 0;
                    // saturation_adjust_meteo(
                    //     T_prime,
                    //     p_prime,
                    //     p_sat,
                    //     qv_prime,
                    //     qc_prime,
                    //     res,
                    //     ref_quant.qref);
                    // y_single_old[qv_idx] += res[qv_idx]/ref_quant.qref;
                    // y_single_old[qc_idx] += res[qc_idx]/ref_quant.qref;
                    // y_single_old[T_idx] += res[T_idx]/ref_quant.Tref;
                    // y_single_old[S_idx] = convert_qv_to_S(
                    //     y_single_old[p_idx].getValue()*ref_quant.pref,
                    //     y_single_old[T_idx].getValue()*ref_quant.Tref,
                    //     y_single_old[qv_idx].getValue()*ref_quant.qref);
                }
//////////////// Add any different scheme and model here
    // I did not check if those two methods still work with CODIPACK
    // #if defined(EXPLICIT_EULER)
    //             euler_step(y_single_new, y_single_old, num_comp, ref_quant, cc);
    // #endif

    // #if defined(IMPLICIT_EULER)
    //             implicit_euler_step(y_single_new, y_single_old, num_comp, ref_quant, cc);
    // #endif


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
                if(last_step)
                {
                    // codi::RealReverse T_prime = y_single_new[T_idx]*ref_quant.Tref;
                    // codi::RealReverse p_prime = y_single_new[p_idx]*ref_quant.pref;
                    // codi::RealReverse qv_prime = y_single_new[qv_idx]*ref_quant.qref;
                    // codi::RealReverse qc_prime = y_single_new[qc_idx]*ref_quant.qref;
                    // codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime);
                    // std::vector<codi::RealReverse> res(7);
                    // for(auto& r: res) r = 0;
                    // saturation_adjust_meteo(
                    //     T_prime,
                    //     p_prime,
                    //     p_sat,
                    //     qv_prime,
                    //     qc_prime,
                    //     res,
                    //     ref_quant.qref);
                    // y_single_new[qv_idx] += res[qv_idx]/ref_quant.qref;
                    // y_single_new[qc_idx] += res[qc_idx]/ref_quant.qref;
                    // y_single_new[T_idx] += res[T_idx]/ref_quant.Tref;
                    // y_single_new[S_idx] = convert_qv_to_S(
                    //     y_single_new[p_idx].getValue()*ref_quant.pref,
                    //     y_single_new[T_idx].getValue()*ref_quant.Tref,
                    //     y_single_new[qv_idx].getValue()*ref_quant.qref);
// #ifdef TRACE_QV
//                     std::cout << "sat ad dqv " << res[qv_idx]
//                         << ", dqc " << res[qc_idx] << "\n";
// #endif
                }
                // CODIPACK: BEGIN
                get_gradients(y_single_new, y_diff, cc, tape);
                // CODIPACK: END

                // Time update
                time_new = (sub + t*cc.num_sub_steps)*cc.dt;

                // Output *if needed* (checks within function)
                write_output(cc, nc_params, y_single_new, y_diff, sub, t,
                    time_new, traj_id, input.write_index,
                    input.snapshot_index, last_step);

                // Interchange old and new for next step
                time_old = time_new;
                y_single_old.swap(y_single_new);
                // While debugging, the bar is not useful.
#if !defined(TRACE_SAT) && !defined(TRACE_ENV) && !defined(TRACE_QV) && !defined(TRACE_QC) && !defined(TRACE_QR) && !defined(TRACE_QS) && !defined(TRACE_QI) && !defined(TRACE_QG) && !defined(TRACE_QH)
                pbar.progress(sub + t*cc.num_sub_steps);
#endif
                // In case that the sub timestep%timestep is not 0
                if( last_step )
                    break;
            } // End substep
#ifdef TRACE_QG
        if(trace)
            std::cout << "\nSediment total q: " << sediment_q_total
                    << "\nSediment total N: " << sediment_n_total << "\n";
        sediment_n_total = 0;
        sediment_q_total = 0;
#endif
        }
        pbar.finish();
        nc_close(ncid);
    } catch(netCDF::exceptions::NcException& e)
    {
        pbar.finish();
        std::cout << e.what() << std::endl;
        std::cout << "ABORTING." << std::endl;
        return 1;
    }
    outfile.close();
    for(int ii = 0 ; ii < num_comp ; ii++)
        out_diff[ii].close();
    std::cout << "-------------------FINISHED-------------------\n\n\n";
}
