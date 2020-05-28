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

#include <netcdf>

#ifdef RK4
#include "include/microphysics/rk4.h"
#endif

#ifdef RK4NOICE
#include "include/microphysics/rk4.h"
#endif

#ifdef RK4ICE
#include "include/microphysics/rk4.h"
#endif

#ifdef EXPLICIT_EULER
#include "include/microphysics/euler.h"
#endif

#ifdef IMPLICIT_EULER
#include "include/microphysics/implicit_euler.h"
#endif

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

    bool need_to_abort = false;
    int opt;
    global_args_t global_args;
    init_global_args(global_args);

    if(argc < 2)
    {
        need_to_abort = true;
        display_usage();
    }else
    {
        opt = getopt(argc, argv, optString);

        while(-1 != opt)
        {
            switch(opt)
            {
                case 'f':
                {
                    global_args.final_time_flag = 1;
                    global_args.final_time_string = optarg;
                    break;
                }
                case 'd':
                {
                    global_args.timestep_flag = 1;
                    global_args.timestep_string = optarg;
                    break;
                }
                case 'i':
                {
                    global_args.snapshot_index_flag = 1;
                    global_args.snapshot_index_string = optarg;
                    break;
                }
                case 'b':
                {
                    global_args.scaling_fact_flag = 1;
                    global_args.scaling_fact_string = optarg;
                    break;
                }
                case 'o':
                {
                    global_args.output_flag = 1;
                    global_args.output_string = optarg;
                    break;
                }
                case 'l':
                {
                    global_args.input_flag = 1;
                    global_args.input_file = optarg;
                    break;
                }
                case 's':
                {
                    global_args.start_over_flag = 1;
                    global_args.start_over_string = optarg;
                    break;
                }
                case 't':
                {
                    global_args.fixed_iteration_flag = 1;
                    global_args.fixed_iteration_string = optarg;
                    break;
                }
                case 'a':
                {
                    global_args.auto_type_flag = 1;
                    global_args.auto_type_string = optarg;
                    break;
                }
                case 'r':
                {
                    global_args.traj_flag = 1;
                    global_args.traj_string = optarg;
                    break;
                }
                case 'w':
                {
                    global_args.write_flag = 1;
                    global_args.write_string = optarg;
                    break;
                }
                case '?':
                {
                    need_to_abort = true;
                    display_usage();
                    break;
                }
                default:
                {
                    need_to_abort = true;
                    display_error_on_command_line();
                    display_usage();
                    break;
                }
            }

            opt = getopt(argc, argv, optString);
        }
    }

    if(need_to_abort){
        std::cout << "ABORTING." << std::endl;
        return 1;			// Report error
    }

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
        global_args.input_file, input.traj, cc) != 0)
    {
        return 1;
    }

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
#if defined(TRACE_QR) || defined(TRACE_QV) || defined(TRACE_QC) || defined(TRACE_QI) || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
                std::cout << "timestep : " << (sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) << "\n";
#endif
                // Set the coefficients from the last timestep and from
                // the input files
                // *Should* only be necessary when parameters from the
                // trajectory are used as start point
                setCoefficients(y_single_old, cc, ref_quant);
                // CODIPACK: BEGIN

                tape.setActive();

                //	 Add the inflow
                y_single_old[qi_idx] += inflow[qi_in_idx]/cc.num_sub_steps;
                y_single_old[qs_idx] += inflow[qs_in_idx]/cc.num_sub_steps;
                y_single_old[qr_idx] += inflow[qr_in_idx]/cc.num_sub_steps;
                y_single_old[qg_idx] += inflow[qg_in_idx]/cc.num_sub_steps;
                y_single_old[Ni_idx] += inflow[Ni_in_idx]/cc.num_sub_steps;
                y_single_old[Ns_idx] += inflow[Ns_in_idx]/cc.num_sub_steps;
                y_single_old[Nr_idx] += inflow[Nr_in_idx]/cc.num_sub_steps;
                y_single_old[Ng_idx] += inflow[Ng_in_idx]/cc.num_sub_steps;
#ifdef TRACE_QR
                std::cout << "Adding qr " << inflow[qr_in_idx]/cc.num_sub_steps
                    << ", Nr " << inflow[Nr_in_idx]/cc.num_sub_steps << "\n";
#endif

                // Dimensional coefficients
                tape.registerInput(cc.a1_prime);    // Autoconversion
                tape.registerInput(cc.a2_prime);    // Accretion
                tape.registerInput(cc.e1_prime);    // Evaporation
                tape.registerInput(cc.e2_prime);    // Evaporation
                tape.registerInput(cc.d_prime);     // Sedimentation
                tape.registerInput(cc.Nc_prime);    // Concentration of cloud droplets

                // Exponents
                tape.registerInput(cc.gamma);       // Autoconversion
                tape.registerInput(cc.betac);       // Accretion
                tape.registerInput(cc.betar);       // Accretion
                tape.registerInput(cc.delta1);      // Evaporation
                tape.registerInput(cc.delta2);      // Evaporation
                tape.registerInput(cc.zeta);        // Sedimentation

                // ICON parameters
                tape.registerInput(rain_gfak);
                tape.registerInput(cloud_k_au);
                tape.registerInput(cloud_k_sc);
                tape.registerInput(kc_autocon);
                tape.registerInput(inv_z);
                cc.rain.register_input(tape);
                cc.cloud.register_input(tape);
                cc.hail.register_input(tape);
                cc.ice.register_input(tape);
                cc.snow.register_input(tape);
                cc.graupel.register_input(tape);
                // CODIPACK: END
//////////////// Add any different scheme and model here
    // I did not check if those two methods still work with CODIPACK
    // #if defined(EXPLICIT_EULER)
    //             euler_step(y_single_new, y_single_old, num_comp, ref_quant, cc);
    // #endif

    // #if defined(IMPLICIT_EULER)
    //             implicit_euler_step(y_single_new, y_single_old, num_comp, ref_quant, cc);
    // #endif

#if defined(RK4)
                RK4_step(y_single_new, y_single_old, ref_quant, cc,
                    nc_params, input.fixed_iteration);
#endif
#if defined(RK4NOICE)
                RK4_step_2_no_ice(y_single_new, y_single_old, ref_quant, cc,
                    nc_params, input.fixed_iteration);
#endif
#if defined(RK4ICE)
                RK4_step_2_sb_ice(y_single_new, y_single_old, ref_quant, cc,
                    nc_params, input.fixed_iteration);
#endif
                // CODIPACK: BEGIN
                for(int ii = 0 ; ii < num_comp ; ii++)
                    tape.registerOutput(y_single_new[ii]);

                tape.setPassive();
                for(int ii = 0 ; ii < num_comp ; ii++)
                {
                    y_single_new[ii].setGradient(1.0);
                    tape.evaluate();

                    y_diff[ii][0] = cc.a1_prime.getGradient();
                    y_diff[ii][1] = cc.a2_prime.getGradient();
                    y_diff[ii][2] = cc.e1_prime.getGradient();
                    y_diff[ii][3] = cc.e2_prime.getGradient();
                    y_diff[ii][4] = cc.d_prime.getGradient();

                    y_diff[ii][5] = cc.gamma.getGradient();
                    y_diff[ii][6] = cc.betac.getGradient();
                    y_diff[ii][7] = cc.betar.getGradient();
                    y_diff[ii][8] = cc.delta1.getGradient();
                    y_diff[ii][9] = cc.delta2.getGradient();
                    y_diff[ii][10] = cc.zeta.getGradient();

                    y_diff[ii][11] = cc.Nc_prime.getGradient();

                    y_diff[ii][12] = rain_gfak.getGradient();
                    y_diff[ii][13] = cloud_k_au.getGradient();
                    y_diff[ii][14] = cloud_k_sc.getGradient();
                    y_diff[ii][15] = kc_autocon.getGradient();
                    y_diff[ii][16] = inv_z.getGradient();

                    uint64_t idx = 17;
                    cc.rain.get_gradient(y_diff[ii], idx);
                    cc.cloud.get_gradient(y_diff[ii], idx);
                    cc.graupel.get_gradient(y_diff[ii], idx);
                    cc.hail.get_gradient(y_diff[ii], idx);
                    cc.ice.get_gradient(y_diff[ii], idx);
                    cc.snow.get_gradient(y_diff[ii], idx);
                    tape.clearAdjoints();
                }
                tape.reset();
                // CODIPACK: END

                // Time update
                time_new = (sub + t*cc.num_sub_steps)*cc.dt;

                // Output *if needed* (checks within function)
                write_output(cc, nc_params, y_single_new, y_diff, sub, t,
                    time_new, traj_id, input.write_index,
                    input.snapshot_index);

                // Interchange old and new for next step
                time_old = time_new;
                    y_single_old.swap(y_single_new);
            } // End substep
        }
        nc_close(ncid);
    } catch(netCDF::exceptions::NcException& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "ABORTING." << std::endl;
        return 1;
    }
    outfile.close();
    for(int ii = 0 ; ii < num_comp ; ii++)
        out_diff[ii].close();
    std::cout << "-------------------FINISHED-------------------\n\n\n";
}
