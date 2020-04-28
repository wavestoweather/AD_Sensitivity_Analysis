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
    const uint64_t n_lookup = 2000;
    const uint64_t n_lookup_highres = 10000;
    const uint64_t n_lookup_hr_dummy = 10;

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

    // Scaling factor from input
    cc.scaling_fact = input.scaling_fact;

    // Accomodation coefficient
    cc.alpha_d = 1.0;


    // // Performance constants for warm cloud; COSMO
    cc.a1_scale = 1.0e-3;
    cc.a2_scale = 1.72 / pow(Ra , 7./8.);
    cc.e1_scale = 1.0 / sqrt(Ra);
    cc.e2_scale = 9.1 / pow(Ra , 11./16.);
    cc.d_scale = ( 130.0*tgamma(4.5) )/( 6.0*(1.0e3)*pow(M_PI*(8.0e6)*Ra , 1.0/8.0) );

    // Performance constants for warm cloud; IFS
    // The file constants.h also defines some constants as nar, ...
    const double Nc = 50; 	// 50 over ocean; 300 over land
    const double F_aut = 1.5;
    const double F_acc = 2.0;
    const double lambdar_pp = pow(cc.nar * cc.ar * tgamma(cc.br + 1.0) , cc.alphar);

    // cc.a1_scale = (1350. * F_aut)/pow(Nc , 1.79);
    // cc.a2_scale = 67.0 * F_acc;
    // cc.e1_scale = 2.0 * M_PI * cc.nar * ( (0.78 * tgamma(2.0 - cc.nbr))/(lambdar_pp*lambdar_pp) );
    // cc.e2_scale = cc.scaling_fact * 2.0 * M_PI * cc.nar * 0.31
    //     * pow(cc.cr/cc.mu, 0.5) * pow(cc.Sc, 1.0/3.0) * pow(cc.rho0, 0.25)
    //     * (tgamma(cc.epsilonr + cc.nbr)/pow(lambdar_pp ,cc.epsilonr));
    // cc.d_scale = 4.0e-3;

    // Inflow from above
    cc.B_prime = 0.0; //1.0e-7;

    // // Exponents of the cloud model
    // // COSMO
    cc.gamma = 1.0;
    cc.betac = 1.0;
    cc.betar = 7./8.;
    cc.delta1 = 0.5;
    cc.delta2 = 11./16.;
    cc.zeta = 9./8.;

    // Exponents of the cloud model
    // IFS
    // cc.gamma = 2.47;
    // cc.betac = 1.15;
    // cc.betar = 1.15;
    // cc.delta1 = 2.0/( cc.br + 1.0 - cc.nbr );
    // cc.delta2 = ( 0.5*cc.dr + 2.5 - cc.nbr )/( cc.br + 1.0 - cc.nbr );
    // cc.zeta = 1.0;

    // Numerics
    cc.t_end_prime = input.t_end_prime;
    cc.t_end = input.t_end_prime/ref_quant.tref;
    // Time of the substeps

    cc.dt = input.dt_prime/ref_quant.tref;
    cc.dt_prime = input.dt_prime;
    cc.dt_traject_prime = cc.dt_traject * ref_quant.tref;
    // The trajectories are calculated with 20 s timesteps.
    cc.num_sub_steps = (floor( 20.0/cc.dt ) < 1) ? 1 : floor( 20.0/cc.dt );

    cc.snapshot_index = input.snapshot_index;
    // Evaluate the general performance constants
    cc.dt_half = cc.dt*0.5;
    cc.dt_third = cc.dt/3.0;
    cc.dt_sixth = cc.dt/6.0;

    // ==================================================
    // Set rain constants
    // See init_2mom_scheme_once in mo_2mom_mcrph_main.f90
    // ==================================================
    // Cosmo5 although nue1nue1 would be okay too, I guess
    //// Cloud
    cc.cloud.nu = 0.0;
    cc.cloud.mu = 1.0/3.0;
    cc.cloud.max_x = 2.6e-10;
    cc.cloud.min_x = 4.2e-15;
    cc.cloud.min_x_act = 4.2e-15;
    cc.cloud.min_x_nuc_homo = 4.2e-15;
    cc.cloud.min_x_nuc_hetero = 4.2e-15;
    cc.cloud.min_x_melt = 4.2e-15;
    cc.cloud.min_x_evap = 4.2e-15;
    cc.cloud.min_x_freezing = 4.2e-15;
    cc.cloud.min_x_depo = 4.2e-15;
    cc.cloud.min_x_collision = 4.2e-15;
    cc.cloud.min_x_collection = 4.2e-15;
    cc.cloud.min_x_conversion = 4.2e-15;
    cc.cloud.min_x_sedimentation = 4.2e-15;
    cc.cloud.a_geo = 1.24e-1;
    cc.cloud.b_geo = 0.333333;
    cc.cloud.a_vel = 3.75e5;
    cc.cloud.b_vel = 0.666667;
    cc.cloud.a_ven = 0.78;
    cc.cloud.b_ven = 0.308;
    cc.cloud.cap = 2.0;
    cc.cloud.vsedi_max = 1.0;
    cc.cloud.vsedi_min = 0.0;
    cc.cloud.c_s = 1.0 / cc.cloud.cap;
    cc.cloud.a_f = vent_coeff_a(cc.cloud, 1);
    cc.cloud.b_f = vent_coeff_b(cc.cloud, 1) * pow(N_sc, n_f) / sqrt(kin_visc_air);
    cc.cloud.c_z = moment_gamma(cc.cloud, 2);

    setup_cloud_autoconversion(cc.cloud);
    setup_bulk_sedi(cc.cloud);

    //// Rain
    cc.rain.nu = 0.0;
    cc.rain.mu = 0.333333;
    cc.rain.max_x = 3.0e-6;
    cc.rain.min_x = 2.6e-10;
    cc.rain.min_x_act = 2.6e-10;
    cc.rain.min_x_nuc_homo = 2.6e-10;
    cc.rain.min_x_nuc_hetero = 2.6e-10;
    cc.rain.min_x_melt = 2.6e-10;
    cc.rain.min_x_evap = 2.6e-10;
    cc.rain.min_x_freezing = 2.6e-10;
    cc.rain.min_x_depo = 2.6e-10;
    cc.rain.min_x_collision = 2.6e-10;
    cc.rain.min_x_collection = 2.6e-10;
    cc.rain.min_x_conversion = 2.6e-10;
    cc.rain.min_x_sedimentation = 2.6e-10;
    cc.rain.a_geo = 1.24e-1;
    cc.rain.b_geo = 0.333333;
    cc.rain.a_vel = 114.0137;
    cc.rain.b_vel = 0.234370;
    cc.rain.cap = 2.0;

    // From rainSBBcoeffs
    cc.rain.alpha = 9.292;
    cc.rain.beta = 9.623;
    cc.rain.gamma = 6.222e2;
    cc.rain.cmu0 = 6.0;
    cc.rain.cmu1 = 3.0e1;
    cc.rain.cmu2 = 1.0e3;
    cc.rain.cmu3 = 1.1e-3;
    cc.rain.cmu4 = 1.0;
    cc.rain.cmu5 = 2.0;
    rain_gfak = 1.0;

    cc.rain.nm1 = (cc.rain.nu+1.0)/cc.rain.mu;
    cc.rain.nm2 = (cc.rain.nu+2.0)/cc.rain.mu;
    cc.rain.nm3 = (cc.rain.nu+3.0)/cc.rain.mu;
    cc.rain.vsedi_max = 20.0;
    cc.rain.vsedi_min = 0.1;
    init_gamma_table(table_r1, n_lookup, n_lookup_hr_dummy, cc.rain.nm1.getValue());
    init_gamma_table(table_r2, n_lookup, n_lookup_hr_dummy, cc.rain.nm2.getValue());
    init_gamma_table(table_r3, n_lookup, n_lookup_hr_dummy, cc.rain.nm3.getValue());
    cc.rain.g1 = table_r1.igf[table_r1.n_bins-1];
    cc.rain.g2 = table_r2.igf[table_r2.n_bins-1];
    cc.rain.c_s = 1.0 / cc.rain.cap;
    cc.rain.a_f = vent_coeff_a(cc.rain, 1);
    cc.rain.b_f = vent_coeff_b(cc.rain, 1) * pow(N_sc, n_f) / sqrt(kin_visc_air);
    cc.rain.c_z = moment_gamma(cc.rain, 2);
    setup_bulk_sedi(cc.rain);

    //// Graupel
    cc.graupel.nu = 1.0;
    cc.graupel.mu = 1.0/3.0; // Not used actually?
    cc.graupel.max_x = 5.0e-4;
    cc.graupel.min_x = 1.0e-9;
    cc.graupel.min_x_act = 1.0e-9;
    cc.graupel.min_x_nuc_homo = 1.0e-9;
    cc.graupel.min_x_nuc_hetero = 1.0e-9;
    cc.graupel.min_x_melt = 1.0e-9;
    cc.graupel.min_x_evap = 1.0e-9;
    cc.graupel.min_x_freezing = 1.0e-9;
    cc.graupel.min_x_depo = 1.0e-9;
    cc.graupel.min_x_collision = 1.0e-9;
    cc.graupel.min_x_collection = 1.0e-9;
    cc.graupel.min_x_conversion = 1.0e-9;
    cc.graupel.min_x_sedimentation = 1.0e-9;
    cc.graupel.a_geo = 1.42e-1;
    cc.graupel.b_geo = 0.314;
    cc.graupel.a_vel = 86.89371;
    cc.graupel.b_vel = 0.268325;
    cc.graupel.a_ven = 0.78;
    cc.graupel.b_ven = 0.308;
    cc.graupel.cap = 2.0;
    cc.graupel.vsedi_max = 30.0;
    cc.graupel.vsedi_min = 0.1;
    cc.graupel.sc_coll_n = 1.0;
    cc.graupel.d_crit_c = 100.0e-6;
    cc.graupel.q_crit_c = 1.0e-6;
    cc.graupel.s_vel = 0.0;

    cc.graupel.nm1 = (cc.graupel.nu+1.0)/cc.graupel.mu;
    cc.graupel.nm2 = (cc.graupel.nu+2.0)/cc.graupel.mu;
    codi::RealReverse a = (cc.graupel.nu+1.0)/cc.graupel.mu;
    init_gamma_table(table_g1, n_lookup, n_lookup_hr_dummy, cc.graupel.nm1.getValue());
    a = (cc.graupel.nu+2.0)/cc.graupel.mu;
    init_gamma_table(table_g2, n_lookup, n_lookup_hr_dummy, cc.graupel.nm2.getValue());
    cc.graupel.g1 = table_g1.igf[table_g1.n_bins-1];
    cc.graupel.g2 = table_g2.igf[table_g2.n_bins-1];
    cc.graupel.c_s = 1.0 / cc.graupel.cap;
    cc.graupel.a_f = vent_coeff_a(cc.graupel, 1);
    cc.graupel.b_f = vent_coeff_b(cc.graupel, 1) * pow(N_sc, n_f) / sqrt(kin_visc_air);
    cc.graupel.c_z = moment_gamma(cc.graupel, 2);
    setup_bulk_sedi(cc.graupel);

    //// Hail
    cc.hail.nu = 1.0;
    cc.hail.mu = 1.0/3.0; // Not used actually?
    cc.hail.max_x = 5.0e-4;
    cc.hail.min_x = 2.6e-9;
    cc.hail.min_x_act= 2.6e-9;
    cc.hail.min_x_nuc_homo = 2.6e-9;
    cc.hail.min_x_nuc_hetero = 2.6e-9;
    cc.hail.min_x_melt = 2.6e-9;
    cc.hail.min_x_evap = 2.6e-9;
    cc.hail.min_x_freezing = 2.6e-9;
    cc.hail.min_x_depo = 2.6e-9;
    cc.hail.min_x_collision = 2.6e-9;
    cc.hail.min_x_collection = 2.6e-9;
    cc.hail.min_x_conversion = 2.6e-9;
    cc.hail.min_x_sedimentation = 2.6e-9;
    cc.hail.a_geo = 0.1366;
    cc.hail.b_geo = 1.0/3.0;
    cc.hail.a_vel = 39.3;
    cc.hail.b_vel = 0.166667;
    cc.hail.a_ven = 0.78;
    cc.hail.b_ven = 0.308;
    cc.hail.cap = 2.0;
    cc.hail.vsedi_max = 30.0;
    cc.hail.vsedi_min = 0.1;
    cc.hail.sc_coll_n = 1.0;
    cc.hail.d_crit_c = 100.0e-6;
    cc.hail.q_crit_c = 1.0e-6;
    cc.hail.s_vel = 0.0;
    cc.hail.c_s = 1.0 / cc.hail.cap;
    cc.hail.a_f = vent_coeff_a(cc.hail, 1);
    cc.hail.b_f = vent_coeff_b(cc.hail, 1) * pow(N_sc, n_f) / sqrt(kin_visc_air);
    cc.hail.c_z = moment_gamma(cc.hail, 2);
    setup_bulk_sedi(cc.hail);

    //// Ice
    cc.ice.nu = 0.0;
    cc.ice.mu = 1.0/3.0; // Not used actually?
    cc.ice.max_x = 1.0e-5;
    cc.ice.min_x = 1.0e-12;
    cc.ice.min_x_act = 1.0e-12;
    cc.ice.min_x_nuc_homo = 1.0e-12;
    cc.ice.min_x_nuc_hetero = 1.0e-12;
    cc.ice.min_x_melt = 1.0e-12;
    cc.ice.min_x_evap = 1.0e-12;
    cc.ice.min_x_freezing = 1.0e-12;
    cc.ice.min_x_depo = 1.0e-12;
    cc.ice.min_x_collision = 1.0e-12;
    cc.ice.min_x_collection = 1.0e-12;
    cc.ice.min_x_conversion = 1.0e-12;
    cc.ice.min_x_sedimentation = 1.0e-12;
    cc.ice.a_geo = 0.835;
    cc.ice.b_geo = 0.39;
    cc.ice.a_vel = 2.77e1;
    cc.ice.b_vel = 0.21579;
    cc.ice.a_ven = 0.78;
    cc.ice.b_ven = 0.308;
    cc.ice.cap = 2.0;
    cc.ice.vsedi_max = 3.0;
    cc.ice.vsedi_min = 0.0;
    cc.ice.sc_coll_n = 0.8;
    cc.ice.d_crit_c = 150.0e-6;
    cc.ice.q_crit_c = 1.0e-5;
    cc.ice.s_vel = 0.05;
    cc.ice.c_s = 1.0 / cc.ice.cap;
    cc.ice.a_f = vent_coeff_a(cc.ice, 1);
    cc.ice.b_f = vent_coeff_b(cc.ice, 1) * pow(N_sc, n_f) / sqrt(kin_visc_air);
    cc.ice.c_z = moment_gamma(cc.ice, 2);
    setup_bulk_sedi(cc.ice);

    //// Snow
    cc.snow.nu = 0.0;
    cc.snow.mu = 0.5; // Not used actually?
    cc.snow.max_x = 2.0e-5;
    cc.snow.min_x = 1.0e-10;
    cc.snow.min_x_act = 1.0e-10;
    cc.snow.min_x_nuc_homo = 1.0e-10;
    cc.snow.min_x_nuc_hetero = 1.0e-10;
    cc.snow.min_x_melt = 1.0e-10;
    cc.snow.min_x_evap = 1.0e-10;
    cc.snow.min_x_freezing = 1.0e-10;
    cc.snow.min_x_depo = 1.0e-10;
    cc.snow.min_x_collision = 1.0e-10;
    cc.snow.min_x_collection = 1.0e-10;
    cc.snow.min_x_conversion = 1.0e-10;
    cc.snow.min_x_sedimentation = 1.0e-10;
    cc.snow.a_geo = 2.4;
    cc.snow.b_geo = 0.455;
    cc.snow.a_vel = 8.8;
    cc.snow.b_vel = 0.15;
    cc.snow.a_ven = 0.78;
    cc.snow.b_ven = 0.308;
    cc.snow.cap = 2.0;
    cc.snow.vsedi_max = 3.0;
    cc.snow.vsedi_min = 0.1;
    cc.snow.sc_coll_n = 0.8;
    cc.snow.d_crit_c = 150.0e-6;
    cc.snow.q_crit_c = 1.0e-5;
    cc.snow.s_vel = 0.25;
    cc.snow.c_s = 1.0 / cc.snow.cap;
    cc.snow.a_f = vent_coeff_a(cc.snow, 1);
    cc.snow.b_f = vent_coeff_b(cc.snow, 1) * pow(N_sc, n_f) / sqrt(kin_visc_air);
    cc.snow.c_z = moment_gamma(cc.snow, 2);
    setup_bulk_sedi(cc.snow);

    init_particle_collection_1(cc.snow, cc.cloud, cc.coeffs_scr);
    init_particle_collection_2(cc.snow, cc.rain, cc.coeffs_srr);
    init_particle_collection_2(cc.ice, cc.rain, cc.coeffs_irr);
    init_particle_collection_1(cc.ice, cc.cloud, cc.coeffs_icr);
    init_particle_collection_1(cc.hail, cc.rain, cc.coeffs_hrr);
    init_particle_collection_1(cc.graupel, cc.rain, cc.coeffs_grr);
    init_particle_collection_1(cc.hail, cc.cloud, cc.coeffs_hcr);
    init_particle_collection_1(cc.graupel, cc.cloud, cc.coeffs_gcr);
    init_particle_collection_1(cc.snow, cc.ice, cc.coeffs_sic);
    init_particle_collection_1(cc.hail, cc.ice, cc.coeffs_hic);
    init_particle_collection_1(cc.graupel, cc.ice, cc.coeffs_gic);
    init_particle_collection_1(cc.hail, cc.snow, cc.coeffs_hsc);
    init_particle_collection_1(cc.graupel, cc.snow, cc.coeffs_gsc);

    // ==================================================
    // Allocate the memory
    // ==================================================

    // Current time
    double time_old = 0.0;
    double time_new = 0.0;

    // Number of components of the solution
    // const int num_comp = 7;
    // const int num_comp = 24;

    // Collect the initial values in a separated vector
    // Storage:
    // See constants.h for storage order
    std::vector<double> y_init(num_comp);
    nc_parameters_t nc_params;

    try
    {
        int dimid, ncid;
        size_t lenp, n_timesteps;
        // Get the amount of trajectories
        nc_open(global_args.input_file, NC_NOWRITE, &ncid);
#ifdef WCB
        nc_inq_dimid(ncid, "ntra", &dimid);
#else
        nc_inq_dimid(ncid, "id", &dimid);
#endif
        nc_inq_dimlen(ncid, dimid, &lenp);
        std::cout << "Number of trajectories in netCDF file: " << lenp << "\n";
        if(lenp <= input.traj)
        {
            std::cout << "You asked for trajectory with index " << input.traj
                      << " which does not exist. ABORTING.\n";
            return 1;
        }
        // Get the amount of timesteps
#ifdef WCB
        nc_inq_dimid(ncid, "ntim", &dimid);
#else
        nc_inq_dimid(ncid, "time", &dimid);
#endif
        nc_inq_dimlen(ncid, dimid, &n_timesteps);
        uint64_t n_timesteps_input = ceil(cc.t_end/20.0);

        cc.num_steps = (n_timesteps-1 > n_timesteps_input) ? n_timesteps_input : n_timesteps-1;
        init_nc_parameters(nc_params, lenp, n_timesteps);
        netCDF::NcFile datafile(global_args.input_file, netCDF::NcFile::read);
        load_nc_parameters_var(nc_params, datafile);

        std::vector<size_t> startp, countp;
        // wcb files have a different ordering
#if defined WCB || defined WCB2
        startp.push_back(0); // time
        startp.push_back(input.traj); // trajectory id
#else
        startp.push_back(input.traj); // trajectory id
        startp.push_back(1); // time (where time == 0 only has zeros)
#endif
        countp.push_back(1);
        countp.push_back(1);
        load_nc_parameters(nc_params, startp, countp,
                           ref_quant, cc.num_sub_steps);

        y_init[p_idx] = nc_params.p;
        y_init[T_idx] = nc_params.t;

        y_init[S_idx] = nc_params.S;
        y_init[qc_idx] = nc_params.qc;
        y_init[qr_idx] = nc_params.qr;
        y_init[qv_idx] = nc_params.qv;
        y_init[qi_idx] = nc_params.qi;
        y_init[qs_idx] = nc_params.qs;
#ifdef WCB
        y_init[w_idx] = 0;
        y_init[qg_idx] = 0;
#else
        y_init[w_idx] = nc_params.w[0];
        y_init[qg_idx] = nc_params.qg;
#endif

        y_init[qh_idx] = 0.0; // hail that is not in the trajectory

        y_init[qh_out_idx] = 0.0;
        y_init[Nh_out_idx] = 0.0;
#ifdef WCB2
        y_init[qi_out_idx] = nc_params.QIout;
        y_init[qs_out_idx] = nc_params.QSout;
        y_init[qr_out_idx] = nc_params.QRout;
        y_init[qg_out_idx] = nc_params.QGout;

        y_init[Ni_out_idx] = nc_params.NIout;
        y_init[Ns_out_idx] = nc_params.NSout;
        y_init[Nr_out_idx] = nc_params.NRout;
        y_init[Ng_out_idx] = nc_params.NGout;

        y_init[Ni_idx] = nc_params.Ni;
        y_init[Ns_idx] = nc_params.Ns;
        y_init[Nr_idx] = nc_params.Nr;
        y_init[Ng_idx] = nc_params.Ng;
        y_init[Nc_idx] = nc_params.Nc;
#else
        // We initialize the sedimentation with 0 for the stepper
        y_init[qi_out_idx] = 0.0;
        y_init[qs_out_idx] = 0.0;
        y_init[qr_out_idx] = 0.0;
        y_init[qg_out_idx] = 0.0;

        y_init[Ni_out_idx] = 0;
        y_init[Ns_out_idx] = 0;
        y_init[Nr_out_idx] = 0;
        y_init[Ng_out_idx] = 0;

        y_init[Ni_idx] = 0;
        y_init[Ns_idx] = 0;
        y_init[Nr_idx] = 0;
        y_init[Ng_idx] = 0;
        y_init[Nc_idx] = 0;
#endif
        y_init[Nv_idx] = 0;
        y_init[z_idx] = nc_params.z[0];

        y_init[n_inact_idx] = 0;
        y_init[depo_idx] = 0;
        y_init[sub_idx] = 0;

    } catch(netCDF::exceptions::NcException& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "ABORTING." << std::endl;
        return 1;
    }

    // Malloc memory
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

    // Storage
    // y_single_old = [p, T, w, S, qc, qr, qv]
    for(int ii = 0 ; ii < num_comp ; ii++)
        y_single_old[ii] = y_init[ii];

    // ==================================================
    // Arrange the Output
    // ==================================================

    std::string out_filename = input.OUTPUT_FILENAME;
    std::string suffix = ".txt";
    std::string full_filename;
    full_filename = out_filename;
    full_filename += suffix;

    std::ofstream outfile;
    std::stringstream out_tmp;
    outfile.open(full_filename);
    outfile.precision(10);

    if( !outfile.is_open() )
    {
        std::cout << "ERROR while opening the outputfile. Aborting." << std::endl;
        return 1;
    }

    std::ofstream outfile_refs;

    outfile_refs.open(input.OUTPUT_FILENAME + "_reference_values.txt");
    outfile_refs.precision(10);

    if( !outfile_refs.is_open() )
    {
        std::cout << "ERROR while opening the outputfile. Aborting." << std::endl;
        return 1;
    }

    // Write the reference quantities
    // Write the reference quantities
    outfile_refs << ref_quant.Tref << " "
	       << ref_quant.pref << " "
	       << ref_quant.qref << " "
	       << ref_quant.Nref << " "
	       << ref_quant.wref << " "
	       << ref_quant.tref << " "
           << ref_quant.zref << "\n";
    print_reference_quantities(ref_quant);
    outfile_refs.close();

    // Append the initial values and write headers
    out_tmp << "timestep,trajectory,LONGITUDE,LATITUDE,"
#if defined WCB || defined WCB2
        << "MAP,"
#endif
        << "p,T,w,S,qc,qr,qv,Nc,Nr,Nv,qi,Ni,vi,"
        << "qs,Ns,qg,Ng,qh,Nh,qiout,qsout,qrout,qgout,qhout,"
        << "latent_heat,latent_cool,Niout,Nsout,Nrout,Ngout,Nhout,z,Inactive,deposition,sublimination\n";

    // CODIPACK: BEGIN
    std::string basename = "_diff_";
    std::string fname;
    std::ofstream out_diff[num_comp];
    std::stringstream out_diff_tmp[num_comp];


    for(int ii = 0 ; ii < num_comp ; ii++)
    {
        fname = input.OUTPUT_FILENAME;
        fname += basename;
        fname += std::to_string(ii);
        fname += suffix;

        out_diff[ii].open(fname);
        out_diff[ii].precision(10);
        if( !out_diff[ii].is_open() )
        {
            std::cout << "ERROR while opening outputfile. Aborting." << std::endl;
            return 1;
        }
        out_diff_tmp[ii]
            << "timestep,"
            << "trajectory,"
            << "Output Parameter,"
            << "LONGITUDE,"
            << "LATITUDE,"
#if defined WCB || defined WCB2
            << "MAP,"
#endif
            << "da_1,"
            << "da_2,"
            << "de_1,"
            << "de_2,"
            << "dd,"
            << "dN_c,"
            << "dgamma,"
            << "dbeta_c,"
            << "dbeta_r,"
            << "ddelta1,"
            << "ddelta2,"
            << "dzeta,"
            << "drain_gfak,"
            << "dcloud_k_au,"
            << "dcloud_k_sc,"
            << "dkc_autocon,"
            << "dinv_z,"
            // Rain
            << "drain_a_geo,"
            << "drain_b_geo,"
            << "drain_min_x,"
            << "drain_min_x_act,"
            << "drain_min_x_nuc_homo,"
            << "drain_min_x_nuc_hetero,"
            << "drain_min_x_melt,"
            << "drain_min_x_evap,"
            << "drain_min_x_freezing,"
            << "drain_min_x_depo,"
            << "drain_min_x_collision,"
            << "drain_min_x_collection,"
            << "drain_min_x_conversion,"
            << "drain_min_x_sedimentation,"
            << "drain_min_x_riming,"
            << "drain_max_x,"
            << "drain_sc_theta_q,"
            << "drain_sc_delta_q,"
            << "drain_sc_theta_n,"
            << "drain_sc_delta_n,"
            << "drain_s_vel,"
            << "drain_a_vel,"
            << "drain_b_vel,"
            << "drain_rho_v,"
            << "drain_c_z,"
            << "drain_sc_coll_n,"
            << "drain_cmu0,"
            << "drain_cmu1,"
            << "drain_cmu2,"
            << "drain_cmu3,"
            << "drain_cmu4,"
            << "drain_cmu5,"
            << "drain_alpha,"
            << "drain_beta,"
            << "drain_gamma,"
            << "drain_nu,"
            << "drain_g1,"
            << "drain_g2,"
            << "drain_mu,"
            << "drain_nm1,"
            << "drain_nm2,"
            << "drain_nm3,"
            << "drain_q_crit_c,"
            << "drain_d_crit_c,"
            << "drain_ecoll_c,"
            << "drain_cap,"
            << "drain_a_ven,"
            << "drain_b_ven,"
            << "drain_c_s,"
            << "drain_a_f,"
            << "drain_b_f,"
            << "drain_alfa_n,"
            << "drain_alfa_q,"
            << "drain_lambda,"
            << "drain_vsedi_min,"
            << "drain_vsedi_max,"
            // Cloud
            << "dcloud_a_geo,"
            << "dcloud_b_geo,"
            << "dcloud_min_x,"
            << "dcloud_min_x_act,"
            << "dcloud_min_x_nuc_homo,"
            << "dcloud_min_x_nuc_hetero,"
            << "dcloud_min_x_melt,"
            << "dcloud_min_x_evap,"
            << "dcloud_min_x_freezing,"
            << "dcloud_min_x_depo,"
            << "dcloud_min_x_collision,"
            << "dcloud_min_x_collection,"
            << "dcloud_min_x_conversion,"
            << "dcloud_min_x_sedimentation,"
            << "dcloud_min_x_riming,"
            << "dcloud_max_x,"
            << "dcloud_sc_theta_q,"
            << "dcloud_sc_delta_q,"
            << "dcloud_sc_theta_n,"
            << "dcloud_sc_delta_n,"
            << "dcloud_s_vel,"
            << "dcloud_a_vel,"
            << "dcloud_b_vel,"
            << "dcloud_rho_v,"
            << "dcloud_c_z,"
            << "dcloud_sc_coll_n,"
            << "dcloud_cmu0,"
            << "dcloud_cmu1,"
            << "dcloud_cmu2,"
            << "dcloud_cmu3,"
            << "dcloud_cmu4,"
            << "dcloud_cmu5,"
            << "dcloud_alpha,"
            << "dcloud_beta,"
            << "dcloud_gamma,"
            << "dcloud_nu,"
            << "dcloud_g1,"
            << "dcloud_g2,"
            << "dcloud_mu,"
            << "dcloud_nm1,"
            << "dcloud_nm2,"
            << "dcloud_nm3,"
            << "dcloud_q_crit_c,"
            << "dcloud_d_crit_c,"
            << "dcloud_ecoll_c,"
            << "dcloud_cap,"
            << "dcloud_a_ven,"
            << "dcloud_b_ven,"
            << "dcloud_c_s,"
            << "dcloud_a_f,"
            << "dcloud_b_f,"
            << "dcloud_alfa_n,"
            << "dcloud_alfa_q,"
            << "dcloud_lambda,"
            << "dcloud_vsedi_min,"
            << "dcloud_vsedi_max,"
            // Graupel
            << "dgraupel_a_geo,"
            << "dgraupel_b_geo,"
            << "dgraupel_min_x,"
            << "dgraupel_min_x_act,"
            << "dgraupel_min_x_nuc_homo,"
            << "dgraupel_min_x_nuc_hetero,"
            << "dgraupel_min_x_melt,"
            << "dgraupel_min_x_evap,"
            << "dgraupel_min_x_freezing,"
            << "dgraupel_min_x_depo,"
            << "dgraupel_min_x_collision,"
            << "dgraupel_min_x_collection,"
            << "dgraupel_min_x_conversion,"
            << "dgraupel_min_x_sedimentation,"
            << "dgraupel_min_x_riming,"
            << "dgraupel_max_x,"
            << "dgraupel_sc_theta_q,"
            << "dgraupel_sc_delta_q,"
            << "dgraupel_sc_theta_n,"
            << "dgraupel_sc_delta_n,"
            << "dgraupel_s_vel,"
            << "dgraupel_a_vel,"
            << "dgraupel_b_vel,"
            << "dgraupel_rho_v,"
            << "dgraupel_c_z,"
            << "dgraupel_sc_coll_n,"
            << "dgraupel_cmu0,"
            << "dgraupel_cmu1,"
            << "dgraupel_cmu2,"
            << "dgraupel_cmu3,"
            << "dgraupel_cmu4,"
            << "dgraupel_cmu5,"
            << "dgraupel_alpha,"
            << "dgraupel_beta,"
            << "dgraupel_gamma,"
            << "dgraupel_nu,"
            << "dgraupel_g1,"
            << "dgraupel_g2,"
            << "dgraupel_mu,"
            << "dgraupel_nm1,"
            << "dgraupel_nm2,"
            << "dgraupel_nm3,"
            << "dgraupel_q_crit_c,"
            << "dgraupel_d_crit_c,"
            << "dgraupel_ecoll_c,"
            << "dgraupel_cap,"
            << "dgraupel_a_ven,"
            << "dgraupel_b_ven,"
            << "dgraupel_c_s,"
            << "dgraupel_a_f,"
            << "dgraupel_b_f,"
            << "dgraupel_alfa_n,"
            << "dgraupel_alfa_q,"
            << "dgraupel_lambda,"
            << "dgraupel_vsedi_min,"
            << "dgraupel_vsedi_max,"
            // Hail
            << "dhail_a_geo,"
            << "dhail_b_geo,"
            << "dhail_min_x,"
            << "dhail_min_x_act,"
            << "dhail_min_x_nuc_homo,"
            << "dhail_min_x_nuc_hetero,"
            << "dhail_min_x_melt,"
            << "dhail_min_x_evap,"
            << "dhail_min_x_freezing,"
            << "dhail_min_x_depo,"
            << "dhail_min_x_collision,"
            << "dhail_min_x_collection,"
            << "dhail_min_x_conversion,"
            << "dhail_min_x_sedimentation,"
            << "dhail_min_x_riming,"
            << "dhail_max_x,"
            << "dhail_sc_theta_q,"
            << "dhail_sc_delta_q,"
            << "dhail_sc_theta_n,"
            << "dhail_sc_delta_n,"
            << "dhail_s_vel,"
            << "dhail_a_vel,"
            << "dhail_b_vel,"
            << "dhail_rho_v,"
            << "dhail_c_z,"
            << "dhail_sc_coll_n,"
            << "dhail_cmu0,"
            << "dhail_cmu1,"
            << "dhail_cmu2,"
            << "dhail_cmu3,"
            << "dhail_cmu4,"
            << "dhail_cmu5,"
            << "dhail_alpha,"
            << "dhail_beta,"
            << "dhail_gamma,"
            << "dhail_nu,"
            << "dhail_g1,"
            << "dhail_g2,"
            << "dhail_mu,"
            << "dhail_nm1,"
            << "dhail_nm2,"
            << "dhail_nm3,"
            << "dhail_q_crit_c,"
            << "dhail_d_crit_c,"
            << "dhail_ecoll_c,"
            << "dhail_cap,"
            << "dhail_a_ven,"
            << "dhail_b_ven,"
            << "dhail_c_s,"
            << "dhail_a_f,"
            << "dhail_b_f,"
            << "dhail_alfa_n,"
            << "dhail_alfa_q,"
            << "dhail_lambda,"
            << "dhail_vsedi_min,"
            << "dhail_vsedi_max,"
            // Ice
            << "dice_a_geo,"
            << "dice_b_geo,"
            << "dice_min_x,"
            << "dice_min_x_act,"
            << "dice_min_x_nuc_homo,"
            << "dice_min_x_nuc_hetero,"
            << "dice_min_x_melt,"
            << "dice_min_x_evap,"
            << "dice_min_x_freezing,"
            << "dice_min_x_depo,"
            << "dice_min_x_collision,"
            << "dice_min_x_collection,"
            << "dice_min_x_conversion,"
            << "dice_min_x_sedimentation,"
            << "dice_min_x_riming,"
            << "dice_max_x,"
            << "dice_sc_theta_q,"
            << "dice_sc_delta_q,"
            << "dice_sc_theta_n,"
            << "dice_sc_delta_n,"
            << "dice_s_vel,"
            << "dice_a_vel,"
            << "dice_b_vel,"
            << "dice_rho_v,"
            << "dice_c_z,"
            << "dice_sc_coll_n,"
            << "dice_cmu0,"
            << "dice_cmu1,"
            << "dice_cmu2,"
            << "dice_cmu3,"
            << "dice_cmu4,"
            << "dice_cmu5,"
            << "dice_alpha,"
            << "dice_beta,"
            << "dice_gamma,"
            << "dice_nu,"
            << "dice_g1,"
            << "dice_g2,"
            << "dice_mu,"
            << "dice_nm1,"
            << "dice_nm2,"
            << "dice_nm3,"
            << "dice_q_crit_c,"
            << "dice_d_crit_c,"
            << "dice_ecoll_c,"
            << "dice_cap,"
            << "dice_a_ven,"
            << "dice_b_ven,"
            << "dice_c_s,"
            << "dice_a_f,"
            << "dice_b_f,"
            << "dice_alfa_n,"
            << "dice_alfa_q,"
            << "dice_lambda,"
            << "dice_vsedi_min,"
            << "dice_vsedi_max,"
            // Snow
            << "dsnow_a_geo,"
            << "dsnow_b_geo,"
            << "dsnow_min_x,"
            << "dsnow_min_x_act,"
            << "dsnow_min_x_nuc_homo,"
            << "dsnow_min_x_nuc_hetero,"
            << "dsnow_min_x_melt,"
            << "dsnow_min_x_evap,"
            << "dsnow_min_x_freezing,"
            << "dsnow_min_x_depo,"
            << "dsnow_min_x_collision,"
            << "dsnow_min_x_collection,"
            << "dsnow_min_x_conversion,"
            << "dsnow_min_x_sedimentation,"
            << "dsnow_min_x_riming,"
            << "dsnow_max_x,"
            << "dsnow_sc_theta_q,"
            << "dsnow_sc_delta_q,"
            << "dsnow_sc_theta_n,"
            << "dsnow_sc_delta_n,"
            << "dsnow_s_vel,"
            << "dsnow_a_vel,"
            << "dsnow_b_vel,"
            << "dsnow_rho_v,"
            << "dsnow_c_z,"
            << "dsnow_sc_coll_n,"
            << "dsnow_cmu0,"
            << "dsnow_cmu1,"
            << "dsnow_cmu2,"
            << "dsnow_cmu3,"
            << "dsnow_cmu4,"
            << "dsnow_cmu5,"
            << "dsnow_alpha,"
            << "dsnow_beta,"
            << "dsnow_gamma,"
            << "dsnow_nu,"
            << "dsnow_g1,"
            << "dsnow_g2,"
            << "dsnow_mu,"
            << "dsnow_nm1,"
            << "dsnow_nm2,"
            << "dsnow_nm3,"
            << "dsnow_q_crit_c,"
            << "dsnow_d_crit_c,"
            << "dsnow_ecoll_c,"
            << "dsnow_cap,"
            << "dsnow_a_ven,"
            << "dsnow_b_ven,"
            << "dsnow_c_s,"
            << "dsnow_a_f,"
            << "dsnow_b_f,"
            << "dsnow_alfa_n,"
            << "dsnow_alfa_q,"
            << "dsnow_lambda,"
            << "dsnow_vsedi_min,"
            << "dsnow_vsedi_max"
            << "\n";
    } // End loop over all components

    // CODIPACK: END
    std::cout << "Various options:\n"
              << "----------------\n"
              << "Output name:\t" << input.OUTPUT_FILENAME << "\n"
              << "start_over:\t" << input.start_over << "\n"
              << "fixed_iter:\t" << input.fixed_iteration << "\n\n";

    // Loop for timestepping: BEGIN
    // ==================================================
    // Read the trajectory file
    // ==================================================
    try
    {

        int ncid;
        nc_open(global_args.input_file, NC_NOWRITE, &ncid);

        // Loop over every timestep that is fixed to 20 s
        std::vector<size_t> startp, countp;
#if defined WCB || defined WCB2
        startp.push_back(1);          // time
        startp.push_back(input.traj); // trajectory
#else
        startp.push_back(input.traj); // trajectory
        startp.push_back(1);          // time
#endif

        countp.push_back(1);
        countp.push_back(1);
        codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
        for(uint32_t t=0; t<cc.num_steps; ++t) //
        {
            // bool written = false;
#if defined WCB || defined WCB2
            startp[0] = t;
#else
            startp[1] = t+1; //  * nc_params.n_trajectories
#endif

            netCDF::NcFile datafile(global_args.input_file, netCDF::NcFile::read);
            load_nc_parameters_var(nc_params, datafile);
            load_nc_parameters(nc_params, startp, countp,
                               ref_quant, cc.num_sub_steps);

            // Set values from a given trajectory
            if(t==0 || input.start_over)
            {
                y_single_old[p_idx]  = nc_params.p;     // p
                y_single_old[T_idx]  = nc_params.t;     // T
                y_single_old[S_idx]  = nc_params.S;     // S
                y_single_old[qc_idx] = nc_params.qc;    // qc
                y_single_old[qr_idx] = nc_params.qr;    // qr
                y_single_old[qv_idx] = nc_params.qv;    // qv
                y_single_old[qi_idx] = nc_params.qi;    // qi
                y_single_old[qs_idx] = nc_params.qs;    // qs
#if !defined(WCB)
                y_single_old[qg_idx] = nc_params.qg;    // qg
#else
                if(t==0)
                    y_single_old[qg_idx] = 0;
#endif

                if(t==0)
                {
                    y_single_old[qh_idx] = 0.0; // qh. We don't have hail in the trajectoris
                    y_single_old[Nh_idx] = 0.0; // Nh. We don't have hail in the trajectoris
                }
                codi::RealReverse denom = 0;
#ifdef WCB2
                y_single_old[Nc_idx] = nc_params.Nc;
                y_single_old[Nr_idx] = nc_params.Nr;
                y_single_old[Ng_idx] = nc_params.Ng;
                y_single_old[Ni_idx] = nc_params.Ni;
                y_single_old[Ns_idx] = nc_params.Ns;

                y_single_old[Nr_out_idx] = nc_params.NRout;
                y_single_old[Ng_out_idx] = nc_params.NGout;
                y_single_old[Ni_out_idx] = nc_params.NIout;
                y_single_old[Ns_out_idx] = nc_params.NSout;

#else
                denom = (cc.cloud.max_x - cc.cloud.min_x) / 2.0 + cc.cloud.min_x;
                y_single_old[Nc_idx] = y_single_old[qc_idx] * ref_quant.qref / (denom); //*10e2);  // Nc
                denom = (cc.rain.max_x - cc.rain.min_x) / 2 + cc.rain.min_x;
                y_single_old[Nr_idx] = y_single_old[qr_idx] * ref_quant.qref / (denom); //*10e2);  // Nr
                denom = cc.cloud.min_x / 2.0;
                y_single_old[Nv_idx] = y_single_old[qv_idx] * ref_quant.qref / (denom); //*10e2);  // Nv
                denom = (cc.ice.max_x - cc.ice.min_x) / 2.0 + cc.ice.min_x;
                y_single_old[Ni_idx] = y_single_old[qi_idx] * ref_quant.qref / (denom); //*10e2); // Ni
                denom = (cc.snow.max_x - cc.snow.min_x) / 2.0 + cc.snow.min_x;
                y_single_old[Ns_idx] = y_single_old[qs_idx] * ref_quant.qref / (denom); //*10e2); // Ns
                denom = (cc.graupel.max_x - cc.graupel.min_x) / 2.0 + cc.graupel.min_x;
                y_single_old[Ng_idx] = y_single_old[qg_idx] * ref_quant.qref / (denom); //*10e2); // Ng
#endif
                cc.Nc_prime = y_single_old[Nc_idx];

                cc.rho_a_prime = compute_rhoa(nc_params.p*ref_quant.pref,//*100,
                    nc_params.t*ref_quant.Tref, nc_params.S);
                y_single_old[w_idx]  = nc_params.w[0]; // w
                cc.dw = nc_params.dw / (cc.dt*cc.num_sub_steps);

                denom = cc.cloud.min_x / 2.0;
                y_single_old[Nv_idx] = y_single_old[qv_idx] * ref_quant.qref / (denom); //*10e2);  // Nv

                y_single_old[z_idx] = nc_params.z[0];

#if defined WCB || defined WCB2
                out_tmp << (t*cc.num_sub_steps)*cc.dt << "," << input.traj << ","
                        << nc_params.lon[0] << "," << nc_params.lat[0] << ","
                        << nc_params.ascent_flag << ",";
#else
                out_tmp << (t*cc.num_sub_steps)*cc.dt << "," << input.traj << ","
                        << nc_params.lon[0] << "," << nc_params.lat[0] << ",";
#endif
                for(int ii = 0 ; ii < num_comp; ii++)
                    out_tmp << y_single_old[ii] <<
                        ((ii == num_comp-1) ? "\n" : ",");


                for(int ii = 0 ; ii < num_comp ; ii++)
                {
#if defined WCB || defined WCB2
                    out_diff_tmp[ii] << t*cc.num_sub_steps*cc.dt << ","
                                    << input.traj << ","
                                    << output_par_idx[ii] << ","
                                    << nc_params.lon[0] << ","
                                    << nc_params.lat[0] << ","
                                    << nc_params.ascent_flag << ",";
#else
                    out_diff_tmp[ii] << t*cc.num_sub_steps*cc.dt << ","
                                    << input.traj << ","
                                    << output_par_idx[ii] << ","
                                    << nc_params.lon[0] << ","
                                    << nc_params.lat[0] << ",";
#endif
                    for(int jj = 0 ; jj < num_par ; jj++)
                        out_diff_tmp[ii] << 0.0
                            << ((jj==num_par-1) ? "\n" : ",");
                }

#if defined(FLUX) && !defined(WCB)
                inflow[qi_in_idx] = nc_params.QIin;
                inflow[qs_in_idx] = nc_params.QSin;
                inflow[qr_in_idx] = nc_params.QRin;
                inflow[qg_in_idx] = nc_params.QGin;
#else
                inflow[qi_in_idx] = 0;
                inflow[qs_in_idx] = 0;
                inflow[qr_in_idx] = 0;
                inflow[qg_in_idx] = 0;
#endif
#if defined(FLUX) && defined(WCB2)
                inflow[Ni_in_idx] = nc_params.NIin;
                inflow[Ns_in_idx] = nc_params.NSin;
                inflow[Nr_in_idx] = nc_params.NRin;
                inflow[Ng_in_idx] = nc_params.NGin;
#else
                inflow[Ni_in_idx] = 0;
                inflow[Ns_in_idx] = 0;
                inflow[Nr_in_idx] = 0;
                inflow[Ng_in_idx] = 0;
#endif
            }
            // Iterate over each substep
            for(uint32_t sub=1; sub<cc.num_sub_steps; ++sub) // cc.num_sub_steps
            {
#if defined(TRACE_QR) || defined(TRACE_QV) || defined(TRACE_QC) || defined(TRACE_QI) || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
                std::cout << "\n\ntimestep : " << (sub*cc.dt_prime + t*cc.num_sub_steps*cc.dt_prime) << "\n";
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
                                                // hidden in coefficient c
                                                // Needed for condensation

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

                // ==================================================
                // Output if needed
                // ==================================================
                if( (0 == (sub + t*cc.num_sub_steps) % input.snapshot_index)
                   || ( t == cc.num_steps-1 && sub == cc.num_sub_steps-1 ) )
                {
                    // Write the results to the output file
#if defined WCB || defined WCB2
                    out_tmp << time_new << "," << input.traj << ","
                            << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                            << (nc_params.lat[0] + sub*nc_params.dlat) << ","
                            << nc_params.ascent_flag << ",";
#else
                    out_tmp << time_new << "," << input.traj << ","
                            << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                            << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#endif
                    for(int ii = 0 ; ii < num_comp; ii++)
                        out_tmp << y_single_new[ii]
                            << ((ii == num_comp-1) ? "\n" : ",");

                    // CODIPACK: BEGIN
                    for(int ii = 0 ; ii < num_comp ; ii++)
                    {
#if defined WCB || defined WCB2
                        out_diff_tmp[ii] << time_new << "," << input.traj << ","
                                     << output_par_idx[ii] << ","
                                     << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                                     << (nc_params.lat[0] + sub*nc_params.dlat) << ","
                                     << nc_params.ascent_flag << ",";
#else
                        out_diff_tmp[ii] << time_new << "," << input.traj << ","
                                     << output_par_idx[ii] << ","
                                     << (nc_params.lon[0] + sub*nc_params.dlon) << ","
                                     << (nc_params.lat[0] + sub*nc_params.dlat) << ",";
#endif
                        for(int jj = 0 ; jj < num_par ; jj++)
                            out_diff_tmp[ii] << y_diff[ii][jj]
                                << ((jj==num_par-1) ? "\n" : ",");
                    }
                    // CODIPACK: END
                }
                if( (0 == (sub + t*cc.num_sub_steps) % input.write_index)
                    || ( t == cc.num_steps-1 && sub == cc.num_sub_steps-1 ) )
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

                // ==================================================
                // Interchange old and new for next step
                // ==================================================
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

    // CODIPACK: BEGIN
    for(int ii = 0 ; ii < num_comp ; ii++)
        out_diff[ii].close();
    // CODIPACK: END

    std::cout << "-------------------FINISHED-------------------\n\n\n";
}
