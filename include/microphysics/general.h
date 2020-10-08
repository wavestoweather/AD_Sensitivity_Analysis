#ifndef GENERAL_H
#define GENERAL_H

#include <stdio.h>
#include "codi.hpp"
#include "constants.h"


/** @defgroup general General Purpose Functions
 * This file contains several general purpose functions.
 * @{
 */

/**
 * This function prints the contents of a vector y with length
 * num_comp to stdout.
 *
 * @param y Contents to print of size num_comp
 */
void printy(const double y[])
{
  std::cout << "\n";
  for(int ii = 0 ; ii < num_comp ; ii++){
    std::cout << "y[" << std::to_string(ii) << "] = " << y[ii] << "\n";
  }
  std::cout << std::endl;
}

/**
 * This function computes
 * \f[ \vec{y} = \vec{v1} + a*\vec{v2} \f]
 * for vectors \f$\vec{y}, \vec{v1}, \vec{v2}\f$ and a scalar \f$a\f$.
 *
 * @param y On out: the result of \f4\vec{y} = \vec{v1} + a*\vec{v2}\f$
 * @param v1 \f$\vec{v1}\f$
 * @param v2 \f$\vec{v2}\f$
 * @param a \f$a\f$
 */
void v1pav2(double y[],
	    const double v1[],
	    const double v2[],
	    const double a)
{

  for(int ii = 0 ; ii < num_comp ; ii++){
    y[ii] = v1[ii] + a*v2[ii];
  }

}

/**
 * Print parameters for given particle.
 *
 * @param pc Particle constants
 * @param title Name of particle
 */
void print_particle_params(
    particle_model_constants_t &pc,
    std::string title)
{
    std::cout << title << "\n"
        << "a_geo = " << pc.a_geo << "\n"
        << "b_geo = " << pc.b_geo << "\n"
        << "min_x = " << pc.min_x << "\n"
        << "max_x = " << pc.max_x << "\n"
        << "sc_theta_q = " << pc.sc_theta_q << "\n"
        << "sc_delta_q = " << pc.sc_delta_q << "\n"
        << "sc_theta_n = " << pc.sc_theta_n << "\n"
        << "sc_delta_n = " << pc.sc_delta_n << "\n"
        << "s_vel = " << pc.s_vel << "\n"
        << "a_vel = " << pc.a_vel << "\n"
        << "b_vel = " << pc.b_vel << "\n"
        << "rho_v = " << pc.rho_v << "\n"
        << "c_z = " << pc.c_z << "\n"
        << "sc_coll_n = " << pc.sc_coll_n << "\n"
        << "nu = " << pc.nu << "\n"
        << "mu = " << pc.mu << "\n"
        << "q_crit_c = " << pc.q_crit_c << "\n"
        << "d_crit_c = " << pc.d_crit_c << "\n"
        << "ecoll_c = " << pc.ecoll_c << "\n"
        << "cap = " << pc.cap << "\n"
        << "a_ven = " << pc.a_ven << "\n"
        << "b_ven = " << pc.b_ven << "\n"
        << "c_s = " << pc.c_s << "\n"
        << "a_f = " << pc.a_f << "\n"
        << "b_f = " << pc.b_f << "\n"
        << "b_f (COSMO variant) = " << pc.b_f / pow(N_sc, n_f) * sqrt(kin_visc_air) << "\n"
        << "alfa_n = " << pc.alfa_n << "\n"
        << "alfa_q = " << pc.alfa_q << "\n"
        << "lambda = " << pc.lambda << "\n"
        << "vsedi_min = " << pc.vsedi_min << "\n"
        << "vsedi_max = " << pc.vsedi_max << "\n\n";
}

/**
 * Print reference quantities temperature, pressure,
 * mixing-ratio, vertical velocity, time.
 *
 * @param ref Struct with reference quantities.
 */
void print_reference_quantities(reference_quantities_t &ref)
{
  std::cout << "\nReference quantities\n"
	    << "--------------------\n"
        << "Temperature: " << ref.Tref << " Kelvin\n"
        << "Pressure: " << ref.pref << " Pascal\n"
        << "Mixing-ratio: " << ref.qref << "\n"
        << "Vertical velocity: " << ref.wref << " meter per second\n"
        << "Time: " << ref.tref << " Second\n"
        << std::endl << std::flush;
}

/**
 * Print constants given a model constants structure, namely integration
 * time, number of steps, scaling factors.
 *
 * @param cc A structure with model constants.
 */
void print_constants(model_constants_t &cc)
{

  std::cout << "\nModel constants:\n"
	    << "----------------\n"
	    << "Final integration time: " << cc.t_end_prime << " seconds\n"
        << "Nondimensional final integration time: " << cc.t_end << "\n"
	    << "Timestep: " << cc.dt_prime << " seconds\n"
	    << "Nondimensional timestep: " << cc.dt << "\n"
	    << "Number of iterations: " << cc.num_steps << "\n"
        << "Number of substeps: " << cc.num_sub_steps << "\n"
	    << "a1_scale: " << cc.a1_scale << "\n"
        << "a2_scale: " << cc.a2_scale << "\n"
        << "e1_scale: " << cc.e1_scale << "\n"
        << "e2_scale: " << cc.e2_scale << "\n"
        << "d_scale: " << cc.d_scale << "\n"
	    << "Scaling factor: " << cc.scaling_fact << "\n"
	    << std::endl << std::flush;
}

/**
 * Print all given input parameters
 *
 * @param in Struct where the input parameters are stored.
 */
void print_input_parameters(input_parameters_t &in)
{
  std::cout << "\n"
	    << "Technical input parameters:\n"
	    << "---------------------------\n"
        << "Time to integrate: " << in.t_end_prime << " Seconds\n"
        << "Timestep: " << in.dt_prime << " Seconds\n"
#ifdef MET3D
        << "Start time (relative to ascend): " << in.start_time << " Seconds\n"
#endif
	    << "Snapshot index: " << in.snapshot_index << "\n"
        << "Write index: " << in.write_index << "\n"
        << "Progressbar index: " << in.progress_index << "\n"
        << "Name of output file: " << in.OUTPUT_FILENAME << "\n"
	    << "Scaling factor: " << in.scaling_fact << "\n"
        << "Name of input file: " << in.INPUT_FILENAME << "\n"
        << "Start over mixing ratios and particle numbers at each timestep of a trajectory?: " << in.start_over << "\n"
        << "Start over pressure, temperature and ascent at each timestep of a trajectory?: " << in.start_over_env << "\n"
        << "Fix temperature and pressure at each substep?: " << in.fixed_iteration << "\n"
        << "Auto type for rain evaporation (1, 2, 3): " << in.auto_type << "\n"
        << "Trajectory used: " << in.traj
        << std::endl << std::flush;
}

/**
 * Display a help message on how to use this program.
 */
void display_usage()
{

  std::cout << "\n"
	    << "USAGE of the program in water cloud mode:\n"
	    << "Invoke the program on the command line with\n"
	    << "$ ./main PARAMETERS\n"
	    << "where PARAMETERS are the following parameters:\n"
	    << "-f: Time to integrate in seconds. "
        << "If higher than number of entries from input file, then it stops early\n"
	    << "-i: Snapshot index.\n"
	    << "-b: Scaling factor.\n"
	    << "-d: Timestep in seconds for a substep between each new trajectory input.\n"
	    << "-o: Name of the output file.\n"
        << "-l: Path and name of input file.\n"
        << "-s: Start over mixing ratios and particle numbers at each timestep of a trajectory.\n"
        << "-e: Start over temperature, pressure and ascent at each timestep of a trajectory.\n"
        << "-t: Set pressure, temperature and vertical velocity fixed at each substep.\n"
        << "-a: Set auto type (1=KB, 2=KK, 3=SB).\n"
        << "-r: Set index of trajectory.\n"
        << "-w: Write index for the snapshots.\n"
        << "-p: Index for updating the progressbar.\n"
#ifdef MET3D
        << "-n: No simulation until given time (relative to ascend).\n"
#endif
	    << "-?: This help message.\n"
	    << std::endl;
}

/**
 * Display an error message when command line arguments are faulty.
 */
void display_error_on_command_line()
{
  std::cout << "==> ERROR: An error occured while dealing with the command line arguments!"
	    << std::endl;
}
/** @} */ // end of group general

#endif
