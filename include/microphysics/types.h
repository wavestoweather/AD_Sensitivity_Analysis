#pragma once
#include "codi.hpp"
#include <netcdf>
using namespace netCDF;

////////////////////////////////////////////////////////////////////////////////
// Collect all reference quantities into this structure.
////////////////////////////////////////////////////////////////////////////////
struct reference_quantities_t{
  double Tref;			// Reference temperature
  double pref;			// Reference pressure
  double qref;			// Reference mixing-ratio
  double wref;			// Reference vertical velocity
  double tref;			// Reference time

  double Nref;			// DUMMY
};

struct particle_model_constants_t{
  codi::RealReverse a_geo;
  codi::RealReverse b_geo;
  codi::RealReverse min_x;
  codi::RealReverse max_x;
  codi::RealReverse sc_theta_q;
  codi::RealReverse sc_delta_q;
  codi::RealReverse sc_theta_n;
  codi::RealReverse sc_delta_n;
  codi::RealReverse s_vel; // sigma_vel
  codi::RealReverse a_vel;
  codi::RealReverse b_vel;
  codi::RealReverse rho_v;
  codi::RealReverse c_z;
  codi::RealReverse sc_coll_n; // e_coll
  codi::RealReverse cmu0, cmu1, cmu2, cmu3, cmu4, cmu5, alpha, beta, gamma;
  codi::RealReverse nu;

  codi::RealReverse g1, g2, mu, nm1, nm2, nm3;
  // Parameters for riming
  codi::RealReverse q_crit_c;
  codi::RealReverse d_crit_c;
  codi::RealReverse ecoll_c;
  codi::RealReverse cap; // coefficient for capacity of particle
  codi::RealReverse a_ven;
  codi::RealReverse b_ven;

  codi::RealReverse c_s; // c_i
  codi::RealReverse a_f;
  codi::RealReverse b_f;

  // Bulk sedimentation velocity
  codi::RealReverse alfa_n;
  codi::RealReverse alfa_q;
  codi::RealReverse lambda;
  codi::RealReverse vsedi_min;
  codi::RealReverse vsedi_max;

  void register_input(codi::RealReverse::TapeType &tape)
  {
      tape.registerInput(this->a_geo);
      tape.registerInput(this->b_geo);
      tape.registerInput(this->min_x);
      tape.registerInput(this->max_x);
      tape.registerInput(this->sc_theta_q);
      tape.registerInput(this->sc_delta_q);
      tape.registerInput(this->sc_theta_n);
      tape.registerInput(this->sc_delta_n);
      tape.registerInput(this->s_vel);
      tape.registerInput(this->a_vel);
      tape.registerInput(this->b_vel);
      tape.registerInput(this->rho_v);
      tape.registerInput(this->c_z);
      tape.registerInput(this->sc_coll_n);
      tape.registerInput(this->cmu0);
      tape.registerInput(this->cmu1);
      tape.registerInput(this->cmu2);
      tape.registerInput(this->cmu3);
      tape.registerInput(this->cmu4);
      tape.registerInput(this->cmu5);
      tape.registerInput(this->alpha);
      tape.registerInput(this->beta);
      tape.registerInput(this->gamma);
      tape.registerInput(this->nu);
      tape.registerInput(this->g1);
      tape.registerInput(this->g2);
      tape.registerInput(this->mu);
      tape.registerInput(this->nm1);
      tape.registerInput(this->nm2);
      tape.registerInput(this->nm3);
      tape.registerInput(this->q_crit_c);
      tape.registerInput(this->d_crit_c);
      tape.registerInput(this->ecoll_c);
      tape.registerInput(this->cap);
      tape.registerInput(this->a_ven);
      tape.registerInput(this->b_ven);
      tape.registerInput(this->c_s);
      tape.registerInput(this->a_f);
      tape.registerInput(this->b_f);
      tape.registerInput(this->alfa_n);
      tape.registerInput(this->alfa_q);
      tape.registerInput(this->lambda);
      tape.registerInput(this->vsedi_min);
      tape.registerInput(this->vsedi_max);
  }

  template<class T>
  void get_gradient(T &out_vec, uint64_t &idx)
  {
      out_vec[idx] = this->a_geo.getGradient();
      idx++;
      out_vec[idx] = this->b_geo.getGradient();
      idx++;
      out_vec[idx] = this->min_x.getGradient();
      idx++;
      out_vec[idx] = this->max_x.getGradient();
      idx++;
      out_vec[idx] = this->sc_theta_q.getGradient();
      idx++;
      out_vec[idx] = this->sc_delta_q.getGradient();
      idx++;
      out_vec[idx] = this->sc_theta_n.getGradient();
      idx++;
      out_vec[idx] = this->sc_delta_n.getGradient();
      idx++;
      out_vec[idx] = this->s_vel.getGradient();
      idx++;
      out_vec[idx] = this->a_vel.getGradient();
      idx++;
      out_vec[idx] = this->b_vel.getGradient();
      idx++;
      out_vec[idx] = this->rho_v.getGradient();
      idx++;
      out_vec[idx] = this->c_z.getGradient();
      idx++;
      out_vec[idx] = this->sc_coll_n.getGradient();
      idx++;
      out_vec[idx] = this->cmu0.getGradient();
      idx++;
      out_vec[idx] = this->cmu1.getGradient();
      idx++;
      out_vec[idx] = this->cmu2.getGradient();
      idx++;
      out_vec[idx] = this->cmu3.getGradient();
      idx++;
      out_vec[idx] = this->cmu4.getGradient();
      idx++;
      out_vec[idx] = this->cmu5.getGradient();
      idx++;
      out_vec[idx] = this->alpha.getGradient();
      idx++;
      out_vec[idx] = this->beta.getGradient();
      idx++;
      out_vec[idx] = this->gamma.getGradient();
      idx++;
      out_vec[idx] = this->nu.getGradient();
      idx++;
      out_vec[idx] = this->g1.getGradient();
      idx++;
      out_vec[idx] = this->g2.getGradient();
      idx++;
      out_vec[idx] = this->mu.getGradient();
      idx++;
      out_vec[idx] = this->nm1.getGradient();
      idx++;
      out_vec[idx] = this->nm2.getGradient();
      idx++;
      out_vec[idx] = this->nm3.getGradient();
      idx++;
      out_vec[idx] = this->q_crit_c.getGradient();
      idx++;
      out_vec[idx] = this->d_crit_c.getGradient();
      idx++;
      out_vec[idx] = this->ecoll_c.getGradient();
      idx++;
      out_vec[idx] = this->cap.getGradient();
      idx++;
      out_vec[idx] = this->a_ven.getGradient();
      idx++;
      out_vec[idx] = this->b_ven.getGradient();
      idx++;
      out_vec[idx] = this->c_s.getGradient();
      idx++;
      out_vec[idx] = this->a_f.getGradient();
      idx++;
      out_vec[idx] = this->b_f.getGradient();
      idx++;
      out_vec[idx] = this->alfa_n.getGradient();
      idx++;
      out_vec[idx] = this->alfa_q.getGradient();
      idx++;
      out_vec[idx] = this->lambda.getGradient();
      idx++;
      out_vec[idx] = this->vsedi_min.getGradient();
      idx++;
      out_vec[idx] = this->vsedi_max.getGradient();
      idx++;
  }
};

struct collection_model_constants_t{
    codi::RealReverse delta_n_aa;
    codi::RealReverse delta_n_ab;
    codi::RealReverse delta_n_bb;
    codi::RealReverse delta_q_aa;
    codi::RealReverse delta_q_ab;
    codi::RealReverse delta_q_bb;
    codi::RealReverse delta_q_ba;

    codi::RealReverse theta_n_aa;
    codi::RealReverse theta_n_ab;
    codi::RealReverse theta_n_bb;
    codi::RealReverse theta_q_aa;
    codi::RealReverse theta_q_ab;
    codi::RealReverse theta_q_bb;
    codi::RealReverse theta_q_ba;
};

// A 4D lookup table
struct table_t{
    uint64_t n1, n2, n3, n4; // Number of grid points in every direction
    std::vector<codi::RealReverse> x1, x2, x3, x4; // Grid vectors
    codi::RealReverse dx1, dx2, dx3, dx4; // Grid distances for every vector
    codi::RealReverse odx1, odx2, odx3, odx4; // One over dx 1, 2, 3, 4
    std::vector<codi::RealReverse> table;

    codi::RealReverse get(
        uint64_t i,
        uint64_t j,
        uint64_t k,
        uint64_t l) const
    {
        return table[i*n3*n2*n1 + j*n2*n1 + k*n1 + l];
    }
};

struct gamma_table_t{
    uint64_t n_bins;
    uint64_t n_bins_highres;
    double a;
    std::vector<double> x; // always starts at 0 and has distant dx
    std::vector<double> x_highres; // as x but with higher resolution
    double dx;
    double dx_highres;
    double odx; // 1/dx
    double odx_highres;
    std::vector<double> igf; // values of the inc. gamma function at (a,x)
    std::vector<double> igf_highres;

    template <class A>
    A look_lo(A x)
    {
        A xt = max(min(x, this->x[this->n_bins-1]), 0.0);
        double tmp = xt.getValue()*this->odx;
        tmp = floor(tmp) + 1;
        uint64_t iu = std::min((uint64_t) tmp, this->n_bins-2);
        uint64_t io = iu + 1;
        return this->igf[iu] + (this->igf[io] - this->igf[iu]) * this->odx*(xt-this->x[iu]);
    }
};

////////////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////////////
struct model_constants_t{
  //
  // Physical constants warm cloud
  //
  // Accomodation coefficient
  double alpha_d;

  // Density of dry air
  double rho_a_prime;

  // Number concentration of cloud droplets
  codi::RealReverse Nc_prime;

  // Dimensional coefficients
  codi::RealReverse a1_prime;
  codi::RealReverse a2_prime;
  codi::RealReverse e1_prime;
  codi::RealReverse e2_prime;
  codi::RealReverse B_prime;
  codi::RealReverse d_prime;

  // Change in buoancy
  codi::RealReverse dw;

  // Exponents
  codi::RealReverse gamma;
  codi::RealReverse betac;
  codi::RealReverse betar;
  codi::RealReverse delta1;
  codi::RealReverse delta2;
  codi::RealReverse zeta;

  // 2 moment stuff
  particle_model_constants_t hail;
  particle_model_constants_t ice;
  particle_model_constants_t snow;
  particle_model_constants_t cloud;
  particle_model_constants_t rain;
  particle_model_constants_t graupel;
  collection_model_constants_t coeffs_scr;
  collection_model_constants_t coeffs_srr;
  collection_model_constants_t coeffs_irr;
  collection_model_constants_t coeffs_icr;
  collection_model_constants_t coeffs_hrr;
  collection_model_constants_t coeffs_grr;
  collection_model_constants_t coeffs_hcr;
  collection_model_constants_t coeffs_gcr;
  collection_model_constants_t coeffs_sic;
  collection_model_constants_t coeffs_hic;
  collection_model_constants_t coeffs_gic;
  collection_model_constants_t coeffs_hsc;
  collection_model_constants_t coeffs_gsc;

  //
  // Technical constants
  //
  double T_end_prime;
  double T_end;
  double dt_prime;
  double dt;
  double dt_traject_prime;
  double dt_traject;
  uint64_t num_steps;
  uint64_t num_sub_steps;
  uint64_t snapshot_index;

  //
  // General performance constants
  //
  double dt_half;
  double dt_sixth;
  double dt_third;

  //
  // Performance constants warm cloud
  //
  codi::RealReverse a1_scale;
  codi::RealReverse a2_scale;
  codi::RealReverse e1_scale;
  codi::RealReverse e2_scale;
  codi::RealReverse d_scale;

  //
  // More constants for the IFS model
  //
  const double nar = 0.22;
  const double nbr = 2.2;
  const double ar = M_PI / 6.0;
  const double br = 3.0;
  const double cr = 386.8;
  const double dr = 0.67;
  const double Sc = 0.6;
  const double mu = 16.0e-6;
  const double rho0 = 1.0; // 1.225?

  const double alphar = 1.0/(br + 1.0 - nbr);
  const double epsilonr = 0.5*dr + 2.5 - nbr;

  // Scaling factor
  double scaling_fact;

};

////////////////////////////////////////////////////////////////////////////////
// Structure to collect all nc parameters.
////////////////////////////////////////////////////////////////////////////////
struct nc_parameters_t{

    uint32_t n_trajectories = 30;
    uint32_t n_timesteps = 7922;
    std::vector<double> w, z;
    double lat, lon, t, p, time_rel,
            qc, qr, qi, qs, qg, qv, S, dw,
            QIin, QSin, QRin, QGin, QIout, QSout, QRout, QGout;
    NcVar lat_var, lon_var, z_var, t_var, p_var, w_var, time_rel_var,
            qc_var, qr_var, qi_var, qs_var, qg_var, qv_var, S_var,
            QIin_var, QSin_var, QRin_var, QGin_var, QIout_var, QSout_var,
            QRout_var, QGout_var;
};