#ifndef USER_FUNCTIONS_H
#define USER_FUNCTIONS_H

#include "codi.hpp"
#include <math.h>
#include <cmath>
#include <tgmath.h>
#include <vector>
#include "constants.h"
#include "physical_parameterizations.h"
#include "types.h"


////////////////////////////////////////////////////////////////////////////////
// =============================================================================
// This file provides the functions to be defined by the user,
// specifically the definition of the right-hand side function RHS of
// the ODE
// y' = RHS(y)
// which is solved.
// =============================================================================
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// This function sets the constants for the cloud model from given
// environmental conditions.
//
void setCoefficients(
    std::vector<codi::RealReverse> & y,
    model_constants_t& cc,
    reference_quantities_t& ref,
    nc_parameters_t& nc)
{
    codi::RealReverse p_prime = y[p_idx]*ref.pref;
    codi::RealReverse T_prime = y[T_idx]*ref.Tref;

    codi::RealReverse rho_prime = p_prime /( Ra * T_prime );
    codi::RealReverse L_vap_prime = latent_heat_water(T_prime);
    codi::RealReverse Ka_prime = thermal_conductivity_dry_air(T_prime);
    codi::RealReverse psat_prime = saturation_pressure_water(T_prime);
    codi::RealReverse A_pp = (L_vap_prime/(Ka_prime*T_prime))*((L_vap_prime/(Rv*T_prime)) - 1.0);
    codi::RealReverse B_pp = (Rv*T_prime)/((2.21/p_prime)*psat_prime);


    cc.e1_prime = cc.e1_scale * ( pow(rho_prime, 2.0*cc.alphar-2.0)/(A_pp + B_pp) );
    cc.e2_prime = cc.e2_scale * ( pow(rho_prime, cc.alphar*cc.epsilonr - (7.0/4.0))/(A_pp + B_pp) );

    cc.a1_prime = cc.a1_scale;	// Constant coefficient
    cc.a2_prime = cc.a2_scale;	// Constant coefficient
    cc.d_prime = cc.d_scale;	// Constant coefficient
}

void setCoefficients(codi::RealReverse p_prime,
		     codi::RealReverse T_prime,
		     model_constants_t &cc)
{
  codi::RealReverse rho_prime = p_prime /( Ra * T_prime );
  codi::RealReverse L_vap_prime = latent_heat_water(T_prime);
  codi::RealReverse Ka_prime = thermal_conductivity_dry_air(T_prime);
  codi::RealReverse psat_prime = saturation_pressure_water(T_prime);
  codi::RealReverse A_pp = (L_vap_prime/(Ka_prime*T_prime))*((L_vap_prime/(Rv*T_prime)) - 1.0);
  codi::RealReverse B_pp = (Rv*T_prime)/((2.21/p_prime)*psat_prime);

  cc.a1_prime = cc.a1_scale;	// Constant coefficient
  cc.a2_prime = cc.a2_scale;	// Constant coefficient

  cc.e1_prime = cc.e1_scale * ( pow(rho_prime, 2.0*cc.alphar-2.0)/(A_pp + B_pp) );
  cc.e2_prime = cc.e2_scale * ( pow(rho_prime, cc.alphar*cc.epsilonr - (7.0/4.0))/(A_pp + B_pp) );

  cc.d_prime = cc.d_scale;	// Constant coefficient
}

////////////////////////////////////////////////////////////////////////////////
// This function evaluates the RHS function of the ODE.
//
void RHS(std::vector<codi::RealReverse> &res,
	 std::vector<codi::RealReverse> &y,
	 const reference_quantities_t &ref,
	 model_constants_t &cc,
     nc_parameters_t& nc)
{


  // ============================================================
  // This is the RHS for the water cloud
  // ============================================================

  // Storage:
  // y = [p, T, w, S, qc, qr, qv]

  // Decrypt the variables
  codi::RealReverse p = y[p_idx];
  codi::RealReverse T = y[T_idx];
  codi::RealReverse w = y[w_idx];

  codi::RealReverse S = y[S_idx];
  codi::RealReverse qc = y[qc_idx];
  codi::RealReverse qr = y[qr_idx];
  codi::RealReverse qv = y[qv_idx];

  // Safety measure: ensure positiveness
  if(0. > qc){
    qc = 0.0;
  }
  if(0. > qr){
    qr = 0.0;
  }

  if(0. > qv)
  {
    qv = 0.0;
  }
  // Change to dimensional variables
  codi::RealReverse p_prime = ref.pref * p;
  codi::RealReverse T_prime = ref.Tref * T;
  codi::RealReverse w_prime = ref.wref * w;
  codi::RealReverse qc_prime = ref.qref * qc;
  codi::RealReverse qr_prime = ref.qref * qr;
  codi::RealReverse qv_prime = ref.qref * qv;

  // Compute parameters
  codi::RealReverse psat_prime = saturation_pressure_water(T_prime);
  codi::RealReverse qsat_prime = Epsilon*( psat_prime/(p_prime - psat_prime) );
  codi::RealReverse p_div_T_prime = p_prime / T_prime;
  codi::RealReverse cpv_prime = specific_heat_water_vapor(T_prime);
  codi::RealReverse cpa_prime = specific_heat_dry_air(T_prime);
  codi::RealReverse cpl_prime = specific_heat_water(T_prime);
  codi::RealReverse rhow_prime = density_water(T_prime);
  codi::RealReverse L_prime = latent_heat_water(T_prime);
  codi::RealReverse H_prime =
    1.0/(((L_prime/(Rv*T_prime)) - 1.0)
    *(L_prime/(thermal_conductivity_dry_air(T_prime)*T_prime))
    + ((Rv*T_prime)/(cc.alpha_d*diffusivity(T_prime,p_prime)*psat_prime)));

  codi::RealReverse c_prime =
    4.0*M_PI*H_prime*pow((3.0/(4.0*M_PI*rhow_prime))
    *cc.Nc_prime*cc.Nc_prime , 1.0/3.0);

  codi::RealReverse qc_third = pow(qc , 1.0/3.0);
  codi::RealReverse qr_delta1 = pow(qr , cc.delta1);
  codi::RealReverse qr_delta2 = pow(qr , cc.delta2);
  codi::RealReverse qc_gamma = pow(qc , cc.gamma);
  codi::RealReverse qc_betac = pow(qc , cc.betac);
  codi::RealReverse qr_betar = pow(qr , cc.betar);
  codi::RealReverse qr_zeta = pow(qr , cc.zeta);


  // Compute nondimensional coefficients
  codi::RealReverse C1 = (ref.tref*gravity_acc*ref.wref)/(Ra*ref.Tref);
  codi::RealReverse C2 = ((1.0-Epsilon)*ref.qref)/Epsilon;
  codi::RealReverse C3 = (cpv_prime*ref.qref)/cpa_prime;
  codi::RealReverse C4 = (cpl_prime*ref.qref)/cpa_prime;
  codi::RealReverse C5 = (gravity_acc*ref.wref*ref.tref)/(ref.Tref*cpa_prime);
  codi::RealReverse C6 = (ref.tref*L_prime*c_prime*pow(ref.qref,1.0/3.0))/(ref.Tref*cpa_prime);
  codi::RealReverse C7 = (ref.tref*L_prime*cc.e1_prime
    *pow(ref.qref,cc.delta1))/(ref.Tref*cpa_prime);
  codi::RealReverse C8 = (ref.tref*L_prime*cc.e2_prime
    *pow(ref.qref,cc.delta2))/(ref.Tref*cpa_prime);
  codi::RealReverse C9 = (ref.tref*c_prime)/pow(ref.qref , 2.0/3.0);
  codi::RealReverse C10 = ref.tref*cc.a1_prime*pow(ref.qref , cc.gamma-1.0);
  codi::RealReverse C11 = ref.tref*cc.a2_prime*pow(ref.qref , cc.betac+cc.betar-1.0);
  codi::RealReverse C12 = ref.tref*cc.e1_prime*pow(ref.qref , cc.delta1-1.0);
  codi::RealReverse C13 = ref.tref*cc.e2_prime*pow(ref.qref , cc.delta2-1.0);
  codi::RealReverse C14 = ref.tref*cc.d_prime*pow(ref.qref , cc.zeta-1.0);
  codi::RealReverse C15 = Epsilon/ref.qref;
  codi::RealReverse C16 = L_prime/(Rv*ref.Tref);
  codi::RealReverse B = ref.tref*nc.QRin;

  //
  // Pressure
  //
  res[p_idx] = -( C1/(1.0 + C2*(qv/(1.0 + qv_prime))) )*( (p*w)/T );

  //
  // Temperature
  //
  res[T_idx] = ( 1.0/(1.0 + C3*qv + C4*(qc + qr)) )*( -C5*w + C6*qc_third*(S-1.0)
    + (C7*qr_delta1 + C8*qr_delta2)*min(S-1.0,0.0) );

  //
  // Vertical velocity. We do not really change anything here
  //
  res[w_idx] = nc.dw;

  //
  // Cloud droplets
  //
  res[qc_idx] = C9*qc_third*(S-1.0) - C10*qc_gamma - C11*qc_betac*qr_betar;

  //
  // Rain drops
  //
  res[qr_idx] = C10*qc_gamma + C11*qc_betac*qr_betar + (C12*qr_delta1 + C13*qr_delta2)
    *min(S-1.0,0.0) - B - C14*qr_zeta; // + B_out;

  //
  // Saturation ratio
  //
  res[S_idx] = (S/p)*res[p_idx] - (S/qv)*(1.0 - (qv/(C15+qv)))*(C9*qc_third*(S-1.0)
    + (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0,0.0)) - C16*(S/(T*T))*res[T_idx];

  // specific humidity
  res[qv_idx] = -C9*(S-1)*qc_third - (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0, 0.0);

//   std::cout << "S: " << S << ", C9: " << C9 << ", qc_third: " << qc_third
//     << ", C12: " << C12 << ", qr_delta1: " << qr_delta1
//     << ", C13: " << C13 << ", qr_delta2: " << qr_delta2 << "\n";

//   std::cout << "derivatives: p: " << res[p_idx] << ", T: " << res[T_idx]
//     << ", w: " << res[w_idx] << ", " << w << ", qc: " << res[qc_idx]
//     << ", qr: " << res[qr_idx] << ", S: " << res[S_idx]
//     << ", qv: " << res[qv_idx] << "\n";

} // End of RHS
////////////////////////////////////////////////////////////////////////////////

void Press_Temp(
    std::vector<codi::RealReverse> &res,
    std::vector<codi::RealReverse> &y,
    const reference_quantities_t &ref,
    model_constants_t &cc,
    nc_parameters_t& nc,
    const uint32_t traj,
    bool set_S=false)
{
  codi::RealReverse p = y[p_idx];
  codi::RealReverse T = y[T_idx];
  codi::RealReverse w = y[w_idx];

  codi::RealReverse S = y[S_idx];
  codi::RealReverse qc = y[qc_idx];
  codi::RealReverse qr = y[qr_idx];
  codi::RealReverse qv = y[qv_idx];

  if(0. > qv)
  {
    qv = 0.0;
  }
  // Change to dimensional variables
  codi::RealReverse p_prime = ref.pref * p;
  codi::RealReverse T_prime = ref.Tref * T;
  codi::RealReverse w_prime = ref.wref * w;
  codi::RealReverse qv_prime = ref.qref * qv;

  // Compute parameters
  codi::RealReverse psat_prime = saturation_pressure_water(T_prime);
  codi::RealReverse qsat_prime = Epsilon*( psat_prime/(p_prime - psat_prime) );
  codi::RealReverse p_div_T_prime = p_prime / T_prime;
  codi::RealReverse cpv_prime = specific_heat_water_vapor(T_prime);
  codi::RealReverse cpa_prime = specific_heat_dry_air(T_prime);
  codi::RealReverse cpl_prime = specific_heat_water(T_prime);
  codi::RealReverse rhow_prime = density_water(T_prime);
  codi::RealReverse L_prime = latent_heat_water(T_prime);
  codi::RealReverse H_prime =
    1.0/(((L_prime/(Rv*T_prime)) - 1.0)
    *(L_prime/(thermal_conductivity_dry_air(T_prime)*T_prime))
    + ((Rv*T_prime)/(cc.alpha_d*diffusivity(T_prime,p_prime)*psat_prime)));

  codi::RealReverse c_prime =
    4.0*M_PI*H_prime*pow((3.0/(4.0*M_PI*rhow_prime))
    *cc.Nc_prime*cc.Nc_prime , 1.0/3.0);

  codi::RealReverse qc_third = pow(qc , 1.0/3.0);
  codi::RealReverse qr_delta1 = pow(qr , cc.delta1);
  codi::RealReverse qr_delta2 = pow(qr , cc.delta2);

  // Compute nondimensional coefficients
  codi::RealReverse C1 = (ref.tref*gravity_acc*ref.wref)/(Ra*ref.Tref);
  codi::RealReverse C2 = ((1.0-Epsilon)*ref.qref)/Epsilon;
  codi::RealReverse C3 = (cpv_prime*ref.qref)/cpa_prime;
  codi::RealReverse C4 = (cpl_prime*ref.qref)/cpa_prime;
  codi::RealReverse C5 = (gravity_acc*ref.wref*ref.tref)/(ref.Tref*cpa_prime);
  codi::RealReverse C6 = (ref.tref*L_prime*c_prime*pow(ref.qref,1.0/3.0))/(ref.Tref*cpa_prime);
  codi::RealReverse C7 = (ref.tref*L_prime*cc.e1_prime
    *pow(ref.qref,cc.delta1))/(ref.Tref*cpa_prime);
  codi::RealReverse C8 = (ref.tref*L_prime*cc.e2_prime
    *pow(ref.qref,cc.delta2))/(ref.Tref*cpa_prime);
  codi::RealReverse B = ref.tref*nc.QRin; //cc.B_prime;

  //
  // Pressure
  //
  res[p_idx] = -( C1/(1.0 + C2*(qv/(1.0 + qv_prime))) )*( (p*w)/T );

  //
  // Temperature
  //
  res[T_idx] = ( 1.0/(1.0 + C3*qv + C4*(qc + qr)) )*( -C5*w + C6*qc_third*(S-1.0)
    + (C7*qr_delta1 + C8*qr_delta2)*min(S-1.0,0.0) );
}


// 2 moment cloud scheme https://doi.org/10.1007/s00703-005-0112-4
// After Seifert and Beheng (2006)
// Based on mo_art_2mom_driver.f90 and
// mo_art_2mom_main.f90 of ICON
void RHS_SB(std::vector<codi::RealReverse> &res,
    std::vector<codi::RealReverse> &y,
    const reference_quantities_t &ref,
    model_constants_t &cc,
    nc_parameters_t &nc,
    const double &dt,
    bool fixed=false,
    const int counter=0)
{
    // // Decrypt the variables
    codi::RealReverse p = y[p_idx];
    codi::RealReverse T = y[T_idx];
    codi::RealReverse w = y[w_idx];

    codi::RealReverse S = y[S_idx];
    codi::RealReverse qc = y[qc_idx];
    codi::RealReverse qr = y[qr_idx];
    codi::RealReverse qv = y[qv_idx];
    codi::RealReverse Nc = y[Nc_idx];
    codi::RealReverse Nr = y[Nr_idx];
    codi::RealReverse Nv = y[Nv_idx];
    codi::RealReverse qi = y[qi_idx];
    codi::RealReverse Ni = y[Ni_idx];
    codi::RealReverse vi = y[vi_idx];
    codi::RealReverse qs = y[qs_idx];
    codi::RealReverse Ns = y[Ns_idx];
    codi::RealReverse qg = y[qg_idx];
    codi::RealReverse Ng = y[Ng_idx];
    codi::RealReverse qh = y[qh_idx];
    codi::RealReverse Nh = y[Nh_idx];

    for(auto &r: res) r = 0;

    // Safety measure: ensure positiveness
    if(0. > qc)
        qc = 0.0;

    if(0. > qr)
        qr = 0.0;

    if(0. > qv)
        qv = 0.0;

    if(0. > Nc)
        Nc = 0.0;

    if(0. > qi)
        qi = 0.0;

    if(0. > qs)
        qs = 0.0;

    if(0. > qg)
        qg = 0.0;

    if(0. > qh)
        qh = 0.0;

    if(0. > Nr)
        Nr = 0.0;

    if(0. > Nv)
        Nv = 0.0;

    if(0. > Ni)
        Ni = 0.0;

    if(0. > Ns)
        Ns = 0.0;

    if(0. > Ng)
        Ng = 0.0;

    if(0. > Nh)
        Nh = 0.0;

    codi::RealReverse dep_rate_ice = 0.0;
    codi::RealReverse dep_rate_snow = 0.0;
    // Change to dimensional variables
    codi::RealReverse p_prime = ref.pref * p;
    codi::RealReverse T_prime = ref.Tref * T;
    codi::RealReverse w_prime = ref.wref * w;
    codi::RealReverse qc_prime = ref.qref * qc;
    codi::RealReverse qr_prime = ref.qref * qr;
    codi::RealReverse qv_prime = ref.qref * qv;
    codi::RealReverse qi_prime = ref.qref * qi;
    codi::RealReverse qs_prime = ref.qref * qs;
    codi::RealReverse qg_prime = ref.qref * qg;
    codi::RealReverse qh_prime = ref.qref * qh;

    codi::RealReverse T_c = T_prime - tmelt;
    codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime);
    codi::RealReverse p_sat_ice = saturation_pressure_ice(T_prime);
    codi::RealReverse ssi = qv_prime * Rv * T_prime / p_sat_ice;
    codi::RealReverse D_vtp = diffusivity(T_prime, p_prime);
    codi::RealReverse e_d = qv_prime * Rv * T_prime; // Could use R_v as well. The difference is minor
    codi::RealReverse s_sw = e_d / p_sat - 1.0; // super saturation over water
    codi::RealReverse s_si = e_d / p_sat_ice - 1.0; // super saturation over ice
    const double EPSILON = 1.0e-20;

    // Nucleation types
    // 0: force constant cloud drop number
    // <6: ccn_activation_hdcp2 (Hande et al)
    // else: ccn_activation_sk (Segal & Khain)
    const int nuc_typ = 1;
    ////////////// ccn_activation_hdcp2
    // Hande et al 2015
    // non maritime case
    std::vector<double> a_ccn = {183230691.161, 0.10147358938,
                                       -0.2922395814, 229189886.226};
    std::vector<double> b_ccn = {0.0001984051994, 4.473190485e-05,
                                       0.0001843225275, 0.0001986158191};
    std::vector<double> c_ccn = {16.2420263911, 3.22011836758,
                                       13.8499423719, 16.2461600644};
    std::vector<double> d_ccn = {287736034.13, 0.6258809883,
                                       0.8907491812, 360848977.55};

    if(qc_prime > EPSILON && w_prime > 0.0)
    {
        codi::RealReverse acoeff = a_ccn[0] * atan(b_ccn[0] * p_prime - c_ccn[0]) + d_ccn[0];
        codi::RealReverse bcoeff = a_ccn[1] * atan(b_ccn[1] * p_prime - c_ccn[1]) + d_ccn[1];
        codi::RealReverse ccoeff = a_ccn[2] * atan(b_ccn[2] * p_prime - c_ccn[2]) + d_ccn[2];
        codi::RealReverse dcoeff = a_ccn[3] * atan(b_ccn[3] * p_prime - c_ccn[3]) + d_ccn[3];
        // concentration of ccn
        codi::RealReverse nuc_n = acoeff * atan(bcoeff * log(w_prime) + ccoeff) + dcoeff;
        // std::cout << "acoeff " << acoeff << ", bcoeff " << bcoeff
        //           << ", ccoeff " << ccoeff << ", dcoeff " << dcoeff
        //           << ", log(w) " << log(w_prime) << ", atan(...) "
        //           << atan(bcoeff * log(w_prime) + ccoeff)
        //           << ", --- " << acoeff * atan(bcoeff * log(w_prime) + ccoeff) <<"\n";
        // this is per m^-3. Change to per cm^-3?
        // nuc_n /= 10e6;
        // codi::RealReverse rho = compute_rhoa(p_prime, T_prime, S);
        // codi::RealReverse tst = (cc.cloud.max_x - cc.cloud.min_x)/2 + cc.cloud.min_x;
        // std::cout << "nuc_n " << nuc_n << ", Nc " << Nc // << ", Nc_avg " << qc_prime/tst
                //   << ", Nc_tmpmin " << qc_prime / cc.cloud.min_x
                //   << ", Nc_tmpmax " << qc_prime / cc.cloud.max_x
                //   << ", p " << p_prime
                //   << ", min_x " << cc.cloud.min_x << ", qv "
                //   << qv_prime << ", dqv_avg " << nuc_n * tst
                //   << ", dqvmin " << nuc_n * cc.cloud.min_x
                //   << ",dqv actually " << (nuc_n-Nc) * cc.cloud.min_x
                //   << ", rhoa " << compute_rhoa(p_prime, T_prime, S)
                //   << "\n";
                //   << ", w " << w_prime << "\n";
        // we need to substract the already "used" ccn in cloud droplets
        // the rest can create new cloud droplets
        // codi::RealReverse Nc_tmp = qv_prime / cc.cloud.min_x;
        // TODO: This looks weird. It should be the line directly under this comment
        codi::RealReverse delta_n = max(max(nuc_n, 10.0e-6) - Nc, 0.0);
        // codi::RealReverse delta_n = max(max(nuc_n, 10.0e7) , 0.0);
        // codi::RealReverse denom = (cc.cloud.max_x - cc.cloud.min_x) / 2.0 + cc.cloud.min_x;
        codi::RealReverse delta_q = min(delta_n * cc.cloud.min_x, qv_prime);
        delta_n = delta_q / cc.cloud.min_x;
        // codi::RealReverse delta_q = min(delta_n * denom, qv_prime);
        // delta_n = delta_q /denom;


        res[Nc_idx] += delta_n;
        res[qc_idx] += delta_q;
        // DEBUG STATEMENT
        // if(delta_q/ref.qref >= 0.003)
        //     std::cout << "line 452, dqc " << delta_q << ", nuc_n " << nuc_n
        //               << ", delta_n " << delta_n << ", Nc " << Nc << ", P "
        //               << p_prime << ", w " << w_prime
        //               << ", qv " << qv_prime << " vs "
        //               << delta_n * cc.cloud.min_x
        //               << ", qc_prime " << qc_prime << "\n";
        res[qv_idx] -= delta_q;

        codi::RealReverse delta_e = latent_heat_evap(T_prime) * delta_q / specific_heat_water_vapor(T_prime);
        // Evaporation
        if(delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        else
            res[lat_heat_idx] += delta_q;

        // std::cout << delta_q << ", q_prime " << qc_prime << "\n";
        // cond += delta_q
    }

    ////////////// ice_nucleation_homhet
    // Homogeneous and heterogeneous ice nucleation based on
    // "A parametrization of cirrus cloud formation: Homogenous
    // freezing of supercooled aerosols" by B. Kaercher and
    // U. Lohmann 2002 (KL02 hereafter)
    //
    // "Physically based parameterization of cirrus cloud formation
    // for use in global atmospheric models" by B. Kaercher, J. Hendricks
    // and U. Lohmann 2006 (KHL06 hereafter)
    //
    // and Phillips et al. (2008) with extensions
    //
    // implementation by Carmen Koehler and AS
    // modified for C++ and Codipack by Maicon Hieronymus

    // heterogeneous nucleation using Hande et al.
    // ice_nucleation_het_hdcp2
    bool use_hdcp2_het = true; // Depends on nucleation type
    codi::RealReverse delta_n_a = 0.0;
    bool ndiag_mask = false;
    if(use_hdcp2_het)
    {
        // TODO: This comes from the macro physical part from which I don't have any data
        codi::RealReverse n_inact = 0.0;
        if(T_prime < T_nuc && T_prime > 180.0 && ssi > 1.0 && n_inact < ni_het_max)
        {

            const double EPSILON = 1.0e-20;
            codi::RealReverse ndiag = 0.0;
            if(qc_prime > EPSILON)
            {
                T_prime = max(T_prime, 237.1501);
                if(T_prime < 261.15)
                {
                    ndiag = nim_imm * exp(-alf_imm
                        * exp(bet_imm*log(T_prime-237.15)));
                }
            } else
            {
                // Hande et al. scheme, Eq. (3) with (2) and (1)
                codi::RealReverse T_tmp = max(T_prime, 220.001);
                if(T_tmp < 253.0)
                {
                    ndiag = nin_dep * exp(-alf_dep
                        * exp(bet_dep*log(T_tmp - 220.0)));
                    ndiag = ndiag * (a_dep * atan(b_dep*(ssi-1.0)+c_dep)+d_dep);
                }
            }
            codi::RealReverse delta_n = max(ndiag - n_inact, 0.0);
            codi::RealReverse delta_q = min(delta_n * cc.ice.min_x, qv_prime);
            delta_n = delta_q/cc.ice.min_x;

            res[qi_idx] += delta_q;
            res[Ni_idx] += delta_n;
            res[qv_idx] -= delta_q;
            n_inact += delta_n;
            // dep += delta_q;
            // lwrite_n_inpot = use_prog_in && ndiag > 1.0e-12;
            // ndiag_mask = lwrite_n_inpot
            delta_n_a = delta_n;

            // latent heating and cooling
            codi::RealReverse delta_e = latent_heat_melt(T_prime) * delta_q / specific_heat_ice(T_prime);
            // Sublimation, cooling
            if(delta_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] += delta_e;
        }
    }  else // end Hande
    {
        // heterogeneous nucleation using Phillips et al.
        // ice_nucleation_het_philips

    } // end Phillips
    // use_prog_in?
    bool use_prog_in = false;
    if(use_prog_in)
    {
        // Stuff that is being done with n_inpot
        // We leave that one out since this is only important for the
        // macrophysical part
    }

    // (optional) homogenous nucleation using KHL06
    codi::RealReverse s_crit = 2.349 - T_prime * (1.0/259.0);
    const double r_0 = 0.25e-6;         // aerosol particle radius prior to freezing
    const double alpha_d = 0.5;         // deposition coefficient (KL02; Spichtinger & Gierens 2009)
    const double M_w = 18.01528e-3;     // molecular mass of water [kg/mol]
    const double M_a = 28.96e-3;        // molecular mass of air [kg/mol]
    const double ma_w = M_w/N_avo;      // mass of water molecule [kg]
    const double svol = ma_w / rho_ice; // specific volume of a water molecule in ice

    if(ssi > s_crit && T_prime < 235.0 && Ni < ni_hom_max)
    {
        codi::RealReverse x_i = particle_mean_mass(qi_prime, Ni,
            cc.ice.min_x, cc.ice.max_x);
        codi::RealReverse r_i = pow(x_i/(4.0/3.0*M_PI*rho_ice), 1.0/3.0);
        codi::RealReverse v_th = sqrt(8.0*k_b*T_prime/(M_PI*ma_w));
        codi::RealReverse flux = alpha_d * v_th/4.0;
        codi::RealReverse n_sat = p_sat_ice/(k_b*T_prime);

        // coeffs of supersaturation equation
        std::vector<codi::RealReverse> acoeff(3);
        acoeff[0] = (L_ed * grav) / (cp * Rv * T_prime*T_prime) - grav/(Ra * T_prime);
        acoeff[1] = 1.0/n_sat;
        acoeff[2] = (L_ed*L_ed * M_w * ma_w)/(cp * p_prime * T_prime * M_a);

        // coeffs of depositional growth equation
        std::vector<codi::RealReverse> bcoeff(2);
        bcoeff[0] = flux * svol * n_sat * (ssi - 1.0);
        bcoeff[1] = flux/diffusivity(T_prime, p_prime);

        // pre-existing ice crystals included as reduced updraft speed
        codi::RealReverse ri_dot = bcoeff[0] / (1.0 + bcoeff[1] * r_i);
        codi::RealReverse R_ik = (4.0*M_PI) / svol * Ni * r_i*r_i * ri_dot;
        codi::RealReverse w_pre = max(0.0, (acoeff[1] + acoeff[2] * ssi)/(acoeff[0]*ssi)*R_ik); // KHL06 Eq. 19

        // homogenous nucleation event
        if(w_prime > w_pre)
        {
            codi::RealReverse cool = grav / cp*w_prime;
            codi::RealReverse ctau = T_prime * (0.004*T_prime - 2.0) + 304.4;
            codi::RealReverse tau = 1.0/(ctau*cool);
            codi::RealReverse delta = bcoeff[1] * r_0;
            codi::RealReverse phi = acoeff[0]*ssi / (acoeff[1]+acoeff[2]*ssi) * (w_prime - w_pre);

            // monodisperse approximation following KHL06
            codi::RealReverse kappa = 2.0 * bcoeff[0]*bcoeff[1]*tau/((1.0+delta)*(1.0+delta));
            codi::RealReverse sqrtkap = sqrt(kappa);
            codi::RealReverse ren = 3.0*sqrtkap / (2.0 + sqrt(1.0+9.0*kappa/M_PI));
            codi::RealReverse R_imfc = 4.0 * M_PI * bcoeff[0]/(bcoeff[1]*bcoeff[1]) / svol;
            codi::RealReverse R_im = R_imfc / (1.0+delta) * (delta*delta - 1.0
                + (1.0+0.5*kappa*(1.0+delta)*(1.0+delta)) * ren/sqrtkap);

            // number concentration and radius of ice particles
            codi::RealReverse ni_hom = phi/R_im;
            codi::RealReverse ri_0 = 1.0 + 0.5*sqrtkap * ren;
            codi::RealReverse ri_hom = (ri_0 * (1.0+delta) - 1.0) / bcoeff[1];
            codi::RealReverse mi_hom = (4.0/3.0 * M_PI * rho_ice) * ni_hom * ri_hom*ri_hom*ri_hom;
            mi_hom = max(mi_hom, cc.ice.min_x);

            codi::RealReverse delta_n = max(min(ni_hom, ni_hom_max), 0.0);
            codi::RealReverse delta_q = min(delta_n * mi_hom, qv_prime);

            res[Ni_idx] += delta_n;
            res[qi_idx] += delta_q;
            res[qv_idx] -= delta_q;
            // depo += nuc_q;

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * delta_q / specific_heat_ice(T_prime);
            // Sublimation, cooling
            if(delta_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] += delta_e;
        }
    }

    ////////////// cloud_freeze
    // Homogeneous freezing of cloud droplets
    // Freeze only if cloud droplets to freeze are available
    // and the temperature is low enough
    if(qc_prime > 0.0 && T_c < -30.0)
    {
        codi::RealReverse x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x, cc.cloud.max_x);

        codi::RealReverse delta_qi;
        codi::RealReverse delta_ni;
        codi::RealReverse j_hom;
        // instantaneous freezing for temperatures below -50 °C
        if(T_c < -50.0)
        {
            delta_qi = qc;
            delta_ni = Nc;
        } else
        {
            if(T_c > -30.0)
                j_hom = 1.0e6 / rho_w * pow(10,
                    -7.63-2.996*(T_c+30.0));
            else
                j_hom = 1.0e6 / rho_w * pow(10,
                    - 243.4
                    - 14.75 * T_c
                    - 0.307 * T_c * T_c
                    - 0.00287 * T_c * T_c * T_c
                    - 0.0000102 * pow(T_c, 4));
            delta_ni = j_hom * qc;
            delta_qi = j_hom * qc * x_c * cc.cloud.c_z;
            delta_ni = min(delta_ni, Nc);
            delta_qi = min(delta_qi, qc_prime);
        }
        res[qi_idx] += delta_qi;
        res[Ni_idx] += delta_ni;
        // Remove cloud droplets
        res[qc_idx] -= delta_qi;
        // DEBUG STATEMENT
        if(delta_qi/ref.qref >= 0.003)
            std::cout << "line 666, dqc " << delta_qi << "\n";
        res[Nc_idx] -= delta_ni;

        codi::RealReverse delta_e = latent_heat_melt(T_prime) * delta_qi / specific_heat_ice(T_prime);
        // Melting, cooling
        if(delta_qi < 0.0)
            res[lat_cool_idx] += delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] += delta_e;

    }
    ////////////// vapor_dep_relaxation
    // Depositional growth of all ice particles
    // Deposition and sublimation, where deposition rate of ice
    // and snow are being stored
    codi::RealReverse dep_ice = 0.0;
    codi::RealReverse dep_snow = 0.0;
    codi::RealReverse dep_graupel = 0.0;
    codi::RealReverse dep_hail = 0.0;
    codi::RealReverse g_i = 0.0;
    codi::RealReverse s_ice = 0.0; // supersaturation over ice

    auto vapor_deposition = [&](
        codi::RealReverse &q,
        codi::RealReverse &N,
        particle_model_constants_t &pc,
        codi::RealReverse &dep)
    {
        if(q == 0.0)
        {
            dep = 0.0;
        } else
        {
            codi::RealReverse x = particle_mean_mass(q, N, pc.min_x, pc.max_x);
            codi::RealReverse d = particle_diameter(x, pc.a_geo, pc.b_geo);
            codi::RealReverse v = particle_velocity(x, pc.a_vel, pc.b_vel) * pc.rho_v;
            codi::RealReverse f_v = pc.a_f + pc.b_f * sqrt(d*v);
            f_v = max(f_v, pc.a_f/pc.a_ven);
            dep = g_i * N * pc.c_s * d * f_v * s_ice;
        }
    };


    if(T_prime < tmelt)
    {
        codi::RealReverse e_d = qv_prime * Rv * T_prime;
        s_ice = e_d/p_sat_ice - 1.0; // supersaturation over ice
        g_i = 4.0*M_PI / (L_ed*L_ed / (K_T*Rv*T_prime*T_prime) + Rv*T_prime / (D_vtp*p_sat_ice));
    }

    vapor_deposition(qi_prime, Ni, cc.ice, dep_ice);
    vapor_deposition(qs_prime, Ns, cc.snow, dep_snow);
    vapor_deposition(qg_prime, Ng, cc.graupel, dep_graupel);
    vapor_deposition(qh_prime, Nh, cc.hail, dep_hail);

    if(T_prime < tmelt)
    {
        // Depositional growth based on
        // "A New Double-Moment Microphysics Parameterization for Application in Cloud and
        // Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov
        codi::RealReverse qvsidiff = qv_prime - p_sat_ice /(Rv*T_prime);
        if(abs(qvsidiff) > EPSILON)
        {
            codi::RealReverse tau_i_i = qvsidiff*dep_ice;
            codi::RealReverse tau_s_i = qvsidiff*dep_snow;
            codi::RealReverse tau_g_i = qvsidiff*dep_graupel;
            codi::RealReverse tau_h_i = qvsidiff*dep_hail;

            codi::RealReverse xi_i = tau_i_i + tau_s_i + tau_g_i + tau_h_i;
            // TODO: Check wether dt is needed here or not
            codi::RealReverse xfac = (xi_i < EPSILON) ? (codi::RealReverse) 0.0 : qvsidiff/xi_i * (1.0-exp(-xi_i));

            dep_ice     = xfac * tau_i_i;
            dep_snow    = xfac * tau_s_i;
            dep_graupel = xfac * tau_g_i;
            dep_hail    = xfac * tau_h_i;

            // Is that even necessary?
            if(qvsidiff < 0.0)
            {
                dep_ice     = max(dep_ice,      -qi_prime);
                dep_snow    = max(dep_snow,     -qs_prime);
                dep_graupel = max(dep_graupel,  -qg_prime);
                dep_hail    = max(dep_hail,     -qh_prime);
            }

            codi::RealReverse dep_sum = dep_ice + dep_graupel + dep_snow + dep_hail;

            res[qi_idx] += dep_ice;
            res[qs_idx] += dep_snow;
            res[qg_idx] += dep_graupel;
            res[qh_idx] += dep_hail;
            res[qv_idx] -= dep_sum;
            // If I want to store the deposition/sublimination
            // if(dep_ice > 0)
            // {
            //     depo += dep_ice;
            // } else
            // {
            //     subl += dep_ice;
            // }
            // if(dep_snow > 0)
            // {
            //     depo += dep_snow;
            // } else
            // {
            //     subl += dep_snow;
            // }
            // if(dep_graupel > 0)
            // {
            //     depo += dep_graupel;
            // } else
            // {
            //     subl += dep_graupel;
            // }
            // if(dep_hail > 0)
            // {
            //     depo += dep_hail;
            // } else
            // {
            //     subl += dep_hail;
            // }
            dep_rate_ice += dep_ice;
            dep_rate_snow += dep_snow;

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * dep_sum / specific_heat_ice(T_prime);
            // Sublimation, cooling
            if(dep_sum > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] -= delta_e;
        }
    }

    ////////////// ice-ice collisions
    codi::RealReverse x_i = particle_mean_mass(qi_prime, Ni, cc.ice.min_x, cc.ice.max_x);
    codi::RealReverse D_i = particle_diameter(x_i, cc.ice.a_geo, cc.ice.b_geo);

    //// ice self collection
    if(Ni > 0.0 && qi_prime > q_crit_i && D_i > D_crit_i)
    {
        codi::RealReverse x_conv_i = pow(D_conv_i/cc.snow.a_geo, 1.0/cc.snow.b_geo);
        // efficiency depends on temperature here (Cotton et al 1986)
        // also Straka 1989, page 53
        codi::RealReverse e_coll = min(pow(10, 0.035*T_c-0.7), 0.2);
        codi::RealReverse vel_i = pow(cc.ice.a_vel * x_i, cc.ice.b_vel * cc.ice.rho_v);

        codi::RealReverse delta_n = M_PI/4.0 * e_coll * cc.ice.sc_delta_n
            * Ni * Ni * D_i * D_i * sqrt(
                cc.ice.sc_theta_n * vel_i * vel_i + 2.0 * cc.ice.s_vel * cc.ice.s_vel);
        codi::RealReverse delta_q = M_PI/4.0 * e_coll * cc.ice.sc_delta_q
            * Ni * qi_prime * D_i * D_i * sqrt(
                cc.ice.sc_theta_q * vel_i * vel_i + 2.0 * cc.ice.s_vel * cc.ice.s_vel);

        delta_q = min(delta_q, qi_prime);
        delta_n = min(min( delta_n, delta_q/x_conv_i), Ni);

        res[qi_idx] -= delta_q;
        res[qs_idx] += delta_q;
        res[Ni_idx] -= delta_n;
        res[Ns_idx] += delta_n/2.0;
    }

    //// snow self collection
    if(qs_prime > q_crit)
    {
        // temperature dependent sticking efficiency Lin (1983)
        codi::RealReverse e_coll = max(0.1, min(exp(0.09*(T_prime-tmelt)), 1.0));
        codi::RealReverse x_s = particle_mean_mass(qs_prime, Ns, cc.snow.min_x, cc.snow.max_x);
        codi::RealReverse D_s = particle_diameter(x_s, cc.snow.a_geo, cc.snow.b_geo);
        codi::RealReverse vel_s = particle_velocity(x_s, cc.snow.a_vel, cc.snow.b_vel) * cc.snow.rho_v;

        res[Ns_idx] -= M_PI/8.0 * e_coll * Ns * Ns * cc.snow.sc_delta_n * D_s * D_s
            * sqrt(cc.snow.sc_theta_n * vel_s * vel_s + 2.0 * cc.snow.s_vel * cc.snow.s_vel);
    }

    // particle particle collection

    // graupel+ice  -> graupel
    // graupel+snow -> graupel
    // hail+ice  -> hail
    // hail+snow -> hail
    // snow+ice  -> snow
    auto particle_collection = [&] (
        codi::RealReverse &q1,
        codi::RealReverse &q2,
        codi::RealReverse &N1,
        codi::RealReverse &N2,
        collection_model_constants_t &coeffs,
        particle_model_constants_t &pc1,
        particle_model_constants_t &pc2)
    {

        codi::RealReverse e_coll = min(exp(0.09*T_c), 1.0);
        codi::RealReverse x_1 = particle_mean_mass(q1, N1, pc1.min_x, pc1.max_x);
        codi::RealReverse d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);
        codi::RealReverse v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;

        codi::RealReverse x_2 = particle_mean_mass(q2, N2, pc1.min_x, pc1.max_x);
        codi::RealReverse d_2 = particle_diameter(x_2, pc1.a_geo, pc1.b_geo);
        codi::RealReverse v_2 = particle_velocity(x_2, pc1.a_vel, pc1.b_vel) * pc1.rho_v;

        codi::RealReverse coll_n = M_PI/4 * N1 * N2 * e_coll
            * (coeffs.delta_n_aa * d_2 * d_2
                + coeffs.delta_n_ab * d_2 * d_1
                + coeffs.delta_n_bb * d_1 * d_1)
            * sqrt(coeffs.theta_n_aa * v_2 * v_2
                - coeffs.theta_n_ab * v_2 * v_1
                + coeffs.theta_n_bb * v_1 * v_1
                + pc1.s_vel * pc1.s_vel);
        codi::RealReverse coll_q = M_PI/4 * N1 * N2 * e_coll
            * (coeffs.delta_q_aa * d_2 * d_2
                + coeffs.delta_q_ab * d_2 * d_1
                + coeffs.delta_q_bb * d_1 * d_1)
            * sqrt(coeffs.theta_q_aa * v_2 * v_2
                - coeffs.theta_q_ab * v_2 * v_1
                + coeffs.theta_q_bb * v_1 * v_1
                + pc1.s_vel * pc1.s_vel);
        coll_n = min(N1, coll_n);
        coll_q = min(q1, coll_q);
        std::vector<codi::RealReverse> r(2);
        r[p_idx] = coll_n;
        r[T_idx] = coll_q;
        return r;
    };

    //// particle_collection snow
    if(qi_prime > q_crit && qs_prime > q_crit)
    {
        std::vector<codi::RealReverse> delta = particle_collection(
            qi_prime, qs_prime, Ni, Ns, cc.coeffs_sic, cc.ice, cc.snow);
        res[qs_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
    }

    //// graupel self collection
    if(qg_prime > q_crit)
    {
        codi::RealReverse x_g = particle_mean_mass(qg_prime, Ng,
            cc.graupel.min_x, cc.graupel.max_x);
        codi::RealReverse d_g = particle_diameter(x_g,
            cc.graupel.a_geo, cc.graupel.b_geo);
        codi::RealReverse v_g = particle_velocity(x_g,
            cc.graupel.a_vel, cc.graupel.b_vel);
        codi::RealReverse delta_n = cc.graupel.sc_coll_n
            * Ng * Ng * d_g * d_g * v_g;
        // sticking efficiency does only distinguish dry and wet
        delta_n *= (T_prime > tmelt) ? ecoll_gg_wet : ecoll_gg;
        delta_n = min(delta_n, Ng);
        res[Ng_idx] -= delta_n;
    }

    // particle particle collection
    // ice and graupel collision
    if(qi_prime > q_crit && qg_prime > q_crit)
    {
        std::vector<codi::RealReverse> delta = particle_collection(
            qi_prime, qg_prime, Ni, Ng, cc.coeffs_gic, cc.ice, cc.graupel);
        res[qg_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
    }

    // particle particle collection
    // snow and graupel collision
    if(qs_prime > q_crit && qg_prime > q_crit)
    {
        std::vector<codi::RealReverse> delta = particle_collection(
            qs_prime, qg_prime, Ns, Ng, cc.coeffs_gsc, cc.snow, cc.graupel);
        res[qg_idx] += delta[1];
        res[qs_idx] -= delta[1];
        res[Ns_idx] -= delta[0];
    }

    ////////////// graupel_hail_conv_wet_gamlook
    // conversion graupel to hail and hail collisions
    codi::RealReverse x_g = particle_mean_mass(qg_prime, Ng,
            cc.graupel.min_x, cc.graupel.max_x);
    codi::RealReverse d_g = particle_diameter(x_g,
        cc.graupel.a_geo, cc.graupel.b_geo);
    codi::RealReverse Ng_tmp = qg_prime/x_g;
    // supercooled liquid water = rain + cloud water
    codi::RealReverse qwa_prime = qr_prime + qc_prime;

    if(qwa_prime > 1e-3 && T_c < 0 && qg_prime > q_crit_c)
    {
        codi::RealReverse qis = qi + qs;

        codi::RealReverse d_sep = wet_growth_diam(p_prime, T_prime, qwa_prime,
            qi_prime, ltabdminwgg);
        if(d_sep > 0.0 && d_sep < 10.0*d_g)
        {
            codi::RealReverse xmin = pow(d_sep/cc.graupel.a_geo, 1.0/cc.graupel.b_geo);
            codi::RealReverse lam = pow(cc.graupel.g2/(cc.graupel.g1*x_g), cc.graupel.mu);
            xmin = pow(xmin, cc.graupel.mu);
            codi::RealReverse n_0 = cc.graupel.mu * Ng * pow(cc.graupel.nm1, cc.graupel.g1);

            codi::RealReverse conv_n = n_0;
            codi::RealReverse conv_q = n_0;

            conv_n = min(conv_n, Ng);
            conv_q = min(conv_q, qg);
            // Graupel q, n
            res[qg_idx] += qg_prime - conv_q;
            res[Ng_idx] += Ng - conv_n;
            // Hail q, n
            res[qh_idx] += qh_prime + conv_q;
            res[Nh_idx] += Nh + conv_n;
        }
    }

    ////////////// hail collisions
    // ice and hail
    if(qi_prime > q_crit && qh_prime > q_crit)
    {
        std::vector<codi::RealReverse> delta = particle_collection(
            qi_prime, qh_prime, Ni, Nh, cc.coeffs_hic, cc.ice, cc.hail);
        res[qh_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
    }

    // snow and hail collision
    if(qs_prime > q_crit && qg_prime > q_crit)
    {
        std::vector<codi::RealReverse> delta = particle_collection(
            qs_prime, qh_prime, Ns, Nh, cc.coeffs_hsc, cc.snow, cc.hail);
        res[qh_idx] += delta[1];
        res[qs_idx] -= delta[1];
        res[Ns_idx] -= delta[0];
    }

    ////////////// Riming of ice with cloud and rain droplets and conversion to graupel
    codi::RealReverse rime_rate_qc, rime_rate_qr, rime_rate_qi, rime_rate_qs;
    codi::RealReverse rime_rate_nc, rime_rate_nr;
    // riming cloud core
    // rate of ice or snow collecting cloud droplets where values are hard
    // coded for *rain* droplets according to COSMO comments
    // TODO: Check if this is actually for cloud or rain droplets
    auto riming_cloud_core = [&](
        codi::RealReverse &q1,
        codi::RealReverse &N1,
        particle_model_constants_t &pc1,
        collection_model_constants_t &coeffs,
        codi::RealReverse &rime_rate_qb,
        codi::RealReverse &rime_rate_nb)
    {
        codi::RealReverse x_1 = particle_mean_mass(q1, N1, pc1.min_x, pc1.max_x);
        codi::RealReverse d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);
        codi::RealReverse x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x, cc.cloud.max_x);
        codi::RealReverse d_c = particle_diameter(x_c, cc.cloud.a_geo, cc.cloud.b_geo);

        codi::RealReverse const1 = const0 * pc1.ecoll_c;

        if(qc_prime > q_crit_c && q1 > pc1.q_crit_c
            && d_c > D_crit_c && d_1 > pc1.d_crit_c)
        {
            codi::RealReverse v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;
            codi::RealReverse v_c = particle_velocity(x_c, cc.cloud.a_vel, cc.cloud.b_vel) * cc.cloud.rho_v;
            codi::RealReverse tmp = const1*(d_c - D_crit_c);
            codi::RealReverse e_coll = min(pc1.ecoll_c, max(tmp, ecoll_min));

            // times dt ?
            rime_rate_qb = M_PI/4.0 * e_coll * N1 * qc_prime
                * (coeffs.delta_q_aa * d_1*d_1
                    + coeffs.delta_q_ab * d_1*d_c
                    + coeffs.delta_q_bb * d_c*d_c)
                * sqrt(coeffs.theta_q_aa * v_1*v_1
                    - coeffs.theta_q_ab * v_1*v_c
                    + coeffs.theta_q_bb * v_c*v_c
                    + pc1.s_vel*pc1.s_vel);

            rime_rate_nb = M_PI/4.0 * e_coll * N1 * Nc
                * (coeffs.delta_n_aa * d_1*d_1
                    + coeffs.delta_n_ab * d_1*d_c
                    + coeffs.delta_n_bb * d_c*d_c)
                * sqrt(coeffs.theta_n_aa * v_1*v_1
                    - coeffs.theta_n_ab * v_1*v_c
                    + coeffs.theta_n_bb * v_c*v_c
                    + pc1.s_vel*pc1.s_vel);
        } else
        {
            rime_rate_qb = 0.0;
            rime_rate_nb = 0.0;
        }
    };

    auto riming_rain_core = [&](
        codi::RealReverse &q1,
        codi::RealReverse &N1,
        particle_model_constants_t &pc1,
        collection_model_constants_t &coeffs,
        codi::RealReverse &rime_rate_qa,
        codi::RealReverse &rime_rate_qb,
        codi::RealReverse &rime_rate_nb)
    {
        codi::RealReverse x_1 = particle_mean_mass(q1, N1, pc1.min_x, pc1.max_x);
        codi::RealReverse d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);

        if(qr_prime > q_crit && q1 > q_crit_r && d_1 > D_crit_r)
        {
            codi::RealReverse x_r = particle_mean_mass(qr_prime, Nc, cc.rain.min_x, cc.rain.max_x);
            codi::RealReverse d_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);

            codi::RealReverse v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;
            codi::RealReverse v_r = particle_velocity(x_r, cc.rain.a_vel, cc.rain.b_vel) * cc.rain.rho_v;

            // times dt ?
            rime_rate_qb = M_PI/4.0 * N1 * qr_prime
                * (coeffs.delta_n_aa * d_1*d_1
                    + coeffs.delta_q_ab * d_1*d_r
                    + coeffs.delta_q_bb * d_r*d_r)
                * sqrt(coeffs.theta_n_aa * v_1*v_1
                    - coeffs.theta_q_ab * v_1*v_r
                    + coeffs.theta_q_bb * v_r*v_r
                    + pc1.s_vel*pc1.s_vel);

            rime_rate_qa = M_PI/4.0 * Nr * q1
                * (coeffs.delta_q_aa * d_1*d_1
                    + coeffs.delta_q_ba * d_1*d_r
                    + coeffs.delta_n_bb * d_r*d_r)
                * sqrt(coeffs.theta_q_aa * v_1*v_1
                    - coeffs.theta_q_ba * v_1*v_r
                    + coeffs.theta_n_bb * v_r*v_r
                    + pc1.s_vel*pc1.s_vel);

            rime_rate_nb = M_PI/4.0 * N1 * Nc
                * (coeffs.delta_n_aa * d_1*d_1
                    + coeffs.delta_n_ab * d_1*d_r
                    + coeffs.delta_n_bb * d_r*d_r)
                * sqrt(coeffs.theta_n_aa * v_1*v_1
                    - coeffs.theta_n_ab * v_1*v_r
                    + coeffs.theta_n_bb * v_r*v_r
                    + pc1.s_vel*pc1.s_vel);
        } else
        {
            rime_rate_qa = 0.0;
            rime_rate_qb = 0.0;
            rime_rate_nb = 0.0;
        }
    };

    riming_cloud_core(qi_prime, Ni, cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc);
    riming_rain_core(qi_prime, Ni, cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr);
    // ice riming
    if(dep_rate_ice > 0.0 && dep_rate_ice >= rime_rate_qc+rime_rate_qr)
    {
        // Depositional growth is stronger than riming growth, therefore ice stays ice
        // ice cloud riming
        if(rime_rate_qc > 0.0)
        {
            codi::RealReverse rime_q = min(qc_prime, rime_rate_qc);
            codi::RealReverse rime_n = min(Nc, rime_rate_nc);
            // Ice
            res[qi_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;
            // DEBUG STATEMENT
            if(rime_q/ref.qref >= 0.003)
                std::cout << "line 1128, dqc " << rime_q << "\n";
            // Cloud N
            res[Nc_idx] -= rime_n;

            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                // Ice N
                res[Ni_idx] += C_mult * mult_1 * mult_2 * rime_q;
            }
        }
        // ice rain riming
        if(rime_rate_qr > 0.0)
        {
            codi::RealReverse rime_q = min(rime_rate_qr, qr_prime);
            codi::RealReverse rime_n = min(Nr, rime_rate_nr);
            // Snow
            res[qs_idx] += rime_q;
            // Rain
            res[qr_idx] -= rime_q;
            // Rain N
            res[Nr_idx] -= rime_n;
            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                // Ice N
                res[Ni_idx] += C_mult * mult_1 * mult_2 * rime_q;
            }

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;
        }
    } else
    {
        // Depositional growth negative or smaller than riming growth,
        // therefore ice is allowed to convert to graupel and / or hail
        // ice cloud riming
        if(rime_rate_qc > 0.0)
        {
            codi::RealReverse x_i = particle_mean_mass(qi_prime, Ni,
                cc.ice.min_x, cc.ice.max_x);
            codi::RealReverse d_i = particle_diameter(x_i,
                cc.ice.a_geo, cc.ice.b_geo);
            codi::RealReverse rime_q = min(rime_rate_qc, qc_prime);
            codi::RealReverse rime_n = min(rime_rate_nc, Nc);
            // Ice
            res[qi_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;
            // DEBUG STATEMENT
            if(rime_q/ref.qref >= 0.003)
                std::cout << "line 1190, dqc " << rime_q << "\n";
            // Cloud N
            res[Nc_idx] -= rime_n;

            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                // Ice N
                res[Ni_idx] += C_mult * mult_1 * mult_2 * rime_q;
            }

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            // Conversion ice -> graupel
            // Technically I had to recalculate x_i given the new qi, Ni
            if(d_i > D_conv_ig)
            {
                codi::RealReverse conv_q = rime_q
                    / (const5*(M_PI/6.0 * rho_ice * d_i*d_i*d_i/x_i -1.0));
                conv_q = min(qi_prime, conv_q);
                codi::RealReverse qi_tmp = qi_prime+dt*res[qi_idx]/ref.qref;
                x_i = particle_mean_mass(qi_tmp, Ni,
                    cc.ice.min_x, cc.ice.max_x);
                codi::RealReverse tmp = conv_q / max(x_i, x_conv);
                codi::RealReverse conv_n = min(tmp, Ni);

                // Ice
                res[qi_idx] -= conv_q;
                // Graupel
                res[qg_idx] += conv_q;
                // Ice N
                res[Ni_idx] -= conv_n;
                // Graupel N
                res[Ng_idx] += conv_n;
            }
        }

        // ice rain riming
        if(rime_rate_qi > 0.0)
        {
            codi::RealReverse rime_qi = min(rime_rate_qi, qi_prime);
            codi::RealReverse rime_qr = min(rime_rate_qr, qr_prime);
            codi::RealReverse rime_n = min(min(rime_rate_nr, Nr), Ni);

            // Ice
            res[qi_idx] -= rime_qi;
            // Rain
            res[qr_idx] -= rime_qr;
            // Ice N
            res[Ni_idx] -= rime_n;
            // Rain N
            res[Nr_idx] -= rime_n;

            codi::RealReverse mult_q = 0.0;
            codi::RealReverse mult_n = 0.0;
            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                mult_n = C_mult * mult_1 * mult_2 * rime_qr;
                codi::RealReverse tmp = mult_n*cc.ice.min_x;
                mult_q = min(rime_qr, tmp);
            }

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_qi / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_qi > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] -= delta_e;

            if(T_prime >= tmelt)
            {
                codi::RealReverse qr_tmp = qr_prime+dt*res[qr_idx]/ref.qref;
                codi::RealReverse Nr_tmp = Nr*res[Nr_idx]*dt;
                codi::RealReverse x_r = particle_mean_mass(
                    qr_tmp, Nr_tmp,
                    cc.rain.min_x, cc.rain.max_x);
                // Ice
                res[qi_idx] += rime_qi;
                // Rain
                res[qr_idx] += rime_qr;
                // Ice N
                res[Ni_idx] += rime_n;
                // Rain N
                res[Nr_idx] += rime_qr/x_r;

                delta_e = latent_heat_melt(T_prime) * rime_qi / specific_heat_ice(T_prime);
                // Melting, cooling
                if(rime_qi > 0.0)
                    res[lat_cool_idx] -= delta_e;
                // Freezing, heating
                else
                    res[lat_heat_idx] -= delta_e;
            } else
            {
                // from multiplication
                // Ice
                res[qi_idx] += mult_q;
                // Ice N
                res[Ni_idx] += mult_n;
                // riming to graupel
                // Graupel
                res[qg_idx] += rime_qi + rime_qr - mult_q;
                // Graupel N
                res[Ng_idx] += rime_n;
            }
        }
    }

    // snow riming
    riming_cloud_core(qs_prime, Ns, cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc);
    riming_rain_core(qs_prime, Ns, cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr);
    if(dep_rate_snow > 0.0 && dep_rate_snow >= rime_rate_qs+rime_rate_qr)
    {

        // Depositional growth is stronger than riming growth, therefore ice stays ice
        // ice cloud riming
        if(rime_rate_qc > 0.0)
        {
            codi::RealReverse rime_q = min(qc_prime, rime_rate_qc);
            codi::RealReverse rime_n = min(Nc, rime_rate_nc);
            // Snow
            res[qs_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;
            // DEBUG STATEMENT
            if(rime_q/ref.qref >= 0.003)
                std::cout << "line 1330, dqc " << rime_q << "\n";
            // Cloud N
            res[Nc_idx] -= rime_n;

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                codi::RealReverse mult_n = C_mult * mult_1 * mult_2 * rime_q;
                codi::RealReverse mult_q = mult_n * cc.ice.min_x;
                mult_q = min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
            }
        }
        // snow rain riming
        if(rime_rate_qr > 0.0)
        {
            codi::RealReverse rime_q = min(rime_rate_qr, qr_prime);
            codi::RealReverse rime_n = min(Nr, rime_rate_nr);
            // Snow
            res[qs_idx] += rime_q;
            // Rain
            res[qr_idx] -= rime_q;
            // Rain N
            res[Nr_idx] -= rime_n;

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                codi::RealReverse mult_n = C_mult * mult_1 * mult_2 * rime_q;
                codi::RealReverse mult_q = mult_n * cc.ice.min_x;
                mult_q = min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
            }
        }
    } else
    {
        // Depositional growth negative or smaller than riming growth,
        // therefore snow is allowed to convert to graupel and / or hail
        // snow cloud riming
        if(rime_rate_qc > 0.0)
        {
            codi::RealReverse x_s = particle_mean_mass(qs_prime, Ns,
                cc.snow.min_x, cc.snow.max_x);
            codi::RealReverse d_s = particle_diameter(x_s,
                cc.snow.a_geo, cc.snow.b_geo);
            codi::RealReverse rime_q = min(rime_rate_qc, qc_prime);
            codi::RealReverse rime_n = min(rime_rate_nc, Nc);
            // Snow
            res[qs_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;
            // DEBUG STATEMENT
            if(rime_q/ref.qref >= 0.003)
                std::cout << "line 1417, dqc " << rime_q << "\n";
            // Cloud N
            res[Nc_idx] -= rime_n;

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            codi::RealReverse mult_q = 0.0;
            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                codi::RealReverse mult_n = C_mult * mult_1 * mult_2 * rime_q;
                mult_q = mult_n * cc.ice.min_x;
                mult_q = min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
            }

            // Conversion snow -> graupel
            if(d_s > D_conv_sg)
            {
                codi::RealReverse conv_q = (rime_q - mult_q)
                    / (const5*(M_PI/6.0 * rho_ice * d_s*d_s*d_s/x_s -1.0));
                conv_q = min(qs_prime, conv_q);
                codi::RealReverse qs_tmp = qs_prime+dt*res[qs_idx]/ref.qref;
                x_s = particle_mean_mass(qs_tmp, Ns,
                    cc.snow.min_x, cc.snow.max_x);
                codi::RealReverse tmp = conv_q / max(x_s, x_conv);
                codi::RealReverse conv_n = min(tmp, Ns);

                // Snow
                res[qs_idx] -= conv_q;
                // Graupel
                res[qg_idx] += conv_q;
                // Snow N
                res[Ns_idx] -= conv_n;
                // Graupel N
                res[Ng_idx] += conv_n;
            }
        }

        // Snow rain riming
        if(rime_rate_qs > 0.0)
        {
            codi::RealReverse rime_qs = min(rime_rate_qs, qs_prime);
            codi::RealReverse rime_qr = min(rime_rate_qr, qr_prime);
            codi::RealReverse rime_n = min(min(rime_rate_nr, Nr), Ns);

            // Snow
            res[qs_idx] -= rime_qs;
            // Rain
            res[qr_idx] -= rime_qr;
            // Snow N
            res[Ns_idx] -= rime_n;
            // Rain N
            res[Nr_idx] -= rime_n;

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_qr / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_qr < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            codi::RealReverse mult_q = 0.0;
            codi::RealReverse mult_n = 0.0;
            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                mult_n = C_mult * mult_1 * mult_2 * rime_qr;
                codi::RealReverse tmp = mult_n*cc.ice.min_x;
                mult_q = min(rime_qr, tmp);
            }
            if(T_prime >= tmelt)
            {
                codi::RealReverse qr_tmp = qr_prime+dt*res[qr_idx]/ref.qref;
                codi::RealReverse Nr_tmp = Nr*res[Nr_idx]*dt;
                codi::RealReverse x_r = particle_mean_mass(
                    qr_tmp, Nr_tmp,
                    cc.rain.min_x, cc.rain.max_x);

                // Snow
                res[qs_idx] += rime_qs;
                // Rain
                res[qr_idx] += rime_qr;
                // Snow N
                res[Ns_idx] += rime_n;
                // Rain N
                res[Nr_idx] += rime_qr/x_r;

                codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_qr / specific_heat_ice(T_prime);
                // Melting, cooling
                if(rime_qr > 0.0)
                    res[lat_cool_idx] -= delta_e;
                // Freezing, heating
                else
                    res[lat_heat_idx] -= delta_e;
            } else
            {
                // from multiplication
                // Ice
                res[qi_idx] += mult_q;
                // Ice N
                res[Ni_idx] += mult_n;
                // riming to graupel
                // Graupel
                res[qg_idx] += rime_qs + rime_qr - mult_q;
                // Graupel N
                res[Ng_idx] += rime_n;
            }
        }
    }

    auto particle_cloud_riming = [&](
        codi::RealReverse &q1,
        codi::RealReverse &N1,
        codi::RealReverse &resq,
        codi::RealReverse &resn,
        collection_model_constants_t &coeffs,
        particle_model_constants_t &pc1)
    {
        codi::RealReverse const1 = const0 * pc1.sc_coll_n;
        codi::RealReverse x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x, cc.cloud.max_x);
        codi::RealReverse d_c = particle_diameter(x_c, cc.cloud.a_geo, cc.cloud.b_geo);
        codi::RealReverse x_1 = particle_mean_mass(q1, N1, pc1.min_x, pc1.max_x);
        codi::RealReverse d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);

        if(qc_prime > q_crit_c && q1 > pc1.q_crit_c
            && d_1 > pc1.d_crit_c && d_c > D_crit_c)
        {
            codi::RealReverse v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;
            codi::RealReverse v_c = particle_velocity(x_c, cc.cloud.a_vel, cc.cloud.b_vel) * cc.cloud.rho_v;

            codi::RealReverse e_coll_n = min(pc1.ecoll_c, max(const1*(d_c-D_crit_c), ecoll_min));
            codi::RealReverse e_coll_q = e_coll_n;

            codi::RealReverse rime_n = M_PI/4.0 * e_coll_n * N1 * Nc
                * (coeffs.delta_n_aa * d_1 * d_1
                    + coeffs.delta_n_ab * d_1 * d_c
                    + coeffs.delta_n_bb * d_c * d_c)
                * sqrt(coeffs.theta_n_aa * v_1 * v_1
                    - coeffs.theta_n_ab * v_1 * v_c
                    + coeffs.theta_n_bb * v_c * v_c);
            codi::RealReverse rime_q = M_PI/4.0 * e_coll_q * N1 * qc_prime
                * (coeffs.delta_q_aa * d_1 * d_1
                    + coeffs.delta_q_ab * d_1 * d_c
                    + coeffs.delta_q_bb * d_c * d_c)
                * sqrt(coeffs.theta_q_aa * v_1 * v_1
                    - coeffs.theta_q_ab * v_1 * v_c
                    + coeffs.theta_q_bb * v_c * v_c);
            rime_q = min(qc_prime, rime_q);
            rime_n = min(Nc, rime_n);

            resq += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;
            // DEBUG STATEMENT
            if(rime_q/ref.qref >= 0.003)
                std::cout << "line 1592, dqc " << rime_q << "\n";
            // Cloud N
            res[Nc_idx] -= rime_n;

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Sublimination, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] += delta_e;

            // ice multiplication based on Hallet and Mossop
            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                codi::RealReverse mult_n = C_mult * mult_1 * mult_2 * rime_q;
                codi::RealReverse mult_q = mult_n * cc.ice.min_x;
                mult_q = min(rime_q, mult_q);

                // Ice
                res[qi_idx] += mult_q;
                // Ice N
                res[Ni_idx] += mult_n;
                resq -= mult_q;
            }

            // Enhancement of melting
            if(T_prime > tmelt && enhanced_melting)
            {
                codi::RealReverse melt_q = (T_prime-tmelt)*const5*rime_q;
                codi::RealReverse melt_n = melt_q/x_1;

                melt_q = min(q1, melt_q);
                melt_n = min(N1, melt_n);

                resq -= melt_q;
                resn -= melt_n;
                // Rain
                res[qr_idx] += melt_q;
                // Rain N
                res[Nr_idx] += melt_n;

                codi::RealReverse delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
                // Melting, cooling
                if(melt_q > 0.0)
                    res[lat_cool_idx] -= delta_e;
                // Freezing, heating
                else
                    res[lat_heat_idx] -= delta_e;
            }
        }
    };

    auto particle_rain_riming = [&](
        codi::RealReverse &q1,
        codi::RealReverse &N1,
        codi::RealReverse &resq,
        codi::RealReverse &resn,
        collection_model_constants_t &coeffs,
        particle_model_constants_t &pc1)
    {
        if(qr_prime > q_crit && q1 > q_crit)
        {

            codi::RealReverse x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x, cc.rain.max_x);
            codi::RealReverse d_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);
            codi::RealReverse x_1 = particle_mean_mass(q1, N1, pc1.min_x, pc1.max_x);
            codi::RealReverse d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);

            codi::RealReverse v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;
            codi::RealReverse v_r = particle_velocity(x_r, cc.rain.a_vel, cc.rain.b_vel) * cc.rain.rho_v;

            codi::RealReverse rime_n = M_PI/4.0 * N1 * Nr
                * (coeffs.delta_n_aa * d_1 * d_1
                    + coeffs.delta_n_ab * d_1 * d_r
                    + coeffs.delta_n_bb * d_r * d_r)
                * sqrt(coeffs.theta_n_aa * v_1 * v_1
                    - coeffs.theta_n_ab * v_1 * v_r
                    + coeffs.theta_n_bb * v_r * v_r);
            codi::RealReverse rime_q = M_PI/4.0 * N1 * qr_prime
                * (coeffs.delta_q_aa * d_1 * d_1
                    + coeffs.delta_q_ab * d_1 * d_r
                    + coeffs.delta_q_bb * d_r * d_r)
                * sqrt(coeffs.theta_q_aa * v_1 * v_1
                    - coeffs.theta_q_ab * v_1 * v_r
                    + coeffs.theta_q_bb * v_r * v_r);
            rime_q = min(qr_prime, rime_q);
            rime_n = min(Nr, rime_n);

            resq += rime_q;
            // Rain
            res[qr_idx] -= rime_q;
            // Rain N
            res[Nr_idx] -= rime_n;

            codi::RealReverse delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            // ice multiplication based on Hallet and Mossop
            if(T_prime < tmelt && ice_multiplication)
            {
                codi::RealReverse mult_1 = (T_prime - T_mult_min)*const3;
                codi::RealReverse mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                codi::RealReverse mult_n = C_mult * mult_1 * mult_2 * rime_q;
                codi::RealReverse mult_q = mult_n * cc.ice.min_x;
                mult_q = min(rime_q, mult_q);

                // Ice
                res[qi_idx] += mult_q;
                // Ice N
                res[Ni_idx] += mult_n;
                resq -= mult_q;
            }

            // Enhancement of melting
            if(T_prime > tmelt && enhanced_melting)
            {
                codi::RealReverse melt_q = (T_prime-tmelt)*const5*rime_q;
                codi::RealReverse melt_n = melt_q/x_1;

                melt_q = min(q1, melt_q);
                melt_n = min(N1, melt_n);

                resq -= melt_q;
                resn -= melt_n;
                // Rain
                res[qr_idx] += melt_q;
                // Rain N
                res[Nr_idx] += melt_n;

                codi::RealReverse delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
                // Melting, cooling
                if(melt_q > 0.0)
                    res[lat_cool_idx] -= delta_e;
                // Freezing, heating
                else
                    res[lat_heat_idx] -= delta_e;
            }
        }
    };
    //// hail cloud rimin
    particle_cloud_riming(qh_prime, Nh, res[qh_idx], res[Nh_idx], cc.coeffs_hcr, cc.rain);

    //// hail rain riming
    particle_rain_riming(qh_prime, Nh, res[qh_idx], res[Nh_idx], cc.coeffs_hrr, cc.rain);

    //// graupel cloud riming
    particle_cloud_riming(qg_prime, Ng, res[qg_idx], res[Ng_idx], cc.coeffs_gcr, cc.rain);

    //// graupel rain riming
    particle_rain_riming(qg_prime, Ng, res[qg_idx], res[Ng_idx], cc.coeffs_grr, cc.rain);

    ////////////// rain_freeze_gamlook
    // Freezing of rain and conversion to ice, graupel, hail
    codi::RealReverse xmax_ice = pow( pow(D_rainfrz_ig/cc.rain.a_geo, 1.0/cc.rain.b_geo), cc.rain.mu );
    codi::RealReverse xmax_gr  = pow( pow(D_rainfrz_gh/cc.rain.a_geo, 1.0/cc.rain.b_geo), cc.rain.mu );
    if(T_prime < T_freeze)
    {
        codi::RealReverse fr_q, fr_n, fr_n_i, fr_q_i, fr_n_g, fr_q_g, fr_n_h,
            fr_q_h, fr_n_tmp, fr_q_tmp;
        fr_q = fr_n = fr_n_i = fr_q_i = fr_n_g = fr_q_g = fr_n_h =
            fr_q_h = fr_n_tmp = fr_q_tmp = 0.0;
        codi::RealReverse Nr_tmp = Nr;
        if(qr_prime <= q_crit_fr)
        {
            if(T_prime < T_f) // instantaneous freezing
            {
                fr_q = fr_q_i = qr_prime;
                fr_n = fr_n_i = Nr;
                fr_n_tmp = fr_q_tmp = 1.0;
            }
        } else
        {
            codi::RealReverse x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x, cc.rain.max_x);
            Nr_tmp = qr_prime/x_r;
            if(T_prime < T_f)
            {
                // Not sure if this branch makes too much sense
                fr_q = qr_prime;
                fr_n = Nr_tmp;
                codi::RealReverse lam = pow(rain_g1/rain_g2*x_r, -cc.rain.mu);
                codi::RealReverse N0 = cc.rain.mu * Nr_tmp * pow(lam, rain_nm1) / rain_g1;
                codi::RealReverse tmp = lam*xmax_ice;
                fr_n_i = N0/pow(cc.rain.mu*lam, rain_nm1) * table_r1.look_lo(tmp);
                fr_q_i = N0/pow(cc.rain.mu*lam, rain_nm2) * table_r2.look_lo(tmp);
                tmp = lam*xmax_gr;
                fr_n_i = N0/pow(cc.rain.mu*lam, rain_nm1) * table_r1.look_lo(tmp);
                fr_n_i = N0/pow(cc.rain.mu*lam, rain_nm2) * table_r2.look_lo(tmp);

                fr_n_h = fr_n - fr_n_g;
                fr_q_h = fr_q - fr_q_g;
                fr_n_g = fr_n_g - fr_n_i;
                fr_q_g = fr_q_g - fr_q_i;
                fr_n_tmp = Nr_tmp/max(fr_n, Nr_tmp);
                fr_q_tmp = qr_prime/max(fr_q, qr_prime);
            } else
            {
                // heterogeneous freezing
                codi::RealReverse j_het = max(b_HET * ( exp(a_HET * (T_prime-tmelt)) - 1.0 ), 0.0) / rho_w;

                if(j_het >= 1.0e-20/dt)
                {
                    fr_n = j_het * qr_prime;
                    fr_q = j_het * qr_prime * x_r * cc.rain.c_z; // rain_coeffs
                    codi::RealReverse lam = pow(cc.rain.g1/cc.rain.g2 * x_r, -cc.rain.mu);
                    codi::RealReverse N0 = cc.rain.mu * Nr_tmp * pow(lam, cc.rain.nm1) / cc.rain.g1;
                    // codi::RealReverse arg1 = lam*xmax_ice;
                    codi::RealReverse tmp = lam*xmax_ice;
                    fr_n_i = j_het * N0/pow(cc.rain.mu*lam, cc.rain.nm2) * table_r2.look_lo(tmp);
                    fr_q_i = j_het * N0/pow(cc.rain.mu*lam, cc.rain.nm3) * table_r3.look_lo(tmp);
                    // arg1 = lam*xmax_gr;
                    tmp = lam*xmax_gr;
                    fr_n_i = j_het * N0/pow(cc.rain.mu*lam, cc.rain.nm2) * table_r2.look_lo(tmp);
                    fr_n_i = j_het * N0/pow(cc.rain.mu*lam, cc.rain.nm3) * table_r3.look_lo(tmp);

                    fr_n_h = fr_n - fr_n_g;
                    fr_q_h = fr_q - fr_q_g;
                    fr_n_g = fr_n_g - fr_n_i;
                    fr_q_g = fr_q_g - fr_q_i;
                    fr_n_tmp = Nr_tmp/max(fr_n, Nr_tmp);
                    fr_q_tmp = qr_prime/max(fr_q,qr_prime);
                } else
                {
                    fr_n = fr_q = fr_n_i = fr_q_i = fr_n_g = fr_q_g
                        = fr_n_h = fr_q_h = fr_n_tmp = fr_q_tmp = 0.0;
                }
            }

            fr_n = fr_n * fr_n_tmp;
            fr_q = fr_q * fr_q_tmp;
            fr_n_i = fr_n_i * fr_n_tmp;
            fr_n_g = fr_n_g * fr_n_tmp;
            fr_n_h = fr_n_h * fr_n_tmp;
            fr_q_i = fr_q_i * fr_q_tmp;
            fr_q_g = fr_q_g * fr_q_tmp;
            fr_q_h = fr_q_h * fr_q_tmp;
        }
        // Rain
        res[qr_idx] -= fr_q;
        // Rain N
        res[Nr_idx] = Nr_tmp - fr_n;

        // Snow
        res[qs_idx] += fr_q_i;
        // Snow N
        res[Ns_idx] += fr_n_i;

        // Graupel
        res[qg_idx] += fr_q_g;
        // Graupel N
        res[Ng_idx] += fr_n_g;

        // Hail
        res[qh_idx] += fr_q_h;
        // Hail N
        res[Nh_idx] += fr_n_h;

        codi::RealReverse delta_e = latent_heat_melt(T_prime) * fr_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(fr_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] += delta_e;

    }

    /////////// ice melting
    if(T_prime > tmelt && qi_prime > 0.0)
    {
        codi::RealReverse x_i = particle_mean_mass(qi_prime, Ni, cc.ice.min_x, cc.ice.max_x);
        // Complete melting
        codi::RealReverse melt_q = qi_prime;
        codi::RealReverse melt_n = Ni;
        qi = 0.0;
        qi_prime = 0.0;
        Ni = 0;

        // Ice
        res[qi_idx] = 0.0;
        y[qi_idx] = 0.0;
        // Ice N
        res[Ni_idx] = 0.0;
        y[Ni_idx] = 0.0;
        // melt into cloud or rain
        if(x_i > cc.cloud.max_x)
        {
            // Rain
            res[qr_idx] += melt_q;
            // Rain N
            res[Nr_idx] += melt_n;
        } else
        {
            // Cloud
            res[qc_idx] += melt_q;
            // DEBUG STATEMENT
            if(melt_q/ref.qref >= 0.003)
                std::cout << "line 1900, dqc " << melt_q << "\n";
            // Cloud N
            res[Nc_idx] += melt_n;
        }

        codi::RealReverse delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }

    /////////// snow melting
    if(T_prime > tmelt && qs_prime > 0.0)
    {
        codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime);

        codi::RealReverse x_s = particle_mean_mass(qs_prime, Ns, cc.snow.min_x, cc.snow.max_x);
        codi::RealReverse d_s = particle_diameter(x_s, cc.snow.a_geo, cc.snow.b_geo);
        codi::RealReverse v_s = particle_velocity(x_s, cc.snow.a_vel, cc.snow.b_vel) * cc.snow.rho_v;

        codi::RealReverse fv_q = cc.snow.a_f + cc.snow.b_f * sqrt(v_s*d_s);
        // From Rasmussen and Heymsfield (1987)
        codi::RealReverse fh_q = 1.05 * fv_q;
        codi::RealReverse melt = 2.0*M_PI/L_ew * d_s * Ns;
        codi::RealReverse melt_h = melt * K_T * (T_prime - tmelt);
        codi::RealReverse melt_v = melt * dv0 * L_wd/Rv * (p_sat/T_prime - p_sat_melt/tmelt);
        codi::RealReverse melt_q = (melt_h * fh_q + melt_v * fv_q);

        codi::RealReverse melt_n = min(max( (melt_q-qs_prime)/x_s + Ns, 0.0), Ns);
        melt_q = min(qs_prime, max(melt_q, 0.0));
        melt_n = min(Ns, max(melt_n, 0.0));

        // Spontaneous melting at 10 C
        if(T_prime - tmelt > 10.0)
        {
            melt_q = qs_prime;
            melt_n = Ns;
        }

        // Snow
        res[qs_idx] -= melt_q;
        // Snow N
        res[Ns_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;

        codi::RealReverse delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
        // TODO: Check if we need that
        // snow%n(i,k) = MAX(snow%n(i,k), snow%q(i,k)/snow%x_max)
    }

    ////////////// Melting of graupel and hail (either simple or LWF-based)
    //// graupel melting
    if(T_prime > tmelt && qg_prime > 0.0)
    {
        codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime);
        codi::RealReverse x_g = particle_mean_mass(qg_prime, Ng, cc.graupel.min_x, cc.graupel.max_x);
        codi::RealReverse d_g = particle_diameter(x_g, cc.graupel.a_geo, cc.graupel.b_geo);
        codi::RealReverse v_g = particle_velocity(x_g, cc.graupel.a_vel, cc.graupel.b_vel) * cc.graupel.rho_v;

        codi::RealReverse fv_q = cc.graupel.a_f + cc.graupel.b_f * sqrt(v_g*d_g);
        codi::RealReverse fh_q = 1.05 * fv_q;
        codi::RealReverse melt = 2.0*M_PI/L_ew * d_g * Ng;
        codi::RealReverse melt_h = melt * K_T * (T_prime - tmelt);
        codi::RealReverse melt_v = melt * dv0 * L_wd/Rv * (p_sat/T_prime - p_sat_melt/tmelt);
        codi::RealReverse melt_q = (melt_h * fh_q + melt_v * fv_q);

        codi::RealReverse melt_n = min(max( (melt_q-qg_prime)/x_g + Ng, 0.0), Ng);
        melt_q = max(0.0, min(melt_q, qg_prime));
        melt_n = max(0.0, max(melt_n, Ng));

        // Graupel
        res[qg_idx] -= melt_q;
        // Graupel N
        res[Ng_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;

        codi::RealReverse delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
    //// hail melting (simple)
    if(T_prime > tmelt && qh_prime > 0.0)
    {
        codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime);
        codi::RealReverse x_h = particle_mean_mass(qh_prime, Nh, cc.hail.min_x, cc.hail.max_x);
        codi::RealReverse d_h = particle_diameter(x_h, cc.hail.a_geo, cc.hail.b_geo);
        codi::RealReverse v_h = particle_velocity(x_h, cc.hail.a_vel, cc.hail.b_vel) * cc.hail.rho_v;

        codi::RealReverse fv_q = cc.hail.a_f + cc.hail.b_f * sqrt(v_h*d_h);
        codi::RealReverse fh_q = 1.05 * fv_q;
        codi::RealReverse melt = 2.0*M_PI/L_ew * d_h * Nh;
        codi::RealReverse melt_h = melt * K_T * (T_prime - tmelt);
        codi::RealReverse melt_v = melt * dv0 * L_wd/Rv * (p_sat/T_prime - p_sat_melt/tmelt);
        codi::RealReverse melt_q = (melt_h * fh_q + melt_v * fv_q);

        codi::RealReverse melt_n = min(max( (melt_q-qh_prime)/x_h + Nh, 0.0), Nh);
        melt_q = max(0.0, min(melt_q, qh_prime));
        melt_n = max(0.0, max(melt_n, Nh));

        // Hail
        res[qh_idx] -= melt_q;
        // Hail N
        res[Nh_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;

        codi::RealReverse delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }

    ////////////// Evaporation from melting ice particles
    // TODO: Check if deposition rates are set to zero here
    auto evaporation = [&](
        codi::RealReverse &q1,
        codi::RealReverse &N1,
        codi::RealReverse &resq,
        particle_model_constants_t &pc1)
    {
        if(q1 > 0.0 && T_prime > tmelt)
        {
            codi::RealReverse p_d = qv_prime * Rv * T_prime;
            codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime);
            codi::RealReverse s_sat = p_d/p_sat - 1.0;
            // TODO: Check D_v
            codi::RealReverse g_d = 4.0*M_PI / (L_wd*L_wd
                / (K_T * Rv * tmelt*tmelt*tmelt) + Rv*tmelt /(D_v * p_sat));

            codi::RealReverse x_1 = particle_mean_mass(q1, N1, pc1.min_x, pc1.max_x);
            codi::RealReverse d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);
            codi::RealReverse v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;

            codi::RealReverse f_v = pc1.a_f + pc1.b_f * sqrt(v_1*d_1);

            codi::RealReverse delta_q = g_d * N1 * pc1.c_s * d_1 * f_v * s_sat;
            delta_q = min(q1, max(-delta_q, 0.0));

            // Vapor
            res[qv_idx] += delta_q;
            resq -= delta_q;
            // subl? Actually not used here but for the broader simulation
            // We could use it later if we want to see how much evaporates
            // sublq -= delta_q;

            codi::RealReverse delta_e = latent_heat_ice(T_prime) * delta_q / latent_heat_melt(T_prime);
            // Sublimination, cooling
            if(delta_q > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] -= delta_e;
        }
    };
    evaporation(qs_prime, Ns, res[qs_idx], cc.snow);
    evaporation(qg_prime, Ng, res[qg_idx], cc.graupel);
    evaporation(qi_prime, Ni, res[Ni_idx], cc.ice);

    ////////////// Warm rain, ie autoconversion, accretion and rain rain collision
    if(auto_type == 1)
    {
        // autoconversionKB
        codi::RealReverse k_a = 6.0 + 25 * pow(9.59, -1.7);
        codi::RealReverse x_s_i = 1.0/cc.cloud.max_x;
        codi::RealReverse x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x, cc.cloud.max_x);
        // Using Beheng 1994
        codi::RealReverse au = k_a * pow(x_c*1e3, 3.3) * pow(qc_prime*1e3, 1.4) * 1e3;
        au = min(qc_prime/cc.dt, au);
        res[Nr_idx] += au*x_s_i;
        res[qr_idx] += au;
        res[Nc_idx] -= au*x_s_i*2.0;
        res[qc_idx] -= au;

        // DEBUG STATEMENT
        if(au/ref.qref >= 0.003)
            std::cout << "line 2099, dqc " << au << "\n";

        // accretionKB
        if(qc_prime > q_crit_i && qr_prime > q_crit_i)
        {
            codi::RealReverse ac = 6.0 * qc_prime * qr_prime; // k_r = 6.0 from Beheng (1994)
            ac = min(qc_prime/cc.dt, ac);
            res[qr_idx] += ac;
            res[qc_idx] -= ac;
            // DEBUG STATEMENT
            if(ac/ref.qref >= 0.003)
                std::cout << "line 2110, dqc " << ac << "\n";
        }
    } else if(auto_type == 2)
    {
        // Not implemented since it appears to be not very interesting
    } else
    {
        const double EPSILON = 1.0e-25;
        // autoconversionSB
        if(qc_prime > q_crit)
        {
            const double k_1 = 6.0e2;
            const double k_2 = 0.68;
            codi::RealReverse x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x, cc.cloud.max_x);
            codi::RealReverse au = cloud_k_au * qc_prime*qc_prime * x_c*x_c * cc.cloud.rho_v;
            codi::RealReverse tau = min(max(1.0-qc_prime/(qc_prime+qr_prime+EPSILON), EPSILON), 0.9);
            codi::RealReverse phi = k_1 * pow(tau, k_2) * pow(1.0-pow(tau, k_2), 3);
            au *= (1.0 + phi/pow(1.0-tau, 2));
            au = max(min(qc_prime, au), 0.0);

            codi::RealReverse sc = cloud_k_sc * qc_prime*qc_prime * cc.cloud.rho_v;

            res[qr_idx] += au;
            res[Nr_idx] += au / cc.cloud.max_x;
            res[qc_idx] -= au;
            res[Nc_idx] -= min(Nc, sc);
        }

        // accretionSB
        if(qc_prime > 0.0 && qr_prime > 0.0)
        {
            const double k_1 = 5.0e-4;
            const double k_r = 5.78;
            codi::RealReverse tau = min(max(1.0-qc_prime/(qc_prime+qr_prime+EPSILON), EPSILON), 1.0);
            codi::RealReverse phi = pow(tau/(tau+k_1), 4);
            codi::RealReverse ac = k_r * qc_prime * qr_prime * phi;
            ac = min(qc_prime, ac);
            codi::RealReverse x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x, cc.cloud.max_x);
            res[qr_idx] += ac;
            res[qc_idx] -= ac;
            // DEBUG STATEMENT
            if(ac/ref.qref >= 0.003)
                std::cout << "line 2152, dqc " << ac << "\n";
            res[Nc_idx] -= min(Nc, x_c);
        }
    }

    // rain_selfcollectionSB
    // Seifert and Beheng (2001)
    if(qr_prime > 0)
    {
        codi::RealReverse x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x, cc.rain.max_x);
        codi::RealReverse D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);
        // Parameters based on Seifert (2008)
        codi::RealReverse sc = 4.33 * Nr * qr_prime * cc.rain.rho_v; // rhain%rho_v(i, k)
        // Breakup Seifert (2008), Eq. A13
        codi::RealReverse breakup = 0.0;
        if(D_r > 0.30e-3)
            breakup = sc * (1.0e+3 * (D_r - 1.10e-3) + 1.0);
        res[Nr_idx] -= min(Nr, sc-breakup);
    }

    // rain_evaporation
    // Seifert (2008)
    if(s_sw < 0.0 && qr_prime > 0.0 && qc_prime < q_crit)
    {
        codi::RealReverse D_v = diffusivity(T_prime, p_prime);
        // Equation (A2) of 10.1175/2008JAS2586.1
        codi::RealReverse g_d = 2*M_PI /
            ( (R_v*T_prime)/(D_v*p_sat) + (L_wd*L_wd)/(K_T*R_v*T_prime*T_prime) );
        codi::RealReverse x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x, cc.rain.max_x);
        codi::RealReverse D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);

        codi::RealReverse mue;
        // Equation 20 of Seifert (2008)
        if(D_v <= cc.rain.cmu3)
            mue = cc.rain.cmu0 * tanh( pow(4.0*cc.rain.cmu2*(D_v-cc.rain.cmu3),
                cc.rain.cmu5)) + cc.rain.cmu4;
        else
            mue = cc.rain.cmu1 * tanh( pow(cc.rain.cmu2*(D_v-cc.rain.cmu3),
                cc.rain.cmu5)) + cc.rain.cmu4;
        // Equation A8
        codi::RealReverse lambda = pow(M_PI/6.0*rho_w*(mue+3.0)*(mue+2.0)*(mue+1.0)/x_r, 1.0/3.0);

        // Approximation of Gamma(mue+5/2) / Gamma(mue+2)
        codi::RealReverse gamma_approx = 0.1357940435E+01
            + mue * (0.3033273220E+00
            + mue * (-0.1299313363E-01
            + mue * (0.4002257774E-03
            - mue * 0.4856703981E-05)));
        // ventilation factor (A7) with (A5) and (A9)
        // Approximation for terminal fall velocity of raindrops
        codi::RealReverse f_v = a_v + b_v * pow(N_Sc, 1.0/3.0) * gamma_approx
            * sqrt(a_prime/kin_visc_air * cc.rain.rho_v/lambda)
            * (1.0 - 0.5 * (b_prime/a_prime) * pow(lambda/(c_prime+lambda), (mue+2.5))
                - 0.125 * pow(b_prime/a_prime, 2.0) * pow(lambda/(2*c_prime+lambda), (mue+2.5))
                - 1.0/16.0 * pow(b_prime/a_prime, 3.0) * pow(lambda/(3*c_prime+lambda), (mue+2.5))
                - 5.0/127.0 * pow(b_prime/a_prime, 4.0) * pow(lambda/(4*c_prime+lambda), (mue+2.5))
            );
        codi::RealReverse gamma_eva;
        // D_br = 1.1e-3 from ICON
        if(rain_gfak > 0.0)
            gamma_eva = rain_gfak * (1.1e-3/D_r) * exp(-0.2*mue);
        else
            gamma_eva = 1.0;

        // Equation A5 with A9
        codi::RealReverse delta_qv = g_d * Nr * (mue+1.0) / lambda * f_v * s_sw; // (mue+1.0) / lambda *
        codi::RealReverse delta_nv = gamma_eva * delta_qv/x_r;

        delta_qv = max(-delta_qv, 0.0);
        delta_nv = max(-delta_nv, 0.0);
        delta_qv = min(delta_qv, qv);
        delta_nv = min(delta_nv, Nv);

        res[qv_idx] += delta_qv;
        res[qr_idx] -= delta_qv;
        res[Nr_idx] -= delta_nv;
        // evap -= delta_qv

        codi::RealReverse delta_e = latent_heat_evap(T_prime) * delta_qv / specific_heat_water_vapor(T_prime);
        // Evaporation, cooling
        if(delta_qv > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
    ///////// Sedimentation with Seifert Beheng (2006, 2008)
    ///////// sedimentation_explicit
    // from mo_2mom_mcrph_driver.f90
    // using an explicit flux-form semi-lagrangian scheme (actually used after the microphysics)
    //// sedi_icon_rain
    codi::RealReverse cmax = 0.0;
    codi::RealReverse p_prime_pascal = p_prime*100;

    // using compute_rhoa(p_prime,T_prime,S) instead of cc.rho_a_prime makes
    // next to no difference
    codi::RealReverse rhocorr = pow(cc.rho_a_prime/rho_0, -rho_vel);
    codi::RealReverse v_n_sedi = 0.0;
    codi::RealReverse v_q_sedi = 0.0;

    auto sedi_icon_core = [&](
        codi::RealReverse &q,
        codi::RealReverse &N,
        codi::RealReverse &v_q_sedi,
        codi::RealReverse &v_n_sedi,
        codi::RealReverse &resQ,
        codi::RealReverse &resN,
        codi::RealReverse &resOut)
    {
        codi::RealReverse cmax_tmp = cmax;
        codi::RealReverse q_flux = 0.0;
        codi::RealReverse n_flux = 0.0;
        codi::RealReverse v_nv = v_n_sedi;
        codi::RealReverse v_qv = v_q_sedi;  // percentage how much trickles down
        // Assuming v_nv, v_qv is negative
        codi::RealReverse c_nv = -v_nv * inv_z;//  * dt; // times inverse of layer thickness...
        codi::RealReverse c_qv = -v_qv * inv_z;//  * dt;
        cmax_tmp = max(cmax_tmp, c_qv);

        // TODO: Check why we should multiply with dt here. That *should* be
        // handled outside in rk4.h or at least I thought so.
        codi::RealReverse s_nv = 0.0;
        if(c_nv <= 1.0)
            s_nv = v_nv*N*inv_z;

        codi::RealReverse s_qv = 0.0;
        if(c_qv <= 1.0)
            s_qv = v_qv*q*inv_z;

        // abs is used for paranoia reasons and should never be needed
        resN -= abs(s_nv);
        resQ -= abs(s_qv);
        resOut -= abs(s_qv);
        cmax = cmax_tmp;
        // precrate = -abs(s_sq);
    };

    auto sedi_icon_sphere = [&](
        codi::RealReverse &q,
        codi::RealReverse &N,
        codi::RealReverse &resQ,
        codi::RealReverse &resN,
        codi::RealReverse &resOut,
        particle_model_constants_t &pc)
    {
        codi::RealReverse v_n_sedi = 0.0;
        codi::RealReverse v_q_sedi = 0.0;

        if(q > q_crit)
        {
            codi::RealReverse x = particle_mean_mass(q, N, pc.min_x, pc.max_x);
            codi::RealReverse lam = pow(pc.lambda*x, pc.b_vel);
            codi::RealReverse v_n = max(pc.alfa_n * lam, pc.vsedi_min);
            codi::RealReverse v_q = max(pc.alfa_q * lam, pc.vsedi_min);
            v_n = min(v_n, pc.vsedi_max);
            v_q = min(v_q, pc.vsedi_max);
            v_n *= rhocorr;
            v_q *= rhocorr;

            v_n_sedi -= v_n;
            v_q_sedi -= v_q;
        }
        sedi_icon_core(q, N, v_q_sedi, v_n_sedi, resQ, resN, resOut);

        resN = max(min(N, q/pc.min_x), q/pc.max_x);
    };

    auto sedi_icon_sphere_lwf = [&](
        )
    {

    };

    if(qr_prime > q_crit)
    {
        codi::RealReverse x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x, cc.rain.max_x);
        codi::RealReverse D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);
        codi::RealReverse mue = (qc_prime >= q_crit)
            ? (cc.rain.nu+1.0)/cc.rain.b_geo - 1.0
            : rain_mue_dm_relation(D_r, cc.rain.cmu0, cc.rain.cmu1,
                cc.rain.cmu2, cc.rain.cmu3, cc.rain.cmu4);
        // inverse of lambda in Eq. (A10) SB (2008)
        codi::RealReverse D_p = D_r * pow( (mue+3)*(mue+2)*(mue+1), -1.0/3.0 );

        // SB (2008), Eq. (A10)
        codi::RealReverse v_nr = cc.rain.alpha - cc.rain.beta
            * pow( 1.0+cc.rain.gamma*D_p, -mue-1.0);
        codi::RealReverse v_qr = cc.rain.alpha - cc.rain.beta
            * pow( 1.0+cc.rain.gamma*D_p, -mue-4.0);

        // Seifert (2008), Eq. A10
        // rhocorr: height dependency of particle fall speed
        v_nr *= rhocorr;
        v_qr *= rhocorr;

        v_n_sedi -= v_nr;
        v_q_sedi -= v_qr;
    }

    sedi_icon_core(qr_prime, Nr, v_q_sedi, v_n_sedi, res[qr_idx], res[Nr_idx], res[qr_out_idx]);

    // sedi_icon_sphere ice
    if(qi_prime > 0.0)
        sedi_icon_sphere(qi_prime, Ni, res[qi_idx], res[Ni_idx], res[qi_out_idx], cc.ice);

    // sedi_icon_sphere snow
    if(qs_prime > 0.0)
        sedi_icon_sphere(qs_prime, Ns, res[qs_idx], res[Ns_idx], res[qs_out_idx], cc.snow);

    // sedi_icon_sphere graupel
    const bool lprogmelt = false; // TODO: Check if true or false
    if(qg_prime > 0.0)
    {
        if(lprogmelt)
            sedi_icon_sphere_lwf();
        else
            sedi_icon_sphere(qg_prime, Ng, res[qg_idx], res[Ng_idx], res[qg_out_idx], cc.graupel);
    }

    // sedi_icon_sphere hail
    if(qh_prime > 0.0)
    {
        if(lprogmelt)
            sedi_icon_sphere_lwf();
        else
            sedi_icon_sphere(qh_prime, Nh, res[qh_idx], res[Nh_idx], res[qh_out_idx], cc.hail);
    }
    //////// end sedimentation_explicit

    ////// Get back to non-prime
    // DEBUG STATEMENT
    if(res[qc_idx]/ref.qref >= 0.003)
        std::cout << "line 2384, dqc " << res[qc_idx] << "\n";
    res[qc_idx] /= ref.qref;

    res[qr_idx] /= ref.qref;
    res[qv_idx] /= ref.qref;
    res[qi_idx] /= ref.qref;
    res[qs_idx] /= ref.qref;
    res[qg_idx] /= ref.qref;
    res[qr_out_idx] /= ref.qref;

     // Compute parameters
    codi::RealReverse psat_prime = saturation_pressure_water(T_prime);
    codi::RealReverse qsat_prime = Epsilon*( psat_prime/(p_prime - psat_prime) );
    codi::RealReverse p_div_T_prime = p_prime / T_prime;
    codi::RealReverse cpv_prime = specific_heat_water_vapor(T_prime);
    codi::RealReverse cpa_prime = specific_heat_dry_air(T_prime);
    codi::RealReverse cpl_prime = specific_heat_water(T_prime);
    codi::RealReverse rhow_prime = density_water(T_prime);
    codi::RealReverse L_prime = latent_heat_water(T_prime);
    codi::RealReverse H_prime =
        1.0/(((L_prime/(Rv*T_prime)) - 1.0)
        *(L_prime/(thermal_conductivity_dry_air(T_prime)*T_prime))
        + ((Rv*T_prime)/(cc.alpha_d*diffusivity(T_prime,p_prime)*psat_prime)));
    codi::RealReverse c_prime =
        4.0*M_PI*H_prime*pow((3.0/(4.0*M_PI*rhow_prime))
        *cc.Nc_prime*cc.Nc_prime , 1.0/3.0);
    codi::RealReverse qc_third = pow(qc , 1.0/3.0);
    codi::RealReverse qr_delta1 = pow(qr , cc.delta1);
    codi::RealReverse qr_delta2 = pow(qr , cc.delta2);

    // Compute nondimensional coefficients
    codi::RealReverse C1 = (ref.tref*gravity_acc*ref.wref)/(Ra*ref.Tref);
    codi::RealReverse C2 = ((1.0-Epsilon)*ref.qref)/Epsilon;
    codi::RealReverse C3 = (cpv_prime*ref.qref)/cpa_prime;
    codi::RealReverse C4 = (cpl_prime*ref.qref)/cpa_prime;
    codi::RealReverse C5 = (gravity_acc*ref.wref*ref.tref)/(ref.Tref*cpa_prime);
    codi::RealReverse C6 = (ref.tref*L_prime*c_prime*pow(ref.qref,1.0/3.0))/(ref.Tref*cpa_prime);
    codi::RealReverse C7 = (ref.tref*L_prime*cc.e1_prime
        *pow(ref.qref,cc.delta1))/(ref.Tref*cpa_prime);
    codi::RealReverse C8 = (ref.tref*L_prime*cc.e2_prime
        *pow(ref.qref,cc.delta2))/(ref.Tref*cpa_prime);
    // codi::RealReverse B = ref.tref*nc.QRin;
    codi::RealReverse C9 = (ref.tref*c_prime)/pow(ref.qref , 2.0/3.0);
    codi::RealReverse C12 = ref.tref*cc.e1_prime*pow(ref.qref , cc.delta1-1.0);
    codi::RealReverse C13 = ref.tref*cc.e2_prime*pow(ref.qref , cc.delta2-1.0);
    codi::RealReverse C15 = Epsilon/ref.qref;
    codi::RealReverse C16 = L_prime/(Rv*ref.Tref);

    if(fixed)
    {
        res[p_idx] = 0;
        res[T_idx] = 0;
        res[w_idx] = 0;
    } else
    {
        res[p_idx] = -( C1/(1.0 + C2*(qv/(1.0 + qv_prime))) )*( (p*w)/T );
        res[T_idx] = ( 1.0/(1.0 + C3*qv + C4*(qc + qr)) )*( -C5*w + C6*qc_third*(S-1.0)
            + (C7*qr_delta1 + C8*qr_delta2)*min(S-1.0,0.0) );
        res[w_idx] += cc.dw;
        res[lat_heat_idx] = ( 1.0/(1.0 + C3*qv + C4*(qc + qr)) ) * C6*qc_third*(S-1.0);
    }
    res[S_idx] = (S/p)*res[p_idx] - (S/qv)*( 1.0 - (qv/(C15+qv)) )*( C9*qc_third*(S-1.0)
        + (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0, 0.0) ) - C16*(S/(T*T))*res[T_idx];
}

void RHS_SB_no_ice(std::vector<codi::RealReverse> &res,
	 std::vector<codi::RealReverse> &y,
	 const reference_quantities_t &ref,
	 model_constants_t &cc,
     nc_parameters_t& nc,
     bool fixed=false )
{

    // Decrypt the variables
    codi::RealReverse p = y[p_idx];
    codi::RealReverse T = y[T_idx];
    codi::RealReverse w = y[w_idx];

    codi::RealReverse S = y[S_idx];
    codi::RealReverse qc = y[qc_idx];
    codi::RealReverse qr = y[qr_idx];
    codi::RealReverse qv = y[qv_idx];
    codi::RealReverse Nc = y[Nc_idx];
    codi::RealReverse Nr = y[Nr_idx];
    codi::RealReverse Nv = y[Nv_idx];
    codi::RealReverse qi = y[qi_idx];
    codi::RealReverse Ni = y[Ni_idx];
    codi::RealReverse vi = y[vi_idx];
    codi::RealReverse qs = y[qs_idx];
    codi::RealReverse Ns = y[Ns_idx];
    codi::RealReverse qg = y[qg_idx];
    codi::RealReverse Ng = y[Ng_idx];
    codi::RealReverse qh = y[qh_idx];
    codi::RealReverse Nh = y[Nh_idx];

    for(auto &r: res) r = 0;

    // Safety measure: ensure positiveness
    if(0. > qc)
        qc = 0.0;

    if(0. > qr)
        qr = 0.0;

    if(0. > qv)
        qv = 0.0;

    if(0. > Nc)
        Nc = 0.0;

    if(0. > Nr)
        Nr = 0.0;

    if(0. > Nv)
        Nv = 0.0;

    // Change to dimensional variables
    codi::RealReverse p_prime = ref.pref * p;
    codi::RealReverse T_prime = ref.Tref * T;
    codi::RealReverse w_prime = ref.wref * w;
    codi::RealReverse qc_prime = ref.qref * qc;
    codi::RealReverse qr_prime = ref.qref * qr;
    codi::RealReverse qv_prime = ref.qref * qv;

    codi::RealReverse T_c = T_prime - tmelt;

    ////////////// Warm rain, ie autoconversion, accretion and rain rain collision
    // autoconversionKB
    codi::RealReverse k_a = 6.0 + 25 * pow(9.59, -1.7);
    codi::RealReverse x_s_i = 1.0/cc.cloud.max_x;
    codi::RealReverse x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x, cc.cloud.max_x);
    // Using Beheng 1994
    codi::RealReverse au = k_a * pow(x_c*1e3, 3.3) * pow(qc_prime*1e3, 1.4) * 1e3;
    au = min(qc_prime/cc.dt, au);
    res[Nr_idx] += au*x_s_i;
    res[qr_idx] += au;
    res[Nv_idx] -= au*x_s_i*2.0;
    res[qc_idx] -= au;

    // accretionKB (= growth of rain hydrometeor by collision with cloud drops)
    if(qc_prime > q_crit_i && qr_prime > q_crit_i)
    {
        codi::RealReverse ac = 6.0 * qc_prime * qr_prime; // k_r = 6.0 from Beheng (1994)
        ac = min(qc_prime/cc.dt, ac);
        res[qr_idx] += ac;
        res[qc_idx] -= ac;
    }

    // rain_selfcollectionSB
    // Seifert and Beheng (2001)
    if(qr > 0)
    {
        codi::RealReverse x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x, cc.rain.max_x);
        codi::RealReverse D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);
        // Parameters based on Seifert (2008)
        codi::RealReverse sc = 4.33 * Nr * qr_prime * cc.rain.rho_v; // rhain%rho_v(i, k)
        // Breakup Seifert (2008), Eq. A13
        codi::RealReverse breakup = 0.0;
        if(D_r > 0.30e-3)
            breakup = sc * (1.0e+3 * (D_r - 1.10e-3) + 1.0);
        res[Nr_idx] -= min(Nr, sc-breakup);
    }

    // rain_evaporation
    // Seifert (2008)
    codi::RealReverse e_d = qv_prime * R_v * T_prime;
    codi::RealReverse e_sat = saturation_pressure_water_icon(T_prime);
    codi::RealReverse s_sw = e_d / e_sat - 1.0;
    if(s_sw < 0.0 && qr_prime > 0.0 && qc_prime < q_crit_i)
    {
        codi::RealReverse D_v = diffusivity(T_prime, p_prime);
        // Equation (A2) of 10.1175/2008JAS2586.1
        codi::RealReverse g_d = 2*M_PI /
            ( (R_v*T_prime)/(D_v*e_sat) + (L_wd*L_wd)/(K_T*R_v*T_prime*T_prime) );
        codi::RealReverse x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x, cc.rain.max_x);
        codi::RealReverse D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);

        codi::RealReverse mue;
        // Equation 20 of Seifert (2008)
        if(D_v <= cc.rain.cmu3)
            mue = cc.rain.cmu0 * tanh( pow(4.0*cc.rain.cmu2*(D_v-cc.rain.cmu3),
                cc.rain.cmu5)) + cc.rain.cmu4;
        else
            mue = cc.rain.cmu1 * tanh( pow(cc.rain.cmu2*(D_v-cc.rain.cmu3),
                cc.rain.cmu5)) + cc.rain.cmu4;
        // Equation A8
        codi::RealReverse lambda = pow(M_PI/6.0*rho_w*(mue+3.0)*(mue+2.0)*(mue+1.0)/x_r, 1.0/3.0);

        // Approximation of Gamma(mue+5/2) / Gamma(mue+2)
        codi::RealReverse gamma_approx = 0.1357940435E+01
            + mue * (0.3033273220E+00
            + mue * (-0.1299313363E-01
            + mue * (0.4002257774E-03
            - mue * 0.4856703981E-05)));
        // ventilation factor (A7) with (A5) and (A9)
        codi::RealReverse f_v = a_v + b_v * pow(N_Sc, 1.0/3.0) * gamma_approx
            * sqrt(a_prime/kin_visc_air * cc.rain.rho_v/lambda)
            * (1.0 - 0.5 * (b_prime/a_prime) * pow(lambda/(c_prime+lambda), (mue+2.5))
                - 0.125 * pow(b_prime/a_prime, 2.0) * pow(lambda/(2*c_prime+lambda), (mue+2.5))
                - 1.0/16.0 * pow(b_prime/a_prime, 3.0) * pow(lambda/(3*c_prime+lambda), (mue+2.5))
                - 5.0/127.0 * pow(b_prime/a_prime, 4.0) * pow(lambda/(4*c_prime+lambda), (mue+2.5))
            );
        codi::RealReverse gamma_eva;
        // D_br = 1.1e-3 from ICON
        if(rain_gfak > 0.0)
            gamma_eva = rain_gfak * (1.1e-3/D_r) * exp(-0.2*mue);
        else
            gamma_eva = 1.0;

        // Equation A5 with A9
        codi::RealReverse delta_qv = g_d * Nr * (mue+1.0) / lambda * f_v * s_sw;
        codi::RealReverse delta_nv = gamma_eva * delta_qv/x_r;

        delta_qv = max(-delta_qv, 0.0);
        delta_nv = max(-delta_nv, 0.0);
        delta_qv = min(delta_qv, qv);
        delta_nv = min(delta_nv, Nv);

        res[qv_idx] += delta_qv;
        res[qr_idx] -= delta_qv;
        res[Nr_idx] -= delta_nv;

        // res[qc_idx] *= delta_qv;
        // How about?
        // res[] += delta_nv;
        //  evap(i,k) = evap(i,k)-eva_q
    }

     // Compute parameters
    codi::RealReverse psat_prime = saturation_pressure_water(T_prime);
    codi::RealReverse qsat_prime = Epsilon*( psat_prime/(p_prime - psat_prime) );
    codi::RealReverse p_div_T_prime = p_prime / T_prime;
    codi::RealReverse cpv_prime = specific_heat_water_vapor(T_prime);
    codi::RealReverse cpa_prime = specific_heat_dry_air(T_prime);
    codi::RealReverse cpl_prime = specific_heat_water(T_prime);
    codi::RealReverse rhow_prime = density_water(T_prime);
    codi::RealReverse L_prime = latent_heat_water(T_prime);
    codi::RealReverse H_prime =
        1.0/(((L_prime/(Rv*T_prime)) - 1.0)
        *(L_prime/(thermal_conductivity_dry_air(T_prime)*T_prime))
        + ((Rv*T_prime)/(cc.alpha_d*diffusivity(T_prime,p_prime)*psat_prime)));
    codi::RealReverse c_prime =
        4.0*M_PI*H_prime*pow((3.0/(4.0*M_PI*rhow_prime))
        *cc.Nc_prime*cc.Nc_prime , 1.0/3.0);
    codi::RealReverse qc_third = pow(qc , 1.0/3.0);
    codi::RealReverse qr_delta1 = pow(qr , cc.delta1);
    codi::RealReverse qr_delta2 = pow(qr , cc.delta2);

    // Compute nondimensional coefficients
    codi::RealReverse C1 = (ref.tref*gravity_acc*ref.wref)/(Ra*ref.Tref);
    codi::RealReverse C2 = ((1.0-Epsilon)*ref.qref)/Epsilon;
    codi::RealReverse C3 = (cpv_prime*ref.qref)/cpa_prime;
    codi::RealReverse C4 = (cpl_prime*ref.qref)/cpa_prime;
    codi::RealReverse C5 = (gravity_acc*ref.wref*ref.tref)/(ref.Tref*cpa_prime);
    codi::RealReverse C6 = (ref.tref*L_prime*c_prime*pow(ref.qref,1.0/3.0))/(ref.Tref*cpa_prime);
    codi::RealReverse C7 = (ref.tref*L_prime*cc.e1_prime
        *pow(ref.qref,cc.delta1))/(ref.Tref*cpa_prime);
    codi::RealReverse C8 = (ref.tref*L_prime*cc.e2_prime
        *pow(ref.qref,cc.delta2))/(ref.Tref*cpa_prime);
    codi::RealReverse B = ref.tref*nc.QRin;

    // Pressure
    res[p_idx] = -( C1/(1.0 + C2*(qv/(1.0 + qv_prime))) )*( (p*w)/T );

    // Temperature
     res[T_idx] = ( 1.0/(1.0 + C3*qv + C4*(qc + qr)) )*( -C5*w + C6*qc_third*(S-1.0)
        + (C7*qr_delta1 + C8*qr_delta2)*min(S-1.0,0.0) );

    // Saturation
    // double psat_prime = saturation_pressure_water(nc.t[traj*countp[1] + i]);
    res[S_idx] = res[qv_idx] * p / (psat_prime*qv+Epsilon*psat_prime)
        + ( res[p_idx] * (psat_prime*qv+Epsilon*psat_prime) + p*psat_prime*res[qv_idx] )
        / ( (psat_prime * qv + Epsilon * psat_prime)*(psat_prime * qv + Epsilon * psat_prime) );
    if(fixed)
    {
        res[p_idx] = 0;
        res[T_idx] = 0;
    }

    // sedimentation_explicit
    // from mo_2mom_mcrph_driver.f90
    // using an explicit flux-form semi-lagrangian scheme (actually used after the microphysics)
    // sedi_icon_rain
    if(qr_prime > q_crit_i)
    {
        codi::RealReverse x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x, cc.rain.max_x);
        codi::RealReverse D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);
        codi::RealReverse mue = (qc_prime >= q_crit_i)
            ? (cc.rain.nu+1.0)/cc.rain.b_geo - 1.0
            : rain_mue_dm_relation(D_r, cc.rain.cmu0, cc.rain.cmu1,
                cc.rain.cmu2, cc.rain.cmu3, cc.rain.cmu4);
        codi::RealReverse D_p = D_r * pow( (mue+3)*(mue+2)*(mue+1), -1.0/3.0 );
        codi::RealReverse v_nr = cc.rain.alpha - cc.rain.beta
            * pow( 1.0+cc.rain.gamma*D_p, -mue-1.0);
        codi::RealReverse v_qr = cc.rain.alpha - cc.rain.beta
            * pow( 1.0+cc.rain.gamma*D_p, -mue-4.0);

        codi::RealReverse rhocorr = pow(compute_rhoa(p_prime, T_prime, S)/1.225, -0.4);
        // codi::RealReverse rhocorr = pow(cc.rho_a_prime/1.225, -0.4);
        v_nr *= rhocorr;
        // TODO: Check the scaling here. I'm not convinced that this is okay as it is!
        v_qr *= rhocorr*ref.qref*ref.qref;

        res[qr_idx] -= v_qr;
        res[Nr_idx] -= v_nr;
    }

    // Get back to non-prime
    res[qc_idx] /= ref.qref;
    res[qr_idx] /= ref.qref;
    res[qv_idx] /= ref.qref;
    res[qi_idx] /= ref.qref;
    res[qs_idx] /= ref.qref;
    res[qg_idx] /= ref.qref;
}



#endif


// mo_satad.f90
// #ifdef __COSMO__
//    ! diagnostically (re)compute density of moist air for time-level nnow
//    ! ... using specific moisture quantities
//    CALL calrho( te, ppe,qve,qce,qle+qie, p0e, rhotot, idim, kdim, r_d, rvd_m_o )
// #endif

//     DO k = klo, kup
//       DO i = ilo , iup

//         ! total content of the species which are changed by the adjustment:
//         qw(i,k) = qve(i,k) + qce(i,k)

//         ! check, which points will still be subsaturated even
//         ! if all the cloud water would have been evaporated.
//         ! At such points, the Newton iteration is not necessary and the
//         ! adjusted values of T, p, qv and qc can be obtained directly.

//         lwdocvd(i,k) = ( lwd + (cp_v - cl)*(te(i,k)-tmelt) - r_v*te(i,k) )/ cvd

//         Ttest(i,k) = te(i,k) - lwdocvd(i,k)*qce(i,k)

//         qtest(i,k) = qsat_rho(Ttest(i,k), rhotot(i,k))

//       END DO
//     END DO

//     DO k = klo, kup
//       DO i = ilo , iup

//        ll_satad = .TRUE.
//        IF (lqtvar) THEN
//          IF ( qtvar(i,k) > 0.0001_ireals * (qve(i,k)+qce(i,k))**2 ) THEN ! satad not in EDMF boundary layer
//            ll_satad = .false.                                            ! sqrt(qtvar) > 0.001
//          ENDIF
//        ENDIF

//        IF ( ll_satad ) THEN

//         IF (qw(i,k) <= qtest(i,k) ) THEN
//           ! In this case, all the cloud water evaporates and there is still (sub)saturation.
//           ! The resulting state depends only on the available cloud water and is
//           ! not saturated, which enables direct computation of the adjusted variables:
//           evap(i,k) = evap(i,k)-qce(i,k)
//           qve(i,k)  = qw(i,k)
//           qce(i,k)  = 0.0_ireals
//           te(i,k)   = Ttest(i,k)
//         ELSE
//           ! In this case, the Newton interation is needed
//           nsat       = nsat+1
//           iwrk(nsat) = i
//           kwrk(nsat) = k
//           ! Field for the iterated temperature, here set the starting value for the below iteration
//           ! to the "old" temperature (as an alternative, the arithmetic mean
//           ! between "old" temperature and dew point has been tested, but did
//           ! not significantly increase the convergence speed of the iteration):
//           twork(nsat)  = te(i,k)
//           ! And this is the storage variable for the "old" values in the below iteration:
//           ! Add some nonesense increment to the starting, which is sufficient to trigger the
//           ! iteration below:
//           tworkold(nsat) = twork(nsat) + 10.0_ireals
//           qv_tmp(nsat) = qve(i,k)
//         END IF

//        END IF

//       END DO
//     END DO