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

/** @defgroup ufunc User Functions
 * Functions to be defined by the user, specifically the definition of
 * the right-hand side function RHS of the ODE
 * \f[ y' = \text{RHS}(y) \f]
 * which is solved.
 * @{
 */

codi::RealReverse trunc(codi::RealReverse x)
{
    int left = x.getValue();
    codi::RealReverse r = left;
    r.setGradient(x.getGradient());
    return r;
}

/**
 * This function evaluates the RHS function of the ODE. It uses the 1 moment
 * cloud scheme.
 *
 * @param res On out: system state difference
 * @param y Old system state
 * @param ref Reference quantities to transform between units
 * @param cc Pointer to constants from the model
 * @param nc Pointer to parameters from the netCDF file
 */
void RHS(std::vector<codi::RealReverse> &res,
	 std::vector<codi::RealReverse> &y,
	 const reference_quantities_t &ref,
	 model_constants_t &cc,
     nc_parameters_t& nc)
{
  // ============================================================
  // This is the RHS for the water cloud
  // ============================================================

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
  codi::RealReverse psat_prime = saturation_pressure_water_icon(T_prime);
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

} // End of RHS

/**
 * This function evaluates the RHS function of the ODE. It uses the 1 moment
 * cloud scheme. Updates only pressure and temperature.
 *
 * @param res On out: system state difference
 * @param y Old system state
 * @param ref Pointer to reference quantities to transform between units
 * @param cc Pointer to constants from the model
 */
void Press_Temp(
    std::vector<codi::RealReverse> &res,
    std::vector<codi::RealReverse> &y,
    const reference_quantities_t &ref,
    model_constants_t &cc)
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
  codi::RealReverse psat_prime = saturation_pressure_water_icon(T_prime);
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

/**
 * Calculates saturation adjustement as in COSMO.
 * It corrects the temperature, specific humidity (qv),
 * the cloud water content and pressure for condensation and
 * evaporation.
 * From ICON:
 * Saturation adjustment at constant total density (adjustment of T and p accordingly)
 * assuming chemical equilibrium of water and vapor. For the heat capacity of
 * of the total system (dry air, vapor, and hydrometeors) the value of dry air
 * is taken, which is a common approximation and introduces only a small error.
 * From COMSO:
 * Diese Routine fuehrt in JEDEM Fall eine Adjustierung der H2O-Komponenten und
 * der Temperatur auf ihre Gleichgewichtswerte hin durch. Dabei wird angenommen,
 * dass sofort bei 100% Saettigung Kondensation eintritt und dass die Phasenum-
 * wandlungen vollstaendig reversibel ablaufende thermodynamische Prozesse sind
 *
 * Die Saettigungsadjustierung ist immer der LETZTE Schritt eines Wolkenmodells
 * Alle anderen Prozesse, wie Advektion, Diffusion und Mikrophysik werden davor
 * berechnet. Die Saettigungsadjustierung bestimmt dann nur noch die Nukleation
 * die max. zum Erreichen von 100% rel. Feuchte im Wolkenraum erforderlich ist.
 *
 * Die tatsaechliche rel. Feuchte in Wolken kann mit der Subroutine saturation_
 * evaluation ermittelt werden. Da das direkte Adjustierungsverfahren eine Nae-
 * herung darstellt, wird nicht immer exakt 100% Feuchte erreicht.  Erfahrungs-
 * gemaess liegen die Werte aber ungefaehr im Intervall
 *
 *                        [ 99.85 % <= RH <= 100.03 % ]
 *
 * oder sogar noch naeher bei 100%. Da auch in der Realitaet ganz leichte Unter
 * und Uebersaettigungen in der Wolke beobachtet werden, ist dieses Resultat OK
 *
 * Die hier programmierte Saettigungsadjustierung setzt die Dampfdruckformel in
 * der Form nach Tetens (1930) voraus. Siehe UP clouds sowie die Referenz Te30!
 */
template<class float_t>
void saturation_adjust(
    const float_t &T_prime,
    const float_t &p_prime,
    const float_t &p_sat,
    const float_t &qv_prime,
    const float_t &qc_prime,
    std::vector<float_t> &res,
    const double &qref)
{
    if(T_prime > T_f || always_sat_adj)
    {
        float_t lrcp = L_wd / cp;
        float_t q_sat = water_vapor_sat_ratio(p_prime, T_prime);
        // Specific saturation humidity after Dotzek (4.33)
        // float_t q_sat = (r_const/r1_const) / (p_prime/p_sat + (r_const/r1_const)-1);
        // Adjust according to Dotzek page 35
        float_t a1 = p_sat_const_a * (T_freeze-p_sat_const_b)
            / ((T_prime - p_sat_const_b)*(T_prime - p_sat_const_b));
        float_t b1 = a1 * lrcp * q_sat;
        // Dotzek page 36
        float_t delta_q = -(qv_prime - q_sat)/(1 + b1);
        // Vaporisation until nothing is left or saturation is reached
        delta_q = max(min(delta_q, qc_prime), -qv_prime);

        res[qv_idx] += delta_q;
        res[qc_idx] -= delta_q;

#ifdef TRACE_QV
        if(trace)
        {
            std::cout << "saturation_adjust: dqv " << delta_q << "\n";
            float_t alt_delta_q = -(qv_prime - q_sat * qref)/(1 + b1);
            alt_delta_q = max(min(alt_delta_q, qc_prime), -qv_prime);
            std::cout << "q_sat " << q_sat << "\n"
                      << "qv_prime " << qv_prime << "\n"
                      << "delta_q with ref " << alt_delta_q << "\n";
        }
#endif
#ifdef TRACE_QC
        if(trace)
            std::cout << "saturation_adjust: dqc " << -delta_q << "\n";
#endif
    }
}

// TODO: Try meteo_utilities.f90 satad
template<class float_t>
void saturation_adjust_meteo(
    const float_t &T_prime,
    const float_t &p_prime,
    const float_t &p_sat,
    const float_t &qv_prime,
    const float_t &qc_prime,
    std::vector<float_t> &res,
    const double &qref)
{
    float_t qv = qv_prime;
    float_t qc = qc_prime;
    float_t T = T_prime;
    // exchanged with convert_S_to_qv(p, T, S)
    // This method calculates the saturated amount of qv using p_sat according
    // to Dotzek. We ignore that and calculate saturation differently.
    auto q_sat_f = [](
        const float_t &p_sat,
        const float_t &p)
    {
        return Epsilon * p_sat / (p - (1-Epsilon)*p_sat);
    };

    auto temp_p_dt = [](
        const float_t &T,
        const float_t p)
    {
        return p_sat_const_a*(T_sat_low_temp-p_sat_const_b)
            * ( 1+(Epsilon-1)*p ) * p / ( (T-p_sat_const_b)*(T-p_sat_const_b) );
    };

    T = T - L_wd*qc_prime/cp;
    qv = qc_prime + qv_prime;
    qc = 0;
    // After Dotzek...
    // float_t p_sat_2 = saturation_pressure_water_icon(T);
    // float_t delta_q = q_sat_f(p_sat_2, p_prime) - qv;
    // ... our version
    float_t delta_q = convert_S_to_qv(p_prime, T, float_t(1)) - qv;
    if(delta_q < 0)
    {

        // Handle over saturation
        // Adjust temperature
        float_t Th = cp*T + L_wd*qv;
        // After Dotzek...
        // p_sat_2 = saturation_pressure_water_icon(T_prime);
        // float_t T_qd0 = q_sat_f(p_sat_2, p_prime);
        // ... our version
        float_t T_qd0 = convert_S_to_qv(p_prime, T_prime, float_t(1));
        float_t T_dt0;
        // Do the Newton
        for(uint32_t i=0; i<1; ++i)
        {
            T_dt0 = temp_p_dt(T_prime, T_qd0);
            T = ( Th-L_wd*(T_qd0-T_dt0*T_prime) )
                / ( cp+L_wd*T_dt0 );
            // After Dotzek ...
            // p_sat_2 = saturation_pressure_water_icon(T);
            // T_qd0 = q_sat_f(p_sat_2, p_prime);
            // ... our version
            T_qd0 = convert_S_to_qv(p_prime, T, float_t(1));
        }

        // Get ratio mixing "back"
        T_dt0 = temp_p_dt(T, T_qd0);
        float_t T_gn = ( Th - L_wd*(T_qd0-T_dt0*T) )
            / ( cp+L_wd*T_dt0 );
        T_qd0 += T_dt0*(T_gn-T);

        qv = T_qd0;
        qc = max(qc_prime+qv_prime - T_qd0, 1.0e-20);
        delta_q = qv - qv_prime;
        res[qv_idx] += delta_q;
        res[qc_idx] -= delta_q;
        // T <- T_gn
        res[T_idx] +=  T_gn - T_prime;
#ifdef TRACE_QV
        if(trace)
            std::cout << "saturation_adjust_meteo (over saturated for new qv = qv+qc) dqv " << delta_q << "\n";

#endif
#ifdef TRACE_ENV
        if(trace)
            std::cout << "saturation_adjust_meteo dT " << T_gn - T_prime
                  << "\n";
#endif
#ifdef TRACE_QC
        if(trace)
            std::cout << "saturation_adjust_meteo (over saturated for new qv = qv+qc) dqc " << -delta_q << "\n";
#endif
    } else
    {
        // Throw everything from qc at qv -> qv will be as saturated as possible
        res[qv_idx] += qc_prime;
        res[qc_idx] -= qc_prime;
        res[T_idx] += T-T_prime;
#ifdef TRACE_QV
        if(trace)
            std::cout << "saturation_adjust_meteo (under saturated) dqv " << qc_prime << "\n";
#endif
#ifdef TRACE_QC
        if(trace)
            std::cout << "saturation_adjust_meteo (under saturated) dqc " << -qc_prime << "\n";
#endif
#ifdef TRACE_ENV
        if(trace)
            std::cout << "saturation_adjust_meteo (under saturated) dT " << T-T_prime << "\n";
#endif
    }
}


/**
 * CCN activation after Hande et al 2015.
 */
template<class float_t>
void ccn_act_hande(
    float_t &p_prime,
    float_t &w_prime,
    float_t &T_prime,
    float_t &qv_prime,
    float_t &qc_prime,
    float_t &Nc,
    const double &EPSILON,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    // non maritime case
    if(qc_prime > EPSILON && w_prime > 0.0)
    {
        float_t acoeff = a_ccn[0] * atan(b_ccn[0] * p_prime - c_ccn[0]) + d_ccn[0];
        float_t bcoeff = a_ccn[1] * atan(b_ccn[1] * p_prime - c_ccn[1]) + d_ccn[1];
        float_t ccoeff = a_ccn[2] * atan(b_ccn[2] * p_prime - c_ccn[2]) + d_ccn[2];
        float_t dcoeff = a_ccn[3] * atan(b_ccn[3] * p_prime - c_ccn[3]) + d_ccn[3];
        // concentration of ccn
        float_t nuc_n = acoeff * atan(bcoeff * log(w_prime) + ccoeff) + dcoeff;

        // we need to substract the already "used" ccn in cloud droplets
        // the rest can create new cloud droplets
        // float_t Nc_tmp = qv_prime / cc.cloud.min_x;
        float_t delta_n = max(nuc_n, 10.0e-6) - Nc;
        /////// DT PROBLEM
        // Problem here: delta_n shall be the difference between *should be*
        // nuc_n and *is* Nc but this needs to be scaled with our timestep
        // if the timestep is bigger than 1s or else we overshoot.
        if(cc.dt_prime > 1) delta_n /= cc.dt_prime;
        float_t delta_q;
        // if(cc.dt_prime > 1)
        delta_q = min(delta_n * cc.cloud.min_x_act, qv_prime/cc.dt_prime);
        // else
        //     delta_q = min(delta_n * cc.cloud.min_x_act, qv_prime);
        // Or: CCN activation is instantanuous, so scale it in any case
        // delta_n /= cc.dt_prime;
        /////// DT PROBLEM SOLVED
        delta_n = max(delta_n, 0.0);

        delta_n = delta_q / cc.cloud.min_x_act;

        res[Nc_idx] += delta_n;
        res[qc_idx] += delta_q;
        res[qv_idx] -= delta_q;
#ifdef TRACE_QC
        if(trace)
            // if(abs(delta_q) > 0)
            std::cout << "Ascent dqc " << delta_q << ", dNc " << delta_n
                << ", nuc_n " << nuc_n << ", Nc " << Nc << ", rest " << 10.0e-6 << "\n";
#endif
#ifdef TRACE_QV
        if(trace)
            std::cout << "Ascent dqv " << -delta_q << "\n";
#endif

        float_t delta_e = latent_heat_evap(T_prime) * delta_q / specific_heat_water_vapor(T_prime);
        // Evaporation
        if(delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        else
            res[lat_heat_idx] += delta_e;
    }
}


/**
 * CCN activation after Seifert & Beheng (2006):
 * cloud_nucleation(..) from Cosmo 5.2
 * There are different cases. For this sumulation, we choose
 * HUCM continental case (Texas CCN)HUCM continental case (Texas CCN)
 * From COSMO 5.2 S = 100 * s_w(...)
 * where s_w = e_v/e_ws(T_2mom) - 1.0
 * e_v = q * T_2mom * r_v
 * q = qv * rho
 * rho = total density of moist air (kg/m3) = (p * p_pertubation) / (r_d * T * (1.0+ (r_v/r_d - 1)*q_v - q_c - (q_r + q_s)))
 * Using gas law rho = p/r_d*tv, r_d gas constant of dry air, pp pressure pertubation
 * T_2mom = T    (I think)
 * r_v = gas constant for water vapour
 * r_d = gas constant for dry air
 * e_ws = saturation pressure over water = e_3 * EXP (A_w * (t_ - T_3) / (t_ - B_w))
 */
template<class float_t>
void ccn_act_seifert(
    float_t &p_prime,
    float_t &p,
    float_t &T_prime,
    float_t &T,
    float_t &qv_prime,
    float_t &qv,
    float_t &qc_prime,
    float_t &qc,
    float_t &Nc,
    float_t &qr,
    float_t &z_prime,
    const double &dt,
    float_t &w,
    float_t &S,
    float_t &psat_prime,
    const reference_quantities_t &ref,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    double N_ccn = 1260e06;
    double N_max = 3000e06;
    double N_min =  300e02; // e06?
    double S_max = 20.000; // in percentage
    double k_ccn = 0.308;
    // Parameter for exp decrease of CCN with height
    double z_nccn = 99999;

    // float_t S_percentage = s_sw*100;
    //     saturation_pressure_ice(T_prime)/saturation_pressure_water_icon(T_prime) * (S+1.0);
    // S_percentage *= 100;
    // oversaturation in percentage
    float_t S_percentage = (S-1)*100;
    // float_t S_percentage = s_sw * 100;

    // In contrast to Seifert & Beheng, we do not approximate dS/dt via
    // w dS/dz but calculate dS/dt
    float_t qsat_prime = Epsilon*( psat_prime/(p_prime - psat_prime) );
    float_t p_div_T_prime = p_prime / T_prime;
    float_t cpv_prime = specific_heat_water_vapor(T_prime);
    float_t cpa_prime = specific_heat_dry_air(T_prime);
    float_t cpl_prime = specific_heat_water(T_prime);
    float_t rhow_prime = density_water(T_prime);
    float_t L_prime = latent_heat_water(T_prime);
    float_t H_prime =
        1.0/(((L_prime/(Rv*T_prime)) - 1.0)
        *(L_prime/(thermal_conductivity_dry_air(T_prime)*T_prime))
        + ((Rv*T_prime)/(cc.alpha_d*diffusivity(T_prime,p_prime)*psat_prime)));
    float_t c_prime =
        4.0*M_PI*H_prime*pow((3.0/(4.0*M_PI*rhow_prime))*Nc*Nc , 1.0/3.0);
    float_t qc_third = pow(qc , 1.0/3.0);
    float_t qr_delta1 = pow(qr , cc.delta1);
    float_t qr_delta2 = pow(qr , cc.delta2);

    float_t C1 = (ref.tref*gravity_acc*ref.wref)/(Ra*ref.Tref);
    float_t C2 = ((1.0-Epsilon)*ref.qref)/Epsilon;
    float_t C3 = (cpv_prime*ref.qref)/cpa_prime;
    float_t C4 = (cpl_prime*ref.qref)/cpa_prime;
    float_t C5 = (gravity_acc*ref.wref*ref.tref)/(ref.Tref*cpa_prime);
    float_t C6 = (ref.tref*L_prime*c_prime*pow(ref.qref,1.0/3.0))/(ref.Tref*cpa_prime);
    float_t C7 = (ref.tref*L_prime*cc.e1_prime
        *pow(ref.qref,cc.delta1))/(ref.Tref*cpa_prime);
    float_t C8 = (ref.tref*L_prime*cc.e2_prime
        *pow(ref.qref,cc.delta2))/(ref.Tref*cpa_prime);
    // float_t B = ref.tref*nc.QRin;
    float_t C9 = (ref.tref*c_prime)/pow(ref.qref , 2.0/3.0);
    float_t C12 = ref.tref*cc.e1_prime*pow(ref.qref , cc.delta1-1.0);
    float_t C13 = ref.tref*cc.e2_prime*pow(ref.qref , cc.delta2-1.0);
    float_t C15 = Epsilon/ref.qref;
    float_t C16 = L_prime/(Rv*ref.Tref);

    float_t dp = -( C1/(1.0 + C2*(qv/(1.0 + qv_prime))) )*( (p*w)/T );
    float_t dT = ( 1.0/(1.0 + C3*qv + C4*(qc + qr)) )*( -C5*w + C6*qc_third*(S-1.0)
        + (C7*qr_delta1 + C8*qr_delta2)*min(S-1.0,0.0) );
    float_t dS = (S/p)*dp - (S/qv)*( 1.0 - (qv/(C15+qv)) )*( C9*qc_third*(S-1.0)
        + (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0, 0.0) ) - C16*(S/(T*T))*dT;
    dS = (S/p)*dp - C16*(S/(T*T))*dT;
#ifdef TRACE_SAT
    if(trace)
        std::cout << "Saturation (CCN activation) dS: " << dS
                    << ", times dt: " << dS * dt
                    << ", Evaporation: " << (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0, 0.0)
                    << ", Condensation: " << C9*qc_third*(S-1.0)
                    << ", Last part: " << - C16*(S/(T*T))*res[T_idx]
                    << ", First part: " << (S/p)*res[p_idx]
                    << ", Middle part: " << - (S/qv)*( 1.0 - (qv/(C15+qv)) )
                    << ", S_percentage: " << S_percentage
                    << ", S: " << S
                    << ", e_d: " << qv_prime * Rv * T_prime
                    << ", qv_prime: " << qv_prime
                    << ", Rv: " << Rv
                    << ", T_prime " << T_prime
                    << ", p_sat: " << saturation_pressure_water_icon(T_prime)
                    << ", s_sw = e_d/p_sat - 1: " << qv_prime * Rv * T_prime/saturation_pressure_water_icon(T_prime)-1
                    << ", S calculated: " << qv_prime * Rv * T_prime/saturation_pressure_water_icon(T_prime)
                    << "\n";
#endif
    // dS *= dt;
    // float_t dSdz_w = w_prime *
#ifdef TRACE_QC
    // if(trace)
    // {
    //     std::cout << "Saturation (CCN activation) dS: " << dS
    //                 << ", times dt: " << dS * dt
    //                 << ", Evaporation: " << (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0, 0.0)
    //                 << ", Condensation: " << C9*qc_third*(S-1.0)
    //                 << ", Last part: " << - C16*(S/(T*T))*res[T_idx]
    //                 << ", First part: " << (S/p)*res[p_idx]
    //                 << ", Middle part: " << - (S/qv)*( 1.0 - (qv/(C15+qv)) )
    //                 << ", S_percentage: " << S_percentage
    //                 << ", S: " << S
    //                 << ", e_d: " << qv_prime * Rv * T_prime
    //                 << ", qv_prime: " << qv_prime
    //                 << ", Rv: " << Rv
    //                 << ", T_prime " << T_prime
    //                 << ", p_sat: " << saturation_pressure_water_icon(T_prime)
    //                 << ", s_sw = e_d/p_sat - 1: " << qv_prime * Rv * T_prime/saturation_pressure_water_icon(T_prime)-1
    //                 << ", S calculated: " << qv_prime * Rv * T_prime/saturation_pressure_water_icon(T_prime)
    //                 << "\n";
    //     std::cout << "dS " << dS << ", s_perc " << S_percentage << ", S " << S
    //               << ", S_max " << S_max << ", T " << T << ", T_prime "
    //               << T_prime << ", T_f " << T_f << "\n";
    //     std::cout << "N_ccn " << N_ccn << ", k_ccn " << k_ccn
    //               << ", pow(S, k-1) " << pow(S_percentage, k_ccn-1.0)
    //               << ", N*k*s^k-1 " << N_ccn * k_ccn * pow(S_percentage, k_ccn-1.0)
    //               << "\n";
    //     std::cout << "Nc " << Nc << "\n";
    //     std::cout << "supersat: " << saturation_pressure_ice(T_prime)/saturation_pressure_water_icon(T_prime) * (S-1)
    //               << ", e_d/e_sat " << qv_prime * R_v * T_prime / saturation_pressure_water_icon(T_prime) - 1.0
    //               << ", e_d/e_sat 2: " << qv_prime * Rv * T_prime / saturation_pressure_water_icon(T_prime) - 1.0
    //               << ", qv to S " << convert_qv_to_S(p_prime, T_prime, qv_prime) << "\n";
    //     std::cout << "z_prime " << z_prime << ", z_nccn " << z_nccn << "\n";
    //     float_t S_cosmo = min(10, 100*S_percentage);
    //     float_t n_ccn_texas = 1260e6 * pow(S_cosmo, 0.308);
    //     std::cout << "over saturation in perc with cosmo " << S_cosmo
    //               << ", n_ccn_texas " << n_ccn_texas << "\n";
    // }
#endif
    // dS can be considered above zero all the time w is above zero
    if(S_percentage > 0.0 && dS >= 0.0 && S_percentage <= S_max && T_prime >= T_f)
    {
        float_t delta_n = 0;
        float_t delta_q = 0;

        if(z_prime <= z_nccn)
        {
            delta_n = N_ccn * k_ccn * pow(S_percentage, k_ccn-1.0) * dS;
        } else
        {
            delta_n = ( N_ccn * k_ccn * pow(S_percentage, k_ccn-1.0) * dS
                - N_ccn/z_nccn * pow(S_percentage, k_ccn) )
                * exp( (z_nccn - z_prime)/z_nccn );
        }
#ifdef TRACE_QV
        // std::cout << "SB ccn activation no bound, no max dqv " << -delta_q << "\n";
#endif
#ifdef TRACE_QC
        // std::cout << "SB ccn activation no bound, no max dqc " << delta_q << ", dNc " << delta_n << "\n";
#endif
        delta_n = max(min(delta_n, N_max-Nc), 0.0);
        delta_q = delta_n * cc.cloud.min_x_act;

        if(delta_q > qv_prime)
        {
            delta_q = qv_prime;
            delta_n = delta_q/cc.cloud.min_x_act;
        }
        res[Nc_idx] += delta_n;
        res[qc_idx] += delta_q;
        res[qv_idx] -= delta_q;
#ifdef TRACE_QV
        if(trace)
            std::cout << "SB ccn activation dqv " << -delta_q << "\n";
#endif
#ifdef TRACE_QC
        if(trace)
            std::cout << "SB ccn activation no bound dqc " << delta_q << ", dNc " << delta_n
                << "\n with bound? " << (res[Nc_idx]*dt + Nc > N_max) << " dNc "
                << (N_max - Nc)/dt
                << "\nelse ? " << (res[Nc_idx]*dt + Nc < N_min) << " dNc "
                << (N_min - Nc)/dt << "\n";


        // std::cout << "SB ccn activation before Nc " << Nc << " res[Nc] " << res[Nc_idx]
        //           << "\nres*dt + Nc " <<  res[Nc_idx]*dt + Nc  << " > "
        //           << N_max << " ? " <<  (res[Nc_idx]*dt + Nc > N_max)
        //           << "\nN_min: " << N_min << ", with true: "
        //           << (N_min - Nc)*dt << "\n";
#endif
        // Hard upper and lower limits
        if(res[Nc_idx]*dt + Nc > N_max)
        {
            res[Nc_idx] = (N_max - Nc)/dt;
        } else if(res[Nc_idx]*dt + Nc < N_min)
        {
            res[Nc_idx] = (N_min - Nc)/dt;
        }
#ifdef TRACE_QV
        // if(trace)
            // std::cout << "SB ccn activation dqv " << -delta_q << "\n";
#endif
        float_t delta_e = latent_heat_evap(T_prime) * delta_q / specific_heat_water_vapor(T_prime);
        // Evaporation
        if(delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        else
            res[lat_heat_idx] += delta_e;
    }
}

/**
 * (optional) homogenous nucleation using KHL06
 */
template<class float_t>
void ice_nuc_hom(
    float_t &T_prime,
    float_t &w_prime,
    float_t &p_prime,
    float_t &qv_prime,
    float_t &qi_prime,
    float_t &Ni,
    float_t &S_i,
    float_t &p_sat_ice,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    float_t s_crit = 2.349 - T_prime * (1.0/259.0);
    const double r_0 = 0.25e-6;         // aerosol particle radius prior to freezing
    const double alpha_d = 0.5;         // deposition coefficient (KL02; Spichtinger & Gierens 2009)
    const double M_w = 18.01528e-3;     // molecular mass of water [kg/mol]
    const double M_a = 28.96e-3;        // molecular mass of air [kg/mol]
    const double ma_w = M_w/N_avo;      // mass of water molecule [kg]
    const double svol = ma_w / rho_ice; // specific volume of a water molecule in ice

    if(S_i > s_crit && T_prime < 235.0 && Ni < ni_hom_max)
    {
        float_t x_i = particle_mean_mass(qi_prime, Ni,
            cc.ice.min_x_nuc_homo, cc.ice.max_x);
        float_t r_i = pow(x_i/(4.0/3.0*M_PI*rho_ice), 1.0/3.0);
        float_t v_th = sqrt(8.0*k_b*T_prime/(M_PI*ma_w));
        float_t flux = alpha_d * v_th/4.0;
        float_t n_sat = p_sat_ice/(k_b*T_prime);

        // coeffs of supersaturation equation
        std::vector<float_t> acoeff(3);
        acoeff[0] = (L_ed * grav) / (cp * Rv * T_prime*T_prime) - grav/(Ra * T_prime);
        acoeff[1] = 1.0/n_sat;
        acoeff[2] = (L_ed*L_ed * M_w * ma_w)/(cp * p_prime * T_prime * M_a);

        // coeffs of depositional growth equation
        std::vector<float_t> bcoeff(2);
        bcoeff[0] = flux * svol * n_sat * (S_i - 1.0);
        bcoeff[1] = flux/diffusivity(T_prime, p_prime);

        // pre-existing ice crystals included as reduced updraft speed
        float_t ri_dot = bcoeff[0] / (1.0 + bcoeff[1] * r_i);
        float_t R_ik = (4.0*M_PI) / svol * Ni * r_i*r_i * ri_dot;
        float_t w_pre = max(0.0, (acoeff[1] + acoeff[2] * S_i)/(acoeff[0]*S_i)*R_ik); // KHL06 Eq. 19

        // homogenous nucleation event
        if(w_prime > w_pre)
        {
            float_t cool = grav / cp*w_prime;
            float_t ctau = T_prime * (0.004*T_prime - 2.0) + 304.4;
            float_t tau = 1.0/(ctau*cool);
            float_t delta = bcoeff[1] * r_0;
            float_t phi = acoeff[0]*S_i / (acoeff[1]+acoeff[2]*S_i) * (w_prime - w_pre);

            // monodisperse approximation following KHL06
            float_t kappa = 2.0 * bcoeff[0]*bcoeff[1]*tau/((1.0+delta)*(1.0+delta));
            float_t sqrtkap = sqrt(kappa);
            float_t ren = 3.0*sqrtkap / (2.0 + sqrt(1.0+9.0*kappa/M_PI));
            float_t R_imfc = 4.0 * M_PI * bcoeff[0]/(bcoeff[1]*bcoeff[1]) / svol;
            float_t R_im = R_imfc / (1.0+delta) * (delta*delta - 1.0
                + (1.0+0.5*kappa*(1.0+delta)*(1.0+delta)) * ren/sqrtkap);

            // number concentration and radius of ice particles
            float_t ni_hom = phi/R_im;
            float_t ri_0 = 1.0 + 0.5*sqrtkap * ren;
            float_t ri_hom = (ri_0 * (1.0+delta) - 1.0) / bcoeff[1];
            float_t mi_hom = (4.0/3.0 * M_PI * rho_ice) * ni_hom * ri_hom*ri_hom*ri_hom;
            mi_hom = max(mi_hom, cc.ice.min_x_nuc_homo);
            /////// DT PROBLEM
            float_t delta_n;
            // if(cc.dt_prime > 1)
                delta_n = max(min(ni_hom, ni_hom_max/cc.dt_prime), 0.0);
            // else:
            //     delta_n = max(min(ni_hom, ni_hom_max), 0.0);
            float_t delta_q;
            // if(cc.dt_prime > 1)
                delta_q = min(delta_n * mi_hom, qv_prime/cc.dt_prime);
            // else
            //     delta_q = min(delta_n * mi_hom, qv_prime);
            /////// DT PROBLEM SOLVED
            // float_t delta_n = max(min(ni_hom, ni_hom_max), 0.0);
            // float_t delta_q = min(delta_n * mi_hom, qv_prime);

            res[Ni_idx] += delta_n;
            res[qi_idx] += delta_q;
            res[qv_idx] -= delta_q;
#ifdef TRACE_QV
            if(trace)
                std::cout << "KHL06 dqv " << -delta_q << "\n";
#endif
#ifdef TRACE_QI
            if(trace)
                std::cout << "KHL06 dqi " << delta_q << ", dNi " << delta_n << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime) * delta_q / specific_heat_ice(T_prime);
            // Sublimation, cooling
            if(delta_e < 0.0)
                res[lat_cool_idx] += delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] += delta_e;
        }
    }
}


/**
 * Heterogeneous nucleation using Hande et al.
 */
template<class float_t>
void ice_activation_hande(
    float_t &qc_prime,
    float_t &qv_prime,
    float_t &T_prime,
    float_t &S_i,
    float_t &n_inact,
    bool &ndiag_mask,
    std::vector<float_t> &res,
    model_constants_t &cc
)
{
    if(T_prime < T_nuc && T_prime > 180.0 && S_i > 1.0 && n_inact < ni_het_max)
    {

        const double EPSILON = 1.0e-20;
        float_t ndiag = 0.0;
        if(qc_prime > EPSILON)
        {
            float_t T_tmp = max(T_prime, 237.1501);
            if(T_tmp < 261.15)
            {
                ndiag = nim_imm * exp(-alf_imm
                    * exp(bet_imm*log(T_tmp-237.15)));
            }
        } else
        {
            // Hande et al. scheme, Eq. (3) with (2) and (1)
            float_t T_tmp = max(T_prime, 220.001);
            if(T_tmp < 253.0)
            {
                ndiag = nin_dep * exp(-alf_dep
                    * exp(bet_dep*log(T_tmp - 220.0)));
                ndiag = ndiag * (a_dep * atan(b_dep*(S_i-1.0)+c_dep)+d_dep);
            }
        }
        /////// DT PROBLEM
        float_t delta_n = max(ndiag - n_inact, 0.0)/cc.dt_prime;
        float_t delta_q;
        // if(cc.dt_prime > 1)
            delta_q = min(delta_n * cc.ice.min_x_nuc_hetero, qv_prime/cc.dt_prime);
        // else:
        //     delta_q = min(delta_n * cc.ice.min_x_nuc_hetero, qv_prime);
        /////// DT PROBLEM SOLVED
        // float_t delta_n = max(ndiag - n_inact, 0.0);
        // float_t delta_q = min(delta_n * cc.ice.min_x_nuc_hetero, qv_prime);
        delta_n = delta_q/cc.ice.min_x_nuc_hetero;

        res[qi_idx] += delta_q;
        res[Ni_idx] += delta_n;
        res[qv_idx] -= delta_q;
#ifdef TRACE_QV
        if(trace)
            std::cout << "heterogeneous nucleation dqv " << -delta_q << "\n";
#endif
#ifdef TRACE_QI
        if(trace)
            std::cout << "heterogeneous nucleation dqi " << delta_q << ", dNi " << delta_n << "\n";
#endif
        n_inact += delta_n;

        // latent heating and cooling
        float_t delta_e = latent_heat_melt(T_prime) * delta_q / specific_heat_ice(T_prime);
        // Sublimation, cooling
        if(delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] += delta_e;
    }
}


/**
 * Heterogeneous nucleation using Phillips et al. (2008)
 * Implementation by Carmen Koehler and AS
 * modified for C++ and Codipack by Maicon Hieronymus
 */
template<class float_t>
void ice_activation_phillips(
    float_t &qc_prime,
    float_t &qv_prime,
    float_t &T_prime,
    // float_t &p_sat_ice,
    float_t &S_i,
    float_t &n_inact,
    bool &use_prog_in,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    // initial number density of dust/soot/organics [1/m3]
    double na_dust;
    double na_soot;
    double na_orga;
    if(nuc_type == 6)
    {
        na_dust = 160e4;
        na_soot = 30e6;
        na_orga = 0.0;
    } else if(nuc_type == 7)
    {
        // Standard values
        na_dust = 160e4;
        na_soot = 25e6;
        na_orga = 30e6;
    } else if(nuc_type == 8)
    {
        na_dust = 70e4;
        na_soot = 0.0;
        na_orga = 0.0;
    }

    const double EPSILON = 1.0e-20;
#ifdef TRACE_QI
    if(trace)
        std::cout << "T_prime " << T_prime << "\n"
                  << "T_nuc " << T_nuc << "\n"
                  << "S_i " << S_i << "\n"
                  << "n_inact " << n_inact << "\n"
                  << "ni_het_max " << ni_het_max << "\n";
#endif
    if(T_prime < T_nuc && T_prime > 180.0 && S_i > 1.0
        && n_inact < ni_het_max)
    {
        float_t x_t = (274.0 - T_prime) / t_tstep;
        x_t = min(x_t, t_tmax-1);
        int tt = (int) x_t.getValue() - 1;

        std::vector<float_t> infrac(3);
        if(qc_prime > EPSILON)
        {
            // Immersion freezing at water saturation
            infrac[0] = ( trunc(x_t)+1.0-x_t ) * afrac_dust[98][tt]
                + ( x_t-trunc(x_t) ) * afrac_dust[98][tt+1];
            infrac[1] = ( trunc(x_t)+1.0-x_t ) * afrac_soot[98][tt]
                + ( x_t-trunc(x_t) ) * afrac_soot[98][tt+1];
            infrac[2] = ( trunc(x_t)+1.0-x_t ) * afrac_orga[98][tt]
                + ( x_t-trunc(x_t) ) * afrac_orga[98][tt+1];
        } else
        {
            // deposition nucleation below water saturation
            // Indices for 2D look-up tables
            float_t x_s = 100.0*(S_i-1.0) / s_sstep;
            x_s = min(x_s, s_smax-1);
            int ss = std::max(0, (int) (x_s.getValue()-1) );
            float_t S_sr = max(1.0, trunc(x_s));
            infrac[0] = ( trunc(x_t)+1.0-x_t ) * ( S_sr+1.0-x_s )
                * afrac_dust[ss][tt]
                + ( x_t-trunc(x_t) ) * ( S_sr+1.0-x_s )
                * afrac_dust[ss][tt+1]
                + ( x_t-trunc(x_t)+1.0-x_t ) * ( x_s-S_sr )
                * afrac_dust[ss+1][tt]
                + ( x_t-trunc(x_t) ) * ( x_s-S_sr )
                * afrac_dust[ss+1][tt+1];
            infrac[1] = ( trunc(x_t)+1.0-x_t ) * ( S_sr+1.0-x_s )
                * afrac_soot[ss][tt]
                + ( x_t-trunc(x_t) ) * ( S_sr+1.0-x_s )
                * afrac_dust[ss][tt+1]
                + ( x_t-trunc(x_t)+1.0-x_t ) * ( x_s-S_sr )
                * afrac_soot[ss+1][tt]
                + ( x_t-trunc(x_t) ) * ( x_s-S_sr )
                * afrac_soot[ss+1][tt+1];
            infrac[2] = ( trunc(x_t)+1.0-x_t ) * ( S_sr+1.0-x_s )
                * afrac_orga[ss][tt]
                + ( x_t-trunc(x_t) ) * ( S_sr+1.0-x_s )
                * afrac_orga[ss][tt+1]
                + ( x_t-trunc(x_t)+1.0-x_t ) * ( x_s-S_sr )
                * afrac_orga[ss+1][tt]
                + ( x_t-trunc(x_t) ) * ( x_s-S_sr )
                * afrac_dust[ss+1][tt+1];
        }
        float_t ndiag = infrac[0]*na_dust + na_soot*infrac[1] + na_orga*infrac[2];
#ifdef TRACE_QI
        if(trace)
            std::cout << "\n\nqc>EPSILON: " << (qc_prime > EPSILON)
                      << "\ninfac[0]: " << infrac[0]
                      << "\ninfrac[1]: " << infrac[1]
                      << "\ninfrac[2]: " << infrac[2]
                      << "\nna_soot: " << na_soot
                      << "\nna_orga: " << na_orga
                      << "\ntt: " << tt
                      << "\nafrac_dust[98][tt]: " << afrac_dust[98][tt]
                      << "\nqv: " << qv_prime
                      << "\nqc: " << qc_prime
                      << "\nT: " << T_prime
                      << "\nn_inact: " << n_inact
                      << "\nS_i: " << S_i
                      << "\nndiag: " << ndiag << "\nni_het_max: " << ni_het_max
                      << "\nn_inact: " << n_inact << "\nna_dust: " << na_dust << "\n";
#endif
        if(use_prog_in)
        {
            // Not implemented
            // n_inpot replaces n_dust
            // ndiag = n_inpot * ndiag;
            // ndiag_dust = n_inpot*infrac[0];
            // ndiag_all = ndiag;
        } else
        {
            // ndiag = na_dust * ndiag;
        }

        /////// DT PROBLEM
        float_t delta_n;
        float_t delta_q;
        // if(cc.dt_prime > 1)
        // {
            ndiag = min(ndiag, ni_het_max);
            delta_n = max(ndiag-n_inact, 0.0)/cc.dt_prime;
            delta_q = min(delta_n*cc.ice.min_x, qv_prime/cc.dt_prime);
        // } else
        // {
        //     ndiag = min(ndiag, ni_het_max);
        //     delta_n = max(ndiag-n_inact, 0.0);
        //     delta_q = min(delta_n*cc.ice.min_x, qv_prime);
        // }
        /////// DT PROBLEM SOLVED

        delta_n = delta_q/cc.ice.min_x;
        res[Ni_idx] += delta_n;
        res[qi_idx] += delta_q;
        res[qv_idx] -= delta_q;
        res[n_inact_idx] += delta_n;
        res[depo_idx] += delta_q;
#ifdef TRACE_QV
        if(trace)
            std::cout << "Phillips nucleation dqv: " << -delta_q << "\n";
#endif
#ifdef TRACE_QI
        if(trace)
            std::cout << "Phillips nucleation dqi: " << delta_q << ", dNi: " << delta_n << "\n";
        if(trace)
            std::cout << "Phillips nucleation dNi (without dt) "
                << min(max(ndiag-n_inact, 0.0)*cc.ice.min_x, qv_prime/cc.dt_prime) /cc.ice.min_x
                << "\n";

#endif
        // latent heating and cooling
        codi::RealReverse delta_e = latent_heat_melt(T_prime) * delta_q / specific_heat_ice(T_prime);
        // Sublimation, cooling
        if(delta_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] += delta_e;
    }
}

/**
 * Homogeneous freezing of cloud droplets.
 * Freeze only if cloud droplets to freeze are available
 * and the temperature is low enough.
 */
template<class float_t>
void cloud_freeze_hom(
    float_t &qc_prime,
    float_t &Nc,
    float_t &T_prime,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(qc_prime > 0.0 && T_c < -30.0)
    {
        float_t x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x_freezing, cc.cloud.max_x);

        float_t delta_qi;
        float_t delta_ni;
        // instantaneous freezing for temperatures below -50 °C
        if(T_c < -50.0)
        {
            delta_qi = qc_prime;
            delta_ni = Nc;
        } else
        {
            float_t j_hom;
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
            delta_ni = j_hom * qc_prime;
            delta_qi = j_hom * qc_prime * x_c * cc.cloud.c_z;
            /////// DT PROBLEM
            // if(cc.dt_prime > 1)
            // {
                delta_ni = min(delta_ni, Nc/cc.dt_prime);
                delta_qi = min(delta_qi, qc_prime/cc.dt_prime);
            // } else
            // {
            //     delta_ni = min(delta_ni, Nc);
            //     delta_qi = min(delta_qi, qc_prime);
            // }
            /////// DT PROBLEM SOLVED
            // delta_ni = min(delta_ni, Nc);
            // delta_qi = min(delta_qi, qc_prime);
        }
        res[qi_idx] += delta_qi;
        res[Ni_idx] += delta_ni;
        // Remove cloud droplets
        res[qc_idx] -= delta_qi;
        res[Nc_idx] -= delta_ni;
#ifdef TRACE_QC
        if(trace)
            if(abs(delta_qi) > 0)
                std::cout << "cloud freeze dqc " << -delta_qi << ", dNc " << -delta_ni << "\n";
#endif
#ifdef TRACE_QI
        if(trace)
            std::cout << "cloud freeze dqi " << delta_qi << ", dNi " << delta_ni << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime) * delta_qi / specific_heat_ice(T_prime);
        // Melting, cooling
        if(delta_qi < 0.0)
            res[lat_cool_idx] += delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] += delta_e;
    }
}


/**
 *
 */
template<class float_t>
void ice_self_collection(
    float_t &qi_prime,
    float_t &Ni,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    codi::RealReverse x_i = particle_mean_mass(qi_prime, Ni, cc.ice.min_x_collision, cc.ice.max_x);
    codi::RealReverse D_i = particle_diameter(x_i, cc.ice.a_geo, cc.ice.b_geo);
    if(Ni > 0.0 && qi_prime > q_crit_i && D_i > D_crit_i)
    {
        float_t x_conv_i = pow(D_conv_i/cc.snow.a_geo, 1.0/cc.snow.b_geo);
        // efficiency depends on temperature here (Cotton et al 1986)
        // also Straka 1989, page 53
        float_t e_coll = min(pow(10, 0.035*T_c-0.7), 0.2);
        float_t vel_i = particle_velocity(x_i, cc.ice.a_vel, cc.ice.b_vel) * cc.ice.rho_v;

        float_t delta_n = M_PI/4.0 * e_coll * cc.ice.sc_delta_n
            * Ni * Ni * D_i * D_i * sqrt(
                cc.ice.sc_theta_n * vel_i * vel_i + 2.0 * cc.ice.s_vel * cc.ice.s_vel);
        float_t delta_q = M_PI/4.0 * e_coll * cc.ice.sc_delta_q
            * Ni * qi_prime * D_i * D_i * sqrt(
                cc.ice.sc_theta_q * vel_i * vel_i + 2.0 * cc.ice.s_vel * cc.ice.s_vel);

        /////// DT PROBLEM
        if(cc.dt_prime > 1)
        {
            delta_q = min(delta_q, qi_prime/cc.dt_prime);
            delta_n = min(min( delta_n, delta_q/x_conv_i), Ni/cc.dt_prime);
        } else
        {
            delta_q = min(delta_q, qi_prime);
            delta_n = min(min( delta_n, delta_q/x_conv_i), Ni);
        }
        /////// DT PROBLEM SOLVED
        // delta_q = min(delta_q, qi_prime);
        // delta_n = min(min( delta_n, delta_q/x_conv_i), Ni);

        res[qi_idx] -= delta_q;
        res[qs_idx] += delta_q;
        res[Ni_idx] -= delta_n;
        res[Ns_idx] += delta_n/2.0;
#ifdef TRACE_QI
        if(trace)
            std::cout << "ice self colletion dqi " << -delta_q << ", dNi " << -delta_n << "\n";
#endif
#ifdef TRACE_QS
        if(trace)
            std::cout << "ice self colletion dqs " << delta_q << ", dNs " << delta_n/2.0 << "\n";
#endif
    }
}

/**
 *
 */
template<class float_t>
void snow_self_collection(
    float_t &qs_prime,
    float_t &Ns,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(qs_prime > q_crit)
    {
        // temperature dependent sticking efficiency Lin (1983)
        float_t e_coll = max(0.1, min(exp(0.09*(T_prime-T_freeze)), 1.0));
        float_t x_s = particle_mean_mass(qs_prime, Ns, cc.snow.min_x_collection, cc.snow.max_x);
        float_t D_s = particle_diameter(x_s, cc.snow.a_geo, cc.snow.b_geo);
        float_t vel_s = particle_velocity(x_s, cc.snow.a_vel, cc.snow.b_vel) * cc.snow.rho_v;

        res[Ns_idx] -= M_PI/8.0 * e_coll * Ns * Ns * cc.snow.sc_delta_n * D_s * D_s
            * sqrt(cc.snow.sc_theta_n * vel_s * vel_s + 2.0 * cc.snow.s_vel * cc.snow.s_vel);
#ifdef TRACE_QS
        if(trace)
            std::cout << "snow self collection dNs " << - (M_PI/8.0 * e_coll * Ns * Ns * cc.snow.sc_delta_n * D_s * D_s
                * sqrt(cc.snow.sc_theta_n * vel_s * vel_s + 2.0 * cc.snow.s_vel * cc.snow.s_vel)) << "\n";
#endif
    }
}


/**
 *
 */
template<class float_t>
void snow_melting(
    float_t &qs_prime,
    float_t &Ns,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(T_prime > T_freeze && qs_prime > 0.0)
    {
        float_t p_sat = saturation_pressure_water_icon(T_prime);

        float_t x_s = particle_mean_mass(qs_prime, Ns, cc.snow.min_x_melt, cc.snow.max_x);
        float_t d_s = particle_diameter(x_s, cc.snow.a_geo, cc.snow.b_geo);
        float_t v_s = particle_velocity(x_s, cc.snow.a_vel, cc.snow.b_vel) * cc.snow.rho_v;

        float_t fv_q = cc.snow.a_f + cc.snow.b_f * sqrt(v_s*d_s);
        // From Rasmussen and Heymsfield (1987)
        float_t fh_q = 1.05 * fv_q;
        float_t melt = 2.0*M_PI/L_ew * d_s * Ns;
        float_t melt_h = melt * K_T * (T_prime - T_freeze);
        float_t melt_v = melt * dv0 * L_wd/Rv * (p_sat/T_prime - p_sat_melt/T_freeze);
        float_t melt_q = (melt_h * fh_q + melt_v * fv_q);
        /////// DT PROBLEM
        float_t melt_n;
        // if(cc.dt_prime > 1)
        // {
            melt_n = min(max( (melt_q-qs_prime)/x_s + Ns, 0.0), Ns/cc.dt_prime);
            melt_q = min(qs_prime/cc.dt_prime, max(melt_q, 0.0));
            melt_n = min(Ns/cc.dt_prime, max(melt_n, 0.0));
            if(T_prime - T_freeze > 10.0)
            {
                melt_q = qs_prime/cc.dt_prime;
                melt_n = Ns/cc.dt_prime;
            }
        // } else
        // {
        //     melt_n = min(max( (melt_q-qs_prime)/x_s + Ns, 0.0), Ns);
        //     melt_q = min(qs_prime, max(melt_q, 0.0));
        //     melt_n = min(Ns, max(melt_n, 0.0));
        //     if(T_prime - T_freeze > 10.0)
        //     {
        //         melt_q = qs_prime;
        //         melt_n = Ns;
        //     }
        // }
        /////// DT PROBLEM SOLVED
        // float_t melt_n = min(max( (melt_q-qs_prime)/x_s + Ns, 0.0), Ns);
        // melt_q = min(qs_prime, max(melt_q, 0.0));
        // melt_n = min(Ns, max(melt_n, 0.0));

        // Spontaneous melting at 10 C
        // if(T_prime - T_freeze > 10.0)
        // {
        //     melt_q = qs_prime;
        //     melt_n = Ns;
        // }

        // Snow
        res[qs_idx] -= melt_q;
        // Snow N
        res[Ns_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;
#ifdef TRACE_QR
        if(trace)
            std::cout << "snow melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
#ifdef TRACE_QS
        if(trace)
            std::cout << "snow melting dqs " << -melt_q << ", dNs " << -melt_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
        // TODO: Check if we need that
        // snow%n(i,k) = MAX(snow%n(i,k), snow%q(i,k)/snow%x_max)
    }
}


/**
 *
 */
template<class float_t>
void graupel_melting(
    float_t &qg_prime,
    float_t &Ng,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(T_prime > T_freeze && qg_prime > 0.0)
    {
        float_t p_sat = saturation_pressure_water_icon(T_prime);
        float_t x_g = particle_mean_mass(qg_prime, Ng, cc.graupel.min_x_melt, cc.graupel.max_x);
        float_t d_g = particle_diameter(x_g, cc.graupel.a_geo, cc.graupel.b_geo);
        float_t v_g = particle_velocity(x_g, cc.graupel.a_vel, cc.graupel.b_vel) * cc.graupel.rho_v;

        float_t fv_q = cc.graupel.a_f + cc.graupel.b_f * sqrt(v_g*d_g);
        float_t fh_q = 1.05 * fv_q;
        float_t melt = 2.0*M_PI/L_ew * d_g * Ng;
        float_t melt_h = melt * K_T * (T_prime - T_freeze);
        float_t melt_v = melt * dv0 * L_wd/Rv * (p_sat/T_prime - p_sat_melt/T_freeze);
        float_t melt_q = (melt_h * fh_q + melt_v * fv_q);
        /////// DT PROBLEM
        float_t melt_n = min(max( (melt_q-qg_prime)/x_g + Ng, 0.0), Ng/cc.dt_prime);
        melt_q = max(0.0, min(melt_q, qg_prime/cc.dt_prime));
        melt_n = max(0.0, max(melt_n, Ng/cc.dt_prime));
        /////// DT PROBLEM SOLVED
        // float_t melt_n = min(max( (melt_q-qg_prime)/x_g + Ng, 0.0), Ng);
        // melt_q = max(0.0, min(melt_q, qg_prime));
        // melt_n = max(0.0, max(melt_n, Ng));

        // Graupel
        res[qg_idx] -= melt_q;
        // Graupel N
        res[Ng_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;
#ifdef TRACE_QR
        if(trace)
            std::cout << "graupel melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
#ifdef TRACE_QG
        if(trace)
            std::cout << "graupel melting dqg " << -melt_q << ", dNg " << -melt_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 *
 */
template<class float_t>
void hail_melting(
    float_t &qh_prime,
    float_t &Nh,
    float_t &T_prime,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(T_prime > T_freeze && qh_prime > 0.0)
    {
        float_t p_sat = saturation_pressure_water_icon(T_prime);
        float_t x_h = particle_mean_mass(qh_prime, Nh, cc.hail.min_x_melt, cc.hail.max_x);
        float_t d_h = particle_diameter(x_h, cc.hail.a_geo, cc.hail.b_geo);
        float_t v_h = particle_velocity(x_h, cc.hail.a_vel, cc.hail.b_vel) * cc.hail.rho_v;

        float_t fv_q = cc.hail.a_f + cc.hail.b_f * sqrt(v_h*d_h);
        float_t fh_q = 1.05 * fv_q;
        float_t melt = 2.0*M_PI/L_ew * d_h * Nh;
        float_t melt_h = melt * K_T * (T_prime - T_freeze);
        float_t melt_v = melt * dv0 * L_wd/Rv * (p_sat/T_prime - p_sat_melt/T_freeze);
        float_t melt_q = (melt_h * fh_q + melt_v * fv_q);
        /////// DT PROBLEM
        float_t melt_n = min(max( (melt_q-qh_prime)/x_h + Nh, 0.0), Nh/cc.dt_prime);
        melt_q = max(0.0, min(melt_q, qh_prime/cc.dt_prime));
        melt_n = max(0.0, max(melt_n, Nh/cc.dt_prime));
        /////// DT PROBLEM SOLVED
        // float_t melt_n = min(max( (melt_q-qh_prime)/x_h + Nh, 0.0), Nh);
        // melt_q = max(0.0, min(melt_q, qh_prime));
        // melt_n = max(0.0, max(melt_n, Nh));

        // Hail
        res[qh_idx] -= melt_q;
        // Hail N
        res[Nh_idx] -= melt_n;
        // Rain
        res[qr_idx] += melt_q;
        // Rain N
        res[Nr_idx] += melt_n;
    #ifdef TRACE_QR
        if(trace)
            std::cout << "hail melting dqr " << melt_q << ", dNr " << melt_n << "\n";
    #endif
    #ifdef TRACE_QH
        if(trace)
            std::cout << "hail melting dqh " << -melt_q << ", dNh " << -melt_n << "\n";
    #endif
        float_t delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 *
 */
template<class float_t>
void auto_conversion_kb(
    float_t &qc_prime,
    float_t &Nc,
    float_t &qr_prime,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    // autoconversionKB
    float_t k_a = 6.0 + 25 * pow(9.59, -1.7);
    float_t x_s_i = 1.0/cc.cloud.max_x;
    float_t x_c = particle_mean_mass(
        qc_prime, Nc, cc.cloud.min_x_conversion, cc.cloud.max_x);
    // Using Beheng 1994
    float_t au = k_a * pow(x_c*1e3, 3.3) * pow(qc_prime*1e3, 1.4) * 1e3;
    au = min(qc_prime/cc.dt_prime, au);
    res[Nr_idx] += au*x_s_i;
    res[qr_idx] += au;
    res[Nc_idx] -= au*x_s_i*2.0;
    res[qc_idx] -= au;
#ifdef TRACE_QC
    if(trace)
        if(abs(au) > 0)
            std::cout << "autoconversion dqc " << -au << ", dNc " << -au*x_s_i*2.0 << "\n";
#endif
#ifdef TRACE_QR
    if(trace)
        std::cout << "autoconversion dqr " << au << ", dNr " << au*x_s_i << "\n";
#endif
    // accretionKB
    if(qc_prime > q_crit_i && qr_prime > q_crit_i)
    {
        // k_r = 6.0 from Beheng (1994)
        float_t ac = 6.0 * qc_prime * qr_prime;
        ac = min(qc_prime/cc.dt_prime, ac);
        res[qr_idx] += ac;
        res[qc_idx] -= ac;
#ifdef TRACE_QR
        if(trace)
            std::cout << "accretionKB dqr " << ac << "\n";
#endif
#ifdef TRACE_QC
        if(trace)
            if(abs(ac) > 0)
                std::cout << "accretionKB dqc " << -ac << "\n";
#endif
    }
}

/**
 * Formation of raindrops by coagulating cloud droplets and growth
 * of raindrops collecting cloud droplets using Seifert and Beheng (2006)
 * section 2.2.1.
 */
template<class float_t>
void auto_conversion_sb(
    float_t &qc_prime,
    float_t &Nc,
    float_t &qr_prime,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    const double EPSILON = 1.0e-25;
    // autoconversionSB
    if(qc_prime > q_crit)
    {

#ifdef SB_CONV
        const double k_1 = 400;
        const double k_2 = 0.7;
#else
        const double k_1 = 600;
        const double k_2 = 0.68;
#endif
        float_t x_c = particle_mean_mass(
            qc_prime, Nc, cc.cloud.min_x_conversion, cc.cloud.max_x);
        float_t au = cc.cloud_k_au * qc_prime*qc_prime
                                * x_c*x_c * cc.cloud.rho_v;
        float_t tau = min(max(1.0-qc_prime/
                                (qc_prime+qr_prime+EPSILON), EPSILON), 0.9);
        float_t phi = k_1 * pow(tau, k_2) * pow(1.0-pow(tau, k_2), 3);
        au *= (1.0 + phi/pow(1.0-tau, 2));
        /////// DT PROBLEM
        au = max(min(qc_prime/cc.dt_prime, au), 0.0);
        /////// DT PROBLEM SOLVED

        float_t sc = cc.cloud_k_sc * qc_prime*qc_prime * cc.cloud.rho_v;

        res[qr_idx] += au;
        // cloud.max_x = rain.min_x
        res[Nr_idx] += au / cc.cloud.max_x;
        res[qc_idx] -= au;
        /////// DT PROBLEM
        res[Nc_idx] -= min(Nc/cc.dt_prime, sc);
        /////// DT PROBLEM SOLVED
#ifdef TRACE_QC
        if(trace)
            if(abs(au) > 0)
                std::cout << "autoconversion dqc " << -au << ", dNc " << -min(Nc, sc) << "\n";
#endif
#ifdef TRACE_QR
        if(trace)
            std::cout << "autoconversion dqr " << au << ", dNr " << au / cc.cloud.max_x << "\n";
#endif
    }

    // accretionSB
    if(qc_prime > 0.0 && qr_prime > 0.0)
    {
        const double k_1 = 5.0e-4;
        const double k_r = 5.78;
        float_t tau = min(max(1.0-qc_prime/
                                (qc_prime+qr_prime+EPSILON), EPSILON), 1.0);
        float_t phi = pow(tau/(tau+k_1), 4);
        float_t ac = k_r * qc_prime * qr_prime * phi;
        /////// DT PROBLEM
        ac = min(qc_prime/cc.dt_prime, ac);
        /////// DT PROBLEM SOLVED
        float_t x_c = particle_mean_mass(
            qc_prime, Nc, cc.cloud.min_x_conversion, cc.cloud.max_x);
        res[qr_idx] += ac;
        res[qc_idx] -= ac;
        /////// DT PROBLEM
        res[Nc_idx] -= min(Nc/cc.dt_prime, x_c);
        /////// DT PROBLEM SOLVED
#ifdef TRACE_QC
        if(trace)
            std::cout << "accretionSB dqc " << -ac << ", dNc " << -min(Nc, x_c) << "\n";
#endif
#ifdef TRACE_QR
        if(trace)
            std::cout << "accretionSB dqr " << ac << "\n";
#endif
    }
}


/**
 * Rain self collection after Seifert and Beheng (2001)
 */
template<class float_t>
void rain_self_collection_sb(
    float_t &qr_prime,
    float_t &Nr,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(qr_prime > 0)
    {
        float_t x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x_collection, cc.rain.max_x);
        float_t D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);
        // Parameters based on Seifert (2008)
        float_t sc = 4.33 * Nr * qr_prime * cc.rain.rho_v; // rhain%rho_v(i, k)
        // Breakup Seifert (2008), Eq. A13
        float_t breakup = 0.0;
        if(D_r > 0.30e-3)
            breakup = sc * (1.0e+3 * (D_r - 1.10e-3) + 1.0);
        /////// DT PROBLEM
        res[Nr_idx] -= min(Nr/cc.dt_prime, sc-breakup);
        /////// DT PROBLEM SOLVED
#ifdef TRACE_QR
        if(trace)
            std::cout << "self collection dNr " << -min(Nr, sc-breakup) << "\n";
        if(trace)
            std::cout << "Nr " << Nr
                      << "\nsc: " << sc
                      << "\nbreakup: " << breakup
                      << "\nrain.rho_v: " << cc.rain.rho_v << "\n";
#endif
    }
}


/**
 * Rain evaporation after Seifert (2008)
 */
template<class float_t>
void rain_evaporation_sb(
    float_t &qr_prime,
    float_t &Nr,
    float_t &qv_prime,
    float_t &qc_prime,
    float_t &T_prime,
    float_t &p_prime,
    float_t &s_sw,
    float_t &p_sat,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(s_sw < 0.0 && qr_prime > 0.0 && qc_prime < q_crit)
    {
        float_t D_v = diffusivity(T_prime, p_prime);
        // Equation (A2) of 10.1175/2008JAS2586.1
        float_t g_d = 2*M_PI /
            ( (R_v*T_prime)/(D_v*p_sat) + (L_wd*L_wd)/(K_T*R_v*T_prime*T_prime) );
        float_t x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x_evap, cc.rain.max_x);
        float_t D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);

        float_t mue;
        // Equation 20 of Seifert (2008)
        if(D_v <= cc.rain.cmu3)
            mue = cc.rain.cmu0 * tanh( pow(4.0*cc.rain.cmu2*(D_v-cc.rain.cmu3),
                cc.rain.cmu5)) + cc.rain.cmu4;
        else
            mue = cc.rain.cmu1 * tanh( pow(cc.rain.cmu2*(D_v-cc.rain.cmu3),
                cc.rain.cmu5)) + cc.rain.cmu4;
        // Equation A8
        float_t lambda = pow(
            M_PI/6.0*rho_w*(mue+3.0)*(mue+2.0)*(mue+1.0)/x_r, 1.0/3.0);

        // Approximation of Gamma(mue+5/2) / Gamma(mue+2)
        float_t gamma_approx = 0.1357940435E+01
            + mue * (0.3033273220E+00
            + mue * (-0.1299313363E-01
            + mue * (0.4002257774E-03
            - mue * 0.4856703981E-05)));
        // ventilation factor (A7) with (A5) and (A9)
        // Approximation for terminal fall velocity of raindrops
        float_t f_v = a_v + b_v * pow(N_Sc, 1.0/3.0) * gamma_approx
            * sqrt(a_prime/kin_visc_air * cc.rain.rho_v/lambda)
            * (1.0 - 0.5 * (b_prime/a_prime) * pow(lambda/(c_prime+lambda), (mue+2.5))
                - 0.125 * pow(b_prime/a_prime, 2.0) * pow(lambda/(2*c_prime+lambda), (mue+2.5))
                - 1.0/16.0 * pow(b_prime/a_prime, 3.0) * pow(lambda/(3*c_prime+lambda), (mue+2.5))
                - 5.0/127.0 * pow(b_prime/a_prime, 4.0) * pow(lambda/(4*c_prime+lambda), (mue+2.5))
            );
        float_t gamma_eva;
        // D_br = 1.1e-3 from ICON
        if(cc.rain_gfak > 0.0)
            gamma_eva = cc.rain_gfak * (1.1e-3/D_r) * exp(-0.2*mue);
        else
            gamma_eva = 1.0;

        // Equation A5 with A9
        float_t delta_qv = g_d * Nr * (mue+1.0) / lambda * f_v * s_sw; // (mue+1.0) / lambda *
        float_t delta_nv = gamma_eva * delta_qv/x_r;

        /////// DT PROBLEM
        delta_qv = max(-delta_qv, 0.0);
        delta_nv = max(-delta_nv, 0.0);
        delta_qv = min(delta_qv, qv_prime/cc.dt_prime);
        /////// DT PROBLEM SOLVED

        res[qv_idx] += delta_qv;
        res[qr_idx] -= delta_qv;
        res[Nr_idx] -= delta_nv;
        // evap -= delta_qv
#ifdef TRACE_QR
        if(trace)
            std::cout << "rain evaporation after Seifert (2008) dqr " << -delta_qv << ", dNr " << -delta_nv << "\n";
#endif
#ifdef TRACE_QV
        if(trace)
            std::cout << "rain evaporation after Seifert (2008) dqv " << delta_qv << "\n";
#endif
        float_t delta_e = latent_heat_evap(T_prime) * delta_qv
            / specific_heat_water_vapor(T_prime);
        // Evaporation, cooling
        if(delta_qv > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 * Sedimentation after Seifert and Beheng (2006, 2008)
 */
template<class float_t>
void sedimentation_explicit(
    float_t &T_prime,
    float_t &S,
    float_t &qc_prime,
    float_t &qr_prime,
    float_t &Nr,
    float_t &qs_prime,
    float_t &Ns,
    float_t &qi_prime,
    float_t &Ni,
    float_t &qh_prime,
    float_t &Nh,
    float_t &qg_prime,
    float_t &Ng,
    float_t &p_prime,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    // from mo_2mom_mcrph_driver.f90
    // using an explicit flux-form semi-lagrangian scheme (actually used after the microphysics)
    //// sedi_icon_rain


    // using compute_rhoa(p_prime,T_prime,S) instead of cc.rho_a_prime makes
    // next to no difference
    float_t rhocorr = pow(compute_rhoa(p_prime,T_prime,S)/rho_0, -rho_vel);
    float_t v_n_sedi = 0.0;
    float_t v_q_sedi = 0.0;
#ifdef TRACE_QI
    // std::cout << "\nrhocorr: " << rhocorr
    //           << "\nrho_a_prime: " << compute_rhoa(p_prime,T_prime,S)
    //           << "\np_prime: " << p_prime
    //           << "\nT_prime: " << T_prime
    //           << "\nS: " << S
    //           << "\nrho_0: " << rho_0
    //           << "\nrho_vel: " << rho_vel << "\n";
#endif

    auto sedi_icon_core = [&](
        float_t &q,
        float_t &N,
        float_t &v_q_sedi,
        float_t &v_n_sedi,
        float_t &resQ,
        float_t &resN,
        float_t &resQOut,
        float_t &resNOut)
    {
    //    std::cout << "\nCore value resN: " << resN << "\n";
        float_t v_nv = v_n_sedi;
        float_t v_qv = v_q_sedi;  // percentage how much trickles down
        // Assuming v_nv, v_qv is negative
        float_t c_nv = -v_nv * cc.inv_z;//  * dt; // times inverse of layer thickness...
        float_t c_qv = -v_qv * cc.inv_z;//  * dt;

        // TODO: Check why we should multiply with dt here. That *should* be
        // handled outside in rk4.h or at least I thought so.
        float_t s_nv = 0.0;
        if(c_nv <= 1.0)
            s_nv = v_nv*N*cc.inv_z;

        float_t s_qv = 0.0;
        if(c_qv <= 1.0)
            s_qv = v_qv*q*cc.inv_z;

        // abs is used for paranoia reasons and should never be needed
        resN -= abs(s_nv);
        resQ -= abs(s_qv);
        resQOut -= abs(s_qv);
        resNOut -= abs(s_nv);
        // precrate = -abs(s_sq);
#ifdef TRACE_QI
        // std::cout << "\n\ns_qv: " << s_qv << "\nc_qv: " << c_qv
        //     << "\nv_qv: " << v_qv << "\nq: " << q << "\ninv_z: "
        //     << cc.inv_z << "\ns_nv: " << s_nv << "\nresN: " << resN << "\n";
#endif
    };

    auto sedi_icon_sphere = [&](
        float_t &q,
        float_t &N,
        float_t &resQ,
        float_t &resN,
        float_t &resQOut,
        float_t &resNOut,
        particle_model_constants_t &pc)
    {
        float_t v_n_sedi = 0.0;
        float_t v_q_sedi = 0.0;

        if(q > q_crit)
        {
            float_t x = particle_mean_mass(q, N, pc.min_x_sedimentation, pc.max_x);
            float_t lam = pow(pc.lambda*x, pc.b_vel);
            /////// DT PROBLEM
            float_t v_n = max(pc.alfa_n * lam, pc.vsedi_min/cc.dt_prime);
            float_t v_q = max(pc.alfa_q * lam, pc.vsedi_min/cc.dt_prime);

            v_n = min(v_n, pc.vsedi_max/cc.dt_prime);
            v_q = min(v_q, pc.vsedi_max/cc.dt_prime);
            /////// DT PROBLEM SOLVED
            v_n *= rhocorr;
            v_q *= rhocorr;

            v_n_sedi -= v_n;
            v_q_sedi -= v_q;
        }
#ifdef TRACE_QI
        // std::cout << "\nWithin sedi_icon_sphere"
        //       << "\nv_q_sedi: " << v_q_sedi << "\nq: " << q
        //       << "\nv_n_sedi: " << v_n_sedi << "\nN: " << N << "\nresN: " << resN
        //       << "\nrhocorr: " << rhocorr
        //       << "\nvsedi_max: " << pc.vsedi_max
        //       << "\nalfa_q: " << pc.alfa_q
        //     //   << "\nlam: " << pow(pc.lambda*particle_mean_mass(q, N, pc.min_x_sedimentation, pc.max_x), pc.b_vel)
        //       << "\nvsedi_min: " << pc.vsedi_min
        //       << "\nlambda: " << pc.lambda
        //       << "\nb_vel: " << pc.b_vel << "\n";
        //     //   << "\nx: " << particle_mean_mass(q, N, pc.min_x_sedimentation, pc.max_x) << "\n";
#endif
        float_t sedi_q = 0;
        float_t sedi_n = 0;
        sedi_icon_core(q, N, v_q_sedi, v_n_sedi, sedi_q, sedi_n, resQOut, resNOut);

        // TODO: Why did I have those three lines?
        // sedi_n = min(max(abs(sedi_n), abs(sedi_q)/pc.max_x), q/pc.min_x_sedimentation);
        // resN -= sedi_n;
        // resQ -= abs(sedi_q);

        // TODO: This breaks everything I think
        // resN gets positive for one even though it should be negative
        // It results in very high values
        // resN = max(min(N, q/pc.min_x_sedimentation), q/pc.max_x);
        // This should be the correct version
        // Do not throw away more particles than ice mass permits aka
        // only bigger particles shall fall down
        // max(.., q/min_x) ?
#ifdef TRACE_QI
        // std::cout << "\nDifferent res"
        //           << "\nresQ: " << resQ
        //           << "\nresN: " << resN
        //           << "\nsedi_n: " << -sedi_n
        //           << "\nsedi_q: " << -sedi_q << "\n";
        //         //   << "\n-max(min(abs(resN), resQ/pc.min_x_sedimentation), resQ/pc.max_x): " << -max(min(abs(resN), abs(resQ)/pc.min_x_sedimentation), abs(resQ)/pc.max_x)
        //         //   << "\n-max(min(abs(resN), q/pc.min_x_sedimentation), q/pc.max_x): " << -max(min(abs(resN), q/pc.min_x_sedimentation), q/pc.max_x)
        //         //   << "\nresQ/pc.max_x: " << resQ/pc.max_x
        //         //   << "\nresQ/pc.min_x: " << resQ/pc.min_x_sedimentation << "\n";
#endif
        // resN = -min(max(abs(resN), abs(resQ)/pc.max_x), q/pc.min_x_sedimentation);
        // resN = -min(max(abs(resN), resQ/pc.max_x), resQ/pc.min_x_sedimentation);

        // This is what ICON does, but it does not make sense.
        // resN = -max(min(abs(resN), resQ/pc.min_x_sedimentation), resQ/pc.max_x);
        // This might work`?
        // resN = -max(min(abs(resN), q/pc.min_x_sedimentation), q/pc.max_x);
    };

    auto sedi_icon_sphere_lwf = [&](
        )
    {

    };

    if(qr_prime > q_crit)
    {
        float_t x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x_sedimentation, cc.rain.max_x);
        float_t D_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);
        float_t mue = (qc_prime >= q_crit)
            ? (cc.rain.nu+1.0)/cc.rain.b_geo - 1.0
            : rain_mue_dm_relation(D_r, cc.rain.cmu0, cc.rain.cmu1,
                cc.rain.cmu2, cc.rain.cmu3, cc.rain.cmu4);
        // inverse of lambda in Eq. (A10) SB (2008)
        float_t D_p = D_r * pow( (mue+3)*(mue+2)*(mue+1), -1.0/3.0 );

        // SB (2008), Eq. (A10)
        float_t v_nr = cc.rain.alpha - cc.rain.beta
            * pow( 1.0+cc.rain.gamma*D_p, -mue-1.0);
        float_t v_qr = cc.rain.alpha - cc.rain.beta
            * pow( 1.0+cc.rain.gamma*D_p, -mue-4.0);

        // Seifert (2008), Eq. A10
        // rhocorr: height dependency of particle fall speed
        v_nr *= rhocorr;
        v_qr *= rhocorr;

        v_n_sedi -= v_nr;
        v_q_sedi -= v_qr;
    }
#ifdef TRACE_QR
    auto beforeN = res[Nr_idx];
    auto beforeq = res[qr_idx];
#endif
    sedi_icon_core(qr_prime, Nr, v_q_sedi, v_n_sedi, res[qr_idx], res[Nr_idx], res[qr_out_idx], res[Nr_out_idx]);
    // Just for fun and giggles
    // res[qr_idx] -= res[qr_out_idx]*2/3;
    // res[Nr_idx] -= res[Nr_out_idx]*2/3;
    // res[qr_out_idx] /= 3;
    // res[Nr_out_idx] /= 3;
#ifdef TRACE_QR
    if(trace)
        std::cout << "sedi_icon_core dqr " << res[qr_idx] - beforeq << ", dNr " << res[Nr_idx] - beforeN << "\n";
#endif

    // sedi_icon_sphere ice
    if(qi_prime > 0.0)
    {
#ifdef TRACE_QI
        auto before_q = res[qi_idx];
        auto before_n = res[Ni_idx];
#endif
        sedi_icon_sphere(qi_prime, Ni, res[qi_idx], res[Ni_idx], res[qi_out_idx], res[Ni_out_idx], cc.ice);
#ifdef TRACE_QI
        if(trace)
        {
            std::cout << "N before: " << before_n << " and after " << res[Ni_idx] << "\n";
            std::cout << "sedi icon sphere dqi " << res[qi_idx] - before_q << ", dNi " << res[Ni_idx] - before_n << "\n";
        }
#endif
    }

    // sedi_icon_sphere snow
    if(qs_prime > 0.0)
    {
#ifdef TRACE_QS
        auto before_q = res[qs_idx];
        auto before_n = res[Ns_idx];
#endif
        sedi_icon_sphere(qs_prime, Ns, res[qs_idx], res[Ns_idx], res[qs_out_idx], res[Ns_out_idx], cc.snow);
#ifdef TRACE_QS
        if(trace)
            std::cout << "sedi icon sphere dqs " << res[qs_idx] - before_q << ", dNs " << res[Ns_idx] - before_n << "\n";
#endif
    }
    // sedi_icon_sphere graupel
    const bool lprogmelt = false; // TODO: Check if true or false
    if(qg_prime > 0.0)
    {
        if(lprogmelt)
            sedi_icon_sphere_lwf();
        else
        {
#ifdef TRACE_QG

            auto before_q = res[qg_idx];
            auto before_n = res[Ng_idx];
#endif
            sedi_icon_sphere(qg_prime, Ng, res[qg_idx], res[Ng_idx], res[qg_out_idx], res[Ng_out_idx], cc.graupel);
#ifdef TRACE_QG

            // auto differ = (res[qg_idx].getValue() - before_q.getValue());
            // differ = differ/1.5;
            // res[qg_idx] = res[qg_idx]+differ;
            // differ = (res[Ng_idx].getValue() - before_n.getValue());
            // differ = differ/1.5;
            // res[Ng_idx] = res[Ng_idx]+differ;
            sediment_q = sediment_q + res[qg_idx].getValue() - before_q.getValue();
            sediment_n = sediment_n + res[Ng_idx].getValue() - before_n.getValue();
            if(trace)
                std::cout << "sedi icon sphere dqg " << res[qg_idx] - before_q
                      << ", dNg " << res[Ng_idx] - before_n << "\n";
#endif
        }
    }

    // sedi_icon_sphere hail
    if(qh_prime > 0.0)
    {
        if(lprogmelt)
            sedi_icon_sphere_lwf();
        else
        {
#ifdef TRACE_QH
            auto before_q = res[qh_idx];
            auto before_n = res[Nh_idx];
#endif
            sedi_icon_sphere(qh_prime, Nh, res[qh_idx], res[Nh_idx], res[qh_out_idx], res[Nh_out_idx], cc.hail);
#ifdef TRACE_QH
            if(trace)
                std::cout << "sedi icon sphere dqh " << res[qh_idx] - before_q << ", dNh " << res[Nh_idx] - before_n << "\n";
#endif
        }
    }
}


/**
 * Evaporation from melting ice particles.
 */
template<class float_t>
void evaporation(
    float_t &qv_prime,
    float_t &e_d,
    float_t &p_sat,
    float_t &s_sw,
    float_t &T_prime,
    float_t &q1,
    float_t &N1,
    float_t &resq,
    particle_model_constants_t &pc1,
    std::vector<float_t> &res,
    const double &dt_prime)
{
    if(q1 > 0.0 && T_prime > T_freeze)
    {
        // TODO: Check D_v
        float_t g_d = 4.0*M_PI / (L_wd*L_wd
            / (K_T * Rv * T_freeze*T_freeze*T_freeze) + Rv*T_freeze /(D_v * p_sat));

        float_t x_1 = particle_mean_mass(q1, N1, pc1.min_x_evap, pc1.max_x);
        float_t d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);
        float_t v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;

        float_t f_v = pc1.a_f + pc1.b_f * sqrt(v_1*d_1);

        float_t delta_q = g_d * N1 * pc1.c_s * d_1 * f_v * s_sw;
        /////// DT PROBLEM
        delta_q = min(q1/dt_prime, max(-delta_q, 0.0));
        /////// DT PROBLEM SOLVED

        // Vapor
        res[qv_idx] += delta_q;
#ifdef TRACE_QV
        if(trace)
            std::cout << "Evaporation from melting ice particles dqv " << delta_q << "\n";
#endif
        resq -= delta_q;
        // subl? Actually not used here but for the broader simulation
        // We could use it later if we want to see how much evaporates
        // sublq -= delta_q;

        float_t delta_e = latent_heat_ice(T_prime) * delta_q / latent_heat_melt(T_prime);
        // Sublimination, cooling
        if(delta_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}


/**
 * Depositional growth of all ice particles.
 * Deposition and sublimation, where deposition rate of ice and snow are
 * being stored.
 */
template<class float_t>
void vapor_dep_relaxation(
    float_t &qv_prime,
    float_t &qi_prime,
    float_t &Ni,
    float_t &qs_prime,
    float_t &Ns,
    float_t &qg_prime,
    float_t &Ng,
    float_t &qh_prime,
    float_t &Nh,
    float_t &s_si,
    float_t &p_sat_ice,
    float_t &T_prime,
    const double &EPSILON,
    float_t &dep_rate_ice,
    float_t &dep_rate_snow,
    float_t &D_vtp,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(T_prime < T_freeze)
    {
        float_t dep_ice = 0.0;
        float_t dep_snow = 0.0;
        float_t dep_graupel = 0.0;
        float_t dep_hail = 0.0;
        float_t g_i = 4.0*M_PI / (L_ed*L_ed / (K_T*Rv*T_prime*T_prime) + Rv*T_prime / (D_vtp*p_sat_ice));

        auto vapor_deposition = [&](
            float_t &q,
            float_t &N,
            particle_model_constants_t &pc,
            float_t &dep)
        {
            if(q <= 0.0)
            {
                dep = 0.0;
            } else
            {
                float_t x = particle_mean_mass(q, N, pc.min_x_depo, pc.max_x);
                float_t d = particle_diameter(x, pc.a_geo, pc.b_geo);
                float_t v = particle_velocity(x, pc.a_vel, pc.b_vel) * pc.rho_v;
                float_t f_v = pc.a_f + pc.b_f * sqrt(d*v);
                f_v = max(f_v, pc.a_f/pc.a_ven);
                dep = g_i * N * pc.c_s * d * f_v * s_si;
            }
        };

        vapor_deposition(qi_prime, Ni, cc.ice, dep_ice);
        vapor_deposition(qs_prime, Ns, cc.snow, dep_snow);
        vapor_deposition(qg_prime, Ng, cc.graupel, dep_graupel);
        vapor_deposition(qh_prime, Nh, cc.hail, dep_hail);

        // Depositional growth based on
        // "A New Double-Moment Microphysics Parameterization for Application in Cloud and
        // Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov
        float_t qvsidiff = qv_prime - p_sat_ice /(Rv*T_prime);
        if(abs(qvsidiff) > EPSILON)
        {
            float_t tau_i_i = 1.0/qvsidiff*dep_ice;
            float_t tau_s_i = 1.0/qvsidiff*dep_snow;
            float_t tau_g_i = 1.0/qvsidiff*dep_graupel;
            float_t tau_h_i = 1.0/qvsidiff*dep_hail;

            float_t xi_i = tau_i_i + tau_s_i + tau_g_i + tau_h_i;
            // TODO: Check wether dt is needed here or not
            float_t xfac = (xi_i < EPSILON) ?
                (float_t) 0.0 : qvsidiff/xi_i * (1.0-exp(-xi_i));

            dep_ice     = xfac * tau_i_i;
            dep_snow    = xfac * tau_s_i;
            dep_graupel = xfac * tau_g_i;
            dep_hail    = xfac * tau_h_i;

            // Is that even necessary?
            if(qvsidiff < 0.0)
            {
                /////// DT PROBLEM
                dep_ice     = max(dep_ice,      -qi_prime/cc.dt_prime);
                dep_snow    = max(dep_snow,     -qs_prime/cc.dt_prime);
                dep_graupel = max(dep_graupel,  -qg_prime/cc.dt_prime);
                dep_hail    = max(dep_hail,     -qh_prime/cc.dt_prime);
                /////// DT PROBLEM SOLVED
            }
            // DEBUG
            // dep_graupel = abs(dep_graupel)*1e6;
            // dep_hail = abs(dep_hail);
            // dep_snow = abs(dep_snow);
            // dep_ice = abs(dep_ice);
            float_t dep_sum = dep_ice + dep_graupel + dep_snow + dep_hail;

            res[qi_idx] += dep_ice;
            res[qs_idx] += dep_snow;
            res[qg_idx] += dep_graupel;
            res[qh_idx] += dep_hail;
            res[qv_idx] -= dep_sum;
#ifdef TRACE_QI
            if(trace)
            {
                std::cout << "Depos growth dqi " << dep_ice << "\n";
                std::cout << "qvsidiff " << qvsidiff << "\n";
                std::cout << "qi_prime " << qi_prime << "\n";
            }
#endif
#ifdef TRACE_QS
            if(trace)
                std::cout << "Depos growth dqs " << dep_snow << "\n";
#endif
#ifdef TRACE_QG
            if(trace)
            {
                std::cout << "Depos growth dqg " << dep_graupel << "\n";
                std::cout << "also dqv " << -dep_sum << "\n";
                std::cout << "xfac " << xfac
                          << "\ntau_g_i " << tau_g_i
                          << "\n-qi_prime " << -qi_prime << "\n";
            }
#endif
#ifdef TRACE_QH
            if(trace)
                std::cout << "Depos growth dqh " << dep_hail << "\n";
#endif
#ifdef TRACE_QV
            if(trace)
                std::cout << "Depos growth dqv " << -dep_sum << "\n";
#endif

            dep_rate_ice += dep_ice;
            dep_rate_snow += dep_snow;

            if(dep_ice > 0)
            {
                res[depo_idx] += dep_ice;
            } else
            {
                res[sub_idx] += dep_ice;
            }
            if(dep_snow > 0)
            {
                res[depo_idx] += dep_snow;
            } else
            {
                res[sub_idx] += dep_snow;
            }
            if(dep_graupel > 0)
            {
                res[depo_idx] += dep_graupel;
            } else
            {
                res[sub_idx] += dep_graupel;
            }
            if(dep_hail > 0)
            {
                res[depo_idx] += dep_hail;
            } else
            {
                res[sub_idx] += dep_hail;
            }
            float_t delta_e = latent_heat_melt(T_prime) * dep_sum / specific_heat_ice(T_prime);
            // Sublimation, cooling
            if(dep_sum > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Deposition, heating
            else
                res[lat_heat_idx] -= delta_e;
        }
    }
}


/**
 * Particle collection for ice particles such as:
 * graupel+ice  -> graupel
 * graupel+snow -> graupel
 * hail+ice  -> hail
 * hail+snow -> hail
 * snow+ice  -> snow
 * This function does only one of these.
 *
 * @return Vector of two float_t with the change in particle number (0)
 * and change in mass mixing ratio (1)
 */
template<class float_t>
std::vector<float_t> particle_collection(
    float_t &q1,
    float_t &q2,
    float_t &N1,
    float_t &N2,
    float_t &T_c,
    collection_model_constants_t &coeffs,
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    const double &dt_prime)
{
    float_t e_coll = min(exp(0.09*T_c), 1.0);
    float_t x_1 = particle_mean_mass(q1, N1, pc1.min_x_collection, pc1.max_x);
    float_t d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);
    float_t v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;

    float_t x_2 = particle_mean_mass(q2, N2, pc2.min_x_collection, pc2.max_x);
    float_t d_2 = particle_diameter(x_2, pc2.a_geo, pc2.b_geo);
    float_t v_2 = particle_velocity(x_2, pc2.a_vel, pc2.b_vel) * pc2.rho_v;

    float_t coll_n = M_PI/4 * N1 * N2 * e_coll
        * (coeffs.delta_n_aa * d_2 * d_2
            + coeffs.delta_n_ab * d_2 * d_1
            + coeffs.delta_n_bb * d_1 * d_1)
        * sqrt(coeffs.theta_n_aa * v_2 * v_2
            - coeffs.theta_n_ab * v_2 * v_1
            + coeffs.theta_n_bb * v_1 * v_1
            + pc1.s_vel * pc1.s_vel);
    float_t coll_q = M_PI/4 * q1 * N2 * e_coll
        * (coeffs.delta_q_aa * d_2 * d_2
            + coeffs.delta_q_ab * d_2 * d_1
            + coeffs.delta_q_bb * d_1 * d_1)
        * sqrt(coeffs.theta_q_aa * v_2 * v_2
            - coeffs.theta_q_ab * v_2 * v_1
            + coeffs.theta_q_bb * v_1 * v_1
            + pc1.s_vel * pc1.s_vel);
    /////// DT PROBLEM
    coll_n = min(N1/dt_prime, coll_n);
    coll_q = min(q1/dt_prime, coll_q);
    /////// DT PROBLEM SOLVED
    std::vector<float_t> r(2);
    r[0] = coll_n;
    r[1] = coll_q;
    return r;
}


/**
 * Ice particle collections for the following processes:
 * snow+ice         -> snow
 * graupel+graupel  -> graupel
 * graupel+ice      -> graupel
 * graupel+snow     -> graupel
 * All those processes are done in this function.
 * Technically, this could be used for hail instead of graupel as well.
 */
template<class float_t>
void particle_particle_collection(
    float_t &qi_prime,
    float_t &Ni,
    float_t &qs_prime,
    float_t &Ns,
    float_t &qg_prime,
    float_t &Ng,
    float_t &T_prime,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    //// particle_collection snow
    if(qi_prime > q_crit && qs_prime > q_crit)
    {
        std::vector<float_t> delta = particle_collection(
            qi_prime, qs_prime, Ni, Ns, T_c, cc.coeffs_sic, cc.ice, cc.snow, cc.dt_prime);
        res[qs_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
#ifdef TRACE_QI
        if(trace)
            std::cout << "ice-snow collision dqi " << -delta[1] << ", dNi " << -delta[0] << "\n";
#endif
#ifdef TRACE_QS
        if(trace)
            std::cout << "ice-snow collision dqs " << delta[1] << "\n";
#endif
    }

    //// graupel self collection
    if(qg_prime > q_crit)
    {
        float_t x_g = particle_mean_mass(qg_prime, Ng,
            cc.graupel.min_x_collection, cc.graupel.max_x);
        float_t d_g = particle_diameter(x_g,
            cc.graupel.a_geo, cc.graupel.b_geo);
        float_t v_g = particle_velocity(x_g,
            cc.graupel.a_vel, cc.graupel.b_vel);
        float_t delta_n = cc.graupel.sc_coll_n
            * Ng * Ng * d_g * d_g * v_g;
        // sticking efficiency does only distinguish dry and wet
        delta_n *= (T_prime > T_freeze) ? ecoll_gg_wet : ecoll_gg;
        /////// DT PROBLEM
        delta_n = min(delta_n, Ng/cc.dt_prime);
        /////// DT PROBLEM SOLVED
        res[Ng_idx] -= delta_n;
#ifdef TRACE_QG
        if(trace)
            std::cout << "graupel self collection dNg " << -delta_n << "\n";
#endif
    }

    // particle particle collection
    // ice and graupel collision
    if(qi_prime > q_crit && qg_prime > q_crit)
    {
        std::vector<float_t> delta = particle_collection(
            qi_prime, qg_prime, Ni, Ng, T_c, cc.coeffs_gic, cc.ice, cc.graupel, cc.dt_prime);
        res[qg_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
#ifdef TRACE_QG
        if(trace)
            std::cout << "ice-graupel collision dqg " << delta[1] << "\n";
#endif
#ifdef TRACE_QI
        if(trace)
            std::cout << "ice-graupel collision dqi " << -delta[1] << ", dNi " << -delta[0] << "\n";
#endif
    }

    // particle particle collection
    // snow and graupel collision
    if(qs_prime > q_crit && qg_prime > q_crit)
    {
        std::vector<float_t> delta = particle_collection(
            qs_prime, qg_prime, Ns, Ng, T_c, cc.coeffs_gsc, cc.snow, cc.graupel, cc.dt_prime);
        res[qg_idx] += delta[1];
        res[qs_idx] -= delta[1];
        res[Ns_idx] -= delta[0];
#ifdef TRACE_QG
        if(trace)
            std::cout << "snow-graupel collision dqg " << delta[1] << "\n";
#endif
#ifdef TRACE_QS
        if(trace)
            std::cout << "snow-graupel collision dqs " << -delta[1] << ", dNs " << -delta[0] << "\n";
#endif
    }
}


/**
 * Conversion graupel to hail and hail collisions.
 */
template<class float_t>
void graupel_hail_conv(
    float_t &qc_prime,
    float_t &qr_prime,
    float_t &qi_prime,
    float_t &qg_prime,
    float_t &Ng,
    float_t &qh_prime,
    float_t &Nh,
    float_t &p_prime,
    float_t &T_prime,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    float_t x_g = particle_mean_mass(qg_prime, Ng,
            cc.graupel.min_x_conversion, cc.graupel.max_x);
    float_t d_g = particle_diameter(x_g,
        cc.graupel.a_geo, cc.graupel.b_geo);
    float_t Ng_tmp = qg_prime/x_g;
    // supercooled liquid water = rain + cloud water
    float_t qwa_prime = qr_prime + qc_prime;

    if(qwa_prime > 1e-3 && T_c < 0 && qg_prime > q_crit_c)
    {
        // float_t qis = qi + qs; // TODO: Check where this comes from

        float_t d_sep = wet_growth_diam(p_prime, T_prime, qwa_prime,
            qi_prime, ltabdminwgg);

        if(d_sep > 0.0 && d_sep < 10.0*d_g)
        {
            float_t xmin = pow(d_sep/cc.graupel.a_geo, 1.0/cc.graupel.b_geo);
            float_t lam = pow(cc.graupel.g2/(cc.graupel.g1*x_g), cc.graupel.mu);
            xmin = pow(xmin, cc.graupel.mu);
            float_t n_0 = cc.graupel.mu * Ng * pow(lam, cc.graupel.nm1)/cc.graupel.g1;
            float_t lam_xmin = lam*xmin;

            float_t conv_n = n_0 / (cc.graupel.mu * pow(lam, cc.graupel.nm1))
                * table_g1.look_up(lam_xmin);
            float_t conv_q = n_0 / (cc.graupel.mu * pow(lam, cc.graupel.nm2))
                * table_g2.look_up(lam_xmin);
            /////// DT PROBLEM
            conv_n = min(conv_n, Ng/cc.dt_prime);
            conv_q = min(conv_q, qg_prime/cc.dt_prime);
            /////// DT PROBLEM SOLVED
            // Graupel q, n
            res[qg_idx] -= conv_q;
            res[Ng_idx] -= conv_n;
            // Hail q, n
            res[qh_idx] += conv_q;
            res[Nh_idx] += conv_n;

#ifdef TRACE_QG
            if(trace)
                std::cout << "conversion graupel->hail dqg " << - conv_q << ", dNq " << - conv_n << "\n";
#endif
#ifdef TRACE_QH
            if(trace)
                std::cout << "conversion graupel->hail dqh " <<  conv_q << ", dNh " << conv_n << "\n";
#endif
        }
    }
}


/**
 * Hail collision with ice and snow.
 */
template<class float_t>
void hail_collision(
    float_t &qh_prime,
    float_t &Nh,
    float_t &qs_prime,
    float_t &Ns,
    float_t &qi_prime,
    float_t &Ni,
    float_t &T_c,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    // ice and hail
    if(qi_prime > q_crit && qh_prime > q_crit)
    {
        std::vector<float_t> delta = particle_collection(
            qi_prime, qh_prime, Ni, Nh, T_c, cc.coeffs_hic, cc.ice, cc.hail, cc.dt_prime);
        res[qh_idx] += delta[1];
        res[qi_idx] -= delta[1];
        res[Ni_idx] -= delta[0];
#ifdef TRACE_QH
        if(trace)
            std::cout << "ice-hail collision dqh " << delta[1] << "\n";
#endif
#ifdef TRACE_QI
        if(trace)
            std::cout << "ice-hail collision dqi " << -delta[1] << ", dNi " << -delta[0] << "\n";
#endif
    }

    // snow and hail collision
    if(qs_prime > q_crit && qh_prime > q_crit)
    {
        std::vector<float_t> delta = particle_collection(
            qs_prime, qh_prime, Ns, Nh, T_c, cc.coeffs_hsc, cc.snow, cc.hail, cc.dt_prime);
        res[qh_idx] += delta[1];
        res[qs_idx] -= delta[1];
        res[Ns_idx] -= delta[0];
#ifdef TRACE_QH
        if(trace)
            std::cout << "snow-hail collision dqh " << delta[1] << "\n";
#endif
#ifdef TRACE_QS
        if(trace)
            std::cout << "snow-hail collision dqs " << -delta[1] << ", dNs " << -delta[0] << "\n";
#endif
    }
}


/**
 * Rate of ice or snow collecting cloud droplets where values are hard
 * coded for *rain* droplets according to COSMO comments.
 */
template<class float_t>
void riming_cloud_core(
    float_t &qc_prime,
    float_t &Nc,
    float_t &q1,
    float_t &N1,
    particle_model_constants_t &pc1,
    collection_model_constants_t &coeffs,
    float_t &rime_rate_qb,
    float_t &rime_rate_nb,
    model_constants_t &cc)
{
    float_t x_1 = particle_mean_mass(q1, N1, pc1.min_x_riming, pc1.max_x);
    float_t d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);
    float_t x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x_riming, cc.cloud.max_x);
    float_t d_c = particle_diameter(x_c, cc.cloud.a_geo, cc.cloud.b_geo);

    float_t const1 = const0 * pc1.ecoll_c;

    if(qc_prime > q_crit_c && q1 > pc1.q_crit_c
        && d_c > D_crit_c && d_1 > pc1.d_crit_c)
    {
        float_t v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;
        float_t v_c = particle_velocity(x_c, cc.cloud.a_vel, cc.cloud.b_vel) * cc.cloud.rho_v;
        float_t tmp = const1*(d_c - D_crit_c);
        /////// DT PROBLEM
        float_t e_coll = min(pc1.ecoll_c/cc.dt_prime, max(tmp, ecoll_min/cc.dt_prime));
        /////// DT PROBLEM SOLVED

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
}


/**
 * Riming of rain droplets with ice or snow particles.
 */
template<class float_t>
void riming_rain_core(
    float_t &qr_prime,
    float_t &Nr,
    float_t &q1,
    float_t &N1,
    particle_model_constants_t &pc1,
    collection_model_constants_t &coeffs,
    float_t &rime_rate_qa,
    float_t &rime_rate_qb,
    float_t &rime_rate_nb,
    model_constants_t &cc)
{
    float_t x_1 = particle_mean_mass(q1, N1, pc1.min_x_riming, pc1.max_x);
    float_t d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);

    if(qr_prime > q_crit && q1 > q_crit_r && d_1 > D_crit_r)
    {
        float_t x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x_riming, cc.rain.max_x);
        float_t d_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);

        float_t v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;
        float_t v_r = particle_velocity(x_r, cc.rain.a_vel, cc.rain.b_vel) * cc.rain.rho_v;

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

        rime_rate_nb = M_PI/4.0 * N1 * Nr
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
}


/**
 *
 */
template<class float_t>
void ice_riming(
    float_t &qc_prime,
    float_t &Nc,
    float_t &qr_prime,
    float_t &Nr,
    float_t &qi_prime,
    float_t &Ni,
    float_t &dep_rate_ice,
    float_t &rime_rate_qc,
    float_t &rime_rate_nc,
    float_t &rime_rate_qr,
    float_t &rime_rate_nr,
    float_t &rime_rate_qi,
    float_t &T_prime,
    const double &dt,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(dep_rate_ice > 0.0 && dep_rate_ice >= rime_rate_qc+rime_rate_qr)
    {
        // Depositional growth is stronger than riming growth, therefore ice stays ice
        // ice cloud riming
        if(rime_rate_qc > 0.0)
        {
            /////// DT PROBLEM
            float_t rime_q = min(qc_prime/cc.dt_prime, rime_rate_qc);
            float_t rime_n = min(Nc/cc.dt_prime, rime_rate_nc);
            /////// DT PROBLEM SOLVED
            // Ice
            res[qi_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;

            // Cloud N
            res[Nc_idx] -= rime_n;

#ifdef TRACE_QC
            if(trace)
                std::cout << "ice riming dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
#ifdef TRACE_QI
            if(trace)
                std::cout << "ice riming dqi " << rime_q << "\n";
#endif
            if(T_prime < T_freeze && ice_multiplication)
            {
                float_t mult_1 = (T_prime - T_mult_min)*const3;
                float_t mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                // Ice N
                res[Ni_idx] += C_mult * mult_1 * mult_2 * rime_q;
#ifdef TRACE_QI
                if(trace)
                    std::cout << "ice riming with mult dNi " <<  C_mult * mult_1 * mult_2 * rime_q << "\n";
#endif
            }
        }
        // ice rain riming
        if(rime_rate_qr > 0.0)
        {
            /////// DT PROBLEM
            float_t rime_q = min(rime_rate_qr, qr_prime/cc.dt_prime);
            float_t rime_n = min(Nr/cc.dt_prime, rime_rate_nr);
            /////// DT PROBLEM SOLVED
            // Ice
            res[qi_idx] += rime_q;
            // Rain
            res[qr_idx] -= rime_q;
            // Rain N
            res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
            if(trace)
                std::cout << "ice rain riming dqr " << -rime_q << ", dNr " << -rime_n << "\n";
#endif
#ifdef TRACE_QI
            if(trace)
                std::cout << "ice rain riming A dqi " << rime_q << "\n";
#endif
            if(T_prime < T_freeze && ice_multiplication)
            {
                float_t mult_1 = (T_prime - T_mult_min)*const3;
                float_t mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                // Ice N
                res[Ni_idx] += C_mult * mult_1 * mult_2 * rime_q;
#ifdef TRACE_QI
                if(trace)
                    std::cout << "ice rain rimingwith mult dNi " << C_mult * mult_1 * mult_2 * rime_q << "\n";
#endif
            }

            float_t delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
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
            float_t x_i = particle_mean_mass(qi_prime, Ni,
                cc.ice.min_x_riming, cc.ice.max_x);
            float_t d_i = particle_diameter(x_i,
                cc.ice.a_geo, cc.ice.b_geo);
            /////// DT PROBLEM
            float_t rime_q = min(rime_rate_qc, qc_prime/cc.dt_prime);
            float_t rime_n = min(rime_rate_nc, Nc/cc.dt_prime);
            /////// DT PROBLEM
            // Ice
            res[qi_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;

            // Cloud N
            res[Nc_idx] -= rime_n;
#ifdef TRACE_QC
            if(trace)
                if(abs(rime_q) > 0)
                    std::cout << "ice cloud riming dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
#ifdef TRACE_QI
            if(trace)
                std::cout << "ice cloud riming dqi " << rime_q << "\n";
#endif
            if(T_prime < T_freeze && ice_multiplication)
            {
                float_t mult_1 = (T_prime - T_mult_min)*const3;
                float_t mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                // Ice N
                res[Ni_idx] += C_mult * mult_1 * mult_2 * rime_q;
#ifdef TRACE_QI
                if(trace)
                    std::cout << "ice cloud rimingwith mult dNi " << C_mult * mult_1 * mult_2 * rime_q << "\n";
#endif
            }

            float_t delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
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
                float_t conv_q = rime_q
                    / (const5*(M_PI/6.0 * rho_ice * d_i*d_i*d_i/x_i -1.0));
                conv_q = min(qi_prime, conv_q);
                float_t qi_tmp = qi_prime+dt*res[qi_idx];
                x_i = particle_mean_mass(qi_tmp, Ni,
                    cc.ice.min_x_conversion, cc.ice.max_x);
                float_t tmp = conv_q / max(x_i, x_conv);
                /////// DT PROBLEM
                float_t conv_n = min(tmp, Ni/cc.dt_prime);
                /////// DT PROBLEM

                // Ice
                res[qi_idx] -= conv_q;
                // Graupel
                res[qg_idx] += conv_q;
                // Ice N
                res[Ni_idx] -= conv_n;
                // Graupel N
                res[Ng_idx] += conv_n;
#ifdef TRACE_QI
                if(trace)
                    std::cout << "conv ice->graupel dqi " << -conv_q << ", dNi " << -conv_n << "\n";
#endif
#ifdef TRACE_QG
                if(trace)
                    std::cout << "conv ice->graupel dqg " << conv_q << ", dNg " << conv_n << "\n";
#endif
            }
        }

        // ice rain riming
        if(rime_rate_qi > 0.0)
        {
            /////// DT PROBLEM
            float_t rime_qi = min(rime_rate_qi, qi_prime/cc.dt_prime);
            float_t rime_qr = min(rime_rate_qr, qr_prime/cc.dt_prime);
            float_t rime_n = min(min(rime_rate_nr, Nr/cc.dt_prime), Ni/cc.dt_prime);
            /////// DT PROBLEM SOLVED

            // Ice
            res[qi_idx] -= rime_qi;
            // Rain
            res[qr_idx] -= rime_qr;
            // Ice N
            res[Ni_idx] -= rime_n;
            // Rain N
            res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
            if(trace)
                std::cout << "ice rain riming dqr " << -rime_qr << ", dNr " << -rime_n << "\n";
#endif
#ifdef TRACE_QI
            if(trace)
                std::cout << "ice rain riming B dqi " << -rime_qi << ", dNi " << -rime_n << "\n";
#endif
            float_t mult_q = 0.0;
            float_t mult_n = 0.0;
            if(T_prime < T_freeze && ice_multiplication)
            {
                float_t mult_1 = (T_prime - T_mult_min)*const3;
                float_t mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                mult_n = C_mult * mult_1 * mult_2 * rime_qr;
                float_t tmp = mult_n*cc.ice.min_x_riming;
                mult_q = min(rime_qr, tmp);
            }

            float_t delta_e = latent_heat_melt(T_prime) * rime_qi / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_qi > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] -= delta_e;

            if(T_prime >= T_freeze)
            {
                float_t qr_tmp = qr_prime+dt*res[qr_idx];
                float_t Nr_tmp = Nr+res[Nr_idx]*dt;
                float_t x_r = particle_mean_mass(
                    qr_tmp, Nr_tmp,
                    cc.rain.min_x_riming, cc.rain.max_x);
                // Ice
                res[qi_idx] += rime_qi;
                // Rain
                res[qr_idx] += rime_qr;
                // Ice N
                res[Ni_idx] += rime_n;
                // Rain N
                res[Nr_idx] += rime_qr/x_r;
#ifdef TRACE_QR
                if(trace)
                    std::cout << "Melting dqr " << rime_qr << ", dNr " << rime_qr/x_r << "\n";
#endif
#ifdef TRACE_QI
                if(trace)
                    std::cout << "Melting dqi " << rime_qi << ", dNi " << rime_n << "\n";
#endif
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
#ifdef TRACE_QG
                if(trace)
                    std::cout << "Melting with mult dqg " << rime_qi + rime_qr - mult_q << ", dNg " << rime_n << "\n";
#endif
#ifdef TRACE_QI
                if(trace)
                    std::cout << "Melting with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
            }
        }
    }
}


/**
 *
 */
template<class float_t>
void snow_riming(
    float_t &qc_prime,
    float_t &Nc,
    float_t &qr_prime,
    float_t &Nr,
    float_t &qs_prime,
    float_t &Ns,
    float_t &dep_rate_snow,
    float_t &rime_rate_qc,
    float_t &rime_rate_nc,
    float_t &rime_rate_qr,
    float_t &rime_rate_nr,
    float_t &rime_rate_qs,
    float_t &T_prime,
    const double &dt,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(dep_rate_snow > 0.0 && dep_rate_snow >= rime_rate_qc+rime_rate_qr)
    {

        // Depositional growth is stronger than riming growth, therefore ice stays ice
        // ice cloud riming
        if(rime_rate_qc > 0.0)
        {
            /////// DT PROBLEM
            float_t rime_q = min(qc_prime/cc.dt_prime, rime_rate_qc);
            float_t rime_n = min(Nc/cc.dt_prime, rime_rate_nc);
            /////// DT PROBLEM SOLVED
            // Snow
            res[qs_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;

            // Cloud N
            res[Nc_idx] -= rime_n;
#ifdef TRACE_QC
            if(trace)
                if(abs(rime_q) > 0)
                    std::cout << "Snow riming dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
#ifdef TRACE_QS
            if(trace)
                std::cout << "Snow riming dqs " << rime_q << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            if(T_prime < T_freeze && ice_multiplication)
            {
                float_t mult_1 = (T_prime - T_mult_min)*const3;
                float_t mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                float_t mult_n = C_mult * mult_1 * mult_2 * rime_q;
                float_t mult_q = mult_n * cc.ice.min_x_riming;
                mult_q = min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
#ifdef TRACE_QI
                if(trace)
                    std::cout << "Snow riming with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
#ifdef TRACE_QS
                if(trace)
                    std::cout << "Snow riming with mult dqs " << -mult_q << "\n";
#endif
            }
        }
        // snow rain riming
        if(rime_rate_qr > 0.0)
        {
            /////// DT PROBLEM
            float_t rime_q = min(rime_rate_qr, qr_prime/cc.dt_prime);
            float_t rime_n = min(Nr/cc.dt_prime, rime_rate_nr);
            /////// DT PROBLEM SOLVED
            // Snow
            res[qs_idx] += rime_q;
            // Rain
            res[qr_idx] -= rime_q;
            // Rain N
            res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
            if(trace)
                std::cout << "snow rain riming dqr " << -rime_q << ", dNr " << -rime_n << "\n";
#endif
#ifdef TRACE_QS
            if(trace)
                std::cout << "Snow rain riming dqs " << rime_q << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            if(T_prime < T_freeze && ice_multiplication)
            {
                float_t mult_1 = (T_prime - T_mult_min)*const3;
                float_t mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                float_t mult_n = C_mult * mult_1 * mult_2 * rime_q;
                float_t mult_q = mult_n * cc.ice.min_x_riming;
                mult_q = min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
#ifdef TRACE_QI
                if(trace)
                    std::cout << "snow rain riming with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
#ifdef TRACE_QS
                if(trace)
                    std::cout << "snow rain riming with mult dqs " << -mult_q << "\n";
#endif
            }
        }
    } else
    {
        // Depositional growth negative or smaller than riming growth,
        // therefore snow is allowed to convert to graupel and / or hail
        // snow cloud riming
        if(rime_rate_qc > 0.0)
        {
            float_t x_s = particle_mean_mass(qs_prime, Ns,
                cc.snow.min_x_riming, cc.snow.max_x);
            float_t d_s = particle_diameter(x_s,
                cc.snow.a_geo, cc.snow.b_geo);
            /////// DT PROBLEM
            float_t rime_q = min(rime_rate_qc, qc_prime/cc.dt_prime);
            float_t rime_n = min(rime_rate_nc, Nc/cc.dt_prime);
            /////// DT PROBLEM SOLVED
            // Snow
            res[qs_idx] += rime_q;
            // Cloud
            res[qc_idx] -= rime_q;

            // Cloud N
            res[Nc_idx] -= rime_n;
#ifdef TRACE_QC
            if(trace)
                if(abs(rime_q) > 0)
                    std::cout << "snow depos dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
#ifdef TRACE_QS
            if(trace)
                std::cout << "snow depos dqs " << rime_q << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_q < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            float_t mult_q = 0.0;
            if(T_prime < T_freeze && ice_multiplication)
            {
                float_t mult_1 = (T_prime - T_mult_min)*const3;
                float_t mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                float_t mult_n = C_mult * mult_1 * mult_2 * rime_q;
                mult_q = mult_n * cc.ice.min_x_riming;
                mult_q = min(rime_q, mult_q);

                // Ice N
                res[Ni_idx] += mult_n;
                // Ice
                res[qi_idx] += mult_q;
                // Snow
                res[qs_idx] -= mult_q;
#ifdef TRACE_QI
                if(trace)
                    std::cout << "snow depos with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
#ifdef TRACE_QS
                if(trace)
                    std::cout << "snow depos with mult dqs " << -mult_q << "\n";
#endif
            }

            // Conversion snow -> graupel
            if(d_s > D_conv_sg)
            {
                float_t conv_q = (rime_q - mult_q)
                    / (const5*(M_PI/6.0 * rho_ice * d_s*d_s*d_s/x_s -1.0));
                /////// DT PROBLEM
                conv_q = min(qs_prime/cc.dt_prime, conv_q);
                /////// DT PROBLEM SOLVED
                // float_t qs_tmp = qs_prime+dt*res[qs_idx];
                x_s = particle_mean_mass(qs_prime, Ns,
                    cc.snow.min_x_riming, cc.snow.max_x);
                float_t tmp = conv_q / max(x_s, x_conv);
                /////// DT PROBLEM
                float_t conv_n = min(tmp, Ns/cc.dt_prime);
                /////// DT PROBLEM SOLVED

                // Snow
                res[qs_idx] -= conv_q;
                // Graupel
                res[qg_idx] += conv_q;
                // Snow N
                res[Ns_idx] -= conv_n;
                // Graupel N
                res[Ng_idx] += conv_n;
#ifdef TRACE_QS
                if(trace)
                    std::cout << "conversion snow->graupel dqs " << -conv_q << ", dNs " << -conv_n << "\n";
#endif
#ifdef TRACE_QG
                if(trace)
                    std::cout << "conversion snow->graupel dqg " << conv_q << ", dNg " << conv_n << "\n";
#endif
            }
        }

        // Snow rain riming
        if(rime_rate_qs > 0.0)
        {
            /////// DT PROBLEM
            float_t rime_qs = min(rime_rate_qs, qs_prime/cc.dt_prime);
            float_t rime_qr = min(rime_rate_qr, qr_prime/cc.dt_prime);
            float_t rime_n = min(min(rime_rate_nr, Nr/cc.dt_prime), Ns/cc.dt_prime);
            /////// DT PROBLEM SOLVED

            // Snow
            res[qs_idx] -= rime_qs;
            // Rain
            res[qr_idx] -= rime_qr;
            // Snow N
            res[Ns_idx] -= rime_n;
            // Rain N
            res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
            if(trace)
                std::cout << "snow rain riming 2 dqr " << -rime_qr << ", dNr " << -rime_n << "\n";
#endif
#ifdef TRACE_QS
            if(trace)
                std::cout << "snow rain riming 2 dqs " << -rime_qs << ", dNs " << -rime_n << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime) * rime_qr / specific_heat_ice(T_prime);
            // Melting, cooling
            if(rime_qr < 0.0)
                res[lat_cool_idx] += delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] += delta_e;

            float_t mult_q = 0.0;
            float_t mult_n = 0.0;
            if(T_prime < T_freeze && ice_multiplication)
            {
                float_t mult_1 = (T_prime - T_mult_min)*const3;
                float_t mult_2 = (T_prime - T_mult_max)*const4;
                mult_1 = max( 0.0, min(mult_1, 1.0));
                mult_2 = max( 0.0, min(mult_2, 1.0));
                mult_n = C_mult * mult_1 * mult_2 * rime_qr;
                float_t tmp = mult_n*cc.ice.min_x_riming;
                mult_q = min(rime_qr, tmp);
            }
            if(T_prime >= T_freeze)
            {
                float_t qr_tmp = qr_prime+dt*res[qr_idx];
                float_t Nr_tmp = Nr+res[Nr_idx]*dt;
                float_t x_r = particle_mean_mass(
                    qr_tmp, Nr_tmp,
                    cc.rain.min_x_riming, cc.rain.max_x);

                // Snow
                res[qs_idx] += rime_qs;
                // Rain
                res[qr_idx] += rime_qr;
                // Snow N
                res[Ns_idx] += rime_n;
                // Rain N
                res[Nr_idx] += rime_qr/x_r;
#ifdef TRACE_QR
                if(trace)
                    std::cout << "More melting dqr " << rime_qr << ", dNr " << rime_qr/x_r << "\n";
#endif
#ifdef TRACE_QS
                if(trace)
                    std::cout << "More melting dqs " << rime_qs << ", dNs " << rime_n << "\n";
#endif
                float_t delta_e = latent_heat_melt(T_prime) * rime_qr
                                            / specific_heat_ice(T_prime);
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
#ifdef TRACE_QI
                if(trace)
                    std::cout << "More melting with mult dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
#ifdef TRACE_QG
                if(trace)
                    std::cout << "More melting with mult dqg " << rime_qs + rime_qr - mult_q << ", dNg " << rime_n << "\n";
#endif
            }
        }
    }
}


/**
 *
 */
template<class float_t>
void particle_cloud_riming(
    float_t &qc_prime,
    float_t &Nc,
    float_t &T_prime,
    float_t &q1,
    float_t &N1,
    float_t &resq,
    float_t &resn,
    collection_model_constants_t &coeffs,
    particle_model_constants_t &pc1,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    float_t const1 = const0 * pc1.sc_coll_n;
    float_t x_c = particle_mean_mass(qc_prime, Nc, cc.cloud.min_x_riming, cc.cloud.max_x);
    float_t d_c = particle_diameter(x_c, cc.cloud.a_geo, cc.cloud.b_geo);
    float_t x_1 = particle_mean_mass(q1, N1, pc1.min_x_riming, pc1.max_x);
    float_t d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);

    if(qc_prime > q_crit_c && q1 > pc1.q_crit_c
        && d_1 > pc1.d_crit_c && d_c > D_crit_c)
    {
        float_t v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;
        float_t v_c = particle_velocity(x_c, cc.cloud.a_vel, cc.cloud.b_vel) * cc.cloud.rho_v;
        /////// DT PROBLEM
        float_t e_coll_n = min(pc1.ecoll_c/cc.dt_prime, max(const1*(d_c-D_crit_c), ecoll_min/cc.dt_prime));
        /////// DT PROBLEM SOLVED
        float_t e_coll_q = e_coll_n;

        float_t rime_n = M_PI/4.0 * e_coll_n * N1 * Nc
            * (coeffs.delta_n_aa * d_1 * d_1
                + coeffs.delta_n_ab * d_1 * d_c
                + coeffs.delta_n_bb * d_c * d_c)
            * sqrt(coeffs.theta_n_aa * v_1 * v_1
                - coeffs.theta_n_ab * v_1 * v_c
                + coeffs.theta_n_bb * v_c * v_c);
        float_t rime_q = M_PI/4.0 * e_coll_q * N1 * qc_prime
            * (coeffs.delta_q_aa * d_1 * d_1
                + coeffs.delta_q_ab * d_1 * d_c
                + coeffs.delta_q_bb * d_c * d_c)
            * sqrt(coeffs.theta_q_aa * v_1 * v_1
                - coeffs.theta_q_ab * v_1 * v_c
                + coeffs.theta_q_bb * v_c * v_c);
        /////// DT PROBLEM
        rime_q = min(qc_prime/cc.dt_prime, rime_q);
        rime_n = min(Nc/cc.dt_prime, rime_n);
        /////// DT PROBLEM SOLVED
        resq += rime_q;
        // Cloud
        res[qc_idx] -= rime_q;

        // Cloud N
        res[Nc_idx] -= rime_n;
#ifdef TRACE_QC
        if(trace)
            if(abs(rime_q) > 0)
                std::cout << "Particle cloud riming dqc " << -rime_q << ", dNc " << -rime_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
        // Sublimination, cooling
        if(rime_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Deposition, heating
        else
            res[lat_heat_idx] += delta_e;

        // ice multiplication based on Hallet and Mossop
        if(T_prime < T_freeze && ice_multiplication)
        {
            float_t mult_1 = (T_prime - T_mult_min)*const3;
            float_t mult_2 = (T_prime - T_mult_max)*const4;
            mult_1 = max( 0.0, min(mult_1, 1.0));
            mult_2 = max( 0.0, min(mult_2, 1.0));
            float_t mult_n = C_mult * mult_1 * mult_2 * rime_q;
            float_t mult_q = mult_n * cc.ice.min_x_riming;
            mult_q = min(rime_q, mult_q);

            // Ice
            res[qi_idx] += mult_q;
            // Ice N
            res[Ni_idx] += mult_n;
            resq -= mult_q;
#ifdef TRACE_QI
            if(trace)
                std::cout << "Particle cloud riming dqi " << mult_q << ", dNi " << mult_n << "\n";
#endif
        }

        // Enhancement of melting
        if(T_prime > T_freeze && enhanced_melting)
        {
            float_t melt_q = (T_prime-T_freeze)*const5*rime_q;
            float_t melt_n = melt_q/x_1;
            /////// DT PROBLEM
            melt_q = min(q1/cc.dt_prime, melt_q);
            melt_n = min(N1/cc.dt_prime, melt_n);
            /////// DT PROBLEM SOLVED
            resq -= melt_q;
            resn -= melt_n;
            // Rain
            res[qr_idx] += melt_q;
            // Rain N
            res[Nr_idx] += melt_n;
#ifdef TRACE_QR
            if(trace)
                std::cout << "enhancement of melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime) * melt_q
                                        / specific_heat_ice(T_prime);
            // Melting, cooling
            if(melt_q > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] -= delta_e;
        }
    }
}


/**
 *
 */
template<class float_t>
void particle_rain_riming(
    float_t &qr_prime,
    float_t &Nr,
    float_t &T_prime,
    float_t &q1,
    float_t &N1,
    float_t &resq,
    float_t &resn,
    collection_model_constants_t &coeffs,
    particle_model_constants_t &pc1,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(qr_prime > q_crit && q1 > q_crit)
    {

        float_t x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x_riming, cc.rain.max_x);
        float_t d_r = particle_diameter(x_r, cc.rain.a_geo, cc.rain.b_geo);
        float_t x_1 = particle_mean_mass(q1, N1, pc1.min_x_riming, pc1.max_x);
        float_t d_1 = particle_diameter(x_1, pc1.a_geo, pc1.b_geo);

        float_t v_1 = particle_velocity(x_1, pc1.a_vel, pc1.b_vel) * pc1.rho_v;
        float_t v_r = particle_velocity(x_r, cc.rain.a_vel, cc.rain.b_vel) * cc.rain.rho_v;

        float_t rime_n = M_PI/4.0 * N1 * Nr
            * (coeffs.delta_n_aa * d_1 * d_1
                + coeffs.delta_n_ab * d_1 * d_r
                + coeffs.delta_n_bb * d_r * d_r)
            * sqrt(coeffs.theta_n_aa * v_1 * v_1
                - coeffs.theta_n_ab * v_1 * v_r
                + coeffs.theta_n_bb * v_r * v_r);
        float_t rime_q = M_PI/4.0 * N1 * qr_prime
            * (coeffs.delta_n_aa * d_1 * d_1
                + coeffs.delta_q_ab * d_1 * d_r
                + coeffs.delta_q_bb * d_r * d_r)
            * sqrt(coeffs.theta_n_aa * v_1 * v_1
                - coeffs.theta_q_ab * v_1 * v_r
                + coeffs.theta_q_bb * v_r * v_r);

        /////// DT PROBLEM
        rime_q = min(qr_prime/cc.dt_prime, rime_q);
        rime_n = min(Nr/cc.dt_prime, rime_n);
        /////// DT PROBLEM SOLVED
        resq += rime_q;
        // Rain
        res[qr_idx] -= rime_q;
        // Rain N
        res[Nr_idx] -= rime_n;
#ifdef TRACE_QR
        if(trace)
            std::cout << "particle rain riming dqr " << -rime_q << ", dNr " << -rime_n << "\n";
#endif
        float_t delta_e = latent_heat_melt(T_prime) * rime_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(rime_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] += delta_e;

        // ice multiplication based on Hallet and Mossop
        if(T_prime < T_freeze && ice_multiplication)
        {
            float_t mult_1 = (T_prime - T_mult_min)*const3;
            float_t mult_2 = (T_prime - T_mult_max)*const4;
            mult_1 = max( 0.0, min(mult_1, 1.0));
            mult_2 = max( 0.0, min(mult_2, 1.0));
            float_t mult_n = C_mult * mult_1 * mult_2 * rime_q;
            float_t mult_q = mult_n * cc.ice.min_x_riming;
            mult_q = min(rime_q, mult_q);

            // Ice
            res[qi_idx] += mult_q;
            // Ice N
            res[Ni_idx] += mult_n;
#ifdef TRACE_QI
            if(trace)
                std::cout << "particle rain riming dqi " << mult_q << ", dNi " << mult_n << "\n";
        // std::cout << "particle rain riming mult_q: " << mult_n * cc.ice.min_x_riming.getValue()
        //           << ", rime_q: " << rime_q.getValue()
        //           << ", ice_min_x: " << cc.ice.min_x_riming.getValue() << "\n";
#endif
            resq -= mult_q;
        }

        // Enhancement of melting
        if(T_prime > T_freeze && enhanced_melting)
        {
            float_t melt_q = (T_prime-T_freeze)*const5*rime_q;
            float_t melt_n = melt_q/x_1;
            /////// DT PROBLEM
            melt_q = min(q1/cc.dt_prime, melt_q);
            melt_n = min(N1/cc.dt_prime, melt_n);
            /////// DT PROBLEM SOLVED
            resq -= melt_q;
            resn -= melt_n;
            // Rain
            res[qr_idx] += melt_q;
            // Rain N
            res[Nr_idx] += melt_n;
#ifdef TRACE_QR
            if(trace)
                std::cout << "particle rain riming enhancement dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
            float_t delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
            // Melting, cooling
            if(melt_q > 0.0)
                res[lat_cool_idx] -= delta_e;
            // Freezing, heating
            else
                res[lat_heat_idx] -= delta_e;
        }

        // TODO: Add shedding
    }
}


/**
 * Freezing of rain and conversion to ice, graupel, hail
 */
template<class float_t>
void rain_freeze(
    float_t &qr_prime,
    float_t &Nr,
    float_t &T_prime,
    const double &dt,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(T_prime < T_freeze)
    {
        float_t xmax_ice = pow( pow(D_rainfrz_ig/cc.rain.a_geo, 1.0/cc.rain.b_geo), cc.rain.mu );
        float_t xmax_gr  = pow( pow(D_rainfrz_gh/cc.rain.a_geo, 1.0/cc.rain.b_geo), cc.rain.mu );
        float_t fr_q, fr_n, fr_n_i, fr_q_i, fr_n_g, fr_q_g, fr_n_h,
            fr_q_h, fr_n_tmp, fr_q_tmp;
        fr_q = fr_n = fr_n_i = fr_q_i = fr_n_g = fr_q_g = fr_n_h =
            fr_q_h = fr_n_tmp = fr_q_tmp = 0.0;
        float_t Nr_tmp = Nr;
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
            float_t x_r = particle_mean_mass(qr_prime, Nr, cc.rain.min_x_freezing, cc.rain.max_x);
            // Nr_tmp = qr_prime/x_r;
            if(T_prime < T_f)
            {
                // Not sure if this branch makes too much sense
                fr_q = qr_prime;
                fr_n = Nr_tmp;
                float_t lam = pow(cc.rain.g1/cc.rain.g2*x_r, -cc.rain.mu);
                float_t N0 = cc.rain.mu * Nr_tmp * pow(lam, cc.rain.nm1) / cc.rain.g1;
                float_t tmp = lam*xmax_ice;
                fr_n_i = N0/(cc.rain.mu*pow(lam, cc.rain.nm1)) * table_r1.look_lo(tmp);
                fr_q_i = N0/(cc.rain.mu*pow(lam, cc.rain.nm2)) * table_r2.look_lo(tmp);
                tmp = lam*xmax_gr;
                fr_n_g = N0/(cc.rain.mu*pow(lam, cc.rain.nm1)) * table_r1.look_lo(tmp);
                fr_q_g = N0/(cc.rain.mu*pow(lam, cc.rain.nm2)) * table_r2.look_lo(tmp);

                fr_n_h = fr_n - fr_n_g;
                fr_q_h = fr_q - fr_q_g;
                fr_n_g = fr_n_g - fr_n_i;
                fr_q_g = fr_q_g - fr_q_i;
                fr_n_tmp = Nr_tmp/max(fr_n, Nr_tmp);
                fr_q_tmp = qr_prime/max(fr_q, qr_prime);
            } else
            {
                // heterogeneous freezing
                float_t j_het = max(b_HET * ( exp(a_HET * (T_freeze-T_prime)) - 1.0 ), 0.0) / rho_w;

                if(j_het >= 1.0e-20/dt)
                {
                    fr_n = j_het * qr_prime;
                    fr_q = j_het * qr_prime * x_r * cc.rain.c_z; // rain_coeffs
                    float_t lam = pow(cc.rain.g1/cc.rain.g2 * x_r, -cc.rain.mu);
                    float_t N0 = cc.rain.mu * Nr_tmp * pow(lam, cc.rain.nm1) / cc.rain.g1;

                    float_t tmp = lam*xmax_ice;
                    fr_n_i = j_het * N0/(cc.rain.mu*pow(lam, cc.rain.nm2)) * table_r2.look_lo(tmp);
                    fr_q_i = j_het * N0/(cc.rain.mu*pow(lam, cc.rain.nm3)) * table_r3.look_lo(tmp);

                    tmp = lam*xmax_gr;
                    fr_n_g = j_het * N0/(cc.rain.mu*pow(lam, cc.rain.nm2)) * table_r2.look_lo(tmp);
                    fr_q_g = j_het * N0/(cc.rain.mu*pow(lam, cc.rain.nm3)) * table_r3.look_lo(tmp);

                    fr_n_h = fr_n - fr_n_g;
                    fr_q_h = fr_q - fr_q_g;
                    fr_n_g = fr_n_g - fr_n_i;
                    fr_q_g = fr_q_g - fr_q_i;
                    fr_n_tmp = Nr_tmp/max(fr_n, Nr_tmp);
                    fr_q_tmp = qr_prime/max(fr_q ,qr_prime);
                } else
                {
                    fr_n = fr_q = fr_n_i = fr_q_i = fr_n_g = fr_q_g
                        = fr_n_h = fr_q_h = fr_n_tmp = fr_q_tmp = 0.0;
                }
#if defined(TRACE_QR) || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
                if(trace)
                    std::cout << "Hetero j_het: " << j_het
                              << ", val: " << b_HET * ( exp(a_HET * (T_freeze-T_prime)) - 1.0 )
                              << "\n";
#endif
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
        res[Nr_idx] -= fr_n;
        // res[Nr_idx] = Nr_tmp - fr_n;
#if defined(TRACE_QR) || defined(TRACE_QS) || defined(TRACE_QG) || defined(TRACE_QH)
        if(trace)
            std::cout << "qr_prime <= q_crit_fr " << (qr_prime <= q_crit_fr)
                      << "\nT_prime < T_freeze " << (T_prime < T_freeze)
                      << "\nT_prime < T_f " << (T_prime < T_f)
                      << "\nfr_q " << -fr_q
                      << "\nfr_q_i " << fr_q_i
                      << "\nfr_q_tmp " << fr_q_tmp
                      << "\n";
#endif
#ifdef TRACE_QR
        if(trace)
            std::cout << "Freezing dqr " << -fr_q << ", dNr " << -fr_n << "\n";
#endif

        // Snow
        res[qs_idx] += fr_q_i;
        // Snow N
        res[Ns_idx] += fr_n_i;
#ifdef TRACE_QS
        if(trace)
            std::cout << "Freezing dqs " << fr_q_i << ", dNs " << fr_n_i << "\n";
#endif

        // Graupel
        res[qg_idx] += fr_q_g;
        // Graupel N
        res[Ng_idx] += fr_n_g;
#ifdef TRACE_QG
        if(trace)
            std::cout << "Freezing dqg " << fr_q_g << ", dNg " << fr_n_g << "\n";
#endif

        // Hail
        res[qh_idx] += fr_q_h;
        // Hail N
        res[Nh_idx] += fr_n_h;
#ifdef TRACE_QH
        if(trace)
            std::cout << "Freezing dqh " << fr_q_h << ", dNh " << fr_n_h << "\n";
#endif

        float_t delta_e = latent_heat_melt(T_prime) * fr_q
                                    / specific_heat_ice(T_prime);
        // Melting, cooling
        if(fr_q < 0.0)
            res[lat_cool_idx] += delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] += delta_e;
    }
}


/**
 * Ice melts instantanuous to rain or cloud droplets.
 */
template<class float_t>
void ice_melting(
    float_t &qi_prime,
    float_t &qi,
    float_t &Ni,
    float_t &T_prime,
    const double &dt,
    std::vector<float_t> &res,
    model_constants_t &cc)
{
    if(T_prime > T_freeze && qi_prime > 0.0)
    {
        float_t x_i = particle_mean_mass(qi_prime, Ni, cc.ice.min_x_melt, cc.ice.max_x);
        // Complete melting
        /////// DT PROBLEM
        float_t melt_q = qi_prime / cc.dt_prime; // Instantanuous, hence "/dt"
        float_t melt_n = Ni / cc.dt_prime;
        /////// DT PROBLEM SOLVED
        // Remove all ice and all changes that came before due to ice processes.
        // Makes sure, no other ice processes are called subsequently.
        qi = 0.0;
        qi_prime = 0.0;
        Ni = 0;

        // Ice
        res[qi_idx] = -melt_q;
        // Ice N
        res[Ni_idx] = -melt_n;
#ifdef TRACE_QI
        if(trace)
        {
            std::cout << "ice melting: qi and Ni set to 0.0\n";
            std::cout << "ice_melting dqi " << -melt_q << ", dNi " << -melt_n << "\n";
        }
#endif
        // melt into cloud or rain
        if(x_i > cc.cloud.max_x)
        {
            // Rain
            res[qr_idx] += melt_q;
            // Rain N
            res[Nr_idx] += melt_n;
#ifdef TRACE_QR
            if(trace)
                std::cout << "ice melting dqr " << melt_q << ", dNr " << melt_n << "\n";
#endif
        } else
        {
            // Cloud
            res[qc_idx] += melt_q;

            // Cloud N
            res[Nc_idx] += melt_n;
#ifdef TRACE_QC
            if(trace)
                if(abs(melt_q) > 0)
                    std::cout << "ice melting dqc " << melt_q << ", dNc " << melt_n << "\n";
#endif
        }

        float_t delta_e = latent_heat_melt(T_prime) * melt_q / specific_heat_ice(T_prime);
        // Melting, cooling
        if(melt_q > 0.0)
            res[lat_cool_idx] -= delta_e;
        // Freezing, heating
        else
            res[lat_heat_idx] -= delta_e;
    }
}

/**
 * This function evaluates the RHS function of the ODE. It uses the 2 moment
 * cloud scheme after Seifert and Beheng (2006),
 * see https://doi.org/10.1007/s00703-005-0112-4
 * Based on mo_art_2mom_driver.f90 and mo_art_2mom_main.f90 of ICON.
 *
 * @param res On out: system state difference
 * @param y Old system state
 * @param ref Pointer to reference quantities to transform between units
 * @param cc Pointer to constants from the model
 * @param nc Pointer to parameters from the netCDF file
 * @param dt Timestep size
 * @param fixed If True: Reset change of pressure, temperature and ascent (w)
 *              at the end to zero
 */
void RHS_SB(std::vector<codi::RealReverse> &res,
    std::vector<codi::RealReverse> &y,
    const reference_quantities_t &ref,
    model_constants_t &cc,
    nc_parameters_t &nc,
    const double &dt,
    bool fixed=false)
{
    // Limit Nx for every particle
    // y[Nc_idx] = min(max(y[Nc_idx], y[qc_idx]/cc.cloud.max_x), y[qc_idx]/cc.cloud.min_x);
    // y[Nr_idx] = min(max(y[Nr_idx], y[qr_idx]/cc.rain.max_x), y[qr_idx]/cc.rain.min_x);
    // y[Ni_idx] = min(max(y[Ni_idx], y[qi_idx]/cc.ice.max_x), y[qi_idx]/cc.ice.min_x);
    // y[Ns_idx] = min(max(y[Ns_idx], y[qs_idx]/cc.snow.max_x), y[qs_idx]/cc.snow.min_x);
    // y[Ng_idx] = min(max(y[Ng_idx], y[qg_idx]/cc.graupel.max_x), y[qg_idx]/cc.graupel.min_x);
    // y[Nh_idx] = min(max(y[Nh_idx], y[qh_idx]/cc.hail.max_x), y[qh_idx]/cc.hail.min_x);

    // // Decrypt the variables
    codi::RealReverse p = y[p_idx];
    codi::RealReverse T = y[T_idx];
    codi::RealReverse w = y[w_idx];
    codi::RealReverse z = y[z_idx];

    codi::RealReverse S = y[S_idx];
    codi::RealReverse qc = y[qc_idx];
    codi::RealReverse qr = y[qr_idx];
    codi::RealReverse qv = y[qv_idx];
    codi::RealReverse Nc = y[Nc_idx];
    codi::RealReverse Nr = y[Nr_idx];
    codi::RealReverse qi = y[qi_idx];
    codi::RealReverse Ni = y[Ni_idx];
    codi::RealReverse qs = y[qs_idx];
    codi::RealReverse Ns = y[Ns_idx];
    codi::RealReverse qg = y[qg_idx];
    codi::RealReverse Ng = y[Ng_idx];
    codi::RealReverse qh = y[qh_idx];
    codi::RealReverse Nh = y[Nh_idx];
    codi::RealReverse n_inact = y[n_inact_idx];
    codi::RealReverse depo = y[depo_idx];
    codi::RealReverse sub = y[sub_idx];

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

    if(0. > Ni)
        Ni = 0.0;

    if(0. > Ns)
        Ns = 0.0;

    if(0. > Ng)
        Ng = 0.0;

    if(0. > Nh)
        Nh = 0.0;


    // Safety measure: ensure no nonsense

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
    codi::RealReverse z_prime = ref.zref * z;
    // Additional variables such as super saturation
    codi::RealReverse T_c = T_prime - T_freeze;
    codi::RealReverse p_sat = saturation_pressure_water_icon(T_prime);
    codi::RealReverse p_sat_ice = saturation_pressure_ice(T_prime);
    codi::RealReverse S_i = qv_prime * Rv * T_prime / p_sat_ice;
    codi::RealReverse D_vtp = diffusivity(T_prime, p_prime);
    codi::RealReverse e_d = qv_prime * Rv * T_prime; // Could use R_v as well. The difference is minor


    S = convert_qv_to_S(
        p_prime,
        T_prime,
        qv_prime);


    codi::RealReverse s_sw = S - 1.0;   // super saturation over water
    // codi::RealReverse s_sw = e_d / p_sat - 1.0; // super saturation over water
    codi::RealReverse s_si = e_d / p_sat_ice - 1.0; // super saturation over ice
#if defined(TRACE_SAT) ||defined(TRACE_QI) || defined(TRACE_QS) || defined(TRACE_QV) || defined(TRACE_QR) || defined(TRACE_QC) || defined(TRACE_QG) || defined(TRACE_QH)
    if(trace)
    {
        codi::RealReverse x = particle_mean_mass(qr_prime, Nr, cc.rain.min_x_depo, cc.rain.max_x);
        codi::RealReverse q_sat = Epsilon*( p_sat/(p_prime - p_sat) );
        codi::RealReverse q_sat_2 = p_sat/(Rv*T_prime);
        codi::RealReverse p_sat_cosmo = saturation_pressure_water(T_prime);
        codi::RealReverse q_sat_cosmo = Epsilon*( p_sat_cosmo/(p_prime - p_sat_cosmo) );
        codi::RealReverse q_sat_2_cosmo = p_sat_cosmo/(Rv*T_prime);
        std::cout << "\np_sat: " << p_sat.getValue()
                  << "\np_sat_ice: " << p_sat_ice.getValue()
                  << "\nS_i: " << S_i.getValue()
                  << "\ns_sw: " << s_sw.getValue()
                  << "\ns_si: " << s_si.getValue()
                  << "\nS: " << S.getValue()
                  << "\nconvert_qv_to_S: " << convert_qv_to_S(p_prime.getValue(), T_prime.getValue(), qv_prime.getValue())

                //   << "\nS_calc : " << e_d / p_sat
                //   << "\nS_with_q_sat2: " << q_sat_2 * Rv * T_prime / p_sat
                //   << "\nS_calc_cosmo: " << e_d/p_sat_cosmo
                //   << "\nS_with_q_sat_cosmo: " << q_sat_cosmo * Rv * T_prime / p_sat_cosmo
                //   << "\nS_with_q_sat2_cosmo: " << q_sat_2_cosmo * Rv * T_prime / p_sat_cosmo
                //   << "\nDotzek : " << (r_const/r1_const) / (p_prime/p_sat + (r_const/r1_const)-1)
                  << "\nqv: " << qv_prime.getValue()
                //   << "\nqv_sat: " <<  q_sat
                //   << "\nqv_sat1: " << p_sat/(Rv*T_prime)
                //   << "\nq_sat2: " << q_sat_2
                //   << "\nS_qv_sat: " << q_sat * Rv * T_prime / p_sat
                  << "\nqc: " << qc_prime.getValue()
                  << "\nqr: " << qr_prime.getValue()
                  << "\nqi: " << qi_prime.getValue()
                  << "\nqg: " << qg_prime.getValue()
                  << "\nqs: " << qs_prime.getValue()
                  << "\nNc: " << Nc.getValue()
                  << "\nNr: " << Nr.getValue()
                  << "\nmeanNr: " << qr_prime/x
                  << "\nmaxNr: " << qr_prime/cc.rain.min_x_depo
                  << "\nminNr: " << qr_prime/cc.rain.max_x
                  << "\nNi: " << Ni.getValue()
                  << "\nNg: " << Ng.getValue()
                  << "\nNs: " << Ns.getValue()
                  << "\nconvert_S_to_qv: " << convert_S_to_qv(p_prime.getValue(), T_prime.getValue(), S.getValue())
                  << "\nconvert_s_to_qv(S=1): " << convert_S_to_qv(p_prime.getValue(), T_prime.getValue(), codi::RealReverse(1).getValue())
                  << "\n\n";
    }
#endif

    codi::RealReverse rime_rate_qc, rime_rate_qr, rime_rate_qi, rime_rate_qs;
    codi::RealReverse rime_rate_nc, rime_rate_nr;

    // Update vertical velocity parameters
    // rho_v = density correction of terminal fall velocity
    // codi::RealReverse virt = 1 + (Rv/Ra - 1) * qv_prime - qc_prime;
    // codi::RealReverse rho_inter = log(p_prime/(Ra*T_prime*virt)/rho_0);

    codi::RealReverse rho_inter = log(compute_rhoh(p_prime, T_prime, S)/rho_0);
    cc.cloud.rho_v = exp(-rho_vel_c * rho_inter);
    cc.rain.rho_v = exp(-rho_vel * rho_inter);
    cc.graupel.rho_v = exp(-rho_vel * rho_inter);
    cc.hail.rho_v = exp(-rho_vel * rho_inter);
    cc.ice.rho_v = exp(-rho_vel * rho_inter);
    cc.snow.rho_v = exp(-rho_vel * rho_inter);
#ifdef TRACE_QR
    // if(trace)
    // {
    //     codi::RealReverse virt = 1 + (Rv/Ra - 1) * qv_prime - qc_prime;
    //     auto visc = [] (double T)
    //     {
    //         return 2.92332e-6 + 5.20712e-8 * T;
    //     };
    //     auto rho_v_cosmo = visc(293.16)/visc(T.getValue());
    //     std::cout << "rho_0 :" << rho_0
    //               << "\np: " << p_prime
    //               << "\nT: " << T_prime
    //               << "\nS: " << S
    //               << "\nrho: " << compute_rhoh(p_prime, T_prime, S)
    //               << "\nrho_icon: " << p_prime/(Ra*T_prime*virt)
    //               << "\nrho_inter: " << rho_inter
    //               << "\nrho_v_cosmo: " << rho_v_cosmo
    //               << "\nrho_v no pow: " << rho_0/compute_rhoh(p_prime, T_prime, S)
    //               << "\nrho_v: " << cc.cloud.rho_v
    //               << "\nrain rho_v: " << cc.rain.rho_v << "\n";
    // }
#endif

    const double EPSILON = 1.0e-20;

    ////////////// ccn_activation_hdcp2
    if(nuc_type == 0)
    {
        // Not implemented
    }
    else if(nuc_type < 6)
    {
        // Hande et al 2015
        ccn_act_hande(p_prime, w_prime, T_prime, qv_prime, qc_prime, Nc,
            EPSILON, res, cc);

    } else if(nuc_type == 6)
    {
        // Not implemented
    } else
    {
        ccn_act_hande(p_prime, w_prime, T_prime, qv_prime, qc_prime, Nc,
            EPSILON, res, cc);
        // Seifert & Beheng (2006)
        // ccn_act_seifert(p_prime, p, T_prime, T, qv_prime, qv, qc_prime, qc,
        //     Nc, qr, z_prime, dt, w, S, p_sat, ref, res, cc);
    }
    // set cloud droplet number to something inbetween their maximum and
    // minimum sizes
    if(nuc_type > 0)
    {
        double sign = (res[Nc_idx] < 0) ? -1 : 1;
        res[Nc_idx] = sign * min(max(abs(res[Nc_idx]), abs(res[qc_idx])/cc.cloud.max_x),
            abs(res[qc_idx])/cc.cloud.min_x);
    }

    ////////////// ice_nucleation_homhet
//     Homogeneous and heterogeneous ice nucleation based on
//     "A parametrization of cirrus cloud formation: Homogenous
//     freezing of supercooled aerosols" by B. Kaercher and
//     U. Lohmann 2002

//     "Physically based parameterization of cirrus cloud formation
//     for use in global atmospheric models" by B. Kaercher, J. Hendricks
//     and U. Lohmann 2006

    bool ndiag_mask = false;

    // use_prog_in?
    bool use_prog_in = false;
    if(use_hdcp2_het)
    {
        ice_activation_hande(qc_prime, qv_prime, T_prime, S_i,
            n_inact, ndiag_mask, res, cc);
    }  else
    {
        ice_activation_phillips(qc_prime, qv_prime, T_prime,
             S_i, n_inact, use_prog_in, res, cc); // p_sat_ice,
    }

    if(use_prog_in)
    {
        // Stuff that is being done with n_inpot
        // We leave that one out since this is only important for the
        // macrophysical part
    }

    // (optional) homogeneous nucleation using KHL06
    ice_nuc_hom(T_prime, w_prime, p_prime, qv_prime,
        qi_prime, Ni, S_i, p_sat_ice, res, cc);

    cloud_freeze_hom(qc_prime, Nc, T_prime, T_c, res, cc);

    vapor_dep_relaxation(qv_prime, qi_prime, Ni,
        qs_prime, Ns, qg_prime, Ng,
        qh_prime, Nh, s_si, p_sat_ice,
        T_prime, EPSILON, dep_rate_ice, dep_rate_snow, D_vtp, res, cc);

    ////////////// ice-ice collisions
    ice_self_collection(qi_prime, Ni, T_c, res, cc);

    snow_self_collection(qs_prime, Ns, T_prime, res, cc);

    particle_particle_collection(qi_prime, Ni, qs_prime, Ns,
        qg_prime, Ng, T_prime, T_c, res, cc);

    graupel_hail_conv(qc_prime, qr_prime, qi_prime, qg_prime, Ng, qh_prime, Nh,
        p_prime, T_prime, T_c, res, cc);

    hail_collision(qh_prime, Nh, qs_prime, Ns, qi_prime, Ni, T_c, res, cc);

    ////////////// Riming of ice with cloud and rain droplets and conversion to graupel
    riming_cloud_core(qc_prime, Nc, qi_prime, Ni,
        cc.ice, cc.coeffs_icr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(qr_prime, Nr, qi_prime, Ni,
        cc.ice, cc.coeffs_irr, rime_rate_qi, rime_rate_qr, rime_rate_nr, cc);
    ice_riming(qc_prime, Nc, qr_prime, Nr, qi_prime, Ni,
        dep_rate_ice, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
        rime_rate_qi, T_prime, dt, res, cc);

    // snow riming
    riming_cloud_core(qc_prime, Nc, qs_prime, Ns,
        cc.snow, cc.coeffs_scr, rime_rate_qc, rime_rate_nc, cc);
    riming_rain_core(qr_prime, Nr, qs_prime, Ns,
        cc.snow, cc.coeffs_srr, rime_rate_qs, rime_rate_qr, rime_rate_nr, cc);
    snow_riming(qc_prime, Nc, qr_prime, Nr, qs_prime, Ns,
        dep_rate_snow, rime_rate_qc, rime_rate_nc, rime_rate_qr, rime_rate_nr,
        rime_rate_qs, T_prime, dt, res, cc);

#ifdef TRACE_QH
    auto qh_before = res[qh_idx];
    auto Nh_before = res[Nh_idx];
#endif
    //// hail cloud riming
    particle_cloud_riming(qc_prime, Nc, T_prime, qh_prime, Nh,
        res[qh_idx], res[Nh_idx], cc.coeffs_hcr, cc.hail, res, cc);

#ifdef TRACE_QH
    if(trace)
        std::cout << "hail cloud riming dqh " << res[qh_idx]-qh_before << ", dNh " << res[Nh_idx]-Nh_before << "\n";
    qh_before = res[qh_idx];
    Nh_before = res[Nh_idx];
#endif
    //// hail rain riming
    particle_rain_riming(qr_prime, Nr, T_prime, qh_prime, Nh,
        res[qh_idx], res[Nh_idx], cc.coeffs_hrr, cc.hail, res, cc);
#ifdef TRACE_QH
    if(trace)
        std::cout << "hail rain riming dqh " << res[qh_idx]-qh_before << ", dNh " << res[Nh_idx]-Nh_before << "\n";
#endif
#ifdef TRACE_QG
    auto qg_before = res[qg_idx];
    auto Ng_before = res[Ng_idx];
#endif
    //// graupel cloud riming
    particle_cloud_riming(qc_prime, Nc, T_prime, qg_prime, Ng,
        res[qg_idx], res[Ng_idx], cc.coeffs_gcr, cc.graupel, res, cc);

#ifdef TRACE_QG
    if(trace)
        std::cout << "graupel cloud riming dqg " << res[qg_idx]-qg_before << ", dNg " << res[Ng_idx]-Ng_before << "\n";
    qg_before = res[qg_idx];
    Ng_before = res[Ng_idx];
#endif
    //// graupel rain riming
    particle_rain_riming(qr_prime, Nr, T_prime, qg_prime, Ng,
        res[qg_idx], res[Ng_idx], cc.coeffs_grr, cc.graupel, res, cc);
#ifdef TRACE_QG
    if(trace)
        std::cout << "graupel rain riming dqg " << res[qg_idx]-qg_before << ", dNg " << res[Ng_idx]-Ng_before << "\n";
#endif

    rain_freeze(qr_prime, Nr, T_prime, dt, res, cc);

    ice_melting(qi_prime, qi, Ni, T_prime, dt, res, cc);

    snow_melting(qs_prime, Ns, T_prime, res, cc);

    graupel_melting(qg_prime, Ng, T_prime, res, cc);

    hail_melting(qh_prime, Nh, T_prime, res, cc);

    ////////////// Evaporation from melting ice particles
    // TODO: Check if deposition rates are set to zero here

#ifdef TRACE_QS
    auto qs_before = res[qs_idx];
#endif
    evaporation(qv_prime, e_d, p_sat, s_sw, T_prime,
        qs_prime, Ns, res[qs_idx], cc.snow, res, cc.dt_prime);

#ifdef TRACE_QS
    if(trace)
        std::cout << "evaporation dqs " << res[qs_idx]-qs_before << "\n";
#endif

#ifdef TRACE_QG
    qg_before = res[qg_idx];
#endif
    evaporation(qv_prime, e_d, p_sat, s_sw, T_prime,
        qg_prime, Ng, res[qg_idx], cc.graupel, res, cc.dt_prime);
#ifdef TRACE_QG
    if(trace)
        std::cout << "evaporation dqg " << res[qg_idx]-qg_before << "\n";
#endif

#ifdef TRACE_QI
    auto qi_before = res[qi_idx];
#endif
    evaporation(qv_prime, e_d, p_sat, s_sw, T_prime,
        qi_prime, Ni, res[qi_idx], cc.ice, res, cc.dt_prime);
#ifdef TRACE_QI
    if(trace)
        std::cout << "evaporation dqi " << res[qi_idx]-qi_before << "\n";
#endif

    ////////////// Warm rain, ie autoconversion, accretion and rain rain collision
    if(auto_type == 1)
    {
        auto_conversion_kb(qc_prime, Nc, qr_prime, res, cc);

    } else if(auto_type == 2)
    {
        // Not implemented since it appears to be not very interesting
    } else
    {
        auto_conversion_sb(qc_prime, Nc, qr_prime, res, cc);
    }

    rain_self_collection_sb(qr_prime, Nr, res, cc);

    rain_evaporation_sb(qr_prime, Nr, qv_prime, qc_prime,
        T_prime, p_prime, s_sw, p_sat, res, cc);

    sedimentation_explicit(T_prime, S,
        qc_prime, qr_prime, Nr,
        qs_prime, Ns, qi_prime, Ni,
        qh_prime, Nh, qg_prime, Ng,
        p_prime, res, cc);

#ifdef TRACE_QI
    // std::cout << "before limits dqi " << res[qi_idx] << ", dNi " << res[Ni_idx] << "\n";
#endif
#ifdef TRACE_QS
    // std::cout << "before limits dqs " << res[qs_idx] << ", dNs " << res[Ns_idx] << "\n";
#endif
#ifdef TRACE_QR
    // std::cout << "before limits dqr " << res[qr_idx] << ", dNr " << res[Nr_idx] << "\n";
#endif
#if defined TRACE_QC && defined IN_SAT_ADJ
    if(trace)
        std::cout << "before saturation adj dqc " << res[qc_idx] << ", dNc " << res[Nc_idx] << "\n";
#endif
#if defined TRACE_QV && defined IN_SAT_ADJ
    if(trace)
        std::cout << "before saturation adj dqv " << res[qv_idx] << "\n";
#endif
#ifdef TRACE_QG
    // std::cout << "before limits dqg " << res[qg_idx] << ", dNg " << res[Ng_idx] << "\n";
#endif
#ifdef TRACE_QH
    // std::cout << "before limits dqh " << res[qh_idx] << ", dNh " << res[Nh_idx] << "\n";
#endif
    // TODO: Try sat_adj at every 20s for everything and adjust temperature. Then do microphysics.
    // saturation_adjust(T_prime, p_prime, p_sat, qv_prime, qc_prime, res, ref.qref);
    // saturation_adjust_qv(T_prime, p_prime, p_sat, qv_prime, qc_prime, res, ref.qref);
#ifdef IN_SAT_ADJ
    saturation_adjust_meteo(T_prime, p_prime, p_sat, qv_prime, qc_prime, res, ref.qref);
#endif
#ifdef TRACE_QC
    // std::cout << "before limits dqc " << res[qc_idx] << ", dNc " << res[Nc_idx] << "\n";
#endif
#ifdef TRACE_QV
    if(trace)
        std::cout << "end dqv " << res[qv_idx] << "\n";
#endif
    // Set hard limits for particle number
    // double sign = (res[qc_idx] < 0) ? -1 : 1;
    // if(nuc_type > 0)
    //     res[Nc_idx] = sign * min(min(max(abs(res[Nc_idx]), abs(res[qc_idx])/cc.cloud.max_x),
    //         abs(res[qc_idx])/cc.cloud.min_x), 5e9);
    // sign = (res[Nr_idx] < 0) ? -1 : 1;
    // res[Nr_idx] = sign * min(max(abs(res[Nr_idx]), abs(res[qr_idx])/cc.rain.max_x),
    //     abs(res[qr_idx])/cc.rain.min_x);
    // sign = (res[Ni_idx] < 0) ? -1 : 1;
    // res[Ni_idx] = sign * min(max(abs(res[Ni_idx]), abs(res[qi_idx])/cc.ice.max_x),
    //     abs(res[qi_idx])/cc.ice.min_x);
    // sign = (res[Ns_idx] < 0) ? -1 : 1;
    // res[Ns_idx] = sign * min(max(abs(res[Ns_idx]), abs(res[qs_idx])/cc.snow.max_x),
    //     abs(res[qs_idx])/cc.snow.min_x);
    // sign = (res[Ng_idx] < 0) ? -1 : 1;
    // res[Ng_idx] = sign * min(max(abs(res[Ng_idx]), abs(res[qg_idx])/cc.graupel.max_x),
    //     abs(res[qg_idx])/cc.graupel.min_x);
    // sign = (res[Nh_idx] < 0) ? -1 : 1;
    // res[Nh_idx] = sign * min(max(abs(res[Nh_idx]), abs(res[qh_idx])/cc.hail.max_x),
    //     abs(res[qh_idx])/cc.hail.min_x);

#ifdef TRACE_QI
    if(trace)
        std::cout << "end dqi " << res[qi_idx] << ", dNi " << res[Ni_idx] << "\n";
#endif
#ifdef TRACE_QS
    if(trace)
        std::cout << "end dqs " << res[qs_idx] << ", dNs " << res[Ns_idx] << "\n";
#endif
#ifdef TRACE_QR
    if(trace)
        std::cout << "end dqr " << res[qr_idx] << ", dNr " << res[Nr_idx] << "\n";
#endif
#ifdef TRACE_QC
    if(trace)
        std::cout << "end dqc " << res[qc_idx] << ", dNc " << res[Nc_idx] << "\n";
#endif
#ifdef TRACE_QG
    if(trace)
        std::cout << "end dqg " << res[qg_idx] << ", dNg " << res[Ng_idx] << "\n";
#endif
#ifdef TRACE_QH
    if(trace)
        std::cout << "end dqh " << res[qh_idx] << ", dNh " << res[Nh_idx] << "\n";
#endif

    ////// Get back to non-prime
    res[qc_idx] /= ref.qref;
    res[qr_idx] /= ref.qref;
    res[qv_idx] /= ref.qref;
    res[qi_idx] /= ref.qref;
    res[qs_idx] /= ref.qref;
    res[qg_idx] /= ref.qref;
    res[qh_idx] /= ref.qref;
    res[qr_out_idx] /= ref.qref;
    res[qi_out_idx] /= ref.qref;
    res[qs_out_idx] /= ref.qref;
    res[qg_out_idx] /= ref.qref;
    res[qh_out_idx] /= ref.qref;

     // Compute parameters
    codi::RealReverse psat_prime = saturation_pressure_water_icon(T_prime);
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
    // codi::RealReverse c_prime =
    //     4.0*M_PI*H_prime*pow((3.0/(4.0*M_PI*rhow_prime))
    //     *cc.Nc_prime*cc.Nc_prime , 1.0/3.0);
    codi::RealReverse c_prime =
        4.0*M_PI*H_prime*pow((3.0/(4.0*M_PI*rhow_prime))*Nc*Nc , 1.0/3.0);
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
        res[z_idx] = 0;
    } else
    {
        // Calculate pressure and temperature change in non-prime directly
        // Pressure change from ascent (C1), change in partial pressure from water vapor.
        res[p_idx] = -( C1/(1.0 + C2*(qv/(1.0 + qv_prime))) )*( (p*w)/T );
        // This would be a formulation based on warm physics
        // res[T_idx] = ( 1.0/(1.0 + C3*qv + C4*(qc + qr)) )*( -C5*w + C6*qc_third*(S-1.0)
        //     + (C7*qr_delta1 + C8*qr_delta2)*min(S-1.0,0.0) );
        // This uses the latent heat and cooling from each process
        res[T_idx] = (res[lat_heat_idx]/specific_heat_dry_air(T_prime)
            + res[lat_cool_idx]/specific_heat_dry_air(T_prime)
            - gravity_acc * w_prime/specific_heat_dry_air(T_prime)) / ref.Tref;

        res[w_idx] += cc.dw/ref.wref;
        // res[lat_heat_idx] = ( 1.0/(1.0 + C3*qv + C4*(qc + qr)) ) * C6*qc_third*(S-1.0);
        res[z_idx] += w_prime/ref.zref;
    }
    res[S_idx] = (S/p)*res[p_idx] - (S/qv)*( 1.0 - (qv/(C15+qv)) )*( C9*qc_third*(S-1.0)
        + (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0, 0.0) ) - C16*(S/(T*T))*res[T_idx];
    // res[S_idx] = convert_qv_to_S(
    //     p_prime.getValue() + res[p_idx].getValue()*ref.qref,
    //     T_prime.getValue() + res[T_idx].getValue()*ref.Tref,
    //     qv_prime.getValue() + res[qv_idx].getValue()*ref.qref) - S;
#ifdef TRACE_ENV
    if(trace)
    {
        std::cout << "\nPressure " << res[p_idx] << "\n"
                  << "Temperat " << res[T_idx] << "\n"
                  << "ascent   " << res[w_idx] << "\n"
                  << "height   " << res[z_idx] << "\n";
        std::cout << "\nlat_heat " << res[lat_heat_idx]
                  << "\nlat_cool " << res[lat_cool_idx]
                  << "\nlat_heat_air " << res[lat_heat_idx]/specific_heat_dry_air(T_prime)
                  << "\nlat_cool " << res[lat_cool_idx]/specific_heat_dry_air(T_prime)
                  << "\ndT ascent " << - gravity_acc * w_prime / specific_heat_dry_air(T_prime) << "\n";
        std::cout << "w_prime " << w_prime << "\n"
                  << "z_prime " << z_prime << "\n"
                  << "zref " << ref.zref << "\n";
        std::cout << "qsat_prime       " << qsat_prime << "\n"
                  << "qsat_prime COSMO " << (Epsilon - psat_prime)/(p_prime-(1-Epsilon)*psat_prime) << "\n";
        std::cout << "qv_prime " << qv_prime << "\n"
                  << "qc_prime " << qc_prime << "\n";
        std::cout << "T_prime " << T_prime << "\n";
    }
#endif
#ifdef TRACE_SAT
    if(trace)
        std::cout << "End: dS " << res[S_idx]
                  << "\ndS primed " << (S/p_prime)*res[p_idx]*ref.pref
                  - (S/qv_prime)*( 1.0 - (qv_prime/(C15+qv_prime)) )*( C9*qc_third*(S-1.0)
        + (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0, 0.0) ) - C16*(S/(T_prime*T_prime))*res[T_idx]*ref.Tref
                  << "\nS " << S
                  << "\nS_calc " << convert_qv_to_S(
                        p_prime.getValue(),
                        T_prime.getValue(),
                        qv_prime.getValue())
                  << "\nS given new T, p and qv: " << convert_qv_to_S(
                        p_prime.getValue() + res[p_idx].getValue()*ref.qref,
                        T_prime.getValue() + res[T_idx].getValue()*ref.Tref,
                        qv_prime.getValue() + res[qv_idx].getValue()*ref.qref)
                  << ", times dt: " << res[S_idx] * dt
                  << ", Evaporation: " << (C12*qr_delta1 + C13*qr_delta2)*min(S-1.0, 0.0)
                  << ", Condensation: " << C9*qc_third*(S-1.0)
                  << ", heating: " << - C16*(S/(T*T))*res[T_idx]
                  << ", pressure change: " << (S/p)*res[p_idx]
                  << ", qv change part: " << - (S/qv)*( 1.0 - (qv/(C15+qv)) )
                  << "\n\t\t\t\t\tEND\n";
#endif
}


/**
 * This function evaluates the RHS function of the ODE. It uses the 1 moment
 * cloud scheme without ice after Seifert and Beheng (2006),
 * see https://doi.org/10.1007/s00703-005-0112-4
 * and 10.1175/2008JAS2586.1
 *
 * @param res On out: system state difference
 * @param y Old system state
 * @param ref Pointer to reference quantities to transform between units
 * @param cc Pointer to constants from the model
 * @param nc Pointer to parameters from the netCDF file
 * @param fixed If True: Reset change of pressure, temperature and ascent (w)
 *              at the end to zero
 */
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
    codi::RealReverse qi = y[qi_idx];
    codi::RealReverse Ni = y[Ni_idx];
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

    // Change to dimensional variables
    codi::RealReverse p_prime = ref.pref * p;
    codi::RealReverse T_prime = ref.Tref * T;
    codi::RealReverse w_prime = ref.wref * w;
    codi::RealReverse qc_prime = ref.qref * qc;
    codi::RealReverse qr_prime = ref.qref * qr;
    codi::RealReverse qv_prime = ref.qref * qv;

    codi::RealReverse T_c = T_prime - T_freeze;

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
    res[Nc_idx] -= au*x_s_i*2.0;
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
        if(cc.rain_gfak > 0.0)
            gamma_eva = cc.rain_gfak * (1.1e-3/D_r) * exp(-0.2*mue);
        else
            gamma_eva = 1.0;

        // Equation A5 with A9
        codi::RealReverse delta_qv = g_d * Nr * (mue+1.0) / lambda * f_v * s_sw;
        codi::RealReverse delta_nv = gamma_eva * delta_qv/x_r;

        delta_qv = max(-delta_qv, 0.0);
        delta_nv = max(-delta_nv, 0.0);
        delta_qv = min(delta_qv, qv);
        delta_nv = min(delta_nv, Nr);

        res[qv_idx] += delta_qv;
        res[qr_idx] -= delta_qv;
        res[Nr_idx] -= delta_nv;
    }

     // Compute parameters
    codi::RealReverse psat_prime = saturation_pressure_water_icon(T_prime);
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
        *Nc*Nc , 1.0/3.0);
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
/** @} */ // end of group ufunc

#endif
/*
SUBROUTINE cu_cond (                                               &
           pt     , pqv     , p1op   , plflag , pcflag , peflag,   &
           idim   , isc     , iec                                  )

!-----------------------------------------------------------------------------!
! Description:
!
!   The module procedure cu_cond does a saturation adjustment for
!   temperature and specific humidity.
!
!   Method:    Thermodynamic adjustment by instantaneous condensation
!              at constant pressure using a double iteration method.
!              Release of latent heat of condendation and of deposition
!              is considered depending on temperature.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------
! Input data
! ----------
  INTEGER (KIND=iintegers), INTENT (IN) ::  &
     idim ,        & ! array dimension in zonal direction
     isc  ,        & ! start index for first  array computation
     iec             ! end   index for first  array computation

  REAL     (KIND=wp   ),     INTENT (IN) ::  &
     p1op  (idim)    ! reciprocal of pressure, 1.0/p

  LOGICAL                 ,  INTENT (IN) ::  &
     plflag (idim),& ! switch for points where adjustment shall be made
     pcflag,       & ! condensation only (.TRUE)
     peflag          ! evaporation only  (.TRUE)

! Input/Output data
! -----------
  REAL     (KIND=wp   ),     INTENT (INOUT) ::  &
     pt    (idim), & ! temperature on input, adjusted on output
     pqv   (idim)    ! specific humidity on input, adjusted on ouput

! Local scalars and automatic arrays:
! ----------------------------------
  INTEGER (KIND=iintegers) ::  &
    i                ! loop indix

  REAL    (KIND=wp   )     ::  &
    zcond(idim)      ! condensation amount

  REAL    (KIND=wp   )     ::  &
    zhldcp, zcond1, zfacc, zface  ! local storage

!DM+AS>
!_nu
!_nu ! Local Parameters
!_nu REAL (KIND=wp),     PARAMETER :: &
!_nu   Tmpmin = 236.15_wp       , & ! Minium temperature of the mixed-phase temperature range [K]
!_nu   Tmpmax = 267.15_wp       , & ! Maximum temperature of the mixed-phase temperature range [K]
!_nu   exp_mp = 1._wp               ! Exponent in the interpolation formula
!_nu                                    ! for the mixed-phase water fraction [-]
!_nu ! Local Scalars
!_nu REAL (KIND=wp) :: &
!_nu   fr_wat            , & ! Water fraction for the water-ice mixed phase [-]
!_nu   qs_w              , & ! Saturation specific humidity over water [-]
!_nu   qs_i              , & ! Saturation specific humidity over ice [-]
!_nu   qs_m              , & ! Saturation specific humidity for mixed phase
!_nu   qdvdt_w           , & ! First derivative of the saturation specific humidity over water
!_nu                         ! with respect to temperature [K^{-1}]
!_nu   qdvdt_i           , & ! First derivative of the saturation specific humidity over ice
!_nu                         ! with respect to temperature [K^{-1}]
!_nu   qdvdt_m               ! First derivative of the saturation specific humidity
!_nu                         ! with respect to temperature for the mixed phase [K^{-1}]
!_nu
!<DM+AS

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine cu_cond
!------------------------------------------------------------------------------
    zfacc = 0.0_wp
    zface = 0.0_wp
    IF (pcflag) zfacc = 1.0_wp
    IF (peflag) zface = 1.0_wp

    DO i = isc, iec
      zcond(i) = 0.0_wp     ! Initialize condensation variable
      IF(plflag(i)) THEN        ! only, if ascent still continues
!DM+AS>
        ! Water fraction for the mixed water-ice phase as dependent on temperature
        IF (pT(i).LE.Tmpmin) THEN
          fr_wat = 0._wp
        ELSE IF (pT(i).GE.Tmpmax) THEN
          fr_wat = 1._wp
        ELSE
          fr_wat = ((pT(i)-Tmpmin)/(Tmpmax-Tmpmin))**exp_mp
        ENDIF
        ! saturation over water and ice
        qs_w = cc2*EXP( b2w*(pt(i)-b3)/(pt(i)-b4w) )*p1op(i)
        qs_i = cc2*EXP( b2i*(pt(i)-b3)/(pt(i)-b4i) )*p1op(i)
        ! Effective saturation for mixed phase region
        qs_m = fr_wat*qs_w + (1._wp-fr_wat)*qs_i
        qs_w = MIN( 0.5_wp, qs_w )
        qs_i = MIN( 0.5_wp, qs_i )
        qs_m = MIN( 0.5_wp, qs_m )
        qs_w = qs_w / (1.0_wp-rvd_m_o*qs_w)
        qs_i = qs_i / (1.0_wp-rvd_m_o*qs_i)
        qs_m = qs_m / (1.0_wp-rvd_m_o*qs_m)
        ! Effective latent heat of evaporation/sublimation for the mixed phase
        zhldcp = fr_wat*chlcdcp + (1._wp-fr_wat)*chlsdcp
        ! The amount of condensate resulting from the saturation adjustment
        qdvdt_w = c5hlccp * qs_w/(1.0_wp-rvd_m_o*qs_w) / (pt(i)-b4w)**2
        qdvdt_i = c5hlscp * qs_i/(1.0_wp-rvd_m_o*qs_i) / (pt(i)-b4i)**2
        qdvdt_m = fr_wat*qdvdt_w + (1.0_wp-fr_wat)*qdvdt_i
        zcond1  = (pqv(i)-qs_m)/(1.0_wp+qdvdt_m)
        ! switches for evaporation vs condensation
        zcond(i) = zfacc*MAX( zcond1, 0.0_wp )  + &
                   zface*MIN( zcond1, 0.0_wp )
        ! integrate T and qv
        pt(i)    = pt(i) + zhldcp*zcond(i)
        pqv(i)   = pqv(i) - zcond(i)
!<DM+AS
      END IF
    END DO
    !Second iteration
    DO i = isc, iec
      IF( plflag(i) .AND. zcond(i).NE.0.0_wp) THEN  !saturation adjustment
!DM+AS>
        ! Water fraction for the mixed water-ice phase as dependent on temperature
        IF (pT(i).LE.Tmpmin) THEN
          fr_wat = 0._wp
        ELSE IF (pT(i).GE.Tmpmax) THEN
          fr_wat = 1._wp
        ELSE
          fr_wat = ((pT(i)-Tmpmin)/(Tmpmax-Tmpmin))**exp_mp
        ENDIF
        ! saturation over water and ice
        qs_w = cc2*EXP( b2w*(pt(i)-b3)/(pt(i)-b4w) )*p1op(i)
        qs_i = cc2*EXP( b2i*(pt(i)-b3)/(pt(i)-b4i) )*p1op(i)
        ! Effective saturation for mixed phase region
        qs_m = fr_wat*qs_w + (1._wp-fr_wat)*qs_i
        qs_w = MIN( 0.5_wp, qs_w )
        qs_i = MIN( 0.5_wp, qs_i )
        qs_m = MIN( 0.5_wp, qs_m )
        qs_w = qs_w / (1.0_wp-rvd_m_o*qs_w)
        qs_i = qs_i / (1.0_wp-rvd_m_o*qs_i)
        qs_m = qs_m / (1.0_wp-rvd_m_o*qs_m)
        ! Effective latent heat of evaporation/sublimation for the mixed phase
        zhldcp = fr_wat*chlcdcp + (1.0_wp-fr_wat)*chlsdcp
        ! The amount of condensate resulting from the saturation adjustment
        qdvdt_w = c5hlccp * qs_w/(1.0_wp-rvd_m_o*qs_w) / (pt(i)-b4w)**2
        qdvdt_i = c5hlscp * qs_i/(1.0_wp-rvd_m_o*qs_i) / (pt(i)-b4i)**2
        qdvdt_m = fr_wat*qdvdt_w + (1.0_wp-fr_wat)*qdvdt_i
!_cdm "pcflag" and "peflag" are not used during the 2nd iteration. Is that OK?
        zcond1  = (pqv(i)-qs_m)/(1.0_wp+qdvdt_m)
        ! integrate T and qv
        pt(i)    = pt(i) + zhldcp*zcond1
        pqv(i)   = pqv(i) - zcond1
!<DM+AS
      END IF
    END DO

!-----------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------
END SUBROUTINE cu_cond
*/
