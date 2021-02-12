#include "include/types/nc_parameters_t.h"

void nc_parameters_t::load_vars(
    NcFile &datafile)
{
#ifdef WCB2
    lat_var      = datafile.getVar("latitude");
    lon_var      = datafile.getVar("longitude");
#else
    lat_var      = datafile.getVar("lat");
    lon_var      = datafile.getVar("lon");
#endif
    z_var        = datafile.getVar("z");
#if defined WCB || defined WCB2
    p_var        = datafile.getVar("P");
    t_var        = datafile.getVar("T");
    qc_var       = datafile.getVar("QC");
    qv_var       = datafile.getVar("QV");
    qr_var       = datafile.getVar("QR");
    qi_var       = datafile.getVar("QI");
    qs_var       = datafile.getVar("QS");
    time_rel_var = datafile.getVar("time");
#elif defined MET3D
    p_var        = datafile.getVar("pressure");
    t_var        = datafile.getVar("T");
    qc_var       = datafile.getVar("QC");
    qv_var       = datafile.getVar("QV");
    qr_var       = datafile.getVar("QR");
    qi_var       = datafile.getVar("QI");
    qs_var       = datafile.getVar("QS");
    time_rel_var = datafile.getVar("time_after_ascent");
    time_abs_var = datafile.getVar("time");
    w_var        = datafile.getVar("w");
    S_var        = datafile.getVar("S");
    type_var     = datafile.getVar("type");
#else
    p_var        = datafile.getVar("p");
    t_var        = datafile.getVar("t");
    w_var        = datafile.getVar("w");
    time_rel_var = datafile.getVar("time_rel");
    qc_var       = datafile.getVar("qc");
    qr_var       = datafile.getVar("qr");
    qi_var       = datafile.getVar("qi");
    qs_var       = datafile.getVar("qs");
    qg_var       = datafile.getVar("qg");
    qv_var       = datafile.getVar("qv");
    QIin_var     = datafile.getVar("QIin");
    QSin_var     = datafile.getVar("QSin");
    QRin_var     = datafile.getVar("QRin");
    QGin_var     = datafile.getVar("QGin");
    QIout_var    = datafile.getVar("QIout");
    QSout_var    = datafile.getVar("QSout");
    QRout_var    = datafile.getVar("QRout");
    QGout_var    = datafile.getVar("QGout");
#endif
#ifdef WCB
    // specific humidity
    S_var        = datafile.getVar("RELHUM");
    // Flag wether an effective ascent region is reached
    ascent_flag_var = datafile.getVar("MAP");
    // Potential vorticity (German: Wirbelstaerke)
    // pot_vortic   = datafile.getVar("POT_VORTIC")
#endif

#if defined WCB2 || defined MET3D
    qg_var       = datafile.getVar("QG");
    QIin_var     = datafile.getVar("QI_IN");
    QSin_var     = datafile.getVar("QS_IN");
    QRin_var     = datafile.getVar("QR_IN");
    QGin_var     = datafile.getVar("QG_IN");
    QIout_var    = datafile.getVar("QI_OUT");
    QSout_var    = datafile.getVar("QS_OUT");
    QRout_var    = datafile.getVar("QR_OUT");
    QGout_var    = datafile.getVar("QG_OUT");

    NIin_var     = datafile.getVar("NI_IN");
    NSin_var     = datafile.getVar("NS_IN");
    NRin_var     = datafile.getVar("NR_IN");
    NGin_var     = datafile.getVar("NG_IN");
    NIout_var    = datafile.getVar("NI_OUT");
    NSout_var    = datafile.getVar("NS_OUT");
    NRout_var    = datafile.getVar("NR_OUT");
    NGout_var    = datafile.getVar("NG_OUT");

    Nc_var       = datafile.getVar("NCCLOUD");
    Nr_var       = datafile.getVar("NCRAIN");
    Ni_var       = datafile.getVar("NCICE");
    Ns_var       = datafile.getVar("NCSNOW");
    Ng_var       = datafile.getVar("NCGRAUPEL");

    conv_400_var = datafile.getVar("conv_400");
    conv_600_var = datafile.getVar("conv_600");
    slan_400_var = datafile.getVar("slan_400");
    slan_600_var = datafile.getVar("slan_600");
#endif
#if defined WCB2
    // Flag wether an effective ascent region is reached
    ascent_flag_var = datafile.getVar("WCB_flag");
    // 2h ascent rate after Oertel et al. (2019)
    dp2h_var     = datafile.getVar("dp2h");
#endif
#if defined MET3D && defined TURBULENCE
    qturb_var    = datafile.getVar("Q_TURBULENCE");
#endif
}

void nc_parameters_t::load_params(
    std::vector<size_t> &startp,
    std::vector<size_t> &countp,
    const reference_quantities_t &ref_quant,
    uint64_t num_sub_steps)
{
    // startp[0] <- ensemble (if available)
    // startp[1] <- trajectory id
    // startp[2] <- timestep


#if defined WCB || defined WCB2
    countp[0]++;
    countp[0]++;
    this->z_var.getVar(startp, countp, this->z.data());
    this->lat_var.getVar(startp, countp, this->lat.data());
    this->lon_var.getVar(startp, countp, this->lon.data());
    countp[0]--;
    countp[0]--;
#elif defined MET3D
    this->z_var.getVar(startp, countp, &(this->z));
    // std::cout << "after z\n";
    countp[2]++;
    this->lat_var.getVar(startp, countp, this->lat.data());
    // std::cout << "after lat\n";
    this->lon_var.getVar(startp, countp, this->lon.data());
    // std::cout << "after lon\n";
    this->w_var.getVar(startp, countp, this->w.data());
    // std::cout << "after w\n";
    countp[2]--;
#else
    countp[1]++;
    countp[1]++;
    this->z_var.getVar(startp, countp, this->z.data());
    this->lat_var.getVar(startp, countp, this->lat.data());
    this->lon_var.getVar(startp, countp, this->lon.data());
    countp[1]--;
    countp[1]--;
#endif

#if defined WCB || defined WCB2
    int map = 0;
    this->ascent_flag_var.getVar(startp, countp, &map);
    this->ascent_flag = (map > 0) ? true : false;
#endif
#if defined WCB2
    this->dp2h_var.getVar(startp, countp, &map);
    this->dp2h = (map > 0) ? true : false;
#endif
#if defined MET3D
    int map = 0;
#endif
#if defined MET3D || defined WCB2
    this->conv_400_var.getVar(startp, countp, &map);
    this->conv_400 = (map > 0) ? true : false;
    this->conv_600_var.getVar(startp, countp, &map);
    this->conv_600 = (map > 0) ? true : false;
    this->slan_400_var.getVar(startp, countp, &map);
    this->slan_400 = (map > 0) ? true : false;
    this->slan_600_var.getVar(startp, countp, &map);
    this->slan_600 = (map > 0) ? true : false;
#endif
    // std::cout << "after flags\n";
    this->t_var.getVar(startp, countp, &(this->t));
    this->p_var.getVar(startp, countp, &(this->p));
    this->time_rel_var.getVar(startp, countp, &(this->time_rel));
    this->qc_var.getVar(startp, countp, &(this->qc));
    this->qr_var.getVar(startp, countp, &(this->qr));
    this->qi_var.getVar(startp, countp, &(this->qi));
    this->qs_var.getVar(startp, countp, &(this->qs));
    this->qv_var.getVar(startp, countp, &(this->qv));
    // std::cout << "after env\n";
#if !defined WCB2 && !defined MET3D
    // We are reading in hPa. Convert to Pa
    this->p        *= 100;
#endif
    this->p        /= ref_quant.pref;
    this->t        /= ref_quant.Tref;
    this->qc       /= ref_quant.qref;
    this->qr       /= ref_quant.qref;
    this->qv       /= ref_quant.qref;
    this->qi       /= ref_quant.qref;
    this->qs       /= ref_quant.qref;

#if defined WCB2 || defined MET3D
    this->NRin_var.getVar(startp, countp, &(this->NRin));
    this->NIin_var.getVar(startp, countp, &(this->NIin));
    this->NSin_var.getVar(startp, countp, &(this->NSin));
    this->NGin_var.getVar(startp, countp, &(this->NGin));

    this->NRout_var.getVar(startp, countp, &(this->NRout));
    this->NIout_var.getVar(startp, countp, &(this->NIout));
    this->NSout_var.getVar(startp, countp, &(this->NSout));
    this->NGout_var.getVar(startp, countp, &(this->NGout));

    this->Nc_var.getVar(startp, countp, &(this->Nc));
    this->Nr_var.getVar(startp, countp, &(this->Nr));
    this->Ni_var.getVar(startp, countp, &(this->Ni));
    this->Ns_var.getVar(startp, countp, &(this->Ns));
    this->Ng_var.getVar(startp, countp, &(this->Ng));

    this->NRin     /= ref_quant.Nref;
    this->NIin     /= ref_quant.Nref;
    this->NSin     /= ref_quant.Nref;
    this->NGin     /= ref_quant.Nref;

    this->NRout    /= ref_quant.Nref;
    this->NIout    /= ref_quant.Nref;
    this->NSout    /= ref_quant.Nref;
    this->NGout    /= ref_quant.Nref;

    this->Nc       /= ref_quant.Nref;
    this->Nr       /= ref_quant.Nref;
    this->Ni       /= ref_quant.Nref;
    this->Ns       /= ref_quant.Nref;
    this->Ng       /= ref_quant.Nref;

    this->NRin     = abs(this->NRin);
    this->NIin     = abs(this->NIin);
    this->NSin     = abs(this->NSin);
    this->NGin     = abs(this->NGin);
#endif

#if defined MET3D && defined TURBULENCE
    this->qturb_var.getVar(startp, countp, &(this->qturb));
    this->qturb /= ref_quant.qref;
#endif
#if !defined WCB
    this->S = 1.0;
    this->qg_var.getVar(startp, countp, &(this->qg));
    this->QIin_var.getVar(startp, countp, &(this->QIin));
    this->QSin_var.getVar(startp, countp, &(this->QSin));
    this->QRin_var.getVar(startp, countp, &(this->QRin));
    this->QGin_var.getVar(startp, countp, &(this->QGin));
    this->QIout_var.getVar(startp, countp, &(this->QIout));
    this->QSout_var.getVar(startp, countp, &(this->QSout));
    this->QRout_var.getVar(startp, countp, &(this->QRout));
    this->QGout_var.getVar(startp, countp, &(this->QGout));

    this->qg       /= ref_quant.qref;
    this->QRin     /= ref_quant.qref;
    this->QRout    /= ref_quant.qref;
    this->QIin     /= ref_quant.qref;
    this->QIout    /= ref_quant.qref;
    this->QSin     /= ref_quant.qref;
    this->QSout    /= ref_quant.qref;
    this->QGin     /= ref_quant.qref;
    this->QGout    /= ref_quant.qref;

    this->QRin     = abs(this->QRin);
    this->QRout    = abs(this->QRout);
    this->QIin     = abs(this->QIin);
    this->QIout    = abs(this->QIout);
    this->QSin     = abs(this->QSin);
    this->QSout    = abs(this->QSout);
    this->QGin     = abs(this->QGin);
    this->QGout    = abs(this->QGout);
#endif
#ifdef MET3D
    this->S_var.getVar(startp, countp, &(this->S));
    this->S /= 100; // from percentage
    this->time_rel_var.getVar(startp, countp, &(this->time_rel));
    // std::cout << "after time_rel\n";
    // there is only a single value for that since each type
    // is divided into different files in the input
    if(std::strcmp(this->type[0], "") == 0)
    {
        this->type_var.getVar(this->type);
    }
    // std::cout << "after type\n";
    this->w[0]     /= ref_quant.wref;
    this->w[1]     /= ref_quant.wref;
#elif !defined WCB && !defined WCB2
    this->w_var.getVar(startp, countp, this->w.data());
    countp[1]++;
    this->w_var.getVar(startp, countp, this->w.data());
    countp[1]--;

    this->w[0]     /= ref_quant.wref;
    this->w[1]     /= ref_quant.wref;
#elif defined WCB
    this->qc /= 1.0e6;
    this->qr /= 1.0e6;
    this->qv /= 1.0e6;
    this->qi /= 1.0e6;
    this->qs /= 1.0e6;
    // Calculate w by getting the z-coordinates
    // and divide it by the amount of substeps
    this->w[0] = (this->z[1] - this->z[0]) / 20.0;
    this->w[1] = (this->z[2] - this->z[1]) / 20.0;
    this->S_var.getVar(startp, countp, &(this->S));
    this->S /= 100; // from percentage
#else
    // WCB 2
    // Calculate w by getting the z-coordinates
    // and divide it by the amount of substeps
    this->w[0] = (this->z[1] - this->z[0]) / 20.0;
    this->w[1] = (this->z[2] - this->z[1]) / 20.0;
#endif
    this->dlat = (this->lat[1] - this->lat[0]) / num_sub_steps;
    this->dlon = (this->lon[1] - this->lon[0]) / num_sub_steps;
}

void nc_parameters_t::init_params(
    uint32_t n,
    uint32_t n_time)
{
    this->n_trajectories = n;
    this->n_timesteps = n_time;
#if defined MET3D
    this->w.resize(2);
    this->time_abs.resize(n_time);
    this->type[0] = (char*) "";
#else
    this->w.resize(2);
    this->z.resize(4);
#endif
    this->lat.resize(2);
    this->lon.resize(2);
}