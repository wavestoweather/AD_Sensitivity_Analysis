#include "include/microphysics/physical_parameterizations.h"

/**
 * Calculate the separation diameter for wet growth of supercooled
 * water to ice.
 *
 * @param p Pressure in Pascal
 * @param T Temperature in Kelvin
 * @param qw Mixing ratio of rain and cloud water
 * @param qi Mixing ratio of ice
 * @param table Lookup-table (usually ltabdminwgg)
 * @return Separation diameter in m
 */
template <class A>
A wet_growth_diam(
    const A &p,
    const A &T,
    const A &qw,
    const A &qi,
    const table_t &table)
{
    A dmin_loc;
    if(T >= table.x2[table.n2-1])
    {
        dmin_loc = 0.0;
    } else if(T > table.x2[0])
    {
        dmin_loc = 999.99;
    } else
    {
        A p_lok = min( max( p, table.x1[0] ), table.x1[table.n1-1] );
        A tmp = (p_lok - table.x1[0]) * table.odx1;
        uint64_t tmp_d = floor(tmp)+1;
        uint64_t iu =  std::min(tmp_d, table.n1-1 );
        A T_lok = min( max( T, table.x2[0] ), table.x2[table.n2-1] );
        tmp = (T_lok - table.x2[0]) * table.odx2;
        tmp_d = floor(tmp)+1;
        uint64_t ju = std::min(tmp_d, table.n2-1);
        A qw_lok = min( max( qw, table.x3[0] ), table.x3[table.n3-1] );
        tmp = (qw_lok - table.x3[0]) * table.odx3;
        tmp_d = floor(tmp)+1;
        uint64_t ku = std::min(tmp_d, table.n3-1 );
        A qi_lok = min( max( qi, table.x4[0] ), table.x4[table.n4-1] );
        tmp = (qi_lok - table.x4[0]) * table.odx4;
        tmp_d = floor(tmp)+1;
        uint64_t lu = std::min(tmp_d, table.n4-1);

        std::vector<A> h1(16);
        std::vector<A> h2(8);
        std::vector<A> h3(4);
        std::vector<A> h4(2);
        // Tetra linear interpolation by Dmin
        for(uint64_t i=0; i<16; ++i)
            h1[i] = table.get(iu + i/8, ju + (i%8)/4, ku + (i%4)/2, lu+i%2);

        for(uint64_t i=0; i<8; ++i)
        {
            h2[i] = h1[i] + (h1[8+i]-h1[i]) * table.odx1 * (p_lok-table.x1[iu]);
        }
        for(uint64_t i=0; i<4; ++i)
        {
            h3[i] = h2[i] + (h2[4+i]-h2[i]) * table.odx2 * (T_lok-table.x2[iu]);
        }

        h4[0] = h3[0] + (h3[2]-h3[0])   * table.odx3 * (qw_lok-table.x3[ku]);
        h4[1] = h3[1] + (h3[3] - h3[1]) * table.odx3 * (qw_lok-table.x3[ku]);
        dmin_loc = h4[0] + (h4[1]-h4[0]) * table.odx4 * (qi_lok-table.x4[lu]);
#ifdef TRACE_GROWTH
        std::cout << "\nh4_0: " << h4[0]
                  << "\nh4_1: " << h4[1]
                  << "\ntable.odx1: " << table.odx1
                  << "\ntable.odx2: " << table.odx2
                  << "\ntable.odx3: " << table.odx3
                  << "\ntable.odx4: " << table.odx4
                  << "\nku: " << ku
                  << "\nlu: " << lu
                  << "\niu: " << iu
                  << "\nlo: " << lo
                  << "\nju: " << ju
                  << "\njo: " << jo << "\n";
#endif

    }
    return dmin_loc;
}
template codi::RealReverse wet_growth_diam<codi::RealReverse>(
    const codi::RealReverse&, const codi::RealReverse&,
    const codi::RealReverse&, const codi::RealReverse&,
    const table_t&);

/**
 * Setup for bulk sedimentation velocity.
 *
 * @param pc Model constants for a certain particle type.
 */
void setup_bulk_sedi(
    particle_model_constants_t &pc)
{
    get_at(pc.constants, Particle_cons_idx::alfa_n) =
        get_at(pc.constants, Particle_cons_idx::a_vel)
        * tgamma( (get_at(pc.constants, Particle_cons_idx::nu)
        + get_at(pc.constants, Particle_cons_idx::b_vel)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu));
    get_at(pc.constants, Particle_cons_idx::alfa_q) =
        get_at(pc.constants, Particle_cons_idx::a_vel)
        * tgamma( (get_at(pc.constants, Particle_cons_idx::nu)
        + get_at(pc.constants, Particle_cons_idx::b_vel)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu) );
    get_at(pc.constants, Particle_cons_idx::lambda) =
        tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+1.0)
        / get_at(pc.constants, Particle_cons_idx::mu) )
        / tgamma( (get_at(pc.constants, Particle_cons_idx::nu)+2.0)
        / get_at(pc.constants, Particle_cons_idx::mu) );
}

/**
 * Initialize the constants for the particle collection using
 * Seifert & Beheng (2006) Eq. 90-93. This function is for all collections
 * except snow - rain and ice - rain collection.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param c Model constants for particle collection
 */
void init_particle_collection_1(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    collection_model_constants_t &c)
{
    c.delta_n_aa = coll_delta_11(pc1, pc2, 0);
    c.delta_n_ab = coll_delta_12(pc1, pc2, 0);
    c.delta_n_bb = coll_delta_22(pc1, pc2, 0);
    c.delta_q_aa = coll_delta_11(pc1, pc2, 0);
    c.delta_q_ab = coll_delta_12(pc1, pc2, 1);
    c.delta_q_bb = coll_delta_22(pc1, pc2, 1);

    c.theta_n_aa = coll_theta_11(pc1, pc2, 0);
    c.theta_n_ab = coll_theta_12(pc1, pc2, 0);
    c.theta_n_bb = coll_theta_22(pc1, pc2, 0);
    c.theta_q_aa = coll_theta_11(pc1, pc2, 0);
    c.theta_q_ab = coll_theta_12(pc1, pc2, 1);
    c.theta_q_bb = coll_theta_22(pc1, pc2, 1);
}


/**
 * Initialize the constants for the particle collection using
 * Seifert & Beheng (2006) Eq. 90-93. This function is for snow - rain and
 * ice - rain collection.
 *
 * @param pc1 Model constants for a particle type that is being collected
 * @param pc2 Model constants for a particle type that collects
 * @param c Model constants for particle collection
 */
void init_particle_collection_2(
    particle_model_constants_t &pc1,
    particle_model_constants_t &pc2,
    collection_model_constants_t &c)
{
    c.delta_n_aa = coll_delta_11(pc1, pc2, 0);
    c.delta_n_ab = coll_delta_12(pc1, pc2, 0);
    c.delta_n_bb = coll_delta_22(pc1, pc2, 0);
    c.delta_q_aa = coll_delta_11(pc1, pc2, 1);
    c.delta_q_ab = coll_delta_12(pc1, pc2, 1);
    c.delta_q_ba = coll_delta_12(pc2, pc1, 1);
    c.delta_q_bb = coll_delta_22(pc1, pc2, 1);

    c.theta_n_aa = coll_theta_11(pc1, pc2, 0);
    c.theta_n_ab = coll_theta_12(pc1, pc2, 0);
    c.theta_n_bb = coll_theta_22(pc1, pc2, 0);
    c.theta_q_aa = coll_theta_11(pc1, pc2, 1);
    c.theta_q_ab = coll_theta_12(pc1, pc2, 1);
    c.theta_q_ba = coll_theta_12(pc2, pc1, 1);
    c.theta_q_bb = coll_theta_22(pc1, pc2, 1);
}

