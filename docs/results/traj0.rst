Example Results for a Single Trajectory
=======================================

If we use the data from ``2mom_wcb/traj_t000000_p001.nc_wcb`` for the
trajectory with index 0 and :math:`0 \text{s} \leq t \leq 187920 \text{s}`, we get data from within
a region that satisfies the WCB-criterion (see Section 2.3 of Oertel2019_).
One can execute the simulation with the following command, executed at the
root folder of the repository:

.. code-block:: bash

    build/apps/src/microphysics/./trajectories \
    -a 3 \
    -t 0 \
    -s 1 \
    -f 187920 \
    -d 0.01 \
    -i 250 \
    -b 1.0 \
    -o data/wcb_sim \
    -l 2mom_wcb/traj_t000000_p001.nc_wcb \
    -r 0

See the section on Executing_ for a description of these commands. In short,
we are resetting variables such as temperature, saturation, pressure and
ascent every :math:`20 \text{s}`, the timestep size of the netCDF file. Those variables,
except for ascent, are not fixed while calculating the microphysics.
The timestep for the microphysics is :math:`0.01 \text{s}`. We capture the results and
derivatives every 250 timesteps (:math:`=2.5 \text{s}`).
We can plot all the following plots with:

.. code-block: bash

    python plot_many_traj.py -i data/ -t 0


Ice Mixing Ratio
----------------

Let's take a look at the ice mixing ratio :math:`q_{\text{ice}}`:

.. figure:: ../gfx/results/traj_qi.png
    :align: center
    :figclass: align-center

.. figure:: ../gfx/results/sim_qi_tr0.png
    :align: center
    :figclass: align-center

.. figure:: ../gfx/results/sim_qi.png
    :align: center
    :figclass: align-center

    Ice crystal mixing ratio for multiple trajectories.

.. figure:: ../gfx/results/qi_1.png
    :align: center
    :figclass: align-center

    Ice mixing ratio derivatives with the highest impact. The colored background
    stands for timesteps where the WCB-criterion is satisfied.

We can see that the minimum size for ice particles has an impact of :math:`\mathcal{O}(10^{-11})`.


.. figure:: ../gfx/results/qi_2.png
    :align: center
    :figclass: align-center

    Ice mixing ratio derivatives with the second highest impact. The colored background
    stands for timesteps where the WCB-criterion is satisfied.

The parameters with the second highest impact of :math:`\mathcal{O}(10^{-21})`
have largely negative derivatives. :math:`z^{-1}` is the inverse height size
of the air parcel, which is set to :math:`150 \text{m}` for this simulation.
The minimum size of a cloud particle is present as well. :math:`\text{vel}_{b, \text{ice}}`
is a parameter for the vertical velocity of ice particles:

.. math::

    \text{vel}_i = a_{\text{vel}} \cdot x^{b_{\text{vel}}} \cdot \rho_{\text{ice}, v}

where :math:`x` is the particle (mean) mass in kg and :math:`a_{\text{vel}}` is
another particle model constant. Those :math:`a_{\text{vel}}` and :math:`b_{\text{vel}}`
are defined for every particle type separately. This particular :math:`\text{vel}_{b, \text{ice}}`
is used in the ice self collection (see Cotton1986_ or Straka1989_ page 53 for
temperature dependency).
The collection of ice crystals after Seifert2006_ Eq. 62 is implemented as follows:

.. math::

    \frac{\partial q_i}{\partial t} |_{\text{vel}_{b, \text{ice}}} =
    - \frac{\pi}{4} e_{\text{coll}} \delta_{sc, q, \text{ice}}
        N_i q_i D_i^2 \sqrt{\theta_{sc, q, \text{ice}} \text{vel}_{\text{ice}}^2 + 2\cdot s_{\text{vel}, \text{ice}}^2}

with :math:`e_{\text{coll}} = \text{min}(10^{0.035 \cdot T_c - 0.7}, 0.2)`
the mean efficiency, :math:`T_c` is the temperature in degree Celcius, :math:`N_i`
is the amount of ice particles and :math:`q_i` the ice mixing ratio.
:math:`\delta_{sc, q}` and :math:`\theta_{sc, q}` are dimensionless constants
(depending on chosen size distributions and velocity- and diameter-mass relations
of the colliding particles).
:math:`D_i` is the mean particle diameter, which is calculated with the mean
particle mass :math:`x` and some model constants:

.. math::

    D_i = a_{\text{ice}, \text{geo}} \cdot x^{b_{\text{ice}, \text{geo}}}


:math:`s_{\text{vel}, i} = 0.2 \text{m}\text{s}^{-1}` is the constant variance
of the particle velocity.


The parameter :math:`\lambda_{\text{ice}}` and :math:`\text{vel}_{b, \text{ice}}`
are used in the sedimentation process:

.. math::

    \frac{\partial q_i}{\partial t } |_{\lambda_{\text{ice}}, \text{vel}_{b, \text{ice}}} =
    - \rho_{\text{corr}} \cdot \alpha_{\text{ice}, q} \cdot \lambda_{\text{ice}}
    \cdot x^{\text{vel}_{b, \text{ice}}} \cdot q_i \cdot z^{-1}

where :math:`\rho_{\text{corr}} = ( \frac{\rho_{a}'}{\rho_0} )^{-\rho_{\text{vel}}}` is
used for density correction.
The air density :math:`\rho_0` is set to :math:`1.225`, the exponent :math:`\rho_{\text{vel}}`
is set to :math:`0.4`. The density of dry air :math:`\rho_{a}'` can be calculated
using pressure, :math:`p` temperature :math:`T`, saturation :math:`S` and
the gas constant for dry air :math:`R_a` via

.. math::

    \rho_{a}' = \frac{ p - S \cdot p_{S, v} }{R_a \cdot T}

The saturation vapor pressure :math:`p_{S, v}` can be calculated with Murphy2005_.

.. figure:: ../gfx/results/qi_3.png
    :align: center
    :figclass: align-center

    Ice mixing ratio derivatives with the third highest impact. The colored background
    stands for timesteps where the WCB-criterion is satisfied.

There are more derivatives with impact of :math:`\mathcal{O}(10^{-24})`.
The minimum size of graupel and snow particles come to have an effect.
The parameters :math:`\text{geo}_{a, \text{ice}}` and :math:`\text{geo}_{b, \text{ice}}`
are used in calculating the particle diameter and hence for ice and ice
collisions (see above), ice riming from cloud and rain particles, conversion from ice to
graupel, vapor deposition, particle collection with snow and graupel.

Formula for cloud riming rate with a particle x (Seifert2006_):

.. math::
    \begin{split}
    \frac{\partial q_i}{\partial t} |_{\text{cloud riming} \\
    = \frac{\Pi}{4} \cdot e_{\text{coll}} \cdot N_x \cdot q_c \cdot \\
    &\quad ( \delta_{q, aa, icr} \cdot D_x^2 + \delta_{q, ab, icr} \cdot D_x \cdot D_c \\
    &\quad + \delta_{q, bb, icr} \cdot D_c^2 ) \\
    &\quad \cdot \sqrt{ \theta_{q, aa, icr} \cdot v_x^2 \\
    &\quad - \theta_{q, ab, icr} \cdot v_x \cdot v_c \\
    &\quad + \theta_{q, bb, icr} \cdot v_c^2 + s^2_{\text{vel}, \text{ice}} }
    \end{split}

The subscript "icr" stands for ice cloud riming model constants and :math:`v`
is the particle velocity. The rain riming rate can be calculated in a similar
way with model constants used from "irr" (ice rain riming) and properties from
rain instead of cloud droplets.

The parameter :math:`\delta_{sc, \text{ice}}` is part of ice ice collisions as above.

The parameter :math:`\text{ven}_{a, \text{ice}}` is used in vapor deposition
(Seifert2006_, Section 3.3, Equations 37 and 38).

.. math::

    Another formula, maybe I add that later.

The parameters :math:`s_{c, \text{ice}}` and :math:`f_{a, \text{ice}}` are used in depositional
growth of ice particles (as above) and evaporation (Seifert2006_, Equations 76, 77):

.. math::
    Another formula

The parameter :math:`q_{\alpha, \text{ice}}` is part of (TODO):

.. math::
    Another formula



Cloud Droplet Mixing Ratio
--------------------------

.. figure:: ../gfx/results/traj_qc.png
    :align: center
    :figclass: align-center

.. figure:: ../gfx/results/sim_qc_tr0.png
    :align: center
    :figclass: align-center

.. figure:: ../gfx/results/sim_qc.png
    :align: center
    :figclass: align-center

    Cloud droplet mixing ratio for multiple trajectories.

.. figure:: ../gfx/results/qc_1.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qc_2.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qc_3.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qc_4.png
    :align: center
    :figclass: align-center


Rain Droplet Mixing Ratio
--------------------------

.. figure:: ../gfx/results/traj_qr.png
    :align: center
    :figclass: align-center

.. figure:: ../gfx/results/sim_qr_tr0.png
    :align: center
    :figclass: align-center

.. figure:: ../gfx/results/sim_qr.png
    :align: center
    :figclass: align-center

    Rain droplet mixing ratio of multiple trajectories.

.. figure:: ../gfx/results/qr_1.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qr_2.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qr_3.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qr_4.png
    :align: center
    :figclass: align-center


Water Vapor Mixing Ratio
--------------------------

.. figure:: ../gfx/results/traj_qv.png
    :align: center
    :figclass: align-center

.. figure:: ../gfx/results/sim_qv_tr0.png
    :align: center
    :figclass: align-center

.. figure:: ../gfx/results/sim_qv.png
    :align: center
    :figclass: align-center

    Water vapor mixing ratio for multiple trajectories.

.. figure:: ../gfx/results/qv_1.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qv_2.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qv_3.png
    :align: center
    :figclass: align-center


.. figure:: ../gfx/results/qv_4.png
    :align: center
    :figclass: align-center



References
----------

.. [Oertel2019] Oertel, A., Boettcher, M., Joos, H., Sprenger, M., and Wernli, H.,
    "Potential vorticity structure of embedded convection in a warm conveyor
    belt and its relevance for the large-scale dynamics", Weather Clim.
    Dynam. Discuss., https://doi.org/10.5194/wcd-2019-3, in review, 2019.

.. [Cotton1986] F. Mimouni, J. L. Ballard, E. T. Ballard, and R. T. Cotton,
    “Necrotizing Tracheobronchitis: Case Report,” Pediatrics, vol. 77, no. 3, p. 366, Mar. 1986.

.. [Straka1989] J. M. Straka,
    “Hail growth in a highly glaciated central High Plains multi-cellular hailstorm,”
    Ph.D.Diss., University of Wisconsin, Madison, 1989.

.. [Seifert2006] A. Seifert and K. D. Beheng,
    “A two-moment cloud microphysics parameterization for mixed-phase clouds.
    Part 1: Model description,”
    Meteorol. Atmos. Phys., vol. 92, no. 1, pp. 45–66, Feb. 2006, doi: 10.1007/s00703-005-0112-4.

.. [Murphy2005] D. M. Murphy and T. Koop,
    “Review of the vapour pressures of ice and supercooled water for
    atmospheric applications,”
    Quarterly Journal of the Royal Meteorological Society, vol. 131, no. 608,
    pp. 1539–1565, 2005, doi: 10.1256/qj.04.94.
