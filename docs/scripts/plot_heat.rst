****************
:mod:`plot_heat`
****************

Load data given a prefix and suffix of the filenames and plot various
figures, mostly heatmaps.
Usage:

- ``-p``, ``--prefix``: Prefix of the files to load such as ``data/sb_ice_traj``.
- ``-s``, ``--suffix``: Suffix of the files to load without datatype such
as ``_start_over_20160922_00`` or ``_start_over`` for derivatives.
- ``-i``, ``--image``: Which plot to generate. Can be either ``results`` or
``derivatives`` to plot lineplots of either the output parameters or the derivatives.
- ``-e``. ``--epsilon``: If ``filter`` is set to true, use this as filter threshold.
- ``-f``, ``--filter``: If you load derivatives, use this to filter derivatives smaller than ``epsilon``.

.. automodule:: scripts.plot_heat
   :members: