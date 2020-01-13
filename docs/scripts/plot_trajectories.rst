************************
:mod:`plot_trajectories`
************************
Plot data from a filtered csv dataset and from a netcdf file.
Usage:

- ``-c``, ``--csv``: File with derivatives.
- ``-n``, ``--netcdf``: File with netCDF info from a Cosmo or Icon simulation.
- ``-o``, ``--outdir``: Path to store plots and transformed data if available.
- ``-r``, ``--rotate``: Use this if longitudes and latitudes have to be rotated on the csv file.
- ``-i``, ``--input_config``: Config file that tells which plots to generate. Default plots are S and T with dzeta each.
- ``-t``, ``--trajectory``: A list of trajectories to consider here. Using a selected number of trajectories makes plots easier to read. If -1 is used, all trajectories will be plotted.
- ``--resolution``: Resolution for the background maps. Default is ``f``. You may need to install basemap-data-hires first for this or choose a lower resolution.
- ``-m``, ``--max_time``: Maximum time to plot in seconds.

.. automodule:: scripts.plot_trajectories
   :members: