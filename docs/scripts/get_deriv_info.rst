*********************
:mod:`get_deriv_info`
*********************

Load derivatives of a simulation and store that to a csv
file after being filtered and transformed if needed.
Usage:

- ``-i``, ``--input``: Path to folder with all files to load.
- ``-o``, ``--output``: Directory and name to save the normed csv.
- ``-e``, ``--epsilon``: Value to filter values with. Default is 0.0.
- ``-f``, ``--filter``: Filter values smaller than epsilon if used.
- ``-t``, ``--trajectory``: A list of trajectories to consider here. Using a selected number of trajectories helps reducing the amount of used RAM.
- ``-c``, ``--csv``: Load an already filtered csv file and just norm the derivatives. Deactivates ``-t``, ``-f`` and ``-e``.
- ``-n``, ``--netcdf``: Path to a netCDF file that is used to include longitude and latitude.
- ``-r``, ``--rotate``: Rotate longitude and latitude coordinates.

.. automodule:: scripts.get_deriv_info
   :members:

Usage
=====
Some info might be useful.