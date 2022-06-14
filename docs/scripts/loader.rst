*************
:mod:`loader`
*************

This is a helper file to load dataframes, netCDF and csv files. When executed
as main-file, it needs three arguments:

- prefix to a netCDF file to load such as ``/data/project/wcb/sb_ice_wcb272280_traj``
- suffix of a netCDF file to load such as ``start_over_20160922_00``
- an identifier for storing the csv file such as ``filt_zero_30``

It loads the specified netCDF file with the function ``load_mult_derivates_big(..)``
and stores it as a csv file with the name
``"sb_ice_wcb272280" + "_" + ident + "_" + suffix + ".csv"``.

.. automodule:: scripts.loader
   :members: