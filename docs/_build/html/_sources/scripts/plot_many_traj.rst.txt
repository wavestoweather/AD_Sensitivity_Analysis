*********************
:mod:`plot_many_traj`
*********************

Load derivatives of a simulation or a netCDF file and plot trajectories over satellite images.
Usage:

- ``-c``, ``--csv``: Load an already processed csv file with derivatives.
- ``-n``, ``--netcdf``: Path to a netCDF file. If you want to plot derivatives from a csv file that lacks ``LONGITUDE`` and ``LATITUDE``, you need to provide a netCDF file with that information.
- ``-o``, ``--output``: Directory and name to save the image.
- ``-t``, ``--trajectory``: A list of trajectories to consider here. Using a selected number of trajectories helps reducing the amount of used RAM and makes it easier to see anything in the image.
- ``-r``, ``--rotate``: Rotate longitude and latitude coordinates.
- ``--out_transformed``: In case the csv file had to be transformed because of missing coordinates, store this as a new csv file. Set to "None" if you don't want to store it.
- ``--input_param``: A list of input parameters of the model to plot such as "dzeta" or "dbeta_r".
- ``--output_param``: A list of output parameters of the model to plot such as "T", "qr" or "p".
- ``--max_time``: Maximum time in seconds to plot.

.. automodule:: scripts.plot_many_traj
   :members: