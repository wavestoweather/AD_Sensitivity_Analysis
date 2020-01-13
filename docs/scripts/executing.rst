*****************
Executing Scripts
*****************

Makefile
=========
The makefile at root supports different variables:

- ``TIMESTEPPER``: Either ``RK4ICE`` for Runge-Kutta 4 method with two-moment scheme and ice after Seifert & Beheng, ``RK4NOICE`` for the same without ice and ``RK4`` for the one-moment scheme without ice
- ``SEASON``: Set to ``SPRING``, ``SUMMER``, ``AUTUMN``, ``WINTER`` or ``SPRING95`` for different HDCP2 simulations (nucleation after Hande et al.)
- ``FLUX``: Wether to use incoming rain, ice, snow and graupel
- ``SOURCE``: If ``WCB``, then data from the netCDF file is read differently

You may execute it using ``make release``.

The makefile under ``docs/`` is used for generating the documentation. Just type ``make html`` in this directory.


Execution
=========
The script ``execute.sh`` can be executed. It sets all necessary variables that
are:

- ``SCALING FACTOR``: The scaling factor. Just leave it at 1.0 if you have no idea
- ``AUTO_TYPE``: Different autoconversion types. Use 1 for KB and Beheng (1994), 3 for Seifert & Beheng. 2 is currently not supported
- ``OUTPUT_FILENAME``: Filename for the results of AD. will be stored at ``data/xx_OUTPUT_FILENAME``, where ``xx`` depends on the used method that is being set via the Makefile
- ``INPUT_FILENAME``: Filepath to the netCDF file from an ICON simulation
- ``TIMESTEP``: Timestep size in seconds for AD
- ``SNAPSHOT_INDEX``: Store a snapshot in ``OUTPUT_FILENAME`` every ``SNAPSHOT_INDEX`` iterations
- ``TRAJ_IDX``: Index of trajectory to load from the netCDF file
- ``TARGET_TIME``: End time in seconds of the simulation
- ``START_OVER``: Set to 1 if the simulation should start over with values from the trajectory everytime the simulation time equals the next datapoint of the netCDF file
- ``FIXED_ITERATION``: Set to 1 if the simulation shall fix pressure, temperature and vertical velocity to values from the netCDF file


Slurm
=====
The file ``make.job`` can be used to execute the script on Mogon I, although Mogon II should work as well but hasn't been tested. The flags set here are explained above under "Execution".