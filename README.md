Algorithmic Differentiation as Sensitivity Analysis in Cloud Microphysics
=========================================================================

Microphysical processes in convective clouds have an impact on radiation in climate models, precipitation and dynamical features. These processes are usually parameterized, i.e. in numerical weather prediction (NWP) models, due to the small scale in which those processes operate and the lack of understanding of some. Understanding the impact of uncertainties in parametrization to uncertainty in cloud responses is crucial for tuning such models. In addition, parameter values can depend on discretization details such as the grid interval, time resolution and other choices made by modeling physical processes. Typical sensitivity analysis investigate few model parameters selected by experts and perturb those based on (joint) proposal distributions that may be refined by experts. In this paper we propose algorithmic differentiation (AD) as a tool to detect the magnitude and time step at which a model state parameter is sensitive to any of hundreds of model parameters. We implement a two-moment scheme and identify the 20 most important parameters to each hydrometeor and validate those parameters by comparing perturbed ensembles and the resulting divergence with predicted deviations. Finally, we use the output of AD for detecting the time step at which perturbing a parameter is necessary with a random forest and a precision and recall of about 75% for each model parameter and hydrometeor.

This repository consists of an implementation of a two-moment scheme (similar to [ICON](https://www.dwd.de/EN/research/weatherforecasting/num_modelling/01_num_weather_prediction_modells/icon_description.html)) with [AD using CoDiPack](https://github.com/scicompkl/codipack) where environment and intial variables are read from [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) files. Python scripts are used for post-processing and data analysis. For the C++ code we follow the
[doxygen standard](http://www.doxygen.nl/manual/docblocks.html). \
We recommend an [anaconda](https://www.anaconda.com/) environment for which you may
use `requirements_conda.txt`. Either way, `requirements_pip.txt` with pip is needed
as well to create a documentation.

We also provide a docker image `mahieronymus/ad_sensitivity:2.3` with all necessary prerequisites installed. 



Contents
---------

- **scripts:** Python scripts for plots and training a random forest.
- **src:** C++ code.
- **include:** Header files for the C++ code.
- **docs:** Documentation. Go in this folder and type `make html` to create a HTML documentation (recommended) or `make latexpdf` to create a pdf documentation with Latex. This documentation is outdated and may only be used to browse the code.
- **exe_scripts:** Bash scripts to compile and run the simulation.
- **configs:** Configuration files for running ensembles during simulation or for specifying parameters for a sensitivity analysis.
- **data:** Example NetCDF files with environment variables and initial datapoints. Output data can be stored here. The results of an ensemble simulations are stored as `data/vladiana_ensembles_postprocess/predictions.nc`.
- **pics:** Default folder to save plots.


Python Prerequisites for Post-Processing and Plotting
------------------------------------------------------
- [Python3](https://www.python.org/) (Tested with v3.7.6)
- [pandas](https://pandas.pydata.org/) (Tested with v1.0.1)
- [progressbar2](https://pypi.org/project/progressbar2/) (Tested with v3.37.1)
- [seaborn](https://seaborn.pydata.org/) (Tested with v0.10.0)
- [scipy](https://www.scipy.org/) (Tested with v1.4.1)
- [scikit-learn](https://scikit-learn.org/stable/index.html) (Tested with v0.24.1)
- [scikit-multilearn](http://scikit.ml/) (Tested with v0.2.0)
- [netCDF4 for Python](https://unidata.github.io/netcdf4-python/netCDF4/index.html) (Tested with v1.5.3)
- [pillow, formerly known as PIL](https://pillow.readthedocs.io/en/stable/) (Tested with v7.0.0)
- [xarray, formerly known as xray](http://xarray.pydata.org/en/stable/) (Tested with v0.15.0)
- [joblib](https://joblib.readthedocs.io/en/latest/) (Tested with v0.16.0)
- [numpy](https://numpy.org/) (Tested with v1.19.1)
- [holoviews](http://holoviews.org/) (Tested with v1.13.2)
- [hvplot](https://hvplot.holoviz.org/) (Tested with v0.5.2)
- [datashader](https://datashader.org/index.html) (Tested with v0.10.0)
- [matplotlib](https://matplotlib.org/) (Tested with v3.1.3)
- [bokeh](https://bokeh.org/) (Tested with v2.0.2)
- [cartopy](https://scitools.org.uk/cartopy/docs/latest/) (Tested with v0.17.0)


Docs Prerequisites
-------------------------
- [Sphinx](http://www.sphinx-doc.org/en/master/) (Tested with v3.0.3)
- [Read the Docs Sphinx Theme](https://sphinx-rtd-theme.readthedocs.io/en/stable/) (Tested with v0.1.8-1)
- [recommonmark](https://www.sphinx-doc.org/en/master/usage/markdown.html) (To use Markdown; Tested with v0.4.0)
- [Doxygen](http://www.doxygen.nl/index.html) (Tested with v1.8.13-4)
- [Exhale](https://exhale.readthedocs.io/en/latest/overview.html) (Tested with v0.2.3)
- [Breathe](https://breathe.readthedocs.io/en/latest/) (Tested with v4.18.1)


C++ Prerequisites
-----------------
- [GCC](https://gcc.gnu.org/) (Tested with v10.2.1-6)
- [OpenMPI](https://www.open-mpi.org/) (Tested with v4.1.2)
- [GSL](https://www.gnu.org/software/gsl/) (Tested with v2.7.1)
- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) (Tested with v4.8.1)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) (Tested with v1.12.1; problems may occur with versions below v1.12)
- [Boost](https://www.boost.org/) (Tested with v1.6.2)
- [CoDiPack](https://www.scicomp.uni-kl.de/software/codi/) (Tested with v2.0; must be at least v2.0)
- [CMake](https://cmake.org/) (v3.7.2 or above; tested with v3.18.4)
- [nlohmann/json](https://github.com/nlohmann/json) (v3.9.1 or above; tested with v3.9.1)
- [(optional) PnetCDF](https://parallel-netcdf.github.io/) (tested with v1.12.2; only if you want to use classic NetCDF-files)
- [(optional) zlib](https://zlib.net/) (tested with v1.2.11; in case you want to compress data with zlib)


Compiling code
---------------
In order to compile the code, create a folder `build` and in this folder type
```
cmake .. -DCMAKE_BUILD_TYPE=release -DTRUSTED_DATA:BOOL=ON -DTARGET=simulation
make -j 6
```
With the option `TRUSTED_DATA` the program just stops simulations for every trajectory with more than 10
consecutive NaNs for every other input column. Without this option, the program
assumes the dataset is broken and terminates directly, giving information which
trajectory and time index starts with consecutive NaNs.

You may change the number for `make` to the number of cores available.
In case the CMake does not find a specific library, you can specify them like this:
```
cmake .. -DNETCDF_INCLUDE_DIR=<path/to/netcdf/include/> -DCODIPACK_INCLUDEDIR=<path/to/codipack/include/> -DHDF5_DIR=<path/to/hdf5/> -DCMAKE_BUILD_TYPE=release
```
You may change the target for timing the microphysics with `-DTARGET=timing` if you wish to check the penalty cost of executing AD at every time step for varying amounts of model state variables and model parameters.

You can add `-DCOMPRESS_OUTPUT:BOOL=ON` to write compressed output with zlib. You need for this HDF5 v1.10.2 or higher. Compression is only supported for sensitivity analysis but not for ensemble simulations at the moment. Default compression level is 6, but you may change it using `-DCOMPRESSION_LEVEL=X` where `X` is the desired compression level. Beware, that level 9 yields only small improvements for trajectories but takes much longer to write.

Running a simulation
---------------------
We provide three bash scripts to execute the code in `exe_scripts`. You may change the value of variable
`NTASKS` to the number of cores available in the script if available.
`perturbed_ensembles.sh` generates the data for estimating the correlation between AD-estimated deviations and
ensemble-estimated deviations. The data is stored under `data/vladiana_ensembles/`.
It runs perturbed ensembles every 30 minutes with a single parameter being perturbed for each member.
For postprocessing, you have to execute the following command in `scripts`:
```
python segment_identifier.py \
    --data_path ../data/vladiana_ensembles/ \
    --verbosity 3 \
    --store_appended_data ../data/vladiana_ensembles_postprocess/predictions.nc \
    --only_append
```
You may now delete the files in `data/vladiana_ensembles/` since we do not need them anymore. You can visualize the results, for example a comparison of the simulated
rain mass density with the gradients:
```
python plot_mse.py \
    --data_path ../data/vladiana_ensembles_postprocess/predictions.nc \
    --backend matplotlib \
    --store_path ../pics/qr_AD_example \
    --xlabel Time\ After\ Ascent\ Begins\ [min] \
    --ylabel "" \
    --title Results\ of\ AD\ Over\ Time \
    --width 1600 \
    --height 900 \
    --plot_variant time_plot \
    --traj 1 \
    --twinlabel Predicted\ Deviation \
    --n_model_params 5 \
    --max_time 7500 \
    --out_parameter QR
```
Or for snow:
```
python plot_mse.py \
    --data_path ../data/vladiana_ensembles_postprocess/predictions.nc \
    --backend matplotlib \
    --store_path ../pics/qs_AD_example \
    --xlabel Time\ After\ Ascent\ Begins\ [min] \
    --ylabel "" \
    --title Results\ of\ AD\ Over\ Time \
    --width 1600 \
    --height 900 \
    --plot_variant time_plot \
    --traj 1 \
    --twinlabel Predicted\ Deviation \
    --n_model_params 1 \
    --max_time 7500 \
    --min_time 2500 \
    --out_parameter QS
```
Plots that show the correlation between AD-estimated and ensemble-estimated
devation for all model state variables at once can be generated with:
```
python plot_mse.py \
    --data_path ../data/vladiana_ensembles_postprocess/predictions.nc \
    --add_zero_sens \
    --plot_types \
    --backend bokeh \
    --store_path ../pics/correlation \
    --xlabel AD-Estimated\ Log\ MSD \
    --ylabel Ensemble-Estimated\ Log\ MSD \
    --title Ensemble\ vs\ AD\ Estimation \
    --width 900 \
    --height 900 \
    --plot_variant correlation \
    --out_parameter all_at_once
```
or only for water vapor mass density:
```
python plot_mse.py \
	--data_path /data/project/wcb/netcdf/vladiana_keyparams_errors/predictions.nc \
	--add_zero_sens \
	--plot_types \
	--backend bokeh \
	--store_path ../pics/correlation \
	--confidence 0.90 \
	--xlabel AD-Estimated\ Log\ MSD \
	--ylabel Ensemble-Estimated\ Log\ MSD \
	--title Ens.\ vs\ AD\ Est. \
	--width 900 \
	--height 900 \
	--plot_variant correlation \
	--out_parameter QV
```
You can generate the plot for other model state variables by exchanging `QV`
with `QS` (snow mass density), `QR` (rain mass density), etc..

`sensitivity_example.sh` is used to run a simulation for comparing the results to the original data.
These can be visualized in folder `scripts`:
```
python plot_simulation_example.py \
    --data_cosmo_path ../data/vladiana_input/ \
    --data_sim_path ../data/comparison/ \
    --width 1200 \
    --height 900 \
    --traj 0 \
    --verbosity 3 \
    --store_path ../pics/comparison \
    --backend bokeh
```

`timing.sh` generates a comma separated output on the terminal with the run time for a
simulation and the number of tracked parameters and model state variables. Be aware
that you have to have compiled the program with target `timing` first!

At last you can plot all input data, which may take a while the first time, with
```
python plot_cosmo.py \
    --pollon 160 \
    --pollat 51 \
    --traj_type conv \
    --store_path ../pics/cosmo_ \
    --verbosity 3 \
    --store_data ../data/all_conv.h5 \
    --data_path ../data/vladiana_complete/
```
This creates a file `data/all_conv.h5` which may be used for additional plots
that can be generated faster.
You can then exchange
```
--store_data ../data/all_conv.h5
```
with
```
--data_path ../data/all_conv.h5
```


Adding New Physics
------------------
Any user-defined functions and processes are under `include/microphysics/user_functions.h`.
Just add your microphysics there to solve the ODE. We use one function to solve the
right-hand side of the ODE, called `RHS_SB`, that calls all necessary functions
to process the microphysics. One may create a similar function that uses different
microphysics.
The datatype `codi::RealReverse` is from CODIPACK and can mostly
be handled like `double`. \
Next you may open `include/microphysics/rk4.h` and either copy `RK4_step(..)` as
`my_RK4_step(..)` and make it use your function,
or simply add pragmas around the function calls like that:
```C++
#ifdef MY_MICROPHYSICS
    my_physics(k, yold, ref, cc, nc);
#else
    RHS(k, yold, ref, cc, nc);
#endif
```
where `MY_MICROPHYSICS` is a name you give to make use of your own function during
compile time (set `TIMESTEPPER=-DMY_MICROPHYSICS` in CMake).
Note, you only need to create your own `RK4_step(..)` if the given arguments are
not sufficient for your own needs! \
At last, go to `src/microphysics/trajectories.cpp` and search for:
```
//////////////// Add any different scheme and model here
```
Add pragmas as before to make use of your timestepper method. If you already
added pragmas in `rk4.h`, it is sufficient to change:
```C++
#elif defined(RK4ICE)
```
to
```C++
#elif defined(RK4) || defined(MY_MICROPHYSICS)
```
Otherwise you need to add additional lines of code, such as:
```C++
#elif defined(MY_MICROPHYSICS)
    my_RK4_step(..);
#elif defined(RK4ICE)
```
If you are using new model parameters, then you have to add them to the enum
`Cons_idx` or if you need model parameters specific to certain
hydrometeors, you have to add them to the enum `Particle_cons_idx` located
in `include/microphysics/constants.h`. The same file has to maps `table_particle_param`
for hydrometeor specific parameters and `table_param` for all other parameters
to map the parameters to a name used for the output.
If your parameters depend on others, you may add such functions to the class
`model_constants_t` which should be called by `setup_model_constants(..)`.
To avoid any confusion, you might want to decorate your version with pragmas.
