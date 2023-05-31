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
cmake .. -DCMAKE_BUILD_TYPE=release -DTRUSTED_DATA:BOOL=ON -DB_EIGHT:BOOL=ON -DTARGET=simulation
make -j 6
```

You may change the number for `make` to the number of cores available.
In case the CMake does not find a specific library, you can specify them like this:
```
cmake .. -DNETCDF_INCLUDE_DIR=<path/to/netcdf/include/> -DCODIPACK_INCLUDEDIR=<path/to/codipack/include/> -DHDF5_DIR=<path/to/hdf5/> -DCMAKE_BUILD_TYPE=release
```
Options are:
- `-DCOMPRESS_OUTPUT:BOOL=ON` Write compressed output with zlib. You need for this HDF5 v1.10.2 or higher. 
Compression is only supported for sensitivity analysis but not for ensemble simulations at the moment. 
Default compression level is 6.
- `-DCOMPRESSION_LEVEL=X` Use compression level `X`, where level 9 yields only small improvements for 
trajectories but takes much longer to write.
- `-DB_EIGHT:BOOL=ON` Output data contains hail. Data may use different time steps.
- `-DTRUSTED_DATA:BOOL=ON` Stop the simulation for every trajectory with more than 10
  consecutive NaNs for every other input column, i.e., trajectories with data from different domains.
  Without this option, the program assumes the dataset is broken and terminates with the first NaN, printing
  information which trajectory and time index starts with consecutive NaNs.
- `-DCCN_AKM:BOOL=ON` Use a a different CCN activation scheme based on Hande et al. (2016) by Annette K. Miltenberger
- `-RELOPTLEVEL=X` Set the optimization level for the compiler. Default optimization is 3.
- `-DSILENT_MODE:BOOL=ON` Suppress prints where possible.

Running a simulation
---------------------
We provide an example bash scripts to execute the code in `exe_scripts`. You may change the value of variable
`NTASKS` to the number of cores available.
The data is stored under `data/simulation/*.nc` and can be analyzed using the various Jupyter notebooks.

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
If you are using new model parameters, then you have to add them to `model_constants_t` and specifically set the 
values in `model_constant_t.setup_model_constants(..)` or in `particle_model_constants_t` if these parameters are
different for different hydrometeors. In addition, please include the default value of these parameters in 
`include/microphysics/constants.h`
If your parameters depend on others, you may add such functions to the class
`model_constants_t` which should be called by `setup_model_constants(..)`.
To avoid any confusion, you might want to decorate your version with pragmas.
