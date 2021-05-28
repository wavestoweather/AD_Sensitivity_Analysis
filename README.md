Algorithmic Differentiation as Sensitivity Analysis in Microphysics
====================================================================

Microphysical processes in convective clouds have an impact on radiation in climate models, precipitation and dynamical features. These processes are usually parametrized, i.e. in numerical weather prediction (NWP) models, due to the small scale in which those processes operate and the lack of understanding of some. Understanding the impact of uncertainties in parametrization to uncertainty in cloud responses is crucial for tuning such models. In addition, parameter values can depend on discretization details such as the grid interval, time resolution and other choices made by modeling physical processes. Typical sensitivity analysis investigate few model parameters selected by experts and perturb those based on (joint) proposal distributions that may be refined by experts. In this paper we propose algorithmic differentiation (AD) as a tool to detect the magnitude and time step at which a model state parameter is sensitive to any of hundreds of model parameters. We implement a two-moment scheme and identify the 20 most important parameters to each hydrometeor and validate those parameters by comparing perturbed ensembles and the resulting divergence with predicted deviations. Finally, we use the output of AD for detecting the time step at which perturbing a parameter is necessary with a random forest and a precision and recall of about 75% for each model parameter and hydrometeor.

This repository consists of an implementation of a two-moment scheme (similar to [ICON](https://www.dwd.de/EN/research/weatherforecasting/num_modelling/01_num_weather_prediction_modells/icon_description.html)) with [AD using CoDiPack](https://github.com/scicompkl/codipack) where environment and intial variables are read from [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) files. Python scripts and Jupyter notebooks are used for post-processing, data analysis and for training a random forest for detecting starts of segments where perturbing a parameter has the highest impact. For the C++ code we follow the
[doxygen standard](http://www.doxygen.nl/manual/docblocks.html). \
We recommend an [anaconda](https://www.anaconda.com/) environment for which you may
use `requirements_conda.txt`. Either way, `requirements_pip.txt` with pip is needed
as well for the documentation.


Contents
---------

- **scripts:** Python scripts for plots and training a random forest.
- **src:** C++ code. The two moment scheme is under "microphysics".
- **include:** Header files.
- **docs:** Documentation. Go in this folder and type `make html` to create a HTML documentation (recommended) and `make latexpdf` to create a pdf documentation with Latex.
- ***.ipynb:** Jupyter Notebooks to get started. See below for more explanation.
- **exe_scripts:** Bash scripts to make and run the simulation. Examples for SLURM scripts are there as well.
- **configs:** Configuration files for running ensembles during simulation.
- **data:** Example NetCDF files with environment variables and initial datapoints. Output data can be stored here.
- **pics:** Default folder to save plots.


Python Prerequisites for Post-Processing and Plotting
------------------------------------------------------
- [Python3](https://www.python.org/) (Tested with v3.7.6)
- [pandas](https://pandas.pydata.org/) (Tested with v1.0.1)
- [DASK](https://dask.org/) (Tested with v2.16.0)
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
- [Read the Docs Sphinx Theme](https://sphinx-rtd-theme.readthedocs.io/en/stable/) (Tested with v0.4.3)
- [recommonmark](https://www.sphinx-doc.org/en/master/usage/markdown.html) (To use Markdown; Tested with v0.6.0)
- [Doxygen](http://www.doxygen.nl/index.html) (Tested with v1.8.13-4)
- [Exhale](https://exhale.readthedocs.io/en/latest/overview.html) (Tested with v0.2.3)
- [Breathe](https://breathe.readthedocs.io/en/latest/) (Tested with v4.18.1)


C++ Prerequisites
------------------
- [GCC](https://gcc.gnu.org/) (Tested with v6.3.0)
- [GSL](https://www.gnu.org/software/gsl/) (Tested with v2.3)
- [NetCDF C++](https://github.com/Unidata/netcdf-cxx4/releases) (Tested with v4.7.3)
- [Boost](https://www.boost.org/) (Tested with v1.3.4)
- [CoDiPack](https://www.scicomp.uni-kl.de/software/codi/) (Tested with v1.8.0)


Optional Prerequisites
-----------------------
- [GNU Parallel](https://www.gnu.org/software/parallel/) (Executing multiple processes with ease; Tested with v20161222)
- [JupyterLab](https://jupyter.org/) (Working with Python and Plotting but in an easier way; Tested with v1.2.6)
- JupyterLab extension for geoviews: `jupyter labextension install @pyviz/jupyterlab_pyviz`
- [Panel](https://panel.holoviz.org/) (Tested with v0.9.5)


Compiling code
---------------
You can alter the Makefile to your choice. \
The variable `SEASON` can be set to
`SPRING`, `SUMMER`, `AUTUMN`, `WINTER` and `SPRING95` and sets variables used
in Hande et al. nucleation (not used by default, so you can ignore it). \
With `FLUX` the microphysics takes incoming particles from above into account (default: on).
Leave it empty, if you do not want to use it. \
`SOURCE` is used to toggle different input files (default is `WCB2`; just leave it with that).
`SAT_CALC` is set in `SOURCE` as well to calculate the saturation at every step using `qv*Rv*T/p_sat`
with `Rv` the gas constant for water vapor and `p_sat` the
saturation vapor pressure over a flat surface of liquid water (see `physical_parametrizations.h`).
If it is not on, then saturation of 1 is assumed. \
Variable `TIMESTEPPER` can be used to set the method to use (default: Runge Kutta)
and the microphysical process (default: two moment scheme with ice phase including hail and graupel).

To compile the code, simply type
```
make release
```
and your binary is under `build/apps/src/microphysics/trajectories` or you
type
```
make scan
```
to create a binary for test cases that scan different microphysical processes
under `build/apps/src/scratch/scan`.


Running a simulation
---------------------
You can use `./execute.sh` in the folder `exe_scripts/` to compile the code, execute a single simulation and do post-processing
to make the data ready for analysis and to plot the simulation results.
All necessary commandline parameters are explained there.
In order to run multiple ensembles, use `./execute_ensembles.sh`.


Using Jupyter Notebooks
-----------------------
- **Plot_MSE.ipynb:** Get the most important parameters for each output parameter, calculate and plot correlation and
predicted errors over time.
- **Plot_Segment_identifier.ipynb:** An interactive notebook that shows different ways to identify segments, comparisons to sensitivities and plots confusion matrix for different input parameters and detection methods.
- **Plot_physics.ipynb:** Plot single microphysical processes with a given range
of input parameters. Helpful if you want to see "what happens" in each processes and for debugging.
- **Plot_sensitivities.ipynb:** Use this as a starting point to plot results
from a simulation and to plot the sensitivities.


Adding New Physics
------------------
Any user-defined functions and processes are under `include/microphysics/user_functions.h`.
Just add your microphysics there and solve the ODE. A simple example is `RHS(..)` that
you may copy and use for your own `my_physics(..)`.
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
compile time (set `TIMESTEPPER=-DMY_MICROPHYSICS` in the Makefile).
Note, you only need to create your own `RK4_step(..)` if the given arguments are
not sufficient for your own needs! \
At last, go to `src/microphysics/trajectories.cpp` and search for:
```
//////////////// Add any different scheme and model here
```
Add pragmas as before to make use of your timestepper method. If you already
added pragmas in `rk4.h`, it is sufficient to change:
```C++
#if defined(RK4)
```
to
```C++
#if defined RK4 || defined MY_MICROPHYSICS
```
Otherwise you need to add additional lines of code, such as:
```C++
#if defined(MY_MICROPHYSICS)
    my_RK4_step(..);
#elif defined(RK4)
```
You can store your model constants and add new ones in the struct
`model_constants_t`, which is located in `include/microphysics/types.h`.
Make sure to use the type `codi::RealReverse` if you want to get
gradients/sensitivites for this parameter. You can set some standard values there
or, in case you need more sophisticated calculations or if you want to change
them during simulation, use `setup_model_constants(..)` in
`include/microphysics/physical_parametrizations.h` and change it to your needs.
Decorate your version with pragmas as explained before.


Getting Gradients for New Physics
---------------------------------
You successfully added parameters to `model_constants_t` and set them up.
Go to `include/microphysics/gradient_handle.h`. The `register_everything(..)`
registers all input parameters on a tape. Add
```C++
tape.registerInput(parameter);
```
for every parameter in your model that is stored in a `model_constants_t`. \
In `get_gradients(..)`, add
```C++
y_diff[ii][idx] == cc.parameter.getGradient();
```
for every parameter in your model with `idx` a running variable for the parameters. \
In `include/microphysics/constants.h` add with your own pragmas the number
of sensitivities your looking for with `num_par` and add the number of output
parameters with `num_comp`. You also need to add indizes of all the output
parameters. \
Change the contents of `output_par_idx` and `output_grad_idx` to your output
parameters and input parameters. These vectors are used to write the headers
of your output files.