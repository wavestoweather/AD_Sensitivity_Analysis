AD with ICON and Processing Results
===================================

Here are various scripts to plot and process data from [netCDF](https://www.unidata.ucar.edu/software/netcdf/) files and my implementation of two moment scheme (similar to [ICON](https://www.dwd.de/EN/research/weatherforecasting/num_modelling/01_num_weather_prediction_modells/icon_description.html)) with [AD using CoDiPack](https://github.com/scicompkl/codipack). For the C++ code we follow the
[doxygen standard](http://www.doxygen.nl/manual/docblocks.html).

Contents
---------

- **scripts:** Python scripts to plot stuff and read files for a quick glance or to create parquet files from simulation output.
- **src:** C++ code. The two moment scheme is under "microphysics".
- **include:** Header files.
- **doc:** Documentation. Go in this folder and type `make html` to create a HTML documentation (recommended) and `make latexpdf` to create a pdf documentation with Latex.
- ***.ipynb:** Jupyter Notebooks to get started. See below for more explanation.


Python and Plotting Prerequisites
---------------------

- [Iris](https://github.com/SciTools/iris) (Recommended to install via conda or from source)
- [pandas](https://pandas.pydata.org/)
- [DASK](https://dask.org/)
- [dask_ml](https://dask-ml.readthedocs.io/en/latest/)
- [progressbar2](https://pypi.org/project/progressbar2/)
- [seaborn](https://seaborn.pydata.org/)
- [matplotlib](https://matplotlib.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib toolkits](https://matplotlib.org/1.4.3/mpl_toolkits/index.html)
- [netCDF4](https://unidata.github.io/netcdf4-python/netCDF4/index.html)
- [pillow, formerly known as PIL](https://pillow.readthedocs.io/en/stable/)
- [xarray, formerly known as xray](http://xarray.pydata.org/en/stable/)
- [GEOS](https://github.com/libgeos/geos)
- [cartopy](https://scitools.org.uk/cartopy/docs/latest/)
- [scikit-learn](https://scikit-learn.org/stable/)
- [datashader](https://datashader.org/index.html)
- [holoviews](http://holoviews.org/)
- [hvplot](https://hvplot.holoviz.org/)
- [fastparquet](https://github.com/dask/fastparquet)

Docs Prerequisites
-------------------------
- [Sphinx](http://www.sphinx-doc.org/en/master/)
- [Read the Docs Sphinx Theme](https://sphinx-rtd-theme.readthedocs.io/en/stable/)
- [recommonmark](https://www.sphinx-doc.org/en/master/usage/markdown.html) to use Markdown
- [Doxygen](http://www.doxygen.nl/index.html)
- [Exhale](https://exhale.readthedocs.io/en/latest/overview.html)
- [Breathe](https://breathe.readthedocs.io/en/latest/)


C++ Prerequisites
-----------------
- [GSL](https://www.gnu.org/software/gsl/)
- [NetCDF C++](https://github.com/Unidata/netcdf-cxx4/releases)
- [Boost](https://www.boost.org/)
- [CoDiPack](https://www.scicomp.uni-kl.de/software/codi/)

Optional Prerequisites
----------------------
- [GNU Parallel](https://www.gnu.org/software/parallel/)


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
and the microphysical process (default: one moment scheme for warm clouds).

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
You can use `./execute.sh` to execute your simulation and do post-processing
to make the data ready for analysis.
All necessary commandline parameters are explained there.

Using Jupyter Notebooks
-----------------------
- **Plot_physics.ipynb:** Plot single microphysical processes with a given range
of input parameters. Helpful if you want to see "what happens" and for debugging.
- **Plot_sensitivities.ipynb:** Use this as a starting point to plot results
from a simulation and to plot the sensitivities.
- **Plot_simulation.ipynb:** Plot the data from netcdf files to get an idea, what's inside.
Or you just use [Panoply](https://www.giss.nasa.gov/tools/panoply/) if you want something more elaborate.

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
