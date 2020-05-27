AD with ICON and Processing Results
===================================

Here are various scripts to plot and process data from [netCDF](https://www.unidata.ucar.edu/software/netcdf/) files and my implementation of two moment scheme (similar to [ICON](https://www.dwd.de/EN/research/weatherforecasting/num_modelling/01_num_weather_prediction_modells/icon_description.html)) with [AD using CoDiPack](https://github.com/scicompkl/codipack). For the C++ code we follow the
[doxygen standard](http://www.doxygen.nl/manual/docblocks.html).

Contents
---------

- **scripts:** Python scripts to plot stuff and read files for a quick glance.
- **src:** C++ code. The two moment scheme is under "microphysics".
- **include:** Header files.
- **doc:** Documentation. Type `make html` to create a HTML documentation (recommended) and `make latexpdf` to create a pdf documentation with Latex.
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

Jupyter Notebooks
-----------------
- **Plot_physics.ipynb:** Plot single microphysical processes with a given range
of input parameters. Helpful if you want to see "what happens" and for debugging.
- **Get_Started.ipynb:** Use this to run a simulation, do some post processing
on the output and plot the results.
- **Plot_netcdf.ipynb:** Plot the data from netcdf files to get an idea, what's inside.
Or you just use Panoply.

I don't want to use Jupyter
---------------------------
Well, that's a weird decision, but here is a way to run everything by yourself:
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
If it is not on, then saturation of 1 is assumed.

To compile the code, simply type
```
make release
```
and your binary is under `build/apps/src/microphysics/trajectories`.

You can use `./execute.sh` to execute your simulation and do post-processing
to make the data ready for analysis.
All necessary commandline parameters are explained there.


