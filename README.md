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

