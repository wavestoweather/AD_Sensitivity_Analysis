AD with ICON and Processing Results
===================================

Here are various scripts to plot and process data from [netCDF](https://www.unidata.ucar.edu/software/netcdf/) files and my implementation of two moment schemes (similar to [ICON](https://www.dwd.de/EN/research/weatherforecasting/num_modelling/01_num_weather_prediction_modells/icon_description.html)) with [AD using CoDiPack](https://github.com/scicompkl/codipack).
The two moment schemes is also available. Documentation is done using [Sphinx](http://www.sphinx-doc.org/en/master/).

Contents
---------

- **scripts:** Python scripts to plot stuff and read files for a quick glance.
- **src:** C++ code. The two moment scheme is under "microphysics".
- **include:** Header files.
- **doc:** Documentation. Type `make html` to create a HTML documentation (recommended) and `make latexpdf` to create a pdf documentation with Latex.


Python Prerequisites
---------------------

- [Iris](https://github.com/SciTools/iris) (Recommended to install via conda or from source)
- [pandas](https://pandas.pydata.org/)
- [progressbar2](https://pypi.org/project/progressbar2/)
- [seaborn](https://seaborn.pydata.org/)
- [matplotlib](https://matplotlib.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib toolkits](https://matplotlib.org/1.4.3/mpl_toolkits/index.html)
- [netCDF4](https://unidata.github.io/netcdf4-python/netCDF4/index.html)
- [pillow, formerly known as PIL](https://pillow.readthedocs.io/en/stable/)

Sphinx Docs Prerequisites
-------------------------
- [Sphinx](http://www.sphinx-doc.org/en/master/)
- [Read the Docs Sphinx Theme](https://sphinx-rtd-theme.readthedocs.io/en/stable/)
- [recommonmark](https://www.sphinx-doc.org/en/master/usage/markdown.html) to use Markdown

C++ Prerequisites
-----------------
- [GSL](https://www.gnu.org/software/gsl/)
- [NetCDF C++](https://github.com/Unidata/netcdf-cxx4/releases)
- [Boost](https://www.boost.org/)
- [CoDiPack](https://www.scicomp.uni-kl.de/software/codi/)
