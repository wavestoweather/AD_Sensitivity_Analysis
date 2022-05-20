# Pre- and Post-Processing
The scripts here are used for either preprocessing data, i.e., changing the format of trajectories such that they can be used by our simulation software, or for postprocessing data, e.g., creating new NetCDF-Files with summarized data, plotting the data, or calculating some statistics.

## Scripts Overview

 - `convert_to_met3d.py`: Loads files in wcb2 style and makes them compatible with Met3D.
 - `dask_loader.py`: Deprecated. Will be removed. File with routines for loading parquet files and filtering or transforming the data to pandas.Dataframe.
 - `Deriv_dask.py`: Deprecated. Will be removed. Class for loading parquet or NetCDF-files of trajectories. The class provides routines to calculate ratios of gradients and plotting routines.
 - `get_stats.py`: Create latex tables and provide statistics for a post-processed file from `segment_identifier.py`. Create histograms and provide statistics for raw output of a sensitivity analysis (a simulation without perturbed ensembles and with gradients).
 - `latexify.py`: Provides routines and lists/dictionaries to map variable names from NetCDF-files to strings formatted to latex. Contains explanations of variables.
 - `merge.py`: Deprecated. Will be removed. Merges multiple NetCDF-files where each file contains only one trajectory.
 - `merge_all.py`: Merge multiple files from ensemble simulations to a single file, where each file contains an ensemble with multiple members where the same parameter had been perturbed.
 - `plot_cosmo.py`: Using trajectories (of an ICON or COSMO simulation), create scatter plots for a quick overview of the data.
 - `plot_mse.py`: Using post-processed data from `segment_identifier.py`, create plots of the mean squared perturbation from ensemble simulations or algorithmic differentiation. The plot can be either a correlation plot or a plot with the predicted squared error over time.
 - `plot_simulation_example.py`: Load trajectories from an ICON or COSMO simulation and the resimulation with our software and create scatter plots of various model state variables for comparison.
 - `segment_identifier.py`: Post-processing of perturbed ensemble simulations with sensitivity analysis. Creates a file with mean predicted errors and mean errors from ensemble simulation. Can also transform data for different machine learning algorithms.
 - `test_output.py`: Tests a given NetCDF-file for suspicious data. This includes checks for NaNs and hydrometeor sizes.

For a thorough description, you may use
```
python script.py --help
```