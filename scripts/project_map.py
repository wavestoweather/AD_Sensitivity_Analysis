import os
from tqdm.auto import tqdm
import xarray as xr

if __name__ == "__main__":
    import argparse
    import textwrap

    from latexify import *

    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Load data from a sensitivity simulation and calculate column-wise mean, maximum, variance of
            any given parameter or sensitivity.
            Can optionally set minimum and maximum height via --min_pressure or --max_pressure to exclude certain heights.
            Can include a density map of trajectories, phases and filter by phases. 
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file",
        default="../data/vladiana_ensembles/",
        help=textwrap.dedent(
            """\
            Path to a folder with many files from a sensitivity analysis simulation.
            """
        ),
    )
    parser.add_argument(
        "--out_file",
        default="../pics/plots.png",
        help=textwrap.dedent(
            """\
            Path and name to store plots.
            """
        ),
    )
    parser.add_argument(
        "--width",
        default=24,
        type=float,
        help=textwrap.dedent(
            """\
            Width in inches for histogram plots.
            """
        ),
    )
    parser.add_argument(
        "--height",
        default=12,
        type=float,
        help=textwrap.dedent(
            """\
            Height in inches for histogram plots.
            """
        ),
    )
    parser.add_argument(
        "--only_asc600",
        action="store_true",
        help=textwrap.dedent(
            """\
            Consider only time steps during the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--inoutflow_time",
        default=-1,
        type=int,
        help=textwrap.dedent(
            """\
            Consider only time steps during the fastest ascent and within the given range before (inflow) 
            and after (outflow) of the fastest ascent.
            """
        ),
    )
    parser.add_argument(
        "--select_phase",
        default=None,
        type=str,
        help=textwrap.dedent(
            """\
            Consider only the given phase. Options are:
            'warm', 'mixed', 'ice', 'all'. The latter creates a different plot for each phase.
            """
        ),
    )
    parser.add_argument(
        "--plot_type",
        default=["density"],
        type=str,
        nargs="+",
        choices=["density", "mean", "min", "max", "sd"],
        help=textwrap.dedent(
            """\
            Define which kind of plot shall be created. Multiple options can be selected to create multiple plots.
            Options are: 
            'density': Plot the amount of datapoints in each grid point as heatmap. 
            'mean': Plot the mean in each grid point as a heatmap.
            'min': Plot the minimum value in each grid point as a heatmap.
            'max': Plot the maximum value in each grid point as a heatmap.
            'sd': Plot the standard deviation in each grid point as a heatmap.
            """
        ),
    )
    parser.add_argument(
        "--var",
        default=None,
        type=str,
        help=textwrap.dedent(
            """\
            Plot the given variable on a map. If the variable is a sensitivity, then set --sensitivity_of 
            to a model state variable for which the sensitivity is for. 
            Otherwise a plot for each model state variable will be made.
            """
        ),
    )
    parser.add_argument(
        "--sensitivity_of",
        default=None,
        type=str,
        help=textwrap.dedent(
            """\
            If --var is a sensitivity (model parameter), then you may define the model state variable here 
            to plot only sensitivities regarding this model state variable.
            """
        ),
    )
    parser.add_argument(
        "--min_pressure",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            Filter such that only values at this pressure or higher (= at this height or lower) are considered.
            """
        ),
    )
    parser.add_argument(
        "--max_pressure",
        default=None,
        type=float,
        help=textwrap.dedent(
            """\
            Filter such that only values at this pressure or lower (= at this height or higher) are considered.
            """
        ),
    )
    parser.add_argument(
        "--store_calculated",
        default=None,
        type=str,
        help=textwrap.dedent(
            """\
            Store the processed data such that it can be loaded and plotted faster again. 
            """
        ),
    )
    parser.add_argument(
        "--load_calculated",
        default=None,
        type=str,
        help=textwrap.dedent(
            """\
            Load previously processed data to plot it again.
            """
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help=textwrap.dedent(
            """\
            More output regarding the calculations.
            """
        ),
    )
    args = parser.parse_args()
