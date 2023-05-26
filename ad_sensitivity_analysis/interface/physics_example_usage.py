"""Plot physics interface results.

Uses physics_t to make some calls to the C++-library and plots the results.
"""
import textwrap

import argparse
import holoviews as hv

from ad_sensitivity_analysis.interface.physics_t import PhysicsT
from ad_sensitivity_analysis.plot import aux_functions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
            """\
            Load the C-library which provides an interface to microphysics processes. Plot some examples of these 
            processes.
            """
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--library_path",
        type=str,
        required=True,
        help=textwrap.dedent(
            """\
            Path to the C-library. You may find it at cmake-build-python_interface/lib/libpython_interface.so.
            """
        ),
    )
    parser.add_argument(
        "--timestep",
        type=float,
        default=30.0,
        help=textwrap.dedent(
            """\
            The timestep size in seconds.
            """
        ),
    )
    parser.add_argument(
        "--uncertainty",
        type=float,
        default=0.01,
        help=textwrap.dedent(
            """\
            The uncertainty that is used for all model parameters.
            """
        ),
    )
    parser.add_argument(
        "--store_path",
        default="../pics/microphysics",
        type=str,
        help=textwrap.dedent(
            """\
            Path to store the generated images.
            """
        ),
    )
    args = parser.parse_args()
    phys = PhysicsT(lib_path=args.library_path)
    phys.setup_model_constants(30.0, 0.1)
    renderer = hv.Store.renderers["matplotlib"].instance(fig="png", dpi=300)
    plot, _ = phys.graupel_melting_plot(
        param_dict={
            "T": 280,
            "qg": 2,
            "Ng": 8,
        },
        x="qg",
        y="final qr",
        hue="T",
        logy=True,
    )
    aux_functions.save_plot_renderer(plot, args.store_path, renderer)
    plot, _ = phys.ccn_act_hande_akm_plot(
        param_dict={
            "w": 10,
            "T": "270",
            "qv": 1e2,
            "qc": -6,
        },
        y="qc",
        logy=True,
    )
    aux_functions.save_plot_renderer(plot, args.store_path, renderer)

    plot, _ = phys.riming_ice_plot(
        param_dict={
            "qc": [0, 1e-3],
            "NC": 1e4,
            "qi": [0, 1e-3],
            "Ni": 0,
            "qr": 1e-6,
            "Nr": 1e3,
            "T": [260, 290],
        },
        x="qc",
        y="qi",
        hue="T",
        logy=False,
        out_param="QV",
    )
    aux_functions.save_plot_renderer(plot, args.store_path, renderer)
