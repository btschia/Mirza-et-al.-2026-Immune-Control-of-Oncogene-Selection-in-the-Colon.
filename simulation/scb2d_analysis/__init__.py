"""
Analysis and figure code for the SCB 2D crypt simulation.

Modules
-------
analytical
    Analytical clone size distribution of the one-dimensional stem cell ring (annihilating random walk).
datasets
    Input paths and the assignment of simulations to conditions. The only
    module containing machine-specific paths.
io_simulation
    Readers for the three raw simulation output formats.
monoclonality
    Statistics for P(Monoclonal | Visible).
fixation
    Fixation and extinction probability curves over simulation time.
replacement
    Replacement probability of a labeled stem cell.
experimental_data
    Loader for the measured mouse data the simulations are compared against.
plot_style
    Palette, line styles and the shared matplotlib configuration.
plots
    Plotting functions, each drawing into the current axes.
"""

from . import (
    analytical,
    datasets,
    experimental_data,
    fixation,
    io_simulation,
    monoclonality,
    plot_style,
    plots,
    replacement,
)

__all__ = [
    "analytical",
    "datasets",
    "experimental_data",
    "fixation",
    "io_simulation",
    "monoclonality",
    "plot_style",
    "plots",
    "replacement",
]
