"""
Plotting functions for the SCB 2D simulation figures.

Every function draws into the *current* matplotlib axes, so that the calling
script retains control over figure size, layout and export.
`plot_style.setup_plot_style` is to be called after the figure has been created
and before plotting into it.
"""

import matplotlib.pyplot as plt
import numpy as np

from .plot_style import LINE_STYLES


def plot_probability_over_days(mean_by_day, std_by_day, line_style, color, label,
                               band_alpha=0.3):
    """
    Plot P(Monoclonal | Visible) against the simulation day, with a ±1 std band.

    Parameters
    ----------
    mean_by_day, std_by_day : dict
        Output of `monoclonality.average_probability_per_day`.
    line_style, color : str, color
        Style of the mean curve; the band uses the same color.
    label : str
        Legend entry; pass "" to keep the curve out of the legend.
    band_alpha : float
        Opacity of the band; 0 draws the mean curve without a band.
    """
    days = sorted(mean_by_day.keys())
    means = np.array([mean_by_day[day] for day in days])
    stds = np.array([std_by_day[day] for day in days])

    # The simulation counts from day 0, the experiment from day 1.
    days = np.array(days) + 1

    plt.plot(days, means, linestyle=line_style, color=color, label=label)
    if band_alpha:
        plt.fill_between(days, means - stds, means + stds, color=color, alpha=band_alpha)


def plot_mean_with_std_band(days, means, stds,
                            band_color="skyblue", line_color="steelblue",
                            label="Mean", line_style=LINE_STYLES[2], band_alpha=0.4):
    """
    Plot a cumulative probability curve and shade the band between mean ± std.

    Parameters
    ----------
    days, means, stds : sequence
        Output of `fixation.average_probability_curves`.
    band_color, line_color : color
        Fill color of the band and color of the mean curve.
    label : str
        Legend entry.
    line_style : str
        Line style of the mean curve.
    band_alpha : float
        Opacity of the band.
    """
    # The curve is anchored at day 0 with probability 0, so that it starts at
    # the origin.
    days = np.array([0] + list(days))
    means = np.array([0] + list(means))
    stds = np.array([0] + list(stds))

    plt.plot(days, means, color=line_color, linestyle=line_style, label=label)
    plt.fill_between(days, means - stds, means + stds,
                     color=band_color, alpha=band_alpha)


def plot_replacement_probability_histogram(run_means, bar_color,
                                           label=None, num_bins=20,
                                           mark_neutral_line=False,
                                           x_limits=(0, 1), y_limits=(0, 10)):
    """
    Plot the distribution of per-run replacement probabilities as a density histogram.

    Parameters
    ----------
    run_means : sequence
        Mean replacement probability of each simulation run.
    bar_color : color
        Fill color of the bars.
    label : str
        Legend entry.
    num_bins : int
        Number of histogram bins.
    mark_neutral_line : bool
        Draw a dashed line at 0.5, i.e. neutral competition between the labeled
        and the unlabeled lineage.
    x_limits, y_limits : tuple
        Axis limits of the panel.
    """
    plt.hist(np.array(run_means), bins=num_bins, density=True,
             color=bar_color, edgecolor="black", label=label, alpha=0.7)

    if mark_neutral_line:
        plt.axvline(x=0.5, color="grey", linestyle="dashed", linewidth=2)

    plt.xlabel("Replacement Probability")
    plt.ylabel("Density")
    plt.xlim(*x_limits)
    plt.ylim(*y_limits)
    plt.grid(False)


def plot_experimental_points(days, averages_percent, std_percent, color):
    """
    Overlay the measured fraction of fixed crypts as error bars.

    `averages_percent` and `std_percent` are given in percent and are divided by
    100 here, so that they share the axis with the simulated probabilities.
    """
    plt.errorbar(days, np.asarray(averages_percent) / 100,
                 yerr=np.asarray(std_percent) / 100,
                 fmt="o", capsize=4, color=color, alpha=1)


def reversed_legend(**legend_kwargs):
    """
    Draw the legend of the current axes in reverse order.

    Matplotlib lists entries in plotting order; reversing them makes the legend
    match the top-to-bottom order of the curves in these panels.
    """
    handles, labels = plt.gca().get_legend_handles_labels()
    return plt.legend(handles[::-1], labels[::-1], **legend_kwargs)
