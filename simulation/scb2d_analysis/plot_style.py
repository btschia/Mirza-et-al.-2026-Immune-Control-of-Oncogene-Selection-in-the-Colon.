"""
Shared figure styling for the SCB 2D simulation figures.

The palette and the rcParams are defined in one place, so that every panel of
the manuscript uses the same colors, line styles and font sizes.
"""

import os

import matplotlib.pyplot as plt
from cycler import cycler

# Conversion factor from centimetres to inches, for `figsize=(7.65 * CM, ...)`.
CM = 1 / 2.54

# Colorblind-friendly qualitative palette (adapted from Okabe-Ito).
COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"]
LINE_STYLES = ["-", "--", "-.", ":", "-"]

# Colors used for the published panels.
PAPER_RED = (209 / 255, 87 / 255, 87 / 255)     # BRAF
PAPER_GREEN = (135 / 255, 196 / 255, 185 / 255)
PAPER_BLUE = COLORS[3]

# Cyclers for figures containing many conditions at once.
LINE_CYCLER = (
    cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2",
                  "#D55E00", "#CC79A7", "#F0E442"])
    + cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."])
)
MARKER_CYCLER = (
    cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2",
                  "#D55E00", "#CC79A7", "#F0E442"])
    + cycler(linestyle=["none"] * 7)
    + cycler(marker=["4", "2", "3", "1", "+", "x", "."])
)

# Shade families for the sensitivity figures. One family is assigned per
# condition: the base shade denotes the reference run, the lighter and the
# darker shade the two variants (e.g. kr doubled and halved).
ORANGE_SHADES = (COLORS[1], "#F2C55C", "#B87C00")
RED_SHADES = (PAPER_RED, "#E88D8D", "#A63C3C")
BLUE_SHADES = (COLORS[0], "#6FB8E0", "#004C77")
GREEN_SHADES = (COLORS[2], "#66C7AC", "#00674B")

SHADE_FAMILIES = [ORANGE_SHADES, RED_SHADES, BLUE_SHADES, GREEN_SHADES]


def shade_family(index):
    """
    Return the (base, light, dark) shades of the `index`-th condition.

    The families are cycled if a figure contains more conditions than there are
    families.
    """
    return SHADE_FAMILIES[index % len(SHADE_FAMILIES)]


# Axis label used for the monoclonality figures.
P_MONOCLONAL_GIVEN_VISIBLE_LABEL = r"$P(Monoclonal \mid Visible)$"


# Formats a panel is written in when a notebook exports its figures.
FIGURE_FORMATS = ("pdf",)


def save_figure(name, output_directory, formats=FIGURE_FORMATS):
    """
    Write the current figure to `output_directory`, one file per format.

    Used by the notebooks to export their panels, so that the file a panel ends
    up in is documented in the notebook itself rather than in a path edited by
    hand.

    Parameters
    ----------
    name : str
        File name without the extension; identifies the panel.
    output_directory : str or Path
        Target directory, created if it does not exist.
    formats : iterable of str
        File formats to write, PDF by default.

    Returns
    -------
    list of str
        The files written.
    """
    os.makedirs(output_directory, exist_ok=True)

    written = []
    for file_format in formats:
        path = os.path.join(str(output_directory), f"{name}.{file_format}")
        plt.savefig(path, bbox_inches="tight")
        written.append(path)

    return written


def setup_plot_style(font_size=14, font_family="Arial", axes=None):
    """
    Apply the manuscript style to one or several panels.

    The grid is removed and the top and right spines are hidden. To be called
    directly after a figure has been created and before plotting into it.

    The font settings are written to the rcParams, which matplotlib reads when a
    figure is created. They therefore apply to the figures created *after* this
    call, not to the one this call styles; the notebooks set the font once in
    their first cell, so that every panel is rendered at the same size.

    Parameters
    ----------
    font_size, font_family : float, str
        Font of the figures created after this call.
    axes : matplotlib axes or iterable of axes, optional
        Panels to style; by default the current axes. Both panels of a figure
        with two plots beside each other are styled by passing them as a pair,
        as only one of them is the current axes.

    Returns
    -------
    The styled axes, or the list of them if several were given.
    """
    plt.rcParams.update({"font.size": font_size, "font.family": font_family})
    plt.rcParams["axes.linewidth"] = 1

    if axes is None:
        panels = [plt.gca()]
    elif hasattr(axes, "plot"):  # a single axes rather than a sequence
        panels = [axes]
    else:
        panels = list(axes)

    for panel in panels:
        panel.grid(False)
        panel.spines["top"].set_visible(False)
        panel.spines["right"].set_visible(False)

    return panels[0] if len(panels) == 1 else panels
