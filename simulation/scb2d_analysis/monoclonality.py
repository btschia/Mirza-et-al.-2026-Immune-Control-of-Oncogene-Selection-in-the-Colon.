"""
Statistics for P(Monoclonal | Visible), the probability that a crypt is
monoclonal given that it carries a visible label.
"""

import numpy as np


def average_probability_per_day(probabilities_by_day):
    """
    Average P(Monoclonal | Visible) over the simulation repetitions of each day.

    Parameters
    ----------
    probabilities_by_day : dict
        Output of `io_simulation.read_visible_monoclonal_probabilities_from_directory`,
        i.e. simulation day -> list of values, one per repetition.

    Returns
    -------
    mean_by_day : dict
        day -> mean over repetitions.
    std_by_day : dict
        day -> sample standard deviation (ddof=1) over repetitions.
    """
    mean_by_day = {}
    std_by_day = {}

    for day, probabilities in probabilities_by_day.items():
        mean_by_day[day] = np.mean(probabilities)
        std_by_day[day] = np.std(probabilities, ddof=1)

    return mean_by_day, std_by_day
