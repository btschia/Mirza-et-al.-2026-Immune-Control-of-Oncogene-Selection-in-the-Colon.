"""
Fixation and extinction of labeled lineages over time.

Terminology used throughout this module:

run
    One simulation repetition, covering 100 crypts simulated in parallel.
event
    One crypt becoming monoclonal, recorded as ``[day, lineage_id]``.
day
    Simulation day, shifted by +1 so that the first day is day 1. The
    simulation itself writes 0 for the first day.
fixation
    A crypt colonized by the *labeled* lineage.
extinction
    A crypt colonized by *unlabeled* cells, i.e. the labeled
    lineage in that crypt is pushed out.
"""

import numpy as np
import scipy.stats as stats

from .io_simulation import UNLABELED_LINEAGE_ID

# Number of crypts simulated in parallel within one run; the denominator of the
# per-run probabilities computed below.
CRYPTS_PER_RUN = 100


def _group_events_by_day(simulation_runs, max_day, keep_unlabeled):
    """
    Bucket the conversion events of every run by simulation day.

    Parameters
    ----------
    simulation_runs : list
        Output of `io_simulation.read_fixation_events_from_directory`.
    max_day : int
        Last day to keep; events after this day are dropped.
    keep_unlabeled : bool
        True  -> keep only the crypts taken over by an unlabeled cell (extinction).
        False -> keep only the crypts taken over by a labeled lineage (fixation).

    Returns
    -------
    list of dict
        One dict per run, mapping day (1 ... max_day) -> list of lineage ids.
    """
    events_per_run = []

    for run in simulation_runs:
        # Every day is pre-filled, so that days without an event remain empty.
        events_by_day = {day: [] for day in range(1, max_day + 1)}

        for day_lineage_pair in run:
            day = day_lineage_pair[0] + 1  # the simulation counts from day 0
            lineage_id = day_lineage_pair[1]

            is_unlabeled = lineage_id == UNLABELED_LINEAGE_ID
            if is_unlabeled == keep_unlabeled and 1 <= day <= max_day:
                events_by_day[day].append(lineage_id)

        events_per_run.append(events_by_day)

    return events_per_run


def group_fixation_events_by_day(simulation_runs, max_day=100):
    """Bucket the fixation events (labeled lineage takes over a crypt) by day."""
    return _group_events_by_day(simulation_runs, max_day, keep_unlabeled=False)


def group_extinction_events_by_day(simulation_runs, max_day=100):
    """Bucket the extinction events (unlabeled cell takes over a crypt) by day."""
    return _group_events_by_day(simulation_runs, max_day, keep_unlabeled=True)


def compute_cumulative_probability_per_run(events_by_day_per_run,
                                           crypts_per_run=CRYPTS_PER_RUN):
    """
    Turn per-day event counts into a cumulative probability curve for each run.

    The probability on a given day is the number of crypts that had the event up
    to and including that day, divided by the number of crypts in the run.

    Parameters
    ----------
    events_by_day_per_run : list of dict
        Output of `group_fixation_events_by_day` or `group_extinction_events_by_day`.
    crypts_per_run : int
        Number of crypts simulated in parallel per run (the denominator).

    Returns
    -------
    list of dict
        One dict per run, mapping day -> cumulative probability.
    """
    cumulative_curves = []

    for events_by_day in events_by_day_per_run:
        days = sorted(events_by_day.keys())

        # Per-day (non-cumulative) probability.
        daily_probability = [len(events_by_day[day]) / crypts_per_run for day in days]

        # Accumulated over days: once a crypt has had the event, it remains counted.
        cumulative_probability = np.cumsum(daily_probability)

        cumulative_curves.append(dict(zip(days, cumulative_probability)))

    return cumulative_curves


def average_probability_curves(cumulative_curves):
    """
    Average the per-run cumulative curves day by day.

    Returns
    -------
    days : list
        Sorted simulation days.
    means : list
        Mean probability across runs for each day.
    stds : list
        Population standard deviation (ddof=0) across runs for each day.
    """
    days = []
    means = []
    stds = []

    for day in sorted(cumulative_curves[0].keys()):
        values_across_runs = [curve.get(day, 0) for curve in cumulative_curves]
        days.append(day)
        means.append(np.mean(values_across_runs))
        stds.append(np.std(values_across_runs))

    return days, means, stds


def average_probability_curves_with_ci(cumulative_curves, confidence=0.95):
    """
    Like `average_probability_curves`, but with a confidence interval of the mean
    and the *sample* standard deviation (ddof=1).

    Returns
    -------
    days : list
        Sorted simulation days.
    means : list
        Mean probability across runs for each day.
    stds : list
        Sample standard deviation across runs for each day.
    confidence_intervals : list of tuple
        (lower, upper) bound of the confidence interval of the mean, per day.
    """
    days = []
    means = []
    stds = []
    confidence_intervals = []

    num_runs = len(cumulative_curves)

    for day in sorted(cumulative_curves[0].keys()):
        values_across_runs = [curve.get(day, 0) for curve in cumulative_curves]

        mean = np.mean(values_across_runs)
        std = np.std(values_across_runs, ddof=1)
        standard_error = std / np.sqrt(num_runs)

        days.append(day)
        means.append(mean)
        stds.append(std)
        confidence_intervals.append(
            stats.t.interval(confidence, df=num_runs - 1, loc=mean, scale=standard_error)
        )

    return days, means, stds, confidence_intervals


def format_probability_at_day(days, means, stds, confidence_intervals, day):
    """
    Format the mean, standard deviation and CI of a single day as one line.

    An empty string is returned if `day` is not part of the simulated period, so
    that reporting a day beyond `max_day` cannot raise.
    """
    if day not in days:
        return ""

    index = days.index(day)
    lower, upper = confidence_intervals[index]
    return (f"day {day}: mean {means[index]:.3f}, std {stds[index]:.3f}, "
            f"CI ({lower:.3f}, {upper:.3f})")
