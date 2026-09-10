"""
Replacement probability of a labeled stem cell.

For every crypt the simulation records how often a labeled stem cell was
replaced by a labeled or the unlabeled cell replaced the labeled cell. The replacement probability is

    labeled / (labeled + unlabeled)

A value of 0.5 means labeled and unlabeled cells compete neutrally; values above
0.5 indicate an advantage of the labeled (mutant) lineage.

The summaries are computed on two levels:

run means
    Averaged within each simulation run first, then across runs. This is the
    level reported in the manuscript, as the runs are the independent units.
pooled
    Every individual crypt is treated as one observation.
"""

import numpy as np
import scipy.stats as stats

from .io_simulation import NOT_MONOCLONAL_MARKER


def compute_replacement_probabilities(simulation_runs, max_crypts_per_run=None):
    """
    Compute the labeled replacement probability of every crypt, per run.

    Parameters
    ----------
    simulation_runs : list
        Output of `io_simulation.read_replacement_counts_from_directory`.
    max_crypts_per_run : int or None
        If set, only the first N crypts of each run are used.

    Returns
    -------
    probabilities_all_crypts : dict
        run index -> list of probabilities, for every crypt (monoclonal or not).
    probabilities_monoclonal_crypts : dict
        run index -> list of probabilities, restricted to the crypts that
        actually became monoclonal.
    """
    probabilities_all_crypts = {}
    probabilities_monoclonal_crypts = {}

    for run_index, run in enumerate(simulation_runs):
        probabilities_all_crypts[run_index] = []
        probabilities_monoclonal_crypts[run_index] = []

        crypt_records = run[:max_crypts_per_run] if max_crypts_per_run else run

        for record in crypt_records:
            if len(record) == 2:
                # [labeled, unlabeled] -> the crypt became monoclonal
                labeled, unlabeled = record
                became_monoclonal = True
            elif len(record) == 3 and record[0] == NOT_MONOCLONAL_MARKER:
                # [marker, labeled, unlabeled] -> the crypt did not become monoclonal
                labeled, unlabeled = record[1], record[2]
                became_monoclonal = False
            else:
                # Records of unexpected length are ignored.
                continue

            total_replacements = labeled + unlabeled
            if total_replacements == 0:
                continue

            probability = labeled / total_replacements
            probabilities_all_crypts[run_index].append(probability)
            if became_monoclonal:
                probabilities_monoclonal_crypts[run_index].append(probability)

    return probabilities_all_crypts, probabilities_monoclonal_crypts


def summarize_run_means(probabilities_by_run):
    """
    Summarise replacement probabilities without confidence intervals.

    Parameters
    ----------
    probabilities_by_run : dict
        run index -> list of per-crypt replacement probabilities.

    Returns
    -------
    run_means : list
        Mean replacement probability of each run.
    mean_of_run_means : float
    std_of_run_means : float
        Population standard deviation (ddof=0) of the run means.
    pooled_mean : float
    pooled_std : float
        Mean and population standard deviation over all individual crypts.
    """
    # One mean per run; runs without a valid crypt are skipped.
    run_means = [np.mean(probabilities)
                 for probabilities in probabilities_by_run.values() if probabilities]

    # All individual crypt probabilities.
    all_probabilities = [probability
                         for probabilities in probabilities_by_run.values()
                         for probability in probabilities]

    return (run_means,
            np.mean(run_means), np.std(run_means),
            np.mean(all_probabilities), np.std(all_probabilities))


def summarize_run_means_with_ci(probabilities_by_run, confidence=0.95):
    """
    Summarise replacement probabilities including confidence intervals.

    Parameters
    ----------
    probabilities_by_run : dict
        run index -> list of per-crypt replacement probabilities.
    confidence : float
        Confidence level of the intervals.

    Returns
    -------
    run_means : list
        Mean replacement probability of each run.
    mean_of_run_means : float
    std_of_run_means : float
        Sample standard deviation (ddof=1) of the run means.
    ci_of_run_means : tuple
        Confidence interval of `mean_of_run_means`.
    pooled_mean : float
    pooled_std : float
    pooled_ci : tuple
        The same three quantities computed over all individual crypts.
    """
    run_means = [np.mean(probabilities)
                 for probabilities in probabilities_by_run.values() if probabilities]
    num_runs = len(run_means)

    all_probabilities = [probability
                         for probabilities in probabilities_by_run.values()
                         for probability in probabilities]
    num_crypts = len(all_probabilities)

    mean_of_run_means = np.mean(run_means)
    std_of_run_means = np.std(run_means, ddof=1)
    pooled_mean = np.mean(all_probabilities)
    pooled_std = np.std(all_probabilities, ddof=1)

    ci_of_run_means = stats.t.interval(
        confidence, df=num_runs - 1,
        loc=mean_of_run_means, scale=std_of_run_means / np.sqrt(num_runs))
    pooled_ci = stats.t.interval(
        confidence, df=num_crypts - 1,
        loc=pooled_mean, scale=pooled_std / np.sqrt(num_crypts))

    return (run_means, mean_of_run_means, std_of_run_means, ci_of_run_means,
            pooled_mean, pooled_std, pooled_ci)


def format_run_mean_summary(mean_of_run_means, std_of_run_means, ci_of_run_means):
    """Format the run-level summary of one condition as a single line."""
    return (f"mean of run means: {mean_of_run_means:.3f}, "
            f"std: {std_of_run_means:.3f}, "
            f"95% CI: ({ci_of_run_means[0]:.3f}, {ci_of_run_means[1]:.3f})")
