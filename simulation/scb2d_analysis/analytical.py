"""
Analytical clone size distribution of the one-dimensional stem cell ring (annihilating random walk).

The stem cells of a crypt are described as a ring of `num_stem_cells` positions
in which a cell is replaced by one of its two neighbours at rate
`replacement_rate`. The number of labeled cells performs a random walk on
0 ... `num_stem_cells`, with the two absorbing boundaries

extinction
    0 labeled cells, i.e. the labeled clone has been lost.
fixation
    `num_stem_cells` labeled cells, i.e. the crypt has become monoclonal for the
    labeled clone.

Both functions below return, for every requested time point,

p_n
    Probability of observing exactly n labeled cells, for n = 1 ... N - 1.
p_0
    Extinction probability.
p_ns
    Fixation probability.
p_persisting
    Clone size distribution conditioned on survival, p_n / (1 - p_0). This is
    the quantity the experiment observes, as a crypt without a labeled cell is
    not detected.

`compute_clone_distribution` treats the two replacement directions as equally
likely. `compute_clone_distribution_asymmetric` introduces the bias `delta`
between them, by which a labeled cell replaces a neighbour more often than the
reverse.

The rates are per day, so that the time argument is in days and the results
share the axis with the simulated and the measured curves.
"""

import numpy as np


def compute_clone_distribution(num_stem_cells, replacement_rate, time):
    """
    Clone size distribution of the symmetric ring, in which both replacement
    directions are equally likely.

    Parameters
    ----------
    num_stem_cells : int
        Number of stem cell positions of the ring, N.
    replacement_rate : float
        Replacement rate per day, lambda.
    time : float or sequence
        Time points in days. A scalar yields scalar results, a sequence yields
        one array per quantity.

    Returns
    -------
    p_n : dict
        n -> probability of exactly n labeled cells, for n = 1 ... N - 1.
    p_0 : ndarray or float
        Extinction probability.
    p_ns : ndarray or float
        Fixation probability.
    p_persisting : dict
        n -> p_n / (1 - p_0), the distribution conditioned on survival, for
        n = 1 ... N.
    """
    time_points = np.atleast_1d(time)

    p_0_history = []
    p_ns_history = []
    p_n_history = {n: [] for n in range(1, num_stem_cells)}
    p_persisting_history = {n: [] for n in range(1, num_stem_cells + 1)}

    for time_point in time_points:
        mode = np.arange(1, num_stem_cells)

        # Eigenvalues of the discrete Laplacian of the ring and their decay.
        eigenvalues = 4 * (np.sin(np.pi * mode / (2 * num_stem_cells)) ** 2)
        decay = np.exp(-eigenvalues * replacement_rate * time_point)
        cosine_squared = np.cos(np.pi * mode / (2 * num_stem_cells)) ** 2

        # Extinction probability.
        p_0 = (2 / num_stem_cells) * np.sum(cosine_squared * (1 - decay))
        p_0_history.append(p_0)

        # Fixation probability; the alternating sign is the second boundary.
        alternating_sign = (-1) ** (mode + 1)
        p_ns = (2 / num_stem_cells) * np.sum(
            alternating_sign * cosine_squared * (1 - decay))
        p_ns_history.append(p_ns)

        # Probability of the intermediate clone sizes.
        for n in range(1, num_stem_cells):
            p_n = (2 / num_stem_cells) * np.sum(
                np.sin(np.pi * mode / num_stem_cells)
                * np.sin(np.pi * mode * n / num_stem_cells)
                * decay)

            p_n_history[n].append(p_n)
            p_persisting_history[n].append(p_n / (1 - p_0) if (1 - p_0) > 0 else 0.0)

        p_persisting_history[num_stem_cells].append(
            p_ns / (1 - p_0) if (1 - p_0) > 0 else 0.0)

    return _as_arrays(p_n_history, p_0_history, p_ns_history,
                      p_persisting_history, scalar_input=np.isscalar(time))


def compute_clone_distribution_asymmetric(num_stem_cells, replacement_rate, delta, time):
    """
    Clone size distribution of the ring with a bias between the two replacement
    directions.

    Parameters
    ----------
    num_stem_cells : int
        Number of stem cell positions of the ring, N.
    replacement_rate : float
        Replacement rate per day, lambda.
    delta : float
        Bias between the two replacement directions, -1 < delta < 1. A positive
        value favours the labeled clone; delta = 0 reproduces
        `compute_clone_distribution`.
    time : float or sequence
        Time points in days.

    Returns
    -------
    p_n, p_0, p_ns, p_persisting
        As in `compute_clone_distribution`.
    """
    time_points = np.atleast_1d(time)

    # The two replacement directions have the rates lambda * (1 + delta) and
    # lambda * (1 - delta). The solution depends on them only through their
    # ratio, which weights the clone sizes, and their geometric mean, which
    # enters the decay rates below.
    bias_ratio = np.sqrt((1 + delta) / (1 - delta))
    bias_scale = np.sqrt(1 - delta ** 2)

    p_0_history = []
    p_ns_history = []
    p_n_history = {n: [] for n in range(1, num_stem_cells)}
    p_persisting_history = {n: [] for n in range(1, num_stem_cells + 1)}

    for time_point in time_points:
        mode = np.arange(1, num_stem_cells)

        # Eigenvalues of the biased walk, replacing those of the symmetric ring.
        eigenvalues = (2 * (1 / bias_scale - 1)
                       + 4 * (np.sin(np.pi * mode / (2 * num_stem_cells)) ** 2))
        decay = np.exp(-bias_scale * replacement_rate * time_point * eigenvalues)
        sine_squared = np.sin(np.pi * mode / num_stem_cells) ** 2

        # Extinction probability.
        p_0 = (2 / (num_stem_cells * bias_ratio)) * np.sum(
            (sine_squared / eigenvalues) * (1 - decay))
        p_0_history.append(p_0)

        # Fixation probability.
        alternating_sign = (-1) ** (mode + 1)
        p_ns = (2 / num_stem_cells) * (bias_ratio ** (num_stem_cells - 1)) * np.sum(
            (alternating_sign * sine_squared / eigenvalues) * (1 - decay))
        p_ns_history.append(p_ns)

        # Probability of the intermediate clone sizes.
        for n in range(1, num_stem_cells):
            p_n = (2 / num_stem_cells) * (bias_ratio ** (n - 1)) * np.sum(
                np.sin(np.pi * mode / num_stem_cells)
                * np.sin(np.pi * mode * n / num_stem_cells)
                * decay)

            p_n_history[n].append(p_n)
            p_persisting_history[n].append(p_n / (1 - p_0) if (1 - p_0) > 0 else 0.0)

        p_persisting_history[num_stem_cells].append(
            p_ns / (1 - p_0) if (1 - p_0) > 0 else 0.0)

    return _as_arrays(p_n_history, p_0_history, p_ns_history,
                      p_persisting_history, scalar_input=np.isscalar(time))


def _as_arrays(p_n_history, p_0_history, p_ns_history, p_persisting_history,
               scalar_input):
    """
    Convert the per-time-point lists into arrays, or into scalars if the caller
    passed a single time point.
    """
    p_0 = np.array(p_0_history)
    p_ns = np.array(p_ns_history)
    p_n = {n: np.array(values) for n, values in p_n_history.items()}
    p_persisting = {n: np.array(values) for n, values in p_persisting_history.items()}

    if scalar_input:
        return ({n: values[0] for n, values in p_n.items()}, p_0[0], p_ns[0],
                {n: values[0] for n, values in p_persisting.items()})

    return p_n, p_0, p_ns, p_persisting
