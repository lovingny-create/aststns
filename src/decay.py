"""Decay simulation utilities."""
from __future__ import annotations

import math
from typing import List, Tuple


def simulate_decay(T12: float, N0: float, dt: float, T_total: float) -> Tuple[List[float], List[float]]:
    """Simulate the decay of a single nuclide population.

    Parameters
    ----------
    T12 : float
        Half-life of the nuclide in seconds.
    N0 : float
        Initial number of nuclei.
    dt : float
        Time interval between samples in seconds. Must be positive.
    T_total : float
        Total duration of the simulation in seconds. Must be non-negative.

    Returns
    -------
    Tuple[List[float], List[float]]
        A tuple ``(t, N)`` containing the sampled time points and the
        corresponding number of nuclei at each time point.

    Notes
    -----
    The simulation assumes ideal exponential decay with decay constant
    ``lambda = ln(2) / T12``. The numerical solution samples the analytic
    expression ``N(t) = N0 * exp(-lambda * t)`` at regular intervals defined
    by ``dt``. The returned sequence is monotonically non-increasing.
    """

    if T12 <= 0:
        raise ValueError("Half-life T12 must be positive.")
    if dt <= 0:
        raise ValueError("Time step dt must be positive.")
    if T_total < 0:
        raise ValueError("Total simulation time must be non-negative.")

    decay_constant = math.log(2.0) / T12

    num_steps = int(math.floor(T_total / dt))
    times = [i * dt for i in range(num_steps + 1)]
    if not times or not math.isclose(times[-1], T_total):
        times.append(T_total)

    populations = [N0 * math.exp(-decay_constant * t) for t in times]
    return times, populations
