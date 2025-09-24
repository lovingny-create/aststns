import math

from src.decay import simulate_decay


def test_half_life_value_close_to_expected():
    T12 = 10.0
    N0 = 1000.0
    dt = 0.5
    T_total = 10.0

    times, populations = simulate_decay(T12, N0, dt, T_total)

    # Identify the sample closest to t = 10 seconds.
    closest_index = min(range(len(times)), key=lambda i: abs(times[i] - 10.0))
    expected = N0 * math.exp(-math.log(2.0))
    assert math.isclose(populations[closest_index], expected, rel_tol=1e-3)


def test_monotonic_decrease():
    T12 = 5.0
    N0 = 500.0
    dt = 0.25
    T_total = 20.0

    _, populations = simulate_decay(T12, N0, dt, T_total)
    differences = [populations[i + 1] - populations[i] for i in range(len(populations) - 1)]
    assert all(diff <= 0 for diff in differences)
