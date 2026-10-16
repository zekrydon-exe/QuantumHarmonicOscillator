import numpy as np

from src.harmonic_oscillator import solve_harmonic_oscillator


def test_consecutive_energy_levels_are_equally_spaced():
    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=5,
    )

    energy_gaps = np.diff(energies)

    assert np.allclose(
        energy_gaps,
        1.0,
        atol=1e-3,
        rtol=0.0,
    )


def test_wavefunctions_have_correct_parity():
    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=5,
    )

    for n in range(5):
        expected = (-1) ** n * wavefunctions[::-1, n]

        assert np.allclose(
            wavefunctions[:, n],
            expected,
            atol=1e-10,
            rtol=1e-8,
        )


def test_probability_density_is_nonnegative_and_symmetric():
    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=5,
    )

    probability_densities = np.abs(wavefunctions) ** 2

    assert np.all(probability_densities >= 0.0)

    for n in range(5):
        assert np.allclose(
            probability_densities[:, n],
            probability_densities[::-1, n],
            atol=1e-10,
            rtol=1e-8,
        )