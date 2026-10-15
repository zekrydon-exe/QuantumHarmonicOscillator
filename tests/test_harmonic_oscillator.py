import numpy as np

from src.harmonic_oscillator import solve_harmonic_oscillator
from src.analytical_solution import analytical_energy


def test_ground_state_energy():

    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=1,
    )

    expected = analytical_energy(0)

    assert abs(energies[0] - expected) < 1e-4


def test_first_five_energy_levels():

    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=5,
    )

    expected = np.arange(5) + 0.5

    assert np.allclose(
        energies,
        expected,
        atol=2.1e-3,
        rtol=0.0,
    )


def test_wavefunction_normalisation():

    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=5,
    )

    for n in range(5):

        norm = np.trapezoid(
            np.abs(wavefunctions[:, n])**2,
            q,
        )

        assert abs(norm - 1.0) < 1e-8


def test_wavefunction_orthogonality():

    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=5,
    )

    for i in range(5):

        for j in range(i + 1, 5):

            overlap = np.trapezoid(
                wavefunctions[:, i] * wavefunctions[:, j],
                q,
            )

            assert abs(overlap) < 1e-3


def test_eigenpair_residuals():

    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=5,
    )

    from src.hamiltonian import build_hamiltonian

    H = build_hamiltonian(q, dq)

    # Convert continuum-scaled wavefunctions back to
    # the discrete eigenvectors used by the eigensolver.
    eigenvectors = wavefunctions * np.sqrt(dq)

    residuals = (
        H @ eigenvectors
        - eigenvectors * energies
    )

    residual_norms = np.linalg.norm(residuals, axis=0)

    assert np.max(residual_norms) < 1e-10