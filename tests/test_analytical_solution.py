import numpy as np
import pytest

from src.analytical_solution import (
    analytical_energy,
    analytical_wavefunction,
)


def test_analytical_energies():

    expected = np.array([0.5, 1.5, 2.5, 3.5, 4.5])

    calculated = np.array([
        analytical_energy(n)
        for n in range(5)
    ])

    assert np.allclose(calculated, expected)


def test_ground_state_wavefunction():

    q = np.array([0.0])

    psi = analytical_wavefunction(q, 0)

    expected = np.pi**(-0.25)

    assert np.isclose(psi[0], expected)


def test_wavefunction_parity():

    q = np.linspace(0.1, 3.0, 100)

    for n in range(5):

        psi_positive = analytical_wavefunction(q, n)

        psi_negative = analytical_wavefunction(-q, n)

        expected = (-1)**n * psi_positive

        assert np.allclose(psi_negative, expected)


def test_negative_quantum_number():

    with pytest.raises(ValueError):

        analytical_energy(-1)


def test_non_integer_quantum_number():

    with pytest.raises(TypeError):

        analytical_wavefunction(
            np.array([0.0, 1.0]),
            1.5
        )